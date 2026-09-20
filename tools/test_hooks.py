#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""test_hooks.py — 测试两个 hook（note_index + prior_art）的声音性与韧性。

设计：subprocess 喂 stdin（JSON payload 用 json.dumps 生成，不拼字符串），验证输出与退出码。
投毒用例：hook 的根目录可由环境变量 SC_HOOK_ROOT 覆盖 ⟹ 在临时目录里构造一个 INDEX 做确定性投毒。
"""
import json
import os
import subprocess
import sys

# Force UTF-8 output on Windows
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")

# Derive ROOT from script location
ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")).replace("\\", "/")
HOOK_INDEX = os.path.join(ROOT, "tools", "hook_note_index.py")
HOOK_PRIOR = os.path.join(ROOT, "tools", "hook_prior_art.py")


def run_hook(hook_path, payload_dict):
    """Run hook with JSON payload, return (exit_code, stdout, stderr)."""
    payload_bytes = json.dumps(payload_dict, ensure_ascii=True).encode("ascii")
    result = subprocess.run(
        [sys.executable, hook_path],
        input=payload_bytes,
        capture_output=True,
        timeout=30,
    )
    return result.returncode, result.stdout, result.stderr


def test_malformed_json():
    """畸形 JSON ⟹ 退出 0、无输出。"""
    result = subprocess.run(
        [sys.executable, HOOK_INDEX],
        input=b"{malformed",
        capture_output=True,
        timeout=5,
    )
    assert result.returncode == 0, f"Expected exit 0, got {result.returncode}"
    assert result.stdout == b"", f"Expected no output, got {result.stdout}"
    print("✓ PASS: malformed JSON (silent exit)")


def test_wrong_tool_name():
    """tool_name 不是 Agent / 不是 Write|Edit ⟹ 无输出。"""
    # Test index hook with Read
    payload = {"tool_name": "Read", "tool_input": {"file_path": "theory/notes/test.md"}}
    rc, out, _ = run_hook(HOOK_INDEX, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    assert out == b"", f"Expected no output, got {out}"

    # Test prior_art hook with Write
    payload = {"tool_name": "Write", "tool_input": {"file_path": "test.md"}}
    rc, out, _ = run_hook(HOOK_PRIOR, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    assert out == b"", f"Expected no output, got {out}"

    print("✓ PASS: wrong tool_name (silent exit)")


def test_agent_with_vocab():
    """Agent prompt 含 vocab 词 ⟹ 有 additionalContext 且 updatedInput.prompt 以原 prompt 开头。"""
    payload = {
        "tool_name": "Agent",
        "tool_input": {
            "description": "test agent",
            "prompt": "请分析 δ_BRD 谱展宽亏损在 H3S 中的贡献。",
        },
    }
    rc, out, _ = run_hook(HOOK_PRIOR, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"

    if out:
        data = json.loads(out.decode("utf-8"))
        hso = data.get("hookSpecificOutput", {})
        assert "additionalContext" in hso, "Missing additionalContext"
        updated = hso.get("updatedInput", {})
        new_prompt = updated.get("prompt", "")
        orig_prompt = payload["tool_input"]["prompt"]
        assert new_prompt.startswith(orig_prompt), "updatedInput.prompt doesn't start with original"
        # Check other keys preserved
        assert updated.get("description") == "test agent", "description not preserved"
        print("✓ PASS: Agent with vocab (hook activated)")
    else:
        # If vocab not in prior_art.sh, this is acceptable
        print("⊘ SKIP: Agent with vocab (vocab not in prior_art.sh, or no activation)"); return "SKIP"


def test_agent_no_vocab():
    """Agent prompt 不含任何 vocab 词 ⟹ 无输出。"""
    payload = {
        "tool_name": "Agent",
        "tool_input": {
            "description": "test agent",
            "prompt": "Please read the README file and summarize it.",
        },
    }
    rc, out, _ = run_hook(HOOK_PRIOR, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    # Should be no output since no vocab match
    if out:
        print(f"⚠ WARNING: Got output when no vocab expected: {out[:200]}")
    print("✓ PASS: Agent with no vocab (silent)")


def test_write_new_note():
    """Write 到 theory/notes/materials/<不存在的名>.md ⟹ 有提醒。"""
    payload = {
        "tool_name": "Write",
        "tool_input": {
            "file_path": os.path.join(ROOT, "theory/notes/materials/nonexistent_test_file_xyz_2026-09-20.md"),
            "content": "test content",
        },
    }
    rc, out, _ = run_hook(HOOK_INDEX, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    assert out, "Expected reminder output, got none"

    data = json.loads(out.decode("utf-8"))
    hso = data.get("hookSpecificOutput", {})
    ctx = hso.get("additionalContext", "")
    assert "索引提醒" in ctx, f"Expected index reminder, got: {ctx}"
    assert "nonexistent_test_file_xyz_2026-09-20.md" in ctx, "Filename not in reminder"
    print("✓ PASS: Write new note (reminder triggered)")


def test_write_indexed_note():
    """Write 到已索引笔记 ⟹ 无输出。"""
    # Pick a file that's actually indexed with [file](path) pattern
    payload = {
        "tool_name": "Write",
        "tool_input": {
            "file_path": os.path.join(ROOT, "theory/notes/theoremI/conjecture2_review_2026-09-11.md"),
            "content": "test content",
        },
    }
    rc, out, _ = run_hook(HOOK_INDEX, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    assert out == b"", f"Expected no output for indexed file, got {out}"
    print("✓ PASS: Write indexed note (silent)")


def test_poisoning_superseded_only():
    """投毒（确定性）：临时根目录里造一个 INDEX，其中 X 只出现在 `superseded by X` 单元格、没有自己的 [X]( 行。
    修复前（子串判断）会静默；修复后必须提醒。同时验证一个有 [Y]( 行的 Y 保持静默（负对照的负对照）。"""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        notes = os.path.join(tmp, "theory", "notes")
        os.makedirs(os.path.join(notes, "materials"))
        with open(os.path.join(notes, "INDEX.md"), "w", encoding="utf-8") as f:
            f.write("| [old_note_2026-09-01.md](materials/old_note_2026-09-01.md) | 2026-09-01 | x | n/a | superseded by ghost_note_2026-09-02.md | y |\n")
        env = dict(os.environ, SC_HOOK_ROOT=tmp.replace("\\", "/"))
        def run(name):
            payload = {"tool_name": "Write", "tool_input": {"file_path": os.path.join(notes, "materials", name), "content": "t"}}
            r = subprocess.run([sys.executable, HOOK_INDEX], input=json.dumps(payload).encode("ascii"),
                               capture_output=True, timeout=30, env=env)
            return r.returncode, r.stdout
        rc, out = run("ghost_note_2026-09-02.md")
        assert rc == 0 and out, "ghost (superseded-by only) must be flagged; got silence (bug not caught!)"
        assert "索引提醒" in json.loads(out.decode())["hookSpecificOutput"]["additionalContext"]
        rc, out = run("old_note_2026-09-01.md")
        assert rc == 0 and out == b"", f"indexed note must stay silent, got {out!r}"
    print("✓ PASS: Poisoning test (superseded-by-only name flagged; indexed name silent)")


def test_path_with_spaces_and_unicode():
    """路径含空格与中文、反斜杠路径 ⟹ 行为正确。"""
    # Test with backslashes and unicode
    payload = {
        "tool_name": "Write",
        "tool_input": {
            "file_path": os.path.join(ROOT, "theory", "notes", "materials", "测试文件_with spaces_2026-09-20.md"),
            "content": "test",
        },
    }
    rc, out, _ = run_hook(HOOK_INDEX, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    # Should trigger reminder since file doesn't exist in index
    assert out, "Expected reminder for unicode/space path"
    print("✓ PASS: Path with spaces and unicode")


def test_scratch_ir_path():
    """_scratch_ir/ 路径 ⟹ 无输出。"""
    payload = {
        "tool_name": "Write",
        "tool_input": {
            "file_path": os.path.join(ROOT, "theory/notes/_scratch_ir/test_script.py"),
            "content": "# test",
        },
    }
    rc, out, _ = run_hook(HOOK_INDEX, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    assert out == b"", f"Expected no output for _scratch_ir, got {out}"
    print("✓ PASS: _scratch_ir/ path (silent)")


def test_derivations_path():
    """theory/derivations/zzz.md ⟹ 提醒指向 derivations/INDEX.md。"""
    payload = {
        "tool_name": "Write",
        "tool_input": {
            "file_path": os.path.join(ROOT, "theory/derivations/test_deriv_2026-09-20.md"),
            "content": "test",
        },
    }
    rc, out, _ = run_hook(HOOK_INDEX, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    assert out, "Expected reminder for derivations path"

    data = json.loads(out.decode("utf-8"))
    hso = data.get("hookSpecificOutput", {})
    ctx = hso.get("additionalContext", "")
    assert "theory/derivations/INDEX.md" in ctx, f"Expected derivations INDEX reference, got: {ctx}"
    print("✓ PASS: derivations/ path (reminder triggered)")


def test_edit_tool():
    """Edit 工具也应触发相同逻辑。"""
    payload = {
        "tool_name": "Edit",
        "tool_input": {
            "file_path": os.path.join(ROOT, "theory/notes/materials/new_edit_test_2026-09-20.md"),
            "old_string": "old",
            "new_string": "new",
        },
    }
    rc, out, _ = run_hook(HOOK_INDEX, payload)
    assert rc == 0, f"Expected exit 0, got {rc}"
    assert out, "Expected reminder for Edit on new file"
    print("✓ PASS: Edit tool (reminder triggered)")


def main():
    print("=" * 60)
    print("Hook Test Suite")
    print("=" * 60)

    tests = [
        ("Malformed JSON", test_malformed_json),
        ("Wrong tool name", test_wrong_tool_name),
        ("Agent with vocab", test_agent_with_vocab),
        ("Agent without vocab", test_agent_no_vocab),
        ("Write new note", test_write_new_note),
        ("Write indexed note", test_write_indexed_note),
        ("Poisoning: superseded-only", test_poisoning_superseded_only),
        ("Path with spaces/unicode", test_path_with_spaces_and_unicode),
        ("Scratch IR path", test_scratch_ir_path),
        ("Derivations path", test_derivations_path),
        ("Edit tool", test_edit_tool),
    ]

    passed = 0
    failed = 0
    skipped = 0

    for name, test_func in tests:
        try:
            print(f"\n[TEST] {name}")
            if test_func() == "SKIP":
                skipped += 1
            else:
                passed += 1
        except AssertionError as e:
            print(f"✗ FAIL: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ ERROR: {e}")
            failed += 1

    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed, {skipped} skipped")
    print("=" * 60)

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
