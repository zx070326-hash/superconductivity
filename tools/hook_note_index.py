#!/usr/bin/env python
"""hook_note_index.py — PostToolUse(Write|Edit) hook：写了笔记却没补索引时提醒。

触发：theory/notes/**/*.md 或 theory/derivations/**/*.md 被 Write/Edit，且文件名不在对应 INDEX.md 里。
行为：只提醒（additionalContext 给当前会话），不阻断、不改文件。
事故背书：2026-09-20 一天 3 篇未索引笔记 + deficit-real-solve 两篇（其一已被自己的后一次交付取代）。
任何异常 ⟹ 静默放行。
"""
import json
import os
import sys

# Derive ROOT from script location; SC_HOOK_ROOT env var overrides (for testing)
ROOT = os.environ.get('SC_HOOK_ROOT') or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
).replace("\\", "/")
INDEXES = {
    "theory/notes/": ROOT + "/theory/notes/INDEX.md",
    "theory/derivations/": ROOT + "/theory/derivations/INDEX.md",
}


def main() -> None:
    raw = sys.stdin.buffer.read()
    try:
        d = json.loads(raw.decode("utf-8"))
    except Exception:
        return
    if d.get("tool_name") not in ("Write", "Edit"):
        return
    fp = (d.get("tool_input") or {}).get("file_path") or ""
    fp = fp.replace("\\", "/")
    if not fp.endswith(".md"):
        return
    base = os.path.basename(fp)
    if base in ("INDEX.md", "CLAUDE.md"):
        return
    for prefix, index in INDEXES.items():
        if ("/" + prefix) in fp or fp.startswith(prefix):
            if "_scratch_ir" in fp:
                return
            try:
                with open(index, encoding="utf-8") as f:
                    idx = f.read()
            except Exception:
                return
            # Check for link pattern [base]( to avoid false positives from "superseded by base"
            if f"[{base}](" in idx:
                return
            msg = (
                f"【索引提醒】{base} 不在 {os.path.relpath(index, ROOT).replace(os.sep, '/')} 里。"
                "路由规则要求落盘同一次编辑里补索引行（含状态列；若取代旧笔记，旧行写 superseded by <文件名>）。"
            )
            payload = {
                "hookSpecificOutput": {
                    "hookEventName": "PostToolUse",
                    "additionalContext": msg,
                }
            }
            sys.stdout.buffer.write(json.dumps(payload, ensure_ascii=True).encode("ascii"))
            return


if __name__ == "__main__":
    try:
        main()
    except Exception:
        pass
    sys.exit(0)
