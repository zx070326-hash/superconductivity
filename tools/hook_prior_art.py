#!/usr/bin/env python
"""hook_prior_art.py — PreToolUse(Agent) hook：派活前自动跑先验检索。

行为（2026-09-20 本机实测，见 tools/hook_probe.log 已删的三次探针）：
  · additionalContext → 只有**派活方（lead）**看得到；子 agent 看不到。
  · updatedInput.prompt → 子 agent 收到的就是改写后的 prompt（实测首行原样可见）。
  ⟹ 两边都给：lead 看到检索结果以核对 D5；执行方在 prompt 末尾看到同一块（验收门的自动化部分）。
纪律：
  · 只**追加**，绝不改动派活方写的正文；追加块显式标明来源，执行方知道那不是派活方的判断。
  · 只报信息，不阻断，不判决（prior_art.sh 的约定）。任何异常 ⟹ 静默放行（exit 0、无输出）。
  · 输出 ensure_ascii，避免 Windows 控制台编码把 JSON 打坏。
"""
import json
import os
import subprocess
import sys

# Derive ROOT from script location; SC_HOOK_ROOT env var overrides (for testing)
ROOT = os.environ.get('SC_HOOK_ROOT') or os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
).replace("\\", "/")


def main() -> None:
    raw = sys.stdin.buffer.read()
    try:
        d = json.loads(raw.decode("utf-8"))
    except Exception:
        return
    if d.get("tool_name") != "Agent":
        return
    ti = d.get("tool_input") or {}
    text = (ti.get("description") or "") + "\n" + (ti.get("prompt") or "")
    try:
        r = subprocess.run(
            ["sh", ROOT + "/tools/prior_art.sh"],
            input=text.encode("utf-8"),
            capture_output=True,
            timeout=25,
        )
        out = r.stdout.decode("utf-8", "replace").strip()
    except Exception:
        return
    if not out:
        return  # 没激活任何承重量 ⟹ 零噪声
    new_ti = dict(ti)
    new_ti["prompt"] = (
        (ti.get("prompt") or "")
        + "\n\n---\n【自动附加 · 先验检索 — 由 hook 生成，不是派活方的判断；读这些文件的当前版本，不信摘要】\n"
        + out
    )
    payload = {
        "hookSpecificOutput": {
            "hookEventName": "PreToolUse",
            "permissionDecision": "allow",
            "additionalContext": out,
            "updatedInput": new_ti,
        }
    }
    sys.stdout.buffer.write(json.dumps(payload, ensure_ascii=True).encode("ascii"))


if __name__ == "__main__":
    try:
        main()
    except Exception:
        pass
    sys.exit(0)
