#!/bin/sh
# 毕业规则自查（2026-09-18）：凡被 notes/ 或 derivations/ 的 .md 点名引用为证据的脚本，
# 必须是 git 跟踪的；否则入库的笔记会引用一个在 clone 里根本不存在的文件。
# 用法： sh tools/check_cited_scripts_tracked.sh      （退出码 1 = 有漏网）
set -e
cd "$(git rev-parse --show-toplevel)"
miss=0
grep -rhoE '[A-Za-z0-9_]+\.py' theory/notes/*.md theory/derivations/*.md 2>/dev/null \
  | sort -u \
  | while read -r f; do
      p="theory/notes/_scratch_ir/$f"
      [ -f "$p" ] || continue
      if ! git ls-files --error-unmatch "$p" >/dev/null 2>&1; then
        echo "UNTRACKED-BUT-CITED  $p"
        miss=1
      fi
    done > /tmp/_cited_check.$$ 2>/dev/null || true
n=$(wc -l < /tmp/_cited_check.$$ 2>/dev/null || echo 0)
cat /tmp/_cited_check.$$ 2>/dev/null || true
rm -f /tmp/_cited_check.$$
if [ "${n:-0}" -gt 0 ]; then
  echo "FAIL: $n 个被引用脚本未入库。对每个跑： git add -f <path>"
  exit 1
fi
echo "PASS: 所有被笔记引用的 _scratch_ir 脚本均已跟踪。"
