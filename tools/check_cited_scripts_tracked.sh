#!/bin/sh
# 毕业规则自查（2026-09-18）：凡被 notes/ 或 derivations/ 的 .md 点名引用为证据的脚本，
# 必须是 git 跟踪的；否则入库的笔记会引用一个在 clone 里根本不存在的文件。
# 用法： sh tools/check_cited_scripts_tracked.sh      （退出码 1 = 有漏网）
#
# 2026-09-19 修订（两处假 PASS）：
#   ① 原 `theory/notes/*.md` 不递归——09-18 笔记全部迁入 12 个子目录后，gate 实际扫不到任何笔记；
#   ② 原正则 `[A-Za-z0-9_]+\.py` 不含 `-`——按命名规则带日期的脚本（`*_2026-09-19.py`）被截成 `19.py` 而漏网。
#   现改为递归扫描（排除 _scratch_ir 自身的 .md 草稿日志），文件名字符集加入 `-` 与 `.`。
#
# 2026-09-20 追加：**打印覆盖面**。上面两处缺陷都不是「检测失灵」而是「输入没到检测器」——
#   投毒负对照只能证明 `P(报警 | 输入到达检测器)`，证明不了 `P(输入到达检测器)`。
#   实测：修复前该 gate 只扫到 1/184 篇笔记、28/434 个可解析引用（覆盖率 6.5%），
#   而它在那 6.5% 上的检测率是 100% —— 于是 PASS 一整天。
#   ⟹ 不打印覆盖面的 PASS 是不可证伪的。任何人改这个脚本，请保留下面两行计数。
set -e
cd "$(git rev-parse --show-toplevel)"
tmp=$(mktemp)
mds=$(mktemp)
find theory/notes theory/derivations -name '*.md' -not -path 'theory/notes/_scratch_ir/*' -print0 > "$mds"
n_md=$(tr -cd '\0' < "$mds" | wc -c)
n_cited=$(xargs -0 grep -hoE '[A-Za-z0-9_.-]+\.py' < "$mds" 2>/dev/null | sort -u \
  | while read -r f; do [ -f "theory/notes/_scratch_ir/$f" ] && echo "$f"; done | wc -l)
xargs -0 grep -hoE '[A-Za-z0-9_.-]+\.py' < "$mds" 2>/dev/null \
  | sort -u \
  | while read -r f; do
      p="theory/notes/_scratch_ir/$f"
      [ -f "$p" ] || continue
      if ! git ls-files --error-unmatch "$p" >/dev/null 2>&1; then
        echo "UNTRACKED-BUT-CITED  $p"
      fi
    done > "$tmp" || true
n=$(grep -c . "$tmp" || true)
cat "$tmp"
rm -f "$tmp" "$mds"
echo "COVERAGE: 扫描 $n_md 篇 .md，解析出 $n_cited 个落在 _scratch_ir 的被引脚本。"
if [ "${n_md:-0}" -lt 100 ] || [ "${n_cited:-0}" -lt 100 ]; then
  echo "FAIL: 覆盖面异常偏低（笔记 $n_md / 被引脚本 $n_cited）——检测器多半没看到输入，PASS 不可信。"
  exit 1
fi
if [ "${n:-0}" -gt 0 ]; then
  echo "FAIL: $n 个被引用脚本未入库。对每个跑： git add -f <path>"
  exit 1
fi
echo "PASS: 所有被笔记引用的 _scratch_ir 脚本均已跟踪。"
