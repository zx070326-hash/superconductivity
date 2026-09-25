#!/bin/sh
# INDEX 结构自查（2026-09-25 新建）。
# 用法： sh tools/check_index_structure.sh            （退出码 1 = 有结构坏）
#        sh tools/check_index_structure.sh --poison   （投毒负对照：确认 FAIL 路径真的会报警）
#
# 缘起：`theory/notes/INDEX.md` 的**表格结构**从来没有任何门在查——只有「被引脚本是否入库」
#   这类**内容**门。2026-09-25 一次手工全表扫描抓到五类坏，**全部是静默的**：
#     ① 少一列（6 列写成 5 列）——markdown 照样渲染，只是右移一格，肉眼看不出；
#     ② 两条记录合并在同一物理行上——后一条**整行不可见**，且链接检查照样 PASS（链接本身是好的）；
#     ③ 单元格内有未转义的 `|`（`|D|/|Λ|`、`max|Δ|=0`）——凭空多出 2–7 列。
#        ⚠ 反引号**不保护** `|`：表格切分发生在 markdown 解析**之前**，代码段里的 `|` 一样断列，必须写 `\|`；
#     ④ 重复行——同一篇笔记被两个会话各写一行，**列数全对、链接全活、覆盖率 100%**。
#        危害不是冗余而是**判词分叉**：两行的`严格性`/`状态`列内容不同 ⟹ 读者落在哪一行决定他信哪个状态；
#     ⑤ 目录总览表的计数陈旧（写 126，实测 250）。
#   ①②③④ 由本脚本覆盖；⑤ 不可自动化（要人判「陈旧」），只打印实测值供比对。
#
# ⚠ 覆盖面自述（沿用 check_cited_scripts_tracked.sh 2026-09-20 的教训）：
#   投毒负对照只能证明 `P(报警 | 输入到达检测器)`，证明不了 `P(输入到达检测器)`。
#   所以下面每一项都**打印它数到了多少行**，不打印计数的 PASS 是不可证伪的。
#
# ⚠ 查重必须区分**行锚**与**单元格内引用**：本脚本只取**每行第一个** `](*.md)` 作为行锚。
#   我的第一版查重器 grep 整行**所有**链接，于是把「依据列里引用了另一篇笔记」也算成重复，12 个里 3 个是假阳性。

cd "$(git rev-parse --show-toplevel)"

# ⚠ 本机 `python3` 是 Windows Store 的占位程序（静默退出 49）——挑一个真解释器。
PY_BIN=""
for c in python python3 py; do
  if "$c" -c "import sys" >/dev/null 2>&1; then PY_BIN="$c"; break; fi
done
if [ -z "$PY_BIN" ]; then echo "找不到可用的 python 解释器"; exit 2; fi

TARGET="theory/notes/INDEX.md"
POISON=0
[ "$1" = "--poison" ] && POISON=1

if [ "$POISON" = "1" ]; then
  TARGET=$(mktemp)
  # 投毒：四类坏各造一条，且都用**当前命名结构**（带子目录、带连字符日期名）——
  # 2026-09-19 的假 PASS 正是因为投毒样本用了旧结构，恰好绕开缺陷。
  {
    printf '## poison\n\n'
    printf '| [a_2026-09-25.md](misc/a_2026-09-25.md) | 2026-09-25 | ok | analytic | live | x |\n'
    printf '| [b_2026-09-25.md](misc/b_2026-09-25.md) | 2026-09-25 | 少一列 | live | x |\n'
    printf '| [c_2026-09-25.md](misc/c_2026-09-25.md) | 2026-09-25 | 未转义 `|D|/|Λ|` | analytic | live | x |\n'
    printf '| [a_2026-09-25.md](misc/a_2026-09-25.md) | 2026-09-25 | 重复行 | analytic | live | x |\n'
    printf '| [gone_2026-09-25.md](misc/gone_2026-09-25.md) | 2026-09-25 | 死链 | analytic | live | x |\n'
  } > "$TARGET"
  echo "POISON MODE: 用一份人造 INDEX（4 类坏 + 1 死链）跑检测器，期望 FAIL。"
fi

OUT=$(mktemp)
INDEX="$TARGET" PYTHONIOENCODING=utf-8 "$PY_BIN" - <<'PY' 2>&1 | tee "$OUT"
import os, re, sys, collections
# Windows 控制台默认 GBK，打印 `⟹` 会 UnicodeEncodeError 把检测器打崩。
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
p = os.environ["INDEX"]
txt = open(p, encoding="utf-8").read()
lines = txt.split("\n")
BS = chr(92)

seg_dist = collections.Counter()
bad_cols, anchors = [], {}
for i, l in enumerate(lines, 1):
    s = l.strip()
    if not (s.startswith("|") and s.endswith("|")):
        continue
    # 先把 \| 抠掉再数，否则转义过的管道会被当成列分隔符
    n = len(s.replace(BS + "|", "\x00").split("|"))
    seg_dist[n] += 1
    if n not in (5, 8):          # 5 = 目录总览表(3 列)，8 = 笔记行(6 列)
        bad_cols.append((i, n, s[:110]))
    m = re.search(r"\]\(([^)]+\.md)\)", s)      # 只取第一个链接 = 行锚
    if m:
        anchors.setdefault(m.group(1), []).append(i)

dups = {k: v for k, v in anchors.items() if len(v) > 1}

links = re.findall(r"\]\(([^)]+\.md)\)", txt)
base = os.path.dirname(p) if os.path.dirname(p) else "."
if os.path.abspath(p) == os.path.abspath("theory/notes/INDEX.md"):
    base = "theory/notes"
dangling = [l for l in links if not os.path.exists(os.path.join(base, l))]

disk = set()
if os.path.isdir("theory/notes"):
    for r, ds, fs in os.walk("theory/notes"):
        ds[:] = [d for d in ds if d != "_scratch_ir"]
        for f in fs:
            if f.endswith(".md") and f not in ("INDEX.md", "CLAUDE.md"):
                disk.add(os.path.relpath(os.path.join(r, f), "theory/notes").replace(os.sep, "/"))
uncovered = sorted(disk - set(links)) if base == "theory/notes" else []

fail = 0
print("[COVERAGE] 本次检测器看到 %d 行表格行、%d 个行锚、%d 个链接、磁盘 %d 篇笔记。"
      % (sum(seg_dist.values()), sum(len(v) for v in anchors.values()), len(links), len(disk)))

print("\n[1/4 列数] 分段数分布: %s" % dict(sorted(seg_dist.items())))
if bad_cols:
    fail = 1
    print("  FAIL: %d 行列数异常（应为 5 或 8 段）" % len(bad_cols))
    for i, n, s in bad_cols[:10]:
        print("    L%-5d %d 段  %s" % (i, n, s))
else:
    print("  PASS")

print("\n[2/4 重复行锚]")
if dups:
    fail = 1
    print("  FAIL: %d 个锚出现多次（同一篇笔记两行 ⟹ 状态列可能分叉）" % len(dups))
    for k in sorted(dups):
        print("    x%d  L%s  %s" % (len(dups[k]), dups[k], k))
else:
    print("  PASS（%d 个行锚全部唯一）" % len(anchors))

print("\n[3/4 死链]")
if dangling:
    fail = 1
    print("  FAIL: %d 条链接指向不存在的 .md" % len(dangling))
    for d in dangling[:10]:
        print("    %s" % d)
else:
    print("  PASS（%d 条链接全部命中）" % len(links))

print("\n[4/4 覆盖]")
if uncovered:
    fail = 1
    print("  FAIL: %d 篇磁盘上的笔记未被索引收录" % len(uncovered))
    for u in uncovered[:10]:
        print("    %s" % u)
else:
    print("  PASS（磁盘 %d 篇，未收录 0 篇）" % len(disk))

print("\n[参考·不判定] 各子目录实测笔记数（目录总览表的计数请人工比对；陈旧计数无法自动判）：")
cnt = collections.Counter(f.split("/")[0] for f in disk)
print("    " + "  ".join("%s=%d" % (k, v) for k, v in sorted(cnt.items())) + "  合计=%d" % len(disk))

sys.exit(1 if fail else 0)
PY
rc=${PIPESTATUS[0]:-$?}

if [ "$POISON" = "1" ]; then
  rm -f "$TARGET"
  echo ""
  if [ "$rc" = "0" ]; then
    echo "POISON FAILED: 检测器对投毒样本给了 PASS —— 这个门是坏的。"
    exit 1
  fi
  # ⚠ 不能只看「非零」：检测器**崩溃**也是非零。实测踩过一次（GBK UnicodeEncodeError
  #   把检测器打崩，退出非零，于是被误判成「正确报 FAIL」）。要求退出码**恰为 1**。
  if [ "$rc" != "1" ]; then
    echo "POISON FAILED: 退出码 $rc != 1 —— 检测器多半是**崩了**而不是报 FAIL。"
    echo "               （非零 != 报警；崩溃与判 FAIL 必须分开。）"
    exit 1
  fi
  # ⚠ 「四项判词齐全」以前是我**印出来但没查**的一句话 —— 崩溃时前 1 项就停了，
  #   那句话照印。既然要断言，就得真查：四个检查头必须都出现，且四项都要有判词。
  miss=""
  for tag in "1/4 列数" "2/4 重复行锚" "3/4 死链" "4/4 覆盖"; do
    grep -q "$tag" "$OUT" || miss="$miss [$tag]"
  done
  n_verdict=$(grep -c -E "^  (PASS|  FAIL|FAIL)" "$OUT")
  rm -f "$OUT"
  if [ -n "$miss" ]; then
    echo "POISON FAILED: 检测器中途停了，缺判词:$miss"
    echo "               （退出码虽为 1，但它没跑完四项 ⟹ 未被跑到的检查未受投毒检验。）"
    exit 1
  fi
  if [ "$n_verdict" -lt 4 ]; then
    echo "POISON FAILED: 只数到 $n_verdict 条 PASS/FAIL 判词，应为 4 条。"
    exit 1
  fi
  echo "POISON OK: 退出码 1，四个检查头齐全，$n_verdict 条判词（已逐项核，非自述）。"
  echo "           覆盖面自述：投毒样本 4 类坏 + 1 死链，检测器 1/2/3 项均报 FAIL。"
  exit 0
fi

exit $rc
