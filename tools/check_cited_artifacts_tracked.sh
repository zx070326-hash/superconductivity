#!/bin/sh
# 毕业规则自查（非 .py 版），2026-09-25。
#
# 消费点：CLAUDE.md §2「脚本毕业规则」的理由——「入库的笔记不该指向任何 clone 里
# 都不存在的文件」——与扩展名无关。既有的 tools/check_cited_scripts_tracked.sh
# 只解析 `*.py`；实测（2026-09-25）被 .md 点名的非 .py 产物（.png/.txt/.json/.csv/...）
# 有 46 个，其中 21 个未跟踪、4 个盘上不存在，而 .py 门在这一整天报 PASS。
# 本脚本补上这一半，**不改** check_cited_scripts_tracked.sh（它硬编码的基线不能动）。
#
# 覆盖的扩展名：.png .txt .json .npy .csv .md（不含 .py——那是另一个门管的）。
# 正则允许连字符/点（09-19 事故：`[A-Za-z0-9_]+` 把带日期的文件名截断）。
# 引用前缀允许 `_scratch_ir/`、`../_scratch_ir/`、`theory/notes/_scratch_ir/`，
# 一律归一为仓库相对路径 `theory/notes/_scratch_ir/<rest>`（<rest> 可含子目录，
# 已核实 _scratch_ir 下确有子目录，故字符集含 `/`）。
#
# 三种状态：
#   (a) MISSING     — 盘上不存在
#   (b) UNTRACKED   — 存在但未被 git 跟踪，且不在允许清单里
#   (c) OK          — 已跟踪，或在允许清单里的 UNTRACKED
#   (d) LOST-REG    — 盘上不存在，但**双钥匙**齐全：路径在 tools/lost_artifacts_registry.txt，
#                     且引用它的那篇 .md 自身带「产物缺失登记」标记。缺任一钥匙仍判 MISSING。
#
# 允许清单：tools/cited_artifacts_allowlist.txt，每行 `<仓库相对路径><TAB><理由>`。
# 用途仅限两类：版权 PDF（*.pdf 是有意的 gitignore）、受版权论文的全文文本抽取层。
# 清单里的条目仍会被打印（标 ALLOWLISTED），但不计入失败。
#
# 用法： sh tools/check_cited_artifacts_tracked.sh   （退出码非 0 = 有 (a) 或未清单化的 (b)）
set -e
cd "$(git rev-parse --show-toplevel)"

allowlist="tools/cited_artifacts_allowlist.txt"
lostreg="tools/lost_artifacts_registry.txt"
mds=$(mktemp)
raw=$(mktemp)
norm=$(mktemp)
result=$(mktemp)
trap 'rm -f "$mds" "$raw" "$norm" "$result"' EXIT

find theory/notes theory/derivations -name '*.md' -not -path 'theory/notes/_scratch_ir/*' -print0 > "$mds"
n_md=$(tr -cd '\0' < "$mds" | wc -c)

# 抽取：形如 (前缀)_scratch_ir/(路径).(ext)，附带来源文件:行号。
# 用 grep -n 逐文件跑以拿到行号（xargs -0 + -n 需要逐文件，故用 while 循环）。
: > "$raw"
find theory/notes theory/derivations -name '*.md' -not -path 'theory/notes/_scratch_ir/*' -print0 \
  | xargs -0 -I{} sh -c '
      grep -noE "_scratch_ir/[A-Za-z0-9_./-]+\.(png|txt|json|npy|csv|md)" "{}" \
        | while IFS=: read -r lineno match; do
            printf "%s\t%s\t%s\n" "{}" "$lineno" "$match"
          done
    ' >> "$raw" 2>/dev/null || true

n_cited_raw=$(wc -l < "$raw" | tr -d ' ')

# 归一化：把 match（可能带 ../ 或 theory/notes/ 前缀，已被上面正则天然剥掉，因为
# 正则只从 `_scratch_ir/` 开始匹配）映射到仓库相对路径。
: > "$norm"
while IFS="$(printf '\t')" read -r srcfile lineno match; do
  rel="theory/notes/${match}"
  printf "%s\t%s\t%s\n" "$rel" "$srcfile" "$lineno" >> "$norm"
done < "$raw"

n_paths=$(cut -f1 "$norm" | sort -u | wc -l | tr -d ' ')

: > "$result"
n_missing=0
n_untracked_flagged=0
n_untracked_allowlisted=0
n_ok=0
sort -u -t"$(printf '\t')" -k1,1 "$norm" | cut -f1 | sort -u | while read -r rel; do
  # 该路径的一条引用来源（供报告；同一路径可能被多处引用，取第一条）
  srcref=$(awk -F'\t' -v r="$rel" '$1==r{print $2":"$3; exit}' "$norm")
  if [ ! -e "$rel" ]; then
    # (d) 已登记缺失：**双钥匙**，缺一不可——
    #     钥匙① 路径在 tools/lost_artifacts_registry.txt 里；
    #     钥匙② 引用它的那篇 .md 自身带「产物缺失登记」标记（即损失是公开写下的，不是被清单悄悄藏掉的）。
    #     单靠钥匙① 无法让一个真缺失静默通过——这正是 CLAUDE.md「记录丢失时不得反向重建」
    #     需要的出口：不重建、不删引用、但也不让门永久 FAIL。
    lostreason=""
    if [ -f "$lostreg" ]; then
      lostreason=$(awk -F'\t' -v r="$rel" '$1==r{print $2; exit}' "$lostreg")
    fi
    if [ -n "$lostreason" ]; then
      srcfile=$(awk -F'\t' -v r="$rel" '$1==r{print $2; exit}' "$norm")
      if grep -q '产物缺失登记' "$srcfile" 2>/dev/null; then
        echo "LOST-REG    $rel | 登记: $lostreason | 引用: $srcref" >> "$result"
        continue
      fi
      echo "MISSING     $rel | ⚠ 在缺失登记册里但引用处 $srcfile 无「产物缺失登记」标记（钥匙②缺失） | 引用: $srcref" >> "$result"
      continue
    fi
    echo "MISSING     $rel | 引用: $srcref" >> "$result"
    continue
  fi
  if git ls-files --error-unmatch "$rel" >/dev/null 2>&1; then
    echo "OK          $rel | 引用: $srcref" >> "$result"
    continue
  fi
  reason=""
  if [ -f "$allowlist" ]; then
    reason=$(awk -F'\t' -v r="$rel" '$1==r{print $2; exit}' "$allowlist")
  fi
  if [ -n "$reason" ]; then
    echo "ALLOWLISTED $rel | 理由: $reason | 引用: $srcref" >> "$result"
  else
    echo "UNTRACKED   $rel | 引用: $srcref" >> "$result"
  fi
done

cat "$result"

n_missing=$(grep -c '^MISSING' "$result" || true)
n_untracked=$(grep -c '^UNTRACKED   ' "$result" || true)
n_allowlisted=$(grep -c '^ALLOWLISTED' "$result" || true)
n_lostreg=$(grep -c '^LOST-REG' "$result" || true)
n_ok=$(grep -c '^OK ' "$result" || true)

echo "COVERAGE: 扫描 $n_md 篇 .md，解析出 $n_cited_raw 条原始引用、去重后 $n_paths 个被引产物路径。"
echo "COUNTS: OK=$n_ok  ALLOWLISTED=$n_allowlisted  UNTRACKED(未清单化)=$n_untracked  LOST-REG=$n_lostreg  MISSING=$n_missing"

if [ "${n_md:-0}" -lt 100 ]; then
  echo "FAIL: 覆盖面异常偏低（笔记 $n_md 篇）——检测器多半没看到输入，PASS 不可信。"
  exit 1
fi

if [ "${n_missing:-0}" -gt 0 ] || [ "${n_untracked:-0}" -gt 0 ]; then
  echo "FAIL: 存在盘上不存在或未跟踪且未清单化的被引产物（见上）。"
  echo "  盘上不存在 ⟹ 检查引用是否打错字，或该产物是否本该重新生成/入库。"
  echo "  未跟踪 ⟹ git add -f <path>，或若确属版权豁免，加入 $allowlist。"
  echo "  盘上不存在且确属不可重建 ⟹ 按双钥匙登记：写入 $lostreg，**并**在引用处写下「产物缺失登记」。"
  exit 1
fi

echo "PASS: 所有被笔记引用的非 .py 产物均已跟踪、已清单化，或已按双钥匙登记为缺失。"
