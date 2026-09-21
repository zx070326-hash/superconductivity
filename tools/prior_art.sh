#!/bin/sh
# prior_art.sh — 先验检索：这个量／材料，仓库里是不是已经算过？是不是已被撤回？
#
# 用法：
#   手工：  sh tools/prior_art.sh "𝔣 的 λ 曲率 q 在 H3S 上"
#   管道：  echo "..." | sh tools/prior_art.sh
#   hook：  由 hook_prior_art.sh 调用（它负责 JSON 化）
#
# 设计约束（每条都有事故背书）：
#  · 只报信息，不做判决，绝不阻断。重算常常是正当的独立二次验证——
#    今天的事故不是「我决定重算」，是「我不知道它已经存在」。
#  · 状态复用 theory/notes/INDEX.md 已有的状态列，不另立账本（避免两处真源漂移 + 台账腐烂）。
#  · 状态解析不依赖列位置——INDEX.md 的单元格里含公式管道符，按列取会在承重行上错位（已实测）。
#    改为行内关键词匹配，按严重度优先：retracted > superseded > exploratory > unclear > live。
#    方向是保守的：宁可过度报警，不可漏报（单边界方向规则）。
#  · 匹配走 grep -F 固定串。**不要**把多字节字符放进正则字符类 `[...]`，会按字节拆成乱码（已实测）。
#  · 输出必须短。它进的是主会话上下文，长了就是污染。

set -u
# Derive ROOT from script location
ROOT=$(cd "$(dirname "$0")/.." && pwd)
NOTES="$ROOT/theory/notes"
VOCAB="$ROOT/tools/vocab.tsv"
MAXLIVE=6        # live 命中最多印几条
MAXSTALE=5       # 非 live 命中最多印几条
TAB=$(printf '\t')

# 输入：参数优先，否则读 stdin
if [ $# -gt 0 ]; then IN="$*"; else IN=$(cat); fi
[ -n "$IN" ] || exit 0
[ -f "$VOCAB" ] || exit 0

# ---- 1. 哪些承重量被这段文本激活了 ----
ACTIVE=""      # "slug|显示名" 每行一条
PATSALL=""     # 所有激活的匹配串，供后续 grep
while IFS="$TAB" read -r slug disp pats; do
    case "$slug" in \#*|"") continue ;; esac
    [ -n "${pats:-}" ] || continue
    hit=0
    OIFS=$IFS; IFS='|'
    for p in $pats; do
        [ -n "$p" ] || continue
        if printf '%s' "$IN" | grep -qF -- "$p" 2>/dev/null; then hit=1; break; fi
    done
    IFS=$OIFS
    if [ "$hit" -eq 1 ]; then
        ACTIVE="$ACTIVE$slug$TAB$disp
"
        PATSALL="$PATSALL$pats|"
    fi
done < "$VOCAB"

[ -n "$ACTIVE" ] || exit 0   # 没激活任何承重量 ⟹ 安静退出，零噪声

# ---- 2. 文件 -> 状态（鲁棒解析，不依赖列位置）----
STATUSMAP=$(awk '/\.md\]\(/ {
    if (match($0, /\[[^]]+\.md\]/)) fn=substr($0,RSTART+1,RLENGTH-2); else next;
    st="?";
    if ($0 ~ /retracted/)         st="retracted";
    else if ($0 ~ /superseded/)   st="superseded";
    else if ($0 ~ /exploratory/)  st="exploratory";
    else if ($0 ~ /unclear/)      st="unclear";
    else if ($0 ~ /live/)         st="live";
    printf "%s\t%s\n", fn, st;
}' "$NOTES/INDEX.md" 2>/dev/null)

lookup_status() {   # $1 = 裸文件名
    printf '%s\n' "$STATUSMAP" | awk -F"$TAB" -v f="$1" '$1==f {print $2; found=1; exit} END{if(!found) print "未索引"}'
}

# ---- 3. 检索命中文件 ----
# 判据不是「命中任一词」——那会让 H3S 这种宽词淹没结果（实测 70 篇）。
# 判据是**同时命中多少个不同的激活承重量**：命中数越高越相关。
# 逐 slug 取并集（同一 slug 的多个别名只算一次），再按 slug 命中数排序。
HITS=""
while IFS="$TAB" read -r slug disp; do
    [ -n "$slug" ] || continue
    sfiles=""
    pats=$(awk -F"$TAB" -v s="$slug" '$1==s {print $3}' "$VOCAB")
    OIFS=$IFS; IFS='|'
    for p in $pats; do
        [ -n "$p" ] || continue
        f=$(grep -rlF --include='*.md' --exclude-dir=_scratch_ir --exclude-dir=probes -- "$p" "$NOTES" 2>/dev/null)
        [ -n "$f" ] && sfiles="$sfiles$f
"
    done
    IFS=$OIFS
    # 同一 slug 内去重，再累加到总表（于是每个文件每个 slug 最多计 1）
    sfiles=$(printf '%s\n' "$sfiles" | grep -v '^$' | sort -u)
    [ -n "$sfiles" ] && HITS="$HITS$sfiles
"
done <<EOF
$(printf '%s' "$ACTIVE" | grep -v '^$')
EOF

NACTIVE=$(printf '%s' "$ACTIVE" | grep -c .)
# RANKED: "<命中slug数> <路径>"，按命中数降序、再按文件名尾部日期降序
RANKED=$(printf '%s\n' "$HITS" | grep -v '^$' | grep -v '/INDEX\.md$' \
    | sort | uniq -c | sort -rn \
    | awk '{n=$1; $1=""; sub(/^ /,""); path=$0;
            if (match(path, /[0-9]{4}-[0-9]{2}-[0-9]{2}/)) d=substr(path,RSTART,RLENGTH); else d="0000-00-00";
            printf "%d\t%s\t%s\n", n, d, path}' \
    | sort -t"$TAB" -k1,1nr -k2,2r)
[ -n "$RANKED" ] || exit 0
FILES=$(printf '%s\n' "$RANKED" | cut -f3)

# ---- 4. 分流：非 live 排最前（陈旧前提风险最高）；组内保持相关度序 ----
STALE=""; LIVE=""
while IFS="$TAB" read -r n d f; do
    [ -n "${f:-}" ] || continue
    base=$(basename "$f")
    rel=$(printf '%s' "$f" | sed "s|$NOTES/||")
    st=$(lookup_status "$base")
    case "$st" in
        live) LIVE="$LIVE$n$TAB$rel$TAB$st
" ;;
        *)    STALE="$STALE$n$TAB$rel$TAB$st
" ;;
    esac
done <<EOF
$RANKED
EOF

# ---- 5. 行级撤回标记（文件级状态装不下「一篇 live 里某个量被撤回」）----
INLINE=$(printf '%s\n' "$FILES" | while read -r f; do
    [ -n "$f" ] || continue
    if grep -qE '撤回|已作废|订正为|本节.*不成立' "$f" 2>/dev/null; then
        printf '%s\n' "$(basename "$f")"
    fi
done | head -4 | tr '\n' ' ')

# ---- 6. 输出（紧凑）----
echo "【先验检索】激活承重量："
printf '%s' "$ACTIVE" | grep -v '^$' | awk -F"$TAB" '{printf "  · %s (%s)\n", $2, $1}'

if [ -n "$STALE" ]; then
    echo "⚠ 非 live / 未索引命中（陈旧前提风险，先读这些）："
    printf '%s' "$STALE" | grep -v '^$' | head -$MAXSTALE \
        | awk -F"$TAB" '{printf "  ! %s  [%s]  (命中 %s/'"$NACTIVE"' 个量)\n", $2, $3, $1}'
fi

if [ -n "$LIVE" ]; then
    n=$(printf '%s' "$LIVE" | grep -c . )
    echo "✓ live 命中 $n 篇（按「同时命中几个承重量」降序，取前 $MAXLIVE）："
    printf '%s' "$LIVE" | grep -v '^$' | head -$MAXLIVE \
        | awk -F"$TAB" '{printf "  - %s  (命中 %s/'"$NACTIVE"' 个量)\n", $2, $1}'
fi

[ -n "$INLINE" ] && echo "⚠ 含行内撤回/订正标记：$INLINE"

echo "→ 这不是判决。重算作为独立二次验证是正当的；本条只是确保「已存在」不被忽略。"
exit 0
