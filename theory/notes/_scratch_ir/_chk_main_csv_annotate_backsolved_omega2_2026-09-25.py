# -*- coding: utf-8 -*-
"""把 §3 反解出的 ω̄₂ 写进 hydrides_theory.csv 的 notes 字段（**不**写进 omega_2_meV 列）。

为什么不填 omega_2_meV：该列其余行装的是一手 α²F 积分（H3S/D3S 四行）或论文印出值
（YH6/YH9，Ding Table I）。反解值是**派生量**，信息类不同；填进去会让下游
「非 not-reported 即可用」的筛选静默吃进派生值。⟹ 列保持 not-reported（事实如此），
值写进 notes，前缀 `omega_2_BACKSOLVED_meV=` 便于 grep。

运行：OMP_NUM_THREADS=1 python theory/notes/_scratch_ir/_chk_main_csv_annotate_backsolved_omega2_2026-09-25.py
"""
import io
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

CSV = "theory/data/hydrides_theory.csv"

# (material, pressure) -> (读法一 meV, 读法二 meV)
BS = {
    ("LaH10", "129"): (100.96, 100.76),
    ("LaH10", "163"): (108.95, 120.88),
    ("LaH10", "214"): (142.51, 141.22),
    ("LaH10", "264"): (156.58, 155.35),
    ("LaD10", "159"): (79.23, 78.88),
    ("LaD10", "210"): (97.35, 97.39),
    ("LaD10", "260"): (109.10, 107.61),
}
TAG = "omega_2_BACKSOLVED_meV="

with open(CSV, encoding="utf-8") as fh:
    lines = fh.read().split("\n")

n = 0
out = []
for ln in lines:
    parts = ln.split(",")
    key = (parts[0], parts[2]) if len(parts) > 3 else None
    if key in BS and TAG not in ln:
        a, b = BS[key]
        note = (f" ★2026-09-25 主会话反解：{TAG}{a:.2f}（读法一，照印出 Tc^Mc）"
                f" / {b:.2f}（读法二，用公式值 Tc^Mc）"
                f"——由 f2=(Tc^AD/Tc^Mc)/f1 反解，**派生量、非论文印出值**，"
                f"故 omega_2_meV 列仍留 not-reported；"
                f"依据与方向审计见 theory/notes/ceiling/errea_f2_backsolve_omega2_2026-09-25.md")
        ln = ln + note
        n += 1
    out.append(ln)

assert n == 7, f"应改 7 行，实改 {n} 行"
with open(CSV, "w", encoding="utf-8", newline="") as fh:
    fh.write("\n".join(out))
print(f"已标注 {n} 行；omega_2_meV 列未改动（仍为 not-reported）。")

# 自检：重读，确认 omega_2_meV 列没被碰
with open(CSV, encoding="utf-8") as fh:
    for ln in fh:
        p = ln.split(",")
        if len(p) > 9 and (p[0], p[2]) in BS:
            assert p[8] == "not-reported", f"omega_2_meV 被改了: {p[0]}@{p[2]} -> {p[8]}"
print("自检 PASS：7 行的 omega_2_meV 仍为 not-reported。")
