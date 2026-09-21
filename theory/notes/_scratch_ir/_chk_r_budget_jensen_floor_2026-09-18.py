#!/usr/bin/env python3
"""r=ω₂/ω_log 的预算：把 r 退化到 Jensen 下确界 r=1，四个 LaH10 行的判决翻转吗？

背景（主会话 2026-09-18 溯源）：`P4_eta_over_M.md` §4.A 写「r 只有 LaH10 四点有
（p4_tc_table.py LAH10_POINTS，P3 一手反演）」。追到 `P3_tc_max_pressure.md` §2，
那个「反演」是**从 Errea Ext. Data Table I 同给的 Tc^Mc 与 Tc^AD 两列，用 Allen–Dynes
f₁f₂ 的确切函数形式解二次方程**得到的，不是 ω₂ 的直接测量/列表（P3 明写「一手 ω₂
列表确认不存在」）。⟹ r 继承了「Errea 用的就是这一支 AD f₂ 形式 + 前因子 1.20」这个前提。

方向：r 被**高估** ⟹ η/M=λ(rω_log)² 高估 ⟹ 界 A 高估 ⟹「界A > Tc」是 fail-open。
故必须算预算：r 退到解析下确界 r≥1（Jensen，⟨ω²⟩≥exp⟨lnω²⟩）后判决还在吗？
"""
import numpy as np
C_INF = 0.1827262477890
MEV_TO_K = 1e-3 / 8.617333262e-5

ROWS = [  # (name, λ, ω_log, r_inverted, Tc^ME)
    ("LaH10 129", 3.62,  76.40, 91.22 / 76.40, 255.3),
    ("LaH10 163", 2.67,  96.40, 103.63 / 96.40, 242.8),
    ("LaH10 214", 2.06, 115.50, 127.63 / 115.50, 237.9),
    ("LaH10 264", 1.73, 126.60, 138.63 / 126.60, 216.9),
]

def boundA(lam, w2):
    return C_INF * np.sqrt(lam * w2 * w2) * MEV_TO_K

print(f"{'row':<11}{'r_inv':>7}{'界A(r_inv)':>11}{'界A(r=1)':>10}{'Tc^ME':>8}"
      f"{'紧度(r_inv)':>12}{'紧度(r=1)':>11}{'判决@r=1':>10}")
allpass = True
for name, lam, wlog, r, tc in ROWS:
    bA = boundA(lam, wlog * r)
    bA1 = boundA(lam, wlog)          # r 退到 Jensen 下确界
    ok = bA1 > tc
    allpass &= ok
    print(f"{name:<11}{r:7.3f}{bA:11.1f}{bA1:10.1f}{tc:8.1f}"
          f"{bA/tc:12.3f}{bA1/tc:11.3f}{'PASS' if ok else 'FAIL':>10}")

print()
print("结论：" + ("四行在 r=1（解析下确界）处**全部仍 PASS** ⟹ "
      "「界A > Tc^ME」这个判决不依赖反演出来的 r，"
      "r 只影响『紧度 1.44–1.78×』这类比较性数字。"
      if allpass else "有行翻转 ⟹ r 是承重输入，反演前提必须升级为承重引用。"))
print("⟹ r 的预算对**判决**是 ∞（删掉照样过），对**紧度数字**是有限的：")
for name, lam, wlog, r, tc in ROWS:
    print(f"   {name}: 紧度从 {boundA(lam, wlog*r)/tc:.3f} 掉到 {boundA(lam, wlog)/tc:.3f}"
          f"（差 {100*(boundA(lam,wlog*r)/tc/(boundA(lam,wlog)/tc)-1):.1f}%）")
