# -*- coding: utf-8 -*-
"""
独立复核 §11.1 的订正界 B_i —— **逐箭头**，不复用主脚本的任何函数。

动机（CLAUDE.md「数值一致性检验必须逐箭头对齐」）：
旧界 `4·M2·ln(i+1)/i² + 2(λ+π²M2/6)/(i+1)` 的 GATE K3 报 0 违反，
是因为它只测 `R_i ≤ B_i`（两个箭头的**和**），而第二项压倒第一项把洞盖住。
⟹ 本脚本把证明拆成 5 个箭头，**每个单独测**，并对旧界做负对照（必须抓到）。

R_i := Σ_{j≥0}[G(|i−j|)+G(i+j+1)]/(2j+1) − S_i/(2i+1)，G(k)=Σ_j λ_j ω_j²/(ω_j²+k²)。
全部 float64（numerical-evidence）——B_i 的成立由解析推导承重，本脚本只查笔误。
"""
import numpy as np

PI2_6 = np.pi**2 / 6.0


def G(lams, wbars, k):
    lams = np.atleast_1d(np.asarray(lams, float))
    wbars = np.atleast_1d(np.asarray(wbars, float))
    w2 = wbars**2
    return float(np.sum(lams * w2 / (w2 + float(k) ** 2)))


def R_exact(lams, wbars, i, Jmax=40000):
    """R_i 的直接求和（Jmax 足够大；尾巴 ~M2/(i·J) 已在 1e-9 以下）。"""
    lam = float(np.sum(np.atleast_1d(lams)))
    s = 0.0
    for j in range(0, Jmax + 1):
        s += (G(lams, wbars, abs(i - j)) + G(lams, wbars, i + j + 1)) / (2 * j + 1)
    S_i = lam + 2.0 * sum(G(lams, wbars, k) for k in range(1, i + 1))
    return s - S_i / (2 * i + 1)


def B_new(lam, M2, i):
    return (8.0 + 4.0 * np.log(i + 1)) * M2 / i**2 + (lam + np.pi**2 * M2 / 2) / i


def B_old(lam, M2, i):
    return 4.0 * M2 * np.log(i + 1) / i**2 + 2 * (lam + np.pi**2 * M2 / 6) / (i + 1)


FAMILIES = [
    ("Einstein λ=3.2953", [3.2953], [0.664060982]),
    ("2atom x=.4 ρ=3", [0.6 * 3.2953, 0.4 * 3.2953], [0.664060982, 3 * 0.664060982]),
    ("2atom x=.1 ρ=20", [0.9 * 3.2953, 0.1 * 3.2953], [0.664060982, 20 * 0.664060982]),
    ("软端 ρ=0.05", [1.0, 2.0], [0.05, 1.3]),
    ("硬端 ρ=200", [2.0, 1.0], [0.3, 60.0]),
]

print("=" * 96)
print("箭头 A1：Σ_{j=0..h} 1/(2j+1) ≤ 1 + ½ln(i+1)，h=⌊i/2⌋  —— 旧界漏掉的就是这条的 j=0 项")
print("=" * 96)
bad_new = bad_old = 0
worst_new = (0.0, None)
for i in range(1, 5001):
    h = i // 2
    lhs = sum(1.0 / (2 * j + 1) for j in range(0, h + 1))
    rhs_new = 1.0 + 0.5 * np.log(i + 1)
    rhs_old = 0.5 * np.log(i + 1)          # 旧界隐含用的那条
    if lhs > rhs_new * (1 + 1e-14):
        bad_new += 1
    if lhs > rhs_old * (1 + 1e-14):
        bad_old += 1
    r = lhs / rhs_new
    if r > worst_new[0]:
        worst_new = (r, i)
print(f"  新版 `1+½ln(i+1)`：违反 {bad_new}/5000，最劣比 {worst_new[0]:.4f} (i={worst_new[1]})")
print(f"  旧版 `½ln(i+1)`  ：违反 {bad_old}/5000  ⟸ **负对照：必须大量违反，否则本脚本无效**")
print(f"    例：i=2 ⟹ 左={1+1/3:.4f}  旧右={0.5*np.log(3):.4f}  新右={1+0.5*np.log(3):.4f}")

print()
print("=" * 96)
print("箭头 A2：j≤h 段 ≤ (8·M2/i²)(1+½ln(i+1))    [单独测，不与 A3 相加]")
print("箭头 A3：j>h 段 ≤ (λ + π²M2/2)/i           [单独测]")
print("箭头 A4：包络 G(k) ≤ M2/k² (k≥1)，G(0)=λ")
print("=" * 96)
JT = 40000
for name, lams, wbars in FAMILIES:
    lam = float(np.sum(lams))
    M2 = float(np.sum(np.asarray(lams, float) * np.asarray(wbars, float) ** 2))
    v2 = v3 = v4 = 0
    w2 = w3 = 0.0
    for i in list(range(1, 61)) + [100, 200, 500, 1000]:
        h = i // 2
        seg1 = sum((G(lams, wbars, i - j) + G(lams, wbars, i + j + 1)) / (2 * j + 1)
                   for j in range(0, h + 1))
        b1 = (8.0 * M2 / i**2) * (1.0 + 0.5 * np.log(i + 1))
        seg2 = sum((G(lams, wbars, abs(i - j)) + G(lams, wbars, i + j + 1)) / (2 * j + 1)
                   for j in range(h + 1, JT + 1))
        b2 = (lam + np.pi**2 * M2 / 2) / i
        if seg1 > b1 * (1 + 1e-12):
            v2 += 1
        if seg2 > b2 * (1 + 1e-12):
            v3 += 1
        w2 = max(w2, seg1 / b1)
        w3 = max(w3, seg2 / b2)
    for k in range(1, 3000):
        if G(lams, wbars, k) > M2 / k**2 * (1 + 1e-12):
            v4 += 1
    okG0 = abs(G(lams, wbars, 0) - lam) < 1e-13
    print(f"  {name:22s} λ={lam:7.4f} M2={M2:9.4f} | A2 违反 {v2:2d} (最劣比 {w2:.4f}) | "
          f"A3 违反 {v3:2d} (最劣比 {w3:.4f}) | A4 违反 {v4} | G(0)=λ {okG0}")

print()
print("=" * 96)
print("箭头 A5：B_i 单调递减（sup_{i>I}R_i ≤ B_{I+1} 靠它）—— 解析：")
print("  d/di[(8+4ln(i+1))/i²] = [4i/(i+1) − 16 − 8ln(i+1)]/i³ < 0（分子 ≤ 4−16 < 0）；1/i 递减 ⟹ B_i 递减")
print("=" * 96)
for name, lams, wbars in FAMILIES[:3]:
    lam = float(np.sum(lams))
    M2 = float(np.sum(np.asarray(lams, float) * np.asarray(wbars, float) ** 2))
    bs = [B_new(lam, M2, i) for i in range(1, 20001)]
    inc = [i + 1 for i in range(len(bs) - 1) if bs[i + 1] > bs[i]]
    print(f"  {name:22s} i=1..20000 递增处 {len(inc)} 个；B_1={bs[0]:.4f} B_20000={bs[-1]:.3e}")

print()
print("=" * 96)
print("合门 A6：R_i ≤ B_i（新）对上；并对**旧界**做负对照 —— 旧界在小 i 上应当被抓到")
print("=" * 96)
for name, lams, wbars in FAMILIES:
    lam = float(np.sum(lams))
    M2 = float(np.sum(np.asarray(lams, float) * np.asarray(wbars, float) ** 2))
    vn = vo = 0
    wn = wo = 0.0
    argo = None
    for i in range(1, 41):
        r = R_exact(lams, wbars, i, JT)
        bn, bo = B_new(lam, M2, i), B_old(lam, M2, i)
        if r > bn:
            vn += 1
        if r > bo:
            vo += 1
            if argo is None:
                argo = i
        wn = max(wn, r / bn)
        if r / bo > wo:
            wo, argo2 = r / bo, i
    print(f"  {name:22s} | 新界违反 {vn:2d}/40 最劣 {wn:.4f} | 旧界违反 {vo:2d}/40 最劣 {wo:.4f}")

print()
print("读法：A1 的旧版负对照必须大量违反（证明本脚本有分辨力）；A2/A3/A4 各自 0 违反")
print("      才算「两个箭头都被单独测过」。A6 的旧界在总和口径上仍可能 0 违反——")
print("      **那正是旧缺陷当初没被 GATE K3 抓到的原因**，不能用它为旧界背书。")
