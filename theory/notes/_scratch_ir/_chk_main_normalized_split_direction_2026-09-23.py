# -*- coding: utf-8 -*-
"""
复核 `extremality-hunt` 的收尾报告：他的「细扫」结论是否又一次栽进 §5 的退化成员劫持？

━━ 被复核的陈述 ━━
(H1)「细扫（λ=5, μ*=0.10, r=10，eps∈[0.05,0.60]×delta∈{0.005..0.3}）真正的 max Δ 在**最小**的
     (eps,delta)=(0.05,0.005)，Δ=−1.80e−4；此前粗网格里 (0.35,0.02)『看起来最优』是网格太粗的假象。」
(H2)「危险方向不在『最铺开的 split』，而是**紧贴简并点本身**。」
(H3)「粗扫 432 点 Δ>0 计数 = 0，max Δ = −3.59e−3。」

━━ 我的判词设计 ━━
`Δ := λ₁(split,T) − λ₁(Einstein,T)` 在 split → Einstein 时**连续趋于 0**。
该族里 Δ 全负 ⟹ **`max Δ` 必然被最接近 Einstein 的那个成员取到**，这是**连续性的恒等式**，
不是「那里最危险」。本篇笔记 §5 已记过同形事故（门 4 首版两列逐行恰差一负号，
极值被近简并谱 (0.999,1.001) 平凡取到，自由度实为 1），修法是**按方向长度归一**：

    Δ̂ := Δ / ‖ΔG‖₁ ,  ‖ΔG‖₁ := Σ_{k≥1} |G_E[k] − G_split[k]|   （G[0] ≡ λ 逐谱相同，不计）

⟹ 判据：(H1) 在**未归一**下为真但**无信息**；真问题是 Δ̂ 的 argmax 落在哪。
   若 Δ̂ 的 argmax 搬离简并角 ⟹ (H2)「危险方向 = 紧贴简并点」**为假**，
   且他"跑到头了"的依据不成立（他按一个平凡统计量选的方向）。

━━ 方向审计 ━━
本脚本**不产出**任何阳性结论。两种结果都只支持「未找到反例」：
· 若 Δ̂ argmax 仍在简并角 ⟹ 他的定性判断侥幸对，但理由仍要换。
· 若 Δ̂ argmax 搬走 ⟹ 他的收尾判据无效，搜索**未**跑到头。
无论哪种，Δ<0 都只支持「未找到反例」，**不支持**「已证极值性」。

严格性：numerical-evidence（float64，N=400，有限网格）。
"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

KB = 8.617333262e-2
WB = 150.0
N = 400


def Gvec(T, aw, kmax=2 * N + 2):
    x2 = (2.0 * np.pi * T * np.arange(kmax, dtype=float)) ** 2
    G = np.zeros(kmax)
    for a, w in aw:
        G += 2.0 * a * w / (w * w + x2)
    return G


def lam1_G(T, G, wc, mus):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N); L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= G[0] * u * u
    Nc = int(np.floor(wc / (2.0 * np.pi * T)))
    if Nc > 0:
        uc = np.where(n < Nc, u, 0.0)
        M -= 2.0 * mus * np.outer(uc, uc)
    return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])


def T_gate_lo(wc):
    return wc / (2.0 * np.pi * (N + 1)) / KB


def Tc(aw, wc, mus, hi=900.0, nb=44):
    f = lambda T: lam1_G(T, Gvec(T, aw), wc, mus) - 1.0
    lo, hi = T_gate_lo(wc) * KB, hi * KB
    if not (f(lo) > 0 > f(hi)):
        return np.nan
    for _ in range(nb):
        m = np.sqrt(lo * hi)
        if f(m) > 0: lo = m
        else: hi = m
    return np.sqrt(lo * hi) / KB


def two_atom(lam, w1, w2, wbar=WB):
    p, q = lam / 2.0, lam * wbar ** 2 / 2.0
    d = w2 ** 2 - w1 ** 2
    if d <= 1e-12: return None
    a1 = w1 * (p * w2 ** 2 - q) / d
    a2 = w2 * (q - p * w1 ** 2) / d
    if a1 < -1e-14 or a2 < -1e-14: return None
    return [(max(a1, 0.0), w1), (max(a2, 0.0), w2)]


LAM, MUS, R = 5.0, 0.10, 10.0
wc = R * WB
E = [(LAM * WB / 2.0, WB)]
tE = Tc(E, wc, MUS); T = tE * KB
GE = Gvec(T, E)
lE = lam1_G(T, GE, wc, MUS)
x = wc / (2.0 * np.pi * T)
print("=" * 104)
print(f"参照点：λ={LAM}, μ*={MUS}, r={R}, N={N}   Tc_E = {tE:.6f} K")
print(f"  N_c = {int(np.floor(x))}  (N_c/N = {np.floor(x)/N:.3f}, 合法)   "
      f"dist(ω_c/2πTc → 整数) = {min(x%1, 1-x%1):.4f}")
print(f"  λ₁(Einstein, T) − 1 = {lE-1:+.3e}   （基准实测，不是名义 1）")
print("=" * 104)

# 他的细扫网格
EPS = np.round(np.arange(0.05, 0.601, 0.05), 3)
DEL = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3)

rows, t0 = [], time.time()
for e in EPS:
    for d in DEL:
        w1, w2 = (1.0 - e) * WB, (1.0 + d) * WB
        a = two_atom(LAM, w1, w2)
        if a is None:
            continue
        G = Gvec(T, a)
        dn = float(np.abs(G[1:] - GE[1:]).sum())      # G[0] ≡ λ 逐谱相同，不计
        if dn < 1e-12:
            continue
        dl = lam1_G(T, G, wc, MUS) - lE
        rows.append((e, d, dl, dn, dl / dn))
print(f"（{len(rows)} 点，用时 {time.time()-t0:.0f}s）\n")

raw = max(rows, key=lambda z: z[2])
nrm = max(rows, key=lambda z: z[4])
print(f"  {'统计量':<28}{'argmax (eps,delta)':>22}{'值':>16}{'‖ΔG‖₁':>14}")
print(f"  {'未归一 Δ（他报的）':<26}{f'({raw[0]:.2f}, {raw[1]:.3f})':>22}{raw[2]:>+16.3e}{raw[3]:>14.3e}")
print(f"  {'归一 Δ̂ = Δ/‖ΔG‖₁':<26}{f'({nrm[0]:.2f}, {nrm[1]:.3f})':>22}{nrm[4]:>+16.3e}{nrm[3]:>14.3e}")
print(f"\n  Δ>0 计数 = {sum(1 for r in rows if r[2] > 0)} / {len(rows)}")

print("\n  按未归一 Δ 排序（前 5）—— 他的口径：")
for r in sorted(rows, key=lambda z: -z[2])[:5]:
    print(f"    eps={r[0]:.2f} delta={r[1]:.3f}   Δ={r[2]:+.3e}   ‖ΔG‖₁={r[3]:.3e}   Δ̂={r[4]:+.3e}")
print("\n  按归一 Δ̂ 排序（前 5）—— 修正口径：")
for r in sorted(rows, key=lambda z: -z[4])[:5]:
    print(f"    eps={r[0]:.2f} delta={r[1]:.3f}   Δ={r[2]:+.3e}   ‖ΔG‖₁={r[3]:.3e}   Δ̂={r[4]:+.3e}")

# 简并极限：Δ 与 ‖ΔG‖₁ 是否同阶趋零？（若 Δ/‖ΔG‖₁ → 有限非零，argmax 落简并角就是平凡的）
print("\n  简并极限诊断（沿 eps=delta→0 的对角线）：")
print(f"  {'s':>10}{'‖ΔG‖₁':>14}{'Δ':>14}{'Δ̂=Δ/‖ΔG‖₁':>16}")
for s in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002):
    a = two_atom(LAM, (1.0 - s) * WB, (1.0 + s) * WB)
    if a is None: continue
    G = Gvec(T, a); dn = float(np.abs(G[1:] - GE[1:]).sum())
    dl = lam1_G(T, G, wc, MUS) - lE
    print(f"  {s:10.3f}{dn:14.3e}{dl:+14.3e}{dl/dn:+16.3e}")
print("  ⟹ 若 Δ 与 ‖ΔG‖₁ 同阶趋 0，则『max Δ 在最小 split』是**连续性恒等式**，零信息。")

print()
print("=" * 104)
print("判词（单边界方向规则；两种结果都只支持「未找到反例」）：")
print("  · (H3) 的 0 违反：复核方向一致，不受影响。")
print("  · (H1) 未归一 max Δ 的 argmax 落在最小 split —— 为真但**零信息**（Δ→0⁻ 的连续性）。")
print("  · (H2)「危险方向 = 紧贴简并点」—— 由上表的 Δ̂ argmax 判定。")
print("=" * 104)
