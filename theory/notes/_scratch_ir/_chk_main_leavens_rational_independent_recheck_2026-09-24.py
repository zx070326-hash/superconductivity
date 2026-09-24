# -*- coding: utf-8 -*-
"""GATE L — 我从零复核 `leavens_c_rational_exact_per_mustar_2026-09-24.md` 的四行。

**不 import 它的任何脚本、不抄它的 w 取点、不抄它的数。** 我自己：
  (1) float64 逐齿求 `c(μ*) = max_w θ_max(w)/(π w)` 的 argmax（齿区间 [k/r,(k+1)/r) 左闭）；
  (2) 在 argmax 附近取一个**有理** w，用 `Fraction` 精确算广义 Rayleigh 商
      `θ_lo = N(y)/D(y) ≤ θ_max`，`y_n := u_n v_n`（把 u 吸进向量 ⟹ 全部系数有理）；
  (3) `c_lo = θ_lo/(π_hi · w)`，`π_hi = 355/113 > π`（方向：分母取大 ⟹ c 取小 ⟹ 下界）；
  (4) `λ_area_up = 1/(π_lo² · c_lo² · M2)`，`π_lo = 3.14159265358979 < π`
      （`λ_area = 1/(π²c²M2)` 对 c 严格递减 ⟹ c 的下界给 λ_area 的**上界**）；
  (5) `ratio_lo = λ_ceil / λ_area_up`（λ_ceil = 我自己解的 Einstein 阈值，float64）。

方向审计（每一步都朝「ratio 变小」退化 ⟹ `ratio_lo > 1` 只支持 PASS）：
  θ_lo ≤ θ_max ⟹ c_lo ≤ c ⟹ λ_area_up ≥ λ_area ⟹ ratio_lo ≤ ratio。单边界方向规则。

严格性：c 的下界一列是 **rational-exact**（格数与有理数个数逐行印出）；
λ_ceil 与 λ_area 整体仍是 **numerical-evidence**（M2 与 λ_ceil 未区间化）。
"""
import os
import sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")   # GBK 终端下 U+27F9 会崩
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh

R = 10.0
M2 = 1.453151466                      # η/M=9.172, T=300K（与头条表同一常数）
PI_HI = F(355, 113)                   # > π
PI_LO = F(314159265358979, 10**14)    # < π
ND = 40                               # 有理格数
NF = 400                              # float64 求 argmax / 试探向量用的维数


# ───────── float64：θ_max(w)、c(w)、逐齿 argmax ─────────
def khat_B(w, mus, Nc, n):
    k = np.arange(2*n + 2); G = w*w/(w*w + k*k)
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    u = 1.0/np.sqrt(2.0*np.arange(n) + 1.0)
    S = G[0] + 2.0*(np.cumsum(G[:n]) - G[0])                  # 闭式 S_n
    K = (G[np.abs(i-j)] + G[i+j+1])*np.outer(u, u) - np.diag(u*u*S)
    uc = np.where(np.arange(n) < Nc, u, 0.0)
    B = np.eye(n) + 2.0*mus*np.outer(uc, uc)
    return K, B


def theta_max(w, mus, n=NF, vec=False):
    Nc = int(np.floor(R*w))
    K, B = khat_B(w, mus, Nc, n)
    e, V = eigh(K, B)
    return (float(e[-1]), V[:, -1]) if vec else float(e[-1])


def c_of_w(w, mus):
    return theta_max(w, mus)/(np.pi*w)


def per_tooth_max(mus, kmin=4, kmax=30, nper=13, nref=25):
    best = (-1.0, None)
    for k in range(kmin, kmax + 1):
        lo, hi = k/R, (k + 1)/R
        ws = np.linspace(lo, hi - 1e-12, nper)              # 齿左闭右开
        cs = [c_of_w(w, mus) for w in ws]
        a = int(np.argmax(cs))
        l2 = max(lo, ws[a] - (hi-lo)/nper); h2 = min(hi - 1e-12, ws[a] + (hi-lo)/nper)
        for w in np.linspace(l2, h2, nref):
            cc = c_of_w(w, mus)
            if cc > best[0]:
                best = (cc, w)
    return best


# ───────── rational-exact：θ_max 的有理下界 ─────────
def theta_lower_rational(wF, mus_F, nd=ND):
    p, q = wF.numerator, wF.denominator
    Nc = int(np.floor(R*float(wF)))
    Gf = [F(p*p, p*p + q*q*k*k) for k in range(2*nd + 2)]
    Sf = [Gf[0] + 2*sum(Gf[1:n+1], F(0)) for n in range(nd)]
    _, V = theta_max(float(wF), float(mus_F), n=nd, vec=True)
    DEN = 10**6
    y = [F(round(float(V[k])/np.sqrt(2*k + 1)*DEN), DEN) for k in range(nd)]
    num = sum((Gf[abs(i-j)] + Gf[i+j+1])*y[i]*y[j]
              for i in range(nd) for j in range(nd)) - sum(Sf[i]*y[i]*y[i] for i in range(nd))
    den = sum(y[i]*y[i]*F(2*i + 1) for i in range(nd)) + 2*mus_F*(sum(y[:Nc], F(0)))**2
    nrat = len(y) + len(Gf) + len(Sf)
    return num/den, Nc, nrat


# ───────── λ_ceil：sup λ₁ = 1 的 λ（M2 固定），float64 ─────────
def lam1_einstein(lam, mus):
    w = np.sqrt(M2/lam); Nc = int(np.floor(R*w))
    K, B = khat_B(w, mus, Nc, NF)
    # λ₁(λK̂ − 2μ*u_cu_c^T) ≥ 1 ⟺ λ ≥ 1/θ_max  (w-参数化引理)
    return lam*theta_max(w, mus)


def lam_ceil(mus, lo=0.5, hi=12.0, ngrid=1201):
    """λ_ceil := inf{λ : λ₁(λ) ≥ 1}。

    ⚠ **先扫变号次数，再二分**：`N_c=⌊r·√(M2/λ)⌋` 的台阶使 `λ₁(λ)−1` 可以多次变号
    （μ*=0.20 实测 3 次），裸二分不保证落到 **inf**-根，而物理陈述
    「λ ≤ λ_ceil ⟹ 达不到 300 K」要求的正是 inf-根（取更大的根 = fail-open）。
    本函数早先是裸二分，在 μ*=0.20 上**碰巧**落对；2026-09-24 改为显式取第一个变号区间。
    见 `_chk_main_lamceil_mu020_step_root_audit_2026-09-24.py`。
    """
    ls = np.linspace(lo, hi, ngrid)
    vs = np.array([lam1_einstein(L, mus) - 1.0 for L in ls])
    sc = np.where(np.diff(np.sign(vs)) != 0)[0]
    if len(sc) == 0:
        return float("nan")
    s = sc[0]                                   # ← 第一个变号区间 = inf-根所在
    a, b, va = ls[s], ls[s + 1], vs[s]
    for _ in range(80):
        m = 0.5*(a + b)
        if (lam1_einstein(m, mus) - 1.0)*va > 0:
            a = m
        else:
            b = m
    return 0.5*(a + b), len(sc)


M2F = F(1453151466, 10**9)
print("GATE L — 我从零复核 leavens-cert 的四行（不 import 其脚本、不抄其 w）")
print(f"{'mus':>5} {'w*(float)':>11} {'c_max(float)':>13} | {'w(有理)':>10} {'Nc':>3} "
      f"{'c_lo(rat-exact)':>17} {'亏损':>10} | {'λ_area_up':>10} {'λ_ceil':>9} {'ratio_lo':>9}")
NOTE = {0.00: 1.0142, 0.10: 1.1841, 0.13: 1.2961, 0.20: 1.7948}   # 0.20 已订正（inf-根）
for mus in (0.00, 0.10, 0.13, 0.20):
    cf, wf = per_tooth_max(mus)
    wF = F(round(wf*10000), 10000)
    thF, Nc, nrat = theta_lower_rational(wF, F(round(mus*100), 100))
    c_lo = thF/(PI_HI*wF)
    la_up = 1/(PI_LO*PI_LO*c_lo*c_lo*M2F)
    lc, nsc = lam_ceil(mus)
    ratio = F(round(lc*10**12), 10**12)/la_up
    print(f"{mus:5.2f} {wf:11.6f} {cf:13.9f} | {float(wF):10.4f} {Nc:3d} "
          f"{float(c_lo):17.9f} {cf-float(c_lo):10.2e} | {float(la_up):10.6f} {lc:9.6f} "
          f"{float(ratio):9.4f} [变号{nsc}次]  (笔记 {NOTE[mus]:.4f}, Δ={float(ratio)-NOTE[mus]:+.4f})"
          f"  [{ND} 格 / {nrat} 个有理数]")
print("\n方向：θ_lo≤θ_max ⟹ c_lo≤c ⟹ λ_area_up≥λ_area ⟹ ratio_lo≤ratio ⟹ 'ratio_lo>1' 只支持 PASS。")
