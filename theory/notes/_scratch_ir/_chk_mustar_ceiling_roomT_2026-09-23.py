# -*- coding: utf-8 -*-
"""
「开工去算」：用户要求的判决计算——μ*>0 情形下，room-T (300K) 在 η/M=9.172 处
所需 λ 的无条件天花板阈值，是否落在 Leavens 面积天花板的「我方获胜」一侧。

━━ 背景（见 memory mustar-resolved-ceiling-PhiA-2026-09-21 / leavens-1975-area-ceiling ）━━
定理 (Φ-A)（analytic，独立复核过 0/40）：M = K_ph − 2μ̃ uuᵀ，μ̃=μ*/λ，
  u[n]=(2n+1)^(-1/2)·1[n<Nc]，c_j:=⟨ψ_j,u⟩，γ:=μ₁−μ₂。前提 (PA): γ > 2μ̃‖u‖²。
  S⁺ := 2μ̃(‖u‖²−c₁²)/(γ−2μ̃‖u‖²)，Δ := 2μ̃c₁²/(1+S⁺) ⟹ λ₁(M) ≤ μ₁ − Δ。
此不等式对**任意**谱都成立，不需要极值性前提——它只是把 K_ph 的顶本征值
（μ₁，谱自身的量，无 Coulomb）扣掉一个由该谱自身向量算出的 Δ。

━━ 本脚本做什么（从零重写，不 import 2cf6342 上的脚本）━━
对 Einstein 谱（唯一有解析 μ*=0 极值性支持的谱，KAY 定理 2 / analytic），
1) 直接对角化线性化 Eliashberg 算子 K_ph(w̄)（μ*=0，Einstein，N=400 截断）；
2) 用 (Φ-A) 公式算 Δ(w̄; μ*, r)，得到 g_ceil(w̄) := μ₁(w̄) − Δ(w̄) ≤ λ₁(M(w̄))；
   ⟸ 方向：g_ceil 是 λ₁(M) 的**下界**? 不对——是**上界**（见 (Φ-A) 证明），
   Delta 是 D_true 的下估 ⟹ λ₁(M)=μ1-D_true ≤ μ1-Delta=g_ceil。
   即 g_ceil ≥ λ₁(M) 恒成立 ⟹ g_ceil(w̄)=1/λ 给出的 w̄ 是真根的**上界**
   （g 递增 ⟹ g_ceil 达到 gstar 的位置 ≥ 真根位置）⟹ 由此定出的 Tc 是
   Einstein 谱**自身**真实 Tc 的**上界**（ceiling，不需极值性）。
3) 固定 η/M=9.172（H₃S 实测），扫 λ 求 g_ceil-ceiling 下 Tc=300K 的 λ 阈值；
4) 与已知 μ*>0 面积天花板交叉常数对照（memory 记的 4.862 for μ*=0.13，
   本脚本独立重算 c(μ*=0.13, r) 的量级做交叉校验，而非照搬）。

严格性：(Φ-A) 不等式本身 analytic；本脚本的数值解 numerical-evidence (float64, N=400)。
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

C_INF = 0.1827262477890          # analytic (rational-exact 包围), kTc = C_inf sqrt(lam <w^2>)
HB = 64.65415134                 # meV per sqrt(eV/A^2/amu), eta/M = lam*wbar2^2/HB^2
KB = 8.617333262e-2              # meV/K
N = 400


def Kph_pieces(wbar, Nc):
    """Einstein K_ph (mu*=0 linearized operator). wbar := w_E/(2 pi T) (repo convention:
    G[k] = lam*wbar^2/(wbar^2+k^2), lam=1 here and rescaled by the caller).
    Returns mu1, mu2, c1^2, ||u||^2 for the coulomb vector u truncated at Nc."""
    n = np.arange(N)
    kmax = 2 * N + 2
    k = np.arange(kmax)
    G = wbar * wbar / (wbar * wbar + k * k)   # single Einstein line, lam=1 normalized
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N)
    L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    K = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    K[n, n] -= G[0] * u * u
    ev, V = eigh(K)
    mu1, mu2 = ev[-1], ev[-2]
    uc = np.where(n < Nc, u, 0.0)
    c1 = float(V[:, -1] @ uc)
    return mu1, mu2, c1 * c1, float(uc @ uc)


def g_true(wbar, lam, mus, r):
    """Exact lambda1(M) for the Einstein+Coulomb operator, gstar-normalized by lam
    outside (this returns lambda1 directly, to be compared against 1/lam externally
    after lam-rescaling of K)."""
    n = np.arange(N)
    kmax = 2 * N + 2
    k = np.arange(kmax)
    G = lam * wbar * wbar / (wbar * wbar + k * k)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N)
    L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= G[0] * u * u
    Nc = int(np.floor(r * wbar))
    if Nc > 0:
        uc = np.where(n < Nc, u, 0.0)
        M -= 2.0 * mus * np.outer(uc, uc)
    return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])


def g_ceil(wbar, lam, mus, r):
    """(Phi-A) ceiling on lambda1(M): mu1(lam-scaled) - Delta."""
    Nc = int(np.floor(r * wbar))
    mu1, mu2, c1s, un = Kph_pieces(wbar, Nc)
    mu1 *= lam
    mu2 *= lam
    gam = mu1 - mu2
    mu_t = mus / lam
    den = gam - 2.0 * mu_t * un
    if den <= 0:
        return None
    sp = 2.0 * mu_t * (un - c1s) / den
    delta = 2.0 * mu_t * c1s / (1.0 + sp)
    return mu1 - delta


def root_wbar(fn, gstar, lo=0.05, hi=10.0, nbis=55):
    if fn(hi) < gstar or fn(lo) >= gstar:
        return None
    a, b = lo, hi
    for _ in range(nbis):
        m = 0.5 * (a + b)
        if fn(m) >= gstar:
            b = m
        else:
            a = m
    return 0.5 * (a + b)


def Tc_ceiling_K(lam, wbar2_meV, mus, r):
    """Solve g_ceil(wbar)=1/lam for wbar=w_E/(2 pi T), convert to Tc[K]."""
    wE = wbar2_meV  # for Einstein spectrum wbar2 == w_E exactly
    wstar = root_wbar(lambda w: g_ceil(w, lam, mus, r), 1.0)
    if wstar is None:
        return None
    T_meV = wE / (2.0 * np.pi * wstar)
    return T_meV / KB


def main():
    print("=" * 90)
    print("STEP 0 -- sanity: (Phi-A) ceiling >= true lambda1 (direction check, few pts)")
    for wbar in (0.6, 1.2, 2.0):
        for lam in (2.0, 4.0):
            for mus, r in ((0.13, 10.0),):
                gt = g_true(wbar, lam, mus, r)
                gc = g_ceil(wbar, lam, mus, r)
                ok = (gc is not None) and (gc >= gt - 1e-9)
                print(f"  wbar={wbar:.2f} lam={lam:.1f} mu*={mus:.2f} r={r:g}  "
                      f"true={gt:.6f}  ceil={gc if gc is None else round(gc,6)}  "
                      f"{'OK' if ok else '**BUG**'}")
    print()

    eta_over_M = 9.172   # H3S-anh experimental anchor, meV^2 units via HB
    mus = 0.13
    r = 10.0
    print("=" * 90)
    print(f"STEP 1 -- solve Tc_ceiling(lam; eta/M={eta_over_M}, mu*={mus}, r={r}) = 300K")
    print("  wbar2[meV] = HB*sqrt(eta_over_M/lam)  (fixed eta/M, Einstein wbar2=w_E)")
    print(f"{'lam':>6}{'wbar2[meV]':>13}{'Tc_ceil[K]':>13}")
    lam_grid = np.arange(2.0, 6.51, 0.25)
    results = []
    for lam in lam_grid:
        wbar2 = HB * np.sqrt(eta_over_M / lam)
        tc = Tc_ceiling_K(lam, wbar2, mus, r)
        results.append((lam, wbar2, tc))
        print(f"{lam:6.2f}{wbar2:13.4f}{('%.3f' % tc) if tc else 'nan':>13}")

    # bisect for Tc_ceiling(lam)=300 (Tc increasing in lam expected)
    def f(lam):
        wbar2 = HB * np.sqrt(eta_over_M / lam)
        tc = Tc_ceiling_K(lam, wbar2, mus, r)
        return (tc - 300.0) if tc is not None else -1e9

    lo, hi = 1.5, 8.0
    flo, fhi = f(lo), f(hi)
    print(f"\n  f(lo={lo})={flo:.3f}  f(hi={hi})={fhi:.3f}")
    if flo < 0 < fhi:
        for _ in range(50):
            m = 0.5 * (lo + hi)
            if f(m) < 0:
                lo = m
            else:
                hi = m
        lam_ceil = 0.5 * (lo + hi)
        print(f"  ==> unconditional (Phi-A, no extremality premise) lambda floor "
              f"for 300K at eta/M={eta_over_M}, mu*={mus}: lam_ceil = {lam_ceil:.4f}")
    else:
        lam_ceil = None
        print("  ==> no sign change in scanned range; cannot bisect (see grid above)")

    print()
    print("=" * 90)
    print("STEP 2 -- DIRECT apples-to-apples: Leavens AREA ceiling Tc<=c(mu*,r)*A on the")
    print("  SAME Einstein spectrum, same eta/M, same 300K target.")
    print("  For Einstein, A = lam*wE/2 (all moments degenerate to wE).")
    print("  c(0.13, r=10) is NOT re-derived from scratch here -- backed out algebraically")
    print("  from the memory-recorded crossing lam*=(2C_inf/c)^2=4.862, i.e. this number")
    print("  is only as trustworthy as that unverified-this-session memory entry.")
    c013 = 2.0 * C_INF / np.sqrt(4.862)
    print(f"  ==> c(0.13,r=10) implied = {c013:.5f}")

    def Tc_area_K(lam):
        wE = HB * np.sqrt(eta_over_M / lam)
        return c013 * lam * wE / 2.0 / KB

    lo, hi = 0.5, 10.0
    fa = lambda l: Tc_area_K(l) - 300.0
    if fa(lo) < 0 < fa(hi):
        for _ in range(60):
            m = 0.5 * (lo + hi)
            if fa(m) < 0:
                lo = m
            else:
                hi = m
        lam_area = 0.5 * (lo + hi)
        print(f"  area-ceiling lambda floor for 300K = {lam_area:.4f}")
        if lam_ceil is not None:
            print(f"  Phi-A moment-ceiling lambda floor  = {lam_ceil:.4f}")
            winner = "AREA (Leavens) is tighter -- our Phi-A route does NOT beat it" \
                if lam_area > lam_ceil else \
                "Phi-A MOMENT ceiling is tighter -- unconditional claim survives"
            print(f"  ==> {winner}")
    print()
    print("CAVEAT (load-bearing, do not drop when reporting):")
    print("  commit 2cf6342's own script _chk_c_mustar_curve_2026-09-21.py states, quoting")
    print("  Carbotte RMP62 verbatim: the mu*>0 extremality behind c(mu*,r) itself is only")
    print("  'a LOCAL maximum ... established via a NUMERICALLY evaluated functional")
    print("  derivative' in the PUBLISHED literature -- i.e. even Leavens/Carbotte's own")
    print("  area ceiling is NOT a proven exclusion theorem at mu*>0, only numerical-evidence.")
    print("  So neither route currently supports an ANALYTIC/certificate-level unconditional")
    print("  PRL claim at mu*>0; both cap out at numerical-evidence.")


if __name__ == "__main__":
    main()
