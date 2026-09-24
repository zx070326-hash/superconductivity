# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ①：尾部机制能否从 Einstein 迁移到两原子核？——从零重建。

━━━━━━━━ 本脚本要判的三件事 ━━━━━━━━

【I】**包络迁移**。`P4_coulomb_global_tc_certified.md` 的 Kato–Temple 尾部走
     `schur_coupling`，而整套 Schur 机制只依赖一个**幂律包络** `|G(k)| <= C/(k^2+A^2)^p`
     （`p4_cii_ci_interval.py` 的 `spec=(C,A,p)`、`env_at`、`row_env`、`col_env`）。
     Einstein: G(k)=lam*w^2/(w^2+k^2) 恰是 (C,A,p)=(lam*w^2, w, 1)。
     两原子: G(k)=sum_j lam_j w_j^2/(w_j^2+k^2)。**断言**（一行 analytic）：
         lam_j w_j^2/(w_j^2+k^2) <= lam_j w_j^2/(w_min^2+k^2)   （因 w_j >= w_min）
       ⟹ G(k) <= M2/(k^2 + w_min^2),   M2 := sum_j lam_j w_j^2 = lam*<w^2>
     **且 M2 在固定 (lam, wbar2) 的整族上是常数** ⟹ 包络的 C **族内不变**，
     只有 A=w_min 随成员变（最坏退化到 A=0，给 G(k) <= M2/k^2，对全族一致）。
     GATE I 数值查该包络不等式（含 A=0 的一致版）。

【II】**S_n 有闭式，不需要截断**（本轮发现，且它订正了我自己的脚本）。
     S_n := sum_{m>=0} [G(|n-m|) - G(n+m+1)]
          = [sum_{j=1}^{n}G(j) + G(0) + sum_{j=1}^{inf}G(j)] - sum_{j=n+1}^{inf}G(j)
          = **G(0) + 2*sum_{j=1}^{n} G(j)**      ← 有限和，精确，与 N 无关
     我此前全部脚本（`_chk_main_phiA_general_spectrum_sup_2026-09-23.py` 等）用的是
     **N x N 块的行和** S_n^(N) = sum_{m<N}[...] < S_n^inf ⟹ 对角偏大（欠负）⟹ lam1 偏大。
     GATE II 量出两者之差，并**用精确 S_n 重跑 (A6) 的族内 sup**
     （CLAUDE.md：链条上对象的性质改了，每个用到它的箭头都要重走）。

【III】**预算**（CLAUDE.md 预算优先，先算再派活）。承认认证包围有相对松 delta：
     `lam1 <= (1+delta) * lam1_true`。问：lam_ceil / Tc_ceil 会移动多少？
     消费点 = 室温阈值 lam_ceil。仓库既有认证宽度 4.2e-3（r=3, Einstein）⟹ 若
     delta~1e-3 只让阈值移动零点几个百分点，则**步骤②是可付的**；若杠杆很大则不是。

严格性：【I】【II】的代数为 analytic；本脚本一切数值为 numerical-evidence（float64）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh


# ───────────────────────── 算子（两种 S_n 口径） ─────────────────────────
def build(lams, wbars, mus, Nc, N, exact_S=True):
    """symmetrized M。exact_S=True 用闭式 S_n = G(0)+2*sum_{j<=n}G(j)（精确，无截断）；
    False 用 N x N 块行和（我此前脚本的口径）。"""
    lams = np.asarray(lams, float); wbars = np.asarray(wbars, float)
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    kmax = 2*N+2
    k = np.arange(kmax)
    G = np.zeros(kmax)
    for L, W in zip(lams, wbars):
        G += L*W*W/(W*W + k*k)
    i, j = np.meshgrid(n, n, indexing="ij")
    Kp = G[np.abs(i-j)]; Km = G[i+j+1]
    if exact_S:
        S = G[0] + 2.0*np.cumsum(G[1:N+1])   # S[n] = G(0)+2*sum_{j=1..n}G(j)
        S = np.concatenate([[G[0]], S[:N-1]])
    else:
        S = (Kp - Km).sum(axis=1)
    K = (Kp + Km)*np.outer(u, u) - np.diag(u*u*S)
    uc = np.where(n < Nc, u, 0.0)
    return K - 2.0*mus*np.outer(uc, uc)


def lam1(lams, wbars, mus, Nc, N, exact_S=True):
    M = build(lams, wbars, mus, Nc, N, exact_S)
    return float(eigh(M, subset_by_index=[N-1, N-1], eigvals_only=True)[0])


def members(lam, wbarE, x, rho):
    w1 = wbarE/np.sqrt((1.0-x) + x*rho*rho)
    return np.array([(1.0-x)*lam, x*lam]), np.array([w1, rho*w1])


def einstein_root(lam, mus, r, N, exact_S=True, lo=0.05, hi=14.0, ng=40, nb=50):
    def f(w):
        return lam1([lam], [w], mus, int(np.floor(r*w)), N, exact_S) - 1.0
    ws = np.linspace(lo, hi, ng); v = np.array([f(w) for w in ws])
    idx = [t for t in range(ng-1) if np.sign(v[t]) != np.sign(v[t+1])]
    if not idx:
        return None
    a, b, fa = ws[idx[0]], ws[idx[0]+1], v[idx[0]]
    for _ in range(nb):
        m = 0.5*(a+b); fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5*(a+b)


def main():
    N, r = 300, 10.0
    rng = np.random.default_rng(20260923)

    # ───────────────── GATE I：包络迁移 ─────────────────
    print("=" * 100)
    print("GATE I  两原子核的幂律包络： G(k) <= M2/(k^2+w_min^2) <= M2/k^2,  M2 = sum_j lam_j w_j^2")
    print("        （M2 在固定 (lam,wbar2) 的整族上**是常数** ⟹ 包络的 C 族内不变）")
    worst_a, worst_0, nviol, M2spread = -1e9, -1e9, 0, []
    for _ in range(4000):
        na = int(rng.integers(2, 5))
        L = rng.uniform(0.02, 3.0, na)
        W = np.exp(rng.uniform(np.log(0.02), np.log(30.0), na))
        M2 = float(L @ (W*W)); wmin = float(W.min())
        M2spread.append(M2)
        kk = np.arange(1, 4001)
        G = np.zeros_like(kk, float)
        for a, b in zip(L, W):
            G += a*b*b/(b*b + kk*kk)
        worst_a = max(worst_a, float(np.max(G - M2/(kk*kk + wmin*wmin))))
        worst_0 = max(worst_0, float(np.max(G - M2/(kk*kk))))
        nviol += (G > M2/(kk*kk + wmin*wmin) + 1e-12).any()
    print(f"        4000 组随机多原子 x 4000 个 k：违反 A=w_min 版 {nviol}/4000，"
          f"最劣越界 {worst_a:+.3e}；A=0 版最劣越界 {worst_0:+.3e}")
    # 族内 M2 不变性（定向：固定 (lam, wbar2) 的两原子族）
    lamT, wT = 3.2953, 0.7
    M2s = [float(members(lamT, wT, x, rho)[0] @ (members(lamT, wT, x, rho)[1]**2))
           for x in np.linspace(0.02, 0.98, 13) for rho in np.geomspace(1.0, 40.0, 13)]
    print(f"        族内 M2：min={min(M2s):.12f} max={max(M2s):.12f} "
          f"相对散布={(max(M2s)-min(M2s))/max(M2s):.3e}  (理论值 lam*wbar2^2={lamT*wT*wT:.12f})")
    print()

    # ───────────────── GATE II：S_n 闭式 vs N 截断行和 ─────────────────
    print("=" * 100)
    print("GATE II  S_n 闭式（精确，无截断） vs 我此前脚本的 N x N 块行和")
    for lam, mus in ((2.0, 0.13), (3.2953, 0.13), (6.0, 0.20)):
        wE_old = einstein_root(lam, mus, r, N, exact_S=False)
        wE_new = einstein_root(lam, mus, r, N, exact_S=True)
        Nc_o, Nc_n = int(np.floor(r*wE_old)), int(np.floor(r*wE_new))
        # 同一成员上两口径的 lam1 差
        d = lam1([lam], [wE_new], mus, Nc_n, N, False) - 1.0
        print(f"  lam={lam:<7.4f} mu*={mus}:  wbarE*(截断行和)={wE_old:.6f}  "
              f"wbarE*(闭式)={wE_new:.6f}  相对差={abs(wE_new-wE_old)/wE_new:.3e}  "
              f"(N_c {Nc_o}->{Nc_n});  同点 lam1 差(截断-闭式)={d:+.3e}")
    print()

    # ───────────────── GATE II-b：用**精确 S_n** 重跑 (A6) 族内 sup ─────────────────
    print("=" * 100)
    print("GATE II-b  用**精确 S_n** 重跑 (A6)：极点族上 Einstein 还是不是 sup？")
    from itertools import product
    xs = np.linspace(0.01, 0.99, 25)
    rhos = np.concatenate([[1.0], np.geomspace(1.005, 60.0, 24)])
    for mus in (0.0, 0.13, 0.20):
        for lam in (2.0, 3.2953, 6.0):
            wE = einstein_root(lam, mus, r, N, exact_S=True)
            Nc = int(np.floor(r*wE))
            l1E = lam1([lam], [wE], mus, Nc, N, True)
            best, bat = -1e9, None
            for x, rho in product(xs, rhos):
                L, W = members(lam, wE, x, rho)
                v = lam1(L, W, mus, Nc, N, True)
                if v > best:
                    best, bat = v, (x, rho)
            flag = "Einstein 是 sup" if best <= l1E + 1e-12 else "** Einstein 不是 sup **"
            print(f"  mu*={mus:<5} lam={lam:<7.4f}  lam1_E={l1E:.12f}  max={best:.12f}  "
                  f"超出={best-l1E:+.2e} at (x={bat[0]:.3f},rho={bat[1]:.3f})  [{len(xs)*len(rhos)}]  {flag}")
    print()

    # ───────────────── GATE III：预算——lam1 的相对松 -> lam_ceil 的位移 ─────────────────
    print("=" * 100)
    print("GATE III  预算：认证上界的相对松 delta（lam1 <= (1+delta)*lam1_true）")
    print("          -> 室温阈值 lam_ceil 位移多少？  消费点 = 300K 阈值")
    HB, KB = 64.65415134, 8.617333262e-2
    etaM = 9.172           # H3S 锚
    Ttar = 300.0
    wbar2 = HB*np.sqrt(etaM/3.2953)   # 占位，下面按 lam 自洽重算
    def Tc_of(lam, mus, delta):
        """给定 lam，解 (1+delta)*lam1(wbar)=1 得 wbar*，再换算 Tc（固定 eta/M）。"""
        def f(w):
            return (1.0+delta)*lam1([lam], [w], mus, int(np.floor(r*w)), N, True) - 1.0
        ws = np.linspace(0.05, 14.0, 40); v = np.array([f(w) for w in ws])
        idx = [t for t in range(39) if np.sign(v[t]) != np.sign(v[t+1])]
        if not idx:
            return None
        a, b, fa = ws[idx[0]], ws[idx[0]+1], v[idx[0]]
        for _ in range(50):
            m = 0.5*(a+b); fm = f(m)
            if (fm > 0) == (fa > 0):
                a, fa = m, fm
            else:
                b = m
        wstar = 0.5*(a+b)
        w2 = HB*np.sqrt(etaM/lam)      # wbar2 in meV，由固定 eta/M = lam*<w^2>/HB^2
        return w2/(2.0*np.pi*wstar)/KB  # Tc in K
    mus = 0.13
    print(f"  {'delta':>10}{'lam_ceil(300K)':>18}{'相对位移':>14}")
    base = None
    for delta in (0.0, 1e-4, 1e-3, 4.2e-3, 1e-2, 3e-2):
        lo, hi = 1.0, 12.0
        for _ in range(45):
            mid = 0.5*(lo+hi)
            t = Tc_of(mid, mus, delta)
            if t is None or t < Ttar:
                lo = mid
            else:
                hi = mid
        lc = 0.5*(lo+hi)
        if base is None:
            base = lc
        print(f"  {delta:>10.1e}{lc:>18.4f}{(lc-base)/base:>13.2%}")
    print("\n  读法：delta 是认证上界对 lam1 的相对松；lam_ceil 越小越强。")
    print("        delta 使 lam1 偏大 ⟹ 天花板偏松 ⟹ lam_ceil 偏小 ⟹ **方向是 fail-open**，")
    print("        必须用 lam1 的**上界**（KT）而不是下界，本表只量杠杆大小。")


if __name__ == "__main__":
    main()
