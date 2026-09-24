# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ②(b) 第二半：把 `λ₂` 的剥指标界**接进 KT**，并补掉两处尾部。

前一脚本 `_chk_main_lam2_peeling_certificate_2026-09-23.py` 的结论：
  λ₂(M) ≤ λ₂(𝔎_ph) ≤ R(p=0, w=u) ≈ 0.097–0.113   （预算线 0.99）
其中 **w_j = u_j = (2j+1)^{−1/2} 是显式解析权**，不需要数值 Perron 向量。

本脚本做三件：
【K1】端到端：把 KT 的 γ 从「大 N 点算 λ₂」（fail-open）换成 `γ = R[w=u]`（fail-closed），
      重新量 δ = KT/λ₁ − 1。这是 CLAUDE.md「换了界就要重解」——不是把新界代进旧解，
      KT 的 v、μ、res² 都随 Nv 重算，只有 γ 这一项换掉。
【K2】**列尾**（j>J）：R_i 的无穷列和用闭式尾界
        Σ_{j>J} [G(|i−j|)+G(i+j+1)]/(2j+1) ≤ (2·M2/(2J+1))·(1/(J−i))      (J>i)
      —— 用 G(k) ≤ M2/k²（GATE I 已证的族内一致包络）。数值查它确是上界。
【K3】**行尾**（i>I）：手推的显式衰减界（只依赖两条矩 λ 与 M2 ⟹ **族内一致**）
        R_i ≤ T_i ≤ B_i := 4·M2·ln(i+1)/i²  +  2(λ + π²M2/6)/(i+1)      (i ≥ 2)
      推法（analytic）：令 h=⌊i/2⌋。
        · j ≤ h ⟹ |i−j| ≥ i/2 且 i+j+1 > i/2 ⟹ 两个 G 都 ≤ 4M2/i²；
          且 Σ_{j=1..h} 1/(2j+1) ≤ ½ln(2h+1) ≤ ½ln(i+1)。
        · j > h ⟹ 2j+1 ≥ 2h+3 ≥ i+1；且 Σ_j[G(|i−j|)+G(i+j+1)] ≤ 2Σ_{k≥0}G(k)
          ≤ 2(λ + π²M2/6)（G(0)=λ、Σ_{k≥1}G(k) ≤ M2·π²/6）。
        · 对角 −D_ii ≤ 0 **直接丢掉**（方向：丢负项使界变大 ⟹ 仍是上界，fail-closed）。
      数值查 B_i ≥ R_i，并给出使 B_I ≤ γ_target 的 I。

严格性：【K2】【K3】的不等式代数为 `analytic`；一切数值为 `numerical-evidence`（float64）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import numpy as np
from scipy.linalg import eigh

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")
PL = _load("_chk_main_lam2_peeling_certificate_2026-09-23")


def R_rows(lams, wbars, N, p=0):
    """返回 R_i 的逐行值（w=u，去掉指标 p），以及 (λ, M2)。"""
    B, D = PL.BD(lams, wbars, N)
    idx = np.array([t for t in range(N) if t != p])
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    Bs = B[np.ix_(idx, idx)]; Ds = D[idx]; ws = u[idx]
    return idx, (Bs @ ws)/ws - Ds


def main():
    r = 10.0
    CASES = [(2.0, 0.13), (3.2953, 0.13), (6.0, 0.20)]
    N_for_R = 1600

    # ───────── K1：把 γ 换成认证式的 R[w=u]，端到端重量 δ ─────────
    print("=" * 112)
    print("GATE K1  端到端：γ 从『点算 λ₂（fail-open）』换成『剥指标界 R[w=u]（fail-closed）』")
    print(f"  {'case':<24}{'Nv':>5}{'γ_point':>12}{'γ_cert=R[u]':>14}"
          f"{'KT(γ_point)':>14}{'KT(γ_cert)':>13}{'λ₁(ref)':>13}{'δ_cert':>11}{'方向':>8}")
    worst = -1.0
    for lam, mus in CASES:
        wE = KT.einstein_root(lam, mus, r, 600)
        Nc = int(np.floor(r*wE))
        for tag, (x, rho) in (("Einstein", (0.0, 1.0)),
                              ("2atom x=.4 ρ=3", (0.4, 3.0)),
                              ("2atom x=.1 ρ=20", (0.1, 20.0))):
            if tag == "Einstein":
                L_, W_ = np.array([lam]), np.array([wE])
            else:
                L_, W_ = KT.members(lam, wE, x, rho)
            ref = KT.lam_k(L_, W_, mus, Nc, 1500)
            g_pt = KT.lam_k(L_, W_, 0.0, 0, 1500, k=2)
            g_ct, _ = PL.R_of(L_, W_, N_for_R, 0, PL.weights("u", N_for_R))
            for Nv in (40, 120):
                kt_p, _, _, _ = KT.kato_temple(L_, W_, mus, Nc, Nv, g_pt)
                kt_c, mu, res2, _ = KT.kato_temple(L_, W_, mus, Nc, Nv, g_ct)
                if kt_c is None:
                    print(f"  {tag:<24}{Nv:>5}  μ={mu:.4f} ≤ γ_cert={g_ct:.4f}  ⟹ KT 不适用")
                    continue
                d = kt_c/ref - 1.0
                worst = max(worst, d)
                ok = "OK" if kt_c >= ref - 1e-12 else "**FAIL**"
                print(f"  {tag+f' λ={lam:.4g}':<24}{Nv:>5}{g_pt:>12.3e}{g_ct:>14.6f}"
                      f"{kt_p:>14.9f}{kt_c:>13.9f}{ref:>13.9f}{d:>11.3e}{ok:>8}")
        print()
    print(f"  汇总：δ_cert 最大值 = {worst:.3e}   （预算表：δ=1e−3 ⟹ λ_ceil −0.26%）")
    print()

    # ───────── K2：列尾界 ─────────
    print("=" * 112)
    print("GATE K2  列尾（j>J）：Σ_{j>J}[G(|i−j|)+G(i+j+1)]/(2j+1) ≤ (2·M2/(2J+1))/(J−i)")
    nv = ntot = 0; worstfac = 0.0
    for lam, mus in CASES:
        wE = KT.einstein_root(lam, mus, r, 600)
        for x, rho in ((0.0, 1.0), (0.4, 3.0), (0.1, 20.0), (0.8, 1.5)):
            L_, W_ = (np.array([lam]), np.array([wE])) if rho == 1.0 and x == 0.0 \
                else KT.members(lam, wE, x, rho)
            M2 = float(np.atleast_1d(L_) @ (np.atleast_1d(W_)**2))
            for J in (200, 500, 1000):
                jj = np.arange(J+1, 20000)
                for i in (0, 1, 5, 20, 100):
                    true = float(np.sum((KT.Gk(L_, W_, np.abs(i-jj))
                                         + KT.Gk(L_, W_, i+jj+1))/(2.0*jj+1.0)))
                    bnd = (2.0*M2/(2.0*J+1.0))/(J-i)
                    ntot += 1; nv += (true > bnd + 1e-18)
                    worstfac = max(worstfac, true/bnd)
    print(f"          违反 {nv}/{ntot}，最劣 true/bound = {worstfac:.4f}（<1 ⟹ 是上界）")
    print()

    # ───────── K3：行尾界 ─────────
    print("=" * 112)
    print("GATE K3  行尾（i>I）：R_i ≤ B_i := 4·M2·ln(i+1)/i² + 2(λ+π²M2/6)/(i+1)，i≥2")
    print("          （只含 λ 与 M2 两条矩 ⟹ **族内一致**；丢掉 −D_ii 是 fail-closed 方向）")
    print(f"  {'case':<26}{'λ':>8}{'M2':>10}{'max R_i':>11}{'max R_i/B_i':>14}"
          f"{'违反':>8}{'I: B_I≤0.15':>14}{'I: B_I≤0.99':>14}")
    for lam, mus in CASES:
        wE = KT.einstein_root(lam, mus, r, 600)
        for x, rho in ((0.0, 1.0), (0.4, 3.0), (0.1, 20.0)):
            L_, W_ = (np.array([lam]), np.array([wE])) if (x == 0.0 and rho == 1.0) \
                else KT.members(lam, wE, x, rho)
            M2 = float(np.atleast_1d(L_) @ (np.atleast_1d(W_)**2))
            lam_tot = float(np.sum(np.atleast_1d(L_)))
            idx, Ri = R_rows(L_, W_, N_for_R, 0)
            ii = idx.astype(float)
            m = idx >= 2
            Bi = 4.0*M2*np.log(ii[m]+1.0)/ii[m]**2 + 2.0*(lam_tot + np.pi**2*M2/6.0)/(ii[m]+1.0)
            ratio = Ri[m]/Bi
            bad = int(np.sum(Ri[m] > Bi + 1e-14))
            big = np.arange(2, 20000, dtype=float)
            Bbig = 4.0*M2*np.log(big+1.0)/big**2 + 2.0*(lam_tot + np.pi**2*M2/6.0)/(big+1.0)
            I15 = int(big[np.argmax(Bbig <= 0.15)]) if (Bbig <= 0.15).any() else -1
            I99 = int(big[np.argmax(Bbig <= 0.99)]) if (Bbig <= 0.99).any() else -1
            tag = f"λ={lam:.4g} x={x:.1f} ρ={rho:.0f}"
            print(f"  {tag:<26}{lam_tot:>8.4f}{M2:>10.4f}{Ri.max():>11.4f}"
                  f"{ratio.max():>14.4f}{bad:>8}{I15:>14}{I99:>14}")
    print()
    print("  读法：max R_i/B_i < 1 ⟹ B_i 确是上界；I 列给出『必须逐行精确算到多少』——")
    print("        i ≤ I 逐行精确（有限）、i > I 用 B_i ⟹ sup_i R_i 的完整包围只差区间算术。")


if __name__ == "__main__":
    main()
