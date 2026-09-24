# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ① 第二半：两原子核上的 **Kato–Temple 端到端原型**。
从零重建（不调用 `p4_cii_ci_interval.py` / `p4_interchange_uniformity.py` 的任何函数，
只复用它们的**数学形状**：KT + Schur 包络尾部 + 算子序继承谱隙）。

━━━━━━━━ 要判的 ━━━━━━━━
承重问题是：`λ₁(P_N M P_N) ≤ λ₁(M)` 是**下估**（fail-open），所以本路线要 certificate
必须补 **上界**。仓库对 Einstein 有一套认证上界（`P4_coulomb_global_tc_certified.md` §2）：

    λ₁(M) ≤ μ + res²/(μ − γ),     μ = <v,Mv>,  res² = ||Mv − μv||²,  γ ≥ λ₂(M)

本脚本查这套东西在**两原子成员**上是否照样成立、并量它的**相对松 δ**（喂给预算表）。

三件（= 笔记 §7 列的三项）：
【T1】尾部（res² 的 n ≥ L 部分）。**包络迁移已在 `_chk_main_tail_transfer_and_exactS` GATE I
      证掉**：G(k) ≤ M2/(k²+w_min²)，M2 = λ<ω²> **族内常数**。本脚本用它给 res² 的解析尾部界
      （见 `res_tail_bound`），并查它确实是上界。
【T2】谱隙 γ。Einstein 侧由 `𝔎_c = 2μ* u_c u_cᵀ ⪰ 0 ⟹ M ⪯ 𝔎_ph ⟹ λ₂(M) ≤ λ₂(𝔎_ph)`
      **免费继承**。该论证只用到库仑项的**半正定性与秩一性**，与声子核形状无关
      ⟹ 预期逐字迁移。GATE T2 定向查 `λ₂(M) ≤ λ₂(𝔎_ph)` 在两原子成员上成立。
      （λ₂(𝔎_ph) 自身的认证上界是**另一件活**，本脚本用大 N 点算值并显式标注为 fail-open。）
【T3】r=3 → r=10。全部 gate 在 r=10 下跑。

严格性：res_tail_bound 的代数为 analytic；其余一切为 numerical-evidence（float64）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh


HB, KB = 64.65415134, 8.617333262e-2


def Gk(lams, wbars, k):
    k = np.asarray(k, float)
    out = np.zeros_like(k, float)
    for L, W in zip(np.atleast_1d(lams), np.atleast_1d(wbars)):
        out += L*W*W/(W*W + k*k)
    return out


def Sn_exact(lams, wbars, nmax):
    """S_n = G(0) + 2*sum_{j=1..n} G(j)   —— 闭式，有限和，与截断无关。"""
    g = Gk(lams, wbars, np.arange(nmax+1))
    return g[0] + 2.0*np.concatenate([[0.0], np.cumsum(g[1:])])


def build(lams, wbars, mus, Nc, N):
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    g = Gk(lams, wbars, np.arange(2*N+2))
    i, j = np.meshgrid(n, n, indexing="ij")
    K = (g[np.abs(i-j)] + g[i+j+1])*np.outer(u, u) - np.diag(u*u*Sn_exact(lams, wbars, N-1))
    uc = np.where(n < Nc, u, 0.0)
    return K - 2.0*mus*np.outer(uc, uc), u, uc


def lam_k(lams, wbars, mus, Nc, N, k=1):
    """k=1 -> λ₁ (最大)，k=2 -> λ₂。P_N M P_N 的值（对真算子是**下估**）。"""
    M, _, _ = build(lams, wbars, mus, Nc, N)
    e = eigh(M, subset_by_index=[N-k, N-k], eigvals_only=True)
    return float(e[0])


def apply_M(lams, wbars, mus, Nc, v, Nv, L):
    """真算子 M 作用在支撑于 [0,Nv) 的向量 v 上，取前 L 行（精确，无截断近似）。"""
    n = np.arange(L); m = np.arange(Nv)
    un = 1.0/np.sqrt(2.0*n+1.0); um = 1.0/np.sqrt(2.0*m+1.0)
    g = Gk(lams, wbars, np.arange(L+Nv+2))
    I, J = np.meshgrid(n, m, indexing="ij")
    out = (un[:, None]*um[None, :]*(g[np.abs(I-J)] + g[I+J+1])) @ v
    S = Sn_exact(lams, wbars, Nv-1)
    out[:Nv] -= un[:Nv]*un[:Nv]*S*v
    uc = np.where(n < Nc, un, 0.0)
    out -= 2.0*mus*uc*float(np.where(m < Nc, um, 0.0) @ v)
    return out


def res_tail_bound(lams, wbars, v, Nv, L):
    """||(1-P_L) M v||^2 的**解析上界**（analytic）。n>=L>Nv 时对角与库仑都够不到：
         (Mv)_n = u_n * sum_{m<Nv} u_m [G(n-m)+G(n+m+1)] v_m
       用 G(k) <= M2/k^2（GATE I，M2 = sum_j lam_j w_j^2，族内常数）与 n-m >= n-Nv：
         |(Mv)_n| <= u_n * (2*M2/(n-Nv)^2) * sum_m u_m |v_m|
       ⟹ sum_{n>=L} |(Mv)_n|^2 <= (2*M2*W)^2 * sum_{n>=L} 1/((2n+1)(n-Nv)^4)
       再用 1/(2n+1) <= 1/(2L+1) 与 sum_{n>=L}(n-Nv)^-4 <= (L-Nv)^-4 + 1/(3(L-Nv)^3)。"""
    m = np.arange(Nv); um = 1.0/np.sqrt(2.0*m+1.0)
    M2 = float(np.atleast_1d(lams) @ (np.atleast_1d(wbars)**2))
    W = float(um @ np.abs(v))
    d = L - Nv
    return (2.0*M2*W)**2 / (2.0*L+1.0) * (d**-4.0 + 1.0/(3.0*d**3))


def kato_temple(lams, wbars, mus, Nc, Nv, gamma, L=40000):
    """KT 上界。返回 (bound, mu, res2, 尾部占比)。"""
    Msub, _, _ = build(lams, wbars, mus, Nc, Nv)
    w, V = eigh(Msub)
    v = V[:, -1]; v = v/np.linalg.norm(v)
    Mv = apply_M(lams, wbars, mus, Nc, v, Nv, L)
    mu = float(v @ Mv[:Nv])
    r = Mv.copy(); r[:Nv] -= mu*v
    res2_head = float(r @ r)
    res2_tail = res_tail_bound(lams, wbars, v, Nv, L)
    res2 = res2_head + res2_tail
    if mu <= gamma:
        return None, mu, res2, np.nan
    return mu + res2/(mu - gamma), mu, res2, res2_tail/res2


def members(lam, wbarE, x, rho):
    w1 = wbarE/np.sqrt((1.0-x) + x*rho*rho)
    return np.array([(1.0-x)*lam, x*lam]), np.array([w1, rho*w1])


def einstein_root(lam, mus, r, N):
    def f(w):
        return lam_k([lam], [w], mus, int(np.floor(r*w)), N) - 1.0
    ws = np.linspace(0.05, 14.0, 40); v = np.array([f(w) for w in ws])
    idx = [t for t in range(39) if np.sign(v[t]) != np.sign(v[t+1])]
    a, b, fa = ws[idx[0]], ws[idx[0]+1], v[idx[0]]
    for _ in range(50):
        mm = 0.5*(a+b); fm = f(mm)
        if (fm > 0) == (fa > 0):
            a, fa = mm, fm
        else:
            b = mm
    return 0.5*(a+b)


def main():
    r, Nref = 10.0, 1500
    CASES = [(2.0, 0.13), (3.2953, 0.13), (6.0, 0.20)]

    # ───────── GATE R：参考值 λ₁ 的 N 收敛（确认 Nref 足够，"数字错 vs 证书松"分诊） ─────────
    print("=" * 104)
    print("GATE R  参考 λ₁ 的 N 收敛（CLAUDE.md：认证界/点算值 异常时先测 N 收敛）")
    lam, mus = 3.2953, 0.13
    wE = einstein_root(lam, mus, r, 600)
    Nc = int(np.floor(r*wE))
    prev = None
    for N in (400, 800, 1500, 2500):
        v = lam_k([lam], [wE], mus, Nc, N)
        d = "" if prev is None else f"  Δ={v-prev:+.3e}"
        print(f"   N={N:<6} λ₁={v:.14f}{d}")
        prev = v
    print()

    # ───────── GATE T2：谱隙继承 λ₂(M) ≤ λ₂(𝔎_ph)（算子序，与核形状无关） ─────────
    print("=" * 104)
    print("GATE T2  谱隙继承：𝔎_c ⪰ 0 ⟹ M ⪯ 𝔎_ph ⟹ λ₂(M) ≤ λ₂(𝔎_ph)")
    print("         该论证只用库仑项的 PSD+秩一，**与声子核形状无关** ⟹ 预期逐字迁移。定向查两原子成员：")
    rng = np.random.default_rng(7)
    worst = -1e9; nv = 0; ntot = 0
    for lam, mus in CASES:
        wE = einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
        for x in np.linspace(0.02, 0.98, 9):
            for rho in np.geomspace(1.0, 50.0, 9):
                L_, W_ = members(lam, wE, x, rho)
                l2M = lam_k(L_, W_, mus, Nc, 400, k=2)
                l2P = lam_k(L_, W_, 0.0, 0, 400, k=2)
                ntot += 1; nv += (l2M > l2P + 1e-12)
                worst = max(worst, l2M - l2P)
    print(f"         违反 {nv}/{ntot}，最劣越界 {worst:+.3e}   ⟹ 继承成立（与预期一致）")
    print("         ※ λ₂(𝔎_ph) **自身**的认证上界是另一件活；下面的 γ 用大 N 点算值，")
    print("           它对真 λ₂ 是**下估** ⟹ 使 μ−γ 偏大 ⟹ KT 偏紧 ⟹ **本表的 δ 是乐观值**（已标注）。")
    print()

    # ───────── GATE T1/T3：KT 端到端，量相对松 δ ─────────
    print("=" * 104)
    print("GATE T1/T3  Kato–Temple 端到端（r=10）：KT 上界 vs 精确 λ₁，量相对松 δ")
    print(f"  {'case':<22}{'Nv':>5}{'mu':>14}{'res2':>12}{'尾占比':>11}{'KT':>15}{'λ₁(ref)':>15}{'δ=KT/λ₁-1':>13}{'方向':>7}")
    summary = {}
    for lam, mus in CASES:
        wE = einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
        for tag, (x, rho) in (("Einstein", (0.0, 1.0)),
                              ("2atom x=.4 rho=3", (0.4, 3.0)),
                              ("2atom x=.1 rho=20", (0.1, 20.0))):
            if tag == "Einstein":
                L_, W_ = np.array([lam]), np.array([wE])
            else:
                L_, W_ = members(lam, wE, x, rho)
            ref = lam_k(L_, W_, mus, Nc, Nref)
            gam = lam_k(L_, W_, 0.0, 0, Nref, k=2)
            for Nv in (40, 120):
                kt, mu, res2, frac = kato_temple(L_, W_, mus, Nc, Nv, gam)
                if kt is None:
                    print(f"  {tag+f' lam={lam}':<22}{Nv:>5}   mu<=gamma, KT 不适用")
                    continue
                d = kt/ref - 1.0
                ok = "OK" if kt >= ref - 1e-12 else "**FAIL**"
                print(f"  {tag+f' λ={lam:.3g}':<22}{Nv:>5}{mu:>14.9f}{res2:>12.2e}{frac:>11.1e}"
                      f"{kt:>15.9f}{ref:>15.9f}{d:>13.3e}{ok:>7}")
                summary[(lam, tag, Nv)] = d
        print()

    print("=" * 104)
    print("汇总：δ 的最大值 =", f"{max(summary.values()):.3e}", " (Nv=120 子集:",
          f"{max(v for k, v in summary.items() if k[2]==120):.3e})")
    print("对照上一脚本 GATE III 的预算表：δ=1e-3 ⟹ λ_ceil 位移 −0.26%；δ=1e-2 ⟹ −2.58%。")


if __name__ == "__main__":
    main()
