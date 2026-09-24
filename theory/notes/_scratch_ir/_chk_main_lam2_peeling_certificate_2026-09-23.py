# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ②(b)：`λ₂(𝔎_ph)` 的**认证上界**原型 —— 从零重建。

━━━━━━━━ 为什么做这件、预算是多少（CLAUDE.md「预算优先」） ━━━━━━━━
KT 上界 `λ₁ ≤ μ + res²/(μ−γ)` 需要一个 **γ ≥ λ₂(M)**。原型里 γ 用的是大 N 点算 λ₂，
那是**下估** ⟹ fail-open，是端到端链条里唯一剩下的不声音件。
`_chk_main_kt_gamma_sensitivity_2026-09-23.py` 已量过预算：把 γ 朝判据变难的方向退化
`γ → γ + s(μ−γ)`，`s=0.99`（即 γ≈0.99、μ≈1）时 δ 仍 ≤1.1e−04
⟹ **消费点只要求 `γ_cert ≤ ~0.99`；准到一位有效数字即可。**
（真值 λ₂≈7e−04，谱隙 μ/λ₂≈1400 ⟹ 允许 3 个数量级的松。）

━━━━━━━━ 工具：Cauchy 交错 + 加权 Schur 检验（= 仓库「引理 Λ」的形状） ━━━━━━━━
写 `𝔎_ph = B − D`，
    B_ij = u_i u_j [G(|i−j|) + G(i+j+1)]  ≥ 0  （**含对角**，因 G ≥ 0）
    D_ii = u_i² S_i,  S_i = G(0) + 2Σ_{j=1..i} G(j)   （闭式，有限和）
两步都是单边、方向朝「判据变难」：
  (Λ-1) **Cauchy 交错**：去掉任意一个指标 p，主子阵 A' 满足 `λ₂(A) ≤ λ₁(A')`。
  (Λ-2) **加权 Schur / Collatz–Wielandt**：B 逐元非负 ⟹ 对任意正权 w，
        `λ₁(B'−D') ≤ sup_{i≠p} [ (1/w_i) Σ_{j≠p} B_ij w_j − D_ii ] =: R(p, w)`。
⟹ `λ₂(𝔎_ph) ≤ R(p,w)` 对**任意** p 与任意正 w 成立。

**估计子的固有松 = 0**（CLAUDE.md「目标常数必须可达」的正面用法）：B−D 是 Metzler
（非对角非负）且不可约 ⟹ Perron–Frobenius ⟹ 其最大特征向量**为正**；取 w = 该向量时
(Λ-2) **取等**。所以 R 的唯一误差来源是 w 的数值质量与尾部处理，不是估计子本身。

━━━━━━━━ 本脚本的四道门 ━━━━━━━━
【Λ0】固有松实测：R(p,w) / λ₂_true —— 在派任何「收紧」活之前先知道这个数。
【Λ1】权族比较：w=u_j、w=1、w=u^q、w=Perron（数值）——哪个够用。
【Λ2】**方向门**：R ≥ λ₂ 必须处处成立（单边界方向规则；出现 R < λ₂ 即是 bug）。
【Λ3】N 收敛 + 两原子族一致性（族内 sup R 是否仍 ≪ 0.99）。

严格性：(Λ-1)(Λ-2) 的代数为 `analytic`；本脚本一切数值为 `numerical-evidence`（float64）。
**尾部（j>N 与 i>N）本轮未认证**，故本脚本产出的不是证书，只是「目标可达」的实测。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh


def Gk(lams, wbars, k):
    k = np.asarray(k, float)
    out = np.zeros_like(k, float)
    for L, W in zip(np.atleast_1d(lams), np.atleast_1d(wbars)):
        out += L*W*W/(W*W + k*k)
    return out


def Sn_exact(lams, wbars, nmax):
    """S_n = G(0) + 2*sum_{j=1..n} G(j)  —— 闭式，有限和，与截断无关。"""
    g = Gk(lams, wbars, np.arange(nmax+1))
    return g[0] + 2.0*np.concatenate([[0.0], np.cumsum(g[1:])])


def BD(lams, wbars, N):
    """返回 (B, D_diag)：𝔎_ph = B − diag(D)，B 逐元非负（含对角）。"""
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    g = Gk(lams, wbars, np.arange(2*N+2))
    i, j = np.meshgrid(n, n, indexing="ij")
    B = (g[np.abs(i-j)] + g[i+j+1])*np.outer(u, u)
    D = u*u*Sn_exact(lams, wbars, N-1)
    return B, D


def kph(lams, wbars, N):
    B, D = BD(lams, wbars, N)
    return B - np.diag(D)


def lam_k(A, k=1):
    N = A.shape[0]
    return float(eigh(A, subset_by_index=[N-k, N-k], eigvals_only=True)[0])


def R_of(lams, wbars, N, p, w):
    """R(p,w) = sup_{i≠p} [ (1/w_i) Σ_{j≠p} B_ij w_j − D_ii ]，正权 w（长度 N）。"""
    B, D = BD(lams, wbars, N)
    idx = np.array([t for t in range(N) if t != p])
    Bs = B[np.ix_(idx, idx)]; Ds = D[idx]; ws = np.asarray(w, float)[idx]
    row = (Bs @ ws)/ws - Ds
    return float(row.max()), idx[int(np.argmax(row))]


def weights(kind, N, lams=None, wbars=None, p=0, q=0.5):
    n = np.arange(N)
    if kind == "u":
        return 1.0/np.sqrt(2.0*n+1.0)
    if kind == "one":
        return np.ones(N)
    if kind == "u^q":
        return (2.0*n+1.0)**(-q)
    if kind == "perron":
        A = kph(lams, wbars, N)
        idx = np.array([t for t in range(N) if t != p])
        _, V = eigh(A[np.ix_(idx, idx)], subset_by_index=[N-2, N-2])
        v = np.abs(V[:, 0].ravel())
        v = np.maximum(v, v.max()*1e-14)
        out = np.empty(N); out[idx] = v; out[p] = v.max()
        return out
    raise ValueError(kind)


def members(lam, wbarE, x, rho):
    w1 = wbarE/np.sqrt((1.0-x) + x*rho*rho)
    return np.array([(1.0-x)*lam, x*lam]), np.array([w1, rho*w1])


def einstein_root(lam, mus, r, N):
    """解 λ₁(M)=1 的 ω̄_E（与前两个脚本同口径：闭式 S_n）。"""
    def f(w):
        n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
        Nc = int(np.floor(r*w))
        uc = np.where(n < Nc, u, 0.0)
        return lam_k(kph([lam], [w], N) - 2.0*mus*np.outer(uc, uc)) - 1.0
    ws = np.linspace(0.05, 14.0, 40); v = np.array([f(t) for t in ws])
    idx = [t for t in range(39) if np.sign(v[t]) != np.sign(v[t+1])]
    a, b, fa = ws[idx[0]], ws[idx[0]+1], v[idx[0]]
    for _ in range(50):
        m = 0.5*(a+b); fm = f(m)
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5*(a+b)


def main():
    r = 10.0
    CASES = [(2.0, 0.13), (3.2953, 0.13), (6.0, 0.20)]
    BAR = 0.99          # 预算给的消费点上限（见文件头）

    print("="*110)
    print("GATE Λ0/Λ1  估计子固有松实测：R(p,w) vs 真 λ₂(𝔎_ph)，以及 vs 预算上限 0.99")
    print("            （Cauchy 交错去掉指标 p=0，再加权 Schur；两步都是 fail-closed 方向）")
    print(f"  {'case':<20}{'N':>6}{'λ₁(𝔎_ph)':>13}{'λ₂(𝔎_ph)':>13}"
          f"{'R[w=u]':>12}{'R[w=1]':>12}{'R[w=u^.25]':>13}{'R[Perron]':>12}{'R/λ₂':>11}{'≤0.99?':>9}")
    N = 800
    store = {}
    for lam, mus in CASES:
        wE = einstein_root(lam, mus, r, 400)
        A = kph([lam], [wE], N)
        l1, l2 = lam_k(A, 1), lam_k(A, 2)
        Rs = {}
        for kind, kw in (("u", {}), ("one", {}), ("u^q", {"q": 0.25}), ("perron", {})):
            w = weights(kind, N, [lam], [wE], 0, **kw)
            Rs[kind], _ = R_of([lam], [wE], N, 0, w)
        best = min(Rs.values())
        store[(lam, mus)] = (wE, l1, l2, Rs)
        print(f"  λ={lam:<6.4g} μ*={mus:<5}{N:>6}{l1:>13.6f}{l2:>13.3e}"
              f"{Rs['u']:>12.4f}{Rs['one']:>12.4f}{Rs['u^q']:>13.4f}{Rs['perron']:>12.3e}"
              f"{best/l2:>11.3f}{'YES' if best <= BAR else '**NO**':>9}")
    print()

    print("="*110)
    print("GATE Λ2  方向门：R(p,w) ≥ λ₂ 必须处处成立（单边界方向规则；R<λ₂ 即 bug）")
    nv = ntot = 0; worst = 1e9
    rng = np.random.default_rng(23)
    for lam, mus in CASES:
        wE = store[(lam, mus)][0]
        for x in np.linspace(0.02, 0.98, 7):
            for rho in np.geomspace(1.0, 40.0, 7):
                L_, W_ = members(lam, wE, x, rho)
                A = kph(L_, W_, 400)
                l2 = lam_k(A, 2)
                for kind, kw in (("u", {}), ("perron", {})):
                    w = weights(kind, 400, L_, W_, 0, **kw)
                    R, _ = R_of(L_, W_, 400, 0, w)
                    ntot += 1; nv += (R < l2 - 1e-12)
                    worst = min(worst, R - l2)
    print(f"          违反 {nv}/{ntot}，最劣余量 min(R−λ₂) = {worst:+.3e}")
    print()

    print("="*110)
    print("GATE Λ3  两原子族一致性：sup_族 R 是否仍 ≪ 0.99？（w=u 与 Perron 两族权）")
    print(f"  {'case':<20}{'sup R[w=u]':>14}{'sup R[Perron]':>16}{'sup λ₂':>12}{'≤0.99?':>9}{'argmax(x,ρ)':>18}")
    for lam, mus in CASES:
        wE = store[(lam, mus)][0]
        sR_u = sR_p = -1e9; sl2 = -1e9; at = None
        for x in np.linspace(0.02, 0.98, 11):
            for rho in np.geomspace(1.0, 40.0, 11):
                L_, W_ = members(lam, wE, x, rho)
                A = kph(L_, W_, 400); sl2 = max(sl2, lam_k(A, 2))
                Ru, _ = R_of(L_, W_, 400, 0, weights("u", 400))
                Rp, _ = R_of(L_, W_, 400, 0, weights("perron", 400, L_, W_, 0))
                if Ru > sR_u:
                    sR_u, at = Ru, (x, rho)
                sR_p = max(sR_p, Rp)
        ok = "YES" if min(sR_u, sR_p) <= BAR else "**NO**"
        print(f"  λ={lam:<6.4g} μ*={mus:<5}{sR_u:>14.4f}{sR_p:>16.3e}{sl2:>12.3e}{ok:>9}"
              f"{f'({at[0]:.2f},{at[1]:.2f})':>18}")
    print()

    print("="*110)
    print("GATE Λ4  N 收敛（R 与 λ₂ 各自随 N）——分诊「数字错 vs 证书松」")
    lam, mus = 3.2953, 0.13
    wE = store[(lam, mus)][0]
    prev = None
    for N_ in (200, 400, 800, 1600):
        A = kph([lam], [wE], N_)
        l2 = lam_k(A, 2)
        Ru, iu = R_of([lam], [wE], N_, 0, weights("u", N_))
        Rp, ip = R_of([lam], [wE], N_, 0, weights("perron", N_, [lam], [wE], 0))
        d = "" if prev is None else f"  ΔR[u]={Ru-prev:+.2e}"
        print(f"   N={N_:<6} λ₂={l2:.10e}   R[w=u]={Ru:.8f} (argmax i={iu})"
              f"   R[Perron]={Rp:.6e} (i={ip}){d}")
        prev = Ru
    print()
    print("  读法：两步都是上界方向 ⟹ 本表只能产出「γ_cert 够小 ⟹ 可用」的 PASS 判决；")
    print("        若某格 R>0.99，那**不是** FAIL，只说明该权 w 选得不好（估计子固有松为 0）。")


if __name__ == "__main__":
    main()
