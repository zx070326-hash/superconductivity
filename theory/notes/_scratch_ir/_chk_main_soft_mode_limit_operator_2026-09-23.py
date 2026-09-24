# -*- coding: utf-8 -*-
"""②-a：`t=1/ρ → 0` 的极限算子 —— 并检验 §1 (A3) 紧性论证的一个漏洞。

ρ→∞（x 固定）时 ω₁ = ω̄_E/√((1−x)+xρ²) → 0，ω₂ = ρω₁ → ω̄_E/√x。
**核**有明确极限：G₁(k) = (1−x)λ ω₁²/(ω₁²+k²) → (1−x)λ·[k=0]，即
    G_∞(k) = (1−x)λ·δ_{k,0} + xλ ω₂²/(ω₂²+k²),   ω₂ = ω̄_E/√x
**但测度没有**：ν₁ = (λ₁ω₁/2)δ_{ω₁}，质量 → 0（弱*收敛到 0），
而 `∫ω⁻¹dν` 却保住 (1−x)λ/2（ω⁻¹ 在 0 附近无界，不是弱*连续泛函）。
⟹ 𝒞 = {∫ω⁻¹dν=λ/2, ∫ωdν=λ⟨ω²⟩/2} **弱*预紧但不闭**；λ₁ 也不弱*连续
（软模的质量没了，它的 G(0) 贡献还在）。

本脚本查三件：
【L1】λ₁(ρ→∞) 是否收敛到**显式极限算子** G_∞ 给的值（而**不是**单原子 (xλ, ω₂) 的值）。
      若是，则「ρ=∞ 端」是一个良定义的、可显式写下的边界成员 ⟹ ②-a 的紧化可做。
【L2】两条矩在 ρ→∞ 的行为：Σλ_j 守恒、Σλ_jω_j² 守恒、但 ν 的**总质量** Σλ_jω_j/2 掉多少。
【L3】极限成员的 λ₁ 是否仍 < λ₁_E（B1 给 0.998，此处用显式极限算子复核）。
numerical-evidence（float64）。
"""
import os; os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util, numpy as np
from scipy.linalg import eigh
p = "theory/notes/_scratch_ir/_chk_main_kato_temple_two_atom_prototype_2026-09-23.py"
s_ = importlib.util.spec_from_file_location("kt", p); m_ = importlib.util.module_from_spec(s_); s_.loader.exec_module(m_)


def lam1_from_G(g, mus, Nc, N):
    """直接用给定的 G 数组（长度 >= 2N+2）建 M 并取 λ₁，S_n 用闭式。"""
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    S = g[0] + 2.0*np.concatenate([[0.0], np.cumsum(g[1:N])])
    K = (g[np.abs(i-j)] + g[i+j+1])*np.outer(u, u) - np.diag(u*u*S)
    uc = np.where(n < Nc, u, 0.0)
    return float(eigh(K - 2.0*mus*np.outer(uc, uc), subset_by_index=[N-1, N-1], eigvals_only=True)[0])


r, N = 10.0, 500
for lam, mus in ((3.2953, 0.13), (6.0, 0.20)):
    wE = m_.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    l1E = m_.lam_k([lam], [wE], mus, Nc, N)
    print("=" * 98)
    print(f"λ={lam}  μ*={mus}  ω̄_E={wE:.6f}  N_c={Nc}  λ₁_E={l1E:.12f}")
    for x in (0.2, 0.5, 0.8):
        w2 = wE/np.sqrt(x)
        k = np.arange(2*N+2)
        gL = x*lam*w2*w2/(w2*w2 + k*k); gL[0] += (1.0-x)*lam        # 显式极限算子 G_∞
        l_lim = lam1_from_G(gL, mus, Nc, N)
        l_single = m_.lam_k([x*lam], [w2], mus, Nc, N)               # 对照：丢掉软模
        print(f"  x={x}:  极限算子 λ₁={l_lim:.9f}   单原子(xλ,ω₂) λ₁={l_single:.9f}   "
              f"差={l_lim-l_single:+.3e}")
        print(f"          ρ →     " + "".join(f"{v:>12}" for v in (20, 100, 1000, 10000)))
        vals = [m_.lam_k(*m_.members(lam, wE, x, rho), mus, Nc, N) for rho in (20, 100, 1000, 10000)]
        print(f"          λ₁      " + "".join(f"{v:>12.9f}" for v in vals))
        print(f"          −极限   " + "".join(f"{v-l_lim:>12.2e}" for v in vals)
              + f"    | 极限/λ₁_E = {l_lim/l1E:.6f}")
        # L2：矩与质量
        L_, W_ = m_.members(lam, wE, x, 10000.0)
        print(f"          Σλ_j={L_.sum():.9f}(应={lam})  Σλ_jω_j²={float(L_@(W_**2)):.9f}"
              f"(应={lam*wE*wE:.9f})  ν 总质量 Σλ_jω_j/2={float(L_@W_)/2:.6f}"
              f"  (ρ=1 时={lam*wE/2:.6f})")


# ───────────────────── 【L4】G(0) 在 M 里**恒等地消失**（投毒门） ─────────────────────
# 代数（analytic，一行）：
#   diag(K)_ii = u_i²[G(0)+G(2i+1)] − u_i²S_i，而 S_i = G(0) + 2Σ_{j=1}^{i}G(j)
#            ⟹ = u_i²[G(2i+1) − 2Σ_{j=1}^{i}G(j)]          ← G(0) 精确抵消
#   offdiag (i≠j) = u_iu_j[G(|i−j|)+G(i+j+1)]，|i−j|≥1 且 i+j+1≥1 ⟹ G(0) 从不出现
#   ⟹ **M 完全不依赖 G(0)**；只有 G(k≥1) 进入。
# 推论：ω ≪ 2πT 的软模对 λ₁（从而 Tc）贡献**恰好为零**，因为它的全部权重都塌在 G(0) 上。
# 这同时修好 §1 (A3) 的紧性漏洞：`∫ω⁻¹dν` 在 ω→0 处不是弱*连续泛函，但**那一端根本不进算子**。
print()
print("=" * 98)
print("L4  投毒门：把 G(0) 任意改动，λ₁ 应当**一点不动**（若动了，上面的代数就是错的）")
for lam, mus in ((3.2953, 0.13), (6.0, 0.20)):
    wE = m_.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    k = np.arange(2*N+2)
    g0 = lam*wE*wE/(wE*wE + k*k)
    base = lam1_from_G(g0.copy(), mus, Nc, N)
    row = []
    for bump in (0.0, 1.0, 1e3, -1e3, 1e6):
        g = g0.copy(); g[0] = bump
        row.append(lam1_from_G(g, mus, Nc, N) - base)
    print(f"  λ={lam:<8} λ₁={base:.12f}   G(0) 置为 {{0,1,1e3,−1e3,1e6}} 时 Δλ₁ = "
          + ", ".join(f"{v:+.1e}" for v in row))
print("  （G(0) 的真值 = Σλ_j = λ；此门说明该值本身不承重，承重的只有 G(k≥1)）")
