# -*- coding: utf-8 -*-
"""任务 ① 之后的接续：步骤 ②（二维区间细分认证）**是否有限**？

① 已判为工程量 ⟹ 承重瓶颈搬到 ②。② 的可行性只看两件：
【B1】**紧化**：ρ∈[1,∞) 不紧。若 max_{ρ≥ρ*} λ₁ 随 ρ* 单调衰减并远离 1，则 ρ>ρ* 的整个尾区
      可被**一个粗界**一次盖住，硬区只剩紧致盒 (x,ρ)∈[0,1]×[1,ρ*]。
【B2】**ρ=1 附近的退化**：ρ=1 整条线上 λ₁ 与 x 无关（Einstein 简并）⟹ 不可能证严格不等式。
      但天花板**不需要**严格：目标是 `sup λ₁ ≤ λ₁_E + ε`，ε 的代价已由预算表给出
      （ε=1e−3 ⟹ λ_ceil −0.26%）。本门量 λ₁_E − λ₁(x,ρ) 在 ρ→1 的**阶**：
      若是 O((ρ−1)²)，则「ε-邻域」宽 ~√ε，盒数 ~1/√ε，有限且小。
numerical-evidence（float64）。
"""
import os; os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util, numpy as np
p = "theory/notes/_scratch_ir/_chk_main_kato_temple_two_atom_prototype_2026-09-23.py"
s_ = importlib.util.spec_from_file_location("kt", p); m_ = importlib.util.module_from_spec(s_); s_.loader.exec_module(m_)

r, N = 10.0, 500
for lam, mus in ((3.2953, 0.13), (6.0, 0.20)):
    wE = m_.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    l1E = m_.lam_k([lam], [wE], mus, Nc, N)
    print("=" * 96)
    print(f"λ={lam}  μ*={mus}  ω̄_E={wE:.6f}  N_c={Nc}  λ₁_E={l1E:.12f}")

    print("  B1 紧化： max_{ρ>=ρ*, x} λ₁   （x 网格 41 点，ρ 网格对数 60 点到 400）")
    xs = np.linspace(0.005, 0.995, 41)
    rs = np.geomspace(1.0, 400.0, 60)
    Lmat = np.array([[m_.lam_k(*m_.members(lam, wE, x, rho), mus, Nc, N) for rho in rs] for x in xs])
    for rstar in (1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0):
        sel = rs >= rstar
        print(f"     ρ*={rstar:<5}  max λ₁ = {Lmat[:, sel].max():.9f}   (= λ₁_E {Lmat[:, sel].max()/l1E:.4f}×)")

    print("  B2 ρ→1 的阶： max_x [λ₁_E − λ₁(x,ρ)]  vs (ρ−1)")
    prev = None
    for eps in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3):
        rho = 1.0 + eps
        gap = max(l1E - m_.lam_k(*m_.members(lam, wE, x, rho), mus, Nc, N) for x in xs)
        q = "" if prev is None else f"   比值={prev[1]/gap:>8.2f}  (ρ−1 比 {prev[0]/eps:.2f})"
        print(f"     ρ−1={eps:<8.0e}  max 亏 = {gap:.3e}{q}")
        prev = (eps, gap)
    print()
