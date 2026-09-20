# -*- coding: utf-8 -*-
"""验收门 FAIL 的归因：R_ME 2.008 (主会话) vs 2.2811 (gapratio-budget)。

★ 这个门的目标数是**我自己**定的 ⟹ FAIL 有两种可能，**方向相反**：
   (a) 执行方实现有别；(b) **我的目标数本身就错**。
   不许默认 (a)。本脚本把两边各自**完整重跑**，逐项隔离差异来源。

两处已从代码读出的实现分歧（不是猜的）：
  D1. 能隙求解温度：主会话 `0.03*Tc`，执行方 `0.02*Tc`。
  D2. 增量的定义：
      主会话 = **有限增量** `½[O(λ+Δλ) − O(λ−Δλ)]`（中心差分，跨度就是实际扰动 Δλ）；
      执行方 = **局部导数 × Δλ**（`dO/dλ` 用小步长算出后乘以 0.15λ）。
      二者只在 O(Δ³) 上等价 ⟹ **曲率大的通道会显著分叉**，而 μ* 通道恰是指数敏感的。

★ 判据不是「谁更接近 2.008」，而是**消费点要哪个量**：
  R_ME 的消费点是「λ 有 15% 不确定、μ* 有 0.05 不确定时，观测量各动多少」——
  这两个不确定度是**有限的**，不是无穷小 ⟹ **有限增量是对的对象**，
  「局部导数 × 有限步」是一次线性化，误差就是曲率项。

严格性：numerical-evidence（float64 Eliashberg 自洽解）。
"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

LAM0, WE, MU0 = 1.876, 106.9, 0.13
N = 1600
DMU, DLAM_FRAC = 0.05, 0.15
DLAM = DLAM_FRAC * LAM0


def obs(lam, mu, tfrac, WCUT):
    spec = Spectrum.einstein(lam=lam, w0=WE)
    solver = Solver(spectrum=spec, mu_star=mu, w_cut=WCUT, n_matsubara=N)
    tc = solver.tc()
    wn, delta, z, conv = solver.solve_gap(tfrac * tc, mixing=0.25, max_iter=4000)
    assert conv, f"未收敛 lam={lam} mu={mu} tfrac={tfrac}"
    return tc, 2.0 * delta[0] / (KB_MEV_PER_K * tc)


# ★ D3（执行方未列出、但代码里就有）：库仑截断
#   主会话脚本 WCUT = 500.0 meV 写死；执行方用 10·ω_E = 1069 meV。
#   memory `h3s-2p3pct-gap-is-omegac-convention`：μ*>0 时 Tc 对 ω_c **非单调**
#   ⟹ 这不是小修正，且 μ* 通道首当其冲。
for WCUT, tfrac in ((500.0, 0.03), (500.0, 0.02), (1069.0, 0.03), (1069.0, 0.02)):
    print("=" * 84)
    print(f"w_cut = {WCUT} meV   T_low = {tfrac}·Tc")
    print("=" * 84)
    tc0, o0 = obs(LAM0, MU0, tfrac, WCUT)
    tlp, olp = obs(LAM0 + DLAM, MU0, tfrac, WCUT)
    tlm, olm = obs(LAM0 - DLAM, MU0, tfrac, WCUT)
    tmp, omp = obs(LAM0, MU0 + DMU, tfrac, WCUT)
    tmm, omm = obs(LAM0, MU0 - DMU, tfrac, WCUT)

    # (A) 有限增量（主会话定义）
    inc_lam = 0.5 * (olp - olm)
    inc_mu = 0.5 * (omp - omm)
    R_inc = abs(inc_lam) / abs(inc_mu)
    dTc_lam = 0.5 * (tlp - tlm)
    dTc_mu = 0.5 * (tmp - tmm)
    R_tc = abs(dTc_lam) / abs(dTc_mu)

    # (B) 局部导数 × Δ（执行方定义）——导数用 1/10 步长的中心差分
    hL, hM = 0.1 * DLAM, 0.1 * DMU
    _, oLp = obs(LAM0 + hL, MU0, tfrac, WCUT); _, oLm = obs(LAM0 - hL, MU0, tfrac, WCUT)
    _, oMp = obs(LAM0, MU0 + hM, tfrac, WCUT); _, oMm = obs(LAM0, MU0 - hM, tfrac, WCUT)
    dO_dlam = (oLp - oLm) / (2 * hL)
    dO_dmu = (oMp - oMm) / (2 * hM)
    R_der = abs(dO_dlam * DLAM) / abs(dO_dmu * DMU)

    print(f"  Tc0 = {tc0:.3f} K    O0 = {o0:.4f}")
    print(f"  (A) 有限增量 : ΔO_λ = {inc_lam:+.5f}   ΔO_μ = {inc_mu:+.5f}   "
          f"R_ME = {R_inc:.4f}   R_Tc = {R_tc:.4f}   净增益 = {R_inc/R_tc:.4f}")
    print(f"  (B) 导数×Δ   : dO/dλ = {dO_dlam:+.4f}  dO/dμ* = {dO_dmu:+.4f}   "
          f"R_ME = {R_der:.4f}")
    print(f"  (A)/(B) = {R_inc/R_der:.4f}   ⟹ 曲率造成的分叉 {100*(R_inc/R_der-1):+.1f}%")
print("=" * 84)
