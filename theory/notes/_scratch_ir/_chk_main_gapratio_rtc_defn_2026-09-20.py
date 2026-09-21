# -*- coding: utf-8 -*-
"""R_Tc 的两种定义（有限增量 vs 导数×Δ）差多少 —— 补完门 FAIL 的归因。

主诊断已定案：R_ME 的分歧 **100% 出在 `w_cut`**（500 vs 1069 meV），
而「有限增量 vs 导数×Δ」在 R_ME 上只有 0.1–0.2% ⟹ **我事先的曲率猜测在 R_ME 上被证伪**。
但 `净增益 = R_ME/R_Tc` 两边仍差 5.9%（我 1.4971 vs 执行方 1.4135）
⟹ 剩余分歧必在 **R_Tc**。Tc 对 λ 是指数型 ⟹ 曲率在这里**可能**真的承重。
本脚本只算 Tc（不解能隙），四格全跑。

严格性：numerical-evidence。
"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from theory.solvers.eliashberg_iso import Spectrum, Solver

LAM0, WE, MU0, N = 1.876, 106.9, 0.13, 1600
DMU, DLAM = 0.05, 0.15 * 1.876


def tc(lam, mu, wcut):
    return Solver(spectrum=Spectrum.einstein(lam=lam, w0=WE),
                  mu_star=mu, w_cut=wcut, n_matsubara=N).tc()


for wcut in (500.0, 1069.0):
    t0 = tc(LAM0, MU0, wcut)
    tlp, tlm = tc(LAM0 + DLAM, MU0, wcut), tc(LAM0 - DLAM, MU0, wcut)
    tmp, tmm = tc(LAM0, MU0 + DMU, wcut), tc(LAM0, MU0 - DMU, wcut)
    hL, hM = 0.1 * DLAM, 0.1 * DMU
    dl = (tc(LAM0 + hL, MU0, wcut) - tc(LAM0 - hL, MU0, wcut)) / (2 * hL)
    dm = (tc(LAM0, MU0 + hM, wcut) - tc(LAM0, MU0 - hM, wcut)) / (2 * hM)

    incA = abs(0.5 * (tlp - tlm)) / abs(0.5 * (tmp - tmm))
    incB = abs(dl * DLAM) / abs(dm * DMU)
    print(f"w_cut = {wcut:>6.1f}   Tc0 = {t0:8.3f} K")
    print(f"   (A) 有限增量 : ΔTc_λ = {0.5*(tlp-tlm):+8.3f}  ΔTc_μ = {0.5*(tmp-tmm):+8.3f}"
          f"   R_Tc = {incA:.4f}")
    print(f"   (B) 导数×Δ   : dTc/dλ·Δλ = {dl*DLAM:+8.3f}  dTc/dμ*·Δμ = {dm*DMU:+8.3f}"
          f"   R_Tc = {incB:.4f}")
    print(f"   (A)/(B) = {incA/incB:.4f}   ⟹ 曲率分叉 {100*(incA/incB-1):+.1f}%")
    print(f"   [参考] dTc/dλ = {dl:+.2f} K,  dTc/dμ* = {dm:+.2f} K")
