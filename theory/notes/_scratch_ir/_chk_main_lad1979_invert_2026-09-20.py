# -*- coding: utf-8 -*-
"""LAD1979 Table 1 的反解：它把渐近常数钉到多少？（_chk_main_lad1979_table1_2026-09-20.py 的 E 节）"""
import os, sys
import numpy as np
from scipy import optimize
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

C_INF = 0.1827262477; W_MAX = 100.0; NGRID = 20001; QGRID = 3.0; NMATS = 512

def build(n, lam):
    if np.isinf(n):
        return Spectrum.einstein(lam, W_MAX)
    C = lam * n / (2.0 * W_MAX ** n)
    t = np.linspace(1.0 / NGRID, 1.0, NGRID)
    w = W_MAX * t ** QGRID
    return Spectrum.from_table(w, C * w ** n, name=f"pow{n}")

def frak_f(n, lam_nominal, nmats=NMATS):
    sp = build(n, lam_nominal)
    lam_true = sp.lambda_total()
    lw2 = lam_true * sp.omega_2() ** 2
    tc_k = Solver(spectrum=sp, mu_star=0.0, w_cut=W_MAX, n_matsubara=nmats).tc()
    return tc_k * KB_MEV_PER_K / (C_INF * np.sqrt(lw2))

ROWS = (("1/2", 0.5, 18.8), ("1", 1.0, 12.3), ("2", 2.0, 9.2), ("inf", np.inf, 7.0))

print("E. 反解：LAD Table 1 经我们的求解器，把渐近常数钉到多少？")
print("   c_implied = C_inf * f_ours(lam_LAD) / 0.95")
print("        n    lam_LAD   f_ours(lam_LAD)   c_implied     偏 C_inf")
cs = []
for nm, n, lad5 in ROWS:
    f_at = frak_f(n, lad5)
    c = C_INF * f_at / 0.95
    cs.append(c)
    print(f"        {nm:>4}  {lad5:>7}     {f_at:.6f}       {c:.7f}    {(c/C_INF-1)*100:+.4f}%")
cs = np.array(cs)
print(f"        => c_implied = {cs.mean():.7f} +- {cs.std(ddof=1):.7f} (4 行 1sigma)")
print(f"           与 C_inf 相对差 {(cs.mean()/C_INF-1)*100:+.4f}%")

print()
print("   [E2] 灵敏度（先算灵敏度再谈吻合）")
for nm, n in (("1", 1.0), ("inf", np.inf)):
    f10, f30 = frak_f(n, 10.0), frak_f(n, 30.0)
    s = (np.log(f30) - np.log(f10)) / np.log(3.0)
    print(f"        n={nm:>4}  dln f/dln lam = {s:.4f}  =>  c 错 1% 会把 lam_(5%) 挪 "
          f"{np.expm1(0.01/s)*100:.0f}%")
    if nm == "1":
        s1 = s
print(f"        LAD 印 3 位（±0.05/12.3 = ±0.41% in lam）")
print(f"        => 本检验对 c 的分辨力 = {0.0041*s1*100:.4f}%  （即 ±{0.0041*s1*C_INF:.2e} 绝对）")

print()
print("   [E3] 反读法：若「5%」的分母是 Tc_exact（f = 1/1.05 = 0.952381）？")
print("        n      lam(f=1/1.05)   LAD lam_(5%)   比值")
for nm, n, lad5 in ROWS:
    L = optimize.brentq(lambda L: frak_f(n, L) - 1.0/1.05, 2.0, 400.0, xtol=1e-3)
    print(f"        {nm:>4}   {L:>12.3f}    {lad5:>8}      {L/lad5:.3f}")
print("        >>> 这一读法 4/4 偏高；f=0.95 读法 4/4 在 0.4% 内 ⟹ 分母是**渐近公式**，由数据定死。")
