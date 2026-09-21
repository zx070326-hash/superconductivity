# -*- coding: utf-8 -*-
"""
主会话从零复核：LAD1979 Eq.(7)-(8) 给出的是 Tc 的**下**界，换成 f 就是**地板**。

一手（页图 _lad1979_p-2.png 逐字）：
  Eq.(7)  lam<w^2>/(w_max^2+(2 pi kB Tc)^2)  <=  2 int_0^{wmax} dw a2F(w) w/(w^2+(2 pi kB Tc)^2)  <=  1
          "The left inequality is obvious; the right one is proved in ref. 2."   [ref.2 = AD1975]
  Eq.(8)  lam <w^2>/w_max^2  <=  1 + (2 pi kB Tc / w_max)^2
  p.595 左栏："For simplicity, we consider explicitly only the case in which the
              Coulomb pseudopotential, mu*, is zero."      <== 全节只在 mu*=0

我的推导（纯代数，方向逐步审）：
  (7) 的两端相连  =>  lam<w^2> <= w_max^2 + (2 pi kB Tc)^2
                  =>  (2 pi kB Tc)^2 >= lam<w^2> - w_max^2
                  =>  kB Tc >= (1/2pi) sqrt( lam<w^2> - w_max^2 )        [右端为负时空洞]
  这是 Tc 的**下**界。与我们的天花板 kB Tc < C_inf sqrt(lam<w^2>) 方向相反 => 两侧夹逼。
  除以 C_inf sqrt(lam<w^2>)：
      f := kB Tc/(C_inf sqrt(lam<w^2>))  >=  (1/sqrt(g2)) * sqrt(1 - w_max^2/(lam<w^2>))
  地板常数 1/sqrt(g2)。

本脚本查三件事（都可能推翻上面）：
  F1 纯算术：1/sqrt(g2) 的值；地板 < 1 是否恒成立（否则与天花板冲突 => 我推错了）。
  F2 在 LAD 的四条谱上直接测 f 与地板，看有没有被违反（违反 = 我的代数或 AD 的 (7) 右端有问题）。
  F3 **前提直测**：AD 的 (7) 右端 lam(n=1) <= 1 在我们求解器的 Tc 处成立吗？
     这是整条链唯一非平凡的输入，必须单独测，不能靠 f 的通过率掩盖
     （CLAUDE.md：数值一致性检验必须逐箭头对齐）。
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K  # noqa: E402

C_INF = 0.1827262477
G2_LO, G2_HI = 1.318140211762671, 1.318140212373478   # CERT-D 包围
W_MAX = 100.0
NGRID, QGRID, NMATS = 20001, 3.0, 512


def powerlaw(n, lam):
    C = lam * n / (2.0 * W_MAX ** n)
    t = np.linspace(1.0 / NGRID, 1.0, NGRID)
    w = W_MAX * t ** QGRID
    return Spectrum.from_table(w, C * w ** n, name=f"pow{n}")


def build(n, lam):
    return Spectrum.einstein(lam, W_MAX) if np.isinf(n) else powerlaw(n, lam)


print("=" * 78)
print("F1. 地板常数（纯算术）")
print(f"   1/sqrt(g2)  in  [{1/np.sqrt(G2_HI):.10f}, {1/np.sqrt(G2_LO):.10f}]")
print(f"   => 地板常数 = 0.8709...  < 1 = 天花板  => 与天花板不冲突（自洽性起码过）")
print(f"   检查 C_inf*2pi = sqrt(g2)? {C_INF*2*np.pi:.10f} vs {np.sqrt(G2_LO):.10f} "
      f"（差 {C_INF*2*np.pi-np.sqrt(G2_LO):+.2e}）")

print()
print("=" * 78)
print("F2. 四条谱上直测：f 是否恒 >= 地板？（地板 = 0.87096*sqrt(1 - w_max^2/(lam<w^2>))）")
print("    注意：LAD 的 w_max 对幂律谱是谱的上端点；Einstein 谱 w_max = w_E。")
print()
print(f"    {'n':>5} {'lam':>7} {'lam<w2>/wmax^2':>15} {'f (实测)':>12} {'地板':>10} {'余量':>9}  判决")
floor_c = 1.0 / np.sqrt(G2_HI)     # 取下端 => 地板本身是下估 => 只支持「未违反」
minmargin = np.inf
for nm, n in (("1/2", 0.5), ("1", 1.0), ("2", 2.0), ("inf", np.inf)):
    for lam in (1.0, 2.0, 5.0, 10.0, 30.0, 100.0):
        sp = build(n, lam)
        lt = sp.lambda_total()
        lw2 = lt * sp.omega_2() ** 2
        tc = Solver(spectrum=sp, mu_star=0.0, w_cut=W_MAX, n_matsubara=NMATS).tc()
        f = tc * KB_MEV_PER_K / (C_INF * np.sqrt(lw2))
        ratio = lw2 / W_MAX ** 2
        if ratio <= 1.0:
            print(f"    {nm:>5} {lam:>7.1f} {ratio:>15.4f} {f:>12.6f} {'空洞':>10} {'—':>9}  (8) 不给信息")
            continue
        fl = floor_c * np.sqrt(1.0 - 1.0 / ratio)
        mg = f - fl
        minmargin = min(minmargin, mg)
        print(f"    {nm:>5} {lam:>7.1f} {ratio:>15.4f} {f:>12.6f} {fl:>10.6f} {mg:>+9.6f}  "
              f"{'OK' if mg >= 0 else '**违反**'}")
print(f"    最小余量 = {minmargin:+.6f}  => {'未被违反' if minmargin >= 0 else '被违反（我的代数或前提有问题）'}")

print()
print("=" * 78)
print("F3. 前提直测（唯一非平凡输入）：AD 的 (7) 右端  lam(n=1) <= 1  在 Tc 处成立？")
print("    lam(n) = 2 int dw a2F(w) w/(w^2+(2 pi n kB Tc)^2)  —— 用求解器自己的 lambda_of_nu")
print()
print(f"    {'n':>5} {'lam':>7} {'Tc/K':>10} {'lam(nu=1)':>12}  判决")
worst = -np.inf
for nm, n in (("1/2", 0.5), ("1", 1.0), ("2", 2.0), ("inf", np.inf)):
    for lam in (1.0, 5.0, 30.0, 100.0):
        sp = build(n, lam)
        tc = Solver(spectrum=sp, mu_star=0.0, w_cut=W_MAX, n_matsubara=NMATS).tc()
        nu1 = 2.0 * np.pi * (tc * KB_MEV_PER_K)      # nu_1 = 2 pi T (n=1 的 Matsubara 玻色频)
        l1 = sp.lambda_of_nu(nu1)
        worst = max(worst, l1)
        print(f"    {nm:>5} {lam:>7.1f} {tc:>10.3f} {l1:>12.6f}  {'OK' if l1 <= 1 else '**>1**'}")
print(f"    max lam(nu=1) over all = {worst:.6f}  => AD 的 (7) 右端 "
      f"{'在我们的求解器上成立' if worst <= 1 else '**不成立**'}")
print()
print("    ⚠ 方向说明：F2/F3 全过只说明「没被证伪」。(7) 右端的**证明**在 AD1975，")
print("      本脚本不重证它；地板的严格性以 AD1975 那一步为准（我们未独立复核该证明）。")
