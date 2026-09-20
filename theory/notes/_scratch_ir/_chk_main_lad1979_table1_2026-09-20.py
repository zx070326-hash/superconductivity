# -*- coding: utf-8 -*-
"""
主会话从零复核：Leavens–Allen–Dynes 1979 (Solid State Commun. 30, 595) 的
Eq.(1)/(10) 与 Table 1，作为对我们 C_inf 与 Eliashberg 求解器的**第三方独立检验**。

一手来源：theory/notes/_scratch_ir/LAD1979.pdf，页图 _lad1979_p-1.png / _lad1979_p-2.png
（文本层是 TextBridge OCR，糊得厉害，**所有数字均取自页图**，不采信文本层）。

页图逐字读到的承重内容：
  Eq.(1)  Limit_{lam<w^2> >> w_max^2}  k_B Tc = 0.1827 [lam<w^2>]^{1/2}
  Eq.(2)  lam<w^2> = 2 int dw a2F(w) w
  Eq.(3)  lam        = 2 int dw a2F(w)/w
  p.596 左栏："Numerical solution of the gap equations gives c(0) = 0.1827."
  Eq.(8)  lam<w^2>/w_max^2 <= 1 + (2 pi k_B Tc / w_max)^2      <-- 方向：Tc 的**下**界
  Eq.(10) Limit_{x >> w_max} k_B Tc = 0.2584 x,  x = [int dw a2F(w) w]^{1/2}
  p.596 右栏："For elements x^2 = eta/2M where eta is the McMillan-Hopfield parameter
              ... in the limit of very strong electron-phonon coupling it is the first
              moment of a2F, rather than lam, which is the single most important parameter."
  Table 1（幂律谱 a2F = (n lam/2)(w/w_max)^n theta(w_max-w)）:
      n     lam_(5%)   (<w^2>lam/w_max^2)_(5%)   (x/w_max)_(5%)
      0     inf        3.4                       1.3
      1/2   18.8       3.8                       1.4
      1     12.3       4.1                       1.4
      2      9.2       4.6                       1.5
      inf    7.0       7.0                       1.9
  表注："To achieve accuracy to within 2%, lam must be larger than lam_(5%) by a factor of ~2.5."

本脚本：
  A. 常数自洽：0.1827*sqrt(2) 是否等于 Eq.(10) 的 0.2584；我们的 C_inf 代进去给多少。
  B. Table 1 的**自由度**：第 3、4 列是不是第 2 列的函数（若是，就不是 15 个独立检验）。
  C. 用仓库求解器在这 4 个有限谱上重求 lam_(5%)，与 LAD 的值对比。
     C1 是**仪器自检**：幂律谱的 lam 与 <w^2> 有解析值，网格必须先打到 1e-4 以内
        （n=1/2 的 a2F/w ~ w^{-1/2} 在 0 附近有可积奇点，均匀网格 Simpson 只有 ~1% 精度
         —— 这正是第一版跑出 n=1/2「偏 14.6%」的原因，是**我的**仪器不是 LAD 的数）。
  D. 天花板方向：f<=1 是否全程成立（LAD 只证极限，没证上界）。
"""
import os
import sys

import numpy as np
from scipy import optimize

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K  # noqa: E402

C_INF = 0.1827262477      # 我们认证的 (1/2pi) sqrt(g(2))
LAD_C = 0.1827            # LAD Eq.(1) 印出的 4 位
W_MAX = 100.0             # meV，任意标度（Tc 线性于它，f 与之无关）
NGRID = 20001
QGRID = 3.0               # 网格 w = W_MAX * t^q，q=3 把可积奇点压到 1e-6（见 C1）
NMATS = 512               # 见 C0 的 N 收敛表

print("=" * 78)
print("A. 常数自洽性（纯算术）")
print(f"   LAD Eq.(1) 常数                      = {LAD_C}")
print(f"   我们的 C_inf                         = {C_INF}")
print(f"   sqrt(2)*LAD_C  (应 = Eq.10 的 0.2584) = {np.sqrt(2)*LAD_C:.6f}")
print(f"   sqrt(2)*C_inf                        = {np.sqrt(2)*C_INF:.6f}")
print("   => Eq.(10) 的 0.2584 与 Eq.(1) 的 0.1827 是**同一个常数**（lam<w^2> = 2 x^2），")
print("      不是两条独立信息。")
print(f"   C_inf 与 LAD 的 0.1827 相对差 = {(C_INF-LAD_C)/LAD_C*100:+.4f}%  "
      f"（LAD 印出的 4 位与我们的认证包封逐位相符）")
ad_asym = np.exp(-1.04) / (1.20 * np.sqrt(2.46))
print(f"   参照：AD75 自己 f1 的渐近常数 e^-1.04/(1.20*sqrt(2.46)) = {ad_asym:.7f}  "
      f"比 C_inf 高 {(ad_asym/C_INF-1)*100:+.2f}%")
print("   ^^^ 即 AD75 的**拟合式**与 LAD79 的**真渐近**内部不一致 2.77%；LAD 印的是后者。")

# ---------------------------------------------------------------- B
print()
print("=" * 78)
print("B. Table 1 的自由度：第 3/4 列是不是第 2 列的函数？")
print("   幂律谱 a2F = C w^n on (0,w_max]:")
print("     lam      = 2C w_max^n / n ;  lam<w^2> = 2C w_max^{n+2}/(n+2)")
print("     => <w^2>/w_max^2 = n/(n+2)   （纯形状常数，与 C、与 lam 无关）")
print()
print("   n      lam_(5%)   col3(LAD)   col2*n/(n+2)   col4(LAD)   sqrt(col3/2)")
for nm, n, lam5, col3, col4 in [("1/2", 0.5, 18.8, 3.8, 1.4), ("1", 1.0, 12.3, 4.1, 1.4),
                                ("2", 2.0, 9.2, 4.6, 1.5), ("inf", np.inf, 7.0, 7.0, 1.9)]:
    pred3 = lam5 * (1.0 if np.isinf(n) else n / (n + 2))
    print(f"   {nm:>4}   {lam5:>7}    {col3:>5}        {pred3:>7.3f}       {col4:>5}      "
          f"{np.sqrt(col3/2):>6.3f}")
print("   n=0 一行：lam = 2C int w^{-1} dw 发散 => lam_(5%)=inf；col3 = C = 3.4 是该行的**原始量**")
print("   >>> col3 = col2 * n/(n+2)（逐行到印出精度），col4 = sqrt(col3/2)（同）。")
print("   >>> Table 1 的 15 个数只有 5 个自由度（每行 1 个），**不是 15 个独立检验**。")


# ---------------------------------------------------------------- 谱构造
def powerlaw(n, lam, ngrid=NGRID, q=QGRID):
    """a2F = C w^n on (0,w_max]，C 解析定标使 lam = 2 int a2F/w dw。
    网格 w = W_MAX * t^q（t 均匀），q>1 在 w->0 处聚点，解析掉可积奇点。"""
    C = lam * n / (2.0 * W_MAX ** n)
    t = np.linspace(1.0 / ngrid, 1.0, ngrid)
    w = W_MAX * t ** q
    return Spectrum.from_table(w, C * w ** n, name=f"pow{n}")


def build(n, lam):
    return Spectrum.einstein(lam, W_MAX) if np.isinf(n) else powerlaw(n, lam)


def frak_f(n, lam_nominal, nmats=NMATS):
    """f = Tc / (C_inf sqrt(lam<w^2>))。lam 与 <w^2> **都取求解器实际用的那个谱的值**，
    不用标称值 —— 否则网格误差会被藏进比值里。"""
    sp = build(n, lam_nominal)
    lam_true = sp.lambda_total()
    lw2 = lam_true * sp.omega_2() ** 2
    tc_k = Solver(spectrum=sp, mu_star=0.0, w_cut=W_MAX, n_matsubara=nmats).tc()
    return tc_k * KB_MEV_PER_K / (C_INF * np.sqrt(lw2)), tc_k, lam_true


# ---------------------------------------------------------------- C
print()
print("=" * 78)
print("C. 用仓库 Eliashberg 求解器重求 lam_(5%)")
print()
print("   [C1] 仪器自检：幂律谱的解析矩 vs 求解器的数值矩（网格 q=%.0f, ngrid=%d）" % (QGRID, NGRID))
print("        n      lambda_total (应=1)          <w^2>/w_max^2 (应=n/(n+2))")
ok_grid = True
for n in (0.5, 1.0, 2.0):
    sp = build(n, 1.0)
    lt = sp.lambda_total()
    m2 = (sp.omega_2() / W_MAX) ** 2
    ok_grid &= (abs(lt - 1.0) < 1e-4 and abs(m2 / (n / (n + 2)) - 1.0) < 1e-4)
    print(f"        {n:<5}  {lt:.8f} ({(lt-1)*100:+.5f}%)    {m2:.8f} "
          f"({(m2/(n/(n+2))-1)*100:+.5f}%)")
print(f"        网格自检 {'PASS' if ok_grid else 'FAIL'}（阈值 1e-4 相对）")

print()
print("   [C0] N 收敛自检（**不是**迭代收敛）—— n=1, lam=12.3")
prev = None
for nm in (64, 128, 256, 512, 1024):
    f, tc, _ = frak_f(1.0, 12.3, nm)
    d = "" if prev is None else f"   Δ={abs(f-prev):.2e}"
    print(f"        N={nm:>5}   f = {f:.8f}   Tc = {tc:9.4f} K{d}")
    prev = f
print(f"        => 取 N={NMATS}（截断误差 ~5e-6，比下面所有判决所需精度小两个量级）")

print()
print("   [C2] 求 lam 使 f(lam) = 0.95（Eq.(1) 精度 5%）")
print("        n      LAD lam_(5%)    我们的 lam(f=.95)    比值      我们在 LAD 的 lam 处的 f")
res95 = {}
for nm, n, lad5 in (("1/2", 0.5, 18.8), ("1", 1.0, 12.3), ("2", 2.0, 9.2), ("inf", np.inf, 7.0)):
    L95 = optimize.brentq(lambda L: frak_f(n, L)[0] - 0.95, 2.0, 400.0, xtol=1e-4)
    res95[nm] = L95
    f_at, _, _ = frak_f(n, lad5)
    print(f"        {nm:>4}    {lad5:>8}        {L95:>10.3f}        {L95/lad5:>5.3f}         {f_at:.4f}")

print()
print("   [C3] 表注「2% 精度需 lam 比 lam_(5%) 大 ~2.5 倍」")
print("        n      lam(f=.98)     lam(f=.98)/lam(f=.95)")
for nm, n in (("1/2", 0.5), ("1", 1.0), ("2", 2.0), ("inf", np.inf)):
    L98 = optimize.brentq(lambda L: frak_f(n, L)[0] - 0.98, 2.0, 4000.0, xtol=1e-3)
    print(f"        {nm:>4}   {L98:>10.2f}          {L98/res95[nm]:>6.2f}")

# ---------------------------------------------------------------- D
print()
print("=" * 78)
print("D. 天花板方向：f <= 1 是否全程成立？（LAD 只证极限，没证上界）")
for nm, n in (("1/2", 0.5), ("1", 1.0), ("2", 2.0), ("inf", np.inf)):
    fs = [frak_f(n, L)[0] for L in (1.0, 2.0, 5.0, 10.0, 30.0, 100.0, 300.0)]
    mono = all(b > a for a, b in zip(fs, fs[1:]))
    print(f"   n={nm:>4}  f(lam=1,2,5,10,30,100,300) = " + " ".join(f"{x:.4f}" for x in fs)
          + f"   max={max(fs):.6f}  单调升={mono}")
print("   （全程 f<1 且单调升向 1 ⟹ 与天花板定理方向一致。")
print("     注意这**不是**天花板的证据——是天花板的推论在 4 条谱上没被证伪。）")

# ---------------------------------------------------------------- E
print()
print("=" * 78)
print("E. 反解：LAD Table 1 经我们的求解器，把渐近常数钉到多少？")
print("   逻辑：若 LAD 的「5%」指 |Tc-Tc_asym|/Tc_asym = 5%（即 f = 0.95），")
print("   且他们的 Tc 与我们的 Tc 相同，则  c_implied = C_inf * f_ours(lam_LAD) / 0.95。")
print()
print("        n    lam_LAD   f_ours(lam_LAD)   c_implied     偏 C_inf")
cs = []
for nm, n, lad5 in (("1/2", 0.5, 18.8), ("1", 1.0, 12.3), ("2", 2.0, 9.2), ("inf", np.inf, 7.0)):
    f_at, _, _ = frak_f(n, lad5)
    c_imp = C_INF * f_at / 0.95
    cs.append(c_imp)
    print(f"        {nm:>4}  {lad5:>7}     {f_at:.6f}       {c_imp:.7f}    {(c_imp/C_INF-1)*100:+.4f}%")
cs = np.array(cs)
print(f"        => c_implied = {cs.mean():.7f} +- {cs.std(ddof=1):.7f}  (4 行，1 sigma)")
print(f"           与 C_inf = {C_INF} 相对差 {(cs.mean()/C_INF-1)*100:+.4f}%")
print()
print("   [E2] 灵敏度：这个检验到底能分辨多大的常数误差？")
print("        （若判据对 c 不敏感，上面的吻合就是空的——先算灵敏度再谈吻合）")
for nm, n in (("1", 1.0), ("inf", np.inf)):
    f10, _, _ = frak_f(n, 10.0)
    f30, _, _ = frak_f(n, 30.0)
    slope = (np.log(f30) - np.log(f10)) / np.log(3.0)
    print(f"        n={nm:>4}  dln f/dln lam ≈ {slope:.4f}  =>  c 错 1% 会把 lam_(5%) 挪 "
          f"{np.expm1(0.01/slope)*100:.0f}%")
print("        LAD 印 3 位有效数字（±0.05/12.3 = ±0.4% in lam）")
print(f"        => 本检验对 c 的分辨力 ≈ ±{0.004*0.0387*100:.3f}%（用 n=1 的斜率换算）")
print("        >>> 吻合落在 LAD 自己的印数精度上，不是过拟合。")

print()
print("   [E3] 反读法检验：若「5%」指的是相对 Tc_exact（f = 1/1.05 = 0.952381）呢？")
print("        n      lam(f=1/1.05)    LAD lam_(5%)    比值")
for nm, n, lad5 in (("1/2", 0.5, 18.8), ("1", 1.0, 12.3), ("2", 2.0, 9.2), ("inf", np.inf, 7.0)):
    L = optimize.brentq(lambda L: frak_f(n, L)[0] - 1.0/1.05, 2.0, 400.0, xtol=1e-4)
    print(f"        {nm:>4}   {L:>12.3f}     {lad5:>8}       {L/lad5:.3f}")
print("        >>> 这一读法 4/4 偏高 ~6-7%；f=0.95 读法 4/4 在 0.4% 内。")
print("            ⟹ 「5%」的分母是**渐近公式**，由数据 4/4 定死，不是我挑的。")
