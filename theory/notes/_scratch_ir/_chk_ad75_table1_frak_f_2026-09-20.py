# -*- coding: utf-8 -*-
"""
AD1975 Table I (p.913) 的亏损因子 f := Tc / (C_inf * sqrt(lambda * wbar2^2))
                                    = Tc / (C_inf * sqrt(lambda) * wbar2)

2026-09-20 主会话从零计算。

★ 数据来源纪律（CLAUDE.md：文本层与页图互为对照）：
  本表**由页图逐行转录**（pdftoppm -f 9 -r 200，见 _ad75_p9-09.png），
  **不用 pdftotext**。原因：AD1975.txt 的文本层把 lambda 列整体**上移一行**
  （Pb 行读成 0.805，实为 In 的 lambda；Pb 真值 1.55）。若信文本层，
  14 行 f 全错，且结论符号反转（0.41 vs 0.52）。

列（页图表头逐字）：Material | w_log(K) | wbar1(K) | wbar2(K) | lambda | w_ph(K) | mu*(w_ph) | Tc(K) | eta(eV/A^2)

非循环性核验：
  - alpha^2F 来自隧道电导反演（McMillan-Rowell），是**实验量**，非 DFT、非拟合公式。
  - Tc 列逐个对上教科书实验值：Pb 7.2(7.19) / Sn 3.75(3.72) / In 3.40(3.41)
    / Hg 4.19(4.15) / Ta 4.48(4.47) / Tl 2.36(2.38) ==> **实验直测**。
  ==> lambda, wbar2, Tc 三者均非由任何 Tc 公式生成 ==> f 非循环。
     （对比：theory/data/ambient_hydride_candidates.csv 的 filling_factor_f
       用 Tc_AD（AD 公式值）==> 闭式循环量，已判定零材料信息。）
"""
import numpy as np

C_INF = 0.1827262477  # P4_g2_analytic.md, CERT-B 下端

# (material, w_log, wbar1, wbar2, lambda, w_ph, mu*, Tc, eta) —— 页图逐行转录，29 行
T1 = [
    ("Pb",              56,  60,  65, 1.55, 110, 0.105, 7.2,  2.4),
    ("In",              68,  79,  89, 0.805,179, 0.097, 3.40, 1.3),
    ("Sn",              99, 110, 121, 0.72, 209, 0.092, 3.75, 2.2),
    ("Hg",              29,  38,  49, 1.6,  162, 0.098, 4.19, 1.4),
    ("Tl",              52,  58,  64, 0.795,127, 0.111, 2.36, 1.2),
    ("Ta",             132, 140, 148, 0.69, 228, 0.093, 4.48, 4.9),
    ("a-Ga",            55,  77, 101, 1.62, 291, 0.095, 8.56, 2.1),
    ("b-Ga",            87, 108, 129, 0.97, 285, 0.092, 5.90, 2.0),
    ("Tl0.9Bi0.1",      48,  55,  62, 0.78, 120, 0.099, 2.30, 1.1),
    ("Pb0.4Tl0.6",      48,  56,  62, 1.15, 121, 0.094, 4.60, 1.6),
    ("Pb0.6Tl0.4",      50,  57,  62, 1.38, 119, 0.103, 5.90, 2.0),
    ("Pb0.8Tl0.2",      50,  56,  61, 1.53, 116, 0.101, 6.80, 2.1),
    ("Pb0.6Tl0.2Bi0.2", 48,  53,  58, 1.81, 112, 0.111, 7.26, 2.3),
    ("Pb0.9Bi0.1",      50,  56,  60, 1.66, 108, 0.081, 7.65, 2.2),
    ("Pb0.8Bi0.2",      46,  52,  57, 1.88, 109, 0.093, 7.95, 2.3),
    ("Pb0.7Bi0.3",      47,  52,  57, 2.01, 110, 0.092, 8.45, 2.4),
    ("Pb0.65Bi0.35",    45,  50,  55, 2.13, 110, 0.093, 8.95, 2.4),
    ("In0.9Tl0.1",      63,  75,  86, 0.85, 176, 0.103, 3.28, 1.4),
    ("In0.73Tl0.27",    55,  67,  77, 0.93, 166, 0.110, 3.36, 1.4),
    ("In0.67Tl0.33",    57,  68,  79, 0.90, 167, 0.110, 3.26, 1.4),
    ("In0.57Tl0.43",    53,  64,  74, 0.85, 165, 0.117, 2.60, 1.3),
    ("In0.5Tl0.5",      53,  64,  73, 0.83, 163, 0.110, 2.52, 1.3),
    ("In0.27Tl0.73",    42,  53,  63, 1.09, 151, 0.094, 3.64, 1.4),
    ("In0.17Tl0.83",    45,  55,  63, 0.98, 144, 0.101, 3.19, 1.3),
    ("In0.07Tl0.93",    49,  56,  63, 0.89, 131, 0.107, 2.77, 1.3),
    ("In2Bi",           46,  57,  67, 1.40, 174, 0.096, 5.6,  1.6),
    ("Sb2Tl7",          37,  48,  58, 1.43, 134, 0.102, 5.2,  1.6),
    ("Bi2Tl",           47,  53,  59, 1.63, 120, 0.101, 6.4,  2.1),
    ("a-Pb0.45Bi0.55",  29,  38,  47, 2.59, 128, 0.116, 7.0,  2.1),
]

print(f"n rows = {len(T1)}")
print()
print("Material           w_log  wbar2  lambda   mu*    Tc(K)    f      Tc/w_log")
print("=" * 78)
rows = []
for m, wl, w1, w2, lam, wph, mus, tc, eta in T1:
    f = tc / (C_INF * np.sqrt(lam) * w2)
    rows.append((m, lam, mus, tc, f, tc / wl))
    print(f"{m:17s} {wl:5.0f}  {w2:5.0f}  {lam:5.2f}  {mus:.3f}  {tc:5.2f}   {f:.4f}   {tc/wl:.4f}")

f = np.array([r[4] for r in rows])
lam = np.array([r[1] for r in rows])
print("=" * 78)
print(f"n = {len(f)}")
print(f"mean  = {f.mean():.4f}    std = {f.std(ddof=1):.4f}   (rel. spread {f.std(ddof=1)/f.mean()*100:.1f}%)")
print(f"min   = {f.min():.4f} ({rows[int(f.argmin())][0]})")
print(f"max   = {f.max():.4f} ({rows[int(f.argmax())][0]})")
print(f"median= {np.median(f):.4f}")
print()
print("*** CEILING CHECK: theorem says f <= 1 for all. ***")
print(f"    max f = {f.max():.4f}  ==> {'PASS (no violation)' if f.max() <= 1 else 'VIOLATION!'}")
print()
print("Compare with H3S/D3S family (f ~ 0.49-0.51, mu* = 0.13):")
print(f"    tunneling-family mean {f.mean():.4f} +- {f.std(ddof=1):.4f}")
print(f"    fraction of rows in [0.45, 0.55]: {np.mean((f >= 0.45) & (f <= 0.55))*100:.0f}%")
print()
print("Correlation of f with lambda (is f a constant, or a function of lambda?):")
print(f"    Pearson r(f, lambda)        = {np.corrcoef(f, lam)[0,1]:+.4f}")
print(f"    Pearson r(f, log lambda)    = {np.corrcoef(f, np.log(lam))[0,1]:+.4f}")
print(f"    f for lambda<1.0  (n={np.sum(lam<1.0):2d}): mean {f[lam<1.0].mean():.4f}  std {f[lam<1.0].std(ddof=1):.4f}")
print(f"    f for lambda>=1.0 (n={np.sum(lam>=1.0):2d}): mean {f[lam>=1.0].mean():.4f}  std {f[lam>=1.0].std(ddof=1):.4f}")

# ---------------------------------------------------------------
# 追加：f 是常数，还是 F(lambda, 形状, mu*)？
# 形状信息由 w_log/wbar2 <= 1 携带（对数平均 <= 二阶矩均方根）
# ---------------------------------------------------------------
print()
print("=" * 78)
print("DECOMPOSITION PROBE: is f constant, or F(lambda, shape, mu*)?")
print("=" * 78)
wlg = np.array([r[1] for r in T1], dtype=float)
wb2 = np.array([r[3] for r in T1], dtype=float)
mus = np.array([r[6] for r in T1], dtype=float)
shape = wlg / wb2
print(f"shape = w_log/wbar2 : min {shape.min():.4f}  max {shape.max():.4f}  mean {shape.mean():.4f}")
print(f"   all <= 1 ? {bool(np.all(shape <= 1))}   (Cauchy-Schwarz direction check)")
print(f"   r(f, shape) = {np.corrcoef(f, shape)[0,1]:+.4f}    r(f, mu*) = {np.corrcoef(f, mus)[0,1]:+.4f}")

X = np.column_stack([np.ones(len(f)), np.log(lam), np.log(shape)])
beta, *_ = np.linalg.lstsq(X, np.log(f), rcond=None)
res = np.log(f) - X @ beta
ss_tot = np.sum((np.log(f) - np.log(f).mean()) ** 2)
print()
print(f"log f = {beta[0]:+.4f} {beta[1]:+.4f}*log(lam) {beta[2]:+.4f}*log(shape)")
print(f"   R^2 = {1 - np.sum(res**2)/ss_tot:.4f}   resid scatter {np.std(res, ddof=3)*100:.2f}%")
X1 = np.column_stack([np.ones(len(f)), np.log(lam)])
b1, *_ = np.linalg.lstsq(X1, np.log(f), rcond=None)
r1 = np.log(f) - X1 @ b1
print(f"   lambda only : R^2 = {1 - np.sum(r1**2)/ss_tot:.4f}   resid scatter {np.std(r1, ddof=2)*100:.2f}%")
print()
print("VERDICT on 'f ~ 0.5 universal':")
print(f"   f spans {f.min():.3f}-{f.max():.3f} = {f.max()/f.min():.2f}x. Not a constant.")
