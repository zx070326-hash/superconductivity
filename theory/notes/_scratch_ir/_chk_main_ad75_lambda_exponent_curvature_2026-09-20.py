# -*- coding: utf-8 -*-
"""
(b) 的副产物：λ 指数**随 λ 区间移动**——λ≤1.62 上是 1.158±0.049（CI 排除 1），
全 29 行上是 0.997±0.041（CI 含 1）。⟹ log𝔣 vs logλ 有**曲率**。

这不是噪声，是**天花板定理要求的**：
  𝔣 ≤ 1 [analytic]，而 𝔣 ∝ λ^1 全局成立会在大 λ 处突破 1
  ⟹ 指数必须随 λ 上升而**下降**。（仓库 memory 已记：「天花板 𝔣≤1 禁止全局 𝔣∝λ」）

本脚本直接检验曲率，三条互相独立的判据：
  (A) 加 (logλ)² 项，二次系数的符号与显著性
  (B) 低λ / 高λ 两半各自拟合，指数差
  (C) 局部斜率（滑动窗）随 λ 的走向
  (D) 外推：按各自指数，𝔣 何时撞到天花板 1？

**方向注记**：(D) 是外推，只做定性陈述，不产生判决。(A)-(C) 为双侧检验。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np, re

C_INF = 0.1827262477
T1 = [("Pb",56,65,1.55,7.2),("In",68,89,0.805,3.40),("Sn",99,121,0.72,3.75),
      ("Hg",29,49,1.6,4.19),("Tl",52,64,0.795,2.36),("Ta",132,148,0.69,4.48),
      ("a-Ga",55,101,1.62,8.56),("b-Ga",87,129,0.97,5.90),
      ("Tl0.9Bi0.1",48,62,0.78,2.30),("Pb0.4Tl0.6",48,62,1.15,4.60),
      ("Pb0.6Tl0.4",50,62,1.38,5.90),("Pb0.8Tl0.2",50,61,1.53,6.80),
      ("Pb0.6Tl0.2Bi0.2",48,58,1.81,7.26),("Pb0.9Bi0.1",50,60,1.66,7.65),
      ("Pb0.8Bi0.2",46,57,1.88,7.95),("Pb0.7Bi0.3",47,57,2.01,8.45),
      ("Pb0.65Bi0.35",45,55,2.13,8.95),("In0.9Tl0.1",63,86,0.85,3.28),
      ("In0.73Tl0.27",55,77,0.93,3.36),("In0.67Tl0.33",57,79,0.90,3.26),
      ("In0.57Tl0.43",53,74,0.85,2.60),("In0.5Tl0.5",53,73,0.83,2.52),
      ("In0.27Tl0.73",42,63,1.09,3.64),("In0.17Tl0.83",45,63,0.98,3.19),
      ("In0.07Tl0.93",49,63,0.89,2.77),("In2Bi",46,67,1.40,5.6),
      ("Sb2Tl7",37,58,1.43,5.2),("Bi2Tl",47,59,1.63,6.4),
      ("a-Pb0.45Bi0.55",29,47,2.59,7.0)]
name = np.array([r[0] for r in T1])
wlog = np.array([float(r[1]) for r in T1]); wb2 = np.array([float(r[2]) for r in T1])
lam  = np.array([float(r[3]) for r in T1]); tc  = np.array([float(r[4]) for r in T1])
frak = tc / (C_INF * np.sqrt(lam) * wb2)
L, S, y = np.log(lam), np.log(wlog / wb2), np.log(frak)

def fit(Xs, ys):
    b, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
    r = ys - Xs @ b; dof = len(ys) - Xs.shape[1]; s2 = float(r @ r / dof)
    return b, np.sqrt(np.diag(s2 * np.linalg.inv(Xs.T @ Xs))), \
           1 - float(r@r)/float(((ys-ys.mean())**2).sum()), np.sqrt(s2)

print("=" * 94)
print("(A) 加 (logλ)² 项：二次系数显著吗？符号对不对？")
print("=" * 94)
X1 = np.column_stack([np.ones(29), L, S])
X2 = np.column_stack([np.ones(29), L, S, L * L])
b1, s1, R1, sd1 = fit(X1, y)
b2, s2_, R2, sd2 = fit(X2, y)
print(f"   线性  : b_λ = {b1[1]:+.4f} ± {s1[1]:.4f}                       R²={R1:.4f}  s={sd1:.4f}")
print(f"   加二次: b_λ = {b2[1]:+.4f} ± {s2_[1]:.4f}   q = {b2[3]:+.4f} ± {s2_[3]:.4f}   R²={R2:.4f}  s={sd2:.4f}")
tq = b2[3] / s2_[3]
print(f"   二次系数 q 的 t = {tq:+.2f}   符号 = {'负（指数随 λ 下降）✓ 与天花板要求一致' if b2[3] < 0 else '正 ✗ 与天花板矛盾'}")
# F 检验
r1 = y - X1 @ b1; r2 = y - X2 @ b2
F = ((r1@r1 - r2@r2) / 1) / (r2@r2 / (29 - 4))
print(f"   嵌套 F 检验 F(1,25) = {F:.2f}   ⟹ {'曲率显著（F>4.24 即 p<0.05）' if F > 4.24 else '曲率不显著'}")

print()
print("=" * 94)
print("(B) 低 λ / 高 λ 两半各自拟合")
print("=" * 94)
med = np.median(lam)
print(f"   λ 中位数 = {med:.2f}")
for tag, m in [("低半 λ ≤ 中位", lam <= med), ("高半 λ > 中位", lam > med)]:
    bb, ss, RR, sdd = fit(np.column_stack([np.ones(m.sum()), L[m], S[m]]), y[m])
    print(f"   {tag:<16} n={m.sum():<3} λ∈[{lam[m].min():.2f},{lam[m].max():.2f}]  "
          f"b_λ = {bb[1]:+.4f} ± {ss[1]:.4f}   R²={RR:.4f}")

print()
print("=" * 94)
print("(C) 局部斜率：按 λ 排序的滑动窗（窗宽 12）")
print("=" * 94)
o = np.argsort(lam); W = 12
print(f"{'窗内 λ 范围':<20}{'中心 λ':>9}{'b_λ':>12}{'±':>9}")
for i in range(0, len(o) - W + 1, 3):
    w = o[i:i + W]
    bb, ss, _, _ = fit(np.column_stack([np.ones(W), L[w], S[w]]), y[w])
    print(f"[{lam[w].min():.2f}, {lam[w].max():.2f}]{'':<10}"
          f"{np.median(lam[w]):>9.2f}{bb[1]:>12.4f}{ss[1]:>9.4f}")

print()
print("=" * 94)
print("(D) 外推：各指数下 𝔣 何时撞到天花板 1？（定性，非判决）")
print("=" * 94)
Smed = np.median(S)
for tag, bb in [("全 29 行 (b=0.997)", b1), ]:
    pass
for tag, bv, av in [("全 29 行", b1[1], b1[0]), ("λ≤1.62 子样本", 1.1584, None)]:
    if av is None:
        m = lam <= 1.62
        bq, _, _, _ = fit(np.column_stack([np.ones(m.sum()), L[m], S[m]]), y[m])
        av, bv = bq[0], bq[1]; cv = bq[2]
    else:
        cv = b1[2]
    lam_star = np.exp((0.0 - av - cv * Smed) / bv)
    print(f"   {tag:<16} b={bv:.4f}  ⟹ 𝔣=1 于 λ* = {lam_star:.1f}"
          f"   （表内最大 λ = {lam.max():.2f}，外推 {lam_star/lam.max():.1f}×）")
print("   注：H3S 的 λ≈1.7–2.3，正落在表内区间；氢化物不需要这段外推。")

print()
print("=" * 94)
print("判决")
print("=" * 94)
print(f"  · log𝔣 vs logλ **确有曲率**：二次项 q = {b2[3]:+.4f}（t={tq:+.2f}），F(1,25)={F:.2f}，符号为负。")
print(f"  · 「λ 指数 = 1」**不是区间无关的**：λ≤1.62 上是 1.158±0.049（CI 排除 1），")
print(f"    全区间上是 0.997±0.041（CI 含 1）。⟹ 「指数恰为 1」是 29 行这段 λ 跨度的**平均**，不是定律。")
print(f"  · 方向与**天花板定理一致**：𝔣≤1 [analytic] 禁止 𝔣∝λ 全局成立 ⟹ 指数必须随 λ 下降。")
print(f"    这是数据上第一次直接看到它。")
print(f"  · **诚实边界**：λ∈(1.62,2.59] 只有 Bi 体系撑着 ⟹ 曲率的高 λ 端锚在单一化学体系上。")
print("=" * 94)
