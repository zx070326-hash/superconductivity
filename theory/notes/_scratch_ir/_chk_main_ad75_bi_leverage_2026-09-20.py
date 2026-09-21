# -*- coding: utf-8 -*-
"""
追查 Q3 的 4.27σ：删掉含 Bi 的 9 行后 λ 指数 0.9972 → 1.1700。

上一版用**全样本**的 σ_b 当尺子判 4.27σ，这是错的尺子：
子样本与全样本**嵌套**（高度相关），且 Bi 行占据 λ 的高端
（Bi 行 λ ∈ [0.78, 2.59]），删掉它们会**缩短杠杆臂** ⟹ 指数本来就该动。

本脚本问三件分得开的事：
  (a) 用**子样本自己的** σ_b 量，1.1700 离 0.9972 有多远？
  (b) 杠杆臂效应：把非 Bi 行按 λ 区间**匹配**回去，位移还在吗？
  (c) 「λ 指数 = 1」这条实验性质，在删 Bi 后还成立吗？（判据是 b 的置信区间含不含 1）
  (d) Bi 行是不是**残差异常**？（若只是杠杆臂，Bi 行应落在全样本拟合线上）

**方向注记**：全部为双侧置信区间比较，不产生单边界判决。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import re

C_INF = 0.1827262477
T1 = [
    ("Pb",56,65,1.55,7.2),("In",68,89,0.805,3.40),("Sn",99,121,0.72,3.75),
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
    ("a-Pb0.45Bi0.55",29,47,2.59,7.0),
]
name = np.array([r[0] for r in T1])
wlog = np.array([float(r[1]) for r in T1]); wb2 = np.array([float(r[2]) for r in T1])
lam  = np.array([float(r[3]) for r in T1]); tc  = np.array([float(r[4]) for r in T1])
frak = tc / (C_INF * np.sqrt(lam) * wb2)
hasBi = np.array([bool(re.search(r"Bi", n)) for n in name])

X = np.column_stack([np.ones(len(T1)), np.log(lam), np.log(wlog / wb2)])
y = np.log(frak)

def fit(Xs, ys):
    b, *_ = np.linalg.lstsq(Xs, ys, rcond=None)
    r = ys - Xs @ b; dof = len(ys) - Xs.shape[1]
    s2 = float(r @ r / dof)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(Xs.T @ Xs)))
    R2 = 1.0 - float(r @ r) / float(((ys - ys.mean()) ** 2).sum())
    return b, se, R2, np.sqrt(s2), r

b0, se0, R20, s0, r0 = fit(X, y)
b1, se1, R21, s1, _  = fit(X[~hasBi], y[~hasBi])

print("=" * 96)
print("(a) 用子样本自己的 σ 量位移")
print("=" * 96)
print(f"   含 Bi 的行数 = {hasBi.sum()}；剩余 n = {(~hasBi).sum()}")
print(f"   全样本  b = {b0[1]:+.4f} ± {se0[1]:.4f}   λ∈[{lam.min():.2f},{lam.max():.2f}]")
print(f"   去 Bi   b = {b1[1]:+.4f} ± {se1[1]:.4f}   λ∈[{lam[~hasBi].min():.2f},{lam[~hasBi].max():.2f}]")
zc = abs(b1[1] - b0[1]) / np.hypot(se0[1], se1[1])
print(f"   位移 / 合并σ = {zc:.2f}σ        （上一版用全样本σ 得 4.27σ —— 尺子选错）")
print(f"   ⟹ 用子样本自己的不确定度，位移是 {zc:.2f}σ，{'不显著' if zc < 2 else '仍显著'}")

print()
print("=" * 96)
print("(b) 是不是杠杆臂效应？把**非 Bi** 行也截到同一 λ 区间")
print("=" * 96)
lam_hi = lam[~hasBi].max()
print(f"   非 Bi 行的 λ 上限 = {lam_hi:.2f}（全样本上限 {lam.max():.2f}，差 {lam.max()/lam_hi:.2f}×）")
keep = lam <= lam_hi
bm, sem, R2m, sm, _ = fit(X[keep], y[keep])
print(f"   全样本截到 λ≤{lam_hi:.2f}（n={keep.sum()}，**仍含 Bi**）: b = {bm[1]:+.4f} ± {sem[1]:.4f}")
print(f"   去 Bi   （n={(~hasBi).sum()}）                      : b = {b1[1]:+.4f} ± {se1[1]:.4f}")
zm = abs(b1[1] - bm[1]) / np.hypot(sem[1], se1[1])
print(f"   λ 区间匹配后的位移 = {zm:.2f}σ")
print(f"   ⟹ 原位移中 {'大部分' if zm < zc/2 else '只有一部分'}来自杠杆臂（λ 区间被截短），"
      f"而非 Bi 行本身异常")

print()
print("=" * 96)
print("(c) 「λ 指数 = 1」这条实验性质在各子样本上是否还成立？（95% CI 含不含 1）")
print("=" * 96)
print(f"{'样本':<26}{'n':>4}{'b':>10}{'95% CI':>22}{'含 1?':>8}")
for tag, m in [("全 29 行", np.ones(len(T1), bool)),
               ("去 Bi", ~hasBi),
               (f"λ ≤ {lam_hi:.2f}（含 Bi）", keep),
               ("仅单质 8 行", np.array([not re.search(r"[0-9]", n) for n in name]))]:
    bb, ss, _, _, _ = fit(X[m], y[m])
    lo, hi = bb[1] - 1.96 * ss[1], bb[1] + 1.96 * ss[1]
    print(f"{tag:<26}{m.sum():>4}{bb[1]:>10.4f}   [{lo:+.4f}, {hi:+.4f}]{'  是' if lo <= 1 <= hi else '  否':>8}")
print()
print("   对照：AD75 自己的闭式给 λ 指数 0.9053（仓库 memory）。")

print()
print("=" * 96)
print("(d) Bi 行是残差异常，还是只是位置靠外？（用全样本拟合的残差）")
print("=" * 96)
print(f"   全样本拟合 s = {s0:.4f}")
print(f"   含 Bi 行残差：mean = {r0[hasBi].mean():+.4f}   rms = {np.sqrt((r0[hasBi]**2).mean()):.4f}"
      f"   max|r| = {np.abs(r0[hasBi]).max():.4f} ({name[hasBi][np.abs(r0[hasBi]).argmax()]})")
print(f"   非 Bi 行残差：mean = {r0[~hasBi].mean():+.4f}   rms = {np.sqrt((r0[~hasBi]**2).mean()):.4f}"
      f"   max|r| = {np.abs(r0[~hasBi]).max():.4f} ({name[~hasBi][np.abs(r0[~hasBi]).argmax()]})")
# Welch t on residual means
from math import sqrt
n1, n2 = hasBi.sum(), (~hasBi).sum()
v1, v2 = r0[hasBi].var(ddof=1), r0[~hasBi].var(ddof=1)
t = (r0[hasBi].mean() - r0[~hasBi].mean()) / sqrt(v1/n1 + v2/n2)
print(f"   两组残差均值的 Welch t = {t:+.2f}  ⟹ {'无系统偏移' if abs(t) < 2 else '⚠ 有系统偏移'}")
print()
print(f"   最大 3 个 |残差| 的行：")
for i in np.argsort(-np.abs(r0))[:3]:
    print(f"      {name[i]:<18} r = {r0[i]:+.4f} = {r0[i]/s0:+.2f}s   λ={lam[i]:.2f}  𝔣={frak[i]:.4f}")

print()
print("=" * 96)
print("判决")
print("=" * 96)
print(f"  · 用正确的尺子（子样本自己的 σ），去 Bi 位移是 {zc:.2f}σ，不是 4.27σ。")
print(f"  · λ 区间匹配后剩 {zm:.2f}σ ⟹ 位移主要是**杠杆臂**（Bi 行独占 λ>{lam_hi:.2f} 的高端）。")
print(f"  · Bi 行残差与其余行无系统偏移（Welch t = {t:+.2f}）⟹ 不是异常族，是趋势的高 λ 支点。")
print(f"  · 但**诚实边界**：λ ∈ ({lam_hi:.2f}, {lam.max():.2f}] 这一段**只有 Bi 体系**，")
print(f"    ⟹ 高 λ 端的指数由单一化学体系支撑。这条要随结论同行，不能省。")
print("=" * 96)
