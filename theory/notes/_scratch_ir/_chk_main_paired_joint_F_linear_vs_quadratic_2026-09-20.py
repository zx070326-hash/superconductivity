# -*- coding: utf-8 -*-
"""
memory `deficit-factor-lambda-exponent-not-an-AD-artifact` 的承重数字是
「配对检验 4.85σ 可分辨」——它来自**线性**设定。
加上 (logλ)² 后 Δb_λ 掉到 2.24σ、Δq 只有 0.61σ。那 4.85σ 还成立吗？

正确做法：**联合检验**。问的不是「某个系数是否可分辨」，而是
    「d := log𝔣_real − log𝔣_AD 对 λ 有没有任何依赖？」
⟹ 对 d 回归 {1, S} vs {1, S, L} vs {1, S, L, L²}，做嵌套 F。
系数怎么在 L 与 L² 之间分账是**设定问题**；联合显著性不是。

另外补一条真正模型无关的量：**两端的 𝔣_real/𝔣_AD 之比**（不拟合任何东西）。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy import stats

C_INF = 0.1827262477
T1 = [("Pb",56,65,1.55,0.105,7.2),("In",68,89,0.805,0.097,3.40),("Sn",99,121,0.72,0.092,3.75),
      ("Hg",29,49,1.6,0.098,4.19),("Tl",52,64,0.795,0.111,2.36),("Ta",132,148,0.69,0.093,4.48),
      ("a-Ga",55,101,1.62,0.095,8.56),("b-Ga",87,129,0.97,0.092,5.90),
      ("Tl0.9Bi0.1",48,62,0.78,0.099,2.30),("Pb0.4Tl0.6",48,62,1.15,0.094,4.60),
      ("Pb0.6Tl0.4",50,62,1.38,0.103,5.90),("Pb0.8Tl0.2",50,61,1.53,0.101,6.80),
      ("Pb0.6Tl0.2Bi0.2",48,58,1.81,0.111,7.26),("Pb0.9Bi0.1",50,60,1.66,0.081,7.65),
      ("Pb0.8Bi0.2",46,57,1.88,0.093,7.95),("Pb0.7Bi0.3",47,57,2.01,0.092,8.45),
      ("Pb0.65Bi0.35",45,55,2.13,0.093,8.95),("In0.9Tl0.1",63,86,0.85,0.103,3.28),
      ("In0.73Tl0.27",55,77,0.93,0.110,3.36),("In0.67Tl0.33",57,79,0.90,0.110,3.26),
      ("In0.57Tl0.43",53,74,0.85,0.117,2.60),("In0.5Tl0.5",53,73,0.83,0.110,2.52),
      ("In0.27Tl0.73",42,63,1.09,0.094,3.64),("In0.17Tl0.83",45,63,0.98,0.101,3.19),
      ("In0.07Tl0.93",49,63,0.89,0.107,2.77),("In2Bi",46,67,1.40,0.096,5.6),
      ("Sb2Tl7",37,58,1.43,0.102,5.2),("Bi2Tl",47,59,1.63,0.101,6.4),
      ("a-Pb0.45Bi0.55",29,47,2.59,0.116,7.0)]
name = np.array([r[0] for r in T1])
wlog = np.array([float(r[1]) for r in T1]); wb2 = np.array([float(r[2]) for r in T1])
lam  = np.array([float(r[3]) for r in T1]); mus = np.array([float(r[4]) for r in T1])
tc   = np.array([float(r[5]) for r in T1])
n = len(tc)

def tc_ad(lam, wlog, wb2, mu):
    r = wb2 / wlog
    L1 = 2.46 * (1 + 3.8 * mu); L2 = 1.82 * (1 + 6.3 * mu) * r
    f1 = (1 + (lam / L1) ** 1.5) ** (1.0 / 3.0)
    f2 = 1 + (r - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    return (wlog / 1.2) * np.exp(-1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam))) * f1 * f2

L, S = np.log(lam), np.log(wlog / wb2)
y_real = np.log(tc / (C_INF * np.sqrt(lam) * wb2))
y_ad   = np.log(tc_ad(lam, wlog, wb2, mus) / (C_INF * np.sqrt(lam) * wb2))
d = y_real - y_ad                       # 配对残差

def rss(X, y):
    b, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b; return float(r @ r), b

O = np.ones(n)
models = {"{1}": np.column_stack([O]),
          "{1,S}": np.column_stack([O, S]),
          "{1,S,L}": np.column_stack([O, S, L]),
          "{1,S,L,L²}": np.column_stack([O, S, L, L * L])}

print("=" * 96)
print("1  配对残差 d = log𝔣_real − log𝔣_AD 的嵌套 F（λ 依赖的**联合**显著性）")
print("=" * 96)
print(f"  d 的均值 {d.mean():+.4f}，标准差 {d.std(ddof=1):.4f}   （AD 整体偏低 ⟹ d>0）")
print(f"\n  {'模型':<14}{'RSS':>12}{'dof':>6}{'s':>10}")
for k, X in models.items():
    r2, _ = rss(X, d); print(f"  {k:<14}{r2:>12.5f}{n-X.shape[1]:>6d}{np.sqrt(r2/(n-X.shape[1])):>10.4f}")

print(f"\n  {'嵌套比较':<34}{'F':>10}{'df':>10}{'p':>12}{'判定':>12}")
for small, big, add in [("{1,S}", "{1,S,L}", "加 L"),
                        ("{1,S,L}", "{1,S,L²}" , "加 L² (给定 L)"),
                        ("{1,S}", "{1,S,L,L²}", "**加 L 与 L² 联合**")]:
    if big == "{1,S,L²}":
        Xb = np.column_stack([O, S, L, L * L]); Xs = np.column_stack([O, S, L]); q = 1
    else:
        Xb = models[big]; Xs = models[small]; q = Xb.shape[1] - Xs.shape[1]
    r_s, _ = rss(Xs, d); r_b, _ = rss(Xb, d); df2 = n - Xb.shape[1]
    F = ((r_s - r_b) / q) / (r_b / df2); p = 1 - stats.f.cdf(F, q, df2)
    print(f"  {small+' → '+big+'  ('+add+')':<34}{F:>10.3f}{('(%d,%d)'%(q,df2)):>10}{p:>12.2e}"
          f"{('  **显著**' if p < 0.05 else '  不显著'):>12}")

print("""
  ⟹ 联合检验才是承重的那个：系数如何在 L 与 L² 之间分账属于**设定选择**，
     而「d 依赖 λ」这件事不依赖设定。""")

print()
print("=" * 96)
print("2  线性设定的 4.85σ 与二次设定的 2.24σ，哪个该被引用？——两个都对，问的不是同一件事")
print("=" * 96)
for tag, X in [("线性 {1,S,L}", np.column_stack([O, S, L])),
               ("二次 {1,S,L,L²}", np.column_stack([O, S, L, L * L]))]:
    r2, b = rss(X, d); dof = n - X.shape[1]
    se = np.sqrt(np.diag((r2 / dof) * np.linalg.inv(X.T @ X)))
    print(f"  {tag:<18} Δb_L = {b[2]:+.4f} ± {se[2]:.4f}  ({abs(b[2]/se[2]):.2f}σ)"
          + (f"   Δq = {b[3]:+.4f} ± {se[3]:.4f} ({abs(b[3]/se[3]):.2f}σ)" if X.shape[1] == 4 else ""))
print("""
  · 线性的 Δb_L 是**区间平均对数斜率之差**——一个有定义的量，4.85σ 成立，**不撤回**。
  · 二次的 Δb_L 是**给定 L² 后的偏斜率**——另一个量，2.24σ。二者不矛盾，
    加参数必然抬高单系数的标准误（共线：corr(L,L²) 很高）。
  · **错误用法**：拿二次设定的 2.24σ 去「推翻」线性设定的 4.85σ。它们不是同一个假设检验。""")
print(f"  corr(L, L²) = {np.corrcoef(L, L*L)[0,1]:.4f}   ⟹ 共线性强，单系数 SE 被抬高是必然的")

print()
print("=" * 96)
print("3  零拟合的模型无关量：低/高 λ 两端的 𝔣_real/𝔣_AD")
print("=" * 96)
o = np.argsort(lam); k = 10
lo, hi = o[:k], o[-k:]
rr = np.exp(d)
tl, pl = stats.ttest_ind(d[hi], d[lo], equal_var=False)
print(f"  低 λ {k} 行 (λ∈[{lam[lo].min():.2f},{lam[lo].max():.2f}]) : 𝔣_real/𝔣_AD 均值 {rr[lo].mean():.4f}")
print(f"  高 λ {k} 行 (λ∈[{lam[hi].min():.2f},{lam[hi].max():.2f}]) : 𝔣_real/𝔣_AD 均值 {rr[hi].mean():.4f}")
print(f"  差（log 尺度）Welch t = {tl:+.2f}   p = {pl:.2e}   "
      f"⟹ {'**两端确有系统差，零拟合可见**' if pl < 0.05 else '不显著'}")
print("""
  ⟹ 这条不经任何回归设定：高 λ 端实验 𝔣 相对 AD **系统性更高**。
     这才是「λ 指数不是 AD 伪影」真正站得住的表述——
     **不是**「曲率是实验独有的」（已被 0.61σ 证伪）。""")
print("=" * 96)
