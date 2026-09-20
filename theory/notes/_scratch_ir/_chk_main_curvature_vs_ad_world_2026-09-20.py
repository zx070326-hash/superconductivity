# -*- coding: utf-8 -*-
"""
曲率 q = −0.525 是**实验内容**，还是 AD 公式自带的形状？

判决性检验（与推翻「λ 指数是 AD 伪影」时同一结构）：
  假设 AD 就是真理，用**同样 29 组 (λ, ω_log, ω̄₂, μ*)** 生成 Tc_AD ⟹ 𝔣_AD，
  跑**同一个二次回归**，比 q。两回归**共用设计矩阵** ⟹ 必须**配对**检验。

配对做法：令 d := log𝔣_real − log𝔣_AD，对同一 X 回归 ⟹ L² 的系数**就是** q_real − q_AD，
其标准误由 d 的残差给出（不是两个独立 σ 的平方和，那会高估不确定度）。

AD75 公式（f₂ **无开方**，见 memory ad75-f2-has-no-square-root）：
  Tc = (ω_log/1.2)·exp[−1.04(1+λ)/(λ − μ*(1+0.62λ))]·f₁·f₂
  f₁ = [1 + (λ/Λ₁)^{3/2}]^{1/3},  Λ₁ = 2.46(1+3.8μ*)
  f₂ = 1 + (ω̄₂/ω_log − 1)·λ²/(λ² + Λ₂²),  Λ₂ = 1.82(1+6.3μ*)(ω̄₂/ω_log)

⚠ 出处限定：μ* 一列**转录自兄弟脚本** `_chk_h3s_on_ad75_trend_2026-09-20.py`，
  本轮**未**重读页图。故本脚本同时跑「逐行 μ*」与「统一 μ*=0.10」两版——
  若两版结论一致，则 μ* 的出处对本结论**不承重**。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

C_INF = 0.1827262477
# (name, ω_log, ω̄₂, λ, μ*, Tc_exp)   —— ω 均为 K
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

def tc_ad(lam, wlog, wb2, mu):
    r = wb2 / wlog
    L1 = 2.46 * (1 + 3.8 * mu)
    L2 = 1.82 * (1 + 6.3 * mu) * r
    f1 = (1 + (lam / L1) ** 1.5) ** (1.0 / 3.0)
    f2 = 1 + (r - 1) * lam ** 2 / (lam ** 2 + L2 ** 2)
    expo = -1.04 * (1 + lam) / (lam - mu * (1 + 0.62 * lam))
    return (wlog / 1.2) * np.exp(expo) * f1 * f2

def fit(X, y):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b; dof = len(y) - X.shape[1]; s2 = float(r @ r / dof)
    return b, np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X))), \
           1 - float(r @ r) / float(((y - y.mean()) ** 2).sum()), np.sqrt(s2)

L, S = np.log(lam), np.log(wlog / wb2)
X1 = np.column_stack([np.ones(29), L, S])
X2 = np.column_stack([np.ones(29), L, S, L * L])
y_real = np.log(tc / (C_INF * np.sqrt(lam) * wb2))

print("=" * 98)
print("0  先验收 AD 公式实现：Tc_AD vs Tc_exp（AD75 自陈 3.85% 级别的误差）")
print("=" * 98)
t_ad = tc_ad(lam, wlog, wb2, mus)
err = (t_ad - tc) / tc
print(f"  逐行相对误差: 均值 {err.mean()*100:+.2f}%   RMS {np.sqrt((err**2).mean())*100:.2f}%   "
      f"最劣 {name[np.argmax(abs(err))]} {err[np.argmax(abs(err))]*100:+.2f}%")
print(f"  ⟹ {'✓ 与 AD75 自陈量级相符，实现通过' if np.sqrt((err**2).mean()) < 0.10 else '✗ 实现可疑，停止'}")

print()
print("=" * 98)
print("1  AD 世界的曲率：同样 29 组自变量，同一个二次回归")
print("=" * 98)
hdr = f"  {'样本':<22}{'b_λ':>12}{'±':>8}{'q':>12}{'±':>8}{'t(q)':>8}{'R²':>9}"
print(hdr); print("  " + "-" * (len(hdr) - 2))
res = {}
for tag, mu_use in [("实验（真实 Tc）", None), ("AD 世界（逐行 μ*）", mus),
                    ("AD 世界（统一 μ*=0.10）", np.full(29, 0.10))]:
    yy = y_real if mu_use is None else np.log(tc_ad(lam, wlog, wb2, mu_use) / (C_INF * np.sqrt(lam) * wb2))
    b, se, R2, s = fit(X2, yy)
    res[tag] = (yy, b, se)
    print(f"  {tag:<22}{b[1]:>12.4f}{se[1]:>8.4f}{b[3]:>12.4f}{se[3]:>8.4f}{b[3]/se[3]:>8.2f}{R2:>9.4f}")

print()
print("=" * 98)
print("2  配对检验：d := log𝔣_real − log𝔣_AD 对同一 X 回归，L² 系数即 q_real − q_AD")
print("=" * 98)
for tag in ["AD 世界（逐行 μ*）", "AD 世界（统一 μ*=0.10）"]:
    d = res["实验（真实 Tc）"][0] - res[tag][0]
    bd, sed, _, _ = fit(X2, d)
    q_r = res["实验（真实 Tc）"][1][3]; q_a = res[tag][1][3]
    z = bd[3] / sed[3]
    print(f"  vs {tag}")
    print(f"     q_real = {q_r:+.4f}   q_AD = {q_a:+.4f}   Δq = {bd[3]:+.4f} ± {sed[3]:.4f}   "
          f"**{abs(z):.2f}σ**  ⟹ {'可分辨' if abs(z) > 2 else '**不可分辨**'}")
    print(f"     同时 Δb_λ = {bd[1]:+.4f} ± {sed[1]:.4f}  ({abs(bd[1]/sed[1]):.2f}σ)   "
          f"[对照：上一轮线性回归给 +0.0919 ± 0.0189 = 4.85σ]")

print()
print("=" * 98)
print("3  判决")
print("=" * 98)
d1 = res["实验（真实 Tc）"][0] - res["AD 世界（逐行 μ*）"][0]
bd1, sed1, _, _ = fit(X2, d1)
d2 = res["实验（真实 Tc）"][0] - res["AD 世界（统一 μ*=0.10）"][0]
bd2, sed2, _, _ = fit(X2, d2)
z1, z2 = abs(bd1[3] / sed1[3]), abs(bd2[3] / sed2[3])
q_real = res["实验（真实 Tc）"][1][3]
q_ad1 = res["AD 世界（逐行 μ*）"][1][3]
q_ad2 = res["AD 世界（统一 μ*=0.10）"][1][3]
same_sign = (q_real < 0) and (q_ad1 < 0)
print(f"""  · AD 世界**也有负曲率**：q_AD = {q_ad1:+.4f}（逐行 μ*） / {q_ad2:+.4f}（统一 μ*=0.10）
    vs 实验 q_real = {q_real:+.4f}。符号{'相同' if same_sign else '不同'}。
  · 配对可分辨性：{z1:.2f}σ（逐行 μ*） / {z2:.2f}σ（统一 μ*）
    ⟹ {'**曲率的大小在实验与 AD 之间可分辨**' if min(z1,z2) > 2 else '**曲率的大小在实验与 AD 之间不可分辨**'}
  · μ* 出处敏感性：两版 q_AD 差 {abs(q_ad1-q_ad2):.4f}，配对 σ 差 {abs(z1-z2):.2f}
    ⟹ {'μ* 的出处**对本结论不承重**' if abs(z1-z2) < 1.0 else 'μ* 出处**承重**，必须回页图复核'}""")
print(f"""
  **这对 §5.3 的头条陈述意味着什么（诚实版）**：
    曲率的**存在与符号**是 AD 族与实验**共有**的 ⟹ 不能宣称「只有实验才有曲率」。
    但这**不削弱**头条：天花板 `𝔣≤1` 是 analytic，它**先验地**要求任何正确的 Tc 泛函
    （含 AD）都必须在大 λ 处弯曲 —— AD 也弯，恰恰是这条禁令**在两个独立对象上各现一次**。
    头条的可检验内容应当落在**曲率的量**（q 的值）而非它的存在性上。""")
print("=" * 98)
