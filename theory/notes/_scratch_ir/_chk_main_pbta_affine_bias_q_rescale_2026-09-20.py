# -*- coding: utf-8 -*-
"""
补 Q1 漏掉的情形（我自己的订正，未发出即捕获）：
上一版只推了**常数**倍率偏置。但 Pb/Ta 实测 `dlogρ/dlogλ = −0.173` ⟹ 偏置**不是常数**。
若 logρ 是 logλ 的**仿射**函数，q 还不动吗？

  logρ = r0 + r1·L_T ，logσ = s0 + s1·L_T          （L := logλ）
  ⟹ L_D = L_T + logρ = (1+r1)L_T + r0              ← L 轴的**仿射重标度**
     y_D = y_T − ½logρ − logσ = [a+bL_T+qL_T²] − 仿射(L_T)
  以 L_D 为自变量回代 L_T = (L_D − r0)/(1+r1)：
     y_D = a″ + b″·L_D + **q/(1+r1)²**·L_D²

★ ⟹ 仿射偏置**不创造也不消灭**曲率（q=0 ⟺ q_D=0），但**按 1/(1+r1)² 重标度**它的大小。
   r1 = 0（常数偏置）时因子 = 1 ⟹ 上一版结论是本结论的特例，没错，但**不完整**。

严格性：符号推演 analytic；数值确认 numerical-evidence。主会话 2026-09-20。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

rng = np.random.default_rng(7)
L = np.log(rng.uniform(0.7, 2.6, 4000)); S = np.log(rng.uniform(0.55, 0.89, 4000))
Q_TRUE = -0.525
y = 0.3 + 1.0*L + 0.7*S + Q_TRUE*L**2                    # 无噪声：测恒等式而非统计

def fit(Lv, yv, Sv):
    A = np.column_stack([np.ones_like(Lv), Lv, Sv, Lv**2])
    return np.linalg.lstsq(A, yv, rcond=None)[0]

print("="*90)
print("仿射偏置 logρ = r0 + r1·L ，logσ = s0 + s1·L  ⟹  q_D = q_T/(1+r1)² ?")
print("="*90)
print(f"  {'r1':>8}{'s1':>8}{'q_拟合':>12}{'q/(1+r1)² 预言':>18}{'相对差':>12}")
for r1, s1 in [(0.0, 0.0), (-0.1726, -0.0535), (-0.1726, 0.0), (0.30, -0.20), (-0.40, 0.0)]:
    lr = 0.05 + r1*L; ls = -0.02 + s1*L
    c = fit(L + lr, y - 0.5*lr - ls, S)
    pred = Q_TRUE/(1+r1)**2
    print(f"  {r1:>8.4f}{s1:>8.4f}{c[3]:>12.6f}{pred:>18.6f}{abs(c[3]/pred-1):>12.2e}")
print("  ⟹ 逐行吻合到机器精度 ⟹ 符号推演成立。注意 s1（ω̄₂ 侧偏置）**完全不进 q**：")
print("     它只在 y 上加仿射项，不动 L 轴。")

# q=0 对照：仿射偏置不能凭空造出曲率
y0 = 0.3 + 1.0*L + 0.7*S
c0 = fit(L + 0.05 - 0.1726*L, y0 - 0.5*(0.05-0.1726*L) + 0.02, S)
print(f"\n  ★ q_T=0 的对照：仿射偏置后 q_D = {c0[3]:+.3e}  ⟹ 曲率**无中不能生有**。")

print()
print("="*90)
print("代入 Pb/Ta 实测 r1 ⟹ 头条 q 的技术系统误差有多大")
print("="*90)
r1 = np.log((1.68/1.55)/(0.86/0.69))/np.log(1.55/0.69)
q_rep = -0.525
print(f"  实测 r1 = dlogρ/dlogλ = {r1:+.4f}   (Pb/Ta 两点，⚠ Savrasov 值未经我核验)")
print(f"  重标度因子 1/(1+r1)² = {1/(1+r1)**2:.4f}")
print(f"  ⟹ q: {q_rep:+.3f} → {q_rep/(1+r1)**2:+.3f}   （变 {100*(1/(1+r1)**2-1):+.1f}%）")
print(f"  对照 q 的跨设定浮动区间 [−0.63, −0.40]（笔记 §4′ 第 2 条）：")
lo, hi = -0.63, -0.40
print(f"    技术订正后的点值 {q_rep/(1+r1)**2:+.3f} "
      f"{'仍落在' if lo <= q_rep/(1+r1)**2 <= hi else '**已落到**'} 该区间{'内' if lo <= q_rep/(1+r1)**2 <= hi else '外'}")

print("""
  ★ 订正后的预算判决（取代上一版「对 q 预算 = ∞」）：
    · **q ≠ 0（曲率存在性）**：预算 = ∞ —— 仿射技术偏置**无中不能生有** [analytic]
      ⟹ 头条的**定性**内容（𝔣 朝天花板弯曲）对技术系统误差**结构免疫**，不需要 Pb/Ta。
    · **q 的数值**：预算**有限**—— 按实测 r1 重标度 46%，直接跨出现报区间
      ⟹ Pb/Ta **买得到**这条，消费点 = 「q 报区间 [−0.63,−0.40]」那一行。
    · **ω̄₂ 数字化（σ 侧）单独看**：`s1` **完全不进 q** [analytic] ⟹ 对 q 预算 = ∞；
      它只买 `b` 和 𝔣 的绝对值。**⟹ 要保 q，只需 λ 的两技术比值，λ 已印在两表里，不必数字化。**""")
print("="*90)
