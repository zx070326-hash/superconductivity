# -*- coding: utf-8 -*-
"""q 的理论侧先验：弱耦合渐近式能不能预言 q<0？（主会话，2026-09-20，从零）

★ 目的：在 f-sample-expansion 的全 Eliashberg λ-扫描回报**之前**，用一个
   完全独立、零求解器的路径先算出 (b,q) 的量级与符号，作为我自己的对照先验。
   这样我不会被它的答案锚定（CLAUDE.md：承重数值必须自己从零重算）。

仓库既有的弱耦合渐近（memory: p4-weak-coupling-secular-series，前因子 0.2529996837）：
    Tc = 0.2530·ω_E·e^{-1/λ}·(1 + 0.358λ)
Einstein 谱上 ω̄₂ = ω_E ⟹
    𝔣(λ) = Tc/(C_∞·√λ·ω̄₂) = (0.2530/C_INF)·e^{-1/λ}·(1+0.358λ)/√λ

⚠ 严格性与适用域（必须随数字同行）：
   该级数是 **λ→0** 的渐近式，在 λ≈2.6 处**没有**合法性保证；
   本脚本只用它取 (b,q) 的**符号与量级**，**不是** q 的理论值。
   真正的理论 q 必须由全 Eliashberg 解给出（已派活）。
   级别：numerical-evidence（且前提=渐近式外推，比一般 numerical-evidence 更弱）。
"""
import sys; sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

C_INF, PREF, B1 = 0.1827262477, 0.2529996837, 0.358
LO, HI = 0.69, 2.59            # AD75 支撑域

def frak_f(lam):
    return PREF/C_INF*np.exp(-1.0/lam)*(1.0+B1*lam)/np.sqrt(lam)

# 解析斜率 dlog𝔣/dlogλ = 1/λ + B1λ/(1+B1λ) − 1/2
def slope(lam):
    return 1.0/lam + B1*lam/(1.0+B1*lam) - 0.5

print("  λ 处的局部斜率 b_loc = dlog𝔣/dlogλ （解析导数）：")
for lam in (0.69, 1.0, 1.5, 2.0, 2.59):
    print(f"    λ={lam:5.2f}   𝔣={frak_f(lam):7.4f}   b_loc={slope(lam):+7.4f}")
print("  ⟹ b_loc 单调**下降** ⟹ log-log 上凹 ⟹ q < 0，**符号被预言**")

# 在支撑域上做与数据同形的二次拟合
lam = np.exp(np.linspace(np.log(LO), np.log(HI), 400))
x, y = np.log(lam), np.log(frak_f(lam))
X = np.vstack([np.ones_like(x), x, x**2]).T
a, b, q = np.linalg.lstsq(X, y, rcond=None)[0]
print(f"\n  支撑域 [{LO}, {HI}] 上二次拟合： b = {b:+.4f}   q = {q:+.4f}")
print(f"  对照 29 行实测：            b = +0.9972 ± 0.0405   q = −0.525")
print(f"  ⟹ 符号一致；q 量级比 {q/-0.525:.2f}×")

# 𝔣>1 报警（CLAUDE.md：当场报警，不许平滑掉）
mx = frak_f(lam).max()
print(f"\n  支撑域内 max 𝔣 = {mx:.4f}  "
      f"{'★ 警报：𝔣>1，天花板被破——但这是渐近式外推的伪影，见下' if mx > 1 else '（≤1，未触发）'}")
print("  ⟹ 若 >1：正是该渐近式在大 λ 失效的证据（它没有 𝔣→1 的饱和机制），")
print("     **不是**天花板定理出问题；这也说明本脚本只能取符号，不能取数值。")

# --- 限制到渐近式自身仍自洽的区间（𝔣<1）再拟一次 -------------------------
from scipy.optimize import brentq
lam_star = brentq(lambda L: frak_f(L)-1.0, 0.5, 5.0)
print(f"\n  渐近式破天花板的位置： λ* = {lam_star:.4f}  ⟹ 它只在 λ ≲ {lam_star:.2f} 自洽")
l2 = np.exp(np.linspace(np.log(LO), np.log(lam_star), 400))
x2, y2 = np.log(l2), np.log(frak_f(l2))
X2 = np.vstack([np.ones_like(x2), x2, x2**2]).T
a2, b2, q2 = np.linalg.lstsq(X2, y2, rcond=None)[0]
print(f"  自洽区间 [{LO}, {lam_star:.2f}] 上： b = {b2:+.4f}   q = {q2:+.4f}")
print(f"\n  ★ 方向性启发（**不是**界，只是方向）：真曲线必须从下方趋近 𝔣→1，")
print(f"     而渐近式在 λ*={lam_star:.2f} 就冲过 1 且继续上升 ⟹ 真曲线在 log-log 上")
print(f"     弯得**更狠** ⟹ 预期 |q_true| > {abs(q2):.3f}。实测 0.525 与此**同向**。")
print(f"     ⚠ 这是启发不是证明：曲率可能集中在窗口之外。")
