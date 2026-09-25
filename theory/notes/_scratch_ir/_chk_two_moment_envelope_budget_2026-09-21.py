#!/usr/bin/env python
"""预算：两矩包络 sup{Tc : int a2F dw = A, int a2F w dw = M1} 比 min(Leavens, moment) 紧多少?

背景（2026-09-21 查清的优先权链）:
  Leavens 1975 SSC 17,1499 : kTc <= 0.2309 A            (0 阶矩; 等号 A delta(w-1.750A))
  Leavens 1976 SSC 19,395  : 宣称 AD 的 sqrt(lam<w^2>) 与上式不兼容、"in fact incorrect"
  LAD    1979 SSC 30,595   : kTc ~= 0.1827 sqrt(lam<w^2>)   (1 阶矩; M1 = lam<w^2>/2)
两者都只用**一个**矩。用 (I₀,I₁) 这一对的包络未见于以下三篇
（Moussa-Cohen 2006 PRB 74,094520 的「两个界」是掺杂限与晶格失稳限, 非矩界
 —— 2026-09-24 一手全文确认本句**成立**, 见下方第二次订正）。

⚠⚠ **2026-09-24 第二次订正（团队长）——上一条标注本身是错的，现予反注销。**
2026-09-24 取回 MC2006 一手全文后，标题/摘要/正文三处互证：
  标题 `Two Bounds on the Maximum Phonon-Mediated Superconducting Transition Temperature`
  摘要 `The first bound is set by the total electron-phonon coupling available … arbitrarily **dope**
        the material. … The second bound is set by the **softening of phonons to instability**.`
  §I末 `two fundamental bounds …, one set by **lattice stability** and the other by the total strength
        of available electron-phonon coupling`
⟹ **第 9 行那句原话（掺杂限 + 晶格失稳限）是对的**，2026-09-23 的注销为假（它 grep 的是**仓库**，
   不是论文——当时手上没有全文）。下面那段 2026-09-24 的第一版标注跟着错了，是同一链的第三层。
   eq.(7)/(9)/(10)/(11) 的转录本身准确，但它们是**中段引理**：§V 逐字
   `using η(E), Eq. (19), and Eq. (10) to tune to the largest Tc` ⟹ eq.(10) 是第一条界的**工具**。

**对本脚本新颖性前提的重判（这条才承重）**：
  本脚本的包络用的是 **(I₀, I₁) = (∫α²F dω, ∫α²F·ω dω)**；
  MC eq.(9)(10) 用的是 **(I₋₁, I₊₁, ω_max) = (λ, λ⟨ω²⟩, 支撑上端)**。
  ⟹ **不同的信息类，MC 不构成对本脚本的反例。**
  但必须同时承认：MC eq.(9) 是「固定 λ、⟨ω²⟩、ω_max 下**最小化** Tc 的谱」的显式解，
  ⟹ 「约束不止一个矩」这个**做法**是他们的，不是我们的。
  ⟹ 本脚本可写的新颖性只剩窄义一条：**(I₀,I₁) 这一对**的包络未见于 MC2006 / Leavens / LAD。
     「无人做过」是**否定存在命题**，只有穷尽检索才能证；当前状态 = **「未见反例」，不是「成立」**。
  （2026-09-24 起本文件被 `theory/notes/ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md`
   §0/§6 点名引用 ⟹ 按毕业规则已 `git add -f` 入库。上一版此处写「未入库、未被引用」，现已过期。）

~~⚠⚠ **2026-09-24 注销上面那句括号**（团队长）。它是 2026-09-21 写下的、**无出处**的陈述，
并已在 `theory/notes/literature/einstein_extremality_and_eta_criterion_priority_2026-09-23.md`
问题 3 里被明确注销……~~ **（整段作废，见上。）**

本脚本只回答预算问题：包络比 min 紧多少？若 ~0 则不派活（CLAUDE.md 预算优先）。

算法（比嵌套 Nelder-Mead-in-bisection 快约 20x）:
  对固定 T, 在约束集上**最大化 top eigenvalue**; 它对 T 递减（每个固定谱的 top_eig 对 T
  递减, 递减函数族的上确界仍递减）==> 外层对 T 单次二分即得包络 Tc。

方向审计：试探族只取有限原子数 ==> 得到的是 sup 的**下估** env_lb <= sup <= min(两界)。
  因此本脚本报的 gap% = 100*(1 - env_lb/min) 是包络**可能收益的上限**：
  真收益 = 100*(1 - sup/min) <= gap%。
  ==> gap 小 ==> 收益必小 ==> 不派活（判决可靠）。
      gap 大 ==> 值得派活, 但**不构成定理**（需要 sup 的上界证明）。
  这是单边界方向规则的正确用法：本脚本只能产出「不派活」的硬判决。
等级：numerical-evidence（float64，有限 N，有限原子数，局部优化 + 二分）。
"""
import math

import numpy as np
from scipy.optimize import minimize

C_INF = 0.1827262477890
C_LEAVENS_PUB = 0.2309        # Leavens 1975 印刷值（外部校验用, 不作为输入）


def top_eig(A_i, w_i, T, N):
    """AD75 对称矩阵, mu*=0, 离散谱 sum A_i delta(w - w_i).
    lam(l) = sum_i 2 A_i w_i / (w_i^2 + (2 pi T l)^2)."""
    nu = 2.0 * math.pi * T * np.arange(0, 2 * N + 2)
    nu2 = nu * nu
    lam_l = np.zeros_like(nu)
    for a, w in zip(A_i, w_i):
        lam_l += 2.0 * a * w / (w * w + nu2)
    n = np.arange(N)
    M = lam_l[np.abs(n[:, None] - n[None, :])] + lam_l[n[:, None] + n[None, :] + 1]
    cum = np.concatenate(([0.0], np.cumsum(lam_l[1:N])))
    M[np.diag_indices(N)] -= (2 * n + 1) + lam_l[0] + 2 * cum
    return np.linalg.eigvalsh(M)[-1]


def tc(A_i, w_i, N=400, hi=5.0):
    lo = 1e-7
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        if top_eig(A_i, w_i, mid, N) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ================================================================ 0. solver blind test
print("=" * 80)
print("0. 求解器盲测：单原子谱必须同时重现 Leavens 0.2309/1.750 与 LAD 0.1827")
print("=" * 80)
best = (None, -1.0)
for wE in np.arange(1.2, 2.4, 0.05):
    t = tc([1.0], [wE], 250)
    if t > best[1]:
        best = (wE, t)
lo, hi = best[0] - 0.05, best[0] + 0.05
for _ in range(30):
    a, b = lo + (hi - lo) / 3, hi - (hi - lo) / 3
    if tc([1.0], [a], 250) < tc([1.0], [b], 250):
        lo = a
    else:
        hi = b
wE_opt = 0.5 * (lo + hi)
c0 = tc([1.0], [wE_opt], 700)
print(f"   max_wE Tc/A = {c0:.6f}  at w_E/A = {wE_opt:.4f}"
      f"   [Leavens 1975 印: 0.2309 at 1.750]")
print(f"   deviation   : c {100*abs(c0-C_LEAVENS_PUB)/C_LEAVENS_PUB:.3f}%"
      f"   w_E {100*abs(wE_opt-1.750)/1.750:.3f}%")
for lam in (5.0, 20.0, 80.0):
    wE = math.sqrt(1.0 / lam)                     # M1 = lam wE^2/2 = 1/2
    print(f"   lam={lam:5.1f}: Tc/sqrt(lam<w^2>) = "
          f"{tc([lam*wE/2], [wE], 900, hi=1.0):.6f}   [LAD 1979: 0.1827]")
lam_star = (2 * C_INF / c0) ** 2
print(f"\n   crossover  lam* = (2 C_inf / c0)^2 = {lam_star:.4f}")
print()

# ================================================================ 1. envelope budget
print("=" * 80)
print("1. 预算：固定 (A=1, M1=m)，k-原子试探族包络  vs  min(Leavens, moment)")
print("=" * 80)
print(f"{'m=M1/A':>7s} {'lam_1atom':>9s} {'1-atom':>8s} {'2-atom':>8s} {'3-atom':>8s} "
      f"{'Leavens':>8s} {'moment':>8s} {'min':>8s} {'gap%':>7s}")


def env_tc(m, k, N=260):
    """max Tc over k-atom spectra with sum A_i = 1, sum A_i w_i = m."""
    def max_eig_at(T):
        def neg(p):
            a = np.append(p[: k - 1], 0.0)
            a[-1] = 1.0 - a[: k - 1].sum()
            w = np.empty(k)
            w[: k - 1] = p[k - 1:]
            if (a <= 1e-9).any() or (w[: k - 1] <= 1e-6).any():
                return 1e3
            w[-1] = (m - (a[: k - 1] * w[: k - 1]).sum()) / a[-1]
            if w[-1] <= 1e-6:
                return 1e3
            return -top_eig(a, w, T, N)
        if k == 1:
            return top_eig([1.0], [m], T, N)
        starts = ([0.5, 0.6 * m], [0.25, 0.35 * m], [0.75, 1.6 * m]) if k == 2 else \
                 ([0.34, 0.33, 0.4 * m, 1.0 * m], [0.2, 0.5, 0.25 * m, 0.9 * m])
        return -min(minimize(neg, x0, method="Nelder-Mead",
                             options={"xatol": 2e-4, "fatol": 1e-9, "maxiter": 300}).fun
                    for x0 in starts)
    lo, hi = 1e-7, 3.0
    for _ in range(34):
        mid = 0.5 * (lo + hi)
        if max_eig_at(mid) > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


for m in (0.4, 0.6, 0.7986, 1.0, 1.5, 2.5, 5.0):
    t1 = env_tc(m, 1)
    t2 = env_tc(m, 2)
    t3 = env_tc(m, 3)
    leav, mom = c0, C_INF * math.sqrt(2.0 * m)
    mn = min(leav, mom)
    best_t = max(t1, t2, t3)
    print(f"{m:7.4f} {2/m:9.3f} {t1:8.5f} {t2:8.5f} {t3:8.5f} "
          f"{leav:8.5f} {mom:8.5f} {mn:8.5f} {100*(1-best_t/mn):7.2f}")

print()
print("   gap% = 100*(1 - max(trial)/min(bound))  = 两矩包络**可能**收益的上限。")
print("   3-atom <= 2-atom 说明极值谱是 2 原子（moment-problem 的典型结构）。")
