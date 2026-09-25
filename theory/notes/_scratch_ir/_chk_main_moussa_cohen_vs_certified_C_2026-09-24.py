# -*- coding: utf-8 -*-
"""优先权判决 — Moussa-Cohen 2006 eq.(11) 的 f(λ,μ*) vs 我们认证的 C(λ,μ*)。

一手来源：arXiv:cond-mat/0607832 = PRB 74, 094520 (2006)，本机 `_mc2006.pdf` 第 2 页**页图**。
  ⚠ `pdftotext` 对该页**不可用**：它吞掉全部根号（eq.(7) 印成 `0.183 λ ω2`）、
     并把 eq.(11) 的分式塌成竖排 ⟹ 层次不可判。逐字式只能由页图定。

eq.(7)  0 ≤ Tc ≤ 0.183 √(λ⟨ω²⟩)
eq.(10) f(λ⟨ω²⟩/ω²_max, μ*)·ω_max ≤ Tc ≤ f(λ,μ*)·√⟨ω²⟩
eq.(11) f(λ,μ*) = 0.69 exp(−(1+λ)/(λ−μ*)) · √((1+0.52λ)/(1+5.6μ*))   for λ>μ*, else 0
        「empirically fit」；最小二乘区间 0.4≤λ≤20, 0≤μ*≤0.2；RMS 偏差 0.01、最大 0.03。

代数桥：Tc ≤ f(λ,μ*)√⟨ω²⟩ 而 √⟨ω²⟩ = √(λ⟨ω²⟩)/√λ = √(M2)/√λ
   ⟹ Tc ≤ [f(λ,μ*)/√λ]·√(M2)  ⟹  **C_MC(λ,μ*) := f(λ,μ*)/√λ** 与我们的 C(λ,μ*) 同一对象。

方向（承重）：两列都是**上界系数**，越小越紧。
  C_MC < C_ours ⟹ 他们数值上更紧（但 f 是拟合，不是界）
  C_MC > C_ours ⟹ 我们既更紧又更严格
严格性：我们那列 `certificate`（转引 C_lambda_certified_curve_2026-09-24.md §1 的 C_upper）；
       C_MC 列是**对已发表拟合式的求值**，是 numerical-evidence，且它本身不是声音上界。
"""
import os
import sys
import math
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

C_INF = 0.1827262477

# 我们的认证 C_upper（逐字抄自 C_lambda_certified_curve_2026-09-24.md §1）
C_UP = {
    0.00: {2.0: 0.151236028, 3.2953: 0.164718304, 4.0: 0.168304891, 6.0: 0.173820789},
    0.13: {1.0: 0.085109595, 1.5: 0.101992146, 2.0: 0.116687205, 2.5: 0.125402058,
           3.0: 0.130991149, 3.2953: 0.133361200, 4.0: 0.137762123, 5.0: 0.141888517,
           6.0: 0.144677171, 8.0: 0.148030265, 10.0: 0.150564456},
    0.20: {2.0: 0.104372372, 3.2953: 0.122242127, 4.0: 0.126754407, 6.0: 0.133713798},
}


def f_mc(lam, mus):
    """Moussa-Cohen eq.(11)，逐字。"""
    if lam <= mus:
        return 0.0
    return 0.69 * math.exp(-(1.0 + lam) / (lam - mus)) * math.sqrt((1.0 + 0.52 * lam) / (1.0 + 5.6 * mus))


def c_mc(lam, mus):
    return f_mc(lam, mus) / math.sqrt(lam)


# ── 自查 0：他们自己的两式必须互洽（eq.(11) 的 λ→∞ 应回到 eq.(7) 的 0.183）──
print("=" * 104)
print("自查 0 — eq.(11) 的强耦合极限是否回到 eq.(7) 的 0.183？（读对了层次的充分条件）")
for lam in (1e2, 1e4, 1e6, 1e8):
    print(f"   λ={lam:9.0e}  μ*=0:  f/√λ = {c_mc(lam, 0.0):.9f}"
          f"   |差 0.183| = {abs(c_mc(lam,0.0)-0.183):.2e}")
lim = 0.69 * math.exp(-1.0) * math.sqrt(0.52)
print(f"   解析极限 0.69·e⁻¹·√0.52 = {lim:.9f}   (他们印的 eq.(7) 系数 = 0.183)")
print(f"   ⟹ 读法自洽 ✓（若把 √ 读在别处，此处必不等于 0.183）")

# ── 自查 1：他们的 eq.(11) 与我们的 analytic C_∞ 谁大 ──
print("\n" + "=" * 104)
print("自查 1 — C_MC(λ,0) 与我们 analytic 的 C_∞ = (1/2π)√𝔤(2) 的关系")
print(f"   我们的 C_∞ = {C_INF:.10f}（analytic，λ 与 μ* 全盲）")
print(f"   他们的 0.183（eq.7）= C_MC 的 λ→∞ 极限 = {lim:.10f}")
print(f"   相对差 = {(lim-C_INF)/C_INF*100:+.3f}%   ⟹ **同一个常数**，他们印到 3 位有效数字。")

# ── 主表：逐格比较 ──
print("\n" + "=" * 104)
print("主表 — C_MC(λ,μ*) = f/√λ  vs  我们认证的 C_upper(λ,μ*)      [两列都是上界系数，越小越紧]")
for mus in (0.00, 0.13, 0.20):
    print(f"\n  μ* = {mus:.2f}")
    print(f"   {'λ':>8}{'C_MC (拟合)':>16}{'C_ours (认证上端)':>20}{'C_MC/C_ours':>14}   判决")
    for lam in sorted(C_UP[mus]):
        cm, co = c_mc(lam, mus), C_UP[mus][lam]
        r = cm / co
        if r < 1.0:
            verd = f"他们紧 {1/r:.3f}×  ⚠（但 f 是拟合，非界）"
        else:
            verd = f"**我们紧 {r:.3f}×**，且我们是 certificate"
        print(f"   {lam:8.4f}{cm:16.9f}{co:20.9f}{r:14.4f}   {verd}")

# ── 他们的拟合误差能否吃掉差距？──
print("\n" + "=" * 104)
print("误差预算 — 他们自陈 f 的 RMS 偏差 0.01、最大偏差 0.03（绝对，作用在 f 上）")
print("   把最大偏差朝**对他们有利**（f 变小 ⟹ 界变紧）的方向退化：C_MC^best = (f−0.03)/√λ")
print(f"   {'λ':>8}{'μ*':>7}{'C_MC':>13}{'C_MC^best':>13}{'C_ours':>13}   最有利情形下谁紧")
for mus in (0.00, 0.13):
    for lam in sorted(C_UP[mus]):
        cm = c_mc(lam, mus)
        best = max(f_mc(lam, mus) - 0.03, 0.0) / math.sqrt(lam)
        co = C_UP[mus][lam]
        who = "他们" if best < co else "**我们**"
        print(f"   {lam:8.4f}{mus:7.2f}{cm:13.9f}{best:13.9f}{co:13.9f}   {who}")

print("\n" + "=" * 104)
print("★ 判词（承重，两条分开写）")
print("  (1) **对象相同**：C_MC = f/√λ 与我们的 C(λ,μ*) 是同一个量（Einstein 谱在固定 (λ,⟨ω²⟩) 下的 Tc）。")
print("      他们把该极值性归给 ref.3 = Allen-Dynes 1975（原文：'It is already known³ that the")
print("      Einstein spectrum maximizes Tc for a given λ and ⟨ω²⟩'）⟹ **定理内容不是他们的，也不是我们的**。")
print("  (2) **严格性不同**：他们的 f 是 0.4≤λ≤20 上的最小二乘拟合（自陈 RMS 0.01/max 0.03），")
print("      拟合值**不构成上界**——偏差可双向。我们的 C_upper 是区间算术的外向包封。")
print("      ⟹ 我们剩下的增量 = 『把一个拟合换成认证包封』+『μ*>0 时极值性前提本身的处理』，")
print("         **不是**『首次给出该上界』。")
