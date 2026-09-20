# -*- coding: utf-8 -*-
"""
主会话从零复核 `f-sample-expansion` 的四材料 `r1` 估计（2026-09-20 晚）。

被复核对象：r1 = 0.119 ± 0.309，R²=0.069，q: -0.525 -> -0.419，建议扩至 >=8 材料。

本脚本回答四问：
  Q1  他们的 OLS 算术对不对（我从零重算，不转引）。
  Q2  n=4 的 r1 有没有内容？（LOO + Ta 双源敏感性 + 1σ 带宽换算成 q 的带宽）
  Q3  扩到 n=8 买得到什么？（把「q 可引用」写成对 σ_r1 的显式要求，反解所需 n）
  Q4  ★ 新产出：`q` 的**符号**对仿射技术偏置是**无条件**免疫的（(1+r1)² 是完全平方），
      而**线性系数 `b` 不免疫、且 `s1` 会进 `b`** ⟹ 我先前「s1 不进 q ⟹ 数字化 ω̄₂
      预算 = ∞」那条裁决**不得外推到 b**。用无噪声恒等式设定数值确认。

另附 Q5：Pb 的 ω̄₂=65 K 撞车问题——他们给的「第三列是 (ω̄₂/ω_log)² 」解释自相矛盾，本脚本算给看。

方向：Q1–Q3 是双侧估计，不产生单边界判决；Q4 是恒等式。
严格性：Q1/Q2/Q3 numerical-evidence；Q4 analytic（代数）+ 无噪声数值确认。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

C_INF = 0.1827262477

# ---------------------------------------------------------------- Q1
# AD75 Table I（隧道）vs Savrasov Table II（DFT），agent 逐项对页图后给出
PAIRS = {           # name: (lam_tun, lam_dft)
    "Pb": (1.55, 1.68),
    "Ta": (0.69, 0.86),
    "Al": (0.43, 0.44),
    "Mo": (0.60, 0.42),
}

def ols(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x)
    dx = x - x.mean(); dy = y - y.mean()
    b = (dx*dy).sum()/(dx*dx).sum()
    a = y.mean() - b*x.mean()
    sse = ((y - (a + b*x))**2).sum()
    r2 = 1.0 - sse/(dy*dy).sum()
    sb = np.sqrt(sse/max(n-2, 1)/(dx*dx).sum()) if n > 2 else np.nan
    return a, b, sb, r2

def fit_r1(pairs):
    L = [np.log(v[0]) for v in pairs.values()]
    R = [np.log(v[1]/v[0]) for v in pairs.values()]
    return ols(L, R)

print("="*86)
print("Q1  四点 OLS：log(rho) = r0 + r1 log(lam_tun)   —— 我从零重算")
print("="*86)
print(f"  {'材料':>6}{'lam_tun':>10}{'lam_dft':>10}{'rho':>10}{'log lam':>12}{'log rho':>12}")
for k, (lt, ld) in PAIRS.items():
    print(f"  {k:>6}{lt:>10.2f}{ld:>10.2f}{ld/lt:>10.4f}{np.log(lt):>12.6f}{np.log(ld/lt):>12.6f}")
a4, r1_4, s4, r2_4 = fit_r1(PAIRS)
print(f"\n  r1 = {r1_4:+.6f}  +/- {s4:.6f}    R^2 = {r2_4:.4f}")
print(f"  他们报 r1 = +0.119 +/- 0.309, R^2 = 0.069  ->  "
      f"{'逐位相符 OK' if abs(r1_4-0.119)<5e-4 and abs(s4-0.309)<5e-4 else '不符!!'}")

# ---------------------------------------------------------------- Q2
print()
print("="*86)
print("Q2  n=4 的 r1 有没有内容？")
print("="*86)
print("  (a) 留一出（LOO）——每去掉一个材料，r1 变多少：")
names = list(PAIRS)
loo = {}
for drop in names:
    sub = {k: v for k, v in PAIRS.items() if k != drop}
    _, b, _, _ = fit_r1(sub)
    loo[drop] = b
    print(f"      去掉 {drop:>3}  ->  r1 = {b:+.6f}")
print(f"      我先前的两点（Pb+Ta）      ->  r1 = {loo['Al'] if False else fit_r1({k:PAIRS[k] for k in ('Pb','Ta')})[1]:+.6f}")
lo, hi = min(loo.values()), max(loo.values())
print(f"      LOO 跨度 [{lo:+.4f}, {hi:+.4f}]，宽 {hi-lo:.4f} = {(hi-lo)/abs(r1_4):.2f} 倍点估计")

print("\n  (b) Ta 的隧道值双源（AD75 0.69 vs Savrasov 引 Ref.3 的 0.78）：")
alt = dict(PAIRS); alt["Ta"] = (0.78, 0.86)
_, r1_alt, s_alt, _ = fit_r1(alt)
print(f"      Ta=0.69 -> r1 = {r1_4:+.6f}   |   Ta=0.78 -> r1 = {r1_alt:+.6f}")
print(f"      该源头分歧对 r1 的影响 = {abs(r1_alt-r1_4):.4f}"
      f" = {abs(r1_alt-r1_4)/s4*100:.1f}% 的 1σ")
print("      ⟹ 「隧道内部 13% 分歧」在单材料 rho 上属实，但对 **r1** 的传播是次主导；")
print("         主导项是四点自身的散布（R^2=0.07）。这是两个不同的陈述。")

print("\n  (c) 把 r1 的 1σ 带换算成 q 的带（消费点：q_D = q_T/(1+r1)^2）：")
Q_T = -0.525
REPORT = (-0.63, -0.40)         # 跨设定浮动区间
for lab, rr in [("点估计", r1_4), ("-1σ", r1_4-s4), ("+1σ", r1_4+s4)]:
    fac = 1.0/(1.0+rr)**2
    print(f"      {lab:>7}: r1={rr:+.4f}  (1+r1)^-2 = {fac:.4f}   q_D = {Q_T*fac:+.4f}")
f_lo = 1.0/(1.0+r1_4+s4)**2; f_hi = 1.0/(1.0+r1_4-s4)**2
print(f"      ⟹ q_D ∈ [{Q_T*f_hi:+.4f}, {Q_T*f_lo:+.4f}]，跨度比 {f_hi/f_lo:.2f}×")
print(f"      对照：报告区间 [{REPORT[0]}, {REPORT[1]}] 跨度比 {REPORT[0]/REPORT[1]:.2f}×")
print(f"      ★ 技术订正的不确定度比它要去修正的那个跨设定浮动还宽 "
      f"{(f_hi/f_lo)/(REPORT[0]/REPORT[1]):.2f} 倍 ⟹ **q 的数值仍不可引用**（与 n=2 时同一判决）")

# ---------------------------------------------------------------- Q3
print()
print("="*86)
print("Q3  扩到 n=8 买得到什么？（把「q 可引用」写成对 σ_r1 的显式要求）")
print("="*86)
span_req = REPORT[0]/REPORT[1]           # 要求技术带 <= 跨设定带
# 解 ((1+r+s)/(1+r-s))^2 <= span_req
rt = np.sqrt(span_req)
s_req = (rt-1.0)*(1.0+r1_4)/(rt+1.0)
print(f"  要求：((1+r1+σ)/(1+r1-σ))^2 <= {span_req:.3f}  ⟹  σ_r1 <= {s_req:.4f}")
print(f"  现有 σ_r1 = {s4:.4f} ⟹ 需缩 {s4/s_req:.2f}×")
print(f"  若散布是纯抽样噪声（σ ∝ n^-1/2、λ 覆盖不变）：n >= 4×({s4/s_req:.2f})^2 = "
      f"**{4*(s4/s_req)**2:.0f} 个同材料双技术对**")
for n in (6, 8, 12, 16, 22, 30):
    sn = s4*np.sqrt(4.0/n)
    sp = ((1+r1_4+sn)/(1+r1_4-sn))**2
    ok = "可引用" if sp <= span_req else "仍不可引用"
    print(f"      n={n:>3}:  σ≈{sn:.4f}   q 带跨度比 {sp:5.2f}×   {ok}")
print("  ⚠ 上表假设散布**全部**是抽样噪声。Ta 的 0.69/0.78 双源说明 rho 本身带")
print("    材料级固有散布 ⟹ 那一份**不随 n 缩小**，上表的 n 是**乐观下界**。")

# ---------------------------------------------------------------- Q4
print()
print("="*86)
print("Q4  ★ q 的符号无条件免疫；b 不免疫、且 s1 会进 b")
print("="*86)
print("  q_D = q_T/(1+r1)^2。分母是**完全平方** ⟹ sign(q_D) = sign(q_T)，对任何 r1 != -1。")
print("  ⟹ 不只是「无中不能生有曲率」，**弯曲的方向也改不了** [analytic]。")
print(f"     （r1 = -1 是轴塌缩的奇点；实测 r1 = {r1_4:+.4f} ± {s4:.4f} 距它 "
      f"{abs(-1-r1_4)/s4:.1f}σ）")

rng = np.random.default_rng(20260920)
n = 40
L = rng.uniform(-0.9, 1.0, n)          # log lam
S = rng.uniform(-0.6, 0.6, n)          # log S
A0, B0, C0 = 0.3, 1.0, 0.7
y = A0 + B0*L + C0*S + Q_T*L**2        # 无噪声：测恒等式而非统计

def refit(r0, r1, s0, s1):
    lr = r0 + r1*L
    ls = s0 + s1*L
    LD = L + lr                         # L_D = (1+r1)L + r0
    yD = y - 0.5*lr - ls
    M = np.column_stack([np.ones(n), LD, S, LD**2])
    return np.linalg.lstsq(M, yD, rcond=None)[0]   # [A, b_D, c, q_D]

print(f"\n  {'r0':>7}{'r1':>8}{'s0':>7}{'s1':>8}{'q_D':>12}{'q预测':>12}{'b_D':>12}{'b预测':>12}")
for r0, r1, s0, s1 in [(0.0, 0.0, 0.0, 0.0),
                       (0.05, -0.1726, -0.02, -0.0535),
                       (0.05, -0.1726, -0.02, 0.0),
                       (0.03, r1_4, -0.01, -0.05),
                       (0.03, r1_4, -0.01, 0.0),
                       (0.10, 0.40, 0.05, 0.30)]:
    A, bD, c, qD = refit(r0, r1, s0, s1)
    u = 1.0/(1.0+r1)
    q_pred = Q_T*u*u
    b_pred = (B0 - 0.5*r1 - s1)*u - 2.0*Q_T*u*u*r0
    print(f"  {r0:>7.2f}{r1:>8.4f}{s0:>7.2f}{s1:>8.4f}"
          f"{qD:>12.6f}{q_pred:>12.6f}{bD:>12.6f}{b_pred:>12.6f}")
    assert abs(qD-q_pred) < 1e-10 and abs(bD-b_pred) < 1e-10, "闭式与拟合不符"
print("  （六行全部逐位吻合闭式，断言门 1e-10 通过）")
print("\n  ★ 对照第 2/3 行与第 4/5 行：s1 改变时 **q_D 逐位不动、b_D 动**。")
print("     ⟹ 「s1 不进 q ⟹ 数字化 ω̄₂ 的预算 = ∞」这条裁决**只对 q 成立**，")
print("        **不得外推到 b（=0.997±0.041，「𝔣 ∝ λ^1」那条）**。b 要用作承重量，")
print("        就必须同时估 r0、r1、s1 三个偏置参数 —— 那才是数字化 ω̄₂ 的真消费点。")

# ---------------------------------------------------------------- Q5
print()
print("="*86)
print("Q5  他们对「Pb 的 65 K 撞车」的解释：(ω̄₂/ω_log)² = 65 ?")
print("="*86)
wlog, wb2_claim, wb2_repo, lam, tc = 56.0, 60.0, 65.0, 1.55, 7.2
print(f"  他们称 AD75 Pb 行为 ω_log={wlog:.0f}, ω̄₂={wb2_claim:.0f}, 且「第三列 = (ω̄₂/ω_log)² = 65」")
print(f"  实算 (60/56)^2 = {(wb2_claim/wlog)**2:.4f}  —— **不是 65**，差 {65/((wb2_claim/wlog)**2):.0f} 倍")
print(f"  且表头逐字为 `w_log(K) | wbar1(K) | wbar2(K) | lambda | ...` ⟹ 第三列是**频率**，量纲为 K")
print(f"\n  仓库既有转录 ω̄₂ = {wb2_repo:.0f} K 的验收门（𝔣_Pb 必须复现 0.4869）：")
for w in (wb2_repo, wb2_claim):
    f = tc/(C_INF*np.sqrt(lam)*w)
    print(f"      ω̄₂ = {w:.0f} K  ->  𝔣_Pb = {f:.4f}"
          f"   {'<= 复现仓库值 OK' if abs(f-0.4869)<5e-4 else '<-- 与仓库值不符'}")
print("\n  ⟹ 撞车的**正确**排除方式是他们自己第二节已经做到的那条：")
print("     DFT ω_log 对隧道 ω_log 为 Pb +16.1% / Ta +21.2%，**同向同量级** ⟹ 无串列；")
print("     65 K 的相撞是巧合（Savrasov 的 Pb ω_log 与 AD75 的 Pb ω̄₂ 是两个不同的量）。")
print("="*86)
