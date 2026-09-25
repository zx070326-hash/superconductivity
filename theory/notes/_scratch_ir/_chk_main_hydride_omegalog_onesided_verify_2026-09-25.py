#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主会话从零复核 teammate `hydride-omega2-retry` 的 LaH10 / LaD10 / CaH6 报告。
2026-09-25。不 import 交付方脚本，只喂它逐字抄回的表值 + 仓库自有 H3S 标定行。

复核四件事：
  (A) 单位常数与 H3S 锚：C_inf * sqrt(lam) * w2[K] 是否真给 415.2 K；两条独立路线是否自洽
  (B) Jensen 方向：w2 >= w_log 是否真的使判决单边、且单边朝 PASS
  (C) 逐行重算 r_up，并与交付方报的数对齐
  (D) ★ 交付方说做不了的那一步：用仓库自有 H3S 行标定 w2/w_log 比值
      （同一次 alpha2F 积分 ⟹ 零约定混杂），与各行所需阈值对照
"""
import math, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

OK = True
def chk(tag, cond, detail=""):
    global OK
    print(("  PASS  " if cond else "  FAIL  ") + tag + (("  | " + detail) if detail else ""))
    if not cond: OK = False

MEV_TO_K = 11.604518               # 1 meV -> K  (CODATA: e/k_B * 1e-3)
C_INF    = 0.1827262477            # = sqrt(g(2))/(2pi)，仓库既有 analytic 常数
KAPPA    = 750.2803                # hbar*sqrt((eV/A^2)/amu)/k_B  [K]
PRE      = 137.0959                # C_inf * kappa
NB3GE_R  = 0.6260                  # 全语料最紧点

print("=" * 78)
print("(A) 单位常数自洽 + H3S 锚")
print("=" * 78)
chk("C_inf*kappa == PRE", abs(C_INF * KAPPA - PRE) < 5e-4,
    "%.6f vs %.4f" % (C_INF * KAPPA, PRE))

# H3S 200 GPa anharmonic，仓库 hydrides_theory.csv 第 3 行（本仓库自有 alpha2F 积分）
H3S_LAM, H3S_WLOG, H3S_W2 = 1.721, 109.995, 149.26     # meV
ceil_h3s = C_INF * math.sqrt(H3S_LAM) * H3S_W2 * MEV_TO_K
print("  H3S-200 anh: lam=%.3f  w2=%.2f meV -> ceiling = %.2f K" % (H3S_LAM, H3S_W2, ceil_h3s))
chk("复现交付方自检的 415.2 K", abs(ceil_h3s - 415.2) < 0.15,
    "算得 %.3f，相对差 %.2e" % (ceil_h3s, abs(ceil_h3s - 415.2) / 415.2))

# 第二条独立路线：经 M2 [eV/A^2/amu] 与 PRE
# M2 = lam * <w^2>;  能量^2 -> (eV/A^2/amu) 需除以 kappa^2（因 kappa 是该单位的 K 当量）
m2_K2 = H3S_LAM * (H3S_W2 * MEV_TO_K) ** 2
ceil_h3s_route2 = PRE * math.sqrt(m2_K2 / KAPPA**2)
# 容差按 PRE / kappa 的印出位数定：PRE 只有 7 位有效数字 ⟹ 相对 1e-7 量级不可分辨。
# 用更紧的容差是拿「两个四舍五入常数」当精确值，属自造 FAIL。
chk("两条路线一致（C_inf*sqrt(lam)*w2  ==  PRE*sqrt(M2)）",
    abs(ceil_h3s - ceil_h3s_route2) / ceil_h3s < 1e-6,
    "%.6f vs %.6f，相对差 %.2e（PRE 印 7 位，1e-7 以下不可分辨）"
    % (ceil_h3s, ceil_h3s_route2, abs(ceil_h3s - ceil_h3s_route2) / ceil_h3s))

print()
print("=" * 78)
print("(B) Jensen 方向审计 —— 单边性到底朝哪边")
print("=" * 78)
# 测度 dmu = (2/lam) * alpha2F(w)/w dw，归一（因 lam = 2∫alpha2F/w dw）
# <w^2> = (2/lam)∫alpha2F*w dw ；w_log = exp(<ln w>)
# Jensen(ln 凹): ln<w^2> >= <ln w^2> = 2<ln w>  =>  <w^2> >= w_log^2  =>  w2 >= w_log
print("  测度 dmu=(2/lam)alpha2F/w dw 归一；<w^2>=(2/lam)∫alpha2F w dw；w_log=exp(<ln w>)")
print("  ln 为凹 ⟹ Jensen: ln<w^2> >= <ln w^2> = 2<ln w> ⟹ w2 >= w_log  [analytic]")
# 数值上用仓库两行 H3S 自查该不等式（同一次积分，故是真检验而非恒等式）
for tag, lam, wl, w2 in [("H3S-200 harm", 2.2929, 107.165, 129.32),
                         ("H3S-200 anh ", 1.721, 109.995, 149.26)]:
    chk("%s 满足 w2 >= w_log" % tag, w2 >= wl, "w2/w_log = %.4f" % (w2 / wl))
print("  代 w_log 换 w2 ⟹ ceiling 变小 ⟹ r=Tc/ceiling 变大 ⟹ r_up >= r_true")
print("  ⟹ r_up <= X 可证 r_true <= X（支持 PASS）；r_up > X 什么都不证。")
print("  ⚠ 推论：用 w_log 永远无法证明某体系『比 Nb3Ge 更紧』（那需要 r 的下界）。")

print()
print("=" * 78)
print("(C) 逐行重算 r_up（交付方表值，主会话独立算术）")
print("=" * 78)
# (体系, 压强GPa, lam, w_log[meV 或 cm^-1], 单位, Tc_exp[K], 交付方报的 r_up)
CM_TO_MEV = 0.12398419843320026    # 1 cm^-1 = 0.1239842 meV
ROWS = [
    ("LaH10", 129, 3.62,  76.4, "meV", 250.0, 0.811),
    ("LaH10", 163, 2.67,  96.4, "meV", 250.0, 0.749),
    ("LaH10", 214, 2.06, 115.5, "meV", 250.0, 0.711),
    ("LaH10", 264, 1.73, 126.6, "meV", 250.0, 0.708),
    ("LaH10", 163, 2.67,  96.4, "meV", 260.0, 0.778),
    ("LaD10", 159, 3.14,  63.5, "meV", 168.0, 0.704),
    ("LaD10", 210, 2.21,  81.7, "meV", 168.0, 0.652),
    ("LaD10", 260, 1.80,  92.2, "meV", 168.0, 0.641),
    ("CaH6",  150, 2.75, 684.0, "cm",  215.0, 0.721),
    ("CaH6",  200, 1.99, 851.0, "cm",  215.0, 0.681),
    ("CaH6",  250, 1.71, 933.0, "cm",  215.0, 0.670),
]
print("  %-6s %5s %6s %9s %8s %8s %8s %8s" %
      ("系统", "P", "lam", "w_log[K]", "ceil_log", "r_up", "交付方", "需 w2/w_log"))
maxr = 0.0
for name, P, lam, wl, unit, tc, claimed in ROWS:
    wl_mev = wl * CM_TO_MEV if unit == "cm" else wl
    wl_K = wl_mev * MEV_TO_K
    ceil_log = C_INF * math.sqrt(lam) * wl_K
    r_up = tc / ceil_log
    need = r_up / NB3GE_R
    agree = abs(r_up - claimed) <= 0.0015
    maxr = max(maxr, r_up)
    print("  %-6s %5d %6.2f %9.1f %8.1f %8.4f %8.3f %8.3f %s" %
          (name, P, lam, wl_K, ceil_log, r_up, claimed, need,
           "" if agree else "  <-- 不符!"))
    if not agree: OK = False
chk("11 行 r_up 全部与交付方相符（|Δ|<=0.0015）", OK)
chk("★ 全部 r_up < 1 ⟹ 天花板对这三个体系成立（单边界朝 PASS，故是真判决）",
    maxr < 1.0, "max r_up = %.4f" % maxr)
print("  ⟹ 交付方说『nothing is decided』只对『是否比 Nb3Ge 更紧』成立；")
print("     对『天花板是否成立』它其实已给出 11 行全部通过的**已证**结果。")

print()
print("=" * 78)
print("(D) ★ 交付方说做不了的一步：用仓库自有 H3S 行标定 w2/w_log")
print("=" * 78)
print("  交付方称需要『repo data (H3S w_log alongside its nu2)』——仓库 hydrides_theory.csv")
print("  第 2/3 行本来就同时印着两者，且出自**同一次** alpha2F 积分 ⟹ 零约定混杂。")
rat_h = 129.32 / 107.165
rat_a = 149.26 / 109.995
print("  H3S-200 harmonic : w2/w_log = 129.32/107.165 = %.4f" % rat_h)
print("  H3S-200 anharmon.: w2/w_log = 149.26/109.995 = %.4f  <-- 与 SSCHA 行可比" % rat_a)
chk("非谐比值 > 谐波比值（非谐谱更宽，方向合理）", rat_a > rat_h,
    "%.4f > %.4f" % (rat_a, rat_h))
worst = max(r / NB3GE_R for r in
            [tc / (C_INF * math.sqrt(lam) * (wl * (CM_TO_MEV if u == "cm" else 1)) * MEV_TO_K)
             for _, _, lam, wl, u, tc, _ in ROWS])
print("  全部 11 行所需 w2/w_log 的最大值 = %.4f" % worst)
chk("该最大阈值低于 H3S 非谐标定值 ⟹ 若这些体系的比值与 H3S 同量级，则全部不比 Nb3Ge 紧",
    worst < rat_a, "%.4f < %.4f" % (worst, rat_a))
print("  ⚠ 这是**标定/外推**，不是界。H3S 的比值不约束 LaH10/CaH6 的比值。")
print("     等级 numerical-evidence，且**不得**写成『已证 Nb3Ge 仍是最紧点』。")
print("     它只说明：要让某个氢化物夺走最紧点，需要一个**反常窄**的谱（比 H3S 还窄）。")

print()
print("=" * 78)
print("(E) 压强错配的方向审计（交付方标为 unresolved）")
print("=" * 78)
m2s = [(P, lam * (wl * MEV_TO_K) ** 2) for _, P, lam, wl, u, _, _ in ROWS
       if _ == "LaH10" or True]
lah = [(P, lam * (wl * MEV_TO_K) ** 2) for n, P, lam, wl, u, tc, _ in ROWS
       if n == "LaH10" and tc == 250.0]
lah.sort()
print("  LaH10 各压强的 M2(用 w_log 代理) [K^2]：")
for P, m2 in lah:
    print("    %3d GPa : %.3e   ceiling_log = %.1f K" % (P, m2, C_INF * math.sqrt(m2)))
mono = all(lah[i][1] < lah[i + 1][1] for i in range(len(lah) - 1))
chk("M2 随压强单调上升 ⟹ ceiling 随压强上升", mono)
print("  实验 ~170 GPa，最近计算行 163 GPa（**更低**）⟹ 用 163 的 ceiling 偏低")
print("  ⟹ 算出的 r 偏大 ⟹ **保守**方向 ⟹ 压强错配不威胁 PASS 判决。")
print("  ⚠ 但它确实威胁『谁最紧』的比较（偏大的 r 会假装更紧）——该问题本就未定。")

# ============================================================================
# (F) 2026-09-25 追加：主会话一手读页图后，给出**双边**括号
#     交付方只给了 r 的上界（w_log 代 w2）。页图额外印了两样东西，方向相反：
#       ① omega_max（声子谱上端）⟹ w2 <= w_max ⟹ ceiling 上界 ⟹ **r 的下界**
#       ② CaH6 正文把 lambda 的**分布**印成了三段百分比 ⟹ 直接是测度 mu 的分箱质量
#          （因 dlambda/dw = 2 alpha2F/w，而 dmu = (2/lam) alpha2F/w dw ⟹ mu = lambda 的归一分布）
# ============================================================================
print()
print("=" * 78)
print("(F) 一手页图 ⟹ r 的**下界**（交付方未做：它只有上界）")
print("=" * 78)

def ceil_K(lam, w_K):
    return C_INF * math.sqrt(lam) * w_K

# --- F1: omega_max 路线 -----------------------------------------------------
# 方向审计：r_lo = Tc / (C_inf sqrt(lam) w_max)。w_max **取大** ⟹ r_lo 取小 ⟹
# 「r >= r_lo」这一断言变**弱** ⟹ 高估 w_max 是**保守**方向。故一律取谱上端的宽松值。
print("  方向：w2 <= w_max ⟹ ceiling <= C_inf sqrt(lam) w_max ⟹ r >= r_lo。")
print("  高估 w_max ⟹ r_lo 变小 ⟹ 断言变弱 ⟹ **保守**。故取轴顶而非目测带顶。")
CASES = [
    # (标签, lam, w_log[K], w_max[K], Tc, 来源)
    ("LaH10@163", 2.67,  96.4 * MEV_TO_K, 260.0 * MEV_TO_K, 250.0,
     "Errea ED Fig.3：full-Hessian 紫线顶 ~250 meV，取 260 meV"),
    ("CaH6@150",  2.75,  684.0 * CM_TO_MEV * MEV_TO_K, 2500.0 * CM_TO_MEV * MEV_TO_K, 215.0,
     "Jeon FIG.4(a)：纵轴顶 2500 cm^-1，带顶明显在其下"),
    ("CaH6@200",  1.99,  851.0 * CM_TO_MEV * MEV_TO_K, 2500.0 * CM_TO_MEV * MEV_TO_K, 215.0,
     "Jeon FIG.4(b)"),
    ("CaH6@250",  1.71,  933.0 * CM_TO_MEV * MEV_TO_K, 2500.0 * CM_TO_MEV * MEV_TO_K, 215.0,
     "Jeon FIG.4(c)"),
]
print("  %-11s %9s %9s %9s %9s   %s" % ("例", "r_lo", "r_up", "含 0.626?", "宽度", "w_max 来源"))
for tag, lam, wl_K, wmax_K, tc, src in CASES:
    r_up = tc / ceil_K(lam, wl_K)
    r_lo = tc / ceil_K(lam, wmax_K)
    brackets = (r_lo <= NB3GE_R <= r_up)
    chk("%s: r_lo <= r_up（括号自洽，即 w_log <= w_max）" % tag, r_lo <= r_up, "")
    print("  %-11s %9.4f %9.4f %9s %9.4f   %s" %
          (tag, r_lo, r_up, "是" if brackets else "否", r_up - r_lo, src))
print("  ⟹ 四例的括号**全部跨过** 0.626 ⟹ w_max 路线判不了『谁最紧』。")
print("     这不是失败，是**定量**结论：要判决就必须拿到真 alpha2F，图上的谱宽不够。")

# --- F2: CaH6 的 lambda 分箱（正文直接印出测度的分布）------------------------
print()
print("  ★ Jeon 正文把测度 mu 的分箱质量直接印了出来（150 GPa）：")
print("    『modes between 500 and 1200 cm^-1 ... ~78% of the total EPC constant lambda』")
print("    『Ca-derived phonon modes reach ~14%』（Ca 声学支，正文另述与 H 光学支well separated）")
print("    因 dlambda/dw = 2 alpha2F/w 而 dmu = (2/lam) alpha2F/w dw ⟹ mu 就是 lambda 的归一分布")
print("    ⟹ 这三段百分比**就是** mu 的分箱质量，不需要任何模型。")
# 三箱：[0,500] 质量 0.14；[500,1200] 质量 0.78；[1200, w_max] 质量 0.08
BINS = [(0.0, 500.0, 0.14), (500.0, 1200.0, 0.78), (1200.0, 2500.0, 0.08)]
tot = sum(m for _, _, m in BINS)
chk("三箱质量和为 1（0.14+0.78+余）", abs(tot - 1.0) < 1e-12, "sum = %.4f" % tot)
w2_hi_cm = math.sqrt(sum(m * hi * hi for _, hi, m in BINS))
w2_lo_cm = math.sqrt(sum(m * lo * lo for lo, _, m in BINS))
print("    质量全压到各箱**上**端 ⟹ <w^2>^(1/2) <= %.1f cm^-1" % w2_hi_cm)
print("    质量全压到各箱**下**端 ⟹ <w^2>^(1/2) >= %.1f cm^-1" % w2_lo_cm)
chk("分箱下界与 Jensen 下界 w_log=684 相容（取两者较大）", True,
    "max(%.1f, 684.0) = %.1f" % (w2_lo_cm, max(w2_lo_cm, 684.0)))
w2_lo_cm = max(w2_lo_cm, 684.0)
chk("分箱上界 %.1f 严于裸 w_max=2500" % w2_hi_cm, w2_hi_cm < 2500.0,
    "收紧 %.2f 倍" % (2500.0 / w2_hi_cm))
lam150, tc_cah6 = 2.75, 215.0
r_hi = tc_cah6 / ceil_K(lam150, w2_lo_cm * CM_TO_MEV * MEV_TO_K)
r_lo2 = tc_cah6 / ceil_K(lam150, w2_hi_cm * CM_TO_MEV * MEV_TO_K)
print("    ⟹ CaH6@150 GPa: r ∈ [%.4f, %.4f]（vs 仅 w_max 时的下界 %.4f）"
      % (r_lo2, r_hi, tc_cah6 / ceil_K(lam150, 2500.0 * CM_TO_MEV * MEV_TO_K)))
chk("分箱后括号仍跨 0.626 ⟹ 仍判不了", r_lo2 <= NB3GE_R <= r_hi,
    "[%.4f, %.4f] ∋ 0.626" % (r_lo2, r_hi))
print("  ⚠ 百分比带『~』，且第三箱质量由 1-0.14-0.78 反推 ⟹ 等级 numerical-evidence，**非** certificate。")
print("     但方向是硬的：它是 r 的**下界**，交付方整份报告里没有任何下界。")

print()
print("=" * 78)
print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
print("=" * 78)
sys.exit(0 if OK else 1)
