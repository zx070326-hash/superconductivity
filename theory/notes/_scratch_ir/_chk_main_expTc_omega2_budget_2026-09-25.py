#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""预算：「实验 Tc + 该材料 ω̄₂」这条路线能不能动 Nb₃Ge 的 0.6260？

消费点（必须先写下来，CLAUDE.md「预算优先」）：
  r_true = Tc_exp / (C_INF * sqrt(lam) * w2[K]) 与 Nb3Ge 的 0.6260 比大小。
  这个判词是**排名**判词（谁最紧），不是 PASS/FAIL 判词。

方向（单边界方向规则的应用）：
  预算要朝**容易翻转**的方向退化 —— Tc_exp 取上端、lam 取下端、w2 取下端。
  若连这样都翻不过去 ⟹ 预算 ∞ ⟹ 不派活。
  若翻得过去 ⟹ 预算有限 ⟹ 派活，但必须同时交付「需要多准」这个消费点。

本脚本不产生新的物理结论，只产生**预算**与**派活规格**。
"""
import math, sys

MEV_K = 1.602176634e-22 / 1.380649e-23        # K per meV
D_LAM = 0.005                                  # lam 印 2 位小数
C_INF = 0.1827262477
NB3GE = 0.6260

assert abs(MEV_K - 11.6045181) < 1e-6, MEV_K

FAIL = []
def chk(cond, msg):
    if not cond:
        FAIL.append(msg); print("  [FAIL]", msg)

# ---------------------------------------------------------------- 数据
# Errea 1907.11916 ED Table 1，**已剔除 @163**（Mc 复现 -3.30%，兄弟行 <=0.251%，离群 13.1x）
# r 区间来自 errea_ed_table1_from_zero_recheck_2026-09-25.md §3（输入侧印刷精度已传播）
# (P, lam, wlog_meV, r_lo, r_hi)
LAH = [(129, 3.62,  76.4, 1.3180, 1.3250),
       (214, 2.06, 115.5, 1.2275, 1.2405),
       (264, 1.73, 126.6, 1.2281, 1.2460)]
LAH163 = (163, 2.67, 96.4, 1.1272, 1.1332)     # 弃用行，仅作诊断对象
LAD = [(159, 3.14,  63.5, 1.2437, 1.2518),
       (210, 2.21,  81.7, 1.1858, 1.1974),
       (260, 1.80,  92.2, 1.1755, 1.1914)]

# 实验（theory/data/hydrides_experimental.csv，一手来源见该表 notes 列）
# LaH10: 170-188 GPa, Tc = 250 K (Drozdov Nature 569:528) / 260 K (Somayazulu PRL 122:027001)
EXP_P_LO, EXP_P_HI = 170.0, 188.0
EXP_TC = [("Drozdov2019", 250.0), ("Somayazulu2019", 260.0)]

def w2K(wlog_meV, r):
    return r * wlog_meV * MEV_K

def r_true(tc, lam, w2):
    return tc / (C_INF * math.sqrt(lam) * w2)

# ---------------------------------------------------------------- §0 自检
print("=" * 78)
print("§0 常数与单调性自检")
print("=" * 78)
chk(abs(MEV_K - 11.6045181) < 1e-6, "MEV_K")
# r_true 对 lam / w2 单调下降、对 tc 单调上升 => 角点扫描合法
chk(r_true(250, 2.0, 1400) > r_true(250, 2.5, 1400), "r_true 对 lam 递减")
chk(r_true(250, 2.0, 1400) > r_true(250, 2.0, 1500), "r_true 对 w2 递减")
chk(r_true(260, 2.0, 1400) > r_true(250, 2.0, 1400), "r_true 对 tc 递增")
print("  MEV_K = %.7f K/meV ; r_true 三方向单调性 OK" % MEV_K)
print("  => 取 (max Tc, min lam, min w2) 即给 r_true 的上包封（角点即极值）")

# ---------------------------------------------------------------- §A 包封预算
print()
print("=" * 78)
print("§A 包封预算：最松的可采纳组合能不能越 0.6260？")
print("=" * 78)
print("  窗口取 [129, 214] GPa（把实验的 170-188 夹在内，且两端都是**存活**行），")
print("  只用 lam(P) 递减 / w2(P) 递增 ⟹ lam >= lam(214), w2 >= w2(129)。")
print("  这是**无形状假设**的包封；它只回答「预算是不是 oo」。")

def lin(x, x0, y0, x1, y1):
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

w2_129 = w2K(LAH[0][2] - 0.05, LAH[0][3])     # wlog/r 均取下端 => w2 下界
w2_214 = w2K(LAH[1][2] - 0.05, LAH[1][3])
lam_win_min = LAH[1][1]                        # 2.06 @214
w2_win_min  = w2_129                           # @129 印刷下端
chk(LAH[0][1] > LAH[1][1] > LAH[2][1], "lam 随 P 递减")
chk(w2_129 < w2_214 < w2K(LAH[2][2], LAH[2][3]), "w2 随 P 递增")
print("  包封: lam_min = %.2f (@214), w2_min = %.1f K (@129)" % (lam_win_min, w2_win_min))
env = {}
for nm, tc in EXP_TC:
    v = r_true(tc, lam_win_min, w2_win_min)
    env[nm] = v
    print("    Tc_exp = %.0f K (%s) -> r_true <= %.4f   %s"
          % (tc, nm, v, "★越线" if v > NB3GE else "未越线"))
chk(max(env.values()) > NB3GE, "包封应当越线（否则预算 oo，本脚本结论要改写）")
print("  => 包封**越线** ⟹ 预算**不是 oo** ⟹ 这条线不能用「反正翻不了」打发掉。")

# ---------------------------------------------------------------- §B 压强匹配估计
print()
print("=" * 78)
print("§B 压强匹配估计：把 lam/w2 插到实验压强上（这才是物理上共存的组合）")
print("=" * 78)
print("  自由度诚实声明：LaH10 存活 3 行，线性插值用 2 参数 ⟹ 残差自由度 1。")
print("  lam(P) 与 w2(P) 各自用 @129/@214 两点定直线，@264 作**留一验证**。")

lam_lin = lambda P: lin(P, 129, LAH[0][1], 214, LAH[1][1])
w2_lin  = lambda P: lin(P, 129, w2_129, 214, w2_214)
pred_lam264, pred_w2264 = lam_lin(264), w2_lin(264)
true_lam264, true_w2264 = LAH[2][1], w2K(LAH[2][2], LAH[2][3])
print("  留一验证 @264（**外**插）: lam 预测 %.3f vs 印出 %.2f (%+.1f%%) ; "
      "w2 预测 %.1f vs 印出 %.1f (%+.1f%%)"
      % (pred_lam264, true_lam264, 100*(pred_lam264-true_lam264)/true_lam264,
         pred_w2264, true_w2264, 100*(pred_w2264-true_w2264)/true_w2264))

# ---- 方向审计（订正 2026-09-25：先前把**外插**的方向搬给了**内插**，反了）----
print()
print("  ★ 方向审计（本节判词的成败全在这里）")
sl_lam = [(LAH[1][1]-LAH[0][1])/(214-129), (LAH[2][1]-LAH[1][1])/(264-214)]
sl_w2  = [(w2_214-w2_129)/(214-129), (true_w2264-w2_214)/(264-214)]
print("    lam 分段斜率 %+.5f -> %+.5f /GPa  ⟹ 递减且**凸**" % tuple(sl_lam))
print("    w2  分段斜率 %+.5f -> %+.5f /GPa  ⟹ 递增且**凹**" % tuple(sl_w2))
chk(sl_lam[1] > sl_lam[0], "lam 凸（斜率上升）")
chk(sl_w2[1] < sl_w2[0], "w2 凹（斜率下降）")
print("    凸函数的弦在区间**内部**位于函数**上**方 ⟹ 内插**高估** lam ⟹ 高估 ceiling ⟹ **低估** r_true")
print("    凹函数的弦在区间**内部**位于函数**下**方 ⟹ 内插**低估** w2  ⟹ 低估 ceiling ⟹ **高估** r_true")
print("    ⟹ 两个偏置**方向相反**，净方向**不定**。")
print("    ⟹ §B 的数是 numerical-evidence 级**点估计**，")
print("       它**不支持**「不越线」这个判词（那需要 r_true 的**上**界）。")
print("    （外插 @264 的方向与此相反——凸函数的弦在区间**外**位于函数**下**方，")
print("      这正是留一验证低估 lam 34% 的原因；两者不是同一个方向命题。）")

print()
print("  %-10s %6s %8s %10s %10s %10s" % ("Tc_exp", "P", "lam", "w2(K)", "r_true", "判"))
matched = []
for P in (170.0, 180.0, 188.0):
    lm, w2 = lam_lin(P), w2_lin(P)
    for nm, tc in EXP_TC:
        v = r_true(tc, lm, w2)
        matched.append((P, tc, v))
        print("  %-10.0f %6.0f %8.3f %10.1f %10.4f %10s"
              % (tc, P, lm, w2, v, "★越线" if v > NB3GE else "未越线"))
mx = max(v for _, _, v in matched)
print("  压强匹配点估计最大 r_true = %.4f" % mx)
print("  ⚠ 上一段已证此数**方向不定** ⟹ 它不是判词，只是量级。")

# ---- §B2 真正的上界：只用单调性 + 形状（凸/凹），不用线性模型 ----
print()
print("  §B2 r_true 的**上**界（这才能下判词）")
print("     (a) 只用单调性（最弱假设，无形状）：")
lam_lo_mono = LAH[1][1]          # lam(P<=214) >= lam(214)
w2_lo_mono  = w2_129             # w2(P>=129) >= w2(129)
v_mono = r_true(260.0, lam_lo_mono, w2_lo_mono)
print("        lam >= %.2f, w2 >= %.1f K  ->  r_true <= %.4f  %s"
      % (lam_lo_mono, w2_lo_mono, v_mono, "★仍越线" if v_mono > NB3GE else "未越线"))
print("     (b) 追加形状（lam 凸 / w2 凹，由存活三行实测，见上）：")
print("        lam 凸 ⟹ 在 [214,264] 割线**外**延处 lam(P) >= 割线")
print("        w2  凹 ⟹ 在 [129,214] 弦**内**部处   w2(P) >= 弦")
sec_lam = lambda P: lin(P, 214, LAH[1][1] - D_LAM, 264, LAH[2][1] - D_LAM)
print("        （lam 印 2 位小数 ⟹ 两端各减 %.3f，使割线成为**下**界）" % D_LAM)
shape = []
for P in (170.0, 180.0, 188.0):
    lo_l, lo_w = sec_lam(P), w2_lin(P)
    for nm, tc in EXP_TC:
        v = r_true(tc, lo_l, lo_w)
        shape.append(v)
        print("        Tc=%.0f @%.0f: lam >= %.4f, w2 >= %.1f  ->  r_true <= %.4f  %s"
              % (tc, P, lo_l, lo_w, v, "★仍越线" if v > NB3GE else "未越线"))
ub_shape = max(shape)
print("     ⟹ 形状上界 max r_true = %.4f  %s 0.6260" % (ub_shape, ">" if ub_shape > NB3GE else "<"))
chk(ub_shape < v_mono, "形状上界应当严于单调性上界")
chk(all(sec_lam(P) > LAH[1][1] - D_LAM for P in (170, 180, 188)), "割线外延确实高于 lam(214)")
print("     ⚠ 形状假设只在**三行**上实测过（残差自由度 1）⟹ 本条是")
print("        conditional bound，不是 certificate；但它**不依赖**被弃用的 @163。")

# ---------------------------------------------------------------- §C 消费点
print()
print("=" * 78)
print("§C 消费点：要判死这件事，w2 需要准到多少？")
print("=" * 78)
print("  判词 r_true > 0.6260 等价于  sqrt(lam)*w2 < Tc_exp/(0.6260*C_INF)")
for nm, tc in EXP_TC:
    thr = tc / (NB3GE * C_INF)
    print("  Tc_exp = %.0f K: 需要 sqrt(lam)*w2 < %.1f K" % (tc, thr))
    for P in (170.0, 180.0):
        lm = lam_lin(P)
        w2_need = thr / math.sqrt(lm)
        w2_have = w2_lin(P)
        print("     @%.0f GPa (lam=%.3f): w2 临界 = %.1f K = %.1f meV ; "
              "插值给 %.1f K ⟹ 需比插值**低** %.1f%% 才越线"
              % (P, lm, w2_need, w2_need / MEV_K, w2_have,
                 100 * (w2_have - w2_need) / w2_have))

# ---------------------------------------------------------------- §D 第三条 @163 诊断
print()
print("=" * 78)
print("§D 顺带：@163 的**第三条**独立诊断——印出的 lam 合不合 P 趋势？")
print("=" * 78)
print("  已有两条诊断：(1) Tc^Mc 零参数复现 -3.30%；(2) r=1.1302 低于 Quan 下端 1.2500。")
print("  本条用**存活三行的 lam(P) 趋势**外推 @163，与印出的 2.67 比。")
# 用 @129/@214/@264 三点拟合 lam = a + b/P（两参数，残差自由度 1）
import itertools
def fit_ab(pts):
    n = len(pts)
    sx = sum(1/P for P, L in pts); sy = sum(L for P, L in pts)
    sxx = sum(1/P**2 for P, L in pts); sxy = sum(L/P for P, L in pts)
    d = n*sxx - sx*sx
    b = (n*sxy - sx*sy)/d; a = (sy - b*sx)/n
    return a, b
pts = [(P, L) for P, L, *_ in LAH]
a, b = fit_ab(pts)
resid = [L - (a + b/P) for P, L in pts]
print("  拟合 lam = a + b/P: a=%.4f b=%.1f ; 三行残差 = %s"
      % (a, b, ", ".join("%+.4f" % x for x in resid)))
lam163_fit = a + b/163
print("  外推 @163: lam_fit = %.4f  vs 印出 %.2f  (印出偏低 %.1f%%)"
      % (lam163_fit, LAH163[1], 100*(lam163_fit - LAH163[1])/lam163_fit))
print("  对照：§B 单改 lam 复现 Tc^Mc 需要 lam = 2.8630")
print("  ⟹ 趋势外推 %.4f 与「复现所需 %.4f」同侧且接近；印出的 2.67 在两者之外"
      % (lam163_fit, 2.8630))
out163 = abs(LAH163[1] - lam163_fit)
scat = max(abs(x) for x in resid)
print("  离群比 = |2.67 - %.4f| / max|存活残差| = %.4f / %.4f = **%.1fx**"
      % (lam163_fit, out163, scat, out163/scat))
print("  ⚠ 对比：Tc^Mc 复现门给的离群比是 **13.1x**。本条只有 %.1fx" % (out163/scat))
print("     ⟹ 本条是**弱**诊断，不单独承重、不改任何已下判词；它的作用只是**同向**。")
print("  ⚠ 限度：三点拟两参数，残差自由度仅 1；lam(P) 未必是 a+b/P。")
chk(abs(max(resid)) < 0.15, "a+b/P 在存活三行上残差应当小（否则本诊断无效）")

# ---------------------------------------------------------------- §E LaD10 对照
print()
print("=" * 78)
print("§E LaD10 对照（同位素侧，Tc_exp 溯源最弱）")
print("=" * 78)
print("  csv notes: LaD10 Tc=180K@~150GPa 是**转引自 Durajski 理论论文**，一手实验源未确认")
print("  ⟹ 本节只算、**不进任何判词**。")
lamD = lambda P: lin(P, 159, LAD[0][1], 210, LAD[1][1])
w2D_159, w2D_210 = w2K(LAD[0][2], LAD[0][3]), w2K(LAD[1][2], LAD[1][3])
w2D = lambda P: lin(P, 159, w2D_159, 210, w2D_210)
for P in (150.0, 159.0, 170.0):
    v = r_true(180.0, lamD(P), w2D(P))
    print("  Tc_exp=180K @%.0f GPa: lam=%.3f w2=%.1fK -> r_true=%.4f %s"
          % (P, lamD(P), w2D(P), v, "★越线" if v > NB3GE else "未越线"))
print("  ⚠ P=150 在 [159,210] 之外，是外插。")

# ---------------------------------------------------------------- §F 投毒
print()
print("=" * 78)
print("§F 投毒负对照（检测器必须真的走过判决路径）")
print("=" * 78)
base_env = max(env.values())
base_matched = mx
base_lam163 = lam163_fit
poison_fail = []
def pchk(name, cond, detail):
    print("  %-34s %s  %s" % (name, "OK" if cond else "MISS", detail))
    if not cond: poison_fail.append(name)

# P1: w2 单位当成 meV 忘乘 MEV_K => r_true 爆炸
v = 250.0 / (C_INF * math.sqrt(lam_lin(170)) * (w2_lin(170)/MEV_K))
pchk("P1 w2 漏乘 K/meV", v > 5, "r_true = %.2f（应荒谬）" % v)
# P2: 把弃用的 @163 放回来（它 w2 最低 => 最偏向越线）
w2_163_printed = w2K(LAH163[2], LAH163[3])
v = r_true(260.0, LAH163[1], w2_163_printed)
pchk("P2 复用弃用行 @163（正对照）", v > NB3GE,
     "r_true = %.4f 越线 ⟹ 弃用它确实是判决的**承重**一步" % v)
# P3: Tc_exp 换成 Errea 算出的 ME_ani（同压强）=> 应回落，且这是信息类混用
v = r_true(237.9, LAH[1][1], w2_214)
pchk("P3 用计算 Tc 冒充实验 Tc", v < base_matched,
     "r_true = %.4f < 压强匹配上端 %.4f" % (v, base_matched))
# P4: lam 取窗口**上**端（方向审反）=> 包封应失效
v = r_true(260.0, LAH[0][1], w2_win_min)          # lam = 3.62 窗口上端
pchk("P4 方向审反（lam 取上端）", v < base_env,
     "r_true = %.4f < 包封 %.4f ⟹ 方向确实在起作用" % (v, base_env))
# P5: @163 的 lam 诊断——它必须是四行里残差最大的那行（**秩**统计量，不拍阈值）
pts5 = pts + [(163, LAH163[1])]
a5, b5 = fit_ab(pts5)
r5 = {P: L - (a5 + b5/P) for P, L in pts5}
worst = min(r5, key=lambda P: -abs(r5[P]))
pchk("P5 @163 是 lam(P) 拟合的最大残差行", worst == 163,
     "各行残差 %s ⟹ 最大在 @%d"
     % (", ".join("@%d:%+.4f" % (P, r5[P]) for P in sorted(r5)), worst))
# P6: 内插方向不定的**正对照**——若强行把点估计当上界，会把「判不了」误报成「不越线」
pchk("P6 点估计冒充上界会翻错判词", (mx < NB3GE) and (ub_shape > NB3GE),
     "点估计 %.4f<0.6260 但真上界 %.4f>0.6260 ⟹ 误报「已判死」" % (mx, ub_shape))

print()
print("  投毒 %d/6 通过" % (6 - len(poison_fail)))
chk(not poison_fail, "投毒未全通过: %s" % poison_fail)

# ---------------------------------------------------------------- 结论
print()
print("=" * 78)
print("预算结论")
print("=" * 78)
print("  无形状包封  r_true <= %.4f  > 0.6260" % base_env)
print("  形状上界    r_true <= %.4f  > 0.6260  ⟹ **判不了**，不是「翻不了」" % ub_shape)
print("  点估计      r_true ~= %.4f  < 0.6260  ⟹ **方向不定，不构成判词**" % mx)
print()
print("  ⟹ 预算**有限**（不是 oo）⟹ 这条线该派活。")
print("  ⟹ 消费点：需要 170-188 GPa 处 LaH10 的 (lam, w2)。")
print("     判死需要的精度：w2 的上界须落在临界值 %.1f-%.1f K 之外（见 §C）"
      % (min(t/(NB3GE*C_INF)/math.sqrt(lam_lin(P)) for _, t in EXP_TC for P in (170.,180.)),
         max(t/(NB3GE*C_INF)/math.sqrt(lam_lin(P)) for _, t in EXP_TC for P in (170.,180.))))
print("     即相对精度需优于 ~4%（Tc=260 侧），这是**可达**的（Quan 印出值精度远好于此）。")

print()
if FAIL:
    print("SOME CHECKS FAILED:"); [print("  -", m) for m in FAIL]; sys.exit(1)
print("ALL CHECKS PASSED")
