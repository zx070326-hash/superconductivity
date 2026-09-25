#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主会话从零复核：Savrasov-Savrasov-Andersen, PRL 72, 372 (1994) TABLE I (p.374)
                以及 teammate `eta-second-method` §5 的三处缺陷。

本脚本**不 import** 任何既有脚本。两台仪器都跑过：
  文本层 theory/notes/_scratch_ir/_savrasov_prl72_utf8.txt (pdftotext -enc UTF-8 -layout)
  页 图 library/pdf/eta_hunt_2026-09-25/savrasov_prl72_tableI-3.png (pdftoppm -r 250, 第 3 页)
文本层是 OCR 产物，把 λ 印成 `A`/`p`（行标签 `pcalc`/`pexp`/`pctarlc` 实为
λ^calc / λ^exp / λ_tr^calc）；**数值与列归属以页图为准**，两者在数值上逐字一致。

检验三件事（互相独立，分开判）：

  (D1) 列名审计：TABLE I 的列**只有 Al / Nb / Mo**，没有 V 行。
       ⟹ teammate 报的 `λ_Savr(V) = 1.19` 不出自 PRL 72,372。

       ⚠⚠ 本脚本 v1 在这里下了一个**假判词**，已订正（2026-09-25 同日，见下）：
       v1 写「1.19 恰是 Radevych TABLE I 里 V 的 λ ⟹ 疑似串源、该行作废」。
       **这是错的。** 1.19 确实同时出现在 Radevych TABLE I 的 λ 列，但那只是**数值巧合**，
       不是出处证据 —— 我拿「两个数相等」当了「同一个来源」，
       正是 CLAUDE.md W1 热区那条：**两个数一致时先问「它们是同一个约定下的量吗」，
       再问「谁错了」**，我把它用反了；且我没有先 grep 仓库里还有没有别的 1.19 出处。
       **真出处**：Savrasov & Savrasov, **PRB 54, 16487 (1996)** TABLE II（期刊页 16 494），
       `λ_calc` 行 = `Al 0.44 | Pb 1.68 | V 1.19 | Nb 1.26 | Ta 0.86 | Mo 0.42 | Cu 0.14 | Pd 0.35`。
       我本轮亲自读了页图 `library/pdf/eta_hunt_2026-09-25/savr_rmta_p8-08.png` 确认。
       ⟹ V 那行**不作废**，它的真缺陷是**两篇论文揉进一张表却没有逐行标出处**
       （16470 与 16487 是两篇不同的论文），以及 (D3) 的混杂同样适用于它。

  (D2) 方向审计（单边界方向规则）：ceiling ∝ √η ⟹ r = Tc/ceiling ∝ 1/√η。
       η 变大 ⟹ r 变**小**。teammate 写的是 r_new = r × √ratio（乘），方向反了。
       正确是 r_new = r / √ratio。
       ⚠ 但它的算术**恰好偏保守**（给出比真值更大的 r），
       ⟹ 它的 PASS 判决**存活**，垮掉的是理由不是结论。这两件事要分开写。

  (D3) 混杂审计：λ = η/(M⟨ω²⟩)。Pap1977 的 λ 用 ⟨ω²⟩ = Θ_D²/2（德拜代餐），
       Savrasov 的 λ 用**算出来的声子谱**。
       ⟹ λ_Savr/λ_Pap **不是** η 的方法偏差，它混了 ⟨ω²⟩ 的差。
       本脚本给出这个混杂的显式分解，并算出「要让它成为 η 比值」所需的前提。

严格性：numerical-evidence（两端都是他人发表的值；无 analytic 新命题、无证书）。
"""
import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

# ---- Savrasov PRL 72, 372 (1994) TABLE I，逐字（页图 savrasov_prl72_tableI-3.png）----
# 列头只有三个：Al | Nb | Mo
SAVRASOV_TABLE_I_COLUMNS = ["Al", "Nb", "Mo"]
SAVRASOV_TABLE_I = {
    #        lam_calc, lam_exp, lam_tr_calc
    "Al": (0.44, 0.42, 0.36),
    "Nb": (1.26, 1.33, 1.17),
    "Mo": (0.42, 0.44, 0.35),
}
# 正文逐字 (p.374 左栏)：`In order to check previous conclusions [7,22] about the
# inapplicability of the RMTA for sp metals, we also performed such a calculation
# and indeed found lambda = 0.14.`   ⟹ 这是**唯一**一个干净的 RMTA-vs-线性响应数对。
SAVRASOV_AL_RMTA_LAMBDA = 0.14
# 过渡金属侧只有散文，没有数：
#   Nb (p.374 左栏): `We obtain similar results with the RMTA, in full agreement
#                     with earlier calculations [17,24]`
#   Mo (p.374 右栏): `Our linear-response calculations of alpha^2 F(omega) are found
#                     to be close to our RMTA calculations and to earlier ones [30]`
SAVRASOV_TM_RMTA_NUMBERS = {}          # 空集合：这是本轮的核心否定发现

# ---- teammate `eta-second-method` §5 报的表 ----
AGENT_LAMBDA_SAVR = {"V": 1.19, "Nb": 1.26, "Mo": 0.42}
AGENT_LAMBDA_PAP = {"V": 1.055, "Nb": 1.222, "Mo": 0.336}
AGENT_RNEW = {"Nb": 0.636, "Mo": 0.700}     # 它写的 r_new
# Radevych et al. arXiv:2508.21281 TABLE I 的 lambda 列（见同目录 radevych 复核脚本）
RADEVYCH_LAMBDA = {"V": 1.19, "Cu": 0.09, "Nb": 1.45, "Mo": 0.35,
                   "Rh": 0.24, "Pd": 0.31, "Ag": 0.06}
# Savrasov & Savrasov, PRB 54, 16487 (1996) TABLE II (期刊页 16 494) 的 lambda_calc 行，
# 逐字读自页图 library/pdf/eta_hunt_2026-09-25/savr_rmta_p8-08.png（主会话本轮亲自读）。
# **这才是 teammate 那个 V=1.19 的真出处** —— 见 (D1) 的订正说明。
# ⚠ 16470 与 16487 是**两篇不同的论文**：16470 是线性响应 MTO 方法论文（只有 Si + NbC 声子频率），
#   16487 才是「金属的电声耦合与相关物理性质」那篇（有 TABLE II）。
SAVRASOV_PRB54_16487_TABLE_II_LAMCALC = {
    "Al": 0.44, "Pb": 1.68, "V": 1.19, "Nb": 1.26,
    "Ta": 0.86, "Mo": 0.42, "Cu": 0.14, "Pd": 0.35,
}

NB3GE_R = 0.6260          # 仓库当前最紧点（Klein 1979 A15, M2_min 档）

fails = []
sep = "=" * 78


def hdr(t):
    print()
    print(sep)
    print(t)
    print(sep)


# ====================== (D1) 列名审计 ======================
hdr("(D1) TABLE I 的列名审计 —— teammate 的 V 行从哪来？（v1 判词已订正）")
print(f"PRL 72,372 页图读出的列头：{SAVRASOV_TABLE_I_COLUMNS}")
for el in AGENT_LAMBDA_SAVR:
    val = AGENT_LAMBDA_SAVR[el]
    if el in SAVRASOV_TABLE_I:
        ok = abs(val - SAVRASOV_TABLE_I[el][0]) <= 1e-9
        print(f"  {el:<3} agent 报 {val:.2f} ; PRL72,372 λ^calc = "
              f"{SAVRASOV_TABLE_I[el][0]:.2f}  ⟹ {'一致' if ok else '不一致'}")
        if not ok:
            fails.append(f"(D1) {el}: agent {val} 与 PRL72,372 的 {SAVRASOV_TABLE_I[el][0]} 不符")
        continue

    # 这一行不在 PRL 72,372 里 —— 穷举**所有已知出处**，不要只看一个就下判词
    print(f"  {el:<3} agent 报 {val:.2f} ; **PRL72,372 无 {el} 行** ⟹ 逐个出处排查：")
    srcs = []
    if el in SAVRASOV_PRB54_16487_TABLE_II_LAMCALC and \
            abs(SAVRASOV_PRB54_16487_TABLE_II_LAMCALC[el] - val) <= 1e-9:
        srcs.append("PRB 54,16487 TABLE II λ_calc（页图 savr_rmta_p8-08.png，期刊页 16 494）")
    if el in RADEVYCH_LAMBDA and abs(RADEVYCH_LAMBDA[el] - val) <= 1e-9:
        srcs.append("Radevych arXiv:2508.21281 TABLE I 的 λ 列（半唯象，基于 RMTA η + Θ_D）")
    for s in srcs:
        print(f"        · 命中：{s}")
    if not srcs:
        fails.append(f"(D1) {el}={val} 找不到任何已知出处 ⟹ 真·无出处，须追")
    print(f"        ⟹ 出处**不唯一**（{len(srcs)} 个）。")
    print( "           ⚠ 『两个数相等』**不是**出处证据。本脚本 v1 只查了 Radevych 一处，")
    print( "             据此判『串源、该行作废』——**那个判词是假的**。")
    print( "             真出处是 PRB 54,16487（与 16470、与 PRL72,372 都是不同的论文）。")
print()
print("⟹ 订正后的判词：V 行**不作废**（1.19 是 PRB 54,16487 TABLE II 的已发表值）。")
print("   它的真缺陷是 (a) 多篇论文揉进一张表却**没有逐行标出处**，")
print("   以及 (b) (D3) 的 ⟨ω²⟩ 混杂对它同样适用 ⟹ 不承重。")
if "V" in SAVRASOV_TABLE_I:
    fails.append("(D1) 内部矛盾：PRL72,372 不该有 V 行")
if abs(SAVRASOV_PRB54_16487_TABLE_II_LAMCALC["V"] - AGENT_LAMBDA_SAVR["V"]) > 1e-9:
    fails.append("(D1) V 的真出处核对失败")

# ====================== (D2) 方向审计 ======================
hdr("(D2) 方向审计 —— r ∝ 1/√η，agent 乘了应该除的那个因子")
print(f"{'元素':<4} {'λ_Savr':>8} {'λ_Pap':>8} {'比值':>7} "
      f"{'agent r_new':>12} {'正确 r_new':>11} {'差':>8}  agent 方向")
for el in ("Nb", "Mo"):
    lam_s = SAVRASOV_TABLE_I[el][0]      # 用我自己读的，不用 agent 的
    lam_p = AGENT_LAMBDA_PAP[el]
    ratio = lam_s / lam_p
    r_agent = AGENT_RNEW[el]
    r_correct = NB3GE_R / math.sqrt(ratio)
    r_agent_recomputed = NB3GE_R * math.sqrt(ratio)
    # 先确认我复现了 agent 的算术（否则是我读错了它写的什么）
    if abs(r_agent_recomputed - r_agent) > 0.006:
        fails.append(f"(D2) 未能复现 agent 对 {el} 的算术："
                     f"重算 {r_agent_recomputed:.3f} vs 它写的 {r_agent}")
    conservative = r_agent > r_correct
    print(f"{el:<4} {lam_s:>8.2f} {lam_p:>8.3f} {ratio:>7.4f} "
          f"{r_agent:>12.3f} {r_correct:>11.4f} {r_agent-r_correct:>+8.4f}  "
          f"{'偏保守(高估 r)' if conservative else '⚠ 偏激进(低估 r)'}")
    if not conservative:
        fails.append(f"(D2) {el}: agent 的方向错**且**结果偏激进 ⟹ PASS 判决不存活")
print()
print("单边界方向规则的应用：判据是 r ≤ 1。**高估** r 只支持 PASS，永不支持 FAIL。")
print("⟹ agent 的算术方向反了，但反的方向恰好使 r 被高估 ⟹ 它的 PASS 结论**存活**，")
print("   垮掉的是**理由**。这两件事必须分两段写（不能说『它错了所以结论也错』）。")

# ====================== (D3) 混杂审计 ======================
hdr("(D3) 混杂审计 —— λ 比值不是 η 比值")
print("λ = η / (M⟨ω²⟩)  ⟹  λ_Savr/λ_Pap = (η_Savr/η_Pap) · (⟨ω²⟩_Pap/⟨ω²⟩_Savr)")
print("Pap1977 取 ⟨ω²⟩ = Θ_D²/2（德拜代餐）；Savrasov 用算出的声子谱 ⟹ 两者不是同一个量。")
print()
print(f"{'元素':<4} {'λ 比值':>8}  {'要让 λ 比值 = η 比值，需':>28}")
for el in ("Nb", "Mo"):
    ratio = SAVRASOV_TABLE_I[el][0] / AGENT_LAMBDA_PAP[el]
    print(f"{el:<4} {ratio:>8.4f}  {'⟨ω²⟩_Pap/⟨ω²⟩_Savr = 1（未验证）':>28}")
print()
print("⟹ 这是一个**两因子**的乘积，我们只观测到乘积。")
print("   把它整个归给 η ⟹ 误差是**双边**的（⟨ω²⟩ 的比值可大于也可小于 1）")
print("   ⟹ 单边界方向规则的豁免**不适用**（同 sandwich-kernel 事故形状）。")
print("   判词只能是：λ 比值**不构成** RMTA 的 η 方法偏差证据。")

# ====================== (D4) 真正干净的那一个数 ======================
hdr("(D4) 本文唯一干净的 RMTA 方法偏差数：Al")
al_lr = SAVRASOV_TABLE_I["Al"][0]
al_rmta = SAVRASOV_AL_RMTA_LAMBDA
print(f"Al：λ(线性响应，同一篇同一套势) = {al_lr:.2f}；λ(RMTA，同一篇) = {al_rmta:.2f}")
print(f"    比值 λ_RMTA/λ_LR = {al_rmta/al_lr:.4f}  ⟹ RMTA 低估 {al_lr/al_rmta:.2f}×")
print("    ⚠ 这个数干净是因为**同一篇、同一套势、同一个 ⟨ω²⟩ 约定** ⟹ (D3) 的混杂不出现。")
print(f"    但它是 **sp 金属**，且作者明说这正是 RMTA 的已知失效域")
print(f"    （逐字：`the inapplicability of the RMTA for sp metals`）。")
if abs(al_rmta/al_lr - 0.3182) > 0.001:
    fails.append("(D4) Al 比值与仓库记录的 0.318 不符")
print()
print(f"过渡金属侧本文给出的 RMTA 数值个数 = {len(SAVRASOV_TM_RMTA_NUMBERS)}")
print("  Nb 逐字：`We obtain similar results with the RMTA, in full agreement with")
print("           earlier calculations [17,24]`   —— 散文，无数")
print("  Mo 逐字：`Our linear-response calculations of α²F(ω) are found to be close")
print("           to our RMTA calculations and to earlier ones [30]`  —— 散文，无数")
print("⟹ 「过渡金属 RMTA 方法偏差」在本文**仍是未定**，本轮未关闭。")
if SAVRASOV_TM_RMTA_NUMBERS:
    fails.append("(D4) 过渡金属 RMTA 数值集合应为空")

# ====================== (D5) 预算：这个缺口值不值得再派活 ======================
hdr("(D5) 预算 —— 要翻掉 Nb₃Ge 的判决需要多大的 RMTA 高估？")
need = (1.0 / NB3GE_R) ** 2
print(f"r = {NB3GE_R} ⟹ 要推过 1 需 RMTA 对 Nb 基**高估** η 达 (1/{NB3GE_R})² = {need:.3f}×")
print(f"现有证据：Al 的 RMTA 是**低估** {al_lr/al_rmta:.2f}×（方向相反）；")
print( "          Radevych 的代际漂移最坏 1.021×（且 Nb/V 方向使 r 下降）。")
print(f"⟹ 已知任何因子都远小于 {need:.2f}× 且多数方向相反 ⟹ 该缺口**大概率不翻判决**。")
print("   但它目前只由散文支撑 ⟹ 判词必须写「未定」，不得写「已关闭」。")

print()
print(sep)
if fails:
    print(f"FAILED ({len(fails)} 项)：")
    for f in fails:
        print("  - " + f)
    raise SystemExit(1)
print("ALL CHECKS PASSED")
print(sep)
