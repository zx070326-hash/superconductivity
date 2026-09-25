#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主会话从零复核：「RMTA 在 Nb/Mo 上可用」这条判词的证据基础，以及由它派生的
                「Nb₃Ge r=0.626 落在 RMTA 可靠侧」这个头条加固。

**订正对象**：本机记忆 `rmta-underestimates-eta-for-sp-metals-by-3x` 写
  「RMTA 不是全局失效——作者明言 Nb『performed our own RMTA calculations and obtained
    complete agreement』、Mo『linear-response … close to our RMTA calculations』。
    失效是 sp 金属专属。」
并据此得出后果 3：A15 全是 Nb/V 基 ⟹ 头条落在 RMTA 可靠的那一侧。

**缺陷**：Nb 那句 `complete agreement` 比的**不是** RMTA vs 线性响应，
而是**他们自己的 RMTA vs 更早的 RMTA 计算**（implementation vs implementation）。
这是本轮第三起同形（前两起：Radevych 的 FLAPW-RMTA vs PS+PW-RMTA；
Savrasov PRL 72,372 的 `in full agreement with earlier calculations [17,24]`）。
⟹ Nb 那一半的证据**不存在**；只有 Mo 那句是真的跨方法陈述（且是散文、无数）。

**但同一段里有一个此前被漏掉的、真正的定量把手**，两篇论文各自独立给出：
文献中 RMTA 对 Nb 的 λ **跨 1.12–1.86**，而线性响应给 1.26。
⟹ 可以据此对「RMTA 高估 Nb 的 λ」取一个包封，并审它对头条的方向。

严格性：numerical-evidence（两端都是他人发表的值）。
本脚本**不 import** 任何既有脚本。

仪器（两台，逐字一致）：
  PRB 54, 16487 文本层 theory/notes/_scratch_ir/_savr_text.txt 第 576-577 行
                页 图 library/pdf/eta_hunt_2026-09-25/savr_rmta_p8-08.png（期刊页 16 494）
  PRL 72,  372 文本层 theory/notes/_scratch_ir/_savrasov_prl72_utf8.txt 第 252 行
                页 图 library/pdf/eta_hunt_2026-09-25/savrasov_prl72_tableI-3.png（p.374）
"""
import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

# ---- 逐字引文（两篇独立论文，同一个区间）----
QUOTE_PRB54 = (
    "Unfortunately, some of the results which have been reported in the literature are "
    "controversial. First, the RMTA calculations give the values of lambda varying from "
    "1.12 to 1.86."
)
QUOTE_PRB54_AGREEMENT = (
    "The discrepancy found by us has already been reported in the past RMTA-based "
    "calculations.[41,42] To check the consistency of our results with the earlier ones, "
    "we have performed our own RMTA calculations and obtained complete agreement between them."
)
QUOTE_PRL72 = (
    "Many theoretical RMTA based calculations exist in the literature [17,22,24] and a "
    "large variety of lambda (1.12-1.86) have been obtained."
)
QUOTE_PRL72_AGREEMENT = (
    "We obtain similar results with the RMTA, in full agreement with earlier "
    "calculations [17,24]"
)
QUOTE_MO = (
    "Our linear-response calculations of alpha^2 F(omega) are found to be close to our "
    "RMTA calculations and to earlier ones [30]."
)

RMTA_NB_LAMBDA_RANGE = (1.12, 1.86)      # 两篇都这么写
LR_NB_LAMBDA = 1.26                      # 线性响应（两篇 TABLE 的 lambda_calc 一致）
NB3GE_R = 0.6260                         # 仓库头条最紧点
FLIP_FACTOR = (1.0 / NB3GE_R) ** 2       # 要翻判决所需的 eta 高估倍数

fails = []
sep = "=" * 78


def hdr(t):
    print()
    print(sep)
    print(t)
    print(sep)


# ============ (A) 「Nb 上 RMTA 可用」的证据审计 ============
hdr("(A) 那句 `complete agreement` 的两端各是什么？")
print("PRB 54,16487 逐字（页图 savr_rmta_p8-08.png / 文本层第 576-577 行）：")
print("  " + QUOTE_PRB54_AGREEMENT)
print()
print("解析：`the consistency of our results with **the earlier ones**`，")
print("      而 the earlier ones 在前一句被明确定义为 `the past **RMTA-based** calculations`。")
print("      ⟹ 比较集合 = {他们自己的 RMTA, 更早的 RMTA} ⟹ **实现 vs 实现**。")
print()
print("PRL 72,372 逐字（同作者，1994）：")
print("  " + QUOTE_PRL72_AGREEMENT)
print("解析：`in full agreement with **earlier calculations**` —— 同一句式、同一集合。")
print()
print("⟹ 判词：**Nb 侧不存在 RMTA-vs-线性响应的一致性陈述。**")
print("   记忆里「Nb 上 RMTA 可用」这一半**没有证据支撑**，须撤回。")
print("   这是本轮第三起「实现一致 ≠ 方法正确」。")
print()
print("Mo 侧逐字（PRB 54,16487 / PRL 72,372 皆有）：")
print("  " + QUOTE_MO)
print("解析：`Our **linear-response** calculations ... close to our **RMTA** calculations`")
print("      ⟹ 这一句**确实**是跨方法陈述，两端正确。")
print("      但它是**散文、无数**（`close` 没有量化）⟹ 只支持定性、不进任何表。")
n_tm_quantified = 0
print(f"⟹ 过渡金属侧被**量化**的跨方法对照个数 = {n_tm_quantified}")
if n_tm_quantified != 0:
    fails.append("(A) 过渡金属量化对照个数应为 0")

# ============ (B) 此前漏掉的定量把手 ============
hdr("(B) 真正的定量把手：文献 RMTA 对 Nb 的 λ 跨 1.12–1.86")
lo, hi = RMTA_NB_LAMBDA_RANGE
print("PRB 54,16487 逐字：" + QUOTE_PRB54)
print("PRL 72,  372 逐字：" + QUOTE_PRL72)
print("⟹ **两篇独立论文给出同一个区间**，两台仪器（文本层 + 页图）各自确认。")
print()
print(f"线性响应（同作者、同代码）λ_calc(Nb) = {LR_NB_LAMBDA}")
print(f"文献 RMTA 区间 = [{lo}, {hi}]，宽度比 {hi/lo:.3f}×")
print(f"  最大**高估**比 = {hi}/{LR_NB_LAMBDA} = {hi/LR_NB_LAMBDA:.4f}×")
print(f"  最大**低估**比 = {LR_NB_LAMBDA}/{lo} = {LR_NB_LAMBDA/lo:.4f}×")
over = hi / LR_NB_LAMBDA
print()
print("⚠ 这是**区间宽度**，不是偏差：1.12–1.86 跨的是不同实现/不同势，")
print("   它给的是「RMTA 对 Nb 最多能偏多少」的包封，**不是**「RMTA 系统性偏多少」。")
print("   且 λ 到 η 仍带 ⟨ω²⟩ 混杂（除非同篇同约定）⟹ 本条只能当**指示性包封**。")

# ============ (C) 方向审计 + 对头条的影响 ============
hdr("(C) 方向审计：高估 η 才是危险方向")
print("ceiling = 137.0959·√(η/M)，r = Tc_exp/ceiling ∝ 1/√η。")
print("若 RMTA **高估** η 达 κ 倍 ⟹ 真 η 更小 ⟹ 真天花板更低 ⟹ 真 r = r_报 × √κ（**上升**）。")
print("⟹ 高估才是使判据变难的方向。这与记忆里「低估 ⟹ 保守 PASS」是**两个相反的方向**，")
print("   必须分开审：sp 金属那条是低估（安全），本条是高估（危险）。")
print()
r_worst = NB3GE_R * math.sqrt(over)
print(f"用 (B) 的最坏高估包封 κ = {over:.4f}：")
print(f"  r(Nb₃Ge) ≤ {NB3GE_R} × √{over:.4f} = {r_worst:.4f}")
print(f"  判据 r ≤ 1 ⟹ {'PASS' if r_worst < 1 else 'FAIL'}（余量 {1/r_worst:.3f}×）")
if r_worst >= 1:
    fails.append("(C) 在最坏 RMTA 高估包封下头条越过 1")
print()
print(f"翻判决所需的 η 高估倍数 = (1/{NB3GE_R})² = {FLIP_FACTOR:.4f}×")
print(f"  (B) 的包封 {over:.4f}× 小于它 {FLIP_FACTOR/over:.3f} 倍 ⟹ 仍有余量。")
if over >= FLIP_FACTOR:
    fails.append("(C) 包封已达到翻判决所需倍数")

# ============ (D) 净变化：记忆该怎么改 ============
hdr("(D) 净变化")
print("撤回：「RMTA 在 Nb 上可用（作者自陈 complete agreement）」")
print("      —— 那句比的是 RMTA vs 更早的 RMTA，**不是** RMTA vs 真值。")
print("撤回：后果 3 中「A15 全是 Nb/V 基 ⟹ 头条落在 RMTA 可靠侧」的**理由**。")
print("      （Nb 的可靠性没有证据；V 侧连散文都没有。）")
print()
print("新陈述（凭什么对）：头条**结论**仍然存活，但换成一条更强的理由——")
print(f"  不再依赖「RMTA 在 Nb 上准」，而是用**文献 RMTA 区间上端**做最坏包封：")
print(f"  即使取 λ 区间上端 {hi} 相对线性响应 {LR_NB_LAMBDA} 的全部高估（{over:.3f}×），")
print(f"  r(Nb₃Ge) ≤ {r_worst:.4f} < 1 仍 PASS。")
print("  ⟹ 从「假定它准」升级为「假定它最差也只差这么多」。这是更弱的前提、更强的结论。")
print()
print("保留不动：sp 金属 Al 的 η_RMTA/η_true = 0.14/0.44 = 0.318（低估 3.14×，无混杂）。")
print("          那条是同篇、同代码、同 ⟨ω²⟩ ⟹ 仍是唯一干净的方法偏差数。")

print()
print(sep)
if fails:
    print(f"FAILED ({len(fails)} 项)：")
    for f in fails:
        print("  - " + f)
    raise SystemExit(1)
print("ALL CHECKS PASSED")
print(sep)
