#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主会话从零复核：**第一篇通过验收门的文献** —— 同一篇论文里同时印出
过渡金属的 η_RMTA 与 η_非RMTA，带单位。

文献：S. K. Bose, "Linear response results for phonons and electron-phonon coupling
      in hcp Sc — spin fluctuations and implications for superconductivity",
      J. Phys.: Condens. Matter 20, 045209 (2008)；预印本 arXiv:0711.1543v1。
取得：curl https://arxiv.org/pdf/0711.1543 (200, 314106 B, 11 页)
      ⚠ TLS 现状已变：arxiv.org 本机现由 `Certainly / Certainly Intermediate R1`
        签发、python ssl 默认上下文**校验通过**，不再是记忆里记录的劫持态。
        全程**未使用** -k / --insecure。

两台仪器（CLAUDE.md §2 强制）：
  文本层 theory/notes/_scratch_ir/_bose0711_utf8.txt（pdftotext -enc UTF-8 -layout）
  页 图 library/pdf/eta_hunt_2026-09-25/bose_p7-07.png（pdftoppm -r 200，PDF 第 7 页）
  ⚠ 实测文本层有**行标错位**：ω_ln / ω_m / ω̄ / ω̃ 四行的**标签**相对数值下移一行。
    η 两行不受影响（两台仪器逐位一致），但这正是「两台仪器」规则的又一次兑现。

严格性：numerical-evidence（两端都是他人发表的计算值）。
本脚本**不 import** 任何既有脚本。
"""
import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

# ---------------- TABLE III 逐字（两台仪器一致）----------------
LATTICE_A_BOHR = [6.255, 6.15, 6.05, 5.85, 5.75]          # a (a_B)，随压力增大而减小
ETA_LR = [0.0381, 0.0433, 0.0519, 0.0666, 0.0771]         # η          Ry/bohr^2  (FP-LMTO 线性响应)
ETA_RMT = [0.0439, 0.0509, 0.0591, 0.0826, 0.1077]        # η(RMT/RAS) Ry/bohr^2  (LMTO-ASA 内的刚性原子球)
I2_LR = [0.0026, 0.0031, 0.0040, 0.0056, 0.0065]          # <I^2>          (Ry/bohr)^2
I2_RMT = [0.0030, 0.0037, 0.0045, 0.0065, 0.0081]         # <I^2>(RMT/RAS) (Ry/bohr)^2

QUOTE_CAPTION = (
    "TABLE III: Hopfield parameters from the linear response calculation eta and the "
    "rigid muffin-tin (atomic sphere) approximation eta (RMT/RAS), the Fermi surface "
    "average of the square of the electron-phonon matrix element <I^2> and "
    "<I^2>(RMT/RAS) from the linear response and the RMT results, respectively; ..."
)
QUOTE_DIRECTION = (
    "As the Hopfield parameter for transition metals has often been calculated using the "
    "rigid muffin-tin approximation of Gaspari and Gyorffy,31 in Table III we have compared "
    "the values obtained via the FP-LMTO linear response method (eta) and those obtained by "
    "using the RMT scheme implemented within the LMTO-ASA method,32,33 known as the rigid "
    "atomic sphere (RAS) method. The latter values, labeled as eta(RMT/RAS), are in general "
    "agreement with the FP-LMTO results, but are somewhat OVERESTIMATED, with the level of "
    "the overestimation increasing with pressure."
)
QUOTE_FE = (
    "This is in contrast with the results for the late transition metal hcp Fe, where RAS "
    "was shown to consistently UNDERESTIMATE the Hopfield parameter, with the level of the "
    "underestimation increasing with pressure.12"
)

NB3GE_R = 0.6260                      # 仓库头条最紧点（Klein 1979 A15，M2_min 档）
FLIP = (1.0 / NB3GE_R) ** 2           # 翻判决所需的 η 高估倍数
AL_RMTA_RATIO = 0.14 / 0.44           # 既有的唯一干净数：sp 金属 Al，低估

fails = []
sep = "=" * 78


def hdr(t):
    print()
    print(sep)
    print(t)
    print(sep)


# ================= (A) 验收门逐条比对 =================
hdr("(A) 这篇是否真的通过验收门？")
gate = [
    ("同一篇论文内", True, "TABLE III 一张表内并列两列"),
    ("同一种过渡金属", True, "hcp Sc（最轻的 3d 过渡金属），五个晶格常数"),
    ("印的是 η 本身而非 λ", True, "行标 `eta` 与 `eta (RMT/RAS)`，非 lambda"),
    ("带单位", True, "Ry/bohr^2（两行同单位）"),
    ("一端是 RMTA", True, "`the RMT scheme ... known as the rigid atomic sphere (RAS)`"),
    ("另一端非 RMTA", True, "`FP-LMTO linear response method`"),
]
for name, ok, why in gate:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name:<22} — {why}")
    if not ok:
        fails.append(f"(A) 验收门项失败：{name}")
print()
print("⟹ **这是本课题至今第一篇通过验收门的文献。**")
print("   此前过渡金属侧被量化的跨方法 η 对照个数 = 0（只有散文）。现在 = 1。")

# ================= (B) 逐列比值 =================
hdr("(B) 逐列 η 比值（RMTA / 线性响应）")
n = len(ETA_LR)
for lst, nm in ((ETA_RMT, "ETA_RMT"), (I2_LR, "I2_LR"), (I2_RMT, "I2_RMT"),
                (LATTICE_A_BOHR, "LATTICE_A_BOHR")):
    if len(lst) != n:
        fails.append(f"(B) 列数不一致：{nm}")
print("  %8s %10s %12s %9s %12s %14s" % (
    "a[a_B]", "eta_LR", "eta_RMT/RAS", "ratio", "<I2> ratio", "N(0)_LR/N(0)_RMT"))
ratios = []
for a, el, er, il, ir in zip(LATTICE_A_BOHR, ETA_LR, ETA_RMT, I2_LR, I2_RMT):
    ratio = er / el
    ratios.append(ratio)
    i2r = ir / il
    # eta = N(0)<I^2>  =>  N(0) = eta/<I^2>；两端 N(0) 若相同则 ratio == i2r
    n0_l, n0_r = el / il, er / ir
    print("  %8.3f %10.4f %12.4f %9.4f %12.4f %14.4f"
          % (a, el, er, ratio, i2r, n0_l / n0_r))
rmax, rmin = max(ratios), min(ratios)
print()
print(f"  比值范围 [{rmin:.4f}, {rmax:.4f}]；**全部 > 1 ⟹ RMTA 一致地高估**。")
if not all(r > 1.0 for r in ratios):
    fails.append("(B) 并非逐列高估")

# --- 作者自陈 `the level of the overestimation increasing with pressure` 的逐列核对 ---
# 我最初把它当成「逐列严格单调」来检验 —— 那是**我**给论文加的强度，论文只说趋势。
# 如实测：端点确实上升，但中间有一处真实的回落。
mono = all(ratios[i] < ratios[i + 1] for i in range(n - 1))
endpoints_up = ratios[-1] > ratios[0]
dips = [i for i in range(n - 1) if ratios[i + 1] < ratios[i]]
print(f"  作者自陈：`the level of the overestimation increasing with pressure`")
print(f"    端点（a 6.255→5.75）：{ratios[0]:.4f} → {ratios[-1]:.4f}  上升 {endpoints_up}")
print(f"    **逐列**严格单调：{mono}")
for i in dips:
    print(f"    ⟹ 回落处：a={LATTICE_A_BOHR[i]} 的 {ratios[i]:.4f}"
          f" → a={LATTICE_A_BOHR[i+1]} 的 {ratios[i+1]:.4f}")
print("    舍入检验（印出值 3 位有效数字的最坏区间）：")
print("      a=6.15 → [1.17301, 1.17803]；a=6.05 → [1.13667, 1.14079]，**两区间不相交**")
print("      ⟹ 这处回落**不是舍入伪影**，是真的。")
print("    ⟹ 作者那句是**粗趋势**陈述，逐列并不成立。此事不影响本轮承重量")
print("       （消费点只用 max(ratio)），但不得把它转述成「逐列单调」。")
if not endpoints_up:
    fails.append("(B) 端点未上升，与作者自陈的趋势方向都不符")
print()
print("  ⚠ **一个必须说清的混杂**：η 的比值 ≠ <I²> 的比值（上表末列 ≠ 1）")
print("    ⟹ 两端的 N(0) 也不同 ⟹ 这不是纯粹的「RMT vs 线性响应」，")
print("       而是 (RMT + LMTO-ASA) vs (线性响应 + FP-LMTO) 的**复合**差。")
print("    但**对本预算这恰好是想要的量**：我们要问的正是「若把 Pap1977 的 RMTA η")
print("    换成一次现代非 RMTA 计算，数字会动多少」——基组差本来就该算在内。")
print("    （对比：λ_Savr/λ_Pap 的 ⟨ω²⟩ 混杂是**有害**的，因为它把误差变成双边且")
print("     ⟨ω²⟩ 不可观测；这里的基组混杂是**无害**的，因为复合量就是消费点要的量。）")

# ================= (C) 方向审计 =================
hdr("(C) 方向审计：这一条是危险方向")
print("ceiling = 137.0959·√(η/M)，r = Tc_exp/ceiling ∝ 1/√η。")
print("RMTA **高估** η 达 κ 倍 ⟹ 真 η 小 ⟹ 真天花板低 ⟹ 真 r = r_报 × √κ（**上升**）。")
print()
print("三个已知的方法偏差，**方向并不统一**（不得取任何一个当全局判词）：")
print(f"  Al  (sp, Savrasov PRB54,16487 同篇同代码) η_RMTA/η_true = {AL_RMTA_RATIO:.4f}"
      f"  ⟹ 低估 {1/AL_RMTA_RATIO:.2f}×  → r 下降（安全）")
print(f"  Sc  (早 3d, 本文 TABLE III)                比值 ∈ [{rmin:.3f}, {rmax:.3f}]"
      f"  ⟹ 高估 ≤{rmax:.3f}×      → r 上升（**危险**）")
print("  Fe  (晚 3d, 本文引 [12]=PRB 67,214518)     作者自陈 consistently UNDERestimate"
      "  → r 下降（安全）")
print()
print("⟹ **不存在统一的 RMTA 偏差方向**；「低估所以 PASS 保守」只能逐材料说，不可推广。")
print("   本条的价值不是给方向，而是给**危险方向上的一个实测量级**。")

# ================= (D) 预算 =================
hdr("(D) 代入头条：Nb₃Ge")
r_worst = NB3GE_R * math.sqrt(rmax)
print(f"取 Sc 上实测的最大高估 κ = {rmax:.4f}（五列里最坏的一列，a=5.75 a_B）：")
print(f"  r(Nb₃Ge) ≤ {NB3GE_R} × √{rmax:.4f} = {r_worst:.4f}")
print(f"  判据 r ≤ 1 ⟹ {'PASS' if r_worst < 1 else 'FAIL'}（余量 {1/r_worst:.3f}×）")
if r_worst >= 1:
    fails.append("(D) 在 Sc 实测最坏高估下头条越过 1")
print()
print(f"翻判决所需 κ = (1/{NB3GE_R})² = {FLIP:.4f}×")
print(f"  实测包封 {rmax:.4f}× 比它小 {FLIP/rmax:.3f} 倍 ⟹ 判决不翻。")
if rmax >= FLIP:
    fails.append("(D) 实测包封已达翻判决倍数")
print()
print("外推限制（如实写）：Sc 不是 Nb/V，更不是 A15 化合物。")
print("  本条把「过渡金属侧完全无数」升级为「有一个同篇 η 级实测点」，")
print("  **不**等于「已测得 Nb 的 RMTA 偏差」。仍属指示性包封。")

# ================= (E) 顺带订正两处文件身份误判 =================
hdr("(E) 订正：两份被判为「误抓取」的 PDF 实为本课题最核心的两篇")
print("旧陈述错在哪：peer 报 `PhysRevLett.28.801.pdf` 是「液晶声模」、")
print("  `PhysRev.167.331.pdf` 是「TiO₂:V⁴⁺ 自旋回波」，判为误抓取、与课题无关。")
print("  成因：只读了 PDF 第 1 页开头几行。APS 按**期刊页**取件，页首是**前一篇文章**")
print("  在同一页上的**续尾**；本篇标题在同页稍下方开始。**页首 ≠ 本篇**。")
print()
print("新陈述凭什么对（全文 grep，非首行）：")
print("  PhysRevLett.28.801.pdf 第 75-77 行逐字：")
print("    `Electron-Phonon Interactions, d Resonances,`")
print("    `and Superconductivity in Transition Metals`")
print("    `G. D. Gaspari*f and B. L. Gyorffy`")
print("    ⟹ = Gaspari & Gyorffy, PRL 28, 801 (1972)，**RMTA/η 的奠基论文**，")
print("       即 Tier-0 语料里每一个 η 的方法出处。第 110 行并印出 λ=n(e_F)<I²>/M<ω²>。")
print("  PhysRev.167.331.pdf 第 25 行逐字：")
print("    `Transition Te1nperature of Strong-Coupled Superconductors`（OCR 的 m→1n）")
print("    ⟹ = McMillan, Phys. Rev. 167, 331 (1968)，**McMillan 方程**本身。")
print()
print("  ⟹ 两份都**不是**误抓取，是本课题最核心的两篇一手文献。不得删除。")
print("  教训同 [audit-the-artifact-not-the-stale-build]：判一份产物的身份，")
print("     要扫它的**全文**，不能用首行代表整份。")

# ================= (F) 引用链交叉核对 =================
hdr("(F) 顺带：本文参考文献独立印证了仓库的一处书目钉死")
print("  Bose ref [7] = S.Y. Savrasov, Phys. Rev. B54, **16470** (1996)")
print("  Bose ref [8] = S.Y. Savrasov and D.Y. Savrasov, Phys. Rev. B54, **16487** (1996)")
print("  ⟹ 第三方独立确认 16470 与 16487 是**两篇不同论文**（我此前的钉死成立）。")
print("  Bose ref [31] = G.D. Gaspari and B.L. Gyorffy, Phys. Rev. Lett. 28, 801 (1972)")
print("  Bose ref [32] = D. Glotzel, D. Rainer, H.R. Schober, Z. Phys. B 35, 317 (1979) [RAS]")
print("  Bose ref [12] = S.K. Bose, O.V. Dolgov, J. Kortus, O. Jepsen, O.K. Andersen,")
print("                  Phys. Rev. B67, 214518 (2003) [hcp Fe，反方向那篇]")

print()
print(sep)
if fails:
    print(f"FAILED ({len(fails)} 项)：")
    for f in fails:
        print("  - " + f)
    raise SystemExit(1)
print("ALL CHECKS PASSED")
print(sep)
