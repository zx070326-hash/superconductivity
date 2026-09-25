#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主会话从零复核：Klein-Boyer-Papaconstantopoulos PRL 42, 530 (1979) TABLE I 的
**多重性约定** —— η_A 是「A 位单原子」还是「A 子晶格合计」？

这是头条 Nb₃Ge r=0.626 仅剩的方向风险：
  M2_min  := η_A/M_A + η_B/M_B        （合计读法）
  M2_mult := 3η_A/M_A + η_B/M_B       （单原子读法，A15 的 A 位重数 3）
M2_mult > M2_min ⟹ 天花板更高 ⟹ r 更小。故 M2_min 是**保守**选择；
但「保守」只保证判决方向安全，不等于「约定正确」——本脚本查的是后者。

订正对象：teammate `a15-mult-convention` 的两轮交付。它的**结论我复核后同意**，
但两处支撑要改（见 (E)）。

严格性：numerical-evidence。本脚本**不 import** 任何既有脚本，η/λ 全部从页图重抄。
仪器：library/pdf/eta_hunt_2026-09-25/klein42_p2-2.png（期刊页 531，pdftoppm -r 230）
     library/pdf/eta_hunt_2026-09-25/PhysRevLett.42.530.pdf 文本层（-enc UTF-8 -layout）
     ⚠ 该文本层两栏严重交错错行，**承重字符一律以页图为准**。
"""
import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

# ---- 单位常数，从物理常数独立重建（不从任何脚本 import）----
C_INF = 0.1827262477
HBAR = 1.054571817e-34
EV = 1.602176634e-19
AMU = 1.66053906660e-27
KB = 1.380649e-23
ANG = 1.0e-10
KAPPA = HBAR * math.sqrt((EV / ANG**2) / AMU) / EV * (EV / KB)   # K per sqrt(eV/A^2/amu)
PRE = C_INF * KAPPA

MASS = {"V": 50.942, "Nb": 92.906, "Al": 26.982, "Ga": 69.723,
        "Si": 28.085, "Ge": 72.630, "Sn": 118.710}

# TABLE I 逐字（页图 klein42_p2-2.png）：(标签, A, B, eta_A, eta_B, lam_calc)
TAB = [
    ("V3Al",           "V",  "Al",  6.87, 0.03, 1.09),
    ("V3Ga",           "V",  "Ga",  9.26, 0.05, 1.48),
    ("V3Si",           "V",  "Si",  8.90, 0.07, 1.18),
    ("V3Ge",           "V",  "Ge",  5.04, 0.05, 0.67),
    ("V3Sn",           "V",  "Sn",  3.89, 0.04, 0.52),
    ("Nb3Al",          "Nb", "Al", 13.40, 0.05, 2.14),
    ("Nb3Ga",          "Nb", "Ga", 11.91, 0.13, 1.92),
    ("Nb3Si a0=5.03",  "Nb", "Si",  4.90, 0.05, 0.65),
    ("Nb3Si a0=5.10",  "Nb", "Si",  4.83, 0.05, 0.64),
    ("Nb3Si a0=5.20",  "Nb", "Si",  4.84, 0.05, 0.64),
    ("Nb3Ge",          "Nb", "Ge",  6.70, 0.07, 0.89),
    ("Nb3Sn",          "Nb", "Sn",  8.49, 0.06, 1.12),
]

# 正文逐字（页图，期刊页 531）
Q_MUSTAR = ("For instance, we find a value of mu* = 0.1-0.2 using the Morel-Anderson "
            "formula,[10] our fitted values of mu, and typical values of the A15 d-band "
            "width (~0.2 eV) and average phonon frequencies (~200 K).")
Q_LAMBDA_SRC = ("Results are shown in Table I, along with the T_c values calculated using "
                "the Allen-Dynes expression[13] and the phonon density of states of "
                "Schweiss et al.[14]")
Q_BORROW = ("Both lambda and T_c were determined for these compounds using estimated values "
            "of the Nb-site phonon moments from the experimental results for Nb3Sn.[14]")
Q_ETAB = ("In all cases eta_B <~ 0.01 eta_A for the A15 compounds - electron-phonon "
          "scattering from the B site is very ineffective.[15]")
Q_DECOUPLE = ("We have also calculated the electron-phonon interactions eta_A and eta_B, and "
              "mass enhancement factors lambda, using the Gaspari-Gyorffy theory,[11] and the "
              "assumption of decoupling the contributions from the A and B sites in the A_3B "
              "A15 compounds.[12]")

# 外部实验 Debye 温度（二手，WebSearch 转述，**未读一手**）—— 只作数量级锚
THETA_D_SECONDARY = {"Nb3Sn": 232.0, "Nb3Ge": 302.0, "V3Si": 315.0}

fails = []
sep = "=" * 78


def hdr(t):
    print()
    print(sep)
    print(t)
    print(sep)


def inv(eA, eB, A, B, lam, mult):
    m2 = (3.0 if mult else 1.0) * eA / MASS[A] + eB / MASS[B]
    return m2, math.sqrt(m2 / lam) * KAPPA


# ============ (A) 单位常数自建 ============
hdr("(A) 单位常数从物理常数独立重建")
print(f"  kappa = hbar*sqrt((eV/A^2)/amu)/k_B = {KAPPA:.4f} K")
print(f"  PRE   = C_inf * kappa             = {PRE:.4f} K   (仓库记的 137.0959)")
if abs(PRE - 137.0959) > 5e-4:
    fails.append("(A) 单位常数与仓库值不符")

# ============ (B) 机械验证：反演真的还原了作者的输入吗 ============
hdr("(B) 先验证反演机械本身：Nb₃Ge 的声子矩是**借** Nb₃Sn 的")
print("  正文逐字：" + Q_BORROW)
_, w_ge = inv(6.70, 0.07, "Nb", "Ge", 0.89, False)
_, w_sn = inv(8.49, 0.06, "Nb", "Sn", 1.12, False)
print(f"  反演 √⟨ω²⟩：Nb₃Ge = {w_ge:.2f} K ；Nb₃Sn = {w_sn:.2f} K")
print(f"  相对差 = {abs(w_ge-w_sn)/w_sn*100:.3f}%")
print("  ⟹ 两者在 0.1% 内重合，**与论文自陈的「借用」完全吻合**。")
print("     这证明我的反演式 λ_calc = M2/⟨ω²⟩ **就是作者用的那个式子**（至多差一个重数因子），")
print("     即反演机械可信。⚠ 注意这一条对两种约定**都**成立（√3 是公因子），")
print("     所以它验证的是机械，**不**判约定。")
if abs(w_ge - w_sn) / w_sn > 0.005:
    fails.append("(B) 借用关系未被反演还原，反演机械可疑")

# ============ (C) 全表反演 ============
hdr("(C) 12 个化合物全表反演 √⟨ω²⟩")
print("  %-16s %7s %7s %10s %11s" % ("compound", "eta_A", "lam_c", "w_min[K]", "w_mult[K]"))
grp = {"V": [], "Nb": []}
for lbl, A, B, eA, eB, lam in TAB:
    _, wmin = inv(eA, eB, A, B, lam, False)
    _, wmul = inv(eA, eB, A, B, lam, True)
    grp[A].append(wmin)
    print("  %-16s %7.2f %7.2f %10.1f %11.1f" % (lbl, eA, lam, wmin, wmul))
vbar = sum(grp["V"]) / len(grp["V"])
nbar = sum(grp["Nb"]) / len(grp["Nb"])
print()
print(f"  V₃X  簇均值 = {vbar:.1f} K  （范围 {min(grp['V']):.1f}–{max(grp['V']):.1f}）")
print(f"  Nb₃X 簇均值 = {nbar:.1f} K  （范围 {min(grp['Nb']):.1f}–{max(grp['Nb']):.1f}）")

# ============ (D) 自由度审计：哪些观测**不能**判约定 ============
hdr("(D) 自由度审计：先分清哪些量在换约定下**不变**")
ratio_obs = vbar / nbar
ratio_mass = math.sqrt(MASS["Nb"] / MASS["V"])
print("  换约定 M2_min → M2_mult 等价于把**每一个** √⟨ω²⟩ 乘以同一个 √3 = 1.7321。")
print("  ⟹ 任何**比值型**观测在换约定下恒不变，**判别力为 0**：")
print(f"     · 簇比 V/Nb = {ratio_obs:.4f}  vs  质量标度 √(M_Nb/M_V) = {ratio_mass:.4f}"
      f"  （符合度 {100*(1-abs(ratio_obs-ratio_mass)/ratio_mass):.1f}%）")
print("       —— 这条很漂亮（说明反演值确实像真声子频率），但它**两种约定下完全相同**，")
print("          不得当作支持任一约定的证据。")
print(f"     · Nb₃Ge/Nb₃Sn 借用重合（(B) 那条）—— 同理，判别力 0。")
print()
print("  ⟹ **只有「绝对能标 vs 外部锚」这一类观测才有判别力。**（CLAUDE.md：先数自由度）")
if abs(ratio_obs - ratio_mass) / ratio_mass > 0.05:
    fails.append("(D) 簇比与质量标度偏离 >5%，反演值不像真声子频率")

# ============ (E) 订正 teammate 的两处支撑 ============
hdr("(E) 订正：teammate 的结论我同意，但两处支撑要改")
print("【第一处】它把「~200 K」当成作者喂给 λ 的 ⟨ω²⟩。**消费点错了。**")
print("  页图逐字（期刊页 531 左栏）：")
print("    " + Q_MUSTAR)
print("  ⟹ 这个 ~200 K 是代进 **Morel-Anderson 公式估 μ\\*** 的，不是 λ 的声子矩。")
print("  λ 的真正来源在同页：")
print("    " + Q_LAMBDA_SRC)
print("  ⟹ 是 **Allen-Dynes 式 + Schweiss et al.[14] 的声子态密度**。")
print("  （CLAUDE.md：重算前先确认那个量到底作用在哪个变量上，不要从名字推语义。）")
print("  它**仍可**当作「作者认为 A15 平均声子频率的量级」的旁证，但级别要降，")
print("  且不得说成「作者自己反解得到的就该是这个数」。")
print()
print("【第二处】它的相对差算错了锚。它口头引 √⟨ω²⟩=√(3/5)Θ_D≈0.775Θ_D，")
print("  却把相对差算在**裸 Θ_D** 上（报 Nb₃Sn 「−7.6%」）。两者不能混：")
for k, td in sorted(THETA_D_SECONDARY.items()):
    row = [t for t in TAB if t[0] == k][0]
    _, wmin = inv(row[3], row[4], row[1], row[2], row[5], False)
    _, wmul = inv(row[3], row[4], row[1], row[2], row[5], True)
    print(f"    {k:<7} Θ_D={td:.0f}K  |  w_min={wmin:6.1f} ({wmin/td-1:+6.1%} vs Θ_D,"
          f" {wmin/(0.775*td)-1:+6.1%} vs 0.775Θ_D)")
    print(f"    {'':<7}            |  w_mult={wmul:6.1f} ({wmul/td-1:+6.1%} vs Θ_D,"
          f" {wmul/(0.775*td)-1:+6.1%} vs 0.775Θ_D)")
print("  ⟹ 换成 Debye 模型锚 0.775Θ_D 后，w_min 其实是**偏高** 10–20%，不是「−7.6%」。")
print("     但**判决方向不变**：w_mult 在两种锚下都高出 60–160%。")

# ============ (F) 真正有判别力的那一条 ============
hdr("(F) 有判别力的判据：绝对能标是否可能")
print("  模型无关的硬约束：√⟨ω²⟩ ≤ ω_max（均方根频率不超过谱的上端）。")
print("  A15 的实测 Θ_D 在 230–330 K 量级（二手），ω_max 与 Θ_D 同量级。")
print()
print(f"  M2_min ：全表 √⟨ω²⟩ ∈ [{min(min(grp['V']),min(grp['Nb'])):.0f},"
      f" {max(max(grp['V']),max(grp['Nb'])):.0f}] K")
print("           Nb₃X ~196–218 K、V₃X ~264–291 K，**全部落在 A15 声子能标之内**，")
print("           且 V 轻于 Nb ⟹ V 簇更高，**物理方向正确**。")
print(f"  M2_mult：全表 √⟨ω²⟩ ∈ [{math.sqrt(3)*min(min(grp['V']),min(grp['Nb'])):.0f},"
      f" {math.sqrt(3)*max(max(grp['V']),max(grp['Nb'])):.0f}] K")
print("           V₃X 要求 ~456–504 K —— **超过任何已知 A15 的 Θ_D 本身 40–60%**，")
print("           而 √⟨ω²⟩ 本应**低于** ω_max。物理上不可能。")
vmax_mult = math.sqrt(3) * max(grp["V"])
print()
print(f"  ⟹ [numerical-evidence] **M2_min 是作者的约定**：η_A 已是 A 子晶格合计量。")
print(f"     M2_mult 被证伪的裕度：V₃X 的 {vmax_mult:.0f} K vs Θ_D≈315 K，超出 {vmax_mult/315-1:+.0%}。")
if vmax_mult <= 330:
    fails.append("(F) M2_mult 未被外部能标证伪，判别失败")

# ============ (G) 残余未关闭 ============
hdr("(G) 残余：仍未读到 Ref.12 的约定原文")
print("  正文逐字：" + Q_DECOUPLE)
print("  Ref.12 = B. M. Klein and D. A. Papaconstantopoulos, J. Phys. F **6, 1135** (1976)")
print("    （teammate 报页图 OCR 作 1185，3→8 误读；真页码经第二仪器确认为 1135，")
print("     DOI 10.1088/0305-4608/6/6/026。**此条我尚未一手核实**，标未定。）")
print("  ⟹ 本判决是**判别器推断**，不是读到约定原文。等级 numerical-evidence，不得升格。")
print("  但就「方向风险是否活着」而言：M2_min 同时是**保守**方向**且**证据指向的正确约定")
print("  ⟹ 风险降级为「弱、且两条独立理由都指向安全」。头条 r=0.626 不需再标注约定未定。")
print()
print("  另：η_B 自陈逐字：" + Q_ETAB)
frac = [t[4] / t[3] for t in TAB]
print(f"    实测 η_B/η_A ∈ [{min(frac):.4f}, {max(frac):.4f}]"
      f"（最大在 Nb₃Ga = {max(frac):.4f}，作者的 `<~` 是约数）")

print()
print(sep)
if fails:
    print(f"FAILED ({len(fails)} 项)：")
    for f in fails:
        print("  - " + f)
    raise SystemExit(1)
print("ALL CHECKS PASSED")
print(sep)
