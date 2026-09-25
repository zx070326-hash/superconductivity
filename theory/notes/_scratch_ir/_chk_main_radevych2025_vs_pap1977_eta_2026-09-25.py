#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
主会话从零复核：Radevych et al. (arXiv:2508.21281) TABLE I 的 eta_ref[14]
                与我们从 Papaconstantopoulos 1977 页图读出的 eta_tot(I) 对照。

两件事，必须分开写（它们是两个独立命题）：

  (A) 外部见证：Radevych 的 eta_ref[14] 引的就是 Pap1977（参考文献逐字为
      D. A. Papaconstantopoulos, L. L. Boyer, B. M. Klein, A. R. Williams,
      V. L. Morruzzi, J. F. Janak），于是它是对**我们抽取 Pap1977 数值**的
      一次独立见证 —— 特别是关掉「OCR 把 eta_tot(I) 列标签印成 &tot(II)、
      我们可能读错列」这个风险。

  (B) eta_new / eta_ref 的比值**不是** RMTA 的方法误差。
      两列都是 RMTA：eta 是他们的 PS+PW(赝势平面波) RMTA，
      eta_ref 是 Pap1977 的 APW+MJW RMTA。
      作者自己的话：`Some discrepancies with older references likely reflect
      the evolution of potentials and exchange-correlation functionals over time.`
      ⟹ 这个比值量的是**势/交换关联泛函的代际漂移**，不是 RMTA 相对全势的偏差。
      （RMTA-vs-FLAPW 的对照在该文别处，本脚本不碰。）

数据来源（两台仪器都跑过，逐字一致）：
  文本层 theory/notes/_scratch_ir/_radevych2508_utf8.txt (pdftotext -enc UTF-8 -layout)
  页 图 library/pdf/eta_hunt_2026-09-25/radevych_p8-08.png (pdftoppm -r 230, 第 8 页)

严格性：numerical-evidence（两端都是他人发表的 DFT 值；无 analytic 新命题、无证书）。

本脚本**不 import** 任何既有脚本；Pap1977 那一列是我从
_chk_main_tier0_eta_ceiling_papaconstantopoulos_2026-09-25.py 手抄过来的，
抄写本身由下面的 MANUAL_TRANSCRIPTION_NOTE 记录，比对是为了查「我们读的列对不对」。
"""
import math
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

MANUAL_TRANSCRIPTION_NOTE = (
    "Pap1977 列手抄自仓库既有脚本第 55-74 行的 eta_tot(I)；"
    "Radevych 列手抄自页图 radevych_p8-08.png 的 TABLE I。"
)

# --- Radevych et al. arXiv:2508.21281, TABLE I (p.8) ---
# 列：Structure | r_MT (bohr) | eta (eV/Ang^2) | eta_ref[14] (eV/Ang^2) | theta_D (K) | lambda
RADEVYCH_TABLE_I = [
    ("V",  "BCC", 2.48, 8.01, 6.89, 385, 1.19),
    ("Cu", "FCC", 2.41, 0.57, 0.46, 343, 0.09),
    ("Nb", "BCC", 2.70, 9.07, 7.63, 275, 1.45),
    ("Mo", "BCC", 2.57, 5.56, 5.80, 430, 0.35),
    ("Rh", "FCC", 2.56, 4.95, 4.72, 480, 0.24),
    ("Pd", "FCC", 2.66, 2.65, 1.99, 299, 0.31),
    ("Ag", "FCC", 2.73, 0.28, 0.29, 225, 0.06),
]

# --- 我们从 Pap1977 TABLE III 页图读出的 eta_tot(I) 与 eta_tot(II) (eV/Ang^2) ---
OURS_PAP1977 = {
    "V":  (6.894, 4.370),
    "Cu": (0.456, 0.380),
    "Nb": (7.627, 5.410),
    "Mo": (5.797, 3.440),
    "Rh": (4.722, 2.490),
    "Pd": (1.992, 1.050),
    "Ag": (0.288, 0.250),
}

C_INF = 0.1827262477          # = sqrt(g(2))/(2*pi)，analytic
NB3GE_R = 0.6260              # 仓库当前最紧点（Klein 1979 A15, M2_min 档）

fails = []


def close(a, b, tol):
    return abs(a - b) <= tol


print("=" * 78)
print("(A) 外部见证：Radevych eta_ref[14]  vs  我们读的 Pap1977 eta_tot(I)")
print("=" * 78)
print(f"{'元素':<4} {'ours_I':>8} {'ours_II':>8} {'their_ref':>10} "
      f"{'|I-ref|':>8} {'|II-ref|':>8}  判词")
n_match_I = 0
for (el, _st, _rmt, _eta_new, eta_ref, _td, _lam) in RADEVYCH_TABLE_I:
    oI, oII = OURS_PAP1977[el]
    dI, dII = abs(oI - eta_ref), abs(oII - eta_ref)
    # 他们印 3 位有效数字 ⟹ 舍入容差取 0.005（末位半个单位）
    hit_I = dI <= 0.005
    hit_II = dII <= 0.005
    if hit_I:
        n_match_I += 1
    verdict = "(I) 命中" if hit_I and not hit_II else (
        "(II) 命中" if hit_II and not hit_I else (
            "两列都命中(不可分辨)" if hit_I and hit_II else "都不命中"))
    print(f"{el:<4} {oI:>8.3f} {oII:>8.3f} {eta_ref:>10.2f} "
          f"{dI:>8.3f} {dII:>8.3f}  {verdict}")
    if not hit_I:
        fails.append(f"(A) {el}: 我们的 eta_tot(I)={oI} 与他们的 eta_ref={eta_ref} 差 {dI:.3f} > 0.005")
    if hit_II:
        fails.append(f"(A) {el}: eta_tot(II) 也命中 ⟹ 该元素不能用来分辨列")

print()
print(f"⟹ eta_tot(I) 命中 {n_match_I}/{len(RADEVYCH_TABLE_I)}；"
      f"eta_tot(II) 命中 0/{len(RADEVYCH_TABLE_I)}")
print("⟹ 我们读的是 (I) 列，外部见证确认。OCR 的 `&tot(II)` 标签错位未污染数值。")

print()
print("=" * 78)
print("(B) eta_new/eta_ref —— 势/泛函代际漂移，**不是** RMTA 方法误差")
print("=" * 78)
print(f"{'元素':<4} {'eta_new':>8} {'eta_ref':>8} {'比值':>7} "
      f"{'r 的乘子':>9}  方向")
ratios = {}
worst_mult, worst_el = 0.0, None
for (el, _st, _rmt, eta_new, eta_ref, _td, _lam) in RADEVYCH_TABLE_I:
    ratio = eta_new / eta_ref
    # ceiling ∝ sqrt(eta) ⟹ r = Tc/ceiling ∝ 1/sqrt(eta) ⟹ r 的乘子 = 1/sqrt(ratio)
    mult = 1.0 / math.sqrt(ratio)
    ratios[el] = (ratio, mult)
    direction = "r 降（安全）" if mult < 1 else "r 升（⚠ 方向不利）"
    print(f"{el:<4} {eta_new:>8.2f} {eta_ref:>8.2f} {ratio:>7.3f} "
          f"{mult:>9.3f}  {direction}")
    if mult > worst_mult:
        worst_mult, worst_el = mult, el

print()
print(f"最不利元素：{worst_el}，r 的乘子 = {worst_mult:.4f} "
      f"（即 r 最多上升 {100*(worst_mult-1):.1f}%）")

# 方向门：新势下 r 最多涨多少，会不会把任何一个点推过 1？
# 我们的元素语料最紧点是 In，r=0.4460（In 不在这 7 个里，但取最坏乘子做包封）
IN_R = 0.4460
print(f"包封检验：元素栏最紧 In r={IN_R} × 最坏乘子 {worst_mult:.4f} "
      f"= {IN_R*worst_mult:.4f} < 1  ⟹ {'PASS' if IN_R*worst_mult < 1 else 'FAIL'}")
if IN_R * worst_mult >= 1:
    fails.append("(B) 元素栏在最坏乘子下越过 1")

# Nb 是 A15 语料的主元素，单独报
nb_ratio, nb_mult = ratios["Nb"]
v_ratio, v_mult = ratios["V"]
print()
print(f"A15 相关：Nb 比值 {nb_ratio:.3f}（r 乘子 {nb_mult:.3f}）、"
      f"V 比值 {v_ratio:.3f}（r 乘子 {v_mult:.3f}）")
print(f"  若把元素 Nb 的漂移**指示性地**搬到 Nb3Ge："
      f"r = {NB3GE_R} × {nb_mult:.3f} = {NB3GE_R*nb_mult:.4f}")
print("  ⚠ 这只是指示：元素 Nb 的势漂移不能严格外推到化合物 Nb3Ge。")
print("     它给的是**方向**（新势下 eta 更大 ⟹ 天花板更高 ⟹ r 更小），不是新的判决。")

print()
print("=" * 78)
print("(C) 单位链自检：C_INF -> 137.0959")
print("=" * 78)
# kTc <= C_inf * sqrt(M2)，M2 = eta/M，eta 单位 eV/Ang^2，M 单位 amu
# sqrt(eta/M) 的量纲是 sqrt(eV/Ang^2/amu) = 能量/hbar 的量纲 -> 转成 K
EV = 1.602176634e-19          # J
ANG = 1.0e-10                 # m
AMU = 1.66053906892e-27       # kg
HBAR = 1.054571817e-34        # J s
KB = 1.380649e-23             # J/K
# sqrt(eV/Ang^2 / amu) 给出角频率 [rad/s]
omega_unit = math.sqrt(EV / (ANG**2) / AMU)
# kTc <= C_inf * hbar * omega  ⟹ Tc[K] <= C_inf * hbar*omega/kB
coeff = C_INF * HBAR * omega_unit / KB
print(f"sqrt(1 eV/Ang^2 / 1 amu) = {omega_unit:.6e} rad/s")
print(f"C_INF * hbar * that / kB = {coeff:.4f}  K   (仓库记录值 137.0959)")
if not close(coeff, 137.0959, 0.02):
    fails.append(f"(C) 单位链不符：重算 {coeff:.4f} vs 记录 137.0959")
else:
    print("⟹ 单位链 PASS（与仓库记录的 137.0959 一致）")

print()
print("=" * 78)
if fails:
    print(f"FAILED ({len(fails)} 项)：")
    for f in fails:
        print("  - " + f)
    raise SystemExit(1)
print("ALL CHECKS PASSED")
print("=" * 78)
print()
print("判词（严格性：numerical-evidence）")
print("-" * 78)
print("1. 我们对 Pap1977 eta_tot(I) 的抽取被外部独立见证确认，7/7 逐元素命中，")
print("   且 eta_tot(II) 0/7 ⟹ 列归属无歧义。OCR 标签错位未污染数值。")
print("2. Radevych TABLE I 的 eta 与 eta_ref **两列都是 RMTA** ⟹ 它们的比值")
print("   量的是势/泛函代际漂移，**不能**当作 RMTA 的方法误差报。")
print("   ⟹ 「过渡金属 RMTA 方法偏差」这个缺口**本脚本没有关闭**，仍是未定。")
print("3. 代际漂移的方向在 7 个元素里 5 升 2 降；最不利是 Mo/Ag，r 最多 +2.1%。")
print("   包封后元素栏仍 PASS。A15 主元素 Nb/V 的漂移方向使 r **下降**。")
