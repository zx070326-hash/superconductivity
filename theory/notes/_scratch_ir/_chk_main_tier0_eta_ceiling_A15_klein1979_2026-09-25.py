#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tier-0 天花板，A15 化合物（强得多的检验：Tc 高一个量级）。

数据：B. M. Klein, L. L. Boyer, D. A. Papaconstantopoulos,
      PRL 42, 530 (1979), TABLE I（p.531 页图逐字抄）。
      N(E_F) [states/Ry-cell, 2 spins]; eta_A, eta_B [eV/A^2];
      lambda_calc, lambda_expt; Tc_calc, Tc_expt [K]。mu* = 0.13 全程。
取得：curl harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.42.530/fulltext (200)

⚠ **多位点约定是本脚本唯一的方向风险**，必须显式审：
  Hopfield/McMillan 对化合物：lambda = sum_j n_j * eta_j/(M_j <w_j^2>)，n_j = 晶胞内该种原子数。
  若各 <w_j^2> 同 => M2 = lambda<w^2> = sum_j n_j * eta_j/M_j。
  A15 是 A_3B => n_A = 3, n_B = 1。
  - 记 M2_min := eta_A/M_A + eta_B/M_B        （**不带**重数）
  - 记 M2_mult:= 3*eta_A/M_A + eta_B/M_B      （带重数）
  M2_min < M2_mult => 天花板更低 => **更难通过**。
  天花板是上界，低估 M2 只支持 FAIL、不支持 PASS —— 但反过来：
  **M2_min 下 PASS ⟹ M2_mult 下必 PASS**（单边界方向规则）。
  故主判决取 **M2_min**（保守/更锐），M2_mult 只作参考。

严格性：天花板 analytic；本判决 numerical-evidence（eta 是作者 1979 年 APW+Gaspari-Gyorffy 计算）。
"""
import math

C_INF = 0.1827262477
HBAR = 1.054571817e-34; EV = 1.602176634e-19
AMU = 1.66053906660e-27; KB = 1.380649e-23; ANG = 1.0e-10
HW_PER_SQRT = HBAR * math.sqrt((EV / ANG**2) / AMU) / EV
K_PER_EV = EV / KB
PRE = C_INF * HW_PER_SQRT * K_PER_EV          # = 137.0959 K per sqrt(eV/A^2/amu)

MASS = {"V": 50.942, "Nb": 92.906, "Al": 26.982, "Ga": 69.723,
        "Si": 28.085, "Ge": 72.630, "Sn": 118.710}

# (标签, A 元素, B 元素, eta_A, eta_B, lam_calc, lam_expt, Tc_calc, Tc_expt)
# Tc_expt = None 表示表中空白；">10" 记为 10.0 并标注（下限，对上界检验是保守的）
A15 = [
    ("V3Al",          "V",  "Al",  6.87, 0.03, 1.09, 0.87, 14.8,  9.6),
    ("V3Ga",          "V",  "Ga",  9.26, 0.05, 1.48, 1.17, 23.0, 16.5),
    ("V3Si",          "V",  "Si",  8.90, 0.07, 1.18, 1.12, 18.4, 17.1),
    ("V3Ge",          "V",  "Ge",  5.04, 0.05, 0.67, 0.70,  5.3,  6.1),
    ("V3Sn",          "V",  "Sn",  3.89, 0.04, 0.52, 0.61,  1.9,  3.8),
    ("Nb3Al",         "Nb", "Al", 13.40, 0.05, 2.14, 1.60, 25.1, 18.6),
    ("Nb3Ga",         "Nb", "Ga", 11.91, 0.13, 1.92, 1.74, 22.5, 20.3),
    ("Nb3Si a0=5.03", "Nb", "Si",  4.90, 0.05, 0.65, None,  3.5, None),
    ("Nb3Si a0=5.10", "Nb", "Si",  4.83, 0.05, 0.64, None,  3.4, 10.0),   # 表中 ">10"
    ("Nb3Si a0=5.20", "Nb", "Si",  4.84, 0.05, 0.64, None,  3.4, None),
    ("Nb3Ge",         "Nb", "Ge",  6.70, 0.07, 0.89, 1.80,  8.2, 23.2),
    ("Nb3Sn",         "Nb", "Sn",  8.49, 0.06, 1.12, 1.44, 12.6, 18.0),
]
LOWER_BOUND_TC = {"Nb3Si a0=5.10"}    # Tc_expt 是 ">10"，非等号


def main():
    print("Tc_ceiling[K] = %.4f * sqrt(M2[eV/A^2/amu])   C_inf=%.10f" % (PRE, C_INF))
    print("M2_min  = eta_A/M_A + eta_B/M_B        (不带重数，**保守**)")
    print("M2_mult = 3*eta_A/M_A + eta_B/M_B      (带 A15 重数 A_3B)")
    print()
    hdr = "%-15s %6s %6s %9s %9s %8s %8s %8s" % (
        "compound", "eta_A", "eta_B", "ceil_min", "ceil_mult", "Tc_exp", "r_min", "r_mult")
    print(hdr); print("-" * len(hdr))

    viol_min = viol_mult = n = 0
    worst = (0.0, "")
    etaB_frac = []
    for lbl, A, B, eA, eB, lc, le, tcc, tce in A15:
        MA, MB = MASS[A], MASS[B]
        m2_min = eA / MA + eB / MB
        m2_mul = 3.0 * eA / MA + eB / MB
        c_min, c_mul = PRE * math.sqrt(m2_min), PRE * math.sqrt(m2_mul)
        etaB_frac.append(eB / eA)
        if tce is None:
            print("%-15s %6.2f %6.2f %9.2f %9.2f %8s %8s %8s"
                  % (lbl, eA, eB, c_min, c_mul, "-", "-", "-"))
            continue
        r_min, r_mul = tce / c_min, tce / c_mul
        n += 1
        viol_min += (r_min > 1.0); viol_mult += (r_mul > 1.0)
        if r_min > worst[0]:
            worst = (r_min, lbl)
        mark = " (Tc 是下限 >10)" if lbl in LOWER_BOUND_TC else ""
        print("%-15s %6.2f %6.2f %9.2f %9.2f %8.1f %8.4f %8.4f%s"
              % (lbl, eA, eB, c_min, c_mul, tce, r_min, r_mul, mark))

    print()
    print("=== 判决 ===")
    print("  M2_min  [保守/更锐]: %d/%d 违反；最劣 %s r=%.4f" % (viol_min, n, worst[1], worst[0]))
    print("  M2_mult [带重数]   : %d/%d 违反" % (viol_mult, n))
    print("  方向：M2_min < M2_mult 逐行 => ceil_min < ceil_mult =>")
    print("        **M2_min 下 PASS 蕴含 M2_mult 下 PASS**（单边界方向规则）。")
    print()
    print("[作者自陈核对] 正文称 'In all cases eta_B <~ 0.01 eta_A'。")
    print("  实测 eta_B/eta_A：min=%.4f  max=%.4f  => %s"
          % (min(etaB_frac), max(etaB_frac),
             "与自陈相符" if max(etaB_frac) <= 0.0115 else
             "⚠ 最大 %.4f 超过 0.01（Nb3Ga），作者的 '<~' 是约数" % max(etaB_frac)))

    # 约定判别器：两种重数约定给出的 sqrt(<w^2>)，看哪个落在 A15 声子能标内
    print()
    print("=== 重数约定判别器：sqrt(M2/lambda_expt) = sqrt(<w^2>) ===")
    print("  物理约束：sqrt(<w^2>) <= w_max ~ Theta_D。A15 的 Theta_D 在 230-500 K 量级。")
    print("  %-15s %8s %12s %12s" % ("compound", "lam_expt", "w_min[K]", "w_mult[K]"))
    hi_min = hi_mul = 0.0
    for lbl, A, B, eA, eB, lc, le, tcc, tce in A15:
        if le is None:
            continue
        m2_min = eA / MASS[A] + eB / MASS[B]
        m2_mul = 3.0 * eA / MASS[A] + eB / MASS[B]
        w_min = math.sqrt(m2_min) * HW_PER_SQRT * K_PER_EV / math.sqrt(le)
        w_mul = math.sqrt(m2_mul) * HW_PER_SQRT * K_PER_EV / math.sqrt(le)
        hi_min = max(hi_min, w_min); hi_mul = max(hi_mul, w_mul)
        print("  %-15s %8.2f %12.1f %12.1f" % (lbl, le, w_min, w_mul))
    print("  最大值：w_min=%.0f K   w_mult=%.0f K" % (hi_min, hi_mul))
    print()
    print("  读法（**不要越级**）：")
    print("   - 带重数时 Nb3Sn 给 sqrt(<w^2>)=328 K，而 Nb3Sn 的 Theta_D 约 290 K")
    print("     => sqrt(<w^2>) > Theta_D，物理上**不可能**（sqrt(<w^2>) <= w_max ~ Theta_D）。")
    print("   - 不带重数时同一量为 189 K，与已知 A15 声子能标相容。")
    print("   => **弱证据**支持「eta_A 已是 A 子晶格的合计量、不再乘 3」，")
    print("      即 M2_min 不只是保守选择，很可能就是正确约定。")
    print("   ⚠ 但这**判不死**：lambda_expt 是作者用 McMillan 式**反解 Tc_expt** 得到的")
    print("     （表头逐字：'measured (by inverting the McMillan-type equation of Ref. 13)'），")
    print("     与 Tc_expt 有生成式关系 => 不是独立量。Theta_D 也未取自本文。")
    print("     约定问题的**定论**需读 Ref. 12（解耦假设的出处），本轮标记为**未定**。")
    print("   ⟹ 主判决仍取 M2_min，理由是**方向**（保守），不是上面这条弱证据。")


if __name__ == "__main__":
    main()
