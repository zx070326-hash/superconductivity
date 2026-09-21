"""
判别性预算比值：|signal(机制)| / Sigma_nu |partial O/partial nu| * Delta_nu
所有输入数字均从仓库既有笔记(numerical-evidence 或更高)的表格直接抄录，
本脚本只做减法/比值算术并打印，不重新求解任何 Eliashberg/算子问题。
来源见每个数字旁的笔记路径注释。级别：算术本身是 analytic(纯减法)，
被减的输入数字都是 numerical-evidence（除 C_inf/ceiling 常数是 rational-exact）。
"""

def show(name, num, den_terms, note=""):
    den = sum(den_terms.values())
    print(f"\n=== {name} ===")
    print(f"  numerator (mechanism swing) = {num:.4f}")
    for k, v in den_terms.items():
        print(f"    nu={k:30s} contributes {v:.4f}")
    print(f"  denominator (linear sum)    = {den:.4f}")
    if den > 0:
        print(f"  ratio (linear)              = {num/den:.3f}")
    import math
    denq = math.sqrt(sum(v * v for v in den_terms.values()))
    if denq > 0:
        print(f"  denominator (quadrature)    = {denq:.4f}")
        print(f"  ratio (quadrature)          = {num/denq:.3f}")
    if note:
        print(f"  note: {note}")


# ---------------------------------------------------------------------------
# Row 1: alpha, PdH -> PdD step, mechanism axis = Coulomb (mu*) amplification
# source: isotope_pdh_discriminating_test_2026-09-18.md Table 5.4 (alpha_budget
# at mu*=0 vs mu*=0.085), Table 5.3.1 (2D envelope width for proxy s / w_B),
# Section 3.3(b) (anchor r scan).
mech_swing_HD = abs(-0.399 - (-0.076))  # alpha_budget(mu*=0.085) - alpha_budget(mu*=0)
den_HD = {
    "proxy s / two-block envelope (w_B in [0.75,0.9023], Ω_ac scan, §5.3.1)": 0.390 - 0.274,
    "mu* anchor r (r0 in {5,10,20}, §3.3b)": 0.421 - 0.378,
    "bare vs reduced mass convention (§5.3.1a)": 0.0006,
}
show("alpha (PdH->PdD), mechanism=mu*/Coulomb", mech_swing_HD, den_HD,
     note="H->D 是本笔记里唯一给出完整二维包络 + 独立锚点扫描的对；"
          "分母的两项来自不同物理入口（谱形状代理 vs 库仑截断约定），"
          "线性求和是保守上界，二次方合成是乐观估计。")

# ---------------------------------------------------------------------------
# Row 2: alpha, PdD -> PdT step -- flagged UNDETERMINED (see main report)
mech_swing_DT = abs(0.3633 - 0.0418)  # alpha_ceiling(mu*=0) vs alpha_budget(mu*=.085), D->T
den_DT_partial = {
    "rounding band on Errea's 2-sig-fig Tc_D/Tc_T (deficit-budget quote, §4.2)": 0.2192 - 0.0638,
}
show("alpha (PdD->PdT), mechanism=mu*/Coulomb [UNDETERMINED]", mech_swing_DT, den_DT_partial,
     note="delta_BRD/形状项的 Delta_nu 在这一步没有独立来源——0.0637 的归因是"
          "'消去法'在单一谐波对照行上做的，不是一个带误差棒的 Delta_nu；"
          "本行只给出可解释的一部分分母，比值不可信，见主报告。")

# ---------------------------------------------------------------------------
# Row 3: filling factor f, deficit-budget (H3S/D3S), mechanism = harmonic vs
# anharmonic alpha^2F (a real microscopic contrast, already independently
# favored because harmonic implies unphysical mu*_fit~0.26-0.27)
# source: deficit_budget_fpred_vs_fobs_2026-09-18.md §1 table, §2 table
mech_swing_f = abs(0.6045 - 0.4792)  # f_pred(harmonic) - f_pred(anharmonic), H3S-200
den_f = {
    "pressure-anchor convention (legacy vs corrected Tc_exp, §2)": abs(0.4889 - 0.4432),
    "experimental digitization noise (sigma(ratio_obs)=3.81% propagated to f, §2)": 0.017,
    "mu*=0.13 fixed value vs literature spread +-0.02 (order-of-magnitude, f_of_mustar table)": 0.015,
}
show("filling factor f (H3S deficit-budget, NOT the circular hydride-table f)",
     mech_swing_f, den_f,
     note="samples: 只有 1 个真正独立自由度（H3S/D3S 共享系统误差，"
          "见 deficit_budget_fpred_vs_fobs 第3节）——比值再高，统计功效仍受限于样本量=1。")

print("\n" + "=" * 70)
print("circularity check: ambient_hydride_candidates.csv filling_factor_f")
print("  is IDENTICAL (48/48, <1.4% residual) to a closed-form F(lambda,r)")
print("  built purely from the Allen-Dynes formula at mu*=0.10 -- it carries")
print("  ZERO material information (hydride_csv_independent_audit_2026-09-17.md §3).")
print("  This is a DIFFERENT object from the H3S/D3S f above; the two must")
print("  never be pooled into one 'f' row.")
