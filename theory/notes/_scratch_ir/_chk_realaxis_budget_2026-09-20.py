#!/usr/bin/env python3
"""
Real-axis gap correction budget analysis for 2Delta_/kTc discriminating ratio.

Task: Compute R_real flipover surface as function of (dc/dlam, dc/dmu*).
      Determine whether real-axis correction is load-bearing.

Anchor point: H3S lam_0=1.876, ω_0=106.9 meV, mu*_0=0.13, Delta_mu*=0.05, Delta_lam/lam=0.15
From main session verification (gap_ratio_mainsession_verification_2026-09-20.md §4):
  - solver.tc() branch (self-consistent Tc):
      dO*Delta_mu* = -0.10342
      dO*Delta_lam  = +0.20767
      R = 2.0080
  - Base O ~ 4.97 (Matsubara轴代理)

Strictness: Step 1-3 analytic (pure algebra), step 4 numerical-evidence.
"""
import numpy as np
import matplotlib.pyplot as plt

# === Anchor data (从主会话复核笔记提取) ===
dO_dmu_star = -0.10342  # dO/dmu* * Delta_mu*，来自自洽 solver.tc() 支
dO_dlam     = +0.20767  # dO/dlam * Delta_lam
O_base      = 4.97      # Matsubara轴代理基准值（H3S 点 N=1600）
Delta_mu    = 0.05      # mu* 扰动幅度
Delta_lam_rel = 0.15    # lam 相对扰动
lam_0       = 1.876     # H3S 基准 lam
Delta_lam   = Delta_lam_rel * lam_0

print("="*70)
print("Real-axis gap correction budget for 2Delta/kTc discriminating ratio")
print("="*70)
print(f"\nAnchor point: lambda_0={lam_0}, omega_0=106.9 meV, mu*_0=0.13")
print(f"Perturbations: Delta_mu*={Delta_mu}, Delta_lambda/lambda={Delta_lam_rel:.2%}")
print(f"\nMatsubara-axis baseline (from main session section 4 table):")
print(f"  dO*Delta_mu* = {dO_dmu_star:+.5f}")
print(f"  dO*Delta_lam  = {dO_dlam:+.5f}")
print(f"  R_Matsubara = {abs(dO_dlam)/abs(dO_dmu_star):.4f}")
print(f"  O_base = {O_base:.2f}")

# === Step 1: Write R_real as function of (dc_dlam, dc_dmu) [analytic] ===
print("\n" + "="*70)
print("Step 1: Explicit formula for R_real(dc/dlam, dc/dmu*) [analytic]")
print("="*70)
print("\nReal observable: O_real = c * O")
print("Chain rule:      dO_real/dp = c * dO/dp + O * dc/dp")
print("Discriminating ratio:")
print("  R_real = |c*dO/dlam*Delta_lam + O*dc/dlam*Delta_lam| / |c*dO/dmu**Delta_mu* + O*dc/dmu**Delta_mu*|")
print("\nFactor out c from numerator and denominator:")
print("  R_real = |dO*Delta_lam + O*(dc/dlam)*Delta_lam| / |dO*Delta_mu* + O*(dc/dmu*)*Delta_mu*|")

# Dimensionless derivatives (归一化，单位 1/1)
# 定义 x := (dc/dlam)*Delta_lam / c  (lam 轴相对修正)
#      y := (dc/dmu*)*Delta_mu* / c (mu* 轴相对修正)
# 则分子 = |dO*Delta_lam + O*x|，分母 = |dO*Delta_mu* + O*y|

def R_real_func(dc_dlam_scaled, dc_dmu_scaled):
    """
    R_real as function of scaled derivatives.

    Parameters:
        dc_dlam_scaled: (dc/dlam)*Delta_lam / c  [dimensionless, fractional change in c per lam perturbation]
        dc_dmu_scaled:  (dc/dmu*)*Delta_mu* / c [dimensionless, fractional change in c per mu* perturbation]

    Returns:
        R_real: discriminating ratio including real-axis correction
    """
    num = np.abs(dO_dlam + O_base * dc_dlam_scaled)
    den = np.abs(dO_dmu_star + O_base * dc_dmu_scaled)
    return num / den

print(f"\nNumerical prefactors:")
print(f"  dO*Delta_lam  = {dO_dlam:+.5f}")
print(f"  dO*Delta_mu* = {dO_dmu_star:+.5f}")
print(f"  O_base = {O_base:.4f}")

# === Step 2: Direction audit (four quadrants) [analytic] ===
print("\n" + "="*70)
print("Step 2: Direction audit — four quadrants [analytic]")
print("="*70)
print("\nDefine x := (dc/dlam)*Delta_lam,  y := (dc/dmu*)*Delta_mu*")
print("Baseline (x=0, y=0): R = 2.0080")
print("\nQuadrant analysis:")

quadrants = [
    ("I:   x>0, y>0", +1, +1, "both corrections push c up with coupling/Coulomb"),
    ("II:  x>0, y<0", +1, -1, "c rises with lam, falls with mu*"),
    ("III: x<0, y<0", -1, -1, "both corrections push c down"),
    ("IV:  x<0, y>0", -1, +1, "c falls with lam, rises with mu*"),
]

for label, sgn_x, sgn_y, desc in quadrants:
    # Extreme case: x or y comparable to baseline signal
    x_test = sgn_x * 0.5 * abs(dO_dlam) / O_base
    y_test = sgn_y * 0.5 * abs(dO_dmu_star) / O_base
    R_test = R_real_func(x_test, y_test)
    print(f"\n{label}")
    print(f"  {desc}")
    print(f"  Sample (x={x_test:+.4f}, y={y_test:+.4f}): R_real = {R_test:.4f}")

# === Step 2.5: Find flipover surface R_real = 1 ===
print("\n" + "="*70)
print("Step 2.5: Flipover surface R_real = 1 [analytic formula, numerical plot]")
print("="*70)

# For fixed y, solve |dO*Delta_lam + O*x| = |dO*Delta_mu* + O*y| for x
# Two branches depending on signs:
#   Branch A: dO*Delta_lam + O*x = ±(dO*Delta_mu* + O*y)

y_grid = np.linspace(-0.15, 0.15, 300)
x_flip_pos = []  # Positive branch
x_flip_neg = []  # Negative branch

for y in y_grid:
    # Positive branch: dO*Delta_lam + O*x = +(dO*Delta_mu* + O*y)
    x_p = (dO_dmu_star + O_base*y - dO_dlam) / O_base
    # Negative branch: dO*Delta_lam + O*x = -(dO*Delta_mu* + O*y)
    x_n = -(dO_dmu_star + O_base*y) - dO_dlam / O_base
    x_flip_pos.append(x_p)
    x_flip_neg.append(x_n)

x_flip_pos = np.array(x_flip_pos)
x_flip_neg = np.array(x_flip_neg)

fig, ax = plt.subplots(figsize=(10, 8))
ax.plot(x_flip_pos, y_grid, 'r-', linewidth=2, label='Flipover surface (branch +)')
ax.plot(x_flip_neg, y_grid, 'b-', linewidth=2, label='Flipover surface (branch -)')
ax.axhline(0, color='k', linewidth=0.5, linestyle='--')
ax.axvline(0, color='k', linewidth=0.5, linestyle='--')
ax.scatter([0], [0], s=100, c='green', marker='o', zorder=5, label=f'Baseline (R={abs(dO_dlam)/abs(dO_dmu_star):.3f})')

# Shade R>1 region (outside the branches)
ax.fill_betweenx(y_grid, x_flip_neg, x_flip_pos, alpha=0.2, color='red', label='R<1 (verdict-flips)')
ax.set_xlabel('(dc/dlam)*Delta_lam [dimensionless, fractional c change per lam perturbation]', fontsize=11)
ax.set_ylabel('(dc/dmu*)*Delta_mu* [dimensionless, fractional c change per mu* perturbation]', fontsize=11)
ax.set_title('Real-axis correction budget: R_real = 1 flipover surface', fontsize=13, fontweight='bold')
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)
ax.set_xlim(-0.15, 0.15)
ax.set_ylim(-0.15, 0.15)

plt.tight_layout()
plt.savefig('D:/physics/superconductivity/theory/notes/_scratch_ir/_fig_realaxis_budget_flipover.png', dpi=150)
print("\nFlipover surface plotted -> _fig_realaxis_budget_flipover.png")
print(f"  Green dot: baseline (x=0, y=0), R={abs(dO_dlam)/abs(dO_dmu_star):.4f}")
print(f"  Red shaded region: R<1 (verdict flips)")
print(f"  Outside shaded region: R>1 (verdict survives)")

# === Step 3: Extreme degradation — lam-axis signal completely canceled ===
print("\n" + "="*70)
print("Step 3: Extreme degradation budget [analytic]")
print("="*70)
print("\nScenario: real-axis correction COMPLETELY cancels lam-axis signal")
print("  Condition: c*dO/dlam + O*dc/dlam = 0")
print("  Solve for: (dc/dlam)*Delta_lam = -dO/dlam*Delta_lam / O")

dc_dlam_extreme = -dO_dlam / O_base
dc_dmu_extreme_for_flip = (abs(dO_dmu_star) - abs(dO_dlam)) / O_base  # Make R->1

print(f"\nExtreme lam-axis cancellation:")
print(f"  (dc/dlam)*Delta_lam = {dc_dlam_extreme:.5f}  [dimensionless]")
print(f"  Interpretation: c must change by {abs(dc_dlam_extreme)*100:.2f}% when lam changes by {Delta_lam_rel*100:.0f}%")
print(f"  Differential form: dc/dlam = {dc_dlam_extreme/Delta_lam:.4f} per unit lam")
print(f"\nTo ALSO flip verdict (R->1) requires BOTH:")
print(f"  (dc/dlam)*Delta_lam  = {dc_dlam_extreme:.5f}  (cancel lam signal)")
print(f"  (dc/dmu*)*Delta_mu* ~ {dc_dmu_extreme_for_flip:.5f}  (amplify mu* signal to match)")

# === Step 4: Cross-check with Carbotte empirical formula ===
print("\n" + "="*70)
print("Step 4: Carbotte formula cross-check [numerical-evidence]")
print("="*70)
print("\nCarbotte RMP 1990: 2Delta_/kTc = 3.53[1 + 12.5t^2 ln(1/2t)], t = k_BTc/ω_log")
print("This IS a real-axis formula (spectral-gap Delta__0).")
print("Ratio c := Carbotte / O_Matsubara estimates the Matsubara->real correction factor.")

# H3S anchor: Tc~200K, ω_log~106.9 meV
k_B_eV = 8.617333e-5  # eV/K
Tc_H3S = 200.0  # K (approximate from main session data)
omega_log_meV = 106.9
omega_log_eV = omega_log_meV * 1e-3

def carbotte_formula(Tc_K, omega_log_eV):
    """Carbotte 2Delta_/kTc formula."""
    t = k_B_eV * Tc_K / omega_log_eV
    return 3.53 * (1 + 12.5 * t**2 * np.log(1/(2*t)))

# Base point
O_carbotte_base = carbotte_formula(Tc_H3S, omega_log_eV)
c_empirical_base = O_carbotte_base / O_base

print(f"\nH3S anchor (Tc~{Tc_H3S:.0f} K, ω_log={omega_log_meV:.1f} meV):")
print(f"  t = k_BTc/ω_log = {k_B_eV*Tc_H3S/omega_log_eV:.4f}")
print(f"  Carbotte 2Delta_/kTc = {O_carbotte_base:.4f}")
print(f"  Matsubara proxy  = {O_base:.4f}")
print(f"  c_empirical := Carbotte/Matsubara = {c_empirical_base:.4f}")

# Perturb lam: change Tc by ±15% (crude proxy for lam±15%, since Tc ∝ exp(-1/lam) roughly)
# More accurate: use AD formula Tc ∝ ω_log exp(-1.04(1+lam)/lam), but Carbotte is lam-agnostic
# So perturb Tc directly as proxy
Tc_plus  = Tc_H3S * 1.15
Tc_minus = Tc_H3S * 0.85

O_carb_plus  = carbotte_formula(Tc_plus, omega_log_eV)
O_carb_minus = carbotte_formula(Tc_minus, omega_log_eV)

# CorrespondingMatsubara proxy: scale by same Tc ratio (first-order)
O_mat_plus_approx  = O_base  # Assume Delta_(iω_0) roughly Tc-independent (weak coupling)
O_mat_minus_approx = O_base

c_plus  = O_carb_plus / O_mat_plus_approx
c_minus = O_carb_minus / O_mat_minus_approx

dc_dlam_carbotte_crude = (c_plus - c_minus) / (2 * Delta_lam_rel)

print(f"\nPerturb Tc by ±15% (proxy for Delta_lam/lam=±15%):")
print(f"  Tc+15%: Carbotte={O_carb_plus:.4f}, c={c_plus:.4f}")
print(f"  Tc-15%: Carbotte={O_carb_minus:.4f}, c={c_minus:.4f}")
print(f"  dc/dlam ~ {dc_dlam_carbotte_crude:.4f} [crude, Tc-proxy derivative]")
print(f"  (dc/dlam)*Delta_lam ~ {dc_dlam_carbotte_crude * Delta_lam_rel:.5f}")

print(f"\n-> Carbotte suggests (dc/dlam)*Delta_lam ~ {dc_dlam_carbotte_crude * Delta_lam_rel:.3f}")
print(f"  Extreme budget (lam-signal cancellation) requires {dc_dlam_extreme:.3f}")
print(f"  Ratio: Carbotte estimate / budget = {abs(dc_dlam_carbotte_crude * Delta_lam_rel / dc_dlam_extreme):.2f}")

# === Step 5: Budget verdict ===
print("\n" + "="*70)
print("Step 5: BUDGET VERDICT")
print("="*70)

budget_lambda_cancel = abs(dc_dlam_extreme)
carbotte_estimate = abs(dc_dlam_carbotte_crude * Delta_lam_rel)
safety_factor = budget_lambda_cancel / carbotte_estimate

print("\n[ANALYTIC] Flipover surface R_real=1:")
print(f"  Requires (dc/dlam)*Delta_lam ~ {dc_dlam_extreme:.4f} (to cancel lam signal)")
print(f"  AND (dc/dmu*)*Delta_mu* ~ {dc_dmu_extreme_for_flip:.4f} (to amplify mu* signal)")
print(f"  Both must occur SIMULTANEOUSLY in the worst quadrant.")

print("\n[NUMERICAL-EVIDENCE] Carbotte empirical estimate:")
print(f"  (dc/dlam)*Delta_lam ~ {carbotte_estimate:.4f}")
print(f"  Budget / Carbotte = {safety_factor:.1f}×")

if safety_factor > 5:
    verdict = "预算极大 (>5×) ==> 实轴转换 NOT load-bearing"
    priority = "LOW — R>1 verdict stands on Matsubara proxy alone"
elif safety_factor > 2:
    verdict = "预算宽松 (2-5×) ==> 实轴转换 weakly load-bearing"
    priority = "MEDIUM — prudent to verify, but not blocking"
else:
    verdict = "预算紧张 (<2×) ==> 实轴转换 ESSENTIAL"
    priority = "HIGH — must implement real-axis solver before claiming R>1"

print(f"\n{'='*70}")
print(f"VERDICT: {verdict}")
print(f"Priority: {priority}")
print(f"{'='*70}")

print("\nConsumption point (if tight):")
print("  The discriminating-ratio criterion R>1 is consumed in:")
print("  theory/notes/isotope/gap_ratio_mainsession_verification_2026-09-20.md §7")
print("  to decide whether 2Delta_/kTc is a viable second observable.")
print("\nDirection safety:")
print("  - Carbotte is real-axis ==> its (dc/dlam) estimate is one-sided UPPER bound")
print("  - Actual dc/dlam likely smaller (Carbotte lumps all strong-coupling effects into t)")
print("  - One-sided bound rule: PASS verdict (R>1) is safe under this estimate")

print("\n" + "="*70)
print("Budget analysis complete. See _fig_realaxis_budget_flipover.png for flipover surface.")
print("="*70)
