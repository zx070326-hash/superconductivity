"""Two-scale vs two-band investigation: quantify the connection between
the μ*>0 extremality counter-example shape (w=0.99/sp=1e5) and two-band
materials like MgB2.

Tasks:
1. Define two-scale parameter space (λ allocation ratio vs frequency ratio)
2. Map counter-examples and real materials onto this space
3. Scan deficit budget for two-scale spectra
4. Check if failure boundary correlates with two-scale-ness

Strictness: numerical-evidence (all float64, finite N).
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import from h3s_deficit_budget
from h3s_deficit_budget import Tc_P, Q_P, CINF_ as CINF, MEV_K

# ============================================================================
# 1. TWO-SCALE PARAMETERIZATION
# ============================================================================

def two_scale_spectrum(lam_total, om2_total, w, freq_ratio):
    """Construct a two-atom spectrum with specified λ allocation and frequency ratio.

    Args:
        lam_total: total λ
        om2_total: total ⟨ω²⟩
        w: weight on mode 1 (∈[0,1])
        freq_ratio: ω₂/ω₁ (>1, mode 2 is harder)

    Returns:
        (omegas, a2f_values) where P(dω) = (2/λ) α²F/ω dω

    λ allocation: λ₁ = w*lam_total to mode 1
    Frequency constraint: w*ω₁² + (1-w)*ω₂² = om2_total with ω₂ = freq_ratio*ω₁
    """
    # Solve for ω₁ from ⟨ω²⟩ constraint
    # w*ω₁² + (1-w)*(r*ω₁)² = om2_total
    # ω₁²[w + (1-w)*r²] = om2_total
    r = freq_ratio
    om1_sq = om2_total / (w + (1 - w) * r**2)
    om1 = math.sqrt(om1_sq)
    om2 = r * om1

    # λ allocation: λᵢ = (2/ωᵢ) ∫ α²F dω for δ spectrum
    # For δ_{ωᵢ}, α²F = (λᵢ*ωᵢ/2) δ(ω-ωᵢ)
    lam1 = w * lam_total
    lam2 = (1 - w) * lam_total

    # Return as (ω, α²F) pairs suitable for Tc_P
    omegas = np.array([om1, om2])
    # α²F values: for δ spectrum, integrate to λᵢ*ωᵢ/2
    a2f = np.array([lam1 * om1 / 2.0, lam2 * om2 / 2.0])

    return omegas, a2f


def two_scale_metrics(w, freq_ratio):
    """Dimensionless metrics characterizing "two-scale-ness".

    Returns:
        λ_asymmetry: |λ₁ - λ₂|/λ_total = |2w - 1|
        freq_separation: log₁₀(ω₂/ω₁)
        variance_carrier: fraction of ⟨ω²⟩ in minority mode
    """
    lam_asym = abs(2*w - 1)  # 0 = equal split, 1 = all in one mode
    freq_sep = math.log10(freq_ratio)

    # Variance carrier: if w is large (mode 1 has most λ),
    # then mode 2 (1-w weight, freq_ratio² scaling) carries variance
    r = freq_ratio
    om2_fraction = (1 - w) * r**2 / (w + (1 - w) * r**2)

    return lam_asym, freq_sep, om2_fraction


# ============================================================================
# 2. MAP COUNTER-EXAMPLES AND MATERIALS
# ============================================================================

print("="*80)
print("PART 1: Two-scale characterization of counter-examples and materials")
print("="*80)

# Counter-example from extremality_failure_map
ce_w = 0.99
ce_sp = 1e5
ce_lam_asym, ce_freq_sep, ce_var_frac = two_scale_metrics(ce_w, math.sqrt(ce_sp))
print(f"\nCounter-example (w={ce_w}, spread={ce_sp}):")
print(f"  λ asymmetry: {ce_lam_asym:.3f}")
print(f"  freq separation: {ce_freq_sep:.2f} decades")
print(f"  variance carrier (minority fraction): {ce_var_frac:.4f}")

# Real materials: need to extract from actual α²F data
# For now, characterize Einstein spectrum (reference)
print(f"\nEinstein spectrum (single frequency):")
print(f"  λ asymmetry: 1.000 (singular)")
print(f"  freq separation: undefined (single peak)")
print(f"  variance carrier: 0.000 (no split)")

# MgB2 two-band: approximate from literature
# σ-band: ~60% of λ at ~60 meV, π-band: ~40% at ~30 meV (rough estimates)
mgb2_w = 0.6  # σ-band weight
mgb2_om_ratio = 60.0 / 30.0  # ~2
mgb2_lam_asym, mgb2_freq_sep, mgb2_var_frac = two_scale_metrics(mgb2_w, mgb2_om_ratio)
print(f"\nMgB2 two-band estimate (w={mgb2_w}, ω_σ/ω_π≈{mgb2_om_ratio}):")
print(f"  λ asymmetry: {mgb2_lam_asym:.3f}")
print(f"  freq separation: {mgb2_freq_sep:.2f} decades")
print(f"  variance carrier: {mgb2_var_frac:.4f}")

# ============================================================================
# 3. DEFICIT BUDGET FOR TWO-SCALE SPECTRA
# ============================================================================

print("\n" + "="*80)
print("PART 2: Deficit budget scan over two-scale parameter space")
print("="*80)

# Fix λ and ⟨ω²⟩ to isolate shape effects
lam_ref = 1.721  # H3S-like
om2_ref = 149.26**2  # meV², H3S-like
sq_om2 = math.sqrt(om2_ref)

# Einstein reference at same (lambda, <om^2>)
print(f"\nReference: lam={lam_ref:.3f}, sqrt<om^2>={sq_om2:.2f} meV")
t_einstein_mu0 = Tc_P(np.array([sq_om2]), np.array([lam_ref * sq_om2 / 2.0]),
                       lam_ref, 0.0, 10.0, N=400, sqom2=sq_om2)
t_einstein_mu13 = Tc_P(np.array([sq_om2]), np.array([lam_ref * sq_om2 / 2.0]),
                        lam_ref, 0.13, 10.0, N=400, sqom2=sq_om2)
print(f"Einstein Tc(mu*=0): {t_einstein_mu0 * MEV_K:.1f} K")
print(f"Einstein Tc(mu*=0.13): {t_einstein_mu13 * MEV_K:.1f} K")
print(f"Einstein Coulomb deficit: {(1 - t_einstein_mu13/t_einstein_mu0)*100:.1f}%")

print(f"\n{'w':>5} {'freq_ratio':>11} {'lam_asym':>8} {'Tc(0) K':>9} "
      f"{'d_BRD%':>8} {'Tc(.13) K':>10} {'d_COU%':>8} {'fill':>6}")

# Scan parameter space
w_values = [0.5, 0.7, 0.8, 0.9, 0.95, 0.99]
freq_ratios = [1.5, 2.0, 5.0, 10.0, 100.0, 316.0]  # 316 ≈ √1e5

results = []
for w in w_values:
    for fr in freq_ratios:
        omegas, a2f = two_scale_spectrum(lam_ref, om2_ref, w, fr)

        # Compute Tc at mu*=0 and mu*=0.13
        try:
            tc_mu0 = Tc_P(omegas, a2f, lam_ref, 0.0, 10.0, N=400, sqom2=sq_om2)
            tc_mu13 = Tc_P(omegas, a2f, lam_ref, 0.13, 10.0, N=400, sqom2=sq_om2)

            if tc_mu0 and tc_mu13:
                # Deficits
                delta_brd = (1 - tc_mu0 / t_einstein_mu0) * 100  # broadening
                delta_cou = (1 - tc_mu13 / tc_mu0) * 100  # Coulomb

                # Fill factor
                ceiling = CINF * math.sqrt(lam_ref * om2_ref) * MEV_K
                fill = tc_mu13 * MEV_K / ceiling

                lam_asym = abs(2*w - 1)

                print(f"{w:>5.2f} {fr:>11.1f} {lam_asym:>8.3f} "
                      f"{tc_mu0*MEV_K:>9.1f} {delta_brd:>7.1f}% "
                      f"{tc_mu13*MEV_K:>10.1f} {delta_cou:>7.1f}% {fill:>6.3f}")

                results.append((w, fr, lam_asym, tc_mu0, tc_mu13, delta_brd, delta_cou, fill))
        except Exception as e:
            print(f"{w:>5.2f} {fr:>11.1f}  FAILED: {e}")

# ============================================================================
# 4. ANALYZE TRENDS
# ============================================================================

print("\n" + "="*80)
print("PART 3: Analysis of two-scale effects")
print("="*80)

if results:
    results = np.array(results)

    # Find maximum Coulomb deficit
    max_cou_idx = np.argmax(results[:, 6])
    print(f"\nMaximum Coulomb deficit:")
    print(f"  w={results[max_cou_idx, 0]:.2f}, freq_ratio={results[max_cou_idx, 1]:.1f}")
    print(f"  d_COU = {results[max_cou_idx, 6]:.1f}%")

    # Compare to Einstein
    ein_cou_def = (1 - t_einstein_mu13/t_einstein_mu0) * 100
    print(f"  Einstein d_COU = {ein_cou_def:.1f}%")
    print(f"  Two-scale enhancement: {results[max_cou_idx, 6] - ein_cou_def:+.1f} pp")

    # Find counter-example-like point (w=0.99, large freq_ratio)
    ce_like = results[(results[:, 0] == 0.99) & (results[:, 1] >= 100)]
    if len(ce_like) > 0:
        print(f"\nCounter-example-like point (w=0.99, freq_ratio>=100):")
        print(f"  Tc(mu*=0.13) = {ce_like[0, 4]*MEV_K:.1f} K")
        print(f"  d_COU = {ce_like[0, 6]:.1f}%")
        print(f"  fill factor = {ce_like[0, 7]:.3f}")

print("\n" + "="*80)
print("PART 4: Falsifiable statement")
print("="*80)
print("""
Hypothesis: Two-scale spectral structure INCREASES Coulomb deficit beyond
Einstein spectrum levels.

Quantitative prediction:
For fixed (lambda, <om^2>), a spectrum with:
  - lambda asymmetry > 0.8 (>90% in one mode)
  - Frequency separation > 1 decade
should exhibit d_COU at least 2-5 percentage points HIGHER than Einstein.

Test case: MgB2
  - Experimental: Tc = 39 K, lambda ~ 0.62, sqrt<om^2> ~ 60 meV
  - Two-band structure: sigma/pi with moderate asymmetry (~0.2) and separation (~0.3 decades)
  - Prediction: d_COU should be intermediate between Einstein and extreme two-scale
  - This can be falsified by computing Tc from actual MgB2 a2F with mu*=0.13
    and comparing to Einstein spectrum at same (lambda, <om^2>)

Alternative (stronger): If two-scale structure enhances deficit, then materials
closest to extremality failure boundary should systematically show two-peak a2F
with large lambda asymmetry. Can be falsified by examining MgB2 (safety ratio 1.8)
a2F shape versus LaH10 (safety ratio >9) shape.
""")

print("="*80)
print("Script complete. All results are [numerical-evidence] level.")
print("="*80)
