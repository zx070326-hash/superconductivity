#!/usr/bin/env python3
"""
Noise-weighted gap-ratio check for parameter separation.

Computes the 2x2 Jacobian for observables {Tc, 2D/kTc} and {Tc, alpha_isotope}
to determine whether gap ratio provides better lam/mu* discrimination than
isotope exponent under realistic measurement uncertainties.

Strictness: numerical-evidence (float64 Eliashberg solver)
"""
import os
import sys
import numpy as np

# Set thread count before importing numpy-dependent modules
os.environ['OMP_NUM_THREADS'] = '1'

# Add solver path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

# Anchor parameters (hardcoded per task spec)
LAM = 1.876
WE = 106.9  # meV
MU = 0.13
WCUT = 500.0  # meV (hardcoded, NOT 10*WE)
N_MATSUBARA = 1600
T_LOW_FACTOR = 0.03

# Perturbation sizes for central difference
DLAM = 0.15 * LAM
DMU = 0.05

def solve_at_point(lam, w_e, mu):
    """Solve for Tc and gap ratio at given parameters."""
    spec = Spectrum.einstein(lam=lam, w0=w_e)
    solver = Solver(spectrum=spec, mu_star=mu, w_cut=WCUT, n_matsubara=N_MATSUBARA)

    tc = solver.tc()
    t_low = T_LOW_FACTOR * tc

    wn, delta, z, conv = solver.solve_gap(t_low, mixing=0.25, max_iter=4000)

    if not conv:
        raise RuntimeError(f"Gap equation did not converge at lam={lam}, w_e={w_e}, mu={mu}")

    gap_ratio = 2.0 * delta[0] / (KB_MEV_PER_K * tc)

    return tc, gap_ratio

def central_diff_lam(w_e, mu):
    """Central difference d/dlam for Tc and gap ratio."""
    tc_p, gr_p = solve_at_point(LAM + DLAM, w_e, mu)
    tc_m, gr_m = solve_at_point(LAM - DLAM, w_e, mu)

    dtc_dlam = (tc_p - tc_m) / (2 * DLAM)
    dgr_dlam = (gr_p - gr_m) / (2 * DLAM)

    return dtc_dlam, dgr_dlam

def central_diff_mu(w_e):
    """Central difference d/dmu* for Tc and gap ratio."""
    tc_p, gr_p = solve_at_point(LAM, w_e, MU + DMU)
    tc_m, gr_m = solve_at_point(LAM, w_e, MU - DMU)

    dtc_dmu = (tc_p - tc_m) / (2 * DMU)
    dgr_dmu = (gr_p - gr_m) / (2 * DMU)

    return dtc_dmu, dgr_dmu

def central_diff_we(mu):
    """Central difference d/dw_E for Tc (for isotope exponent)."""
    dwe = 0.01 * WE  # 1% perturbation

    tc_p, _ = solve_at_point(LAM, WE + dwe, mu)
    tc_m, _ = solve_at_point(LAM, WE - dwe, mu)

    dtc_dwe = (tc_p - tc_m) / (2 * dwe)

    # alpha = +½ · dlnTc/dlnw_E
    tc_center, _ = solve_at_point(LAM, WE, mu)
    alpha = 0.5 * (dtc_dwe * WE / tc_center)

    return alpha

def confidence_ellipse_area(jac, sigma_vec):
    """
    Compute confidence region area for 2D observable vector.

    Area prop 1/sqrtdet(J^T Sigma^{-1} J) where Sigma = diag(σ_1^2, σ_2^2)

    For comparing two sets with different observables, the ratio
    Area_A / Area_B = sqrt[det(J_B^T Sigma_B^{-1} J_B) / det(J_A^T Sigma_A^{-1} J_A)]
    """
    # Sigma^{-1} = diag(1/σ_1^2, 1/σ_2^2)
    sigma_inv = np.diag(1.0 / sigma_vec**2)

    # Fisher information matrix: J^T Sigma^{-1} J
    fisher = jac.T @ sigma_inv @ jac

    det_fisher = np.linalg.det(fisher)

    if det_fisher <= 0:
        return np.inf  # Singular or indefinite

    # Area prop 1/sqrtdet(Fisher)
    return 1.0 / np.sqrt(det_fisher)

print("=" * 70)
print("VALIDATION GATES")
print("=" * 70)

# Gate 1: Reproduce R_ME, R_Tc, net gain at anchor
print("\n[Gate 1] Reproducing shape ratios at anchor point...")
tc_0, gr_0 = solve_at_point(LAM, WE, MU)
print(f"  Anchor: Tc = {tc_0:.4f} K, 2D/kTc = {gr_0:.6f}")

dtc_dlam, dgr_dlam = central_diff_lam(WE, MU)
dtc_dmu, dgr_dmu = central_diff_mu(WE)

print(f"\n  dTc/dlam = {dtc_dlam:.4f} K")
print(f"  dTc/dmu = {dtc_dmu:.4f} K")
print(f"  d(2D/kTc)/dlam = {dgr_dlam:.6f}")
print(f"  d(2D/kTc)/dmu = {dgr_dmu:.6f}")

# Shape ratios (noise-independent)
delta_tc_lam = dtc_dlam * DLAM
delta_tc_mu = dtc_dmu * DMU
delta_gr_lam = dgr_dlam * DLAM
delta_gr_mu = dgr_dmu * DMU

R_Tc = abs(delta_tc_lam / delta_tc_mu)
R_ME = abs(delta_gr_lam / delta_gr_mu)
net_gain = R_ME / R_Tc

print(f"\n  |DTc_lam / DTc_mu| = {R_Tc:.4f}")
print(f"  |D(2D/kTc)_lam / D(2D/kTc)_mu| = {R_ME:.4f}")
print(f"  Net gain = {net_gain:.4f}")

# Check against target: R_ME = 2.0080, R_Tc = 1.3710, gain = 1.4646
tol = 0.02  # ±2%
target_R_ME = 2.0080
target_R_Tc = 1.3710
target_gain = 1.4646

pass_R_ME = abs(R_ME - target_R_ME) / target_R_ME < tol
pass_R_Tc = abs(R_Tc - target_R_Tc) / target_R_Tc < tol
pass_gain = abs(net_gain - target_gain) / target_gain < tol

print(f"\n  Validation (±2%):")
print(f"    R_ME: {'PASS' if pass_R_ME else 'FAIL'} (target {target_R_ME}, got {R_ME:.4f})")
print(f"    R_Tc: {'PASS' if pass_R_Tc else 'FAIL'} (target {target_R_Tc}, got {R_Tc:.4f})")
print(f"    Gain: {'PASS' if pass_gain else 'FAIL'} (target {target_gain}, got {net_gain:.4f})")

if not (pass_R_ME and pass_R_Tc and pass_gain):
    print("\n  *** GATE 1 FAILED: Cannot proceed. ***")
    sys.exit(1)

print("\n  OK Gate 1 PASSED")

# Gate 2: alpha(mu*=0) = 0.5
print("\n[Gate 2] Isotope exponent at mu*=0 must equal 0.5...")
alpha_mu0 = central_diff_we(mu=0.0)
print(f"  alpha(mu*=0) = {alpha_mu0:.6f}")

pass_alpha = abs(alpha_mu0 - 0.5) < 0.01
print(f"  Validation: {'PASS' if pass_alpha else 'FAIL'} (target 0.5, tol ±0.01)")

if not pass_alpha:
    print("\n  *** GATE 2 FAILED: Isotope exponent implementation is wrong. ***")
    sys.exit(1)

print("  OK Gate 2 PASSED")

# Now compute isotope exponent at mu*=0.13
print("\n[Isotope exponent at anchor mu*=0.13]")
alpha_anchor = central_diff_we(mu=MU)
print(f"  alpha(mu*={MU}) = {alpha_anchor:.6f}")

print("\n" + "=" * 70)
print("JACOBIAN COMPUTATION")
print("=" * 70)

# Jacobian for Set A: {Tc, 2D/kTc}
J_A = np.array([
    [dtc_dlam, dtc_dmu],       # dTc/d(lam, mu*)
    [dgr_dlam, dgr_dmu]        # d(2D/kTc)/d(lam, mu*)
])

print("\n[Set A] Observables: {Tc, 2D/kTc}")
print(f"  Jacobian J_A = d(O_1,O_2)/d(lam,mu*):")
print(f"    [{J_A[0,0]:+.6f}  {J_A[0,1]:+.6f}]")
print(f"    [{J_A[1,0]:+.6f}  {J_A[1,1]:+.6f}]")
print(f"  det(J_A) = {np.linalg.det(J_A):.6e}")

# Jacobian for Set B: {Tc, alpha}
# dalpha/dlam: alpha only depends on lam,mu* through Tc's dependence, need to compute
print("\n[Set B] Observables: {Tc, alpha_isotope}")
print("  Computing dalpha/dlam and dalpha/dmu* via finite differences...")

# Perturb lam
alpha_p_lam, _ = solve_at_point(LAM + DLAM, WE, MU)
alpha_m_lam, _ = solve_at_point(LAM - DLAM, WE, MU)
# Recompute alpha at perturbed lam
alpha_p_lam_val = central_diff_we(mu=MU)  # Need to be careful here

# Actually, alpha depends on (lam, w_E, mu*), and we're varying lam,mu* at fixed w_E
# So dalpha/dlam at fixed w_E, fixed mu*
def alpha_at_params(lam, mu):
    """Compute alpha at given lam, mu*, fixed w_E."""
    return central_diff_we(mu=mu)  # This computes alpha from Tc(lam, w_E, mu*)

# Wait, central_diff_we uses global LAM. Need to be more careful.
def compute_alpha_at(lam, w_e, mu):
    """Compute isotope exponent at given parameters."""
    dwe = 0.01 * w_e
    tc_p, _ = solve_at_point(lam, w_e + dwe, mu)
    tc_m, _ = solve_at_point(lam, w_e - dwe, mu)
    dtc_dwe = (tc_p - tc_m) / (2 * dwe)
    tc_center, _ = solve_at_point(lam, w_e, mu)
    alpha = 0.5 * (dtc_dwe * w_e / tc_center)
    return alpha

# Now compute derivatives of alpha
alpha_anchor_check = compute_alpha_at(LAM, WE, MU)
print(f"  alpha at anchor (recomputed): {alpha_anchor_check:.6f}")

# dalpha/dlam
alpha_p_lam = compute_alpha_at(LAM + DLAM, WE, MU)
alpha_m_lam = compute_alpha_at(LAM - DLAM, WE, MU)
dalpha_dlam = (alpha_p_lam - alpha_m_lam) / (2 * DLAM)

# dalpha/dmu*
alpha_p_mu = compute_alpha_at(LAM, WE, MU + DMU)
alpha_m_mu = compute_alpha_at(LAM, WE, MU - DMU)
dalpha_dmu = (alpha_p_mu - alpha_m_mu) / (2 * DMU)

print(f"  dalpha/dlam = {dalpha_dlam:.6f}")
print(f"  dalpha/dmu* = {dalpha_dmu:.6f}")

J_B = np.array([
    [dtc_dlam, dtc_dmu],       # dTc/d(lam, mu*)
    [dalpha_dlam, dalpha_dmu]  # dalpha/d(lam, mu*)
])

print(f"\n  Jacobian J_B = d(Tc,alpha)/d(lam,mu*):")
print(f"    [{J_B[0,0]:+.6f}  {J_B[0,1]:+.6f}]")
print(f"    [{J_B[1,0]:+.6f}  {J_B[1,1]:+.6f}]")
print(f"  det(J_B) = {np.linalg.det(J_B):.6e}")

print("\n" + "=" * 70)
print("HEAD-TO-HEAD COMPARISON")
print("=" * 70)

# We need to compute area ratios under different error models
# The flip threshold is where Area_A / Area_B = 1

# For a 2x2 system, if observables are O = (O_1, O_2) with errors (σ_1, σ_2):
# Fisher = J^T Sigma^{-1} J where Sigma = diag(σ_1^2, σ_2^2)
# Area prop 1/sqrtdet(Fisher)

# For Set A: σ_A = (σ_Tc, σ_gr)
# For Set B: σ_B = (σ_Tc, σ_alpha)

# The ratio Area_A/Area_B depends on the ratio σ_gr/σ_Tc and σ_alpha/σ_Tc

# Let's parameterize: let σ_Tc = σ₀ (reference), σ_gr = r_gr · σ₀, σ_alpha = r_alpha · σ₀
# Then we can study how the area ratio depends on r_gr

# But actually, we want to find the threshold r_gr such that Area_A = Area_B

# Let's work in relative units: σ_Tc = 1 (arbitrary units)
# Then σ_gr = r_gr, σ_alpha = r_alpha

# For typical literature values:
# - Tc measured to ~1% or better (but with systematic ~10% for parameter inference)
# - 2D/kTc typically has scatter of several percent
# - alpha typically has scatter of order 0.05 in absolute terms

# Let's compute the flip threshold analytically

print("\n[Flip threshold analysis]")
print("  Finding σ(2D/kTc)/σ(Tc) threshold where Set A and Set B have equal area...")

# Set σ_Tc = 1 (relative units)
sigma_tc = 1.0

# For Set B, we need typical σ_alpha. But alpha is dimensionless ~0.3-0.5
# Typical scatter in alpha is ~0.05 absolute, or ~10-15% relative
# Let's use σ_alpha as given by typical literature scatter

# But we need to be careful: the area ratio is
# Area_A/Area_B = sqrt[det(J_B^T Sigma_B^{-1} J_B) / det(J_A^T Sigma_A^{-1} J_A)]

# Let σ_A = (σ_Tc, σ_gr), σ_B = (σ_Tc, σ_alpha)
# Both share the same σ_Tc, so the ratio simplifies

# Let r = σ_gr/σ_Tc be the unknown threshold
# We want to find r such that Area_A(r) = Area_B(r_alpha) for some assumed σ_alpha

# Actually, let's compute it parametrically and find the crossing point

# For Set A with σ_A = (1, r):
def area_ratio_vs_r(r, r_alpha=None):
    """
    Compute Area_A/Area_B as function of r = σ_gr/σ_Tc.

    If r_alpha is None, we compute the threshold assuming typical σ_alpha.
    """
    sigma_A = np.array([sigma_tc, r * sigma_tc])

    # For Set B, we need σ_alpha/σ_Tc
    # Typical absolute σ_alpha ~ 0.05, and alpha ~ 0.4, so relative ~12%
    # But Tc ~ 50-200 K with σ_Tc ~ 1-2 K gives relative ~1-2%
    # So σ_alpha/σ_Tc in absolute units is harder to define

    # Let's instead think: if σ_Tc/Tc ~ 0.01 (1%) and σ_alpha/alpha ~ 0.12 (12%)
    # with alpha ~ 0.4, then σ_alpha ~ 0.05 (absolute)
    # This is independent of Tc in absolute units

    # But for the Fisher information, what matters is the relative error
    # in the observable space

    # Let me reconsider: the Jacobian relates parameter perturbations to
    # observable perturbations: δO = J δtheta
    # The Fisher information is F = J^T Sigma^{-1} J
    # where Sigma is the covariance of observables

    # For Tc in Kelvin with σ_Tc ~ 2 K (H3S example)
    # For alpha dimensionless with σ_alpha ~ 0.05 (typical scatter)
    # These are in different units!

    # The correct way: assume σ_Tc = 2 K, σ_gr = ? (dimensionless), σ_alpha = 0.05
    # But this mixes units...

    # Actually, the Fisher info doesn't care about units - it's asking
    # "how much information about theta do we get from measuring O?"
    # The key is the SNR: J δtheta / σ_O

    # Let me restart with dimensional analysis:
    # dTc/dlam has units [K]
    # dalpha/dlam has units [dimensionless]
    # σ_Tc has units [K]
    # σ_alpha is dimensionless

    # The Fisher matrix element is (dO/dtheta)^2 / σ_O^2
    # For Tc: (dTc/dlam)^2 / σ_Tc^2  has units [K^2/K^2] = dimensionless OK
    # For alpha: (dalpha/dlam)^2 / σ_alpha^2  has units [1/1] = dimensionless OK

    # So we can work with actual units. Let me assume:
    # σ_Tc = 2 K (from Drozdov, but actually represents systematic uncertainty)
    # σ_alpha = 0.05 (typical scatter in alpha measurements)

    if r_alpha is None:
        # Assume typical σ_alpha = 0.05 (absolute, dimensionless)
        # and σ_Tc = 2 K
        sigma_tc_actual = 2.0  # K
        sigma_alpha_actual = 0.05  # dimensionless
        sigma_B = np.array([sigma_tc_actual, sigma_alpha_actual])
    else:
        sigma_B = np.array([sigma_tc, r_alpha * sigma_tc])

    area_A = confidence_ellipse_area(J_A, sigma_A)
    area_B = confidence_ellipse_area(J_B, sigma_B)

    return area_A / area_B

# But we still have a units problem. Let me think differently.
# The actual observables are:
# - Tc in Kelvin: ~50-200 K, measured with σ ~ 1-2 K (statistical) or ~10% (systematic)
# - 2D/kTc dimensionless: ~3-4, measured with σ ~ 0.1-0.3 (several percent)
# - alpha dimensionless: ~0.3-0.5, measured with σ ~ 0.05 (absolute)

# In our calculation:
# - Tc ~ 50 K (from anchor)
# - 2D/kTc ~ 3.5 (typical)
# - alpha ~ 0.4 (from our calc)

print(f"\n  Reference values at anchor:")
print(f"    Tc = {tc_0:.4f} K")
print(f"    2D/kTc = {gr_0:.6f} (dimensionless)")
print(f"    alpha = {alpha_anchor:.6f} (dimensionless)")

# Now let's set up two error models

print("\n" + "=" * 70)
print("ERROR MODEL 1: Statistical (measurement precision)")
print("=" * 70)

# Typical statistical uncertainties:
# Tc: ~1% of value (for clean samples with good thermometry)
# 2D/kTc: ~3-5% (from tunneling spectroscopy, see Semenok 2020)
# alpha: ~10% relative (~0.05 absolute for alpha~0.5)

sigma_Tc_stat = 0.01 * tc_0  # 1% statistical
sigma_gr_stat = 0.04 * gr_0  # 4% (mid-range of 3-5%)
sigma_alpha_stat = 0.05  # 0.05 absolute

sigma_A_stat = np.array([sigma_Tc_stat, sigma_gr_stat])
sigma_B_stat = np.array([sigma_Tc_stat, sigma_alpha_stat])

area_A_stat = confidence_ellipse_area(J_A, sigma_A_stat)
area_B_stat = confidence_ellipse_area(J_B, sigma_B_stat)
ratio_stat = area_A_stat / area_B_stat

print(f"\n  Assumed statistical uncertainties:")
print(f"    σ(Tc) = {sigma_Tc_stat:.4f} K  ({100*sigma_Tc_stat/tc_0:.2f}%)")
print(f"    σ(2D/kTc) = {sigma_gr_stat:.6f}  ({100*sigma_gr_stat/gr_0:.2f}%)")
print(f"    σ(alpha) = {sigma_alpha_stat:.6f}  ({100*sigma_alpha_stat/alpha_anchor:.2f}%)")

print(f"\n  Confidence region areas (relative units):")
print(f"    Area(Set A) = {area_A_stat:.6e}")
print(f"    Area(Set B) = {area_B_stat:.6e}")
print(f"    Ratio A/B = {ratio_stat:.4f}")

if ratio_stat < 1.0:
    verdict_stat = "Set A WINS (gap ratio provides better discrimination)"
    improvement_stat = 1.0 / ratio_stat
    print(f"\n  -> {verdict_stat}")
    print(f"  -> Improvement factor: {improvement_stat:.2f}×")
else:
    verdict_stat = "Set B WINS (isotope exponent provides better discrimination)"
    degradation_stat = ratio_stat
    print(f"\n  -> {verdict_stat}")
    print(f"  -> Set A is worse by factor: {degradation_stat:.2f}×")

print("\n" + "=" * 70)
print("ERROR MODEL 2: Total (including systematics)")
print("=" * 70)

# Total uncertainties including sample-to-sample variation and systematics:
# Tc: ~10% (from pressure mismatch, annealing state, etc.) [H3S example: 184 vs 203 K]
# 2D/kTc: ~5-10% (sample quality, gap extraction method)
# alpha: ~15-20% relative (sample variation, mass-rescaling issues)

sigma_Tc_total = 0.10 * tc_0  # 10% total
sigma_gr_total = 0.075 * gr_0  # 7.5% (mid-range)
sigma_alpha_total = 0.15 * alpha_anchor  # 15% relative

sigma_A_total = np.array([sigma_Tc_total, sigma_gr_total])
sigma_B_total = np.array([sigma_Tc_total, sigma_alpha_total])

area_A_total = confidence_ellipse_area(J_A, sigma_A_total)
area_B_total = confidence_ellipse_area(J_B, sigma_B_total)
ratio_total = area_A_total / area_B_total

print(f"\n  Assumed total uncertainties:")
print(f"    σ(Tc) = {sigma_Tc_total:.4f} K  ({100*sigma_Tc_total/tc_0:.2f}%)")
print(f"    σ(2D/kTc) = {sigma_gr_total:.6f}  ({100*sigma_gr_total/gr_0:.2f}%)")
print(f"    σ(alpha) = {sigma_alpha_total:.6f}  ({100*sigma_alpha_total/alpha_anchor:.2f}%)")

print(f"\n  Confidence region areas (relative units):")
print(f"    Area(Set A) = {area_A_total:.6e}")
print(f"    Area(Set B) = {area_B_total:.6e}")
print(f"    Ratio A/B = {ratio_total:.4f}")

if ratio_total < 1.0:
    verdict_total = "Set A WINS (gap ratio provides better discrimination)"
    improvement_total = 1.0 / ratio_total
    print(f"\n  -> {verdict_total}")
    print(f"  -> Improvement factor: {improvement_total:.2f}×")
else:
    verdict_total = "Set B WINS (isotope exponent provides better discrimination)"
    degradation_total = ratio_total
    print(f"\n  -> {verdict_total}")
    print(f"  -> Set A is worse by factor: {degradation_total:.2f}×")

print("\n" + "=" * 70)
print("FLIP THRESHOLD ANALYSIS")
print("=" * 70)

# Find the threshold r* = σ(2D/kTc)/σ(Tc) where Area_A = Area_B
# This means Area_A/Area_B = 1

# For a given σ(Tc) and σ(alpha), we want to find σ(2D/kTc) such that areas are equal

# Using the same σ(Tc) for both sets:
# Area_A prop 1/sqrtdet(J_A^T diag(1/σ_Tc^2, 1/σ_gr^2) J_A)
# Area_B prop 1/sqrtdet(J_B^T diag(1/σ_Tc^2, 1/σ_alpha^2) J_B)

# Setting equal and solving for σ_gr:

def find_flip_threshold(sigma_tc, sigma_alpha):
    """
    Find σ_gr such that Area_A = Area_B.
    Returns σ_gr and the ratio σ_gr/σ_tc.
    """
    # We want det(Fisher_A) = det(Fisher_B)
    # Fisher = J^T Sigma^{-1} J

    # For Set B (fixed):
    sigma_B = np.array([sigma_tc, sigma_alpha])
    sigma_inv_B = np.diag(1.0 / sigma_B**2)
    fisher_B = J_B.T @ sigma_inv_B @ J_B
    det_B = np.linalg.det(fisher_B)

    # For Set A, we have σ_A = (σ_tc, σ_gr)
    # det(Fisher_A) = det(J_A^T diag(1/σ_tc^2, 1/σ_gr^2) J_A)
    # We want det(Fisher_A) = det(Fisher_B)

    # Let's parameterize σ_gr and search for the zero crossing
    from scipy.optimize import brentq

    def det_diff(sigma_gr):
        sigma_A = np.array([sigma_tc, sigma_gr])
        sigma_inv_A = np.diag(1.0 / sigma_A**2)
        fisher_A = J_A.T @ sigma_inv_A @ J_A
        det_A = np.linalg.det(fisher_A)
        return det_A - det_B

    # Search over a reasonable range
    # σ_gr should be positive and probably not more than 10× σ_tc
    try:
        sigma_gr_threshold = brentq(det_diff, 1e-6 * sigma_tc, 100 * sigma_tc)
        ratio_threshold = sigma_gr_threshold / sigma_tc
        return sigma_gr_threshold, ratio_threshold
    except ValueError:
        # No zero crossing in range
        return None, None

# For Model 1 (statistical)
sigma_gr_thresh_stat, ratio_thresh_stat = find_flip_threshold(sigma_Tc_stat, sigma_alpha_stat)

print(f"\n[Model 1: Statistical errors]")
print(f"  Given: σ(Tc) = {sigma_Tc_stat:.4f} K, σ(alpha) = {sigma_alpha_stat:.6f}")
if sigma_gr_thresh_stat is not None:
    print(f"  Flip threshold: σ(2D/kTc) = {sigma_gr_thresh_stat:.6f}")
    print(f"  Or equivalently: σ(2D/kTc)/σ(Tc) = {ratio_thresh_stat:.4f}")
    print(f"  Current ratio: σ(2D/kTc)/σ(Tc) = {sigma_gr_stat/sigma_Tc_stat:.4f}")

    if sigma_gr_stat < sigma_gr_thresh_stat:
        print(f"  -> Below threshold: Set A wins (gap ratio is better)")
    else:
        print(f"  -> Above threshold: Set B wins (isotope exponent is better)")
else:
    print(f"  No flip threshold found (one set always dominates)")

# For Model 2 (total)
sigma_gr_thresh_total, ratio_thresh_total = find_flip_threshold(sigma_Tc_total, sigma_alpha_total)

print(f"\n[Model 2: Total errors]")
print(f"  Given: σ(Tc) = {sigma_Tc_total:.4f} K, σ(alpha) = {sigma_alpha_total:.6f}")
if sigma_gr_thresh_total is not None:
    print(f"  Flip threshold: σ(2D/kTc) = {sigma_gr_thresh_total:.6f}")
    print(f"  Or equivalently: σ(2D/kTc)/σ(Tc) = {ratio_thresh_total:.4f}")
    print(f"  Current ratio: σ(2D/kTc)/σ(Tc) = {sigma_gr_total/sigma_Tc_total:.4f}")

    if sigma_gr_total < sigma_gr_thresh_total:
        print(f"  -> Below threshold: Set A wins (gap ratio is better)")
    else:
        print(f"  -> Above threshold: Set B wins (isotope exponent is better)")
else:
    print(f"  No flip threshold found (one set always dominates)")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\nShape-only analysis (noise-independent):")
print(f"  |DTc_lam / DTc_mu| = {R_Tc:.4f}")
print(f"  |D(2D/kTc)_lam / D(2D/kTc)_mu| = {R_ME:.4f}")
print(f"  Net gain factor = {net_gain:.4f}")

print(f"\nNoise-weighted comparison:")
print(f"  Model 1 (statistical): Area_A/Area_B = {ratio_stat:.4f} -> {verdict_stat}")
print(f"  Model 2 (total):       Area_A/Area_B = {ratio_total:.4f} -> {verdict_total}")

print(f"\nFlip thresholds σ(2D/kTc)/σ(Tc):")
if ratio_thresh_stat is not None:
    print(f"  Model 1: {ratio_thresh_stat:.4f}")
else:
    print(f"  Model 1: (no crossing)")
if ratio_thresh_total is not None:
    print(f"  Model 2: {ratio_thresh_total:.4f}")
else:
    print(f"  Model 2: (no crossing)")

print(f"\nStrictness level: numerical-evidence (float64 Eliashberg)")
print("\n" + "=" * 70)
