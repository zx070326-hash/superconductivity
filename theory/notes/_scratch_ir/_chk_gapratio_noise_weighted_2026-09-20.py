#!/usr/bin/env python3
"""
Noise-weighted gap-ratio check for parameter separation.

Optimized version with fewer Eliashberg solves.
Strictness: numerical-evidence (float64 Eliashberg solver)
"""
import os
import sys
import numpy as np
from scipy.optimize import brentq

os.environ['OMP_NUM_THREADS'] = '1'

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

# Anchor parameters
LAM = 1.876
WE = 106.9  # meV
MU = 0.13
WCUT = 500.0  # meV
N_MATSUBARA = 1600
T_LOW_FACTOR = 0.03

DLAM = 0.15 * LAM
DMU = 0.05

def solve_at_point(lam, w_e, mu, verbose=False):
    """Solve for Tc and gap ratio at given parameters."""
    if verbose:
        print(f"    Solving at lam={lam:.4f}, w_e={w_e:.4f}, mu={mu:.4f}...", flush=True)

    spec = Spectrum.einstein(lam=lam, w0=w_e)
    solver = Solver(spectrum=spec, mu_star=mu, w_cut=WCUT, n_matsubara=N_MATSUBARA)

    tc = solver.tc()
    t_low = T_LOW_FACTOR * tc

    wn, delta, z, conv = solver.solve_gap(t_low, mixing=0.25, max_iter=4000)

    if not conv:
        raise RuntimeError(f"Gap equation did not converge at lam={lam}, w_e={w_e}, mu={mu}")

    gap_ratio = 2.0 * delta[0] / (KB_MEV_PER_K * tc)

    return tc, gap_ratio

print("=" * 70)
print("VALIDATION GATES")
print("=" * 70)

# Gate 1: Reproduce R_ME, R_Tc, net gain
print("\n[Gate 1] Reproducing shape ratios at anchor point...")
print("  (6 Eliashberg solves total)")

tc_0, gr_0 = solve_at_point(LAM, WE, MU, verbose=True)
print(f"  Anchor: Tc = {tc_0:.4f} K, 2D/kTc = {gr_0:.6f}")

# d/dlam derivatives
print("\n  Computing d/dlam derivatives...")
tc_p, gr_p = solve_at_point(LAM + DLAM, WE, MU, verbose=True)
tc_m, gr_m = solve_at_point(LAM - DLAM, WE, MU, verbose=True)

dtc_dlam = (tc_p - tc_m) / (2 * DLAM)
dgr_dlam = (gr_p - gr_m) / (2 * DLAM)

# d/dmu derivatives
print("\n  Computing d/dmu derivatives...")
tc_p, gr_p = solve_at_point(LAM, WE, MU + DMU, verbose=True)
tc_m, gr_m = solve_at_point(LAM, WE, MU - DMU, verbose=True)

dtc_dmu = (tc_p - tc_m) / (2 * DMU)
dgr_dmu = (gr_p - gr_m) / (2 * DMU)

print(f"\n  dTc/dlam = {dtc_dlam:.4f} K")
print(f"  dTc/dmu = {dtc_dmu:.4f} K")
print(f"  d(2D/kTc)/dlam = {dgr_dlam:.6f}")
print(f"  d(2D/kTc)/dmu = {dgr_dmu:.6f}")

# Shape ratios
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

# Validation
tol = 0.02
target_R_ME = 2.0080
target_R_Tc = 1.3710
target_gain = 1.4646

pass_R_ME = abs(R_ME - target_R_ME) / target_R_ME < tol
pass_R_Tc = abs(R_Tc - target_R_Tc) / target_R_Tc < tol
pass_gain = abs(net_gain - target_gain) / target_gain < tol

print(f"\n  Validation (+/-2%):")
print(f"    R_ME: {'PASS' if pass_R_ME else 'FAIL'} (target {target_R_ME}, got {R_ME:.4f})")
print(f"    R_Tc: {'PASS' if pass_R_Tc else 'FAIL'} (target {target_R_Tc}, got {R_Tc:.4f})")
print(f"    Gain: {'PASS' if pass_gain else 'FAIL'} (target {target_gain}, got {net_gain:.4f})")

if not (pass_R_ME and pass_R_Tc and pass_gain):
    print("\n  *** GATE 1 FAILED ***")
    sys.exit(1)

print("\n  OK Gate 1 PASSED")

# Gate 2: alpha(mu*=0) = 0.5
print("\n[Gate 2] Isotope exponent at mu*=0 must equal 0.5...")
print("  (3 Eliashberg solves)")

dwe = 0.01 * WE
tc_p, _ = solve_at_point(LAM, WE + dwe, 0.0, verbose=True)
tc_m, _ = solve_at_point(LAM, WE - dwe, 0.0, verbose=True)
tc_c, _ = solve_at_point(LAM, WE, 0.0, verbose=True)

dtc_dwe = (tc_p - tc_m) / (2 * dwe)
alpha_mu0 = 0.5 * (dtc_dwe * WE / tc_c)

print(f"  alpha(mu*=0) = {alpha_mu0:.6f}")

pass_alpha = abs(alpha_mu0 - 0.5) < 0.01
print(f"  Validation: {'PASS' if pass_alpha else 'FAIL'} (target 0.5, tol +/-0.01)")

if not pass_alpha:
    print("\n  *** GATE 2 FAILED ***")
    sys.exit(1)

print("  OK Gate 2 PASSED")

# Compute alpha at anchor
print("\n[Isotope exponent at anchor mu*=0.13]")
print("  (2 additional Eliashberg solves, reusing anchor)")

tc_p_we, _ = solve_at_point(LAM, WE + dwe, MU, verbose=True)
tc_m_we, _ = solve_at_point(LAM, WE - dwe, MU, verbose=True)

dtc_dwe_anchor = (tc_p_we - tc_m_we) / (2 * dwe)
alpha_anchor = 0.5 * (dtc_dwe_anchor * WE / tc_0)

print(f"  alpha(mu*={MU}) = {alpha_anchor:.6f}")

print("\n" + "=" * 70)
print("JACOBIAN COMPUTATION")
print("=" * 70)

# Jacobian for Set A
J_A = np.array([
    [dtc_dlam, dtc_dmu],
    [dgr_dlam, dgr_dmu]
])

print("\n[Set A] Observables: {Tc, 2D/kTc}")
print(f"  Jacobian J_A:")
print(f"    [{J_A[0,0]:+.6f}  {J_A[0,1]:+.6f}]")
print(f"    [{J_A[1,0]:+.6f}  {J_A[1,1]:+.6f}]")
print(f"  det(J_A) = {np.linalg.det(J_A):.6e}")

# Jacobian for Set B - need dalpha/dlam and dalpha/dmu
print("\n[Set B] Observables: {Tc, alpha_isotope}")
print("  Computing dalpha/dlam and dalpha/dmu...")
print("  (4 additional Eliashberg solves)")

# dalpha/dlam
def compute_alpha_at(lam, w_e, mu, verbose=False):
    dwe_local = 0.01 * w_e
    tc_p, _ = solve_at_point(lam, w_e + dwe_local, mu, verbose=verbose)
    tc_m, _ = solve_at_point(lam, w_e - dwe_local, mu, verbose=verbose)
    tc_c, _ = solve_at_point(lam, w_e, mu, verbose=verbose)
    dtc_dwe_local = (tc_p - tc_m) / (2 * dwe_local)
    return 0.5 * (dtc_dwe_local * w_e / tc_c)

print("\n  For dalpha/dlam:")
alpha_p_lam = compute_alpha_at(LAM + DLAM, WE, MU, verbose=True)
alpha_m_lam = compute_alpha_at(LAM - DLAM, WE, MU, verbose=True)
dalpha_dlam = (alpha_p_lam - alpha_m_lam) / (2 * DLAM)

print("\n  For dalpha/dmu:")
alpha_p_mu = compute_alpha_at(LAM, WE, MU + DMU, verbose=True)
alpha_m_mu = compute_alpha_at(LAM, WE, MU - DMU, verbose=True)
dalpha_dmu = (alpha_p_mu - alpha_m_mu) / (2 * DMU)

print(f"\n  dalpha/dlam = {dalpha_dlam:.6f}")
print(f"  dalpha/dmu = {dalpha_dmu:.6f}")

J_B = np.array([
    [dtc_dlam, dtc_dmu],
    [dalpha_dlam, dalpha_dmu]
])

print(f"\n  Jacobian J_B:")
print(f"    [{J_B[0,0]:+.6f}  {J_B[0,1]:+.6f}]")
print(f"    [{J_B[1,0]:+.6f}  {J_B[1,1]:+.6f}]")
print(f"  det(J_B) = {np.linalg.det(J_B):.6e}")

print(f"\n  Total Eliashberg solves: 6 + 3 + 2 + 6 + 6 = 23")

print("\n" + "=" * 70)
print("HEAD-TO-HEAD COMPARISON")
print("=" * 70)

def confidence_ellipse_area(jac, sigma_vec):
    """Compute confidence region area (relative units)."""
    sigma_inv = np.diag(1.0 / sigma_vec**2)
    fisher = jac.T @ sigma_inv @ jac
    det_fisher = np.linalg.det(fisher)
    if det_fisher <= 0:
        return np.inf
    return 1.0 / np.sqrt(det_fisher)

# Error Model 1: Statistical
print("\n" + "=" * 70)
print("ERROR MODEL 1: Statistical (measurement precision)")
print("=" * 70)

sigma_Tc_stat = 0.01 * tc_0  # 1%
sigma_gr_stat = 0.04 * gr_0  # 4%
sigma_alpha_stat = 0.05      # absolute

sigma_A_stat = np.array([sigma_Tc_stat, sigma_gr_stat])
sigma_B_stat = np.array([sigma_Tc_stat, sigma_alpha_stat])

area_A_stat = confidence_ellipse_area(J_A, sigma_A_stat)
area_B_stat = confidence_ellipse_area(J_B, sigma_B_stat)
ratio_stat = area_A_stat / area_B_stat

print(f"\n  Assumed uncertainties:")
print(f"    sigma(Tc) = {sigma_Tc_stat:.4f} K  ({100*sigma_Tc_stat/tc_0:.2f}%)")
print(f"    sigma(2D/kTc) = {sigma_gr_stat:.6f}  ({100*sigma_gr_stat/gr_0:.2f}%)")
print(f"    sigma(alpha) = {sigma_alpha_stat:.6f}  ({100*sigma_alpha_stat/alpha_anchor:.2f}%)")

print(f"\n  Confidence region areas (relative):")
print(f"    Area(Set A) = {area_A_stat:.6e}")
print(f"    Area(Set B) = {area_B_stat:.6e}")
print(f"    Ratio A/B = {ratio_stat:.4f}")

if ratio_stat < 1.0:
    print(f"\n  -> Set A WINS (gap ratio better)")
    print(f"  -> Improvement: {1.0/ratio_stat:.2f}x")
else:
    print(f"\n  -> Set B WINS (isotope exponent better)")
    print(f"  -> Set A worse by: {ratio_stat:.2f}x")

# Error Model 2: Total
print("\n" + "=" * 70)
print("ERROR MODEL 2: Total (including systematics)")
print("=" * 70)

sigma_Tc_total = 0.10 * tc_0    # 10%
sigma_gr_total = 0.075 * gr_0   # 7.5%
sigma_alpha_total = 0.15 * alpha_anchor  # 15% relative

sigma_A_total = np.array([sigma_Tc_total, sigma_gr_total])
sigma_B_total = np.array([sigma_Tc_total, sigma_alpha_total])

area_A_total = confidence_ellipse_area(J_A, sigma_A_total)
area_B_total = confidence_ellipse_area(J_B, sigma_B_total)
ratio_total = area_A_total / area_B_total

print(f"\n  Assumed uncertainties:")
print(f"    sigma(Tc) = {sigma_Tc_total:.4f} K  ({100*sigma_Tc_total/tc_0:.2f}%)")
print(f"    sigma(2D/kTc) = {sigma_gr_total:.6f}  ({100*sigma_gr_total/gr_0:.2f}%)")
print(f"    sigma(alpha) = {sigma_alpha_total:.6f}  ({100*sigma_alpha_total/alpha_anchor:.2f}%)")

print(f"\n  Confidence region areas (relative):")
print(f"    Area(Set A) = {area_A_total:.6e}")
print(f"    Area(Set B) = {area_B_total:.6e}")
print(f"    Ratio A/B = {ratio_total:.4f}")

if ratio_total < 1.0:
    print(f"\n  -> Set A WINS (gap ratio better)")
    print(f"  -> Improvement: {1.0/ratio_total:.2f}x")
else:
    print(f"\n  -> Set B WINS (isotope exponent better)")
    print(f"  -> Set A worse by: {ratio_total:.2f}x")

print("\n" + "=" * 70)
print("FLIP THRESHOLD ANALYSIS")
print("=" * 70)

def find_flip_threshold(sigma_tc, sigma_alpha):
    """Find sigma_gr such that Area_A = Area_B."""
    sigma_B = np.array([sigma_tc, sigma_alpha])
    sigma_inv_B = np.diag(1.0 / sigma_B**2)
    fisher_B = J_B.T @ sigma_inv_B @ J_B
    det_B = np.linalg.det(fisher_B)

    def det_diff(sigma_gr):
        sigma_A = np.array([sigma_tc, sigma_gr])
        sigma_inv_A = np.diag(1.0 / sigma_A**2)
        fisher_A = J_A.T @ sigma_inv_A @ J_A
        return np.linalg.det(fisher_A) - det_B

    try:
        sigma_gr_threshold = brentq(det_diff, 1e-6 * sigma_tc, 100 * sigma_tc)
        return sigma_gr_threshold, sigma_gr_threshold / sigma_tc
    except ValueError:
        return None, None

# Model 1 threshold
sigma_gr_thresh_stat, ratio_thresh_stat = find_flip_threshold(sigma_Tc_stat, sigma_alpha_stat)

print(f"\n[Model 1: Statistical]")
print(f"  Given: sigma(Tc) = {sigma_Tc_stat:.4f} K, sigma(alpha) = {sigma_alpha_stat:.6f}")
if sigma_gr_thresh_stat is not None:
    print(f"  Flip threshold: sigma(2D/kTc) = {sigma_gr_thresh_stat:.6f}")
    print(f"  Or: sigma(2D/kTc)/sigma(Tc) = {ratio_thresh_stat:.4f}")
    print(f"  Current: sigma(2D/kTc)/sigma(Tc) = {sigma_gr_stat/sigma_Tc_stat:.4f}")

    if sigma_gr_stat < sigma_gr_thresh_stat:
        print(f"  -> Below threshold: Set A wins")
    else:
        print(f"  -> Above threshold: Set B wins")
else:
    print(f"  No flip threshold found")

# Model 2 threshold
sigma_gr_thresh_total, ratio_thresh_total = find_flip_threshold(sigma_Tc_total, sigma_alpha_total)

print(f"\n[Model 2: Total]")
print(f"  Given: sigma(Tc) = {sigma_Tc_total:.4f} K, sigma(alpha) = {sigma_alpha_total:.6f}")
if sigma_gr_thresh_total is not None:
    print(f"  Flip threshold: sigma(2D/kTc) = {sigma_gr_thresh_total:.6f}")
    print(f"  Or: sigma(2D/kTc)/sigma(Tc) = {ratio_thresh_total:.4f}")
    print(f"  Current: sigma(2D/kTc)/sigma(Tc) = {sigma_gr_total/sigma_Tc_total:.4f}")

    if sigma_gr_total < sigma_gr_thresh_total:
        print(f"  -> Below threshold: Set A wins")
    else:
        print(f"  -> Above threshold: Set B wins")
else:
    print(f"  No flip threshold found")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\nShape-only (noise-independent):")
print(f"  |DTc_lam / DTc_mu| = {R_Tc:.4f}")
print(f"  |D(2D/kTc)_lam / D(2D/kTc)_mu| = {R_ME:.4f}")
print(f"  Net gain = {net_gain:.4f}")

print(f"\nNoise-weighted:")
print(f"  Model 1: Area_A/Area_B = {ratio_stat:.4f}")
print(f"  Model 2: Area_A/Area_B = {ratio_total:.4f}")

print(f"\nFlip thresholds sigma(2D/kTc)/sigma(Tc):")
if ratio_thresh_stat is not None:
    print(f"  Model 1: {ratio_thresh_stat:.4f}")
else:
    print(f"  Model 1: (no crossing)")
if ratio_thresh_total is not None:
    print(f"  Model 2: {ratio_thresh_total:.4f}")
else:
    print(f"  Model 2: (no crossing)")

print(f"\nStrictness: numerical-evidence (float64 Eliashberg)")
print("=" * 70)
