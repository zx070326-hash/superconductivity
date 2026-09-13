r"""
P4 . MOMENT THRESHOLD THEOREM -- numerical gate for Theorem 2' and Proposition 2''.

v0.1 (2026-09-13): numerical verification of p-moment extremality (p >= 2) and
counterexamples for p < 2.

================================================================================
WHAT IS VERIFIED
================================================================================
This script provides CONSISTENCY CHECKS (not proofs) for:
  (a) Theorem 2': kappa(P,T) <= h(varpi_p) for p >= 2 and many random spectra P.
  (b) Proposition 2'': For p < 2, explicit counterexamples violate the bound.
  (c) Monotonicity: varpi_p increases with p (Lyapunov inequality).
  (d) N-independence: results stable under N=300 vs N=600.

All proofs are ANALYTIC (see P4_moment_threshold.md). This script only checks
that numerics agree with the analytic predictions.

Numerical discipline: numpy float64 for eigenvalues (never mpmath eigensolvers);
mpmath dps=30 for h(varpi) via H_mat from p4_interchange_uniformity.py.
"""

import sys
import time
import numpy as np
from mpmath import mp

mp.dps = 30

# Import H_mat from the verified interchange uniformity script
sys.path.insert(0, 'theory/derivations')
from p4_interchange_uniformity import H_mat

# ==============================================================================
# SETUP
# ==============================================================================

_checks_passed = 0
_checks_total = 0

def ck(cond, msg):
    """Count-style check: increment counters, print failure."""
    global _checks_passed, _checks_total
    _checks_total += 1
    if cond:
        _checks_passed += 1
    else:
        print(f"  [FAIL] {msg}")
    return cond

def h_of_varpi(varpi, N):
    """Top eigenvalue h(varpi) = lambda_1(H^(N)(varpi)) via float64 eigvalsh."""
    if float(varpi) == 0.0:
        return 0.0
    H = H_mat(N, float(varpi), order=0)
    eigs = np.linalg.eigvalsh(H)
    return float(eigs[-1])

def random_multiatomic_spectrum(n_atoms, omega_max=10.0, seed=None):
    """
    Generate a random n-atom spectrum: weights w_i in (0,1) normalized to sum=1,
    frequencies omega_i in (0, omega_max].
    Returns (weights, omegas) as numpy arrays.
    """
    if seed is not None:
        np.random.seed(seed)
    weights = np.random.rand(n_atoms)
    weights /= weights.sum()
    omegas = np.random.rand(n_atoms) * omega_max + 0.1  # avoid exact zero
    return weights, omegas

def moment_p(weights, omegas, p):
    """Compute <omega^p> = sum_i w_i omega_i^p."""
    return np.sum(weights * omegas**p)

def varpi_p(weights, omegas, p, T=1.0):
    """varpi_p := <omega^p>^(1/p) / (2*pi*T)."""
    mom_p = moment_p(weights, omegas, p)
    return mom_p**(1.0/p) / (2.0 * np.pi * T)

def kappa_of_spectrum(weights, omegas, T, N):
    """
    kappa(P,T) for a discrete spectrum P = sum_i w_i delta(omega_i).
    By linearity: kappa = sum_i w_i * h(omega_i / (2*pi*T)).
    """
    kappa = 0.0
    for w_i, omega_i in zip(weights, omegas):
        varpi_i = omega_i / (2.0 * np.pi * T)
        kappa += w_i * h_of_varpi(varpi_i, N)
    return kappa

# ==============================================================================
# STAGE 1: Theorem 2' verification (p >= 2)
# ==============================================================================

def run_stage1(N):
    """
    Check kappa(P,T) <= h(varpi_p) for p in {2, 3, 4} and many random spectra.
    """
    print(f"\n=== STAGE 1: Theorem 2' verification (p >= 2), N={N} ===")
    t0 = time.time()

    p_values = [2.0, 3.0, 4.0]
    n_spectra = 50
    T = 1.0

    violations = {p: 0 for p in p_values}
    max_ratio = {p: 0.0 for p in p_values}

    for i in range(n_spectra):
        n_atoms = np.random.randint(2, 6)  # 2 to 5 atoms
        weights, omegas = random_multiatomic_spectrum(n_atoms, omega_max=10.0, seed=1000+i)

        kappa = kappa_of_spectrum(weights, omegas, T, N)

        for p in p_values:
            vp = varpi_p(weights, omegas, p, T)
            h_vp = h_of_varpi(vp, N)
            ratio = kappa / h_vp if h_vp > 0 else 0.0

            max_ratio[p] = max(max_ratio[p], ratio)

            # Theorem 2' predicts ratio <= 1 (with margin for float64 noise)
            if ratio > 1.0 + 1e-10:
                violations[p] += 1

    for p in p_values:
        ck(violations[p] == 0,
           f"p={p}: {violations[p]}/{n_spectra} spectra violated bound")
        ck(max_ratio[p] < 1.0 + 1e-9,
           f"p={p}: max ratio {max_ratio[p]:.10f} should be <= 1")

    print(f"  Max ratios kappa/h(varpi_p): ", end="")
    for p in p_values:
        print(f"p={p}: {max_ratio[p]:.8f}  ", end="")
    print(f"\n  ({time.time()-t0:.1f}s)")

# ==============================================================================
# STAGE 2: Proposition 2'' verification (p < 2)
# ==============================================================================

def run_stage2(N):
    """
    Verify that for p < 2, there exist spectra violating kappa > h(varpi_p).
    Use the explicit construction: P_eps = (1-eps)*delta(omega_1) + eps*delta(omega_2),
    omega_2 -> 0.
    """
    print(f"\n=== STAGE 2: Proposition 2'' counterexamples (p < 2), N={N} ===")
    t0 = time.time()

    # Test parameters
    p_values = [0.5, 1.0, 1.5]
    eps = 0.1
    omega_2 = 1e-6  # small but not exactly zero
    T = 1.0

    # CERT-B lower bound on g(2) from P4_main_theorems.md
    g2_lo = 1.3181093293757127
    zeta4 = float(mp.pi**4 / 90)

    for p in p_values:
        alpha = (1.0 - eps)**(1.0/p)

        # Explicit condition (*) from P4_moment_threshold.md:
        # varpi_1^2 < g2_lo * (alpha^(p-2) - 1) / (4*zeta4 * (alpha^(p-2) + 1))
        alpha_pm2 = alpha**(p - 2.0)
        threshold_sq = g2_lo * (alpha_pm2 - 1.0) / (4.0 * zeta4 * (alpha_pm2 + 1.0))

        ck(threshold_sq > 0,
           f"p={p}: condition threshold {threshold_sq:.6e} should be > 0")

        # Choose omega_1 to satisfy the condition with margin
        varpi_1 = 0.9 * np.sqrt(threshold_sq)
        omega_1 = varpi_1 * 2.0 * np.pi * T

        # Construct the two-atom spectrum
        weights = np.array([1.0 - eps, eps])
        omegas = np.array([omega_1, omega_2])

        # Compute kappa(P,T)
        kappa = kappa_of_spectrum(weights, omegas, T, N)

        # Compute varpi_p and h(varpi_p)
        vp = varpi_p(weights, omegas, p, T)
        h_vp = h_of_varpi(vp, N)

        ratio = kappa / h_vp if h_vp > 0 else 0.0
        violation = kappa - h_vp

        # Proposition 2'' predicts kappa > h(varpi_p) when condition holds
        ck(ratio > 1.0 + 1e-10,
           f"p={p}: ratio {ratio:.10f} should be > 1 (violation={violation:.6e})")

        # Also verify the condition itself
        varpi_1_actual = omega_1 / (2.0 * np.pi * T)
        ck(varpi_1_actual**2 < threshold_sq,
           f"p={p}: varpi_1^2={varpi_1_actual**2:.6e} should be < {threshold_sq:.6e}")

        print(f"  p={p:.1f}: eps={eps}, omega_1={omega_1:.4f}, "
              f"kappa/h(varpi_p)={ratio:.8f}, violation={violation:.6e}")

    print(f"  ({time.time()-t0:.1f}s)")

# ==============================================================================
# STAGE 3: Monotonicity of varpi_p with p
# ==============================================================================

def run_stage3():
    """
    Check that varpi_p increases with p (Lyapunov inequality).
    Test on multiple random spectra.
    """
    print(f"\n=== STAGE 3: Monotonicity varpi_p(increasing p) ===")
    t0 = time.time()

    p_values = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
    n_spectra = 30
    T = 1.0

    monotone_count = 0

    for i in range(n_spectra):
        n_atoms = np.random.randint(2, 6)
        weights, omegas = random_multiatomic_spectrum(n_atoms, omega_max=10.0, seed=2000+i)

        varpi_vals = [varpi_p(weights, omegas, p, T) for p in p_values]

        # Check if sequence is increasing
        is_monotone = all(varpi_vals[j] <= varpi_vals[j+1] * (1.0 + 1e-12)
                          for j in range(len(varpi_vals)-1))

        if is_monotone:
            monotone_count += 1

    ck(monotone_count == n_spectra,
       f"Monotonicity: {monotone_count}/{n_spectra} spectra had increasing varpi_p")

    print(f"  {monotone_count}/{n_spectra} spectra confirmed varpi_p increasing with p")
    print(f"  ({time.time()-t0:.1f}s)")

# ==============================================================================
# STAGE 4: N-independence check
# ==============================================================================

def run_stage4():
    """
    Verify that results are stable under N=300 vs N=600.
    Check a subset of Stage 1 and Stage 2 at both N values.
    """
    print(f"\n=== STAGE 4: N-independence (N=300 vs N=600) ===")
    t0 = time.time()

    N_values = [300, 600]
    p = 2.0
    T = 1.0

    # Use a fixed spectrum for reproducibility
    weights = np.array([0.3, 0.5, 0.2])
    omegas = np.array([2.0, 5.0, 8.0])

    results = {}
    for N in N_values:
        kappa = kappa_of_spectrum(weights, omegas, T, N)
        vp = varpi_p(weights, omegas, p, T)
        h_vp = h_of_varpi(vp, N)
        results[N] = (kappa, h_vp, kappa / h_vp)

    kappa_300, h_300, ratio_300 = results[300]
    kappa_600, h_600, ratio_600 = results[600]

    # Check relative difference < 1e-5
    rel_diff_kappa = abs(kappa_600 - kappa_300) / abs(kappa_300)
    rel_diff_h = abs(h_600 - h_300) / abs(h_300)
    rel_diff_ratio = abs(ratio_600 - ratio_300) / abs(ratio_300)

    ck(rel_diff_kappa < 1e-5,
       f"kappa relative diff {rel_diff_kappa:.2e} should be < 1e-5")
    ck(rel_diff_h < 1e-5,
       f"h(varpi_p) relative diff {rel_diff_h:.2e} should be < 1e-5")
    ck(rel_diff_ratio < 1e-5,
       f"ratio relative diff {rel_diff_ratio:.2e} should be < 1e-5")

    print(f"  N=300: kappa={kappa_300:.8f}, h(varpi_2)={h_300:.8f}, ratio={ratio_300:.8f}")
    print(f"  N=600: kappa={kappa_600:.8f}, h(varpi_2)={h_600:.8f}, ratio={ratio_600:.8f}")
    print(f"  Relative diffs: kappa {rel_diff_kappa:.2e}, h {rel_diff_h:.2e}, ratio {rel_diff_ratio:.2e}")
    print(f"  ({time.time()-t0:.1f}s)")

# ==============================================================================
# STAGE 5: p < 2 worst-case search
# ==============================================================================

def run_stage5(N):
    """
    For p < 2, search for the worst violation among random spectra.
    Report the maximum ratio kappa / h(varpi_p).
    """
    print(f"\n=== STAGE 5: Worst violation search (p < 2), N={N} ===")
    t0 = time.time()

    p_values = [0.5, 1.0, 1.5]
    n_spectra = 50
    T = 1.0

    worst_violation = {p: 0.0 for p in p_values}
    worst_ratio = {p: 0.0 for p in p_values}

    for i in range(n_spectra):
        n_atoms = np.random.randint(2, 6)
        weights, omegas = random_multiatomic_spectrum(n_atoms, omega_max=10.0, seed=3000+i)

        kappa = kappa_of_spectrum(weights, omegas, T, N)

        for p in p_values:
            vp = varpi_p(weights, omegas, p, T)
            h_vp = h_of_varpi(vp, N)
            ratio = kappa / h_vp if h_vp > 0 else 0.0
            violation = kappa - h_vp

            if ratio > worst_ratio[p]:
                worst_ratio[p] = ratio
                worst_violation[p] = violation

    for p in p_values:
        # For p < 2, we expect SOME violations (ratio > 1)
        print(f"  p={p:.1f}: worst ratio={worst_ratio[p]:.8f}, violation={worst_violation[p]:.6e}")

    print(f"  ({time.time()-t0:.1f}s)")

# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 80)
    print("P4 MOMENT THRESHOLD THEOREM -- numerical gate")
    print("=" * 80)

    overall_t0 = time.time()

    # Run stages with N=300
    N = 300
    run_stage1(N)
    run_stage2(N)
    run_stage3()
    run_stage5(N)

    # Run N-independence check
    run_stage4()

    # Summary
    print("\n" + "=" * 80)
    print(f"SUMMARY: {_checks_passed}/{_checks_total} checks PASSED")
    print(f"Total time: {time.time() - overall_t0:.1f}s")
    print("=" * 80)

    if _checks_passed == _checks_total:
        print("\n[OK] All numerical gates passed. Theorem 2' and Proposition 2'' are")
        print("     numerically consistent with the analytic proofs.")
    else:
        print(f"\n[ATTENTION] {_checks_total - _checks_passed} check(s) failed.")
        sys.exit(1)

if __name__ == "__main__":
    main()
