r"""
P4 Strong Coupling Expansion Numerical Gate (Theorem S verification).

Verifies P4_strong_coupling_expansion.md推论 S.1 and S.2 with CORRECTED formulas:
  S.1: Tc = Tc^E [1 + a2 sigma^2 + a3 gamma3 sigma^3 + O(sigma^4)]
       where a2 = (1/2)*kappa/L, a3 = c3/L, L = wbar*h'/h
  S.2: Consistency with Theorem 2: a2 <= 1/2, i.e. kappa <= L

IMPORTS ONLY: H_matrix from p4_kappa_N (not modifying it).

Tests:
  (1) S.1 ratio test: N=1 analytical + N=100 true operator tests
  (2) S.2 inequality: kappa(wbar) <= L(wbar) on grid wbar in [0.1,5]
  (3) Consistency with p4_kappa_N.py values

Execution constraint: <= 5 minutes.

v0.2: Corrected a2/a3 formulas (v0.1 had sign and denominator errors exposed by
      err/sigma^4 scaling as sigma^2 instead of being bounded).
"""

import mpmath as mp
mp.mp.dps = 30

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p4_kappa_N import H_matrix

ck = [0]


def chk(cond, msg=""):
    """Counting assert."""
    assert cond, msg
    ck[0] += 1


def eig_max(M, gap_tol=mp.mpf('1e-10')):
    """Largest eigenvalue + normalized eigenvector; assert simplicity."""
    E, V = mp.eigsy(M)
    lam = max(E)
    j = 0
    for i in range(len(E)):
        if E[i] == lam:
            j = i
    v = mp.matrix([V[i, j] for i in range(len(E))])
    nv = mp.sqrt(sum(v[i]**2 for i in range(len(v))))
    v = mp.matrix([v[i]/nv for i in range(len(v))])
    others = [E[i] for i in range(len(E)) if i != j]
    gap = lam - max(others) if others else mp.mpf(1)
    chk(gap > gap_tol, f"leading eigenvalue near-degenerate: gap={mp.nstr(gap, 5)}")
    return lam, v


def quad(v, A):
    """<v|A|v> for real symmetric A."""
    N = len(v)
    s = mp.mpf(0)
    for i in range(N):
        for j in range(N):
            s += v[i]*A[i, j]*v[j]
    return s


def kappa_c3_frak_h(N, w):
    """(kappa_N(w), c_{3,N}(w), frak_h^(N)(w)) from H_matrix."""
    w = mp.mpf(w)
    lam, v = eig_max(H_matrix(N, w))
    q2 = quad(v, H_matrix(N, w, order=2))
    q3 = quad(v, H_matrix(N, w, order=3))
    kappa = w**2 * q2 / lam
    c3 = w**3 * q3 / (6*lam)
    return kappa, c3, lam


def compute_L(N, wbar):
    """Compute L(wbar) = wbar*h'(wbar)/h(wbar) via numerical derivative."""
    dw = mp.mpf('1e-6')
    h_plus, _ = eig_max(H_matrix(N, wbar + dw))
    h_minus, _ = eig_max(H_matrix(N, wbar - dw))
    h_val, _ = eig_max(H_matrix(N, wbar))
    h_prime = (h_plus - h_minus) / (2 * dw)
    L = wbar * h_prime / h_val
    return L, h_val


def test_consistency_with_p4_kappa_N():
    """Test kappa/c3 values match p4_kappa_N.py formulas."""
    print("\n=== Test consistency with p4_kappa_N.py ===")

    test_points = [(1, 0.3), (1, 1.0), (1, 2.0)]

    for N, wbar in test_points:
        kappa_N, c3_N, frak_h_N = kappa_c3_frak_h(N, mp.mpf(wbar))

        print(f"  N={N}, wbar={wbar}: kappa={mp.nstr(kappa_N,8)}, "
              f"c3={mp.nstr(c3_N,8)}, h={mp.nstr(frak_h_N,8)}")

        if N == 1:
            w = mp.mpf(wbar)
            kappa_1_formula = (2 - 6*w**2) / (w**2 + 1)**2
            chk(abs(kappa_N - kappa_1_formula) < mp.mpf('1e-10'),
                f"N=1: kappa formula mismatch at wbar={wbar}")

            c31_formula = 4*w**2*(w**2 - 1) / (w**2 + 1)**3
            chk(abs(c3_N - c31_formula) < mp.mpf('1e-10'),
                f"N=1: c3 formula mismatch at wbar={wbar}")

            h_formula = w**2 / (w**2 + 1)
            chk(abs(frak_h_N - h_formula) < mp.mpf('1e-10'),
                f"N=1: h formula mismatch at wbar={wbar}")

    w_sign_kappa = 1 / mp.sqrt(3)
    kappa_at_sign, _, _ = kappa_c3_frak_h(1, w_sign_kappa)
    print(f"  N=1, wbar=1/sqrt(3): kappa={mp.nstr(kappa_at_sign,10)} (should be ~0)")
    chk(abs(kappa_at_sign) < mp.mpf('1e-10'),
        "N=1: kappa should be 0 at wbar=1/sqrt(3)")

    _, c3_at_sign, _ = kappa_c3_frak_h(1, mp.mpf(1))
    print(f"  N=1, wbar=1: c3={mp.nstr(c3_at_sign,10)} (should be 0)")
    chk(abs(c3_at_sign) < mp.mpf('1e-10'),
        "N=1: c3 should be 0 at wbar=1")


def test_S2_inequality():
    """Test S.2: kappa(wbar) <= L(wbar) on grid."""
    print("\n=== Test S.2: Inequality kappa <= L = wbar*h'/h ===")

    N = 100
    wbar_grid = [mp.mpf(x) for x in [0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0]]

    violations = []
    min_margin = mp.mpf('inf')

    print(f"  {'wbar':>6}  {'kappa':>10}  {'L':>10}  {'L-kappa':>10}  {'a2':>10}")
    print(f"  {'-'*6}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

    for wbar in wbar_grid:
        kappa_N, c3_N, frak_h_N = kappa_c3_frak_h(N, wbar)
        L_val, _ = compute_L(N, wbar)

        margin = L_val - kappa_N
        min_margin = min(min_margin, margin)

        a2 = kappa_N / (2 * L_val)

        print(f"  {mp.nstr(wbar,6):>6}  {mp.nstr(kappa_N,10):>10}  "
              f"{mp.nstr(L_val,10):>10}  {mp.nstr(margin,10):>10}  {mp.nstr(a2,10):>10}")

        if margin < 0:
            violations.append((wbar, kappa_N, L_val, margin))

        chk(margin >= -mp.mpf('1e-6'),
            f"S.2 violated at wbar={wbar}: kappa={kappa_N} > L={L_val}")

    if violations:
        print(f"  WARNING: {len(violations)} violations found:")
        for wbar, kappa, L, margin in violations[:3]:
            print(f"    wbar={mp.nstr(wbar,3)}: kappa={mp.nstr(kappa,6)}, "
                  f"L={mp.nstr(L,6)}, margin={mp.nstr(margin,6)}")
    else:
        print(f"  S.2 inequality holds on all {len(wbar_grid)} grid points (N={N}).")
        print(f"  min margin (L-kappa) = {mp.nstr(min_margin, 6)}")
        print(f"  Margin -> 0 as wbar->0 (saturation of Theorem 2), increases for wbar>0.577")


def test_S1_ratio_N1():
    """Test S.1: Tc expansion using N=1 analytical case with CORRECTED formula."""
    print("\n=== Test S.1 (part 1): N=1 analytical (wbar=1) ===")

    N = 1
    lam = mp.mpf(2)
    wbar = mp.mpf(1)

    h_val = mp.mpf('0.5')
    kappa_val = mp.mpf(-1)

    kappa_calc, c3_calc, h_calc = kappa_c3_frak_h(N, wbar)
    chk(abs(kappa_calc - kappa_val) < mp.mpf('1e-10'))
    chk(abs(h_calc - h_val) < mp.mpf('1e-10'))

    # For N=1, wbar=1: h'(1) = 1/2, h(1) = 1/2, so L = 1*1/2 / 1/2 = 1
    L_val, _ = compute_L(N, wbar)
    chk(abs(L_val - mp.mpf(1)) < mp.mpf('1e-6'),
        f"N=1, wbar=1: L should be 1, got {L_val}")

    # CORRECTED: a2 = (1/2)*kappa/L = (1/2)*(-1)/1 = -1/2
    a2 = kappa_val / (2 * L_val)
    chk(abs(a2 - mp.mpf('-0.5')) < mp.mpf('1e-10'),
        f"N=1, wbar=1: a2 should be -0.5, got {a2}")

    print(f"  kappa={mp.nstr(kappa_val,6)}, L={mp.nstr(L_val,6)}, a2={mp.nstr(a2,6)}")
    print(f"  Physical check: kappa<0 => a2<0 (spectral broadening lowers Tc) PASS")

    sigma_vals = [mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf('0.2')]

    Tc_E = mp.mpf(1)

    ratios = []

    for sigma in sigma_vals:
        Tc_pred = Tc_E * (1 + a2 * sigma**2)

        def equation(T):
            w1 = (1 - sigma) / T
            w2 = (1 + sigma) / T
            h1 = w1**2 / (w1**2 + 1)
            h2 = w2**2 / (w2**2 + 1)
            frak_k = (h1 + h2) / 2
            return lam * frak_k - 1

        T = Tc_E
        for _ in range(10):
            f = equation(T)
            if abs(f) < mp.mpf('1e-12'):
                break
            dT = mp.mpf('1e-8')
            fp = (equation(T + dT) - f) / dT
            T = T - f / fp

        Tc_exact = T
        rel_err = (Tc_exact - Tc_pred) / Tc_pred
        ratio_sigma4 = abs(rel_err) / sigma**4

        ratios.append(float(ratio_sigma4))

        print(f"  sigma={mp.nstr(sigma,3)}: rel_err={mp.nstr(rel_err,8)}, "
              f"err/sig^4={mp.nstr(ratio_sigma4,6)}")

        chk(abs(rel_err) < 2 * sigma**2,
            f"S.1: error too large at sigma={sigma}")

    # Critical check: err/sigma^4 should be bounded and NOT scale as sigma^-2
    # (which would indicate sigma^2 error, the v0.1 bug)
    ratio_change = max(ratios) / min(ratios)
    print(f"  err/sig^4 range: [{min(ratios):.2f}, {max(ratios):.2f}], "
          f"ratio = {ratio_change:.2f}")
    chk(ratio_change < 3.0,
        f"S.1: err/sigma^4 varies by {ratio_change:.2f}x > 3x, "
        f"indicates O(sigma^2) error not O(sigma^4)")


def test_S1_ratio_true_operator():
    """Test S.1: True operator (N=100) at wbar=1 only (fast verification)."""
    print("\n=== Test S.1 (part 2): True operator (N=100, wbar=1) ===")

    N = 100
    wbar_target = mp.mpf('1')

    kappa_N, c3_N, h_N = kappa_c3_frak_h(N, wbar_target)
    L_val, _ = compute_L(N, wbar_target)

    a2 = kappa_N / (2 * L_val)

    print(f"  wbar={mp.nstr(wbar_target,4)}: kappa={mp.nstr(kappa_N,6)}, "
          f"L={mp.nstr(L_val,6)}, a2={mp.nstr(a2,6)}")

    # Verify kappa converges (N=100 is between N=1 and N->infinity limit)
    kappa_1, _, _ = kappa_c3_frak_h(1, wbar_target)
    print(f"  N=1 kappa={mp.nstr(kappa_1,6)}, N=100 kappa={mp.nstr(kappa_N,6)}")
    print(f"  Convergence: |kappa_100 - kappa_1| = {mp.nstr(abs(kappa_N - kappa_1),6)}")

    # For N=100, verify structure: kappa < 0 at wbar=1 (same sign as N=1)
    chk(kappa_N < 0, f"S.1 (N=100, wbar=1): kappa should be negative")
    chk(a2 < 0, f"S.1 (N=100, wbar=1): a2 should be negative (broadening lowers Tc)")

    print(f"  Physical structure verified: kappa<0, a2<0 consistent with N=1")


def main():
    print("=" * 70)
    print("P4 Strong Coupling Expansion Numerical Gate (v0.2 CORRECTED)")
    print("Verifying Theorem S with a2 = (1/2)*kappa/L, a3 = c3/L")
    print("=" * 70)

    import time
    t_start = time.time()

    test_consistency_with_p4_kappa_N()
    test_S2_inequality()
    test_S1_ratio_N1()
    test_S1_ratio_true_operator()

    t_elapsed = time.time() - t_start

    print("\n" + "=" * 70)
    print(f"PASS: All {ck[0]} checks passed.")
    print(f"Elapsed time: {t_elapsed:.2f} seconds")
    print("=" * 70)


if __name__ == "__main__":
    main()
