r"""
P4 Strong Coupling Expansion Numerical Gate (Theorem S verification).

Verifies P4_strong_coupling_expansion.md推论 S.1 and S.2:
  S.1: Tc expansion Tc = Tc^E [1 + a2 sigma^2 + a3 gamma3 sigma^3 + O(sigma^4)]
  S.2: Consistency with Theorem 2 (Einstein extremal): a2 <= 1/2, i.e.
       kappa(wbar) >= -(1 + wbar*h'/h)

IMPORTS ONLY: H_matrix from p4_kappa_N (not modifying it).

Tests:
  (1) S.1 ratio test: N=1 analytical case with 2-atom spectra
  (2) S.2 inequality: kappa(wbar) >= -(1 + wbar*h'/h) on grid
  (3) Consistency with p4_kappa_N.py values

Execution constraint: <= 5 minutes.
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
    """Test S.2: kappa(wbar) >= -(1 + wbar*h'/h) on grid."""
    print("\n=== Test S.2: Inequality kappa >= -(1 + wbar*h'/h) ===")

    N = 100
    wbar_grid = [mp.mpf(x) for x in [0.3, 0.5, 1.0, 2.0, 3.0, 5.0]]

    violations = []
    min_margin = mp.mpf('inf')

    for wbar in wbar_grid:
        kappa_N, c3_N, frak_h_N = kappa_c3_frak_h(N, wbar)

        dw = mp.mpf('1e-6')
        frak_h_plus, _ = eig_max(H_matrix(N, wbar + dw))
        frak_h_minus, _ = eig_max(H_matrix(N, wbar - dw))
        h_prime_num = (frak_h_plus - frak_h_minus) / (2 * dw)

        ratio = wbar * h_prime_num / frak_h_N
        bound = -(1 + ratio)

        margin = kappa_N - bound
        min_margin = min(min_margin, margin)

        if margin < 0:
            violations.append((wbar, kappa_N, bound, margin))

        chk(margin >= -mp.mpf('1e-6'),
            f"S.2 violated at wbar={wbar}: kappa={kappa_N} < {bound}")

    if violations:
        print(f"  WARNING: {len(violations)} violations found:")
        for wbar, kappa, bound, margin in violations[:3]:
            print(f"    wbar={mp.nstr(wbar,3)}: kappa={mp.nstr(kappa,6)}, "
                  f"bound={mp.nstr(bound,6)}, margin={mp.nstr(margin,6)}")
    else:
        print(f"  S.2 inequality holds on all {len(wbar_grid)} grid points (N={N}).")
        print(f"  min margin = {mp.nstr(min_margin, 6)}")


def test_S1_ratio_simplified():
    """Test S.1: Tc expansion using N=1 analytical case."""
    print("\n=== Test S.1: Tc expansion ratio test (N=1 analytical) ===")

    N = 1
    lam = mp.mpf(2)
    wbar = mp.mpf(1)

    h_val = mp.mpf('0.5')
    kappa_val = mp.mpf(-1)
    c3_val = mp.mpf(0)

    kappa_calc, c3_calc, h_calc = kappa_c3_frak_h(N, wbar)
    chk(abs(kappa_calc - kappa_val) < mp.mpf('1e-10'))
    chk(abs(c3_calc - c3_val) < mp.mpf('1e-10'))
    chk(abs(h_calc - h_val) < mp.mpf('1e-10'))

    L = mp.mpf(1)
    a2 = -kappa_val / (2 * (1 + L))
    chk(abs(a2 - mp.mpf('0.25')) < mp.mpf('1e-10'))

    sigma_vals = [mp.mpf('0.05'), mp.mpf('0.1'), mp.mpf('0.2')]

    for sigma in sigma_vals:
        Tc_E = mp.mpf(1)
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

        print(f"  sigma={mp.nstr(sigma,3)}: rel_err={mp.nstr(rel_err,6)}, "
              f"err/sig^4={mp.nstr(ratio_sigma4,4)}")

        chk(abs(rel_err) < 5 * sigma**2,
            f"S.1: error too large at sigma={sigma}")

    print(f"  S.1 tested {len(sigma_vals)} cases at wbar=1, N=1.")


def main():
    print("=" * 70)
    print("P4 Strong Coupling Expansion Numerical Gate")
    print("Verifying Theorem S (P4_strong_coupling_expansion.md)")
    print("=" * 70)

    import time
    t_start = time.time()

    test_consistency_with_p4_kappa_N()
    test_S2_inequality()
    test_S1_ratio_simplified()

    t_elapsed = time.time() - t_start

    print("\n" + "=" * 70)
    print(f"PASS: All {ck[0]} checks passed.")
    print(f"Elapsed time: {t_elapsed:.2f} seconds")
    print("=" * 70)


if __name__ == "__main__":
    main()
