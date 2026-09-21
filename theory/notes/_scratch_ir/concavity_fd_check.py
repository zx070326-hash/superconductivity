"""Independent verification of H'' = -2 Psi + R, from zero:
H(x) := lam_1( O[f_x] - 2 mu~ u_c u_c^T ) at FIXED N_c; H'' by 5-point central FD
on H itself (no perturbation theory at all).  This is the check that catches the
known trap (using A's eigenvectors instead of M's understates |C| by ~6500x).
numerical-evidence (float64, N=400; FD step s and its halving shown)."""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concavity_exact_resolvent import build
from concavity_region_scan import hpp

N = 400


def H_of_x(x, Nc, mut):
    n, u, A, Mp, Opsi = build(N, x)
    uc = np.where(n < Nc, u, 0.0)
    return float(eigh(A - 2.0 * mut * np.outer(uc, uc), eigvals_only=True)[-1])


def fd2(x, Nc, mut, s):
    f = [H_of_x(x + k * s, Nc, mut) for k in (-2, -1, 0, 1, 2)]
    return (-f[0] + 16 * f[1] - 30 * f[2] + 16 * f[3] - f[4]) / (12 * s * s)


if __name__ == "__main__":
    print(f"{'x':>6} {'Nc':>4} {'mu~':>6} {'H\"(theory)':>14} {'H\"(FD,s)':>14} "
          f"{'H\"(FD,s/2)':>14} {'rel err':>10}")
    worst = 0.0
    for x, Nc in ((0.09, 6), (0.04, 6), (0.25, 6), (0.25, 45), (1.0, 6),
                  (1.0, 10), (4.0, 20), (0.49, 7)):
        for mut in (0.0, 0.05, 0.13, 0.5):
            th = hpp(N, x, Nc, mut)[0]
            s = 1e-3 * max(x, 1e-2)
            a, b = fd2(x, Nc, mut, s), fd2(x, Nc, mut, s / 2)
            rel = abs(b - th) / max(abs(th), 1e-300)
            worst = max(worst, rel)
            print(f"{x:>6} {Nc:>4} {mut:>6.2f} {th:>14.6e} {a:>14.6e} {b:>14.6e} "
                  f"{rel:>10.2e}" + ("   **" if rel > 1e-4 else ""))
    print(f"\nworst rel |H''(theory) - H''(FD,s/2)| = {worst:.2e}   "
          f"(32 points; FD is fully independent of the resolvent identity)")
