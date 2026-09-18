"""EXHAUSTIVE test of Einstein extremality at mu*>0.

Reduction that makes this exhaustive rather than a sample:
    M(P) = O[F_P] - 2 mu_t |u_c><u_c|  is AFFINE in P (O is linear in P, the
    Coulomb term is P-independent), so kappa(P) = lam1(M(P)) is a CONVEX
    function of P (sup of affine functionals).  A convex function on the
    convex set  {P >= 0 : mass = 1, mean = xbar}  attains its MAXIMUM at an
    extreme point, and the extreme points of a set cut out by 2 linear
    constraints are measures supported on at most 2 points.
    ==> scanning all 2-atom measures with matched (mass, mean) is exhaustive.

Theorem E at mu*>0 is TRUE at (xbar, mu_t, Nc) iff
    max over 2-atom P of  kappa(P)  =  kappa(delta_xbar) .

A positive entry in the "best excess" column is a COUNTEREXAMPLE.
Strictness: numerical-evidence (float64, finite basis N, finite grid).
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigvalsh

N = 400
VPMAX = 40.0          # basis validity: need N >~ 5*varpi
_n = np.arange(N)
_u = 1.0 / np.sqrt(2.0 * _n + 1.0)
_ks = np.arange(2 * N + 1, dtype=float)
_absd = np.abs(_n[:, None] - _n[None, :])
_sumd = _n[:, None] + _n[None, :] + 1
_uu = _u[:, None] * _u[None, :]


def kappa(atoms, mu_t, Nc):
    """atoms: [(w, x)], x = varpi^2."""
    f = np.zeros(2 * N + 1)
    for w, x in atoms:
        f += w * x / (_ks ** 2 + x)
    f[0] = sum(w for w, _ in atoms)
    A = (f[_absd] + f[_sumd]) * _uu
    L1 = np.zeros(N)
    L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * _n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * _n + 1] * _u ** 2)
    if mu_t:
        uc = np.where(_n < (Nc if Nc else N), _u, 0.0)
        A = A - 2.0 * mu_t * np.outer(uc, uc)
    return float(eigvalsh(A)[-1])


def scan(xbar, mu_t, Nc, nw=25, nx=40):
    kE = kappa([(1.0, xbar)], mu_t, Nc)
    best, arg = -1e9, None
    for w in np.linspace(0.02, 0.98, nw):
        x1max = min(xbar / w, VPMAX ** 2)
        for x1 in np.geomspace(1e-4 * xbar, x1max * 0.999, nx):
            x2 = (xbar - w * x1) / (1.0 - w)
            if x2 <= 0 or x2 > VPMAX ** 2:
                continue
            e = kappa([(w, x1), (1 - w, x2)], mu_t, Nc) - kE
            if e > best:
                best, arg = e, (w, x1, x2)
    return kE, best, arg


if __name__ == "__main__":
    print(f"N={N}  (exhaustive over 2-atom extreme points; see docstring)")
    print(f"{'xbar':>8} {'varpi':>7} {'mu_t':>6} {'Nc':>5} {'kappa_E':>10} "
          f"{'best excess':>13} {'w*':>6} {'x1*':>9} {'x2*':>9}")
    for mu_t, Nc in ((0.0, None), (0.05, 6), (0.13, 6), (0.26, 6),
                     (0.05, 45), (0.13, 45), (0.13, 200)):
        for xbar in (0.01, 0.09, 0.1089, 0.25, 1.0, 9.0, 100.0):
            kE, best, arg = scan(xbar, mu_t, Nc)
            w, x1, x2 = arg if arg else (float('nan'),) * 3
            flag = "   <<< COUNTEREXAMPLE" if best > 1e-12 else ""
            print(f"{xbar:>8.4f} {math.sqrt(xbar):>7.4f} {mu_t:>6.2f} {str(Nc):>5} "
                  f"{kE:>10.6f} {best:>13.3e} {w:>6.3f} {x1:>9.4f} {x2:>9.3f}{flag}")
        print()
