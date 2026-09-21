"""Verify the CONVEXITY BOUND -- a one-line, unconditional, all-mu* upper bound
that bypasses Lemma M entirely.

    M(P) = O[F_P] - 2 mu_t |u_c><u_c|

O[.] is LINEAR in the channel profile F, F_P is LINEAR in P, and the Coulomb
term is P-INDEPENDENT.  So P |-> M(P) is affine, and

    kappa(P) := lam1(M(P)) = sup_{||v||=1} <v, M(P) v>

is a supremum of affine functionals of P, hence CONVEX in P.  Therefore for
P = sum_i w_i delta_{x_i} with sum w_i = 1:

    kappa(P)  <=  sum_i w_i kappa(delta_{x_i})  =  <h_mu>_P            (CVX)

i.e. the eigenvalue of ANY spectrum is at most the P-average of the Einstein
eigenvalues at the SAME mu_t and N_c.  This is analytic and unconditional --
it needs no cone, no Perron structure, no kernel monotonicity, and it holds
for every mu* >= 0 and every varpi.

Relation to Theorem E:  (CVX) then  <h_mu>_P <= h_mu(sqrt(xbar))  [Jensen,
needs only concavity of H_mu(x) = h_mu(sqrt x)]  ==>  Theorem E.
So (CVX) is the mu*-robust half of a proof route that replaces E1+E4.

This script checks (CVX) numerically and measures how much it beats / loses to
the Theorem E bound h_mu(sqrt(xbar)).
Strictness: numerical-evidence (float64, finite basis).
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigvalsh

N = 400
_n = np.arange(N)
_u = 1.0 / np.sqrt(2.0 * _n + 1.0)
_ks = np.arange(2 * N + 1, dtype=float)
_absd = np.abs(_n[:, None] - _n[None, :])
_sumd = _n[:, None] + _n[None, :] + 1
_uu = _u[:, None] * _u[None, :]


def kappa(atoms, mu_t, Nc):
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


if __name__ == "__main__":
    rng = np.random.default_rng(20260916)
    print("(CVX) check: kappa(P) <= <h_mu>_P ?   violation => bound is WRONG")
    print(f"{'mu_t':>6} {'Nc':>5} {'natoms':>7} {'xbar':>9} {'kappa(P)':>10} "
          f"{'<h>_P':>10} {'h(sqrt xbar)':>13} {'CVX slack':>11} {'ThmE slack':>11}")
    worst_cvx = 1e9
    for mu_t, Nc in ((0.0, None), (0.05, 6), (0.13, 6), (0.26, 6),
                     (0.13, 45), (0.13, 200)):
        for trial in range(6):
            na = int(rng.integers(2, 6))
            xs = np.exp(rng.uniform(math.log(1e-3), math.log(400.0), na))
            ws = rng.dirichlet(np.ones(na))
            atoms = list(zip(ws, xs))
            xbar = float(ws @ xs)
            kP = kappa(atoms, mu_t, Nc)
            hbar = float(sum(w * kappa([(1.0, x)], mu_t, Nc) for w, x in atoms))
            hE = kappa([(1.0, xbar)], mu_t, Nc)
            s_cvx, s_thmE = hbar - kP, hE - kP
            worst_cvx = min(worst_cvx, s_cvx)
            bad = "  <<< CVX VIOLATED" if s_cvx < -1e-12 else ""
            bad += "  <<< ThmE VIOLATED" if s_thmE < -1e-12 else ""
            print(f"{mu_t:>6.2f} {str(Nc):>5} {na:>7d} {xbar:>9.3f} {kP:>10.6f} "
                  f"{hbar:>10.6f} {hE:>13.6f} {s_cvx:>11.3e} {s_thmE:>11.3e}{bad}")
    print(f"\nworst (CVX) slack over all trials: {worst_cvx:.3e}  "
          f"(must be >= 0)")
