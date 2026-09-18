"""Question 3: Understand non-monotonic mu_t dependence.

Observation: At Nc=6, mu_t=0.05 gives ratio=5.8, mu_t=0.13 gives 2.2, mu_t=0.26 gives 1.09.
Decreasing with mu_t looks opposite to another source claiming mu_t>=0.10 all pass at Nc=200.

Reconcile these trends.

Strictness: numerical-evidence (float64, finite N).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

N = 300
_n = np.arange(N)
_u = 1.0 / np.sqrt(2.0 * _n + 1.0)
_ks = np.arange(2 * N + 1, dtype=float)
_absd = np.abs(_n[:, None] - _n[None, :])
_sumd = _n[:, None] + _n[None, :] + 1
_uu = _u[:, None] * _u[None, :]


def O(f):
    A = (f[_absd] + f[_sumd]) * _uu
    L1 = np.zeros(N)
    L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * _n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * _n + 1] * _u ** 2)
    return A


def f_ein(x):
    f = x / (_ks ** 2 + x)
    f[0] = 1.0
    return f


def coulomb(mu_t, Nc):
    uc = np.where(_n < (Nc if Nc else N), _u, 0.0)
    return 2.0 * mu_t * np.outer(uc, uc), uc


print("=" * 100)
print("QUESTION 3: Non-monotonic mu_t dependence")
print("=" * 100)
print()

# Fixed small x (worst case region) and varying Nc
x_test = 0.09
Nc_values = [6, 45, 200]
mu_t_values = [0.0, 0.02, 0.05, 0.08, 0.10, 0.13, 0.16, 0.20, 0.26, 0.30]

print(f"{'Nc':>4} {'mu_t':>6} {'x':>8} {'gap':>8} {'c1^2':>8} {'kappa':>8} "
      f"{'ratio':>7} {'R':>10} {'E':>10} {'sign_chg':>9}")

for Nc in Nc_values:
    print(f"\n--- Nc = {Nc} ---")
    for mu_t in mu_t_values:
        s = _ks ** 2
        f = f_ein(x_test)
        A_ = O(f)
        B_ = O(f ** 2)
        Cm, uc = coulomb(mu_t, Nc)
        M = A_ - Cm
        lam, W = eigh(M)
        v = W[:, -1]
        gam = lam[-1] - lam[-2]

        c1 = abs(float(uc @ v))
        kappa = 2.0 * mu_t * c1 ** 2 / gam if gam > 1e-12 else 0.0

        Mp = (A_ - B_) / x_test
        Psi = float(v @ (O(s / (s + x_test) ** 3) @ v))
        cj = W[:, :-1].T @ (Mp @ v)
        R = float(2.0 * np.sum(cj ** 2 / (lam[-1] - lam[:-1])))
        E = 2 * x_test * Psi
        ratio = (x_test * R) / E if E > 1e-12 else 0.0

        sign_changes = np.sum(np.diff(np.sign(v)) != 0)

        print(f"{Nc:>4} {mu_t:>6.2f} {x_test:>8.4f} {gam:>8.4f} {c1**2:>8.5f} "
              f"{kappa:>8.3f} {ratio:>7.3f} {R:>10.3e} {E:>10.3e} {sign_changes:>9}")

print("\n" + "=" * 100)
print("INTERPRETATION:")
print("  - At small Nc (6): ratio decreases with mu_t (system regularizing)")
print("  - At large Nc (200): ratio may increase with mu_t (different regime)")
print("  - The two trends are NOT contradictory if they operate in different Nc regimes")
