"""Question 2: Understand non-monotonic Nc dependence.

Observation: Nc=6 shows worst ratio; Nc=45 recovers. Why?

Hypothesis: Nc controls overlap structure — c1 = <uc, v> and <uc, w2>.
Check if ratio is organized by kappa := 2*mu_t*c1^2/gamma (repo's near-crossing criterion).

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


print("=" * 110)
print("QUESTION 2: Non-monotonic Nc dependence at fixed mu_t=0.13, varying x")
print("=" * 110)
print()

mu_t = 0.13
x_values = [0.01, 0.09, 0.25, 1.0]
Nc_values = [2, 3, 4, 6, 8, 12, 20, 30, 45, 64, 100, 200]

print(f"{'Nc':>4} {'x':>8} {'varpi':>6} {'gap':>8} {'c1':>8} {'c1^2':>8} "
      f"{'kappa':>8} {'<uc,w2>':>9} {'ratio':>7} {'R':>10} {'E':>10} {'sign_changes':>14}")

for x in x_values:
    print(f"\n--- x = {x:.4f}, varpi = {np.sqrt(x):.3f} ---")
    for Nc in Nc_values:
        s = _ks ** 2
        f = f_ein(x)
        A_ = O(f)
        B_ = O(f ** 2)
        Cm, uc = coulomb(mu_t, Nc)
        M = A_ - Cm
        lam, W = eigh(M)
        v = W[:, -1]
        w2 = W[:, -2]
        gam = lam[-1] - lam[-2]

        # Overlaps
        c1 = abs(float(uc @ v))
        uc_w2 = float(uc @ w2)
        kappa = 2.0 * mu_t * c1 ** 2 / gam if gam > 1e-12 else 0.0

        # Ratio computation
        Mp = (A_ - B_) / x
        Psi = float(v @ (O(s / (s + x) ** 3) @ v))
        cj = W[:, :-1].T @ (Mp @ v)
        R = float(2.0 * np.sum(cj ** 2 / (lam[-1] - lam[:-1])))
        E = 2 * x * Psi
        ratio = (x * R) / E if E > 1e-12 else 0.0

        # Sign changes in v (cone exit indicator)
        sign_changes = np.sum(np.diff(np.sign(v)) != 0)

        print(f"{Nc:>4} {x:>8.4f} {np.sqrt(x):>6.3f} {gam:>8.4f} {c1:>8.5f} "
              f"{c1**2:>8.5f} {kappa:>8.3f} {uc_w2:>9.5f} {ratio:>7.3f} "
              f"{R:>10.3e} {E:>10.3e} {sign_changes:>14}")

print("\n" + "=" * 110)
print("Look for:")
print("  - Peak in ratio at intermediate Nc (expected ~6-8)")
print("  - Correlation between kappa and ratio")
print("  - Correlation between sign_changes and ratio")
print("  - How c1^2 and gap evolve with Nc")
