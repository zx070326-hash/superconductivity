"""Question 4: Search for ratio < 1 sufficiency criterion.

Try to organize ratio as a function of (mu_t, Nc, x) or kappa.
Look for empirical scaling law.

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
print("QUESTION 4: Search for ratio < 1 criterion")
print("=" * 110)
print()

# Systematic scan
data = []
for mu_t in [0.05, 0.08, 0.10, 0.13, 0.16, 0.20, 0.26]:
    for Nc in [6, 8, 12, 20, 30, 45, 64, 100]:
        for x in [0.01, 0.04, 0.09, 0.16, 0.25, 0.49, 1.0, 4.0]:
            s = _ks ** 2
            f = f_ein(x)
            A_ = O(f)
            B_ = O(f ** 2)
            Cm, uc = coulomb(mu_t, Nc)
            M = A_ - Cm
            lam, W = eigh(M)
            v = W[:, -1]
            gam = lam[-1] - lam[-2]

            c1 = abs(float(uc @ v))
            kappa = 2.0 * mu_t * c1 ** 2 / gam if gam > 1e-12 else 0.0

            Mp = (A_ - B_) / x
            Psi = float(v @ (O(s / (s + x) ** 3) @ v))
            cj = W[:, :-1].T @ (Mp @ v)
            R = float(2.0 * np.sum(cj ** 2 / (lam[-1] - lam[:-1])))
            E = 2 * x * Psi
            ratio = (x * R) / E if E > 1e-12 else 0.0

            varpi = np.sqrt(x)
            data.append((mu_t, Nc, x, varpi, gam, c1**2, kappa, ratio))

# Sort by ratio descending to find worst cases
data.sort(key=lambda row: row[-1], reverse=True)

print("Top 30 worst cases (ratio > 1):")
print(f"{'mu_t':>6} {'Nc':>4} {'x':>8} {'varpi':>6} {'gap':>8} {'c1^2':>8} "
      f"{'kappa':>8} {'ratio':>7}")
for row in data[:30]:
    print(f"{row[0]:>6.2f} {row[1]:>4} {row[2]:>8.4f} {row[3]:>6.3f} "
          f"{row[4]:>8.4f} {row[5]:>8.5f} {row[6]:>8.3f} {row[7]:>7.3f}")

print("\n" + "=" * 110)
print("\nLooking for patterns:")

# Check if varpi threshold works
pass_data = [row for row in data if row[-1] < 1.0]
fail_data = [row for row in data if row[-1] >= 1.0]

if pass_data:
    min_varpi_pass = min(row[3] for row in pass_data)
    print(f"  Minimum varpi among passing cases: {min_varpi_pass:.3f}")
if fail_data:
    max_varpi_fail = max(row[3] for row in fail_data)
    print(f"  Maximum varpi among failing cases: {max_varpi_fail:.3f}")

# Check if Nc threshold works at fixed small varpi
small_varpi_data = [row for row in data if row[3] < 0.35]
if small_varpi_data:
    fail_small = [row for row in small_varpi_data if row[-1] >= 1.0]
    pass_small = [row for row in small_varpi_data if row[-1] < 1.0]
    if fail_small:
        max_Nc_fail = max(row[1] for row in fail_small)
        print(f"  At varpi < 0.35: Maximum Nc among failures: {max_Nc_fail}")
    if pass_small:
        min_Nc_pass = min(row[1] for row in pass_small)
        print(f"  At varpi < 0.35: Minimum Nc among passes: {min_Nc_pass}")

# Try simple criterion: varpi >= 0.5 OR Nc >= 45
simple_rule_pass = sum(1 for row in data if (row[3] >= 0.5 or row[1] >= 45) and row[-1] < 1.0)
simple_rule_fail = sum(1 for row in data if (row[3] >= 0.5 or row[1] >= 45) and row[-1] >= 1.0)
simple_rule_total = simple_rule_pass + simple_rule_fail
print(f"\nSimple rule: varpi >= 0.5 OR Nc >= 45")
print(f"  Passes: {simple_rule_pass}/{simple_rule_total} "
      f"({100*simple_rule_pass/simple_rule_total:.1f}% if total > 0)")

# Try kappa-based criterion
kappa_threshold_candidates = [0.2, 0.3, 0.5, 0.7, 1.0]
print(f"\nKappa-based criteria:")
for kappa_th in kappa_threshold_candidates:
    below_pass = sum(1 for row in data if row[6] < kappa_th and row[-1] < 1.0)
    below_fail = sum(1 for row in data if row[6] < kappa_th and row[-1] >= 1.0)
    below_total = below_pass + below_fail
    if below_total > 0:
        print(f"  kappa < {kappa_th:.1f}: {below_pass}/{below_total} pass "
              f"({100*below_pass/below_total:.1f}%)")
