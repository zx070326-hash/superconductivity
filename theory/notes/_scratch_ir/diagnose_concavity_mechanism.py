"""Diagnose mechanism of concavity failure at small varpi + small Nc.

Task: decompose R's growth, understand non-monotonic Nc/mu_t dependencies,
      search for ratio < 1 criteria, physical interpretation.

WHY R GROWS
-----------
R = 2 sum_{j>=2} c_j^2 / (h - lam_j),  c_j = <w_j, M' v>.
Three potential sources:
  (i)   denominator collapse:  gap gamma = lam_1 - lam_2
  (ii)  numerator amplification: c_j^2 increases
  (iii) weight redistribution: sum dominated by j=2 term

KEY QUESTION: Is R ~= 2 c_2^2 / gamma  (j=2 dominance)?
If yes, mechanism collapses to two scalars => analytically tractable.

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
    return 2.0 * mu_t * np.outer(uc, uc)


print("=" * 100)
print("QUESTION 1: Decompose R's growth — j=2 dominance vs full sum")
print("=" * 100)
print()

# Pick the worst case from the table: mu_t=0.05, Nc=6, x=0.01 (ratio=5.815)
# and a baseline for comparison: mu_t=0, x=0.01 (ratio=0.000)
print(f"{'mu_t':>5} {'Nc':>4} {'x':>8} {'varpi':>6} {'gap':>8} {'R_total':>10} "
      f"{'c2^2':>10} {'2c2^2/g':>10} {'j=2 %':>7} {'sum_j c_j^2':>12} {'top5 j':>20}")

for mu_t, Nc in [(0.0, None), (0.05, 6), (0.13, 6), (0.13, 45), (0.26, 6)]:
    Cm = coulomb(mu_t, Nc) if mu_t else np.zeros((N, N))
    for x in [0.01, 0.09, 0.1089, 0.16, 0.25, 1.0]:
        s = _ks ** 2
        f = f_ein(x)
        A_ = O(f)
        B_ = O(f ** 2)
        M = A_ - Cm
        lam, W = eigh(M)
        v = W[:, -1]
        w2 = W[:, -2]
        gam = lam[-1] - lam[-2]
        Mp = (A_ - B_) / x
        cj = W[:, :-1].T @ (Mp @ v)
        c2 = cj[-1]  # j=2 is the second-to-last column
        denom_j = lam[-1] - lam[:-1]
        terms = 2.0 * cj ** 2 / denom_j
        R_total = float(np.sum(terms))
        R_j2 = 2.0 * c2 ** 2 / gam
        j2_frac = R_j2 / R_total if R_total > 1e-12 else 0.0
        sum_cj2 = float(np.sum(cj ** 2))

        # Top 5 contributing j indices (counting from j=2 as index 0)
        top5_idx = np.argsort(terms)[-5:][::-1]
        top5_str = ",".join([f"{N-1-i}" for i in top5_idx])

        print(f"{mu_t:>5.2f} {str(Nc):>4} {x:>8.4f} {np.sqrt(x):>6.3f} "
              f"{gam:>8.4f} {R_total:>10.3e} {c2**2:>10.3e} {R_j2:>10.3e} "
              f"{j2_frac:>7.3f} {sum_cj2:>12.3e} {top5_str:>20}")
    print()

print()
print("INTERPRETATION:")
print("  - If 'j=2 %' ~= 1.0 in failure region => R is j=2-dominated => R ~= 2c2^2/gamma")
print("  - If 'j=2 %' << 1.0 => multi-level conspiracy, harder to isolate")
print()
