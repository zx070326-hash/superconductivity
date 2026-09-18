"""psi_0 = k^2/(k^2+x)^3 = f^2 (1-f) / x^2, hence

    Psi := <v, O[psi_0] v> = ( <v,O[f^2]v> - <v,O[f^3]v> ) / x^2 .

Derivation (exact, one line each):
    f = x/(k^2+x)                  =>  df/dx  = k^2/(k^2+x)^2 = f(1-f)/x
    d2f/dx2 = d/dx [f(1-f)/x] = -2 f^2 (1-f) / x^2 = -2 k^2/(k^2+x)^3
So psi_0 = -(1/2) d2f/dx2 = f^2(1-f)/x^2.  Note (f^2-f^3)[0] = 1-1 = 0,
matching psi_0(0) = 0, so no k=0 special-casing is needed.

This identifies the extremality functional Psi with the repo's condition (A)
<v,O[f^3]v> < <v,O[f^2]v> recorded in the 2026-09-14 concavity note.
Strictness: numerical-evidence (float64, finite basis).
"""
import os, math
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


print(f"{'x':>9} {'mu_t':>6} {'Nc':>5} {'Psi direct':>14} "
      f"{'(f2-f3)/x^2':>14} {'rel err':>11}")
worst = 0.0
for mu_t, Nc in ((0.0, None), (0.05, 6), (0.13, 45), (0.26, 6)):
    C = coulomb(mu_t, Nc) if mu_t else 0.0
    for x in (0.01, 0.1089, 1.0, 9.0, 100.0):
        s = _ks ** 2
        f = f_ein(x)
        v = eigh(O(f) - C)[1][:, -1]
        Psi = float(v @ (O(s / (s + x) ** 3) @ v))
        alt = float(v @ ((O(f ** 2) - O(f ** 3)) @ v)) / x ** 2
        rel = abs(Psi - alt) / abs(Psi)
        worst = max(worst, rel)
        print(f"{x:>9.4f} {mu_t:>6.2f} {str(Nc):>5} {Psi:>14.6e} "
              f"{alt:>14.6e} {rel:>11.2e}")
    print()
print(f"worst relative deviation: {worst:.2e}   "
      f"({'IDENTITY HOLDS' if worst < 1e-10 else 'IDENTITY FAILS'})")
