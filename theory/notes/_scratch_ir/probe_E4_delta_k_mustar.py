"""Decisive test for Theorem E at mu*>0, replacing the H''(x) route.

WHY H''(x) IS THE WRONG DIAGNOSTIC
----------------------------------
Theorem E's chain is (E1) per-channel concavity + (E4) kernel monotonicity.
Per channel, f_x(k) = x/(k^2+x) has d2f/dx2 = -2k^2/(k^2+x)^3 < 0 for k>=1
and f_x(0) = 1 is x-independent.  So for ANY measure P with mean x,
    F_P(k) <= f_x(k)  for every k,   with equality at k=0.          (E1)
This is entirely mu*-independent -- the Coulomb term does not touch f.
Therefore extremality at mu*>0 needs ONLY:

    delta >= 0 with delta(0) = 0   ==>   lam1(M + O[delta]) >= lam1(M)   (E4)

where M = O[f] - 2 mu_t |u_c><u_c|.  By the variational principle it suffices
that, for v = top eigenvector of M,

    <v, O[delta_k] v>  >=  0    for every k >= 1.                    (*)

The repo proves (*) via cone invariance (v >= 0 and nonincreasing), which is
what fails at mu*>0.  But the cone is SUFFICIENT, not necessary.  This script
tests (*) directly on the actual Perron vector of the mu*-perturbed operator.

O[delta_k] entries (delta_k = indicator of channel k, k >= 1):
    off-diag n!=m : (1[|n-m|=k] + 1[n+m+1=k]) u_n u_m
    diag    n=m   : -(2/(2n+1)) 1[k<=n] + 1[2n+1=k] u_n^2
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh


def build_from_f(f, N):
    """Operator O[f] for an arbitrary channel profile f (length >= 2N)."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    A = (f[np.abs(n[:, None] - n[None, :])] + f[n[:, None] + n[None, :] + 1]) \
        * (u[:, None] * u[None, :])
    L1 = np.zeros(N)
    L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * n + 1] * u ** 2)
    return A, u


def f_einstein(x, N):
    ks = np.arange(2 * N + 1, dtype=float)
    f = x / (ks ** 2 + x)
    f[0] = 1.0
    return f


def top_vec(x, mu_t, Nc, N):
    f = f_einstein(x, N)
    A, u = build_from_f(f, N)
    if mu_t:
        uc = np.where(np.arange(N) < (Nc if Nc else N), u, 0.0)
        A = A - 2.0 * mu_t * np.outer(uc, uc)
    w, V = eigh(A)
    v = V[:, -1]
    return w[-1], v * np.sign(v[np.argmax(np.abs(v))])


def O_delta(k, N):
    e = np.zeros(2 * N + 1)
    e[k] = 1.0
    return build_from_f(e, N)[0]


if __name__ == "__main__":
    N = 200
    print("(*) min over k>=1 of <v, O[delta_k] v> for the ACTUAL Perron vector")
    print("    negative => E4 fails at that point => Theorem E route broken there\n")
    print(f"{'x':>8} {'varpi':>7} {'mu_t':>6} {'Nc':>5} {'lam1':>10} "
          f"{'min_k <v,Od_k,v>':>18} {'argmin k':>9} {'sign chg':>8} {'v>=0?':>6}")
    Ods = {k: O_delta(k, N) for k in range(1, 61)}
    for mu_t, Nc in ((0.0, None), (0.05, 6), (0.05, 45), (0.13, 6),
                     (0.13, 45), (0.26, 6), (0.05, 200), (0.13, 200)):
        for x in (0.01, 0.09, 0.1089, 0.25, 1.0, 9.0, 100.0):
            lam, v = top_vec(x, mu_t, Nc, N)
            vals = {k: float(v @ (Ods[k] @ v)) for k in Ods}
            kmin = min(vals, key=vals.get)
            s = np.sign(v[np.abs(v) > 1e-14])
            nchg = int(np.sum(s[1:] != s[:-1]))
            allpos = bool(np.all(v >= -1e-14))
            flag = "  <<< NEGATIVE" if vals[kmin] < 0 else ""
            print(f"{x:>8.4f} {math.sqrt(x):>7.4f} {mu_t:>6.2f} {str(Nc):>5} "
                  f"{lam:>10.6f} {vals[kmin]:>18.6e} {kmin:>9d} {nchg:>8d} "
                  f"{str(allpos):>6}{flag}")
        print()
