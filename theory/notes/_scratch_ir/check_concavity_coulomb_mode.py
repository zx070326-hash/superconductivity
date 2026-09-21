"""Identify the single mode that carries 94-100% of R (check_concavity_R_Nstability.py):
hypothesis = the COULOMB MODE, i.e. the isolated most-negative eigenvalue that
the rank-one repulsion -2mu~ u_c u_c^T pushes out of the bulk, with eigenvector
w_c ~ u_c/||u_c|| and lam_c ~ -2 mu~ ||u_c||^2.  If so
      R ~ 2 <w_c, M'v>^2 / (h - lam_c) ~ 2 <u_c,M'v>^2 / (||u_c||^2 (h + 2mu~||u_c||^2))
and the concavity ratio R/(2Psi) should scale like 1/mu~ at fixed (x,N_c) once
the mode is isolated.  Strictness: numerical-evidence (N=400).
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

N = 400
n = np.arange(N); u = 1.0 / np.sqrt(2.0 * n + 1.0)
ks = np.arange(2 * N + 1, dtype=float); s = ks ** 2
absd = np.abs(n[:, None] - n[None, :]); sumd = n[:, None] + n[None, :] + 1
uu = u[:, None] * u[None, :]


def O(f):
    A = (f[absd] + f[sumd]) * uu
    L1 = np.zeros(N); L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * n + 1] * u ** 2)
    return A


print(f"{'x':>5} {'Nc':>3} {'mu~':>5} | {'lam_min':>9} {'-2mu~|uc|^2':>11} {'|<w_min,uc^>|':>13} "
      f"{'lam_2nd_min':>11} | {'R':>9} {'R_1mode':>9} {'R_approx':>9} {'Psi':>8} {'ratio':>6} {'ratio*mu~':>9}")
for x, Nc in ((0.09, 6), (0.04, 6), (0.25, 6), (0.09, 20), (0.25, 45), (1.0, 6)):
    f = x / (s + x); f[0] = 1.0
    A = O(f)
    f1 = s / (s + x) ** 2; f1[0] = 0.0
    psi0 = s / (s + x) ** 3; psi0[0] = 0.0
    Mp = O(f1); Opsi = O(psi0)
    uc = np.where(n < Nc, u, 0.0); nuc2 = float(uc @ uc); uch = uc / math.sqrt(nuc2)
    for mut in (0.02, 0.05, 0.13, 0.26, 0.5):
        M = A - 2.0 * mut * np.outer(uc, uc)
        lam, W = eigh(M)
        h = lam[-1]; v = W[:, -1]
        Psi = float(v @ (Opsi @ v))
        c = W[:, :-1].T @ (Mp @ v)
        terms = 2.0 * c ** 2 / (h - lam[:-1]); R = float(terms.sum())
        R1 = float(terms[0])
        ov = abs(float(W[:, 0] @ uch))
        Rapx = 2.0 * float(uc @ (Mp @ v)) ** 2 / (nuc2 * (h + 2 * mut * nuc2))
        print(f"{x:>5} {Nc:>3} {mut:>5.2f} | {lam[0]:>9.4f} {-2*mut*nuc2:>11.4f} {ov:>13.4f} "
              f"{lam[1]:>11.4f} | {R:>9.3e} {R1:>9.3e} {Rapx:>9.3e} {Psi:>8.3f} {R/(2*Psi):>6.3f} {R/(2*Psi)*mut:>9.4f}")
    print()
