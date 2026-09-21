"""Item 4a: exact rank-one resolvent form of the concavity remainder R.

H_mu''(x) = -2 Psi + R,   Psi = <v, O[psi0] v>,
R = 2 sum_{k != top} <w_k, M'v>^2/(h - lam_k),  M = A - 2 mu~ u_c u_c^T.

CLAIM (analytic, derived this session).  Let h = lam_1(M), v its top vector,
y_perp := (I - v v^T) M'v.  Then, with A's OWN resolvent,

        R = 2 < y_perp , (h - A)^{-1} y_perp >                         (*)

i.e. M need not be diagonalised at all.  Proof: G z (z _|_ v) is the unique
w _|_ v with (h-M)w = z.  Writing (h-M) = (h-A) + 2mu~ u_c u_c^T and using
u_c.(h-A)^{-1}u_c = sum_j c_j^2/(h-mu_j) = -1/(2mu~)  [secular equation],
the Sherman-Morrison correction has coefficient 1/nu + 2mu~ c_v = 0, so
w = (h-A)^{-1}z + beta v for ANY beta; pairing with y_perp _|_ v kills beta.

In A's eigenbasis, with a_j := <psi_j, y_perp> and h in (mu_2, mu_1):

        R = 2 [ sum_{j>=2} a_j^2/(h - mu_j)  -  a_1^2/(mu_1 - h) ]     (**)

-- the j=1 term is NEGATIVE and is a genuine, keepable gain.  Sound bound:

        R <= 2 ( ||y_perp||^2 - a_1^2 ) / (h - mu_2)                   (B0)
        R <= 2 [ sum_{j>=2} a_j^2/(h-mu_j) ] - 2 a_1^2/(mu_1-h)        (B1 = exact)

This script verifies (*) and (**) against a direct eigh of M, and reports how
much the psi_1 term and the exact-vs-worst-gap distribution are worth.
Strictness: identity (*) is analytic; all numbers here are numerical-evidence
(float64, N as stated).
"""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh


def build(N, x):
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    ks = np.arange(2*N+1, dtype=float); s = ks**2
    absd = np.abs(n[:, None]-n[None, :]); sumd = n[:, None]+n[None, :]+1
    uu = u[:, None]*u[None, :]

    def O(f):
        Am = (f[absd]+f[sumd])*uu
        L1 = np.zeros(N); L1[1:] = np.cumsum(2.0*f[1:N])[:N-1]/(2*n[1:]+1)
        np.fill_diagonal(Am, -L1 + f[2*n+1]*u**2)
        return Am
    f = x/(s+x); f[0] = 1.0
    f1 = s/(s+x)**2; f1[0] = 0.0
    psi0 = s/(s+x)**3; psi0[0] = 0.0
    return n, u, O(f), O(f1), O(psi0)


CASES = [(0.09, 6), (0.04, 6), (0.25, 6), (0.09, 20), (0.25, 45), (1.0, 6)]
MUTS = (0.02, 0.05, 0.13, 0.26, 0.5)

if __name__ == "__main__":
    for N in (300, 600):
        print(f"\n===== N = {N} "
              f"(R_dir: eigh(M);  R_id: identity (*);  R_B0: worst-gap bound) =====")
        print(f"{'x':>5} {'Nc':>3} {'mu~':>5} | {'h-mu2':>9} {'mu1-h':>9} | "
              f"{'R_dir':>10} {'R_id':>10} {'relerr':>9} | {'psi1 gain':>9} "
              f"{'R_B0':>10} {'B0/R':>6} | {'Psi':>7} {'H\"':>9}")
        for x, Nc in CASES:
            n, u, A, Mp, Opsi = build(N, x)
            mu, PSI = eigh(A)
            uc = np.where(n < Nc, u, 0.0)
            for mut in MUTS:
                M = A - 2.0*mut*np.outer(uc, uc)
                lam, W = eigh(M)
                h = lam[-1]; v = W[:, -1]
                Psi = float(v @ (Opsi @ v))
                # direct
                c = W[:, :-1].T @ (Mp @ v)
                R_dir = float(np.sum(2.0*c**2/(h - lam[:-1])))
                # identity (*): A's own resolvent, M never diagonalised
                y = Mp @ v; yp = y - float(v @ y)*v
                a = PSI.T @ yp
                R_id = float(np.sum(2.0*a**2/(h - mu)))
                # split
                neg = 2.0*a[-1]**2/(mu[-1] - h)          # the psi_1 term, subtracted
                pos = R_id + neg
                R_B0 = 2.0*(float(yp @ yp) - a[-1]**2)/(h - mu[-2])
                rel = abs(R_id - R_dir)/max(abs(R_dir), 1e-300)
                print(f"{x:>5} {Nc:>3} {mut:>5.2f} | {h-mu[-2]:>9.2e} {mu[-1]-h:>9.2e} | "
                      f"{R_dir:>10.4e} {R_id:>10.4e} {rel:>9.2e} | {neg/pos:>9.4f} "
                      f"{R_B0:>10.4e} {R_B0/R_dir:>6.2f} | {Psi:>7.4f} {-2*Psi+R_dir:>9.3e}")
