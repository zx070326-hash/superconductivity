"""[numerical-evidence] Explicit end-to-end counterexample: under (Z-dom) the ceiling
CANNOT use the pairing moment M2^D = M2^ph - M2^sf; it must use a Z-channel moment.

Construction (B=1, so no off-diagonal blocks and team-lead's M2' = M2^D exactly):
  a2F^ph = A_ph delta(w - w_ph)  with w_ph SMALL   -> large lam^ph = 2A/w, small M2 = 2Aw
  a2F^sf = A_sf delta(w - w_sf)  with w_sf LARGE   -> small lam^sf,        large M2
Choose A_sf w_sf = A_ph w_ph  ==>  M2^D = 0 exactly, while the pairing kernel
Lam^D(k) = Lam^ph(k) - Lam^sf(k) is still strongly attractive at small k ==> lam1(OO) > 0.
(Z-dom) holds automatically: |Lam^ph - Lam^sf| <= Lam^ph + Lam^sf = Lam^Z.

So  "Tc <= C_inf sqrt(rho(M2^D + 2 offdiag(M2^D)))"  gives a ZERO ceiling at a positive lam1
==> violated by an unbounded factor.  My M2G = (M2^ph + 3 M2^sf) + 2 offdiag(M2^ph + M2^sf)
stays finite and is respected.

DIRECTION: lam1(P_N OO P_N) <= lam1(OO); the ratio is a LOWER estimate, so ratio > 1 is a
sound violation.  kappa = (2 pi T)^2 = 1 throughout (the scaled variable of the repo).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import _chk_ml_build as B
import _chk_ml_absM2 as A

G2 = 1.318140212132


def einstein(lamm, varpi, N, Bmat):
    return B.lam_factorized(Bmat * lamm, varpi, N)


def offd(M):
    return M - np.diag(np.diag(M))


print(f"{'N':>5} {'lam1(OO)':>11} {'M2^D':>9} {'M2G':>9} {'ratio_theirs':>13} {'ratio_mine':>11}")
w_ph, w_sf = 0.10, 5.0
lam_ph = 10.0                      # dimensionless coupling of the phonon piece
M2_ph = lam_ph * w_ph ** 2         # = 0.1
lam_sf = M2_ph / w_sf ** 2         # so that M2_sf = M2_ph  ==>  M2^D = 0
M2_sf = lam_sf * w_sf ** 2
one = np.ones((1, 1))
for N in (100, 200, 400, 800, 1600):
    lph = einstein(lam_ph, w_ph, N, one)
    lsf = einstein(lam_sf, w_sf, N, one)
    OO = A.multiband_O2(lph + lsf, lph - lsf, N, np.array([1.0]))
    l1 = B.lam1(OO)
    M2D = np.array([[M2_ph - M2_sf]])
    M2G = np.array([[M2_ph + 3 * M2_sf]])
    rt = np.inf if abs(M2D[0, 0]) < 1e-14 else l1 / (G2 * abs(M2D[0, 0]))
    print(f"{N:5d} {l1:11.6f} {M2D[0,0]:9.2e} {M2G[0,0]:9.4f} {rt:13.4g} "
          f"{l1/(G2*M2G[0,0]):11.6f}")

# a less extreme, fully finite version: M2^D small but nonzero
print("\nfinite-ratio version (M2^D small but nonzero):")
print(f"{'w_sf':>6} {'lam1':>10} {'M2^D':>9} {'ratio_theirs':>13} {'ratio_mine':>11}")
for w_sf in (1.0, 2.0, 3.0, 4.0, 6.0, 9.0):
    lam_sf = 0.9 * M2_ph / w_sf ** 2        # M2_sf = 0.9 M2_ph  ==>  M2^D = 0.1 M2_ph
    M2_sf = lam_sf * w_sf ** 2
    N = 800
    lph = einstein(lam_ph, w_ph, N, one)
    lsf = einstein(lam_sf, w_sf, N, one)
    l1 = B.lam1(A.multiband_O2(lph + lsf, lph - lsf, N, np.array([1.0])))
    M2D = M2_ph - M2_sf
    M2G = M2_ph + 3 * M2_sf
    print(f"{w_sf:6.1f} {l1:10.6f} {M2D:9.4f} {l1/(G2*M2D):13.4f} {l1/(G2*M2G):11.6f}")

# control: pure phonon (M2^sf = 0) -- both formulas coincide and the ceiling must hold
N = 800
lph = einstein(lam_ph, w_ph, N, one)
l1 = B.lam1(A.multiband_O2(lph, lph, N, np.array([1.0])))
print(f"\ncontrol pure phonon: lam1={l1:.6f}, ceiling {G2*M2_ph:.6f}, ratio {l1/(G2*M2_ph):.6f}")
