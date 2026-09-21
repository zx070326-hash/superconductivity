"""[numerical-evidence] Which moment should the (Z-dom) ceiling be written on?

Candidates for the matrix whose spectral radius appears under the square root:
  (a) |M2^D| = |M2^ph - M2^sf|          -- what the current notes say.  FALSE: zero ceiling
                                            at positive lam1 (see _chk_ml_G3_counterex.py).
  (b) M2G = (M2^ph+3M2^sf) + 2 offdiag(M2^ph+M2^sf)   -- my G3-route repair, safe but loose.
  (c) M2^+ := M2^ph  = moment of the ATTRACTIVE (positive) part of a2F^D alone.
      Reason it should work in the CONE route: Lam^D(k) <= Lam^ph(k) <= M2^ph/k^2 elementwise,
      so (S3)/(S4) can take Ghat = M2^ph/k^2; the repulsive piece only lowers lam1 twice over
      (repulsion in Delta, extra mass renormalisation in Z).
This scan asks whether (c) is violated anywhere.  If it survives, the ceiling under (Z-dom)
is the SAME constant as the phonon-only ceiling, with M2 read off the attractive part only.

DIRECTION: lam1(P_N OO P_N) <= lam1(OO) ==> ratio at finite N is a LOWER estimate
==> ratio > 1 is a sound violation; ratio <= 1 is support, not proof.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import _chk_ml_build as B
import _chk_ml_absM2 as A

rng = np.random.default_rng(20260918_7)
G2 = 1.318140212132
N = 150


def rho(M):
    return float(max(abs(np.linalg.eigvals(M))))


def offd(M):
    return M - np.diag(np.diag(M))


nv = 0
worst = 0.0
cfg = None
T = 2000
for _ in range(T):
    Bb = int(rng.integers(1, 5))
    ph = rng.uniform(0.0, 1.5, size=(Bb, Bb)); ph = 0.5 * (ph + ph.T)
    sf = rng.uniform(0.0, 1.5, size=(Bb, Bb)); sf = 0.5 * (sf + sf.T)
    wph = float(rng.uniform(0.05, 8.0)); wsf = float(rng.uniform(0.05, 8.0))
    w = np.concatenate([[1.0], rng.uniform(0.2, 5.0, size=Bb - 1)]) if Bb > 1 else np.array([1.0])
    lph = B.lam_factorized(ph, wph, N); lsf = B.lam_factorized(sf, wsf, N)
    l1 = B.lam1(A.multiband_O2(lph + lsf, lph - lsf, N, w))
    M2ph = wph ** 2 * ph
    if rho(M2ph) <= 0:
        continue
    r = l1 / (G2 * rho(M2ph))
    nv += r > 1.0
    if r > worst:
        worst, cfg = r, (Bb, round(wph, 3), round(wsf, 3))
print(f"(c) ceiling on the ATTRACTIVE moment only: {nv}/{T} violations, "
      f"worst ratio {worst:.6f} at (B,wph,wsf)={cfg}")

# adversarial: sf tuned to be as unhelpful as possible (tiny sf, or sf at extreme frequencies)
nv2 = 0; worst2 = 0.0
for _ in range(2000):
    Bb = int(rng.integers(1, 4))
    ph = rng.uniform(0.0, 2.0, size=(Bb, Bb)); ph = 0.5 * (ph + ph.T)
    scale = 10.0 ** rng.uniform(-3, 0.3)
    sf = scale * ph                                  # proportional sf, worst for cancellation
    wph = 10.0 ** rng.uniform(-1.2, 0.9); wsf = 10.0 ** rng.uniform(-1.2, 1.2)
    w = np.concatenate([[1.0], rng.uniform(0.2, 5.0, size=Bb - 1)]) if Bb > 1 else np.array([1.0])
    lph = B.lam_factorized(ph, wph, N); lsf = B.lam_factorized(sf, wsf, N)
    l1 = B.lam1(A.multiband_O2(lph + lsf, lph - lsf, N, w))
    M2ph = wph ** 2 * ph
    if rho(M2ph) <= 0:
        continue
    r = l1 / (G2 * rho(M2ph))
    nv2 += r > 1.0
    worst2 = max(worst2, r)
print(f"    adversarial sf proportional to ph: {nv2}/2000 violations, worst {worst2:.6f}")

# the counterexample of _chk_ml_G3_counterex.py, scored against (c)
one = np.ones((1, 1))
lam_ph, w_ph, w_sf = 10.0, 0.10, 5.0
M2_ph = lam_ph * w_ph ** 2
lam_sf = M2_ph / w_sf ** 2
for Nn in (200, 800, 1600):
    lph = B.lam_factorized(one * lam_ph, w_ph, Nn)
    lsf = B.lam_factorized(one * lam_sf, w_sf, Nn)
    l1 = B.lam1(A.multiband_O2(lph + lsf, lph - lsf, Nn, np.array([1.0])))
    print(f"    counterexample N={Nn}: lam1={l1:.6f}, (c) ceiling={G2*M2_ph:.6f}, "
          f"ratio={l1/(G2*M2_ph):.6f}")
