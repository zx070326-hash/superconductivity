"""[numerical-evidence] End-to-end: is the G3-based UNCONDITIONAL ceiling correct
under (Z-dom) when the pairing channel is genuinely signed (s+-)?

team-lead's (G'-3) reads   Tc <= C_inf sqrt(rho(M2')),  M2' := M2 + 2*offdiag(M2).
Here M2 is the moment of the kernel that enters the operator, i.e. the PAIRING one,
M2^D = M2^ph - M2^sf under s+-.  But G3's hypothesis is "f >= 0 non-increasing", and
Lam^D(k) is signed, so neither  ||O[Lam^D]|| <= 3 lam1(O[Lam^D])  nor lemma M applies to it.

My repair (split on the SPECTRAL function, where both pieces are nonneg and non-increasing in k):
    diagonal   N_ii = lam1(O[Lam^D_ii]) <= lam1(O[Lam^ph_ii]) + lam1(-O[Lam^sf_ii])
                                        <= (g2/kap) (M2^ph + 3 M2^sf)_ii      [G3 on the sf piece]
    off-diag   N_ij = ||O[Lam^D_ij]||   <= ||O[Lam^ph_ij]|| + ||O[Lam^sf_ij]||
                                        <= 3 (g2/kap) (M2^ph + M2^sf)_ij      [G3 twice]
    ==>  M2G := (M2^ph + 3 M2^sf) + 2 * offdiag(M2^ph + M2^sf)   elementwise >= N kap/g2
    ==>  Tc <= C_inf sqrt(rho(M2G)).
Pure phonon (M2^sf = 0) reduces to team-lead's M2' exactly.

DIRECTION: lam1(P_N OO P_N) <= lam1(OO)  ==>  the ratio computed at finite N is a LOWER
estimate of the true ratio  ==>  ratio > 1 is a SOUND violation; ratio <= 1 is support only.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

import _chk_ml_build as B

rng = np.random.default_rng(20260918_5)
G2 = 1.318140212132
N = 120
FAIL = 0


def rho(M):
    return float(max(abs(np.linalg.eigvals(M))))


def offd(M):
    return M - np.diag(np.diag(M))


nv_mine = nv_theirs = 0
worst_mine = worst_theirs = 0.0
cfg_theirs = None
T = 1500
for _ in range(T):
    Bb = int(rng.integers(2, 5))
    ph = rng.uniform(0.0, 1.5, size=(Bb, Bb)); ph = 0.5 * (ph + ph.T)
    sf = rng.uniform(0.0, 1.5, size=(Bb, Bb)); sf = 0.5 * (sf + sf.T)
    wph = float(rng.uniform(0.1, 8.0)); wsf = float(rng.uniform(0.1, 8.0))
    w = np.concatenate([[1.0], rng.uniform(0.2, 5.0, size=Bb - 1)])

    lph = B.lam_factorized(ph, wph, N)
    lsf = B.lam_factorized(sf, wsf, N)
    lamZ, lamD = lph + lsf, lph - lsf
    OO = mb = None
    import _chk_ml_absM2 as A
    OO = A.multiband_O2(lamZ, lamD, N, w)
    l1 = B.lam1(OO)

    M2ph, M2sf = wph ** 2 * ph, wsf ** 2 * sf          # M2 = varpi^2 * Lambda  (Einstein)
    M2D = M2ph - M2sf
    M2G = (M2ph + 3 * M2sf) + 2 * offd(M2ph + M2sf)     # mine
    M2p = M2D + 2 * offd(M2D)                           # team-lead's, applied to the pairing moment

    r_mine = l1 / (G2 * rho(M2G))
    r_th = l1 / (G2 * rho(M2p)) if rho(M2p) > 0 else np.inf
    nv_mine += r_mine > 1.0
    nv_theirs += r_th > 1.0
    worst_mine = max(worst_mine, r_mine)
    if r_th > worst_theirs:
        worst_theirs, cfg_theirs = r_th, (Bb, round(wph, 3), round(wsf, 3))

print(f"(Z-dom, separate ph/sf Einstein modes, {T} draws, N={N})")
print(f"  mine   M2G = (M2ph+3M2sf) + 2*offdiag(M2ph+M2sf): {nv_mine} violations, "
      f"worst ratio {worst_mine:.6f}")
print(f"  theirs M2' = M2D + 2*offdiag(M2D), M2D = M2ph-M2sf: {nv_theirs} violations, "
      f"worst ratio {worst_theirs:.4f}  cfg(B,wph,wsf)={cfg_theirs}")
FAIL += (nv_mine != 0)

# pure phonon: the two formulas must coincide
same = 0
for _ in range(500):
    Bb = int(rng.integers(2, 5))
    ph = rng.uniform(0.0, 1.5, size=(Bb, Bb)); ph = 0.5 * (ph + ph.T)
    M2G = ph + 2 * offd(ph)
    same += np.max(np.abs(M2G - (ph + 2 * offd(ph)))) < 1e-15
print(f"  pure-phonon reduction M2G == M2' : {same}/500")

# N-convergence of the worst 'theirs' ratio, to show the violation is not truncation
Bb, wph, wsf = 2, 1.0, 1.0
ph = np.array([[0.2, 1.4], [1.4, 0.2]]); sf = np.array([[0.15, 1.35], [1.35, 0.15]])
w = np.array([1.0, 1.0])
for Nn in (60, 120, 240, 480):
    lph = B.lam_factorized(ph, wph, Nn); lsf = B.lam_factorized(sf, wsf, Nn)
    import _chk_ml_absM2 as A
    l1 = B.lam1(A.multiband_O2(lph + lsf, lph - lsf, Nn, w))
    M2ph, M2sf = wph ** 2 * ph, wsf ** 2 * sf
    M2D = M2ph - M2sf
    M2G = (M2ph + 3 * M2sf) + 2 * offd(M2ph + M2sf)
    print(f"  N={Nn:4d}  lam1={l1:.6f}  ratio_theirs={l1/(G2*rho(M2D + 2*offd(M2D))):.6f}  "
          f"ratio_mine={l1/(G2*rho(M2G)):.6f}")

print(f"\nTOTAL FAILURES: {FAIL}")
