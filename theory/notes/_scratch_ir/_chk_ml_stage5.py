"""Stage 5: viability of the two candidate proof routes for the non-factorized case.
Level: numerical-evidence."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chk_ml_build import (u_vec, brk, multiband_O, lam_from_measures, lam1)

N = 120
rng = np.random.default_rng(2026)

print("=" * 78)
print("(1) MATRIX LEMMA M:  Lt(k) <= Gh(k) elementwise  =>  lam1(OO[Lt]) <= lam1(OO[Gh]) ?")
print("    (both built from elementwise-nonneg matrix measures, same w)")
print("=" * 78)
bad = 0
worst = None
for trial in range(1500):
    B = int(rng.integers(2, 5))
    w = np.exp(rng.uniform(-1.2, 1.2, B))
    P1, P2 = [], []
    for _ in range(int(rng.integers(1, 4))):
        W = rng.random((B, B)); W = (W + W.T) / 2 * rng.uniform(0.05, 3.0)
        vp = float(np.exp(rng.uniform(np.log(0.05), np.log(20.0))))
        P1.append((W, vp))
        P2.append((W, vp))
    # Gh = Lt + extra nonneg measure  => Gh(k) >= Lt(k) elementwise for every k
    for _ in range(int(rng.integers(1, 3))):
        W = rng.random((B, B)); W = (W + W.T) / 2 * rng.uniform(0.01, 2.0)
        P2.append((W, float(np.exp(rng.uniform(np.log(0.05), np.log(20.0))))))
    l1 = lam1(multiband_O(lam_from_measures(P1, N), N, w=w))
    l2 = lam1(multiband_O(lam_from_measures(P2, N), N, w=w))
    if l2 < l1 - 1e-12:
        bad += 1
        if worst is None or (l1 - l2) > worst[0]:
            worst = (l1 - l2, B, l1, l2, [f"{v:.3f}" for _, v in P1])
print(f"   violations: {bad}/1500")
if worst:
    print(f"   worst deficit {worst[0]:.6e}  B={worst[1]}  lam1(Lt)={worst[2]:.6f} "
          f"lam1(Gh)={worst[3]:.6f}  varpis={worst[4]}")

print()
print("=" * 78)
print("(2) cone positivity:  psi in C (all bands nonneg & nonincreasing), d elementwise>=0")
print("    is  <psi, OO[d] psi>  >= 0 ?   (needed for a cone proof of matrix Lemma M)")
print("=" * 78)
bad2 = 0
minval = 1e9
for trial in range(2000):
    B = int(rng.integers(2, 5))
    w = np.exp(rng.uniform(-1.2, 1.2, B))
    Pd = []
    for _ in range(int(rng.integers(1, 3))):
        W = rng.random((B, B)); W = (W + W.T) / 2 * rng.uniform(0.05, 3.0)
        Pd.append((W, float(np.exp(rng.uniform(np.log(0.05), np.log(20.0))))))
    Od = multiband_O(lam_from_measures(Pd, N), N, w=w)
    # random cone vector: each band nonneg, nonincreasing
    psi = np.cumsum(rng.random((N, B))[::-1], axis=0)[::-1].copy()
    psi *= rng.random(B)[None, :]
    psi = psi.reshape(-1)
    nrm = psi @ psi
    val = (psi @ Od @ psi) / nrm
    minval = min(minval, val)
    if val < -1e-12:
        bad2 += 1
print(f"   violations: {bad2}/2000   min Rayleigh value = {minval:.6e}")

print()
print("=" * 78)
print("(3) Metzler in the step(x)band basis: coefficient of e_L (x) e_i in OO(e_L' (x) e_j)")
print("=" * 78)
Ns = 40
Bn = 2
w = np.array([1.0, 1.7])
P = [(np.array([[1.0, 0.4], [0.4, 0.6]]), 2.0)]
M = multiband_O(lam_from_measures(P, Ns), Ns, w=w).reshape(Ns, Bn, Ns, Bn)
neg = []
for Lp in range(1, Ns):
    for j in range(Bn):
        v = np.zeros((Ns, Bn)); v[:Lp, j] = 1.0          # e_{L'} (x) e_j
        out = np.einsum("nimj,mj->ni", M, v)
        for L in range(1, Ns):
            for i in range(Bn):
                coef = out[L - 1, i] - (out[L, i] if L < Ns else 0.0)
                if (L, i) != (Lp, j) and coef < -1e-12:
                    neg.append((coef, L, i, Lp, j))
neg.sort()
print(f"   negative off-diagonal step-basis coefficients: {len(neg)}")
for t in neg[:5]:
    print(f"     coef={t[0]: .6f}  (L={t[1]},i={t[2]}) <- (L'={t[3]},j={t[4]})")
same_band = [t for t in neg if t[2] == t[4]]
cross_band = [t for t in neg if t[2] != t[4]]
print(f"   of which same-band: {len(same_band)}, cross-band: {len(cross_band)}")
