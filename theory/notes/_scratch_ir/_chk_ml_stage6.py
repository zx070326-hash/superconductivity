"""Stage 6: the relaxed operator OO_pair -- does IT have the monotonicity/cone
structure that the full OO lacks?  This decides whether the general (non-factorized)
case closes.  Level: numerical-evidence."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chk_ml_build import u_vec, brk, lam_from_measures, lam1, g2

G2 = g2(1600)
rng = np.random.default_rng(777)


def OO_pair(lam, N):
    """sum_{k>=1} O[delta_k] (x) Lt(k).  Independent of w (no Z row sums)."""
    B = lam.shape[1]
    u = u_vec(N)
    n = np.arange(N)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    uu = np.outer(u, u)
    OO = np.zeros((N, B, N, B))
    pair = np.where(D == 0, 0.0, 1.0)[:, :, None, None] * lam[D] + lam[S]
    OO += uu[:, None, :, None] * np.transpose(pair, (0, 2, 1, 3))
    for i in range(B):
        for j in range(B):
            csum = np.concatenate([[0.0], np.cumsum(lam[1:N, i, j])])
            OO[n, i, n, j] -= 2.0 * u ** 2 * csum
    return OO.reshape(N * B, N * B)


N = 120
print("=" * 78)
print("(1') MONOTONICITY of OO_pair:  d elementwise>=0  =>  lam1 increases?")
print("=" * 78)
bad = 0
worst = None
for trial in range(1500):
    B = int(rng.integers(2, 5))
    P1 = []
    for _ in range(int(rng.integers(1, 4))):
        W = rng.random((B, B)); W = (W + W.T) / 2 * rng.uniform(0.05, 3.0)
        P1.append((W, float(np.exp(rng.uniform(np.log(0.05), np.log(20.0))))))
    P2 = list(P1)
    for _ in range(int(rng.integers(1, 3))):
        W = rng.random((B, B)); W = (W + W.T) / 2 * rng.uniform(0.01, 2.0)
        P2.append((W, float(np.exp(rng.uniform(np.log(0.05), np.log(20.0))))))
    l1 = lam1(OO_pair(lam_from_measures(P1, N), N))
    l2 = lam1(OO_pair(lam_from_measures(P2, N), N))
    if l2 < l1 - 1e-11:
        bad += 1
        if worst is None or (l1 - l2) > worst[0]:
            worst = (l1 - l2, B, l1, l2)
print(f"   violations: {bad}/1500" + (f"   worst deficit {worst[0]:.3e} "
      f"(lam1 {worst[2]:.5f} -> {worst[3]:.5f}, B={worst[1]})" if worst else ""))

print()
print("=" * 78)
print("(2') cone positivity for OO_pair:  psi in C, d elementwise>=0 => <psi,OO_pair[d]psi> >=0")
print("     (this one MUST hold: Lemma A gives <psi_i, O[delta_k] psi_j> >= 0 termwise)")
print("=" * 78)
bad2 = 0
minval = 1e9
for trial in range(2000):
    B = int(rng.integers(2, 5))
    Pd = []
    for _ in range(int(rng.integers(1, 3))):
        W = rng.random((B, B)); W = (W + W.T) / 2 * rng.uniform(0.05, 3.0)
        Pd.append((W, float(np.exp(rng.uniform(np.log(0.05), np.log(20.0))))))
    Od = OO_pair(lam_from_measures(Pd, N), N)
    psi = np.cumsum(rng.random((N, B))[::-1], axis=0)[::-1].copy()
    psi *= rng.random(B)[None, :]
    psi = psi.reshape(-1)
    val = (psi @ Od @ psi) / (psi @ psi)
    minval = min(minval, val)
    if val < -1e-12:
        bad2 += 1
print(f"   violations: {bad2}/2000   min Rayleigh = {minval:.6e}")

print()
print("=" * 78)
print("(3') is the TOP EIGENVECTOR of OO_pair in the cone C (all bands >=0, nonincreasing)?")
print("=" * 78)
vs = vm = 0
tot = 0
for trial in range(400):
    B = int(rng.integers(2, 5))
    P = []
    for _ in range(int(rng.integers(1, 3))):
        W = rng.random((B, B)); W = (W + W.T) / 2 * rng.uniform(0.05, 3.0)
        P.append((W, float(np.exp(rng.uniform(np.log(0.1), np.log(20.0))))))
    M = OO_pair(lam_from_measures(P, N), N)
    ev, V = np.linalg.eigh(M)
    v = V[:, -1].reshape(N, B)
    if v.sum() < 0:
        v = -v
    tot += 1
    if v.min() < -1e-9:
        vs += 1
    if np.diff(v, axis=0).max() > 1e-9:
        vm += 1
print(f"   sign violations {vs}/{tot}   monotonicity violations {vm}/{tot}")

print()
print("=" * 78)
print("(4') direct: lam1(OO_pair[Lt]) <= g(2)*lam_max(M2) ?  [the general-case ceiling]")
print("=" * 78)
worst = []
for trial in range(2000):
    B = int(rng.integers(2, 6))
    P = []
    M2 = np.zeros((B, B))
    for _ in range(int(rng.integers(1, 5))):
        W = rng.random((B, B)) ** rng.uniform(0.3, 6.0)
        W = (W + W.T) / 2 * rng.uniform(0.02, 4.0)
        vp = float(np.exp(rng.uniform(np.log(0.02), np.log(50.0))))
        P.append((W, vp)); M2 += W * vp ** 2
    l = lam1(OO_pair(lam_from_measures(P, N), N))
    worst.append(l / (G2 * np.linalg.eigvalsh(M2)[-1]))
worst.sort(reverse=True)
print(f"   max ratio {worst[0]:.9f}  (2nd {worst[1]:.9f})   "
      f"{'OK' if worst[0] < 1 + 1e-9 else '*** VIOLATION ***'}")
