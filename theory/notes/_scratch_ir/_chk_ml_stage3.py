"""Stage 3: (i) lam_min of the step kernels, (ii) sharp constant in Lemma N,
(iii) exact saturation condition, (iv) B=3..6 physical scan, (v) MgB2-like slack.
Level: numerical-evidence."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chk_ml_build import (brk, scalar_O, multiband_O, lam_factorized,
                           lam_from_measures, lam1, lammin, h_of, g2)

G2 = g2(1600)
np.set_printoptions(precision=6, suppress=True)

print("=" * 78)
print("(i) lam_min( O[1_{k<=K}] ) : is it >= -1 ?  (would give lam_min(O[F])>=-F(1))")
print("=" * 78)
N = 600
for K in (1, 2, 3, 5, 10, 30, 100, 300, 599):
    f = np.zeros(2 * N + 2)
    f[1:K + 1] = 1.0
    M = scalar_O(f, N)
    print(f"   K={K:4d}   lam_min={lammin(M): .6f}   lam1={lam1(M): .6f}   "
          f"|lmin|/lam1={abs(lammin(M))/lam1(M):.4f}")

print()
print("=" * 78)
print("(ii) Lemma N sharp constant:  sup |lam_min(O[F])| / lam1(O[F]) over class")
print("=" * 78)
N = 500
ks = np.arange(2 * N + 2)
cands = []
# Einstein sweep
for vp in np.exp(np.linspace(np.log(1e-3), np.log(1e4), 60)):
    cands.append((brk(ks, vp), f"Einstein vp={vp:.4g}"))
# pure 1/k^2 (the varpi->0 limit shape)
f = np.zeros(2 * N + 2); f[1:] = 1.0 / ks[1:] ** 2
cands.append((f, "1/k^2 (g(2) kernel)"))
# heavy mixtures skewed to small varpi
rng = np.random.default_rng(17)
for _ in range(400):
    F = np.zeros(2 * N + 2)
    for _ in range(int(rng.integers(1, 6))):
        vp = float(np.exp(rng.uniform(np.log(1e-3), np.log(1e3))))
        F += rng.uniform(0.01, 5.0) * brk(ks, vp)
    cands.append((F, "random mix"))
best = []
for F, tag in cands:
    M = scalar_O(F, N)
    a, b = lam1(M), lammin(M)
    if a <= 0:
        continue
    best.append((abs(b) / a, tag, a, b))
best.sort(reverse=True)
for t in best[:8]:
    print(f"   ratio={t[0]:.9f}   {t[1]:28s} lam1={t[2]: .6f} lmin={t[3]: .6f}")
print(f"   |lam_min(G(2))|/g(2) = {abs(lammin(scalar_O(f, N)))/lam1(scalar_O(f, N)):.9f}")

print()
print("=" * 78)
print("(iii) saturation:  R_i = (Lt w)_i / w_i constant  <=>  ratio == 1 exactly?")
print("=" * 78)
N = 240
rng = np.random.default_rng(23)
for trial in range(6):
    B = int(rng.integers(2, 5))
    Lt = rng.random((B, B)); Lt = (Lt + Lt.T) / 2
    ev, V = np.linalg.eigh(Lt)
    perron = np.abs(V[:, -1])
    vp = float(np.exp(rng.uniform(np.log(0.2), np.log(8.0))))
    h = h_of(vp, N)
    # w == Perron vector  =>  R_i == lam_max(Lt) constant
    r_sat = lam1(multiband_O(lam_factorized(Lt, vp, N), N, w=perron)) / (h * ev[-1])
    w_rand = np.exp(rng.uniform(-1, 1, B))
    r_gen = lam1(multiband_O(lam_factorized(Lt, vp, N), N, w=w_rand)) / (h * ev[-1])
    R_sat = (Lt @ perron) / perron
    R_gen = (Lt @ w_rand) / w_rand
    print(f"   B={B} vp={vp:.3f}  w=Perron: ratio={r_sat:.12f} (R={R_sat}) ")
    print(f"              w=random: ratio={r_gen:.12f} (R spread={R_gen.max()-R_gen.min():.4f})")

print()
print("=" * 78)
print("(iv) B=2..6, physical (elementwise>=0) random matrix measures, ceiling test")
print("     ratio = lam1(OO) / [ g(2) * lam_max(M2) ]     >1 would be a VIOLATION")
print("=" * 78)
N = 200
rng = np.random.default_rng(31)
worst = []
for trial in range(4000):
    B = int(rng.integers(2, 7))
    nm = int(rng.integers(1, 5))
    P = []
    M2 = np.zeros((B, B))
    for _ in range(nm):
        W = rng.random((B, B)) ** rng.uniform(0.3, 6.0)
        W = (W + W.T) / 2 * rng.uniform(0.02, 4.0)
        if rng.random() < 0.3:      # sometimes make it near rank-one / near-uniform
            v = rng.random(B)
            W = np.outer(v, v) * rng.uniform(0.1, 3.0)
        vp = float(np.exp(rng.uniform(np.log(0.02), np.log(50.0))))
        P.append((W, vp))
        M2 += W * vp ** 2
    w = np.exp(rng.uniform(-2.0, 2.0, B))
    M = multiband_O(lam_from_measures(P, N), N, w=w)
    bnd = G2 * np.linalg.eigvalsh(M2)[-1]
    worst.append((lam1(M) / bnd, B, nm, [f"{v:.3f}" for _, v in P]))
worst.sort(reverse=True)
print("   ratio       B nm  varpis")
for t in worst[:8]:
    print(f"   {t[0]:.9f}  {t[1]} {t[2]}   {t[3]}")
print(f"   max over {len(worst)}: {worst[0][0]:.9f}   median {np.median([t[0] for t in worst]):.4f}")

print()
print("=" * 78)
print("(v) MgB2-like two-band lambda matrix: how much slack does the ceiling leave?")
print("=" * 78)
# Golubov/Choi-type: lam_sigsig=1.017, lam_sigpi=0.213, lam_pisig=0.155, lam_pipi=0.448
# (lam_ij = N_j V_ij; detailed balance => symmetrized Lt_ij = sqrt(lam_ij lam_ji))
lam_mat = np.array([[1.017, 0.213], [0.155, 0.448]])
Lt = np.array([[lam_mat[0, 0], np.sqrt(lam_mat[0, 1] * lam_mat[1, 0])],
               [np.sqrt(lam_mat[0, 1] * lam_mat[1, 0]), lam_mat[1, 1]]])
w = np.array([1.0, np.sqrt(lam_mat[1, 0] / lam_mat[0, 1])])
R = (Lt @ w) / w
print(f"   Lt =\n{Lt}\n   w={w}  R={R}  (rows of the unsymmetrized lam: "
      f"{lam_mat.sum(axis=1)})")
print(f"   lam_max(Lt) = {np.linalg.eigvalsh(Lt)[-1]:.6f}   "
      f"total scalar lambda (Perron) vs max row sum {R.max():.6f}")
N = 400
# Einstein spectrum at 60 meV for both channels; find Tc from lam1(OO)=1
wE = 0.060  # eV
for label, shape in (("common Einstein 60meV", None),):
    lo, hi = 1e-4, 0.05
    for _ in range(60):
        T = 0.5 * (lo + hi)
        vp = wE / (2 * np.pi * T)
        v = lam1(multiband_O(lam_factorized(Lt, vp, N), N, w=w))
        if v > 1.0:
            lo = T
        else:
            hi = T
    Tc = 0.5 * (lo + hi)
    M2 = Lt * wE ** 2
    Tceil = (1 / (2 * np.pi)) * np.sqrt(G2 * np.linalg.eigvalsh(M2)[-1])
    K = 11604.5
    print(f"   {label}: Tc = {Tc*K:.2f} K ; ceiling = {Tceil*K:.2f} K ; "
          f"ratio Tc/ceiling = {Tc/Tceil:.4f}")
