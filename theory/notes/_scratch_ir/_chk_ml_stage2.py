"""Stage 2: (a) sharpen the c>=0 boundary, (b) N-convergence of the worst cases,
(c) test the two structural sub-claims a proof would need.
Level: numerical-evidence."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chk_ml_build import (u_vec, brk, scalar_O, multiband_O, lam_factorized,
                           lam_from_measures, lam1, lammin, h_of, g2)

np.set_printoptions(precision=6, suppress=True)
G2 = g2(1600)


def pair_only(lam, N, w):
    """OO_pair = sum_{k>=1} O[delta_k] (x) Lt(k): same as OO but with the Z-channel
    row sums R_i replaced by the matrix Lt itself (i.e. Laplacian term dropped)."""
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
            OO[n, i, n, j] += u ** 2 * lam[0][i, j]
            csum = np.concatenate([[0.0], np.cumsum(lam[1:N, i, j])])
            OO[n, i, n, j] -= u ** 2 * (lam[0][i, j] + 2.0 * csum)
    return OO.reshape(N * B, N * B)


print("=" * 78)
print("(a) c>=0 boundary, full precision.  ratio = lam1(OO) / [h(varpi)*lam_max(Lt)]")
print("=" * 78)
N = 240
for vp in (0.5, 2.0, 8.0):
    h = h_of(vp, N)
    print(f"  varpi={vp}")
    for c in (0.0, 1e-6, 1e-3, 0.05, 0.2, 0.5, 1.0):
        for (a, b) in ((1.0, 1.0), (1.0, 0.3), (2.0, 0.5)):
            Lt = np.array([[a, c], [c, b]])
            lm = np.linalg.eigvalsh(Lt)[-1]
            for w2 in (1.0, 3.0):
                M = multiband_O(lam_factorized(Lt, vp, N), N, w=np.array([1.0, w2]))
                r = lam1(M) / (h * lm)
                print(f"    a={a} b={b} c={c:<8g} w2={w2}  ratio={r:.12f}")

print()
print("=" * 78)
print("(b) N-convergence of the worst PHYSICAL (elementwise>=0) cases")
print("=" * 78)
cases = []
rng = np.random.default_rng(7)
# reproduce the tightest general-measure cases: small varpi, one mode
for vp in (0.05, 0.08, 0.15, 0.4):
    W = np.array([[1.0, 0.9], [0.9, 1.0]])
    cases.append(([(W, vp)], np.array([1.0, 1.0]), f"1-mode vp={vp} near-rank-1 W"))
    W2 = np.array([[1.0, 0.02], [0.02, 1.0]])
    cases.append(([(W2, vp)], np.array([1.0, 4.0]), f"1-mode vp={vp} weak coupling w2=4"))
for P, w, tag in cases:
    M2 = sum(W * v ** 2 for W, v in P)
    bnd = G2 * np.linalg.eigvalsh(M2)[-1]
    row = []
    for Nn in (120, 240, 480, 960):
        M = multiband_O(lam_from_measures(P, Nn), Nn, w=w)
        row.append(lam1(M) / bnd)
    print(f"  {tag:36s} " + "  ".join(f"N={n}:{r:.9f}" for n, r in
                                      zip((120, 240, 480, 960), row)))

print()
print("=" * 78)
print("(c1) Lemma N candidate:  |lam_min(O[F])| <= lam1(O[F])  for F in the class?")
print("     (needed to turn lam1(O (x) Lt) into lam1(O)*lam_max(Lt))")
print("=" * 78)
N = 400
rng = np.random.default_rng(3)
worst = []
for trial in range(600):
    nm = int(rng.integers(1, 5))
    ks = np.arange(2 * N + 2)
    F = np.zeros(2 * N + 2)
    for _ in range(nm):
        vp = float(np.exp(rng.uniform(np.log(0.02), np.log(200.0))))
        F += rng.uniform(0.05, 3.0) * brk(ks, vp)
    F /= F[1]
    M = scalar_O(F, N)
    a, b = lam1(M), lammin(M)
    worst.append((abs(b) / a, a, b, nm))
worst.sort(reverse=True)
print("   |lmin|/lmax   lam1        lmin      nmodes")
for t in worst[:6]:
    print(f"   {t[0]:.6f}     {t[1]: .5f}  {t[2]: .5f}   {t[3]}")
print(f"   max over 600 random F in the class: {worst[0][0]:.6f}")

print()
print("=" * 78)
print("(c2) is the RELAXATION OO_pair (Laplacian dropped) already enough?")
print("     ratio_pair = lam1(OO_pair)/[h*lam_max(Lt)]   vs   ratio_full")
print("=" * 78)
N = 240
for vp in (0.5, 2.0, 8.0):
    h = h_of(vp, N)
    for (a, b, c) in ((1.0, 1.0, 0.9), (1.0, 1.0, 0.1), (2.0, 0.3, 0.7), (1.0, 1.0, 1.0)):
        Lt = np.array([[a, c], [c, b]])
        lm = np.linalg.eigvalsh(Lt)[-1]
        for w2 in (1.0, 3.0):
            w = np.array([1.0, w2])
            lamk = lam_factorized(Lt, vp, N)
            rf = lam1(multiband_O(lamk, N, w=w)) / (h * lm)
            rp = lam1(pair_only(lamk, N, w)) / (h * lm)
            print(f"  vp={vp:4.1f} a={a} b={b} c={c} w2={w2}: full={rf:.9f}  pair={rp:.9f}")

print()
print("=" * 78)
print("(c3) does the top eigenvector of OO lie in the cone  {psi_i >= 0, nonincreasing}?")
print("=" * 78)
N = 200
rng = np.random.default_rng(5)
viol_sign = viol_mono = 0
tot = 0
worst_desc = None
for trial in range(400):
    B = int(rng.integers(2, 4))
    nm = int(rng.integers(1, 3))
    P = []
    for _ in range(nm):
        W = rng.random((B, B)) ** rng.uniform(0.5, 4.0)
        W = (W + W.T) / 2 * rng.uniform(0.05, 3.0)
        P.append((W, float(np.exp(rng.uniform(np.log(0.1), np.log(20.0))))))
    w = np.exp(rng.uniform(-1.0, 1.0, B))
    M = multiband_O(lam_from_measures(P, N), N, w=w)
    ev, V = np.linalg.eigh(M)
    v = V[:, -1].reshape(N, B)
    if v.sum() < 0:
        v = -v
    tot += 1
    if v.min() < -1e-10:
        viol_sign += 1
    d = np.diff(v, axis=0)
    if d.max() > 1e-10:
        viol_mono += 1
        if worst_desc is None:
            worst_desc = (d.max(), [f"{x:.3f}" for _, x in P])
print(f"  sign violations (some psi_i(n) < 0): {viol_sign}/{tot}")
print(f"  monotonicity violations (psi_i not nonincreasing): {viol_mono}/{tot}")
if worst_desc:
    print(f"  first mono violation: max increment {worst_desc[0]:.3e}, varpis {worst_desc[1]}")
