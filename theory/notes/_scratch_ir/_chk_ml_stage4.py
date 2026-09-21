"""Stage 4: (a) mechanism of the c<0 counterexample, (b) is elementwise-nonneg
NECESSARY or does R_i>=0 suffice, (c) Coulomb mu* matrix.
Level: numerical-evidence."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chk_ml_build import (u_vec, brk, multiband_O, lam_factorized,
                           lam_from_measures, lam1, lammin, h_of, g2)

G2 = g2(1600)
N = 240
np.set_printoptions(precision=6, suppress=True)

print("=" * 78)
print("(a) mechanism of the c<0 violation: what happens to R_i (the Z-channel)?")
print("=" * 78)
for (a, b, c, w2) in ((0.2, 0.2, -3.0, 5.0), (1.0, 1.0, -0.5, 1.0),
                      (1.0, 1.0, -0.9, 1.0), (1.0, 1.0, -0.2, 1.0)):
    Lt = np.array([[a, c], [c, b]])
    w = np.array([1.0, w2])
    R = (Lt @ w) / w
    ev = np.linalg.eigvalsh(Lt)
    vp = 0.5
    M = multiband_O(lam_factorized(Lt, vp, N), N, w=w)
    l1 = lam1(M)
    _, V = np.linalg.eigh(M)
    psi = V[:, -1].reshape(N, 2)
    bandsign = np.sign(psi[0])
    r = l1 / (h_of(vp, N) * ev[-1]) if ev[-1] > 0 else float('nan')
    print(f"   a={a} b={b} c={c:5.1f} w2={w2}: R={R}  lam_max(Lt)={ev[-1]: .4f}"
          f"  ratio={r:.4f}  top-vec band signs at n=0: {bandsign}")

print()
print("=" * 78)
print("(b) is elementwise>=0 NECESSARY?  scan c<0 but REQUIRE R_i >= 0 (Z_i >= 1)")
print("=" * 78)
worst = []
for a in (0.3, 1.0, 3.0):
    for b in (0.3, 1.0, 3.0):
        for w2 in (0.5, 1.0, 2.0):
            for c in -np.linspace(0.005, 3.0, 40):
                Lt = np.array([[a, c], [c, b]])
                w = np.array([1.0, w2])
                R = (Lt @ w) / w
                if R.min() < 0:
                    continue          # Z_i < 1 -- excluded by hand
                ev = np.linalg.eigvalsh(Lt)
                if ev[-1] <= 0:
                    continue
                for vp in (0.3, 1.0, 4.0):
                    M = multiband_O(lam_factorized(Lt, vp, N), N, w=w)
                    worst.append((lam1(M) / (h_of(vp, N) * ev[-1]), a, b, c, w2, vp,
                                  R.min()))
worst.sort(reverse=True)
print(f"   {len(worst)} configs with c<0 and R_i>=0")
print("   ratio       a    b     c      w2   varpi   min R_i")
for t in worst[:8]:
    print(f"   {t[0]:.6f}  {t[1]:4.1f} {t[2]:4.1f} {t[3]:7.3f} {t[4]:4.1f} {t[5]:5.1f}"
          f"  {t[6]:.4f}")

print()
print("=" * 78)
print("(c) Coulomb mu* matrix:  OO_mu = OO_ph - mu_ij * K_c    (K_c = 2 u_c u_c^T)")
print("    does the mu*=0 ceiling survive?  ratio = lam1 / [g2*lam_max(M2)]")
print("=" * 78)


def coulomb_block(mu, N, Nc):
    B = mu.shape[0]
    u = u_vec(N)
    uc = u.copy()
    uc[Nc:] = 0.0
    Kc = 2.0 * np.outer(uc, uc)
    out = np.zeros((N, B, N, B))
    out += Kc[:, None, :, None] * mu[None, :, None, :]
    return out.reshape(N * B, N * B)


rng = np.random.default_rng(41)
res = {"psd": [], "nonneg_nonpsd": [], "mixed_sign": []}
for trial in range(1200):
    B = 2
    Lt = rng.random((B, B)); Lt = (Lt + Lt.T) / 2 * rng.uniform(0.3, 3.0)
    w = np.exp(rng.uniform(-1, 1, B))
    vp = float(np.exp(rng.uniform(np.log(0.3), np.log(15.0))))
    kind = ["psd", "nonneg_nonpsd", "mixed_sign"][trial % 3]
    if kind == "psd":
        v = rng.random(B); mu = np.outer(v, v) * rng.uniform(0.05, 0.4) + \
            np.diag(rng.random(B) * 0.2)
    elif kind == "nonneg_nonpsd":
        d = rng.random(B) * 0.05
        off = rng.uniform(0.15, 0.5)
        mu = np.array([[d[0], off], [off, d[1]]])
    else:
        mu = np.array([[rng.uniform(0.05, 0.3), -rng.uniform(0.05, 0.5)],
                       [-rng.uniform(0.05, 0.5), rng.uniform(0.05, 0.3)]])
        mu = (mu + mu.T) / 2
    Nc = int(rng.integers(3, 40))
    M = multiband_O(lam_factorized(Lt, vp, N), N, w=w) - coulomb_block(mu, N, Nc)
    M2 = Lt * vp ** 2
    bnd = G2 * np.linalg.eigvalsh(M2)[-1]
    res[kind].append((lam1(M) / bnd, float(np.linalg.eigvalsh(mu)[0])))
for kind, vals in res.items():
    mx = max(v[0] for v in vals)
    mineig = min(v[1] for v in vals)
    print(f"   mu* {kind:16s}: n={len(vals):4d}  max ratio={mx:.6f}  "
          f"min lam_min(mu*)={mineig:+.4f}   {'OK' if mx < 1 else '*** VIOLATION ***'}")
