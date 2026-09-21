"""[numerical-evidence] Where do the (G'-abs) violations live, and does the
(ML3-abs) premise |LD| <= LZ elementwise alone suffice?  Plus lemmaG's budget
claim that the degradation factor is 1+(c-1)s, not c."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
exec(open("theory/notes/_scratch_ir/_chk_absM2_premise_indep.py").read().split(
     "# ---------- gate 0")[0].split('"""')[2])

rng = np.random.default_rng(7)
N = 90

# (1) do violations require a NEGATIVE Z row sum (i.e. unphysical mass renorm)?
rows = []
for _ in range(3000):
    a, b = rng.uniform(0, 2, 2)
    c = rng.uniform(-2 * np.sqrt(a * b), 2 * np.sqrt(a * b))
    L = np.array([[a, c], [c, b]]); vp = 10 ** rng.uniform(-1.0, 0.6)
    r = ceiling_ratio(L, L, vp, N) if False else (
        np.linalg.eigvalsh(build(L, L, vp, N))[-1]
        / (G2 * np.abs(np.linalg.eigvals(np.abs(L))).max() * vp ** 2))
    rows.append((r, L.sum(axis=1).min()))
R, MINROW = map(np.array, zip(*rows))
print(f"(1) violations vs sign of min_i (LZ row sum):")
for tag, m in (("row sum >= 0", MINROW >= 0), ("row sum <  0", MINROW < 0)):
    print(f"    {tag}: {int((R[m] > 1).sum()):4d}/{int(m.sum()):4d} violate, "
          f"worst {R[m].max():8.3f}")

# (2) does |LD| <= LZ elementwise (LZ nonneg) alone kill it?
worst, nv = -np.inf, 0
for _ in range(3000):
    Z = rng.uniform(0, 2, (2, 2)); Z = (Z + Z.T) / 2
    D = np.clip(Z * rng.uniform(-1.6, 1.6, (2, 2)), -Z, Z); D = (D + D.T) / 2
    D = np.clip(D, -Z, Z)
    vp = 10 ** rng.uniform(-1.0, 0.6)
    r = (np.linalg.eigvalsh(build(Z, D, vp, N))[-1]
         / (G2 * np.abs(np.linalg.eigvals(np.abs(D))).max() * vp ** 2))
    nv += r > 1; worst = max(worst, r)
print(f"\n(2) premise |LD|<=LZ elementwise, LZ>=0 : {nv}/3000 violate, worst {worst:.4f}")

# (3) lemmaG budget: degradation = 1 + (c-1)*s, s = off-diagonal Perron share
print("\n(3) budget 1+(c-1)s  (c=3 from G3):")
for M2 in (np.array([[1.0, 0.0], [0.0, 1.0]]),
           np.array([[1.0, 0.3], [0.3, 1.0]]),
           np.array([[0.5, 1.2], [1.2, 0.5]]),
           np.array([[1.0, 1.0], [1.0, 1.0]])):
    ev, V = np.linalg.eigh(M2); rho = ev[-1]; w = np.abs(V[:, -1])
    s = max(((M2 - np.diag(np.diag(M2))) @ w / w / rho))
    print(f"    M2={M2.tolist()}  rho={rho:.4f}  s={s:.4f}"
          f"  -> factor {1 + 2 * s:.4f}  (sqrt {np.sqrt(1 + 2 * s):.4f})")
