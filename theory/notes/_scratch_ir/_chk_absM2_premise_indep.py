"""[numerical-evidence] Independent check of universality-lift's claim that
(G'-abs) as I stated it is missing a premise, plus lemmaG's Theorem G3.

I rebuild the multiband operator MYSELF from (ML1) as I re-derived it:

  O[(i,n),(j,m)] =  u_n u_m [ LD_ij f(|n-m|)(1-delta_nm) + LD_ij f(n+m+1) ]
                  - delta_ij delta_nm (2/(2n+1)) sum_{k=1}^{n} SZ_i(k)
                  - delta_nm u_n^2 [ delta_ij SZ_i(0) - LD_ij f(0) ]

with SZ_i(k) := sum_j LZ_ij f(k).  The last line is the residual of the
Z/pairing cancellation that is exact in the scalar case -- L(0) = diag(SZ) - LD.

Sanity gate: B=1 with LZ=LD must reproduce my validated scalar O[f].
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

G2 = 1.318140212132
rng = np.random.default_rng(20260918)


def scalar_O(f, N):
    n = np.arange(N); u = 1 / np.sqrt(2 * n + 1.)
    fv = np.array([f(k) for k in range(2 * N + 2)], dtype=float)
    T = fv[np.abs(n[:, None] - n[None, :])].copy(); np.fill_diagonal(T, 0.0)
    M = (T + fv[n[:, None] + n[None, :] + 1]) * np.outer(u, u)
    M[n, n] -= 2.0 * np.concatenate(([0.], np.cumsum(fv[1:N]))) / (2 * n + 1.)
    return M


def build(LZ, LD, varpi, N):
    B = len(LZ)
    n = np.arange(N); u = 1 / np.sqrt(2 * n + 1.)
    kk = np.arange(0, 2 * N + 2, dtype=float)
    fv = varpi ** 2 / (kk ** 2 + varpi ** 2)          # f(0) = 1
    T = fv[np.abs(n[:, None] - n[None, :])].copy(); np.fill_diagonal(T, 0.0)
    P = (T + fv[n[:, None] + n[None, :] + 1]) * np.outer(u, u)
    d = 2.0 * np.concatenate(([0.], np.cumsum(fv[1:N]))) / (2 * n + 1.)
    rsZ = LZ.sum(axis=1)
    M = np.zeros((B * N, B * N))
    for i in range(B):
        for j in range(B):
            blk = LD[i, j] * P
            if i == j:
                blk = blk - np.diag(rsZ[i] * d) - np.diag(u * u * (rsZ[i] - LD[i, i]))
            else:
                blk = blk + np.diag(u * u * LD[i, j])
            M[i * N:(i + 1) * N, j * N:(j + 1) * N] = blk
    return M


# ---------- gate 0: scalar reduction ----------
N = 200
w = 7.0
one = np.array([[1.0]])
A = build(one, one, w, N)
Bm = scalar_O(lambda k: w ** 2 / (k * k + w ** 2), N)
print(f"gate0  builder vs scalar O : max|diff| = {np.max(np.abs(A - Bm)):.3e}")
print(f"gate0  lam1(H(7))          = {np.linalg.eigvalsh(A)[-1]:.6f}"
      f"   (repo-implied 2.5284)")
print(f"gate0  ceiling/lam1        = "
      f"{G2 * w**2 / np.linalg.eigvalsh(A)[-1]:.4f}   (repo 25.5454)\n")


def ceiling_ratio(LZ, LD, varpi, N):
    lam1 = np.linalg.eigvalsh(build(LZ, LD, varpi, N))[-1]
    rho = np.abs(np.linalg.eigvals(np.abs(LD))).max()
    return lam1 / (G2 * rho * varpi ** 2)


# ---------- R1: ONE signed matrix drives BOTH channels ----------
worst, arg = -np.inf, None
nviol = 0
for _ in range(1500):
    a, b = rng.uniform(0, 2, 2)
    c = rng.uniform(-np.sqrt(a * b) * 2.0, np.sqrt(a * b) * 2.0)
    L = np.array([[a, c], [c, b]])
    vp = 10 ** rng.uniform(-1.0, 0.6)
    r = ceiling_ratio(L, L, vp, 90)
    nviol += r > 1.0
    if r > worst:
        worst, arg = r, (a, b, c, vp)
print(f"R1  LZ = LD = signed matrix : {nviol}/1500 violate the ceiling, "
      f"worst ratio = {worst:.3f}")
print(f"    worst at a={arg[0]:.3f} b={arg[1]:.3f} c={arg[2]:.3f} varpi={arg[3]:.4f}")
for NN in (90, 180, 360):
    print(f"      N={NN:4d} -> ratio {ceiling_ratio(np.array([[arg[0],arg[2]],[arg[2],arg[1]]]),
                                                    np.array([[arg[0],arg[2]],[arg[2],arg[1]]]),
                                                    arg[3], NN):.6f}")

# ---------- R2: physical s+-, LZ = Lph+Lsf, LD = Lph-Lsf ----------
worst2, nviol2 = -np.inf, 0
for _ in range(1500):
    Aph = rng.uniform(0, 1.5, (2, 2)); Aph = (Aph + Aph.T) / 2
    Asf = rng.uniform(0, 1.5, (2, 2)); Asf = (Asf + Asf.T) / 2
    vp = 10 ** rng.uniform(-1.0, 0.6)
    r = ceiling_ratio(Aph + Asf, Aph - Asf, vp, 90)
    nviol2 += r > 1.0
    worst2 = max(worst2, r)
print(f"\nR2  LZ=Lph+Lsf, LD=Lph-Lsf  : {nviol2}/1500 violate, worst ratio = {worst2:.4f}")

Asf = rng.uniform(0, 1.5, (2, 2)); Asf = (Asf + Asf.T) / 2
worst3, nviol3 = -np.inf, 0
for _ in range(600):
    Asf = rng.uniform(0, 1.5, (2, 2)); Asf = (Asf + Asf.T) / 2
    vp = 10 ** rng.uniform(-1.0, 0.6)
    r = ceiling_ratio(Asf, -Asf, vp, 90)          # pure spin fluctuation
    nviol3 += r > 1.0
    worst3 = max(worst3, r)
print(f"    pure spin-fluctuation (Lph=0): {nviol3}/600 violate, worst = {worst3:.4f}")

# ---------- (ML3-abs) as a pure quadratic-form statement ----------
bad = 0; mn = np.inf
for _ in range(20000):
    B = rng.integers(2, 7)
    Z = rng.uniform(0, 2, (B, B)); Z = (Z + Z.T) / 2
    D = Z * rng.uniform(-1, 1, (B, B)); D = (D + D.T) / 2
    D = np.clip(D, -Z, Z)                     # enforce |D| <= Z elementwise
    wv = rng.uniform(0.2, 3.0, B)
    R = (Z @ wv) / wv
    L = np.diag(R) - D
    e = np.linalg.eigvalsh((L + L.T) / 2).min() / max(np.linalg.norm(L), 1e-12)
    mn = min(mn, e); bad += e < -1e-9
print(f"\n(ML3-abs)  |D|<=Z elementwise : {bad}/20000 fail, min norm-eig = {mn:.3e}")
bad = 0; mn = np.inf
for _ in range(5000):                          # control: break dominance
    B = rng.integers(2, 7)
    Z = rng.uniform(0, 2, (B, B)); Z = (Z + Z.T) / 2
    D = Z * rng.uniform(-3, 3, (B, B)); D = (D + D.T) / 2
    wv = rng.uniform(0.2, 3.0, B)
    L = np.diag((Z @ wv) / wv) - D
    e = np.linalg.eigvalsh((L + L.T) / 2).min() / max(np.linalg.norm(L), 1e-12)
    mn = min(mn, e); bad += e < -1e-9
print(f"           control (no dominance): {bad}/5000 fail, min norm-eig = {mn:.3f}")

# ---------- lemmaG Theorem G3:  |lmin(O[f])| <= lam1 + 2 f(1) <= 3 lam1 ----------
print()
worstc = -np.inf
for _ in range(400):
    kind = rng.integers(0, 3)
    if kind == 0:
        K = int(rng.integers(1, 40)); f = lambda k, K=K: 1.0 * (1 <= k <= K)
    elif kind == 1:
        g = rng.uniform(0.3, 3.0); f = lambda k, g=g: 0.0 if k == 0 else k ** (-g)
    else:
        v = 10 ** rng.uniform(-1, 2); f = lambda k, v=v: v * v / (k * k + v * v)
    M = scalar_O(f, 140)
    ev = np.linalg.eigvalsh(M); l1, lm = ev[-1], ev[0]
    if l1 <= 0:
        continue
    assert abs(lm) <= l1 + 2 * f(1) + 1e-9, (kind, abs(lm), l1, f(1))
    worstc = max(worstc, abs(lm) / l1)
print(f"G3  |lmin| <= lam1 + 2f(1) : 0/400 violations; "
      f"worst observed |lmin|/lam1 = {worstc:.6f}  (G3 gives c=3)")
print(f"G3  step check  lam1 >= f(1)? and lam_max(L1) <= f(1)?")
for f, tag in ((lambda k: 0.0 if k == 0 else 1.0 / k ** 2, "1/k^2"),
               (lambda k: 1.0 * (1 <= k <= 1), "step k<=1"),
               (lambda k: 1.0 * (1 <= k <= 100), "step k<=100")):
    M = scalar_O(f, 400); n = np.arange(400)
    fv = np.array([f(k) for k in range(802)])
    L1d = 2.0 * np.concatenate(([0.], np.cumsum(fv[1:400]))) / (2 * n + 1.)
    l1 = np.linalg.eigvalsh(M)[-1]
    print(f"    {tag:12s} f(1)={f(1):.4f}  lam1={l1:.6f}  max(L1 diag)={L1d.max():.6f}"
          f"   -> lam1>=f(1): {l1 >= f(1)},  L1<=f(1): {L1d.max() <= f(1)}")
