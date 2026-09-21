"""rev_lm: independent float64 check of Lemma M (kernel monotonicity).
Self-built O[f] from the elementwise definition (P4_conjecture2.md sec 6.2):
  O[f]_{nm} = -delta_{nm} 2 u_n^2 sum_{k<=n} f(k) + (1-delta_{nm}) f(|n-m|) u_n u_m + f(n+m+1) u_n u_m,
  u_n = (2n+1)^{-1/2}.  No import from derivations/.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, time
rng = np.random.default_rng(20260912)

def O_of(f, N):
    """f: callable on positive integer array -> kernel values."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    A = np.abs(n[:, None] - n[None, :]).astype(float)
    S = (n[:, None] + n[None, :] + 1).astype(float)
    off = f(np.where(A == 0, 1.0, A))  # value at 0 never used
    np.fill_diagonal(off, 0.0)
    M = (off + f(S)) * np.outer(u, u)
    ks = np.arange(1, N, dtype=float)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * f(ks))))
    M[np.diag_indices(N)] -= cum * u * u
    return M

def top(M):
    w, V = np.linalg.eigh(M)
    v = V[:, -1]
    if v[0] < 0: v = -v
    return w[-1], w[-2], v

def in_cone(v, tol=1e-10):
    return (v.min() >= -tol) and (np.diff(v).max() <= tol)

def qf(M, v): return float(v @ M @ v)

brk = lambda w: (lambda k: w * w / (k * k + w * w))
inv2 = lambda k: 1.0 / (k * k)

# ---------------- 0. g(2) = lambda_1(O[1/k^2]) ----------------
print("== 0. g(2) from O[1/k^2] ==")
for N in (300, 600, 1200):
    l1, l2, v = top(O_of(inv2, N))
    print(f"  N={N:5d}  lam1={l1:.15f}  lam2={l2:.3e}  lam_min={np.linalg.eigvalsh(O_of(inv2,N))[0]:.6f}  cone={in_cone(v)}")
G_LO, G_HI = 1.318140211762671, 1.318140212136387
l1, _, _ = top(O_of(inv2, 1200))
print(f"  certified enclosure [{G_LO},{G_HI}]  N=1200 value inside: {G_LO-1e-9 <= l1 <= G_HI+1e-9}")

# ---------------- 1. Application 1 chain ----------------
print("\n== 1. Application 1: h(w) < g(2) w^2, explicit chain, N=600 ==")
N = 600
G2 = O_of(inv2, N)
g2 = top(G2)[0]
ok = True
for w in (0.05, 0.3, 1.0, 3.0, 10.0, 30.0):
    F = O_of(brk(w), N)
    Gk = w * w * G2                      # O[w^2/k^2] = w^2 O[1/k^2] (linearity)
    Gk2 = O_of(lambda k: w * w / (k * k), N)
    lin = np.abs(Gk - Gk2).max()
    D = O_of(lambda k: w**4 / (k * k * (k * k + w * w)), N)
    lin2 = np.abs((Gk - F) - D).max()    # O[G]-O[F] = O[d] incl. -L1 sign
    h, h2, v = top(F)
    cone = in_cone(v)
    qG = qf(Gk, v); qD = qf(D, v)
    ident = abs(qG - h - qD)
    chain = (h <= qG) and (qD > 0) and (qG <= g2 * w * w) and (h < g2 * w * w)
    ok &= chain and cone and lin < 1e-12 and lin2 < 1e-12
    print(f"  w={w:5.2f} h={h:.10f} <v,O[G]v>={qG:.10f} <v,O[d]v>={qD:.3e} g2w^2={g2*w*w:.10f} "
          f"gap={h-h2:.3e} cone={cone} lin={lin:.1e}/{lin2:.1e} ident={ident:.1e} chain={chain}")
print("  APP1 all OK:", ok)

# bonus: r(w)=h/w^2 monotone decreasing (kernel 1/(k^2+w^2) decreasing in w)
print("\n== 1b. r(w)=h/w^2 monotone (Lemma M with F=1/(k^2+w'^2), G=1/(k^2+w^2), w<w') ==")
ws = np.array([0.02, 0.05, 0.1, 0.3, 0.7, 1, 1.5, 2.4, 4, 8, 15, 30, 60])
rs = np.array([top(O_of(brk(w), N))[0] / w**2 for w in ws])
print("  r:", np.array2string(rs, precision=8))
print("  strictly decreasing:", bool(np.all(np.diff(rs) < 0)), " r(0.02)<g2:", rs[0] < g2)

# ---------------- 2. Application 2 ----------------
print("\n== 2. Application 2: k(P,T) <= h(w_rms), random multi-atom P (x = w^2 atoms) ==")
def app2(xs, ps):
    FT = O_of(lambda k: sum(p * x / (k * k + x) for x, p in zip(xs, ps)), N)
    xm = float(np.dot(ps, xs))
    G = O_of(lambda k: xm / (k * k + xm), N)
    kap, kap2, v = top(FT)
    hrms = top(G)[0]
    qG = qf(G, v); qD = qf(G - FT, v)
    # kernel-level Jensen check d>=0
    ks = np.arange(1, 2 * N, dtype=float)
    d = xm / (ks**2 + xm) - sum(p * x / (ks**2 + x) for x, p in zip(xs, ps))
    return kap, qG, hrms, qD, d.min(), in_cone(v), kap - kap2

worst = 1.0; viol = 0; ncase = 0; mind = 1.0; minqD = 1.0; conefail = 0
t0 = time.time()
cases = []
for i in range(150):
    m = rng.integers(2, 6)
    # log-uniform atoms in w in [0.03, 40]; near-degenerate variants
    if i % 5 == 0:
        w0 = np.exp(rng.uniform(np.log(0.05), np.log(30)))
        w = w0 * (1 + 1e-3 * rng.standard_normal(m))   # near point mass
    elif i % 5 == 1:
        w = np.exp(rng.uniform(np.log(0.03), np.log(40), m))
        ps = np.full(m, 1e-4); ps[0] = 1 - 1e-4 * (m - 1)  # near point mass in weight
    else:
        w = np.exp(rng.uniform(np.log(0.03), np.log(40), m))
    if i % 5 != 1:
        ps = rng.dirichlet(np.ones(m))
    xs = w**2
    kap, qG, hrms, qD, dmin, cone, gap = app2(xs, ps)
    ncase += 1
    if not (kap <= qG + 1e-13 and qG <= hrms + 1e-13 and kap <= hrms): viol += 1
    worst = min(worst, hrms - kap); mind = min(mind, dmin); minqD = min(minqD, qD)
    conefail += (not cone)
    cases.append((kap, hrms, qD, m))
print(f"  cases={ncase}  violations={viol}  min(h_rms-k)={worst:.3e}  min d(k)={mind:.2e}  min<v,O[d]v>={minqD:.2e}  "
      f"cone failures={conefail}  ({time.time()-t0:.1f}s)")
# extreme: two atoms w=0.01 and 100, weights 0.5/0.5 and 0.999/0.001
for xs, ps in (((0.01**2, 100.0**2), (0.5, 0.5)), ((0.01**2, 100.0**2), (0.999, 0.001)), ((0.3**2, 0.3**2 * 1.0000001), (0.5, 0.5))):
    kap, qG, hrms, qD, dmin, cone, gap = app2(np.array(xs), np.array(ps))
    print(f"  extreme xs={xs} ps={ps}: k={kap:.8f} <v,Gv>={qG:.8f} h_rms={hrms:.8f} <v,O[d]v>={qD:.2e} dmin={dmin:.1e} cone={cone}")

# ---------------- 3. Corollary A' for general d>=0 (non-monotone) with random cone vectors ----------------
print("\n== 3. <v,O[d]v> >= 0 for random d>=0 (non-monotone, incl. sparse) and random v in C ==")
mn = 1.0
for i in range(300):
    kind = i % 3
    if kind == 0:
        dk = rng.uniform(0, 1, 2 * N) / np.arange(1, 2 * N + 1)**rng.uniform(0, 2.5)
    elif kind == 1:
        dk = np.zeros(2 * N); dk[rng.integers(0, 2 * N, 5)] = rng.uniform(0, 1, 5)
    else:
        dk = rng.uniform(0, 1, 2 * N) * (rng.uniform(0, 1, 2 * N) < 0.05)
    D = O_of(lambda k: dk[k.astype(int) - 1], N)
    a = rng.exponential(1.0, N) * (rng.uniform(0, 1, N) < 0.1)  # sparse Abel coeffs
    v = np.cumsum(a[::-1])[::-1]  # v_n = sum_{L>n} a_L
    if v.max() == 0: continue
    v /= np.linalg.norm(v)
    mn = min(mn, qf(D, v))
print(f"  min <v,O[d]v> over 300 = {mn:.3e} (>= -1e-12: {mn >= -1e-12})")
# control: non-cone vector gives negative for d = delta_3
D3 = O_of(lambda k: (k == 3).astype(float), N)
e5 = np.zeros(N); e5[4] = 1
print(f"  control e_5 with d=delta_3: {qf(D3, e5):.4f} (negative expected)")

# ---------------- 4. Necessity of monotone F: non-monotone F <= G can violate lam1 ordering? ----------------
print("\n== 4. probe: is lam1(O[F]) <= lam1(O[G]) without F monotone? ==")
bad = 0
for i in range(200):
    Fk = rng.uniform(0, 1, 2 * N) / np.arange(1, 2 * N + 1)**2 * (rng.uniform(0,1,2*N) < 0.3)
    dk = rng.uniform(0, 1, 2 * N) / np.arange(1, 2 * N + 1)**2 * (rng.uniform(0,1,2*N) < 0.3)
    lF = top(O_of(lambda k: Fk[k.astype(int) - 1], N))[0]
    lG = top(O_of(lambda k: (Fk + dk)[k.astype(int) - 1], N))[0]
    if lG < lF - 1e-12: bad += 1
print(f"  non-monotone F: violations {bad}/200")
# explicit: F = delta_2 (non-monotone), G = delta_2 + c*delta_1?  and F = delta_3, G = delta_3+delta_1
for (Fk, Gk, name) in (((3,), (3, 1), "F=d3, G=d3+d1"), ((2,), (2, 1), "F=d2, G=d2+d1"), ((5,), (5, 2), "F=d5,G=d5+d2")):
    lF = top(O_of(lambda k: np.isin(k, Fk).astype(float), N))[0]
    lG = top(O_of(lambda k: np.isin(k, Gk).astype(float), N))[0]
    print(f"  {name}: lam1F={lF:.6f} lam1G={lG:.6f}  ordered={lG >= lF}")
print("DONE")
