"""[numerical-evidence] Consistency checks for the off-diagonal GF factor.

Claims being exercised (proofs are in lemma_G_offdiag_GF_factor_2026-09-18.md):

  (g)   g(v) := (1 - pi v / sinh(pi v)) / v^2  =  2 sum_{n>=1} (-1)^{n+1} / (v^2 + n^2),
        strictly decreasing, g(0+) = pi^2/6, and 2 C_{Einstein(v)} = v^2 g(v).
  (GFm) f = int Einstein(v) P(dv), P >= 0  ==>
        |lmin(O[f])| <= max( f(1), int v^2 g(v) P )  <= (pi^2/6) * int v^2 P  = (pi^2/6) M2/kappa.
  (MB)  block matrix N (N_ii = lam1(O[L_ii]), N_ij = ||O[L_ij]||) satisfies
        lam1(multiband O) <= lmax(N) <= rho(Chat)/kappa,
        Chat = g2 * diag(M2) + (pi^2/6) * offdiag(M2).
  (MgB2) cost factor sqrt(rho(M2 + (c-1) offdiag)/rho(M2)) for c = 3 vs c = pi^2/(6 g2).

Everything here is float64 / finite N: navigation + consistency only, never proof.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

G2_LO, G2_HI = 1.3181093293757127, 1.3181615784372800      # repo CERT-B (rational-exact)
G2 = 1.318140212132                                          # CERT-D midpoint, navigation only
PI26 = np.pi ** 2 / 6
rng = np.random.default_rng(20260918)


# ---------------------------------------------------------------- scalar operator (own builder)
def scalar_O(fk, N):
    """fk[k] = f(k) for k = 0..2N+1 (fk[0] unused).  O = -L1 + L2 + L3."""
    n = np.arange(N); u = 1 / np.sqrt(2 * n + 1.)
    T = fk[np.abs(n[:, None] - n[None, :])].copy(); np.fill_diagonal(T, 0.0)
    M = (T + fk[n[:, None] + n[None, :] + 1]) * np.outer(u, u)
    M[n, n] -= 2.0 * np.concatenate(([0.], np.cumsum(fk[1:N]))) / (2 * n + 1.)
    return M


def einstein(v, N):
    k = np.arange(0, 2 * N + 2, dtype=float)
    return v * v / (k * k + v * v)


def g_closed(v):
    return (1 - np.pi * v / np.sinh(np.pi * v)) / v ** 2


def g_series(v, M=200000):
    n = np.arange(1, M + 1, dtype=float)
    s = 2 * np.sum((-1.0) ** (n + 1) / (v * v + n * n))
    # alternating tail |remainder| <= first omitted term
    return s, 2.0 / (v * v + (M + 1) ** 2)


def C_f(fk, Mfft=1 << 16):
    """C_f = -inf_theta sum_{k>=1} f(k) cos k theta on an FFT grid (grid UNDER-estimates the sup)."""
    h = np.zeros(Mfft); L = min(len(fk) - 1, Mfft - 1)
    h[1:L + 1] = fk[1:L + 1]
    return float((-np.fft.rfft(h).real).max())


# ---------------------------------------------------------------- (g)
print("== (g) closed form vs alternating series vs FFT symbol ==")
for v in (0.05, 0.3, 0.45, 1.0, 3.0, 10.0):
    gc = g_closed(v); gs, tail = g_series(v)
    fk = einstein(v, 4000)
    twoC = 2 * C_f(fk)
    print(f"  v={v:5.2f}  g_closed={gc:.10f}  g_series={gs:.10f} (tail<{tail:.1e})"
          f"  2C_fft/v^2={twoC / v**2:.10f}   diff={abs(gc - gs):.1e}")
vs = np.linspace(1e-3, 50, 20001)
gv = g_closed(vs)
print(f"  g(0+) -> {g_closed(1e-3):.8f}  (pi^2/6 = {PI26:.8f});  monotone decreasing on grid: "
      f"{bool(np.all(np.diff(gv) < 0))};  max g = {gv.max():.8f}")

# threshold v*: g(v*) = g2  (where the off-diagonal factor becomes exactly 1)
from scipy.optimize import brentq
for tag, g2v in (("g2_lo", G2_LO), ("g2_hi", G2_HI), ("g2_mid", G2)):
    vstar = brentq(lambda v: g_closed(v) - g2v, 0.05, 2.0)
    print(f"  v* with g(v*) = {tag} = {g2v:.10f}:  v* = {vstar:.6f}   "
          f"(omega* = 2 pi T v* = {2*np.pi*vstar:.4f} k_B T)")
print(f"  worst-case factor c = (pi^2/6)/g2_lo = {PI26 / G2_LO:.6f},  sqrt = {np.sqrt(PI26 / G2_LO):.6f}")

# ---------------------------------------------------------------- (GFm) random Einstein mixtures
print("\n== (GFm) random positive mixtures of Einstein kernels, N=160 ==")
N = 160
worst_a = worst_b = worst_c = -np.inf
viol_a = viol_b = 0
for trial in range(600):
    m = int(rng.integers(1, 6))
    vv = 10 ** rng.uniform(-1.3, 1.0, m)          # v in [0.05, 10]
    pp = rng.uniform(0.1, 1.0, m)
    fk = sum(p * einstein(v, N) for p, v in zip(pp, vv))
    O = scalar_O(fk, N)
    ev = np.linalg.eigvalsh(O); lmin, l1 = ev[0], ev[-1]
    M2k = float(np.sum(pp * vv ** 2))             # M2/kappa = int v^2 P
    bound_res = max(fk[1], float(np.sum(pp * vv ** 2 * g_closed(vv))))   # resolved
    bound_sup = PI26 * M2k                                                 # sup form
    ra = abs(lmin) / bound_res;  rb = abs(lmin) / bound_sup
    rc = max(l1, abs(lmin)) / (max(G2, PI26) * M2k)
    viol_a += ra > 1; viol_b += rb > 1
    worst_a = max(worst_a, ra); worst_b = max(worst_b, rb); worst_c = max(worst_c, rc)
print(f"  |lmin| <= max(f(1), int v^2 g P):  {viol_a}/600 violations, worst ratio {worst_a:.4f}")
print(f"  |lmin| <= (pi^2/6) M2/kappa      :  {viol_b}/600 violations, worst ratio {worst_b:.4f}")
print(f"  ||O|| <= (pi^2/6) M2/kappa       :  worst ratio {worst_c:.4f}")
# single Einstein, v -> 0: where the sup bound is loosest
for v in (0.05, 0.2, 0.45, 1.0, 3.0):
    fk = einstein(v, 400); ev = np.linalg.eigvalsh(scalar_O(fk, 400))
    print(f"  Einstein v={v:4.2f}:  |lmin|/(v^2 g)={abs(ev[0])/(v*v*g_closed(v)):.4f}"
          f"   |lmin|/(g2 v^2)={abs(ev[0])/(G2*v*v):.4f}   lam1/(g2 v^2)={ev[-1]/(G2*v*v):.4f}"
          f"   GF factor c(v)=max(1,g/g2)={max(1.0, g_closed(v)/G2):.4f}")

# ---------------------------------------------------------------- (MB) multiband end-to-end
print("\n== (MB) multiband Einstein blocks, per-block v_ij, my own ML1 builder ==")


def build_mb(Lam, V, N):
    """Physical multiband operator from (ML1) as re-derived in the route note:
       O[(i,n),(j,m)] = u_n u_m Lam_ij [ f_ij(|n-m|)(1-d_nm) + f_ij(n+m+1) ]
                        - d_ij d_nm (2/(2n+1)) sum_{k<=n} sum_l Lam_il f_il(k)
                        - d_nm u_n^2 [ d_ij sum_l Lam_il f_il(0) - Lam_ij f_ij(0) ],   f_ij(0)=1."""
    B = len(Lam); n = np.arange(N); u = 1 / np.sqrt(2 * n + 1.)
    F = {(i, j): einstein(V[i, j], N) for i in range(B) for j in range(B)}
    rowS = [sum(Lam[i, l] * np.concatenate(([0.], np.cumsum(F[(i, l)][1:N]))) for l in range(B))
            for i in range(B)]
    rowZ0 = [sum(Lam[i, l] * F[(i, l)][0] for l in range(B)) for i in range(B)]
    M = np.zeros((B * N, B * N))
    for i in range(B):
        for j in range(B):
            fk = F[(i, j)]
            T = fk[np.abs(n[:, None] - n[None, :])].copy(); np.fill_diagonal(T, 0.0)
            blk = Lam[i, j] * (T + fk[n[:, None] + n[None, :] + 1]) * np.outer(u, u)
            if i == j:
                blk[n, n] -= 2.0 * rowS[i] / (2 * n + 1.) + u * u * (rowZ0[i] - Lam[i, i] * fk[0])
            else:
                blk[n, n] += u * u * Lam[i, j] * fk[0]
            M[i * N:(i + 1) * N, j * N:(j + 1) * N] = blk
    return M


# gate: B=1 reduces to scalar_O
N = 120
A = build_mb(np.array([[1.0]]), np.array([[7.0]]), N)
print(f"  gate B=1 vs scalar_O: max|diff| = {np.max(np.abs(A - scalar_O(einstein(7.0, N), N))):.2e}")

viol1 = viol2 = 0; worst2 = -np.inf; worstT = -np.inf
for trial in range(300):
    B = int(rng.integers(2, 4))
    Lam = rng.uniform(0.05, 1.5, (B, B)); Lam = (Lam + Lam.T) / 2
    V = 10 ** rng.uniform(-1.0, 0.6, (B, B)); V = (V + V.T) / 2     # v in [0.1, 4], symmetric
    O = build_mb(Lam, V, N)
    l1 = np.linalg.eigvalsh(O)[-1]
    Nm = np.zeros((B, B)); Chat = np.zeros((B, B))
    for i in range(B):
        for j in range(B):
            ev = np.linalg.eigvalsh(scalar_O(Lam[i, j] * einstein(V[i, j], N), N))
            Nm[i, j] = ev[-1] if i == j else max(ev[-1], abs(ev[0]))
            Chat[i, j] = (G2 if i == j else PI26) * Lam[i, j] * V[i, j] ** 2
    lN = np.linalg.eigvalsh(Nm)[-1]
    rC = np.abs(np.linalg.eigvals(Chat)).max()
    viol1 += l1 > lN * (1 + 1e-9); viol2 += lN > rC * (1 + 1e-9)
    worst2 = max(worst2, lN / rC); worstT = max(worstT, l1 / rC)
print(f"  (G'-1)+(ML3): lam1(O) <= lmax(N)      {viol1}/300 violations")
print(f"  (GF-M):       lmax(N) <= rho(Chat)    {viol2}/300 violations, worst lmax(N)/rho(Chat) = {worst2:.4f}")
print(f"  end-to-end:   worst lam1(O)/rho(Chat) = {worstT:.4f}")

# ---------------------------------------------------------------- (MgB2) cost table
print("\n== cost factor sqrt(rho(M2 + (c-1) offdiag)/rho(M2)) ==")


def cost(M2, c):
    Mp = M2 + (c - 1) * (M2 - np.diag(np.diag(M2)))
    r = lambda X: np.abs(np.linalg.eigvals(X)).max()
    return np.sqrt(r(Mp) / r(M2)), r(M2), r(Mp)


golubov = {"MgB2 Golubov, conv. B (geom. mean)": np.array([[1.017, np.sqrt(0.213 * 0.155)],
                                                             [np.sqrt(0.213 * 0.155), 0.448]]),
           "MgB2 Golubov, conv. A (arith. mean)": np.array([[1.017, 0.184], [0.184, 0.448]]),
           "all-equal [[1,1],[1,1]]": np.ones((2, 2)),
           "pure off-diagonal [[0,1],[1,0]]": np.array([[0., 1.], [1., 0.]])}
cGF = PI26 / G2_LO
for tag, M2 in golubov.items():
    s3, r0, r3 = cost(M2, 3.0); sg, _, rg = cost(M2, cGF)
    print(f"  {tag:38s} rho(M2)={r0:.4f}  c=3: rho'={r3:.4f} factor={s3:.4f}"
          f"   c=pi^2/(6 g2_lo)={cGF:.4f}: rho''={rg:.4f} factor={sg:.4f}")
