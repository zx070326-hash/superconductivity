"""[numerical-evidence] My own re-derivation + check of theorem GF and GF-mix, so that
M2''_Z (c_GF = (pi^2/6)/g2) can be adopted as MY verified statement, not team-lead's.

GF (re-derived, analytic):  f >= 0 non-increasing, summable.  z_n := u_n x_n.
  <x,O x> = sum_{n!=m} f(|n-m|) z_n z_m + sum_{n,m} f(n+m+1) z_n z_m - 2 sum_n S(n) z_n^2,
  S(n) = sum_{k<=n} f(k).  Reflect  zt_{-1-n} := z_n  (n>=0).  Then on Z:
  Tt := sum_{n,m in Z} f(|n-m|) zt_n zt_m
      = 2[ sum_{n,m>=0} f(|n-m|) z_n z_m + sum_{n,m>=0} f(n+m+1) z_n z_m ]
  (four quadrants: (+,+),(-,-) give the Toeplitz block twice; (+,-),(-,+) give the
  Hankel block twice since |n-(-1-b)| = n+b+1).  Hence
  <x,O x> = Tt/2 - f(0)||z||^2 - 2 sum S(n) z_n^2.
  Parseval: Tt = (1/2pi) int F |zhat|^2,  F(th) = f(0) + 2 sum_{k>=1} f(k) cos k th,
  (1/2pi) int |zhat|^2 = ||zt||^2 = 2||z||^2  ==>  Tt/2 >= (f(0) - 2C_f)||z||^2,
  C_f := -min_th sum_{k>=1} f(k) cos k th  (>= 0 since the sum has zero mean in th).
  f(0) cancels EXACTLY.  With S(n) <= n f(1):
  <x,O x> >= -sum_n z_n^2 [2C_f + 2n f(1)] >= -M sum_n (2n+1) z_n^2 = -M||x||^2
  iff M(2n+1) = M + 2nM >= 2C_f + 2n f(1) for all n  <==  M := max(2C_f, f(1)).   QED

Einstein closed form (analytic):  sum_{k>=1} cos(k th)/(k^2+a^2)
  = (pi/2a) cosh(a(pi-th))/sinh(a pi) - 1/(2a^2)  on [0,2pi]  -> decreasing on [0,pi]
  ==> argmin at th = pi  ==>  2C_[[w]] = 1 - pi w / sinh(pi w) = w^2 g(w),
  g(w) = 2 sum_{n>=1} (-1)^{n+1}/(w^2+n^2) = 2 sum_j (4j-1)/[(w^2+(2j-1)^2)(w^2+(2j)^2)]
  paired terms positive & decreasing in w  ==>  g decreasing, sup g = g(0+) = 2 eta(2) = pi^2/6.

GF-mix (analytic):  Lam = int [[w]] dP, P >= 0, int w^2 dP = M2/kappa.  C_f is a sup of
  linear functionals of f  ==> convex/subadditive  ==> 2C_Lam <= int w^2 g(w) dP
  <= (pi^2/6) M2/kappa.  f(1) = int w^2/(1+w^2) dP <= M2/kappa.  Lemma M: lam1 <= g2 M2/kappa.
  ==> ||O[Lam]|| <= max(g2, pi^2/6) M2/kappa = (pi^2/6) M2/kappa,  c_GF = (pi^2/6)/g2.

M2''_Z := diag(M2ph + c_GF M2sf) + c_GF offdiag(M2ph + M2sf):
  offdiag N_ij <= ||O[Lph_ij]|| + ||O[Lsf_ij]|| <= (pi^2/6)(M2ph_ij + M2sf_ij)/kappa
  diag    N_ii <= lam1(O[Lph_ii]) + (-lmin(O[Lsf_ii])) <= [g2 M2ph_ii + (pi^2/6) M2sf_ii]/kappa

DIRECTION: lam1(P_N O P_N) <= lam1(O); lmin(P_N O P_N) >= lmin(O).  For the GF check
  -lmin is thus UNDER-estimated at finite N: a finite-N violation is sound, a pass is
  evidence only (N-trend shown).  For the ceiling, ratio > 1 is a sound violation.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh, eigvals
import _chk_ml_build as Bld
from _chk_ml_absM2 import multiband_O2

G2 = 1.318140212132
PI2_6 = np.pi ** 2 / 6
C_GF = PI2_6 / G2
rng = np.random.default_rng(20260918_3)
fail = 0


def C_f(fv, K=2 ** 14):
    """-min_theta sum_{k>=1} f(k) cos(k theta), dense theta grid (K points) -- a grid
    minimum OVER-estimates the true min => UNDER-estimates C_f => the GF bound built
    from it is slightly too SMALL => any pass here is conservative in the wrong direction
    only by the grid error; we also check against the Einstein closed form."""
    th = np.linspace(0, np.pi, K)
    k = np.arange(1, len(fv))
    s = np.cos(np.outer(th, k)) @ fv[1:]
    return -s.min()


# (1) Einstein closed form vs grid
print("(1) 2C_[[w]] vs 1 - pi w/sinh(pi w)")
for w in (0.1, 0.4589, 1.0, 3.0):
    fv = Bld.brk(np.arange(4001), w)
    num = 2 * C_f(fv); closed = 1 - np.pi * w / np.sinh(np.pi * w)
    print(f"    w={w:7.4f}  grid={num:.10f}  closed={closed:.10f}")
    fail += abs(num - closed) > 1e-6

# (2) g decreasing (mpmath, dps 50 -- float64 cannot resolve the w->0 cancellation)
import mpmath as mp
mp.mp.dps = 50
g = lambda w: (1 - mp.pi * w / mp.sinh(mp.pi * w)) / w ** 2
ws = [mp.mpf(10) ** e for e in np.linspace(-4, 2, 200)]
gv = [g(w) for w in ws]
mono = all(gv[i + 1] < gv[i] for i in range(len(gv) - 1))
print(f"(2) g strictly decreasing on [1e-4,1e2] (dps50): {mono};  g(1e-4)={float(gv[0]):.10f}"
      f"  pi^2/6={PI2_6:.10f}")
fail += not mono
wstar = mp.findroot(lambda w: g(w) - G2, 0.46)
print(f"    w* with g(w*)=g2 : {float(wstar):.8f}   2 pi w* = {float(2*mp.pi*wstar):.4f}"
      f"   c_GF = {C_GF:.7f}")

# (3) GF on 400 kernels, N=200, with N trend on the tightest
print("(3) GF  -lmin(O[f]) <= max(f(1), 2C_f)")
worst, warg = -np.inf, None
for _ in range(400):
    kind = rng.integers(0, 4)
    ks = np.arange(2 * 200 + 2, dtype=float)
    if kind == 0:
        K = int(rng.integers(1, 60)); fv = ((ks >= 1) & (ks <= K)).astype(float)
    elif kind == 1:
        gam = rng.uniform(0.5, 3.0); fv = np.where(ks == 0, 0.0, ks ** -gam)
    elif kind == 2:
        w = 10 ** rng.uniform(-1.5, 1.7); fv = Bld.brk(ks, w)
    else:
        a = 10 ** rng.uniform(-0.5, 1.0); fv = (1 + (ks / a) ** 2) ** -2.0
    fv = fv.astype(float)
    lm = eigvalsh(Bld.scalar_O(fv, 200))[0]
    bound = max(fv[1], 2 * C_f(fv))
    r = -lm / bound
    fail += r > 1 + 1e-9
    if r > worst:
        worst, warg = r, (kind, fv)
print(f"    400 kernels: worst  -lmin/bound = {worst:.6f}  (kind {warg[0]})")
for NN in (100, 200, 400, 800):
    ks = np.arange(2 * NN + 2, dtype=float)
    fv = warg[1][:2 * NN + 2] if len(warg[1]) >= 2 * NN + 2 else None
    if fv is None:
        break
    lm = eigvalsh(Bld.scalar_O(fv, NN))[0]
    print(f"      N={NN:4d}  -lmin={-lm:.8f}  bound={max(fv[1], 2*C_f(fv)):.8f}")

# (4) GF-mix on random positive mixtures: -lmin <= (pi^2/6) M2/kappa, lam1 <= g2 M2/kappa
print("(4) GF-mix on random positive Einstein mixtures (M2/kappa = sum c w^2)")
w1 = w2 = -np.inf
for _ in range(300):
    m = int(rng.integers(1, 5)); ks = np.arange(2 * 150 + 2, dtype=float)
    cs = rng.uniform(0.1, 2, m); wsv = 10 ** rng.uniform(-1.5, 1.2, m)
    fv = sum(c * Bld.brk(ks, w) for c, w in zip(cs, wsv)); m2 = float(np.sum(cs * wsv ** 2))
    ev = eigvalsh(Bld.scalar_O(fv, 150))
    w1 = max(w1, -ev[0] / (PI2_6 * m2)); w2 = max(w2, ev[-1] / (G2 * m2))
    fail += (-ev[0] > PI2_6 * m2 * (1 + 1e-9)) + (ev[-1] > G2 * m2 * (1 + 1e-9))
print(f"    worst -lmin/((pi^2/6)M2) = {w1:.6f}    worst lam1/(g2 M2) = {w2:.6f}")

# (5) M2''_Z end-to-end on the DIRECTED family r = M2sf/M2ph -> 1 (where |M2^D| died)
print("(5) M2''_Z ceiling, directed r = M2^sf/M2^ph, two frequencies, B in {1,2,3}")
off = lambda X: X - np.diag(np.diag(X)); dg = lambda X: np.diag(np.diag(X))
M2pp = lambda a, b: dg(a + C_GF * b) + C_GF * off(a + b)
M2G = lambda a, b: dg(a + 3 * b) + 3 * off(a + b)
rho = lambda M: float(np.max(np.abs(eigvals(M))))
N = 90
for r in (1.0, 0.999, 0.99, 0.9, 0.5):
    v = 0; vg = 0; wr = -np.inf; cnt = 0
    for _ in range(60):
        Bb = int(rng.integers(1, 4))
        Lph = rng.uniform(0.2, 1.5, (Bb, Bb)); Lph = (Lph + Lph.T) / 2
        wph = 10 ** rng.uniform(-1.0, 0.0); wsf = 10 ** rng.uniform(0.0, 0.8)
        Lsf = r * Lph * (wph / wsf) ** 2                 # M2sf = r * M2ph elementwise
        ks = np.arange(2 * N + 2)
        fp = Bld.brk(ks, wph)[:, None, None] * Lph[None]
        fs = Bld.brk(ks, wsf)[:, None, None] * Lsf[None]
        l1 = eigvalsh(multiband_O2(fp + fs, fp - fs, N, np.ones(Bb)))[-1]
        a = Lph * wph ** 2; b = Lsf * wsf ** 2
        rr = l1 / (G2 * rho(M2pp(a, b))); v += rr > 1; wr = max(wr, rr)
        vg += l1 / (G2 * rho(M2G(a, b))) > 1; cnt += 1
        fail += rr > 1
    print(f"    r={r:6.3f}: M2''_Z {v}/{cnt} viol, worst {wr:.4f};   M2^G(c=3) {vg}/{cnt} viol")

# single-point counterexample
wph, wsf, lph = 0.10, 5.0, 10.0; lsf = lph * wph ** 2 / wsf ** 2
ks = np.arange(2 * 400 + 2)
fp = (Bld.brk(ks, wph) * lph)[:, None, None]; fs = (Bld.brk(ks, wsf) * lsf)[:, None, None]
l1 = eigvalsh(multiband_O2(fp + fs, fp - fs, 400, np.ones(1)))[-1]
a = np.array([[lph * wph ** 2]]); b = np.array([[lsf * wsf ** 2]])
print(f"    counterexample: lam1={l1:.6f}  ratio M2''_Z={l1/(G2*rho(M2pp(a,b))):.4f}"
      f"  ratio M2^G={l1/(G2*rho(M2G(a,b))):.4f}")
print("\nTOTAL FAILURES:", int(fail))
