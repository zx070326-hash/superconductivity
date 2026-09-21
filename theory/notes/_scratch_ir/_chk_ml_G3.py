"""[numerical-evidence, consistency only]  My own check of theorem G3, re-derived
independently before running anything.

G3.  f >= 0 non-increasing, lam1 := lam1(O[f]) > 0.  Then |lmin(O[f])| <= lam1 + 2 f(1) <= 3 lam1.

My re-derivation against THIS repo's builder convention
    O[f]_{nm} = u_n u_m ( f(|n-m|)[n!=m] + f(n+m+1) )  -  delta_{nm} 2 u_n^2 sum_{k=1}^{n} f(k)
i.e. O = -L1 + P,  L1 = diag(2/(2n+1) sum_{k<=n} f(k)),  P = L2+L3 elementwise >= 0.

 1. 0 <= L1_nn = (2/(2n+1)) sum_{k=1}^n f(k) <= 2n f(1)/(2n+1) < f(1).      [f non-incr]
 2. lam1 >= <e0,O e0> = O_00 = f(1)   (u_0 = 1, L1_00 = empty sum = 0, L2 has no diagonal).
 3. P elementwise >=0 symmetric ==> |<x,Px>| <= <|x|,P|x|> <= lmax(P)||x||^2 ==> ||P|| = lmax(P).
 4. P = O + L1, Weyl ==> lmax(P) <= lam1 + lmax(L1) <= lam1 + f(1).
 5. lmin(O) >= lmin(P) - lmax(L1) >= -lmax(P) - f(1) >= -(lam1 + 2 f(1)) >= -3 lam1  [by 2].

All five steps are elementary and hold for EVERY N (truncation-independent), so the
statement is analytic; the numbers below only look for a transcription slip.

DIRECTION: each printed quantity is compared against its own bound; a violation would be
a genuine FAIL (the bound is claimed for the truncated matrix too, step by step).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import importlib.util
import numpy as np

spec = importlib.util.spec_from_file_location(
    "B", os.path.join(os.path.dirname(os.path.abspath(__file__)), "_chk_ml_build.py"))
B = importlib.util.module_from_spec(spec); spec.loader.exec_module(B)

rng = np.random.default_rng(20260918)
fails = 0
N = 200
ks = np.arange(0, 2 * N + 2, dtype=float)


def kernels():
    yield "1/k^2", np.concatenate(([0.0], 1.0 / ks[1:] ** 2))
    for K in (1, 2, 5, 100):
        yield f"1_[k<={K}]", (ks <= K).astype(float)
    for g in (0.3, 1.0, 2.5, 4.0):
        yield f"k^-{g}", np.concatenate(([0.0], ks[1:] ** (-g)))
    for w in (0.15, 0.8, 4.0):
        yield f"Einstein w={w}", w ** 2 / (ks ** 2 + w ** 2)
    for _ in range(300):                     # random non-increasing f >= 0
        d = rng.exponential(size=len(ks))
        f = np.cumsum(d[::-1])[::-1]
        f = f / f[1]
        yield "random non-incr", f


worst = 0.0
worst_name = ""
n_seen = 0
for name, f in kernels():
    O = B.scalar_O(f, N)
    ev = np.linalg.eigvalsh(O)
    lam1, lmin = ev[-1], ev[0]
    f1 = f[1]
    L1 = 2.0 / (2 * np.arange(N) + 1.0) * np.concatenate(([0.0], np.cumsum(f[1:N])))
    P = O + np.diag(L1)
    n_seen += 1
    # step-by-step
    if L1.max() > f1 + 1e-12:
        fails += 1; print(f"  [{name}] step1 FAIL: max L1={L1.max():.6g} > f(1)={f1:.6g}")
    if lam1 < f1 - 1e-12:
        fails += 1; print(f"  [{name}] step2 FAIL: lam1={lam1:.6g} < f(1)={f1:.6g}")
    if P.min() < -1e-14:
        fails += 1; print(f"  [{name}] step3 FAIL: P has a negative entry {P.min():.3g}")
    pev = np.linalg.eigvalsh(P)
    if abs(pev).max() > pev[-1] + 1e-10:
        fails += 1; print(f"  [{name}] step3 FAIL: ||P|| != lmax(P)")
    if pev[-1] > lam1 + f1 + 1e-10:
        fails += 1; print(f"  [{name}] step4 FAIL")
    if abs(lmin) > lam1 + 2 * f1 + 1e-10:
        fails += 1; print(f"  [{name}] G3 FAIL: |lmin|={abs(lmin):.6g} > lam1+2f(1)={lam1 + 2 * f1:.6g}")
    r = abs(lmin) / lam1
    if r > worst:
        worst, worst_name = r, name

print(f"kernels tested: {n_seen} (N={N});  worst |lmin|/lam1 = {worst:.6f}  at {worst_name}")
print("  (lemma (G) conjectured sup 0.814262 at 1_[k<=1]; G3's constant 3 is far above it)")

# N-independence of the worst case (G3 is per-N, but check the ratio is not a truncation effect)
f = (ks <= 1).astype(float)
for Nn in (50, 100, 200, 400, 800):
    kk = np.arange(0, 2 * Nn + 2, dtype=float)
    ff = (kk <= 1).astype(float)
    ev = np.linalg.eigvalsh(B.scalar_O(ff, Nn))
    print(f"  1_[k<=1]: N={Nn:4d}  lam1={ev[-1]:.6f}  lmin={ev[0]:.6f}  ratio={abs(ev[0])/ev[-1]:.6f}")

# 7.2: N <= (g2/kappa) M2' elementwise with M2' = M2 + 2*offdiag(M2) ==> rho(N) <= rho(M2')
bad = 0
for _ in range(2000):
    Bb = rng.integers(2, 5)
    A = rng.exponential(size=(Bb, Bb)); M2 = (A + A.T) / 2
    M2p = M2 + 2 * (M2 - np.diag(np.diag(M2)))
    Nmat = np.diag(np.diag(M2)) + rng.uniform(0, 3, size=(Bb, Bb)) * (M2 - np.diag(np.diag(M2)))
    Nmat = (Nmat + Nmat.T) / 2
    if max(abs(np.linalg.eigvals(Nmat))) > max(abs(np.linalg.eigvals(M2p))) + 1e-10:
        bad += 1
print(f"rho monotone on elementwise-dominated nonneg matrices: {bad}/2000 violations")
fails += bad

# MgB2-style cost factor, my own arithmetic
for lbl, M in (("Golubov-like lam", np.array([[1.017, 0.213], [0.155, 0.448]])),
               ("identical coupling", np.array([[1.0, 1.0], [1.0, 1.0]])),
               ("pure offdiag", np.array([[0.0, 1.0], [1.0, 0.0]]))):
    Ms = (M + M.T) / 2
    Mp = Ms + 2 * (Ms - np.diag(np.diag(Ms)))
    r1 = max(abs(np.linalg.eigvals(Ms))); r2 = max(abs(np.linalg.eigvals(Mp)))
    print(f"  {lbl:20s}: rho={r1:.4f} rho'={r2:.4f} ceiling x{np.sqrt(r2/r1):.4f}")

print(f"\nTOTAL FAILURES: {fails}")
