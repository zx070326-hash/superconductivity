"""rev_cu_5: (a) deeper N convergence at wbar=51.2/13.5; (b) secular engine vs dense
eigh at the near-degenerate counterexample; (c) independent mutilde_crit; (d) the
minimum of mutilde_crit/mutilde_* over a grid (how close does the theorem come to
claiming a false point?); (e) Lemma-1 hypotheses: c2, simplicity of lam_1(A), psi_1 in C."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import RCell, uc, myH, myOf

import sys
if "skipA" in sys.argv:
    pass
print("(a) deeper N-convergence, the two large-N_c anchors")
for (w, Nc) in ((13.5, 135), (51.2, 153)):
    prev = None
    for NN in (400, 800, 1200, 1600):
        t0=time.time(); C = RCell(w, Nc, NN)
        d = "" if prev is None else f"  d(mut_*)={C.thr-prev:+.2e}"
        print(f"   w={w} Nc={Nc} N={NN:>5d}: mut_*={C.thr:.7f} gam={C.gam:.7f} "
              f"q1={C.q1:.7f} m={C.m:.7f}{d}  ({time.time()-t0:.0f}s)")
        prev = C.thr
    print()

print("(b) secular reduction vs dense eigh at the near-degenerate point (0.74, N_c=4)")
NN = 600
C = RCell(0.74, 4, NN)
print(f"   mu1={C.mu1:.10e} mu2={C.mu2:.10e} mu3={C.mu3:.10e} gam={C.gam:.6f}")
print(f"   c1={C.c1:.6f} c2={C.c2:.4e} c3={abs(C.c[-3]):.4e}")
u = uc(NN, 4)
for mut in (0.13, 0.26, 0.320, 0.325, 0.330, 0.335, 0.340):
    e, V = np.linalg.eigh(C.A - 2*mut*np.outer(u, u))
    v = V[:, -1]
    q_dense = float(v @ (C.Of @ v))
    # secular root by bisection on (mu2, mu1)
    c2 = C.c**2
    f = lambda l: 2*mut*float(np.sum(c2/(C.e - l))) - 1.0
    lo, hi = C.mu2 + 1e-18, C.mu1 - 1e-18
    for _ in range(300):
        mid = 0.5*(lo+hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    lam_s = 0.5*(lo+hi)
    y = C.c/(C.e - lam_s); y = y/np.linalg.norm(y)
    q_sec = float(y @ (C.G @ y))
    print(f"   mut={mut:.3f}: lam_eigh={e[-1]:.12e} lam_sec={lam_s:.12e} "
          f"|d|={abs(e[-1]-lam_s):.2e}  q_dense={q_dense:+.4e} q_sec={q_sec:+.4e}"
          f"  lam1(M)-mu2={e[-1]-C.mu2:.3e}")
print(f"   mutilde_*(0.74,4) at N={NN}: {C.thr:.5f}")

print("\n(c) independent mutilde_crit(0.74,4) by DENSE eigh bisection")
def q_dense_of(C, Nc, mut):
    u = uc(C.N, Nc)
    e, V = np.linalg.eigh(C.A - 2*mut*np.outer(u, u))
    v = V[:, -1]
    return float(v @ (C.Of @ v))
lo, hi = 0.30, 0.325
assert q_dense_of(C, 4, lo) > 0 and q_dense_of(C, 4, hi) < 0
for _ in range(60):
    mid = 0.5*(lo+hi)
    if q_dense_of(C, 4, mid) > 0: lo = mid
    else: hi = mid
print(f"   mutilde_crit (first zero of q, dense) = {0.5*(lo+hi):.6f}   "
      f"(note claims 0.32208; probe window start ~0.322 / N=1200 [0.325,0.330])")
# upper end of the bad window
lo2, hi2 = 0.331, 0.36
for _ in range(60):
    mid = 0.5*(lo2+hi2)
    if q_dense_of(C, 4, mid) <= 0: lo2 = mid
    else: hi2 = mid
print(f"   upper end of the non-monotone window   = {0.5*(lo2+hi2):.6f}")

print("\n(e) Lemma-1 hypotheses over a grid: min |c2|, lam_1(A) simplicity, psi_1 in C")
NN = 250
worst_c2 = (np.inf, None); worst_gap = (np.inf, None); cone_bad = []
minpsi = (np.inf, None)
for w in (0.05,0.1,0.3,0.5,0.74,1.0,2.0,5.0,10.0,20.0,50.0):
    A = myH(NN, w); e, V = np.linalg.eigh(A)
    p1 = V[:, -1]*(1.0 if V[0,-1] >= 0 else -1.0)
    if p1.min() < -1e-12 or (p1[:-1]-p1[1:]).min() < -1e-12:
        cone_bad.append(w)
    if e[-1]-e[-2] < worst_gap[0]: worst_gap = (e[-1]-e[-2], w)
    if p1[0] < minpsi[0]: minpsi = (p1[0], w)
    for Nc in (1,2,3,4,5,8,12,20,40,80,120):
        c = V.T @ uc(NN, Nc)
        if abs(c[-2]) < worst_c2[0]: worst_c2 = (abs(c[-2]), (w, Nc))
print(f"   min gap gam over wbar grid: {worst_gap[0]:.4e} at wbar={worst_gap[1]}  (all > 0 => lam_1(A) simple numerically)")
print(f"   min |c_2| = <u_c,psi_2> over (wbar,N_c) grid: {worst_c2[0]:.4e} at {worst_c2[1]}")
print(f"   psi_1 outside the cone C at wbar = {cone_bad}   (empty = always in C)")
print(f"   min psi_1[0] = {minpsi[0]:.4f} at wbar={minpsi[1]}  => c_1 >= psi_1[0] > 0")
