"""4a refinement: locate the concavity boundary varpi_c(r) (largest varpi at which
H'' > 0 for SOME mu~ in [1e-3,1]) and check N-convergence there.
Also records that the failure set in mu~ need not be an interval.
numerical-evidence (float64)."""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concavity_region_scan import hpp

GRID = np.geomspace(1e-3, 1.0, 121)


def nfail(N, vp, r):
    x = vp*vp; Nc = max(1, int(r*vp))
    s = np.array([hpp(N, x, Nc, m)[0] > 0 for m in GRID])
    worst = max(hpp(N, x, Nc, m)[0] for m in GRID)
    return int(s.sum()), worst, Nc, s


print("boundary varpi_c(r): largest varpi with ANY H''>0 on mu~ in [1e-3,1] (121 pts)")
for r in (3, 5, 10, 20):
    print(f"\n  r = {r}")
    print(f"  {'varpi':>6} {'Nc':>4} {'N':>5} {'#fail/121':>10} {'max H\"':>11} "
          f"{'mu~ argmax':>10}   contiguous?")
    for vp in (0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.6):
        Ncq = max(1, int(r*vp)); N = max(400, 8*Ncq)
        k, w, Nc, s = nfail(N, vp, r)
        if k:
            i = np.argmax([hpp(N, vp*vp, Nc, m)[0] for m in GRID])
            am = f"{GRID[i]:.4f}"
            idx = np.flatnonzero(s)
            cont = "yes" if idx[-1]-idx[0]+1 == len(idx) else "NO (not an interval)"
        else:
            am, cont = "-", "-"
        print(f"  {vp:>6.1f} {Nc:>4} {N:>5} {k:>10} {w:>11.3e} {am:>10}   {cont}")

print("\nN-convergence at the tightest boundary rows")
print(f"{'r':>3} {'varpi':>6} {'Nc':>4} {'N':>6} {'#fail/121':>10} {'max H\"':>12}")
for r, vp in ((10, 0.7), (10, 0.8), (10, 0.9), (10, 1.0), (3, 2.0), (3, 6.0)):
    for N in (400, 800, 1600):
        k, w, Nc, _ = nfail(N, vp, r)
        print(f"{r:>3} {vp:>6.1f} {Nc:>4} {N:>6} {k:>10} {w:>12.4e}")
    print()
