"""rev_cu_7: (a) r=10 below wbar=0.34 -- are the jumps really all upward?
(b) the 'six digits from N=200' claim of S6/H1, checked at the anchors S6 omits."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import RCell, uc, myH

NN = 300
def lam1(w, Nc, mt):
    A = myH(NN, w); u = uc(NN, Nc)
    return float(np.linalg.eigvalsh(A - 2*mt*np.outer(u, u))[-1])

print("(a) r=10 stepped path BELOW the note's scan window [0.34,6]")
for mt in (0.065, 0.10, 0.13, 0.26):
    ws = np.arange(0.05, 0.40, 0.002)
    vals = np.array([lam1(w, max(0, int(np.floor(10*w))), mt) for w in ws])
    run = np.maximum.accumulate(vals)
    bad = np.where(vals < run - 1e-12)[0]
    if len(bad) == 0:
        print(f"   mut={mt}: monotone on [0.05,0.40]")
    else:
        print(f"   mut={mt}: {len(bad)} depressed grid points, worst dip "
              f"{(vals-run)[bad].min():.3e}, first at wbar={ws[bad[0]]:.3f} "
              f"(lam1 {run[bad[0]]:.6f} -> {vals[bad].min():.6f}"
              f"  => NON-UNIQUE for lambda in ({1/run[bad[0]]:.2f}, {1/vals[bad].min():.2f}))")

print("\n(b) S6's '6 significant digits from N=200' -- checked at the anchors S6 omits")
print(f"   {'anchor':>16} {'N=200':>11} {'N=400':>11} {'N=800':>11} {'N=1600':>11} {'rel |400-1600|':>15}")
for (w, Nc, tag) in ((1.0, 10, "S6 cell (1,10)"), (0.74, 4, "S6 cell (.74,4)"),
                     (4.56, 45, "lam=0.7 r=10"), (13.5, 135, "lam=0.5 r=10"),
                     (15.6, 46, "lam=0.5 r=3"), (51.2, 153, "lam=0.4 r=3")):
    vs = []
    for N2 in (200, 400, 800, 1600):
        vs.append(RCell(w, Nc, N2).thr if Nc <= N2//2 else np.nan)
    rel = abs(vs[1]-vs[3])/vs[3]
    print(f"   {tag:>16} " + " ".join(f"{v:>11.6f}" for v in vs) + f" {rel:>15.2e}")
