"""rev_cu_2: (A) is mutilde_* evaluated at the RIGHT point?  The theorem claims
window-wide monotonicity from a SINGLE wbar.  Compute inf over each anchor's
actual N_c-window [Nc/r, (Nc+1)/r)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import RCell

N = 300
anchors = [(3.0,0.13,0.60,10),(2.0,0.13,0.98,10),(2.0,0.10,0.98,10),(1.0,0.13,2.20,10),
           (1.0,0.10,2.00,10),(0.7,0.13,4.56,10),(0.5,0.13,13.5,10),(0.5,0.13,15.6,3),
           (2.0,0.13,3.00,3),(1.0,0.13,2.50,3)]
print(f"{'lam':>5} {'mu*':>5} {'r':>3} {'w_c':>7} {'Nc':>5} {'mut':>7} "
      f"{'mut_*(w_c)':>11} {'window':>17} {'inf mut_* on win':>17} {'verdict':>9}")
for (lam, mus, wc, r) in anchors:
    mut = mus/lam; Nc = max(1, int(np.floor(r*wc)))
    if Nc > N//2: continue
    thr_c = RCell(wc, Nc, N).thr
    wl, wr = Nc/r, (Nc+1)/r
    ws = np.linspace(wl, wr - 1e-9, 9)
    thrs = [RCell(w, Nc, N).thr for w in ws]
    imin = int(np.argmin(thrs))
    ok = "OK" if mut < min(thrs) else "**FAILS**"
    print(f"{lam:>5} {mus:>5} {r:>3d} {wc:>7.2f} {Nc:>5d} {mut:>7.4f} {thr_c:>11.5f} "
          f"[{wl:>6.3f},{wr:>6.3f}) {min(thrs):>10.5f}@{ws[imin]:>5.3f} {ok:>9}")
