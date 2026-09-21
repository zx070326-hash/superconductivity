"""rev_cu_3: HOSTILE SCAN.  Is there ANY (wbar,N_c,mutilde) with mutilde <
mutilde_*(wbar,N_c) but <vhat,O[f]vhat> <= 0 ?  That would be fatal.
q is evaluated by a DENSE eigh of the pencil (not by the secular engine)."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import RCell, uc

N = 250
t0 = time.time()
ws = [0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.74,0.9,1.0,1.3,1.7,2.0,2.5,3.0,4.0,5.0,7.0,10.0,20.0,50.0]
Ncs = [1,2,3,4,5,6,7,8,10,12,15,20,30,40,60,80,120]
worst = (np.inf, None)
worst_ratio = (np.inf, None)   # min over cells of (min q at mut<=thr)/... just track q
bad = []
lowest_thr = (np.inf, None)
rows = 0
for w in ws:
    for Nc in Ncs:
        C = RCell(w, Nc, N)
        thr = C.thr
        if not np.isfinite(thr) or thr <= 0:
            bad.append((w, Nc, "thr non-positive", thr)); continue
        if thr < lowest_thr[0]:
            lowest_thr = (thr, (w, Nc))
        muts = np.linspace(1e-4, thr*(1-1e-9), 25)
        A = C.A; u = uc(N, Nc); Of = C.Of
        qs = []
        for mt in muts:
            e, V = np.linalg.eigh(A - 2.0*mt*np.outer(u, u))
            v = V[:, -1]
            qs.append(float(v @ (Of @ v)))
        qs = np.array(qs)
        rows += 1
        if qs.min() <= 0:
            bad.append((w, Nc, "q<=0 inside the theorem's hypothesis", float(qs.min())))
        if qs.min() < worst[0]:
            worst = (float(qs.min()), (w, Nc))
        # also: does the Lemma-2 bound actually hold on the whole allowed range?
print(f"scanned {rows} (wbar,N_c) cells x 25 mutilde each, N={N}, {time.time()-t0:.0f}s")
print(f"min over the WHOLE scan of <vhat,O[f]vhat> inside mutilde<mutilde_*: {worst[0]:.4e} at {worst[1]}")
print(f"smallest mutilde_* on the grid: {lowest_thr[0]:.5f} at {lowest_thr[1]}")
print(f"VIOLATIONS: {len(bad)}")
for b in bad[:40]:
    print("   ", b)
