"""rev_cu_4: (a) N-convergence at the two TIGHTEST anchors (large N_c);
(b) load-bearing sensitivity of the six float64 constants."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import RCell

def thr_of(gam, q1, b, m, c1, up):
    tmax = np.inf if m == 0 else (-b + np.sqrt(b*b + m*q1))/m
    if np.isinf(tmax):
        return gam/(2*c1*c1)
    return tmax*gam/(2*c1*up + 2*tmax*c1*c1)

print("(a) N-convergence of mutilde_* at the two anchors with the thinnest margin")
for (tag, w, Nc, mut) in (("lam=0.5,r=10,w_c=13.5", 13.5, 135, 0.26),
                          ("lam=0.4,r=3 ,w_c=51.2", 51.2, 153, 0.325),
                          ("lam=0.7,r=10,w_c=4.56", 4.56, 45, 0.13/0.7)):
    print(f"  {tag}   mutilde={mut:.4f}")
    for NN in (300, 400, 600, 800):
        if Nc > NN//2: continue
        t0=time.time(); C = RCell(w, Nc, NN)
        print(f"     N={NN:>4d}: mut_*={C.thr:.6f}  gam={C.gam:.7f} q1={C.q1:.7f} "
              f"b={C.b:.7f} m={C.m:.7f} c1={C.c1:.6f} up={C.up:.6f} "
              f"margin={C.thr/mut:.4f}  ({time.time()-t0:.0f}s)")
    print()

print("(b) elasticities d ln(mut_*)/d ln(X) and the relative error in X that flips each anchor")
anchors = [(3.0,0.13,0.60,10),(2.0,0.13,0.98,10),(1.0,0.13,2.20,10),(1.0,0.10,2.00,10),
           (0.7,0.13,4.56,10),(0.5,0.13,13.5,10),(0.5,0.13,15.6,3),(2.0,0.13,3.00,3),
           (1.0,0.13,2.50,3)]
names = ["gam","q1","b","m","c1","up"]
print(f"  {'anchor':>22} {'margin':>7} " + " ".join(f"{n+'(el/flip%)':>16}" for n in names))
NN = 400
worst = {n: (np.inf, None) for n in names}
for (lam, mus, wc, r) in anchors:
    mut = mus/lam; Nc = max(1, int(np.floor(r*wc)))
    if Nc > NN//2: continue
    C = RCell(wc, Nc, NN)
    base = [C.gam, C.q1, C.b, C.m, C.c1, C.up]
    t0 = thr_of(*base)
    out = []
    for i, n in enumerate(names):
        h = 1e-5
        p = list(base); p[i] = base[i]*(1+h)
        el = (thr_of(*p) - t0)/t0/h
        # relative change in X needed to bring mut_* down to mut  (bisect in log)
        lo, hi = -0.999, 5.0
        target = mut
        def f(s):
            p = list(base); p[i] = base[i]*(1+s); return thr_of(*p) - target
        # direction: decrease threshold
        s = None
        for _ in range(200):
            mid = 0.5*(lo+hi)
            if (f(mid) > 0) == (f(0.0) > 0):
                pass
            lo = lo  # placeholder
            break
        # simple scan
        grid = np.concatenate([-np.geomspace(1e-4, 0.999, 400)[::-1], np.geomspace(1e-4, 5.0, 400)])
        vals = np.array([f(g) for g in grid])
        cand = grid[vals <= 0]
        flip = np.min(np.abs(cand)) if cand.size else np.inf
        out.append(f"{el:>+6.2f}/{100*flip:>7.2f}%")
        if flip < worst[n][0]:
            worst[n] = (flip, f"lam={lam},r={r}")
    print(f"  {f'lam={lam},mu*={mus},r={r}':>22} {t0/mut:>7.4f} " + " ".join(f"{o:>16}" for o in out))
print("\n  smallest relative error in each constant that would flip SOME anchor:")
for n in names:
    print(f"    {n:>4}: {100*worst[n][0]:>8.2f}%   ({worst[n][1]})")
