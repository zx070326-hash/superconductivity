"""Post-fix verification of the secular-route threshold (W1 discipline):
(a) c_2 not degenerate (bracket endpoint of the secular root is genuine),
(b) N-convergence of thr_SEC at the two tightest anchors."""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.argv = ["x"]
import p4d_coulomb_extremal_threshold as ce
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from ce_threshold_resolvent_v2 import constants, thr, tan_sec, tan_old, tmax, ANCH

print("(a) c_j magnitudes near the top (N=400)")
print(f"{'anchor':>10} {'r':>2} {'c1':>9} {'c2':>10} {'c3':>10} {'gam':>9}")
for name, lam, mus, vr in ANCH:
    vr = vr or 0.8709/math.sqrt(lam); oc = 2*math.pi*vr
    for r in (3, 10):
        Nc = max(1, int(r*vr))
        C = ce.CellExtremal(np.array([.5,.5]), np.array([oc/1.05, oc*1.05]), 1., Nc, NN=400)
        K = constants(C)
        print(f"{name:>10} {r:>2} {K['c1']:>9.5f} {K['cj'][-1]:>10.3e} "
              f"{K['cj'][-2]:>10.3e} {K['gam']:>9.5f}")

print("\n(b) N-convergence of thr_SEC / thr_OLD at the tight anchors")
print(f"{'anchor':>10} {'r':>2} {'N':>5} {'mu~':>7} {'thr_OLD':>9} {'thr_SEC':>9} "
      f"{'m_SEC':>7} {'tmax':>7} {'tan_sec':>8}")
for name, lam, mus, vr in [a for a in ANCH if a[0] in ("MgB2", "H3S 200", "Nb")]:
    vr = vr or 0.8709/math.sqrt(lam); mut = mus/lam; oc = 2*math.pi*vr
    for r in (3, 10):
        Nc = max(1, int(r*vr))
        for NN in (400, 800, 1600):
            C = ce.CellExtremal(np.array([.5,.5]), np.array([oc/1.05, oc*1.05]), 1., Nc, NN=NN)
            K = constants(C)
            ts = thr(K, tan_sec, cap_kappa=False)
            print(f"{name:>10} {r:>2} {NN:>5} {mut:>7.4f} {thr(K,tan_old):>9.5f} "
                  f"{ts:>9.5f} {ts/mut:>7.3f} {tmax(K):>7.4f} {tan_sec(K,mut):>8.4f}")
