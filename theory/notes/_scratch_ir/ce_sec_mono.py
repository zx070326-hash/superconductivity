"""tan Theta_P(mu~) strictly increasing => the crossing with t^d is unique,
so thr_SEC is a genuine threshold (not just 'a' root).  Analytic reason:
secular eq gives dlam/dmu~ < 0, and both (mu_1 - lam) and S_2(lam) decrease in lam."""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
sys.argv = ["x"]
import p4d_coulomb_extremal_threshold as ce
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ce_threshold_resolvent_v2 import constants, tan_sec, lam_secular, ANCH

bad = 0; worst_lam = -np.inf
for name, lam, mus, vr in ANCH:
    vr = vr or 0.8709/math.sqrt(lam); oc = 2*math.pi*vr
    for r in (3, 10):
        Nc = max(1, int(r*vr))
        C = ce.CellExtremal(np.array([.5,.5]), np.array([oc/1.05, oc*1.05]), 1., Nc, NN=400)
        K = constants(C)
        ms = np.geomspace(1e-4, 1.0, 400)
        ts = np.array([tan_sec(K, m) for m in ms])
        ls = np.array([lam_secular(K, m) for m in ms])
        dt, dl = np.diff(ts), np.diff(ls)
        if dt.min() <= 0: bad += 1; print(f"  NONMONO tan {name} r={r}: {dt.min():.2e}")
        if dl.max() >= 0: bad += 1; print(f"  NONMONO lam {name} r={r}: {dl.max():.2e}")
        worst_lam = max(worst_lam, dl.max())
print(f"14 rows x 400 mu~ in [1e-4,1]: monotonicity violations = {bad}; "
      f"worst d(lam) = {worst_lam:.2e} (must be < 0)")
