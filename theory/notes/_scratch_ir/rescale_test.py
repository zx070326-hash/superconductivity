"""RED FLAG CHECK: is the 'anharmonic' a2F in the Mishra set an independent
calculation, or a uniform frequency rescaling of the harmonic one?
Hypothesis H: a2F_anh(om) = (1/s^2) a2F_har(om/s)  =>  preserves 2 int a2F om dom
exactly, gives lam_anh/lam_har = 1/s^2 and sq<om2>_anh/sq<om2>_har = s."""
import os, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_fill_factor import moments, ROOT

for tag in ("H3S-200GPa", "D3S-200GPa"):
    ga = [x for x in glob.glob(f"{ROOT}/{tag}/EPW_anharmonic/*.a2f") if not x.endswith(".a2fv")]
    gh = [x for x in glob.glob(f"{ROOT}/{tag}/EPW_harmonic/*.a2f") if not x.endswith(".a2fv")]
    la, o2a, _, wa, aa = moments(ga[0])
    lh, o2h, _, wh, ah = moments(gh[0])
    s_from_lam = math.sqrt(lh/la)
    s_from_om  = math.sqrt(o2a/o2h)
    print(f"\n=== {tag} ===  files:\n  anh {os.path.basename(ga[0])}  n={len(wa)} range [{wa[0]:.4f},{wa[-1]:.4f}]"
          f"\n  har {os.path.basename(gh[0])}  n={len(wh)} range [{wh[0]:.4f},{wh[-1]:.4f}]")
    print(f"  s from lam ratio  sqrt(lh/la) = {s_from_lam:.6f}")
    print(f"  s from <om2> ratio            = {s_from_om:.6f}   rel diff {abs(s_from_lam-s_from_om)/s_from_om:.2e}")
    # direct pointwise test of H on the shared grid
    if len(wa) == len(wh) and np.allclose(wa, wh):
        s = s_from_om
        pred = np.interp(wa/s, wh, ah)/s**2
        m = aa > 0.01*aa.max()
        if m.sum():
            r = np.abs(pred[m]-aa[m])/aa[m]
            print(f"  SAME grid. pointwise a2F_anh vs (1/s^2)a2F_har(om/s): "
                  f"median rel {np.median(r):.3e}  max rel {r.max():.3e}  (n={m.sum()})")
    else:
        print(f"  grids differ: ratio of endpoints {wa[-1]/wh[-1]:.6f} vs s={s_from_om:.6f}")
        # test whether the anh grid IS the har grid stretched
        if len(wa) == len(wh):
            rr = wa/wh
            print(f"  om_anh/om_har pointwise: min {rr.min():.6f} max {rr.max():.6f} "
                  f"median {np.median(rr):.6f}")
            ratio_a = aa[ah > 0]/ah[ah > 0] if len(aa)==len(ah) else None
            if ratio_a is not None and len(ratio_a):
                print(f"  a2F_anh/a2F_har  min {ratio_a.min():.6f} max {ratio_a.max():.6f} "
                      f"median {np.median(ratio_a):.6f}   (1/s^2 = {1/s_from_om**2:.6f})")
