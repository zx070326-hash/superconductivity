"""[numerical-evidence] How large may C_R be before the 12-cell criterion degrades?

The lambda_2 bound actually used is  min(P_beta, Lam_full(C_R), THETA*, 0.6h).
THETA* = 0.63493 (the repo's pre-existing Schur bound, which does NOT depend on
C_R) sits inside that min, so Lam_full can blow up without the criterion failing:
the worst case can never exceed the OLD (C_R-free) column.  This script measures
that explicitly, so the T1 target is sized against the real requirement instead
of against a number I picked by hand.

Run: OMP_NUM_THREADS=1 python _chk_CR_saturation.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
def _load(name, alias):
    sp = importlib.util.spec_from_file_location(alias, os.path.join(HERE, name))
    m = importlib.util.module_from_spec(sp)
    with contextlib.redirect_stdout(io.StringIO()):
        sp.loader.exec_module(m)
    return m

gc = _load("_chk_grid_coarsen.py", "gc")
pf = _load("_chk_padR_full.py", "pf")
THETA_STAR = gc.THETA_STAR
CELLS = pf.CELLS
CRS = [1.6, 20.0, 25.0, 30.0, 50.0, 1e2, 1e3, 1e6, np.inf]

print(f"THETA_STAR = {THETA_STAR:.8f}")
WS = {}
for lo, hi, al in CELLS:
    for w in np.geomspace(lo, hi, 11):
        WS[float(w)] = None
print(f"precomputing {len(WS)} grid points ...", flush=True)
cache = {}
for w in WS:
    bs = gc.base(w)
    cache[w] = (bs, float(gc.peel_all_betas(w).min()))
print("done.\n")

hdr = "".join(f"{('C_R=inf' if np.isinf(c) else f'C_R={c:g}'):>11}" for c in CRS)
print(f"{'cell':>15} {'alpha':>5} |{hdr}   |{'old':>8}")
worst = {c: -np.inf for c in CRS}; worst_old = -np.inf
for lo, hi, al in CELLS:
    row = {c: -np.inf for c in CRS}; ro = -np.inf
    for w in np.geomspace(lo, hi, 11):
        w = float(w); bs, Pb = cache[w]; h = bs["h"]
        ro = max(ro, gc.ratio(bs, al, 0.6 * h if w <= 2 else THETA_STAR))
        for c in CRS:
            Lf = np.inf if np.isinf(c) else pf.lam_full(w, c)
            row[c] = max(row[c], gc.ratio(bs, al, min(Pb, Lf, THETA_STAR, 0.6 * h)))
    for c in CRS:
        worst[c] = max(worst[c], row[c])
    worst_old = max(worst_old, ro)
    print(f"[{lo:5.2f},{hi:5.2f}] {al:5.2f} |"
          + "".join(f"{row[c]:11.4f}" for c in CRS) + f"   |{ro:8.4f}")
print(f"{'12-cell worst':>15} {'':>5} |"
      + "".join(f"{worst[c]:11.4f}" for c in CRS) + f"   |{worst_old:8.4f}")
print(f"\nTARGET (grid_coarsen) = {gc.TARGET};  pass requires ratio < 1.")
