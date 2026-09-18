"""[numerical-evidence] Which term of min(P_beta, Lam_full(C_R), THETA*, 0.6h)
actually BINDS at the worst point of each cell, for C_R in {1.6, 20, inf}?

Purpose: size the remaining proof obligations.  _chk_CR_saturation.py showed the
12-cell criterion already passes (worst 0.2809 < 1) with NO C_R bound at all, so
C_R is a margin refinement, not load-bearing.  This says what IS load-bearing.

Run: OMP_NUM_THREADS=1 python _chk_CR_attribution.py
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
TS = gc.THETA_STAR

print(f"{'cell':>15} |" + "".join(f"{f'C_R={c}':>22}" for c in ("1.6", "20", "inf")))
for lo, hi, al in pf.CELLS:
    out = []
    for c in (1.6, 20.0, np.inf):
        best, who, ww = -np.inf, "", 0.0
        for w in np.geomspace(lo, hi, 11):
            w = float(w); bs = gc.base(w); h = bs["h"]
            Pb = float(gc.peel_all_betas(w).min())
            Lf = np.inf if np.isinf(c) else pf.lam_full(w, c)
            cand = {"P_b": Pb, "Lam": Lf, "Th*": TS, "0.6h": 0.6 * h}
            k = min(cand, key=cand.get)
            r = gc.ratio(bs, al, cand[k])
            if r > best:
                best, who, ww = r, k, w
        out.append(f"{best:8.4f} {who:>5}@{ww:5.2f}")
    print(f"[{lo:5.2f},{hi:5.2f}] |" + "".join(f"{o:>22}" for o in out))
