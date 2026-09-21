"""[numerical-evidence] Dense re-check of the single-cell claim from
_chk_grid_coarsen.py: alpha=1.40, beta=0.0 over the whole range [0.30, 20.0].

The greedy search used a 49-point geometric sample; a wide cell could hide a
spike between sample points.  Here: 161 points, plus a local refinement around
the worst one.  Also prints the alpha-sweep at the worst varpi so the choice of
a single alpha can be seen to be non-fragile.
"""
import os, io, contextlib, importlib
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import importlib.util

spec = importlib.util.spec_from_file_location(
    "gc", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "_chk_grid_coarsen.py"))
gc = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(gc)

THETA_STAR = gc.THETA_STAR


def r_at(w, al, beta):
    bs = gc.base(w)
    P = gc.peel_all_betas(w)[list(gc.BETAS).index(beta)]
    return gc.ratio(bs, al, min(P, THETA_STAR, 0.6 * bs["h"])), P, bs["h"]


if __name__ == "__main__":
    AL, BE = 1.40, 0.0
    WS = np.geomspace(0.30, 20.0, 161)
    vals = []
    for w in WS:
        r, P, h = r_at(w, AL, BE)
        vals.append((r, w, P, h))
    vals.sort(reverse=True)
    print(f"dense sweep, alpha={AL}, beta={BE}, 161 points on [0.30,20.0]")
    print(f"{'rank':>4} {'varpi':>8} {'ratio':>9} {'P':>9} {'h':>9} {'P/h':>8}")
    for i, (r, w, P, h) in enumerate(vals[:8]):
        print(f"{i+1:4d} {w:8.4f} {r:9.4f} {P:9.5f} {h:9.5f} {P/h:8.4f}")
    print(f"\nglobal worst over 161 points: {vals[0][0]:.4f} at varpi={vals[0][1]:.4f}")

    w0 = vals[0][1]
    print(f"\nlocal refinement around varpi={w0:.4f} (31 points, +-12%):")
    loc = [(r_at(w, AL, BE)[0], w) for w in np.geomspace(w0 * 0.88, w0 * 1.12, 31)]
    print(f"  refined worst = {max(loc)[0]:.4f} at varpi={max(loc)[1]:.4f}")

    print(f"\nalpha sweep at varpi={w0:.4f} (is the single-alpha choice fragile?):")
    for al in (1.20, 1.30, 1.35, 1.40, 1.45, 1.50, 1.60, 1.80, 2.10):
        print(f"  alpha={al:4.2f}  ratio={r_at(w0, al, BE)[0]:.4f}")
