"""[numerical-evidence] Does route R' stay closed under the NEW lambda_2 bound?

routeRpM1 closed R' on the grounds "exact-R buys only single-digit % while
rat/cell x2.25", using gains 3.9-7.7% computed with the OLD analytic gap
    gap_an_old(w,h) = 0.4h (w<=2) else h - Theta*.
The new bound is
    gap_an_new(w,h) = h - min(P_beta(w), Lam_full(w), Theta*, 0.6h)
which moved the 12-cell worst 0.5751 -> 0.2255 and the bottleneck cell from
[3.40,4.80] to [13.60,20.00].  g_an enters `ratio` THREE times, so the RELATIVE
gain of the expansion route is not invariant -- it must be recomputed before the
closure verdict is carried over.

Reuses _chk_expansion_12cell.py's own functions (build/trial_vector/
sinTheta_bound/F_lo_of/best_ub) unchanged; only g_an is swapped.

Run: OMP_NUM_THREADS=1 python _chk_Rprime_under_new_lam2.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, HERE)


def _load(name, alias):
    sp = importlib.util.spec_from_file_location(alias, os.path.join(HERE, name))
    m = importlib.util.module_from_spec(sp)
    with contextlib.redirect_stdout(io.StringIO()):
        sp.loader.exec_module(m)
    return m


ex = _load("_chk_expansion_12cell.py", "ex")      # runs its own loop; output eaten
gc = _load("_chk_grid_coarsen.py", "gc")
pf = _load("_chk_padR_full.py", "pf")             # also runs its own loop
THETA_STAR = gc.THETA_STAR
C_R = 1.60

CELLS = ex.CELLS


def lam2_new(w):
    Pb = gc.peel_all_betas(w).min()
    Lf = pf.lam_full(w, C_R)
    return min(Pb, Lf, THETA_STAR)


if __name__ == "__main__":
    print(f"{'cell':>15} {'alpha':>6} | {'dk_old':>8} {'exp_old':>8} {'g_old':>7}"
          f" | {'dk_new':>8} {'exp_new':>8} {'g_new':>7} | {'gain_old':>8} {'gain_new':>8}")
    wo_dk = wo_ex = wn_dk = wn_ex = 0.0
    rows = []
    for wa, wb, alpha in CELLS:
        best = {}
        for tag in ("old", "new"):
            worst = None
            for w in np.geomspace(wa, wb, 11):
                pr = ex.build(w)
                h = pr["h"]
                g_an = (ex.gap_an(w, h) if tag == "old"
                        else h - min(lam2_new(w), 0.6 * h))
                if g_an <= 0:
                    continue
                t = ex.trial_vector(pr, alpha)
                s = ex.sinTheta_bound(pr, t, g_an)
                x = pr["x"]
                Ophi = ex.O_of(lambda k: k * k / (k * k + x) ** 3, pr["N"])
                nphi = np.abs(np.linalg.eigvalsh(Ophi)).max()
                F_lo = ex.F_lo_of(pr, t, s, Ophi, nphi)
                if F_lo <= 0:
                    continue
                dk_n, ub_n = ex.best_ub(pr, t, s, g_an, m2=True)
                r_dk, r_ex = x * x * dk_n / F_lo, x * x * ub_n / F_lo
                if worst is None or r_ex > worst[1]:
                    worst = (w, r_ex, r_dk, g_an)
            best[tag] = worst
        (_, eo, do, go) = best["old"]
        (_, en, dn, gn) = best["new"]
        wo_dk, wo_ex = max(wo_dk, do), max(wo_ex, eo)
        wn_dk, wn_ex = max(wn_dk, dn), max(wn_ex, en)
        print(f"[{wa:5.2f},{wb:5.2f}] {alpha:6.2f} | {do:8.4f} {eo:8.4f} {go:7.4f}"
              f" | {dn:8.4f} {en:8.4f} {gn:7.4f} |"
              f" {100 * (1 - eo / do):7.1f}% {100 * (1 - en / dn):7.1f}%")
    print(f"\n{'12-cell worst':>15} {'':>6} | {wo_dk:8.4f} {wo_ex:8.4f} {'':>7}"
          f" | {wn_dk:8.4f} {wn_ex:8.4f}")
    print(f"   exact-R gain at the 12-cell worst:  old {100*(1-wo_ex/wo_dk):.1f}%"
          f"   new {100*(1-wn_ex/wn_dk):.1f}%")
