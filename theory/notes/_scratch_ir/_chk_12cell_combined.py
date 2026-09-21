"""[numerical-evidence] 12-cell DK criterion with the combined lambda_2 bound

    lam2 <= min( P_beta(w) ,  Stil(w) ,  0.6 h )

where both ingredients are Lemma-Lambda instances:
  P_beta : Lambda applied to H = B - Delta directly, S={0}, q_n=(2n+1)^{-beta/2}
           (good at small varpi, blows up at large varpi)
  Stil   : min-max (repo 1.F.3(c))  lam2(H) <= lam1(R - D),  then Lambda with
           k=0, q=p on (R,D):  sup_n [ sum_m w_m Rk(s_n,s_m) - D_nn ]
           -> continuum constant Thetatil* = 0.0737378 (vs the repo's Theta*=0.63493)
           (flat ~0.074-0.087 in varpi, so it is the large-varpi workhorse)

Baseline for regression: old gap = 0.4h (w<=2) else h - Theta*, worst cell 0.5751.

Run: OMP_NUM_THREADS=1 python _chk_12cell_combined.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("gc", os.path.join(HERE, "_chk_grid_coarsen.py"))
gc = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(gc)
THETA_STAR = gc.THETA_STAR

spec2 = importlib.util.spec_from_file_location("tt", os.path.join(HERE, "_chk_theta_tilde.py"))
tt = importlib.util.module_from_spec(spec2)
with contextlib.redirect_stdout(io.StringIO()):
    spec2.loader.exec_module(tt)

CELLS = [(0.30, 0.42, 5.0), (0.42, 0.60, 3.6), (0.60, 0.85, 2.7), (0.85, 1.20, 2.1),
         (1.20, 1.70, 1.8), (1.70, 2.40, 1.6), (2.40, 3.40, 1.5), (3.40, 4.80, 1.42),
         (4.80, 6.80, 1.35), (6.80, 9.60, 1.30), (9.60, 13.60, 1.27), (13.60, 20.0, 1.24)]
NPT = 11


def bounds_at(w):
    bs = gc.base(w)
    h = bs["h"]
    Pb = gc.peel_all_betas(w)           # array over gc.BETAS
    St = tt.discrete_sup(w)[0]
    return bs, h, Pb, St


if __name__ == "__main__":
    print(f"{'cell':>15} {'alpha':>5} {'b*':>5} | {'old':>7} {'P_beta':>7}"
          f" {'+Stil':>7} | {'which':>6} {'gain':>6}")
    wo = wp = wn = -np.inf
    for lo, hi, al in CELLS:
        ws = np.geomspace(lo, hi, NPT)
        ro = rp = rn = -np.inf
        who = ""
        for w in ws:
            bs, h, Pb, St = bounds_at(w)
            l2_old = 0.6 * h if w <= 2 else THETA_STAR
            ro = max(ro, gc.ratio(bs, al, l2_old))
            bp = min(Pb.min(), THETA_STAR, 0.6 * h)
            rp = max(rp, gc.ratio(bs, al, bp))
            bn = min(Pb.min(), St, THETA_STAR, 0.6 * h)
            r = gc.ratio(bs, al, bn)
            if r > rn:
                rn, who = r, ("Stil" if bn == St else
                              "P_b" if bn == Pb.min() else "0.6h")
        wo, wp, wn = max(wo, ro), max(wp, rp), max(wn, rn)
        print(f"[{lo:5.2f},{hi:5.2f}] {al:5.2f} {'':>5} | {ro:7.4f} {rp:7.4f}"
              f" {rn:7.4f} | {who:>6} {100 * (1 - rn / ro):5.1f}%")
    print(f"{'12-cell worst':>15} {'':>5} {'':>5} | {wo:7.4f} {wp:7.4f} {wn:7.4f}"
          f" | {'':>6} {100 * (1 - wn / wo):5.1f}%")
