"""[numerical-evidence] What does the 12-cell criterion give if the lambda_2 bound
is restricted to ANALYTIC FORM, i.e. no float64 eigenvalue and no float64 discrete
row sum -- only

    Lam_an(w) := sup_{n>=0} [ Thetatil(s_n) + pad_D(s_n) ] ,  s_n = (n+1/2)/w
    Thetatil(s) = Theta(s) - arctan(s)/s ,  sup_s Thetatil = Thetatil* = 0.0737378
    pad_D(s)    = (D^2/8) * V2[G] / (2 s) = 0.162380 / (w^2 s),  V2[G]=3sqrt(3)/2

pad_D is the midpoint-rule slack on D_nn, which is SUBTRACTED, so it must be
added back (we need a LOWER bound on D_nn -> an UPPER bound on the whole thing).

HONEST OMISSION: the row-sum side still needs its own upward pad,
sum_m w_m Rk(s_n,s_m) <= Theta(s_n) + pad_R(s_n)  (repo M2's q*V2[u] family,
which I have NOT independently re-verified).  pad_R is NOT included below, so
Lam_an here is an under-estimate of the true analytic bound.  An under-estimated
upper bound supports only PASS -> the numbers below are NOT a certificate and a
"pass" here must not be read as one.  They answer only: "how much does insisting
on analytic form cost, relative to the discrete sup?"

Run: OMP_NUM_THREADS=1 python _chk_12cell_analytic_form.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    sp = importlib.util.spec_from_file_location(name[:4], os.path.join(HERE, name))
    m = importlib.util.module_from_spec(sp)
    with contextlib.redirect_stdout(io.StringIO()):
        sp.loader.exec_module(m)
    return m


gc = _load("_chk_grid_coarsen.py")
tt = _load("_chk_theta_tilde.py")
THETA_STAR = gc.THETA_STAR
V2G_OVER_16 = 3.0 * np.sqrt(3.0) / 2.0 / 16.0        # = 0.1623798...

CELLS = [(0.30, 0.42, 5.0), (0.42, 0.60, 3.6), (0.60, 0.85, 2.7), (0.85, 1.20, 2.1),
         (1.20, 1.70, 1.8), (1.70, 2.40, 1.6), (2.40, 3.40, 1.5), (3.40, 4.80, 1.42),
         (4.80, 6.80, 1.35), (6.80, 9.60, 1.30), (9.60, 13.60, 1.27), (13.60, 20.0, 1.24)]


def theta(s):
    return (np.log1p(s * s) + 2.0 * s * np.arctan(s)) / (2.0 * (1.0 + s * s))


def lam_an(w, NMAX=200000):
    n = np.arange(NMAX, dtype=float)
    s = (n + 0.5) / w
    return float((theta(s) - np.arctan(s) / s + V2G_OVER_16 / (w * w * s)).max())


if __name__ == "__main__":
    print("== Lam_an(w) vs discrete sup Stil(w) vs Thetatil* ==")
    print(f"{'varpi':>7} {'Lam_an':>10} {'Stil(disc)':>11} {'Thetatil*':>10}"
          f" {'Lam_an/Stil':>12}")
    for w in (0.30, 0.60, 1.20, 1.70, 2.40, 3.40, 4.80, 9.60, 20.0, 45.0):
        St = tt.discrete_sup(w)[0]
        L = lam_an(w)
        print(f"{w:7.2f} {L:10.6f} {St:11.6f} {0.0737378:10.6f} {L / St:12.3f}")

    print("\n== 12-cell, analytic-form bound  min(P_beta, Lam_an, 0.6h) ==")
    print(f"{'cell':>15} {'alpha':>5} | {'old':>7} {'disc Stil':>10}"
          f" {'Lam_an':>8} | {'which':>6}")
    wo = wd = wa = -np.inf
    for lo, hi, al in CELLS:
        ro = rd = ra = -np.inf
        who = ""
        for w in np.geomspace(lo, hi, 11):
            bs = gc.base(w); h = bs["h"]
            Pb = gc.peel_all_betas(w).min()
            St = tt.discrete_sup(w)[0]
            L = lam_an(w)
            ro = max(ro, gc.ratio(bs, al, 0.6 * h if w <= 2 else THETA_STAR))
            rd = max(rd, gc.ratio(bs, al, min(Pb, St, THETA_STAR, 0.6 * h)))
            ban = min(Pb, L, THETA_STAR, 0.6 * h)
            r = gc.ratio(bs, al, ban)
            if r > ra:
                ra, who = r, ("Lam_an" if ban == L else
                              "P_b" if ban == Pb else
                              "Th*" if ban == THETA_STAR else "0.6h")
        wo, wd, wa = max(wo, ro), max(wd, rd), max(wa, ra)
        print(f"[{lo:5.2f},{hi:5.2f}] {al:5.2f} | {ro:7.4f} {rd:10.4f}"
              f" {ra:8.4f} | {who:>6}")
    print(f"{'12-cell worst':>15} {'':>5} | {wo:7.4f} {wd:10.4f} {wa:8.4f}")
