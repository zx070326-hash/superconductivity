"""[numerical-evidence] Close debt (b): the FULL analytic lambda_2 bound, both pads.

From _chk_padR_rowsum.py:  row sum = D * midpoint sum of Phi_s, so

    sum_m w_m Rk(s,s_m) <= Theta(s) + pad_R(s),   pad_R(s) = (D^2/8) V2[Phi_s]
    D_nn                >= arctan(s)/s - pad_D(s), pad_D(s) = (D^2/8) V2[G]/(2s)

Both pads ADD (row sum needs an upper bound; D_nn is subtracted so needs a lower
bound).  Hence the fully analytic-form bound

    Lam_full(w) := sup_{n>=0} [ Thetatil(s_n) + pad_R(s_n) + pad_D(s_n) ]

Empirically V2[Phi_s] ~ C_R/s, so write pad_R(s) <= (D^2/8) * C_R / s with
C_R := sup_{s>0} s*V2[Phi_s].  A GRID sup under-estimates C_R (fail-open), so we
scan finely, report the grid sup, and then ROUND C_R UP to a safe constant.
Two independent schemes are compared: (i) per-n exact V2[Phi_{s_n}], (ii) the
closed C_R/s envelope.  Scheme (ii) must dominate (i) everywhere it is used.

Run: OMP_NUM_THREADS=1 python _chk_padR_full.py
"""
import os, io, contextlib, importlib.util
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, alias):
    sp = importlib.util.spec_from_file_location(alias, os.path.join(HERE, name))
    m = importlib.util.module_from_spec(sp)
    with contextlib.redirect_stdout(io.StringIO()):
        sp.loader.exec_module(m)
    return m


pr = _load("_chk_padR_rowsum.py", "pr")
gc = _load("_chk_grid_coarsen.py", "gc")
tt = _load("_chk_theta_tilde.py", "tt")
THETA_STAR = gc.THETA_STAR
V2G = 3.0 * np.sqrt(3.0) / 2.0

CELLS = [(0.30, 0.42, 5.0), (0.42, 0.60, 3.6), (0.60, 0.85, 2.7), (0.85, 1.20, 2.1),
         (1.20, 1.70, 1.8), (1.70, 2.40, 1.6), (2.40, 3.40, 1.5), (3.40, 4.80, 1.42),
         (4.80, 6.80, 1.35), (6.80, 9.60, 1.30), (9.60, 13.60, 1.27), (13.60, 20.0, 1.24)]


def theta(s):
    return (np.log1p(s * s) + 2.0 * s * np.arctan(s)) / (2.0 * (1.0 + s * s))


def thetatil(s):
    return theta(s) - np.arctan(s) / s


def lam_full(w, C_R, NMAX=200000):
    """sup_n [ Thetatil(s_n) + (D^2/8)(C_R + V2G/2)/s_n ],  D=1/w."""
    n = np.arange(NMAX, dtype=float)
    s = (n + 0.5) / w
    pad = (1.0 / (w * w) / 8.0) * (C_R + V2G / 2.0) / s
    return float((thetatil(s) + pad).max())


if __name__ == "__main__":
    print("== 1. scan C_R(s) = s * V2[Phi_s] for its supremum ==")
    SS = np.concatenate([np.geomspace(0.05, 0.9, 12), np.linspace(1.0, 3.0, 11),
                         np.geomspace(3.5, 200.0, 14)])
    rec = []
    for s in SS:
        v, k = pr.V2_Phi(float(s), ngrid=4001)
        rec.append((float(s), v, float(s) * v, k))
    print(f"{'s':>10} {'V2[Phi_s]':>12} {'s*V2':>10} {'#crit':>6}")
    for s, v, sv, k in rec:
        print(f"{s:10.4f} {v:12.6f} {sv:10.6f} {k:6d}")
    C_R_grid = max(r[2] for r in rec)
    s_at = [r[0] for r in rec if r[2] == C_R_grid][0]
    C_R = 1.60                      # rounded UP from the grid sup (see below)
    print(f"\n   grid sup  s*V2 = {C_R_grid:.6f}  at s={s_at:.4f}")
    print(f"   grid sup is a LOWER estimate of the true sup -> round UP:"
          f"  C_R := {C_R}")
    print(f"   compare V2[G]/2 = {V2G / 2:.6f}  (the pad_D constant, same form)")

    print("\n== 2. envelope must dominate the exact per-s V2 (fail-open guard) ==")
    bad = [(s, v) for s, v, sv, k in rec if C_R / s < v]
    print(f"   points where C_R/s < V2[Phi_s]: {len(bad)}  "
          f"({'PASS' if not bad else 'FAIL ' + str(bad[:3])})")

    print("\n== 3. Lam_full vs Lam_an (pad_D only) vs discrete Stil ==")
    print(f"{'varpi':>7} {'Lam_full':>10} {'Lam_an':>10} {'Stil(disc)':>11}"
          f" {'full/Stil':>10}")
    V2G16 = V2G / 16.0
    for w in (0.30, 0.60, 1.20, 1.70, 2.40, 3.40, 4.80, 9.60, 20.0, 45.0):
        n = np.arange(200000, dtype=float); s = (n + 0.5) / w
        Lan = float((thetatil(s) + V2G16 / (w * w * s)).max())
        Lf = lam_full(w, C_R)
        St = tt.discrete_sup(w)[0]
        print(f"{w:7.2f} {Lf:10.6f} {Lan:10.6f} {St:11.6f} {Lf / St:10.3f}")

    print("\n== 4. 12-cell with the FULL analytic bound min(P_beta, Lam_full, 0.6h) ==")
    print(f"{'cell':>15} {'alpha':>5} | {'old':>7} {'disc':>8} {'Lam_an':>8}"
          f" {'Lam_full':>9} | {'which':>8}")
    wo = wd = wa = wf = -np.inf
    for lo, hi, al in CELLS:
        ro = rd = ra = rf = -np.inf
        who = ""
        for w in np.geomspace(lo, hi, 11):
            bs = gc.base(w); h = bs["h"]
            Pb = gc.peel_all_betas(w).min()
            St = tt.discrete_sup(w)[0]
            nn = np.arange(200000, dtype=float); ss = (nn + 0.5) / w
            Lan = float((thetatil(ss) + V2G16 / (w * w * ss)).max())
            Lf = lam_full(w, C_R)
            ro = max(ro, gc.ratio(bs, al, 0.6 * h if w <= 2 else THETA_STAR))
            rd = max(rd, gc.ratio(bs, al, min(Pb, St, THETA_STAR, 0.6 * h)))
            ra = max(ra, gc.ratio(bs, al, min(Pb, Lan, THETA_STAR, 0.6 * h)))
            bf = min(Pb, Lf, THETA_STAR, 0.6 * h)
            r = gc.ratio(bs, al, bf)
            if r > rf:
                rf, who = r, ("Lam_full" if bf == Lf else "P_b" if bf == Pb else
                              "Th*" if bf == THETA_STAR else "0.6h")
        wo, wd, wa, wf = max(wo, ro), max(wd, rd), max(wa, ra), max(wf, rf)
        print(f"[{lo:5.2f},{hi:5.2f}] {al:5.2f} | {ro:7.4f} {rd:8.4f} {ra:8.4f}"
              f" {rf:9.4f} | {who:>8}")
    print(f"{'12-cell worst':>15} {'':>5} | {wo:7.4f} {wd:8.4f} {wa:8.4f} {wf:9.4f}")
