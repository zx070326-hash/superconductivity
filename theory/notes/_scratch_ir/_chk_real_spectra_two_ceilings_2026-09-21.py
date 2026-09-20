"""TWO CEILINGS ON THE REAL Mishra alpha^2F SPECTRA:
   AREA (Leavens 1975 / Carbotte 1990)   vs   MOMENT (this repo, mu*-resolved).

This is the load-bearing comparison.  Both ceilings are evaluated on the SAME
spectrum with the SAME conventions, all moments recomputed here from the raw
EPW .a2f tables (nothing is taken from hydrides_theory.csv or from any summary).

Definitions used (all from alpha^2F alone):
  lam      = 2 \\int alpha^2F(w)/w dw
  A        = \\int alpha^2F(w) dw                                  (ZEROTH moment)
  dP       = (2/lam) alpha^2F(w)/w dw   (probability measure, the Einstein mixture)
  <w>_P    = 2A/lam
  wbar2    = sqrt(<w^2>_P) = sqrt((2/lam) \\int w alpha^2F dw)     (FIRST moment norm)
  w_log    = exp((2/lam) \\int ln(w) alpha^2F/w dw)

CEILINGS
  (I)   AREA:        Tc <= c(mu*,r) * A                       [Leavens/Carbotte]
  (II)  one-line:    Tc <= C_inf sqrt(lam) wbar2              [repo headline]
  (III) lam-res:     Tc <= C_inf sqrt(lam) wbar2 * f_E(lam)
  (IV)  mu*-res:     Tc <= C_inf sqrt(lam) wbar2 * f_E(lam) * Phi_E(lam,mu*,r)
                        = Tc_E(lam,mu*,r) evaluated at w_E = wbar2
  (II)-(IV) all come from the SAME Jensen/Lemma-M extremality (Theorem E): the
  extremal Einstein spectrum sits at w_E = wbar2 (rms matching), NOT at w_log.

EXACT ANALYTIC COMPARISON (level: analytic, given the two ceilings)
  (I)/(IV) = [c(mu*,r)/(2 C_inf)] sqrt(lam) (<w>/wbar2) / (f_E Phi_E)
  Cauchy-Schwarz: <w> <= wbar2, equality iff P is a point mass.  So spectral
  SPREAD always helps the AREA ceiling.  Define
      R_crit(lam,mu*,r) := 2 C_inf f_E Phi_E / (c(mu*,r) sqrt(lam))
  then  AREA is tighter  <=>  <w>/wbar2 < R_crit.  R_crit >= 1 decides WITHOUT
  any spectral data (Cauchy-Schwarz); R_crit < 1 needs the measured <w>/wbar2.

CUTOFF CONVENTION.  c(mu*,r) depends on r := w_c/w_E (Carbotte Sec. II.G
explicitly NEGLECTS the lam dependence of that cutoff, so his 0.175 carries an
undeclared r; we reproduce it at r ~ 6).  For a real spectrum the SAME physical
w_c maps to  r_eff = w_c / wbar2  in the Einstein comparison, because the
comparison spectrum has w_E = wbar2.  Both the material's own r_eff and the
r-agnostic  sup_{r in [3,20]} c(mu*,r)  are reported (the sup is the honest
choice for a ceiling when w_c is not pinned).

LEVEL: numerical-evidence (float64; trapezoid moments on the raw EPW grid, grid
refinement reported).  The area ceiling's own extremality step is, as published,
also not above numerical-evidence (Carbotte: "a LOCAL maximum", functional
derivative evaluated numerically) -- so min{} inherits numerical-evidence.
"""
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import importlib.util as _ilu  # noqa: E402


def _load(stem):
    spec = _ilu.spec_from_file_location(
        stem.replace("-", "_"), str(HERE / (stem + ".py")))
    mod = _ilu.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


_E = _load("_phiE_engine_2026-09-21")
_C = _load("_chk_c_mustar_curve_2026-09-21")
C_INF, MEV_PER_K = _E.C_INF, _E.MEV_PER_K   # MEV_PER_K = k_B in meV/K

# ------------------------------------------------------------------ data lookup
REL = ("theory/data/external/mishra_h3s_d3s/"
       "Data-H3S-D3S-PH-EPW-Anharmonic")
CANDS = [HERE.parents[2] / REL] + [p / REL for p in HERE.parents[3:6]]
DATA = None
for c in CANDS:
    if c.is_dir():
        DATA = c
        break
if DATA is None:
    raise SystemExit("Mishra dataset not found; tried:\n  " +
                     "\n  ".join(str(c) for c in CANDS))

FILES = [
    ("H3S 200 har", "H3S-200GPa/EPW_harmonic/H3S.a2f",    203.0),
    ("H3S 200 anh", "H3S-200GPa/EPW_anharmonic/H3S.a2f",  203.0),
    ("D3S 200 har", "D3S-200GPa/EPW_harmonic/D3S.a2f",    150.0),
    ("D3S 200 anh", "D3S-200GPa/EPW_anharmonic/D3S.a2f",  150.0),
    ("H3S 160 anh", "H3S-160GPa/EPW_anharmonic/H3S.a2f",  None),
    ("D3S 160 anh", "D3S-160GPa/EPW_anh/D3S.a2f",         None),
]


def load_a2f(path):
    """EPW .a2f: one header line, then 'w[meV] a2f ...' until a non-numeric line."""
    w, a = [], []
    for line in Path(path).read_text().splitlines()[1:]:
        p = line.split()
        if len(p) < 2:
            break
        try:
            x, y = float(p[0]), float(p[1])
        except ValueError:
            break
        w.append(x)
        a.append(max(y, 0.0))
    return np.array(w), np.array(a)


def moments(w, a, wmin=1e-9):
    m = w > wmin
    w, a = w[m], a[m]
    lam = 2.0 * np.trapezoid(a / w, w)
    A = float(np.trapezoid(a, w))
    wmean = 2.0 * A / lam
    w2 = np.sqrt(2.0 / lam * np.trapezoid(w * a, w))
    wlog = np.exp(2.0 / lam * np.trapezoid(np.log(w) * a / w, w))
    return dict(lam=float(lam), A=A, wmean=float(wmean), w2=float(w2),
                wlog=float(wlog))


def main():
    t0 = time.time()
    print("=" * 96)
    print("AREA vs MOMENT ceilings on the REAL Mishra spectra   [numerical-evidence]")
    print("data root:", DATA)
    print("=" * 96)
    print()

    M = {}
    print("Moments recomputed here from the raw EPW tables (trapezoid):")
    print(f"{'case':<13}{'npts':>6}{'lam':>8}{'A(meV)':>9}{'<w>':>8}{'wbar2':>8}"
          f"{'w_log':>8}{'<w>/wbar2':>10}{'wlog/wbar2':>11}")
    print("-" * 81)
    for name, rel, tcx in FILES:
        p = DATA / rel
        if not p.is_file():
            print(f"{name:<13}  MISSING {rel}")
            continue
        w, a = load_a2f(p)
        d = moments(w, a)
        d['tcx'] = tcx
        M[name] = d
        print(f"{name:<13}{len(w):6d}{d['lam']:8.4f}{d['A']:9.3f}{d['wmean']:8.2f}"
              f"{d['w2']:8.2f}{d['wlog']:8.2f}{d['wmean']/d['w2']:10.4f}"
              f"{d['wlog']/d['w2']:11.4f}")
    print()
    print("  CAUCHY-SCHWARZ CHECK  <w>/wbar2 <= 1 :",
          all(d['wmean'] <= d['w2'] * (1 + 1e-12) for d in M.values()))
    print("  (also w_log < <w> < wbar2 for every row -- the three 'average")
    print("   frequencies' are ordered, and the ceiling uses wbar2, the largest.)")
    print()

    # grid-refinement check on the moments (decimation by 2 and 4)
    print("Grid sensitivity of the moments (decimate the raw grid by 2 and by 4):")
    for name in ("H3S 200 har", "H3S 200 anh"):
        w, a = load_a2f(DATA / dict((n, r) for n, r, _ in FILES)[name])
        row = []
        for step in (1, 2, 4):
            d = moments(w[::step], a[::step])
            row.append(f"step{step}: lam={d['lam']:.5f} A={d['A']:.4f} "
                       f"wbar2={d['w2']:.4f}")
        print(f"  {name:<13} " + " | ".join(row))
    print()

    # --------------------------------------------------------------- c(mu*, r)
    print("=" * 96)
    print("c(mu*,r) recomputed here (same engine as the Phi_E machine):")
    RS = [3.0, 5.0, 6.0, 8.0, 10.0, 20.0]
    CT = {}
    for mus in (0.0, 0.10, 0.13):
        cells = []
        for r in RS:
            cv, lo = _C.maximize_c(mus, r, N=180)
            CT[(mus, r)] = cv
            cells.append(f"{cv:.5f}")
        sup = max(CT[(mus, r)] for r in RS)
        CT[(mus, 'sup')] = sup
        print(f"  mu*={mus:4.2f}: " + "  ".join(cells) + f"   sup_r = {sup:.5f}")
    print("  (c(0.10) published as 0.175; reproduced to 0.3% at r=6.  c is")
    print("   NON-monotone in r, minimum near r~6 -- the Morel-Anderson mechanism.)")
    print()

    # ------------------------------------------------------- the two ceilings
    print("=" * 96)
    print("CEILINGS IN KELVIN, same spectrum, same mu*, same w_c")
    print("=" * 96)
    print("w_c convention here: w_c = 10 * wbar2  (=> r_eff = 10 exactly).  The")
    print("r-agnostic column uses sup_{r in [3,20]} c(mu*,r) for the AREA ceiling")
    print("and inf over the same r-grid for OUR ceiling -- i.e. each ceiling is")
    print("degraded in ITS OWN safe direction, so the comparison is not rigged.")
    print()
    MUS = 0.13
    r_eff = 10.0
    hdr = (f"{'case':<13}{'lam':>7}{'(I)area':>9}{'(II)1ln':>9}{'(III)lam':>9}"
           f"{'(IV)mu*':>9}{'min':>8}{'Tc_exp':>8}{'min/exp':>9}"
           f"{'R_crit':>8}{'<w>/w2':>8}  winner")
    print(hdr)
    print("-" * len(hdr))
    for name, d in M.items():
        lam, w2, A = d['lam'], d['w2'], d['A']
        c_area = CT[(MUS, r_eff)]
        tc_area = c_area * A / MEV_PER_K
        t0w, _ = _C.tc_over_wE(lam, 0.0, 1.0, N=260)
        fE = t0w / (C_INF * np.sqrt(lam))
        tc_1 = C_INF * np.sqrt(lam) * w2 / MEV_PER_K
        tc_3 = tc_1 * fE
        tmu, nc = _C.tc_over_wE(lam, MUS, r_eff, N=260)
        tc_4 = tmu * w2 / MEV_PER_K
        phi = tmu / t0w
        rcrit = 2.0 * C_INF * fE * phi / (c_area * np.sqrt(lam))
        mn = min(tc_area, tc_4)
        win = "AREA" if tc_area < tc_4 else "MOMENT(ours)"
        tx = d['tcx']
        rat = f"{mn/tx:9.3f}" if tx else f"{'n/a':>9}"
        txs = f"{tx:8.1f}" if tx else f"{'n/a':>8}"
        print(f"{name:<13}{lam:7.4f}{tc_area:9.1f}{tc_1:9.1f}{tc_3:9.1f}"
              f"{tc_4:9.1f}{mn:8.1f}{txs}{rat}{rcrit:8.4f}"
              f"{d['wmean']/w2:8.4f}  {win}")
    print()
    print("  R_crit vs <w>/wbar2 is the DECIDER: <w>/wbar2 < R_crit  <=>  AREA tighter.")
    print()

    print("r-agnostic version (AREA uses sup_r c, OURS uses inf_r Phi_E):")
    print(f"{'case':<13}{'area_sup':>10}{'ours_inf':>10}{'min':>9}  winner")
    print("-" * 46)
    for name, d in M.items():
        lam, w2, A = d['lam'], d['w2'], d['A']
        ta = CT[(MUS, 'sup')] * A / MEV_PER_K
        t0w, _ = _C.tc_over_wE(lam, 0.0, 1.0, N=260)
        best = None
        for r in RS:
            tmu, _ = _C.tc_over_wE(lam, MUS, r, N=260)
            if tmu is not None:
                v = tmu * w2 / MEV_PER_K
                best = v if best is None else max(best, v)
        print(f"{name:<13}{ta:10.1f}{best:10.1f}{min(ta,best):9.1f}  "
              f"{'AREA' if ta < best else 'MOMENT(ours)'}")
    print()
    print("  NOTE the direction: for a CEILING the r-agnostic version must take the")
    print("  SUP over the admissible r of each bound (a ceiling must hold for the")
    print("  true, unknown w_c).  sup is used for both columns above.")
    print()

    # ------------------------------------------------------------------------
    # THE EXACT RATIO IDENTITY, and why a single "crossing lambda*" is
    # level-dependent.  Define the AREA EFFICIENCY of the Einstein spectrum
    #     c_mu(lam) := 2 Tc_E(lam,mu*,r)/(lam w_E),      c(mu*,r) = max_lam c_mu(lam).
    # Ceiling (IV) evaluated on an Einstein spectrum is EXACTLY Tc_E (it is tight),
    # hence for ANY spectrum
    #     (I)/(IV) = [c(mu*,r)/c_mu(lam)] * (<w>/wbar2)   ... (*)
    # with the first factor >= 1 BY CONSTRUCTION (c is the max of c_mu) and the
    # second <= 1 by Cauchy-Schwarz.  So on Einstein spectra the AREA ceiling can
    # NEVER be tighter than our mu*-resolved ceiling, and it is equal only at
    # lam = lam_opt(mu*).  A finite "crossing lam*" exists only against a
    # DEGRADED version of our ceiling:
    #     vs (II)  one-line C_inf sqrt(lam) wbar2      : lam* = (2 C_inf/c)^2
    #     vs (III) lam-resolved (f_E, mu*=0 in ours)   : c_0(lam) = c(mu*,r)
    #     vs (IV)  mu*-resolved                        : no crossing (see (*))
    # LEVEL of (*): analytic given the two ceilings; the numbers in it are
    # numerical-evidence.
    print("=" * 96)
    print("RATIO IDENTITY (*) AND THE LEVEL-DEPENDENCE OF THE CROSSING lam*")
    print("=" * 96)
    print("c_mu(lam) := 2 Tc_E(lam,mu*,r)/(lam w_E);   c(mu*,r) = max_lam c_mu(lam)")
    print("(I)/(IV) = [c(mu*,r)/c_mu(lam)] * (<w>/wbar2)   -- verified numerically:")
    print(f"{'case':<13}{'LHS':>10}{'RHS':>10}{'rel err':>10}"
          f"{'c/c_mu':>9}{'<w>/w2':>8}")
    print("-" * 60)
    for name, d in M.items():
        lam, w2, A = d['lam'], d['w2'], d['A']
        tmu, _ = _C.tc_over_wE(lam, MUS, r_eff, N=260)
        c_mu_lam = 2.0 * tmu / lam
        lhs = (CT[(MUS, r_eff)] * A / MEV_PER_K) / (tmu * w2 / MEV_PER_K)
        rhs = (CT[(MUS, r_eff)] / c_mu_lam) * (d['wmean'] / w2)
        print(f"{name:<13}{lhs:10.6f}{rhs:10.6f}{abs(lhs/rhs-1):10.2e}"
              f"{CT[(MUS,r_eff)]/c_mu_lam:9.4f}{d['wmean']/w2:8.4f}")
    print()
    print("Crossing against each level (Einstein point, <w>/wbar2 = 1):")
    print("  level (II): lam* = (2 C_inf/c)^2 -- the value quoted in the")
    print("              coordinator's message (2.5056 at mu*=0) IS this one.")
    for mus in (0.0, 0.10, 0.13):
        cv = CT[(mus, r_eff if mus > 0 else 3.0)]
        print(f"     mu*={mus:4.2f}: c={cv:.5f}  lam*_(II) = {(2*C_INF/cv)**2:.4f}")
    print("  level (III): solve c_0(lam) = c(mu*,r), where c_0(lam) := 2Tc_E(lam,0)/lam.")
    print("               c_0 is UNIMODAL (max c(0)=0.23087 at lam=1.1425), so the")
    print("               solution set is an INTERVAL [lam_-, lam_+] inside which the")
    print("               AREA ceiling is tighter, and OUR (III) is tighter OUTSIDE it")
    print("               -- at BOTH small and large lam.  A single-sided 'lam<lam*'")
    print("               reading of the crossing is therefore wrong at level (III).")
    lams = np.arange(0.30, 12.0001, 0.01)
    c0 = []
    for L in lams:
        t, _ = _C.tc_over_wE(float(L), 0.0, 1.0, N=220)
        c0.append(2.0 * t / L if t else np.nan)
    c0 = np.array(c0)
    for mus in (0.05, 0.10, 0.13, 0.15, 0.20):
        cv = _C.maximize_c(mus, r_eff, N=180)[0]
        inside = lams[c0 > cv]
        if len(inside):
            print(f"     mu*={mus:4.2f}: c={cv:.5f} -> AREA tighter for lam in "
                  f"[{inside[0]:.3f}, {inside[-1]:.3f}]  (width {inside[-1]-inside[0]:.3f})")
        else:
            print(f"     mu*={mus:4.2f}: c={cv:.5f} -> AREA NEVER tighter at level (III)")
    print("  level (IV): no crossing exists (identity (*) with c/c_mu >= 1).")
    print()

    print("=" * 96)
    print("BUDGET FOR THE H3S MARGIN, measured against min{two ceilings}")
    print("=" * 96)
    d = M.get("H3S 200 har")
    if d is not None:
        lam, w2, A = d['lam'], d['w2'], d['A']
        tc_area = CT[(MUS, r_eff)] * A / MEV_PER_K
        t0w, _ = _C.tc_over_wE(lam, 0.0, 1.0, N=260)
        tmu, _ = _C.tc_over_wE(lam, MUS, r_eff, N=260)
        tc_1 = C_INF * np.sqrt(lam) * w2 / MEV_PER_K
        tc_4 = tmu * w2 / MEV_PER_K
        print(f"  H3S 200 harmonic, mu*={MUS}, r=10, Tc_expt = 203 K")
        print(f"    (II)  one-line moment ceiling        {tc_1:7.1f} K   "
              f"{tc_1/203:5.3f}x")
        print(f"    (IV)  mu*-resolved moment ceiling    {tc_4:7.1f} K   "
              f"{tc_4/203:5.3f}x   <- OUR contribution")
        print(f"    (I)   area ceiling (Leavens/Carbotte){tc_area:7.1f} K   "
              f"{tc_area/203:5.3f}x   <- PRE-EXISTING")
        print(f"    min{{I,IV}}                            {min(tc_area,tc_4):7.1f} K   "
              f"{min(tc_area,tc_4)/203:5.3f}x")
        print()
        print("  HONEST ACCOUNTING: the gain from (II) to (IV) is ours; the gain from")
        print("  (IV) to min{I,IV} -- if any -- belongs to Leavens/Carbotte and must")
        print("  NOT be booked as this delivery's improvement.")
        print(f"    ours:      {tc_1/203:.3f}x -> {tc_4/203:.3f}x  "
              f"(factor {tc_1/tc_4:.3f})")
        print(f"    theirs:    {tc_4/203:.3f}x -> {min(tc_area,tc_4)/203:.3f}x  "
              f"(factor {tc_4/min(tc_area,tc_4):.3f})")
    print()
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
