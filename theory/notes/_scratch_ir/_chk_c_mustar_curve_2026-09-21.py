"""Leavens(1975)/Carbotte(1990) AREA ceiling constant c(mu*) -- independent recompute,
and the crossing with the MOMENT ceiling  Tc <= C_inf sqrt(lam <w^2>).

WHAT IS PRE-EXISTING (verified by reading the repo's own full text CARB.txt, not
quoted from anyone's summary):
  CARB.txt ~l.1990-2010  (Carbotte, RMP 62, 1027 (1990), Sec. II.H):
    "It follows that the inequality  T_c <= A c(mu*)   (2.50)  holds.  The equality
     applies to a delta function and the inequality to any other spectral shape.
     This result was first established by Leavens (1975) using somewhat different
     arguments."
    "for mu*=0.1, c(mu*)=0.175 and omega_E(mu*)=1.3 A ... the optimum value of mass
     renormalization is quite finite, namely lam=1.54"
    Fig.6 caption: "Constant c(mu*) in the relation k_B T_c = c(mu*) A for the
     maximum T_c associated with a given A as a function of mu* ... Placed on the
     same figure (solid dots) are the results for T_c/A obtained in the case of many
     strong-coupling superconductors ... The solid points all fall below the maximum
     curve, as they must."
  Sec. II.G: "Equations (2.45) and (2.46) depend only on the single material
     parameter B_E [= omega_E/A] and on mu*, provided as did Leavens (1975), we
     NEGLECT the lam dependence that appears in the cutoff on the Coulomb repulsion
     term."   <== the cutoff convention is declared neglected, i.e. r := omega_c/omega_E
     is held FIXED.  So c(mu*) is really c(mu*, r) and the published 0.175 carries an
     undeclared r.  This script measures which r reproduces it.

  A := integral of alpha^2 F(w) dw  (the AREA).  For the Einstein spectrum
  alpha^2F = (lam w_E/2) delta(w-w_E):  A = lam w_E / 2.
  So  c(mu*, r) = max_lam  2 (T_c(lam,mu*,r)/w_E) / lam,  and the optimum obeys
  w_E/A = 2/lam_opt.

WHY THE TWO CEILINGS ARE INDEPENDENT (level: analytic):
  Let P be the normalized spectral measure (dP ~ 2 alpha^2F(w)/(lam w) dw), so
  <w>_P = 2A/lam  and  <w^2>_P = wbar2^2.  Then
      area ceiling  / moment ceiling  =  c(mu*) A / (C_inf sqrt(lam) wbar2)
                                      =  [c(mu*)/(2 C_inf)] sqrt(lam) * (<w>/wbar2).
  Cauchy-Schwarz gives <w> <= wbar2 with equality iff P is a point mass.  Hence
  spectral SPREAD always moves the comparison in favour of the AREA ceiling, and the
  Einstein crossing lambda* = (2 C_inf/c(mu*))^2 is the SMALLEST crossing over the
  whole class.  Neither ceiling implies the other: our Jensen route proves
  Tc(P) <= Tc_E(lam, w_rms) at FIXED lam, and A_E = lam w_rms/2 >= A_P, so the chain
  produces c(mu*) A_E, never c(mu*) A_P.

LEVEL: numerical-evidence (float64, N truncation; convergence in N reported).
"""
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util as _ilu  # noqa: E402


def _load(stem):
    spec = _ilu.spec_from_file_location(
        stem.replace("-", "_"), os.path.join(HERE, stem + ".py"))
    mod = _ilu.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


_E = _load("_phiE_engine_2026-09-21")
C_INF, TWO_PI, MEV_PER_K = _E.C_INF, _E.TWO_PI, _E.MEV_PER_K
wstar_sup = _E.wstar_sup

MUS = [0.0, 0.05, 0.10, 0.13, 0.15, 0.20]
RS = [3.0, 5.0, 6.0, 8.0, 10.0, 20.0]


def _w0(lam, N):
    w, _, _ = wstar_sup(lam, 0.0, 1.0, N=N, sub=10, nbis=50)
    return w


_W0 = {}


def w0_cached(lam, N=180):
    key = (round(lam, 6), N)
    if key not in _W0:
        _W0[key] = _w0(lam, N)
    return _W0[key]


def tc_over_wE(lam, mus, r, N=180):
    """Tc/omega_E via the sup rule, using the PROVEN free bound as the window cut."""
    wf = w0_cached(lam, N)
    if wf is None:
        return None, None
    if mus == 0.0:
        return 1.0 / (TWO_PI * wf), 0
    w, nc, _ = wstar_sup(lam, mus, r, N=N, sub=10, nbis=50, w_free=wf)
    if w is None:
        return None, None
    return 1.0 / (TWO_PI * w), nc


def c_of(lam, mus, r, N=180):
    """2 Tc/(lam w_E) = Tc/A for the Einstein spectrum."""
    t, nc = tc_over_wE(lam, mus, r, N)
    return (None, None) if t is None else (2.0 * t / lam, nc)


def maximize_c(mus, r, N=180, lo=0.55, hi=4.0, coarse=0.025, fine=0.0025):
    lams = np.arange(lo, hi + 1e-12, coarse)
    vals = []
    for lam in lams:
        v, nc = c_of(float(lam), mus, r, N)
        vals.append(-1.0 if v is None else v)
    i = int(np.argmax(vals))
    a = max(lo, lams[i] - coarse)
    b = min(hi, lams[i] + coarse)
    best = (vals[i], float(lams[i]))
    for lam in np.arange(a, b + 1e-12, fine):
        v, nc = c_of(float(lam), mus, r, N)
        if v is not None and v > best[0]:
            best = (v, float(lam))
    return best


def main():
    t0 = time.time()
    print("=" * 78)
    print("c(mu*, r) = max_lam 2 Tc_E(lam,mu*,r)/(lam w_E)   [numerical-evidence]")
    print("Carbotte RMP 1990 / Leavens 1975 quote c(0.1)=0.175, w_E/A=1.3, lam=1.54")
    print("=" * 78)

    # --- N convergence of the mu*=0 endpoint (the one number with no cutoff amb.) ---
    print("N-convergence of c(0) (no Coulomb => r irrelevant):")
    for N in (140, 180, 260, 360):
        cv, lo = maximize_c(0.0, 1.0, N=N)
        print(f"  N={N:4d}  c(0) = {cv:.7f}   lam_opt = {lo:.4f}   "
              f"w_E/A = {2.0/lo:.5f}   lam*=(2C_inf/c)^2 = {(2*C_INF/cv)**2:.5f}")
    print()

    N = 180
    print(f"{'mu*':>5} | " + " | ".join(f"r={r:g}" for r in RS))
    print("  c(mu*,r)   [row 1], lam_opt [row 2], w_E/A=2/lam_opt [row 3], "
          "Nc at opt [row 4]")
    print("-" * 78)
    CTAB = {}
    for mus in MUS:
        rows = [[], [], [], []]
        for r in RS:
            cv, lo = maximize_c(mus, r, N=N)
            _, nc = c_of(lo, mus, r, N)
            CTAB[(mus, r)] = (cv, lo, nc)
            rows[0].append(f"{cv:.5f}")
            rows[1].append(f"{lo:6.3f}")
            rows[2].append(f"{2.0/lo:6.3f}")
            rows[3].append(f"{(-1 if nc is None else nc):6d}")
        print(f"{mus:5.2f} c  " + "  ".join(rows[0]))
        print(f"      lam" + "  ".join(rows[1]))
        print(f"      wEA" + "  ".join(rows[2]))
        print(f"      Nc " + "  ".join(rows[3]))
    print()

    print("COMPARISON WITH THE PUBLISHED NUMBERS (mu*=0.10):")
    for r in RS:
        cv, lo, nc = CTAB[(0.10, r)]
        print("  r=%5.1f:  c=%.5f (pub 0.175, ratio %.3f)   w_E/A=%.4f (pub 1.3)"
              "   lam_opt=%.3f (pub 1.54)" % (r, cv, cv / 0.175, 2.0 / lo, lo))
    print()
    print("CROSSING WITH THE MOMENT CEILING (Einstein case; lam* = (2C_inf/c)^2):")
    print("  lam < lam*: AREA ceiling is tighter.   lam > lam*: MOMENT ceiling wins.")
    print(f"{'mu*':>5} | " + " | ".join(f"r={r:g}" for r in RS))
    for mus in MUS:
        cells = []
        for r in RS:
            cv, lo, nc = CTAB[(mus, r)]
            cells.append(f"{(2*C_INF/cv)**2:7.4f}")
        print(f"{mus:5.2f}   " + "  ".join(cells))
    print()
    print("  ==> DOES mu* MOVE lam*?  c(mu*) DECREASES with mu* => lam* INCREASES:")
    for r in (3.0, 10.0):
        c0 = CTAB[(0.0, r)][0]
        c1 = CTAB[(0.13, r)][0]
        print(f"    r={r:g}: lam*(0)={(2*C_INF/c0)**2:.4f} -> "
              f"lam*(0.13)={(2*C_INF/c1)**2:.4f}  "
              f"(factor {((2*C_INF/c1)**2)/((2*C_INF/c0)**2):.3f})")
    print("    DIRECTION: mu* pushes the crossing to LARGER lam, i.e. it ENLARGES the")
    print("    region where the area ceiling is the tighter of the two.  Hydrides with")
    print("    lam = 1.7-3.6 therefore sit MORE firmly on the area side once mu*>0.")
    print()

    # ------------------------------------------------------- hydride margin table
    print("=" * 78)
    print("HYDRIDE MARGINS against min{AREA, MOMENT} ceilings   [numerical-evidence]")
    print("=" * 78)
    print("NOTE ON RIGOR ASYMMETRY (do not average these two):")
    print("  * moment ceiling: extremality is the repo's Theorem E (Jensen on")
    print("    x/(x+j^2), analytic) + C_inf certificate.")
    print("  * area ceiling: Carbotte's own wording is 'We can prove that (2.48)")
    print("    represents a LOCAL maximum', established via a NUMERICALLY evaluated")
    print("    functional derivative (his Fig.5 lower frame).  Leavens (1975) not")
    print("    read here, so no claim is made about ITS level.  As published, the")
    print("    area ceiling's extremality step is not above numerical-evidence.")
    print("  => min{} inherits the LOWER level.")
    print()
    # lam, wbar2 (meV), w_log (meV), mu*, Tc_expt, label
    # A is NOT in hydrides_theory.csv (neither A nor <w> is tabulated), so the AREA
    # ceiling CANNOT be evaluated for these rows from repo data.  Instead we report
    # the EXACT threshold on the one missing ratio:
    #     area ceiling tighter  <=>  (c/(2 C_inf)) sqrt(lam) * (<w>/wbar2) < 1
    #                           <=>  <w>/wbar2  <  2 C_inf/(c sqrt(lam)) =: R_crit.
    # Since Cauchy-Schwarz forces <w>/wbar2 <= 1, R_crit >= 1 means the area ceiling
    # is tighter for EVERY spectrum with that (lam, mu*) -- no data needed.
    MATS = [
        ("H3S 200 harm",  2.2929, 129.32, 107.165, 0.13, 203.0),
        ("H3S 200 anh",   1.7210, 149.26, 132.60, 0.13, 203.0),
        ("LaH10 129",     3.6200, None,   76.400, 0.10, 250.0),
        ("LaH10 163",     2.6700, None,   96.400, 0.10, 250.0),
        ("LaH10 214",     2.0600, None,  115.500, 0.10, 250.0),
        ("LaH10 264",     1.7300, None,  126.600, 0.10, 250.0),
    ]
    print("  (w_log is NOT <w> and NOT wbar2; rows with no wbar2 in")
    print("   theory/data/hydrides_theory.csv are left n/a, never filled with w_log.)")
    print()
    hdr = (f"{'material':<15}{'lam':>6}{'mu*':>6}{'Tc_ex':>7}"
           f"{'1-line':>9}{'lam-res':>9}{'mu*-res':>9}{'margin':>8}{'R_crit':>8}")
    print(hdr)
    print("-" * len(hdr))
    for name, lam, wb2, wlog, mus, tcx in MATS:
        cv = CTAB[(0.13 if mus == 0.13 else 0.10, 10.0)][0]
        rcrit = 2.0 * C_INF / (cv * np.sqrt(lam))
        if wb2 is None:
            print(f"{name:<15}{lam:6.3f}{mus:6.2f}{tcx:7.1f}"
                  f"{'n/a':>9}{'n/a':>9}{'n/a':>9}{'n/a':>8}{rcrit:8.3f}")
            continue
        t0_, _ = tc_over_wE(lam, 0.0, 1.0, N=260)
        fE = t0_ / (C_INF * np.sqrt(lam))
        c1 = C_INF * np.sqrt(lam) * wb2 / MEV_PER_K
        c2 = c1 * fE
        tmu, nc = tc_over_wE(lam, mus, 10.0, N=260)
        c3 = tmu * wb2 / MEV_PER_K
        print(f"{name:<15}{lam:6.3f}{mus:6.2f}{tcx:7.1f}"
              f"{c1:9.1f}{c2:9.1f}{c3:9.1f}{c3/tcx:8.3f}x{rcrit:8.3f}")
    print()
    print("  R_crit >= 1 means: for THAT (lam,mu*) the area ceiling beats even the")
    print("  ONE-LINE moment ceiling for every admissible spectrum (Cauchy-Schwarz),")
    print("  so no spectral data are needed to decide.  R_crit < 1 means the verdict")
    print("  depends on the spectrum's <w>/wbar2, which the repo CSV does not carry.")
    print()
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
