"""THE DECIDER CURVE between the AREA ceiling and our mu*-RESOLVED moment ceiling.

CORRECTION OF FRAMING (two separate claims, kept separate per CLAUDE.md):

  (A) WHAT IS WRONG WITH "lam* = (2 C_inf/c(mu*))^2 is THE crossing".
      That number is the crossing of the AREA ceiling against the ONE-LINE
      ceiling  Tc <= C_inf sqrt(lam <w^2>)  ONLY.  The one-line form is the
      lam -> infinity asymptote of our own machine and is loose by exactly
      1/f_E(lam) (f_E = 0.63 at lam=1, 0.93 at lam=5).  Against the
      lam-resolved and mu*-resolved versions the crossing is a DIFFERENT
      number, and at the mu*-resolved level there is NO crossing in lam at
      all at the Einstein point (see (B)).  So a statement of the form
      "lam < lam* => area tighter" is level-dependent and must name its level.

  (B) WHY THE REPLACEMENT IS RIGHT.  Let
          c_mu(lam) := 2 Tc_E(lam,mu*,r)/(lam w_E)   (area efficiency of the
                                                      Einstein spectrum)
          c(mu*,r)  := max_lam c_mu(lam)             (the Leavens constant)
      Our mu*-resolved ceiling is TIGHT on Einstein spectra (it IS Tc_E), so for
      ANY admissible spectrum P with lam, <w>_P, wbar2:
          AREA/OURS = [c(mu*,r)/c_mu(lam)] * (<w>_P/wbar2)          ... (*)
      First factor >= 1 by construction; second <= 1 by Cauchy-Schwarz.  Hence
      the decider is ONE inequality between two measurable ratios:
          AREA is tighter  <==>  <w>_P/wbar2  <  rho_crit(lam,mu*,r)
                                                := c_mu(lam)/c(mu*,r)  <= 1.
      At the Einstein point <w>/wbar2 = 1 > rho_crit, so OURS wins there for
      every lam -- confirming that the finite crossing seen at level (II) is an
      artifact of the degraded one-line form, not a property of the two
      functionals.  The real boundary is a CURVE in the (lam, <w>/wbar2) plane,
      not a point on the lam axis.

  LEVEL of (*): analytic, given the two ceilings and the definition of c.
  LEVEL of the tabulated rho_crit: numerical-evidence (float64).

CONSUMPTION POINT (declared, per CLAUDE.md budget rule): rho_crit(lam,mu*,r) is
consumed by the sentence "for THIS material, use THAT bound" -- i.e. it is the
selector inside  Tc <= min{c(mu*)A, C_inf sqrt(lam <w^2>) f_E Phi_E}.  If a
material's <w>/wbar2 lies above the curve, our factor Phi_E changes the verdict;
if below, it only buys margin and must NOT be booked as a gain.
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
_C = _load("_chk_c_mustar_curve_2026-09-21")
C_INF = _E.C_INF

LAMS = [1.0, 1.25, 1.5, 1.75, 2.0, 2.2929, 2.5, 2.67, 3.0, 3.62, 4.0, 5.0, 6.0, 8.0]
MUS = [0.0, 0.05, 0.10, 0.13, 0.15, 0.20]
R = 10.0


def main():
    t0 = time.time()
    print("=" * 92)
    print("rho_crit(lam,mu*,r) = c_mu(lam)/c(mu*,r):  AREA tighter <=> <w>/wbar2 < rho_crit")
    print(f"r = w_c/w_E = {R:g}   [numerical-evidence]")
    print("=" * 92)
    print()
    CC = {}
    for mus in MUS:
        CC[mus] = _C.maximize_c(mus, R if mus > 0 else 1.0, N=180)
        print(f"  c({mus:4.2f}) = {CC[mus][0]:.6f}   lam_opt = {CC[mus][1]:.4f}")
    print()
    print("Also reported: lam*_(II) = (2 C_inf/c)^2, the crossing against the ONE-LINE")
    print("ceiling, so the two readings can be compared side by side.")
    for mus in MUS:
        print(f"  mu*={mus:4.2f}:  lam*_(II) = {(2*C_INF/CC[mus][0])**2:.4f}")
    print()
    hdr = f"{'lam':>7}" + "".join(f"{('mu*='+format(m,'.2f')):>10}" for m in MUS)
    print(hdr)
    print("-" * len(hdr))
    for lam in LAMS:
        cells = []
        for mus in MUS:
            t, nc = _C.tc_over_wE(lam, mus, R, N=260)
            if t is None:
                cells.append(f"{'n/a':>10}")
                continue
            c_mu = 2.0 * t / lam
            cells.append(f"{c_mu/CC[mus][0]:10.4f}")
        print(f"{lam:7.4f}" + "".join(cells))
    print()
    print("READING: every entry is <= 1 (by construction, c is the max of c_mu), and")
    print("it is the LARGEST <w>/wbar2 for which the AREA ceiling still wins.  A")
    print("material is on OUR side iff its measured <w>/wbar2 EXCEEDS the entry.")
    print("The entry DECREASES away from lam_opt, so strong coupling helps us -- the")
    print("same conclusion as the level-(II) crossing, but now with the spectral")
    print("shape ratio made explicit instead of hidden in 'sqrt(lam)'.")
    print()
    print("Reference: an Einstein spectrum has <w>/wbar2 = 1 exactly; a wide spectrum")
    print("has <w>/wbar2 < 1.  Two analytic reference shapes:")
    print("  * two deltas, equal weight, ratio rho:  <w>/wbar2 = (1+rho)/sqrt(2(1+rho^2))")
    for rho in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0):
        print(f"      rho={rho:5.1f}: {(1+rho)/np.sqrt(2*(1+rho*rho)):.4f}")
    print("  * uniform P on [0,W]: <w>/wbar2 = (W/2)/(W/sqrt(3)) = sqrt(3)/2 = 0.8660")
    print("  So ANY spectrum broader than a 2:1 two-delta (0.9487) or a box (0.8660)")
    print("  sits near the boundary -- the verdict is genuinely material-specific.")
    print()
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
