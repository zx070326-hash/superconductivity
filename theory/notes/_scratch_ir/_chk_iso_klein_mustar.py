#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_iso_klein_mustar.py  --  2026-09-18, main session, written from zero.

QUESTION
--------
Klein, PRL 39, 574 (1977), p.576 states verbatim:

  "our quantitative prediction for the inverse isotope effect is independent of
   the exact values of mu* that we have used, provided only that mu* is the same
   for both PdH_x and PdD_x."

Our deficit budget says the opposite and loudly: alpha_budget(mu*=0) = -0.076 vs
alpha_budget(mu*=0.085) = -0.399, i.e. ~81% of the negative isotope exponent is
the Coulomb channel, with d(alpha)/d(mu*) ~ -4.3.

Is that a defect of OUR machine, or is Klein's claim simply false on HIS OWN
tooling?  This script settles it on the McMillan/Allen-Dynes exponent, which is
what a 1977 PdH paper used.  No repo imports, no external data beyond the two
numbers I read off Klein's Table I page image myself.

INPUT (level L1 -- main session read the harvest.aps.org page image directly)
    Klein Table I, x = 1.00 row:  lambda(PdH) = 0.54, lambda(PdD) = 0.62
                                  Tc(PdH)     = 7.9 K, Tc(PdD)    = 9.6 K
                                  mu*(x)      = 0.26 N/(1+N) -> 0.085 at x=1.0

*** PREMISE CORRECTION, 2026-09-18, BEFORE THIS SCRIPT WAS EVER CITED ***

    The first version of this file claimed to test Klein "on his own tooling"
    using the McMillan exponent.  THAT PREMISE IS FALSE.  Klein p.576 verbatim:

      "The variables just described, along with the band-structure potentials,
       were used to determine the quantities in Eqs. (1)-(3), which were then
       used to obtain a solution of the ELIASHBERG EQUATIONS for T_c(x)..."

    He solved Eliashberg, not McMillan.  Consequences, stated separately from
    the old error per CLAUDE.md ("old statement wrong" and "new statement right"
    are two independent propositions):

      (a) RETRACTED: "Theta_H/Theta_D = 1.2287 is Klein's phonon input."  It is
          a McMillan-inversion artifact with no counterpart in his calculation.
      (b) RETRACTED: the closure test printed in [1] is a TAUTOLOGY -- x is
          pinned from the same Tc pair it is then checked against, so 5.5e-17
          was guaranteed by construction and validated nothing.
      (c) SURVIVES [analytic]: within the McMillan/AD exponent, alpha is
          independent of mu* IFF lambda_D == lambda_H.  Closed form, exact.
      (d) SURVIVES [numerical-evidence]: with Klein's lambda pair, the McMillan
          exponent gives d(alpha)/d(mu*) = -3.62, same sign and within 16% of
          our Eliashberg machine's -4.3.  Two different Tc functionals agree
          that the mu*-dependence is strong => it is NOT an artifact of our
          Einstein proxy or step-Coulomb cutoff.
      (e) NOT RESOLVED: whether Klein's OWN Eliashberg solution is mu*-flat.
          Settling that needs his alpha^2F, which we do not have.

METHOD
------
McMillan(1968)/Allen-Dynes(1975) share the exponent

    ln Tc = ln(Theta/A) - D(lambda, mu*),   D := 1.04(1+L) / (L - mu*(1+0.62 L))

The prefactor constant A (1.45 vs 1.2) CANCELS in the isotope ratio, so the
result below is identical for McMillan and for Allen-Dynes.  Only the ratio
Theta_H/Theta_D survives, and I pin it from Klein's OWN Tc pair at his OWN
mu*=0.085 -- i.e. the counterfactual "hold the phonon input fixed, move mu*"
is exactly the one his sentence asserts to be flat.

    alpha = [ln Tc_H - ln Tc_D] / ln(M_D/M_H)
          = [ x + D(L_D,mu*) - D(L_H,mu*) ] / ln(M_D/M_H),   x := ln(Theta_H/Theta_D)

so alpha splits EXACTLY into a mu*-free piece and a piece carrying all of the
mu* dependence:

    alpha = alpha_Theta + alpha_lambda(mu*)
    alpha_Theta  = x / ln(M_D/M_H)                       (mu*-independent)
    alpha_lambda = [D(L_D,mu*) - D(L_H,mu*)] / ln(M_D/M_H)

DIRECTION AUDIT (CLAUDE.md one-sided rule): nothing here is a bound.  Every
number is an exact evaluation of a closed-form algebraic expression at given
inputs; the only inequality used is the sanity check on Theta_H/Theta_D, which
is reported, not consumed by any verdict.  Level: numerical-evidence (float64
evaluation of a closed form; the closed form itself is analytic).
"""

import math

# ---------------------------------------------------------------- constants
# CODATA atomic masses (u).  alpha is defined with BARE masses -- see
# isotope_pdh_discriminating_test_2026-09-18.md Sec 5.3.1a for why not reduced.
M_H = 1.00782503
M_D = 2.01410178
LN_MASS = math.log(M_D / M_H)          # 0.692379...

# Klein PRL 39,574 (1977) Table I, x=1.00 row  [L1: read off page image]
LAM_H, LAM_D = 0.54, 0.62
TC_H,  TC_D  = 7.9, 9.6
MU_KLEIN = 0.085


def D(lam, mu):
    """McMillan/Allen-Dynes exponent.  Returns +inf if the denominator dies."""
    den = lam - mu * (1.0 + 0.62 * lam)
    if den <= 0.0:
        return float('inf')
    return 1.04 * (1.0 + lam) / den


def alpha_of(x, mu, lam_h=LAM_H, lam_d=LAM_D):
    """alpha := -dlnTc/dlnM with x := ln(Theta_H/Theta_D) held fixed."""
    return (x + D(lam_d, mu) - D(lam_h, mu)) / LN_MASS


def main():
    print("=" * 74)
    print("Klein 1977 mu*-independence claim, tested on the McMillan/AD exponent")
    print("=" * 74)

    # ---- 0. Klein's own alpha, straight from his Tc pair -----------------
    a_klein = math.log(TC_H / TC_D) / LN_MASS
    print("\n[0] Klein's own numbers")
    print("    lambda_H/lambda_D = %.2f / %.2f   (lambda RISES with mass)" % (LAM_H, LAM_D))
    print("    Tc_H/Tc_D         = %.1f / %.1f K" % (TC_H, TC_D))
    print("    alpha(Klein)      = %+.4f" % a_klein)

    # ---- 1. pin Theta ratio from Klein's own pair at his own mu* ---------
    x = math.log(TC_H / TC_D) + D(LAM_H, MU_KLEIN) - D(LAM_D, MU_KLEIN)
    theta_ratio = math.exp(x)
    print("\n[1] Phonon input pinned from Klein's pair at mu*=%.3f" % MU_KLEIN)
    print("    ln(Theta_H/Theta_D) = %+.6f  ->  Theta_H/Theta_D = %.4f" % (x, theta_ratio))
    print("    sanity: uniform-scaling value sqrt(M_D/M_H) = %.4f" % math.sqrt(M_D / M_H))
    print("            Errea harmonic omega_log ratio       = %.4f" % (205.0 / 150.0))
    print("    -> %.4f is physically sensible (partial substitution, < sqrt ratio): %s"
          % (theta_ratio, "OK" if 1.0 < theta_ratio < math.sqrt(M_D / M_H) else "OUT OF RANGE"))
    # closure check: TAUTOLOGY -- x was pinned from this very pair.  Printed
    # only to show the algebra is self-consistent; it is NOT a validation.
    print("    closure (TAUTOLOGY, x pinned from this pair -- validates nothing):")
    print("      alpha(x, mu*=0.085) = %+.6f  vs Klein %+.6f  (diff %.2e)"
          % (alpha_of(x, MU_KLEIN), a_klein, abs(alpha_of(x, MU_KLEIN) - a_klein)))
    print("    !! Klein solved the ELIASHBERG equations (p.576), not McMillan.")
    print("    !! So Theta_H/Theta_D above is a McMillan-inversion artifact,")
    print("    !! NOT Klein's phonon input.  See header (a),(b).")

    # ---- 2. THE TEST: sweep mu* with the phonon input held fixed ---------
    print("\n[2] Klein's stated counterfactual: same mu* for both isotopes, vary its VALUE")
    print("    (this is precisely the sweep his sentence claims to be flat)\n")
    print("    %-8s %10s %10s %10s %10s" % ("mu*", "Tc_H/Tc_D", "alpha", "alpha_Th", "alpha_lam"))
    print("    " + "-" * 52)
    a_theta = x / LN_MASS
    rows = []
    for mu in (0.0, 0.02, 0.04, 0.06, 0.085, 0.10, 0.12, 0.13, 0.15):
        a = alpha_of(x, mu)
        a_lam = (D(LAM_D, mu) - D(LAM_H, mu)) / LN_MASS
        ratio = math.exp(a * LN_MASS)
        rows.append((mu, a))
        mark = "  <- Klein" if abs(mu - MU_KLEIN) < 1e-12 else ""
        print("    %-8.3f %10.4f %+10.4f %+10.4f %+10.4f%s"
              % (mu, ratio, a, a_theta, a_lam, mark))

    a0 = alpha_of(x, 0.0)
    a85 = alpha_of(x, MU_KLEIN)
    a13 = alpha_of(x, 0.13)
    print("\n    span over mu* in [0, 0.13]:  %+.4f -> %+.4f   (Delta = %.4f)"
          % (a0, a13, abs(a13 - a0)))
    print("    Coulomb share of Klein's own alpha at mu*=0.085:")
    print("        alpha(mu*=0)     = %+.4f" % a0)
    print("        alpha(mu*=0.085) = %+.4f" % a85)
    print("        Coulomb channel  = %+.4f   (%.1f%% of the total)"
          % (a85 - a0, 100.0 * abs(a85 - a0) / abs(a85)))

    # ---- 3. slope, analytic + finite difference --------------------------
    def dD_dmu(lam, mu):
        den = lam - mu * (1.0 + 0.62 * lam)
        return 1.04 * (1.0 + lam) * (1.0 + 0.62 * lam) / den ** 2

    slope_an = (dD_dmu(LAM_D, MU_KLEIN) - dD_dmu(LAM_H, MU_KLEIN)) / LN_MASS
    h = 1e-6
    slope_fd = (alpha_of(x, MU_KLEIN + h) - alpha_of(x, MU_KLEIN - h)) / (2 * h)
    print("\n[3] d(alpha)/d(mu*) at mu*=0.085")
    print("    analytic        = %+.4f" % slope_an)
    print("    finite-difference= %+.4f   (agree to %.2e)"
          % (slope_fd, abs(slope_an - slope_fd)))
    print("    our machine (deficit budget, note Sec 5.4) ~ -4.3")

    # ---- 4. WHEN is Klein's claim exactly true? --------------------------
    print("\n[4] Exact condition for mu*-independence")
    print("    alpha_lambda(mu*) = [D(L_D,mu*) - D(L_H,mu*)] / ln(M_D/M_H)")
    print("    This vanishes identically in mu* IFF L_D == L_H.")
    a_eq = alpha_of(x, 0.0, lam_h=LAM_H, lam_d=LAM_H)
    print("    counterfactual L_D := L_H = %.2f (lambda mass-independent):" % LAM_H)
    for mu in (0.0, 0.085, 0.13):
        print("        mu* = %-6.3f -> alpha = %+.6f" % (mu, alpha_of(x, mu, LAM_H, LAM_H)))
    print("    -> flat to machine precision, and equal to alpha_Theta = %+.4f" % a_eq)
    print("    So Klein's sentence is EXACT only in the limit his own Table I violates:")
    print("    his lambda moves 0.54 -> 0.62, which is the entire source of the effect.")

    # ---- 5. alpha_exp provenance: Klein's Table I vs Errea's quotes -----
    print("\n[5] alpha_exp provenance (read off Klein Table I page image, x=1.00)")
    print("    Klein's OWN experimental columns disagree with the pairs the repo note uses.\n")
    print("    %-40s %11s %9s" % ("source (H / D in K)", "pair", "alpha"))
    print("    " + "-" * 62)
    exp_rows = [
        ("Klein Tab.I ref12 Schirber-Northrup", 8.0, 9.8),
        ("Klein Tab.I ref13 Miller-Satterthwaite", 9.1, 10.3),
        ("Errea ref[8] Schirber-Northrup (quoted)", 8.0, 10.0),
        ("Errea ref[7] Stritzker-Buckel (quoted)", 9.0, 11.0),
    ]
    vals = []
    for nm, th, td in exp_rows:
        a = math.log(th / td) / LN_MASS
        vals.append(a)
        print("    %-40s %5.1f/%-5.1f %+9.4f" % (nm, th, td, a))
    lo, hi = min(vals), max(vals)
    print("\n    SAME EXPERIMENT (Schirber-Northrup), two secondary quotes:")
    print("        Klein 8.0/9.8 -> %+.4f   Errea 8.0/10.0 -> %+.4f   (differ %.4f)"
          % (vals[0], vals[2], abs(vals[0] - vals[2])))
    print("    experimental band actually in evidence: [%+.4f, %+.4f], width %.4f"
          % (lo, hi, hi - lo))
    print("    band assumed by the repo note:          [-0.3223, -0.2898], width 0.0325")
    print("    -> the assumed target band is %.1fx too narrow." % ((hi - lo) / 0.0325))
    print("    Direction audit: adding Miller-Satterthwaite widens the band toward")
    print("    LESS negative, so the NEAREST-edge overshoot of our -0.399 is")
    print("    unchanged at %.4f (set by Errea's SN quote); what collapses is the"
          % abs(-0.3992 - lo))
    print("    claim that the target is known to +-0.016.")

    # ---- 6. verdict ------------------------------------------------------
    print("\n" + "=" * 74)
    print("VERDICT (rewritten after the premise correction)")
    print("=" * 74)
    print("""
  [RETRACTED] "Klein's alpha is not mu*-independent on his own tooling."
      Not established.  His tooling is the Eliashberg equations (p.576); this
      script exercises the McMillan exponent.  Deciding his actual solution
      needs his alpha^2F, which we do not have.

  [analytic] Within the McMillan/AD exponent, alpha is independent of mu*
      IFF lambda is mass-independent.  Klein's own Table I has lambda rising
      0.54 -> 0.62, so within THIS functional form his stated independence
      cannot hold.  Exact, from the closed form -- no computer needed.

  [numerical-evidence] With Klein's lambda pair, McMillan gives
      d(alpha)/d(mu*) = %+.2f, versus our Eliashberg machine's ~ -4.3:
      same sign, same order, within 16%%.  TWO DIFFERENT Tc FUNCTIONALS AGREE
      that alpha is strongly mu*-dependent once lambda(M) rises.
      => our -4.3 is NOT an artifact of the Einstein proxy or step cutoff.

  [reconciliation -- the substantive result]
      Klein attributes the inverse effect to the OPTIC-MODE contribution
      (p.576 + Fig. 2), our budget attributes 81%% of it to the COULOMB
      channel.  These are NOT in conflict; they answer different questions:

          DRIVER    = lambda(M) rising, from optic-mode anharmonicity.
                      Klein: 20%% enhancement of <w^2>_H over harmonic (p.575).
                      Ours : note Sec 5.4, "lambda_D -> 0.40 sends alpha to +0.414".
                      BOTH FRAMEWORKS AGREE.
          AMPLIFIER = the Coulomb term's nonlinearity, which multiplies that
                      driver by ~5x (alpha: -0.076 -> -0.399 in our machine,
                      -0.062 -> -0.282 in McMillan at Klein's mu*).
                      Klein says this is irrelevant; McMillan and our machine
                      both say it is strong and strongly mu*-dependent.

      So the disagreement is narrow and specific: not "what causes the inverse
      isotope effect", but "how much does mu* amplify it".  That is the one
      question Klein's sentence answers differently from every functional form
      we can evaluate.

  [data provenance -- unplanned find, see section 5]
      The experimental target band is [%+.4f, %+.4f] (width %.4f), not the
      [-0.322, -0.290] the repo note assumes.  The same Schirber-Northrup
      experiment is quoted as 8.0/9.8 by Klein and 8.0/10.0 by Errea, worth
      0.029 in alpha.  Miller-Satterthwaite (%+.4f) is absent from the note
      entirely.  Minimum overshoot is unchanged (%.3f); its stated PRECISION
      is not defensible.
""" % (slope_an, lo, hi, hi - lo, vals[1], abs(-0.3992 - lo)))

    neg = all(a < 0 for _, a in rows)
    print("  sign check: alpha < 0 at every mu* in the swept table: %s" % neg)
    print("  (this IS robust, and is the reading of Klein's sentence that survives)")


if __name__ == "__main__":
    main()
