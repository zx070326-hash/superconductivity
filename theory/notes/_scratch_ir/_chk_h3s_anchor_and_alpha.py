#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_h3s_anchor_and_alpha.py -- 2026-09-18, main session, written from zero.

WHY
---
h3s-pressure-audit (subagent) correctly transcribed Errea PRL 114,157004
Table II and found that the pressure-matched experimental anchor Errea
themselves used at 200 GPa is 190 K, not the 203 K the repo's deficit-budget
note compares against.  Its RECOMMENDATION -- "switch the anchor to 190 K,
the residual flips -2.0% -> +4.7%" -- does not survive, because of what the
D-side row reveals about where Errea's experimental column came from.

L1 EVIDENCE (main session read the page images / text layer itself)
-------------------------------------------------------------------
Errea PRL 114,157004 (2015) TABLE II, verbatim, p.157004-3:
    Compound        lam_har  wlog_har  lam_anh  wlog_anh  Tc_har  Tc_anh  Tc(expt)
    H3S (200 GPa)    2.64      90.4     1.84     92.86     250    194.0     190
    H3S (250 GPa)    1.96     109.1     1.71    101.3      226    190       --
    D3S (200 GPa)    2.64      68.5     1.87     73.3      183    152.0      90
  (mu* = 0.16 throughout; Tc in K, frequencies in meV)

Drozdov et al., Nature 525, 73 (2015), Fig. 2a caption, verbatim:
    "(Note here that the only point for sulfur deuteride presented in
     arXiv:1412.0460 was determined WITHOUT SAMPLE ANNEALING, and Tc would
     increase after annealing at room temperature.)"
  and Fig. 2a body: at 145 GPa, un-annealed Tc = 170 K -> after annealing 195 K.

Drozdov p.2, verbatim:
    "Comparison Tc`s in the pressure range P>170 GPa (Fig. 2c) gives alpha_H
     approx 0.3."

Drozdov p.3, verbatim (ref 23 = Errea):
    "Inclusion of zero point vibrations and anharmonicity in the calculations23
     corrected the calculated Tc to 190 K, and the isotope coefficient from
     alpha=0.5 to alpha=0.35 - both in agreement with the present work."

=> Errea's Tc(expt) column is sourced from the PREPRINT arXiv:1412.0460, whose
   deuteride point the experimentalists later explicitly disavowed as
   un-annealed, and whose hydride points were likewise pre-annealing.  The bias
   DIRECTION is known and one-sided: annealing RAISES Tc.  So 190 K is not a
   clean pressure-matched anchor; it is a pressure-matched LOW-BIASED anchor.

DIRECTION AUDIT (CLAUDE.md one-sided rule)
------------------------------------------
Nothing below is a bound.  Every number is an exact evaluation of a closed
form at stated inputs.  The one-sided statement used is qualitative and its
direction is stated explicitly: annealing raises Tc, so 190 K underestimates
the annealed 200 GPa value => using it makes our prediction look HIGH.
Level: numerical-evidence (float64 arithmetic on L1-transcribed inputs).
"""

import math

M_H = 1.00782503          # CODATA, u
M_D = 2.01410178
L = math.log(M_D / M_H)   # 0.6923787...


def alpha(tc_h, tc_d):
    """Isotope exponent alpha := -dlnTc/dlnM, from an (H, D) Tc pair."""
    return math.log(tc_h / tc_d) / L


def pct(pred, ref):
    """Signed residual of a prediction against a reference, in percent."""
    return 100.0 * (pred / ref - 1.0)


def main():
    print("=" * 72)
    print("H3S anchor audit + isotope exponent, from L1-transcribed inputs")
    print("=" * 72)
    print("ln(M_D/M_H) = %.7f" % L)

    # ---------------- 1. the Tc anchor question -----------------------
    PRED_H = 199.0        # repo deficit budget, mu*=0.13 (repo's own output,
    PRED_D = 155.0        # NOT recomputed here -- see caveat in the report)
    print("\n[1] Tc anchor for the 200 GPa H3S prediction (%.1f K)" % PRED_H)
    for name, ref, note in [
        ("Drozdov 203 K @155 GPa (published, SQUID-confirmed)", 203.0,
         "well-provenanced but PRESSURE-MISMATCHED (dome peak)"),
        ("Errea Tc(expt) 190 K @200 GPa (from arXiv:1412.0460)", 190.0,
         "pressure-matched but UN-ANNEALED => biased LOW"),
    ]:
        print("    %-52s residual %+6.2f%%" % (name, pct(PRED_H, ref)))
        print("    %-52s   ^ %s" % ("", note))
    print("    => the two candidate anchors straddle the prediction.")
    print("       Sign of the residual is NOT DETERMINED by current evidence.")

    # ---------------- 2. the isotope exponent -------------------------
    print("\n[2] Isotope exponent alpha -- far more anchor-robust")
    rows = [
        ("Errea harmonic   (Tc 250 / 183)", 250.0, 183.0, "theory"),
        ("Errea anharmonic (Tc 194.0/152.0)", 194.0, 152.0, "theory"),
        ("repo deficit budget (199 / 155)", PRED_H, PRED_D, "theory"),
        ("Errea Tc(expt) column (190 / 90)", 190.0, 90.0, "un-annealed, disavowed"),
    ]
    print("    %-36s %9s   %s" % ("source", "alpha", "status"))
    print("    " + "-" * 64)
    for nm, th, td, st in rows:
        print("    %-36s %9.4f   %s" % (nm, alpha(th, td), st))
    print("    %-36s %9.4f   %s" % ("Drozdov, annealed, P>170 GPa", 0.30,
                                    "L1 verbatim: 'gives alpha_H approx 0.3'"))
    print("    %-36s %9.4f   %s" % ("Errea's own quoted value", 0.35,
                                    "L1 verbatim, matches its table to 0.6%"))

    a_errea = alpha(194.0, 152.0)
    a_repo = alpha(PRED_H, PRED_D)
    print("\n    cross-check: Errea's table -> %.4f vs their quoted 0.35  (%.1f%%)"
          % (a_errea, abs(pct(a_errea, 0.35))))
    print("    repo vs Errea (two INDEPENDENT anharmonic treatments):")
    print("        %.4f vs %.4f   -> differ by %.1f%%" % (a_repo, a_errea,
                                                          abs(pct(a_repo, a_errea))))
    print("    both vs published experiment 0.30:")
    print("        repo  overshoots by %+.1f%%" % pct(a_repo, 0.30))
    print("        Errea overshoots by %+.1f%%" % pct(a_errea, 0.30))

    # ---------------- 3. why alpha is the robust observable -----------
    print("\n[3] Why alpha survives the anchor problem but Tc does not")
    print("    Tc anchor moves 203 -> 190 K: %.1f%% swing in the residual"
          % abs(pct(PRED_H, 190.0) - pct(PRED_H, 203.0)))
    print("    alpha is a RATIO: a common multiplicative shift of both")
    print("    isotopes' Tc cancels exactly.  Check -- scale both by 0.95:")
    print("        alpha(199,155) = %.6f ; alpha(0.95*199, 0.95*155) = %.6f"
          % (a_repo, alpha(0.95 * PRED_H, 0.95 * PRED_D)))
    print("    identical to machine precision => alpha is blind to any error")
    print("    that multiplies both isotopes equally (incl. a uniform pressure")
    print("    offset along the dome, to first order).")

    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    print("""
  [L1] Errea Table II Tc(expt) = 190 K @200 GPa, 90 K for D3S.  The subagent's
       transcription is CORRECT -- independently re-read from the page image.

  [L1] But that column is sourced from arXiv:1412.0460, whose deuteride point
       Drozdov et al. explicitly disavow in print as un-annealed, and annealing
       RAISES Tc (170 -> 195 K at 145 GPa, their Fig. 2a).  So 190 K is
       pressure-matched but biased LOW by a known-signed amount.

  [numerical-evidence] Therefore the recommendation "switch to 190 K, residual
       flips to +4.7%" does NOT follow.  203 K (mismatched pressure) and 190 K
       (mismatched sample prep) bracket the prediction from opposite sides.
       The residual's SIGN is undetermined; settling it needs the annealed
       Tc(P) curve of Drozdov Fig. 2c read at 200 GPa.

  [numerical-evidence -- the substantive find] On the ISOTOPE exponent the
       anchor problem largely evaporates, because alpha is a ratio:
         repo   0.3610   |  Errea anharmonic 0.3523  |  experiment ~0.30
       Two independent anharmonic treatments agree to 2.5% and BOTH overshoot
       the published experimental alpha by ~17-20%, in the SAME direction.
       This is the H3S instance of the pattern already recorded for PdH
       (our -0.399 and Errea SSCHA -0.38 both overshoot experiment).
""")


if __name__ == "__main__":
    main()
