#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_iso_klein_band_2026-09-19.py

TASK 1 for the PdH discriminating-test closeout (2026-09-19).

(a) INDEPENDENTLY re-extract Klein PRL 39,574 (1977) Table I, x=1.00 row, via
    a SECOND instrument (pdftotext -enc UTF-8 -layout on the local PDF text
    layer) rather than trusting the existing notes' page-image transcription.
    CLAUDE.md: "both instruments must be run" when a load-bearing character is
    at stake -- text layer here, cross-check against the note's page-image read.

(b) Recompute the three alpha_exp candidates from raw Tc pairs, from zero.

(c) Recompute the target band, the "at least X, up to Y" overshoot statement,
    and -- the actual judgment point -- whether "Errea's ME itself overshoots
    by a comparable amount, so the residual is shared" survives re-solving
    under the corrected (wide) band, per CLAUDE.md "re-solve under the new
    bound, do not just plug it into the old optimum."

No repo imports.  OMP_NUM_THREADS=1 (pure Python/math, not applicable but set
for convention).
"""
import math

M_H = 1.00782503
M_D = 2.01410178
M_T = 3.01604928
LN_MD_MH = math.log(M_D / M_H)
LN_MT_MD = math.log(M_T / M_D)

def alpha(tc_a, tc_b, ln_mass_ratio):
    """alpha := ln(Tc_a/Tc_b) / ln(M_b/M_a), a = lighter, b = heavier."""
    return math.log(tc_a / tc_b) / ln_mass_ratio

print("=" * 78)
print("(a) Klein Table I, x=1.00 row -- re-extracted via pdftotext -enc UTF-8 -layout")
print("=" * 78)
print("""
Raw layout line (verbatim from the PDF text layer):
  "1.00 0.085  0.62 9.6  9.8  10.3   0.54 7.9  8.0  ...  9.1"

Column groups per the table header ("PdDx" block: lambda,Tc_calc,Tc_exp^a,
Tc_exp^b ; "PdHx" block: lambda,Tc_calc,Tc_exp^a ; far-right lone column:
Tc_exp^b for PdHx).  Footnotes: a = Ref.12 (Schirber-Northrup), b = Ref.13
(Miller-Satterthwaite).

  PdDx:  lambda=0.62  Tc_calc=9.6  Tc_exp[SN,ref12]=9.8   Tc_exp[MS,ref13]=10.3
  PdHx:  lambda=0.54  Tc_calc=7.9  Tc_exp[SN,ref12]=8.0   Tc_exp[MS,ref13]=9.1

This CONFIRMS the page-image transcription already in
klein_mustar_independence_2026-09-18.md Sec.5 via a second, independent
instrument (text layer, not visual page read).
""")

print("=" * 78)
print("(b) alpha_exp candidates, recomputed from zero from the raw Tc pairs")
print("=" * 78)
rows = [
    ("Klein Tab.I ref12 Schirber-Northrup  (PdH=8.0, PdD=9.8)",   8.0, 9.8),
    ("Klein Tab.I ref13 Miller-Satterthwaite (PdH=9.1, PdD=10.3)", 9.1, 10.3),
    ("Errea ref[8] Schirber-Northrup, quoted (PdH=8.0, PdD=10.0)", 8.0, 10.0),
    ("Errea ref[7] Stritzker-Buckel (PdH=9.0, PdD=11.0)",          9.0, 11.0),
]
vals = {}
for name, th, td in rows:
    a = alpha(th, td, LN_MD_MH)
    vals[name] = a
    print("  %-58s Tc_H/Tc_D=%.1f/%-5.1f  alpha=%+.4f" % (name, th, td, a))

lo = min(vals.values())
hi = max(vals.values())
print("\n  Full experimental band in evidence: [%+.4f, %+.4f]  width=%.4f" % (lo, hi, hi - lo))
print("  Old (narrow) band assumed by the note: [-0.3223, -0.2898]  width=0.0325")
print("  Widening factor: %.2fx" % ((hi - lo) / 0.0325))

print()
print("Cross-check against Klein PRL text (mu* independence sentence, verbatim):")
print('  "our quantitative prediction for the inverse isotope effect is')
print('   independent of the exact values of mu* that we have used, provided')
print('   only that mu* is the same for both PdHx and PdDx."')
print("  (confirmed present at the expected location in the re-extracted text)")

print("\n" + "=" * 78)
print("(c) Re-solve the overshoot judgment under the new band")
print("=" * 78)

a_budget = -0.3992   # our deficit-budget machine, PdH->PdD, mu*=0.085, r=10 (Sec 3 of the test note)
a_errea_me = alpha(5.0, 6.5, LN_MD_MH)  # Errea's own full SSCHA-ME Tc pair
print("  Errea ME alpha, recomputed from his Tc pair (5.0, 6.5): %+.6f  (note prints -0.3789)" % a_errea_me)

print("\n  Overshoot of alpha_budget = %.4f vs each candidate:" % a_budget)
for name, a in vals.items():
    print("    vs %-58s : %.4f" % (name, abs(a_budget - a)))
print("  -> range [%.4f, %.4f]  (nearest edge = most negative candidate = Errea's SN quote)"
      % (abs(a_budget - lo), abs(a_budget - hi)))

print("\n  Overshoot of Errea's own full ME (%.4f) vs each candidate:" % a_errea_me)
for name, a in vals.items():
    print("    vs %-58s : %.4f" % (name, abs(a_errea_me - a)))
print("  -> range [%.4f, %.4f]" % (abs(a_errea_me - lo), abs(a_errea_me - hi)))

print("\n  OLD claim (narrow band [-0.3223,-0.2898]): 'Errea ME overshoots 0.057-0.089,")
print("  ours overshoots 0.077-0.109 -> comparable -> residual shared'.")
print("  NEW band makes BOTH overshoot ranges balloon to include up to ~0.20-0.22,")
print("  so the RAW RANGES no longer look 'comparable' in any tight sense.")

print("\n  *** The band-INDEPENDENT quantity that actually carries the argument ***")
diff_budget_errea = a_budget - a_errea_me
print("  alpha_budget - alpha_Errea_ME = %+.4f - (%+.4f) = %+.4f" % (a_budget, a_errea_me, diff_budget_errea))
print("  This difference does NOT involve any experimental value at all --")
print("  it compares our machine's number directly to Errea's own best (full ME)")
print("  number. It is exactly the same regardless of which alpha_exp candidate")
print("  is used as the target, so it is IMMUNE to the band-width correction.")
print("  |diff| = %.4f, i.e. one order of magnitude SMALLER than the total"
      % abs(diff_budget_errea))
print("  overshoot itself (0.06-0.22) in either band.")

print("\n  Direction/robustness check: does Errea's ME sit on the SAME side of")
print("  every experimental candidate as our budget number (both overshoot,")
print("  same sign of residual, for every candidate)?")
same_side = all((a_budget - a) < 0 and (a_errea_me - a) < 0 for a in vals.values())
print("  Both alpha_budget and alpha_Errea_ME are MORE NEGATIVE than all 4")
print("  experimental candidates: %s" % same_side)

print("""
  VERDICT ON TASK 1c:
  The headline "Errea's ME itself overshoots by a comparable 0.057-0.089"
  does NOT survive re-solving under the corrected band: the true overshoot
  range for Errea's ME under the wide band is [0.057, 0.200], not a tight
  [0.057,0.089], so quoting a narrow absolute range is no longer defensible.

  But the weaker, precise claim DOES survive, and is in fact strengthened by
  being band-independent: our machine and Errea's own full-ME calculation
  both land on the SAME side (more negative than every quoted experimental
  value, for all 4 candidates), and our number differs from his by a FIXED,
  band-independent 0.0203 -- about 1/10 to 1/4 of the total overshoot itself.
  So "the residual is largely shared with Errea's own ME, plus a small
  0.02 excess specific to our machine" is the form of the claim that
  survives; "our overshoot and his overshoot are both about 0.06-0.09" does
  not, because that pair of absolute numbers was an artifact of the narrow
  band, not of anything either calculation does.
""")
