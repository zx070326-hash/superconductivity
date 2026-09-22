#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_fese_signless_amplitude_test_2026-09-22.py

Task (team-lead, 2026-09-22, third pass on this line of work):
  test condition (ii) lambda^ph_ij >= |lambda^sf_ij| for i != j on FeSe,
  using ONLY 1411.2121-printed numbers, in the worst case where SF
  off-diagonal sign is unknown/adversarial (sign-free amplitude bound,
  not Eq.(37)'s elementwise-positivity claim).

DATA (L1, read from a fresh 300dpi render of p.3, this session):
  Phonon matrix (Fig.3a, printed ON the curves):
      lam_++=0.028, lam_--=0.029, lam_+-=0.042, lam_-+=0.048
  SF matrix: paper prints ONLY the matrix's largest eigenvalue lambda_max,
  at five (alpha/alpha_c) points (Fig.3c inset labels, this session's read
  added 0.16 and 0.37 which prior notes had not recorded) plus one z_Se
  point (Fig.4a):
      alpha/alpha_c = 0.90 -> lambda_max = 0.07
      alpha/alpha_c = 0.92 -> lambda_max = 0.16
      alpha/alpha_c = 0.94 -> lambda_max = 0.37
      alpha/alpha_c = 0.98 -> lambda_max = 0.48
      z_Se = 0.22           -> lambda_max = 1.30

BLOCKER (reported, not resolved): the paper NEVER prints the individual
off-diagonal SF entries lambda^sf_{+-}, lambda^sf_{-+} -- only lambda_max
(the Perron/largest eigenvalue of the |.|-integrated, entrywise-nonnegative
SF matrix).  For a general nonnegative asymmetric 2x2 matrix [[a,b],[c,d]],
Perron-Frobenius gives

      lambda_max = (a+d)/2 + sqrt(((a-d)/2)^2 + b*c)   >= sqrt(b*c)

i.e. lambda_max bounds only the GEOMETRIC MEAN of the two off-diagonal
entries, not either one individually, UNLESS the matrix is symmetric
(b=c) or a=d=0.  The paper's own text ("in this limit SF are dominated
by off-diagonal components") is the closest thing to a justification for
treating lambda_max as a proxy for the off-diagonal entries in the
alpha->alpha_c limit, but this is NOT a printed number and NOT exact.

Given this, we report TWO things:
  (A) the RIGOROUS bound: lambda_max >= sqrt(lam^sf_+- * lam^sf_-+), so
      IF lambda^ph_+- * lambda^ph_-+ >= lambda_max^2 the individual
      product test would still be inconclusive (does not pin each factor).
      What IS rigorous: if lambda_max were itself an upper bound on each
      off-diagonal entry (true whenever the matrix is PSD, since for PSD
      M, |M_ij| <= sqrt(M_ii M_jj) <= lambda_max), the comparison below
      is valid.  The |.|-integrated printed matrix need NOT be PSD (the
      entrywise |.| specifically DESTROYS the PSD structure the raw
      Im Lambda^SF(omega) enjoys -- see _chk_essenberger_lambda_ij_psd.py
      DIRECTION AUDIT).  So "lambda_max as an upper bound on each entry"
      is a MODELING ASSUMPTION, not derived here -- flagged explicitly.
  (B) under that explicit, flagged assumption, run condition (ii) at
      every available data point, under BOTH 1.5x-convention branches.
"""

import math

LAM_PH = {("+", "+"): 0.028, ("-", "-"): 0.029,
          ("+", "-"): 0.042, ("-", "+"): 0.048}

# lambda_max(SF) at each reported point (alpha/alpha_c or z_Se), L1 this session
SF_LMAX = [
    ("alpha/alpha_c=0.90", 0.07),
    ("alpha/alpha_c=0.92", 0.16),
    ("alpha/alpha_c=0.94", 0.37),
    ("alpha/alpha_c=0.98", 0.48),
    ("z_Se=0.22          ", 1.30),
]


def main():
    print("=" * 78)
    print("Condition (ii): lambda^ph_ij >= |lambda^sf_ij|  for i != j (worst-case sign)")
    print("Proxy used for |lambda^sf_ij| (i!=j): lambda_max(SF matrix), FLAGGED")
    print("as an explicit modeling assumption (see docstring), not a printed number.")
    print("=" * 78)

    lam_ph_pm = LAM_PH[("+", "-")]
    lam_ph_mp = LAM_PH[("-", "+")]
    print("\nPrinted phonon off-diagonal: lam_ph_+- = %.3f, lam_ph_-+ = %.3f"
          % (lam_ph_pm, lam_ph_mp))

    for convention, factor in [("no 3/2 correction needed (as-printed)", 1.0),
                                ("divide by 1.5 (if 3/2 not yet removed)", 1.0 / 1.5),
                                ("multiply by 1.5 (opposite convention)", 1.5)]:
        print("\n--- convention: %s ---" % convention)
        all_fail = True
        for label, lmax_raw in SF_LMAX:
            lmax = lmax_raw * factor
            fail_pm = lam_ph_pm < lmax
            fail_mp = lam_ph_mp < lmax
            verdict = "FAIL" if (fail_pm or fail_mp) else "pass (at this point)"
            print("  %s: lambda_max(SF)=%.4f  vs  lam_ph_+-=%.3f, lam_ph_-+=%.3f  -> %s"
                  % (label, lmax, lam_ph_pm, lam_ph_mp, verdict))
            if not (fail_pm or fail_mp):
                all_fail = False
        print("  ALL FIVE POINTS FAIL under this convention: %s" % all_fail)

    print("\n" + "=" * 78)
    print("Rigorous PSD-based per-entry bound (Cauchy-Schwarz on a PSD matrix):")
    print("  if SF matrix M is PSD, |M_ij| <= sqrt(M_ii M_jj) <= lambda_max(M).")
    print("  The PRINTED lambda_ij matrix (this paper's convention) is built by")
    print("  entrywise |.| BEFORE omega-integration -- this is explicitly NOT the")
    print("  same object as the PSD Im Lambda^SF(omega); PSD-ness of the printed,")
    print("  |.|-inflated matrix is NOT established here (flagged, not resolved).")
    print("  Absent that, the only rigorous Perron-Frobenius fact is")
    print("      lambda_max >= sqrt(lam^sf_+- * lam^sf_-+)")
    print("  which bounds the GEOMETRIC MEAN, not either entry individually.")
    lhs = lam_ph_pm * lam_ph_mp
    for label, lmax in SF_LMAX:
        print("  %s: lam_ph_+- * lam_ph_-+ = %.6f  vs  lambda_max^2 = %.6f  -> %s"
              % (label, lhs, lmax ** 2,
                 "product test also fails" if lhs < lmax ** 2 else "product test holds"))


    print("\n" + "=" * 78)
    print("SELF-CORRECTION (found while writing up the above): the per-point")
    print("FAIL verdicts printed above rest on an UNPROVEN premise --")
    print("'lambda_max(SF) bounds each individual off-diagonal entry'.")
    print("That premise is FALSE in general for a nonnegative, non-symmetric")
    print("2x2 matrix. Counterexample: M=[[0,100],[0,0]] is entrywise")
    print("nonnegative, its eigenvalues are both 0 (char. poly = lambda^2),")
    print("so lambda_max = 0 while the (1,2) entry is 100 -- lambda_max gives")
    print("ZERO information about an individual off-diagonal entry when the")
    print("matrix is asymmetric and the paired entry is small.")
    print("The printed FeSe phonon matrix IS mildly asymmetric (0.042 vs")
    print("0.048, ratio 1.14), so nothing rules out the SF matrix being much")
    print("more asymmetric -- the paper gives no information to exclude it.")
    print("")
    print("The 'product test' (lambda_max^2 >= lam^sf_+- * lam^sf_-+, true by")
    print("Perron-Frobenius for ANY nonnegative matrix, no symmetry needed)")
    print("is rigorous, but it is the WRONG DIRECTION for refuting condition")
    print("(ii): (ii) holding requires lam_ph_+- >= lam^sf_+- AND")
    print("lam_ph_-+ >= lam^sf_-+, which by AM-GM implies the PRODUCT")
    print("lam_ph_+- * lam_ph_-+ >= lam^sf_+- * lam^sf_-+ -- a NECESSARY")
    print("condition, not sufficient. Failing the product test (product too")
    print("small) WOULD refute (ii) validly -- but only if we knew")
    print("lam^sf_+- * lam^sf_-+ itself, which we don't (only its upper bound")
    print("lambda_max^2, which goes the WRONG way for this purpose: knowing")
    print("bc <= lambda_max^2 does not tell us bc is actually large).")
    print("")
    print("VERDICT (corrected): with ONLY lambda_max(SF) printed and NO other")
    print("structural information (symmetry, diagonal dominance) established")
    print("for the SF off-diagonal entries, condition (ii) is NOT DECIDABLE")
    print("from 1411.2121's printed numbers -- not 'sound FAIL', but genuinely")
    print("UNDETERMINED. The per-point 'FAIL' verdicts printed above in the")
    print("main loop are RETRACTED as unjustified; kept in this script only")
    print("as a documented false start (see W1 discipline).")


if __name__ == "__main__":
    main()
