#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_ummarino_arsenide_lambda_ij.py -- 2026-09-18, main session, from zero.

Independent re-verification of source-hunt's claim that the multiband
criterion  lam_ii lam_jj >= lam_ij lam_ji  fails on the published
three-band coupling matrix for IRON ARSENIDES.

SOURCE, read by the main session from the TEXT LAYER of
    theory/notes/_scratch_ir/_ummarino0904.1808.pdf
  (fetched independently from arxiv.org/pdf/0904.1808v2; 4 pp, text PDF)
  = G. A. Ummarino, M. Tortello, D. Daghero, R. S. Gonnelli,
    "Three-band s+- Eliashberg theory and the superconducting gaps of iron
    pnictides", Phys. Rev. B 80, 172503 (2009).

L1 QUOTES (verbatim, text layer)
--------------------------------
p.2, after the parameter-reduction argument:
  "We include these features in the most simple three-band model by posing:
   lam^ph_ii = lam^ph_ij = 0, lam^sp_ii = 0 and mu*_ii(w_c) = mu*_ij(w_c) = 0.
   In addition, we set Gamma^N_ij = Gamma^M_ij = 0 ...
   Under these approximations, the electron-boson coupling-constant matrix is
   then [8]:
        [ 0        0        lam31 nu1 ]
        [ 0        0        lam32 nu2 ]
        [ lam31    lam32    0         ]
   where nu1 = N1(0)/N3(0), nu2 = N2(0)/N3(0) and Ni(0) is the normal density
   of states at the Fermi level for the i-band"

p.2/p.3 (sentence spans the page break -- BOTH halves must be read):
  "In the 122 case (Tc = 37 K) we know that nu1 = 1 and nu2 = 2 [8] while in
   the 1111 case (Tc = 52 K) we have nu1 = 0.4 and nu2 = 0.5 [18]."
  >>> The second half alone ("nu1 = 0.4 and nu2 = 0.5") is the 1111 compound.
      Table I is the 122 compound, so nu1 = 1, nu2 = 2 is the right pair there.
      source-hunt's values are CORRECT; reading only the visible fragment on
      p.3 would have given the wrong ones.

p.3, TABLE I (verbatim):
  "lam31 = 3.866 ; lam32 = 0.471 ; mu*_ij = 0 , lam_ii = 0.4 ;
   D1 = 10.30 meV  D2 = 5.62 meV  |D3| = 10.24 meV
   TABLE I: The effect of a small contribution of intraband coupling,
   lam_ii = 0.4 for the case of Ba0.6K0.4Fe2As2 at T = 2 K."
p.4, TABLE II (verbatim):
  "lam31 = 2.730 ; lam32 = 0.758 ; mu*_ij = 0.1 , lam_ii = 0.4 ;
   D1 = 7.49 meV  D2 = 5.72 meV  |D3| = 7.98 meV
   TABLE II: The effect of the Coulomb interaction, mu*_ij for the case
   shown in Table I."          <- "the case shown in Table I" => same 122 compound
p.3, Fig. 3 parameter sets (verbatim) -- TWO MORE published sets:
  "The parameters used for the 122 compound are w0 = 10 meV, lam31 = 4.267 and
   lam32 = 0.569; for the 1111 compound we used w0 = 10 mev, lam31 = 14.520 and
   lam32 = 1.708."
  >>> CORRECTED 8th round.  An earlier version of this comment claimed these two
      sets "close the open item '1111's lam31/lam32 are only in Fig. 2'" and
      implied they were new relative to source-hunt.  BOTH WRONG:
      (a) its appendix C.3 (6th round) already lists 4.267/0.569 and
          14.520/1.708 with the provenance "p.3 left column body text";
      (b) the "only in Fig. 2" open item was MY OWN earlier comment -- it never
          appears in its note (repo-wide grep: this file and my own note only).
      What IS new here: TABLE II (p.4, which it had not read at the time), the
      assumption-free rescue threshold below, and the lam_ii = 0.4 substitution
      on the Fig. 3 sets.  No number in this script changes.
p.4: "In the case of Ba0.6K0.4Fe2As2 we used lam_ii = 0.4 since we know
   indeed that this coupling cannot be very high [5]."
p.4, the authors' own characterisation of the physics:
  "the newly discovered iron pnictides very likely represent a case of dominant
   negative interband-channel pairing superconductivity ... the coupling
   constants are very high (lam31 > 4)"
  >>> interband dominance is the AUTHORS' claim, not our reading.

THE ASSUMPTION-FREE QUANTITY: THE RESCUE THRESHOLD
--------------------------------------------------
For the pair (i,3), with the printed off-diagonals lam_i3 = lam3i nu_i and
lam_3i = lam3i, the criterion lam_ii lam_33 >= lam_i3 lam_3i reads

        sqrt( lam_ii lam_33 )  >=  lam3i * sqrt(nu_i)                    (*)

The right-hand side contains ONLY published numbers -- no assumption about the
diagonal at all.  So instead of asking "does it fail at lam_ii = 0.4", ask
"how large would the intraband coupling have to be to rescue it".  That number
is reported for every published parameter set below, and compared against the
authors' own 0.4 ("cannot be very high").

DIRECTION AUDIT.  lhs of (*) is INCREASING in the diagonal, rhs is fixed.
=> the verdict at any assumed diagonal d is conservative for every d' <= d.
Reporting the threshold makes the whole one-sided question moot: the criterion
holds iff the true geometric-mean diagonal exceeds it.

CONVENTION [analytic] -- OPPOSITE to Sanna PRB 85,184514
--------------------------------------------------------
Printed: lam_13 / lam_31 = nu1 = N1/N3.
  * lam_ij = N_i V_ij (DOS on the LEFT index) gives lam_13/lam_31 = N1/N3  OK
  * lam_ij = N_j V_ij (DOS on the RIGHT index, = Sanna Eq.(11)) would give
    N3/N1, i.e. the reciprocal                                            NO
Detailed balance confirms: N_3 lam_13 = N3 lam31 N1/N3 = N1 lam31 = N_1 lam_31.
=> Two papers in the same field use OPPOSITE conventions.  This is what makes
   the convention-invariance theorem (lam_ij = f(i)g(j)V_ij => the ratio
   lam_ij lam_ji / (lam_ii lam_jj) is invariant) LOAD-BEARING for cross-paper
   comparison, not a decoration.  Both readings are evaluated below.

WHAT THIS IS AND IS NOT
-----------------------
This is a PHENOMENOLOGICAL FIT with only TWO free parameters (lam31, lam32)
tuned to reproduce Tc and the small gap; nu1, nu2 are inputs taken from [8]/
[18].  It is NOT a first-principles band-resolved lambda_ij table.  The zero
diagonal is a MODELLING CHOICE ("by posing"), adopted because the authors drop
the phonon channel entirely -- note they actually assert lam^ph_ii >> lam^ph_ij,
i.e. phonons are INTRAband dominated.  So the interband dominance here is a
property of the retained SF channel, not of the full coupling.
That is exactly why the load-bearing quantity is the rescue threshold above,
not the lam_ii = 0 row.
"""

import itertools
import math

# ---- every published (lam31, lam32) parameter set in the paper ----------
#      (tag, compound, lam31, lam32, nu1, nu2, provenance)
SETS = [
    ("T1", "122  Ba0.6K0.4Fe2As2", 3.866, 0.471, 1.0, 2.0,
     "TABLE I  (mu*_ij = 0, lam_ii = 0.4, T = 2 K)"),
    ("T2", "122  Ba0.6K0.4Fe2As2", 2.730, 0.758, 1.0, 2.0,
     "TABLE II (mu*_ij = 0.1, lam_ii = 0.4) -- authors' DISFAVOURED scenario"),
    ("F3a", "122  Ba0.6K0.4Fe2As2", 4.267, 0.569, 1.0, 2.0,
     "Fig. 3 text, w0 = 10 meV (primary model, lam_ii = 0)"),
    ("F3b", "1111 SmFeAsO0.8F0.2", 14.520, 1.708, 0.4, 0.5,
     "Fig. 3 text, w0 = 10 meV (primary model, lam_ii = 0)"),
]
LAM_II = 0.4                    # the authors' own plausible intraband value


def matrix(diag, lam31, lam32, nu1, nu2):
    return [[diag, 0.0, lam31 * nu1],
            [0.0, diag, lam32 * nu2],
            [lam31, lam32, diag]]


def check(M, tag, labels=("1", "2", "3")):
    print("\n  %s" % tag)
    print("    [[%.4f, %.4f, %.4f], [%.4f, %.4f, %.4f], [%.4f, %.4f, %.4f]]"
          % tuple(M[i][j] for i in range(3) for j in range(3)))
    nviol = 0
    for i, j in itertools.combinations(range(3), 2):
        lhs = M[i][i] * M[j][j]
        rhs = M[i][j] * M[j][i]
        if rhs == 0.0 and lhs == 0.0:
            verdict, r = "degenerate (0 vs 0)", float("nan")
        elif lhs == 0.0:
            verdict, r = "FAIL (lhs = 0)", float("inf")
        elif rhs == 0.0:
            verdict, r = "pass (rhs = 0 exactly: decoupled pair)", 0.0
        else:
            r = rhs / lhs
            verdict = "FAIL %.3f x" % r if r > 1 else "pass (margin %.3f x)" % (1 / r)
        if r > 1:
            nviol += 1
        print("      (%s,%s)  lam_ii lam_jj = %-10.5f  lam_ij lam_ji = %-10.5f  %s"
              % (labels[i], labels[j], lhs, rhs, verdict))
    return nviol


def main():
    print("=" * 74)
    print("Ummarino et al., PRB 80, 172503 (2009) -- iron ARSENIDES")
    print("criterion  lam_ii lam_jj >= lam_ij lam_ji")
    print("=" * 74)

    print("\n--- (A) the paper's PRIMARY model: lam_ii = 0 'by posing' ---")
    check(matrix(0.0, 3.866, 0.471, 1.0, 2.0), "T1 couplings, diag = 0")
    print("    -> trivially fails wherever an off-diagonal pair is nonzero;")
    print("       this is the published realisation of the one-line analytic")
    print("       counterexample g = [[0,g12],[g12,0]].  But the zero diagonal")
    print("       is a MODELLING CHOICE, so this row alone proves little.")

    print("\n--- (B) ASSUMPTION-FREE: the rescue threshold, all published sets ---")
    print("    criterion for pair (i,3)  <=>  sqrt(lam_ii lam_33) >= lam3i sqrt(nu_i)")
    print("    the rhs uses ONLY published numbers.\n")
    print("    %-4s %-22s %-10s %-10s %s"
          % ("set", "compound", "thr (1,3)", "thr (2,3)", "vs authors' lam_ii = 0.4"))
    print("    " + "-" * 82)
    for tag, comp, l31, l32, nu1, nu2, prov in SETS:
        t13 = l31 * math.sqrt(nu1)
        t23 = l32 * math.sqrt(nu2)
        print("    %-4s %-22s %-10.4f %-10.4f  (1,3) needs %5.1f x , (2,3) needs %4.1f x"
              % (tag, comp, t13, t23, t13 / LAM_II, t23 / LAM_II))
    print("\n    -> EVERY published set needs an intraband coupling ABOVE the")
    print("       authors' own 0.4, on BOTH pairs (smallest threshold 0.666).")
    print("       For (1,3) the threshold runs 2.73-9.18, i.e. 0.63-1.00 x the")
    print("       paper's own INTERband lam31 -- rescuing the criterion needs an")
    print("       intraband coupling comparable to the interband one, which")
    print("       contradicts the paper's premise ('cannot be very high [5]')")
    print("       and its conclusion ('dominant interband-channel pairing').")
    for tag, comp, l31, l32, nu1, nu2, prov in SETS:
        print("         %-4s thr(1,3)/lam31 = %.4f / %.3f = %.2f"
              % (tag, l31 * math.sqrt(nu1), l31, math.sqrt(nu1)))

    print("\n--- (C) the same, expressed as violation factors at lam_ii = 0.4 ---")
    for tag, comp, l31, l32, nu1, nu2, prov in SETS:
        n = check(matrix(LAM_II, l31, l32, nu1, nu2),
                  "%s  %s  [%s]" % (tag, comp, prov))
        print("      -> %d / 3 pairs violate" % n)
    print("\n    CAVEAT: lam_ii = 0.4 is published only for T1/T2.  For F3a/F3b")
    print("    the published diagonal is 0, so imposing 0.4 there is MY choice,")
    print("    made to be conservative -- it is the largest value the authors")
    print("    endorse anywhere.  The (B) thresholds carry no such choice.")

    print("\n--- (D) convention invariance: transpose / opposite DOS placement ---")
    M = matrix(LAM_II, 3.866, 0.471, 1.0, 2.0)
    MT = [[M[j][i] for j in range(3)] for i in range(3)]
    check(MT, "T1 matrix TRANSPOSED (= the opposite DOS convention)")
    print("    -> identical verdicts, as the invariance theorem requires.")
    print("       This matters here: Ummarino puts the DOS on the LEFT index,")
    print("       Sanna PRB 85,184514 Eq.(11) puts it on the RIGHT.  Without")
    print("       the invariance theorem the two papers could not be compared.")

    print("\n--- (E) how many independent material-side data points? ---")
    print("    T1 / T2 / F3a are the SAME compound (122) fitted under three")
    print("    different assumptions => ONE data point, three robustness rows.")
    print("    F3b is a DIFFERENT compound (1111) => a second data point.")
    print("    Both are phenomenological 2-parameter fits, not first principles.")

    nu_sensitivity()


def nu_sensitivity():
    """How much does the 1111 verdict depend on an unsourceable input?

    source-hunt (8th round) found, and I re-read from the p.4 reference list of
    the same PDF, verbatim:

        "[18] I.I. Mazin, private communication."

    Ref. [18] is the ONLY source for the 1111 DOS ratios nu1 = 0.4, nu2 = 0.5.
    So the F3b row rests on an input nobody can trace.  The right response is
    not to drop the data point but to ask how far nu would have to move before
    the verdict changes.  The rescue threshold makes that a one-liner:

        threshold(i,3) = lam3i * sqrt(nu_i)  >=  lam_ii   <=>  nu_i >= (lam_ii/lam3i)^2

    DIRECTION AUDIT: the threshold is INCREASING in nu_i, so the smallest nu_i
    that still fails is the binding one; any larger nu_i fails harder.

    LIMIT, stated honestly: lam31/lam32 for 1111 were FITTED GIVEN those nu, so
    varying nu without refitting is a sensitivity of the PUBLISHED numbers, not
    a refit of the model.  It bounds how wrong [18] would have to be for the
    published parameter set to stop violating the criterion -- nothing more.
    """
    print("\n--- (F) how much does 1111 depend on the untraceable Ref.[18]? ---")
    l31, l32, nu1, nu2 = 14.520, 1.708, 0.4, 0.5
    for lab, lam, nu in (("(1,3)", l31, nu1), ("(2,3)", l32, nu2)):
        nu_crit = (LAM_II / lam) ** 2
        print("    %s  published nu = %.3f ;  criterion would be rescued at"
              " lam_ii = 0.4 only if nu < %.3e  (%.0f x smaller)"
              % (lab, nu, nu_crit, nu / nu_crit))
    print("    -> even a 9x error in nu2 (0.5 -> 0.055) leaves (2,3) failing,")
    print("       and (1,3) survives a ~500x error.  The 1111 verdict is")
    print("       insensitive to the unsourceable input; what IS sourced-limited")
    print("       is the exact violation factor (527x), not the verdict.")



if __name__ == "__main__":
    main()
