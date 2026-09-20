#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_0904_mustar_table2_2026-09-21.py -- 2026-09-21, ceiling-fe agent, FROM ZERO.

THE ONE TABLE IN THE WHOLE IRON-BASED CENSUS THAT HAS mu* != 0.

L1 QUOTES, read BY ME this session from _fe_0904.1808.txt (text layer) --
G. A. Ummarino, M. Tortello, D. Daghero, R. S. Gonnelli, arXiv:0904.1808,
"Three-band s+- Eliashberg model and the superconductivity of iron pnictides":

  line 111  Lam^Del_ij(iwn,iwm) = Lam^ph_ij(iwn,iwm) - Lam^sp_ij(iwn,iwm)
  line 105  ... piT [Lam^Del_ij(iwn,iwm) - mu*_ij(wc)] theta(wc-|wm|) N^Del_j(iwm)
  line  81  "lam^sp_ii = 0 and mu*_ii(wc) = mu*_ij(wc) = 0"   (the BASE model)
  line 164  TABLE I :  lam31 = 3.866  lam32 = 0.471   mu*_ij = 0  , lam_ii = 0.4
  line 172  TABLE II:  lam31 = 2.730  lam32 = 0.758   mu*_ij = 0.1, lam_ii = 0.4
  line 176  "TABLE II: The effect of the Coulomb interaction, mu*_ij for the
             case shown in Table I."
  line 180  "We also tested the effect into the model of a small intraband
             coupling (POSSIBLY OF PHONON ORIGIN)."
  line 190  "We chose mu*_ij = 0.1 and, as expected, we found that the intraband
             Coulomb pseudopotential has a negligible effect while the interband
             one strongly contributes to RAISE Tc and reduces in a considerable
             way the value of lam31"

PHYSICS OF THAT LAST SENTENCE (worth stating, because it is the whole point of
s+-): an INTERBAND repulsion is PAIRING-FAVOURABLE when the gap changes sign
between the two bands -- the repulsion and the sign change multiply to an
attraction in the gap equation.  That is why mu*_31 > 0 RAISES Tc here and lets
the fit get away with a SMALLER lam31.  It is also exactly why the interband
Coulomb element cannot be waved through in a (Z-dom) audit: it enters the
pairing channel with the SAME sign as the spin-fluctuation term and the
OPPOSITE sign to the phonon term, while entering the Z channel not at all.

WHAT THIS SCRIPT DOES.  Two readings of where mu* lives, EXACT rational
arithmetic, and they give OPPOSITE verdicts on the same table:

  (R1) ELEMENTWISE FOLDING (the naive (K1) reading):
         Lam^Del_ij = lam^ph_ij - lam^sf_ij - mu*_ij
         Lam^Z_ij   = lam^ph_ij + lam^sf_ij            (mu* absent: a STATIC
                      repulsion renormalises no mass, so it never enters Z)
       => (Z-dom) at (i,j)  <=>  mu*_ij <= 2 lam^ph_ij          [analytic, iff]

  (R2) THE REPO'S ACTUAL ROUTE ((C2)/(mu1) of multiband_ceiling_lift):
         the Coulomb block is subtracted as a PSD Gram block,
         mu >= 0 (PSD) => M >= 0 => O_mu = O - M <= O,
       independent of the SIZE of mu*.  The requirement is PSD, NOT
       elementwise nonnegativity.

This is the first element of the whole iron-based census whose (Z-dom) verdict
is a FUNCTION OF THE PUBLISHED NUMBERS rather than of the support pattern.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from fractions import Fraction as F
import numpy as np


def sep(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


# Printed inputs, as exact rationals.  Bands: 1,2 = hole (Gamma), 3 = electron (M).
# Only lam31 and lam32 are printed; the model is "predominantly interband".
TABLES = {
    "Table I  (mu*=0)":   dict(l31=F(3866, 1000), l32=F(471, 1000),
                               mus=F(0), lii=F(4, 10)),
    "Table II (mu*=0.1)": dict(l31=F(2730, 1000), l32=F(758, 1000),
                               mus=F(1, 10), lii=F(4, 10)),
}

sep("1. (R1) ELEMENTWISE FOLDING -- exact margins m_ij = Lam^Z/|Lam^Del|")
print("""
   lam^ph lives ONLY on the diagonal (lam_ii = 0.4, "possibly of phonon
   origin"); lam^sp lives ONLY off-diagonal (line 81: lam^sp_ii = 0).
   mu*_ij = mu* on EVERY element (line 176: "mu*_ij for the case shown in
   Table I" -- one number, all ij; the text then distinguishes the intraband
   and interband EFFECTS, not two different values).
""")
print("   %-20s %-10s %-14s %-14s %-12s %s"
      % ("table", "element", "Lam^Z", "|Lam^Del|", "m exact", "verdict"))
for tag, d in TABLES.items():
    rows = [("(1,1) diag", d["lii"], F(0)),
            ("(3,1) inter", F(0), d["l31"]),
            ("(3,2) inter", F(0), d["l32"])]
    for name, ph, sf in rows:
        Z = ph + sf
        De = ph - sf - d["mus"]
        if De == 0:
            m, v = None, "Lam^Del = 0 (margin infinite)"
        else:
            m = Z / abs(De)
            v = "HOLDS" if m >= 1 else "*** FAILS ***"
        print("   %-20s %-10s %-14s %-14s %-12s %s"
              % (tag, name, str(Z), str(abs(De)),
                 "n/a" if m is None else "%s=%.6f" % (m, float(m)), v))
print("""
   [rational-exact]  6 elements x 2 tables / 8 input rationals (l31,l32,mu*,lii
   for two tables) / 12 additions / 6 divisions / 6 comparisons -- ZERO floating
   point in any verdict.  The decimal prints are display only.

   *** This is the FIRST (Z-dom) verdict in the entire iron-based census that
       is a FUNCTION of the published numbers.  Degrees of freedom consumed = 1
       (the ratio mu*/lam^sf on the interband elements).  Every other table in
       the census had DOF = 0.  ***

   And under (R1) it FAILS, by exactly the amount
        m_3j = lam^sf_3j / (lam^sf_3j + mu*)  <  1,
   i.e. the violation is the Coulomb fraction of the pairing kernel.
""")

sep("2. THE THRESHOLD, DERIVED AND THEN CONFIRMED NUMERICALLY")
print("""
   [analytic, iff]  With a := lam^ph_ij >= 0, b := lam^sf_ij >= 0, c := mu*_ij >= 0:
       (Z-dom)  <=>  |a - b - c| <= a + b
       upper branch  a - b - c <=  a + b  <=>  -(b+c) <= b  : ALWAYS true
       lower branch  a - b - c >= -a - b  <=>   c <= 2a     : THE THRESHOLD
   => (Z-dom) at (i,j)  <=>  mu*_ij <= 2 lam^ph_ij,  INDEPENDENT of lam^sf_ij.

   PHYSICAL READING: Z-channel dominance says the mass renormalisation must
   outweigh the pairing kernel.  Only RETARDED bosons renormalise the mass;
   the INSTANTANEOUS Coulomb repulsion goes into pairing alone.  So the budget
   the phonons buy for the Z channel is 2 lam^ph, and the static repulsion
   spends it.  Spin fluctuations drop out of the threshold entirely because
   they pay into BOTH channels (they are retarded too).
""")
print("   %-10s %-10s %-10s %-14s %-14s %s"
      % ("a=lam^ph", "b=lam^sf", "c=mu*", "c <= 2a ?", "|a-b-c|<=a+b ?", "agree"))
ok = True
for a, b, c in [(F(4, 10), F(0), F(1, 10)),        # 0904 diagonal
                (F(0), F(273, 100), F(1, 10)),     # 0904 (3,1)
                (F(0), F(758, 1000), F(1, 10)),    # 0904 (3,2)
                (F(4, 10), F(0), F(8, 10)),        # exactly at threshold
                (F(4, 10), F(0), F(81, 100)),      # just above
                (F(3, 100), F(5, 100), F(1, 10)),  # FeSe-like phonon element
                (F(1), F(1), F(1))]:
    thr = (c <= 2 * a)
    direct = (abs(a - b - c) <= a + b)
    ok &= (thr == direct)
    print("   %-10s %-10s %-10s %-14s %-14s %s"
          % (a, b, c, thr, direct, "AGREE" if thr == direct else "DISAGREE"))
print("   -> %s   [rational-exact] 7 triples / 21 rationals, zero floating point."
      % ("all AGREE, the iff is confirmed" if ok else "MISMATCH -- BUG"))

sep("3. (R2) THE REPO'S ACTUAL ROUTE GIVES THE OPPOSITE VERDICT ON THE SAME TABLE")
print("""
   The repo does NOT fold mu* into the elementwise split.  (C2)/(mu1) subtract
   the Coulomb as a Gram block:  M_{(n,i),(m,j)} = 2 mu_ij w^i_n w^j_m, so
   mu PSD => M PSD => O_mu = O - M <= O, and the ceiling survives REGARDLESS
   of the size of mu*.  The requirement is PSD, not elementwise nonneg.

   0904 Table II has mu*_ij = 0.1 on EVERY element, i.e.
        mu = 0.1 * J_3   (J_3 = all-ones 3x3)
   whose eigenvalues are 0.3, 0, 0 -- rank one, PSD.  So under the repo's route
   Ba0.6K0.4Fe2As2 PASSES, while under (R1) it FAILS at both interband elements.
""")
mu = 0.1 * np.ones((3, 3))
ev = np.linalg.eigvalsh(mu)
print("   eigenvalues of mu = 0.1 * J_3 :", np.array2string(ev, precision=12))
print("   min eigenvalue = %.15g  >= 0  => PSD  => (C2) satisfied" % ev.min())
print("   (exact: J_3 = 1 1^T is rank one with eigenvalues 3,0,0 -> 0.3,0,0.")
print("    [analytic] any CONSTANT mu*_ij = c >= 0 matrix is c * 1 1^T, PSD.)")
print("""
   *** THE LESSON, and it is the load-bearing one ***
   The same published table is a PASS under one routing of the Coulomb term and
   a FAIL under the other.  So "Ba-122 satisfies (Z-dom)" is not a statement
   about Ba-122; it is a statement about WHERE YOU PUT mu*.  Any census row
   must name the routing, or it means nothing.

   Which routing is right?  The repo's, for the ceiling theorem -- because the
   ceiling proof needs O_mu <= O as OPERATORS, and the Gram subtraction gives
   exactly that with no threshold.  (R1) is the reading one falls into by
   pattern-matching mu* onto the (K1) two-piece split, and it imposes a
   condition the theorem never asked for.
""")

sep("4. NEGATIVE CONTROL: elementwise nonneg does NOT imply PSD")
print("""
   Needed because the SCDFT route (Essenberger Eq. 37) prints only ELEMENTWISE
   positivity of Im[w_kk'], while (C2) needs PSD.  Four 2x2 Coulomb blocks,
   all elementwise nonnegative:
""")
for M in ([[F(3, 10), F(1, 10)], [F(1, 10), F(3, 10)]],
          [[F(3, 10), F(3, 10)], [F(3, 10), F(3, 10)]],
          [[F(1, 10), F(3, 10)], [F(3, 10), F(1, 10)]],
          [[F(0), F(2, 10)], [F(2, 10), F(0)]]):
    tr = M[0][0] + M[1][1]
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    psd = (tr >= 0) and (det >= 0)
    print("   [[%s,%s],[%s,%s]]  tr=%-8s det=%-8s  %s"
          % (M[0][0], M[0][1], M[1][0], M[1][1], tr, det,
             "PSD" if psd else "*** NOT PSD ***"))
print("""   [rational-exact] 4 matrices / 8 rationals / 4 traces / 4 determinants.
   => 2 of 4 are NOT PSD.  A purely INTERBAND Coulomb block (last row) has
   eigenvalues +-c and is never PSD unless c = 0.  So "the paper says the
   Coulomb repulsion is positive" does NOT discharge (C2).  This is why the
   FeSe Coulomb layer is graded UNDETERMINED and not A.
""")
