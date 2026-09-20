#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_essenberger_partition_2026-09-21.py -- 2026-09-21, ceiling-fe agent, FROM ZERO.

SETTLES THE BINARY the coordinator asked for on Essenberger et al.,
"Superconducting pairing mediated by spin fluctuations from first principles"
(FeSe), arXiv:1411.2121.

QUESTION AS POSED: is the 5-patch SF decomposition a REFINEMENT of the 2-band
(+/-) partition (route alive), or an INCOMPATIBLE partition (route dead)?

ANSWER: the question rests on a FALSE PREMISE.  There is no "lam^ph on 2 bands
vs lam^sf on 5 patches" mismatch in this paper at all.

L1 EVIDENCE, read BY ME this session from BOTH instruments
(text layer  pdftotext -layout -enc UTF-8  ->  _fe_1411_utf8.txt , AND
 300/200-dpi page crops  _fe1411_p2b-2.png , _fe1411_p2c-2.png , _fe1411_p4b-4.png
 rendered by me with pdftoppm -- the two instruments agree verbatim):

 p.2 (txt L127-134, image _fe1411_p2b-2.png):
   "First, to keep the formalism complexity at minimum, we consider a multiband
    isotropic approximation, meaning that states and pairings are averaged over
    Brillouin zone volumes on isoenergy surfaces. ...  This way we can group
    together the holes (see Fig. 1) around the Gamma-point (labeled with n = +)
    and the electrons at the M-point (labeled with n = -). ...  We will release
    this averaging approximation at the end of the paper."

 p.2 (txt L135-140, image _fe1411_p2c-2.png):
   "Under this approximation we account for: 1) the phononic interaction by a
    band-resolved Eliashberg function a2F_ij ... ; 2) the Coulomb repulsion by
    the average of the RPA screened matrix elements ... ; 3) the SF term
    containing an effective interaction mediated by paramagnon excitations
    Lam^SF as introduced in Ref. 51"
   => ALL THREE channels are built under the SAME 2-band (+/-) averaging.

 p.3 Fig. 3 caption (txt L179-186):
   "Two bands resolved superconducting pairing functions ...  a) electron-phonon
    Eliashberg function (PH); b) Screened Coulomb matrix elements ...; c) Spin
    fluctuation Eliashberg function (SF)"
   => the SAME "two bands resolved" label covers PH, Coulomb AND SF.

 p.4 (txt L253-262, image _fe1411_p4b-4.png):
   "We conclude the paper with an investigation of the gap function in k-space,
    going beyond the two-band isotropic approximation used so far (Eq. 2).  We
    divide the xy-plane of the Brillouin zone into sectors as shown in Fig. 5,
    and compute the superconducting gap corresponding to each sub-band. ...
    The large dip in the Delta_- is related to the CB fluctuation (M1 - M2) in
    this band ... not seen in Delta_+ (Gamma_1 - Gamma_2)."
   => the sectors appear ONLY in the final paragraph, as a post-hoc refinement
      of the WHOLE calculation (Tc moves "by a few percent"), and the sector
      labels are Gamma_1, Gamma_2, Gamma_{1,3} (inside the + band) and M_1, M_2
      (inside the - band) -- i.e. the sectors SUBDIVIDE the Gamma/M grouping.
      It is a refinement, not a rival partition, and it is not SF-specific.

 VERDICT: branch (b), a fortiori -- the partitions are not merely reconcilable,
 they are THE SAME partition in the body of the paper.  grep confirms that the
 strings "patch"/"sector"/"sub-band" occur NOWHERE before line 253 of 386.

THE NUMBERS THE PAPER ACTUALLY PRINTS (all read by me, L1):
 footnote [70] (txt L347-350):
   "The fact that phonons reduce the Tc in this material is a direct consequence
    of the stronger interband than the intraband coupling lam_ij (~0.03
    intraband and ~0.05 interband), effectively acting against the s+- ordering."
 body (txt L179-182): el-ph "integrated coupling lambda being less that 0.1"
 body (txt L195-205): SF "The maximum eigenvalue of the lam_ij matrix being
   lam_max = 0.07" at alpha/alpha_c = 0.9; "In the limit of alpha -> alpha_c,
   lam_max rises up to 0.48"; "SF are dominated by off diagonal components"
 body (txt L207-211): Tc SF-only 11 K; +Coulomb 11->4 K; +phonons 11->10 K;
   combined estimate 3 K; experiment 8 K.
 NOT PRINTED ANYWHERE: the individual lam^sf_ij matrix elements, the DOS N_i,
   and the SCDFT kernel definitions (deferred to Ref. 51).
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


# =============================================== 1. the margin, now NON-DEGENERATE
sep("1. FeSe IS THE FIRST TABLE IN THE CENSUS WITH OVERLAPPING SUPPORTS")
print("""
   On every previously examined s+- table (arXiv:0904.1808, arXiv:1304.3638)
   the two channels sat on DISJOINT matrix elements (ph only on the diagonal,
   sf only interband), so the (Z-dom) margin

       m_ij = Lam^Z_ij / |Lam^Del_ij| = (lam^ph_ij + lam^sf_ij)/|lam^ph_ij - lam^sf_ij|

   collapsed to the exact rational 1 and consumed ZERO degrees of freedom.

   FeSe is different: footnote [70] gives lam^ph on BOTH the diagonal (0.03)
   and the off-diagonal (0.05), while the SF matrix is "dominated by off
   diagonal components" with lam_max = 0.07 .. 0.48.  Both channels are
   therefore present on the SAME (+,-) element => r := lam^sf/lam^ph is a
   genuine finite positive number and the margin is a MATERIAL quantity.
""")
# lam^sf = [[a, s],[s, a]] with a ~ 0 (off-diagonal dominated) => lam_max = a+s
PH_INTRA, PH_INTER = F(3, 100), F(5, 100)
print("   %-26s %-12s %-12s %-12s %-12s"
      % ("SF regime", "lam^sf_+-", "r = sf/ph", "m_+- exact", "m_+-"))
for tag, sfoff in (("alpha/alpha_c = 0.9", F(7, 100)),
                   ("alpha -> alpha_c",    F(48, 100)),
                   ("zSe = 0.22 (Fig.4a)", F(130, 100))):
    Z = PH_INTER + sfoff
    D = PH_INTER - sfoff
    m = Z / abs(D)
    r = sfoff / PH_INTER
    print("   %-26s %-12s %-12s %-12s %-12.4f"
          % (tag, sfoff, r, m, float(m)))
print("""
   [rational-exact]  6 rationals / 3 additions / 3 subtractions / 6 divisions,
   zero floating point in the decision (floats above are display only).

   DEGREES OF FREEDOM CONSUMED HERE = 1 (the single ratio r = lam^sf/lam^ph on
   the interband element), versus 0 on every earlier table.  This is the first
   non-degenerate instance in the whole census.

   AND THE DOS CONVENTION CANCELS TOO: whichever index carries N_j, BOTH
   Lam^ph_ij and Lam^sf_ij carry the SAME factor, so it divides out of m_ij.
   No need to resolve the 0904-vs-1304 convention clash here either.  [analytic]
""")
print("   the diagonal element is still knife-edge, because SF is stated to be")
print("   off-diagonal dominated:  lam^ph_++ = 0.03, lam^sf_++ ~ 0  =>  m_++ = 1.")

# ============================== 2. BUT THE CRITERION IS STILL A TAUTOLOGY
sep("2. THE DECISIVE POINT: MORE DATA CANNOT RESCUE (P4)")
print("""
   Even with genuinely overlapping supports, (Z-dom) CANNOT FAIL under (K1).
   Proof (already in section 2.1 of the note, restated here because it is what
   kills the route):

       m = (1+r)/|1-r| >= 1  for every r in [0, +inf],   with equality iff
       r in {0, +inf}, and m -> +inf as r -> 1.

   Equivalently: |a - b| <= a + b for any a,b >= 0.  That is the triangle
   inequality and nothing else.  [analytic]

   So the arrow (P4) -- "does this a2F^ph / a2F^sf decomposition describe a real
   material?" -- is NOT untested for want of data.  It is UNTESTABLE IN
   PRINCIPLE under (K1): no measurement of lam^ph_ij and lam^sf_ij, however
   accurate and however ab initio, can produce a violation.  The canonical note
   already says this in one line ((K1) => "automatic, triangle inequality");
   what this session adds is that the literature search for a decisive material
   was therefore searching for something that does not exist.
""")
rng = np.random.default_rng(20260921)
worst = np.inf
for _ in range(200000):
    a, b = rng.random(2) * 10
    worst = min(worst, (a + b) / abs(a - b))
print("   brute-force confirmation, 200000 random (a,b) in [0,10]^2 :")
print("      min over samples of (a+b)/|a-b| = %.6f   (must be >= 1)" % worst)
print("   NOTE this random scan is NOT evidence for the theorem -- it cannot")
print("   reach the equality face {a=0 or b=0}, which is ATOMIC and has zero")
print("   probability under a continuous draw.  The equality face is exactly")
print("   what the published tables sit on.  Reported per CLAUDE.md sampling")
print("   discipline: the scan covers the INTERIOR only.  [numerical-evidence]")

# ================================= 3. WHERE mu* ACTUALLY ENTERS -- routing audit
sep("3. ROUTING AUDIT: mu* DOES NOT GO INTO THE (Z-dom) SPLIT")
print("""
   Section 2D of the note folded mu* into Lam^Del elementwise and checked the
   direction on the (1,1) element.  That was the NAIVE route.  Before writing
   any more about it I checked WHICH OPERATOR the repo's proof actually puts
   mu* into -- CLAUDE.md W1 lesson: check which variable the symbol acts on
   before recomputing.

   multiband_ceiling_lift_2026-09-18.md lines 242-247 and 341-342, and
   multiband_ceiling_CANONICAL_2026-09-18.md line 444:
       the Coulomb block enters as a Gram form  M_{(n,i),(m,j)} = 2 mu_ij
       w^(i)_n w^(j)_m ,  and  mu >= 0 (PSD)  =>  M >= 0  =>  O_mu = O - M <= O.
       "the ceiling holds verbatim, INDEPENDENTLY OF THE SIZE of mu*."
       "the condition is REVERSED relative to the boson layer: the boson layer
        needs a2F_ij elementwise nonnegative (PSD is not enough), the Coulomb
        layer needs mu* PSD (elementwise nonnegative is neither necessary nor
        sufficient)."

   => mu* is SUBTRACTED as a PSD block, OUTSIDE the (K1) split.  It never
   appears inside Lam^Del for the purposes of (Z-dom).  My section 2D is a
   robustness remark about a route the repo does not take, not a gap.

   WHY THE PSD ROUTE IS NECESSARY AND NOT COSMETIC -- the naive elementwise
   route really does break, and here is the exact threshold:
""")
print("""     naive route:  Lam^Del_ij = lam^ph_ij - lam^sf_ij - mu*_ij ,
                    Lam^Z_ij   = lam^ph_ij + lam^sf_ij     (mu* is static,
                                 it renormalises no mass, so it is absent from Z)
     |Lam^Del| <= Lam^Z
        <=>  -(ph+sf) <= ph - sf - mu* <= ph + sf
        right ineq:  -sf - mu* <= sf      <=>  mu* >= -2 sf   (always true)
        left  ineq:  ph - sf - mu* >= -ph - sf   <=>   mu* <= 2 lam^ph
     =>  (Z-dom), IF mu* WERE FOLDED IN ELEMENTWISE, holds at (i,j)
         IF AND ONLY IF   mu*_ij <= 2 lam^ph_ij .              [analytic]
""")
print("   numerical confirmation of the iff (boundary included), and the")
print("   contrast with the PSD route:")
for ph, sf, mus in ((0.05, 0.07, 0.00), (0.05, 0.07, 0.10),
                    (0.05, 0.07, 0.11), (0.03, 0.00, 0.06),
                    (0.03, 0.00, 0.07)):
    Z, D = ph + sf, ph - sf - mus
    ok = abs(D) <= Z + 1e-15
    pred = mus <= 2 * ph + 1e-15
    print("      ph=%.2f sf=%.2f mu*=%.2f -> |Del|/Z = %7.4f  naive %s "
          "| predicted by mu*<=2ph: %s  %s"
          % (ph, sf, mus, abs(D) / Z, "OK  " if ok else "FAIL",
             "OK  " if pred else "FAIL", "AGREE" if ok == pred else "MISMATCH"))
print("""
   PHYSICAL MEANING of the threshold: Z-channel dominance says the mass
   renormalisation must outweigh the pairing kernel.  Only RETARDED bosons
   renormalise the mass; the INSTANTANEOUS Coulomb repulsion contributes to
   pairing and nothing to Z.  So mu* is precisely the term that can break
   Z-dominance, and it breaks it once it exceeds twice the phonon attraction
   on that element.  For a weak-phonon s+- material (FeSe: lam^ph ~ 0.03-0.05,
   while a typical mu* is 0.1) the naive route would FAIL -- which is exactly
   why the repo does not take it.  The PSD route pays nothing for the SIZE of
   mu* and asks only for its SHAPE.
""")
print("   ONE CAVEAT I MUST STATE, because the PSD route is not free either:")
print("   a PURELY INTERBAND mu* is NOT PSD.  mu = [[0,c],[c,0]] has eigenvalues")
print("   +-c, so it fails (C2) for any c>0:")
for c in (0.05, 0.1, 0.3):
    ev = np.linalg.eigvalsh(np.array([[0.0, c], [c, 0.0]]))
    print("      mu = [[0,%.2f],[%.2f,0]] -> eigenvalues %+.3f, %+.3f -> NOT PSD"
          % (c, c, ev[0], ev[1]))
print("   A physical mu_ij = <rho_i, V rho_j> with Vhat >= 0 is PSD by Parseval")
print("   (lift note (mu2)) and is diagonally dominated -- and Essenberger's")
print("   Fig. 3b is read by the authors as 'diagonally dominated (++ and --")
print("   components)', consistent with PSD.  But the paper prints no numbers,")
print("   so PSD-ness of THIS mu is NOT VERIFIED, only plausible.  [not obtained]")

# ================================= 4. what would actually be needed
sep("4. WHAT IS STILL MISSING FROM 1411.2121, STATED HONESTLY")
print("""
   OBTAINED (L1, two instruments):
     - both channels live on one and the same 2-band (+/-) partition;
     - lam^ph_ij ~ 0.03 (intra) / 0.05 (inter), the FIRST non-fitted,
       non-circular, ab-initio lam^ph_ij matrix in this census;
     - lam^sf_max = 0.07 (alpha/alpha_c=0.9) .. 0.48 (alpha->alpha_c),
       1.3 at zSe=0.22, off-diagonal dominated;
     - Coulomb diagonally dominated;
     - Tc budget 11 K (SF only) -> 4 K (+Coulomb) -> 10 K (+phonons) -> 3 K.

   NOT OBTAINED:
     - the individual lam^sf_ij elements (only lam_max is printed; the matrix
       itself is only PLOTTED, Fig. 3c);
     - the per-band DOS N_i (needed for the ceiling's symmetrisation);
     - the SCDFT kernel definitions, which 1411.2121 defers to Ref. 51.
       *** SUPERSEDED WITHIN THIS SAME SCRIPT: Ref.51 = arXiv:1409.7968 is in
       the repo and DOES print them (Eq. 36).  See section 5.  I leave this
       line standing, struck through, rather than deleting it, because it
       records that "the paper does not print its equations" was true of
       1411.2121 alone and false of the paper it points to -- which is where
       the earlier census went wrong. ***
     - alpha is admitted by the authors to be non-first-principles, so the SF
       column carries one adjustable parameter after all.

   CONSEQUENCE FOR THE ROUTE: the partition obstruction is GONE (that reason
   was wrong), but the route dies anyway one level up, analytically, because
   (K1) makes (Z-dom) a tautology.  Getting Ref. 51 and digitising Fig. 3c
   would buy a sharper MARGIN, not a TEST.  Per CLAUDE.md's budget rule that
   is a purchase of slack, not of a verdict => do not dispatch that work.
""")

# ============================= 5. Ref. [51] = arXiv:1409.7968 -- I READ THE KERNELS
sep("5. REF.[51] DOES PRINT THE EQUATIONS.  MY INDEPENDENT READING OF Eq.(36)")
print("""
   The coordinator found that Ref.[51] of 1411.2121 is arXiv:1409.7968, the
   formalism sister paper (same authors, electron-gas model), and that its
   Eq.(36) prints both kernels.  I re-extracted that PDF MYSELF
   (pdftotext -layout -enc UTF-8 -> _ess1409_mine_utf8.txt, 1000 lines) and
   rendered MY OWN page crops at 200 dpi:
       _ess1409_eq36w-08.png  (Eq. 36, full-width display, p.8)
       _ess1409_stmt1-08.png  (numbered statement 1, p.8 right column)
       _ess1409_eq37-09.png   (Eq. 37, p.9 left column)

   WHAT THE PAGE IMAGE ACTUALLY SHOWS (and it differs from the transcription
   I was sent -- see the correction below):

     Delta^xc_k = -Delta^xc_k Z^D_k - SUM_k' K^C_kk' [tanh(beta E_k'/2)/(2E_k')] Delta^xc_k'

     Z^D_k    = (1/pi) SUM_k' INT_0^inf dw { Im[ (3/2) Lam^SF_kk'(w) ] - Lam^Ph_kk'(w) }
                        / (2 tanh(beta zeta_k/2))  *  (d/d zeta_k) J^+
     K^C_kk'  = (2/pi)        INT_0^inf dw { Im[ w_kk'(w) + (3/2) Lam^SF_kk'(w) ] + Lam^Ph_kk'(w) }
                        / (tanh * tanh)            *  J^-

   *** CORRECTION OF THE TRANSCRIPTION I WAS SENT: the Im[...] bracket CLOSES
   BEFORE Lam^Ph in BOTH kernels.  I was sent Im[(3/2)Lam^SF - Lam^Ph] and
   Im[w + (3/2)Lam^SF + Lam^Ph], i.e. with Lam^Ph INSIDE the bracket.  The page
   image puts it OUTSIDE.  This is not pedantry: it is what makes the paper's
   own sign statements self-consistent, since statement 1 asserts
   "Lam^Ph_kk'(w) < 0" (about Lam^Ph ITSELF, no Im) while Eq.(37) asserts
   "Im[Lam^SF] > 0 and Im[w] > 0" (about the imaginary parts).  With the
   bracket where the transcription put it, statement 1 would be about a
   quantity that never appears.  ***

   *** BUT THE ALGEBRA IS UNCHANGED, so the transcription error is NOT
   consequential for the verdict.  I verify that explicitly: ***
""")
print("""     define the three NONNEGATIVE quantities fixed by the authors themselves
     (statement 1, verified on _ess1409_stmt1-08.png; Eq.37, verified on
      _ess1409_eq37-09.png; both with FREE indices k,k' => ELEMENTWISE):
         a := -Lam^Ph_kk'(w)      >= 0     (statement 1: Lam^Ph < 0, attractive)
         s := (3/2) Im Lam^SF_kk'(w) >= 0  (Eq. 37)
         c := Im w_kk'(w)         >= 0     (Eq. 37)
     then
         Z^D numerator  = Im[(3/2)Lam^SF] - Lam^Ph  =  s + a
         K^C numerator  = Im[w + (3/2)Lam^SF] + Lam^Ph = c + s - a
     -- identical to what the (wrong-bracket) transcription gives.  Confirmed.

     (the 3/2 carries the SAME coefficient in both kernels, so it merely
      rescales s and cannot affect any ratio.  Confirmed on the page image.)
""")
print("""   *** ALSO VERIFIED ON THE PAGE IMAGE, AND LOAD-BEARING: w_kk' APPEARS IN
   K^C AND NOT IN Z^D. ***  This is the first-hand confirmation of the physical
   statement I used in section 3 above: the Coulomb term renormalises no mass
   and therefore contributes to the pairing channel only.  I had asserted it
   from the structure of the Ummarino equations; here it is visible directly in
   the SCDFT kernels of a different formalism.  [fact-check, two instruments]

   AND in 1411.2121 p.3 (image _fe1411_fig3a-3.png) the authors state
   "For the Coulomb interaction we use a STATIC screening", matching the third
   line of their approximation scheme  w_ij(E,E',w) ~ w_ij(E,E',0)  =>  for
   FeSe the Coulomb block is STATIC, i.e. exactly the object the repo's (C2)
   layer and the Anderson-Morel lemma (mu3) were derived for.
""")

sep("6. WHERE I DISAGREE WITH THE RECLASSIFICATION, AND BY HOW MUCH")
print("""
   WHAT I CONFIRM.  Taking only the two BOSONIC pieces, Eq.(36) gives
       Z channel   :  a + s        (two nonneg pieces ADDED)
       pair channel:      s - a    (phonon SUBTRACTED)
   which IS the (K1) add/subtract convention with Lam^ph := a, Lam^sf := s,
   both elementwise nonnegative by the authors' own statements.  The repo
   already has (K1) => (Z-dom) automatic (triangle inequality, analytic).
   So for the BOSON layer the reclassification is right and I verify it.
   GRADE: convention reading = fact-check (two instruments); the implication
   = analytic (pre-existing).  NOT to be abbreviated as "analytic".

   WHAT I DO NOT CONFIRM.  The message says the pairing channel is (K1) "only
   with one extra Coulomb term c (equally nonnegative, equally entering only
   the pairing channel)".  That parenthesis is exactly the step that must not
   be waved through: c is the ONE term that can break (Z-dom), because it adds
   to the pairing channel WITHOUT adding to Z.  Fold it in elementwise and the
   repo's Lam^Del becomes  a - s - c , so

       (Z-dom):  |a - s - c| <= a + s
                 upper branch: always true
                 lower branch: <=>  c <= 2a          [analytic, iff]

   This is NOT automatic, and it is the same threshold derived in section 3
   above from the Ummarino-side equations -- two independent routes, one
   number.  For FeSe it is very likely VIOLATED: the authors report that the
   Coulomb term cuts Tc from 11 K to 4 K (the single largest effect in their
   budget) while the phonon coupling is "less that 0.1" and cuts 11 K to 10 K.
""")
print("   quantified with the paper's own Tc budget as a proxy for relative size:")
for name, dTc in (("Coulomb: 11 K -> 4 K", 11 - 4), ("phonons: 11 K -> 10 K", 11 - 10)):
    print("      %-24s  Tc reduction %2d K" % (name, dTc))
print("      ratio of the two reductions = %.1f  =>  c/a is plausibly >> 2"
      % (7.0 / 1.0))
print("      (this is a PROXY, not a measurement of c/a: Tc depends on the")
print("       kernels non-linearly.  Graded numerical-evidence, indicative only.")
print("       The paper prints no numbers for w_ij, so c/a is NOT OBTAINED.)")
print("""
   THE RESOLUTION -- and it is the repo's, not mine.  The repo does NOT fold
   the Coulomb block into the elementwise split at all.  It subtracts it as a
   PSD Gram block ((mu1)/(C2) of the lift note): mu >= 0 => O_mu = O - M <= O,
   "independently of the SIZE of mu*".  Under that route c never meets the
   2a threshold and the boson-layer (Z-dom) verdict stands untouched.  What
   the Coulomb layer then needs instead is that Im[w_ij] be PSD as a matrix in
   the band indices -- and Eq.(37) gives only ELEMENTWISE positivity, which the
   lift note states explicitly is "neither necessary nor sufficient" for PSD.

   => MY VERDICT ON FeSe / Essenberger:
        boson layer   : A  (via (K1); convention read first-hand, two instruments)
        Coulomb layer : UNDETERMINED (needs Im[w_ij] PSD; only elementwise
                        positivity is printed, and no numbers at all)
      NOT plain "A (conditional)".  The two layers have different consumption
      points and must be graded separately.
""")
print("   the gap between elementwise-positive and PSD, made concrete on 2x2:")
for lab, M in (("diag dominated (authors' Fig.3b reading)", [[0.30, 0.10], [0.10, 0.30]]),
               ("borderline",                                [[0.30, 0.30], [0.30, 0.30]]),
               ("off-diag dominated",                        [[0.10, 0.30], [0.30, 0.10]]),
               ("purely interband",                          [[0.00, 0.20], [0.20, 0.00]])):
    A = np.array(M)
    ev = np.linalg.eigvalsh(A)
    print("      %-42s elementwise>=0: yes | eig %+.3f %+.3f -> %s"
          % (lab, ev[0], ev[1], "PSD" if ev[0] >= -1e-15 else "NOT PSD"))
print("   => all four are elementwise nonnegative; two of them are not PSD.")
print("      Elementwise positivity therefore cannot certify the Coulomb layer.")
print("      Essenberger's own words 'diagonally dominated (++ and -- components)'")
print("      put FeSe in the first row, which IS PSD -- but that is a reading of")
print("      a figure by its authors, not a number.  Hence UNDETERMINED.")
