#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_lifeas_intel3_2026-09-21.py -- 2026-09-21, ceiling-fe agent, FROM ZERO.

Answers the coordinator's THIRD intel batch on arXiv:1304.3638 (LiFeAs).
Everything below was re-read BY ME from the UTF-8 text layer
(_fe_1304_utf8.txt, produced with  pdftotext -layout -enc UTF-8).

VERBATIM L1 QUOTES I READ THIS SESSION (line numbers of _fe_1304_utf8.txt):

 L86-97  Eq.(1):  wn Zi(iwn) = wn + piT sum_{m,j} Lam^Z_ij(iwn,iwm) N^Z_j(iwm)
                                    + sum_j [Gam_ij + Gam^M_ij] N^Z_j(iwn)
         *** NO Theta in the Z equation ***
 L99-103 Eq.(2):  Zi(iwn) Del_i(iwn) = piT sum_{m,j} [Lam^Del_ij(iwn,iwm)
                        - mu*_ij(wc)] x Theta(wc - |wm|) N^Del_j(iwm) + ...
         *** the Theta multiplies the WHOLE bracket, Lam^Del INCLUDED ***
 L118    Lam^{ph,sf}_ij(iwn,iwm) = 2 int_0^{+infinity} dW W a2F^{ph,sf}_ij(W)
                                      / [(wn-wm)^2 + W^2]
         *** the SPECTRAL integral runs to +infinity, it is NOT truncated ***
 L141-3  "we should take: lam^ph_ii = lam^ph_ij = 0,  mu*_ii(wc) = mu*_ij(wc) = 0
          i.e. the electron-phonon coupling constant and the Coulomb
          pseudopotential COMPENSATE EACH OTHER, in first approximation"
 L190    "The cut-off energy is wc = 18 W0^sf and the maximum QUASIPARTICLE
          energy is wmax = 21 W0^sf"
 L186-9  "In all the calculations we set W_ij = W^sf_ij = W^sf_0 = 8 meV,
          and Y_ij = Y^sf_ij = W^sf_ij/2"
 L72-84  Table 1 "Fermi Surface resolved Kohn Sham properties":
          N(0) in states/spin/eV:  N1=0.556 N2=0.646 N3=0.616 N4=0.370 N5=0.039
 L194-6  "nu13 = 0.9019, nu14 = 1.5010, nu23 = 1.0483, nu24 = 1.7447"
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from fractions import Fraction as F
import numpy as np
from scipy.integrate import quad


def sep(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


# ============================================================ 1. GENERATING FORMULA
sep("1. THE GENERATING FORMULA OF THE nu MATRIX IS PRINTED IN THE SAME PAPER")
print("""
   The coordinator verified rank-one through the CONSISTENCY relation
   nu13*nu24 == nu14*nu23.  That is a necessary condition.  Table 1 gives
   something strictly stronger: the four N_i(0) themselves, i.e. the actual
   GENERATING FORMULA.  If nu_ij = N_i/N_j really generates the printed nu,
   then the four printed nu are not four measurements but ONE vector of
   four DOS -- and, being ratios, only THREE independent numbers.
""")
N = {1: F(556, 1000), 2: F(646, 1000), 3: F(616, 1000), 4: F(370, 1000)}
PRINTED = {(1, 3): F(9019, 10000), (1, 4): F(15010, 10000),
           (2, 3): F(10483, 10000), (2, 4): F(17447, 10000)}
print("   %-8s %-14s %-14s %-12s %-12s" %
      ("nu_ij", "from Table 1", "printed", "rel.diff", "rounding hw"))
ok = True
for (i, j), p in PRINTED.items():
    gen = N[i] / N[j]
    rel = abs(gen - p) / p
    # N printed to 3 decimals => absolute half-width 5e-4 on each N
    hw = F(5, 10000) / N[i] + F(5, 10000) / N[j]
    flag = "OK" if rel <= hw else "EXCEEDS"
    ok &= (rel <= hw)
    print("   nu%d%d    %-14.6f %-14.4f %-12.2e %-12.2e %s"
          % (i, j, float(gen), float(p), float(rel), float(hw), flag))
print("   VERDICT: all four printed nu are reproduced by N_i/N_j from Table 1")
print("            within the 3-decimal rounding of N.  -> %s" % ("CONFIRMED" if ok else "FAILED"))
print("   [rational-exact]  8 rationals (4 DOS + 4 nu) / 4 divisions / 4")
print("   comparisons / 8 divisions for the half-widths -- zero floating point")
print("   in the decision.  The float prints above are display only.")
print()
print("   => the nu matrix is rank-one BY CONSTRUCTION, not by coincidence,")
print("      and the material content of the four nu is the DOS ratio vector")
print("      N1:N2:N3:N4 = %s : %s : %s : %s  (3 independent numbers)."
      % tuple(str(N[k] / N[1]) for k in (1, 2, 3, 4)))
print("   This is INDEPENDENT of, and stronger than, the nu13*nu24=nu14*nu23")
print("   check: that test can only falsify rank-one, this one exhibits the")
print("   generator.  Exact cross-check of the consistency relation itself:")
a, b = PRINTED[(1, 3)] * PRINTED[(2, 4)], PRINTED[(1, 4)] * PRINTED[(2, 3)]
print("      nu13*nu24 - nu14*nu23 = %s  (rel %.3e)"
      % (a - b, float(abs(a - b) / ((a + b) / 2))))
print("      exact from Table 1     = %s   (identically zero, as it must be)"
      % (N[1] / N[3] * (N[2] / N[4]) - N[1] / N[4] * (N[2] / N[3])))

# ====================================================== 2. HOMOGENEITY DEGREE AUDIT
sep("2. HOMOGENEITY AUDIT -- correcting the coordinator's claim")
print("""
   The coordinator wrote: "M2^sf_ij = <w^2>_sf * lam^sf_ij, the energy scale
   <w^2>_sf is a common factor => in any criterion HOMOGENEOUS in M2 it
   cancels EXACTLY."

   The first half is right and I confirm it below.  The second half is NOT:
   it silently assumes degree ZERO.  Count the degrees.

     ceiling  = C_inf sqrt(rho(M2''_Z))           -- degree +1/2 in M2
                                                     = degree +1 in the energy
     Tc       (solution of lam_1(O(T)) = 1)       -- degree +1 in the energy
     Tc/ceiling                                   -- degree  0   => CANCELS

   So <w^2>_sf does NOT cancel from the CEILING (it cannot: the ceiling is an
   energy and must carry an energy scale).  It cancels from the RATIO, and it
   cancels there only because Tc is itself homogeneous of degree 1 -- which
   requires that the model contain NO OTHER energy scale.

   *** And THAT is exactly what the coordinator's intel item (1) buys. ***
   mu*_ij == 0 removes the one dimensionless-coupling-with-its-own-cutoff that
   would have broken the homogeneity, and wc = 18 W0^sf, wmax = 21 W0^sf are
   TIED TO W0 rather than being independent energies.  With mu* = 0 and both
   cutoffs proportional to W0, the ONLY energy in the whole model is W0.
   => Tc/ceiling is EXACTLY scale free, not approximately.        [analytic]

   PHYSICAL MEANING: the LiFeAs fit cannot be used to test the ceiling,
   because the ceiling and Tc are the same energy scale times two different
   pure functions of (lam_ij, Y/W0).  The test has one equation and zero
   material inputs on the energy axis.
""")


def shape(W, W0, Y):
    return 1.0 / ((W - W0) ** 2 + Y ** 2) - 1.0 / ((W + W0) ** 2 + Y ** 2)


print("   numerical confirmation that the shape factors out (same (W0,Y) for")
print("   every (i,j) => M2_ij / lam_ij is the SAME number for every element):")
W0, Y = 8.0, 4.0
num = quad(lambda W: 2 * W * shape(W, W0, Y), 0, np.inf, limit=800)[0]
den = quad(lambda W: 2 * shape(W, W0, Y) / W, 0, np.inf, limit=800)[0]
print("      M2/lam = %.10f meV^2   vs   W0^2+Y^2 = %.10f meV^2  (rel %.2e)"
      % (num / den, W0 ** 2 + Y ** 2,
         abs(num / den - (W0 ** 2 + Y ** 2)) / (W0 ** 2 + Y ** 2)))

# ============================================== 3. IS M2 TRUNCATION-DOMINATED?
sep("3. THE COORDINATOR'S THIRD PIT: is M2 dominated by the cutoff?")
print("""
   The worry: "a Lorentzian's second moment is cutoff-dominated (tail ~1/w^2)".

   TWO SEPARATE QUESTIONS, and CLAUDE.md's W1 lesson #4 says to answer them in
   this order -- first WHICH VARIABLE the cutoff acts on, only then recompute.

   (a) WHERE DOES wmax ACT?  L1, line 118 of the text layer: the spectral
       integral defining Lam is printed  int_0^{+infinity} dW ... -- it is NOT
       truncated.  Line 190 calls wmax "the maximum QUASIPARTICLE energy", and
       line 99 shows the only Theta in the paper sits in Eq.(2) as
       Theta(wc - |wm|), whose argument is a MATSUBARA FREQUENCY wm, not W.
       => neither wc nor wmax is a spectral cutoff.  The premise of the worry
       does not hold for this paper.  (This is verbatim the accident recorded
       in CLAUDE.md: reading "w_c" as a SPECTRAL cutoff -- do not infer
       semantics from a name.)

   (b) EVEN SO, HOW BIG WOULD IT BE?  Answer it anyway, because "the premise is
       false" is a weaker deliverable than "and here is the number".
       The key point the worry misses: a2F here is a DIFFERENCE of two
       Lorentzians.  L(W-W0) - L(W+W0) ~ 4 W0 / W^3 at large W, one power
       FASTER than a single Lorentzian's W^-2.  So int W a2F dW converges like
       int W^-2, and the tail beyond wmax is a FINITE, SMALL fraction.
""")
tot_num = quad(lambda W: W * shape(W, W0, Y), 0, np.inf, limit=800)[0]
tot_den = quad(lambda W: shape(W, W0, Y) / W, 0, np.inf, limit=800)[0]
print("   %-14s %-14s %-14s %-12s %-12s"
      % ("W_cut (meV)", "M2/lam", "vs untrunc.", "lam deficit", "ceiling x"))
for wc in (84.0, 144.0, 168.0, 336.0, 1e4):
    n = quad(lambda W: W * shape(W, W0, Y), 0, wc, limit=800)[0]
    d = quad(lambda W: shape(W, W0, Y) / W, 0, wc, limit=800)[0]
    lab = "%.0f" % wc if wc < 1e3 else "1e4 (~inf)"
    print("   %-14s %-14.6f %-14.6f %-12.2e %-12.6f"
          % (lab, n / d, (n / d) / (tot_num / tot_den),
             1 - d / tot_den, np.sqrt((n / d) / (tot_num / tot_den))))
print("""
   READING: the untruncated value is exactly W0^2+Y^2 = 80.000 meV^2.  At the
   paper's own wmax = 168 meV the second moment is 77.574 = 96.97% of it.
   Relative to THAT reference value, HALVING the cutoff to 84 meV moves M2 by
   -3.1% and DOUBLING it to 336 meV by +1.6%.  The ceiling moves by the square
   root, i.e. -1.6% / +0.8%.  The lambda integral is insensitive at the 1e-5
   level (its integrand decays like W^-4).
   => M2 is NOT truncation-dominated here.  The worry would have been correct
      for a SINGLE Lorentzian; the difference form removes it.        [checked]
   DIRECTION: truncation LOWERS M2 => LOWERS the ceiling => makes the test
   HARDER.  So the untruncated ceiling I report is the PERMISSIVE choice and
   the reported "consistent" verdicts stay sound under either convention.
""")
# same check for a SINGLE Lorentzian, to show the contrast (negative control)
print("   NEGATIVE CONTROL -- single Lorentzian L(W-W0,Y) alone, which is what")
print("   the worry actually describes:")
sng = lambda W: 1.0 / ((W - W0) ** 2 + Y ** 2)
for wc in (84.0, 168.0, 336.0, 3360.0):
    n = quad(lambda W: W * sng(W), 0, wc, limit=800)[0]
    d = quad(lambda W: sng(W) / W, 0, wc, limit=800)[0]
    print("      W_cut=%7.0f : M2/lam = %12.3f meV^2  (grows without bound"
          " ~ 2 W_cut Y/pi)" % (wc, n / d))
print("   => the contrast is the whole point: the odd difference converges,")
print("      the single Lorentzian does not.  The paper uses the difference.")

# ========================================== 4. (Z-dom-k) AS EXPLICIT ALGEBRA
sep("4. (Z-dom-k) WRITTEN OUT AS EXPLICIT ALGEBRA + THE CANCELLATION LIST")
print("""
   Inputs actually printed by the paper, all 9 of them:
       fitted   : lam11, lam13, lam23, lam14, lam24          (5)
       1st-princ: nu13, nu14, nu23, nu24 = N_i/N_j           (4 printed, 3 free)
       hand-set : lam^ph_ij = 0 (i != j), lam^sf_ii = 0, mu* = 0,
                  W_ij = 8 meV and Y_ij = 4 meV for ALL (i,j)

   Because every (i,j) shares one shape, Lam^X_ij(k) = lam^X_ij * J(k) with
       J(k) := [int 2W D(W)/(W^2+(2 pi T k)^2) dW] / [int 2 D(W)/W dW],
   ONE scalar function of k, the SAME for every band pair.        [analytic]

   Therefore
       Lam^Z_ij(k)   = (lam^ph_ij + lam^sf_ij) J(k)
       Lam^Del_ij(k) = (lam^ph_ij - lam^sf_ij) J(k)
   and the (Z-dom-k) margin is

       m_ij(k) = |Lam^Z_ij(k)| / |Lam^Del_ij(k)|
               = (lam^ph_ij + lam^sf_ij) / |lam^ph_ij - lam^sf_ij|

   *** J(k) CANCELS IDENTICALLY -- k, T, W0, Y all gone. ***
   The criterion cannot see: the energy scale W0, the width Y, the temperature,
   the Matsubara index, the cutoffs wc and wmax.  [analytic]  This is the
   coordinator's predicted cancellation, CONFIRMED, and it is wider than
   predicted (it also eats k and T, not just the energy scale).

   Then the supports: lam^ph lives ONLY at (1,1), lam^sf ONLY off-diagonal.
   At every element one of the two is exactly 0, so

       m_ij(k) = x/|x| = 1   for every i,j,k and for ANY values of the 9 numbers.

   *** FULL CANCELLATION LIST (what the (Z-dom-k) verdict is a function of) ***
       lam11  : cancels      lam13 : cancels      lam23 : cancels
       lam14  : cancels      lam24 : cancels
       nu13   : cancels      nu14  : cancels      nu23  : cancels
       nu24   : cancels
       W0, Y, T, k, wc, wmax : cancel
       -> the verdict is a function of NOTHING.  Degrees of freedom = 0.
   The only input it does depend on is the SUPPORT PATTERN, i.e. the modelling
   choice of which channel is allowed on which matrix element.
""")

# --------- the one place where the channel assignment is itself questionable
sep("5. THE ONE SUBTLETY IN THE CHANNEL ASSIGNMENT, AND ITS DIRECTION")
print("""
   L1 line 141-143: the paper sets lam^ph = 0 AND mu* = 0 together, because
   "the electron-phonon coupling constant and the Coulomb pseudopotential
   COMPENSATE EACH OTHER, in first approximation".  When it later reinstates
   lam11 != 0, that lam11 is therefore a NET effective intraband coupling,
   not a bare phonon coupling: lam11 ~ lam^ph_11 - (something mu*-like).

   DIRECTION AUDIT of what that does to the margin at (1,1):
     if lam11 = lam^ph_11 - mu*_11 with mu*_11 >= 0, then
        Lam^Del_11 = lam^ph_11 - mu*_11 = lam11        (mu* enters Del with a -)
        Lam^Z_11   = lam^ph_11          = lam11 + mu*_11   (mu* does NOT enter Z)
     => m_11 = (lam11 + mu*_11)/lam11 >= 1, INCREASING in mu*_11.
   So the compensation reading can only make (Z-dom) MORE comfortable at (1,1),
   never less.  The direction is safe, and the verdict is unchanged.  But note
   the consequence for the OTHER direction: it means the "phonon" label on
   lam11 is an effective one, so this table cannot be used as evidence about
   lam^ph_11 itself either.
""")
for mu11 in (0.0, 0.05, 0.10, 0.20):
    for l11 in (0.86, 0.90):
        print("      lam11 = %.2f, mu*_11 = %.2f  ->  m_11 = %.4f"
              % (l11, mu11, (l11 + mu11) / l11))
