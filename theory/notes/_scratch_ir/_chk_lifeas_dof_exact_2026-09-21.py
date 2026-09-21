#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_lifeas_dof_exact_2026-09-21.py -- 2026-09-21, ceiling-fe agent, FROM ZERO.

TASK B, exact-rational core: how many degrees of freedom does the (Z-dom)
verdict actually consume on the LiFeAs 4-band table?

L1 SOURCE, read first-hand this session (page text layer re-extracted with
   pdftotext -layout -enc UTF-8  ->  _fe_1304_utf8.txt ; the pre-existing
   _umm1304.txt in the main repo has the Greek letters DROPPED, e.g. "nu13"
   prints as bare "13" -- exactly the CLAUDE.md two-instrument failure mode):

G. A. Ummarino et al., "A phenomenological multiband Eliashberg model for
LiFeAs", arXiv:1304.3638.

  Eq.(4)   lam_ij =  [[ l11,  0  , l13      , l14      ],
                      [ 0  ,  0  , l23      , l24      ],
                      [ l13*nu13, l23*nu23, 0, 0       ],
                      [ l14*nu14, l24*nu24, 0, 0       ]]
           "where nu_ij = N_i(0)/N_j(0)"
  nu13 = 0.9019, nu14 = 1.5010, nu23 = 1.0483, nu24 = 1.7447   (band structure)
  Eq.(5)   a2F_ij = C_ij { L(W+W_ij,Y_ij) - L(W-W_ij,Y_ij) },  L = 1/((W+-W)^2+Y^2)
  "In all the calculations we set  W_ij = W^sf_ij = W^sf_0 = 8 meV, and
   Y_ij = Y^sf_ij = W^sf_ij / 2"           => W0 = 8 meV, Y = 4 meV
  "the cut-off energy is w_c = 18 W^sf_0 and the maximum quasiparticle energy
   is w_max = 21 W^sf_0.  We put Gam_ij = Gam^M_ij = 0"
  "the free parameters are reduced to the five coupling constants
   l13, l23, l14, l24 and l11"
  "we fixed the free parameters in order to reproduce the gap values at low T"
  "we have not considered the possibility to have l22, l33, l44 != 0"  (hand-set 0)

  Table 2 (Tc in K, gaps in meV):
    row          l11   ltot  l13   l23   l14   l24   D1   D2   D3   D4   Tc
    Ex.           -     -     -     -     -     -    5.0  2.6  3.6  2.9  18.0
    sf           0.00  1.80  1.78  0.66  0.45  0.52  3.7  2.6  3.6  2.9  15.9
    sf,*         2.10  2.00  1.15  0.80  0.45  0.30  5.0  2.6  3.6  2.9  18.6
    sf,ph(1)     0.86  1.62  1.06  0.79  0.42  0.30  5.1  2.6  3.7  2.9  20.0
    sf,ph(2)     0.90  1.63  1.15  0.80  0.45  0.30  5.0  2.6  3.6  2.9  20.1

  Channel assignment (verbatim): "lam^ph_ii = lam^ph_ij = 0 ... and lam^sf_ii
  = 0, i.e. SF produce only interband coupling"; then "it is necessary to
  introduce at least an intraband coupling in the first band, then lam11 != 0",
  and for rows sf,ph(1)/(2) the (1,1) spectral function is the CALCULATED
  phonon DOS G(W) / the calculated a2F_tot(W) -- i.e. lam11 is the PHONON
  element.  => supports are DISJOINT: ph only at (1,1), sf only interband.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from fractions import Fraction as F
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

KB = 0.08617333262
G2 = 1.318140212132
C_GF = (np.pi ** 2 / 6.0) / G2
C_INF = np.sqrt(G2) / (2 * np.pi)

NU = {"13": F(9019, 10000), "14": F(15010, 10000),
      "23": F(10483, 10000), "24": F(17447, 10000)}

ROWS = [   # tag, l11, l13, l23, l14, l24, Tc_calc
    ("sf",       F(0),          F(178, 100), F(66, 100), F(45, 100), F(52, 100), 15.9),
    ("sf,*",     F(210, 100),   F(115, 100), F(80, 100), F(45, 100), F(30, 100), 18.6),
    ("sf,ph(1)", F(86, 100),    F(106, 100), F(79, 100), F(42, 100), F(30, 100), 20.0),
    ("sf,ph(2)", F(90, 100),    F(115, 100), F(80, 100), F(45, 100), F(30, 100), 20.1),
]
TC_EXP = 18.0


def sep(t):
    print()
    print("=" * 78)
    print(t)
    print("=" * 78)


# ---------------------------------------------------------------- 1. nu rank-1
sep("1. Is the nu matrix rank-one?  (detailed balance nu_ij = N_i/N_j forces")
print("   nu13 * nu24 == nu14 * nu23.)   EXACT RATIONAL ARITHMETIC.")
a = NU["13"] * NU["24"]
b = NU["14"] * NU["23"]
print("   nu13*nu24 = %s = %.8f" % (a, float(a)))
print("   nu14*nu23 = %s = %.8f" % (b, float(b)))
d = abs(a - b) / ((a + b) / 2)
print("   exact relative difference = %s = %.3e" % (d, float(d)))
# printing-precision propagation: each nu is printed to 4 decimals => half-width
# 0.5e-4 absolute; relative half-width of a product of two = sum of the two
hw = sum(F(5, 100000) / NU[k] for k in ("13", "24"))
hw2 = sum(F(5, 100000) / NU[k] for k in ("14", "23"))
tot = hw + hw2
print("   propagated relative half-width of the comparison = %s = %.3e"
      % (tot, float(tot)))
print("   VERDICT: |rel diff| %s propagated half-width  ->  %s"
      % ("<=" if d <= tot else ">", "consistent with rank-one"
         if d <= tot else "INCONSISTENT"))
print("   [rational-exact]  4 rationals (the four printed nu) / 2 multiplications")
print("   / 1 comparison / 4 divisions for the half-width -- zero floating point")
print("   in the decision.  => the four printed nu carry only THREE independent")
print("   numbers (the ratios N1:N2:N3:N4).")
print()
print("   Two independent routes to N4/N3 (must agree if rank-one):")
print("      nu13/nu14 = %s = %.6f" % (NU["13"] / NU["14"], float(NU["13"] / NU["14"])))
print("      nu23/nu24 = %s = %.6f" % (NU["23"] / NU["24"], float(NU["23"] / NU["24"])))


# -------------------------------------------------- 2. the (Z-dom) margin, exact
sep("2. THE (Z-dom) VERDICT ON THIS TABLE -- exact, element by element")
print("""
   (K1) split as the paper defines it:
        Lam^Z_ij  = Lam^ph_ij + Lam^sf_ij     (both elementwise >= 0)
        Lam^Dl_ij = Lam^ph_ij - Lam^sf_ij
   (Z-dom) asks:  |Lam^Dl_ij(k)| <= Lam^Z_ij(k)  for every i,j,k.
   Elementwise margin  m_ij(k) := Lam^Z_ij(k)/|Lam^Dl_ij(k)|.

   The paper's own channel assignment puts
        ph  ONLY at (1,1)          [phonon DOS G(W) or a2F_tot(W)]
        sf  ONLY at the interband positions (1,3),(1,4),(2,3),(2,4) + mirrors
   i.e. the supports are DISJOINT, so at every element exactly one channel is
   present and  m_ij(k) = (x+0)/|x-0| = 1  identically.
""")


def margins_exact(l11, l13, l23, l14, l24):
    """Return dict (i,j) -> exact margin Fraction, or None where Lam^Dl = 0."""
    ph = [[F(0)] * 4 for _ in range(4)]
    sf = [[F(0)] * 4 for _ in range(4)]
    ph[0][0] = l11
    sf[0][2], sf[2][0] = l13, l13 * NU["13"]
    sf[0][3], sf[3][0] = l14, l14 * NU["14"]
    sf[1][2], sf[2][1] = l23, l23 * NU["23"]
    sf[1][3], sf[3][1] = l24, l24 * NU["24"]
    out = {}
    for i in range(4):
        for j in range(4):
            Z = ph[i][j] + sf[i][j]
            D = ph[i][j] - sf[i][j]
            out[(i, j)] = None if D == 0 else Z / abs(D)
    return out


for tag, l11, l13, l23, l14, l24, tcc in ROWS:
    m = margins_exact(l11, l13, l23, l14, l24)
    nz = {k: v for k, v in m.items() if v is not None}
    vals = sorted(set(nz.values()))
    print("   row %-9s : %d elements with Lam^Dl != 0 ; distinct exact margins ="
          " %s" % (tag, len(nz), [str(v) for v in vals]))
print("   => every margin is the EXACT rational 1.  (Z-dom) holds with EQUALITY.")

sep("3. DEGREES OF FREEDOM CONSUMED BY THE VERDICT  (the load-bearing count)")
print("""
   The coordinator's question was: does the criterion eat the 5 FITTED
   couplings, or the 3 first-principles DOS ratios?

   ANSWER: NEITHER.  It eats ZERO of the 8.  Proof, and then a stress test.

   [analytic]  m_ij = (1 + r_ij)/|1 - r_ij| with r_ij := Lam^sf_ij/Lam^ph_ij.
   Disjoint supports => r_ij in {0, +inf} at every element => m_ij = 1
   identically, for ANY values of l11, l13, l23, l14, l24, nu13, nu14, nu23,
   nu24, and for any W0, Y, T.  The verdict is a TAUTOLOGY of the paper's
   MODELLING CHOICE (which channel sits on which matrix element), not a
   measurement of LiFeAs.
""")
print("   EXACT STRESS TEST: replace all 8 published numbers by arbitrary")
print("   rationals (including adversarial ones) and recompute the margins.")
import random
random.seed(20260921)
bad = 0
ntr = 20000
for _ in range(ntr):
    vals = [F(random.randint(1, 10**6), random.randint(1, 10**6))
            for _ in range(5)]
    saveNU = dict(NU)
    for k in NU:
        NU[k] = F(random.randint(1, 10**6), random.randint(1, 10**6))
    m = margins_exact(*vals)
    NU.update(saveNU)
    if any(v is not None and v != 1 for v in m.values()):
        bad += 1
print("   rows with any margin != exactly 1 : %d / %d      => DOF = 0" % (bad, ntr))
print("   [rational-exact] 20000 x 9 random rationals, zero floating point.")
print()
print("   NEGATIVE CONTROL (the test is not vacuous): let the two channels")
print("   OVERLAP by putting a phonon piece p on the (1,3) element as well.")
for p in (F(1, 10), F(1, 2), F(1), F(178, 100)):
    l13 = F(178, 100)
    Z, D = p + l13, p - l13
    print("      ph_13 = %-8s  sf_13 = 1.78  ->  margin_13 = %s"
          % (p, "INFINITE (exact compensation, Lam^Dl_13 = 0)" if D == 0
             else "%s = %.6f" % (Z / abs(D), float(Z / abs(D)))))
print("   => the margin does move once supports overlap; it is pinned at 1")
print("      only because the paper assigned the channels disjointly.")

# ------------------------------------------------------------- 4. moment closed form
sep("4. SECOND MOMENT OF THE LORENTZIAN-DIFFERENCE SPECTRUM -- closed form")


def shape(W, W0, Y):
    return 1.0 / ((W - W0) ** 2 + Y ** 2) - 1.0 / ((W + W0) ** 2 + Y ** 2)


print("""   [analytic] D is ODD (D(-W) = -D(W)), so with u-shifts and
      int_R du/(u^2+Y^2) = pi/Y :
         int_0^inf W D dW = (1/2) int_R W D dW = pi W0 / Y
      and by partial fractions  PV int_R dW/(W((W-a)^2+Y^2)) = a pi/(Y(a^2+Y^2)),
      applied at a = +W0 and a = -W0 :
         int_0^inf D/W dW = pi W0 / (Y (W0^2 + Y^2))
      =>  <W^2> := (int W D)/(int D/W) = W0^2 + Y^2   EXACTLY.""")
for (W0, Y) in ((10.0, 2.0), (8.0, 4.0), (8.0, 1.3), (25.0, 12.5)):
    Ip = quad(lambda W: W * shape(W, W0, Y), 0, np.inf, limit=500)[0]
    Im = quad(lambda W: shape(W, W0, Y) / W, 0, np.inf, limit=500)[0]
    print("      W0=%5.1f Y=%5.1f : int W D = %11.6f (pi W0/Y = %11.6f, rel %.1e)"
          "  <W^2> = %9.4f (W0^2+Y^2 = %9.4f, rel %.1e)"
          % (W0, Y, Ip, np.pi * W0 / Y, abs(Ip - np.pi * W0 / Y) / Ip,
             Ip / Im, W0 ** 2 + Y ** 2,
             abs(Ip / Im - (W0 ** 2 + Y ** 2)) / (W0 ** 2 + Y ** 2)))
print("   PHYSICAL MEANING: peak energy and half-width add in quadrature.")
print("   A broader Lorentzian carries MORE second moment at fixed lambda,")
print("   hence a HIGHER ceiling -- width is not a cosmetic parameter.")

# ------------------------------------------------ 5. symmetrisation, convention-free
sep("5. SYMMETRISATION IS CONVENTION-FREE:  Lam~_ij = sqrt(lam_ij lam_ji)")
print("""
   NOTE a genuine convention clash between two papers of the same group:
     arXiv:0904.1808 prints element (1,3) = lam31 * nu1 with nu1 = N1/N3
        => lam_13/lam_31 = N1/N3   (DOS on the LEFT index)
     arXiv:1304.3638 prints element (3,1) = lam13 * nu13 with nu13 = N1/N3
        => lam_31/lam_13 = N1/N3   (DOS on the RIGHT index -- the standard
           Eliashberg convention, a2F_ij ~ N_j |g_ij|^2)
   These are OPPOSITE.  I do NOT need to resolve which paper mislabels its
   nu, because the symmetrised (detailed-balance) matrix is the GEOMETRIC MEAN

        Lam~_ij = sqrt(lam_ij lam_ji) = sqrt(N_i N_j) V_ij

   which is identical under both readings.  [analytic]  Every number below
   uses only sqrt(lam_ij lam_ji), so no ceiling result depends on the clash.
""")
for tag, l11, l13, l23, l14, l24, tcc in ROWS:
    t13 = float(l13) * float(NU["13"]) ** 0.5
    t14 = float(l14) * float(NU["14"]) ** 0.5
    t23 = float(l23) * float(NU["23"]) ** 0.5
    t24 = float(l24) * float(NU["24"]) ** 0.5
    print("   row %-9s  Lam~_13=%.5f  Lam~_14=%.5f  Lam~_23=%.5f  Lam~_24=%.5f"
          % (tag, t13, t14, t23, t24))


# --------------------------------------------------------- 6. operator + ceiling
def u_vec(N):
    n = np.arange(N)
    return 1.0 / np.sqrt(2.0 * n + 1.0)


def build_O(LamD, LamZ, N):
    B = LamD(0).shape[0]
    u = u_vec(N)
    LD = [LamD(k) for k in range(2 * N + 2)]
    LZ = [LamZ(k) for k in range(N + 2)]
    RZ = np.zeros((N, B))
    for i in range(B):
        acc = LZ[0][i].sum()
        RZ[0, i] = acc
        for n in range(1, N):
            acc += 2 * LZ[n][i].sum()
            RZ[n, i] = acc
    O = np.zeros((N * B, N * B))
    for n in range(N):
        for m in range(N):
            blk = LD[n + m + 1].copy()
            if n != m:
                blk = blk + LD[abs(n - m)]
            blk = blk * (u[n] * u[m])
            if n == m:
                blk = blk + u[n] ** 2 * LD[0] - np.diag(u[n] ** 2 * RZ[n])
            O[n * B:(n + 1) * B, m * B:(m + 1) * B] = blk
    return 0.5 * (O + O.T)


def sf_matrix(l13, l23, l14, l24):
    M = np.zeros((4, 4))
    M[0, 2] = M[2, 0] = float(l13) * float(NU["13"]) ** 0.5
    M[0, 3] = M[3, 0] = float(l14) * float(NU["14"]) ** 0.5
    M[1, 2] = M[2, 1] = float(l23) * float(NU["23"]) ** 0.5
    M[1, 3] = M[3, 1] = float(l24) * float(NU["24"]) ** 0.5
    return M


sep("6. GATE: solve the PURE-SF row from my own operator and compare to the")
print("   paper's own Tc_calc = 15.9 K.  This row (lam11 = 0) is COMPLETELY")
print("   specified by the paper -- one spectrum (W0=8, Y=4), one matrix -- so")
print("   it is a true end-to-end check of inputs + convention + operator.")
W0sf, Ysf = 8.0, 4.0
Im1 = quad(lambda W: shape(W, W0sf, Ysf) / W, 0, np.inf, limit=500)[0]


def Jfac(k, T, W0=W0sf, Y=Ysf):
    a = (2 * np.pi * KB * T * k) ** 2
    return quad(lambda W: W * shape(W, W0, Y) / (W ** 2 + a), 0, np.inf,
                limit=500)[0] / Im1


tag, l11, l13, l23, l14, l24, tcc = ROWS[0]
SF = sf_matrix(l13, l23, l14, l24)


def f(T):
    cache = {}

    def J(k):
        if k not in cache:
            cache[k] = Jfac(k, T)
        return cache[k]
    return np.linalg.eigvalsh(build_O(lambda k: -SF * J(k),
                                      lambda k: SF * J(k), 160))[-1] - 1.0


Tc_mine = brentq(f, 1.0, 200.0, xtol=1e-3)
print("   paper Tc_calc = 15.9 K   |   my Tc (no cutoff) = %.2f K   (diff %+.1f%%)"
      % (Tc_mine, 100 * (Tc_mine - 15.9) / 15.9))
print("   DIRECTION: the paper solves WITH theta(w_c-|w_m|), w_c = 18 W0 = 144")
print("   meV, which removes pairing weight => lowers Tc.  My no-cutoff value")
print("   is therefore expected to sit slightly ABOVE the paper's.")

sep("7. CEILING TEST FOR LiFeAs   Tc <= C_inf sqrt(rho(M2''_Z))")
print("""
   DIRECTION AUDIT.  I do NOT have the calculated phonon spectrum of Fig.2
   digitised, so <W^2>_ph for the (1,1) element is UNKNOWN.  rho of a
   nonnegative matrix is monotone increasing in every entry (Perron-Frobenius),
   so setting M2^ph_11 := 0 gives a STRICT LOWER BOUND on the ceiling.
   A lower bound on an upper bound supports only the verdict "consistent";
   it can never be used to claim a violation.  Every row below is therefore a
   SOUND consistency statement and NOT a sound falsification test.
""")
W2sf = W0sf ** 2 + Ysf ** 2
print("   <W^2>_sf = W0^2 + Y^2 = %.1f meV^2   (W0 = 8, Y = 4 meV)" % W2sf)
print("   %-9s %9s %9s %9s %10s %10s %9s %9s"
      % ("row", "lam_pap", "lam_rowmx", "rho(lam~)", "ceil_lo(K)",
         "Tc_exp/ceil", "Tc_calc", "calc/ceil"))
LAM_PAPER = {"sf": 1.80, "sf,*": 2.00, "sf,ph(1)": 1.62, "sf,ph(2)": 1.63}
for tag, l11, l13, l23, l14, l24, tcc in ROWS:
    SF = sf_matrix(l13, l23, l14, l24)
    M2 = SF * W2sf                     # phonon (1,1) DROPPED -> lower bound
    Mpp = C_GF * M2                    # purely off-diagonal => all get c_GF
    rho = np.abs(np.linalg.eigvals(Mpp)).max()
    ceil = C_INF * np.sqrt(rho) / KB
    lam_rowmx = float(l11) + SF.sum(axis=1).max()
    print("   %-9s %9.2f %9.3f %9.4f %10.2f %10.4f %9.1f %9.4f"
          % (tag, LAM_PAPER[tag], lam_rowmx, np.abs(np.linalg.eigvals(SF)).max(),
             ceil, TC_EXP / ceil, tcc, tcc / ceil))
print()
print("   ! 'lam_pap' is the paper's own lam_tot column; 'lam_rowmx' is MY")
print("     explicitly-defined quantity lam11 + max_i sum_j Lam~_ij.  They")
print("     differ (e.g. 1.80 vs 2.366) because the paper does not state the")
print("     averaging convention behind its lam_tot.  I do NOT reconcile them;")
print("     I report both and use only my own explicitly-defined column.")
print("   lam is reported next to every ratio, as required: these rows sit at")
print("   lam ~ 1.6-4.0, i.e. genuinely strong coupling, NOT the vacuous")
print("   lam <= 0.5 regime where the repo measured Tc/ceiling ~ 0.003-0.022.")
