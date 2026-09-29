#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_cac6_fese_noise_and_dof_2026-09-19.py -- 2026-09-19, spinfluc-quant worktree,
written from zero.

Question this script answers
----------------------------
For the two FIRST-PRINCIPLES material-side failures of the convention-independent
multiband criterion

        lam_ii lam_jj  >=  lam_ij lam_ji

(CaC6, Sanna PRB 85 184514 Table I; FeSe, Essenberger arXiv:1411.2121 Fig.3a),
is the failure

  (a) a real property of the computed coupling matrix, or
  (b) an artifact of how the published table was produced
      -- the MgB2 mu*_ij precedent, where four printed numbers carried only
         three degrees of freedom and the criterion value was a CONSTANT
         (2.23 x 2.48) independent of everything the paper actually fitted?

Three independent things are computed here:

  SECTION 1  Completeness of the criterion as an invariant.  [analytic, checked]
  SECTION 2  CaC6: rigorous rounding interval on the criterion ratio,
             in EXACT RATIONAL ARITHMETIC, plus the Metropolis-noise estimate.
  SECTION 3  FeSe: rounding interval + the degrees-of-freedom count.

ALL INPUT NUMBERS RE-EXTRACTED BY ME (2026-09-19) from the PDF text layer /
page image, NOT copied from the existing repo scripts or notes.
  CaC6 : pdftotext -layout -enc UTF-8 -f 3 -l 3 _sanna_prb85_184514.pdf
  FeSe : pdftoppm -png -r 200 -f 3 -l 3 essenberger1411.2121.pdf  (page image)
My CaC6 extraction agrees entry-for-entry with _chk_cac6_lambda_ij_criterion.py;
my FeSe extraction agrees with the four labels in
ceiling/lambda_ij_criterion_fails_on_FeSe_2026-09-18.md.

DIRECTION AUDIT (CLAUDE.md one-sided-bound rule) -- read before using any number
-------------------------------------------------------------------------------
The criterion verdict here is FAIL (a violation).  A FAIL may NOT be supported by
a one-sided estimate that happens to push the ratio up.  So:

  * To certify a violation as REAL against printing roundoff I must use
    r_min = the SMALLEST ratio consistent with the printed digits.
    r_min > 1  ==>  violation survives  (this is the sound direction for FAIL).
  * r_max is reported only for context; NO verdict is ever read off r_max.
  * The Metropolis-noise column is an ESTIMATE, not a bound.  It can only
    WEAKEN a violation claim ("not resolved"), never strengthen one.
    Accordingly it is used in exactly one direction: to move pairs OUT of the
    "real" column, never into it.
"""

from fractions import Fraction as F
import itertools
import math

# =====================================================================
# INPUT DATA -- my own extraction
# =====================================================================

# --- CaC6, Sanna PRB 85, 184514 (2012), Table I, p.184514-3 ----------
# Reading A (rows = first index), the reading that makes the table
# self-consistent given Eq.(11).  The criterion is transpose-invariant,
# so no verdict below depends on this choice.
T3 = [["0.301", "0.136", "0.257"],
      ["0.546", "0.239", "0.479"],
      ["0.427", "0.198", "0.367"]]
N3 = ["0.412", "0.104", "0.249"]
L3 = ["1", "2", "3"]

T6 = [["0.163", "0.126", "0.099", "0.033", "0.201", "0.046"],
      ["0.179", "0.140", "0.105", "0.035", "0.221", "0.050"],
      ["0.331", "0.245", "0.151", "0.084", "0.384", "0.096"],
      ["0.271", "0.202", "0.206", "0.047", "0.400", "0.080"],
      ["0.252", "0.194", "0.145", "0.061", "0.309", "0.073"],
      ["0.206", "0.157", "0.128", "0.044", "0.259", "0.060"]]
N6 = ["0.250", "0.176", "0.075", "0.031", "0.200", "0.056"]
L6 = ["1a", "1b", "2a", "2b", "3a", "3b"]

# --- FeSe, Essenberger arXiv:1411.2121v1, Fig. 3(a) typeset labels ---
# lam_{++}=0.028 (purple, positive peak), lam_{--}=0.029 (green, positive peak),
# lam_{+-}=0.042 (red, NEGATIVE trough), lam_{-+}=0.048 (black, NEGATIVE trough)
FESE = [["0.028", "0.042"],
        ["0.048", "0.029"]]
FESE_L = ["+", "-"]


def frac(s):
    """Exact rational value of a printed decimal string."""
    return F(s)


def halfwidth(s):
    """Half-width of the rounding interval of a printed decimal string:
    '0.301' -> 1/2000.  Exact."""
    ndec = len(s.split(".")[1])
    return F(1, 2 * 10 ** ndec)


# =====================================================================
# SECTION 1 -- the criterion is a COMPLETE invariant  [analytic]
# =====================================================================
def section1():
    print("=" * 78)
    print("SECTION 1  Is the criterion blind to a degree of freedom?  [analytic]")
    print("=" * 78)
    print("""
  Claim (mine, 2026-09-19).  Let V be symmetric d x d with strictly positive
  diagonal, and let the 'convention group' act by V -> D V D with D = diag(d_i),
  d_i > 0  (this is the group that the repo's convention-independence theorem
  quotients out: lam_ij = f(i) g(j) V_ij).  Then the normalisation

        d_i := 1 / sqrt(V_ii)          (uniquely determined)

  sends V to the correlation matrix C with C_ii = 1 and

        C_ij = V_ij / sqrt(V_ii V_jj).

  Hence C is a COMPLETE invariant, and since the criterion ratio is

        r_ij = V_ij V_ji / (V_ii V_jj) = C_ij^2,

  the set { r_ij }_{i<j} together with the signs of V_ij determines C, i.e.
  determines the entire convention-invariant content of the table.

  CONSEQUENCE, and this is the point:
    There is NO residual degree of freedom for the criterion to be blind to.
    The MgB2 mu* pathology was NOT 'the criterion missed a d.o.f.' -- it was
    that the single invariant the criterion sees had been IMPORTED from
    first-principles refs [15]/[16] rather than produced by that paper's fit.
    So the artifact question for CaC6/FeSe is not 'what is hidden from the
    criterion' but 'is the invariant computed, or is it fixed by an ansatz?'
    -> that is a question about the GENERATIVE PROCESS, settled by reading
       the Computational Details, not by counting printed numbers.

  Numerical check that r_ij = C_ij^2 and that both are D-invariant:
""")
    import random
    random.seed(20260919)
    worst = 0.0
    for trial in range(2000):
        d = random.randrange(2, 7)
        V = [[0.0] * d for _ in range(d)]
        for i in range(d):
            for j in range(i, d):
                V[i][j] = V[j][i] = random.uniform(-1.5, 2.5)
            V[i][i] = abs(V[i][i]) + 0.3          # positive diagonal
        D = [random.uniform(0.05, 20.0) for _ in range(d)]
        W = [[D[i] * V[i][j] * D[j] for j in range(d)] for i in range(d)]
        for i, j in itertools.combinations(range(d), 2):
            rV = V[i][j] * V[j][i] / (V[i][i] * V[j][j])
            rW = W[i][j] * W[j][i] / (W[i][i] * W[j][j])
            C = V[i][j] / math.sqrt(V[i][i] * V[j][j])
            worst = max(worst, abs(rV - rW), abs(rV - C * C))
    print("    worst |r(V) - r(DVD)| and |r - C^2| over 2000 random cases: %.3e" % worst)
    print("    -> criterion ratio IS the squared correlation entry, and is exactly")
    print("       D-invariant.  Complete-invariant claim checked.\n")


# =====================================================================
# SECTION 2 -- CaC6
# =====================================================================
def detailed_balance(T, N, lab, tag):
    """N_i lam_ij = N_j lam_ji is an EXACT consequence of Sanna Eq.(11)
    (given <|g|^2> symmetric).  Its residual is therefore a direct empirical
    error bar on the table: Metropolis sampling noise + printing roundoff.
    Nothing in the physics can make it nonzero."""
    print("  [detailed balance]  %s  -- exact identity, so residual = noise" % tag)
    worst, wp = 0.0, None
    for i, j in itertools.combinations(range(len(N)), 2):
        a = float(N[i]) * float(T[i][j])
        b = float(N[j]) * float(T[j][i])
        rel = abs(a - b) / max(a, b)
        if rel > worst:
            worst, wp = rel, (lab[i], lab[j])
    print("      worst relative residual = %.2f%%  at %s-%s" % (100 * worst, wp[0], wp[1]))
    return worst


def criterion_exact(T, N, lab, tag, dbnoise):
    """Exact-rational criterion ratio + RIGOROUS rounding interval.

    DIRECTION: the verdict of interest is FAIL (r > 1).  The sound quantity for
    a FAIL is r_min, the smallest ratio compatible with the printed digits:
        r_min = (lam_ij - h)(lam_ji - h) / ((lam_ii + h)(lam_jj + h)).
    r_min > 1  ==>  the violation cannot be an artifact of printing roundoff.
    r_max is context only; no verdict is read from it."""
    print("\n  [criterion, exact rational + rounding interval]  %s" % tag)
    print("      %-9s %-11s %-11s %-11s %s" % ("pair", "r (exact)", "r_min", "r_max", "verdict"))
    n = len(N)
    rows = []
    nrat = 0
    for i, j in itertools.combinations(range(n), 2):
        lii, ljj = frac(T[i][i]), frac(T[j][j])
        lij, lji = frac(T[i][j]), frac(T[j][i])
        hii, hjj = halfwidth(T[i][i]), halfwidth(T[j][j])
        hij, hji = halfwidth(T[i][j]), halfwidth(T[j][i])
        nrat += 8
        r = (lij * lji) / (lii * ljj)
        r_min = ((lij - hij) * (lji - hji)) / ((lii + hii) * (ljj + hjj))
        r_max = ((lij + hij) * (lji + hji)) / ((lii - hii) * (ljj - hjj))
        # Metropolis-noise estimate: 4 entries in quadrature, each at the
        # detailed-balance residual level.  ESTIMATE, NOT A BOUND.
        sig = 2.0 * dbnoise            # sqrt(4) * per-entry relative noise
        if float(r) <= 1.0:
            verdict = "ok"
        elif float(r_min) <= 1.0:
            verdict = "FAIL, but roundoff alone can erase it"
        elif (float(r) - 1.0) < sig:
            verdict = "FAIL, but < 1 sigma of MC noise -- UNRESOLVED"
        elif (float(r) - 1.0) < 3.0 * sig:
            verdict = "FAIL, 1-3 sigma -- weak"
        else:
            verdict = "FAIL, > 3 sigma -- REAL"
        rows.append((float(r), lab[i], lab[j], verdict))
        print("      %-3s-%-3s   %-11.5f %-11.5f %-11.5f %s"
              % (lab[i], lab[j], float(r), float(r_min), float(r_max), verdict))
    print("      [rational-exact: %d pairs / %d rationals in this block]"
          % (n * (n - 1) // 2, nrat))
    print("      [MC-noise sigma used: %.2f%% per entry -> %.2f%% on the ratio]"
          % (100 * dbnoise, 100 * sig))
    return rows


def section2():
    print("=" * 78)
    print("SECTION 2  CaC6 -- is the violation resolved by the data?")
    print("=" * 78)
    print("""
  GENERATIVE PROCESS  [L1, my own pdftotext -enc UTF-8 of p.184514-3]:

    "The double Brillouin zone integration appearing in the definition of the
     band-resolved Eliashberg functions in Eq. (11) is evaluated with a
     METROPOLIS INTEGRATION SCHEME.  Random k points are generated on the
     Brillouin zone, then each k point is accepted or rejected with a
     probability depending on xi_k ... We use a set of about 2 x 10^4 ACCEPTED
     K POINTS PER BAND.  Then, electron-phonon matrix elements on this random
     mesh are obtained via INTERPOLATION from those calculated on the regular
     grid."
    (ESPRESSO, GGA-PW, ultrasoft PP, 30/300 Ry, 6x6x6 k-grid, e-ph on 10x10x10.)

  => Each lam_JJ' is an INDEPENDENT stochastic double-BZ integral over genuinely
     k-resolved DFPT matrix elements.  There is NO fitted scale factor and NO
     DOS-ratio reduction anywhere in the chain.  Structurally unlike MgB2 mu*.
     ==> artifact hypothesis (b) is REFUTED for CaC6.

  BUT the same sentence says the entries carry MONTE CARLO NOISE.  That is a
  new, quantifiable caveat, and it is what this section measures.
""")
    for T, N, lab, tag in ((T3, N3, L3, "three-band"), (T6, N6, L6, "six-band")):
        db = detailed_balance(T, N, lab, tag)
        criterion_exact(T, N, lab, tag, db)
        print()


# =====================================================================
# SECTION 3 -- FeSe
# =====================================================================
def section3():
    print("=" * 78)
    print("SECTION 3  FeSe -- degrees of freedom and roundoff")
    print("=" * 78)
    print("""
  GENERATIVE PROCESS  [L1, my own text-layer extraction]:
    p.1: "As neither the functional, nor the theoretical framework, contain any
          ADJUSTABLE PARAMETERS, this scheme promises the prediction of Tc ..."
    p.2: "we account for: 1) the phononic interaction by a band-resolved
          Eliashberg function alpha^2 F_ij (detailed definitions can be found in
          Ref. 45 ...)"
    Sanna is a co-author of BOTH papers; the band-resolved alpha^2 F_ij is the
    same Eq.(11)-type object as CaC6's.
    The single free parameter in the paper (alpha, the xc scaling) enters the
    SPIN-FLUCTUATION susceptibility only -- Fig.3c shows Lambda^SF vs alpha/alpha_c
    -- NOT the phonon alpha^2 F_ij of Fig.3a.
  ==> the four phonon lam_ij are parameter-free first-principles output.
      Artifact hypothesis (b) REFUTED for FeSe too.

  DEGREES OF FREEDOM (the MgB2 discipline, applied honestly):
    4 printed numbers.  Convention group f(i)g(j) is 3-dimensional (4 params,
    1 redundant overall scale).  4 - 3 = 1 invariant, and by SECTION 1 that
    invariant IS the criterion ratio.  So the criterion consumes 100% of the
    convention-invariant content -- there is no spare d.o.f. it could be
    blind to, and no printed number is inert.
    Contrast MgB2 mu*: also 2x2, also 1 invariant -- but there the invariant
    was p11 p22/(p12 p21) = 2.23 x 2.48, IMPORTED from refs [15]/[16].  Here
    nothing is imported.
""")
    dbnoise_fese = None
    # detailed balance cannot be run: the paper prints no N_+ , N_- .
    print("  [detailed balance]  NOT AVAILABLE -- the paper prints no DOS for the")
    print("      two FS portions, so the exact identity N_i lam_ij = N_j lam_ji")
    print("      cannot be used as an error bar.  FeSe has NO independent noise")
    print("      estimate.  Marked UNTESTED, not 'passed'.\n")
    print("  [criterion, exact rational + rounding interval]")
    lii, ljj = frac(FESE[0][0]), frac(FESE[1][1])
    lij, lji = frac(FESE[0][1]), frac(FESE[1][0])
    h = F(1, 2000)
    r = (lij * lji) / (lii * ljj)
    r_min = ((lij - h) * (lji - h)) / ((lii + h) * (ljj + h))
    r_max = ((lij + h) * (lji + h)) / ((lii - h) * (ljj - h))
    print("      r       = %s = %.5f" % (r, float(r)))
    print("      r_min   = %s = %.5f   <-- sound direction for a FAIL" % (r_min, float(r_min)))
    print("      r_max   = %s = %.5f   (context only)" % (r_max, float(r_max)))
    print("      [rational-exact: 1 pair / 8 rationals]")
    if float(r_min) > 1.0:
        print("      VERDICT: violation survives printing roundoff by %.2fx" % float(r_min))
    else:
        print("      VERDICT: roundoff alone can erase the violation")
    print("""
  Author's own text-layer corroboration [L1, footnote 70, my extraction]:
    "The fact that phonons reduce the Tc in this material is a direct
     consequence of the stronger interband than the intraband coupling lam_ij
     (~0.03 intraband and ~0.05 interband), effectively acting against the
     s+- ordering."
  Using the author's own round numbers 0.03 / 0.05:""")
    r2 = (F("0.05") * F("0.05")) / (F("0.03") * F("0.03"))
    print("      r = (0.05)^2/(0.03)^2 = %s = %.4f  -> still FAIL" % (r2, float(r2)))
    print("      => the verdict does not depend on reading 3 significant figures")
    print("         off a figure; it survives on the author's prose alone.")


# =====================================================================
# SECTION 4 -- what the sign structure does and does not touch
# =====================================================================
def section4():
    print()
    print("=" * 78)
    print("SECTION 4  Sign structure: which repo condition is actually at risk")
    print("=" * 78)
    print("""
  [L1, my own reading of the 200-dpi page image of arXiv:1411.2121v1 p.3]
  Fig. 3(a), y-axis labelled alpha^2 F(omega):
      lam_{--}=0.029 (green dashed) and lam_{++}=0.028 (purple): POSITIVE peaks
      lam_{+-}=0.042 (red) and lam_{-+}=0.048 (black): NEGATIVE troughs, min ~ -0.2
  This independently closes the premise that the 2026-09-18 note flagged as
  "read from the same figure, not independently verified".

  [analytic]  A bare Eq.(11) object
      alpha^2 F_ij = (1/N_i) SUM_{k in i, k' in j} delta delta SUM_nu
                     delta(omega - omega_{q,nu}) |g^nu_{k,k'}|^2
  is a sum of SQUARED MODULI and is therefore >= 0 POINTWISE.
  ==> the negative off-diagonal traces in Fig. 3(a) CANNOT be the bare Eq.(11)
      alpha^2 F_ij.  Some sign convention or kernel structure is folded into the
      plotted object.  WHAT that convention is, I do NOT determine.  OPEN.

  Why this matters, and why it does NOT reopen the ceiling chain:

  (i) For the CRITERION the signs are irrelevant -- both off-diagonal entries
      carry the same sign, so the product lam_ij lam_ji is unchanged.  The
      2026-09-18 verdict (2.4828x) stands either way.  Already established there.

  (ii) For the PHONON-SIDE condition the signs are NOT irrelevant, because
      (ML3) consumes ELEMENTWISE NONNEGATIVITY of Lambda-tilde, and the signed
      form consumes (Z-dom)  |Lambda^Delta_ij| <= Lambda^Z_ij.

  (iii) [analytic]  But for a PURE PHONON channel (Z-dom) holds with EQUALITY.
      In the standard multiband Eliashberg structure the mass-renormalisation
      (Z) kernel and the pairing (Delta) kernel are built from the SAME
      alpha^2 F_ij: the Z-channel sum and the Delta-channel sum differ by which
      gap function multiplies them, not by the kernel.  Hence
          |Lambda^{ph,Delta}_ij| = Lambda^{ph,Z}_ij   for every i,j,
      and (Z-dom) is SATURATED, never violated, by phonons alone.
      ==> FeSe's phonon sign structure BY ITSELF cannot break (Z-dom).

  (iv) What could break (Z-dom) is the SF contribution, where the two channels
      need not share a kernel.  That needs the Lambda^ph and Lambda^SF tables,
      which are still not in hand (the paper prints only lam^max for SF).
      ==> UNCHANGED from 2026-09-18.  This section does not close it and does
          not claim to.
""")


def section5():
    print()
    print("=" * 78)
    print("SECTION 5  Sign GAUGE: why FeSe's negative entries cannot break (Z-dom)")
    print("=" * 78)
    print("""
  [L1, my own 9x-zoom render of the Fig.3c inset, arXiv:1411.2121v1 p.3]
      Lambda^SF_{+-} (red)  and Lambda^SF_{-+} (black dashed): POSITIVE
      Lambda^SF_{++} (purple) and Lambda^SF_{--} (green dashed): NEGATIVE (to ~ -0.03)
  Combined with Fig.3a (phonon: diagonal POSITIVE, off-diagonal NEGATIVE), the
  two plotted matrices have EXACTLY COMPLEMENTARY sign patterns:

      plotted ph :  [[+, -], [-, +]]        plotted sf :  [[-, +], [+, -]]

  Both patterns are of the form  sign = s_i s_j  with s = (+1, -1).
  That is the s+- gauge, and it is the ONLY nontrivial symmetric sign pattern a
  2x2 admits of that form.

  [analytic]  (Z-dom) IS INVARIANT UNDER ANY DIAGONAL SIGN GAUGE.
    Let S = diag(s_i), s_i in {+1,-1}, and suppose the decomposition the ceiling
    consumes is, up to that gauge,
        Lambda^Delta = S (Lambda^ph - Lambda^sf) S ,
        Lambda^Z     =    Lambda^ph + Lambda^sf ,
    with Lambda^ph, Lambda^sf ELEMENTWISE NONNEGATIVE.  Then for every i,j
        |Lambda^Delta_ij| = |s_i s_j| * |Lambda^ph_ij - Lambda^sf_ij|
                          = |Lambda^ph_ij - Lambda^sf_ij|
                          <= Lambda^ph_ij + Lambda^sf_ij  =  Lambda^Z_ij .
    The gauge cancels because |s_i s_j| = 1, and the last step is the (K1)
    triangle inequality.  So (Z-dom) holds REGARDLESS of the gauge.

  CONSEQUENCE for the task question:
    (Z-dom) is stated with ABSOLUTE VALUES.  The signs I read off Figs. 3a/3c are
    therefore IRRELEVANT to it.  FeSe's negative off-diagonal phonon entries do
    NOT threaten (Z-dom), and this holds without my having to resolve what the
    plotted sign convention actually means.
    What (Z-dom) does still need is that the MAGNITUDES |alpha^2 F_ij| and
    |Lambda^SF_ij| are the two (K1) pieces.  That is NOT established here --
    it needs the kernel definitions in Essenberger Refs. 45 / 51.  OPEN.

  Numerical check of the gauge-invariance step:""")
    import random
    random.seed(20260919)
    worst = -1.0
    bad = 0
    for _ in range(20000):
        d = random.randrange(2, 6)
        ph = [[0.0] * d for _ in range(d)]
        sf = [[0.0] * d for _ in range(d)]
        for i in range(d):
            for j in range(i, d):
                ph[i][j] = ph[j][i] = random.uniform(0, 2)      # elementwise >= 0
                sf[i][j] = sf[j][i] = random.uniform(0, 2)      # elementwise >= 0
        s = [random.choice((1, -1)) for _ in range(d)]
        for i in range(d):
            for j in range(d):
                lhs = abs(s[i] * s[j] * (ph[i][j] - sf[i][j]))
                rhs = ph[i][j] + sf[i][j]
                worst = max(worst, lhs - rhs)
                if lhs > rhs + 1e-12:
                    bad += 1
    print("      20000 random cases, d=2..5, random sign gauge:")
    print("      violations of |Lambda^Delta_ij| <= Lambda^Z_ij : %d" % bad)
    print("      worst (lhs - rhs) = %.3e   (<= 0 means (Z-dom) holds)" % worst)
    print("""
  And the CONTRAST that shows the check is not vacuous: if the two pieces are
  NOT elementwise nonnegative (i.e. (K1)'s premise genuinely fails), (Z-dom)
  does break.  Same driver, pieces allowed to go negative:""")
    random.seed(20260919)
    bad2 = 0
    tot = 0
    for _ in range(20000):
        d = random.randrange(2, 6)
        ph = [[0.0] * d for _ in range(d)]
        sf = [[0.0] * d for _ in range(d)]
        for i in range(d):
            for j in range(i, d):
                ph[i][j] = ph[j][i] = random.uniform(-2, 2)     # signed
                sf[i][j] = sf[j][i] = random.uniform(-2, 2)     # signed
        s = [random.choice((1, -1)) for _ in range(d)]
        for i in range(d):
            for j in range(d):
                tot += 1
                lhs = abs(s[i] * s[j] * (ph[i][j] - sf[i][j]))
                rhs = ph[i][j] + sf[i][j]
                if lhs > rhs + 1e-12:
                    bad2 += 1
    print("      violations: %d / %d entries (%.1f%%)  -> the premise is load-bearing"
          % (bad2, tot, 100.0 * bad2 / tot))


if __name__ == "__main__":
    section1()
    section2()
    section3()
    section4()
    section5()
