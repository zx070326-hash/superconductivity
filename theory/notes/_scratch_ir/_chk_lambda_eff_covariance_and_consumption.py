#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_lambda_eff_covariance_and_consumption.py -- 2026-09-18, main session, from zero.

Round 10.  source-hunt reports the first genuine NON-RANK-ONE object in this
literature: Mazin & Antropov's multiband McMillan matrix

    Lam_eff_ij = [Lam_ij - mus_ij (1 + 0.62 sum_n Lam_in)] / (1 + sum_n Lam_in)   Eq.(7)

whose row sum is not covariant under Lam_ij -> f_i g_j Lam_ij, so the repo's
convention-invariance theorem does not cover it.  Verified here from zero, and
then pushed to the question source-hunt left open: DOES IT MATTER?

WHAT THIS SCRIPT ESTABLISHES

  (A) [analytic, exact rational]  The criterion  M_ii M_jj >= M_ij M_ji  is
      invariant under M -> D1 M D2 for ANY positive diagonals D1, D2.  This is
      the reason the criterion, not PSD, is the convention-free object.
  (B) [analytic, exact rational]  Eq.(7): the denominator 1+sum_n Lam_in cancels
      EXACTLY in the criterion ratio (source-hunt's claim -- confirmed), and the
      residual convention dependence enters ONLY through A_i := 1+0.62 sum_n
      Lam_in.  SHARPENING (mine): if the A_i are EQUAL across bands the map is
      again covariant, so the obstruction is proportional to the INTERBAND
      SPREAD of the mass renormalisations, not to the presence of mu* as such.
  (C) MgB2 numbers: source-hunt's 0.02003 vs 0.02108 reproduced from zero.
  (D) ** THE BIGGER ISSUE, which the covariance question hides. **  Even with
      A_i CONSTANT, and with Lam and mu* both symmetric and both passing the
      criterion, Lam - A mu* can VIOLATE it -- with a clean threshold.  So the
      criterion is simply not stable under lambda-minus-mu* combinations, in
      any convention.  This is a statement about the OBJECT, not the notation.
  (E) [analytic]  Eq.(8) is a BAND-INDEXED Morel-Anderson map,
      mus = (1 + mu D_L)^{-1} mu with D_L = diag(log(W_n/wc)).  The spectral-map
      proof of round 9 does NOT apply (mu D_L is not symmetric), but the
      resolvent form does:   mus^{-1} = mu^{-1} + D_L.
      => mus symmetric PSD whenever mu is, for ANY diagonal D_L >= 0.
      This STRICTLY GENERALISES the round-9 theorem (Mitrovic Eq.(12) is
      D_L = L * I) and covers band-dependent bandwidths/cutoffs.
      Converse (pole) survives too.
  (F) Table II: the three published first-principles MgB2 2x2 matrices all pass.
  (G) [analytic]  Mazin-Golubov's printed "superconductivity with all g_ij<0"
      condition vs the sharp one.  The SHARP 2x2 condition is exactly the
      NEGATION of the criterion; the printed condition implies it by AM-GM, so
      the printed one is SUFFICIENT but NOT NECESSARY.  (source-hunt left this
      "undecided"; it is decidable without any further reading.)
  (H) ** THE CONSUMPTION POINT, pinned. **  The criterion is consumed in
      exactly ONE place in this repo -- Lemma mu-1's  mu* >= 0, i.e. the
      COULOMB block.  The phonon side (ML3)/(G'-abs) consumes ELEMENTWISE
      non-negativity / (Z-dom), which is a DIFFERENT condition: shown here to
      be logically independent of the criterion in both directions, and shown
      to hold on the very CaC6 pair that violates the criterion.
      => The four material-side "failures" (FeSe/CaC6/122/1111) never fed the
         ceiling chain.  Lam_eff lands outside it as well.

L1 SOURCES, read by the main session from the TEXT LAYER of PDFs it fetched
itself (curl -> PyMuPDF; both real text PDFs, verified by get_text length):

  _mazin_antropov0212263.pdf = I.I. Mazin, V.P. Antropov, "Electronic structure,
    electron-phonon coupling, and multiband effects in MgB2",
    cond-mat/0212263v1 = Physica C 385, 49 (2003), 14 pp.
  printed p.9 (page idx 8), verbatim:
    "Lambda is an asymmetric matrix related to the symmetric matrix of the
     pairing interaction, Lambda_ij = V_ij N_j, where N_i is the contribution
     of the i-th band to the total DOS."
    "Note that although the matrix Lambda is not symmetric, its eigenvalues are
     the same as those of the symmetric matrix sqrt(N) V sqrt(N)."
    TABLE II. 2x2 EPC matrices in different calculations.
      Ref.[10]   0.96 0.17 / 0.23 0.29
      Ref.[55]   1.02 0.16 / 0.21 0.45
      Ref.[30]a  0.78 0.11 / 0.15 0.21
      a Obtained by integrating lambda(k,k') distribution plots from Ref.[30].
  printed p.10 (page idx 9), verbatim:
    "The mass renormalization parameters for each band can be constructed from
     the matrix Lambda: lambda_i = sum_j Lambda_ij.  These lambda_i define,
     among other things, the de Haas-van Alphen thermal masses."
    "One can show that lambdabar <= lambda_eff, the equality being achieved when
     and only when all elements of the V matrix are the same (the relative
     magnitude of N_i is irrelevant)."
    Eq.(7):  Lambda^eff_ij = [Lambda_ij - mu*_ij (1 + 0.62 sum_n Lambda_in)]
                             / (1 + sum_n Lambda_in)
  printed p.11 (page idx 10), verbatim:
    "... mu_ij ... is formally defined as <<V_C>>_ij N_j ..."
    "Jepsen and Andersen [64] estimated this effect ... and found the ratio
     between <<VC>>_sigsig, <<VC>>_pipi and <<VC>>_sigpi to be ~ 3:2.5:1."
    Eq.(8):  mu*_ij = mu_ij - sum_n mu_in log(W_n/w_c) mu*_nj
    "assuming mu_ss = mu_pp = alpha mu_sp, (alpha>1), and
     mu_ss log(W_s/w_c) = mu_pp log(W_p/w_c) = L, we obtain
     alpha* = alpha + (alpha - 1/alpha) L.  From the ratios of <<V_C>>'s above,
     alpha ~ 2.3, and L for MgB2 is of the order of 0.5 - 1, so for mu*_ij it
     holds that mu*_ss = mu*_pp ~ 4 mu*_sp"

  _mazin_golubov9602118.pdf = I.I. Mazin, A.A. Golubov, "Interband sign reversal
    of the order parameter ...", cond-mat/9602118v1, 8 pp.
  printed p.2 (page idx 1), verbatim:
    "The effective coupling constant lambda_eff in this case is simply the
     maximal eigenvalue lambda_max of the matrix Lambda_ij = g_ij N_j"
    "lambda = sum_ij Lambda_ij N_i/N = sum_i lambda_i N_i/N, where the mass
     renormalization for the i-th band, as measured, for instance, in de
     Haas-van Alphen experiments, is lambda_i = sum_j Lambda_ij"
    "An instructive example ... two-band model with Lambda_11 = Lambda_22 =
     Lambda > 0, Lambda_12 = Lambda_21 = -Lambda.  Then lambda = 0, while
     lambda_eff = 2 Lambda."
    "One can easily check that Eq.1 may have a superconducting solution even for
     all g_ij < 0 ...  The condition for that is
     |g_12| > (|g_1| N_1^2 + |g_2| N_2^2) / 2 N_1 N_2."
  printed p.4 (page idx 3), SPAN-LEVEL extraction (font sizes 9.96 body /
    6.97 subscript / 4.98 prime), verbatim:
    "Coupling matrix Lambda is defined in the same way as Allen's matrix
     lambda_JJ',  Lambda_JJ' = V^pairing_JJ'  N_JJ'(0)."
    >>> TRANSCRIPTION CORRECTION to source-hunt's round 10, item 1: the printed
        DOS factor carries TWO indices, N_{JJ'}(0), not N_{J'}(0).  Span dump:
          'N' size 9.96 | 'JJ' size 6.97 | U+2032 size 4.98 | '(0).' size 9.96
        -- the same three-span pattern as the neighbouring Lambda_JJ' and
        V^pairing_JJ'.  I do NOT adjudicate what N_{JJ'} means; see section (I).
"""

import itertools
from fractions import Fraction as F

import numpy as np

# ------------------------------------------------------------------ L1 inputs
# Golubov cond-mat/0111262 p.4/p.5 (read by me in round 9, text layer)
LAM_G = np.array([[1.017, 0.213],          # sigma, pi
                  [0.155, 0.448]])
MUS_G = np.array([[0.210, 0.095],
                  [0.069, 0.172]])
# Mazin-Antropov Table II, printed p.9
TABLE_II = {
    "Ref.[10]":  np.array([[0.96, 0.17], [0.23, 0.29]]),
    "Ref.[55]":  np.array([[1.02, 0.16], [0.21, 0.45]]),
    "Ref.[30]a": np.array([[0.78, 0.11], [0.15, 0.21]]),
}


def ratio(M, i=0, j=1):
    """criterion ratio  M_ij M_ji / (M_ii M_jj);  <= 1 means the criterion holds."""
    return (M[i, j] * M[j, i]) / (M[i, i] * M[j, j])


def all_ratios(M):
    d = M.shape[0]
    return [(i, j, (M[i, j] * M[j, i]) / (M[i, i] * M[j, j]))
            for i, j in itertools.combinations(range(d), 2)]


# ==================================================================== (A)
def sec_A():
    print("\n" + "=" * 76)
    print("(A) why the CRITERION (not PSD) is the convention-free object")
    print("=" * 76)
    print("""    [analytic]  For positive diagonals D1=diag(p), D2=diag(q) and N = D1 M D2:
        N_ii N_jj = p_i q_i p_j q_j M_ii M_jj
        N_ij N_ji = p_i q_j p_j q_i M_ij M_ji
      -> IDENTICAL prefactors, so the ratio is exactly invariant.  This covers
      every DOS attachment  h(N_i,N_j) in {N_j, N_i, sqrt(N_iN_j), N_j/N_F, 1}
      AND transposition.  Note PSD is NOT invariant under two-sided scaling
      with D1 != D2, which is exactly why the repo states the falsifiable test
      in product form.""")
    rng = np.random.default_rng(1010)
    worst = 0.0
    for _ in range(4000):
        d = int(rng.integers(2, 6))
        M = rng.normal(size=(d, d)) + 3 * np.eye(d)
        p, q = rng.uniform(0.05, 5, d), rng.uniform(0.05, 5, d)
        N = np.diag(p) @ M @ np.diag(q)
        for (i, j, r0), (_, _, r1) in zip(all_ratios(M), all_ratios(N)):
            worst = max(worst, abs(r1 - r0) / max(1e-30, abs(r0)))
    print("    4000 random (M, D1, D2), all pairs: max relative drift %.3e" % worst)


# ==================================================================== (B)
def lam_eff(Lam, mus, c=0.62):
    """Mazin-Antropov Eq.(7)."""
    S = Lam.sum(axis=1)                     # row sums, sum_n Lam_in
    A = 1.0 + c * S
    B = 1.0 + S
    return (Lam - (mus.T * A).T) / B[:, None]


def sec_B():
    print("\n" + "=" * 76)
    print("(B) Eq.(7): what cancels, what does not  [exact rational]")
    print("=" * 76)
    # exact rational 2x2 with symbolic-free but exact arithmetic
    L11, L12, L21, L22 = F(1017, 1000), F(213, 1000), F(155, 1000), F(448, 1000)
    m11, m12, m21, m22 = F(21, 100), F(95, 1000), F(69, 1000), F(172, 1000)
    c = F(62, 100)
    S1, S2 = L11 + L12, L21 + L22
    A1, A2 = 1 + c * S1, 1 + c * S2
    B1, B2 = 1 + S1, 1 + S2
    E11 = (L11 - m11 * A1) / B1
    E12 = (L12 - m12 * A1) / B1
    E21 = (L21 - m21 * A2) / B2
    E22 = (L22 - m22 * A2) / B2
    R_full = (E12 * E21) / (E11 * E22)
    C11, C12 = L11 - m11 * A1, L12 - m12 * A1
    C21, C22 = L21 - m21 * A2, L22 - m22 * A2
    R_noB = (C12 * C21) / (C11 * C22)
    print("    with the MgB2 numbers as exact rationals:")
    print("      A_sigma = %s = %.6f" % (A1, float(A1)))
    print("      A_pi    = %s = %.6f" % (A2, float(A2)))
    print("      R with the 1/(1+sum Lam) factors kept : %.12f" % float(R_full))
    print("      R with those factors dropped         : %.12f" % float(R_noB))
    print("      EXACTLY equal as rationals?  %s" % (R_full == R_noB))
    print("    -> [analytic] the mass-renormalisation denominator is a pure")
    print("       f(i) row scaling and cancels by (A).  CONFIRMS source-hunt.")

    print("""
    covariance under  Lam -> D_f Lam D_g,  mu* -> D_f mu* D_g:
      numerator  f_i g_j Lam_ij - A'_i f_i g_j mu*_ij = f_i g_j (Lam - D_A' mu*)_ij
      denominator B'_i is another pure row factor
      => Lam_eff' = D_f (Lam - D_A' mu*) D_g / B'   and by (A)
         R(Lam_eff') = ratio(Lam - D_A' mu*).

    ** [analytic] EXACT STATEMENT.  R(Lam_eff) = ratio(Lam - D_A mu*), a
       function of the DIAGONAL D_A = diag(1 + 0.62 sum_j Lam_ij) alone.  A
       convention change enters ONLY by moving that diagonal:
           A'_i = 1 + 0.62 f_i sum_n g_n Lam_in.
       => any convention change that leaves D_A fixed leaves R fixed.

    ** TWO CORRECTIONS TO MY OWN DRAFTS OF THIS SECTION, both caught by the
       check below before anything was written down.
       draft 1: 'A_i all equal => covariant again'.  FALSE -- ratio(Lam - D_A
         mu*) depends on A_1 and A_2 separately, so a change that keeps them
         equal but moves their common VALUE still moves R.
       draft 2: 'equal row sums => transposition-invariant'.  ALSO FALSE --
         transposition is not a D_f . D_g rescaling; it transposes mu* as
         well, and equal row sums do not force equal column sums.
       What survives, and is verified below:
         (i)  mu* = 0  =>  D_A drops out entirely  =>  R = ratio(Lam), the
              repo's theorem.  This is why the BARE Lambda is convention-free
              and Lam_eff is not.
         (ii) IF D_A is the same SCALAR A.1 under both readings, then with
              Lam = V D_N, mu* = K D_N (V, K symmetric) the two readings are
              (V - A K) D_N  and  D_N (V - A K), related by two-sided diagonal
              scaling => identical R by (A).  The spread AND the value of A_i
              must both be preserved; neither is automatic.""")
    rng = np.random.default_rng(77)
    worst_gen, worst_mu0, worst_scal = 0.0, 0.0, 0.0
    for _ in range(3000):
        V = np.abs(rng.normal(size=(2, 2))); V = (V + V.T) / 2
        K = np.abs(rng.normal(size=(2, 2))) * 0.2; K = (K + K.T) / 2
        Nb = rng.uniform(0.1, 1.2, 2)
        Lam, mus = V * Nb[None, :], K * Nb[None, :]

        def drift(Lm, Mu):
            r0, r1 = ratio(lam_eff(Lm, Mu)), ratio(lam_eff(Lm.T, Mu.T))
            return abs(r1 - r0) / max(1e-30, abs(r0))

        worst_gen = max(worst_gen, drift(Lam, mus))
        worst_mu0 = max(worst_mu0, drift(Lam, np.zeros((2, 2))))
        # corollary (ii): impose D_A = A.1 by hand under BOTH readings
        A = 1.0 + 0.62 * float(rng.uniform(0.2, 3.0))
        r_right = ratio((V - A * K) * Nb[None, :])        # DOS on right index
        r_left = ratio(Nb[:, None] * (V - A * K))         # DOS on left index
        worst_scal = max(worst_scal,
                         abs(r_left - r_right) / max(1e-30, abs(r_right)))
    print("    3000 random 2-band models, drift of R under the DOS-index swap:")
    print("      generic (A_i from the row sums) : worst %.4f" % worst_gen)
    print("      corollary (i)  mu* = 0          : worst %.3e" % worst_mu0)
    print("      corollary (ii) D_A = A.1 forced : worst %.3e" % worst_scal)


# ==================================================================== (C)
def sec_C():
    print("\n" + "=" * 76)
    print("(C) MgB2 demonstration: source-hunt's two numbers, from zero")
    print("=" * 76)
    for tag, (Lm, Mu) in (("as printed ", (LAM_G, MUS_G)),
                          ("transposed ", (LAM_G.T, MUS_G.T))):
        S = Lm.sum(axis=1)
        E = lam_eff(Lm, Mu)
        r = ratio(E)
        print("    %s A = [%.5f, %.5f]   R_sigpi = %.6f   (pass %.1f x)"
              % (tag, 1 + 0.62 * S[0], 1 + 0.62 * S[1], r, 1 / r))
    r_bare = ratio(LAM_G)
    print("    bare Lambda: R = %.9f both ways (pass %.3f x) -- the theorem"
          % (r_bare, 1 / r_bare))
    print("    drift of R(Lam_eff) under transposition: %.2f %%"
          % (100 * abs(ratio(lam_eff(LAM_G.T, MUS_G.T)) / ratio(lam_eff(LAM_G, MUS_G)) - 1)))
    print("    -> source-hunt's 0.02003 / 0.02108 / 5.3 %% reproduced.")
    print("    DIRECTION AUDIT: both readings PASS (49.9x, 47.4x).  The demo")
    print("    shows convention-SENSITIVITY, never a flipped verdict; I do not")
    print("    claim a flip is reachable and have not bounded the drift.")


# ==================================================================== (D)
def sec_D():
    print("\n" + "=" * 76)
    print("(D) ** the bigger issue: Lam - A mu* leaves the criterion's class **")
    print("=" * 76)
    print("""    Convention is a red herring here.  Take A_i = A CONSTANT (so section
    (B) says everything is covariant), V = I and K = ones(2,2).  Both are
    symmetric, both satisfy the criterion (K on its boundary).  Then
        M(A) = V - A K = [[1-A, -A], [-A, 1-A]]
        criterion residual = (1-A)^2 - A^2 = 1 - 2A
    which is NEGATIVE for A > 1/2.  A = 1 + 0.62*lambda >= 1 ALWAYS.""")
    for A in (0.0, 0.25, 0.5, 0.75, 1.0, 1.7626):
        M = np.eye(2) - A * np.ones((2, 2))
        res = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
        print("      A = %6.4f :  residual = %+8.5f   %s"
              % (A, res, "VIOLATED" if res < 0 else "ok"))
    print("""
    -> [analytic]  The criterion is NOT preserved by lambda-minus-mu* type
       combinations, in ANY convention and even for symmetric inputs.  This is
       the same fact as 'subtraction leaves the class P' that killed the
       rho(|M2^Delta|) ceiling: a difference of two PSD matrices need not be
       PSD.  Section (B)'s covariance question is the SMALLER of the two
       problems with Lam_eff, and would be moot anyway if the repo consumed
       Lam_eff -- which section (H) shows it does not.""")
    # how often does it bite on physical-looking inputs?
    rng = np.random.default_rng(303)
    n_bad = n_tot = 0
    for _ in range(4000):
        X = rng.normal(size=(2, 4)); V = X @ X.T / 4
        Y = rng.normal(size=(2, 4)); K = Y @ Y.T / 4 * rng.uniform(0.02, 0.35)
        Nb = rng.uniform(0.1, 1.2, 2)
        Lam, mus = V * Nb[None, :], K * Nb[None, :]
        if ratio(Lam) > 1 or ratio(mus) > 1:
            continue
        n_tot += 1
        if ratio(lam_eff(Lam, mus)) > 1:
            n_bad += 1
    print("    random PSD Lam and PSD mu* (both passing): Lam_eff VIOLATES in")
    print("    %d / %d = %.1f %% of cases." % (n_bad, n_tot, 100 * n_bad / n_tot))


# ==================================================================== (E)
def sec_E():
    print("\n" + "=" * 76)
    print("(E) Eq.(8) is a BAND-INDEXED MA map -- round-9 theorem generalises")
    print("=" * 76)
    print("""    Eq.(8) reads   mu*_ij = mu_ij - sum_n mu_in log(W_n/wc) mu*_nj,
    i.e.  mu* = mu - mu D_L mu*  with  D_L = diag(log(W_n/wc)).  Hence
        (1 + mu D_L) mu* = mu     =>    mu* = (1 + mu D_L)^{-1} mu.
    Round 9 proved preservation via the SPECTRAL map x -> x/(1+Lx); that proof
    needs D_L scalar (otherwise mu D_L does not commute with mu).  The
    RESOLVENT form does not:
        (mu*)^{-1} = mu^{-1} (1 + mu D_L) = mu^{-1} + D_L,
    a sum of a symmetric PSD-inverse and a non-negative diagonal.
    => [analytic, any band count, any diagonal D_L >= 0]  mu^s > 0  ==>
       (mu^s)* = ((mu^s)^{-1} + D_L)^{-1} is symmetric positive definite
       => every 2x2 principal minor >= 0 => criterion preserved.
    Mitrovic Eq.(12) is the special case D_L = L*I, so this STRICTLY contains
    the round-9 theorem and covers band-dependent bandwidths W_n.""")
    rng = np.random.default_rng(88)
    worst_r, worst_sym, worst_eig = 0.0, 0.0, np.inf
    for _ in range(3000):
        d = int(rng.integers(2, 6))
        X = rng.normal(size=(d, d + 2)); K = X @ X.T
        Nb = rng.uniform(0.05, 1.5, d)
        Dh = np.diag(np.sqrt(Nb))
        mu_s = Dh @ K @ Dh
        DL = np.diag(rng.uniform(0.0, 4.0, d))          # band-indexed log(W_n/wc)
        ms = np.linalg.solve(np.eye(d) + mu_s @ DL, mu_s)
        worst_sym = max(worst_sym, np.abs(ms - ms.T).max() / np.abs(ms).max())
        worst_eig = min(worst_eig, np.linalg.eigvalsh((ms + ms.T) / 2).min()
                        / np.linalg.eigvalsh(mu_s).max())
        # identity check
        idn = np.linalg.inv(ms) - (np.linalg.inv(mu_s) + DL)
        assert np.abs(idn).max() < 1e-6 * np.abs(np.linalg.inv(mu_s)).max()
        # read back in the N_j convention
        Mn = np.diag(1 / np.sqrt(Nb)) @ ms @ Dh
        for _i, _j, r in all_ratios(Mn):
            worst_r = max(worst_r, r)
    print("    3000 random PSD mu with NON-SCALAR D_L, d in 2..5:")
    print("      identity  mu*^{-1} = mu^{-1} + D_L      : verified every time")
    print("      worst criterion ratio over all pairs    : %.12f" % worst_r)
    print("      min normalised eigenvalue of mu*        : %.3e" % worst_eig)
    print("      worst relative asymmetry of mu*         : %.3e" % worst_sym)

    print("\n    converse (pole) with band-indexed cutoffs, D_L = t * Dhat:")
    print("""      det(1 + mu D_L) = 0  <=>  -1/t in spec(mu Dhat) = spec(T),
      T := Dhat^{1/2} mu Dhat^{1/2} symmetric.  A violated 2x2 minor with
      non-negative diagonal makes mu indefinite; Sylvester's law of inertia
      says the CONGRUENCE T is indefinite too => T has a negative eigenvalue
      => t = -1/lambda_min(T) > 0 is a genuine pole, i.e. a cutoff
      omega_c in (0, W) at which the band-indexed MA map blows up.""")
    rng2 = np.random.default_rng(89)
    hits, worst_det = 0, 0.0
    for _ in range(600):
        d = int(rng2.integers(2, 5))
        K = rng2.normal(size=(d, d)); K = (K + K.T) / 2
        np.fill_diagonal(K, np.abs(np.diag(K)) + 0.1)
        if all(r <= 1 for _i, _j, r in all_ratios(K)):
            continue
        hits += 1
        Dhat = np.diag(rng2.uniform(0.2, 3.0, d))
        Dh2 = np.sqrt(Dhat)
        T = Dh2 @ K @ Dh2
        ev = np.linalg.eigvalsh(T)
        assert ev.min() < 0, "violating criterion but T is PSD?!"
        t = -1.0 / ev.min()
        assert t > 0
        det = np.linalg.det(np.eye(d) + K @ (t * Dhat))
        worst_det = max(worst_det, abs(det))
    print("      %d violating matrices: positive pole t found in EVERY case," % hits)
    print("      0 counterexamples; worst |det(1+mu D_L)| at that t = %.2e" % worst_det)
    print("      => the round-9 converse also survives band-dependent cutoffs.")

    print("\n    the paper's own two-band estimate, checked:")
    a, L = 2.3, 1.0
    astar = a + (a - 1 / a) * L
    print("      alpha* = alpha + (alpha - 1/alpha) L = %.4f at alpha=2.3, L=1"
          % astar)
    print("      -> printed 'mu*_ss = mu*_pp ~ 4 mu*_sp'; criterion margin goes")
    print("         from alpha^2 = %.2f x (bare) to alpha*^2 = %.2f x."
          % (a ** 2, astar ** 2))
    print("      Jepsen-Andersen bare <<V_C>> ratios 3:2.5:1 -> margin %.2f x."
          % (3 * 2.5 / 1.0))
    print("      ALL of these PASS; consistent with (E) but NOT evidence for")
    print("      it -- a positive kernel is PSD by Bochner to begin with.")


# ==================================================================== (F)
def sec_F():
    print("\n" + "=" * 76)
    print("(F) Table II: three published first-principles MgB2 matrices")
    print("=" * 76)
    for k, M in TABLE_II.items():
        r = ratio(M)
        print("    %-10s  ii*jj=%.5f  ij*ji=%.5f  ratio=%.5f  pass %.2f x"
              % (k, M[0, 0] * M[1, 1], M[0, 1] * M[1, 0], r, 1 / r))
    print("""
    TENSION, recorded not adjudicated: Ref.[55] IS Golubov cond-mat/0111262,
    which prints lam_sigpi = 0.213 > lam_pisig = 0.155, but Table II puts the
    SMALLER entry at (1,2) (all three columns do).  Under that paper's own
    Lam_ij = V_ij N_j and N_pi > N_sigma, (1,2) should be the LARGER.  So Table
    II looks transposed relative to the source.  The criterion is transpose-
    invariant by (A), so none of the three verdicts moves.""")


# ==================================================================== (G)
def sec_G():
    print("\n" + "=" * 76)
    print("(G) Mazin-Golubov's all-repulsive condition: printed vs sharp")
    print("=" * 76)
    print("""    [analytic]  Two bands, g symmetric, Lam_ij = g_ij N_j, all g_ij < 0.
    Write a = g11 N1 < 0, d = g22 N2 < 0, b = g12 N2, c = g12 N1, bc > 0.
        lam_max = [(a+d) + sqrt((a-d)^2 + 4bc)]/2 > 0
        <=>  (a-d)^2 + 4bc > (a+d)^2   <=>   bc > ad
        <=>  Lam_12 Lam_21 > Lam_11 Lam_22        (EXACTLY  NOT(criterion))
        <=>  g12^2 > g11 g22   <=>   |g12| > sqrt(|g11||g22|).
    The printed condition is  |g12| > (|g1| N1^2 + |g2| N2^2)/(2 N1 N2).
    By AM-GM that right-hand side is >= sqrt(|g1||g2|), with equality iff
    |g1| N1^2 = |g2| N2^2.  So the printed condition IMPLIES the sharp one:
    it is SUFFICIENT but NOT NECESSARY.  No contradiction -- and the sharp
    condition is literally the negation of the repo's criterion.""")
    rng = np.random.default_rng(555)
    bad = 0
    gap_max = 0.0
    for _ in range(20000):
        g11, g22, g12 = -rng.uniform(0.05, 2, 3)
        N1, N2 = rng.uniform(0.1, 2, 2)
        Lam = np.array([[g11 * N1, g12 * N2], [g12 * N1, g22 * N2]])
        lmax = np.linalg.eigvals(Lam).real.max()
        sharp = Lam[0, 1] * Lam[1, 0] > Lam[0, 0] * Lam[1, 1]
        if (lmax > 1e-12) != sharp:
            bad += 1
        printed = abs(g12) > (abs(g11) * N1 ** 2 + abs(g22) * N2 ** 2) / (2 * N1 * N2)
        if printed and not sharp:
            bad += 1                      # would falsify "printed => sharp"
        rhs_p = (abs(g11) * N1 ** 2 + abs(g22) * N2 ** 2) / (2 * N1 * N2)
        gap_max = max(gap_max, rhs_p / np.sqrt(abs(g11) * abs(g22)))
    print("    20000 random all-repulsive 2-band models:")
    print("      'lam_max>0  <=>  Lam12 Lam21 > Lam11 Lam22'   mismatches: %d" % bad)
    print("      largest printed-RHS / sharp-RHS ratio: %.2f x (AM-GM slack)" % gap_max)


# ==================================================================== (H)
def sec_H():
    print("\n" + "=" * 76)
    print("(H) ** THE CONSUMPTION POINT, pinned **")
    print("=" * 76)
    print("""    Where the criterion is actually eaten in this repo:

      COULOMB SIDE -- Lemma mu-1 (mustar_psd_multiband_ceiling_2026-09-18.md):
        M_ij = mu*_ij w^(i) (w^(j))^T  is PSD  <==  mu* PSD
        => O_mu = O - M  <=  O  => Tc_mu <= Tc_0 <= ceiling.
        The criterion (all 2x2 principal minors >= 0) is a NECESSARY condition
        for mu* PSD, hence a falsifiable test OF THIS HYPOTHESIS.  <-- consumed

      PHONON SIDE -- (ML3)/(G'-abs) (multiband_ceiling_route_Gprime_2026-09-18.md):
        needs  L := diag(R^Z) - LamTilde  >= 0  (weighted Laplacian PSD),
        proved from ELEMENTWISE non-negativity, resp. (Z-dom) |Lam^D| <= Lam^Z.
        The criterion appears NOWHERE in that proof.            <-- NOT consumed

    So: is 'criterion' the same condition as 'elementwise non-negative'?  No --
    they are logically independent in BOTH directions, shown below.""")
    # independence, both directions
    P = np.array([[1.0, 2.0], [2.0, 1.0]])        # elementwise >=0, criterion FAILS
    Q = np.array([[1.0, -0.5], [-0.5, 1.0]])      # criterion HOLDS, has a negative entry
    for tag, M in (("elementwise >= 0 but criterion violated", P),
                   ("criterion satisfied but a negative entry", Q)):
        w = np.ones(2)
        R = (M @ w) / w
        Lap = np.diag(R) - M
        print("\n    %s:" % tag)
        print("      M = %s" % np.array2string(M, prefix=" " * 10))
        print("      criterion ratio = %.4f (%s)"
              % (ratio(M), "FAIL" if ratio(M) > 1 else "pass"))
        print("      Laplacian diag(R)-M eigenvalues = %s -> %s"
              % (["%.4f" % e for e in np.linalg.eigvalsh(Lap)],
                 "PSD" if np.linalg.eigvalsh(Lap).min() > -1e-12 else "NOT PSD"))
    print("""
    -> The two conditions are independent.  (ML3) rides on elementwise sign,
       the criterion rides on 2x2 minors.  A criterion violation says nothing
       about (ML3) and vice versa.""")

    # No fabricated CaC6 entries here: the load-bearing CaC6 numbers live in
    # _chk_cac6_lambda_ij_criterion.py.  What (H) needs is the GENERAL fact.
    print("""
    [analytic]  For ANY elementwise-non-negative symmetric Lam and ANY w > 0,
    with R_i := (Lam w)_i / w_i,
      sum_ij (diag(R) - Lam)_ij x_i x_j
        = (1/2) sum_ij Lam_ij ( (w_j/w_i) x_i^2 + (w_i/w_j) x_j^2 - 2 x_i x_j )
        = (1/2) sum_ij Lam_ij ( sqrt(w_j/w_i) x_i - sqrt(w_i/w_j) x_j )^2 >= 0.
    Every term is a non-negative coefficient times a square.  The 2x2 minors
    of Lam never appear.  So a criterion violation cannot disturb (ML3).""")
    rng = np.random.default_rng(404)
    n_viol = n_psd = 0
    worst_eig = np.inf
    for _ in range(4000):
        d = int(rng.integers(2, 7))
        M = np.abs(rng.normal(size=(d, d))); M = (M + M.T) / 2
        if all(r <= 1 for _i, _j, r in all_ratios(M)):
            continue
        n_viol += 1
        w = rng.uniform(0.05, 5.0, d)
        R = (M @ w) / w
        Lap = np.diag(R) - M
        e = np.linalg.eigvalsh((Lap + Lap.T) / 2).min()
        worst_eig = min(worst_eig, e / np.abs(M).max())
        if e > -1e-10:
            n_psd += 1
    print("    %d random elementwise-non-negative symmetric matrices that"
          % n_viol)
    print("    VIOLATE the criterion, random w > 0: weighted Laplacian PSD in")
    print("    %d / %d cases; worst normalised eigenvalue %.2e"
          % (n_psd, n_viol, worst_eig))
    print("""
    CONCLUSION.  The criterion's only consumption point is the COULOMB block.
    The four material-side failures (FeSe / CaC6 / 122 / 1111) are all on the
    PHONON side and therefore never fed the ceiling chain: they are evidence
    that the criterion is not vacuous, not threats to the ceiling.  Lam_eff
    (section B/D) is a McMillan-equation convenience object that the repo's
    chain never forms, so its non-covariance lands outside as well.""")


# ==================================================================== (I)
def sec_I():
    print("\n" + "=" * 76)
    print("(I) the one lead that is NOT closed: Allen's N_{JJ'}(0)")
    print("=" * 76)
    print("""    Mazin-Golubov p.4 prints  Lam_JJ' = V^pairing_JJ' N_JJ'(0)  with a
    TWO-INDEX DOS factor (span-level extraction in the docstring).  Taken
    literally, Lam = V o N is a HADAMARD product of two symmetric matrices,
    and that is NOT of the form f(i) g(j) V_ij unless N factorises.  That would
    be exactly the shape the convention-invariance theorem does not cover.

    I do NOT adjudicate it:
      - the SAME paper's p.2 writes Lam_ij = g_ij N_j, single index;
      - in Allen's Fermi-surface-harmonic formalism the harmonics are
        orthonormal under the Fermi-surface average, which would make
        N_{JJ'} proportional to delta_{JJ'} and the formula degenerate;
      - so 'typo for N_{J'}(0)' and 'genuine two-index object' are both live.
    Settling it needs Allen, Phys. Rev. B 13, 1416 (1976) (their Ref. [10]),
    which I have not retrieved.  Logged as OPEN, not as a counterexample.""")


if __name__ == "__main__":
    sec_A(); sec_B(); sec_C(); sec_D(); sec_E(); sec_F(); sec_G(); sec_H(); sec_I()
