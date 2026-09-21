#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_mustar_MA_matrix_criterion.py -- 2026-09-18, main session, from zero.

QUESTION (handed to me by source-hunt, 9th round, explicitly undecided by it):
  Mitrovic's multiband Morel-Anderson renormalisation
        mu*(wc) = (1 + mu^s ln(E_F/wc))^{-1} mu^s                      Eq.(12)
  is NOT a rank-one rescaling mu_ij -> f(i)g(j)mu_ij, so the repo's
  convention-invariance theorem does NOT cover it.  Does the multiband
  criterion   M_ii M_jj >= M_ij M_ji   survive this map?

ANSWER PROVED HERE:
  (i)  [analytic, any number of bands]  YES, given mu^s PSD -- which the repo's
       Bochner/Gram argument supplies.  Reason: mu^s symmetric PSD => the map
       is the SPECTRAL map  x -> x/(1+Lx)  on mu^s, so mu*(wc) is symmetric PSD
       for every wc, hence every 2x2 principal minor is >= 0.
  (ii) [analytic, any number of bands]  CONVERSELY, if the criterion is
       VIOLATED on any pair, then Eq.(12) itself has a POLE at a finite,
       physical cutoff  wc = E_F exp(-1/|mu_min|)  in (0, E_F).
       => On the Coulomb side the criterion is not an extra assumption: it is
       exactly the condition for the standard multiband MA machinery to be
       well defined.  (Two bands: equivalence is exact.  d>=3: only this
       direction; the converse fails, see (E).)
  (iii)[analytic] The RATIO mu*_ij mu*_ji/(mu*_ii mu*_jj) DOES move with wc
       (source-hunt's observation is right) -- it just never crosses 1.

L1 SOURCES, read by the main session from the TEXT LAYER of PDFs it fetched
itself (both real text PDFs, not scans):

  _mitrovic0402473.pdf = B. Mitrovic, "Functional derivatives of Tc for a
      two-band superconductor: application to MgB2", cond-mat/0402473v1, 6 pp.
    p.2 verbatim:
      "one first takes the cutoff wc in the Matsubara frequency sums to be
       large enough so that mu*_ij(wc) = mu_ij = V^c_ij N_j, where V^c_ij is
       the Fermi surface averaged screened Coulomb matrix element between the
       states in the bands i and j; clearly V^c_ji = V^c_ij.  The electron-
       phonon coupling functions can be written as
       alpha^2 F_ij(Om) = alpha^2 f_ij(Om) N_j  with  a^2f_ij = a^2f_ji."
      >>> THIRD-PARTY L1 CONFIRMATION of the repo's  mu_ij = N_j K_ij  with K
          symmetric, AND of lambda_ij = N_j V_ij (DOS on the RIGHT index),
          i.e. the Sanna convention, stated for BOTH channels in one sentence.
      Eq.(10) verbatim:   mu^s_ij = sqrt(N_i/N_j) mu_ij
      >>> = sqrt(N_i N_j) V^c_ij  =  D^{1/2} V^c D^{1/2}.  The paper performs
          the symmetrisation itself; my "similarity to a symmetric matrix"
          argument is its own Eqs.(7)-(10).
      Eq.(12) verbatim:  hat-mu*(wc) = [hat-1 + hat-mu^s ln(E_F/wc)]^{-1} hat-mu^s
      Eqs.(13)-(16), two bands, verbatim (with "mu^s_pi-sigma = mu^s_sigma-pi"
      printed in the lead-in parenthesis):
        mu*_ss = [a + (ac - b^2) L]/D ,  mu*_sp = b/D ,
        mu*_pp = [c + (ac - b^2) L]/D ,  L := ln(E_F/wc)
        "where D is the determinant of hat-1 + hat-mu^s ln(E_F/wc)"
        D = 1 + (a+c) L  ?  (ac-b^2) L^2                                  (16)
      >>> TRANSCRIPTION NOTE: in the text layer the operator before the L^2
          term of Eq.(16) comes out as a multiplication sign.  It must be "+":
          the line above Eq.(16) DEFINES D as det(1 + mu^s L), and
          det(1+A L) = 1 + tr(A) L + det(A) L^2 for 2x2.  Checked symbolically
          in section (A) below rather than read off the glyph.

  _golubov_cond-mat0111262.pdf = Golubov, Kortus, Kong, Dolgov, Jepsen, Shulga,
      Drechsler, Kremer, "Specific heat of MgB2 ... " (= J. Phys. Cond. Mat.
      14, 1353), cond-mat/0111262v2, 9 pp.
    printed p.4 (page idx 3) verbatim:
      "Their respective densities of states at the Fermi energy have values of
       N_sigma(0) = 0.300 states/cell*eV and N_pi(0) = 0.410 states/cell*eV"
      "lambda_ss = 1.017, lambda_pp = 0.448, lambda_sp = 0.213, and
       lambda_ps = 0.155."
    printed p.5 (page idx 4) verbatim:
      "The sigma-sigma-, pi-pi- and sigma-pi-values were in the ratio
       2.23/2.48/1."
      "mu*_ss(wc)=0.210, mu*_sp(wc)=0.095, mu*_ps(wc)=0.069, mu*_pp(wc)=0.172."
    >>> This CLOSES the repo's long-standing "MgB2 source unread" marker: the
        four lambda_ij the repo quotes are verbatim correct.  Upgrades that row
        from "transcribed, source unread" to L1.
"""

import itertools

import numpy as np

# ----------------------------------------------------------------- L1 inputs
N_S, N_P = 0.300, 0.410                      # states/cell*eV, Golubov p.4
LAM = np.array([[1.017, 0.213],              # rows/cols order: sigma, pi
                [0.155, 0.448]])
MUSTAR = np.array([[0.210, 0.095],
                   [0.069, 0.172]])


def minors(M, labels=None):
    """all 2x2 principal minors  M_ii M_jj - M_ij M_ji  (the criterion)."""
    d = M.shape[0]
    labels = labels or [str(i + 1) for i in range(d)]
    out = []
    for i, j in itertools.combinations(range(d), 2):
        lhs, rhs = M[i, i] * M[j, j], M[i, j] * M[j, i]
        out.append((labels[i], labels[j], lhs, rhs, lhs - rhs))
    return out


def report(M, tag, labels=None):
    print("  %s" % tag)
    worst = None
    for a, b, lhs, rhs, m in minors(M, labels):
        r = rhs / lhs if lhs else float("inf")
        print("    (%s,%s)  ii*jj=%-12.6g  ij*ji=%-12.6g  minor=%-12.6g  %s"
              % (a, b, lhs, rhs, m,
                 "FAIL %.4g x" % r if m < 0 else "pass (margin %.4g x)" % (1 / r)))
        worst = r if worst is None else max(worst, r)
    return worst


# ==================================================================== (A)
def sec_A():
    print("\n" + "=" * 74)
    print("(A) Eq.(16): is the disputed operator '+' ?  [symbolic, exact]")
    print("=" * 74)
    from fractions import Fraction as F
    a, c, b, L = F(7, 3), F(5, 2), F(4, 5), F(9, 7)      # arbitrary exact rationals
    A = [[a, b], [b, c]]
    M = [[1 + A[0][0] * L, A[0][1] * L], [A[1][0] * L, 1 + A[1][1] * L]]
    det = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    plus = 1 + (a + c) * L + (a * c - b * b) * L * L
    times = 1 + (a + c) * L * (a * c - b * b) * L * L
    print("    det(1 + mu^s L)                 = %s" % det)
    print("    1 + (a+c)L + (ac-b^2)L^2        = %s   %s" % (plus, "MATCH" if plus == det else "no"))
    print("    1 + (a+c)L x (ac-b^2)L^2        = %s   %s" % (times, "MATCH" if times == det else "no"))
    print("    -> the text-layer glyph is a mis-rendered '+'.  Read as '+'.")

    # and the closed form of the criterion residual, exactly
    D = det
    m_ss = (a + (a * c - b * b) * L) / D
    m_sp = b / D
    m_pp = (c + (a * c - b * b) * L) / D
    resid = m_ss * m_pp - m_sp * m_sp
    print("\n    criterion residual of mu* (Eqs.13-15), exact rational:")
    print("      mu*_ss mu*_pp - mu*_sp^2 = %s" % resid)
    print("      Delta / D                = %s   %s"
          % ((a * c - b * b) / D, "MATCH" if resid == (a * c - b * b) / D else "no"))
    print("    -> [analytic]  det mu* = Delta/D, and for 2x2 the criterion")
    print("       residual EQUALS det.  So sign(criterion for mu*) =")
    print("       sign(Delta) whenever D > 0, i.e. the criterion is carried")
    print("       over from mu^s unchanged.")


# ==================================================================== (B)
def sec_B():
    print("\n" + "=" * 74)
    print("(B) MgB2, Golubov's printed numbers: criterion on BOTH channels")
    print("=" * 74)
    lab = ["sigma", "pi"]
    print("\n  electron-phonon lambda_ij (first-principles):")
    wl = report(LAM, "", lab)
    print("\n  Coulomb mu*_ij(wc) (fitted to Tc = 39.4 K):")
    wm = report(MUSTAR, "", lab)
    print("\n    -> MgB2 PASSES on both channels: lambda margin %.3f x," % (1 / wl))
    print("       mu* margin %.3f x.  (Repo previously recorded 13.8 x for" % (1 / wm))
    print("       lambda with the source UNREAD; now verbatim, and the mu*")
    print("       row is new.)")

    print("\n  convention self-check  M_ij/M_ji  vs  N_j/N_i  (= Sanna/Mitrovic")
    print("  convention M_ij = N_j V_ij with V symmetric):")
    print("    N_pi/N_sigma                 = %.6f" % (N_P / N_S))
    print("    lambda_sp/lambda_ps          = %.6f   (%.2f %% off)"
          % (LAM[0, 1] / LAM[1, 0], 100 * abs(LAM[0, 1] / LAM[1, 0] / (N_P / N_S) - 1)))
    print("    mu*_sp/mu*_ps                = %.6f   (%.2f %% off)"
          % (MUSTAR[0, 1] / MUSTAR[1, 0], 100 * abs(MUSTAR[0, 1] / MUSTAR[1, 0] / (N_P / N_S) - 1)))
    print("    -> both off-diagonal ratios reproduce N_pi/N_sigma to < 1 %,")
    print("       from numbers printed one page apart.  This is an internal")
    print("       consistency check of the transcription AND an independent")
    print("       confirmation of the DOS-on-the-right-index convention --")
    print("       which Mitrovic p.2 states in words for this very paper.")


# ==================================================================== (C)
def ma_map(mu_s, L):
    """Mitrovic Eq.(12) on the symmetrised matrix."""
    d = mu_s.shape[0]
    return np.linalg.solve(np.eye(d) + mu_s * L, mu_s)


def sec_C():
    print("\n" + "=" * 74)
    print("(C) the theorem, numerically: PSD mu^s => criterion at every wc")
    print("=" * 74)
    print("    mu^s symmetric PSD => the map is the SPECTRAL map x -> x/(1+Lx)")
    print("    on mu^s (same eigenvectors), so every eigenvalue stays >= 0 and")
    print("    mu*(wc) is symmetric PSD for every L >= 0.  Checked below, and")
    print("    the verdict is then read back in the N_j convention.")
    rng = np.random.default_rng(20260918)
    worst_over_all, worst_sym_err = -np.inf, 0.0
    worst_eig, rmin, rmax, drift = np.inf, np.inf, -np.inf, 0.0
    for d in (2, 3, 4, 6):
        for _ in range(400):
            # PSD K by construction (Gram form -- the repo's Bochner input)
            X = rng.normal(size=(d, d + 2))
            K = X @ X.T
            Nb = rng.uniform(0.05, 1.5, size=d)
            Dh = np.diag(np.sqrt(Nb))
            mu_s = Dh @ K @ Dh                       # Eq.(10)
            seen = []
            for L in (0.0, 0.3, 1.0, 3.0, 10.0, 100.0):
                MS = ma_map(mu_s, L)
                worst_sym_err = max(worst_sym_err, np.abs(MS - MS.T).max())
                worst_eig = min(worst_eig, np.linalg.eigvalsh((MS + MS.T) / 2).min())
                # back to the N_j convention, as Golubov prints it
                Mn = np.diag(1 / np.sqrt(Nb)) @ MS @ Dh
                for _i, _j, lhs, rhs, m in minors(Mn):
                    worst_over_all = max(worst_over_all, rhs / lhs)
                r = MS[0, 1] * MS[1, 0] / (MS[0, 0] * MS[1, 1])
                seen.append(r)
            rmin, rmax = min(rmin, min(seen)), max(rmax, max(seen))
            drift = max(drift, max(seen) - min(seen))
    print("\n    1600 random PSD Coulomb matrices, d in {2,3,4,6}, 6 cutoffs each")
    print("    worst  M_ij M_ji / (M_ii M_jj)  over ALL pairs, d, wc :  %.12f" % worst_over_all)
    print("    min eigenvalue of mu*(wc) over everything             :  %.3e" % worst_eig)
    print("    max |mu* - mu*^T| (symmetry of the map's output)      :  %.3e" % worst_sym_err)
    print("    (1,2) ratio: range over all runs [%.4f, %.4f], largest" % (rmin, rmax))
    print("    within-matrix swing across the six cutoffs            :  %.4f" % drift)
    print("    -> criterion never violated (ratio <= 1), as the spectral-map")
    print("       argument requires; the RATIO does move with wc (source-hunt's")
    print("       observation is correct, swing up to %.2f here) but is pinned" % drift)
    print("       below 1.  Verdict preserved in BOTH conventions (checked")
    print("       after mapping back to M_ij = N_j V_ij).")


# ==================================================================== (D)
def sec_D():
    print("\n" + "=" * 74)
    print("(D) converse: violated criterion => Eq.(12) has a pole at finite wc")
    print("=" * 74)
    print("    [analytic]  a violated 2x2 principal minor with non-negative")
    print("    diagonal makes that 2x2 submatrix indefinite; by Cauchy")
    print("    interlacing mu^s then has an eigenvalue m < 0, and")
    print("        det(1 + mu^s L) = prod_k (1 + m_k L)")
    print("    vanishes at  L = 1/|m| > 0, i.e. at the physical cutoff")
    print("        wc = E_F exp(-1/|m|)  in (0, E_F).")
    rng = np.random.default_rng(7)
    hits = 0
    for _ in range(300):
        d = rng.integers(2, 5)
        K = rng.normal(size=(d, d)); K = (K + K.T) / 2
        np.fill_diagonal(K, np.abs(np.diag(K)) + 0.1)
        Nb = rng.uniform(0.1, 1.0, size=d)
        Dh = np.diag(np.sqrt(Nb))
        mu_s = Dh @ K @ Dh
        viol = any(m < 0 for *_x, m in minors(mu_s))
        if not viol:
            continue
        hits += 1
        ev = np.linalg.eigvalsh(mu_s)
        assert ev.min() < 0, "violated criterion but PSD?!"
        Lpole = 1.0 / abs(ev.min())
        assert abs(np.linalg.det(np.eye(d) + mu_s * Lpole)) < 1e-8 * max(1, abs(np.linalg.det(mu_s)))
    print("    checked on %d random VIOLATING Coulomb matrices: in every case" % hits)
    print("    mu^s has a negative eigenvalue and det(1+mu^s L) vanishes at")
    print("    L = 1/|m_min|, to machine precision.  0 counterexamples.")
    print("    => the Coulomb-side criterion is not an extra assumption; it is")
    print("       exactly the well-definedness of the standard multiband MA")
    print("       renormalisation for all cutoffs.")


# ==================================================================== (E)
def sec_E():
    print("\n" + "=" * 74)
    print("(E) the honest limitation: at d >= 3 the criterion alone is WEAKER")
    print("    than PSD, and the MA map can then BREAK it")
    print("=" * 74)
    K = np.array([[1.0, 0.9, -0.9],
                  [0.9, 1.0, 0.9],
                  [-0.9, 0.9, 1.0]])
    print("    K =\n%s" % np.array2string(K, prefix="        "))
    print("    all 2x2 principal minors of K: %s"
          % ["%.4f" % m for *_x, m in minors(K)])
    print("    eigenvalues of K             : %s"
          % ["%.4f" % e for e in np.linalg.eigvalsh(K)])
    print("    -> criterion HOLDS pairwise, but K is INDEFINITE.")
    for L in (0.0, 0.2, 0.5, 0.8):
        MS = ma_map(K, L)
        w = max(rhs / lhs for *_x, lhs, rhs, _m in
                [(a, b, l_, r_, m_) for a, b, l_, r_, m_ in minors(MS)])
        print("      L = %4.2f :  worst ij*ji/(ii*jj) = %8.4f   %s"
              % (L, w, "FAIL" if w > 1 else "pass"))
    print("    => the preservation theorem needs mu^s PSD, NOT merely 'mu^s")
    print("       satisfies the criterion'.  The repo's Bochner/Gram argument")
    print("       does deliver PSD, so the theorem applies -- but the stronger")
    print("       hypothesis must be carried with the statement.")


# ==================================================================== (F)
def sec_F():
    print("\n" + "=" * 74)
    print("(F) what Golubov's printed mu* actually IS  --  and a false start")
    print("=" * 74)
    print("""
    >>> A FALSE START OF MINE, recorded because the shape matters (CLAUDE.md
        W1: the sentence that says 'the other party's numbers do not work' is
        itself the one that skips verification).
        I first read 'the ratios 2.23/2.48/1' as ratios of the BARE, symmetric
        V^c, and tested whether Mitrovic Eq.(12) maps them onto the printed
        mu*.  Best fit over (scale, L) was 25.7 % off, with sigma-sigma and
        pi-pi INVERTED relative to print -- and the inversion is structural:
        for mu_ij = N_j V_ij the map sends mu*_ss/mu*_pp monotonically to 1
        from whichever side it starts, so NO (s, L) can reproduce
        mu*_ss > mu*_pp when N_s V_ss < N_p V_pp.  I was one line away from
        writing 'Golubov's mu* does not come out of the MA formula'.
        THE HYPOTHESIS WAS MINE, NOT THE PAPER'S.  Verbatim, p.5:
          "we can approximately calculate the ratios for the MU-MATRIX [37].
           ... This allows one to express mu*_ij(wc) by these ratios and one
           single free parameter which is fixed to get the experimental Tc"
        -- the ratios are ratios OF mu* itself, fixed by one scale.  There is
        no bare -> renormalised step in this paper to test at all.
    """)
    ss, sp, ps, pp = MUSTAR[0, 0], MUSTAR[0, 1], MUSTAR[1, 0], MUSTAR[1, 1]
    print("    which pairing do the printed numbers support?")
    print("      mu*_ss/mu*_sp = %.4f   vs 2.23  -> %.2f %%   (row-wise reading)"
          % (ss / sp, 100 * abs(ss / sp / 2.23 - 1)))
    print("      mu*_pp/mu*_ps = %.4f   vs 2.48  -> %.2f %%   (row-wise reading)"
          % (pp / ps, 100 * abs(pp / ps / 2.48 - 1)))
    print("      mu*_pp/mu*_sp = %.4f   vs 2.48  -> %.2f %%   (single-reference)"
          % (pp / sp, 100 * abs(pp / sp / 2.48 - 1)))
    print("    -> ROW-WISE: each diagonal is normalised to the off-diagonal in")
    print("       ITS OWN row.  0.87 % / 0.51 % -- inside the rounding of")
    print("       3-decimal printed values.  The single-reference reading is")
    print("       27 % off, so the ambiguity in the sentence is settled by the")
    print("       numbers, not by grammar.")
    print("\n    consequence [arithmetic, inputs L1] -- transcription-free form:")
    print("      margin = mu*_ss mu*_pp/(mu*_sp mu*_ps) = (mu*_ss/mu*_sp)(mu*_pp/mu*_ps)")
    print("             = 2.23 x 2.48 = %.4f      (from the printed numbers: %.4f)"
          % (2.23 * 2.48, ss * pp / (sp * ps)))
    print("    -> MgB2's Coulomb criterion margin IS the product of the two")
    print("       printed ratios.  It passes for the trivial reason that each")
    print("       band's own intraband Coulomb element exceeds the interband")
    print("       one IN THAT ROW -- no cancellation, no DOS bookkeeping.")
    print("\n    implied symmetric V^c (from mu*_ij = c N_j V_ij):")
    print("      V_ss/V_sp = 2.23 N_p/N_s = %.4f ,  V_pp/V_sp = 2.48 N_s/N_p = %.4f"
          % (2.23 * N_P / N_S, 2.48 * N_S / N_P))
    print("      => V_ss/V_pp = %.3f : the sigma (in-plane bonding, more"
          % ((2.23 * N_P / N_S) / (2.48 * N_S / N_P)))
    print("         localised) states carry the larger Coulomb matrix element,")
    print("         and the PRINTED mu*_ss > mu*_pp survives despite N_s < N_p.")
    print("      => V_ss V_pp / V_sp^2 = %.3f  = the same margin, as the"
          % ((2.23 * N_P / N_S) * (2.48 * N_S / N_P)))
    print("         convention-invariance theorem requires (N's cancel).")
    print("\n    NOT CLAIMED: that this material confirms Eq.(12).  It cannot --")
    print("    Golubov never applies Eq.(12); the mu* matrix is an ansatz with")
    print("    one fitted scale.  Eq.(12) is tested here only as algebra, in")
    print("    sections (A)-(E).")


if __name__ == "__main__":
    sec_A(); sec_B(); sec_C(); sec_D(); sec_E(); sec_F()
