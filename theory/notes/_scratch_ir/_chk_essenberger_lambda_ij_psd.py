#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_essenberger_lambda_ij_psd.py -- 2026-09-18, main session, from zero.

WHAT
----
First material-side evaluation of the convention-independent multiband
criterion

        lam_11 * lam_22  >=  lam_12 * lam_21          (equivalently det >= 0)

on a PUBLISHED band-resolved lambda_ij matrix.

L1 SOURCE (main session read the page image itself; the four numbers are
TYPESET LABELS inside the figure, not curve readings)
-------------------------------------------------------------------------
Essenberger, Sanna, Linscheid, Tandetzky, Profeta, Cudazzo, Gross,
"Superconducting pairing mediated by spin fluctuations from first
principles", arXiv:1411.2121v1 (= PRB 94, 014503 (2016)), Fig. 3(a),
electron-phonon Eliashberg function alpha^2 F_ij for FeSe, two bands {+,-}:

        lam_++ = 0.028      (positive peak)
        lam_-- = 0.029      (positive peak)
        lam_+- = 0.042      (curve dips NEGATIVE, valley ~ -0.2)
        lam_-+ = 0.048      (curve dips NEGATIVE)

Definition, PRB 94,014503 footnote [86], verbatim (via source-hunt, L1):
        lam_ij = 2 * INT |A_ij(omega)| / omega  d omega
        A_ij = alpha^2 F_ij (phonon) or Lambda^SF_ij (spin fluctuation)
        lam  = (1/N) SUM_ij lam_ij N_i,   N = N_+ + N_-
NOTE the absolute value is on the INTEGRAND.  arXiv v1 carries no definition
at all (its footnote [69] is a bare citation to Suhl-Matthias-Walker); the
definition is new in the published version.

DOS CONVENTION -- the criterion does not depend on it (CORRECTED 2026-09-18).
    (a) WRONG, first version of this file: "the scalar collapse
        lam = (1/N) SUM_ij lam_ij N_i is sensible ONLY if lam_ij = N_j V_ij."
        False.  (1/N) SUM_ij lam_ij N_i = SUM_i (N_i/N) lam_i with
        lam_i = SUM_j lam_ij -- a DOS-weighted average of per-band TOTAL
        couplings.  That is sensible under ANY convention for lam_ij, so the
        collapse formula pins nothing.
    (b) RIGHT, and stronger: the criterion's VERDICT is invariant under
        lam_ij = f(i) g(j) V_ij for ANY positive f, g, with V symmetric:
              lam_11 lam_22 = f1 g1 V11 * f2 g2 V22 = (f1 f2 g1 g2) V11 V22
              lam_12 lam_21 = f1 g2 V12 * f2 g1 V21 = (f1 f2 g1 g2) V12 V21
        The prefactors are IDENTICAL, so
              lam_12 lam_21 / (lam_11 lam_22) = V12 V21 / (V11 V22)
        exactly, for every one of: no DOS factor (f=g=1), lam_ij=N_j V_ij
        (f=1,g=N), lam_ij=N_i V_ij (f=N,g=1), lam_ij=sqrt(N_i N_j) V_ij.
        => whether the published table carries a DOS factor, and which index
        it sits on, CANNOT change the 2.48x.  No retrieval of arXiv:1409.7968
        is needed -- for this reason, not the one in (a).
    (c) Decoration, not load-bearing: lam_+- = 0.042 != lam_-+ = 0.048 shows
        lam_ij is NOT symmetric.  Since |g_ij|^2 = |g_ji|^2 by hermiticity
        makes V symmetric, some index-asymmetric factor must be present; so
        "no DOS factor inside lam_ij AND V symmetric" is contradicted by the
        data.  (Main session has NOT read footnote [86] itself; source-hunt
        paraphrases it as asserting no DOS factor.  The tension is flagged,
        not resolved -- and by (b) it does not matter.)

DIRECTION AUDIT (CLAUDE.md one-sided rule)
------------------------------------------
The published numbers are |.|-integrated, hence >= |signed integral| for every
entry, with EQUALITY iff that A_ij has a definite sign on the support.
So BOTH sides of the criterion are inflated and the induced error does NOT
have a determined direction.  Consequently:
  * the FAIL computed below is a fact about the PUBLISHED (|.|) matrix;
  * it does NOT by itself establish that the SIGNED matrix fails;
  * it DOES refute "the criterion is automatically satisfied", which is the
    only thing it is used for here.
Level: numerical-evidence on L1-transcribed inputs.
"""

import math

# ----- L1: Essenberger Fig. 3(a), FeSe, phonon channel -----------------
LAM = {("+", "+"): 0.028,
       ("-", "-"): 0.029,
       ("+", "-"): 0.042,
       ("-", "+"): 0.048}


def report(name, l11, l22, l12, l21):
    lhs = l11 * l22
    rhs = l12 * l21
    det = lhs - rhs
    tr = l11 + l22
    disc = tr * tr - 4.0 * det
    r1 = (tr + math.sqrt(disc)) / 2.0
    r2 = (tr - math.sqrt(disc)) / 2.0
    print("\n  %s" % name)
    print("    matrix  [[%.3f, %.3f], [%.3f, %.3f]]" % (l11, l12, l21, l22))
    print("    lam11*lam22 = %.6f" % lhs)
    print("    lam12*lam21 = %.6f" % rhs)
    print("    det         = %+.6f   -> criterion %s" %
          (det, "HOLDS" if det >= 0 else "FAILS"))
    if det < 0:
        print("    violation factor rhs/lhs = %.4f x" % (rhs / lhs))
    print("    eigenvalues = %+.5f , %+.5f  -> matrix is %s" %
          (r1, r2, "indefinite" if r1 * r2 < 0 else "semidefinite"))
    return det


def main():
    print("=" * 70)
    print("Multiband criterion lam11*lam22 >= lam12*lam21 on published lam_ij")
    print("FeSe, electron-phonon channel, Essenberger et al. Fig. 3(a)")
    print("=" * 70)

    det = report("published |.|-integrated values (L1 figure labels)",
                 LAM[("+", "+")], LAM[("-", "-")],
                 LAM[("+", "-")], LAM[("-", "+")])

    # --- what would it take to rescue the criterion? ----------------------
    print("\n  [rescue analysis] how much would the off-diagonals have to")
    print("  shrink under sign cancellation for the criterion to be restored?")
    lhs = LAM[("+", "+")] * LAM[("-", "-")]
    rhs = LAM[("+", "-")] * LAM[("-", "+")]
    f = math.sqrt(lhs / rhs)
    print("    need |signed lam_offdiag| <= %.4f x the published value" % f)
    print("    i.e. lam_+- <= %.4f and lam_-+ <= %.4f"
          % (f * LAM[("+", "-")], f * LAM[("-", "+")]))
    print("    (and that is only sufficient if the DIAGONALS do not shrink")
    print("     at all -- Fig. 3a shows the diagonal traces are the positive")
    print("     peaks, so their sign cancellation is plausibly the smaller")
    print("     effect, but this is NOT established here.)")

    # --- the definite-sign branch ----------------------------------------
    print("\n  [definite-sign branch] Fig. 3(a) shows the two DIAGONAL traces")
    print("  as the positive peaks and the two OFF-DIAGONAL traces as the")
    print("  negative valleys.  IF each A_ij has a definite sign on its")
    print("  support, then |.| is a PURE SIGN FLIP with no magnitude")
    print("  inflation, and the signed matrix is exactly:")
    report("signed matrix under the definite-sign premise",
           +LAM[("+", "+")], +LAM[("-", "-")],
           -LAM[("+", "-")], -LAM[("-", "+")])
    print("    -> identical: the two minus signs multiply out.  Under this")
    print("       premise the FAIL is EXACT, not merely an artefact of |.|.")
    print("    !! premise read from the SAME figure; not independently checked.")

    # --- convention invariance, demonstrated numerically -----------------
    print("\n  [convention invariance] lam_ij = f(i) g(j) V_ij, V symmetric.")
    print("  Recover V from the published lam under several conventions and")
    print("  confirm the criterion's verdict is identical in all of them.")
    l11, l22 = LAM[("+", "+")], LAM[("-", "-")]
    l12, l21 = LAM[("+", "-")], LAM[("-", "+")]
    base = (l12 * l21) / (l11 * l22)
    # pick a deliberately lopsided DOS pair to make any leakage visible
    N1, N2 = 0.37, 1.93
    for name, f, g in [
            ("no DOS factor        (f=1,        g=1)", (1.0, 1.0), (1.0, 1.0)),
            ("lam_ij = N_j V_ij    (f=1,        g=N)", (1.0, 1.0), (N1, N2)),
            ("lam_ij = N_i V_ij    (f=N,        g=1)", (N1, N2), (1.0, 1.0)),
            ("lam_ij = sqrt(NiNj)V (f=sqrt N,   g=sqrt N)",
             (N1 ** .5, N2 ** .5), (N1 ** .5, N2 ** .5))]:
        V11 = l11 / (f[0] * g[0])
        V22 = l22 / (f[1] * g[1])
        V12 = l12 / (f[0] * g[1])
        V21 = l21 / (f[1] * g[0])
        ratio = (V12 * V21) / (V11 * V22)
        print("    %-44s V12V21/(V11V22) = %.10f  sym(V)=%s"
              % (name, ratio, "yes" if abs(V12 - V21) < 1e-12 else "NO"))
    print("    published lam ratio                          = %.10f" % base)
    print("    -> all identical to the published ratio: the DOS convention")
    print("       cannot change the verdict.  (V is only symmetric in the")
    print("       convention that actually holds; the RATIO is invariant")
    print("       regardless, which is what the criterion consumes.)")

    # --- the structural point --------------------------------------------
    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    print("""
  [numerical-evidence, L1 inputs] The published FeSe electron-phonon lam_ij
  matrix is INDEFINITE: det = %+.6f, i.e. lam12*lam21 exceeds lam11*lam22
  by %.2fx.  The criterion FAILS on it.

  [what this does establish] It refutes "a first-principles PHONON lam_ij
  matrix automatically satisfies the criterion".  That claim was made --
  by me -- earlier today on a false Bochner/Gram analogy with mu_ij, and
  retracted on an analytic counterexample (pure interband coupling
  g = [[0,g12],[g12,0]] gives lam11*lam22 = 0 < lam12*lam21).  FeSe is the
  physical realisation of exactly that regime: interband-dominated coupling.
  The retraction is now confirmed from the materials side, independently.

  [what this does NOT establish] That the SIGNED lam_ij matrix fails.  The
  published entries are |.|-integrated; both sides of the criterion are
  inflated and the net direction is undetermined (see DIRECTION AUDIT).
  Settling it needs alpha^2F_ij(omega) itself, not the integrated table.

  [caveat on scope] Which repo lemma consumes this criterion has NOT been
  re-derived here.  Do not propagate this to a ceiling verdict without
  first pinning the consumption point.
""" % (det, rhs / lhs))


if __name__ == "__main__":
    main()
