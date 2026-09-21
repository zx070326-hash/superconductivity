#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_chk_cac6_lambda_ij_criterion.py -- 2026-09-18, main session, from zero.

Independent re-verification of source-hunt's claim that the convention-
independent multiband criterion

        lam_ii * lam_jj  >=  lam_ij * lam_ji

fails on the CaC6 band-resolved electron-phonon lambda_ij tables.

ALL NUMBERS BELOW ARE READ BY THE MAIN SESSION FROM THE PDF **TEXT LAYER** of
    A. Sanna, S. Pittalis, J. K. Dewhurst, M. Monni, S. Sharma,
    G. Ummarino, S. Massidda, E. K. U. Gross,
    Phys. Rev. B 85, 184514 (2012),  Table I, p.184514-3
fetched independently from harvest.aps.org to
    theory/notes/_scratch_ir/_sanna_prb85_184514.pdf   (8 pp, LaTeX text PDF)
NOT copied from source-hunt's note.  Two entries differ from its transcription
(flagged below).

DEFINITIONS, verbatim from the same text layer [L1]
---------------------------------------------------
Eq. (11), p.184514-2:
    alpha^2 F_JJ'(omega) = (1/N_J) SUM_{k in J} SUM_{k' in J'} SUM_nu
                           delta(xi_k) delta(xi_k') delta(omega - omega_{q,nu})
                           |g^nu_{k,k'}|^2
    "where N_J is the density of states per spin at the Jth FS portion."

Eq. (15), p.184514-3:
    lambda_JJ' = 2 INT alpha^2 F_JJ'(omega) / omega  d omega
    -> this settles source-hunt's B.7 open item ("table values are 2 int
       a2F/w dw: INFERRED"): it is now L1, not inferred.

Table I caption, verbatim:
    "Electron-phonon coupling lambda_JJ' and DOS for the three-band FS (top)
     and six-band FS (bottom) divisions.  DOS is given in states/eV/spin.
     The first index runs over the columns.  Isotropic DOS and lambda are
     0.787 and 0.870, respectively."

CONVENTION, derived from Eq. (11) [analytic]
--------------------------------------------
    <|g|^2>_JJ' := [SUM_{k in J,k' in J'} dd |g|^2] / (N_J N_J')   is SYMMETRIC
    => alpha^2 F_JJ' = (1/N_J) * N_J N_J' <|g|^2> = N_J' <|g|^2>
    => lam_JJ' = N_J' V_JJ' with V symmetric          (DOS on the RIGHT index)
    => detailed balance   N_J lam_JJ' = N_J' lam_J'J.

*** WHERE I DIVERGE FROM source-hunt ***
It claims the detailed-balance test "proves which index the DOS weight sits
on" (its B.3).  It does NOT, on its own: the test only pins the COMBINATION
(row/column reading) x (convention).  Reading A (rows = first index) with
lam = N_j V, and Reading B (columns = first index, i.e. what the caption
literally says) with lam = N_i V, give numerically IDENTICAL detailed-balance
residuals -- they are the same array transposed with the convention also
transposed.  Both are checked below and both pass.
What actually pins the convention is Eq. (11) alone (the 1/N_J sits on the
first index).  Given that, detailed balance then pins the LAYOUT to Reading A
-- i.e. the printed caption's "first index runs over the columns" is not the
reading that makes the table self-consistent.
The criterion is transpose-invariant either way, so no verdict depends on it.
"""

import itertools

# ---------------------------------------------------------------- Table I ---
# three-band, as printed (rows 1,2,3 ; columns 1,2,3)
T3 = [[0.301, 0.136, 0.257],
      [0.546, 0.239, 0.479],
      [0.427, 0.198, 0.367]]
N3 = [0.412, 0.104, 0.249]
L3 = ["1", "2", "3"]

# six-band, as printed
# !! source-hunt transcribed (1b,2a) as 0.055 and flagged it as a suspected
#    misread; the text layer prints 0.105 -- its guess was right.
# !! source-hunt transcribed (3a,2a) as 0.148; the text layer prints 0.145.
#    This second discrepancy was NOT flagged by it.
T6 = [[0.163, 0.126, 0.099, 0.033, 0.201, 0.046],
      [0.179, 0.140, 0.105, 0.035, 0.221, 0.050],
      [0.331, 0.245, 0.151, 0.084, 0.384, 0.096],
      [0.271, 0.202, 0.206, 0.047, 0.400, 0.080],
      [0.252, 0.194, 0.145, 0.061, 0.309, 0.073],
      [0.206, 0.157, 0.128, 0.044, 0.259, 0.060]]
N6 = [0.250, 0.176, 0.075, 0.031, 0.200, 0.056]
L6 = ["1a", "1b", "2a", "2b", "3a", "3b"]


def detailed_balance(T, N, lab, tag):
    """N_i lam_ij = N_j lam_ji, under Reading A (rows = first index)."""
    print("\n  [detailed balance] %s   (Reading A: rows = first index)" % tag)
    worst = 0.0
    worst_pair = None
    for i, j in itertools.combinations(range(len(N)), 2):
        a = N[i] * T[i][j]
        b = N[j] * T[j][i]
        rel = abs(a - b) / max(abs(a), abs(b))
        if rel > worst:
            worst, worst_pair = rel, (lab[i], lab[j])
        flag = "  <== FAIL" if rel > 0.05 else ""
        print("    %-3s-%-3s  N_i l_ij = %.6f   N_j l_ji = %.6f   rel = %6.2f%%%s"
              % (lab[i], lab[j], a, b, 100 * rel, flag))
    print("    worst = %.2f%% at %s" % (100 * worst, worst_pair))
    return worst


def detailed_balance_B(T, N, lab, tag):
    """Reading B (columns = first index, i.e. the printed caption) combined
    with the OPPOSITE convention lam_ij = N_i V_ij.  Should be identical."""
    worst = 0.0
    for i, j in itertools.combinations(range(len(N)), 2):
        # lam_ij = T[j][i] under Reading B ; test N_j lam_ij = N_i lam_ji
        a = N[j] * T[j][i]
        b = N[i] * T[i][j]
        rel = abs(a - b) / max(abs(a), abs(b))
        worst = max(worst, rel)
    print("    [Reading B + lam=N_i V] worst = %.2f%%  -> %s"
          % (100 * worst, "identical residuals, test cannot separate them"))
    return worst


def criterion(T, N, lab, tag):
    print("\n  [criterion  lam_ii lam_jj >= lam_ij lam_ji]  %s" % tag)
    n = len(N)
    viol = []
    for i, j in itertools.combinations(range(n), 2):
        lhs = T[i][i] * T[j][j]
        rhs = T[i][j] * T[j][i]
        r = rhs / lhs
        if r > 1.0:
            viol.append((r, lab[i], lab[j], lhs, rhs))
        print("    %-3s-%-3s  lhs = %.7f   rhs = %.7f   rhs/lhs = %7.4f  %s"
              % (lab[i], lab[j], lhs, rhs, r, "VIOLATION" if r > 1 else "ok"))
    npair = n * (n - 1) // 2
    viol.sort(reverse=True)
    print("    -> %d / %d pairs violate" % (len(viol), npair))
    for r, a, b, lhs, rhs in viol[:4]:
        print("       worst: %s-%s  %.4f x  (%.7f vs %.7f)" % (a, b, r, lhs, rhs))
    return viol, npair


def main():
    print("=" * 74)
    print("CaC6 lambda_ij (Sanna et al., PRB 85, 184514 (2012), Table I)")
    print("main-session transcription from the PDF TEXT LAYER")
    print("=" * 74)

    for T, N, lab, tag in ((T3, N3, L3, "three-band"), (T6, N6, L6, "six-band")):
        detailed_balance(T, N, lab, tag)
        detailed_balance_B(T, N, lab, tag)
        criterion(T, N, lab, tag)

    # --- are the two tables independent data points? ----------------------
    print("\n" + "=" * 74)
    print("  [independence check] is the three-band table a coarsening of the")
    print("  six-band one?  Compare summed DOS:")
    for k, (a, b) in enumerate([(0, 1), (2, 3), (4, 5)]):
        s = N6[a] + N6[b]
        print("    N_%s + N_%s = %.3f   vs  three-band N_%d = %.3f   (%+.1f%%)"
              % (L6[a], L6[b], s, k + 1, N3[k], 100 * (s / N3[k] - 1)))
    print("  -> same FS partitioned twice, DOS agree to ~3%: the two tables are")
    print("     NOT independent measurements.  Treat CaC6 as ONE material-side")
    print("     data point with two partitions, not two.")

    # --- isotropic cross-check --------------------------------------------
    print("\n  [isotropic cross-check, caption: N=0.787, lam=0.870]")
    for T, N, lab, tag in ((T3, N3, L3, "three-band"), (T6, N6, L6, "six-band")):
        Ntot = sum(N)
        lam = sum(N[i] * T[i][j] for i in range(len(N)) for j in range(len(N))) / Ntot
        print("    %-10s  sum N = %.3f (vs 0.787)   (1/N) sum_ij N_i lam_ij = %.4f"
              " (vs 0.870)" % (tag, Ntot, lam))
    print("  -> this is the scalar collapse lam = (1/N) sum_ij lam_ij N_i with")
    print("     the DOS on the LEFT index in the sum; agreement at the few-%%")
    print("     level is a transcription check on BOTH the table and the DOS row.")
    print("     NB this does NOT pin the convention (see my earlier correction):")
    print("     it equals sum_i (N_i/N) lam_i under any convention.  It is a")
    print("     TRANSCRIPTION check, nothing more.")

    separability()


def separability():
    """[analytic] What the criterion really tests, and why violations are
    neither surprising nor pathological.

    Write g^nu_{k,k'} for the e-ph matrix element of mode nu.  Two facts:

    (S1) SEPARABLE, SINGLE MODE  g_{k,k'} = gamma_k conj(gamma_k')
         => |g|^2 = |gamma_k|^2 |gamma_k'|^2
         => SUM_{k in J,k' in J'} dd |g|^2 = A_J A_J' with A_J := SUM_{k in J} d |gamma_k|^2
         => lam_JJ' ∝ A_J A_J' / N_J
         => lam_ii lam_jj = (A_i^2/N_i)(A_j^2/N_j) = (A_i A_j/N_i)(A_j A_i/N_j)
                          = lam_ij lam_ji                       EXACTLY.
         The criterion is SATURATED -- for ANY gamma and ANY partition.

    (S2) SEPARABLE, MANY MODES   g^nu_{k,k'} = gamma^nu_k conj(gamma^nu_k')
         => V = SUM_nu (2/omega_nu) a^nu (a^nu)^T ,  a^nu_J := A^nu_J / N_J
         a sum of rank-one terms with POSITIVE weights => V is PSD
         => the criterion HOLDS (with slack, generically strict).

    => Violating the criterion requires g^nu_{k,k'} to be NON-separable for
       individual modes.  And that is exactly what momentum conservation
       delivers: g^nu_{k,k'} depends on the transfer q = k - k', i.e. it is a
       TOEPLITZ kernel, not a product.  |g(q)|^2 >= 0 pointwise in q does NOT
       make the kernel |g(k-k')|^2 positive semidefinite in (k,k') -- that is
       Bochner's condition on the real-space transform, a different statement.

    This is the SAME Hadamard-vs-Gram distinction that killed the earlier
    "phonon lam_ij is automatically PSD" claim, now visible in the momentum
    structure rather than in the band structure.  CaC6 is not a pathology; it
    is the generic consequence of g depending on k - k'.
    """
    import random
    print("\n" + "=" * 74)
    print("  [separability, numerical check of (S1) and (S2)]")
    random.seed(20260918)
    nk, nb = 240, 4
    part = [random.randrange(nb) for _ in range(nk)]

    def blocks(absg2):
        """absg2[k][kp] -> lam_JJ' = (1/N_J) sum_{k in J,k' in J'} |g|^2."""
        N = [part.count(J) for J in range(nb)]
        lam = [[0.0] * nb for _ in range(nb)]
        for k in range(nk):
            for kp in range(nk):
                lam[part[k]][part[kp]] += absg2[k][kp]
        for J in range(nb):
            for Jp in range(nb):
                lam[J][Jp] /= N[J]
        return lam

    def worst_ratio(lam):
        w = 0.0
        for i, j in itertools.combinations(range(nb), 2):
            w = max(w, (lam[i][j] * lam[j][i]) / (lam[i][i] * lam[j][j]))
        return w

    # (S1) single separable mode
    gam = [random.uniform(0.1, 2.0) for _ in range(nk)]
    a2 = [[(gam[k] * gam[kp]) ** 2 for kp in range(nk)] for k in range(nk)]
    print("    (S1) single separable mode      worst rhs/lhs = %.12f  (expect 1)"
          % worst_ratio(blocks(a2)))

    # (S2) five separable modes, positive weights
    a2 = [[0.0] * nk for _ in range(nk)]
    for _ in range(5):
        gam = [random.uniform(0.1, 2.0) for _ in range(nk)]
        w = random.uniform(0.5, 3.0)
        for k in range(nk):
            for kp in range(nk):
                a2[k][kp] += w * (gam[k] * gam[kp]) ** 2
    print("    (S2) 5 separable modes          worst rhs/lhs = %.6f  (expect <=1)"
          % worst_ratio(blocks(a2)))

    # (S3) Toeplitz |g(k-k')|^2, i.e. momentum conservation -- generic case
    nv = 6
    kv = [2.0 * 3.141592653589793 * k / nk for k in range(nk)]
    import math
    a2 = [[0.0] * nk for _ in range(nk)]
    for _ in range(nv):
        m = random.randrange(1, 5)
        ph = random.uniform(0, 6.283185307179586)
        w = random.uniform(0.5, 3.0)
        for k in range(nk):
            for kp in range(nk):
                g = math.cos(m * (kv[k] - kv[kp]) + ph)
                a2[k][kp] += w * g * g          # |g(k-k')|^2 >= 0 pointwise
    print("    (S3) Toeplitz |g(k-k')|^2       worst rhs/lhs = %.6f  (>1 = violates)"
          % worst_ratio(blocks(a2)))
    print("    -> nonnegativity of |g|^2 in the TRANSFER variable does not give")
    print("       positive semidefiniteness of the kernel in (k,k').  Same")
    print("       Hadamard-vs-Gram distinction, seen in momentum space.")


if __name__ == "__main__":
    main()
