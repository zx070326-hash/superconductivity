r"""
P4 . N -> infinity INTERCHANGE UNIFORMITY -- computer-assisted certificate for
the three pillars of `P4_interchange_uniformity.md`.

v2 (2026-09-10): rewritten after the adversarial review
`theory/notes/2026-09-10_interchange_review.md`.  The v1 certificate was
rigorous but only AT 7 wbar NODES; v2 closes it over the whole interval [0.3,5]
and removes the one float64-dependent branch.  Changes, all reviewer-driven:
  (1) PILLAR 3 REWRITTEN via Poincare separation + Cauchy interlacing.  The Weyl
      route (N_0, rho_N, and a float64 eigvalsh direct check for N < N_0 -- audit
      point R2) is REDUNDANT: interlacing gives inf_N gap_N >= mu_2 - alpha
      outright, with no N_0, no float64 eigenvalue, and a bound 1.3-1.9x larger.
  (2) STAGE 3B ADDED (audit point R1): the node certificate is extended to the
      interval [0.3,5] on a 246-point geometric grid with segment envelopes.
      A Lipschitz argument on the ORIGINAL 7-node spacing provably does NOT
      close (it goes negative near wbar~0.36), so refinement was mandatory.
  (3) `cert_sup_on_J` de-regressed: its Lipschitz constant used to be a grid
      maximum of another certified quantity; it is now a closed-form segment
      Schur bound (finding 5).
  (4) `hs_offhead` weakened to the form its own derivation actually proves
      (finding 4b), and mu/alpha padded for float64 entry rounding (finding 7).
  (5) Stage 5 re-scoped in wording: it is a finite-N consistency check at
      N=400, not a direct verification of the true-operator bound (finding 6).

================================================================================
WHICH OPERATOR IS BEING CERTIFIED (read this first -- anti-confusion note)
================================================================================
This file certifies the KAY **II** operator at a FIXED frequency parameter wbar,

    H(w) = -H1(w) + H2(w) + H3(w)   on  l^2(N_0),     [[k]](w) := w^2/(k^2+w^2),
      (H1)_{nn} = (1/(2n+1)) sum_{k=1}^n 2[[k]](w)                (diagonal)
      (H2)_{nm} = (1-delta_{nm}) [[|n-m|]](w) / sqrt((2n+1)(2m+1))
      (H3)_{nm} = [[n+m+1]](w)     / sqrt((2n+1)(2m+1)),

entry formulas taken VERBATIM from the verified `p4_remainder_interchange.py`
(`H_matrix`, T5, commit 520443d).  It is NOT the KAY **I** gamma-model operator
G(gamma) whose top eigenvalue frak_g(2)=1.3181405 was enclosed in P4-B and whose
lambda_2(G) <= 2e-2 was certified in P4-C.  The two are linked only through the
small-w BRIDGE  lim_{w->0} rho(H^(N)(w))/w^2 = frak_g^(N)(2)  (P4-A), which is
re-checked here in Stage 0 as an ANCHOR ASSERTION ONLY.  No P4-B/P4-C numerical
value ever enters the derivation of any bound below (anti-circularity).

================================================================================
WHAT IS CERTIFIED
================================================================================
Because the entry formulas do not depend on the truncation size,
  H^(N)(w) = P_N H(w) P_N          (leading principal block / compression)
so gap_N := lam_1(H^(N)) - lam_2(H^(N)) is the gap of a compression of the
compact self-adjoint H(w).  The capstone needs inf_N gap_N > 0.

PILLAR 1 (Stage 1).  A certified upper bound on
    rho_N(w) := 2 ||(I-P_N) H(w)||_op   >=  ||H - P_N H P_N||_op.
Two independent rigorous routes, both implemented; the certificate uses the min:
  (a) HS route, as sketched in the derivation note: H1 by its max diagonal
      element (closed form C1(w) = pi w coth(pi w) - 1), H2+H3 by a
      Hilbert-Schmidt row tail with a zeta/integral bound.  Decays like N^-1/2.
  (b) SCHUR/Holmgren route (added here): ||A|| <= sqrt(rowsum_max * colsum_max).
      Decays like N^-1; this is the route that carries every number, with (a)
      kept as an independent cross-check on the same quantity.
  Since the Pillar-3 rewrite, rho_N is DEMOTED: it no longer establishes "the
  gap does not collapse" (interlacing does that), and serves only the corollary
  (7) convergence rate and the Davis-Kahan step giving kappa_N -> kappa.

PILLAR 2 (Stage 2).  A certified positive lower bound delta(w) <= lam_1(H) -
lam_2(H) for the LIMIT operator:
  lam_1 >= mu   : Rayleigh quotient of a finitely-supported trial vector,
                  evaluated in mpmath (exact: finite support => <v,Hv> uses only
                  the head block).  Rigorous with no truncation error.
  lam_2 <= alpha: Cholesky-deflation certificate, same machine as P4-C
                  (`p4_gap_conjecture2.py` / `p4_gamma_general.py`
                  `alpha_certificate`), transplanted from G(gamma) to H(w):
                  alpha = lam_2(head A) + ||B|| + ||D||, where the head bound
                  comes from an mpmath dps=30 Cholesky test on  c*I - P A P
                  (P = deflation off the head top eigenvector) and ||B||, ||D||
                  come from the Pillar-1 Schur envelopes.  Every link of this
                  chain was re-derived by the reviewer and CONFIRMED.
  delta := mu - alpha.

PILLAR 3 (Stage 3).  inf_{N>=2} gap_N(w) >= g_alt(w) := mu_2(w) - alpha(w),
by Poincare separation (lam_2(A_N) <= lam_2(A) <= alpha for every N) plus Cauchy
interlacing (lam_1(A_N) nondecreasing in N, so >= lam_1(A_2) >= mu_2), with
mu_2 the EXACT closed-form top eigenvalue of the 2x2 block.  No N_0, no rho_N,
no float64 eigenvalue enters this bound.

PILLAR 3B (Stage 3B).  The same bound over the whole INTERVAL [0.3,5]: mu_2 is
bounded below per cell EXACTLY (a fixed trial vector makes the Rayleigh quotient
a positive combination of the increasing functions [[k]](w)), ||B||,||D|| by
one segment-envelope evaluation per cell, and lam_2(head) by endpoint Cholesky
plus a closed-form Lipschitz slack.

PILLAR 4 (Stage 4).  Compression monotonicity ||H^(N)(j)|| <= ||H^(N')(j)||
(N<N') verified, and REGRESS-FREE certified upper bounds B_2^inf, B_4^inf on the
LIMIT operator derivative norms, sup over the support-image interval J.

ASSEMBLY (Stage 5).  K_*(w) = w^4 B_4^inf eta_max^2/24
+ (1/2) w^4 (B_2^inf)^2 / g_alt, then a finite-N consistency check at N=400:
non-vacuity, the sigma^4 rate, and N-independence of the measured coefficient.

NUMERICAL DISCIPLINE (mandated): numpy float64 `eigvalsh` only (never generic
`eig`, never mpmath eigensolvers); mpmath dps=30 restricted to the Cholesky
certificate, the Rayleigh quotient and the closed-form tail integrals; all loops
have fixed bounds, no open-ended search.

USAGE: `python p4_interchange_uniformity.py` runs everything; pass stage names
(0 1 2 3 3B 4 5) to run a subset.
"""

import sys
import time

import numpy as np
from mpmath import mp

mp.dps = 30

PASS = [0]
_STAGE0 = [0]


def _check(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name} {detail}")
    PASS[0] += 1


def _stage_summary(tag, note=""):
    k = PASS[0] - _STAGE0[0]
    _STAGE0[0] = PASS[0]
    print(f"  [{tag}] {k} checks passed. {note}")


# ============================================================================
# STAGE 0 -- builders.  `_brk_der` / `H_matrix` / `G_matrix` are copied VERBATIM
# from the verified `p4_remainder_interchange.py` (T5) and are NOT re-derived.
# A vectorized twin `H_mat` is added for speed and checked to agree exactly.
# ============================================================================

def _brk_der(k, w, order):
    """d^order/dw^order of [[k]](w)=w^2/(k^2+w^2).  VERBATIM from T5."""
    k2 = float(k) ** 2
    w = float(w)
    if k == 0:
        return 1.0 if order == 0 else 0.0
    d = k2 + w * w
    if order == 0:
        return w * w / d
    if order == 1:
        return 2 * k2 * w / d**2
    if order == 2:
        return 2 * k2 * (k2 - 3 * w * w) / d**3
    if order == 3:
        return 24 * k2 * w * (w * w - k2) / d**4
    if order == 4:
        return 24 * k2 * (-5 * w**4 + 10 * k2 * w * w - k2 * k2) / d**5
    if order == 5:
        return 240 * k2 * w * (3 * w**4 - 10 * k2 * w * w + 3 * k2 * k2) / d**6
    raise ValueError(order)


def H_matrix(N, w, order=0):
    """KAY II H^(N)(w) or its elementwise w-derivative.  VERBATIM from T5."""
    M = np.zeros((N, N), dtype=np.float64)
    for n in range(N):
        for m in range(N):
            nrm = np.sqrt((2 * n + 1) * (2 * m + 1))
            val = _brk_der(n + m + 1, w, order) / nrm
            if n != m:
                val += _brk_der(abs(n - m), w, order) / nrm
            else:
                s = sum(2 * _brk_der(k, w, order) for k in range(1, n + 1))
                val -= s / (2 * n + 1)
            M[n, m] = val
    return M


def G_matrix(N, gamma):
    """KAY I G^(N)(gamma).  VERBATIM from T5 (anchor assertions only)."""
    M = np.zeros((N, N), dtype=np.float64)
    for n in range(N):
        for m in range(N):
            nrm = np.sqrt((2 * n + 1) * (2 * m + 1))
            val = 1.0 / (nrm * (n + m + 1) ** gamma)
            if n != m:
                val += 1.0 / (nrm * abs(n - m) ** gamma)
            else:
                s = sum(2.0 / k**gamma for k in range(1, n + 1))
                val -= s / (2 * n + 1)
            M[n, m] = val
    return M


def _brk_vec(K, w, order):
    """Vectorized `_brk_der` over an integer array K (same closed forms)."""
    K = np.asarray(K, dtype=np.float64)
    k2 = K * K
    w = float(w)
    d = k2 + w * w
    if order == 0:
        out = w * w / d
    elif order == 1:
        out = 2 * k2 * w / d**2
    elif order == 2:
        out = 2 * k2 * (k2 - 3 * w * w) / d**3
    elif order == 3:
        out = 24 * k2 * w * (w * w - k2) / d**4
    elif order == 4:
        out = 24 * k2 * (-5 * w**4 + 10 * k2 * w * w - k2 * k2) / d**5
    elif order == 5:
        out = 240 * k2 * w * (3 * w**4 - 10 * k2 * w * w + 3 * k2 * k2) / d**6
    else:
        raise ValueError(order)
    zero = (K == 0)
    if np.any(zero):
        out = np.where(zero, 1.0 if order == 0 else 0.0, out)
    return out


def H_mat(N, w, order=0):
    """Vectorized twin of `H_matrix` (checked bit-comparable in Stage 0)."""
    n = np.arange(N)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)
    A = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    off = _brk_vec(A, w, order)
    np.fill_diagonal(off, 0.0)
    M = (off + _brk_vec(S, w, order)) * nrm
    cum = np.concatenate(([0.0], np.cumsum(2.0 * _brk_vec(np.arange(1, N), w, order))))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M


def gap_of(A):
    """(lam1, lam2, gap) of a symmetric matrix via float64 eigvalsh."""
    E = np.linalg.eigvalsh(A)
    if E.size == 1:
        return float(E[-1]), -np.inf, np.inf
    return float(E[-1]), float(E[-2]), float(E[-1] - E[-2])


def opnorm(A):
    return float(np.max(np.abs(np.linalg.eigvalsh(A))))


def run_stage0():
    print("=== STAGE 0: builders (verbatim T5 port) + anchor assertions ===")
    t0 = time.time()
    # (a) vectorized twin == verbatim builder, all orders used downstream
    for N in (5, 17, 40):
        for w in (0.3, 1.0, 2.0, 5.0):
            for order in (0, 1, 2, 3, 4, 5):
                d = np.max(np.abs(H_mat(N, w, order) - H_matrix(N, w, order)))
                sc = max(1.0, np.max(np.abs(H_matrix(N, w, order))))
                _check("H_mat==H_matrix", d <= 1e-12 * sc, f"N={N} w={w} o={order} d={d:.3e}")
    # (b) N=1 closed form
    for w in (0.3, 1.0, 2.0, 5.0):
        _check("H^(1)=w^2/(w^2+1)", abs(H_mat(1, w)[0, 0] - w**2 / (w**2 + 1)) < 1e-13)
    # (c) ANCHOR ONLY (never used in any bound): KAY I fingerprint + small-w bridge
    g2 = float(np.max(np.linalg.eigvalsh(G_matrix(2, 2))))
    _check("KAY I g^(2)(2) fingerprint", abs(g2 - (20 + np.sqrt(13819)) / 108) < 1e-12)
    for N in (1, 2, 3, 4):
        gN = float(np.max(np.linalg.eigvalsh(G_matrix(N, 2))))
        ratio = float(np.max(np.linalg.eigvalsh(H_mat(N, 1e-3)))) / 1e-6
        _check("bridge rho(H^(N)(w))/w^2 -> g^(N)(2)", abs(ratio - gN) < 1e-4,
               f"N={N} {ratio} vs {gN}")
    gG = float(np.max(np.linalg.eigvalsh(G_matrix(300, 2))))
    _check("P4-B anchor g(2) in [1.31714,1.31820] (ASSERTION ONLY)",
           1.31714 <= gG <= 1.31820, f"g^(300)(2)={gG:.7f}")
    print(f"  bridge/anchor: g^(300)(2) = {gG:.7f}  (P4-B enclosure [1.31714,1.31820];"
          " used as a self-test gate ONLY, never inside a bound)")
    _stage_summary("Stage 0", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE 1 -- PILLAR 1: certified tail rate rho_N(w).
#
# Envelope machinery.  Every block of H (and of every w-derivative of H) has
# entries of the form  psi(k) / sqrt((2n+1)(2m+1))  with k = |n-m| >= 1 (H2),
# k = n+m+1 >= 1 (H3), or the H1 diagonal (1/(2n+1)) sum_{k=1}^n 2 psi(k), where
# psi is POSITIVE and DECREASING on [1/2, inf).  For order 0, psi(t)=[[t]](w).
# For the derivatives we use the closed envelopes (each verified numerically in
# Stage 1a against the exact `_brk_der`):
#     |[[k]]'|   <= 2 w / D,        |[[k]]''|  <= 6 / D,
#     |[[k]]'''| <= 24 w / D^2,     |[[k]]''''|<= 240 / D^2,
#     |[[k]]^(5)| <= 2400 w / D^3,          D := k^2 + w^2.
# (Proofs: e.g. |[[k]]''| = 2k^2|k^2-3w^2|/D^3 <= 2k^2 * 3D / D^3 <= 6/D since
#  |k^2-3w^2| <= 3D and k^2 <= D.  Same one-line argument for the others, using
#  5w^4+10k^2w^2+k^4 <= 10 D^2 and 3w^4+10k^2w^2+3k^4 <= 10 D^2.)
# ============================================================================

# psi(t) = PSI_C[order](w) / D^PSI_P[order],  D = t^2 + w^2
PSI_P = {0: 1, 1: 1, 2: 1, 3: 2, 4: 2, 5: 3}


def _psi_c(order, w):
    w = float(w)
    return {0: w * w, 1: 2 * w, 2: 6.0, 3: 24 * w, 4: 240.0, 5: 2400 * w}[order]


def psi(t, w, order=0):
    """Positive decreasing envelope for |d^order [[t]](w)| (see block comment)."""
    return _psi_c(order, w) / (float(t) ** 2 + float(w) ** 2) ** PSI_P[order]


# ---------------------------------------------------------------------------
# SEGMENT specs.  Every envelope below is driven by a triple
#     spec = (C, A, p)   meaning   psi(t) = C / (t^2 + A^2)^p,
# so the SAME machinery serves two uses:
#   pointwise   spec_point(w, j)   = (c_j(w), w, p_j)      -- psi at a single w;
#   segment     spec_seg(a, b, j)  = (c_j(b), a, p_j)      -- a valid upper
#       envelope for |d^j [[t]](w)| SIMULTANEOUSLY for every w in [a,b], because
#       c_j(w) in {w^2, 2w, 6, 24w, 240, 2400w} is nondecreasing in w (so
#       c_j(w) <= c_j(b)) and (t^2+w^2)^-p is increasing in w (so it is
#       <= (t^2+a^2)^-p).  This is the reviewer's `L_seg` construction
#       (`_scratch_ir/e8_lipclose.py`) and it has NO regress: the bound is a
#       closed form in (a,b), not a grid maximum of another certified quantity.
# CAUTION (why the naive "evaluate at the right endpoint" shortcut is WRONG for
# j>=1): psi_j(t;w) is increasing in w ONLY for j=0.  For j=2,4 it DECREASES in
# w, and for j=1,3,5 it peaks at an interior w (e.g. psi_1 = 2w/(t^2+w^2) peaks
# at w=t).  So sup_{[a,b]} psi_j = psi_j(.;b) holds for order 0 only, and the
# mixed numerator-at-b / denominator-at-a form above is what covers all orders.
# ---------------------------------------------------------------------------

def spec_point(w, order=0):
    return (mp.mpf(_psi_c(order, w)), mp.mpf(float(w)), PSI_P[order])


def spec_seg(a, b, order=0):
    return (mp.mpf(_psi_c(order, float(b))), mp.mpf(float(a)), PSI_P[order])


_ENV_CACHE = {}


def envelope(spec):
    """(S1, IntTail) for psi(t)=C/(t^2+A^2)^p: S1 = sum_{k>=1} psi(k) (upper
    bound), IntTail(x) = int_x^inf psi.

    S1 is a partial sum to k=2000 plus the integral tail (psi decreasing =>
    sum_{k>K} psi(k) <= int_K^inf psi), so S1 is a rigorous UPPER bound on the
    series.  Everything is evaluated in mpmath (closed forms below, elementary).
    """
    C, A, p = spec
    key = (str(C), str(A), p)
    if key in _ENV_CACHE:
        return _ENV_CACHE[key]

    def int_tail(x):
        x = mp.mpf(x)
        if p == 1:                      # int_x^inf dt/(t^2+A^2)
            return C * (mp.pi / 2 - mp.atan(x / A)) / A
        if p == 2:                      # int_x^inf dt/(t^2+A^2)^2
            return C * ((mp.pi / 2 - mp.atan(x / A)) / (2 * A**3)
                        - x / (2 * A**2 * (x**2 + A**2)))
        if p == 3:                      # int_x^inf dt/(t^2+A^2)^3
            return C * (3 * (mp.pi / 2 - mp.atan(x / A)) / (8 * A**5)
                        - x / (4 * A**2 * (x**2 + A**2) ** 2)
                        - 3 * x / (8 * A**4 * (x**2 + A**2)))
        raise ValueError(p)

    K = 2000
    ks = np.arange(1, K + 1, dtype=np.float64)
    Cf, Af = float(C), float(A)
    S1 = mp.mpf(float(np.sum(Cf / (ks**2 + Af**2) ** p))) + int_tail(K)
    _ENV_CACHE[key] = (S1, int_tail)
    return S1, int_tail


def psi_spec(t, spec):
    C, A, p = spec
    return C / (mp.mpf(float(t)) ** 2 + A**2) ** p


def row_env(n, spec):
    """R_bar(n): a MANIFESTLY DECREASING upper envelope for the full row sum
    sum_{m>=0} |H^{(order)}_{nm}|, valid for every n >= 1.

    Derivation (h := ceil(n/2); all steps elementary):
      H1 diagonal          : (1/(2n+1)) sum_{k=1}^n 2 psi(k)  <= 2 S1/(2n+1)
      H2 row, m >= h       : 1/sqrt(2m+1) <= 1/sqrt(n+1), sum_{m!=n} psi(|n-m|)
                             <= 2 S1        ==> 2 S1 / sqrt((n+1)(2n+1))
      H2 row, m <  h       : |n-m| >= n-h+1 > n/2 so psi(|n-m|) <= psi(n/2), and
                             sum_{m<h} (2m+1)^-1/2 <= sqrt(2h-1) <= sqrt(n);
                             sqrt(n)/sqrt(2n+1) <= 1/sqrt(2)  ==> psi(n/2)/sqrt(2)
      H3 row, m >= h       : sum_{m>=0} psi(n+m+1) = sum_{k>n} psi(k)
                             <= IntTail(n)  ==> IntTail(n)/sqrt((n+1)(2n+1))
      H3 row, m <  h       : psi(n+m+1) <= psi(n)   ==> psi(n)/sqrt(2)
    Each of the five terms is decreasing in n (psi and IntTail are decreasing),
    hence max_{n>=N} R_bar(n) = R_bar(N) with no monotonicity hand-waving.
    """
    S1, int_tail = envelope(spec)
    n = mp.mpf(n)
    r = mp.sqrt((n + 1) * (2 * n + 1))
    return (2 * S1 / (2 * n + 1)
            + 2 * S1 / r
            + psi_spec(float(n) / 2.0, spec) / mp.sqrt(2)
            + int_tail(n) / r
            + psi_spec(float(n), spec) / mp.sqrt(2))


def row_env0(spec):
    """Row n=0 (excluded from `row_env`, whose split assumes n>=1): the H1
    diagonal vanishes and both off-diagonal halves are bounded by S1, since
    1/sqrt((2*0+1)(2m+1)) <= 1 and sum_{m>=1} psi(m) = sum_{m>=0} psi(m+1) = S1."""
    S1, _ = envelope(spec)
    return 2 * S1


def col_env(m, N, spec):
    """C_bar(m) for m < N: upper bound on sum_{n>=N} |H^{(order)}_{nm}|.

      sum_{n>=N} psi(n-m)   = sum_{k>=N-m} psi(k) <= IntTail(N-m-1)   (N-m>=1)
      sum_{n>=N} psi(n+m+1) = sum_{k>=N+m+1} psi(k) <= IntTail(N+m)
      1/sqrt((2n+1)(2m+1)) <= 1/sqrt((2N+1)(2m+1))
    """
    _, int_tail = envelope(spec)
    return (int_tail(N - m - 1) + int_tail(N + m)) / mp.sqrt((2 * m + 1) * (2 * N + 1))


def schur_offhead(N, w, order=0, spec=None):
    """||(I-P_N) H^{(order)}||_op bound by the Schur/Holmgren test
        ||A|| <= sqrt( max_n sum_m |A_nm| * max_m sum_n |A_nm| ).
    Rows n >= N: sum over all m <= R_bar(n) <= R_bar(N).
    Columns: m >= N gives sum_{n>=N} <= full row sum at m <= R_bar(N);
             m <  N is the finite exact max of C_bar(m).
    `spec` overrides (w, order) -- pass spec_seg(a,b,order) for a bound valid on
    the whole segment [a,b]."""
    if spec is None:
        spec = spec_point(w, order)
    rmax = row_env(N, spec)
    cmax = rmax
    for m in range(N):
        c = col_env(m, N, spec)
        if c > cmax:
            cmax = c
    return mp.sqrt(rmax * cmax)


def schur_tailblock(N, w, order=0, spec=None):
    """||H^{(order)}|_{n,m >= N}||_op <= R_bar(N) (Schur; restricting the column
    range only shrinks both Schur maxima)."""
    if spec is None:
        spec = spec_point(w, order)
    return row_env(N, spec)


def schur_coupling(Na, w, order=0, spec=None):
    """||B||_op for the head-tail coupling B = (H_nm)_{n<Na, m>=Na}.
    row sums (over m>=Na, n<Na) are exactly C_bar(n); column sums (over n<Na,
    m>=Na) are bounded by the full row sum at m, i.e. R_bar(m) <= R_bar(Na)."""
    if spec is None:
        spec = spec_point(w, order)
    cmax = mp.mpf(0)
    for n in range(Na):
        c = col_env(n, Na, spec)
        if c > cmax:
            cmax = c
    return mp.sqrt(cmax * row_env(Na, spec))


def schur_full_seg(a, b, order, Mrows=200, Mcols=8000):
    """RIGOROUS sup_{w in [a,b]} ||H^{(order)}(w)||_op, with NO regress.

    The symmetric envelope matrix E (E_nm >= |H^{(order)}_nm(w)| entrywise for
    every w in [a,b], built from spec_seg) dominates H^{(order)} in the sense
    ||H^{(order)}|| <= ||E||, and for a symmetric nonnegative E the Schur test
    gives ||E|| <= max_n (row sum of E).  Rows n < Mrows are summed explicitly
    (exact for m < Mcols, closed-form integral tail beyond); rows n >= Mrows use
    the decreasing envelope R_bar(Mrows).  Row 0 uses `row_env0`.

    This is the replacement for the "L = grid max of the next certified norm"
    slack used earlier: that L was itself only a grid maximum, so the Lipschitz
    argument regressed (reviewer finding 5).  Here L is a closed form in (a,b).
    """
    spec = spec_seg(a, b, order)
    C, A, p = float(spec[0]), float(spec[1]), spec[2]
    _, int_tail = envelope(spec)

    m = np.arange(Mcols, dtype=np.float64)
    invm = 1.0 / np.sqrt(2.0 * m + 1.0)
    kk = np.arange(Mrows + Mcols + 2, dtype=np.float64)
    psi_k = C / (kk**2 + A**2) ** p
    cum = np.concatenate(([0.0], np.cumsum(psi_k[1:])))     # cum[n] = sum_{k=1}^n psi(k)
    inv_sqrt_last = 1.0 / np.sqrt(2.0 * Mcols + 1.0)

    best = float(row_env0(spec))
    for n in range(1, Mrows):
        nrm_n = 1.0 / np.sqrt(2.0 * n + 1.0)
        r = 2.0 * cum[n] / (2.0 * n + 1.0)                   # H1 diagonal
        d = np.abs(m - n)
        h2 = psi_k[d.astype(np.int64)] * invm
        h2[n] = 0.0                                          # no diagonal in H2
        h3 = psi_k[(n + m + 1).astype(np.int64)] * invm
        tail = (float(int_tail(Mcols - n - 1)) + float(int_tail(Mcols + n))) * inv_sqrt_last
        r += nrm_n * (float(h2.sum()) + float(h3.sum()) + tail)
        best = max(best, r)
    return max(best, float(row_env(Mrows, spec)))


def hs_offhead(N, w, M=200000):
    """||(I-P_N) H(w)||_op via the HS route of the derivation note (order 0 only).

      ||(I-P_N)H|| <= max_{n>=N}(H1)_nn + ||(I-P_N)(H2+H3)||_HS
                   <= 2 S1/(2N+1)      + sqrt(S(N)),
      S(N) = sum_{n>=N} sum_m ((H2+H3)_nm)^2 <= sum_{n>=N} 2 (T(n)+U(n))/(2n+1),
      using (a+b)^2 <= 2a^2+2b^2 and, with Z2 := 2 sum_{k>=1}[[k]]^2,
        T(n) = Z2/(n+1) + [[n/2]]^2 (3/2 + ln(n)/2)        (H2 inner sum)
        U(n) = [[n]]^2 (3/2 + ln(n)/2) + Z2/(2(n+1))       (H3 inner sum)
      (same m>=ceil(n/2) vs m<ceil(n/2) split as `row_env`, with
       sum_{m<h} 1/(2m+1) <= 3/2 + ln(n)/2).
    [The U(n) tail term is Z2/(2(n+1)), NOT the tighter Z2/(2(2n+1)) an earlier
    revision used: the m>=h half only yields 1/(2m+1) <= 1/(n+1), so 1/(2(n+1))
    is what the stated split actually proves.  Reviewer finding (2026-09-10,
    `_scratch_ir/e3_hs.py`): the tighter form held numerically (ratio <=0.34) but
    was NOT derivable from the comment.  We weaken the CODE to the provable form
    rather than only the comment, so every number this function returns is
    derived.  It is never load-bearing anyway (see `rho_bound`).]
    The n-sum is exact for N <= n < M and closed-form beyond, using
    [[n/2]]^2+[[n]]^2 <= 17 w^4/n^4 and int_{M-1}^inf (3/2+ln x/2)/x^5 dx; the
    Z2 part of the beyond-M tail is sum_{n>=M} 3 Z2/((n+1)(2n+1)) <= 3Z2/(2(M-1)).
    Decays like N^-1/2 -- strictly weaker than `schur_offhead`; kept as an
    INDEPENDENT cross-check of the same quantity."""
    wf = float(w)
    ks = np.arange(1, 20001, dtype=np.float64)
    brk = wf**2 / (ks**2 + wf**2)
    Z2 = 2.0 * (np.sum(brk**2) + wf**4 / (3.0 * 20000.0**3))  # + integral tail
    n = np.arange(N, M, dtype=np.float64)
    b_half = wf**2 / (n**2 / 4.0 + wf**2)
    b_n = wf**2 / (n**2 + wf**2)
    lg = 1.5 + 0.5 * np.log(n)
    T = Z2 / (n + 1.0) + b_half**2 * lg
    U = b_n**2 * lg + Z2 / (2.0 * (n + 1.0))
    S = float(np.sum(2.0 * (T + U) / (2.0 * n + 1.0)))
    Mm = float(M - 1)
    S += 1.5 * Z2 / Mm + 17.0 * wf**4 * ((1.5 + 0.5 * np.log(Mm)) / (4 * Mm**4)
                                         + 0.5 / (16 * Mm**4))
    S1, _ = envelope(spec_point(w, 0))
    return float(2 * S1 / (2 * N + 1)) + np.sqrt(S)


def rho_bound(N, w):
    """rho_N(w) = 2||(I-P_N)H(w)|| certified upper bound (min of both routes)."""
    return float(2.0 * min(float(schur_offhead(N, w)), hs_offhead(N, w)))


WBARS = [0.3, 0.5, 1.0 / np.sqrt(3.0), 1.0, 2.0, 3.0, 5.0]


def run_stage1():
    print("=== STAGE 1: PILLAR 1 -- certified tail rate rho_N(wbar) ===")
    t0 = time.time()

    # 1a. the psi envelopes really dominate the exact derivatives
    for w in WBARS:
        for order in range(6):
            ks = np.arange(1, 400)
            ex = np.abs(_brk_vec(ks, w, order))
            en = np.array([psi(k, w, order) for k in ks])
            _check("psi envelope dominates |[[k]]^(order)|", np.all(ex <= en * (1 + 1e-12)),
                   f"w={w} order={order} worst={np.max(ex/en):.4f}")
        # S1 is a partial sum + integral tail, hence an UPPER bound on the exact
        # series C1(w)/2 = sum_{k>=1}[[k]]; it must dominate it and be tight.
        S1, _ = envelope(spec_point(w, 0))
        C1 = np.pi * w / np.tanh(np.pi * w) - 1.0
        _check("2*S1 >= C1(w)=pi w coth(pi w)-1 and tight",
               C1 <= float(2 * S1) <= C1 + 1e-5 * max(1.0, C1),
               f"w={w} {float(2*S1)} vs {C1}")

    # 1b. SOUNDNESS: the row/col envelopes must dominate the true finite sums
    Mbig = 1500
    for w in (0.3, 1.0, 5.0):
        Hb = np.abs(H_mat(Mbig, w))
        rows = Hb.sum(axis=1)
        for N in (10, 40, 120):
            _check("row_env dominates true row sums", float(row_env(N, spec_point(w, 0))) >= rows[N] - 1e-12,
                   f"w={w} N={N}: env={float(row_env(N, spec_point(w, 0))):.6f} true={rows[N]:.6f}")
            cols = Hb[N:, :].sum(axis=0)
            for m in (0, N // 3, N - 1):
                _check("col_env dominates true col sums",
                       float(col_env(m, N, spec_point(w, 0))) >= cols[m] - 1e-12,
                       f"w={w} N={N} m={m}")
            # the certified operator-norm bound must dominate a large finite
            # numerical evaluation of ||(I-P_N) H P_Mbig||
            blk = H_mat(Mbig, w)[N:, :]
            true_op = float(np.sqrt(np.max(np.linalg.eigvalsh(blk @ blk.T))))
            for nm, val in (("schur", float(schur_offhead(N, w))), ("hs", hs_offhead(N, w))):
                _check(f"{nm} bound >= true ||(I-P_N)H||", val >= true_op,
                       f"w={w} N={N}: {nm}={val:.6f} true>={true_op:.6f}")

    # 1c. SOUNDNESS at the DERIVATIVE orders.  1b covers order 0 only, but
    # orders 1-5 feed cert_opnorm_inf -> B_2^inf/B_4^inf and the segment
    # Lipschitz constants.  (Reviewer补测, `_scratch_ir/e4_env.py`.)
    for w in (0.3, 1.0, 5.0):
        for order in (1, 2, 3, 4, 5):
            Hb = np.abs(H_mat(400, w, order))
            rows = Hb.sum(axis=1)
            sp = spec_point(w, order)
            for N in (2, 10, 40, 120):
                _check("row_env dominates true row sums (derivative orders)",
                       float(row_env(N, sp)) >= rows[N],
                       f"w={w} order={order} N={N}")
                cols = Hb[N:, :].sum(axis=0)
                _check("col_env dominates true col sums (derivative orders)",
                       float(col_env(0, N, sp)) >= cols[0],
                       f"w={w} order={order} N={N}")
                blk = H_mat(400, w, order)[N:, :]
                true_op = float(np.sqrt(np.max(np.linalg.eigvalsh(blk @ blk.T))))
                _check("schur_offhead >= true ||(I-P_N)H^(j)||",
                       float(schur_offhead(N, w, order)) >= true_op,
                       f"w={w} order={order} N={N}")

    # 1d. SOUNDNESS of the SEGMENT envelope and of schur_full_seg: the segment
    # bound must dominate the pointwise truth at every w inside the segment,
    # for every order (this is what licenses the R1 closure and cert_sup_on_J).
    for a, b in ((0.30, 0.31), (0.95, 1.05), (2.0, 2.5), (4.5, 5.0)):
        for order in (0, 1, 2, 3, 4, 5):
            sp = spec_seg(a, b, order)
            Lseg = schur_full_seg(a, b, order, Mrows=60, Mcols=3000)
            for w in np.linspace(a, b, 5):
                ks = np.arange(1, 300)
                _check("segment psi envelope dominates |[[k]]^(j)(w)| on [a,b]",
                       np.all(np.abs(_brk_vec(ks, w, order))
                              <= np.array([float(psi_spec(k, sp)) for k in ks]) * (1 + 1e-12)),
                       f"[{a},{b}] order={order} w={w}")
                _check("schur_full_seg >= true ||H^(j)(w)|| on [a,b]",
                       Lseg >= opnorm(H_mat(220, w, order)),
                       f"[{a},{b}] order={order} w={w}: {Lseg:.5f}")

    print(f"  {'wbar':>7} " + "".join(f"{('N='+str(N)):>11}" for N in (10, 20, 40, 80, 160)))
    for w in WBARS:
        print(f"  {w:>7.4f} " + "".join(f"{rho_bound(N, w):>11.5f}" for N in (10, 20, 40, 80, 160)))
    print("  (rho_N = certified bound on 2||(I-P_N)H||; Schur route ~1/N dominates)")
    print(f"  {'wbar':>7}  {'schur(N=40)':>12}  {'HS(N=40)':>12}   (both rigorous; min used)")
    for w in WBARS:
        print(f"  {w:>7.4f}  {2*float(schur_offhead(40, w)):>12.5f}  {2*hs_offhead(40, w):>12.5f}")
    _stage_summary("Stage 1", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE 2 -- PILLAR 2: certified limit spectral gap delta(w) = mu - alpha.
# ============================================================================

NA_HEAD = 70      # head size for the Cholesky-deflation certificate
NV_RAY = 140      # support size of the Rayleigh trial vector

# Rounding pad (reviewer finding, precision hygiene).  `rayleigh_mu` /
# `alpha_certificate_H` do EXACT mpmath arithmetic, but on matrix entries that
# were produced in float64 (elementwise rounding <= 1.1e-16, crude operator-norm
# gap ||A_f - A|| <= 8e-15 at Na=70/Nv=140).  We therefore deflate mu and
# inflate alpha by 1e-12 -- three orders of magnitude above that residual and
# nine below the smallest certified quantity (delta >= 0.1), so it costs nothing
# and makes every returned number a bound on the EXACT operator, not on its
# float64 image.
ROUND_PAD = mp.mpf("1e-12")


def rayleigh_mu(w, Nv=NV_RAY):
    """Certified lower bound mu <= lam_1(H(w)).

    v is the float64 top eigenvector of H^(Nv)(w), zero-padded to l^2; since v
    has FINITE SUPPORT, <v,H v> involves only the head block and is computed
    EXACTLY (mpmath dps=30) -- no truncation error, no eigenvalue solver enters
    the bound (the eigenvector is merely a trial vector; any vector would give a
    valid, just weaker, bound).  Minus ROUND_PAD for float64 entry rounding."""
    A = H_mat(Nv, w)
    _, V = np.linalg.eigh(A)
    v = V[:, -1]
    vm = mp.matrix([mp.mpf(float(x)) for x in v])
    Am = mp.matrix(Nv, Nv)
    for i in range(Nv):
        for j in range(Nv):
            Am[i, j] = mp.mpf(float(A[i, j]))
    Av = Am * vm
    num = sum(vm[i] * Av[i] for i in range(Nv))
    den = sum(vm[i] ** 2 for i in range(Nv))
    return num / den - ROUND_PAD


def _brk_mp(k, w):
    return w * w / (mp.mpf(k) ** 2 + w * w)


def lam1_H2_exact(w):
    """lam_1(H^(2)(w)) in CLOSED FORM (mpmath), no eigensolver:
        H^(2) = [[ [[1]],  ([[1]]+[[2]])/sqrt3 ],
                 [ .    ,  ([[3]] - 2[[1]])/3  ]]
    (n=m=0: H1 diagonal is empty, H3 gives [[1]]; n=0,m=1: H2 gives [[1]] and H3
    gives [[2]], both over sqrt(1*3); n=m=1: H3 gives [[3]]/3, H1 gives -2[[1]]/3.)
    lam_1 = (tr + sqrt(tr^2 - 4 det))/2."""
    w = mp.mpf(w)
    b1, b2, b3 = _brk_mp(1, w), _brk_mp(2, w), _brk_mp(3, w)
    h00 = b1
    h01 = (b1 + b2) / mp.sqrt(3)
    h11 = (b3 - 2 * b1) / 3
    tr = h00 + h11
    det = h00 * h11 - h01 * h01
    return (tr + mp.sqrt(tr * tr - 4 * det)) / 2, (h00, h01, h11)


def mu2_cell_lb(a, b):
    """RIGOROUS lower bound on lam_1(H^(2)(w)) valid for EVERY w in [a,b].

    Fix v = top eigenvector of H^(2)(a) (any fixed vector gives a valid bound).
    Then Q(w) := <v,H^(2)(w)v>/||v||^2 <= lam_1(H^(2)(w)) for all w, and

        Q(w) = k1 [[1]](w) + k2 [[2]](w) + k3 [[3]](w),
        k1 = (v0^2 + 2 v0 v1/sqrt3 - 2 v1^2/3)/||v||^2,
        k2 = (2 v0 v1/sqrt3)/||v||^2,   k3 = (v1^2/3)/||v||^2.

    Every [[k]](w) = w^2/(k^2+w^2) is strictly INCREASING in w, so if
    k1,k2,k3 >= 0 then Q is increasing and  Q(w) >= Q(a) = lam_1(H^(2)(a))
    for all w in [a,b] -- an EXACT segment bound with zero Lipschitz slack.
    The sign of k1 is asserted, not assumed (it is ~1.34 at small w and ~1.44
    as w -> inf; k2,k3 >= 0 automatically since v0,v1 > 0).
    Returns (bound, k1, k2, k3)."""
    lam, (h00, h01, h11) = lam1_H2_exact(a)
    v0, v1 = h01, lam - h00          # eigenvector of a 2x2; both > 0 when h01>0
    nn = v0 * v0 + v1 * v1
    k1 = (v0 * v0 + 2 * v0 * v1 / mp.sqrt(3) - 2 * v1 * v1 / 3) / nn
    k2 = (2 * v0 * v1 / mp.sqrt(3)) / nn
    k3 = (v1 * v1 / 3) / nn
    return lam - ROUND_PAD, k1, k2, k3


def alpha_certificate_H(w, Na=NA_HEAD, seg=None):
    """Certified upper bound alpha >= lam_2(H(w)).

    Structure copied from P4-C's `alpha_certificate` (`p4_gamma_general.py`),
    with G(gamma) replaced by H(w) and the gamma-specific tail bounds replaced by
    the Stage-1 Schur envelopes:

      H = [[A, B],[B^T, D]],  A = head (Na x Na), D = tail block.
      lam_2(H) <= lam_2(diag(A,D)) + ||B||          (Weyl, off-diagonal part)
               <= max(lam_2(A), lam_1(D)) + ||B||
               <= lam_2(A) + ||D|| + ||B||          (both terms >= 0 below)
      lam_2(A) <= c   is CERTIFIED by an mpmath dps=30 Cholesky test on
                      c*I - P A P  (P = deflation off the head top eigenvector
                      u): if that matrix is SPD then lam_1(P A P) < c, and since
                      lam_2(A) <= lam_1(P A P) + |residual| we anchor c on the
                      float64 top eigenvalue of PAP plus a 1e-7 pad, exactly as
                      in P4-C (u is only float64-accurate, so PAP retains a
                      ~1e-15 residual in the deflated direction).
    Each link of this chain was re-derived independently by the reviewer
    (2026-09-10 report Sec. 3) and confirmed.

    `seg=(a,b)` makes ||B|| and ||D|| valid SIMULTANEOUSLY on [a,b] by using the
    segment envelope; the c part still refers to the single w passed in (the
    caller adds the Lipschitz slack, see `alpha_cell_ub`).
    Returns (alpha, chol_ok, lam2A_up, normB, normD)."""
    A_np = H_mat(Na, w)
    _, V = np.linalg.eigh(A_np)
    u = V[:, -1]
    Au_f = A_np @ u
    uAu_f = float(u @ Au_f)
    PAP_f = A_np - np.outer(u, Au_f) - np.outer(Au_f, u) + uAu_f * np.outer(u, u)
    c = mp.mpf(float(np.linalg.eigvalsh(PAP_f)[-1])) + mp.mpf("1e-7")

    Am = mp.matrix(Na, Na)
    for i in range(Na):
        for j in range(Na):
            Am[i, j] = mp.mpf(float(A_np[i, j]))
    um = mp.matrix([mp.mpf(float(x)) for x in u])
    Au = Am * um
    uAu = sum(um[i] * Au[i] for i in range(Na))
    M = mp.matrix(Na, Na)
    for i in range(Na):
        for j in range(Na):
            pap = Am[i, j] - um[i] * Au[j] - Au[i] * um[j] + um[i] * um[j] * uAu
            M[i, j] = (c if i == j else mp.mpf(0)) - pap
    try:
        mp.cholesky(M)
        chol_ok = True
    except Exception:  # noqa: BLE001
        chol_ok = False
    if not chol_ok:
        return None, False, None, None, None
    sp = spec_seg(seg[0], seg[1], 0) if seg is not None else spec_point(w, 0)
    nB = schur_coupling(Na, w, 0, spec=sp)
    nD = schur_tailblock(Na, w, 0, spec=sp)
    return c + nB + nD + ROUND_PAD, True, c, nB, nD


DELTA = {}
ALPHA = {}
MU2 = {}
GALT = {}


def run_stage2():
    print("=== STAGE 2: PILLAR 2 -- certified limit gap delta(wbar) = mu - alpha ===")
    t0 = time.time()
    print(f"  head Na={NA_HEAD} (Cholesky-deflation, mpmath dps={mp.dps}), "
          f"Rayleigh support Nv={NV_RAY}, rounding pad {mp.nstr(ROUND_PAD, 3)}")
    print(f"  {'wbar':>7}  {'mu (<=lam1)':>12}  {'lam2(A)<=c':>11}  {'||B||':>9}  "
          f"{'||D||':>9}  {'alpha>=lam2':>12}  {'delta':>10}")
    for w in WBARS:
        mu = rayleigh_mu(w)
        alpha, ok, c, nB, nD = alpha_certificate_H(w)
        _check("Cholesky deflation certificate SPD", ok, f"w={w}")
        delta = mu - alpha
        DELTA[w] = float(delta)
        ALPHA[w] = alpha
        _check("delta(wbar) > 0", delta > 0, f"w={w} delta={float(delta):.6f}")
        # soundness vs float64 truth at large N (assertion only, not used in bound)
        l1, l2, _ = gap_of(H_mat(400, w))
        _check("mu <= lam1 numeric", float(mu) <= l1 + 1e-9, f"w={w}")
        _check("alpha >= lam2 numeric", float(alpha) >= l2 - 1e-9, f"w={w}")
        _check("delta <= true gap numeric", float(delta) <= l1 - l2 + 1e-9, f"w={w}")
        print(f"  {w:>7.4f}  {float(mu):>12.7f}  {float(c):>11.3e}  {float(nB):>9.5f}  "
              f"{float(nD):>9.5f}  {float(alpha):>12.7f}  {float(delta):>10.6f}")
    _stage_summary("Stage 2", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE 3 -- PILLAR 3, REWRITTEN VIA POINCARE SEPARATION (reviewer, Sec. 2.2).
#
# The earlier version routed the uniform gap bound through Weyl:
#     gap_N >= delta - 2 rho_N,  N_0 := min{N : 2 rho_N < delta/2},
#     g_* := min( min_{2<=N<N_0} gap_N , delta/2 ),
# whose first branch was a float64 eigvalsh direct check (audit point R2).  That
# whole branch is REDUNDANT.  Two elementary facts about a compression
# A_N = P_N A P_N of a compact self-adjoint A do the job outright:
#
#   (i)  POINCARE SEPARATION:  lam_k(A_N) <= lam_k(A) for every N and k.
#        In particular lam_2(A_N) <= lam_2(A) <= alpha -- no rho_N, no Weyl.
#   (ii) CAUCHY INTERLACING on nested principal blocks: lam_1(A_N) is
#        nondecreasing in N, so lam_1(A_N) >= lam_1(A_2) for every N >= 2.
#
# Hence, with mu_2 <= lam_1(H^(2)) any lower bound on the 2x2 top eigenvalue,
#
#     inf_{N>=2} gap_N(w)  >=  mu_2(w) - alpha(w)  =:  g_alt(w),          (3')
#
# which is fully rigorous, has NO N_0, NO float64 eigenvalue, and no discussion
# of the zero-padded operator's lam_2 sign -- and is 1.3-1.9x LARGER than the old
# g_* = delta/2 at every node.  N=2 is the binding case precisely because the
# bound must hold for the SMALLEST N in the infimum; using a larger trial
# dimension Nv would only bound inf_{N>=Nv}, reopening the finite direct check.
#
# rho_N SURVIVES but is demoted: it no longer serves "the gap does not collapse",
# only the convergence RATE of corollary (7) and the Davis-Kahan sin-Theta step
# that gives kappa_N -> kappa.  N_0 is still reported for that purpose.
# ============================================================================

NMAX_SCAN = 400
GSTAR = {}
N0 = {}


def run_stage3():
    print("=== STAGE 3: PILLAR 3 -- inf_N gap_N >= mu_2 - alpha (Poincare route) ===")
    t0 = time.time()
    print("  (i) Poincare separation: lam_2(A_N) <= lam_2(A) <= alpha, every N.")
    print("  (ii) Cauchy interlacing: lam_1(A_N) nondecreasing => >= lam_1(A_2) >= mu_2.")
    print("  => inf_{N>=2} gap_N >= g_alt := mu_2 - alpha.  No N_0, no float64 gap.")
    print(f"  {'wbar':>7}  {'mu_2 (exact 2x2)':>17}  {'alpha':>10}  {'g_alt':>10}  "
          f"{'old g_*=delta/2':>15}  {'gain':>6}  {'N_0 (rate only)':>15}")
    tables = {}
    for w in WBARS:
        mu2, k1, k2, k3 = mu2_cell_lb(w, w)
        _check("mu_2 monotonicity coefficients k1,k2,k3 >= 0", k1 >= 0 and k2 >= 0 and k3 >= 0,
               f"w={w} k=({float(k1):.4f},{float(k2):.4f},{float(k3):.4f})")
        MU2[w] = mu2
        alpha = ALPHA[w]
        galt = float(mu2 - alpha)
        GALT[w] = galt
        GSTAR[w] = galt
        _check("g_alt > 0", galt > 0, f"w={w} g_alt={galt:.6f}")
        # mu_2 must really be lam_1 of the 2x2 (float64 assertion only)
        l1_2 = float(np.linalg.eigvalsh(H_mat(2, w))[-1])
        _check("mu_2 == lam_1(H^(2)) numeric", abs(float(mu2) - l1_2) < 1e-11, f"w={w}")
        # N_0 retained ONLY for the corollary-(7) rate / Davis-Kahan
        n0 = None
        for N in range(2, 2001):
            if 2.0 * rho_bound(N, w) < DELTA[w] / 2.0:
                n0 = N
                break
        _check("N_0 exists within scan budget", n0 is not None, f"w={w}")
        N0[w] = n0
        Hbig = H_mat(NMAX_SCAN, w)
        gaps = np.empty(NMAX_SCAN + 1)
        gaps[:2] = np.nan
        for N in range(2, NMAX_SCAN + 1):
            E = np.linalg.eigvalsh(Hbig[:N, :N])
            gaps[N] = E[-1] - E[-2]
        tables[w] = gaps
        gall = float(np.min(gaps[2:NMAX_SCAN + 1]))
        # EVIDENCE (not part of the certificate): the certified bound must of
        # course sit below every directly computed gap_N.
        _check("g_alt <= every directly computed gap_N (N=2..400)", gall >= galt - 1e-12,
               f"w={w} min={gall:.6f} g_alt={galt:.6f}")
        _check("g_alt beats the old Weyl g_*=delta/2", galt > DELTA[w] / 2.0, f"w={w}")
        print(f"  {w:>7.4f}  {float(mu2):>17.7f}  {float(alpha):>10.6f}  {galt:>10.6f}  "
              f"{DELTA[w]/2:>15.6f}  {galt/(DELTA[w]/2):>6.2f}x  {n0:>15}")

    print("\n  gap_N(wbar) vs N  [EVIDENCE ONLY -- the certificate above uses no "
          "float64 gap]:")
    Ns = [2, 4, 8, 16, 32, 64, 128, 256, 400]
    print(f"    {'wbar':>7} " + "".join(f"{('N='+str(N)):>10}" for N in Ns)
          + f"{'g_alt':>10}")
    for w in WBARS:
        print(f"    {w:>7.4f} " + "".join(f"{tables[w][N]:>10.5f}" for N in Ns)
              + f"{GALT[w]:>10.5f}")
    g1 = tables[1.0]
    _check("T5 cross-check gap_2(1)~1.137", abs(g1[2] - 1.137) < 5e-3, f"{g1[2]:.5f}")
    _check("T5 cross-check gap_32(1)~0.717", abs(g1[32] - 0.717) < 5e-3, f"{g1[32]:.5f}")
    print(f"  T5 cross-check at wbar=1: gap_2={g1[2]:.5f} (T5: 1.137), "
          f"gap_32={g1[32]:.5f} (T5: 0.717), gap_400={g1[400]:.5f} -> stable, no collapse.")
    _stage_summary("Stage 3", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE 3B -- R1 CLOSURE: from 7 nodes to the whole interval [0.3, 5].
#
# The node certificate says nothing between nodes, and the reviewer showed
# (`_scratch_ir/e8_lipclose.py`) that a Lipschitz argument on the 7-node spacing
# does NOT close: cell [0.3,0.5] has a NEGATIVE closure bound near w~0.36 even
# with the most optimistic Lipschitz constant.  So we refine the grid AND make
# every ingredient a genuine segment bound:
#
#   mu_2 over [a,b] : EXACT, no slack -- lam_1(H^(2)(.)) is bounded below by
#                     lam_1(H^(2)(a)) via the fixed-trial-vector monotonicity
#                     argument of `mu_2_cell_lb`.
#   ||B||, ||D||    : segment Schur envelope spec_seg(a,b,0) -- one evaluation,
#                     valid for all w in [a,b] at once.
#   lam_2(head A)   : the Cholesky certificate is NOT continuous in w (c is
#                     anchored on a discrete float64 eigenvalue), so instead of a
#                     Lipschitz bound on alpha as a whole we bound its three
#                     components separately, as dispatched: run the certificate
#                     at BOTH endpoints and add L*(b-a)/2, where
#                     L = schur_full_seg(a,b,1) is a rigorous, regress-free sup
#                     of ||H'|| on the cell (|lam_2(A(w)) - lam_2(A(p))| <=
#                     ||A(w)-A(p)|| <= L|w-p| by Weyl, and A' = P_Na H' P_Na is a
#                     compression of H' so ||A'|| <= ||H'||).
# ============================================================================

NGRID = 241       # GEOMETRIC grid points over [0.3, 5]


def wbar_grid(npts=NGRID):
    """Geometric (log-uniform) grid on [0.3, 5], with the 7 physical nodes forced
    in.  Geometric rather than uniform because the cell slack is L*h/2 and the
    rigorous L = sup||H'|| GROWS with wbar (1.9 at 0.3, 6.0 at 5): a uniform grid
    therefore wastes resolution at the left end and starves it at the right,
    which is exactly what capped the [3,5] retention at 0.73 in the uniform run.
    Constant h/wbar equalizes the relative slack across the interval."""
    pts = np.geomspace(WBARS[0], WBARS[-1], npts)
    pts = np.unique(np.concatenate([pts, np.array(WBARS)]))
    return pts


def alpha_cell_ub(a, b, cA, cB, L):
    """Upper bound on lam_2(H(w)) valid for every w in [a,b].

    sup_w min(cA + L(w-a), cB + L(b-w)) = (cA+cB)/2 + L(b-a)/2 <= max(cA,cB) +
    L(b-a)/2 for the head part; the envelope parts are segment bounds already."""
    nB = schur_coupling(NA_HEAD, b, 0, spec=spec_seg(a, b, 0))
    nD = schur_tailblock(NA_HEAD, b, 0, spec=spec_seg(a, b, 0))
    head = max(cA, cB) + mp.mpf(L) * (mp.mpf(b) - mp.mpf(a)) / 2
    return head + nB + nD + ROUND_PAD


G_INTERVAL = [None]


def run_stage3b():
    print("=== STAGE 3B: R1 CLOSURE -- inf over the INTERVAL [0.3, 5] ===")
    t0 = time.time()
    grid = wbar_grid()
    print(f"  grid: {len(grid)} geometric points on [0.3, 5] (h/wbar ~ "
          f"{grid[1]/grid[0]-1:.4f} constant); mu_2 exact per cell, ||B||,||D|| by "
          f"segment envelope, lam_2(A) by endpoint Cholesky + rigorous Lipschitz slack")

    # Cholesky c at every grid point (endpoints of all cells)
    cs = []
    for w in grid:
        _, ok, c, _, _ = alpha_certificate_H(float(w))
        _check("Cholesky SPD on the refined grid", ok, f"w={w}")
        cs.append(c)

    worst = None
    rows = []
    for i in range(len(grid) - 1):
        a, b = float(grid[i]), float(grid[i + 1])
        L = schur_full_seg(a, b, 1)
        mu2, k1, k2, k3 = mu2_cell_lb(a, b)
        _check("mu_2 segment monotonicity k1,k2,k3 >= 0",
               k1 >= 0 and k2 >= 0 and k3 >= 0, f"cell [{a},{b}]")
        al = alpha_cell_ub(a, b, cs[i], cs[i + 1], L)
        g = float(mu2 - al)
        _check("cell g_alt > 0", g > 0, f"cell [{a:.4f},{b:.4f}] g={g:.6f}")
        rows.append((a, b, L, float(mu2), float(al), g))
        if worst is None or g < worst[5]:
            worst = rows[-1]

    G_INTERVAL[0] = min(r[5] for r in rows)
    print(f"  {'segment':>18}  {'h':>7}  {'L_seg max':>9}  {'cell g_alt min':>14}  "
          f"{'nodal g_alt':>11}  {'retain':>7}")
    for i in range(len(WBARS) - 1):
        lo, hi = WBARS[i], WBARS[i + 1]
        sel = [r for r in rows if r[0] >= lo - 1e-9 and r[1] <= hi + 1e-9]
        gmin = min(r[5] for r in sel)
        nodal = min(GALT[lo], GALT[hi])
        print(f"  [{lo:6.4f},{hi:6.4f}]  {sel[0][1]-sel[0][0]:>7.4f}  "
              f"{max(r[2] for r in sel):>9.3f}  {gmin:>14.6f}  {nodal:>11.6f}  "
              f"{gmin/nodal:>7.3f}")
    print(f"  worst cell: [{worst[0]:.4f},{worst[1]:.4f}]  L_seg={worst[2]:.3f}  "
          f"mu_2={worst[3]:.6f}  alpha={worst[4]:.6f}  g_alt={worst[5]:.6f}")
    print(f"  ==> inf_{{wbar in [0.3,5]}} inf_{{N>=2}} gap_N(wbar) >= "
          f"{G_INTERVAL[0]:.6f}  (INTERVAL, not just the 7 nodes)")
    _check("interval infimum positive", G_INTERVAL[0] > 0)
    _check("interval infimum retains >=80% of the worst nodal g_alt",
           G_INTERVAL[0] >= 0.8 * min(GALT.values()),
           f"{G_INTERVAL[0]:.6f} vs {min(GALT.values()):.6f}")
    _stage_summary("Stage 3B", f"({time.time()-t0:.1f}s)")




# ============================================================================
# STAGE 4 -- PILLAR 4: compression monotonicity + certified B_2^inf, B_4^inf.
# ============================================================================

ETA_MAX = 2.0
SIGMA0 = 0.25
NB = 200          # head size for the limit-derivative-norm certificate


def cert_opnorm_inf(w, order, N=NB):
    """Certified upper bound on ||H^{(order)}(w)||_op for the LIMIT operator:
        ||H^(j)|| <= ||P_N H^(j) P_N|| + 2||(I-P_N) H^(j)||   (eq. (2) of the note)
    with the second term from the Stage-1 Schur envelope."""
    return opnorm(H_mat(N, w, order)) + 2.0 * float(schur_offhead(N, w, order))


def cert_sup_on_J(w, order, ncell=96, N=NB):
    """CERTIFIED sup over J = [wbar(1-sigma0 eta_max), wbar(1+sigma0 eta_max)] of
    ||H^{(order)}(.)||_op for the LIMIT operator -- REGRESS-FREE.

    Earlier version (reviewer finding 5, DOWNGRADE): sup_J was bounded by
    "grid max + L*h/2" with L itself only the grid MAXIMUM of the next certified
    norm.  That is not a sup, so the Lipschitz argument regressed onto another
    un-closed grid bound -- and closing it in that style at order 4 would have
    needed a 6th derivative, which `_brk_der` does not have.

    Now: split J into cells and on each cell [p,q] use
        sup_{[p,q]} ||H^(j)||  <=  max(||A(p)||, ||A(q)||) + L*(q-p)/2
                                   + 2 * ||(I-P_N) H^(j)||_seg,
    where A = P_N H^(j) P_N (float64 eigvalsh on the head, exact),
      L = schur_full_seg(p, q, j+1)  is a CLOSED-FORM rigorous sup of
          ||H^(j+1)|| on the cell (Schur on the segment envelope), and
          |d/dw ||A(w)|| | <= ||A'(w)|| = ||P_N H^(j+1) P_N|| <= ||H^(j+1)||,
      the tail term uses spec_seg(p,q,j), valid on the whole cell at once.
    Order 4 therefore needs only the order-5 envelope, which exists.  No grid
    maximum of a certified quantity is ever used as a Lipschitz constant.

    ncell is set high (96) because the rigorous Schur L is 2-7x larger than the
    eigvalsh values it replaces, and the slack L*h/2 shrinks like 1/ncell; the
    one-shot alternative sup_J <= schur_full_seg(a,b,j) with no cells at all is
    valid but far worse (e.g. 19.5 vs 2.7 for B_2 at wbar=0.3).  Head norms are
    evaluated once per cell EDGE and shared between adjacent cells."""
    a = max(w * (1.0 - SIGMA0 * ETA_MAX), 1e-6)
    b = w * (1.0 + SIGMA0 * ETA_MAX)
    edges = np.linspace(a, b, ncell + 1)
    heads = [opnorm(H_mat(N, float(x), order)) for x in edges]
    best = 0.0
    for i in range(ncell):
        p, q = float(edges[i]), float(edges[i + 1])
        L = schur_full_seg(p, q, order + 1)
        tail = 2.0 * float(schur_offhead(N, q, order, spec=spec_seg(p, q, order)))
        best = max(best, max(heads[i], heads[i + 1]) + L * (q - p) / 2.0 + tail)
    return best


BINF = {}


def run_stage4():
    print("=== STAGE 4: PILLAR 4 -- compression monotonicity + B_2^inf, B_4^inf ===")
    t0 = time.time()
    # (a) compression monotonicity of the derivative norms
    for w in WBARS:
        for order in (2, 4):
            prev = -np.inf
            for N in (2, 4, 8, 16, 32, 64, 128, 256):
                v = opnorm(H_mat(N, w, order))
                _check("||H^(N)(j)|| monotone in N (compression)", v >= prev - 1e-12,
                       f"w={w} order={order} N={N}: {v} < {prev}")
                prev = v
            _check("||H^(N)(j)|| <= certified limit bound",
                   prev <= cert_opnorm_inf(w, order) + 1e-12, f"w={w} order={order}")
    print(f"  J = [wbar*{1-SIGMA0*ETA_MAX:.2f}, wbar*{1+SIGMA0*ETA_MAX:.2f}] "
          f"(eta_max={ETA_MAX}, sigma0={SIGMA0}); head N={NB} + certified Schur tail")
    print(f"  {'wbar':>7}  {'|H(2)|_N=256':>13}  {'B2^inf':>10}  "
          f"{'|H(4)|_N=256':>15}  {'B4^inf':>10}")
    for w in WBARS:
        B2 = cert_sup_on_J(w, 2)
        B4 = cert_sup_on_J(w, 4)
        BINF[w] = (B2, B4)
        _check("B2^inf finite/positive", np.isfinite(B2) and B2 > 0, f"w={w}")
        _check("B4^inf finite/positive", np.isfinite(B4) and B4 > 0, f"w={w}")
        print(f"  {w:>7.4f}  {opnorm(H_mat(256, w, 2)):>13.5f}  {B2:>10.5f}  "
              f"{opnorm(H_mat(256, w, 4)):>15.4f}  {B4:>10.4f}")
    _stage_summary("Stage 4", f"({time.time()-t0:.1f}s)")


# ============================================================================
# STAGE 5 -- assembly of K_*(wbar) and INDEPENDENT test of the true-operator
# bound (6):  |frak_k - frak_h (1 + kappa sigma^2/2 + c3 gamma3 sigma^3)|
#                <= K_*(wbar) sigma^4.
# ============================================================================

NTRUE = 400       # truncation used as the proxy for the limit operator
SMAX = {}         # per-wbar admissible relative width (Stage 5)


def kappa_c3(N, w):
    """(kappa_N, c_{3,N}, frak_h^(N)) -- VERBATIM formulas from T5."""
    A = H_mat(N, w)
    _, V = np.linalg.eigh(A)
    v = V[:, -1]
    lam = float(v @ A @ v)
    q2 = float(v @ H_mat(N, w, 2) @ v)
    q3 = float(v @ H_mat(N, w, 3) @ v)
    return w**2 * q2 / lam, w**3 * q3 / (6 * lam), lam


def frak_k_of_P(N, wbar, xis, ps):
    """frak_k^(N)(P,T) = rho( sum_j p_j H^(N)(wbar(1+xi_j)) )."""
    M = np.zeros((N, N))
    for xi, p in zip(xis, ps):
        M += p * H_mat(N, wbar * (1.0 + xi))
    return float(np.max(np.linalg.eigvalsh(M)))


def _fixed_shapes():
    """Fixed zero-mean unit-variance discrete shapes eta with |eta| <= ETA_MAX."""
    out = [("2pk-sym", np.array([-1.0, 1.0]), np.array([0.5, 0.5]))]
    a = np.sqrt(2.0)
    out.append(("3pk-sym", np.array([-a, 0.0, a]), np.array([0.25, 0.5, 0.25])))
    e = np.array([-0.8, -0.2, 1.6]); p = np.array([0.4, 0.4, 0.2])
    e = e - p @ e
    e = e / np.sqrt(p @ e**2)
    out.append(("3pk-asym", e, p))
    e = np.array([-1.5, -0.5, 0.5, 1.5]); p = np.array([0.15, 0.35, 0.35, 0.15])
    e = e - p @ e
    e = e / np.sqrt(p @ e**2)
    out.append(("4pk-sym", e, p))
    return out


def run_stage5():
    print("=== STAGE 5: K_*(wbar) assembly + finite-N consistency check ===")
    t0 = time.time()
    print("  SCOPE (reviewer finding 6, wording): the numeric test below evaluates")
    print("  frak_k at N=400, so what it checks is T5's FINITE-N theorem at N=400")
    print("  together with K_400 <= K_* -- it is a consistency + non-vacuity +")
    print("  sigma^4-rate check, NOT a direct verification of the true-operator")
    print("  bound (6).  (6) for the true operator follows from the limit argument")
    print("  of the note's Sec. 3; the N=400 proxy error is NOT quantified by rho_N")
    print("  (rho_400 ~ 4e-3..1.4e-1 >> K_* sigma^4 ~ 7e-7), it is supported by the")
    print("  N-independence evidence printed at the end of this stage.")
    print(f"  K_*(w) = w^4 B4^inf eta_max^2/24 + (1/2) w^4 (B2^inf)^2 / g_alt")
    print(f"  {'wbar':>7}  {'g_alt':>9}  {'B2^inf':>9}  {'B4^inf':>9}  {'K_tail':>11}  "
          f"{'K_2nd':>11}  {'K_*':>11}  {'s_max':>9}")
    KSTAR = {}
    for w in WBARS:
        B2, B4 = BINF[w]
        g = GSTAR[w]
        # The K_2nd form (1/2) w^4 B2^2 / g comes from simplifying the exact
        # denominator g - 2||E|| >= g/2, which needs ||E|| <= (1/2) w^2 B2 s^2 to
        # satisfy  w^2 B2 s^2 < g/2.  With the now-rigorous (hence looser) B2,
        # the global s=SIGMA0=0.25 VIOLATES this at wbar=5, so we certify a
        # per-wbar admissible width instead of printing a caveat and moving on:
        #     s_max(w) = 0.95 * sqrt( g / (2 w^2 B2) ),  capped at SIGMA0.
        # (The 0.95 keeps the inequality strict.)  B2/B4 are NOT recomputed on
        # the shrunken support-image interval -- they were taken over the larger
        # J built from SIGMA0, so they remain valid sups; this keeps s_max free
        # of any circular dependence on itself.
        s_max = min(SIGMA0, 0.95 * np.sqrt(g / (2.0 * w**2 * B2)))
        SMAX[w] = s_max
        K_tail = w**4 * B4 * ETA_MAX**2 / 24.0
        K_2nd = 0.5 * w**4 * B2**2 / g
        K = K_tail + K_2nd
        KSTAR[w] = K
        cond = w**2 * B2 * s_max**2 < g / 2.0
        _check("K_* finite positive", np.isfinite(K) and K > 0, f"w={w}")
        _check("uniform-denominator condition holds at the certified s_max(wbar)",
               cond, f"w={w}: w^2 B2 s^2={w**2*B2*s_max**2:.4f} vs g/2={g/2:.4f}")
        print(f"  {w:>7.4f}  {g:>9.6f}  {B2:>9.4f}  {B4:>9.3f}  {K_tail:>11.4e}  "
              f"{K_2nd:>11.4e}  {K:>11.4e}  {s_max:>9.4f}")
    print(f"  (s_max < {SIGMA0} only at wbar=5, where the rigorous B2 is large "
          "enough that sigma=0.25 would break g-2||E|| >= g/2.)")

    print(f"\n  INDEPENDENT test of (6) at N={NTRUE} (limit proxy), fixed shapes, "
          f"sigma in [0.01, {SIGMA0}]:")
    print(f"    {'wbar':>7}  {'shape':>9}  {'max |rem|/s^4':>14}  {'K_*':>11}  "
          f"{'ratio':>7}  {'r(.01)/r(.25)':>13}  {'settle':>8}")
    worst_ratio = 0.0
    worst_spread = 0.0
    worst_settle = 0.0
    for w in WBARS:
        kap, c3, hN = kappa_c3(NTRUE, w)
        K = KSTAR[w]
        # only widths inside the certified regime for this wbar
        sigmas = tuple(s for s in (0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25)
                       if s <= SMAX[w])
        for name, eta, p in _fixed_shapes():
            g3 = float(p @ eta**3)
            per_sigma = []
            for sigma in sigmas:
                exact = frak_k_of_P(NTRUE, w, sigma * eta, p)
                pred = hN * (1.0 + 0.5 * kap * sigma**2 + c3 * g3 * sigma**3)
                per_sigma.append(abs(exact - pred) / sigma**4)
            worst = max(per_sigma)
            r = worst / K
            worst_ratio = max(worst_ratio, r)
            _check("finite-N remainder <= K_* sigma^4 at N=400", r <= 1.0,
                   f"w={w} {name}: {worst:.4e} > {K:.4e}")
            # The sigma^4 RATE is genuine iff |rem|/sigma^4 CONVERGES as sigma->0.
            # (i) it must not diverge: a merely-sigma^2 remainder would inflate the
            #     ratio by (0.25/0.01)^2 = 625 between the two ends;
            # (ii) it must settle: the two smallest sigmas must agree to <25%.
            spread = per_sigma[0] / per_sigma[-1]
            settle = abs(per_sigma[0] - per_sigma[1]) / per_sigma[1]
            worst_spread = max(worst_spread, spread)
            worst_settle = max(worst_settle, settle)
            _check("|rem|/sigma^4 does not diverge as sigma->0", spread < 30.0,
                   f"w={w} {name}: r(lo)/r(hi)={spread:.2f} (sigma^2 rate would give ~625)")
            _check("|rem|/sigma^4 settles as sigma->0", settle < 0.25,
                   f"w={w} {name}: |r(.01)-r(.02)|/r(.02)={settle:.3f}")
            print(f"    {w:>7.4f}  {name:>9}  {worst:>14.5e}  {K:>11.4e}  {r:>7.4f}"
                  f"  {spread:>13.2f}  {settle:>8.4f}")
    _check("K_* not vacuous (some case fills >1e-4 of it)", worst_ratio > 1e-4,
           f"worst ratio {worst_ratio:.2e}")
    print(f"  worst |rem|/(K_* sigma^4) over all (wbar, shape) = {worst_ratio:.4f}  "
          f"(<1 => holds; not ~0 => not vacuous)")
    print(f"  worst r(.01)/r(.25) = {worst_spread:.2f} (a sigma^2 remainder would give "
          f"~625) and worst")
    print(f"  |r(.01)-r(.02)|/r(.02) = {worst_settle:.4f}: |rem|/sigma^4 CONVERGES to a "
          "finite limit,")
    print("  i.e. the sigma^4 rate is real and K_* bounds its coefficient.")

    # N-independence of the measured remainder: this is what actually supports
    # treating N=400 as the limit proxy (rho_400 cannot, see the scope note).
    print("\n  N-independence of |rem|/sigma^4 (worst shape 3pk-asym, sigma=0.01):")
    print(f"    {'wbar':>7} " + "".join(f"{('N='+str(N)):>13}" for N in (50, 100, 200, 400, 600)))
    name, eta, p = [s for s in _fixed_shapes() if s[0] == "3pk-asym"][0]
    g3 = float(p @ eta**3)
    for w in WBARS:
        vals = []
        for N in (50, 100, 200, 400, 600):
            kap, c3, hN = kappa_c3(N, w)
            exact = frak_k_of_P(N, w, 0.01 * eta, p)
            pred = hN * (1.0 + 0.5 * kap * 1e-4 + c3 * g3 * 1e-6)
            vals.append(abs(exact - pred) / 1e-8)
        rel = (max(vals) - min(vals)) / max(vals)
        _check("|rem|/sigma^4 is N-independent for N in [50,600]", rel < 1e-3,
               f"w={w} spread={rel:.2e}")
        print(f"    {w:>7.4f} " + "".join(f"{v:>13.6f}" for v in vals))
    print("    => the measured coefficient is already converged well below N=400,")
    print("       so N=400 is a sound proxy even though rho_400 does not certify it.")
    _stage_summary("Stage 5", f"({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    stages = sys.argv[1:] or ["0", "1", "2", "3", "3B", "4", "5"]
    stages = [s.upper() for s in stages]
    T0 = time.time()
    if "0" in stages:
        run_stage0()
    if "1" in stages:
        run_stage1()
    if "2" in stages:
        run_stage2()
    if "3" in stages:
        run_stage3()
    if "3B" in stages:
        run_stage3b()
    if "4" in stages:
        run_stage4()
    if "5" in stages:
        run_stage5()
    print(f"\nALL {PASS[0]} checks passed  ({time.time()-T0:.1f}s).")
