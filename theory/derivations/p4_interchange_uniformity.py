r"""
P4 . N -> infinity INTERCHANGE UNIFORMITY -- computer-assisted certificate for
the three pillars of `P4_interchange_uniformity.md` (v0.1, commit 9b0bacf).

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
      Decays like N^-1 and is what actually makes N_0 ~ 30-80 instead of ~500.
      The derivation note predicted "N_0 ~ 30-40 at wbar ~ 1"; only route (b)
      reaches that, so route (b) is the one that carries the certificate and
      route (a) is kept as an independent cross-check on the same quantity.

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
                  come from the Pillar-1 Schur envelopes.
  delta := mu - alpha.

PILLAR 3 (Stage 3).  N_0(w) := min{N : 2 rho_N(w) < delta(w)/2} and, by Weyl
(eq. (4) of the note)  gap_N >= delta - 2 rho_N,
    inf_{N>=2} gap_N(w) >= g_*(w) := min( min_{2<=N<N_0} gap_N(w), delta(w)/2 ).
The finite direct check is run over ALL N in [2, 400] (not only N < N_0).

PILLAR 4 (Stage 4).  Compression monotonicity ||H^(N)(j)|| <= ||H^(N')(j)||
(N<N') verified, and certified upper bounds B_2^inf, B_4^inf on the LIMIT
operator derivative norms, sup over the support-image interval J.

ASSEMBLY + INDEPENDENT TEST (Stage 5).  K_*(w) = w^4 B_4^inf eta_max^2/24
+ (1/2) w^4 (B_2^inf)^2 / g_*, then a direct test of eq. (6): for fixed-shape
compactly supported P, compare the true frak_k (N=400 proxy for the limit)
against frak_h (1 + kappa sigma^2/2 + c_3 gamma_3 sigma^3) and check the
difference divided by sigma^4 stays <= K_* and does not collapse to 0 (i.e. the
bound is not vacuous).

NUMERICAL DISCIPLINE (mandated): numpy float64 `eigvalsh` only (never generic
`eig`, never mpmath eigensolvers); mpmath dps=30 restricted to the Cholesky
certificate, the Rayleigh quotient and the closed-form tail integrals; all loops
have fixed bounds, no open-ended search.
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


_ENV_CACHE = {}


def envelope(w, order=0):
    """(S1, IntTail) for psi: S1 = sum_{k>=1} psi(k),  IntTail(a) = int_a^inf psi.

    S1 is a partial sum to k=2000 plus the integral tail (psi decreasing =>
    sum_{k>K} psi(k) <= int_K^inf psi), so S1 is a rigorous UPPER bound on the
    series.  IntTail is evaluated in mpmath (closed forms below, all elementary).
    """
    key = (round(float(w), 12), order)
    if key in _ENV_CACHE:
        return _ENV_CACHE[key]
    wm = mp.mpf(float(w))
    c = mp.mpf(_psi_c(order, w))
    p = PSI_P[order]

    def int_tail(a):
        a = mp.mpf(a)
        if p == 1:                      # int_a^inf dt/(t^2+w^2) = (pi/2 - atan(a/w))/w
            return c * (mp.pi / 2 - mp.atan(a / wm)) / wm
        if p == 2:                      # int_a^inf dt/(t^2+w^2)^2
            return c * ((mp.pi / 2 - mp.atan(a / wm)) / (2 * wm**3)
                        - a / (2 * wm**2 * (a**2 + wm**2)))
        if p == 3:                      # int_a^inf dt/(t^2+w^2)^3
            return c * (3 * (mp.pi / 2 - mp.atan(a / wm)) / (8 * wm**5)
                        - a / (4 * wm**2 * (a**2 + wm**2) ** 2)
                        - 3 * a / (8 * wm**4 * (a**2 + wm**2)))
        raise ValueError(p)

    K = 2000
    ks = np.arange(1, K + 1, dtype=np.float64)
    S1 = mp.mpf(float(np.sum(_psi_c(order, w) / (ks**2 + float(w) ** 2) ** p))) + int_tail(K)
    _ENV_CACHE[key] = (S1, int_tail)
    return S1, int_tail


def row_env(n, w, order=0):
    """R_bar(n): a MANIFESTLY DECREASING upper envelope for the full row sum
    sum_{m>=0} |H^{(order)}_{nm}(w)|, valid for every n >= 1.

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
    S1, int_tail = envelope(w, order)
    n = mp.mpf(n)
    r = mp.sqrt((n + 1) * (2 * n + 1))
    return (2 * S1 / (2 * n + 1)
            + 2 * S1 / r
            + mp.mpf(psi(float(n) / 2.0, w, order)) / mp.sqrt(2)
            + int_tail(n) / r
            + mp.mpf(psi(float(n), w, order)) / mp.sqrt(2))


def col_env(m, N, w, order=0):
    """C_bar(m) for m < N: upper bound on sum_{n>=N} |H^{(order)}_{nm}(w)|.

      sum_{n>=N} psi(n-m)   = sum_{k>=N-m} psi(k) <= IntTail(N-m-1)   (N-m>=1)
      sum_{n>=N} psi(n+m+1) = sum_{k>=N+m+1} psi(k) <= IntTail(N+m)
      1/sqrt((2n+1)(2m+1)) <= 1/sqrt((2N+1)(2m+1))
    """
    _, int_tail = envelope(w, order)
    return (int_tail(N - m - 1) + int_tail(N + m)) / mp.sqrt((2 * m + 1) * (2 * N + 1))


def schur_offhead(N, w, order=0):
    """||(I-P_N) H^{(order)}(w)||_op bound by the Schur/Holmgren test
        ||A|| <= sqrt( max_n sum_m |A_nm| * max_m sum_n |A_nm| ).
    Rows n >= N: sum over all m <= R_bar(n) <= R_bar(N).
    Columns: m >= N gives sum_{n>=N} <= full row sum at m <= R_bar(N);
             m <  N is the finite exact max of C_bar(m).
    """
    rmax = row_env(N, w, order)
    cmax = rmax
    for m in range(N):
        c = col_env(m, N, w, order)
        if c > cmax:
            cmax = c
    return mp.sqrt(rmax * cmax)


def schur_tailblock(N, w, order=0):
    """||H^{(order)}|_{n,m >= N}||_op <= R_bar(N) (Schur; restricting the column
    range only shrinks both Schur maxima)."""
    return row_env(N, w, order)


def schur_coupling(Na, w, order=0):
    """||B||_op for the head-tail coupling B = (H_nm)_{n<Na, m>=Na}.
    row sums (over m>=Na, n<Na) are exactly C_bar(n); column sums (over n<Na,
    m>=Na) are bounded by the full row sum at m, i.e. R_bar(m) <= R_bar(Na)."""
    cmax = mp.mpf(0)
    for n in range(Na):
        c = col_env(n, Na, w, order)
        if c > cmax:
            cmax = c
    return mp.sqrt(cmax * row_env(Na, w, order))


def hs_offhead(N, w, M=200000):
    """||(I-P_N) H(w)||_op via the HS route of the derivation note (order 0 only).

      ||(I-P_N)H|| <= max_{n>=N}(H1)_nn + ||(I-P_N)(H2+H3)||_HS
                   <= 2 S1/(2N+1)      + sqrt(S(N)),
      S(N) = sum_{n>=N} sum_m ((H2+H3)_nm)^2 <= sum_{n>=N} 2 (T(n)+U(n))/(2n+1),
      using (a+b)^2 <= 2a^2+2b^2 and, with Z2 := 2 sum_{k>=1}[[k]]^2,
        T(n) = Z2/(n+1) + [[n/2]]^2 (3/2 + ln(n)/2)        (H2 inner sum)
        U(n) = [[n]]^2 (3/2 + ln(n)/2) + Z2/(2(2n+1))      (H3 inner sum)
      (same m>=ceil(n/2) vs m<ceil(n/2) split as `row_env`, with
       sum_{m<h} 1/(2m+1) <= 3/2 + ln(n)/2).
    The n-sum is exact for N <= n < M and closed-form beyond, using
    [[n/2]]^2+[[n]]^2 <= 17 w^4/n^4 and int_{M-1}^inf (3/2+ln x/2)/x^5 dx.
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
    U = b_n**2 * lg + Z2 / (2.0 * (2.0 * n + 1.0))
    S = float(np.sum(2.0 * (T + U) / (2.0 * n + 1.0)))
    Mm = float(M - 1)
    S += 1.25 * Z2 / Mm + 17.0 * wf**4 * ((1.5 + 0.5 * np.log(Mm)) / (4 * Mm**4)
                                          + 0.5 / (16 * Mm**4))
    S1, _ = envelope(w, 0)
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
        S1, _ = envelope(w, 0)
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
            _check("row_env dominates true row sums", float(row_env(N, w)) >= rows[N] - 1e-12,
                   f"w={w} N={N}: env={float(row_env(N,w)):.6f} true={rows[N]:.6f}")
            cols = Hb[N:, :].sum(axis=0)
            for m in (0, N // 3, N - 1):
                _check("col_env dominates true col sums",
                       float(col_env(m, N, w)) >= cols[m] - 1e-12,
                       f"w={w} N={N} m={m}")
            # the certified operator-norm bound must dominate a large finite
            # numerical evaluation of ||(I-P_N) H P_Mbig||
            blk = H_mat(Mbig, w)[N:, :]
            true_op = float(np.sqrt(np.max(np.linalg.eigvalsh(blk @ blk.T))))
            for nm, val in (("schur", float(schur_offhead(N, w))), ("hs", hs_offhead(N, w))):
                _check(f"{nm} bound >= true ||(I-P_N)H||", val >= true_op,
                       f"w={w} N={N}: {nm}={val:.6f} true>={true_op:.6f}")

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


def rayleigh_mu(w, Nv=NV_RAY):
    """Certified lower bound mu <= lam_1(H(w)).

    v is the float64 top eigenvector of H^(Nv)(w), zero-padded to l^2; since v
    has FINITE SUPPORT, <v,H v> involves only the head block and is computed
    EXACTLY (mpmath dps=30) -- no truncation error, no eigenvalue solver enters
    the bound (the eigenvector is merely a trial vector; any vector would give a
    valid, just weaker, bound)."""
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
    return num / den


def alpha_certificate_H(w, Na=NA_HEAD):
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
    nB = schur_coupling(Na, w)
    nD = schur_tailblock(Na, w)
    return c + nB + nD, True, c, nB, nD


DELTA = {}


def run_stage2():
    print("=== STAGE 2: PILLAR 2 -- certified limit gap delta(wbar) = mu - alpha ===")
    t0 = time.time()
    print(f"  head Na={NA_HEAD} (Cholesky-deflation, mpmath dps={mp.dps}), "
          f"Rayleigh support Nv={NV_RAY}")
    print(f"  {'wbar':>7}  {'mu (<=lam1)':>12}  {'lam2(A)<=c':>11}  {'||B||':>9}  "
          f"{'||D||':>9}  {'alpha>=lam2':>12}  {'delta':>10}")
    for w in WBARS:
        mu = rayleigh_mu(w)
        alpha, ok, c, nB, nD = alpha_certificate_H(w)
        _check("Cholesky deflation certificate SPD", ok, f"w={w}")
        delta = mu - alpha
        DELTA[w] = float(delta)
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
# STAGE 3 -- PILLAR 3: N_0(w), the gap_N table, and g_*(w).
# ============================================================================

NMAX_SCAN = 400
GSTAR = {}
N0 = {}


def run_stage3():
    print("=== STAGE 3: PILLAR 3 -- N_0(wbar), gap_N vs N, and g_*(wbar) ===")
    t0 = time.time()
    print("  Weyl (eq.4): gap_N >= delta - 2 rho_N.  Note P_N H P_N (zero-padded)")
    print("  has lam_2 >= lam_2(H^(N)), so the bound on the PADDED operator's gap")
    print("  is a fortiori a lower bound on the matrix gap gap_N -- sign is safe.")
    print(f"  {'wbar':>7}  {'delta':>9}  {'N_0':>5}  {'2rho_{N_0}':>11}  "
          f"{'min gap 2..399':>14}  {'argmin':>6}  {'g_*':>9}")
    tables = {}
    for w in WBARS:
        delta = DELTA[w]
        n0 = None
        for N in range(2, 2001):
            if 2.0 * rho_bound(N, w) < delta / 2.0:
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
        gmin_pre = float(np.min(gaps[2:n0])) if n0 > 2 else np.inf
        gstar = min(gmin_pre, delta / 2.0)
        GSTAR[w] = gstar
        _check("g_* > 0", gstar > 0, f"w={w}")
        gall = float(np.min(gaps[2:NMAX_SCAN + 1]))
        amin = int(np.nanargmin(gaps[2:NMAX_SCAN + 1])) + 2
        _check("no intermediate N in [2,400] dips below g_*", gall >= gstar - 1e-12,
               f"w={w} min={gall:.6f} g_*={gstar:.6f} at N={amin}")
        _check("Weyl bound consistent at N_0", gaps[n0] >= delta - 2 * rho_bound(n0, w) - 1e-12,
               f"w={w}")
        print(f"  {w:>7.4f}  {delta:>9.6f}  {n0:>5}  {2*rho_bound(n0, w):>11.6f}  "
              f"{gall:>14.6f}  {amin:>6}  {gstar:>9.6f}")

    print("\n  gap_N(wbar) vs N  (converges from ABOVE to the limit gap; never collapses):")
    Ns = [2, 4, 8, 16, 32, 64, 128, 256, 400]
    print(f"    {'wbar':>7} " + "".join(f"{('N='+str(N)):>10}" for N in Ns) + f"{'delta(LB)':>11}")
    for w in WBARS:
        print(f"    {w:>7.4f} " + "".join(f"{tables[w][N]:>10.5f}" for N in Ns)
              + f"{DELTA[w]:>11.5f}")
    # cross-check against T5's reported observation at wbar=1
    g1 = tables[1.0]
    _check("T5 cross-check gap_2(1)~1.137", abs(g1[2] - 1.137) < 5e-3, f"{g1[2]:.5f}")
    _check("T5 cross-check gap_32(1)~0.717", abs(g1[32] - 0.717) < 5e-3, f"{g1[32]:.5f}")
    print(f"  T5 cross-check at wbar=1: gap_2={g1[2]:.5f} (T5: 1.137), "
          f"gap_32={g1[32]:.5f} (T5: 0.717), gap_400={g1[400]:.5f} -> stable, no collapse.")
    _stage_summary("Stage 3", f"({time.time()-t0:.1f}s)")


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


def cert_sup_on_J(w, order, ngrid=201):
    """Certified sup over J = [wbar(1-sigma0 eta_max), wbar(1+sigma0 eta_max)] of
    the certified limit norm: grid max + Lipschitz slack L h/2 with
    L = grid max of the (order+1) certified norm (since d/dw ||H^(j)(w)||
    <= ||H^(j+1)(w)||)."""
    a = max(w * (1.0 - SIGMA0 * ETA_MAX), 1e-6)
    b = w * (1.0 + SIGMA0 * ETA_MAX)
    ws = np.linspace(a, b, ngrid)
    vals = [cert_opnorm_inf(x, order) for x in ws]
    L = max(cert_opnorm_inf(x, order + 1) for x in ws)
    return max(vals) + L * (b - a) / (ngrid - 1) / 2.0


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
    print("=== STAGE 5: K_*(wbar) assembly + INDEPENDENT test of bound (6) ===")
    t0 = time.time()
    print(f"  K_*(w) = w^4 B4^inf eta_max^2/24 + (1/2) w^4 (B2^inf)^2 / g_*")
    print(f"  {'wbar':>7}  {'g_*':>9}  {'B2^inf':>9}  {'B4^inf':>9}  {'K_tail':>11}  "
          f"{'K_2nd':>11}  {'K_*':>11}  {'sigma0 ok':>9}")
    KSTAR = {}
    for w in WBARS:
        B2, B4 = BINF[w]
        g = GSTAR[w]
        K_tail = w**4 * B4 * ETA_MAX**2 / 24.0
        K_2nd = 0.5 * w**4 * B2**2 / g
        K = K_tail + K_2nd
        KSTAR[w] = K
        # the note's uniform-denominator condition: wbar^2 B2 sigma0^2 < g_*/2
        cond = w**2 * B2 * SIGMA0**2 < g / 2.0
        _check("K_* finite positive", np.isfinite(K) and K > 0, f"w={w}")
        print(f"  {w:>7.4f}  {g:>9.6f}  {B2:>9.4f}  {B4:>9.3f}  {K_tail:>11.4e}  "
              f"{K_2nd:>11.4e}  {K:>11.4e}  {('yes' if cond else 'NO'):>9}")
        if not cond:
            print(f"          ^ wbar^2 B2 sigma0^2 = {w**2*B2*SIGMA0**2:.4f} >= g_*/2 = "
                  f"{g/2:.4f}: bound (6) then holds only for smaller sigma0 at this wbar")

    print(f"\n  INDEPENDENT test of (6) at N={NTRUE} (limit proxy), fixed shapes, "
          f"sigma in [0.01, {SIGMA0}]:")
    print(f"    {'wbar':>7}  {'shape':>9}  {'max |rem|/s^4':>14}  {'K_*':>11}  "
          f"{'ratio':>7}  {'r(.01)/r(.25)':>13}  {'settle':>8}")
    sigmas = (0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25)
    worst_ratio = 0.0
    worst_spread = 0.0
    worst_settle = 0.0
    for w in WBARS:
        kap, c3, hN = kappa_c3(NTRUE, w)
        K = KSTAR[w]
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
            _check("bound (6) holds: |rem|/sigma^4 <= K_*", r <= 1.0,
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
                   f"w={w} {name}: r(.01)/r(.25)={spread:.2f} (sigma^2 rate would give ~625)")
            _check("|rem|/sigma^4 settles as sigma->0", settle < 0.25,
                   f"w={w} {name}: |r(.01)-r(.02)|/r(.02)={settle:.3f}")
            print(f"    {w:>7.4f}  {name:>9}  {worst:>14.5e}  {K:>11.4e}  {r:>7.4f}"
                  f"  {spread:>13.2f}  {settle:>8.4f}")
    _check("bound (6) not vacuous (some case fills >1e-4 of it)", worst_ratio > 1e-4,
           f"worst ratio {worst_ratio:.2e}")
    print(f"  worst |rem|/(K_* sigma^4) over all (wbar, shape) = {worst_ratio:.4f}  "
          f"(<1 => holds; not ~0 => not vacuous)")
    print(f"  worst r(.01)/r(.25) = {worst_spread:.2f} (a sigma^2 remainder would give "
          f"~625) and worst")
    print(f"  |r(.01)-r(.02)|/r(.02) = {worst_settle:.4f}: |rem|/sigma^4 CONVERGES to a "
          "finite limit,")
    print("  i.e. the sigma^4 rate is real and K_* bounds its coefficient.")
    _stage_summary("Stage 5", f"({time.time()-t0:.1f}s)")


if __name__ == "__main__":
    stages = sys.argv[1:] or ["0", "1", "2", "3", "4", "5"]
    T0 = time.time()
    if "0" in stages:
        run_stage0()
    if "1" in stages:
        run_stage1()
    if "2" in stages:
        run_stage2()
    if "3" in stages:
        run_stage3()
    if "4" in stages:
        run_stage4()
    if "5" in stages:
        run_stage5()
    print(f"\nALL {PASS[0]} checks passed  ({time.time()-T0:.1f}s).")
