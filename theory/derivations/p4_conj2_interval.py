r"""
P4 . KAY CONJECTURE 2 -- Lemma B (h(wbar) <= frak_g(2) wbar^2) on the MIDDLE and
LARGE wbar segments, plus the MONOTONICITY of h beyond KAY's sqrt(2).

RUN:  python theory/derivations/p4_conj2_interval.py            (S0 S3 S2)
      python theory/derivations/p4_conj2_interval.py M          (the monotonicity
                                                                 stage, separate
                                                                 because it is the
                                                                 expensive one)
      python theory/derivations/p4_conj2_interval.py ALL

SCOPE (this file owns Segment II, Segment III and Sec.3 of
`theory/derivations/P4_conjecture2.md`; Segment I -- the small-wbar analytic
perturbation on (0, wbar_0] -- is a separate deliverable).

    Lemma B :  h(wbar) := lam_1(H(wbar))  <  frak_g(2) wbar^2   for all wbar > 0.

  Segment II  (Stage S2):  [0.3, WBAR2]  by a per-cell interval certificate.
  Segment III (Stage S3):  [wbar_1, inf) by a CLOSED-FORM analytic bound, with
                           wbar_1 <= WBAR2 so the two segments overlap.
  Sec. 3      (Stage M ):  h'(wbar) > 0, certified per cell on [sqrt(2), W].

NOTHING IS RE-DERIVED HERE.  Every envelope, Cholesky certificate and segment
spec is IMPORTED UNCHANGED from the certified modules
  * `p4_interchange_uniformity` (1668 checks): `H_mat`, the Schur segment
    envelopes `schur_coupling` / `schur_tailblock` / `schur_full_seg`,
    `spec_seg` / `spec_point` / `envelope` / `col_env` / `psi_spec`,
    `mu2_cell_lb`, `alpha_certificate_H`, `ROUND_PAD`;
  * `p4_interval_extension`: `schur_full_seg_w0` (the weighted-row-0 Lipschitz
    constant, reviewer-adopted), `BASE_RATIO`, `head_size`;
  * `p4d_kappa_pencil`: `lambda1_upper_pencil` -- the NON-deflated Cholesky head
    certificate lam_1(A_head) <= c (its mutilde=0 branch is bit-identical to the
    phonon operator, Stage A gate of `p4d_coulomb_certify`).
Two small envelopes are NEW and derived in full below: `coupling_restricted`
(Sec. M.2) and the analytic row bound of Sec. S3.2.  Neither touches a tracked
base file.

NUMERICAL DISCIPLINE (inherited, mandated): every quantity that enters a bound is
mpmath dps=30 arithmetic on float64-rounded matrix entries, padded by ROUND_PAD =
1e-12.  Sizing of that pad at the head sizes used here: the entry rounding is
<= 1.1e-16 per entry and every entry of H, H' is O(1), so the operator-norm
residual of the float64 image of any block used below is bounded by its Frobenius
norm <= 1.1e-16 * sqrt(rows*cols) <= 6.6e-14 (the largest block is the 1192 x 298
residual block of Stage M); the mpmath arithmetic itself contributes ~1e-24.
ROUND_PAD is therefore 15-30x above the worst residual and 1e-9 x below the
smallest certified quantity, so it costs nothing and makes every number returned
a bound on the EXACT operator rather than on its float64 image.  numpy float64
`eigvalsh` / `eigh` appear ONLY as trial-vector generators and as EVIDENCE
checks.  All loops have fixed bounds; the one adaptive refinement (Stage M) has a
hard cap.

frak_g(2) IS USED ONLY THROUGH ITS CERTIFIED LOWER BOUND  frak_g(2) >= 1.31714
(P4-C two-sided enclosure [1.31714, 1.31820], `P4_rigorous_tc_bound.md` Sec.P4-C).
Since the claim is h < frak_g(2) wbar^2, only the LOWER bound is load-bearing.
"""

import sys
import os
import time

import numpy as np
from mpmath import mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_interchange_uniformity as piu                                  # noqa: E402
from p4_interchange_uniformity import (H_mat, H_matrix, ROUND_PAD,       # noqa: E402
                                       schur_coupling, schur_tailblock,
                                       schur_full_seg, spec_seg, spec_point,
                                       envelope, col_env, psi_spec,
                                       mu2_cell_lb, alpha_certificate_H,
                                       _brk_vec)
from p4_interval_extension import (schur_full_seg_w0, BASE_RATIO,        # noqa: E402
                                   head_size)
from p4d_kappa_pencil import lambda1_upper_pencil                        # noqa: E402

mp.dps = 30

PASS = [0]
_LAST = [0]


def chk(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name} {detail}")
    PASS[0] += 1


def _summary(tag, note=""):
    k = PASS[0] - _LAST[0]
    _LAST[0] = PASS[0]
    print(f"  [{tag}] {k} checks passed. {note}")


# The ONLY external numerical constant.  P4-C certified enclosure of the KAY I
# top eigenvalue frak_g(2) = rho(G(2)):  1.31714 <= frak_g(2) <= 1.31820.
G2_LO = mp.mpf("1.31714")

SQRT2 = mp.sqrt(2)


# ===========================================================================
# Shared helpers
# ===========================================================================

def H_rect(Nrow, Ncol, w, order=0):
    """The rectangular block (H^{(order)}(w))_{n<Nrow, m<Ncol}, same entry
    formulas as `piu.H_mat` (checked bit-comparable against it on the square
    overlap in Stage S0)."""
    n = np.arange(Nrow)
    m = np.arange(Ncol)
    invn = 1.0 / np.sqrt(2.0 * n + 1.0)
    invm = 1.0 / np.sqrt(2.0 * m + 1.0)
    A = np.abs(n[:, None] - m[None, :])
    S = n[:, None] + m[None, :] + 1
    off = _brk_vec(A, w, order)
    off[A == 0] = 0.0                     # H2 has no diagonal
    M = (off + _brk_vec(S, w, order)) * np.outer(invn, invm)
    d = min(Nrow, Ncol)
    if d > 1:
        cum = np.concatenate(([0.0], np.cumsum(2.0 * _brk_vec(np.arange(1, d), w, order))))
    else:
        cum = np.zeros(1)
    idx = np.arange(d)
    M[idx, idx] -= cum / (2.0 * idx + 1.0)
    return M


def _mp_vec(v):
    return [mp.mpf(float(x)) for x in v]


def _mp_matvec(Af, vm):
    """(A v) in exact mpmath arithmetic on the float64 entries of A."""
    out = []
    nc = len(vm)
    for i in range(Af.shape[0]):
        row = Af[i]
        s = mp.zero
        for j in range(nc):
            s += mp.mpf(float(row[j])) * vm[j]
        out.append(s)
    return out


def _mp_quad(Af, vm):
    """<v, A v> in exact mpmath arithmetic on the float64 entries of A."""
    y = _mp_matvec(Af, vm)
    return sum(vm[i] * y[i] for i in range(len(vm)))


# ===========================================================================
# STAGE S0 -- self-check gate: the imported machinery is the machinery the base
# files certified, and the two new local builders agree with it.
# ===========================================================================

def run_S0():
    print("=== STAGE S0: self-check gate (imports + local builders) ===")
    t0 = time.time()
    # (a) H_rect == H_mat == H_matrix on the square overlap, all orders used here
    for N in (7, 31) :
        for w in (0.3, 1.0, 2.5, 7.0):
            for order in (0, 1, 2):
                R = H_rect(N, N, w, order)
                d1 = np.max(np.abs(R - H_mat(N, w, order)))
                d2 = np.max(np.abs(H_mat(N, w, order) - H_matrix(N, w, order)))
                chk("H_rect == H_mat (square overlap)", d1 == 0.0, f"N={N} w={w} o={order} {d1:.3e}")
                chk("H_mat == H_matrix (verbatim T5)", d2 <= 1e-13, f"N={N} w={w} o={order}")
    # (b) rectangular block == the corresponding sub-block of a bigger square one
    for w in (0.5, 3.0):
        big = H_mat(60, w, 0)
        chk("H_rect block == sub-block of H_mat",
            np.max(np.abs(H_rect(60, 12, w, 0) - big[:, :12])) == 0.0, f"w={w}")
    # (c) the imported non-deflated Cholesky head certificate really bounds
    #     lam_1(A_head) (float64 evidence) and is >= the float64 value
    for w in (0.3, 1.0, 2.5):
        ub, ok, c, nB, nD = lambda1_upper_pencil(w, 0.0, 0, Na=70)
        chk("imported Cholesky head certificate SPD", ok, f"w={w}")
        l1 = float(np.linalg.eigvalsh(H_mat(70, w))[-1])
        chk("c >= lam_1(A_head) (float64 evidence)", float(c) >= l1, f"w={w}")
        chk("c + ||B|| + ||D|| >= lam_1(H) at N=1500 (float64 evidence)",
            float(ub) >= float(np.linalg.eigvalsh(H_mat(1500, w))[-1]), f"w={w}")
    # (d) anchor assertion ONLY (never used in a bound): P4-C's frak_g(2) enclosure
    gG = float(np.max(np.linalg.eigvalsh(piu.G_matrix(300, 2))))
    chk("P4-C anchor frak_g(2) in [1.31714,1.31820] (ASSERTION ONLY)",
        1.31714 <= gG <= 1.31820, f"g^(300)(2)={gG:.7f}")
    chk("G2_LO is the P4-C lower end", G2_LO == mp.mpf("1.31714"))
    _summary("Stage S0", f"({time.time()-t0:.1f}s)")


# ===========================================================================
# STAGE S3 -- SEGMENT III: a closed-form analytic bound for large wbar.
#
# S3.1  THE SERIES BOUND.  S1(wbar) := sum_{k>=1} [[k]](wbar),
#       [[k]] = wbar^2/(k^2+wbar^2).  t -> wbar^2/(t^2+wbar^2) is POSITIVE and
#       DECREASING on [0,inf), so the integral comparison for a decreasing
#       function gives
#           S1 = sum_{k>=1} f(k) <= int_0^inf f(t) dt
#                                 = wbar^2 * (1/wbar) * (pi/2) = pi wbar / 2.
#       (Exact closed form, used only for reporting: S1 = (pi w coth(pi w)-1)/2.)
#
# S3.2  THE SCHUR ROW BOUND.  H is symmetric, so the Schur/Holmgren test
#       ||A|| <= sqrt(rowmax * colmax) reduces to  ||H|| <= max_n R_n,
#       R_n := sum_m |H_nm|.  Split R_n = T1 + T2 + T3 exactly as `row_env`
#       does (H1 diagonal, H2 off-diagonal, H3), with s := (2n+1)^{-1/2}:
#           T1 = 2 sigma_n/(2n+1),  sigma_n := sum_{k=1}^n [[k]] <= S1
#           T2 = s sum_{m!=n} [[|n-m|]]/sqrt(2m+1)
#           T3 = s sum_{m>=0} [[n+m+1]]/sqrt(2m+1)
#
#       ROW n = 0  (T1 is empty: the sum k=1..0 has no terms).
#           T2 = sum_{m>=1} [[m]]/sqrt(2m+1) <= (1/sqrt3) S1
#           T3 = [[1]] + sum_{m>=1}[[m+1]]/sqrt(2m+1) <= [[1]] + (S1-[[1]])/sqrt3
#         ==> R_0 <= 2 S1/sqrt3 + [[1]](1 - 1/sqrt3) <= (1 + 1/sqrt3) S1 = 1.5774 S1
#         (last step: [[1]] <= S1).
#
#       ROW n = 1.
#           T1 = 2[[1]]/3
#           T2 = (1/sqrt3)( [[1]]  +  sum_{k>=1} [[k]]/sqrt(2k+3) )
#              <= [[1]]/sqrt3 + S1/sqrt15
#           T3 = (1/sqrt3)( [[2]]  +  sum_{m>=1}[[m+2]]/sqrt(2m+1) )
#              <= [[2]]/sqrt3 + (1/3)(S1 - [[1]] - [[2]])
#         ==> R_1 <= [[1]](1/3 + 1/sqrt3) + [[2]](1/sqrt3 - 1/3) + S1(1/3 + 1/sqrt15)
#                 <= (1/3 + 2/sqrt3 + 1/sqrt15) S1 = 1.746233... S1
#         (both bracket coefficients are >= 0, so [[1]],[[2]] <= S1 may be used).
#
#       ROW n >= 2.
#           T1 <= 2 S1/(2n+1)
#           T2, m > n : sum_{m>n}[[m-n]]/sqrt(2m+1) <= S1/sqrt(2n+3)
#           T2, m < n : sum_{j=1}^{n} [[j]] * (2(n-j)+1)^{-1/2}.  The first factor
#                       DECREASES in j, the second INCREASES in j, so CHEBYSHEV's
#                       sum inequality for oppositely ordered sequences gives
#                         <= (1/n) (sum_j [[j]]) (sum_{m=0}^{n-1}(2m+1)^{-1/2})
#                         <= (sigma_n/n) sqrt(2n-1) <= (S1/n) sqrt(2n-1),
#                       using sum_{m=0}^{n-1}(2m+1)^{-1/2} <= sqrt(2n-1)
#                       (induction: sqrt(2n+1)-sqrt(2n-1) = 2/(sqrt(2n+1)+sqrt(2n-1))
#                        >= 1/sqrt(2n+1)).
#           T3 <= s sum_{k>n}[[k]] <= S1/sqrt(2n+1)
#         ==> R_n <= S1 * phi(n),
#             phi(n) = 2/(2n+1) + sqrt(2n-1)/(n sqrt(2n+1))
#                      + 1/sqrt((2n+1)(2n+3)) + 1/sqrt(2n+1).
#         phi is DECREASING for n >= 1: the 1st, 3rd, 4th terms visibly are, and
#         the 2nd is 1/(n sqrt((2n+1)/(2n-1))) with n^2(2n+1)/(2n-1) increasing
#         (derivative numerator 2n(4n^2-2n-1) > 0).  Hence max_{n>=2} phi = phi(2)
#         = 1.403543...
#
#       ==> ||H(wbar)|| <= C S1(wbar),   C := 1/3 + 2/sqrt3 + 1/sqrt15 = 1.746233,
#       the maximum of the three row constants (row 1 is the binding row).
#
# S3.3  THE CLOSURE.  h(wbar) <= ||H(wbar)|| <= C S1 <= C pi wbar/2, while the
#       target is frak_g(2) wbar^2 >= 1.31714 wbar^2.  Hence for
#           wbar > wbar_1 := C pi / (2 * 1.31714) = 2.08260...
#       we get h(wbar) <= C pi wbar/2 < 1.31714 wbar^2 <= frak_g(2) wbar^2,
#       STRICTLY, with no computation at all.  Segment II below is certified up
#       to WBAR2 = 2.5 > wbar_1, so the two segments overlap.
# ===========================================================================

# the three proven row constants
C_ROW0 = 1 + 1 / mp.sqrt(3)
C_ROW1 = mp.mpf(1) / 3 + 2 / mp.sqrt(3) + 1 / mp.sqrt(15)
C_ROWN = (mp.mpf(2) / 5 + mp.sqrt(3) / (2 * mp.sqrt(5))
          + 1 / mp.sqrt(35) + 1 / mp.sqrt(5))            # = phi(2)
C_SCHUR = max(C_ROW0, C_ROW1, C_ROWN)
WBAR1 = C_SCHUR * mp.pi / (2 * G2_LO)
WBAR2 = 2.5          # right end of Segment II (> WBAR1, with margin)


def phi_row(n):
    n = mp.mpf(n)
    return (2 / (2 * n + 1) + mp.sqrt(2 * n - 1) / (n * mp.sqrt(2 * n + 1))
            + 1 / mp.sqrt((2 * n + 1) * (2 * n + 3)) + 1 / mp.sqrt(2 * n + 1))


def S1_closed(w):
    """S1(w) = sum_{k>=1} w^2/(k^2+w^2) = (pi w coth(pi w) - 1)/2 (classical)."""
    w = mp.mpf(w)
    return (mp.pi * w / mp.tanh(mp.pi * w) - 1) / 2


def row_bound_analytic(n, w):
    """The proven per-row upper bound on sum_m |H_nm(w)| (the sharp form, i.e.
    BEFORE replacing [[1]],[[2]] by S1)."""
    S1 = S1_closed(w)
    b1 = mp.mpf(w) ** 2 / (1 + mp.mpf(w) ** 2)
    b2 = mp.mpf(w) ** 2 / (4 + mp.mpf(w) ** 2)
    if n == 0:
        return 2 * S1 / mp.sqrt(3) + b1 * (1 - 1 / mp.sqrt(3))
    if n == 1:
        return (b1 * (mp.mpf(1) / 3 + 1 / mp.sqrt(3))
                + b2 * (1 / mp.sqrt(3) - mp.mpf(1) / 3)
                + S1 * (mp.mpf(1) / 3 + 1 / mp.sqrt(15)))
    return S1 * phi_row(n)


def row_sum_exact_ub(n, w, M=20000):
    """Rigorous upper bound on the EXACT row sum sum_{m>=0}|H_nm(w)| : explicit
    for m < M, closed-form integral tail beyond (the tail of row n is the same
    expression `col_env` bounds for column n, by the symmetry of H)."""
    spec = spec_point(w, 0)
    row = np.abs(H_rect(n + 1, M, float(w), 0)[n])
    return mp.mpf(float(row.sum())) + col_env(n, M, spec)


def run_S3():
    print("=== STAGE S3: SEGMENT III -- closed-form bound ||H|| <= C S1 <= C pi w/2 ===")
    t0 = time.time()
    print(f"  row constants:  n=0 : {mp.nstr(C_ROW0, 10)}   n=1 : {mp.nstr(C_ROW1, 10)}"
          f"   n>=2 : phi(2) = {mp.nstr(C_ROWN, 10)}")
    print(f"  ==> C = 1/3 + 2/sqrt3 + 1/sqrt15 = {mp.nstr(C_SCHUR, 12)}  (row 1 binds)")

    # --- (1) S1 <= pi w / 2, and S1_closed is the true series -------------
    for w in (0.05, 0.3, 1.0, 2.0, 5.0, 20.0, 100.0):
        S1c = S1_closed(w)
        S1e = envelope(spec_point(w, 0))[0]          # certified UPPER bound on the series
        chk("S1 closed form <= certified series envelope", S1c <= S1e + mp.mpf("1e-20"),
            f"w={w} {mp.nstr(S1c,12)} vs {mp.nstr(S1e,12)}")
        chk("S1 envelope - S1 closed form is tiny (both are the same series; the"
            " envelope replaces the k>2000 tail by its integral)",
            S1e - S1c < mp.mpf("1e-4") * (1 + S1c), f"w={w}")
        chk("integral comparison S1 <= pi w / 2", S1c <= mp.pi * mp.mpf(w) / 2, f"w={w}")

    # --- (2) phi decreasing, and the three row constants dominate ----------
    prev = phi_row(1)
    for n in range(2, 400):
        cur = phi_row(n)
        chk("phi(n) decreasing", cur < prev, f"n={n}")
        chk("phi(n) <= phi(2) for n>=2", cur <= C_ROWN + mp.mpf("1e-25"), f"n={n}")
        prev = cur
    chk("C = max of the three row constants",
        C_SCHUR == C_ROW1 and C_ROW1 > C_ROW0 and C_ROW1 > C_ROWN)

    # --- (3) the analytic row bound really dominates the exact row sums -----
    #     (this is the check of the DERIVATION of S3.2, not evidence: the exact
    #      row sum is itself computed with a rigorous integral tail)
    for w in (0.2, 0.7, 1.5, 3.0, 8.0, 25.0):
        for n in list(range(0, 12)) + [20, 45, 90, 175]:
            ex = row_sum_exact_ub(n, w)
            an = row_bound_analytic(n, w)
            chk("analytic row bound >= exact row sum", an >= ex,
                f"w={w} n={n} {mp.nstr(an,10)} < {mp.nstr(ex,10)}")
            chk("sharp row bound <= C * S1", an <= C_SCHUR * S1_closed(w) + mp.mpf("1e-25"),
                f"w={w} n={n}")

    # --- (4) wbar_1 and the closure --------------------------------------
    print(f"  wbar_1 = C pi / (2 * frak_g(2)_lo) = {mp.nstr(WBAR1, 10)}   "
          f"(frak_g(2)_lo = {mp.nstr(G2_LO,6)})")
    chk("wbar_1 < WBAR2 (Segment II reaches past it)", WBAR1 < mp.mpf(WBAR2),
        f"{mp.nstr(WBAR1,8)} vs {WBAR2}")
    for w in (WBAR1 * mp.mpf("1.0000001"), mp.mpf(WBAR2), mp.mpf(3), mp.mpf(10),
              mp.mpf(60), mp.mpf(1000)):
        chk("C pi w/2 < g2_lo w^2 for w > wbar_1",
            C_SCHUR * mp.pi * w / 2 < G2_LO * w * w, f"w={mp.nstr(w,8)}")
    chk("equality at wbar_1 (the threshold is exactly C pi/(2 g2_lo))",
        abs(C_SCHUR * mp.pi * WBAR1 / 2 - G2_LO * WBAR1 ** 2) < mp.mpf("1e-25"))

    # --- (5) float64 EVIDENCE: the analytic bound dominates ||H^(1500)|| ----
    print(f"  {'wbar':>7}  {'||H^(1500)||':>12}  {'C*S1':>10}  {'C pi w/2':>10}  "
          f"{'g2_lo w^2':>10}  {'ratio h/(C S1)':>14}")
    for w in (2.1, 2.5, 3.0, 5.0, 10.0, 20.0, 40.0):
        nrm = float(np.max(np.abs(np.linalg.eigvalsh(H_mat(1500, w)))))
        cs1 = C_SCHUR * S1_closed(w)
        cpi = C_SCHUR * mp.pi * mp.mpf(w) / 2
        chk("evidence: ||H^(1500)(w)|| <= C S1(w)", mp.mpf(nrm) <= cs1, f"w={w}")
        chk("evidence: ||H^(1500)(w)|| <= C pi w/2", mp.mpf(nrm) <= cpi, f"w={w}")
        chk("evidence: lam_1(H^(1500)) < g2_lo w^2 on segment III",
            nrm < float(G2_LO) * w * w, f"w={w}")
        print(f"  {w:>7.2f}  {nrm:>12.6f}  {float(cs1):>10.4f}  {float(cpi):>10.4f}  "
              f"{float(G2_LO)*w*w:>10.4f}  {nrm/float(cs1):>14.4f}")
    # the sharper (non-load-bearing) threshold if S1 is kept instead of pi w/2
    lo, hi = mp.mpf("1.0"), mp.mpf(WBAR2)
    for _ in range(200):
        mid = (lo + hi) / 2
        if C_SCHUR * S1_closed(mid) < G2_LO * mid * mid:
            hi = mid
        else:
            lo = mid
    print(f"  (keeping S1 itself instead of pi w/2 would already close from "
          f"wbar = {mp.nstr(hi, 8)}; we quote the closed form wbar_1 = "
          f"{mp.nstr(WBAR1,8)} as the theorem's threshold)")
    _summary("Stage S3", f"({time.time()-t0:.1f}s)")
    return float(WBAR1)


# ===========================================================================
# STAGE S2 -- SEGMENT II: per-cell interval certificate on [0.3, WBAR2].
#
# Per cell [a,b], write H = [[A, B],[B^T, D]] with A = P_Na H P_Na the head:
#     lam_1(H) <= lam_1(diag(A,D)) + ||B||          (Weyl; the off-diagonal part
#                                                    has operator norm ||B||)
#              =  max(lam_1(A), lam_1(D)) + ||B||
#              <= max(c, ||D||) + ||B||             (lam_1(D) <= ||D||)
# NOTE the `max`, not the `c + ||D||` of the base `alpha_cell_ub`: there the head
# term is lam_2(A) ~ 1e-3, comparable to ||D||, so the sum costs nothing; here the
# head term is lam_1(A) ~ h itself, two orders ABOVE ||D||, and keeping the max
# recovers the whole ||D|| (4.97e-3 of a 8.79e-3 total budget at wbar = 0.3).
# with
#   lam_1(A(w)) <= c(w)   CERTIFIED by an mpmath dps=30 Cholesky test on
#                         c I - A(w)  (SPD  ==>  lam_1(A) < c).  NO deflation is
#                         needed: we want the TOP eigenvalue of the head, not the
#                         second.  This is exactly `p4d_kappa_pencil.
#                         lambda1_upper_pencil` at mutilde = 0, imported verbatim.
#   sup_{[a,b]} lam_1(A)  <= max(c(a), c(b)) + L (b-a)/2, because
#                         |lam_1(A(w)) - lam_1(A(p))| <= ||A(w)-A(p)|| <= L|w-p|
#                         (Weyl + the mean value theorem entrywise) and
#                         sup_w min(c_a + L(w-a), c_b + L(b-w))
#                              = (c_a+c_b)/2 + L(b-a)/2 <= max(c_a,c_b)+L(b-a)/2.
#                         L := `schur_full_seg_w0(a,b,1)` is a rigorous,
#                         regress-free sup of ||H'|| on [a,b]; A' = P_Na H' P_Na
#                         is a compression of H', so ||A'|| <= ||H'|| <= L.
#   ||B||, ||D||          the SEGMENT Schur envelopes `schur_coupling` /
#                         `schur_tailblock` evaluated with `spec_seg(a,b,0)`,
#                         hence valid simultaneously for every w in [a,b].
# The cell check is  h_ub < frak_g(2)_lo * a^2  (a = LEFT end, because w^2 is
# increasing, so inf_{[a,b]} frak_g(2) w^2 = frak_g(2) a^2).
# ===========================================================================

def seg2_grid(lo, hi, ratio=BASE_RATIO):
    g = [float(lo)]
    while g[-1] * ratio < hi - 1e-12:
        g.append(g[-1] * ratio)
    g.append(float(hi))
    return g


def seg2_head_size(b):
    """Head size for the Segment-II lam_1 certificate, keyed on the cell's RIGHT
    end.  The base `p4_interval_extension.head_size` bands (70/140/280) are keyed
    to the ALPHA certificate, whose budget is set by lam_2 ~ 1e-3; here the
    budget is the gap frak_g(2) a^2 - h(a), which is only 7.4% at wbar = 0.3 and
    shrinks like 1.19 a^2 as a -> 0, so the small-wbar end needs the bigger head.
    Na = 140 below 1.2 (slack 2.33e-3 at wbar = 0.3, vs 2.6e-4 at Na = 70) and
    the base Na = 70 above it, where the margin is >= 25%."""
    return 140 if b <= 1.2 + 1e-12 else 70


def h_cell_ub(a, b, cA, cB, L, Na):
    """Certified sup_{w in [a,b]} lam_1(H(w))."""
    sp = spec_seg(a, b, 0)
    nB = schur_coupling(Na, b, 0, spec=sp)
    nD = schur_tailblock(Na, b, 0, spec=sp)
    head = max(cA, cB) + mp.mpf(L) * (mp.mpf(b) - mp.mpf(a)) / 2
    return max(head, nD) + nB + ROUND_PAD, nB, nD


SEG2 = {}


def run_S2():
    print(f"=== STAGE S2: SEGMENT II -- interval certificate on [0.3, {WBAR2}] ===")
    t0 = time.time()
    grid = seg2_grid(0.3, WBAR2)
    nas = [seg2_head_size(grid[min(i + 1, len(grid) - 1)]) for i in range(len(grid))]
    print(f"  grid: {len(grid)} geometric points, ratio = BASE_RATIO = {BASE_RATIO:.6f} "
          f"(the ratio the base interval theorem certified on [0.3,5]); "
          f"head Na = 140 up to wbar 1.2, then the base 70")

    ccache = {}

    def c_at(w, Na):
        key = (round(float(w), 12), Na)
        if key not in ccache:
            _, ok, c, _, _ = lambda1_upper_pencil(float(w), 0.0, 0, Na=Na)
            chk("non-deflated Cholesky head certificate SPD", ok, f"w={w:.5f} Na={Na}")
            chk("head bound c > 0", c > 0, f"w={w:.5f}")
            ccache[key] = c
        return ccache[key]

    rows = []
    for i in range(len(grid) - 1):
        a, b = grid[i], grid[i + 1]
        Na = min(nas[i], nas[i + 1])
        # c certified for the Na-head is NOT automatically valid for a LARGER
        # head (Cauchy interlacing makes lam_1 of a bigger head larger), so both
        # endpoints are certified at the SAME Na -- the smaller of the two bands.
        cA, cB = c_at(a, Na), c_at(b, Na)
        Lw = schur_full_seg_w0(a, b, 1)
        Lb = schur_full_seg(a, b, 1)
        chk("weighted-row-0 L <= base L (subset of the same derivation)",
            Lw <= Lb + 1e-15, f"cell [{a:.4f},{b:.4f}]")
        hub, nB, nD = h_cell_ub(a, b, cA, cB, Lw, Na)
        tgt = G2_LO * mp.mpf(a) ** 2
        slack = tgt - hub
        chk("cell: h_ub < g2_lo * a^2", slack > 0,
            f"cell [{a:.5f},{b:.5f}] h_ub={mp.nstr(hub,8)} tgt={mp.nstr(tgt,8)}")
        rows.append((a, b, Lw, float(hub), float(tgt), float(slack),
                     float(hub / tgt), float(nB), float(nD)))

    worst = min(rows, key=lambda r: r[5])
    worst_rel = max(rows, key=lambda r: r[6])
    SEG2["rows"] = rows
    SEG2["inf_slack"] = worst[5]
    SEG2["hi"] = WBAR2

    print(f"  {'segment':>16}  {'cells':>6}  {'L max':>7}  {'min(tgt-h_ub)':>14}  "
          f"{'max h_ub/tgt':>13}")
    bands = [(0.3, 0.5), (0.5, 0.7), (0.7, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, WBAR2)]
    for lo, hi in bands:
        sel = [r for r in rows if r[0] >= lo - 1e-9 and r[1] <= hi + 1e-9]
        if not sel:
            continue
        print(f"  [{lo:5.3f},{hi:5.3f}]  {len(sel):>6d}  {max(r[2] for r in sel):>7.3f}  "
              f"{min(r[5] for r in sel):>14.6f}  {max(r[6] for r in sel):>13.4f}")
    print(f"  binding cell (smallest absolute slack): [{worst[0]:.5f},{worst[1]:.5f}]  "
          f"L={worst[2]:.4f}  h_ub={worst[3]:.6f}  g2_lo*a^2={worst[4]:.6f}  "
          f"slack={worst[5]:.6f}  (||B||={worst[7]:.5f}, ||D||={worst[8]:.5f})")
    print(f"  binding cell (largest ratio)          : [{worst_rel[0]:.5f},"
          f"{worst_rel[1]:.5f}]  h_ub/(g2_lo a^2) = {worst_rel[6]:.6f}")
    print(f"  ==> inf over cells of (frak_g(2)_lo a^2 - h_ub) = {worst[5]:.6f} > 0 "
          f"on [0.3, {WBAR2}]")
    chk("Segment II infimum of the slack is positive", worst[5] > 0)
    chk("Segment II reaches past wbar_1 (overlap with Segment III)",
        mp.mpf(WBAR2) > WBAR1)

    # --- float64 EVIDENCE: the certificate is not vacuous ------------------
    step = max(1, len(rows) // 14)
    for r in rows[::step] + [rows[-1]]:
        a, b = r[0], r[1]
        mid = 0.5 * (a + b)
        l1 = float(np.linalg.eigvalsh(H_mat(1500, mid))[-1])
        chk("evidence: h_ub >= lam_1(H^(1500)) at the cell midpoint", r[3] >= l1,
            f"cell [{a:.4f},{b:.4f}] h_ub={r[3]:.6f} lam1={l1:.6f}")
        chk("evidence: lam_1(H^(1500)) < g2_lo mid^2 (the true margin)",
            l1 < float(G2_LO) * mid * mid, f"mid={mid:.4f}")
    _summary("Stage S2", f"({time.time()-t0:.1f}s)")


# ===========================================================================
# STAGE M -- Sec. 3: MONOTONICITY of h on [sqrt(2), W].
#
# M.1  HELLMANN-FEYNMAN, justified.
#   (i)  w -> H(w) is C^1 in OPERATOR NORM on (0,inf): the entries are rational
#        in w, and by the entrywise Taylor remainder
#          |H_nm(w+d) - H_nm(w) - d H'_nm(w)| <= (d^2/2) sup |H''_nm|,
#        so  ||H(w+d) - H(w) - d H'(w)|| <= (d^2/2) * schur_full_seg_w0(.,.,2),
#        because the Schur test applied to the ENVELOPE matrix of |H''| bounds
#        the operator norm of any matrix it dominates entrywise.  The same
#        estimate makes H' itself norm-continuous.
#   (ii) On a cell where the certificate below produces g > 0, lam_1(H(w)) is a
#        SIMPLE ISOLATED eigenvalue (lam_2 <= alpha < theta <= lam_1), so the
#        Riesz projection P(w) = (2 pi i)^{-1} \oint (z-H)^{-1} dz is C^1, hence
#        so are h(w) = tr(H P)/tr(P) and a normalized eigenvector v(w) chosen
#        with <v,v'> = 0 (Kato, Ch.II Sec.5 / Rellich).  Differentiating
#        H v = h v and pairing with v gives  h'(w) = <v(w), H'(w) v(w)>.
#   (iii) v' = (h - H)^{+} H' v with (h-H)^{+} the reduced resolvent on {v}^perp,
#        whose norm is 1/dist(lam_1, spec \ {lam_1}) = 1/gap.  Hence
#        ||v'|| <= ||H' v||/gap <= ||H'||/gap.
#
# M.2  THE PER-CELL LOWER BOUND.  Let p = (a+b)/2, hw = (b-a)/2, and let v_hat be
#   the float64 top eigenvector of H^(Ns)(p), zero-padded (FINITE SUPPORT).
#     q_hat := <v_hat, H'(p) v_hat>/||v_hat||^2   -- EXACT in mpmath: finite
#              support means only the Ns x Ns head of H'(p) is touched.
#     theta := <v_hat, H(p) v_hat>/||v_hat||^2 <= lam_1(H(p))   (Rayleigh)
#     r     := H(p) v_hat/||v_hat|| - theta v_hat/||v_hat||.  Its rows m < Nr are
#              computed EXACTLY (mpmath on the float64 rectangular block
#              H_rect(Nr,Ns,p)); rows m >= Nr are the block (m >= Nr, n < Ns) of H
#              applied to v_hat, bounded by the NEW envelope
#                  ||(H)_{m>=Nr, n<Ns}|| <= sqrt(rowmax * colmax),
#                  colmax = max_{n<Ns} col_env(n, Nr, spec)          [IMPORTED]
#                  rowmax <= 2 psi(Nr-Ns+1) sqrt(2Ns-1)/sqrt(2Nr+1)  [derived:
#                     sum_{n<Ns} (psi(m-n)+psi(m+n+1))/sqrt((2m+1)(2n+1))
#                     <= (2 psi(m-Ns+1)/sqrt(2m+1)) sum_{n<Ns}(2n+1)^{-1/2}
#                     <= 2 psi(m-Ns+1) sqrt(2Ns-1)/sqrt(2m+1),  decreasing in m,
#                     so its max over m >= Nr is at m = Nr].
#              THIS TAIL TERM IS MANDATORY -- omitting it is exactly the bug the
#              main session found in the first P4-D kappa certificate
#              (`P4D_coulomb_operator.md` Sec.2.4b: head-only residual 1e-16 vs
#              true full residual 1e-4).
#     Davis-Kahan sinTheta:  sin angle(v_hat, v(p)) <= ||r|| / gamma,
#              gamma := theta - alpha_cell  (alpha_cell >= lam_2 >= sup of the
#              rest of the spectrum, so dist(theta, spec \ {lam_1}) >= gamma).
#     |<v,H'v> - q_hat| <= 2 sqrt2 ||H'|| sin angle   (||v - v_hat|| <= sqrt2 sin).
#     Lipschitz over the cell:  Phi(w) := <v(w),H'(w)v(w)>,
#              Phi' = 2<v',H'v> + <v,H''v>,  |Phi'| <= 2 L1^2/g + L2,
#              g := a cell-wide lower bound on the gap:
#                   g = max( theta - L1 hw - alpha_cell ,  mu_2 - alpha_cell ).
#   ==>  for every w in [a,b],
#          h'(w) >= q_hat - 2 sqrt2 L1 ||r||/gamma - (2 L1^2/g + L2) hw - 2 ROUND_PAD.
#   L1 = schur_full_seg_w0(a,b,1), L2 = schur_full_seg_w0(a,b,2) (both asserted
#   <= their base `schur_full_seg` twins), alpha_cell = the Segment-II-style
#   segment upper bound on lam_2 built from `alpha_certificate_H` (WITH deflation
#   -- here we do want the second eigenvalue) on a COARSE grid; a coarse-cell
#   alpha/L1/L2 is valid on every sub-cell because spec_seg(a',b',j) is
#   dominated by spec_seg(a,b,j) entrywise when [a',b'] subset [a,b].
# ===========================================================================

COARSE_RATIO = 1.05        # alpha / L1 / L2 grid (expensive Cholesky)
MAX_SPLIT_DOUBLINGS = 5    # hard cap on the adaptive refinement


def coupling_restricted(Ns, Nr, spec):
    """||(H^{(0)})_{m >= Nr, n < Ns}||_op, Ns <= Nr.  See M.2 for the derivation."""
    cmax = mp.mpf(0)
    for n in range(Ns):
        c = col_env(n, Nr, spec)
        if c > cmax:
            cmax = c
    rmax = 2 * psi_spec(Nr - Ns + 1, spec) * mp.sqrt(2 * Ns - 1) / mp.sqrt(2 * Nr + 1)
    return mp.sqrt(rmax * cmax)


def support_size(w):
    """Trial-vector support: big enough that the true top eigenvector's tail
    beyond Ns (which sets the Davis-Kahan residual) is small; empirically the
    residual scales like w/Ns, so Ns ~ 10 w."""
    return max(140, int(10 * float(w)) // 2 * 2)


def alpha_cell_seg(a, b, cA, cB, L, Na):
    """Segment upper bound on lam_2(H(w)) for w in [a,b] -- `piu.alpha_cell_ub`
    with an explicit head size (identical algebra)."""
    sp = spec_seg(a, b, 0)
    nB = schur_coupling(Na, b, 0, spec=sp)
    nD = schur_tailblock(Na, b, 0, spec=sp)
    return max(cA, cB) + mp.mpf(L) * (mp.mpf(b) - mp.mpf(a)) / 2 + nB + nD + ROUND_PAD


def hprime_probe(p, Ns):
    """float64 ESTIMATE of (h'(p), lam_1(p)) -- used ONLY to choose the
    refinement depth.  Every cell it proposes is verified afterwards by
    `hprime_cell_lb`, so nothing it returns enters a bound."""
    A = H_mat(Ns, p)
    e, V = np.linalg.eigh(A)
    v = V[:, -1]
    return float(v @ H_mat(Ns, p, 1) @ v), float(e[-1])


def hprime_cell_lb(a, b, L1, L2, alpha, mu2):
    """Certified lower bound on h'(w) valid for EVERY w in [a,b] (M.2)."""
    p = 0.5 * (a + b)
    hw = mp.mpf(b - a) / 2
    Ns = support_size(p)
    Nr = 4 * Ns
    A_head = H_mat(Ns, p, 0)
    v = np.linalg.eigh(A_head)[1][:, -1]
    vm = _mp_vec(v)
    nv2 = sum(x * x for x in vm)
    nv = mp.sqrt(nv2)

    Rblk = H_rect(Nr, Ns, p, 0)
    y = _mp_matvec(Rblk, vm)
    theta = sum(vm[i] * y[i] for i in range(Ns)) / nv2
    r2 = mp.zero
    for i in range(Ns):
        d = (y[i] - theta * vm[i]) / nv
        r2 += d * d
    for i in range(Ns, Nr):
        d = y[i] / nv
        r2 += d * d
    rtail = coupling_restricted(Ns, Nr, spec_point(p, 0))
    rnorm = mp.sqrt(r2 + rtail * rtail) + ROUND_PAD

    qhat = _mp_quad(H_mat(Ns, p, 1), vm) / nv2 - ROUND_PAD
    theta_lo = theta - ROUND_PAD

    gamma = theta_lo - alpha
    g = max(theta_lo - mp.mpf(L1) * hw - alpha, mu2 - alpha)
    if gamma <= 0 or g <= 0:
        return None, dict(reason="no certified gap", gamma=float(gamma), g=float(g))
    sin_b = rnorm / gamma
    dk = 2 * mp.sqrt(2) * mp.mpf(L1) * sin_b
    lipr = 2 * mp.mpf(L1) ** 2 / g + mp.mpf(L2)
    lip = lipr * hw
    lb = qhat - dk - lip
    return lb, dict(p=p, Ns=Ns, Nr=Nr, qhat=float(qhat), theta=float(theta_lo),
                    alpha=float(alpha), gamma=float(gamma), g=float(g),
                    rnorm=float(rnorm), rtail=float(rtail), sin=float(sin_b),
                    dk=float(dk), lip=float(lip), liprate=float(lipr), hw=float(hw))


MONO = {}


def run_M(W=30.0):
    print(f"=== STAGE M: Sec.3 -- h'(wbar) > 0 certified on [sqrt2, {W}] ===")
    t0 = time.time()
    lo = float(SQRT2)
    coarse = seg2_grid(lo, W, COARSE_RATIO)
    print(f"  coarse alpha/L grid: {len(coarse)} points, ratio {COARSE_RATIO} "
          f"(Cholesky-deflation alpha, banded Na); fine cells adaptive "
          f"(cap {MAX_SPLIT_DOUBLINGS} doublings)")

    ca = []
    for w in coarse:
        Na = head_size(w)
        al, ok, c, _, _ = alpha_certificate_H(float(w), Na=Na)
        chk("Cholesky-deflation alpha certificate SPD (coarse grid)", ok, f"w={w:.4f}")
        ca.append((c, Na))

    fine_rows = []
    worst = None
    n_fail = 0
    stop_at = None
    for i in range(len(coarse) - 1):
        A, B = coarse[i], coarse[i + 1]
        Na = min(ca[i][1], ca[i + 1][1])
        cA, cB = ca[i][0], ca[i + 1][0]
        if ca[i][1] != ca[i + 1][1]:
            _, okA, cA, _, _ = alpha_certificate_H(A, Na=Na)
            _, okB, cB, _, _ = alpha_certificate_H(B, Na=Na)
            chk("band-boundary alpha recertification SPD", okA and okB, f"[{A},{B}]")
        L1 = schur_full_seg_w0(A, B, 1)
        L2 = schur_full_seg_w0(A, B, 2)
        chk("L1 (weighted row 0) <= base schur_full_seg order 1",
            L1 <= schur_full_seg(A, B, 1) + 1e-15, f"[{A:.4f},{B:.4f}]")
        chk("L2 (weighted row 0) <= base schur_full_seg order 2",
            L2 <= schur_full_seg(A, B, 2) + 1e-15, f"[{A:.4f},{B:.4f}]")
        alpha = alpha_cell_seg(A, B, cA, cB, L1, Na)
        mu2, k1, k2, k3 = mu2_cell_lb(A, B)
        chk("mu_2 segment monotonicity k1,k2,k3 >= 0", k1 >= 0 and k2 >= 0 and k3 >= 0,
            f"[{A:.4f},{B:.4f}]")

        # predict the number of sub-cells from a float64 probe, then VERIFY
        pmid = 0.5 * (A + B)
        est, lam1_est = hprime_probe(pmid, support_size(pmid))
        g_pred = max(float(mu2 - alpha),
                     lam1_est - L1 * (B - A) / 2 - float(alpha), 1e-9)
        rate = 2.0 * L1 ** 2 / g_pred + L2
        k = max(1, int(np.ceil((B - A) * rate / (2 * max(est, 1e-6) * 0.45))))
        ok_cell = False
        for _ in range(MAX_SPLIT_DOUBLINGS + 1):
            sub = np.linspace(A, B, k + 1)
            got = []
            bad = False
            for j in range(k):
                lb, info = hprime_cell_lb(float(sub[j]), float(sub[j + 1]),
                                          L1, L2, alpha, mu2)
                if lb is None or lb <= 0:
                    bad = True
                    break
                got.append((float(sub[j]), float(sub[j + 1]), float(lb), info))
            if not bad:
                ok_cell = True
                fine_rows.extend(got)
                for gr in got:
                    chk("fine cell: certified h'(w) > 0 for all w in the cell",
                        gr[2] > 0, f"[{gr[0]:.5f},{gr[1]:.5f}]")
                    if worst is None or gr[2] < worst[2]:
                        worst = gr
                break
            k *= 2
        if not ok_cell:
            n_fail += 1
            stop_at = A
            print(f"  !! coarse cell [{A:.4f},{B:.4f}] did NOT close after "
                  f"{MAX_SPLIT_DOUBLINGS} doublings (k={k//2}); stopping here.")
            break

    Wdone = fine_rows[-1][1] if fine_rows else lo
    MONO["W"] = Wdone
    MONO["cells"] = len(fine_rows)
    MONO["worst"] = worst

    print(f"  {'band':>14}  {'cells':>6}  {'min hp_lb':>10}  {'q_hat':>8}  "
          f"{'DK term':>9}  {'Lip term':>9}  {'Ns':>5}")
    edges = [lo, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0, W + 1e-9]
    for e0, e1 in zip(edges[:-1], edges[1:]):
        sel = [r for r in fine_rows if r[0] >= e0 - 1e-9 and r[1] <= e1 + 1e-9]
        if not sel:
            continue
        wc = min(sel, key=lambda r: r[2])
        print(f"  [{e0:5.3f},{e1:6.3f}]  {len(sel):>6d}  {wc[2]:>10.5f}  "
              f"{wc[3]['qhat']:>8.5f}  {wc[3]['dk']:>9.5f}  {wc[3]['lip']:>9.5f}  "
              f"{wc[3]['Ns']:>5d}")
    if worst is not None:
        iw = worst[3]
        bind = max((("Davis-Kahan angle", iw['dk']), ("Lipschitz slack", iw['lip'])),
                   key=lambda t: t[1])
        print(f"  binding fine cell [{worst[0]:.5f},{worst[1]:.5f}]: "
              f"h'_lb = {worst[2]:.6f}, q_hat = {iw['qhat']:.6f}, "
              f"DK = {iw['dk']:.6f} (sin<={iw['sin']:.3e}, ||r||={iw['rnorm']:.3e}, "
              f"tail={iw['rtail']:.3e}, gamma={iw['gamma']:.4f}), "
              f"Lip = {iw['lip']:.6f} (rate {iw['liprate']:.3f}, hw={iw['hw']:.5f})")
        print(f"  ==> the LARGER slack term on the binding cell is: {bind[0]}")
    print(f"  ==> h'(wbar) > 0 CERTIFIED on [sqrt2, {Wdone:.4f}]  "
          f"({len(fine_rows)} fine cells)")
    chk("monotonicity interval is non-trivial", Wdone > 2.0, f"W={Wdone}")

    # --- float64 EVIDENCE: central differences of lam_1 at N=1500 ---------
    for w in (1.5, 2.5, 4.0, 7.0, 12.0, 18.0):
        if w > Wdone:
            continue
        d = 1e-3
        fd = (float(np.linalg.eigvalsh(H_mat(1500, w + d))[-1])
              - float(np.linalg.eigvalsh(H_mat(1500, w - d))[-1])) / (2 * d)
        near = [r for r in fine_rows if r[0] <= w <= r[1]]
        chk("evidence: finite-difference h'(w) > 0 at N=1500", fd > 0, f"w={w}")
        if near:
            chk("evidence: certified lower bound <= finite-difference h'",
                near[0][2] <= fd + 1e-9, f"w={w} lb={near[0][2]:.6f} fd={fd:.6f}")
        # Hellmann-Feynman itself, as float64 evidence
        A = H_mat(1500, w)
        v = np.linalg.eigh(A)[1][:, -1]
        hf = float(v @ H_mat(1500, w, 1) @ v)
        chk("evidence: <v,H'v> == d lam_1/dw at N=1500", abs(hf - fd) < 1e-5,
            f"w={w} hf={hf:.8f} fd={fd:.8f}")
    _summary("Stage M", f"({time.time()-t0:.1f}s)")


def run_verdict():
    print("\n=== VERDICT ===")
    print(f"  Segment II  : h(wbar) < frak_g(2) wbar^2 CERTIFIED on [0.3, {WBAR2}] "
          f"({len(SEG2['rows'])} cells, inf slack {SEG2['inf_slack']:.6f})")
    print(f"  Segment III : ||H(wbar)|| <= C S1(wbar) <= C pi wbar/2, "
          f"C = {mp.nstr(C_SCHUR,10)}  ==>  closed for wbar > wbar_1 = "
          f"{mp.nstr(WBAR1,8)}")
    print(f"  overlap     : [{mp.nstr(WBAR1,8)}, {WBAR2}] is covered by BOTH "
          f"segments, so Lemma B holds on [0.3, inf) once Segment I closes "
          f"(0, wbar_0] with wbar_0 >= 0.3.")
    if MONO:
        print(f"  Sec.3       : h'(wbar) > 0 CERTIFIED on [sqrt2, {MONO['W']:.4f}] "
              f"({MONO['cells']} fine cells)")


if __name__ == "__main__":
    want = [s.upper() for s in sys.argv[1:]] or ["S0", "S3", "S2"]
    if "ALL" in want:
        want = ["S0", "S3", "S2", "M"]
    T0 = time.time()
    if "S0" in want:
        run_S0()
    if "S3" in want:
        run_S3()
    if "S2" in want:
        run_S2()
    if "M" in want:
        W = float(os.environ.get("CONJ2_W", "30.0"))
        if "S0" not in want:
            run_S0()
        run_M(W)
    if {"S2", "S3"} <= set(want):
        run_verdict()
    print(f"\nALL {PASS[0]} checks passed  ({time.time()-T0:.1f}s).")
