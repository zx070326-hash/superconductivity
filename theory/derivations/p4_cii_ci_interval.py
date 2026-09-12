r"""
P4 . EINSTEIN EXTREMAL -- Segment C-I, INTERVAL-ARITHMETIC CERTIFICATE.

    H(x) := h(sqrt x) = lam_1(H(sqrt x))  is STRICTLY CONCAVE on (0, X],
    certified cell by cell with mpmath.iv (OUTWARD-ROUNDED INTERVAL ARITHMETIC),
    dps = 30.

RUN
    python theory/derivations/p4_cii_ci_interval.py              gate + Stage S + ladder to 3
    python theory/derivations/p4_cii_ci_interval.py EXT          ladder to 45
    python theory/derivations/p4_cii_ci_interval.py EXT 120      ladder to 120
    python theory/derivations/p4_cii_ci_interval.py EXT 45 -j 22 22-process pool
    python theory/derivations/p4_cii_ci_interval.py GATE         Stage A' only

SCOPE AND RELATION TO `p4_extremal_concavity.py`
    That file OWNS the mathematics of Segment C-I: the derivation in its
    docstring Sec.1--Sec.2 ((P2)/(P3), the identities (I1)-(I3), the
    cancellation-free form (C**), the per-cell machine).  NOTHING there is
    re-derived here and NOTHING there is modified.  This file re-implements the
    same certificate so that EVERY inequality is supported by outward-rounded
    interval operations, closing residual assumption #1 of
    `P4_main_theorems.md` Sec.5 for Segment C-I:

        "mpmath dps >= 30, explicit rounding pads (1e-12 / 1e-25), Cholesky
         success treated as a positive-definiteness certificate; NON-INTERVAL
         ARITHMETIC."

    The plan, the assertion inventory and the cost analysis are in
    `theory/notes/2026-09-12_ci_intervalization_plan.md`.

WHAT "INTERVAL" MEANS HERE (the three things that were not true before)

 (a) NO ROUNDING PADS.  `ROUND_PAD = 1e-12` and the `+1e-25` tolerances are
     GONE; `PAD` below is identically zero and is asserted to be zero.  Their
     role -- absorbing the (unquantified) floating-point error of the mpf
     evaluation -- is taken over by outward rounding, which is quantified by
     construction.  Precisely: the ONLY arithmetic primitives used on certified
     quantities are  + - * /  , integer powers, sqrt and min/max, all taken from
     `mpmath.iv`, each of which returns an interval containing the exact real
     result (mpmath's `iv` implements these by directed rounding of the endpoint
     mpf operations).  ONE transcendental is needed (arctan, for the closed-form
     tail integral) and it is NOT taken from mpmath: `iv.atan2` is not outward
     rounded in mpmath 1.3.0 (`mpf_atan2` evaluates `mpf_atan(mpf_div(y,x))` at
     prec+4 with round-to-NEAREST and only directs the LAST rounding, so the
     returned interval can fail to enclose -- measured: ~9% of random operands,
     worst 0.24 ulp at prec 103).  `iv_atan_nn`/`iv_pi_half` below rebuild it
     from +,-,*,/,sqrt and integer powers only, with an explicit alternating-
     series remainder, and are cross-checked in Stage A' against three rigorous
     identities (atan t + atan 1/t = pi/2, the doubling formula, and Machin).

 (b) POSITIVE DEFINITENESS IS CERTIFIED, NOT ASSUMED.  The old chain ran
     `mp.cholesky` on an mpf copy of a FLOAT64 head block and treated "no
     exception" as "SPD".  Two gaps: the arithmetic was inexact, and the matrix
     was the float64 one, not the exact operator's head.  Here:

     * the head block is built in `iv` directly from the exact kernel (integers
       k and the exactly-representable binary cell endpoints), so the interval
       matrix [A] provably contains the exact head;
     * positive definiteness of every member of [cI - PAP] is certified by an
       interval Cholesky, justified by:

       THEOREM (interval Cholesky => uniform positive definiteness;
       Alefeld 1980, see also Rump, BIT 46 (2006) 433-452, Sec.2).
       Let [S] be a symmetric interval matrix and run Cholesky with interval
       arithmetic (every operation outward-rounded).  If the algorithm
       completes and every pivot interval [d_i] satisfies inf[d_i] > 0, then
       EVERY symmetric real S in [S] is positive definite.

       PROOF.  Fix a symmetric S in [S] and run EXACT Cholesky on S, producing
       L^S as far as it is defined.  Claim: for every i, the exact partial
       quantities satisfy L^S_{ij} in [L_{ij}] for all j <= i.  Induction on i.
       Suppose it holds for all rows < i.  For j < i the exact off-diagonal
       entry is L^S_{ij} = (S_{ij} - sum_{k<j} L^S_{ik} L^S_{jk}) / L^S_{jj};
       by the inductive hypothesis and the INCLUSION PROPERTY of interval
       arithmetic (an interval operation on intervals containing the exact
       operands returns an interval containing the exact result), this lies in
       the interval computed by the algorithm, provided the divisor interval
       [L_{jj}] does not contain 0 -- which holds because inf[L_{jj}] =
       sqrt(inf[d_j]) > 0.  For j = i the exact radicand is
       d^S_i = S_{ii} - sum_{k<i} (L^S_{ik})^2, which lies in [d_i], hence
       d^S_i >= inf[d_i] > 0: the exact Cholesky step is defined and its pivot
       is POSITIVE.  Induction closes.  Therefore exact Cholesky runs to
       completion on S with all pivots positive, i.e. S = L^S (L^S)^T with L^S
       lower triangular and invertible, so S is positive definite.  []

       (Only the inclusion property and Sylvester's criterion are used; no
       reference to the numerical stability of Cholesky is needed, which is
       exactly why this replaces the old "no exception was raised" argument.)

 (c) EXACT GRID.  All cell endpoints are integer multiples of 2^-16 and are
     verified to be so (`_dyadic_ok`), so consecutive cells abut EXACTLY and the
     union of the certified cells provably covers (0, WMAX].  The float64
     `eigh` that produces trial vectors and the float64 heuristic that CHOOSES
     the step size are NOT bounds: any fixed vector gives a valid Rayleigh /
     Davis-Kahan bound, and any exactly-representable tiling is a valid tiling.
     Both are re-certified downstream in `iv`.

ONE-SIDED BOUND DISCIPLINE (plan Sec.5).  Every elementary operation happens in
`iv` (see (a) for the exact primitive list); every NAMED certified quantity is
then collapsed to a single mpf endpoint
(`LB` or `UB`) with the monotonicity that justifies the choice stated inline.
This prevents the two classic holes ("an intermediate width blew up but a later
max swallowed it" -- impossible when only one endpoint flows; "an iv value was
compared to an iv value with an ambiguous outcome" -- impossible when the
comparison operands are mpf endpoints).  Reference:

    E_hi = ehat_psi_hi + s*(N_hi + |qhat_psi|_hi)     increasing in all  -> UB
    F_lo = qhat_phi_lo - s^2*(a_phi_hi + N_hi) - 2*s*ehat_phi_hi
                                                     decreasing in s,a,N,ehat;
                                                     increasing in qhat_lo -> LB
    g_lo = max(theta_lo - L1_hi*hw, mu2_lo) - alpha_hi
    pos  = b^4 * E_hi^2 / g_lo                       increasing E_hi, decr g_lo
    hpp_ub = -2*F_lo + 2*pos     <  0   is the CONCLUSION.

REMAINING NON-INTERVAL INGREDIENTS (honest list; none is a bound)
  1. `numpy.linalg.eigh` on a float64 matrix -> the trial vector vhat and the
     head deflation vector u.  Legitimate: for ANY fixed vector the Rayleigh
     quotient is a lower bound on lam_1 and `lam_2(A) <= lam_max(PAP)` holds for
     ANY unit u (see `alpha_head_iv`).  Their float64 entries are treated as the
     exact binary rationals they are, and normalisation happens inside `iv`.
  2. The choice of the coarse step, of the number of fine subdivisions, and of
     Ns/Na (`support_size`, `na_band`).  Legitimate: these select WHICH
     inequalities to certify, not whether they hold.
  3. Stage A'-only float64 cross-checks against `p4_extremal_concavity.py`
     (typo insurance, explicitly labelled EVIDENCE, never used in a bound).
"""

import os
import sys
import time

import numpy as np
from mpmath import mp, iv

mp.dps = 30
iv.dps = 30

PASS = [0]
_LAST = [0]
T0 = time.time()


def chk(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name} {detail}")
    PASS[0] += 1


def _summary(tag, note=""):
    k = PASS[0] - _LAST[0]
    _LAST[0] = PASS[0]
    print(f"  [{tag}] {k} checks passed. {note}")


# ===========================================================================
# SEC. 0.  INTERVAL PRIMITIVES
# ===========================================================================

IVT = type(iv.mpf(0))
PAD = mp.zero                    # the old ROUND_PAD; MUST stay exactly zero
assert PAD == 0


def ivf(x):
    """Exact point interval.  int / float / mpf / decimal string -> [x, x] when x
    is exactly representable in binary (int, float, mpf: always), else the
    tightest enclosing interval (decimal strings)."""
    if isinstance(x, IVT):
        return x
    return iv.mpf(x)


def ivlh(lo, hi):
    return iv.mpf([lo, hi])


def LB(x):
    """Exact mpf lower endpoint of an interval."""
    return mp.make_mpf(x._mpi_[0])


def UB(x):
    """Exact mpf upper endpoint of an interval."""
    return mp.make_mpf(x._mpi_[1])


IZERO = iv.mpf(0)
IONE = iv.mpf(1)


def ivdot(xs, ys):
    """sum_i xs[i]*ys[i] in interval arithmetic (outward rounded each step)."""
    s = IZERO
    for x, y in zip(xs, ys):
        s = s + x * y
    return s


def ivsumsq(xs):
    s = IZERO
    for x in xs:
        s = s + x ** 2
    return s


def ivsum(xs):
    s = IZERO
    for x in xs:
        s = s + x
    return s


def iv_sqrt_nn(x):
    """sqrt of a quantity KNOWN to be >= 0.  If the interval straddles 0 (which
    can only happen through outward rounding of a genuinely nonnegative
    quantity) the lower endpoint is clamped to 0 -- still an enclosure of the
    exact value, since the exact value is >= 0 by hypothesis."""
    lo, hi = LB(x), UB(x)
    if hi < 0:
        raise AssertionError(f"iv_sqrt_nn of a negative interval [{lo}, {hi}]")
    if lo < 0:
        lo = mp.zero
    return iv.sqrt(ivlh(lo, hi))


def absmax_ub(x):
    """UB(|x|) for an interval x, i.e. max(|inf|, |sup|) -- a valid upper bound
    on |t| for every t in x, whatever the sign of x."""
    return max(abs(LB(x)), abs(UB(x)))


def ub_pos(x):
    """UB(max(x, 0))."""
    u = UB(x)
    return u if u > 0 else mp.zero


def as_tuple(x):
    """Picklable exact representation of an mpf (for multiprocessing)."""
    return x._mpf_


def from_tuple(t):
    return mp.make_mpf(t)


def sub_lo(x, y):
    """LB(x - y) for mpf x, y -- OUTWARD-rounded.  Plain mpf subtraction rounds
    to NEAREST and would silently produce a value slightly ABOVE the exact
    difference, i.e. an invalid lower bound; every mpf-level combination of two
    already-certified endpoints in this file goes through sub_lo/add_hi/mul_hi/
    div_hi so that the one-sidedness survives the last arithmetic step too."""
    return LB(ivf(x) - ivf(y))


def add_hi(x, y):
    return UB(ivf(x) + ivf(y))


def mul_hi(x, y):
    return UB(ivf(x) * ivf(y))


def div_hi(x, y):
    return UB(ivf(x) / ivf(y))


# ===========================================================================
# SEC. 1.  KERNEL TABLES  (built in iv from integers k and exact binary w)
# ===========================================================================
# index 0 is never used by O[.] (see p4_extremal_concavity.py Sec.1.3); it is
# set to 0 so that c_0 = 0 multiplies a finite number.

_K2 = {}


def k2_table(kmax):
    """[k^2] for k = 0..kmax as EXACT point intervals (integers up to 2^51 are
    exact in binary, and kmax <= 1e4 here)."""
    t = _K2.get("t")
    if t is not None and len(t) > kmax:
        return t
    t = [ivf(k * k) for k in range(kmax + 1)]
    _K2["t"] = t
    return t


def tab_H_iv(w, kmax):
    """f_w(k) = w^2/(k^2+w^2)."""
    w2 = ivf(w) ** 2
    k2 = k2_table(kmax)
    return [IZERO] + [w2 / (k2[k] + w2) for k in range(1, kmax + 1)]


def tab_phi_iv(w, kmax):
    """phi_x(k) = k^2/(k^2+w^2)^3,  x = w^2  --  A_xx = -2 O[phi_x]."""
    w2 = ivf(w) ** 2
    k2 = k2_table(kmax)
    return [IZERO] + [k2[k] / (k2[k] + w2) ** 3 for k in range(1, kmax + 1)]


def tab_psi_iv(w, kmax):
    """psi_x(k) = 1/(k^2+w^2)^2."""
    w2 = ivf(w) ** 2
    k2 = k2_table(kmax)
    return [IZERO] + [IONE / (k2[k] + w2) ** 2 for k in range(1, kmax + 1)]


def tab_kappa_iv(w, kmax):
    """kappa_w(k) = 1/(k^2(k^2+w^2))."""
    w2 = ivf(w) ** 2
    k2 = k2_table(kmax)
    return [IZERO] + [IONE / (k2[k] * (k2[k] + w2)) for k in range(1, kmax + 1)]


def tab_gamma_iv(gamma, kmax):
    """k^{-gamma} for INTEGER gamma (exact power, no log/exp)."""
    g = int(gamma)
    assert g == gamma
    return [IZERO] + [IONE / (ivf(k) ** g) for k in range(1, kmax + 1)]


def cell_table(ft_b, ft_a):
    """[f_b(k), f_a(k)] -- the cell-uniform enclosure of f_w(k) for w in [a,b].

    Every kernel used here (f_w, phi_w, psi_w, kappa_w) is NON-INCREASING in w
    pointwise in k, so f_b(k) <= f_w(k) <= f_a(k) for a <= w <= b.  Feeding this
    interval table to the structure-coefficient contraction reproduces EXACTLY
    the sign split of p4_extremal_concavity.py Sec.2.1: interval multiplication
    of a (possibly sign-straddling) c_k by [f_b, f_a] yields
    [c_k^+ f_b - c_k^- f_a, c_k^+ f_a - c_k^- f_b] term by term."""
    n = min(len(ft_a), len(ft_b))
    return [ivlh(LB(ft_b[k]), UB(ft_a[k])) for k in range(n)]


_ISQ = {}


def isq_table(n):
    """[1/sqrt(2j+1)] as intervals, j = 0 .. n-1.  Monotone cache: a LONGER
    table is returned when one is already built (every consumer either indexes
    j < n or zips against a length-n list, so extra entries are inert)."""
    t = _ISQ.get("t")
    if t is not None and len(t) >= n:
        return t
    t = [IONE / iv.sqrt(ivf(2 * j + 1)) for j in range(n)]
    _ISQ["t"] = t
    return t


# ===========================================================================
# SEC. 2.  ENVELOPES   spec = (C, A, p)  =>  env(k) = C/(k^2+A^2)^p
# ===========================================================================
# `key` is a hashable exact tag used only for caching.

class Spec:
    __slots__ = ("C", "A", "p", "key")

    def __init__(self, C, A, p, key):
        self.C = ivf(C)
        self.A = ivf(A)
        self.p = int(p)
        self.key = key

    def __iter__(self):
        return iter((self.C, self.A, self.p))


# psi_j(t; w) = c_j(w) / (t^2 + w^2)^{PSI_P[j]}  is an upper envelope for
# |d_w^j [[t]](w)|, [[t]](w) = w^2/(t^2+w^2)  (verbatim from
# p4_interchange_uniformity.py, block comment above its PSI_P).  E.g. j = 1:
# |d_w [[t]]| = 2w t^2/(t^2+w^2)^2 <= 2w/(t^2+w^2);  j = 2:
# |d_ww [[t]]| = 2t^2|t^2-3w^2|/(t^2+w^2)^3 <= 6/(t^2+w^2).
PSI_P = {0: 1, 1: 1, 2: 1, 3: 2, 4: 2, 5: 3}


def _psi_c_iv(order, w):
    """c_j(w) in {w^2, 2w, 6, 24w, 240, 2400w} -- the numerator of the
    order-j elementwise envelope (verbatim from piu._psi_c)."""
    w = ivf(w)
    return {0: w ** 2, 1: 2 * w, 2: ivf(6), 3: 24 * w,
            4: ivf(240), 5: 2400 * w}[order]


def _kk(x):
    """Exact hashable tag for a cache key (float() would raise on a non-point
    interval, and repr(float(.)) would silently merge distinct mpf values)."""
    if isinstance(x, IVT):
        return x._mpi_
    return mp.mpf(x)._mpf_


def spec_point_iv(w, order=0):
    return Spec(_psi_c_iv(order, w), w, PSI_P[order], ("pt", _kk(w), order))


def spec_seg_iv(a, b, order=0):
    """Valid for EVERY w in [a,b] simultaneously: c_j is nondecreasing in w and
    (t^2+w^2)^-p is increasing in w (piu.spec_seg)."""
    return Spec(_psi_c_iv(order, b), a, PSI_P[order],
                ("sg", _kk(a), _kk(b), order))


def spec_maj_iv(a):
    """(1, a, 2): for every wbar >= a >= 0 and k >= 1,
         phi_w(k) = k^2/(k^2+w^2)^3 <= 1/(k^2+w^2)^2 <= 1/(k^2+a^2)^2,
         psi_w(k) = 1/(k^2+w^2)^2   <= 1/(k^2+a^2)^2,
         kappa_w(k) = 1/(k^2(k^2+w^2)) <= 1/k^4      (a = 0 case)."""
    return Spec(IONE, a, 2, ("mj", _kk(a)))


def spec_d3_iv(a, coef):
    return Spec(coef, a, 3, ("d3", _kk(a), _kk(coef)))


def env_at(spec, k):
    return spec.C / (ivf(k) ** 2 + spec.A ** 2) ** spec.p


# --- certified arctangent ---------------------------------------------------
# mpmath's `iv.atan2` / `iv.atan` are NOT outward rounded (mpmath 1.3.0):
# `mpf_atan2` evaluates `mpf_atan(mpf_div(y, x, prec + 4), prec + 4)` with
# ROUND-TO-NEAREST and only directs the FINAL rounding, so the interval it
# returns can MISS the exact value by a fraction of an ulp (measured on this
# build: ~9% of random operands, worst 0.24 ulp at prec 103, and 14 misses in
# 20000 calls when fed through `tail_int`).  An unquantified floating-point
# error of exactly the kind this module exists to remove -- and `tail_int` is
# on the load-bearing path of every Stage L cell.  So arctan is done here with
# nothing but +, -, *, /, sqrt and integer powers, the operations for which
# mpmath.iv IS provably outward rounded.  (Adversarial review M-1,
# `notes/2026-09-12_ci_interval_review.md`.)
ATAN_HALVINGS = 80
ATAN_TERMS = 20
_ATAN = {}


def iv_atan_nn(t):
    r"""Enclosure of atan(t) for an interval t >= 0.

    Argument halving:  atan(t) = 2 * atan( t / (1 + sqrt(1 + t^2)) ).  This is
    an exact identity for every t >= 0 (put t = tan(2u)).  One application maps
    [0, inf) into [0, 1), and from there the value drops below 1/16 in at most
    five more steps, so the loop always terminates -- asserted.

    Then the alternating series  atan(u) = sum_j (-1)^j u^(2j+1)/(2j+1), whose
    terms decrease monotonically for 0 <= u < 1; hence the remainder after N
    terms is bounded IN ABSOLUTE VALUE by the first omitted term
    u^(2N+1)/(2N+1).  With u <= 1/16 and N = 20 that is <= 16^-41 < 2^-164,
    i.e. far below one ulp at prec 103, so this costs nothing in width.
    """
    t = t if isinstance(t, IVT) else ivf(t)
    key = t._mpi_
    if key in _ATAN:
        return _ATAN[key]
    assert LB(t) >= 0, "iv_atan_nn: needs t >= 0"
    u, m = t, 0
    while UB(u) > mp.mpf(0.0625) and m < ATAN_HALVINGS:
        u = u / (IONE + iv_sqrt_nn(IONE + u ** 2))
        m += 1
    assert UB(u) <= mp.mpf(0.0625), "iv_atan_nn: halving did not converge"
    u2, s, tp = u ** 2, IZERO, u
    for j in range(ATAN_TERMS):
        c = tp / ivf(2 * j + 1)
        s = s + c if (j % 2 == 0) else s - c
        tp = tp * u2                       # tp = u^(2j+3) at the end of step j
    rem = UB(tp / ivf(2 * ATAN_TERMS + 1))  # >= |first omitted term|, tp >= 0
    out = (s + ivlh(-rem, rem)) * (2 ** m)
    _ATAN[key] = out
    return out


def iv_pi_half():
    """pi/2 = 2*atan(1), from the certified series above -- so that not even
    `iv.pi` (whose outward rounding rests on a comment in mpmath's `pi_fixed`)
    is trusted."""
    return 2 * iv_atan_nn(IONE)


def tail_int(spec, x):
    """UPPER interval enclosure of  int_x^inf C dt/(t^2+A^2)^p,  x >= 0.

    A = 0 branch: C/((2p-1) x^(2p-1)).
    A > 0 branch: written with  at := atan(A/x) = pi/2 - atan(x/A)  (identity
    atan(t) + atan(1/t) = pi/2 for t > 0).  Using atan(A/x) DIRECTLY rather
    than the difference pi/2 - atan(x/A) avoids 2-3 digits of cancellation at
    large x/A -- irrelevant at dps 30 but free, and it keeps every interval
    width at its natural size.  At x = 0 the value is atan(+inf) = pi/2.
    `iv_atan_nn`, NOT `iv.atan2`: see its docstring."""
    C, A, p = spec.C, spec.A, spec.p
    x = ivf(x)
    if UB(A) == 0 and LB(A) == 0:
        assert LB(x) > 0, "tail_int: A = 0 branch needs x > 0"
        return C / ((2 * p - 1) * x ** (2 * p - 1))
    if UB(x) == 0 and LB(x) == 0:
        at = iv_pi_half()
    else:
        at = iv_atan_nn(A / x)
    x2A2 = x ** 2 + A ** 2
    if p == 1:
        return C * at / A
    if p == 2:
        return C * (at / (2 * A ** 3) - x / (2 * A ** 2 * x2A2))
    if p == 3:
        return C * (3 * at / (8 * A ** 5) - x / (4 * A ** 2 * x2A2 ** 2)
                    - 3 * x / (8 * A ** 4 * x2A2))
    raise ValueError(p)


_S1 = {}
S1_K = 3000


def S1_of(spec):
    """UPPER bound on sum_{k>=1} env(k): partial sum to K plus the integral
    tail (env decreasing => sum_{k>K} env(k) <= int_K^inf env)."""
    if spec.key in _S1:
        return _S1[spec.key]
    s = ivsum([env_at(spec, k) for k in range(1, S1_K + 1)]) + tail_int(spec, S1_K)
    _S1[spec.key] = s
    return s


def tail_sum(spec, J):
    """UPPER bound on sum_{k>=J} env(k)  (env decreasing => <= int_{J-1}^inf)."""
    if J <= 1:
        return S1_of(spec)
    return tail_int(spec, J - 1)


def row_env(n, spec):
    """Manifestly decreasing upper envelope for the full row sum
    sum_{m>=0} |O[f]_{nm}|, valid for n >= 1.  Verbatim algebra of piu.row_env
    (five decreasing terms; see its docstring for the derivation)."""
    assert n >= 1
    S1 = S1_of(spec)
    nn = ivf(n)
    r = iv.sqrt((nn + 1) * (2 * nn + 1))
    half = spec.C / ((nn / 2) ** 2 + spec.A ** 2) ** spec.p
    at_n = spec.C / (nn ** 2 + spec.A ** 2) ** spec.p
    s2 = iv.sqrt(ivf(2))
    return (2 * S1 / (2 * nn + 1)
            + 2 * S1 / r
            + half / s2
            + tail_int(spec, nn) / r
            + at_n / s2)


def row_env0(spec):
    """Row n = 0 (piu.row_env0): 2*S1."""
    return 2 * S1_of(spec)


def col_env(m, N, spec):
    """Upper bound on sum_{n>=N} |O[f]_{nm}| for m < N (piu.col_env)."""
    return ((tail_int(spec, N - m - 1) + tail_int(spec, N + m))
            / iv.sqrt(ivf((2 * m + 1) * (2 * N + 1))))


def schur_coupling(Na, spec):
    """||B||, B = (O[f]_{nm})_{n<Na, m>=Na}  (piu.schur_coupling)."""
    cmax = IZERO
    best = mp.mpf(-1)
    for n in range(Na):
        c = col_env(n, Na, spec)
        if UB(c) > best:
            best = UB(c)
            cmax = c
    return iv.sqrt(cmax * row_env(Na, spec))


def schur_tailblock(Na, spec):
    """||D||, D = (O[f]_{nm})_{n,m >= Na}  <= row_env(Na)  (piu.schur_tailblock)."""
    return row_env(Na, spec)


def coupling_restricted(Ns, Nr, spec):
    """||(O[f])_{m >= Nr, n < Ns}||, Ns <= Nr.  Stage M.2 of p4_conj2_interval."""
    best = mp.mpf(-1)
    cmax = IZERO
    for n in range(Ns):
        c = col_env(n, Nr, spec)
        if UB(c) > best:
            best = UB(c)
            cmax = c
    psi_r = spec.C / (ivf(Nr - Ns + 1) ** 2 + spec.A ** 2) ** spec.p
    rmax = 2 * psi_r * iv.sqrt(ivf(2 * Ns - 1)) / iv.sqrt(ivf(2 * Nr + 1))
    return iv.sqrt(rmax * cmax)


_RSUM = {}
RSUM_N0 = 60
RSUM_M = 900


def rsum(spec, n0=RSUM_N0, M=RSUM_M):
    """RIGOROUS  ||O[f]|| <= max_n sum_m |O[f]_{nm}|  (Schur/Holmgren) for every
    kernel with |f(k)| <= env(k).  Rows n < n0 explicitly (m < M exact plus the
    closed-form remainders), rows n >= n0 by the decreasing bulk bound.  Same
    algebra as p4_extremal_concavity.rsum, in iv."""
    key = (spec.key, n0, M)
    if key in _RSUM:
        return _RSUM[key]
    ev = [IZERO] + [env_at(spec, k) for k in range(1, M + n0 + 3)]
    cum = [IZERO] * (n0 + 1)
    acc = IZERO
    for n in range(1, n0 + 1):
        acc = acc + ev[n]
        cum[n] = acc
    invm = isq_table(M)
    S1 = S1_of(spec)
    isqM = IONE / iv.sqrt(ivf(2 * M + 1))
    best = IZERO
    bub = mp.mpf(-1)
    for n in range(n0):
        row = 2 * cum[n] / (2 * n + 1)
        if n < M:
            r2 = ev[n:0:-1] + [IZERO] + ev[1:M - n]
        else:
            r2 = ev[n:n - M:-1]
        acc2 = ivdot(r2[:M], invm) + ivdot(ev[n + 1:n + M + 1], invm)
        acc2 = acc2 + (tail_sum(spec, M - n) + tail_sum(spec, M + n + 1)) * isqM
        row = row + acc2 / iv.sqrt(ivf(2 * n + 1))
        if UB(row) > bub:
            bub = UB(row)
            best = row
    bulk = (2 * S1 / (2 * n0 + 1) + 2 * S1 / iv.sqrt(ivf(2 * n0 + 1))
            + tail_sum(spec, n0 + 1) / iv.sqrt(ivf(2 * n0 + 1)))
    out = ivlh(mp.zero, max(UB(best), UB(bulk)))
    _RSUM[key] = out
    return out


SEG_MROWS = 200
SEG_MCOLS = 1200
_SEG = {}


def rsum_seg(a, b, order, Mrows=SEG_MROWS, Mcols=SEG_MCOLS):
    """RIGOROUS  sup_{w in [a,b]} ||H^{(order)}(w)||_op  -- the iv twin of
    `p4_interval_extension.schur_full_seg_w0` (rows summed from n = 0 with their
    column weights).  Identical algebra; only the column cut differs (Mcols is
    1200 instead of 8000, the remainder being absorbed by the SAME closed-form
    integral tail that the base function already uses beyond Mcols, so the bound
    is still rigorous -- just looser by ~1e-3 relative, which is irrelevant:
    L1 enters only through theta - L1*hw and alpha's Lipschitz slack).

    Row n:  2 cum[n]/(2n+1)                                    (H1 diagonal)
          + (1/sqrt(2n+1)) [ sum_{m<Mcols, m!=n} psi(|n-m|)/sqrt(2m+1)
                           + sum_{m<Mcols}      psi(n+m+1)/sqrt(2m+1)
                           + (IntTail(Mcols-n-1)+IntTail(Mcols+n))/sqrt(2Mcols+1) ]
    Rows n >= Mrows: row_env(Mrows) (decreasing), which is a valid bound for
    every such row because row_env is decreasing in n."""
    key = (repr(float(a)), repr(float(b)), order, Mrows, Mcols)
    if key in _SEG:
        return _SEG[key]
    spec = spec_seg_iv(a, b, order)
    ev = [IZERO] + [env_at(spec, k) for k in range(1, Mrows + Mcols + 3)]
    cum = [IZERO] * (Mrows + 1)
    acc = IZERO
    for n in range(1, Mrows + 1):
        acc = acc + ev[n]
        cum[n] = acc
    invm = isq_table(Mcols)
    isqM = IONE / iv.sqrt(ivf(2 * Mcols + 1))
    bub = mp.mpf(-1)
    for n in range(Mrows):
        row = 2 * cum[n] / (2 * n + 1)
        r2 = ev[n:0:-1] + [IZERO] + ev[1:Mcols - n]
        acc2 = ivdot(r2[:Mcols], invm) + ivdot(ev[n + 1:n + Mcols + 1], invm)
        acc2 = acc2 + (tail_sum(spec, Mcols - n) + tail_sum(spec, Mcols + n + 1)) * isqM
        row = row + acc2 / iv.sqrt(ivf(2 * n + 1))
        if UB(row) > bub:
            bub = UB(row)
    out = ivlh(mp.zero, max(bub, UB(row_env(Mrows, spec))))
    _SEG.clear()
    _SEG[key] = out
    return out


# ===========================================================================
# SEC. 3.  TRIAL VECTOR, STRUCTURE COEFFICIENTS, MATVEC   (all in iv)
# ===========================================================================

class TrialIV:
    """Finite-support unit trial vector vhat together with u_n = vhat_n/sqrt(2n+1),
    the structure coefficients c_k = T_k + U_k - P_k of
    p4_extremal_concavity.py Sec.2.1, and Cu = sum|u_n|.

    vnp is float64 -- its entries are EXACT binary rationals and are treated as
    such.  Normalisation is carried out in iv, so ||vhat|| = 1 only up to the
    interval width; every use of vhat below is therefore written so that it
    remains valid for any vector in the enclosure:
      * `qform` divides by nv2 = <v,v> explicitly (Rayleigh quotient);
      * `_ehat` / the residual use ||.||, and the certified norm is taken as an
        upper endpoint;
      * the (P2) chain only needs SOME unit vector, and vhat/||vhat|| is one.
    """

    def __init__(self, vnp):
        Nv = len(vnp)
        v0 = [ivf(float(t)) for t in vnp]
        nv2 = ivsumsq(v0)
        nrm = iv.sqrt(nv2)
        v = [t / nrm for t in v0]
        self.Nv = Nv
        self.v = v
        isq = isq_table(Nv)
        u = [v[n] * isq[n] for n in range(Nv)]
        self.u = u
        self.Cu = ivsum([abs(t) for t in u])
        # P_k = 2 sum_{n>=k} u_n^2
        P = [IZERO] * (2 * Nv + 2)
        acc = IZERO
        for n in range(Nv - 1, 0, -1):
            acc = acc + u[n] ** 2
            P[n] = 2 * acc
        P[0] = 2 * (acc + u[0] ** 2)
        # T_k = 2 sum_n u_n u_{n+k}
        T = [IZERO] * (2 * Nv + 2)
        for j in range(1, Nv):
            T[j] = 2 * ivdot(u[:Nv - j], u[j:])
        # U_s = sum_{n+m+1=s} u_n u_m
        U = [IZERO] * (2 * Nv + 2)
        for s in range(1, 2 * Nv):
            lo = max(0, s - Nv)
            hi = min(s - 1, Nv - 1)
            if hi < lo:
                continue
            U[s] = ivdot(u[lo:hi + 1], u[s - 1 - hi:s - lo][::-1])
        c = [T[k] + U[k] - P[k] for k in range(2 * Nv)]
        c[0] = IZERO
        self.c = c
        self.all_nonneg = all(LB(t) >= 0 for t in c)
        self.min_c = min(LB(t) for t in c[1:])

    def qform(self, ft):
        """<vhat, O[f] vhat> as an interval (vhat already unit-normalised in iv)."""
        n = min(len(ft), len(self.c))
        return ivdot(self.c[:n], ft[:n])

    def qform_cell(self, ftcell):
        """The cell-uniform enclosure of qhat(w) for w in [a,b]: feed the
        interval table [f_b(k), f_a(k)] (see `cell_table`).  Equivalent to the
        sign split qhat_lo / qhat_hi of p4_extremal_concavity.py Sec.2.1."""
        n = min(len(ftcell), len(self.c))
        return ivdot(self.c[:n], ftcell[:n])

    def matvec(self, ft, M):
        """y_n = (O[f] vhat)_n, n < M.  Needs ft up to index M + Nv."""
        Nv, u, v = self.Nv, self.u, self.v
        isq = isq_table(max(M, Nv))
        y = []
        run = IZERO
        for n in range(M):
            if n >= 1:
                run = run + 2 * ft[n]
            if n < Nv:
                row2 = ft[n:0:-1] + [IZERO] + ft[1:Nv - n]
            else:
                row2 = ft[n:n - Nv:-1]
            s = ivdot(row2, u) + ivdot(ft[n + 1:n + Nv + 1], u)
            yn = s * isq[n]
            if n < Nv:
                yn = yn - (run * isq[n] ** 2) * v[n]
            y.append(yn)
        return y

    def tail_sq(self, spec, M):
        """Upper bound on sum_{n>=M} y_n^2 for |f| <= env(spec), env non-increasing
        (p4_extremal_concavity.Trial.tail_sq, same derivation):
          sum_{n>=M} y_n^2 <= (4 Cu^2/(2M+1)) sum_{j>=M-Nv+1} F(j)^2,
          sum_{j>=J} F(j)^2 <= C^2/((4p-1)(J-1)^{4p-1})   (integral test, J>=2)."""
        C, A, p = spec.C, spec.A, spec.p
        J = M - self.Nv + 1
        chk("matvec tail: row cut beyond the support", J >= 2, f"M={M} Nv={self.Nv}")
        s2 = C ** 2 / ((4 * p - 1) * ivf(J - 1) ** (4 * p - 1))
        return 4 * self.Cu ** 2 / (2 * ivf(M) + 1) * s2


def ehat_iv(tr, ft, M, spec):
    """UB of ||(I - vhat vhat*) O[f] vhat||: exact iv head + closed-form row tail."""
    y = tr.matvec(ft, M)
    q = ivdot(tr.v, y[:tr.Nv])
    s2 = IZERO
    for n in range(M):
        dd = y[n] - (q * tr.v[n] if n < tr.Nv else IZERO)
        s2 = s2 + dd ** 2
    s2 = s2 + tr.tail_sq(spec, M)
    return UB(iv_sqrt_nn(s2))


# ===========================================================================
# SEC. 4.  INTERVAL CHOLESKY AND THE lam_2 CERTIFICATE
# ===========================================================================

def chol_pd(S, n):
    """Interval Cholesky.  S is an n x n list-of-lists of intervals containing a
    symmetric real matrix.  Returns (ok, min_pivot_lb).

    By the theorem in the module docstring, ok=True (all pivots have strictly
    positive lower endpoint) PROVES that every symmetric real matrix in [S] is
    positive definite.  Row-wise formulation:
        L[i][j] = (S[i][j] - sum_{k<j} L[i][k] L[j][k]) / L[j][j],   j < i
        L[i][i] = sqrt(S[i][i] - sum_{k<i} L[i][k]^2)."""
    L = []
    worst = None
    for i in range(n):
        Si = S[i]
        Li = []
        for j in range(i):
            Lj = L[j]
            t = Si[j] - ivdot(Li, Lj)
            Li.append(t / Lj[j])
        d = Si[i] - ivsumsq(Li)
        dl = LB(d)
        if worst is None or dl < worst:
            worst = dl
        if dl <= 0:
            return False, dl
        Li.append(iv.sqrt(ivlh(dl, UB(d))))
        L.append(Li)
    return True, worst


def head_iv(ft, Na):
    """The exact head block A = (O[f]_{nm})_{n,m < Na} as an interval matrix,
    built DIRECTLY from the iv kernel table (this is the step the base
    `alpha_certificate_H` skipped: it built the head in float64 and then copied
    the float64 numbers into mpf, so its Cholesky certified the float64 matrix,
    not the operator's head).  ft must reach index 2Na."""
    isq = isq_table(Na)
    diag = []
    run = IZERO
    for n in range(Na):
        if n >= 1:
            run = run + 2 * ft[n]
        diag.append(run * isq[n] ** 2)
    A = []
    for n in range(Na):
        isn = isq[n]
        row = []
        ftn = ft
        for m in range(Na):
            w = isn * isq[m]
            val = ftn[n + m + 1] * w
            if n != m:
                val = val + ftn[abs(n - m)] * w
            else:
                val = val - diag[n]
            row.append(val)
        A.append(row)
    return A


CHOL_C0 = mp.mpf("1e-12")
CHOL_GROW = mp.mpf(4)
CHOL_TRIES = 64


def lam2_head_iv(ft, Na, uf, c0=None):
    """CERTIFIED  c  with  lam_2(A) <= c  for the exact head block A = O[f]|_{<Na}.

    For ANY unit vector u and P = I - u u^T,
        lam_2(A) = min_{dim S = 1} max_{x perp S} <x,Ax>
                <= max_{x perp u} <x,Ax> = lam_max(P A P|_{u^perp})
                <= lam_max(P A P),
    so a positive-definiteness certificate for cI - PAP gives lam_2(A) <= c.
    u is supplied in float64 (the top eigenvector of the float64 head) purely as
    a deflation direction -- ANY unit vector works, so no bound depends on its
    accuracy.  It is normalised inside iv.

    c0 (optional, float64 lam_2 of the float64 head) is only a STARTING GUESS:
    c is raised by factors of CHOL_GROW until the interval Cholesky succeeds, and
    a larger c is a weaker (still valid) conclusion.  The seed keeps the
    overshoot below CHOL_GROW; without it the escalation can overshoot the true
    lam_2 by orders of magnitude, which costs margin in alpha (lam_2 of the head
    is genuinely positive, ~2.4e-4 at wbar = 1, ~3.6e-5 at wbar = 0.3125).
    Returns (c, pivot_lb, tries)."""
    A = head_iv(ft, Na)
    u0 = [ivf(float(t)) for t in uf]
    un = iv.sqrt(ivsumsq(u0))
    u = [t / un for t in u0]
    Au = [ivdot(A[i], u) for i in range(Na)]
    uAu = ivdot(u, Au)
    # PAP = A - u (Au)^T - (Au) u^T + (u^T A u) u u^T   (A symmetric)
    negPAP = []
    for i in range(Na):
        ui, Aui = u[i], Au[i]
        Ai = A[i]
        row = [-(Ai[j] - ui * Au[j] - Aui * u[j] + (ui * u[j]) * uAu)
               for j in range(Na)]
        negPAP.append(row)
    c = CHOL_C0
    if c0 is not None and float(c0) > 0:
        c = max(c, add_hi(mul_hi(mp.mpf(float(c0)), mp.mpf("1.000001")), CHOL_C0))
    for t in range(CHOL_TRIES):
        civ = ivf(c)
        S = []
        for i in range(Na):
            row = list(negPAP[i])
            row[i] = row[i] + civ
            S.append(row)
        ok, piv = chol_pd(S, Na)
        if ok:
            return c, piv, t + 1
        c = c * CHOL_GROW
    raise AssertionError(f"lam2_head_iv: interval Cholesky failed for Na={Na} "
                         f"up to c={c}")


def alpha_cell_iv(A, B, cA, cB, L1_ub, Na):
    """Cell-uniform  alpha >= sup_{w in [A,B]} lam_2(H(w)):
       lam_2(H(w)) <= max(lam_2(head), ||D||_seg) + ||B||_seg
                   <= max(max(cA,cB) + L1 (B-A)/2, ||D||_seg) + ||B||_seg.
    The Lipschitz term covers w in [A,B] from the nearer endpoint via Weyl and
    ||d_w H|| <= L1.  (The base `alpha_cell` uses c + ||D|| + ||B||; max(.,.) is
    the tighter and equally valid form of the same block bound, since the
    spectrum of blockdiag(head, D) is spec(head) u spec(D) and hence its second
    largest value is <= max(lam_2(head), lam_1(D)).)  Returns UB as mpf."""
    sp = spec_seg_iv(A, B, 0)
    nB = schur_coupling(Na, sp)
    nD = schur_tailblock(Na, sp)
    chead = ivf(max(cA, cB)) + ivf(L1_ub) * (ivf(B) - ivf(A)) / 2
    m = ivlh(mp.zero, max(UB(chead), UB(nD)))
    return UB(m + nB)


def mu2_cell_iv(a, b):
    """RIGOROUS lower bound on lam_1(H^(2)(w)) -- hence on lam_1(H(w)) by
    min-max on the 2-dimensional subspace -- valid for EVERY w in [a,b].
    piu.mu2_cell_lb in iv:  v = top eigenvector of the 2x2 head H^(2)(a) in
    closed form; Q(w) = k1[[1]](w)+k2[[2]](w)+k3[[3]](w) <= lam_1(H^(2)(w));
    every [[k]](w) = w^2/(k^2+w^2) is increasing in w, so k1,k2,k3 >= 0 implies
    Q(w) >= Q(a) = lam_1(H^(2)(a)).  The signs of k1,k2,k3 are ASSERTED.
    Returns (LB(lam), ok)."""
    w = ivf(a)
    w2 = w ** 2
    b1 = w2 / (ivf(1) + w2)
    b2 = w2 / (ivf(4) + w2)
    b3 = w2 / (ivf(9) + w2)
    s3 = iv.sqrt(ivf(3))
    h00 = b1
    h01 = (b1 + b2) / s3
    h11 = (b3 - 2 * b1) / 3
    tr = h00 + h11
    det = h00 * h11 - h01 ** 2
    disc = tr ** 2 - 4 * det
    lam = (tr + iv_sqrt_nn(disc)) / 2
    v0 = h01
    v1 = lam - h00
    if LB(v0) <= 0 or LB(v1) <= 0:
        return mp.zero, False
    nn = v0 ** 2 + v1 ** 2
    k1 = (v0 ** 2 + 2 * v0 * v1 / s3 - 2 * v1 ** 2 / 3) / nn
    k2 = (2 * v0 * v1 / s3) / nn
    k3 = (v1 ** 2 / 3) / nn
    ok = LB(k1) >= 0 and LB(k2) >= 0 and LB(k3) >= 0
    return LB(lam), ok


# ---------------------------------------------------------------------------
# the G(2) analogue used by Stage S
# ---------------------------------------------------------------------------

def tailblock_G2(Na, spec):
    """||D|| for the tail block of O[f], f <= env(spec) (psw.tailblock_bound):
       L1 <= 2 S1/(2Na+1); L2 <= 2 S1/(2Na+1); L3 <= tail_sum(2Na+1)/(2Na+1)."""
    return (4 * S1_of(spec) + tail_sum(spec, 2 * Na + 1)) / (2 * ivf(Na) + 1)


def coupling_G2(Na, spec):
    """||B|| for the head x tail block (psw.coupling_bound), A = 0 safe."""
    sqNa = iv.sqrt(2 * ivf(Na) + 1)
    rmax = IZERO
    rb = mp.mpf(-1)
    for n in range(Na):
        r = (tail_sum(spec, Na - n) + tail_sum(spec, Na + n + 1)) \
            / (iv.sqrt(ivf(2 * n + 1)) * sqNa)
        if UB(r) > rb:
            rb = UB(r)
            rmax = r
    cmax = IZERO
    for n in range(Na):
        cmax = cmax + (env_at(spec, Na - n) + env_at(spec, Na + n + 1)) \
            / iv.sqrt(ivf(2 * n + 1))
    cmax = cmax / sqNa
    return iv.sqrt(rmax * cmax)


def alpha_G2_iv(Na, ftab):
    """alpha >= lam_2(G(2)), fully in iv.  Head deflation as in lam2_head_iv,
    tail/coupling by the k^-2 envelope spec (1, 0, 1).
    Returns (alpha_ub, c, nD_ub, nB_ub, pivot)."""
    A_np = op_np(Na, np_tab(2 * Na + 2, lambda k: k ** -2.0))
    ev, vec = np.linalg.eigh(A_np)
    c, piv, _ = lam2_head_iv(ftab, Na, vec[:, -1], c0=float(ev[-2]))
    spec = Spec(IONE, 0, 1, ("g2env",))
    nD = tailblock_G2(Na, spec)
    nB = coupling_G2(Na, spec)
    return add_hi(max(c, UB(nD)), UB(nB)), c, UB(nD), UB(nB), piv


# ===========================================================================
# SEC. 5.  float64 helpers (TRIAL VECTORS AND EVIDENCE ONLY -- never a bound)
# ===========================================================================

def np_tab(kmax, fun):
    k = np.arange(0, kmax + 1, dtype=np.float64)
    out = np.zeros_like(k)
    out[1:] = fun(k[1:])
    return out


def op_np(N, ftab):
    n = np.arange(N)
    sq = np.sqrt(2.0 * n + 1.0)
    diag = np.zeros(N)
    run = 0.0
    for k in range(N):
        if k >= 1:
            run += 2.0 * ftab[k]
        diag[k] = run / (2 * k + 1)
    Dm = np.abs(n[:, None] - n[None, :])
    L2 = np.where(Dm > 0, ftab[np.where(Dm > 0, Dm, 1)] / (sq[:, None] * sq[None, :]), 0.0)
    Sm = n[:, None] + n[None, :] + 1
    L3 = ftab[Sm] / (sq[:, None] * sq[None, :])
    return -np.diag(diag) + L2 + L3


def trial_vec_np(w, Ns):
    """float64 top eigenvector of the head of H(w) -- a TRIAL VECTOR generator."""
    tb = np_tab(2 * Ns + 2, lambda k: w * w / (k * k + w * w))
    vnp = np.linalg.eigh(op_np(Ns, tb))[1][:, -1]
    return vnp * np.sign(vnp[0])


# ===========================================================================
# SEC. 6.  THE PER-CELL CERTIFICATES
# ===========================================================================

def support_size(w):
    return max(140, int(10 * float(w)) // 2 * 2)


NA_RULE = ["band"]


def na_band(w):
    w = float(w)
    if NA_RULE[0] == "ext":
        if w > 20.0:
            return min(1000, -(-int(14 * w) // 20) * 20)
    if w <= 5.0:
        return 70
    if w <= 12.0:
        return 140
    return 280


_TRIAL = {}


def trial_at(p, Ns):
    key = (repr(float(p)), Ns)
    if key in _TRIAL:
        return _TRIAL[key]
    tr = TrialIV(trial_vec_np(float(p), Ns))
    _TRIAL.clear()
    _TRIAL[key] = tr
    return tr


TARGET = 0.75


def cell_L_iv(a, b, alpha_ub, L1_ub, mu2_lb):
    """Certified  sup_{x in [a^2,b^2]} H''(x)  on a cell [a,b], 0 < a <= b.
    Every quantity is computed in iv and collapsed to the single endpoint the
    (C**) assembly needs; the monotonicity justifying each choice is in the
    module docstring.  Returns (hpp_ub as float or None, info)."""
    a, b = float(a), float(b)
    # The endpoints are integer multiples of 2^-16 and |a|,|b| < 2^8, so a+b and
    # b-a need <= 25 significant bits: p and hw are EXACT in float64 and hence
    # exact point intervals.  This matters: it makes "distance from p to any
    # point of [a,b] <= hw" an exact statement rather than one needing slack.
    p = 0.5 * (a + b)
    hw_f = 0.5 * (b - a)
    assert p - a == hw_f and b - p == hw_f, "cell midpoint is not exact"
    aiv, biv = ivf(a), ivf(b)
    piv = ivf(p)
    hw = ivf(hw_f)
    Ns = support_size(p)
    Nr = 4 * Ns
    tr = trial_at(p, Ns)
    KT = Nr + Ns + 8
    fH = tab_H_iv(piv, KT)

    # --- Rayleigh, residual (incl. the MANDATORY tail coupling), angle at p ---
    theta_i = tr.qform(fH)                 # <= lam_1(H(p)); use LB
    theta_lo = LB(theta_i)
    y = tr.matvec(fH, Nr)
    r2 = IZERO
    for n in range(Nr):
        dd = y[n] - (ivf(theta_lo) * tr.v[n] if n < Ns else IZERO)
        r2 = r2 + dd ** 2
    rtail = coupling_restricted(Ns, Nr, spec_point_iv(piv, 0))
    rnorm_ub = UB(iv_sqrt_nn(r2 + rtail ** 2))
    gamma_lo = sub_lo(theta_lo, alpha_ub)
    lam1_lo = max(LB(ivf(theta_lo) - ivf(L1_ub) * hw), mu2_lb)
    g_lo = sub_lo(lam1_lo, alpha_ub)
    if gamma_lo <= 0 or g_lo <= 0:
        return None, dict(reason="no certified gap", gamma=float(gamma_lo),
                          g=float(g_lo))
    s_p = div_hi(rnorm_ub, gamma_lo)

    # --- phi / psi pieces ----------------------------------------------------
    MK = Ns + 300
    sp2 = spec_maj_iv(a)
    R3_ub = UB(rsum(spec_d3_iv(a, 1)))
    N_ub = UB(rsum(sp2))
    fphi_cell = cell_table(tab_phi_iv(biv, KT), tab_phi_iv(aiv, KT))
    fpsi_cell = cell_table(tab_psi_iv(biv, KT), tab_psi_iv(aiv, KT))
    qphi = tr.qform_cell(fphi_cell)
    qpsi = tr.qform_cell(fpsi_cell)
    qphi_lo = LB(qphi)
    aphi_ub = absmax_ub(qphi)
    apsi_ub = absmax_ub(qpsi)
    eph = ehat_iv(tr, tab_phi_iv(piv, KT), MK, sp2)
    eps_ = ehat_iv(tr, tab_psi_iv(piv, KT), MK, sp2)
    # |d_w ehat_f| <= ||O[d_w f]||;  |d_w psi| <= 4b/(k^2+a^2)^3,
    # |d_w phi| <= 6b k^2/(k^2+a^2)^4 <= 6b/(k^2+a^2)^3.
    eph_hi = UB(ivf(eph) + hw * 6 * biv * ivf(R3_ub))
    eps_hi = UB(ivf(eps_) + hw * 4 * biv * ivf(R3_ub))

    # --- cell-uniform angle: the monotone bootstrap (V1)-(V3) ----------------
    # s |-> s_p + hw*2b^3*(eps_hi + s*(N+|qpsi|))/g_lo is an increasing affine
    # map sending VALID upper bounds on sin angle to VALID upper bounds; every
    # iterate is independently legitimate.  Start from E <= N.
    b3 = biv ** 3
    s = min(mp.one, UB(ivf(s_p) + hw * 2 * b3 * ivf(N_ub) / ivf(g_lo)))
    for _ in range(12):
        E_hi = UB(ivf(eps_hi) + ivf(s) * (ivf(N_ub) + ivf(apsi_ub)))
        s_new = UB(ivf(s_p) + hw * 2 * b3 * ivf(E_hi) / ivf(g_lo))
        if s_new >= s:
            break
        s = s_new
    if s >= 1:
        return None, dict(reason="Davis-Kahan angle >= 1", s=float(s))
    E_hi = UB(ivf(eps_hi) + ivf(s) * (ivf(N_ub) + ivf(apsi_ub)))
    F_lo = LB(ivf(qphi_lo) - ivf(s) ** 2 * (ivf(aphi_ub) + ivf(N_ub))
              - 2 * ivf(s) * ivf(eph_hi))
    pos_ub = UB(biv ** 4 * ivf(E_hi) ** 2 / ivf(g_lo))
    hpp_ub = UB(-2 * ivf(F_lo) + 2 * ivf(pos_ub))
    ratio = float(div_hi(pos_ub, F_lo)) if F_lo > 0 else 9e9
    info = dict(a=float(a), b=float(b), p=float(p), Ns=Ns, Nr=Nr,
                F_lo=float(F_lo), pos=float(pos_ub), s=float(s), s_p=float(s_p),
                g_lo=float(g_lo), E_hi=float(E_hi), theta=float(theta_lo),
                alpha=float(alpha_ub), rnorm=float(rnorm_ub), rtail=float(UB(rtail)),
                ehat_phi=float(eph), ehat_psi=float(eps_), N=float(N_ub),
                hpp_ub=float(hpp_ub), ratio=ratio, regime="L",
                cnn=bool(tr.all_nonneg), min_c=float(tr.min_c),
                lip_psi=float(eps_hi - eps_),
                ang_psi=float(mul_hi(s, add_hi(N_ub, apsi_ub))),
                gamma=float(gamma_lo))
    return hpp_ub, info


# ---------------------------------------------------------------------------
# Stage S: the first cell (0, w0], rescaled operator h~ = G(2) - w^2 K(w)
# ---------------------------------------------------------------------------

NV_S = 200
M_RES_S = 2400
NA_G2 = 200
SDATA = {}


def stage_S_setup():
    if "tr" in SDATA:
        return SDATA
    KG = 2 * NV_S + M_RES_S + 8
    fG2 = tab_gamma_iv(2, KG)
    vnp = np.linalg.eigh(op_np(NV_S, np_tab(2 * NV_S + 2, lambda k: k ** -2.0)))[1][:, -1]
    vnp = vnp * np.sign(vnp[0])
    chk("Stage S trial vector is entrywise positive (Perron structure)",
        float(np.min(vnp)) > 0, f"min={float(np.min(vnp)):.3e}")
    tr = TrialIV(vnp)
    lamhat_lo = LB(tr.qform(fG2))                  # <= frak_g(2)  (Rayleigh)
    y = tr.matvec(fG2, M_RES_S)
    r2 = IZERO
    for n in range(M_RES_S):
        dd = y[n] - (ivf(lamhat_lo) * tr.v[n] if n < NV_S else IZERO)
        r2 = r2 + dd ** 2
    r2 = r2 + tr.tail_sq(Spec(IONE, 0, 1, ("g2env",)), M_RES_S)
    rG2_ub = UB(iv_sqrt_nn(r2))
    fG2a = tab_gamma_iv(2, 2 * NA_G2 + 4)
    alpha_G2, cG2, nDG2, nBG2, pivG2 = alpha_G2_iv(NA_G2, fG2a)
    chk("alpha_G2 < lamhat (Davis-Kahan hypothesis on G(2))", alpha_G2 < lamhat_lo,
        f"{mp.nstr(alpha_G2,8)} vs {mp.nstr(lamhat_lo,8)}")
    chk("alpha_G2 beats the P4-C published 2.0e-2 bound on lam_2(G(2))",
        alpha_G2 <= mp.mpf("2.1e-2"), mp.nstr(alpha_G2, 8))
    chk("interval Cholesky pivot for lam_2(G(2)) is strictly positive", pivG2 > 0,
        mp.nstr(pivG2, 6))
    SDATA.update(tr=tr, lamhat=lamhat_lo, rG2=rG2_ub, alpha_G2=alpha_G2, fG2=fG2)
    print(f"  G(2) trial: lamhat = {mp.nstr(lamhat_lo,16)} <= frak_g(2),  "
          f"||r_G2|| = {mp.nstr(rG2_ub,6)},  alpha_G2 = {mp.nstr(alpha_G2,8)} "
          f"(c = {mp.nstr(cG2,4)}, ||D|| = {mp.nstr(nDG2,4)}, ||B|| = {mp.nstr(nBG2,4)})")
    return SDATA


def cell_S_iv(a, b):
    """Certified sup_{x in [a^2,b^2]} H''(x) for 0 <= a <= b small, via the
    RESCALED operator h~(w) = H(w)/w^2 = G(2) - w^2 K(w)
    (p4_extremal_concavity.py Sec.2.4)."""
    D = stage_S_setup()
    tr, lamhat_lo, rG2_ub, alG2 = D["tr"], D["lamhat"], D["rG2"], D["alpha_G2"]
    aiv, biv = ivf(a), ivf(b)
    a, b = float(a), float(b)
    p = 0.5 * (a + b)
    hw_f = 0.5 * (b - a)
    assert p - a == hw_f and b - p == hw_f, "cell midpoint is not exact"
    piv = ivf(p)
    hw = ivf(hw_f)
    KT = 2 * tr.Nv + 600
    MK = tr.Nv + 300
    sp0 = spec_maj_iv(0.0)                      # 1/k^4 dominates phi, psi, kappa
    sp3 = spec_d3_iv(0.0, 1)

    qk = tr.qform_cell(cell_table(tab_kappa_iv(biv, KT), tab_kappa_iv(aiv, KT)))
    NK_ub = UB(rsum(sp0))
    theta_lo = LB(ivf(lamhat_lo) - biv ** 2 * ivf(ub_pos(qk)))
    alpha_t = UB(ivf(alG2) + biv ** 2 * ivf(NK_ub))
    g_lo = sub_lo(theta_lo, alpha_t)
    if g_lo <= 0:
        return None, dict(reason="no certified gap in the rescaled operator",
                          g=float(g_lo))
    rt_ub = UB(ivf(rG2_ub) + biv ** 2 * (ivf(NK_ub) + ivf(absmax_ub(qk))))
    s = div_hi(rt_ub, g_lo)
    if s >= 1:
        return None, dict(reason="Davis-Kahan angle >= 1", s=float(s))

    N_ub = UB(rsum(sp0))
    R3_ub = UB(rsum(sp3))
    qphi = tr.qform_cell(cell_table(tab_phi_iv(biv, KT), tab_phi_iv(aiv, KT)))
    qpsi = tr.qform_cell(cell_table(tab_psi_iv(biv, KT), tab_psi_iv(aiv, KT)))
    qphi_lo = LB(qphi)
    aphi_ub = absmax_ub(qphi)
    apsi_ub = absmax_ub(qpsi)
    eph = ehat_iv(tr, tab_phi_iv(piv, KT), MK, sp0)
    eps_ = ehat_iv(tr, tab_psi_iv(piv, KT), MK, sp0)
    eph_hi = UB(ivf(eph) + hw * 6 * biv * ivf(R3_ub))
    eps_hi = UB(ivf(eps_) + hw * 4 * biv * ivf(R3_ub))

    F_lo = LB(ivf(qphi_lo) - ivf(s) ** 2 * (ivf(aphi_ub) + ivf(N_ub))
              - 2 * ivf(s) * ivf(eph_hi))
    E_hi = UB(ivf(eps_hi) + ivf(s) * (ivf(N_ub) + ivf(apsi_ub)))
    # g(w) = w^2 g~(w) >= w^2 g_lo  =>  x^2 E^2/g = w^2 E^2/g~ <= b^2 E_hi^2/g_lo
    pos_ub = UB(biv ** 2 * ivf(E_hi) ** 2 / ivf(g_lo))
    hpp_ub = UB(-2 * ivf(F_lo) + 2 * ivf(pos_ub))
    info = dict(a=float(a), b=float(b), F_lo=float(F_lo), pos=float(pos_ub),
                s=float(s), g_lo=float(g_lo), E_hi=float(E_hi), theta=float(theta_lo),
                alpha=float(alpha_t), qk=float(UB(qk)), ehat_phi=float(eph),
                ehat_psi=float(eps_), N=float(N_ub), hpp_ub=float(hpp_ub),
                ratio=float(div_hi(pos_ub, F_lo)) if F_lo > 0 else 9e9, regime="S")
    return hpp_ub, info


# ===========================================================================
# SEC. 7.  STAGE A' -- the interval gate
# ===========================================================================
# Everything here is either (i) an EXACT interval identity check (the enclosure
# must CONTAIN 0 -- strictly stronger than the base script's "float64 difference
# < 1e-12"), or (ii) an explicitly labelled float64 EVIDENCE cross-check whose
# only job is to catch typos.  No bound produced by Secs 1-6 depends on anything
# in this section.

def _contains0(x):
    return LB(x) <= 0 <= UB(x)


def run_gate(full=True):
    print("=== STAGE A': interval gate (primitives, identities, cross-checks) ===")
    t0 = time.time()

    # ---- (a) interval primitives -----------------------------------------
    chk("PAD (the base script's ROUND_PAD) is identically zero", PAD == 0)
    x = ivf(0.3)
    chk("ivf(float) is a POINT interval (a float64 IS an exact binary rational)",
        LB(x) == UB(x) and LB(x) == mp.mpf(0.3))
    y = ivlh(-1, 2)
    chk("x**2 of a sign-straddling interval is [0,4], not [-2,4]",
        LB(y ** 2) == 0 and UB(y ** 2) == 4)
    z = IONE / ivf(3)
    with mp.extraprec(20):                # 3*endpoint is then exact
        strict = 3 * LB(z) < 1 < 3 * UB(z)
    chk("1/3 is enclosed STRICTLY on both sides (outward rounding really "
        "happens; the nearest-rounded mpf would sit on one endpoint)",
        LB(z) < UB(z) and strict)
    chk("LB/UB round-trip is exact", LB(ivlh(LB(z), UB(z))) == LB(z))
    chk("absmax_ub of a straddling interval", absmax_ub(ivlh(-5, 2)) == 5)
    chk("ub_pos clamps a negative interval to 0", ub_pos(ivlh(-5, -2)) == 0)
    # The four one-sided mpf combinators.  These tests are DISCRIMINATING: the
    # operands are chosen so that the exact result is NOT representable at
    # mp.prec and so that round-to-NEAREST goes the WRONG way.  A naive
    # implementation (`x - y` etc. in plain mpf) fails every one of them --
    # verified by mutation testing.  Reference values are computed in
    # extraprec, where they ARE exact (mpf comparison is exact across precs).
    _one, _tiny = mp.mpf(1), mp.ldexp(1, -200)
    _e100 = mp.mpf(1) + mp.ldexp(1, -100)
    with mp.extraprec(400):
        _ex_sub, _ex_add = _one - _tiny, _one + _tiny
        _ex_mul = _e100 * _e100
    chk("sub_lo rounds OUTWARD: 1 - 2^-200 lands strictly below 1 "
        "(plain mpf returns exactly 1, i.e. ABOVE the exact difference)",
        sub_lo(_one, _tiny) <= _ex_sub and (_one - _tiny) > _ex_sub)
    chk("add_hi rounds OUTWARD: 1 + 2^-200 lands strictly above 1 "
        "(plain mpf returns exactly 1, i.e. BELOW the exact sum)",
        add_hi(_one, _tiny) >= _ex_add and (_one + _tiny) < _ex_add)
    chk("mul_hi rounds OUTWARD: (1+2^-100)^2 = 1+2^-99+2^-200 needs 201 bits "
        "(plain mpf rounds it DOWN to 1+2^-99)",
        mul_hi(_e100, _e100) >= _ex_mul and (_e100 * _e100) < _ex_mul)
    _third, _dh = _one / 3, div_hi(_one, mp.mpf(3))   # nearest vs outward
    with mp.extraprec(60):              # 3*(a <=103-bit mpf) is then EXACT
        _ok_dh, _ok_naive = 3 * _dh > 1, 3 * _third < 1
    chk("div_hi rounds OUTWARD: 3*div_hi(1,3) > 1 "
        "(1/3 rounds DOWN at prec 103, so 3*(1/3) < 1 in plain mpf)",
        _ok_dh and _ok_naive)
    chk("iv_sqrt_nn clamps a slightly negative lower endpoint to 0",
        LB(iv_sqrt_nn(ivlh(mp.mpf("-1e-40"), 4))) == 0)

    # ---- (b) the interval Cholesky itself ---------------------------------
    def _mat(rows):
        return [[ivf(t) for t in r] for r in rows]
    ok, piv = chol_pd(_mat([[2, 1, 0], [1, 2, 1], [0, 1, 2]]), 3)
    chk("chol_pd accepts a tridiagonal SPD matrix", ok and piv > 0)
    ok, _ = chol_pd(_mat([[1, 1], [1, 1]]), 2)
    chk("chol_pd REJECTS a singular (psd, not pd) matrix", not ok)
    ok, _ = chol_pd(_mat([[1, 2], [2, 1]]), 2)
    chk("chol_pd REJECTS an indefinite matrix", not ok)
    ok, _ = chol_pd([[ivf(1), ivf(0)], [ivf(0), ivlh(-1e-30, 1e-30)]], 2)
    chk("chol_pd REJECTS when a pivot interval straddles 0 (inf[d] > 0 test)",
        not ok)
    H6 = [[IONE / ivf(i + j + 1) for j in range(6)] for i in range(6)]
    ok, piv = chol_pd(H6, 6)
    chk("chol_pd certifies the 6x6 Hilbert matrix (cond ~1.5e7)", ok and piv > 0,
        f"min pivot {mp.nstr(piv, 6)}")
    ok, piv = chol_pd([[ivf(1e-12), ivf(0)], [ivf(0), ivf(1)]], 2)
    chk("chol_pd certifies a matrix whose smallest eigenvalue is 1e-12 "
        "(the regime of cI - PAP)", ok and piv > 0)

    # ---- (c) IDENTITIES (I1), (I2), (I) -- EXACT interval containments -----
    # O[.] is LINEAR entrywise in the kernel, so the scalar kernel identities
    # imply the operator ones; the operator-level check below confirms it too.
    for wtxt in ("0.05", "0.2", "0.9", "2.5", "7", "20", "44.1"):
        w = ivf(wtxt)
        w2 = w ** 2
        for k in (1, 2, 3, 7, 40, 999, 5000):
            k2 = ivf(k) ** 2
            den = k2 + w2
            brk1 = 2 * w * k2 / den ** 2                   # d_w [[k]](w)
            brk2 = 2 * k2 * (k2 - 3 * w2) / den ** 3       # d_ww [[k]](w)
            phi = k2 / den ** 3
            psi = IONE / den ** 2
            kap = IONE / (k2 * den)
            chk("(I1) [[k]]'' - [[k]]'/w + 8 w^2 phi_x = 0  (iv, exact)",
                _contains0(brk2 - brk1 / w + 8 * w2 * phi), f"w={wtxt} k={k}")
            chk("(I1) A_xx kernel equals -2 phi_x", _contains0(
                (brk2 - brk1 / w) / (4 * w2) + 2 * phi), f"w={wtxt} k={k}")
            chk("(I2) [[k]]'/(2w) = [[k]]/w^2 - w^2 psi_x  (iv, exact)",
                _contains0(brk1 / (2 * w) - (IONE / den - w2 * psi)),
                f"w={wtxt} k={k}")
            chk("(I) [[k]]/w^2 = 1/k^2 - w^2 kappa_w  (iv, exact)",
                _contains0(IONE / den - (IONE / k2 - w2 * kap)), f"w={wtxt} k={k}")
    for k in (1, 5, 300):
        chk("phi_0 = psi_0 = kappa_0 = 1/k^4 (so A_xx(0) = -2 G(4))", all(
            _contains0(t[k] - IONE / ivf(k) ** 4)
            for t in (tab_phi_iv(0, 400), tab_psi_iv(0, 400),
                      tab_kappa_iv(0, 400), tab_gamma_iv(4, 400))), f"k={k}")
    for N in (5, 12):
        for wtxt in ("0.3", "3.0"):
            w = ivf(wtxt)
            KTg = 2 * N + 4
            k2t = [ivf(k) ** 2 for k in range(KTg + 1)]
            Hxx = head_iv([IZERO] + [2 * k2t[k] * (k2t[k] - 3 * w ** 2)
                                     / (k2t[k] + w ** 2) ** 3
                                     for k in range(1, KTg + 1)], N)
            Hx = head_iv([IZERO] + [2 * w * k2t[k] / (k2t[k] + w ** 2) ** 2
                                    for k in range(1, KTg + 1)], N)
            Phi = head_iv(tab_phi_iv(w, KTg), N)
            bad = mp.zero
            allzero = True
            for i in range(N):
                for j in range(N):
                    d = Hxx[i][j] - Hx[i][j] / w + 8 * w ** 2 * Phi[i][j]
                    bad = max(bad, absmax_ub(d))
                    allzero = allzero and _contains0(d)
            chk("(I1) at OPERATOR level: A_xx + 2 O[phi_x] ENCLOSES 0 entrywise "
                "(a containment, not a tolerance)",
                allzero, f"N={N} w={wtxt} worst half-width={mp.nstr(bad, 4)}")

    # ---- (d) structure coefficients vs. the matvec -------------------------
    for N, wtxt in ((40, "0.5"), (120, "3.0")):
        w = ivf(wtxt)
        wf = float(wtxt)
        tr = TrialIV(trial_vec_np(wf, N))
        for name, ft in (("H", tab_H_iv(w, 2 * N + 4)),
                         ("phi", tab_phi_iv(w, 2 * N + 4)),
                         ("psi", tab_psi_iv(w, 2 * N + 4))):
            q1 = tr.qform(ft)
            q2 = ivdot(tr.v, tr.matvec(ft, N))
            chk("structure-coefficient qform and the matvec quadratic form have "
                "a nonempty intersection (both enclose the exact value)",
                LB(q1) <= UB(q2) and LB(q2) <= UB(q1), f"{name} N={N} w={wtxt}")
            chk("both enclosures are tight (width < 1e-22)",
                UB(q1) - LB(q1) < mp.mpf("1e-22")
                and UB(q2) - LB(q2) < mp.mpf("1e-22"), f"{name} N={N} w={wtxt}")
        a, b = 0.75 * wf, 1.25 * wf
        for nm, tabf in (("phi", tab_phi_iv), ("psi", tab_psi_iv),
                         ("kappa", tab_kappa_iv)):
            cellq = tr.qform_cell(cell_table(tabf(b, 2 * N + 4), tabf(a, 2 * N + 4)))
            for t in (0.0, 0.13, 0.5, 0.87, 1.0):
                wt = a + t * (b - a)
                q = tr.qform(tabf(wt, 2 * N + 4))
                chk("the cell-uniform interval table brackets the exact "
                    "quadratic form at every interior w",
                    LB(cellq) <= LB(q) and UB(q) <= UB(cellq),
                    f"{nm} N={N} w={wtxt} t={t}")
        lam = float(np.linalg.eigvalsh(op_np(N, np_tab(
            2 * N + 2, lambda k: wf ** 2 / (k * k + wf ** 2))))[-1])
        qq = tr.qform(tab_H_iv(w, 2 * N + 4))
        chk("evidence: the trial Rayleigh quotient matches the float64 lam_1",
            abs(float(LB(qq)) - lam) < 1e-11, f"N={N} w={wtxt}")
        chk("the Rayleigh quotient is a genuine LOWER bound (<= float64 lam_1 "
            "up to float64 error)", float(LB(qq)) <= lam + 1e-12,
            f"N={N} w={wtxt}")

    # ---- (e) envelopes -----------------------------------------------------
    # (e0) the envelope family itself: psi_j must dominate |d_w^j [[k]]|, both
    # pointwise and in the segment form.  (A wrong exponent here is exactly the
    # kind of silent hole this gate exists to catch.)
    for wtxt in ("0.05", "0.9", "7", "44.1"):
        w = ivf(wtxt)
        k2t = None
        for k in (1, 2, 5, 50, 1000):
            k2 = ivf(k) ** 2
            den = k2 + w ** 2
            ders = ((0, w ** 2 / den),
                    (1, 2 * w * k2 / den ** 2),
                    (2, 2 * k2 * (k2 - 3 * w ** 2) / den ** 3))
            for j, dv in ders:
                chk("psi_j dominates |d_w^j [[k]]| pointwise (the envelope "
                    "family, hence spec_point/spec_seg, is correct)",
                    UB(env_at(spec_point_iv(w, j), k)) >= absmax_ub(dv),
                    f"j={j} w={wtxt} k={k}")
        del k2t
    for (a, b) in ((0.25, 0.5), (1.0, 1.25), (10.0, 10.5), (40.0, 44.0)):
        for j in (0, 1, 2):
            sp = spec_seg_iv(a, b, j)
            for t in (0.0, 0.31, 0.77, 1.0):
                wv = ivf(a + t * (b - a))
                for k in (1, 3, 40, 700):
                    k2 = ivf(k) ** 2
                    den = k2 + wv ** 2
                    dv = {0: wv ** 2 / den,
                          1: 2 * wv * k2 / den ** 2,
                          2: 2 * k2 * (k2 - 3 * wv ** 2) / den ** 3}[j]
                    chk("the SEGMENT envelope spec_seg(a,b,j) dominates "
                        "|d_w^j [[k]](w)| simultaneously for every w in [a,b]",
                        UB(env_at(sp, k)) >= absmax_ub(dv),
                        f"j={j} [{a},{b}] t={t} k={k}")
    # -- the CERTIFIED arctangent (adversarial review M-1) -------------------
    # Three RIGOROUS self-consistency identities that use nothing but
    # iv_atan_nn itself, followed by one evidence-grade comparison against a
    # 600-bit reference.  The identities pin the halving loop AND the series
    # against each other at wildly different argument scales, which is what
    # catches a wrong constant / a wrong term count.
    pi2 = iv_pi_half()
    for ttxt in ("0.001", "0.3", "1", "2.5", "1000", "1e8"):
        tt = ivf(mp.mpf(ttxt))
        chk("atan(t) + atan(1/t) = pi/2  (t > 0) -- relates the small-argument "
            "series to the halving loop at the reciprocal scale",
            _contains0(iv_atan_nn(tt) + iv_atan_nn(IONE / tt) - pi2), f"t={ttxt}")
    for ttxt in ("0.05", "0.2", "0.4"):
        tt = ivf(mp.mpf(ttxt))
        chk("atan(2t/(1-t^2)) = 2 atan(t)  (|t| < 1) -- the doubling identity",
            _contains0(iv_atan_nn(2 * tt / (IONE - tt ** 2)) - 2 * iv_atan_nn(tt)),
            f"t={ttxt}")
    machin = 4 * (4 * iv_atan_nn(IONE / ivf(5)) - iv_atan_nn(IONE / ivf(239)))
    chk("Machin: pi = 16 atan(1/5) - 4 atan(1/239) agrees with 2*pi/2",
        _contains0(machin - 2 * pi2))
    chk("the certified atan is NARROW (width <= 1e-27 relative), i.e. the "
        "halving + alternating-series remainder costs nothing",
        UB(iv_atan_nn(ivf(mp.mpf("0.7")))) - LB(iv_atan_nn(ivf(mp.mpf("0.7"))))
        <= mp.mpf("1e-27"))
    _bad_at = 0
    for _i in range(1, 241):                     # EVIDENCE (typo insurance)
        _t = mp.mpf(_i) / 40 * mp.mpf(10) ** ((_i % 11) - 5)
        _g = iv_atan_nn(ivf(_t))
        with mp.workprec(600):
            _r = mp.atan(_t)
        if not (LB(_g) <= _r <= UB(_g)):
            _bad_at += 1
    chk("EVIDENCE: iv_atan_nn encloses a 600-bit reference on 240 operands "
        "spanning 10 decades (mpmath's own iv.atan2 fails this: its mpf_atan "
        "runs at prec+4 with round-to-NEAREST and only directs the LAST step)",
        _bad_at == 0, f"{_bad_at} misses")

    # -- REGRESSION GUARD: no mpmath transcendental on the certified path -----
    # M-1 was exactly "a library primitive that is not outward rounded slipped
    # into a load-bearing bound".  The only durable defence is a source-level
    # one: tokenise THIS file (so that prose in comments/docstrings cannot
    # trigger or mask it) and assert that Secs 0-6 -- everything that produces a
    # bound -- calls no transcendental from `mp` or `iv`.  Stage A' (Sec.7) and
    # the drivers are exempt: they only ever use them as EVIDENCE references.
    import io
    import tokenize as _tk
    _src = open(os.path.abspath(__file__), "r", encoding="utf-8").read()
    _lines = _src.splitlines()
    _s0 = next(i for i, L in enumerate(_lines) if "SEC. 0." in L)
    _s7 = next(i for i, L in enumerate(_lines) if "SEC. 7." in L)
    _BANNED = {"pi", "e", "atan", "atan2", "asin", "acos", "sin", "cos", "tan",
               "sinh", "cosh", "tanh", "log", "log10", "exp", "cbrt",
               "power", "quad", "zeta", "gamma", "expm1", "log1p", "hypot",
               "besselj", "polylog", "nsum", "findroot", "cholesky"}
    # `iv.sqrt` is the ONE library transcendental that is genuinely outward
    # rounded: mpmath's `mpi_sqrt` calls `mpf_sqrt(a, prec, round_floor)` and
    # `mpf_sqrt(b, prec, round_ceiling)`, i.e. the direction is applied to the
    # square root itself, not to a nearest-rounded intermediate (contrast
    # `mpf_atan2`).  It is reached only through `iv_sqrt_nn`.  `mp.sqrt` (plain
    # mpf, round-to-nearest) stays banned.
    _BAN = {"mp": _BANNED | {"sqrt"}, "iv": _BANNED}
    _toks = [t for t in _tk.generate_tokens(io.StringIO(_src).readline)
             if t.type in (_tk.NAME, _tk.OP)]
    _hits = []
    for _i in range(len(_toks) - 2):
        _a, _b, _c = _toks[_i], _toks[_i + 1], _toks[_i + 2]
        if (_a.type == _tk.NAME and _a.string in ("mp", "iv")
                and _b.string == "." and _c.type == _tk.NAME
                and _c.string in _BAN[_a.string]
                and _s0 < _a.start[0] - 1 < _s7):
            _hits.append(f"{_a.string}.{_c.string}@{_a.start[0]}")
    chk("no mpmath transcendental (pi/atan/log/exp/mp.sqrt/quad/zeta/cholesky) "
        "is called anywhere in Secs 0-6, i.e. on the certified path: the only "
        "primitives producing a bound are +,-,*,/,**int, the outward-rounded "
        "iv.sqrt (via iv_sqrt_nn) and the hand-rolled iv_atan_nn/iv_pi_half",
        not _hits, "; ".join(_hits))
    _bad_sq = 0
    for _i in range(1, 301):                     # the same test that killed atan2
        _t = mp.mpf(_i) / 7 * mp.mpf(10) ** ((_i % 13) - 6)
        _g = iv_sqrt_nn(ivf(_t))
        with mp.workprec(600):
            _r = mp.sqrt(_t)
        if not (LB(_g) <= _r <= UB(_g)):
            _bad_sq += 1
    chk("iv.sqrt (the only library primitive left on the certified path) "
        "encloses a 600-bit reference on 300 operands spanning 13 decades",
        _bad_sq == 0, f"{_bad_sq} misses")

    # -- the two derivative-envelope CONSTANTS 4b and 6b (review blind spot 1)
    # cell_L_iv bounds  |d_w psi_w(k)| <= 4b/(k^2+a^2)^3  and
    #                   |d_w phi_w(k)| <= 6b/(k^2+a^2)^3  for all w in [a,b].
    # Both are sharp at the corner w = b, a = b, so a wrong constant would slip
    # through any loose test; check them by exact interval domination.
    for (af, bf) in ((0.3125, 0.375), (1.0, 1.125), (9.9375, 10.0625),
                     (43.43066, 43.62891)):
        a_i, b_i = ivf(af), ivf(bf)
        wcell = ivlh(a_i, b_i)
        k2c = ivf(0) ** 2
        for k in (1, 2, 3, 5, 17, 140, 3000):
            k2 = ivf(k) ** 2
            den = k2 + wcell ** 2
            dpsi = 4 * wcell / den ** 3                  # |d_w (k^2+w^2)^-2|
            dphi = 6 * wcell * k2 / den ** 4             # |d_w k^2(k^2+w^2)^-3|
            majo = IONE / (k2 + a_i ** 2) ** 3
            chk("|d_w psi_w(k)| <= 4b/(k^2+a^2)^3 uniformly on the cell",
                UB(4 * b_i * majo) >= absmax_ub(dpsi), f"[{af},{bf}] k={k}")
            chk("|d_w phi_w(k)| <= 6b/(k^2+a^2)^3 uniformly on the cell",
                UB(6 * b_i * majo) >= absmax_ub(dphi), f"[{af},{bf}] k={k}")
        del k2c
    chk("tail_int(A=0, p=2, x=5) = 1/375",
        _contains0(tail_int(Spec(IONE, 0, 2, ("t1",)), 5) - IONE / ivf(375)))
    chk("tail_int(A=0, p=1, x=7) = 1/7",
        _contains0(tail_int(Spec(IONE, 0, 1, ("t2",)), 7) - IONE / ivf(7)))
    s1 = S1_of(Spec(IONE, 0, 2, ("z4gate",)))
    z4 = mp.zeta(4)
    chk("S1_of((1,0,2)) is an UPPER bound on zeta(4) and is tight",
        UB(s1) >= z4 and UB(s1) <= z4 + mp.mpf("1e-10"),
        f"{mp.nstr(UB(s1), 20)} vs zeta(4) = {mp.nstr(z4, 20)}")
    s1b = S1_of(Spec(IONE, 0, 1, ("z2gate",)))
    chk("S1_of((1,0,1)) is an UPPER bound on zeta(2)", UB(s1b) >= mp.zeta(2),
        f"{mp.nstr(UB(s1b), 20)}")
    for (C, A, p) in ((1, 2, 2), (3, mp.mpf("0.7"), 3), (9, 11, 1)):
        sp = Spec(C, A, p, ("q", str(C), str(A), p))
        for xx in (1, 7, 200):
            ref = mp.quad(lambda t: mp.mpf(C) / (t ** 2 + mp.mpf(A) ** 2) ** p,
                          [xx, mp.inf])
            v = tail_int(sp, xx)
            chk("tail_int closed form = high-precision quadrature",
                LB(v) <= ref * (1 + mp.mpf("1e-20"))
                and UB(v) >= ref * (1 - mp.mpf("1e-20")),
                f"C={C} A={A} p={p} x={xx}")
    for a in ((0.0, 1.0, 20.0) if full else (1.0,)):
        R = rsum(spec_maj_iv(a))
        R2 = rsum(spec_maj_iv(a), n0=120, M=2000)
        chk("rsum is stable under refining (n0, M) -- the cut is not load-bearing",
            abs(UB(R) - UB(R2)) <= UB(R) * mp.mpf("1e-5"),
            f"a={a} {mp.nstr(UB(R), 10)} vs {mp.nstr(UB(R2), 10)}")
        wv = max(a, 0.05)
        Nb = 500
        for nm, fn in (("psi", lambda k: 1.0 / (k * k + wv * wv) ** 2),
                       ("phi", lambda k: k * k / (k * k + wv * wv) ** 3)):
            nn = float(np.max(np.abs(np.linalg.eigvalsh(
                op_np(Nb, np_tab(2 * Nb + 2, fn))))))
            chk("evidence: rsum(majorant) >= ||O[f]|| at N=500",
                float(UB(R)) >= nn, f"a={a} {nm} rsum={float(UB(R)):.4e} "
                f"norm={nn:.4e}")
    for (a, b) in ((0.3125, 0.375), (10.0, 10.25)):
        L1 = rsum_seg(a, b, 1)
        L1f = rsum_seg(a, b, 1, Mrows=240, Mcols=1600)
        chk("rsum_seg is stable under refining (Mrows, Mcols)",
            abs(UB(L1) - UB(L1f)) <= UB(L1) * mp.mpf("2e-3"),
            f"[{a},{b}] {mp.nstr(UB(L1), 8)} vs {mp.nstr(UB(L1f), 8)}")
        Nb = 400
        h = 1e-6
        dH = (op_np(Nb, np_tab(2 * Nb + 2, lambda k: (b + h) ** 2 / (k * k + (b + h) ** 2)))
              - op_np(Nb, np_tab(2 * Nb + 2, lambda k: (b - h) ** 2 / (k * k + (b - h) ** 2)))) / (2 * h)
        nn = float(np.max(np.abs(np.linalg.eigvalsh(dH))))
        chk("evidence: L1 >= ||d_w H(b)|| (float64 central difference, N=400)",
            float(UB(L1)) >= nn, f"[{a},{b}] {float(UB(L1)):.5f} vs {nn:.5f}")
    # the grid heuristic must track the certified envelope (a wrong exponent in
    # _L1_hint would silently produce coarse cells whose Lipschitz slack ruins
    # alpha -- it cannot make a bound unsound, but it must not be allowed to
    # drift unnoticed)
    for w in (0.3125, 1.0, 5.0, 44.0):
        hint = _L1_hint(w)
        cert = float(UB(rsum_seg(w, w * 1.02, 1)))
        chk("the float64 grid heuristic _L1_hint tracks the certified "
            "rsum_seg(w, 1.02w, 1) to within a factor of 1.5",
            0.5 * cert <= hint <= 1.5 * cert,
            f"w={w} hint={hint:.5f} certified={cert:.5f}")

    # ---- (f) mu_2 and the head certificate ---------------------------------
    for (a, b) in ((0.3125, 0.375), (1.0, 1.125), (10.0, 10.25), (44.0, 44.75)):
        m2, ok = mu2_cell_iv(a, b)
        chk("mu_2 segment monotonicity k1, k2, k3 >= 0", ok, f"[{a},{b}]")
        lam_true = float(np.linalg.eigvalsh(op_np(400, np_tab(
            802, lambda k: a * a / (k * k + a * a))))[-1])
        chk("evidence: mu_2 <= lam_1(H(a)) (N=400 float64)",
            float(m2) <= lam_true, f"[{a},{b}] {float(m2):.6f} vs {lam_true:.6f}")

    # ---- (g) cross-check against the mpf certificate ------------------------
    # EVIDENCE ONLY (typo insurance).  cell_L_iv is fed EXACTLY the same alpha,
    # L1, mu2 as the base script so that the comparison isolates the (C**)
    # assembly itself; the interval alpha is compared separately, where the only
    # expected difference is the tighter max(.,.) block form.
    if full:
        try:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            import p4_extremal_concavity as C0
            from p4_interchange_uniformity import alpha_certificate_H, mu2_cell_lb
            from p4_interval_extension import schur_full_seg_w0
            have = True
        except Exception as exc:                                   # noqa: BLE001
            print(f"  (cross-check against the mpf script skipped: {exc})")
            have = False
    else:
        have = False
    if have:
        for (A, B) in ((0.3125, 0.375), (1.0, 1.125), (9.9375, 10.0625)):
            Na = na_band(A)
            L1_ref = schur_full_seg_w0(A, B, 1)
            L1_iv = UB(rsum_seg(A, B, 1))
            chk("evidence: the iv L1 is >= the numpy L1 (same algebra, looser "
                "column cut)", L1_iv >= L1_ref * mp.mpf("0.999"),
                f"[{A},{B}] {mp.nstr(L1_iv, 8)} vs {mp.nstr(L1_ref, 8)}")
            chk("evidence: the iv L1 is within 1% of the Mcols=8000 value",
                L1_iv <= L1_ref * mp.mpf("1.01"), f"[{A},{B}]")
            al0, ok0, c0, nB0, nD0 = alpha_certificate_H(A, Na=Na)
            ftab = tab_H_iv(A, 2 * Na + 4)
            evh, vech = np.linalg.eigh(op_np(Na, np_tab(
                2 * Na + 2, lambda k: A * A / (k * k + A * A))))
            c_iv, pivot, tries = lam2_head_iv(ftab, Na, vech[:, -1],
                                              c0=float(evh[-2]))
            chk("interval Cholesky certifies the EXACT head at a coarse node",
                pivot > 0, f"w={A} Na={Na} pivot={mp.nstr(pivot, 6)} "
                f"tries={tries}")
            chk("the interval c really bounds lam_2 of the float64 head from "
                "above (it is a certified bound for the EXACT head; this only "
                "confirms no gross transcription error)", c_iv >= evh[-2],
                f"w={A} c_iv={mp.nstr(c_iv, 6)} lam2_f64={evh[-2]:.6e}")
            chk("evidence: the interval c does not overshoot the mpf c by more "
                "than the escalation factor", c_iv <= 4 * c0 + mp.mpf("1e-11"),
                f"w={A} c_iv={mp.nstr(c_iv, 4)} c_mpf={mp.nstr(c0, 4)}")
            mu0, _, _, _ = mu2_cell_lb(A, B)
            mu_iv, okm = mu2_cell_iv(A, B)
            chk("evidence: the iv mu_2 agrees with the mpf mu_2 to the ROUND_PAD",
                okm and abs(mu_iv - mu0) < mp.mpf("1e-11"),
                f"[{A},{B}] {mp.nstr(mu_iv, 12)} vs {mp.nstr(mu0, 12)}")
            alpha_0 = C0.alpha_cell(A, B, c0, c0, L1_ref, Na)
            alpha_iv = alpha_cell_iv(A, B, c_iv, c_iv, L1_iv, Na)
            chk("the interval alpha (max-block form) is <= the mpf alpha "
                "(sum-block form), as it must be", alpha_iv <= alpha_0,
                f"[{A},{B}] {mp.nstr(alpha_iv, 10)} vs {mp.nstr(alpha_0, 10)}")
            chk("evidence: the two alphas are of the same order (the only "
                "structural difference is max(.,.) vs the sum of the SAME two "
                "blocks)", alpha_iv >= alpha_0 / 3,
                f"[{A},{B}] {mp.nstr(alpha_iv, 8)} vs {mp.nstr(alpha_0, 8)}")
            # identical inputs => the assembly must agree to ~ROUND_PAD
            ub_iv, info_iv = cell_L_iv(A, B, alpha_0, L1_ref, mu0)
            ub_0, info_0 = C0.cell_L(A, B, alpha_0, L1_ref, mu0)
            chk("the iv cell certificate closes wherever the mpf one does",
                ub_iv is not None and ub_iv < 0 and ub_0 < 0, f"[{A},{B}]")
            for kname in ("F_lo", "E_hi", "g_lo", "s", "s_p", "theta", "rnorm",
                          "ehat_phi", "ehat_psi", "N", "pos", "hpp_ub"):
                r0, r1 = info_0[kname], info_iv[kname]
                chk("evidence: iv vs mpf, quantity by quantity, same inputs",
                    abs(r1 - r0) <= 1e-6 * abs(r0) + 1e-11,
                    f"[{A},{B}] {kname} mpf={r0:.12e} iv={r1:.12e}")
            # and the certificate that will actually be used (tighter alpha)
            ub_use, info_use = cell_L_iv(A, B, alpha_iv, L1_iv, mu_iv)
            chk("the production iv certificate (tighter alpha, iv L1) closes",
                ub_use is not None and ub_use < 0, f"[{A},{B}]")
            if A == 1.0:
                # ---- safety monotonicity in the three INPUTS (blind spot 5) --
                # cell_L_iv is fed alpha, L1, mu2 from OTHER certificates.  If
                # it were not monotone the right way, a slightly-too-generous
                # upstream bound could make hpp_ub LOOK better.  Degrading each
                # input must never improve the conclusion.
                for (lbl, aa, ll, mm) in (
                        ("alpha up", alpha_0 * mp.mpf("1.05"), L1_ref, mu0),
                        ("L1 up", alpha_0, L1_ref * mp.mpf("1.05"), mu0),
                        ("mu_2 down", alpha_0, L1_ref, mu0 * mp.mpf("0.95"))):
                    ub_d, _ = cell_L_iv(A, B, aa, ll, mm)
                    chk("degrading an INPUT bound never improves hpp_ub "
                        "(cell_L_iv is monotone in the safe direction)",
                        ub_d is None or ub_d >= ub_iv,
                        f"[{A},{B}] {lbl}: {ub_d} vs {ub_iv}")
            print(f"  cross-check [{A},{B}]: hpp_ub mpf {ub_0:.8e} / iv "
                  f"{float(ub_iv):.8e} (ratio {float(ub_iv) / ub_0:.8f}); "
                  f"production iv {float(ub_use):.8e}, pos/F "
                  f"{info_0['ratio']:.4f} -> {info_use['ratio']:.4f}")
    if have:
        # ---- (g2) the SAME cross-check for Stage S (rescaled operator) ------
        # cell_S_iv takes no external inputs, so this compares the two Stage S
        # pipelines end to end (G(2) trial vector, Davis-Kahan on h~, the
        # kappa/phi/psi quadratic forms).  EVIDENCE ONLY.
        for (a, b) in ((0.0, 0.3125), (0.0, 0.25)):
            ubS_iv, iS_iv = cell_S_iv(a, b)
            ubS_0, iS_0 = C0.cell_S(a, b)
            chk("the iv Stage S certificate closes wherever the mpf one does",
                ubS_iv is not None and ubS_iv < 0 and ubS_0 < 0, f"[{a},{b}]")
            for kname in ("F_lo", "E_hi", "g_lo", "s", "theta", "alpha", "qk",
                          "ehat_phi", "ehat_psi", "N", "pos", "hpp_ub"):
                r0, r1 = iS_0[kname], iS_iv[kname]
                chk("evidence: Stage S iv vs mpf, quantity by quantity",
                    abs(r1 - r0) <= 1e-6 * abs(r0) + 1e-11,
                    f"[{a},{b}] {kname} mpf={r0:.12e} iv={r1:.12e}")
            chk("the iv Stage S bound is CONSERVATIVE w.r.t. the mpf one "
                "(interval arithmetic may only give up, never gain)",
                float(ubS_iv) >= ubS_0 - 1e-11,
                f"[{a},{b}] iv={float(ubS_iv):.12e} mpf={ubS_0:.12e}")
            print(f"  cross-check S [{a},{b}]: hpp_ub mpf {ubS_0:.8e} / iv "
                  f"{float(ubS_iv):.8e} (ratio {float(ubS_iv) / ubS_0:.8f})")
    _summary("Stage A'", f"({time.time() - t0:.1f}s)")


# ===========================================================================
# SEC. 8.  THE GRID: EXACT DYADIC ENDPOINTS AND AN EXPLICIT TILING CHECK
# ===========================================================================

COARSE_SLACK = 0.35
COARSE_REL_MAX = 0.20
MAX_SPLIT = 6
DYAD_COARSE = 1024          # coarse nodes are integer multiples of 2^-10
DYAD_ALL = 65536            # every endpoint is an integer multiple of 2^-16


def is_dyadic(x, q=DYAD_ALL):
    """True iff x is an integer multiple of 1/q (q a power of two) and small
    enough that the grid arithmetic is exact in float64."""
    y = float(x) * q
    return y == float(int(y)) and abs(y) < 2.0 ** 50


def _L1_hint(w):
    """float64 heuristic for sup||d_w H|| near w, used ONLY to CHOOSE the coarse
    step.  It is NOT a bound: the L1 that enters alpha is recomputed on the
    chosen cell by `rsum_seg` (interval), and the tiling check does not look at
    this number at all.  A faithful float64 transcription of the same row sums
    (minus the integral tails, which only make the hint slightly smaller and
    hence the step slightly larger -- absorbed by adaptive subdivision).
    The envelope triple is READ OFF `spec_seg_iv` rather than re-typed, so the
    hint cannot drift out of sync with the certified envelope."""
    a, b = float(w), float(w) * 1.02
    sp = spec_seg_iv(a, b, 1)
    C, A, p = float(UB(sp.C)), float(UB(sp.A)), sp.p
    Mrows, Mcols = 60, 900
    m = np.arange(Mcols, dtype=np.float64)
    invm = 1.0 / np.sqrt(2.0 * m + 1.0)
    kk = np.arange(Mrows + Mcols + 2, dtype=np.float64)
    psi_k = np.zeros_like(kk)
    psi_k[1:] = C / (kk[1:] ** 2 + A ** 2) ** p
    cum = np.concatenate(([0.0], np.cumsum(psi_k[1:])))
    mi = m.astype(np.int64)
    best = 0.0
    for n in range(Mrows):
        d = np.abs(mi - n)
        h2 = psi_k[d] * invm
        h2[n] = 0.0
        h3 = psi_k[mi + n + 1] * invm
        r = (2.0 * cum[n] / (2 * n + 1)
             + (h2.sum() + h3.sum()) / np.sqrt(2.0 * n + 1.0))
        best = max(best, r)
    return best


def coarse_grid(W0, WMAX):
    """Coarse nodes are integer multiples of 2^-10, so consecutive coarse cells
    abut EXACTLY and the grid is bit-for-bit reproducible: the float64 heuristic
    only has to land in the same 2^-10 bucket, and nothing downstream trusts
    it."""
    chk("coarse grid: the endpoints of the ladder are dyadic",
        is_dyadic(W0, DYAD_COARSE) and is_dyadic(WMAX, DYAD_COARSE),
        f"{W0!r} {WMAX!r}")
    g = [float(W0)]
    while g[-1] < WMAX:
        w = g[-1]
        d = min(2 * COARSE_SLACK / max(_L1_hint(w), 1e-9), COARSE_REL_MAX * w)
        dq = max(1, int(d * DYAD_COARSE)) / DYAD_COARSE       # floor to 2^-10
        nxt = min(WMAX, w + dq)
        chk("coarse grid: the step is strictly positive (no stall)", nxt > w,
            f"w={w!r} d={d!r}")
        g.append(nxt)
    for t in g:
        chk("coarse grid: every node is an exact multiple of 2^-10",
            is_dyadic(t, DYAD_COARSE), f"{t!r}")
    return g


def subdivide(A, B, k):
    """[A,B] into k equal subcells with EXACT shared endpoints.  A and B are
    multiples of 2^-10 and k is a power of two <= 2^MAX_SPLIT, so every endpoint
    is a multiple of 2^-16 and np.linspace reproduces it exactly (it also pins
    the last entry to B)."""
    assert k & (k - 1) == 0 and k <= 2 ** MAX_SPLIT
    sub = [float(t) for t in np.linspace(A, B, k + 1)]
    chk("fine grid: the subdivision spans the coarse cell exactly",
        sub[0] == A and sub[-1] == B, f"[{A!r},{B!r}] k={k}")
    for t in sub:
        chk("fine grid: every endpoint is an exact multiple of 2^-16",
            is_dyadic(t), f"{t!r}")
    return sub


def check_tiling(cells, W0, WMAX):
    """The certified cells must tile [W0, WMAX] with no gap and no float slop."""
    chk("tiling: the ladder starts exactly at wbar_0", cells[0][0] == W0,
        f"{cells[0][0]!r} vs {W0!r}")
    chk("tiling: the ladder ends exactly at the reported W", cells[-1][1] == WMAX,
        f"{cells[-1][1]!r} vs {WMAX!r}")
    for i in range(len(cells) - 1):
        chk("tiling: consecutive fine cells share a bit-identical endpoint",
            cells[i][1] == cells[i + 1][0],
            f"{cells[i][1]!r} vs {cells[i + 1][0]!r}")
    for (a, b) in cells:
        chk("tiling: every cell is nondegenerate with dyadic endpoints",
            is_dyadic(a) and is_dyadic(b) and a < b, f"[{a!r},{b!r}]")


# ===========================================================================
# SEC. 9.  DRIVERS
# ===========================================================================

def _do_node(arg):
    """Certify c >= lam_2(head of H(w)) at one coarse node (worker entry)."""
    w, Na = arg
    n0 = PASS[0]
    ftab = tab_H_iv(w, 2 * Na + 4)
    ev, vec = np.linalg.eigh(op_np(Na, np_tab(
        2 * Na + 2, lambda k: w * w / (k * k + w * w))))
    c, pivot, tries = lam2_head_iv(ftab, Na, vec[:, -1], c0=float(ev[-2]))
    chk("interval Cholesky: cI - PAP is positive definite (coarse node)",
        pivot > 0, f"w={w} Na={Na}")
    nch = PASS[0] - n0
    # The PARENT adds nch (see run_L).  Resetting the worker's own counter makes
    # the serial path (same process) and the pooled path (child process, counter
    # discarded) report EXACTLY the same number of checks.
    PASS[0] = n0
    return (w, Na, as_tuple(c), as_tuple(pivot), tries, nch)


def _do_cell(arg):
    """Certify one coarse cell: L1, alpha, mu_2, then adaptive fine cells."""
    A, B, Na, ctA, ctB = arg
    n0 = PASS[0]
    cA, cB = from_tuple(ctA), from_tuple(ctB)
    L1 = UB(rsum_seg(A, B, 1))
    chk("segment envelope L1 = sup||d_w H|| is finite and positive",
        0 < L1 < mp.mpf("1e6"), f"[{A},{B}]")
    alpha = alpha_cell_iv(A, B, cA, cB, L1, Na)
    mu2, okm = mu2_cell_iv(A, B)
    chk("mu_2 segment monotonicity k1, k2, k3 >= 0", okm, f"[{A},{B}]")
    k, got = 1, None
    while k <= 2 ** MAX_SPLIT:
        sub = subdivide(A, B, k)
        acc, bad = [], False
        for j in range(k):
            ub, info = cell_L_iv(sub[j], sub[j + 1], alpha, L1, mu2)
            if ub is None or ub >= 0 or info["ratio"] > TARGET:
                bad = True
                break
            acc.append(info)
        if not bad:
            got = acc
            break
        k *= 2
    if got is not None:
        for info in got:
            chk("fine cell: H''(x) < 0 certified for all x in [a^2, b^2]",
                info["hpp_ub"] < 0, f"[{info['a']:.6f},{info['b']:.6f}]")
    nch = PASS[0] - n0
    PASS[0] = n0                            # see _do_node
    return (A, B, got, nch, float(L1), float(alpha), float(mu2))


def _pool_map(fn, args, jobs):
    if jobs <= 1 or len(args) <= 1:
        return [fn(a) for a in args]
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        return list(ex.map(fn, args, chunksize=1))


LAD = {}


def run_L(W0, WMAX, jobs=1, tag="Stage L"):
    print(f"=== {tag}: interval concavity ladder on [{W0}, {WMAX}] "
          f"(jobs={jobs}) ===")
    t0 = time.time()
    coarse = coarse_grid(W0, WMAX)
    nas = [min(na_band(coarse[i]), na_band(coarse[i + 1]))
           for i in range(len(coarse) - 1)]
    need = sorted({(coarse[i], nas[i]) for i in range(len(nas))}
                  | {(coarse[i + 1], nas[i]) for i in range(len(nas))})
    print(f"  coarse grid: {len(coarse)} dyadic nodes, {len(need)} (node, Na) "
          f"head certificates, Na in {sorted(set(nas))}; fine cells adaptive "
          f"(cap {MAX_SPLIT} doublings, target pos/F <= {TARGET})")

    t1 = time.time()
    res = _pool_map(_do_node, need, jobs)
    cmap, pivmin, tries_max = {}, None, 0
    for (w, Na, ct, pt, tries, nch) in res:
        cmap[(w, Na)] = ct
        PASS[0] += nch
        p = from_tuple(pt)
        pivmin = p if pivmin is None else min(pivmin, p)
        tries_max = max(tries_max, tries)
    print(f"  head certificates: {len(need)} interval Choleskys in "
          f"{time.time() - t1:.1f}s; smallest certified pivot = "
          f"{mp.nstr(pivmin, 6)} > 0; max ladder escalations = {tries_max}")

    t2 = time.time()
    cargs = [(coarse[i], coarse[i + 1], nas[i],
              cmap[(coarse[i], nas[i])], cmap[(coarse[i + 1], nas[i])])
             for i in range(len(nas))]
    out = _pool_map(_do_cell, cargs, jobs)
    rows, cells, worst, stopped = [], [], None, None
    for (A, B, got, nch, L1, alpha, mu2) in out:
        PASS[0] += nch
        if got is None:
            stopped = A
            print(f"  !! coarse cell [{A:.5f},{B:.5f}] did NOT close after "
                  f"{MAX_SPLIT} doublings; the ladder stops here.")
            break
        for info in got:
            rows.append(info)
            cells.append((info["a"], info["b"]))
            if worst is None or info["ratio"] > worst["ratio"]:
                worst = info
    print(f"  fine cells: {len(rows)} certified in {time.time() - t2:.1f}s")
    if not rows:
        raise AssertionError("Stage L closed no cell at all")
    Wdone = cells[-1][1]
    check_tiling(cells, float(W0), Wdone)
    # The two checks inside check_tiling that involve Wdone are internal
    # consistency only (Wdone is DEFINED as cells[-1][1]).  The claim that
    # actually matters -- that the ladder covered what it was asked to cover --
    # is this one, and it must abort rather than silently report a smaller W:
    chk("the ladder reached the REQUESTED W_max (no coarse cell was abandoned "
        "and the last certified endpoint is the target)",
        stopped is None and Wdone == float(WMAX),
        f"stopped at {stopped}, reached {Wdone!r}, requested {float(WMAX)!r}")
    LAD.update(rows=rows, W=Wdone, worst=worst, stopped=stopped, cells=cells,
               coarse=coarse, pivmin=pivmin)

    print(f"  {'band':>18} {'cells':>6} {'Ns':>5} {'max pos/F':>10} {'F_lo':>11} "
          f"{'pos':>11} {'g_lo':>6} | E_hi split %: {'ehat':>5} {'width':>5} "
          f"{'angle':>5} | {'s':>9} {'s_p':>9}")
    edges = sorted(set([float(W0)]
                       + [e for e in (0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0,
                                      30.0, 45.0, 60.0, 80.0, 100.0, 120.0)
                          if float(W0) < e < Wdone] + [Wdone]))
    for e0, e1 in zip(edges[:-1], edges[1:]):
        sel = [r for r in rows if r["a"] >= e0 - 1e-9 and r["b"] <= e1 + 1e-9]
        if not sel:
            continue
        wc = max(sel, key=lambda r: r["ratio"])
        tot = max(wc["E_hi"], 1e-300)
        print(f"  [{e0:6.3f},{e1:8.3f}] {len(sel):>6d} {wc['Ns']:>5d} "
              f"{wc['ratio']:>10.4f} {wc['F_lo']:>11.4e} {wc['pos']:>11.4e} "
              f"{wc['g_lo']:>6.3f} |               "
              f"{100 * wc['ehat_psi'] / tot:5.1f} "
              f"{100 * wc['lip_psi'] / tot:5.1f} {100 * wc['ang_psi'] / tot:5.1f}"
              f" | {wc['s']:>9.3e} {wc['s_p']:>9.3e}")
    if worst is not None:
        print(f"  binding fine cell [{worst['a']:.5f},{worst['b']:.5f}]: "
              f"pos/F_lo = {worst['ratio']:.4f}, F_lo = {worst['F_lo']:.5e}, "
              f"pos = {worst['pos']:.5e}, hpp_ub = {worst['hpp_ub']:.5e}")
        print(f"    E_hi = {worst['E_hi']:.4e} = ehat_psi {worst['ehat_psi']:.4e}"
              f" + width {worst['lip_psi']:.4e} + angle {worst['ang_psi']:.4e}; "
              f"s = {worst['s']:.3e} (s_p = {worst['s_p']:.3e}), g_lo = "
              f"{worst['g_lo']:.4f}, alpha = {worst['alpha']:.4f}")
    nn = sum(1 for r in rows if r.get("cnn"))
    print(f"  ==> H''(x) < 0 CERTIFIED (interval arithmetic) for x in "
          f"(0, {Wdone ** 2:.4f}]  (wbar up to {Wdone}); {len(rows)} fine "
          f"cells; c_k >= 0 on {nn}/{len(rows)}")
    _summary(tag, f"({time.time() - t0:.1f}s)")
    return Wdone


# dyadic replacements for the base script's (0.5, 0.45, 0.4, 0.35, 0.3, ...):
S_CAND = (0.5, 0.4375, 0.375, 0.3125, 0.25, 0.1875, 0.125)
SFIRST = {}


def run_S():
    print("=== STAGE S: the first cell (0, wbar_0], reaching x = 0 exactly ===")
    t0 = time.time()
    stage_S_setup()
    chosen = None
    for w0 in S_CAND:
        chk("Stage S candidate wbar_0 is dyadic", is_dyadic(w0, DYAD_COARSE),
            f"{w0!r}")
        ub, info = cell_S_iv(0.0, w0)
        if ub is not None and ub < 0 and info["ratio"] <= 0.85:
            chosen = (w0, ub, info)
            break
    chk("the first cell (0, wbar_0] closes for some dyadic wbar_0",
        chosen is not None)
    w0, ub, info = chosen
    SFIRST.update(w0=w0, ub=float(ub), info=info)
    print(f"  wbar_0 = {w0}  ==>  sup_(x in (0,{w0 ** 2:.5f}]) H''(x) <= "
          f"{float(ub):.8f} < 0")
    print(f"    F_lo = {info['F_lo']:.8f}, positive term = {info['pos']:.8f} "
          f"(pos/F = {info['ratio']:.4f}), sin angle <= {info['s']:.4e}, "
          f"g~_lo = {info['g_lo']:.8f}")
    chk("first cell: H''(x) < 0 on (0, wbar_0^2]", ub < 0)
    ub0, info0 = cell_S_iv(0.0, 0.0)
    chk("the degenerate cell at x = 0 closes", ub0 is not None and ub0 < 0)
    chk("the positive term vanishes identically at x = 0", info0["pos"] == 0.0)
    print(f"  x -> 0 limit: H''(0) <= {float(ub0):.10f}   (c_2 >= "
          f"{info0['F_lo']:.10f})")
    chk("H''(0) <= -2*1.193 (matches c_2 from Segment I)",
        ub0 <= mp.mpf("-2.386"), f"ub0={float(ub0)}")
    chk("the c_2 lower bound agrees with Segment I's "
        "[1.19312762268, 1.19320906584]",
        1.1931 <= info0["F_lo"] <= 1.19321, f"{info0['F_lo']}")
    _summary("Stage S", f"({time.time() - t0:.1f}s)")
    return w0


def main():
    args = sys.argv[1:]
    jobs = 1
    if "-j" in args:
        i = args.index("-j")
        jobs = int(args[i + 1])
        del args[i:i + 2]
    ext = "EXT" in args
    gate_only = "GATE" in args
    fast = "FAST" in args
    args = [a for a in args if a not in ("EXT", "GATE", "FAST")]
    WMAX = float(args[0]) if args else (45.0 if ext else 3.0)
    if WMAX > 45.0:
        NA_RULE[0] = "ext"
    print("P4 Segment C-I -- INTERVAL-ARITHMETIC certificate "
          f"(mpmath.iv, iv.dps = {iv.dps}, mp.prec = {mp.prec})")
    print(f"  target W_max = {WMAX},  Na rule = {NA_RULE[0]},  jobs = {jobs}")
    run_gate(full=not fast)
    if gate_only:
        print(f"\nGATE ONLY: {PASS[0]} checks passed in {time.time() - T0:.1f}s")
        return
    w0 = run_S()
    W = run_L(w0, WMAX, jobs=jobs)
    print()
    print(f"RESULT: H(x) = h(sqrt x) is strictly concave on (0, {W ** 2:.4f}] "
          f"(wbar in (0, {W}]), certified with outward-rounded interval "
          f"arithmetic end to end.")
    print(f"ALL {PASS[0]} checks passed  ({time.time() - T0:.1f} s)")


if __name__ == "__main__":
    main()
