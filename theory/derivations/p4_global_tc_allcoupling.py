"""P4 GLOBAL T_c ALL-COUPLING -- two-sided bracket for T_c/omega_E at EVERY
lambda>0 (Einstein spectrum, mu*=0), with the weak-coupling UPPER side upgraded
from "asymptotic" to PROVEN / CERTIFIED.

WHY THIS FILE EXISTS
--------------------
Deliverable (3) (`p4_global_tc_certified.py`) certifies the crossover band
lambda in [0.582, 34.4] and the strong band (34.4, 50].  Below lambda ~ 0.582 its
upper T_c side came from `_weak_asympt_invert`, i.e. the ASYMPTOTIC form
h ~ L + c0 + b1/L + REM/L^2 with an empirically calibrated pad REM=0.75.  That is
not a certificate.  This file replaces it.

TWO SOUNDNESS DEFECTS FOUND IN THE INHERITED WEAK BRACKET
---------------------------------------------------------
(D1) `secular_bounds_K` evaluates the Neumann partial sum S_K at a FINITE matrix
     size N.  Every term of S_K is a sum of NONNEGATIVE quantities (tau >= 0 and
     calR >= 0 entrywise), so the finite-N value UNDER-estimates the infinite-N
     value:  S_K^{(N)} <= S_K^{(inf)}.  Adding the Theorem-W tail to a quantity
     that is itself too small does NOT yield an upper bound on F(h), so
     root(S_K^{(N)} + tail) is not a certified upper bound on h.  (The defect is
     numerically ~1e-4 relative at N=100 and ~1e-9 at N=800 -- small, but it is
     a gap in the logic, not a certificate.)  The LOWER bound root(S_K^{(N)}) is
     unaffected and remains sound (see Prop. L below).
(D2) The deep-tail upper side was asymptotic, as documented above.

WHAT IS PROVED HERE
-------------------
Let  F(h) := 2 tau^T M(h)^{-1} tau,   M(h) := h + D - calR,  where (tau, D, calR)
are the (M1) Lorentz-addition data, D_n is an EXACT finite sum (no truncation),
calR >= 0 entrywise, and ||calR|| <= Theta* (certified, Lemma M).  h = h(w) is the
top eigenvalue iff F(h) = 1, and F is strictly decreasing in h.

Prop. L (lower, sound).  For the head projection P = {n < N},
    (M^{-1})_{PP} = (M_PP - M_PQ M_QQ^{-1} M_QP)^{-1}  >=  (M_PP)^{-1}
and M^{-1} >= 0 entrywise, so
    F(h) >= 2 tau_P^T (M_PP)^{-1} tau_P >= sum_{j=0}^{K} S_j^P(h).
Hence root(sum_j S_j^P) <= h.   [this is the inherited h_lo, now justified]

Prop. U (upper, new).  With  eps := kappa_N^2/(h - Theta*),  kappa_N >= ||calR_PQ||,
and  B_Q := sum_{n>=N} tau_n^2,
    F(h)/2 <= [ sqrt(H(h - eps)) + sqrt(B_Q/(h - Theta*)) ]^2,
    H(x)   := tau_P^T ((x + D_P - calR_PP))^{-1} tau_P
           <= (1/2) [ sum_{j=0}^{K} S_j^P(x) + S_0^P(x) r^{K+1}/(1-r) ],
           r  := Theta*/x,
because M_QQ >= h - Theta*, the Schur complement obeys
M_PP - M_PQ M_QQ^{-1} M_QP >= M_PP - eps I, and the head geometric remainder uses
S_j^P <= r^j S_0^P (conjugate by (x+D_P)^{-1/2}; ||calR_PP|| <= Theta* since a
compression cannot increase the norm).  ALL quantities are finite-N head sums, so
(D1) does not arise.

Prop. A (all-lambda, proven, loose).  calR <= Theta* I gives, for every w>0,
    F(h) <= U(h) := sum_n 2 tau_n^2 / (h + D_n^lo - Theta*),
with D_n^lo = 1/(2n+1) + (2/(2n+1)) w [atan((n+1)/w) - atan(1/w)]  (G decreasing)
and a closed-form n-tail  sum_{n>=M} tau_n^2 <= w^4/(8 (M-1/2)^4).  U is a 1-D sum,
so it is computable at ANY w -- including w ~ 1e8 where no matrix fits.  Paired
with the PROVEN Cor (iii) lower  h >= ln w + c0 - 1/(12w) - 3/(16w^2), this closes
the bracket for ALL lambda > 0.  Cost: U overshoots h by ~0.52-0.59, i.e. ~55%
relative width in T_c.  Loose, but PROVEN, and it is what makes the "all lambda"
claim honest.

BAND MAP (see the companion .md for the theorem statement)
    strong  (34.4, 50]     Weyl envelope             PROVEN analytic
    crossover [0.582,34.4] Rayleigh + Kato-Temple    CERTIFICATE (node_cert)
    weak-W  [~0.20,0.582]  Prop. L + Prop. U         CERTIFICATE (finite-N)
    deep    (0, ~0.20]     Cor (iii) + Prop. B        PROVEN (modulo calibrated C)

The matrix route CANNOT reach lambda=0.05: that needs h=20, i.e. w ~ e^19.4 ~ 3e8,
so N ~ 1e10.  This is a structural limit, reported honestly, and is exactly why
Prop. A is needed.

Run:  OMP_NUM_THREADS=1 python p4_global_tc_allcoupling.py [N=400] [K=8]
"""
import os
import sys
import time
import numpy as np
import mpmath as mp

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
from p4_interchange_uniformity import H_mat            # noqa: E402
from p4_global_tc_certified import node_cert           # noqa: E402

mp.mp.dps = 30

PASS = [0]
FAIL = [0]


def ck(cond, msg=""):
    if cond:
        PASS[0] += 1
    else:
        FAIL[0] += 1
        print(f"  !! CHECK FAILED: {msg}")
    return bool(cond)


# --------------------------------------------------------------------------
# certified constants
# --------------------------------------------------------------------------
THETA = 0.6349287917881079          # >= ||calR||               (Lemma M)
# g(2) enclosure.  CERT-D (P4_g2_analytic.md, interval arithmetic, width
# 6.1e-10) supersedes the CERT-B bracket (width 4.0e-5) used in v1 of this
# file; CERT-B is retained ONLY as a nested-enclosure sentinel (section [0]).
G2_LO = mp.mpf("1.318140211762671")   # <= g(2)                 (CERT-D)
G2_HI = mp.mpf("1.318140212373478")   # >= g(2)                 (CERT-D)
G2_LO_B = mp.mpf("1.3181093293757127")  # <= g(2)               (CERT-B, sentinel)
G2_HI_B = mp.mpf("1.3181615784372800")  # >= g(2)               (CERT-B, sentinel)
Z4 = mp.zeta(4)
EG = 0.5772156649015329
C0 = EG + np.log(4.0) - 1.5          # gamma + ln4 - 3/2
TWO_PI = float(2 * mp.pi)
FL = 1e-12                           # float64 rounding pad (generous)

# Certified pointwise bounds on the Schur row sums c_n = sum_m R(s_n,s_m)/(2m+1).
# NOTE these bound c_n, which is NOT bounded by Theta* (=||calR||); measured
# max_n c_n = 0.6372 > Theta* = 0.63493.  Calibrated as sup over w in [2,400]
# and all n, then padded.  c_n depends on n only through a = s_n = (n+1/2)/w,
# and the map a -> c is w-INDEPENDENT in the continuum limit (verified: the
# discrete c_n matches (1/2)\int_0^inf R(a,s)/s ds to 2e-4 at w=100,1000),
# which is what makes a single a-envelope valid at every w.
# --- CLOSED FORM (2026-09-15): the three sups are now PROVEN, not calibrated ---
# The continuum row-sum has an exact closed form.  With a = s_n and the midpoint
# nodes s_m=(m+1/2)/w, 1/(2m+1) = 1/(2 w s_m), so c_n = (1/w) sum_m g(s_m) with
# g(b) = R(a,b)/(2b) on cells of width 1/w -- a midpoint rule for
#        C(a) := int_0^inf R(a,b)/(2b) db
#              = [ln(1+a^2)/2 + a arctan a] / (1+a^2).
# Verified against mp.quad to 40 digits at a=0.1..200.  This is EXACTLY the
# Lemma-M function Theta(s) already in p4_cii_m1.py:534 and p4_b1_closed_form.py,
# so the three envelope constants inherit closed-form sups:
#   sup_a C(a)                 = Theta* = 0.63492879178810795  at a*=1.500904283
#                                (root of 2a + (1-a^2)arctan a - a ln(1+a^2) = 0)
#   sup_a C(a)(1+a^2)/a^2      = 3/2  EXACTLY  (a->0 limit; decreasing in a)
#   sup_a C(a)*a               = 1.68503693513504 at a=10.0752247 (interior max;
#                                the a->inf limit is only pi/2 = 1.5708)
# Each constant below therefore = a PROVEN sup + an explicit pad.  The pad's job
# is to absorb the discrete-minus-continuum (midpoint-rule) excess, measured as
# +0.0197/w^2 at the a~1.5 peak, stable over w=3..200.  Pad 2.571e-3 covers that
# for w >= 2.77; c_shape/C_K1/C_K2/C_SUP are RETIRED from the shipped chain
# (h_hi_dyadic uses Prop. C envelope c_env_blocksup, proven, zero calibrated
# constants); kept here as legacy audit constants for gates [1]/[1b].
# HONEST BOUNDARY: the O(1/w^2) midpoint excess is measured, not proven, so the
# label is proven-sup + calibrated-pad.  It is NOT valid for small w -- at w=0.3
# max_n c_n = 0.979 and at w=1.0 it is 0.317/0.304 = 1.04x over c_shape.  Those w
# are dispatched to strong-Weyl / node_cert / head-K, which never call c_shape.
C_K1 = 1.5190      # >= 3/2 (proven sup) + 1.900e-02 pad
C_SUP = 0.6375     # >= Theta* = 0.6349287917881 (proven sup) + 2.571e-03 pad
C_K2 = 1.6870      # >= 1.6850369351350 (proven sup) + 1.963e-03 pad


# Location of the interior maximum of c_shape (where C_K1 a^2/(1+a^2) meets
# C_SUP): a_peak = sqrt(C_SUP/(C_K1-C_SUP)).  Used for block-wise sup.
A_PEAK = float(np.sqrt(C_SUP / (C_K1 - C_SUP)))




def _coth(x):
    """coth(x) for x>0, overflow-safe."""
    if x > 350.0:
        return 1.0
    return 1.0 / np.tanh(x)

def C_cont(a):
    """EXACT continuum row-sum  C(a) = int_0^inf R(a,b)/(2b) db.

    Closed form  [ln(1+a^2)/2 + a arctan a]/(1+a^2).  Identical to the Lemma-M
    Theta(s) of p4_cii_m1.py:534, whose max is the certified Theta*.  c_shape
    must dominate THIS function for the Prop A-prime/B Schur test to be sound.
    """
    a = np.asarray(a, dtype=float)
    return (np.log1p(a * a) / 2.0 + a * np.arctan(a)) / (1.0 + a * a)


def c_shape(a):
    """Certified pointwise upper bound on c_n as a function of a = (n+1/2)/w.

    These are three SCALAR bounds on the same scalar quantity c(a), so taking
    their pointwise minimum is legitimate (unlike minimising two competing
    matrix dominations, which is not -- see schur_c).
    """
    a = np.asarray(a, dtype=float)
    return np.minimum(np.minimum(C_K1 * a * a / (1.0 + a * a), C_SUP),
                      C_K2 / np.maximum(a, 1e-300))


# ==========================================================================
# (M1) head data.  D_n is EXACT (finite sum), tau >= 0, calR >= 0 entrywise.
# ==========================================================================
def pieces_M1(N, w):
    n = np.arange(N)
    wn = 1.0 / (2 * n + 1.0)
    p = np.sqrt(wn)
    s = (n + 0.5) / w
    G = 1.0 / (1.0 + s**2)
    tau = p * G
    S, Sg = s[:, None], s[None, :]
    R = 2 * S**2 * Sg**2 * (5 + S**2 + Sg**2)
    R /= (1 + (S - Sg)**2) * (1 + (S + Sg)**2) * (1 + S**2) * (1 + Sg**2)
    calR = np.outer(p, p) * R
    k = np.arange(1, N)
    Gk = 1.0 / (1.0 + (k / w)**2)
    lG = np.concatenate(([0.0], np.cumsum(Gk))) * (2.0 / (2 * n + 1.0))
    return tau, calR, wn + lG


def bisect(f, a, b, it=200):
    """Root of decreasing-in-x  f; returns midpoint of final bracket."""
    fa = f(a)
    for _ in range(it):
        m = 0.5 * (a + b)
        if (f(m) > 0) == (fa > 0):
            a = m
        else:
            b = m
    return 0.5 * (a + b)


# ==========================================================================
# Prop. L -- certified LOWER bound on h  (head-only Neumann, sound)
# ==========================================================================
def h_lo_headK(w, N, K):
    """root of sum_{j<=K} S_j^P = 1.  Each S_j^P >= 0 and F >= sum S_j^P."""
    tau, calR, D = pieces_M1(N, w)

    def g(h):
        Ai = 1.0 / (h + D)
        s = 2.0 * float(np.dot(tau**2, Ai))
        v = Ai * tau
        for _ in range(K):
            v = Ai * (calR @ v)
            s += 2.0 * float(np.dot(tau, v))
        return s - 1.0

    hi = 2.0 * float(np.dot(tau, tau)) + 1.0
    return bisect(g, THETA + 1e-9, hi) - FL


# ==========================================================================
# Prop. U -- certified UPPER bound on h  (Schur + geometric head remainder)
# ==========================================================================
def _kappa_PQ(N, w, Mq):
    """kappa >= ||calR_PQ||_2 via Frobenius on the P x (N..N+Mq) block."""
    tau, _, _ = pieces_M1(N + Mq, w)
    n = np.arange(N + Mq)
    p = np.sqrt(1.0 / (2 * n + 1.0))
    s = (n + 0.5) / w
    a = s[:N][:, None]
    b = s[N:][None, :]
    R = 2 * a**2 * b**2 * (5 + a**2 + b**2)
    R /= (1 + (a - b)**2) * (1 + (a + b)**2) * (1 + a**2) * (1 + b**2)
    blk = np.outer(p[:N], p[N:]) * R
    fro = float(np.linalg.norm(blk))
    # beyond N+Mq: R <= 5/(1+(b-a)^2) <= 5, tau-weights give a crude but valid pad
    pad = 5.0 * np.sqrt(float(np.sum(1.0 / (2 * np.arange(N) + 1.0)))
                        * w**4 / (8.0 * (N + Mq - 0.5)**4))
    return fro + pad


def h_hi_headK(w, N, K, Mq=4000):
    """Certified upper bound on h(w) (Prop. U).  All sums are finite-N head sums."""
    tau, calR, D = pieces_M1(N, w)
    tau_P2 = float(np.dot(tau, tau))
    # tail mass B_Q = sum_{n>=N} tau_n^2 <= sum_{n>=N} (w/(n+1/2))^4/(2n+1)
    B_Q = w**4 / (8.0 * (N - 0.5)**4)
    kap = _kappa_PQ(N, w, Mq)

    def H_up(x):
        """upper bound on tau_P^T (x + D_P - calR_PP)^{-1} tau_P."""
        Ai = 1.0 / (x + D)
        S0 = float(np.dot(tau**2, Ai))
        s = S0
        v = Ai * tau
        for _ in range(K):
            v = Ai * (calR @ v)
            s += float(np.dot(tau, v))
        r = THETA / x
        return s + S0 * r**(K + 1) / (1.0 - r)

    def g(h):
        eps = kap**2 / (h - THETA)
        x = h - eps
        if x <= THETA + 1e-12:
            return 1.0
        val = 2.0 * (np.sqrt(H_up(x)) + np.sqrt(B_Q / (h - THETA)))**2
        return val - 1.0

    hi = 2.0 * tau_P2 + 2.0 * B_Q + THETA + 2.0
    return bisect(g, THETA + 1e-6, hi) + FL


# ==========================================================================
# Prop. A -- PROVEN analytic (Theta*-fill), valid at EVERY w (no matrix)
# ==========================================================================
def h_hi_analytic(w, M=400000):
    """PROVEN upper bound (Theta*-fill, no calibrated envelope).  Root of U(h)=1.

    U(h) = sum_n 2 tau_n^2/(h + D_n^lo - Theta*)  with
      D_n^lo = 1/(2n+1) + (2/(2n+1)) w [atan((n+1)/w) - atan(1/w)]   <= D_n
              (G(k/w) is decreasing in k, so sum_{k=1}^{n} G >= integral_1^{n+1} G)
      tau_n^2 = G((n+1/2)/w)^2/(2n+1) <= (w/(n+1/2))^4/(2n+1)
      => sum_{n>=M} tau_n^2 <= w^4/(8 (M-1/2)^4).
    """
    n = np.arange(M, dtype=float)
    s = (n + 0.5) / w
    G = 1.0 / (1.0 + s**2)
    tau2 = G**2 / (2 * n + 1.0)
    Dlo = 1.0 / (2 * n + 1.0) + (2.0 / (2 * n + 1.0)) * w * (
        np.arctan((n + 1.0) / w) - np.arctan(1.0 / w))
    Btail = w**4 / (8.0 * (M - 0.5)**4)

    def g(h):
        d = h + Dlo - THETA
        return 2.0 * float(np.sum(tau2 / d)) + 2.0 * Btail / (h - THETA) - 1.0

    hi = 2.0 * float(np.sum(tau2)) + 2.0 * Btail + THETA + 1.0
    return bisect(g, THETA + 1e-9, hi) + FL


def h_lo_analytic(w):
    """PROVEN lower bound on h(w): Cor (iii) of P4_main_theorems.md."""
    return np.log(w) + C0 - 1.0 / (12.0 * w) - 3.0 / (16.0 * w * w) - FL


# ==========================================================================
# Prop. A' -- Schur-weighted DIAGONAL domination (PROVEN, ~60x tighter than A)
# ==========================================================================
# Schur test (phi-weighted Cauchy-Schwarz).  Let calR >= 0 entrywise, symmetric,
# phi > 0, and c_n := (calR phi)_n / phi_n.  For any v,
#   <v,calR v> <= sum_{n,m} calR_nm |v_n||v_m|
#              <= (1/2) sum_{n,m} calR_nm (phi_m/phi_n v_n^2 + phi_n/phi_m v_m^2)
#               = sum_n c_n v_n^2,
# i.e.  calR <= diag(c).   With phi_m = (m+1/2)^{-1/2} the weight COLLAPSES:
#   calR_nm phi_m/phi_n = [1/sqrt((2n+1)(2m+1))] sqrt((2n+1)/(2m+1)) R = R/(2m+1),
# so  c_n = sum_m R(s_n, s_m)/(2m+1)  -- independent of n's own weight.
# (verified against the explicit matrix product to 1e-16).
#
# Why this is 60x better than calR <= Theta* I:  c_n is TINY for small n
# (c_0 ~ 1e-2 at w=10, ~5e-4 at w=100) exactly where tau_n is LARGE, while
# max_n c_n ~ Theta* is attained near n ~ w where tau is already small.
def schur_c(w, Nc, M):
    """Upper bound on c_n for n < Nc (plus the uniform Theta* for n >= Nc).

    Sum over m < M exactly; tail m >= M uses R(a,b) <= CT * a^2/((1+a^2) b^2)
    for b >= 2a (verified numerically with margin; CT=6 is generous).
    """
    n = np.arange(M, dtype=float)
    s = (n + 0.5) / w
    inv = 1.0 / (2 * n + 1.0)
    Nc = min(Nc, M)
    c = np.empty(Nc)
    CT = 6.0                      # verified: max_{b>=2a} R (1+a^2) b^2/a^2 = 4.445
    sM = (M + 0.5) / w
    # tail m >= M:  R <= CT a^2/((1+a^2) b^2), b=(m+1/2)/w, so
    #   sum_{m>=M} R/(2m+1) <= CT a^2/(1+a^2) w^2 sum_{m>=M} 1/((m+1/2)^2 (2m+1))
    #                        <= CT a^2/(1+a^2) w^2 / (4 (M-1/2)^2)     [verified]
    tail_fac = CT * w**2 / (4.0 * (M - 0.5)**2)
    CHUNK = 256                   # blocked matmul: vectorized, bounded memory
    for i0 in range(0, Nc, CHUNK):
        i1 = min(i0 + CHUNK, Nc)
        a = s[i0:i1][:, None]
        b = s[None, :]
        R = 2 * a**2 * b**2 * (5 + a**2 + b**2)
        R /= (1 + (a - b)**2) * (1 + (a + b)**2) * (1 + a**2) * (1 + b**2)
        head = R @ inv
        aa = s[i0:i1]
        # Where the b>=2a hypothesis fails (sM < 2a) we cannot use the tail
        # majorant.  Fall back on C_SUP, a certified pointwise bound on the FULL
        # row sum c_n (C_SUP >= sup_n c_n over all w, calibrated in [2,400] and
        # consistent with the w-independent continuum limit).  Using Theta* here
        # would be UNSOUND: Theta* bounds ||calR||, not the row sums c_n, and the
        # measured max_n c_n exceeds Theta* at w=3,100,1000.
        tail = np.where(sM >= 2 * aa, tail_fac * aa**2 / (1 + aa**2), C_SUP)
        c[i0:i1] = head + tail
    # NO elementwise min with Theta* here.  calR <= diag(c) (Schur test) and
    # calR <= Theta* I (Lemma M) are two DIFFERENT matrix dominations; the
    # entrywise minimum of the two diagonals is NOT a valid dominator.
    # Counterexample: A=ones(2,2) <= diag(2,2) and <= diag(1.5,4), but
    # diag(1.5,2)-A has an eigenvalue -0.28.  Raising a valid c is always safe
    # (c' >= c => diag(c') >= diag(c) >= calR), so c is returned uncapped.
    return c


def h_hi_schur(w, Nc=3000, M=200000):
    """Proven-sup + calibrated-pad upper bound (Prop. A'); UNUSED in dispatch."""
    n = np.arange(M, dtype=float)
    s = (n + 0.5) / w
    G = 1.0 / (1.0 + s**2)
    tau2 = G**2 / (2 * n + 1.0)
    Dlo = 1.0 / (2 * n + 1.0) + (2.0 / (2 * n + 1.0)) * w * (
        np.arctan((n + 1.0) / w) - np.arctan(1.0 / w))
    # n >= Nc: use the certified pointwise envelope c_shape(s_n), NOT Theta*.
    # For n >= Nc (i.e. a >= Nc/w) c_shape decays like C_K2/a, so this is both
    # SOUND (Theta* does not bound c_n) and much tighter than a constant fill.
    cv = c_shape(s)
    cv[:Nc] = np.minimum(schur_c(w, Nc, M), cv[:Nc])
    Btail = w**4 / (8.0 * (M - 0.5)**4)

    def g(h):
        d = h + Dlo - cv
        if np.min(d) <= 0:
            return 1.0
        return 2.0 * float(np.sum(tau2 / d)) + 2.0 * Btail / (h - THETA) - 1.0

    return bisect(g, THETA + 1e-9, 60.0) + FL


# ==========================================================================
# Prop. B -- dyadic-block bound: PROVEN (Prop. C envelope, zero calibrated
# constants), valid on the dispatched band w>25.  Replaces the former
# proven-sup + calibrated-pad route: the per-block c-sup now comes from the
# closed-form Prop. C envelope c_env_blocksup, not from c_shape/C_SUP/C_K2.
# ==========================================================================
def h_hi_dyadic(w, ratio=1.002, nmax_mult=400.0):
    """PROVEN upper bound on h(w) via the Prop.-C row-sum envelope.

    Valid for ALL w > 0.75 analytically (the Prop. C envelope dominates the
    row sum everywhere); dispatched on w > 25.  Returns None for w <= 0.75
    (bisection cap would be hit vacuously).

    Partition n >= 1 into blocks [a,b).  On each block:
        sum_{n in [a,b)} tau_n^2 <= G(s_a)^2 * sum_{n in [a,b)} 1/(2n+1)
                                   = G(s_a)^2 * (1/2)[psi(b+1/2) - psi(a+1/2)]
    (G decreasing, 1/(2n+1) exact via digamma),  and
        D_n >= D_blk := (1/(2b-1)) + (2/(2b-1)) w [atan((a+1)/w) - atan(1/w)]
    for all n in [a,b) (D_n increasing, use the min).  c_n <= c_blk := RIGOROUS
    closed-form sup over [a,b) of the Prop. C envelope (each factor bounded at
    its own monotone extreme -- see c_env_blocksup).  Beyond n >= Nmax,
    algebraic tau and c tail bounds.  F(h) <= U(h) := sum_blk
    2 S_blk/(h + D_blk - c_blk) + tail, and root(U=1) >= h(w).
    Cost: O(log w) blocks with ratio ~ 1.0002.
    """
    if w <= 0.75:
        return None  # bisection cap 200.0 would be hit vacuously; dispatched only w>25
    from scipy.special import psi  # digamma function
    Nmax = max(64.0, nmax_mult * w)
    blocks = []
    a = 1.0
    at1 = np.arctan(1.0 / w)
    while a < Nmax:
        b = min(np.ceil(a * ratio) + 1.0, Nmax)
        Ga = 1.0 / (1.0 + ((a + 0.5) / w)**2)
        # EXACT block sum via digamma: sum_{n=a}^{b-1} 1/(2n+1) = (1/2)[psi(b+1/2)-psi(a+1/2)]
        S = Ga**2 * 0.5 * (psi(b + 0.5) - psi(a + 0.5))
        # D lower bound: use n=b-1 (largest n in block, smallest D_n)
        Dblk = 1.0 / (2.0 * b - 1.0) + (2.0 / (2.0 * b - 1.0)) * w * (
            np.arctan((a + 1.0) / w) - at1)
        # c upper bound: RIGOROUS closed-form sup of the Prop.-C proven envelope
        # over a in [(a+.5)/w, (b-.5)/w] -- no calibrated constant enters.
        cblk = c_env_blocksup((a + 0.5) / w, (b - 0.5) / w, w)
        blocks.append((S, max(Dblk, 0.0), cblk))
        a = b
    # n=0 block (always separate: D_0=1 exactly, tau_0^2 <= G(0.5/w)^2)
    S0 = 1.0 / (1.0 + (0.5 / w)**2)**2
    blocks.append((S0, 1.0, c_env_blocksup(0.5 / w, 0.5 / w, w)))
    Btail = w**4 / (8.0 * (Nmax - 0.5)**4)
    Sarr = np.array([x[0] for x in blocks])
    Darr = np.array([x[1] for x in blocks])
    Carr = np.array([x[2] for x in blocks])
    # tail n >= Nmax: envelope is termwise-decreasing there (a >= 400), so the
    # rigorous closed-form sup sits at the left endpoint a = (Nmax+1/2)/w.
    c_tail = c_env_blocksup((Nmax + 0.5) / w, np.inf, w)

    def g(h):
        d = h + Darr - Carr
        if np.min(d) <= 0.0 or h - c_tail <= 0.0:
            return 1.0
        return (2.0 * float(np.sum(Sarr / d))
                + 2.0 * Btail / (h - c_tail) - 1.0)

    # Bisection ENDPOINTS, not bounds: g is a decreasing majorant of F-1, so its
    # root is the returned upper bound wherever it is bracketed.  The left
    # endpoint is the certified THETA (Lemma M), not the retired calibrated
    # C_SUP -- verified numerically identical (0.000e+00 at w=26..1000), swapped
    # so that no calibrated symbol appears anywhere in the shipped deep route.
    return bisect(g, THETA + 1e-9, 200.0) + FL


# ==========================================================================
# strong band -- Weyl envelope (PROVEN), taken from deliverable (3)
# ==========================================================================
def strong_bracket(w):
    """[h_lo, h_hi] for small w (strong coupling): h = g(2) w^2 (1 + O(w^2))."""
    wm = mp.mpf(float(w))
    lo = G2_LO * wm**2 * (1 - 4 * Z4 * wm**2 / G2_LO)
    hi = G2_HI * wm**2
    return float(lo) - FL, float(hi) + FL


def strong_invert(u):
    """w* with h(w*)=u, u small.  h increasing in w => invert both envelopes."""
    if u >= 0.30:
        return None
    w_lo = bisect(lambda w: u - strong_bracket(w)[1], 1e-9, 0.60)
    w_hi = bisect(lambda w: u - strong_bracket(w)[0], 1e-9, 0.60)
    if w_hi < w_lo:
        w_lo, w_hi = w_hi, w_lo
    return w_lo, w_hi


# ==========================================================================
# INVERSION  h-envelope -> T_c/omega_E bracket
# ==========================================================================
def tc_from_w(w_lo, w_hi):
    """T_c/omega_E = 1/(2 pi w*).  Larger w* => SMALLER T_c."""
    return 1.0 / (TWO_PI * w_hi), 1.0 / (TWO_PI * w_lo)


def invert_monotone(u, hlo_f, hhi_f, wa, wb):
    """Enclose w* solving h(w*)=u, given hlo_f <= h <= hhi_f both increasing.

    hhi(w) <= u  ==>  h(w) <= u  ==>  w <= w*   (h increasing)  => w_lo
    hlo(w) >= u  ==>  h(w) >= u  ==>  w >= w*                   => w_hi
    """
    if hhi_f(wb) < u or hlo_f(wa) > u:
        return None
    w_lo = bisect(lambda w: hhi_f(w) - u, wb, wa)
    w_hi = bisect(lambda w: hlo_f(w) - u, wb, wa)
    if w_hi < w_lo:
        w_lo, w_hi = w_hi, w_lo
    return w_lo, w_hi


def invert_grid(ws, hlo, hhi, u):
    """Cell-wise inversion on a monotone tabulated envelope (deliverable-3 style)."""
    w_lo = w_hi = None
    for i in range(len(ws)):
        if hhi[i] <= u:
            w_lo = ws[i]
    for i in range(len(ws) - 1, -1, -1):
        if hlo[i] >= u:
            w_hi = ws[i]
    if w_lo is None or w_hi is None or w_hi < w_lo:
        return None
    return float(w_lo), float(w_hi)




def _psi_complex(x, y, shift=24):
    """psi(x+iy) for x>0 (array x, scalar y), via recurrence + Stirling series.

    psi(z) = psi(z+k) - sum_{j=0}^{k-1} 1/(z+j); apply until Re >= shift, then
    the standard asymptotic psi(z) ~ ln z - 1/(2z) - sum B_2n/(2n z^2n).
    Accurate to ~1e-14 for shift=24, which is far tighter than anything the
    bracket needs (the PROVEN bound c_proven_upper does not rely on this
    routine's error; it is only used for reporting/gate comparison of c_exact).
    """
    x = np.asarray(x, dtype=float)
    z = x + 1j*y
    acc = np.zeros_like(z)
    for _ in range(4000):
        m = z.real < shift
        if not np.any(m):
            break
        acc[m] += 1.0/z[m]
        z = np.where(m, z+1.0, z)
    zi = 1.0/z
    z2 = zi*zi
    ser = z2*(1.0/12.0 + z2*(-1.0/120.0 + z2*(1.0/252.0
          + z2*(-1.0/240.0 + z2*(1.0/132.0 + z2*(-691.0/32760.0 + z2*(1.0/12.0)))))))
    return np.log(z) - 0.5*zi - ser - acc


def _re_psi_half(w):
    """Re psi(1/2 + i w)."""
    return float(np.real(_psi_complex(np.array([0.5]), w)[0]))


def c_exact(n, w):
    """EXACT Schur row sum c_n via (C1).  Identity, no calibration."""
    n = np.asarray(n, dtype=float)
    a = (n+0.5)/w
    x = n+1.0
    P = _psi_complex(x, w)
    num = (-_re_psi_half(w) + P.real - a*P.imag
           + a*(np.pi/2.0)*_coth(np.pi*w))
    return num/(1.0+a*a)


def c_proven_upper(n, w):
    """PROVEN upper bound on c_n: (C1) with B1 on Re, B3 on Im, B2 on psi(1/2).

    No calibrated constant enters.  Verified ratio >= 1 at every tested
    (n, w); the bound is tight (ratio -> 1 as n grows, worst 1.000003 at
    n=5000).  Small n is loose (factor ~4) but c_n is ~1e-5 there, so the
    looseness is far inside the pad it replaces.
    """
    n = np.asarray(n, dtype=float)
    a = (n+0.5)/w
    x = n+1.0
    q = x*x + w*w
    ReP_up = 0.5*np.log(q)
    ImP_lo = np.arctan(w/x) - w/(2.0*q)
    qh = 0.25 + w*w
    Repsi_half_lo = 0.5*np.log(qh) - 0.5/qh
    num = (-Repsi_half_lo + ReP_up - a*ImP_lo
           + a*(np.pi/2.0)*_coth(np.pi*w))
    return num/(1.0+a*a)


def c_env_proven(a, w):
    """PROVEN upper envelope on the row sum, as a function of a=(n+1/2)/w.

    From the identity (C1) with the three standard digamma inequalities:
      Re psi(n+1+iw) <= (1/2)ln q                         (B1)
      Re psi(1/2+iw) >= (1/2)ln q_h - 1/(2 q_h)           (B2)
      -a Im psi(n+1+iw) <= -a atan(w/x) + a w/(2q)              (B3)
    with x = a w + 1/2, q = x^2 + w^2, q_h = 1/4 + w^2.  Every term is closed
    form in (a, w): NO calibrated constant, no measured sup.  Unlike c_shape
    this is valid for all a>0 and all w>0, and it is what Prop. B now uses.

    The -a*atan(w/x) term is RETAINED (not discarded for a cheaper bound): it
    costs nothing analytically and makes the envelope essentially exact -- the
    measured min ratio env/c_exact is 1.000000 at large w, and sup_a env
    reproduces Theta* = 0.63492879 from the closed form instead of assuming it.
    """
    a = np.asarray(a, dtype=float)
    x = a*w + 0.5
    q = x*x + w*w
    qh = 0.25 + w*w
    T1 = 0.5*np.log(q/qh) + 0.5/qh
    T2 = a*((np.pi/2.0)*_coth(np.pi*w) - np.arctan(w/x) + w/(2.0*q))
    return (T1+T2)/(1.0+a*a)


def c_env_blocksup(a_lo, a_hi, w):
    """RIGOROUS closed-form sup of the Prop.-C envelope over a in [a_lo,a_hi].

    E(a,w) = [ T1(a) + a*B(a) ] / (1+a^2),  T1 = (1/2)ln(q/q_h) + 1/(2 q_h),
    B(a) = kappa(w) - atan(w/x) + w/(2q),  x = a w + 1/2,  q = x^2 + w^2,
    q_h = 1/4 + w^2,  kappa(w) = (pi/2)coth(pi w).

    No grid and no pad: each factor is bounded by its own monotone extreme.
      * T1 is INCREASING in a (q increasing)              -> evaluate at a_hi
      * -atan(w/x) is INCREASING in a                     -> evaluate at a_hi
      * +w/(2q)   is DECREASING in a                      -> evaluate at a_lo
      * 1/(1+a^2) is DECREASING in a                      -> evaluate at a_lo
      * a*B: B has a sign, so use a_hi if B>=0 else a_lo (never mixes signs)
    Each replacement only increases the value, so the result is a true upper
    bound for EVERY a in the block -- valid without unimodality.  This replaces
    an earlier 96-point grid sup inflated by a Lipschitz-style slack whose
    justification did not survive checking (the a*kappa term breaks the
    termwise |a dE/da| <= 2 max E estimate the slack relied on).  A measured
    grid sup plus a guessed pad is calibration, not proof; this is closed form.

    Cost on the dispatched band w>25: ratio <= 1.041 vs the true block max
    (tail block 1.002), so the bracket does not measurably loosen.
    """
    kap = 0.5*np.pi*_coth(np.pi*w)
    qh = 0.25 + w*w
    if not np.isfinite(a_hi):
        if a_lo < 2.0:
            raise ValueError("tail block needs a_lo>=2, got %r" % (a_lo,))
        Cw = 0.5*np.log(((w+0.5)**2 + w*w)/qh) + 0.5/qh
        return float((np.log(a_lo)+Cw)/(1.0+a_lo*a_lo)
                     + a_lo*kap/(1.0+a_lo*a_lo)
                     + 1.0/(2.0*a_lo*w*(1.0+a_lo*a_lo)))
    if a_hi < a_lo:
        a_lo, a_hi = a_hi, a_lo
    x_hi = a_hi*w + 0.5
    q_hi = x_hi*x_hi + w*w
    x_lo = a_lo*w + 0.5
    q_lo = x_lo*x_lo + w*w
    T1 = 0.5*np.log(q_hi/qh) + 0.5/qh
    B = kap - np.arctan(w/x_hi) + w/(2.0*q_lo)
    aB = a_hi*B if B >= 0.0 else a_lo*B
    return float((T1+aB)/(1.0+a_lo*a_lo))


# ==========================================================================
# BAND DISPATCH -- one certified [h_lo, h_hi] per w, with a method label
# ==========================================================================
W_STRONG = 0.14        # below this the Weyl envelope is the tightest proven route
W_CROSS = 3.0          # node_cert (Rayleigh+Kato-Temple) band; head-K needs h>Theta*
W_HEAD = 25.0          # above this the O(N^2) head route costs >5s/node
N_HEAD_MAX = 2400      # head matrix cap (cost ~ N^2 per bisection step)


def h_bracket(w, K=8, N=None, mode="auto"):
    if mode is None:
        mode = "auto"
    """Certified [h_lo, h_hi, method, rigor] for h(w) at a single w.

    mode="auto" picks the tightest AFFORDABLE route:
      w <= 0.6            -> strong Weyl envelope        (PROVEN analytic)
      0.6 < w <= W_HEAD   -> Prop. L / Prop. U head-K    (CERTIFICATE, finite-N)
      w > W_HEAD          -> Cor (iii) / Prop. B dyadic  (PROVEN analytic)
    """
    if mode == "strong" or (mode == "auto" and w <= W_STRONG):
        lo, hi = strong_bracket(w)
        return lo, hi, "strong-Weyl", "proven"
    if mode == "dyadic" or (mode == "auto" and w > W_HEAD):
        hi = h_hi_dyadic(w)
        if hi is None:
            # w <= 0.75: bisection cap unreachable; fall back to proven Weyl envelope
            lo2, hi2 = strong_bracket(w)
            return lo2, hi2, "strong-Weyl", "proven"
        return (h_lo_analytic(w), hi, "deep-dyadic", "proven")
    if mode == "cross" or (mode == "auto" and w <= W_CROSS):
        # Crossover band.  The head-K route needs h(w) > Theta* (the Neumann
        # series and the Schur-complement step both divide by h - Theta*), and
        # h(w) < Theta* for w <~ 0.93.  Deliverable (3)'s node_cert (eigenvector
        # Rayleigh + Kato-Temple) has no such restriction and is ~0.1 s/node.
        lo, hi = node_cert(w, Nv=140, Na=70)[:2]
        if hi is None:
            lo2, hi2 = strong_bracket(w)
            return lo2, hi2, "strong-Weyl", "proven"
        return float(lo), float(hi), "crossover-KT", "certificate"
    if N is None:
        N = int(min(N_HEAD_MAX, max(300, 60 * w)))
    return (h_lo_headK(w, N, K), h_hi_headK(w, N, K), "head-K%d" % K, "certificate")


# ==========================================================================
# MAIN GATE
# ==========================================================================
def main(K=8, quick=False):
    t0 = time.time()
    print("=" * 78)
    print("P4 GLOBAL T_c ALL-COUPLING -- unified two-sided bracket, ALL lambda>0")
    print("=" * 78)

    print("\n[0] CERTIFIED-CONSTANT HYGIENE (g(2) CERT-D nested in CERT-B)")
    # Guard against a transcription slip in the tighter CERT-D digits: CERT-D
    # must nest strictly inside the independently-derived CERT-B bracket.
    ck(G2_LO_B < G2_LO and G2_HI < G2_HI_B,
       "CERT-D nests strictly inside CERT-B")
    ck(G2_LO < G2_HI and float(G2_HI - G2_LO) < 1e-8,
       "CERT-D bracket ordered, width < 1e-8")
    print("    CERT-D width %.3e  nested in CERT-B width %.3e"
          % (float(G2_HI - G2_LO), float(G2_HI_B - G2_LO_B)))

    print("\n[1a] CLOSED-FORM SUPS of the envelope constants (reviewer T2)")
    print("     C(a)=[ln(1+a^2)/2+a atan a]/(1+a^2) == Lemma-M Theta(a); sups PROVEN")
    for aq in (0.1, 1.0, 1.5, 5.0, 20.0):
        am = mp.mpf(aq)
        Rk = lambda b: (2 * am**2 * b**2 * (5 + am**2 + b**2) /
                        ((1 + (am - b)**2) * (1 + (am + b)**2) *
                         (1 + am**2) * (1 + b**2)))
        q = mp.quad(lambda b: Rk(b) / (2 * b), [0, am, 2 * am, mp.inf])
        ck(abs(float(q) - float(C_cont(aq))) < 1e-13,
           "closed form = quad(R/2b) at a=%s" % aq)
    THETA_STAR = 0.6349287917881079
    A_STAR = 1.500904283161906
    ck(abs(float(C_cont(A_STAR)) - THETA_STAR) < 1e-13,
       "C(a*) = certified Lemma-M Theta* at a*=1.5009042832")
    Nf = lambda a: 2 * a + (1 - a * a) * np.arctan(a) - a * np.log1p(a * a)
    ck(abs(Nf(A_STAR)) < 1e-12, "a* is a root of the exact stationarity numerator")
    ck(Nf(1.0) > 0 and Nf(2.0) < 0, "N changes sign around a* (interior max)")
    ag = np.geomspace(1e-6, 1e6, 2000001)
    Cg = C_cont(ag)
    s_C = float(Cg.max())
    s_K1 = float((Cg * (1 + ag**2) / ag**2).max())
    s_K2 = float((Cg * ag).max())
    ck(s_C <= THETA_STAR + 1e-9, "sup_a C(a) = Theta* (%.13f)" % s_C)
    ck(abs(s_K1 - 1.5) < 1e-6, "sup_a C(a)(1+a^2)/a^2 = 3/2 exactly (%.9f)" % s_K1)
    ck(abs(s_K2 - 1.6850369351350385) < 1e-6,
       "sup_a C(a)*a = 1.68503693514 interior (%.11f)" % s_K2)
    ck(s_K2 > float(mp.pi / 2), "the a*C sup is INTERIOR, not the a->oo limit pi/2")
    for nm, cal, sup in (("C_SUP", C_SUP, s_C), ("C_K1", C_K1, 1.5),
                         ("C_K2", C_K2, s_K2)):
        ck(cal > sup, "%s=%.4f exceeds its proven sup %.13f (pad %.3e)"
           % (nm, cal, sup, cal - sup))
        print("     %-6s = %.4f  >=  proven sup %.13f   pad %.3e"
              % (nm, cal, sup, cal - sup))
    rr = float((Cg / c_shape(ag)).max())
    ck(rr <= 1.0, "c_shape dominates C(a) for ALL a>0 (max ratio %.8f)" % rr)
    print("     max_a C(a)/c_shape(a) = %.8f  (ALL a>0, closed form)" % rr)
    bad = np.minimum(np.minimum(C_K1 * ag**2 / (1 + ag**2), THETA_STAR - 1e-6),
                     C_K2 / ag)
    ck(float((Cg / bad).max()) > 1.0,
       "teeth: C_SUP shaved under Theta* stops dominating C (control)")

    print("\n[1b] LEGACY audit: discrete-vs-continuum excess (RETIRED calibrated-pad scheme)")
    print("     Shipped route uses Prop. C proven envelope (pad-free); this audits the old c_shape.")
    K_EXC = 0.0197
    w_ok = (K_EXC / (C_SUP - THETA_STAR))**0.5
    ck(w_ok < W_HEAD, "pad-coverage threshold w=%.3f is below the live band %.0f"
       % (w_ok, W_HEAD))
    print("     excess ~ %.4f/w^2 at the peak; pad %.3e covers it for w >= %.3f"
          % (K_EXC, C_SUP - THETA_STAR, w_ok))
    print("     c_shape/C_SUP/C_K2 RETIRED from shipped chain (Prop. C envelope is proven,"
          " pad-free); kept as LEGACY audit below.")
    ck(K_EXC / W_HEAD**2 < 0.05 * (C_SUP - THETA_STAR),
       "at w=%.0f the excess is >20x inside the pad" % W_HEAD)
    for wsm in (0.3, 1.0):
        Nc = max(80, int(4 * wsm))
        cd = schur_c(wsm, Nc, 200000)
        asm = (np.arange(Nc) + 0.5) / wsm
        rsm = float((cd / c_shape(asm)).max())
        ck(rsm > 1.0, "teeth: c_shape FAILS at w=%s (ratio %.4f) -- band excluded"
           % (wsm, rsm))
        print("     w=%4.1f: max c_n/c_shape = %.4f  EXCEEDS 1 (not dispatched here)"
              % (wsm, rsm))

    print("\n[1] LEGACY audit: Schur envelope c_shape soundness (RETIRED from shipped chain)")
    print("    Shipped route uses Prop. C proven envelope (proven, pad-free).")
    print("    This checks the OLD c_shape envelope on former dispatch band w>25:")
    for w in (26.0, 50.0, 100.0, 200.0):
        Nc = max(80, int(4 * w))
        c = schur_c(w, Nc, int(max(60 * w, 60000)))
        a = (np.arange(Nc) + 0.5) / w
        r = float((c / c_shape(a)).max())
        ck(r <= 1.0 + 1e-9, "c_n <= c_shape at w=%s (ratio %.6f)" % (w, r))
        print("    w=%7.1f: max c_n/c_shape(s_n) = %.6f" % (w, r))
    print("    NEGATIVE CONTROL (w < 2.8 boundary, c_shape FAILS; NOT dispatched):")
    for w in (0.7, 1.0, 1.8):
        Nc = max(80, int(4 * w))
        c = schur_c(w, Nc, 200000)
        a = (np.arange(Nc) + 0.5) / w
        r = float((c / c_shape(a)).max())
        ck(r > 1.0, "teeth: c_shape FAILS at w=%s (%.4f>1, band excluded)" % (w, r))
        print("    w=%7.1f: max c_n/c_shape(s_n) = %.6f  EXCEEDS (ok, not used)" % (w, r))

    print("\n[2] ENCLOSURE vs an INDEPENDENT builder (dense H_mat eigenvalue)")
    print("    NOTE: H_mat is built from the KAY entry formulas, NOT from the")
    print("    (M1) tau/calR/D decomposition -- a genuinely independent route.")
    print("    %7s %11s %11s %11s %12s %5s" % ("w", "h_lo", "h_true", "h_hi", "method", "ok"))
    for w in (0.3, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 26.0, 50.0, 200.0):
        lo, hi, meth, rig = h_bracket(w, K=K)
        ht = float(np.linalg.eigvalsh(H_mat(2500, w))[-1])
        ok = (lo <= ht <= hi)
        ck(ok, "bracket encloses h_true at w=%s" % w)
        print("    %7.2f %11.6f %11.6f %11.6f %12s %5s" % (w, lo, ht, hi, meth, ok))

    print("\n[3] INVERSION DIRECTION (larger h_hi => smaller w_lo => larger F_hi)")
    u = 1.0 / 0.3
    wlo_a = bisect(lambda x: h_hi_dyadic(x) - u, 1e5, 1.0)
    wlo_b = bisect(lambda x: 1.02 * h_hi_dyadic(x) - u, 1e5, 1.0)
    F_a = 1.0 / (TWO_PI * wlo_a)
    F_b = 1.0 / (TWO_PI * wlo_b)
    ck(F_b >= F_a, "inflating h_hi does not shrink the F_hi bound")
    print("    lambda=0.3: F_hi(base)=%.6e  F_hi(h_hi x1.02)=%.6e  monotone=%s"
          % (F_a, F_b, F_b >= F_a))

    print("\n[4] SEAM continuity: adjacent bands must OVERLAP (both contain h)")
    # Each pair is evaluated at a w where BOTH routes are valid.  (The head-K
    # route is NOT valid for w <~ 0.93, where h(w) < Theta*, so the strong/cross
    # seam is tested against node_cert, not against head-K.)
    for w, m1, m2 in ((0.14, "strong", "cross"), (0.30, "strong", "cross"),
                      (3.0, "cross", "head"), (5.0, "cross", "head"),
                      (24.0, "head", "dyadic"), (30.0, "head", "dyadic")):
        la, ha, ma, _ = h_bracket(w, K=K, mode=m1)
        lb, hb, mb, _ = h_bracket(w, K=K, mode=m2)
        ov = (la <= hb) and (lb <= ha)
        ck(ov, "seam overlap at w=%s" % w)
        print("    w=%6.2f: %12s=[%.6f,%.6f]  %12s=[%.6f,%.6f]  overlap=%s"
              % (w, ma, la, ha, mb, lb, hb, ov))

    print("\n[5] UNIFIED CERTIFIED T_c(lambda) TABLE  (F_lo <= T_c/w_E <= F_hi)")
    lam_grid = [0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.582,
                0.7, 1.0, 2.0, 5.0, 10.0, 20.0, 34.0, 40.0, 50.0]
    if quick:
        lam_grid = [0.05, 0.1, 0.3, 0.582, 1.0, 10.0, 50.0]
    # Precompute ONE monotone certified envelope on a log-w grid, then invert by
    # table lookup (deliverable-3 style).  Nested bisection over h_bracket would
    # rebuild an NxN matrix 200x200 times per lambda -- hours; this is seconds.
    # Grid density is the DOMINANT width source: the staircase inversion loses
    # ~one cell, so rel width in T_c ~ the grid step in w.  The head-K band is
    # the expensive one (O(N^2) per node), so put most nodes in the cheap
    # dyadic/crossover regions by using a piecewise-geometric grid.
    if quick:
        ws_grid = np.unique(np.concatenate([
            np.geomspace(0.10, W_CROSS, 30),      # crossover: node_cert, ~0.1s
            np.geomspace(W_CROSS, W_HEAD, 10),    # head-K: expensive
            np.geomspace(W_HEAD, 1e13, 70)]))     # dyadic: ~0.01s
    else:
        ws_grid = np.unique(np.concatenate([
            np.geomspace(0.10, W_CROSS, 120),
            np.geomspace(W_CROSS, W_HEAD, 26),
            np.geomspace(W_HEAD, 1e13, 320)]))
    NW = len(ws_grid)
    lo_grid = np.zeros(NW)
    hi_grid = np.zeros(NW)
    meth_grid = []
    print("    building certified w-grid (%d nodes, w in [%.2g, 1e13]) ..." % (NW, ws_grid[0]))
    for i, wv in enumerate(ws_grid):
        lo_grid[i], hi_grid[i], m, _rg = h_bracket(float(wv), K=K)
        meth_grid.append(m)
    # MONOTONE HARDENING (sound, and strictly tightening).  h is PROVEN
    # strictly increasing, so for i' > i:  h(w_i') > h(w_i) >= lo_grid[i],
    # i.e. a lower bound at ANY smaller node is also a lower bound here.
    # Hence the running maximum of lo_grid is still a valid lower-bound array,
    # and symmetrically the reverse running minimum of hi_grid is still a valid
    # upper-bound array.  This is not cosmetic: at the head-K -> dyadic seam
    # (w=25) the head route is tighter than dyadic just above it, so without
    # this the raw arrays are non-monotone (a legitimate method-switch artifact,
    # NOT a soundness bug -- invert_grid scans for extremes and is correct
    # either way).  Hardening both fixes the artifact and inherits the tighter
    # bound across the seam for free.
    lo_raw = lo_grid.copy()
    hi_raw = hi_grid.copy()
    lo_grid = np.maximum.accumulate(lo_grid)
    hi_grid = np.minimum.accumulate(hi_grid[::-1])[::-1]
    n_lift = int(np.sum(lo_grid > lo_raw + 0.0))
    n_drop = int(np.sum(hi_grid < hi_raw - 0.0))
    print("    monotone hardening: %d lo nodes lifted, %d hi nodes lowered"
          % (n_lift, n_drop))
    # The hardening must never invert the bracket, and must remain an enclosure.
    ck(bool(np.all(hi_grid - lo_grid > 0.0)),
       "hardened bracket still strictly positive-width at every node")
    ck(bool(np.all(lo_grid >= lo_raw - 1e-15)) and bool(np.all(hi_grid <= hi_raw + 1e-15)),
       "hardening only tightens (lo up, hi down), never loosens")
    # After clamping the arrays are non-DEcreasing: a lifted node ties its
    # predecessor by construction, so ">0" is not the right assertion (and a
    # method switch legitimately breaks strictness -- see the w=25 seam).
    ck(bool(np.all(np.diff(lo_grid) >= 0.0)), "grid h_lo non-decreasing (hardened)")
    ck(bool(np.all(np.diff(hi_grid) >= 0.0)), "grid h_hi non-decreasing (hardened)")
    # Strictness IS required wherever the method does not change: within a
    # single band the bound is a smooth function of w and must resolve h's
    # proven strict growth.  Ties/decreases there would signal a real bug.
    same = np.array([meth_grid[i] == meth_grid[i + 1] for i in range(NW - 1)])
    d_lo = np.diff(lo_raw)[same]
    d_hi = np.diff(hi_raw)[same]
    ck(bool(np.all(d_lo > 0.0)),
       "raw h_lo strictly increasing within every band (%d intra-band steps)" % len(d_lo))
    ck(bool(np.all(d_hi > 0.0)),
       "raw h_hi strictly increasing within every band (%d intra-band steps)" % len(d_hi))
    n_seam_break = int(np.sum(np.diff(lo_raw)[~same] <= 0.0))
    print("    intra-band strict: %d/%d steps; method-switch steps with a "
          "non-monotone lo: %d (clamped above)"
          % (len(d_lo), NW - 1, n_seam_break))

    print("    %8s %13s %13s %9s %13s %12s"
          % ("lambda", "F_lo", "F_hi", "rel_wid", "method", "rigor"))
    rows = []
    for lam in lam_grid:
        u = 1.0 / lam
        if u < 0.02:
            r = strong_invert(u)
            if r is None:
                print("    %8.3f  (not covered)" % lam)
                ck(False, "lambda=%s uncovered" % lam)
                continue
            wlo, whi = r
            meth, rig = "strong-Weyl", "proven"
        else:
            r = invert_grid(ws_grid, lo_grid, hi_grid, u)
            if r is None:
                print("    %8.3f  (not covered)" % lam)
                ck(False, "lambda=%s uncovered" % lam)
                continue
            wlo, whi = r
            j = int(np.argmin(np.abs(ws_grid - 0.5 * (wlo + whi))))
            meth = meth_grid[j]
            rig = ("proven" if meth in ("strong-Weyl", "deep-dyadic")
                   else "certificate")
        Flo, Fhi = tc_from_w(wlo, whi)
        rw = (Fhi - Flo) / Flo
        rows.append((lam, Flo, Fhi, rw, meth, rig))
        ck(Flo <= Fhi, "F_lo<=F_hi at lambda=%s" % lam)
        print("    %8.3f %13.6e %13.6e %9.2e %13s %12s"
              % (lam, Flo, Fhi, rw, meth, rig))

    print("\n[6] WORST-CASE RELATIVE WIDTH by band")
    by = {}
    for (lam, Flo, Fhi, rw, meth, rig) in rows:
        by.setdefault(meth, []).append((rw, lam))
    for meth in sorted(by):
        w_, l_ = max(by[meth])
        print("    %13s: worst rel width %.3e (at lambda=%s)" % (meth, w_, l_))
    if rows:
        gw, gl = max((r[3], r[0]) for r in rows)
        print("    %13s: worst rel width %.3e (at lambda=%s)" % ("GLOBAL", gw, gl))

    print("\n[7] ASYMPTOTIC LIMITS")
    pref = float(mp.e**mp.mpf(C0) / TWO_PI)
    ck(abs(pref - 0.2529996837) < 1e-8, "weak prefactor = %s" % pref)
    print("    weak  lambda->0 : (1/2pi)exp(c0) = %.10f  (target 0.2529996837)" % pref)
    # CORRECTED 2026-09-15 (v1 had this reciprocal-flipped, see .md section 5).
    # h(w) -> g(2) w^2 as w -> 0, so h(w*)=1/lambda gives w* = 1/sqrt(g(2) lambda)
    # and T_c/w_E = 1/(2 pi w*) = sqrt(g(2) lambda)/(2 pi).  The prefactor is
    # therefore sqrt(g(2))/(2 pi), NOT 1/(2 pi sqrt(g(2))): the two differ by
    # exactly the factor g(2)=1.3181, a reciprocal-sqrt slip, and for an
    # Einstein spectrum sqrt(<w^2>) = w_E so no normalisation can absorb it.
    s_lo = float(mp.sqrt(G2_LO) / TWO_PI)
    s_hi = float(mp.sqrt(G2_HI) / TWO_PI)
    ck(s_lo <= s_hi and abs(s_hi - s_lo) < 1e-8,
       "strong prefactor bracket is tight and ordered (CERT-D)")
    # NOTE the CERT-D bracket has width 4.2e-11, i.e. TIGHTER than the 10-digit
    # literal 0.1827262478, and it straddles the ...77|...78 rounding boundary
    # ([0.18272624774670, 0.18272624778904]).  So NEITHER 10-dp literal is
    # enclosed, and asserting enclosure of one would be a false check.  The two
    # sound checks are (a) nesting inside the INDEPENDENT CERT-B derivation
    # (M=12 rational arithmetic -- a different proof, not a re-run of CERT-D),
    # and (b) agreement with the published headline to the published precision.
    b_lo = float(mp.sqrt(G2_LO_B) / TWO_PI)
    b_hi = float(mp.sqrt(G2_HI_B) / TWO_PI)
    mid = 0.5 * (s_lo + s_hi)
    ck(b_lo < s_lo and s_hi < b_hi,
       "CERT-D prefactor nests strictly in the independent CERT-B prefactor")
    ck(abs(mid - 0.1827262477) < 1e-9 and abs(mid - 0.1827262478) < 1e-9,
       "prefactor agrees with KAY f(0)=0.1827262477 to 1e-9 (straddles 77|78)")
    # Claim C publishes 0.1827 to FOUR decimals, so the correct agreement test is
    # a 4-dp round-trip (half-ulp 5e-5), not a 5e-6 window -- the true gap is
    # 0.18272624-0.1827 = 2.6e-5, which is agreement at the published precision
    # and would be a spurious FAIL against a tolerance tighter than the headline.
    ck(round(mid, 4) == 0.1827 and abs(mid - 0.1827) < 5e-5,
       "prefactor rounds to the merged Claim C headline 0.1827 at 4 dp")
    # the WRONG v1 constant must be OUTSIDE the bracket -- teeth test
    ck(not (s_lo <= 0.1386242875 <= s_hi),
       "the v1 value 1/(2pi sqrt(g2))=0.13862 is correctly EXCLUDED")
    ck(abs(0.1827262478 / 0.1386242875 - float(0.5 * (G2_LO + G2_HI))) < 1e-6,
       "the v1/v2 ratio is exactly g(2) -- confirms a reciprocal-sqrt slip")
    print("    strong lambda->oo: sqrt(g2)/(2pi) in [%.14f,%.14f] (CERT-D, w %.1e)"
          % (s_lo, s_hi, s_hi - s_lo))
    print("      nested in CERT-B [%.14f,%.14f] (independent derivation)"
          % (b_lo, b_hi))
    print("      reconciles with the merged Claim C (T_c < 0.1827 sqrt(eta/M))")
    print("      and KAY Thm 1(c) f(0)=0.1827262477..; same w_E scale.")
    # end-to-end: the certified inversion itself must land on this constant
    for lam in (1.0e4, 1.0e6):
        r = strong_invert(1.0 / lam)
        tlo, thi = tc_from_w(r[0], r[1])
        rl, rh = tlo / mp.sqrt(lam), thi / mp.sqrt(lam)
        # the inverted bracket is WIDER than the prefactor bracket (it carries the
        # Weyl O(w^2) correction too), so here enclosure IS the right assertion:
        # it must contain the whole CERT-D prefactor interval, hence the limit.
        ck(float(rl) <= s_lo and s_hi <= float(rh),
           "inverted bracket encloses the CERT-D prefactor at lambda=%g" % lam)
        print("      lambda=%-9g T_c/(w_E sqrt(lambda)) in [%.9f,%.9f]"
              % (lam, float(rl), float(rh)))

    # --- INDEPENDENT confirmation of the strong constant (added 2026-09-15) ---
    # The two rows above both come from strong_invert, i.e. from the SAME Weyl
    # envelope whose constant is under test -- self-consistent, not independent.
    # Strong coupling means SMALL w*, which is cheap for the KAY-entry builder
    # H_mat, so we can solve lam*h(w*)=1 directly with no envelope anywhere and
    # Richardson-extrapolate the O(w^2) Weyl correction away.  This is the check
    # that would have caught the v1 reciprocal-sqrt slip on its own.
    def _h_true(w, Nq):
        return float(np.linalg.eigvalsh(H_mat(Nq, w))[-1])

    def _solve_w(lam, Nq):
        u, a, b = 1.0 / lam, 1e-9, 1.0
        for _ in range(60):
            m = 0.5 * (a + b)
            if _h_true(m, Nq) > u:
                b = m
            else:
                a = m
        return 0.5 * (a + b)

    pf = {}
    for lam in (1600.0, 6400.0):
        pf[lam] = 1.0 / (2.0 * np.pi * _solve_w(lam, 300) * np.sqrt(lam))
    rich = (4.0 * pf[6400.0] - pf[1600.0]) / 3.0
    print("      independent H_mat solve of lam*h(w*)=1 (no envelope):")
    print("        prefactor(1600)=%.9f  prefactor(6400)=%.9f" % (pf[1600.0], pf[6400.0]))
    print("        Richardson limit=%.9f  vs CERT-D mid %.9f" % (rich, mid))
    ck(abs(rich - mid) < 1e-8,
       "INDEPENDENT builder reproduces sqrt(g2)/2pi to 1e-8 (envelope-free)")
    ck(abs(rich - 0.1386242875) > 1e-3,
       "INDEPENDENT builder EXCLUDES the v1 constant 0.13862 (teeth test)")

    # ------------------------------------------------------------------
    print("\n[8] PROP. C -- exact closed form (C1) and the PROVEN envelope")
    print("    %6s %5s %16s %16s %11s" % ("w", "n", "c_Richardson", "c_exact_C1", "rel_err"))
    # The brute row sum truncates at M with remainder ~1/M (terms ~ w^2/2m^3);
    # the 1/M expansion behind the Richardson step c_R=2c(2M)-c(M) only holds
    # for M >> w, so M scales with w.  The residual of the combination is then
    # the float64 cancellation floor (~2e-8), flat in w -- set tolerance above.
    worst_id = 0.0
    for w in (0.5, 2.0, 7.0, 30.0, 200.0):
        M_w = int(max(100000.0, 4000.0 * w))
        cM = schur_c(w, 12, M_w)
        c2M = schur_c(w, 12, 2 * M_w)
        for n in (0, 1, 5, 11):
            ce = c_exact(n, w)
            cR = 2.0 * c2M[n] - cM[n]
            rel = abs(cR - ce) / max(abs(ce), 1e-300)
            worst_id = max(worst_id, rel)
            if n in (0, 5):
                print("    %6.4g %5d %16.12f %16.12f %11.2e"
                      % (w, n, cR, ce, rel))
    ck(worst_id < 1e-7,
       "(C1) closed form == Richardson row sum, worst rel %.2e" % worst_id)
    # teeth: the raw (un-extrapolated) sum must MISS at this tolerance, else
    # the check above would pass for reasons unrelated to the identity.
    _cm = schur_c(30.0, 12, 100000)
    ck(abs(_cm[5] - c_exact(5, 30.0)) / c_exact(5, 30.0) > 1e-9,
       "teeth: un-extrapolated M=1e5 sum is truncation-limited (>1e-9)")

    # Per-node proven bound must dominate the exact row sum everywhere.
    worst_dom = 1.0e30
    for w in (0.3, 1.0, 3.0, 12.0, 60.0, 400.0, 5000.0):
        for n in range(0, 24):
            worst_dom = min(worst_dom, c_proven_upper(n, w) - c_exact(n, w))
    ck(worst_dom >= 0.0,
       "PROVEN envelope dominates exact c_n, worst slack %.3e" % worst_dom)
    # Teeth for the domination check: shave the envelope at its TIGHTEST point
    # (small n carries big slack, so a global shave proves nothing).
    tight = min(((c_proven_upper(n, w) - c_exact(n, w)) / c_proven_upper(n, w), n, w)
                for w in (0.3, 1.0, 3.0, 12.0, 60.0, 400.0, 5000.0)
                for n in range(0, 24))
    _rel_slack, n_t, w_t = tight
    print("    tightest envelope point: n=%d w=%.4g rel slack %.3e"
          % (n_t, w_t, _rel_slack))
    shave = 1.0 - min(_rel_slack * 2.0, 0.5)
    ck(c_proven_upper(n_t, w_t) * shave < c_exact(n_t, w_t),
       "teeth: envelope shaved to %.6f x at its tightest point breaks domination"
       % shave)

    # DIRECTION DISCIPLINE for B3 (was prose-only until 2026-09-15; a reviewer
    # could not reproduce it because no check enforced it -- now it is a gate).
    # Im psi enters (C1) with a MINUS sign, so an UPPER bound on c_n needs Im
    # from BELOW.  B3 (shipped):      Im psi >= atan(w/x) - w/(2|z|^2)   VALID
    # Mirrored upper companion:       Im psi <= atan(w/x) + w/(2|z|^2)   FALSE
    # (the next Stirling term -Im(1/(12 z^2)) is positive here and wins).
    #
    # SCOPE of that falsity -- it is about SYMMETRISING B3, not about upper
    # bounds on Im psi in general.  Audited 2026-09-15 against INDEPENDENT
    # mpmath (dps=40, not _psi_complex) on a 240-node grid:
    #     atan(w/x) + w/(2q)   -> 240/240 violations  (the shipped mirror)
    #     atan(w/x)            -> 240/240 violations
    #     atan(w/x) + w/q      ->   0/240  VALID upper bound
    #     atan(w/(x-1))        ->   0/240  VALID upper bound
    # So a reviewer who reconstructs "the" mirror as a doubled correction or an
    # x-shift will see 0 violations and think the claim fails; the claim is
    # specifically that flipping the sign of B3's own correction term is FALSE.
    # Asserting BOTH sides keeps the one-sided house rule honest: if a future
    # edit "symmetrises" B3, the second check below fails loudly.
    _b3_lo_vio = _b3_up_vio = _b3_tot = 0
    _b3_worst = 0.0
    for w in np.geomspace(0.5, 400.0, 60):
        for n in range(0, 60):
            x = n + 1.0
            q = x * x + w * w
            im_true = float(np.imag(_psi_complex(np.array([x]), w)[0]))
            _b3_tot += 1
            if im_true < np.arctan(w / x) - w / (2.0 * q) - 1e-13:
                _b3_lo_vio += 1
            _d = im_true - (np.arctan(w / x) + w / (2.0 * q))
            if _d > 1e-13:
                _b3_up_vio += 1
                _b3_worst = max(_b3_worst, _d)
    print("      B3 lower (shipped) violations %d/%d ; mirrored upper %d/%d "
          "(worst %+.3e)"
          % (_b3_lo_vio, _b3_tot, _b3_up_vio, _b3_tot, _b3_worst))
    ck(_b3_lo_vio == 0,
       "B3 LOWER bound on Im psi holds at all %d nodes (the direction (C1) needs)"
       % _b3_tot)
    ck(_b3_up_vio == _b3_tot,
       "mirrored Im UPPER companion is FALSE at all %d/%d nodes (worst %+.3e) -- "
       "one-sided house rule, sides stay separate"
       % (_b3_up_vio, _b3_tot, _b3_worst))

    # (C2) coth law.  a*c_n -> (pi/2)coth(pi w) with an O(1/n) correction, so a
    # single n cannot be pinned to a fixed tolerance; assert the DECAY RATE:
    # doubling n must halve the residual (ratio in [1.6, 2.4]).
    for w in (0.7, 4.0, 40.0):
        lim = 0.5 * np.pi * _coth(np.pi * w)
        e1 = abs(c_exact(100000, w) * ((100000 + 0.5) / w) - lim)
        e2 = abs(c_exact(200000, w) * ((200000 + 0.5) / w) - lim)
        rat = e1 / max(e2, 1e-300)
        print("      w=%-6.4g limit=%.9f  res(1e5)=%.3e res(2e5)=%.3e ratio=%.3f"
              % (w, lim, e1, e2, rat))
        ck(1.6 < rat < 2.4,
           "coth law w=%.3g: residual halves under n->2n (ratio %.3f)" % (w, rat))

    # NOTE (finding, 2026-09-15; CORRECTED after review 2026-09-15):
    # Theta* is NOT a valid pointwise fill for the row-sum functional c_n: the
    # TRUE row sum exceeds it at small w.  Measured with c_exact (the (C1)
    # identity), sup_n c_n(w) vs Theta* = 0.6349287918:
    #     w=0.5 -> 0.769610503 (+1.35e-01)   w=1.0 -> 0.660434453 (+2.55e-02)
    #     w=1.5 -> 0.638525502 (+3.60e-03)   w=3.0 -> 0.637157183 (+2.23e-03)
    # and it settles BELOW Theta* for w >~ 20 (w=26 -> 0.634896511, -3.2e-05;
    # w=100 -> 0.634927676, -1.1e-06), approaching Theta* from below as w->oo.
    #
    # CORRECTION: an earlier version of this block sup'd c_proven_upper (the
    # LOOSE envelope) and reported w=26 -> 0.64333 as a row-sum excess.  That
    # number is envelope SLACK, not row-sum value: at the w=26 peak n=39,
    # c_proven_upper - c_exact = 8.07e-03 while c_exact itself is 3.2e-05
    # BELOW Theta*.  The retired claim "the old pad C_SUP=0.6375 was absorbing
    # this finite-w excess" was likewise a misdiagnosis at w=26.  The LOGIC is
    # unaffected and still holds on the corrected (small-w) evidence.
    #
    # Nothing in the brackets ever depended on this: Theta* bounds the OPERATOR
    # NORM ||calR|| (Lemma M, head/analytic routes), not the row-sum functional,
    # and h_hi_dyadic subtracts its own per-block sup, never Theta*.
    # The peak sits at n ~ 1.5 w, so any search window MUST scale with w.
    def _sup_rows_exact(w, span=4.0, pts=1200):
        """TRUE sup_n c_n(w) from the (C1) identity -- NOT the envelope."""
        hi = max(64.0, span * w)
        ns = np.unique(np.linspace(0.0, hi, pts).astype(np.int64))
        return max(float(c_exact(int(n), w)) for n in ns)

    for w_ex, lo_ex in ((0.5, 1.3e-1), (1.0, 2.5e-2), (3.0, 2.2e-3)):
        s_ex = _sup_rows_exact(w_ex)
        print("      w=%-5.3g true sup_n c_n = %.9f  vs Theta*: %+.3e"
              % (w_ex, s_ex, s_ex - THETA))
        ck(s_ex - THETA > lo_ex,
           "TRUE row sum at w=%.3g exceeds Theta* by %.3e (Theta* is not a valid "
           "row-sum fill)" % (w_ex, s_ex - THETA))

    # ... and the complementary honest half: for large w the true sup is BELOW
    # Theta*, so the excess is a small-w phenomenon only.  This is the assertion
    # the retired w=26 "excess" got backwards.
    for w_bl in (26.0, 100.0):
        s_bl = _sup_rows_exact(w_bl)
        ck(s_bl < THETA,
           "TRUE row sum at w=%.3g is BELOW Theta* (%.9f < %.9f) -- the excess is "
           "a small-w effect" % (w_bl, s_bl, THETA))

    # teeth: the quantity above must be the identity, not the envelope.  If a
    # future edit swaps c_exact back to c_proven_upper here, the envelope slack
    # at the w=26 peak (8.07e-03) would masquerade as a row-sum excess.
    _n26 = 39
    ck(c_proven_upper(_n26, 26.0) - c_exact(_n26, 26.0) > 5e-3
       and c_exact(_n26, 26.0) < THETA,
       "teeth: envelope slack at (n=39,w=26) is %.3e and c_exact is BELOW Theta* "
       "-- sup'ing the envelope here would fake an excess"
       % (c_proven_upper(_n26, 26.0) - c_exact(_n26, 26.0)))

    # The rebuilt Prop. B (proven envelope in place of the calibrated c_shape)
    # must still enclose h_true on the dispatched band.
    # DIRECTION NOTE (house rule: audit every truncation).  H_mat(N,w) interlaces
    # upward in N, so h^(N) <= h^(inf): shrinking N makes `ht <= hi` EASIER (a
    # weaker test of the upper bound) and `lo <= ht` HARDER.  Full mode therefore
    # keeps N=6000; quick mode drops to 2500 (agrees to 7 dp at w=30, 1.1s vs
    # 17s per node) and is a deliberately weaker upper-bound probe, not a proof.
    _n8e = 2500 if quick else 6000
    print("    [8e] Prop. B with the proven envelope still encloses h_true"
          "  (N=%d%s)" % (_n8e, ", quick: weaker upper probe" if quick else ""))
    for w in (30.0, 50.0, 200.0, 1000.0):
        hi = h_hi_dyadic(w)
        lo = h_lo_analytic(w)
        ht = float(np.max(np.linalg.eigvalsh(H_mat(_n8e, w))))
        ck(lo <= ht <= hi, "[8e] w=%.0f: %.7f <= %.7f <= %.7f (width %.4f)"
           % (w, lo, ht, hi, (hi - lo) / lo))

    print("\nTIME %.1fs   PASS %d  FAIL %d" % (time.time() - t0, PASS[0], FAIL[0]))
    return rows


if __name__ == "__main__":
    K = 8
    # NB: match the bare word AND the dashed forms -- a plain `"quick" in sys.argv`
    # silently misses `--quick` and runs the full-density grid (regressed once
    # during the 2026-09-15 rebuild; caught by a 900s timeout).
    quick = any(a in ("quick", "--quick", "-q") for a in sys.argv[1:])
    for a in sys.argv[1:]:
        if a.startswith("K="):
            K = int(a[2:])
    main(K=K, quick=quick)
