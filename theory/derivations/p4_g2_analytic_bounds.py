#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
p4_g2_analytic_bounds.py  --  regression gate for the ANALYTIC two-sided bound on

        g(2) := lambda_1( G(2) ),      G(2) := O[k -> 1/k^2],

the constant that carries Theorem 1(c) of `P4_main_theorems.md`
(Tc < Tctilde = (1/2pi) sqrt( g(2) <w^2> lambda )).

The companion write-up is `P4_g2_analytic.md`.  This script is NOT the proof;
it is the regression gate: every closed form, every inequality and every
constant of the write-up is recomputed here from the *elementwise* definition
of O[f] given in `P4_main_theorems.md` §0, at several truncations N, and
compared grid-by-grid.  The certificate arithmetic itself (Stage 4/5) is done
in exact rationals (`fractions.Fraction`) -- no floating point, no rounding
pads, no interval arithmetic needed.

Run:  OMP_NUM_THREADS=1 ... python p4_g2_analytic_bounds.py [quick|full]

====================================================================
                        SETTING
====================================================================
u_n := (2n+1)^{-1/2},  n >= 0.   For f : N_{>=1} -> R,

  O[f]_{nm} = - delta_{nm} * 2 u_n^2 * sum_{k=1}^{n} f(k)
              + (1-delta_{nm}) * f(|n-m|) * u_n u_m
              + f(n+m+1) * u_n u_m.

A := G(2) = O[1/k^2]:

  A_{nn} = u_n^2 ( 1/(2n+1)^2 - 2 H_n^{(2)} ),     H_n^{(2)} := sum_{k=1}^n 1/k^2,
  A_{nm} = u_n u_m ( 1/(n-m)^2 + 1/(n+m+1)^2 )     (n != m).

====================================================================
        LEMMA 0  (A is Hilbert-Schmidt, hence compact)
====================================================================
||A||_HS^2 = sum_{n,m} A_{nm}^2 < inf, so A is a compact self-adjoint operator on
l^2 and g(2) := lambda_1(A) = max spec A is an attained eigenvalue.

PROOF.  Diagonal: |1/(2n+1)^2 - 2H_n^{(2)}| <= 2 zeta(2) for every n >= 0
(the bracket lies in [-2 zeta(2), 1]), so A_{nn}^2 <= 4 zeta(2)^2/(2n+1)^2 and
  sum_n A_{nn}^2 <= 4 zeta(2)^2 * pi^2/8 = zeta(2)^2 pi^2/2.
Off-diagonal, with (x+y)^2 <= 2x^2+2y^2 and 2/((2n+1)(2m+1)) <= 1/(2n+1)^2 + 1/(2m+1)^2,
  sum_{n!=m} A_{nm}^2 <= sum_{n!=m} (1/(2n+1)^2 + 1/(2m+1)^2)(1/(n-m)^4 + 1/(n+m+1)^4)
       <= 2 sum_n 1/(2n+1)^2 * ( 2 zeta(4) + zeta(4) ) = 3 zeta(4) pi^2/4 .
Both are finite.  (The numerical value is ||A||_HS^2 = 3.8496...; the bound above
is crude, only finiteness is used.)   []

Throughout, for v in l^2 we write

        p_n := u_n v_n         (equivalently v_n = sqrt(2n+1) p_n),

so that  ||v||^2 = sum_n (2n+1) p_n^2 , and every bilinear expression below is
a *rational* function of the p's (all the irrational u_n's cancel).

====================================================================
        LEMMA 1  (exact rearrangement identity)
====================================================================
For every finitely supported real v (p_n = u_n v_n):

  <v,Av> = <v,Bv> - (1/2) D[p],                                        (L1)

  <v,Bv> := sum_{n,m>=0} p_n p_m/(n+m+1)^2 + sum_{n>=0} psi'(n+1) p_n^2 ,
  D[p]   := sum_{n != m >= 0} (p_n - p_m)^2/(n-m)^2  >= 0,
  psi'(n+1) = sum_{k>n} 1/k^2 = zeta(2) - H_n^{(2)} .

i.e. B is the operator with entries

  B_{nm} = u_n u_m/(n+m+1)^2 + delta_{nm} u_n^2 psi'(n+1).

PROOF.  Substituting u_n v_n = p_n in the elementwise definition,

  <v,Av> = sum_{n!=m} p_n p_m/(n-m)^2 + sum_{n,m} p_n p_m/(n+m+1)^2
           - 2 sum_n H_n^{(2)} p_n^2 .

For the first sum use the polarisation  2 p_n p_m = p_n^2 + p_m^2 - (p_n-p_m)^2
and  c_n := sum_{m>=0, m!=n} 1/(n-m)^2 = zeta(2) + H_n^{(2)}
(the m>n half gives zeta(2), the m<n half gives H_n^{(2)}):

  sum_{n!=m} p_n p_m/(n-m)^2 = sum_n c_n p_n^2 - (1/2) D[p].

Then  c_n - 2 H_n^{(2)} = zeta(2) - H_n^{(2)} = psi'(n+1).  All sums converge
absolutely for finitely supported p.   []

CONSEQUENCE (L1'):   <v,Av> <= <v,Bv>  for all finitely supported v, hence -- both
forms being continuous (A is bounded, B is trace class, see Lemma 3) -- for all
v in l^2.  Equivalently  A <= B  in the quadratic-form order.

====================================================================
        LEMMA 2  (B >= 0)
====================================================================
B >= 0.   PROOF.  1/(n+m+1)^2 = int_0^inf t e^{-(n+m+1)t} dt, hence for finitely
supported p

  sum_{n,m} p_n p_m/(n+m+1)^2 = int_0^inf t e^{-t} ( sum_n p_n e^{-nt} )^2 dt >= 0,

and sum_n psi'(n+1) p_n^2 >= 0.  Density + boundedness of B finishes.   []

(So the Hankel kernel 1/(n+m+1)^2 is a moment matrix; this is the only place
where positivity of L3 is used.)

====================================================================
        LEMMA 3  (trace of B, and the gap constant beta)
====================================================================
  tr B = sum_{n>=0} B_{nn} = sum_{n>=0} [ 1/(2n+1)^3 + psi'(n+1)/(2n+1) ] < inf,

so B is trace class; being >= 0 all its eigenvalues are >= 0 and sum to tr B.
Since  lambda_1(B) >= B_{00} = 1 + zeta(2),

  lambda_2(B) <= sum_{i>=2} lambda_i(B) = tr B - lambda_1(B) <= tr B - B_{00}
               = sum_{n>=1} [ 1/(2n+1)^3 + psi'(n+1)/(2n+1) ] =: beta_0 .

Two classical closed forms and one elementary bound finish it.  Everything below
uses only the MIDPOINT RULE (Lemma 4b): for f convex on (c,inf) and a-1/2 > c,
f(a) <= int_{a-1/2}^{a+1/2} f, hence sum_{k>=K} f(k) <= int_{K-1/2}^inf f.

  sum_{n>=1} 1/(2n+1)^3 = (7/8) zeta(3) - 1                    [classical]
  psi'(n+1) = sum_{k>n} 1/k^2 <= int_{n+1/2}^inf dx/x^2 = 1/(n+1/2) = 2/(2n+1),
  hence sum_{n>=1} psi'(n+1)/(2n+1) <= sum_{n>=1} 2/(2n+1)^2 = 2(pi^2/8 - 1) = pi^2/4 - 2.

Therefore
        beta_0 <= (7/8) zeta(3) + pi^2/4 - 3 = (7/8) zeta(3) + (3/2) zeta(2) - 3 .

To certify that this is < 13/25 WITHOUT quoting decimal expansions of zeta(2),
zeta(3), bound both zeta values by their 10-term partial sum plus a midpoint tail:

  zeta(2) <= H_10^{(2)} + 1/(10+1/2) = 2089289/1270080          =: R2
  zeta(3) <= H_10^{(3)} + 1/(2(10+1/2)^2) = 19236689947/16003008000  =: R3
  (7/8) R3 + (3/2) R2 - 3 = 9497876347/18289152000 = 0.51931748...  <  13/25 .

This is a single inequality between explicit rationals, checkable by hand.

Combining with Lemma 1' and Courant-Fischer (lambda_k(A) <= lambda_k(B) for the
min-max values of two compact self-adjoint operators with A <= B):

  (L3')   lambda_2(A) <= lambda_2(B) <= beta := 13/25 .

SPECTRAL SEPARATION (the hypothesis Temple actually needs).  Let P be the spectral
projection of A onto (beta, inf).  A is compact (Lemma 0) and beta > 0, so P has
finite rank = #{eigenvalues > beta}.  If that rank were >= 2, choosing any 1-dim
subspace S inside ran P leaves a vector v in ran P with v _|_ S, ||v||=1 and
<v,Av> > beta; taking the min over S still gives lambda_2(A) > beta, contradicting
(L3').  So the rank is <= 1, and lambda_1 >= rho > beta forces it to be exactly 1:

        spec(A) \cap (13/25, inf) = {lambda_1},   lambda_1 simple.

(This is the step that turns "lambda_2 <= beta" into "everything except lambda_1
is <= beta"; it is not automatic and is where compactness is used.)   []

====================================================================
        LEMMA 4  (kernel comparison and two tail sums)
====================================================================
(4a)  For integers 0 <= m < N1 <= n:   (n-m)(N1+1) >= (n+1)(N1-m).
      PROOF. (n-m)(N1+1) - (n+1)(N1-m) = (n-M... ) expand:
        nN1 + n - mN1 - m - ( nN1 - nm + N1 - m ) = n + nm - mN1 - N1
        = (n - N1)(1 + m) >= 0.    []
      Hence   1/(n-m)^2 <= (N1+1)^2 / ( (n+1)^2 (N1-m)^2 ).

(4b)  MIDPOINT RULE.  If f is convex on an open interval containing
      [a-1/2, a+1/2] then f(a) <= int_{a-1/2}^{a+1/2} f(x) dx.
      PROOF.  Convexity gives f(x) >= f(a) + f'_+(a)(x-a); integrating over the
      interval, which is symmetric about a, kills the linear term.   []
      COROLLARY.  If f is convex on (c,inf) and K-1/2 > c then
              sum_{k>=K} f(k) <= int_{K-1/2}^inf f(x) dx.

(4c)  T(N1) := sum_{n>=N1} 1/( (n+1)^4 (2n+1) )  <=  1/( N1 (2N1+1)^3 )   (N1 >= 1).
      PROOF.  Let h(x) = 1/((x+1)^4 (2x+1)) on (-1/2, inf).  Then
        (log h)'' = 4/(x+1)^2 + 4/(2x+1)^2 > 0,
      so h is log-convex, hence convex.  By (4b),
        T(N1) <= int_{N1-1/2}^inf h(x) dx = int_a^inf dt/(2t^5 - t^4),   a := N1+1/2
      (substituting t = x+1 and using (x+1)^4(2x+1) = 2t^5 - t^4).  For t >= a,
        2t^5 - t^4 = 2t^5 (1 - 1/(2t)) >= 2t^5 (1 - 1/(2a)) = 2t^5 (2a-1)/(2a),
      so
        T(N1) <= (2a/(2a-1)) * (1/2) int_a^inf t^{-5} dt = (2a/(2a-1)) / (8 a^4)
               = 1/( 4 a^3 (2a-1) ) = 8 / ( 4 (2N1+1)^3 * 2N1 ) = 1/( N1 (2N1+1)^3 ).  []
      (Sharp to O(1/N1): the ratio T(N1)*N1*(2N1+1)^3 increases to 1.  Stage 3
      checks it for 1 <= N1 <= 400 and at N1 = 1000, 5000.)

====================================================================
        THEOREM (the certificate)
====================================================================
Fix an integer M >= 1, an integer N1 >= M and positive rationals p_0,...,p_{M-1}.
Let v be the l^2 vector v_n = sqrt(2n+1) p_n for n < M, v_n = 0 for n >= M.
Define the *rational* quantities

  ||v||^2 = sum_{n<M} (2n+1) p_n^2,
  Q       = sum_{n,m<M} p_n p_m K(n,m) - 2 sum_{n<M} H_n^{(2)} p_n^2,
              K(n,m) := 1/(n+m+1)^2 + [n != m]/(n-m)^2,
  rho     = Q / ||v||^2,
  g_n     = sum_{m<M} p_m K(n,m) - 2 H_n^{(2)} p_n [n<M]           (so (Av)_n = u_n g_n),
  R       = sum_{n<M} (g_n - rho (2n+1) p_n)^2/(2n+1)
            + sum_{M<=n<N1} g_n^2/(2n+1)
            + W^2 / ( N1 (2N1+1)^3 ),
  W       = sum_{m<M} p_m ( (N1+1)^2/(N1-m)^2 + 1 ),
  r^2     = R/||v||^2.

Then, provided rho > 13/25,

      rho  <=  g(2)  <=  rho + r^2/(rho - 13/25).                     (*)

PROOF.  Lower half: Rayleigh, g(2) = max spec A >= <v,Av>/||v||^2 = rho.

Upper half:  ||Av - rho v||^2 = sum_n ( u_n g_n - rho sqrt(2n+1) p_n )^2
= sum_n (g_n - rho(2n+1)p_n)^2/(2n+1).  For n >= M the trial vector vanishes, so
the summand is g_n^2/(2n+1) with g_n = sum_{m<M} p_m K(n,m) > 0; by Lemma 4a and
1/(n+m+1)^2 <= 1/(n+1)^2 we get, for n >= N1,
        0 < g_n <= W/(n+1)^2,
hence  sum_{n>=N1} g_n^2/(2n+1) <= W^2 T(N1) <= W^2 / (N1 (2N1+1)^3) (Lemma 4c).
Therefore ||Av - rho v||^2 <= R and r^2 >= ||Av-rho v||^2/||v||^2.

Now Temple's inequality.  Let lambda_1 = g(2) and beta = 13/25.  By (L3')
the spectral measure of A on the unit vector w := v/||v|| is supported in
{lambda_1} u (-inf, beta], so

  0 <= <w, (A-lambda_1)(A-beta) w> = ||Aw||^2 - (lambda_1+beta) rho + lambda_1 beta
     = r^2 + rho^2 - (lambda_1+beta) rho + lambda_1 beta
     = r^2 + rho(rho-beta) - lambda_1 (rho - beta),
(using ||Aw||^2 = r^2 + rho^2 which is exactly the definition of the residual),
and rho - beta > 0 gives lambda_1 <= rho + r^2/(rho-beta).   []

[Each factor is non-positive/non-negative as required: for a spectral point mu
!= lambda_1 we have mu <= beta < rho <= lambda_1, so (mu-lambda_1)(mu-beta) >= 0;
at mu = lambda_1 the first factor vanishes.]

Note the two ingredients that make (*) sharp:
  * Temple is *second order* in the trial-vector error (width = r^2/(rho-beta)),
    whereas the Collatz-Wielandt / weighted Schur bound sup_n (Aw)_n/w_n is only
    first order.  See Stage 7 for the numerical comparison.
  * beta = 13/25 is crude (the true lambda_2(A) is 5.354e-4) but costs only the
    factor (rho-lambda_2)/(rho-beta) = 1.65 in the width.

====================================================================
        CERTIFICATES SHIPPED
====================================================================
All trial vectors are POSITIVE INTEGER vectors (the certificate is invariant
under p -> c p, so integers suffice; an overall scale is pure bookkeeping).

CERT-A  (hand-checkable) : M=3,   N1=12,  p = (200, 45, 13)
CERT-A' (hand-checkable) : M=4,   N1=16,  p = (1000, 228, 69, 25)
CERT-B  (headline)       : M=12,  N1=40,  p = (10000, 2286, 704, 275, 129, 70,
                                               41, 27, 18, 13, 9, 7)
CERT-C  (extended)       : M=30,  N1=80   (integer vector generated in-script)
CERT-D  (limit demo)     : M=200, N1=500  (integer vector generated in-script)

The width scales like M^{-4}; see the scaling table in Stage 5.
See `P4_g2_analytic.md` for the write-up and the honesty-boundary section.
"""

import os
import sys
import time
from fractions import Fraction as F

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import mpmath as mp

mp.mp.dps = 40

# certified enclosure of g(2) from P4_main_theorems.md §4.1 (mpmath Rayleigh + Kato-Temple)
G2_LO_REF = mp.mpf("1.318140211762671")
G2_HI_REF = mp.mpf("1.318140212136387")
G2_REF = mp.mpf("1.3181402121321")          # best float64/mpmath value, N -> inf
G2_LO_REF_F = F(1318140211762671, 10 ** 15)
G2_HI_REF_F = F(1318140212136387, 10 ** 15)
G2_REF_F = F(13181402121321, 10 ** 13)

NCHECK = 0
NFAIL = 0
_T0 = time.time()


def chk(cond, label, extra=""):
    global NCHECK, NFAIL
    NCHECK += 1
    if not cond:
        NFAIL += 1
        print(f"    !! FAIL  [{NCHECK:4d}] {label}  {extra}")
    return bool(cond)


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# ---------------------------------------------------------------------------
# Stage 0 : elementwise definition, two independent constructions
# ---------------------------------------------------------------------------

def build_A_elementwise(N):
    """Literal transcription of O[f]_{nm} with f(k)=1/k^2 (double loop, slow, N small)."""
    A = np.zeros((N, N))
    u = np.array([1.0 / np.sqrt(2 * n + 1.0) for n in range(N)])
    H = np.zeros(N)
    for n in range(1, N):
        H[n] = H[n - 1] + 1.0 / n ** 2
    for n in range(N):
        for m in range(N):
            val = 0.0
            if n == m:
                val -= 2 * u[n] ** 2 * H[n]
            else:
                val += u[n] * u[m] / (n - m) ** 2
            val += u[n] * u[m] / (n + m + 1) ** 2
            A[n, m] = val
    return A, u, H


def build_A(N):
    """Vectorised construction (used for the large-N gates)."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    off = np.zeros((N, N))
    nz = D > 0
    off[nz] = 1.0 / D[nz].astype(float) ** 2
    A = (off + 1.0 / S.astype(float) ** 2) * u[:, None] * u[None, :]
    H = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, N) ** 2)))
    A[n, n] -= 2 * u ** 2 * H
    return A, u, H


def build_B(N, H):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    S = n[:, None] + n[None, :] + 1
    B = (1.0 / S.astype(float) ** 2) * u[:, None] * u[None, :]
    zeta2 = float(mp.zeta(2))
    B[n, n] += (zeta2 - H[:N]) * u ** 2
    return B


def stage0():
    hdr("Stage 0 -- elementwise definition of G(2); lambda_1 at N = 400/800/1600")
    A1, u1, H1 = build_A_elementwise(60)
    A2, u2, H2 = build_A(60)
    chk(np.max(np.abs(A1 - A2)) < 1e-15, "elementwise vs vectorised build (N=60)",
        f"maxdiff={np.max(np.abs(A1 - A2)):.2e}")
    chk(np.max(np.abs(A2 - A2.T)) < 1e-16, "A symmetric")
    chk(abs(A2[0, 0] - 1.0) < 1e-15, "A_00 = 1", f"{A2[0,0]}")
    chk(abs(A2[1, 1] - (1.0 / 3) * (1.0 / 9 - 2)) < 1e-15, "A_11 = (1/3)(1/9-2)")

    out = {}
    hs_bound = float(mp.zeta(2) ** 2 * mp.pi ** 2 / 2 + 3 * mp.zeta(4) * mp.pi ** 2 / 4)
    print(f"    Lemma 0 analytic HS bound: ||A||_HS^2 <= zeta(2)^2 pi^2/2 + 3 zeta(4) pi^2/4"
          f" = {hs_bound:.4f}")
    for N in (400, 800, 1600):
        A, u, H = build_A(N)
        w = np.linalg.eigvalsh(A)
        out[N] = (w[-1], w[-2], w[0], float(np.sum(A * A)))
        print(f"    N={N:5d}  lam1={w[-1]:.13f}  lam2={w[-2]:.6e}  lam_min={w[0]:.9f}"
              f"  ||A||_HS^2={np.sum(A*A):.6f}")
        chk(abs(w[-1] - float(G2_REF)) < 3e-9, f"lambda_1(A_N) ~ g(2), N={N}",
            f"{w[-1]:.13f}")
        chk(w[0] > -0.99, f"lambda_min(A_N) > -0.99, N={N}", f"{w[0]:.9f}")
        chk(w[-2] < 6e-4, f"lambda_2(A_N) < 6e-4, N={N}", f"{w[-2]:.3e}")
        chk(float(np.sum(A * A)) <= hs_bound, f"Lemma 0: ||A_N||_HS^2 <= analytic bound, N={N}",
            f"{np.sum(A*A):.4f} <= {hs_bound:.4f}")
    return out


# ---------------------------------------------------------------------------
# Stage 1 : Lemma 1, the rearrangement identity
# ---------------------------------------------------------------------------

def psi1(n):
    """psi'(n) = sum_{k>=n} 1/k^2, exact mpmath."""
    return mp.polygamma(1, n)


def stage1():
    hdr("Stage 1 -- Lemma 1: <v,Av> = <v,Bv> - (1/2) D[p]   (exact, mpmath dps=40)")
    rng = np.random.default_rng(20260912)
    zeta2 = mp.zeta(2)
    Amat, _, Hf = build_A(2200)
    for trial, M in enumerate([1, 2, 3, 5, 8, 13, 21, 34]):
        praw = rng.normal(size=M)
        p = [mp.mpf(float(x)) for x in praw]
        # LHS, from the elementwise definition (exact: v has support < M)
        v = np.array([float(p[n]) * np.sqrt(2 * n + 1.0) for n in range(M)])
        lhs_mat = float(v @ Amat[:M, :M] @ v)
        # LHS, from the p-form
        H = [mp.mpf(0)] * M
        for n in range(1, M):
            H[n] = H[n - 1] + mp.mpf(1) / n ** 2
        lhs = mp.mpf(0)
        for n in range(M):
            for m in range(M):
                k = mp.mpf(1) / (n + m + 1) ** 2
                if n != m:
                    k += mp.mpf(1) / (n - m) ** 2
                lhs += p[n] * p[m] * k
        for n in range(M):
            lhs -= 2 * H[n] * p[n] ** 2
        chk(abs(lhs - lhs_mat) < 1e-11 * max(1.0, abs(lhs_mat)),
            f"p-form of <v,Av> vs matrix (M={M})", f"{float(lhs)} vs {lhs_mat}")
        # RHS
        Bform = mp.mpf(0)
        for n in range(M):
            for m in range(M):
                Bform += p[n] * p[m] / mp.mpf((n + m + 1) ** 2)
        for n in range(M):
            Bform += (zeta2 - H[n]) * p[n] ** 2
        Dir = mp.mpf(0)
        for n in range(M):
            for m in range(M):
                if n != m:
                    Dir += (p[n] - p[m]) ** 2 / mp.mpf((n - m) ** 2)
        for n in range(M):                      # m >= M half, twice
            Dir += 2 * p[n] ** 2 * psi1(M - n)
        rhs = Bform - Dir / 2
        chk(abs(lhs - rhs) < mp.mpf("1e-30") * max(1, abs(lhs)),
            f"Lemma 1 identity (M={M})", f"diff={float(lhs-rhs):.3e}")
        chk(Dir >= 0, f"D[p] >= 0 (M={M})")
        chk(lhs <= Bform + mp.mpf("1e-30"), f"Lemma 1' : <v,Av> <= <v,Bv> (M={M})")
    # non-negative trial vectors too
    for M in (4, 9, 16):
        p = [mp.mpf(float(abs(x))) + mp.mpf("0.01") for x in rng.normal(size=M)]
        H = [mp.mpf(0)] * M
        for n in range(1, M):
            H[n] = H[n - 1] + mp.mpf(1) / n ** 2
        lhs = sum(p[n] * p[m] * ((mp.mpf(1) / (n + m + 1) ** 2) +
                                 (mp.mpf(1) / (n - m) ** 2 if n != m else 0))
                  for n in range(M) for m in range(M)) \
            - 2 * sum(H[n] * p[n] ** 2 for n in range(M))
        Bform = sum(p[n] * p[m] / mp.mpf((n + m + 1) ** 2) for n in range(M) for m in range(M)) \
            + sum((zeta2 - H[n]) * p[n] ** 2 for n in range(M))
        chk(lhs <= Bform, f"A <= B on a positive vector (M={M})")


# ---------------------------------------------------------------------------
# Stage 2 : Lemma 2 (B >= 0) and Lemma 3 (trace / beta)
# ---------------------------------------------------------------------------

BETA = F(13, 25)


def stage2():
    hdr("Stage 2 -- Lemma 2 (B >= 0), Lemma 3 (trace bound, beta = 13/25)")
    rng = np.random.default_rng(7)
    # (a) moment representation of the Hankel kernel
    for M in (3, 6, 12):
        p = rng.normal(size=M)
        direct = sum(p[n] * p[m] / (n + m + 1) ** 2 for n in range(M) for m in range(M))
        f = lambda t: t * mp.e ** (-t) * (sum(p[n] * mp.e ** (-n * t) for n in range(M))) ** 2
        integ = float(mp.quad(f, [0, mp.inf]))
        chk(abs(direct - integ) < 1e-12 * max(1, abs(direct)),
            f"Hankel moment representation (M={M})", f"{direct} vs {integ}")
        chk(direct >= -1e-14, f"Hankel form >= 0 (M={M})")
    # (b) spectrum of B
    for N in (400, 800, 1600):
        _, _, H = build_A(N)
        B = build_B(N, H)
        wB = np.linalg.eigvalsh(B)
        chk(wB[0] > -1e-12, f"lambda_min(B_N) >= 0 (N={N})", f"{wB[0]:.3e}")
        A, _, _ = build_A(N)
        wD = np.linalg.eigvalsh(B - A)
        chk(wD[0] > -1e-9, f"B - A >= 0 (N={N})", f"{wD[0]:.3e}")
        wA = np.linalg.eigvalsh(A)
        chk(wA[-2] <= wB[-2] + 1e-12, f"lambda_2(A) <= lambda_2(B) (N={N})",
            f"{wA[-2]:.6e} vs {wB[-2]:.6e}")
        chk(wB[-2] <= float(BETA), f"lambda_2(B_N) <= 13/25 (N={N})", f"{wB[-2]:.6f}")
        chk(wA[-2] <= float(BETA), f"lambda_2(A_N) <= 13/25 (N={N})", f"{wA[-2]:.3e}")
        print(f"    N={N:5d}  lam1(B)={wB[-1]:.9f} lam2(B)={wB[-2]:.9f} "
              f"trB_N={np.trace(B):.9f}")
    # (c) psi'(n+1) <= 2/(2n+1)  (== 1/(n+1/2)); mp.zeta(2, n+1) is the Hurwitz form
    bad = 0
    for n in list(range(0, 3001)) + [10 ** 4, 10 ** 5, 10 ** 6, 10 ** 9]:
        if mp.zeta(2, n + 1) > mp.mpf(2) / (2 * n + 1):
            bad += 1
    chk(bad == 0, "psi'(n+1) <= 2/(2n+1) for n <= 3000 and n in {1e4,1e5,1e6,1e9}",
        f"violations={bad}")
    # (d) beta_0 and the closed forms
    beta0 = mp.nsum(lambda n: mp.mpf(1) / (2 * n + 1) ** 3 + mp.zeta(2, n + 1) / (2 * n + 1),
                    [1, mp.inf])
    closed = mp.mpf(7) / 8 * mp.zeta(3) + mp.pi ** 2 / 4 - 3
    print(f"    beta_0 = sum_{{n>=1}} B_nn = {mp.nstr(beta0, 15)}")
    print(f"    (7/8)zeta(3) + pi^2/4 - 3   = {mp.nstr(closed, 15)}")
    print(f"    tr B = {mp.nstr(beta0 + 1 + mp.zeta(2), 15)}   B_00 = 1+zeta(2) = "
          f"{mp.nstr(1 + mp.zeta(2), 15)}")
    chk(beta0 <= closed, "beta_0 <= (7/8)zeta(3)+pi^2/4-3",
        f"{mp.nstr(beta0,12)} vs {mp.nstr(closed,12)}")
    chk(closed < mp.mpf(BETA.numerator) / BETA.denominator, "(7/8)zeta(3)+pi^2/4-3 < 13/25",
        f"{mp.nstr(closed,12)}")
    # rational majorisation actually used in the write-up:
    #   zeta(2) <= H_10^(2) + 1/(10+1/2),  zeta(3) <= H_10^(3) + 1/(2(10+1/2)^2)
    # (midpoint rule, Lemma 4b) -- no decimal expansion of zeta is quoted.
    NZ = 10
    R2 = sum(F(1, k * k) for k in range(1, NZ + 1)) + F(2, 2 * NZ + 1)
    R3 = sum(F(1, k ** 3) for k in range(1, NZ + 1)) + F(2, (2 * NZ + 1) ** 2)
    chk(mp.zeta(2) <= mp.mpf(R2.numerator) / R2.denominator,
        "zeta(2) <= H_10^(2) + 2/21 (midpoint tail)", f"{R2}")
    chk(mp.zeta(3) <= mp.mpf(R3.numerator) / R3.denominator,
        "zeta(3) <= H_10^(3) + 2/441 (midpoint tail)", f"{R3}")
    rat = F(7, 8) * R3 + F(3, 2) * R2 - 3
    chk(rat < BETA, "(7/8)R3 + (3/2)R2 - 3 < 13/25 (exact rationals)",
        f"{rat} = {float(rat):.10f}")
    print(f"    rational certificate: (7/8)R3+(3/2)R2-3 = {rat} = {float(rat):.10f} < 13/25")
    # tail bounds behind R2, R3 (midpoint rule), checked against mpmath
    chk(mp.zeta(2, NZ + 1) <= mp.mpf(2) / (2 * NZ + 1),
        "sum_{k>10} 1/k^2 <= 1/(10+1/2)")
    chk(mp.zeta(3, NZ + 1) <= mp.mpf(2) / (2 * NZ + 1) ** 2,
        "sum_{k>10} 1/k^3 <= 1/(2(10+1/2)^2)")
    # closed forms sum_{n>=1} 1/(2n+1)^3 and sum_{n>=1} 1/(2n+1)^2
    s3 = mp.nsum(lambda n: mp.mpf(1) / (2 * n + 1) ** 3, [1, mp.inf])
    s2 = mp.nsum(lambda n: mp.mpf(1) / (2 * n + 1) ** 2, [1, mp.inf])
    chk(abs(s3 - (mp.mpf(7) / 8 * mp.zeta(3) - 1)) < mp.mpf("1e-30"),
        "sum_{n>=1} 1/(2n+1)^3 = 7zeta(3)/8 - 1")
    chk(abs(s2 - (mp.pi ** 2 / 8 - 1)) < mp.mpf("1e-30"),
        "sum_{n>=1} 1/(2n+1)^2 = pi^2/8 - 1")


# ---------------------------------------------------------------------------
# Stage 3 : Lemma 4, the elementary tail inequalities
# ---------------------------------------------------------------------------

def stage3():
    hdr("Stage 3 -- Lemma 4: kernel comparison (4a), midpoint rule (4b), T(N1) (4c)")
    bad = 0
    for N1 in range(1, 61):
        for m in range(0, N1):
            for n in range(N1, N1 + 200):
                if (n - m) * (N1 + 1) < (n + 1) * (N1 - m):
                    bad += 1
    chk(bad == 0, "(4a) (n-m)(N1+1) >= (n+1)(N1-m) for 0<=m<N1<=n (N1<=60, n<N1+200)",
        f"violations={bad}")
    # (4a) symbolic identity  (n-m)(N1+1) - (n+1)(N1-m) = (n-N1)(1+m)
    import sympy as sp
    nn, mm, LL = sp.symbols("n m L", integer=True)
    chk(sp.simplify((nn - mm) * (LL + 1) - (nn + 1) * (LL - mm) - (nn - LL) * (1 + mm)) == 0,
        "(4a) sympy: (n-m)(L+1)-(n+1)(L-m) = (n-L)(1+m)")
    # (4b) midpoint rule, on the three convex functions actually used
    for lab, f in (("x^-2", lambda x: 1 / x ** 2),
                   ("x^-3", lambda x: 1 / x ** 3),
                   ("h(x)=1/((x+1)^4(2x+1))", lambda x: 1 / ((x + 1) ** 4 * (2 * x + 1)))):
        bad = 0
        for a in (1, 2, 3, 5, 8, 13, 21, 50):
            if f(mp.mpf(a)) > mp.quad(f, [a - mp.mpf(1) / 2, a + mp.mpf(1) / 2]):
                bad += 1
        chk(bad == 0, f"(4b) midpoint rule holds for {lab}", f"violations={bad}")
    # (4b) log-convexity of h:  (log h)'' = 4/(x+1)^2 + 4/(2x+1)^2 > 0
    xs = sp.symbols("x", positive=True)
    hh = 1 / ((xs + 1) ** 4 * (2 * xs + 1))
    chk(sp.simplify(sp.diff(sp.log(hh), xs, 2)
                    - (4 / (xs + 1) ** 2 + 4 / (2 * xs + 1) ** 2)) == 0,
        "(4b) sympy: (log h)'' = 4/(x+1)^2 + 4/(2x+1)^2  (> 0 => h convex)")
    # (4c)  T(N1) = sum_{n>=N1} 1/((n+1)^4 (2n+1)) <= 1/(N1 (2N1+1)^3)
    KBIG = 2_000_000
    kk = np.arange(1, KBIG + 1, dtype=np.float64)
    term = 1.0 / (kk ** 4 * (2 * kk - 1.0))          # k = n+1
    csum = np.concatenate(([0.0], np.cumsum(term)))
    total = csum[-1] + 1.0 / (4 * KBIG ** 4)          # crude but > remainder
    bad = 0
    worst = 0.0
    for N1 in range(1, 401):
        T = total - csum[N1]                          # sum_{k>=N1+1}
        ub = 1.0 / (N1 * (2.0 * N1 + 1) ** 3)
        if T > ub:
            bad += 1
        worst = max(worst, T / ub)
    chk(bad == 0, "(4c) T(N1) <= 1/(N1(2N1+1)^3) for 1<=N1<=400", f"violations={bad}")
    print(f"    (4c) tightness: max T(N1)/bound = {worst:.8f}  (sup is 1, attained as N1->inf)")
    # mpmath spot check of (4c), including the N1 values actually shipped
    for N1 in (1, 3, 7, 12, 16, 40, 80, 500, 1000, 5000):
        T = mp.nsum(lambda n: mp.mpf(1) / ((n + 1) ** 4 * (2 * n + 1)), [N1, mp.inf])
        chk(T <= mp.mpf(1) / (N1 * (2 * N1 + 1) ** 3),
            f"(4c) mpmath spot check N1={N1}", f"{mp.nstr(T,10)}")


# ---------------------------------------------------------------------------
# Stage 4/5 : the certificate, exact rational arithmetic
# ---------------------------------------------------------------------------

_H2_CACHE = {0: F(0)}


def H2(n):
    if n not in _H2_CACHE:
        _H2_CACHE[n] = H2(n - 1) + F(1, n * n)
    return _H2_CACHE[n]


def Kf(n, m):
    s = F(1, (n + m + 1) ** 2)
    if n != m:
        s += F(1, (n - m) ** 2)
    return s


def certificate(p, N1, beta=BETA):
    """Exact-rational Theorem (*).  p : list of Fractions (length M), N1 >= M integer."""
    M = len(p)
    assert N1 >= M
    nrm2 = sum((2 * n + 1) * p[n] ** 2 for n in range(M))
    quad = sum(p[n] * p[m] * Kf(n, m) for n in range(M) for m in range(M)) \
        - 2 * sum(H2(n) * p[n] ** 2 for n in range(M))
    rho = quad / nrm2

    def g(n):
        val = sum(p[m] * Kf(n, m) for m in range(M))
        if n < M:
            val -= 2 * H2(n) * p[n]
        return val

    gs = {n: g(n) for n in range(N1)}
    res = sum((gs[n] - rho * (2 * n + 1) * p[n]) ** 2 / F(2 * n + 1) for n in range(M))
    res += sum(gs[n] ** 2 / F(2 * n + 1) for n in range(M, N1))
    W = sum(p[m] * (F((N1 + 1) ** 2, (N1 - m) ** 2) + 1) for m in range(M))
    T = F(1, N1 * (2 * N1 + 1) ** 3)            # Lemma 4c
    R = res + W * W * T
    r2 = R / nrm2
    assert rho > beta, "Temple needs rho > beta"
    hi = rho + r2 / (rho - beta)
    return dict(M=M, N1=N1, p=p, nrm2=nrm2, rho=rho, r2=r2, W=W, T=T,
                res_exact=res, tail=W * W * T, lo=rho, hi=hi, gs=gs)


def fdown(x, d=12):
    """decimal string rounding a Fraction DOWN to d places (safe lower bound)"""
    s = 10 ** d
    k = (x.numerator * s) // x.denominator
    return f"{k // s}.{abs(k) % s:0{d}d}"


def fup(x, d=12):
    s = 10 ** d
    k = -((-x.numerator * s) // x.denominator)
    return f"{k // s}.{abs(k) % s:0{d}d}"


def pvec_float(M, H):
    n = np.arange(M)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    off = np.zeros((M, M))
    nz = D > 0
    off[nz] = 1.0 / D[nz].astype(float) ** 2
    A = (off + 1.0 / S.astype(float) ** 2) * u[:, None] * u[None, :]
    A[n, n] -= 2 * u ** 2 * H[:M]
    w, V = np.linalg.eigh(A)
    v = V[:, -1] * np.sign(V[0, -1])
    p = v * u
    return p / p[0]


def pint(M, H, scale):
    """Positive INTEGER trial vector: nearest integers to scale * (top eigenvector of A_M)."""
    pf = pvec_float(M, H)
    p = [F(int(round(float(x) * scale))) for x in pf]
    assert all(x > 0 for x in p), "scale too small, zero entry"
    return p


def report(cert, name, big_A=None):
    lo, hi, rho, r2 = cert["lo"], cert["hi"], cert["rho"], cert["r2"]
    wid = hi - lo
    print(f"\n  [{name}]  M={cert['M']}  N1={cert['N1']}")
    print(f"     p = " + ", ".join(str(x) for x in cert["p"][:min(6, cert['M'])])
          + (" ..." if cert["M"] > 6 else ""))
    print(f"     rho  = {fdown(rho,15)}   (exact rational, numerator {len(str(rho.numerator))} digits)")
    print(f"     r^2  = {fup(r2,15)}    (exact head {fup(cert['res_exact']/cert['nrm2'],15)}"
          f" + analytic tail {fup(cert['tail']/cert['nrm2'],15)})")
    print(f"     g_lo = {fdown(lo,15)}")
    print(f"     g_hi = {fup(hi,15)}")
    print(f"     width = {float(wid):.6e}   relative = {float(wid/lo):.3e}")
    # the real gate: the analytic bracket must contain the true lambda_1.
    chk(lo <= G2_REF_F <= hi, f"{name}: encloses the reference value g(2)=1.3181402121321")
    chk(lo <= G2_HI_REF_F and hi >= G2_LO_REF_F,
        f"{name}: intersects the certified reference interval")
    if float(wid) > 2 * float(G2_HI_REF_F - G2_LO_REF_F):
        chk(lo <= G2_LO_REF_F and hi >= G2_HI_REF_F,
            f"{name}: contains the certified reference interval")
    else:
        print("     (certificate is comparable to / sharper than the reference interval;"
              " containment not required)")
    if big_A is not None:
        for N, lam in big_A.items():
            chk(float(lo) <= lam[0] <= float(hi), f"{name}: encloses lambda_1(A_N), N={N}",
                f"{lam[0]:.12f}")
    return cert


# integer trial vectors printed literally in the write-up
P_CERT_A = [F(200), F(45), F(13)]
P_CERT_A2 = [F(1000), F(228), F(69), F(25)]
P_CERT_B = [F(x) for x in (10000, 2286, 704, 275, 129, 70, 41, 27, 18, 13, 9, 7)]

# CERT-C / CERT-D: frozen integer literals.  They were produced once by
#   pint(30, H, 10**8)  and  pint(200, H, 10**14),
# i.e. nearest integers to scale * (top eigenvector of the M x M truncation, float64).
# They are FROZEN here so that the certificate -- and the exact rationals quoted in
# P4_g2_analytic.md -- do not depend on bit-for-bit LAPACK reproducibility.  Stage 5
# still re-derives them and checks agreement, but a mismatch is reported as a
# provenance note, not as a failure of the certificate (Theorem (*) holds for ANY
# positive rational p).
P_CERT_C = [F(x) for x in (
    100000000, 22857428, 7040522, 2750396, 1291917, 695581,
    413939, 265377, 180060, 127702, 93829, 70958,
    54962, 43439, 34928, 28506, 23568, 19709,
    16649, 14192, 12195, 10557, 9199, 8064,
    7108, 6297, 5603, 5005, 4482, 3995,
)]

P_CERT_D = [F(x) for x in (
    100000000000000, 22857432119997, 7040527826434, 2750401415493, 1291922387591,
    695585481598, 413942558793, 265380906783, 180063539928, 127704891952,
    93831771449, 70960718652, 54964469410, 43441671888, 34930973580,
    28508472364, 23570593365, 19711655708, 16652118172, 14195082013,
    12199126227, 10560888104, 9203625820, 8069547534, 7114566882,
    6304654320, 5613258034, 5019454146, 4506602113, 4061355221,
    3672924019, 3332522214, 3032945681, 2768249668, 2533499151,
    2324574197, 2138017071, 1970911264, 1820785148, 1685534745,
    1563361451, 1452721539, 1352284980, 1260901699, 1177573794,
    1101432542, 1031719299, 967769561, 908999606, 854895263,
    805002431, 758919044, 716288246, 676792565, 640148937,
    606104428, 574432556, 544930124, 517414468, 491721077,
    467701525, 445221663, 424160041, 404406530, 385861108,
    368432793, 352038704, 336603234, 322057312, 308337758,
    295386699, 283151062, 271582112, 260635044, 250268619,
    240444840, 231128651, 222287683, 213892011, 205913945,
    198327836, 191109901, 184238069, 177691840, 171452151,
    165501265, 159822662, 154400944, 149221745, 144271654,
    139538139, 135009482, 130674718, 126523580, 122546444,
    118734288, 115078644, 111571560, 108205566, 104973635,
    101869158, 98885914, 96018042, 93260016, 90606627,
    88052959, 85594371, 83226478, 80945135, 78746423,
    76626636, 74582262, 72609977, 70706631, 68869239,
    67094966, 65381127, 63725170, 62124670, 60577326,
    59080949, 57633458, 56232871, 54877305, 53564965,
    52294141, 51063204, 49870602, 48714853, 47594545,
    46508330, 45454920, 44433086, 43441654, 42479502,
    41545555, 40638789, 39758221, 38902912, 38071963,
    37264512, 36479735, 35716842, 34975073, 34253705,
    33552039, 32869408, 32205170, 31558710, 30929438,
    30316785, 29720207, 29139179, 28573199, 28021781,
    27484462, 26960792, 26450342, 25952695, 25467455,
    24994235, 24532667, 24082393, 23643070, 23214368,
    22795966, 22387557, 21988845, 21599543, 21219375,
    20848074, 20485384, 20131055, 19784847, 19446530,
    19115879, 18792677, 18476716, 18167793, 17865712,
    17570285, 17281327, 16998661, 16722115, 16451523,
    16186721, 15927554, 15673866, 15425509, 15182335,
    14944197, 14710951, 14482445, 14258521, 14038997,
    13823642, 13612096, 13403624, 13195949, 12974371,
)]


def stage45(eig, full=True):
    hdr("Stage 4 -- CERT-A/A' (hand-checkable) and CERT-B (headline), exact rationals")
    _, _, H = build_A(4000)

    certA = report(certificate(P_CERT_A, 12), "CERT-A", eig)
    certA2 = report(certificate(P_CERT_A2, 16), "CERT-A'", eig)
    certB = report(certificate(P_CERT_B, 40), "CERT-B", eig)
    certs = {"CERT-A": certA, "CERT-A'": certA2, "CERT-B": certB}

    # independent reconstruction of the shipped integer vectors from A_M's top
    # eigenvector, so the "where did p come from" step is itself regression-tested.
    # CERT-A and CERT-B are the raw roundings; CERT-A' was additionally tuned by a
    # +-2 integer search (this changes nothing about the validity of (*), only the
    # width) -- the check records exactly that.
    chk(pint(3, H, 200) == P_CERT_A, "CERT-A  p = round(200 * v_1(A_3))")
    chk(pint(12, H, 10000) == P_CERT_B, "CERT-B  p = round(10^4 * v_1(A_12))")
    raw4 = pint(4, H, 1000)
    chk(all(abs(a - b) <= 2 for a, b in zip(raw4, P_CERT_A2)),
        "CERT-A' p within +-2 of round(1000 * v_1(A_4))",
        f"raw={[str(x) for x in raw4]}")

    if full:
        hdr("Stage 5 -- CERT-C / CERT-D and the M^{-4} scaling table")
        certs["CERT-C"] = report(certificate(P_CERT_C, 80), "CERT-C", eig)
        certs["CERT-D"] = report(certificate(P_CERT_D, 500), "CERT-D", eig)
        # provenance of the frozen literals: re-derive and compare.  A mismatch in
        # the last digits only means this machine's LAPACK differs; it does NOT
        # invalidate the certificate, so it is reported, not asserted.
        for nm, Pfrozen, M, sc in (("CERT-C", P_CERT_C, 30, 10 ** 8),
                                   ("CERT-D", P_CERT_D, 200, 10 ** 14)):
            raw = pint(M, H, sc)
            worst = max(abs(float((a - b) / b)) for a, b in zip(raw, Pfrozen))
            print(f"    provenance {nm}: max relative deviation of frozen literal "
                  f"from round({sc} * v_1(A_{M})) on this machine = {worst:.2e}")
            chk(worst < 1e-6, f"{nm} frozen p reproduces round(scale * v_1(A_M))",
                "(informational: certificate is valid for ANY positive rational p)")
        print("\n    width scaling of the certificate (M -> relative width):")
        prev = None
        for M, N1, sc in ((3, 12, 10 ** 3), (4, 16, 10 ** 4), (6, 20, 10 ** 5),
                          (8, 24, 10 ** 6), (12, 40, 10 ** 6), (16, 40, 10 ** 6),
                          (20, 50, 10 ** 7), (30, 80, 10 ** 8), (40, 100, 10 ** 9),
                          (60, 150, 10 ** 10), (80, 200, 10 ** 11), (120, 300, 10 ** 12),
                          (160, 400, 10 ** 13), (200, 500, 10 ** 14)):
            c = certificate(pint(M, H, sc), N1)
            rel = float((c["hi"] - c["lo"]) / c["lo"])
            ratio = f"{prev/rel:6.2f}" if prev else "   ---"
            print(f"      M={M:4d} N1={N1:4d}   rel={rel:.4e}   gain vs previous: {ratio}")
            chk(c["lo"] <= G2_REF_F <= c["hi"], f"scaling table M={M} brackets g(2)")
            prev = rel
    return certs


# ---------------------------------------------------------------------------
# Stage 6 : regression of every certificate quantity against the matrix
# ---------------------------------------------------------------------------

def stage6(certs, eig):
    hdr("Stage 6 -- regression: exact certificate quantities vs the N x N matrix")
    for N in (400, 800, 1600):
        A, u, H = build_A(N)
        lam1N = eig[N][0]
        nn = np.arange(N)
        for name, c in certs.items():
            M, N1 = c["M"], c["N1"]
            if N1 + 50 > N:        # need a genuine tail beyond N1 inside the matrix
                continue
            p = np.array([float(x) for x in c["p"]])
            v = np.zeros(N)
            v[:M] = p * np.sqrt(2 * np.arange(M) + 1.0)
            nrm2 = v @ v
            chk(abs(nrm2 - float(c["nrm2"])) < 1e-12 * nrm2,
                f"{name} ||v||^2 vs matrix (N={N})")
            rho_m = (v @ A @ v) / nrm2
            chk(abs(rho_m - float(c["rho"])) < 1e-12 * abs(rho_m),
                f"{name} rho vs v^T A v/||v||^2 (N={N})",
                f"{rho_m:.15f} vs {float(c['rho']):.15f}")
            Av = A @ v
            gm = Av / u                      # (Av)_n = u_n g_n
            probe = sorted(set([0, 1, 2, M // 2, M - 1, M, M + 1, M + 3,
                                (M + N1) // 2, N1 - 1]))
            for n in probe:
                if not (0 <= n < N1):
                    continue
                chk(abs(gm[n] - float(c["gs"][n])) < 1e-11 * max(1.0, abs(gm[n])),
                    f"{name} g_{n} vs (Av)_n/u_n (N={N})",
                    f"{gm[n]:.12f} vs {float(c['gs'][n]):.12f}")
            resid = Av - rho_m * v
            r2_m = (resid @ resid) / nrm2
            chk(r2_m <= float(c["r2"]) * (1 + 1e-12),
                f"{name} certified r^2 >= true residual (N={N})",
                f"cert={float(c['r2']):.6e} true={r2_m:.6e}")
            env = float(c["W"]) / (nn[N1:] + 1.0) ** 2
            chk(bool(np.all(gm[N1:] > 0)), f"{name} g_n > 0 for n >= N1 (N={N})")
            chk(bool(np.all(gm[N1:] <= env * (1 + 1e-12))),
                f"{name} g_n <= W/(n+1)^2 for N1<=n<{N} (N={N})",
                f"max ratio={np.max(gm[N1:]/env):.6f}")
            tail_true = float(np.sum(gm[N1:] ** 2 / (2 * nn[N1:] + 1.0)))
            chk(tail_true <= float(c["tail"]) * (1 + 1e-12),
                f"{name} analytic tail >= true tail (N={N})",
                f"{float(c['tail']):.4e} vs {tail_true:.4e}")
            chk(float(c["lo"]) <= lam1N <= float(c["hi"]),
                f"{name} brackets lambda_1(A_N) (N={N})")


# ---------------------------------------------------------------------------
# Stage 7 : Collatz-Wielandt cross-check (numerical, NOT part of the proof)
# ---------------------------------------------------------------------------

def stage7():
    hdr("Stage 7 -- Collatz-Wielandt / weighted-Schur cross-check (numerical only)")
    NM = 20000
    n = np.arange(NM)
    H = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, NM) ** 2)))

    def cw(p, Ncheck):
        out = np.empty(Ncheck)
        for k in range(Ncheck):
            d = k - n
            t = np.where(d != 0, 1.0 / np.where(d == 0, 1, d) ** 2, 0.0)
            out[k] = np.dot(p, t + 1.0 / (k + n + 1.0) ** 2) - 2 * H[k] * p[k]
        return out / ((2 * np.arange(Ncheck) + 1.0) * p[:Ncheck])

    for name, p in (("p_n = 1/(n+1)^2", 1.0 / (n + 1.0) ** 2),
                    ("p_n = 1/((n+1)(n+2)(n+3))", 1.0 / ((n + 1.0) * (n + 2.0) * (n + 3.0)))):
        r = cw(p, 200)
        s = r.max()
        print(f"    {name:30s}  sup_n (Aw)_n/w_n <= {s:.7f}  (arg n={r.argmax()})")
        chk(s >= float(G2_REF), f"CW bound is an upper bound ({name})", f"{s:.7f}")
        chk(s < 1.45, f"CW bound < 1.45 ({name})", f"{s:.7f}")
    # head (top eigenvector of A_M) + continuity-matched c/(n+1)^3 tail
    _, _, Hf = build_A(4000)
    print("    head = top eigenvector of A_M, tail = c/(n+1)^3 matched at n=M-1:")
    for M in (12, 40, 100, 200, 400, 800):
        pf = pvec_float(M, Hf)
        p = np.zeros(NM)
        p[:M] = pf
        p[M:] = pf[M - 1] * float(M) ** 3 / (n[M:] + 1.0) ** 3
        s = cw(p, min(2 * M, 900)).max()
        print(f"      M={M:4d}   sup_n (Aw)_n/w_n = {s:.9f}   excess {s-float(G2_REF):.3e}")
        chk(s >= float(G2_REF), f"CW head+tail is an upper bound (M={M})")
    print("    With this (simple, explicitly analysable) weight family the CW excess decays")
    print("    like ~0.8/M -- the sup is attained at the matching point n=M-1.  Even at")
    print("    M=800 it is 1.0e-03, i.e. 20x WORSE than Temple at M=12 (CERT-B, 5.2e-05),")
    print("    and 10^3 x worse than Temple at M=30 (CERT-C, 1.2e-06).")
    print("    -> CW is viable but FIRST order in the trial-vector error and would need")
    print("       uniform-in-n analytic bounds on THREE infinite sums (the H_n^(2) block,")
    print("       the |n-m|^-2 block and the (n+m+1)^-2 block) for every n, whereas Temple")
    print("       is SECOND order and needs ONE tail bound (Lemma 4c).  Hence the demotion;")
    print("       see the route-comparison section of P4_g2_analytic.md.")


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "full"
    eig = stage0()
    stage1()
    stage2()
    stage3()
    certs = stage45(eig, full=(mode != "quick"))
    stage6(certs, eig)
    if mode != "quick":
        stage7()
    hdr(f"TOTAL: {NCHECK} checks, {NFAIL} failures, {time.time()-_T0:.1f} s")
    best = min(certs.values(), key=lambda c: c["hi"] - c["lo"])
    hb = certs["CERT-B"]
    print(f"\n  headline (hand-tabulated, M={hb['M']}, N1={hb['N1']}):")
    print(f"      {fdown(hb['lo'],13)}  <=  g(2)  <=  {fup(hb['hi'],13)}")
    print(f"      relative width {float((hb['hi']-hb['lo'])/hb['lo']):.3e}")
    print(f"\n  sharpest shipped (M={best['M']}, N1={best['N1']}):")
    print(f"      {fdown(best['lo'],15)}  <=  g(2)  <=  {fup(best['hi'],15)}")
    print(f"      relative width {float((best['hi']-best['lo'])/best['lo']):.3e}")
    print(f"\n  physics: (1/2pi) sqrt(g(2)) in "
          f"[{float(mp.sqrt(float(best['lo']))/(2*mp.pi)):.9f}, "
          f"{float(mp.sqrt(float(best['hi']))/(2*mp.pi)):.9f}]"
          f"   (Allen-Dynes 1975 strong-coupling constant 0.1827)")
    return 0 if NFAIL == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
