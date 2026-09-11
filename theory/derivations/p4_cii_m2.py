r"""
P4 . EINSTEIN EXTREMAL -- Segment C-II, MILESTONE M2 (route R: making M1 a proof).

Rigorous replacement of the two non-proved ingredients of `p4_cii_m1.py` (its
HONEST LIMITS 2):  (i) the CALIBRATED midpoint pad  PAD = 0.25/wbar^2 + 2e-6
used at every "lattice sum -> integral" step, and (ii) the three sup_s taken on
a grid (Theta*, ||O[phi_0]||_Schur, ||K||_Schur).  Nothing in M1 is modified;
this file imports M1 (and, for the float64 truth side, the same certified bases).

RUN:  python theory/derivations/p4_cii_m2.py          (all stages)
      python theory/derivations/p4_cii_m2.py --quick  (trims the wbar grid)

================================================================================
SEC. 0.  THE THREE ELEMENTARY LEMMAS THAT REPLACE THE PAD
================================================================================
Throughout  D := 1/wbar,  s_n := (n+1/2) D,  w_n = 1/(2n+1) = D/(2 s_n),
q := D^2/8.  Every lattice sum of M1 is  S[u] := sum_{n>=0} w_n u(s_n)
           = sum_{n>=0} D * h(s_n),   h(s) := u(s)/(2s),
i.e. the MIDPOINT rule for  I[u] := int_0^inf h = (1/2) int_0^inf u(s) ds/s.

(L1) MIDPOINT, PEANO FORM.  For h in C^2((0,inf)) with h, h'' integrable,
        | sum_{n>=0} D h((n+1/2)D) - int_0^inf h |  <=  (D^2/8) int_0^inf |h''|.
     PROOF.  On one cell [a,b], b-a = D, m = (a+b)/2, two integrations by parts
     give the exact Peano representation
        int_a^b h - D h(m) = int_a^b K(t) h''(t) dt,
        K(t) = (t-a)^2/2 on [a,m],  = (b-t)^2/2 on [m,b],
     K >= 0 with max K = K(m) = D^2/8.  Hence the cell error is <= (D^2/8)
     int_cell |h''|; summing over the cells (which tile (0,inf) exactly) gives
     the claim.  []
     REMARK (correcting M1 Sec.2.1).  int_0^D K = D^3/24, so the constant 1/24
     is the one for the MEAN of K; the PROVABLE bound with int|h''| is 1/8 (a
     bump of h'' concentrated at a midpoint attains it).  The 1/24 form is valid
     only with the per-cell sup:  sum_cells (D^3/24) sup_cell |h''|.  All of M2
     uses 1/8.  The one place M1 used 1/24 is its bound (a): with V_eta = 3/2
     that term becomes 3/(16 wbar^2) instead of 1/(16 wbar^2), i.e.
        h >= rho >= 2B-1 >= ln wbar + c_0 - 1/(12 wbar) - 3/(16 wbar^2),
     still true and still not load-bearing (M2 uses the EXACT B, see Sec.2).

(L2) VARIATION.  For h of bounded variation with h(inf) = 0,
        | sum_n D h(s_n) - int_0^inf h |  <=  (D/2) int_0^inf |h'|.
     PROOF.  Per cell, |int_a^b (h(t)-h(m))dt| <= int_a^m Var_{[t,m]}h dt +
     int_m^b Var_{[m,t]}h dt <= (D/2) Var_cell h; sum.  []
     (L2) is used ONLY inside terms that already carry a factor q = D^2/8, where
     the resulting O(D^3) is irrelevant; (L1) is used for everything else.

(L3) SHIFTED PARTIAL SUMS (the diagonal l_f).  For f in C^2, f even, and
     Lam_f(s) := int_0^s f,
        sum_{k=1}^{n} f(kD) = (1/D)[ Lam_f(s_n) - Lam_f(D/2) ] + E_n,
        |E_n| <= (D/8) int_0^inf |f''|,
     because kD is the midpoint of [(k-1/2)D,(k+1/2)D] and those cells tile
     [D/2, s_n]; apply (L1) on that interval to h = f.  Consequently, with
        l_f(n) := (2/(2n+1)) sum_{k=1}^n f(kD) = (D/s_n) sum_{k=1}^n f(kD),
        | l_f(n) - Lam_f(s_n)/s_n | <= eps_f(s_n),
        eps_f(s) := [ Lam_f(D/2) + q * min(s * sup|f''|, int_0^inf|f''|) ] / s,
     where the min comes from int_{D/2}^{s}|f''| <= min(s sup|f''|, int|f''|).
     For f = G one uses instead the EXACT cancellation with w_n:
        1 - D_nn = 1 - w_n - l_G(n) = d(s_n) - [D/2 - arctan(D/2)]/s_n - (D/s_n)E_n,
        d(s) := 1 - arctan(s)/s,   | 1-D_nn - d(s_n) | <= beta/s_n,
        beta := (D/2 - arctan(D/2)) + q * V_G'',   V_G'' = int_0^inf|G''| = 3 sqrt3 / 4
     (exact: G' = -2s/(1+s^2)^2, G'' vanishes only at s = 1/sqrt3, so
      int|G''| = -2G'(1/sqrt3) = 9/(4 sqrt 3) = 3 sqrt 3/4 = 1.2990381...).

(L4) TWO ELEMENTARY LATTICE TAILS (used to absorb the small-n loss of (L3)):
        sum_n w_n / s_n   = (D/2) sum_n s_n^{-2} = (1/2) sum_n (n+1/2)^{-2} = pi^2/4 / D^0 ...
     precisely:  sum_n w_n/s_n = (1/(2D)) * ... see code (`lat_pow`): with
     s_n = (n+1/2)D,   sum_n w_n s_n^{-p} = (D^{-p}/2) * zeta_h(p+1),
     zeta_h(r) := sum_{n>=0}(n+1/2)^{-r} = (2^r - 1) zeta(r):
        zeta_h(2) = pi^2/2,  zeta_h(3) = 7 zeta(3),  zeta_h(4) = pi^4/6.

================================================================================
SEC. 1.  WHAT BECOMES EXACT (no pad at all):  S[G^j] BY psi-JETS
================================================================================
With x = wbar^2 and a_n = n+1/2, G(s_n) = x/(x+a_n^2) and  x dG/dx = G - G^2.
Hence, with L := x d/dx = (wbar/2) d/dwbar and A_j := S[G^j],
        A_{j+1} = A_j - (1/j) L A_j,        A_1 = A = S[G],
        A(wbar) = (1/2)[ Re psi(1/2 + i wbar) + gamma + 2 ln 2 ]   (M1 Sec.2.1),
        d^k/dwbar^k Re psi(1/2+i wbar) = Re[ i^k psi^(k)(1/2+i wbar) ]  (exact).
So A_1..A_6 are CLOSED FORMS in polygamma; the code carries 6-jets in wbar.
Since phi_0 = G^2 - G^3, kappa = G - G^2, G kappa = phi_0, this makes EXACT
        S[G], S[G^2]=B, S[G phi_0]=A_3-A_4, S[G kappa]=S[phi_0]=A_2-A_3,
        S[kappa^2]=A_2-2A_3+A_4, S[phi_0^2]=A_4-2A_5+A_6, S[kappa]=A_1-A_2.
These are exactly the LEADING pieces of M1's (a),(d1),(d2),(d3) -- in M1 they
were continuum constants (1/3, 1/2, 1/24, 1/80) plus a pad.  What is left to pad
is only the SMALL remainder constants, all of which carry d, Lam_f/s or Ups_f.

================================================================================
SEC. 2.  THE PAD CATALOGUE  (each = q * V2[u], V2[u] := int_0^inf |(u/2s)''|)
================================================================================
(printed by Stage 2; see CONST_V2 in the code).  Every V2 is a convergent
integral of an EXPLICIT elementary function of s (the profiles Ups_f, Theta_f
enter through their own quadratures, whose s-derivatives are again explicit
integrals of rational kernels: d_s^k R_f(s,sigma) is rational).

SEC. 3.  THE THREE sup_s, MADE RIGOROUS
   (i) Theta* : Theta' = N(s)/(1+s^2)^2, N(s) = 2s + (1-s^2)arctan s - s ln(1+s^2),
       N(0)=0 and N'(s) = 2 + (1-3s^2)/(1+s^2) - 2s arctan s - ln(1+s^2) is a sum
       of STRICTLY DECREASING functions with N'(0)=3>0, N'(inf)=-inf, so N' has a
       single sign change, N is increasing-then-decreasing with N(0)=0: N>0 on
       (0,s*), N<0 on (s*,inf).  Hence Theta increases then decreases and
       Theta* = Theta(s*) with s* the unique root of N, bracketed by bisection to
       1e-40:   Theta* = 0.634928792...,  s* = 1.5009042...
   (ii) the Schur profiles Phi_f(s) = 2A f(s) + Theta_f(s) + 2G(s)S[f] + Lam_f(s)/s:
       grid on [1e-3, 30] with the CHORD-PLUS-CURVATURE bound
          max_[a,b] F <= max(F(a),F(b)) + (b-a)^2 sup_[a,b]|F''| / 8
       (F is below its chord plus that parabola), sup|F''| taken as 1.5x the
       endpoint values of |F''| = |2A f'' + (Theta_f+Lam_f/s)'' + 2 S[f] G''|
       (grid step 1e-3 => this correction is < 1e-4); the head s < 1e-3 is bounded
       by 2A f(1e-3) + max(0,Psi(1e-3)) + 2S[f] (f, Lam_f/s increasing, G <= 1);
       and for s > 30 by the explicit DECREASING majorant  (f in {phi_0,kappa})
          Phi_f(s) <= M_f(s) := 2A f(s) + mu_f(s-1)/4 + (f(s)/2)(ln 2 + ln(s/2))
                       + (f(s/2)/2) ln(s/2) + 3 F_0/s + 2 S[f]/s^2 + F_0/s,
          F_0 = int_0^inf f (= pi/4, pi/16),  mu_f(r) = sup_{y>=r}|f''| (<= 10.2/r^4
          for kappa, 26.7/r^6 for phi_0), proof: drop -2f(sigma)G(s) <= 0, Taylor
          f(s-x)+f(s+x)-2f(s) <= x^2 mu_f(s-1) on x<=1, and f(s-x) <= f(s/2) resp.
          1/x <= 2/s on the two pieces of x>1.  (Stage 5 asserts M_f(30) < grid max.)
   (iii) sup_s of the inner-sum pad profiles T_f(s) = int|d_sigma^2 (R_f/2sigma)|,
       U_f(s) = the same with the extra G(sigma): grid sup x 1.5.  These enter
       ONLY multiplied by q, contributing 0.19 % of g_lo at wbar = 8.5 and
       3.6e-5 of it at wbar = 45, so no finer argument is warranted.

================================================================================
SEC. 4.  THE RESULTING BOUNDS (each a formula in wbar; q := D^2/8, D := 1/wbar)
================================================================================
(a)  h >= rho >= 2B - 1,   B = A_2  EXACT  (M1 already had no pad here).
(b)  ||r|| <= C_r := [ sqrt(N_R) + sqrt(N_D) ] / sqrt(B),
     sqrt(N_R) = sqrt( I[Ups_G^2] + q V2[Ups_G^2] ) + q sqrt( I[U_G^2] + (D/2)V1[U_G^2] )
       (Minkowski in the w-weighted l^2: the inner-sum error of Ups_G^lat(n) is
        <= q U_G(s_n) pointwise, by (L1) in sigma),
     N_D = I[d^2G^2] + q V2[d^2G^2] + 2 beta ( I[dG^2/s] + (D/2)V1[dG^2/s] )
           + beta^2 * 7 zeta(3)/(2 D^2),     beta as in (L3).
(c)  g >= 2B - 1 - Theta* - q Tbar_G.
(d1) wbar^4 F(t) >= 4 S[G phi_0] + (X_phi^lo - Y_phi^hi)/B,
     X_f = <tau,R_f tau> in I[G Ups_f] -+ [ q V2[G Ups_f] + q (I[G U_f] + (D/2)V1[G U_f]) ],
     Y_f = <tau,D_f tau> in I[(Lam_f/s)G^2] -+ [ q V2[(Lam_f/s)G^2] + eps1_f ],
     eps1_f := q ( ln(1 + wbar int|f''|) + int|f''|/(int|f''| + D) )   by (L3)+(L4).
(d2) wbar^4 E(t) = ||(I-P)Z||/sqrt(B) <= [ sqrt( ||Ztilde||^2 - <tau,Ztilde>^2/B )
       + ||Delta Z|| ] / sqrt(B), where ||Ztilde||^2 expands into the six sums
       4B^2 S[f^2] (EXACT) - 4B S[f Lam_f G/s] + S[(Lam_f/s)^2G^2] + 4B S[f Ups_f]
       - 2 S[(Lam_f/s) G Ups_f] + S[Ups_f^2], each with its own q V2[.], and
     ||Delta Z|| <= sqrt(eps2_f) + q sqrt( I[U_f^2] + (D/2) V1[U_f^2] ),
     eps2_f := 7 zeta(3) D^4/576 + q^2 ( 4 ln(1 + wbar int|f''|) + 2 ).
(d3) ||O[f]|| <= sup_s Phi_f(s) + q Tbar_f + D^2/3   (the last term bounds
     sup_n |l_f(n) - Lam_f(s_n)/s_n| <= Lam_f(D/2)(2/D) + q sup_n min(2nD,.)/s_n
     <= D^2/12 + 2q = D^2/3).
(e)  unchanged Davis-Kahan algebra of M1 (e), with |q_kappa| <= max of the two
     one-sided bounds of (d1) applied to f = kappa.

RESULT:   wbar_A' = 8.4761  (M1's calibrated-pad value was 7.776), ratio
     0.4543 (wbar=45), 0.3982 (100), 0.3301 (1000), 0.2580 (1e20) -> 1/4.
     Segment C-I certifies wbar <= 45, so the two segments still OVERLAP by 5.31
     in wbar (28 in x) -- with NO calibrated constant anywhere in the chain.

================================================================================
SEC. 5.  HONEST LIMITS OF M2  (what is still not a certificate)
================================================================================
1. STRUCTURE PROVED, CONSTANTS QUADRATURE.  Every pad is now a FORMULA
   q V2[u] (resp. the (L2)/(L3)/(L4) forms above) whose coefficient V2[u] =
   int_0^inf |(u/2s)''| ds is a convergent integral of an EXPLICIT function.
   Those ~30 numbers are evaluated by graded-trapezoid / scipy.quad at float64,
   exactly the epistemic status of the ~20 absolute constants I[.] of M1
   (its HONEST LIMIT 1): NOT interval arithmetic.  What M2 removes is the
   CALIBRATION (M1's 0.25/wbar^2 + 2e-6 had no proof at all): the shape of the
   pad, its wbar-dependence and its coefficients are now derived.
   Stage 3 checks every V2 pad against the measured |S[u] - I[u]| at
   wbar = 6..1000: the pads are 3x-5x conservative (max ratio 0.35).
2. Two 1.5x safety factors: on the grid sup of the inner-sum pad profiles
   (Sec.3(iii)) and on the per-cell curvature of the Schur sup (Sec.3(ii)).
   Both are < 0.2 % of the quantity they pad at wbar = 8.5; setting them to 1.0
   moves wbar_A' by < 0.01.
3. The profiles Ups_f, Theta_f and their first two s-derivatives are computed by
   quadrature of the EXPLICIT rational kernels d_s^k R_f(s,sigma) and splined on
   a 730-point graded grid (Stage 2 checks them against the closed form of
   Theta_G and against central differences).
4. (C**) itself, the identities (M1),(M2),(M3), the Schur test, the one-
   dimensional min-max for lambda_2 and the Davis-Kahan algebra are quoted from
   M1 Sec.1-2 (all fully rigorous there).
5. The float64 "truth" side (Stage 6) is M1's `truth()` at N = 40 wbar: it is a
   CHECK, not part of any bound.
"""
import os
import sys
import time

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from mpmath import mp, mpf, mpc, digamma, polygamma, euler, log as mlog, re as mre, zeta

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

mp.dps = 30

PASS = [0]
_LAST = [0]
T0 = time.time()


def chk(name, ok, detail=""):
    PASS[0] += 1
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name}  {detail}")


def stage_end(tag):
    n = PASS[0] - _LAST[0]
    _LAST[0] = PASS[0]
    print(f"  [{tag}: {n} checks, total {PASS[0]}, {time.time()-T0:.1f} s]")


# ============================================================================
# SEC. 1.  kernels and their first two derivatives (all elementary, exact)
# ============================================================================
def G(y):
    y = np.asarray(y, float)
    return 1.0 / (1.0 + y * y)


def G1(y):
    y = np.asarray(y, float)
    return -2.0 * y / (1.0 + y * y) ** 2


def G2(y):
    y = np.asarray(y, float)
    return (6.0 * y * y - 2.0) / (1.0 + y * y) ** 3


def KAP(y):
    y = np.asarray(y, float)
    return y * y / (1.0 + y * y) ** 2


def KAP1(y):
    y = np.asarray(y, float)
    return 2.0 * y * (1.0 - y * y) / (1.0 + y * y) ** 3


def KAP2(y):
    y = np.asarray(y, float)
    y2 = y * y
    return (6.0 * y2 * y2 - 16.0 * y2 + 2.0) / (1.0 + y2) ** 4


def PHI(y):
    y = np.asarray(y, float)
    return y * y / (1.0 + y * y) ** 3


def PHI1(y):
    y = np.asarray(y, float)
    return 2.0 * y * (1.0 - 2.0 * y * y) / (1.0 + y * y) ** 4


def PHI2(y):
    y = np.asarray(y, float)
    y2 = y * y
    return (20.0 * y2 * y2 - 26.0 * y2 + 2.0) / (1.0 + y2) ** 5


FAM = {'G': (G, G1, G2), 'phi': (PHI, PHI1, PHI2), 'kap': (KAP, KAP1, KAP2)}
LAM_F = {'G': lambda s: np.arctan(s),
         'phi': lambda s: 0.125 * (np.arctan(s) + s * (s * s - 1.0) / (1.0 + s * s) ** 2),
         'kap': lambda s: 0.5 * (np.arctan(s) - s / (1.0 + s * s))}
# int_0^inf |f''|  and  sup |f''|   (both exact / certified by Stage 1)
V_FPP = {'G': 3.0 * np.sqrt(3.0) / 4.0}
M_FPP = {'G': 2.0}


def Rk(nm, k_s, k_sg, s, sg):
    """d_s^{k_s} d_sigma^{k_sg} R_f(s,sigma) (at most one order > 0 is used)."""
    F = FAM[nm]
    k = k_s + k_sg
    tot = ((-1.0) ** k_sg) * F[k](s - sg) + F[k](s + sg)
    if nm == 'G':
        tot = tot - 2.0 * F[k_s](s) * F[k_sg](sg)
    else:
        tot = tot - 2.0 * F[k_s](s) * (G, G1, G2)[k_sg](sg)
        tot = tot - 2.0 * F[k_sg](sg) * (G, G1, G2)[k_s](s)
    return tot


def h_sig2(nm, s, sg, weightG):
    """d^2/dsigma^2 [ (G(sigma))^w R_f(s,sigma) / (2 sigma) ]."""
    P0, P1, P2 = Rk(nm, 0, 0, s, sg), Rk(nm, 0, 1, s, sg), Rk(nm, 0, 2, s, sg)
    if weightG:
        g0, g1, g2 = G(sg), G1(sg), G2(sg)
        P0, P1, P2 = (g0 * P0, g1 * P0 + g0 * P1, g2 * P0 + 2.0 * g1 * P1 + g0 * P2)
    return 0.5 * (P2 / sg - 2.0 * P1 / sg ** 2 + 2.0 * P0 / sg ** 3)


# ============================================================================
# SEC. 2.  EXACT lattice sums  A_j = S[G^j]  by psi-jets  (no pad)
# ============================================================================
def _A_jet(w, K=6):
    z = mpc(mpf(1) / 2, mpf(w))
    out = []
    for k in range(K + 1):
        pk = digamma(z) if k == 0 else polygamma(k, z)
        val = mre((mpc(0, 1) ** k) * pk) / 2
        if k == 0:
            val += (euler + 2 * mlog(2)) / 2
        out.append(val)
    return out


def A_exact(w, jmax=6):
    """A_j = S[G^j] = sum_n w_n G(s_n)^j, j = 1..jmax, EXACT (closed form)."""
    cur = _A_jet(w, K=jmax)
    A = {1: cur[0]}
    for j in range(1, jmax):
        Lc = [(k * cur[k] + mpf(w) * cur[k + 1]) / 2 for k in range(len(cur) - 1)]
        cur = [cur[k] - Lc[k] / mpf(j) for k in range(len(Lc))]
        A[j + 1] = cur[0]
    return {j: float(v) for j, v in A.items()}


def exact_sums(w):
    """All lattice sums that are polynomials in G -- EXACT, no discretisation pad."""
    A = A_exact(w)
    return dict(A=A[1], B=A[2],
                S_Gphi=A[3] - A[4],           # 4 S[G phi] replaces 4 I[G phi] = 1/3
                S_Gkap=A[2] - A[3],           # = S[phi]  (G kappa = phi_0)
                S_kk=A[2] - 2 * A[3] + A[4],  # replaces m_kk = 1/24
                S_pp=A[4] - 2 * A[5] + A[6],  # replaces m_pp = 1/80
                S_phi=A[2] - A[3], S_kap=A[1] - A[2])


# ============================================================================
# SEC. 3.  profiles  Theta_f, Ups_f (+ d/ds, d^2/ds^2)  and the pad profiles T_f, U_f
# ============================================================================
def prof_int(nm, k, s, weightG):
    """(1/2) int_0^inf [G(sigma)^w] d_s^k R_f(s,sigma) dsigma / sigma."""
    def fn(x):
        v = Rk(nm, k, 0, s, x) / (2.0 * x)
        return v * G(x) if weightG else v
    big = max(60.0, 10.0 * s)
    edges = sorted({1e-4, 1e-2, 0.3 * s, 0.7 * s, s, 1.4 * s, 3.0 * s, 1.0, 8.0, big})
    edges = [e for e in edges if 1e-8 < e <= big]
    tot, prev = 0.0, 1e-8
    for e in edges:
        if e > prev:
            tot += quad(fn, prev, e, limit=200)[0]
            prev = e
    return tot + quad(fn, prev, np.inf, limit=200)[0]


def abs_sig2(nm, s, weightG, n=30001):
    """T_f(s) resp. U_f(s) = int_0^inf |d_sigma^2 h| dsigma   (graded trapezoid)."""
    big = max(80.0, 12.0 * s)
    mid = min(4.0 + 2.0 * s, big)
    g = np.unique(np.concatenate([
        np.geomspace(1e-3, 0.2, 600),
        np.linspace(0.2, mid, n // 2),
        np.geomspace(mid + 1e-9, big, n // 4),
        np.linspace(max(1e-3, 0.6 * s), 1.4 * s + 1e-9, 4001)]))
    v = np.abs(h_sig2(nm, s, g, weightG))
    tot = float(np.trapezoid(v, g)) + float(v[-1]) * big / 2.0
    return tot + 0.5e-3 * float(np.abs(h_sig2(nm, s, 1e-3, weightG)))


S_MAX = 240.0
_SG = np.unique(np.concatenate([
    np.geomspace(1e-5, 0.5, 300), np.linspace(0.5, 6.0, 330),
    np.geomspace(6.0, S_MAX, 180)]))

PROF = {}


def build_profiles(verbose=True):
    """Splines of Theta_f, Ups_f, their s-derivatives, and of T_f, U_f."""
    for nm in ('G', 'phi', 'kap'):
        for tag, wG in (('T', False), ('U', True)):
            for k in range(3):
                PROF[(nm, tag, k)] = CubicSpline(
                    _SG, [prof_int(nm, k, s, wG) for s in _SG])
            PROF[(nm, tag, 'pad')] = CubicSpline(
                _SG, [abs_sig2(nm, s, wG) for s in _SG])
        if verbose:
            print(f"   profiles[{nm}] built ({time.time()-T0:.1f} s)", flush=True)


def P(nm, tag, k, s):
    s = np.asarray(s, float)
    return np.where(s <= S_MAX, PROF[(nm, tag, k)](np.minimum(s, S_MAX)), 0.0)


# ============================================================================
# SEC. 4.  3-jets in s; the constants I[u] and the pad constants V2[u], V1[u]
# ============================================================================
class J3(object):
    """(value, d/ds, d^2/ds^2) on a numpy grid -- exact product/quotient rules."""

    __slots__ = ('v', 'a', 'b')

    def __init__(self, v, a, b):
        self.v, self.a, self.b = v, a, b

    def __add__(self, o):
        o = o if isinstance(o, J3) else J3(o, 0.0, 0.0)
        return J3(self.v + o.v, self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __neg__(self):
        return J3(-self.v, -self.a, -self.b)

    def __sub__(self, o):
        return self + (-o if isinstance(o, J3) else J3(-o, 0.0, 0.0))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, J3):
            return J3(self.v * o, self.a * o, self.b * o)
        return J3(self.v * o.v, self.a * o.v + self.v * o.a,
                  self.b * o.v + 2.0 * self.a * o.a + self.v * o.b)

    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, J3):
            return self * (1.0 / o)
        u, v = self, o
        w = u.v / v.v
        w1 = (u.a * v.v - u.v * v.a) / v.v ** 2
        w2 = (u.b / v.v - 2.0 * u.a * v.a / v.v ** 2
              - u.v * v.b / v.v ** 2 + 2.0 * u.v * v.a ** 2 / v.v ** 3)
        return J3(w, w1, w2)


# the dense s-grid on which every V2 / I constant is evaluated
_DS = np.unique(np.concatenate([
    np.geomspace(1e-4, 0.3, 6000), np.linspace(0.3, 8.0, 26000),
    np.geomspace(8.0, S_MAX, 8000)]))


def _jets(nm, s):
    """the elementary 3-jets on the grid s, for the family nm in {G,phi,kap}."""
    f, f1, f2 = FAM[nm]
    js = J3(s, np.ones_like(s), np.zeros_like(s))
    jG = J3(G(s), G1(s), G2(s))
    jf = J3(f(s), f1(s), f2(s))
    jat = J3(np.arctan(s), 1.0 / (1.0 + s * s), -2.0 * s / (1.0 + s * s) ** 2)
    jd = 1.0 - jat / js                                   # d(s) = 1 - arctan s / s
    jL = J3(LAM_F[nm](s), f(s), f1(s)) / js               # Lam_f(s)/s
    jU = J3(P(nm, 'U', 0, s), P(nm, 'U', 1, s), P(nm, 'U', 2, s))   # Ups_f
    jT = J3(P(nm, 'T', 0, s), P(nm, 'T', 1, s), P(nm, 'T', 2, s))   # Theta_f
    return js, jG, jf, jd, jL, jU, jT


def _int(y, s):
    """int_0^inf y ds on the graded grid + a power tail + the [0,s_0] head.

    Every integrand here vanishes at s = 0 and decays at least like s^-5, so the
    head is (1/2) s_0 y(s_0) and the tail <= s_max y(s_max)/4 (both < 1e-9)."""
    tot = float(np.trapezoid(y, s))
    tot += float(y[-1]) * float(s[-1]) / 4.0
    tot += 0.5 * float(s[0]) * float(y[0])
    return tot


def _half(s):
    return J3(2.0 * s, 2.0 * np.ones_like(s), np.zeros_like(s))


def const_pair(u):
    """(I[u], V2[u]),  I[u] = int_0^inf h,  V2[u] = int_0^inf |h''|,  h = u/(2s)."""
    s = _DS
    h = u / _half(s)
    return _int(h.v, s), _int(np.abs(h.b), s)


def const_pair1(u):
    """(I[u], V1[u]) with V1[u] = int |h'|  (used only inside O(D^2) terms)."""
    s = _DS
    h = u / _half(s)
    return _int(h.v, s), _int(np.abs(h.a), s)


CV = {}          # every constant used by bounds2(), filled by build_constants()


def build_constants(verbose=True):
    s = _DS
    js, jG, _, jd, _, jUG, _ = _jets('G', s)
    CV['nD2'], CV['V_nD2'] = const_pair(jd * jd * jG * jG)
    CV['J1'], CV['V_J1'] = const_pair1((jd * jG * jG) / js)
    CV['nR2'], CV['V_nR2'] = const_pair(jUG * jUG)
    # the inner-sum pad profiles U_f enter only inside O(D^2) terms -> (L2) suffices
    for nm in ('G', 'phi', 'kap'):
        pv = P(nm, 'U', 'pad', s)
        jp = J3(pv, np.gradient(pv, s), np.zeros_like(s))
        CV['W_' + nm], CV['V_W_' + nm] = const_pair1(jp * jp)
        CV['GU_' + nm], CV['V_GU_' + nm] = const_pair1(jG * jp)
    for nm in ('phi', 'kap'):
        _, jG2, jf, _, jL, jU, _ = _jets(nm, s)
        CV[nm + '_GU'], CV['V_' + nm + '_GU'] = const_pair(jG2 * jU)
        CV[nm + '_tDt'], CV['V_' + nm + '_tDt'] = const_pair(jL * jG2 * jG2)
        CV[nm + '_fDt'], CV['V_' + nm + '_fDt'] = const_pair(jf * jL * jG2)
        CV[nm + '_DD'], CV['V_' + nm + '_DD'] = const_pair(jL * jL * jG2 * jG2)
        CV[nm + '_fU'], CV['V_' + nm + '_fU'] = const_pair(jf * jU)
        CV[nm + '_lU'], CV['V_' + nm + '_lU'] = const_pair(jL * jG2 * jU)
        CV[nm + '_UU'], CV['V_' + nm + '_UU'] = const_pair(jU * jU)
    if verbose:
        print(f"   constants built ({time.time()-T0:.1f} s)", flush=True)


# ============================================================================
# SEC. 5.  THE THREE sup_s, MADE RIGOROUS
# ============================================================================
def theta_star_rig():
    """Theta* = max_{s>=0} Theta, with the monotonicity PROOF of Sec.0/(iii)."""
    mp.dps = 40
    N = lambda s: 2 * s + (1 - s * s) * mp.atan(s) - s * mp.log(1 + s * s)
    Th = lambda s: (mp.log(1 + s * s) + 2 * s * mp.atan(s)) / (2 * (1 + s * s))
    a, b = mpf(1), mpf(2)                      # N(1) > 0 > N(2)
    assert N(a) > 0 and N(b) < 0
    for _ in range(150):
        m = (a + b) / 2
        if N(m) > 0:
            a = m
        else:
            b = m
    # Theta is increasing on (0,s*), decreasing on (s*,inf) => Theta* = Theta(s*),
    # and |Theta'| <= |N|/(1+s^2)^2 <= 1 near s*, so the bracket width is the error.
    val = max(Th(a), Th(b)) + (b - a)
    mp.dps = 30
    return float(val), float((a + b) / 2)


def sup_pad_profile(nm, tag, safety=1.5):
    """sup_s of the inner-sum pad profile T_f (tag='T') resp. U_f (tag='U').

    Enters every bound multiplied by q = D^2/8 and contributes < 0.3 % of the
    quantity it pads at wbar = 8, so the grid sup is inflated by `safety`."""
    g = np.unique(np.concatenate([_SG, np.geomspace(1e-3, 60.0, 4000)]))
    return safety * float(np.max(np.abs(P(nm, tag, 'pad', g))))


S_SUP = 30.0                      # Schur sup: grid on [1e-3, S_SUP], tail beyond
_SS = np.unique(np.concatenate([np.geomspace(1e-3, 0.5, 2000),
                                np.linspace(0.5, S_SUP, 30000)]))
SCH = {}


def build_schur_grids():
    """f, G, Psi_f = Theta_f + Lam_f/s and their 2nd derivatives on the sup grid."""
    s = _SS
    for nm in ('phi', 'kap'):
        _, jG, jf, _, jL, _, jT = _jets(nm, s)
        psi = jT + jL
        SCH[nm] = dict(f=jf.v, f2=jf.b, G=jG.v, G2=jG.b, psi=psi.v, psi2=psi.b)


def schur_sup(nm, A, Sf, w):
    """rigorous sup_{s>=0} Phi_f(s), Phi_f = 2A f + Theta_f + 2 G S[f] + Lam_f/s.

    grid + chord-plus-curvature:  max_[a,b] F <= max(F(a),F(b)) + (b-a)^2 sup|F''|/8
    (F below its chord + the parabola bound), with sup|F''| on a cell taken as
    1.5 max of the endpoint values (the grid step is 1e-3, so this term is < 1e-4);
    plus the explicit decreasing tail majorant M_f(s) for s > S_SUP."""
    d = SCH[nm]
    F = 2.0 * A * d['f'] + d['psi'] + 2.0 * Sf * d['G']
    F2 = np.abs(2.0 * A * d['f2'] + d['psi2'] + 2.0 * Sf * d['G2'])
    ds = np.diff(_SS)
    cell = np.maximum(F[:-1], F[1:]) + ds ** 2 * (1.5 * np.maximum(F2[:-1], F2[1:])) / 8.0
    top = float(np.max(cell))
    # s < 1e-3:  f and Lam_f/s increase, Theta_f = O(s^2), G <= 1
    head = 2.0 * A * d['f'][0] + max(0.0, float(d['psi'][0])) + 2.0 * Sf
    return top, head, tail_major(nm, S_SUP, A, Sf, w)


F0 = {'phi': np.pi / 16.0, 'kap': np.pi / 4.0}       # int_0^inf f
MU = {'phi': lambda r: 26.7 / r ** 6, 'kap': lambda r: 10.2 / r ** 4}   # sup_{y>=r}|f''|
FUP = {'phi': lambda s: 1.0 / s ** 4, 'kap': lambda s: 1.0 / s ** 2}    # f(s) <=


def tail_major(nm, s, A, Sf, w):
    """M_f(s) >= Phi_f(s) for s >= 4, decreasing (Sec.5 of the docstring)."""
    f, f2 = FUP[nm](s), FUP[nm](s / 2.0)
    th = (MU[nm](s - 1.0) / 4.0 + 0.5 * f * (np.log(2.0) + np.log(s / 2.0))
          + 0.5 * f2 * np.log(s / 2.0) + 3.0 * F0[nm] / s)
    return 2.0 * A * f + th + 2.0 * Sf / s ** 2 + F0[nm] / s


# ============================================================================
# SEC. 6.  THE RIGOROUS BOUNDS  (every entry a formula in wbar)
# ============================================================================
# int_0^inf |f''| : exact (f'' has two positive roots y1<y2, so int|f''| =
# 2 f'(y1) - 2 f'(y2));  sup |f''| = |f''(0)| = 2 for all three kernels.
V_FPP['phi'] = 0.96059156149540
V_FPP['kap'] = 1.09092291179862
M_FPP['phi'] = M_FPP['kap'] = 2.0
Z3H = float(7 * zeta(3))            # sum_n (n+1/2)^-3
Z2H = float(mp.pi ** 2 / 2)         # sum_n (n+1/2)^-2
THETA_STAR, S_STAR = None, None
TBAR = {}


def build_sups(verbose=True):
    global THETA_STAR, S_STAR
    THETA_STAR, S_STAR = theta_star_rig()
    for nm in ('G', 'phi', 'kap'):
        TBAR[nm] = sup_pad_profile(nm, 'T')
        TBAR['U' + nm] = sup_pad_profile(nm, 'U')
    build_schur_grids()
    if verbose:
        print(f"   sups built: Theta*={THETA_STAR:.9f} at s*={S_STAR:.6f}, "
              f"Tbar={ {k: round(v,4) for k,v in TBAR.items()} } "
              f"({time.time()-T0:.1f} s)", flush=True)


def eps1(nm, w):
    """sum_n |l_f(n) - Lam_f(s_n)/s_n| w_n G(s_n)^2  <=  ... (L3)+(L4)."""
    D, V = 1.0 / w, V_FPP[nm]
    q = D * D / 8.0
    return q * (np.log(1.0 + V * w) + V / (V + D))


def eps2(nm, w):
    """sum_n |l_f(n) - Lam_f(s_n)/s_n|^2 w_n G(s_n)^2  <=  ..."""
    D, V = 1.0 / w, V_FPP[nm]
    q = D * D / 8.0
    return Z3H * D ** 4 / 576.0 + q * q * (4.0 * np.log(1.0 + V * w) + 2.0)


def padX(nm, w):
    """two-sided pad of  <tau, R_f tau> = sum_n w_n G Ups_f^lat  vs  I[G Ups_f]."""
    D = 1.0 / w
    q = D * D / 8.0
    return q * CV['V_' + nm + '_GU'] + q * (CV['GU_' + nm] + (D / 2.0) * CV['V_GU_' + nm])


def bounds2(w, detail=False):
    """ALL bounds of M1 Sec.2, with every pad replaced by (L1)-(L4).  Formulas."""
    D = 1.0 / w
    q = D * D / 8.0
    ex = exact_sums(w)
    A, B = ex['A'], ex['B']
    rho = 2.0 * B - 1.0                                   # (a)  EXACT
    g = rho - THETA_STAR - q * TBAR['G']                  # (c)
    # ---- (b) residual ------------------------------------------------------
    beta = (D / 2.0 - np.arctan(D / 2.0)) + q * V_FPP['G']
    ND_hi = (CV['nD2'] + q * CV['V_nD2']
             + 2.0 * beta * (CV['J1'] + (D / 2.0) * CV['V_J1'])
             + beta ** 2 * Z3H / (2.0 * D * D))
    NR_hi = (np.sqrt(CV['nR2'] + q * CV['V_nR2'])
             + q * np.sqrt(CV['W_G'] + (D / 2.0) * CV['V_W_G']))
    Cr = (NR_hi + np.sqrt(ND_hi)) / np.sqrt(B)
    sn = Cr / g
    # ---- (d1) F(t) and (e) q_kappa ----------------------------------------
    out = {}
    for nm, Sgf in (('phi', ex['S_Gphi']), ('kap', ex['S_Gkap'])):
        pX = padX(nm, w)
        Yhi = CV[nm + '_tDt'] + q * CV['V_' + nm + '_tDt'] + eps1(nm, w)
        Ylo = max(0.0, CV[nm + '_tDt'] - q * CV['V_' + nm + '_tDt'] - eps1(nm, w))
        Xlo, Xhi = CV[nm + '_GU'] - pX, CV[nm + '_GU'] + pX
        out[nm + '_lo'] = 4.0 * Sgf + (Xlo - Yhi) / B
        out[nm + '_hi'] = 4.0 * Sgf + (Xhi - Ylo) / B
    F_t = out['phi_lo']
    q_k = max(abs(out['kap_lo']), abs(out['kap_hi']))
    # ---- (d2) e_phi, e_kappa ----------------------------------------------
    e = {}
    for nm, Sff, Sgf in (('phi', ex['S_pp'], ex['S_Gphi']),
                         ('kap', ex['S_kk'], ex['S_Gkap'])):
        V = lambda k: q * CV['V_' + nm + '_' + k]
        nZ2 = (4.0 * B * B * Sff
               - 4.0 * B * (CV[nm + '_fDt'] - V('fDt'))
               + (CV[nm + '_DD'] + V('DD'))
               + 4.0 * B * (CV[nm + '_fU'] + V('fU'))
               - 2.0 * (CV[nm + '_lU'] - V('lU'))
               + (CV[nm + '_UU'] + V('UU')))
        tZ = (2.0 * B * Sgf - (CV[nm + '_tDt'] + V('tDt'))
              + (CV[nm + '_GU'] - V('GU')))
        dZ = (np.sqrt(eps2(nm, w))
              + q * np.sqrt(CV['W_' + nm] + (D / 2.0) * CV['V_W_' + nm]))
        e[nm] = (np.sqrt(max(0.0, nZ2 - max(0.0, tZ) ** 2 / B)) + dZ) / np.sqrt(B)
        e[nm + '_tZ'] = tZ
    # ---- (d3) the two Schur norms -----------------------------------------
    NS = {}
    for nm, Sf in (('phi', ex['S_phi']), ('kap', ex['S_kap'])):
        top, head, tail = schur_sup(nm, A, Sf, w)
        NS[nm] = top + q * TBAR[nm] + D * D / 3.0
        NS[nm + '_head'], NS[nm + '_tail'], NS[nm + '_top'] = head, tail, top
    # ---- (e) Davis-Kahan assembly -----------------------------------------
    F_lo = (1.0 - sn * sn) * F_t - 2.0 * sn * e['phi'] - sn * sn * NS['phi']
    E_hi = e['kap'] + sn * (q_k + NS['kap'])
    ratio = E_hi ** 2 / (g * F_lo) if F_lo > 0 else np.inf
    r = dict(A=A, B=B, rho=rho, Cr=Cr, g=g, s=sn, F_t=F_t, q_k=q_k,
             e_phi=e['phi'], e_kap=e['kap'], NP=NS['phi'], NK=NS['kap'],
             F_lo=F_lo, E_hi=E_hi, ratio=ratio)
    if detail:
        r.update(ND_hi=ND_hi, NR_hi=NR_hi, beta=beta, NS=NS, e=e, out=out)
    return r


# ============================================================================
# SEC. 7.  STAGES -- every lemma and every pad verified against the lattice
# ============================================================================
def lattice(N, w):
    n = np.arange(N)
    wn = 1.0 / (2.0 * n + 1.0)
    return n, wn, (n + 0.5) / w


def S_lat(u, w, smax=400.0):
    """the genuine lattice sum S[u] = sum_n w_n u(s_n) (truth side of (L1))."""
    N = int(smax * w) + 10
    _, wn, sn = lattice(N, w)
    return float(np.sum(wn * u(sn)))


def run_stage1():
    print("\n=== Stage 1: the lemmas (L1),(L3) and the psi-jets ===")
    # (L1) on explicit test functions with known int h and int |h''|
    for D in (0.25, 0.1, 0.02):
        n = np.arange(int(200.0 / D))
        sn = (n + 0.5) * D
        for name, h, Ih, Vh in [
                ("e^-s", lambda s: np.exp(-s), 1.0, 2.0),
                ("s e^-s", lambda s: s * np.exp(-s), 1.0, 4.0 / np.e),
                ("1/(1+s^2)^2", lambda s: 1.0 / (1.0 + s * s) ** 2, np.pi / 4.0,
                 2.0 * 2.0)]:
            err = abs(float(np.sum(D * h(sn))) - Ih)
            chk(f"(L1) {name} D={D}", err <= D * D / 8.0 * Vh + 1e-12,
                f"{err:.3e} vs {D*D/8*Vh:.3e}")
    # eta of M1 (a):  V_eta = int|eta''| = 3/2 exactly, int eta = 1/2
    eta = lambda u: u / (1.0 + u * u) ** 2
    for D in (0.2, 0.05):
        n = np.arange(int(400.0 / D))
        err = abs(float(np.sum(D * eta((n + 0.5) * D))) - 0.5)
        chk(f"(L1) eta D={D}", err <= D * D / 8.0 * 1.5, f"{err:.3e}")
    # V_f'' constants
    for nm, V in V_FPP.items():
        gg = np.geomspace(1e-6, 200.0, 400001)
        num = float(np.trapezoid(np.abs(FAM[nm][2](gg)), gg))
        chk(f"int|f''| {nm}", abs(num - V) < 2e-5, f"{num:.8f} vs {V:.8f}")
        chk(f"sup|f''| {nm}", np.max(np.abs(FAM[nm][2](gg))) <= M_FPP.get(nm, 2.0) + 1e-9)
    # psi-jets:  A_j = S[G^j] exactly
    for w in (3.0, 8.0, 45.0, 1000.0):
        A = A_exact(w)
        for j in range(1, 7):
            lat = S_lat(lambda s, j=j: G(s) ** j, w, smax=4000.0)
            tol = 3e-7 if j == 1 else 1e-9
            chk(f"A_{j} exact w={w}", abs(A[j] - lat) < tol, f"{A[j]:.12f} {lat:.12f}")
        ex = exact_sums(w)
        chk(f"S[G kappa]=S[phi] w={w}", abs(ex['S_Gkap'] - ex['S_phi']) < 1e-14)
    # (L3): the partial sums of f(k D)
    for w in (8.0, 45.0):
        D = 1.0 / w
        for nm in ('G', 'phi', 'kap'):
            f = FAM[nm][0]
            N = int(60 * w)
            k = np.arange(1, N)
            cs = np.concatenate(([0.0], np.cumsum(f(k * D))))
            sn = (np.arange(N) + 0.5) * D
            lhs = np.abs(cs - (LAM_F[nm](sn) - LAM_F[nm](D / 2)) / D)
            chk(f"(L3) {nm} w={w}", np.all(lhs <= D / 8.0 * V_FPP[nm] + 1e-12),
                f"{lhs.max():.3e} vs {D/8*V_FPP[nm]:.3e}")
    stage_end("Stage 1")


def run_stage2():
    print("\n=== Stage 2: profiles Theta_f, Ups_f, their s-derivatives, T_f, U_f ===")
    for s in (0.05, 0.3, 1.0, 1.5009, 3.0, 10.0, 40.0):
        ex = (np.log1p(s * s) + 2 * s * np.arctan(s)) / (2 * (1 + s * s))
        chk(f"Theta_G closed form s={s}", abs(float(P('G', 'T', 0, s)) - ex) < 5e-7,
            f"{float(P('G','T',0,s)):.9f} vs {ex:.9f}")
    # d/ds and d2/ds2 of the profiles vs central differences of the profile itself
    for nm in ('G', 'phi', 'kap'):
        for tag in ('T', 'U'):
            for s in (0.4, 1.0, 2.5, 6.0):
                hh = 1e-3 * s
                f0 = prof_int(nm, 0, s, tag == 'U')
                fp = prof_int(nm, 0, s + hh, tag == 'U')
                fm = prof_int(nm, 0, s - hh, tag == 'U')
                d1, d2 = (fp - fm) / (2 * hh), (fp - 2 * f0 + fm) / hh ** 2
                chk(f"d/ds {nm}{tag} s={s}", abs(d1 - float(P(nm, tag, 1, s))) < 2e-5,
                    f"{d1:.7f} vs {float(P(nm,tag,1,s)):.7f}")
                chk(f"d2/ds2 {nm}{tag} s={s}", abs(d2 - float(P(nm, tag, 2, s))) < 3e-3,
                    f"{d2:.6f} vs {float(P(nm,tag,2,s)):.6f}")
    # T_f, U_f really majorise the inner-sum discretisation error (L1 in sigma)
    print("   inner-sum pad q*T_f, q*U_f vs measured |lattice - continuum|:")
    for w in (8.0, 45.0):
        D, q = 1.0 / w, 1.0 / (8.0 * w * w)
        N = int(300 * w)
        _, wn, sn = lattice(N, w)
        for nm in ('G', 'phi', 'kap'):
            for tag, wt in (('T', 1.0), ('U', None)):
                worst = 0.0
                for s in (0.3, 1.0, 2.0, 5.0):
                    wm = wn if wt else wn * G(sn)
                    lat = float(np.sum(wm * Rk(nm, 0, 0, s, sn)))
                    con = float(P(nm, tag, 0, s))
                    pad = q * float(P(nm, tag, 'pad', s))
                    worst = max(worst, abs(lat - con) / pad)
                    chk(f"inner pad {nm}{tag} w={w} s={s}", abs(lat - con) <= pad,
                        f"{abs(lat-con):.3e} vs {pad:.3e}")
                print(f"     w={w:6g} {nm}-{tag}: max |lat-cont| / pad = {worst:.3f}")
    stage_end("Stage 2")


def run_stage3():
    print("\n=== Stage 3: the constants and their V2 pads vs the true lattice error ===")
    import p4_cii_m1 as m1
    for k, v in m1.CONST.items():
        if k in CV:
            chk(f"const {k} == M1", abs(CV[k] - v) < 3e-6, f"{CV[k]:.9f} vs {v:.9f}")
    chk("phi_tDt closed form", abs(CV['phi_tDt'] - (7 / 48 - 3 * np.pi ** 2 / 256)) < 3e-8,
        f"{CV['phi_tDt']:.9f}")
    chk("kap_tDt closed form", abs(CV['kap_tDt'] - (0.5 - 3 * np.pi ** 2 / 64)) < 3e-8,
        f"{CV['kap_tDt']:.9f}")
    # every V2 pad tested against the measured midpoint error of its own sum
    s_ = _DS
    tests = []
    jsG, jG, _, jd, _, jUG, _ = _jets('G', s_)
    U_of = {nm: (lambda ss, nm=nm: P(nm, 'U', 0, ss)) for nm in ('G', 'phi', 'kap')}
    tests.append(('nD2', lambda ss: (1 - np.arctan(ss) / ss) ** 2 * G(ss) ** 2))
    tests.append(('nR2', lambda ss: U_of['G'](ss) ** 2))
    for nm in ('phi', 'kap'):
        f = FAM[nm][0]
        Lo = lambda ss, nm=nm: LAM_F[nm](ss) / ss
        tests += [(nm + '_GU', lambda ss, nm=nm: G(ss) * U_of[nm](ss)),
                  (nm + '_tDt', lambda ss, Lo=Lo: Lo(ss) * G(ss) ** 2),
                  (nm + '_fDt', lambda ss, f=f, Lo=Lo: f(ss) * Lo(ss) * G(ss)),
                  (nm + '_DD', lambda ss, Lo=Lo: Lo(ss) ** 2 * G(ss) ** 2),
                  (nm + '_fU', lambda ss, f=f, nm=nm: f(ss) * U_of[nm](ss)),
                  (nm + '_lU', lambda ss, Lo=Lo, nm=nm: Lo(ss) * G(ss) * U_of[nm](ss)),
                  (nm + '_UU', lambda ss, nm=nm: U_of[nm](ss) ** 2)]
    print("   key = |S[u]-I[u]| / (q V2[u])   (must be <= 1 for every wbar)")
    for w in (6.0, 8.0, 10.0, 20.0, 45.0, 100.0, 1000.0):
        q = 1.0 / (8.0 * w * w)
        row = []
        for key, u in tests:
            err = abs(S_lat(u, w) - CV[key])
            pad = q * CV['V_' + key]
            chk(f"V2 pad {key} w={w}", err <= pad, f"{err:.3e} vs {pad:.3e}")
            row.append(err / pad)
        print(f"     wbar={w:6g}:  max = {max(row):.3f}  (at {tests[int(np.argmax(row))][0]})")
    stage_end("Stage 3")


def run_stage4(ws):
    print("\n=== Stage 4: the (L3)/(L4) lattice bounds and the Schur row sums ===")
    for w in ws[:4]:
        D, q = 1.0 / w, 1.0 / (8.0 * w * w)
        N = int(300 * w)
        n, wn, sn = lattice(N, w)
        # 1 - D_nn  vs  d(s_n)
        lG = 2.0 * np.concatenate(([0.0], np.cumsum(G(np.arange(1, N) / w)))) / (2 * n + 1)
        beta = (D / 2.0 - np.arctan(D / 2.0)) + q * V_FPP['G']
        dd = 1.0 - np.arctan(sn) / sn
        err = np.abs((1.0 - wn - lG) - dd)
        chk(f"(L3)+w_n cancellation w={w}", np.all(err <= beta / sn + 1e-14),
            f"{np.max(err*sn):.3e} vs {beta:.3e}")
        # l_f(n) vs Lam_f/s_n
        for nm in ('phi', 'kap'):
            lf = 2.0 * np.concatenate(
                ([0.0], np.cumsum(FAM[nm][0](np.arange(1, N) / w)))) / (2 * n + 1)
            eps = (LAM_F[nm](D / 2) + q * np.minimum(2.0 * n * D, V_FPP[nm])) / sn
            e2 = np.abs(lf - LAM_F[nm](sn) / sn)
            chk(f"eps_f {nm} w={w}", np.all(e2 <= eps + 1e-15),
                f"{np.max(e2/eps):.3f}")
            # and the two summed versions
            s1 = float(np.sum(e2 * wn * G(sn) ** 2))
            s2 = float(np.sum(e2 ** 2 * wn * G(sn) ** 2))
            chk(f"eps1 {nm} w={w}", s1 <= eps1(nm, w), f"{s1:.3e} vs {eps1(nm,w):.3e}")
            chk(f"eps2 {nm} w={w}", s2 <= eps2(nm, w), f"{s2:.3e} vs {eps2(nm,w):.3e}")
        # Schur row sums of R_G:  sup_n sum_m R(s_n,s_m) w_m <= Theta* + q Tbar_G
        idx = np.unique(np.concatenate([np.arange(0, 30),
                                        np.geomspace(30, N - 1, 300).astype(int)]))
        Sig = max(float(np.sum(Rk('G', 0, 0, sn[i], sn) * wn)) for i in idx)
        chk(f"row sum <= Theta*+q Tbar w={w}", Sig <= THETA_STAR + q * TBAR['G'],
            f"{Sig:.8f} vs {THETA_STAR + q*TBAR['G']:.8f}")
        # the full Schur bound for O[f]
        for nm in ('phi', 'kap'):
            f = FAM[nm][0]
            lf = 2.0 * np.concatenate(
                ([0.0], np.cumsum(f(np.arange(1, N) / w)))) / (2 * n + 1)
            rows = []
            for i in idx:
                rows.append(float(np.sum(wn * (f(np.abs(sn[i] - sn)) + f(sn[i] + sn))))
                            + lf[i])
            ex = exact_sums(w)
            b = bounds2(w)
            NS = b['NP'] if nm == 'phi' else b['NK']
            chk(f"Schur sup {nm} w={w}", max(rows) <= NS, f"{max(rows):.6f} vs {NS:.6f}")
        print(f"   wbar={w:6g}: row-sum(R_G) {Sig:.7f} <= {THETA_STAR + q*TBAR['G']:.7f}"
              f" ; Schur(phi) {max(rows):.5f}")
    stage_end("Stage 4")


def run_stage5():
    print("\n=== Stage 5: Theta* and the Schur sup (monotonicity + grid + tail) ===")
    s = np.linspace(1e-6, 400.0, 400001)
    Np = (2.0 + (1.0 - 3.0 * s * s) / (1.0 + s * s) - 2.0 * s * np.arctan(s)
          - np.log1p(s * s))
    chk("N' strictly decreasing", np.all(np.diff(Np) < 0), "")
    chk("N'(0+)>0", Np[0] > 2.9)
    Nf = 2.0 * s + (1.0 - s * s) * np.arctan(s) - s * np.log1p(s * s)
    i = int(np.argmax(Nf < 0))
    chk("N has a single sign change", np.all(Nf[:i] > 0) and np.all(Nf[i:] < 0),
        f"root near {s[i]:.6f}")
    chk("Theta* value", abs(THETA_STAR - 0.63492879) < 1e-7, f"{THETA_STAR:.9f}")
    chk("s* value", abs(S_STAR - 1.500904) < 1e-5, f"{S_STAR:.6f}")
    TH = lambda x: (np.log1p(x * x) + 2 * x * np.arctan(x)) / (2 * (1 + x * x))
    chk("Theta <= Theta* on a dense grid", np.max(TH(s)) <= THETA_STAR, "")
    for w in (8.0, 45.0, 1000.0):
        ex = exact_sums(w)
        for nm, Sf in (('phi', ex['S_phi']), ('kap', ex['S_kap'])):
            top, head, tail = schur_sup(nm, ex['A'], Sf, w)
            chk(f"Schur head < top {nm} w={w}", head <= top, f"{head:.5f} {top:.5f}")
            chk(f"Schur tail < top {nm} w={w}", tail <= top, f"{tail:.5f} {top:.5f}")
            # the tail majorant is decreasing beyond S_SUP
            ss = np.linspace(S_SUP, 4000.0, 40000)
            tj = tail_major(nm, ss, ex['A'], Sf, w)
            chk(f"tail majorant decreasing {nm} w={w}", np.all(np.diff(tj) < 0))
    print(f"   Theta* = {THETA_STAR:.9f} at s* = {S_STAR:.6f} (bisection width 1e-40)")
    stage_end("Stage 5")


def run_stage6(ws):
    print("\n=== Stage 6: the rigorous bounds vs float64 truth ===")
    import p4_cii_m1 as m1
    print("  wbar |  rho_lo     h     | C_r    ||r||  | F_lo(t)  F(t)  |"
          " e_kap  exact | q_k   true | sin_lo  sin")
    T = {}
    for w in ws:
        b, t = bounds2(w), m1.truth(w)
        T[w] = (b, t)
        chk(f"(a) rho_lo<=h w={w}", b['rho'] <= t['h'], f"{b['rho']:.6f} {t['h']:.6f}")
        chk(f"(b) C_r>=||r|| w={w}", b['Cr'] >= t['nr'], f"{b['Cr']:.6f} {t['nr']:.6f}")
        chk(f"(c) g_lo<=true gap w={w}", b['g'] <= t['h'] + 1e-9)
        chk(f"(d1) F_lo(t)<=F(t) w={w}", b['F_t'] <= t['F_t'], f"{b['F_t']:.7f} {t['F_t']:.7f}")
        chk(f"(d2) e_kap>= w={w}", b['e_kap'] >= t['e_kap'], f"{b['e_kap']:.7f} {t['e_kap']:.7f}")
        chk(f"(d2) e_phi>= w={w}", b['e_phi'] >= t['e_phi'], f"{b['e_phi']:.7f} {t['e_phi']:.7f}")
        chk(f"(e) |q_k| bound w={w}", b['q_k'] >= abs(t['q_k']), f"{b['q_k']:.6f} {t['q_k']:.6f}")
        chk(f"(d3) NP>=||O[phi]that|| w={w}", b['NP'] >= t['nP'], f"{b['NP']:.5f} {t['nP']:.5f}")
        chk(f"(d3) NK>=||K that|| w={w}", b['NK'] >= t['nK'], f"{b['NK']:.5f} {t['nK']:.5f}")
        chk(f"(e) sin bound w={w}", b['s'] >= t['sin'], f"{b['s']:.5f} {t['sin']:.5f}")
        chk(f"(e) E_hi>=E(v) w={w}", b['E_hi'] >= t['E_v'], f"{b['E_hi']:.5f} {t['E_v']:.5f}")
        chk(f"(e) F_lo<=F(v) w={w}", b['F_lo'] <= t['F_v'], f"{b['F_lo']:.5f} {t['F_v']:.5f}")
        print(f" {w:6g} | {b['rho']:9.5f} {t['h']:9.5f} | {b['Cr']:.4f} {t['nr']:.4f} |"
              f" {b['F_t']:.5f} {t['F_t']:.5f} | {b['e_kap']:.5f} {t['e_kap']:.5f} |"
              f" {b['q_k']:.4f} {t['q_k']:.4f} | {b['s']:.5f} {t['sin']:.5f}")
    stage_end("Stage 6")
    return T


def wbar_A(fn=None):
    fn = fn or (lambda x: bounds2(x)['ratio'])
    lo, hi = 2.0, 400.0
    for _ in range(70):
        mid = 0.5 * (lo + hi)
        try:
            bad = fn(mid) > 1.0
        except Exception:
            bad = True
        if bad:
            lo = mid
        else:
            hi = mid
    return hi


def pad_report(w):
    """the complete pad catalogue at one wbar:  M1's blanket PAD vs M2's formulas."""
    D, q = 1.0 / w, 1.0 / (8.0 * w * w)
    beta = (D / 2.0 - np.arctan(D / 2.0)) + q * V_FPP['G']
    rows = [
        ("g_lo      sup_n R-row sum", "q Tbar_G", q * TBAR['G']),
        ("C_r       ||(I-D)tau||^2 ", "q V2[d^2G^2] + 2b(J1+..) + b^2 7z(3)/2D^2",
         q * CV['V_nD2'] + 2 * beta * (CV['J1'] + D / 2 * CV['V_J1'])
         + beta ** 2 * Z3H / (2 * D * D)),
        ("C_r       ||R tau||^2    ", "q V2[Ups_G^2]", q * CV['V_nR2']),
        ("C_r       inner (Ups)    ", "q sqrt(W_G+..)",
         q * np.sqrt(CV['W_G'] + D / 2 * CV['V_W_G'])),
        ("F(t)      <tau,R_phi tau>", "padX_phi", padX('phi', w)),
        ("F(t)      <tau,D_phi tau>", "q V2[(Lam/s)G^2] + eps1",
         q * CV['V_phi_tDt'] + eps1('phi', w)),
        ("q_kappa   <tau,R_kap tau>", "padX_kap", padX('kap', w)),
        ("q_kappa   <tau,D_kap tau>", "q V2[(Lam/s)G^2] + eps1",
         q * CV['V_kap_tDt'] + eps1('kap', w)),
        ("e_kappa   ||Z||^2 (5 sums)", "q sum V2[.]",
         q * sum(CV['V_kap_' + k] for k in ('fDt', 'DD', 'fU', 'lU', 'UU'))),
        ("e_kappa   <tau,Z>        ", "q(V2[tDt]+V2[GU])",
         q * (CV['V_kap_tDt'] + CV['V_kap_GU'])),
        ("e_kappa   Delta Z        ", "sqrt(eps2)+q sqrt(W_kap+..)",
         np.sqrt(eps2('kap', w)) + q * np.sqrt(CV['W_kap'] + D / 2 * CV['V_W_kap'])),
        ("||K||     inner + l_kap  ", "q Tbar_kap + D^2/3",
         q * TBAR['kap'] + D * D / 3.0),
        ("||O[phi]|| inner + l_phi ", "q Tbar_phi + D^2/3",
         q * TBAR['phi'] + D * D / 3.0),
    ]
    print(f"\n   PAD CATALOGUE at wbar = {w:g}   (M1 used the single calibrated"
          f" PAD = 0.25/wbar^2 + 2e-6 = {0.25/w/w + 2e-6:.3e})")
    for a, b, v in rows:
        print(f"     {a} :  {v:.3e}   = {b}")
    return rows


def run_stage7(ws, T):
    print("\n=== Stage 7: wbar_A' , the criterion table, and the M1 comparison ===")
    B1 = {45.0: None, 100.0: -0.493280384, 200.0: -0.494738431, 1000.0: -0.496742993}
    print("  wbar |  g_lo   sin(th)  |  F_lo    E_hi   | RATIO  | true ratio | x^2H''<= ")
    for w in ws:
        b, t = T[w]
        x2H = -2.0 * b['F_lo'] + 2.0 * b['E_hi'] ** 2 / b['g']
        if B1.get(w):
            chk(f"x^2H'' >= truth w={w}", x2H >= B1[w], f"{x2H:.6f} {B1[w]:.6f}")
            chk(f"x^2H'' non-vacuous w={w}", x2H < 0.55 * B1[w], f"{x2H:.6f}")
        print(f" {w:6g} | {b['g']:7.4f} {b['s']:.5f} | {b['F_lo']:.5f} {b['E_hi']:.5f} |"
              f" {b['ratio']:.4f} | {t['ratio']:.4f}     | {x2H:9.6f}")
    wA = wbar_A()
    pad_report(45.0)
    chk("wbar_A' found", 2.0 < wA < 45.0, f"{wA}")
    chk("closes above wbar_A'", bounds2(wA * 1.001)['ratio'] < 1.0)
    chk("fails below wbar_A'", bounds2(wA * 0.999)['ratio'] > 1.0)
    chk("wbar_A' <= 45 (overlap with Segment C-I)", wA <= 45.0, f"{wA:.4f}")
    grid = np.geomspace(wA * 1.0001, 1e20, 200)
    rr = np.array([bounds2(x)['ratio'] for x in grid])
    chk("ratio<1 on [wbar_A',1e20]", float(rr.max()) < 1.0, f"{rr.max():.4f}")
    chk("ratio -> 1/4", abs(rr[-1] - 0.25) < 0.02, f"{rr[-1]:.5f}")
    chk("ratio monotone at large wbar", rr[-1] < rr[len(rr) // 2] < rr[2])
    print(f"\n   ***  wbar_A' = {wA:.4f}  (RIGOROUS pads; M1's calibrated value was"
          f" 7.776)  ***")
    print(f"   ratio: {bounds2(45.0)['ratio']:.4f} (wbar=45), "
          f"{bounds2(1000.0)['ratio']:.4f} (wbar=1000), {rr[-1]:.5f} (wbar=1e20) -> 1/4")
    print(f"   Segment C-I certifies wbar <= 45  ==>  overlap factor {45/wA:.2f} in wbar"
          f" ({45**2/wA**2:.1f} in x).")
    # side by side with M1
    import p4_cii_m1 as m1
    m1.SPL = m1.build_profiles()
    print("\n   M1 (calibrated pad) vs M2 (rigorous pads):")
    print("    wbar |  g_lo M1 / M2    |  F_lo M1 / M2   |  E_hi M1 / M2   | ratio M1 / M2")
    for w in ws:
        a, b = m1.bounds(w), bounds2(w)
        chk(f"M2 ~ M1 g_lo w={w}", abs(b['g'] - a['g']) < 0.02 * a['g'],
            f"{a['g']:.6f} {b['g']:.6f}")
        chk(f"M2 ratio >= M1 - 1e-3 w={w}", b['ratio'] >= a['ratio'] - 1e-3,
            f"{a['ratio']:.6f} {b['ratio']:.6f}")
        print(f"  {w:6g} | {a['g']:8.4f} {b['g']:8.4f} | {a['F_lo']:.5f} {b['F_lo']:.5f}"
              f" | {a['E_hi']:.5f} {b['E_hi']:.5f} | {a['ratio']:.4f} {b['ratio']:.4f}")
    stage_end("Stage 7")
    return wA


def main():
    quick = "--quick" in sys.argv
    ws = (8.0, 10.0, 20.0, 45.0) if quick else (8.0, 10.0, 20.0, 45.0, 100.0, 1000.0)
    print("P4 Segment C-II Milestone M2 -- the rigorous pads (M1 HONEST LIMIT 2 closed)")
    build_profiles()
    build_constants()
    build_sups()
    run_stage1()
    run_stage2()
    run_stage3()
    run_stage4(ws)
    run_stage5()
    T = run_stage6(ws)
    wA = run_stage7(ws, T)
    print(f"\n=== ALL {PASS[0]} CHECKS PASSED in {time.time()-T0:.1f} s ===")
    print(f"=== wbar_A' = {wA:.4f} ; ratio {bounds2(45.0)['ratio']:.4f} (45), "
          f"{bounds2(1000.0)['ratio']:.4f} (1000) ===")


if __name__ == "__main__":
    main()

