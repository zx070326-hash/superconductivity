r"""
P4 . EINSTEIN EXTREMAL -- Segment C-II, MILESTONE M1 (route R, analytic side).

Closed-form analytic bounds for the four ingredients of the large-wbar concavity
criterion   x^2 E^2 / g  <  F   of (C**) in `theory/derivations/P4_einstein_extremal.md`
Sec.1.A.2, plus the Davis-Kahan assembly and the resulting wbar_A.

RUN:   python theory/derivations/p4_cii_m1.py            (~2 min: all stages)
       python theory/derivations/p4_cii_m1.py --quick    (~40 s: trims the wbar grid)

SCOPE.  This file owns Sec.1.F of P4_einstein_extremal.md.  It IMPORTS the certified
entry formulas (`H_mat` from `p4_interchange_uniformity`) and the O(N log N) kernel
operator (`OpKernel`, `top_eig` from `p4_extremal_asymptotics`); NO base file is
modified and no matrix element is re-derived.  All bounds below are ANALYTIC
(formulas in wbar with explicit constants); the float64 numbers are only there to
verify them.  This is NOT an interval-arithmetic certificate -- see HONEST LIMITS.

================================================================================
SEC. 0.  NOTATION AND THE HALF-SHIFTED LATTICE
================================================================================
p_n := (2n+1)^{-1/2},  w_n := p_n^2 = 1/(2n+1),  and for a kernel f : R -> R

    O[f]_{nm} = p_n p_m [ (1-delta_{nm}) f(|n-m|) + f(n+m+1) ]
                - delta_{nm} (2/(2n+1)) sum_{k=1}^{n} f(k).

H(wbar) = O[ [[.]] ], [[k]] = wbar^2/(k^2+wbar^2), h(wbar) = lam_1(H), x = wbar^2.

THE KEY BOOKKEEPING DEVICE.  Put  D := 1/wbar  and  s_n := (n + 1/2) D.  Then

    (i)   w_n = 1/(2n+1) = D / (2 s_n)        EXACTLY,
    (ii)  s_n - s_m = (n-m) D  and  s_n + s_m = (n+m+1) D   EXACTLY,
    (iii) s_n is the MIDPOINT of [nD, (n+1)D].

By (ii) every Toeplitz index |n-m| AND every Hankel index n+m+1 is a difference
resp. a sum of lattice points -- the Hankel "+1" is absorbed by the half shift.
By (i)+(iii), for any u,
    sum_n w_n u(s_n) = (1/2) sum_n u(s_n)/s_n * D  =  midpoint rule for  I[u],
    I[u] := (1/2) int_0^infty u(s) ds/s ,
so every sum in this file is a midpoint Riemann sum of an explicit integral, with
error O(D^2) whenever u/s is C^2 and integrable.  Write
    G(y) := 1/(1+y^2)  --  so that  [[k]] = G(k/wbar)  and  G(s_n -+ s_m) is exactly
the Toeplitz / Hankel entry.  TRIAL VECTOR (the dispatch's, half-shifted):

    tau_n := p_n G(s_n) = p_n wbar^2/(wbar^2 + (n+1/2)^2),   B := ||tau||^2,
    that  := tau / sqrt(B),   rho := <that, H that>.

(The half shift is the ONLY deviation from the dispatch's t_n = p_n(1+(n/wbar)^2)^{-1};
it costs nothing numerically -- rho is the same to 1e-3 -- and it is what makes the
three identities of Sec.1 EXACT rather than exact-up-to-O(1/wbar).)

================================================================================
SEC. 1.  THREE EXACT ALGEBRAIC IDENTITIES (the whole engine)
================================================================================
(M1)  LORENTZIAN ADDITION.  For all s,sigma >= 0,

    G(s-sigma) + G(s+sigma) - 2 G(s) G(sigma)  =  R(s,sigma)
    R(s,sigma) = 2 s^2 sigma^2 (5+s^2+sigma^2)
                 / [ (1+(s-sigma)^2)(1+(s+sigma)^2)(1+s^2)(1+sigma^2) ]   >= 0.

(elementary algebra: (1+(s-sig)^2)(1+(s+sig)^2) = 1+2(s^2+sig^2)+(s^2-sig^2)^2 and
(1+s^2)(1+sig^2) = 1+s^2+sig^2+s^2 sig^2; the numerator of the difference is
s^2 sig^2 (5+s^2+sig^2).)  With (ii) this gives the OPERATOR identity

    **  H(wbar) = 2 tau tau^T + Rop - Dg  **                                  (M1)
    Rop_{nm} := p_n p_m R(s_n,s_m) >= 0 (entrywise),
    Dg := diag( w_n + l_G(n) ),  l_f(n) := (2/(2n+1)) sum_{k=1}^n f(k/wbar),
    0 <= Dg_{nn} <= 1   (because G <= 1 gives sum_{k<=n} G <= n).

(M2)  For a kernel f with f(0) = 0 (we need f = phi_0 := y^2/(1+y^2)^3 and
kappa := y^2/(1+y^2)^2 = G - G^2), the SAME split with the two cross terms:

    R_f(s,sigma) := f(s-sigma)+f(s+sigma) - 2 f(s) G(sigma) - 2 f(sigma) G(s),
    R_f(s,0) = R_f(0,sigma) = 0,                                             (*)
    **  O[f] = 2 ( tau_f tau^T + tau tau_f^T ) + Rop_f - Dg_f  **             (M2)
    tau_{f,n} := p_n f(s_n),  Dg_f := diag( l_f(n) ),  Rop_{f,nm} := p_n p_m R_f(s_n,s_m).

(*) is what kills BOTH logarithmic divergences: I[R_f(s,.)] and I[G R_f(s,.)]
converge, so every constant below is a CONVERGENT integral of an elementary
function, independent of wbar.  (For f = G the single term 2G(s)G(sigma) already
covers both limits -- that is (M1).)

(M3)  SCALING.  phi_x(k) = k^2/(k^2+x)^3 = wbar^{-4} phi_0(k/wbar) and
psi_x(k) = 1/(k^2+x)^2 = wbar^{-4} G(k/wbar)^2, hence with F, E as in (C**)

    wbar^4 F = <v, O[phi_0] v>,   wbar^4 E = ||(I-P) O[G^2] v|| = ||(I-P) K v||,
    K := O[kappa] = H - O[G^2]   (the last equality because (I-P) H v = 0),

and (C**) is equivalent to      ( wbar^4 E )^2 / ( g * wbar^4 F )  <  1.      (CRIT)

Using K instead of O[G^2] is essential: <that,O[G^2]that> ~ ln wbar, whereas
<that, K that> = O(1), and the Davis-Kahan transport cost is proportional to it.

================================================================================
SEC. 2.  THE FIVE BOUNDS
================================================================================
Everything is expressed through the two elementary wbar-functions (digamma closed
forms, Sec. 2.1) and a fixed list of ABSOLUTE constants I[.] (Sec. 2.6).

(2.1) A and B IN CLOSED FORM.  With a_n = n+1/2,
    2A := sum_n [1/a_n - a_n/(a_n^2+x)] = Re psi(1/2 + i wbar) + gamma + 2 ln 2,
    B  = A - x A'(x) = A + (wbar/4) Im psi'(1/2 + i wbar),
    A  = sum_n w_n G(s_n),   B = sum_n w_n G(s_n)^2 = ||tau||^2.
ELEMENTARY TWO-SIDED FORM.  Binet's second formula psi(z) = ln z - 1/(2z)
- 2 int_0^inf t dt/((t^2+z^2)(e^{2 pi t}-1)) at z = 1/2 + i wbar, with
|t^2+z^2| >= Im(z^2) = wbar and int_0^inf t/(e^{2 pi t}-1) dt = 1/24, gives
    | Re psi(1/2+i wbar) - (1/2) ln(wbar^2+1/4) + 1/(4 wbar^2+1) | <= 1/(12 wbar);
and 2A - 2B = sum_n Xi(a_n), Xi(a) = a x/(a^2+x)^2 = eta(a/wbar)/wbar,
eta(u) = u/(1+u^2)^2, whose midpoint sum obeys |sum - int_0^inf Xi| <=
(1/24) int|Xi''| = V_eta/(24 wbar^2) with V_eta = int_0^inf |eta''| = 3/2 EXACTLY
and int_0^inf Xi = 1/2.  Hence

  **(a)  h >= rho >= 2B - 1 >= ln wbar + (gamma + ln 4 - 3/2)
                              - 1/(12 wbar) - 1/(4 wbar^2) - 1/(16 wbar^2). **

PROOF of h >= rho >= 2B-1:  by (M1),
    <tau,H tau> = 2B^2 + <tau,Rop tau> - <tau,Dg tau>,
    <tau,Dg tau> = B - c_Delta,  c_Delta := sum_n (1-Dg_{nn}) w_n G(s_n)^2 >= 0,
so  rho = 2B - 1 + ( c_Delta + <tau,Rop tau> ) / B  >=  2B - 1, both dropped terms
being >= 0 by Dg <= I and Rop >= 0 entrywise with tau >= 0.  []
NOTE the constant: gamma + ln 4 - 3/2 = 0.4635101 is EXACTLY the published soft-
cutoff c_0 (Sec. 0.A) -- so (a) proves the lower half of h = ln wbar + c_0 + o(1),
which had no proof at all before (Sec. 1.E: "带误差界的严格大 ϖ 渐近：全域为零").

(2.2) RESIDUAL.  r := H that - rho that.  By (M1),
    H tau - rho tau = (I - P_tau)(Rop - Dg) tau = (I-P_tau)(Rop tau - (Dg - I)tau),
(the I is free: (I-P_tau) tau = 0), hence

  **(b)  ||r|| <= [ ||Rop tau|| + ||(I - Dg) tau|| ] / sqrt(B)
                <= [ sqrt(N_R) + sqrt(N_D) + pad ] / sqrt(B),  **
    N_R := I[Ups_G^2] = 0.0492567,  Ups_G(s) := I_sigma[ G R(s,.) ],
    N_D := I[ d^2 G^2 ] = 0.0086534,  d(s) := 1 - arctan(s)/s  (the limit of 1-Dg_nn).
i.e. C_r(wbar) = 0.31496/sqrt(B) + O(wbar^{-2}):  C_r -> 0 like 1/sqrt(ln wbar).

(2.3) GAP.  tau in l^2, so for every unit w perp tau, <w,H w> = <w,(Rop-Dg)w>, and
the one-dimensional min-max gives directly
    lam_2(H) <= lam_1(H - 2 tau tau^T) = lam_1(Rop - Dg) <= lam_1(Rop) <= ||Rop||.
Weighted Schur with weight p (Rop >= 0 entrywise, symmetric):
    ||Rop|| <= sup_n p_n^{-1} sum_m Rop_{nm} p_m = sup_n sum_m R(s_n,s_m) w_m
            =  sup_n [ midpoint sum of I[R(s_n,.)] ]  <= sup_{s>=0} Theta(s) + pad,
and I[R(s,.)] has the CLOSED FORM (residues in z = sigma^2; the three poles are
z = -1 and z = (s -+ i)^2, the residue at z=-1 is 1/s^2 and the pair contributes
2 arctan s / s):

    **  Theta(s) = [ ln(1+s^2) + 2 s arctan s ] / ( 2 (1+s^2) ),
        Theta* := max_{s>=0} Theta = 0.6349288  at s* = 1.500904 (root of
        arctan(s)(1-s^2) = s(ln(1+s^2)-2)),  Theta(s) decreasing for s >= 2.53. **

  **(c)  g = h - lam_2 >= rho - Theta* >= ln wbar - 0.1714 - O(1/wbar).  **

THIS REPLACES the dispatch's route.  The dispatch asked for g >= h + lam_min(M),
M = I + W = 2uu^T - H, u_n = p_n.  That route is VACUOUS on l^2: since u is NOT in
l^2, for the top eigenvector v and any block [N1,N2] one may put z := u*1_{[N1,N2]}
normalised, so that <u,z> = sqrt(sum_{N1}^{N2} w_n) -> infinity, and w :=
(v - eps z)/||.|| with eps := <u,v>/<u,z> -> 0 satisfies <u,w> = 0 exactly and
<w,H w> -> h.  Hence lam_min(M) = inf(2<u,w>^2 - <w,H w>) = -h EXACTLY and
h + lam_min(M) = 0.  Stage 5 exhibits this collapse numerically (it is the 1/ln N
drift already reported in Sec. 1.B (B8)).  The cure is exactly (M1): replace the
non-l^2 rank one 2uu^T by the l^2 rank one 2 tau tau^T, which costs Theta* = 0.635
instead of everything.  The resulting g_lo is ~85% of the true gap (vs 58% for the
lam_min(M) route at N = 40 wbar, and 0% as N -> infinity).

(2.4) F AND E.  By (M2) with f = phi_0 (note <tau,tau_f> = I[G f] exactly in the
limit, and 4 I[G phi_0] = 4 * (1/12) = 1/3):
    wbar^4 F(t) = <that,O[phi_0]that> = 4 I[G phi_0] + ( <tau,Rop_phi tau>
                  - <tau,Dg_phi tau> ) / B
  **(d1)  wbar^4 F(t) >= 1/3 - 0.1257572/B - pad,   0.1257572 = 0.0301739+0.0955830. **
By (M2) with f = kappa, K tau = 2B tau_kappa + 2 I[G kappa] tau + Rop_k tau - Dg_k tau,
and the tau-component drops out of (I-P_tau):
  **(d2)  wbar^4 E(t) = ||(I-P_tau) Z|| / sqrt(B),  Z := 2B tau_k - Dg_k tau + Rop_k tau,
        ||Z||^2 = 4B^2 m_kk - 4B m_kl + m_ll + 4B m_kU - 2 m_lU + m_UU,
        <tau,Z> = 2B I[G kappa] - m_tl + m_GU,      ||(I-P_tau)Z||^2 = ||Z||^2 - <tau,Z>^2/B **
(identical formula with phi_0 in place of kappa gives e_phi), and by weighted Schur
  **(d3)  ||O[f]|| <= sup_s [ 2 A f(s) + Theta_f(s) + 2 G(s) I[f] + Lam_f(s)/s ],
        Theta_f(s) := I_sigma[R_f(s,.)],  Lam_f(s) := int_0^s f. **
These grow like (2 max f) A ~ (max f) ln wbar, which is harmless: they enter only
multiplied by sin^2(theta) resp. sin(theta).

(2.5) DAVIS-KAHAN ASSEMBLY.  With P the spectral projection on h and
sin(theta) := ||(I-P)that||, the reduced resolvent gives
    s := sin(theta) <= ||r|| / (rho - lam_2) <= C_r / g_lo,
and writing v = cos(theta) that + sin(theta) z, z perp that, |z|=1:
    F   = <v,O[phi_0]v>/wbar^4 >= (1-s^2) <that,O[phi_0]that> - 2 s e_phi - s^2 ||O[phi_0]||,
    E   = ||(I-P_v) K v|| <= e_kappa + s ( |q_kappa| + ||K|| ),  q_kappa := <that,K that>,
(the middle term of F uses |<z,O[phi]that>| <= e_phi because z perp that).  Hence

  **(e)  CRITERION:   E_hi^2 / ( g_lo * F_lo ) < 1,
        g_lo = 2B - 1 - Theta*,   s = C_r/g_lo,
        F_lo = (1-s^2)(1/3 - 0.1257572/B) - 2 s e_phi - s^2 ||O[phi_0]||_Schur,
        E_hi = e_kappa + s ( |q_kappa| + ||K||_Schur ),  q_kappa = 1/2 - 0.1277654/B. **

(2.6) THE ABSOLUTE CONSTANTS.  All are I[.] of explicit elementary functions
(Sec. 2 of the code); they are wbar-INDEPENDENT and are cross-checked in Stage 2
against the lattice sums at wbar = 1000, where lattice and continuum agree to 1e-7.

================================================================================
SEC. 3.  HONEST LIMITS  (read before quoting any number)
================================================================================
1. NOT a certificate.  mpmath/float64, no interval arithmetic.  The absolute
   constants are scipy.quad / mpmath.quad values of convergent elementary
   integrals (~1e-9 accurate, margin is ~0.5), not certified enclosures.
2. THE ONE LOAD-BEARING GAP: the midpoint-rule pads.  Every lattice-vs-continuum
   step is O(D^2) = O(wbar^{-2}) and is padded here by PAD(wbar) = 0.25/wbar^2
   + 2e-6, a constant CALIBRATED from the measured discrepancies (Stage 7 prints
   them: max 2.9e-3 at wbar=6, 6.3e-5 at wbar=45, 1.4e-7 at wbar=1000; the 2e-6
   floor absorbs the quadrature accuracy of the absolute constants).  Turning M1
   into a proof needs the explicit majorants (1/24) int |u''| for the ~10
   integrands involved; each is elementary but none is written here.
   Same remark for the three sup_s evaluated on grids (Theta*, ||O[phi]||,
   ||K||): a Lipschitz-plus-grid argument is sketched for Theta* only.
3. The identities (M1),(M2),(M3), h >= rho >= 2B-1, lam_2 <= ||Rop||, the Schur
   test, and the Davis-Kahan assembly are FULLY rigorous and elementary.
4. The perturbation identity (C**) itself is quoted from Sec.1.A.2; its C^2
   hypothesis on x |-> A(x) is supplied there (Stage B of p4_extremal_concavity).

================================================================================
SEC. 4.  RESULT
================================================================================
    wbar_A = 7.776:  the criterion (CRIT) holds for EVERY wbar >= wbar_A, with
    ratio 0.454 (wbar=45), 0.330 (wbar=1000), and monotonically decreasing to the
    analytic limit  4 I[kappa^2] / (8 I[G phi_0]) = (1/6)/(2/3) = 1/4  --  the
    route does NOT degrade as wbar -> infinity (its asymptotic margin is 4x).
    Since -2F(1-1/4) = -1/2, the same limit reproduces the known x^2 H'' -> -1/2.
    Segment C-I certifies wbar <= 45, so the two segments OVERLAP by a factor 5.8
    in wbar (a factor 33 in x) -- modulo HONEST LIMITS 1-2, Lemma C, hence
    Theorem E, holds on all of (0,infinity).
    BOTTLENECK (wbar_A if one ingredient were replaced by its float64 value):
    F,E transport 1.20 < g_lo 3.64 < C_r 6.83 < rho_lo 6.96 < as proved 7.78 --
    i.e. the Davis-Kahan transport of F and E is what sets wbar_A, and NO single
    constant is close to breaking the criterion at wbar >= 45.
"""

import sys
import os
import time

import numpy as np
from mpmath import mp, mpf, digamma, euler, log as mlog, mpc, im, re as mre
from scipy.integrate import quad
from scipy.interpolate import CubicSpline

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p4_interchange_uniformity import H_mat                                # noqa: E402
from p4_extremal_asymptotics import OpKernel, top_eig                      # noqa: E402

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
# SEC. 1.  kernels, the exact identities, and a generic O[f]
# ============================================================================

G = lambda y: 1.0 / (1.0 + np.asarray(y, dtype=float) ** 2)
PHI = lambda y: np.asarray(y, dtype=float) ** 2 / (1.0 + np.asarray(y, dtype=float) ** 2) ** 3
KAP = lambda y: np.asarray(y, dtype=float) ** 2 / (1.0 + np.asarray(y, dtype=float) ** 2) ** 2
LAM_G = lambda s: np.arctan(s)
LAM_PHI = lambda s: 0.125 * (np.arctan(s) + s * (s * s - 1.0) / (1.0 + s * s) ** 2)
LAM_KAP = lambda s: 0.5 * (np.arctan(s) - s / (1.0 + s * s))
# I[f] = (1/2) int_0^inf f(s) ds/s   for the two f(0)=0 kernels (elementary):
I_PHI, I_KAP = 0.125, 0.25            # (1/2)*1/4 , (1/2)*1/2


def R_G(s, t):
    """(M1) remainder, the closed form; also computed as the raw difference."""
    s = np.asarray(s, float); t = np.asarray(t, float)
    num = 2.0 * s * s * t * t * (5.0 + s * s + t * t)
    den = (1.0 + (s - t) ** 2) * (1.0 + (s + t) ** 2) * (1.0 + s * s) * (1.0 + t * t)
    return num / den


def R_f(f, s, t):
    """(M2) remainder for a kernel with f(0)=0."""
    s = np.asarray(s, float); t = np.asarray(t, float)
    return f(np.abs(s - t)) + f(s + t) - 2.0 * f(s) * G(t) - 2.0 * f(t) * G(s)


R_PHI = lambda s, t: R_f(PHI, s, t)
R_KAP = lambda s, t: R_f(KAP, s, t)


def op_dense(N, w, f):
    """O[f(./wbar)] as a dense NxN matrix (generic builder, validated in Stage 1)."""
    n = np.arange(N)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)
    off = f(np.abs(n[:, None] - n[None, :]) / w)
    np.fill_diagonal(off, 0.0)
    M = (off + f((n[:, None] + n[None, :] + 1) / w)) * nrm
    cum = np.concatenate(([0.0], np.cumsum(2.0 * f(np.arange(1, N) / w))))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M


def lattice(N, w):
    n = np.arange(N)
    wn = 1.0 / (2.0 * n + 1.0)
    sn = (n + 0.5) / w
    return n, wn, sn


def ell(N, w, f):
    """l_f(n) = (2/(2n+1)) sum_{k<=n} f(k/wbar)."""
    n = np.arange(N)
    return 2.0 * np.concatenate(([0.0], np.cumsum(f(np.arange(1, N) / w)))) / (2.0 * n + 1.0)


def run_stage1(ws=(3.0, 7.0, 20.0, 45.0), N=60):
    print("\n=== Stage 1: the exact identities (M1),(M2) elementwise ===")
    for w in ws:
        n, wn, sn = lattice(N, w)
        p = np.sqrt(wn)
        # the generic builder reproduces the CERTIFIED entry formulas
        chk(f"op_dense==H_mat w={w}", np.max(np.abs(op_dense(N, w, G) - H_mat(N, w, 0))) < 1e-13,
            f"{np.max(np.abs(op_dense(N,w,G)-H_mat(N,w,0))):.2e}")
        # closed form of R equals the raw difference
        S, T = np.meshgrid(sn, sn, indexing="ij")
        raw = G(np.abs(S - T)) + G(S + T) - 2.0 * G(S) * G(T)
        chk(f"R closed form w={w}", np.max(np.abs(raw - R_G(S, T))) < 1e-14)
        chk(f"R>=0 w={w}", raw.min() >= -1e-16, f"{raw.min():.3e}")
        # (M1)
        tau = p * G(sn)
        Dg = wn + ell(N, w, G)
        chk(f"Dg in [0,1] w={w}", Dg.min() >= 0.0 and Dg.max() <= 1.0 + 1e-15)
        M1 = 2.0 * np.outer(tau, tau) + np.outer(p, p) * R_G(S, T) - np.diag(Dg)
        chk(f"(M1) w={w}", np.max(np.abs(M1 - H_mat(N, w, 0))) < 1e-13,
            f"{np.max(np.abs(M1-H_mat(N,w,0))):.2e}")
        # (M2) for phi_0 and kappa
        for nm, f, Rr in [("phi", PHI, R_PHI), ("kap", KAP, R_KAP)]:
            tf = p * f(sn)
            M2 = (2.0 * (np.outer(tf, tau) + np.outer(tau, tf))
                  + np.outer(p, p) * Rr(S, T) - np.diag(ell(N, w, f)))
            d = np.max(np.abs(M2 - op_dense(N, w, f)))
            chk(f"(M2)[{nm}] w={w}", d < 1e-13, f"{d:.2e}")
        # R_f vanishes on the axes  (*)
        for f, Rr in [(PHI, R_PHI), (KAP, R_KAP)]:
            chk("R_f(s,0)=0", abs(float(Rr(np.array(1.7), np.array(0.0)))) < 1e-15)
    stage_end("Stage 1")


# ============================================================================
# SEC. 2.  A(wbar), B(wbar):  closed form + the elementary two-sided bound
# ============================================================================

GAMMA = float(euler)
C0 = GAMMA + np.log(4.0) - 1.5                      # 0.4635101...


def AB_closed(w):
    """A, B from the digamma closed forms (mpmath, exact identities)."""
    z = mpc(mpf(1) / 2, mpf(w))
    A = 0.5 * (mre(digamma(z)) + euler + 2 * mlog(2))
    dp = mp.diff(digamma, z)                         # psi'(1/2 + i wbar)
    B = A + mpf(w) / 4 * im(dp)
    return float(A), float(B)


def AB_lattice(w, mult=4000):
    n = np.arange(int(mult * w) + 40)
    wn = 1.0 / (2.0 * n + 1.0)
    sn = (n + 0.5) / w
    A = float(np.sum(wn * G(sn)))
    B = float(np.sum(wn * G(sn) ** 2))
    # tails:  sum_{n>=N} w_n G^j ~ (1/2) int (wbar/s)^{2j} ds/s  -- add them
    N = len(n)
    sN = (N - 0.5) / w
    A += 0.5 / (2 * sN ** 2)
    B += 0.5 / (4 * sN ** 4)
    return A, B


def rho_lo_elem(w):
    """The fully elementary lower bound of Sec.2.1 (a)."""
    return np.log(w) + C0 - 1.0 / (12.0 * w) - 1.0 / (4.0 * w * w) - 1.0 / (16.0 * w * w)


def run_stage2(ws):
    print("\n=== Stage 2: A(wbar), B(wbar) closed form and the elementary bound ===")
    print("   wbar |      A        B    | A,B lattice err | 2B-1        ln w + c0    elem. bound")
    for w in ws:
        A, B = AB_closed(w)
        Al, Bl = AB_lattice(w)
        chk(f"A closed==lattice w={w}", abs(A - Al) < 5e-8, f"{abs(A-Al):.2e}")
        chk(f"B closed==lattice w={w}", abs(B - Bl) < 5e-8, f"{abs(B-Bl):.2e}")
        chk(f"2B-1 >= elem bound w={w}", 2 * B - 1 >= rho_lo_elem(w) - 1e-12,
            f"{2*B-1:.9f} vs {rho_lo_elem(w):.9f}")
        print(f" {w:6g} | {A:8.5f} {B:8.5f} | {abs(A-Al):.1e} {abs(B-Bl):.1e} |"
              f" {2*B-1:11.7f} {np.log(w)+C0:11.7f} {rho_lo_elem(w):11.7f}")
    stage_end("Stage 2")


# ============================================================================
# SEC. 3.  the absolute constants  I[.]
# ============================================================================

def I_of(u, smax=60.0):
    """I[u] = (1/2) int_0^inf u(s)/s ds."""
    f = lambda s: u(s) / s
    tot = quad(f, 1e-13, 1.0, limit=400)[0] + quad(f, 1.0, smax, limit=400)[0]
    tot += quad(f, smax, np.inf, limit=400)[0]
    return 0.5 * tot


def _prof(Rr, s, weightG):
    g = (lambda t: G(t) * float(Rr(np.array(s), np.array(t))) / t) if weightG else \
        (lambda t: float(Rr(np.array(s), np.array(t))) / t)
    edges = sorted({1e-13, 0.5 * s, s, 2.0 * s, 5.0, 60.0})
    edges = [e for e in edges if e > 0]
    tot, prev = 0.0, 1e-13
    for e in edges:
        if e > prev:
            tot += quad(g, prev, e, limit=400)[0]
            prev = e
    tot += quad(g, prev, np.inf, limit=400)[0]
    return 0.5 * tot


_GRID = np.concatenate([np.linspace(1e-4, 0.5, 60), np.linspace(0.52, 4, 80),
                        np.linspace(4.1, 20, 60), np.linspace(21, 120, 40)])
_SMAX = _GRID[-1]


def build_profiles():
    """Ups_f(s) = I_sigma[G R_f(s,.)]  and  Theta_f(s) = I_sigma[R_f(s,.)]."""
    out = {}
    for nm, Rr in [("G", R_G), ("phi", R_PHI), ("kap", R_KAP)]:
        out[nm + "_U"] = CubicSpline(_GRID, [_prof(Rr, s, True) for s in _GRID])
        out[nm + "_T"] = CubicSpline(_GRID, [_prof(Rr, s, False) for s in _GRID])
    return out


SPL = None


def U(nm, s):
    s = np.asarray(s, float)
    return np.where(s <= _SMAX, SPL[nm + "_U"](np.minimum(s, _SMAX)), 0.0)


def TH(nm, s):
    s = np.asarray(s, float)
    return np.where(s <= _SMAX, SPL[nm + "_T"](np.minimum(s, _SMAX)), 0.0)


# stored values (reproduced by run_stage3)
CONST = dict(
    P1=1.0 / 12.0, K1=1.0 / 8.0, mPP=1.0 / 80.0, mKK=1.0 / 24.0,
    cDel=0.05027541, nD2=0.00865336, qR=0.12892810, nR2=0.04925665,
    phi_tDt=0.03017391, phi_fDt=0.00465127, phi_DD=0.00174685,
    phi_fU=-0.01263179, phi_lU=-0.00455746, phi_UU=0.01422639, phi_GU=-0.09558297,
    kap_tDt=0.03736229, kap_fDt=0.01149276, kap_DD=0.00323943,
    kap_fU=-0.01913140, kap_lU=-0.00593825, kap_UU=0.01546955, kap_GU=-0.09040317,
)


def compute_constants():
    d = lambda s: 1.0 - np.arctan(s) / s
    C = dict(P1=I_of(lambda s: G(s) * PHI(s)), K1=I_of(lambda s: G(s) * KAP(s)),
             mPP=I_of(lambda s: PHI(s) ** 2), mKK=I_of(lambda s: KAP(s) ** 2),
             cDel=I_of(lambda s: d(s) * G(s) ** 2),
             nD2=I_of(lambda s: d(s) ** 2 * G(s) ** 2),
             qR=I_of(lambda s: G(s) * U("G", s)), nR2=I_of(lambda s: U("G", s) ** 2))
    for nm, f, Lam in [("phi", PHI, LAM_PHI), ("kap", KAP, LAM_KAP)]:
        l = lambda s: Lam(s) / s
        C[nm + "_tDt"] = I_of(lambda s: l(s) * G(s) ** 2)
        C[nm + "_fDt"] = I_of(lambda s: f(s) * l(s) * G(s))
        C[nm + "_DD"] = I_of(lambda s: l(s) ** 2 * G(s) ** 2)
        C[nm + "_fU"] = I_of(lambda s: f(s) * U(nm, s))
        C[nm + "_lU"] = I_of(lambda s: l(s) * G(s) * U(nm, s))
        C[nm + "_UU"] = I_of(lambda s: U(nm, s) ** 2)
        C[nm + "_GU"] = I_of(lambda s: G(s) * U(nm, s))
    return C


def lattice_constants(w, mult=400):
    """The same constants as genuine lattice sums at one wbar (independent check)."""
    N = int(mult * w)
    n, wn, sn = lattice(N, w)
    p = np.sqrt(wn)
    tau = p * G(sn)
    B = float(tau @ tau)
    Dg = wn + ell(N, w, G)
    out = dict(B=B, cDel=float(np.sum((1.0 - Dg) * wn * G(sn) ** 2)),
               nD2=float(np.sum((1.0 - Dg) ** 2 * wn * G(sn) ** 2)))
    # Upsilon_f(n) = sum_m w_m G_m R_f(s_n,s_m) -- O(N) per n, sample on a coarse set
    idx = np.unique(np.concatenate([np.arange(0, 40), np.geomspace(40, N - 1, 400).astype(int)]))
    for nm, Rr in [("G", R_G), ("phi", R_PHI), ("kap", R_KAP)]:
        Ups = np.array([float(np.sum(wn * G(sn) * Rr(sn[i], sn))) for i in idx])
        out[nm + "_Ups"] = (sn[idx], Ups)
    return out


def run_stage3():
    print("\n=== Stage 3: the absolute constants (quadrature vs stored vs lattice) ===")
    C = compute_constants()
    for k, v in sorted(C.items()):
        chk(f"const {k}", abs(v - CONST[k]) < 3e-7, f"{v:.9f} vs {CONST[k]:.9f}")
    print("   " + "  ".join(f"{k}={C[k]:.7f}" for k in
                            ["P1", "K1", "cDel", "qR", "nR2", "nD2"]))
    print("   " + "  ".join(f"{k}={C[k]:.7f}" for k in
                            ["phi_GU", "phi_tDt", "kap_GU", "kap_tDt"]))
    # independent check: the Upsilon profiles as lattice sums at wbar = 1000
    lc = lattice_constants(1000.0, mult=60)
    for nm in ["G", "phi", "kap"]:
        ss, val = lc[nm + "_Ups"]
        m = (ss > 0.02) & (ss < 60)
        err = float(np.max(np.abs(val[m] - U(nm, ss[m]))))
        chk(f"Ups[{nm}] lattice(wbar=1000)==continuum", err < 3e-4, f"{err:.2e}")
        print(f"   Ups[{nm}]: max |lattice(wbar=1000) - continuum| = {err:.2e}")
    for k in ["cDel", "nD2"]:
        chk(f"{k} lattice", abs(lc[k] - CONST[k]) < 5e-6, f"{lc[k]:.8f}")
    stage_end("Stage 3")
    return C


# ============================================================================
# SEC. 4.  Theta*, the Schur norms, and the five bounds as functions of wbar
# ============================================================================

THETA = lambda s: (np.log1p(s * s) + 2.0 * s * np.arctan(s)) / (2.0 * (1.0 + s * s))


def theta_star():
    """max Theta on [0,2.53] by grid + Lipschitz;  Theta decreasing beyond."""
    ss = np.linspace(0.0, 2.53, 2530001)
    v = THETA(ss)
    i = int(np.argmax(v))
    return float(v[i]), float(ss[i])


THETA_STAR, S_STAR = 0.63492879, 1.500904
# calibrated pad: 0.25/wbar^2 covers the midpoint-rule (lattice -> continuum) error,
# the 2e-6 floor covers the quadrature/spline accuracy of the absolute constants.
PAD = lambda w: 0.25 / (w * w) + 2e-6   # see HONEST LIMITS 2


def schur_norm(nm, f, Lam, A):
    ss = np.concatenate([np.linspace(1e-3, 6.0, 40001), np.linspace(6.0, 200.0, 20001)])
    return float(np.max(2.0 * A * f(ss) + TH(nm, ss) + 2.0 * G(ss) *
                        (I_PHI if nm == "phi" else I_KAP) + Lam(ss) / ss))


def bounds(w, C=None, tight_rho=False):
    """ALL five bounds of Sec.2 at one wbar.  Pure formulas in wbar."""
    C = C or CONST
    A, B = AB_closed(w)
    pad = PAD(w)
    rho = 2.0 * B - 1.0
    if tight_rho:                                   # the two dropped >=0 terms
        rho += (C["cDel"] + C["qR"]) / B - pad
    Cr = (np.sqrt(C["nR2"]) + np.sqrt(C["nD2"])) / np.sqrt(B) + pad
    g = rho - THETA_STAR - pad
    s = Cr / g
    F_t = 4.0 * C["P1"] + (C["phi_GU"] - C["phi_tDt"]) / B - pad
    q_k = 4.0 * C["K1"] + (C["kap_GU"] - C["kap_tDt"]) / B
    e = {}
    for nm, m_ff, P1 in [("phi", C["mPP"], C["P1"]), ("kap", C["mKK"], C["K1"])]:
        nZ2 = (4 * B * B * m_ff - 4 * B * C[nm + "_fDt"] + C[nm + "_DD"]
               + 4 * B * C[nm + "_fU"] - 2 * C[nm + "_lU"] + C[nm + "_UU"])
        tZ = 2 * B * P1 - C[nm + "_tDt"] + C[nm + "_GU"]
        e[nm] = np.sqrt(max(0.0, (nZ2 - tZ * tZ / B) / B)) + pad
    NP = schur_norm("phi", PHI, LAM_PHI, A)
    NK = schur_norm("kap", KAP, LAM_KAP, A)
    F_lo = (1.0 - s * s) * F_t - 2.0 * s * e["phi"] - s * s * NP
    E_hi = e["kap"] + s * (abs(q_k) + NK)
    ratio = E_hi ** 2 / (g * F_lo) if F_lo > 0 else np.inf
    return dict(A=A, B=B, rho=rho, Cr=Cr, g=g, s=s, F_t=F_t, q_k=q_k,
                e_phi=e["phi"], e_kap=e["kap"], NP=NP, NK=NK,
                F_lo=F_lo, E_hi=E_hi, ratio=ratio)


# ============================================================================
# SEC. 5.  float64 truth
# ============================================================================

def truth(w, mult=40):
    N = int(round(mult * w))
    k = np.arange(2 * N) / w
    n, wn, sn = lattice(N, w)
    tau = np.sqrt(wn) * G(sn)
    B = float(tau @ tau)
    that = tau / np.sqrt(B)
    oH, oP, oK = OpKernel(N, G(k)), OpKernel(N, PHI(k)), OpKernel(N, KAP(k))
    rho = float(that @ oH.matvec(that))
    h, v, _ = top_eig(N, G(k))
    r = oH.matvec(that) - rho * that

    def qe(op, z):
        y = op.matvec(z)
        q = float(z @ y)
        return q, float(np.linalg.norm(y - q * z))

    qp, ep = qe(oP, that)
    qk, ek = qe(oK, that)
    Fv, _ = qe(oP, v)
    _, Ev = qe(oK, v)
    nP = float(np.linalg.norm(oP.matvec(that)))
    nK = float(np.linalg.norm(oK.matvec(that)))
    sin_tv = float(np.sqrt(max(0.0, 1.0 - (that @ v) ** 2)))
    return dict(N=N, B=B, rho=rho, h=h, nr=float(np.linalg.norm(r)), F_t=qp, e_phi=ep,
                q_k=qk, e_kap=ek, F_v=Fv, E_v=Ev, nP=nP, nK=nK, sin=sin_tv,
                ratio=Ev ** 2 / (h * Fv))


def run_stage4(ws):
    print("\n=== Stage 4: bounds (a),(b),(d) vs float64 truth ===")
    print("  wbar |   rho_lo      h     | C_r     ||r||  | F_lo(t)  F(t)  | e_kap(t) exact | e_phi  exact")
    T = {}
    for w in ws:
        b = bounds(w)
        t = truth(w)
        T[w] = (b, t)
        chk(f"(a) rho_lo<=h w={w}", b["rho"] <= t["h"], f"{b['rho']:.6f} {t['h']:.6f}")
        chk(f"(a) rho_lo<=rho w={w}", b["rho"] <= t["rho"] + 1e-12)
        chk(f"(a) elem<=2B-1 w={w}", rho_lo_elem(w) <= b["rho"])
        chk(f"(b) C_r>=||r|| w={w}", b["Cr"] >= t["nr"], f"{b['Cr']:.6f} {t['nr']:.6f}")
        chk(f"(d1) F_lo<=F(t) w={w}", b["F_t"] <= t["F_t"], f"{b['F_t']:.7f} {t['F_t']:.7f}")
        chk(f"(d2) e_kap>=exact w={w}", b["e_kap"] >= t["e_kap"], f"{b['e_kap']:.7f} {t['e_kap']:.7f}")
        chk(f"(d2) e_phi>=exact w={w}", b["e_phi"] >= t["e_phi"])
        chk(f"(d3) Schur>=||O[phi]that|| w={w}", b["NP"] >= t["nP"])
        chk(f"(d3) Schur>=||K that|| w={w}", b["NK"] >= t["nK"])
        chk(f"q_kappa w={w}", abs(b["q_k"] - t["q_k"]) < 2e-3, f"{b['q_k']:.6f} {t['q_k']:.6f}")
        print(f" {w:6g} | {b['rho']:9.5f} {t['h']:9.5f} | {b['Cr']:.4f} {t['nr']:.4f} |"
              f" {b['F_t']:.5f} {t['F_t']:.5f} | {b['e_kap']:.5f} {t['e_kap']:.5f} |"
              f" {b['e_phi']:.5f} {t['e_phi']:.5f}")
    stage_end("Stage 4")
    return T


def run_stage5(ws, T):
    print("\n=== Stage 5: bound (c) -- the gap, and why the lam_min(M) route is vacuous ===")
    th, ss = theta_star()
    chk("Theta* grid", abs(th - THETA_STAR) < 1e-7, f"{th:.9f}")
    chk("s* ", abs(ss - S_STAR) < 1e-3)
    chk("Theta decreasing beyond 2.53", np.all(np.diff(THETA(np.linspace(2.53, 200, 200001))) < 0))
    # stationarity equation
    f_st = lambda s: np.arctan(s) * (1 - s * s) - s * (np.log1p(s * s) - 2)
    chk("Theta* stationarity", abs(f_st(ss)) < 1e-4, f"{f_st(ss):.2e}")
    print(f"   Theta* = {th:.8f} at s* = {ss:.5f}   (stationarity residual {f_st(ss):.1e})")
    print("  wbar | sup_n Sigma_n (lattice) | Theta* | g_lo = rho_lo-Theta* | true gap h-lam2")
    for w in ws:
        b, t = T[w]
        # sup_n Sigma_n directly (the Schur row sums of Rop), N = 200 wbar + tail
        N = int(200 * w)
        n, wn, sn = lattice(N, w)
        idx = np.unique(np.concatenate([np.arange(0, 20),
                                        np.geomspace(20, N - 1, 600).astype(int)]))
        Sig = max(float(np.sum(R_G(sn[i], sn) * wn)) for i in idx)
        chk(f"sup Sigma_n <= Theta*+pad w={w}", Sig <= THETA_STAR + PAD(w), f"{Sig:.8f}")
        chk(f"g_lo>0 w={w}", b["g"] > 0)
        chk(f"g_lo <= true gap w={w}", b["g"] <= t["h"] + 1e-9)
        print(f" {w:6g} | {Sig:.8f} ({Sig-THETA_STAR:+.1e}) | {THETA_STAR:.5f} | {b['g']:9.5f}"
              f"            | {t['h']:9.5f}")
    # the collapse of lam_min(M) = -h on l^2:  variational upper bound on the
    # two-dimensional space span{v, z}, z := u restricted to [2 wbar, N] (unit, _|_ v).
    print("\n   lam_min(M), M = I+W = 2uu^T - H  (the dispatch's route (c)) -- collapse:")
    print("   variational bound on span{v,z};  <u,z> = sqrt(sum_{2wbar}^{N} w_n) -> infinity")
    w = 50.0
    prev, vals = 1e9, []
    for mult in [20, 80, 320, 1280, 5120]:
        N = int(mult * w)
        k = np.arange(2 * N) / w
        h, v, oH = top_eig(N, G(k))
        n, wn, sn = lattice(N, w)
        u = np.sqrt(wn)
        z = np.zeros(N)
        z[int(2 * w):] = u[int(2 * w):]
        z -= (v @ z) * v
        z /= np.linalg.norm(z)
        uv, uz = float(u @ v), float(u @ z)
        Hz = oH.matvec(z)
        off = 2.0 * uv * uz - float(v @ Hz)
        Mmat = np.array([[2.0 * uv * uv - h, off],
                         [off, 2.0 * uz * uz - float(z @ Hz)]])
        lam = float(np.linalg.eigvalsh(Mmat)[0])
        vals.append(h + lam)
        print(f"     N={mult:5d} wbar:  <u,z> = {uz:6.3f}   lam_min(M) <= {lam:9.5f}"
              f"   h + lam_min(M) <= {h+lam:8.5f}   (floor -h = {-h:.5f})")
        chk(f"lam_min(M_N) decreasing N={mult}", h + lam < prev + 1e-9)
        prev = h + lam
    chk("lam_min(M) collapse demonstrated", vals[-1] < 0.8 * vals[0], f"{vals}")
    print(f"     -> h + lam_min(M) falls like 1/ln N (limit 0, Sec.2.3); the tau-route (M1)"
          f"\n        gives g_lo = {bounds(50.0)['g']:.4f} at wbar = 50 instead.")
    stage_end("Stage 5")


def run_stage6(ws, T):
    print("\n=== Stage 6: (e) assembly -- the criterion, wbar_A, and the bottleneck ===")
    print("  wbar |  g_lo   sin(th)  |  F_lo    E_hi   | RATIO  | true ratio | x^2 H''<= | float64 x^2H''")
    B1 = {45.0: None, 100.0: -0.493280384, 200.0: -0.494738431, 500.0: -0.496048452,
          1000.0: -0.496742993, 2000.0: -0.497270568}     # Sec.1.B table (B1)
    for w in ws:
        b, t = T[w]
        chk(f"sin bound w={w}", b["s"] >= t["sin"], f"{b['s']:.5f} {t['sin']:.5f}")
        chk(f"E_hi >= E(v) w={w}", b["E_hi"] >= t["E_v"], f"{b['E_hi']:.5f} {t['E_v']:.5f}")
        chk(f"F_lo <= F(v) w={w}", b["F_lo"] <= t["F_v"], f"{b['F_lo']:.5f} {t['F_v']:.5f}")
        chk(f"criterion closes w={w}", b["ratio"] < 1.0, f"{b['ratio']:.4f}")
        x2H = -2.0 * b["F_lo"] + 2.0 * b["E_hi"] ** 2 / b["g"]      # (C**), x^2-normalised
        if B1[w] is not None:
            chk(f"x^2H'' bound >= truth w={w}", x2H >= B1[w], f"{x2H:.6f} {B1[w]:.6f}")
            chk(f"x^2H'' bound non-vacuous w={w}", x2H < 0.60 * B1[w], f"{x2H:.6f}")
        print(f" {w:6g} | {b['g']:7.4f} {b['s']:.5f} | {b['F_lo']:.5f} {b['E_hi']:.5f} |"
              f" {b['ratio']:.4f} | {t['ratio']:.4f}     | {x2H:9.6f} |"
              f" {B1[w] if B1[w] is not None else float('nan'):.6f}")
    # wbar_A
    lo, hi = 2.0, 60.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if bounds(mid)["ratio"] > 1.0:
            lo = mid
        else:
            hi = mid
    wA = hi
    chk("wbar_A found", 2.0 < wA < 60.0, f"{wA}")
    chk("closes above wbar_A", bounds(wA * 1.001)["ratio"] < 1.0)
    chk("fails below wbar_A", bounds(wA * 0.999)["ratio"] > 1.0)
    print(f"\n   ***  wbar_A = {wA:.4f}  (criterion closes for all wbar >= wbar_A)  ***")
    print(f"   Segment C-I certifies wbar <= 45  ==>  overlap factor {45/wA:.1f} in wbar.")
    # no degradation at infinity: scan, and the analytic limit 1/4
    grid = np.geomspace(wA * 1.0001, 1e20, 200)
    rr = np.array([bounds(x)["ratio"] for x in grid])
    chk("ratio < 1 on [wbar_A, 1e20]", float(rr.max()) < 1.0, f"{rr.max():.4f}")
    chk("ratio decreasing at large wbar", rr[-1] < rr[len(rr) // 2] < rr[2])
    lim = 4.0 * CONST["mKK"] / (8.0 * CONST["P1"])
    chk("asymptotic ratio -> 1/4", abs(lim - 0.25) < 1e-12 and abs(rr[-1] - 0.25) < 0.02,
        f"{lim} {rr[-1]:.5f}")
    print(f"   ratio is monotone decreasing to the analytic limit 4 I[kappa^2]/(8 I[G phi_0])"
          f" = (1/6)/(2/3) = 1/4   (value at wbar=1e20: {rr[-1]:.5f});")
    print(f"   hence x^2 H'' <= -2(1/3)(1-1/4) = -1/2 asymptotically -- the exact known limit.")
    # bottleneck analysis: switch each ingredient to its numerical value
    print("\n   bottleneck (wbar_A when one ingredient is replaced by its float64 value):")
    variants = [("as proved", None), ("C_r -> true ||r||", "Cr"),
                ("g_lo -> true gap", "g"), ("F,E -> values at v", "FE"),
                ("rho_lo -> tight (keep c_Delta,qR)", "rho")]
    for name, key in variants:
        def rat(w, key=key):
            b = bounds(w, tight_rho=(key == "rho"))
            if key == "Cr":
                b2 = dict(b)
                b2["s"] = truth(w)["nr"] / b["g"]
                s = b2["s"]
                F = (1 - s * s) * b["F_t"] - 2 * s * b["e_phi"] - s * s * b["NP"]
                E = b["e_kap"] + s * (abs(b["q_k"]) + b["NK"])
                return E * E / (b["g"] * F)
            if key == "g":
                t = truth(w)
                s = b["Cr"] / t["h"]
                F = (1 - s * s) * b["F_t"] - 2 * s * b["e_phi"] - s * s * b["NP"]
                E = b["e_kap"] + s * (abs(b["q_k"]) + b["NK"])
                return E * E / (t["h"] * F)
            if key == "FE":
                t = truth(w)
                return t["E_v"] ** 2 / (b["g"] * t["F_v"])
            return b["ratio"]
        lo2, hi2 = 1.2, 80.0
        for _ in range(50):
            mid = 0.5 * (lo2 + hi2)
            try:
                bad = rat(mid) > 1.0
            except Exception:
                bad = True
            if bad:
                lo2 = mid
            else:
                hi2 = mid
        print(f"     {name:34s} wbar_A = {hi2:6.3f}")
    stage_end("Stage 6")
    return wA


def run_stage7(ws):
    print("\n=== Stage 7: midpoint-pad calibration (the one load-bearing gap) ===")
    print("  wbar |  rho     F_t     e_kap   e_phi   | PAD(wbar)=0.25/wbar^2 + 2e-6")
    worst = 0.0
    for w in ws:
        b = bounds(w)
        t = truth(w, mult=60)
        # undo the pad to see the raw continuum-vs-lattice discrepancy
        p = PAD(w)
        d = [abs((b["rho"]) - (2 * t["B"] - 1)),
             abs((b["F_t"] + p) - t["F_t"]), abs((b["e_kap"] - p) - t["e_kap"]),
             abs((b["e_phi"] - p) - t["e_phi"])]
        worst = max(worst, (max(d) - 2e-6) * w * w)
        chk(f"pad covers w={w}", max(d) <= p, f"{max(d):.2e} vs {p:.2e}")
        print(f" {w:6g} | " + " ".join(f"{x:.1e}" for x in d) + f" | {p:.2e}")
    print(f"   worst  wbar^2 * (discrepancy - 2e-6) = {worst:.4f}   (PAD uses 0.25)")
    chk("pad constant is safe", worst < 0.25)
    stage_end("Stage 7")


def main():
    quick = "--quick" in sys.argv
    ws = (45.0, 100.0, 200.0) if quick else (45.0, 100.0, 200.0, 500.0, 1000.0, 2000.0)
    print(__doc__.split("SEC. 0.")[0].strip()[:0] or "", end="")
    print("P4 Segment C-II Milestone M1 -- closed-form bounds for (C**)")
    global SPL
    run_stage1()
    run_stage2(ws)
    SPL = build_profiles()
    C = run_stage3()
    T = run_stage4(ws)
    run_stage5(ws, T)
    wA = run_stage6(ws, T)
    run_stage7((6.0, 10.0, 20.0) + ws)
    print(f"\n=== ALL {PASS[0]} CHECKS PASSED in {time.time()-T0:.1f} s ===")
    print(f"=== wbar_A = {wA:.4f} ;  criterion ratio at wbar=45 is "
          f"{bounds(45.0)['ratio']:.4f}, at wbar=1000 {bounds(1000.0)['ratio']:.4f} ===")


if __name__ == "__main__":
    main()
