r"""
P4 . EINSTEIN EXTREMAL -- Segment C-I:  H(x) := h(sqrt x) = lam_1(H(sqrt x)) is
STRICTLY CONCAVE on (0, X], certified cell by cell (mpmath dps=30).

RUN:  python theory/derivations/p4_extremal_concavity.py         (FAST default, ~35 s:
                                                                  stages A B S + the
                                                                  ladder up to wbar=3)
      python theory/derivations/p4_extremal_concavity.py EXT     (EXTENDED, ~10 min:
                                                                  the ladder to wbar=45)
      python theory/derivations/p4_extremal_concavity.py EXT 60  (extended, W_max=60)

SCOPE: this file owns Sec.1 Segment C-I of `theory/derivations/P4_einstein_extremal.md`
(Lemma C on (0, X]).  Segment C-II (the analytic large-x half) is NOT here.

NOTHING FROM THE CERTIFIED BASE IS RE-DERIVED: `H_mat`/`H_matrix`/`G_matrix`,
`spec_point`/`spec_seg`/`envelope`, `schur_coupling`, `schur_tailblock`,
`schur_full_seg`, `alpha_certificate_H`, `mu2_cell_lb`, `ROUND_PAD` are
IMPORTED from `p4_interchange_uniformity`; `schur_full_seg_w0` from
`p4_interval_extension`; `coupling_restricted` (the restricted head-tail
coupling envelope, Stage M Sec.M.2) from `p4_conj2_interval`; `alpha_G2` (the
lam_2(G(2)) Cholesky-deflation certificate) and `schur_rowsum_bound` from
`p4_conj2_small_w`.  No base file is modified.

================================================================================
SEC. 1.  THE INEQUALITY (C*) AND ITS CANCELLATION-FREE FORM (C**)
================================================================================

Write  x := wbar^2,  A(x) := H(sqrt x)  (KAY's operator at wbar = sqrt x), and
H(x) := lam_1(A(x)).

(1.1) SECOND-ORDER PERTURBATION.  Let x |-> A(x) be a C^2 family of bounded
self-adjoint operators on l^2 on an open interval J, and suppose that for x in J
the top of the spectrum lam_1(x) is a SIMPLE ISOLATED eigenvalue with unit
eigenvector v(x) and gap g(x) = lam_1(x) - sup(spec A(x) \ {lam_1(x)}) > 0.
Then (Kato II.Sec.5 / Rellich) lam_1, v are C^2 resp. C^1 on J with the
normalisation <v,v'> = 0, and

    lam_1'  = <v, A_x v>,
    lam_1'' = <v, A_xx v> + 2 <A_x v, (lam_1 - A)^+ A_x v>,                 (P2)

where (lam_1-A)^+ is the reduced resolvent on {v}^perp.  Since
spec(A|_{v^perp}) <= lam_1 - g, the operator (lam_1-A)^+ is positive with norm
<= 1/g, so the second term of (P2) is >= 0 and

    lam_1''(x) <= <v, A_xx v> + (2/g) ||(I-P) A_x v||^2,   P := |v><v|.      (P3)

(In eigen-coordinates (P2) reads 2 sum_{j>=2} |<v_j,A_x v>|^2/(lam_1-lam_j),
the form quoted in P4_einstein_extremal.md Sec.1; (P3) is the same with
lam_1-lam_j >= g.  We use the resolvent form so that continuous spectrum is
covered without change.)

(1.2) CHAIN RULE.  With wbar = sqrt x,
    A_x  = H'(wbar)/(2 wbar),
    A_xx = (H''(wbar) - H'(wbar)/wbar)/(4 wbar^2),
so multiplying (P3) by 4 wbar^2 > 0 gives the note's

    4 wbar^2 H''(x)  <=  <v,H''v> - <v,H'v>/wbar + 2 ||(I-P) H' v||^2 / g.   (C*)

(1.3) WHY (C*) MUST NOT BE CERTIFIED TERM BY TERM.  As wbar -> 0,
<v,H''v> -> 2 frak_g(2) and <v,H'v>/wbar -> 2 frak_g(2) while their difference
is O(wbar^2): a term-by-term certificate has an error that does not vanish and
(C*) would close only for wbar bounded away from 0 -- the same "both sides
cancel" problem Segment I of P4_conjecture2.md met.  The cure is the same:
an EXACT identity.  With the kernel-generic operator
    O[f] := -L1[f] + L2[f] + L3[f],
    L1[f]_nn = (1/(2n+1)) sum_{k=1}^{n} 2 f(k),
    L2[f]_nm = (1-delta_nm) f(|n-m|)/sqrt((2n+1)(2m+1)),
    L3[f]_nm = f(n+m+1)/sqrt((2n+1)(2m+1)),
one has H(wbar) = O[f_wbar], f_wbar(k) = wbar^2/(k^2+wbar^2), and ELEMENTWISE

    [[k]]''(wbar) - [[k]]'(wbar)/wbar
        = 2k^2(k^2-3wbar^2)/(k^2+wbar^2)^3 - 2k^2/(k^2+wbar^2)^2
        = -8 wbar^2 k^2/(k^2+wbar^2)^3,

hence, with  phi_x(k) := k^2/(k^2+x)^3  (x = wbar^2),

    **  A_xx = -2 O[phi_x]  **                                              (I1)

-- one single POSITIVE kernel, no cancellation, and phi_0(k) = 1/k^4 so
A_xx(0) = -2 G(4).  Likewise, from  k^2/(k^2+wbar^2)^2 = 1/(k^2+wbar^2)
- wbar^2/(k^2+wbar^2)^2  and  H(wbar)/wbar^2 = O[1/(k^2+wbar^2)],

    **  A_x = H(wbar)/wbar^2 - x O[psi_x],   psi_x(k) := 1/(k^2+x)^2  **     (I2)

and since  (I-P) H(wbar) v = 0  (v is the eigenvector!),

    **  (I-P) A_x v = -x (I-P) O[psi_x] v.  **                               (I3)

(I1)-(I3) turn (P3) into the form this file certifies:

    H''(x)  <=  -2 F(wbar) + 2 x^2 E(wbar)^2 / g(wbar),                     (C**)
        F(wbar) := <v, O[phi_x] v>,   E(wbar) := ||(I-P) O[psi_x] v||.

Both O[phi_x] and O[psi_x] are O(wbar^{-4}) at large wbar and O(1) at wbar -> 0;
NOTHING cancels.  (I1)/(I2) are verified elementwise against `H_mat(N,w,1)`,
`H_mat(N,w,2)` in Stage A.

(1.4) THE x -> 0 LIMIT (no wbar_min gap).  phi_0 = 1/k^4 gives F(0) = <v_0,G(4)v_0>
= c_2 (Segment I: c_2 in [1.19312762268, 1.19320906584]); psi_0 = 1/k^4 gives
E(0) = ||(I-P_0) G(4) v_0|| < infinity; and g(wbar) = wbar^2 g~(wbar) with
g~(0) = frak_g(2) - lam_2(G(2)) >= 1.30443, so the positive term of (C**) is
2 wbar^4 E^2/(wbar^2 g~) = 2 wbar^2 E^2/g~ -> 0.  Hence

    H''(0) = -2 c_2 = -2.3862...,   and the RHS of (C*) is  -8 c_2 x + O(x^2),

i.e. 4 wbar^2 H'' ~ -9.5450 wbar^2.  The first cell of the certificate reaches
x = 0 exactly (Stage S), with no wbar_min gap: it works in the RESCALED operator
h~(wbar) := H(wbar)/wbar^2 = G(2) - wbar^2 K(wbar),  K(wbar) = O[1/(k^2(k^2+wbar^2))]
(Segment I's identity (I)), whose top eigenvalue is simple and isolated
UNIFORMLY down to wbar = 0, unlike H(wbar) itself (whose whole spectrum
collapses to 0).

(1.5) C^2 IN OPERATOR NORM.  For wbar in [a,b] with a>0, the elementwise Taylor
remainder of order 3 is dominated by the order-3 segment envelope matrix
`spec_seg(a,b,3)`, whose Schur norm `schur_full_seg_w0(a,b,3)` is finite; hence
||H(w+d) - H(w) - d H'(w) - d^2 H''(w)/2|| <= (|d|^3/6) schur_full_seg_w0(a,b,3)
and ||H''(w+d)-H''(w)|| <= |d| schur_full_seg_w0(a,b,3): wbar |-> H(wbar) is C^2
in operator norm on every [a,b] with 0<a<b<inf, so x |-> A(x) is C^2 on
(0,inf).  Near x=0 use instead A(x) = x G(2) - x^2 K~(x),
K~(x) := O[1/(k^2(k^2+x))], whose kernel is real-analytic in x on (-1,inf)
with an x-uniform Schur majorant: A is ANALYTIC in x on a neighbourhood of 0
(Stage B reports both).

================================================================================
SEC. 2.  THE PER-CELL CERTIFICATE
================================================================================

Fix a cell [a,b], 0 <= a <= b.  All kernels used are non-increasing in wbar
pointwise in k:
    phi_wbar(k) = k^2/(k^2+wbar^2)^3,  psi_wbar(k) = 1/(k^2+wbar^2)^2,
    kappa_wbar(k) = 1/(k^2(k^2+wbar^2)),
and all are dominated on the cell by the NON-INCREASING (in k) majorant
    m_a(k) := 1/(k^2+a^2)^2       (phi <= psi <= m_a,  kappa <= m_a).

(2.1) QUADRATIC FORMS ARE EXACT.  For a finite-support trial vector vhat
(support Nv, u_n := vhat_n/sqrt(2n+1)),
    <vhat,O[f]vhat>/||vhat||^2 = sum_{k>=1} f(k) c_k,
    c_k = T_k + U_k - P_k,  T_k = 2 sum_n u_n u_{n+k},
    U_k = sum_{n+m+1=k} u_n u_m,  P_k = 2 sum_{n>=k} u_n^2,     (all /||vhat||^2)
-- an identity, verified against the direct double sum in Stage A.  Since
f_wbar(k) is monotone in wbar, splitting c_k into c^+ - c^- gives the EXACT
cell-uniform envelope
    qhat_lo(a,b) = sum_k [c_k^+ f_b(k) - c_k^- f_a(k)]  <=  qhat(wbar)
                <= sum_k [c_k^+ f_a(k) - c_k^- f_b(k)] = qhat_hi(a,b).
(Empirically c_k >= 0 for every cell here, so qhat_lo = qhat(b) is ATTAINED --
there is no cell-width penalty at all in this term; the sign split is used
anyway, it costs nothing and needs no assumption.)

(2.2) THE ANGLE.  theta := <vhat,H(p)vhat> <= lam_1(H(p)) (Rayleigh, p the cell
midpoint); r := H(p)vhat - theta vhat with the m >= Nr rows bounded by
`coupling_restricted` (the tail-coupling term is MANDATORY, cf.
P4D_coulomb_operator.md Sec.2.4b); alpha >= lam_2 from the coarse grid.  Then
Davis-Kahan sin(angle(vhat,v(p))) <= ||r||/(theta-alpha) =: s_p.  For the whole
cell, v(.) is a C^1 curve on the sphere with
    ||v'|| <= ||(I-P)H'v||/g = 2 wbar^3 E(wbar)/g(wbar)                      (V1)
(by (I3) and H' = 2 wbar A_x), so
    sin angle(vhat, v(wbar)) <= s_p + hw * sup_{[a,b]} ||v'||  =: s          (V2)
(using sin(A+B) <= sin A + sin B).  (V1)-(V2) and
    E <= ehat_psi + s (N + |qhat_psi|)                                       (V3)
form a monotone system: starting from the crude E <= N (operator norm bound)
and iterating (V1)->(V3) gives a DECREASING sequence of valid bounds.  This is
the one place where the naive bound ||v'|| <= ||H'||/g would be ~20x too weak
(it would force ~20x more cells); (V1) is sharp because ||(I-P)H'v|| is exactly
the quantity the certificate is already bounding.

(2.3) THE TWO SIDES.  With  N >= sup_{[a,b]} ||O[f_wbar]||  (Schur row sum of
the majorant m_a) and R3 >= Schur row sum of 1/(k^2+a^2)^3:
    |d/dwbar ehat_f| <= ||O[d_wbar f]||,  |d_wbar psi| <= 4b/(k^2+a^2)^3,
    |d_wbar phi| <= 6b k^2/(k^2+a^2)^4 <= 6b/(k^2+a^2)^3,
so  ehat_f(wbar) <= ehat_f(p) + hw * (4b or 6b) * R3, and with v = c vhat + s u:
    F(wbar) = <v,O[phi]v> >= qhat_phi_lo - s^2 (|qhat_phi_hi| + N) - 2 s ehat_phi
    E(wbar) = ||(I-P)O[psi]v|| = min_t ||O[psi]v - t v||
                              <= ehat_psi + s (N + |qhat_psi_hi|).
    g(wbar) >= max(theta - L1 hw, mu_2) - alpha =: g_lo,  L1 = schur_full_seg_w0(a,b,1).
CERTIFICATE:   b^4 E_hi^2 / g_lo  <  F_lo   ==>   sup_{x in [a^2,b^2]} H''(x) < 0.

(2.4) THE FIRST CELL (0, wbar_0] (Stage S) uses h~ = G(2) - wbar^2 K(wbar)
instead of H: thetatilde >= lamhat_G2 - wbar_0^2 qhat_kappa_hi,
alphatilde = alpha_G2 + wbar_0^2 Lam, ||rtilde|| <= ||r_G2|| + wbar_0^2(Lam+|qhat_kappa|),
s = ||rtilde||/(thetatilde-alphatilde), g(wbar) = wbar^2 g~ >= wbar^2 (thetatilde-alphatilde),
so the positive term of (C**) is 2 wbar^2 E^2/g~ <= 2 wbar_0^2 E_hi^2/g~_lo:
it vanishes at wbar=0 and the cell reaches x = 0 with no gap.

NUMERICAL DISCIPLINE.  Every number that enters a bound is mpmath dps=30 and is
padded by ROUND_PAD = 1e-12 in the safe direction.  UNLIKE the base files the
kernel tables here are built in mpmath from exact inputs, so there is no
float64 entry rounding to absorb at all except inside the IMPORTED certificates
(`alpha_certificate_H`, whose head is float64 -- that is exactly what its own
ROUND_PAD covers).  numpy float64 `eigh` appears ONLY as a trial-vector
generator and in the evidence tables.  All loops have fixed bounds; the one
adaptive refinement has a hard cap.
"""

import os
import sys
import time

import numpy as np
from mpmath import mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_interchange_uniformity as piu                                   # noqa: E402
from p4_interchange_uniformity import (H_mat, H_matrix, G_matrix,         # noqa: E402
                                       ROUND_PAD, spec_point, spec_seg,
                                       envelope, schur_coupling,
                                       schur_tailblock, schur_full_seg,
                                       alpha_certificate_H, mu2_cell_lb)
from p4_interval_extension import schur_full_seg_w0                       # noqa: E402
from p4_conj2_interval import coupling_restricted                         # noqa: E402
import p4_conj2_small_w as psw                                            # noqa: E402

mp.dps = 30

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


ROUND = ROUND_PAD                      # 1e-12, the base-file pad convention
SQ2 = mp.sqrt(2)


# ===========================================================================
# KERNELS  (mpmath tables, index 0..kmax; index 0 is never used by O[.])
# ===========================================================================

def tab_H(w, kmax):
    """f(k) = w^2/(k^2+w^2)  -- the kernel of H(wbar)."""
    w2 = mp.mpf(w) ** 2
    return [mp.zero] + [w2 / (mp.mpf(k) ** 2 + w2) for k in range(1, kmax + 1)]


def tab_phi(w, kmax):
    """phi_x(k) = k^2/(k^2+w^2)^3,  x = w^2  --  A_xx = -2 O[phi_x]."""
    w2 = mp.mpf(w) ** 2
    out = [mp.zero]
    for k in range(1, kmax + 1):
        k2 = mp.mpf(k) ** 2
        out.append(k2 / (k2 + w2) ** 3)
    return out


def tab_psi(w, kmax):
    """psi_x(k) = 1/(k^2+w^2)^2  --  (I-P)A_x v = -x (I-P) O[psi_x] v."""
    w2 = mp.mpf(w) ** 2
    return [mp.zero] + [1 / (mp.mpf(k) ** 2 + w2) ** 2 for k in range(1, kmax + 1)]


def tab_kappa(w, kmax):
    """kappa_wbar(k) = 1/(k^2(k^2+w^2))  --  H(wbar)/wbar^2 = G(2) - wbar^2 K(wbar)."""
    w2 = mp.mpf(w) ** 2
    out = [mp.zero]
    for k in range(1, kmax + 1):
        k2 = mp.mpf(k) ** 2
        out.append(1 / (k2 * (k2 + w2)))
    return out


def tab_gamma(gamma, kmax):
    """k^{-gamma}  --  G(gamma) = O[k^{-gamma}]."""
    g = mp.mpf(gamma)
    return [mp.zero] + [1 / mp.mpf(k) ** g for k in range(1, kmax + 1)]


def np_tab(kmax, fun):
    k = np.arange(0, kmax + 1, dtype=np.float64)
    out = np.zeros_like(k)
    out[1:] = fun(k[1:])
    return out


def op_np(N, ftab):
    """float64 O[f] from a float64 kernel table (trial vectors / evidence ONLY)."""
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


def op_mp(N, ft):
    """exact mpmath O[f] head block of size N (Stage A gate only)."""
    sq = [mp.sqrt(2 * n + 1) for n in range(N)]
    diag = []
    run = mp.zero
    for k in range(N):
        if k >= 1:
            run += 2 * ft[k]
        diag.append(run / (2 * k + 1))
    M = mp.matrix(N)
    for n in range(N):
        for m in range(N):
            v = ft[n + m + 1] / (sq[n] * sq[m])
            if n != m:
                v += ft[abs(n - m)] / (sq[n] * sq[m])
            else:
                v -= diag[n]
            M[n, m] = v
    return M


# ===========================================================================
# TRIAL VECTOR, STRUCTURE COEFFICIENTS, MATVEC
# ===========================================================================

class Trial:
    """A finite-support unit trial vector vhat together with everything that can
    be precomputed from it: u_n = vhat_n/sqrt(2n+1), the structure coefficients
    c_k of Sec.2.1, and Cu = sum|u_n| (the matvec tail constant)."""

    def __init__(self, vnp):
        Nv = len(vnp)
        nrm = float(np.linalg.norm(vnp))
        v = [mp.mpf(float(x) / nrm) for x in vnp]
        nv2 = mp.fdot(v, v)
        v = [x / mp.sqrt(nv2) for x in v]              # exactly normalised in mpmath
        self.Nv = Nv
        self.v = v
        self.u = [v[n] / mp.sqrt(2 * n + 1) for n in range(Nv)]
        self.Cu = mp.fsum([abs(t) for t in self.u])
        u = self.u
        P = [mp.zero] * (2 * Nv + 2)
        acc = mp.zero
        for n in range(Nv - 1, 0, -1):
            acc += u[n] * u[n]
            P[n] = 2 * acc
        P[0] = 2 * (acc + u[0] * u[0])
        T = [mp.zero] * (2 * Nv + 2)
        for j in range(1, Nv):
            T[j] = 2 * mp.fdot(u[:Nv - j], u[j:])
        U = [mp.zero] * (2 * Nv + 2)
        for s in range(1, 2 * Nv):
            lo = max(0, s - Nv)
            hi = min(s - 1, Nv - 1)
            if hi < lo:
                continue
            U[s] = mp.fdot(u[lo:hi + 1], u[s - 1 - hi:s - lo][::-1])
        self.c = [T[k] + U[k] - P[k] for k in range(2 * Nv)]
        self.c[0] = mp.zero
        self.cp = [t if t > 0 else mp.zero for t in self.c]
        self.cm = [-t if t < 0 else mp.zero for t in self.c]
        self.all_nonneg = all(t >= 0 for t in self.c)

    def qform(self, ft):
        """<vhat, O[f] vhat> exactly, via the structure coefficients."""
        n = min(len(ft), len(self.c))
        return mp.fdot(self.c[:n], ft[:n])

    def qform_lo(self, ft_b, ft_a):
        n = min(len(ft_a), len(ft_b), len(self.c))
        return mp.fdot(self.cp[:n], ft_b[:n]) - mp.fdot(self.cm[:n], ft_a[:n])

    def qform_hi(self, ft_b, ft_a):
        n = min(len(ft_a), len(ft_b), len(self.c))
        return mp.fdot(self.cp[:n], ft_a[:n]) - mp.fdot(self.cm[:n], ft_b[:n])

    def matvec(self, ft, M):
        """y_n = (O[f] vhat)_n, n < M, exact mpmath.  Needs ft up to index M+Nv."""
        Nv, u, v = self.Nv, self.u, self.v
        zero = mp.zero
        y = []
        run = zero
        for n in range(M):
            if n >= 1:
                run += 2 * ft[n]
            if n < Nv:
                row2 = ft[n:0:-1] + [zero] + ft[1:Nv - n]
            else:
                row2 = ft[n:n - Nv:-1]
            s = mp.fdot(row2, u) + mp.fdot(ft[n + 1:n + Nv + 1], u)
            yn = s / mp.sqrt(2 * n + 1)
            if n < Nv:
                yn -= (run / (2 * n + 1)) * v[n]
            y.append(yn)
        return y

    def tail_sq(self, spec, M):
        """Upper bound on sum_{n>=M} y_n^2 for |f| <= env(spec), env non-increasing.

        For n >= Nv the L1 term drops (vhat_n = 0) and
          |y_n| <= (1/sqrt(2n+1)) sum_m [F(n-m)+F(n+m+1)] |u_m|
                <= (1/sqrt(2n+1)) 2 F(n-Nv+1) Cu
        (F non-increasing, n-m >= n-Nv+1 and n+m+1 >= n+1 >= n-Nv+1), hence
          sum_{n>=M} y_n^2 <= (4 Cu^2/(2M+1)) sum_{j>=M-Nv+1} F(j)^2.
        With F(k) = C/(k^2+A^2)^p <= C/k^{2p}:  sum_{j>=J} F(j)^2 <=
        C^2/((4p-1)(J-1)^{4p-1})  (integral test, J>=2)."""
        C, A, p = spec
        J = M - self.Nv + 1
        chk("matvec tail: row cut beyond the support", J >= 2, f"M={M} Nv={self.Nv}")
        s2 = C ** 2 / ((4 * p - 1) * mp.mpf(J - 1) ** (4 * p - 1))
        return 4 * self.Cu ** 2 / (2 * mp.mpf(M) + 1) * s2


# ===========================================================================
# ENVELOPES:  spec = (C, A, p)  means  env(k) = C/(k^2+A^2)^p,  A >= 0
# ===========================================================================

def env_at(spec, t):
    C, A, p = spec
    return C / (mp.mpf(t) ** 2 + A ** 2) ** p


def tail_int(spec, x):
    """int_x^inf C dt/(t^2+A^2)^p  (closed form; the A=0 branch is elementary)."""
    C, A, p = spec
    x = mp.mpf(x)
    if A == 0:
        return C / ((2 * p - 1) * x ** (2 * p - 1))
    at = mp.pi / 2 - mp.atan(x / A)
    if p == 1:
        return C * at / A
    if p == 2:
        return C * (at / (2 * A ** 3) - x / (2 * A ** 2 * (x * x + A * A)))
    if p == 3:
        return C * (3 * at / (8 * A ** 5) - x / (4 * A ** 2 * (x * x + A * A) ** 2)
                    - 3 * x / (8 * A ** 4 * (x * x + A * A)))
    raise ValueError(p)


_S1 = {}


def S1_of(spec, K=3000):
    """Upper bound on sum_{k>=1} env(k):  partial sum + integral tail."""
    key = (mp.nstr(spec[0], 25), mp.nstr(spec[1], 25), spec[2])
    if key in _S1:
        return _S1[key]
    s = mp.fsum([env_at(spec, k) for k in range(1, K + 1)]) + tail_int(spec, K)
    _S1[key] = s
    return s


def tail_sum(spec, J):
    """Upper bound on sum_{k>=J} env(k)  (env decreasing => <= int_{J-1}^inf)."""
    if J <= 1:
        return S1_of(spec)
    return tail_int(spec, J - 1)


_RSUM = {}


def rsum(spec, n0=60, M=900):
    """RIGOROUS  ||O[f]|| <= max_n sum_m |O[f]_{nm}|  (Schur/Holmgren) for every
    kernel with |f(k)| <= env(k).  Same algebra as
    `p4_conj2_small_w.schur_rowsum_bound`, spec-driven (and A=0-safe):
      row_n <= (2/(2n+1)) sum_{k<=n} env(k)
             + sum_{m!=n} env(|n-m|)/(sqrt(2n+1)sqrt(2m+1))
             + sum_{m>=0} env(n+m+1)/(sqrt(2n+1)sqrt(2m+1));
    rows n<n0 exactly for m<M plus the closed-form remainders; rows n>=n0 by
    the decreasing bulk bound 2S1/(2n0+1) + 2S1/sqrt(2n0+1)
    + tail_sum(n0+1)/sqrt(2n0+1)."""
    key = (mp.nstr(spec[0], 25), mp.nstr(spec[1], 25), spec[2], n0, M)
    if key in _RSUM:
        return _RSUM[key]
    ev = [mp.zero] + [env_at(spec, k) for k in range(1, M + n0 + 3)]
    cum = [mp.zero] * (n0 + 1)
    acc = mp.zero
    for n in range(1, n0 + 1):
        acc += ev[n]
        cum[n] = acc
    invm = [1 / mp.sqrt(2 * m + 1) for m in range(M)]
    S1 = S1_of(spec)
    best = mp.zero
    for n in range(n0):
        row = 2 * cum[n] / (2 * n + 1)
        r2 = ev[n:0:-1] + [mp.zero] + ev[1:M - n] if n < M else ev[n:n - M:-1]
        acc2 = mp.fdot(r2[:M], invm) + mp.fdot(ev[n + 1:n + M + 1], invm)
        acc2 += (tail_sum(spec, M - n) + tail_sum(spec, M + n + 1)) / mp.sqrt(2 * M + 1)
        row += acc2 / mp.sqrt(2 * n + 1)
        if row > best:
            best = row
    bulk = (2 * S1 / (2 * n0 + 1) + 2 * S1 / mp.sqrt(2 * n0 + 1)
            + tail_sum(spec, n0 + 1) / mp.sqrt(2 * n0 + 1))
    out = max(best, bulk) + ROUND
    _RSUM[key] = out
    return out


def spec_maj(a):
    """(1, a, 2), i.e. env(k) = 1/(k^2+a^2)^2.  For every wbar >= a and k >= 1
        phi_wbar(k) = k^2/(k^2+wbar^2)^3 <= 1/(k^2+wbar^2)^2 <= 1/(k^2+a^2)^2,
        psi_wbar(k) = 1/(k^2+wbar^2)^2                       <= 1/(k^2+a^2)^2,
    so `rsum(spec_maj(a))` bounds ||O[phi_wbar]|| and ||O[psi_wbar]|| on [a,b].
    (For a = 0 it also dominates kappa_wbar(k) = 1/(k^2(k^2+wbar^2)) <= 1/k^4,
    which is how Stage S uses it; for a > 0 kappa is not needed.)"""
    return (mp.one, mp.mpf(a), 2)


def spec_d3(a, coef):
    """(coef, a, 3): coef/(k^2+a^2)^3."""
    return (mp.mpf(coef), mp.mpf(a), 3)


# ===========================================================================
# STAGE A -- self-check gate: local builders == certified builders; (I1)-(I3)
# ===========================================================================

def run_A():
    print("=== STAGE A: gate -- builders, identities (I1)-(I3), envelopes ===")
    t0 = time.time()

    # (a) local float64 / mpmath O[f] builders == the VERBATIM base builders
    for N in (6, 25, 60):
        for w in (0.1, 0.7, 2.4, 9.0):
            base = H_matrix(N, w, 0)
            base2 = H_mat(N, w, 0)
            mine = op_np(N, np_tab(2 * N + 2, lambda k: w * w / (k * k + w * w)))
            chk("op_np == H_matrix (verbatim T5 builder)",
                float(np.max(np.abs(mine - base))) < 1e-14, f"N={N} w={w}")
            chk("op_np == H_mat (vectorised twin)",
                float(np.max(np.abs(mine - base2))) < 1e-14, f"N={N} w={w}")
    for N in (6, 25):
        for gam in (2, 4):
            base = G_matrix(N, gam)
            mine = op_np(N, np_tab(2 * N + 2, lambda k: k ** (-float(gam))))
            chk("op_np == G_matrix (verbatim T5 builder)",
                float(np.max(np.abs(mine - base))) < 1e-14, f"N={N} gamma={gam}")
    for N in (5, 18):
        for w in (0.4, 3.0):
            A = op_mp(N, tab_H(w, 2 * N + 2))
            B = H_mat(N, w, 0)
            d = max(abs(A[i, j] - mp.mpf(float(B[i, j])))
                    for i in range(N) for j in range(N))
            chk("op_mp (exact mpmath) == H_mat to float64 rounding", d < mp.mpf("1e-15"),
                f"N={N} w={w} d={mp.nstr(d,4)}")

    # (b) IDENTITY (I1):  A_xx = (H'' - H'/w)/(4w^2) = -2 O[phi_x]
    for N in (8, 30):
        for w in (0.2, 0.9, 2.5, 7.0, 20.0):
            lhs = (H_mat(N, w, 2) - H_mat(N, w, 1) / w) / (4 * w * w)
            rhs = -2.0 * op_np(N, np_tab(2 * N + 2,
                                         lambda k: k * k / (k * k + w * w) ** 3))
            sc = max(1.0, float(np.max(np.abs(rhs))))
            chk("(I1) A_xx == -2 O[phi_x] elementwise",
                float(np.max(np.abs(lhs - rhs))) / sc < 1e-12, f"N={N} w={w}")

    # (c) IDENTITY (I2):  A_x = H'/(2w) = H/w^2 - w^2 O[psi_x]
    for N in (8, 30):
        for w in (0.2, 0.9, 2.5, 7.0, 20.0):
            lhs = H_mat(N, w, 1) / (2 * w)
            rhs = (H_mat(N, w, 0) / (w * w)
                   - w * w * op_np(N, np_tab(2 * N + 2,
                                             lambda k: 1.0 / (k * k + w * w) ** 2)))
            sc = max(1.0, float(np.max(np.abs(rhs))))
            chk("(I2) A_x == H/w^2 - x O[psi_x] elementwise",
                float(np.max(np.abs(lhs - rhs))) / sc < 1e-12, f"N={N} w={w}")

    # (d) Segment I's identity (I):  H(w)/w^2 = G(2) - w^2 K(w)
    for N in (8, 30):
        for w in (0.3, 1.5, 6.0):
            lhs = H_mat(N, w, 0) / (w * w)
            rhs = (op_np(N, np_tab(2 * N + 2, lambda k: k ** -2.0))
                   - w * w * op_np(N, np_tab(2 * N + 2,
                                             lambda k: 1.0 / (k * k * (k * k + w * w)))))
            chk("(I) H/w^2 == G(2) - w^2 K(w) elementwise",
                float(np.max(np.abs(lhs - rhs))) < 1e-13, f"N={N} w={w}")
    # phi_0 = psi_0 = kappa_0 = k^-4  =>  A_xx(0) = -2 G(4)
    for N in (12,):
        d = float(np.max(np.abs(op_np(N, np_tab(2 * N + 2, lambda k: k ** -4.0))
                                - G_matrix(N, 4))))
        chk("phi_0 = 1/k^4  =>  A_xx(0) = -2 G(4)", d < 1e-14, f"d={d}")

    # (e) structure coefficients == direct double sum
    for N in (40, 120):
        for w in (0.5, 3.0):
            vnp = np.linalg.eigh(H_mat(N, w))[1][:, -1]
            vnp = vnp * np.sign(vnp[0])
            tr = Trial(vnp)
            for name, ft in (("H", tab_H(w, 2 * N + 4)),
                             ("phi", tab_phi(w, 2 * N + 4)),
                             ("psi", tab_psi(w, 2 * N + 4))):
                q1 = tr.qform(ft)
                y = tr.matvec(ft, N)
                q2 = mp.fdot(tr.v, y)
                chk("structure-coefficient qform == <vhat, O[f] vhat> from the matvec",
                    abs(q1 - q2) < mp.mpf("1e-25"), f"{name} N={N} w={w} "
                    f"{mp.nstr(q1,18)} vs {mp.nstr(q2,18)}")
            chk("structure coefficients: sign split is exhaustive (c = c+ - c-)",
                all(abs(tr.c[k] - (tr.cp[k] - tr.cm[k])) < mp.mpf("1e-30")
                    for k in range(len(tr.c))))
            # Rayleigh quotient of the trial vector == float64 eigenvalue
            lam = float(np.linalg.eigvalsh(H_mat(N, w))[-1])
            chk("trial Rayleigh quotient == float64 lam_1",
                abs(float(tr.qform(tab_H(w, 2 * N + 4))) - lam) < 1e-12,
                f"N={N} w={w}")
            # the cell-uniform sign-split bracket really brackets qhat(wbar)
            a, b = 0.8 * w, 1.25 * w
            for tabf in (tab_phi, tab_psi, tab_kappa):
                lo = tr.qform_lo(tabf(b, 2 * N + 4), tabf(a, 2 * N + 4))
                hi = tr.qform_hi(tabf(b, 2 * N + 4), tabf(a, 2 * N + 4))
                for t in (0.0, 0.13, 0.5, 0.87, 1.0):
                    wt = a + t * (b - a)
                    q = tr.qform(tabf(wt, 2 * N + 4))
                    chk("cell-uniform sign-split brackets the exact quadratic form",
                        lo <= q <= hi, f"N={N} w={w} t={t} {mp.nstr(lo,12)} "
                        f"<= {mp.nstr(q,12)} <= {mp.nstr(hi,12)}")
            # _ehat == the float64 residual norm ||(I-P)O[f]vhat||
            ftb = np_tab(2 * N + 4, lambda k: 1.0 / (k * k + w * w) ** 2)
            Of = op_np(N, ftb)
            vv = np.array([float(t) for t in tr.v])
            yy = Of @ vv
            e64 = float(np.linalg.norm(yy - float(vv @ yy) * vv))
            chk("_ehat (mpmath, incl. row tail) >= the float64 head residual",
                float(_ehat(tr, tab_psi(w, 4 * N + 400), N + 300,
                            (mp.one, mp.mpf(0.8 * w), 2))) >= e64 - 1e-13,
                f"N={N} w={w}")

    # (f) my envelope closed forms == the certified `piu.envelope` ones (A>0)
    for (C, A, p) in ((mp.one, mp.mpf(2), 2), (mp.mpf(3), mp.mpf(0.7), 3),
                      (mp.mpf(9), mp.mpf(11), 1)):
        S1p, int_tail = envelope((C, A, p))
        for x in (1, 7, 200):
            chk("tail_int == piu.envelope int_tail",
                abs(tail_int((C, A, p), x) - int_tail(x)) < mp.mpf("1e-25"),
                f"C={C} A={A} p={p} x={x}")
        # both are "partial sum + integral tail" upper bounds on the same series;
        # ours cuts at K=3000 vs piu's K=2000, so ours is the tighter one.
        chk("S1_of is a (tighter) upper bound consistent with piu.envelope S1",
            S1_of((C, A, p)) <= S1p + mp.mpf("1e-15")
            and abs(S1_of((C, A, p)) - S1p) < mp.mpf("1e-5"),
            f"{mp.nstr(S1_of((C,A,p)),20)} vs {mp.nstr(S1p,20)}")
    # A = 0 branch
    chk("tail_int(A=0,p=2) == 1/(3 x^3)",
        abs(tail_int((mp.one, mp.zero, 2), 5) - mp.mpf(1) / 375) < mp.mpf("1e-28"))
    chk("S1_of((1,0,2)) == zeta(4)", abs(S1_of((mp.one, mp.zero, 2)) - mp.zeta(4))
        < mp.mpf("1e-12"), mp.nstr(S1_of((mp.one, mp.zero, 2)), 20))

    # (g) rsum really dominates the operator norm (float64 evidence) and is
    #     monotone in the exactness parameters
    for a in (0.0, 1.0, 5.0, 20.0):
        R = rsum(spec_maj(a))
        R2 = rsum(spec_maj(a), n0=120, M=2500)
        chk("rsum is stable under refining (n0,M)", abs(R - R2) <= abs(R) * mp.mpf("1e-6"),
            f"a={a} {mp.nstr(R,10)} vs {mp.nstr(R2,10)}")
        w = max(a, 0.05)
        Nb = 700
        for tb, nm in ((np_tab(2 * Nb + 2, lambda k: 1.0 / (k * k + w * w) ** 2), "psi"),
                       (np_tab(2 * Nb + 2, lambda k: k * k / (k * k + w * w) ** 3), "phi")):
            nn = float(np.max(np.abs(np.linalg.eigvalsh(op_np(Nb, tb)))))
            chk("evidence: rsum(majorant) >= ||O[f]|| at N=700", float(R) >= nn,
                f"a={a} {nm} rsum={float(R):.4e} norm={nn:.4e}")
    # cross-check against Segment I's own row-sum bound for the shared kernel 1/k^4
    Rs = psw.schur_rowsum_bound(psw.ker_gamma(4), n0=60, M=900)
    chk("rsum((1,0,2)) agrees with p4_conj2_small_w.schur_rowsum_bound[k^-4]",
        abs(rsum(spec_maj(0.0)) - Rs) < mp.mpf("1e-9"),
        f"{mp.nstr(rsum(spec_maj(0.0)),12)} vs {mp.nstr(Rs,12)}")

    # (h) float64 EVIDENCE that (P2)/(C**) is the right formula:
    #     H''(x) from central differences == -2F + 2x^2 sum_j .../(lam1-lam_j)
    Nb = 900
    for w in (0.3, 1.0, 3.0, 10.0):
        x = w * w
        A = H_mat(Nb, w)
        e, V = np.linalg.eigh(A)
        v = V[:, -1]
        Phi = op_np(Nb, np_tab(2 * Nb + 2, lambda k: k * k / (k * k + x) ** 3))
        Psi = op_np(Nb, np_tab(2 * Nb + 2, lambda k: 1.0 / (k * k + x) ** 2))
        F = float(v @ Phi @ v)
        cj = V.T @ (Psi @ v)
        pos = 2.0 * x * x * float(np.sum(cj[:-1] ** 2 / (e[-1] - e[:-1])))
        hpp = -2 * F + pos
        dx = max(1e-5, 1e-4 * x)
        f = lambda xx: float(np.linalg.eigvalsh(H_mat(Nb, np.sqrt(xx)))[-1])
        fd = (f(x + dx) - 2 * f(x) + f(x - dx)) / dx ** 2
        chk("evidence: (P2)+(I1)+(I3) reproduce the finite-difference H''(x)",
            abs(hpp - fd) < 1e-6 * max(1.0, abs(fd)) + 1e-9,
            f"w={w} formula={hpp:.10e} fd={fd:.10e}")
        chk("evidence: H''(x) < 0 at N=900", fd < 0, f"w={w}")
    _summary("Stage A", f"({time.time()-t0:.1f}s)")


# ===========================================================================
# STAGE B -- C^2 in operator norm (Sec.1.5)
# ===========================================================================

def run_B():
    print("=== STAGE B: x |-> A(x) is C^2 in operator norm (Sec.1.5) ===")
    t0 = time.time()
    print("  %16s  %12s  %13s  %14s" % ("cell", "L1=sup|H'|", "L2=sup|H''|",
                                        "L3=sup|H'''|"))
    for (a, b) in ((0.2, 0.3), (1.0, 1.1), (5.0, 5.2), (20.0, 20.5), (40.0, 41.0)):
        L1 = schur_full_seg_w0(a, b, 1)
        L2 = schur_full_seg_w0(a, b, 2)
        L3 = schur_full_seg_w0(a, b, 3)
        for j, L in ((1, L1), (2, L2), (3, L3)):
            chk("weighted-row-0 segment envelope <= base schur_full_seg",
                L <= schur_full_seg(a, b, j) + 1e-15, f"[{a},{b}] order {j}")
            chk("segment envelope is finite and positive", 0 < L < 1e6,
                f"[{a},{b}] order {j}")
        print(f"  [{a:6.2f},{b:6.2f}]  {L1:12.6f}  {L2:13.6f}  {L3:14.6f}")
        # float64 EVIDENCE: the order-3 Taylor remainder really obeys the bound
        N = 400
        d = 0.02 * (b - a)
        w = 0.5 * (a + b)
        R = (H_mat(N, w + d) - H_mat(N, w) - d * H_mat(N, w, 1)
             - d * d / 2 * H_mat(N, w, 2))
        nr = float(np.max(np.abs(np.linalg.eigvalsh(R))))
        chk("evidence: ||3rd-order Taylor remainder|| <= |d|^3 L3/6 at N=400",
            nr <= float(L3) * d ** 3 / 6 + 1e-15, f"[{a},{b}] {nr:.3e}")
    # near x = 0: A(x) = x G(2) - x^2 K~(x), kernel analytic in x, uniform majorant
    R0 = rsum(spec_maj(0.0))
    print(f"  near x=0: A(x) = x G(2) - x^2 K~(x) with |d^j_x kernel| <= j!/k^(2j+2),")
    print(f"            Schur majorant for every derivative order <= {float(R0):.6f} "
          f"=> A is ANALYTIC in x on |x| < 1 (operator norm).")
    chk("uniform Schur majorant at x=0 is finite", R0 < 10)
    _summary("Stage B", f"({time.time()-t0:.1f}s)")


# ===========================================================================
# STAGE S -- the first cell (0, wbar_0]:  rescaled operator h~ = G(2) - w^2 K(w)
# ===========================================================================

NV_S = 200            # trial-vector support for G(2)
M_RES_S = 2400        # row cut for the G(2) residual (kernel decays only like k^-2)
NA_G2 = 200           # head cut for the lam_2(G(2)) Cholesky-deflation certificate

SDATA = {}


def _stage_S_setup():
    if "tr" in SDATA:
        return SDATA
    ftab_G2 = tab_gamma(2, 2 * NV_S + M_RES_S + 8)
    vnp = np.linalg.eigh(op_np(NV_S, np_tab(2 * NV_S + 2, lambda k: k ** -2.0)))[1][:, -1]
    vnp = vnp * np.sign(vnp[0])
    chk("Stage S trial vector is entrywise positive (Perron structure)",
        float(np.min(vnp)) > 0, f"min={float(np.min(vnp)):.3e}")
    tr = Trial(vnp)
    lamhat = tr.qform(ftab_G2)                       # <= frak_g(2)  (Rayleigh)
    # residual of G(2) at vhat, with the row tail bounded by the k^-2 majorant
    y = tr.matvec(ftab_G2, M_RES_S)
    r2 = mp.zero
    for n in range(M_RES_S):
        dd = y[n] - (lamhat * tr.v[n] if n < NV_S else mp.zero)
        r2 += dd * dd
    r2 += tr.tail_sq((mp.one, mp.zero, 1), M_RES_S)
    rG2 = mp.sqrt(r2) + ROUND
    alpha_G2, cG2, nDG2, nBG2 = psw.alpha_G2(NA_G2)
    chk("alpha_G2 < lamhat (Davis-Kahan hypothesis on G(2))", alpha_G2 < lamhat)
    chk("alpha_G2 beats the P4-C published 2.0e-2 bound on lam_2(G(2))",
        alpha_G2 <= mp.mpf("2.1e-2"), mp.nstr(alpha_G2, 8))
    SDATA.update(tr=tr, lamhat=lamhat, rG2=rG2, alpha_G2=alpha_G2, fG2=ftab_G2)
    print(f"  G(2) trial: lamhat = {mp.nstr(lamhat,16)} <= frak_g(2),  "
          f"||r_G2|| = {mp.nstr(rG2,6)},  alpha_G2 = {mp.nstr(alpha_G2,8)}")
    return SDATA


def cell_S(a, b, verbose=False):
    """Certified sup_{x in [a^2,b^2]} H''(x) for 0 <= a <= b <= (small), via the
    RESCALED operator h~(w) = H(w)/w^2 = G(2) - w^2 K(w)  (Sec.2.4)."""
    D = _stage_S_setup()
    tr, lamhat, rG2, alG2 = D["tr"], D["lamhat"], D["rG2"], D["alpha_G2"]
    a = mp.mpf(a)
    b = mp.mpf(b)
    KT = 2 * tr.Nv + 600
    MK = tr.Nv + 300
    fk_a, fk_b = tab_kappa(a, KT), tab_kappa(b, KT)
    fphi_a, fphi_b = tab_phi(a, KT), tab_phi(b, KT)
    fpsi_a, fpsi_b = tab_psi(a, KT), tab_psi(b, KT)

    # --- kappa quadratic form and the norm of K(w) on the cell ---------------
    qk_lo = tr.qform_lo(fk_b, fk_a) - ROUND
    qk_hi = tr.qform_hi(fk_b, fk_a) + ROUND
    # kappa_w(k) = 1/(k^2(k^2+w^2)) <= 1/k^4 for every w >= 0
    NK = rsum((mp.one, mp.zero, 2))
    theta_lo = lamhat - b * b * max(qk_hi, mp.zero) - ROUND
    alpha_t = alG2 + b * b * NK + ROUND
    g_lo = theta_lo - alpha_t
    if g_lo <= 0:
        return None, dict(reason="no certified gap in the rescaled operator",
                          g=float(g_lo))
    rt = rG2 + b * b * (NK + max(abs(qk_lo), abs(qk_hi))) + ROUND
    s = rt / g_lo
    if s >= 1:
        return None, dict(reason="Davis-Kahan angle >= 1", s=float(s))

    # --- phi / psi on the cell ------------------------------------------------
    Nmaj = rsum((mp.one, mp.zero, 2))
    R3 = rsum((mp.one, mp.zero, 3))
    qphi_lo = tr.qform_lo(fphi_b, fphi_a) - ROUND
    qphi_hi = tr.qform_hi(fphi_b, fphi_a) + ROUND
    qpsi_lo = tr.qform_lo(fpsi_b, fpsi_a) - ROUND
    qpsi_hi = tr.qform_hi(fpsi_b, fpsi_a) + ROUND
    # |qhat_f(wbar)| <= max(|lo|,|hi|) for every wbar in the cell (the sign of the
    # quadratic form is never assumed)
    aphi = max(abs(qphi_lo), abs(qphi_hi))
    apsi = max(abs(qpsi_lo), abs(qpsi_hi))
    p = (a + b) / 2
    fphi_p, fpsi_p = tab_phi(p, KT), tab_psi(p, KT)
    eph = _ehat(tr, fphi_p, MK, (mp.one, mp.zero, 2))
    eps_ = _ehat(tr, fpsi_p, MK, (mp.one, mp.zero, 2))
    hw = (b - a) / 2
    eph_hi = eph + hw * 6 * b * R3 + ROUND
    eps_hi = eps_ + hw * 4 * b * R3 + ROUND

    F_lo = qphi_lo - s * s * (aphi + Nmaj) - 2 * s * eph_hi
    E_hi = eps_hi + s * (Nmaj + apsi)
    pos = b * b * E_hi * E_hi / g_lo            # = sup 2x^2E^2/g / 2, since g=w^2 g~
    hpp_ub = -2 * F_lo + 2 * pos
    info = dict(a=float(a), b=float(b), F_lo=float(F_lo), pos=float(pos),
                s=float(s), g_lo=float(g_lo), E_hi=float(E_hi), theta=float(theta_lo),
                alpha=float(alpha_t), qk=float(qk_hi), ehat_phi=float(eph),
                ehat_psi=float(eps_), N=float(Nmaj), hpp_ub=float(hpp_ub),
                ratio=float(pos / F_lo) if F_lo > 0 else 9e9, regime="S")
    return float(hpp_ub), info


def _ehat(tr, ft, M, spec):
    """||(I - vhat vhat*) O[f] vhat||, exact mpmath head + closed-form row tail."""
    y = tr.matvec(ft, M)
    q = mp.fdot(tr.v, y[:tr.Nv])
    s2 = mp.zero
    for n in range(M):
        dd = y[n] - (q * tr.v[n] if n < tr.Nv else mp.zero)
        s2 += dd * dd
    s2 += tr.tail_sq(spec, M)
    return mp.sqrt(s2) + ROUND


# ===========================================================================
# STAGE L -- the ladder on [wbar_0, X]:  operator H(wbar) itself
# ===========================================================================

def support_size(w):
    """Trial-vector support.  The residual ||H(p)vhat - theta vhat|| scales like
    w^2/Ns^2, so Ns ~ 10 w keeps the Davis-Kahan angle ~5e-3 (same rule as
    Stage M of `p4_conj2_interval`)."""
    return max(140, int(10 * float(w)) // 2 * 2)


def na_band(w):
    """Head size for the lam_2 Cholesky-deflation certificate (cost is O(Na^3),
    the envelope is O(w/Na); this banding keeps alpha well below theta while
    the coarse grid stays affordable)."""
    if w <= 5.0:
        return 70
    if w <= 12.0:
        return 140
    return 280


_TRIAL = {}


def trial_at(p, Ns):
    key = (float(p), Ns)
    if key in _TRIAL:
        return _TRIAL[key]
    vnp = np.linalg.eigh(H_mat(Ns, float(p)))[1][:, -1]
    vnp = vnp * np.sign(vnp[0])
    tr = Trial(vnp)
    _TRIAL.clear()
    _TRIAL[key] = tr
    return tr


def cell_L(a, b, alpha, L1, mu2):
    """Certified sup_{x in [a^2,b^2]} H''(x) on a cell [a,b] with 0 < a <= b."""
    a, b = mp.mpf(a), mp.mpf(b)
    p = (a + b) / 2
    hw = (b - a) / 2
    Ns = support_size(p)
    Nr = 4 * Ns
    tr = trial_at(p, Ns)
    KT = Nr + Ns + 8
    fH = tab_H(p, KT)

    # --- Rayleigh, residual (with the MANDATORY tail coupling), angle at p ----
    theta = tr.qform(fH) - ROUND
    y = tr.matvec(fH, Nr)
    r2 = mp.zero
    for n in range(Nr):
        dd = y[n] - (theta * tr.v[n] if n < Ns else mp.zero)
        r2 += dd * dd
    rtail = coupling_restricted(Ns, Nr, spec_point(float(p), 0))
    rnorm = mp.sqrt(r2 + rtail * rtail) + ROUND
    gamma = theta - alpha
    g_lo = max(theta - mp.mpf(L1) * hw, mu2) - alpha
    if gamma <= 0 or g_lo <= 0:
        return None, dict(reason="no certified gap", gamma=float(gamma), g=float(g_lo))
    s_p = rnorm / gamma

    # --- phi / psi pieces ----------------------------------------------------
    MK = Ns + 300
    sp2 = spec_maj(a)
    R3 = rsum(spec_d3(a, 1))
    Nmaj = rsum(sp2)
    fphi_a, fphi_b, fphi_p = tab_phi(a, KT), tab_phi(b, KT), tab_phi(p, KT)
    fpsi_a, fpsi_b, fpsi_p = tab_psi(a, KT), tab_psi(b, KT), tab_psi(p, KT)
    qphi_lo = tr.qform_lo(fphi_b, fphi_a) - ROUND
    qphi_hi = tr.qform_hi(fphi_b, fphi_a) + ROUND
    qpsi_lo = tr.qform_lo(fpsi_b, fpsi_a) - ROUND
    qpsi_hi = tr.qform_hi(fpsi_b, fpsi_a) + ROUND
    aphi = max(abs(qphi_lo), abs(qphi_hi))
    apsi = max(abs(qpsi_lo), abs(qpsi_hi))
    eph = _ehat(tr, fphi_p, MK, sp2)
    eps_ = _ehat(tr, fpsi_p, MK, sp2)
    eph_hi = eph + hw * 6 * b * R3 + ROUND
    eps_hi = eps_ + hw * 4 * b * R3 + ROUND

    # --- cell-uniform angle by the monotone bootstrap (V1)-(V3) --------------
    s = min(mp.one, s_p + hw * 2 * b ** 3 * Nmaj / g_lo)
    for _ in range(12):
        E_hi = eps_hi + s * (Nmaj + apsi)
        s_new = s_p + hw * 2 * b ** 3 * E_hi / g_lo
        if s_new >= s:
            break
        s = s_new
    if s >= 1:
        return None, dict(reason="Davis-Kahan angle >= 1", s=float(s))
    E_hi = eps_hi + s * (Nmaj + apsi)
    F_lo = qphi_lo - s * s * (aphi + Nmaj) - 2 * s * eph_hi
    pos = b ** 4 * E_hi * E_hi / g_lo
    hpp_ub = -2 * F_lo + 2 * pos
    info = dict(a=float(a), b=float(b), p=float(p), Ns=Ns, Nr=Nr,
                F_lo=float(F_lo), pos=float(pos), s=float(s), s_p=float(s_p),
                g_lo=float(g_lo), E_hi=float(E_hi), theta=float(theta),
                alpha=float(alpha), rnorm=float(rnorm), rtail=float(rtail),
                ehat_phi=float(eph), ehat_psi=float(eps_), N=float(Nmaj),
                hpp_ub=float(hpp_ub),
                ratio=float(pos / F_lo) if F_lo > 0 else 9e9, regime="L",
                cnn=bool(tr.all_nonneg),
                lip_psi=float(hw * 4 * b * R3), ang_psi=float(s * (Nmaj + apsi)),
                gamma=float(gamma))
    return float(hpp_ub), info


def alpha_cell(A, B, cA, cB, L1, Na):
    """Cell-uniform  alpha >= sup_{w in [A,B]} lam_2(H(w))  --  the Stage M /
    Segment II algebra with an explicit head size:
       lam_2(H(w)) <= lam_2(head) + ||B||_seg + ||D||_seg
                   <= max(c_A,c_B) + L1 (B-A)/2 + ||B||_seg + ||D||_seg."""
    sp = spec_seg(A, B, 0)
    nB = schur_coupling(Na, B, 0, spec=sp)
    nD = schur_tailblock(Na, B, 0, spec=sp)
    return max(cA, cB) + mp.mpf(L1) * (mp.mpf(B) - mp.mpf(A)) / 2 + nB + nD + ROUND


COARSE_SLACK = 0.35        # target for the L1*(coarse width)/2 term of alpha
COARSE_REL_MAX = 0.20      # and at most this relative coarse width
MAX_SPLIT = 6              # hard cap on adaptive doublings per coarse cell
TARGET = 0.75              # accept a fine cell when pos/F_lo <= TARGET


def coarse_grid(W0, WMAX):
    g = [float(W0)]
    while g[-1] < WMAX:
        w = g[-1]
        L1w = float(schur_full_seg_w0(w, w * 1.02, 1))
        d = min(2 * COARSE_SLACK / max(L1w, 1e-9), COARSE_REL_MAX * w)
        g.append(min(WMAX, w + d))
    return g


LAD = {"rows": [], "W": None}


def run_L(W0, WMAX, tag="Stage L"):
    print(f"=== {tag}: concavity ladder on [{W0}, {WMAX}] ===")
    t0 = time.time()
    coarse = coarse_grid(W0, WMAX)
    print(f"  coarse alpha grid: {len(coarse)} nodes (Cholesky-deflation, banded Na); "
          f"fine cells adaptive (cap {MAX_SPLIT} doublings, target pos/F <= {TARGET})")
    ca = []
    for w in coarse:
        Na = na_band(w)
        al, ok, c, _, _ = alpha_certificate_H(float(w), Na=Na)
        chk("Cholesky-deflation lam_2 certificate SPD (coarse node)", ok, f"w={w:.4f}")
        ca.append((c, Na))

    rows = []
    worst = None
    stopped = None
    for i in range(len(coarse) - 1):
        A, B = coarse[i], coarse[i + 1]
        Na = min(ca[i][1], ca[i + 1][1])
        cA, cB = ca[i][0], ca[i + 1][0]
        if ca[i][1] != ca[i + 1][1]:
            _, okA, cA, _, _ = alpha_certificate_H(A, Na=Na)
            _, okB, cB, _, _ = alpha_certificate_H(B, Na=Na)
            chk("band-boundary lam_2 recertification SPD", okA and okB, f"[{A},{B}]")
        L1 = schur_full_seg_w0(A, B, 1)
        chk("L1 (weighted row 0) <= base schur_full_seg order 1",
            L1 <= schur_full_seg(A, B, 1) + 1e-15, f"[{A:.4f},{B:.4f}]")
        alpha = alpha_cell(A, B, cA, cB, L1, Na)
        mu2, k1, k2, k3 = mu2_cell_lb(A, B)
        chk("mu_2 segment monotonicity k1,k2,k3 >= 0", k1 >= 0 and k2 >= 0 and k3 >= 0,
            f"[{A:.4f},{B:.4f}]")

        k = 1
        got = None
        for _ in range(MAX_SPLIT + 1):
            sub = np.linspace(A, B, k + 1)
            acc = []
            bad = False
            for j in range(k):
                ub, info = cell_L(float(sub[j]), float(sub[j + 1]), alpha, L1, mu2)
                if ub is None or ub >= 0 or info["ratio"] > TARGET:
                    bad = True
                    break
                acc.append(info)
            if not bad:
                got = acc
                break
            k *= 2
        if got is None:
            stopped = A
            print(f"  !! coarse cell [{A:.4f},{B:.4f}] did NOT close after {MAX_SPLIT} "
                  f"doublings (k={k//2}); stopping the ladder here.")
            break
        for info in got:
            chk("fine cell: certified H''(x) < 0 for all x in [a^2,b^2]",
                info["hpp_ub"] < 0, f"[{info['a']:.5f},{info['b']:.5f}]")
            rows.append(info)
            if worst is None or info["ratio"] > worst["ratio"]:
                worst = info

    Wdone = rows[-1]["b"] if rows else float(W0)
    LAD["rows"] = rows
    LAD["W"] = Wdone
    LAD["worst"] = worst
    LAD["stopped"] = stopped

    print(f"  {'band':>16} {'cells':>6} {'Ns':>5} {'max pos/F':>10} {'F_lo':>11} "
          f"{'pos':>11} {'g_lo':>6} | E_hi split %: {'ehat':>5} {'width':>5} {'angle':>5}"
          f" | {'s':>9} {'s_p':>9}")
    edges = [float(W0), 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 45.0, 1e9]
    edges = [e for e in edges if e >= float(W0)]
    edges = sorted(set([float(W0)] + [e for e in edges if e <= Wdone + 1e-9] + [Wdone]))
    for e0, e1 in zip(edges[:-1], edges[1:]):
        sel = [r for r in rows if r["a"] >= e0 - 1e-9 and r["b"] <= e1 + 1e-9]
        if not sel:
            continue
        wc = max(sel, key=lambda r: r["ratio"])
        tot = wc["E_hi"]
        print(f"  [{e0:6.3f},{e1:7.3f}] {len(sel):>6d} {wc['Ns']:>5d} {wc['ratio']:>10.4f} "
              f"{wc['F_lo']:>11.4e} {wc['pos']:>11.4e} {wc['g_lo']:>6.3f} | "
              f"               {100*wc['ehat_psi']/tot:5.1f} {100*wc['lip_psi']/tot:5.1f} "
              f"{100*wc['ang_psi']/tot:5.1f} | {wc['s']:>9.3e} {wc['s_p']:>9.3e}")
    if worst is not None:
        bind = max((("Davis-Kahan angle (s_p, i.e. the trial-vector residual)",
                     worst["ang_psi"]),
                    ("cell-width Lipschitz of ehat_psi", worst["lip_psi"]),
                    ("ehat_psi itself (structural)", worst["ehat_psi"])),
                   key=lambda t: t[1])
        print(f"  binding fine cell [{worst['a']:.5f},{worst['b']:.5f}]: pos/F_lo = "
              f"{worst['ratio']:.4f}, F_lo = {worst['F_lo']:.5e}, pos = {worst['pos']:.5e}")
        print(f"    E_hi = {worst['E_hi']:.4e} = ehat_psi {worst['ehat_psi']:.4e} "
              f"+ width {worst['lip_psi']:.4e} + angle {worst['ang_psi']:.4e}; "
              f"s = {worst['s']:.3e} (s_p = {worst['s_p']:.3e}), g_lo = {worst['g_lo']:.4f}")
        print(f"    ==> the LARGEST contribution to E_hi is: {bind[0]}")
    print(f"  ==> H''(x) < 0 CERTIFIED on x in (0, {Wdone**2:.4f}] "
          f"(wbar up to {Wdone:.4f}), {len(rows)} fine cells")
    nn = sum(1 for r in rows if r.get("cnn"))
    print(f"      (structure coefficients c_k >= 0 on {nn}/{len(rows)} cells -- the "
          f"sign split of Sec.2.1 is used regardless, so this is a remark, not a "
          f"hypothesis)")
    _summary(tag, f"({time.time()-t0:.1f}s)")
    return Wdone


# ===========================================================================
# STAGE S driver -- pick the largest first cell (0, wbar_0] that closes
# ===========================================================================

SFIRST = {}


def run_S(candidates=(0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15)):
    print("=== STAGE S: the first cell (0, wbar_0] -- reaches x = 0 exactly ===")
    t0 = time.time()
    _stage_S_setup()
    chosen = None
    for w0 in candidates:
        ub, info = cell_S(0.0, w0)
        if ub is not None and ub < 0 and info["ratio"] <= 0.85:
            chosen = (w0, ub, info)
            break
    chk("first cell (0, wbar_0] closes for some wbar_0", chosen is not None)
    w0, ub, info = chosen
    SFIRST.update(w0=w0, ub=ub, info=info)
    print(f"  wbar_0 = {w0}   ==>  sup_{{x in (0,{w0**2:.4f}]}} H''(x) <= {ub:.6f} < 0")
    print(f"    F_lo = {info['F_lo']:.6f}, positive term = {info['pos']:.6f} "
          f"(pos/F = {info['ratio']:.4f}), sin angle <= {info['s']:.4e}, "
          f"g~_lo = {info['g_lo']:.6f}")
    chk("first cell: H''(x) < 0 on (0, wbar_0^2]", ub < 0)
    # the x -> 0 limit itself:  H''(0) = -2 c_2
    ub0, info0 = cell_S(0.0, 0.0)
    chk("degenerate cell at x=0 closes", ub0 is not None and ub0 < 0)
    print(f"  x -> 0 limit: H''(0) <= {ub0:.9f}  (positive term vanishes identically; "
          f"c_2 >= {info0['F_lo']:.9f})")
    chk("H''(0) <= -2*1.193 (matches c_2 from Segment I)", ub0 <= -2 * 1.193,
        f"ub0={ub0}")
    chk("c_2 lower bound agrees with Segment I's [1.19312762268, 1.19320906584]",
        1.1931 <= info0["F_lo"] <= 1.19321, f"{info0['F_lo']}")
    print(f"  ==> RHS of (C*) = 4 wbar^2 H''  ~  -8 c_2 wbar^2 = "
          f"{-8*info0['F_lo']:.6f} wbar^2   as wbar -> 0  (no wbar_min gap)")
    _summary("Stage S", f"({time.time()-t0:.1f}s)")
    return w0


# ===========================================================================
# STAGE E -- evidence table: certified upper bound vs float64 truth
# ===========================================================================

def _point_bound(w):
    """Certified upper bound on 4 wbar^2 H''(wbar^2) at the single point wbar."""
    if w <= SFIRST.get("w0", 0.3):
        ub, info = cell_S(w, w)
        return ub, info
    Na = na_band(w)
    al, ok, c, _, _ = alpha_certificate_H(float(w), Na=Na)
    chk("point-evaluation Cholesky-deflation SPD", ok, f"w={w}")
    d = 1e-6 * w
    A, B = w - d, w + d
    L1 = schur_full_seg_w0(A, B, 1)
    alpha = alpha_cell(A, B, c, c, L1, Na)
    mu2, _, _, _ = mu2_cell_lb(A, B)
    return cell_L(A, B, alpha, L1, mu2)


def run_E(points):
    print("=== STAGE E: certified upper bound on 4 wbar^2 H''  vs  float64 truth ===")
    t0 = time.time()
    print(f"  {'wbar':>8} {'x':>10} {'4w^2 H'' UB (cert)':>20} "
          f"{'4w^2 H'' (fd, N=1500)':>22} {'UB/truth':>9} {'regime':>7}")
    N = 1500
    for w in points:
        ub, info = _point_bound(w)
        if ub is None:
            print(f"  {w:8.3f}   -- no certificate ({info.get('reason')})")
            continue
        x = w * w
        dx = max(1e-6, 1e-4 * x)
        f = lambda xx: float(np.linalg.eigvalsh(H_mat(N, np.sqrt(xx)))[-1])
        fd = (f(x + dx) - 2 * f(x) + f(x - dx)) / dx ** 2
        cert = 4 * x * ub
        truth = 4 * x * fd
        chk("evidence: float64 4w^2 H'' < 0 at N=1500", truth < 0, f"w={w}")
        chk("evidence: certified upper bound >= the float64 truth",
            cert >= truth - 1e-9 * abs(truth) - 1e-12,
            f"w={w} cert={cert:.6e} truth={truth:.6e}")
        chk("certified upper bound is negative", cert < 0, f"w={w}")
        print(f"  {w:8.3f} {x:10.3f} {cert:20.8e} {truth:22.8e} "
              f"{cert/truth:9.4f} {info['regime']:>7}")
    _summary("Stage E", f"({time.time()-t0:.1f}s)")


# ===========================================================================
def main():
    args = [a.upper() for a in sys.argv[1:]]
    ext = "EXT" in args or "ALL" in args
    wmax = 45.0
    for a in args:
        try:
            wmax = float(a)
        except ValueError:
            pass
    run_A()
    run_B()
    w0 = run_S()
    W = run_L(w0, wmax if ext else 3.0, "Stage L" if ext else "Stage L (fast, wbar<=3)")
    pts = [0.1, 0.5, 1.0, 2.0]
    if W >= 5.0:
        pts += [5.0]
    if W >= 10.0:
        pts += [10.0]
    pts += [round(W, 4)]
    run_E([p for p in pts if p <= W + 1e-9])

    print("\n" + "=" * 78)
    print("SYNTHESIS -- P4 Einstein extremal, Segment C-I")
    print(f"  CERTIFIED:  H''(x) < 0  for all x in (0, {W**2:.4f}]   "
          f"(wbar in (0, {W:.4f}])")
    print(f"              ==> H(x) = h(sqrt x) is STRICTLY CONCAVE on (0, {W**2:.4f}],")
    print(f"                  and concave on [0, {W**2:.4f}] (H''(0) = -2 c_2 < 0).")
    print(f"              first cell (0, {SFIRST['w0']}] reaches x = 0 exactly: "
          f"NO wbar_min gap.")
    if LAD["worst"] is not None:
        wst = LAD["worst"]
        print(f"  binding cell [{wst['a']:.5f},{wst['b']:.5f}]: pos/F_lo = "
              f"{wst['ratio']:.4f}")
    if LAD["stopped"] is not None:
        print(f"  ladder stopped at wbar = {LAD['stopped']:.4f} (adaptive cap reached)")
    print(f"  ALL {PASS[0]} checks passed  ({time.time()-T0:.1f}s).")
    print("=" * 78)


if __name__ == "__main__":
    main()
