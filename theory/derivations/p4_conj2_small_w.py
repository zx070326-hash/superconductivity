r"""
P4 . KAY Conjecture 2 -- SEGMENT I of Lemma B (`theory/derivations/P4_conjecture2.md` §2)
=========================================================================================
Computer-assisted certificate (mpmath, dps=40) for the pointwise inequality

    h(wbar) := lam_1( H(wbar) )  <  frak_g(2) * wbar^2        for all wbar in (0, wbar_0],

with wbar_0 reported by Stage 4.  H(wbar) is the KAY II operator of
`p4_interchange_uniformity.H_matrix` / `H_mat` ([[k]](wbar) = wbar^2/(k^2+wbar^2)) and
frak_g(2) = rho(G(2)) is the top eigenvalue of the KAY I gamma-model operator
`G_matrix(N,2)`.  NOTHING here modifies a base certificate script: the verbatim builders
`H_matrix` / `H_mat` / `G_matrix` are IMPORTED and used as the correctness gate on the
(kernel-generic, vectorized) builders written here.

--------------------------------------------------------------------------------------
0.  THE ONE STRUCTURAL OBSERVATION
--------------------------------------------------------------------------------------
Every operator in this family is  O[f] := -L1[f] + L2[f] + L3[f]  for a positive kernel
f : {1,2,...} -> (0,inf):

    L1[f]_{nm} = delta_{nm} (1/(2n+1)) sum_{k=1}^{n} 2 f(k)
    L2[f]_{nm} = (1-delta_{nm}) f(|n-m|) / ( sqrt(2n+1) sqrt(2m+1) )
    L3[f]_{nm} =               f(n+m+1) / ( sqrt(2n+1) sqrt(2m+1) )

and f |-> O[f] is LINEAR, with O[f] entrywise-dominated by Oabs[f] := L1[f]+L2[f]+L3[f]
and Oabs monotone in f.  With

    G(gamma) = O[ k -> k^-gamma ],    H(wbar) = O[ k -> wbar^2/(k^2+wbar^2) ],
    K(wbar)  = O[ k -> 1/(k^2 (k^2+wbar^2)) ],

the elementary kernel identity  1/(k^2+w^2) = 1/k^2 - w^2/(k^2(k^2+w^2))  gives the EXACT
operator identity (valid for ALL wbar>0 -- no series, no |wbar|<1 restriction)

    (I)      H(wbar)/wbar^2  =  G(2)  -  wbar^2 K(wbar).

Expanding 1/(k^2+w^2) as a geometric series instead gives the note's series form
    (I')     H(wbar)/wbar^2 = sum_{j>=0} (-wbar^2)^j G(2j+2),      wbar < 1,
which Stage 0 verifies numerically and Stage 3 uses for a (weaker) cross-check route.
The delivered certificate of Stage 4 rests on (I) alone.

--------------------------------------------------------------------------------------
1.  THE PERTURBATION LEMMA  (self-contained proof; sharpened form of the remainder (*)
    at `p4_remainder_interchange.py` line ~40)
--------------------------------------------------------------------------------------
LEMMA.  Let A be bounded self-adjoint with lam_1 := max spec(A) a SIMPLE isolated
eigenvalue, unit eigenvector v0, P0 = |v0><v0|, and
    g := lam_1 - sup spec( A|_{v0^perp} )  > 0.

[Simplicity/isolation is not assumed but CERTIFIED, and without any compactness input:
 Stage 2 produces alpha with  alpha >= min_u sup_{x perp u, ||x||=1} <x,Ax>  and a trial
 vector with lamhat = <vhat,Avhat>/<vhat,vhat> > alpha.  If the spectral projection
 Q := 1_{(alpha,inf)}(A) had rank >= 2 then for EVERY unit u the subspace ran(Q) ∩ u^perp
 would be nonzero and would contain a unit x with <x,Ax> > alpha, i.e.
 min_u sup_{x perp u} <x,Ax> > alpha -- contradiction.  So rank Q <= 1; and Q != 0 since
 sup spec(A) >= lamhat > alpha.  Hence spec(A) ∩ (alpha,inf) is a single simple eigenvalue,
 necessarily lam_1, with eigenvector v0 = ran Q, and sup spec(A|_{v0^perp}) <= alpha.]

Let E be bounded self-adjoint and put
    e := || (I-P0) E v0 ||,     Etil := sup_{u perp v0, ||u||=1} <u,E u>.
For every unit w write w = c v0 + s u (u perp v0 unit, c^2+s^2=1, c>=0, real case).  Then

    <w,A w> = c^2 lam_1 + <(I-P0)w, A (I-P0)w>          <= lam_1 - s^2 g,
    <w,E w> = <v0,Ev0> + s^2(<u,Eu> - <v0,Ev0>) + 2cs<v0,Eu>
            <= <v0,Ev0> + s^2( Etil - <v0,Ev0> ) + 2 s e,

the last step by <v0,Eu> = <(I-P0)Ev0,u> <= e and c <= 1.  Taking the sup over w,

    lam_1(A+E) <= lam_1(A) + <v0,Ev0> + max_{s in [0,1]} [ s^2 beta + 2 s e ],      (*)
    beta := Etil - <v0,Ev0> - g.

No attainment of the sup is needed, so (*) holds for any bounded self-adjoint A+E.
With the coarse choices e <= ||E||, Etil <= ||E||, -<v0,Ev0> <= ||E|| (so beta <= 2||E||-g)
and beta < 0, the inner max is e^2/(-beta) and (*) degenerates EXACTLY to the quoted
    0 <= lam_1(A+E) - lam_1(A) - <v0,Ev0> <= ||E||^2/(g - 2||E||),
so (*) is a strict sharpening of the (*) of `p4_remainder_interchange.py` -- and that
sharpening is what buys wbar_0 >> 0.7:  e = ||(I-P0)Ev0|| is ~8x smaller than ||E||, and
Etil (which is  -wbar^2 lam_min(K),  NOT  +wbar^2||K||) is ~1.5x smaller than ||E||.

APPLICATION.  A := G(2),  E := -wbar^2 K(wbar)  (identity (I)).  Put
    kappa(wbar) := <v0, K(wbar) v0>  (>0),    Lam(wbar) := -lam_min(K(wbar))  (>0),
    ehat(wbar)  := || (I-P0) K(wbar) v0 ||.
Then <v0,Ev0> = -wbar^2 kappa, Etil <= wbar^2 Lam, e <= wbar^2 ehat, and (*) reads

    lam_1(H(wbar)/wbar^2) - frak_g(2)  <=  wbar^2 * max_{s in [0,1]} Psi(s),
    Psi(s) := -(1-s^2) kappa  +  s^2 ( Lam - g/wbar^2 )  +  2 s ehat.           (**)

*** frak_g(2) CANCELS.  The target h(wbar) < frak_g(2) wbar^2 is EXACTLY
    max_{s in [0,1]} Psi(s) < 0, which involves only kappa, Lam, ehat and a LOWER bound
    on the gap g.  The WIDTH of the frak_g(2) enclosure never enters, so there is NO
    small-wbar cutoff wbar_min: the s~0 end is governed by -kappa < 0 (degree-2
    homogeneous, so it survives wbar -> 0) and the s~1 end by wbar^2(Lam+2ehat) < g.
    This answers plan item (iv): we do NOT compare a certified upper bound of h/wbar^2
    against frak_g(2)_lo; we cancel the TRUE frak_g(2) on both sides, which is strictly
    stronger and removes the (0,wbar_min) gap entirely.  (As a by-product Stage 2 also
    tightens the frak_g(2) enclosure from P4-C's width 1.06e-3 to ~1e-11, which would
    have closed the small-wbar end of the naive formulation as well.)

--------------------------------------------------------------------------------------
2.  CERTIFIED INGREDIENTS  (mpmath everywhere; float64 appears ONLY as (a) trial vectors,
    (b) Cholesky shift anchors, (c) evidence tables -- never inside a bound)
--------------------------------------------------------------------------------------
 (a) SCHUR/HOLMGREN ENVELOPE  ||O[f]|| <= Rsum[f] := sup_n sum_m |O[f]_{nm}|, as exact
     mpmath partial sums for n<n0, m<M plus closed-form integral-test remainders from a
     majorant f(k) <= C_f/k^p.  Rsum is monotone in f, hence Rsum[k^-gamma] is
     NON-INCREASING in gamma: ONE evaluation at gamma0 bounds ||G(gamma)|| for ALL
     gamma >= gamma0.  That is the uniform majorant plan item (ii) asks for, and it makes
     the series tail sum_{j>=2} wbar^{2j}||G(2j+2)|| <= S8 wbar^6/(1-wbar^2) closed-form.
 (b) HEAD/COUPLING/TAIL SPLIT at cut Na:  O[f] = [[A,B],[B^T,D]] and
        lam_max(O) <= max(lam_max(A), ||D||) + ||B||,
       -lam_min(O) <= max(-lam_min(A), ||D||) + ||B||,
        ||O||      <= max(|lam_max(A)|, |lam_min(A)|, ||D||) + ||B||,
     with lam_max(A) / -lam_min(A) certified by mpmath Cholesky of (cI -+ A) and ||D||,
     ||B|| by closed-form Schur envelopes.
 (c) TRIAL VECTOR + DAVIS-KAHAN.  vhat = float64 top eigenvector of the NV-head of G(2)
     (EVIDENCE).  lamhat = <vhat,G(2)vhat>/<vhat,vhat> exact (finite support), so
     frak_g(2) >= lamhat by Rayleigh.  alpha >= lam_2(G(2)) by head Cholesky-deflation +
     Schur tail (same machine as P4-C).  g_lo := lamhat - alpha <= g.  The residual
     r = G(2)vhat - lamhat vhat is evaluated over rows n < M_RES exactly PLUS a
     closed-form bound on rows n >= M_RES -- i.e. the TAIL COUPLING of vhat to rows
     beyond its support is included.  That is the P4-D §2.4b Davis-Kahan fix (there the
     head-only residual was ~1e-16 and wrong); here the tail term is ~10^10 x the head
     term, so the fix is load-bearing and Stage 2 asserts it.
     Then sin angle(vhat,v0) <= s_hi := ||r|| / (lamhat - alpha)   (proof: with
     vhat = c v0 + s u, r = c(lam_1-lamhat)v0 + s(A-lamhat)u and ||(A-lamhat)u|| >=
     dist(lamhat, spec\{lam_1}) >= lamhat - alpha, so ||r|| >= s (lamhat-alpha)).
 (d) CELL-UNIFORM BOUNDS.  For wbar in [a,b]:  f_b <= f_wbar <= f_a pointwise, so with
     the difference kernel d := f_a - f_b >= 0,
        | K(wbar) - K(b) | = | O[f_wbar - f_b] | <= Oabs[d]   entrywise,
     and every quantity at wbar is its value at b plus an Oabs[d] penalty.  Combined with
     monotonicity of the nonneg-coefficient quadratic forms <vhat,L_i[f]vhat> this yields
     rigorous kappa_lo, Lam_hi, ehat_hi valid on the WHOLE cell.

--------------------------------------------------------------------------------------
3.  HONEST BOUNDARIES
--------------------------------------------------------------------------------------
 * mpmath is NOT interval arithmetic.  Every certified quantity is padded by
   ROUND = 1e-25 (absolute) in the safe direction; at dps=40 the accumulated rounding of
   the largest sum used here (~10^6 terms of O(1) magnitude) is <= 1e-33, eight orders
   below the pad.  The pad is invisible at the reported precision.
 * A successful floating-point Cholesky at dps=40 is taken as a positive-definiteness
   certificate (the convention of P4-B / P4-C / P4-D), with an explicit shift pad
   CHOL_PAD = 1e-12.
 * v0, the top eigenvector of the LIMIT operator G(2) on l^2, is NEVER computed.  Only
   vhat (finite support, exact float64 entries) and the certified angle s_hi enter.
 * The certificate is a THEOREM on (0, wbar_0] only.  wbar > wbar_0 is Segment II/III.

Run:  python theory/derivations/p4_conj2_small_w.py            (all stages)
      python theory/derivations/p4_conj2_small_w.py 0 1 2      (selected stages)
"""
import os
import sys
import time

import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# ---- base certificate scripts: IMPORTED ONLY (never modified) ---------------------
from p4_interchange_uniformity import H_mat, H_matrix, G_matrix          # noqa: E402
from p4_gamma_general import prop5_rhs                                   # noqa: E402

mp.mp.dps = 40

ROUND = mp.mpf("1e-25")        # safety pad on every certified quantity
CHOL_PAD = mp.mpf("1e-12")     # shift pad for the Cholesky SPD certificates

PASS = [0]
_STAGE0 = [0]
T_START = time.time()


def _check(name, ok, detail=""):
    if not ok:
        raise AssertionError(f"CHECK FAILED: {name}  {detail}")
    PASS[0] += 1


def _stage_summary(tag, note=""):
    k = PASS[0] - _STAGE0[0]
    _STAGE0[0] = PASS[0]
    print(f"  [{tag}] {k} checks passed.  {note}  (t={time.time()-T_START:.1f}s)")


# =====================================================================================
# STAGE 0 -- kernel-generic builders, gated against the verbatim base builders
# =====================================================================================

class Kernel:
    """Positive non-increasing kernel f(k), k>=1, with a rigorous majorant f(k)<=C/k^p
    and (where available) a CLOSED FORM for Sf = sum_{k>=1} f(k).

    Closed-form integral-test remainders (k^-p <= int_{k-1}^k x^-p dx, valid J>=2):
        sum_{k>=J} f(k)   <= C   / ((p-1)  (J-1)^(p-1)),
        sum_{k>=J} f(k)^2 <= C^2 / ((2p-1) (J-1)^(2p-1)).
    """

    def __init__(self, name, fun, C, p, fun_np=None, total_fn=None):
        self.name = name
        self._f = fun
        self.C = mp.mpf(C)
        self.p = mp.mpf(p)
        self._np = fun_np
        self._total_fn = total_fn
        self._cache = {}
        self._total = None

    def f(self, k):
        k = int(k)
        r = self._cache.get(k)
        if r is None:
            r = self._f(k)
            self._cache[k] = r
        return r

    def table(self, kmax):
        return [mp.mpf(0)] + [self.f(k) for k in range(1, kmax + 1)]

    def np_table(self, kmax):
        k = np.arange(0, kmax + 1, dtype=np.float64)
        out = np.zeros_like(k)
        out[1:] = self._np(k[1:])
        return out

    def total(self):
        """Upper bound on Sf = sum_{k>=1} f(k) (closed form when supplied)."""
        if self._total is None:
            if self._total_fn is not None:
                self._total = self._total_fn() + ROUND
            else:
                K = 4000
                s = mp.mpf(0)
                for k in range(1, K):
                    s += self.f(k)
                self._total = s + self.C / ((self.p - 1) * mp.mpf(K - 1) ** (self.p - 1)) + ROUND
        return self._total

    def tail_sum(self, J):
        J = int(J)
        if J <= 1:
            return self.total()
        return self.C / ((self.p - 1) * mp.mpf(J - 1) ** (self.p - 1)) + ROUND

    def tail_sq_sum(self, J):
        J = max(int(J), 2)
        return self.C ** 2 / ((2 * self.p - 1) * mp.mpf(J - 1) ** (2 * self.p - 1)) + ROUND


def _sum_inv_k2_plus_w2(w):
    """sum_{k>=1} 1/(k^2+w^2) = (pi/(2w)) coth(pi w) - 1/(2 w^2)   (exact, w>0)."""
    w = mp.mpf(w)
    return mp.pi / (2 * w) * mp.coth(mp.pi * w) - 1 / (2 * w ** 2)


def _sum_K(w):
    """sum_{k>=1} 1/(k^2(k^2+w^2)) = (zeta(2) - sum 1/(k^2+w^2))/w^2  (exact), w>=0."""
    w = mp.mpf(w)
    if w == 0:
        return mp.zeta(4)
    return (mp.zeta(2) - _sum_inv_k2_plus_w2(w)) / w ** 2


def ker_gamma(gamma):
    g = mp.mpf(gamma)
    gf = float(gamma)
    return Kernel(f"k^-{gamma}", lambda k: 1 / mp.mpf(k) ** g, 1, gamma,
                  fun_np=lambda k: np.exp(-gf * np.log(k)),
                  total_fn=lambda: mp.zeta(g))


def ker_K(w):
    """f(k) = 1/(k^2 (k^2+w^2)) -- the kernel of K(wbar).  f(k) <= k^-4."""
    wm = mp.mpf(w)
    w2 = wm ** 2
    wf = float(w)
    return Kernel(f"K({w})", lambda k: 1 / (mp.mpf(k) ** 2 * (mp.mpf(k) ** 2 + w2)), 1, 4,
                  fun_np=lambda k: 1.0 / (k ** 2 * (k ** 2 + wf * wf)),
                  total_fn=lambda: _sum_K(wm))


def ker_H(w):
    """f(k) = w^2/(k^2+w^2) -- the kernel of H(wbar).  f(k) <= w^2 k^-2."""
    wm = mp.mpf(w)
    w2 = wm ** 2
    wf = float(w)
    return Kernel(f"H({w})", lambda k: w2 / (mp.mpf(k) ** 2 + w2), w2, 2,
                  fun_np=lambda k: wf * wf / (k ** 2 + wf * wf),
                  total_fn=lambda: w2 * _sum_inv_k2_plus_w2(wm))


def ker_diff(a, b):
    """d(k) = f_a(k) - f_b(k) = (b^2-a^2)/(k^2 (k^2+a^2)(k^2+b^2)) >= 0,  d(k) <= (b^2-a^2)k^-6."""
    am, bm = mp.mpf(a), mp.mpf(b)
    a2, b2 = am ** 2, bm ** 2
    af, bf = float(a), float(b)
    dd = b2 - a2

    def fun(k):
        k2 = mp.mpf(k) ** 2
        return dd / (k2 * (k2 + a2) * (k2 + b2))

    def fun_np(k):
        k2 = k ** 2
        return (bf * bf - af * af) / (k2 * (k2 + af * af) * (k2 + bf * bf))

    return Kernel(f"d[{a},{b}]", fun, dd, 6, fun_np=fun_np,
                  total_fn=lambda: _sum_K(am) - _sum_K(bm))


def op_np(N, ftab, signs=(-1.0, 1.0, 1.0)):
    """float64 O[f] (or Oabs[f] with signs=(1,1,1)) from a float64 kernel table."""
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
    return signs[0] * np.diag(diag) + signs[1] * L2 + signs[2] * L3


def op_mp(N, ker, signs=(-1, 1, 1)):
    """exact mpmath O[f] head block of size N."""
    ft = ker.table(2 * N + 2)
    sq = [mp.sqrt(2 * n + 1) for n in range(N)]
    diag = []
    run = mp.mpf(0)
    for k in range(N):
        if k >= 1:
            run += 2 * ft[k]
        diag.append(run / (2 * k + 1))
    Mx = mp.matrix(N)
    s0, s1, s2 = [mp.mpf(s) for s in signs]
    for n in range(N):
        for m in range(N):
            v = s2 * ft[n + m + 1] / (sq[n] * sq[m])
            if n != m:
                v += s1 * ft[abs(n - m)] / (sq[n] * sq[m])
            else:
                v += s0 * diag[n]
            Mx[n, m] = v
    return Mx


def run_stage0():
    print("=== STAGE 0: kernel-generic builders, gated on the verbatim base builders ===")
    # (a) op_np with the gamma kernel reproduces the VERBATIM G_matrix
    for N in (4, 11, 30):
        for gam in (2, 3, 4, 6):
            ft = ker_gamma(gam).np_table(2 * N + 2)
            d = np.max(np.abs(op_np(N, ft) - G_matrix(N, gam)))
            _check("op_np==G_matrix", d <= 1e-14, f"N={N} gamma={gam} d={d:.3e}")
    # (b) op_np with the H kernel reproduces the VERBATIM H_matrix / H_mat
    for N in (4, 11, 30):
        for w in (0.1, 0.5, 1.0, 2.4, 7.0):
            ft = ker_H(w).np_table(2 * N + 2)
            d = np.max(np.abs(op_np(N, ft) - H_matrix(N, w)))
            _check("op_np==H_matrix", d <= 1e-14, f"N={N} w={w} d={d:.3e}")
            d2 = np.max(np.abs(op_np(N, ft) - H_mat(N, w)))
            _check("op_np==H_mat", d2 <= 1e-14, f"N={N} w={w} d={d2:.3e}")
    # (c) mpmath builder == float64 builder
    for N in (5, 25):
        for ker in (ker_gamma(2), ker_gamma(4), ker_K(0.7), ker_diff(0.4, 0.7)):
            A = op_mp(N, ker)
            B = op_np(N, ker.np_table(2 * N + 2))
            d = max(abs(float(A[i, j]) - B[i, j]) for i in range(N) for j in range(N))
            _check("op_mp==op_np", d <= 1e-13, f"N={N} ker={ker.name} d={d:.3e}")
    # (c2) closed-form kernel totals agree with brute-force partial sums
    KBR = 20000
    for ker in (ker_gamma(2), ker_gamma(4), ker_K("0.9"), ker_H("1.7"), ker_diff("0.5", "0.8")):
        brute = mp.mpf(0)
        for k in range(1, KBR):
            brute += ker.f(k)
        tot = ker.total()
        # brute (exact partial) <= Sf <= brute + integral-test tail; tot must lie in there
        _check("closed-form Sf is a valid upper bound", tot >= brute - ROUND,
               f"{ker.name}: {mp.nstr(tot,20)} < {mp.nstr(brute,20)}")
        _check("closed-form Sf matches the partial sum + tail",
               tot <= brute + ker.tail_sum(KBR) + ROUND,
               f"{ker.name}: {mp.nstr(tot,20)} vs {mp.nstr(brute,20)}"
               f"+{mp.nstr(ker.tail_sum(KBR),6)}")

    # (d) IDENTITY (I):  H(w)/w^2 = G(2) - w^2 K(w)   -- exact, all w>0
    N = 200
    for w in (0.05, 0.3, 0.7, 1.0, 2.4, 6.0, 30.0):
        Hn = op_np(N, ker_H(w).np_table(2 * N + 2)) / w ** 2
        Gn = op_np(N, ker_gamma(2).np_table(2 * N + 2))
        Kn = op_np(N, ker_K(w).np_table(2 * N + 2))
        d = np.max(np.abs(Hn - (Gn - w * w * Kn)))
        _check("identity (I) H/w^2 = G(2)-w^2 K(w)", d <= 1e-13, f"w={w} d={d:.3e}")
    for w in ("0.3", "1.7"):
        n0 = 12
        Hm, Gm2, Km = op_mp(n0, ker_H(w)), op_mp(n0, ker_gamma(2)), op_mp(n0, ker_K(w))
        w2 = mp.mpf(w) ** 2
        d = max(abs(Hm[i, j] / w2 - (Gm2[i, j] - w2 * Km[i, j]))
                for i in range(n0) for j in range(n0))
        _check("identity (I) exact in mpmath", d <= mp.mpf("1e-35"), f"w={w} d={mp.nstr(d,3)}")

    # (e) IDENTITY (I'):  H(w)/w^2 = sum_j (-w^2)^j G(2j+2)   for w<1  (plan item (i))
    N = 200
    Gcache = {}
    for w in (0.1, 0.2, 0.3, 0.4, 0.5):
        Hn = op_np(N, ker_H(w).np_table(2 * N + 2)) / w ** 2
        Sm = np.zeros_like(Hn)
        jmax = int(np.ceil(36.0 / (-np.log10(w ** 2)))) + 4
        for j in range(jmax):
            gam = 2 * j + 2
            if gam not in Gcache:
                Gcache[gam] = op_np(N, ker_gamma(gam).np_table(2 * N + 2))
            Sm = Sm + ((-w * w) ** j) * Gcache[gam]
        d = np.max(np.abs(Hn - Sm))
        _check("identity (I') series", d <= 1e-14, f"w={w} jmax={jmax} d={d:.3e}")
        print(f"    (I') at w={w}: max |H/w^2 - sum_j (-w^2)^j G(2j+2)| = {d:.2e}"
              f"   (N={N}, {jmax} terms)")

    # (f) G(gamma) -> G(inf) structure sanity: only |n-m|=1 and n=m=0 survive
    Ginf = op_np(6, ker_gamma(60).np_table(14))
    _check("G(gamma->inf) subdiagonal", abs(Ginf[0, 1] - 1 / np.sqrt(3)) < 1e-12)
    _check("G(gamma->inf) (0,0) entry", abs(Ginf[0, 0] - 1.0) < 1e-12)
    _stage_summary("STAGE 0")


# =====================================================================================
# STAGE 1 -- rigorous norm machinery:  Schur envelope + head/coupling/tail split
# =====================================================================================

def schur_rowsum_bound(ker, n0=100, M=1200):
    """Rigorous  ||O[f]|| <= Rsum[f] = sup_n sum_m |O[f]_{nm}|  (Schur/Holmgren).

    Using |O_nn| <= L1_nn + L3_nn,
      row_n <= (2/(2n+1)) sum_{k<=n} f(k)
               + sum_{m!=n} f(|n-m|)/(sqrt(2n+1)sqrt(2m+1))
               + sum_{m>=0} f(n+m+1)/(sqrt(2n+1)sqrt(2m+1)).
    n < n0 : exact m-sums for m < M plus the closed-form remainders
      sum_{m>=M} f(m-n)/(sqrt(2n+1)sqrt(2m+1))  <= tail_sum(M-n)  /(sqrt(2n+1)sqrt(2M+1)),
      sum_{m>=M} f(m+n+1)/(sqrt(2n+1)sqrt(2m+1))<= tail_sum(M+n+1)/(sqrt(2n+1)sqrt(2M+1)).
    n >= n0 (bulk), using sqrt(2m+1) >= 1 and sum_{m!=n} f(|n-m|) <= 2 Sf :
      row_n <= 2 Sf/(2n0+1) + 2 Sf/sqrt(2n0+1) + tail_sum(n0+1)/sqrt(2n0+1).
    Monotone in f (f <= h pointwise => Rsum[f] <= Rsum[h]); in particular
    Rsum[k^-gamma] is non-increasing in gamma.
    """
    ft = ker.table(M + n0 + 4)
    sq = [mp.sqrt(2 * m + 1) for m in range(M + 1)]
    Sf = ker.total()
    best = mp.mpf(0)
    run = mp.mpf(0)
    for n in range(n0):
        if n >= 1:
            run += 2 * ft[n]
        row = run / (2 * n + 1)
        acc = mp.mpf(0)
        for m in range(M):
            if m != n:
                acc += ft[abs(n - m)] / sq[m]
            acc += ft[n + m + 1] / sq[m]
        acc += (ker.tail_sum(M - n) + ker.tail_sum(M + n + 1)) / sq[M]
        row += acc / sq[n]
        if row > best:
            best = row
    bulk = (2 * Sf / (2 * n0 + 1) + 2 * Sf / mp.sqrt(2 * n0 + 1)
            + ker.tail_sum(n0 + 1) / mp.sqrt(2 * n0 + 1))
    return max(best, bulk) + ROUND


def tailblock_bound(ker, Na):
    """||D|| for the tail block D = O[f] restricted to n,m >= Na.
    Row n >= Na:  L1 <= 2Sf/(2Na+1);  L2 <= 2Sf/(2Na+1) (both sqrt's >= sqrt(2Na+1));
    L3 <= tail_sum(2Na+1)/(2Na+1).  =>  ||D|| <= (4 Sf + tail_sum(2Na+1))/(2Na+1)."""
    return (4 * ker.total() + ker.tail_sum(2 * Na + 1)) / (2 * mp.mpf(Na) + 1) + ROUND


def coupling_bound(ker, Na):
    """||B|| for the head(<Na) x tail(>=Na) block, via ||B|| <= sqrt(rowmax*colmax).
      row n < Na : sum_{m>=Na}[f(m-n)+f(m+n+1)]/(sqrt(2n+1)sqrt(2m+1))
                   <= [tail_sum(Na-n)+tail_sum(Na+n+1)]/(sqrt(2n+1)sqrt(2Na+1))
      col m >= Na: every term is non-increasing in m (f non-increasing, 1/sqrt(2m+1)
                   decreasing), so the max over m>=Na is at m=Na and is computed there."""
    sqNa = mp.sqrt(2 * mp.mpf(Na) + 1)
    rowmax = mp.mpf(0)
    for n in range(Na):
        r = (ker.tail_sum(Na - n) + ker.tail_sum(Na + n + 1)) / (mp.sqrt(2 * n + 1) * sqNa)
        if r > rowmax:
            rowmax = r
    colmax = mp.mpf(0)
    for n in range(Na):
        colmax += (ker.f(Na - n) + ker.f(Na + n + 1)) / mp.sqrt(2 * n + 1)
    colmax /= sqNa
    return mp.sqrt(rowmax * colmax) + ROUND


def _chol_ok(Mx):
    try:
        mp.cholesky(Mx)
        return True
    except Exception:          # noqa: BLE001
        return False


def head_spec_bounds(ker, Na):
    """(lam_max(A) upper bound, -lam_min(A) upper bound) for the head block A of O[f],
    certified by mpmath Cholesky:  cholesky(cI - A) OK => lam_max(A) < c ;
                                   cholesky(cI + A) OK => lam_min(A) > -c."""
    A_np = op_np(Na, ker.np_table(2 * Na + 2))
    ev = np.linalg.eigvalsh(A_np)
    A = op_mp(Na, ker)
    out = []
    for sgn, anchor in ((+1, float(ev[-1])), (-1, -float(ev[0]))):
        c = mp.mpf(anchor) + CHOL_PAD
        for _ in range(60):
            Mx = mp.matrix(Na)
            for i in range(Na):
                for j in range(Na):
                    Mx[i, j] = (c if i == j else mp.mpf(0)) - sgn * A[i, j]
            if _chol_ok(Mx):
                break
            c = c + abs(c) * mp.mpf("1e-9") + mp.mpf("1e-12")
        else:
            raise AssertionError(f"Cholesky escalation failed for {ker.name} sgn={sgn}")
        out.append(c + ROUND)
    return out[0], out[1]


def op_spec_bounds(ker, Na):
    """Certified (lam_max_hi, -lam_min_hi, norm_hi, ||D||, ||B||) for the FULL O[f]."""
    lmaxA, mlminA = head_spec_bounds(ker, Na)
    nD = tailblock_bound(ker, Na)
    nB = coupling_bound(ker, Na)
    return (max(lmaxA, nD) + nB, max(mlminA, nD) + nB,
            max(lmaxA, mlminA, nD) + nB, nD, nB)


NA_NORM = 120          # head cut for the K(wbar) / G(gamma) norm certificates
G_NORM = {}


def run_stage1():
    print("=== STAGE 1: rigorous operator-norm machinery (plan item (ii)) ===")
    for gam in (2, 4, 6, 8, 10):
        R = schur_rowsum_bound(ker_gamma(gam), n0=80, M=900)
        Gn = op_np(500, ker_gamma(gam).np_table(1100))
        true = float(np.max(np.abs(np.linalg.eigvalsh(Gn))))
        _check("Rsum[k^-g] >= float64 ||G(gamma)||", float(R) > true,
               f"gamma={gam} R={mp.nstr(R,8)} true={true:.6f}")
        G_NORM[gam] = R
        print(f"    Rsum[k^-{gam}] = {mp.nstr(R, 9)}    (float64 ||G({gam})|| = {true:.6f})")
    for g1, g2 in ((2, 4), (4, 6), (6, 8), (8, 10)):
        _check("Rsum non-increasing in gamma", G_NORM[g1] > G_NORM[g2],
               f"{mp.nstr(G_NORM[g1],8)} vs {mp.nstr(G_NORM[g2],8)}")
    for gam in (12, 20, 40, 80):
        R = schur_rowsum_bound(ker_gamma(gam), n0=40, M=400)
        _check("Rsum[k^-8] majorizes all gamma>=8", R <= G_NORM[8] + ROUND,
               f"gamma={gam}: {mp.nstr(R,8)} vs S8={mp.nstr(G_NORM[8],8)}")
    p5, eps = prop5_rhs(2)
    kay_chain = 2 * p5 + 1
    print(f"    KAY Prop5/6 chain at gamma=2: rho(G(2)) <= 2*Prop5+1 = {mp.nstr(kay_chain,8)}"
          f"  (eps={mp.nstr(eps,4)});  Schur envelope here = {mp.nstr(G_NORM[2],8)}")
    _check("Schur envelope beats the KAY chain at gamma=2", G_NORM[2] < kay_chain)

    for gam in (4, 6):
        lmax, mlmin, nrm, nD, nB = op_spec_bounds(ker_gamma(gam), NA_NORM)
        Gn = op_np(600, ker_gamma(gam).np_table(1300))
        ev = np.linalg.eigvalsh(Gn)
        _check("head/tail lam_max(G) valid", float(lmax) > float(ev[-1]))
        _check("head/tail -lam_min(G) valid", float(mlmin) > -float(ev[0]))
        _check("head/tail ||G|| valid", float(nrm) > float(max(abs(ev[0]), abs(ev[-1]))))
        _check("head/tail ||G|| beats the Schur envelope", nrm < G_NORM[gam])
        G_NORM[f"cert{gam}"] = nrm
        print(f"    ||G({gam})|| <= {mp.nstr(nrm, 10)}  (float64 "
              f"{max(abs(ev[0]),abs(ev[-1])):.7f});  ||D||={mp.nstr(nD,3)} ||B||={mp.nstr(nB,3)}")
    _stage_summary("STAGE 1", f"S6={mp.nstr(G_NORM[6],7)}  S8={mp.nstr(G_NORM[8],7)}")


# =====================================================================================
# STAGE 2 -- frak_g(2) lower bound, alpha, certified gap, trial vector, Davis-Kahan
# =====================================================================================

NV = 200               # trial-vector support
M_RES = 2400           # row cut for the G(2) residual (its kernel decays only like k^-2)
M_VEC = NV + 250       # row cut for the K(wbar) matvecs (kernel decays like k^-4)
NA_ALPHA = 300         # head cut for the lam_2(G(2)) Cholesky-deflation certificate

S = {}


def qform_parts(ker):
    """(<vhat,L1[f]vhat>, <vhat,(L2+L3)[f]vhat>)/<vhat,vhat>, EXACT via the precomputed
    structure coefficients.  Both have NONNEGATIVE coefficients, hence both are monotone
    increasing in the kernel f -- the fact that makes the cell-uniform bounds work."""
    P, T, U = S["P"], S["T"], S["U"]
    neg = mp.mpf(0)
    for k in range(1, NV):
        neg += ker.f(k) * P[k]
    pos = mp.mpf(0)
    for j in range(1, NV):
        pos += ker.f(j) * T[j]
    for s in range(1, 2 * NV):
        pos += ker.f(s) * U[s]
    return neg / S["nrm2"], pos / S["nrm2"]


def qform(ker):
    neg, pos = qform_parts(ker)
    return pos - neg


def matvec(ker, M, signs=(-1, 1, 1)):
    """y_n = (O[f] vhat)_n for n < M (for the UNNORMALISED vhat) and a rigorous bound on
    sum_{n>=M} (y_n/nu)^2.  For n >= NV the L1 term drops (vhat_n = 0) and
      |y_n| <= (1/sqrt(2n+1)) * 2 f(n-NV+1) * sum_m u_m
    since f(|n-m|) <= f(n-NV+1) and f(n+m+1) <= f(n+1) <= f(n-NV+1) for m < NV; hence
      sum_{n>=M} y_n^2 <= (4 Cu^2/(2M+1)) sum_{j>=M-NV+1} f(j)^2."""
    v, u = S["v"], S["u"]
    ft = ker.table(M + NV + 4)
    s0, s1, s2 = [mp.mpf(s) for s in signs]
    y = []
    run = mp.mpf(0)
    for n in range(M):
        if n >= 1:
            run += 2 * ft[n]
        acc = mp.mpf(0)
        for m in range(NV):
            if m != n:
                acc += s1 * ft[abs(n - m)] * u[m]
            acc += s2 * ft[n + m + 1] * u[m]
        yn = acc / mp.sqrt(2 * n + 1)
        if n < NV:
            yn += s0 * (run / (2 * n + 1)) * v[n]
        y.append(yn)
    tail2 = (4 * S["Cu"] ** 2 / (2 * mp.mpf(M) + 1)) * ker.tail_sq_sum(M - NV + 1)
    return y, tail2 / S["nrm2"] + ROUND


def alpha_G2(Na):
    """alpha >= lam_2(G(2)).  For ANY unit u,
         lam_2(A) = min_{dim S=1} max_{x perp S} <x,Ax> <= max_{x perp u} <x,Ax>
                  = lam_max(PAP|_{u^perp}) <= lam_max(PAP),   P = I - uu^T,
    so a Cholesky certificate for cI - PAP gives lam_2(A) < c on the HEAD block A.
    Blockwise lam_2(G) <= lam_2(blockdiag(A,D)) + ||B|| <= max(lam_2(A),||D||) + ||B||."""
    ker = ker_gamma(2)
    A_np = op_np(Na, ker.np_table(2 * Na + 2))
    u = np.linalg.eigh(A_np)[1][:, -1]
    Au_f = A_np @ u
    PAP_f = A_np - np.outer(u, Au_f) - np.outer(Au_f, u) + float(u @ Au_f) * np.outer(u, u)
    c_anchor = float(np.linalg.eigvalsh(PAP_f)[-1])

    A = op_mp(Na, ker)
    umb = mp.matrix([mp.mpf(float(x)) for x in u])
    nrm = mp.sqrt(sum(umb[i] ** 2 for i in range(Na)))
    for i in range(Na):
        umb[i] = umb[i] / nrm
    Au = A * umb
    uAu = sum(umb[i] * Au[i] for i in range(Na))
    c = mp.mpf(c_anchor) + CHOL_PAD
    for _ in range(80):
        Mx = mp.matrix(Na)
        for i in range(Na):
            for j in range(Na):
                pap = A[i, j] - umb[i] * Au[j] - Au[i] * umb[j] + umb[i] * umb[j] * uAu
                Mx[i, j] = (c if i == j else mp.mpf(0)) - pap
        if _chol_ok(Mx):
            break
        c = c + mp.mpf("1e-10")
    else:
        raise AssertionError("alpha_G2: Cholesky escalation failed")
    nD = tailblock_bound(ker, Na)
    nB = coupling_bound(ker, Na)
    return max(c, nD) + nB + ROUND, c, nD, nB


def run_stage2():
    print("=== STAGE 2: frak_g(2) lower bound, alpha, gap, trial vector, Davis-Kahan ===")
    ker2 = ker_gamma(2)

    A_np = op_np(NV, ker2.np_table(2 * NV + 2))
    vh = np.linalg.eigh(A_np)[1][:, -1]
    vh = vh * np.sign(vh[0])
    _check("trial vector is entrywise positive (Perron structure)", float(np.min(vh)) > 0,
           f"min={np.min(vh):.3e}")
    v = [mp.mpf(float(x)) for x in vh]
    nrm2 = sum(x * x for x in v)
    S["v"], S["nrm2"], S["nu"] = v, nrm2, mp.sqrt(nrm2)

    u = [v[n] / mp.sqrt(2 * n + 1) for n in range(NV)]
    S["u"] = u
    S["Cu"] = sum(u)
    P = [mp.mpf(0)] * (NV + 1)
    acc = mp.mpf(0)
    for n in range(NV - 1, 0, -1):
        acc += v[n] ** 2 / (2 * n + 1)
        P[n] = 2 * acc
    T = [mp.mpf(0)] * (NV + 1)
    for j in range(1, NV):
        s = mp.mpf(0)
        for n in range(NV - j):
            s += u[n] * u[n + j]
        T[j] = 2 * s
    U = [mp.mpf(0)] * (2 * NV + 1)
    for n in range(NV):
        for m in range(NV):
            U[n + m + 1] += u[n] * u[m]
    S["P"], S["T"], S["U"] = P, T, U
    for arr, nm in ((P, "P"), (T, "T"), (U, "U")):
        _check(f"structure coefficients {nm} >= 0", all(x >= 0 for x in arr))

    for ker in (ker_gamma(2), ker_gamma(4), ker_K("1.3")):
        q_fast = qform(ker)
        Amp = op_mp(NV, ker)
        q_slow = mp.mpf(0)
        for n in range(NV):
            row = mp.mpf(0)
            for m in range(NV):
                row += Amp[n, m] * v[m]
            q_slow += v[n] * row
        q_slow /= nrm2
        _check("structure-coefficient quadratic form == direct double sum",
               abs(q_fast - q_slow) < mp.mpf("1e-30"),
               f"{ker.name}: {mp.nstr(q_fast,20)} vs {mp.nstr(q_slow,20)}")

    lamhat = qform(ker2)
    S["lamhat"] = lamhat

    t0 = time.time()
    alpha, c_defl, nD2, nB2 = alpha_G2(NA_ALPHA)
    S["alpha"] = alpha
    print(f"    alpha >= lam_2(G(2)):  {mp.nstr(alpha, 10)}   (deflated head c="
          f"{mp.nstr(c_defl,6)}, ||D||={mp.nstr(nD2,4)}, ||B||={mp.nstr(nB2,4)}; "
          f"Na={NA_ALPHA}, {time.time()-t0:.1f}s)")
    _check("alpha <= P4-C published 2.0e-2", alpha <= mp.mpf("2.0e-2"),
           f"alpha={mp.nstr(alpha,8)}")
    _check("alpha < lamhat (Kato-Temple / Davis-Kahan hypothesis)", alpha < lamhat)

    # rho(G(2)) = lam_1(G(2)):  certify that the negative end cannot win
    _, mlmin2, _, _, _ = op_spec_bounds(ker2, NA_NORM)
    S["mlmin_G2"] = mlmin2
    _check("rho(G(2)) = lam_1(G(2)) (i.e. -lam_min < lamhat)", mlmin2 < lamhat,
           f"-lam_min(G2) <= {mp.nstr(mlmin2,8)} vs lamhat={mp.nstr(lamhat,8)}")

    y, ytail2 = matvec(ker2, M_RES)
    r2 = mp.mpf(0)
    for n in range(M_RES):
        rn = y[n] - lamhat * (v[n] if n < NV else mp.mpf(0))
        r2 += rn * rn
    r2 = r2 / nrm2 + ytail2 + ROUND
    rnorm = mp.sqrt(r2)
    r_head2 = mp.mpf(0)
    for n in range(NV):
        rn = y[n] - lamhat * v[n]
        r_head2 += rn * rn
    r_head2 /= nrm2
    S["rnorm"] = rnorm
    print(f"    ||r_full|| = {mp.nstr(rnorm, 8)}   (head-only ||r_head|| = "
          f"{mp.nstr(mp.sqrt(r_head2),4)};  the P4-D §2.4b tail-coupling term is "
          f"{mp.nstr(rnorm/mp.sqrt(r_head2+mp.mpf('1e-80')),4)}x larger)")
    _check("tail coupling dominates the head residual (P4-D §2.4b fix is live)",
           r_head2 < r2 / 4)

    g_lo = lamhat - alpha
    s_hi = rnorm / (lamhat - alpha) + ROUND
    S["g_lo"], S["s_hi"] = g_lo, s_hi
    _check("s_hi < 1e-3 (Davis-Kahan angle is tight)", s_hi < mp.mpf("1e-3"))

    kt_up = lamhat + r2 / (lamhat - alpha) + ROUND
    S["g2_lo"], S["g2_hi"] = lamhat, kt_up
    print(f"    frak_g(2) in [{mp.nstr(lamhat, 16)}, {mp.nstr(kt_up, 16)}]    width "
          f"{mp.nstr(kt_up-lamhat, 4)}")
    # P4-B anchor: C_inf = 0.1827262477... = sqrt(frak_g(2))/(2 pi)  =>  frak_g(2) =
    # 1.31814021108...  (NOTE: the string "1.3181405" that appears in
    # P4_rigorous_tc_bound.md §P4-C and in p4_gamma_general.G_TRUE_REF_2 is a rounding
    # slip -- the 7th digit is 2, not 5; p4_C_infinity_enclosure.py line ~657 itself
    # writes "g(2) = 1.31814021...".  Reported, not silently absorbed.)
    g2_from_Cinf = (2 * mp.pi * mp.mpf("0.1827262477")) ** 2
    _check("enclosure contains the P4-B anchor (2 pi C_inf)^2",
           lamhat <= g2_from_Cinf + mp.mpf("1e-9") and kt_up >= g2_from_Cinf,
           f"[{mp.nstr(lamhat,16)},{mp.nstr(kt_up,16)}] vs {mp.nstr(g2_from_Cinf,16)}")
    _check("enclosure inside the P4-C bracket [1.31714,1.31820]",
           lamhat >= mp.mpf("1.31714") and kt_up <= mp.mpf("1.31820"))
    _check("enclosure width beats P4-C's 1.06e-3", kt_up - lamhat < mp.mpf("1e-6"))
    print(f"    g_lo = lamhat - alpha = {mp.nstr(g_lo, 12)}     "
          f"s_hi = sin angle(vhat,v0) <= {mp.nstr(s_hi, 6)}")
    _stage_summary("STAGE 2")


# =====================================================================================
# STAGE 3/4 shared -- the cell-uniform certified bundle
# =====================================================================================

NODE = {}


def node(w):
    """All wbar-local certified quantities for K(wbar) at a single wbar >= 0."""
    key = mp.nstr(mp.mpf(w), 30)
    if key in NODE:
        return NODE[key]
    ker = ker_K(w)
    qneg, qpos = qform_parts(ker)
    q = qpos - qneg
    y, tail2 = matvec(ker, M_VEC)
    v, nrm2 = S["v"], S["nrm2"]
    ny2 = eh2 = mp.mpf(0)
    for n in range(M_VEC):
        vn = v[n] if n < NV else mp.mpf(0)
        ny2 += y[n] * y[n]
        d = y[n] - q * vn
        eh2 += d * d
    ny = mp.sqrt(ny2 / nrm2 + tail2 + ROUND)
    ehv = mp.sqrt(eh2 / nrm2 + tail2 + ROUND)
    lmax, mlmin, nrm, nD, nB = op_spec_bounds(ker, NA_NORM)
    out = dict(qneg=qneg, qpos=qpos, q=q, ny=ny, ehv=ehv, Lam=mlmin, NK=nrm,
               lmax=lmax, nD=nD, nB=nB)
    NODE[key] = out
    return out


def cell_data(a, b):
    """Cell-uniform certified bundle valid for every wbar in [a,b]  (a == b allowed).

    kappa  = <v0,K v0>,  with v0 = c vhat + s utilde (utilde unit perp vhat, s <= s_hi):
        kappa = c^2 qhat + 2 c s <vhat,K utilde> + s^2 <utilde,K utilde>
      => qhat_lo - pad <= kappa <= qhat_hi + pad,   pad := s^2(|qhat_hi|+N) + 2 s ||K vhat||
    ehat = ||(I-P0) K v0||:   K v0 - kappa v0 = c(Kvhat - qhat vhat) + c(qhat-kappa)vhat
                                              + s(K utilde - kappa utilde), and
      |qhat - kappa| <= pad, so  ehat <= ||(I-P0hat)K vhat|| + pad + s(N + kappa_hi).
    Cell-uniformity: f_b <= f_wbar <= f_a, d := f_a - f_b, |K(wbar)-K(b)| <= Oabs[d],
      so  ||K(wbar)vhat|| <= ||K(b)vhat|| + ||Oabs[d]vhat||,
          ||(I-P0hat)K(wbar)vhat|| <= ||(I-P0hat)K(b)vhat|| + ||Oabs[d]vhat||,
          ||K(wbar)|| <= ||K(b)|| + Rsum[d],   -lam_min(K(wbar)) <= -lam_min(K(b)) + Rsum[d];
      and the L1 / (L2+L3) quadratic forms are monotone in f, so qhat_lo takes L1 at f_a
      and (L2+L3) at f_b (and vice versa for qhat_hi)."""
    nb = node(b)
    na = node(a) if mp.mpf(a) != mp.mpf(b) else nb
    s_hi = S["s_hi"]
    if mp.mpf(a) != mp.mpf(b):
        kd = ker_diff(a, b)
        yW, tW = matvec(kd, M_VEC, signs=(1, 1, 1))
        wv2 = mp.mpf(0)
        for n in range(M_VEC):
            wv2 += yW[n] * yW[n]
        Wv = mp.sqrt(wv2 / S["nrm2"] + tW + ROUND)
        Wn = schur_rowsum_bound(kd, n0=60, M=600)
    else:
        Wv = Wn = mp.mpf(0)
    N_hi = nb["NK"] + Wn
    ny_hi = nb["ny"] + Wv
    qhat_lo = nb["qpos"] - na["qneg"]
    qhat_hi = na["qpos"] - nb["qneg"]
    pad = s_hi * s_hi * (abs(qhat_hi) + N_hi) + 2 * s_hi * ny_hi
    kap_lo = qhat_lo - pad - ROUND
    kap_hi = qhat_hi + pad + ROUND
    eh_hi = nb["ehv"] + Wv + pad + s_hi * (N_hi + kap_hi) + ROUND
    Lam_hi = nb["Lam"] + Wn
    return dict(kap_lo=kap_lo, kap_hi=kap_hi, qhat_lo=qhat_lo, qhat_hi=qhat_hi,
                eh_hi=eh_hi, Lam_hi=Lam_hi, N_hi=N_hi, ny_hi=ny_hi, Wv=Wv, Wn=Wn)


def cell_bound(a, b):
    """Certified  sup_{wbar in [a,b]} max_{s in [0,1]} Psi(s)   (see (**)).
    < 0 certifies  h(wbar) < frak_g(2) wbar^2  on [a,b].
    Psi(s) <= Psibar(s) = -kap_lo + s^2 (kap_lo + Lam_hi - g_lo/b^2) + 2 s eh_hi,
    using wbar <= b (so g/wbar^2 >= g_lo/b^2) and 1-s^2 >= 0."""
    dd = cell_data(a, b)
    b2 = mp.mpf(b) ** 2
    bt = dd["kap_lo"] + dd["Lam_hi"] - S["g_lo"] / b2
    if bt >= 0:
        val, branch = -dd["kap_lo"] + bt + 2 * dd["eh_hi"], "s=1"
    else:
        st = dd["eh_hi"] / (-bt)
        if st <= 1:
            val, branch = -dd["kap_lo"] + dd["eh_hi"] ** 2 / (-bt), "interior"
        else:
            val, branch = -dd["kap_lo"] + bt + 2 * dd["eh_hi"], "s=1"
    dd["bt"] = bt
    dd["branch"] = branch
    return val, dd


# =====================================================================================
# STAGE 3 -- ROUTE A: the note's series formulation (cross-check, weaker wbar_0)
# =====================================================================================

def route_A_F(w, c2_lo, N4, N6, S8, g_lo):
    """F(w) := -c2_lo + w^2 N4^2/(g_lo - 2 w^2 N4) + w^2 N6 + S8 w^4/(1-w^2).
    From (I') with A = G(2), E = -w^2 G(4), R = sum_{j>=2}(-w^2)^j G(2j+2):
      h/w^2 <= lam_1(G(2) - w^2 G(4)) + ||R||
            <= frak_g(2) - w^2 c2 + w^4 N4^2/(g - 2 w^2 N4) + w^4 N6 + S8 w^6/(1-w^2),
    so h/w^2 - frak_g(2) <= w^2 F(w) and F(w) < 0 certifies the target on (0,w]."""
    w2 = mp.mpf(w) ** 2
    den = g_lo - 2 * w2 * N4
    if den <= 0 or w2 >= 1:
        return mp.mpf("+inf")
    return -c2_lo + w2 * N4 ** 2 / den + w2 * N6 + S8 * w2 ** 2 / (1 - w2)


def run_stage3():
    print("=== STAGE 3: ROUTE A -- the note's series/Weyl formulation (cross-check) ===")
    _, _, N4, _, _ = op_spec_bounds(ker_gamma(4), NA_NORM)
    N6, S8, g_lo = G_NORM["cert6"], G_NORM[8], S["g_lo"]
    dd = cell_data(mp.mpf(0), mp.mpf(0))          # K(0) = G(4)
    c2_lo, c2_hi, qhat4 = dd["kap_lo"], dd["kap_hi"], dd["qhat_lo"]
    S["c2_lo"] = c2_lo
    print(f"    <vhat,G(4)vhat> = {mp.nstr(qhat4, 12)};  c2 = <v0,G(4)v0> in "
          f"[{mp.nstr(c2_lo,12)}, {mp.nstr(c2_hi,12)}]")
    _check("c2 ~ 1.19317 (the note's float64 probe value)",
           abs(qhat4 - mp.mpf("1.19317")) < mp.mpf("1e-4"))
    _check("c2_lo > 1.19", c2_lo > mp.mpf("1.19"))
    print(f"    ||G(4)|| <= {mp.nstr(N4,10)};  ||G(6)|| <= {mp.nstr(N6,10)};  "
          f"S8 = sup_{{gamma>=8}}||G(gamma)|| <= {mp.nstr(S8,10)};  g_lo = {mp.nstr(g_lo,10)}")

    lo, hi = mp.mpf("0.01"), mp.mpf("0.99")
    _check("route A: F(0.01) < 0", route_A_F(lo, c2_lo, N4, N6, S8, g_lo) < 0)
    if route_A_F(hi, c2_lo, N4, N6, S8, g_lo) < 0:
        w0A = hi
    else:
        for _ in range(80):
            mid = (lo + hi) / 2
            if route_A_F(mid, c2_lo, N4, N6, S8, g_lo) < 0:
                lo = mid
            else:
                hi = mid
        w0A = lo
    S["w0A"] = w0A
    prev = mp.mpf("-inf")
    for i in range(1, 41):
        wv = w0A * i / 40
        Fv = route_A_F(wv, c2_lo, N4, N6, S8, g_lo)
        _check("route A: F < 0 on the grid", Fv < 0, f"w={mp.nstr(wv,6)} F={mp.nstr(Fv,6)}")
        _check("route A: F increasing in wbar", Fv > prev)
        prev = Fv
    print(f"    ROUTE A certifies  h(wbar) < frak_g(2) wbar^2  on (0, {mp.nstr(w0A, 8)}]")
    print(f"      budget at wbar_0^A: 2nd-order remainder "
          f"{mp.nstr(w0A**2*N4**2/(g_lo-2*w0A**2*N4),6)} + ||G(6)|| term "
          f"{mp.nstr(w0A**2*N6,6)} + series tail {mp.nstr(S8*w0A**4/(1-w0A**2),6)}"
          f"  vs  c2_lo = {mp.nstr(c2_lo,6)}")
    _stage_summary("STAGE 3", f"wbar_0^A = {mp.nstr(w0A,6)}")


# =====================================================================================
# STAGE 4 -- ROUTE B: the exact-K cell ladder (the delivered certificate)
# =====================================================================================

W_TARGET = mp.mpf("6.0")
D_MAX = mp.mpf("0.5")
D_MIN = mp.mpf("0.004")
R_MIN = mp.mpf("0.05")     # a cell is only accepted if its margin is >= 5% of kappa_lo


def _next_width(d, margin, kap_lo):
    """Adapt the next cell width from the achieved relative margin.  Purely a scheduling
    heuristic (the certificate is whatever the cells prove); tuned to keep every cell's
    relative margin above ~5%, so that no cell in the ladder is marginal."""
    r = margin / kap_lo
    if r > mp.mpf("0.35"):
        return d * mp.mpf("1.8")
    if r > mp.mpf("0.18"):
        return d * mp.mpf("1.3")
    if r > mp.mpf("0.10"):
        return d * mp.mpf("1.05")
    if r > mp.mpf("0.06"):
        return d
    return d * mp.mpf("0.75")


def run_stage4():
    print("=== STAGE 4: ROUTE B -- exact-K cell ladder (the delivered certificate) ===")
    print("      cell                kappa_lo      Lam_hi      ehat_hi   branch     margin")
    w, d = mp.mpf(0), D_MAX
    cells = []
    while w < W_TARGET:
        b = min(w + d, W_TARGET)
        val, dg = cell_bound(w, b)
        if val < 0 and -val >= R_MIN * dg["kap_lo"]:
            cells.append((w, b, val, dg))
            print(f"    [{mp.nstr(w,5):>7},{mp.nstr(b,5):>7}] {mp.nstr(dg['kap_lo'],7):>13}"
                  f" {mp.nstr(dg['Lam_hi'],7):>11} {mp.nstr(dg['eh_hi'],6):>10}"
                  f" {dg['branch']:>9} {mp.nstr(-val,6):>10}")
            _check("cell certificate margin > 0", val < 0)
            _check("cell kappa_lo > 0", dg["kap_lo"] > 0)
            _check("cell ehat_hi >= 0", dg["eh_hi"] >= 0)
            _check("cell Lam_hi > 0", dg["Lam_hi"] > 0)
            w = b
            d = min(D_MAX, _next_width(d, -val, dg["kap_lo"]))
        else:
            d = d / 2
            if d < D_MIN:
                print(f"    LADDER STOPS at wbar={mp.nstr(w,8)}  (cell width < "
                      f"{mp.nstr(D_MIN,3)}; last trial value {mp.nstr(val,6)})")
                break
    S["w0"], S["cells"] = w, cells
    _check("wbar_0 >= 0.5 (task target)", w >= mp.mpf("0.5"))
    _check("wbar_0 >= 0.3 (Segment II hand-off floor)", w >= mp.mpf("0.3"))
    if "w0A" in S:
        _check("route B reaches at least as far as route A", w >= S["w0A"])
    a_t, b_t, v_t, dg_t = cells[-1]
    i_tight = min(range(len(cells)), key=lambda i: -cells[i][2])
    a_g, b_g, v_g, dg_g = cells[i_tight]
    # the s=1 branch condition  wbar^2 (Lam + 2 ehat) < g  is the structural constraint;
    # report the worst cell for it (it is what ultimately terminates the ladder).
    worst = max(cells, key=lambda c: c[1] ** 2 * (c[3]["Lam_hi"] + 2 * c[3]["eh_hi"]))
    lhs_w = worst[1] ** 2 * (worst[3]["Lam_hi"] + 2 * worst[3]["eh_hi"])
    print(f"\n    ROUTE B certifies  h(wbar) < frak_g(2) wbar^2  on (0, {mp.nstr(w, 8)}]  "
          f"with {len(cells)} cells.")
    print(f"    tightest cell [{mp.nstr(a_g,6)},{mp.nstr(b_g,6)}]: margin {mp.nstr(-v_g,6)}"
          f"  (= {mp.nstr(-v_g/dg_g['kap_lo']*100,4)}% of kappa_lo), branch {dg_g['branch']}")
    print(f"    BINDING TERM (the structural one, s=1 branch = directions ORTHOGONAL to v0):")
    print(f"      max over cells of  wbar^2 (Lam + 2 ehat)  =  {mp.nstr(lhs_w,8)}   at cell "
          f"[{mp.nstr(worst[0],6)},{mp.nstr(worst[1],6)}]")
    print(f"      must stay < g_lo = {mp.nstr(S['g_lo'],8)}      [slack "
          f"{mp.nstr(S['g_lo']-lhs_w,6)} = {mp.nstr((S['g_lo']-lhs_w)/S['g_lo']*100,4)}%]")
    _check("worst-cell s=1 slack > 3% of g_lo",
           (S["g_lo"] - lhs_w) / S["g_lo"] > mp.mpf("0.03"),
           f"slack={mp.nstr((S['g_lo']-lhs_w)/S['g_lo']*100,5)}%")
    print(f"    top cell [{mp.nstr(a_t,6)},{mp.nstr(b_t,6)}]: branch={dg_t['branch']}, "
          f"kappa_lo={mp.nstr(dg_t['kap_lo'],8)}, Lam_hi={mp.nstr(dg_t['Lam_hi'],8)}, "
          f"ehat_hi={mp.nstr(dg_t['eh_hi'],8)}, margin={mp.nstr(-v_t,6)}")
    S["bind_lhs"] = lhs_w

    # ---- independent float64 validation sweep (EVIDENCE, not part of any bound) -----
    # Recompute kappa, Lam, ehat DIRECTLY from a large float64 truncation with the
    # float64 top eigenvector of G(2) at N=NBIG, and confirm every certified bound sits
    # on the correct side.  This is the gate that would catch a sign/derivation slip.
    NBIG = 700
    G2n = op_np(NBIG, ker_gamma(2).np_table(2 * NBIG + 2))
    evn, Vn = np.linalg.eigh(G2n)
    v0n = Vn[:, -1] * np.sign(Vn[0, -1])
    print("\n    independent float64 validation (N=700, true v0; EVIDENCE only):")
    print("      wbar     kappa_lo <= kappa      Lam_hi >= Lam      ehat_hi >= ehat")
    for wv in ("0.3", "0.8", "1.5", "2.4", "4.0", "6.0"):
        wf = float(wv)
        if mp.mpf(wv) > S["w0"]:
            continue
        Kn = op_np(NBIG, ker_K(wf).np_table(2 * NBIG + 2))
        kap_t = float(v0n @ Kn @ v0n)
        Kv = Kn @ v0n
        eh_t = float(np.sqrt(max(0.0, Kv @ Kv - kap_t * kap_t)))
        Lam_t = -float(np.linalg.eigvalsh(Kn)[0])
        dd = cell_data(mp.mpf(wv), mp.mpf(wv))
        _check("kappa_lo <= float64 kappa", float(dd["kap_lo"]) <= kap_t,
               f"w={wv}: {mp.nstr(dd['kap_lo'],10)} vs {kap_t:.10f}")
        _check("kappa_hi >= float64 kappa", float(dd["kap_hi"]) >= kap_t, f"w={wv}")
        _check("Lam_hi >= float64 -lam_min(K)", float(dd["Lam_hi"]) >= Lam_t,
               f"w={wv}: {mp.nstr(dd['Lam_hi'],10)} vs {Lam_t:.10f}")
        _check("ehat_hi >= float64 ehat", float(dd["eh_hi"]) >= eh_t,
               f"w={wv}: {mp.nstr(dd['eh_hi'],10)} vs {eh_t:.10f}")
        print(f"    {wv:>6}   {mp.nstr(dd['kap_lo'],9):>10} <= {kap_t:<12.9f} "
              f"{mp.nstr(dd['Lam_hi'],9):>10} >= {Lam_t:<12.9f} "
              f"{mp.nstr(dd['eh_hi'],9):>10} >= {eh_t:<12.9f}")
    _stage_summary("STAGE 4", f"wbar_0 = {mp.nstr(w,8)}")


# =====================================================================================
# STAGE 5 -- evidence table: certified upper bound on h/wbar^2 vs float64 truth
# =====================================================================================

def run_stage5():
    print("=== STAGE 5: evidence table (plan item (v)) ===")
    g2hi, g2lo = S["g2_hi"], S["g2_lo"]
    print("     wbar   certified UB on h/wbar^2    float64 h/wbar^2 (N=400)     UB-truth"
          "     UB/frak_g(2)")
    Nbig = 400
    for wv in ("0.05", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "1.0",
               "1.5", "2.0", "2.4", "3.0", "4.0", "5.0", "6.0"):
        w = mp.mpf(wv)
        if w > S["w0"]:
            continue
        val, _ = cell_bound(w, w)
        ub = g2hi + w ** 2 * val
        Hn = op_np(Nbig, ker_H(float(wv)).np_table(2 * Nbig + 2))
        truth = float(np.linalg.eigvalsh(Hn)[-1]) / float(w) ** 2
        _check("certified UB >= float64 truth", float(ub) > truth,
               f"w={wv} ub={mp.nstr(ub,10)} truth={truth:.8f}")
        _check("certified UB < frak_g(2)_lo", ub < g2lo, f"w={wv}")
        print(f"    {wv:>5}   {mp.nstr(ub, 12):>22}   {truth:>22.10f}   "
              f"{float(ub)-truth:>10.6f}   {float(ub)/float(g2lo):>12.8f}")
    for wv in ("1e-3", "1e-6"):
        w = mp.mpf(wv)
        val, _ = cell_bound(w, w)
        _check("Psi_max < 0 down to wbar=1e-6", val < 0, f"w={wv} val={mp.nstr(val,6)}")
        if "c2_lo" in S:
            _check("Psi_max -> -c2 as wbar->0", abs(val + S["c2_lo"]) < mp.mpf("1e-4"),
                   f"w={wv} val={mp.nstr(val,10)} c2_lo={mp.nstr(S['c2_lo'],10)}")
    if "c2_lo" in S:
        print(f"    wbar -> 0 limit of max_s Psi(s)  =  -c2  =  -{mp.nstr(S['c2_lo'],10)}"
              f"   (< 0 : no wbar_min gap)")
    _stage_summary("STAGE 5")


# =====================================================================================
if __name__ == "__main__":
    stages = sys.argv[1:] or ["0", "1", "2", "3", "4", "5"]
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
    print("\n" + "=" * 78)
    print("SYNTHESIS -- P4 KAY conjecture 2, Segment I")
    if "w0" in S:
        print(f"  CERTIFIED:  h(wbar) < frak_g(2) wbar^2  for ALL wbar in "
              f"(0, {mp.nstr(S['w0'],8)}]")
        print("              no wbar_min gap -- the frak_g(2) enclosure WIDTH never enters,")
        print("              frak_g(2) cancels on both sides of (**).")
    if "w0A" in S:
        print(f"  cross-check, route A (the note's series/Weyl plan): (0, "
              f"{mp.nstr(S['w0A'],6)}]")
    if "g2_lo" in S:
        print(f"  by-product: frak_g(2) in [{mp.nstr(S['g2_lo'],16)}, "
              f"{mp.nstr(S['g2_hi'],16)}]   (P4-C: [1.31714, 1.31820])")
    print(f"  ALL {PASS[0]} checks passed  ({time.time()-T_START:.1f}s).")
    print("=" * 78)
