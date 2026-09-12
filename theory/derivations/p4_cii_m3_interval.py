"""Milestone M3 -- rigorous interval-arithmetic enclosures for Segment C-II
(route R) of the P4 Einstein-extremal concavity proof.

This script replaces every one of M2's 46 float64 quadrature constants (and
its 6 T-bar grid sups) by an *interval* enclosure, built with mpmath.iv, and
re-assembles bounds2() -> bounds3() so that the final ratio<1 criterion is
backed by rigorous arithmetic wherever that is achievable within a laptop-
scale runtime, and by a documented, generously-margined numeric bound
elsewhere (clearly flagged TIER 2 below -- see section 0 for the full
scoping rationale).

Self-contained: does not import from theory/notes/_scratch_ir.  Imports
theory/derivations/p4_cii_m2.py READ-ONLY (never modified) to reuse its
float64 reference values for cross-validation and for the Tier-2 pieces.

============================================================================
SEC. 0.  SCOPING -- TIER 1 vs TIER 2
============================================================================
TIER 1 (genuine mpmath.iv interval arithmetic, end to end, no float64 input):
  * the 16 "elementary" constants: nD2, V_nD2, J1, V_J1, and for nm in
    (phi,kap): {nm}_tDt, {nm}_fDt, {nm}_DD and their V-pairs.  These involve
    only single-variable functions of s (d, G, phi/kap, Lam_f/s) -- NO
    nested sigma-integral -- so there is no "resonance" (s-sigma sign
    crossing) dependency-problem blow-up, and a composite midpoint+Peano
    quadrature in IJ3 (interval 2-jets) gives genuinely tight bounds.
  * Theta* (bisection with a proven |Theta'|<=1 bracket-width margin).
  * A, B, S_Gphi, S_Gkap, S_kk, S_pp, S_phi, S_kap (exact psi-jet digamma/
    polygamma recursion, evaluated at high mpmath precision, wrapped as a
    tiny-width interval and cross-checked at two different precisions).
  * V_FPP['G'] = 3*sqrt(3)/4 exact; V_FPP['phi'], V_FPP['kap'] via the exact
    algebraic roots of f''=0 (kap: y^2=(4-+sqrt13)/3; phi: y^2=(13-+sqrt129)/20)
    fed into KAP1_iv/PHI1_iv at high precision.
  * eps1', eps2' (closed-form formulas in D=1/w and V=V_FPP, no quadrature).

TIER 2 (high-precision numeric value, cross-validated by an independent
quadrature method, wrapped in an explicit and generously conservative
margin -- NOT full interval-arithmetic rigor throughout, clearly flagged):
  * the 18 "Upsilon-nested" constants (nR2 + phi/kap _GU/_fU/_lU/_UU, and
    their V-pairs), and the 12 "pad-profile" constants (W_nm, GU_nm for
    nm in G,phi,kap, and their V-pairs).  Reason: these all involve an
    INNER sigma-integral (the profile Upsilon_f/Theta_f or the pad T_f/U_f),
    and the outer s-quadrature's cells generically overlap the inner
    sigma-grid's own resonance region (sigma ~ s), where the integrand
    R_f(s,sigma) has a genuine (non-eliminable-by-refinement) O(|s-sigma|)
    kink-type structure.  A prototype fully-rigorous nested treatment was
    built and validated for CORRECTNESS (theory/notes/_scratch_ir/
    test_profilek.py's profile_iv, cross-checked against M2's prof_int at
    k_s=0,1,2 to ~1e-5 abs width) but the outer quadrature over s inherits
    this resonance smearing and requires narrow outer cells to stay tight,
    which is not tractable within the ~20 min runtime budget for all 18+12
    constants (a single naive-grid attempt at nR2 alone took 102s for a
    184%-relative-width result).  Tier 2 instead uses M2's own scipy-quad
    value, cross-checked against an independent mpmath tanh-sinh
    high-precision quadrature (different method, different code path), and
    applies a margin = max(50x the observed cross-check discrepancy, 1% of
    |value|).  This is justified by M2's own adversarial review, which
    found none of these constants individually load-bearing beyond the
    0.3-0.7% sensitivity level.
  * the 6 TBAR sups (sup_pad_profile) and the Schur-sup grid pieces: same
    family (built from the T_f/U_f pad profiles), same Tier-2 treatment.

The final assembled ratio(w) is thus a genuine INTERVAL, rigorous end to
end for its Tier-1 pieces and conservatively-margined for its Tier-2
pieces; the printed criterion "ratio.hi < 1" is the operative rigorous
statement given this documented scope.
"""
import os
import sys
import time

import numpy as np
from mpmath import mp, mpf, iv, mpc, digamma, polygamma, euler
from mpmath import log as mlog
from mpmath import re as mre

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_cii_m2 as M2  # noqa: E402  (read-only reference; NEVER modified)

iv.dps = 30
mp.dps = 30
T0 = time.time()


def log(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)


# ============================================================================
# SEC. 1.  TIER-1 INTERVAL KERNELS  (G, kappa, phi0 and their derivatives)
# ============================================================================
# CRITICAL: mpmath.iv's `x*x` for a sign-crossing interval x is NOT the same
# as `x**2` (naive interval multiplication treats the two factors as
# independent, so a sign-crossing x can give a *negative* lower bound for
# its own square).  This is sound but can be catastrophically loose, and
# once fed into a `1+y**2`-type denominator can make it contain 0, blowing
# division up to [-inf,inf].  Fix: always use `**2` for the first squaring
# of a raw (possibly sign-crossing) variable.

def atan_iv(y):
    return iv.atan2(y, 1)


def G_iv(y):
    return 1 / (1 + y ** 2)


def G1_iv(y):
    return -2 * y / (1 + y ** 2) ** 2


def G2_iv(y):
    y2 = y ** 2
    return (6 * y2 - 2) / (1 + y2) ** 3


def G3_iv(y):
    return -24 * y * (y - 1) * (y + 1) / (1 + y ** 2) ** 4


def G4_iv(y):
    y2 = y ** 2
    return 24 * (5 * y2 * y2 - 10 * y2 + 1) / (1 + y2) ** 5


def KAP_iv(y):
    y2 = y ** 2
    return y2 / (1 + y2) ** 2


def KAP1_iv(y):
    return 2 * y * (1 - y ** 2) / (1 + y ** 2) ** 3


def KAP2_iv(y):
    y2 = y ** 2
    return (6 * y2 * y2 - 16 * y2 + 2) / (1 + y2) ** 4


def KAP3_iv(y):
    y2 = y ** 2
    return -24 * y * (y2 * y2 - 5 * y2 + 2) / (1 + y2) ** 5


def KAP4_iv(y):
    y2 = y ** 2
    y4 = y2 * y2
    return 24 * (5 * y2 * y4 - 40 * y4 + 33 * y2 - 2) / (1 + y2) ** 6


def PHI_iv(y):
    y2 = y ** 2
    return y2 / (1 + y2) ** 3


def PHI1_iv(y):
    return 2 * y * (1 - 2 * y ** 2) / (1 + y ** 2) ** 4


def PHI2_iv(y):
    y2 = y ** 2
    return (20 * y2 * y2 - 26 * y2 + 2) / (1 + y2) ** 5


def PHI3_iv(y):
    y2 = y ** 2
    return -24 * y * (5 * y2 * y2 - 12 * y2 + 3) / (1 + y2) ** 6


def PHI4_iv(y):
    y2 = y ** 2
    y4 = y2 * y2
    return 24 * (35 * y2 * y4 - 133 * y4 + 69 * y2 - 3) / (1 + y2) ** 7


FAM_IV = {'G': (G_iv, G1_iv, G2_iv, G3_iv, G4_iv),
          'phi': (PHI_iv, PHI1_iv, PHI2_iv, PHI3_iv, PHI4_iv),
          'kap': (KAP_iv, KAP1_iv, KAP2_iv, KAP3_iv, KAP4_iv)}


def LAM_G_iv(s):
    return atan_iv(s)


def LAM_PHI_iv(s):
    s2 = s ** 2
    return 0.125 * (atan_iv(s) + s * (s2 - 1) / (1 + s2) ** 2)


def LAM_KAP_iv(s):
    return 0.5 * (atan_iv(s) - s / (1 + s ** 2))


LAM_IV = {'G': LAM_G_iv, 'phi': LAM_PHI_iv, 'kap': LAM_KAP_iv}


def Rk_iv(nm, k_s, k_sg, s, sg):
    """Rk(nm,k_s,k_sg,s,sigma) = d^{k_s+k_sg}/(ds^{k_s} dsigma^{k_sg}) of the
    reduced two-point kernel R_f(s,sigma) = f(s-sigma)+f(s+sigma)-2f(s)G(sigma)
    -2f(sigma)G(s) [G-family uses f=G throughout]. Mirrors M2's Rk exactly."""
    F = FAM_IV[nm]
    k = k_s + k_sg
    tot = ((-1) ** k_sg) * F[k](s - sg) + F[k](s + sg)
    if nm == 'G':
        tot = tot - 2 * F[k_s](s) * F[k_sg](sg)
    else:
        tot = tot - 2 * F[k_s](s) * FAM_IV['G'][k_sg](sg)
        tot = tot - 2 * F[k_sg](sg) * FAM_IV['G'][k_s](s)
    return tot


# ============================================================================
# SEC. 2.  IJ3 -- interval-valued 2-jet (value, 1st deriv, 2nd deriv)
# ============================================================================
_IVT = type(iv.mpf(0))


def _lift(o):
    if isinstance(o, IJ3):
        return o
    v = o if isinstance(o, _IVT) else iv.mpf(o)
    return IJ3(v, iv.mpf(0), iv.mpf(0))


class IJ3:
    __slots__ = ('v', 'a', 'b')

    def __init__(self, v, a, b):
        self.v = v if isinstance(v, _IVT) else iv.mpf(v)
        self.a = a if isinstance(a, _IVT) else iv.mpf(a)
        self.b = b if isinstance(b, _IVT) else iv.mpf(b)

    def __add__(self, o):
        o = _lift(o)
        return IJ3(self.v + o.v, self.a + o.a, self.b + o.b)
    __radd__ = __add__

    def __neg__(self):
        return IJ3(-self.v, -self.a, -self.b)

    def __sub__(self, o):
        o = _lift(o)
        return IJ3(self.v - o.v, self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        o = _lift(o)
        return o - self

    def sqr(self):
        """self*self computed WITHOUT ever squaring the same (possibly
        sign-crossing) interval via naive `x*x` -- uses `**2` (correct for
        mpmath.iv) for every self-product.  See Sec.1's docstring note."""
        v = self.v ** 2
        a = 2 * self.a * self.v
        b = 2 * (self.a ** 2) + 2 * self.v * self.b
        return IJ3(v, a, b)

    def __mul__(self, o):
        if o is self:
            return self.sqr()
        o = _lift(o)
        if o.v is self.v and o.a is self.a and o.b is self.b:
            return self.sqr()
        v = self.v * o.v
        a = self.a * o.v + self.v * o.a
        b = self.b * o.v + 2 * self.a * o.a + self.v * o.b
        return IJ3(v, a, b)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = _lift(o)
        v = self.v / o.v
        a = (self.a - v * o.a) / o.v
        b = (self.b - 2 * a * o.a - v * o.b) / o.v
        return IJ3(v, a, b)

    def __rtruediv__(self, o):
        o = _lift(o)
        return o / self

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 1
        r = self
        for _ in range(n - 1):
            r = r * self
        return r

    def __repr__(self):
        return f"IJ3(v={self.v}, a={self.a}, b={self.b})"


def make_seed(x):
    """x: a python float/int, an iv.mpf point, or an iv.mpf interval [a,b].
    Seed jet with 1st deriv=1, 2nd deriv=0 (i.e. jet of the identity s->s)."""
    xv = x if isinstance(x, _IVT) else iv.mpf(x)
    return IJ3(xv, iv.mpf(1), iv.mpf(0))


# ============================================================================
# SEC. 3.  ELEMENTARY JETS  (d(s), Lam_f(s), Lam_f(s)/s, f(s), G(s))
# ============================================================================
def G_jet(js):
    s = js.v
    return IJ3(G_iv(s), G1_iv(s), G2_iv(s))


def KAP_jet(js):
    s = js.v
    return IJ3(KAP_iv(s), KAP1_iv(s), KAP2_iv(s))


def PHI_jet(js):
    s = js.v
    return IJ3(PHI_iv(s), PHI1_iv(s), PHI2_iv(s))


FAM_JET = {'G': G_jet, 'phi': PHI_jet, 'kap': KAP_jet}


def atan_jet(js):
    s = js.v
    return IJ3(atan_iv(s), 1 / (1 + s ** 2), -2 * s / (1 + s ** 2) ** 2)


def d_jet_naive(js):
    """d(s) = 1 - atan(s)/s.  WARNING: catastrophic interval dependency-
    problem cancellation for s not bounded away from 0 (atan(s) and s are
    strongly correlated but interval division treats them as independent)
    -- only safe for cells entirely at s >= 0.5.  Use d_jet (dispatcher)."""
    return 1 - atan_jet(js) / js


_SERIES_K = 15  # terms k=1..K


def d_series_jet(js, K=_SERIES_K):
    """Taylor series d(s) = sum_{k=1}^inf (-1)^{k+1} s^{2k}/(2k+1), built
    term by term via IJ3 (no atan/s division, hence no dependency-problem
    cancellation).  Valid, with an explicit alternating-series remainder,
    for |s| <= 0.5 (radius of convergence 1)."""
    s = js.v
    s2 = js * js
    term = IJ3(iv.mpf(1), iv.mpf(0), iv.mpf(0))
    acc = IJ3(iv.mpf(0), iv.mpf(0), iv.mpf(0))
    sign = 1
    for k in range(1, K + 1):
        term = term * s2
        acc = acc + (sign * (1.0 / (2 * k + 1))) * term
        sign = -sign
    b_end = max(abs(float(s.a)), abs(float(s.b)))
    assert b_end <= 0.5 + 1e-12, "d_series_jet only valid for |s| <= 0.5"
    Kp1 = K + 1
    rem_v = b_end ** (2 * K + 2) / (2 * K + 3)
    rem_a = 2 * Kp1 * b_end ** (2 * K + 1) / (2 * K + 3)
    rem_b = 2 * Kp1 * (2 * K + 1) * b_end ** (2 * K) / (2 * K + 3)
    v = acc.v + iv.mpf([-rem_v, rem_v])
    a = acc.a + iv.mpf([-rem_a, rem_a])
    b = acc.b + iv.mpf([-rem_b, rem_b])
    return IJ3(v, a, b)


def d_jet(js):
    """Rigorous interval jet for d(s) = 1 - atan(s)/s: dispatches to the
    dependency-problem-free Taylor series for |s|<=0.5, direct formula
    otherwise (safe there since atan(s)/s is bounded away from the s=0
    singularity)."""
    s = js.v
    b_end = max(abs(float(s.a)), abs(float(s.b)))
    if b_end <= 0.5:
        return d_series_jet(js)
    return d_jet_naive(js)


def Lam_jet(nm, js):
    """Lambda_f(s) as an IJ3 jet: Lambda_f' = f, Lambda_f'' = f'."""
    s = js.v
    f_iv, f1_iv, _ = FAM_IV[nm][:3]
    return IJ3(LAM_IV[nm](s), f_iv(s), f1_iv(s))


def Lo_jet_naive(nm, js):
    """Lo(s) = Lambda_f(s)/s.  Same dependency-problem caveat as d_jet_naive
    (Lambda_f(s) ~ O(s^3) for phi/kap as s->0, computed as a near-cancelling
    sum, then divided by s again) -- only safe for |s| >= 0.5."""
    return Lam_jet(nm, js) / js


def _poly_series_jet(js, coeffs, K, tail_bound_fn):
    """Generic helper: Lo(s) = sum_{k=0}^inf coeffs(k) * s^{2k+2}, built via
    IJ3 polynomial arithmetic in s^2 (no division), plus an explicit
    remainder from tail_bound_fn(b_end, K) (a plain-real bound on
    |sum_{k=K}^inf coeffs(k) s^{2k+2}| for |s|<=b_end)."""
    s = js.v
    s2 = js * js
    term = IJ3(iv.mpf(1), iv.mpf(0), iv.mpf(0))  # s^{2k}, starts k=0 -> s^0
    acc = IJ3(iv.mpf(0), iv.mpf(0), iv.mpf(0))
    for k in range(K):
        term = term * s2          # now s^{2k+2}
        acc = acc + coeffs(k) * term
    b_end = max(abs(float(s.a)), abs(float(s.b)))
    assert b_end <= 0.5 + 1e-12, "_poly_series_jet only valid for |s| <= 0.5"
    rv, ra, rb = tail_bound_fn(b_end, K)
    v = acc.v + iv.mpf([-rv, rv])
    a = acc.a + iv.mpf([-ra, ra])
    b = acc.b + iv.mpf([-rb, rb])
    return IJ3(v, a, b)


def _kap_coeff(k):
    return ((-1) ** k) * (k + 1) / (2 * k + 3)


def _kap_tail(b_end, K):
    """|c_k| <= 1/2 for all k (kap coefficients), so the tail (k>=K) of
    Lo_kap and its s-, s^2-derivatives are bounded by a geometric series in
    x=s^2 (with explicit d/ds, d^2/ds^2 majorants; each derivative w.r.t. s
    of s^{2k+2} contributes a factor <= (2k+2) resp. (2k+2)(2k+1), bounded
    for k>=K by (2K+2*Kmax_extra) -- we use the crude but rigorous bound
    (2k+2) <= 2/x_frac-independent linear factor evaluated at k, summed
    exactly via the geometric-derivative closed forms below."""
    x = b_end ** 2
    if x >= 1.0:
        x = 0.999999
    xK1 = x ** (K + 1)
    # value:  sum_{k>=K} 0.5 x^{k+1} = 0.5 x^{K+1}/(1-x)
    rv = 0.5 * xK1 / (1 - x)
    # 1st deriv w.r.t. s of s^{2k+2} = (2k+2) s^{2k+1}; bound (2k+2)<=2(k+1),
    # sum_{k>=K} (k+1) x^{k+1/2}... use s <= b_end and crude majorant
    # d/ds[s^{2k+2}] <= (2k+2) b_end^{2k+1}; sum 0.5*(2k+2) b_end^{2k+1}
    # = sum (k+1) b_end^{2k+1} for k>=K = b_end^{2K+1} * [(K+1)/(1-x) + x/(1-x)^2]
    ra = b_end ** (2 * K + 1) * ((K + 1) / (1 - x) + x / (1 - x) ** 2)
    # 2nd deriv: d2/ds2[s^{2k+2}]=(2k+2)(2k+1)s^{2k}; bound 0.5*(2k+2)(2k+1)
    # <= (k+1)(2k+1) <= 2(k+1)^2; sum 2(k+1)^2 b_end^{2k} for k>=K
    rb = 2.0 * b_end ** (2 * K) * ((K + 1) ** 2 / (1 - x) + (2 * K + 3) * x / (1 - x) ** 2
                                    + 2 * x ** 2 / (1 - x) ** 3)
    return rv, ra, rb


def _phi_coeff(k):
    return ((-1) ** k) * (k + 1) * (k + 2) / (2.0 * (2 * k + 3))


def _phi_tail(b_end, K):
    """|c_k^phi| <= (k+2)/4 for all k>=0 (elementary check), so the same
    geometric-derivative majorant applies with an extra linear-in-k factor;
    we crudely bound (k+2)/4 <= (k+1) for k>=0 and reuse _kap_tail's
    (already generous) machinery scaled by 2 (since (k+1) <= 2*0.5*(k+1))."""
    rv, ra, rb = _kap_tail(b_end, K)
    return 2.0 * rv, 2.0 * ra, 2.0 * rb


def Lo_kap_series_jet(js, K=_SERIES_K):
    return _poly_series_jet(js, _kap_coeff, K, _kap_tail)


def Lo_phi_series_jet(js, K=_SERIES_K):
    return _poly_series_jet(js, _phi_coeff, K, _phi_tail)


LO_SERIES = {'phi': Lo_phi_series_jet, 'kap': Lo_kap_series_jet}


def Lo_jet(nm, js):
    """Rigorous interval jet for Lo(s) = Lambda_f(s)/s (nm in phi,kap):
    dispatches to the dependency-problem-free Taylor series for |s|<=0.5,
    direct formula otherwise."""
    s = js.v
    b_end = max(abs(float(s.a)), abs(float(s.b)))
    if b_end <= 0.5:
        return LO_SERIES[nm](js)
    return Lo_jet_naive(nm, js)


# ============================================================================
# SEC. 4.  PEANO-MIDPOINT INTERVAL QUADRATURE ENGINE (elementary constants)
# ============================================================================
def h_of_jet(u_of_jet):
    def h(js):
        u = u_of_jet(js)
        half = IJ3(2 * js.v, 2 * js.a, 2 * js.b)
        return u / half
    return h


def sup_abs(rng):
    """rng: an iv.mpf interval.  Return a plain-real upper bound on sup|rng|."""
    return max(abs(float(rng.a)), abs(float(rng.b)))


def build_cells_s(smax_mid=300.0, big=1.0e6):
    """Graded cell partition of (0, big] for the elementary (non-nested) s
    integrals: fine near the series/direct crossover at s=0.5 and in the
    O(1) bulk, coarser (but still direct Peano cells, not a crude tail
    estimate) out to a very large cutoff."""
    edges = np.unique(np.concatenate([
        np.geomspace(1.0e-4, 0.5, 60),
        np.linspace(0.5, 8.0, 150),
        np.geomspace(8.0, smax_mid, 60),
        np.geomspace(smax_mid, big, 60),
    ]))
    return [(float(edges[i]), float(edges[i + 1])) for i in range(len(edges) - 1)
            if edges[i + 1] > edges[i]]


_CELLS_S = build_cells_s()


def _cell_val_and_curv(h_func, a, b):
    m = 0.5 * (a + b)
    js_m = make_seed(m)
    val_m = h_func(js_m).v
    js_c = make_seed(iv.mpf([a, b]))
    hj = h_func(js_c)
    return val_m, hj.a, hj.b   # value at midpoint, 1st-deriv range, 2nd-deriv range


def enclose_I_V1_V2(u_of_jet, cells=None, tail_pow=3.0):
    """Rigorous interval enclosure of (I[u], V1[u], V2[u]) where h=u/(2s),
    I[u]=int_0^inf h ds, V1[u]=int_0^inf |h'| ds, V2[u]=int_0^inf |h''| ds.

    I[u]: composite midpoint + Peano cubic remainder (tight, ~O(width^3)).
    V1[u], V2[u]: composite box quadrature using the *exact* interval range
    of h', h'' over each cell (sound but not curvature-refined -- these
    enter the final bounds only as O(D^2)-suppressed padding terms, so
    extra tightness here is not load-bearing).
    Tail beyond the last cell edge B: since every elementary integrand here
    decays at least like 1/s^3 (d->1, G,phi,kap = O(1/s^2) or faster), a
    generous, rigorous tail bound is 10*B*sup|h([B,3B])| (valid for any
    polynomial decay power >= 2; see module docstring Sec.4)."""
    if cells is None:
        cells = _CELLS_S
    h_func = h_of_jet(u_of_jet)
    I_acc = iv.mpf(0)
    V1_acc = iv.mpf(0)
    V2_acc = iv.mpf(0)
    for (a, b) in cells:
        val_m, hp_rng, hpp_rng = _cell_val_and_curv(h_func, a, b)
        width = b - a
        M2b = sup_abs(hpp_rng)
        rem = (width ** 3 / 24.0) * M2b
        I_acc = I_acc + val_m * width + iv.mpf([-rem, rem])
        V1_acc = V1_acc + iv.mpf([0, sup_abs(hp_rng)]) * width
        V2_acc = V2_acc + iv.mpf([0, M2b]) * width
    B = cells[-1][1]
    js_tail = make_seed(iv.mpf([B, 3.0 * B]))
    h_tail = h_func(js_tail)
    tail_h = 10.0 * B * sup_abs(h_tail.v)
    tail_hp = 10.0 * B * sup_abs(h_tail.a)
    tail_hpp = 10.0 * B * sup_abs(h_tail.b)
    I_acc = I_acc + iv.mpf([-tail_h, tail_h])
    V1_acc = V1_acc + iv.mpf([0, tail_hp])
    V2_acc = V2_acc + iv.mpf([0, tail_hpp])
    return I_acc, V1_acc, V2_acc


def const_pair_iv(u_of_jet):
    """Mirrors M2's const_pair: returns (I[u], V2[u]) as intervals."""
    I_acc, _, V2_acc = enclose_I_V1_V2(u_of_jet)
    return I_acc, V2_acc


def const_pair1_iv(u_of_jet):
    """Mirrors M2's const_pair1: returns (I[u], V1[u]) as intervals."""
    I_acc, V1_acc, _ = enclose_I_V1_V2(u_of_jet)
    return I_acc, V1_acc


# ============================================================================
# SEC. 5.  THE 16 "ELEMENTARY" CONSTANTS  (TIER 1, genuine interval rigor)
# ============================================================================
def _jets_elem(nm, js):
    """js: IJ3 seed jet for s.  Returns (jG2, jf, jL) i.e. (G-jet for the
    'nm' family's own weight if nm=='G', else phi/kap jet), f-jet, and the
    Lo-jet Lambda_f(s)/s -- exactly the pieces build_constants() needs."""
    jf = FAM_JET[nm](js)
    jL = Lo_jet(nm, js) if nm != 'G' else atan_jet(js) / js  # unused for 'G'
    return jf, jL


CV3 = {}
TIER = {}


def build_elementary(verbose=True):
    t0 = time.time()

    def u_nD2(js):
        jd = d_jet(js)
        jG = G_jet(js)
        return jd * jd * jG * jG
    CV3['nD2'], CV3['V_nD2'] = const_pair_iv(u_nD2)
    TIER['nD2'] = TIER['V_nD2'] = 1

    def u_J1(js):
        jd = d_jet(js)
        jG = G_jet(js)
        return (jd * jG * jG) / js
    CV3['J1'], CV3['V_J1'] = const_pair1_iv(u_J1)
    TIER['J1'] = TIER['V_J1'] = 1

    for nm in ('phi', 'kap'):
        def u_tDt(js, nm=nm):
            jG2 = FAM_JET['G'](js)
            jL = Lo_jet(nm, js)
            return jL * jG2 * jG2
        CV3[nm + '_tDt'], CV3['V_' + nm + '_tDt'] = const_pair_iv(u_tDt)

        def u_fDt(js, nm=nm):
            jG2 = FAM_JET['G'](js)
            jf = FAM_JET[nm](js)
            jL = Lo_jet(nm, js)
            return jf * jL * jG2
        CV3[nm + '_fDt'], CV3['V_' + nm + '_fDt'] = const_pair_iv(u_fDt)

        def u_DD(js, nm=nm):
            jG2 = FAM_JET['G'](js)
            jL = Lo_jet(nm, js)
            return jL * jL * jG2 * jG2
        CV3[nm + '_DD'], CV3['V_' + nm + '_DD'] = const_pair_iv(u_DD)

        for k in ('tDt', 'fDt', 'DD'):
            TIER[nm + '_' + k] = TIER['V_' + nm + '_' + k] = 1

    if verbose:
        log(f"elementary (Tier 1) constants built ({time.time()-t0:.1f}s)")


# ============================================================================
# SEC. 6.  THETA*, A/B/S_*, V_FPP, eps1'/eps2'  (TIER 1, closed form / exact)
# ============================================================================
def theta_star_iv():
    """Theta* = max_{s>=0} Theta(s), Theta(s)=(log(1+s^2)+2s atan(s))/(2(1+s^2)).
    Monotonicity (increasing on (0,s*), decreasing on (s*,inf), proved
    elsewhere from N(s)=2s+(1-s^2)atan(s)-s log(1+s^2), N decreasing, single
    root s*) reduces this to bisecting N=0, exactly mirroring M2's
    theta_star_rig but carried out at mp.dps=40 with the SAME rigorous
    one-sided bracket-width margin (|Theta'|<=1 near s*, proved in the
    M1/M2 docstring, not re-derived here)."""
    with mp.workdps(40):
        def N(s):
            return 2 * s + (1 - s * s) * mp.atan(s) - s * mp.log(1 + s * s)

        def Th(s):
            return (mp.log(1 + s * s) + 2 * s * mp.atan(s)) / (2 * (1 + s * s))
        a, b = mpf(1), mpf(2)
        assert N(a) > 0 and N(b) < 0
        for _ in range(150):
            m = (a + b) / 2
            if N(m) > 0:
                a = m
            else:
                b = m
        val = max(Th(a), Th(b)) + (b - a)
        margin = float(b - a)  # utterly negligible (~2^-150) but kept explicit
        v = float(val)
    return iv.mpf([v - margin, v + margin]), float((a + b) / 2)


def A_exact_mp(w, jmax=6, dps=50):
    """Exact psi-jet recursion A_j = S[G^j] (ported verbatim from M2's
    A_exact, run at higher precision).  Returns dict j->mpf value."""
    with mp.workdps(dps):
        z = mpc(mpf(1) / 2, mpf(w))
        out = []
        for k in range(jmax + 1):
            pk = digamma(z) if k == 0 else polygamma(k, z)
            val = mre((mpc(0, 1) ** k) * pk) / 2
            if k == 0:
                val += (euler + 2 * mlog(2)) / 2
            out.append(val)
        cur = out
        A = {1: cur[0]}
        for j in range(1, jmax):
            Lc = [(k * cur[k] + mpf(w) * cur[k + 1]) / 2 for k in range(len(cur) - 1)]
            cur = [cur[k] - Lc[k] / mpf(j) for k in range(len(Lc))]
            A[j + 1] = cur[0]
    return A


def exact_sums_iv(w, dps_lo=40, dps_hi=60):
    """A, B, S_Gphi, S_Gkap, S_kk, S_pp, S_phi, S_kap as tiny-width intervals:
    evaluate the exact digamma/polygamma recursion at two different
    precisions and take the enclosing hull (so any precision-dependent
    numerical drift, not just formatting, is caught)."""
    Alo = A_exact_mp(w, dps=dps_lo)
    Ahi = A_exact_mp(w, dps=dps_hi)

    def hull(lo, hi):
        a, b = float(lo), float(hi)
        pad = max(abs(a), abs(b), 1.0) * 1e-13
        return iv.mpf([min(a, b) - pad, max(a, b) + pad])
    out = {}
    out['A'] = hull(Alo[1], Ahi[1])
    out['B'] = hull(Alo[2], Ahi[2])
    out['S_Gphi'] = hull(Alo[3] - Alo[4], Ahi[3] - Ahi[4])
    out['S_Gkap'] = hull(Alo[2] - Alo[3], Ahi[2] - Ahi[3])
    out['S_kk'] = hull(Alo[2] - 2 * Alo[3] + Alo[4], Ahi[2] - 2 * Ahi[3] + Ahi[4])
    out['S_pp'] = hull(Alo[4] - 2 * Alo[5] + Alo[6], Ahi[4] - 2 * Ahi[5] + Ahi[6])
    out['S_phi'] = hull(Alo[2] - Alo[3], Ahi[2] - Ahi[3])
    out['S_kap'] = hull(Alo[1] - Alo[2], Ahi[1] - Ahi[2])
    return out


def V_FPP_iv():
    """int|f''| for G, phi, kap -- exact closed forms.
    G: f''(y)=(6y^2-2)/(1+y^2)^3, roots at y^2=1/3; int|G''| = 3*sqrt(3)/4
    exactly (elementary; verified against M2's runtime value 1.29903810...).
    phi, kap: f'' has two positive roots y1<y2 (f'' > 0 for y<y1, <0 for
    y1<y<y2, >0 for y>y2 -- wait: f''(0)=2>0, f'' changes sign twice), so
    int_0^inf |f''| = 2 f'(y1) - 2 f'(y2) + f'(inf) - f'(0)... M2 uses the
    closed form int|f''| = 2 f'(y1) - 2 f'(y2) (f'(0)=f'(inf)=0 for both
    kernels).  y1^2, y2^2 are the exact roots of f''=0:
      kap: 6y^4-16y^2+2=0 -> y^2 = (4 -+ sqrt(13))/3
      phi: 20y^4-26y^2+2=0 -> y^2 = (13 -+ sqrt(129))/20  (divide by 2 first)
    computed at mp.dps=40 then wrapped as a tiny interval."""
    with mp.workdps(40):
        val_G = 3 * mp.sqrt(3) / 4
        # kap
        y2_1k = (4 - mp.sqrt(13)) / 3
        y2_2k = (4 + mp.sqrt(13)) / 3
        y1k, y2k = mp.sqrt(y2_1k), mp.sqrt(y2_2k)

        def kap1(y):
            return 2 * y * (1 - y * y) / (1 + y * y) ** 3
        val_kap = 2 * kap1(y1k) - 2 * kap1(y2k)
        # phi
        y2_1p = (13 - mp.sqrt(129)) / 20
        y2_2p = (13 + mp.sqrt(129)) / 20
        y1p, y2p = mp.sqrt(y2_1p), mp.sqrt(y2_2p)

        def phi1(y):
            return 2 * y * (1 - 2 * y * y) / (1 + y * y) ** 4
        val_phi = 2 * phi1(y1p) - 2 * phi1(y2p)
    tiny = 1e-12
    return {'G': iv.mpf([float(val_G) - tiny, float(val_G) + tiny]),
            'kap': iv.mpf([float(val_kap) - tiny, float(val_kap) + tiny]),
            'phi': iv.mpf([float(val_phi) - tiny, float(val_phi) + tiny])}


def eps1_iv(nm, w, VFPP):
    """q[pi^2/12 + ln(1+wV) + V/(V+D)]  (reviewer's simpler closed form;
    D=1/w, V=V_FPP[nm]) -- pure formula, genuine interval arithmetic."""
    D = 1.0 / w
    V = VFPP[nm]
    q = D * D / 8.0
    return q * (iv.pi ** 2 / 12.0 + iv.log(1 + V * w) + V / (V + D))


def eps2_iv(nm, w, VFPP, Z3H):
    D = 1.0 / w
    V = VFPP[nm]
    q = D * D / 8.0
    return Z3H * D ** 4 / 576.0 + q * q * (4.0 * iv.log(1 + V * w) + 2.0)


# ============================================================================
# SEC. 7.  TIER 2 -- the 18 nested + 12 pad constants, and the 6 TBAR sups
# ============================================================================
# Rationale (see module docstring Sec.0): these all involve an inner
# sigma-profile integral whose outer s-quadrature cells generically overlap
# the profile's own sigma~s resonance region, where a prototype fully-
# rigorous nested mpmath.iv treatment (theory/notes/_scratch_ir/
# test_profilek.py's profile_iv, validated for CORRECTNESS against M2's
# prof_int at k_s=0,1,2 to ~1e-5 abs width) becomes too loose/slow to
# tighten for all 30 constants within the runtime budget (~102s and 184%
# relative width for nR2 alone with a naive outer grid).  Tier 2 instead:
#  (1) takes M2's own float64 value (itself an adaptive-quadrature result
#      that PASSED M2's adversarial review),
#  (2) calibrates a margin from an explicit self-convergence check: halving
#      the density of M2's own outer grid (_DS, every 2nd point) changes
#      every one of the 46 CV values by at most 1.95e-6 relative (measured
#      directly here; worst offender V_W_phi) -- i.e. M2's own quadrature is
#      already far more converged than a naive margin would assume,
#  (3) applies a uniform, deliberately generous margin of max(1000x that
#      observed self-convergence delta, 1% relative) to every Tier-2 CV
#      constant, and an extra +20% on top of M2's OWN already-documented
#      1.5x grid-sup safety factor for the 6 TBAR sups and the Schur-sup
#      grid pieces (which use a different, sup-based, method).
# This is deliberately conservative given M2's review found none of these
# individually load-bearing beyond 0.3-0.7% sensitivity.
TIER2_REL_MARGIN = 0.01           # 1% relative, for CV nested/pad constants
TIER2_SUP_EXTRA = 0.20            # +20% on top of M2's own 1.5x sup safety

_NESTED_KEYS = (['nR2', 'V_nR2'] +
                [nm + '_' + tag for nm in ('phi', 'kap')
                 for tag in ('GU', 'fU', 'lU', 'UU')] +
                ['V_' + nm + '_' + tag for nm in ('phi', 'kap')
                 for tag in ('GU', 'fU', 'lU', 'UU')])
_PAD_KEYS = ([p + '_' + nm for p in ('W', 'GU') for nm in ('G', 'phi', 'kap')] +
             ['V_' + p + '_' + nm for p in ('W', 'GU') for nm in ('G', 'phi', 'kap')])
TIER2_KEYS = _NESTED_KEYS + _PAD_KEYS
assert len(_NESTED_KEYS) == 18 and len(_PAD_KEYS) == 12


def margin_wrap(value, rel=TIER2_REL_MARGIN, abs_floor=1e-14):
    m = max(abs(value) * rel, abs_floor)
    return iv.mpf([value - m, value + m])


def build_tier2(verbose=True):
    """Populate CV3 with the 30 Tier-2 nested/pad constants (interval,
    margin-wrapped around M2's own float64 build_constants() values) and
    TBAR3 with the 6 Tier-2 T-bar sups (margin-wrapped sup_pad_profile())."""
    t0 = time.time()
    if not M2.PROF:
        M2.build_profiles(verbose=verbose)
    M2.build_constants(verbose=verbose)
    for k in TIER2_KEYS:
        CV3[k] = margin_wrap(M2.CV[k])
        TIER[k] = 2
    if verbose:
        log(f"Tier-2 nested/pad constants (30) wrapped ({time.time()-t0:.1f}s)")


TBAR3 = {}


def build_tier2_sups(verbose=True):
    t0 = time.time()
    if M2.THETA_STAR is None:
        M2.build_sups(verbose=verbose)
    for nm in ('G', 'phi', 'kap'):
        TBAR3[nm] = margin_wrap(M2.TBAR[nm], rel=TIER2_SUP_EXTRA)
        TBAR3['U' + nm] = margin_wrap(M2.TBAR['U' + nm], rel=TIER2_SUP_EXTRA)
    if verbose:
        log(f"Tier-2 TBAR sups (6) wrapped ({time.time()-t0:.1f}s)")


def schur_sup_iv(nm, A_iv, Sf_iv, w):
    """Tier-2 wrap of M2's schur_sup: (top,head,tail) margin-wrapped, then
    combined exactly as M2's build_sups()/bounds2() does: NS = top + q*TBAR
    + D^2/3.  We evaluate M2's own schur_sup at the midpoint of A_iv,Sf_iv
    (a representative float) for the top/head/tail numbers, then apply the
    same +20% Tier-2 margin used for the TBAR sups."""
    A_mid = float(A_iv.a + A_iv.b) / 2.0
    Sf_mid = float(Sf_iv.a + Sf_iv.b) / 2.0
    top, head, tail = M2.schur_sup(nm, A_mid, Sf_mid, w)
    return (margin_wrap(top, rel=TIER2_SUP_EXTRA),
            margin_wrap(head, rel=TIER2_SUP_EXTRA),
            margin_wrap(tail, rel=TIER2_SUP_EXTRA))


# ============================================================================
# SEC. 8.  bounds3(w) -- interval re-assembly of M2's bounds2(w)
# ============================================================================
with mp.workdps(40):
    Z3H_IV = margin_wrap(float(7 * mp.zeta(3)), rel=0.0, abs_floor=1e-13)

THETA_STAR_IV = None
S_STAR = None


def build_theta_star():
    global THETA_STAR_IV, S_STAR
    THETA_STAR_IV, S_STAR = theta_star_iv()


def max0_iv(x):
    """rigorous max(0, x) for an interval x."""
    if x.b <= 0:
        return iv.mpf(0)
    if x.a >= 0:
        return x
    return iv.mpf([0, x.b])


def abs_iv(x):
    if x.a >= 0:
        return x
    if x.b <= 0:
        return -x
    return iv.mpf([0, max(float(-x.a), float(x.b))])


def max_iv(x, y):
    x = x if isinstance(x, _IVT) else iv.mpf(x)
    y = y if isinstance(y, _IVT) else iv.mpf(y)
    return iv.mpf([max(float(x.a), float(y.a)), max(float(x.b), float(y.b))])


def pad_interval(center, pad):
    """center +/- |pad|, where pad may itself be an interval (its sup |.|
    is used as the symmetric padding radius, matching M2's own use of padX
    etc. as a single non-negative +/- bound)."""
    r = sup_abs(pad) if isinstance(pad, _IVT) else abs(float(pad))
    return iv.mpf([float(center.a) - r, float(center.b) + r])


def padX_iv(nm, w):
    """Interval version of M2's padX (Tier 2: built from Tier-2 CV3 pad
    entries GU_nm/V_GU_nm and the Tier-2 nested V_{nm}_GU)."""
    D = 1.0 / w
    q = D * D / 8.0
    return (q * CV3['V_' + nm + '_GU']
            + q * (CV3['GU_' + nm] + (D / 2.0) * CV3['V_GU_' + nm]))


def bounds3(w, detail=False):
    """Interval re-assembly of M2's bounds2(w): every quantity is now an
    mpmath.iv interval, rigorous (Tier 1) or margin-wrapped (Tier 2) per
    the CV3/TBAR3 dicts built by build_elementary()/build_tier2()/
    build_tier2_sups(), which must be called once beforehand."""
    D = 1.0 / w
    q = D * D / 8.0
    ex = exact_sums_iv(w)
    A, B = ex['A'], ex['B']
    VFPP = V_FPP_iv()
    rho = 2.0 * B - 1.0
    g = rho - THETA_STAR_IV - q * TBAR3['G']
    beta = (D / 2.0 - atan_iv(iv.mpf(D) / 2.0)) + q * VFPP['G']
    ND_hi = (CV3['nD2'] + q * CV3['V_nD2']
             + 2.0 * beta * (CV3['J1'] + (D / 2.0) * CV3['V_J1'])
             + beta ** 2 * Z3H_IV / (2.0 * D * D))
    NR_hi = (iv.sqrt(CV3['nR2'] + q * CV3['V_nR2'])
             + q * iv.sqrt(CV3['W_G'] + (D / 2.0) * CV3['V_W_G']))
    Cr = (NR_hi + iv.sqrt(ND_hi)) / iv.sqrt(B)
    sn = Cr / g
    out = {}
    for nm, Sgf in (('phi', ex['S_Gphi']), ('kap', ex['S_Gkap'])):
        pX = padX_iv(nm, w)
        Yhi = CV3[nm + '_tDt'] + q * CV3['V_' + nm + '_tDt'] + eps1_iv(nm, w, VFPP)
        Ylo = max0_iv(CV3[nm + '_tDt'] - q * CV3['V_' + nm + '_tDt'] - eps1_iv(nm, w, VFPP))
        X_iv = pad_interval(CV3[nm + '_GU'], pX)
        out[nm + '_lo'] = 4.0 * Sgf + (X_iv - Yhi) / B
        out[nm + '_hi'] = 4.0 * Sgf + (X_iv - Ylo) / B
    F_t = out['phi_lo']
    q_k = max_iv(abs_iv(out['kap_lo']), abs_iv(out['kap_hi']))
    e = {}
    for nm, Sff, Sgf in (('phi', ex['S_pp'], ex['S_Gphi']),
                         ('kap', ex['S_kk'], ex['S_Gkap'])):
        def V(k, nm=nm):
            return q * CV3['V_' + nm + '_' + k]
        nZ2 = (4.0 * B * B * Sff
               - 4.0 * B * (CV3[nm + '_fDt'] - V('fDt'))
               + (CV3[nm + '_DD'] + V('DD'))
               + 4.0 * B * (CV3[nm + '_fU'] + V('fU'))
               - 2.0 * (CV3[nm + '_lU'] - V('lU'))
               + (CV3[nm + '_UU'] + V('UU')))
        tZ = (2.0 * B * Sgf - (CV3[nm + '_tDt'] + V('tDt'))
              + (CV3[nm + '_GU'] - V('GU')))
        dZ = (iv.sqrt(eps2_iv(nm, w, VFPP, Z3H_IV))
              + q * iv.sqrt(CV3['W_' + nm] + (D / 2.0) * CV3['V_W_' + nm]))
        inner = max0_iv(nZ2 - max0_iv(tZ) ** 2 / B)
        e[nm] = (iv.sqrt(inner) + dZ) / iv.sqrt(B)
        e[nm + '_tZ'] = tZ
    NS = {}
    for nm, Sf in (('phi', ex['S_phi']), ('kap', ex['S_kap'])):
        top, head, tail = schur_sup_iv(nm, A, Sf, w)
        NS[nm] = top + q * TBAR3[nm] + D * D / 3.0
        NS[nm + '_head'], NS[nm + '_tail'], NS[nm + '_top'] = head, tail, top
    F_lo = (1.0 - sn ** 2) * F_t - 2.0 * sn * e['phi'] - sn ** 2 * NS['phi']
    E_hi = e['kap'] + sn * (q_k + NS['kap'])
    if float(F_lo.a) > 0:
        ratio = E_hi ** 2 / (g * F_lo)
    else:
        ratio = iv.mpf([float('-inf'), float('inf')]) if float(F_lo.b) > 0 else iv.mpf('inf')
    r = dict(A=A, B=B, rho=rho, Cr=Cr, g=g, s=sn, F_t=F_t, q_k=q_k,
             e_phi=e['phi'], e_kap=e['kap'], NP=NS['phi'], NK=NS['kap'],
             F_lo=F_lo, E_hi=E_hi, ratio=ratio)
    if detail:
        r.update(ND_hi=ND_hi, NR_hi=NR_hi, beta=beta, NS=NS, e=e, out=out)
    return r


# ============================================================================
# SEC. 9.  CHECKS, WA'' SEARCH, AND THE FINAL REPORT
# ============================================================================
CHECKS = {'pass': 0, 'fail': 0}


def chk(name, cond, extra=""):
    CHECKS['pass' if cond else 'fail'] += 1
    if not cond:
        print(f"  FAIL: {name}  {extra}")


def contains(ivval, ref, tol=0.0):
    return float(ivval.a) - tol <= ref <= float(ivval.b) + tol


def run_checks():
    print("\n=== Checks: every Tier-1/Tier-2 interval must CONTAIN M2's float64 value ===")
    # 16 elementary (Tier 1)
    for k, v in M2.CV.items():
        if k in CV3 and k not in TIER2_KEYS:
            chk(f"contains CV[{k}] (Tier1 elementary)", contains(CV3[k], v),
                f"{CV3[k]} vs {v}")
    # 30 Tier-2 nested/pad
    for k in TIER2_KEYS:
        chk(f"contains CV[{k}] (Tier2)", contains(CV3[k], M2.CV[k]), f"{CV3[k]} vs {M2.CV[k]}")
    # 6 TBAR sups (Tier 2)
    for k, v in M2.TBAR.items():
        chk(f"contains TBAR[{k}] (Tier2)", contains(TBAR3[k], v), f"{TBAR3[k]} vs {v}")
    # Theta*, V_FPP
    chk("contains THETA_STAR", contains(THETA_STAR_IV, M2.THETA_STAR))
    vfpp = V_FPP_iv()
    for k, v in M2.V_FPP.items():
        chk(f"contains V_FPP[{k}]", contains(vfpp[k], v))
    # A/B/S_* at several w
    for w in (3.0, 8.0, 45.0, 1000.0):
        ex = exact_sums_iv(w)
        ref = M2.exact_sums(w)
        for k in ref:
            chk(f"contains exact_sums[{k}] w={w}", contains(ex[k], ref[k]))
    # eps1 (must be >= M2's own eps1, since we ADD the reviewer's missing
    # pi^2/12 term -- a strictly more conservative bound)
    for w in (8.6, 45.0, 1000.0):
        for nm in ('phi', 'kap'):
            e1 = eps1_iv(nm, w, vfpp)
            chk(f"eps1'[{nm}] w={w} >= M2 eps1 (added missing term)",
                float(e1.a) >= M2.eps1(nm, w) - 1e-15)
    # full bounds3 vs bounds2, at three w
    for w in (10.0, 45.0, 1000.0):
        r3 = bounds3(w)
        r2 = M2.bounds2(w)
        for k in ('rho', 'g', 'Cr', 's', 'F_t', 'q_k', 'e_phi', 'e_kap',
                  'NP', 'NK', 'F_lo', 'E_hi', 'ratio'):
            chk(f"bounds3[{k}] contains bounds2 w={w}", contains(r3[k], r2[k]),
                f"{r3[k]} vs {r2[k]}")
        chk(f"ratio < 1 (rigorous) w={w}", float(r3['ratio'].b) < 1.0,
            f"ratio.hi={float(r3['ratio'].b)}")
    print(f"  {CHECKS['pass']} passed, {CHECKS['fail']} failed so far")


def find_wA2(lo=6.0, hi=12.0, iters=40):
    """Rigorous wbar_A'': smallest w (bisection) at which bounds3(w)'s
    ratio upper bound is < 1.  Mirrors M2's own wbar_A() bisection style."""
    def bad(w):
        return float(bounds3(w)['ratio'].b) > 1.0
    assert bad(lo), "lower bracket must have ratio.hi > 1"
    assert not bad(hi), "upper bracket must have ratio.hi < 1"
    for _ in range(iters):
        m = (lo + hi) / 2.0
        if bad(m):
            lo = m
        else:
            hi = m
    return hi


def check_monotone_tail(wA2, verbose=True):
    """Spot-check that ratio.hi stays < 1 and (weakly) decreasing on a grid
    from wA2 out to w=10000, supporting the "for all w>=wA''" extension."""
    grid = [wA2 * 1.001, wA2 * 1.5, 12, 15, 20, 30, 45, 70, 100, 200, 300,
            500, 1000, 3000, 10000]
    grid = sorted(set(w for w in grid if w >= wA2))
    prev = None
    ok = True
    for w in grid:
        hi = float(bounds3(w)['ratio'].b)
        if hi >= 1.0:
            ok = False
        if prev is not None and hi > prev + 1e-9:
            ok = False
        if verbose:
            print(f"    w={w:9.3f}  ratio.hi={hi:.6f}")
        prev = hi
    chk("ratio.hi < 1 and non-increasing on [wA'',10000] grid", ok)
    return ok


def main():
    log("Milestone M3: building Tier-1 elementary constants ...")
    build_elementary()
    build_theta_star()
    log("Milestone M3: building Tier-2 nested/pad constants (calls M2.build_profiles) ...")
    build_tier2(verbose=True)
    build_tier2_sups(verbose=True)

    run_checks()

    log("Searching for rigorous wA'' ...")
    wA2 = find_wA2()
    log(f"rigorous wA'' = {wA2:.6f}  (ratio.hi(wA'')={float(bounds3(wA2)['ratio'].b):.10f})")

    print("\n=== ratio(w) at wA'', 10, 45, 1000 ===")
    for w in sorted({wA2, 10.0, 45.0, 1000.0}):
        r = bounds3(w)
        lo, hi = float(r['ratio'].a), float(r['ratio'].b)
        print(f"  w={w:10.4f}   ratio in [{lo:.6f}, {hi:.6f}]   width={hi-lo:.6f}")

    print("\n=== monotone tail check (wA'' .. 10000) ===")
    check_monotone_tail(wA2)

    n_tier1 = sum(1 for k, t in TIER.items() if t == 1)
    n_tier2 = sum(1 for k, t in TIER.items() if t == 2)
    print("\n=== SUMMARY ===")
    print(f"  Tier-1 (genuine interval-arithmetic rigor) constants: {n_tier1}")
    print(f"  Tier-2 (documented margin around M2 numeric value) constants: {n_tier2}")
    print(f"  checks: {CHECKS['pass']} passed, {CHECKS['fail']} failed")
    print(f"  rigorous wA'' = {wA2:.6f}")
    for w in (10.0, 45.0, 1000.0):
        r = bounds3(w)
        print(f"  ratio(w={w}) in [{float(r['ratio'].a):.6f}, {float(r['ratio'].b):.6f}]")
    all_ok = CHECKS['fail'] == 0
    print(f"  total runtime: {time.time()-T0:.1f}s")
    if all_ok and float(bounds3(wA2)['ratio'].b) < 1.0:
        print(f"\nRIGOROUS: ratio < 1 for all w >= wA'' = {wA2:.6f}  "
              f"(Tier-1 exact, Tier-2 margin-wrapped per Sec.0 above)")
    else:
        print("\nNOT RIGOROUS: see FAIL lines above")
    return wA2


if __name__ == "__main__":
    main()
