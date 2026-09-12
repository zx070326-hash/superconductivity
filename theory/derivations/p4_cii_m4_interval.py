r"""Milestone M4 -- CLOSING THE INNER sigma-INTEGRAL IN CLOSED FORM, and the
resulting Tier-1 interval enclosures for Segment C-II (route R).

Builds on `p4_cii_m3_interval.py` (imported READ-ONLY): M3 rigorised 16
elementary constants + Theta* + A/B/S_* + V_FPP + eps1'/eps2' (Tier 1) but
left 36 quantities at Tier 2 (M2's float64 value wrapped in a documented
margin): 18 "Upsilon-nested" constants, 12 "pad-profile" constants, 6 T-bar
sups.  M3's stated obstruction was that all 36 carry an INNER sigma-integral
whose outer s-quadrature cells overlap the sigma~s resonance, making a naive
nested interval treatment both loose and slow.

M4 removes that obstruction for the 18 Upsilon-nested constants by
eliminating the inner integral ENTIRELY: Theta_f and Ups_f have elementary
closed forms.

================================================================================
SEC. 0.  THE CLOSED FORM  (the whole content of this milestone)
================================================================================
Write g_n(y) := 1/(1+y^2)^n, so that

        G = g_1,      kappa = g_1 - g_2,      phi_0 = g_2 - g_3.

M2's reduced two-point kernel and the two profiles are

  R_G(s,sg)  = G(s-sg) + G(s+sg) - 2 G(s) G(sg)
  R_f(s,sg)  = f(s-sg) + f(s+sg) - 2 f(s) G(sg) - 2 f(sg) G(s)   (f = kappa,phi_0)
  Theta_f(s) = (1/2) int_0^inf R_f(s,sg) dsg/sg
  Ups_f(s)   = (1/2) int_0^inf G(sg) R_f(s,sg) dsg/sg

Introduce the ONE-PARAMETER family  P(y;u) := 1/(u + y^2)  (so P(.;1) = G,
and g_n = (-1)^{n-1}/(n-1)! d_u^{n-1} P |_{u=1}), and the two u-blocks

  Chat(s;u) := int_0^inf [ P(sg-s;u) + P(sg+s;u) - 2 P(s;u) G(sg) ] dsg/sg
  Ehat(s;u) := int_0^inf G(sg)[ same bracket ] dsg/sg

Both are elementary.  With  K := u + s^2,  a := sqrt(u),  Del := (K-1)^2 + 4 s^2:

  Chat(s;u) = [ ln K + (2 s/a) atan(s/a) ] / K
  Ehat(s;u) = [ ln K + (2 s/a) atan(s/a) + 1 ] / K
              - (K-1)/Del * ln K  -  2 s (K+1)/(a Del) * atan(s/a)

(At u = 1 the first reduces to M2's own known Theta_G closed form, which is
how the family was normalised; both were checked against 30-digit tanh-sinh
quadrature of their defining integrals over a (s,u) grid -- max relative
error 4.3e-29, see theory/notes/_scratch_ir/m4_verify_family.py.)

Differentiating under the integral sign in u (legitimate: the bracket is
O(sg) at sg->0 and O(sg^-2) at infinity, locally uniformly in u near 1, for
every u-derivative order) and setting u = 1 gives the blocks

  C_n := (-1)^{n-1}/(n-1)! d_u^{n-1} Chat |_{u=1},   E_n likewise from Ehat,

and, using int_0^inf kappa/sg = 1/2, int phi_0/sg = 1/4, int G kappa/sg = 1/4,
int G phi_0/sg = 1/6 for the -2 f(sg) G(s) piece of R_f,

  Theta_G   = C_1/2                  Ups_G   = E_1/2
  Theta_kap = (C_1-C_2)/2 - G(s)/2   Ups_kap = (E_1-E_2)/2 - G(s)/4
  Theta_phi = (C_2-C_3)/2 - G(s)/4   Ups_phi = (E_2-E_3)/2 - G(s)/6

Every one of the six lands in the SAME three-term basis

        F(s) = A(s) + B(s) * log(1+s^2) + C(s) * atan(s),     A,B,C in Q(s),

which is closed under d/ds (since L' = 2s/(1+s^2) and T' = 1/(1+s^2) are
rational).  The exact integer coefficient tables are hard-coded in Sec. 1.
All six F, and their first and second s-derivatives, were checked against
35-digit quadrature of the ORIGINAL nested definitions (and of the
differentiated kernels d_s^k R_f) at 15 resp. 5 values of s -- max relative
error 3.5e-27 resp. 3.0e-31 (theory/notes/_scratch_ir/m4_verify_prof.py).

--------------------------------------------------------------------------------
SEC. 0b.  THE TWO REPRESENTATIONS  (why one formula is not enough)
--------------------------------------------------------------------------------
The direct form A + B L + C T is catastrophically ill-conditioned as s -> 0:
for Ups_kap, B ~ s^-2 and C ~ s^-3 while Ups_kap = O(s^2), so the three terms
cancel to ~10 digits at s = 1e-4 and ~19 digits at s = 1e-5 (measured).  In
INTERVAL arithmetic over a cell of nonzero width that cancellation is not
recovered, and the resulting h'' enclosure would destroy the Peano remainder.

Fix: substitute the two cancellation-free primitives

  m(s) := 1 - log(1+s^2)/s^2 = sum_{k>=1} (-1)^{k+1} s^{2k}/(k+1)
  d(s) := 1 - atan(s)/s      = sum_{k>=1} (-1)^{k+1} s^{2k}/(2k+1)

(the second is exactly M3's d_jet), i.e. L = s^2 (1-m), T = s (1-d), split
m = M_p + m_p and d = D_p + d_p into the degree-p truncation plus tail, and
absorb the polynomial parts into the rational block:

  F = R_p(s) - (B s^2) m_p(s) - (C s) d_p(s),
  R_p := A + B s^2 + C s - (B s^2) M_p - (C s) D_p    (exact, sympy `cancel`)

Since m_p, d_p = O(s^{2p+2}), the two correction blocks are multiplied by a
high power of s; writing them in "hatted" form

  Bs2h := (B s^2) s^{2p+2},  Csh := (C s) s^{2p+2},
  mhat_p := m_p / s^{2p+2},  dhat_p := d_p / s^{2p+2}

makes EVERY rational block regular at s = 0 (checked symbolically: with
p = 2 all 18 hatted/R denominators are nonvanishing at s=0, whereas the raw
C s of Ups_phi still carries a double pole).  Measured condition number
max|term|/|F| of this representation is <= 1.006 for all s <= 0.5 and all six
profiles, versus up to 8.2e18 for the direct form.  Conversely the direct
form has condition number <= 3.4 for all s >= 0.5.  So:

        s <= 0.5  ->  cancellation-free form (mhat/dhat as truncated
                      alternating series with explicit remainder)
        s >  0.5  ->  direct form A + B L + C T

exactly mirroring M3's own d_jet / Lo_jet dispatch at 0.5.

================================================================================
SEC. 0c.  WHAT M4 RIGORISES  (18 of M3's 36 Tier-2 quantities)
================================================================================
PROMOTED TO TIER 1 (genuine mpmath.iv end to end, no float64 input):
    nR2, V_nR2                                                    (2)
    {phi,kap}_{GU,fU,lU,UU}  and their V_ pairs                  (16)
These are exactly M3's `_NESTED_KEYS`.  With Ups_f now an elementary
function of s, each is a single-variable integral of an explicit closed-form
integrand.  M3's own Tier-1 Peano-midpoint scheme is kept for I[u]; for the
V constants its one-sided O(width) box bound is replaced by the sign-definite
exact rule of Sec. 3c, on an adaptive grid.  Result (46 checks, 0 failures,
189 s): all 18 intervals contain M2's independent float64 value STRICTLY,
and every one is tighter -- in the direction bounds3 uses -- than the Tier-2
margin it replaces, so the assembled threshold improves,
        wbar_A''  =  9.831910  (M3)   ->   9.711925  (M4).

STILL TIER 2 (unchanged from M3, and honestly so):
    W_{G,phi,kap}, GU_{G,phi,kap} and their V_ pairs              (12)
    TBAR[{G,phi,kap}], TBAR[U{G,phi,kap}]                          (6)
These are built from the PAD profiles
    T_f(s) = int_0^inf |d_sg^2 (R_f/2sg)| dsg,   U_f(s) = same with G(sg),
i.e. integrals of an ABSOLUTE VALUE of the sigma-derivative, which admits no
closed form; the u-family trick does not apply.  They remain at M3's
documented margin.  (They are the least load-bearing of the 36: each enters
only inside a factor q = D^2/8 = 1/(8 wbar^2); per M2's pad catalogue
(P4_einstein_extremal.md Sec. 1.G.3) each such term is 1e-5..4e-4 of the O(1)
quantity it pads at wbar = 45.)  Note that d_sg^2(R_f/2sg) is an elementary
RATIONAL function of sg, so its sign changes can be isolated by interval
arithmetic and the |.|-integral done piecewise with an elementary
antiderivative -- i.e. the 12 pad constants ARE closable by a different
route; not attempted here.

NOTE on M3's stated obstruction.  M3 recorded the nested constants as blocked
by a "fundamental" resonance between the outer s-cells and the inner
sigma-profile at sigma ~ s.  That diagnosis was correct for a NESTED interval
quadrature, but the obstruction is not fundamental: there is no inner
integral to enclose once Sec. 0's closed form is used.  What did turn out to
matter, and what M3's grid would have hit anyway, are two separate
conditioning problems -- the s -> 0 cancellation of A + B L + C T (Sec. 0b)
and the s -> 0 cancellation in forming u/(2s) (Sec. 3b) -- both of which are
removed by rearrangement rather than by refinement.

RESIDUAL, EXPLICITLY UNCLOSED (inherited from M3, quantified here):
M3's cell grid starts at s = 1e-4 and the head (0, 1e-4) is not enclosed by
any term.  For the M4 nested constants the integrand h is O(s) at worst
there, so the omission is < 1e-8 in M3.  M4 extends its own grid down to
s = 1e-8, reducing the residual to < 1e-16 for I and < 3e-15 for V2 (worst
case phi_GU/kap_GU: sup|h| = 7.6e-9, sup|h''| = 2.2e-7 on the head) --
immaterial against constants of size 1e-2..1, but it is an omission, not an
enclosure, and is reported as such.  The tail bound 10*B*sup|h([B,3B])|
beyond B = 1e6 (inherited from M3) rests on the power-law decay of the
integrands (measured exponents >= 5), not on a self-contained enclosure.
Apart from these two items every number entering CV4 is produced by
outward-rounded mpmath.iv arithmetic (adversarial review 2026-09-12:
SOUND; the two bit-level float leaks it flagged -- float Peano/V radii and
a float midpoint -- were closed afterwards, see _cell_record).

RUN:  python theory/derivations/p4_cii_m4_interval.py
"""
import os
import sys
import time

import numpy as np
from mpmath import mp, iv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_cii_m3_interval as M3  # noqa: E402  (read-only; NEVER modified)
import p4_cii_m2 as M2           # noqa: E402  (read-only reference)

from p4_cii_m3_interval import IJ3, make_seed, G_jet, FAM_JET, Lo_jet  # noqa: E402

iv.dps = 30
mp.dps = 30
T0 = time.time()


def log(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)


# ============================================================================
# SEC. 1.  EXACT COEFFICIENT TABLES  (ascending powers of s, integer coeffs)
# ============================================================================
# Generated once by theory/notes/_scratch_ir/m4_tab.py from the Sec.0 closed
# forms; every entry is an EXACT integer.  For each profile F = Theta_f
# (tag 'T') resp. Ups_f (tag 'U'):
#     direct form           F = A + B*log(1+s^2) + C*atan(s)
#     cancellation-free     F = R - Bs2h*mhat_p - Csh*dhat_p        (p = 2)
P_TRUNC = 2
SHIFT = 2 * P_TRUNC + 2          # m_p, d_p vanish to order s^SHIFT

PROF_TAB = {
    ('G', 'T'): {
        'A_n': (0,),
        'A_d': (1,),
        'B_n': (1,),
        'B_d': (2, 0, 2),
        'C_n': (0, 1),
        'C_d': (1, 0, 1),
        'R_n': (0, 0, 90, 0, -35, 0, 22),
        'R_d': (60, 0, 60),
        'Bs2h_n': (0, 0, 0, 0, 0, 0, 0, 0, 1),
        'Bs2h_d': (2, 0, 2),
        'Csh_n': (0, 0, 0, 0, 0, 0, 0, 0, 1),
        'Csh_d': (1, 0, 1),
    },
    ('G', 'U'): {
        'A_n': (1,),
        'A_d': (2, 0, 2),
        'B_n': (3,),
        'B_d': (8, 0, 10, 0, 2),
        'C_n': (-2, 0, 1),
        'C_d': (0, 4, 0, 5, 0, 1),
        'R_n': (0, 0, 220, 0, -89, 0, 42),
        'R_d': (240, 0, 300, 0, 60),
        'Bs2h_n': (0, 0, 0, 0, 0, 0, 0, 0, 3),
        'Bs2h_d': (8, 0, 10, 0, 2),
        'Csh_n': (0, 0, 0, 0, 0, 0, -2, 0, 1),
        'Csh_d': (4, 0, 5, 0, 1),
    },
    ('kap', 'T'): {
        'A_n': (0, 0, -1),
        'A_d': (1, 0, 2, 0, 1),
        'B_n': (0, 0, 1),
        'B_d': (2, 0, 4, 0, 2),
        'C_n': (0, -1, 0, 1),
        'C_d': (2, 0, 4, 0, 2),
        'R_n': (0, 0, -90, 0, 70, 0, -31, 0, 16),
        'R_d': (60, 0, 120, 0, 60),
        'Bs2h_n': (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
        'Bs2h_d': (2, 0, 4, 0, 2),
        'Csh_n': (0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1),
        'Csh_d': (2, 0, 4, 0, 2),
    },
    ('kap', 'U'): {
        'A_n': (6, 0, 1, 0, 1),
        'A_d': (16, 0, 36, 0, 24, 0, 4),
        'B_n': (-4, 0, -11, 0, 5, 0, 3),
        'B_d': (0, 0, 32, 0, 80, 0, 66, 0, 20, 0, 2),
        'C_n': (-8, 0, -36, 0, -9, 0, 1),
        'C_d': (0, 32, 0, 80, 0, 66, 0, 20, 0, 2),
        'R_n': (0, 0, -1120, 0, 392, 0, -176, 0, -59, 0, 36),
        'R_d': (960, 0, 2400, 0, 1980, 0, 600, 0, 60),
        'Bs2h_n': (0, 0, 0, 0, 0, 0, -4, 0, -11, 0, 5, 0, 3),
        'Bs2h_d': (32, 0, 80, 0, 66, 0, 20, 0, 2),
        'Csh_n': (0, 0, 0, 0, 0, 0, -8, 0, -36, 0, -9, 0, 1),
        'Csh_d': (32, 0, 80, 0, 66, 0, 20, 0, 2),
    },
    ('phi', 'T'): {
        'A_n': (0, 0, -13, 0, -1),
        'A_d': (8, 0, 24, 0, 24, 0, 8),
        'B_n': (0, 0, 1),
        'B_d': (2, 0, 6, 0, 6, 0, 2),
        'C_n': (0, -3, 0, 6, 0, 1),
        'C_d': (8, 0, 24, 0, 24, 0, 8),
        'R_n': (0, 0, -80, 0, 50, 0, -18, 0, 11, 0, 1),
        'R_d': (40, 0, 120, 0, 120, 0, 40),
        'Bs2h_n': (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
        'Bs2h_d': (2, 0, 6, 0, 6, 0, 2),
        'Csh_n': (0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 6, 0, 1),
        'Csh_d': (8, 0, 24, 0, 24, 0, 8),
    },
    ('phi', 'U'): {
        'A_n': (48, 0, 164, 0, -220, 0, -141, 0, -25, 0, -4),
        'A_d': (0, 0, 384, 0, 1344, 0, 1752, 0, 1032, 0, 264, 0, 24),
        'B_n': (4, 0, 13, 0, 78, 0, 52, 0, 10),
        'B_d': (0, 0, 128, 0, 480, 0, 696, 0, 490, 0, 174, 0, 30, 0, 2),
        'C_n': (-64, 0, -272, 0, -664, 0, -198, 0, 73, 0, 32, 0, 1),
        'C_d': (0, 0, 0, 512, 0, 1920, 0, 2784, 0, 1960, 0, 696, 0, 120, 0, 8),
        'R_n': (0, 0, -11712, 0, -16, 0, -72, 0, -84, 0, 794, 0, 291, 0, 3),
        'R_d': (7680, 0, 28800, 0, 41760, 0, 29400, 0, 10440, 0, 1800, 0, 120),
        'Bs2h_n': (0, 0, 0, 0, 0, 0, 4, 0, 13, 0, 78, 0, 52, 0, 10),
        'Bs2h_d': (128, 0, 480, 0, 696, 0, 490, 0, 174, 0, 30, 0, 2),
        'Csh_n': (0, 0, 0, 0, -64, 0, -272, 0, -664, 0, -198, 0, 73, 0, 32, 0, 1),
        'Csh_d': (512, 0, 1920, 0, 2784, 0, 1960, 0, 696, 0, 120, 0, 8),
    },
}


# ============================================================================
# SEC. 2.  INTERVAL POLYNOMIAL / RATIONAL / SERIES JETS
# ============================================================================
def _horner(coefs, x):
    acc = iv.mpf(0)
    for c in reversed(coefs):
        acc = acc * x + c
    return acc


def _derivs(coefs):
    """(coefs, coefs', coefs'') as ascending-power tuples."""
    return (tuple(coefs),
            tuple(k * coefs[k] for k in range(1, len(coefs))),
            tuple(k * (k - 1) * coefs[k] for k in range(2, len(coefs))))


_PJ = {}     # (table-key, field) -> (c0, c1, c2), built eagerly at import


def poly_jet(pj, js):
    """IJ3 jet of a polynomial at the SEED jet js, from a precomputed
    (c0, c1, c2) derivative triple.

    Value and both derivatives each get their own Horner scheme on the
    interval js.v -- tighter than accumulating IJ3 products, and it never
    squares a sign-crossing interval."""
    c0, c1, c2 = pj
    s = js.v
    return IJ3(_horner(c0, s),
               _horner(c1, s) if c1 else iv.mpf(0),
               _horner(c2, s) if c2 else iv.mpf(0))


def rat_jet(key, field, js):
    return poly_jet(_PJ[(key, field + '_n')], js) / poly_jet(_PJ[(key, field + '_d')], js)


def log1ps2_jet(js):
    """L(s) = log(1+s^2);  L' = 2s/(1+s^2);  L'' = 2(1-s^2)/(1+s^2)^2."""
    s = js.v
    s2 = s ** 2
    return IJ3(iv.log(1 + s2), 2 * s / (1 + s2), 2 * (1 - s2) / (1 + s2) ** 2)


_SER_K = 50          # terms j = 0..K of mhat/dhat


def _series_coeffs(which):
    """ascending coefficients (in s) of  sum_{j=0}^{K} (-1)^j s^{2j} / den(j),
    den(j) = j+P+2 for mhat, 2j+2P+3 for dhat  (see Sec.0b, P = P_TRUNC).
    Each 1/den is a correctly-rounded INTERVAL, so the truncation is enclosed
    exactly and the only remaining error is the tail bounded in _series_jet."""
    out = [iv.mpf(0)] * (2 * _SER_K + 1)
    for j in range(_SER_K + 1):
        den = (j + P_TRUNC + 2) if which == 'm' else (2 * j + 2 * P_TRUNC + 3)
        out[2 * j] = iv.mpf((-1) ** j) / den
    return tuple(out)


for _k, _tab in PROF_TAB.items():
    for _f, _c in _tab.items():
        _PJ[(_k, _f)] = _derivs(_c)
    # Shifted numerators, giving R/s, Bs2h/s, Csh/s with NO division by s.
    # Every one of these numerators has a zero constant term (in fact a
    # double zero), so the shift is exact polynomial division.  See Sec. 3b.
    for _f in ('R_n', 'Bs2h_n', 'Csh_n'):
        _c = _tab[_f]
        assert _c[0] == 0, f"{_k} {_f} is not divisible by s"
        _PJ[(_k, _f + '_s')] = _derivs(_c[1:])
_PJ[('series', 'm')] = _derivs(_series_coeffs('m'))
_PJ[('series', 'd')] = _derivs(_series_coeffs('d'))


def _series_jet(which, js):
    """IJ3 jet of mhat_p (which='m') resp. dhat_p (which='d'): the K-term
    truncation plus a rigorous remainder.  Valid for |s| <= 1/2, where the
    omitted terms are dominated termwise by C0 * b^{2j}, C0 = 1/(P+2) the
    largest |coefficient|, giving the geometric-derivative majorants below:
        sum_{j>J} b^{2j}            <= b^{2J+2}/(1-x)                (x = b^2)
        sum_{j>J} 2j b^{2j-1}       <= 2(J+1) b^{2J+1}/(1-x)^2   [code: 2(J+2)]
        sum_{j>J} 2j(2j-1) b^{2j-2} <= 4 sum_{j>J} j^2 x^{j-1}
                                    <= 8 (J+2)^2 b^{2J}/(1-x)^3
    (the last from sum_{i>=0}(m+i)^2 x^i <= 2(m+1)^2/(1-x)^3, m = J+1)."""
    s = js.v
    b = max(abs(float(s.a)), abs(float(s.b)))
    assert b <= 0.5 + 1e-12, "_series_jet only valid for |s| <= 0.5"
    acc = poly_jet(_PJ[('series', which)], js)
    # majorants evaluated in interval arithmetic; the upper endpoint is taken
    bb = iv.mpf(b)
    x = bb ** 2
    J = _SER_K
    C0 = iv.mpf(1) / (P_TRUNC + 2)
    rv = (C0 * bb ** (2 * J + 2) / (1 - x)).b
    ra = (C0 * (2 * J + 4) * bb ** (2 * J + 1) / (1 - x) ** 2).b
    rb = (C0 * 8 * (J + 2) ** 2 * bb ** (2 * J) / (1 - x) ** 3).b
    return IJ3(acc.v + iv.mpf([-rv, rv]),
               acc.a + iv.mpf([-ra, ra]),
               acc.b + iv.mpf([-rb, rb]))


# ============================================================================
# SEC. 3.  THE PROFILE JETS  Theta_f(s), Ups_f(s)  -- rigorous, closed form
# ============================================================================
SWITCH = 0.5          # cancellation-free below, direct above (see Sec. 0b)


def prof_jet(nm, tag, js):
    """Rigorous interval 2-jet of Theta_f (tag='T') resp. Ups_f (tag='U')."""
    key = (nm, tag)
    s = js.v
    b = max(abs(float(s.a)), abs(float(s.b)))
    if b <= SWITCH:
        return (rat_jet(key, 'R', js)
                - rat_jet(key, 'Bs2h', js) * _series_jet('m', js)
                - rat_jet(key, 'Csh', js) * _series_jet('d', js))
    return (rat_jet(key, 'A', js)
            + rat_jet(key, 'B', js) * log1ps2_jet(js)
            + rat_jet(key, 'C', js) * M3.atan_jet(js))


def ups_jet(nm, js):
    return prof_jet(nm, 'U', js)


# ---- SEC. 3b.  Ups_f(s)/s, built WITHOUT ever dividing by s ---------------
# Every integrand of Sec. 4 is h = u/(2s) with u = (bounded factors) * Ups_f,
# and Ups_f = O(s^2) as s -> 0.  Forming u first and dividing by 2s afterwards
# (M3's h_of_jet) is a disaster in interval arithmetic: the three terms of
#     h'' = u''/(2s) - u'/s^2 + u/s^3
# are each O(1/s) and cancel down to O(s), so the enclosure of h'' has width
# ~1/s while |h''| ~ s -- a measured 10^6 blow-up at s = 4e-5 for phi_GU,
# and one that bisection CANNOT fix (halving the cell halves its width but
# doubles the 1/s prefactor, leaving the product constant).
#
# The fix is to never divide by s at all.  In the cancellation-free
# representation F = R - Bs2h*mhat - Csh*dhat, all three numerators carry an
# explicit factor s (asserted at import above), so F/s is obtained by simply
# shifting those numerators down one power -- exact, and with the same
# conditioning as F itself.  Above the switch there is no small-s
# cancellation to worry about and plain division is used.
def prof_over_s_jet(nm, tag, js):
    """Rigorous interval 2-jet of Theta_f(s)/s resp. Ups_f(s)/s."""
    key = (nm, tag)
    s = js.v
    b = max(abs(float(s.a)), abs(float(s.b)))
    if b <= SWITCH:
        return (poly_jet(_PJ[(key, 'R_n_s')], js) / poly_jet(_PJ[(key, 'R_d')], js)
                - (poly_jet(_PJ[(key, 'Bs2h_n_s')], js)
                   / poly_jet(_PJ[(key, 'Bs2h_d')], js)) * _series_jet('m', js)
                - (poly_jet(_PJ[(key, 'Csh_n_s')], js)
                   / poly_jet(_PJ[(key, 'Csh_d')], js)) * _series_jet('d', js))
    return prof_jet(nm, tag, js) / js


def ups_over_s_jet(nm, js):
    return prof_over_s_jet(nm, 'U', js)


# ---- SEC. 3c.  cell grid and adaptive quadrature --------------------------
# M3's engine uses a FIXED graded grid and, per cell, the box bounds
#     I  :  h(midpoint)*width  +-  (width^3/24) sup_cell|h''|      (Peano)
#     V2 :  [0, width * sup_cell|h''|]
# The Peano term is O(width^3) and harmless, but the V2 box bound is only
# O(width) accurate AND one-sided, so on M3's grid it overestimates int|h''|
# by up to ~20x for the Ups-integrands (measured), which would more than
# cancel the benefit of the Tier-1 promotion.  M4 keeps the Peano treatment
# of I and replaces the V treatment by the sign-definite exact rule of
# _cell_record (int|h''| = |h'(b)-h'(a)| wherever h'' provably does not
# change sign on the cell), with adaptive bisection to chase the remaining
# sign-change cells and the Peano budget.
def build_cells_m4(smax_mid=300.0, big=1.0e6):
    """Base grid, extended down to 1e-8 so the unenclosed head shrinks from
    M3's <1e-8 to <1e-16 (Sec. 0c).  Adaptive refinement takes it from here."""
    edges = np.unique(np.concatenate([
        np.geomspace(1.0e-8, 0.5, 120),
        np.linspace(0.5, 8.0, 200),
        np.geomspace(8.0, smax_mid, 90),
        np.geomspace(smax_mid, big, 60),
    ]))
    return [(float(edges[i]), float(edges[i + 1])) for i in range(len(edges) - 1)
            if edges[i + 1] > edges[i]]


_CELLS_M4 = build_cells_m4()


def _abs_range(rng):
    """(inf, sup) of |x| over the interval rng -- both exact mpf (no float
    rounding: abs/min/max of the endpoints are exact)."""
    p, q = rng.a, rng.b
    lo = mp.mpf(0) if (p <= 0 <= q) else min(abs(p), abs(q))
    return lo, max(abs(p), abs(q))


def _definite(rng):
    return not (rng.a <= 0 <= rng.b)


def _cell_record(h_func, a, b, want_v1=False):
    """Per-cell data: midpoint value + Peano radius for I, and a two-sided
    enclosure of int_cell|h'| and int_cell|h''|.  Every stored quantity is
    an mpf obtained by outward-rounded interval arithmetic (the cell width
    w = b - a and the midpoint are formed as intervals, so neither float
    subtraction nor the float midpoint leaks into the result).

    KEY POINT for the V integrals.  The box bound width*[inf|.|, sup|.|] is
    only O(width) accurate, which is what makes M3's V constants so loose.
    But wherever the interval enclosure of h'' over the cell does not contain
    0, h'' is PROVABLY sign-definite there, and then

        int_cell |h''| = | int_cell h'' | = | h'(b) - h'(a) |

    exactly -- evaluated at the two endpoints as point intervals, i.e. to
    roundoff.  Same for int|h'| via |h(b)-h(a)|.  Since h', h'' are smooth
    with finitely many sign changes, all but a handful of cells take this
    exact branch, and the adaptive loop only has to chase the sign-change
    cells.  This is what makes the V constants tight enough to beat M3's
    Tier-2 margin at a tractable cell count."""
    ia, ib = iv.mpf(a), iv.mpf(b)
    wi = ib - ia
    jm = h_func(make_seed((ia + ib) / 2))
    hj = h_func(make_seed(iv.mpf([a, b])))
    _, d2hi = _abs_range(hj.b)
    _ends = {}

    def ends():
        if not _ends:
            _ends['a'] = h_func(make_seed(ia))
            _ends['b'] = h_func(make_seed(ib))
        return _ends['a'], _ends['b']

    if want_v1 and _definite(hj.a):
        ja, jb = ends()
        v1lo, v1hi = _abs_range(jb.v - ja.v)
    else:
        d1lo, d1hi = _abs_range(hj.a)
        v1lo, v1hi = (wi * d1lo).a, (wi * d1hi).b

    if _definite(hj.b):
        ja, jb = ends()
        v2lo, v2hi = _abs_range(jb.a - ja.a)
    else:
        d2lo, _ = _abs_range(hj.b)
        v2lo, v2hi = (wi * d2lo).a, (wi * d2hi).b

    return {'a': a, 'b': b, 'wi': wi, 'w': float(wi.b), 'val': jm.v,
            'peano': (wi ** 3 / 24 * d2hi).b,
            'v1lo': v1lo, 'v1hi': v1hi, 'v2lo': v2lo, 'v2hi': v2hi}


def enclose_m4(h_func, rel_I=2.0e-5, rel_V=2.0e-3, max_cells=12000,
               cells=None, want_v1=False):
    """Rigorous (I[u], V1[u], V2[u]) enclosure by adaptive bisection, where
    h_func is the INTEGRAND h = u/(2s) itself as an IJ3 jet function (built
    cancellation-free via ups_over_s_jet -- see Sec. 3b -- never as u/(2s)).

    rel_I: target for the summed Peano radius, relative to |I|.
    rel_V: target for the summed V box slack, relative to the V upper bound.
           Only the UPPER bound of V1/V2 is load-bearing downstream (every V
           enters bounds3 as a one-sided padding), so that is what the cells
           are spent on.
    want_v1: V1 is only needed by the const_pair1-style constants; the 18
           nested ones all use V2, so V1 is left at its (sound but loose) box
           bound by default rather than consuming cells.
    Tail beyond the last edge B: M3's bound 10*B*sup|h([B,3B])|, valid for
    any polynomial decay of order >= 2 (every integrand here decays >= s^-3).
    """
    recs = [_cell_record(h_func, a, b, want_v1) for (a, b) in (cells or _CELLS_M4)]

    # Refinement decisions below use plain floats; they only choose WHICH
    # cells to split and never enter the returned intervals.
    while len(recs) < max_cells:
        I_abs = abs(sum(float(r['val'].a) * r['w'] for r in recs))
        pe = sum(float(r['peano']) for r in recs)
        v2hi = sum(float(r['v2hi']) for r in recs)
        sl2 = sum(float(r['v2hi'] - r['v2lo']) for r in recs)
        tI = rel_I * max(I_abs, 1e-12)
        tV = rel_V * max(v2hi, 1e-12)
        need = (pe > tI) or (sl2 > tV)
        if want_v1:
            v1hi = sum(float(r['v1hi']) for r in recs)
            sl1 = sum(float(r['v1hi'] - r['v1lo']) for r in recs)
            tV1 = rel_V * max(v1hi, 1e-12)
            need = need or (sl1 > tV1)
        else:
            tV1 = None
        if not need:
            break

        def score(r):
            sc = float(r['peano']) / tI + float(r['v2hi'] - r['v2lo']) / tV
            if tV1 is not None:
                sc += float(r['v1hi'] - r['v1lo']) / tV1
            return sc
        n_split = max(1, len(recs) // 6)
        order = sorted(range(len(recs)), key=lambda i: -score(recs[i]))
        hot = set(order[:n_split])
        nxt = []
        for i, r in enumerate(recs):
            if i in hot:
                mid = 0.5 * (r['a'] + r['b'])
                if r['a'] < mid < r['b']:
                    nxt.append(_cell_record(h_func, r['a'], mid, want_v1))
                    nxt.append(_cell_record(h_func, mid, r['b'], want_v1))
                    continue
            nxt.append(r)
        if len(nxt) == len(recs):
            break
        recs = nxt

    I_acc = iv.mpf(0)
    V1_acc = iv.mpf(0)
    V2_acc = iv.mpf(0)
    for r in recs:
        I_acc = I_acc + r['val'] * r['wi'] + iv.mpf([-r['peano'], r['peano']])
        V1_acc = V1_acc + iv.mpf([r['v1lo'], r['v1hi']])
        V2_acc = V2_acc + iv.mpf([r['v2lo'], r['v2hi']])
    B = recs[-1]['b']
    h_tail = h_func(make_seed(iv.mpf([B, 3.0 * B])))
    tenB = 10 * iv.mpf(B)
    t0 = (tenB * _abs_range(h_tail.v)[1]).b
    t1 = (tenB * _abs_range(h_tail.a)[1]).b
    t2 = (tenB * _abs_range(h_tail.b)[1]).b
    I_acc = I_acc + iv.mpf([-t0, t0])
    V1_acc = V1_acc + iv.mpf([0, t1])
    V2_acc = V2_acc + iv.mpf([0, t2])
    return I_acc, V1_acc, V2_acc, len(recs)


NCELLS = {}


def const_pair_m4(name, h_func):
    """(I[u], V2[u]) as intervals -- mirrors M2's const_pair / M3's
    const_pair_iv, but on M4's adaptive grid and with h supplied directly."""
    I_acc, _, V2_acc, n = enclose_m4(h_func)
    NCELLS[name] = n
    return I_acc, V2_acc


# ============================================================================
# SEC. 4.  THE 18 UPSILON-NESTED CONSTANTS, NOW TIER 1
# ============================================================================
NESTED_KEYS = (['nR2', 'V_nR2'] +
               [nm + '_' + tag for nm in ('phi', 'kap')
                for tag in ('GU', 'fU', 'lU', 'UU')] +
               ['V_' + nm + '_' + tag for nm in ('phi', 'kap')
                for tag in ('GU', 'fU', 'lU', 'UU')])
assert len(NESTED_KEYS) == 18
assert set(NESTED_KEYS) == set(M3._NESTED_KEYS)

CV4 = {}
_HALF = 0.5


def build_nested_tier1(verbose=True):
    """Rebuild M2's 18 Upsilon-nested constants as genuine interval
    enclosures.  The integrands mirror M2.build_constants() exactly --
    h = u/(2s) with u as listed there -- except that M2's spline
    P(nm,'U',k,s) is replaced by the closed-form Ups jet, and the 1/(2s) is
    absorbed into ups_over_s_jet rather than applied afterwards (Sec. 3b)."""
    t0 = time.time()

    def h_nR2(js):                       # u = Ups_G^2
        return _HALF * ups_jet('G', js) * ups_over_s_jet('G', js)
    CV4['nR2'], CV4['V_nR2'] = const_pair_m4('nR2', h_nR2)
    if verbose:
        log(f"  nR2 done ({time.time()-t0:.1f}s, {NCELLS['nR2']} cells)")

    for nm in ('phi', 'kap'):
        def h_GU(js, nm=nm):             # u = G * Ups_f
            return _HALF * G_jet(js) * ups_over_s_jet(nm, js)
        CV4[nm + '_GU'], CV4['V_' + nm + '_GU'] = const_pair_m4(nm + '_GU', h_GU)

        def h_fU(js, nm=nm):             # u = f * Ups_f
            return _HALF * FAM_JET[nm](js) * ups_over_s_jet(nm, js)
        CV4[nm + '_fU'], CV4['V_' + nm + '_fU'] = const_pair_m4(nm + '_fU', h_fU)

        def h_lU(js, nm=nm):             # u = (Lam_f/s) * G * Ups_f
            return (_HALF * Lo_jet(nm, js) * G_jet(js)
                    * ups_over_s_jet(nm, js))
        CV4[nm + '_lU'], CV4['V_' + nm + '_lU'] = const_pair_m4(nm + '_lU', h_lU)

        def h_UU(js, nm=nm):             # u = Ups_f^2
            return _HALF * ups_jet(nm, js) * ups_over_s_jet(nm, js)
        CV4[nm + '_UU'], CV4['V_' + nm + '_UU'] = const_pair_m4(nm + '_UU', h_UU)
        if verbose:
            log(f"  {nm} GU/fU/lU/UU done ({time.time()-t0:.1f}s, cells "
                + "/".join(str(NCELLS[nm + '_' + k]) for k in ('GU', 'fU', 'lU', 'UU')) + ")")


# ============================================================================
# SEC. 5.  CHECKS
# ============================================================================
CHECKS = {'pass': 0, 'fail': 0}


def chk(name, cond, extra=""):
    CHECKS['pass' if cond else 'fail'] += 1
    if not cond:
        print(f"  FAIL: {name}  {extra}")


def check_profiles_vs_m2():
    """The closed-form profiles must reproduce M2's spline profiles (which
    come from an INDEPENDENT adaptive nested quadrature) at every k = 0,1,2."""
    print("\n=== closed-form Theta_f/Ups_f vs M2's nested-quadrature splines ===")
    pts = [0.02, 0.1, 0.25, 0.5, 0.9, 1.5009, 2.5, 5.0, 12.0, 40.0, 150.0]
    for nm in ('G', 'phi', 'kap'):
        for tag in ('T', 'U'):
            worst = [0.0, 0.0, 0.0]
            for x in pts:
                jet = prof_jet(nm, tag, make_seed(iv.mpf(x)))
                for k, comp in enumerate((jet.v, jet.a, jet.b)):
                    ref = float(M2.P(nm, tag, k, np.array([x]))[0])
                    mid = 0.5 * (float(comp.a) + float(comp.b))
                    worst[k] = max(worst[k], abs(mid - ref) /
                                   max(abs(ref), 1e-4))
            ok = worst[0] < 2e-5 and worst[1] < 2e-4 and worst[2] < 5e-3
            chk(f"{nm}/{tag} closed form matches M2 spline", ok,
                f"rel(k=0,1,2) = {worst[0]:.2e}, {worst[1]:.2e}, {worst[2]:.2e}")
            print(f"  {nm:4s} {tag}: rel vs M2 spline  k0={worst[0]:.2e} "
                  f"k1={worst[1]:.2e} k2={worst[2]:.2e}")


def check_constants_vs_m2():
    """Two things must hold for each of the 18:

    (1) CONTAINMENT: the Tier-1 interval must contain M2's float64 value (up
        to M2's own quadrature error -- M3 measured <= 1.95e-6 relative by
        grid halving; we allow 1e-4 relative + 1e-9 absolute).  This is the
        correctness test: two completely independent routes (M2's adaptive
        nested quadrature vs M4's closed form + interval quadrature) must
        agree.
    (2) NO REGRESSION: the M4 interval must not be WIDER, in the direction
        that bounds3 actually uses, than the Tier-2 margin it replaces.  For
        the eight I-constants both ends are load-bearing, so we compare full
        widths; for the ten V-constants only the upper end is (every V enters
        as a one-sided padding), so we compare upper bounds.
    """
    print("\n=== 18 Tier-1 nested constants: containment and tightness ===")
    print(f"  {'key':14s} {'M2 float64':>13s} {'M4 lo':>13s} {'M4 hi':>13s} "
          f"{'M3 Tier2 hi':>13s}  verdict")
    for k in NESTED_KEYS:
        ref = M2.CV[k]
        lo, hi = float(CV4[k].a), float(CV4[k].b)
        tol = 1e-4 * abs(ref) + 1e-9
        inside = (lo - tol) <= ref <= (hi + tol)
        chk(f"CV4[{k}] contains M2 value", inside,
            f"[{lo:.10g}, {hi:.10g}] vs {ref:.10g}")
        m3_hi = ref + abs(ref) * M3.TIER2_REL_MARGIN
        if k.startswith('V_'):
            better = hi <= m3_hi
            chk(f"{k}: M4 upper bound no looser than M3 Tier-2", better,
                f"{hi:.6g} vs {m3_hi:.6g}")
        else:
            better = (hi - lo) <= 2.0 * M3.TIER2_REL_MARGIN * abs(ref)
            chk(f"{k}: M4 width no looser than M3 Tier-2", better,
                f"{hi-lo:.3e} vs {2*M3.TIER2_REL_MARGIN*abs(ref):.3e}")
        verd = ("ok" if (inside and better) else
                ("CONTAINMENT FAIL" if not inside else "wider than M3"))
        print(f"  {k:14s} {ref:13.7g} {lo:13.7g} {hi:13.7g} {m3_hi:13.7g}  {verd}")


def check_cells():
    print("\n=== adaptive cell counts ===")
    print("  " + "  ".join(f"{k}={v}" for k, v in sorted(NCELLS.items())))


# ============================================================================
# SEC. 6.  RE-ASSEMBLY:  M3's bounds3() WITH THE 18 CONSTANTS UPGRADED
# ============================================================================
def install_into_m3():
    """Overwrite M3's CV3 entries for the 18 nested constants with the M4
    Tier-1 enclosures and flip their TIER flag.  Nothing else in M3 changes,
    so bounds3() is re-used verbatim."""
    for k in NESTED_KEYS:
        M3.CV3[k] = CV4[k]
        M3.TIER[k] = 1


def report_ratio():
    print("\n=== ratio(w) from M3's bounds3 with the 18 constants at Tier 1 ===")
    wA2 = M3.find_wA2()
    print(f"  rigorous wA'' = {wA2:.6f}")
    for w in sorted({wA2, 10.0, 45.0, 1000.0}):
        r = M3.bounds3(w)
        lo, hi = float(r['ratio'].a), float(r['ratio'].b)
        print(f"  w={w:10.4f}   ratio in [{lo:.6f}, {hi:.6f}]   width={hi-lo:.6f}")
        chk(f"ratio.hi < 1 at w={w}", hi < 1.0)
    return wA2


def main():
    log("M4: building M3's Tier-1 elementary constants and Theta* ...")
    M3.build_elementary(verbose=False)
    M3.build_theta_star()
    log("M4: building M2 profiles (needed only for the 12 pad + 6 sup "
        "Tier-2 quantities, and for cross-validation) ...")
    M3.build_tier2(verbose=False)
    M3.build_tier2_sups(verbose=False)

    check_profiles_vs_m2()

    log("M4: enclosing the 18 Upsilon-nested constants (Tier 1) ...")
    build_nested_tier1()

    check_constants_vs_m2()
    check_cells()

    wA2_m3 = M3.find_wA2()
    install_into_m3()
    wA2 = report_ratio()

    n1 = sum(1 for t in M3.TIER.values() if t == 1)
    n2 = sum(1 for t in M3.TIER.values() if t == 2)
    print("\n=== SUMMARY (M4) ===")
    print(f"  constants promoted Tier 2 -> Tier 1 by M4: {len(NESTED_KEYS)}")
    print(f"  of the 46 CV constants: Tier-1 {n1}, Tier-2 {n2} (the 12 pad "
          f"constants)")
    print(f"  plus {len(M3.TBAR3)} T-bar sups, still Tier 2 (tracked outside CV)")
    print(f"  wA''  M3 = {wA2_m3:.6f}   ->   M4 = {wA2:.6f}")
    print(f"  checks: {CHECKS['pass']} passed, {CHECKS['fail']} failed")
    print(f"  total runtime: {time.time()-T0:.1f}s")
    if CHECKS['fail'] == 0:
        print("\nOK: the 18 Upsilon-nested constants are now backed by genuine "
              "interval arithmetic on an elementary closed form.")
    else:
        print("\nNOT OK: see FAIL lines above")
    return wA2


if __name__ == "__main__":
    main()
