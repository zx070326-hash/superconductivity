r"""
P4 . b_1 CLOSED FORM  (Theorem W' sub-leading weak-coupling constant)
================================================================================

Companion numerical gate for `theory/derivations/P4_b1_closed_form.md`.

PART A -- closed form of

    b_1 = lim_{wbar->oo} [ 2<tau,(I-D)tau> + 2<tau,R tau> ]

for the KAY II Einstein operator  H(wbar) = 2 tau tau^T + calR - D  (identity
(M1) of `P4_einstein_extremal.md` 1.F.2), with

    D := 1/wbar,  s_n := (n+1/2)D,  w_n := 1/(2n+1),  p_n := sqrt(w_n),
    G(y) := 1/(1+y^2),  tau_n := p_n G(s_n),
    calR_{nm} := p_n p_m R(s_n,s_m),
    R(s,sg) := 2 s^2 sg^2 (5+s^2+sg^2)/[(1+(s-sg)^2)(1+(s+sg)^2)(1+s^2)(1+sg^2)],
    D_{nn} := w_n + ell_G(n),  ell_G(n) := (2/(2n+1)) sum_{k=1}^n G(k/wbar).

RESULT (this file):

    2<tau,(I-D)tau>  ->  J   := 3 pi^2/16 - 7/4   = 0.100550825204254...
    2<tau,calR tau>  ->  I_R := 11 pi^2/72 - 5/4  = 0.257856227944207...
    ------------------------------------------------------------------
    b_1 = J + I_R = 49 pi^2/144 - 3 = (7 pi/12)^2 - 3 = 0.3584070531484623...

and therefore

    T_c = (w_E/2pi) e^{c_0} e^{-1/lam} [1 + b_1 lam + O(lam^2)],
    c_0 = gamma + ln 4 - 3/2,   (1/2pi) e^{c_0} = (2/pi) e^{gamma-3/2} = 0.2529746...

PART B -- three-way cross-check of T_c(lam) for the Einstein spectrum, mu*=0,
constant DOS:
  (i)   `theory/solvers/eliashberg_iso.py` (independent linearised Eliashberg
        solver, dense non-symmetric eigenvalue problem),
  (ii)  the KAY operator lam_1(H^(N)(wbar)) via the FFT matvec of
        `p4_extremal_asymptotics.top_eig`,
  (iii) the asymptotic formula T_c^asym = (w_E/2pi) e^{c_0} e^{-1/lam}.

USAGE
    python p4_b1_closed_form.py            # ~4 min, all default checks
    python p4_b1_closed_form.py --hp       # + mpmath dps=60 / PSLQ block
    python p4_b1_closed_form.py --slow     # + solver at lam=0.15 (~10 min more)

Imports only (unchanged) `p4_interchange_uniformity.H_mat`,
`p4_extremal_asymptotics.{OpKernel,k_brk,top_eig}` and the solver package.
Counter `ck` records every assertion.
"""
import os
import sys
import time
from pathlib import Path

import numpy as np
import mpmath as mp
from scipy import integrate
from scipy.optimize import brentq
from scipy.sparse.linalg import LinearOperator, eigsh

HERE = Path(os.path.abspath(__file__)).parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "solvers"))
from p4_interchange_uniformity import H_mat                      # noqa: E402
from p4_extremal_asymptotics import OpKernel, k_brk, top_eig     # noqa: E402
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K        # noqa: E402

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HP = "--hp" in sys.argv
SLOW = "--slow" in sys.argv

ck = 0
T_START = time.time()

EG = 0.5772156649015329
C0 = EG + np.log(4.0) - 1.5                  # 0.46351001...
THETA_STAR = 0.6349287917881079              # certified (Lemma M)

J_CLOSED = 3 * np.pi ** 2 / 16 - 7 / 4       # 2<tau,(I-D)tau> limit
IR_CLOSED = 11 * np.pi ** 2 / 72 - 5 / 4     # 2<tau,R tau>   limit
B1_CLOSED = 49 * np.pi ** 2 / 144 - 3        # = J + I_R = (7 pi/12)^2 - 3


def _chk(name, ok, detail=""):
    global ck
    ck += 1
    if not ok:
        raise AssertionError(f"[FAIL #{ck}] {name}: {detail}")
    print(f"  ok  [{ck:3d}] {name}" + (f"   {detail}" if detail else ""))


def hdr(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# ===========================================================================
# closed forms used throughout
# ===========================================================================
def R_kernel(S, Sg):
    """R(s,sigma) of (M1), broadcasting."""
    num = 2 * S ** 2 * Sg ** 2 * (5 + S ** 2 + Sg ** 2)
    den = (1 + (S - Sg) ** 2) * (1 + (S + Sg) ** 2) * (1 + S ** 2) * (1 + Sg ** 2)
    return num / den


def Ups_G(s):
    """Upsilon_G(s) := (1/2) int_0^oo G(sg) R(s,sg) dsg/sg   [closed form, A.3]."""
    return ((3 * np.log(1 + s * s) + 2 * (s * s - 2) * np.arctan(s) / s + (4 + s * s))
            / (2 * (1 + s * s) * (4 + s * s)))


def Theta(s):
    """Theta(s) := (1/2) int_0^oo R(s,sg) dsg/sg   [1.F.3(c) of the extremal note]."""
    return (np.log(1 + s * s) + 2 * s * np.arctan(s)) / (2 * (1 + s * s))


def lattice_terms(N, w):
    """Finite-wbar lattice values of the two b_1 pieces, O(N log N) via FFT.

    Returns (two_tIDt, two_tRt, B, quad_H, ell_G) where
      two_tIDt = 2<tau,(I-D)tau>,  two_tRt = 2<tau,calR tau>,
      B = ||tau||^2,  quad_H = <tau, H^(N) tau>.

    The calR quadratic form is evaluated through (M1) backwards:
      2<tau,calR tau> = 2 [ a^T (T+H) a + sum_n w_n^2 G_n^2 ] - 4 B^2,
    a_n := w_n G(s_n), (T+H) the Toeplitz [[|n-m|]] (diagonal removed) plus
    Hankel [[n+m+1]] pair carried by OpKernel -- no O(N^2) matrix is formed.
    """
    n = np.arange(N)
    wn = 1.0 / (2 * n + 1.0)
    s = (n + 0.5) / w
    G = 1.0 / (1.0 + s ** 2)
    tau = np.sqrt(wn) * G
    B = float(tau @ tau)
    op = OpKernel(N, k_brk(N, w, 0))
    quad_H = op.quad(tau)
    ell_G = -op.diag                                     # OpKernel diag = -ell_G
    TH = quad_H + float((tau ** 2) @ ell_G)
    two_tRt = 2.0 * (TH + float(np.sum(wn ** 2 * G ** 2))) - 4.0 * B ** 2
    kap = np.arange(1, N) / w
    kap = kap ** 2 / (1.0 + kap ** 2)                    # 1 - G(k/wbar)
    one_minus_D = 2 * wn * np.concatenate(([0.0], np.cumsum(kap)))
    two_tIDt = 2.0 * float(np.sum(wn * G ** 2 * one_minus_D))
    return two_tIDt, two_tRt, B, quad_H, ell_G


# ===========================================================================
# A.1  the six elementary integral blocks of the derivation
# ===========================================================================
hdr("A.1  elementary integral blocks (the whole derivation is built from these)")
mp.mp.dps = 30
EPS = mp.mpf(10) ** -6

blocks = [
    # (name, numeric with cutoff eps, closed form, tolerance)
    ("X    = int_eps^oo ds/(s(1+s^2)^2)          = -ln eps - 1/2",
     mp.quad(lambda s: 1 / (s * (1 + s ** 2) ** 2), [EPS, 1, mp.inf]),
     -mp.log(EPS) - mp.mpf(1) / 2, mp.mpf(10) ** -11),
    ("A_2  = int_eps^oo arctan(s)/s^2 ds         = 1 - ln eps",
     mp.quad(lambda s: mp.atan(s) / s ** 2, [EPS, 1, mp.inf]),
     1 - mp.log(EPS), mp.mpf(10) ** -11),
    ("W(1) = int_0^oo arctan(s)/(1+s^2) ds       = pi^2/8",
     mp.quad(lambda s: mp.atan(s) / (1 + s ** 2), [0, 1, mp.inf]),
     mp.pi ** 2 / 8, mp.mpf(10) ** -22),
    ("V    = int_0^oo arctan(s)/(1+s^2)^2 ds     = pi^2/16 - 1/4",
     mp.quad(lambda s: mp.atan(s) / (1 + s ** 2) ** 2, [0, 1, mp.inf]),
     mp.pi ** 2 / 16 - mp.mpf(1) / 4, mp.mpf(10) ** -22),
    ("Psi  = int_eps^oo arctan(s)/(s^2(1+s^2)^2) = 5/4 - ln eps - 3 pi^2/16",
     mp.quad(lambda s: mp.atan(s) / (s ** 2 * (1 + s ** 2) ** 2), [EPS, 1, mp.inf]),
     mp.mpf(5) / 4 - mp.log(EPS) - 3 * mp.pi ** 2 / 16, mp.mpf(10) ** -11),
    ("Lam(1)= int_0^oo ln(1+s^2)/(s(1+s^2)) ds   = pi^2/12",
     mp.quad(lambda s: mp.log(1 + s ** 2) / (s * (1 + s ** 2)), [0, 1, mp.inf]),
     mp.pi ** 2 / 12, mp.mpf(10) ** -22),
    ("Lam2 = int_0^oo ln(1+s^2)/(s(1+s^2)^2) ds  = pi^2/12 - 1/2",
     mp.quad(lambda s: mp.log(1 + s ** 2) / (s * (1 + s ** 2) ** 2), [0, 1, mp.inf]),
     mp.pi ** 2 / 12 - mp.mpf(1) / 2, mp.mpf(10) ** -22),
    ("h(2) = int_0^2 ln x/(x-1) dx               = pi^2/4",
     mp.quad(lambda x: mp.log(x) / (x - 1), [0, 1, 2]),
     mp.pi ** 2 / 4, mp.mpf(10) ** -22),
]
for name, num, cls, tol in blocks:
    _chk(name, abs(num - cls) < tol, "err=%.2e" % float(abs(num - cls)))

# the ONE non-elementary-looking pair, killed by differentiation under the integral
W2 = mp.quad(lambda s: mp.atan(s) / (s ** 2 + 4), [0, 1, mp.inf])
Lam2a = mp.quad(lambda s: mp.log(1 + s ** 2) / (s * (s ** 2 + 4)), [0, 1, mp.inf])
_chk("W(2)+Lam(2) = pi^2/8   (the identity that removes all dilogarithms)",
     abs(W2 + Lam2a - mp.pi ** 2 / 8) < mp.mpf(10) ** -22,
     "W(2)=%s Lam(2)=%s err=%.2e" % (mp.nstr(W2, 12), mp.nstr(Lam2a, 12),
                                     float(abs(W2 + Lam2a - mp.pi ** 2 / 8))))

# ===========================================================================
# A.2  inner sigma-integral: Upsilon_G closed form, and Theta as a control
# ===========================================================================
hdr("A.2  inner integral  Upsilon_G(s) = (1/2) int_0^oo G(sg) R(s,sg) dsg/sg")
print("      claimed:  [3 ln(1+s^2) + 2(s^2-2) arctan(s)/s + (4+s^2)]"
      " / [2(1+s^2)(4+s^2)]")
for s0 in (0.05, 0.2, 0.5, 1.0, 1.500904, 2.5, 5.0, 20.0):
    fu = lambda x: R_kernel(s0, np.tan(x)) / (1 + np.tan(x) ** 2) / np.tan(x) / np.cos(x) ** 2
    vu, _ = integrate.quad(fu, 0, np.pi / 2, epsabs=1e-14, epsrel=1e-13, limit=400)
    ft = lambda x: R_kernel(s0, np.tan(x)) / np.tan(x) / np.cos(x) ** 2
    vt, _ = integrate.quad(ft, 0, np.pi / 2, epsabs=1e-14, epsrel=1e-13, limit=400)
    du, dt = abs(0.5 * vu - Ups_G(s0)), abs(0.5 * vt - Theta(s0))
    _chk("Upsilon_G & Theta at s=%.6f" % s0, du < 1e-11 and dt < 1e-11,
         "Ups=%.12f (err %.1e)  Theta=%.12f (err %.1e)" % (Ups_G(s0), du, Theta(s0), dt))
_chk("Theta(s*) reproduces the certified Theta* = 0.6349287917881079",
     abs(Theta(1.5009040) - THETA_STAR) < 1e-12,
     "Theta(1.5009040)=%.16f" % Theta(1.5009040))

# ===========================================================================
# A.3  I_R = 11 pi^2/72 - 5/4   (and the factor 1/2, not 1/4)
# ===========================================================================
hdr("A.3  outer integral  I_R = (1/2) int int G(s)G(sg)R(s,sg) ds dsg/(s sg)")
mp.mp.dps = 50 if not HP else 60


def Ups_mp(s):
    return ((3 * mp.log(1 + s ** 2) + 2 * (s ** 2 - 2) * mp.atan(s) / s + (4 + s ** 2))
            / (2 * (1 + s ** 2) * (4 + s ** 2)))


IR_mp = mp.quad(lambda s: Ups_mp(s) / (s * (1 + s ** 2)),
                [0, mp.mpf(1) / 2, 1, 2, 10, mp.inf])
IRc_mp = 11 * mp.pi ** 2 / 72 - mp.mpf(5) / 4
print("  I_R (mpmath dps=%d, 1-D reduced) = %s" % (mp.mp.dps, mp.nstr(IR_mp, 35)))
print("  11 pi^2/72 - 5/4                 = %s" % mp.nstr(IRc_mp, 35))
_chk("I_R = 11 pi^2/72 - 5/4", abs(IR_mp - IRc_mp) < mp.mpf(10) ** -30,
     "err=%.2e" % float(abs(IR_mp - IRc_mp)))

# genuinely 2-D control (no use of the Upsilon closed form at all)
f2 = lambda y, x: (R_kernel(np.tan(x), np.tan(y)) / (1 + np.tan(x) ** 2)
                   / (1 + np.tan(y) ** 2) / (np.tan(x) * np.tan(y))
                   / np.cos(x) ** 2 / np.cos(y) ** 2)
v2, e2 = integrate.dblquad(f2, 0, np.pi / 2, 0, np.pi / 2, epsabs=1e-12, epsrel=1e-11)
_chk("2-D dblquad control: (1/2)*I2 = I_R  (factor is 1/2, NOT 1/4)",
     abs(0.5 * v2 - IR_CLOSED) < 1e-10,
     "(1/2)I2=%.12f  (1/4)I2=%.12f  closed=%.12f" % (0.5 * v2, 0.25 * v2, IR_CLOSED))

# ===========================================================================
# A.4  (I-D) part:  exact discrete identity, then J = 3 pi^2/16 - 7/4
# ===========================================================================
hdr("A.4  (I-D) part:  1 - D_nn = 2 w_n sum_{k<=n} (1 - G(k/wbar))   [EXACT]")
for w0, N0 in ((7.0, 400), (50.0, 2000), (137.0, 5000)):
    n = np.arange(N0)
    wn = 1.0 / (2 * n + 1.0)
    Gk = 1.0 / (1.0 + (np.arange(1, N0) / w0) ** 2)
    D_nn = wn + np.concatenate(([0.0], np.cumsum(Gk))) * (2.0 / (2 * n + 1.0))
    alt = 2 * wn * np.concatenate(([0.0], np.cumsum(1.0 - Gk)))
    e = float(np.abs((1 - D_nn) - alt).max())
    _chk("exact identity at wbar=%g, N=%d" % (w0, N0), e < 1e-14, "max err=%.2e" % e)

_Jf = lambda s: (1.0 - np.arctan(s) / s) / (s * (1.0 + s * s) ** 2)
Jn1, _ = integrate.quad(_Jf, 0, 1, epsabs=1e-15, epsrel=1e-13, limit=400)
Jn2, _ = integrate.quad(_Jf, 1, np.inf, epsabs=1e-15, epsrel=1e-13, limit=400)
Jnum = Jn1 + Jn2
_chk("J = int_0^oo G^2 (1 - arctan s/s) ds/s = 3 pi^2/16 - 7/4",
     abs(Jnum - J_CLOSED) < 1e-12, "num=%.14f closed=%.14f" % (Jnum, J_CLOSED))

# Euler-Maclaurin coefficient: Delta_n = d(s_n) - (D^2/24) f'(s_n)/s_n + O(D^4),
# f(y) = 1-G(y),  f'(y)/y = 2 G(y)^2
print("\n  Euler-Maclaurin check  [Delta_n - d(s_n)] * 24 wbar^2 / (2 G(s_n)^2) -> -1")
for w0 in (40.0, 160.0, 640.0):
    N0 = int(20 * w0)
    n = np.arange(N0)
    wn = 1.0 / (2 * n + 1.0)
    s = (n + 0.5) / w0
    Gk = 1.0 / (1.0 + (np.arange(1, N0) / w0) ** 2)
    Delta = 2 * wn * np.concatenate(([0.0], np.cumsum(1.0 - Gk)))
    d = 1.0 - np.arctan(s) / s
    G = 1.0 / (1 + s ** 2)
    probe = (Delta - d) * 24 * w0 ** 2 / (2 * G ** 2)
    sl = slice(int(0.5 * w0), int(3 * w0))
    _chk("EM coefficient at wbar=%g" % w0, abs(probe[sl].mean() + 1.0) < 0.02,
         "mean over n in [wbar/2, 3wbar] = %.5f (target -1)" % probe[sl].mean())

# ===========================================================================
# A.5  lattice sums at finite wbar  vs  closed forms
# ===========================================================================
hdr("A.5  finite-wbar lattice sums vs the closed forms  (N = 200 wbar)")
N_check, w_check = 900, 11.0
tau_c = np.sqrt(1.0 / (2 * np.arange(N_check) + 1.0)) / (1 + ((np.arange(N_check) + .5) / w_check) ** 2)
d_c, r_c, B_c, qH_c, lG_c = lattice_terms(N_check, w_check)
H_ref = H_mat(N_check, w_check, 0)
_chk("(M1) quadratic form: <tau,H tau> = 2B^2 + <tau,R tau> - <tau,D tau>",
     abs(qH_c - (2 * B_c ** 2 + r_c / 2 - (B_c - d_c / 2))) < 1e-11,
     "lhs=%.12f rhs=%.12f" % (qH_c, 2 * B_c ** 2 + r_c / 2 - (B_c - d_c / 2)))
_chk("FFT quadratic form == dense H_mat quadratic form",
     abs(qH_c - float(tau_c @ (H_ref @ tau_c))) < 1e-11,
     "diff=%.2e" % abs(qH_c - float(tau_c @ (H_ref @ tau_c))))

print("\n%7s %9s %18s %14s %18s %14s" %
      ("wbar", "N", "2<t,(I-D)t>", "(err)*wbar^2", "2<t,Rt>", "(err)*wbar^2"))
rows_A = []
for w0 in (100, 200, 400):
    N0 = 200 * w0
    d1, r1, B, qH, _ = lattice_terms(N0, float(w0))
    rows_A.append((w0, N0, d1, r1))
    print("%7d %9d %18.12f %14.6f %18.12f %14.6f"
          % (w0, N0, d1, (d1 - J_CLOSED) * w0 ** 2, r1, (r1 - IR_CLOSED) * w0 ** 2))
for w0, N0, d1, r1 in rows_A:
    _chk("wbar=%d: |2<t,(I-D)t> - J| < 1.2 (ln wbar + 3)/(12 wbar^2)" % w0,
         abs(d1 - J_CLOSED) < 1.2 * (np.log(w0) + 3) / (12.0 * w0 ** 2),
         "err=%.3e  predicted ~ %.3e" % (abs(d1 - J_CLOSED), np.log(w0) / (12.0 * w0 ** 2)))
    _chk("wbar=%d: |2<t,Rt> - I_R| < 0.15/wbar^2" % w0,
         abs(r1 - IR_CLOSED) < 0.15 / w0 ** 2, "err=%.3e" % abs(r1 - IR_CLOSED))
sl = [(w, abs(d - J_CLOSED)) for w, _, d, _ in rows_A]
_chk("(I-D) error decays like wbar^-2 (ratio in [3.2,4.2] per doubling)",
     all(3.2 < sl[i][1] / sl[i + 1][1] < 4.2 for i in range(len(sl) - 1)),
     "ratios=" + ", ".join("%.3f" % (sl[i][1] / sl[i + 1][1]) for i in range(len(sl) - 1)))
sr = [(w, abs(r - IR_CLOSED)) for w, _, _, r in rows_A]
_chk("R-part error decays like wbar^-2 (ratio in [3.5,4.5] per doubling)",
     all(3.5 < sr[i][1] / sr[i + 1][1] < 4.5 for i in range(len(sr) - 1)),
     "ratios=" + ", ".join("%.3f" % (sr[i][1] / sr[i + 1][1]) for i in range(len(sr) - 1)))

# --- sign control on the secular-series step 4 of Theorem W'.
# Rayleigh quotient at tau (1.F.3(a)):  rho = 2B - 1 + (c_Delta + <tau,R tau>)/B,
# and the secular solve gives h = rho + O(1/ln^2 wbar).  The SIGN is PLUS; the
# 2026-09-13 draft of P4_weak_coupling_series.md step 4 prints delta = -[...]/(2B).
print("\n  sign control:  h - (2B-1)  vs  +(c_Delta + <tau,R tau>)/B   [must be PLUS]")
for w0 in (45.0, 200.0):
    N0 = int(40 * w0)
    d1, r1, B, qH, _ = lattice_terms(N0, w0)
    hh = top_eig(N0, k_brk(N0, w0, 0))[0]
    delta_num = hh - (2 * B - 1)
    delta_pred = (d1 / 2 + r1 / 2) / B
    print("    wbar=%6.1f:  h-(2B-1) = %+.6f   +(c_D+<t,Rt>)/B = %+.6f   "
          "(-) variant = %+.6f" % (w0, delta_num, delta_pred, -delta_pred))
    _chk("wbar=%g: h-(2B-1) matches the PLUS sign to O(1/ln^2 wbar)" % w0,
         delta_num > 0 and abs(delta_num - delta_pred) < 2.0 / np.log(w0) ** 2,
         "|diff| = %.4f, O(1/ln^2 wbar) budget = %.4f"
         % (abs(delta_num - delta_pred), 2.0 / np.log(w0) ** 2))

# ===========================================================================
# A.6  LEMMA B1: the limit EXISTS, with provable remainder orders
# ===========================================================================
hdr("A.6  Lemma B1 -- explicit remainder constants (existence of the limit)")
print("""  Lemma B1.  With D := 1/wbar, for every wbar >= 1

   (a) 2<tau,(I-D)tau> = J  -  (D^2/6) sum_n w_n G(s_n)^4  +  rho_a(wbar),
       |rho_a| <= (D^2/8) V2[phi] + C_a D^3,     phi(s) := G(s)^2 d(s)/(2s),
       d(s) := 1 - arctan(s)/s,   J = 3 pi^2/16 - 7/4.
       Since 0 <= sum_n w_n G(s_n)^4 <= (ln wbar + gamma + ln4 + 1)/2, the whole
       remainder is O(ln wbar / wbar^2).   NO LATTICE CONSTANT: the w_n inside
       D_nn = w_n + ell_G(n) cancels IDENTICALLY (see A.4), so C_lattice = 0.

   (b) 2<tau,calR tau> = I_R + rho_b(wbar),  |rho_b| <= (D^2/4) V2x[psi],
       psi(s,sg) := G(s)G(sg) R(s,sg)/(4 s sg),  I_R = 11 pi^2/72 - 5/4,
       V2x[psi] := int int |d^2 psi/ds^2| ds dsg  (finite; computed below).

   Hence lim 2<tau,(I+calR-D)tau> EXISTS and equals b_1 = 49 pi^2/144 - 3.

   Ingredients, all elementary:
     (i)  w_n = D/(2 s_n) and s_n = midpoint of [nD,(n+1)D]  (1.F.1 (i),(iii));
     (ii) 4th-order Euler-Maclaurin for sum_{k<=n} f(kD), f := 1-G, f(0)=f'(0)=
          f'''(0)=0, remainder controlled by V4[f] := int_0^oo |f''''| < oo;
          after the half-shift nD = s_n - D/2 the D^1 AND D^2 terms collapse to
          the single term -(D/24) f'(s_n), giving
          Delta_n = d(s_n) - (D^2/12) G(s_n)^2 + O(D^4/s_n);
    (iii) composite midpoint rule with the PROVABLE Peano constant D^2/8 paired
          with int|u''| (NOT D^2/24 -- that pairs with per-cell sup, see the
          M2 correction in P4_einstein_extremal.md 1.G).""")

import sympy as _sp  # noqa: E402

_s, _sg = _sp.symbols('s sigma', positive=True)
_G = 1 / (1 + _s ** 2)
_d = 1 - _sp.atan(_s) / _s
_phi = _G ** 2 * _d / (2 * _s)
_phi2 = _sp.lambdify(_s, _sp.diff(_phi, _s, 2), 'numpy')
V2_phi, _e = integrate.quad(lambda x: abs(_phi2(np.tan(x))) / np.cos(x) ** 2,
                            1e-12, np.pi / 2 - 1e-12, epsabs=1e-12, epsrel=1e-10,
                            limit=600)
_chk("V2[phi] = int_0^oo |phi''| ds is finite", np.isfinite(V2_phi) and V2_phi > 0,
     "V2[phi] = %.8f   (=> (D^2/8) V2[phi] = %.4f/wbar^2)" % (V2_phi, V2_phi / 8))

_Rsym = (2 * _s ** 2 * _sg ** 2 * (5 + _s ** 2 + _sg ** 2)
         / ((1 + (_s - _sg) ** 2) * (1 + (_s + _sg) ** 2) * (1 + _s ** 2) * (1 + _sg ** 2)))
_psi = _Rsym / ((1 + _s ** 2) * (1 + _sg ** 2) * 4 * _s * _sg)
_psi_ss = _sp.lambdify((_s, _sg), _sp.diff(_psi, _s, 2), 'numpy')


def V2x_psi(M):
    """int int |d^2 psi/ds^2| on (0,oo)^2 by composite midpoint on the tan-map."""
    x = (np.arange(M) + 0.5) * (np.pi / 2) / M
    hx = (np.pi / 2) / M
    t, jac = np.tan(x), 1.0 / np.cos(x) ** 2
    tot = 0.0
    for i0 in range(0, M, 256):
        S = t[i0:i0 + 256][:, None]
        tot += float(np.sum(np.abs(_psi_ss(S, t[None, :]))
                            * jac[i0:i0 + 256][:, None] * jac[None, :]))
    return tot * hx * hx


v_lo, v_hi = V2x_psi(1200), V2x_psi(2400)
_chk("V2x[psi] = int int |d^2psi/ds^2| is finite (two resolutions agree to 1%)",
     np.isfinite(v_hi) and abs(v_hi / v_lo - 1) < 0.01,
     "M=1200: %.6f   M=2400: %.6f  => (D^2/4)V2x = %.3f/wbar^2" % (v_lo, v_hi, v_hi / 4))

print("\n  Lemma B1 bounds evaluated against the measured finite-wbar errors:")
print("  %7s %16s %16s %14s %14s" %
      ("wbar", "meas err (I-D)", "after -D^2/6 term", "bound a", "meas/bound"))
for w0, N0, d1, r1 in rows_A:
    n = np.arange(N0)
    G4 = (1.0 / (1.0 + ((n + 0.5) / w0) ** 2)) ** 4
    S4 = float(np.sum(G4 / (2 * n + 1.0)))
    main = -(1.0 / (6.0 * w0 ** 2)) * S4
    resid = (d1 - J_CLOSED) - main
    bnd = V2_phi / (8.0 * w0 ** 2)
    print("  %7d %16.3e %16.3e %14.3e %14.4f" % (w0, d1 - J_CLOSED, resid, bnd, abs(resid) / bnd))
    _chk("Lemma B1(a) bound holds at wbar=%d" % w0, abs(resid) < bnd,
         "|residual|=%.3e <= (D^2/8)V2[phi]=%.3e" % (abs(resid), bnd))
print("  %7s %16s %30s %14s" % ("wbar", "meas err (R)", "", "bound b"))
for w0, N0, d1, r1 in rows_A:
    bnd = v_hi / (4.0 * w0 ** 2)
    print("  %7d %16.3e %30s %14.3e" % (w0, r1 - IR_CLOSED, "", bnd))
    _chk("Lemma B1(b) bound holds at wbar=%d" % w0, abs(r1 - IR_CLOSED) < bnd,
         "|err|=%.3e <= (D^2/4)V2x[psi]=%.3e" % (abs(r1 - IR_CLOSED), bnd))

# ===========================================================================
# A.7  b_1 and its identification
# ===========================================================================
hdr("A.7  b_1 = 49 pi^2/144 - 3 = (7 pi/12)^2 - 3")
mp.mp.dps = 60 if HP else 50
b1_mp = IR_mp + 3 * mp.pi ** 2 / 16 - mp.mpf(7) / 4
b1c_mp = 49 * mp.pi ** 2 / 144 - 3
print("  2<t,(I-D)t> -> 3 pi^2/16 - 7/4  = %s" % mp.nstr(3 * mp.pi ** 2 / 16 - mp.mpf(7) / 4, 30))
print("  2<t,R  t>   -> 11 pi^2/72 - 5/4 = %s" % mp.nstr(IRc_mp, 30))
print("  b_1         -> 49 pi^2/144 - 3  = %s" % mp.nstr(b1c_mp, 30))
_chk("b_1 (integrals) == 49 pi^2/144 - 3", abs(b1_mp - b1c_mp) < mp.mpf(10) ** -30,
     "err=%.2e" % float(abs(b1_mp - b1c_mp)))
ident = mp.identify(b1_mp, ['pi**2'])
print("  mp.identify(b_1, ['pi**2']) = %s" % ident)
_chk("mp.identify returns the rational multiple of pi^2", ident is not None and "49" in str(ident),
     str(ident))
rel = mp.pslq([b1_mp, mp.mpf(1), mp.pi ** 2, mp.pi, mp.log(2), mp.log(2) ** 2,
               mp.catalan, mp.zeta(3)], maxcoeff=10 ** 6, maxsteps=10 ** 6,
              tol=mp.mpf(10) ** -(mp.mp.dps - 8))
print("  PSLQ [b1, 1, pi^2, pi, ln2, ln^2 2, G_Catalan, zeta(3)] = %s" % (rel,))
_chk("PSLQ finds 144 b1 + 432 - 49 pi^2 = 0 and nothing else",
     rel is not None and abs(rel[0]) == 144 and abs(rel[2]) == 49
     and all(c == 0 for c in rel[3:]), str(rel))

# corollary: a two-term dilogarithm identity falls out of W(2)+Lam(2)=pi^2/8
dl = mp.polylog(2, -3) + 2 * mp.polylog(2, -2) + mp.pi ** 2 / 3 + 2 * mp.log(2) * mp.log(3)
_chk("corollary  Li_2(-3) + 2 Li_2(-2) = -pi^2/3 - 2 ln2 ln3",
     abs(dl) < mp.mpf(10) ** -30, "residual=%.2e" % float(abs(dl)))

print("\n  physical translation:")
pref = np.exp(C0) / (2 * np.pi)
print("    T_c = (w_E/2pi) e^{c0} e^{-1/lam} [1 + b1 lam + O(lam^2)]")
print("    (1/2pi) e^{c0} = (2/pi) e^{gamma-3/2} = %.10f" % pref)
print("    b1 = 49 pi^2/144 - 3 = %.13f" % B1_CLOSED)
_chk("prefactor (2/pi) e^{gamma-3/2} = 0.2529997 (the 0.2530 of the literature)",
     abs(pref - 2 / np.pi * np.exp(EG - 1.5)) < 1e-15 and abs(pref - 0.2529997) < 1e-6,
     "%.10f" % pref)

# ===========================================================================
# B.  solver / KAY operator / asymptotic formula
# ===========================================================================
hdr("B.  T_c cross-check: eliashberg_iso solver vs KAY operator vs asymptotics")

W0_MEV = 100.0          # w_E; only sets the unit, every output is T_c/w_E
LAMS = [0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
SOLVER_LAMS = [0.20, 0.25, 0.30, 0.40, 0.50]


def h_N(w, N):
    return top_eig(N, k_brk(N, w, 0))[0]


def wbar_of_lambda(lam, mult=40, nmin=2000):
    """Solve lam * lam_1(H^(N)(wbar)) = 1 for wbar. Two-stage: cheap bracket, then
    secant refinement at the production truncation N = max(nmin, mult*wbar)."""
    target = 1.0 / lam
    f_cheap = lambda x: h_N(np.exp(x), max(1000, int(8 * np.exp(x)))) - target
    x = brentq(f_cheap, target - C0 - 1.5, target - C0 + 1.5, xtol=1e-9)
    x0, x1 = x, x + 1e-3
    f0 = h_N(np.exp(x0), max(nmin, int(mult * np.exp(x0)))) - target
    for _ in range(6):
        f1 = h_N(np.exp(x1), max(nmin, int(mult * np.exp(x1)))) - target
        if abs(f1) < 1e-12 or abs(f1 - f0) < 1e-15:
            break
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x0, f0, x1 = x1, f1, x2
    return float(np.exp(x1)), max(nmin, int(mult * np.exp(x1)))


print("\n  B.1  mu*=0: w_cut is inert except through the automatic N estimate")
sp_t = Spectrum.einstein(lam=0.4, w0=W0_MEV)
e_a = Solver(sp_t, mu_star=0.0, w_cut=3 * W0_MEV, n_matsubara=300)._max_eigen_linearized(3.0, 300)
e_b = Solver(sp_t, mu_star=0.0, w_cut=50 * W0_MEV, n_matsubara=300)._max_eigen_linearized(3.0, 300)
_chk("two w_cut values, same n_matsubara -> bit-identical linearised eigenvalue",
     e_a == e_b, "%.16f vs %.16f" % (e_a, e_b))

print("\n  B.2  identified convention difference (solver Z-sum truncation)")
print("       solver: Z_n = 1 + lam w_n [1 + 2 sum_{k<=n}[[k]] - sum_{j=N-n}^{N+n}[[j]]]")
print("       KAY   : D_nn uses the exact infinite tail, i.e. the last sum absent.")
print("       => the solver eigenproblem is  lam (H^(N) + diag(eps)) x = x  with")
print("          eps_n = w_n sum_{j=N-n}^{N+n} [[j]],  <v,diag(eps)v> = wbar^2/N^2 + ...")


def eps_diag(N, w):
    n = np.arange(N)
    f = k_brk(N + 1, w, 0)[:2 * N + 2]
    cs = np.concatenate(([0.0], np.cumsum(f[1:2 * N + 1])))
    lo = np.maximum(N - n, 1)
    return (cs[N + n] - cs[lo - 1]) / (2 * n + 1.0)


for w0 in (15.6, 40.0):
    for mult in (20, 40):
        N0 = int(mult * w0)
        h0 = np.linalg.eigvalsh(H_mat(N0, w0, 0))[-1]
        h1 = np.linalg.eigvalsh(H_mat(N0, w0, 0) + np.diag(eps_diag(N0, w0)))[-1]
        c = (h1 - h0) * N0 ** 2 / w0 ** 2
        _chk("delta h * N^2/wbar^2 = 1 at wbar=%g, N=%dwbar" % (w0, mult),
             abs(c - 1.0) < 0.01, "coefficient = %.5f" % c)

print("\n  B.3  KAY operator route  (lam * lam_1(H^(N)(wbar)) = 1, N = 40 wbar)")
kay = {}
for lam in LAMS:
    t0 = time.time()
    w, N = wbar_of_lambda(lam)
    hh = h_N(w, N)
    tc_over_wE = 1.0 / (2 * np.pi * w)
    tc_asym = np.exp(C0) * np.exp(-1.0 / lam) / (2 * np.pi)
    r = (tc_over_wE / tc_asym - 1.0) / lam
    kay[lam] = dict(wbar=w, N=N, h=hh, tc=tc_over_wE, tc_asym=tc_asym, r=r,
                    resid=(hh - np.log(w) - C0) * np.log(w))
    print("    lam=%.2f  wbar=%12.5f  N=%7d  h=%.10f  Tc/wE=%.6e  r=%.6f  (%.1fs)"
          % (lam, w, N, hh, tc_over_wE, r, time.time() - t0))
    _chk("KAY root solves lam*h = 1 at lam=%.2f" % lam, abs(lam * hh - 1.0) < 1e-10,
         "lam*h - 1 = %.2e" % (lam * hh - 1.0))

# truncation control: redo lam=0.30 at N = 160 wbar
w30 = kay[0.30]["wbar"]
h30_fine = h_N(w30, int(160 * w30))
_chk("N-truncation of the KAY route is below 1e-6 in h (N=40wbar vs 160wbar)",
     abs(h30_fine - kay[0.30]["h"]) < 1e-6,
     "|dh| = %.2e" % abs(h30_fine - kay[0.30]["h"]))

print("\n  B.4  solver route  (eliashberg_iso, mu*=0, constant DOS)")
print("       (i)  MATCHED CONVENTION: at the SAME N, solver T_c must equal the root of")
print("            lam * lam_1(H^(N) + diag(eps_N)) = 1  -- this is an identity, not a fit.")
print("       (ii) Richardson in 1/N^2 removes the wbar^2/N^2 shift -> N->oo KAY value.")


def h_N_eps(w, N, v0=None):
    """lam_1 of the SOLVER's finite-N operator H^(N)(wbar) + diag(eps_N(wbar))."""
    op = OpKernel(N, k_brk(N, w, 0))
    ep = eps_diag(N, w)
    L = LinearOperator((N, N), matvec=lambda v: op.matvec(v) + ep * v, dtype=np.float64)
    seed = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
    vals = eigsh(L, k=1, which="LA", tol=0, ncv=min(48, N - 1), v0=seed)[0]
    return float(vals[0])


sol_res = {}
for lam in sorted(SOLVER_LAMS):
    w = kay[lam]["wbar"]
    N1 = max(256, int(10 * w))
    N2 = 2 * N1
    tc_asym_K = np.exp(C0) * np.exp(-1.0 / lam) * W0_MEV / (2 * np.pi) / KB_MEV_PER_K
    br = (0.70 * tc_asym_K, 1.60 * tc_asym_K)
    sp = Spectrum.einstein(lam=lam, w0=W0_MEV)
    t0 = time.time()
    tc1 = Solver(sp, mu_star=0.0, w_cut=5 * W0_MEV, n_matsubara=N1).tc(t_bracket_k=br, xtol=1e-11)
    tc2 = Solver(sp, mu_star=0.0, w_cut=5 * W0_MEV, n_matsubara=N2).tc(t_bracket_k=br, xtol=1e-11)
    tc_rich = (4 * tc2 - tc1) / 3.0                      # leading error ~ wbar^2/N^2
    tc_over_wE = tc_rich * KB_MEV_PER_K / W0_MEV
    r = (tc_over_wE / kay[lam]["tc_asym"] - 1.0) / lam
    dt = time.time() - t0
    # matched-convention reference at the SAME N
    match = []
    for NN, tcN in ((N1, tc1), (N2, tc2)):
        g = lambda x: h_N_eps(np.exp(x), NN) - 1.0 / lam
        wm = np.exp(brentq(g, np.log(w) - 0.08, np.log(w) + 0.08, xtol=1e-12))
        tcm = W0_MEV / (2 * np.pi * wm) / KB_MEV_PER_K
        match.append(tcN / tcm - 1.0)
    rel1 = tc1 * KB_MEV_PER_K / W0_MEV / kay[lam]["tc"] - 1
    rel2 = tc2 * KB_MEV_PER_K / W0_MEV / kay[lam]["tc"] - 1
    relR = tc_over_wE / kay[lam]["tc"] - 1
    sol_res[lam] = dict(N1=N1, N2=N2, tc1=tc1, tc2=tc2, tc=tc_over_wE, r=r)
    print("    lam=%.2f  N=(%d,%d)  raw dev vs KAY(N=oo): %+.3e, %+.3e  "
          "(wbar^2/N^2 = %.3e, %.3e)   Richardson: %+.3e   (%.1fs)"
          % (lam, N1, N2, rel1, rel2, (w / N1) ** 2, (w / N2) ** 2, relR, dt))
    _chk("solver T_c == matched-convention KAY root at the SAME N (lam=%.2f)" % lam,
         abs(match[0]) < 1e-7 and abs(match[1]) < 1e-7,
         "rel dev: N=%d %+.2e, N=%d %+.2e" % (N1, match[0], N2, match[1]))
    _chk("raw solver-vs-KAY gap at N=%d is the wbar^2/N^2 shift to 10%%" % N2,
         abs(rel2 / (w / N2) ** 2 - 1.0) < 0.10, "ratio=%.4f" % (rel2 / (w / N2) ** 2))
    _chk("Richardson agrees with KAY(N=oo) to < 5e-4 (task threshold was 1e-3)",
         abs(relR) < 5e-4, "rel dev = %+.3e, expected ~ (wbar/N1)^4 = %.1e"
         % (relR, (w / N1) ** 4))

if SLOW:
    print("\n  B.4b  --slow: lam=0.15 (wbar ~ 470) at a SINGLE truncation N = 3 wbar.")
    print("        No Richardson (N = 6 wbar would cost ~1 h of dense O(N^3) eigvals);")
    print("        the point is only that the matched-convention identity still holds.")
    lam = 0.15
    w = kay[lam]["wbar"]
    N1 = int(3 * w)
    tc_asym_K = np.exp(C0) * np.exp(-1.0 / lam) * W0_MEV / (2 * np.pi) / KB_MEV_PER_K
    t0 = time.time()
    tc1 = Solver(Spectrum.einstein(lam=lam, w0=W0_MEV), mu_star=0.0, w_cut=5 * W0_MEV,
                 n_matsubara=N1).tc(t_bracket_k=(0.70 * tc_asym_K, 1.60 * tc_asym_K), xtol=1e-11)
    g = lambda x: h_N_eps(np.exp(x), N1) - 1.0 / lam
    wm = np.exp(brentq(g, np.log(w) - 0.20, np.log(w) + 0.20, xtol=1e-12))
    tcm = W0_MEV / (2 * np.pi * wm) / KB_MEV_PER_K
    rel1 = tc1 * KB_MEV_PER_K / W0_MEV / kay[lam]["tc"] - 1
    print("        N=%d  Tc_solver=%.8f K  matched-convention KAY=%.8f K  (%.0fs)"
          % (N1, tc1, tcm, time.time() - t0))
    _chk("lam=0.15: solver T_c == matched-convention KAY root at N=%d" % N1,
         abs(tc1 / tcm - 1) < 1e-7, "rel dev = %+.2e" % (tc1 / tcm - 1))
    _chk("lam=0.15: raw solver-vs-KAY(N=oo) gap is the wbar^2/N^2 shift to 30%",
         abs(rel1 / (w / N1) ** 2 - 1.0) < 0.30,
         "gap=%+.3e, wbar^2/N^2=%.3e, ratio=%.4f  (only 30%% here: at N=3wbar the"
         " shift is 11%%, far outside the asymptotic regime)"
         % (rel1, (w / N1) ** 2, rel1 / (w / N1) ** 2))

print("\n  B.5  three-way table   T_c^asym = (w_E/2pi) e^{c0} e^{-1/lam},"
      "   r(lam) = (T_c/T_c^asym - 1)/lam")
print("\n  %5s %12s %14s %14s %14s %10s %10s %12s" %
      ("lam", "wbar", "Tc/wE (KAY)", "Tc/wE (solv)", "Tc/wE (asym)",
       "r (KAY)", "r (solv)", "(h-L-c0)L"))
for lam in LAMS:
    k = kay[lam]
    s = sol_res.get(lam)
    print("  %5.2f %12.4f %14.7e %14s %14.7e %10.6f %10s %12.6f"
          % (lam, k["wbar"], k["tc"],
             ("%.7e" % s["tc"]) if s else "     --      ", k["tc_asym"],
             k["r"], ("%.6f" % s["r"]) if s else "   --   ", k["resid"]))
print("\n  b_1 = 49 pi^2/144 - 3 = %.10f   <-- the lam->0 limit of both r columns"
      % B1_CLOSED)

# quadratic fit of r(lam) on the three smallest lambdas
xs = np.array(LAMS[:3])
ys = np.array([kay[l]["r"] for l in LAMS[:3]])
co = np.polyfit(xs, ys, 2)
print("  quadratic fit r(lam) = %.6f %+.5f lam %+.5f lam^2   (intercept vs b1: %+.2e)"
      % (co[2], co[1], co[0], co[2] - B1_CLOSED))
_chk("quadratic extrapolation of r(lam) hits b_1 to 1e-3",
     abs(co[2] - B1_CLOSED) < 1e-3, "intercept=%.6f b1=%.6f" % (co[2], B1_CLOSED))
# the NATURAL variable is 1/ln wbar, not lam:  R(wbar) := (h - ln wbar - c0) ln wbar
xs2 = np.array([1.0 / np.log(kay[l]["wbar"]) for l in LAMS[:3]])
co2 = np.polyfit(xs2, np.array([kay[l]["resid"] for l in LAMS[:3]]), 2)
print("  operator residual fitted in 1/ln wbar over wbar = %.1f, %.1f, %.1f:"
      % tuple(kay[l]["wbar"] for l in LAMS[:3]))
print("    R = %.6f %+.5f /ln wbar %+.5f /ln^2 wbar   (intercept vs b1: %+.2e;"
      "  by-product estimate b_2 ~ %+.2f, NOT a result)"
      % (co2[2], co2[1], co2[0], co2[2] - B1_CLOSED, co2[1]))
_chk("(h - ln wbar - c0) ln wbar extrapolates to b_1 to 1e-3 in the 1/ln wbar variable",
     abs(co2[2] - B1_CLOSED) < 1e-3, "intercept=%.6f b1=%.6f" % (co2[2], B1_CLOSED))
_chk("r(lam) increases monotonically as lam decreases",
     all(kay[LAMS[i]]["r"] > kay[LAMS[i + 1]]["r"] for i in range(len(LAMS) - 1)),
     "r = " + ", ".join("%.5f" % kay[l]["r"] for l in LAMS))
_chk("operator residual (h - ln wbar - c0) ln wbar -> b_1 from below, monotone",
     all(kay[LAMS[i]]["resid"] > kay[LAMS[i + 1]]["resid"] for i in range(len(LAMS) - 1))
     and kay[0.10]["resid"] < B1_CLOSED,
     "resid(lam=0.10) = %.6f < b1 = %.6f" % (kay[0.10]["resid"], B1_CLOSED))

# ===========================================================================
hdr("SUMMARY")
print("  2<tau,(I-D)tau>  ->  3 pi^2/16 - 7/4   = %.15f" % J_CLOSED)
print("  2<tau,calR tau>  ->  11 pi^2/72 - 5/4  = %.15f" % IR_CLOSED)
print("  b_1              =   49 pi^2/144 - 3   = %.15f   = (7 pi/12)^2 - 3" % B1_CLOSED)
print("  T_c = (2/pi) e^{gamma-3/2} w_E e^{-1/lam} [1 + b_1 lam + O(lam^2)]")
print("      = %.6f w_E e^{-1/lam} [1 + %.6f lam + ...]" % (pref, B1_CLOSED))
if not SLOW:
    print("\n  (no solver row for lam <= 0.15: it needs n_matsubara >~ 10 wbar = 4700 and"
          "\n   the solver's linearised step is a dense non-symmetric eigvals, O(N^3);"
          "\n   `--slow` runs lam=0.15 at N = 3 wbar as a convention check only.)")
print("\nALL %d CHECKS PASSED   (%.1f s)" % (ck, time.time() - T_START))
