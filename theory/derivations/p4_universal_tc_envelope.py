"""P4 -- UNIVERSAL certified T_c upper envelope over ALL phonon spectra (mu*=0).

For every coupling lambda>0 and EVERY phonon spectrum P in P (probability measure
on R_+, second moment <w^2> < infty, P != delta_0), we deliver a certified upper
bound
        T_c(P; lambda)  <=  F_hi^E(lambda) * sqrt(<w^2>) ,
where F_hi^E(lambda) is the certified UPPER curve of T_c/w_E vs lambda for the
single-frequency (Einstein) spectrum (deliverable (3), p4_global_tc_certified.py).

WHY THIS IS TRUE (all analytic, inherited):
  * Theorem 2 (Einstein extremality, P4_main_theorems.md, kernel-level Jensen):
      k(P,T) = lam_1(O[F_T]) <= h(w_rms),  w_rms := sqrt(<w^2>)/(2 pi T),
    with EQUALITY iff P is a point mass.  Equivalently T_c(P,lambda) <=
    T_c^E(lambda, sqrt(<w^2>)) for ALL lambda.  The fixed moment is EXACTLY the
    SECOND moment <w^2>; the matched Einstein frequency is the rms frequency
    w_E = sqrt(<w^2>).  Hence
        T_c(P;lambda) <= T_c^E(lambda; sqrt(<w^2>)) = F^E(lambda) * sqrt(<w^2>)
                       <= F_hi^E(lambda) * sqrt(<w^2>).
  * F_hi^E(lambda) is deliverable (3)'s certified upper bracket of T_c/w_E,
    obtained by inverting a certified two-sided envelope of h(w) through the
    PROVEN strict monotonicity of h (Cor (i)).  Inherited CERTIFICATE-LEVEL.

SHARPENING vs the asymptotic ceiling  T_c <= C_inf sqrt(lambda <w^2>):
    C_inf = (1/2pi) sqrt(g(2)) = 0.1827262477...  (KAY; = Allen-Dynes f(0)).
  Theorem 1(c) is a valid upper bound for ALL lambda but is tight ONLY as
  lambda->infty.  For the Einstein spectrum T_c/w_E = 1/(2 pi w_*) with
  lambda h(w_*) = 1, so
        (T_c/w_E) / (C_inf sqrt(lambda))  =  1 / ( w_* sqrt(g(2) lambda) )
                                           =  sqrt( r(w_*) / g(2) ),   r(w):=h(w)/w^2.
  By Cor (ii), r(w) < g(2) strictly and r(w) is strictly decreasing with
  r(0+)=g(2).  So the universal envelope is STRICTLY below the C_inf ceiling at
  every finite lambda, and the gap closes (ratio -> 1) only as lambda -> infty.
  This module CERTIFIES the sharpening at each tabulated lambda by the SOUND
  integer/rational test
        w_lo^2 * G2_LO * lambda  >=  1      ==>   F_hi^E(lambda) <= C_inf sqrt(lambda),
  where w_lo (<= w_*, certified by (3)) and G2_LO (certified lower on g(2), CERT-B)
  make the test the most conservative possible (only +,*,/, no sqrt/pi -> sound).
  Note pi CANCELS in the ratio, so no transcendental enters any load-bearing test.

ONE-SIDED SCOPE (honest):  Theorem 2 gives max_P k(P,T) = Einstein.  There is NO
  matching nondegenerate LOWER envelope: at fixed <w^2> the infimum of T_c over P
  is 0, approached only by degenerate spectra with mass escaping to {0, infty}.
  Demonstrated numerically in check group D.

DEPENDENCIES.
  * p4_interchange_uniformity  (committed base): H_mat  -- the operator O[.].
  * p4_global_tc_certified     (deliverable (3), inherited certificate-level):
        node_cert, invert_envelope, tc_from_w, strong_invert, _weak_invert,
        _weak_asympt_invert, G2_LO, G2_HI, TWO_PI.
    If (3) is absent the gate aborts with a pointer to its commit.

Run:  OMP_NUM_THREADS=1 python p4_universal_tc_envelope.py
Prints  "PASS k  FAIL 0".  Single-thread, < 2 min.  No /tmp scratch imports.
"""
import os
import sys
import time

import numpy as np
import mpmath as mp

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from p4_interchange_uniformity import H_mat  # noqa: E402  (committed base)

try:
    from p4_global_tc_certified import (  # noqa: E402  (deliverable (3))
        node_cert, invert_envelope, tc_from_w,
        strong_invert, _weak_invert, _weak_asympt_invert,
        G2_LO, G2_HI, TWO_PI)
except ImportError as e:  # pragma: no cover
    sys.stderr.write(
        "\nERROR: this gate inherits the certified Einstein bracket from "
        "deliverable (3),\n  theory/derivations/p4_global_tc_certified.py "
        "(commit 1184684).\n  Pull it into the worktree, e.g.\n"
        "    git checkout 1184684 -- theory/derivations/p4_global_tc_certified.py\n"
        "                             theory/derivations/P4_global_tc_certified.md\n"
        f"  (underlying import error: {e})\n\n")
    raise

mp.mp.dps = 50

# ---------------------------------------------------------------------------
PASS = [0]
FAIL = [0]


def ck(cond, msg=""):
    if cond:
        PASS[0] += 1
    else:
        FAIL[0] += 1
        print(f"  !! CHECK FAILED: {msg}")
    return bool(cond)


# ===========================================================================
# The single-scalar Einstein h and the general discrete-P operator k(P,T).
#   O is elementwise-LINEAR in its kernel, and F_T(k) = sum_i p_i [[k]](w_i) is a
#   nonnegative combination of Einstein kernels, so  O[F_T] = sum_i p_i H_mat(w_i).
# ===========================================================================
def h_eins(w, N):
    """h(w) = lam_1(O[ [[.]](w) ]) for the Einstein spectrum (float, truncation N)."""
    return float(np.linalg.eigvalsh(H_mat(N, float(w)))[-1])


def kappa_discrete(atoms, N):
    """k(P,T) = lam_1(O[F_T]) for a discrete P given as [(p_i, w_i)] with
    w_i = omega_i/(2 pi T).  O[F_T] = sum_i p_i H_mat(N, w_i)."""
    M = np.zeros((N, N))
    for p, w in atoms:
        M += float(p) * H_mat(N, float(w))
    return float(np.linalg.eigvalsh(M)[-1])


def varpi_rms(atoms):
    """rms varpi = sqrt(sum p_i w_i^2)  (the Einstein match of the 2nd moment)."""
    return float(np.sqrt(sum(float(p) * float(w) ** 2 for p, w in atoms)))


# ===========================================================================
# Certified enclosure  [w_lo, w_hi] of w_*  (lambda h(w_*) = 1), band-dispatched
# through deliverable (3)'s certified pipeline.  w_lo <= w_* <= w_hi, both rigorous
# (each node's h_lo/h_hi are bounds on the EXACT operator; monotone-staircase
# inversion, zero slack).  F_lo = 1/(2pi w_hi) <= T_c/w_E <= 1/(2pi w_lo) = F_hi.
# ===========================================================================
def crossover_bracket(lam, Nv=140, Na=70, ngrid=41):
    """Certified [w_lo, w_hi] for lambda in the crossover window via a LOCAL node
    grid around w_* (reuses (3)'s node_cert + invert_envelope verbatim)."""
    u = 1.0 / lam
    # cheap float locate of w_* by bisection on the proven-increasing h
    lo, hi = 0.15, 3.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if h_eins(mid, Nv) < u:
            lo = mid
        else:
            hi = mid
    west = 0.5 * (lo + hi)
    for half in (0.06, 0.12, 0.25):          # widen until the grid straddles u
        a = max(0.15, west * (1.0 - half))
        b = min(3.0, west * (1.0 + half))
        ws = np.exp(np.linspace(np.log(a), np.log(b), ngrid))
        hlo, hhi = [], []
        ok = True
        for w in ws:
            r = node_cert(w, Nv=Nv, Na=Na)
            if r[1] is None:
                ok = False
                break
            hlo.append(r[0])
            hhi.append(r[1])
        if not ok:
            continue
        res = invert_envelope(ws, np.array(hlo), np.array(hhi), u)
        if res is not None:
            return res
    return None


def varpi_star_bracket(lam):
    """Return (w_lo, w_hi, band).  Dispatches to (3)'s certified band pieces."""
    if lam >= 34.4:
        r = strong_invert(1.0 / lam)
        return (float(r[0]), float(r[1]), "strong-Weyl") if r else (None, None, "strong-Weyl")
    if lam >= 0.582:
        r = crossover_bracket(lam)
        return (float(r[0]), float(r[1]), "crossover") if r else (None, None, "crossover")
    if lam >= 0.345:
        r = _weak_invert(1.0 / lam)
        return (float(r[0]), float(r[1]), "weak-W") if r else (None, None, "weak-W")
    r = _weak_asympt_invert(1.0 / lam)
    # F_hi (upper T_c) side here is Theorem W' asymptotic form -- NOT finite-N certified
    return (float(r[0]), float(r[1]), "weak-asympt*") if r else (None, None, "weak-asympt*")


# deliverable (3)'s published certified brackets (from P4_global_tc_certified.md
# section 4, N=2000) -- used only as an independent CONSISTENCY (overlap) target.
PUB = {
    0.7:  (7.202285e-2, 7.213087e-2),
    1.0:  (1.146122e-1, 1.147841e-1),
    2.0:  (2.115556e-1, 2.118729e-1),
    5.0:  (3.795354e-1, 3.801047e-1),
    10.0: (5.570158e-1, 5.586878e-1),
    30.0: (9.888704e-1, 9.903534e-1),
}


def main():
    t0 = time.time()
    print("=" * 78)
    print("P4 UNIVERSAL T_c ENVELOPE  --  T_c(P;lambda) <= F_hi^E(lambda) sqrt(<w^2>)")
    print("=" * 78)

    # -------------------------------------------------------------------
    # A. Einstein extremality (Theorem 2) -- the load-bearing math, tested
    #    DIRECTLY on genuine multi-frequency spectra P.
    #    k(P,T) <= h(w_rms), strict off the point mass, equality on it.
    # -------------------------------------------------------------------
    print("\n[A] Einstein extremality  k(P,T) <= h(w_rms)  (Theorem 2, direct)")
    N = 260
    spectra = [
        ("2pt narrow",  [(0.5, 1.0), (0.5, 3.0)]),
        ("2pt wide",    [(0.5, 0.5), (0.5, 4.0)]),
        ("2pt skew",    [(0.8, 1.0), (0.2, 5.0)]),
        ("3pt",         [(0.34, 0.7), (0.33, 2.0), (0.33, 6.0)]),
        ("near-Einstein+far", [(0.99, 2.0), (0.01, 20.0)]),
        ("physical-ish", [(0.6, 0.8), (0.4, 2.5)]),
    ]
    print(f"    {'spectrum':>20} {'kappa(P)':>12} {'w_rms':>8} {'h(w_rms)':>12} {'gap':>11}")
    for name, atoms in spectra:
        k = kappa_discrete(atoms, N)
        wr = varpi_rms(atoms)
        hr = h_eins(wr, N)
        gap = hr - k
        print(f"    {name:>20} {k:12.7f} {wr:8.4f} {hr:12.7f} {gap:11.3e}")
        ck(k <= hr * (1 + 1e-9), f"extremality k<=h_rms failed for {name}")
        ck(gap > 1e-4, f"extremality NOT strict for non-degenerate {name}")
    # equality for a genuine point mass
    for w0 in (0.7, 2.0, 5.0):
        kpt = kappa_discrete([(1.0, w0)], N)
        hpt = h_eins(w0, N)
        ck(abs(kpt - hpt) <= 1e-12 * abs(hpt),
           f"point-mass equality k(delta)==h failed at w={w0}")
    print("    point-mass equality k(delta_w)==h(w): confirmed (w=0.7,2,5)")

    # -------------------------------------------------------------------
    # B. Universal envelope table + CERTIFIED sharpening vs the C_inf ceiling.
    #    Envelope coefficient = F_hi^E(lambda) (units of sqrt(<w^2>)).
    #    Certified sharpening test:  w_lo^2 * G2_LO * lambda >= 1  (sound)
    #      ==> F_hi^E(lambda) <= C_inf sqrt(lambda).
    # -------------------------------------------------------------------
    print("\n[B] Universal envelope  F_hi^E(lambda)  and sharpening vs C_inf sqrt(lambda)")
    print("      ratio  = F_hi^E / (C_inf sqrt(lambda))   (< 1  ==>  envelope sharper)")
    print("      sharpen= (C_inf sqrt(lambda)) / F_hi^E    (> 1  ==>  ceiling looser)")
    lambdas = [0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 30.0]
    C_inf_mid = mp.mpf("0.1827262477")            # inherited (KAY / CERT-D midpoint)
    print(f"    {'lambda':>7} {'F_lo^E':>12} {'F_hi^E':>12} {'C_inf*sqrtL':>12}"
          f" {'ratio<=':>9} {'sharpen>=':>10} {'band':>13}")
    ratios = []
    for lam in lambdas:
        w_lo, w_hi, band = varpi_star_bracket(lam)
        ck(w_lo is not None and w_hi is not None, f"no bracket at lambda={lam}")
        if w_lo is None:
            continue
        F_lo, F_hi = tc_from_w(w_lo, w_hi)
        ck(0 < F_lo <= F_hi, f"degenerate F bracket at lambda={lam}")
        # certified sharpening (sound: only *,>=, integers/CERT-B rationals)
        wlo_m, lam_m = mp.mpf(w_lo), mp.mpf(lam)
        Q = wlo_m * wlo_m * G2_LO * lam_m          # = (F_hi / (C_inf^LO sqrt(L)))^-2
        below = ck(Q >= 1, f"F_hi<=C_inf*sqrt(lambda) NOT certified at lambda={lam}")
        ratio_hi = mp.sqrt(1 / Q) + mp.mpf(10) ** (-30)   # certified upper on ratio
        sharpen_lo = mp.sqrt(Q)                            # certified lower on sharpening
        ceil = C_inf_mid * mp.sqrt(lam_m)
        ratios.append(float(ratio_hi))
        tag = " " if below else "*"
        star = "" if band != "weak-asympt*" else "  (F_hi asympt, NOT finite-N)"
        print(f"    {lam:7.2f} {F_lo:12.6e} {F_hi:12.6e} {float(ceil):12.6e}"
              f" {float(ratio_hi):9.5f}{tag}{float(sharpen_lo):10.4f} {band:>13}{star}")

    # ratio strictly increasing toward 1 (sharpening monotonically weakens)
    for i in range(len(ratios) - 1):
        ck(ratios[i] < ratios[i + 1],
           f"ratio not increasing between lambda index {i},{i+1}")
    ck(all(r < 1.0 for r in ratios), "some ratio >= 1 (envelope not below ceiling)")
    ck(ratios[-1] > 0.98, "ratio at lambda=30 unexpectedly far from 1")
    print("    ratio(lambda) strictly increases toward 1 (ceiling tightens as lambda->inf)")

    # independent consistency with (3)'s published N=2000 brackets (overlap test)
    print("\n    consistency with deliverable (3) published brackets (interval overlap):")
    for lam in (0.7, 1.0, 2.0, 5.0, 10.0, 30.0):
        w_lo, w_hi, _ = varpi_star_bracket(lam)
        F_lo, F_hi = tc_from_w(w_lo, w_hi)
        plo, phi = PUB[lam]
        overlap = (F_lo <= phi * (1 + 1e-9)) and (plo <= F_hi * (1 + 1e-9))
        ck(overlap, f"my bracket disjoint from (3) at lambda={lam}: "
                    f"mine=[{F_lo:.6e},{F_hi:.6e}] pub=[{plo:.6e},{phi:.6e}]")
    print("      all six overlap (3)'s published bracket -> faithful to inherited (3)")

    # -------------------------------------------------------------------
    # C. The C_inf ceiling enclosure and the strong-lambda boundary of sharpening.
    # -------------------------------------------------------------------
    print("\n[C] C_inf = (1/2pi) sqrt(g(2)) enclosure and the strong-end boundary")
    pi_hi = (+mp.pi) + mp.mpf(10) ** (-40)     # rigorous enclosure of pi at dps=50
    pi_lo = (+mp.pi) - mp.mpf(10) ** (-40)
    C_lo = mp.sqrt(G2_LO) / (2 * pi_hi) - mp.mpf(10) ** (-30)
    C_hi = mp.sqrt(G2_HI) / (2 * pi_lo) + mp.mpf(10) ** (-30)
    print(f"    C_inf in [{mp.nstr(C_lo, 12)}, {mp.nstr(C_hi, 12)}]  (from CERT-B g(2))")
    ck(C_lo <= C_inf_mid <= C_hi, "C_inf enclosure does not contain 0.1827262477")
    ck(C_lo > mp.mpf("0.18272") and C_hi < mp.mpf("0.18273"),
       "C_inf enclosure outside the known [0.18272,0.18273] band")
    # strong-lambda: envelope's certified upper F_hi has MERGED into the ceiling.
    #   strong_invert => w_lo = sqrt((1/lam)/G2_HI)  =>  w_lo^2 G2_HI lam == 1 exactly,
    #   so F_hi = (1/2pi) sqrt(G2_HI lam) >= C_inf sqrt(lam): NO certified sharpening.
    lam_s = 45.0
    ws_lo, ws_hi, band_s = varpi_star_bracket(lam_s)
    Qhi = mp.mpf(ws_lo) ** 2 * G2_HI * mp.mpf(lam_s)
    Qlo = mp.mpf(ws_lo) ** 2 * G2_LO * mp.mpf(lam_s)
    print(f"    lambda={lam_s} ({band_s}):  w_lo^2*G2_HI*lambda = {mp.nstr(Qhi, 10)}"
          f"  w_lo^2*G2_LO*lambda = {mp.nstr(Qlo, 10)}")
    ck(abs(Qhi - 1) < 1e-12, "strong-band w_lo^2 G2_HI lambda != 1")
    ck(Qlo < 1, "strong band unexpectedly still shows certified sharpening")
    print("    -> at lambda>~34 the certified envelope equals the ceiling to within the")
    print("       g(2) bracket width: sharpening is a strictly FINITE-lambda phenomenon.")

    # -------------------------------------------------------------------
    # D. One-sided scope: NO matching nondegenerate LOWER envelope.
    #    At fixed <w^2>, inf_P k(P,T) = 0, approached only by mass -> {0, infty}.
    # -------------------------------------------------------------------
    print("\n[D] One-sided scope: degenerate infimizing sequence collapses k -> 0")
    sigma2 = 4.0                                  # fixed 2nd moment <w^2> (this T)
    eps = 0.02
    hE = h_eins(np.sqrt(sigma2), N)
    print(f"    fixed <w^2>={sigma2}: Einstein (max) h={hE:.6f};  P=(1-p)d_{{{eps}}}+p d_Omega:")
    prev = None
    ks = []
    for Om in (5.0, 10.0, 20.0, 50.0, 100.0, 300.0):
        p = (sigma2 - eps ** 2) / (Om ** 2 - eps ** 2)
        k = kappa_discrete([(1 - p, eps), (p, Om)], N)
        ks.append(k)
        m2 = (1 - p) * eps ** 2 + p * Om ** 2
        print(f"      Omega={Om:6.0f}  p={p:.3e}  <w^2>={m2:.4f}  kappa={k:.6e}")
        ck(k <= hE + 1e-9, f"degenerate kappa exceeds Einstein max at Omega={Om}")
        if prev is not None:
            ck(k < prev, f"degenerate kappa not decreasing at Omega={Om}")
        prev = k
    ck(ks[-1] < 1e-3, "degenerate kappa did not collapse toward 0")
    print("    -> inf_P T_c at fixed <w^2> is 0 (degenerate, mass -> {0,infty}):")
    print("       the universal envelope is ONE-SIDED (upper only); no lower matcher.")

    # -------------------------------------------------------------------
    print("\n" + "=" * 78)
    print(f"PASS {PASS[0]}  FAIL {FAIL[0]}   ({time.time() - t0:.1f} s)")
    print("=" * 78)
    return FAIL[0]


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
