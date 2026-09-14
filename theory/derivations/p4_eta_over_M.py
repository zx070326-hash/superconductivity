r"""
P4 . ETA/M THEOREM ("why hydrogen") -- numerical gate for `P4_eta_over_M.md`.

v0.2 (2026-09-15).  Consistency checks (NOT proofs) for the eta/M theorem:

  (i)   the McMillan-Hopfield chain  lambda<w^2> = eta/M  is self-consistent on
        explicitly constructed alpha^2F objects, and eta is isotope-invariant
        under harmonic rigid-ion scaling (Stage 1);
  (ii)  the bound  Tc <= C_inf sqrt(eta/M),  C_inf = sqrt(g(2))/2pi, holds for
        H3S / LaH10 / YH6 / PdH against SSCHA + experimental Tc (Stage 2);
  (iii) the isotope exponent of the BOUND is exactly 1/2 (numerical derivative
        w.r.t. M), and the partial-substitution exponent equals the closed form
        (1/2)*(eta_H/M_H)/(eta/M) (Stage 3);
  (iv)  the all-lambda version  Tc <= G(lambda) sqrt(eta/M)  is tighter than the
        asymptotic version and consistent with it: G increasing, G < C_inf,
        G -> C_inf, plus the explicit error bound (Stage 4);
  (v)   the room-temperature verdicts respect the ONE-SIDED rule of section 4.A:
        a starred row (w_log standing in for w_2) yields only a LOWER ESTIMATE
        of the bound, so it may return "clears" but never an exclusion.  The
        B column quotes the CERTIFIED bracket F_hi^E imported from the merged
        P4_universal_tc_envelope, not the float64 point value of G (Stage 4b);
  (vi)  the two arithmetic facts behind the u_+ branch exclusion lemma that
        makes (2.2) legitimate: the rational inequality at both CERT-D ends of
        g(2), and the first-basis Rayleigh bound h >= H_00 = w^2/(1+w^2)
        (Stage 4c).

All theorems are proven in `P4_eta_over_M.md`; this file only checks that the
numerics agree.  UNITS ARE THE HISTORICAL FAILURE MODE OF THIS PROJECT
(`p4_tc_sandwich.py` once divided by an extra 2pi and produced an "upper bound"
below the theory value), so Stage 0 asserts every conversion in the chain
meV <-> K <-> varpi and prints a direction sentinel: an upper bound that comes
out BELOW the quantity it bounds is a bug, not a discovery.

Numerical discipline: numpy float64 eigenvalues (never mpmath eigensolvers);
mpmath only for the analytic constants (effective dps is 30: the imported
`p4_interchange_uniformity` resets it; every asserted constant passes at that
precision).  h(varpi) via H_mat imported VERBATIM from
`p4_interchange_uniformity.py` (nothing re-derived here).
NOTE: the per-stage lines print CUMULATIVE counters; the per-stage increments
are 8 / 70 / 3 / 30 / 8 / 20 / 56 / 27 / 16 (sum 238).
Runtime target < 5 min; measured ~120 s single-process.
"""

import os
import sys
import time

import numpy as np
from mpmath import mp
from scipy.optimize import brentq

mp.dps = 40

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from p4_interchange_uniformity import H_mat          # noqa: E402  (verbatim reuse)
import p4_tc_table as TCT                           # noqa: E402  (LaH10 r values)

# ==============================================================================
# check counter (style of p4_moment_threshold.py / p4_extremal_constraints.py)
# ==============================================================================
_passed = 0
_total = 0
_failures = []


def ck(cond, msg):
    """Count-style check: tally, report failures, never abort silently."""
    global _passed, _total
    _total += 1
    cond = bool(cond)
    if cond:
        _passed += 1
    else:
        _failures.append(msg)
        print(f"  [FAIL] {msg}")
    return cond


def stage(tag):
    print(f"\n=== {tag} ===")


# ==============================================================================
# constants -- every one of them cited, none of them fitted
# ==============================================================================
# g(2) analytic enclosure, CERT-D (P4_g2_analytic.md; P4_main_theorems.md sec.4.1).
G2_LO = mp.mpf('1.318140211762671')
G2_HI = mp.mpf('1.318140212373478')

# Published enclosure of (1/2pi)sqrt(g(2)) (P4_main_theorems.md sec.2, note).
DOC_C_LO = mp.mpf('0.182726247746')
DOC_C_HI = mp.mpf('0.182726247790')

# CODATA-2018 (exact SI definition): k_B = 8.617333262e-5 eV/K.
KB_EV_PER_K = mp.mpf('8.617333262e-5')
MEV_TO_K = float(mp.mpf(1) / (mp.mpf(1000) * KB_EV_PER_K))   # 1 meV -> K

# CODATA-2018: hbar*c = 197.3269804 MeV*fm = 1973.269804 eV*Angstrom;
# atomic mass unit in energy units: 931.49410242 MeV.
HBARC_EV_A = mp.mpf('1973.269804')
AMU_EV = mp.mpf('931.49410242e6')
# eta[eV/A^2] = (M[amu] * lambda<w^2>[eV^2]) * amu/(hbar c)^2
ETA_FAC = float(AMU_EV / HBARC_EV_A ** 2)

ZETA4 = float(mp.pi ** 4 / 90)
C0_WEAK = float(mp.euler + mp.log(4) - mp.mpf(3) / 2)   # h ~ ln w + c0

# Bound-side constant: the UPPER end, so the printed bound is a valid bound.
C_INF = float(mp.sqrt(G2_HI) / (2 * (+mp.pi)))          # +mp.pi: force eval
G2_HI_F = float(G2_HI)

N_EIG = 1200          # truncation for h^(N); Stage 5 certifies convergence


# ==============================================================================
# core objects
# ==============================================================================
def h_of(varpi, N=N_EIG):
    """h^(N)(varpi) = lambda_1(H^(N)(varpi)).  h^(N) <= h (compression)."""
    if float(varpi) <= 0.0:
        return 0.0
    return float(np.linalg.eigvalsh(H_mat(N, float(varpi)))[-1])


def varpi_c(lam, N=N_EIG):
    """Solve h^(N)(varpi) = 1/lambda.  h^(N) strictly increasing (Cor. (i))."""
    tgt = 1.0 / lam
    return brentq(lambda w: h_of(w, N) - tgt, 1e-6, 1e7,
                  xtol=1e-14, rtol=8.9e-16, maxiter=100)


def G_of_lambda(lam, N=N_EIG):
    """G(lambda) := 1/(2 pi varpi_c(lambda) sqrt(lambda)), so that
    Tc^E = G(lambda) * sqrt(lambda <w^2>).  Theorem B's constant.

    DIRECTION WARNING.  This is a float64 POINT value at N=N_EIG, and since
    h^(N) <= h (compression) we get varpi_c^(N) >= varpi_c, hence
    G^(N) <= G_true: the value returned is a LOWER ESTIMATE of the true G,
    NOT an upper bound.  It may be used one-sidedly (to certify "clears"),
    never as the legal Tc ceiling.  For a legal upper bound use
    F_HI_E[lam]/sqrt(lam) -- the certified bracket below.
    """
    return 1.0 / (2 * np.pi * varpi_c(lam, N) * np.sqrt(lam))


# ==============================================================================
# CERTIFIED Einstein upper bracket F_hi^E(lambda), imported from the MERGED
# universal envelope P4_universal_tc_envelope.{md,py} (which inherits it from
# deliverable (3), p4_global_tc_certified.py).  This is the step (3) that this
# document's own Theorem B lacks:
#
#     T_c(P; lambda) <= F_hi^E(lambda) * sqrt(<w^2>)
#                     = (F_hi^E(lambda)/sqrt(lambda)) * sqrt(eta/M),
#
# using <w^2> = eta/(lambda M).  So the LEGAL Theorem-B-shaped constant is
# G_cert(lambda) := F_hi^E(lambda)/sqrt(lambda)  >=  G(lambda).
#
# Values below were produced by running the merged gate's certified pipeline
#   varpi_star_bracket(lam) -> tc_from_w(w_lo, w_hi) -> (F_lo^E, F_hi^E)
# at exactly the seven material lambda of MATERIALS (the merged doc's own table
# stops at lambda in {0.3,0.5,1,2,5,10,30}, which does not cover them).
# Reproduce with theory/notes/_scratch_ir/etaM_fhi_compute.py.
# Band: lambda=0.40 -> "weak-W" (finite-N certified); the rest -> "crossover"
# (finite-N Kato-Temple).  NONE falls in the non-certified "weak-asympt*" band
# (lambda < 0.345), so every row here is finite-N certified.
#         lambda: (F_lo^E,         F_hi^E,         band)
F_CERT_E = {
    0.40: (0.022053692056437223, 0.023927460043843720, "weak-W"),
    1.71: (0.187871119698441280, 0.188436259161995200, "crossover"),
    1.73: (0.189585699582598970, 0.190155996713587350, "crossover"),
    1.84: (0.198775736212639550, 0.199373678105523640, "crossover"),
    2.06: (0.216074751742515560, 0.216724731204411100, "crossover"),
    2.67: (0.258347369451086400, 0.259124510152735970, "crossover"),
    3.62: (0.313335817616550600, 0.314278370341867100, "crossover"),
    1.00: (0.114482008870714680, 0.114826384851352980, "crossover"),
}


def G_cert(lam):
    """LEGAL upper-bound analogue of G(lambda): F_hi^E(lambda)/sqrt(lambda).

    Tc <= G_cert(lambda) * sqrt(eta/M) is a genuine ceiling (certificate-level,
    inherited), unlike G_of_lambda which is a one-sided lower estimate.
    """
    return F_CERT_E[round(lam, 2)][1] / np.sqrt(lam)


def G_cert_lo(lam):
    """Certified LOWER estimate of the same constant: F_lo^E(lambda)/sqrt(lambda).
    Used for the "clears 300 K" direction, which needs a lower estimate."""
    return F_CERT_E[round(lam, 2)][0] / np.sqrt(lam)


def eta_over_M_meV2(lam, omega2_meV):
    """lambda<w^2> = eta/M in meV^2 (the theorem's invariant combination)."""
    return lam * omega2_meV ** 2


def eta_eV_per_A2(M_amu, lam, omega2_meV):
    """eta = M * lambda<w^2> in eV/A^2 -- the unit AD Table III/IV uses."""
    return M_amu * lam * (omega2_meV * 1e-3) ** 2 * ETA_FAC


def bound_K(lam, omega2_meV, asymptotic=True):
    """Tc upper bound in K.  asymptotic=True -> C_inf (Theorem A, all lambda
    but weaker); False -> G(lambda) (Theorem B, tighter)."""
    C = C_INF if asymptotic else G_of_lambda(lam)
    return C * np.sqrt(eta_over_M_meV2(lam, omega2_meV)) * MEV_TO_K


# ==============================================================================
# STAGE 0: units.  Every conversion asserted; direction sentinel printed.
# ==============================================================================
def run_stage0():
    stage("STAGE 0: unit chain meV <-> K <-> varpi, and the analytic constants")

    ck(abs(MEV_TO_K - 11.6045181217) < 1e-9, f"meV->K wrong: {MEV_TO_K}")
    print(f"    1 meV = {MEV_TO_K:.10f} K   (k_B = 8.617333262e-5 eV/K, CODATA-2018)")

    # C_inf against the independently published enclosure.
    ck(float(DOC_C_LO) <= C_INF <= float(DOC_C_HI) + 1e-12,
       f"C_inf {C_INF} outside published enclosure")
    ck(abs(C_INF - 0.1827262477890) < 1e-12, f"C_inf digits moved: {C_INF}")
    ck(G2_LO < G2_HI, "g(2) enclosure inverted")
    print(f"    C_inf = sqrt(g(2))/2pi = {C_INF:.13f}  in [{float(DOC_C_LO)}, "
          f"{float(DOC_C_HI)}]  (CERT-D, analytic)")
    print(f"    AD 1975 Eq. (25) first-hand: Tc = 0.182 (lambda<w^2>)^1/2 "
          f"(lambda >~ 10, mu*=0)  -> ratio C_inf/0.182 = {C_INF/0.182:.5f}")

    # varpi is dimensionless and the 2pi sits where the definition says.
    # varpi = omega/(2 pi T): with omega in meV and T in K, T[K] = omega[meV]
    # *MEV_TO_K/(2 pi varpi).  Assert on a case where the answer is known.
    om, T = 100.0, 50.0
    w = om * MEV_TO_K / (2 * np.pi * T)
    ck(abs(T - om * MEV_TO_K / (2 * np.pi * w)) < 1e-12 * T, "varpi round-trip failed")
    print(f"    round-trip: omega=100 meV, T=50 K -> varpi={w:.9f} -> T={om*MEV_TO_K/(2*np.pi*w):.9f} K")

    # THE historical bug: an extra (or missing) 2pi.  Show it would be caught.
    T_bug = om * MEV_TO_K / w          # 2pi dropped
    ck(abs(T_bug / T - 2 * np.pi) < 1e-9, "2pi sentinel miscomputed")
    print(f"    2pi sentinel: dropping the 2pi would give T={T_bug:.4f} K "
          f"= {T_bug/T:.6f}x the truth (= 2pi).  Caught by this assert.")

    # eta unit factor, anchored on AD 1975 Table III (first-hand, p.918):
    # Pb: <w^2>^1/2 = 64 K, lambda = 1.55, eta = 2.4 eV/A^2, M = 207.2 amu.
    for name, w2_K, lam, eta_AD, M in (("Pb", 64.0, 1.55, 2.4, 207.2),
                                       ("a-Pb", 53.0, 1.91, 2.0, 207.2)):
        om2_meV = w2_K / MEV_TO_K
        e = eta_eV_per_A2(M, lam, om2_meV)
        ck(abs(e / eta_AD - 1.0) < 0.05,
           f"{name}: eta {e:.3f} vs AD {eta_AD} eV/A^2 off by >5%")
        print(f"    eta anchor {name:>5}: computed {e:.3f} vs AD Table III {eta_AD} "
              f"eV/A^2  ({100*(e/eta_AD-1):+.1f}%, AD prints 2 digits)")

    print(f"Stage 0: {_passed}/{_total} checks passed")


# ==============================================================================
# STAGE 1: the McMillan-Hopfield chain on explicit alpha^2F objects.
#
# We build alpha^2F from FIRST PRINCIPLES of the rigid-ion model rather than
# assuming the identity we want to test:
#     alpha^2F(w) = sum_j A_j delta(w - w_j),   A_j := eta_j /(2 M w_j)   (*)
# with eta_j MASS- AND FREQUENCY-INDEPENDENT weights (this is (A2)+(A3) of the
# derivation: |g|^2 = hbar |G|^2/(2 M w), G electronic).  Then
#     lambda      = 2 sum_j A_j/w_j
#     lambda<w^2> = 2 sum_j A_j w_j = sum_j eta_j/M     <- must come out exact.
# The test has teeth because (*) never mentions <w^2>: if the derivation's
# factorization were wrong, the two sides would disagree.
# ==============================================================================
def a2F_rigid_ion(etas, omegas, M):
    """Rigid-ion alpha^2F weights A_j = eta_j/(2 M w_j)  (derivation eq. (1.6))."""
    etas = np.asarray(etas, float)
    omegas = np.asarray(omegas, float)
    return etas / (2.0 * M * omegas)


def lam_and_moments(A, omegas):
    """(lambda, <w^2>, <w>_log) from AD 1975 Eqs. (3), (5), (A2)."""
    A = np.asarray(A, float)
    omegas = np.asarray(omegas, float)
    lam = 2.0 * np.sum(A / omegas)
    w2 = (2.0 / lam) * np.sum(A * omegas)           # AD Eq. (5) with n=2
    wlog = np.exp((2.0 / lam) * np.sum(A * np.log(omegas) / omegas))
    return lam, w2, wlog


def run_stage1():
    stage("STAGE 1: lambda<w^2> = eta/M on explicit alpha^2F (and isotope invariance)")
    rng = np.random.default_rng(20260914)

    n_ok = 0
    worst = 0.0
    for trial in range(60):
        n_modes = int(rng.integers(1, 7))
        etas = rng.uniform(0.2, 5.0, n_modes)        # eV/A^2-like, arbitrary scale
        omegas = rng.uniform(5.0, 250.0, n_modes)    # meV
        M = float(rng.uniform(1.0, 210.0))           # amu

        A = a2F_rigid_ion(etas, omegas, M)
        lam, w2, _ = lam_and_moments(A, omegas)

        lhs = lam * w2                 # lambda<w^2>, computed through the moments
        rhs = np.sum(etas) / M         # eta/M, computed from the electronic input
        rel = abs(lhs - rhs) / rhs
        worst = max(worst, rel)
        if ck(rel < 1e-12,
              f"trial {trial}: lambda<w^2>={lhs:.10e} != eta/M={rhs:.10e} (rel {rel:.2e})"):
            n_ok += 1

    print(f"    {n_ok}/60 random multi-mode rigid-ion spectra: "
          f"lambda<w^2> = eta/M exactly (worst rel. dev. {worst:.2e}, float64 noise)")

    # --- isotope invariance: harmonic scaling w -> w/sqrt(s), M -> sM leaves
    #     eta invariant and lambda invariant, and scales lambda<w^2> by 1/s.
    etas = np.array([2.3, 0.9, 4.1])
    omegas = np.array([140.0, 62.0, 210.0])
    M0 = 1.008
    A0 = a2F_rigid_ion(etas, omegas, M0)
    lam0, w20, _ = lam_and_moments(A0, omegas)
    for s in (2.0, 3.0, 1.5):                      # H -> D, H -> T, arbitrary
        Ms, oms = s * M0, omegas / np.sqrt(s)
        As = a2F_rigid_ion(etas, oms, Ms)
        lams, w2s, _ = lam_and_moments(As, oms)
        ck(abs(lams / lam0 - 1) < 1e-13,
           f"s={s}: lambda not isotope-invariant ({lams} vs {lam0})")
        ck(abs((lams * w2s) / (lam0 * w20) - 1.0 / s) < 1e-13,
           f"s={s}: lambda<w^2> did not scale as 1/s")
        ck(abs((Ms * lams * w2s) / (M0 * lam0 * w20) - 1.0) < 1e-13,
           f"s={s}: eta = M lambda<w^2> not invariant")
    print("    harmonic rigid-ion scaling (s = 2, 3, 1.5): lambda invariant, "
          "lambda<w^2> ~ 1/s, eta = M lambda<w^2> invariant")

    # --- the identity is about <w^2>, NOT <w>_log: show they genuinely differ,
    #     so that using w_log in the bound would be unjustified (sandwich lesson).
    r = np.sqrt(w20) / lam_and_moments(A0, omegas)[2]
    ck(r > 1.0, f"w_2/w_log = {r} <= 1 contradicts Cauchy-Schwarz/Lyapunov")
    print(f"    w_2/w_log = {r:.4f} > 1 (Lyapunov): w_log is NOT a substitute for w_2")

    print(f"Stage 1: {_passed}/{_total} checks passed")


# ==============================================================================
# STAGE 1b: multi-species  lambda<w^2> = sum_i eta_i/M_i  -- exact for
# species-pure modes, and NOT exact when modes mix species (the cross-term
# caveat (A6) of the derivation).  The negative control is the point: it shows
# the assumption is load-bearing rather than decorative.
# ==============================================================================
def run_stage1b():
    stage("STAGE 1b: multi-species decomposition sum_i eta_i/M_i (with its negative control)")

    # H3S-like: 3 H modes (M=1.008) + 1 S mode (M=32.06), species-pure.
    M_H, M_S = 1.008, 32.06
    eta_H = np.array([2.6, 1.9, 3.1])
    eta_S = np.array([1.4])
    om_H = np.array([150.0, 120.0, 195.0])          # meV, H-dominated (high)
    om_S = np.array([40.0])                          # meV, host (low)

    A = np.concatenate([a2F_rigid_ion(eta_H, om_H, M_H),
                        a2F_rigid_ion(eta_S, om_S, M_S)])
    om = np.concatenate([om_H, om_S])
    lam, w2, _ = lam_and_moments(A, om)

    lhs = lam * w2
    rhs = eta_H.sum() / M_H + eta_S.sum() / M_S
    ck(abs(lhs / rhs - 1) < 1e-12,
       f"species-pure: lambda<w^2>={lhs:.8e} != sum eta_i/M_i={rhs:.8e}")
    share_H = (eta_H.sum() / M_H) / rhs
    print(f"    species-pure H+S: lambda<w^2> = {lhs:.4f} = sum eta_i/M_i = {rhs:.4f} meV^2")
    print(f"    H share of eta/M = {100*share_H:.2f}%  (M_S/M_H = {M_S/M_H:.1f}x "
          f"amplifies the H channel)")

    # Negative control: one mode carries BOTH species (mixed eigenvector, cross
    # term I_H I_S e_H e_S != 0).  Then the per-species sum is NOT the answer.
    # Mixed mode: |g|^2 = hbar |I_H e_H/sqrt(M_H) + I_S e_S/sqrt(M_S)|^2/(2 w).
    I_H, I_S = 1.5, 1.1                      # deformation potentials
    e_H, e_S = np.sqrt(0.6), np.sqrt(0.4)    # mass-weighted eigenvector, |e|^2=1
    om_mix = 90.0
    amp = I_H * e_H / np.sqrt(M_H) + I_S * e_S / np.sqrt(M_S)
    A_mix = amp ** 2 / (2.0 * om_mix)        # = eta_mix/(2 w), eta_mix := amp^2
    lam_m, w2_m, _ = lam_and_moments(np.array([A_mix]), np.array([om_mix]))
    exact = lam_m * w2_m
    naive = (I_H * e_H) ** 2 / M_H + (I_S * e_S) ** 2 / M_S   # cross term dropped
    cross = 2 * I_H * e_H * I_S * e_S / np.sqrt(M_H * M_S)
    ck(abs(exact - (naive + cross)) < 1e-12 * exact,
       "mixed mode: exact != naive + cross (algebra error)")
    ck(abs(exact / naive - 1) > 0.05,
       f"negative control toothless: cross term only {100*(exact/naive-1):.2f}%")
    print(f"    mixed-species mode (|e_H|^2=0.6): exact = {exact:.4f}, "
          f"sum_i eta_i/M_i alone = {naive:.4f}  -> cross term = {cross:+.4f} "
          f"({100*(exact/naive-1):+.1f}%)")
    print("    => sum_i eta_i/M_i is EXACT only for species-pure modes (A6); the")
    print("       aggregate lambda<w^2> = eta/M needs no such assumption.")

    print(f"Stage 1b: {_passed}/{_total} checks passed")


# ==============================================================================
# STAGE 2: the four materials.
#
# DATA PROVENANCE (nothing here is fitted by this script):
#   lambda, omega_log, Tc_theory : theory/data/isotope_curated.md (SSCHA tables,
#       Errea 2015 / Errea 2020 / Troyan 2021 / Errea 2013), same rows as
#       p4_tc_sandwich.py.
#   r = omega_2/omega_log : p4_tc_table.py LAH10_POINTS (P3 first-hand
#       inversion).  ONLY LaH10 HAS IT.  H3S / YH6 / PdH have omega_log only.
#   Tc_exp : isotope_curated.md (Drozdov 2015 / Drozdov 2019 / Troyan 2021 /
#       Errea 2013 ref. 7-8).
#
# THE omega_2 PROBLEM, HANDLED ONE-SIDEDLY (P4_tc_sandwich.md lesson).
# The bound needs <w^2>^{1/2} = omega_2, not omega_log, and Lyapunov gives
# omega_2 >= omega_log.  So for the three materials lacking omega_2:
#       C sqrt(lambda) omega_log  <=  C sqrt(lambda) omega_2  =  the true bound.
# The omega_log number is therefore a LOWER ESTIMATE OF THE BOUND, never a
# bound on Tc.  But if even that lower estimate already exceeds Tc, the true
# bound exceeds it too -- a valid one-sided verification.  We report it as
# such and never call it "the bound".
# ==============================================================================
#         name        P     M_light lambda  w_log   r=w2/wlog  Tc_thy  Tc_exp  mu*
MATERIALS = [
    ("H3S",    200,   1.008,  1.84,   92.86,  None,     194.0,  203.0, 0.16),
    ("LaH10",  129,   1.008,  3.62,   76.40,  1.194,    255.3,  None,  0.10),
    ("LaH10",  163,   1.008,  2.67,   96.40,  1.075,    242.8,  250.0, 0.10),
    ("LaH10",  214,   1.008,  2.06,  115.50,  1.105,    237.9,  None,  0.10),
    ("LaH10",  264,   1.008,  1.73,  126.60,  1.095,    216.9,  None,  0.10),
    ("YH6",    165,   1.008,  1.71,  114.87,  None,     241.5,  224.0, 0.125),
    ("PdH",      0,   1.008,  0.40,   50.21,  None,       5.0,    9.0, 0.085),
]


def run_stage2():
    stage("STAGE 2: Tc <= C_inf sqrt(eta/M) for H3S / LaH10 / YH6 / PdH")

    # cross-check that the LaH10 r values really are p4_tc_table.py's
    r_ours = {P: r for n, P, _, _, _, r, _, _, _ in MATERIALS if n == "LaH10"}
    for P, r in TCT.LAH10_POINTS:
        ck(abs(r_ours[P] - r) < 1e-12,
           f"LaH10 r mismatch at {P} GPa: ours {r_ours[P]} vs p4_tc_table {r}")
    print(f"    LaH10 r = w_2/w_log cross-checked against p4_tc_table.LAH10_POINTS "
          f"({len(TCT.LAH10_POINTS)}/4 points, first-hand P3 inversion)")

    print()
    print(f"    {'system':<8}{'P':>5} {'lam':>5} {'w_log':>7} {'w_2':>7} "
          f"{'eta/M':>9} {'eta_H^ub':>9} {'bound':>8} {'G-bnd':>8} "
          f"{'Tc_thy':>7} {'Tc_exp':>7} {'b/exp':>6}")
    print(f"    {'':<8}{'GPa':>5} {'':>5} {'meV':>7} {'meV':>7} "
          f"{'meV^2':>9} {'eV/A^2':>9} {'K':>8} {'K':>8} {'K':>7} {'K':>7} {'':>6}")
    print("    " + "-" * 100)

    rows = []
    for name, P, M_l, lam, w_log, r, Tc_thy, Tc_exp, mus in MATERIALS:
        have_w2 = r is not None
        w2 = w_log * r if have_w2 else w_log      # w_log => lower estimate only
        eta_M = eta_over_M_meV2(lam, w2)          # meV^2
        eta_H_ub = eta_eV_per_A2(M_l, lam, w2)    # eV/A^2, upper est. of eta_H
        b = bound_K(lam, w2, asymptotic=True)
        bG = bound_K(lam, w2, asymptotic=False)
        rows.append((name, P, lam, w_log, r, w2, eta_M, b, bG, Tc_thy, Tc_exp))

        # --- the load-bearing assertions -------------------------------------
        # (1) an upper bound must exceed what it bounds.  Direction sentinel.
        ck(b > Tc_thy, f"{name}@{P}: bound {b:.1f} K <= Tc_theory {Tc_thy:.1f} K "
                       f"-- THIS IS A BUG, NOT A DISCOVERY")
        if Tc_exp is not None:
            ck(b > Tc_exp, f"{name}@{P}: bound {b:.1f} K <= Tc_exp {Tc_exp:.1f} K")
        # (2) the all-lambda bound must be tighter than (or equal to) C_inf's.
        ck(bG <= b * (1 + 1e-12), f"{name}@{P}: G-bound {bG:.3f} > C_inf bound {b:.3f}")
        # (3) and still a bound.
        ck(bG > Tc_thy, f"{name}@{P}: G-bound {bG:.1f} K <= Tc_theory {Tc_thy:.1f} K")

        w2s = f"{w2:7.2f}" if have_w2 else "  (none)"
        exps = f"{Tc_exp:7.1f}" if Tc_exp is not None else "      -"
        rexp = f"{b/Tc_exp:6.2f}" if Tc_exp is not None else "     -"
        flag = "" if have_w2 else "  <- w_log only: LOWER ESTIMATE of the bound"
        print(f"    {name:<8}{P:>5} {lam:>5.2f} {w_log:>7.2f} {w2s} "
              f"{eta_M:>9.1f} {eta_H_ub:>9.2f} {b:>8.1f} {bG:>8.1f} "
              f"{Tc_thy:>7.1f} {exps} {rexp}{flag}")

    print()
    print("    eta/M column = lambda<w^2> in meV^2; it is the SAME NUMBER as")
    print("    sum_i eta_i/M_i and needs no mass assignment (Stage 1).")
    print("    eta_H^ub = M_H*lambda<w^2>: the Hopfield parameter one infers IF all")
    print("    coupling sat on hydrogen; since eta_H/M_H <= sum_i eta_i/M_i, this is")
    print("    an UPPER estimate of eta_H.  Both are BACK-INFERRED from (lambda,w_2),")
    print("    not independent measurements.")

    # tightness summary, only over rows where the bound is genuine (LaH10)
    gen = [(n, P, b, t, e) for (n, P, l, wl, r, w2, em, b, bG, t, e) in rows
           if r is not None]
    ratios = [b / t for _, _, b, t, _ in gen]
    print(f"\n    tightness on the {len(gen)} rows with a genuine w_2 (LaH10): "
          f"bound/Tc_theory = {min(ratios):.2f}-{max(ratios):.2f}x")
    ck(all(x > 1 for x in ratios), "some LaH10 ratio <= 1")

    print(f"Stage 2: {_passed}/{_total} checks passed")
    return rows


# ==============================================================================
# STAGE 3: the isotope exponent OF THE BOUND is exactly 1/2 -- and this is a
# different object from the isotope exponent of the true Tc (P1's subject).
# ==============================================================================
def run_stage3():
    stage("STAGE 3: isotope exponent of the BOUND (exactly 1/2) vs of the true Tc")

    # (a) full substitution, all masses scaled by s: numerical d ln(bound)/d ln M.
    eta_tot, M0, lam = 3.4e4, 1.008, 2.67       # eta in meV^2*amu, arbitrary
    def bound_of_M(M):
        # eta fixed (A3) => lambda<w^2> = eta/M; lambda is isotope-invariant.
        return C_INF * np.sqrt(eta_tot / M) * MEV_TO_K

    for hstep in (1e-4, 1e-5, 1e-6):
        d = (np.log(bound_of_M(M0 * np.exp(hstep)))
             - np.log(bound_of_M(M0 * np.exp(-hstep)))) / (2 * hstep)
        ck(abs(-d - 0.5) < 1e-9,
           f"alpha_bound = {-d:.12f} != 1/2 at step {hstep}")
    print(f"    d ln(bound)/d ln M = {d:+.12f}  =>  alpha_bound = {-d:.12f} "
          f"(exact 1/2, all masses scaled)")

    # exact ratio check H -> D: bound must fall by exactly sqrt(2)
    ck(abs(bound_of_M(M0) / bound_of_M(2 * M0) - np.sqrt(2)) < 1e-12,
       "H->D did not scale the bound by 1/sqrt(2)")
    print(f"    bound(M)/bound(2M) = {bound_of_M(M0)/bound_of_M(2*M0):.12f} "
          f"= sqrt(2) exactly  (alpha_bound = 1/2 per fully substituted mass)")

    # (b) PARTIAL substitution: alpha_H^bound = (1/2)*(eta_H/M_H)/(sum_i eta_i/M_i).
    #     Verified against the numerical derivative on the Stage-1b spectrum.
    M_H, M_S = 1.008, 32.06
    eta_H, eta_S = 7.6, 1.4
    def u_of_MH(MH):
        return eta_H / MH + eta_S / M_S
    def bnd_partial(MH):
        return C_INF * np.sqrt(u_of_MH(MH)) * MEV_TO_K
    hstep = 1e-6
    d = (np.log(bnd_partial(M_H * np.exp(hstep)))
         - np.log(bnd_partial(M_H * np.exp(-hstep)))) / (2 * hstep)
    share = (eta_H / M_H) / u_of_MH(M_H)
    closed = 0.5 * share
    ck(abs(-d - closed) < 1e-8,
       f"partial exponent {-d:.10f} != closed form {closed:.10f}")
    print(f"    partial (H only): alpha_H^bound = {-d:.8f} = (1/2)*shareH "
          f"= (1/2)*{share:.6f} = {closed:.8f}  in [0, 1/2]")
    ck(0.0 <= closed <= 0.5, "partial exponent outside [0,1/2]")

    # (c) the distinction that must not be blurred: PdH, first-hand SSCHA.
    #     H: lambda=0.40, w_log=50.21, Tc=5.0 K;  D: lambda=0.46, w_log=37.69,
    #     Tc=6.5 K  (isotope_curated.md, Errea 2013).  Tc RISES with mass.
    lamH, wH, TcH = 0.40, 50.21, 5.0
    lamD, wD, TcD = 0.46, 37.69, 6.5
    uH, uD = lamH * wH ** 2, lamD * wD ** 2
    a_bound = -np.log(np.sqrt(uD) / np.sqrt(uH)) / np.log(2.0)
    a_true = -np.log(TcD / TcH) / np.log(2.0)
    ck(a_bound > 0 > a_true,
       f"PdH sign contrast lost: a_bound={a_bound:.3f}, a_true={a_true:.3f}")
    print(f"\n    PdH/PdD (SSCHA, w_log-based since w_2 is unavailable):")
    print(f"      ceiling  : eta/M = {uH:.1f} -> {uD:.1f} meV^2, "
          f"alpha_bound = {a_bound:+.3f}  (ceiling FALLS, H->D)")
    print(f"      true Tc  : {TcH} K -> {TcD} K, "
          f"alpha_true  = {a_true:+.3f}  (Tc RISES, H->D)")
    print(f"    => OPPOSITE SIGNS.  alpha_bound is the exponent of the CEILING;")
    print(f"       alpha_true is P1's object.  Conflating them is a category error.")
    print(f"       (Here eta itself drifts: eta_D/eta_H = {2*uD/uH:.3f} != 1, i.e.")
    print(f"       anharmonicity breaks (A3); alpha_bound = 1/2 - (1/2)dln(eta)/dln(M).)")
    ck(abs(2 * uD / uH - 1.0) > 0.1, "PdH eta drift unexpectedly small")

    print(f"Stage 3: {_passed}/{_total} checks passed")


# ==============================================================================
# STAGE 4: the all-lambda version.  Tc <= G(lambda) sqrt(eta/M),
# G(lambda) = 1/(2 pi sqrt(lambda * H^{-1}(1/lambda))), G increasing to C_inf.
# ==============================================================================
def run_stage4():
    stage("STAGE 4: all-lambda bound G(lambda) and its consistency with C_inf")
    t0 = time.time()

    lam_grid = [0.30, 0.40, 0.60, 1.00, 1.71, 2.06, 2.67, 3.62, 10.0, 30.0, 100.0]
    Gs, wcs = [], []
    for lam in lam_grid:
        wc = varpi_c(lam)
        wcs.append(wc)
        Gs.append(1.0 / (2 * np.pi * wc * np.sqrt(lam)))

    print(f"    {'lambda':>7} {'varpi_c':>10} {'G(lambda)':>11} {'G/C_inf':>9} "
          f"{'err bnd':>9}")
    for lam, wc, G in zip(lam_grid, wcs, Gs):
        # PROVEN explicit error bound (derivation sec. 3.3):
        #   G(lambda)^2 >= [g2 + sqrt(g2^2 - 16 zeta(4)/lambda)]/(8 pi^2),
        # valid for lambda > 16 zeta(4)/g2^2 ~ 9.97, equivalently
        #   C_inf/G <= 1/sqrt(1 - 8 zeta(4)/(g2*lambda)).
        thr = 16 * ZETA4 / G2_HI_F ** 2
        eb = (1.0 / np.sqrt(1 - 8 * ZETA4 / (G2_HI_F * lam))
              if lam > thr else float('nan'))
        ebs = f"{eb:9.4f}" if eb == eb else "        -"
        print(f"    {lam:>7.2f} {wc:>10.6f} {G:>11.8f} {G/C_INF:>9.6f} {ebs}")
        ck(G < C_INF * (1 + 1e-12),
           f"G({lam})={G:.9f} >= C_inf={C_INF:.9f} (violates h < g(2) varpi^2)")
        if eb == eb:
            ck(C_INF / G <= eb + 1e-9,
               f"lambda={lam}: C_inf/G={C_INF/G:.6f} exceeds proven bound {eb:.6f}")
            ck(G >= C_INF * np.sqrt(1 - 8 * ZETA4 / (G2_HI_F * lam)) - 1e-12,
               f"lambda={lam}: G below its proven lower envelope")

    ck(all(np.diff(Gs) > 0), f"G(lambda) not strictly increasing: {Gs}")
    print(f"\n    G strictly increasing over lambda in [0.30, 100]  "
          f"(Cor. (ii): h/varpi^2 decreasing)")
    ck(abs(Gs[-1] / C_INF - 1) < 5e-3,
       f"G(100)/C_inf = {Gs[-1]/C_INF:.6f} not within 5e-3 of 1")
    print(f"    G(100)/C_inf = {Gs[-1]/C_INF:.6f} -> 1 as lambda -> infinity "
          f"(so C_inf is the OPTIMAL eta/M-only constant)")
    print(f"    proven error envelope C_inf/G <= 1/sqrt(1 - 8 zeta(4)/(g(2) lambda)) "
          f"verified for lambda > {16*ZETA4/G2_HI_F**2:.2f}")

    # weak coupling: G(lambda) ~ exp(c0 - 1/lambda)/(2 pi sqrt(lambda)).
    # The approach is slow (1/ln varpi corrections), so we check the DIRECTION
    # and magnitude, not a tight limit.
    lam_w = 0.30
    G_asym = np.exp(C0_WEAK - 1.0 / lam_w) / (2 * np.pi * np.sqrt(lam_w))
    G_num = Gs[0]
    ck(0.5 < G_num / G_asym < 2.0,
       f"weak-coupling shape off: G={G_num:.4e} vs asym {G_asym:.4e}")
    print(f"    weak coupling lambda=0.30: G = {G_num:.6e}, "
          f"exp(c0-1/lambda)/(2pi sqrt(lambda)) = {G_asym:.6e} "
          f"(ratio {G_num/G_asym:.3f}; 1/ln(varpi) corrections are slow)")
    print(f"    => the all-lambda bound carries the exp(-1/lambda) factor, which")
    print(f"       C_inf sqrt(eta/M) alone cannot (PdH: 13.6 K vs 67.3 K).")

    print(f"    [stage 4 runtime {time.time()-t0:.1f} s]")
    print(f"    observed approach rate: 1 - G/C_inf ~ "
          f"{(1-Gs[-3]/C_INF)*lam_grid[-3]:.2f}/lambda, "
          f"{(1-Gs[-2]/C_INF)*lam_grid[-2]:.2f}/lambda, "
          f"{(1-Gs[-1]/C_INF)*lam_grid[-1]:.2f}/lambda "
          f"(coefficient ~0.34; the PROVEN envelope has coefficient "
          f"{8*ZETA4/G2_HI_F:.2f}, i.e. loose by ~12x but first-order correct)")
    print(f"Stage 4: {_passed}/{_total} checks passed")


# ==============================================================================
# STAGE 4b: the quantitative "room temperature needs what" statement.
# ==============================================================================
def run_stage4b():
    stage("STAGE 4b: room-temperature necessary conditions from Theorems A and B")

    print("    (a) lambda-free ceilings (Theorem A, optimal at lambda -> infinity):")
    for T in (200.0, 250.0, 300.0, 350.0):
        u = (T / (C_INF * MEV_TO_K)) ** 2
        ck(u > 0, "negative threshold")
        print(f"        Tc >= {T:>3.0f} K  requires  eta/M >= {u:>8.0f} meV^2  "
              f"(eta_H >= {eta_eV_per_A2(1.008, 1.0, np.sqrt(u)):5.2f} eV/A^2 "
              f"if all on H)")

    print("\n    (b) per-material verdicts (one-sided where w_2 is missing):")
    eta_map = {(n, P): eta_over_M_meV2(lam, w_log * r if r else w_log)
               for n, P, _, lam, w_log, r, _, _, _ in MATERIALS}
    # 300 K threshold in eta/M (A-level, Theorem A).
    U300 = (300.0 / (C_INF * MEV_TO_K)) ** 2
    for name, P, _, lam, w_log, r, _, _, _ in MATERIALS:
        u = eta_map[(name, P)]
        b = C_INF * np.sqrt(u) * MEV_TO_K
        # B column, both directions kept separate (this is the F3b fix):
        #   "clears 300 K" is a claim that the TRUE ceiling exceeds 300 K, so it
        #       needs a LOWER estimate of it   -> F_lo^E/sqrt(lambda) (certified),
        #   "NO"           is a claim that the TRUE ceiling falls below 300 K, so
        #       it needs an UPPER bound        -> F_hi^E/sqrt(lambda) (certified),
        #       AND a genuine w_2 (a starred row's eta/M is itself only a lower
        #       estimate, so it can never support the "NO" direction).
        bG_hi = G_cert(lam) * np.sqrt(u) * MEV_TO_K          # legal ceiling
        bG_lo = G_cert_lo(lam) * np.sqrt(u) * MEV_TO_K       # legal lower estimate
        bG_pt = G_of_lambda(lam) * np.sqrt(u) * MEV_TO_K     # float64 point value
        # The certified bracket must straddle the float64 point value.
        ck(bG_lo <= bG_pt <= bG_hi,
           f"{name}@{P}: float64 G point value outside certified bracket "
           f"[{bG_lo:.3f}, {bG_hi:.3f}] (got {bG_pt:.3f})")
        # A starred (r is None) row uses w_log in place of w_2, so its number is a
        # LOWER ESTIMATE OF THE BOUND, not the bound.  A lower estimate that sits
        # below 300 K excludes NOTHING (the true bound may be anywhere above it):
        # such a row can only ever return "clears" (valid direction), never "NO".
        # See P4_eta_over_M.md section 4.A and the rule line under the 5.3 table.
        if b > 300:
            verdict_A = "clears 300 K"
        elif r is None:
            verdict_A = "UNDECIDABLE (one-sided lower estimate)"
        else:
            verdict_A = "NO (bound below 300 K)"
        if bG_lo > 300:
            verdict_G = "clears 300 K"
        elif r is None:
            verdict_G = "UNDECIDABLE (one-sided lower estimate)"
        elif bG_hi < 300:
            verdict_G = "NO (certified ceiling below 300 K)"
        else:
            verdict_G = "UNDECIDABLE (300 K inside certified bracket)"
        onside = "" if r is not None else " [one-sided]"
        print(f"        {name:<6}@{P:>3} GPa: eta/M={u:>8.0f} meV^2, "
              f"A-level bound {b:>5.0f} K -> {verdict_A};  "
              f"B-level cert (lam={lam:.2f}) [{bG_lo:>5.0f}, {bG_hi:>5.0f}] K "
              f"-> {verdict_G}{onside}"
              f"   [F_hi^E={F_CERT_E[round(lam,2)][1]:.6e} "
              f"{F_CERT_E[round(lam,2)][2]}; float64 G pt {bG_pt:.0f} K]")
        # No starred row may ever produce an exclusion.
        # INVARIANT (section 4.A): a starred row may never produce an exclusion,
        # in either column.  Only "clears" (valid direction) or "UNDECIDABLE".
        if r is None:
            for col, v in (("A", verdict_A), ("B", verdict_G)):
                ck(v.startswith("clears") or v.startswith("UNDECIDABLE"),
                   f"{name}@{P}: starred row produced verdict '{v}' in column "
                   f"{col} -- a one-sided lower estimate cannot exclude")
        if r is None:
            # r-threshold: the r = w_2/w_log that would carry each level to 300 K.
            # B level uses the CERTIFIED ceiling, so the threshold is the honest
            # one for the legal bound.
            rA = np.sqrt(U300 / lam) / w_log
            rB = 300.0 / (G_cert(lam) * np.sqrt(lam) * w_log * MEV_TO_K)
            ck(rA > 0 and rB > 0, f"{name}@{P}: nonpositive r threshold")
            # consistency: plugging the threshold r back in must land exactly on 300 K
            ck(abs(C_INF * np.sqrt(eta_over_M_meV2(lam, w_log * rA)) * MEV_TO_K
                   - 300.0) < 1e-6, f"{name}@{P}: r_A does not reproduce 300 K")
            ck(abs(G_cert(lam) * np.sqrt(eta_over_M_meV2(lam, w_log * rB))
                   * MEV_TO_K - 300.0) < 1e-6,
               f"{name}@{P}: r_B does not reproduce 300 K")
            print(f"{'':>16}-> crossing 300 K needs r=w_2/w_log >= {rA:.3f} (A), "
                  f">= {rB:.3f} (B cert); empirical LaH10 range [1.075, 1.194]")

    # lambda at which LaH10@163's eta/M would clear 300 K (solve G(lambda)).
    u_La = eta_map[("LaH10", 163)]
    f = lambda L: G_of_lambda(L) * np.sqrt(u_La) * MEV_TO_K - 300.0
    lam_crit = brentq(f, 1.0, 100.0, xtol=1e-10)
    ck(f(lam_crit + 1e-6) > 0 > f(lam_crit - 1e-6) or
       f(lam_crit - 1e-6) > 0 > f(lam_crit + 1e-6), "lambda_crit bracket wrong")
    # The margin quoted here is a LAMBDA margin: how much larger the actual
    # lambda is than the lambda at which the B-level ceiling first reaches
    # 300 K.  (The old label called 1.037 a "margin" but that number is the
    # T_c RATIO bound(2.67)/300 K -- a different quantity.  Both are printed.)
    lam_margin = 2.67 / lam_crit
    tc_ratio = G_of_lambda(2.67) * np.sqrt(u_La) * MEV_TO_K / 300.0
    ck(abs(lam_margin - 2.67 / lam_crit) < 1e-12, "lambda margin miscomputed")
    ck(lam_margin > 1.0, "LaH10@163 lambda margin should exceed 1")
    ck(abs(lam_margin - tc_ratio) > 0.1,
       "lambda margin and Tc ratio must not be conflated (they differ)")
    print(f"\n    LaH10 @163 GPa (eta/M = {u_La:.0f} meV^2) clears 300 K for "
          f"lambda >= {lam_crit:.4f} (actual lambda = 2.67: "
          f"lambda margin = 2.67/{lam_crit:.4f} = {lam_margin:.4f}; "
          f"separately, the Tc ratio bound/300 K = {tc_ratio:.3f})")
    ck(2.67 > lam_crit, "LaH10@163 should clear the room-T necessary condition")

    # (c) LaH10 -> LaD10 harmonic ceiling shift (full H -> D, eta invariant).
    b_H = C_INF * np.sqrt(u_La) * MEV_TO_K
    b_D = C_INF * np.sqrt(u_La / 2) * MEV_TO_K
    ck(abs(b_H / b_D - np.sqrt(2)) < 1e-12, "H->D ceiling not sqrt(2)")
    print(f"    harmonic full H->D: LaH10 ceiling {b_H:.0f} K -> LaD10 "
          f"{b_D:.0f} K (x 1/sqrt(2)); SSCHA LaD10 theory 180.4 K "
          f"(isotope_curated.md) sits below both, as it must.")

    # (d) F3b: the certified import is faithful.  At lambda=1 the merged gate's
    # bracket must straddle this document's float64 point value of G.
    print("\n    (d) certified bracket F_hi^E (merged P4_universal_tc_envelope) "
          "vs float64 G:")
    G1 = G_of_lambda(1.0)
    F_lo_1, F_hi_1, _ = F_CERT_E[1.00]
    ck(F_lo_1 <= G1 <= F_hi_1,
       f"lambda=1: float64 G={G1:.7f} outside certified [{F_lo_1:.7f},{F_hi_1:.7f}]")
    ck(abs(F_hi_1 - 0.1148264) < 1e-6, "lambda=1 F_hi^E not the merged value")
    ck(abs(G1 - 0.1146195) < 1e-6, "lambda=1 float64 G drifted")
    print(f"        lambda=1: certified [{F_lo_1:.7f}, {F_hi_1:.7f}] straddles "
          f"float64 G={G1:.7f}  (F_hi^E >= G: {F_hi_1 >= G1})")
    for lam in sorted(k for k in F_CERT_E if k != 1.00):
        F_lo_c, F_hi_c, band = F_CERT_E[lam]
        Gp = G_of_lambda(lam) * np.sqrt(lam)     # float64 F^E point value
        ck(F_lo_c <= Gp <= F_hi_c,
           f"lambda={lam}: float64 F^E={Gp:.7e} outside certified bracket")
        ck(F_lo_c < F_hi_c, f"lambda={lam}: degenerate certified bracket")
        ck(band in ("crossover", "weak-W"),
           f"lambda={lam}: band {band} is not finite-N certified")
        print(f"        lambda={lam:4.2f}: F^E in [{F_lo_c:.7e}, {F_hi_c:.7e}] "
              f"({band}), float64 {Gp:.7e}, width "
              f"{(F_hi_c-F_lo_c)/F_hi_c:.2e} rel")

    print(f"Stage 4b: {_passed}/{_total} checks passed")


# ==============================================================================
# STAGE 4c: the (2.2) u_+ branch exclusion lemma (F2).  Fully analytic; this
# gate only certifies the two arithmetic facts the proof leans on.
#
#   (i)  the rational inequality  16 zeta(4) - 8 zeta(4) g(2) - g(2)^2 > 0,
#        which is EQUIVALENT to  u_0/(1+u_0) > 1/lambda*,  evaluated at BOTH
#        ends of the CERT-D enclosure of g(2);
#   (ii) the first-basis Rayleigh bound  h(varpi) >= H_00 = varpi^2/(1+varpi^2),
#        including the claim that H_00 equals it EXACTLY and N-independently.
# ==============================================================================
def run_stage4c():
    stage("STAGE 4c: (2.2) u_+ branch exclusion lemma (F2)")

    z4 = mp.zeta(4)
    print("    (i) rational inequality at both CERT-D ends of g(2):")
    for tag, g2 in (("LO", G2_LO), ("HI", G2_HI)):
        expr = 16 * z4 - 8 * z4 * g2 - g2 ** 2
        ck(expr > 0, f"u_+ exclusion rational inequality FAILS at {tag} end")
        u0 = g2 / (8 * z4)
        lhs = u0 / (1 + u0)
        inv_lam_star = g2 ** 2 / (16 * z4)
        ck(lhs > inv_lam_star,
           f"u_0/(1+u_0) > 1/lambda* FAILS at {tag} end")
        print(f"        {tag}: 16z4-8z4*g2-g2^2 = {mp.nstr(expr, 9)} > 0;  "
              f"u_0 = {mp.nstr(u0, 9)},  u_0/(1+u_0) = {mp.nstr(lhs, 9)} > "
              f"1/lambda* = {mp.nstr(inv_lam_star, 9)}")
    lam_star = float(16 * z4 / G2_HI ** 2)
    ck(abs(lam_star - 9.9668) < 1e-3, "lambda* drifted from 9.9668")
    print(f"        => lambda* = 16 zeta(4)/g(2)^2 = {lam_star:.4f}; (2.2) is "
          f"VACUOUS for lambda < lambda* (no real roots), so it says nothing "
          f"in the physical range lambda in (0.5, 3).")

    print("\n    (ii) first-basis Rayleigh bound h(varpi) >= H_00 = "
          "varpi^2/(1+varpi^2), exact and N-independent:")
    for w in (0.05, 0.3, 1.0, 2.5, 8.0):
        pred = w * w / (1.0 + w * w)
        for N in (200, 600):
            H = H_mat(N, float(w))
            ck(abs(H[0, 0] - pred) < 1e-14,
               f"H_00 != varpi^2/(1+varpi^2) at varpi={w}, N={N}")
            lam_max = float(np.linalg.eigvalsh(H)[-1])
            ck(lam_max >= H[0, 0] - 1e-14,
               f"lambda_max < H_00 at varpi={w}, N={N} (Rayleigh violated)")
        print(f"        varpi={w:5.2f}: H_00 = {H[0,0]:.12f} = "
              f"varpi^2/(1+varpi^2) (N=200 and 600 identical); "
              f"lambda_max = {lam_max:.12f} >= H_00")

    # the lemma's punchline, numerically witnessed at the parabola vertex
    u0f = float(G2_HI / (8 * z4))
    w0 = np.sqrt(u0f)
    h_at_vertex = h_of(w0)
    ck(h_at_vertex >= u0f / (1 + u0f),
       "h(sqrt(u_0)) below its own Rayleigh bound")
    ck(h_at_vertex > float(G2_HI ** 2 / (16 * z4)),
       "h(sqrt(u_0)) > 1/lambda* FAILS -- u_+ branch NOT excluded")
    print(f"\n        vertex u_0 = {u0f:.10f}, sqrt(u_0) = {w0:.10f}: "
          f"h = {h_at_vertex:.10f} >= u_0/(1+u_0) = {u0f/(1+u0f):.10f} > "
          f"1/lambda* = {float(G2_HI**2/(16*z4)):.10f}")
    print("        => varpi_c < sqrt(u_0) <= sqrt(u_+): the u_+ branch is "
          "excluded, so varpi_c^2 <= u_- and (2.2) stands.")

    print(f"Stage 4c: {_passed}/{_total} checks passed")


# ==============================================================================
# STAGE 5: the analytic inputs the two theorems rest on, spot-checked.
# ==============================================================================
def run_stage5():
    stage("STAGE 5: analytic inputs (Cor. (i)-(iii)) and h^(N) convergence")

    wg = np.array([0.05, 0.1, 0.3, 0.6, 1.0, 2.0, 5.0, 10.0, 20.0])
    hv = np.array([h_of(w) for w in wg])

    ck(np.all(np.diff(hv) > 0), "h not strictly increasing (Cor. (i))")
    ck(np.all(hv < G2_HI_F * wg ** 2), "h >= g(2) varpi^2 somewhere (Cor. (ii))")
    ratio = hv / wg ** 2
    ck(np.all(np.diff(ratio) < 0), f"h/varpi^2 not decreasing (Cor. (ii)): {ratio}")
    print(f"    h increasing, h < g(2)varpi^2, h/varpi^2 decreasing on "
          f"{len(wg)} grid points (Cor. (i),(ii))")
    print(f"    h/varpi^2: {ratio[0]:.6f} (varpi=0.05) -> {ratio[-1]:.6e} "
          f"(varpi=20), sup = g(2) = {G2_HI_F:.9f}")

    # Cor. (iii) lower half: h(w) >= ln w + gamma + ln4 - 3/2 - 1/(12w) - 3/(16w^2)
    for w, h in zip(wg, hv):
        if w >= 2.0:
            lb = np.log(w) + C0_WEAK - 1.0 / (12 * w) - 3.0 / (16 * w ** 2)
            ck(h >= lb - 1e-12, f"Cor.(iii) lower bound violated at varpi={w}: "
                                f"h={h:.9f} < {lb:.9f}")
    print(f"    Cor. (iii) lower envelope ln(varpi)+gamma+ln4-3/2-1/(12varpi)"
          f"-3/(16varpi^2) holds for varpi >= 2")

    # small-varpi quadratic: h(w)/w^2 -> g(2), and h >= g(2)w^2 - 4 zeta(4) w^4
    for w in (0.02, 0.05, 0.1):
        h = h_of(w)
        ck(h >= G2_HI_F * w ** 2 - 4 * ZETA4 * w ** 4 - 1e-14,
           f"h({w}) below g(2)w^2 - 4zeta(4)w^4")
    h_small = h_of(0.02)
    ck(abs(h_small / (G2_HI_F * 0.02 ** 2) - 1) < 2e-3,
       f"h/g(2)w^2 at varpi=0.02 = {h_small/(G2_HI_F*0.02**2)}")
    print(f"    h(0.02)/(g(2)*0.02^2) = {h_small/(G2_HI_F*0.02**2):.8f} -> 1 "
          f"(quadratic regime; g(2) is the curvature, CERT-D)")

    # h^(N) convergence at the largest varpi actually used (lambda=0.30 -> ~6.8)
    for w in (0.6, 6.8, 20.0):
        h1, h2 = h_of(w, 600), h_of(w, 1200)
        rel = abs(h2 - h1) / h2
        ck(rel < 3e-6, f"h^(N) not converged at varpi={w}: rel {rel:.2e}")
        print(f"    varpi={w:>5}: h^(600)={h1:.10f}, h^(1200)={h2:.10f}, "
              f"rel diff {rel:.2e}  (h^(N) <= h, monotone in N)")

    ck(abs(ZETA4 - 1.0823232337111382) < 1e-14, "zeta(4) wrong")
    ck(abs(C0_WEAK - 0.4635100260) < 1e-9, f"c0 = {C0_WEAK} wrong")
    print(f"    zeta(4) = {ZETA4:.13f}, c0 = gamma+ln4-3/2 = {C0_WEAK:.10f}")

    print(f"Stage 5: {_passed}/{_total} checks passed")


# ==============================================================================
def main():
    t0 = time.time()
    print(__doc__.strip().splitlines()[0])
    print(f"[N_EIG = {N_EIG}, numpy float64 eigenvalues, mpmath dps={mp.dps} "
          f"for constants only]")
    run_stage0()
    run_stage1()
    run_stage1b()
    run_stage2()
    run_stage3()
    run_stage4()
    run_stage4b()
    run_stage4c()
    run_stage5()

    print("\n" + "=" * 78)
    if _failures:
        print(f"RESULT: {_passed}/{_total} checks passed, {len(_failures)} FAILED")
        for f in _failures:
            print(f"  - {f}")
    else:
        print(f"RESULT: ALL {_passed}/{_total} CHECKS PASSED")
    print(f"total runtime {time.time()-t0:.1f} s")
    print("=" * 78)
    return 1 if _failures else 0


if __name__ == "__main__":
    sys.exit(main())
