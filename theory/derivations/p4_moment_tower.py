r"""
P4 - Theorem W'' (moment tower): gate/probe for P4_weak_coupling_series.md sec.2b.

Reproduces, from the SCALING-LIMIT integrals only (mpmath, dps=60), every
number folded into the "Theorem W'' (moment tower)" section:

  Reduction (exact, non-asymptotic):  K := M + I = R + (I - D),
     F(h) = 2<tau,(hI-M)^-1 tau> = 2<tau,((h+1)I - K)^-1 tau> = 1,
     delta := h-(2B-1) = 2 sum_{j>=1} kappa_j / y^j,   y := h+1 = L + C + delta,
     kappa_j := <tau, K^j tau>,   L := ln(varpi),  C := gamma+ln4-1/2 = c0+1.

  Scaling image (both gamma-free):  (K tau)(s) -> Upsilon_G(s) + d(s) G(s) =: Phi_K(s),
     Upsilon_G = (E4) from P4_b1_closed_form.md,  d(s)=1-atan(s)/s,  G(s)=1/(1+s^2).

  Per-order (A3): b1 = 2 kappa1,  b2 = 2 kappa2 - b1*C,  b3 = 2 kappa3 - 4 kappa2 C + 2 kappa1(C^2 - b1).

  kappa1 = 1/2 int_0^inf G Phi_K ds/s      (= b1/2 exactly, closed form 49 pi^2/144 - 3)
  kappa2 = 1/2 int_0^inf Phi_K^2  ds/s     (weight-3 period, NO known Q[pi^2] closed form)
  kappa3 = 1/2 int_0^inf Phi_K (K Phi_K) ds/s,  (K phi)(s)=1/2 int R(s,t) phi(t) dt/t + d(s) phi(s)

Physical T_c translation (careful h<->T_c inversion, extends P4_b1_closed_form.md sec.5):
  h = L + c0 + b1/L + b2/L^2 + ... ,  lam*h = 1,  L = ln(omega_E/2pi T_c):
    T_c = (omega_E/2pi) e^{c0} e^{-1/lam} [ 1 + b1 lam + b2Tc lam^2 + ... ],
    b2Tc = b1 c0 + b1^2/2 + b2 = 2 kappa2 - b1 + b1^2/2   (gamma cancels; still not in Q[pi^2]).

Standalone: imports only mpmath. Counter `ck`. Runtime a few minutes (kappa3 nested 2D).
"""
import mpmath as mp

mp.mp.dps = 60

ck = [0]


def chk(cond, msg=""):
    assert cond, msg
    ck[0] += 1
    return cond


# ---------------------------------------------------------------- kernels
def G(s):
    return 1 / (1 + s**2)


def d(s):
    return 1 - mp.atan(s) / s


def Ups(s):
    # (E4): Upsilon_G(s) = [3 ln(1+s^2) + 2(s^2-2) atan(s)/s + (4+s^2)] / [2(1+s^2)(4+s^2)]
    return (3 * mp.log(1 + s**2) + 2 * (s**2 - 2) * mp.atan(s) / s + (4 + s**2)) \
        / (2 * (1 + s**2) * (4 + s**2))


def PhiK(s):
    return Ups(s) + d(s) * G(s)


def Rk(s, t):
    # R(s,sigma) = 2 s^2 t^2 (5+s^2+t^2) / [(1+(s-t)^2)(1+(s+t)^2)(1+s^2)(1+t^2)]
    return 2 * s**2 * t**2 * (5 + s**2 + t**2) \
        / ((1 + (s - t)**2) * (1 + (s + t)**2) * (1 + s**2) * (1 + t**2))


# split points that straddle the kernel structure (0, poles of atan growth, infinity)
BRK = [0, mp.mpf(1) / 4, mp.mpf(1) / 2, 1, 2, 5, 20, mp.inf]

# reference numbers from the delivered derivation (branch worktree-agent-a098766863e5ef484,
# P4_strat_pi2_tower.md, pi2_tower_final.py) -- this probe RE-DERIVES them independently.
KAPPA2_REF = mp.mpf('0.09873042407289524101111089465604911963035')
B2_REF = mp.mpf('-0.327071467533777298057653195349')


def main():
    print("=== P4 Theorem W'' (moment tower) gate  [mpmath dps=%d] ===" % mp.mp.dps)

    # ---- exact analytic constants -------------------------------------------
    b1 = 49 * mp.pi**2 / 144 - 3                      # = (7 pi/12)^2 - 3
    C = mp.euler + mp.log(4) - mp.mpf(1) / 2          # = gamma + ln4 - 1/2
    c0 = C - 1                                        # = gamma + ln4 - 3/2
    chk(abs(c0 - (mp.euler + mp.log(4) - mp.mpf(3) / 2)) < mp.mpf(10)**-50,
        "c0 identity")
    print("  b1 = 49pi^2/144-3 =", mp.nstr(b1, 30))
    print("  C  = gamma+ln4-1/2=", mp.nstr(C, 30), " (c0 = C-1)")

    # ---- kappa_1 from the scaling integral: must equal b1/2 exactly ---------
    k1 = mp.quad(lambda s: G(s) * PhiK(s) / (2 * s), BRK)
    print("  kappa_1            =", mp.nstr(k1, 40))
    chk(abs(k1 - b1 / 2) < mp.mpf(10)**-45,
        "kappa1 != b1/2 : %s" % mp.nstr(k1 - b1 / 2, 6))
    print("  kappa_1 - b1/2     =", mp.nstr(k1 - b1 / 2, 4), " (reduction exact)  OK")

    # ---- kappa_2 from the scaling integral (weight-3 period) ---------------
    k2 = mp.quad(lambda s: PhiK(s)**2 / (2 * s), BRK)
    print("  kappa_2            =", mp.nstr(k2, 40))
    chk(abs(k2 - KAPPA2_REF) < mp.mpf(10)**-38,
        "kappa2 mismatch vs delivered value: %s" % mp.nstr(k2 - KAPPA2_REF, 6))
    print("  kappa_2 - ref      =", mp.nstr(k2 - KAPPA2_REF, 4), " OK")

    # ---- b2 = 2 kappa2 - b1 C  (h-expansion 1/L^2 coefficient) --------------
    b2 = 2 * k2 - b1 * C
    print("  b2 = 2 kappa2-b1 C =", mp.nstr(b2, 30))
    chk(abs(b2 - B2_REF) < mp.mpf(10)**-28, "b2 mismatch: %s" % mp.nstr(b2 - B2_REF, 6))
    # equivalent gamma-explicit form b2 = 2k2 - b1 gamma - b1 ln4 + b1/2
    b2_gamma = 2 * k2 - b1 * mp.euler - b1 * mp.log(4) + b1 / 2
    chk(abs(b2 - b2_gamma) < mp.mpf(10)**-50, "b2 gamma-explicit form")
    print("  b2 (gamma-explicit form 2k2 - b1 gamma - b1 ln4 + b1/2) agrees  OK")

    # ---- Theorem W'' : Lambda-variable coefficients beta_n = 2 kappa_n (gamma-free)
    #      delta = 2 kappa1/Lambda + 2 kappa2/Lambda^2 + (2 kappa3 - 4 kappa1^2)/Lambda^3 + ...
    beta1 = 2 * k1
    beta2 = 2 * k2
    chk(abs(beta1 - b1) < mp.mpf(10)**-45, "beta1 != b1")
    print("  Lambda-tower: beta1 = 2 kappa1 = b1 =", mp.nstr(beta1, 20),
          " ; beta2 = 2 kappa2 =", mp.nstr(beta2, 20), " (both gamma-free)  OK")

    # ---- physical T_c translation (careful h<->T_c inversion) --------------
    # h = L + c0 + b1/L + b2/L^2,  lam h = 1  =>  T_c bracket lam^2 coeff:
    b2Tc = b1 * c0 + b1**2 / 2 + b2
    b2Tc_alt = 2 * k2 - b1 + b1**2 / 2            # gamma-and-ln4-free reduced form
    chk(abs(b2Tc - b2Tc_alt) < mp.mpf(10)**-45,
        "b2Tc identity b1 c0 + b1^2/2 + b2 = 2k2 - b1 + b1^2/2 failed")
    print("  physical T_c lam^2 coeff b2Tc = b1 c0 + b1^2/2 + b2 =", mp.nstr(b2Tc, 20))
    print("      = 2 kappa2 - b1 + b1^2/2 (gamma, ln4 cancel)     =", mp.nstr(b2Tc_alt, 20), " OK")

    # numeric confirmation of the inversion: solve h(L)=1/lam for small lam,
    # T_c bracket = exp(-L + 1/lam - c0); fit (bracket-1-b1 lam)/lam^2 -> b2Tc
    def hfun(L):
        return L + c0 + b1 / L + b2 / L**2

    def bracket_coeff(lam):
        L = mp.findroot(lambda L: hfun(L) - 1 / lam, 1 / lam - c0)
        br = mp.e**(-L + 1 / lam - c0)
        return (br - 1 - b1 * lam) / lam**2

    lams = [mp.mpf('0.02'), mp.mpf('0.015'), mp.mpf('0.01')]
    cvals = [bracket_coeff(l) for l in lams]
    # linear-in-lam extrapolation lam->0 (truncation of h at b2 leaves O(lam) error)
    b2Tc_fit = cvals[-1] + (cvals[-1] - cvals[-2]) * lams[-1] / (lams[-2] - lams[-1])
    print("  T_c bracket lam^2 fit (lam=0.02,.015,.01 -> 0):", mp.nstr(b2Tc_fit, 8),
          " vs formula", mp.nstr(b2Tc, 8))
    chk(abs(b2Tc_fit - b2Tc) < mp.mpf('2e-3'),
        "T_c inversion fit disagrees: %s vs %s" % (mp.nstr(b2Tc_fit, 8), mp.nstr(b2Tc, 8)))
    # and the O(lam) coefficient of the bracket is exactly b1 (translation preserves b1)
    b1Tc_fit = (mp.e**(-mp.findroot(lambda L: hfun(L) - 1 / mp.mpf('0.005'), 200) +
                        1 / mp.mpf('0.005') - c0) - 1) / mp.mpf('0.005')
    chk(abs(b1Tc_fit - b1) < mp.mpf('5e-3'), "O(lam) bracket coeff != b1")
    print("  O(lam) bracket coeff = b1 (translation trivial at first order)   OK")

    # ---- refutation EVIDENCE (not a proof): PSLQ nulls ----------------------
    r1 = mp.pslq([b2, mp.mpf(1), mp.pi**2, mp.pi**4], tol=mp.mpf(10)**-45, maxcoeff=10**7)
    r2 = mp.pslq([2 * k2, mp.mpf(1), mp.pi**2, mp.pi**3, mp.zeta(3)],
                 tol=mp.mpf(10)**-45, maxcoeff=10**7)
    chk(r1 is None, "PSLQ(b2;1,pi2,pi4) unexpectedly found a relation: %s" % (r1,))
    chk(r2 is None, "PSLQ(2k2;1,pi2,pi3,z3) unexpectedly found a relation: %s" % (r2,))
    print("  PSLQ(b2;1,pi2,pi4)          =", r1, " (60-digit null)")
    print("  PSLQ(2 kappa2;1,pi2,pi3,z3) =", r2, " (weight-3, not in Q[pi^2])")
    # spurious-PSLQ artifact: a large-basis "relation" appears at dps=45 but
    # VANISHES at dps=60 -> it was a numerical artifact, so no closed form is
    # folded in.  A real relation is precision-stable; this one changes.
    basis7 = [mp.mpf(1), mp.pi**2, mp.pi**3, mp.zeta(3),
              mp.pi**2 * mp.log(2), mp.pi * mp.log(2)**2, mp.log(2)**3]
    with mp.workdps(45):
        k2_lo = mp.quad(lambda s: PhiK(s)**2 / (2 * s), BRK)
        spur = mp.pslq([k2_lo] + basis7, tol=mp.mpf(10)**-37, maxcoeff=10**7, maxsteps=40000)
    spur60 = mp.pslq([k2] + basis7, tol=mp.mpf(10)**-52, maxcoeff=10**7, maxsteps=40000)
    chk(spur is not None, "expected a spurious PSLQ relation at dps=45 (basis/tol drift?)")
    chk(spur60 is None, "7-basis relation survived at dps=60 (would be a real closed form!)")
    print("  spurious-PSLQ check: dps=45 7-basis =", spur[:3], "... -> dps=60 =", spur60)
    print("      (artifact changes with precision; kappa2 is a genuine weight-3 period)  OK")

    # ---- certified constants used by the section ---------------------------
    pref = (2 / mp.pi) * mp.e**(mp.euler - mp.mpf(3) / 2)
    chk(abs(pref - mp.mpf('0.2529996837')) < mp.mpf('1e-10'), "prefactor")
    chk(abs((1 / (2 * mp.pi)) * mp.e**c0 - pref) < mp.mpf(10)**-40, "(1/2pi)e^{c0}=(2/pi)e^{gamma-3/2}")
    print("  prefactor (2/pi)e^{gamma-3/2} =", mp.nstr(pref, 16), " = (1/2pi)e^{c0}  OK")

    # Theta* : max over s of the Schur bound function Theta(s) (same residue as E4)
    def Theta(s):
        return (mp.log(1 + s**2) + 2 * s * mp.atan(s)) / (2 * (1 + s**2))
    s_star = mp.findroot(lambda s: mp.diff(Theta, s), mp.mpf('1.5'))
    Theta_star = Theta(s_star)
    print("  Theta* = max_s Theta(s) =", mp.nstr(Theta_star, 20), " at s*=", mp.nstr(s_star, 10))
    chk(abs(Theta_star - mp.mpf('0.63492879178810795')) < mp.mpf('1e-15'), "Theta*")
    print("  Theta* matches certified 0.63492879178810795  OK")

    # ---- (optional) kappa_3 / b_3 : demonstrates all-orders computability ---
    with mp.workdps(30):
        def KPhi(s):
            inner = mp.quad(lambda t: Rk(s, t) * PhiK(t) / (2 * t),
                            [0, mp.mpf(1) / 2, 1, 2, 5, 20, mp.inf])
            return inner + d(s) * PhiK(s)
        k3 = mp.quad(lambda s: PhiK(s) * KPhi(s) / (2 * s),
                     [0, mp.mpf(1) / 2, 1, 2, 5, 20, mp.inf])
        k1_30 = mp.quad(lambda s: G(s) * PhiK(s) / (2 * s),
                        [0, mp.mpf(1) / 2, 1, 2, 5, 20, mp.inf])
        b1_30 = 49 * mp.pi**2 / 144 - 3
        C_30 = mp.euler + mp.log(4) - mp.mpf(1) / 2
        b3 = 2 * k3 - 4 * k2 * C_30 + 2 * k1_30 * (C_30**2 - b1_30)
    print("  kappa_3            =", mp.nstr(k3, 12), " (ref 0.0707480)")
    print("  b3 = 2k3-4k2 C+2k1(C^2-b1) =", mp.nstr(b3, 8), " (ref +0.20273; operator fit +0.21)")
    chk(abs(k3 - mp.mpf('0.0707480')) < mp.mpf('5e-6'), "kappa3 vs ref")
    chk(abs(b3 - mp.mpf('0.20273')) < mp.mpf('2e-3'), "b3 vs ref")
    print("  all-orders demonstrator (kappa3,b3) reproduces delivered values  OK")

    print("\nALL %d checks passed." % ck[0])


if __name__ == '__main__':
    main()
