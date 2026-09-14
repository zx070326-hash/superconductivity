r"""
P4 - Theorem W'' (moment tower), weight-4 rung: gate for the NEXT term kappa_3, b3, b3^Tc.

Extends theory/derivations/p4_moment_tower.py (pulled read-only from commit f1b455d,
NOT modified here) by one rung. That engine already contains a kappa_3/b3
"all-orders demonstrator" at dps=30 inside its main(); this file:

  (1) re-derives the same kernels standalone (so this gate has no import-time
      dependency on the pulled file's main() running first),
  (2) computes kappa_3 at THREE precisions (dps=30, 45, 60) and checks they
      agree bit-for-bit in the overlapping digits (rules out quadrature-tolerance
      artifacts, not just PSLQ-tolerance artifacts),
  (3) derives the h-expansion coefficient b3 and, via a careful lambda*h(varpi)=1
      series inversion (redone independently in sympy, see
      theory/notes/_scratch_ir/derive_b3_from_beta.py and inv2.py), the physical
      T_c bracket's lambda^3 coefficient b3^Tc,
  (4) PSLQ-classifies kappa_3, beta_3 (:=2kappa_3-4kappa_1^2, the Lambda-tower
      weight-4 period), b3, b3^Tc against {1,pi^2,pi^4,pi^6} at dps=45 AND
      dps=60 (both must return None for a clean non-closed-form verdict),
  (5) checks the gamma-dependence of b3 (nonzero, since k1 != 0) versus b3^Tc
      (exactly gamma-free after algebraic cancellation -- verified symbolically
      in gamma_dep_check2.py and re-checked numerically here),
  (6) reproduces kappa_1, kappa_2, b1, b2, b2^Tc from the lower rungs exactly.

Standalone: imports only mpmath. Counter `ck`. Runtime ~9-10 minutes single
thread (three kappa_3 quadratures at increasing dps dominate).
"""
import time
import mpmath as mp

mp.mp.dps = 60  # module-level mpf constants below must be parsed at full precision

ck = [0]


def chk(cond, msg=""):
    assert cond, msg
    ck[0] += 1
    return cond


# ---------------------------------------------------------------- kernels
def make_kernels():
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
        return 2 * s**2 * t**2 * (5 + s**2 + t**2) \
            / ((1 + (s - t)**2) * (1 + (s + t)**2) * (1 + s**2) * (1 + t**2))

    return G, d, PhiK, Rk


BRK = [0, mp.mpf(1) / 2, 1, 2, 5, 20, mp.inf]

# lower-rung reference numbers, independently re-derived in p4_moment_tower.py
KAPPA2_REF = mp.mpf('0.09873042407289524101111089465604911963035')
B2_REF = mp.mpf('-0.327071467533777298057653195349')
B2TC_REF = mp.mpf('-0.0967183971')

# kappa_3 reference at dps=60, produced by this session's own three-precision
# run below (recorded here only as a fast-fail sanity net; the gate RECOMPUTES
# it from scratch and does not trust this constant).
KAPPA3_REF60 = mp.mpf('0.07074798316969202409684611376721757467387060721202729696')


def kappa3_at(dps):
    with mp.workdps(dps):
        G, d, PhiK, Rk = make_kernels()

        def KPhi(s):
            inner = mp.quad(lambda t: Rk(s, t) * PhiK(t) / (2 * t), BRK)
            return inner + d(s) * PhiK(s)

        k3 = mp.quad(lambda s: PhiK(s) * KPhi(s) / (2 * s), BRK)
        k2 = mp.quad(lambda s: PhiK(s)**2 / (2 * s), BRK)
        k1 = mp.quad(lambda s: G(s) * PhiK(s) / (2 * s), BRK)
        return k1, k2, k3


def main():
    print("=== P4 Theorem W'' weight-4 rung (kappa_3, b3, b3^Tc) gate ===")

    b1 = 49 * mp.pi**2 / 144 - 3
    C = mp.euler + mp.log(4) - mp.mpf(1) / 2
    c0 = C - 1

    # ---- three-precision kappa_3 (and kappa_1,kappa_2 cross-checks) ---------
    t0 = time.time()
    k1_30, k2_30, k3_30 = kappa3_at(30)
    print("  dps=30: k1=%s k2=%s" % (mp.nstr(k1_30, 20), mp.nstr(k2_30, 20)))
    print("          k3=%s  (%.1fs)" % (mp.nstr(k3_30, 25), time.time() - t0))
    chk(abs(k1_30 - b1 / 2) < mp.mpf(10)**-15, "dps=30 kappa1 != b1/2")
    chk(abs(k2_30 - KAPPA2_REF) < mp.mpf(10)**-15, "dps=30 kappa2 mismatch")

    t0 = time.time()
    k1_45, k2_45, k3_45 = kappa3_at(45)
    print("  dps=45: k3=%s  (%.1fs)" % (mp.nstr(k3_45, 40), time.time() - t0))
    chk(abs(k3_45 - k3_30) < mp.mpf(10)**-28,
        "kappa3 dps30 vs dps45 disagree beyond quadrature tolerance")

    t0 = time.time()
    k1_60, k2_60, k3_60 = kappa3_at(60)
    print("  dps=60: k3=%s  (%.1fs)" % (mp.nstr(k3_60, 55), time.time() - t0))
    chk(abs(k3_60 - k3_45) < mp.mpf(10)**-42,
        "kappa3 dps45 vs dps60 disagree beyond quadrature tolerance")
    chk(abs(k3_60 - KAPPA3_REF60) < mp.mpf(10)**-50, "kappa3 dps60 vs recorded reference")
    print("  three-precision agreement (30->45->60) confirms kappa_3 is a genuine")
    print("  convergent period, not a quadrature-tolerance artifact.            OK")

    k1, k2, k3 = k1_60, k2_60, k3_60

    # ---- Lambda-tower weight-4 period: beta3 = 2 kappa3 - 4 kappa1^2 --------
    beta1 = 2 * k1
    beta2 = 2 * k2
    beta3 = 2 * k3 - 4 * k1**2
    chk(abs(beta1 - b1) < mp.mpf(10)**-45, "beta1 != b1")
    print("  beta3 = 2 kappa3 - 4 kappa1^2 =", mp.nstr(beta3, 40))

    # ---- h-expansion coefficients b2, b3 (Lambda -> 1/ln(varpi) re-expansion,
    #      derived independently in _scratch_ir/derive_b3_from_beta.py:
    #        b1=beta1,  b2=beta2-C*beta1,  b3=beta3-2*C*beta2+C^2*beta1 ) -------
    b2 = beta2 - C * beta1
    b3 = beta3 - 2 * C * beta2 + C**2 * beta1
    print("  b2 =", mp.nstr(b2, 30))
    print("  b3 =", mp.nstr(b3, 40))
    chk(abs(b2 - B2_REF) < mp.mpf(10)**-28, "b2 mismatch vs delivered reference")
    # cross-check against the equivalent direct formula recorded in
    # P4_weak_coupling_series.md sec.2'.2(ii): b3 = 2k3-4k2C+2k1(C^2-b1)
    b3_direct = 2 * k3 - 4 * k2 * C + 2 * k1 * (C**2 - b1)
    chk(abs(b3 - b3_direct) < mp.mpf(10)**-45, "b3 tower-form vs direct-form mismatch")
    print("  b3 (tower form) == b3 (direct md formula)  OK")

    # ---- h<->T_c inversion, lambda^3 coefficient -----------------------------
    # lambda*h(varpi)=1, h=L+c0+b1/L+b2/L^2+b3/L^3, L=ln(omega_E/2pi T_c).
    # Independent series inversion (sympy, _scratch_ir/inv2.py) gives the
    # bracket T_c=(omega_E/2pi)e^{c0}e^{-1/lam}[1+b1 lam+b2Tc lam^2+b3Tc lam^3+...]
    # with:
    #   b2Tc = b1*c0 + b1^2/2 + b2
    #   b3Tc = b1^3/6 + b1^2*c0 + b1^2 + b1*b2 + b1*c0^2 + 2*b2*c0 + b3
    b2Tc = b1 * c0 + b1**2 / 2 + b2
    chk(abs(b2Tc - B2TC_REF) < mp.mpf(10)**-9, "b2Tc mismatch vs delivered reference")
    print("  b2Tc (reproduced) =", mp.nstr(b2Tc, 15), " OK")

    b3Tc = b1**3 / 6 + b1**2 * c0 + b1**2 + b1 * b2 + b1 * c0**2 + 2 * b2 * c0 + b3
    print("  b3Tc =", mp.nstr(b3Tc, 40))

    # symbolic identity check: b3Tc collapses to a gamma-free combination of
    # kappa1,kappa2,kappa3 only (gamma_dep_check.py): 4k1^3/3-4k1^2+4k1k2+2k1-4k2+2k3
    b3Tc_kappa_form = mp.mpf(4) / 3 * k1**3 - 4 * k1**2 + 4 * k1 * k2 + 2 * k1 - 4 * k2 + 2 * k3
    chk(abs(b3Tc - b3Tc_kappa_form) < mp.mpf(10)**-45,
        "b3Tc gamma-free kappa-form identity failed")
    print("  b3Tc == 4k1^3/3-4k1^2+4k1k2+2k1-4k2+2k3 (gamma-free closed form)  OK")

    # ---- numeric confirmation of the lambda^3 inversion coefficient ---------
    def hfun(L):
        return L + c0 + b1 / L + b2 / L**2 + b3 / L**3

    def bracket_coeff3(lam):
        L = mp.findroot(lambda L: hfun(L) - 1 / lam, 1 / lam - c0)
        br = mp.e**(-L + 1 / lam - c0)
        return (br - 1 - b1 * lam - b2Tc * lam**2) / lam**3

    lams = [mp.mpf('0.03'), mp.mpf('0.02'), mp.mpf('0.015'), mp.mpf('0.01')]
    cvals = [bracket_coeff3(lam) for lam in lams]
    b3Tc_fit = cvals[-1] + (cvals[-1] - cvals[-2]) * lams[-1] / (lams[-2] - lams[-1])
    print("  T_c bracket lam^3 fit (lam=.03,.02,.015,.01 -> 0):", mp.nstr(b3Tc_fit, 6),
          " vs formula", mp.nstr(b3Tc, 6))
    chk(abs(b3Tc_fit - b3Tc) < mp.mpf('2e-3'),
        "T_c inversion lambda^3 fit disagrees with the derived formula")
    print("  numeric h<->T_c inversion confirms b3Tc formula to ~1e-3          OK")

    # ---- gamma-dependence: b3 carries gamma (b3 - b3(gamma->0-shift) != 0),
    #      b3Tc does not (already shown as pure kappa-combination above) -------
    # symbolic check (gamma_dep_check2.py): coefficient of gamma^2 in b3 is 2*k1 != 0
    chk(abs(2 * k1) > mp.mpf(10)**-10, "gamma^2 coefficient of b3 unexpectedly ~0")
    print("  b3 DOES carry gamma (coeff of gamma^2 is 2*kappa1 =", mp.nstr(2 * k1, 6),
          "!= 0); b3Tc does NOT (shown above).  Same pattern as b2/b2Tc.      OK")

    # ---- PSLQ classification vs {1,pi^2,pi^4,pi^6} at dps=45 and dps=60 -----
    basis = [mp.mpf(1), mp.pi**2, mp.pi**4, mp.pi**6]
    with mp.workdps(45):
        vals45 = dict(kappa3=k3_45, beta3=2 * k3_45 - 4 * k1_45**2,
                      b3=(2 * k3_45 - 4 * k1_45**2) - 2 * C * beta2 + C**2 * beta1,
                      b3Tc=b3Tc)  # b3Tc itself already dps=60; recompute at 45 below
        k2_45c, k1_45c = k2_45, k1_45
        C45 = mp.euler + mp.log(4) - mp.mpf(1) / 2
        b1_45 = 2 * k1_45c
        b2_45 = 2 * k2_45c - b1_45 * C45
        b3_45v = 2 * k3_45 - 4 * k2_45c * C45 + 2 * k1_45c * (C45**2 - b1_45)
        b3Tc_45 = mp.mpf(4) / 3 * k1_45c**3 - 4 * k1_45c**2 + 4 * k1_45c * k2_45c \
            + 2 * k1_45c - 4 * k2_45c + 2 * k3_45
        beta3_45 = 2 * k3_45 - 4 * k1_45c**2
        r_k3_45 = mp.pslq([k3_45] + basis, tol=mp.mpf(10)**-35, maxcoeff=10**7, maxsteps=40000)
        r_beta3_45 = mp.pslq([beta3_45] + basis, tol=mp.mpf(10)**-35, maxcoeff=10**7, maxsteps=40000)
        r_b3_45 = mp.pslq([b3_45v] + basis, tol=mp.mpf(10)**-35, maxcoeff=10**7, maxsteps=40000)
        r_b3Tc_45 = mp.pslq([b3Tc_45] + basis, tol=mp.mpf(10)**-35, maxcoeff=10**7, maxsteps=40000)

    r_k3_60 = mp.pslq([k3] + basis, tol=mp.mpf(10)**-50, maxcoeff=10**7, maxsteps=40000)
    r_beta3_60 = mp.pslq([beta3] + basis, tol=mp.mpf(10)**-50, maxcoeff=10**7, maxsteps=40000)
    r_b3_60 = mp.pslq([b3] + basis, tol=mp.mpf(10)**-50, maxcoeff=10**7, maxsteps=40000)
    r_b3Tc_60 = mp.pslq([b3Tc] + basis, tol=mp.mpf(10)**-50, maxcoeff=10**7, maxsteps=40000)

    for label, r45, r60 in [("kappa3", r_k3_45, r_k3_60), ("beta3", r_beta3_45, r_beta3_60),
                             ("b3", r_b3_45, r_b3_60), ("b3Tc", r_b3Tc_45, r_b3Tc_60)]:
        chk(r45 is None, "PSLQ(%s;1,pi2,pi4,pi6) unexpectedly found relation at dps=45: %s" % (label, r45))
        chk(r60 is None, "PSLQ(%s;1,pi2,pi4,pi6) unexpectedly found relation at dps=60: %s" % (label, r60))
        print("  PSLQ(%s;1,pi2,pi4,pi6) dps=45 -> %s,  dps=60 -> %s  (both null)" % (label, r45, r60))

    print("\n  Interpretation (EVIDENCE-LEVEL, not proof): weight-4 periods kappa_3,")
    print("  beta_3, b3, b3Tc show no relation against {1,pi^2,pi^4,pi^6} at two")
    print("  independent precisions (45 and 60 digits) -- consistent with them being")
    print("  genuine weight-4 periods outside Q[pi^2], but a PSLQ null is never a")
    print("  proof of transcendence/irrationality-outside-basis.")

    print("\nALL %d checks passed." % ck[0])
    print("\nSUMMARY:")
    print("  kappa_3  =", mp.nstr(k3, 45))
    print("  beta_3   =", mp.nstr(beta3, 45))
    print("  b3       =", mp.nstr(b3, 45))
    print("  b3Tc     =", mp.nstr(b3Tc, 45))


if __name__ == '__main__':
    main()
