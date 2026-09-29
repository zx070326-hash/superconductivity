"""CERTIFIED Phi_E grid points, via the repo's certified pencil inversion.

Uses `p4_coulomb_global_tc_certified.invert_sup` (certificate level: interval
Rayleigh lower bound + Kato-Temple upper bound on the EXACT operator, inherited
alpha-side gamma, PROVEN free bound / inherited within-window monotonicity for the
lower-window no-root coverage, sup-rule inversion).

Phi_E bracket from two certified wbar* enclosures:
    Tc/omega_E in [1/(2 pi w_hi), 1/(2 pi w_lo)]
    Phi_lo = Tc_lo(mu*)/Tc_hi(0) = w_lo(0)/w_hi(mu*)
    Phi_hi = Tc_hi(mu*)/Tc_lo(0) = w_hi(0)/w_lo(mu*)

DIRECTION AUDIT (one-sided rule): the ceiling is
    Tc <= C_inf sqrt(lam <w^2>) * f_E(lam) * Phi_E,
so only the UPPER end Phi_hi may be quoted in a ceiling.  Phi_lo is reported for
width only, never consumed by a PASS/FAIL.  Likewise the certified product
    f_E(lam)*Phi_E(lam,mu*,r)  =  Tc_E(lam,mu*,r)/(C_inf sqrt(lam) omega_E)
is bracketed above by  Tc_hi(mu*)/(C_inf_lo sqrt(lam) omega_E).

Usage:  python _chk_PhiE_cert_2026-09-21.py [quick|full]
LEVEL: certificate (per row; rows whose status != 'ok' or sup_ok != True are
reported as NOT certified and must not be quoted at certificate level).
"""
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DERIV = os.path.normpath(os.path.join(HERE, "..", "..", "derivations"))
sys.path.insert(0, DERIV)

import p4_coulomb_global_tc_certified as CG  # noqa: E402

TWO_PI = 2.0 * np.pi
# C_inf certified bracket (P4_C_infinity_enclosure / repo headline 0.1827262477890)
C_INF = 0.1827262477890

MODE = sys.argv[1] if len(sys.argv) > 1 else "quick"
if MODE == "full":
    GRID = [(lam, mus, r) for lam in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 5.0)
            for mus in (0.05, 0.10, 0.13, 0.15) for r in (3.0, 5.0, 10.0, 20.0)]
else:
    GRID = [(lam, mus, r) for lam in (1.0, 1.5, 2.0, 2.5, 3.0)
            for mus in (0.10, 0.13) for r in (3.0, 10.0)]


def cert_wstar(lam, mu_star, r, w_max=6.0):
    mu_t = 0.0 if mu_star == 0.0 else mu_star / lam
    return CG.invert_sup(mu_t, 1.0 / lam, r=r, w_min=0.30, w_max=w_max)


def cert_wstar_mu0(lam, w_max=6.0):
    """mu*=0 reference.

    invert_sup cannot certify a root that lies in the k=0 Coulomb window: it
    brackets from w_guess = k_star/r, and for k_star=0 that is wbar=0, where the
    interval machinery divides by A=0 (ZeroDivisionError).  With mu_t=0 the
    operator is K_c-free, so the window index is a pure bookkeeping device and ANY
    r may be used; we pick r=10 so that every root wbar* in [0.36,1.39]
    (lam in [1,5]) sits in a window k>=3, well inside [w_min, w_max].  The
    lower-window coverage then reduces to the PROVEN free bound (M = K_ph exactly,
    h strictly increasing), so it is not a weakening.
    """
    return CG.invert_sup(0.0, 1.0 / lam, r=10.0, w_min=0.30, w_max=w_max)


def main():
    t0 = time.time()
    print("CERTIFIED Phi_E grid   mode=%s   grid points=%d" % (MODE, len(GRID)))
    print("engine: p4_coulomb_global_tc_certified.invert_sup (certificate level)")
    print()
    # mu*=0 references (certified), cached per lambda
    ref = {}
    for lam in sorted({g[0] for g in GRID}):
        d = cert_wstar_mu0(lam)
        ref[lam] = d
        if d['status'] == 'ok':
            print(f"  mu*=0  lam={lam:4.1f}: wbar* in [{d['w_lo']:.6f},{d['w_hi']:.6f}]"
                  f"  Tc/wE in [{1/(TWO_PI*d['w_hi']):.7f},{1/(TWO_PI*d['w_lo']):.7f}]"
                  f"  f_E in [{1/(TWO_PI*d['w_hi'])/(C_INF*np.sqrt(lam)):.5f},"
                  f"{1/(TWO_PI*d['w_lo'])/(C_INF*np.sqrt(lam)):.5f}]")
        else:
            print(f"  mu*=0  lam={lam:4.1f}: status={d['status']} NOT CERTIFIED")
    print()
    print(f"{'lam':>5} {'mu*':>5} {'r':>5} {'Nc':>3} {'status':>18} {'supOK':>6} "
          f"{'Phi_lo':>9} {'Phi_hi':>9} {'relwid':>9} {'roots':>5}")
    print("-" * 92)
    ok_n = 0
    for lam, mus, r in GRID:
        d0 = ref[lam]
        d = cert_wstar(lam, mus, r)
        st, so = d['status'], d.get('sup_ok')
        if st == 'ok' and so is True and d0['status'] == 'ok':
            phi_lo = d0['w_lo'] / d['w_hi']
            phi_hi = d0['w_hi'] / d['w_lo']
            rel = phi_hi / phi_lo - 1.0
            ok_n += 1
            print(f"{lam:5.1f} {mus:5.2f} {r:5.1f} {d['Nc']:3d} {st:>18} {str(so):>6} "
                  f"{phi_lo:9.5f} {phi_hi:9.5f} {rel:9.2e} {d['nroots']:5d}")
        else:
            print(f"{lam:5.1f} {mus:5.2f} {r:5.1f} "
                  f"{(d['Nc'] if d['Nc'] is not None else -1):3d} {st:>18} "
                  f"{str(so):>6} {'---':>9} {'---':>9} {'---':>9} "
                  f"{d['nroots']:5d}")
            if 'lower_cov' in d and so is False:
                print(f"        lower_cov = {d['lower_cov']}")
    print()
    print(f"CERTIFIED rows: {ok_n}/{len(GRID)}   wall {time.time()-t0:.1f}s")
    print()
    print("NOTE: rows with status != 'ok' or sup_ok != True are NOT certified here;")
    print("they stay at numerical-evidence (see _chk_PhiE_table_2026-09-21.py).")
    print()
    print("KNOWN REPORTING DEFECT OF invert_sup AT LARGE r (found here, direction")
    print("audited, NOT a soundness defect):  _crossings() detects only INTERIOR")
    print("up-crossings gL < gstar <= gR.  But the window-boundary jump of g is DOWN")
    print("only while <u_c,v>^2 grows with Nc, and <u_c,v>^2 is NON-monotone in Nc")
    print("(the top eigenvector of M changes sign).  Measured, N-converged to 10")
    print("digits at N=220/400/800:")
    print("    lam=1, mu*=0.13, r=10:  g(Nc=21, wbar=2.2-) = 0.9983828740")
    print("                            g(Nc=22, wbar=2.2 ) = 1.0001528407")
    print("so gstar = 1/lam = 1 is crossed AT the boundary wbar = 2.2 and the scan")
    print("reports nroots=0 -> status 'above_range'.  DIRECTION: 'above_range' makes")
    print("the row NOT CERTIFIED (safe); it must NOT be read as 'no root / no Tc'.")
    print("A root missed this way also cannot yield a FALSE certificate elsewhere,")
    print("because the lower-window coverage test (free bound / within-window")
    print("monotonicity) cannot fire on a window that actually contains a root.")


if __name__ == "__main__":
    main()
