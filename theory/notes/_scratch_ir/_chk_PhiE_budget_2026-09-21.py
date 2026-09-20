"""BUDGET FIRST (CLAUDE.md): what does Phi_E buy on the H3S margin column?

Phi_E(lam, mu*, r) := Tc_E(lam, mu*, r) / Tc_E(lam, 0),   r := omega_c/omega_E.

Consumption point (declared up front): the "material check margin" column
    margin = Tc_ceiling / Tc_expt,
    Tc_ceiling(1-line)      = C_inf sqrt(lam) wbar2
    Tc_ceiling(lam-resolved)= C_inf sqrt(lam) wbar2 * f_E(lam)
    Tc_ceiling(mu*-resolved)= C_inf sqrt(lam) wbar2 * f_E(lam) * Phi_E(lam,mu*,r)
H3S: lam=2.293, wbar2=129.31 meV, Tc_expt=203 K.  Current lam-resolved margin 1.73x.

DEGRADE DIRECTION (one-sided rule): Phi_E <= 1 always (K_c PSD => operator order
=> lam_1(M) <= lam_1(K_ph) => wbar* larger => Tc smaller).  Dropping Phi_E from the
chain (Phi_E -> 1) makes the ceiling LARGER, i.e. makes the PASS criterion EASIER.
So "budget = infinity" here would mean: the ceiling passes anyway without Phi_E.
That is trivially true (a ceiling always "passes" if margin > 1).  The real
consumption point is therefore NOT a pass/fail flip but the MARGIN VALUE itself --
so the honest budget question is: how many of the 1.73x does Phi_E remove, and does
the remainder (shape deficit) still dominate?  This script answers that in float64.

LEVEL: numerical-evidence (float64 point evaluation, N=200/400 truncation).
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DERIV = os.path.normpath(os.path.join(HERE, "..", "..", "derivations"))
sys.path.insert(0, DERIV)

from p4_interchange_uniformity import H_mat  # noqa: E402

TWO_PI = 2.0 * np.pi
C_INF = 0.1827262477890
MEV_PER_K = 8.617333262e-2   # k_B in meV/K


def K_c_np(N, Nc):
    u = np.zeros(N)
    k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return 2.0 * np.outer(u, u)


def g_float(w, mu_t, r, N=200):
    """lam_1(M(w, mutilde, Nc=floor(r w))) in float64 (repo pencil convention)."""
    Nc = int(np.floor(r * w + 1e-9))
    A = H_mat(N, w) - mu_t * K_c_np(N, Nc)
    return float(np.linalg.eigvalsh(A)[-1])


def tc_over_wE(lam, mu_star, r, N=200, w_max=60.0, ngrid=None):
    """T_c/omega_E via the SUP rule: smallest root of g(w)=1/lam (largest T).

    Scans w upward on a fine grid from w_min; the first up-crossing of gstar is the
    sup root.  Down-jumps at w=k/r are handled automatically because we take the
    FIRST crossing.  Returns None if no root on the scanned range.
    """
    if mu_star == 0.0:
        mu_t = 0.0
    else:
        mu_t = mu_star / lam
    gstar = 1.0 / lam
    # coarse scan then bisect on the first up-crossing
    if ngrid is None:
        ngrid = 4000
    ws = np.linspace(0.05, w_max, ngrid)
    prev_w, prev_g = None, None
    for w in ws:
        g = g_float(w, mu_t, r, N)
        if prev_g is not None and prev_g < gstar <= g:
            a, b = prev_w, w
            for _ in range(60):
                m = 0.5 * (a + b)
                if g_float(m, mu_t, r, N) < gstar:
                    a = m
                else:
                    b = m
            wstar = 0.5 * (a + b)
            return 1.0 / (TWO_PI * wstar), wstar
        prev_w, prev_g = w, g
    return None, None


def main():
    print("=" * 78)
    print("BUDGET: Phi_E on the H3S margin column  [numerical-evidence, float64]")
    print("=" * 78)
    lam = 2.293
    wbar2 = 129.31          # meV
    tc_expt = 203.0         # K
    f_E = 0.8439            # f_E(2.293), from _chk_einstein_deficit_curve_fE_2026-09-21.py

    ceil1 = C_INF * np.sqrt(lam) * wbar2 / MEV_PER_K
    ceil2 = ceil1 * f_E
    print(f"lam={lam}  wbar2={wbar2} meV  Tc_expt={tc_expt} K")
    print(f"  one-line ceiling    = {ceil1:9.2f} K   margin {ceil1/tc_expt:6.3f}x")
    print(f"  lam-resolved (f_E)  = {ceil2:9.2f} K   margin {ceil2/tc_expt:6.3f}x")
    print()

    # r conventions:  omega_c = 10*wbar2 -> r=10 ;  w_cut=500 meV -> r=500/129.31
    r_list = [("r=3   (repo physical step)", 3.0),
              ("r=3.867 (w_cut=500 meV)", 500.0 / wbar2),
              ("r=10  (omega_c=10 wbar2)", 10.0),
              ("r=20", 20.0)]

    tc0, w0 = tc_over_wE(lam, 0.0, 3.0)      # mu*=0: r irrelevant
    print(f"  Tc_E(lam,0)/omega_E = {tc0:.6f}   (wbar*={w0:.5f})")
    print(f"  cross-check f_E(lam) from this solver = {tc0/(C_INF*np.sqrt(lam)):.4f} "
          f"(task value {f_E})")
    print()
    print(f"{'convention':<28} {'mu*':>5} {'Tc/wE':>10} {'Phi_E':>8} "
          f"{'ceiling K':>10} {'margin':>8}")
    print("-" * 78)
    rows = []
    for name, r in r_list:
        for mus in (0.10, 0.13):
            tc, ws = tc_over_wE(lam, mus, r)
            if tc is None:
                print(f"{name:<28} {mus:5.2f}  NO ROOT on scanned range")
                continue
            phi = tc / tc0
            ceil3 = ceil2 * phi
            print(f"{name:<28} {mus:5.2f} {tc:10.6f} {phi:8.4f} "
                  f"{ceil3:10.2f} {ceil3/tc_expt:7.3f}x")
            rows.append((name, r, mus, phi, ceil3 / tc_expt))
    print()
    print("READING:")
    worst = min(x[4] for x in rows)
    best = max(x[4] for x in rows)
    print(f"  margin after Phi_E spans [{worst:.3f}x, {best:.3f}x] over the four")
    print(f"  omega_c conventions x two mu*.  lam-resolved margin was "
          f"{ceil2/tc_expt:.3f}x.")
    print(f"  => Phi_E removes at most {(1-worst/(ceil2/tc_expt))*100:.1f}% of the margin.")
    if worst > 1.5:
        print("  VERDICT: margin still > 1.5x => SHAPE deficit (not mu*) dominates.")
    else:
        print("  VERDICT: margin drops below 1.5x => mu* is a material part of the gap;")
        print("           the residual is the honest shape-deficit budget.")
    print()
    print("  NOTE (one-sided direction audit): Phi_E<=1 is the operator-order")
    print("  direction (K_c PSD).  A ceiling multiplied by Phi_E<=1 is a TIGHTER")
    print("  ceiling, so any error that OVERSTATES Phi_E keeps the ceiling valid;")
    print("  an error that UNDERSTATES Phi_E could break it.  Hence the certified")
    print("  version must use the UPPER end of the Phi_E bracket.")


if __name__ == "__main__":
    main()
