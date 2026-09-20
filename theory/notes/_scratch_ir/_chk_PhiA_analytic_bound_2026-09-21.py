"""ANALYTIC upper bound on Phi_E  (rigor-first route; the coordinator asked for a
looser-but-PROVEN Phi in preference to a tighter float64 one).

SETUP.  M(wbar) = K_ph(wbar) - 2 mu_t u u^T,  u = u_c(Nc),  mu_t = mu*/lam.
K_ph symmetric with eigenvalues mu_1 > mu_2 >= ... and eigenvectors psi_j;
c_j := <psi_j, u>.  Write lam1 := lambda_1(M) and D := mu_1 - lam1 >= 0.

LEMMA A (analytic).  Cauchy interlacing for a rank-one PSD downward
perturbation gives  mu_2 <= lam1 <= mu_1,  and Weyl gives
lam1 >= mu_1 - 2 mu_t ||u||^2.

LEMMA B (analytic).  If c_1 != 0 then lam1 is the largest root of the secular
equation  sum_j 2 mu_t c_j^2/(mu_j - lam) = 1  on (mu_2, mu_1).  Splitting off
j=1 and moving the rest to the right:

        2 mu_t c_1^2 / D  =  1 + sum_{j>=2} 2 mu_t c_j^2/(lam1 - mu_j).

THEOREM (Phi-A).  Assume  (PA)  gamma := mu_1 - mu_2  >  2 mu_t ||u||^2.
Then, using Lemma A to bound lam1 - mu_j >= gamma - 2 mu_t ||u||^2 from below,

        S+ := 2 mu_t (||u||^2 - c_1^2) / (gamma - 2 mu_t ||u||^2)

        D  >=  Delta(wbar) := 2 mu_t c_1^2 / (1 + S+)                  (**)

i.e.   lambda_1(M(wbar))  <=  mu_1(wbar) - Delta(wbar).

DIRECTION AUDIT.  (**) is a LOWER bound on the Coulomb-induced eigenvalue
reduction, hence an UPPER bound on lambda_1(M), hence (g is increasing in wbar,
Tc/w_E = 1/(2 pi wbar*)) an UPPER bound on Tc, hence an UPPER bound on Phi_E.
That is the direction a CEILING needs.  A violation of (**) would be a bug.

WHY THIS IS WORTH HAVING.  Every ingredient of Delta -- mu_1, c_1^2, ||u||^2,
gamma -- belongs to the mu*-FREE operator K_ph.  So the mu*-dependence of the
ceiling is reduced to explicit algebra in mu_t on top of a Coulomb-free
eigenproblem that the repo already certifies.  Note ||u_c||^2 = sum_{n<Nc}
1/(2n+1) ~ (1/2) ln Nc: the Morel-Anderson logarithm appears here on its own,
in the DENOMINATOR of the reduction -- i.e. as the mechanism that WEAKENS the
Coulomb suppression as the cutoff grows.

LIMITS.  mu_t -> 0: Delta -> 2 mu_t c_1^2 (linear), so Phi_A -> 1 continuously
and the level-(IV) ceiling collapses to level (III) = C_inf sqrt(lam) wbar2
f_E(lam).  Premise (PA) fails when 2 mu_t ||u||^2 >= gamma; that is the stated
mu* validity interval and it is CHECKED, not assumed, at every grid point.

LEVEL: the inequality (**) is analytic; the numbers below are
numerical-evidence (float64) and would need interval arithmetic on
(mu_1, mu_2, c_1^2, ||u||^2) to become certificate.
"""
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = Path(__file__).resolve().parent
DERIV = (HERE.parents[1] / "derivations").resolve()
sys.path.insert(0, str(DERIV))
sys.path.insert(0, str(HERE))
from p4_interchange_uniformity import H_mat  # noqa: E402

N_DEF = 400


def u_c(N, Nc):
    u = np.zeros(N)
    k = min(int(Nc), N)
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return u


def pieces(wbar, Nc, N=N_DEF):
    A = H_mat(N, float(wbar))
    ev, V = np.linalg.eigh(A)
    mu1, mu2 = ev[-1], ev[-2]
    u = u_c(N, Nc)
    c1 = float(V[:, -1] @ u)
    return mu1, mu2, c1 * c1, float(u @ u), A, u


def delta_bound(wbar, Nc, mu_t, N=N_DEF):
    """Delta of (**), plus the premise flag.  Returns (Delta, premise_ok)."""
    mu1, mu2, c1s, un, _, _ = pieces(wbar, Nc, N)
    gam = mu1 - mu2
    den = gam - 2.0 * mu_t * un
    if den <= 0:
        return 0.0, False
    sp = 2.0 * mu_t * (un - c1s) / den
    return 2.0 * mu_t * c1s / (1.0 + sp), True


def lam1_true(wbar, Nc, mu_t, N=N_DEF):
    A = H_mat(N, float(wbar))
    u = u_c(N, Nc)
    return float(np.linalg.eigvalsh(A - mu_t * 2.0 * np.outer(u, u))[-1])


def root_of(fn, gstar, lo=0.05, hi=8.0, nbis=54):
    """Smallest wbar with fn(wbar) >= gstar  (sup rule; fn increasing in wbar)."""
    if fn(hi) < gstar:
        return None
    if fn(lo) >= gstar:
        return None
    a, b = lo, hi
    for _ in range(nbis):
        m = 0.5 * (a + b)
        if fn(m) >= gstar:
            b = m
        else:
            a = m
    return 0.5 * (a + b)


def main():
    t0 = time.time()
    print("=" * 96)
    print("THEOREM (Phi-A): lambda_1(K_ph - 2 mu_t u u^T) <= mu_1 - Delta,"
          "  Delta = 2 mu_t c_1^2/(1+S+)")
    print("=" * 96)
    print("STEP 1 -- verify the ANALYTIC inequality (**) itself (a violation = bug)")
    print(f"{'wbar':>7}{'Nc':>5}{'mu_t':>8}{'gamma':>9}{'2mt|u|^2':>10}"
          f"{'(PA)':>6}{'D_true':>10}{'Delta':>10}{'D/Delta':>9}")
    nbad = ntot = 0
    worst_ratio = 1e9
    for wbar in (0.4, 0.6, 1.0, 1.4, 2.2, 3.5):
        for Nc in (1, 2, 5, 12, 22, 40):
            for mu_t in (0.01, 0.03, 0.057, 0.10, 0.1351, 0.25):
                mu1, mu2, c1s, un, _, _ = pieces(wbar, Nc)
                gam = mu1 - mu2
                d, ok = delta_bound(wbar, Nc, mu_t)
                dtrue = mu1 - lam1_true(wbar, Nc, mu_t)
                ntot += 1
                if ok:
                    if dtrue < d - 1e-12:
                        nbad += 1
                    worst_ratio = min(worst_ratio, dtrue / d if d > 0 else 1e9)
                if (wbar, Nc) in ((1.0, 12), (2.2, 22), (0.6, 5)) and \
                        mu_t in (0.03, 0.1351, 0.25):
                    print(f"{wbar:7.2f}{Nc:5d}{mu_t:8.4f}{gam:9.4f}"
                          f"{2*mu_t*un:10.4f}{str(ok):>6}{dtrue:10.5f}"
                          f"{d:10.5f}"
                          + (f"{dtrue/d:9.3f}" if d > 0 else f"{'-':>9}"))
    print(f"  (**) violations: {nbad}/{ntot}"
          f"   worst  D_true/Delta = {worst_ratio:.4f}  (must be >= 1)")
    print()

    print("STEP 2 -- premise (PA) on the physically relevant range")
    print("  (PA):  gamma = mu_1 - mu_2  >  2 mu_t ||u_c||^2")
    print(f"{'wbar':>7}{'Nc':>5}{'gamma':>9}{'||u||^2':>9}"
          f"{'mu_t_max = gamma/(2||u||^2)':>28}")
    for wbar, Nc in ((0.4, 4), (0.6, 6), (1.0, 10), (1.39, 13), (2.2, 22),
                     (3.5, 35)):
        mu1, mu2, c1s, un, _, _ = pieces(wbar, Nc)
        gam = mu1 - mu2
        print(f"{wbar:7.2f}{Nc:5d}{gam:9.4f}{un:9.4f}"
              f"{gam/(2*un):28.4f}")
    print()

    print("STEP 3 -- Phi_A (the PROVEN bound) vs Phi_E (float64 truth)")
    print("  Phi_A := wbar*_0 / wbar*_A  with wbar*_A the root of"
          " g_0(wbar) - Delta(wbar) = 1/lam")
    print("  Nc frozen at the value realized at the TRUE root (so both sides")
    print("  use the same physical cutoff); direction: Phi_A >= Phi_E required.")
    print(f"{'lam':>5}{'mu*':>7}{'r':>6}{'Nc':>5}{'Phi_E':>10}{'Phi_A':>10}"
          f"{'Phi_A/Phi_E':>13}{'(PA)':>6}")
    bad = 0
    for lam in (1.0, 1.5, 2.0, 3.0, 5.0):
        for mus in (0.05, 0.10, 0.13):
            for r in (3.0, 10.0):
                mu_t = mus / lam
                gstar = 1.0 / lam
                w0 = root_of(lambda w: lam1_true(w, 0, 0.0), gstar)
                if w0 is None:
                    continue

                def g_true(w):
                    return lam1_true(w, int(np.floor(r * w)), mu_t)

                wt = root_of(g_true, gstar)
                if wt is None:
                    continue
                Nc = int(np.floor(r * wt))

                def g_A(w):
                    d, ok = delta_bound(w, Nc, mu_t)
                    return lam1_true(w, 0, 0.0) - (d if ok else 0.0)

                wa = root_of(g_A, gstar)
                _, ok = delta_bound(wt, Nc, mu_t)
                if wa is None:
                    print(f"{lam:5.1f}{mus:7.2f}{r:6.1f}{Nc:5d}"
                          f"{w0/wt:10.5f}{'no root':>10}")
                    continue
                phiE, phiA = w0 / wt, w0 / wa
                if phiA < phiE - 1e-9:
                    bad += 1
                print(f"{lam:5.1f}{mus:7.2f}{r:6.1f}{Nc:5d}{phiE:10.5f}"
                      f"{phiA:10.5f}{phiA/phiE:13.5f}{str(ok):>6}"
                      + ("  <== DIRECTION BUG" if phiA < phiE - 1e-9 else ""))
    print(f"  direction violations: {bad}  (must be 0)")
    print()

    print("STEP 4 -- mu* -> 0 continuity of Delta (should be linear in mu_t)")
    print(f"{'mu_t':>10}{'Delta':>12}{'Delta/mu_t':>13}{'2 c_1^2':>10}")
    for mu_t in (1e-4, 1e-3, 1e-2, 3e-2, 1e-1):
        d, ok = delta_bound(1.0, 10, mu_t)
        _, _, c1s, _, _, _ = pieces(1.0, 10)
        print(f"{mu_t:10.5f}{d:12.6f}{d/mu_t:13.6f}{2*c1s:10.6f}")
    print("  => Delta -> 2 mu_t c_1^2 as mu_t -> 0, so Phi_A -> 1 CONTINUOUSLY")
    print("     and the level-(IV) ceiling collapses to level (III).")
    print()

    print("STEP 5 -- N convergence of the ingredients (wbar=1.0, Nc=10)")
    print(f"{'N':>6}{'mu_1':>12}{'gamma':>12}{'c_1^2':>12}{'||u||^2':>12}")
    for N in (200, 300, 400, 600, 900):
        mu1, mu2, c1s, un, _, _ = pieces(1.0, 10, N)
        print(f"{N:6d}{mu1:12.8f}{mu1-mu2:12.8f}{c1s:12.8f}{un:12.8f}")
    print()
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
