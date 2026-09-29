"""Global test of the AREA ceiling  k_B Tc <= c(mu*) A  on the two-delta family,
plus N-convergence audit of the (FA) counterexamples found by the sibling script.

WHY THIS SCRIPT EXISTS.  Two separate questions, both now load-bearing:

(1) Carbotte RMP 62 Sec.II.H proves only that the optimal Einstein spectrum is a
    LOCAL maximum of Tc at fixed A ("We can prove that (2.48) represents a local
    maximum for any shape of alpha^2F with a fixed A").  Leavens 1975 declares his
    mu*!=0 extension "approximate".  So the GLOBAL validity of  Tc <= c(mu*) A  is
    not established in the literature at theorem level.  A single two-delta
    spectrum with  Tc > c(mu*) A  would kill it outright.  DIRECTION: the scan can
    only FALSIFY.  Zero violations is NOT a proof; it upgrades nothing.

(2) The sibling script found (FA) [Einstein extremal at fixed lam AND fixed A]
    violated 34/100.  Per repo discipline a numerical refutation must be shown
    N-converged before it is reported, otherwise it may be a truncation artifact.

CONVENTION.  Same as the sibling: physical w_c, Nc = floor(w_c/(2 pi T)) at each
side's own root; sup rule (smallest root).  c(mu*, r) is RECOMPUTED here from the
same engine, never taken from Leavens' printed table.

LEVEL: numerical-evidence (float64).
"""
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import importlib.util as _ilu  # noqa: E402


def _load(stem):
    spec = _ilu.spec_from_file_location(
        stem.replace("-", "_"), str(HERE / (stem + ".py")))
    mod = _ilu.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


_F = _load("_chk_fixed_lam_area_extremality_2026-09-21")
tc_of, H_from_kernel = _F.tc_of, _F.H_from_kernel
kernel_of_mixture = _F.kernel_of_mixture


def tc_of_N(wts, freqs, lam, mu_star, wc, Nfix, nbis=60):
    """Same root solve but with the truncation N PINNED, for convergence audits."""
    wts = np.asarray(wts, float)
    freqs = np.asarray(freqs, float)
    gstar = 1.0 / lam
    mu_t = 0.0 if mu_star == 0 else mu_star / lam

    def g(T):
        wb = freqs / (2.0 * np.pi * T)
        f = kernel_of_mixture(wts, wb, 2 * Nfix + 2)
        A = H_from_kernel(Nfix, f)
        if mu_t:
            nc = min(int(np.floor(wc / (2.0 * np.pi * T))), Nfix)
            if nc > 0:
                u = np.zeros(Nfix)
                u[:nc] = 1.0 / np.sqrt(2.0 * np.arange(nc) + 1.0)
                A -= mu_t * 2.0 * np.outer(u, u)
        return float(np.linalg.eigvalsh(A)[-1])

    a, b = 1e-4, 400.0
    if g(b) >= gstar or g(a) < gstar:
        return None
    for _ in range(nbis):
        m = 0.5 * (a + b)
        if g(m) >= gstar:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


# c(mu*, r) recomputed here, same engine, Einstein spectrum, max over lam.
def c_of(mu_star, r, lam_grid=None):
    if lam_grid is None:
        lam_grid = np.arange(0.60, 3.001, 0.005)
    best, bl = 0.0, None
    for lam in lam_grid:
        wc = r * 1.0                       # w_E = 1
        t = tc_of([1.0], [1.0], float(lam), mu_star, wc)
        if t is None:
            continue
        v = 2.0 * t / lam                  # A = lam * w_E / 2 = lam/2
        if v > best:
            best, bl = v, float(lam)
    return best, bl


def main():
    t0 = time.time()
    print("=" * 94)
    print("PART 1 -- N-CONVERGENCE AUDIT OF THE (FA) COUNTEREXAMPLES")
    print("=" * 94)
    print("A numerical refutation is only reportable if it survives N -> larger.")
    print(f"{'lam':>5}{'rho':>6}{'w2':>6}{'N':>6}{'Tc(P)':>11}{'Tc_E(<w>)':>12}"
          f"{'FA ratio':>11}")
    for (lam, rho, w2) in ((5.0, 20.0, 0.20), (5.0, 8.0, 0.50),
                           (3.0, 8.0, 0.50), (2.0, 4.0, 0.80)):
        wts = [1 - w2, w2]
        fr = np.array([1.0, rho])
        wm = float(np.dot(wts, fr))
        for N in (300, 600, 1200, 2400):
            tp = tc_of_N(wts, fr, lam, 0.0, None, N)
            te = tc_of_N([1.0], [wm], lam, 0.0, None, N)
            print(f"{lam:5.1f}{rho:6.1f}{w2:6.2f}{N:6d}{tp:11.6f}{te:12.6f}"
                  f"{tp/te:11.6f}")
        print()

    # ---------------------------------------------------------------- PART 2
    print("=" * 94)
    print("PART 2 -- IS THE AREA CEILING  k_B Tc <= c(mu*) A  GLOBAL?")
    print("=" * 94)
    for mu_star, r in ((0.0, None), (0.10, 10.0), (0.13, 10.0)):
        c, lopt = c_of(mu_star, r if r else 10.0)
        print(f"  c({mu_star:.2f}) = {c:.6f}  at lam_opt = {lopt:.3f}"
              "   [recomputed here, not quoted]")
        print(f"{'lam':>5}{'rho':>7}{'w2':>6}{'<w>':>9}{'A':>9}"
              f"{'Tc(P)':>10}{'c*A':>10}{'Tc/(cA)':>10}")
        worst, warg, nviol, ntest = 0.0, None, 0, 0
        for lam in (0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 8.0):
            for rho in (1.0, 1.5, 2.0, 4.0, 8.0, 20.0, 50.0):
                for w2 in (0.02, 0.05, 0.2, 0.5, 0.8, 0.95):
                    if rho == 1.0 and w2 != 0.5:
                        continue
                    wts = [1 - w2, w2]
                    fr = np.array([1.0, rho])
                    wm = float(np.dot(wts, fr))
                    A = 0.5 * lam * wm
                    wc = None if mu_star == 0 else r * wm
                    tp = tc_of(wts, fr, lam, mu_star, wc)
                    if tp is None:
                        continue
                    ntest += 1
                    ratio = tp / (c * A)
                    if ratio > 1 + 1e-9:
                        nviol += 1
                    if ratio > worst:
                        worst, warg = ratio, (lam, rho, w2, wm, A, tp)
                    if ratio > 0.97:
                        print(f"{lam:5.1f}{rho:7.1f}{w2:6.2f}{wm:9.3f}{A:9.3f}"
                              f"{tp:10.5f}{c*A:10.5f}{ratio:10.6f}"
                              + ("  <== AREA CEILING VIOLATED"
                                 if ratio > 1 + 1e-9 else ""))
        print(f"  violations: {nviol}/{ntest}   worst Tc/(cA) = {worst:.6f}"
              f"   at (lam,rho,w2)={warg[:3] if warg else None}")
        print()

    print("DIRECTION AUDIT: this scan can only FALSIFY the area ceiling.")
    print("0 violations leaves it exactly where Carbotte left it (local-max")
    print("argument + numerical table); it does NOT make it a theorem.")
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
