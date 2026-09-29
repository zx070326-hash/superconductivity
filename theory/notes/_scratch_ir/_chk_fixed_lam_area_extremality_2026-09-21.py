"""CONJECTURE (FA): at FIXED lam AND FIXED area A, is the Einstein spectrum extremal?

    Tc(P)  <=  Tc_E(lam, w_E = <w>_P, mu*, w_c)      with  <w>_P = 2A/lam   ?

WHY THIS IS THE RIGHT QUESTION.  The two ceilings on the table fix different
things:
  * repo Theorem E fixes lam and <w^2>  -> Einstein at w_E = wbar2   (PROVEN:
    kernel-level Jensen in x = wbar^2, which is concave, + Lemma M)
  * Leavens/Carbotte fixes A only       -> Einstein at w_E = d(mu*)A, then
    maximizes over lam, giving c(mu*)A
(FA) fixes BOTH lam and A.  If true it implies Leavens (maximize its RHS over
lam at fixed A) AND is strictly tighter than ours whenever <w> < wbar2, i.e.
for every non-Einstein spectrum (Cauchy-Schwarz).  So (FA) would strictly
dominate both ceilings simultaneously.

WHY IT IS NOT IMPLIED BY THEOREM E.  Theorem E's Jensen step is concavity of
x -> x/(k^2+x) in x = wbar^2, which produces the SECOND moment.  Getting <w>
instead would require concavity in wbar itself; wbar -> wbar^2/(k^2+wbar^2) is
NOT concave in wbar near 0 (it is convex there: second derivative at wbar=0 is
+2/k^2 > 0).  So the same proof route cannot deliver (FA), and the question is
genuinely open at the analytic level.

WHAT THIS SCRIPT DOES.  Searches for a COUNTEREXAMPLE in the two-delta family
(the smallest family that leaves the Einstein point) and evaluates (FA) on the
six real Mishra spectra.  DIRECTION AUDIT: a scan can only FALSIFY.  Zero
violations means "no counterexample in this family", never "proven".  One
violation, if N-converged, kills (FA).

CONVENTION.  Coulomb cutoff w_c is held FIXED IN PHYSICAL UNITS between P and
its Einstein comparison partner (this is the only honest choice: w_c is a
property of the electronic structure, not of the phonons).  In reduced units
that means Nc = floor(w_c/(2 pi T)) evaluated at each side's OWN root.

LEVEL: numerical-evidence (float64).
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
import importlib.util as _ilu  # noqa: E402


def _load(stem):
    spec = _ilu.spec_from_file_location(
        stem.replace("-", "_"), str(HERE / (stem + ".py")))
    mod = _ilu.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


_E = _load("_phiE_engine_2026-09-21")
_R = _load("_chk_real_spectra_two_ceilings_2026-09-21")
C_INF, MEV_PER_K = _E.C_INF, _E.MEV_PER_K


def n_for(wmax):
    return int(min(1400, max(300, 14.0 * wmax + 80)))


def H_from_kernel(N, f):
    """Build O[f] directly from the kernel array f[0..2N] (f(0) never used).

    Structure (the repo's H_mat, re-derived here so a 500-delta mixture costs one
    matrix instead of 500):
        M[n,m] = ( f(|n-m|)*(1-delta_nm) + f(n+m+1) ) / sqrt((2n+1)(2m+1))
        M[n,n] -= ( sum_{k=1..n} 2 f(k) ) / (2n+1)
    Checked against H_mat for single deltas in main().
    """
    n = np.arange(N)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    off = f[D].copy()
    np.fill_diagonal(off, 0.0)
    M = (off + f[S]) * np.outer(inv, inv)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * f[1:N])))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M


def kernel_of_mixture(wts, wbars, K):
    """f_P(k) = sum_i w_i wbar_i^2/(k^2+wbar_i^2) for k=0..K."""
    k = np.arange(K + 1, dtype=float)
    wb = np.asarray(wbars, float)[:, None]
    return (np.asarray(wts, float)[:, None]
            * (wb**2 / (k[None, :] ** 2 + wb**2))).sum(axis=0)


def lam1_mix(wts, wbars, N, mu_t=0.0, Nc=0):
    f = kernel_of_mixture(wts, wbars, 2 * N + 2)
    A = H_from_kernel(N, f)
    if mu_t and Nc > 0:
        u = np.zeros(N)
        k = min(int(Nc), N)
        u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
        A -= mu_t * 2.0 * np.outer(u, u)
    return float(np.linalg.eigvalsh(A)[-1])


def tc_of(wts, freqs, lam, mu_star, wc, lo=1e-4, hi=400.0, nbis=60):
    """Tc in the SAME units as freqs and wc.  Sup rule: smallest root.

    Solves  lam * lam1( sum_i w_i H(N, f_i/(2 pi T)) - (mu*/lam) Kc(Nc) ) = 1
    with Nc = floor(wc/(2 pi T)) -- the PHYSICAL cutoff, identical for P and E.
    """
    wts = np.asarray(wts, float)
    freqs = np.asarray(freqs, float)
    gstar = 1.0 / lam
    mu_t = 0.0 if mu_star == 0 else mu_star / lam

    def g(T):
        wb = freqs / (2.0 * np.pi * T)
        N = n_for(wb.max())
        nc = int(np.floor(wc / (2.0 * np.pi * T))) if mu_t else 0
        return lam1_mix(wts, wb, N, mu_t, nc)

    a, b = lo, hi
    if g(b) >= gstar:
        return None
    if g(a) < gstar:
        return None
    for _ in range(nbis):
        m = 0.5 * (a + b)
        if g(m) >= gstar:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def main():
    t0 = time.time()
    print("=" * 92)
    print("CONJECTURE (FA): Tc(P) <= Tc_E(lam, <w>_P)  at fixed lam AND fixed A")
    print("=" * 92)
    print("DIRECTION: the scan can only FALSIFY.  0 violations != proof.")
    print()

    print("SELF-CHECK of the kernel-built operator against the repo's H_mat:")
    for N, wb in ((60, 0.7), (60, 3.0), (120, 1.3)):
        A1 = H_from_kernel(N, kernel_of_mixture([1.0], [wb], 2 * N + 2))
        A2 = H_mat(N, wb)
        print(f"  N={N:4d} wbar={wb:4.1f}  max|H_from_kernel - H_mat| = "
              f"{np.abs(A1-A2).max():.3e}")
    print()

    # ------------------------------------------------ two-delta counterexample hunt
    for mu_star, wc_mult in ((0.0, None), (0.13, 10.0)):
        tag = "mu*=0" if mu_star == 0 else f"mu*={mu_star}, w_c = {wc_mult}*<w>"
        print(f"--- two-delta scan, {tag} ---")
        print(f"{'lam':>5}{'rho':>7}{'w2':>6}{'<w>':>8}{'wrms':>8}"
              f"{'Tc(P)':>9}{'Tc_E(<w>)':>10}{'Tc_E(wrms)':>11}"
              f"{'FA ratio':>10}{'E ratio':>9}")
        nviol = 0
        ntest = 0
        worst = 0.0
        for lam in (0.8, 1.2, 2.0, 3.0, 5.0):
            for rho in (1.5, 2.0, 4.0, 8.0, 20.0):
                for w2 in (0.05, 0.2, 0.5, 0.8):
                    wts = [1 - w2, w2]
                    fr = np.array([1.0, rho])
                    wm = float(np.dot(wts, fr))
                    wr = float(np.sqrt(np.dot(wts, fr**2)))
                    wc = None if mu_star == 0 else wc_mult * wm
                    tp = tc_of(wts, fr, lam, mu_star, wc)
                    te_m = tc_of([1.0], [wm], lam, mu_star, wc)
                    te_r = tc_of([1.0], [wr], lam, mu_star, wc)
                    if None in (tp, te_m, te_r):
                        continue
                    ntest += 1
                    fa = tp / te_m
                    er = tp / te_r
                    if fa > 1 + 1e-9:
                        nviol += 1
                        worst = max(worst, fa)
                    if fa > 0.995 or nviol and fa > 1:
                        print(f"{lam:5.1f}{rho:7.1f}{w2:6.2f}{wm:8.3f}{wr:8.3f}"
                              f"{tp:9.5f}{te_m:10.5f}{te_r:11.5f}"
                              f"{fa:10.6f}{er:9.6f}"
                              f"{'  <== (FA) VIOLATED' if fa > 1+1e-9 else ''}")
        print(f"  (FA) violations: {nviol}/{ntest}"
              + (f"   worst ratio {worst:.6f}" if nviol else ""))
        print()

    # ---------------------------------------------------------- real Mishra spectra
    print("=" * 92)
    print("(FA) ON THE SIX REAL SPECTRA (full alpha^2F, not a two-delta proxy)")
    print("=" * 92)
    print("Tc(P) here is the FULL linearized-Eliashberg Tc of the real spectrum,")
    print("solved with the same operator; compare with Tc_E(<w>) and Tc_E(wbar2).")
    print(f"{'case':<13}{'lam':>7}{'Tc(P)':>8}{'TcE(<w>)':>9}{'TcE(w2)':>9}"
          f"{'areaceil':>9}{'FA':>8}{'Tc_exp':>8}")
    print("-" * 71)
    MUS, RMULT = 0.13, 10.0
    for name, rel, tcx in _R.FILES:
        p = _R.DATA / rel
        if not p.is_file():
            continue
        w, a = _R.load_a2f(p)
        d = _R.moments(w, a)
        lam, wm, w2, A = d['lam'], d['wmean'], d['w2'], d['A']
        wc = RMULT * w2
        # bin the spectrum onto a coarse delta-mixture (mass-and-lam preserving):
        # weights of the Einstein mixture are dP = (2/lam) a2F/w dw
        m = w > 1e-9
        ww, aa = w[m], a[m]
        dens = aa / ww
        wts = dens / dens.sum()
        keep = wts > 1e-6
        tp = tc_of(wts[keep] / wts[keep].sum(), ww[keep], lam, MUS, wc)
        te_m = tc_of([1.0], [wm], lam, MUS, wc)
        te_r = tc_of([1.0], [w2], lam, MUS, wc)
        cA = 0.16574 * A          # c(0.13, r=10) recomputed in the sibling script
        if tp is None:
            print(f"{name:<13}  root not found")
            continue
        f = 1.0 / MEV_PER_K
        print(f"{name:<13}{lam:7.4f}{tp*f:8.1f}{te_m*f:9.1f}{te_r*f:9.1f}"
              f"{cA*f:9.1f}{tp/te_m:8.4f}"
              + (f"{tcx:8.1f}" if tcx else f"{'n/a':>8}"))
    print()
    print("  NOTE: Tc(P) is computed from the SAME discretized Einstein mixture that")
    print("  defines lam and the moments, so the comparison is internally consistent;")
    print("  it is NOT a claim about the published Tc of these compounds.")
    print()
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
