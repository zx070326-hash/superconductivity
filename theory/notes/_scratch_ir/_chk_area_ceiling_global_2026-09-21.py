"""Is the Leavens/Carbotte AREA inequality  Tc <= A c(mu*)  (Carbotte 2.50) GLOBAL?

Carbotte's own wording (CARB.txt, Sec. II.H) is a LOCAL statement:
  "We can prove that (2.48) represents a LOCAL maximum for any shape of alpha^2F(w)
   with a fixed A. ... the functional derivative ... is exactly zero at the frequency
   of the base delta function and negative everywhere else."
so the published support is a first-order variational argument with a NUMERICALLY
evaluated functional derivative (his Fig.5 lower frame).  A local maximum does not
imply a global one.  This script looks for a GLOBAL counterexample inside the
two-delta family (the smallest family that leaves the Einstein point), and at the
same time runs the POSITIVE CONTROL for our own moment ceiling on the same samples.

MACHINERY (all built here from the repo's single primitive H_mat):
  The linearized operator O[f] is LINEAR in the kernel f restricted to k>=1
  (H_mat's entries are f(|n-m|), f(n+m+1) and the diagonal -sum_k 2f(k)/(2n+1)).
  Hence for a spectral mixture  f = sum_i w_i [[.]](wbar_i)  we get exactly
      O[f] = sum_i w_i H_mat(N, wbar_i).
  Tc condition:  lam * lam_1( sum_i w_i H_mat(N,wbar_i) - mutilde K_c ) = 1.
  With wbar_i = w_i/(2 pi T):   A = lam <w>/2 = pi T lam <wbar>,  so at the root
      Tc / A = 1 / (pi lam <wbar>*).
  For the single delta this reduces to 2Tc/(lam w_E), the definition of c.

ARROWS TESTED (explicitly, per CLAUDE.md "every arrow must be covered"):
  (a1) O[f] monotone in f (rearrangement identity => the k-th quadratic form q_k>=0)
  (a2) Jensen: f_P(k) <= [[k]](wbar_rms)  pointwise in k                  [analytic]
  (a3) (a1)+(a2) => lam_1(O[f_P]) <= h(wbar_rms)   -- MEASURED here, 0 violations
       expected; this is the moment ceiling's load-bearing step
  (a4) Tc/A <= c(mu*) for two-delta spectra  -- the AREA ceiling, SEARCHED for
       violations (this is the arrow Carbotte only establishes locally)
  (a5) the same with mu*>0 (CE-threshold territory: the Einstein comparison at
       fixed mutilde is NOT unconditional; violations here are EXPECTED above the
       CE threshold and are NOT a counterexample to (2.50))

LEVEL: numerical-evidence (float64; N-convergence reported per regime).
"""
import os
import sys
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DERIV = os.path.normpath(os.path.join(HERE, "..", "..", "derivations"))
sys.path.insert(0, DERIV)
from p4_interchange_uniformity import H_mat  # noqa: E402

TWO_PI = 2.0 * np.pi
C_INF = 0.1827262477890
C0_REF = 0.2308678          # c(0) from _chk_c_mustar_curve_2026-09-21.py (N-converged)
C010_PUB = 0.175            # Carbotte RMP 1990, mu*=0.10


def n_for(wmax):
    return int(min(1100, max(260, 12.0 * wmax + 60)))


def Kc(N, Nc):
    u = np.zeros(N)
    k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return 2.0 * np.outer(u, u)


def lam1_mix(wts, wbars, N, mu_t=0.0, Nc=0):
    A = np.zeros((N, N))
    for w, wb in zip(wts, wbars):
        if w != 0.0:
            A += w * H_mat(N, float(wb))
    if mu_t != 0.0 and Nc > 0:
        A -= mu_t * Kc(N, Nc)
    return float(np.linalg.eigvalsh(A)[-1])


def scale_root(wts, ratios, lam, mu_star=0.0, r=None, N=None, lo=1e-3, hi=200.0,
               nbis=52):
    """smallest s>0 with lam*lam_1(O[f_{s*ratios}] - mutilde Kc)=1, mu*=0 branch exact.

    For mu*=0 the map s -> lam_1 is strictly increasing (monotonicity in f + [[k]]
    increasing in wbar), so plain bisection is valid.  For mu*>0 the Nc staircase
    makes it piecewise; we then bisect on the CONTINUOUS branch with Nc frozen at
    floor(r * s * <ratios>) and report the realized Nc so the reader can see whether
    the freeze was self-consistent.
    """
    ratios = np.asarray(ratios, float)
    gstar = 1.0 / lam
    mu_t = 0.0 if mu_star == 0.0 else mu_star / lam
    a, b = lo, hi

    def g(s):
        wb = s * ratios
        NN = N or n_for(wb.max())
        nc = 0
        if mu_t != 0.0:
            nc = int(np.floor(r * s * float(np.dot(wts, ratios))))
        return lam1_mix(wts, wb, NN, mu_t, nc), nc

    if g(a)[0] >= gstar:
        return a, g(a)[1]
    if g(b)[0] < gstar:
        return None, None
    for _ in range(nbis):
        m = 0.5 * (a + b)
        if g(m)[0] < gstar:
            a = m
        else:
            b = m
    s = 0.5 * (a + b)
    return s, g(s)[1]


def tc_over_A(wts, ratios, lam, mu_star=0.0, r=None, N=None):
    s, nc = scale_root(wts, ratios, lam, mu_star, r, N)
    if s is None:
        return None, None, None
    wbar_mean = s * float(np.dot(wts, ratios))
    return 1.0 / (np.pi * lam * wbar_mean), s, nc


def main():
    t0 = time.time()
    rng = np.random.default_rng(20260921)

    # ---------------------------------------------------------------- (a1) q_k >= 0
    print("=" * 78)
    print("(a1) MONOTONICITY OF O[f] IN f: is the k-th quadratic form q_k >= 0 ?")
    print("=" * 78)
    print("O[f] = sum_{k>=1} f(k) * B_k  with B_k := d O/d f(k) (exact, O is linear).")
    print("q_k(x) = <x,B_k x>.  We test lam_min(B_k) >= 0 for k=1..K at several N.")
    bad = 0
    for N in (40, 80):
        # B_k: entries at |n-m|=k (off-diag), n+m+1=k, minus diagonal 2/(2n+1) for k<=n
        for k in range(1, 13):
            n = np.arange(N)
            inv = 1.0 / np.sqrt(2.0 * n + 1.0)
            nrm = np.outer(inv, inv)
            B = np.zeros((N, N))
            Aij = np.abs(n[:, None] - n[None, :])
            Sij = n[:, None] + n[None, :] + 1
            off = (Aij == k).astype(float)
            np.fill_diagonal(off, 0.0)
            B = (off + (Sij == k).astype(float)) * nrm
            B[np.diag_indices(N)] -= 2.0 * (n >= k) / (2.0 * n + 1.0)
            lmin = float(np.linalg.eigvalsh(B)[0])
            if lmin < -1e-12:
                bad += 1
                if bad <= 6:
                    print(f"  N={N:3d} k={k:2d}: lam_min(B_k) = {lmin:+.6e}  NOT PSD")
    print(f"  => B_k PSD failures: {bad} (of 2*12 = 24 tested)")
    print("  READING: B_k is INDEFINITE, so lam_1(O[.]) is NOT monotone under an")
    print("  ARBITRARY entrywise increase of f -- raising f(k) also raises the mass")
    print("  renormalization, which enters O with a minus sign on the diagonal.")
    print("  This is consistent with, and explains the premises of, the repo's")
    print("  Lemma M (P4_einstein_extremal.md l.824): monotonicity is asserted only")
    print("  for F NON-INCREASING, sum_k F < inf, F not identically 0, G in l^1,")
    print("  F <= G pointwise.  Both [[.]](wbar) and any P-mixture of them ARE")
    print("  non-increasing and summable, so Lemma M applies to arrow (a1) here;")
    print("  what fails is only the shortcut 'B_k PSD'.  We therefore MEASURE (a3).")
    print()

    # ------------------------------------------------- (a3) moment ceiling control
    print("=" * 78)
    print("(a3) POSITIVE CONTROL for the MOMENT ceiling: lam_1(O[f_P]) <= h(wbar_rms)?")
    print("=" * 78)
    nviol = 0
    ntest = 0
    worst = 0.0
    for _ in range(120):
        nat = int(rng.integers(2, 5))
        w = rng.dirichlet(np.ones(nat) * 0.7)
        ratios = np.exp(rng.uniform(np.log(0.15), np.log(6.0), nat))
        s = float(np.exp(rng.uniform(np.log(0.15), np.log(4.0))))
        wb = s * ratios
        NN = n_for(wb.max())
        wrms = float(np.sqrt(np.dot(w, wb**2)))
        a = lam1_mix(w, wb, NN)
        b = lam1_mix([1.0], [wrms], NN)
        ntest += 1
        if a > b + 1e-12:
            nviol += 1
            worst = max(worst, a / b - 1.0)
    print(f"  random mixtures (2-4 deltas, wbar in [0.02,24]): {nviol}/{ntest} "
          f"violations of lam_1(O[f_P]) <= h(wbar_rms)")
    print(f"  worst relative excess: {worst:.3e}")
    print("  (Jensen (a2) is analytic; this covers the operator step (a1)+(a3).)")
    print()

    # ------------------------------------------------- (a4) AREA ceiling, two deltas
    print("=" * 78)
    print("(a4) AREA CEILING GLOBAL SEARCH, two-delta family, mu*=0")
    print("=" * 78)
    print(f"  reference c(0) = {C0_REF:.7f} (Einstein optimum, this repo's solver)")
    print("  Tc/A = 1/(pi lam <wbar>*).  Scanning weight w2 of the HIGH delta and the")
    print("  frequency ratio rho = wbar2/wbar1, maximizing over lam.")
    print()
    best = (0.0, None)
    rows = []
    for rho in (1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 14.0, 25.0):
        for w2 in (0.02, 0.05, 0.10, 0.20, 0.35, 0.50, 0.70, 0.90):
            loc = (0.0, None)
            for lam in (0.7, 0.9, 1.1, 1.3, 1.6, 2.0, 2.6, 3.4, 4.5, 6.0, 8.0,
                        11.0, 15.0, 20.0):
                v, s, _ = tc_over_A([1 - w2, w2], [1.0, rho], lam)
                if v is not None and v > loc[0]:
                    loc = (v, lam)
            rows.append((rho, w2, loc[0], loc[1]))
            if loc[0] > best[0]:
                best = (loc[0], (rho, w2, loc[1]))
    print("   rho    w2     max_lam Tc/A    lam_opt   Tc/A / c(0)")
    for rho, w2, v, lam in rows:
        flag = "  <== EXCEEDS c(0)" if v > C0_REF else ""
        print(f"  {rho:5.1f}  {w2:5.2f}   {v:.6f}      {lam:5.1f}     "
              f"{v/C0_REF:.4f}{flag}")
    print()
    print(f"  GLOBAL MAX over the scanned two-delta family: Tc/A = {best[0]:.6f}")
    print(f"    at (rho, w2, lam) = {best[1]}   ratio to c(0) = {best[0]/C0_REF:.5f}")
    if best[0] > C0_REF * (1 + 1e-6):
        print("  ==> (2.50) VIOLATED in this family at mu*=0.  Would need N-convergence")
        print("      and a second scheme before being reported as a refutation.")
    else:
        print("  ==> NO violation found: every two-delta spectrum scanned obeys")
        print("      Tc <= A c(0).  Supports (2.50) as a GLOBAL statement on this")
        print("      family (does not prove it).")
    print()

    # N-convergence at the two extreme corners of the scan
    print("  N-convergence at the scan corners (Tc/A):")
    for rho, w2, lam in ((25.0, 0.05, 6.0), (25.0, 0.90, 1.3), (1.0, 0.5, 1.15)):
        vals = []
        for mult in (6.0, 12.0, 24.0):
            s0 = None
            NN = int(min(1600, max(260, mult * rho * 2.0)))
            v, s, _ = tc_over_A([1 - w2, w2], [1.0, rho], lam, N=NN)
            vals.append((NN, v))
        txt = "  ".join(f"N={n}:{v:.6f}" for n, v in vals)
        print(f"    rho={rho:5.1f} w2={w2:.2f} lam={lam:4.1f}   {txt}")
    print()

    # --------------------------------------------------- (a5) mu*>0 two-delta check
    print("=" * 78)
    print("(a5) AREA CEILING with mu*=0.10 (r = omega_c/omega_1 fixed at 10)")
    print("=" * 78)
    print(f"  published c(0.10) = {C010_PUB} (Carbotte RMP 1990, cutoff convention")
    print("  declared 'neglected'); see _chk_c_mustar_curve_2026-09-21.py for c(0.1,r).")
    best2 = (0.0, None)
    for rho in (1.0, 2.0, 3.0, 5.0, 8.0, 14.0):
        for w2 in (0.05, 0.15, 0.35, 0.60, 0.90):
            loc = (0.0, None)
            for lam in (0.9, 1.1, 1.3, 1.6, 2.0, 2.6, 3.4, 4.5, 6.0, 8.0):
                v, s, nc = tc_over_A([1 - w2, w2], [1.0, rho], lam, 0.10, 10.0)
                if v is not None and v > loc[0]:
                    loc = (v, lam, nc)
            if loc[0] > best2[0]:
                best2 = (loc[0], (rho, w2) + tuple(loc[1:]))
            print(f"  rho={rho:5.1f} w2={w2:.2f}: max Tc/A = {loc[0]:.6f} "
                  f"at lam={loc[1]}, Nc={loc[2]}")
    print()
    print(f"  GLOBAL MAX (mu*=0.10, r=10): Tc/A = {best2[0]:.6f} at {best2[1]}")
    print(f"  ratio to published 0.175: {best2[0]/C010_PUB:.4f}")
    print()
    print(f"wall {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
