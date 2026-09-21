"""Stage 1: try to FALSIFY the matrix-valued ceiling.
Level: numerical-evidence (float64, finite N).

Two claims under test, for the multiband operator OO:

  (C1) factorized / Einstein:   lam1(OO[Lt, [[.]](varpi)])  <=  h(varpi) * lam_max(Lt)
  (C2) ceiling (general matrix measure):  lam1(OO)  <  g(2) * lam_max(M2),
       M2 := sum_a W_a varpi_a^2   (the matrix second moment, i.e. rho(2 int a2F_ij w dw)/(2piT)^2)

Direction of every one-sided step is recorded: both claims are UPPER bounds on lam1,
so a numerically computed lam1 that EXCEEDS the bound by more than the truncation
drift is a genuine violation (finite-N truncation gives a LOWER estimate of lam1 --
Rayleigh on the truncated subspace -- so a violation at finite N is a violation, not
an artifact; whereas a PASS at finite N is only evidence).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chk_ml_build import (brk, multiband_O, lam_factorized, lam_from_measures,
                           lam1, h_of, g2)

N = 240
G2 = g2(1600)


def scan_factorized():
    print("=" * 78)
    print("(C1) factorized:  lam1(OO) <= h(varpi)*lam_max(Lt) ?   [ratio>1 = VIOLATION]")
    print("=" * 78)
    worst = []
    vps = [0.2, 0.5, 1.0, 2.0, 5.0, 15.0]
    hs = {vp: h_of(vp, N) for vp in vps}
    grid_ab = [0.2, 0.5, 1.0, 2.0, 4.0]
    for w2 in (1.0, 2.0, 5.0):
        w = np.array([1.0, w2])
        for a in grid_ab:
            for b in grid_ab:
                cmax = np.sqrt(a * b)
                for c in np.concatenate([-np.linspace(0.02, 3.0, 14),
                                         np.linspace(0.0, 3.0, 15)]):
                    Lt = np.array([[a, c], [c, b]])
                    lmaxL = np.linalg.eigvalsh(Lt)[-1]
                    if lmaxL <= 0:
                        continue
                    for vp in vps:
                        M = multiband_O(lam_factorized(Lt, vp, N), N, w=w)
                        r = lam1(M) / (hs[vp] * lmaxL)
                        worst.append((r, w2, a, b, c, vp, lam1(M), hs[vp] * lmaxL))
    worst.sort(reverse=True)
    print(f"  {len(worst)} configurations")
    print("   ratio      w2    a     b      c     varpi     lam1      bound")
    for t in worst[:12]:
        print(f"  {t[0]:.6f}  {t[1]:4.1f} {t[2]:4.1f} {t[3]:4.1f} {t[4]:7.3f} {t[5]:6.2f}"
              f"  {t[6]: .6f}  {t[7]: .6f}")
    print("  --- best (most slack) ---")
    for t in worst[-3:]:
        print(f"  {t[0]:.6f}  {t[1]:4.1f} {t[2]:4.1f} {t[3]:4.1f} {t[4]:7.3f} {t[5]:6.2f}")
    # split by sign of c
    neg = [t for t in worst if t[4] < 0]
    pos = [t for t in worst if t[4] >= 0]
    print(f"  max ratio  c>=0 (physical): {max(t[0] for t in pos):.6f}")
    print(f"  max ratio  c< 0 (unphysical repulsive interband): {max(t[0] for t in neg):.6f}")
    return worst


def scan_general():
    print()
    print("=" * 78)
    print("(C2) general matrix measure: lam1(OO) < g(2)*lam_max(M2) ? [ratio>1=VIOLATION]")
    print("=" * 78)
    rng = np.random.default_rng(7)
    worst = []
    for trial in range(3000):
        B = int(rng.integers(2, 4))
        nmode = int(rng.integers(1, 4))
        P = []
        M2 = np.zeros((B, B))
        for _ in range(nmode):
            W = rng.random((B, B)) ** rng.uniform(0.5, 4.0)
            W = (W + W.T) / 2
            W *= rng.uniform(0.05, 3.0)
            vp = float(np.exp(rng.uniform(np.log(0.05), np.log(30.0))))
            P.append((W, vp))
            M2 += W * vp ** 2
        w = np.exp(rng.uniform(-1.5, 1.5, B))
        lam = lam_from_measures(P, N)
        M = multiband_O(lam, N, w=w)
        bnd = G2 * np.linalg.eigvalsh(M2)[-1]
        r = lam1(M) / bnd
        worst.append((r, B, nmode, [f"{v:.3f}" for _, v in P], lam1(M), bnd))
    worst.sort(reverse=True)
    print("   ratio     B  nmode  varpis                lam1        bound")
    for t in worst[:10]:
        print(f"  {t[0]:.6f}  {t[1]}  {t[2]}   {str(t[3]):28s} {t[4]: .5f}  {t[5]: .5f}")
    print(f"  median ratio {np.median([t[0] for t in worst]):.4f}, "
          f"max {worst[0][0]:.6f}")
    return worst


def scan_general_signed():
    print()
    print("=" * 78)
    print("(C2') same, but allowing NEGATIVE interband weights (unphysical for a2F)")
    print("=" * 78)
    rng = np.random.default_rng(11)
    worst = []
    for trial in range(3000):
        B = 2
        nmode = int(rng.integers(1, 3))
        P = []
        M2 = np.zeros((B, B))
        for _ in range(nmode):
            d = rng.random(2) * 2
            off = rng.uniform(-2.0, 0.0)
            W = np.array([[d[0], off], [off, d[1]]])
            vp = float(np.exp(rng.uniform(np.log(0.05), np.log(30.0))))
            P.append((W, vp))
            M2 += W * vp ** 2
        w = np.exp(rng.uniform(-1.0, 1.0, B))
        lam = lam_from_measures(P, N)
        M = multiband_O(lam, N, w=w)
        ev = np.linalg.eigvalsh(M2)
        bnd = G2 * ev[-1]
        if bnd <= 0:
            continue
        r = lam1(M) / bnd
        worst.append((r, [f"{v:.3f}" for _, v in P], P[0][0][0, 1], lam1(M), bnd))
    worst.sort(reverse=True)
    print("   ratio    varpis                  c(mode0)    lam1       bound")
    for t in worst[:10]:
        print(f"  {t[0]:.6f}  {str(t[1]):24s} {t[2]: .4f}  {t[3]: .5f}  {t[4]: .5f}")
    return worst


if __name__ == "__main__":
    print(f"N={N}   g(2)={G2:.9f}")
    w1 = scan_factorized()
    w2 = scan_general()
    w3 = scan_general_signed()
    print()
    print("SUMMARY max ratios:  C1 %.6f   C2 %.6f   C2' %.6f" %
          (w1[0][0], w2[0][0], w3[0][0]))
