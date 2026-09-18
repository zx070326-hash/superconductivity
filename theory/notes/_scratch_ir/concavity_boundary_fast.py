"""Fast H'' = -2 Psi + R via identity (*): one eigh(A) per (N,x), then O(N^2) per mu~.
h = secular root of 2 mu~ sum_j c_j^2/(mu_j - h) = 1 on (mu_2, mu_1);  v_j = c_j/(h-mu_j)/nu;
R = 2 sum_j a_j^2/(h-mu_j), a = PSI^T y_perp.  Cross-checked against concavity_region_scan.hpp
(which diagonalises M directly) at the start.  numerical-evidence (float64)."""
import os, sys, math, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concavity_exact_resolvent import build

_cache = {}


def prep(N, x):
    key = (N, x)
    if key not in _cache:
        n, u, A, Mp, Opsi = build(N, x)
        mu, PSI = eigh(A)
        _cache[key] = (n, u, Mp, Opsi, mu, PSI)
    return _cache[key]


def hpp_fast(N, x, Nc, mut):
    n, u, Mp, Opsi, mu, PSI = prep(N, x)
    uc = np.where(n < Nc, u, 0.0)
    c = PSI.T @ uc
    m1, m2 = mu[-1], mu[-2]
    if mut == 0.0:
        h = m1
        v = PSI[:, -1]
    else:
        w = (m1 - m2) * 1e-13 + 1e-300
        phi = lambda l: 2 * mut * float(np.sum(c ** 2 / (mu - l))) - 1.0
        h = brentq(phi, m2 + w, m1 - w, xtol=1e-16, rtol=1e-15)
        vj = c / (h - mu)
        vj /= np.linalg.norm(vj)
        v = PSI @ vj
    Psi = float(v @ (Opsi @ v))
    y = Mp @ v
    yp = y - float(v @ y) * v
    a = PSI.T @ yp
    d = h - mu
    if mut == 0.0:
        R = float(np.sum(2.0 * a[:-1] ** 2 / d[:-1]))
    else:
        R = float(np.sum(2.0 * a ** 2 / d))
    return -2.0 * Psi + R, R, Psi


def first_fail(N, x, Nc, muts):
    """returns (mu~_c or None, n_fail, worst H'' (max))"""
    vals = np.array([hpp_fast(N, x, Nc, m)[0] for m in muts])
    bad = np.where(vals > 0)[0]
    return (muts[bad[0]] if bad.size else None), int(bad.size), float(vals.max())


if __name__ == "__main__":
    from concavity_region_scan import hpp
    # 0. cross-check vs direct eigh(M)
    worst = 0.0
    for x, Nc, mut in ((0.09, 6, 0.05), (0.25, 6, 0.13), (1.0, 10, 0.5), (0.49, 7, 0.0), (4.0, 20, 0.02)):
        a = hpp(400, x, Nc, mut)[0]
        b = hpp_fast(400, x, Nc, mut)[0]
        worst = max(worst, abs(a - b) / abs(a))
    print(f"cross-check hpp_fast vs hpp (eigh M): worst rel = {worst:.2e}\n")

    muts = np.logspace(-3, 0, 121)
    N = 400
    print("boundary map, N=400, 121-pt log grid mu~ in [1e-3,1], Nc = floor(r*varpi)")
    print(f"{'varpi':>6} " + " ".join(f"{'r='+str(r):>22}" for r in (3, 5, 10, 20)))
    t0 = time.time()
    rows = {}
    for vp in (0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.6):
        x = vp * vp
        line = f"{vp:>6.2f} "
        for r in (3, 5, 10, 20):
            Nc = max(1, int(math.floor(r * vp)))
            mc, nf, wm = first_fail(N, x, Nc, muts)
            rows[(vp, r)] = (Nc, mc, nf, wm)
            s = "all" if mc is None else f"{mc:.4f}"
            line += f" Nc={Nc:>3} {s:>7} nf={nf:>3}"
        print(line, flush=True)
    print(f"  [{time.time()-t0:.0f}s]\n")

    # N-convergence at boundary-ish rows
    print("N-convergence (Nc, mu~_c, n_fail, max H'') at N=400/800/1600")
    for vp, r in ((0.7, 10), (0.9, 10), (1.0, 10), (1.0, 3), (0.8, 20)):
        x = vp * vp
        Nc = max(1, int(math.floor(r * vp)))
        out = []
        for NN in (400, 800, 1600):
            mc, nf, wm = first_fail(NN, x, Nc, muts)
            out.append(f"N={NN}: {('all' if mc is None else f'{mc:.4f}'):>7} nf={nf:>3} maxH''={wm:+.3e}")
        print(f"varpi={vp} r={r} Nc={Nc} | " + " | ".join(out), flush=True)
