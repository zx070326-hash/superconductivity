"""REVIEWER shared: alpha cells (mu*=0 machinery, inherited), cached as floats.
Only PUBLIC functions of the certified mu*=0 file are called; nothing is patched."""
import os, sys, time, pickle
import numpy as np
from mpmath import mp
mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "derivations")))
from p4_interchange_uniformity import (H_mat, alpha_certificate_H, wbar_grid,   # noqa
                                       schur_full_seg, alpha_cell_ub)
CACHE = os.path.join(HERE, "rev_alpha_cells.pkl")

def Kc(N, Nc):
    u = np.zeros(N); k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return 2.0 * np.outer(u, u)

def alpha_cells():
    """grid (246 pts), als[i] = certified sup_{w in cell i} lam_2(H) (float),
    alnode[i] = nodal alpha at grid[i]."""
    if os.path.exists(CACHE):
        with open(CACHE, "rb") as f:
            return pickle.load(f)
    t0 = time.time()
    grid = wbar_grid()
    res = [alpha_certificate_H(float(w)) for w in grid]
    cs = [r[2] for r in res]; alnode = [float(r[0]) for r in res]
    als = []
    for i in range(len(grid) - 1):
        a, b = float(grid[i]), float(grid[i + 1])
        als.append(float(alpha_cell_ub(a, b, cs[i], cs[i + 1], schur_full_seg(a, b, 1))))
    out = (np.array(grid, dtype=float), np.array(als), np.array(alnode))
    with open(CACHE, "wb") as f:
        pickle.dump(out, f)
    print(f"[rev_common] alpha cells built in {time.time()-t0:.0f}s")
    return out

# ---- independent 2x2 pencil pieces (NOT imported from the file under test) ----
U2 = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (1.0, 1.0 / np.sqrt(3.0))}

def h2(w):
    b = lambda k: w * w / (k * k + w * w)
    return b(1), (b(1) + b(2)) / np.sqrt(3.0), (b(3) - 2 * b(1)) / 3.0

def M2(w, mt, case):
    h00, h01, h11 = h2(w); u0, u1 = U2[case]
    return np.array([[h00 - 2 * mt * u0 * u0, h01 - 2 * mt * u0 * u1],
                     [h01 - 2 * mt * u0 * u1, h11 - 2 * mt * u1 * u1]])

def lam1_2x2(w, mt, case):
    return float(np.linalg.eigvalsh(M2(w, mt, case))[-1])

def ksigns(w, mt, case):
    e, V = np.linalg.eigh(M2(w, mt, case)); v = V[:, -1]; nn = v @ v
    k1 = (v[0]**2 + 2*v[0]*v[1]/np.sqrt(3) - 2*v[1]**2/3)/nn
    k2 = (2*v[0]*v[1]/np.sqrt(3))/nn
    k3 = (v[1]**2/3)/nn
    return k1, k2, k3, v

def brk(k, w):
    return w * w / (k * k + w * w)

def cell_bound_signfree(a, b, mt, case):
    """Sign-AGNOSTIC rigorous cell bound: v = top eigvec of M2(a); Q(w) = sum k_i [[i]](w) - C;
    each term bounded below by its worse endpoint (each [[i]] is increasing in w)."""
    k1, k2, k3, v = ksigns(a, mt, case); u0, u1 = U2[case]
    C = 2 * mt * (u0 * v[0] + u1 * v[1])**2 / (v @ v)
    q = 0.0
    for k, i in ((k1, 1), (k2, 2), (k3, 3)):
        q += min(k * brk(i, a), k * brk(i, b))
    return q - C, (k1, k2, k3)

def cell_bound_author_style(a, b, mt):
    """Replicates the author's logic independently: min over cases of lam1(M2(a)),
    k-signs of the BEST case only.  Returns (bound, ksig_best_ok, ksig_all_ok, case)."""
    vals = [(lam1_2x2(a, mt, c), c) for c in (0, 1, 2)]
    lam, case = min(vals)
    kb = ksigns(a, mt, case)[:3]
    ok_best = min(kb) >= 0
    ok_all = all(min(ksigns(a, mt, c)[:3]) >= 0 for c in (0, 1, 2))
    return lam, ok_best, ok_all, case
