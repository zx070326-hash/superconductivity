# Probe: does the block-7 interval machinery (Cholesky alpha + segment Schur + exact 2x2 mu_2)
# close cells outside [0.3,5]?  Lower end [LO,0.3], upper end [5,HI].  Same geometric ratio.
import sys, time, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations"))
from p4_interchange_uniformity import (alpha_certificate_H, alpha_cell_ub, schur_full_seg, mu2_cell_lb, rayleigh_mu)
lo, hi = float(sys.argv[1]), float(sys.argv[2])
r = (5.0 / 0.3) ** (1.0 / 245)
grid = [lo]
while grid[-1] * r < hi - 1e-12:
    grid.append(grid[-1] * r)
grid.append(hi)
print(f"[{lo},{hi}] cells: {len(grid)-1}, ratio {r:.5f}", flush=True)
t0 = time.time()
cs = []; fail = []
for w in grid:
    al, ok, c, nB, nD = alpha_certificate_H(w)
    if not ok: fail.append(w); cs.append(None)
    else: cs.append(c)
print(f"Cholesky failures: {len(fail)} {fail[:5]}  ({time.time()-t0:.1f}s)", flush=True)
ginf = np.inf; worst = None; bad = []
rows = []
for i in range(len(grid) - 1):
    a, b = grid[i], grid[i + 1]
    if cs[i] is None or cs[i+1] is None: bad.append((a,b,'chol')); continue
    L = schur_full_seg(a, b, 1)
    al = float(alpha_cell_ub(a, b, cs[i], cs[i + 1], L))
    mu2, k1, k2, k3 = [float(x) for x in mu2_cell_lb(a, b)[:4]]
    g = mu2 - al
    ks = k1 >= 0 and k2 >= 0 and k3 >= 0
    rows.append((a,b,L,mu2,al,g,ks))
    if g <= 0 or not ks: bad.append((a,b,g,ks))
    if g < ginf: ginf, worst = g, rows[-1]
print(f"inf g_alt = {ginf:.6f}  worst cell [{worst[0]:.4f},{worst[1]:.4f}] L={worst[2]:.3f} mu2={worst[3]:.5f} alpha_cell={worst[4]:.5f} ksign={worst[6]}  bad cells={len(bad)}")
for bcell in bad[:8]: print("  BAD", bcell)
# sample rows
for k in np.linspace(0, len(rows)-1, 8).astype(int):
    a,b,L,mu2,al,g,ks = rows[k]
    print(f"  cell [{a:.4f},{b:.4f}] L={L:.3f} mu2={mu2:.5f} alpha_cell={al:.5f} g={g:.5f} ks={ks}")
# nodal values at ends
for w in (lo, hi):
    al, ok, c, nB, nD = alpha_certificate_H(w)
    mu = float(rayleigh_mu(w)); mu2 = float(mu2_cell_lb(w, w)[0])
    print(f"  node w={w}: chol_ok={ok} c={float(c) if c is not None else None}  nB={float(nB) if nB else None} nD={float(nD) if nD else None} alpha={float(al) if al else None}  mu={mu:.5f} mu2={mu2:.5f} delta={mu-float(al) if al else None} g_alt={mu2-float(al) if al else None}")
print(f"{time.time()-t0:.1f}s")
