"""rev3 item 3c: recompute X4L suffix logic on the base grid (Na=70) for mutilde=0.33
from scratch (Cholesky c per point, L_w0 and L_base per cell, mu2_pencil_cell), own suffix code."""
import numpy as np, time, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from mpmath import mp
import p4_interchange_uniformity as piu
from p4_interchange_uniformity import alpha_certificate_H, schur_full_seg, ROUND_PAD
from p4_interval_extension import schur_full_seg_w0, alpha_cell_ub_Na
from p4d_coulomb_certify import mu2_pencil_cell
mp.dps = 30
t0 = time.time()
grid = [float(w) for w in piu.wbar_grid()]
cs = []
for w in grid:
    _, ok, c, _, _ = alpha_certificate_H(w, Na=70); assert ok; cs.append(c)
print(f"  {len(grid)} points, Cholesky done ({time.time()-t0:.0f}s)", flush=True)
rows = []
for i in range(len(grid)-1):
    a, b = grid[i], grid[i+1]
    Lw = schur_full_seg_w0(a, b, 1); Lb = schur_full_seg(a, b, 1)
    alw = float(alpha_cell_ub_Na(a, b, cs[i], cs[i+1], Lw, 70)); alb = float(alpha_cell_ub_Na(a, b, cs[i], cs[i+1], Lb, 70))
    d = {}
    for mt in (0.26, 0.33):
        m, k1, k2, k3, case = mu2_pencil_cell(a, b, mt); d[mt] = (float(m), bool(k1 >= 0 and k2 >= 0 and k3 >= 0))
    rows.append((a, b, alw, alb, d))
def al_of(r, which):
    return r[2] if which == "w" else r[3]
def suffix(mt, which):
    ok = [(r[4][mt][0] - al_of(r, which) > 0) and r[4][mt][1] for r in rows]
    i = len(ok)
    while i > 0 and ok[i-1]: i -= 1
    kept = rows[i:]
    ginf = min(r[4][mt][0] - al_of(r, which) for r in kept)
    drop = rows[:i]
    nk = sum(1 for r in drop if not r[4][mt][1]); ng = sum(1 for r in drop if r[4][mt][1] and r[4][mt][0] - al_of(r, which) <= 0)
    return kept[0][0], len(kept), ginf, nk, ng, (drop[-1] if drop else None), sum(ok[:i])
for mt in (0.26, 0.33):
    for which, tag in (("w", "weighted L"), ("b", "base L")):
        wmin, n, ginf, nk, ng, last, npass_before = suffix(mt, which)
        lg = last[4][mt][0] - al_of(last, which)
        print(f"  mutilde={mt} {tag:>10}: wbar_min={wmin:.4f} cells {n}/{len(rows)} inf g_alt={ginf:.6f} | dropped: {nk} k-sign fails, {ng} g<=0 with k ok, {npass_before} passing cells inside dropped prefix; last dropped [{last[0]:.4f},{last[1]:.4f}] g={lg:+.5f} ks={last[4][mt][1]}", flush=True)
print(f"({time.time()-t0:.0f}s)")
