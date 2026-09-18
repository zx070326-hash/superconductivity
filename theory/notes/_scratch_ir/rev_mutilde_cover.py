"""REVIEWER (audit item iv): wbar_min(mutilde) beyond the two nodes, incl. the lambda<1
corner (mu*=0.13 at lambda 0.77..0.5 -> mutilde 0.17..0.26) and mu*=0.1/lambda=3 -> 0.033.
Also: is a certificate at mutilde_max automatically valid for all smaller mutilde
(operator order M(mt) <= M(mt') for mt >= mt')?"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rev_common as rc
from p4d_coulomb_certify import mu2_pencil_cell
t0 = time.time()
grid, als, aln = rc.alpha_cells()
ncell = len(grid) - 1
print(f"{'mutilde':>8} {'wbar_min':>9} {'cells':>6} {'inf g_alt':>10} {'g_alt(w=1)':>11} {'g_alt(w=1.55)':>13} {'g_alt(w=3)':>11}  first dropped")
hist = {}
for mt in (0.0, 0.033, 0.05, 0.1, 0.13, 0.15, 0.17, 0.2, 0.26, 0.3, 0.5, 1.0):
    rows = []
    for i in range(ncell):
        a, b = float(grid[i]), float(grid[i + 1])
        mu2, k1, k2, k3, case = mu2_pencil_cell(a, b, mt)
        rows.append((a, b, float(mu2) - als[i], bool(k1 >= 0 and k2 >= 0 and k3 >= 0)))
    start = next((i for i in range(ncell) if all(r[2] > 0 and r[3] for r in rows[i:])), None)
    if start is None:
        print(f"{mt:8.3f}  NO certified suffix"); continue
    kept = rows[start:]
    def gat(x):
        r = next((r for r in kept if r[0] - 1e-9 <= x < r[1]), None)
        return r[2] if r else float('nan')
    ld = rows[start - 1] if start > 0 else None
    why = "-" if ld is None else f"g={ld[2]:+.4f} ksign_ok={ld[3]}"
    hist[mt] = (kept[0][0], min(r[2] for r in kept))
    print(f"{mt:8.3f} {kept[0][0]:9.4f} {len(kept):6d} {min(r[2] for r in kept):10.6f} {gat(1.0):11.6f} {gat(1.55):13.6f} {gat(3.0):11.6f}  {why}")

ks = sorted(hist)
print("\nMonotonicity of the certified quantities in mutilde:")
print("  wbar_min(mt) nondecreasing:", all(hist[a][0] <= hist[b][0] + 1e-12 for a, b in zip(ks, ks[1:])))
print("  inf g_alt(mt) nonincreasing:", all(hist[a][1] >= hist[b][1] - 1e-12 for a, b in zip(ks, ks[1:])))
viol = 0
for w in np.geomspace(0.3, 5, 60):
    for c in (0, 1, 2):
        prev = None
        for mt in (0.0, 0.033, 0.1, 0.13, 0.2, 0.26):
            v = rc.lam1_2x2(w, mt, c)
            if prev is not None and v > prev + 1e-14: viol += 1
            prev = v
print(f"  lam1(M2(w;mt)) pointwise nonincreasing in mt (60 w x 3 cases x 6 mt): violations={viol}")

print("\nWeyl-in-mutilde (numerical): worst-Nc true gap at N=300 vs mutilde")
for w in (0.6427, 1.0, 1.55, 3.0):
    H = rc.H_mat(300, w); vals = []
    for mt in (0.0, 0.033, 0.1, 0.13, 0.2, 0.26):
        g = min(np.diff(np.linalg.eigvalsh(H - mt * rc.Kc(300, Nc))[-2:])[0] for Nc in range(0, 301))
        vals.append(g)
    mono = all(vals[i] >= vals[i + 1] - 1e-12 for i in range(len(vals) - 1))
    print(f"  w={w}: gaps={['%.4f' % v for v in vals]} monotone_nonincreasing={mono}")
print(f"\ndone ({time.time()-t0:.0f}s)")
