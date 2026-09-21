"""rev2: anatomy of alpha_cell at large wbar (claim-side helpers imported on purpose:
this decomposes the CERTIFICATE, not the truth).  Questions:
 (1) at the mutilde=0.26 cutoff cell [16.3172,16.5057]: which term of alpha_cell
     dominates -- Lipschitz slack L(b-a)/2 or the tails ||B||+||D||?  Would halving
     the cell (same Na) already re-close it?
 (2) at wbar=30 on the SAME geometric grid with Na=420: does the CELL close
     (the docstring quotes only the nodal g_alt=1.17)?"""
import sys, time
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from mpmath import mp
import p4_interval_extension as X
from p4_interchange_uniformity import (alpha_certificate_H, schur_full_seg, schur_coupling,
                                       schur_tailblock, spec_seg, mu2_cell_lb, ROUND_PAD)
from p4d_coulomb_certify import mu2_pencil_cell
mp.dps = 30
def anatomy(a, b, Na, mts=(0.0, 0.13, 0.26)):
    t0 = time.time()
    _, okA, cA, _, _ = alpha_certificate_H(a, Na=Na)
    _, okB, cB, _, _ = alpha_certificate_H(b, Na=Na)
    L = schur_full_seg(a, b, 1)
    sp = spec_seg(a, b, 0)
    nB = float(schur_coupling(Na, b, 0, spec=sp)); nD = float(schur_tailblock(Na, b, 0, spec=sp))
    slack = L * (b - a) / 2
    al = float(max(cA, cB)) + slack + nB + nD + float(ROUND_PAD)
    print(f"  cell [{a:.4f},{b:.4f}] h={b-a:.4f} Na={Na} chol={okA and okB}: c={float(max(cA,cB)):.1e}  "
          f"L={L:.3f} slack=L*h/2={slack:.4f}  ||B||={nB:.4f} ||D||={nD:.4f}  => alpha_cell={al:.4f}   ({time.time()-t0:.0f}s)")
    for mt in mts:
        if mt == 0.0:
            mu2 = float(mu2_cell_lb(a, b)[0]); ks = True
        else:
            mu2, k1, k2, k3, case = mu2_pencil_cell(a, b, mt); mu2 = float(mu2); ks = k1 >= 0 and k2 >= 0 and k3 >= 0
        print(f"      mutilde={mt:<5} mu2={mu2:.4f}  g={mu2-al:+.4f}  ksign={ks}")
    return al
print("(1) mutilde=0.26 cutoff cell and its two halves (Na=280, same as X1):")
a, b = 16.317203, None
grid = X.ext_grid(5.0, 20.0)
i = [j for j in range(len(grid)-1) if abs(grid[j]-16.3172) < 1e-3][0]
a, b = grid[i], grid[i+1]
anatomy(a, b, 280)
m = 0.5 * (a + b)
anatomy(a, m, 280); anatomy(m, b, 280)
print("    -> also the last cell [19.8347,20] at Na=280:")
anatomy(grid[-2], grid[-1], 280)
print("\n(2) wbar=30: geometric cell [30, 30*ratio] with Na=420 (docstring: 'Na=420 still closes at wbar=30, g_alt=1.17' -- nodal only?)")
anatomy(30.0, 30.0 * X.BASE_RATIO, 420, mts=(0.0, 0.13))
# nodal (point) certificate at 30 for comparison
al, ok, c, nB, nD = alpha_certificate_H(30.0, Na=420)
mu2 = float(mu2_cell_lb(30.0, 30.0)[0])
print(f"  nodal w=30 Na=420: alpha={float(al):.4f} (B={float(nB):.4f} D={float(nD):.4f}) mu2={mu2:.4f} g_alt(nodal)={mu2-float(al):+.4f}")
