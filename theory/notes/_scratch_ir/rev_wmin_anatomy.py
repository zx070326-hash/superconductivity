"""REVIEWER: which side sets wbar_min?  At the first kept cell and the first dropped cell:
mu2 (2x2 cell bound), alpha_cell, nodal alpha, true lam1(M2) worst case, true lam1(M) worst Nc,
true lam2(M) (N=400).  Plus: multi-root detail of the Tc condition at omega_c = 3 omega_E."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rev_common as rc
grid, als, aln = rc.alpha_cells()
for mt, wmin in ((0.1, 0.5391), (0.13, 0.6427), (0.26, 1.1686)):
    i = int(np.argmin(np.abs(grid - wmin)))
    for j, tag in ((i - 1, "first DROPPED"), (i, "first KEPT")):
        a, b = grid[j], grid[j + 1]
        mu2 = min(rc.lam1_2x2(a, mt, c) for c in (0, 1, 2))
        H = rc.H_mat(400, a)
        Es = [np.linalg.eigvalsh(H - mt * rc.Kc(400, Nc))[-2:] for Nc in range(0, 60)]
        k = int(np.argmin([e[1] for e in Es])); l1, l2 = Es[k][1], Es[k][0]
        l2max = max(e[0] for e in Es)
        print(f"mt={mt} {tag:>13} cell [{a:.4f},{b:.4f}]: mu2(2x2)={mu2:.5f}  alpha_cell={als[j]:.5f}  alpha_node={aln[j]:.5f}  "
              f"g={mu2-als[j]:+.5f} | truth@N=400: min_Nc lam1(M)={l1:.5f} (Nc={k}), lam2 there={l2:.2e}, max_Nc lam2={l2max:.2e}")
print("\nmu*=0 reference: alpha_node vs true lam2(H) at small wbar (looseness of the inherited alpha side)")
for w in (0.3, 0.4, 0.5391, 0.6427, 1.0):
    E = np.linalg.eigvalsh(rc.H_mat(400, w)); j = int(np.argmin(np.abs(grid - w)))
    print(f"  w={w:.4f}: alpha_node={aln[j]:.5f}  true lam2(H)={E[-2]:.2e}  ratio={aln[j]/E[-2]:.0f}x")

print("\nMulti-root detail: f(wbar)=lam1(H - mt Kc(floor(r*wbar))) - 1/lambda, r=3, lambda=1, N=400")
N = 400
for mt in (0.1, 0.13):
    ws = np.linspace(0.3, 5.0, 941)
    f = np.array([np.linalg.eigvalsh(rc.H_mat(N, w) - mt * rc.Kc(N, int(np.floor(3 * w))))[-1] for w in ws]) - 1.0
    s = np.sign(f); idx = np.where(s[1:] != s[:-1])[0]
    print(f"  mt={mt}: sign changes at wbar ~ {[f'{ws[i]:.3f}({int(np.floor(3*ws[i]))}->{int(np.floor(3*ws[i+1]))})' for i in idx]}; "
          f"f at 1.9/2.0/2.1: {f[np.argmin(abs(ws-1.9))]:+.4f} {f[np.argmin(abs(ws-2.0))]:+.4f} {f[np.argmin(abs(ws-2.1))]:+.4f}")
    # also lambda=1.2 and r=5
    for lam, r in ((1.0, 5), (1.2, 3), (0.9, 3)):
        f = np.array([np.linalg.eigvalsh(rc.H_mat(N, w) - mt * rc.Kc(N, int(np.floor(r * w))))[-1] for w in ws]) - 1.0 / lam
        s = np.sign(f); idx = np.where(s[1:] != s[:-1])[0]
        print(f"     lam={lam} r={r}: {len(idx)} sign changes at {[f'{ws[i]:.3f}' for i in idx]}")
