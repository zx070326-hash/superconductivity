"""REVIEWER (audit item ii): (A) k-sign completeness across ALL three window cases on
every retained cell (author checks only the best case); (B) per-cell soundness of the
author's mu2_pencil_cell vs dense 2x2 truth for all cases; (B') discarded prefix: is the
bound actually invalid there, or only unproven?; (C) end-to-end: true gap_N (dense
eigvalsh, all N_c<=N, several N) >= g_alt(cell) on [wbar_min,5], dense near wbar_min."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rev_common as rc
from p4d_coulomb_certify import mu2_pencil_cell   # object under test

t0 = time.time()
grid, als, _ = rc.alpha_cells()
for mt in (0.1, 0.13):
    rows = []
    for i in range(len(grid) - 1):
        a, b = float(grid[i]), float(grid[i + 1])
        mu2, k1, k2, k3, case = mu2_pencil_cell(a, b, mt)
        rows.append((a, b, float(mu2), als[i], float(mu2) - als[i],
                     bool(k1 >= 0 and k2 >= 0 and k3 >= 0), case))
    start = next(i for i in range(len(rows)) if all(r[4] > 0 and r[5] for r in rows[i:]))
    kept = rows[start:]
    print(f"\n=== mutilde={mt}: wbar_min={kept[0][0]:.4f}, {len(kept)} cells kept, inf g={min(r[4] for r in kept):.6f}")
    fd = rows[start - 1]
    binding = "k-sign ONLY" if fd[4] > 0 else "g<=0"
    print(f"  first DROPPED cell [{fd[0]:.4f},{fd[1]:.4f}]: g={fd[4]:+.5f} ksign_best_ok={fd[5]} case={fd[6]}  -> binding filter: {binding}")

    fails = {0: 0, 1: 0, 2: 0}; worst_k = {0: 9, 1: 9, 2: 9}
    for (a, b, mu2, al, g, ks, case) in kept:
        for c in (0, 1, 2):
            k = rc.ksigns(a, mt, c)[:3]
            worst_k[c] = min(worst_k[c], min(k))
            if min(k) < 0: fails[c] += 1
    print(f"  (A) k-sign failures per case on retained cells: {fails}; min k per case: "
          + ", ".join(f"{c}:{v:+.4f}" for c, v in worst_k.items()))
    print(f"      best case chosen by author on retained cells: {sorted(set(r[6] for r in kept))}")

    viol = 0; minslack = 9; where = None
    for (a, b, mu2, al, g, ks, case) in kept:
        for w in np.linspace(a, b, 25):
            for c in (0, 1, 2):
                v = rc.lam1_2x2(w, mt, c)
                if v - mu2 < minslack: minslack, where = v - mu2, (a, w, c)
                if v < mu2 - 1e-13: viol += 1
    print(f"  (B) author's cell bound vs dense 2x2 truth (25 pts x 3 cases x {len(kept)} cells): "
          f"violations={viol}, min slack={minslack:.3e} at cell a={where[0]:.4f} w={where[1]:.4f} case={where[2]}")

    drop = rows[:start]
    bad_bound = 0; bad_cells = []
    for (a, b, mu2, al, g, ks, case) in drop:
        if min(rc.lam1_2x2(w, mt, c) for w in np.linspace(a, b, 9) for c in (0, 1, 2)) < mu2 - 1e-13:
            bad_bound += 1; bad_cells.append(a)
    extra = f"; invalid cells a in [{min(bad_cells):.4f},{max(bad_cells):.4f}]" if bad_cells else ""
    print(f"  (B') discarded prefix: {len(drop)} cells; actually-INVALID bound in {bad_bound} "
          f"(k-sign fail={sum(1 for r in drop if not r[5])}, g<=0={sum(1 for r in drop if r[4] <= 0)}){extra}")

    def cell_g(w):
        for r in kept:
            if r[0] - 1e-12 <= w <= r[1] + 1e-12: return r[4]
        return None
    wmin = kept[0][0]
    ws = list(np.linspace(wmin, wmin * 1.15, 12)) + list(np.geomspace(wmin * 1.15, 5.0, 16))
    Ns = (2, 3, 4, 5, 6, 8, 12, 16, 24, 32, 48, 64, 96, 140, 200)
    viol = 0; minratio = 9; worst = None
    for w in ws:
        g = cell_g(w)
        if g is None: continue
        tg = np.inf; arg = None
        for N in Ns:
            Hm = rc.H_mat(N, w)
            for Nc in range(0, N + 1):
                E = np.linalg.eigvalsh(Hm - mt * rc.Kc(N, Nc))
                gp = E[-1] - E[-2]
                if gp < tg: tg, arg = gp, (N, Nc, E[-1], E[-2])
        if tg < g - 1e-12: viol += 1
        if tg / g < minratio: minratio, worst = tg / g, (w, g, tg, arg)
        if w < wmin * 1.05 or abs(w - 1) < 0.03 or w > 4.9:
            print(f"      w={w:.4f}: true inf gap={tg:.6f} at (N,Nc)=({arg[0]},{arg[1]}) lam1={arg[2]:.5f} lam2={arg[3]:.5f}; cert g_alt={g:.6f}")
    print(f"  (C) end-to-end over {len(ws)} w x {len(Ns)} N x all Nc: violations={viol}; "
          f"min(true/cert)={minratio:.3f} at w={worst[0]:.4f}")
print(f"\ndone ({time.time()-t0:.0f}s)")
