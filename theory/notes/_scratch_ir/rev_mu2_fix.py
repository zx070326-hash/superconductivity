"""REVIEWER (audit item ii, wbar_min location): three variants of the 2x2 cell bound
  (a) author: min-case lam1(M2(a)), k-signs of the BEST case only     [as shipped]
  (b) author + k-signs required for ALL three cases                    [logically needed]
  (c) sign-agnostic endpoint bound (no k-sign test at all)             [proposed fix]
and the two truth boundaries the 2x2 route could ever reach:
  (t1) dense-w  min_c lam1(M2_c(w)) - alpha_cell > 0
  (t2) true operator gap inf_{N<=140} inf_{Nc} gap_N at wbar around the boundaries."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import rev_common as rc
t0 = time.time()
grid, als, _ = rc.alpha_cells()
ncell = len(grid) - 1

def suffix_start(ok):
    for i in range(len(ok)):
        if all(ok[i:]): return i
    return None

def fmt(s):
    return f"{grid[s]:.4f}" if s is not None else "NONE"

print(f"{'mt':>5} | {'(a) author':>10} {'(b) all-k':>10} {'(c) signfree':>12} | {'(t1) dense 2x2':>14} | first-drop reason under (a)")
res = {}
for mt in (0.05, 0.1, 0.13, 0.17, 0.2, 0.26):
    ok_a, ok_b, ok_c, g_a, g_c = [], [], [], [], []
    for i in range(ncell):
        a, b = float(grid[i]), float(grid[i + 1])
        lam, okb, okall, case = rc.cell_bound_author_style(a, b, mt)
        g = lam - als[i]
        ok_a.append(g > 0 and okb); ok_b.append(g > 0 and okall); g_a.append(g)
        bc = min(rc.cell_bound_signfree(a, b, mt, c)[0] for c in (0, 1, 2))
        gc = bc - als[i]; ok_c.append(gc > 0); g_c.append(gc)
    sa, sb, sc = suffix_start(ok_a), suffix_start(ok_b), suffix_start(ok_c)
    ok_t = []
    for i in range(ncell):
        a, b = float(grid[i]), float(grid[i + 1])
        m = min(rc.lam1_2x2(w, mt, c) for w in np.linspace(a, b, 41) for c in (0, 1, 2))
        ok_t.append(m - als[i] > 0)
    st = suffix_start(ok_t)
    fd = sa - 1
    okb_fd = rc.cell_bound_author_style(float(grid[fd]), float(grid[fd + 1]), mt)[1]
    reason = f"g={g_a[fd]:+.5f}, ksign_best_ok={okb_fd}"
    print(f"{mt:5.2f} | {fmt(sa):>10} {fmt(sb):>10} {fmt(sc):>12} | {fmt(st):>14} | {reason}")
    res[mt] = (sa, sb, sc, st)
    print(f"      inf g_alt: author suffix={min(g_a[sa:]):.6f}; signfree suffix={min(g_c[sc:]):.6f}; "
          f"signfree bound restricted to author suffix={min(g_c[sa:]):.6f}")

print("\n(t2) true inf_(N<=140) inf_(Nc<=N) gap_N(M) around the certified boundaries:")
for mt in (0.1, 0.13, 0.26):
    sa, sb, sc, st = res[mt]
    wa, wc = grid[sa], grid[sc]
    ws = sorted(set([0.3, 0.35, 0.4, 0.45, 0.5, round(wc * 0.97, 4), round(wc, 4), round(wa * 0.97, 4), round(wa, 4), 0.7, 1.0]))
    print(f"  mt={mt} [author wmin={wa:.4f}, signfree wmin={wc:.4f}]")
    for w in ws:
        tg = np.inf; arg = None
        for N in (2, 3, 4, 6, 8, 12, 20, 40, 80, 140):
            Hm = rc.H_mat(N, w)
            for Nc in range(0, N + 1):
                E = np.linalg.eigvalsh(Hm - mt * rc.Kc(N, Nc))
                if E[-1] - E[-2] < tg: tg, arg = E[-1] - E[-2], (N, Nc, E[-1])
        print(f"     w={w:.4f}: gap={tg:.5f} (N={arg[0]},Nc={arg[1]},lam1={arg[2]:.4f})")
print(f"\ndone ({time.time()-t0:.0f}s)")
