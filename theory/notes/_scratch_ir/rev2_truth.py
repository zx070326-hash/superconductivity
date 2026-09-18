"""rev2 items 2, 3, 6: truth side with the independent builder (rev2_truth_lib);
claim side = cells dumped by rev2_cells.py."""
import numpy as np, pickle, time, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H, lam12, lam1_2x2
t0 = time.time()
with open(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev2_cells.pkl", "rb") as f:
    D = pickle.load(f)
cells = D["lo"] + D["hi"]
rows0 = D["lo_rows"][0.0] + D["hi_rows"][0.0]
OK, FATAL = "OK", "FATAL"

def containing(w):
    return [(c, r) for c, r in zip(cells, rows0) if c[0] - 1e-12 <= w <= c[1] + 1e-12]

print("== item 6: spot points, N=2000 (19.9 also at N=1000/3000) ==")
print(f"{'w':>6} {'N':>5} {'lam1':>9} {'lam2':>10} {'cell':>18} {'Na':>4} {'alpha_cell':>10} {'mu2':>9} {'alpha-lam2':>11} {'lam1(H2)-mu2':>13}")
fatal = 0
for w in (0.2, 0.25, 7.0, 12.0, 19.9):
    for N in ((1000, 2000, 3000) if w == 19.9 else (2000,)):
        l1, l2 = lam12(N, w)
        for c, r in containing(w):
            ok = (c[4] >= l2) and (r[5] <= lam1_2x2(w))
            fatal += (not ok)
            print(f"{w:>6.2f} {N:>5} {l1:>9.5f} {l2:>10.6f} [{c[0]:.4f},{c[1]:.4f}] {c[2]:>4} {c[4]:>10.5f} {r[5]:>9.5f} {c[4]-l2:>11.5f} {lam1_2x2(w)-r[5]:>13.2e}  {OK if ok else FATAL}")

print("\n== dense sweep: every X0+X1 cell, endpoints+midpoint, N=1000 ==")
pts = {}; worst = (np.inf, None); worst2 = (np.inf, None); nviol = 0
for c, r in zip(cells, rows0):
    a, b = c[0], c[1]
    for w in (a, 0.5 * (a + b), b):
        if w not in pts: pts[w] = lam12(1000, w)
        l1, l2 = pts[w]; m = c[4] - l2; m2 = lam1_2x2(w) - r[5]
        nviol += (m < 0) or (m2 < 0)
        if m < worst[0]: worst = (m, (a, b, w, c[4], l2))
        if m2 < worst2[0]: worst2 = (m2, (a, b, w, r[5]))
print(f"  cells={len(cells)} points={len(pts)} violations={nviol}")
print(f"  min(alpha_cell - lam2) = {worst[0]:.5f} at [{worst[1][0]:.4f},{worst[1][1]:.4f}] w={worst[1][2]:.4f} (alpha={worst[1][3]:.5f}, lam2={worst[1][4]:.6f})")
print(f"  min(lam1(H2) - mu2)   = {worst2[0]:.2e} at [{worst2[1][0]:.4f},{worst2[1][1]:.4f}] w={worst2[1][2]:.4f}")
print("  slack factor alpha_cell/lam2 at midpoints (every 20th hi cell + last):")
for c in D["hi"][::20] + [D["hi"][-1]]:
    w = 0.5 * (c[0] + c[1]); l2 = pts[w][1]
    print(f"    [{c[0]:.3f},{c[1]:.3f}] Na={c[2]} alpha_cell={c[4]:.4f} lam2={l2:.5f} ratio={c[4]/l2:.0f}")

print("\n== item 2: L=schur_full_seg(a,b,1) vs float64 ||H_w(w)|| N=1500 ==")
for c in (D["hi"][-1], D["hi"][60], D["hi"][0], D["lo"][0]):
    a, b = c[0], c[1]
    m = max(np.max(np.abs(np.linalg.eigvalsh(H(1500, w, 1)))) for w in np.linspace(a, b, 5))
    print(f"  [{a:.4f},{b:.4f}] L={c[3]:.4f}  max||H_w||={m:.4f}  ratio={c[3]/m:.2f}  {OK if c[3] >= m else FATAL}")

print("\n== item 3: B2/B4 (ext_full.log) vs float64 ||H^(j)(x)|| N=1500 on J=[0.5w,1.5w] ==")
B = {7.0: (0.2568, 0.6760), 10.0: (0.1677, 0.2465), 15.0: (0.1049, 0.0763), 20.0: (0.0759, 0.0328)}
for w, (B2, B4) in B.items():
    xs = np.linspace(0.5 * w, 1.5 * w, 9)
    m2 = max(np.max(np.abs(np.linalg.eigvalsh(H(1500, x, 2)))) for x in xs)
    m4 = max(np.max(np.abs(np.linalg.eigvalsh(H(1500, x, 4)))) for x in xs)
    print(f"  w={w:>4}: B2={B2:.4f} vs {m2:.4f} ({B2/m2:.2f}x)  B4={B4:.4f} vs {m4:.4f} ({B4/m4:.2f}x)  {OK if B2>=m2 and B4>=m4 else FATAL}")
print(f"\nfatal={fatal}  ({time.time()-t0:.0f}s)")
