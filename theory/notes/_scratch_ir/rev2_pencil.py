"""rev2 item 5: mutilde=0.26 pencil truth (own H, own u).
(a) inf over ALL windows N_c=0..N of the gap of M(w,0.26,N_c) at N=1500 on [5,20]
    via the secular equation of the rank-one downdate, cross-checked vs eigvalsh;
(b) on every retained cell: mu2_cell <= lam1(M^(2)(w,case)) for interior w, 3 cases."""
import numpy as np, pickle, time, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
t0 = time.time()
MT = 0.26; RHO = 2 * MT
def u_vec(N, Nc):
    u = np.zeros(N); k = min(Nc, N); u[:k] = 1 / np.sqrt(2 * np.arange(k) + 1.0); return u
def top2_all_windows(N, w):
    d, Q = np.linalg.eigh(H(N, w)); d = d[::-1]; Q = Q[:, ::-1]
    Upref = np.cumsum(Q * (1 / np.sqrt(2 * np.arange(N) + 1.0))[:, None], axis=0)
    Z2 = np.vstack([np.zeros(N), Upref**2])           # Z2[Nc,i] = (Q^T u(Nc))_i^2
    out = np.empty((N + 1, 2))
    # eigenvalues x of H - rho u u^T solve f(x) = 1 + rho sum_i z_i^2/(d_i - x) = 0;
    # on (d[k+1], d[k]) f is decreasing from +inf to -inf, so f(mid)>0 => root above mid.
    for k in (0, 1):
        lo = np.full(N + 1, d[k + 1]); hi = np.full(N + 1, d[k])
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            f = 1 + RHO * np.sum(Z2 / (d[None, :] - mid[:, None]), axis=1)
            up = f > 0
            lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
        out[:, k] = 0.5 * (lo + hi)
    out[0] = d[:2]
    return out, d
print("== (a) inf over N_c of the true gap of M(w,0.26,N_c), N=1500 ==")
print(f"{'w':>6} {'lam1(H)':>8} {'gap(H)':>8} {'min gap/Nc':>11} {'argmin':>6} {'gap@Nc=1':>9} {'gap@Nc=2':>9} {'gap@Nc=N':>9} {'lam1min':>8} {'xcheck':>8}")
N = 1500
for w in (5.0, 7.0, 10.0, 13.0, 15.0, 16.0, 16.3, 16.5, 17.0, 18.0, 19.0, 20.0):
    out, d = top2_all_windows(N, w)
    gap = out[:, 0] - out[:, 1]; j = int(np.argmin(gap))
    chk = []
    for Nc in (j, int(round(10 * w)), 1, 2):
        u = u_vec(N, Nc); E = np.linalg.eigvalsh(H(N, w) - RHO * np.outer(u, u))
        chk.append(abs((E[-1] - E[-2]) - gap[Nc]))
    print(f"{w:>6.1f} {d[0]:>8.5f} {d[0]-d[1]:>8.5f} {gap.min():>11.6f} {j:>6} {gap[1]:>9.5f} {gap[2]:>9.5f} {gap[N]:>9.5f} {out[:,0].min():>8.5f} {max(chk):>8.1e}")
print("\n== (b) retained cells: mu2_cell <= lam1(M^(2)(w,case)), 11 w per cell, 3 cases ==")
with open(r"D:\physics\superconductivity\theory\notes\_scratch_ir\rev2_cells.pkl", "rb") as f:
    D = pickle.load(f)
rows = D["hi_rows"][0.26]
n_ok = 0
for r in rows:
    if r[6] > 0 and r[7]: n_ok += 1
    else: break
kept = rows[:n_ok]; dropped = rows[n_ok:n_ok + 3]
print(f"  retained prefix: {n_ok} cells, wmax={kept[-1][1]:.4f}")
print("  first dropped: " + "; ".join(f"[{r[0]:.4f},{r[1]:.4f}] g={r[6]:+.5f} ks={r[7]} case={r[8]} alpha={r[4]:.4f} mu2={r[5]:.4f}" for r in dropped))
U2 = {0: np.zeros(2), 1: np.array([1.0, 0.0]), 2: np.array([1.0, 1 / np.sqrt(3)])}
worst = (np.inf, None); nviol = 0; cases = {}
for r in kept:
    a, b = r[0], r[1]; cases[r[8]] = cases.get(r[8], 0) + 1
    for w in np.linspace(a, b, 11):
        for case, u in U2.items():
            m = np.linalg.eigvalsh(H(2, w) - RHO * np.outer(u, u))[-1] - r[5]
            nviol += m < 0
            if m < worst[0]: worst = (m, (a, b, w, case))
print(f"  violations={nviol}; min(lam1(M2)-mu2_cell)={worst[0]:.2e} at [{worst[1][0]:.4f},{worst[1][1]:.4f}] w={worst[1][2]:.4f} case={worst[1][3]}; binding-case histogram={cases}")
a, b = kept[-1][0], kept[-1][1]
print(f"  last kept [{a:.4f},{b:.4f}]: alpha_cell={kept[-1][4]:.5f} mu2={kept[-1][5]:.5f} g={kept[-1][6]:.5f}; " +
      " ".join(f"case{c}:{np.linalg.eigvalsh(H(2,a)-RHO*np.outer(u,u))[-1]:.5f}" for c, u in U2.items()))
print("\n  kept cells: cell g_alt vs true inf_Nc gap (N=1500, midpoint), every 12th + last:")
for r in kept[::12] + [kept[-1]]:
    w = 0.5 * (r[0] + r[1]); out, d = top2_all_windows(N, w); gap = out[:, 0] - out[:, 1]
    print(f"    [{r[0]:.4f},{r[1]:.4f}] g_alt={r[6]:.5f}  true inf_Nc gap={gap.min():.5f}  lam2(H)={d[1]:.5f}  alpha_cell={r[4]:.5f}  true/cert={gap.min()/r[6]:.0f}x")
print(f"({time.time()-t0:.0f}s)")
