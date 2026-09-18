"""rev2 item 5(a) redo: direct eigvalsh of M(w,0.26,N_c)=H-0.52 u u^T at N=1500 over
an explicit window set (small, physical r=3/10, large, full), plus the rank-one
secular sweep over ALL N_c=0..N with the CORRECT sign, cross-checked."""
import numpy as np, time, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
t0 = time.time(); RHO = 0.52; N = 1500
def u_vec(N, Nc):
    u = np.zeros(N); k = min(Nc, N); u[:k] = 1 / np.sqrt(2 * np.arange(k) + 1.0); return u
def all_windows(N, w):
    d, Q = np.linalg.eigh(H(N, w)); d = d[::-1]; Q = Q[:, ::-1]
    Upref = np.cumsum(Q * (1 / np.sqrt(2 * np.arange(N) + 1.0))[:, None], axis=0)
    Z2 = np.vstack([np.zeros(N), Upref**2])
    out = np.empty((N + 1, 2))
    for k in (0, 1):   # root of f(x)=1-rho*sum z_i^2/(d_i-x) in (d[k+1],d[k]); f decreasing there
        lo = np.full(N + 1, d[k + 1]) + 1e-13; hi = np.full(N + 1, d[k]) - 1e-13
        for _ in range(70):
            mid = 0.5 * (lo + hi)
            f = 1 - RHO * np.sum(Z2 / (d[None, :] - mid[:, None]), axis=1)
            up = f > 0
            lo = np.where(up, mid, lo); hi = np.where(up, hi, mid)
        out[:, k] = 0.5 * (lo + hi)
    out[0] = d[:2]
    return out, d
print(f"{'w':>6} {'gap(H)':>8} | {'inf_Nc gap (secular, all Nc)':>28} {'argmin':>6} | direct eigvalsh gaps at Nc = 1, 2, 3, 5, 10, 3w, 10w, 300, 1500 | xcheck")
for w in (5.0, 8.0, 11.0, 13.5, 15.6, 16.3, 16.5, 18.0, 20.0):
    out, d = all_windows(N, w); gap = out[:, 0] - out[:, 1]; j = int(np.argmin(gap))
    Ncs = [1, 2, 3, 5, 10, int(round(3 * w)), int(round(10 * w)), 300, 1500]
    direct = []
    for Nc in Ncs:
        u = u_vec(N, Nc); E = np.linalg.eigvalsh(H(N, w) - RHO * np.outer(u, u)); direct.append(E[-1] - E[-2])
    xc = max(abs(direct[i] - gap[Ncs[i]]) for i in range(len(Ncs)))
    print(f"{w:>6.1f} {d[0]-d[1]:>8.5f} | {gap.min():>28.5f} {j:>6} | " + " ".join(f"{g:.4f}" for g in direct) + f" | {xc:.1e}", flush=True)
print(f"({time.time()-t0:.0f}s)")
