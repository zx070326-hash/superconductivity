"""rev2 item 7 (X3 margin): wbar_c for the tightest IN entry (lambda=0.5, mu*=0.13,
r=3 -> 15.6 vs certified 16.3172) recomputed with the independent builder at
N=600 (reproduce) and N=1500 (truncation sensitivity), plus the alternative
N_c convention floor(r wbar) instead of floor(r wbar + 1/2)."""
import numpy as np, sys, time
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
t0 = time.time()
def u_vec(N, Nc):
    u = np.zeros(N); k = min(Nc, N); u[:k] = 1 / np.sqrt(2 * np.arange(k) + 1.0); return u
lam, mus = 0.5, 0.13; mt = mus / lam
for r in (3, 10):
    for N in (600, 1500):
        for conv in ("floor(rw+1/2)", "floor(rw)"):
            wc = None
            for w in np.arange(12.0, 30.0, 0.1):
                Nc = int(np.floor(r * w + 0.5)) if conv.endswith("2)") else int(np.floor(r * w))
                u = u_vec(N, Nc)
                f = np.linalg.eigvalsh(H(N, w) - 2 * mt * np.outer(u, u))[-1]
                if f >= 1 / lam:
                    wc = w; break
            print(f"  r={r:>2} N={N:>4} Nc={conv:<14}: wbar_c={wc:.2f}  (certified upper end for mutilde=0.26: 16.3172)  {time.time()-t0:.0f}s", flush=True)
