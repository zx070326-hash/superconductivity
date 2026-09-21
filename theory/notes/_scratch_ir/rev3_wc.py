"""rev3 item 6: independent recompute of the four lambda=0.4 anchor entries (self-written H, N=800,
same scan rule as X6: first w on the 0.1 grid from 5 with lam1(H - mt K_c(floor(r w+1/2))) >= 1/lam)."""
import numpy as np, sys, time
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
t0 = time.time(); N = 800; lam = 0.4
def u(N, Nc):
    v = np.zeros(N); k = min(Nc, N); v[:k] = 1/np.sqrt(2*np.arange(k)+1.0); return v
for r, mus, ref in ((10, 0.10, 25.50), (10, 0.13, 39.40), (3, 0.10, 28.90), (3, 0.13, 51.20)):
    mt = mus/lam
    for w in np.arange(5, 80.0, 0.1):
        Nc = int(np.floor(r*w+0.5)); uu = u(N, Nc)
        f = np.linalg.eigvalsh(H(N, w) - 2*mt*np.outer(uu, uu))[-1]
        if f >= 1/lam:
            print(f"r={r} mu*={mus} mt={mt:.3f}: wbar_c={w:.2f} (log {ref})  {'match' if abs(w-ref)<0.051 else 'DIFF'}", flush=True); break
print(f"({time.time()-t0:.0f}s)")
