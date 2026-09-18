"""rev3 item 4: float64 entry rounding at Nv=1920 (w=60) and Na=960 vs ROUND_PAD=1e-12.
Exact entries via mpmath on a random sample of 4000 entries (max relative error), and
||A||_F in float64 -> Frobenius bound on ||A_f - A||."""
import numpy as np, sys, time
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
from mpmath import mp
mp.dps = 30
t0 = time.time()
def exact(n, m, w):
    w = mp.mpf(w); s = 1/mp.sqrt((2*n+1)*(2*m+1))
    b = lambda k: w*w/(mp.mpf(k)**2+w*w)
    v = b(n+m+1)*s
    if n != m: v += b(abs(n-m))*s
    else: v -= sum(2*b(k) for k in range(1, n+1))/(2*n+1)
    return v
rng = np.random.default_rng(0)
for N, w in ((1920, 60.0), (960, 60.0), (960, 59.883)):
    A = H(N, w); fro = np.linalg.norm(A)
    idx = rng.integers(0, N, size=(3000, 2)); idx = np.vstack([idx, np.array([[k, k] for k in rng.integers(0, N, 300)])])
    maxrel = 0.0; maxabs = 0.0
    for n, m in idx:
        e = exact(int(n), int(m), w); d = abs(mp.mpf(float(A[n, m])) - e)
        maxabs = max(maxabs, float(d)); maxrel = max(maxrel, float(d/abs(e)) if e != 0 else 0.0)
    print(f"N={N} w={w}: ||A||_F={fro:.3f} max|entry|={np.max(np.abs(A)):.3f}; sampled {len(idx)} entries: max abs err={maxabs:.2e} max rel err={maxrel:.2e}; crude ||A_f-A||_F <= maxrel*||A||_F = {maxrel*fro:.2e}  vs ROUND_PAD 1e-12", flush=True)
print(f"({time.time()-t0:.0f}s)")
