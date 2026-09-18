"""rev2 item 7: independent recount of the X3 table (IN/OUT of 64 entries) with the
own builder, N=600, bisection on the same 0.02/0.1 grids as X3 (monotone lam1 in wbar
checked by a final forward scan around the bisection result)."""
import numpy as np, sys, time
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from rev2_truth_lib import H
t0 = time.time(); N = 600
ws = np.concatenate([np.arange(0.2, 5, 0.02), np.arange(5, 60.0, 0.1)])
def u_vec(N, Nc):
    u = np.zeros(N); k = min(Nc, N); u[:k] = 1 / np.sqrt(2 * np.arange(k) + 1.0); return u
def lam1(w, mt, r):
    Nc = int(np.floor(r * w + 0.5)); u = u_vec(N, Nc)
    return np.linalg.eigvalsh(H(N, w) - 2 * mt * np.outer(u, u))[-1]
WMIN = {0.0: 0.2, 0.1: 0.5391, 0.13: 0.6427, 0.26: 1.1686}; WMAX = {0.0: 20.0, 0.1: 20.0, 0.13: 20.0, 0.26: 16.3172}
def interval(mt):
    keys = [k for k in WMIN if k >= mt - 1e-12]
    if not keys: return None
    k = min(keys); return WMIN[k], WMAX[k]
now = before = 0; total = 0
for r in (10, 3):
    for lam in (0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0):
        line = []
        for mus in (0.0, 0.05, 0.10, 0.13):
            mt = mus / lam; total += 1
            # bisection on index of ws for the first w with lam1 >= 1/lam (lam1 increasing in w)
            lo, hi = 0, len(ws) - 1
            if lam1(ws[hi], mt, r) < 1 / lam:
                wc = np.nan
            else:
                while hi - lo > 1:
                    mid = (lo + hi) // 2
                    if lam1(ws[mid], mt, r) >= 1 / lam: hi = mid
                    else: lo = mid
                wc = ws[hi]
                # guard against non-monotone N_c steps: forward scan the 6 grid points below
                for j in range(max(0, hi - 6), hi):
                    if lam1(ws[j], mt, r) >= 1 / lam: wc = ws[j]; break
            iv = interval(mt)
            if np.isnan(wc) or iv is None: tag = "--"
            else:
                tag = "IN" if iv[0] <= wc <= iv[1] else "OUT"
                now += tag == "IN"; before += iv[0] <= wc <= 5.0
            line.append(f"{wc:6.2f} {tag:>3}")
        print(f"  r={r:>2} lam={lam:.1f} | " + " | ".join(line), flush=True)
print(f"IN now={now}, IN with old upper end 5={before}, total={total}   ({time.time()-t0:.0f}s)")
