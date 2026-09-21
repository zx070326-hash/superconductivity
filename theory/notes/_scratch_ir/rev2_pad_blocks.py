"""rev2 items 1(c)/4: (i) ||H_float64 - H_exact||_F at Nv=560 (Rayleigh support at the
new nodes) vs ROUND_PAD=1e-12; (ii) direct domination of row_env/col_env and of the
certified ||B||,||D|| (Na=140/280) by float64 block norms at N=3000, w=10/20."""
import sys, time, numpy as np
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
from mpmath import mp
mp.dps = 30
from p4_interchange_uniformity import H_mat, row_env, col_env, schur_coupling, schur_tailblock, spec_point, spec_seg
from rev2_truth_lib import H
t0 = time.time()
print("(i) float64 entry rounding vs ROUND_PAD, Frobenius norm (>= operator norm):")
for N, w in ((560, 20.0), (560, 15.0), (280, 20.0)):
    A = H_mat(N, w); wm = mp.mpf(w)
    b = [mp.mpf(0)] + [wm * wm / (mp.mpf(k) ** 2 + wm * wm) for k in range(1, 2 * N + 1)]
    inv = [1 / mp.sqrt(2 * n + 1) for n in range(N)]
    cum = [mp.mpf(0)]
    for k in range(1, N): cum.append(cum[-1] + 2 * b[k])
    fro2 = mp.mpf(0); mx = mp.mpf(0)
    for n in range(N):
        for m in range(N):
            ex = (b[n + m + 1] + (b[abs(n - m)] if n != m else 0)) * inv[n] * inv[m]
            if n == m: ex -= cum[n] / (2 * n + 1)
            d = abs(mp.mpf(float(A[n, m])) - ex); fro2 += d * d; mx = max(mx, d)
    print(f"  N={N} w={w}: max|entry err|={mp.nstr(mx,3)}  ||A_f-A||_F={mp.nstr(mp.sqrt(fro2),3)}  (pad 1e-12)  {time.time()-t0:.0f}s", flush=True)
print("(ii) envelope domination at the extension head sizes, N=3000 float64 truth:")
for w, Na in ((20.0, 280), (10.0, 140), (20.0, 140)):
    Hb = np.abs(H(3000, w)); rows = Hb.sum(axis=1); cols = Hb[Na:, :].sum(axis=0)
    sp = spec_point(w, 0)
    r_ok = float(row_env(Na, sp)) >= rows[Na:].max()
    c_ok = all(float(col_env(m, Na, sp)) >= cols[m] for m in range(Na))
    Hf = H(3000, w); B = Hf[:Na, Na:]; Dm = Hf[Na:, Na:]
    nB = np.sqrt(np.max(np.linalg.eigvalsh(B @ B.T))); nD = np.max(np.abs(np.linalg.eigvalsh(Dm)))
    cB = float(schur_coupling(Na, w, 0, spec=sp)); cD = float(schur_tailblock(Na, w, 0, spec=sp))
    l2 = np.linalg.eigvalsh(Hf)[-2]
    print(f"  w={w} Na={Na}: row_env>=max row sum n>=Na: {r_ok} ({float(row_env(Na,sp)):.4f} vs {rows[Na:].max():.4f}); col_env>=all cols: {c_ok}; "
          f"||B|| cert {cB:.4f} vs true {nB:.4f} ({cB/nB:.1f}x); ||D|| cert {cD:.4f} vs true {nD:.4f} ({cD/nD:.1f}x); lam2(H)={l2:.5f}, lam1(D)={np.linalg.eigvalsh(Dm)[-1]:.5f}", flush=True)
print(f"({time.time()-t0:.0f}s)")
