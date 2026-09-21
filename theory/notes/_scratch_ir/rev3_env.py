"""rev3 item 1: envelope domination at Na=960/wbar~60, Na=580/35, Na=340/20.1, Na=280/20.
Claim side: row_env/col_env/schur_coupling/schur_tailblock imported (they ARE the claim);
truth side: exact entry sums from the self-written formulas at N=20000 (rows/cols) and
float64 svd/eigvalsh of the B and D blocks at N=3000."""
import numpy as np, time, sys
sys.path.insert(0, r"D:\physics\superconductivity\theory\notes\_scratch_ir")
sys.path.insert(0, r"D:\physics\superconductivity\theory\derivations")
from rev2_truth_lib import H, brk
from mpmath import mp
from p4_interchange_uniformity import row_env, col_env, spec_seg, schur_coupling, schur_tailblock
mp.dps = 30
t0 = time.time()
def entries_row(n, N, w):
    m = np.arange(N); s = 1/np.sqrt(2*n+1.0)/np.sqrt(2*m+1.0)
    d = np.abs(m-n); off = np.where(d==0, 0.0, brk(np.maximum(d,1), w))
    v = (off + brk(n+m+1, w))*s
    v[n] -= np.sum(2*brk(np.arange(1,n+1), w))/(2*n+1.0)
    return v
cells = [(59.883, 60.0, 960), (34.906, 35.107, 580), (20.0, 20.115, 340), (19.8347, 20.0, 280)]
NBIG = 20000; N3 = 3000
for a, b, Na in cells:
    sp = spec_seg(a, b, 0)
    renv = float(row_env(Na, sp)); nB = float(schur_coupling(Na, b, 0, spec=sp)); nD = float(schur_tailblock(Na, b, 0, spec=sp))
    cenv = np.array([float(col_env(m, Na, sp)) for m in range(Na)])
    print(f"cell [{a},{b}] Na={Na}: row_env(Na)={renv:.5f}  max col_env={cenv.max():.5f}  ||B||cert={nB:.5f} ||D||cert={nD:.5f}", flush=True)
    for w in (a, 0.5*(a+b), b):
        rows = list(range(Na, Na+300)) + list(range(Na+300, 6000, 97))
        rs = max(np.sum(np.abs(entries_row(n, NBIG, w))) for n in rows)
        n = np.arange(Na, NBIG); m = np.arange(Na)
        s = (1/np.sqrt(2*n+1.0))[:,None]*(1/np.sqrt(2*m+1.0))[None,:]
        blk = (brk(n[:,None]-m[None,:], w) + brk(n[:,None]+m[None,:]+1, w))*s
        cs = np.sum(np.abs(blk), axis=0)
        worst_c = np.min(cenv - cs)
        A = H(N3, w); B = A[:Na, Na:]; D = A[Na:, Na:]
        nBt = np.linalg.svd(B, compute_uv=False)[0]; ED = np.linalg.eigvalsh(D); nDt = np.max(np.abs(ED))
        l2 = np.linalg.eigvalsh(A)[-2]
        print(f"   w={w:.4f}: max true rowsum(n>=Na)={rs:.5f} {'<=' if rs<=renv else 'VIOL'} row_env | min(col_env-true colsum)={worst_c:+.2e} {'ok' if worst_c>=0 else 'VIOL'} | ||B|| true={nBt:.5f} ({nB/nBt:.1f}x) | ||D|| true={nDt:.5f} lam1(D)={ED[-1]:+.5f} ({nD/nDt:.1f}x) | true lam2(H)={l2:.5f}", flush=True)
print(f"({time.time()-t0:.0f}s)")
