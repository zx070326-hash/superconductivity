"""
Scratch exploration for effective-rank universality note (2026-09-16).
Independent from-scratch numpy build of the KAY Einstein operator H(varpi),
plus toy alternative kernels to probe generalization/counterexamples.

NOT a certificate script. numerical-evidence level only (float64, finite N).
Not committed (per CLAUDE.md scratch convention).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

def H_einstein(N, w):
    """Exact KAY-II Einstein-kernel operator, built independently from the
    formula in P4_interchange_uniformity.md sec 1 (not imported)."""
    n = np.arange(N)
    bracket = lambda k: w*w/(k*k + w*w)
    # H1 diagonal
    H1diag = np.zeros(N)
    cumsum = 0.0
    ks = np.arange(1, N+1)
    br_k = bracket(ks.astype(float))
    csum = np.cumsum(2*br_k)
    for i in range(N):
        if i == 0:
            H1diag[0] = 0.0
        else:
            H1diag[i] = csum[i-1] / (2*i+1)
    nn, mm = np.meshgrid(n, n, indexing='ij')
    denom = np.sqrt((2*nn+1)*(2*mm+1))
    H2 = np.where(nn != mm, bracket(np.abs(nn-mm).astype(float))/denom, 0.0)
    H3 = bracket((nn+mm+1).astype(float))/denom
    H = -np.diag(H1diag) + H2 + H3
    return H

def gap_stats(H):
    evals = np.linalg.eigvalsh(H)
    evals = np.sort(evals)[::-1]
    lam1, lam2 = evals[0], evals[1]
    gap = lam1 - lam2
    tr = evals.sum()
    tr2 = (evals**2).sum()
    P = tr**2/tr2 if tr2 > 0 else np.nan
    return lam1, lam2, gap, gap/lam1 if lam1 != 0 else np.nan, P

if __name__ == "__main__":
    N = 900
    print("=== Einstein kernel: gap/lambda1 across varpi in [0.2,60] ===")
    print(f"{'varpi':>8} {'lam1':>10} {'lam2':>10} {'gap':>10} {'gap/lam1':>10} {'partic.ratio':>12}")
    varpis = [0.2,0.3,0.5,1,1.5,2,3,5,7,10,15,20,30,45,60]
    rows = []
    for w in varpis:
        H = H_einstein(N, w)
        lam1, lam2, gap, ratio, P = gap_stats(H)
        rows.append((w, lam1, lam2, gap, ratio, P))
        print(f"{w:8.2f} {lam1:10.5f} {lam2:10.5f} {gap:10.5f} {ratio:10.5f} {P:12.5f}")

    # cross-check against analytic lower bound: h >= ln(w) + c0 - O(1/w), c0=0.4635101
    print()
    print("=== cross-check lower bound rho=2B-1 ~ ln(w)+c0-1, vs numeric lam1 ===")
    import math
    c0 = 0.4635101
    for w, lam1, lam2, gap, ratio, P in rows:
        lb = math.log(w) + c0
        print(f"varpi={w:7.2f}  lam1_numeric={lam1:9.5f}  ln(w)+c0={lb:9.5f}  diff={lam1-lb:8.5f}")

    print()
    print("=== gap universal-constant check: lam2 <= Theta* = 0.63492879 (M1 bound) ===")
    Theta_star = 0.63492879
    for w, lam1, lam2, gap, ratio, P in rows:
        print(f"varpi={w:7.2f}  lam2={lam2:9.5f}  <= Theta*={Theta_star}?  {lam2 <= Theta_star}")
