"""Direction cross-validation (team-lead constraint #4): the new certified ||C||_schur
bound (c1sq_schur_bound.schur_C_bound_hi, mpmath.iv) must be >= the true ||C|| for
every (varpi, N) tested. A bound BELOW the true value is a bug (same failure mode as
the interlacing mistake in the previous round), not a good result.

True ||C|| here is float64 numerical-evidence: build the (N+1) x (M-N) block with a
large tail truncation M (converged -- doubling M changes the printed digits by <1e-6
relative, checked at varpi=100 the slowest-converging case) and take its numpy SVD
top singular value (= operator 2-norm).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from c1sq_schur_bound import schur_C_bound_hi
from c1sq_lower_bound import build_data

def f(k, w): return w*w/(w*w+k*k)
def u(n): return 1.0/np.sqrt(2.0*n+1.0)

def true_Cnorm(w, N, M):
    n = np.arange(0, N+1)
    m = np.arange(N+1, M+1)
    K = f(np.abs(n[:, None]-m[None, :]), w) + f(n[:, None]+m[None, :]+1.0, w)
    C = u(n)[:, None]*u(m)[None, :]*K
    return np.linalg.norm(C, 2)

grid = [0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0]
bad = []
for w in grid:
    for N in (60, 200):
        M1 = N + 40000
        M2 = N + 80000
        t1 = true_Cnorm(w, N, M1)
        t2 = true_Cnorm(w, N, M2)
        rel = abs(t2-t1)/t2
        data = build_data(w, N)
        cert = float(schur_C_bound_hi(data, N))
        ok = cert >= t2
        tag = "OK " if ok else "*** BUG: bound < true ***"
        print(f"varpi={w:6.1f} N={N:4d} true(M={M1})={t1:.8f} true(M={M2})={t2:.8f} "
              f"conv_rel={rel:.2e}  cert_schur={cert:.8f}  ratio={cert/t2:6.3f}  {tag}")
        if not ok:
            bad.append((w, N))

print()
if bad:
    print(f"FAIL: {len(bad)} cases with bound < true: {bad}")
else:
    print("PASS: certified Schur bound >= true ||C|| in all cases tested (direction confirmed).")
