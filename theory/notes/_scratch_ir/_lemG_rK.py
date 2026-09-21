"""r_K = -lambda_min(O[1_{k<=K}]) and lambda_1, N-convergence scan."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _lemG_reduce import build, ev, step

for K in [1, 2, 3, 4, 5, 8, 12, 20, 40, 80, 160, 320]:
    row = [f"K={K:5d}"]
    for N in [4 * K + 200, 12 * K + 400, 30 * K + 800]:
        B = build(step(K), N)
        l1, lmin = ev(B["O"])
        row.append(f"N={N:6d}: l1={l1:9.6f} r={-lmin:9.6f} rho={-lmin/l1:8.6f}")
    print("  ".join(row))
