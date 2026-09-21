import os, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")
import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

lam = 0.5
spec = Spectrum.einstein(lam, 50.0)
N = 3200
solver = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=N)
t0 = time.time()
tc = solver.tc()
t1 = time.time()
print(f"N={N} tc={tc:.5f} K  ({t1-t0:.1f}s)", flush=True)
t_low = 0.02*tc
wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000, mixing=0.2)
t2 = time.time()
if conv:
    ratio = 2*delta[0]/(KB_MEV_PER_K*tc)
    print(f"N={N} delta[0]={delta[0]:.5f} meV ratio={ratio:.4f} conv={conv} ({t2-t1:.1f}s)", flush=True)
else:
    print("NOT CONVERGED", flush=True)
