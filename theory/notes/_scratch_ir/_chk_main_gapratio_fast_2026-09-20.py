# -*- coding: utf-8 -*-
"""精简版：只判两件事。主会话 2026-09-20。"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, r"D:\physics\superconductivity")
import numpy as np
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

BCS = 2*np.pi/np.exp(np.euler_gamma)
print(f"BCS 2D/kTc = {BCS:.6f}", flush=True)

print("\n=== CHECK 1: lambda=0.3, N sweep (mu*=0, wE=50, wcut=500, t=0.02Tc) ===", flush=True)
spec = Spectrum.einstein(0.3, 50.0)
for N in (400, 800, 1600, 3200):
    t0=time.time()
    s = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=N)
    tc = s.tc(); wn,d,z,c = s.solve_gap(0.02*tc, max_iter=4000, mixing=0.2)
    r = 2*d[0]/(KB_MEV_PER_K*tc)
    print(f"  N={N:5d} Tc={tc:9.5f}K  ratio={r:.4f}  conv={c}  [{time.time()-t0:.0f}s]", flush=True)
