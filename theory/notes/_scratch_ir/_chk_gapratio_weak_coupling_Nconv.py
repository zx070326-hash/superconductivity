"""
Task 1+2 diagnostic for team-lead's gapratio review.

Reproduces gapratio's own pipeline EXACTLY (same Solver class, same
delta[0] usage, same w_cut=500 meV / omega_E=50 meV / mu*=0 / t_low=0.02*Tc
as in gap_ratio_diagnostic.py) at lambda=0.2, 0.3, 0.5, but SWEEPS
n_matsubara to check convergence -- gapratio hard-coded n_matsubara=400
regardless of lambda, and lambda=0.2 has varpi_c~87 (Tc very low), so the
Matsubara grid at t_low=0.02*Tc is far finer than at lambda=0.5. Need to
check whether N=400 is actually converged there.

Strict level: numerical-evidence (float64), independent rebuild.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")

import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

BCS = 3.5275  # 2*pi/e^gamma... actually BCS is 2*pi/exp(gamma) numerically ~3.528; use gapratio's own constant

lam_values = [0.5, 0.3, 0.2]
omega_e = 50.0  # meV
N_list = [400, 800, 1600]

print("="*90)
print("Reproducing gapratio pipeline: delta[0]=Delta(i*omega_0) used as 'Delta(0)'")
print("mu*=0, w_cut=500 meV, omega_E=50 meV, t_low=0.02*Tc, mixing=0.2, max_iter=3000")
print("="*90)

for lam in lam_values:
    spec = Spectrum.einstein(lam, omega_e)
    print(f"\n--- lambda={lam} ---")
    for N in N_list:
        solver = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=N)
        tc = solver.tc()
        t_low = 0.02 * tc
        # Matsubara spacing check
        t_mev = KB_MEV_PER_K * t_low
        spacing = 2*np.pi*t_mev
        w_max_covered = spacing * N
        wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000, mixing=0.2)
        if conv:
            delta0 = delta[0]
            ratio = 2*delta0/(KB_MEV_PER_K*tc)
            print(f"  N={N:5d}  Tc={tc:10.5f} K  w_max_covered={w_max_covered:9.2f} meV "
                  f"(w_cut=500)  delta[0]={delta0:9.5f} meV  ratio={ratio:.4f}  conv={conv}", flush=True)
        else:
            print(f"  N={N:5d}  Tc={tc:10.5f} K  NOT CONVERGED", flush=True)
