"""
Follow-up: use the solver's OWN auto-N estimator (n_matsubara=None) instead
of gapratio's hard-coded n_matsubara=400. This should automatically scale N
to cover w_cut given T, which is exactly the axis gapratio starved by fixing
N=400 regardless of lambda.

Also computes what N gapratio's fixed N=400 would have needed to actually
cover w_cut at t_low=0.02*Tc for each lambda (spacing = 2*pi*T_meV).

Strict level: numerical-evidence (float64).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, time
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")
import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

omega_e = 50.0
for lam in (0.5, 0.3, 0.2):
    spec = Spectrum.einstein(lam, omega_e)
    solver = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=None,
                     matsubara_oversample=3.0)
    t0 = time.time()
    tc = solver.tc()
    t_low = 0.02 * tc
    t_mev = KB_MEV_PER_K * t_low
    n_needed = solver._n_for_T(t_mev)
    spacing = 2*np.pi*t_mev
    t1 = time.time()
    print(f"lambda={lam}: Tc={tc:.5f} K (auto-N tc_time={t1-t0:.1f}s); "
          f"at t_low=0.02*Tc, auto N_for_solve_gap={n_needed}, spacing={spacing:.5f} meV, "
          f"needed to cover w_cut=500meV with gapratio's N=400: "
          f"N_min={500/spacing:.0f}", flush=True)
