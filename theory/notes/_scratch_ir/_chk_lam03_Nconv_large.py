"""
Team-lead's requested check: lambda=0.3 N-convergence extrapolation, at N large
enough to approach the auto-estimated requirement (N_auto=16067 at t_low=0.03*Tc,
matching gap_ratio_key_findings.py's actual t_low convention -- see
_chk_reconcile_teamlead_2026-09-17.py PART B). This is the test that's supposed
to tell us whether the method is sound once N is adequate (Carbotte reference at
lambda=0.3 is 3.5470, close to BCS 3.53).

Uses the pre-established converged Tc=5.7584 K directly (from N=1600 in the
earlier _chk_gapratio_weak_coupling_Nconv.py run, cross-checked against
team-lead's independent Carbotte-script Tc=5.758 K) and calls solve_gap directly
with n_matsubara=N fixed, SKIPPING the tc() bisection (which would otherwise
redo ~10-30 expensive eig calls per N). This is a legitimate shortcut because
Tc itself was already shown to converge by N~1600 (same behavior seen at
lambda=0.5: Tc flat to 5 digits from N=400 up), while only the T=0.03*Tc gap
value needs the larger N.

Strict level: numerical-evidence (float64).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, time
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")
import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

TC_LAM03 = 5.7584  # K, converged (see docstring)
spec = Spectrum.einstein(0.3, 50.0)
t_low = 0.03 * TC_LAM03

for N in (400, 800, 1600, 3200, 8000, 16000):
    solver = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=N)
    t0 = time.time()
    wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000, mixing=0.2)
    t1 = time.time()
    ratio = 2*delta[0]/(KB_MEV_PER_K*TC_LAM03)
    print(f"N={N:6d}  delta[0]={delta[0]:.6f} meV  ratio={ratio:.4f}  conv={conv}  ({t1-t0:.1f}s)", flush=True)
