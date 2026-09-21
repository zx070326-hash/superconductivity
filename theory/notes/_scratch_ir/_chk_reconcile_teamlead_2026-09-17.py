"""
Correction run after team-lead's cross-check flagged two mismatches:

1. My earlier lambda=0.2 "N=24064 needed" was WRONG: it came from calling
   solver.tc(n_matsubara=None) and letting the RuntimeError fire during the
   bisection's *bracket-widening* search at some irrelevant low trial T
   (T=0.1153 K), not at the intended t_low=0.02*Tc. Recomputing N-requirement
   directly from the converged Tc via the closed formula gives ~43700,
   matching team-lead's independent number.

2. My earlier lambda=0.5, N=400 reproduction gave ratio=4.454, not gapratio's
   reported 4.05. Root cause found: gap_ratio_key_findings.py (the script
   that actually produced "4.05") uses t_low=0.03*Tc, not t_low=0.02*Tc
   (which is what gap_ratio_diagnostic.py uses -- I grabbed the wrong
   script's parameter). With t_low=0.03*Tc, N=400 gives ratio=4.0514,
   matching 4.05 almost exactly.

This script: redo the N-convergence table for lambda=0.5 with the CORRECT
t_low=0.03*Tc (matching gap_ratio_key_findings.py exactly), and compute the
lambda=0.3 / lambda=0.2 N-requirement using converged Tc directly (bypassing
solver.tc(n_matsubara=None) entirely to avoid the bracket-widening trap).

Strict level: numerical-evidence (float64), independent rebuild.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, time
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")
import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

print("="*90)
print("PART A: lambda=0.5, CORRECTED t_low=0.03*Tc (matches gap_ratio_key_findings.py)")
print("="*90)
spec05 = Spectrum.einstein(0.5, 50.0)
solver_ref = Solver(spectrum=spec05, mu_star=0.0, w_cut=500.0, n_matsubara=400)
tc05 = solver_ref.tc()
print(f"Tc(lambda=0.5) = {tc05:.5f} K (N=400; Tc itself converges fast, see earlier run: "
      f"22.880/22.878/22.878/22.878 at N=400/800/1600/3200)")
for N in (400, 800, 1600, 3200):
    solver = Solver(spectrum=spec05, mu_star=0.0, w_cut=500.0, n_matsubara=N)
    t_low = 0.03 * tc05
    t0 = time.time()
    wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000, mixing=0.2)
    t1 = time.time()
    ratio = 2*delta[0]/(KB_MEV_PER_K*tc05)
    print(f"  N={N:5d}  delta[0]={delta[0]:.5f} meV  ratio={ratio:.4f}  conv={conv}  ({t1-t0:.1f}s)", flush=True)

print()
print("="*90)
print("PART B: N-requirement recomputed from CONVERGED Tc (not via buggy tc(n=None) call)")
print("="*90)
# Use Tc values already independently confirmed converged (this session + team-lead cross-check):
#   lambda=0.5: Tc=22.878 K (converges by N=400 itself)
#   lambda=0.3: Tc=5.7584 K (converged at N=1600, matches team-lead's Carbotte-script Tc=5.758 K)
#   lambda=0.2: Tc=1.0565 K (team-lead's independent N-convergence, N=800..3200 all agree to 0.008%)
tc_converged = {0.5: 22.878, 0.3: 5.7584, 0.2: 1.0565}
for lam, tc_k in tc_converged.items():
    for frac, label in ((0.02, "0.02*Tc (my earlier table)"), (0.03, "0.03*Tc (matches key_findings.py)")):
        t_low_k = frac * tc_k
        t_mev = KB_MEV_PER_K * t_low_k
        spacing = 2*np.pi*t_mev
        n_min_bare = 500.0/spacing            # bare coverage of w_cut, oversample=1
        n_auto = 3.0*n_min_bare + 30          # matches Solver default matsubara_oversample=3.0, margin=30
        print(f"lambda={lam}  Tc={tc_k} K  t_low={frac}*Tc: spacing={spacing:.6f} meV  "
              f"N_min(bare)={n_min_bare:.0f}  N_auto(oversample=3)={n_auto:.0f}   [{label}]")
