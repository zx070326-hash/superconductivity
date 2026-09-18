"""Diagnostic: Why is 2Delta/kTc overestimated and BCS limit wrong?

Critical issue found:
- H3S theory: 2Delta/kTc = 4.55 vs experiment 3.43 (32% high)
- BCS limit: lambda=0.3, mu*=0.1 gives 3.24 vs expected 3.53 (8% low)

Possible causes:
1. Delta(omega_0) != Delta(0) - need to check frequency dependence
2. T=0.1*Tc is not low enough - need T -> 0 properly
3. Solver has a systematic bias
4. Definition mismatch with experiment

This script: test BCS limit carefully, check Delta frequency structure
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")

import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

print("="*70)
print("DIAGNOSTIC 1: BCS weak-coupling limit")
print("="*70)
print()
print("BCS theory: 2*Delta(0)/kTc = 2*pi/exp(gamma) = 3.5275")
print("where gamma = 0.5772... is Euler-Mascheroni constant")
print()

# Test with very weak coupling and mu*=0
lam_values = [0.1, 0.2, 0.3, 0.5]
print("Testing Einstein model at weak coupling, mu*=0:")
print()
print("lambda   Tc(K)   Delta(0) meV   2*Delta/kTc   Error vs BCS")
print("-" * 65)

for lam in lam_values:
    omega_e = 50.0  # meV
    spec = Spectrum.einstein(lam, omega_e)
    solver = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=400)

    tc = solver.tc()

    # Solve at very low T
    t_low = 0.02 * tc  # Even lower than before
    wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000, mixing=0.2)

    if conv:
        delta0 = delta[0]
        ratio = 2 * delta0 / (KB_MEV_PER_K * tc)
        error = (ratio / 3.5275 - 1) * 100
        print(f"{lam:6.2f}  {tc:6.2f}  {delta0:8.3f}      {ratio:7.3f}      {error:+6.1f}%")
    else:
        print(f"{lam:6.2f}  {tc:6.2f}  NOT CONVERGED")

print()
print("BCS limit: 2*Delta/kTc should approach 3.527 as lambda->0")
print()

print("="*70)
print("DIAGNOSTIC 2: Frequency structure of Delta")
print("="*70)
print()

# Check how Delta varies with Matsubara frequency
lam = 0.3
omega_e = 50.0
spec = Spectrum.einstein(lam, omega_e)
solver = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=200)
tc = solver.tc()
t_low = 0.05 * tc

wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000)

if conv:
    print(f"T = {t_low:.2f} K = {t_low/tc:.3f} Tc")
    print()
    print("First 10 Matsubara frequencies:")
    print("n    omega_n (meV)   Delta_n (meV)   Delta_n/Delta_0")
    print("-" * 60)
    for i in range(min(10, len(wn))):
        ratio = delta[i] / delta[0]
        print(f"{i:2d}   {wn[i]:8.2f}        {delta[i]:8.3f}        {ratio:6.3f}")
    print()

    # Check if Delta is roughly constant in the low-frequency region
    delta_avg_low = np.mean(delta[:5])
    print(f"Average of first 5: {delta_avg_low:.3f} meV")
    print(f"Delta[0] = {delta[0]:.3f} meV")
    print(f"Ratio = {delta_avg_low/delta[0]:.3f}")
    print()

print("="*70)
print("DIAGNOSTIC 3: What is the experimental 2Delta?")
print("="*70)
print()
print("From hydrides_experimental.csv:")
print("  H3S: gap_meV = 30 => 2*Delta = 60 meV")
print("  Tc = 203 K")
print("  2*Delta/(k_B*Tc) = 60/(0.08617*203) = 3.43")
print()
print("This is the FULL gap 2*Delta measured by tunneling.")
print("In BCS: Delta(T=0) is constant in k-space (s-wave).")
print("In Eliashberg: Delta(omega_n) varies, but Delta(0) is")
print("the physical gap at Fermi surface.")
print()
print("So we should compare Delta(omega_0) with experiment,")
print("which is what we're doing. The comparison is correct.")
print()

print("="*70)
print("DIAGNOSTIC 4: H3S with different N_matsubara")
print("="*70)
print()

# Load H3S spectrum
path = "D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-160GPa/EPW_anharmonic/H3S.a2f"
w_list, a2f_list = [], []
with open(path, 'r') as f:
    for line in f:
        p = line.split()
        if len(p) >= 2:
            try:
                w_list.append(float(p[0]))
                a2f_list.append(float(p[1]))
            except ValueError:
                continue

w = np.array(w_list)
a2f = np.array(a2f_list)
mask = w > 0
w, a2f = w[mask], a2f[mask]

spec = Spectrum.from_table(w, a2f)
lam = spec.lambda_total()
wlog = spec.omega_log()
w2 = spec.omega_2()

print(f"H3S: lambda={lam:.3f}, omega_log={wlog:.1f} meV")
print()

print("Testing N_matsubara convergence:")
print("N_mat   Tc(K)   Delta(0) meV   2*Delta/kTc")
print("-" * 50)

mu_star = 0.13
w_cut = 10.0 * w2

for n_mat in [400, 600, 800, 1000]:
    solver = Solver(spectrum=spec, mu_star=mu_star, w_cut=w_cut,
                   n_matsubara=n_mat)
    tc = solver.tc()

    t_low = 0.05 * tc
    wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000, mixing=0.25)

    if conv:
        delta0 = delta[0]
        ratio = 2 * delta0 / (KB_MEV_PER_K * tc)
        print(f"{n_mat:5d}  {tc:6.2f}  {delta0:8.3f}      {ratio:7.3f}")
    else:
        print(f"{n_mat:5d}  {tc:6.2f}  NOT CONVERGED")

print()

print("="*70)
print("CONCLUSION")
print("="*70)
print()
print("If BCS limit is wrong, there's a systematic issue.")
print("If BCS limit is correct, H3S overestimation is real and")
print("could mean:")
print("  1. mu* should be higher than 0.13")
print("  2. Spectrum needs anharmonic corrections")
print("  3. Non-Migdal vertex corrections lower the gap")
print("  4. DOS effects matter")
print()
print("This is the critical diagnostic to understand before proceeding.")
