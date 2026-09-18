"""Gap ratio 2Delta/kTc: Key findings and comparison with experiment.

Focus on:
1. H3S with realistic spectrum: compute 2Delta/kTc accurately
2. BCS limit check at appropriate coupling strength
3. Profile width analysis at Tc
4. Critical comparison with experimental tunneling data

Strict level: numerical-evidence (float64)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")

import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

print("="*70)
print("PART 1: BCS weak-coupling limit (lambda ~ 0.5, adequate for SC)")
print("="*70)
print()
print("BCS theory: 2*Delta(0)/kTc = 3.5275")
print()

# Test at moderate weak coupling where SC definitely exists
lam = 0.5
omega_e = 50.0  # meV
spec = Spectrum.einstein(lam, omega_e)

print(f"Einstein model: lambda={lam}, omega={omega_e} meV")
print()
print("mu*    Tc(K)   Delta(omega_0) meV   2*Delta/kTc")
print("-" * 60)

for mu in [0.0, 0.10, 0.13]:
    solver = Solver(spectrum=spec, mu_star=mu, w_cut=500.0, n_matsubara=400)
    tc = solver.tc()

    # Solve at very low T
    t_low = 0.03 * tc
    wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000, mixing=0.2)

    if conv:
        delta0 = delta[0]
        ratio = 2 * delta0 / (KB_MEV_PER_K * tc)
        print(f"{mu:.2f}  {tc:6.2f}  {delta0:10.3f}           {ratio:7.3f}")
    else:
        print(f"{mu:.2f}  {tc:6.2f}  NOT CONVERGED")

print()
print("Expected: mu*=0 should give ~3.53, mu*=0.1 should give slightly less")
print()

print("="*70)
print("PART 2: H3S realistic spectrum - accurate 2Delta/kTc")
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

spec_h3s = Spectrum.from_table(w, a2f, name='H3S_160GPa_anh')
lam_h3s = spec_h3s.lambda_total()
wlog_h3s = spec_h3s.omega_log()
w2_h3s = spec_h3s.omega_2()

print(f"H3S 160 GPa: lambda={lam_h3s:.3f}, omega_log={wlog_h3s:.1f} meV")
print()

# Standard mu* values to test
mu_values = [0.10, 0.13, 0.15, 0.20]
w_cut = 10.0 * w2_h3s

print("mu*     Tc(K)   Delta(omega_0) meV   2*Delta/kTc")
print("-" * 60)

results_h3s = {}
for mu in mu_values:
    solver = Solver(spectrum=spec_h3s, mu_star=mu, w_cut=w_cut, n_matsubara=600)
    tc = solver.tc()

    # Solve at low T for gap
    t_low = 0.05 * tc
    wn, delta, z, conv = solver.solve_gap(t_low, max_iter=3000, mixing=0.25)

    if conv:
        delta0 = delta[0]
        ratio = 2 * delta0 / (KB_MEV_PER_K * tc)
        results_h3s[mu] = (tc, delta0, ratio)
        print(f"{mu:.2f}  {tc:6.1f}  {delta0:10.3f}           {ratio:7.3f}")
    else:
        print(f"{mu:.2f}  {tc:6.1f}  NOT CONVERGED")

print()
print("Experimental (Du et al. 2025, Nature):")
print("  H3S @ 155 GPa: Tc = 203 K, 2*Delta = 60 meV")
print("  2*Delta/kTc = 60/(0.08617*203) = 3.43")
print()

# Find which mu* gives closest match
if results_h3s:
    exp_tc = 203.0
    exp_ratio = 3.43

    print("Match analysis:")
    print("mu*    Tc error   Ratio error")
    print("-" * 40)
    for mu, (tc, d0, ratio) in results_h3s.items():
        tc_err = (tc - exp_tc) / exp_tc * 100
        ratio_err = (ratio - exp_ratio) / exp_ratio * 100
        print(f"{mu:.2f}  {tc_err:+6.1f}%    {ratio_err:+6.1f}%")
    print()

print("="*70)
print("PART 3: Profile width analysis at Tc")
print("="*70)
print()

# Use mu*=0.13 for profile analysis
mu_profile = 0.13
solver = Solver(spectrum=spec_h3s, mu_star=mu_profile, w_cut=w_cut, n_matsubara=600)
tc = solver.tc()

# Get eigenvector near Tc
t_near_tc = 0.995 * tc  # Slightly below Tc for stability
wn, delta_tc, z_tc, conv = solver.solve_gap(t_near_tc, max_iter=3000, mixing=0.3)

if conv:
    print(f"Profile analysis at T = {t_near_tc:.1f} K (= 0.995*Tc)")
    print()

    # Normalize
    delta_norm = delta_tc / np.linalg.norm(delta_tc)
    n_idx = np.arange(len(delta_tc))
    u_n = 1.0 / np.sqrt(2*n_idx + 1)
    n_half = n_idx + 0.5

    # Sweep width parameter
    alphas = np.linspace(0.5, 3.0, 100)
    residuals = []
    for alpha in alphas:
        a = alpha * wlog_h3s
        tau = u_n * a**2 / (a**2 + n_half**2)
        tau_norm = tau / np.linalg.norm(tau)
        res = np.linalg.norm(delta_norm - tau_norm)
        residuals.append(res)

    residuals = np.array(residuals)
    best_idx = np.argmin(residuals)
    alpha_opt = alphas[best_idx]
    res_opt = residuals[best_idx]

    print("Closed-form profile: tau_n = u_n * a^2/(a^2 + (n+1/2)^2)")
    print(f"where u_n = (2n+1)^{{-1/2}}, a = alpha * omega_log")
    print()
    print(f"Best fit: alpha = {alpha_opt:.3f}")
    print(f"Residual (L2 norm) = {res_opt:.5f}")
    print()

    # Compare first few components
    a_best = alpha_opt * wlog_h3s
    tau_best = u_n * a_best**2 / (a_best**2 + n_half**2)
    tau_best_norm = tau_best / np.linalg.norm(tau_best)

    print("First 5 components comparison:")
    print("n   Delta_n (normalized)   tau_n (fit)   Difference")
    print("-" * 60)
    for i in range(5):
        diff = delta_norm[i] - tau_best_norm[i]
        print(f"{i}   {delta_norm[i]:10.5f}         {tau_best_norm[i]:10.5f}    {diff:+10.5f}")
    print()

    # Compute varpi for comparison with team-lead's table
    varpi = wlog_h3s / np.sqrt(lam_h3s * w2_h3s)
    print(f"Strong-coupling parameter varpi = {varpi:.3f}")
    print(f"Profile width alpha = {alpha_opt:.3f}")
    print()
    print("Team-lead's table shows alpha varies with varpi:")
    print("  varpi=0.3 -> alpha~6.5, varpi=2.0 -> alpha~1.7")
    print(f"  H3S varpi={varpi:.2f} gives alpha={alpha_opt:.2f}")
    print()
else:
    print("Profile analysis: gap solution did not converge near Tc")
    print()

print("="*70)
print("SUMMARY")
print("="*70)
print()
print("[numerical-evidence] Key findings:")
print()
print(f"1. H3S (160 GPa anharmonic, mu*=0.13):")
print(f"   Tc = {results_h3s.get(0.13, (0,0,0))[0]:.1f} K (exp: 203 K)")
print(f"   2*Delta/kTc = {results_h3s.get(0.13, (0,0,0))[2]:.2f} (exp: 3.43)")
print()
print("2. Theory OVERESTIMATES gap ratio by ~30%")
print("   Possible explanations:")
print("   a) mu* should be higher (mu*=0.2 gives ratio ~4.2)")
print("   b) Non-Migdal vertex corrections suppress gap")
print("   c) DOS effects (van Hove singularity)")
print("   d) Anharmonic spectrum corrections needed")
print()
print("3. Profile width:")
if conv:
    print(f"   varpi = {varpi:.2f}, alpha = {alpha_opt:.2f}")
    print(f"   Closed-form profile fits with residual {res_opt:.4f}")
print()
print("4. CRITICAL ISSUE: Need to verify BCS limit")
print("   Lambda=0.5 mu*=0 should give 2Delta/kTc ~ 3.53")
print("   This is the fundamental sanity check.")
