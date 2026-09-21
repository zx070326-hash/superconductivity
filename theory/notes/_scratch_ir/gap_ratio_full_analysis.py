"""Full gap ratio analysis for H3S: T->0 extrapolation, profile ansatz test,
strong-coupling parameter exploration.

Key findings to establish:
1. True 2Delta(0)/kTc via T->0 extrapolation
2. Whether closed-form eigenvector profile reproduces the gap
3. Relationship between profile width alpha(varpi) and Tc/omega_log
4. mu* sensitivity of the gap ratio

Strict level: numerical-evidence (float64)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")

import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

# Load H3S alpha^2F
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

spec = Spectrum.from_table(w, a2f, name='H3S_160GPa_anh')
lam = spec.lambda_total()
wlog = spec.omega_log()
w2 = spec.omega_2()

print("="*70)
print("PART 1: T->0 extrapolation for accurate Delta(0)")
print("="*70)
print(f"Spectrum: lambda={lam:.4f}, omega_log={wlog:.2f} meV")
print()

mu_star = 0.13
w_cut = 10.0 * w2
solver = Solver(spectrum=spec, mu_star=mu_star, w_cut=w_cut, n_matsubara=600)
tc_k = solver.tc()
print(f"Tc = {tc_k:.2f} K")
print()

# Solve at multiple temperatures to extrapolate to T=0
t_ratios = [0.5, 0.3, 0.2, 0.15, 0.1, 0.05]
deltas_0 = []
temps = []

print("T/Tc    T(K)    Delta(omega_0) meV    2*Delta/kTc")
print("-" * 55)
for ratio in t_ratios:
    t = ratio * tc_k
    wn, delta, z, conv = solver.solve_gap(t, max_iter=3000, mixing=0.25)
    if conv:
        delta0 = delta[0]
        ratio_val = 2 * delta0 / (KB_MEV_PER_K * tc_k)
        deltas_0.append(delta0)
        temps.append(t)
        print(f"{ratio:.3f}  {t:6.2f}  {delta0:10.3f}          {ratio_val:.3f}")
    else:
        print(f"{ratio:.3f}  {t:6.2f}  NOT CONVERGED")

print()

# Linear extrapolation from lowest two temperatures
if len(deltas_0) >= 2:
    # Delta(T) ~ Delta(0) * (1 - (T/Tc)^2) near T=0 (BCS-like)
    # But use simple linear extrapolation in Delta vs T
    t1, t2 = temps[-2], temps[-1]
    d1, d2 = deltas_0[-2], deltas_0[-1]
    delta_extrap = d2 - (d2 - d1)/(t2 - t1) * t2
    ratio_extrap = 2 * delta_extrap / (KB_MEV_PER_K * tc_k)

    print(f"Linear extrapolation to T=0:")
    print(f"  Delta(0) ~ {delta_extrap:.3f} meV")
    print(f"  2*Delta(0)/kTc ~ {ratio_extrap:.3f}")
    print()

print("Experimental value: 2*Delta/kTc ~ 3.43")
print()

print("="*70)
print("PART 2: Profile width vs strong-coupling parameter")
print("="*70)
print()

# Solve at Tc (linearized) to get eigenvector
t_near_tc = 0.999 * tc_k
wn, delta_tc, z_tc, conv = solver.solve_gap(t_near_tc, max_iter=3000)

if conv:
    delta_norm = delta_tc / np.linalg.norm(delta_tc)
    n_idx = np.arange(len(delta_tc))
    u_n = 1.0 / np.sqrt(2*n_idx + 1)
    n_half = n_idx + 0.5

    # Sweep width to find best fit
    alphas = np.linspace(0.3, 3.0, 100)
    residuals = []
    for alpha in alphas:
        a = alpha * wlog
        tau = u_n * a**2 / (a**2 + n_half**2)
        tau_norm = tau / np.linalg.norm(tau)
        res = np.linalg.norm(delta_norm - tau_norm)
        residuals.append(res)

    residuals = np.array(residuals)
    best_idx = np.argmin(residuals)
    alpha_opt = alphas[best_idx]
    res_opt = residuals[best_idx]

    print(f"Closed-form profile: tau_n = u_n * a^2/(a^2 + (n+1/2)^2)")
    print(f"where u_n = (2n+1)^{-1/2}, a = alpha * omega_log")
    print()
    print(f"Best fit: alpha = {alpha_opt:.3f}")
    print(f"Residual (L2 norm) = {res_opt:.5f}")
    print()

    # Compute varpi = omega_log / sqrt(lambda * omega_2)
    varpi = wlog / np.sqrt(lam * w2)
    print(f"Strong-coupling parameters:")
    print(f"  lambda = {lam:.3f}")
    print(f"  Tc/omega_log = {tc_k*KB_MEV_PER_K/wlog:.4f}")
    print(f"  varpi = omega_log/sqrt(lambda*omega_2) = {varpi:.3f}")
    print()

    # Compare with team-lead's table
    print("Team-lead's profile width table (Einstein spectrum):")
    print("varpi  alpha_opt  residual")
    table = [
        (0.3, 6.47, 0.0046),
        (0.5, 4.03, 0.0099),
        (1.0, 2.32, 0.0149),
        (2.0, 1.67, 0.0046),
        (4.0, 1.42, 0.0073),
        (8.48, 1.31, 0.0112),
        (15.0, 1.25, 0.0118),
        (20.0, 1.22, 0.0123),
    ]
    for v, a, r in table:
        print(f"{v:5.2f}  {a:7.2f}     {r:.4f}")
    print()

    print(f"H3S (realistic spectrum): varpi={varpi:.3f}, alpha={alpha_opt:.3f}")
    print(f"This is varpi~{varpi:.1f}, between the table entries.")
    print()

print("="*70)
print("PART 3: mu* sensitivity")
print("="*70)
print()

mu_values = [0.10, 0.13, 0.15, 0.18]
print("mu*     Tc(K)   2*Delta/kTc  Tc_deficit")
print("-" * 45)

tc_mu0 = None
for mu in mu_values:
    solver_mu = Solver(spectrum=spec, mu_star=mu, w_cut=w_cut, n_matsubara=600)
    tc_mu = solver_mu.tc()

    # Solve at low T for gap
    t_low = 0.1 * tc_mu
    wn, delta, z, conv = solver_mu.solve_gap(t_low, max_iter=3000, mixing=0.25)

    if conv:
        delta0 = delta[0]
        ratio_val = 2 * delta0 / (KB_MEV_PER_K * tc_mu)

        if mu == 0.10:
            tc_mu0 = tc_mu

        if tc_mu0 is not None:
            deficit = (1 - tc_mu/tc_mu0) * 100
            print(f"{mu:.2f}  {tc_mu:6.1f}  {ratio_val:8.3f}     {deficit:5.1f}%")
        else:
            print(f"{mu:.2f}  {tc_mu:6.1f}  {ratio_val:8.3f}        -")
    else:
        print(f"{mu:.2f}  {tc_mu:6.1f}  NOT CONVERGED")

print()
print("Experimental: Tc=203K with mu*=0.13 gives ~3% deficit vs theory")
print("(210.7K theory vs 203K experiment)")
print()

print("="*70)
print("PART 4: BCS weak-coupling limit check")
print("="*70)
print()

# Einstein model at weak coupling
lam_weak = 0.3
omega_e = 50.0  # meV
spec_weak = Spectrum.einstein(lam_weak, omega_e)

for mu in [0.0, 0.10]:
    solver_weak = Solver(spectrum=spec_weak, mu_star=mu,
                        w_cut=10*omega_e, n_matsubara=400)
    tc_weak = solver_weak.tc()

    t_low = 0.05 * tc_weak
    wn, delta, z, conv = solver_weak.solve_gap(t_low, max_iter=2000)

    if conv:
        delta0 = delta[0]
        ratio_val = 2 * delta0 / (KB_MEV_PER_K * tc_weak)
        print(f"lambda={lam_weak}, mu*={mu:.2f}: 2*Delta/kTc = {ratio_val:.3f}")
    else:
        print(f"lambda={lam_weak}, mu*={mu:.2f}: NOT CONVERGED")

print()
print("BCS weak-coupling prediction: 2*Delta/kTc = 3.527")
print("This is the crucial limit check!")
print()

print("="*70)
print("SUMMARY")
print("="*70)
print()
print("[numerical-evidence] H3S 160 GPa anharmonic spectrum:")
print(f"  lambda = {lam:.3f}")
print(f"  omega_log = {wlog:.1f} meV")
print(f"  Tc = {tc_k:.1f} K (mu*=0.13)")
print()
if len(deltas_0) >= 2:
    print(f"  2*Delta(0)/kTc = {ratio_extrap:.2f} (T->0 extrapolation)")
else:
    print(f"  2*Delta(0)/kTc = {2*deltas_0[0]/(KB_MEV_PER_K*tc_k):.2f} (T=0.05*Tc)")
print(f"  Experiment: 2*Delta/kTc ~ 3.43 (Du et al. 2025)")
print()
print(f"  Profile width: alpha = {alpha_opt:.2f}")
print(f"  Profile residual: {res_opt:.4f}")
print()
print("Next: Test closed-form profile as ansatz for self-consistent equations")
