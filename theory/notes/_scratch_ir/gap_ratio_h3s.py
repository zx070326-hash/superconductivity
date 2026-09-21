"""Gap ratio 2Delta/kTc for H3S: full nonlinear Eliashberg solution.

Roadmap:
1. Load H3S alpha^2F (160 GPa anharmonic, closest to experimental 155 GPa)
2. Solve nonlinear gap equation at T << Tc to get Delta(T->0)
3. Compute Tc via linearized equation
4. Compare 2Delta(0)/kTc with experimental tunneling value ~3.4
5. Test whether closed-form eigenvector profile can reproduce the solution

Strict level: numerical-evidence (float64, finite Matsubara grid)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, "D:/physics/superconductivity/theory/solvers")

import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

# Load H3S alpha^2F (160 GPa anharmonic)
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

print("=" * 70)
print("H3S 160 GPa anharmonic spectrum")
print("=" * 70)
print(f"lambda     = {lam:.4f}")
print(f"omega_log  = {wlog:.2f} meV = {wlog/KB_MEV_PER_K:.1f} K")
print(f"omega_2    = {w2:.2f} meV = {w2/KB_MEV_PER_K:.1f} K")
print(f"w_max      = {spec.w_max:.2f} meV")
print()

# Solver setup: mu* = 0.13 (from deficit budget work), w_cut = 10*sqrt(<om^2>)
mu_star = 0.13
w_cut = 10.0 * w2

solver = Solver(spectrum=spec, mu_star=mu_star, w_cut=w_cut, n_matsubara=600)

print("=" * 70)
print("Tc calculation (linearized)")
print("=" * 70)
print(f"mu*   = {mu_star:.3f}")
print(f"w_cut = {w_cut:.2f} meV = {w_cut/KB_MEV_PER_K:.1f} K")

tc_k = solver.tc()
print(f"Tc    = {tc_k:.2f} K")
print()

# Allen-Dynes for comparison
try:
    tc_ad = spec.allen_dynes(mu_star)
    print(f"Allen-Dynes Tc = {tc_ad:.2f} K (for reference)")
except Exception as e:
    print(f"Allen-Dynes: {e}")
print()

# Nonlinear gap at low temperature (T = 0.1*Tc to approximate T->0)
print("=" * 70)
print("Nonlinear gap solution at T = 0.1*Tc")
print("=" * 70)
t_low = 0.1 * tc_k
print(f"T = {t_low:.2f} K = {t_low/tc_k:.3f} Tc")

wn, delta, z, converged = solver.solve_gap(t_low, max_iter=3000, mixing=0.25)
if not converged:
    print("WARNING: Gap solution did not converge!")
else:
    print("Gap solution converged.")

delta0 = delta[0]  # Gap at first Matsubara frequency (meV)
ratio = 2 * delta0 / (KB_MEV_PER_K * tc_k)

print(f"Delta(omega_0) = {delta0:.3f} meV")
print(f"2*Delta(omega_0) / (k_B*Tc) = {ratio:.3f}")
print()

# BCS weak-coupling limit check
bcs_ratio = 3.527  # 2*Delta/kTc for BCS
print(f"BCS weak-coupling value = {bcs_ratio:.3f}")
print(f"Enhancement factor = {ratio/bcs_ratio:.3f}")
print()

# Experimental comparison
print("=" * 70)
print("Comparison with experiment")
print("=" * 70)
print("Experimental (H3S, 155 GPa, Du et al. 2025):")
print("  Tc ~ 203 K")
print("  2*Delta = 60 meV (tunneling)")
print("  2*Delta/kTc ~ 60/(0.08617*203) ~ 3.43")
print()
print(f"This calculation:")
print(f"  Tc = {tc_k:.1f} K")
print(f"  2*Delta/kTc = {ratio:.2f}")
print()

# Profile analysis: check how well the closed-form eigenvector fits
print("=" * 70)
print("Eigenvector profile analysis")
print("=" * 70)
print("Testing closed-form profile: tau_n = u_n * a^2/(a^2 + (n+1/2)^2)")
print("where u_n = (2n+1)^{-1/2}, a = alpha * omega_bar")
print()

# Normalize delta
delta_norm = delta / np.linalg.norm(delta)

# Test profile with different widths
n_idx = np.arange(len(delta))
u_n = 1.0 / np.sqrt(2*n_idx + 1)
n_half = n_idx + 0.5

# Estimate characteristic frequency scale from delta decay
# Find where delta drops to 1/e of its peak
delta_peak = delta[0]
idx_decay = np.where(delta < delta_peak/np.e)[0]
if len(idx_decay) > 0:
    n_decay = idx_decay[0]
    # omega_n = pi*T*(2n+1), so characteristic scale is pi*T*2*n_decay
    omega_char = np.pi * KB_MEV_PER_K * t_low * (2*n_decay + 1)
    print(f"Gap decay scale: n ~ {n_decay}, omega ~ {omega_char:.1f} meV")
    print()

# Sweep width parameter
alphas = np.linspace(0.5, 5.0, 50)
residuals = []
for alpha in alphas:
    a = alpha * wlog
    tau = u_n * a**2 / (a**2 + n_half**2)
    tau_norm = tau / np.linalg.norm(tau)
    res = np.linalg.norm(delta_norm - tau_norm)
    residuals.append(res)

residuals = np.array(residuals)
best_idx = np.argmin(residuals)
best_alpha = alphas[best_idx]
best_res = residuals[best_idx]

print(f"Best width parameter: alpha = {best_alpha:.3f}")
print(f"omega_log = {wlog:.2f} meV")
print(f"a = alpha * omega_log = {best_alpha * wlog:.2f} meV")
print(f"Residual = {best_res:.5f}")
print()

# Compare with linearized eigenvector at Tc
print("For comparison, testing at T = Tc (linearized regime):")
wn_tc, delta_tc, z_tc, conv_tc = solver.solve_gap(tc_k*0.999, max_iter=3000)
if conv_tc:
    delta_tc_norm = delta_tc / np.linalg.norm(delta_tc)
    n_tc = np.arange(len(delta_tc))
    u_tc = 1.0 / np.sqrt(2*n_tc + 1)
    n_half_tc = n_tc + 0.5

    residuals_tc = []
    for alpha in alphas:
        a = alpha * wlog
        tau = u_tc * a**2 / (a**2 + n_half_tc**2)
        tau_norm = tau / np.linalg.norm(tau)
        res = np.linalg.norm(delta_tc_norm - tau_norm)
        residuals_tc.append(res)

    residuals_tc = np.array(residuals_tc)
    best_idx_tc = np.argmin(residuals_tc)
    best_alpha_tc = alphas[best_idx_tc]
    best_res_tc = residuals_tc[best_idx_tc]

    print(f"  Best alpha = {best_alpha_tc:.3f}")
    print(f"  Residual = {best_res_tc:.5f}")
    print()

print("=" * 70)
print("Summary")
print("=" * 70)
print(f"[numerical-evidence] H3S @ 160 GPa, mu* = {mu_star:.2f}:")
print(f"  Tc = {tc_k:.1f} K")
print(f"  2*Delta(0)/kTc = {ratio:.2f}")
print(f"  Experiment: Tc ~ 203 K, 2*Delta/kTc ~ 3.43")
print(f"  Closed-form profile: alpha = {best_alpha:.2f}, residual = {best_res:.4f}")
print()
print("Next steps:")
print("1. Solve at even lower T (T->0 extrapolation)")
print("2. Test mu* sensitivity")
print("3. Test closed-form profile as ansatz for full nonlinear equation")
