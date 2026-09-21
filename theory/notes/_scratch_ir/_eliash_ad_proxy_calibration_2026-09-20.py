# -*- coding: utf-8 -*-
"""AD proxy calibration on H3S/D3S: compare AD formula vs true Eliashberg on same parameters

Task: For 3 H3S/D3S rows with full spectra, compute:
  1. omega_log from full alpha2F (the only missing quantity)
  2. Four-channel decomposition using AD75 formula (no sqrt bug in f2)
  3. Compare delta_BRD(AD) vs delta_BRD(true) row by row

Reference true values: deficit_budget_h3s_2026-09-17.md lines 36/38
Convention: r = w_cut/omega_2 = 10, multiplicative normalization
Units: omega_2 in meV, Tc in K
"""
import os
os.environ['OMP_NUM_THREADS'] = '1'

import sys
sys.path.insert(0, 'theory/solvers')

import numpy as np
from eliashberg_iso import KB_MEV_PER_K

C_INF = 0.1827262477

# Reference true solutions (from validation gate, already verified)
reference = [
    ("H3S-200 anharm", 1.721, 149.26, 415.2, 327.4, 0.211, 283.3, 0.135, 199.0, 0.298),
    ("D3S-200 anharm", 1.784, 105.43, 298.6, 237.8, 0.204, 214.4, 0.098, 155.4, 0.275),
    ("H3S-160 anharm", 1.876, 137.87, 400.4, 323.1, 0.193, 290.7, 0.100, 211.2, 0.274),
]

# Data file paths
data_paths = {
    "H3S-200 anharm": "theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S-200 anharm": "theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/D3S-200GPa/EPW_anharmonic/D3S.a2f",
    "H3S-160 anharm": "theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-160GPa/EPW_anharmonic/H3S.a2f",
}

print("=" * 80)
print("AD proxy calibration: AD formula vs true Eliashberg on H3S/D3S")
print("=" * 80)
print()

# Load alpha2F and compute omega_log
def load_and_compute_moments(path):
    """Load alpha2F and compute lambda, omega_log, omega_2"""
    omega_list = []
    alpha2f_list = []
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or not line[0].isdigit():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    omega = float(parts[0])
                    a2f = float(parts[1])
                    if omega > 0:
                        omega_list.append(omega)
                        alpha2f_list.append(a2f)
                except ValueError:
                    break

    omega = np.array(omega_list)
    alpha2f = np.array(alpha2f_list)

    # Compute moments using trapezoidal integration
    from scipy import integrate
    lam = 2.0 * integrate.trapezoid(alpha2f / omega, omega)
    omega_log_meV = np.exp(2.0 * integrate.trapezoid(alpha2f / omega * np.log(omega), omega) / lam)
    omega_2_meV = np.sqrt(2.0 * integrate.trapezoid(omega * alpha2f, omega) / lam)

    return lam, omega_log_meV, omega_2_meV

# Allen-Dynes formula (AD75 Eq.34/35, NO sqrt bug in f2)
def allen_dynes_formula(lam, omega_log_meV, omega_2_meV, mu_star):
    """AD75 formula with correct f2 (no sqrt)"""
    r = omega_2_meV / omega_log_meV
    Lambda1 = 2.46 * (1 + 3.8 * mu_star)
    Lambda2 = 1.82 * (1 + 6.3 * mu_star) * r
    f1 = (1 + (lam / Lambda1)**1.5)**(1.0/3.0)
    f2 = 1 + (r - 1) * lam**2 / (lam**2 + Lambda2**2)  # NO sqrt here

    Tc_meV = (omega_log_meV / 1.2) * np.exp(-1.04 * (1 + lam) / (lam - mu_star * (1 + 0.62 * lam))) * f1 * f2
    return Tc_meV / KB_MEV_PER_K

print("Step 1: Computing omega_log from full alpha2F")
print("=" * 80)
print()

computed_moments = {}
for name in ["H3S-200 anharm", "D3S-200 anharm", "H3S-160 anharm"]:
    path = data_paths[name]
    lam_calc, omega_log_calc, omega_2_calc = load_and_compute_moments(path)
    computed_moments[name] = (lam_calc, omega_log_calc, omega_2_calc)

    # Get reference values
    ref_row = [r for r in reference if r[0] == name][0]
    lam_ref = ref_row[1]
    omega_2_ref = ref_row[2]

    print(f"{name}:")
    print(f"  lambda:    computed {lam_calc:.4f}, ref {lam_ref:.4f} (diff {abs(lam_calc-lam_ref)/lam_ref*100:.2f}%)")
    print(f"  omega_2:   computed {omega_2_calc:.2f} meV, ref {omega_2_ref:.2f} meV (diff {abs(omega_2_calc-omega_2_ref)/omega_2_ref*100:.2f}%)")
    print(f"  omega_log: computed {omega_log_calc:.2f} meV (NEW—not in reference table)")
    print()

print()
print("Step 2: Four-channel decomposition using AD formula")
print("=" * 80)
print()

results = []

print(f"{'Material':<20} {'Channel':<12} {'AD proxy':<12} {'True soln':<12} {'Diff':<12}")
print("-" * 80)

for name, lam_ref, omega_2_ref, ceil_ref, TcE0_ref, dFIN_ref, TcP0_ref, dBRD_ref, TcP_ref, dCOU_ref in reference:
    lam_calc, omega_log_calc, omega_2_calc = computed_moments[name]

    # Use computed values (more accurate than reference which may have been rounded)
    lam = lam_calc
    omega_log_meV = omega_log_calc
    omega_2_meV = omega_2_calc

    # Compute ceiling
    ceiling_K = C_INF * np.sqrt(lam) * omega_2_meV / KB_MEV_PER_K

    # AD formula four channels
    # delta_FIN: Einstein (omega_log = omega_2 => f2 = 1)
    Tc_E_AD = allen_dynes_formula(lam, omega_2_meV, omega_2_meV, 0.0)
    delta_FIN_AD = 1 - Tc_E_AD / ceiling_K

    # delta_BRD: Real spectrum (use omega_log) at mu*=0
    Tc_P0_AD = allen_dynes_formula(lam, omega_log_meV, omega_2_meV, 0.0)
    delta_BRD_AD = 1 - Tc_P0_AD / Tc_E_AD

    # delta_COU: Real spectrum at mu*=0.13
    Tc_P_AD = allen_dynes_formula(lam, omega_log_meV, omega_2_meV, 0.13)
    delta_COU_AD = 1 - Tc_P_AD / Tc_P0_AD

    # Print comparison
    print(f"{name:<20} {'delta_FIN':<12} {delta_FIN_AD*100:6.1f}%     {dFIN_ref*100:6.1f}%     {(delta_FIN_AD-dFIN_ref)*100:+6.1f} pp")
    print(f"{'':<20} {'delta_BRD':<12} {delta_BRD_AD*100:6.1f}%     {dBRD_ref*100:6.1f}%     {(delta_BRD_AD-dBRD_ref)*100:+6.1f} pp")
    print(f"{'':<20} {'delta_COU':<12} {delta_COU_AD*100:6.1f}%     {dCOU_ref*100:6.1f}%     {(delta_COU_AD-dCOU_ref)*100:+6.1f} pp")
    print()

    results.append({
        'name': name,
        'lambda': lam,
        'omega_log': omega_log_meV,
        'omega_2': omega_2_meV,
        'delta_FIN_AD': delta_FIN_AD,
        'delta_FIN_true': dFIN_ref,
        'delta_BRD_AD': delta_BRD_AD,
        'delta_BRD_true': dBRD_ref,
        'delta_COU_AD': delta_COU_AD,
        'delta_COU_true': dCOU_ref,
    })

print()
print("=" * 80)
print("SUMMARY: AD proxy bias")
print("=" * 80)
print()

delta_BRD_AD_vals = [r['delta_BRD_AD'] for r in results]
delta_BRD_true_vals = [r['delta_BRD_true'] for r in results]
bias = [(ad - true)*100 for ad, true in zip(delta_BRD_AD_vals, delta_BRD_true_vals)]

print(f"delta_BRD bias (AD - true):")
for i, name in enumerate([r['name'] for r in results]):
    print(f"  {name:<20} {bias[i]:+5.1f} pp")
print()
print(f"Mean bias: {np.mean(bias):+5.1f} pp")
print(f"Std:       {np.std(bias):5.1f} pp")
print()

print("Interpretation:")
if np.mean(bias) > 1.0:
    print(f"  AD proxy SYSTEMATICALLY OVERESTIMATES delta_BRD by ~{np.mean(bias):.1f} pp at lambda~1.7-1.9")
elif np.mean(bias) < -1.0:
    print(f"  AD proxy SYSTEMATICALLY UNDERESTIMATES delta_BRD by ~{abs(np.mean(bias)):.1f} pp at lambda~1.7-1.9")
else:
    print(f"  AD proxy shows NO SYSTEMATIC BIAS (mean {np.mean(bias):+.1f} pp within ±1 pp)")

print()
print("CRITICAL LIMITATION:")
print("  This calibration is on 3 rows, all Im3m H3S family, lambda in [1.72, 1.88].")
print("  AD75's 29 rows span lambda in [0.69, 2.59], 14 disjoint element sets, s-p metals.")
print("  **Proxy accuracy at lambda~1.8 DOES NOT IMPLY accuracy at lambda~1.2**")
print("  This test can FALSIFY the proxy (if mismatch is large) but cannot VALIDATE it")
print("  for the AD75 population.")
