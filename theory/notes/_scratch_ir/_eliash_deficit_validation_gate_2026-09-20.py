# -*- coding: utf-8 -*-
"""Validation gate: reproduce 3-row true Eliashberg decomposition before attempting AD75 rows

Reference: theory/notes/materials/deficit_budget_h3s_2026-09-17.md lines 36/38
Convention: r = w_cut/omega_2 = 10, multiplicative normalization
Units: omega_2 in meV, Tc in K
"""
import os
os.environ['OMP_NUM_THREADS'] = '1'

import sys
sys.path.insert(0, 'theory/solvers')

import numpy as np
from eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

C_INF = 0.1827262477

# Reference table from deficit_budget_h3s_2026-09-17.md
# Dataset | lambda | omega_2(meV) | ceiling(K) | Tc_E(0) | delta_FIN | Tc_P(0) | delta_BRD | Tc_P(.13) | delta_COU
reference = [
    ("H3S-200 anharm", 1.721, 149.26, 415.2, 327.4, 0.211, 283.3, 0.135, 199.0, 0.298),
    ("D3S-200 anharm", 1.784, 105.43, 298.6, 237.8, 0.204, 214.4, 0.098, 155.4, 0.275),
    ("H3S-160 anharm", 1.876, 137.87, 400.4, 323.1, 0.193, 290.7, 0.100, 211.2, 0.274),
]

print("=" * 80)
print("VALIDATION GATE: 3-row true Eliashberg decomposition")
print("Reference: deficit_budget_h3s_2026-09-17.md")
print("=" * 80)
print()

# Load full spectra
print("Loading full alpha2F spectra...")

# Check if data files exist
import os
h3s_200_path = "theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-200GPa/EPW_anharmonic/H3S.a2f"
d3s_200_path = "theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/D3S-200GPa/EPW_anharmonic/D3S.a2f"
h3s_160_path = "theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-160GPa/EPW_anharmonic/H3S.a2f"

if not all(os.path.exists(p) for p in [h3s_200_path, d3s_200_path, h3s_160_path]):
    print("ERROR: Full spectra files not found. Expected:")
    print(f"  {h3s_200_path}")
    print(f"  {d3s_200_path}")
    print(f"  {h3s_160_path}")
    print()
    print("These files are needed for validation. Cannot proceed without them.")
    sys.exit(1)

# Load spectra
def load_spectrum(path):
    """Load alpha2F from EPW .a2f file (3 columns: omega_meV, alpha2F, integrated_lambda)"""
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
                    break  # Hit footer text, stop reading
    return np.array(omega_list), np.array(alpha2f_list)

spectra = {
    "H3S-200 anharm": load_spectrum(h3s_200_path),
    "D3S-200 anharm": load_spectrum(d3s_200_path),
    "H3S-160 anharm": load_spectrum(h3s_160_path),
}

print("Spectra loaded successfully.")
print()

print("=" * 80)
print("Running validation...")
print("=" * 80)
print()

results = []
all_pass = True

for name, lam, omega2_meV, ceil_ref, TcE0_ref, dFIN_ref, TcP0_ref, dBRD_ref, TcP_ref, dCOU_ref in reference:
    print(f"Processing: {name}")
    print(f"  lambda = {lam:.3f}, omega_2 = {omega2_meV:.2f} meV")

    # Construct spectrum
    omega_meV, alpha2f = spectra[name]
    spec = Spectrum.from_table(omega_meV, alpha2f, name=name)

    # Verify moments
    lam_computed = spec.lambda_total()
    omega2_computed = spec.omega_2()

    print(f"  Moment verification:")
    print(f"    lambda: expected {lam:.4f}, computed {lam_computed:.4f} (diff {abs(lam_computed-lam)/lam*100:.2f}%)")
    print(f"    omega_2: expected {omega2_meV:.2f} meV, computed {omega2_computed:.2f} meV (diff {abs(omega2_computed-omega2_meV)/omega2_meV*100:.2f}%)")

    # Convention: r = w_cut/omega_2 = 10
    w_cut = 10.0 * omega2_meV

    # Compute ceiling
    ceiling_K = C_INF * np.sqrt(lam) * omega2_meV / KB_MEV_PER_K

    # delta_FIN: Einstein spectrum at mu*=0
    spec_einstein = Spectrum.einstein(lam=lam, w0=omega2_meV)
    solver_E = Solver(spectrum=spec_einstein, mu_star=0.0, w_cut=w_cut, n_matsubara=1600)
    Tc_E0 = solver_E.tc()
    delta_FIN = 1 - Tc_E0 / ceiling_K

    # delta_BRD: Real spectrum at mu*=0
    solver_P0 = Solver(spectrum=spec, mu_star=0.0, w_cut=w_cut, n_matsubara=1600)
    Tc_P0 = solver_P0.tc()
    delta_BRD = 1 - Tc_P0 / Tc_E0

    # delta_COU: Real spectrum at mu*=0.13
    solver_P = Solver(spectrum=spec, mu_star=0.13, w_cut=w_cut, n_matsubara=1600)
    Tc_P = solver_P.tc()
    delta_COU = 1 - Tc_P / Tc_P0

    # Check multiplicative closure: f_pred = product of (1 - delta_i)
    f_pred = (1 - delta_FIN) * (1 - delta_BRD) * (1 - delta_COU)
    f_obs = Tc_P / ceiling_K

    print(f"  Results:")
    print(f"    Ceiling:     {ceiling_K:6.1f} K  (ref: {ceil_ref:6.1f} K, diff {abs(ceiling_K-ceil_ref):.1f} K)")
    print(f"    Tc_E(0):     {Tc_E0:6.1f} K  (ref: {TcE0_ref:6.1f} K, diff {abs(Tc_E0-TcE0_ref):.1f} K)")
    print(f"    delta_FIN:   {delta_FIN*100:5.1f}%  (ref: {dFIN_ref*100:5.1f}%, diff {abs(delta_FIN-dFIN_ref)*100:.1f} pp)")
    print(f"    Tc_P(0):     {Tc_P0:6.1f} K  (ref: {TcP0_ref:6.1f} K, diff {abs(Tc_P0-TcP0_ref):.1f} K)")
    print(f"    delta_BRD:   {delta_BRD*100:5.1f}%  (ref: {dBRD_ref*100:5.1f}%, diff {abs(delta_BRD-dBRD_ref)*100:.1f} pp)")
    print(f"    Tc_P(0.13):  {Tc_P:6.1f} K  (ref: {TcP_ref:6.1f} K, diff {abs(Tc_P-TcP_ref):.1f} K)")
    print(f"    delta_COU:   {delta_COU*100:5.1f}%  (ref: {dCOU_ref*100:5.1f}%, diff {abs(delta_COU-dCOU_ref)*100:.1f} pp)")
    print(f"    Closure: f_pred = {f_pred:.4f}, f_obs = {f_obs:.4f} (diff {abs(f_pred-f_obs)*100:.2f}%)")

    # Check tolerances (2 K for Tc, 1 pp for deltas)
    Tc_E0_pass = abs(Tc_E0 - TcE0_ref) < 2.0
    Tc_P0_pass = abs(Tc_P0 - TcP0_ref) < 2.0
    Tc_P_pass = abs(Tc_P - TcP_ref) < 2.0
    dFIN_pass = abs(delta_FIN - dFIN_ref) < 0.01
    dBRD_pass = abs(delta_BRD - dBRD_ref) < 0.01
    dCOU_pass = abs(delta_COU - dCOU_ref) < 0.01

    row_pass = all([Tc_E0_pass, Tc_P0_pass, Tc_P_pass, dFIN_pass, dBRD_pass, dCOU_pass])

    if row_pass:
        print(f"  STATUS: PASS")
    else:
        print(f"  STATUS: FAIL")
        if not Tc_E0_pass: print(f"    Tc_E(0) out of tolerance")
        if not Tc_P0_pass: print(f"    Tc_P(0) out of tolerance")
        if not Tc_P_pass: print(f"    Tc_P(0.13) out of tolerance")
        if not dFIN_pass: print(f"    delta_FIN out of tolerance")
        if not dBRD_pass: print(f"    delta_BRD out of tolerance")
        if not dCOU_pass: print(f"    delta_COU out of tolerance")
        all_pass = False

    print()

    results.append({
        'name': name,
        'Tc_E0': Tc_E0,
        'Tc_P0': Tc_P0,
        'Tc_P': Tc_P,
        'delta_FIN': delta_FIN,
        'delta_BRD': delta_BRD,
        'delta_COU': delta_COU,
        'pass': row_pass
    })

print("=" * 80)
print("VALIDATION GATE SUMMARY")
print("=" * 80)

if all_pass:
    print("STATUS: ALL PASSED")
    print()
    print("Pipeline validated. Ready to proceed to AD75 29-row analysis.")
else:
    print("STATUS: VALIDATION FAILED")
    print()
    print("Pipeline has implementation errors. DO NOT proceed to AD75 rows.")
    print("Fix the discrepancies above before continuing.")
    sys.exit(1)
