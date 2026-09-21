#!/usr/bin/env python3
"""
Investigate the 2.4% ME bias between KAY linearization and full Eliashberg.
Focus on:
1. Which dataset (fsr_adia/fsr_na/fbw_adia/fbw_na) was compared
2. μ_c value in ME vs μ*=0.13 in our calculation
3. Z-channel treatment differences
4. Sign direction audit
"""
import numpy as np
import os

# Data directory
DATA_DIR = "D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic"

def extrapolate_tc_sqrtgap(temps, gaps, m=3):
    """
    Extrapolate Tc from Δ(T) using √-method: fit Δ² vs T to linear, find zero.
    m: number of last points to use
    """
    if len(temps) < m:
        m = len(temps)

    T_fit = temps[-m:]
    gap_sq_fit = gaps[-m:]**2

    # Linear fit: Δ² = a + b*T, Tc at Δ²=0 => Tc = -a/b
    coeffs = np.polyfit(T_fit, gap_sq_fit, 1)
    b, a = coeffs
    tc = -a / b

    return tc

def load_me_gap_file(filepath):
    """Load ME .imag_iso_gap0 file: columns are T(K), Δ(meV)"""
    data = np.loadtxt(filepath)
    temps = data[:, 0]
    gaps = data[:, 1]
    return temps, gaps

def get_tc_from_me_file(filepath, windows=[2,3,5]):
    """Get Tc from ME gap file using multiple extrapolation windows"""
    temps, gaps = load_me_gap_file(filepath)

    tcs = []
    for m in windows:
        tc = extrapolate_tc_sqrtgap(temps, gaps, m)
        tcs.append(tc)

    return np.array(tcs)

# Define all four datasets for H3S-200 anharmonic
datasets = {
    'fsr_adia': os.path.join(DATA_DIR, 'H3S-200GPa/EPW_anharmonic/fsr_adia/H3S.imag_iso_gap0'),
    'fsr_na': os.path.join(DATA_DIR, 'H3S-200GPa/EPW_anharmonic/fsr_na/H3S.imag_iso_gap0'),
    'fbw_adia': os.path.join(DATA_DIR, 'H3S-200GPa/EPW_anharmonic/fbw_adia/H3S.imag_iso_gap0'),
    'fbw_na': os.path.join(DATA_DIR, 'H3S-200GPa/EPW_anharmonic/fbw_na/H3S.imag_iso_gap0'),
}

print("="*70)
print("ME Tc extraction: H3S-200 anharmonic, all four groups")
print("="*70)

for name, filepath in datasets.items():
    if os.path.exists(filepath):
        tcs = get_tc_from_me_file(filepath, windows=[2,3,5])
        print(f"\n{name:12s}: windows [2,3,5] => {tcs[0]:.2f}, {tcs[1]:.2f}, {tcs[2]:.2f} K")
        print(f"              mean={np.mean(tcs):.2f} K, std={np.std(tcs):.3f} K")
    else:
        print(f"\n{name:12s}: FILE NOT FOUND")

print("\n" + "="*70)
print("Now check D3S-200, H3S-160, D3S-160 for fsr_adia only")
print("="*70)

# Check other datasets with fsr_adia
other_datasets = {
    'D3S-200': os.path.join(DATA_DIR, 'D3S-200GPa/EPW_anharmonic/fsr_adia/D3S.imag_iso_gap0'),
    'H3S-160': os.path.join(DATA_DIR, 'H3S-160GPa/EPW_anharmonic/fsr_adia/H3S.imag_iso_gap0'),
    'D3S-160': os.path.join(DATA_DIR, 'D3S-160GPa/EPW_anh/fsr_adia/D3S.imag_iso_gap0'),
}

for name, filepath in other_datasets.items():
    if os.path.exists(filepath):
        tcs = get_tc_from_me_file(filepath, windows=[2,3,5])
        print(f"\n{name:12s}: windows [2,3,5] => {tcs[0]:.2f}, {tcs[1]:.2f}, {tcs[2]:.2f} K")
        print(f"              mean={np.mean(tcs):.2f} K, std={np.std(tcs):.3f} K")
    else:
        print(f"\n{name:12s}: FILE NOT FOUND")

print("\n" + "="*70)
print("Our KAY Tc values (from deficit_budget table):")
print("="*70)
print("H3S-200 anharm: 203.9 K (at μ*=0.13)")
print("D3S-200 anharm: 158.8 K (at μ*=0.13)")
print("H3S-160 anharm: 215.4 K (at μ*=0.13)")
print("D3S-160 anharm: 164.1 K (at μ*=0.13)")
