#!/usr/bin/env python3
"""
Check ALL ME gap files including muc_0.12 subdirectory.
This is critical: the default might be muc=0.16, but there's a muc_0.12 subdirectory!
"""
import numpy as np
import os

DATA_DIR = "D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic"

def extrapolate_tc_sqrtgap(temps, gaps, m=3):
    """Extrapolate Tc from Δ(T) using √-method"""
    if len(temps) < m:
        m = len(temps)
    T_fit = temps[-m:]
    gap_sq_fit = gaps[-m:]**2
    coeffs = np.polyfit(T_fit, gap_sq_fit, 1)
    b, a = coeffs
    tc = -a / b
    return tc

def load_me_gap_file(filepath):
    """Load ME .imag_iso_gap0 file"""
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

# Check H3S-200 with BOTH muc values
print("="*70)
print("H3S-200GPa anharmonic: muc comparison")
print("="*70)

groups = ['fsr_adia', 'fsr_na', 'fbw_adia', 'fbw_na']
muc_dirs = ['', 'muc_0.12/']

for muc_label in ['muc=0.16 (default)', 'muc=0.12']:
    muc_dir = muc_dirs[0 if muc_label.startswith('muc=0.16') else 1]
    print(f"\n{muc_label}:")
    print("-" * 70)

    for group in groups:
        filepath = os.path.join(DATA_DIR, f'H3S-200GPa/EPW_anharmonic/{muc_dir}{group}/H3S.imag_iso_gap0')
        if os.path.exists(filepath):
            tcs = get_tc_from_me_file(filepath, windows=[2,3,5])
            mean_tc = np.mean(tcs)
            std_tc = np.std(tcs)
            print(f"  {group:12s}: {mean_tc:7.2f} K  (windows: {tcs[0]:.2f}, {tcs[1]:.2f}, {tcs[2]:.2f}; std={std_tc:.3f})")
        else:
            print(f"  {group:12s}: FILE NOT FOUND")

# Now check D3S-200, H3S-160, D3S-160 with both muc values
print("\n" + "="*70)
print("Other datasets with muc comparison (fsr_adia only)")
print("="*70)

datasets_info = [
    ('D3S-200GPa', 'EPW_anharmonic', 'D3S'),
    ('H3S-160GPa', 'EPW_anharmonic', 'H3S'),
    ('D3S-160GPa', 'EPW_anh', 'D3S'),
]

for compound_dir, epw_subdir, compound_name in datasets_info:
    print(f"\n{compound_dir}:")
    print("-" * 70)

    for muc_label in ['muc=0.16 (default)', 'muc=0.12']:
        muc_dir = '' if muc_label.startswith('muc=0.16') else 'muc_0.12/'
        filepath = os.path.join(DATA_DIR, f'{compound_dir}/{epw_subdir}/{muc_dir}fsr_adia/{compound_name}.imag_iso_gap0')

        if os.path.exists(filepath):
            tcs = get_tc_from_me_file(filepath, windows=[2,3,5])
            mean_tc = np.mean(tcs)
            print(f"  {muc_label:20s}: {mean_tc:7.2f} K")
        else:
            print(f"  {muc_label:20s}: FILE NOT FOUND")

print("\n" + "="*70)
print("Compare to team-lead's reported ME Tc values:")
print("="*70)
print("  H3S-200 anharm: 208.2 K (team-lead)")
print("  D3S-200 anharm: 163.0 K (team-lead)")
print("  H3S-160 anharm: 220.1 K (team-lead)")
print("  D3S-160 anharm: 168.8 K (team-lead)")
print("\nOur KAY Tc (μ*=0.13):")
print("  H3S-200: 203.9 K")
print("  D3S-200: 158.8 K")
print("  H3S-160: 215.4 K")
print("  D3S-160: 164.1 K")
