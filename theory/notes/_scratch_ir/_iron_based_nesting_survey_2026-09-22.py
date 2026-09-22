#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Iron-based superconductor geometric nesting survey - Framework

Team C Phase 2: Simplified 2-orbital tight-binding models for 5 materials

IMPORTANT CAVEAT:
- Exact Kuroki (2008) / Graser (2009,2010) hopping parameters (t_ij) are
  NOT accessible (web search budget exhausted, arXiv blocked, no local PDF).
- This script uses a SIMPLIFIED, QUALITATIVE 2-orbital model
  (d_xz, d_yz only) capturing the essential physics:
    - Disconnected Fermi surface (Kuroki 2008 title)
    - Gamma-centered hole pockets + M-centered electron pockets
- Material-to-material variation is modeled via a SINGLE effective
  "correlation strength" parameter (bandwidth narrowing), calibrated
  qualitatively to known Fermi surface size trends from ARPES literature
  (LaFeAsO/BaFe2As2/LiFeAs have larger, more 3D-like pockets than FeSe;
  FeTe0.5Se0.5 has strong pocket reconstruction).

STRICTNESS LEVEL: numerical-evidence, QUALITATIVE ONLY.
Do NOT quote ξ(q) absolute values as literature-grade for any material
other than FeSe (which reuses Team B's validated 5-orbital model).

Physical basis for 2-orbital model (standard, textbook-level):
  H(k) = [[eps_xz(k) - mu,  eps_xy(k)      ],
          [eps_xy(k),        eps_yz(k) - mu]]

  eps_xz(k) = -2t1*cos(kx) - 2t2*cos(ky) - 4t3*cos(kx)*cos(ky)
  eps_yz(k) = -2t2*cos(kx) - 2t1*cos(ky) - 4t3*cos(kx)*cos(ky)
  eps_xy(k) = -4t4*sin(kx)*sin(ky)

This is the standard Raghu et al. (arXiv:0810.0019) minimal 2-orbital
model form, itself a simplification of Kuroki's 5-orbital model,
capturing Gamma hole pockets + M/(pi,0)+(0,pi) electron pockets.
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from numpy.linalg import eigvalsh
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# 2-orbital tight-binding model (Raghu et al. form)
# ==============================================================================

def build_2orbital_hamiltonian(kx, ky, t1, t2, t3, t4, mu):
    """
    Standard 2-orbital (d_xz, d_yz) tight-binding model for Fe-pnictides.

    Reference form: Raghu, Qi, Liu, Scalapino, Zhang, PRB 77, 220503(R) (2008)
    arXiv:0804.1113 / related to arXiv:0810.0019

    Parameters t1-t4 control:
    - t1: nearest-neighbor same-orbital (x-direction for xz)
    - t2: nearest-neighbor same-orbital (y-direction for xz, weaker)
    - t3: next-nearest-neighbor same-orbital (diagonal)
    - t4: next-nearest-neighbor inter-orbital (diagonal, xz-yz mixing)
    """
    eps_xz = -2*t1*np.cos(kx) - 2*t2*np.cos(ky) - 4*t3*np.cos(kx)*np.cos(ky) - mu
    eps_yz = -2*t2*np.cos(kx) - 2*t1*np.cos(ky) - 4*t3*np.cos(kx)*np.cos(ky) - mu
    eps_xy_offdiag = -4*t4*np.sin(kx)*np.sin(ky)

    H = np.array([[eps_xz, eps_xy_offdiag],
                  [eps_xy_offdiag, eps_yz]])
    return H

# ==============================================================================
# Material parameter sets (QUALITATIVE, calibrated to known FS topology trends)
# ==============================================================================

# Base parameters from Raghu et al. minimal model (eV), used as reference scale
BASE_PARAMS = {
    't1': 1.0, 't2': 1.3, 't3': -0.85, 't4': -0.85, 'mu': 1.45
}

# Material-specific scaling factors (QUALITATIVE, not fitted to data)
# Rationale: larger bandwidth W ~ t_i => less correlated, larger FS pockets
# (consistent with ARPES trend: LaFeAsO/BaFe2As2/LiFeAs > FeSe in pocket size,
#  FeTe0.5Se0.5 has strong pocket reconstruction from Te/Se disorder)
MATERIAL_SCALING = {
    'LaFeAsO':      {'W_scale': 1.00, 'mu_shift': 0.00, 'Tc': 26},
    'BaFe2As2':     {'W_scale': 0.95, 'mu_shift': 0.05, 'Tc': 38},  # optimally doped
    'LiFeAs':       {'W_scale': 0.90, 'mu_shift': -0.05, 'Tc': 18},
    'FeSe':         {'W_scale': 0.70, 'mu_shift': 0.10, 'Tc': 8},   # smaller pockets (ARPES-known)
    'FeTeSe':       {'W_scale': 0.80, 'mu_shift': 0.15, 'Tc': 14},
}

def get_material_hamiltonian_params(material_name):
    """Get scaled tight-binding parameters for given material."""
    scale = MATERIAL_SCALING[material_name]
    params = {
        't1': BASE_PARAMS['t1'] * scale['W_scale'],
        't2': BASE_PARAMS['t2'] * scale['W_scale'],
        't3': BASE_PARAMS['t3'] * scale['W_scale'],
        't4': BASE_PARAMS['t4'] * scale['W_scale'],
        'mu': BASE_PARAMS['mu'] * scale['W_scale'] + scale['mu_shift'],
    }
    return params, scale['Tc']

# ==============================================================================
# Nesting function (reused from Phase 1)
# ==============================================================================

def compute_nesting_function_2orbital(q, params, fermi_window=0.05, k_grid=150):
    """
    Fermi surface nesting function for 2-orbital model.

    xi(q) = sum_k [Fermi surface indicator at k] x [Fermi surface indicator at k+q]
    """
    qx, qy = q
    t1, t2, t3, t4, mu = params['t1'], params['t2'], params['t3'], params['t4'], params['mu']

    kx_arr = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)

    nesting_count = 0

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = build_2orbital_hamiltonian(kx, ky, t1, t2, t3, t4, mu)
            H_kq = build_2orbital_hamiltonian(kx + qx, ky + qy, t1, t2, t3, t4, mu)

            eig_k = eigvalsh(H_k)
            eig_kq = eigvalsh(H_kq)

            on_fermi_k = np.abs(eig_k) < fermi_window  # mu already subtracted
            on_fermi_kq = np.abs(eig_kq) < fermi_window

            nesting_count += np.sum(on_fermi_k[:, None] & on_fermi_kq[None, :])

    return nesting_count / (k_grid**2)

def scan_nesting_bz(params, fermi_window=0.05, k_grid=150, q_grid=30, q_exclude_radius=0.15):
    """
    Scan xi(q) over full BZ to find Q_AFM.

    IMPORTANT: xi(q=0) is a TRIVIAL peak (every Fermi-surface point overlaps
    with itself), and always dominates any finite-window nesting count
    regardless of band structure. It carries no AFM-nesting information.
    We therefore exclude a small disk around q=0 (radius q_exclude_radius,
    in units of pi) from the Q_AFM search -- this is standard practice
    (e.g. Kuroki et al. 2008 identify Q_AFM from the *off-origin* peak of
    xi(q) / Im chi0(q), not from the trivial q->0 intraband response).
    """
    qx_arr = np.linspace(0, np.pi, q_grid)
    qy_arr = np.linspace(0, np.pi, q_grid)

    xi_map = np.zeros((q_grid, q_grid))

    for i, qx in enumerate(qx_arr):
        for j, qy in enumerate(qy_arr):
            xi_map[i, j] = compute_nesting_function_2orbital(
                (qx, qy), params, fermi_window, k_grid)

    # Mask out the trivial q=0 self-overlap peak before taking argmax
    q_exclude = q_exclude_radius * np.pi
    QX, QY = np.meshgrid(qx_arr, qy_arr, indexing='ij')
    mask_trivial = (QX**2 + QY**2) < q_exclude**2

    xi_map_masked = xi_map.copy()
    xi_map_masked[mask_trivial] = -np.inf

    idx_max = np.unravel_index(np.argmax(xi_map_masked), xi_map_masked.shape)
    q_afm = (qx_arr[idx_max[0]], qy_arr[idx_max[1]])
    xi_max = xi_map[idx_max]  # report the TRUE (unmasked) value at that q

    return xi_map, q_afm, xi_max, qx_arr, qy_arr

# ==============================================================================
# Main survey
# ==============================================================================

def run_survey():
    print("=" * 80)
    print("Iron-based Superconductor Geometric Nesting Survey (Phase 2)")
    print("QUALITATIVE 2-orbital model - NOT literature t_ij values")
    print("=" * 80)

    materials = ['LaFeAsO', 'BaFe2As2', 'LiFeAs', 'FeSe', 'FeTeSe']
    results = {}

    for material in materials:
        print(f"\n[{material}]")
        print("-" * 80)

        params, Tc = get_material_hamiltonian_params(material)
        print(f"  Parameters: t1={params['t1']:.3f}, t2={params['t2']:.3f}, "
              f"t3={params['t3']:.3f}, t4={params['t4']:.3f}, mu={params['mu']:.3f}")
        print(f"  Tc (experimental) = {Tc} K")

        # Coarse scan for Q_AFM (using fewer points for speed, then refine)
        xi_map, q_afm, xi_max, qx_arr, qy_arr = scan_nesting_bz(
            params, fermi_window=0.08, k_grid=100, q_grid=15)

        # Normalized nesting (relative to BZ average, EXCLUDING the trivial
        # q=0 self-overlap disk -- including it would inflate the average
        # with a physically meaningless peak and understate the true
        # AFM-channel enhancement)
        QX, QY = np.meshgrid(qx_arr, qy_arr, indexing='ij')
        q_exclude = 0.15 * np.pi
        mask_trivial = (QX**2 + QY**2) < q_exclude**2
        xi_avg = np.mean(xi_map[~mask_trivial])
        xi_normalized = xi_max / xi_avg if xi_avg > 0 else 0

        print(f"  Q_AFM = ({q_afm[0]/np.pi:.3f}pi, {q_afm[1]/np.pi:.3f}pi)")
        print(f"  xi(Q_AFM) = {xi_max:.6f}")
        print(f"  xi_avg(BZ) = {xi_avg:.6f}")
        print(f"  Normalized enhancement = {xi_normalized:.3f}x")

        results[material] = {
            'params': params,
            'Tc': Tc,
            'q_afm': q_afm,
            'xi_max': xi_max,
            'xi_avg': xi_avg,
            'xi_normalized': xi_normalized,
            'xi_map': xi_map,
            'qx_arr': qx_arr,
            'qy_arr': qy_arr
        }

    # ==============================================================================
    # Summary table
    # ==============================================================================

    print("\n" + "=" * 80)
    print("SUMMARY TABLE")
    print("=" * 80)
    print(f"{'Material':<12} {'Tc(K)':<8} {'Q_AFM (pi units)':<20} "
          f"{'xi(Q_AFM)':<12} {'Normalized':<12}")
    print("-" * 80)

    for material in materials:
        r = results[material]
        q_str = f"({r['q_afm'][0]/np.pi:.2f},{r['q_afm'][1]/np.pi:.2f})"
        print(f"{material:<12} {r['Tc']:<8} {q_str:<20} "
              f"{r['xi_max']:<12.6f} {r['xi_normalized']:<12.3f}")

    # ==============================================================================
    # Tc correlation
    # ==============================================================================

    print("\n" + "=" * 80)
    print("Tc CORRELATION ANALYSIS")
    print("=" * 80)

    Tc_array = np.array([results[m]['Tc'] for m in materials])
    xi_array = np.array([results[m]['xi_max'] for m in materials])
    xi_norm_array = np.array([results[m]['xi_normalized'] for m in materials])

    corr_raw = np.corrcoef(Tc_array, xi_array)[0, 1]
    corr_norm = np.corrcoef(Tc_array, xi_norm_array)[0, 1]

    print(f"\nPearson correlation Tc vs xi(Q_AFM): r = {corr_raw:.3f}")
    print(f"Pearson correlation Tc vs normalized xi: r = {corr_norm:.3f}")

    print("\n" + "=" * 80)
    print("CAVEATS (must be reported alongside any numbers above)")
    print("=" * 80)
    print("""
1. Tight-binding parameters are a QUALITATIVE 2-orbital model, NOT the
   literature Kuroki (2008) / Graser (2009) 5-orbital fits.
2. Material-to-material scaling (W_scale, mu_shift) is HAND-CHOSEN to
   reflect known qualitative ARPES trends (pocket size), NOT fitted to
   any specific dataset.
3. Any correlation found here is illustrative of the SURVEY METHOD,
   not a validated physical result. Must be redone with real t_ij
   before any publication claim.
4. FeSe is the only material with a validated 5-orbital model (Team B).
   Its results here (2-orbital) may differ from Phase 1 (5-orbital).
""")

    return results

if __name__ == "__main__":
    import os
    os.environ['OMP_NUM_THREADS'] = '1'
    results = run_survey()
