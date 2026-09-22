#!/usr/bin/env python3
"""
Bare RPA negative-eigenvalue full-BZ survey, REAL Graser et al. arXiv:1003.0133
BaFe2As2 5-orbital tight-binding parameters (kz=0 cut), replacing the toy-model
Hamiltonian used in the 2026-09-22 baseline survey.

Purpose: check whether the toy-model finding (100% of BZ has negative chi0/chi_s
eigenvalues, RPA 0% improvement) is an artifact of the placeholder H(k), or persists
under a real DFT-fitted band structure.

mu fixed at -0.0656 eV (see _graser_find_mu_fermi_surface_2026-09-22.py: T=0 filling
frac 6/10 electrons/Fe, exact fraction below mu = 0.6000 on a 120x120 grid).

Caveat (must propagate to any note using this script as evidence): material is
BaFe2As2 (122-family), not FeSe (11-family) itself -- same qualitative Fermi-surface
topology (per Graser Fig 1-2 cross-check) but quantitatively different bands.
"""

import numpy as np
from numpy.linalg import eigh, inv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import importlib.util, os

_here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location(
    "graser_H", os.path.join(_here, "_graser1003_hamiltonian_2026-09-22.py"))
graser_H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(graser_H)
build_H_graser = graser_H.build_H_graser


def fermi_func(E, mu, T=0.01):
    x = (E - mu) / T
    return 1.0 / (1.0 + np.exp(np.clip(x, -50, 50)))


def lindhard_chi0_5orbital(q, mu, T=0.01, eta=0.01, k_grid=30):
    qx, qy = q
    chi0_orbital = np.zeros((5, 5), dtype=complex)

    kx_arr = np.linspace(-np.pi, np.pi, k_grid)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid)
    dk = (2*np.pi / k_grid)**2
    BZ_area = (2*np.pi)**2
    dk_normalized = dk / BZ_area

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = build_H_graser(kx, ky)
            H_kq = build_H_graser(kx + qx, ky + qy)

            eig_k, U_k = eigh(H_k)
            eig_kq, U_kq = eigh(H_kq)

            f_k = np.array([fermi_func(e, mu, T) for e in eig_k])
            f_kq = np.array([fermi_func(e, mu, T) for e in eig_kq])

            for m in range(5):
                for n in range(5):
                    denom = eig_k[m] - eig_kq[n] + 1j*eta

                    if abs(qx) < 1e-10 and abs(qy) < 1e-10 and m == n:
                        chi0_band_mn = f_k[m] * (1 - f_k[m]) / T
                    else:
                        chi0_band_mn = (f_k[m] - f_kq[n]) / denom

                    for i in range(5):
                        for j in range(5):
                            chi0_orbital[i, j] += (
                                np.conj(U_k[i, m]) * chi0_band_mn * U_kq[j, n]
                                * dk_normalized
                            )

    return chi0_orbital


def rpa_spin_susceptibility(chi0, U_spin):
    I = np.eye(chi0.shape[0])
    try:
        chi_s = chi0 @ inv(I - U_spin * chi0)
        return chi_s
    except np.linalg.LinAlgError:
        return None


print("=" * 70)
print("REAL Graser-parameter (BaFe2As2) bare RPA negative-eigenvalue BZ survey")
print("=" * 70)

mu = -0.0656  # from _graser_find_mu_fermi_surface_2026-09-22.py, filling=0.6000 exact
U_spin = 0.5
T = 0.01
k_grid_chi = 30

q_grid = 16
qx_arr = np.linspace(-np.pi, np.pi, q_grid)
qy_arr = np.linspace(-np.pi, np.pi, q_grid)

chi0_neg_eigs = np.zeros((q_grid, q_grid), dtype=int)
chi_s_neg_eigs = np.zeros((q_grid, q_grid), dtype=int)
chi0_min_eig = np.zeros((q_grid, q_grid))
chi_s_min_eig = np.zeros((q_grid, q_grid))

print(f"\nScanning {q_grid}x{q_grid} = {q_grid**2} q-points (k-grid {k_grid_chi}x{k_grid_chi})...")

for i, qx in enumerate(qx_arr):
    for j, qy in enumerate(qy_arr):
        if (i * q_grid + j) % 50 == 0:
            print(f"  Progress: {i*q_grid + j}/{q_grid**2}", end='\r')

        chi0 = lindhard_chi0_5orbital((qx, qy), mu, T=T, k_grid=k_grid_chi)
        eigs0 = np.linalg.eigvalsh(chi0.real)

        chi_s = rpa_spin_susceptibility(chi0, U_spin)
        if chi_s is not None:
            eigs_s = np.linalg.eigvalsh(chi_s.real)
        else:
            eigs_s = np.array([np.nan] * 5)

        chi0_neg_eigs[i, j] = np.sum(eigs0 < -1e-6)
        chi_s_neg_eigs[i, j] = np.sum(eigs_s < -1e-6)
        chi0_min_eig[i, j] = np.min(eigs0)
        chi_s_min_eig[i, j] = np.min(eigs_s) if not np.isnan(eigs_s[0]) else np.nan

print(f"\nScan complete: {q_grid**2} q-points")

print("\n" + "=" * 70)
print("Results (REAL Graser BaFe2As2 parameters)")
print("=" * 70)

total_points = q_grid**2
chi0_has_neg = np.sum(chi0_neg_eigs > 0)
chi_s_has_neg = np.sum(chi_s_neg_eigs > 0)

print(f"\nchi0 negative eigenvalues:")
print(f"  q-points with negative eigenvalues: {chi0_has_neg}/{total_points} ({100*chi0_has_neg/total_points:.1f}%)")
print(f"  Min eigenvalue range: [{np.min(chi0_min_eig):.4f}, {np.max(chi0_min_eig):.4f}]")

print(f"\nchi_s (RPA) negative eigenvalues:")
print(f"  q-points with negative eigenvalues: {chi_s_has_neg}/{total_points} ({100*chi_s_has_neg/total_points:.1f}%)")
print(f"  Min eigenvalue range: [{np.nanmin(chi_s_min_eig):.4f}, {np.nanmax(chi_s_min_eig):.4f}]")

worsened = np.sum(chi_s_neg_eigs > chi0_neg_eigs)
improved = np.sum(chi_s_neg_eigs < chi0_neg_eigs)
unchanged = np.sum(chi_s_neg_eigs == chi0_neg_eigs)

print(f"\nRPA effect:")
print(f"  Worsened: {worsened}/{total_points} ({100*worsened/total_points:.1f}%)")
print(f"  Improved: {improved}/{total_points} ({100*improved/total_points:.1f}%)")
print(f"  Unchanged: {unchanged}/{total_points} ({100*unchanged/total_points:.1f}%)")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

im0 = axes[0, 0].imshow(chi0_neg_eigs.T, origin='lower', cmap='Reds',
                         extent=[-np.pi, np.pi, -np.pi, np.pi], vmin=0, vmax=5)
axes[0, 0].set_title(r'$\chi_0$ (Graser): # negative eigenvalues', fontsize=12)
plt.colorbar(im0, ax=axes[0, 0])

im1 = axes[0, 1].imshow(chi_s_neg_eigs.T, origin='lower', cmap='Reds',
                         extent=[-np.pi, np.pi, -np.pi, np.pi], vmin=0, vmax=5)
axes[0, 1].set_title(r'$\chi^s$ (Graser, RPA): # negative eigenvalues', fontsize=12)
plt.colorbar(im1, ax=axes[0, 1])

im2 = axes[1, 0].imshow(chi0_min_eig.T, origin='lower', cmap='RdBu_r',
                         extent=[-np.pi, np.pi, -np.pi, np.pi])
axes[1, 0].set_title(r'$\chi_0$ (Graser): min eigenvalue', fontsize=12)
plt.colorbar(im2, ax=axes[1, 0])

im3 = axes[1, 1].imshow(chi_s_min_eig.T, origin='lower', cmap='RdBu_r',
                         extent=[-np.pi, np.pi, -np.pi, np.pi])
axes[1, 1].set_title(r'$\chi^s$ (Graser, RPA): min eigenvalue', fontsize=12)
plt.colorbar(im3, ax=axes[1, 1])

plt.tight_layout()
outpng = os.path.join(_here, '_rpa_negative_eigenvalue_survey_GRASER_2026-09-22.png')
plt.savefig(outpng, dpi=150)
print(f"\nHeatmap saved: {outpng}")
