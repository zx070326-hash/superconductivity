#!/usr/bin/env python3
"""
Full-BZ negative-eigenvalue survey, REAL Graser et al. arXiv:1003.0133 BaFe2As2
5-orbital H(k), using the CORRECT (4-U-insertion, cross-term) chi0 formula --
i.e. the same real-material rerun as
`_rpa_negative_eigenvalue_survey_GRASER_2026-09-22.py`, but with the bug fixed.

That script's own chi0 projection step (`np.conj(U_k[i,m])*chi0_band_mn*U_kq[j,n]`)
carries the identical 2-U-insertion bug confirmed in
`_chk_chi0_correct_vs_buggy_recompute_2026-09-22.py` and
`_chk_chi0_gauge_invariance_2026-09-22.py` (team-lead, independent derivation).
This script reruns the SAME survey (same H(k), same mu/T/eta/k_grid/q_grid/U_spin)
with the correct formula, for a direct "which conclusions survive" comparison.
"""

import numpy as np
from numpy.linalg import eigh, inv
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


def lindhard_chi0_5orbital_correct(q, mu, T=0.01, eta=0.01, k_grid=30):
    """Correct 4-U-insertion cross-term formula (see two_band_chi0_derivation_2026-09-22.md)."""
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

            f_k = fermi_func(eig_k, mu, T)
            f_kq = fermi_func(eig_kq, mu, T)

            for m in range(5):
                for n in range(5):
                    denom = eig_k[m] - eig_kq[n] + 1j*eta
                    if abs(qx) < 1e-10 and abs(qy) < 1e-10 and m == n:
                        chi0_band_mn = f_k[m] * (1 - f_k[m]) / T
                    else:
                        chi0_band_mn = (f_k[m] - f_kq[n]) / denom

                    factor1 = np.conj(U_k[:, m])[:, None] * U_k[:, m][None, :]
                    factor2 = U_kq[:, n][:, None] * np.conj(U_kq[:, n])[None, :]
                    chi0_orbital += factor1 * factor2 * chi0_band_mn * dk_normalized

    return chi0_orbital


def rpa_spin_susceptibility(chi0, U_spin):
    I = np.eye(chi0.shape[0])
    try:
        return chi0 @ inv(I - U_spin * chi0)
    except np.linalg.LinAlgError:
        return None


def main():
    print("=" * 70)
    print("REAL Graser-parameter (BaFe2As2) BZ survey -- CORRECT (4-U) formula")
    print("=" * 70)

    mu = -0.0656
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
            if (i * q_grid + j) % 20 == 0:
                print(f"  Progress: {i*q_grid + j}/{q_grid**2}", end='\r')

            chi0 = lindhard_chi0_5orbital_correct((qx, qy), mu, T=T, k_grid=k_grid_chi)
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

    total_points = q_grid**2
    chi0_has_neg = np.sum(chi0_neg_eigs > 0)
    chi_s_has_neg = np.sum(chi_s_neg_eigs > 0)

    print("\n" + "=" * 70)
    print("Results (REAL Graser BaFe2As2 parameters, CORRECT formula)")
    print("=" * 70)
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

    np.savez(os.path.join(_here, "_rpa_negative_eigenvalue_survey_GRASER_CORRECT_2026-09-22.npz"),
              chi0_neg_eigs=chi0_neg_eigs, chi_s_neg_eigs=chi_s_neg_eigs,
              chi0_min_eig=chi0_min_eig, chi_s_min_eig=chi_s_min_eig)
    print("\nSaved raw arrays to _rpa_negative_eigenvalue_survey_GRASER_CORRECT_2026-09-22.npz")


if __name__ == "__main__":
    main()
