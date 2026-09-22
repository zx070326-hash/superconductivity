#!/usr/bin/env python3
"""
Corrected chi0_ij(q,omega) for the 5-orbital toy model, using the validated 4-U-insertion
formula from theory/notes/meta/two_band_chi0_derivation_2026-09-22.md Sec.2:

  chi_ll'(q,w) = (1/N) sum_k sum_{mu,nu} U_{l',mu}(k) U*_{l,mu}(k)
                                          U_{l,nu}(k+q) U*_{l',nu}(k+q)
                                          [f(E_mu(k)) - f(E_nu(k+q))]
                                          / (w + i eta + E_mu(k) - E_nu(k+q))

vs the OLD (buggy, single-U-rotation) formula used in
_rpa_negative_eigenvalue_survey_2026-09-22.py / _fese_5orbital_FINAL_2026-09-22.py:

  chi_orbital[i,j] += conj(U_k[i,m]) * chi0_band_mn * U_kq[j,n]

Rerun ONLY the two diagnostic points (q=0, q=(pi,0)) reported in
misc/essenberger_chi_sign_structure_2026-09-22.md Sec.3.1, on the SAME toy H(k), to
quantify the impact before touching the full BZ scan or the Graser-parameter rerun.
"""
import numpy as np
from numpy.linalg import eigh


def build_5orbital_hamiltonian(kx, ky):
    """Same toy model as _rpa_negative_eigenvalue_survey_2026-09-22.py."""
    t_xz = -0.5
    t_yz = -0.5
    t_xy = -0.3

    H = np.zeros((5, 5), dtype=complex)
    H[0, 0] = -2*t_xz*np.cos(kx)
    H[1, 1] = -2*t_yz*np.cos(ky)
    H[2, 2] = -2*t_xy*(np.cos(kx) + np.cos(ky))
    H[3, 3] = -0.2*(np.cos(kx) - np.cos(ky))
    H[4, 4] = -0.1

    t_inter = 0.1
    H[0, 1] = H[1, 0] = t_inter * np.sin(kx) * np.sin(ky)
    H[0, 2] = H[2, 0] = t_inter * np.cos(kx)
    H[1, 2] = H[2, 1] = t_inter * np.cos(ky)

    return H


def fermi_func(E, mu, T=0.01):
    x = (E - mu) / T
    return 1.0 / (1.0 + np.exp(np.clip(x, -50, 50)))


def chi0_OLD_buggy(q, mu, T=0.01, eta=0.01, k_grid=30):
    """Exact reproduction of the single-U-rotation formula (Test E's diagnosed bug)."""
    qx, qy = q
    chi0 = np.zeros((5, 5), dtype=complex)
    kx_arr = np.linspace(-np.pi, np.pi, k_grid)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid)
    dk_norm = (2*np.pi/k_grid)**2 / (2*np.pi)**2

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = build_5orbital_hamiltonian(kx, ky)
            H_kq = build_5orbital_hamiltonian(kx+qx, ky+qy)
            eig_k, U_k = eigh(H_k)
            eig_kq, U_kq = eigh(H_kq)
            f_k = fermi_func(eig_k, mu, T)
            f_kq = fermi_func(eig_kq, mu, T)
            for m in range(5):
                for n in range(5):
                    denom = eig_k[m] - eig_kq[n] + 1j*eta
                    if abs(qx) < 1e-10 and abs(qy) < 1e-10 and m == n:
                        band_mn = f_k[m]*(1-f_k[m])/T
                    else:
                        band_mn = (f_k[m]-f_kq[n])/denom
                    for i in range(5):
                        for j in range(5):
                            chi0[i, j] += np.conj(U_k[i, m]) * band_mn * U_kq[j, n] * dk_norm
    return chi0


def chi0_CORRECTED(q, mu, T=0.01, eta=0.01, k_grid=30):
    """Validated 4-U-insertion formula from two_band_chi0_derivation_2026-09-22.md."""
    qx, qy = q
    chi0 = np.zeros((5, 5), dtype=complex)
    kx_arr = np.linspace(-np.pi, np.pi, k_grid)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid)
    dk_norm = (2*np.pi/k_grid)**2 / (2*np.pi)**2

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = build_5orbital_hamiltonian(kx, ky)
            H_kq = build_5orbital_hamiltonian(kx+qx, ky+qy)
            eig_k, U_k = eigh(H_k)     # U_k[a, mu] = <a|mu,k>
            eig_kq, U_kq = eigh(H_kq)  # U_kq[a, nu] = <a|nu,k+q>
            f_k = fermi_func(eig_k, mu, T)
            f_kq = fermi_func(eig_kq, mu, T)
            for mu_ in range(5):
                for nu_ in range(5):
                    denom = eig_k[mu_] - eig_kq[nu_] + 1j*eta
                    if abs(qx) < 1e-10 and abs(qy) < 1e-10 and mu_ == nu_:
                        weight = f_k[mu_]*(1-f_k[mu_])/T
                    else:
                        weight = (f_k[mu_]-f_kq[nu_])/denom
                    # chi_{l,l'} = sum U_{l',mu}(k) U*_{l,mu}(k) U_{l,nu}(k+q) U*_{l',nu}(k+q) * weight
                    for l in range(5):
                        for lp in range(5):
                            chi0[l, lp] += (U_k[lp, mu_] * np.conj(U_k[l, mu_])
                                            * U_kq[l, nu_] * np.conj(U_kq[lp, nu_])
                                            * weight * dk_norm)
    return chi0


if __name__ == "__main__":
    mu, T, eta, k_grid = 0.0, 0.01, 0.01, 30

    for label, q in [("q=(0,0)", (0.0, 0.0)), ("q=(pi,0)", (np.pi, 0.0))]:
        print("=" * 70)
        print(label)
        print("=" * 70)

        chi_old = chi0_OLD_buggy(q, mu, T, eta, k_grid)
        chi_new = chi0_CORRECTED(q, mu, T, eta, k_grid)

        eig_old = np.linalg.eigvalsh(chi_old.real)
        eig_new = np.linalg.eigvalsh(chi_new.real)

        print(f"OLD (buggy, single-U-rotation) eigenvalues: {np.round(eig_old, 4)}")
        print(f"NEW (corrected, 4-U-insertion) eigenvalues: {np.round(eig_new, 4)}")
        print(f"  # negative (OLD): {np.sum(eig_old < -1e-6)}")
        print(f"  # negative (NEW): {np.sum(eig_new < -1e-6)}")
        print(f"  max abs elementwise diff: {np.max(np.abs(chi_old - chi_new)):.4e}")
        print(f"  diagonal OLD: {np.round(np.diag(chi_old).real, 4)}")
        print(f"  diagonal NEW: {np.round(np.diag(chi_new).real, 4)}")
        print()
