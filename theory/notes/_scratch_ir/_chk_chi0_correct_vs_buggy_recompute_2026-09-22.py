"""Recompute Q_AFM and q=0 5x5 chi0/chi_s eigenvalues with the CORRECT
(4-U-insertion, cross-term) formula, using the IDENTICAL toy Hamiltonian
and parameters as the buggy production script, for a direct point-by-point
comparison of "which conclusions change, which stay the same."

Reference (buggy) script: _fese_5orbital_FINAL_2026-09-22.py
  chi0_orbital[i,j] += conj(U_k[i,m]) * chi0_band_mn * U_kq[j,n]   (2 U insertions)

Correct formula (derived + validated on 2-band model in
two_band_chi0_derivation_2026-09-22.md / two_band_chi0_verification_2026-09-22.py,
also independently re-derived by team-lead in _chk_chi0_gauge_invariance_2026-09-22.py):
  chi0_orbital[i,j] += U_k[j,m]*conj(U_k[i,m]) * chi0_band_mn * U_kq[i,n]*conj(U_kq[j,n])
                                                                    (4 U insertions)

Same Hamiltonian, same k-grid, same mu/T/eta/U_spin as the production script,
so any eigenvalue-count difference is attributable ONLY to the formula bug.
"""

import numpy as np
from numpy.linalg import eigh, inv, eigvalsh

NORB = 5


def build_5orbital_hamiltonian(kx, ky):
    """Identical to _fese_5orbital_FINAL_2026-09-22.py -- copied verbatim
    so this script has no import dependency on the (buggy) production file."""
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


def chi0_5orbital(q, mu, T=0.01, eta=0.01, k_grid=30, mode="correct"):
    """mode='buggy'   -> reproduces the exact 2-U-insertion formula (for reference/contrast)
       mode='correct' -> the validated 4-U-insertion cross-term formula
    """
    assert mode in ("buggy", "correct")
    qx, qy = q
    chi0_orbital = np.zeros((NORB, NORB), dtype=complex)

    kx_arr = np.linspace(-np.pi, np.pi, k_grid)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid)
    dk = (2*np.pi / k_grid)**2
    BZ_area = (2*np.pi)**2
    dk_normalized = dk / BZ_area

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = build_5orbital_hamiltonian(kx, ky)
            H_kq = build_5orbital_hamiltonian(kx + qx, ky + qy)

            eig_k, U_k = eigh(H_k)
            eig_kq, U_kq = eigh(H_kq)

            f_k = np.array([fermi_func(e, mu, T) for e in eig_k])
            f_kq = np.array([fermi_func(e, mu, T) for e in eig_kq])

            for m in range(NORB):
                for n in range(NORB):
                    denom = eig_k[m] - eig_kq[n] + 1j*eta
                    if abs(qx) < 1e-10 and abs(qy) < 1e-10 and m == n:
                        chi0_band_mn = f_k[m] * (1 - f_k[m]) / T
                    else:
                        chi0_band_mn = (f_k[m] - f_kq[n]) / denom

                    if mode == "buggy":
                        for i in range(NORB):
                            for j in range(NORB):
                                chi0_orbital[i, j] += (
                                    np.conj(U_k[i, m]) * chi0_band_mn * U_kq[j, n]
                                    * dk_normalized
                                )
                    else:
                        # cross term: factor1[i,j]=U_k[j,m]*conj(U_k[i,m]);
                        #             factor2[i,j]=U_kq[i,n]*conj(U_kq[j,n])
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


def n_neg(mat):
    ev = eigvalsh((mat + mat.conj().T).real / 2)
    return int((ev < -1e-6).sum()), ev


def report(q, label, mu, T, eta, k_grid, U_spin):
    print("=" * 70)
    print(f"q = {label} = {q}")
    print("=" * 70)

    chi0_bug = chi0_5orbital(q, mu, T=T, eta=eta, k_grid=k_grid, mode="buggy")
    chi0_cor = chi0_5orbital(q, mu, T=T, eta=eta, k_grid=k_grid, mode="correct")

    nb, evb = n_neg(chi0_bug)
    nc, evc = n_neg(chi0_cor)
    print(f"chi0 (buggy,   2-U): n_neg = {nb}, eigs = {np.round(evb,4)}")
    print(f"chi0 (correct, 4-U): n_neg = {nc}, eigs = {np.round(evc,4)}")

    diag_diff = np.abs(np.diag(chi0_bug) - np.diag(chi0_cor)).max()
    print(f"max |diag(buggy)-diag(correct)| = {diag_diff:.6e}  "
          f"(nonzero here means DIAGONAL is also corrupted by the bug)")

    chis_bug = rpa_spin_susceptibility(chi0_bug, U_spin)
    chis_cor = rpa_spin_susceptibility(chi0_cor, U_spin)

    if chis_bug is not None:
        nbs, evbs = n_neg(chis_bug)
        print(f"chi_s(RPA) (buggy,   2-U): n_neg = {nbs}, eigs = {np.round(evbs,4)}")
    else:
        nbs = None
        print("chi_s(RPA) (buggy): RPA inversion failed")

    if chis_cor is not None:
        ncs, evcs = n_neg(chis_cor)
        print(f"chi_s(RPA) (correct, 4-U): n_neg = {ncs}, eigs = {np.round(evcs,4)}")
    else:
        ncs = None
        print("chi_s(RPA) (correct): RPA inversion failed")

    print()
    return dict(q=label, n_neg_chi0_buggy=nb, n_neg_chi0_correct=nc,
                n_neg_chis_buggy=nbs, n_neg_chis_correct=ncs,
                diag_diff=diag_diff)


def main():
    mu = 0.0
    U_spin = 0.5
    T = 0.01
    eta = 0.01
    k_grid = 30

    results = []
    results.append(report((0.0, 0.0), "Gamma (0,0)", mu, T, eta, k_grid, U_spin))
    results.append(report((np.pi, 0.0), "Q_AFM (pi,0)", mu, T, eta, k_grid, U_spin))

    print("=" * 70)
    print("SUMMARY TABLE (n_neg = number of eigenvalues < -1e-6)")
    print("=" * 70)
    print(f"{'q':16s} {'chi0 buggy':>12s} {'chi0 correct':>13s} "
          f"{'chi_s buggy':>12s} {'chi_s correct':>14s} {'diag corrupted?':>16s}")
    for r in results:
        print(f"{r['q']:16s} {r['n_neg_chi0_buggy']:>12d} {r['n_neg_chi0_correct']:>13d} "
              f"{str(r['n_neg_chis_buggy']):>12s} {str(r['n_neg_chis_correct']):>14s} "
              f"{'YES' if r['diag_diff']>1e-8 else 'no':>16s}")


if __name__ == "__main__":
    main()
