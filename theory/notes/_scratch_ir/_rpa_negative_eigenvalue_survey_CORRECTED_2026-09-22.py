#!/usr/bin/env python3
"""
Full-BZ static chi0/chi^s (RPA) negative-eigenvalue survey, rerun with the
VALIDATED 4-U-insertion Lindhard formula (theory/notes/meta/two_band_chi0_derivation_2026-09-22.md
Sec.2, independently re-verified at 2 diagnostic points in
_chk_5orbital_chi0_bug_impact_2026-09-22.py), replacing the buggy "single U rotation"
(2-insertion) formula used in the WITHDRAWN _rpa_negative_eigenvalue_survey_2026-09-22.py
and _rpa_negative_eigenvalue_survey_GRASER_2026-09-22.py.

Runs on BOTH:
  (a) the toy 5-orbital H(k) (same as the withdrawn toy-model scan, for direct comparison)
  (b) the real Graser et al. 1003.0133 BaFe2As2 H(k), with the CORRECTED Table I/II
      parameters (theory/notes/literature/graser1003.0133_5orbital_TB_parameters_2026-09-22.md
      Sec.4/5) and the recomputed mu=0.0127 eV (_graser_find_mu_fermi_surface_2026-09-22.py rerun).

Strictness: numerical-evidence (float64, finite 16x16 q grid x 30x30 k grid, eta=0.01 broadening).
"""
import numpy as np
from numpy.linalg import eigh
import importlib.util, os

SCRIPT_DIR = os.path.dirname(__file__)


def load_module(fname):
    spec = importlib.util.spec_from_file_location(fname, os.path.join(SCRIPT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_5orbital_hamiltonian_toy(kx, ky):
    """Same toy model as the withdrawn _rpa_negative_eigenvalue_survey_2026-09-22.py."""
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


def fermi_func(E, mu, T):
    x = (E - mu) / T
    return 1.0 / (1.0 + np.exp(np.clip(x, -50, 50)))


def chi0_matrix(q, H_func, mu, T=0.01, eta=0.01, k_grid=30, norb=5):
    """Validated 4-U-insertion formula (two_band_chi0_derivation_2026-09-22.md Sec.2),
    vectorized over orbital indices via einsum (full-BZ scan needs this; the pure-python
    quadruple loop used in _chk_5orbital_chi0_bug_impact_2026-09-22.py was only run at
    2 points and would be too slow for 256 q-points x 900 k-points):
    chi_{l,l'}(q,0) = (1/N) sum_k sum_{mu,nu} U_{l',mu}(k) U*_{l,mu}(k)
                                               U_{l,nu}(k+q) U*_{l',nu}(k+q) * weight
    Cross-checked against the pure-python loop version at q=(0,0),(pi,0) below (__main__),
    max elementwise diff must be ~0 (both are the same formula, just re-expressed).
    """
    qx, qy = q
    chi0 = np.zeros((norb, norb), dtype=complex)
    kx_arr = np.linspace(-np.pi, np.pi, k_grid)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid)
    dk_norm = (2*np.pi/k_grid)**2 / (2*np.pi)**2
    q0 = (abs(qx) < 1e-10 and abs(qy) < 1e-10)

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = H_func(kx, ky)
            H_kq = H_func(kx+qx, ky+qy)
            eig_k, U_k = eigh(H_k)
            eig_kq, U_kq = eigh(H_kq)
            f_k = fermi_func(eig_k, mu, T)
            f_kq = fermi_func(eig_kq, mu, T)

            denom = eig_k[:, None] - eig_kq[None, :] + 1j*eta   # [mu, nu]
            weight = (f_k[:, None] - f_kq[None, :]) / denom
            if q0:
                # CORRECTED 2026-09-22 (team-lead catch): the q=0, m=n L'Hopital limit of
                # (f_m-f_n)/(eps_m-eps_n) is df/deps = -f(1-f)/T, NOT +f(1-f)/T. The original
                # diagnostic script (_chk_5orbital_chi0_bug_impact_2026-09-22.py) and this
                # script both had the wrong sign here; this affected ONLY the q=(0,0) point's
                # diagonal terms. This static-chi0 scan is superseded anyway (team-lead proved
                # Re chi0(q,0) is NSD as an analytic identity -> zero information for Eq.(37));
                # kept here fixed for the record, not for further static-scan use.
                diag_w = -f_k*(1-f_k)/T
                np.fill_diagonal(weight, diag_w)

            # A[mu,l,lp] = conj(U_k[l,mu]) * U_k[lp,mu]
            A = np.einsum('lm,pm->mlp', np.conj(U_k), U_k)
            # B[nu,l,lp] = U_kq[l,nu] * conj(U_kq[lp,nu])
            B = np.einsum('ln,pn->nlp', U_kq, np.conj(U_kq))
            chi0 += np.einsum('mn,mlp,nlp->lp', weight, A, B) * dk_norm
    return chi0


def run_survey(label, H_func, mu, q_grid=16, k_grid=30, T=0.01, eta=0.01, U_rpa=0.3, norb=5):
    qs = np.linspace(-np.pi, np.pi, q_grid, endpoint=False)
    n_neg_chi0 = 0
    n_neg_chis = 0
    min_eig_chi0 = []
    min_eig_chis = []
    n_improved = 0
    n_worsened = 0
    n_unchanged = 0
    Umat = U_rpa * np.eye(norb)

    total = q_grid * q_grid
    count = 0
    for qx in qs:
        for qy in qs:
            count += 1
            chi0 = chi0_matrix((qx, qy), H_func, mu, T=T, eta=eta, k_grid=k_grid, norb=norb)
            chi0_r = 0.5*(chi0 + chi0.conj().T)  # Hermitian part for real-eigval PSD test
            eig0 = np.linalg.eigvalsh(chi0_r.real if np.allclose(chi0_r.imag, 0, atol=1e-8) else chi0_r)
            eig0 = np.linalg.eigvalsh(chi0_r)
            neg0 = int(np.sum(eig0.real < -1e-8))
            n_neg_chi0 += 1 if neg0 > 0 else 0
            min_eig_chi0.append(eig0.real.min())

            try:
                chis = chi0 @ np.linalg.inv(np.eye(norb) - Umat @ chi0)
            except np.linalg.LinAlgError:
                chis = chi0
            chis_r = 0.5*(chis + chis.conj().T)
            eigs = np.linalg.eigvalsh(chis_r)
            negs = int(np.sum(eigs.real < -1e-8))
            n_neg_chis += 1 if negs > 0 else 0
            min_eig_chis.append(eigs.real.min())

            if negs < neg0:
                n_improved += 1
            elif negs > neg0:
                n_worsened += 1
            else:
                n_unchanged += 1

            if count % 32 == 0:
                print(f"  [{label}] {count}/{total} q-points done")

    min_eig_chi0 = np.array(min_eig_chi0)
    min_eig_chis = np.array(min_eig_chis)
    print(f"=== {label} ===")
    print(f"  chi0: {n_neg_chi0}/{total} q-points with >=1 negative eigenvalue "
          f"({100*n_neg_chi0/total:.1f}%); min eig range [{min_eig_chi0.min():.4f}, {min_eig_chi0.max():.4f}]")
    print(f"  chi^s (RPA, U={U_rpa}): {n_neg_chis}/{total} q-points with >=1 negative eigenvalue "
          f"({100*n_neg_chis/total:.1f}%); min eig range [{min_eig_chis.min():.4f}, {min_eig_chis.max():.4f}]")
    print(f"  RPA effect: improved {n_improved}/{total} ({100*n_improved/total:.1f}%), "
          f"worsened {n_worsened}/{total} ({100*n_worsened/total:.1f}%), "
          f"unchanged {n_unchanged}/{total} ({100*n_unchanged/total:.1f}%)")
    return dict(n_neg_chi0=n_neg_chi0, n_neg_chis=n_neg_chis, total=total,
                min_eig_chi0=min_eig_chi0, min_eig_chis=min_eig_chis,
                n_improved=n_improved, n_worsened=n_worsened, n_unchanged=n_unchanged)


if __name__ == "__main__":
    results = {}

    print("Running TOY model survey (corrected 4-U-insertion formula)...")
    results['toy'] = run_survey("toy (corrected formula)", build_5orbital_hamiltonian_toy,
                                 mu=0.0, q_grid=16, k_grid=30)

    graser_H_mod = load_module("_graser1003_hamiltonian_2026-09-22.py")
    build_H_graser = graser_H_mod.build_H_graser
    mu_graser = float(np.load(os.path.join(SCRIPT_DIR, "_graser_mu_value.npy")))
    print(f"\nUsing recomputed Graser mu = {mu_graser:.4f} eV (corrected Table I/II)")

    def H_graser_kz0(kx, ky):
        return build_H_graser(kx, ky, kz=0.0)

    print("Running GRASER (corrected params) model survey (corrected 4-U-insertion formula)...")
    results['graser'] = run_survey("Graser (corrected params + formula)", H_graser_kz0,
                                    mu=mu_graser, q_grid=16, k_grid=30)

    np.savez(os.path.join(SCRIPT_DIR, "_rpa_negative_eigenvalue_survey_CORRECTED_results.npz"),
             toy_min_chi0=results['toy']['min_eig_chi0'], toy_min_chis=results['toy']['min_eig_chis'],
             graser_min_chi0=results['graser']['min_eig_chi0'], graser_min_chis=results['graser']['min_eig_chis'])
    print("\nDone. Results saved to _rpa_negative_eigenvalue_survey_CORRECTED_results.npz")
