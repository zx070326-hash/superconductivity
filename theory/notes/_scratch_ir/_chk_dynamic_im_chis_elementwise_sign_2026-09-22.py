#!/usr/bin/env python3
"""
Im chi^s_{ll'}(q,omega>0) elementwise sign scan -- the actual open question per
team-lead's 2026-09-22 assignment (static Re chi(q,0) non-PSD is a proven identity,
zero information for Essenberger Eq.(37); Eq.(37) needs elementwise Im Lambda^SF>0,
which traces back through Eq.(36) to Im chi^s(q,omega>0), NOT the static real part).

Uses the gauge-invariant 4-U-insertion formula (validated at machine precision in
_chk_gauge_invariance_gate_dynamic_2026-09-22.py: max diff 1.1e-16 over 4 (q,omega)
points x 8 random complex rephasings each -- GATE PASS, run first, required before
trusting any number here).

RPA dressing: chi^s = chi0 (I - U chi0)^-1, U = U_rpa * I (same simple on-site
Hubbard U as the (now-retracted-as-zero-information) static survey, U_rpa=0.3, kept
for continuity/comparability, not re-derived here).

Reports, per (q,omega) grid point: whether Im chi^s_{ll'} is (a) uniformly one sign
across all l,l' (the FDT-consistent, "boring" case) or (b) mixed sign across
elements (the case that would matter for Eq.(37)'s elementwise-positivity claim).
Diagonal and off-diagonal are reported separately since Eq.(37)'s novel content is
specifically about off-diagonal (inter-orbital) elements.

Strictness: numerical-evidence (float64, finite q/omega/k grids, eta=0.02 broadening).
"""
import numpy as np
import importlib.util, os

SCRIPT_DIR = os.path.dirname(__file__)


def load_module(fname):
    spec = importlib.util.spec_from_file_location(fname, os.path.join(SCRIPT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gate_mod = load_module("_chk_gauge_invariance_gate_dynamic_2026-09-22.py")
chi0_dynamic = gate_mod.chi0_dynamic


def run_dynamic_scan(H_func, mu, q_grid=8, omegas=(0.05, 0.1, 0.2, 0.5, 1.0),
                      T=0.02, eta=0.02, k_grid=20, U_rpa=0.3, norb=5):
    qs = np.linspace(-np.pi, np.pi, q_grid, endpoint=False)
    Umat = U_rpa * np.eye(norb)

    n_diag_neg = 0
    n_diag_pos = 0
    n_offdiag_neg = 0
    n_offdiag_pos = 0
    n_mixed_sign_points = 0
    n_total_points = 0
    n_herm_indefinite = 0  # Hermitian part of Im chi^s has mixed-sign eigenvalues

    for omega in omegas:
        for qx in qs:
            for qy in qs:
                n_total_points += 1
                chi0 = chi0_dynamic((qx, qy), omega, H_func, mu, T=T, eta=eta,
                                     k_grid=k_grid, norb=norb, rephase=False)
                try:
                    chis = chi0 @ np.linalg.inv(np.eye(norb) - Umat @ chi0)
                except np.linalg.LinAlgError:
                    continue
                im_chis = chis.imag

                diag = np.diag(im_chis)
                offdiag = im_chis[~np.eye(norb, dtype=bool)]
                n_diag_neg += int(np.sum(diag < -1e-10))
                n_diag_pos += int(np.sum(diag > 1e-10))
                n_offdiag_neg += int(np.sum(offdiag < -1e-10))
                n_offdiag_pos += int(np.sum(offdiag > 1e-10))

                all_elems = np.concatenate([diag, offdiag])
                has_pos = np.any(all_elems > 1e-10)
                has_neg = np.any(all_elems < -1e-10)
                if has_pos and has_neg:
                    n_mixed_sign_points += 1

                herm_im = 0.5*(im_chis + im_chis.T)  # symmetric part, real eigvals well-defined
                eig_herm = np.linalg.eigvalsh(herm_im)
                if np.any(eig_herm > 1e-10) and np.any(eig_herm < -1e-10):
                    n_herm_indefinite += 1

    print(f"Scanned {n_total_points} (q,omega) points, omega in {omegas}, q_grid={q_grid}x{q_grid}, k_grid={k_grid}")
    print(f"Diagonal Im chi^s_ll: {n_diag_neg} negative-element-instances, {n_diag_pos} positive")
    print(f"Off-diagonal Im chi^s_ll' (l!=l'): {n_offdiag_neg} negative-element-instances, {n_offdiag_pos} positive")
    print(f"(q,omega) points with MIXED sign across all matrix elements: {n_mixed_sign_points}/{n_total_points}")
    print(f"(q,omega) points where Hermitian(Im chi^s) has indefinite eigenvalues (violates the "
          f"'consistent sign, not literal PSD' FDT check): {n_herm_indefinite}/{n_total_points}")
    return dict(n_diag_neg=n_diag_neg, n_diag_pos=n_diag_pos,
                n_offdiag_neg=n_offdiag_neg, n_offdiag_pos=n_offdiag_pos,
                n_mixed_sign_points=n_mixed_sign_points, n_total_points=n_total_points,
                n_herm_indefinite=n_herm_indefinite)


if __name__ == "__main__":
    graser_H_mod = load_module("_graser1003_hamiltonian_2026-09-22.py")
    build_H_graser = graser_H_mod.build_H_graser
    mu_graser = float(np.load(os.path.join(SCRIPT_DIR, "_graser_mu_value.npy")))

    def H_graser_kz0(kx, ky):
        return build_H_graser(kx, ky, kz=0.0)

    print("=== Graser BaFe2As2 (corrected params), Im chi^s(q,omega>0) elementwise sign scan ===")
    result = run_dynamic_scan(H_graser_kz0, mu_graser, q_grid=8,
                               omegas=(0.05, 0.1, 0.2, 0.5, 1.0), T=0.02, eta=0.02, k_grid=20)

    np.savez(os.path.join(SCRIPT_DIR, "_dynamic_im_chis_sign_scan_results.npz"), **result)
    print("\nDone.")
