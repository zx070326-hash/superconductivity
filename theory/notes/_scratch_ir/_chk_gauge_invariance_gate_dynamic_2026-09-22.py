#!/usr/bin/env python3
"""
Gauge-invariance acceptance gate for the validated 4-U-insertion chi0 formula
(theory/notes/meta/two_band_chi0_derivation_2026-09-22.md Sec.2), extended to
general (q, omega) -- needed before trusting any Im chi^s(q,omega>0) elementwise
sign scan (team-lead's assignment, 2026-09-22).

Team-lead's finding (chi0_two_factor_form_is_gauge_dependent_2026-09-22.md): the OLD
"single U rotation" formula is gauge-dependent (LAPACK-phase-convention-dependent).
The four-factor / cross-term formula is claimed gauge-invariant by construction
(each orbital index l appears once conjugated and once not, for both k and k+q,
so a per-band phase e^{i theta_mu(k)} cancels exactly against its conjugate).
This script actually TESTS that claim numerically, at nonzero q and omega, with
COMPLEX per-band rephasing (U(1), not just +-1 -- Graser H(k) is complex Hermitian,
so eigenvectors from eigh are complex and the residual gauge freedom is a full
phase, not just a sign).

Strictness: numerical-evidence (float64), but the property being tested (gauge
invariance) is an algebraic identity -- a single violation found here would be a
genuine falsification of the formula, not evidence needing large-sample stats.
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


def fermi_func(E, mu, T):
    x = (E - mu) / T
    return 1.0 / (1.0 + np.exp(np.clip(x, -50, 50)))


def chi0_dynamic(q, omega, H_func, mu, T=0.01, eta=0.02, k_grid=20, norb=5,
                  rng=None, rephase=False):
    """4-U-insertion chi0(q, omega), optionally with a random per-(k,band) complex
    rephasing applied to U_k and U_kq independently, to test gauge invariance.
    chi_{l,l'}(q,w) = (1/N) sum_k sum_{mu,nu} U_{l',mu}(k) U*_{l,mu}(k)
                                               U_{l,nu}(k+q) U*_{l',nu}(k+q)
                                               [f(E_mu(k))-f(E_nu(k+q))]
                                               / (w+i*eta+E_mu(k)-E_nu(k+q))
    """
    qx, qy = q
    chi0 = np.zeros((norb, norb), dtype=complex)
    kx_arr = np.linspace(-np.pi, np.pi, k_grid)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid)
    dk_norm = (2*np.pi/k_grid)**2 / (2*np.pi)**2

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = H_func(kx, ky)
            H_kq = H_func(kx+qx, ky+qy)
            eig_k, U_k = eigh(H_k)
            eig_kq, U_kq = eigh(H_kq)

            if rephase:
                phase_k = np.exp(1j*rng.uniform(0, 2*np.pi, size=norb))
                phase_kq = np.exp(1j*rng.uniform(0, 2*np.pi, size=norb))
                U_k = U_k * phase_k[None, :]
                U_kq = U_kq * phase_kq[None, :]

            f_k = fermi_func(eig_k, mu, T)
            f_kq = fermi_func(eig_kq, mu, T)

            denom = eig_k[:, None] - eig_kq[None, :] + omega + 1j*eta
            weight = (f_k[:, None] - f_kq[None, :]) / denom

            A = np.einsum('lm,pm->mlp', np.conj(U_k), U_k)
            B = np.einsum('ln,pn->nlp', U_kq, np.conj(U_kq))
            chi0 += np.einsum('mn,mlp,nlp->lp', weight, A, B) * dk_norm
    return chi0


if __name__ == "__main__":
    graser_H_mod = load_module("_graser1003_hamiltonian_2026-09-22.py")
    build_H_graser = graser_H_mod.build_H_graser
    mu_graser = float(np.load(os.path.join(SCRIPT_DIR, "_graser_mu_value.npy")))

    def H_graser_kz0(kx, ky):
        return build_H_graser(kx, ky, kz=0.0)

    rng = np.random.default_rng(20260922)
    k_grid = 16  # small for gate speed; not a headline-number run
    T, eta = 0.02, 0.02

    test_points = [((0.0, 0.0), 0.0), ((np.pi, 0.0), 0.0),
                   ((0.7, 1.3), 0.2), ((np.pi, np.pi), 0.5)]
    n_trials = 8
    max_diff_overall = 0.0

    for (q, omega) in test_points:
        baseline = chi0_dynamic(q, omega, H_graser_kz0, mu_graser, T=T, eta=eta,
                                 k_grid=k_grid, rephase=False)
        diffs = []
        for trial in range(n_trials):
            rephased = chi0_dynamic(q, omega, H_graser_kz0, mu_graser, T=T, eta=eta,
                                     k_grid=k_grid, rng=rng, rephase=True)
            d = np.max(np.abs(baseline - rephased))
            diffs.append(d)
        max_d = max(diffs)
        max_diff_overall = max(max_diff_overall, max_d)
        print(f"q={q}, omega={omega}: max|baseline-rephased| over {n_trials} random "
              f"complex-phase trials = {max_d:.3e}  (baseline max|chi0|={np.max(np.abs(baseline)):.3e})")

    print()
    if max_diff_overall < 1e-9:
        print(f"GATE PASS: four-factor formula is gauge-invariant to machine precision "
              f"(max diff {max_diff_overall:.3e} across all {len(test_points)} (q,omega) points x "
              f"{n_trials} random complex rephasings each).")
    else:
        print(f"GATE FAIL: max diff {max_diff_overall:.3e} -- formula is NOT gauge-invariant, "
              f"do not trust downstream Im chi(q,omega) numbers until this is resolved.")
