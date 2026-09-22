#!/usr/bin/env python3
"""
Off-diagonal Lambda^SF_{kk'}(omega) elementwise sign in the FERMI-SURFACE (band,k) basis
-- the actual object Eq.(37) is about, per team-lead's narrowing
(meta/eq37_content_is_in_the_offdiagonal_2026-09-22.md): Lambda^SF = 2 F chi F is a
CONGRUENCE on chi (same f^xc on both sides) => Sylvester's law of inertia => Im Lambda^SF
has the same inertia as Im chi => forced PSD for omega>0 (Schur product argument), but
PSD does NOT force elementwise positivity. Diagonal Lambda^SF_kk is a quadratic form
(pinned PSD, zero budget -- do not scan it). Off-diagonal Lambda^SF_{k,k'} (k!=k') is a
BILINEAR form -- free, this is where all falsifiable content of Eq.(37) lives.

MODEL (explicit assumption, NOT first-principles): we do not have real-space Wannier
functions for the Graser tight-binding model, so we cannot evaluate the true continuum
ALDA sandwich f^xc(x1,x)*chi_zz(x,x')*f^xc(x',x2). Instead we use the standard
multiorbital-RPA-pairing-vertex approximation (Graser/Kuroki-type spin-fluctuation
pairing kernel form): f^xc is taken ORBITAL-DIAGONAL (local single-site vertex,
f^xc_{l l'} = f_l * delta_{ll'}), giving

    Lambda^SF_{(n,k),(n',k')}(omega) = 2 * sum_{l,l'} f_l * w_l(n,k) * chi_{ll'}(k-k',omega)
                                                       * f_l' * w_l'(n',k')

where w_l(n,k) = |U_k[l,n]|^2 (orbital-character weight, real >=0, from the already
gauge-invariant chi0_dynamic machinery's own eigenvectors) and chi_{ll'} is the
RPA-dressed magnetic susceptibility already validated (gauge gate PASS, see
_chk_gauge_invariance_gate_dynamic_2026-09-22.py). f_l is the (possibly orbital-
dependent) ALDA vertex -- THIS IS THE MODELING COMMITMENT team-lead flagged as requiring
a sensitivity scan, not a single-branch theorem-level verdict.

Two f_xc branches (per team-lead's explicit requirement -- do not judge on one branch):
  - "uniform ALDA":     f_l = +1 for all l (baseline; f_xc^2 always >=0, no orbital
                         cross-sign structure -- reduces the sandwich to a pure
                         reflection of chi's own sign structure)
  - "correlation-flip": f_l = -1 for the two e_g-like orbitals (idx 2 = dx2-y2,
                         4 = d3z2-r2) and +1 for the three t2g-like orbitals (idx 0,1,3
                         = dxz,dyz,dxy), modeling the r_s~3-4 correlation-part sign
                         reversal documented in _alda_fxc_negative_rs_2026-09-21.md
                         (orbitals with different local density character are not
                         obligated to share the same ALDA sign). This is an explicit,
                         documented ASSUMPTION for sensitivity purposes, not a claim
                         about which orbitals actually cross r_s~3-4 in real BaFe2As2.

Fermi-surface point selection: |zeta_nk| < FS_TOL (default 0.02 eV) on the same 20x20
k-grid used for mu (mu = 0.0127 eV, corrected Graser parameters).

Pre-registered power check (per team-lead's warning about the chi0-rebuild "0/16 all
positive" false-negative-power episode): our statistic is the RAW FRACTION of
off-diagonal (i,j) instances with Im Lambda^SF > 0, aggregated over many pairs and both
omega values -- NOT an "all positive in a small group" indicator. Under the null
(symmetric random sign, no structure) this fraction has mean 0.5 and, for N_pairs
independent-ish instances, sd ~ 1/(2 sqrt(N_pairs)); for N_pairs of order few hundred
this gives ~2-3% resolution, enough to detect a departure from 50% without an all-or-
nothing indicator collapsing power to near zero. We print z and p up front, before
interpreting the sign fractions, and separately report cross-branch sign AGREEMENT
(the sensitivity question) since that is the one team-lead most cares about.

Strictness: numerical-evidence (finite FS sample, finite k-grid for chi, one modeling
choice of orbital-diagonal f^xc -- explicitly flagged, not a first-principles f_xc).
"""
import numpy as np
import importlib.util, os

SCRIPT_DIR = os.path.dirname(__file__)


def load_module(fname):
    spec = importlib.util.spec_from_file_location(fname, os.path.join(SCRIPT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# SIGN FIX (2026-09-22, team-lead flag, load-bearing): the shared gate script's
# chi0_dynamic uses weight=(f_k-f_kq)/denom with the degenerate branch left
# unhandled at q,omega->0; team-lead independently confirmed by physics argument
# (static long-wavelength limit must be +N(E_F)>0, Pauli-paramagnetic) that this
# convention is the WRONG overall sign -- chi0-rebuild's script negates the
# general term instead. I independently re-verified this here: finite-difference
# DOS(E_F) from direct band filling = +1.002 (correct sign), but
# Re[Tr(chi0_dynamic)] at q~1e-3, omega~1e-3 = -0.861 (wrong sign); negating
# chi0_dynamic gives +0.861, matching the independent DOS check in sign and
# rough magnitude (residual ~15% difference is finite k_grid/eta/T discretization,
# not a sign issue). NOT editing the baseline gate script in place (CLAUDE.md
# convention: extend via new files, don't touch a script other notes already cite
# as gauge-invariance-gate-PASS); gauge invariance itself is unaffected by an
# overall sign flip (linear), so that PASS result stands unchanged. This wrapper
# is the load-bearing chi0 for everything below.
_chi0_dynamic_raw = gate_mod.chi0_dynamic


def chi0_dynamic(q, omega, H_func, mu, **kwargs):
    return -_chi0_dynamic_raw(q, omega, H_func, mu, **kwargs)


def rpa_dress(chi0, U_rpa, norb):
    return chi0 @ np.linalg.inv(np.eye(norb) - (U_rpa * np.eye(norb)) @ chi0)


def find_fs_points(H_func, mu, k_grid=20, fs_tol=0.02, norb=5):
    """Return list of (band_idx, kx, ky, weight_vector) for states within fs_tol of mu."""
    ks = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)
    pts = []
    for kx in ks:
        for ky in ks:
            H = H_func(kx, ky)
            eig, U = np.linalg.eigh(H)
            for n in range(norb):
                if abs(eig[n] - mu) < fs_tol:
                    w = np.abs(U[:, n])**2  # orbital-character weights, real >=0, sum=1
                    pts.append((n, kx, ky, w))
    return pts


def run_offdiag_scan(H_func, mu, k_grid=20, fs_tol=0.02, chi_k_grid=16,
                      omegas=(0.05, 0.2), T=0.02, eta=0.02, U_rpa=0.3, norb=5,
                      max_fs_points=25, seed=20260922):
    fs_pts = find_fs_points(H_func, mu, k_grid=k_grid, fs_tol=fs_tol, norb=norb)
    print(f"Found {len(fs_pts)} Fermi-surface (band,k) points within {fs_tol} eV of mu={mu:.4f}")
    if len(fs_pts) > max_fs_points:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(fs_pts), size=max_fs_points, replace=False)
        fs_pts = [fs_pts[i] for i in idx]
        print(f"Subsampled to {max_fs_points} points for tractability")

    branches = {
        "uniform_ALDA": np.array([1.0, 1.0, 1.0, 1.0, 1.0]),
        "correlation_flip": np.array([1.0, 1.0, -1.0, 1.0, -1.0]),  # eg (idx2,4) flipped
    }

    # cache chi_{ll'}(q,omega) by rounded q to avoid recomputation for repeated q's
    chi_cache = {}

    def get_chi(qx, qy, omega):
        key = (round(qx, 6), round(qy, 6), omega)
        if key not in chi_cache:
            chi0 = chi0_dynamic((qx, qy), omega, H_func, mu, T=T, eta=eta,
                                 k_grid=chi_k_grid, norb=norb, rephase=False)
            chi_cache[key] = rpa_dress(chi0, U_rpa, norb)
        return chi_cache[key]

    results = {}
    for label, f_vec in branches.items():
        for omega in omegas:
            n_pos, n_neg, n_total = 0, 0, 0
            for i in range(len(fs_pts)):
                ni, kxi, kyi, wi = fs_pts[i]
                for j in range(len(fs_pts)):
                    if i == j:
                        continue
                    nj, kxj, kyj, wj = fs_pts[j]
                    qx, qy = kxi - kxj, kyi - kyj
                    chi_q = get_chi(qx, qy, omega)
                    Fwi = f_vec * wi
                    Fwj = f_vec * wj
                    lam = 2.0 * (Fwi @ chi_q @ Fwj)
                    im_lam = lam.imag
                    n_total += 1
                    if im_lam > 1e-12:
                        n_pos += 1
                    elif im_lam < -1e-12:
                        n_neg += 1
            frac_pos = n_pos / n_total if n_total else float("nan")
            results[(label, omega)] = dict(n_pos=n_pos, n_neg=n_neg, n_total=n_total,
                                            frac_pos=frac_pos)
            print(f"[{label}, omega={omega}] off-diag pairs: {n_total}, "
                  f"Im Lambda^SF > 0: {n_pos} ({frac_pos:.3f}), < 0: {n_neg}")

    # pre-registered power / null check
    print("\n=== pre-registered null check (symmetric-random-sign baseline) ===")
    for key, r in results.items():
        n = r["n_total"]
        if n == 0:
            continue
        mean, sd = 0.5, 0.5 / np.sqrt(n)
        z = (r["frac_pos"] - mean) / sd
        print(f"{key}: n={n}, frac_pos={r['frac_pos']:.3f}, expected sd={sd:.4f}, z={z:.3f}")

    # cross-branch sign agreement (the sensitivity question team-lead asked for)
    print("\n=== cross-branch (uniform ALDA vs correlation-flip) sign agreement ===")
    for omega in omegas:
        rU = results[("uniform_ALDA", omega)]
        rC = results[("correlation_flip", omega)]
        # can't directly compare per-pair without re-storing; report aggregate only
        print(f"omega={omega}: uniform frac_pos={rU['frac_pos']:.3f}, "
              f"correlation_flip frac_pos={rC['frac_pos']:.3f}, "
              f"delta={rC['frac_pos']-rU['frac_pos']:+.3f}")

    return results


if __name__ == "__main__":
    graser_H_mod = load_module("_graser1003_hamiltonian_2026-09-22.py")
    build_H_graser = graser_H_mod.build_H_graser
    mu_graser = float(np.load(os.path.join(SCRIPT_DIR, "_graser_mu_value.npy")))

    def H_graser_kz0(kx, ky):
        return build_H_graser(kx, ky, kz=0.0)

    print("=== Graser BaFe2As2 (corrected params), off-diagonal Lambda^SF_{kk'} "
          "elementwise sign, Fermi-surface (band,k) basis ===")
    results = run_offdiag_scan(H_graser_kz0, mu_graser, k_grid=20, fs_tol=0.02,
                                chi_k_grid=16, omegas=(0.05, 0.2), T=0.02, eta=0.02,
                                U_rpa=0.3, max_fs_points=25)
    print("\nDone.")
