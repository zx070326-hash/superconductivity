#!/usr/bin/env python3
"""
Dynamic Im chi0_{ll'}(q, omega>0) for the REAL Graser et al. 1003.0133 BaFe2As2
5-orbital H(k), with:
  (a) the four-factor cross-term projection (gauge invariant -- see
      chi0_two_factor_form_is_gauge_dependent_2026-09-22.md, and
      two_band_chi0_derivation_2026-09-22.md for the Wick-contraction derivation),
  (b) the SIGN FIX team-lead identified in _chk_q0_vs_qneq0_signature_2026-09-22.py /
      the second teammate-message of 2026-09-22: the NON-degenerate branch
      c_mn = (f_m - f_n)/(e_m - e_n + i*eta) must be NEGATED (chi0 = -sum(...)) to match
      the physical convention already used (correctly) in the degenerate branch
      c_mm = +f(1-f)/T.  Do NOT touch the degenerate branch.

Per team-lead's explicit budget call: static Re chi0(q,0) is now proven to be a
zero-degrees-of-freedom theorem (NSD everywhere, any material) -- rebuilding/
rescanning it has infinite budget and is explicitly out of scope here.  The only
channel with remaining degrees of freedom is Im chi0(q, omega>0); this script
targets exactly that.

MANDATORY GATE (per team-lead, "this door has not been run by anyone yet"):
  gauge_invariance_gate() must print PASS before any (q,omega) production numbers
  are reported.  It resigns eigenvectors with random per-(k,band) +/-1 (Graser H(k)
  is real symmetric at all sampled k -- eigh returns real eigenvectors, so the
  residual gauge is exactly a sign, matching the "real branch" case team-lead
  exhaustively enumerated for the static case) and checks the full k-summed chi0
  output is unchanged to machine precision.

Rigor labels (see report_2026-09-22 note for the full breakdown):
  - four-factor formula gauge invariance: analytic (termwise phase cancellation)
  - this script's gate: numerical-evidence (random resigning, not exhaustive --
    exhaustive is infeasible over a 20x20 k-grid; termwise invariance is what the
    analytic argument guarantees, so per-(k,band) random resigning is a faithful
    stress test of that argument, not an independent unproven claim)
  - elementwise Im chi0 sign statistics, RPA sign change, amplitude ratios:
    numerical-evidence (float64, finite k-grid/eta)
"""

import itertools
import numpy as np
from numpy.linalg import eigh, inv
import importlib.util, os

_here = os.path.dirname(__file__)
spec = importlib.util.spec_from_file_location(
    "graser_H", os.path.join(_here, "_graser1003_hamiltonian_2026-09-22.py"))
graser_H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(graser_H)
build_H_graser = graser_H.build_H_graser

NORB = 5
MU = -0.0656  # _graser_find_mu_fermi_surface_2026-09-22.py, filling=0.6000 exact
T = 0.01


def fermi(e, mu=MU, T=T):
    return 1.0 / (1.0 + np.exp(np.clip((e - mu) / T, -500, 500)))


def chi0_dynamic(q, omega, eta=0.02, k_grid=20, sign_k=None, sign_kq=None):
    """Four-factor cross-term chi0_{ll'}(q, omega+i*eta), sign-fixed non-degenerate
    branch.  sign_k/sign_kq: optional dict {(ik,jk): array of +/-1 length NORB} for
    the gauge-invariance gate; None means "use whatever LAPACK returns" (no resign).
    """
    qx, qy = q
    chi = np.zeros((NORB, NORB), dtype=complex)
    kx_arr = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)
    dk_norm = 1.0 / k_grid ** 2

    for ik, kx in enumerate(kx_arr):
        for jk, ky in enumerate(ky_arr):
            e_k, u_k = eigh(build_H_graser(kx, ky))
            e_q, u_q = eigh(build_H_graser(kx + qx, ky + qy))

            if sign_k is not None:
                u_k = u_k * sign_k[(ik, jk)][None, :]
            if sign_kq is not None:
                u_q = u_q * sign_kq[(ik, jk)][None, :]

            f_k, f_q = fermi(e_k), fermi(e_q)

            for m in range(NORB):
                factor1 = np.conj(u_k[:, m])[:, None] * u_k[:, m][None, :]
                for n in range(NORB):
                    d = e_k[m] - e_q[n]
                    if abs(d) < 1e-10:
                        # degenerate branch -- ALREADY correct sign, do not touch
                        c = fermi(e_k[m]) * (1 - fermi(e_k[m])) / T
                    else:
                        # non-degenerate branch -- SIGN FIX: negate.
                        # standard retarded convention: chi0(q,w+i0+) = sum (f_m-f_n)/(w+i eta+e_m-e_n),
                        # this is the ANALYTIC CONTINUATION of the w=0 static formula used upstream
                        # (_chk_q0_vs_qneq0_signature_2026-09-22.py, d = e_k[m]-e_q[n], denom = d+i*eta
                        # at w=0), so the frequency enters additively: denom = w + d + i*eta.
                        c = -(f_k[m] - f_q[n]) / (omega + d + 1j * eta)
                    factor2 = u_q[:, n][:, None] * np.conj(u_q[:, n])[None, :]
                    chi += c * factor1 * factor2
    return chi * dk_norm


def gauge_invariance_gate(n_trials=60, k_grid=10, q=(np.pi, 0.0), omega=0.2, eta=0.02):
    """Random per-(k,band) +/-1 resign (real branch, matches Graser's real H(k)).
    PASS iff max|Delta chi0| over all trials is at machine precision."""
    print("=" * 70)
    print("GAUGE INVARIANCE GATE (mandatory, per team-lead 2026-09-22)")
    print("=" * 70)
    rng = np.random.default_rng(20260922)
    base = chi0_dynamic(q, omega, eta=eta, k_grid=k_grid)

    worst = 0.0
    for t in range(n_trials):
        sign_k = {(i, j): rng.choice([1.0, -1.0], size=NORB)
                   for i in range(k_grid) for j in range(k_grid)}
        sign_kq = {(i, j): rng.choice([1.0, -1.0], size=NORB)
                   for i in range(k_grid) for j in range(k_grid)}
        resigned = chi0_dynamic(q, omega, eta=eta, k_grid=k_grid,
                                 sign_k=sign_k, sign_kq=sign_kq)
        worst = max(worst, np.abs(base - resigned).max())

    print(f"  {n_trials} random per-(k,band) sign-gauge trials, k_grid={k_grid}, q={q}, omega={omega}")
    print(f"  max|Delta chi0| over all trials = {worst:.3e}")
    passed = worst < 1e-9
    print(f"  GATE: {'PASS' if passed else 'FAIL'}")
    print()
    return passed


def elementwise_sign_report(q, label, omegas, k_grid=20, eta=0.02, U_spin=0.5):
    print("=" * 70)
    print(f"q = {label} = {q}")
    print("=" * 70)
    results = []
    for w in omegas:
        chi0 = chi0_dynamic(q, w, eta=eta, k_grid=k_grid)
        im0 = chi0.imag
        # elementwise off-diagonal sign stats (upper triangle, l<l', 10 pairs)
        offdiag = [(l, lp, im0[l, lp]) for l in range(NORB) for lp in range(l+1, NORB)]
        n_neg_off = sum(1 for _, _, v in offdiag if v < -1e-8)
        n_pos_off = sum(1 for _, _, v in offdiag if v > 1e-8)
        diag_sign = np.sign(np.round(np.diag(im0), 8))

        I = np.eye(NORB)
        try:
            chi_s = chi0 @ inv(I - U_spin * chi0)
            im_s = chi_s.imag
            offdiag_s = [(l, lp, im_s[l, lp]) for l in range(NORB) for lp in range(l+1, NORB)]
            n_neg_off_s = sum(1 for _, _, v in offdiag_s if v < -1e-8)
            n_pos_off_s = sum(1 for _, _, v in offdiag_s if v > 1e-8)
            flips = sum(1 for (l, lp, v0), (_, _, vs) in zip(offdiag, offdiag_s)
                        if np.sign(np.round(v0, 8)) != 0 and np.sign(np.round(vs, 8)) != 0
                        and np.sign(v0) != np.sign(vs))
        except np.linalg.LinAlgError:
            n_neg_off_s = n_pos_off_s = flips = None

        fro = np.linalg.norm(im0, 'fro')
        print(f"  omega={w:.3f}: Im chi0 off-diag(10) neg={n_neg_off} pos={n_pos_off} "
              f"zero={10-n_neg_off-n_pos_off} | diag signs={diag_sign.astype(int)} "
              f"| ||Im chi0||_F={fro:.4f}")
        print(f"             RPA Im chi_s off-diag neg={n_neg_off_s} pos={n_pos_off_s} "
              f"| sign flips vs chi0: {flips}/10")
        results.append(dict(omega=w, n_neg_off=n_neg_off, n_pos_off=n_pos_off,
                             fro=fro, n_neg_off_s=n_neg_off_s, n_pos_off_s=n_pos_off_s,
                             flips=flips))
    print()
    return results


def main():
    gate_ok = gauge_invariance_gate()
    if not gate_ok:
        print("GATE FAILED -- refusing to report any (q,omega) numbers below.")
        return

    omegas = [0.05, 0.2, 0.5, 1.0]
    k_grid = 20

    all_results = {}
    all_results["Gamma (0,0)"] = elementwise_sign_report((0.0, 0.0), "Gamma (0,0)", omegas, k_grid=k_grid)
    all_results["Q_AFM (pi,0)"] = elementwise_sign_report((np.pi, 0.0), "Q_AFM (pi,0)", omegas, k_grid=k_grid)
    all_results["M (pi,pi)"] = elementwise_sign_report((np.pi, np.pi), "M (pi,pi)", omegas, k_grid=k_grid)
    all_results["generic (0.5pi,1.1pi)"] = elementwise_sign_report(
        (0.5*np.pi, 1.1*np.pi), "generic (0.5pi,1.1pi)", omegas, k_grid=k_grid)

    print("=" * 70)
    print("AMPLITUDE RATIO: ||Im chi0(q,omega)||_F relative to Q_AFM")
    print("=" * 70)
    for w_idx, w in enumerate(omegas):
        qafm_fro = all_results["Q_AFM (pi,0)"][w_idx]["fro"]
        print(f"  omega={w:.3f}:")
        for label, res in all_results.items():
            ratio = res[w_idx]["fro"] / qafm_fro if qafm_fro > 0 else float('nan')
            print(f"    {label:24s} ||Im chi0||_F={res[w_idx]['fro']:.4f}  ratio to Q_AFM={ratio:.4f}")


if __name__ == "__main__":
    main()
