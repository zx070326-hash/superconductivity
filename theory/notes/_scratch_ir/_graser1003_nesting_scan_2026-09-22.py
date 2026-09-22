#!/usr/bin/env python3
"""
Fine-grid Fermi-surface nesting function xi(q) for the Graser et al. arXiv:1003.0133
5-orbital BaFe2As2 tight-binding model (kz=0 cut), reusing the corrected H(k) from
_graser1003_hamiltonian_2026-09-22.py (essenberger-chi-reconstruction's build, same
page-image-verified Table I/II values as
theory/notes/literature/graser1003.0133_5orbital_TB_parameters_2026-09-22.md).

Team C priority 1: replace the coarse 0.15pi placeholder-model nesting scan with a
real-parameter, fine (2pi/N, N=128 -> ~0.05pi step) scan, and locate Q_AFM precisely.

MATERIAL CAVEAT (must stay attached to every downstream claim):
  This is BaFe2As2 (122-family), NOT FeSe (11-family) -- Essenberger's actual material.
  Only the FERMI-SURFACE TOPOLOGY class (Gamma-centered hole pockets + M-centered
  electron pockets, same space group / orbital content) is expected to transfer;
  QUANTITATIVE nesting magnitudes, exact peak sharpness, and pocket sizes are
  BaFe2As2-specific and do NOT transfer to FeSe without FeSe-specific DFT input.

Method: xi(q) = sum_k A(k) A(k+q), A(k) = sum_n Lorentzian(eps_n(k) - mu, eta)
        (Lorentzian-broadened joint spectral weight -- standard nesting-function proxy).
Computed via 2D circular autocorrelation (FFT) over the periodic BZ grid -- O(N^2 log N)
instead of the O(N^4) double loop, so a genuinely fine grid is affordable.
"""
import importlib.util
import os
import numpy as np

_here = os.path.dirname(__file__)
_spec = importlib.util.spec_from_file_location(
    "graser1003_hamiltonian", os.path.join(_here, "_graser1003_hamiltonian_2026-09-22.py"))
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
build_H_graser = _mod.build_H_graser

N = 128          # grid -> step 2*pi/128 = 0.0156*2pi = 0.0313pi in each direction (finer than the 0.05pi target)
ETA = 0.03        # Lorentzian broadening (eV), ~ typical experimental/thermal scale
NBAND = 5
FILLING = 0.6 * 2 * NBAND / 10.0  # placeholder, mu found by bisection below on 6 e/Fe = 6/10 of 2*5 states

def find_mu(bands, target_filling_per_site=6.0):
    """bisection on mu s.t. mean occupation (2 spins * 5 orbitals) = target electrons/Fe."""
    flat = bands.reshape(-1, NBAND)
    lo, hi = flat.min() - 1, flat.max() + 1
    for _ in range(60):
        mu = 0.5 * (lo + hi)
        n = 2 * np.mean(flat < mu) * NBAND  # spin-degenerate count of filled (band,k) states
        if n < target_filling_per_site:
            lo = mu
        else:
            hi = mu
    return 0.5 * (lo + hi)

def main():
    kx = np.linspace(-np.pi, np.pi, N, endpoint=False)
    ky = np.linspace(-np.pi, np.pi, N, endpoint=False)
    bands = np.empty((N, N, NBAND))
    for i, kxi in enumerate(kx):
        for j, kyj in enumerate(ky):
            H = build_H_graser(kxi, kyj, 0.0)
            bands[i, j, :] = np.linalg.eigvalsh(H)

    mu = find_mu(bands, target_filling_per_site=6.0)
    print(f"mu (6 e/Fe target) = {mu:.4f} eV  [essenberger's independent solve: -0.0656 eV, cross-check]")

    # A(k) = sum_n Lorentzian(eps_n(k)-mu, eta)
    A = np.sum(ETA / np.pi / ((bands - mu) ** 2 + ETA ** 2), axis=2)

    # circular autocorrelation via FFT: xi(q) = sum_k A(k) A(k+q)
    F = np.fft.fft2(A)
    xi = np.fft.ifft2(F * np.conj(F)).real
    xi = np.fft.fftshift(xi)  # now q=0 at center index (N//2, N//2)

    qx = np.fft.fftshift(np.fft.fftfreq(N, d=1.0/N)) * (2*np.pi/N)  # in [-pi,pi)
    qy = qx.copy()

    # special points (unfolded 1-Fe BZ convention: Gamma=(0,0), M=(pi,pi), X=(pi,0))
    special = {
        "Gamma (0,0)": (0.0, 0.0),
        "X (pi,0)": (np.pi, 0.0),
        "Y (0,pi)": (0.0, np.pi),
        "M (pi,pi)": (np.pi, np.pi),
    }

    def xi_at(qx0, qy0):
        i = np.argmin(np.abs(qx - qx0))
        j = np.argmin(np.abs(qy - qy0))
        return xi[i, j], qx[i], qy[j]

    print(f"grid step = {2*np.pi/N:.4f} rad = {2/N:.4f} pi  (target was 0.05 pi -> {'MET' if 2*np.pi/N <= 0.05*np.pi + 1e-9 else 'NOT MET'})")
    print()
    results = {}
    for name, (qx0, qy0) in special.items():
        val, qxf, qyf = xi_at(qx0, qy0)
        results[name] = val
        print(f"xi({name}) = {val:10.4f}  (grid pt qx={qxf:+.4f}, qy={qyf:+.4f})")

    # global argmax over full grid (exclude q=0 trivial peak, which is just total DOS^2-like)
    xi_noself = xi.copy()
    # mask a small disk around q=0
    ci, cj = N//2, N//2
    rmask = 3
    for di in range(-rmask, rmask+1):
        for dj in range(-rmask, rmask+1):
            xi_noself[(ci+di) % N, (cj+dj) % N] = -np.inf
    imax, jmax = np.unravel_index(np.argmax(xi_noself), xi_noself.shape)
    print()
    print(f"global argmax (excluding q=0 self-peak): q=({qx[imax]:+.4f},{qy[jmax]:+.4f}) "
          f"= ({qx[imax]/np.pi:+.4f},{qy[jmax]/np.pi:+.4f}) pi, xi={xi[imax,jmax]:.4f}")

    avg = np.mean(xi)
    print(f"BZ average xi = {avg:.4f}")
    for name in results:
        print(f"  enhancement at {name}: {results[name]/avg:.3f}x average")

    np.savez(os.path.join(_here, "_graser1003_nesting_scan_2026-09-22.npz"),
             qx=qx, qy=qy, xi=xi, mu=mu, N=N, eta=ETA)

if __name__ == "__main__":
    main()
