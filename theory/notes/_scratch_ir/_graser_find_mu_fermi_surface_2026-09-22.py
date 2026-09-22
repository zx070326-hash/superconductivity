#!/usr/bin/env python3
"""
Graser 5-orbital BaFe2As2: find chemical potential mu for total filling n=6 electrons/Fe
(Fe2+ d6 configuration -> filling fraction 6/10 across the 5 orbitals x 2 spins),
then check qualitative Fermi surface topology (Gamma hole pockets + zone-corner/edge
electron pockets) against Graser Fig. 1-2 before proceeding to Lindhard chi0.
"""
import numpy as np
import importlib.util, os
spec = importlib.util.spec_from_file_location(
    "graser_H", os.path.join(os.path.dirname(__file__), "_graser1003_hamiltonian_2026-09-22.py"))
graser_H = importlib.util.module_from_spec(spec)
spec.loader.exec_module(graser_H)
build_H_graser = graser_H.build_H_graser

N = 120
ks = np.linspace(-np.pi, np.pi, N, endpoint=False)
all_eigs = np.zeros((N, N, 5))
for i, kx in enumerate(ks):
    for j, ky in enumerate(ks):
        H = build_H_graser(kx, ky)
        all_eigs[i, j, :] = np.linalg.eigvalsh(H)

flat = all_eigs.reshape(-1)
flat_sorted = np.sort(flat)

# T=0 filling: want <n> = 6 electrons/Fe out of 10 states (5 orbitals x 2 spins)
# each eigenvalue counted once per k-point represents ONE orbital-band state;
# with spin degeneracy 2, filled fraction of the 5N^2 eigenvalues that must be
# below mu is 6/10 = 0.6 (spin cancels out of the fraction).
target_frac = 0.6
idx = int(round(target_frac * len(flat_sorted))) - 1
mu = flat_sorted[idx]
print(f"N={N}x{N} k-grid, {len(flat_sorted)} eigenvalues total")
print(f"mu (target filling frac=0.6) = {mu:.4f} eV")

# cross-check with actual fraction below this mu
frac_below = np.mean(flat < mu)
print(f"actual fraction below mu: {frac_below:.4f}  (target 0.6000)")

# Fermi surface check: at k=Gamma=(0,0) and various points near it, is there a band
# crossing mu with hole-like (downward) curvature? Near M=(pi,pi) or (pi,0), electron-like?
def bands_at(kx, ky):
    H = build_H_graser(kx, ky)
    return np.sort(np.linalg.eigvalsh(H))

print("\n--- band energies relative to mu at high-symmetry points ---")
for name, (kx, ky) in [("Gamma", (0, 0)), ("X=(pi,0)", (np.pi, 0)),
                        ("M=(pi,pi)", (np.pi, np.pi))]:
    e = bands_at(kx, ky) - mu
    print(f"{name:12s}: {np.round(e, 4)}")

# small-radius scan around Gamma and M to see if bands cross mu (pocket exists)
def pocket_scan(k0, label, radius=0.3, npts=9):
    print(f"\n-- scan around {label} (radius {radius}) --")
    crossings = 0
    for ang in np.linspace(0, 2*np.pi, npts, endpoint=False):
        kx = k0[0] + radius*np.cos(ang)
        ky = k0[1] + radius*np.sin(ang)
        e = bands_at(kx, ky) - mu
        n_below = np.sum(e < 0)
        n_above = np.sum(e > 0)
        crossings += (0 < n_below < 5)  # partial filling = band crosses mu somewhere on this ring
    print(f"  rings with mixed occupied/empty bands: {crossings}/{npts}")

pocket_scan((0, 0), "Gamma")
pocket_scan((np.pi, np.pi), "M=(pi,pi)")
pocket_scan((np.pi, 0), "X=(pi,0)")

np.save(os.path.join(os.path.dirname(__file__), "_graser_mu_value.npy"), mu)
print(f"\nmu saved to _graser_mu_value.npy = {mu:.6f}")
