#!/usr/bin/env python3
"""
检查 k 空间积分的归一化

问题：chi0 ~ 100 量级，但应该 ~ DOS ~ 1
可能原因：dk 归一化不对
"""

import numpy as np
from numpy.linalg import eigh

def build_5orbital_hamiltonian(kx, ky):
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

# 计算总态密度（作为参考）
print("=" * 70)
print("Check k-space integral normalization")
print("=" * 70)

mu = 0.0
T = 0.01

k_grid = 50
kx_arr = np.linspace(-np.pi, np.pi, k_grid)
ky_arr = np.linspace(-np.pi, np.pi, k_grid)

# 方法1：直接算 dk = (2π/N)^2
dk_v1 = (2*np.pi / k_grid)**2
print(f"\nMethod 1: dk = (2π/N)^2 = {dk_v1:.6f}")
print(f"  Total k-space volume: N * dk = {k_grid**2 * dk_v1:.6f}")
print(f"  Expected: (2π)^2 = {(2*np.pi)**2:.6f}")

# 方法2：用 np.gradient
dk_v2 = (kx_arr[1] - kx_arr[0]) * (ky_arr[1] - ky_arr[0])
print(f"\nMethod 2: dk = Δkx * Δky = {dk_v2:.6f}")
print(f"  Total k-space volume: N * dk = {k_grid**2 * dk_v2:.6f}")

# 计算总态密度（应该 ~ 每个轨道 O(1)）
dos_total = 0.0
dos_by_band = np.zeros(5)

for kx in kx_arr:
    for ky in ky_arr:
        H = build_5orbital_hamiltonian(kx, ky)
        eig, _ = eigh(H)

        for n in range(5):
            f = fermi_func(eig[n], mu, T)
            dos_contrib = f * (1 - f) / T  # f'(E) ~ f(1-f)/T
            dos_by_band[n] += dos_contrib * dk_v1

dos_total = dos_by_band.sum()

print(f"\nDOS at Fermi level (sum rule check):")
print(f"  By band: {dos_by_band}")
print(f"  Total: {dos_total:.4f}")
print(f"  Per band average: {dos_total/5:.4f}")
print(f"\nExpected: O(1) per band for tight-binding model")

# 直接计算 chi0(q=0) 对角元
print(f"\n" + "="*70)
print("chi0(q=0) diagonal - direct calculation:")

chi0_diag_sum = 0.0

for kx in kx_arr:
    for ky in ky_arr:
        H = build_5orbital_hamiltonian(kx, ky)
        eig, _ = eigh(H)

        for n in range(5):
            f = fermi_func(eig[n], mu, T)
            chi0_nn = f * (1 - f) / T
            chi0_diag_sum += chi0_nn * dk_v1

print(f"  Sum of chi0_nn over all bands: {chi0_diag_sum:.4f}")
print(f"  Should equal total DOS: {dos_total:.4f}")
print(f"  Ratio: {chi0_diag_sum / dos_total:.6f} (should be 1.0)")

print("\n" + "="*70)
print("DIAGNOSIS:")
if abs(chi0_diag_sum - dos_total) < 1e-6:
    print("  ✓ Normalization correct")
else:
    print("  ✗ Normalization problem!")

if dos_total > 100:
    print(f"  ✗ DOS = {dos_total:.1f} >> 1, likely wrong dk or BZ size")
elif dos_total < 0.01:
    print(f"  ✗ DOS = {dos_total:.6f} << 1, likely missing (2π)^2 factor")
else:
    print(f"  ✓ DOS ~ {dos_total:.2f}, reasonable magnitude")
print("="*70)
