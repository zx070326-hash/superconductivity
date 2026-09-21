#!/usr/bin/env python3
"""
深度调试：逐步检查 Lindhard 计算的每个环节
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

# 单个 k 点测试
print("=" * 70)
print("Single k-point debug: Lindhard at k=(0,0), q=(0,0)")
print("=" * 70)

kx, ky = 0.0, 0.0
qx, qy = 0.0, 0.0
mu = 0.0
T = 0.01
eta = 0.01

# k 点
H_k = build_5orbital_hamiltonian(kx, ky)
eig_k, U_k = eigh(H_k)
f_k = np.array([fermi_func(e, mu, T) for e in eig_k])

print(f"\nAt k = ({kx}, {ky}):")
print(f"  Eigenvalues: {eig_k}")
print(f"  Fermi occupation: {f_k}")

# k+q 点（这里 q=0）
H_kq = build_5orbital_hamiltonian(kx + qx, ky + qy)
eig_kq, U_kq = eigh(H_kq)
f_kq = np.array([fermi_func(e, mu, T) for e in eig_kq])

print(f"\nAt k+q = ({kx+qx}, {ky+qy}):")
print(f"  Eigenvalues: {eig_kq}")
print(f"  Fermi occupation: {f_kq}")

# 带基 Lindhard（对角元 m=n）
print(f"\nBand-basis Lindhard (diagonal, m=n):")
chi0_band_diag = []
for m in range(5):
    denom = eig_k[m] - eig_kq[m] + 1j*eta
    numer = f_k[m] - f_kq[m]

    # q=0 时应该用 L'Hospital
    if abs(qx) < 1e-10 and abs(qy) < 1e-10:
        # χ₀ ~ -∂f/∂ε ~ f'(ε_m) / T
        # 对于 q=0: χ₀_{mm} = -∂f_m/∂ε_m ≈ f(1-f)/T
        chi0_lh = f_k[m] * (1 - f_k[m]) / T
        chi0_band_diag.append(chi0_lh)
        print(f"  m={m}: f={f_k[m]:.4f}, chi0(L'Hospital)={chi0_lh:.6f}")
    else:
        chi0_mn = numer / denom
        chi0_band_diag.append(chi0_mn.real)
        print(f"  m={m}: (f_k-f_kq)/(e_k-e_kq) = {numer:.6f}/{denom:.6f} = {chi0_mn:.6f}")

# 投影到轨道基（只看对角）
print(f"\nProjection to orbital basis (diagonal elements):")
chi0_orb = np.zeros((5, 5), dtype=complex)

for m in range(5):
    for n in range(5):
        # 简化：只用对角带贡献
        if abs(eig_k[m] - eig_kq[n]) < 1e-10:
            chi0_band_mn = f_k[m] * (1 - f_k[m]) / T
        else:
            chi0_band_mn = (f_k[m] - f_kq[n]) / (eig_k[m] - eig_kq[n] + 1j*eta)

        for i in range(5):
            for j in range(5):
                chi0_orb[i, j] += np.conj(U_k[i, m]) * chi0_band_mn * U_kq[j, n]

print(f"chi0_orbital diagonal:")
for i in range(5):
    print(f"  chi0[{i},{i}] = {chi0_orb[i,i].real:.6f}")

print(f"\nchi0_orbital eigenvalues: {np.linalg.eigvalsh(chi0_orb.real)}")

print("\n" + "=" * 70)
print("EXPECTED at q=0:")
print("  - chi0_band diagonal ~ f(1-f)/T (all positive)")
print("  - chi0_orbital should preserve positivity (unitary transform)")
print("  - Typical value ~ 0.01-1 (not -99 to +97!)")
print("=" * 70)
