#!/usr/bin/env python3
"""
FeSe 5-orbital RPA - FINAL CORRECTED VERSION

所有已修复的问题：
1. ✓ Lindhard 公式无多余负号
2. ✓ 轨道投影公式正确
3. ✓ k 空间归一化：除以 (2π)²
4. ✓ 带基计算后投影
"""

import numpy as np
from numpy.linalg import eigh, inv

def build_5orbital_hamiltonian(kx, ky):
    """简化 5-轨道紧束缚模型"""
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
    """费米分布"""
    x = (E - mu) / T
    return 1.0 / (1.0 + np.exp(np.clip(x, -50, 50)))

def lindhard_chi0_5orbital_FINAL(q, mu, T=0.01, eta=0.01, k_grid=30):
    """
    最终修正版 Lindhard χ₀

    关键修正：
    - dk 归一化：dk / (2π)²
    - 公式：χ₀ = Σ_k [f_k - f_{k+q}] / [ε_k - ε_{k+q} + iη]
    - 投影：U†(k) χ₀^band U(k+q)
    """
    qx, qy = q
    chi0_orbital = np.zeros((5, 5), dtype=complex)

    kx_arr = np.linspace(-np.pi, np.pi, k_grid)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid)
    dk = (2*np.pi / k_grid)**2

    # 关键：除以布里渊区面积 (2π)²
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

            # 带基 Lindhard
            for m in range(5):
                for n in range(5):
                    denom = eig_k[m] - eig_kq[n] + 1j*eta

                    # q=0 且 m=n：用 L'Hospital
                    if abs(qx) < 1e-10 and abs(qy) < 1e-10 and m == n:
                        chi0_band_mn = f_k[m] * (1 - f_k[m]) / T
                    else:
                        chi0_band_mn = (f_k[m] - f_kq[n]) / denom

                    # 投影到轨道基
                    for i in range(5):
                        for j in range(5):
                            chi0_orbital[i, j] += (
                                np.conj(U_k[i, m]) * chi0_band_mn * U_kq[j, n]
                                * dk_normalized
                            )

    return chi0_orbital

def rpa_spin_susceptibility(chi0, U_spin):
    """RPA: χ^s = χ₀ (1 - U χ₀)^{-1}"""
    I = np.eye(chi0.shape[0])
    try:
        chi_s = chi0 @ inv(I - U_spin * chi0)
        return chi_s
    except np.linalg.LinAlgError:
        print("WARNING: RPA inversion failed (Stoner instability?)")
        return None

# ==============================================================================
# 测试
# ==============================================================================

print("=" * 70)
print("FeSe 5-orbital RPA - FINAL VERSION")
print("=" * 70)

mu = 0.0
U_spin = 0.5
T = 0.01

# Test 1: q = 0
print("\n[TEST 1] q = (0, 0)")
print("-" * 70)

chi0_q0 = lindhard_chi0_5orbital_FINAL((0.0, 0.0), mu, T=T, k_grid=30)

print("chi0 diagonal elements:")
for i in range(5):
    print(f"  chi0[{i},{i}] = {chi0_q0[i,i].real:.6f}")

eigs_chi0 = np.linalg.eigvalsh(chi0_q0.real)
print(f"\nchi0 eigenvalues: {eigs_chi0}")
print(f"All positive? {np.all(eigs_chi0 > -1e-6)}")

if np.all(eigs_chi0 > -1e-6):
    print("\n✓ chi0 is positive semi-definite")
else:
    print("\n✗ chi0 has negative eigenvalues")

# RPA
chi_s_q0 = rpa_spin_susceptibility(chi0_q0, U_spin)

if chi_s_q0 is not None:
    eigs_s = np.linalg.eigvalsh(chi_s_q0.real)
    print(f"\nchi_s eigenvalues: {eigs_s}")
    print(f"All positive? {np.all(eigs_s > -1e-6)}")

    print(f"\nOff-diagonal elements:")
    for i in range(4):
        print(f"  chi_s[{i},{i+1}] = {chi_s_q0[i,i+1].real:.6f} "
              f"({'POS' if chi_s_q0[i,i+1].real > 0 else 'NEG'})")

# Test 2: q = (π, 0) - AFM wavevector
print("\n" + "=" * 70)
print("[TEST 2] q = (π, 0) - AFM wavevector")
print("-" * 70)

chi0_qm = lindhard_chi0_5orbital_FINAL((np.pi, 0.0), mu, T=T, k_grid=30)

print("chi0 diagonal elements:")
for i in range(5):
    print(f"  chi0[{i},{i}] = {chi0_qm[i,i].real:.6f}")

eigs_chi0_m = np.linalg.eigvalsh(chi0_qm.real)
print(f"\nchi0 eigenvalues: {eigs_chi0_m}")
print(f"All positive? {np.all(eigs_chi0_m > -1e-6)}")

chi_s_qm = rpa_spin_susceptibility(chi0_qm, U_spin)

if chi_s_qm is not None:
    eigs_s_m = np.linalg.eigvalsh(chi_s_qm.real)
    print(f"\nchi_s eigenvalues: {eigs_s_m}")
    print(f"All positive? {np.all(eigs_s_m > -1e-6)}")

    print(f"\nOff-diagonal elements:")
    for i in range(4):
        print(f"  chi_s[{i},{i+1}] = {chi_s_qm[i,i+1].real:.6f} "
              f"({'POS' if chi_s_qm[i,i+1].real > 0 else 'NEG'})")

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("Expected behavior:")
print("  1. chi0 should be PSD (all eigenvalues >= 0)")
print("  2. RPA should enhance but preserve PSD")
print("  3. Off-diagonal signs: system-specific (depends on Fermi surface)")
print("=" * 70)
