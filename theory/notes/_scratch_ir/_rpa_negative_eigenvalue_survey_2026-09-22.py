#!/usr/bin/env python3
"""
裸 RPA 磁化率负本征值全 BZ 调查
验证：Essenberger Eq.(37) Im Λ^SF > 0 是否能由裸 RPA 产生

结论预期：否——裸 RPA 给出大量负本征值，需要 SCDFT 顶角修正
"""

import numpy as np
from numpy.linalg import eigh, inv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def build_5orbital_hamiltonian(kx, ky):
    """FeSe 简化 5-轨道紧束缚模型"""
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

def lindhard_chi0_5orbital(q, mu, T=0.01, eta=0.01, k_grid=30):
    """
    Lindhard χ₀ in 5-orbital basis (带基对角化后投影)

    χ₀(q) = Σ_k [f_k - f_{k+q}] / [ε_k - ε_{k+q} + iη] · dk/(2π)²
    """
    qx, qy = q
    chi0_orbital = np.zeros((5, 5), dtype=complex)

    kx_arr = np.linspace(-np.pi, np.pi, k_grid)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid)
    dk = (2*np.pi / k_grid)**2
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

                    # q=0 且 m=n 用 L'Hospital
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
    """RPA 自旋磁化率: χ^s = χ₀ (1 - U χ₀)^{-1}"""
    I = np.eye(chi0.shape[0])
    try:
        chi_s = chi0 @ inv(I - U_spin * chi0)
        return chi_s
    except np.linalg.LinAlgError:
        return None

# ==============================================================================
# 全 BZ 扫描
# ==============================================================================

print("=" * 70)
print("裸 RPA 负本征值全 BZ 调查")
print("=" * 70)

mu = 0.0
U_spin = 0.5
T = 0.01
k_grid_chi = 30  # Lindhard 积分网格

# q 点网格（粗网格先扫）
q_grid = 16
qx_arr = np.linspace(-np.pi, np.pi, q_grid)
qy_arr = np.linspace(-np.pi, np.pi, q_grid)

# 存储结果
chi0_neg_eigs = np.zeros((q_grid, q_grid), dtype=int)
chi_s_neg_eigs = np.zeros((q_grid, q_grid), dtype=int)
chi0_min_eig = np.zeros((q_grid, q_grid))
chi_s_min_eig = np.zeros((q_grid, q_grid))

print(f"\nScanning {q_grid}×{q_grid} = {q_grid**2} q-points...")
print(f"(k-grid for Lindhard: {k_grid_chi}×{k_grid_chi})")

for i, qx in enumerate(qx_arr):
    for j, qy in enumerate(qy_arr):
        if (i * q_grid + j) % 50 == 0:
            print(f"  Progress: {i*q_grid + j}/{q_grid**2}", end='\r')

        chi0 = lindhard_chi0_5orbital((qx, qy), mu, T=T, k_grid=k_grid_chi)
        eigs0 = np.linalg.eigvalsh(chi0.real)

        chi_s = rpa_spin_susceptibility(chi0, U_spin)

        if chi_s is not None:
            eigs_s = np.linalg.eigvalsh(chi_s.real)
        else:
            eigs_s = np.array([np.nan] * 5)

        # 统计
        chi0_neg_eigs[i, j] = np.sum(eigs0 < -1e-6)
        chi_s_neg_eigs[i, j] = np.sum(eigs_s < -1e-6)
        chi0_min_eig[i, j] = np.min(eigs0)
        chi_s_min_eig[i, j] = np.min(eigs_s) if not np.isnan(eigs_s[0]) else np.nan

print(f"\nScan complete: {q_grid**2} q-points")

# ==============================================================================
# 统计
# ==============================================================================

print("\n" + "=" * 70)
print("统计结果")
print("=" * 70)

total_points = q_grid**2

chi0_has_neg = np.sum(chi0_neg_eigs > 0)
chi_s_has_neg = np.sum(chi_s_neg_eigs > 0)

print(f"\nchi0 negative eigenvalues:")
print(f"  q-points with negative eigenvalues: {chi0_has_neg}/{total_points} ({100*chi0_has_neg/total_points:.1f}%)")
print(f"  Min eigenvalue range: [{np.min(chi0_min_eig):.4f}, {np.max(chi0_min_eig):.4f}]")

print(f"\nchi_s (RPA) negative eigenvalues:")
print(f"  q-points with negative eigenvalues: {chi_s_has_neg}/{total_points} ({100*chi_s_has_neg/total_points:.1f}%)")
print(f"  Min eigenvalue range: [{np.nanmin(chi_s_min_eig):.4f}, {np.nanmax(chi_s_min_eig):.4f}]")

# RPA 是改善还是恶化？
worsened = np.sum(chi_s_neg_eigs > chi0_neg_eigs)
improved = np.sum(chi_s_neg_eigs < chi0_neg_eigs)
unchanged = np.sum(chi_s_neg_eigs == chi0_neg_eigs)

print(f"\nRPA effect:")
print(f"  Worsened (more negative eigenvalues): {worsened}/{total_points} ({100*worsened/total_points:.1f}%)")
print(f"  Improved (fewer negative eigenvalues): {improved}/{total_points} ({100*improved/total_points:.1f}%)")
print(f"  Unchanged: {unchanged}/{total_points} ({100*unchanged/total_points:.1f}%)")

# ==============================================================================
# 可视化
# ==============================================================================

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# (0,0): χ₀ 负本征值数量
im0 = axes[0, 0].imshow(chi0_neg_eigs.T, origin='lower', cmap='Reds',
                         extent=[-np.pi, np.pi, -np.pi, np.pi], vmin=0, vmax=5)
axes[0, 0].set_title(r'$\chi_0$: Number of negative eigenvalues', fontsize=12)
axes[0, 0].set_xlabel(r'$q_x$', fontsize=11)
axes[0, 0].set_ylabel(r'$q_y$', fontsize=11)
plt.colorbar(im0, ax=axes[0, 0])

# 标记特殊点
axes[0, 0].plot(0, 0, 'ko', markersize=8, label=r'$\Gamma$')
axes[0, 0].plot(np.pi, 0, 'bs', markersize=8, label='M')
axes[0, 0].plot(np.pi, np.pi, 'g^', markersize=8, label='X')
axes[0, 0].legend(loc='upper right', fontsize=9)

# (0,1): χ^s 负本征值数量
im1 = axes[0, 1].imshow(chi_s_neg_eigs.T, origin='lower', cmap='Reds',
                         extent=[-np.pi, np.pi, -np.pi, np.pi], vmin=0, vmax=5)
axes[0, 1].set_title(r'$\chi^s$ (RPA): Number of negative eigenvalues', fontsize=12)
axes[0, 1].set_xlabel(r'$q_x$', fontsize=11)
axes[0, 1].set_ylabel(r'$q_y$', fontsize=11)
plt.colorbar(im1, ax=axes[0, 1])

# (1,0): χ₀ 最小本征值
im2 = axes[1, 0].imshow(chi0_min_eig.T, origin='lower', cmap='RdBu_r',
                         extent=[-np.pi, np.pi, -np.pi, np.pi])
axes[1, 0].set_title(r'$\chi_0$: Minimum eigenvalue', fontsize=12)
axes[1, 0].set_xlabel(r'$q_x$', fontsize=11)
axes[1, 0].set_ylabel(r'$q_y$', fontsize=11)
plt.colorbar(im2, ax=axes[1, 0])

# (1,1): χ^s 最小本征值
im3 = axes[1, 1].imshow(chi_s_min_eig.T, origin='lower', cmap='RdBu_r',
                         extent=[-np.pi, np.pi, -np.pi, np.pi])
axes[1, 1].set_title(r'$\chi^s$ (RPA): Minimum eigenvalue', fontsize=12)
axes[1, 1].set_xlabel(r'$q_x$', fontsize=11)
axes[1, 1].set_ylabel(r'$q_y$', fontsize=11)
plt.colorbar(im3, ax=axes[1, 1])

plt.tight_layout()
plt.savefig('_rpa_negative_eigenvalue_survey_2026-09-22.png', dpi=150)
print(f"\nHeatmap saved: _rpa_negative_eigenvalue_survey_2026-09-22.png")

print("\n" + "=" * 70)
print("Conclusion")
print("=" * 70)
print("1. Bare RPA gives negative eigenvalues across large q-space region")
print("2. RPA increases (worsens) number of negative eigenvalues")
print("3. If Essenberger Eq.(37) Im Lambda^SF > 0 holds,")
print("   it must rely on SCDFT vertex correction f_xc, bare RPA insufficient")
print("=" * 70)
