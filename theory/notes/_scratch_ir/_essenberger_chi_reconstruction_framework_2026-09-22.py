#!/usr/bin/env python3
"""
Essenberger FeSe χ^s_{ij} 符号结构复现框架
目标：验证 Eq.(37) 的 Im[χ^s_{kk'}] > 0 声称

Phase 0: 最简 2×2 Lindhard + RPA 框架
- 各向同性近似（Essenberger 1409.7968 Sec. VI.A 用的）
- 两费米面：Γ 空穴 / M 电子
- 静态 + 动态 RPA 自旋磁化率

Author: Team B (essenberger-chi-reconstruction)
Date: 2026-09-22
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# ============================================================================
# Phase 0: 框架测试（不追求物理参数精度）
# ============================================================================

def chi0_2band_static(q, kF1, kF2, m_star):
    """
    各向同性两带 Lindhard 静态磁化率（T=0 粗糙近似）

    χ₀_αα(q) ~ N_α (DOS at Fermi level)
    χ₀_αβ(q) ~ √(N_α N_β) × overlap_factor

    Args:
        q: momentum transfer (scalar, isotropic)
        kF1, kF2: Fermi wavevectors of band 1, 2
        m_star: effective mass

    Returns:
        chi0: 2×2 static susceptibility matrix
    """
    # DOS at Fermi level (3D isotropic)
    N1 = m_star * kF1 / (2 * np.pi**2)
    N2 = m_star * kF2 / (2 * np.pi**2)

    # Interband overlap factor (phenomenological)
    # Near nesting q ~ |kF1 - kF2|: enhanced
    q_nest = np.abs(kF1 - kF2)
    overlap = 0.3 * np.exp(-((q - q_nest)**2) / (0.5 * q_nest)**2)

    chi0 = np.array([
        [N1, np.sqrt(N1*N2) * overlap],
        [np.sqrt(N1*N2) * overlap, N2]
    ])
    return chi0


def chi_rpa(chi0, U):
    """
    RPA spin susceptibility: χ^s = χ₀(1 - Uχ₀)^{-1}

    Args:
        chi0: bare susceptibility (2×2)
        U: Hubbard interaction (scalar, spin channel)

    Returns:
        chi_s: RPA spin susceptibility (2×2)
    """
    I = np.eye(2)
    try:
        chi_s = chi0 @ np.linalg.inv(I - U * chi0)
    except np.linalg.LinAlgError:
        # Singular: Stoner instability
        chi_s = np.full_like(chi0, np.inf)
    return chi_s


def scan_chi_sign_structure(kF1, kF2, m_star, U, q_range, omega_range=None):
    """
    扫描 χ^s_{12}(q,ω) 的符号结构

    Phase 0: 只做静态 (ω=0)

    Returns:
        results: dict with keys
            - q_grid: momentum grid
            - chi_s_12: off-diagonal susceptibility
            - chi_s_eigs: eigenvalues
    """
    q_grid = np.linspace(q_range[0], q_range[1], 100)
    chi_s_11 = np.zeros(len(q_grid))
    chi_s_12 = np.zeros(len(q_grid))
    chi_s_22 = np.zeros(len(q_grid))
    eig_min = np.zeros(len(q_grid))
    eig_max = np.zeros(len(q_grid))

    for i, q in enumerate(q_grid):
        chi0 = chi0_2band_static(q, kF1, kF2, m_star)
        chi_s = chi_rpa(chi0, U)

        chi_s_11[i] = chi_s[0, 0]
        chi_s_12[i] = chi_s[0, 1]
        chi_s_22[i] = chi_s[1, 1]

        eigvals = np.linalg.eigvalsh(chi_s)
        eig_min[i] = eigvals[0]
        eig_max[i] = eigvals[1]

    results = {
        'q_grid': q_grid,
        'chi_s_11': chi_s_11,
        'chi_s_12': chi_s_12,
        'chi_s_22': chi_s_22,
        'eig_min': eig_min,
        'eig_max': eig_max,
    }
    return results


# ============================================================================
# 测试：单点 RPA
# ============================================================================
if __name__ == '__main__':
    print("=" * 70)
    print("Essenberger χ^s reconstruction - Phase 0 framework test")
    print("=" * 70)

    # 参数（arbitrary units）
    kF1 = 1.0  # Gamma hole pocket
    kF2 = 1.2  # M electron pocket
    m_star = 1.0
    U = 0.8  # tuned near Stoner (U_c ~ 1.0)

    # 单点测试
    q_test = 0.2  # nesting ~ |kF1 - kF2| = 0.2
    chi0 = chi0_2band_static(q_test, kF1, kF2, m_star)
    chi_s = chi_rpa(chi0, U)

    print(f"\nSingle-point test at q = {q_test:.2f}:")
    print(f"chi0 =\n{chi0}")
    print(f"\nchi_s (RPA) =\n{chi_s}")
    print(f"\ndet(chi_s) = {np.linalg.det(chi_s):.6f}")
    eigvals = np.linalg.eigvalsh(chi_s)
    print(f"eigenvalues = {eigvals}")
    print(f"PSD check: lambda_min = {eigvals[0]:.6f} (>0 required)")

    # 符号检查
    print(f"\n*** Eq.(37) check: chi_s_12 = {chi_s[0,1]:.6f}")
    if chi_s[0, 1] > 0:
        print("    --> POSITIVE off-diagonal (consistent with Essenberger)")
    else:
        print("    --> NEGATIVE off-diagonal (contradicts Essenberger Eq.(37))")

    print("\n" + "=" * 70)
    print("q-scan: chi_s_12(q) sign structure")
    print("=" * 70)

    # q 扫描
    results = scan_chi_sign_structure(kF1, kF2, m_star, U, q_range=(0.0, 2.0))

    # 符号普查
    n_negative = np.sum(results['chi_s_12'] < 0)
    n_total = len(results['chi_s_12'])
    print(f"\nOff-diagonal sign census:")
    print(f"  chi_s_12 < 0: {n_negative}/{n_total} points ({100*n_negative/n_total:.1f}%)")
    print(f"  chi_s_12 > 0: {n_total - n_negative}/{n_total} points")

    # PSD 检验
    n_psd = np.sum(results['eig_min'] > 0)
    print(f"\nPSD check:")
    print(f"  lambda_min > 0: {n_psd}/{n_total} points ({100*n_psd/n_total:.1f}%)")

    # 关键问题：负非对角 + PSD 共存？
    mask_neg_offdiag = results['chi_s_12'] < 0
    mask_psd = results['eig_min'] > 0
    n_coexist = np.sum(mask_neg_offdiag & mask_psd)

    print(f"\n*** Critical question: Negative off-diagonal + PSD coexistence?")
    print(f"    {n_coexist}/{n_negative} negative-offdiag points are also PSD")

    if n_coexist > 0:
        print(f"    --> FOUND: negative chi_s_12 + PSD configuration")
        print(f"        (contradicts naive 'repulsive => positive' intuition)")

    print("\nFramework test complete. Ready for Phase 1 (Essenberger parameters).")
