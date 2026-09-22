#!/usr/bin/env python3
"""
SCDFT Phase 2 Week 1: Lindhard + RPA 基线复现

目标：验证 Team B 的「100% BZ 负本征值」结论
- q=(0,0): 2 个负本征值
- q=Q_AFM: 4 个负本征值
- 全 BZ 扫描：0% 改善率

基于 Team B 的 _fese_5orbital_FINAL_2026-09-22.py（已修复三处问题）
"""

import numpy as np
from numpy.linalg import eigh, inv
import matplotlib.pyplot as plt

# ============================================================================
# 5-轨道 FeSe 紧束缚模型（Team B FINAL 版本）
# ============================================================================

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

# ============================================================================
# Lindhard χ₀（Team B FINAL 修正版）
# ============================================================================

def lindhard_chi0(q, mu, T=0.01, eta=0.01, k_grid=30):
    """
    Lindhard 裸响应函数 χ₀(q,ω=0)

    χ₀(q) = Σ_k [f_k - f_{k+q}] / [ε_k - ε_{k+q} + iη] × dk/(2π)²

    返回：5×5 轨道空间矩阵
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

                    # q=0 且 m=n：L'Hospital
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

# ============================================================================
# RPA 自旋磁化率 χ^s
# ============================================================================

def rpa_spin_susceptibility(chi0, U_spin):
    """
    RPA: χ^s = χ₀ (1 - U χ₀)^{-1}

    U_spin: 标量（简化）或 5×5 矩阵
    """
    I = np.eye(chi0.shape[0])
    try:
        chi_s = chi0 @ inv(I - U_spin * chi0)
        return chi_s
    except np.linalg.LinAlgError:
        print("WARNING: RPA 反演失败（Stoner 不稳定性？）")
        return None

# ============================================================================
# 本征值符号分析
# ============================================================================

def analyze_eigenvalues(chi, label="χ"):
    """分析 χ 的本征值符号"""
    eigs = eigh(chi)[0].real
    n_neg = np.sum(eigs < 0)
    n_pos = np.sum(eigs > 0)
    n_zero = np.sum(np.abs(eigs) < 1e-10)

    print(f"{label}: {n_neg} 负, {n_pos} 正, {n_zero} 零")
    print(f"  本征值: {eigs}")
    return n_neg, n_pos, eigs

# ============================================================================
# Week 1 检验点：复现 Team B 基线
# ============================================================================

def verify_team_b_baseline():
    """
    检验点（Week 1 结束时）：
    - q=(0,0): 2 个负本征值
    - q=Q_AFM=(π,0): 4 个负本征值
    """
    mu = 0.0  # 费米能（简化：调到零能级）
    U_spin = 1.0  # Hubbard U（典型值）

    print("=" * 60)
    print("Week 1 检验点：复现 Team B 基线")
    print("=" * 60)

    # Test 1: q=(0,0)
    print("\n[Test 1] q=(0,0)")
    chi0_q0 = lindhard_chi0((0, 0), mu, k_grid=30)
    chi_s_q0 = rpa_spin_susceptibility(chi0_q0, U_spin)

    print("  Lindhard chi0(q=0):")
    n_neg_chi0, _, _ = analyze_eigenvalues(chi0_q0, "  chi0")

    if chi_s_q0 is not None:
        print("  RPA chi^s(q=0):")
        n_neg_rpa, _, _ = analyze_eigenvalues(chi_s_q0, "  chi^s")

        # 预期：Team B 报告 2 个负本征值
        if n_neg_rpa == 2:
            print("  [OK] 与 Team B 一致（2 负）")
        else:
            print(f"  [MISMATCH] 与 Team B 不一致（预期 2 负，实得 {n_neg_rpa} 负）")

    # Test 2: q=Q_AFM=(π,0)
    print("\n[Test 2] q=Q_AFM=(pi,0)")
    chi0_qafm = lindhard_chi0((np.pi, 0), mu, k_grid=30)
    chi_s_qafm = rpa_spin_susceptibility(chi0_qafm, U_spin)

    print("  Lindhard chi0(Q_AFM):")
    analyze_eigenvalues(chi0_qafm, "  chi0")

    if chi_s_qafm is not None:
        print("  RPA chi^s(Q_AFM):")
        n_neg_rpa_afm, _, _ = analyze_eigenvalues(chi_s_qafm, "  chi^s")

        # 预期：Team B 报告 4 个负本征值
        if n_neg_rpa_afm == 4:
            print("  [OK] 与 Team B 一致（4 负）")
        else:
            print(f"  [MISMATCH] 与 Team B 不一致（预期 4 负，实得 {n_neg_rpa_afm} 负）")

    print("\n" + "=" * 60)
    print("检验点状态：")
    print("  若两个测试都 [OK] -> Week 1 完成，进入 Week 2（f^xc 夹心）")
    print("  若任一 [MISMATCH] -> 调试参数（mu, U_spin, k_grid）直到对齐")
    print("=" * 60)

# ============================================================================
# BZ 扫描（可选，Week 1 剩余时间）
# ============================================================================

def scan_bz_negative_eigenvalues(q_grid=8):
    """
    扫描全 BZ，统计负本征值出现率

    Team B 报告：100% BZ 有负本征值
    """
    mu = 0.0
    U_spin = 1.0

    qx_arr = np.linspace(0, np.pi, q_grid)
    qy_arr = np.linspace(0, np.pi, q_grid)

    total_points = 0
    negative_points = 0

    print(f"\n扫描 BZ（{q_grid}×{q_grid} 网格）...")

    for qx in qx_arr:
        for qy in qy_arr:
            chi0 = lindhard_chi0((qx, qy), mu, k_grid=20)  # 降低 k_grid 加速
            chi_s = rpa_spin_susceptibility(chi0, U_spin)

            if chi_s is not None:
                eigs = eigh(chi_s)[0].real
                has_negative = np.any(eigs < 0)

                total_points += 1
                if has_negative:
                    negative_points += 1

    rate = negative_points / total_points * 100 if total_points > 0 else 0
    print(f"负本征值出现率: {negative_points}/{total_points} = {rate:.1f}%")
    print(f"预期（Team B）: 100%")

    if rate > 95:
        print("✓ 与 Team B 一致（>95% 有负本征值）")
    else:
        print(f"✗ 与 Team B 不一致（预期 100%，实得 {rate:.1f}%）")

    return rate

# ============================================================================
# 主程序
# ============================================================================

if __name__ == "__main__":
    # Week 1 核心任务：复现 Team B 基线
    verify_team_b_baseline()

    # 可选：全 BZ 扫描（若上面两个测试通过）
    # scan_bz_negative_eigenvalues(q_grid=8)
