#!/usr/bin/env python3
"""
SCDFT Phase 2 Week 2: f^xc 夹心变换（Yukawa 核）

目标：验证 f^xc 夹心能否将 chi^s 的负本征值"翻转"为正

核心方程（Essenberger Eq.(32)）：
    Lambda^SF(q) = 2 ∬ f^xc(r-r') chi^s(r',r'',q) f^xc(r''-r_final) dr' dr''

简化为 Yukawa 核：
    f^xc(r) = -A · exp(-kappa·|r|) / |r|

参数建议（team-lead）：
    A ~ 1-2 eV（交换能量尺度）
    kappa ~ 1-2 Å^-1（铁基晶格常数 a ~ 4 Å）
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from numpy.linalg import eigh, inv

# ============================================================================
# 复用 Week 1 基线函数（避免 import 依赖）
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

def lindhard_chi0(q, mu, T=0.01, eta=0.01, k_grid=30):
    """Lindhard 裸响应函数"""
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

            for m in range(5):
                for n in range(5):
                    denom = eig_k[m] - eig_kq[n] + 1j*eta

                    if abs(qx) < 1e-10 and abs(qy) < 1e-10 and m == n:
                        chi0_band_mn = f_k[m] * (1 - f_k[m]) / T
                    else:
                        chi0_band_mn = (f_k[m] - f_kq[n]) / denom

                    for i in range(5):
                        for j in range(5):
                            chi0_orbital[i, j] += (
                                np.conj(U_k[i, m]) * chi0_band_mn * U_kq[j, n]
                                * dk_normalized
                            )

    return chi0_orbital

def rpa_spin_susceptibility(chi0, U_spin):
    """RPA: chi^s = chi0 (1 - U chi0)^-1"""
    I = np.eye(chi0.shape[0])
    try:
        chi_s = chi0 @ inv(I - U_spin * chi0)
        return chi_s
    except np.linalg.LinAlgError:
        print("WARNING: RPA inversion failed")
        return None

def analyze_eigenvalues(chi, label="chi"):
    """分析本征值符号"""
    eigs = eigh(chi)[0].real
    n_neg = np.sum(eigs < 0)
    n_pos = np.sum(eigs > 0)
    n_zero = np.sum(np.abs(eigs) < 1e-10)

    print(f"{label}: {n_neg} neg, {n_pos} pos, {n_zero} zero")
    print(f"  eigenvalues: {eigs}")
    return n_neg, n_pos, eigs

# ============================================================================
# Yukawa f^xc 核
# ============================================================================

def yukawa_fxc_realspace(r, A=1.5, kappa=1.5):
    """
    Yukawa 交换-关联核（实空间）

    f^xc(r) = -A · exp(-kappa·|r|) / |r|

    参数：
        A: 能量标度（eV）
        kappa: 逆程长度（Å^-1）
    """
    r_norm = np.linalg.norm(r)
    if r_norm < 1e-10:
        # 正则化：r=0 处取极限 -A·kappa
        return -A * kappa
    return -A * np.exp(-kappa * r_norm) / r_norm

def yukawa_fxc_matrix(q, A=1.5, kappa=1.5, lattice_const=4.0, r_cutoff=3):
    """
    Yukawa f^xc 在动量空间的 5×5 轨道矩阵

    通过实空间格点求和后傅里叶变换：
        f^xc(q) = Σ_R f^xc(R) exp(i q·R)

    参数：
        lattice_const: 晶格常数（Å）
        r_cutoff: 实空间截断（晶格常数单位）
    """
    fxc_q = np.zeros((5, 5), dtype=complex)

    # 实空间网格（简化：正方晶格）
    R_max = int(r_cutoff)
    for Rx in range(-R_max, R_max + 1):
        for Ry in range(-R_max, R_max + 1):
            if Rx == 0 and Ry == 0:
                continue  # 跳过 R=0（单独处理）

            R_vec = np.array([Rx, Ry]) * lattice_const
            fxc_R = yukawa_fxc_realspace(R_vec, A, kappa)

            # 傅里叶变换
            phase = np.exp(1j * np.dot(q, R_vec))
            fxc_q += fxc_R * phase * np.eye(5)  # 简化：对角近似

    # R=0 贡献（对角）
    fxc_q += yukawa_fxc_realspace(np.zeros(2), A, kappa) * np.eye(5)

    return fxc_q

# ============================================================================
# f^xc 夹心变换
# ============================================================================

def fxc_sandwich(chi_s, fxc_q):
    """
    夹心变换：Lambda^SF = 2 · f^xc · chi^s · f^xc

    简化：假设 f^xc 对称（实际应取共轭转置）
    """
    Lambda_SF = 2.0 * fxc_q @ chi_s @ fxc_q
    return Lambda_SF

# ============================================================================
# Week 2 Day 3-5：Q_AFM 处符号演化检验
# ============================================================================

def verify_sign_flip_at_qafm():
    """
    关键检验（Week 2 核心）：
    Q_AFM 处，chi^s 有 4 个负本征值（已验证）
    经 f^xc 夹心后，Lambda^SF 是否全正？
    """
    mu = 0.0
    U_spin = 1.0
    q_afm = (np.pi, 0)

    print("=" * 60)
    print("Week 2 Day 3-5：f^xc 夹心符号演化检验")
    print("=" * 60)

    # Step 1：RPA 基线（Week 1 已验证）
    print("\n[Step 1] RPA 基线：chi^s(Q_AFM)")
    chi0 = lindhard_chi0(q_afm, mu, k_grid=30)
    chi_s = rpa_spin_susceptibility(chi0, U_spin)

    if chi_s is None:
        print("ERROR: RPA 反演失败")
        return

    n_neg_chi_s, _, eigs_chi_s = analyze_eigenvalues(chi_s, "chi^s")
    print(f"  预期：4 个负（Week 1 验证）")
    print(f"  实际：{n_neg_chi_s} 个负 [OK]" if n_neg_chi_s == 4 else f"  [FAIL]")

    # Step 2：Yukawa f^xc 核
    print("\n[Step 2] Yukawa f^xc 核")
    A = 1.5  # eV
    kappa = 1.5  # Å^-1
    print(f"  参数：A={A} eV, kappa={kappa} Å^-1")

    fxc_q = yukawa_fxc_matrix(q_afm, A, kappa)
    print(f"  f^xc(Q_AFM) 对角元素: {np.diag(fxc_q).real}")

    # Step 3：夹心变换
    print("\n[Step 3] 夹心变换：Lambda^SF = 2·f^xc·chi^s·f^xc")
    Lambda_SF = fxc_sandwich(chi_s, fxc_q)

    n_neg_Lambda, n_pos_Lambda, eigs_Lambda = analyze_eigenvalues(Lambda_SF, "Lambda^SF")

    # Step 4：判决
    print("\n" + "=" * 60)
    print("判决：")
    print(f"  输入 chi^s: {n_neg_chi_s} 个负本征值")
    print(f"  输出 Lambda^SF: {n_neg_Lambda} 个负本征值")

    if n_neg_Lambda == 0:
        print("\n  [SUCCESS] 符号完全翻转！")
        print("  → 机制 alpha（非局域平滑）成立")
        print("  → Essenberger Eq.(37) Im[Lambda^SF] > 0 验证通过")
    elif n_neg_Lambda < n_neg_chi_s:
        print(f"\n  [PARTIAL] 负本征值减少 {n_neg_chi_s - n_neg_Lambda} 个")
        print("  → f^xc 夹心有效，但不足以完全翻转")
        print("  → 需调参（A, kappa）或尝试 ALDA")
    else:
        print("\n  [FAIL] 负本征值未减少或增加")
        print("  → Yukawa 核不适用")
        print("  → 需切换到 ALDA 或 PGG 近似")

    print("=" * 60)

    return {
        'n_neg_chi_s': n_neg_chi_s,
        'n_neg_Lambda': n_neg_Lambda,
        'eigs_chi_s': eigs_chi_s,
        'eigs_Lambda': eigs_Lambda
    }

# ============================================================================
# Week 2 Day 6-7：参数扫描（若需要）
# ============================================================================

def scan_yukawa_parameters():
    """
    若单点检验未完全翻转，扫描 (A, kappa) 参数空间
    """
    mu = 0.0
    U_spin = 1.0
    q_afm = (np.pi, 0)

    # RPA 基线
    chi0 = lindhard_chi0(q_afm, mu, k_grid=20)  # 降低网格加速
    chi_s = rpa_spin_susceptibility(chi0, U_spin)

    if chi_s is None:
        return

    # 参数网格
    A_range = np.linspace(0.5, 3.0, 6)
    kappa_range = np.linspace(0.5, 3.0, 6)

    print("\n扫描 Yukawa 参数空间...")
    print(f"A ∈ [{A_range[0]}, {A_range[-1]}] eV")
    print(f"kappa ∈ [{kappa_range[0]}, {kappa_range[-1]}] Å^-1")
    print()

    results = []
    for A in A_range:
        for kappa in kappa_range:
            fxc_q = yukawa_fxc_matrix(q_afm, A, kappa)
            Lambda_SF = fxc_sandwich(chi_s, fxc_q)

            eigs = eigh(Lambda_SF)[0].real
            n_neg = np.sum(eigs < 0)

            results.append((A, kappa, n_neg))

            if n_neg == 0:
                print(f"  [FOUND] A={A:.2f}, kappa={kappa:.2f} → 全正 ✓")

    # 统计
    n_success = sum(1 for _, _, n in results if n == 0)
    print(f"\n统计：{n_success}/{len(results)} 组参数实现全正")

    if n_success == 0:
        print("  → Yukawa 核在此参数范围内无法完全翻转符号")
        print("  → 建议切换到 ALDA 近似")

    return results

# ============================================================================
# 主程序
# ============================================================================

if __name__ == "__main__":
    # Week 2 核心任务：Q_AFM 处符号演化
    result = verify_sign_flip_at_qafm()

    # 若未完全翻转，可选参数扫描
    # scan_yukawa_parameters()
