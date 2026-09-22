#!/usr/bin/env python3
"""
SCDFT Phase 2 Week 2 Day 4: 第一层 TDDFT Dyson 重求和

任务（team-lead 派活）：
    检验"对象错位假说"（出路 2）——
    Team B 的 chi^s（Hubbard-U RPA，标量 U，q 无关）
    是否与 Essenberger 的 P（TDDFT 第一层重求和，f^xc(q) 矩阵，q 依赖）
    是不同的对象（符号结构不同）？

核心区别：
    chi^s = chi0 @ inv(I - U*chi0)          # U: 标量, q 无关
    P     = chi0 @ inv(I - fxc(q) @ chi0)   # fxc(q): 矩阵, q 依赖（Yukawa 傅里叶变换）

注意：矩阵求逆是非线性操作，不是 congruence transform，
不受 Sylvester 惯性定律约束——P 的符号结构原则上可以与 chi0/chi^s 不同。

验证点：P(Q_AFM) 的负本征值个数 vs chi^s(Q_AFM) 的负本征值个数（4 个，Week 1 已验证）
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from numpy.linalg import eigh, inv, cond

# ============================================================================
# 复用 Week 1/2 基线函数
# ============================================================================

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

def lindhard_chi0(q, mu, T=0.01, eta=0.01, k_grid=30):
    """裸 KS 响应 chi^KS（Lindhard 函数）"""
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

def rpa_spin_susceptibility_hubbard(chi0, U_spin):
    """Team B 的 chi^s：标量 Hubbard U，q 无关"""
    I = np.eye(chi0.shape[0])
    try:
        chi_s = chi0 @ inv(I - U_spin * chi0)
        return chi_s
    except np.linalg.LinAlgError:
        print("WARNING: Hubbard-RPA inversion failed")
        return None

def analyze_eigenvalues(chi, label="chi"):
    eigs = eigh(chi)[0].real
    n_neg = np.sum(eigs < 0)
    n_pos = np.sum(eigs > 0)
    n_zero = np.sum(np.abs(eigs) < 1e-10)
    print(f"{label}: {n_neg} neg, {n_pos} pos, {n_zero} zero")
    print(f"  eigenvalues: {eigs}")
    return n_neg, n_pos, eigs

# ============================================================================
# q 依赖 Yukawa f^xc（TDDFT 第一层核）
# ============================================================================

def yukawa_fxc_realspace(r, A=1.5, kappa=1.5):
    r_norm = np.linalg.norm(r)
    if r_norm < 1e-10:
        return -A * kappa
    return -A * np.exp(-kappa * r_norm) / r_norm

def yukawa_fxc_matrix(q, A=1.5, kappa=1.5, lattice_const=4.0, r_cutoff=3):
    """
    Yukawa f^xc(q)：q 依赖矩阵（本实现仍是对角/单位阵比例，
    但数值随 q 变化——这是与标量 U 的关键区别）
    """
    fxc_q = np.zeros((5, 5), dtype=complex)

    R_max = int(r_cutoff)
    for Rx in range(-R_max, R_max + 1):
        for Ry in range(-R_max, R_max + 1):
            if Rx == 0 and Ry == 0:
                continue
            R_vec = np.array([Rx, Ry]) * lattice_const
            fxc_R = yukawa_fxc_realspace(R_vec, A, kappa)
            phase = np.exp(1j * np.dot(q, R_vec))
            fxc_q += fxc_R * phase * np.eye(5)

    fxc_q += yukawa_fxc_realspace(np.zeros(2), A, kappa) * np.eye(5)
    return fxc_q

# ============================================================================
# 第一层 TDDFT Dyson 重求和：P(q) = chi0 [1 - fxc(q) chi0]^-1
# ============================================================================

def tddft_dyson_proper_response(chi0, fxc_q):
    """
    Eq.(26)-(29) 的第一层重求和：

    P(q) = chi^KS(q) / (1 - f^xc(q) chi^KS(q))
         = chi0 @ inv(I - fxc_q @ chi0)

    这是矩阵求逆（非线性），不是 congruence transform，
    不受 Sylvester 惯性定律约束。
    """
    I = np.eye(chi0.shape[0])
    denom = I - fxc_q @ chi0
    kappa_cond = cond(denom)

    try:
        P = chi0 @ inv(denom)
        return P, kappa_cond
    except np.linalg.LinAlgError:
        print("WARNING: TDDFT-Dyson inversion failed (Stoner instability in f^xc channel?)")
        return None, kappa_cond

# ============================================================================
# Day 4 核心检验
# ============================================================================

def compare_hubbard_vs_tddft_at_qafm():
    mu = 0.0
    U_hubbard = 1.0
    q_afm = (np.pi, 0)

    print("=" * 70)
    print("Week 2 Day 4: Hubbard-U chi^s vs TDDFT-Dyson P 对比（Q_AFM）")
    print("=" * 70)

    # 裸响应（两者共用同一个 chi0=chi^KS）
    print("\n[Step 0] 裸响应 chi0 = chi^KS(Q_AFM)")
    chi0 = lindhard_chi0(q_afm, mu, k_grid=30)
    analyze_eigenvalues(chi0, "chi0/chi^KS")

    # 路线 A：Team B 的 Hubbard-U RPA（标量 U，q 无关）
    print("\n[路线 A] chi^s = chi0 / (1 - U*chi0)，U=1.0（标量，q 无关）")
    chi_s = rpa_spin_susceptibility_hubbard(chi0, U_hubbard)
    n_neg_hubbard, _, eigs_hubbard = analyze_eigenvalues(chi_s, "chi^s (Hubbard)")

    # 路线 B：TDDFT 第一层重求和（f^xc(q) 矩阵，q 依赖）
    print("\n[路线 B] P = chi0 / (1 - fxc(q)*chi0)，fxc(q)=Yukawa 傅里叶变换")
    A_param, kappa_param = 1.5, 1.5
    fxc_q = yukawa_fxc_matrix(q_afm, A_param, kappa_param)
    print(f"  fxc(Q_AFM) 对角元素: {np.diag(fxc_q).real}")

    P, cond_num = tddft_dyson_proper_response(chi0, fxc_q)
    print(f"  [1 - fxc*chi0] 条件数: {cond_num:.4e}")

    if P is not None:
        n_neg_tddft, _, eigs_tddft = analyze_eigenvalues(P, "P (TDDFT)")
    else:
        n_neg_tddft = None

    # 判决
    print("\n" + "=" * 70)
    print("判决：")
    print(f"  Hubbard-U chi^s: {n_neg_hubbard} 个负本征值")
    if n_neg_tddft is not None:
        print(f"  TDDFT-Dyson P:   {n_neg_tddft} 个负本征值")

        if n_neg_tddft < n_neg_hubbard:
            print(f"\n  [PARTIAL SUCCESS] P 比 chi^s 更正（少 {n_neg_hubbard - n_neg_tddft} 个负本征值）")
            print("  → 出路 2（对象错位假说）部分成立")
            print("  → 第一层 TDDFT Dyson 重求和确实改善了符号结构")
        elif n_neg_tddft == 0:
            print("\n  [SUCCESS] P 完全正定！")
            print("  → 出路 2 完全成立：Team B 的 chi^s 不是正确的 Essenberger 输入对象")
        elif n_neg_tddft == n_neg_hubbard:
            print("\n  [NO CHANGE] P 与 chi^s 符号结构相同")
            print("  → 出路 2 不成立（至少此 Yukawa 参数下）")
            print("  → 需要考虑出路 1（f_σ^F 奇异）或出路 3（频率积分语境）")
        else:
            print("\n  [WORSE] P 比 chi^s 更负")
            print("  → 出路 2 不成立，且方向相反")

    print("=" * 70)

    return {
        'n_neg_hubbard': n_neg_hubbard,
        'n_neg_tddft': n_neg_tddft,
        'eigs_hubbard': eigs_hubbard,
        'eigs_tddft': eigs_tddft if P is not None else None,
        'cond_num': cond_num
    }

# ============================================================================
# 主程序
# ============================================================================

if __name__ == "__main__":
    result = compare_hubbard_vs_tddft_at_qafm()
