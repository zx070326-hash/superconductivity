#!/usr/bin/env python3
"""
SCDFT Phase 2 Week 2 Day 7: 动力学 Im chi(Q_AFM, omega) 最小验证 — 出路 3

任务（team-lead 派活，Day 7，Week 2 最后一天）：
    Essenberger Eq.(36) 的核是频率积分：∫domega [权重] Im chi(q,omega)。
    Day 4-6 的三次尝试全部工作在静态 (omega=0) 对象上，从未真正对齐
    Essenberger 的定义域。本脚本构建动力学响应，检验：

    1) Im chi(Q_AFM, omega) 的正确定义是什么（不是逐元 Im 部分，是谱函数）
    2) 涨落-耗散定理保证的 PSD 性质是否数值确认
    3) 峰值频率 omega_peak 处符号结构是否比 omega=0 处更好
    4) （若 1-3 成功）频率积分构造 Lambda^SF，检查其符号

关键澄清（本脚本的核心技术点）：
    对于多轨道矩阵响应函数，"Im chi >= 0（PSD）"这个 FDT 陈述，
    正确的数学对象不是逐元虚部 Im(chi_ij)（那是一个实反对称矩阵，
    本征值纯虚数，谈不上正负），而是"谱函数"矩阵：

        A(q,omega) := i * [ chi(q,omega) - chi(q,omega)^dagger ]

    这是 Hermitian 矩阵（i(M-M^dagger) 对任意 M 都是 Hermitian），
    其谱表示 A_ij(omega) = 2*pi * sum_n <0|O_i|n><n|O_j|0> delta(omega-E_n)（omega>0 支路）
    是一个 Gram 矩阵之和 ⟹ PSD，这是涨落-耗散定理的矩阵版本。
    对角元 A_ii(omega) 正是通常写的 "Im chi_ii(omega)"（标量情形下的熟悉表述）。
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from numpy.linalg import eigh, inv

# ============================================================================
# 复用基线函数
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

def lindhard_chi0_omega(q, omega, mu, T=0.01, eta=0.05, k_grid=30):
    """
    动力学 Lindhard 响应（推迟，retarded）：

        chi0_{ij}(q,omega) = sum_{k,m,n} conj(U_k[i,m]) U_{k+q}[j,n]
                              * (f_k[m] - f_{k+q}[n]) / (omega + eig_k[m] - eig_{k+q}[n] + i*eta)

    omega=0 时退化为 Week 1 的静态 lindhard_chi0（此处不做 m==n,q=0 的
    L'Hospital 特判，因为本脚本只在 q=Q_AFM≠0 处使用）。
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

            for m in range(5):
                for n in range(5):
                    denom = omega + eig_k[m] - eig_kq[n] + 1j*eta
                    chi0_band_mn = (f_k[m] - f_kq[n]) / denom

                    for i in range(5):
                        for j in range(5):
                            chi0_orbital[i, j] += (
                                np.conj(U_k[i, m]) * chi0_band_mn * U_kq[j, n]
                                * dk_normalized
                            )

    return chi0_orbital

def rpa_spin_susceptibility(chi0, U_spin):
    I = np.eye(chi0.shape[0])
    return chi0 @ inv(I - U_spin * chi0)

def spectral_function(chi):
    """
    A(q,omega) := i*(chi - chi^dagger)   [Hermitian, PSD 若 FDT 成立]
    对角元 = 通常意义的 2*Im(chi_ii)（标量情形下的熟悉量）
    """
    return 1j * (chi - chi.conj().T)

def analyze_hermitian(M, label="M"):
    herm_residual = np.linalg.norm(M - M.conj().T)
    eigs = eigh(M)[0].real
    n_neg = np.sum(eigs < -1e-10)
    n_pos = np.sum(eigs > 1e-10)
    print(f"{label}: Hermitian residual={herm_residual:.2e}, {n_neg} neg, {n_pos} pos")
    print(f"  eigenvalues: {np.round(eigs, 5)}")
    return n_neg, n_pos, eigs

# ============================================================================
# Day 7 优先级 1：Im chi(Q_AFM, omega) 扫描
# ============================================================================

def scan_dynamic_spectral_function():
    mu = 0.0
    U_spin = 1.0
    q_afm = (np.pi, 0)
    k_grid = 20  # Day 7 时间预算内用较粗网格（半天任务），若需要可事后精修

    # 带宽尺度：t_xz=t_yz=0.5 → 典型能级差量级 ~1-2；扫描到 2*bandwidth
    omega_list = [0.0, 0.2, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0]

    print("=" * 70)
    print("Day 7 优先级 1：动力学 Im chi(Q_AFM, omega) 扫描")
    print("=" * 70)

    results = []
    for omega in omega_list:
        chi0_w = lindhard_chi0_omega(q_afm, omega, mu, k_grid=k_grid)
        chi_s_w = rpa_spin_susceptibility(chi0_w, U_spin)

        A_w = spectral_function(chi_s_w)
        n_neg, n_pos, eigs = analyze_hermitian(A_w, f"A(Q_AFM, omega={omega:.2f})")
        tr_A = np.trace(A_w).real

        print(f"  tr[A(omega={omega:.2f})] = {tr_A:.5f}")
        print()

        results.append({
            'omega': omega, 'n_neg': n_neg, 'n_pos': n_pos,
            'eigs': eigs, 'trace': tr_A
        })

    # 峰值定位
    traces = np.array([r['trace'] for r in results])
    idx_peak = np.argmax(traces)
    omega_peak = results[idx_peak]['omega']

    print("=" * 70)
    print("判决：")
    print(f"  omega=0 (静态): {results[0]['n_neg']} 负, {results[0]['n_pos']} 正, tr[A]={results[0]['trace']:.5f}")
    print(f"  peak omega={omega_peak}: {results[idx_peak]['n_neg']} 负, {results[idx_peak]['n_pos']} 正, tr[A]={traces[idx_peak]:.5f}")
    print()

    all_psd = all(r['n_neg'] == 0 for r in results)
    if all_psd:
        print("  [FDT CONFIRMED] 所有 omega>=0 处 A(q,omega) 全部 PSD（数值确认涨落-耗散定理）")
    else:
        n_bad = sum(1 for r in results if r['n_neg'] > 0)
        print(f"  [FDT VIOLATION?] {n_bad}/{len(results)} 个频率点出现负本征值 —— 需要排查数值/物理前提")

    if results[0]['n_neg'] > results[idx_peak]['n_neg']:
        print("  → 峰值处符号结构确实比 omega=0 更好：支持出路 3 的机制链条")
    elif results[0]['n_neg'] == results[idx_peak]['n_neg'] == 0:
        print("  → 全程 PSD（含 omega=0），Sylvester/静态负本征值问题在这个定义下不出现")
        print("    （提醒：这与 Day 4-6 的静态 Re chi(q,0) 负本征值问题是不同对象——")
        print("     那里测的是 chi^s 实部本身的本征值，这里测的是谱函数 A=i(chi-chi^dagger)）")
    print("=" * 70)

    return results

if __name__ == "__main__":
    results = scan_dynamic_spectral_function()
