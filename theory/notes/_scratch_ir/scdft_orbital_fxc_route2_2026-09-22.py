#!/usr/bin/env python3
"""
SCDFT Phase 2 Week 2 Day 5: 真正的轨道结构 f^xc(q) — 重新检验出路 2

背景（订正见 day4_scalar_fxc_artifact_retraction_2026-09-22.md）：
    Day 4 的 yukawa_fxc_matrix() 每个傅里叶项都乘 np.eye(5)，
    导致 f^xc(q) 恒为 标量函数(q) * I，与 Hubbard-U（也是标量）
    数学等价 —— 不构成"对象错位假说"（出路 2）的检验。

本脚本的改动（唯一实质性改动）：
    f^xc_{ij}(q) = sum_R A_{ij} * exp(-kappa_{ij}|R|)/|R| * exp(i q.R)
    其中耦合矩阵 A_{ij}, kappa_{ij} 依轨道对不同：
      - 对角元（同轨道 Hubbard-型 onsite 排斥）：A_ii ~ 1.5-2.0
      - (xz,yz) 非对角（Hund's 型自旋关联，简并 t2g 对）：A_01 ~ 0.6
      - 其余非对角（跨对称性，弱）：A_ij ~ 0.15-0.3
    这给 f^xc(q) 提供了 chi0/chi^s 之外的独立矩阵结构（非对角元不为零、
    且非对角元与对角元的比例本身是 q 的非平凡函数——因为每个 (i,j) 分量
    有自己的 kappa_{ij}，衰减速度不同，傅里叶变换后 q 依赖形状也不同）。

检验：P(Q_AFM) = chi0 @ inv(I - fxc(q) @ chi0) 的负本征值个数
      vs chi^s(Q_AFM)（Hubbard-U，标量，4 个负，Week 1 已验证）
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from numpy.linalg import eigh, inv, cond

# ============================================================================
# 复用基线函数（Week 1/2 一致）
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
    I = np.eye(chi0.shape[0])
    chi_s = chi0 @ inv(I - U_spin * chi0)
    return chi_s

def analyze_eigenvalues(chi, label="chi"):
    eigs = eigh(chi)[0].real
    n_neg = np.sum(eigs < 0)
    n_pos = np.sum(eigs > 0)
    n_zero = np.sum(np.abs(eigs) < 1e-10)
    print(f"{label}: {n_neg} neg, {n_pos} pos, {n_zero} zero")
    print(f"  eigenvalues: {eigs}")
    return n_neg, n_pos, eigs

# ============================================================================
# 轨道结构 f^xc(q)：真正的非对角矩阵
# ============================================================================

# 轨道标号：0=xz, 1=yz, 2=xy, 3=(x2-y2 型), 4=(z2 型 / eg 剩余)
# 耦合矩阵 A_{ij}（对称）：对角=onsite Hubbard 型，非对角=Hund's/轨道混合
A_MATRIX = np.array([
    [1.8, 0.6, 0.25, 0.10, 0.05],   # xz
    [0.6, 1.8, 0.25, 0.10, 0.05],   # yz  (与 xz 简并对，强 Hund's 耦合)
    [0.25, 0.25, 1.2, 0.15, 0.05],  # xy
    [0.10, 0.10, 0.15, 0.8, 0.10],  # 其余
    [0.05, 0.05, 0.05, 0.10, 0.6],  # 其余
])

# 衰减长度矩阵 kappa_{ij}：非对角项衰减更快（短程 Hund's 关联）
KAPPA_MATRIX = np.array([
    [1.2, 2.0, 2.5, 3.0, 3.0],
    [2.0, 1.2, 2.5, 3.0, 3.0],
    [2.5, 2.5, 1.5, 3.0, 3.0],
    [3.0, 3.0, 3.0, 2.0, 3.0],
    [3.0, 3.0, 3.0, 3.0, 2.0],
])

assert np.allclose(A_MATRIX, A_MATRIX.T), "A_MATRIX 必须对称（Hermitian f^xc 前提）"
assert np.allclose(KAPPA_MATRIX, KAPPA_MATRIX.T), "KAPPA_MATRIX 必须对称"

def yukawa_fxc_realspace_orbital(r_norm, i, j):
    """轨道对 (i,j) 各自的 Yukawa 衰减，实空间"""
    A_ij = A_MATRIX[i, j]
    kappa_ij = KAPPA_MATRIX[i, j]
    if r_norm < 1e-10:
        return -A_ij * kappa_ij
    return -A_ij * np.exp(-kappa_ij * r_norm) / r_norm

def orbital_fxc_matrix(q, lattice_const=4.0, r_cutoff=3):
    """
    真正非对角的 f^xc(q)：每个轨道对 (i,j) 有独立的傅里叶求和，
    不再是 标量(q) * I。
    """
    fxc_q = np.zeros((5, 5), dtype=complex)

    R_max = int(r_cutoff)
    R_list = []
    for Rx in range(-R_max, R_max + 1):
        for Ry in range(-R_max, R_max + 1):
            if Rx == 0 and Ry == 0:
                continue
            R_list.append(np.array([Rx, Ry]) * lattice_const)

    for i in range(5):
        for j in range(5):
            acc = 0.0 + 0.0j
            for R_vec in R_list:
                r_norm = np.linalg.norm(R_vec)
                fxc_R = yukawa_fxc_realspace_orbital(r_norm, i, j)
                phase = np.exp(1j * np.dot(q, R_vec))
                acc += fxc_R * phase
            acc += yukawa_fxc_realspace_orbital(0.0, i, j)  # R=0 项
            fxc_q[i, j] = acc

    # 强制 Hermitian（数值上应已近似成立，因 A/kappa 矩阵对称 + phase(-R)=phase(R)^*）
    fxc_q = 0.5 * (fxc_q + fxc_q.conj().T)
    return fxc_q

def tddft_dyson_proper_response(chi0, fxc_q):
    I = np.eye(chi0.shape[0])
    denom = I - fxc_q @ chi0
    kappa_cond = cond(denom)
    P = chi0 @ inv(denom)
    return P, kappa_cond

# ============================================================================
# 检验
# ============================================================================

def compare_hubbard_vs_orbital_tddft_at_qafm():
    mu = 0.0
    U_hubbard = 1.0
    q_afm = (np.pi, 0)

    print("=" * 70)
    print("Week 2 Day 5: 轨道结构 f^xc(q) vs Hubbard-U chi^s（Q_AFM）")
    print("=" * 70)

    print("\n[Step 0] 裸响应 chi0 = chi^KS(Q_AFM)")
    chi0 = lindhard_chi0(q_afm, mu, k_grid=30)
    analyze_eigenvalues(chi0, "chi0/chi^KS")

    print("\n[路线 A] chi^s = chi0 / (1 - U*chi0)，U=1.0（标量基线，Week 1）")
    chi_s = rpa_spin_susceptibility_hubbard(chi0, U_hubbard)
    n_neg_hubbard, _, eigs_hubbard = analyze_eigenvalues(chi_s, "chi^s (Hubbard)")

    print("\n[路线 B'] P = chi0 / (1 - fxc(q)*chi0)，fxc(q) 为真正轨道矩阵")
    fxc_q = orbital_fxc_matrix(q_afm)
    print("  fxc(Q_AFM) 矩阵（实部）：")
    print(np.round(fxc_q.real, 4))
    off_diag_norm = np.linalg.norm(fxc_q - np.diag(np.diag(fxc_q)))
    diag_norm = np.linalg.norm(np.diag(fxc_q))
    print(f"  非对角/对角 Frobenius 范数比 = {off_diag_norm/diag_norm:.4f}  (非零 ⟹ 真实非对角结构)")
    print(f"  Hermitian 残差 (应为 0): {np.linalg.norm(fxc_q - fxc_q.conj().T):.2e}")

    P, cond_num = tddft_dyson_proper_response(chi0, fxc_q)
    print(f"  [1 - fxc*chi0] 条件数: {cond_num:.4e}")
    n_neg_orbital, _, eigs_orbital = analyze_eigenvalues(P, "P (orbital-fxc TDDFT)")

    print("\n" + "=" * 70)
    print("判决：")
    print(f"  Hubbard-U chi^s (标量):        {n_neg_hubbard} 个负本征值")
    print(f"  P (真正轨道结构 fxc, TDDFT):   {n_neg_orbital} 个负本征值")

    if n_neg_orbital < n_neg_hubbard:
        print(f"\n  [ROUTE 2 SUPPORTED] 真正轨道结构使负本征值减少 {n_neg_hubbard - n_neg_orbital} 个")
        print("  → 对象错位假说（出路 2）在此参数化下获得初步支持")
        print("  → 仍需：(a) 排除本次 A/kappa 参数选择的偶然性 (b) 对多个 q 点复核")
    elif n_neg_orbital == n_neg_hubbard:
        print("\n  [NO CHANGE] 轨道结构未改变符号计数（可能仍是幅度效应，需换参数/对角占优检验）")
    else:
        print("\n  [WORSE] 轨道结构使负本征值增多")
        print("  → 出路 2 在此参数化下不成立")

    print("=" * 70)

    return {
        'n_neg_hubbard': n_neg_hubbard,
        'n_neg_orbital': n_neg_orbital,
        'eigs_hubbard': eigs_hubbard,
        'eigs_orbital': eigs_orbital,
        'cond_num': cond_num,
        'off_diag_norm': off_diag_norm,
        'diag_norm': diag_norm,
    }

def robustness_scan_ABC():
    """
    排除"参数偶然性"：在若干独立的 (A,kappa) 矩阵扰动下重复检验，
    看 n_neg_orbital < n_neg_hubbard 是否稳健，还是像 Day 4 一样对参数敏感。
    """
    mu = 0.0
    U_hubbard = 1.0
    q_afm = (np.pi, 0)

    chi0 = lindhard_chi0(q_afm, mu, k_grid=20)  # 降低网格加速扫描
    chi_s = rpa_spin_susceptibility_hubbard(chi0, U_hubbard)
    n_neg_hubbard = np.sum(eigh(chi_s)[0].real < 0)

    print("\n" + "=" * 70)
    print("稳健性扫描：随机扰动 A/kappa 矩阵（保持对称、保持对角>非对角量级关系）")
    print("=" * 70)

    rng = np.random.default_rng(42)
    global A_MATRIX, KAPPA_MATRIX
    A_base = A_MATRIX.copy()
    kappa_base = KAPPA_MATRIX.copy()

    n_trials = 12
    results = []
    for trial in range(n_trials):
        # 保持对称扰动
        dA = rng.uniform(-0.3, 0.3, size=(5, 5))
        dA = 0.5 * (dA + dA.T)
        dK = rng.uniform(-0.5, 0.5, size=(5, 5))
        dK = 0.5 * (dK + dK.T)

        A_MATRIX = np.clip(A_base + dA, 0.02, None)
        KAPPA_MATRIX = np.clip(kappa_base + dK, 0.3, None)

        fxc_q = orbital_fxc_matrix(q_afm)
        P, _ = tddft_dyson_proper_response(chi0, fxc_q)
        n_neg_P = np.sum(eigh(P)[0].real < 0)
        results.append(n_neg_P)
        print(f"  trial {trial:2d}: n_neg(P) = {n_neg_P}   (chi^s baseline = {n_neg_hubbard})")

    A_MATRIX = A_base
    KAPPA_MATRIX = kappa_base

    results = np.array(results)
    print(f"\n统计：n_neg(P) 范围 [{results.min()}, {results.max()}], 均值 {results.mean():.2f}")
    print(f"chi^s (Hubbard) 基线固定于 {n_neg_hubbard}")
    frac_better = np.mean(results < n_neg_hubbard)
    print(f"改善比例（n_neg(P) < n_neg(chi^s)）: {frac_better:.2%}")

    if frac_better > 0.7:
        print("→ 结果对参数扰动稳健：轨道结构本身（不是某个特定数值）驱动改善")
    elif frac_better < 0.3:
        print("→ 结果不稳健：多数扰动下无改善，Day 5 单点结果可能仍是参数巧合")
    else:
        print("→ 结果部分稳健：改善方向存在但幅度/发生率依赖参数，需要进一步机制分析")

    return results

if __name__ == "__main__":
    result = compare_hubbard_vs_orbital_tddft_at_qafm()
    scan_results = robustness_scan_ABC()
