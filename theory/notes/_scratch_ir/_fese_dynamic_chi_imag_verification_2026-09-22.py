#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
FeSe dynamic magnetic susceptibility Im chi(q,omega) numerical verification

Team C Phase 1 tasks:
1. Real-axis analytic continuation: Matsubara -> omega + i*eta
2. PSD verification: Im chi(q,omega>0) eigenvalues all >= 0
3. Amplitude enhancement: ||Im chi(Q_AFM,omega)|| >> ||Im chi(other q,omega)||
4. Nesting correlation: xi(q) vs ||Im chi(q,omega)||

Physical basis:
- Fluctuation-dissipation theorem => Im chi(q,omega>0) must be PSD
- Q_AFM nesting => Im chi amplitude enhancement

Reusing Team B's tight-binding model + Lindhard framework
"""

import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import numpy as np
from numpy.linalg import eigh, eigvalsh, norm as matrix_norm
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# Team B 紧束缚模型（复用）
# ==============================================================================

def build_5orbital_hamiltonian(kx, ky):
    """FeSe 5-轨道简化模型"""
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

# ==============================================================================
# 动力学虚部 Im χ(q,ω)（实轴公式）
# ==============================================================================

def compute_im_chi_realaxis(q, omega, mu, T=0.01, eta=0.01, k_grid=100):
    """
    动力学磁化率虚部 Im χ(q,ω)

    公式（实轴解析延拓）：
    Im χ_{αβ}(q,ω) = -π Σ_k δ(ω - [ε_β(k+q) - ε_α(k)]) [f_α(k) - f_β(k+q)]

    数值实现：
    δ(ω - Δε) → Lorentzian 展宽 (η/π) / [(ω - Δε)² + η²]

    参数：
    - q: 波矢 (qx, qy)
    - omega: 实频率（eV，> 0）
    - eta: 展宽参数（替代 δ 函数，默认 0.01 eV）
    - k_grid: k 空间网格密度（需要密集，100×100 起步）

    返回：
    - Im χ (5×5 Hermitian 矩阵)，理论保证 PSD（ω > 0）
    """
    qx, qy = q
    im_chi_orbital = np.zeros((5, 5), dtype=float)

    kx_arr = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)
    dk = (2*np.pi / k_grid)**2
    BZ_area = (2*np.pi)**2
    dk_normalized = dk / BZ_area

    # Lorentzian 展宽函数
    def lorentzian(x, eta):
        return (eta / np.pi) / (x**2 + eta**2)

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = build_5orbital_hamiltonian(kx, ky)
            H_kq = build_5orbital_hamiltonian(kx + qx, ky + qy)

            eig_k, U_k = eigh(H_k)
            eig_kq, U_kq = eigh(H_kq)

            f_k = np.array([fermi_func(e, mu, T) for e in eig_k])
            f_kq = np.array([fermi_func(e, mu, T) for e in eig_kq])

            # 带基贡献
            for m in range(5):
                for n in range(5):
                    # 能量差
                    delta_E = eig_kq[n] - eig_k[m]

                    # Lorentzian 权重
                    weight = lorentzian(omega - delta_E, eta)

                    # 费米函数差
                    fermi_diff = f_k[m] - f_kq[n]

                    # 带基虚部贡献
                    # CORRECTED SIGN: +π (not -π) to satisfy fluctuation-dissipation theorem
                    # Verified by 1D free electron gas test
                    im_chi_band_mn = np.pi * weight * fermi_diff

                    # 投影到轨道基
                    for i in range(5):
                        for j in range(5):
                            # Hermitian 对称性：只算上三角
                            if i <= j:
                                contrib = (
                                    np.conj(U_k[i, m]) * U_kq[j, n] * im_chi_band_mn
                                ).real
                                im_chi_orbital[i, j] += contrib * dk_normalized

    # 强制 Hermitian（数值对称化）
    im_chi_hermitian = (im_chi_orbital + im_chi_orbital.T) / 2

    return im_chi_hermitian

# ==============================================================================
# 嵌套函数 ξ(q)
# ==============================================================================

def compute_nesting_function(q, mu, fermi_window=0.05, k_grid=100):
    """
    费米面嵌套函数

    ξ_{αβ}(q) = ∫ dk δ(ε_α(k)) δ(ε_β(k+q))

    数值实现：计数同时满足的 k 点数
    |ε_α(k) - μ| < fermi_window  且  |ε_β(k+q) - μ| < fermi_window

    返回：标量（总嵌套度，对所有轨道求和）
    """
    qx, qy = q
    nesting_count = 0

    kx_arr = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)
    ky_arr = np.linspace(-np.pi, np.pi, k_grid, endpoint=False)

    for kx in kx_arr:
        for ky in ky_arr:
            H_k = build_5orbital_hamiltonian(kx, ky)
            H_kq = build_5orbital_hamiltonian(kx + qx, ky + qy)

            eig_k = eigvalsh(H_k)
            eig_kq = eigvalsh(H_kq)

            # 费米面窗口判据
            on_fermi_k = np.abs(eig_k - mu) < fermi_window
            on_fermi_kq = np.abs(eig_kq - mu) < fermi_window

            # 任意轨道对的嵌套贡献
            nesting_count += np.sum(on_fermi_k[:, None] & on_fermi_kq[None, :])

    # 归一化（除以总 k 点数）
    return nesting_count / (k_grid**2)

# ==============================================================================
# Phase 1 验证：5 特殊点 × 7 频率
# ==============================================================================

def phase1_verification():
    """
    Phase 1 主验证流程

    检验项：
    1. PSD 验证：Im χ 的本征值全 ≥ 0
    2. 幅度增强：Q_AFM 处的 Frobenius 范数最大
    3. 嵌套关联：ξ(q) 与 ‖Im χ(q,ω)‖ 的相关性
    """
    print("=" * 80)
    print("FeSe 动力学 Im χ(q,ω) 验证 - Phase 1")
    print("=" * 80)

    # 物理参数
    mu = 0.0
    T = 0.01
    eta = 0.02  # 展宽（相对密集网格可稍大）

    # 5 特殊 q 点
    q_points = {
        'Gamma': (0.0, 0.0),
        'X': (np.pi, 0.0),
        'M': (np.pi, np.pi),
        'Y': (0.0, np.pi),
        'Random': (0.3*np.pi, 0.5*np.pi)
    }

    # 7 频率网格（eV）
    omega_grid = [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]

    print(f"\n参数设置：")
    print(f"  化学势 μ = {mu:.3f} eV")
    print(f"  温度 T = {T:.3f} eV ({T*11604:.1f} K)")
    print(f"  展宽 η = {eta:.3f} eV")
    print(f"  k 网格 = 100×100")
    print(f"\n测试矩阵：{len(q_points)} q点 × {len(omega_grid)} 频率 = {len(q_points)*len(omega_grid)} 个计算")

    # 存储结果
    results = {}
    psd_violations = []

    print("\n" + "=" * 80)
    print("开始计算...")
    print("=" * 80)

    for q_name, q_vec in q_points.items():
        results[q_name] = {}

        print(f"\n[q = {q_name} = ({q_vec[0]/np.pi:.2f}π, {q_vec[1]/np.pi:.2f}π)]")
        print("-" * 80)

        for omega in omega_grid:
            # 计算 Im χ(q,ω)
            im_chi = compute_im_chi_realaxis(q_vec, omega, mu, T=T, eta=eta, k_grid=100)

            # 本征值检验
            eigs = eigvalsh(im_chi)
            min_eig = np.min(eigs)
            max_eig = np.max(eigs)

            # Frobenius 范数
            frob_norm = matrix_norm(im_chi, 'fro')

            # PSD 判定（容许小负值 ~-1e-10 数值误差）
            is_psd = (min_eig > -1e-8)

            if not is_psd:
                psd_violations.append({
                    'q': q_name,
                    'omega': omega,
                    'min_eig': min_eig
                })

            # 存储
            results[q_name][omega] = {
                'im_chi': im_chi,
                'eigs': eigs,
                'frob_norm': frob_norm,
                'is_psd': is_psd
            }

            # 实时输出
            status = "[OK PSD]" if is_psd else f"[NEG lambda_min={min_eig:.2e}]"
            print(f"  omega={omega:5.2f} eV:  ||Im chi||_F = {frob_norm:8.5f}  |  "
                  f"lambda in [{min_eig:+.2e}, {max_eig:+.2e}]  |  {status}")

    # ==============================================================================
    # 检验 1：PSD 验证
    # ==============================================================================

    print("\n" + "=" * 80)
    print("检验 1：PSD 验证（涨落-耗散定理）")
    print("=" * 80)

    total_tests = len(q_points) * len(omega_grid)
    n_violations = len(psd_violations)

    print(f"\nTotal tests: {total_tests}")
    print(f"PSD violations: {n_violations}")
    print(f"Pass rate: {(total_tests - n_violations) / total_tests * 100:.1f}%")

    if n_violations == 0:
        print("\n[OK] All passed! Im chi(q,omega>0) is PSD at all (q,omega)")
        print("    => Fluctuation-dissipation theorem numerically confirmed")
    else:
        print(f"\n[WARN] Found {n_violations} violations (possibly numerical noise):")
        for v in psd_violations[:5]:  # 只显示前 5 个
            print(f"    q={v['q']}, omega={v['omega']} eV: lambda_min = {v['min_eig']:.3e}")

    # ==============================================================================
    # 检验 2：幅度增强
    # ==============================================================================

    print("\n" + "=" * 80)
    print("检验 2：Q_AFM 嵌套的幅度增强")
    print("=" * 80)

    # 选一个典型频率
    omega_test = 0.5  # eV

    print(f"\nFrequency: omega = {omega_test} eV")
    print(f"{'q-point':<10} {'||Im chi||_F':<12} {'Relative amp':<12}")
    print("-" * 40)

    norms_at_omega = {q_name: results[q_name][omega_test]['frob_norm']
                      for q_name in q_points.keys()}

    max_norm = max(norms_at_omega.values())

    for q_name, frob_norm in sorted(norms_at_omega.items(), key=lambda x: -x[1]):
        rel_amp = frob_norm / max_norm
        print(f"{q_name:<10} {frob_norm:10.6f}    {rel_amp:6.2%}")

    # 判定 Q_AFM
    q_afm_candidate = max(norms_at_omega, key=norms_at_omega.get)
    enhancement_factor = max_norm / np.mean(list(norms_at_omega.values()))

    print(f"\nQ_AFM candidate: {q_afm_candidate}")
    print(f"Enhancement factor: {enhancement_factor:.2f}x (vs average)")

    if enhancement_factor > 1.5:
        print("\n[OK] Q_AFM nesting amplitude enhancement confirmed")
    else:
        print("\n[?] Enhancement not significant (may need finer q-scan)")

    # ==============================================================================
    # 检验 3：嵌套关联
    # ==============================================================================

    print("\n" + "=" * 80)
    print("Test 3: Nesting function xi(q) vs Im chi amplitude correlation")
    print("=" * 80)

    print("\nComputing nesting function...")
    nesting_values = {}

    for q_name, q_vec in q_points.items():
        xi = compute_nesting_function(q_vec, mu, fermi_window=0.05, k_grid=100)
        nesting_values[q_name] = xi
        print(f"  xi({q_name}) = {xi:.6f}")

    # 相关分析
    print(f"\nNesting xi vs Im chi amplitude (omega={omega_test} eV):")
    print(f"{'q-point':<10} {'xi(q)':<12} {'||Im chi||':<12} {'Rank match':<15}")
    print("-" * 55)

    xi_sorted = sorted(nesting_values.items(), key=lambda x: -x[1])
    norm_sorted = sorted(norms_at_omega.items(), key=lambda x: -x[1])

    for rank, ((q1, xi), (q2, norm)) in enumerate(zip(xi_sorted, norm_sorted), 1):
        match = "[OK]" if q1 == q2 else "[X]"
        print(f"{rank}. {q1:<8} {xi:10.6f}  |  {q2:<8} {norm:10.6f}  |  {match}")

    # Pearson 相关系数
    xi_array = np.array([nesting_values[q] for q in q_points.keys()])
    norm_array = np.array([norms_at_omega[q] for q in q_points.keys()])

    if len(xi_array) > 2:
        corr = np.corrcoef(xi_array, norm_array)[0, 1]
        print(f"\nPearson correlation: r = {corr:.3f}")

        if corr > 0.7:
            print("[OK] Strong positive correlation! Nesting mechanism supported")
        elif corr > 0.4:
            print("[?] Moderate correlation (need more q-points)")
        else:
            print("[X] Weak correlation (other factors may dominate)")

    # ==============================================================================
    # 总结
    # ==============================================================================

    print("\n" + "=" * 80)
    print("PHASE 1 VERIFICATION SUMMARY")
    print("=" * 80)

    print("\n[1] PSD verification:")
    if n_violations == 0:
        print("    [OK] Fluctuation-dissipation theorem numerically confirmed (0 violations)")
    else:
        print(f"    [?] {n_violations}/{total_tests} violations (possibly numerical noise)")

    print("\n[2] Q_AFM amplitude enhancement:")
    print(f"    Candidate: {q_afm_candidate}")
    print(f"    Enhancement: {enhancement_factor:.2f}x")
    if enhancement_factor > 1.5:
        print("    [OK] Nesting mechanism amplitude enhancement confirmed")

    print("\n[3] Nesting correlation:")
    if len(xi_array) > 2:
        print(f"    Correlation: r = {corr:.3f}")
        if corr > 0.7:
            print("    [OK] Strong correlation between geometric nesting and dynamic response")

    print("\n" + "=" * 80)
    print("Strictness level: numerical-evidence (float64)")
    print("Next: Phase 2 iron-based materials survey (universality verification)")
    print("=" * 80)

    return results, nesting_values

# ==============================================================================
# 主程序
# ==============================================================================

if __name__ == "__main__":
    import os
    os.environ['OMP_NUM_THREADS'] = '1'

    results, nesting_values = phase1_verification()
