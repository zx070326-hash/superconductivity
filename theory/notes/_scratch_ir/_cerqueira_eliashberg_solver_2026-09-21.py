#!/usr/bin/env python3
"""
Cerqueira 6912: 自解 Eliashberg + 定理 2 检验

读取 a2F.dos* 文件（频率单位 Rydberg），用各向同性线性化 Eliashberg 求解器
（μ*=0）解 Tc，计算 𝔣 = Tc/(C_∞√λ·ω̄₂) 与 𝔣/𝔣_E(λ)，检验定理 2。

单位约定：
- 输入：Rydberg（a2F.dos* 表头）
- 内部：meV（k_B·K 单位，Tc 直接是开尔文）
- 输出：开尔文、meV

Rydberg 转换：1 Ry = 13605.693122994 meV = 1178.2662 K
"""
import os, sys, re
import numpy as np
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d

# 单位转换常数
RY_TO_MEV = 13605.693122994  # 1 Ry in meV
MEV_TO_K = 11.604518121550082  # k_B in K/meV
RY_TO_K = RY_TO_MEV * MEV_TO_K  # 1 Ry in K = 157887.51 K

# 天花板常数（仓库）
C_INF = 0.1827262477890
CARBOTTE_C0 = 0.230868  # Carbotte-Leavens 面积界常数

def load_a2f_dos(filepath):
    """
    读取 a2F.dos* 文件，返回 (omega_ry, a2f_total, lambda_file, delta_file)
    频率单位 Rydberg，第 2 列是 a2F_total
    """
    data = []
    lambda_file = None
    delta_file = None

    with open(filepath, 'r') as f:
        for line in f:
            if line.strip().startswith('#'):
                continue
            if 'lambda' in line.lower():
                # 末行：lambda = ... Delta = ...
                parts = line.split()
                for i, p in enumerate(parts):
                    if p == 'lambda' and i+2 < len(parts):
                        lambda_file = float(parts[i+2])
                    if p == 'Delta' and i+2 < len(parts):
                        delta_file = float(parts[i+2])
                continue

            vals = line.split()
            if len(vals) >= 2:
                try:
                    omega = float(vals[0])
                    a2f = float(vals[1])
                    data.append((omega, a2f))
                except ValueError:
                    continue

    if not data:
        return None, None, None, None

    omega_ry = np.array([d[0] for d in data])
    a2f_total = np.array([d[1] for d in data])

    return omega_ry, a2f_total, lambda_file, delta_file

def compute_moments(omega_ry, a2f_total):
    """
    计算矩：λ = 2∫α²F/ω dω, ⟨ω²⟩ = (2/λ)∫α²F·ω dω, A = ∫α²F dω

    全部在 Ry 单位下积分（与用户复核脚本一致）

    返回：
    - lam: 无量纲
    - omega2_mev: meV
    - area_ry: Ry（积分值，用于 Carbotte 界时需乘 RY_K）
    - M1_ry: Ry（一阶矩，用于矩界）
    """
    # 只取正频率
    mask = omega_ry > 0
    omega_ry_pos = omega_ry[mask]
    a2f_pos = a2f_total[mask]

    if len(omega_ry_pos) == 0:
        return None, None, None, None

    # 在 Ry 单位下积分
    lam = 2.0 * trapezoid(a2f_pos / omega_ry_pos, omega_ry_pos)

    if lam <= 0:
        return None, None, None, None

    M1_ry = trapezoid(a2f_pos * omega_ry_pos, omega_ry_pos)
    area_ry = trapezoid(a2f_pos, omega_ry_pos)

    omega2_avg_ry2 = 2.0 * M1_ry / lam
    omega2_mev = np.sqrt(omega2_avg_ry2) * RY_TO_MEV

    return lam, omega2_mev, area_ry, M1_ry

def solve_eliashberg_isotropic(omega_mev, a2f_total, T_K, N=400):
    """
    各向同性线性化 Eliashberg（μ*=0），AD75 完整对称矩阵形式

    返回对称矩阵 M 的最大本征值（Tc 时 = 0）

    参考：用户验证脚本 _chk_verify_t5_area_vs_moment_batchq_2026-09-21.py
    """
    # 只取正频率
    mask = omega_mev > 0
    omega_mev_pos = omega_mev[mask]
    a2f_pos = a2f_total[mask]

    # 转回 Ry 单位（用户的代码在 Ry 下工作）
    omega_ry = omega_mev_pos / RY_TO_MEV
    T_ry = T_K / RY_TO_K

    # 计算所有需要的 lam_l，l = 0, 1, ..., 2N+1
    nu = 2.0 * np.pi * T_ry * np.arange(0, 2 * N + 2)  # Ry
    integ = 2.0 * a2f_pos[None, :] * omega_ry[None, :] / (omega_ry[None, :] ** 2 + nu[:, None] ** 2)
    lam_l = np.trapezoid(integ, omega_ry, axis=1)

    # 构造 AD75 对称矩阵
    n = np.arange(N)
    M = lam_l[np.abs(n[:, None] - n[None, :])] + lam_l[n[:, None] + n[None, :] + 1]

    # 对角元修正
    cum = np.concatenate(([0.0], np.cumsum(lam_l[1:N])))
    M[np.diag_indices(N)] -= (2 * n + 1) + lam_l[0] + 2 * cum

    # 返回最大本征值
    return np.linalg.eigvalsh(M)[-1]

def find_tc_bisection(omega_mev, a2f_total, T_low=None, T_high=None, N=400, tol=1e-6):
    """
    二分法找 Tc（M 的最大本征值=0）
    温度单位：开尔文
    自动设置范围：T_low = 1e-12 K, T_high = 0.3 * omega_max (converted to K)
    """
    if T_low is None:
        T_low = 1e-12  # K，极低温
    if T_high is None:
        # 用户用 0.3 * omega_max(Ry) 作为温度(Ry)，转成 K：
        omega_max_ry = np.max(omega_mev) / RY_TO_MEV
        T_high = 0.3 * omega_max_ry * RY_TO_K  # K

    lam_low = solve_eliashberg_isotropic(omega_mev, a2f_total, T_low, N)
    lam_high = solve_eliashberg_isotropic(omega_mev, a2f_total, T_high, N)

    # 二分法要求：lam_low > 0（T低时本征值正）且 lam_high < 0（T高时本征值负）
    # 如果 lam_high > 0，说明 T_high 不够高
    if lam_high > 0.0:
        # T_high 不够高，翻倍重试
        T_high *= 2.0
        lam_high = solve_eliashberg_isotropic(omega_mev, a2f_total, T_high, N)
        if lam_high > 0.0:
            return np.nan  # 仍然不够

    for _ in range(60):  # 用户用60次迭代
        T_mid = 0.5 * (T_low + T_high)
        lam_mid = solve_eliashberg_isotropic(omega_mev, a2f_total, T_mid, N)

        if lam_mid > 0.0:
            T_low = T_mid
        else:
            T_high = T_mid

    return 0.5 * (T_low + T_high)

def solve_einstein_tc(lam, omega2_mev, N=400, tol=1e-6):
    """
    Einstein 谱的 Tc（用于 𝔣_E(λ)）
    λ(iω_n) = λ·ω₀²/(ω₀²+ω_n²)
    """
    omega0_mev = omega2_mev

    def lambda_max_einstein(T_K):
        omega_n = 2 * np.pi * T_K / MEV_TO_K * (np.arange(N) + 0.5)
        lam_n = lam * omega0_mev**2 / (omega0_mev**2 + omega_n**2)

        # 对角核矩阵
        K_diag = (np.pi * MEV_TO_K / T_K) * lam_n / (2 * omega_n)

        return np.max(K_diag)

    # 二分法，要求 lam_max(Tc) = 1
    T_low, T_high = 1e-12, 50.0
    lam_low = lambda_max_einstein(T_low)
    lam_high = lambda_max_einstein(T_high)

    # 如果 T_high 处仍 > 1，扩大上界
    if lam_high > 1.0:
        T_high = 200.0
        lam_high = lambda_max_einstein(T_high)
        if lam_high > 1.0:
            return np.nan

    while T_high - T_low > tol:
        T_mid = 0.5 * (T_low + T_high)
        lam_mid = lambda_max_einstein(T_mid)

        if lam_mid > 1.0:
            T_low = T_mid
        else:
            T_high = T_mid

    return 0.5 * (T_low + T_high)

def process_material(material_dir, smearing_idx=5, N=400):
    """
    处理单个材料：读 a2F.dos<smearing_idx>，解 Tc，计算 𝔣

    返回 dict 或 None
    """
    a2f_file = os.path.join(material_dir, f'a2F.dos{smearing_idx}')
    if not os.path.exists(a2f_file):
        return None

    omega_ry, a2f_total, lambda_file, delta_file = load_a2f_dos(a2f_file)
    if omega_ry is None:
        return None

    lam, omega2_mev, area_ry, M1_ry = compute_moments(omega_ry, a2f_total)
    if lam is None or lam <= 0:
        return None

    omega_mev = omega_ry * RY_TO_MEV

    # 解真实谱 Tc
    Tc_real = find_tc_bisection(omega_mev, a2f_total, N=N)
    if np.isnan(Tc_real):
        return None  # 数值失败

    # 解 Einstein Tc
    Tc_einstein = solve_einstein_tc(lam, omega2_mev, N=N)
    if np.isnan(Tc_einstein):
        return None

    # 计算天花板和 𝔣
    # 正确公式：k_B Tc ≤ C_∞ · √λ · ω̄₂  （ω̄₂ 不开方）
    Tc_ceiling = C_INF * np.sqrt(lam) * omega2_mev * MEV_TO_K
    f_deficit = Tc_real / Tc_ceiling
    f_ratio = Tc_real / Tc_einstein  # 𝔣/𝔣_E

    # Carbotte-Leavens 面积界
    # 正确公式：k_B Tc ≤ c(0) · A，其中 A 在 Ry 单位
    Tc_carbotte = CARBOTTE_C0 * area_ry * RY_TO_K
    f_carbotte = Tc_real / Tc_carbotte

    # 矩界（独立验证）
    # 正确公式：k_B Tc ≤ C_∞ · √(2M₁)
    Tc_moment = C_INF * np.sqrt(2 * M1_ry) * RY_TO_K

    # 界比较：直接比大小，不要用填充因子
    carbotte_tighter = (Tc_carbotte < Tc_ceiling)
    exceeds_carbotte = (Tc_real > Tc_carbotte)
    exceeds_ceiling = (Tc_real > Tc_ceiling)
    exceeds_moment = (Tc_real > Tc_moment)

    material_name = os.path.basename(material_dir)

    return {
        'material': material_name,
        'smearing': smearing_idx,
        'lambda': lam,
        'omega2_mev': omega2_mev,
        'area_ry': area_ry,
        'M1_ry': M1_ry,
        'lambda_file': lambda_file,
        'Tc_real_K': Tc_real,
        'Tc_einstein_K': Tc_einstein,
        'Tc_ceiling_K': Tc_ceiling,
        'Tc_carbotte_K': Tc_carbotte,
        'f_deficit': f_deficit,
        'f_ratio': f_ratio,
        'f_carbotte': f_carbotte,
        'carbotte_tighter': carbotte_tighter,
        'exceeds_carbotte': exceeds_carbotte,
        'exceeds_ceiling': exceeds_ceiling,
        'exceeds_moment': exceeds_moment,
        'Tc_moment_K': Tc_moment,
        'N': N
    }

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print("Usage: python _cerqueira_eliashberg_solver_2026-09-21.py <material_dir>")
        sys.exit(1)

    material_dir = sys.argv[1]

    try:
        result = process_material(material_dir, smearing_idx=5, N=400)
    except Exception as e:
        print(f"ERROR: {material_dir}")
        print(f"Exception: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    if result is None:
        print(f"FAIL: {material_dir} (process_material returned None)")
        sys.exit(1)

    print(f"Material: {result['material']}")
    print(f"  lambda = {result['lambda']:.6f} (file: {result['lambda_file']:.6f})")
    print(f"  omega2_mev = {result['omega2_mev']:.3f} meV")
    print(f"  area = {result['area']:.6f}")
    print(f"  Tc_real = {result['Tc_real_K']:.4f} K")
    print(f"  Tc_einstein = {result['Tc_einstein_K']:.4f} K")
    print(f"  Tc_ceiling = {result['Tc_ceiling_K']:.4f} K")
    print(f"  Tc_carbotte = {result['Tc_carbotte_K']:.4f} K")
    print(f"  f_deficit = {result['f_deficit']:.6f}")
    print(f"  f_ratio = {result['f_ratio']:.6f} {'VIOLATION' if result['f_ratio'] > 1 else 'PASS'}")
    print(f"  f_carbotte = {result['f_carbotte']:.6f}")
