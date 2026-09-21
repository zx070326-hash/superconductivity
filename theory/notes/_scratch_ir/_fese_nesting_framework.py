#!/usr/bin/env python3
"""
FeSe 费米面嵌套函数 + Stoner 稳定性计算框架
Team C (fese-nesting-mechanism) Phase 1 准备

目标：
1. 计算嵌套函数 ξ_{αβ}(q) = ∫ dk δ(ε_α(k)) δ(ε_β(k+q))
2. 计算 Lindhard χ₀ 矩阵
3. 计算 Stoner 谱半径 ρ(q) = λ_max(χ₀(q)·U)
4. 关联：ρ(q) vs ξ_{αβ}(q) vs χ^s_{ij} 符号

更新于 2026-09-22：整合 Team A Stoner 稳定性机制

参考：
- Essenberger PRB 94, 014503 (2016) Fig. 4c
- 研究计划 theory/notes/meta/rpa_offdiag_positivity_research_plan_2026-09-22.md
- Team A 洞察 theory/notes/meta/rpa_stoner_boundary_insight_2026-09-22.md
- 理论框架 theory/notes/meta/fese_chi_sign_stoner_mechanism_2026-09-22.md
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from typing import Tuple, Callable

# ============================================================================
# 1. 简化紧束缚模型 (待 Team B 给出实际参数)
# ============================================================================

class SimpleFeSeFermiSurface:
    """
    简化 FeSe 费米面模型

    基于 Essenberger 1411.2121 Fig. 1:
    - Band +: Γ 点空穴袋 (hole pocket at Γ)
    - Band −: M 点电子袋 (electron pocket at M)

    这里用简化的圆形/椭圆形费米面作为占位
    实际参数需要 Team B 从文献中提取
    """

    def __init__(self):
        # 占位参数 (需要 Team B 提供实际值)
        self.k_F_hole = 0.2  # Γ 点空穴费米波矢 (2π/a 单位)
        self.k_F_elec = 0.3  # M 点电子费米波矢
        self.M_point = (np.pi, 0.0)  # M 点位置 (约化布里渊区)

    def epsilon_hole(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        """Γ 点空穴带色散 (占位)"""
        return -1.0 + (kx**2 + ky**2) / (2 * self.k_F_hole**2)

    def epsilon_elec(self, kx: np.ndarray, ky: np.ndarray) -> np.ndarray:
        """M 点电子带色散 (占位)"""
        qx = kx - self.M_point[0]
        qy = ky - self.M_point[1]
        return 1.0 - (qx**2 + qy**2) / (2 * self.k_F_elec**2)



# ============================================================================
# 2. Lindhard 函数计算 (新增：Phase 1 扩展)
# ============================================================================

def compute_lindhard_matrix(
    q: Tuple[float, float],
    fermi_surface: SimpleFeSeFermiSurface,
    omega: float = 0.0,
    T: float = 0.01,
    n_k: int = 200
) -> np.ndarray:
    """
    计算 Lindhard 极化函数矩阵 χ₀(q,ω)

    χ₀_{αβ}(q,0) = Σ_k [f(ε_α(k)) - f(ε_β(k+q))] / [ε_β(k+q) - ε_α(k)]

    返回 2×2 矩阵（hole-hole, hole-elec; elec-hole, elec-elec）

    参数：
        q: 波矢 (qx, qy)
        fermi_surface: 费米面模型
        omega: 频率（静态 ω=0）
        T: 温度（T→0 极限，用于 Fermi 函数平滑）
        n_k: k 空间网格密度

    返回：
        chi0: 2×2 Lindhard 矩阵
    """
    qx, qy = q
    kx = np.linspace(-np.pi, np.pi, n_k)
    ky = np.linspace(-np.pi, np.pi, n_k)
    KX, KY = np.meshgrid(kx, ky)

    # 计算两个带的色散
    E_hole = fermi_surface.epsilon_hole(KX, KY)
    E_elec = fermi_surface.epsilon_elec(KX, KY)

    # k+q 处的色散（周期边界条件）
    KXq = np.mod(KX + qx + np.pi, 2*np.pi) - np.pi
    KYq = np.mod(KY + qy + np.pi, 2*np.pi) - np.pi
    E_hole_q = fermi_surface.epsilon_hole(KXq, KYq)
    E_elec_q = fermi_surface.epsilon_elec(KXq, KYq)

    # Fermi 分布（T→0 极限用阶跃函数近似）
    def fermi(E):
        if T < 1e-6:
            return (E < 0).astype(float)
        else:
            return 1.0 / (1.0 + np.exp(E / T))

    f_hole = fermi(E_hole)
    f_elec = fermi(E_elec)
    f_hole_q = fermi(E_hole_q)
    f_elec_q = fermi(E_elec_q)

    # 计算 χ₀ 矩阵元
    chi0 = np.zeros((2, 2))
    dk = (2*np.pi / n_k)**2  # k 空间积分元
    eta = 0.01  # 避免分母为零

    # (0,0): hole-hole
    denom_hh = E_hole_q - E_hole + omega + 1j*eta
    chi0[0, 0] = np.sum((f_hole - f_hole_q) / denom_hh).real * dk

    # (0,1): hole-elec
    denom_he = E_elec_q - E_hole + omega + 1j*eta
    chi0[0, 1] = np.sum((f_hole - f_elec_q) / denom_he).real * dk

    # (1,0): elec-hole
    denom_eh = E_hole_q - E_elec + omega + 1j*eta
    chi0[1, 0] = np.sum((f_elec - f_hole_q) / denom_eh).real * dk

    # (1,1): elec-elec
    denom_ee = E_elec_q - E_elec + omega + 1j*eta
    chi0[1, 1] = np.sum((f_elec - f_elec_q) / denom_ee).real * dk

    return chi0


def compute_spectral_radius(chi0: np.ndarray, U: np.ndarray) -> float:
    """
    计算 Stoner 谱半径 ρ(χ₀U)

    ρ := max |λ_i(χ₀·U)|

    Stoner 稳定性条件：ρ < 1
    """
    M = np.dot(chi0, U)
    eigvals = np.linalg.eigvals(M)
    return np.max(np.abs(eigvals))


def compute_correlation_coefficient(chi: np.ndarray) -> float:
    """
    计算相关系数（对角主导度量的逆）

    r := |χ₁₂| / √(χ₁₁·χ₂₂)

    Team A 定理：ρ < 1 ⟹ r_n → 0（几何衰减）
    """
    if chi.shape != (2, 2):
        raise ValueError("Only 2×2 matrices supported")

    chi11, chi22 = chi[0, 0], chi[1, 1]
    chi12 = chi[0, 1]

    if chi11 <= 0 or chi22 <= 0:
        return np.nan  # 对角元非正，定义失效

    return np.abs(chi12) / np.sqrt(chi11 * chi22)


# ============================================================================
# 3. 嵌套函数计算（原有，保持不变）
# ============================================================================

def compute_nesting_function(
    q: Tuple[float, float],
    fermi_surface: SimpleFeSeFermiSurface,
    band_alpha: str = 'hole',
    band_beta: str = 'elec',
    n_k: int = 200,
    delta_E: float = 0.01
) -> float:
    """
    计算嵌套函数 ξ_{αβ}(q)

    ξ_{αβ}(q) = ∫ dk δ(ε_α(k)) δ(ε_β(k+q))

    物理意义：衡量波矢 q 连接轨道 α 和 β 费米面的程度

    参数：
        q: 波矢 (qx, qy)
        fermi_surface: 费米面模型
        band_alpha, band_beta: 'hole' 或 'elec'
        n_k: k 空间网格密度
        delta_E: 费米面窗口宽度

    返回：
        nesting: 嵌套度 (归一化)
    """
    qx, qy = q

    # k 空间网格 (约化布里渊区 [-π, π])
    kx = np.linspace(-np.pi, np.pi, n_k)
    ky = np.linspace(-np.pi, np.pi, n_k)
    KX, KY = np.meshgrid(kx, ky)

    # 计算 ε_α(k)
    if band_alpha == 'hole':
        E_alpha = fermi_surface.epsilon_hole(KX, KY)
    else:
        E_alpha = fermi_surface.epsilon_elec(KX, KY)

    # 计算 ε_β(k+q) (周期边界条件)
    KXq = np.mod(KX + qx + np.pi, 2*np.pi) - np.pi
    KYq = np.mod(KY + qy + np.pi, 2*np.pi) - np.pi

    if band_beta == 'hole':
        E_beta = fermi_surface.epsilon_hole(KXq, KYq)
    else:
        E_beta = fermi_surface.epsilon_elec(KXq, KYq)

    # δ 函数近似为窗口函数
    on_FS_alpha = np.abs(E_alpha) < delta_E
    on_FS_beta = np.abs(E_beta) < delta_E

    # 嵌套度 = 同时在两个费米面上的 k 点数量
    nesting = np.sum(on_FS_alpha & on_FS_beta)

    # 归一化 (除以费米面 α 的总点数)
    norm = np.sum(on_FS_alpha) + 1e-10

    return nesting / norm



# ============================================================================
# 4. Stoner 稳定性扫描（新增：核心任务）
# ============================================================================

def scan_stoner_stability(
    fermi_surface: SimpleFeSeFermiSurface,
    U: np.ndarray,
    n_q: int = 50
) -> dict:
    """
    扫描布里渊区，绘制 ρ(q) 热图

    这是 Phase 2 的核心任务：验证 Essenberger FeSe 满足 ρ < 1

    返回：
        results: {
            'qx': q_x 网格,
            'qy': q_y 网格,
            'rho': ρ(q) 热图,
            'chi0_norm': ‖χ₀(q)‖ 热图,
            'rho_max': max_q ρ(q),
            'q_max': argmax_q ρ(q)
        }
    """
    qx = np.linspace(-np.pi, np.pi, n_q)
    qy = np.linspace(0, np.pi, n_q)

    rho_map = np.zeros((n_q, n_q))
    chi0_norm_map = np.zeros((n_q, n_q))

    print(f"扫描 Stoner 稳定性...")
    print(f"网格：{n_q}×{n_q} = {n_q**2} 点")

    rho_max = 0.0
    q_max = (0.0, 0.0)

    for i, qxi in enumerate(qx):
        for j, qyj in enumerate(qy):
            chi0 = compute_lindhard_matrix((qxi, qyj), fermi_surface)
            rho = compute_spectral_radius(chi0, U)
            chi0_norm = np.linalg.norm(chi0)

            rho_map[j, i] = rho
            chi0_norm_map[j, i] = chi0_norm

            if rho > rho_max:
                rho_max = rho
                q_max = (qxi, qyj)

        if (i+1) % 10 == 0:
            print(f"  进度：{i+1}/{n_q}")

    results = {
        'qx': qx,
        'qy': qy,
        'rho': rho_map,
        'chi0_norm': chi0_norm_map,
        'rho_max': rho_max,
        'q_max': q_max
    }

    print(f"完成。max_q ρ(q) = {rho_max:.4f} at q = {q_max}")

    if rho_max >= 1.0:
        print("⚠ WARNING: Stoner 不稳定！ρ ≥ 1")
    else:
        print("✓ Stoner 稳定：ρ < 1")

    return results


def check_special_points_stoner(
    fermi_surface: SimpleFeSeFermiSurface,
    U: np.ndarray
):
    """
    检验特殊波矢的 Stoner 稳定性 + 嵌套度

    对照 Team B 的符号检验点
    """
    test_points = {
        'Q_AFM': (np.pi, 0.0),
        'Gamma': (0.0, 0.0),
        'M': (np.pi, np.pi),
        'X_1': (np.pi, 0.0),
        'X_2': (0.0, np.pi),
    }

    print("\n=== 特殊点：Stoner 稳定性 + 嵌套度 ===")
    print(f"{'Point':<10} {'ρ(q)':<10} {'‖χ₀‖':<10} {'r':<10} {'ξ_he':<10}")
    print("-" * 60)

    for name, q in test_points.items():
        # Stoner 量
        chi0 = compute_lindhard_matrix(q, fermi_surface)
        rho = compute_spectral_radius(chi0, U)
        chi0_norm = np.linalg.norm(chi0)
        r = compute_correlation_coefficient(chi0)

        # 嵌套度
        xi_he = compute_nesting_function(q, fermi_surface, 'hole', 'elec')

        status = "✓" if rho < 1.0 else "✗"
        print(f"{name:<10} {rho:<10.4f} {chi0_norm:<10.4f} {r:<10.4f} {xi_he:<10.4f} {status}")


# ============================================================================
# 5. 布里渊区扫描（更新：增加嵌套度扫描）
# ============================================================================

def scan_brillouin_zone(
    fermi_surface: SimpleFeSeFermiSurface,
    band_pairs: list = [('hole', 'elec'), ('hole', 'hole'), ('elec', 'elec')],
    n_q: int = 50
) -> dict:
    """
    扫描布里渊区，计算所有轨道对的嵌套函数

    这是 Team B 符号图的物理对应
    """
    qx = np.linspace(-np.pi, np.pi, n_q)
    qy = np.linspace(0, np.pi, n_q)  # 利用对称性只扫半个 BZ

    results = {}

    for alpha, beta in band_pairs:
        label = f"{alpha}-{beta}"
        nesting_map = np.zeros((n_q, n_q))

        print(f"Computing nesting for {label}...")

        for i, qxi in enumerate(qx):
            for j, qyj in enumerate(qy):
                nesting_map[j, i] = compute_nesting_function(
                    (qxi, qyj), fermi_surface, alpha, beta
                )

        results[label] = {
            'qx': qx,
            'qy': qy,
            'nesting': nesting_map
        }

    return results


# ============================================================================
# 6. 特殊点检验（更新：整合 Stoner）
# ============================================================================

def check_special_points(fermi_surface: SimpleFeSeFermiSurface):
    """
    检验特殊波矢的嵌套度

    关键预期：
    - Q_AFM ≈ (π, 0): 连接 Γ 空穴和 M 电子 → 高嵌套度
    - q = 0: 带内散射 → 中等嵌套度
    - 随机 q: 低嵌套度
    """
    test_points = {
        'Q_AFM': (np.pi, 0.0),
        'Gamma': (0.0, 0.0),
        'M': (np.pi, np.pi),
        'X': (np.pi/2, 0.0),
    }

    print("\n=== 特殊点嵌套度检验 ===")
    print(f"{'Point':<10} {'hole-elec':<12} {'hole-hole':<12} {'elec-elec':<12}")
    print("-" * 50)

    for name, q in test_points.items():
        xi_he = compute_nesting_function(q, fermi_surface, 'hole', 'elec')
        xi_hh = compute_nesting_function(q, fermi_surface, 'hole', 'hole')
        xi_ee = compute_nesting_function(q, fermi_surface, 'elec', 'elec')

        print(f"{name:<10} {xi_he:<12.4f} {xi_hh:<12.4f} {xi_ee:<12.4f}")


# ============================================================================
# 7. 主函数（更新：Phase 1 → Phase 2 准备）
# ============================================================================

if __name__ == "__main__":
    print("FeSe 费米面嵌套 + Stoner 稳定性分析框架")
    print("=" * 70)
    print("状态：Phase 1 准备完成")
    print("更新：整合 Team A Stoner 稳定性机制（2026-09-22）")
    print("依赖：等待 Team B 的 χ^s_{ij}(q,ω) 符号结构 + 紧束缚参数")
    print("=" * 70)

    # 初始化简化模型
    fs = SimpleFeSeFermiSurface()

    # 占位 Hubbard U（需要 Team B 提供实际值）
    U_placeholder = np.array([[2.0, 0.0],
                              [0.0, 2.0]])  # 对角近似

    print("\n使用占位参数（玩具模型）:")
    print(f"  k_F_hole = {fs.k_F_hole}")
    print(f"  k_F_elec = {fs.k_F_elec}")
    print(f"  U = diag({U_placeholder[0,0]}, {U_placeholder[1,1]})")
    print("\n⚠ 这些参数需要 Team B 从文献中提取实际值")

    # Phase 1 检验：特殊点
    check_special_points(fs)

    # Phase 1 新增：Stoner 稳定性检验
    check_special_points_stoner(fs, U_placeholder)

    print("\n" + "=" * 70)
    print("框架准备完成。收到 Team B 数值后，将执行：")
    print()
    print("Phase 2 任务（按优先级）：")
    print()
    print("1. Stoner 稳定性检验（最关键）")
    print("   - scan_stoner_stability(fs, U_real, n_q=50)")
    print("   - 验证 max_q ρ(q) < 1")
    print("   - 绘制 ρ(q) 热图")
    print()
    print("2. 嵌套度 vs ρ 关联")
    print("   - 对 5 个特殊点：计算 (ξ, ρ, ‖χ₀‖, r)")
    print("   - 检验：ξ 大 → ‖χ₀‖ 大 → ρ 接近 1")
    print("   - 但 ρ 仍 < 1（物理约束）")
    print()
    print("3. 符号结构检验")
    print("   - 与 Team B 符号图对照")
    print("   - 若出现负非对角：分析对应 q 点的 ρ 和 ξ")
    print()
    print("预期交付：2-3 天完成 Phase 2 + 物理机制报告")
    print("=" * 70)

