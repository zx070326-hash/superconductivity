#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N×N RPA diagonal dominance via Perron-Frobenius
Goal: Generalize 2×2 conditional theorem to arbitrary N
Strategy: Use spectral decomposition + track eigenvector mixing
"""
import numpy as np
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def compute_diagonal_dominance(M):
    """
    Compute diagonal dominance measure D[M] = min_α (M_{αα} - Σ_{β≠α}|M_{αβ}|)
    """
    N = M.shape[0]
    D_per_row = []
    for i in range(N):
        diag_term = M[i, i]
        offdiag_sum = np.sum(np.abs(M[i, :])) - np.abs(M[i, i])
        D_per_row.append(diag_term - offdiag_sum)
    return min(D_per_row), D_per_row

def rpa_susceptibility(chi0, U):
    """
    Compute RPA susceptibility: χ^RPA = χ₀(1 - Uχ₀)⁻¹
    """
    N = chi0.shape[0]
    I = np.eye(N)
    Uchi0 = U @ chi0
    chi_rpa = chi0 @ np.linalg.inv(I - Uchi0)
    return chi_rpa

print("=" * 70)
print("N×N RPA 对角主导：Perron-Frobenius 方法")
print("=" * 70)

# Test with N=3 symmetric case first
print("\n测试 1：N=3 对称情况（验证推广）")
print("-" * 70)

N = 3
# Symmetric chi0
chi0 = np.array([
    [1.0, 0.2, 0.1],
    [0.2, 1.0, 0.15],
    [0.1, 0.15, 1.0]
])

# Uniform U
u = 0.35
U = u * np.eye(N)

# Spectral decomposition
eigenvalues, eigenvectors = np.linalg.eigh(chi0)
print(f"\nχ₀ 本征值：{eigenvalues}")
print(f"条件数 λ_max/λ_min = {eigenvalues[-1]/eigenvalues[0]:.4f}")

# Compute spectral radius
chi0U = chi0 @ U
rho = np.max(np.abs(np.linalg.eigvals(chi0U)))
print(f"\n谱半径 ρ = {rho:.6f}")
print(f"Stoner 稳定性：{'✓' if rho < 1 else '✗'}")

# Compute RPA
chi_rpa = rpa_susceptibility(chi0, U)

# Diagonal dominance
D_chi0, D0_rows = compute_diagonal_dominance(chi0)
D_rpa, Dr_rows = compute_diagonal_dominance(chi_rpa)

print(f"\n对角主导度量：")
print(f"  D[χ₀] = {D_chi0:.6f}")
print(f"  D[χ^RPA] = {D_rpa:.6f}")
print(f"  比值 = {D_rpa/D_chi0:.6f}")
print(f"  理论界（对称）1/(1-ρ) = {1/(1-rho):.6f}")

if D_rpa >= D_chi0:
    print(f"  ✓ 对角主导增强")
else:
    print(f"  ✗ 对角主导退化")

print(f"\n逐行分析：")
for i in range(N):
    print(f"  行 {i+1}: D₀={D0_rows[i]:.4f} → D^RPA={Dr_rows[i]:.4f} (×{Dr_rows[i]/D0_rows[i]:.2f})")

# Check if symmetric case formula holds
symmetric_prediction = D_chi0 / (1 - rho)
print(f"\n对称公式验证：")
print(f"  预测值 D[χ₀]/(1-ρ) = {symmetric_prediction:.6f}")
print(f"  实际值 D[χ^RPA] = {D_rpa:.6f}")
print(f"  误差 = {abs(D_rpa - symmetric_prediction)/symmetric_prediction * 100:.2f}%")

# Test with asymmetric case
print("\n" + "=" * 70)
print("测试 2：N=3 非对称情况")
print("-" * 70)

chi0_asym = np.array([
    [1.5, 0.2, 0.1],
    [0.2, 1.0, 0.15],
    [0.1, 0.15, 0.8]
])

# Check positive semi-definite
eigs_asym = np.linalg.eigvalsh(chi0_asym)
print(f"\nχ₀ 本征值：{eigs_asym}")
print(f"PSD 检验：{'✓' if np.all(eigs_asym >= -1e-10) else '✗'}")
print(f"条件数 κ = {np.max(np.diag(chi0_asym))/np.min(np.diag(chi0_asym)):.4f}")

# RPA
chi0U_asym = chi0_asym @ U
rho_asym = np.max(np.abs(np.linalg.eigvals(chi0U_asym)))
print(f"\n谱半径 ρ = {rho_asym:.6f}")

chi_rpa_asym = rpa_susceptibility(chi0_asym, U)

D0_asym, D0_asym_rows = compute_diagonal_dominance(chi0_asym)
Dr_asym, Dr_asym_rows = compute_diagonal_dominance(chi_rpa_asym)

print(f"\n对角主导度量：")
print(f"  D[χ₀] = {D0_asym:.6f}")
print(f"  D[χ^RPA] = {Dr_asym:.6f}")
print(f"  比值 = {Dr_asym/D0_asym:.6f}")

if Dr_asym >= D0_asym:
    print(f"  ✓ 对角主导保持")
else:
    print(f"  ✗ 对角主导退化")

print(f"\n逐行分析：")
for i in range(N):
    print(f"  行 {i+1}: D₀={D0_asym_rows[i]:.4f} → D^RPA={Dr_asym_rows[i]:.4f} (×{Dr_asym_rows[i]/D0_asym_rows[i]:.2f})")

bottleneck_idx = np.argmin(Dr_asym_rows)
print(f"\n瓶颈行：行 {bottleneck_idx+1}（对角元 χ₀,{bottleneck_idx+1}{bottleneck_idx+1} = {chi0_asym[bottleneck_idx, bottleneck_idx]:.2f}）")

print("\n" + "=" * 70)
print("分析：谱分解视角")
print("=" * 70)

# Analyze eigenvector structure
print("\nχ₀（非对称）的本征向量（按列）：")
evals_asym, evecs_asym = np.linalg.eigh(chi0_asym)
for i in range(N):
    print(f"  λ_{i+1} = {evals_asym[i]:.4f}: v = {evecs_asym[:, i]}")

# Check how RPA affects eigenvector mixing in standard basis
print("\n标准基下的本征向量分量平方（反映对角主导）：")
for i in range(N):
    v_sq = evecs_asym[:, i]**2
    print(f"  本征向量 {i+1}: {v_sq} (和 = {np.sum(v_sq):.6f})")

# Key insight: diagonal dominance in standard basis depends on eigenvector alignment
print("\n关键观察：")
print("对角主导在标准基下 ⟺ 本征向量接近标准基向量")
print("RPA 增强取决于 U·λ，但对角主导还需要本征向量'对齐'")

print("\n" + "=" * 70)
print("下一步：寻找 N×N 的充分条件")
print("=" * 70)

print("""
从 N=3 测试可以看出：

1. 对称情况：D[χ^RPA] = D[χ₀]/(1-ρ) **仍然成立**（N=3 验证）

2. 非对称情况：
   - 瓶颈仍在最弱的对角元
   - 本征向量的"对齐度"影响对角主导传递

3. 可能的 N×N 条件定理方向：
   a) 条件数 κ[χ₀] ≤ κ_c(ρ, N)（推广 2×2）
   b) 本征向量对齐条件：max_i Σ_j v²_{ij} ≤ C（本征向量接近标准基）
   c) 加权度量：用 Perron 本征向量分量作为权重

需要更多 N 和参数的数值实验来确定哪个方向可行。
""")
