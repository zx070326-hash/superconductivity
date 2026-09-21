#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
N×N symmetric case: why 2×2 formula fails
Goal: Find correct formula for general N
"""
import numpy as np
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def diagonal_dominance(M):
    N = M.shape[0]
    D_rows = [M[i,i] - (np.sum(np.abs(M[i,:])) - abs(M[i,i])) for i in range(N)]
    return min(D_rows), D_rows

def rpa_chi(chi0, u):
    N = chi0.shape[0]
    I = np.eye(N)
    U = u * np.eye(N)
    return chi0 @ np.linalg.inv(I - U @ chi0)

print("=" * 70)
print("为什么 2×2 公式在 N=3 失败？")
print("=" * 70)

# Test symmetric case with uniform off-diagonal
print("\n测试：N=3，对角元全相同，非对角元全相同")
print("-" * 70)

c = 1.0
c_off = 0.2
u = 0.35

chi0 = c * np.eye(3) + c_off * (np.ones((3, 3)) - np.eye(3))
print(f"\nχ₀ =")
print(chi0)

# Spectral analysis
evals = np.linalg.eigvalsh(chi0)
print(f"\n本征值：{evals}")

# RPA
chi0U = chi0 * u  # element-wise since U = u*I
rho = np.max(np.abs(np.linalg.eigvals(chi0U)))
print(f"\n谱半径 ρ = {rho:.6f}")

chi_rpa = rpa_chi(chi0, u)

D0, D0_rows = diagonal_dominance(chi0)
Dr, Dr_rows = diagonal_dominance(chi_rpa)

print(f"\n对角主导：")
print(f"  D[χ₀] = {D0:.6f}")
print(f"  D[χ^RPA] = {Dr:.6f}")
print(f"  比值 = {Dr/D0:.6f}")

# 2×2 formula prediction
formula_2x2 = D0 / (1 - rho)
print(f"\n2×2 公式预测：D[χ₀]/(1-ρ) = {formula_2x2:.6f}")
print(f"实际值：{Dr:.6f}")
print(f"误差：{abs(Dr - formula_2x2)/Dr * 100:.2f}%")

# Key insight: what is ρ really?
print("\n" + "=" * 70)
print("关键问题：2×2 的 ρ 定义")
print("=" * 70)

print(f"""
2×2 情况：
  χ₀ = [[c, c₁₂], [c₁₂, c]]
  U = diag(u, u)

对称性 ⟹ χ₀U 的本征值 = u·λᵢ(χ₀)
χ₀ 的本征值 = c ± |c₁₂|

所以 ρ(χ₀U) = u·max(c+|c₁₂|, c-|c₁₂|) = u(c + |c₁₂|)

而 D[χ₀] = c - |c₁₂|

所以：
  D[χ₀] / (1 - ρ) = (c - |c₁₂|) / (1 - u(c + |c₁₂|))

这在 2×2 恰好等于 D[χ^RPA]！

但在 N=3：
  χ₀ = c·I + c_off·(J - I)，其中 J 是全 1 矩阵

本征值：
  - λ₁ = c - c_off（重数 N-1）
  - λ₂ = c + (N-1)c_off（重数 1）

所以 ρ = u·(c + (N-1)c_off)

而 D[χ₀] = c - (N-1)c_off

公式给出：
  D[χ₀]/(1-ρ) = (c - (N-1)c_off) / (1 - u(c + (N-1)c_off))

但这**不等于** D[χ^RPA]！
""")

print("\n实际计算验证：")
c_minus = c - (3-1)*c_off
c_plus = c + (3-1)*c_off
rho_formula = u * c_plus
D0_formula = c_minus

print(f"  c - (N-1)c_off = {c_minus:.6f}")
print(f"  c + (N-1)c_off = {c_plus:.6f}")
print(f"  ρ（公式）= u(c + (N-1)c_off) = {rho_formula:.6f}")
print(f"  ρ（数值）= {rho:.6f}")
print(f"  匹配：{'✓' if abs(rho_formula - rho) < 1e-6 else '✗'}")

prediction_wrong = D0_formula / (1 - rho_formula)
print(f"\n  公式预测：{prediction_wrong:.6f}")
print(f"  实际 D[χ^RPA]：{Dr:.6f}")
print(f"  误差：{abs(Dr - prediction_wrong)/Dr * 100:.2f}%")

print("\n" + "=" * 70)
print("为什么公式失败？")
print("=" * 70)

print("""
关键发现：2×2 公式依赖特殊结构！

在 2×2：
  非对角元只有一个（对称矩阵上三角 1 个自由度）
  D[χ₀] = c - |c₁₂|
  ρ = u(c + |c₁₂|)

  两者恰好"对偶"：一个是 c - x，一个是 u(c + x)

在 N≥3：
  非对角元有 N-1 个（每行）
  D[χ₀] = c - (N-1)·c_off（假设全相同）
  ρ = u·λ_max = u(c + (N-1)·c_off)

  但 RPA 的增强**不只看 ρ**，还要看矩阵结构！

具体：χ^RPA 的非对角元不是简单地按 1/(1-ρ) 缩放。
""")

# Let's compute the actual RPA off-diagonal elements
print("\n实际 RPA 矩阵元：")
print(f"\nχ^RPA =")
print(chi_rpa)

print(f"\nχ^RPA 的对角元：{np.diag(chi_rpa)}")
print(f"χ^RPA 的非对角元（第 1 行）：{chi_rpa[0, 1:]}")

chi_rpa_offdiag = chi_rpa[0, 1]
chi0_offdiag = chi0[0, 1]
offdiag_ratio = chi_rpa_offdiag / chi0_offdiag

print(f"\n非对角元增强倍数：χ^RPA_{12} / χ₀_{12} = {offdiag_ratio:.4f}")
print(f"对角元增强倍数：χ^RPA_{11} / χ₀_{11} = {chi_rpa[0,0] / chi0[0,0]:.4f}")

print(f"\n关键：非对角元和对角元的增强倍数**不同**！")
print(f"这导致 D[χ^RPA]/D[χ₀] ≠ 简单的 1/(1-ρ) 关系")

print("\n" + "=" * 70)
print("正确的 N×N 对称公式（猜想）")
print("=" * 70)

print("""
需要找到：
  D[χ^RPA] 与 D[χ₀]、ρ、N 的真实关系

可能形式：
  D[χ^RPA] = D[χ₀] · g(ρ, N)

其中 g(ρ, N) 不是简单的 1/(1-ρ)

从 N=2 到 N=3 的数据：
  N=2: g = 1/(1-ρ) ≈ 2.08（ρ=0.52 时）
  N=3: g ≈ 1.13（ρ=0.46 时）

猜测：g 随 N 递减？

需要更多数据点来拟合 g(ρ, N)。
""")
