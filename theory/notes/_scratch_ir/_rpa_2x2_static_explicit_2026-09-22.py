#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RPA 2x2 static limit explicit calculation
Goal: Find exact relation between D[chi^RPA] and D[chi_0]
"""
import sympy as sp
from sympy import symbols, Matrix, simplify, sqrt, Abs, factor, expand
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

# 定义符号变量（全部实数）
c11, c22, c12 = symbols('c11 c22 c12', real=True, positive=True)
u1, u2 = symbols('u1 u2', real=True, positive=True)

# χ₀ 矩阵（2×2 实对称 PSD）
chi0 = Matrix([
    [c11, c12],
    [c12, c22]
])

# U 矩阵（对角）
U = Matrix([
    [u1, 0],
    [0, u2]
])

print("=" * 60)
print("RPA 2×2 静态极限显式计算")
print("=" * 60)

# 计算 χ₀U
chi0U = chi0 * U
print("\nχ₀U =")
print(chi0U)

# 计算 1 - Uχ₀
I = Matrix([[1, 0], [0, 1]])
Uchi0 = U * chi0
one_minus_Uchi0 = I - Uchi0
print("\n1 - Uχ₀ =")
print(one_minus_Uchi0)

# 计算 (1 - Uχ₀)^{-1}
print("\n计算 (1 - Uχ₀)^{-1}...")
inv_term = one_minus_Uchi0.inv()
print("(1 - Uχ₀)^{-1} =")
for i in range(2):
    for j in range(2):
        print(f"  [{i},{j}]: {simplify(inv_term[i,j])}")

# 计算 χ^RPA = χ₀(1 - Uχ₀)^{-1}
print("\n计算 χ^RPA = χ₀(1 - Uχ₀)^{-1}...")
chi_rpa = chi0 * inv_term
print("χ^RPA =")
for i in range(2):
    for j in range(2):
        chi_rpa[i,j] = simplify(chi_rpa[i,j])
        print(f"  [{i},{j}]: {chi_rpa[i,j]}")

# 提取对角元和非对角元
print("\n" + "=" * 60)
print("对角元和非对角元")
print("=" * 60)

rpa11 = chi_rpa[0,0]
rpa22 = chi_rpa[1,1]
rpa12 = chi_rpa[0,1]

print(f"\nχ^RPA_{{11}} = {rpa11}")
print(f"χ^RPA_{{22}} = {rpa22}")
print(f"χ^RPA_{{12}} = {rpa12}")

# 计算对角主导度量
print("\n" + "=" * 60)
print("对角主导度量")
print("=" * 60)

# 原始 χ₀ 的对角主导度量
D_chi0_1 = c11 - Abs(c12)
D_chi0_2 = c22 - Abs(c12)
D_chi0 = sp.Min(D_chi0_1, D_chi0_2)

print(f"\nD[χ₀] = min(c11 - |c12|, c22 - |c12|)")
print(f"      = min({D_chi0_1}, {D_chi0_2})")

# χ^RPA 的对角主导度量
D_rpa_1 = rpa11 - Abs(rpa12)
D_rpa_2 = rpa22 - Abs(rpa12)

print(f"\nD[χ^RPA] 第1行: χ^RPA_{{11}} - |χ^RPA_{{12}}|")
print(f"简化中...")
D_rpa_1_simplified = simplify(D_rpa_1)
print(f"结果: {D_rpa_1_simplified}")

print(f"\nD[χ^RPA] 第2行: χ^RPA_{{22}} - |χ^RPA_{{12}}|")
print(f"简化中...")
D_rpa_2_simplified = simplify(D_rpa_2)
print(f"结果: {D_rpa_2_simplified}")

# 计算谱半径 ρ(χ₀U)
print("\n" + "=" * 60)
print("谱半径 ρ(χ₀U)")
print("=" * 60)

chi0U_mat = chi0 * U
eigenvals = chi0U_mat.eigenvals()
print(f"\nχ₀U 的本征值（符号形式）:")
for ev, mult in eigenvals.items():
    print(f"  λ = {ev}, 重数 = {mult}")

# 对于 2×2 实对称，谱半径是最大本征值的绝对值
# 但 χ₀U 不对称！需要单独计算
# χ₀U = [[c11*u1, c12*u2], [c12*u1, c22*u2]]
# 特征多项式 det(χ₀U - λI) = 0

print("\n特征多项式展开:")
lam = symbols('lambda', real=True)
char_poly = (chi0U_mat - lam*I).det()
print(f"det(χ₀U - λI) = {expand(char_poly)}")

# 解特征方程
eigenvalues_explicit = sp.solve(char_poly, lam)
print(f"\n本征值（显式）:")
for i, ev in enumerate(eigenvalues_explicit):
    print(f"  λ_{i+1} = {ev}")

# 谱半径（取绝对值最大的）
print("\n对于 ρ < 1 的条件，我们需要:")
print("  max(|λ₁|, |λ₂|) < 1")
print("\n这等价于 Stoner 稳定性条件。")

# 数值测试案例
print("\n" + "=" * 60)
print("数值测试案例")
print("=" * 60)

test_cases = [
    # (c11, c22, c12, u1, u2, 描述)
    (1.0, 1.0, 0.3, 0.4, 0.4, "对称 χ₀, 均匀 U"),
    (1.0, 0.8, 0.2, 0.5, 0.5, "弱非对称 χ₀, 均匀 U"),
    (1.0, 1.0, 0.3, 0.3, 0.5, "对称 χ₀, 非均匀 U"),
    (1.5, 0.8, 0.2, 0.4, 0.4, "强非对称 χ₀, 均匀 U"),
]

for c11_val, c22_val, c12_val, u1_val, u2_val, desc in test_cases:
    print(f"\n案例: {desc}")
    print(f"  χ₀ = [[{c11_val}, {c12_val}], [{c12_val}, {c22_val}]]")
    print(f"  U = diag({u1_val}, {u2_val})")

    # 代入数值
    subs_dict = {c11: c11_val, c22: c22_val, c12: c12_val,
                 u1: u1_val, u2: u2_val}

    # 计算数值结果
    chi_rpa_num = chi_rpa.subs(subs_dict)

    rpa11_num = float(chi_rpa_num[0,0])
    rpa22_num = float(chi_rpa_num[1,1])
    rpa12_num = float(chi_rpa_num[0,1])

    # 对角主导度量
    D_chi0_num = min(c11_val - abs(c12_val), c22_val - abs(c12_val))
    D_rpa_1_num = rpa11_num - abs(rpa12_num)
    D_rpa_2_num = rpa22_num - abs(rpa12_num)
    D_rpa_num = min(D_rpa_1_num, D_rpa_2_num)

    # 谱半径
    chi0U_num = chi0U_mat.subs(subs_dict)
    import numpy as np
    chi0U_np = np.array(chi0U_num.tolist(), dtype=float)
    eigenvalues_np = np.linalg.eigvals(chi0U_np)
    rho = max(abs(eigenvalues_np))

    print(f"\n  ρ(χ₀U) = {rho:.6f}")
    print(f"  D[χ₀] = {D_chi0_num:.6f}")
    print(f"  D[χ^RPA] 第1行 = {D_rpa_1_num:.6f}")
    print(f"  D[χ^RPA] 第2行 = {D_rpa_2_num:.6f}")
    print(f"  D[χ^RPA] = {D_rpa_num:.6f}")
    print(f"  比值 D[χ^RPA]/D[χ₀] = {D_rpa_num/D_chi0_num:.6f}")
    print(f"  理论界 1/(1-ρ) = {1/(1-rho):.6f}")

    if D_rpa_num > D_chi0_num:
        print(f"  ✓ 对角主导**增强** (比值 {D_rpa_num/D_chi0_num:.3f}x)")
    else:
        print(f"  ✗ 对角主导退化")

print("\n" + "=" * 60)
print("结论")
print("=" * 60)
print("""
从 2×2 显式计算可以看出：
1. χ^RPA 的表达式可以完全显式化（分母是 det(1 - Uχ₀)）
2. 对角主导度量 D[χ^RPA] 的符号表达式非常复杂
3. 数值测试表明：ρ < 1 时，D[χ^RPA] > D[χ₀]（对角主导增强）
4. 比值 D[χ^RPA]/D[χ₀] 与 1/(1-ρ) 的关系需要进一步分析

下一步：尝试找到 D[χ^RPA] 的更简单的下界表达式。
""")
