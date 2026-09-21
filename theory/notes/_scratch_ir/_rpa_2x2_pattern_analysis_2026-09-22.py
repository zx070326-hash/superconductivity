#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2x2 RPA diagonal dominance pattern analysis
Goal: Extract pattern from numerical tests to guide general proof
"""
import numpy as np
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

print("=" * 60)
print("2×2 RPA 对角主导模式分析")
print("=" * 60)

# Test cases from explicit calculation
test_cases = [
    # (c11, c22, c12, u1, u2, description)
    (1.0, 1.0, 0.3, 0.4, 0.4, "对称 χ₀, 均匀 U"),
    (1.0, 0.8, 0.2, 0.5, 0.5, "弱非对称 χ₀, 均匀 U"),
    (1.0, 1.0, 0.3, 0.3, 0.5, "对称 χ₀, 非均匀 U"),
    (1.5, 0.8, 0.2, 0.4, 0.4, "强非对称 χ₀, 均匀 U"),
]

print("\n模式提取：对角主导比值 vs ρ")
print("-" * 60)

for c11, c22, c12, u1, u2, desc in test_cases:
    # Build matrices
    chi0 = np.array([[c11, c12], [c12, c22]])
    U = np.diag([u1, u2])

    # Compute chi0U and its spectral radius
    chi0U = chi0 @ U
    eigenvalues = np.linalg.eigvals(chi0U)
    rho = np.max(np.abs(eigenvalues))

    # Compute RPA chi
    I = np.eye(2)
    Uchi0 = U @ chi0
    chi_rpa = chi0 @ np.linalg.inv(I - Uchi0)

    # Diagonal dominance measures
    D_chi0_1 = c11 - abs(c12)
    D_chi0_2 = c22 - abs(c12)
    D_chi0 = min(D_chi0_1, D_chi0_2)

    D_rpa_1 = chi_rpa[0,0] - abs(chi_rpa[0,1])
    D_rpa_2 = chi_rpa[1,1] - abs(chi_rpa[0,1])
    D_rpa = min(D_rpa_1, D_rpa_2)

    # Key ratios
    ratio_D = D_rpa / D_chi0
    theoretical_bound = 1 / (1 - rho)

    # Analyze which row is bottleneck
    bottleneck = 1 if D_rpa_1 < D_rpa_2 else 2

    print(f"\n{desc}")
    print(f"  ρ = {rho:.6f}")
    print(f"  D[χ₀] = {D_chi0:.6f}")
    print(f"  D[χ^RPA] = {D_rpa:.6f}")
    print(f"  比值 D[χ^RPA]/D[χ₀] = {ratio_D:.6f}")
    print(f"  理论界 1/(1-ρ) = {theoretical_bound:.6f}")
    print(f"  瓶颈行 = {bottleneck}")

    # Check if diagonal dominance is preserved
    if ratio_D >= 1:
        print(f"  ✓ 对角主导增强")
    else:
        print(f"  ✗ 对角主导退化")

    # Analyze the failure mechanism for case 4
    if ratio_D < 1:
        print(f"\n  失效分析：")
        print(f"    χ₀ 对角元比值 = {c11/c22:.3f}")
        print(f"    χ^RPA 对角元：{chi_rpa[0,0]:.3f} vs {chi_rpa[1,1]:.3f}")
        print(f"    行 1 对角主导：{D_rpa_1:.3f}")
        print(f"    行 2 对角主导：{D_rpa_2:.3f}")
        print(f"    瓶颈在较弱的对角元（行 {bottleneck}）")

print("\n" + "=" * 60)
print("关键模式总结")
print("=" * 60)

print("""
从 4 个测试案例可以看出：

1. **均匀 U + 对称 χ₀**：对角主导总是增强（案例 1）
   - ρ = 0.52 → 比值 1.39（增强 39%）

2. **均匀 U + 弱非对称 χ₀**：对角主导仍增强（案例 2）
   - ρ = 0.56 → 比值 1.26（增强 26%）

3. **非均匀 U + 对称 χ₀**：对角主导轻微增强（案例 3）
   - ρ = 0.55 → 比值 1.04（增强 4%）
   - U 的不均匀性削弱了增强效果

4. **均匀 U + 强非对称 χ₀**：对角主导退化（案例 4）
   - ρ = 0.62 → 比值 0.85（退化 15%）
   - χ₀ 对角元比值 1.5/0.8 = 1.875（差异大）
   - 瓶颈在较弱的对角元（行 2）

**失效机制**：
- RPA 对强耦合轨道增强更多，对弱耦合轨道增强有限
- 当对角元差异大时，瓶颈在弱轨道
- min 操作放大了不对称效应

**可能的定理方向**：
1. **条件定理**：若 χ₀ 对角元接近（κ := max/min ≤ C），则对角主导保持
2. **加权度量**：用加权和代替 min，避免木桶效应
3. **相对度量**：证明 D[χ^RPA]/max(χ^RPA_{αα}) ≥ D[χ₀]/max(χ₀_{αα})
""")

print("\n" + "=" * 60)
print("下一步：探索条件数界")
print("=" * 60)

# Scan condition number vs diagonal dominance preservation
print("\n扫描对角元条件数 κ = c11/c22 的影响：")
print("（固定 c12=0.2, u1=u2=0.4, c22=0.8）")
print()

kappa_values = [1.0, 1.2, 1.5, 1.8, 2.0, 2.5]
c22_fixed = 0.8
c12_fixed = 0.2
u_fixed = 0.4

results = []

for kappa in kappa_values:
    c11 = kappa * c22_fixed
    chi0 = np.array([[c11, c12_fixed], [c12_fixed, c22_fixed]])
    U = np.diag([u_fixed, u_fixed])

    chi0U = chi0 @ U
    rho = np.max(np.abs(np.linalg.eigvals(chi0U)))

    I = np.eye(2)
    Uchi0 = U @ chi0
    chi_rpa = chi0 @ np.linalg.inv(I - Uchi0)

    D_chi0 = min(c11 - abs(c12_fixed), c22_fixed - abs(c12_fixed))
    D_rpa = min(chi_rpa[0,0] - abs(chi_rpa[0,1]),
                chi_rpa[1,1] - abs(chi_rpa[0,1]))

    ratio = D_rpa / D_chi0
    preserved = ratio >= 1.0

    results.append((kappa, rho, ratio, preserved))

    status = "✓ 保持" if preserved else "✗ 退化"
    print(f"κ = {kappa:.2f}:  ρ = {rho:.4f},  比值 = {ratio:.4f}  {status}")

# Find the critical kappa
print("\n寻找临界条件数 κ_c（对角主导刚好保持的边界）：")
preserved_kappas = [k for k, r, ratio, p in results if p]
failed_kappas = [k for k, r, ratio, p in results if not p]

if preserved_kappas and failed_kappas:
    kappa_lower = max(preserved_kappas)
    kappa_upper = min(failed_kappas)
    print(f"  κ_c ∈ [{kappa_lower:.2f}, {kappa_upper:.2f}]")
    print(f"  （在这个参数下）")
else:
    print("  需要更密集的扫描")

print("\n" + "=" * 60)
print("结论")
print("=" * 60)
print("""
1. 对角主导保持与否**不只依赖 ρ**，还依赖 χ₀ 的结构
2. 对角元条件数 κ 是关键参数：κ 大 → 容易失败
3. 可能的条件定理：
   - 若 κ ≤ κ_c(ρ)，则 D[χ^RPA] ≥ D[χ₀]
   - 需要找到 κ_c(ρ) 的解析表达式

下一步：
- 用符号计算展开 2×2 的 D[χ^RPA] 表达式
- 寻找 κ 和 ρ 的显式关系
- 推广到 N×N（若可能）
""")
