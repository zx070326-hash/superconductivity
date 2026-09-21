#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Binary search for critical condition number kappa_c
Goal: Find kappa_c(rho) relation
"""
import numpy as np
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def compute_diagonal_dominance_ratio(c11, c22, c12, u1, u2):
    """
    Compute D[chi_RPA] / D[chi_0] for given parameters
    Returns: (rho, D_chi0, D_rpa, ratio)
    """
    chi0 = np.array([[c11, c12], [c12, c22]])
    U = np.diag([u1, u2])

    # Spectral radius
    chi0U = chi0 @ U
    rho = np.max(np.abs(np.linalg.eigvals(chi0U)))

    # RPA chi
    I = np.eye(2)
    Uchi0 = U @ chi0
    chi_rpa = chi0 @ np.linalg.inv(I - Uchi0)

    # Diagonal dominance
    D_chi0 = min(c11 - abs(c12), c22 - abs(c12))
    D_rpa = min(chi_rpa[0,0] - abs(chi_rpa[0,1]),
                chi_rpa[1,1] - abs(chi_rpa[0,1]))

    ratio = D_rpa / D_chi0

    return rho, D_chi0, D_rpa, ratio

def find_critical_kappa(c22, c12, u, tol=1e-6):
    """
    Binary search for critical kappa where D[chi_RPA] / D[chi_0] = 1
    Fixed: c22, c12, u (uniform)
    Variable: c11 = kappa * c22
    """
    # Start with wide bracket
    kappa_low = 1.0  # Symmetric case (always works)
    kappa_high = 3.0  # Very asymmetric (likely fails)

    # Check endpoints
    c11_low = kappa_low * c22
    _, _, _, ratio_low = compute_diagonal_dominance_ratio(c11_low, c22, c12, u, u)

    c11_high = kappa_high * c22
    _, _, _, ratio_high = compute_diagonal_dominance_ratio(c11_high, c22, c12, u, u)

    if ratio_low < 1.0:
        return None, "Even symmetric case fails"

    if ratio_high >= 1.0:
        return None, "Even kappa=3.0 works, need higher upper bound"

    # Binary search
    while kappa_high - kappa_low > tol:
        kappa_mid = (kappa_low + kappa_high) / 2
        c11_mid = kappa_mid * c22
        rho, _, _, ratio = compute_diagonal_dominance_ratio(c11_mid, c22, c12, u, u)

        if ratio >= 1.0:
            kappa_low = kappa_mid
        else:
            kappa_high = kappa_mid

    kappa_c = (kappa_low + kappa_high) / 2
    c11_c = kappa_c * c22
    rho_c, _, _, ratio_c = compute_diagonal_dominance_ratio(c11_c, c22, c12, u, u)

    return kappa_c, rho_c

print("=" * 60)
print("精确定位临界条件数 κ_c")
print("=" * 60)

# Fixed parameters
c22 = 0.8
c12 = 0.2
u = 0.4

print(f"\n固定参数：c22 = {c22}, c12 = {c12}, u = {u}")
print(f"（搜索 c11 = κ · c22，使得 D[χ^RPA] / D[χ₀] = 1）")

kappa_c, rho_c = find_critical_kappa(c22, c12, u)

if kappa_c is not None:
    print(f"\n临界条件数：κ_c = {kappa_c:.6f}")
    print(f"对应谱半径：ρ_c = {rho_c:.6f}")
    print(f"\n验证：")

    # Test slightly below and above
    for delta in [-0.01, 0.0, 0.01]:
        kappa_test = kappa_c + delta
        c11_test = kappa_test * c22
        rho, D0, Dr, ratio = compute_diagonal_dominance_ratio(c11_test, c22, c12, u, u)
        status = "✓" if ratio >= 1.0 else "✗"
        print(f"  κ = {kappa_test:.4f}: ρ = {rho:.4f}, 比值 = {ratio:.6f}  {status}")

print("\n" + "=" * 60)
print("扫描不同参数下的 κ_c(ρ) 关系")
print("=" * 60)

# Vary u to get different rho values
u_values = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
results = []

print(f"\n固定 c22={c22}, c12={c12}，改变 u：")
print(f"{'u':>6} {'ρ_c':>8} {'κ_c':>8}")
print("-" * 24)

for u_val in u_values:
    kappa_c, rho_c = find_critical_kappa(c22, c12, u_val)
    if kappa_c is not None:
        results.append((u_val, rho_c, kappa_c))
        print(f"{u_val:>6.2f} {rho_c:>8.4f} {kappa_c:>8.4f}")

# Try to fit kappa_c(rho)
if len(results) >= 3:
    print("\n" + "=" * 60)
    print("寻找 κ_c(ρ) 的函数关系")
    print("=" * 60)

    rho_arr = np.array([r[1] for r in results])
    kappa_arr = np.array([r[2] for r in results])

    # Try linear fit: kappa_c = a + b * rho
    A = np.vstack([np.ones_like(rho_arr), rho_arr]).T
    [a, b], residuals, _, _ = np.linalg.lstsq(A, kappa_arr, rcond=None)

    print(f"\n线性拟合：κ_c ≈ {a:.4f} + {b:.4f} · ρ")
    print(f"残差平方和：{residuals[0] if len(residuals) > 0 else 0:.6e}")

    print(f"\n验证：")
    print(f"{'ρ实际':>8} {'κ_c实际':>10} {'κ_c拟合':>10} {'误差%':>8}")
    print("-" * 36)
    for u_val, rho_val, kappa_val in results:
        kappa_fit = a + b * rho_val
        error = abs(kappa_fit - kappa_val) / kappa_val * 100
        print(f"{rho_val:>8.4f} {kappa_val:>10.4f} {kappa_fit:>10.4f} {error:>7.2f}%")

    # Try other functional forms
    print("\n" + "-" * 60)

    # Power law: kappa_c = a * rho^b
    log_rho = np.log(rho_arr)
    log_kappa = np.log(kappa_arr)
    A_log = np.vstack([np.ones_like(log_rho), log_rho]).T
    [log_a, b_pow], res_pow, _, _ = np.linalg.lstsq(A_log, log_kappa, rcond=None)
    a_pow = np.exp(log_a)

    print(f"幂律拟合：κ_c ≈ {a_pow:.4f} · ρ^{b_pow:.4f}")
    print(f"残差平方和（对数）：{res_pow[0] if len(res_pow) > 0 else 0:.6e}")

    # Rational: kappa_c = a / (1 - b*rho)
    y_inv = 1 / kappa_arr
    A_rat = np.vstack([np.ones_like(rho_arr), rho_arr]).T
    [c, d], res_rat, _, _ = np.linalg.lstsq(A_rat, y_inv, rcond=None)
    a_rat = 1 / c
    b_rat = d / c

    print(f"有理式拟合：κ_c ≈ {a_rat:.4f} / (1 - {b_rat:.4f}·ρ)")
    print(f"残差平方和（倒数）：{res_rat[0] if len(res_rat) > 0 else 0:.6e}")

print("\n" + "=" * 60)
print("物理解释")
print("=" * 60)

print("""
关键发现：

1. **临界条件数存在**：对于给定的 (c12, u)，存在临界 κ_c
   - κ < κ_c → 对角主导保持
   - κ > κ_c → 对角主导退化

2. **κ_c 与 ρ 的关系**：κ_c 随 ρ 递减
   - ρ 小（弱耦合）→ κ_c 大（允许更大的非对称）
   - ρ 大（强耦合）→ κ_c 小（要求更接近对称）

3. **物理图像**：
   - ρ 控制 RPA 增强的整体强度
   - κ 控制不同轨道增强的不均匀性
   - 两者竞争：ρ↑ 使增强变强但也使不对称加剧

4. **定理形式**：
   若 κ[χ₀] ≤ κ_c(ρ)，则 D[χ^RPA] ≥ D[χ₀]

   其中 κ_c(ρ) 可以用线性或有理式近似。
""")
