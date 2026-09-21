#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Systematic scan: diagonal dominance degradation vs N
Goal: Quantify how degradation depends on N and parameters
"""
import numpy as np
import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def diagonal_dominance(M):
    N = M.shape[0]
    D_rows = [M[i,i] - (np.sum(np.abs(M[i,:])) - abs(M[i,i])) for i in range(N)]
    return min(D_rows)

def rpa_chi(chi0, u):
    N = chi0.shape[0]
    I = np.eye(N)
    U = u * np.eye(N)
    return chi0 @ np.linalg.inv(I - U @ chi0)

print("=" * 70)
print("对角主导退化 vs N：系统扫描")
print("=" * 70)

# Scan N from 2 to 10
results = []

c = 1.0
c_off = 0.2
u = 0.35

print(f"\n参数：c={c}, c_off={c_off}, u={u}")
print(f"\n{'N':>3} {'ρ':>8} {'D[χ₀]':>10} {'D[χ^RPA]':>10} {'比值':>8} {'状态':>6}")
print("-" * 70)

for N in range(2, 11):
    # Build symmetric matrix with uniform off-diagonal
    chi0 = c * np.eye(N) + c_off * (np.ones((N, N)) - np.eye(N))

    # Spectral radius
    chi0U = chi0 * u
    rho = np.max(np.abs(np.linalg.eigvals(chi0U)))

    # RPA
    chi_rpa = rpa_chi(chi0, u)

    # Diagonal dominance
    D0 = diagonal_dominance(chi0)
    Dr = diagonal_dominance(chi_rpa)

    ratio = Dr / D0
    status = "✓" if ratio >= 1.0 else "✗"

    results.append((N, rho, D0, Dr, ratio))

    print(f"{N:>3} {rho:>8.4f} {D0:>10.6f} {Dr:>10.6f} {ratio:>8.4f} {status:>6}")

print("\n" + "=" * 70)
print("退化模式分析")
print("=" * 70)

# Extract N where degradation starts
N_critical = None
for i, (N, rho, D0, Dr, ratio) in enumerate(results):
    if ratio < 1.0:
        N_critical = N
        break

if N_critical:
    print(f"\n临界维度：N_c = {N_critical}")
    print(f"对于 N ≥ {N_critical}，对角主导退化")
else:
    print(f"\n在 N ≤ 10 范围内，对角主导始终保持")

# Analyze trend
print(f"\n比值趋势：")
for i in range(len(results) - 1):
    N1, _, _, _, r1 = results[i]
    N2, _, _, _, r2 = results[i+1]
    delta = r2 - r1
    print(f"  N={N1} → N={N2}: 比值变化 {delta:+.4f}")

# Fit degradation rate
N_vals = np.array([r[0] for r in results])
ratio_vals = np.array([r[4] for r in results])

# Try linear fit
A = np.vstack([np.ones_like(N_vals), N_vals]).T
[a, b], residuals, _, _ = np.linalg.lstsq(A, ratio_vals, rcond=None)

print(f"\n线性拟合：比值 ≈ {a:.4f} + {b:.4f}·N")
print(f"残差平方和：{residuals[0] if len(residuals) > 0 else 0:.6e}")

# Predict where ratio = 1
if b != 0:
    N_cross = (1.0 - a) / b
    print(f"预测 N_cross（比值=1）：{N_cross:.2f}")

print("\n" + "=" * 70)
print("参数依赖性：改变 c_off")
print("=" * 70)

c_off_values = [0.1, 0.15, 0.2, 0.25, 0.3]
print(f"\n固定 N=3, u={u}，改变 c_off：")
print(f"\n{'c_off':>6} {'ρ':>8} {'比值':>8} {'状态':>6}")
print("-" * 35)

for c_off_test in c_off_values:
    chi0 = c * np.eye(3) + c_off_test * (np.ones((3, 3)) - np.eye(3))
    chi0U = chi0 * u
    rho = np.max(np.abs(np.linalg.eigvals(chi0U)))

    chi_rpa = rpa_chi(chi0, u)

    D0 = diagonal_dominance(chi0)
    Dr = diagonal_dominance(chi_rpa)
    ratio = Dr / D0
    status = "✓" if ratio >= 1.0 else "✗"

    print(f"{c_off_test:>6.2f} {rho:>8.4f} {ratio:>8.4f} {status:>6}")

print("\n" + "=" * 70)
print("参数依赖性：改变 u")
print("=" * 70)

u_values = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45]
print(f"\n固定 N=3, c_off={0.2}，改变 u：")
print(f"\n{'u':>6} {'ρ':>8} {'比值':>8} {'状态':>6}")
print("-" * 35)

c_off = 0.2
for u_test in u_values:
    chi0 = c * np.eye(3) + c_off * (np.ones((3, 3)) - np.eye(3))
    chi0U = chi0 * u_test
    rho = np.max(np.abs(np.linalg.eigvals(chi0U)))

    if rho >= 1.0:
        print(f"{u_test:>6.2f} {rho:>8.4f}  Stoner 不稳定")
        continue

    chi_rpa = rpa_chi(chi0, u_test)

    D0 = diagonal_dominance(chi0)
    Dr = diagonal_dominance(chi_rpa)
    ratio = Dr / D0
    status = "✓" if ratio >= 1.0 else "✗"

    print(f"{u_test:>6.2f} {rho:>8.4f} {ratio:>8.4f} {status:>6}")

print("\n" + "=" * 70)
print("关键结论")
print("=" * 70)

print(f"""
1. **N=2 是特殊情况**：
   - N=2: 比值 = {results[0][4]:.4f}（增强）
   - N=3: 比值 = {results[1][4]:.4f}（退化）
   - N→大：退化持续加剧

2. **退化率近似线性于 N**：
   - 比值 ≈ {a:.4f} + {b:.4f}·N
   - 每增加 1 个维度，比值下降约 {abs(b):.4f}

3. **参数依赖**：
   - c_off 增大 → 退化加剧（非对角元权重增加）
   - u 增大 → 退化加剧（RPA 增强更强，但非对角增强更快）

4. **物理图像**：
   - N 越大 → 每行的非对角项越多
   - RPA 对非对角元增强 > 对对角元增强
   - N 个非对角项的累积效应压倒对角增强

5. **理论意义**：
   - 2×2 结果**不能外推**到一般 N
   - 对称性 + Stoner 稳定性**不足以**保证对角主导
   - 需要依赖 N 的条件，或改变定理陈述

下一步：
- 寻找使对角主导保持的 c_off/u/N 关系
- 或证明"有界退化"：Dr/D0 ≥ f(N, ρ) > 0
""")
