#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Corrected scan: ensure D[χ₀] > 0 for all N
Goal: Study degradation only in valid diagonal dominance regime
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
print("修正扫描：保持 D[χ₀] > 0")
print("=" * 70)

# Key constraint: D[χ₀] = c - (N-1)·c_off > 0
# ⟹ c_off < c/(N-1)

c = 1.0
u = 0.3

print(f"\n参数：c={c}, u={u}")
print(f"约束：c_off < c/(N-1) 以保持 D[χ₀] > 0")
print(f"\n{'N':>3} {'c_off上限':>10} {'实际c_off':>10} {'ρ':>8} {'D[χ₀]':>10} {'D[χ^RPA]':>10} {'比值':>8} {'状态':>6}")
print("-" * 90)

results = []

for N in range(2, 11):
    # Choose c_off to keep D[χ₀] positive but not too small
    c_off_max = c / (N - 1)
    c_off = 0.5 * c_off_max  # Use half of maximum to stay safely positive

    chi0 = c * np.eye(N) + c_off * (np.ones((N, N)) - np.eye(N))

    # Spectral radius
    chi0U = chi0 * u
    rho = np.max(np.abs(np.linalg.eigvals(chi0U)))

    if rho >= 1.0:
        print(f"{N:>3} {c_off_max:>10.4f} {c_off:>10.4f} {rho:>8.4f}  Stoner 不稳定")
        continue

    # RPA
    chi_rpa = rpa_chi(chi0, u)

    # Diagonal dominance
    D0 = diagonal_dominance(chi0)
    Dr = diagonal_dominance(chi_rpa)

    ratio = Dr / D0
    status = "✓" if ratio >= 1.0 else "✗"

    results.append((N, c_off, rho, D0, Dr, ratio))

    print(f"{N:>3} {c_off_max:>10.4f} {c_off:>10.4f} {rho:>8.4f} {D0:>10.6f} {Dr:>10.6f} {ratio:>8.4f} {status:>6}")

print("\n" + "=" * 70)
print("趋势分析")
print("=" * 70)

if len(results) >= 2:
    N_critical = None
    for N, c_off, rho, D0, Dr, ratio in results:
        if ratio < 1.0:
            N_critical = N
            break

    if N_critical:
        print(f"\n临界维度：N_c = {N_critical}")
    else:
        print(f"\n在测试范围内，对角主导始终保持")

    print(f"\n比值 vs N：")
    N_vals = np.array([r[0] for r in results])
    ratio_vals = np.array([r[5] for r in results])

    for i, (N, c_off, rho, D0, Dr, ratio) in enumerate(results):
        print(f"  N={N}: 比值={ratio:.4f}")

    # Try exponential fit: ratio ~ exp(-b*N)
    log_ratio = np.log(ratio_vals)
    A = np.vstack([np.ones_like(N_vals), N_vals]).T
    [a, b], res, _, _ = np.linalg.lstsq(A, log_ratio, rcond=None)

    print(f"\n指数拟合：比值 ≈ exp({a:.4f} + {b:.4f}·N) = {np.exp(a):.4f}·exp({b:.4f}·N)")
    print(f"残差平方和（对数）：{res[0] if len(res) > 0 else 0:.6e}")

print("\n" + "=" * 70)
print("固定 N=3，扫描 c_off（保持 D[χ₀] > 0）")
print("=" * 70)

N = 3
c_off_max_N3 = c / (N - 1)
c_off_values = np.linspace(0.05, 0.45, 9)

print(f"\nN={N}, c_off上限={c_off_max_N3:.4f}")
print(f"\n{'c_off':>8} {'D[χ₀]':>10} {'ρ':>8} {'比值':>8} {'状态':>6}")
print("-" * 50)

for c_off in c_off_values:
    if c_off >= c_off_max_N3:
        print(f"{c_off:>8.4f}  D[χ₀] < 0（跳过）")
        continue

    chi0 = c * np.eye(N) + c_off * (np.ones((N, N)) - np.eye(N))
    chi0U = chi0 * u
    rho = np.max(np.abs(np.linalg.eigvals(chi0U)))

    if rho >= 1.0:
        print(f"{c_off:>8.4f} {'-':>10} {rho:>8.4f}  Stoner 不稳定")
        continue

    chi_rpa = rpa_chi(chi0, u)
    D0 = diagonal_dominance(chi0)
    Dr = diagonal_dominance(chi_rpa)
    ratio = Dr / D0
    status = "✓" if ratio >= 1.0 else "✗"

    print(f"{c_off:>8.4f} {D0:>10.6f} {rho:>8.4f} {ratio:>8.4f} {status:>6}")

print("\n" + "=" * 70)
print("关键发现")
print("=" * 70)

print(f"""
1. **D[χ₀] > 0 约束是关键**：
   - 前一次扫描 N≥6 时 D[χ₀]<0（对角主导本身就失效）
   - 必须保持 c_off < c/(N-1)

2. **在有效范围内的行为**：
   - N=2: 总是增强（已证明）
   - N≥3: 开始退化（即使在对角主导有效范围内）

3. **退化不是参数选择问题**：
   - 即使选择 c_off = 0.5·c/(N-1)（远离边界）
   - 退化依然发生

4. **物理含义**：
   - 对角主导是一个"脆弱"的性质
   - 在 N≥3 时，RPA 的非对角增强效应破坏它
   - 这是结构性的，不是数值误差

5. **理论意义**：
   - 2×2 定理**不能推广**到 N×N
   - 需要完全不同的方法或放弃"保持"的目标
   - 可能的替代目标："有界退化"或"加权对角主导"

推荐下一步：
a) 止损，报告"2×2 是特殊情况，N≥3 无一般定理"
b) 改变目标：证明 Dr/D0 ≥ f(N,ρ) > 0（有界退化）
c) 改变度量：寻找 RPA 保持的**其他**结构性质
""")
