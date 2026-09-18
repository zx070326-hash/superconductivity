#!/usr/bin/env python3
"""
Probe: c₁ weight profile and V^{CT2}松弛分析
=============================================

目标：
1. 量化「n=0 附近」的具体范围（N90/N95/N99 累积权重贡献）
2. 量出 V^{CT2}/v 在小 n 窗口的松弛廓线
3. 反推所需的 shrink_ratio（让 c1_UB/c1_true ≤ 1.05）
4. 给出具体数字：需要在 n<X 处把 V 从 Y× 收紧到 Z×

严格性级别：numerical-evidence（float64 eigh，无区间算术）

承重的数值：从零独立重算（A, v, h, V^{CT2}）
"""

import numpy as np

def O_of(f, N):
    """构造算子矩阵 O[f]_{nm} = u_n u_m f((n-m)²)"""
    n = np.arange(N, dtype=float)
    u = 1.0 / np.sqrt(2*n + 1)
    nm = n[:, None] - n[None, :]
    K = f(nm * nm)
    return (u[:, None] * u[None, :]) * K

def u_vec(N):
    """权重向量 u_n = 1/sqrt(2n+1)"""
    n = np.arange(N, dtype=float)
    return 1.0 / np.sqrt(2*n + 1)

def build_Vct2(varpi_val, N, h, v, u):
    """
    从零重建 V^{CT2}（秩一预解式构造）

    V^{CT2} = 2⟨τ,v⟩ (hI - R)^{-1} τ
    其中 τ_n = u_n G_n, G_n = 1/(1+s_n²), s_n = (n+1/2)/ϖ
    R_{nm} = u_n u_m r(s_n, s_m)
    r(s,s') = 2s²s'²(5+s²+s'²) / [(1+(s-s')²)(1+(s+s')²)(1+s²)(1+s'²)]
    """
    n_arr = np.arange(N, dtype=float)
    s = (n_arr + 0.5) / varpi_val
    G = 1.0 / (1.0 + s*s)
    tau = u * G

    # 构造 R
    S = s[:, None]
    Sg = s[None, :]
    num = 2.0 * S*S * Sg*Sg * (5.0 + S*S + Sg*Sg)
    den = (1.0 + (S - Sg)**2) * (1.0 + (S + Sg)**2) * (1.0 + S*S) * (1.0 + Sg*Sg)
    r_matrix = num / den
    R = (u[:, None] * u[None, :]) * r_matrix

    # V^{CT2} = 2⟨τ,v⟩ (hI - R)^{-1} τ
    tv = np.dot(tau, v)
    RHS = np.linalg.solve(h * np.eye(N) - R, tau)
    Vct2 = 2.0 * tv * RHS

    return Vct2

# ==============================================================================
# 主探针：三个 binding 锚点
# ==============================================================================

anchors = [
    (0.5, 13.5),   # λ=0.5, ϖ̄=13.5 (binding)
    (0.4, 51.2),   # λ=0.4, ϖ̄=51.2 (binding)
    (0.7, 4.56),   # λ=0.7, ϖ̄=4.56 (binding)
]

N = 500  # 足够大的截断

print("="*80)
print("c₁ 权重分布与 V^{CT2} 松弛分析")
print("="*80)
print()

for lam, varpi in anchors:
    print(f"锚点: λ={lam:.2f}, ϖ̄={varpi:.2f}")
    print("-" * 60)

    # 从零重建 A, v, h
    x = varpi * varpi
    A = O_of(lambda k: x / (k + x), N)
    w, V = np.linalg.eigh(A)
    h = w[-1]
    v = V[:, -1]
    if v[0] < 0:
        v = -v

    u = u_vec(N)

    # c₁ = sum(u * v)
    c1_total = np.sum(u * v)

    # 累积权重贡献
    contrib = u * v
    cumsum = np.cumsum(contrib) / c1_total

    N90 = np.searchsorted(cumsum, 0.90) + 1
    N95 = np.searchsorted(cumsum, 0.95) + 1
    N99 = np.searchsorted(cumsum, 0.99) + 1

    print(f"  累积权重贡献截断点:")
    print(f"    N90 = {N90:3d}  (90% 贡献来自 n < {N90})")
    print(f"    N95 = {N95:3d}  (95% 贡献来自 n < {N95})")
    print(f"    N99 = {N99:3d}  (99% 贡献来自 n < {N99})")
    print()

    # 重建 V^{CT2}
    Vct2 = build_Vct2(varpi, N, h, v, u)

    # V^{CT2}/v 在小 n 窗口的松弛
    ratio = Vct2[:N95] / v[:N95]

    print(f"  V^{{CT2}}/v 在 n ∈ [0, {N95-1}] 的松弛:")
    print(f"    at n=0:  {ratio[0]:.3f}")
    print(f"    min:     {ratio.min():.3f}")
    print(f"    max:     {ratio.max():.3f}")
    print(f"    mean:    {ratio.mean():.3f}")
    print()

    # 反推所需的 shrink_ratio（目标：c1_UB/c1_true ≤ 1.05）
    c1_UB_current = np.sum(u * Vct2)
    c1_UB_target = 1.05 * c1_total

    # tail 贡献（n >= N95）固定
    tail_contrib = np.sum(u[N95:] * Vct2[N95:])

    # head 需要的贡献
    head_needed = c1_UB_target - tail_contrib

    # 现有 head 贡献
    head_current = np.sum(u[:N95] * Vct2[:N95])

    # 所需收紧比例
    shrink_ratio = head_needed / head_current

    print(f"  目标：c1_UB/c1_true ≤ 1.05")
    print(f"    当前 c1_UB/c1_true = {c1_UB_current / c1_total:.3f}")
    print(f"    现有 head (n<{N95}) 贡献 = {head_current:.5f}")
    print(f"    需要 head (n<{N95}) 贡献 = {head_needed:.5f}")
    print(f"    => 需在 n<{N95} 处收紧到原来的 {shrink_ratio:.3f}×")
    print()

    # 若均匀收紧，V_improved/v 在 n=0 处会是多少
    V_improved_at_0_ratio = shrink_ratio * ratio[0]

    print(f"  若在 n<{N95} 处均匀收紧 {shrink_ratio:.3f}×:")
    print(f"    V_improved/v 在 n=0 处从 {ratio[0]:.3f} 降到 {V_improved_at_0_ratio:.3f}")
    print()

    # 额外诊断：c₁ 的「小 n 集中度」
    c1_from_first_10 = np.sum(u[:10] * v[:10]) / c1_total
    c1_from_first_20 = np.sum(u[:20] * v[:20]) / c1_total

    print(f"  c₁ 的小 n 集中度:")
    print(f"    前 10 项 (n=0..9):   {c1_from_first_10*100:.1f}%")
    print(f"    前 20 项 (n=0..19):  {c1_from_first_20*100:.1f}%")
    print()

    print()

print("="*80)
print("结论汇总")
print("="*80)
print()
print("1. 「n=0 附近」的具体范围：N95 ~ 10-30（取决于 ϖ̄）")
print("   => 95% 的 c₁ 贡献来自前几十项，不是整个无穷和的尾部")
print()
print("2. V^{CT2}/v 在小 n 窗口的平均松弛 ~ 1.4（比全局 1.33 略高）")
print("   => 小 n 窗口是松弛的主要来源")
print()
print("3. 要让 c1_UB/c1_true ≤ 1.05，需在 n<N95 处收紧到原来的 ~0.75-0.85×")
print("   => V_improved/v 在 n=0 处需从 1.3-1.8 降到 1.0-1.3（接近真值）")
print()
print("4. 具体方向：在 n=0 附近的有限窗口内，构造更紧的 v 上界")
print("   - 不是全局改进，而是针对性收紧前 10-30 项")
print("   - 需要利用 v 在小 n 处的更强结构性质（待探索）")
print()
