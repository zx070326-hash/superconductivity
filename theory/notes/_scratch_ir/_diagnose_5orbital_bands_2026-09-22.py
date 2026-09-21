#!/usr/bin/env python3
"""
诊断：检查 5-orbital 模型的能带结构

如果能带结构本身就错了（比如全在费米面以上/以下），
Lindhard 会给出非物理结果
"""

import numpy as np
from numpy.linalg import eigh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def build_5orbital_hamiltonian(kx, ky):
    """简化 5-轨道模型（与主代码相同）"""
    t_xz = -0.5
    t_yz = -0.5
    t_xy = -0.3

    H = np.zeros((5, 5), dtype=complex)
    H[0, 0] = -2*t_xz*np.cos(kx)
    H[1, 1] = -2*t_yz*np.cos(ky)
    H[2, 2] = -2*t_xy*(np.cos(kx) + np.cos(ky))
    H[3, 3] = -0.2*(np.cos(kx) - np.cos(ky))
    H[4, 4] = -0.1

    t_inter = 0.1
    H[0, 1] = H[1, 0] = t_inter * np.sin(kx) * np.sin(ky)
    H[0, 2] = H[2, 0] = t_inter * np.cos(kx)
    H[1, 2] = H[2, 1] = t_inter * np.cos(ky)

    return H


# 计算沿高对称线的能带
k_points = 100
path_Gamma_X = [(t*np.pi, 0) for t in np.linspace(0, 1, k_points)]
path_X_M = [(np.pi, t*np.pi) for t in np.linspace(0, 1, k_points)]
path_M_Gamma = [(np.pi - t*np.pi, np.pi - t*np.pi) for t in np.linspace(0, 1, k_points)]

full_path = path_Gamma_X + path_X_M + path_M_Gamma
bands = []

for kx, ky in full_path:
    H = build_5orbital_hamiltonian(kx, ky)
    eig, _ = eigh(H)
    bands.append(eig)

bands = np.array(bands)

# 统计
print("=" * 70)
print("5-orbital band structure diagnostics")
print("=" * 70)

for n in range(5):
    band_min = bands[:, n].min()
    band_max = bands[:, n].max()
    band_width = band_max - band_min
    print(f"\nBand {n}:")
    print(f"  Range: [{band_min:.4f}, {band_max:.4f}]")
    print(f"  Width: {band_width:.4f}")

# 总能量范围
E_min = bands.min()
E_max = bands.max()
print(f"\nTotal energy range: [{E_min:.4f}, {E_max:.4f}]")
print(f"Center: {(E_min + E_max)/2:.4f}")

# 推荐的化学势范围
print(f"\nRecommended mu range for partial filling:")
print(f"  mu ∈ [{E_min:.4f}, {E_max:.4f}]")
print(f"  Current mu = 0.0 → {'INSIDE' if E_min < 0 < E_max else 'OUTSIDE'} band range")

# 可视化
fig, ax = plt.subplots(figsize=(8, 6))

for n in range(5):
    ax.plot(bands[:, n], label=f'Band {n}', linewidth=2)

ax.axhline(0, color='red', linestyle='--', linewidth=1, label='mu = 0')
ax.set_xlabel('k-path: Γ → X → M → Γ', fontsize=12)
ax.set_ylabel('Energy', fontsize=12)
ax.set_title('5-orbital band structure (simplified model)', fontsize=14)
ax.legend()
ax.grid(True, alpha=0.3)

# 标记高对称点
n_seg = len(path_Gamma_X)
ax.axvline(0, color='gray', linestyle=':', alpha=0.5)
ax.axvline(n_seg, color='gray', linestyle=':', alpha=0.5)
ax.axvline(2*n_seg, color='gray', linestyle=':', alpha=0.5)
ax.axvline(3*n_seg, color='gray', linestyle=':', alpha=0.5)

ax.text(0, ax.get_ylim()[0], 'Γ', ha='center')
ax.text(n_seg, ax.get_ylim()[0], 'X', ha='center')
ax.text(2*n_seg, ax.get_ylim()[0], 'M', ha='center')
ax.text(3*n_seg, ax.get_ylim()[0], 'Γ', ha='center')

plt.tight_layout()
plt.savefig('_5orbital_bands_diagnostic.png', dpi=150)
print(f"\nBand structure plot saved: _5orbital_bands_diagnostic.png")

print("\n" + "=" * 70)
print("DIAGNOSIS:")
if E_min > 0:
    print("  ERROR: All bands above mu=0 → no occupied states → chi0 = 0")
    print("  FIX: Lower mu or shift bands down")
elif E_max < 0:
    print("  ERROR: All bands below mu=0 → all occupied → no response")
    print("  FIX: Raise mu or shift bands up")
else:
    print("  OK: Fermi level crosses bands → should have Fermi surface")
    print(f"  Crossing bands: {np.sum((bands.min(axis=0) < 0) & (bands.max(axis=0) > 0))}/5")
print("=" * 70)
