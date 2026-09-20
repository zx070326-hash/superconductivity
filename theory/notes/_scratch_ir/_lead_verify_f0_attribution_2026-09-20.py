# -*- coding: utf-8 -*-
"""
lead 复核 team-lead 指出的 `f(0)` 归属细节（2026-09-20），并量化它在双核/单核
λ₁ 差里占多大份额。

两种写法：
  (A) 笔记约定：L1 从 k=1 起求和、L2 去掉对角；
  (B) builder 约定（`_chk_ml_absM2.py`）：L1 含 f(0)、L2 含对角。

对角元之差 (B) − (A) = u_n² · ( f_{L2}(0) − f_{L1}(0) )：
  * 单核（f_{L1} = f_{L2}）⟹ 恒等于 0，两写法等价；
  * 双核（L1←Λ^Z, L2+L3←Λ^Δ）⟹ = u_n²·(Λ^Δ(0) − Λ^Z(0)) = −2u_n²·λ_sf ≠ 0，
    必须按物理来源分别归属。

本脚本不 import 任何仓库脚本；算三个算子：
  a  = 物理归属（L1 整条来自 Λ^Z，含 Λ^Z(0)）        —— 仓库的 0.115511
  b  = 错误归属（行和来自 Λ^Z，但 f(0) 取了 Λ^Δ 的）
  c  = 单核 𝒪[Λ^Δ]                                  —— 线 C 的 0.123471
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh


def O_b(f_rowsum, f_pair, f0_src, N):
    """约定 (B)：L2/L3 含对角；L1 的 f(0) 由 f0_src 提供、k>=1 的行和由 f_rowsum 提供。"""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    M = np.outer(u, u) * (f_pair[D] + f_pair[S])
    zrow = f0_src[0] + 2.0 * np.concatenate(([0.0], np.cumsum(f_rowsum[1:N])))
    M[n, n] -= u ** 2 * zrow
    return M


def einstein(lam, wbar, K):
    k = np.arange(K, dtype=float)
    return lam * wbar ** 2 / (k ** 2 + wbar ** 2)


W_PH, LAM_PH, W_SF = 0.1, 10.0, 5.0
LAM_SF = LAM_PH * W_PH ** 2 / W_SF ** 2

print("=== f(0) 归属对 lam1 的影响（反例参数 w_ph=0.1, lam_ph=10, w_sf=5）===")
print(f"  lam_sf = {LAM_SF:.10g}   M2^Delta = {LAM_PH*W_PH**2 - LAM_SF*W_SF**2:.3e}")
print()
print("   N        a: 物理归属      b: f(0) 错归属     c: 单核 O[Lam^D]   (a-b)        (b-c)")
for N in (400, 800, 1600, 3200):
    K = 2 * N + 2
    ph, sf = einstein(LAM_PH, W_PH, K), einstein(LAM_SF, W_SF, K)
    Z, Dl = ph + sf, ph - sf
    a = float(eigvalsh(O_b(Z, Dl, Z, N))[-1])
    b = float(eigvalsh(O_b(Z, Dl, Dl, N))[-1])
    c = float(eigvalsh(O_b(Dl, Dl, Dl, N))[-1])
    print(f" {N:5d}   {a:.10f}     {b:.10f}     {c:.10f}   {a-b:+.3e}  {b-c:+.3e}")


def O_a_naive(f_rowsum, f_pair, N):
    """约定 (A) 的**直译**：L1 只从 k=1 起（无 f(0)）、L2 去掉对角（丢掉 f(0)）。
    双核下这等于「两处 f(0) 双双丢弃」——应与上面的 b 逐位相同。"""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    M = np.outer(u, u) * (np.where(D == 0, 0.0, f_pair[D]) + f_pair[S])
    M[n, n] -= 2.0 * u ** 2 * np.concatenate(([0.0], np.cumsum(f_rowsum[1:N])))
    return M


print()
print("  中间态 b 的第二个拼法（(1-delta) 简写直译到双核 = 两处 f(0) 双双丢弃）应与 b 逐位相同：")
for N in (400, 3200):
    K = 2 * N + 2
    ph, sf = einstein(LAM_PH, W_PH, K), einstein(LAM_SF, W_SF, K)
    Z, Dl = ph + sf, ph - sf
    b2 = float(eigvalsh(O_a_naive(Z, Dl, N))[-1])
    b = float(eigvalsh(O_b(Z, Dl, Dl, N))[-1])
    print(f"    N={N:5d}  b2 = {b2:.12f}   b = {b:.12f}   |b2-b| = {abs(b2-b):.3e}")

print()
print("  解析核对：对角元之差应为 u_n^2 * (Lam^D(0) - Lam^Z(0)) = -2 u_n^2 lam_sf")
print(f"    Lam^D(0) - Lam^Z(0) = {Dl[0]-Z[0]:.6f}    -2*lam_sf = {-2*LAM_SF:.6f}")
print()
print("  share of the (single-kernel vs dual) gap carried by the f(0) attribution:",
      f"|a-b|/|a-c| = {abs(a-b)/abs(a-c)*100:.1f}%")
