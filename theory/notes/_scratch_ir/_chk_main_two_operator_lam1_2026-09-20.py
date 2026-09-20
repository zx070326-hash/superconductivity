# -*- coding: utf-8 -*-
"""
主会话从零复核 guidance-lead 对线 C 的推翻（2026-09-20）。

争议：仓库三处笔记写 `M2^Delta = 0` 的定向反例给 lam1 = 0.115511；线 C 算得 0.123471
并据此报「差 6.9%、未闭合」；guidance-lead 判定两个数是**两个不同算子**，仓库无需订正。

本脚本不 import 任何仓库脚本（不用 _chk_ml_build / _chk_ml_absM2），只用
P4_main_theorems.md §0 的定义从零重建：
    u_n = (2n+1)^(-1/2),  n = 0..N-1
    L1[f] 对角 = (2n+1)^(-1) * ( f(0) + 2*sum_{k=1..n} f(k) )
    L2[f]      = f(|n-m|) u_n u_m
    L3[f]      = f(n+m+1) u_n u_m
    O[f] = -L1[f] + L2[f] + L3[f]
（注：笔记里「L1 从 k=1 起」+「L2 去掉对角」是同一算子的另一写法——单核时两处的
 f(0) 恰好抵消。两核时不能再抵消，必须按物理来源分别归属：Z 通道给 L1、配对通道给
 L2+L3，这正是下面 dual 的写法。）

校准：O[1/k^2] 的顶本征值应复现仓库 g(2) = 1.3181402121。
"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh


def O_general(f_diag_src, f_pair_src, N, drop_f0=False):
    """f_diag_src, f_pair_src: 长度 >= 2N 的数组，第 k 项是核在整数 k 上的值。

    drop_f0=True 复现「照笔记 (1−δ) 约定直译」的写法：L1 从 k=1 起求和、L2 去掉对角
    ⟹ **两处 f(0) 都被丢弃**。单核时它与 drop_f0=False 等价（f(0) 精确抵消），
    双核时不等价——见文末 Stage C。
    """
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    pair = f_pair_src[D].copy()
    if drop_f0:
        pair[n, n] = 0.0
    M = np.outer(u, u) * (pair + f_pair_src[S])
    zrow = (0.0 if drop_f0 else f_diag_src[0]) \
        + 2.0 * np.concatenate(([0.0], np.cumsum(f_diag_src[1:N])))
    M[n, n] -= u ** 2 * zrow
    return M


def einstein(lam, wbar, K):
    k = np.arange(K, dtype=float)
    return lam * wbar ** 2 / (k ** 2 + wbar ** 2)


def lam1(M):
    return float(eigvalsh(M)[-1])


print("=== 校准：O[1/k^2] 的 lam1 应 = g(2) = 1.3181402121 ===")
for N in (400, 1600, 3200):
    K = 2 * N + 2
    k = np.arange(K, dtype=float)
    f = np.zeros(K)
    f[1:] = 1.0 / k[1:] ** 2
    f[0] = 0.0          # f(0)=1/0^2 发散，但单核时 f(0) 在 -L1 与 L2 之间精确抵消
    print(f"  N={N:5d}  lam1(O[1/k^2]) = {lam1(O_general(f, f, N)):.10f}")

print()
print("=== 反例参数：w_ph=0.1, lam_ph=10 / w_sf=5, lam_sf = lam_ph*w_ph^2/w_sf^2 ===")
W_PH, LAM_PH, W_SF = 0.1, 10.0, 5.0
LAM_SF = LAM_PH * W_PH ** 2 / W_SF ** 2
print(f"  lam_sf = {LAM_SF:.10g}")
print(f"  M2^ph = {LAM_PH*W_PH**2:.10g}   M2^sf = {LAM_SF*W_SF**2:.10g}"
      f"   M2^Delta = {LAM_PH*W_PH**2 - LAM_SF*W_SF**2:.3e}")

print()
print("  N      lam1( O[Lam^Delta] )   lam1( dual: L1<-Lam^Z, L2+L3<-Lam^Delta )")
for N in (100, 200, 400, 800, 1600, 3200):
    K = 2 * N + 2
    ph, sf = einstein(LAM_PH, W_PH, K), einstein(LAM_SF, W_SF, K)
    Z, Dl = ph + sf, ph - sf
    print(f" {N:5d}     {lam1(O_general(Dl, Dl, N)):.10f}          {lam1(O_general(Z, Dl, N)):.10f}")

print()
print("  方向自洽性：Lam^Delta <= Lam^Z 逐元？", bool(np.all(Dl <= Z + 0.0)))
print("  ⟹ 单核版 -L1 扣减更小（更不负）⟹ lam1 更大。")

print()
print("=== Stage C: f(0) 归属占那 6.9% 差异的多少？（复核 guidance-lead 的自我订正）===")
print("  (a) 物理归属：f(0) 的 Z 侧归 Lam^Z、配对侧归 Lam^Delta")
print("  (b) 照笔记 (1-delta) 约定直译到双核：两处 f(0) 双双丢弃")
print("  (c) 单核 O[Lam^Delta]")
print("  N      (a)                (b)                (c)")
for N in (400, 800, 1600, 3200):
    K = 2 * N + 2
    ph, sf = einstein(LAM_PH, W_PH, K), einstein(LAM_SF, W_SF, K)
    Z, Dl = ph + sf, ph - sf
    a = lam1(O_general(Z, Dl, N))
    b = lam1(O_general(Z, Dl, N, drop_f0=True))
    c = lam1(O_general(Dl, Dl, N))
    print(f" {N:5d}  {a:.10f}     {b:.10f}     {c:.10f}")
print(f"  a-b = {a-b:+.4e}  占 (a-c) 的 {100*(a-b)/(a-c):.1f}%")
print(f"  b-c = {b-c:+.4e}  占 (a-c) 的 {100*(b-c)/(a-c):.1f}%")
print(f"  解析检验（对角元之差 = u_n^2 (Lam^D(0)-Lam^Z(0)) = -2 u_n^2 lam_sf）：")
print(f"    Lam^D(0)-Lam^Z(0) = {Dl[0]-Z[0]:.6f}   -2*lam_sf = {-2*LAM_SF:.6f}")
print("  ⟹ 只把行和换成 Lam^Z、不动 f(0)，会停在 (b)=0.1227，离 (a) 还差 6.2%，")
print("     并会误判「约定还是没对上」。f(0) 归属是主项（90.5%），不是脚注。")
