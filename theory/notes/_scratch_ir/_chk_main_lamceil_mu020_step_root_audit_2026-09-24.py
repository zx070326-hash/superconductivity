# -*- coding: utf-8 -*-
"""GATE L2 — μ*=0.20 的 λ_ceil 到底是 5.7975 还是 5.8589：N_c 台阶造出两个根？

我的 GATE L 二分得 5.7975，头条表反解出 5.8589（差 1.06%）。两者的 w 分居
`w=0.5`（N_c: 5↔4）两侧 —— 正是仓库记过的「库仑饱和造赝根 / λ₁ 对 N_c 非单调」的形状。
判词依赖**定义**：物理陈述是「λ ≤ λ_ceil ⟹ 达不到 300K」，要求 λ₁<1 对**所有** λ'≤λ_ceil
⟹ `λ_ceil = inf{λ : λ₁(λ) ≥ 1}` = **最小**根，不是最大根。
本脚本扫 λ₁(λ)，印 N_c、变号次数与每个根，对四个 μ* 都做（看是否只有 0.20 有这毛病）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

R, M2, N = 10.0, 1.453151466, 400


def lam1(lam, mus):
    w = np.sqrt(M2/lam); Nc = int(np.floor(R*w)); k = np.arange(2*N + 2)
    G = w*w/(w*w + k*k); n = np.arange(N); u = 1.0/np.sqrt(2.0*n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    S = G[0] + 2.0*(np.cumsum(G[:N]) - G[0])
    K = (G[np.abs(i-j)] + G[i+j+1])*np.outer(u, u) - np.diag(u*u*S)
    uc = np.where(n < Nc, u, 0.0)
    B = np.eye(N) + 2.0*mus*np.outer(uc, uc)
    return lam*float(eigh(K, B, eigvals_only=True)[-1]), Nc


for mus, lo, hi in ((0.00, 1.1, 1.6), (0.10, 2.3, 3.0), (0.13, 2.9, 3.7), (0.20, 5.2, 6.6)):
    ls = np.linspace(lo, hi, 1401)
    vals = []; ncs = []
    for L in ls:
        v, nc = lam1(L, mus); vals.append(v - 1.0); ncs.append(nc)
    vals = np.array(vals); ncs = np.array(ncs)
    sgn = np.where(np.diff(np.sign(vals)) != 0)[0]
    steps = np.where(np.diff(ncs) != 0)[0]
    roots = []
    for s in sgn:
        a, b = ls[s], ls[s + 1]; va = vals[s]
        for _ in range(60):
            m = 0.5*(a + b)
            if (lam1(m, mus)[0] - 1.0)*va > 0:
                a = m
            else:
                b = m
        roots.append(0.5*(a + b))
    down = int((np.diff(vals) < 0).sum())
    print(f"μ*={mus:.2f}  变号 {len(sgn)} 次  根={['%.6f' % r for r in roots]}  "
          f"N_c 台阶于 λ≈{['%.4f' % ls[s] for s in steps]} (N_c {ncs[0]}→{ncs[-1]})  "
          f"下降步 {down}/{len(vals)-1}")
    if roots:
        r0 = min(roots); w0 = np.sqrt(M2/r0)
        print(f"        ⟹ inf-根 = {r0:.6f}（w={w0:.6f}, N_c={int(np.floor(R*w0))}）"
              f"；max-根 = {max(roots):.6f}")
