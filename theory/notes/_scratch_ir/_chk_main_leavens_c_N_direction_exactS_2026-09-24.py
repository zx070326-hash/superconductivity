"""从零复核 leavens-cert Phase 1 的两条承重方向断言（不 import 任何既有模块）。

断言 A: lam_area ∝ c^{-2} 递减 ⟹ 要 lam_area 上界就要 c 下界；c = max_w c(w) ⟹ 单点即可。
断言 B: c(N) 随 N 单调**下降** ⟹ 有限 N 的 c 是**高**估 ⟹ 认证下界不声音，需解析尾界引理。

本脚本检验：断言 B 是否是「S_n 用了 N 截断行和」的指纹。
S_n 闭式（我从定义重推）: S_n = G(0) + 2*sum_{j=1..n} G(j)，两条无穷尾恒等相消。
若用闭式，N x N 块 = K^inf 的真主子块 ⟹ Cauchy 交错成立 ⟹ theta_max(N) 单调**上升**、有限 N 是下估。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

def S_exact(w, n_dim):
    k = np.arange(n_dim)
    G = w * w / (w * w + k * k)
    return G[0] + 2.0 * (np.cumsum(G) - G[0])          # G(0) + 2*sum_{j=1..n}G(j)

def S_trunc(w, n_dim):
    n = np.arange(n_dim)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * n_dim + 2)
    G = w * w / (w * w + k * k)
    return (G[np.abs(i - j)] - G[i + j + 1]).sum(axis=1)

def Khat(w, n_dim, exact):
    n = np.arange(n_dim)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * n_dim + 2)
    G = w * w / (w * w + k * k)
    S = S_exact(w, n_dim) if exact else S_trunc(w, n_dim)
    return (G[np.abs(i - j)] + G[i + j + 1]) * np.outer(u, u) - np.diag(u * u * S), u

def theta_max(w, mus, r, n_dim, exact):
    K, u = Khat(w, n_dim, exact)
    Nc = int(np.floor(r * w))
    B = np.eye(n_dim)
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(n_dim) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    return float(eigh(K, B, subset_by_index=[n_dim-1, n_dim-1], eigvals_only=True)[0])

def c_of_w(w, mus, r, n_dim, exact):
    return theta_max(w, mus, r, n_dim, exact) / (np.pi * w)

def c_max(mus, r, n_dim, exact, lo=0.15, hi=6.0, ng=220, nref=60):
    ws = np.linspace(lo, hi, ng)
    cs = np.array([c_of_w(w, mus, r, n_dim, exact) for w in ws])
    i = int(np.argmax(cs))
    a, b = ws[max(i-1,0)], ws[min(i+1,ng-1)]
    gr = 0.5*(np.sqrt(5.0)-1.0)
    x1, x2 = b-gr*(b-a), a+gr*(b-a)
    f = lambda w: c_of_w(w, mus, r, n_dim, exact)
    f1, f2 = f(x1), f(x2)
    for _ in range(nref):
        if f1 < f2: a, x1, f1 = x1, x2, f2; x2 = a+gr*(b-a); f2 = f(x2)
        else:       b, x2, f2 = x2, x1, f1; x1 = b-gr*(b-a); f1 = f(x1)
    wm = 0.5*(a+b)
    return f(wm), wm

print("="*96)
print("GATE X0  S_n 闭式 vs N 截断行和（同一 w、同一 n_dim）")
for w in (0.3, 0.6664, 1.0, 2.2):
    for nd in (200, 800):
        se, st = S_exact(w, nd), S_trunc(w, nd)
        d = np.abs(se - st)
        print(f"  w={w:6.4f} n_dim={nd:4d}  max|S_ex-S_tr|={d.max():.3e} "
              f"(at n={int(np.argmax(d)):4d})  S_ex[0]={se[0]:.6f} S_tr[0]={st[0]:.6f}  "
              f"尾端 S_ex[-1]={se[-1]:.6f} S_tr[-1]={st[-1]:.6f}")

print()
print("="*96)
print("GATE X1  c(N) 的 N 方向：闭式 S vs 截断 S。断言 B 说『下降』。")
for mus, r in ((0.0, 10.0), (0.13, 10.0)):
    for exact in (True, False):
        row = []
        for nd in (100, 200, 400, 800):
            c, wm = c_max(mus, r, nd, exact)
            row.append((nd, c, wm))
        tag = "闭式S " if exact else "截断S "
        s = "  ".join(f"N={nd}:{c:.9f}" for nd, c, _ in row)
        dif = [row[k+1][1]-row[k][1] for k in range(len(row)-1)]
        mono = "单调↑(下估,Cauchy交错)" if all(d > 0 for d in dif) else (
               "单调↓(高估)" if all(d < 0 for d in dif) else "非单调")
        print(f"  mu*={mus:4.2f} {tag} {s}")
        print(f"        增量 {['%+.2e'%d for d in dif]}  ==> {mono}   (w*={row[-1][2]:.6f})")
