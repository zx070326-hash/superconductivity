#!/usr/bin/env python3
"""
候选 A · 秩破坏定量理论 —— Step 1：秩一情形，两能级闭式

独立重建 A = H(ϖ̄) = O[Einstein 核]（对照 P4_main_theorems.md §0 定义、
theory/notes/_scratch_ir/coulomb_true_mechanism.py 的 build_H_mat，逐条自查矩阵元）。
M = A - 2*mu_tilde*|u_c><u_c|，u_c[n]=1/sqrt(2n+1) (n<Nc) 否则 0。

产出：
 (a) Cauchy 交错数值核验（标准定理，只验证不重证）。
 (b) 两能级闭式 gap/lambda1(kappa, rho) 的推导实现 + 与 dense eigh 全网格比对。
 (c) 复现 (0.74, 4, 0.322) -> gap/lambda1 ~= 0.14。
 (d) N=400 vs N=800 收敛核对（部分格点）。

层级：numerical-evidence（float64 dense eigh），闭式推导本身是 analytic（标准秩一微扰理论 /
Cauchy 交错，未引入新数学工具）。
"""
import numpy as np
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")

def build_A(N, w):
    """A = O[Einstein kernel](w)，逐条对照 P4_main_theorems.md §0。"""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2*n + 1)
    k_arr = np.arange(1, N+1)
    f_k = w*w / (k_arr**2 + w*w)          # f(k) = varpi^2/(k^2+varpi^2), k=1..N
    # L1 diagonal: (2n+1)^{-1} * sum_{k=1}^n 2 f(k)
    csum = np.cumsum(2*f_k)               # csum[n-1] = sum_{k=1}^n 2f(k)
    L1diag = np.zeros(N)
    L1diag[1:] = csum[:N-1] / (2*n[1:] + 1)
    idx = np.arange(N)
    absdiff = np.abs(idx[:, None] - idx[None, :])
    summ = idx[:, None] + idx[None, :] + 1
    def f(k):
        return w*w / (k.astype(float)**2 + w*w)
    L2 = f(absdiff) * np.outer(u, u)
    np.fill_diagonal(L2, 0.0)
    L3 = f(summ) * np.outer(u, u)
    A = -np.diag(L1diag) + L2 + L3
    A = 0.5*(A + A.T)
    return A

def u_c_vec(N, Nc):
    n = np.arange(N)
    u = 1.0/np.sqrt(2*n+1)
    u[Nc:] = 0.0
    return u

def secular_two_level(mu1, mu2, c1, c2, mu_tilde):
    """两能级闭式：返回 (lambda1, lambda2_block, kappa, rho)。"""
    gamma = mu1 - mu2
    kappa = 2*mu_tilde*c1**2/gamma
    rho = (c2/c1)**2
    B = 1 + kappa*(1+rho)
    disc = B*B - 4*kappa
    if disc < 0:
        return None
    sq = np.sqrt(disc)
    x_minus = (B - sq)/2.0
    x_plus  = (B + sq)/2.0
    lam1 = mu1 - gamma*x_minus
    lam2_block = mu1 - gamma*x_plus
    return lam1, lam2_block, kappa, rho

def exact_top2(A, u_c, mu_tilde):
    M = A - 2*mu_tilde*np.outer(u_c, u_c)
    ev = np.linalg.eigvalsh(M)
    return ev[-1], ev[-2], ev

def exact_A_top3(A):
    ev = np.linalg.eigvalsh(A)
    return ev[-1], ev[-2], ev[-3]

print("="*78)
print("(a) Cauchy interlacing check (standard theorem, numerical spot-check)")
print("="*78)
N = 400
rng_checks = [(0.74,4,0.10),(1.0,10,0.13),(5.0,20,0.26),(13.5,135,0.26)]
for w, Nc, mut in rng_checks:
    A = build_A(N, w)
    uc = u_c_vec(N, Nc)
    mu1, mu2, mu3 = exact_A_top3(A)
    l1, l2, ev = exact_top2(A, uc, mut)
    ok_a = (mu2 <= l1 <= mu1 + 1e-12)
    ok_b = (mu3 <= l2 <= mu2 + 1e-12)
    print(f"  w={w:6.2f} Nc={Nc:4d} mu~={mut:.3f}  mu2<=l1<=mu1: {ok_a} ({mu2:.6f}<={l1:.6f}<={mu1:.6f})"
          f"   mu3<=l2<=mu2: {ok_b} ({mu3:.6f}<={l2:.6f}<={mu2:.6f})")

print()
print("="*78)
print("(b) Two-level closed form vs exact eigh over grid")
print("="*78)
w_grid = [0.3, 0.5, 0.74, 1.0, 1.5, 2.0, 3.0]
Nc_grid = [2, 4, 8, 16, 32, 40]
mut_grid = [0.02, 0.05, 0.1, 0.2, 0.3, 0.45]

N = 400
results = []
for w in w_grid:
    A = build_A(N, w)
    ev_full, evec_full = np.linalg.eigh(A)
    mu1, mu2 = ev_full[-1], ev_full[-2]
    psi1, psi2 = evec_full[:, -1], evec_full[:, -2]
    # 符号归一化（psi1 >= 0，由 Perron-Frobenius，引理0）
    if psi1[0] < 0: psi1 = -psi1
    for Nc in Nc_grid:
        uc = u_c_vec(N, Nc)
        c1 = float(uc @ psi1)
        c2 = float(uc @ psi2)
        if c1 == 0:
            continue
        for mut in mut_grid:
            out = secular_two_level(mu1, mu2, c1, c2, mut)
            if out is None:
                continue
            lam1_cf, lam2_cf, kappa, rho = out
            l1_ex, l2_ex, ev_full_M = exact_top2(A, uc, mut)
            mu3 = ev_full[-3]
            l2_ex_check = max(lam2_cf, mu3)  # 若两能级块第二值低于mu3，被mu3继承覆盖
            rel_l1 = abs(lam1_cf - l1_ex)/abs(l1_ex)
            gap_cf = 1 - max(lam2_cf, mu3)/lam1_cf
            gap_ex = 1 - l2_ex/l1_ex
            rel_gap = abs(gap_cf - gap_ex)/max(abs(gap_ex), 1e-12)
            results.append((w, Nc, mut, kappa, rho, rel_l1, rel_gap, gap_ex, gap_cf))

results = np.array([r for r in results])
rel_l1_arr = results[:,5].astype(float)
rel_gap_arr = results[:,6].astype(float)
kappa_arr = results[:,3].astype(float)

print(f"grid size: {len(results)} points")
print(f"max relative error in lambda1 (2-level vs exact): {rel_l1_arr.max():.4e}")
print(f"max relative error in gap/lambda1 (2-level vs exact): {rel_gap_arr.max():.4e}")

# breakdown by kappa range
for lo, hi in [(0,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.5),(2.5,10)]:
    mask = (kappa_arr>=lo)&(kappa_arr<hi)
    if mask.sum()==0:
        continue
    print(f"  kappa in [{lo},{hi}): n={mask.sum():3d}  max rel err(gap)={rel_gap_arr[mask].max():.3e}  "
          f"median={np.median(rel_gap_arr[mask]):.3e}")

print()
print("="*78)
print("(c) Reproduce (w=0.74, Nc=4, mu~=0.322) -> gap/lambda1 ~= 0.14")
print("="*78)
w, Nc, mut = 0.74, 4, 0.322
for N in [400, 800]:
    A = build_A(N, w)
    ev_full, evec_full = np.linalg.eigh(A)
    mu1, mu2, mu3 = ev_full[-1], ev_full[-2], ev_full[-3]
    psi1, psi2 = evec_full[:,-1], evec_full[:,-2]
    if psi1[0] < 0: psi1 = -psi1
    uc = u_c_vec(N, Nc)
    c1 = float(uc@psi1); c2 = float(uc@psi2)
    out = secular_two_level(mu1, mu2, c1, c2, mut)
    lam1_cf, lam2_cf, kappa, rho = out
    l1_ex, l2_ex, _ = exact_top2(A, uc, mut)
    gap_ex = 1 - l2_ex/l1_ex
    gap_cf = 1 - max(lam2_cf, mu3)/lam1_cf
    print(f"  N={N}: kappa={kappa:.4f} rho={rho:.4f}  gap/l1 exact={gap_ex:.5f}  2-level cf={gap_cf:.5f}"
          f"  l1_ex={l1_ex:.6e} l2_ex={l2_ex:.6e} mu3={mu3:.6e}")

print()
print("Done.")
