#!/usr/bin/env python3
"""
候选 A · 秩破坏定量理论 —— Step 2：一般负核，普适性检验

M = O[phi_+] - b*O[phi_-] =: A - b*B。定义 kappa_eff := b*<psi1|B|psi1> / gamma，
gamma = mu1(A)-mu2(A)，psi1 = A 的顶本征向量（未受扰动）。
检验：不同核族的 (kappa_eff, gap/lambda1(M)) 是否落在同一条坍缩曲线上。
gap/lambda1(M) 用 M 自身的顶两个本征值算（Step 1 已确认这是唯一物理相关的定义，
不是 P4D_coulomb_uniqueness.md 里 kappa 控制的"lambda1(M) 与 mu2(A) 的边界穿越"那个不同的量）。

层级：numerical-evidence，float64 dense eigh，N=500（部分点 N=1000 复核）。
"""
import numpy as np
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")

def build_O(N, fvals):
    """通用 O[f] 算子。fvals(karray) 需对 k=1..2N 可求值（向量化）。"""
    n = np.arange(N)
    u = 1.0/np.sqrt(2*n+1)
    k_diag = np.arange(1, N+1)
    f_diag_k = fvals(k_diag)
    csum = np.cumsum(2*f_diag_k)
    L1diag = np.zeros(N)
    L1diag[1:] = csum[:N-1] / (2*n[1:]+1)
    idx = np.arange(N)
    absdiff = np.abs(idx[:,None]-idx[None,:]).astype(float)
    summ = (idx[:,None]+idx[None,:]+1).astype(float)
    L2 = fvals(absdiff) * np.outer(u,u)
    np.fill_diagonal(L2, 0.0)
    L3 = fvals(summ) * np.outer(u,u)
    O = -np.diag(L1diag) + L2 + L3
    return 0.5*(O+O.T)

def top2_gap(M):
    ev = np.linalg.eigvalsh(M)
    l1, l2 = ev[-1], ev[-2]
    return l1, l2, (l1-l2)/l1 if l1 != 0 else np.nan

def family_scan(name, phi_plus, phi_minus, b_values, N=500):
    A = build_O(N, phi_plus)
    evA, evecA = np.linalg.eigh(A)
    mu1, mu2 = evA[-1], evA[-2]
    gamma = mu1 - mu2
    psi1 = evecA[:,-1]
    if psi1[0] < 0: psi1 = -psi1
    B = build_O(N, phi_minus)
    Bpsi1 = B @ psi1
    overlap = float(psi1 @ Bpsi1)   # <psi1|B|psi1>
    rows = []
    for b in b_values:
        M = A - b*B
        l1, l2, gap = top2_gap(M)
        kappa_eff = b*overlap/gamma
        rows.append((name, b, kappa_eff, gap, l1, l2))
    return rows, gamma, overlap

k = np.arange(1, 2001, dtype=float)  # 足够覆盖 N=500 下 |n-m|, n+m+1 <= ~1000

def lorentz(k, x):
    return x*x/(k*k+x*x) if np.isscalar(x) else x*x/(k**2+x*x)

all_rows = []

# --- 族 1: 墨西哥帽 phi = 1/(1+k^2) - b/(1+(k/c)^2), c in {2,4,8} ---
for c in [2,4,8]:
    def phi_plus(kk):
        return 1.0/(1.0+kk**2)
    def phi_minus(kk, c=c):
        return 1.0/(1.0+(kk/c)**2)
    bs = np.linspace(0.05, 1.2, 14)
    rows, gamma, ov = family_scan(f"mexhat_c{c}", phi_plus, phi_minus, bs)
    all_rows += rows
    print(f"mexhat c={c}: gamma={gamma:.6f} <psi1|B|psi1>={ov:.6f}")

# --- 族 2: 振荡核，phi = phi+ - b*phi- 中 phi+/phi- 取自 cos(alpha k)/(1+k^2) 的正负部 ---
for alpha in [0.5, 1.0]:
    def phi_plus(kk, alpha=alpha):
        g = np.cos(alpha*kk)/(1.0+kk**2)
        return np.maximum(g, 0.0)
    def phi_minus(kk, alpha=alpha):
        g = np.cos(alpha*kk)/(1.0+kk**2)
        return np.maximum(-g, 0.0)
    bs = np.linspace(0.1, 2.0, 14)
    rows, gamma, ov = family_scan(f"osc_alpha{alpha}", phi_plus, phi_minus, bs)
    all_rows += rows
    print(f"osc alpha={alpha}: gamma={gamma:.6f} <psi1|B|psi1>={ov:.6f}")

# --- 族 3: 高斯负肩，phi = 1/(1+k^2) - b*exp(-((k-k0)/sigma)^2) ---
for (k0, sigma) in [(5.0, 3.0), (15.0, 6.0)]:
    def phi_plus(kk):
        return 1.0/(1.0+kk**2)
    def phi_minus(kk, k0=k0, sigma=sigma):
        return np.exp(-((kk-k0)/sigma)**2)
    bs = np.linspace(0.1, 3.0, 14)
    rows, gamma, ov = family_scan(f"gauss_k0{k0}_s{sigma}", phi_plus, phi_minus, bs)
    all_rows += rows
    print(f"gauss k0={k0} sigma={sigma}: gamma={gamma:.6f} <psi1|B|psi1>={ov:.6f}")

# --- 族 4（自设计）: 幂律尾负部，phi = 1/(1+k^2) - b/(1+k)^p ---
for p in [1.2, 1.6]:
    def phi_plus(kk):
        return 1.0/(1.0+kk**2)
    def phi_minus(kk, p=p):
        return 1.0/(1.0+kk)**p
    bs = np.linspace(0.02, 0.6, 14)
    rows, gamma, ov = family_scan(f"powertail_p{p}", phi_plus, phi_minus, bs)
    all_rows += rows
    print(f"powertail p={p}: gamma={gamma:.6f} <psi1|B|psi1>={ov:.6f}")

# --- 族 5（自设计）: 局域洛伦兹凹陷（notch），phi = 1/(1+k^2) - b/(1+((k-k0)/w)^2) ---
for (k0, w) in [(3.0, 1.5), (10.0, 2.0)]:
    def phi_plus(kk):
        return 1.0/(1.0+kk**2)
    def phi_minus(kk, k0=k0, w=w):
        return 1.0/(1.0+((kk-k0)/w)**2)
    bs = np.linspace(0.05, 1.5, 14)
    rows, gamma, ov = family_scan(f"notch_k0{k0}_w{w}", phi_plus, phi_minus, bs)
    all_rows += rows
    print(f"notch k0={k0} w={w}: gamma={gamma:.6f} <psi1|B|psi1>={ov:.6f}")

print()
print("="*90)
print(f"{'family':22s} {'b':>6s} {'kappa_eff':>10s} {'gap/l1':>10s}")
print("="*90)
for r in all_rows:
    name,b,kap,gap,l1,l2 = r
    print(f"{name:22s} {b:6.3f} {kap:10.4f} {gap:10.5f}")

# ---- 普适性检验：在共同 kappa_eff 网格上插值各族，计算跨族离散度 ----
import numpy as np
families = sorted(set(r[0] for r in all_rows))
common_kappa = np.array([0.05,0.1,0.2,0.3,0.5,0.7,1.0,1.3,1.6,2.0,2.5,3.0])
interp_table = {}
for fam in families:
    ks = np.array([r[2] for r in all_rows if r[0]==fam])
    gs = np.array([r[3] for r in all_rows if r[0]==fam])
    order = np.argsort(ks)
    ks, gs = ks[order], gs[order]
    vals = np.full_like(common_kappa, np.nan, dtype=float)
    for i,ck in enumerate(common_kappa):
        if ks.min() <= ck <= ks.max():
            vals[i] = np.interp(ck, ks, gs)
    interp_table[fam] = vals

print()
print("="*90)
print("普适性表：跨族在同一 kappa_eff 处的 gap/l1（NaN=该族未覆盖该 kappa_eff）")
print("="*90)
header = "kappa_eff " + " ".join(f"{f:>16s}" for f in families)
print(header)
for i,ck in enumerate(common_kappa):
    row = f"{ck:9.3f} " + " ".join(f"{interp_table[f][i]:16.4f}" if not np.isnan(interp_table[f][i]) else f"{'--':>16s}" for f in families)
    print(row)

print()
max_spread = 0.0
for i,ck in enumerate(common_kappa):
    vals = [interp_table[f][i] for f in families if not np.isnan(interp_table[f][i])]
    if len(vals) >= 2:
        spread = max(vals)-min(vals)
        max_spread = max(max_spread, spread)
        print(f"kappa_eff={ck:.3f}: n_families={len(vals)} spread(max-min gap/l1)={spread:.4f}  vals={[f'{v:.3f}' for v in vals]}")

print()
print(f"MAX cross-family spread in gap/l1 at matched kappa_eff: {max_spread:.4f}")
