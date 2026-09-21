#!/usr/bin/env python3
"""
Step 2 续：单参数 kappa_eff 坍缩失败之后——2x2 投影块（{psi1,psi2} 有效哈密顿量）
能否给出跨核族普适的 gap/lambda1 预测？

M_block(b) := [[mu1 - b*B11, -b*B12],
               [-b*B12,      mu2 - b*B22]]
B11=<psi1|B|psi1>, B22=<psi2|B|psi2>, B12=<psi1|B|psi2>（A 的顶两个未受扰动本征向量）。
用 M_block 的两个本征值预测 gap/lambda1，和 O[phi_+ - b phi_-] 的真实 dense eigh 比较。

若这个 2x2 块（3 个数 B11,B22,B12，而非 1 个数 kappa_eff）能跨核族给出准确预测，
说明"坍缩曲线"本质上是二维的（(kappa1,kappa2,offdiag) 或等价的 (kappa_eff,rho_eff,mix_eff)），
不是一维的。
"""
import numpy as np
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, "notes/_scratch_ir")
from rank_collapse_step2 import build_O, top2_gap

def family_block_test(name, phi_plus, phi_minus, b_values, N=500):
    A = build_O(N, phi_plus)
    evA, evecA = np.linalg.eigh(A)
    mu1, mu2 = evA[-1], evA[-2]
    gamma = mu1-mu2
    psi1, psi2 = evecA[:,-1], evecA[:,-2]
    if psi1[0] < 0: psi1 = -psi1
    B = build_O(N, phi_minus)
    B11 = float(psi1@B@psi1); B22 = float(psi2@B@psi2); B12 = float(psi1@B@psi2)
    rows = []
    for b in b_values:
        Mb = np.array([[mu1-b*B11, -b*B12],[-b*B12, mu2-b*B22]])
        evb = np.linalg.eigvalsh(Mb)
        l1_blk, l2_blk = evb[-1], evb[-2]
        gap_blk = (l1_blk-l2_blk)/l1_blk if l1_blk!=0 else np.nan

        M = A - b*B
        l1_ex, l2_ex, gap_ex = top2_gap(M)

        kappa_eff = b*B11/gamma
        rho_eff = B22/B11 if B11!=0 else np.nan
        rows.append((name,b,kappa_eff,rho_eff,gap_ex,gap_blk,l1_ex,l1_blk))
    return rows

families_defs = []

for c in [2,4,8]:
    def pp(kk): return 1.0/(1.0+kk**2)
    def pm(kk,c=c): return 1.0/(1.0+(kk/c)**2)
    families_defs.append((f"mexhat_c{c}", pp, pm, np.linspace(0.05,1.2,10)))

for alpha in [0.5,1.0]:
    def pp(kk,alpha=alpha):
        g=np.cos(alpha*kk)/(1.0+kk**2); return np.maximum(g,0.0)
    def pm(kk,alpha=alpha):
        g=np.cos(alpha*kk)/(1.0+kk**2); return np.maximum(-g,0.0)
    families_defs.append((f"osc_alpha{alpha}", pp, pm, np.linspace(0.1,2.0,10)))

for (k0,sigma) in [(5.0,3.0),(15.0,6.0)]:
    def pp(kk): return 1.0/(1.0+kk**2)
    def pm(kk,k0=k0,sigma=sigma): return np.exp(-((kk-k0)/sigma)**2)
    families_defs.append((f"gauss_k0{k0}_s{sigma}", pp, pm, np.linspace(0.1,3.0,10)))

for p in [1.2,1.6]:
    def pp(kk): return 1.0/(1.0+kk**2)
    def pm(kk,p=p): return 1.0/(1.0+kk)**p
    families_defs.append((f"powertail_p{p}", pp, pm, np.linspace(0.02,0.6,10)))

for (k0,w) in [(3.0,1.5),(10.0,2.0)]:
    def pp(kk): return 1.0/(1.0+kk**2)
    def pm(kk,k0=k0,w=w): return 1.0/(1.0+((kk-k0)/w)**2)
    families_defs.append((f"notch_k0{k0}_w{w}", pp, pm, np.linspace(0.05,1.5,10)))

all_rows=[]
for name, pp, pm, bs in families_defs:
    rows = family_block_test(name, pp, pm, bs)
    all_rows += rows

print(f"{'family':22s} {'b':>6s} {'kappa_eff':>9s} {'rho_eff':>8s} {'gap_exact':>10s} {'gap_2x2blk':>11s} {'abserr':>8s}")
abserrs=[]
for r in all_rows:
    name,b,kap,rho,gap_ex,gap_blk,l1e,l1b = r
    err = abs(gap_ex-gap_blk)
    abserrs.append(err)
    print(f"{name:22s} {b:6.3f} {kap:9.4f} {rho:8.4f} {gap_ex:10.5f} {gap_blk:11.5f} {err:8.4f}")

abserrs=np.array(abserrs)
print()
print(f"2x2 block prediction vs exact dense eigh: max |err|={abserrs.max():.4f}  median={np.median(abserrs):.4f}  "
      f"n={len(abserrs)}  frac(err<0.03)={(abserrs<0.03).mean():.3f}")
