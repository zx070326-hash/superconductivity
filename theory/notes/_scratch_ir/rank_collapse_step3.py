#!/usr/bin/env python3
"""
候选 A · 秩破坏定量理论 —— Step 3：物理翻译

对 P4D_coulomb_uniqueness.md Sec.4 的 11 个物理锚点（取自其文件，逐条复制 lambda, mu*, r,
varpi_bar_c, Nc, mu_tilde=mu*/lambda），从零重算：
  - kappa := 2*mu_tilde*c1^2/gamma  （P4D 文件定义的"控制参数"，c1=<u_c,psi1>, gamma=mu1-mu2(A)）
  - 真实 gap/lambda1(M) = (lambda1(M)-lambda2(M))/lambda1(M)（Step 1 已确认这是唯一物理相关的
    "秩坍缩"量，且与 kappa 不是同一件事——见下方独立发现）。

层级：numerical-evidence（float64 dense eigh，N=600，Nc 最大 153 留足够 buffer）。
"""
import numpy as np
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
sys.path.insert(0, "notes/_scratch_ir")
from rank_collapse_step1 import build_A, u_c_vec, exact_top2

anchors = [
    # lambda, mu_star, r, varpi_c, Nc
    (3.0, 0.13, 10, 0.60, 6),
    (2.0, 0.13, 10, 0.98, 9),
    (2.0, 0.10, 10, 0.98, 9),
    (1.0, 0.13, 10, 2.20, 22),
    (1.0, 0.10, 10, 2.00, 20),
    (0.7, 0.13, 10, 4.56, 45),
    (0.5, 0.13, 10, 13.5, 135),
    (0.5, 0.13, 3, 15.6, 46),
    (0.4, 0.13, 3, 51.2, 153),
    (2.0, 0.13, 3, 3.00, 9),
    (1.0, 0.13, 3, 2.50, 7),
]

N = 600
print(f"{'lambda':>7s} {'mu*':>6s} {'r':>4s} {'varpi_c':>8s} {'Nc':>4s} {'mu~':>7s} "
      f"{'gamma':>9s} {'c1':>8s} {'kappa':>8s} {'l1(M)':>10s} {'l2(M)':>10s} {'gap/l1(M)':>10s}")
rows = []
for (lam, mus, r, wc, Nc) in anchors:
    mut = mus/lam
    A = build_A(N, wc)
    ev, evec = np.linalg.eigh(A)
    mu1, mu2 = ev[-1], ev[-2]
    gamma = mu1-mu2
    psi1 = evec[:,-1]
    if psi1[0] < 0: psi1 = -psi1
    uc = u_c_vec(N, Nc)
    c1 = float(uc@psi1)
    kappa = 2*mut*c1**2/gamma
    l1, l2, _ = exact_top2(A, uc, mut)
    gap = (l1-l2)/l1
    rows.append((lam,mus,r,wc,Nc,mut,gamma,c1,kappa,l1,l2,gap))
    print(f"{lam:7.2f} {mus:6.2f} {r:4d} {wc:8.2f} {Nc:4d} {mut:7.4f} "
          f"{gamma:9.5f} {c1:8.4f} {kappa:8.4f} {l1:10.3e} {l2:10.3e} {gap:10.5f}")

kappas = np.array([r[8] for r in rows])
gaps = np.array([r[11] for r in rows])
print()
print(f"kappa range across 11 anchors: [{kappas.min():.4f}, {kappas.max():.4f}]")
print(f"true gap/l1(M) range across 11 anchors: [{gaps.min():.4f}, {gaps.max():.4f}]")

print()
print("="*78)
print("直接检验：kappa 是否是 gap/l1(M) 坍缩的好代理？（本会话独立复算，见 Step1 的发现）")
print("="*78)
# 对最高 kappa 的锚点做精细 mu~ 扫描，找真实 gap/l1(M) 的最小值，并与 kappa 对照
for (lam,mus,r,wc,Nc) in [anchors[8], anchors[6], anchors[9]]:
    A = build_A(N, wc)
    ev, evec = np.linalg.eigh(A)
    mu1, mu2 = ev[-1], ev[-2]
    gamma = mu1-mu2
    psi1 = evec[:,-1]
    if psi1[0]<0: psi1=-psi1
    uc = u_c_vec(N, Nc)
    c1 = float(uc@psi1)
    mut_anchor = mus/lam
    best = None
    for mut in np.linspace(0.001, min(0.6, 3*mut_anchor), 200):
        l1,l2,_ = exact_top2(A, uc, mut)
        gap = (l1-l2)/l1 if l1>0 else np.nan
        kap = 2*mut*c1**2/gamma
        if best is None or gap < best[2]:
            best = (mut, kap, gap)
    print(f"lambda={lam} mu*={mus} r={r} varpi_c={wc} Nc={Nc}: anchor mu~={mut_anchor:.4f}, "
          f"anchor kappa={2*mut_anchor*c1**2/gamma:.4f}, anchor gap/l1={ (lambda: None)(None) if False else None}")
    l1a,l2a,_ = exact_top2(A, uc, mut_anchor)
    gap_a = (l1a-l2a)/l1a
    print(f"    at anchor: gap/l1={gap_a:.5f}  |  min gap/l1 over mu~ in [0,{min(0.6,3*mut_anchor):.2f}]: "
          f"mu~*={best[0]:.4f} kappa*={best[1]:.4f} gap*={best[2]:.5f}")
