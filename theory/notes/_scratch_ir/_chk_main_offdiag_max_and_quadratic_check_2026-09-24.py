# -*- coding: utf-8 -*-
"""补测 #2（team-lead 2026-09-24 订正后）：λ=2.00/2.50 的**非对角**子集 max，
以及「亏损/|σ−q|²」是否 O(1)（判「余量是物理的还是抽样密度指纹」）。

背景：team-lead 指出「全体抽样 max == Einstein」这件事自由度为 0——抽样族含对角点
(q=σ) 就恒等于 Einstein，逐位相符不测任何东西。可报的只有：
  (i) 过滤 |q-σ|>0 后的 max，及其到对角线的距离；
  (ii) 亏损 deficit := λ₁(Einstein) − λ₁(offdiag point) 除以 |q-σ|² 是否稳定 O(1)
       （若稳定，说明这个「亏损」只是抽样离对角线多近的指纹，不是 Einstein 的真实领先幅度——
       因为 G_E − G 对 (σ-q) 恰好二次，见 §12.2）。

两部分：
  A. 复用与 §12 GATE B4b 同参数的逐叶抽样（8 点/叶），过滤 |q-σ|>1e-9，报非对角 max。
  B. 系统径向探测：固定 q=0.5，沿 σ=0.5±ε 扫多个 ε（径向逼近对角线），
     直接测 deficit(ε)/ε² 是否稳定——这是判 O(1) 指纹假说的**直接**实验，
     不依赖抽样恰好落在哪。

标级：全部 numerical-evidence（float64）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import numpy as np

_d = os.path.dirname(os.path.abspath(__file__))
def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

BB = _load("_chk_main_box_branch_bound_certificate_2026-09-23")
KT, QS = BB.KT, BB.QS

r, mus, lam_star = 10.0, 0.13, 3.2953
M2 = lam_star*KT.einstein_root(lam_star, mus, r, 600)**2
rng = np.random.default_rng(20260923)

for lam in (2.00, 2.50):
    wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
    lE800 = KT.lam_k([lam], [wE], mus, Nc, 800)
    print("=" * 110)
    print(f"λ={lam}  ω̄₂={wE:.9f}  N_c={Nc}  λ₁(Einstein, N=800)={lE800:.12f}")

    # ── A. 逐叶抽样，过滤 |q-σ|>1e-9，报非对角 max（N=500，同 B4b 协议） ──
    R = BB.certify(lam, wE, mus, Nc, 1.0, verbose=False, maxbox=40000)
    best_off = -1e18; arg_off = None
    for (qlo, qhi, slo, shi, tool, b) in R["leaves"]:
        pts = [(qlo, slo), (qlo, shi), (qhi, slo), (qhi, shi),
               (0.5*(qlo+qhi), 0.5*(slo+shi))]
        pts += [(rng.uniform(qlo, qhi), rng.uniform(slo, shi)) for _ in range(3)]
        for (q, s) in pts:
            if abs(q - s) <= 1e-9:
                continue
            l1 = KT.lam_k(*QS.mem_qs(lam, wE, q, s), mus, Nc, 500)
            if l1 > best_off:
                best_off = l1; arg_off = (q, s)
    q_a, s_a = arg_off
    l1_800 = KT.lam_k(*QS.mem_qs(lam, wE, q_a, s_a), mus, Nc, 800)
    dist = abs(q_a - s_a)
    deficit = lE800 - l1_800
    ratio = deficit / (dist*dist) if dist > 0 else float("nan")
    print(f"  A. 逐叶抽样(N=500)非对角 max={best_off:.9f} @ (q,σ)=({q_a:.9f},{s_a:.9f})"
          f" |q-σ|={dist:.6e}")
    print(f"     重算(N=800): λ₁={l1_800:.9f}  deficit=λ₁(E)-λ₁(pt)={deficit:.6e}"
          f"  deficit/|q-σ|² = {ratio:.4f}")

    # ── B. 系统径向探测：q=0.5 固定，σ=0.5±ε ──
    print("  B. 径向探测 q=0.5, σ=0.5±ε（N=800）：deficit(ε) 与 deficit/ε² 是否稳定 O(1)")
    print(f"     {'eps':>10}{'sigma':>12}{'lam1(pt)':>15}{'deficit':>14}{'deficit/eps^2':>16}")
    for eps in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001):
        s = 0.5 + eps
        l1 = KT.lam_k(*QS.mem_qs(lam, wE, 0.5, s), mus, Nc, 800)
        d = lE800 - l1
        print(f"     {eps:>10.4f}{s:>12.4f}{l1:>15.9f}{d:>14.6e}{d/eps**2:>16.4f}")
    print()
