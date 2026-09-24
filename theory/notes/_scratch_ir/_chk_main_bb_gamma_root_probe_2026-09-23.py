# -*- coding: utf-8 -*-
"""探针：根盒（整块 (q,σ) 正方形）上的 γ 到底有多大 —— 它决定 KT 的适用范围。

γ 只要是「对根盒成立的 λ₂ 上界」，对任何子盒就都成立（根盒 ⊇ 子盒）⟹ 只算一次。
但若 γ_root 太大（比如 > 0.6），中等大小的成员（λ₁∈(γ,1)）会同时被 KT（前提 μ>γ 不成立）
和 R₀（>1）拒绝 ⟹ 分支定界会在中区无限细分。本探针就是先量这件事，再决定 I 取多少。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import time
import numpy as np
from mpmath import iv

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


BB = _load("_chk_main_box_branch_bound_certificate_2026-09-23")
KT = BB.KT

r, mus, lam_star = 10.0, 0.13, 3.2953
wE_star = KT.einstein_root(lam_star, mus, r, 600)
M2 = lam_star*wE_star*wE_star
lam = 3.29
wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
lam_iv = iv.mpf(lam); M2_iv = iv.mpf(lam)*iv.mpf(wE)*iv.mpf(wE)
print(f"λ={lam} ω̄₂={wE:.9f} N_c={Nc}  M2={M2:.9f}")
print()
print(f"{'I':>5}{'J':>6}{'行尾 B_(I+1)':>16}{'R₀(根盒)':>14}{'γ(根盒)':>14}{'秒':>8}")
for (I, J) in [(25, 250), (40, 300), (60, 400), (90, 500), (140, 700)]:
    t0 = time.time()
    g = BB.G_box(lam, wE, (0.0, 1.0, 0.0, 1.0), I + J + 2)
    r0, gam = BB.schur_bounds(g, lam_iv, M2_iv, I, J)
    tail = float(BB.row_tail_up(lam_iv, M2_iv, I+1).b)
    print(f"{I:>5}{J:>6}{tail:>16.6f}{r0:>14.6f}{gam:>14.6f}{time.time()-t0:>8.1f}")

# ─────────────────────────────────────────────────────────────────────────────
# 两工具的覆盖：R₀<1 要求 λ₁ 小、KT 要求 μ>γ（λ₁ 大）⟹ 必须有重叠，否则中区无限细分。
# 逐盒算 γ（不是根盒 γ）时，γ 随盒子变小而变小 ⟹ 重叠是否成立要实测。
print()
print("逐盒 (R₀, γ) vs 盒内 float64 λ₁ —— 看两工具的适用区间是否重叠（I=60, J=400）")
print(f"  {'盒':<30}{'R₀':>10}{'γ':>10}{'λ₁(中心,float)':>16}{'R₀<1?':>8}{'μ>γ?':>8}{'秒':>7}")
for bx in [(0.30, 0.35, 0.30, 0.35), (0.10, 0.20, 0.70, 0.90), (0.05, 0.10, 0.80, 0.85),
           (0.02, 0.04, 0.90, 0.92), (0.00, 0.05, 0.95, 1.00), (0.20, 0.25, 0.60, 0.65),
           (0.15, 0.20, 0.50, 0.55), (0.30, 0.35, 0.45, 0.50), (0.00, 0.05, 0.00, 0.05),
           (0.95, 1.00, 0.95, 1.00)]:
    t0 = time.time()
    g = BB.G_box(lam, wE, bx, 60 + 400 + 2)
    r0, gam = BB.schur_bounds(g, lam_iv, M2_iv, 60, 400)
    qc, sc = 0.5*(bx[0]+bx[1]), 0.5*(bx[2]+bx[3])
    l1 = KT.lam_k(*BB.QS.mem_qs(lam, wE, qc, sc), mus, Nc, 400)
    print(f"  q[{bx[0]:.2f},{bx[1]:.2f}] σ[{bx[2]:.2f},{bx[3]:.2f}]      "
          f"{r0:>10.4f}{gam:>10.4f}{l1:>16.6f}{'是' if r0 < 1 else '否':>8}"
          f"{'是' if l1 > gam else '否':>8}{time.time()-t0:>7.1f}")
