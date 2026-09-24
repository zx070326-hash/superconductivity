# -*- coding: utf-8 -*-
"""GATE Z — `mu_lo` 到底 fail-open 没有：走**支撑链**，不再拿它和任何「我算的 sup」比。

我此前的判词「`kt_box` 的 `mu_lo` 在 4/7 个下侧端点上高于真 sup ⟹ fail-open」建立在
`mu_lo` vs 我自己扫出来的 sup 之上——而那个 sup 只覆盖两原子族，不覆盖 𝒞̄ 的全部极点片。
那是**同一个对象错配**的第 N 次。本脚本改判**声音性本身**，逐箭头：

  (1) `mu_lo = float(mu.a)`，`mu = ⟨v,M_Nv v⟩/⟨v,v⟩` 的区间包围下端  ⟹ 应 ≤ 同参数下的
      float64 Rayleigh 商 `R_float`；
  (2) `R_float ≤ λ₁(P_Nv)`（变分原理，同一个 Nv×Nv 矩阵）；
  (3) `λ₁(P_Nv) ≤ λ₁(P_600)`（Cauchy 交错，S_n 用闭式 ⟹ 真主子块）。

三条全过 ⟹ `mu_lo` 是 `λ₁(Einstein at that w)` 的合法下界，**不 fail-open**，
我之前的缺陷判词必须撤回。任一条失败 ⟹ 缺陷落在那一步，且位置被钉死。

参数**逐字抄自** `_chk_mustar020_certified_2026-09-24.py` 的 lower_cert()：
q=σ=0.5、Nc=floor(R*w)、Nv=40、Lr=Nv+pad、gamma_up=-1e9。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import math
import numpy as np
from mpmath import iv

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


BB = _load("_chk_main_box_branch_bound_certificate_2026-09-23")
KT = BB.KT
QS = _load("_chk_main_qsigma_coords_and_deficit_identity_2026-09-23")

R = 10.0
NV = 80                        # 交付脚本 lower_cert(Nv=80, pad=200)——**逐字核对过，不是我猜的**
PAD = 200

# (mus, lam, w_lo, 交付日志里印出的 mu_lo)
ROWS = [
    (0.20, 2.0,    1.100033, 1.021190594),
    (0.20, 3.2953, 0.731708, 1.020081705),
    (0.20, 4.0,    0.640491, 1.021883712),
    (0.20, 6.0,    0.495741, 1.025180308),
    (0.13, 2.0,    1.116872, None),
    (0.13, 3.2953, 0.742145, None),
    (0.13, 4.0,    0.649667, None),
]

print("GATE Z — mu_lo 支撑链（全部在同一 (lam,w,mus,Nc) 上）")
print(f"{'mus':>5} {'lam':>7} {'w_lo':>10} {'Nc':>3} | {'mu_lo':>14} {'R_float':>14} "
      f"{'lam1(P_Nv)':>14} {'lam1(P_600)':>14} | {'(1)':>4}{'(2)':>4}{'(3)':>4}")
bad = 0
for mus, lam, w, mlo_ref in ROWS:
    Nc = int(math.floor(R * w))
    Lr = NV + PAD
    lam_iv = iv.mpf(lam); M2_iv = iv.mpf(lam) * iv.mpf(w) * iv.mpf(w)

    v = BB.vstar(lam, w, 0.5, 0.5, mus, Nc, NV)
    g = BB.G_box(lam, w, (0.5, 0.5, 0.5, 0.5), Lr + NV + 2)
    _, mlo, mhi, r2 = BB.kt_box(g, lam_iv, M2_iv, mus, Nc, v, NV, Lr, -1e9)

    lams_a, ws_a = QS.mem_qs(lam, w, 0.5, 0.5)                 # = Einstein
    Msub, _, _ = KT.build(lams_a, ws_a, mus, Nc, NV)
    Rf = float(v @ (Msub @ v) / (v @ v))
    l_nv = KT.lam_k(lams_a, ws_a, mus, Nc, NV, k=1)
    l_600 = KT.lam_k(lams_a, ws_a, mus, Nc, 600, k=1)
    l_big = KT.lam_k(lams_a, ws_a, mus, Nc, 2400, k=1)

    a1 = mlo <= Rf
    a2 = Rf <= l_nv
    a3 = l_nv <= l_600
    bad += (not a1) + (not a2) + (not a3)
    if mlo_ref is not None and abs(mlo - mlo_ref) > 5e-9:
        print(f"    !! mu_lo 复现失败：{mlo:.9f} vs 日志 {mlo_ref:.9f}")
    print(f"{mus:5.2f} {lam:7.4f} {w:10.6f} {Nc:3d} | {mlo:14.9f} {Rf:14.9f} "
          f"{l_nv:14.9f} {l_600:14.9f} | {'ok' if a1 else 'XX':>4}{'ok' if a2 else 'XX':>4}"
          f"{'ok' if a3 else 'XX':>4}   d1={Rf-mlo:+.2e} d2={l_nv-Rf:+.2e} d3={l_600-l_nv:+.2e}"
          f"  lam1(P_2400)={l_big:.9f}  (P2400-mu_lo)={l_big-mlo:+.2e}")

print(f"\n违反总数 = {bad}")
print("0 ⟹ mu_lo ≤ λ₁(同参数 Einstein) 逐箭头成立 ⟹ **fail-open 判词撤回**。")
