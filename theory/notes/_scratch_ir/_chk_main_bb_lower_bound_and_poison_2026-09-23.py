# -*- coding: utf-8 -*-
"""GATE B3：**另一侧**的认证 —— λ₁ 的认证**下**界，给出 λ 阈值的双侧夹逼 + 负对照。

分支定界只会给上界（sup λ₁ ≤ b）。要证明这套工具**不是 fail-open**（不是「什么都能过」），
需要一个方向相反的认证量：Rayleigh 商

    μ(v) = ⟨v, M v⟩ / ⟨v, v⟩ ≤ λ₁(M)     （变分原理，对任意有限支撑 v 成立）

在 iv 下算出的 **μ 的下端 μ_lo 就是 λ₁ 的认证下界**。μ_lo > 1 ⟹ 该成员确实破天花板 ⟹ FAIL。
⚠ 方向：这里的「下估」支持的是 FAIL，不违反单边界方向规则——因为它下估的是 **λ₁ 本身**，
   不是某个上界。（下估一个**上界**只支持 PASS；下估**被界的量**只支持 FAIL。两件事。）

于是：λ₁ 的临界 λ 被夹在 [分支定界仍 PASS 的最大 λ, μ_lo>1 的最小 λ] 之间，两端都是 certificate。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import numpy as np
from mpmath import iv

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


BB = _load("_chk_main_box_branch_bound_certificate_2026-09-23")
KT, QS = BB.KT, BB.QS

r, mus, lam_star = 10.0, 0.13, 3.2953
M2 = lam_star*KT.einstein_root(lam_star, mus, r, 600)**2
print(f"M2 = {M2:.9f}   μ*={mus}   r={r}   （Einstein 根的 float64 定位 λ*≈{lam_star}）")
print(f"mpmath.iv dps={iv.dps}")
print()
print("GATE B3a  Einstein 成员的认证 Rayleigh 下界 μ_lo（Nv=80）", flush=True)
print(f"  {'λ':>7}{'ω̄₂':>14}{'N_c':>5}{'λ₁(float,N=800)':>18}{'μ_lo（认证下界）':>20}"
      f"{'μ_lo>1?':>10}")
rows = []
for lam in [3.29, 3.2953, 3.31, 3.33, 3.36]:
    wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
    Nv = 80
    Lr = Nv + 200
    lf = KT.lam_k([lam], [wE], mus, Nc, 800)
    v = BB.vstar(lam, wE, 0.5, 0.5, mus, Nc, Nv)          # 对角 q=σ ⟹ Einstein
    g = BB.G_box(lam, wE, (0.5, 0.5, 0.5, 0.5), Lr + Nv + 2)
    lam_iv = iv.mpf(lam); M2_iv = iv.mpf(lam)*iv.mpf(wE)*iv.mpf(wE)
    _, mlo, mhi, _ = BB.kt_box(g, lam_iv, M2_iv, mus, Nc, v, Nv, Lr, -1e9)
    rows.append((lam, mlo))
    print(f"  {lam:>7.4f}{wE:>14.9f}{Nc:>5}{lf:>18.9f}{mlo:>20.9f}"
          f"{('是 ⟹ FAIL' if mlo > 1 else '否'):>10}", flush=True)
print()
print("  ⟹ Einstein 成员在上表最小的「μ_lo>1」处已被**认证**破天花板；")
print("     分支定界在更小的 λ 上给出 sup<1 ⟹ 临界 λ 被双侧夹住。")
print()

# ─────────────────────────────────────────────────────────────────────────────
print("GATE B3b  投毒：把 λ 抬到 Einstein 根之上，分支定界必须**不通过**（限额内）")
for lam in [3.36]:
    wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
    lf = KT.lam_k([lam], [wE], mus, Nc, 800)
    print(f"  λ={lam}  λ₁(Einstein,float)={lf:.9f} > 1 ⟹ 真 sup>1，任何细分都不该给出 PASS", flush=True)
    R = BB.certify(lam, wE, mus, Nc, 1.0, maxbox=150, verbose=False)
    print(f"    限额 150 盒内：叶 {R['nleaf']}，失败盒 {len(R['fails'])}"
          f"（{R['fails'][0][0] if R['fails'] else '—'}），{R['secs']:.0f}s"
          f"  ⟹ {'未通过（正确）' if R['fails'] else '**通过 = fail-open 缺陷**'}")
