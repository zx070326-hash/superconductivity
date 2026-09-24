# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ②(c) 前置：**试探向量能不能跨盒子复用？**

这决定 ② 的真实成本。Kato–Temple
    λ₁ ≤ μ + res²/(μ−γ),  μ=⟨v,Mv⟩, res²=‖Mv−μv‖², γ≥λ₂
对**任意**单位向量 v 成立 —— v 不必是特征向量，选错只让界变松、**不破坏声音性**。
⟹ 区间化时可以在盒**中心**算一个**固定的浮点 v**（当作精确有理数看待），
   只把 μ、res²、γ 做区间包围。**不需要区间特征值求解器**——那本是 (c) 看起来最难的一块。

本脚本判：固定 v（取自盒中心成员）用到盒内其他成员上，δ 涨多少？
    ⟹ 直接读出「盒子能开多大」⟹ 盒子数 ⟹ ② 的工作量。

坐标用**紧化后的** (x, t)，t := 1/ρ ∈ (0,1]，t=0 是 ②-a 已证的 Einstein 极限边界。

严格性：KT 对任意 v 成立是 `analytic`；本脚本一切数值为 `numerical-evidence`（float64）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import numpy as np
from scipy.linalg import eigh

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")
PL = _load("_chk_main_lam2_peeling_certificate_2026-09-23")


def mem_xt(lam, wE, x, t):
    """t = 1/ρ；t=0 ⟹ ②-a 的极限成员 (xλ, ω̄_E/√x)。"""
    if t <= 0.0:
        return np.array([x*lam]), np.array([wE/np.sqrt(x)])
    return KT.members(lam, wE, x, 1.0/t)


def kt_with_v(L_, W_, mus, Nc, v, Nv, gamma, L=40000):
    Mv = KT.apply_M(L_, W_, mus, Nc, v, Nv, L)
    mu = float(v @ Mv[:Nv])
    rr = Mv.copy(); rr[:Nv] -= mu*v
    res2 = float(rr @ rr) + KT.res_tail_bound(L_, W_, v, Nv, L)
    if mu <= gamma:
        return None, mu
    return mu + res2/(mu - gamma), mu


def vstar(L_, W_, mus, Nc, Nv):
    Msub, _, _ = KT.build(L_, W_, mus, Nc, Nv)
    _, V = eigh(Msub, subset_by_index=[Nv-1, Nv-1])
    v = V[:, 0].ravel()
    return v/np.linalg.norm(v)


def main():
    r, Nv, N_R, N_ref = 10.0, 120, 1200, 1000
    lam, mus = 3.2953, 0.13
    wE = KT.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    print(f"基准 λ={lam} μ*={mus} r={r} ⟹ ω̄_E*={wE:.9f} N_c={Nc}；坐标 (x, t=1/ρ)")
    print()

    print("=" * 108)
    print("GATE V1  固定 v（盒中心）跨盒复用：δ = KT(v_center)/λ₁ − 1 随盒半径 h 怎么涨？")
    print("         γ 也固定取盒内**最大**的 R[w=u]（保守；区间化时就是区间上端）")
    print(f"  {'box center (x,t)':<20}{'h':>7}{'δ(自身 v*)':>14}{'δ(中心 v, 角点)':>18}"
          f"{'max over 4 角':>14}{'方向':>8}")
    centers = [(0.30, 0.30), (0.50, 0.50), (0.70, 0.20), (0.90, 0.05), (0.50, 0.95)]
    for (xc, tc) in centers:
        Lc, Wc = mem_xt(lam, wE, xc, tc)
        vc = vstar(Lc, Wc, mus, Nc, Nv)
        for h in (0.20, 0.10, 0.05, 0.02):
            corners = [(np.clip(xc+sx*h, 1e-3, 1-1e-3), np.clip(tc+st*h, 0.0, 1.0))
                       for sx in (-1, 1) for st in (-1, 1)]
            gam = max(PL.R_of(*mem_xt(lam, wE, cx, ct), N_R, 0, PL.weights("u", N_R))[0]
                      for cx, ct in corners)
            dself = dmax = -1.0; ok = "OK"
            for cx, ct in corners:
                L_, W_ = mem_xt(lam, wE, cx, ct)
                ref = KT.lam_k(L_, W_, mus, Nc, N_ref)
                b_self, _ = kt_with_v(L_, W_, mus, Nc, vstar(L_, W_, mus, Nc, Nv), Nv, gam)
                b_ctr, _ = kt_with_v(L_, W_, mus, Nc, vc, Nv, gam)
                if b_self is not None:
                    dself = max(dself, b_self/ref - 1.0)
                if b_ctr is None:
                    dmax = np.inf; ok = "μ≤γ"
                else:
                    dmax = max(dmax, b_ctr/ref - 1.0)
                    if b_ctr < ref - 1e-12:
                        ok = "**FAIL**"
            print(f"  ({xc:.2f},{tc:.2f}){'':<9}{h:>7.2f}{dself:>14.2e}{dmax:>18.2e}"
                  f"{dmax:>14.2e}{ok:>8}")
        print()

    print("=" * 108)
    print("GATE V2  成本估算：在 h 给定下，δ 是否 ≤ ε=1e−3？⟹ 盒子数 ≈ (1/h)² per (λ,μ*) 格")
    print(f"  {'h':>7}{'盒子数/格':>12}{'max δ 全域':>14}{'max KT 上界':>14}{'≤1+ε?':>9}")
    for h in (0.25, 0.125, 0.0625):
        nb = int(np.ceil(1.0/(2*h)))**2
        dmax = -1.0; ubmax = -1.0; bad = 0
        xs = np.arange(h, 1.0, 2*h); ts = np.arange(h, 1.0+1e-9, 2*h)
        for xc in xs:
            for tc in ts:
                Lc, Wc = mem_xt(lam, wE, xc, tc)
                vc = vstar(Lc, Wc, mus, Nc, Nv)
                corners = [(np.clip(xc+sx*h, 1e-3, 1-1e-3), np.clip(tc+st*h, 0.0, 1.0))
                           for sx in (-1, 1) for st in (-1, 1)]
                gam = max(PL.R_of(*mem_xt(lam, wE, cx, ct), N_R, 0, PL.weights("u", N_R))[0]
                          for cx, ct in corners)
                for cx, ct in corners:
                    L_, W_ = mem_xt(lam, wE, cx, ct)
                    ref = KT.lam_k(L_, W_, mus, Nc, N_ref)
                    b, _ = kt_with_v(L_, W_, mus, Nc, vc, Nv, gam)
                    if b is None:
                        # 远端小 λ₁ 成员：交给工具 1（粗界 R₀），这里只记不算失败
                        continue
                    dmax = max(dmax, b/ref - 1.0); ubmax = max(ubmax, b)
                    bad += (b < ref - 1e-12)
        print(f"  {h:>7.4f}{nb:>12}{dmax:>14.2e}{ubmax:>14.9f}"
              f"{('YES' if ubmax <= 1.0+1e-3 else 'no'):>9}   方向违反={bad}")
    print()
    print("  读法：δ(中心 v) 若与 δ(自身 v*) 同量级，则**固定 v 可复用** ⟹ (c) 不需要区间特征值求解器，")
    print("        只需把 μ、res²、γ 三个显式表达式做区间包围。盒子数 ≈ (1/2h)² per (λ,μ*) 格。")


if __name__ == "__main__":
    main()
