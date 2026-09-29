# -*- coding: utf-8 -*-
"""
主会话 2026-09-25：新天花板 B(λ, ω̄₁, ω̄₂) 在真实隧穿数据上的紧度（接
`_chk_main_shape_factor_single_atom_and_leavens_2026-09-25.py`）。

B := sup{Tc : ν≥0, λ(ν)≤λ, ⟨ω⟩-面积 A(ν)≤λω̄₁/2, M2(ν)≤λω̄₂²}，μ*=0（⟹ 对一切 μ*≥0 是上界）。
单原子归约 ⟹ Tc_B = ω̄₂/(2π) · max_a min(ρ₁/a, 1/√a)·τ(aλ)，ρ₁=ω̄₁/ω̄₂，τ(x)=t*(x,1)。
B 对 (λ, ω̄₁, ω̄₂) 各自单调不减（约束取 ≤，集合嵌套）⟹ 三者的**上界**代入即给声音上界。

对照四条（同一行、同一输入）：
  r_∞   = Tc / (C_∞ √λ ω̄₂)             仓库头条（λ 无关，analytic）
  r_E   = Tc / (C_E(λ) √λ ω̄₂)          λ 分辨 Einstein（μ*=0，float）
  r_L   = Tc / (0.2309 A),  A=λω̄₁/2    Leavens 1975 面积界（数值局部极大，非严格）
  r_new = Tc / Tc_B                      本文
一致性：r_new ≥ max(r_E, r_L·(0.2309/0.230868)) 必须成立（集合包含）。

数据：AD75 Table I 29 行（页图转录，`_chk_ad75_table1_frak_f_2026-09-20.py`，列 wbar1/wbar2 已三重钉死为频率 K）；
Kihlstrom 1984 Table II 13 行（⟨ω⟩ 的 λ 加权定义一手钉死；⟨ω²⟩ 未定义 ⟹ 另给 Jensen 最坏列）。
方向：隧穿反演是**双边估计** ⟹ 本表只量「紧了多少」，不产生任何 PASS/FAIL 判词。numerical-evidence。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import math
import re
import numpy as np
from multiprocessing import Pool

_d = os.path.dirname(os.path.abspath(__file__))


def _load(fn, name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, fn))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


S = _load("_chk_main_shape_factor_single_atom_and_leavens_2026-09-25.py", "S")
tstar, t1 = S.tstar, S.t1
C_INF = 0.1827262477
LEAV = 0.2309


def _grab(fn, var):
    src = open(os.path.join(_d, fn), encoding="utf-8").read()
    i = src.index(var + " = ["); j = src.index("\n]", i)
    return eval(src[i+len(var)+3:j+2])


T1 = _grab("_chk_ad75_table1_frak_f_2026-09-20.py", "T1")
KH = _grab("_chk_main_kihlstrom1984_nb3ge_tunneling_kappa_2026-09-25.py", "ROWS")


def row(args):
    name, lam, w1, w2, tc = args
    rho = w1/w2
    tE = tstar(lam, [1.], [1.]); tb, a = t1(lam, min(rho, 1.0))
    r_inf = tc/(C_INF*math.sqrt(lam)*w2)
    r_E = tc/(tE*w2/(2*math.pi))
    r_L = tc/(LEAV*lam*w1/2)
    r_new = tc/(tb*w2/(2*math.pi))
    return name, lam, rho, r_inf, r_E, r_L, r_new, a


def main():
    jobs = [(m, lam, w1, w2, tc) for m, wl, w1, w2, lam, wph, mus, tc, eta in T1]
    # Kihlstrom：⟨ω⟩ [meV]、⟨ω²⟩ [meV²] ⟹ 换成 K（×11.6045）；同一单位下比值无关
    MEV2K = 11.604518
    for r in KH:
        lam = float(r[3]); w1 = float(r[6])*MEV2K; w2 = math.sqrt(float(r[8]))*MEV2K
        jobs.append((f"Nb3Ge#{r[0]}({r[1]}%)", lam, w1, w2, float(r[2])))
    with Pool(18) as P: res = P.map(row, jobs)
    print(f"{'材料':22s} {'λ':>5s} {'ρ₁':>6s} {'r_∞':>7s} {'r_E':>7s} {'r_L':>7s} {'r_new':>7s} "
          f"{'r_new/max(r_E,r_L)':>18s}  a*")
    viol = 0; gains = []
    for name, lam, rho, ri, rE, rL, rn, a in res:
        m = max(rE, rL*LEAV/0.230868)
        g = rn/max(rE, rL); gains.append((g, name, rn/rE, rn/ri))
        if rn < m*(1-1e-6): viol += 1
        print(f"{name:22s} {lam:5.2f} {rho:6.3f} {ri:7.4f} {rE:7.4f} {rL:7.4f} {rn:7.4f} {g:18.4f}  {a:.3f}")
    print(f"\n一致性（r_new ≥ max(r_E, r_L) 必成立）违反：{viol}/{len(res)}")
    gains.sort()
    print(f"r_new/max(r_E,r_L)：中位 {np.median([g[0] for g in gains]):.4f}，最大 {gains[-1][0]:.4f}（{gains[-1][1]}）")
    print(f"r_new/r_E（对 λ 分辨 Einstein 的紧化）：中位 {np.median([g[2] for g in gains]):.4f}，"
          f"范围 [{min(g[2] for g in gains):.4f}, {max(g[2] for g in gains):.4f}]")
    print(f"r_new/r_∞（对仓库头条的紧化）：中位 {np.median([g[3] for g in gains]):.4f}，"
          f"范围 [{min(g[3] for g in gains):.4f}, {max(g[3] for g in gains):.4f}]")
    print(f"max r_new = {max(x[6] for x in res):.4f}（{max(res, key=lambda x: x[6])[0]}）；max r_∞ = {max(x[3] for x in res):.4f}")


if __name__ == "__main__":
    main()
