# -*- coding: utf-8 -*-
"""GATE B3b（重跑，修了限额记数缺陷）+ GATE B5（修正本身的复核）。

缺陷：`certify` 原来用 `nleaf + len(stack) > maxbox` 当限额。DFS 的栈深 ≤ maxdepth
⟹ 在**永远收不了**的投毒场景里这个量到不了限额，程序会把深度 maxdepth 的整棵树走完
（实测 λ=3.36 跑了 19 分钟还没停）。改为数**已处理盒数** `nproc`。

方向：该缺陷**只影响什么时候放弃**，不影响任何一次接受判定 ⟹ 三个 PASS 跑（自然清空栈、
0 失败盒）的盒子数不受影响。但按「修正本身要独立复核」，GATE B5 用修好的代码重跑 λ=3.20，
盒子数必须**逐位复现** 161。
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
KT = BB.KT
r, mus, lam_star = 10.0, 0.13, 3.2953
M2 = lam_star*KT.einstein_root(lam_star, mus, r, 600)**2

print("GATE B3b  投毒：λ 抬到 Einstein 根之上 ⟹ 真 sup>1 ⟹ 分支定界必须报失败盒", flush=True)
for lam in [3.31, 3.36]:
    wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
    lf = KT.lam_k([lam], [wE], mus, Nc, 800)
    R = BB.certify(lam, wE, mus, Nc, 1.0, maxbox=400, maxdepth=14, verbose=False)
    kinds = sorted({f[0] for f in R["fails"]})
    print(f"  λ={lam}  λ₁(E,float)={lf:.9f}  处理 {R['nproc']} 盒，收 {R['nleaf']}，"
          f"失败 {len(R['fails'])}{kinds}，{R['secs']:.0f}s  ⟹ "
          f"{'未通过（正确）' if R['fails'] else '**通过 = fail-open 缺陷**'}", flush=True)

print()
print("GATE B5  修正复核：用修好限额记数的代码重跑 λ=3.20，盒子数必须仍是 161", flush=True)
lam = 3.20
wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
R = BB.certify(lam, wE, mus, Nc, 1.0, maxbox=40000, verbose=False)
print(f"  λ=3.20：处理 {R['nproc']} 盒，叶 {R['nleaf']}（R₀ {R['n_R0']} / KT {R['n_KT']}），"
      f"深度 {R['maxdepth']}，最劣 {R['worst']:.9f}，失败 {len(R['fails'])}，{R['secs']:.0f}s"
      f"  ⟹ {'复现（161）' if R['nleaf'] == 161 else '**与 161 不符**'}")
