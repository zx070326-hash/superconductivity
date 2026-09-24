# -*- coding: utf-8 -*-
"""GATE B4：分支定界的两条**未被 B2 测到的箭头**。

B2 只说「每个叶子上的认证上界 < 1」。整条判词 `sup_{ν∈𝒞̄} λ₁ < 1` 还要两件事：

  箭头①（覆盖）：叶子的并 ∪ 其镜像 ⊇ 整块 (q,σ) 正方形。
      对称裁剪（跳过 σ<q）与「中点二分生成的盒族无缝铺满」都在这条箭头上，B2 一个都没测。
  箭头②（逐叶方向）：每个叶子的认证上界 ≥ 该叶内**任意**成员的 λ₁。
      B1 只对 5 个手选盒子测过；这里对**每一个**叶子测（角点 + 随机内点）。

⚠ 箭头② 的参照量是 float64 截断本征值（`numerical-evidence`，且对真算子是**下估**）
   ⟹ 它只能抓「上界比真值还小」的粗缺陷，不能替代区间算术。判词仍由 B2 的区间量给出。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import sys
import time
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
lam = float(sys.argv[1]) if len(sys.argv) > 1 else 3.29
wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
print(f"λ={lam}  ω̄₂={wE:.9f}  N_c={Nc}  M2={M2:.9f}", flush=True)
R = BB.certify(lam, wE, mus, Nc, 1.0, verbose=False, maxbox=40000)
Lv = R["leaves"]
print(f"分支定界：叶 {R['nleaf']}（R₀ {R['n_R0']} / KT {R['n_KT']}），跳过 {R['nskip']}，"
      f"最劣上界 {R['worst']:.9f}，{R['secs']:.0f}s，失败 {len(R['fails'])}", flush=True)
A = np.array([[l[0], l[1], l[2], l[3]] for l in Lv])
B = np.array([l[5] for l in Lv])

# ── 箭头①：覆盖 ────────────────────────────────────────────────────────────────
print()
print("GATE B4a  覆盖：随机 (q,σ) 必须落在某个叶子里，或其镜像 (1−q,1−σ) 落在某个叶子里")
rng = np.random.default_rng(20260923)
npt = 40000
Q = rng.uniform(0, 1, npt); S = rng.uniform(0, 1, npt)


def inside(q, s):
    return np.any((A[:, 0] <= q) & (q <= A[:, 1]) & (A[:, 2] <= s) & (s <= A[:, 3]))


miss = [(q, s) for q, s in zip(Q, S) if not (inside(q, s) or inside(1-q, 1-s))]
# 额外：对角线、四条边、四个角这些「盒子接缝」上的点单独测
edge = [(t, t) for t in np.linspace(0, 1, 501)]
edge += [(0.0, t) for t in np.linspace(0, 1, 201)] + [(t, 1.0) for t in np.linspace(0, 1, 201)]
edge += [(1.0, t) for t in np.linspace(0, 1, 201)] + [(t, 0.0) for t in np.linspace(0, 1, 201)]
miss_e = [(q, s) for q, s in edge if not (inside(q, s) or inside(1-q, 1-s))]
print(f"  随机内点 {npt}：未覆盖 {len(miss)}；接缝点 {len(edge)}（对角线+四边）：未覆盖 {len(miss_e)}")
for p in (miss[:3] + miss_e[:3]):
    print(f"    MISS {p}")

# ── 箭头②：逐叶方向 ──────────────────────────────────────────────────────────
print()
print("GATE B4b  逐叶方向：叶内成员的 float64 λ₁（N=500）必须 ≤ 该叶的认证上界")
t0 = time.time(); bad = 0; worst_gap = 1e18; worst_at = None; mx_all = -1e18
for idx, (qlo, qhi, slo, shi, tool, b) in enumerate(Lv):
    pts = [(qlo, slo), (qlo, shi), (qhi, slo), (qhi, shi),
           (0.5*(qlo+qhi), 0.5*(slo+shi))]
    pts += [(rng.uniform(qlo, qhi), rng.uniform(slo, shi)) for _ in range(3)]
    for (q, s) in pts:
        l1 = KT.lam_k(*QS.mem_qs(lam, wE, q, s), mus, Nc, 500)
        mx_all = max(mx_all, l1)
        if l1 > b + 1e-12:
            bad += 1
            print(f"    **FAIL** 叶 {idx} q[{qlo:.4f},{qhi:.4f}] σ[{slo:.4f},{shi:.4f}] "
                  f"工具{tool} 上界 {b:.9f} < λ₁ {l1:.9f}")
        if b - l1 < worst_gap:
            worst_gap = b - l1; worst_at = (idx, qlo, qhi, slo, shi, tool, b, l1)
print(f"  {len(Lv)} 叶 × 8 点 = {8*len(Lv)} 次：方向违反 {bad}；"
      f"全体 float64 max λ₁ = {mx_all:.9f}；{time.time()-t0:.0f}s")
print(f"  最紧的一叶（上界−λ₁ 最小 = {worst_gap:.3e}）：{worst_at}")
