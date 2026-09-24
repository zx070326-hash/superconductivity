# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ② 成本估算：**自适应细分要多少个盒子？**

GATE M 已判：(x,t) 与 (x,u) 都**不是**能让均匀网格通吃的度规
（λ₁ 沿 x 的快变尺度是 1/ρ²，在 x→0 的 Einstein 退化边附近任选一种固定度规都会有坏盒子）。
⟹ ② 的正确形状是**自适应区间分支定界**（标准做法），不是均匀网格。
本脚本估它的规模：细分到「盒内上界 ≤ 1+ε」为止，数盒子。

⚠ **本脚本给的不是证书，是规模估计。** 判据在**角点**上算，角点不外包盒内部
⟹ 得到的盒子数是真实（区间算术）盒子数的**下估**。真正的 ② 必须把
`G`、`B`、`D`、`μ`、`res²`、`γ` 全部换成区间量，那时每个盒子还会更保守、细分更深。
标这一条是因为 CLAUDE.md：单边界的方向要审——这里方向是「低估工作量」，
所以本表只能用来说「② 是可付的工程量」，**不能**用来说「② 已经做完了」。

两个工具（都 fail-closed，见 `_chk_main_two_regime_coverage_2026-09-23.py`）：
  工具 1：λ₁(M) ≤ λ₁(𝔎_ph) ≤ R₀[w=u]        （无前提，远端用）
  工具 2：Kato–Temple + γ=R[w=u]，需 μ>γ     （近 Einstein 用，δ~1e−9）

━━━━━━━━ 实测结论（2026-09-23 跑完后回填）：**本代理判据是非检验，盒子数 UNCOMPUTED** ━━━━━━━━
基线 (λ=3.2953, μ\*=0.13) 下，**根盒子（深度 0）四角判据就通过**：上界 1.000028420 ≤ 1+ε
⟹ 叶子数 = 1。这不是「② 只要 1 个盒子」，是**这个代理没有内部分辨率**：
它只看四个角，而本问题的 sup 恰好落在边界（x→0 / x→1 / t=1 三条 Einstein 退化边）上，
于是角点判据在基线处**恒真**，被测自由度为 0。
功效对照（投毒：把族的 λ 整体放大）：λ×1.00 过、×1.05 起全部报否（角点上界 1.0676 → 3.0568，
与族内真 max 逐位相符），⟹ 它对**全局 λ 放大**有牙，但那种投毒**同时移动了边界**，
所以这条功效**不构成**「它能看见内部」的证据。
⟹ 判词：**② 的盒子数无法用角点代理估计，报 UNCOMPUTED。**
   要得到真数字，必须先把 `G/B/D/μ/res²/γ` 换成区间量（真正的 ②(c)），本轮未做。
（这是 CLAUDE.md「先数自由度 / 先算功效再跑」的又一起：100% 通过率先查它是不是恒真。）
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


KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")
PL = _load("_chk_main_lam2_peeling_certificate_2026-09-23")
CV = _load("_chk_main_two_regime_coverage_2026-09-23")
TV = _load("_chk_main_trial_vector_reuse_boxes_2026-09-23")


def bound_at(lam, wE, mus, Nc, x, t, vc, gam, N_R, Nv):
    L_, W_ = TV.mem_xt(lam, wE, x, t)
    r0, _ = CV.R0(L_, W_, N_R)
    if r0 <= 1.0 + 1e-3:
        return r0
    b, _ = TV.kt_with_v(L_, W_, mus, Nc, vc, Nv, gam)
    return np.inf if b is None else min(r0, b)


def main():
    EPS, r = 1e-3, 10.0
    N_R, Nv = 600, 120
    MAXDEPTH = 12
    for lam, mus in ((3.2953, 0.13),):
        wE = KT.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
        print(f"λ={lam} μ*={mus} r={r} ⟹ ω̄₂={wE:.9f} N_c={Nc}；坐标 (x, t=1/ρ) ∈ (0,1]×[0,1]")
        print(f"ε={EPS:g}，最大深度={MAXDEPTH}，N_R={N_R}, Nv={Nv}")
        print()
        stack = [(0.5, 0.5, 0.5, 0)]        # (xc, tc, h, depth)
        nbox = nleaf = nfail = 0
        bydepth = {}
        worst = None
        while stack:
            xc, tc, h, dep = stack.pop()
            nbox += 1
            corners = [(np.clip(xc+sx*h, 1e-4, 1-1e-4), np.clip(tc+st*h, 0.0, 1.0))
                       for sx in (-1, 1) for st in (-1, 1)]
            Lc, Wc = TV.mem_xt(lam, wE, xc, tc)
            vc = TV.vstar(Lc, Wc, mus, Nc, Nv)
            gam = max(PL.R_of(*TV.mem_xt(lam, wE, cx, ct), N_R, 0, PL.weights("u", N_R))[0]
                      for cx, ct in corners)
            vals = [bound_at(lam, wE, mus, Nc, cx, ct, vc, gam, N_R, Nv) for cx, ct in corners]
            mx = max(vals)
            if mx <= 1.0 + EPS:
                nleaf += 1
                bydepth[dep] = bydepth.get(dep, 0) + 1
                if worst is None or mx > worst[0]:
                    worst = (mx, xc, tc, h, dep)
                continue
            if dep >= MAXDEPTH:
                nfail += 1
                bydepth[dep] = bydepth.get(dep, 0) + 1
                continue
            for sx in (-1, 1):
                for st in (-1, 1):
                    stack.append((xc + sx*h/2, tc + st*h/2, h/2, dep+1))
        print(f"  访问盒子 {nbox}，**叶子（通过）{nleaf}**，达最大深度仍不过 {nfail}")
        print(f"  叶子按深度分布：{dict(sorted(bydepth.items()))}")
        if worst:
            print(f"  最劣叶子：上界={worst[0]:.9f} 于 (x={worst[1]:.5f}, t={worst[2]:.5f}, h={worst[3]:.2e}, 深度={worst[4]})")
        print()
        print("  读法：叶子数 = ② 需要跑的区间 KT/Schur 次数的**数量级估计**（真实值更大，见文件头）。")
        print("        nfail>0 ⟹ 该区域用角点判据都不收敛，需要更强的工具或更深的细分，必须逐个列出。")


if __name__ == "__main__":
    main()
