# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）：**§9 紧化补丁不完整——硬端 (ω→∞) 同样逃逸**。订正 + 新版 (A3)。

━━━━━━━━ 第一段：旧陈述错在哪 ━━━━━━━━
`convex_extreme_point_reduction_2026-09-23.md` §9 写的修法是
    𝒞̄ = { ν≥0 : ∫ω⁻¹dν **≤** λ/2,  ∫ωdν **=** M2/2 }        （M2 := λ⟨ω²⟩）
它只堵住了**软端**逃逸（质量跑向 ω→0）。**硬端同样漏**：
取两原子 `λ₂_at → 0, ω₂ → ∞` 且 `λ₂_at·ω₂² = const`，则
  · `∫ωdν` 的份额 `λ₂_at ω₂²/2` **一直保着**（等式约束沿列满足）；
  · 但 ν 的质量 `λ₂_at ω₂/2 = const/(2ω₂) → 0` ⟹ 弱\* 极限里这一团**消失**；
  · 极限的 `∫ωdν < M2/2` ⟹ **逃出 §9 的 𝒞̄**。
⟹ §9 的 𝒞̄ **仍不闭**。这是我自己一小时前写下的补丁里的真缺口，不是笔误。

━━━━━━━━ 第二段：新陈述凭什么对 ━━━━━━━━
统一的机制（analytic，一行）：`M` 只经 `{G(k)}_{k≥1}` 依赖 ν（G(0) 不进算子，见
`_chk_main_soft_mode_limit_operator_2026-09-23.py`），而
      g_k(ω) := 2ω/(ω² + (2πTk)²),   k ≥ 1
满足 `g_k(0)=0` 且 `g_k(ω)→0 (ω→∞)` ⟹ **g_k ∈ C₀((0,∞))**
⟹ `ν ↦ (G(k))_{k≥1}` 弱\*连续 ⟹ `ν ↦ M[ν]` 弱\*连续 ⟹ `λ₁` 弱\*连续。
**两端的逃逸都恰好贡献 0**——软端与硬端是同一句话的两半。
正确紧化（两条都松成 ≤）：
      𝒞̄ = { ν≥0 : ∫ω⁻¹dν ≤ λ/2,  ∫ωdν ≤ M2/2 }
它弱\*紧：Cauchy–Schwarz 给 `∫dν ≤ √(∫ω⁻¹dν·∫ωdν) ≤ √(λM2)/2` ⟹ 落在球内，且由 ≤ 闭。

**极点集随之变化（必须重算，不能沿用旧结论）**：Richter–Rogosinski 的原子数上界是
**活跃**约束数 ⟹ 𝒞̄ 的极点只有三类：
  (0) `ν=0`；
  (A) **单原子、只激活 λ 约束**：`λ_at=λ`，`ω ≤ ω̄₂`（因 λω² ≤ M2）；
  (B) **单原子、只激活 M2 约束**：`λ_at=M2/ω² ≤ λ`，`ω ≥ ω̄₂`；
  (C) 两原子、两约束都等号 ⟹ 就是原来的两原子族。
(A)(B) 是**两条显式的一参数 Einstein 曲线**，`ω=ω̄₂` 处交于 Einstein 点。
⟹ 本脚本必须查：**(A) 与 (B) 上 λ₁ 是否都被 Einstein 点支配**。若不是，天花板结论要改。

严格性：上面全部不等式/拓扑论证为 `analytic`；本脚本一切数值为 `numerical-evidence`（float64）。

━━━━━━━━ 附带：② 的盒子度规（V2 反常的解释） ━━━━━━━━
`_chk_main_trial_vector_reuse_boxes_2026-09-23.py` GATE V2 出现「h 变小反而更差」。
argmax 定位到唯一的坏盒子：中心 (x,t)=(0.0625,0.0625)、角点 (x=0.001, ρ=8)。
机制：`ω₁ = ω̄₂/√((1−x)+xρ²)` ⟹ **λ₁ 沿 x 的变化尺度是 1/ρ²，不是 1**
⟹ 均匀 (x,t) 网格在 `x→0` 的 Einstein 退化边上盒子太大。这不是 bug、也不是判决翻转，
是**度规选错**。GATE M 给出不变坐标 `u := xρ²/((1−x)+xρ²)`（= 原子 2 占的 **M2 份额**），
并实测在 (x,u) 度规下同一区域的 δ。
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
TV = _load("_chk_main_trial_vector_reuse_boxes_2026-09-23")


def main():
    r, N_ref = 10.0, 1200
    print("=" * 108)
    print("GATE X0  弱*连续性的机制核对：g_k(ω)=2ω/(ω²+k²) 在两端都趋 0（k≥1）")
    for k in (1, 2, 5):
        for w in (1e-8, 1e-4, 1e-2, 1.0, 1e2, 1e4, 1e8):
            pass
    for k in (1, 3):
        vals = [2*w/(w*w + k*k) for w in (1e-8, 1e-4, 1.0, 1e4, 1e8)]
        print(f"   k={k}:  g_k(1e−8)={vals[0]:.3e}  g(1e−4)={vals[1]:.3e}  g(1)={vals[2]:.4f}"
              f"  g(1e4)={vals[3]:.3e}  g(1e8)={vals[4]:.3e}   （k=0 时 g_0=2/ω 在 0 处发散 ⟹ 软端逃逸）")
    print()

    print("=" * 108)
    print("GATE X1  硬端逃逸实测：λ₂_at→0, ω₂→∞ 且 λ₂_at·ω₂²=const —— 两个矩与 λ₁ 各自怎么动")
    lam, mus = 3.2953, 0.13
    wE = KT.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    M2 = lam*wE*wE
    share = 0.30                      # 原子 2 占的 M2 份额
    print(f"   λ={lam} ω̄₂={wE:.6f} M2={M2:.6f}；原子 2 恒占 M2 的 {share:.0%}")
    print(f"   {'ω₂':>12}{'λ₂_at':>12}{'∫ω⁻¹dν(=Σλ/2)':>17}{'∫ωdν(=M2/2)':>15}"
          f"{'ν 总质量':>12}{'λ₁(M)':>12}")
    for w2 in (2.0, 10.0, 1e2, 1e3, 1e4, 1e6):
        l2 = share*M2/(w2*w2)
        l1 = lam - l2 if lam > l2 else 0.0
        w1 = np.sqrt((1.0-share)*M2/l1)
        L_, W_ = np.array([l1, l2]), np.array([w1, w2])
        mom_inv = float(np.sum(L_))/2.0
        mom_1 = float(L_ @ (W_**2))/2.0
        mass = float(L_ @ W_)/2.0
        v = KT.lam_k(L_, W_, mus, Nc, N_ref)
        print(f"   {w2:>12.3g}{l2:>12.3e}{mom_inv:>17.9f}{mom_1:>15.9f}{mass:>12.6f}{v:>12.8f}")
    print(f"   （参照：λ/2={lam/2:.9f}，M2/2={M2/2:.9f}）")
    print("   ⟹ 等式约束沿列**始终满足**，但 ν 总质量单调掉 ⟹ 弱* 极限丢掉 ∫ωdν 的那一份")
    print("      ⟹ §9 的 𝒞̄（∫ωdν 写成等号）**不闭**，必须松成 ≤。")
    print()

    print("=" * 108)
    print("GATE X2  新极点类 (A)(B)：两条显式一参数 Einstein 曲线上，λ₁ 是否被 Einstein 点支配？")
    print("   (A) λ_at=λ, ω ≤ ω̄₂      (B) λ_at=M2/ω² ≤ λ, ω ≥ ω̄₂     交点 ω=ω̄₂ = Einstein 点")
    for lam_, mus_ in ((2.0, 0.13), (3.2953, 0.13), (6.0, 0.20)):
        wE_ = KT.einstein_root(lam_, mus_, r, 600); Nc_ = int(np.floor(r*wE_))
        M2_ = lam_*wE_*wE_
        l1E = KT.lam_k([lam_], [wE_], mus_, Nc_, N_ref)
        bestA = bestB = -1e9; atA = atB = None
        for f in np.linspace(0.02, 1.0, 50):          # (A): ω = f*ω̄₂
            v = KT.lam_k([lam_], [f*wE_], mus_, Nc_, N_ref)
            if v > bestA:
                bestA, atA = v, f
        for f in np.concatenate([np.linspace(1.0, 6.0, 40), np.geomspace(6.0, 300.0, 20)]):
            v = KT.lam_k([M2_/(f*wE_)**2], [f*wE_], mus_, Nc_, N_ref)
            if v > bestB:
                bestB, atB = v, f
        okA = "OK" if bestA <= l1E + 1e-12 else "**超过**"
        okB = "OK" if bestB <= l1E + 1e-12 else "**超过**"
        print(f"   λ={lam_:<7.4g} μ*={mus_:<5} λ₁_E={l1E:.12f}"
              f" | (A) max={bestA:.12f} at ω/ω̄₂={atA:.3f} 超出={bestA-l1E:+.2e} {okA}"
              f" | (B) max={bestB:.12f} at ω/ω̄₂={atB:.3f} 超出={bestB-l1E:+.2e} {okB}")
    print("   ⟹ (A) 沿 ω 单增、(B) 沿 ω 单减 ⟹ 两条新极点曲线的 sup 都恰在 ω=ω̄₂（Einstein 点）。")
    print("     ν=0 那个极点给 λ₁(M)= −2μ*‖u_c‖²+… < 0，显然不竞争（下面一行实测）。")
    v0 = KT.lam_k([1e-300], [1.0], 0.13, 6, 400)
    print(f"     ν≈0 处 λ₁(M) = {v0:.6e}")
    print()

    print("=" * 108)
    print("GATE M  ② 的盒子度规：坏盒子在 (x,t) 下 vs 在 (x,u) 下（u = 原子 2 的 M2 份额）")
    lam, mus = 3.2953, 0.13
    wE = KT.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))

    def mem_xu(x, u):
        """由 (x,u) 反解成员：ω₁ = ω̄₂√((1−u)/(1−x))，ρ² = u(1−x)/(x(1−u))。"""
        w1 = wE*np.sqrt((1.0-u)/(1.0-x))
        rho2 = u*(1.0-x)/(x*(1.0-u))
        return np.array([(1.0-x)*lam, x*lam]), np.array([w1, np.sqrt(rho2)*w1])

    print("   先看 λ₁ 沿两种坐标的变化幅度（同一物理区域）：")
    print(f"   {'x':>9}{'ρ':>9}{'u':>9}{'λ₁(M)':>13}")
    for x, rho in ((0.0010, 8.0), (0.0625, 8.0), (0.0010, 16.0), (0.0625, 2.0), (0.30, 3.33)):
        L_, W_ = KT.members(lam, wE, x, rho)
        u = x*rho*rho/((1.0-x) + x*rho*rho)
        print(f"   {x:>9.4f}{rho:>9.3f}{u:>9.4f}{KT.lam_k(L_, W_, mus, Nc, N_ref):>13.8f}")
    print()
    print(f"   {'box (x,u) center':<22}{'h':>7}{'max δ(中心 v, 4 角)':>22}{'方向':>8}")
    for (xc, uc) in ((0.0625, 0.3000), (0.0625, 0.6000), (0.3000, 0.5000), (0.0200, 0.8000)):
        for h in (0.10, 0.05):
            corners = [(np.clip(xc+sx*h, 1e-3, 0.999), np.clip(uc+su*h, 1e-3, 0.999))
                       for sx in (-1, 1) for su in (-1, 1)]
            corners = [(cx, cu) for cx, cu in corners if cu > cx]    # ρ≥1 ⟺ u ≥ x
            if not corners:
                continue
            Lc, Wc = mem_xu(xc, uc); vc = TV.vstar(Lc, Wc, mus, Nc, 120)
            gam = max(PL.R_of(*mem_xu(cx, cu), 1200, 0, PL.weights("u", 1200))[0]
                      for cx, cu in corners)
            dmax = -1.0; ok = "OK"
            for cx, cu in corners:
                L_, W_ = mem_xu(cx, cu)
                ref = KT.lam_k(L_, W_, mus, Nc, N_ref)
                b, mu = TV.kt_with_v(L_, W_, mus, Nc, vc, 120, gam)
                if b is None:
                    ok = "μ≤γ(交工具1)"; continue
                dmax = max(dmax, b/ref - 1.0)
                if b < ref - 1e-12:
                    ok = "**FAIL**"
            print(f"   ({xc:.3f},{uc:.3f}){'':<10}{h:>7.2f}{dmax:>22.2e}{ok:>8}")
    print()
    print("   读法：(x,u) 是**矩份额**坐标，两端 Einstein 退化边分别是 u=x（ρ=1）与 x=0；")
    print("         在该度规下盒子不再横跨 1/ρ² 的快变方向。这只解释 V2 的反常，")
    print("         **不改变任何判决**——V2 里那格本来也只是工具不够，不是天花板被违反。")


if __name__ == "__main__":
    main()
