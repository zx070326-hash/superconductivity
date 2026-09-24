# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ②(b) 第三半：**两区制覆盖**——KT 在远端成员上失效，谁来接？

K1 暴露的事实（不是 bug，是工具边界）：`(λ=2, x=.1, ρ=20)` 处 `μ=0.0637 ≤ γ_cert=0.0656`
⟹ Kato–Temple 的 `μ−γ>0` 前提不成立 ⟹ KT **不适用**。
机制清楚：远离 Einstein 的成员 `λ₁(M)` 本身就掉到 ~0.06，而剥指标界 γ_cert≈0.07 是**粗界**，
两者撞在一起。KT 是「把 λ₁ 夹到 1e−7」的精密工具，用在 λ₁≈0.06 的地方是杀鸡用牛刀且够不着。

⟹ 设计上应当是**两区制**，各用各的工具，各自 fail-closed：
  · 工具 1（粗、无前提）：`λ₁(M) ≤ λ₁(𝔎_ph) ≤ R₀ := sup_i [(1/u_i)(B u)_i − D_ii]`
    —— **不剥指标**的加权 Schur 检验（B 逐元非负 ⟹ 合法）。第一步用了
    `M = 𝔎_ph − 2μ*u_cu_cᵀ ⪯ 𝔎_ph`（库仑项 PSD），也是 fail-closed。
  · 工具 2（精密、有前提）：KT + γ_cert，要求 `μ > γ_cert`。
本脚本判的是：**两个工具的适用域并起来盖不盖得住整族？**
判据：对每个成员，`min(可用工具给出的上界) ≤ 1+ε`（ε=1e−3）是否成立。

⚠ 覆盖判决对**每个成员**给的是 PASS/未覆盖；「未覆盖」不等于该成员违反天花板，
只说明这两个工具都不够——那是工具问题，不是结论问题（单边界方向规则）。

严格性：两条不等式链均为 `analytic`；本脚本一切数值为 `numerical-evidence`（float64）。
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


def R0(lams, wbars, N):
    """不剥指标的加权 Schur 界：λ₁(𝔎_ph) ≤ sup_i [(1/u_i)(Bu)_i − D_ii]，w=u。"""
    B, D = PL.BD(lams, wbars, N)
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    row = (B @ u)/u - D
    return float(row.max()), int(np.argmax(row))


def main():
    r, EPS = 10.0, 1e-3
    N_R, N_ref = 1200, 1000
    print("=" * 116)
    print("GATE C1  不剥指标的粗界 R₀ 本身有多大？（Einstein 处必然 >1，因 λ₁(𝔎_ph)>λ₁(M)=1）")
    print(f"  {'case':<22}{'λ₁(M)':>12}{'λ₁(𝔎_ph)':>13}{'R₀':>12}{'R₀/λ₁(𝔎_ph)':>15}{'R₀≤1+ε?':>11}")
    bases = {}
    for lam, mus in ((2.0, 0.13), (3.2953, 0.13), (6.0, 0.20)):
        wE = KT.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
        bases[(lam, mus)] = (wE, Nc)
        l1M = KT.lam_k([lam], [wE], mus, Nc, N_ref)
        l1P = KT.lam_k([lam], [wE], 0.0, 0, N_ref)
        v, i0 = R0([lam], [wE], N_R)
        print(f"  λ={lam:<6.4g} μ*={mus:<5}{l1M:>12.6f}{l1P:>13.6f}{v:>12.4f}"
              f"{v/l1P:>15.4f}{'YES' if v <= 1+EPS else 'no':>11}   argmax i={i0}")
    print("  ⟹ R₀ 在 Einstein 附近（λ₁≈1）**够不着**，必须由 KT 接手；问题只剩远端。")
    print()

    print("=" * 116)
    print("GATE C2  两区制覆盖：对每个成员，(工具1: R₀) 与 (工具2: KT, 需 μ>γ_cert) 谁能给出 ≤1+ε 的上界？")
    print(f"  ε={EPS:g}。每格打印：族内被覆盖的成员数 / 总数，以及两工具各自的覆盖数、重叠数、缺口")
    xs = np.linspace(0.02, 0.98, 13)
    rhos = np.concatenate([[1.0], np.geomspace(1.02, 40.0, 12)])
    for lam, mus in ((2.0, 0.13), (3.2953, 0.13), (6.0, 0.20)):
        wE, Nc = bases[(lam, mus)]
        g_ct, _ = PL.R_of([lam], [wE], N_R, 0, PL.weights("u", N_R))   # Einstein 处的 γ（仅供参考）
        n_t1 = n_t2 = n_cov = n_gap = 0
        gaps = []
        maxlam1 = -1e9
        for x in xs:
            for rho in rhos:
                L_, W_ = KT.members(lam, wE, x, rho)
                ref = KT.lam_k(L_, W_, mus, Nc, N_ref)
                maxlam1 = max(maxlam1, ref)
                r0, _ = R0(L_, W_, N_R)
                gam, _ = PL.R_of(L_, W_, N_R, 0, PL.weights("u", N_R))
                t1 = (r0 <= 1.0 + EPS)
                kt, mu, _, _ = KT.kato_temple(L_, W_, mus, Nc, 120, gam)
                t2 = (kt is not None) and (kt <= 1.0 + EPS)
                n_t1 += t1; n_t2 += t2
                if t1 or t2:
                    n_cov += 1
                else:
                    n_gap += 1
                    gaps.append((x, rho, ref, r0, mu, gam, kt))
        tot = len(xs)*len(rhos)
        print(f"  λ={lam:<6.4g} μ*={mus:<5}  覆盖 {n_cov}/{tot}"
              f"   [工具1(R₀≤1+ε): {n_t1}   工具2(KT≤1+ε): {n_t2}   缺口: {n_gap}]"
              f"   族内 max λ₁(M)={maxlam1:.9f}")
        for g in gaps[:6]:
            x, rho, ref, r0, mu, gam, kt = g
            ks = "N/A(μ≤γ)" if kt is None else f"{kt:.6f}"
            print(f"       缺口 x={x:.3f} ρ={rho:.3f}:  λ₁(M)={ref:.6f}  R₀={r0:.4f}"
                  f"  μ={mu:.6f} γ={gam:.6f} KT={ks}")
        if len(gaps) > 6:
            print(f"       … 另有 {len(gaps)-6} 个缺口")
    print()
    print("  读法：缺口 = 这两个工具都不够，**不是**该成员违反天花板（单边界方向规则）。")
    print("        缺口若集中在 ρ≈1 附近，说明需要的是 KT 而 γ_cert 太粗；")
    print("        缺口若集中在远端，说明需要的是更好的粗界（如剥指标版 R 直接界 λ₁ 的变体）。")


if __name__ == "__main__":
    main()
