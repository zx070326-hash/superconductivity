# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ②(b) 第四半：**投毒门 + ρ→∞ 端**。

前一脚本给出 169/169 零缺口。**100% 通过率必须先证它不是恒等式自证**
（CLAUDE.md：「把一组数字当作 N 个独立检验之前，先数它的自由度」；
 以及「判据是我的抽样族对退化集给了多少质量」）。

【P1 投毒门】把成员人为改造成**应当被判失败**的样子，查两个工具是否真的报不出 ≤1+ε：
    (a) λ 放大 ×(1+s)：真 λ₁(M) > 1 ⟹ 任何**上界**都必须 >1 ⟹ 必须出缺口。
    (b) ω̄ 缩小（等价于 T 降低）：同样使 λ₁>1。
    若投毒后仍报「覆盖」，则 C2 的 169/169 是假 PASS。
    ※ 注意方向：工具给的是**上界**，上界 >1+ε **不等于** λ₁>1+ε；
      所以投毒门查的是「工具会不会盲目说 YES」，这是覆盖判据的**必要**条件。

【P2 ρ→∞ 端】族是 ρ∈[1,∞)，C2 只扫到 ρ=40。②-a 已证极限成员**恰是单原子
    `(xλ, ω̄_E/√x)`**（保 M2、丢 λ）。这里直接查 ρ∈{40,100,1000,10⁴} 与**极限成员本身**
    是否仍被两工具覆盖 —— 紧化后的闭区间端点必须也在覆盖域里，否则 ② 的盒子盖不满。

严格性：一切数值为 `numerical-evidence`（float64）。
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


def cover(L_, W_, mus, Nc, N_R, N_ref, eps):
    """返回 (是否覆盖, λ₁(M)真值, R₀, KT, μ, γ)。"""
    ref = KT.lam_k(L_, W_, mus, Nc, N_ref)
    r0, _ = CV.R0(L_, W_, N_R)
    gam, _ = PL.R_of(L_, W_, N_R, 0, PL.weights("u", N_R))
    kt, mu, _, _ = KT.kato_temple(L_, W_, mus, Nc, 120, gam)
    t1 = r0 <= 1.0 + eps
    t2 = (kt is not None) and (kt <= 1.0 + eps)
    return (t1 or t2), ref, r0, kt, mu, gam


def main():
    r, EPS = 10.0, 1e-3
    N_R, N_ref = 1200, 1000
    lam, mus = 3.2953, 0.13
    wE = KT.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    print(f"基准：λ={lam} μ*={mus} r={r}  ⟹ ω̄_E*={wE:.9f}  N_c={Nc}   ε={EPS:g}")
    print()

    # ───────── P1 投毒门 ─────────
    print("=" * 112)
    print("GATE P1  投毒门：人为造出 λ₁(M)>1 的成员，两工具**必须**报不出 ≤1+ε（否则 C2 是假 PASS）")
    print(f"  {'投毒':<30}{'λ₁(M)真值':>13}{'R₀':>11}{'KT':>14}{'覆盖?':>8}{'期望':>10}{'判':>10}")
    poisons = []
    for s in (0.00, 0.02, 0.05, 0.20, 0.50):
        poisons.append((f"λ ×(1+{s:.2f})  [Einstein]", np.array([lam*(1+s)]), np.array([wE])))
    for s in (0.02, 0.05, 0.20):
        poisons.append((f"ω̄ ×(1−{s:.2f})  [Einstein]", np.array([lam]), np.array([wE*(1-s)])))
    # 双原子投毒：先取一个原本被覆盖的成员，再放大 λ
    L0, W0 = KT.members(lam, wE, 0.4, 3.0)
    for s in (0.00, 0.20, 0.60):
        poisons.append((f"2atom x=.4 ρ=3, λ×(1+{s:.2f})", L0*(1+s), W0))
    nbad = 0
    for tag, L_, W_ in poisons:
        cov, ref, r0, kt, mu, gam = cover(L_, W_, mus, Nc, N_R, N_ref, EPS)
        expect = "应出缺口" if ref > 1.0 + EPS else "可覆盖"
        bad = (ref > 1.0 + EPS) and cov          # 真值超标却报覆盖 = 假 PASS
        nbad += bad
        ks = "N/A(μ≤γ)" if kt is None else f"{kt:.6f}"
        print(f"  {tag:<30}{ref:>13.6f}{r0:>11.4f}{ks:>14}"
              f"{('YES' if cov else 'GAP'):>8}{expect:>10}{('**假PASS**' if bad else 'ok'):>10}")
    print(f"  ⟹ 假 PASS 计数 = {nbad}   （必须为 0；>0 则 C2 的 169/169 不可信）")
    print()

    # ───────── P2 ρ→∞ 端 ─────────
    print("=" * 112)
    print("GATE P2  ρ→∞ 端：紧化后的边界成员（②-a 已证 = 单原子 (xλ, ω̄_E/√x)）是否也被覆盖？")
    print(f"  {'x':>6}{'ρ':>10}{'λ₁(M)':>12}{'R₀':>10}{'KT':>13}{'μ':>10}{'γ':>10}{'覆盖?':>8}")
    ngap = 0
    for x in (0.05, 0.2, 0.5, 0.8, 0.95):
        for rho in (40.0, 100.0, 1000.0, 1e4, np.inf):
            if np.isinf(rho):
                L_, W_ = np.array([x*lam]), np.array([wE/np.sqrt(x)])   # ②-a 的极限成员
            else:
                L_, W_ = KT.members(lam, wE, x, rho)
            cov, ref, r0, kt, mu, gam = cover(L_, W_, mus, Nc, N_R, N_ref, EPS)
            ngap += (not cov)
            ks = "N/A(μ≤γ)" if kt is None else f"{kt:.6f}"
            rs = "∞(极限)" if np.isinf(rho) else f"{rho:.0f}"
            print(f"  {x:>6.2f}{rs:>10}{ref:>12.6f}{r0:>10.4f}{ks:>13}{mu:>10.6f}{gam:>10.6f}"
                  f"{('YES' if cov else '**GAP**'):>8}")
    print(f"  ⟹ 远端缺口数 = {ngap}")
    print()
    print("  读法：P1 若出现「假 PASS」，C2 作废；P2 的缺口若非 0，②的盒子覆盖不到闭包端点。")
    print("  两门都过**不**等于 ② 已完成：C2/P2 都只是**网格点**，")
    print("  ② 真正要的是 (x,1/ρ) 的**区间盒子细分 + 区间算术**，本轮全部为 float64。")


if __name__ == "__main__":
    main()
