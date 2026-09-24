# -*- coding: utf-8 -*-
"""
从零复核 einstein-margin 的「径向探测：deficit/ε² 收敛到常数」，并把它升级成一个
**有自由度**的检验（2026-09-24，主会话）。

━━━ 为什么原报告的判词自由度很低 ━━━
`_chk_main_qsigma_coords_and_deficit_identity_2026-09-23.py` 的 GATE P2 已经写下：
    Δ(k) := G_E(k) − G(k;q,σ) = λ ω̄₂⁴ k² (σ−q)² / [ (ω̄₂²+k²) · u · (W−u) ]
    u := ω̄₂²σ + k²q,  W := ω̄₂²+k²
—— 亏损在 **G 层面**对 (σ−q) 是**恰好二次**的（代数恒等式）。于是
「λ₁ 的亏损 / ε² 在 ε→0 收敛到有限正数」几乎是这条恒等式 + λ₁ 的可微性的直接推论，
测出来「收敛」并不构成对任何物理陈述的检验。（原报告自己观察到的 1.0555→1.0283 漂移
恰恰说明 λ₁ 层不是恰好二次——那是 O(ε³) 项，符合预期。）

━━━ 升级成有自由度的检验：闭式系数 + q 不变量 ━━━
沿 σ=q+ε 令 ε→0：u → qW、W−u → (1−q)W，代入得**领头项闭式**
    Δ(k)/ε²  →  Δ̂(k) := λ ω̄₂⁴ k² / [ q(1−q) (ω̄₂²+k²)³ ]         （k≥0，Δ̂(0)=0）
又 M[ν] = M[E] − 𝒦₀[Δ]（库仑项与 ν 无关，在差里恰好抵消），一阶微扰给
    λ₁(E) − λ₁(q,q+ε)  =  ⟨v_E, 𝒦₀[Δ̂] v_E⟩ · ε²  +  O(ε³)
其中 v_E 是 Einstein 算子的顶本征向量——**v_E 只依赖 (λ, ω̄₂, μ*, N_c)，与 q 无关**。
⟹ 全部 q 依赖都在显式因子 1/(q(1−q)) 里 ⟹ **预言**：

    P(q) := lim_{ε→0} deficit(q,ε)/ε²  ×  q(1−q)     对 q **恒定**，
    且该常数 = q(1−q)·⟨v_E,𝒦₀[Δ̂]v_E⟩ 有闭式（与 q 无关的表达式）。

这条有真自由度：q 在 [0.2,0.85] 上变 4 倍多，1/(q(1−q)) 变 ~2.7 倍，若机制判错就对不上。
它同时解释了原报告表 A 与表 B 的「不一致」：q≈0.83 处 1.8995 vs q=0.5 处 1.028，
比值应为 (0.5·0.5)/(0.83·0.17) = 1.772 —— 不是「系数依赖切入方向」这种含糊说法。

━━━ 方向审计 ━━━
本脚本不产出任何上界/下界判决，只做「机制识别」：全部结论 numerical-evidence（float64），
闭式领头项本身是 analytic（上面一行代数极限）。**不得**用本脚本的任何数字去
声称「Einstein 领先 X」——那正是被撤下的那类陈述。

严格性：闭式 Δ̂ = analytic；一切数值 = numerical-evidence。
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
QS = _load("_chk_main_qsigma_coords_and_deficit_identity_2026-09-23")

M2_BASE = 1.453151466
MUS, R = 0.13, 10.0


def lam1_of(lams, wbars, Nc, N):
    M, _, _ = KT.build(lams, wbars, MUS, Nc, N)
    return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])


def einstein_pair(lam, N):
    """固定 M2 = M2_BASE ⟹ ω̄₂ = √(M2/λ)；返回 (ω̄₂, Nc, λ₁_E, v_E)。"""
    w = np.sqrt(M2_BASE / lam)
    Nc = int(np.floor(R * w))
    M, _, _ = KT.build([lam], [w], MUS, Nc, N)
    e, V = eigh(M, subset_by_index=[N - 1, N - 1])
    return w, Nc, float(e[0]), V[:, 0]


def K0_of_kernel(gvals, N):
    """把一个核 k ↦ g(k)（长度 2N+2 的取值表）送进 𝒦₀：与 KT.build 同构、无库仑项。"""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    S = gvals[0] + 2.0 * np.concatenate([[0.0], np.cumsum(gvals[1:N])])
    return (gvals[np.abs(i - j)] + gvals[i + j + 1]) * np.outer(u, u) - np.diag(u * u * S)


def Delta_hat(lam, w, q, N):
    """闭式领头项 Δ̂(k) = λ ω̄₂⁴ k² / [ q(1−q) (ω̄₂²+k²)³ ]。"""
    k = np.arange(2 * N + 2, dtype=float)
    return lam * w**4 * k * k / (q * (1.0 - q) * (w * w + k * k) ** 3)


def deficit(lam, w, Nc, lam1E, q, eps, N):
    lams, ws = QS.mem_qs(lam, w, q, q + eps)
    return lam1E - lam1_of(lams, ws, Nc, N)


def main():
    N = 800
    print("=" * 94)
    print(f"(I) 从零重跑 einstein-margin 的径向探测（q=0.5, σ=0.5+ε）  N={N}, μ*={MUS}, r={R}")
    print("=" * 94)
    print("      ε        λ=2.00 deficit/ε²     λ=2.50 deficit/ε²")
    rows = {}
    for lam in (2.00, 2.50):
        w, Nc, l1E, _ = einstein_pair(lam, N)
        rows[lam] = (w, Nc, l1E)
        print(f"   [λ={lam:.2f}]  ω̄₂={w:.9f}  Nc={Nc}  λ₁(Einstein)={l1E:.12f}")
    for eps in (0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001):
        cells = []
        for lam in (2.00, 2.50):
            w, Nc, l1E = rows[lam]
            cells.append(deficit(lam, w, Nc, l1E, 0.5, eps, N) / eps**2)
        print(f"   {eps:7.4f}      {cells[0]:14.6f}      {cells[1]:14.6f}")

    print()
    print("=" * 94)
    print("(II) GATE D1  闭式一阶微扰系数 vs ε→0 外推（q=0.5）")
    print("=" * 94)
    for lam in (2.00, 2.50):
        w, Nc, l1E, vE = einstein_pair(lam, N)
        cf = float(vE @ (K0_of_kernel(Delta_hat(lam, w, 0.5, N), N) @ vE))
        num = deficit(lam, w, Nc, l1E, 0.5, 1e-4, N) / 1e-8
        print(f"   λ={lam:.2f}:  闭式 ⟨v_E,𝒦₀[Δ̂]v_E⟩ = {cf:.9f}   "
              f"ε=1e-4 实测 = {num:.9f}   相对差 = {abs(cf-num)/cf:.2e}")

    print()
    print("=" * 94)
    print("(III) GATE D2  q 不变量：P(q) := (deficit/ε²)·q(1−q) 应与 q 无关（ε=1e-4）")
    print("=" * 94)
    print("      q      1/(q(1−q))   deficit/ε²(λ=2.00)   P(q)         闭式 P        相对差")
    for lam in (2.00, 2.50):
        w, Nc, l1E, vE = einstein_pair(lam, N)
        print(f"   ── λ={lam:.2f} ──")
        for q in (0.20, 0.30, 0.50, 0.70, 0.83):
            d = deficit(lam, w, Nc, l1E, q, 1e-4, N) / 1e-8
            P = d * q * (1.0 - q)
            cf = float(vE @ (K0_of_kernel(Delta_hat(lam, w, q, N), N) @ vE)) * q * (1.0 - q)
            print(f"   {q:5.2f}    {1/(q*(1-q)):8.4f}    {d:16.6f}   {P:10.6f}   "
                  f"{cf:10.6f}   {abs(P-cf)/cf:.2e}")

    print()
    print("=" * 94)
    print("(IV) 解释原报告表 A/表 B 的「不一致」：q≈0.83 的 1.8995 vs q=0.5 的 1.028")
    print("=" * 94)
    qA = 0.831503418
    pred = (0.5 * 0.5) / (qA * (1.0 - qA))
    print(f"   1/(q(1−q)) 之比 = (0.5·0.5)/({qA:.4f}·{1-qA:.4f}) = {pred:.4f}")
    w, Nc, l1E, _ = einstein_pair(2.00, N)
    d05 = deficit(2.00, w, Nc, l1E, 0.5, 1e-4, N) / 1e-8
    print(f"   ⟹ 预言 q={qA:.4f} 处的领头系数 = {d05:.4f} × {pred:.4f} = {d05*pred:.4f}")
    epsA = 0.849240402 - 0.831503418
    dA = deficit(2.00, w, Nc, l1E, qA, epsA, N) / epsA**2
    print(f"   实测（用原报告那一点的 ε={epsA:.6e}）= {dA:.4f}"
          f"   （原报告报 1.8995；差异是 O(ε) 高阶项，ε 不小）")
    dA0 = deficit(2.00, w, Nc, l1E, qA, 1e-4, N) / 1e-8
    print(f"   同一 q 上取 ε=1e-4 = {dA0:.4f}  ⟹ 与预言比 = {dA0/(d05*pred):.6f}")

    print()
    print("=" * 94)
    print("(V) N 收敛（GATE D1 的闭式 vs 实测，q=0.5）")
    print("=" * 94)
    for NN in (400, 600, 800):
        line = []
        for lam in (2.00, 2.50):
            w, Nc, l1E, vE = einstein_pair(lam, NN)
            cf = float(vE @ (K0_of_kernel(Delta_hat(lam, w, 0.5, NN), NN) @ vE))
            line.append(f"λ={lam:.2f}: 闭式={cf:.8f}")
        print(f"   N={NN:4d}   " + "   ".join(line))

    print()
    print("=" * 94)
    print("(VI) 把**我自己**记的那个 1.79（λ=3.29 的 N_c 扫描）也换成不变量")
    print("=" * 94)
    print("   P(λ) := (deficit/ε²)·q(1−q) —— 与 q 无关的那个量（这才是可记的）：")
    for lam in (2.00, 2.50, 3.29):
        w, Nc, l1E, vE = einstein_pair(lam, N)
        P = float(vE @ (K0_of_kernel(Delta_hat(lam, w, 0.5, N), N) @ vE)) * 0.25
        print(f"     λ={lam:.2f}  ω̄₂={w:.9f}  Nc={Nc}   P = {P:.9f}")
    w, Nc, l1E, vE = einstein_pair(3.29, N)
    P329 = float(vE @ (K0_of_kernel(Delta_hat(3.29, w, 0.5, N), N) @ vE)) * 0.25
    qstar = 0.5 * (1.0 - np.sqrt(max(0.0, 1.0 - 4.0 * P329 / 1.79)))
    print(f"   ⟹ 我记的 1.79 对应 q(1−q) = P/1.79 = {P329/1.79:.6f} ⟹ q ≈ {qstar:.4f} 或 {1-qstar:.4f}")
    print("     （即那个数只是在说「抽样点离对角线多近、落在 q 的哪一侧」，"
          "不含任何关于 Einstein 领先幅度的信息）")


if __name__ == "__main__":
    main()
