# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）任务 ②(c) 前置：**换坐标 (q,σ) + 亏损恒等式**，并顺手测一个
「如果成立就不用 B&B 了」的强命题。

━━━━━━━━ 动机 ━━━━━━━━
§10.6 的判词是「没有能让均匀网格通吃的固定度规」，那是在 (x,t=1/ρ) 坐标里得到的。
但两原子族有一组**天然坐标**：按两个约束各自的**份额**来参数化——
    q := λ₂/λ                （第二个原子占的 λ 份额，= ∫ω⁻¹dν 的份额）
    σ := λ₂ω₂²/M2            （第二个原子占的 M2 份额，= ∫ωdν 的份额）
于是 ω₂ = ω̄₂√(σ/q)、ω₁ = ω̄₂√((1−σ)/(1−q))，(q,σ) ∈ [0,1]²，**整块正方形**。
这组坐标的好处（本脚本逐条查）：
  · 紧化后的**全部**极点（两原子族 + 单原子曲线 (A)/(B) + ν=0）都落在这个正方形里：
        σ=q  ⟹ ω₁=ω₂=ω̄₂        （Einstein，**整条对角线**）
        q=0  ⟹ 单原子 (λ, ω̄₂√(1−σ))，只饱和第一条约束  = 曲线 (A)
        σ=1  ⟹ 单原子 (qλ, ω̄₂/√q)，只饱和第二条约束  = 曲线 (B)
        (q,σ)=(0,1) ⟹ G≡0 = ν=0
    ⟹ §10.1 里「2 原子族 ∪ 两条曲线 ∪ {0}」合并成**一个闭正方形**，不再有「角落里的曲线」。
  · G(k) 对 (q,σ) 的依赖是**两个逐变量单调的 Möbius 项之和** ⟹ 盒上区间包围可取到（几乎）精确。
  · **亏损恒等式**（本脚本 GATE P2）：
        G_E(k) − G(k;q,σ) = λ ω̄₂⁴ k² (σ−q)² / [ (ω̄₂²+k²) · u · (W−u) ],
        u := ω̄₂²σ + k²q,  W := ω̄₂²+k²,  W−u = ω̄₂²(1−σ) + k²(1−q)
    ⟹ 亏损对「离对角线的距离」是**恰好二次**的 ⟹ 跨对角线的盒子里区间宽度是二阶量，
       这正是 §10.6 那个「一阶包围 vs 二阶真变化」的病根所在，换坐标后消失。

━━━━━━━━ 附带的解析引理（GATE P3） ━━━━━━━━
上面的恒等式只对两原子族成立，但同一结论对**整个 𝒞̄** 有一行解析证明：
  y ↦ 2y/(y+k²) 在 y=ω² 上凹 ⟹ 取 y=ω̄₂² 处的切线 ⟹  g_k(ω)=2ω/(ω²+k²) ≤ α/ω + βω，
  α = 2ω̄₂⁴/(ω̄₂²+k²)² ≥ 0,  β = 2k²/(ω̄₂²+k²)² ≥ 0，两系数非负 ⟹ 可用 𝒞̄ 的**两条 ≤ 约束**：
      G(k;ν) = ∫g_k dν ≤ α·(λ/2) + β·(M2/2) = λω̄₂²/(ω̄₂²+k²) = G_E(k)   （k ≥ 1）
  ⟹ **Einstein 的 G 逐 k 支配全族**（`analytic`）。

⚠ **但这不蕴含 λ₁ ≤ λ₁_E**：算子 M 对 G **不单调**（对角 −S_i 与非对角同向增大）。
GATE P4 就是去测那个「如果成立就不用 B&B」的强命题：
      M[ν] = M[E] − 𝒦₀[Δ],  Δ := G_E − G ≥ 0     （库仑项与 ν 无关，在差里恰好抵消）
      若 𝒦₀[Δ] ⪰ 0 ⟹ M[ν] ⪯ M[E] ⟹ λ₁(ν) ≤ λ₁(E)，整个定理解析化。
先验上**预期为假**（对正测度 ν′，𝒦₀[G′]=𝔎_ph 本身就是不定的：λ₁>0 而对角很负），
但 Δ 不是任何正测度的 G，所以值得实测一次。⟹ 本门是**朝「省掉 B&B」方向的一次尝试**，
失败就照常走 ②(c)；**不允许**用它的失败去动摇任何已有结论。

严格性：GATE P3 的不等式为 `analytic`（上面的切线论证）；本脚本一切数值为
`numerical-evidence`（float64）。GATE P2 的恒等式是代数恒等式（手推），数值只排查笔误。
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


# ───────────────────────── 坐标 ─────────────────────────

def mem_qs(lam, wE, q, s):
    """(q,σ) ⟹ (lams, wbars)。q=λ₂/λ、σ=λ₂ω₂²/M2。退化端用连续延拓（权重为 0 的原子丢掉）。"""
    lams, ws = [], []
    if 1.0 - q > 0 and 1.0 - s > 0:
        lams.append((1.0 - q)*lam); ws.append(wE*np.sqrt((1.0 - s)/(1.0 - q)))
    if q > 0 and s > 0:
        lams.append(q*lam); ws.append(wE*np.sqrt(s/q))
    if not lams:
        return np.array([0.0]), np.array([wE])          # ν = 0
    return np.array(lams), np.array(ws)


def qs_from_xrho(x, rho):
    """旧坐标 (x,ρ)（`KT.members`）⟹ 新坐标 (q,σ)。"""
    return x, x*rho*rho/((1.0 - x) + x*rho*rho)


def G_of(lams, wbars, k):
    return KT.Gk(lams, wbars, k)


def G_E(lam, wE, k):
    k = np.asarray(k, float)
    return lam*wE*wE/(wE*wE + k*k)


def deficit_closed(lam, wE, q, s, k):
    """亏损的闭式：G_E(k) − G(k;q,σ)。"""
    k = np.asarray(k, float)
    W = wE*wE + k*k
    u = wE*wE*s + k*k*q
    return lam*wE**4*k*k*(s - q)**2/(W*u*(W - u))


# ───────────────────────── 𝒦₀：只由 {G(k)}_{k≥1} 决定的算子 ─────────────────────────

def K0(gfun, N):
    """𝒦₀[g]：对角用**抵消后**的写法 u_i²(g(2i+1) − 2Σ_{j=1..i}g(j))，**不含 g(0)**。
    （G(0) 在对角上恰好抵消——见 [[G0-does-not-enter-the-operator]]。）"""
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n + 1.0)
    g = np.concatenate([[np.nan], gfun(np.arange(1, 2*N + 2))])   # g[0] 置 nan：确保没人用到它
    i, j = np.meshgrid(n, n, indexing="ij")
    off = g[np.abs(i - j)] + g[i + j + 1]
    np.fill_diagonal(off, 0.0)
    cum = np.concatenate([[0.0], np.cumsum(g[1:N + 1])])          # Σ_{j=1..i} g(j)
    diag = g[2*n + 1] - 2.0*cum[:N]
    return np.outer(u, u)*(off + np.diag(diag))


def M_op(lam, wE, q, s, mus, Nc, N):
    """M[ν] = 𝒦₀[G] − 2μ* u_c u_cᵀ，用 (q,σ) 成员。"""
    lams, ws = mem_qs(lam, wE, q, s)
    K = K0(lambda kk: G_of(lams, ws, kk), N)
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n + 1.0)
    uc = np.where(n < Nc, u, 0.0)
    return K - 2.0*mus*np.outer(uc, uc)


def main():
    r, mus, lam = 10.0, 0.13, 3.2953
    wE = KT.einstein_root(lam, mus, r, 600); Nc = int(np.floor(r*wE))
    print(f"基准：λ={lam} μ*={mus} r={r} ⟹ ω̄₂={wE:.9f}  N_c={Nc}  M2=λω̄₂²={lam*wE*wE:.9f}")
    print()

    # ───────────────────────── GATE P0 ─────────────────────────
    print("=" * 112)
    print("GATE P0  新坐标 (q,σ) 与旧坐标 (x,ρ) 给出**同一个成员**（λ₁ 逐位相同）")
    print(f"  {'(x,ρ)':<16}{'(q,σ)':<24}{'λ₁ 旧坐标':>16}{'λ₁ 新坐标':>16}{'|差|':>12}")
    worst = 0.0
    for x in (0.1, 0.4, 0.7, 0.9):
        for rho in (1.0, 1.3, 3.0, 20.0):
            q, s = qs_from_xrho(x, rho)
            L1, W1 = KT.members(lam, wE, x, rho)
            L2, W2 = mem_qs(lam, wE, q, s)
            a = KT.lam_k(L1, W1, mus, Nc, 400)
            b = KT.lam_k(L2, W2, mus, Nc, 400)
            worst = max(worst, abs(a - b))
            if rho in (1.0, 3.0):
                print(f"  x={x:<4} ρ={rho:<7.3g}({q:.6f}, {s:.6f})   {a:>15.12f}{b:>16.12f}{abs(a-b):>12.1e}")
    print(f"  ⟹ 16 组最劣 |差| = {worst:.2e}")
    print()

    # ───────────────────────── GATE P0b：边界 = 两条曲线 + ν=0 ─────────────────────────
    print("=" * 112)
    print("GATE P0b  正方形的四条边就是 §10.1 的紧化极点集（逐条核对频率与权重）")
    for (q, s, tag) in ((0.37, 0.37, "对角 σ=q  ⟹ Einstein"),
                        (0.0, 0.55, "左边 q=0   ⟹ 曲线(A) 单原子 (λ, ω̄₂√(1−σ))"),
                        (0.42, 1.0, "上边 σ=1   ⟹ 曲线(B) 单原子 (qλ, ω̄₂/√q)"),
                        (0.0, 1.0, "角 (0,1)   ⟹ ν=0")):
        L_, W_ = mem_qs(lam, wE, q, s)
        mom1 = float(np.sum(L_))                      # 2∫ω⁻¹dν = Σλ_j
        mom2 = float(L_ @ (W_**2))                    # 2∫ωdν   = Σλ_jω_j²
        print(f"  (q,σ)=({q:.2f},{s:.2f})  λ_j={np.array2string(L_, precision=5)}  "
              f"ω_j={np.array2string(W_, precision=5)}")
        print(f"      Σλ_j={mom1:.6f} (≤λ={lam})   Σλ_jω_j²={mom2:.6f} (≤M2={lam*wE*wE:.6f})   {tag}")
    print()

    # ───────────────────────── GATE P1：单调结构（盒上区间包围可取精确） ─────────────────────────
    print("=" * 112)
    print("GATE P1  G(k) = λω̄₂²[F(1−σ,1−q) + F(σ,q)]，F(a,b)=ab/(ω̄₂²a+k²b) 对 a、b 各自**不减**")
    print("  ⟹ 每项在盒上的 sup/inf 都在角点取到（逐变量只出现一次 ⟹ 该项的区间包围是精确的）")
    rng = np.random.default_rng(20260923)
    bad = 0
    for _ in range(4000):
        a, b, k = rng.uniform(1e-3, 1), rng.uniform(1e-3, 1), rng.integers(1, 60)
        F = lambda A, B: A*B/(wE*wE*A + k*k*B)
        if F(a + 1e-6, b) < F(a, b) - 1e-15 or F(a, b + 1e-6) < F(a, b) - 1e-15:
            bad += 1
    print(f"  4000 组随机 (a,b,k) 的单调性违反：{bad}")
    # 与直接构造的成员对照
    err = 0.0
    for q in (0.05, 0.3, 0.8):
        for s in (0.1, 0.5, 0.95):
            L_, W_ = mem_qs(lam, wE, q, s)
            ks = np.arange(1, 200)
            direct = G_of(L_, W_, ks)
            via = lam*wE*wE*((1-s)*(1-q)/(wE*wE*(1-s) + ks*ks*(1-q)) + s*q/(wE*wE*s + ks*ks*q))
            err = max(err, float(np.max(np.abs(direct - via))))
    print(f"  与直接求和的最大偏差：{err:.2e}")
    print()

    # ───────────────────────── GATE P2：亏损恒等式 ─────────────────────────
    print("=" * 112)
    print("GATE P2  亏损恒等式 G_E−G = λω̄₂⁴k²(σ−q)²/[(ω̄₂²+k²)·u·(W−u)]（代数恒等式，数值只查笔误）")
    print(f"  {'(q,σ)':<20}{'max_k |恒等式−直接|':>24}{'max_k 相对':>16}{'G≤G_E?':>10}")
    worst_rel = 0.0
    for (q, s) in ((0.2, 0.2), (0.2, 0.25), (0.05, 0.6), (0.5, 0.9), (0.8, 0.99), (0.01, 0.02)):
        L_, W_ = mem_qs(lam, wE, q, s)
        ks = np.arange(1, 400)
        d_direct = G_E(lam, wE, ks) - G_of(L_, W_, ks)
        d_closed = deficit_closed(lam, wE, q, s, ks)
        ae = float(np.max(np.abs(d_direct - d_closed)))
        re = float(np.max(np.abs(d_direct - d_closed)/np.maximum(np.abs(d_direct), 1e-300)))
        worst_rel = max(worst_rel, re)
        print(f"  ({q:.2f},{s:.2f}){ae:>26.2e}{re:>16.2e}{'YES' if d_direct.min() >= -1e-15 else '**NO**':>10}")
    print(f"  ⟹ 最劣相对偏差 {worst_rel:.2e}")
    print()

    # ───────────────────────── GATE P3：切线引理（对整个 𝒞̄） ─────────────────────────
    print("=" * 112)
    print("GATE P3  切线引理 g_k(ω) ≤ α/ω + βω（α,β ≥ 0，ω=ω̄₂ 处取等）—— 这条对**整个 𝒞̄** 给 G ≤ G_E")
    nbad = 0; worst_gap = 1e18; tight = 0
    for k in range(1, 200):
        al = 2*wE**4/(wE*wE + k*k)**2
        be = 2*k*k/(wE*wE + k*k)**2
        if al < 0 or be < 0:
            nbad += 1
        ws = np.concatenate([np.geomspace(1e-4, 1e4, 4000), [wE]])
        gap = al/ws + be*ws - 2*ws/(ws*ws + k*k)
        nbad += int(np.sum(gap < -1e-14*np.maximum(1.0, al/ws + be*ws)))
        worst_gap = min(worst_gap, float(gap.min()))
        tight += int(abs(gap[-1]) < 1e-14)
    print(f"  k=1..199 × 4001 个 ω：α或β 为负 / 不等式违反 共 {nbad} 处；最劣间隙 {worst_gap:.2e}")
    print(f"  ω=ω̄₂ 处取等的 k 数：{tight}/199（切点必须取等，否则界不紧、G_E 不是该界的值）")
    print("  ⟹ (analytic) 对任何 ν∈𝒞̄ 与 k≥1：G(k;ν) ≤ α·(λ/2)+β·(M2/2) = G_E(k)。")
    print("     注意：**这不蕴含 λ₁ ≤ λ₁_E**（M 对 G 不单调）——GATE P4 去测那个更强的命题。")
    print()

    # ───────────────────────── GATE P4：M[ν] ⪯ M[E] ? ─────────────────────────
    print("=" * 112)
    print("GATE P4  强命题：𝒦₀[Δ] ⪰ 0 ?（若真 ⟹ M[ν] ⪯ M[E] ⟹ λ₁ ≤ λ₁_E，整定理解析化）")
    N = 300
    print(f"  N={N}。印 λ_min(𝒦₀[Δ])、λ_max(𝒦₀[Δ])，以及 λ₁(M[ν]) 与 λ₁(M[E]) 的实际差")
    print(f"  {'(q,σ)':<16}{'λ_min(𝒦₀[Δ])':>16}{'λ_max(𝒦₀[Δ])':>16}{'λ₁(M[E])−λ₁(M[ν])':>22}{'⪰0?':>8}")
    MEs = M_op(lam, wE, 0.3, 0.3, mus, Nc, N)          # 对角线上任一点都是 Einstein
    l1E = float(eigh(MEs, subset_by_index=[N-1, N-1], eigvals_only=True)[0])
    npsd = 0; ntot = 0
    for (q, s) in ((0.3, 0.32), (0.3, 0.5), (0.05, 0.6), (0.5, 0.9), (0.8, 0.99),
                   (0.01, 0.02), (0.1, 0.11), (0.45, 0.55), (0.02, 0.9)):
        dfun = lambda kk: deficit_closed(lam, wE, q, s, kk)
        KD = K0(dfun, N)
        ev = eigh(KD, eigvals_only=True)
        Mv = M_op(lam, wE, q, s, mus, Nc, N)
        l1 = float(eigh(Mv, subset_by_index=[N-1, N-1], eigvals_only=True)[0])
        ntot += 1; npsd += int(ev.min() >= -1e-12)
        print(f"  ({q:.2f},{s:.2f}){ev.min():>16.3e}{ev.max():>16.3e}{l1E - l1:>22.9f}"
              f"{'YES' if ev.min() >= -1e-12 else 'no':>8}")
    print(f"  ⟹ {npsd}/{ntot} 个成员满足 𝒦₀[Δ] ⪰ 0")
    print("  读法：只要有**一个** no，强命题即被证伪 ⟹ 回到 ②(c) 的分支定界（不影响任何已有结论）；")
    print("        若全 YES，则须**另行**证明它（抽样不构成证明），且要先找它可能失败的方向。")


if __name__ == "__main__":
    main()
