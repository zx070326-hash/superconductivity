# -*- coding: utf-8 -*-
"""
主会话从零：细化 (a) —— 把「μ*>0 Einstein 极值性未证」从**口头前提**变成**可算预算**。

判据表（300 K 需要多大 η/M）建在 Einstein 族上，其方向靠
「Einstein 在固定 (λ,ω̄₂) 下使 Tc 最大」：μ*=0 是 analytic（定理 2），**μ*>0 未证**。
口头挂一句「未证」没有信息量。本脚本给两个数：

━━ (a1) 违反有多大？（在判据表**实际用到**的区域搜反例）━━
  固定 (λ,ω̄₂) 的多原子族全扫。**便宜筛 + 昂贵确认**：
    把温度**钉在** T = Tc_E，只算一次 λ₁(M)。λ₁>1 ⟺ 该谱在 Tc_E 仍超导 ⟹ Tc ≥ Tc_E。
    ⟹ 一次 eigh 而非 44 次二分；`λ₁−1` 的**符号**就是判词，候选再做完整 Tc 确认。
  ⚠ **正对照承重**：memory（`mustar-extremality-weak-coupling-failure-2026-09-17`）记
    「固定 (λ,⟨ω²⟩) 下软化升 Tc ⟺ λ < 2μ*_rms，边界 ≈2.5μ*_rms；精确算子在
     μ=0.30,r=5 处 λ=0.50/0.48 违反、0.52/0.55 不违反」。
    ⟹ 我**必须**在那一带看见违反；看不见 ⟹ 本检验对违反是盲的 ⟹ 判据区的「0 违反」作废。
    正对照按 **μ* 扫**（0.13/0.20/0.30）× λ 跨越 2μ*，不硬套别人的 μ*_rms 约定。
  ⚠ **N_c ≤ N 门**：Coulomb 截断格数 N_c=⌊ω_c/2πT⌋ 超过矩阵维数 N 时库仑项覆盖全矩阵，
    结果不可信（memory 记 μ*=0.13 的失效线要 N_c~2500–5000）。逐行印 N_c/N，超界标 ⚠。

━━ (a2) 万一违反，判据表要动多少？（前提的**预算**）━━
  若某谱对 Einstein 有增益 g>1，它用「Einstein 达 300/g K 所需的 η/M」即可到 300 K。
  ⟹ 把目标温度按 g 退化，读 η/M 的下移量。按 CLAUDE.md「预算优先」：
     若 g=1.25 也不翻转判决 ⟹ 该前提**不承重**，只买余量不买通过，不必去证。

━━ 方向审计 ━━
· (a1) 只能**证伪**极值性，不能证明。判词：反例 / 「N 点未找到」，不写「已证」。
· (a2) 是单调退化：g↑ ⟹ 目标温度↓ ⟹ 需求 η/M ↓ ⟹ 判据**变松**。
  ⟹ (a2) 给的是「前提若失效，判据最多松到哪」，**这一侧不依赖极值性、无条件**。
· 两侧 ω_c 都以同一个 meV 数传入 ⟹ 约定在比值里抵消。

━━ 两矩 2 原子族（analytic）━━
  p := λ/2 = Σa_j/w_j ,  q := λω̄₂²/2 = Σa_j w_j
  ⟹ a1 = w1(p·w2²−q)/(w2²−w1²) ,  a2 = w2(q−p·w1²)/(w2²−w1²)；非负 ⟺ w1 ≤ ω̄₂ ≤ w2。
  固定 w1<ω̄₂ 后 w2 仍自由 ⟹ 族是二维的，全扫。（重建须逐位命中，断言门 1e-10。）

严格性：numerical-evidence。
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

KB = 8.617333262e-2
HB = 64.65415134
rng = np.random.default_rng(20260923)


def top_eig_G(T, G, N, mus, wc):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N)
    L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= G[0] * u * u
    if mus > 0.0:
        Nc = int(np.floor(wc / (2.0 * np.pi * T)))
        if Nc > 0:
            uc = np.where(n < Nc, u, 0.0)
            M -= 2.0 * mus * np.outer(uc, uc)
    return eigh(M, eigvals_only=True, subset_by_index=[N - 1, N - 1])[0]


def G_atoms(kmax, T, aw):
    x2 = (2.0 * np.pi * T * np.arange(kmax, dtype=float)) ** 2
    G = np.zeros(kmax)
    for a, w in aw:
        G += 2.0 * a * w / (w * w + x2)
    return G


def ev_atoms(T, aw, N, mus, wc):
    return top_eig_G(T, G_atoms(2 * N + 2, T, aw), N, mus, wc)


def Tc_atoms(aw, N, mus, wc, lo=0.02, hi=1200.0):
    f = lambda T: ev_atoms(T, aw, N, mus, wc) - 1.0
    lo, hi = lo * KB, hi * KB
    if not (f(lo) > 0 > f(hi)):
        return np.nan
    for _ in range(46):
        mid = np.sqrt(lo * hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi) / KB


def two_atom(lam, wbar, w1, wB):
    p, q = lam / 2.0, lam * wbar ** 2 / 2.0
    d = wB ** 2 - w1 ** 2
    if d <= 1e-12:
        return None
    a1 = w1 * (p * wB ** 2 - q) / d
    a2 = wB * (q - p * w1 ** 2) / d
    if a1 < -1e-14 or a2 < -1e-14:
        return None
    return [(max(a1, 0.0), w1), (max(a2, 0.0), wB)]


def moments_atoms(aw):
    lam = 2.0 * sum(a / w for a, w in aw)
    return lam, np.sqrt(2.0 * sum(a * w for a, w in aw) / lam)


print("=" * 104)
print("两矩 2 原子族自检（analytic 闭式，数值须逐位命中）")
print("=" * 104)
for lam, wb, w1, wB in [(2.0, 150.0, 80.0, 300.0), (3.5, 120.0, 119.0, 121.0), (1.5, 200.0, 5.0, 4000.0)]:
    aw = two_atom(lam, wb, w1, wB)
    lt, vt = moments_atoms(aw)
    e = max(abs(lt / lam - 1), abs(vt / wb - 1))
    print(f"  (λ,ω̄₂)=({lam},{wb})  w1={w1:7.1f} w2={wB:7.1f}  a=({aw[0][0]:.5f},{aw[1][0]:.5f})  误差 {e:.2e}")
    assert e < 1e-10, "两矩闭式不自洽，停"
print("  ⟹ 全部 < 1e-10\n")

WB, N = 150.0, 400          # ω̄₂ 标度无关，固定；N 见 N_c 门


def screen(lam, mus, r, N=N, n3=150):
    """钉温筛：T = Tc_E，逐谱算 λ₁−1。返回 (Tc_E, Nc, best_excess, best_param, nviol, ntot)。"""
    wc = WB * r
    E = [(lam * WB / 2.0, WB)]
    tE = Tc_atoms(E, N, mus, wc)
    if not np.isfinite(tE):
        return np.nan, 0, np.nan, None, 0, 0
    T = tE * KB
    Nc = int(np.floor(wc / (2.0 * np.pi * T)))
    base = ev_atoms(T, E, N, mus, wc)                 # 应 ≈ 1（二分残差）
    best, nv, nt = (-np.inf, None), 0, 0
    w1s = WB * np.concatenate([[1e-3, 1e-2, 0.05, 0.1], np.linspace(0.2, 0.95, 7),
                               1.0 - np.array([3e-2, 1e-2, 3e-3, 1e-3, 1e-4])])
    w2s = WB * np.concatenate([1.0 + np.array([1e-4, 1e-3, 3e-3, 1e-2, 3e-2]),
                               np.linspace(1.1, 3.0, 7), [5.0, 10.0, 30.0, 100.0]])
    for w1 in w1s:
        for wB in w2s:
            aw = two_atom(lam, WB, w1, wB)
            if aw is None:
                continue
            x = ev_atoms(T, aw, N, mus, wc) - base
            nt += 1
            if x > 1e-12:
                nv += 1
            if x > best[0]:
                best = (x, (round(w1 / WB, 5), round(wB / WB, 5)))
    p, q = lam / 2.0, lam * WB ** 2 / 2.0
    for _ in range(n3):
        ws = np.sort(WB * np.exp(rng.uniform(-3.0, 3.0, 3)))
        a3 = rng.uniform(0.0, q / ws[2])
        Mx = np.array([[1.0 / ws[0], 1.0 / ws[1]], [ws[0], ws[1]]])
        if abs(np.linalg.det(Mx)) < 1e-12:
            continue
        a12 = np.linalg.solve(Mx, np.array([p - a3 / ws[2], q - a3 * ws[2]]))
        if np.any(a12 < 0):
            continue
        aw = [(a12[0], ws[0]), (a12[1], ws[1]), (a3, ws[2])]
        lt, vt = moments_atoms(aw)
        if max(abs(lt / lam - 1), abs(vt / WB - 1)) > 1e-9:
            continue
        x = ev_atoms(T, aw, N, mus, wc) - base
        nt += 1
        if x > 1e-12:
            nv += 1
        if x > best[0]:
            best = (x, ("3atom",) + tuple(np.round(ws / WB, 3)))
    return tE, Nc, best[0], best[1], nv, nt


# ── 正对照：按 μ* 扫，λ 跨越 2μ* 边界 ────────────────────────────────
print("=" * 104)
print("(a1-正对照)  memory：固定 (λ,⟨ω²⟩) 下软化升 Tc ⟺ λ ≲ 2–2.5 μ*_rms ⟹ 此带**应当**出现违反")
print("=" * 104)
print(f"  {'μ*':>6}{'r':>4}{'λ':>7}{'Tc_E[K]':>10}{'N_c':>8}{'N_c/N':>8}"
      f"{'max(λ₁−1)':>13}{'违反/总':>10}   参数")
pos_fire, pos_rows = 0, 0
for mus, r in [(0.30, 5.0), (0.20, 5.0), (0.13, 5.0)]:
    for lam in [0.9 * 2 * mus, 1.1 * 2 * mus, 1.4 * 2 * mus, 2.0 * 2 * mus]:
        tE, Nc, x, prm, nv, nt = screen(lam, mus, r)
        if not np.isfinite(tE):
            print(f"  {mus:6.2f}{r:4.0f}{lam:7.3f}   —  (Einstein 本身不超导)"); continue
        flag = " ⚠N_c>N" if Nc > N else ""
        print(f"  {mus:6.2f}{r:4.0f}{lam:7.3f}{tE:10.3f}{Nc:8d}{Nc/N:8.2f}"
              f"{x:13.3e}{f'{nv}/{nt}':>10}   {prm}{'  ← 反例' if nv else ''}{flag}", flush=True)
        if Nc <= N:
            pos_rows += 1
            pos_fire += (nv > 0)

# ── 判据区 ─────────────────────────────────────────────────────────
print()
print("=" * 104)
print("(a1-判据区)  300 K 判据表实际用到的 λ∈[1.5,5]，μ* = 0.13")
print("=" * 104)
print(f"  {'μ*':>6}{'r':>4}{'λ':>7}{'Tc_E[K]':>10}{'N_c':>8}{'N_c/N':>8}"
      f"{'max(λ₁−1)':>13}{'违反/总':>10}   参数")
crit_v, crit_t, crit_x, crit_rows = 0, 0, -np.inf, 0
for r in (6.0, 10.0):
    for lam in [1.5, 2.0, 3.0, 4.0, 5.0]:
        tE, Nc, x, prm, nv, nt = screen(lam, 0.13, r)
        flag = " ⚠N_c>N" if Nc > N else ""
        print(f"  {0.13:6.2f}{r:4.0f}{lam:7.3f}{tE:10.2f}{Nc:8d}{Nc/N:8.2f}"
              f"{x:13.3e}{f'{nv}/{nt}':>10}   {prm}{'  ← 反例' if nv else ''}{flag}", flush=True)
        if Nc <= N:
            crit_rows += 1
            crit_v += nv; crit_t += nt; crit_x = max(crit_x, x)

print()
print("-" * 104)
if pos_fire == 0:
    print(f"⚠⚠ VOID：{pos_rows} 个合法正对照行（N_c≤N）**全部未出现违反**。")
    print("   memory 记该带极值性真失效 ⟹ 本检验对违反是盲的 ⟹ 判据区的「0 违反」**不作数**。")
else:
    print(f"✓ 正对照亮了：{pos_fire}/{pos_rows} 个合法行出现违反 ⟹ 本检验对违反**不盲**。")
    print(f"⟹ 判据区可读：{crit_rows} 行、违反 {crit_v}/{crit_t}、max(λ₁−1) = {crit_x:.3e}")
    print("   判词：**未找到反例**（有限族，不写「已证」）。")
print("-" * 104)

# ── (a2) 前提预算 ──────────────────────────────────────────────────
print()
print("=" * 104)
print("(a2) 前提预算：若极值性以增益 g 失效，300 K 的 η/M 需求下移多少？（g↑⟹判据变松；此侧无条件）")
print("=" * 104)


def wE_for_Tc(T0_K, lam, mus, r, Nb=400):
    T = T0_K * KB
    k = np.arange(2 * Nb + 2, dtype=float)
    f = lambda w: top_eig_G(T, lam * w ** 2 / (w ** 2 + (2.0 * np.pi * T * k) ** 2),
                            Nb, mus, r * w) - 1.0
    lo, hi = 1.0, 30000.0
    if f(hi) < 0:
        return np.nan
    for _ in range(46):
        mid = np.sqrt(lo * hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


GS = [1.00, 1.02, 1.05, 1.10, 1.25, 1.50]
print(f"  {'λ':>6}" + "".join(f"{'g=%.2f' % g:>12}" for g in GS) + f"{'g=1.5 降幅':>14}")
for lam in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0]:
    v = []
    for g in GS:
        w = wE_for_Tc(300.0 / g, lam, 0.13, 10.0)
        v.append(lam * w ** 2 / HB ** 2 if np.isfinite(w) else np.nan)
    print(f"  {lam:6.2f}" + "".join(f"{x:12.3f}" for x in v)
          + f"{(v[-1] / v[0] - 1) * 100:13.1f}%", flush=True)
print("\n  单位 η/M [eV/Å²/amu]；对照 H₃S-200 实测 η/M = 9.172。")
print("=" * 104)
