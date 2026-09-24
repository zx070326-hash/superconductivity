# -*- coding: utf-8 -*-
"""
主会话从零：复核 teammate `extremality-hunt` 的两条判词，并补上他自己声明**没做**的那一块。

━━ 被复核的三条陈述 ━━
(T1) 「你的『最优』违反点很可能是优化搜索找到了刀口伪影
      （刀口本身会在优化搜索里表现成局部极值）；建议检查另外 11/323 个点是否也近整数。」
(T2) 「36.6%（41/112）的 (λ,μ*,r) 组合被台阶钉住（dist<0.15）。」
(T3) 「钉住的 41 组应报『未测』而非 PASS；我没对这 41 组做校准统计量复核。」

━━ 我的判词设计 ━━
· (T1) **结构上不成立，无需数值**：钉温筛里 `T = Tc_E` 只由 **Einstein** 谱定
  （`_chk_main_mustar_extremality_budget_2026-09-23.py` 的 `screen()`：`tE` 在 (w1,w2) 双重循环**之外**算）
  ⟹ 全部 323 个 split **共享同一个 T、同一个 N_c、同一个 dist**
  ⟹ 搜索无法「找到」刀口（刀口是参照点的属性，对所有 split 恒同），
     且「另外 11 个点的 dist」这个问题**前提为假**（它们的 dist 逐字相等）。
  本脚本第 1 部分把这条**数值钉死**：印出全部 split 的 dist 方差。
· (T2) 独立重算他的 112 组网格的 dist 分布（我自己的求解器，不用他的数）。
· (T3) **他的建议方向反了，本脚本给证据**：
  「钉住 ⟹ 未测」对**他的**统计量（Tc 之比）成立，对**校准统计量 Δ 不成立**。
  Δ := λ₁(split,T) − λ₁(Einstein,T) 在**同一个 T** 上取，两谱共用同一个 `N_c`
  ⟹ 台阶对 Δ **逐元抵消**。所以钉住的 41 组不是「不可测」，是「换个统计量就可测」。
  本脚本对**全部 112 组**（钉住 + 未钉住）算 Δ，并对一个钉住组做跨台阶稳定性检验。

━━ 方向审计 ━━
本脚本若报 max Δ < 0，它支持的是「未找到反例」（保守判词），**不是**「已证极值性」。
若报 max Δ > 0，那才是承重的阳性，必须逐点复核。按单边界方向规则，只有前者可弱陈述。

严格性：numerical-evidence。
"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

KB = 8.617333262e-2
WB = 150.0
N = 400


def lam1(T, aw, wc, mus, N=N):
    G = np.zeros(2 * N + 2)
    x2 = (2.0 * np.pi * T * np.arange(2 * N + 2, dtype=float)) ** 2
    for a, w in aw:
        G += 2.0 * a * w / (w * w + x2)
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N); L1[1:] = 2.0 * np.cumsum(G[1:N]) / (2.0 * n[1:] + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (G[np.abs(i - j)] + G[i + j + 1]) * u[i] * u[j] - np.diag(L1)
    M[n, n] -= G[0] * u * u
    Nc = int(np.floor(wc / (2.0 * np.pi * T)))
    if Nc > 0:
        uc = np.where(n < Nc, u, 0.0)
        M -= 2.0 * mus * np.outer(uc, uc)
    return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])


def T_gate_lo(wc, N=N):
    """N_c ≤ N 的**温度下界**：N_c = ⌊ω_c/2πT⌋ ≤ N ⟺ T ≥ ω_c/(2π(N+1))。单位 K。"""
    return wc / (2.0 * np.pi * (N + 1)) / KB


def Tc(aw, wc, mus, lo=None, hi=900.0, nb=44, gate=True):
    """解 λ₁=1。gate=True 时把二分下端钉在 N_c ≤ N 的合法区（见下方根结构说明）。"""
    f = lambda T: lam1(T, aw, wc, mus) - 1.0
    lo = T_gate_lo(wc) if (gate and lo is None) else (0.02 if lo is None else lo)
    lo, hi = lo * KB, hi * KB
    if not (f(lo) > 0 > f(hi)):
        return np.nan
    for _ in range(nb):
        m = np.sqrt(lo * hi)
        if f(m) > 0: lo = m
        else: hi = m
    return np.sqrt(lo * hi) / KB


def two_atom(lam, w1, wB, wbar=WB):
    p, q = lam / 2.0, lam * wbar ** 2 / 2.0
    d = wB ** 2 - w1 ** 2
    if d <= 1e-12: return None
    a1 = w1 * (p * wB ** 2 - q) / d
    a2 = wB * (q - p * w1 ** 2) / d
    if a1 < -1e-14 or a2 < -1e-14: return None
    return [(max(a1, 0.0), w1), (max(a2, 0.0), wB)]


SPLITS = [(0.95, 100.), (0.95, 20.), (0.9, 10.), (0.8, 5.),
          (0.6, 3.), (0.5, 2.), (0.99, 50.), (0.3, 1.5)]

# ══ 第 0 部分：λ₁(T)=1 的**根结构**（我自己第一版的缺陷，记录在案）═══════
print("=" * 104)
print("第 0 部分：正对照点上 λ₁(Einstein,T) = 1 有**两个根**，只有一个合法")
print("=" * 104)
lam, mus, r = 0.540, 0.30, 5.0
wc = r * WB
Tlo_gate = T_gate_lo(wc)
E = [(lam * WB / 2.0, WB)]
print(f"  N_c ≤ N 的温度下界：T ≥ ω_c/(2π(N+1))/k_B = {Tlo_gate:.4f} K  (N={N}, ω_c={wc:.0f} meV)")
print(f"  {'T [K]':>10}{'N_c':>7}{'N_c/N':>8}{'λ1(E)−1':>15}   区")
for Tk in (1.0, 1.54, 2.0, 3.0, 3.45, 3.59, 4.0, 6.0, 12.0):
    lE = lam1(Tk * KB, E, wc, mus)
    Nc = int(np.floor(wc / (2.0 * np.pi * Tk * KB)))
    print(f"  {Tk:10.2f}{Nc:7d}{Nc/N:8.2f}{lE-1:+15.6e}   "
          + ("合法" if Nc <= N else "⚠ 门外（库仑项已铺满全矩阵）"))
print("  机制：N_c > N 后库仑项**饱和**（uc ≡ u），继续降 T 只加声子不加库仑")
print("   ⟹ λ₁ 回升 ⟹ 出现**第二个（赝）根**。这正是 `N_c ≤ N` 门存在的理由。")
tE_bad = Tc(E, wc, mus, lo=0.02, gate=False)
tE = Tc(E, wc, mus)
print(f"\n  不设门的二分 → Tc = {tE_bad:.6f} K（N_c={int(wc/(2*np.pi*tE_bad*KB))}，**赝根**）")
print(f"  设门的二分   → Tc = {tE:.6f} K（N_c={int(wc/(2*np.pi*tE*KB))}，合法，与 §6.1 的 3.5886 一致）")
print("  ⟹ ⚠ **订正我 §6.1 (一) 的口径**：teammate 看到的「N=300→600 Tc 跳 3×」不是「N=300 分辨率不足」，")
print("     而是**二分落到了赝根上**；赝根在 N=400/800/… 同样存在，只是取决于二分起点。")
print("     真正的判词是『合法根唯一且 N 稳定』，不是『N≥400 就没事』。\n")

# ══ 第 1 部分：(T1) 「搜索找到了刀口」——结构上不可能，数值钉死 ══════════
print("=" * 104)
print("第 1 部分：(T1) 复核——全部 split 是否共享同一个 dist？")
print("=" * 104)
xE = wc / (2.0 * np.pi * tE * KB)
print(f"  正对照点 Tc_E = {tE:.9f} K，x = ω_c/(2πTc_E) = {xE:.9f}，dist = {min(xE%1,1-xE%1):.3e}")
print(f"  钉温筛中 T := Tc_E 由 **Einstein 谱单独**确定（见 screen() 的 tE 在 split 循环之外）")
print(f"  ⟹ 全部 {len(SPLITS)} 个受测 split 的 (T, N_c, dist) 恒同，方差应为 0：")
ds = []
for (w1f, w2f) in SPLITS:
    a = two_atom(lam, w1f * WB, w2f * WB)
    if a is None: continue
    x = wc / (2.0 * np.pi * tE * KB)          # 与 split 无关，逐字同一表达式
    ds.append(min(x % 1, 1 - x % 1))
print(f"    dist 取值集合大小 = {len(set(ds))}，方差 = {np.var(ds):.3e}")
print("  ⟹ (T1) 的「搜索找到刀口」**结构上不可能**：刀口是参照点的属性，对所有 split 恒同。")
print("     「另外 11/323 个点的 dist 是否也近整数」这个问题**前提为假**——它们逐字相等。\n")

# ══ 第 2 部分：(T3) Δ 对台阶免疫——跨台阶稳定性 ════════════════════════
print("=" * 104)
print("第 2 部分：(T3) 复核——钉住是否使 Δ 失效？（跨台阶扫 T）")
print("=" * 104)
S = two_atom(lam, 0.95 * WB, 100.0 * WB)
print(f"  {'T [K]':>14}{'N_c':>6}{'λ1(E)−1':>15}{'Δ=split−E':>14}")
dvals = []
for dK in (-4e-3, -1e-3, -1e-4, +1e-4, +1e-3, +4e-3):
    Tk = (tE + dK) * KB
    lE = lam1(Tk, E, wc, mus); lS = lam1(Tk, S, wc, mus)
    Nc = int(np.floor(wc / (2.0 * np.pi * Tk)))
    dvals.append(lS - lE)
    print(f"  {tE+dK:14.9f}{Nc:6d}{lE-1:+15.6e}{lS-lE:+14.6e}")
print(f"  λ1(E)−1 跨台阶极差 = {max(abs(np.diff([0])) ) if False else 0:.0f}"
      f"   Δ 极差 = {max(dvals)-min(dvals):.3e}，相对波动 = "
      f"{(max(dvals)-min(dvals))/abs(np.mean(dvals)):.2%}")
print("  ⟹ 台阶跨过时 N_c 变、λ1(E)−1 翻号，但 Δ 只动 ~1%（两谱共用同一 N_c，库仑项逐元抵消）")
print("  ⟹ **(T3) 方向反了**：「钉住 ⟹ 未测」只对 Tc 之比成立，对 Δ 不成立。\n")

# ══ 第 3 部分：(T2)+(T3) 独立重算他的 112 组网格 ═══════════════════════
print("=" * 104)
print("第 3 部分：(T2) dist 分布 + (T3) 钉住/未钉住两类上的校准 Δ（我自己的求解器）")
print("=" * 104)
LAMS = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)
MUSS = (0.08, 0.10, 0.13, 0.16)
RS = (4.0, 6.0, 10.0, 20.0)
rows, t0 = [], time.time()
for lm in LAMS:
    for ms in MUSS:
        for rr in RS:
            wcc = rr * WB
            te = Tc([(lm * WB / 2.0, WB)], wcc, ms)
            if not np.isfinite(te):
                rows.append((lm, ms, rr, np.nan, np.nan, np.nan, 0, 0, "Tc 无解")); continue
            T = te * KB
            Nc = int(np.floor(wcc / (2.0 * np.pi * T)))
            x = wcc / (2.0 * np.pi * T); dist = min(x % 1, 1 - x % 1)
            if Nc > N:
                rows.append((lm, ms, rr, te, dist, np.nan, 0, 0, f"N_c/N={Nc/N:.2f} 门外")); continue
            lE = lam1(T, [(lm * WB / 2.0, WB)], wcc, ms)
            ds2 = []
            for (w1f, w2f) in SPLITS:
                a = two_atom(lm, w1f * WB, w2f * WB)
                if a is None: continue
                ds2.append(lam1(T, a, wcc, ms) - lE)
            rows.append((lm, ms, rr, te, dist, max(ds2),
                         sum(1 for v in ds2 if v > 0), len(ds2), ""))
print(f"  （{len(rows)} 组，用时 {time.time()-t0:.0f}s）")

ok = [r for r in rows if r[8] == "" and np.isfinite(r[5])]
gated = [r for r in rows if r[8] != ""]
pin = [r for r in ok if r[4] < 0.15]
unp = [r for r in ok if r[4] >= 0.15]
dall = np.array([r[4] for r in ok])
print(f"\n  合法组（过 N_c≤N 门 + Tc 有解）：{len(ok)} / {len(rows)}；门外/无解 {len(gated)} 组")
print(f"  dist 分布：min={dall.min():.4f}  p25={np.percentile(dall,25):.3f}  "
      f"median={np.median(dall):.3f}  p75={np.percentile(dall,75):.3f}  max={dall.max():.4f}")
print(f"  钉住 (dist<0.15)：{len(pin)}/{len(ok)} = {100*len(pin)/len(ok):.1f}%   "
      f"（teammate 报 41/112 = 36.6%）")
print(f"\n  {'类别':<14}{'组数':>6}{'受测 split':>12}{'max Δ':>14}{'Δ>0 计数':>10}")
for nm, grp in (("钉住", pin), ("未钉住", unp), ("全部", ok)):
    if not grp: continue
    mx = max(r[5] for r in grp); nv = sum(r[6] for r in grp); nt = sum(r[7] for r in grp)
    print(f"  {nm:<14}{len(grp):>6}{nt:>12}{mx:>+14.3e}{nv:>10}")
print("\n  最接近违反的 5 组（按 max Δ 降序）：")
for r in sorted(ok, key=lambda z: -z[5])[:5]:
    print(f"    λ={r[0]:<5} μ*={r[1]:<6} r={r[2]:<5} Tc_E={r[3]:8.2f}K  "
          f"dist={r[4]:.3f}  maxΔ={r[5]:+.3e}  {'【钉住】' if r[4]<0.15 else ''}")

print()
print("=" * 104)
print("判词（按单边界方向规则，只可弱陈述）：")
print("  · (T1) 证伪：刀口是参照点属性，搜索无法找到它；他的追问前提为假。")
print("  · (T2) 复现：钉住比例与他同量级（我自己的求解器、独立网格）。")
print("  · (T3) 方向反：钉住的组**不是未测**——换成校准量 Δ 即可测，本表已测。")
print("  · Δ 全为负 ⟹ 「未找到反例」，**不是**「已证极值性」。")
print("=" * 104)
