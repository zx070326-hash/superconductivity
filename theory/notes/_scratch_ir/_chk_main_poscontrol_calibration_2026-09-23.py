# -*- coding: utf-8 -*-
"""
主会话从零：§6 钉温筛的**校准缺陷**（自查，起因于 teammate `extremality-hunt` 的收敛质疑）。

━━ 旧陈述（被订正的对象）━━
笔记 `ceiling/mustar_extremality_two_routes_closed_2026-09-23.md` §6：
  「正对照必须亮：μ*=0.30、r=5、λ=0.540 → **12/323 违反**，max(λ₁−1)=+2.385e−04」

━━ 缺陷（本脚本证明）━━
钉温筛把「违反」判为 `λ₁(split) − 1 > 0`。这默认了 `λ₁(Einstein, Tc_E) = 1`。
**在库仑台阶上这个默认为假**：`N_c = ⌊ω_c/2πT⌋` 是 T 的阶梯函数 ⟹ `λ₁(T)` 有跳跃间断；
弱耦合端 Tc 被**钉**在台阶边上（实测 `ω_c/(2πTc) = 386.000000`，到整数距离 5.7e−9），
此时二分返回的 Tc_E 处 `λ₁(Einstein) − 1 = +1.14e−4 ≠ 0`（**残差**）。
⟹ 旧统计量 = 真统计量 + 残差，**偏置为正**。

正确统计量（本脚本）：`Δ := λ₁(split, T) − λ₁(Einstein, T)`，同 T、同 N_c、同库仑项。

━━ 方向审计（承重，决定哪些结论受损）━━
偏置为正 ⟹ 旧筛**高报**违反 ⟹
  · 判据区「0/3172 未找到反例」：筛子本来偏向找到违反却仍为 0 ⟹ **保守，结论不受损**（反而更强）。
  · 正对照「12/323 亮」：它恰恰是靠这个偏置亮的 ⟹ **受损，必须重算**。
⟹ 这是单边界方向规则的一次直接应用：同一个偏置对两个判决的后果相反。

━━ 前提核查（W1 纪律：订正前先确认符号作用在哪个变量上）━━
`N_c` 只依赖 (T, ω_c)，**不依赖谱形** ⟹ 同 T 下两谱的库仑项逐元相同 ⟹ Δ 的比较是干净的。
（若 N_c 依赖谱，Δ 本身也会被污染，那就得换工具。已核：不依赖。）

严格性：numerical-evidence。
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

KB = 8.617333262e-2
WB = 150.0
MUS = 0.30
N = 600


def lam1(T, aw, wc, mus=MUS, N=N):
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


def Tc(aw, wc, lo=0.05, hi=400.0, nb=60):
    f = lambda T: lam1(T, aw, wc) - 1.0
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
    a1 = w1 * (p * wB ** 2 - q) / d
    a2 = wB * (q - p * w1 ** 2) / d
    if a1 < -1e-14 or a2 < -1e-14:
        return None
    return [(max(a1, 0.0), w1), (max(a2, 0.0), wB)]


# ── 第 1 部分：台阶钉住的直接证据 ──────────────────────────────────────
print("=" * 106)
print("第 1 部分：Tc 被库仑台阶钉住（λ₁ 在 Tc 处不等于 1）")
print("=" * 106)
LAM, WC = 0.540, 5.0 * WB
EIN = [(LAM * WB / 2.0, WB)]
SPL = two_atom(LAM, 0.95 * WB, 100.0 * WB)
tE = Tc(EIN, WC); T0 = tE * KB
x = WC / (2.0 * np.pi * T0)
print(f"  Tc_E = {tE:.9f} K    ω_c/(2πTc) = {x:.6f}    到整数距离 = {min(x%1, 1-x%1):.2e}")
print("  （memory `coulomb-ramp-makes-tc-non-c1`：距离 < 0.15 格 ⟹ 该点 Tc 不可微、导数不可报）")
print(f"\n  {'T [K]':>14}{'N_c':>6}{'λ1(E)−1':>15}{'λ1(split)−1':>15}{'Δ=split−E':>14}")
for dK in (-6e-3, -2e-3, 0.0, +2e-4, +2e-3, +6e-3):
    Tk = tE + dK
    lE = lam1(Tk * KB, EIN, WC); lS = lam1(Tk * KB, SPL, WC)
    Nc = int(np.floor(WC / (2.0 * np.pi * Tk * KB)))
    print(f"  {Tk:14.9f}{Nc:6d}{lE-1:+15.6e}{lS-1:+15.6e}{lS-lE:+14.3e}")
print("  ⟹ λ1(E)−1 在台阶左侧为 **+1.1e−4（非 0）**，跨台阶跳到 −4.9e−4；")
print("     而 Δ 跨台阶**稳定在 +2.0e−4** ⟹ Δ 是干净量，λ1−1 不是。\n")

# ── 第 2 部分：旧统计 vs 新统计，逐行重算 ──────────────────────────────
print("=" * 106)
print("第 2 部分：正对照带上，旧统计（λ1s−1）与校准统计（Δ）的违反计数")
print("=" * 106)
W1S, W2S = (0.99, 0.95, 0.9, 0.8, 0.6), (100., 50., 20., 10., 5., 2.)
print(f"  {'r':>5}{'λ':>6}{'dist':>7}{'E残差':>11}{'旧 max':>12}{'旧#':>5}"
      f"{'校准 maxΔ':>12}{'新#':>5}{'点数':>6}")
tot_old = tot_new = tot = 0
fired = []
for r in (5.0, 5.5, 6.5, 7.5):
    wc = r * WB
    for lam in (0.50, 0.54, 0.58, 0.62):
        tEr = Tc([(lam * WB / 2.0, WB)], wc)
        if not np.isfinite(tEr):
            continue
        T = tEr * KB
        xr = wc / (2.0 * np.pi * T); dist = min(xr % 1, 1 - xr % 1)
        if int(np.floor(xr)) > N:
            continue
        lE = lam1(T, [(lam * WB / 2.0, WB)], wc)
        old, new = [], []
        for w1f in W1S:
            for w2f in W2S:
                a = two_atom(lam, w1f * WB, w2f * WB)
                if a is None: continue
                ls = lam1(T, a, wc)
                old.append(ls - 1.0); new.append(ls - lE)
        no = sum(1 for v in old if v > 0); nn = sum(1 for v in new if v > 0)
        tot_old += no; tot_new += nn; tot += len(old)
        if nn > 0:
            fired.append((r, lam, nn, len(new), max(new)))
        print(f"  {r:5.1f}{lam:6.2f}{dist:7.3f}{lE-1:+11.2e}{max(old):+12.3e}{no:5d}"
              f"{max(new):+12.3e}{nn:5d}{len(old):6d}", flush=True)

print(f"\n  合计 {tot} 点：旧统计判违反 {tot_old} 个 ⟹ **校准后 {tot_new} 个**"
      f"（{100*(tot_old-tot_new)/max(tot_old,1):.0f}% 是偏置造成的假违反）")
print(f"  校准后仍亮的行：{fired}")
print()
print("=" * 106)
print("判词：")
print("  · 正对照**仍然亮**（非盲），但强度大幅下降：从「带上 12/323」降为「单点 8/30」。")
print("  · 判据区 0/3172 的方向是**保守**的（旧筛偏向高报），主结论不受损。")
print("  · Tc 层面：该点 Tc 被台阶钉住 ⟹ ΔTc ≡ 0 ⟹ 反例只在 **λ₁ 层面**成立，")
print("    Tc 层面的严格不等号**未证**（台阶是尖截断库仑的离散赝物，不是物理）。")
print("=" * 106)
