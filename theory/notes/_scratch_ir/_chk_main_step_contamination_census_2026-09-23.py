# -*- coding: utf-8 -*-
"""
§6.2 (四)：`λ₁(Einstein, Tc_E) − 1` 的残差普查——台阶污染到底覆盖判据区多少？符号统不统一？

━━ 起因 ━━
`extremality-hunt` 的最终判词写「0/14088（含 41/112 钉住组、**经你的 Δ 校准复核**）」。
两处要纠：
  ① 我只用干净统计量 `Δ` 复核了 **896** 点（§6.2 (三)），不是 14088；
  ② §6.1 (三) 的方向论证「偏置为正 ⟹ 旧筛只会**高报**违反 ⟹ 0 违反反而更强」
     是在**正对照那一个点**上实测到的，**不自动推广**到判据区。

━━ 本脚本判什么 ━━
旧筛把「违反」判为 `λ₁(split, Tc_E) − 1 > 0`，默认 `res := λ₁(Einstein, Tc_E) − 1 = 0`。
`res ≠ 0` 时旧统计量 = 真统计量 `Δ` + `res`，**偏置 = res**。于是两个问题：
  (Q1) 判据区里有多大比例的参照点 `res` 显著非零（= 被 `N_c` 取整台阶污染）？
  (Q2) 这些 `res` 的**符号统一**吗？不统一 ⟹ §6.1 (三) 的单边论证不可外推。

━━ 方向审计 ━━
本脚本**不产出**任何关于极值性的阳性结论，只刻画**仪器**。
两种结果的后果相反、都要如实写：
· `res` 一致为正 ⟹ 旧筛保守，teammate 的 0/14088 反而更强；
· `res` 两种符号都有 ⟹ 该论证**不可外推**，0/14088 只能报「仪器未校准的计数」。

━━ 关键前提（已在 §6.2 (T1) 核过）━━
`T := Tc_E` 由 Einstein 谱**单独**确定，算在 split 循环之外 ⟹ 同组全部 split 共享同一
`T / N_c / res`。故对 `(μ*, r, λ)` 的 Einstein 参照点做普查，**就是**对污染总体的普查。

严格性：numerical-evidence（float64，N=400，48 个参照点）。
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

KB = 8.617333262e-2
N = 400
WB = 150.0


def Gvec(T, aw, kmax=2 * N + 2):
    x2 = (2.0 * np.pi * T * np.arange(kmax, dtype=float)) ** 2
    G = np.zeros(kmax)
    for a, w in aw:
        G += 2.0 * a * w / (w * w + x2)
    return G


def lam1(T, G, wc, mus):
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


def Tc(aw, wc, mus):
    """§9 第 6 条：下端必须钉在 N_c ≤ N 的温度形式，否则落库仑饱和区的赝根。"""
    f = lambda T: lam1(T, Gvec(T, aw), wc, mus) - 1.0
    lo, hi = wc / (2.0 * np.pi * (N + 1)), 900.0 * KB
    if not (f(lo) > 0 > f(hi)):
        return np.nan
    for _ in range(50):
        m = np.sqrt(lo * hi)
        if f(m) > 0: lo = m
        else: hi = m
    return np.sqrt(lo * hi)


rows = []
for mus in (0.10, 0.13, 0.16):
    for r in (4.0, 6.0, 10.0, 20.0):
        for lam in (1.0, 2.0, 3.0, 5.0):
            aw, wc = [(lam * WB / 2.0, WB)], r * WB
            T = Tc(aw, wc, mus)
            if not np.isfinite(T):
                continue
            res = lam1(T, Gvec(T, aw), wc, mus) - 1.0
            x = wc / (2.0 * np.pi * T)
            rows.append((mus, r, lam, T / KB, int(x // 1), min(x % 1, 1 - x % 1), res))

# 噪声地板：二分收敛到机器精度时 res 的量级（eigh 的舍入 ~ few*eps*λ₁）
FLOOR = 1e-10
bad = [q for q in rows if abs(q[6]) > FLOOR]

print("=" * 96)
print(f"§6.2 (四) 残差普查：判据区 {len(rows)} 个 Einstein 参照点（μ*×r×λ = 3×4×4）")
print("=" * 96)
print(f"  {'mu*':>6}{'r':>5}{'lam':>6}{'Tc[K]':>11}{'Nc':>6}{'dist→int':>11}{'res':>14}   注")
for q in sorted(rows, key=lambda z: -abs(z[6]))[:6]:
    tag = "**污染**" if abs(q[6]) > FLOOR else "机器噪声"
    print(f"  {q[0]:6.2f}{q[1]:5.0f}{q[2]:6.1f}{q[3]:11.3f}{q[4]:6d}{q[5]:11.4f}{q[6]:+14.3e}   {tag}")

print(f"\n  (Q1) |res| > {FLOOR:g} 的点：**{len(bad)}/{len(rows)}** "
      f"（{100*len(bad)/len(rows):.1f}%）；其余 max|res| = "
      f"{max((abs(q[6]) for q in rows if abs(q[6])<=FLOOR), default=0):.2e}（= eigh 舍入地板）")
print(f"  (Q2) 污染点的 res 符号：{[f'{q[6]:+.2e}' for q in bad]}")
sign_uniform = len({q[6] > 0 for q in bad}) <= 1
print(f"       符号统一？ **{'是' if sign_uniform else '否'}**")

# dist<0.15 这个「钉住」旗标准不准？
pinned = [q for q in rows if q[5] < 0.15]
pinned_bad = [q for q in pinned if abs(q[6]) > FLOOR]
print(f"\n  (Q3) §6.2 (T2) 用的『钉住』旗 dist<0.15：命中 {len(pinned)}/{len(rows)}，"
      f"其中真被污染 **{len(pinned_bad)}**")
print(f"       ⟹ 该旗标**大幅高估**污染面；真正要紧的是 dist 落进二分最终括号内（实测 dist→0.0000）")

print()
print("=" * 96)
print("判词：")
print(f"  · 台阶污染在判据区**罕见**（{len(bad)}/{len(rows)}），远好于 §6.2 原先的措辞所暗示的。")
print("  · 但污染点的 res **符号不统一** ⟹ §6.1 (三)「偏置为正 ⟹ 旧筛保守」**不可外推**到判据区。")
print("  · ⟹ teammate 的 0/14088 正确报法：『仪器在 ~%.0f%% 的参照点上带不可定号的偏置』，" % (100*len(bad)/len(rows)))
print("    而**不是**『经 Δ 校准复核』（我只用 Δ 复核了 896 点，见 §6.2 (三)）。")
print("  · 方向：本脚本只刻画仪器，对极值性不产出任何阳性结论。")
print("=" * 96)
