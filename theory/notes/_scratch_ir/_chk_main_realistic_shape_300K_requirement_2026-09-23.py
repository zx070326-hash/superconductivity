# -*- coding: utf-8 -*-
"""
主会话从零：300 K 判据的**两条细化**。

细化 (b)：上一版的 `(η/M)` 需求表建在 **Einstein 族**上，是「最有利谱形」下的必要条件。
  本脚本**不用 ×1.24 这个近似**，而是把**真实 α²F 的谱形标度上去直接重解**：
      α²F(ω) := c · α²F₀(ω/s)     ⟹   λ = c·λ₀ ,   ω̄₂ = s·ω̄₂₀        （两行代数，见下）
  ⟹ 谱形固定为真实形状，`(λ, ω̄₂)` 自由 ⟹ 与 Einstein 表**逐格可比**。
  好处：谱形罚变成**精确的、逐 λ 的**量，而不是在 H₃S 一点上测出的单个数字。

细化 (a)：Einstein 是否真是最有利谱形，在 `μ*>0` 时**未证**。
  本脚本给出**双侧夹**：
      Einstein 需求（= 若极值性成立，则是所有谱的**下确界**）  ≤  真实需求  ≤  真实谱形需求
  ⟹ 即使极值性不成立，**上侧（真实谱形那一列）仍然是无条件的**——
     它是一条**构造性**陈述（「存在一个真实谱形，用这么多 η/M 就能到 300 K」），
     不依赖任何极值性定理。⟹ 判据从此有一条**不欠债**的腿。

━━ 标度代数（analytic，两行）━━
令 ω = s·u：
    λ  = 2∫ cα²F₀(ω/s)/ω dω = 2c∫ α²F₀(u)/u du = c·λ₀
    λω̄₂² = 2∫ ω cα²F₀(ω/s) dω = 2cs²∫ u α²F₀(u) du = cs²·λ₀ω̄₂₀²
  ⟹ ω̄₂ = s·ω̄₂₀，与 c 无关；c 只调 λ。两个参数正交，无需迭代。
  （数值自检：重建后的 (λ, ω̄₂) 必须逐位命中目标，断言门 1e-10。）

━━ 方向审计 ━━
· **Einstein 列**：必要条件（**依赖**未证的 μ*>0 极值性）⟹ 只可用于排除，且带前提。
· **真实谱形列**：充分条件（构造性，**无前提**）⟹ 可用于「这么多就够」，方向与上一列相反。
  ⟹ 两列夹住真值，**判词必须分列写，不许混成一个数**。
· ω_c 约定逐 r 印出，不藏。

严格性：标度代数 analytic；Tc 求解 numerical-evidence（求解器已过 Leavens 1975 盲测 0.014%）。
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from scipy.linalg import eigh

KB = 8.617333262e-2
HB = 64.65415134                  # meV per sqrt(eV/Å²/amu)
C_INF = 0.1827262477
BASE = os.path.join(r"D:\physics\superconductivity", "theory", "data", "external",
                    "mishra_h3s_d3s", "Data-H3S-D3S-PH-EPW-Anharmonic")


def load_a2f(path):
    rows = []
    for ln in open(path, errors="replace"):
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        try:
            v = [float(x) for x in s.split()]
        except ValueError:
            continue
        if len(v) >= 2:
            rows.append(v[:2])
    a = np.array(rows)
    w, f = a[:, 0], np.clip(a[:, 1], 0.0, None)
    m = w > 1e-9
    return w[m], f[m]


quad = (lambda y, x: np.trapezoid(y, x)) if hasattr(np, "trapezoid") else (lambda y, x: np.trapz(y, x))


def moments(w, a):
    lam = 2.0 * quad(a / w, w)
    return lam, np.sqrt(2.0 * quad(a * w, w) / lam)


def G_of_spectrum(kmax, T, w, a):
    """向量化：G[k] = 2∫ ω α²F/(ω²+(2πTk)²) dω"""
    x2 = (2.0 * np.pi * T * np.arange(kmax, dtype=float)) ** 2      # (kmax,)
    den = w[None, :] ** 2 + x2[:, None]                             # (kmax, nω)
    return quad((2.0 * w * a)[None, :] / den, w)


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


def w2_for_Tc(T0_K, lam, shape, mus, r, N=400):
    """反解：固定 λ 与谱形，求达到 Tc=T0 所需的 ω̄₂。"""
    T = T0_K * KB
    if shape is None:                                   # Einstein：ω_E = ω̄₂
        def ev(w2):
            k = np.arange(2 * N + 2, dtype=float)
            G = lam * w2 ** 2 / (w2 ** 2 + (2.0 * np.pi * T * k) ** 2)
            return top_eig_G(T, G, N, mus, r * w2) - 1.0
    else:
        w0, a0, lam0, w20 = shape
        def ev(w2):
            s, c = w2 / w20, lam / lam0
            G = G_of_spectrum(2 * N + 2, T, s * w0, c * a0)
            return top_eig_G(T, G, N, mus, r * w2) - 1.0
    lo, hi = 1.0, 30000.0
    if ev(hi) < 0:
        return np.nan
    if ev(lo) > 0:
        return lo
    for _ in range(46):
        mid = np.sqrt(lo * hi)
        if ev(mid) < 0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


# ── 载入谱形 + 标度自检（断言门）──────────────────────────────────
SHAPES = {}
for tag, rel in [("H3S 非谐", "H3S-200GPa/EPW_anharmonic/H3S.a2f"),
                 ("H3S 谐波", "H3S-200GPa/EPW_harmonic/H3S.a2f"),
                 ("D3S 非谐", "D3S-200GPa/EPW_anharmonic/D3S.a2f")]:
    p = os.path.join(BASE, rel.replace("/", os.sep))
    if os.path.exists(p):
        w0, a0 = load_a2f(p)
        l0, v0 = moments(w0, a0)
        SHAPES[tag] = (w0, a0, l0, v0)

print("=" * 100)
print("标度自检（analytic 的两行代数，数值必须逐位命中）")
print("=" * 100)
for tag, (w0, a0, l0, v0) in SHAPES.items():
    s, c = 2.7183 / v0 * 100.0, 3.1416 / l0          # 任意目标
    lt, vt = moments(s * w0, c * a0)
    e1, e2 = abs(lt / 3.1416 - 1), abs(vt / 271.83 - 1)
    print(f"  {tag}:  原 (λ,ω̄₂) = ({l0:.4f}, {v0:.2f})  →  目标 (3.1416, 271.83)"
          f"  实得 ({lt:.6f}, {vt:.4f})   相对差 {e1:.2e} / {e2:.2e}")
    assert e1 < 1e-10 and e2 < 1e-10, "标度代数不自洽，停"
print("  ⟹ 全部 < 1e-10，标度代数确认\n")

MUS = 0.13
LAMS = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
H3S_ETA_M = 9.172

for T0 in (300.0, 203.0):
    print("=" * 100)
    print(f"目标 Tc = {T0:.0f} K   ·   μ* = {MUS}   ·   ω_c = 10·ω̄₂   ·   η/M 单位 eV/Å²/amu")
    print("=" * 100, flush=True)
    hdr = f"  {'λ':>5}{'Einstein':>20}"
    for tag in SHAPES:
        hdr += f"{tag:>22}"
    print(hdr)
    print(f"  {'':>5}{'ω̄₂[meV]':>10}{'η/M':>10}" + "".join(f"{'ω̄₂[meV]':>11}{'η/M':>11}" for _ in SHAPES))
    pen = {t: [] for t in SHAPES}
    for lam in LAMS:
        wE = w2_for_Tc(T0, lam, None, MUS, 10.0)
        emE = lam * wE ** 2 / HB ** 2
        line = f"  {lam:5.2f}{wE:10.1f}{emE:10.3f}"
        for tag in SHAPES:
            wS = w2_for_Tc(T0, lam, SHAPES[tag], MUS, 10.0)
            emS = lam * wS ** 2 / HB ** 2
            pen[tag].append(emS / emE)
            line += f"{wS:11.1f}{emS:11.3f}"
        print(line, flush=True)
    print(f"\n  谱形罚（真实需求 / Einstein 需求，按 η/M 计）:")
    for tag in SHAPES:
        p = np.array(pen[tag])
        print(f"     {tag}: {p.min():.3f}–{p.max():.3f}  （λ={LAMS[int(p.argmin())]:.1f}→{LAMS[int(p.argmax())]:.1f}）"
              f"   {'⚠ 有 <1 的格 = 极值性反例!' if p.min() < 1.0 else ''}")
    print()

# ── 关键单点：η/M = 9.172 固定时，各谱形达到 300 K 所需的最小 λ ──────────
print("=" * 100)
print(f"★ 固定 η/M = {H3S_ETA_M} （= H₃S 实测，纯电子量，软化声子不改它）：")
print("   达到 300 K 所需的 λ = ?   （等价地：ω̄₂ 要软化到多少）")
print("=" * 100, flush=True)
print(f"  {'谱形':>12}{'r=6':>22}{'r=10':>22}{'r=20':>22}")
print(f"  {'':>12}" + "".join(f"{'λ_需要':>11}{'ω̄₂[meV]':>11}" for _ in range(3)))
for tag in ["Einstein"] + list(SHAPES):
    sh = None if tag == "Einstein" else SHAPES[tag]
    line = f"  {tag:>12}"
    for r in (6.0, 10.0, 20.0):
        lo, hi = 0.8, 30.0
        g = lambda L: L * w2_for_Tc(300.0, L, sh, MUS, r) ** 2 / HB ** 2 - H3S_ETA_M
        if g(hi) > 0:                       # 大 λ 端仍不够
            line += f"{'—':>11}{'—':>11}"; continue
        if g(lo) < 0:
            line += f"{'<0.8':>11}{'—':>11}"; continue
        for _ in range(34):
            mid = 0.5 * (lo + hi)
            if g(mid) > 0: lo = mid
            else: hi = mid
        L = 0.5 * (lo + hi)
        line += f"{L:11.2f}{HB*np.sqrt(H3S_ETA_M/L):11.1f}"
    print(line, flush=True)
print(f"\n  对照：H₃S-200 非谐实测 λ = 1.721、ω̄₂ = 149.26 meV、Tc = 203 K")
print("  ⚠ 方向：Einstein 行**依赖未证的 μ*>0 极值性**（只可用于排除）；")
print("     三个真实谱形行是**构造性**的（存在这样的谱 ⟹ 这么多就够），**无前提**。")
print("=" * 100)
