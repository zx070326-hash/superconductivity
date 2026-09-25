# -*- coding: utf-8 -*-
"""预算 — 路线 1+2：把**认证的** C(λ,μ*) 代进真实氢化物，看天花板从 2.05× 收到多少。

问题：仓库已有的 `Tc ≤ C_∞√(λ⟨ω²⟩)`（C_∞=0.1827262477，analytic，λ/μ* 全盲）在
H3S 上松 2.05×（𝔣=0.489）⟹ 当筛子只能剪掉真 Tc < target/2 的材料，很钝。
2026-09-24 认证出的 `C(λ,μ*)` 分辨了 λ 与 μ*。**代进去能收到多少？** 这是预算问题：
若收不动（𝔣 仍 ~0.5），路线 1 的增量为 ∞，不该派活。

方向审计（每一步都朝「天花板变大」退化 ⟹ 「Tc ≤ 天花板」是 PASS 才有意义）：
  (a) **C 对 λ 单调升**（GATE M，certificate，11 点）⟹ 对 λ 取**上邻格**的 `C_upper`
      是 C(λ_true) 的合法上界。**不插值**——插值的是上界列，不是真值。
  (b) **C 对 μ\* 单调降**（下面 GATE μ 实测）⟹ 用比材料真 μ\* **小**的那档是上界。
      `C(λ,0)` 无条件成立（μ\*≥0）；`C(λ,0.13)` 条件于「该材料 μ\*≥0.13」。
  (c) ω̄₂ 缺失时用 ω_log 代（Jensen: ω_log ≤ ω̄₂）⟹ M2 **下估** ⟹ 天花板**下估**
      ⟹ 方向**相反**，只支持「无违反」的 PASS，不支持「违反」的 FAIL。该列单独标记。

严格性：C 的两列是 `certificate`（转引 C_lambda_certified_curve_2026-09-24.md §1 的端点）；
λ/ω̄₂/Tc 来自外部论文 ⟹ 整行结论 `numerical-evidence`。
"""
import os
import sys
import csv
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import math

import numpy as np
from scipy.linalg import eigh

KB = 0.086173332621452  # meV/K
C_INF = 0.1827262477

# 认证端点，逐字抄自 C_lambda_certified_curve_2026-09-24.md §1 表（C_upper 列）
C_UP = {
    0.00: [(2.0, 0.151236028), (3.2953, 0.164718304), (4.0, 0.168304891), (6.0, 0.173820789)],
    0.13: [(1.0, 0.085109595), (1.5, 0.101992146), (2.0, 0.116687205), (2.5, 0.125402058),
           (3.0, 0.130991149), (3.2953, 0.133361200), (4.0, 0.137762123), (5.0, 0.141888517),
           (6.0, 0.144677171), (8.0, 0.148030265), (10.0, 0.150564456)],
    0.20: [(2.0, 0.104372372), (3.2953, 0.122242127), (4.0, 0.126754407), (6.0, 0.133713798)],
}


def c_upper(lam, mus):
    """C(λ,μ*) 的**声音**上界：取 λ 的上邻认证格（依赖 C 对 λ 单调升）。
    λ 超出表的右端 ⟹ 回退到 C_∞（无条件 analytic 上界）。"""
    for L, c in C_UP[mus]:
        if lam <= L:
            return c, L
    return C_INF, math.inf


# ─────────── GATE μ：C 对 μ* 单调降？（(b) 的前提，我自己测，不转引） ───────────
R = 10.0


def theta_max(w, mus, n=400):
    Nc = int(math.floor(R * w))
    k = np.arange(2 * n + 2)
    G = w * w / (w * w + k * k)
    i, j = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    u = 1.0 / np.sqrt(2.0 * np.arange(n) + 1.0)
    S = G[0] + 2.0 * (np.cumsum(G[:n]) - G[0])
    K = (G[np.abs(i - j)] + G[i + j + 1]) * np.outer(u, u) - np.diag(u * u * S)
    uc = np.where(np.arange(n) < Nc, u, 0.0)
    B = np.eye(n) + 2.0 * mus * np.outer(uc, uc)
    return float(eigh(K, B, eigvals_only=True)[-1])


def w_thresh(lam, mus, lo=0.05, hi=8.0):
    """解 λ·θ_max(w) = 1 的 w = ω̄₂^E（无量纲阈值）；C = 1/(2π ω̄₂ √λ)。
    ⚠ 先扫变号再二分（⌊r·w⌋ 台阶，见 coulomb-saturation-creates-a-pseudo-root）。"""
    ws = np.linspace(lo, hi, 601)
    vs = np.array([lam * theta_max(w, mus) - 1.0 for w in ws])
    sc = np.where(np.diff(np.sign(vs)) != 0)[0]
    if len(sc) == 0:
        return float("nan"), 0
    s = sc[-1]                      # θ_max 随 w 升 ⟹ 取最大变号（C 最小 ⟹ 上界方向保守）
    a, b, va = ws[s], ws[s + 1], vs[s]
    for _ in range(70):
        m = 0.5 * (a + b)
        if (lam * theta_max(m, mus) - 1.0) * va > 0:
            a = m
        else:
            b = m
    return 0.5 * (a + b), len(sc)


print("=" * 108)
print("GATE μ — C(λ,μ*) 对 μ* 单调降？（(b) 的前提；float64，numerical-evidence）")
print(f"{'λ':>6} | " + "".join(f"{f'μ*={m:.2f}':>14}" for m in (0.00, 0.05, 0.10, 0.13, 0.16, 0.20)) + "   单调?")
mono_ok = True
for lam in (1.5, 2.0, 3.2953, 4.0, 6.0):
    cs = []
    for mus in (0.00, 0.05, 0.10, 0.13, 0.16, 0.20):
        w, _ = w_thresh(lam, mus)
        cs.append(1.0 / (2 * math.pi * w * math.sqrt(lam)))
    ok = all(cs[i] > cs[i + 1] for i in range(len(cs) - 1))
    mono_ok &= ok
    print(f"{lam:6.3f} | " + "".join(f"{c:14.9f}" for c in cs) + f"   {'降✓' if ok else 'XX'}")
print(f"⟹ 五行全部严格递减：{mono_ok}   （成立 ⟹ 用低档 μ* 的 C 是上界，方向 (b) 合法）")

# 与认证表对照（我自己的 float64 应落在 [C_lower, C_upper] 内）
print("\n对照认证端点（我的 float64 应 ≤ C_upper 且 ≥ C_lower）：")
CHK = [(2.0, 0.13, 0.114376567, 0.116687205), (4.0, 0.13, 0.135034160, 0.137762123),
       (4.0, 0.00, 0.164972121, 0.168304891), (4.0, 0.20, 0.124244418, 0.126754407)]
for lam, mus, clo, chi in CHK:
    w, nsc = w_thresh(lam, mus)
    c = 1.0 / (2 * math.pi * w * math.sqrt(lam))
    print(f"  λ={lam:<7.4f} μ*={mus:.2f}  我={c:.9f}  括号=[{clo:.9f},{chi:.9f}]  "
          f"{'在内✓' if clo <= c <= chi else 'XX 越界'}  [变号{nsc}次]")

# ─────────── 路线 1：真实氢化物逐行 ───────────
rows = list(csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             "..", "..", "data", "hydrides_theory.csv"),
                                encoding="utf-8")))
# 实验 Tc 锚（K）——用于 𝔣；压力失配的行标注
TC_EXP = {("H3S", "H"): (203, "Drozdov2015 @155GPa"), ("D3S", "D"): (150, "Drozdov2015"),
          ("LaH10", "H"): (250, "Drozdov2019 @170GPa"), ("LaD10", "D"): (180, "Drozdov2019"),
          ("YH6", "H"): (224, "Troyan2021 @166GPa"), ("YH9", "H"): (243, "Kong2021 @201GPa"),
          ("CaH6", "H"): (215, "Ma/Li2022 @172GPa")}


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


print("\n" + "=" * 108)
print("路线 1 — 真实氢化物：C_∞（旧）vs 认证 C(λ,μ*)（新）")
print(f"{'material':<8}{'iso':<4}{'P':>5}{'treat':<10}{'λ':>7}{'ω̄₂':>8}{'src':>4} | "
      f"{'M2(meV²)':>10}{'ceil_∞':>9}{'ceil_C0':>9}{'ceil_C13':>10} | "
      f"{'Tc_exp':>7}{'𝔣_∞':>7}{'𝔣_C0':>7}{'𝔣_C13':>7}")
out = []
for d in rows:
    lam = fnum(d["lambda"])
    w2 = fnum(d["omega_2_meV"])
    src = "ω̄₂"
    if w2 is None:                       # 回退 ω_log ⟹ M2 下估 ⟹ 天花板下估（方向相反！）
        w2 = fnum(d["omega_log_meV"])
        src = "ωlog"
    if lam is None or w2 is None:
        continue
    M2 = lam * w2 * w2
    ceil_inf = C_INF * math.sqrt(M2) / KB
    c0, _ = c_upper(lam, 0.00)
    c13, _ = c_upper(lam, 0.13)
    ceil0 = c0 * math.sqrt(M2) / KB
    ceil13 = c13 * math.sqrt(M2) / KB
    key = (d["material"], d["isotope"].split()[0])
    tce = TC_EXP.get(key, (None, ""))[0]
    f = lambda c: (tce / c) if tce else float("nan")
    print(f"{d['material']:<8}{d['isotope'].split()[0]:<4}{d['pressure_GPa']:>5}"
          f"{d['treatment'][:9]:<10}{lam:7.3f}{w2:8.2f}{src:>5} | "
          f"{M2:10.0f}{ceil_inf:9.1f}{ceil0:9.1f}{ceil13:10.1f} | "
          f"{(tce if tce else 0):7d}{f(ceil_inf):7.3f}{f(ceil0):7.3f}{f(ceil13):7.3f}"
          + ("" if src == "ω̄₂" else "   ⚠ωlog代，天花板下估"))
    out.append((d["material"], src, tce, ceil_inf, ceil0, ceil13))

print("\n" + "-" * 108)
have = [o for o in out if o[2] and o[1] == "ω̄₂"]
if have:
    fi = [o[2] / o[3] for o in have]
    f0 = [o[2] / o[4] for o in have]
    f13 = [o[2] / o[5] for o in have]
    print(f"只看 ω̄₂ 真值行（n={len(have)}）：")
    print(f"  𝔣_∞   ∈ [{min(fi):.3f},{max(fi):.3f}]  中位 {np.median(fi):.3f}  ⟹ 松 {1/np.median(fi):.2f}×")
    print(f"  𝔣_C0  ∈ [{min(f0):.3f},{max(f0):.3f}]  中位 {np.median(f0):.3f}  ⟹ 松 {1/np.median(f0):.2f}×  (无条件)")
    print(f"  𝔣_C13 ∈ [{min(f13):.3f},{max(f13):.3f}]  中位 {np.median(f13):.3f}  ⟹ 松 {1/np.median(f13):.2f}×  (条件 μ*≥0.13)")
    print(f"\n  **收紧倍数**：C_∞ → C(λ,0) = {np.median(f0)/np.median(fi):.3f}×，"
          f"C_∞ → C(λ,0.13) = {np.median(f13)/np.median(fi):.3f}×")
print("\n⚠ 0/N 违反 = 天花板与实验相容；不构成「天花板紧」。紧度看 𝔣。")
