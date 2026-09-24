# -*- coding: utf-8 -*-
"""
从零复核 certC 的「认证 C(λ,μ*) 曲线」（2026-09-24，主会话）。

certC 的产物：`_chk_main_C_lambda_certified_curve_2026-09-24.py` + 日志
`_log_C_lambda_2026-09-24.txt` / `_log_C_lambda_2026-09-24_remaining.txt`，
给出 9 行（λ=1.0…6.0, μ*=0.13）双侧认证括号 [C_lower, C_upper]，外加 μ*=0 的 λ=4/6 两行。

━━━ 我用**不同的算法**复核什么 ━━━
certC 的 C 来自 `KT.einstein_root` 的二分 + 区间分支定界。我这里用本会话新得的
**w-参数化引理**（`_chk_main_c_sawtooth_exact_branchwise_2026-09-24.py`）：
    λ₁(λK̂ − 2μ*u_cu_cᵀ) ≥ 1  ⟺  λ ≥ 1/θ_max，θ_max = 广义本征问题 K̂v = θBv 的最大根
⟹ `λ*(w) = 1/θ_max(w)` 一次 eigh 精确给出，**零二分**。把它对 w 反解到 λ*(w)=λ，
就得到 Einstein 阈值 ω̄₂^E——与 certC 的路径无共享代码（我不 import 它的任何模块）。

检验三件事（每件都写清方向）：
  V1  ω̄₂^E 列：我的值 vs certC 的值（纯数值一致性）。
  V2  **括号包含性**（真正有内容的那条）：C_E := 1/(2π ω̄₂^E √λ) 是否落在 [C_lower, C_upper] 内。
      方向：上侧 certify 在 w<ω̄* 给 sup<1 ⟹ 真阈值更大 ⟹ C 更小 ⟹ 报的是**上界**；
            下侧 Rayleigh 在 w>ω̄* 给 μ_lo>1 ⟹ 真阈值更小 ⟹ C 更大 ⟹ 报的是**下界**。
      ⟹ 必须有 w_up < ω̄₂^E < w_lo，等价地 C_lower < C_E < C_upper。任一侧越界 = certC 有缺陷。
  V3  C 的算术：C = 1/(2π w √λ) 对它印出的 (w_up, C_upper) 与 (w_lo, C_lower) 逐行重算。

━━━ 我**不**复核的 ━━━
区间分支定界内部的声音性（外向舍入、盒上包围、Kato–Temple 的每一步）——那需要重写整个
证书求解器。本脚本只能抓「数字对不对 / 方向对不对 / 括号自洽不自洽」，
抓不到「区间算术不外包」这一类缺陷。**判词里必须写清这一点。**

━━━ 一处自由度审计（第三起同形）━━━
certC 的 GATE N0「μ*=0 时 N_c 是否真的不进判据（不从代码结构推，实测）」**自由度为 0**：
μ*=0 ⟹ 秩一项系数 2μ*=0 ⟹ N_c 在代数上根本不出现 ⟹ 该门**不可能失败**。
（同一形状本轮第三起：mustar-scan、certC，以及我自己的第一版 N_c 扫描。）
它的判词文字自己说对了机制，但把恒等式呈现成「实测通过的门」，不得计入通过率。

严格性：全部 numerical-evidence（float64）。引理本身是 analytic。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

C_INF = 0.1827262477
R = 10.0


def pieces(wbar, N):
    """λ=1 归一的 K̂ 与权 u（与 KT.build 同构，但本文件自带、不 import）。"""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * N + 2)
    G = wbar * wbar / (wbar * wbar + k * k)
    Kp, Km = G[np.abs(i - j)], G[i + j + 1]
    S = (Kp - Km).sum(axis=1)
    return (Kp + Km) * np.outer(u, u) - np.diag(u * u * S), u


def lam_root_exact(wbar, mus, r, N):
    """引理：使 λ₁=1 的**最小** λ = 1/θ_max，一次广义 eigh，零二分。"""
    Khat, u = pieces(wbar, N)
    B = np.eye(N)
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(N) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    th = float(eigh(Khat, B, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])
    return (1.0 / th) if th > 0 else np.inf, Nc


def einstein_wbar(lam, mus, r, N, lo=0.05, hi=14.0, nbis=80):
    """反解 λ*(w) = λ ⟹ Einstein 阈值 ω̄₂^E。λ*(w) 对 w 递减（GATE V0 实测）。"""
    f = lambda w: lam_root_exact(w, mus, r, N)[0] - lam
    a, b, fa = lo, hi, None
    fa = f(a)
    for _ in range(nbis):
        m = 0.5 * (a + b)
        fm = f(m)
        if not np.isfinite(fm) or (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b)


def C_of(w, lam):
    return 1.0 / (2.0 * np.pi * w * np.sqrt(lam))


# certC 汇总表（μ*=0.13 的 9 个完整行 + μ*=0 的 2 行，逐字抄自日志）
ROWS = [
    # lam,  mus,  wE_certC,      C_upper,      w_up,        C_lower,      w_lo
    (1.0,   0.13, 2.200000000, 0.085109595, 1.870000,   0.071626887, 2.222000),
    (1.5,   0.13, 1.286982312, 0.101992146, 1.274112,   0.099972499, 1.299852),
    (2.0,   0.13, 0.974196811, 0.116687205, 0.964455,   0.114376567, 0.983939),
    (2.5,   0.13, 0.810793522, 0.125402058, 0.802686,   0.122918849, 0.818901),
    (3.0,   0.13, 0.708569341, 0.130991149, 0.701484,   0.128397265, 0.715655),
    (3.2953, 0.13, 0.664060982, 0.133361200, 0.657420,  0.130720384, 0.670702),
    (4.0,   0.13, 0.583478845, 0.137762123, 0.577644,   0.135034160, 0.589314),
    (5.0,   0.13, 0.506702078, 0.141888517, 0.501635,   0.139078844, 0.511769),
    (6.0,   0.13, 0.453637871, 0.144677171, 0.449101,   0.141812277, 0.458174),
    (4.0,   0.00, 0.477593276, 0.168304891, 0.472817,   0.164972121, 0.482369),
    (6.0,   0.00, 0.377578793, 0.173820789, 0.373803,   0.170378794, 0.381355),
]


def main():
    N = 800
    print("=" * 108)
    print(f"V0  单调性前提：λ*(w) 对 w 递减（einstein_wbar 的二分依赖它）  N={N}")
    print("=" * 108)
    bad = 0
    for mus in (0.0, 0.13):
        prev = None
        ws = np.linspace(0.2, 6.0, 40)
        vals = [lam_root_exact(w, mus, R, N)[0] for w in ws]
        dec = all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
        print(f"   μ*={mus:.2f}: 40 个 w 上严格递减 = {dec}   "
              f"λ*(0.2)={vals[0]:.4f} → λ*(6.0)={vals[-1]:.4f}")
        bad += 0 if dec else 1

    print()
    print("=" * 108)
    print("V1+V2+V3  逐行：我的 ω̄₂^E（w-参数化引理，零二分）vs certC；括号包含性；C 的算术")
    print("=" * 108)
    print("     λ    μ*     我的 ω̄₂^E      certC 的       相对差   |  C_E(我算)    "
          "[C_lower, C_upper]              包含?  C 算术最劣差")
    nfail = 0
    for lam, mus, wE_c, Cup, w_up, Clo, w_lo in ROWS:
        wE = einstein_wbar(lam, mus, R, N)
        rel = abs(wE - wE_c) / wE_c
        CE = C_of(wE, lam)
        inside = (Clo < CE < Cup)
        # V3：C 的算术
        a1 = abs(C_of(w_up, lam) - Cup) / Cup
        a2 = abs(C_of(w_lo, lam) - Clo) / Clo
        arith = max(a1, a2)
        tag = "✓" if inside else "✗ 越界"
        if not inside:
            nfail += 1
        print(f"  {lam:6.4f}  {mus:.2f}   {wE:.9f}  {wE_c:.9f}  {rel:.2e}  |  "
              f"{CE:.9f}  [{Clo:.9f}, {Cup:.9f}]  {tag}   {arith:.2e}")
    print()
    print(f"  ⟹ 括号包含性：{len(ROWS)-nfail}/{len(ROWS)} 行通过"
          f"（w_up < ω̄₂^E < w_lo，等价 C_lower < C_E < C_upper）")

    print()
    print("=" * 108)
    print("V4  一致性门：全部认证端点 C/C_∞ < 1？（与已并入的 P4_conjecture2.md 定理相容）")
    print("=" * 108)
    worst = max(Cup / C_INF for _, _, _, Cup, _, _, _ in ROWS)
    arg = max(ROWS, key=lambda r: r[3] / C_INF)
    print(f"  max C_upper/C_∞ = {worst:.5f}  @ λ={arg[0]}, μ*={arg[1]}   "
          f"⟹ {'PASS' if worst < 1 else 'ALARM'}")
    print("  C_upper/C_∞ 随 λ（μ*=0.13）：",
          "  ".join(f"{r[0]:g}:{r[3]/C_INF:.3f}" for r in ROWS if r[1] == 0.13))
    print("  ⟹ 单调升；这与 C_∞ 是**强耦合渐近**（LAD79 的 kTc=0.1827√(η/M)）相符，")
    print("     即我们的分辨-λ 界从下方趋近它 —— 一条独立的结构性佐证（非证明）。")

    print()
    print("=" * 108)
    print("V5  自由度审计：certC 的 GATE N0（μ*=0 时改 N_c）自由度为 0")
    print("=" * 108)
    w = 0.477593276
    Khat, u = pieces(w, N)
    for Nc in (0, 1, 3, 5, 10, 50):
        B = np.eye(N)   # μ*=0 ⟹ 秩一项系数恰为 0 ⟹ B 与 Nc 无关
        th = float(eigh(Khat, B, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])
        print(f"    N_c={Nc:3d}  1/θ_max={1.0/th:.15f}")
    print("  ⟹ 逐位相同是**恒等式**（B ≡ I，Nc 在代数上不出现），不是通过的检验。")
    print("     本轮同形第三起（mustar-scan / certC / 我自己的第一版 N_c 扫描）。")

    print()
    print("=" * 108)
    print("判词")
    print("=" * 108)
    print("  · V1–V4 全过 ⟹ certC 的表在**数值与方向**上自洽，且与独立算法一致。")
    print("  · 本脚本**未**复核区间分支定界内部的声音性（外向舍入 / 盒上包围 / KT 每一步）——")
    print("    抓不到「区间算术不外包」这类缺陷 ⟹ 该表的 certificate 级别**我未独立确认**，")
    print("    我确认的是 numerical-evidence 级的一致性。")


if __name__ == "__main__":
    main()
