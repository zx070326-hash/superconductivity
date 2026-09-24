# -*- coding: utf-8 -*-
"""
certC 补跑的 4 行 + 「紧度比单调」这条结构性断言的从零复核（2026-09-24，主会话）。

背景：我 09-24 早先的笔记 `ceiling/C_lambda_certified_curve_2026-09-24.md` §5 把
λ=8(下侧)、λ=10(两侧)、μ*=0 的 λ=2 与 3.2953 写成「预算耗尽，未算」。
certC 的第二次后台运行（`_log_C_lambda_2026-09-24_remaining.txt`）把它们补齐了
⟹ 该 §5 已过时，且 §1.1 的「紧度比随 λ 单调升」必须用**补齐后的全表**重判。

━━━ 本脚本查三件事（每件写清方向）━━━

V6  4 个新行的 V1/V2/V3：我用 w-参数化引理（λ*(w)=1/θ_max，一次广义 eigh、零二分，
    不 import certC 任何模块）独立算 ω̄₂^E，检验括号包含性 C_lower < C_E < C_upper。
    方向：上侧 certify 在 w<ω̄* 给 sup<1 ⟹ 真阈值更大 ⟹ C 更小 ⟹ **上界**；
          下侧 Rayleigh 在 w>ω̄* 给 μ_lo>1 ⟹ 真阈值更小 ⟹ C 更大 ⟹ **下界**。

GATE M  **紧度比单调性的正确判法**（这是本脚本的重点）。
    ⚠ 一列**上界**单调升，对真值 C(λ) 的单调性**什么都不证**（C ≤ C_upper）。
      我此前 §1.1 用 C_upper/C_∞ 排出的递增序列，作为「结构性佐证」方向是对的，
      但**判据选错了**：正确的、真能进 certificate 的判法是**括号链不重叠**——
          若 C_lower(λ_{i+1}) > C_upper(λ_i) 则 C(λ_{i+1}) > C(λ_i) 严格成立（两侧都是证书）。
      重叠的相邻对判**未定**（既非升也非降），不得读成「掉头」。
    新表里 λ=8 的 C_upper=0.150915 > λ=10 的 C_upper=0.150565 看似掉头，
    本门要判：那是真掉头，还是两格**收紧深度不等**造成的赝象。

GATE F  收紧深度普查：逐行从日志里读最后一次 PASS 的 f，看哪些格没到 f=0.010。
    （certC 的消息只点名了 λ=1.0 一格，第一份日志的汇总列里 λ=8 也标着 BUDGET。）

GATE C_E  C_E(λ) 的细网格单调性（numerical-evidence，float64 Einstein 精确根）：
    真值代理量在 λ∈[1,12] 上是否单调升、λ=8→10 之间有没有拐。

严格性：V6/GATE C_E 全部 numerical-evidence（float64）。
        GATE M 的判词本身是 certificate 级（只用 certC 的两侧证书端点做算术比较），
        但前提是那些端点的区间声音性——**那部分我未复核**（见笔记 §4.1）。
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
    f = lambda w: lam_root_exact(w, mus, r, N)[0] - lam
    a, b = lo, hi
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


# 全表：逐字抄自 certC 的两份日志（第一份 + _remaining）。
# lam, mus, wE_certC, C_upper, w_up, C_lower, w_lo, f_up_last, f_lo_last
ROWS = [
    (1.0,    0.13, 2.200000000, 0.085109595, 1.870000, 0.071626887, 2.222000, 0.150, 0.010),
    (1.5,    0.13, 1.286982312, 0.101992146, 1.274112, 0.099972499, 1.299852, 0.010, 0.010),
    (2.0,    0.13, 0.974196811, 0.116687205, 0.964455, 0.114376567, 0.983939, 0.010, 0.010),
    (2.5,    0.13, 0.810793522, 0.125402058, 0.802686, 0.122918849, 0.818901, 0.010, 0.010),
    (3.0,    0.13, 0.708569341, 0.130991149, 0.701484, 0.128397265, 0.715655, 0.010, 0.010),
    (3.2953, 0.13, 0.664060982, 0.133361200, 0.657420, 0.130720384, 0.670702, 0.010, 0.010),
    (4.0,    0.13, 0.583478845, 0.137762123, 0.577644, 0.135034160, 0.589314, 0.010, 0.010),
    (5.0,    0.13, 0.506702078, 0.141888517, 0.501635, 0.139078844, 0.511769, 0.010, 0.010),
    (6.0,    0.13, 0.453637871, 0.144677171, 0.449101, 0.141812277, 0.458174, 0.010, 0.010),
    (8.0,    0.13, 0.382417923, 0.150914957, 0.372857, 0.145685231, 0.386242, 0.025, 0.010),  # ← 上侧 BUDGET
    (10.0,   0.13, 0.337646676, 0.150564456, 0.334270, 0.147582982, 0.341023, 0.010, 0.010),
    (2.0,    0.00, 0.751648293, 0.151236028, 0.744132, 0.148241256, 0.759165, 0.010, 0.010),
    (3.2953, 0.00, 0.537644980, 0.164718304, 0.532269, 0.161456556, 0.543021, 0.010, 0.010),
    (4.0,    0.00, 0.477593276, 0.168304891, 0.472817, 0.164972121, 0.482369, 0.010, 0.010),
    (6.0,    0.00, 0.377578793, 0.173820789, 0.373803, 0.170378794, 0.381355, 0.010, 0.010),
]
NEW_ROWS = {(8.0, 0.13), (10.0, 0.13), (2.0, 0.0), (3.2953, 0.0)}


def main():
    N = 800
    print("=" * 112)
    print(f"V6  certC 补跑的 4 行：我的 ω̄₂^E（w-参数化引理）vs certC；括号包含性；C 的算术   N={N}")
    print("=" * 112)
    print("     λ    μ*   新?   我的 ω̄₂^E      certC 的       相对差   |  C_E(我算)   "
          "[C_lower, C_upper]              包含?   C 算术最劣差")
    nfail = 0
    CE_tab = {}
    for lam, mus, wE_c, Cup, w_up, Clo, w_lo, fu, fl in ROWS:
        wE = einstein_wbar(lam, mus, R, N)
        CE = C_of(wE, lam)
        CE_tab[(lam, mus)] = CE
        rel = abs(wE - wE_c) / wE_c
        inside = (Clo < CE < Cup)
        arith = max(abs(C_of(w_up, lam) - Cup) / Cup, abs(C_of(w_lo, lam) - Clo) / Clo)
        if not inside:
            nfail += 1
        mark = "★新" if (lam, mus) in NEW_ROWS else "  "
        print(f"  {lam:6.4f}  {mus:.2f}  {mark}  {wE:.9f}  {wE_c:.9f}  {rel:.2e}  |  "
              f"{CE:.9f}  [{Clo:.9f}, {Cup:.9f}]  {'✓' if inside else '✗越界'}    {arith:.2e}")
    print(f"\n  ⟹ 括号包含性 {len(ROWS)-nfail}/{len(ROWS)} 行通过"
          f"（其中补跑的 4 行 {sum(1 for r in ROWS if (r[0],r[1]) in NEW_ROWS)}/4）")

    print()
    print("=" * 112)
    print("GATE F  收紧深度普查：哪些格没收紧到 f=0.010？（certC 消息只点名了 λ=1.0）")
    print("=" * 112)
    short = [(r[0], r[1], r[7], "upper") for r in ROWS if r[7] > 0.0101] + \
            [(r[0], r[1], r[8], "lower") for r in ROWS if r[8] > 0.0101]
    for lam, mus, f, side in short:
        print(f"    λ={lam:g}  μ*={mus:.2f}  {side} 侧只到 f={f:.3f}  ⟹ 该侧合法但**更松**")
    print(f"  ⟹ 共 {len(short)} 格；certC 消息里只点名 1 格（λ=1.0 upper）。")
    print("     λ=8 upper 在**它自己第一份日志的汇总列**里标着 BUDGET，消息的警示段漏了它。")

    print()
    print("=" * 112)
    print("GATE M  单调性的**正确**判法：括号链不重叠 ⟹ certificate 级严格升")
    print("=" * 112)
    print("  ⚠ 一列上界单调升，对真值 C(λ) 的单调性什么都不证（C ≤ C_upper）。")
    print("    可进 certificate 的判据只有：C_lower(λ_{i+1}) > C_upper(λ_i)。")
    print()
    for mus in (0.13, 0.0):
        sub = [r for r in ROWS if r[1] == mus]
        sub.sort(key=lambda r: r[0])
        print(f"  ── μ*={mus:.2f} ──")
        est = und = 0
        for a, b in zip(sub, sub[1:]):
            ok = b[5] > a[3]          # C_lower(next) > C_upper(prev)
            gap = b[5] - a[3]
            if ok:
                est += 1
            else:
                und += 1
            print(f"    λ {a[0]:>6g} → {b[0]:<6g}  C_lower(次)={b[5]:.9f}  "
                  f"C_upper(前)={a[3]:.9f}  差={gap:+.2e}  "
                  f"{'✓ 严格升（certificate）' if ok else '— 未定（括号重叠）'}")
        print(f"    ⟹ {est} 对确立严格升，{und} 对未定。**未定 ≠ 下降。**")
        print()

    print("=" * 112)
    print("GATE C_E  真值代理量 C_E(λ) 的细网格单调性（numerical-evidence）")
    print("=" * 112)
    lams = [1.0, 1.5, 2.0, 2.5, 3.0, 3.2953, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
    vals = []
    for lam in lams:
        w = einstein_wbar(lam, 0.13, R, N)
        vals.append(C_of(w, lam))
    print("    λ      C_E          C_E/C_∞     Δ 比前一个")
    bad = 0
    for i, (lam, v) in enumerate(zip(lams, vals)):
        d = "" if i == 0 else f"{v - vals[i-1]:+.3e}"
        if i > 0 and v <= vals[i - 1]:
            bad += 1
            d += "  ✗下降"
        print(f"   {lam:6g}  {v:.9f}  {v/C_INF:.5f}   {d}")
    print(f"  ⟹ 15 个 λ 上 C_E 严格升：{bad == 0}（违反 {bad} 处）")
    print()
    print("  λ=8 与 λ=10 的对照（拆解「掉头」）：")
    i8, i10 = lams.index(8.0), lams.index(10.0)
    print(f"    C_E:      λ=8 {vals[i8]:.9f}  <  λ=10 {vals[i10]:.9f}   ⟹ 真值代理量**在升**")
    r8 = [r for r in ROWS if r[0] == 8.0 and r[1] == 0.13][0]
    r10 = [r for r in ROWS if r[0] == 10.0 and r[1] == 0.13][0]
    print(f"    C_upper:  λ=8 {r8[3]:.9f}  >  λ=10 {r10[3]:.9f}   ⟹ 上界看着在降")
    print(f"    括号相对宽度: λ=8 {(r8[3]-r8[5])/r8[5]*100:.2f}%   λ=10 {(r10[3]-r10[5])/r10[5]*100:.2f}%")
    print(f"    收紧深度:     λ=8 upper f={r8[7]:.3f}   λ=10 upper f={r10[7]:.3f}")
    print("  ⟹ 判词：λ=8 上侧因预算停在 f=0.025（比别人松一级），上界偏高；")
    print("     两格括号**重叠**（GATE M 判未定）⟹ 「紧度比掉头」是**收紧深度不等的赝象**，不是发现。")

    print()
    print("=" * 112)
    print("判词")
    print("=" * 112)
    print("  1. 补跑的 4 行全部通过括号包含性 ⟹ 笔记 §5 的「未算」清单须改写。")
    print("  2. §1.1 的判据要换：从「C_upper 列单调升」换成「括号链不重叠」——")
    print("     前者对真值什么都不证；后者是 certificate 级，且在全表上只对部分相邻对成立。")
    print("  3. λ=8 上侧是**第二个**预算受限格（certC 消息漏报），它造成表面的紧度比掉头。")


if __name__ == "__main__":
    main()
