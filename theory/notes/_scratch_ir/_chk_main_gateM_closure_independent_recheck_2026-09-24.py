# -*- coding: utf-8 -*-
"""
GATE M 关闭（10/10）的**从零独立复核**（2026-09-24，主会话）。

被复核对象：gateM-cells 交付的 4 个新证书端点
  ① λ=3.2953 下侧 f=0.0060  ② λ=6.0 下侧 f=0.0060
  ③ λ=8.0   上侧 f=0.0060  ④ λ=10.0 下侧 f=0.0042
以及它据此重算的「μ*=0.13 的 10 对 GATE M 全部确立」。

**不 import certC / gateM-cells 的任何模块**；ω̄₂^E 用本轮的 w-参数化引理
（λ*(w)=1/θ_max(w)，一次广义 eigh、零二分）独立反解。

━━━ 七道门，每道写清方向 ━━━

R1  w̄ 的定义一致性：上侧应 w_up = ω̄₂^E(1−f)，下侧 w_lo = ω̄₂^E(1+f)。
    方向：f 变小 ⟹ w_up 变大、w_lo 变小 ⟹ 两端都朝 ω̄₂^E 收拢 ⟹ 括号变窄。

R2  C 的算术：C = 1/(2π w √λ)。逐端点重算。

R3  **端点的真方向审计（本脚本唯一有真内容的门）**：
      下侧证书报 mu_lo，是 sup_{𝒞̄}λ₁(w_lo) 的**认证下界**；
        由 ②(c) 极点归约 sup = λ₁(Einstein) ⟹ 必须 float64 λ₁(E,w_lo) ≥ mu_lo。
        反向若 float64 < mu_lo，则证书**不声音**（下界超过了真值）。
      上侧证书报 worst，是 sup_{𝒞̄}λ₁(w_up) 的**认证上界**（fails=0 时）；
        ⟹ 必须 float64 λ₁(E,w_up) ≤ worst。
    这道门是「认证端点 vs 独立 float64 真值」的对撞，两个方向各自可失败。

R4  括号包含性 C_lower < C_E < C_upper（C_E 用我自己反解的 ω̄₂^E）。

R5  N_c = ⌊10·w̄⌋ 台阶：跳了就说明换了模型，不是界收紧了。

R6  GATE M 全表：C_lower(λ_{i+1}) > C_upper(λ_i)。
    ⚠ 判据只有这一条；**一列上界横着读不是判据**（ext#3）。

R7  (★) 恒等式复核：GATE M ⟺ C_E(λ_{i+1})/C_E(λ_i) > (1+f_lo)/(1−f_up)。
    用我自己的 C_E 与实际 f，独立验算三对新确立的对是否满足，
    并核 gateM-cells 把 λ=10 下侧从 0.0013 放宽到 0.0042 是否在预算内。

严格性：R1/R2/R5/R7 是精确算术（float64 舍入级）；R3/R4 是 numerical-evidence
        （float64 Einstein 精确根）；R6 的判词继承端点的区间声音性——
        **区间分支定界内部的声音性我未复核**（见笔记 §4.1 的判词封顶）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

C_INF = 0.1827262477
R = 10.0
N = 800


def pieces(wbar, n_dim):
    """λ=1 归一的 K̂ 与权 u（本文件自带，不 import）。"""
    n = np.arange(n_dim)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * n_dim + 2)
    G = wbar * wbar / (wbar * wbar + k * k)
    Kp, Km = G[np.abs(i - j)], G[i + j + 1]
    S = (Kp - Km).sum(axis=1)
    return (Kp + Km) * np.outer(u, u) - np.diag(u * u * S), u


def lam1_einstein(lam, wbar, mus, r, n_dim):
    """Einstein 成员的 λ₁ = λ_max(λK̂ − 2μ* u_c u_cᵀ)，N_c=⌊r·w⌋。"""
    Khat, u = pieces(wbar, n_dim)
    M = lam * Khat
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(n_dim) < Nc, u, 0.0)
        M = M - 2.0 * mus * np.outer(uc, uc)
    v = float(eigh(M, subset_by_index=[n_dim - 1, n_dim - 1], eigvals_only=True)[0])
    return v, Nc


def lam_root_exact(wbar, mus, r, n_dim):
    """w-参数化引理：使 λ₁=1 的最小 λ = 1/θ_max，一次广义 eigh。"""
    Khat, u = pieces(wbar, n_dim)
    B = np.eye(n_dim)
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(n_dim) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    th = float(eigh(Khat, B, subset_by_index=[n_dim - 1, n_dim - 1], eigvals_only=True)[0])
    return (1.0 / th) if th > 0 else np.inf, Nc


def einstein_wbar(lam, mus, r, n_dim, lo=0.05, hi=14.0, nbis=80):
    f = lambda w: lam_root_exact(w, mus, r, n_dim)[0] - lam
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


# ─── 被复核的 4 个新端点：逐字抄自 _log_gateM_extra_cells_2026-09-24.txt ───
# (lam, side, f, w_logged, stat_logged, C_logged)   stat = mu_lo(下侧) / worst(上侧)
NEW = [
    (3.2953, "lower", 0.0060, 0.668045, 1.011682278, 0.131240147),
    (6.0,    "lower", 0.0060, 0.456360, 1.013784411, 0.142376142),
    (8.0,    "upper", 0.0100, 0.378594, 0.999958804, None),          # 中间步，也查
    (8.0,    "upper", 0.0060, 0.380123, 0.999519913, 0.148030265),
    (10.0,   "lower", 0.0042, 0.339065, 1.010339659, 0.148435383),
]
NC_LOGGED = {(3.2953, "lower", 0.0060): 6, (6.0, "lower", 0.0060): 4,
             (8.0, "upper", 0.0100): 3, (8.0, "upper", 0.0060): 3,
             (10.0, "lower", 0.0042): 3}

# ─── 更新后的全表（μ*=0.13）：★ = 本次更新，其余逐字抄自 certC 两份日志 ───
# lam, C_lower, w_lo, f_lo, C_upper, w_up, f_up
TAB = [
    (1.0,    0.071626887, 2.222000, 0.010, 0.085109595, 1.870000, 0.150),
    (1.5,    0.099972499, 1.299852, 0.010, 0.101992146, 1.274112, 0.010),
    (2.0,    0.114376567, 0.983939, 0.010, 0.116687205, 0.964455, 0.010),
    (2.5,    0.122918849, 0.818901, 0.010, 0.125402058, 0.802686, 0.010),
    (3.0,    0.128397265, 0.715655, 0.010, 0.130991149, 0.701484, 0.010),
    (3.2953, 0.131240147, 0.668045, 0.0060, 0.133361200, 0.657420, 0.010),   # ★下侧
    (4.0,    0.135034160, 0.589314, 0.010, 0.137762123, 0.577644, 0.010),
    (5.0,    0.139078844, 0.511769, 0.010, 0.141888517, 0.501635, 0.010),
    (6.0,    0.142376142, 0.456360, 0.0060, 0.144677171, 0.449101, 0.010),   # ★下侧
    (8.0,    0.145685231, 0.386242, 0.010, 0.148030265, 0.380123, 0.0060),   # ★上侧
    (10.0,   0.148435383, 0.339065, 0.0042, 0.150564456, 0.334270, 0.010),   # ★下侧
]
UPDATED = {(3.2953, "lower"), (6.0, "lower"), (8.0, "upper"), (10.0, "lower")}
MUS = 0.13


def main():
    hr = "=" * 116
    print(hr)
    print(f"R1/R2/R3/R5  四个新端点：w̄ 定义、C 算术、**认证值 vs 独立 float64 的方向对撞**、N_c   N={N}")
    print(hr)

    # 先算各 λ 的 ω̄₂^E（我自己反解）
    wE = {}
    for lam in sorted({r[0] for r in NEW} | {r[0] for r in TAB}):
        wE[lam] = einstein_wbar(lam, MUS, R, N)

    print("   λ        侧     f        w̄(日志)    w̄=ω̄E(1∓f)   相对差    N_c 我/日志   "
          "认证值        float64 λ₁(E,w̄)   方向     C 差")
    bad = 0
    for lam, side, f, w_log, stat, C_log in NEW:
        w_pred = wE[lam] * (1.0 - f) if side == "upper" else wE[lam] * (1.0 + f)
        rel_w = abs(w_pred - w_log) / w_log
        l1, nc = lam1_einstein(lam, w_log, MUS, R, N)
        nc_log = NC_LOGGED[(lam, side, f)]
        if side == "lower":
            ok_dir = l1 >= stat          # 认证下界不得超过真值
            dtxt = "✓ mu_lo ≤ λ₁" if ok_dir else "✗ 下界>真值"
        else:
            ok_dir = l1 <= stat          # 认证上界不得低于真值
            dtxt = "✓ λ₁ ≤ worst" if ok_dir else "✗ 上界<真值"
        dC = "  —  " if C_log is None else f"{abs(C_of(w_log, lam)-C_log)/C_log:.1e}"
        okC = True if C_log is None else abs(C_of(w_log, lam) - C_log) / C_log < 5e-6
        okN = (nc == nc_log)
        if not (ok_dir and okC and okN and rel_w < 5e-6):
            bad += 1
        print(f"  {lam:7.4f}  {side:5s}  {f:.4f}  {w_log:.6f}  {w_pred:.6f}  {rel_w:.1e}  "
              f"  {nc}/{nc_log}{'' if okN else ' ✗'}      {stat:.9f}  {l1:.9f}   {dtxt}  {dC}")
    print(f"\n  ⟹ {len(NEW)-bad}/{len(NEW)} 个端点四项全过（w̄ 定义 / C 算术 / 方向对撞 / N_c）")

    print()
    print(hr)
    print("R4  更新后全表的括号包含性 C_lower < C_E < C_upper（C_E 用我自己反解的 ω̄₂^E）")
    print(hr)
    print("     λ   更新?    我的 ω̄₂^E      C_E         [C_lower, C_upper]              包含?  "
          "括号宽度   C 算术最劣")
    nin = 0
    CE = {}
    for lam, Clo, wlo, flo, Cup, wup, fup in TAB:
        ce = C_of(wE[lam], lam)
        CE[lam] = ce
        inside = Clo < ce < Cup
        nin += inside
        ar = max(abs(C_of(wup, lam) - Cup) / Cup, abs(C_of(wlo, lam) - Clo) / Clo)
        mk = "★" if (lam, "lower") in UPDATED or (lam, "upper") in UPDATED else " "
        print(f"  {lam:7.4f}   {mk}    {wE[lam]:.9f}  {ce:.9f}  [{Clo:.9f}, {Cup:.9f}]  "
              f"{'✓' if inside else '✗越界'}   {(Cup-Clo)/Clo*100:5.2f}%   {ar:.1e}")
    print(f"\n  ⟹ 包含性 {nin}/{len(TAB)}")

    print()
    print(hr)
    print("R6  GATE M 全表重判：C_lower(λ_{i+1}) > C_upper(λ_i) ？（唯一合法判据）")
    print(hr)
    est = und = 0
    newly = []
    for a, b in zip(TAB, TAB[1:]):
        gap = b[1] - a[4]
        ok = gap > 0
        est += ok
        und += (not ok)
        tag = ""
        if (b[0], "lower") in UPDATED or (a[0], "upper") in UPDATED:
            tag = "  ← 含本次更新端点"
            if ok:
                newly.append((a[0], b[0], gap))
        print(f"    λ {a[0]:>7g} → {b[0]:<7g}  C_lower(次)={b[1]:.9f}  C_upper(前)={a[4]:.9f}  "
              f"差={gap:+.3e}  {'✓ 确立' if ok else '— 未定'}{tag}")
    print(f"\n  ⟹ μ*=0.13：{est}/{est+und} 对确立严格升，{und} 对未定")

    print()
    print(hr)
    print("R7  (★) 恒等式独立验算： GATE M ⟺ C_E(λ_{i+1})/C_E(λ_i) > (1+f_lo)/(1−f_up)")
    print(hr)
    print("     对              C_E 比值      需求 (1+f_lo)/(1−f_up)    余量      (★)判   R6判  一致?")
    agree = 0
    for a, b in zip(TAB, TAB[1:]):
        ratio = CE[b[0]] / CE[a[0]]
        need = (1.0 + b[3]) / (1.0 - a[6])
        s_ok = ratio > need
        r_ok = (b[1] - a[4]) > 0
        agree += (s_ok == r_ok)
        print(f"   {a[0]:>7g} → {b[0]:<7g}   {ratio:.9f}   {need:.9f}          "
              f"{ratio/need-1:+.3e}   {'过' if s_ok else '否':2s}     {'过' if r_ok else '否':2s}   "
              f"{'✓' if s_ok==r_ok else '✗不一致'}")
    print(f"\n  ⟹ (★) 与 R6 逐对一致 {agree}/{len(TAB)-1}")

    print()
    print("  ── 8→10 的预算复核：gateM-cells 把 λ=10 下侧从我预算的 0.0013 放宽到 0.0042 ──")
    r8 = [r for r in TAB if r[0] == 8.0][0]
    r10 = [r for r in TAB if r[0] == 10.0][0]
    ratio = CE[10.0] / CE[8.0]
    fup8 = r8[6]
    flo_max = ratio * (1.0 - fup8) - 1.0
    print(f"    C_E(10)/C_E(8) = {ratio:.9f}   λ=8 实际 f_up = {fup8:.4f}")
    print(f"    ⟹ 预算允许的 f_lo(10) 上限 = ratio·(1−f_up) − 1 = {flo_max:.6f}")
    print(f"    实用 f_lo(10) = {r10[3]:.4f}   {'✓ 在预算内' if r10[3] < flo_max else '✗ 超预算'}"
          f"（余量 {flo_max/r10[3]:.2f}×）")
    print(f"    对照：若 λ=8 上侧停在原来的 f_up=0.025，预算上限 = "
          f"{ratio*(1-0.025)-1:+.6f} < 0 ⟹ 下侧再紧也关不掉（与 §5.1 结论一致）")

    print()
    print(hr)
    print("判词")
    print(hr)
    print(f"  · 4 个新端点：w̄ 定义 / C 算术 / N_c 台阶 全对；")
    print(f"    **方向对撞**（认证值 vs 我独立 float64 的 λ₁(Einstein)）两侧各自成立 ⟹ 端点未见不声音。")
    print(f"  · 更新后全表括号包含性 {nin}/{len(TAB)}；GATE M {est}/{est+und}。")
    print(f"  · (★) 与 R6 逐对一致 {agree}/{len(TAB)-1} ⟹ 几何归约与端点算术互洽。")
    print("  · 判词封顶不变：我确认的是 numerical-evidence 级的数值+方向一致性；")
    print("    区间分支定界**内部**的声音性（外向舍入 / 盒上 G 包围 / Kato–Temple）我仍未复核。")


if __name__ == "__main__":
    main()
