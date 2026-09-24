# -*- coding: utf-8 -*-
"""
派活前的预算与可达性测算（2026-09-24，主会话）。

用户要求：看四条候选能不能深推，不能就收尾；GATE M 再尝试推进。
`CLAUDE.md` §「预算优先（派活前必做——已三次踩中）」：
  (a) 写下消费点——该量在哪一格、哪个不等式里被用掉；
  (b) 朝判据变难的方向退化，看下游判决在哪里翻转（预算 ∞ ⟹ 不派）；
  (c) **先在 float64 上实测所选估计子的固有松**，确认目标低于它。
      「估计子自带的松不是执行方能靠努力压掉的；目标定在固有松之下 = 派了件不可能完成的活。」

━━━ 消费点 ━━━

GATE M：`C(λ)` 在全表上严格单调升，升到 `certificate` 级。
判据（两端都是证书）   C_lower(λ_{i+1}) > C_upper(λ_i)
标度 C = 1/(2π w √λ)；上侧 w_up = wE(1−f_up) ⟹ C_upper = C_E/(1−f_up)
                      下侧 w_lo = wE(1+f_lo) ⟹ C_lower = C_E/(1+f_lo)
⟹ 判据等价于纯几何条件（analytic 的代数恒等变形）：

        C_E(λ_{i+1}) / C_E(λ_i)  >  (1 + f_lo) / (1 − f_up)                 (★)

(★) 让**需求先于算力**被解出来：给定对侧现有的 f，另一侧需要收紧到多少。

━━━ 可达性（本脚本的重点）━━━

Einstein 成员是该区域的极值成员（②(c) 的极点归约）⟹ sup_{𝒞̄} λ₁ = λ₁(Einstein)。
于是两个估计子的**固有松**可以逐点直接测出来：
  · 下侧 Rayleigh：slack_R := λ₁(E,w_lo) − mu_lo(日志)   余量 head_R := λ₁(E,w_lo) − 1
  · 上侧 B&B：     slack_B := worst(日志) − λ₁(E,w_up)   余量 head_B := 1 − λ₁(E,w_up)
判据要成立必须 slack < head。**f 越小 head 越小，而 slack 由盒子分辨率定、不随 f 同步缩小**
⟹ 存在硬下限 f_min ≈ f·(slack/head)。低于它就是不可能完成的活。

旁证（日志里已可见）：上侧 `worst` 在 f 上**非单调**（λ=2.5：0.060→0.998711、0.025→0.998596、
0.010→0.999984；λ=10：0.998808 / 0.998641 / 0.999526）⟹ 它由 slack_B 主导，不是真值的函数。

严格性：全部 numerical-evidence（float64）。(★) 是 analytic。
本脚本**不产生任何证书**，只产生「该不该派、派什么参数」的判断。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

C_INF = 0.1827262477
R = 10.0
MUS = 0.13
N = 800


def pieces(wbar, N):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * N + 2)
    G = wbar * wbar / (wbar * wbar + k * k)
    Kp, Km = G[np.abs(i - j)], G[i + j + 1]
    S = (Kp - Km).sum(axis=1)
    return (Kp + Km) * np.outer(u, u) - np.diag(u * u * S), u


def lam1_einstein(lam, wbar, mus, r, N):
    """Einstein 成员的 λ₁ = λ_max(λK̂ − 2μ* u_c u_cᵀ)，N_c=⌊r·w⌋。"""
    Khat, u = pieces(wbar, N)
    M = lam * Khat
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(N) < Nc, u, 0.0)
        M = M - 2.0 * mus * np.outer(uc, uc)
    return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0]), Nc


def lam_star(wbar, mus, r, N):
    Khat, u = pieces(wbar, N)
    B = np.eye(N)
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(N) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    th = float(eigh(Khat, B, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])
    return (1.0 / th) if th > 0 else np.inf


def einstein_wbar(lam, mus, r, N, lo=0.05, hi=14.0, nbis=80):
    f = lambda w: lam_star(w, mus, r, N) - lam
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


# ── 锚点：**逐字**抄自 certC 两份日志（w̄ 用日志印出的值，即实际被认证的那个点）──
ANCH_UP = [   # (lam, f, w_up_logged, worst_logged)
    (2.5, 0.150, 0.689174, 0.999039460), (2.5, 0.060, 0.762146, 0.998711047),
    (2.5, 0.025, 0.790524, 0.998595548), (2.5, 0.010, 0.802686, 0.999984463),
    (3.0, 0.150, 0.602284, 0.997145311), (3.0, 0.060, 0.666055, 0.999169315),
    (3.0, 0.025, 0.690855, 0.999439861), (3.0, 0.010, 0.701484, 0.999283061),
    (3.2953, 0.010, 0.657420, 0.999160096),
    (4.0, 0.010, 0.577644, 0.999945446),
    (5.0, 0.010, 0.501635, 0.999140223),
    (6.0, 0.010, 0.449101, 0.998831508),
    (8.0, 0.150, 0.325055, 0.997405701), (8.0, 0.060, 0.359473, 0.998435309),
    (8.0, 0.025, 0.372857, 0.999633996),
    (10.0, 0.150, 0.287000, 0.997905666), (10.0, 0.060, 0.317388, 0.998808239),
    (10.0, 0.025, 0.329206, 0.998641498), (10.0, 0.010, 0.334270, 0.999525712),
]
ANCH_LO = [   # (lam, f, w_lo_logged, mu_lo_logged)
    (2.5, 0.010, 0.818901, 1.017131890),
    (3.0, 0.025, 0.726284, 1.046654408), (3.0, 0.010, 0.715655, 1.018667116),
    (3.2953, 0.025, 0.680663, 1.048693422), (3.2953, 0.010, 0.670702, 1.019472255),
    (4.0, 0.010, 0.589314, 1.020855517),
    (5.0, 0.025, 0.519370, 1.055216591), (5.0, 0.010, 0.511769, 1.022029078),
    (6.0, 0.025, 0.464979, 1.057675251), (6.0, 0.010, 0.458174, 1.022994585),
    (8.0, 0.025, 0.391978, 1.060371776), (8.0, 0.010, 0.386242, 1.024047549),
    (10.0, 0.025, 0.346088, 1.061963784), (10.0, 0.010, 0.341023, 1.024665254),
]

PAIRS = [  # λ_i, λ_{i+1}, 现状 f_up(λ_i), 现状 f_lo(λ_{i+1})
    (3.0, 3.2953, 0.010, 0.010),
    (5.0, 6.0, 0.010, 0.010),
    (8.0, 10.0, 0.025, 0.010),
]


def main():
    print("=" * 116)
    print("STEP 1  需求：(★) C_E(λ')/C_E(λ) > (1+f_lo)/(1−f_up)   —— 零算力，需求先于算力解出来")
    print("=" * 116)
    wE, CE = {}, {}
    for lam in (3.0, 3.2953, 5.0, 6.0, 8.0, 10.0):
        wE[lam] = einstein_wbar(lam, MUS, R, N)
        CE[lam] = C_of(wE[lam], lam)
        l1, Nc = lam1_einstein(lam, wE[lam], MUS, R, N)
        print(f"   λ={lam:<7g} ω̄₂^E={wE[lam]:.9f}  N_c={Nc:2d}  C_E={CE[lam]:.9f}   [自检 λ₁(E)={l1:.9f}]")
    print()
    need = {}
    for li, lj, fup, flo in PAIRS:
        ratio = CE[lj] / CE[li]
        cur = (1.0 + flo) / (1.0 - fup)
        f_sym = (ratio - 1.0) / (ratio + 1.0)
        flo_need = ratio * (1.0 - fup) - 1.0
        fup_need = 1.0 - (1.0 + flo) / ratio
        need[(li, lj)] = dict(ratio=ratio, f_sym=f_sym, flo_need=flo_need, fup_need=fup_need)
        print(f"   λ {li:g}→{lj:g}:  C_E 比={ratio:.7f}   现状 (1+f_lo)/(1−f_up)={cur:.7f}   差 {cur-ratio:+.2e}")
        print(f"        对称需求 f ≤ {f_sym:.6f} | 只收下侧需 f_lo ≤ {flo_need:+.6f}（现{flo:.3f}）"
              f" | 只收上侧需 f_up ≤ {fup_need:+.6f}（现{fup:.3f}）")

    print()
    print("=" * 116)
    print("STEP 2  固有松实测（sup_{𝒞̄}λ₁ = λ₁(Einstein)，②(c) 的极点归约）")
    print("=" * 116)
    print("\n  ── 上侧 B&B：slack_B = worst − λ₁(E,w_up)；余量 head_B = 1 − λ₁(E,w_up) ──")
    print("      λ     f       w_up       λ₁(E,w_up)     worst(日志)   slack_B      head_B      slack/head")
    up_rows = []
    for lam, f, w, logged in ANCH_UP:
        l1, _ = lam1_einstein(lam, w, MUS, R, N)
        sb, hb = logged - l1, 1.0 - l1
        up_rows.append((lam, f, sb, hb))
        flag = "  ← 濒临固有松" if sb / hb > 0.8 else ""
        print(f"   {lam:6g}  {f:.3f}  {w:.6f}  {l1:.9f}  {logged:.9f}  {sb:+.3e}  {hb:.3e}  {sb/hb:7.3f}{flag}")
    print("\n  ── 下侧 Rayleigh：slack_R = λ₁(E,w_lo) − mu_lo；余量 head_R = λ₁(E,w_lo) − 1 ──")
    print("      λ     f       w_lo       λ₁(E,w_lo)     mu_lo(日志)   slack_R      head_R      slack/head")
    lo_rows = []
    for lam, f, w, logged in ANCH_LO:
        l1, _ = lam1_einstein(lam, w, MUS, R, N)
        sr, hr = l1 - logged, l1 - 1.0
        lo_rows.append((lam, f, sr, hr))
        print(f"   {lam:6g}  {f:.3f}  {w:.6f}  {l1:.9f}  {logged:.9f}  {sr:+.3e}  {hr:.3e}  {sr/hr:7.3f}")
    print()
    print("  ⟹ 两个估计子的固有松（全锚点范围）：")
    print(f"     上侧 B&B      slack_B ∈ [{min(r[2] for r in up_rows):+.2e}, {max(r[2] for r in up_rows):+.2e}]")
    print(f"     下侧 Rayleigh slack_R ∈ [{min(r[2] for r in lo_rows):+.2e}, {max(r[2] for r in lo_rows):+.2e}]")

    print()
    print("=" * 116)
    print("STEP 3  硬下限到底存不存在：slack 是常数还是与 head 成比例？")
    print("=" * 116)
    print("""
  ⚠ 自订正（本脚本初稿的 STEP 3 判据错，此处作废并替换）
    ── 旧陈述错在哪 ──
       初稿写 f_min ≈ f·(slack/head)，**预设 slack_B 在 f 缩小时保持不变**（常数松）。
       若真如此，head ∝ f 迟早被常数 slack 吃掉 ⟹ 存在硬下限。
       上面 STEP 2 的数据**直接否证**这个预设：slack_B 从 f=0.150 的 2.8e−01 一路缩到
       f=0.010 的 1.8e−02，与 head_B 同步缩小，比值稳定在 0.95–0.999。
       ⟹ 由此推出的任何「f_min」以及「λ=8 上侧 f=0.010 会 FAIL」都是**假命题**，撤回。
    ── 新陈述凭什么对 ──
       `worst ≈ 1⁻` 是**停机规则的指纹**，不是估计子的松：分支定界一旦把所有盒子证到 <1 就停，
       于是印出的 max 必然贴着 1 从下方逼近，与 f 无关。真正的可行性指标是**叶数能不能在
       maxbox 内收敛**——这是算力问题，不是固有松问题。下面用日志叶数直接测它的增长律。
""")
    print("  ── 比例性检验：slack_B / head_B 是否与 f 无关（若无关 ⟹ 无硬下限）──")
    for lam in (2.5, 3.0, 8.0, 10.0):
        sub = sorted([(f, s / h) for (L, f, s, h) in up_rows if L == lam], reverse=True)
        s = "  ".join(f"f={f:.3f}:{r:.3f}" for f, r in sub)
        print(f"     λ={lam:<5g} {s}")
    ratios = [s / h for (_, _, s, h) in up_rows]
    print(f"     ⟹ 全部 {len(ratios)} 个锚点：min={min(ratios):.3f} max={max(ratios):.3f} "
          f"（对 f 无系统趋势 ⟹ **无硬下限**）")
    print()
    print("  ── 叶数增长律（真正的可行性指标）：日志逐字 ──")
    NLEAF = {  # lam -> [(f, nleaf, 秒)]
        2.5: [(0.150, 51, 32), (0.060, 85, 55), (0.025, 117, 73), (0.010, 148, 91)],
        3.0: [(0.150, 48, 29), (0.060, 81, 50), (0.025, 110, 68), (0.010, 150, 100)],
        3.2953: [(0.150, 47, 29), (0.060, 79, 48), (0.025, 111, 68), (0.010, 146, 91)],
        4.0: [(0.150, 42, 25), (0.060, 80, 48), (0.025, 108, 66), (0.010, 141, 87)],
        5.0: [(0.150, 47, 28), (0.060, 77, 47), (0.025, 109, 67), (0.010, 139, 85)],
        6.0: [(0.150, 47, 28), (0.060, 83, 50), (0.025, 110, 67), (0.010, 139, 85)],
        8.0: [(0.150, 47, 28), (0.060, 87, 52), (0.025, 117, 71)],
        10.0: [(0.150, 52, 30), (0.060, 98, 57), (0.025, 133, 78), (0.010, 153, 90)],
    }
    print("      λ      f=0.150    f=0.060    f=0.025    f=0.010")
    for lam in sorted(NLEAF):
        cells = {f: (n, t) for f, n, t in NLEAF[lam]}
        row = "  ".join(f"{cells[f][0]:3d}叶/{cells[f][1]:3d}s" if f in cells else "   —     "
                        for f in (0.150, 0.060, 0.025, 0.010))
        print(f"   {lam:6g}   {row}")
    have = [c for lam in NLEAF for c in NLEAF[lam] if abs(c[0] - 0.010) < 1e-9]
    print(f"\n     7 个 λ 在 f=0.010 收敛：叶数 {min(c[1] for c in have)}–{max(c[1] for c in have)}，"
          f"耗时 {min(c[2] for c in have)}–{max(c[2] for c in have)}s，**零 FAIL**。")
    print("     λ=8 在 f=0.025 用 117 叶（与同 f 的中位 110 相当，不是异常格）")
    print("     ⟹ λ=8 上侧 f=0.010 预计 ~150 叶 / ~90s，**无理由预期失败**；")
    print("        它此前未算的唯一原因是第一次运行撞到**全局**时间预算（日志汇总列 BUDGET）。")
    fmin_up = {}           # 已证无硬下限 ⟹ 上侧不设 f_min
    fmin_lo = {}
    for lam, f, sr, hr in lo_rows:
        if abs(f - 0.010) < 1e-9:
            fmin_lo[lam] = 0.010 * sr / hr   # 下侧 slack 是**真常数级**(1e−6)，此式适用
    print()
    print("  ── 下侧：slack_R ~ 1e−06 是**真常数**（不随 f 缩）⟹ 此式对下侧适用 ──")
    for lam in sorted(fmin_lo):
        print(f"     λ={lam:<7g} f_min(下侧) ≈ {fmin_lo[lam]:.2e}")

    print()
    print("=" * 116)
    print("STEP 4  派活参数：带安全裕度 δ 的目标 f，以及每个目标的消费点")
    print("=" * 116)
    DELTA = 0.001   # 要求 (1+f_lo)/(1−f_up) ≤ ratio/(1+δ)，即 (★) 留 0.1% 余量
    print(f"   安全裕度 δ={DELTA}：要求 (1+f_lo)/(1−f_up) ≤ C_E 比 /(1+δ)")
    print("   （C_E 是 float64 真值代理；证书端点由实际认证的 w̄ 定，留裕度防 6 位截断与台阶）")
    orders = []
    for li, lj, fup, flo in PAIRS:
        d = need[(li, lj)]
        tgt = d["ratio"] / (1.0 + DELTA)
        print(f"\n   ── λ {li:g} → {lj:g}   C_E 比={d['ratio']:.7f}   带裕度目标 ≤ {tgt:.7f} ──")
        fup_max = 1.0 - 1.0 / tgt          # 上侧无论下侧多紧都必须满足的硬约束
        flo_keep = tgt * (1.0 - fup) - 1.0  # ← 用**现状** f_up，不是假想的 0.010
        if flo_keep > 0:
            pick = 0.008 if flo_keep > 0.0085 else (0.006 if flo_keep > 0.0065 else round(flo_keep * 0.8, 4))
            print(f"      上侧 λ={li:g} 保持现状 f_up={fup:.3f}（已有证书，**不必重跑**）")
            print(f"      下侧 λ={lj:g} 需 f_lo ≤ {flo_keep:.6f}  ⟹ **取 f_lo = {pick:.4f}**"
                  f"（Rayleigh 0.2s，硬下限 {fmin_lo.get(lj, 0):.1e}）")
            orders.append((lj, "lower", pick, f"关 λ{li:g}→{lj:g}"))
        else:
            print(f"      上侧 λ={li:g} 现状 f_up={fup:.3f} ⟹ (1+f_lo) ≤ {tgt*(1-fup):.6f} < 1，")
            print(f"        **无论下侧收多紧都关不掉** ⟹ 上侧必须重跑到 f_up ≤ {fup_max:.6f}")
            for fu in (0.010, 0.006):
                fl = tgt * (1.0 - fu) - 1.0
                if fl > 0:
                    print(f"        · 上侧到 f_up={fu:.3f} ⟹ 下侧 λ={lj:g} 需 f_lo ≤ {fl:.6f}"
                          f"  ⟹ 取 f_lo = {round(fl*0.7, 4):.4f}")
            fl10 = tgt * (1.0 - 0.010) - 1.0
            orders.append((li, "upper", 0.010, f"关 λ{li:g}→{lj:g}（硬需求 f_up≤{fup_max:.4f}）"))
            orders.append((lj, "lower", round(fl10 * 0.7, 4), f"关 λ{li:g}→{lj:g}"))
    print("\n   ── N_c 台阶自检（收紧后 N_c=⌊10w⌋ 有没有跳档？跳了就是换了模型，不是收紧）──")
    for lam, side, f, why in orders:
        w = wE[lam] * (1.0 - f) if side == "upper" else wE[lam] * (1.0 + f)
        nc_new, nc_old = int(np.floor(R * w)), int(np.floor(R * wE[lam]))
        print(f"      λ={lam:<7g} {side:<6s} f={f:<7.4f} w={w:.6f}  N_c {nc_old} → {nc_new}"
              f"  {'✓ 不跳' if nc_new == nc_old else '✗ 跳档！'}")
    print("\n   ── 派活清单（共 %d 格）──" % len(orders))
    for lam, side, f, why in orders:
        cost = "~90s B&B" if side == "upper" else "~0.2s Rayleigh"
        print(f"      λ={lam:<7g} μ*=0.13 {side:<6s} f={f:<7.4f}  [{cost}]   消费点：{why}")

    print()
    print("=" * 116)
    print("STEP 5  其余三条候选的预算（退化到判据变难的方向，看判决在哪翻转）")
    print("=" * 116)
    print("""
   ② 认证 μ*=0.20
      消费点：头条表「优势随 μ* 增大」的第四列（比值 1.8156，现 numerical-evidence）。  【2026-09-24 订正：该列现为 1.7967，原 1.8156 的 λ_ceil 取了最大根】
      退化实验：整列删掉 ⟹ 仍余 μ*=0(1.0152)/0.10(1.1853)/0.13(1.2982) 三列，
                「随 μ* 单调升」的定性结论**不翻转**，只是覆盖的 μ* 上限从 0.20 缩到 0.13。
      ⟹ 预算**有限、不翻判决**：买的是覆盖范围（物理上限 μ*≈0.2），不是通过。
      可达性：B&B 在 μ*=0 与 0.13 都跑通，0.20 只是参数 ⟹ 无固有松障碍。
      判词：**可深推，优先级低于 ③**。

   ③ 逐分支区间认证 c(μ*,r)（Leavens 侧，27 分支）
      消费点：头条比值 λ_ceil/λ_area > 1 的**分母**。
      方向审计：要证 ratio>1，需 λ_ceil 的**下界** + λ_area 的**上界**；
                λ_area 由 c 单调决定，c 已被 w-参数化引理化成显式分支 ⟹ 每支可闭式外向舍入。
      退化实验：删掉它 ⟹ 分母只有 float64 ⟹ **整个头条比较掉回 numerical-evidence**。
      ⟹ 预算**有限且承重**——四条里唯一一条「删掉会掉级」的。
      判词：**可深推，优先级最高**。

   ④ 收紧临界 λ 的 [3.278824, 3.328253]（宽 1.5%）
      消费点：报出的天花板数值本身。
      退化实验：退化到最不利端 3.278824，与 λ_area=2.5382 比 ⟹ 1.292 > 1，**判决不翻转**；
                宽度翻倍也不翻。
      ⟹ **预算 ≈ ∞（只买余量与报数位数）⟹ 按 CLAUDE.md「预算为 ∞ ⟹ 不要派活」收尾。**
""")


if __name__ == "__main__":
    main()
