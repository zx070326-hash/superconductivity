# -*- coding: utf-8 -*-
"""主会话从零复核 + 新结果：从 Errea LaH10 的 Tc^Mc / Tc^AD 两列**反解** ω̄₂/ω_log。

起因：agent `lah10-omega2` 交付阴性结论，其中一条断言为
    「也没有同时报 Tc^AD 与 Tc^McMillan 供反解」
而我一手读过的 Extended Data Table I 表头是
    λ | ω_log(meV) | Tc^Mc_{μ*=0.1}(K) | Tc^AD_{μ*=0.1}(K) | Tc^ME_ani(K) | Tc^SCDFT(K)
两列都在、且 **μ* 相同**。⟹ 该断言为假，反解路线是开的。

因为 Tc^AD = f1 f2 · Tc^Mc（同一 ω_log、同一 μ*），而 f1 只依赖 (λ, μ*)：
    f1 f2 = Tc^AD / Tc^Mc   ⟹   f2 = (Tc^AD/Tc^Mc)/f1   ⟹   反解 r := ω̄₂/ω_log。

本脚本四节，每节都是**零参数**检验（不拟合任何东西）：
  (A) 用 AD75 Table I 的 29 行钉死 f2 的形式（仓库里有三个互相矛盾的实现）
  (B) 用 (λ, ω_log, μ*=0.1) 零参数重算 Errea 的 Tc^Mc 列 —— 钉死它的前因子约定
  (C) 反解 r = ω̄₂/ω_log（含双根的物理甄别）
  (D) 代回 Tier-0，判 LaH10 是否夺走 Nb3Ge(0.6260) 的最紧点

方向审计写在各节。运行：
  OMP_NUM_THREADS=1 python theory/notes/_scratch_ir/_chk_main_errea_f2_backsolve_omega2_2026-09-25.py
"""
import math
import sys
import io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

OK = True


def chk(label, cond, extra=""):
    global OK
    if not cond:
        OK = False
    print(f"  [{'PASS' if cond else 'FAIL'}] {label}" + (f"   {extra}" if extra else ""))


# ---------------------------------------------------------------- 常数
C_INF = 0.1827262477          # = sqrt(g(2))/(2pi)，analytic
MEV_TO_K = 11.604518          # meV -> K
NB3GE_R = 0.6260              # 现头条最紧点

# AD75 Eq.(35) 的辅助常数（仓库权威：theory/derivations/p1_ad_sensitivity.py:19-20）
def Lam1(mus):
    return 2.46 * (1 + 3.8 * mus)


def Lam2(mus, r):
    return 1.82 * (1 + 6.3 * mus) * r


def f1_of(lam, mus):
    return (1 + (lam / Lam1(mus)) ** 1.5) ** (1 / 3)


def f2_of(lam, mus, r, form):
    """form: 'r-1' (AD75 原式) / 'r2-1' (p1_ad_sensitivity.py 的写法) / 'sqrt' (已知 bug)"""
    L2 = Lam2(mus, r)
    if form == "r-1":
        return 1 + (r - 1) * lam**2 / (lam**2 + L2**2)
    if form == "r2-1":
        return 1 + (r**2 - 1) * lam**2 / (lam**2 + L2**2)
    if form == "sqrt":
        return 1 + math.sqrt(abs((r - 1) * lam**2 / (lam**2 + L2**2)))
    raise ValueError(form)


def tc_mcmillan(lam, mus, wlog_K):
    """AD75 改写形式的 McMillan 式：Tc = (ω_log/1.2) exp[-1.04(1+λ)/(λ-μ*(1+0.62λ))]"""
    return (wlog_K / 1.20) * math.exp(-1.04 * (1 + lam) / (lam - mus * (1 + 0.62 * lam)))


# ================================================================ (A)
print("=" * 78)
print("(A) 用 AD75 Table I 的 29 行钉死 f2 的形式（零参数：λ,ω_log,ω̄₂,μ* 全是表里的）")
print("=" * 78)
print("  仓库里三个互相矛盾的实现：")
print("    memory `ad75-f2-has-no-square-root` + _chk_main_ad_f2_sqrt_bug : (r-1)")
print("    theory/derivations/p1_ad_sensitivity.py:23                     : (r**2-1)")
print("    theory/derivations/P4_strong_coupling_expansion.md:123         : 带开方")
print("  ⟹ 不靠「哪个实现多」定案，靠**它们对 AD 自己那张表的复现误差**定案。\n")

# AD1975 Table I，29 行，页图逐行转录
# (material, w_log(K), wbar1, wbar2(K), lambda, w_ph, mu*, Tc(K), eta)
T1 = [
    ("Pb", 56, 60, 65, 1.55, 110, 0.105, 7.2, 2.4),
    ("In", 68, 79, 89, 0.805, 179, 0.097, 3.40, 1.3),
    ("Sn", 99, 110, 121, 0.72, 209, 0.092, 3.75, 2.2),
    ("Hg", 29, 38, 49, 1.6, 162, 0.098, 4.19, 1.4),
    ("Tl", 52, 58, 64, 0.795, 127, 0.111, 2.36, 1.2),
    ("Ta", 132, 140, 148, 0.69, 228, 0.093, 4.48, 4.9),
    ("a-Ga", 55, 77, 101, 1.62, 291, 0.095, 8.56, 2.1),
    ("b-Ga", 87, 108, 129, 0.97, 285, 0.092, 5.90, 2.0),
    ("Tl0.9Bi0.1", 48, 55, 62, 0.78, 120, 0.099, 2.30, 1.1),
    ("Pb0.4Tl0.6", 48, 56, 62, 1.15, 121, 0.094, 4.60, 1.6),
    ("Pb0.6Tl0.4", 50, 57, 62, 1.38, 119, 0.103, 5.90, 2.0),
    ("Pb0.8Tl0.2", 50, 56, 61, 1.53, 116, 0.101, 6.80, 2.1),
    ("Pb0.6Tl0.2Bi0.2", 48, 53, 58, 1.81, 112, 0.111, 7.26, 2.3),
    ("Pb0.9Bi0.1", 50, 56, 60, 1.66, 108, 0.081, 7.65, 2.2),
    ("Pb0.8Bi0.2", 46, 52, 57, 1.88, 109, 0.093, 7.95, 2.3),
    ("Pb0.7Bi0.3", 47, 52, 57, 2.01, 110, 0.092, 8.45, 2.4),
    ("Pb0.65Bi0.35", 45, 50, 55, 2.13, 110, 0.093, 8.95, 2.4),
    ("In0.9Tl0.1", 63, 75, 86, 0.85, 176, 0.103, 3.28, 1.4),
    ("In0.73Tl0.27", 55, 67, 77, 0.93, 166, 0.110, 3.36, 1.4),
    ("In0.67Tl0.33", 57, 68, 79, 0.90, 167, 0.110, 3.26, 1.4),
    ("In0.57Tl0.43", 53, 64, 74, 0.85, 165, 0.117, 2.60, 1.3),
    ("In0.5Tl0.5", 53, 64, 73, 0.83, 163, 0.110, 2.52, 1.3),
    ("In0.27Tl0.73", 42, 53, 63, 1.09, 151, 0.094, 3.64, 1.4),
    ("In0.17Tl0.83", 45, 55, 63, 0.98, 144, 0.101, 3.19, 1.3),
    ("In0.07Tl0.93", 49, 56, 63, 0.89, 131, 0.107, 2.77, 1.3),
    ("In2Bi", 46, 57, 67, 1.40, 174, 0.096, 5.6, 1.6),
    ("Sb2Tl7", 37, 48, 58, 1.43, 134, 0.102, 5.2, 1.6),
    ("Bi2Tl", 47, 53, 59, 1.63, 120, 0.101, 6.4, 2.1),
    ("a-Pb0.45Bi0.55", 29, 38, 47, 2.59, 128, 0.116, 7.0, 2.1),
]
chk("AD75 Table I 行数 = 29", len(T1) == 29, f"实得 {len(T1)}")

res = {}
for form in ("r-1", "r2-1", "sqrt"):
    errs = []
    for (nm, wlog, _w1, wb2, lam, _wph, mus, Tc, _eta) in T1:
        r = wb2 / wlog
        tc = f1_of(lam, mus) * f2_of(lam, mus, r, form) * tc_mcmillan(lam, mus, wlog)
        errs.append(abs(tc / Tc - 1))
    res[form] = sum(errs) / len(errs)
    print(f"  f2 用 {form:>6} : |Tc_AD/Tc_exp-1| 均值 = {res[form]*100:6.2f}%   最大 {max(errs)*100:6.2f}%")

print()
print("  memory 记载的独立指纹：AD75 原式 3.85%、带开方 bug 13.15%")
chk("(r-1) 复现 AD75 自己的表到 3.85%（与 memory 指纹一致，容差 0.05%）",
    abs(res["r-1"] * 100 - 3.85) < 0.05, f"实得 {res['r-1']*100:.2f}%")
chk("带开方变体复现 13.15%（负对照：确认我调的是同一套数据/同一条路径）",
    abs(res["sqrt"] * 100 - 13.15) < 0.05, f"实得 {res['sqrt']*100:.2f}%")
chk("(r-1) 严格优于 (r**2-1) ⟹ p1_ad_sensitivity.py:23 的 r**2-1 是该文件的局部 bug",
    res["r-1"] < res["r2-1"], f"{res['r-1']*100:.2f}% vs {res['r2-1']*100:.2f}%")
print("  ⟹ 承重形式钉死：f2 = 1 + (r-1)λ²/(λ²+Λ₂²)，Λ₂ = 1.82(1+6.3μ*)r。")
print("     [方向] 本节是双边检验（复现误差），不是单边界，无方向风险。")

# ================================================================ (B)
print()
print("=" * 78)
print("(B) 零参数重算 Errea Extended Data Table I 的 Tc^Mc 列 —— 钉死前因子约定")
print("=" * 78)
print("  待检前提：Errea 的「McMillan equation」用的是 AD75 改写形式 (ω_log/1.2)·exp[...]，")
print("  而不是 McMillan 原文的 (Θ_D/1.45)·exp[...]。表里**只印 ω_log、不印 Θ_D**，")
print("  且正文 line 505-507 明说两列同为 μ*=0.1 ⟹ 合理猜测。但猜测必须被检验：")
print("  Tc^Mc 由 (λ, ω_log, μ*=0.1) **唯一确定、零自由度** ⟹ 能直接对表。\n")

MUS = 0.10
# Errea 1907.11916 Extended Data Table I，7 行，一手页图 errea_lah10_tab-11.png 逐字转录
# (system, P(GPa), lambda, w_log(meV), Tc^Mc(K), Tc^AD(K), Tc^ME_ani(K), Tc^SCDFT(K))
ERREA = [
    ("LaH10", 129, 3.62, 76.4, 171.8, 252.6, 255.3, 230),
    ("LaH10", 163, 2.67, 96.4, 197.1, 247.0, 242.8, 225),
    ("LaH10", 214, 2.06, 115.5, 196.3, 235.9, 237.9, 210),
    ("LaH10", 264, 1.73, 126.6, 189.5, 219.2, 216.9, 201),
    ("LaD10", 159, 3.14, 63.5, 135.0, 184.2, 180.4, 171),
    ("LaD10", 210, 2.21, 81.7, 145.5, 176.5, 172.9, 158),
    ("LaD10", 260, 1.80, 92.2, 142.2, 164.6, 157.9, 151),
]
print(f"  {'系统':<6}{'P':>5}{'λ':>7}{'ω_log':>8}{'Tc^Mc印':>9}{'Tc^Mc算':>9}{'相对差':>9}")
mc_rel = []
for (nm, P, lam, wlog_meV, tcmc, tcad, _a, _s) in ERREA:
    pred = tc_mcmillan(lam, MUS, wlog_meV * MEV_TO_K)
    d = pred / tcmc - 1
    mc_rel.append(abs(d))
    print(f"  {nm:<6}{P:5d}{lam:7.2f}{wlog_meV:8.1f}{tcmc:9.1f}{pred:9.2f}{d*100:8.2f}%")
print()
worst = max(mc_rel)
nout = sum(1 for d in mc_rel if d < 0.005)
print(f"  相对差 < 0.5% 的行数 = {nout}/7，最大相对差 {worst*100:.2f}%")
chk("6/7 行零参数重算 Tc^Mc 落在 0.5% 内 ⟹ 前因子约定 = AD75 改写式，**已证不是假设**",
    nout >= 6, f"{nout}/7")
print("     [前提] 这一节把『Tc^Mc 与 Tc^AD 共用同一个 ω_log 和同一个 μ*』从猜测升为实测：")
print("     若 Errea 用 Θ_D/1.45 或不同 μ*，不可能有 6 行同时落在 0.25% 内（λ 跨 1.73–3.62）。")
print()
print("  ⚠ 异常行：LaH10@163 差 -3.30%，而其余 6 行全 ≤ 0.25%。**不许当噪声放过**——")
print("     它恰好是 (D) 里唯一给出「夺冠」的那一行（见 (E) 的专项追查）。")
chk("异常是**孤立**的（恰 1 行 > 0.5%）⟹ 不是系统性约定错，而是单行数据问题",
    sum(1 for d in mc_rel if d >= 0.005) == 1,
    f"超差行数 {sum(1 for d in mc_rel if d >= 0.005)}")
MC_PRED = {(nm, P): tc_mcmillan(lam, MUS, wl * MEV_TO_K)
           for (nm, P, lam, wl, _a, _b, _c, _d) in ERREA}

# ================================================================ (C)
print()
print("=" * 78
      )
print("(C) 反解 r = ω̄₂/ω_log")
print("=" * 78)
print("  f1 f2 = Tc^AD/Tc^Mc，f1 = [1+(λ/Λ₁)^{3/2}]^{1/3} 只含 (λ,μ*) ⟹ f2 已知。")
print("  解 f2 = 1 + (r-1)λ²/(λ²+c r²)，c := [1.82(1+6.3μ*)]² ⟹ 二次方程")
print("      (f2-1)c·r² - λ²·r + [λ²+(f2-1)λ²] = 0 ... 整理见代码")
print("  ⚠ f2(r) **非单调**：在 r* = 1+√(1+λ²/c) 处取极大 ⟹ 一般有两个根，必须甄别。\n")

cc = (1.82 * (1 + 6.3 * MUS)) ** 2
print(f"  c = [1.82(1+6.3·0.10)]² = {cc:.6f}")
print()
hdr = f"  {'系统':<6}{'P':>5}{'λ':>6}{'AD/Mc':>8}{'f1':>7}{'f2':>8}{'r*峰':>7}{'小根r':>8}{'大根r':>9}"
print(hdr)
BACK = []
for (nm, P, lam, wlog_meV, tcmc, tcad, _a, _s) in ERREA:
    ratio = tcad / tcmc
    f1 = f1_of(lam, MUS)
    f2 = ratio / f1
    # (f2-1)(λ² + c r²) = (r-1)λ²  ⟹  (f2-1)c r² - λ² r + λ²[1 + (f2-1)] = 0
    A = (f2 - 1) * cc
    B = -lam**2
    Cq = lam**2 * f2
    disc = B * B - 4 * A * Cq
    rstar = 1 + math.sqrt(1 + lam**2 / cc)
    if disc < 0:
        print(f"  {nm:<6}{P:5d}{lam:6.2f}{ratio:8.4f}{f1:7.4f}{f2:8.4f}{rstar:7.3f}   无实根(disc<0)")
        BACK.append((nm, P, lam, wlog_meV, None, None))
        continue
    sq = math.sqrt(disc)
    r_lo = (-B - sq) / (2 * A)
    r_hi = (-B + sq) / (2 * A)
    print(f"  {nm:<6}{P:5d}{lam:6.2f}{ratio:8.4f}{f1:7.4f}{f2:8.4f}{rstar:7.3f}{r_lo:8.4f}{r_hi:9.2f}")
    BACK.append((nm, P, lam, wlog_meV, r_lo, r_hi))

print()
print("  甄别双根：")
print("   · 大根甄别不用拍阈值，用**物理不可能**：ω̄₂ ≤ ω_max，而任何声子谱的上端都远低于")
print("     500 meV（自由 H₂ 伸缩 ~4400 cm⁻¹ = 546 meV 已是刚性上限，晶体里只会更低）。")
print("     逐行检查大根对应的 ω̄₂ = r_hi·ω_log：")
for (nm, P, lam, wl, r_lo, r_hi) in BACK:
    if r_hi is not None:
        print(f"       {nm}@{P}: r_hi={r_hi:6.2f} ⟹ ω̄₂ = {r_hi*wl:7.1f} meV  "
              f"{'(> 500 ⟹ 非物理，弃)' if r_hi*wl > 500 else '(<= 500 ⟹ 无法排除!)'}")
chk("7/7 大根对应 ω̄₂ > 500 meV ⟹ 非物理，双根甄别是**物理判定**而非拍阈值",
    all(b[5] * b[3] > 500 for b in BACK if b[5] is not None),
    f"min = {min(b[5]*b[3] for b in BACK if b[5] is not None):.1f} meV")
print("   · 小根须满足 Jensen 下界 r ≥ 1（[analytic]，ln 凹）—— 这是**独立于反解**的约束，")
print("     所以它是一次真检验，不是自证。")
lo_all = [b[4] for b in BACK if b[4] is not None]
chk("7/7 小根满足 Jensen r ≥ 1（独立解析约束，反解没有被喂进这个信息）",
    all(x >= 1.0 for x in lo_all), f"min = {min(lo_all):.4f}")
chk("7/7 小根 ≤ 3.3（氢化物 ω_max/ω_log 的量级上限）",
    all(x <= 3.3 for x in lo_all), f"max = {max(lo_all):.4f}")

print()
print("  与独立标定对照（H₃S/D₃S 四点真 ω̄₂/ω_log，来自 theory/data/hydrides_theory.csv）：")
CAL = [("H3S harm", 107.165, 129.32), ("H3S anh ", 109.995, 149.26),
       ("D3S harm", 79.562, 92.97), ("D3S anh ", 85.535, 105.43)]
cal_r = [w2 / wl for (_n, wl, w2) in CAL]
for (n, wl, w2), rr in zip(CAL, cal_r):
    print(f"    {n}:  ω_log={wl:7.3f}  ω̄₂={w2:7.2f}  r={rr:.4f}")
print(f"    标定范围 [{min(cal_r):.4f}, {max(cal_r):.4f}]")
print(f"    LaH10/LaD10 反解范围 [{min(lo_all):.4f}, {max(lo_all):.4f}]")
print("    ⟹ 两族**量级一致**（都在 1.1–1.4）。这是一次**外部一致性**检验：")
print("       反解完全不知道 H₃S 的任何数，却落进同一个窄带。")

# ================================================================ (D)
print()
print("=" * 78)
print("(D) 代回 Tier-0：LaH10 是否夺走 Nb3Ge(0.6260) 的最紧点？")
print("=" * 78)
print("  Tier-0: kTc ≤ C_∞·√(λ⟨ω²⟩) ⟹ Tc_ceiling[K] = C_∞·√λ·ω̄₂[meV]·11.604518")
print("  r_up 用 ω_log 顶替 ω̄₂（Jensen 单边、朝 PASS）；现在有了 ω̄₂ ⟹ r_true = r_up / r。\n")

# Tc_exp 取值（与 hydride_omegalog_corpus_extension_2026-09-25.md 同一约定，保守取值）
TCEXP = {("LaH10", 129): 250.0, ("LaH10", 163): 260.0, ("LaH10", 214): 260.0,
         ("LaH10", 264): 260.0, ("LaD10", 159): 180.0, ("LaD10", 210): 180.0,
         ("LaD10", 260): 180.0}

print(f"  {'系统':<6}{'P':>5}{'Tc_exp':>8}{'ceil_up':>9}{'r_up':>8}{'r反解':>8}{'ceil_true':>10}{'r_true':>8}  判决")
rows = []
for (nm, P, lam, wlog_meV, r_lo, _r_hi) in BACK:
    tce = TCEXP[(nm, P)]
    ceil_up = C_INF * math.sqrt(lam) * wlog_meV * MEV_TO_K
    r_up = tce / ceil_up
    ceil_true = ceil_up * r_lo
    r_true = tce / ceil_true
    verdict = "夺冠(比 Nb3Ge 紧)" if r_true > NB3GE_R else "未夺冠"
    rows.append((nm, P, r_up, r_lo, r_true, verdict))
    print(f"  {nm:<6}{P:5d}{tce:8.1f}{ceil_up:9.2f}{r_up:8.4f}{r_lo:8.4f}{ceil_true:10.2f}{r_true:8.4f}  {verdict}")

print()
r_ups = [x[2] for x in rows]
r_trues = [x[4] for x in rows]
chk("(a) 天花板成立：7/7 行 r_true < 1（本节结论不依赖 Jensen，用的是反解出的真 ω̄₂）",
    all(x < 1.0 for x in r_trues), f"max r_true = {max(r_trues):.4f}")
chk("方向自查：r_true ≤ r_up 逐行成立（因 r ≥ 1）⟹ 与 Jensen 单边性一致，无矛盾",
    all(t <= u + 1e-12 for t, u in zip(r_trues, r_ups)))

mx = max(rows, key=lambda x: x[4])
print()
print(f"  【读法一：照印出的 Tc^Mc 用】最紧行 {mx[0]}@{mx[1]} GPa，r_true = {mx[4]:.4f}"
      f"   （Nb3Ge = {NB3GE_R}）⟹ {'易主' if mx[4] > NB3GE_R else '不易主'}")
rows_ok = [x for x in rows if not (x[0] == "LaH10" and x[1] == 163)]
mx_ok = max(rows_ok, key=lambda x: x[4])
print(f"  【剔除 (B) 的异常行后】最紧行 {mx_ok[0]}@{mx_ok[1]} GPa，"
      f"r_true = {mx_ok[4]:.4f} ⟹ {'易主' if mx_ok[4] > NB3GE_R else '不易主'}")
print("  ⟹ 判决**完全由那一行不可复现的 Tc^Mc 决定**。转 (E) 追查，不在此处下结论。")

# ================================================================ (E)
print()
print("=" * 78)
print("(E) 专项追查 LaH10@163：印出的 Tc^Mc=197.1 与公式值 190.60 差 -3.30%")
print("=" * 78)
print("  [方向审计——先写清，再算] 印出的 Tc^Mc 若**偏高**：")
print("    ratio = Tc^AD/Tc^Mc 偏低 ⟹ f2 偏低 ⟹ 反解 r 偏低 ⟹ ceiling_true 偏低")
print("    ⟹ r_true 偏**高** ⟹ **偏向「夺冠」**。")
print("  ⟹ 唯一支持易主的那一行，其未解释的偏差方向**正好朝易主** ⟹ 不得照用。\n")

print("  假设 H：印出的 Tc^Mc 是笔误，Tc^AD 列本身是对的（由公式值算出）。")
print("  H 的**可检验推论**：用公式值 Tc^Mc 反解出的 r，应当与相邻压强行连成一条光滑曲线。")
print(f"  {'系统':<6}{'P':>5}{'r(用印出Mc)':>13}{'r(用公式Mc)':>13}")
r_alt = {}
for (nm, P, lam, wlog_meV, r_lo, _r_hi) in BACK:
    tcad = [e[5] for e in ERREA if e[0] == nm and e[1] == P][0]
    ratio2 = tcad / MC_PRED[(nm, P)]
    f1 = f1_of(lam, MUS)
    f2 = ratio2 / f1
    A = (f2 - 1) * cc
    B = -lam**2
    Cq = lam**2 * f2
    disc = B * B - 4 * A * Cq
    rr = (-B - math.sqrt(disc)) / (2 * A) if disc >= 0 else float("nan")
    r_alt[(nm, P)] = rr
    print(f"  {nm:<6}{P:5d}{r_lo:13.4f}{rr:13.4f}")

print()
la_alt = [r_alt[("LaH10", p)] for p in (129, 163, 214, 264)]
print(f"  LaH10 四行（用公式 Mc）r = {['%.4f' % x for x in la_alt]}")
print(f"  其中 @163 的邻居是 @129={la_alt[0]:.4f} 与 @214={la_alt[2]:.4f}")
inside = min(la_alt[0], la_alt[2]) - 0.02 <= la_alt[1] <= max(la_alt[0], la_alt[2]) + 0.02
chk("H 的推论成立：用公式 Mc 后 @163 的 r 落进两个邻居之间（±0.02）⟹ 支持「印出值是笔误」",
    inside, f"@163 = {la_alt[1]:.4f}，邻居区间 [{min(la_alt[0],la_alt[2]):.4f}, {max(la_alt[0],la_alt[2]):.4f}]")
print(f"  对照：照印出值反解时 @163 给 r = {BACK[1][4]:.4f}，**低于两个邻居**（"
      f"{BACK[0][4]:.4f}, {BACK[2][4]:.4f}）⟹ 单点离群，与 (B) 的 -3.30% 同一根因。")
print("  ⚠ 这是**一致性证据，不是证明**：我没有 Errea 的勘误表，也没有他们的 α²F。")
print("     两种读法都保留，判决按**更不利于我方主张**的那个报（单边界方向规则）。")

print()
print("  两种读法下的全语料最紧点：")
rows2 = []
for (nm, P, lam, wlog_meV, r_lo, _r_hi) in BACK:
    tce = TCEXP[(nm, P)]
    ceil_up = C_INF * math.sqrt(lam) * wlog_meV * MEV_TO_K
    r_up = tce / ceil_up
    rows2.append((nm, P, r_up / r_lo, r_up / r_alt[(nm, P)]))
m1 = max(x[2] for x in rows2)
m2 = max(x[3] for x in rows2)
a1 = max(rows2, key=lambda x: x[2])
a2 = max(rows2, key=lambda x: x[3])
print(f"    读法一（印出 Mc）：max r_true = {m1:.4f} ({a1[0]}@{a1[1]}) ⟹ "
      f"{'易主' if m1 > NB3GE_R else '不易主'}")
print(f"    读法二（公式 Mc）：max r_true = {m2:.4f} ({a2[0]}@{a2[1]}) ⟹ "
      f"{'易主' if m2 > NB3GE_R else '不易主'}")
print(f"    区间 [{min(m1,m2):.4f}, {max(m1,m2):.4f}] vs 阈值 {NB3GE_R}")
straddle = min(m1, m2) < NB3GE_R < max(m1, m2)
print(f"    ⟹ {'**跨过**阈值 ⟹ 仍判不死' if straddle else '同侧 ⟹ 可判'}")
chk("两读法区间跨过 0.626 ⟹ 诚实判词是「未定」，**不许挑有利的一端报**",
    straddle, f"[{min(m1,m2):.4f}, {max(m1,m2):.4f}]")
print()
print("  ★ 交付量：反解出的 ω̄₂（meV）—— 这是 hydrides_theory.csv 里原本 not-reported 的格子")
print(f"  {'系统':<6}{'P':>5}{'ω_log':>8}{'ω̄₂(读法一)':>13}{'ω̄₂(读法二)':>13}")
for (nm, P, lam, wl, r_lo, _r_hi) in BACK:
    print(f"  {nm:<6}{P:5d}{wl:8.1f}{r_lo*wl:13.2f}{r_alt[(nm,P)]*wl:13.2f}")
print("  （两读法只在 LaH10@163 上有实质差别：108.95 vs 120.88 meV）")

print()
print("  **但「未定」的性质彻底变了**（这是本轮的净收益）：")
print("   · 旧「未定」：缺 ω̄₂ 这个量本身，需要采购真 α²F —— 无界的数据需求。")
print("   · 新「未定」：ω̄₂ 已从论文自己的两列反解出来，7 行全落 1.13–1.32、与 H₃S/D₃S")
print("     四点标定 [1.1685,1.3570] 同带；剩下的唯一不确定性是**一行印刷值的 3.3% 异常**。")
print("   · 且两读法**都**给 max r_true < 0.69 ⟹ (a) 天花板成立这一条**不受影响**。")

print()
print("  ⚠ 严格性定级（承重）：")
print("   · (A)(B) 是零参数复现，analytic 公式 + float64 点算 ⟹ numerical-evidence")
print("   · (C) 反解依赖『Errea 的 Tc^AD 列确实由 AD75 Eq.(35) 原式算出』。(B) 已证前因子约定，")
print("     但 f2 的**具体形式**是我从 AD75 补的；若 Errea 用了某个变体，反解出的 r 会偏。")
print("     ⟹ (C)(D) 只能是 **numerical-evidence**，不得升级。")
print("   · 方向：r 若被**低估**，ceiling_true 偏小、r_true 偏大 ⟹ 偏向「夺冠」。")
print("     ⟹ 「夺冠」这个判词需要 r 的**上估**才安全，而我手上正是可能被低估的那一侧")
print("     （(B) 的异常行就是这个形状）⟹ **不得单独报读法一的 0.6887**。")
print("     反之「(a) 天花板成立 / r_true < 1」只需 r 的下估即可 ⟹ 单边安全、可报。")

print()
print("=" * 78)
print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
print("=" * 78)
sys.exit(0 if OK else 1)
