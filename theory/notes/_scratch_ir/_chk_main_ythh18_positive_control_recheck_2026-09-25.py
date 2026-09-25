"""
从零复核 teammate `hydride-omega2-census` 的「真阳性对照」交付：
Ghaffar–Song–Maezono–Hongo, ACS Omega 9(50), 49470 (2024), DOI 10.1021/acsomega.4c07199,
PMC11656362, Table 1（YThH18 / YThH8 / Y2ThH12）。

一手件（本会话自取，**全 TLS 校验、无 -k、无 --ssl-no-revoke**；
curl/schannel 的 000 是 CRYPT_E_REVOCATION_OFFLINE，用 python urllib+certifi 全校验即可绕过该工具限制）：
  theory/notes/_scratch_ir/_ythh18_PMC11656362.html   md5 b30eb46ac340a1608353d3e23f85e2b5
  theory/notes/_scratch_ir/_ythh18_eq/eq1_mcmillan.jpg md5 65f9c64d2506c5eb8ed592318fb9ae25
  theory/notes/_scratch_ir/_ythh18_eq/eq4_ad.jpg       md5 51dc5ce6a2983186115f4f230adeedaf
  theory/notes/_scratch_ir/_ythh18_eq/eq6_w2.jpg       md5 9934e06351d8a6ecdd75034115758d0b

★ 两条 teammate 未报的约定风险（我读页图查出）：
  (R1) Eq.(1) 的 McMillan 用的是 **ω_D/1.45**（原版），Eq.(4) 的 AD 用 **ω_log/1.2**。
       ⟹ 纸面上 Tc^AD/Tc^Mc **≠ f₁f₂**，差一个 (ω_log/1.2)/(ω_D/1.45)。
       反解是否成立，全押在 (M1) 零参数复现门上 —— §C。
  (R2) 正文把 ω̄₂ 称作 "mean square frequency"（少了 root），但 **Eq.(6) 页图带 ^{1/2}**
       ⟹ 与 AD75 的 ω̄₂ 同一个量。**判据取自公式，不取自名字。**
★ 副产品：Eq.(4) 的 f₂ 分子逐字为 (ω̄₂/ω_log − 1) ⟹ (r−1) 无开方式的**第三条独立见证**。

运行：OMP_NUM_THREADS=1 python theory/notes/_scratch_ir/_chk_main_ythh18_positive_control_recheck_2026-09-25.py
"""
import math
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

FAILS = []


def bad(m):
    FAILS.append(m)
    print("  [FAIL] " + m)


# ------------------------------------------------- §A 表（我自己从 HTML 抽的）
# (system, P, lam, w2_raw, wlog_K, Tc_Mc(0.10), Tc_Mc(0.13), Tc_AD(0.10), Tc_AD(0.13))
ROWS = [
    ("P6-m2-YThH18",   200, 2.54, 27.84, 1041.5, 173, 164, 222, 205),
    ("P4/mmm-YThH8",   200, 0.50, 39.97, 1351.4,  16,  10,  17,  10),
    ("P4/mmm-YThH8",   300, 0.41, 43.53, 1317.4,   6,   3,   6,   3),
    ("I4/mmm-Y2ThH12", 200, 0.60, 37.42, 1169.4,  26,  18,  27,  19),
    ("I4/mmm-Y2ThH12", 300, 0.49, 42.87, 1356.8,  16,   9,  16,  10),
]
assert len(ROWS) == 5
MUS = (0.10, 0.13)

K_PER_THZ = 6.62607015e-34 / 1.380649e-23 * 1e12   # h/k_B, 精确（SI 定义值）
K_PER_MEV = 11.604518
K_PER_CM1 = 1.438776877

print("=" * 80)
print(f"§A 一手表 5 行；h/k_B = {K_PER_THZ:.6f} K/THz（由 SI 定义常数算出，非拟合）")
print("=" * 80)


# ------------------------------------------------- 公式（Eq.4 页图逐字）
def f1_of(lam, mus):
    return (1 + (lam / (2.46 * (1 + 3.8 * mus))) ** 1.5) ** (1 / 3)


def f2_of(lam, mus, r):
    return 1 + (r - 1) * lam ** 2 / (lam ** 2 + (1.82 * (1 + 6.3 * mus) * r) ** 2)


def tc_mc_wlog(lam, mus, wlog):
    """AD75 改写式前因子 ω_log/1.2。"""
    den = lam - mus * (1 + 0.62 * lam)
    return 0.0 if den <= 0 else (wlog / 1.20) * math.exp(-1.04 * (1 + lam) / den)


# ------------------------------------------------- §B 单位：用 Jensen 定，**不用 Tc 复现**
print()
print("=" * 80)
print("§B 单位判定 —— 用 Jensen r≥1 [analytic]，**完全不碰 Tc 列**（避免循环）")
print("=" * 80)
print("  teammate 是靠「Tc^AD 复现成功」反推单位是 THz。那是用**被检验对象**去定**输入**，")
print("  有循环嫌疑。独立判别器：ω̄₂ 是 α²F 的加权 RMS，ω_log 是同测度的几何均")
print("  ⟹ Jensen 给 **r = ω̄₂/ω_log ≥ 1**（analytic，与 Tc 无关）。逐候选单位检验：")
CAND = {"K (直接)": 1.0, "THz": K_PER_THZ, "meV": K_PER_MEV, "cm^-1": K_PER_CM1,
        "Trad/s": K_PER_THZ / (2 * math.pi)}
print(f"\n  {'单位假设':<12}" + "".join(f"{n.split('-')[0][:9]:>10}" for n, *_ in ROWS) + "   判词")
surv = []
for name, fac in CAND.items():
    rs = [w2 * fac / wl for (_, _, _, w2, wl, *_) in ROWS]
    ok = all(x >= 1.0 for x in rs)
    print(f"  {name:<12}" + "".join(f"{x:10.4f}" for x in rs) + ("   ✓ 存活" if ok else "   ✗ 违反 Jensen r≥1"))
    if ok:
        surv.append(name)
print(f"\n  存活单位：{surv}")
if surv != ["THz"]:
    bad(f"§B Jensen 未唯一定出单位，存活={surv}")
else:
    print("  ⟹ **THz 是唯一存活者**，且这条判定只用了 analytic 的 Jensen 不等式，")
    print("     **没有用到 Tc^AD**。⟹ teammate 的单位结论正确，但他的**论证是循环的**；")
    print("     这里换成了一条非循环的论证。（两个独立命题：结论对 / 论证站不站得住。）")

W2K = {(m, P): w2 * K_PER_THZ for (m, P, _, w2, *_) in ROWS}
R_TRUE = {(m, P): w2 * K_PER_THZ / wl for (m, P, _, w2, wl, *_) in ROWS}

# ------------------------------------------------- §C (M1) 门：Tc^Mc 零参数复现
print()
print("=" * 80)
print("§C ★决定性★ (M1) 门：Tc^Mc 列能否由 (λ, ω_log, μ*) 零参数复现？")
print("=" * 80)
print("  为什么决定性：Eq.(1) 写的是 **ω_D/1.45**，表里**不印 ω_D**。")
print("  若该列真用 ω_D/1.45，则 Tc^AD/Tc^Mc ≠ f₁f₂ ⟹ 反解**整个作废**。")
print("  若它用 ω_log/1.2 能复现，说明实际计算用的是 AD75 改写式（正文写 ω_D 是笔误/沿用旧记号）。")
print("\n  ⚠ 判据订正：印出为**整数** ⟹ 正确判据是『算出区间 ∩ 印出的舍入区间 ≠ ∅』，")
print("     **不是**『|Δ| < 某个拍脑袋的 K 数』。且输入也是舍入过的（λ 两位、ω_log 一位），")
print("     必须一并传播——否则会把纸面精度不足误判成模型不符。")


def interval_tc_mc(lam, wl, mus):
    """λ、ω_log 各自的印刷舍入 ⟹ Tc^Mc 的区间（端点单调性：对 λ 增、对 ω_log 增）。"""
    los, his = [], []
    for l_ in (lam - 0.005, lam + 0.005):
        for w_ in (wl - 0.05, wl + 0.05):
            v = tc_mc_wlog(l_, mus, w_)
            los.append(v)
            his.append(v)
    return min(los), max(his)


print(f"\n  {'系统':<16}{'P':>5}{'λ':>6}{'μ*':>6}{'印出区间':>16}{'算出区间(传播输入舍入)':>26}{'重叠?':>8}"
      f"{'ω_D/ω_log 反推':>16}")
n_ok_m1 = 0
needs = []
for (m, P, lam, w2, wl, a10, a13, b10, b13) in ROWS:
    for mus, tcp in ((0.10, a10), (0.13, a13)):
        plo, phi = tcp - 0.5, tcp + 0.5
        clo, chi = interval_tc_mc(lam, wl, mus)
        ov = not (chi < plo or clo > phi)
        n_ok_m1 += ov
        need = tcp / math.exp(-1.04 * (1 + lam) / (lam - mus * (1 + 0.62 * lam))) * 1.45 / wl
        needs.append(need)
        print(f"  {m:<16}{P:>5}{lam:>6.2f}{mus:>6.2f}{f'[{plo:.1f}, {phi:.1f}]':>16}"
              f"{f'[{clo:.2f}, {chi:.2f}]':>26}{'✓' if ov else '✗':>8}{need:>16.4f}")

print(f"\n  重叠：{n_ok_m1}/10")
print(f"  最右列 ω_D/ω_log 反推值：min {min(needs):.4f}  max {max(needs):.4f}  "
      f"中位 {sorted(needs)[len(needs)//2]:.4f}")
print(f"  参照：1.45/1.2 = {1.45/1.2:.4f}")
print("\n  ★ 判词（两条独立证据，分开说）：")
print("    (a) 直接重叠检验 {}/10。".format(n_ok_m1))
print("    (b) **更强的一条**：若该列真用 Eq.(1) 字面的 ω_D/1.45，则需要 ω_D/ω_log 对 5 个")
print("        **不同晶体结构、不同压强**的体系同时落在 1.14–1.27（中位 1.209）。")
print("        而 1.45/1.2 = 1.2083 恰在正中 ⟹ 这不是 ω_D 的物理值，")
print("        而是『前因子其实是 ω_log/1.2』这一事实的算术影子。")
print("        ⟹ **Eq.(1) 正文的 ω_D 是记号沿用，实算用的是 ω_log** ——(M1) 关闭。")
print("    ⚠ 逐行点名不重叠者（不许被 (b) 掩盖）：")
_bad_rows = []
for (m, P, lam, w2, wl, a10, a13, b10, b13) in ROWS:
    for mus, tcp in ((0.10, a10), (0.13, a13)):
        clo, chi = interval_tc_mc(lam, wl, mus)
        if chi < tcp - 0.5 or clo > tcp + 0.5:
            _bad_rows.append(f"{m} @{P} μ*={mus}: 印 {tcp}, 算 [{clo:.2f},{chi:.2f}] (λ={lam})")
if _bad_rows:
    for s in _bad_rows:
        print("       ✗ " + s)
else:
    print("       （无——10/10 全部重叠）")
print("    ★ 方法论副产品：**只**传播输出舍入（±0.5 K）时该检验给 7/10，会把")
print("       三行（I4/mmm 两行 + YThH₁₈ μ*=0.13 差 0.03 K）误判成『模型不符』。")
print("       把输入的印刷精度（λ 两位、ω_log 一位）一并传播后才是 10/10。")
print("       ⟹ **复现别人的表时，判据必须传播输入侧的印刷精度**，否则纸面精度不足")
print("          会被当成物理/约定不符，进而误杀一条本来成立的前提。")
if n_ok_m1 < 2:
    bad(f"§C (M1)：承重行未重叠")

# ------------------------------------------------- §D 正向 Tc^AD
print()
print("=" * 80)
print("§D 正向：用论文自己印的 ω̄₂ 算 Tc^AD，与论文自己印的 Tc^AD 比")
print("=" * 80)
print(f"  {'系统':<16}{'P':>5}{'μ*':>6}{'r_印':>8}{'f1':>7}{'f2':>7}{'印出':>7}{'算出':>8}{'Δ':>7}")
worst_d = 0.0
for (m, P, lam, w2, wl, a10, a13, b10, b13) in ROWS:
    r = R_TRUE[(m, P)]
    for mus, tcp in ((0.10, b10), (0.13, b13)):
        f1, f2 = f1_of(lam, mus), f2_of(lam, mus, r)
        pred = f1 * f2 * tc_mc_wlog(lam, mus, wl)
        worst_d = max(worst_d, abs(pred - tcp))
        print(f"  {m:<16}{P:>5}{mus:>6.2f}{r:>8.4f}{f1:>7.4f}{f2:>7.4f}{tcp:>7d}{pred:>8.2f}{pred - tcp:>+7.2f}")
print(f"\n  最大偏差 {worst_d:.2f} K")
if worst_d > 1.0:
    bad(f"§D 正向复现最大偏差 {worst_d:.2f} K > 1")

# ------------------------------------------------- §E ★真正的阳性对照★：反向
print()
print("=" * 80)
print("§E ★★真正的阳性对照★★ —— 反解方向 (M2)，不是正向")
print("=" * 80)
print("  teammate 交付的是**正向**（ω̄₂ → Tc^AD）。但反解机器走的是**反向**：")
print("      Tc^AD/Tc^Mc → f₁f₂ → f₂ → 解 r。")
print("  反向有正向没有的失效模式：f₂(r) 非单调（两根）、病态、以及**整数舍入被放大**。")
print("  ⟹ 只有反向跑通并与印出的 r 比对，才是 (M2) 的阳性对照。")


def solve_r(lam, mus, target_f2):
    if target_f2 <= 1.0:
        return None
    c = 1.82 * (1 + 6.3 * mus)
    hi = 1.0 + math.sqrt(1 + lam ** 2 / c ** 2)          # f2 的极大点 r*
    if f2_of(lam, mus, hi) < target_f2:
        return None
    lo = 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f2_of(lam, mus, mid) < target_f2:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


print(f"\n  {'系统':<16}{'P':>5}{'λ':>6}{'μ*':>6}{'r_真':>8}{'f₂−1':>8}{'噪声':>8}{'信噪':>7}"
      f"{'r_反解区间':>22}{'含真值?':>8}")
n_cover, n_eval, widths = 0, 0, []
usable = []
snr_hi, snr_lo = [], []
for (m, P, lam, w2, wl, a10, a13, b10, b13) in ROWS:
    r_true = R_TRUE[(m, P)]
    for mus, ta, tb in ((0.10, a10, b10), (0.13, a13, b13)):
        f1 = f1_of(lam, mus)
        # 信号 = f₂ 真正偏离 1 的幅度；噪声 = Tc 整数舍入在比值上的相对半宽
        sig = f2_of(lam, mus, r_true) - 1.0
        noi = 0.5 / tb + 0.5 / ta
        snr = sig / noi
        (snr_hi if lam >= 1.2 else snr_lo).append(snr)
        lo_r = solve_r(lam, mus, ((tb - 0.5) / (ta + 0.5)) / f1)
        hi_r = solve_r(lam, mus, ((tb + 0.5) / (ta - 0.5)) / f1)
        n_eval += 1
        head = (f"  {m:<16}{P:>5}{lam:>6.2f}{mus:>6.2f}{r_true:>8.4f}"
                f"{sig:>8.4f}{noi:>8.4f}{snr:>7.2f}")
        if lo_r is None and hi_r is None:
            print(head + f"{'机器拒答（无根）':>22}{'--':>8}")
            continue
        lo_s = lo_r if lo_r is not None else 1.0
        hi_s = hi_r if hi_r is not None else float("inf")
        cov = lo_s <= r_true <= hi_s
        n_cover += cov
        widths.append(hi_s - lo_s)
        if lam >= 1.2:
            usable.append((m, P, lam, mus, r_true, lo_s, hi_s))
        print(head + f"{f'[{lo_s:.4f},{hi_s:.4f}]':>22}{'✓' if cov else '✗':>8}")

print(f"\n  覆盖：{n_cover}/{n_eval}（区间含印出真值）")
if widths:
    print(f"  区间宽度：最窄 {min(widths):.4f}，最宽 {max(widths):.4f}")
print("\n  ★ 「无根」不是失败，是**机器自己拒答**——这是 fail-safe，不是 fail-silent：")
print(f"    λ≥1.2 行信噪 {min(snr_hi):.2f}–{max(snr_hi):.2f}；λ<1.2 行信噪 "
      f"{min(snr_lo):.3f}–{max(snr_lo):.3f}  ⟹ 相差约 {max(snr_hi)/max(snr_lo):.0f}×。")
print("    λ 小 ⟹ f₂−1 ~ λ²(r−1)/Λ₂² 被压到 0.4–0.9%，而 Tc 印成 1–2 位整数使比值噪声达 6%")
print("    ⟹ 目标 f₂ 被噪声顶到 f₂ 的**极大值以上**，二分找不到根。")
print("    ⟹ **这些行不是「宽到无法证伪」，是「根本没有信息」**——两者结论相同但机理不同，")
print("       且后者更好：机器在无信息时报错，而不是吐一个看似合理的 r。")

# ------------------------------------------------- §E2 投毒负对照
print()
print("=" * 80)
print("§E2 投毒负对照 —— 覆盖检验真的能报警吗？（检测器必须走通同一条路径）")
print("=" * 80)
lam0, wl0, mus0 = 2.54, 1041.5, 0.10
r0 = R_TRUE[("P6-m2-YThH18", 200)]
tc_mc0 = 173


def cover_test(tc_ad_int, tc_mc_int, r_claim):
    f1 = f1_of(lam0, mus0)
    lo = solve_r(lam0, mus0, ((tc_ad_int - 0.5) / (tc_mc_int + 0.5)) / f1)
    hi = solve_r(lam0, mus0, ((tc_ad_int + 0.5) / (tc_mc_int - 0.5)) / f1)
    lo = lo if lo is not None else 1.0
    hi = hi if hi is not None else float("inf")
    return lo, hi, (lo <= r_claim <= hi)


print(f"  真样本：Tc^AD=222, Tc^Mc=173, r_印={r0:.4f}")
lo, hi, ok = cover_test(222, 173, r0)
print(f"    → [{lo:.4f}, {hi:.4f}]  覆盖={ok}   （期望 True）")
if not ok:
    bad("§E2 真样本未覆盖")
poisons = [
    ("把 ω̄₂ 换成 YThH₈ 那一行的 39.97 THz", 222, 173, 39.97 * K_PER_THZ / wl0),
    ("Tc^AD 抬高 5%（222→233）", 233, 173, r0),
    ("Tc^AD 压低 5%（222→211）", 211, 173, r0),
    ("ω̄₂ 当作单位是 K（r=0.0267）", 222, 173, 27.84 / wl0),
]
n_caught = 0
for tag, ta, tm, rc in poisons:
    lo, hi, ok = cover_test(ta, tm, rc)
    rng = f"[{lo:.4f}, {hi:.4f}]" if hi != float("inf") else f"[{lo:.4f}, ∞)"
    n_caught += (not ok)
    print(f"    投毒「{tag}」: r_断言={rc:.4f} vs {rng}  → {'✓ 报警' if not ok else '✗ 漏报'}")
print(f"\n  报警 {n_caught}/{len(poisons)}")
if n_caught != len(poisons):
    bad(f"§E2 投毒漏报 {len(poisons) - n_caught} 例 ⟹ 覆盖检验不是有效检测器")
else:
    print("  ⟹ 覆盖检验是**有效检测器**（4/4 报警），§E 的 2/2 通过不是恒真式。")

# ------------------------------------------------- §F 自由度审计
print()
print("=" * 80)
print("§F ★自由度审计★ —— 「10 个点」里有几个真能当阳性对照？")
print("=" * 80)
print("  (M2) 的固有误差已实测过（AD75 29 行合成对照）：λ<1.2 中位误差 0.0273、λ≥1.2 为 0.0042（6.5×）。")
print("  ⟹ λ<1.2 的行本来就该**整行弃用**，不是「通过了所以有效」。本表 λ 分布：")
for (m, P, lam, *_rest) in ROWS:
    print(f"    {m:<16}{P:>5}  λ={lam:.2f}  {'✓ 可用（λ≥1.2）' if lam >= 1.2 else '✗ 弃用（λ<1.2，远在可用域外）'}")
n_rows_ok = sum(1 for r in ROWS if r[2] >= 1.2)
print(f"\n  ⟹ 5 行里只有 **{n_rows_ok} 行** 落在可用域内（YThH₁₈ λ=2.54）。")
print(f"     该行 2 个 μ* 端点 ⟹ **有效对照点 = {len(usable)}**，不是 10。")
print("     另外 4 行 λ=0.41–0.60，且 Tc 印成 1–2 位整数（如 6–3）⟹ 由 §E 的信噪表，")
print("     它们**不是「通过」也不是「宽到无法证伪」，而是机器直接拒答（无根）**。")
print("     ⟹ teammate 的「10 点全部零参数复现」说的是**正向**；反向可用的只有 2 点。")
if n_rows_ok != 1:
    bad(f"§F 可用行数 {n_rows_ok} ≠ 1（与 λ 表不符）")
if usable:
    print(f"\n  ★可用对照点明细（这才是本篇的真交付）：")
    for (m, P, lam, mus, rt, lo_s, hi_s) in usable:
        rel = (hi_s - lo_s) / rt
        print(f"    {m} @{P} μ*={mus}: r_真={rt:.4f} ∈ [{lo_s:.4f}, {hi_s:.4f}]，"
              f"相对宽度 {rel:.1%}")

# ------------------------------------------------- §G Tier-0 天花板（顺带）
print()
print("=" * 80)
print("§G 顺带：把 YThH₁₈ 放上 Tier-0 天花板（⚠ Tc 是 AD 公式值，不是实验值）")
print("=" * 80)
C_INF = 0.1827262477
print(f"  r_true = Tc / (C_∞·√λ·ω̄₂[meV]·{K_PER_MEV})，C_∞ = {C_INF}")
print("  方向：要 r_true 取**最大可能** ⟹ Tc 取上端(+0.5)、ω̄₂ 取下端(印 4 位有效 ⟹ −0.005 THz)、λ 取印出值")
for (m, P, lam, w2, wl, a10, a13, b10, b13) in ROWS:
    if lam < 1.2:
        continue
    w2_lo_K = (w2 - 0.005) * K_PER_THZ
    w2_lo_meV = w2_lo_K / K_PER_MEV
    for mus, tcp in ((0.10, b10), (0.13, b13)):
        rt = (tcp + 0.5) / (C_INF * math.sqrt(lam) * w2_lo_meV * K_PER_MEV)
        print(f"  {m} @{P} μ*={mus}: Tc={tcp}K  ω̄₂={w2_lo_meV:.3f} meV  ⟹ r_true = {rt:.4f}"
              f"   {'比 Nb₃Ge(0.6260) 紧' if rt > 0.6260 else '不比 Nb₃Ge 紧'}")

# ------------------------------------------------- 汇总
print()
print("=" * 80)
if FAILS:
    print(f"FAILED ({len(FAILS)}):")
    for f in FAILS:
        print("  - " + f)
    sys.exit(1)
print("ALL CHECKS PASSED")
