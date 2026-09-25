# -*- coding: utf-8 -*-
"""
从零独立复核 arXiv:2603.11590 (Shahid et al., La-Zr-H hydrides) 能否当反解机器的阳性对照。

起因：teammate 会话 `mcad-omega2-positive-control` 报该文「同印 Tc^Mc/Tc^AD + 直接印 ω̄₂」
⟹ 疑似 6 个对照点（3 行 λ≥1.2 × 2 个 μ* 端点）。它的 WebFetch 对所有域名前置校验失败，
全部证据自标 L2 ⟹ 本脚本不采信任何一条，只用我自己经 ar5iv 全 TLS 校验取得的一手 LaTeXML 源。

一手产物：theory/notes/_scratch_ir/_lazrh_2603.11590_ar5iv.html

严格性：全篇 numerical-evidence（float64 复现一张印刷表）。无 certificate 级结论。

节：
  §0 产物与逐字公式
  §A 转录 + Eq.(10) 自洽（零参数，校验 ω_ph / ω_max^Mats 的转录）
  §B 单位与 Jensen r≥1、ω̄₂ ≤ ω_max
  §C (M1a) Tc^Mc 列零参数复现（ω_log/1.2, f1f2=1），传播输入侧印刷精度
  §D (M1b) Tc^AD 列零参数复现（AD75 正确 f2）
  §E 变体模型甄别（预先列全，不许挑）
  §F 反解可行性：所需 f2 vs f2 的极大值 ⟹ 机器给不给根
  §G Tier-0 天花板：对 Tc(Eliashberg) 与 Tc(AD)
  §H 投毒阴性对照
  §I 自由度审计
"""
import hashlib
import math
import os
import sys

FAIL = []


def chk(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAIL.append(name)
    print(f"  [{tag}] {name}" + (f"  {detail}" if detail else ""))


HERE = os.path.dirname(os.path.abspath(__file__))
HTML = os.path.join(HERE, "_lazrh_2603.11590_ar5iv.html")

# ---------------------------------------------------------------- §0
print("=" * 78)
print("§0  一手产物与逐字公式")
print("=" * 78)
h = hashlib.md5(open(HTML, "rb").read()).hexdigest()
n = os.path.getsize(HTML)
print(f"  {os.path.basename(HTML)}  size={n}  md5={h}")
print("""
  取得路线：https://ar5iv.labs.arxiv.org/html/2603.11590   (http=200, ssl_verify_result=0)
  ★ 未放松任何校验位（无 -k / --insecure / --ssl-no-revoke / unverified_context）。
     本机 arxiv.org 被 TLS 劫持，ar5iv 是已验干净路线。

  逐字（LaTeXML 源，非转述）：
  Eq.(1)  T_c = f_1 f_2 \\frac{\\omega_{log}}{1.2}
                  \\exp(\\frac{-1.04(1+\\lambda)}{\\lambda-\\mu^{*}(1+0.62\\lambda)})
          ⟹ 前因子是 ω_log/1.2（AD 改写式），**不是** McMillan 原文的 Θ_D/1.45。
             与 PMC11656362 恰好相反（那篇 Eq.(1) 逐字写 ω_D/1.45）⟹ (R1) 风险在本篇不存在。
  Eq.(2)  f_1 = [1+(\\frac{\\lambda}{2.46(1+3.8\\mu^{*})})^{3/2}]^{1/3}
  Eq.(3)  f_2 = 1 + \\frac{(\\frac{\\bar{\\omega_2}}{\\omega_{log}}-1)\\lambda^2}
                        {\\lambda^2+[1.82(1+6.3\\mu^{*})\\frac{\\bar{\\omega_2}}{\\omega_{log}}]^2}
          ⟹ 分子 (r-1)，无平方、无开方，`1+` 在分数外，Λ₂ **含** r 因子。
             = AD75 Eq.(35) 逐字 ⟹ **第四条独立第三方见证**。
  Eq.(6)  \\bar{\\omega}_{2} = \\sqrt{\\frac{1}{\\lambda}\\int_0^{\\omega_{max}}
                              [\\frac{2\\alpha^2 F(\\omega)}{\\omega}]\\omega^2 d\\omega}
          ⟹ 带 \\sqrt，被积函数 [2α²F/ω]·ω² = 2α²F·ω ⟹ 与 AD75 **逐字同一个量**。
             (R2) 关闭：正文措辞 "second moment"（少 root）不算证据，判据取自公式。
  Eq.(10) 1/\\mu^{*}_{Eliashberg} = 1/\\mu^{*}_{McMillan} + \\ln(\\omega_{ph}/\\omega^{Matsub}_{max})
""")

# ---------------------------------------------------------------- §A
print("=" * 78)
print("§A  转录 + Eq.(10) 自洽（零参数）")
print("=" * 78)

# 表列：名称, Tc_Mc(hi,lo), Tc_AD(hi,lo), Tc_ME(hi,lo), mus_E(hi,lo),
#       w_ph[THz], w_max_Mats(hi,lo)[THz], lam, wlog[K], w2[K]
ROWS = [
    ("R3m-Zr2H17@300",      (129, 115), (188, 158), (209, 178), (0.13, 0.21),
     70.29, (590.96, 496.62), 2.11,  841, 1613),
    ("P6b-m2-LaZrH18@300",  (151, 133), (193, 162), (206, 173), (0.13, 0.22),
     67.27, (606.67, 509.19), 1.84, 1111, 1690),
    ("P6/mmm-LaZr2H24@200", (143, 129), (191, 163), (202, 175), (0.13, 0.22),
     62.79, (600.38, 512.33), 2.20,  941, 1348),
    ("I4/mmm-La2ZrH12@250", (90, 69),   (98, 74),   (97, 75),   (0.12, 0.18),
     67.27, (308.01, 232.51), 1.00, 1283, 1788),
]
MUS = (0.10, 0.15)   # μ*_McMillan 行：全四列同一个 "0.10 - 0.15"，Mc 与 AD 两列共用
K_PER_THZ = 6.62607015e-34 / 1.380649e-23 * 1e12   # h/k_B, SI 定义值
print(f"  h/k_B = {K_PER_THZ:.6f} K/THz")
print("""
  ★ 共用 μ* 已由一手源钉死（teammate 把这条列为「最大风险」，实为假）：
    表里 μ*_McMillan 是**单独一行、四个材料列取同一个 "0.10 - 0.15"**，
    正文 "For the McMillan and Allen-Dynes approach, standard Coulomb pseudopotential
    values of mu*_McMillan = 0.1 ..." ⟹ Mc 与 AD **同 μ***。
    0.10/0.15 是 Tc 区间的两个端点（μ* 大 ⟹ Tc 小 ⟹ 对应 lo）。
    μ*_Eliashberg(0.12-0.22) 是经 Eq.(10) **导出**的另一个量，只供 Eliashberg 列。
""")
print(f"  {'material':22s} {'mu*':>5s} {'1/muE calc':>10s} {'muE calc':>9s} {'printed':>8s}")
for (nm, _tmc, _tad, _tme, musE, wph, wmax, lam, wl, w2) in ROWS:
    for j, ms in enumerate(MUS):
        val = 1.0 / ms + math.log(wph / wmax[j])
        got = 1.0 / val
        print(f"  {nm:22s} {ms:5.2f} {val:10.5f} {got:9.5f} {musE[j]:8.2f}")
        chk(f"Eq.(10) {nm} mu*={ms}", abs(got - musE[j]) <= 0.005 + 1e-9,
            f"|{got:.5f}-{musE[j]:.2f}|={abs(got-musE[j]):.5f} (印两位 ⟹ ±0.005)")
print("""
  ⟹ Eq.(10) 在 8/8 上自洽到印刷精度。这一步的价值不是「验证作者会做除法」，
     而是**零参数地校验我对 ω_ph、ω_max^Matsub、μ*_E 三组共 20 个数的转录**
     ——它们通过一个我没参与拟合的恒等式互相锁死。
""")

# ---------------------------------------------------------------- §B
print("=" * 78)
print("§B  单位与 Jensen")
print("=" * 78)
print("  ω̄₂ 与 ω_log **都印为 K**（表头逐字 `\\bar{\\omega}_{2}, K` / `\\omega_{\\log}, K`）")
print("  ⟹ r := ω̄₂/ω_log 无量纲，**不需要任何单位推断**（与 PMC11656362 相反）。")
print(f"\n  {'material':22s} {'lam':>5s} {'wlog[K]':>8s} {'w2[K]':>7s} {'r':>8s} {'w_max[K]':>9s}")
for (nm, _a, _b, _c, _d, wph, _w, lam, wl, w2) in ROWS:
    r = w2 / wl
    wmaxK = wph * K_PER_THZ
    print(f"  {nm:22s} {lam:5.2f} {wl:8d} {w2:7d} {r:8.4f} {wmaxK:9.1f}")
    chk(f"Jensen r>=1 {nm}", r >= 1.0, f"r={r:.4f}")
    chk(f"w2<=w_max {nm}", w2 <= wmaxK, f"{w2} <= {wmaxK:.1f} K")
print("""
  Jensen [analytic]：ω̄₂ 是 α²F 归一权重下的**均方根**，ω_log 是同一测度的**几何平均**
  ⟹ r ≥ 1 恒成立。4/4 满足。⚠ 这不是一次强检验（本篇单位已印出，无需推断）；
     它在这里只起转录自查作用。
  ★ r = 1.9180（Zr₂H₁₇）**超过** AD75 Table I 29 种真实材料的最大值 1.8364
     ⟹ 若属实，是语料里形状最「铺开」的谱之一；也提高了该行印刷值出错的先验。
""")


# ------------------------------------------------ AD 公式族
def f1_ad(lam, mus):
    return (1.0 + (lam / (2.46 * (1.0 + 3.8 * mus))) ** 1.5) ** (1.0 / 3.0)


def f2_ad(lam, mus, r):
    """AD75 Eq.(35) = 本文 Eq.(3) 逐字。"""
    L2 = 1.82 * (1.0 + 6.3 * mus) * r
    return 1.0 + (r - 1.0) * lam ** 2 / (lam ** 2 + L2 ** 2)


def tc_base(lam, mus, wl):
    """Eq.(1) 去掉 f1f2 的部分 = 假设的 Tc^Mc。"""
    den = lam - mus * (1.0 + 0.62 * lam)
    if den <= 0:
        return 0.0
    return (wl / 1.20) * math.exp(-1.04 * (1.0 + lam) / den)


def interval(fn, lam, wl, mus, dlam=0.005, dwl=0.5, extra=None):
    """传播输入侧印刷精度（角点扫；被复现函数对 λ、ω_log 单调 ⟹ 角点即够）。"""
    vals = []
    ex = extra if extra is not None else [None]
    for l_ in (lam - dlam, lam + dlam):
        for w_ in (wl - dwl, wl + dwl):
            for e_ in ex:
                vals.append(fn(l_, mus, w_) if e_ is None else fn(l_, mus, w_, e_))
    return min(vals), max(vals)


# ---------------------------------------------------------------- §C
print("=" * 78)
print("§C  (M1a) Tc^Mc 列零参数复现 —— 前因子 ω_log/1.2、f1f2=1")
print("=" * 78)
print("  判据：[算出区间] ∩ [印出±0.5] ≠ ∅ ，且**输入侧印刷精度也传播**")
print("        （λ 印两位 ⟹ ±0.005；ω_log 印整数 ⟹ ±0.5）")
print("        ——只传输出舍入会给假 FAIL，见 memory 的同名规则。\n")
print(f"  {'material':22s} {'mu*':>5s} {'calc':>8s} {'[lo':>8s} {'hi]':>8s} {'print':>6s} {'rel%':>7s} ov")
mc_bad = []
for (nm, tmc, _b, _c, _d, _e, _f, lam, wl, w2) in ROWS:
    for j, ms in enumerate(MUS):
        c = tc_base(lam, ms, wl)
        lo, hi = interval(tc_base, lam, wl, ms)
        p = tmc[j]
        ov = not (hi < p - 0.5 or lo > p + 0.5)
        rel = 100.0 * (p - c) / c
        print(f"  {nm:22s} {ms:5.2f} {c:8.2f} {lo:8.2f} {hi:8.2f} {p:6d} {rel:+7.2f} {'Y' if ov else 'N'}")
        if not ov:
            mc_bad.append((nm, ms, c, p, rel))
print(f"\n  重叠 {8-len(mc_bad)}/8")
if mc_bad:
    print("  ✗ 不重叠的行：")
    for (nm, ms, c, p, rel) in mc_bad:
        print(f"      {nm} mu*={ms}: 算 {c:.2f} vs 印 {p}  ({rel:+.2f}%)")
    print("""
  ★ 两类失效必须分清（幅度 + 双端点符号一致性）：
    · **结构性**：Zr₂H₁₇ 两个 μ* 端点 **同号 +2.7/+2.9%**，是容差（≈±0.45%）的 6 倍
      ⟹ 不可能由印刷精度产生 ⟹ 该行印出的某个数与其余列不自洽。
      方向审：印出 **高于** 算出 ⟹ 比值 Tc^AD/Tc^Mc 被**压低** ⟹ 反解 f₂ 偏小
      ⟹ r 偏小 ⟹ r_true = r_up/r 偏**大** ⟹ **偏向「夺冠」** ⟹ 该行绝不可照用。
      （与 LaH₁₀@163 的 Mc 偏低 −3.30% 同形反号，但下游偏向相同。）
    · **边缘性**：LaZr2H24 只 μ*=0.10 端点差 **−0.69%**（≈1.3 倍容差），
      另一端点通过，且**符号与上一类相反**（算出高于印出）
      ⟹ 与「作者内部用了比印刷位数更多的有效位」相容，不构成该行不自洽的证据。
      方向审：印出偏低 ⟹ 比值偏高 ⟹ r 偏大 ⟹ r_true 偏**小** ⟹ **保守侧**。""")
    # 单一数字修正假说：ω_log misprint？
    print("\n  单数修正假说（**假说，不是结论**）：若该行 ω_log 被误印，反解所需 ω_log =")
    for (nm, ms, c, p, rel) in mc_bad:
        row = [x for x in ROWS if x[0] == nm][0]
        need = row[8] * p / c
        print(f"      {nm} mu*={ms}:  ω_log ≈ {need:.1f} K （印出 {row[8]}）")
    print("      两端点互相一致 ⟹ 是**乘性**偏差，与 μ* 无关 ⟹ 指向前因子/ω_log 而非 μ* 约定。")
    print("      ⚠ 我无法从这篇判定它错在哪（α²F 只有 Fig.5 图形，非表格）⟹ 只报「该行不自洽」。")

# ---------------------------------------------------------------- §D
print("\n" + "=" * 78)
print("§D  (M1b) Tc^AD 列零参数复现 —— AD75 正确 f2（Eq.3 逐字）")
print("=" * 78)


def tc_ad_correct(lam, mus, wl, w2):
    r = w2 / wl
    return tc_base(lam, mus, wl) * f1_ad(lam, mus) * f2_ad(lam, mus, r)


print(f"  {'material':22s} {'mu*':>5s} {'calc':>8s} {'[lo':>8s} {'hi]':>8s} {'print':>6s} {'rel%':>7s} ov")
ad_bad = []
for (nm, _a, tad, _c, _d, _e, _f, lam, wl, w2) in ROWS:
    for j, ms in enumerate(MUS):
        c = tc_ad_correct(lam, ms, wl, w2)
        vals = []
        for l_ in (lam - 0.005, lam + 0.005):
            for w_ in (wl - 0.5, wl + 0.5):
                for x_ in (w2 - 0.5, w2 + 0.5):
                    vals.append(tc_ad_correct(l_, ms, w_, x_))
        lo, hi = min(vals), max(vals)
        p = tad[j]
        ov = not (hi < p - 0.5 or lo > p + 0.5)
        rel = 100.0 * (p - c) / c
        print(f"  {nm:22s} {ms:5.2f} {c:8.2f} {lo:8.2f} {hi:8.2f} {p:6d} {rel:+7.2f} {'Y' if ov else 'N'}")
        if not ov:
            ad_bad.append((nm, ms, c, p, rel))
print(f"\n  重叠 {8-len(ad_bad)}/8")
print("""
  ★★ 这是本篇的判决：**AD 列与 Mc 列不由 AD75 的 f₁f₂ 相联**。
     印出的 AD 值**系统性高于**用同一 (λ, ω_log, ω̄₂, μ*) 算出的 AD 值。
     注意这一条与 §C 独立：§C 只测 Mc 列，这里测的是 f₁f₂ 这个箭头。
""")

# ---------------------------------------------------------------- §E
print("=" * 78)
print("§E  变体模型甄别（候选**预先列全**，逐个报，不许挑）")
print("=" * 78)


def variants():
    c0 = lambda ms: 1.82 * (1.0 + 6.3 * ms)
    V = {}
    V["(a) AD75 正确: L2=1.82(1+6.3u)r"] = lambda lam, ms, r: 1 + (r - 1) * lam**2 / (lam**2 + (c0(ms) * r)**2)
    V["(b) L2 丢掉 r 因子"] = lambda lam, ms, r: 1 + (r - 1) * lam**2 / (lam**2 + c0(ms)**2)
    V["(c) 分子 (r^2-1)"] = lambda lam, ms, r: 1 + (r*r - 1) * lam**2 / (lam**2 + (c0(ms) * r)**2)
    V["(d) 整体开方"] = lambda lam, ms, r: 1 + math.sqrt(max(0.0, (r - 1)) * lam**2 / (lam**2 + (c0(ms) * r)**2))
    V["(e) L2 用 1/r"] = lambda lam, ms, r: 1 + (r - 1) * lam**2 / (lam**2 + (c0(ms) / r)**2)
    V["(f) L2 无 mu* 项: L2=1.82r"] = lambda lam, ms, r: 1 + (r - 1) * lam**2 / (lam**2 + (1.82 * r)**2)
    V["(g) f2 = r (完全饱和)"] = lambda lam, ms, r: r
    V["(h) f2 = 1 (无形状修正)"] = lambda lam, ms, r: 1.0
    return V


V = variants()
print(f"  {'变体':34s} {'max|rel|%':>10s} {'中位|rel|%':>11s} {'<=0.7% 的点':>12s}")
scores = []
for key, f2v in V.items():
    errs = []
    for (nm, _a, tad, _c, _d, _e, _f, lam, wl, w2) in ROWS:
        r = w2 / wl
        for j, ms in enumerate(MUS):
            c = tc_base(lam, ms, wl) * f1_ad(lam, ms) * f2v(lam, ms, r)
            errs.append(abs(tad[j] - c) / c * 100.0)
    errs.sort()
    med = 0.5 * (errs[3] + errs[4])
    ngood = sum(1 for e in errs if e <= 0.7)
    scores.append((max(errs), med, ngood, key))
    print(f"  {key:34s} {max(errs):10.2f} {med:11.2f} {ngood:12d}")
scores.sort()
print(f"\n  最优变体：{scores[0][3]}   max|rel| = {scores[0][0]:.2f}%")
chk("变体甄别有唯一赢家（次优的 max|rel| 至少 2 倍）",
    scores[1][0] >= 2.0 * scores[0][0],
    f"best={scores[0][0]:.2f}% next={scores[1][0]:.2f}% ({scores[1][3]})")
print("""
  ★ 逐行看赢家 (b)：""")
best = V[scores[0][3]]
print(f"  {'material':22s} {'mu*':>5s} {'calc':>8s} {'print':>6s} {'rel%':>7s}")
for (nm, _a, tad, _c, _d, _e, _f, lam, wl, w2) in ROWS:
    r = w2 / wl
    for j, ms in enumerate(MUS):
        c = tc_base(lam, ms, wl) * f1_ad(lam, ms) * best(lam, ms, r)
        print(f"  {nm:22s} {ms:5.2f} {c:8.2f} {tad[j]:6d} {100.0*(tad[j]-c)/c:+7.2f}")
print("""
  判词 [numerical-evidence]：**印出的公式不是跑出那张表的公式。**
  Eq.(3) 逐字含 Λ₂ = 1.82(1+6.3μ*)·r；而复现那张表的是 Λ₂ = 1.82(1+6.3μ*)（丢了 r 因子）。
  方向自洽：丢掉 r>1 ⟹ Λ₂ 变小 ⟹ f₂ 变大 ⟹ Tc 变高 = §D 观测到的偏差方向。
  ⚠ 这是对第三方实现的**推断**，不是证明；我只能说「(b) 在预先列出的 8 个候选里唯一拟合」。
     Zr₂H₁₇ 行仍是残差最大的一行（与 §C 同一行）⟹ **两个互不相干的诊断指向同一行**。
""")

# ---------------------------------------------------------------- §F
print("=" * 78)
print("§F  反解可行性：所需 f2 vs f2 在 r 上的极大值")
print("=" * 78)
print("  f2(r) 非单调，极大在 r* = 1 + sqrt(1 + λ²/c²), c = 1.82(1+6.3μ*)")
print("  ⟹ 若「所需 f2」> f2(r*)，二分**无根**，机器**拒答**（fail-safe，不是 fail-silent）\n")
print(f"  {'material':22s} {'mu*':>5s} {'req f1f2':>9s} {'req f2':>8s} {'f2max':>8s} {'r*':>7s} 可解")
nsolv = 0
for (nm, tmc, tad, _c, _d, _e, _f, lam, wl, w2) in ROWS:
    for j, ms in enumerate(MUS):
        c0 = 1.82 * (1.0 + 6.3 * ms)
        rstar = 1.0 + math.sqrt(1.0 + lam ** 2 / c0 ** 2)
        f2max = f2_ad(lam, ms, rstar)
        req = tad[j] / tmc[j]
        reqf2 = req / f1_ad(lam, ms)
        ok = reqf2 <= f2max
        if ok:
            nsolv += 1
        print(f"  {nm:22s} {ms:5.2f} {req:9.5f} {reqf2:8.5f} {f2max:8.5f} {rstar:7.4f} {'Y' if ok else 'N'}")
print(f"\n  可解 {nsolv}/8")
chk("机器对本篇 8/8 拒答（与 §D 独立地拒收）", nsolv == 0,
    f"可解 {nsolv} 个")
print("""
  ★★ 这一节的意义：反解机器**自己**就把这篇挡住了，不需要我先做 §D 的诊断。
     「所需 f₂ 超过 f₂ 的解析极大值」是一个**内建的前提自检**
     ⟹ 机器不会对这篇吐出一个看起来合理的假 r。
     ⟹ 与 PMC11656362 的 λ<1.2 行是**同一个保护机制**在不同失效模式下触发
        （那边：信噪太低；这边：两列不共享 f₁f₂）。
""")

# ---------------------------------------------------------------- §G
print("=" * 78)
print("§G  Tier-0 天花板 —— 本篇真正可用的东西")
print("=" * 78)
C_INF = 0.1827262477
print(f"  天花板（μ*=0）：Tc <= C_inf * sqrt(lam) * ω̄₂ ,  C_inf = {C_INF:.10f}")
print("  ⟹ r_true := Tc / (C_inf sqrt(lam) ω̄₂)  必须 <= 1")
print("  方向（取最大可能 r_true）：Tc 取印出上端+0.5，λ、ω̄₂ 取印出下端")
print("  ★ μ*>0 只会把 Tc 压低 ⟹ 用 μ*=0 的天花板检验 μ*>0 的数据是**保守**的（单边，合法）\n")
print(f"  {'material':22s} {'column':>12s} {'Tc[K]':>7s} {'r_true_up':>10s}")
rt = []
for (nm, tmc, tad, tme, _d, _e, _f, lam, wl, w2) in ROWS:
    for lbl, col in (("Eliashberg", tme), ("AllenDynes", tad)):
        tc = col[0] + 0.5
        ceil = C_INF * math.sqrt(lam - 0.005) * (w2 - 0.5)
        v = tc / ceil
        rt.append((v, nm, lbl))
        print(f"  {nm:22s} {lbl:>12s} {col[0]:7d} {v:10.4f}")
        chk(f"天花板 {nm}/{lbl}", v <= 1.0, f"r_true_up={v:.4f}")
rt.sort(reverse=True)
print(f"\n  最紧一点：{rt[0][1]} / {rt[0][2]}  r_true_up = {rt[0][0]:.4f}")
print(f"  语料现任最紧：Nb3Ge 0.6260  ⟹ {'夺冠' if rt[0][0] > 0.6260 else '未夺冠'}")
print("""
  ⚠⚠ 两条限制，必须与数字一起报：
  1. **Tc 不是实验值**。Tc(AD) 是公式值 ⟹ 与天花板同源（都是 λ、ω̄₂ 的函数）⟹ 自由度低。
     Tc(Eliashberg) 是**真解 Eliashberg 方程**得到的 ⟹ 它才是我们定理界的那个对象，
     检验它是**有内容的**（不是恒等式）；但它仍是计算值、非实验测量
     ⟹ **两列都不得计入 0/41 实验语料计数**（同 YThH₁₈ 的处置）。
  2. Zr₂H₁₇ 行在 §C/§E 两处诊断为不自洽 ⟹ 它的任何判决都带保留。
""")

# ---------------------------------------------------------------- §H
print("=" * 78)
print("§H  投毒阴性对照（检测器必须走通检测路径）")
print("=" * 78)


def run_gates(rows, mus=MUS):
    """返回 (Mc 重叠数, AD 重叠数, 反解可解数, 最优变体, 最优 max|rel|%, 最优变体逐点重叠数)。"""
    nmc = nad = nsol = 0
    for (nm, tmc, tad, _c, _d, _e, _f, lam, wl, w2) in rows:
        for j, ms in enumerate(mus):
            lo, hi = interval(tc_base, lam, wl, ms)
            if not (hi < tmc[j] - 0.5 or lo > tmc[j] + 0.5):
                nmc += 1
            vals = []
            for l_ in (lam - 0.005, lam + 0.005):
                for w_ in (wl - 0.5, wl + 0.5):
                    for x_ in (w2 - 0.5, w2 + 0.5):
                        vals.append(tc_ad_correct(l_, ms, w_, x_))
            if not (max(vals) < tad[j] - 0.5 or min(vals) > tad[j] + 0.5):
                nad += 1
            c0 = 1.82 * (1.0 + 6.3 * ms)
            rstar = 1.0 + math.sqrt(1.0 + lam ** 2 / c0 ** 2)
            if (tad[j] / tmc[j]) / f1_ad(lam, ms) <= f2_ad(lam, ms, rstar):
                nsol += 1
    sc = []
    for key, f2v in V.items():
        e = []
        for (nm, _a, tad, _c, _d, _e2, _f, lam, wl, w2) in rows:
            r = w2 / wl
            for j, ms in enumerate(mus):
                c = tc_base(lam, ms, wl) * f1_ad(lam, ms) * f2v(lam, ms, r)
                e.append(abs(tad[j] - c) / c * 100.0)
        sc.append((max(e), key))
    sc.sort()
    bestkey = sc[0][1]
    # ★ 正经判据：赢家变体的**逐点重叠数**（传播输入印刷精度），不是 max|rel| 的倍数
    f2b = V[bestkey]

    def tc_var(l_, ms, w_, x_):
        return tc_base(l_, ms, w_) * f1_ad(l_, ms) * f2b(l_, ms, x_ / w_)

    nov = 0
    for (nm, _a, tad, _c, _d, _e2, _f, lam, wl, w2) in rows:
        for j, ms in enumerate(mus):
            vals = []
            for l_ in (lam - 0.005, lam + 0.005):
                for w_ in (wl - 0.5, wl + 0.5):
                    for x_ in (w2 - 0.5, w2 + 0.5):
                        vals.append(tc_var(l_, ms, w_, x_))
            if not (max(vals) < tad[j] - 0.5 or min(vals) > tad[j] + 0.5):
                nov += 1
    return nmc, nad, nsol, bestkey, sc[0][0], nov


base = run_gates(ROWS)
print(f"  基线：Mc 重叠 {base[0]}/8 · AD(正确式) 重叠 {base[1]}/8 · 反解可解 {base[2]}/8")
print(f"        最优变体 = {base[3]}  (max|rel| {base[4]:.2f}%)，其**逐点重叠 {base[5]}/8**")
print("""
  ★ 报警判据用「赢家变体的逐点重叠数 / 赢家身份 / Mc 门」三者的**恶化**，
    **不用** max|rel| 的倍数——倍数判据对「把两行的 ω̄₂ 互换」这类投毒不灵敏
    （实测：1.56%→1.94% 不到 2 倍，但逐点重叠会掉）。
    这是我第一版判据的缺陷，已订正；订正本身也要有方向：逐点重叠只在
    「算出区间与印出区间脱开」时下降 ⟹ 它不会因为投毒变**宽松**。
""")

poisons = []
# P1: ω̄₂ 按 THz 读（本篇印的是 K）——单位投毒
p = [list(r) for r in ROWS]
for r_ in p:
    r_[9] = int(round(r_[9] / K_PER_THZ))
poisons.append(("P1 ω̄₂ 误当 THz 读", [tuple(x) for x in p]))
# P2: Tc^AD 列整体 +5%
p = [list(r) for r in ROWS]
for r_ in p:
    r_[2] = (int(round(r_[2][0] * 1.05)), int(round(r_[2][1] * 1.05)))
poisons.append(("P2 Tc^AD 列 +5%", [tuple(x) for x in p]))
# P3: 交换两行的 ω̄₂
p = [list(r) for r in ROWS]
p[0][9], p[1][9] = p[1][9], p[0][9]
poisons.append(("P3 交换前两行 ω̄₂", [tuple(x) for x in p]))
# P4: μ* 端点反序（0.15 配 hi）
poisons.append(("P4 μ* 端点反序", ROWS))
# P5: λ 全部 ×1.1
p = [list(r) for r in ROWS]
for r_ in p:
    r_[7] = round(r_[7] * 1.1, 2)
poisons.append(("P5 λ 全部 ×1.1", [tuple(x) for x in p]))

print(f"\n  {'投毒':24s} {'Mc':>4s} {'AD':>4s} {'解':>4s} {'重叠':>5s} {'最优变体':>28s} {'max|rel|%':>10s} 报警")
nalarm = 0
for (lbl, rows) in poisons:
    mm = (MUS[1], MUS[0]) if lbl.startswith("P4") else MUS
    g = run_gates(rows, mm)
    alarm = (g[0] < base[0]) or (g[3] != base[3]) or (g[5] < base[5])
    if alarm:
        nalarm += 1
    print(f"  {lbl:24s} {g[0]:4d} {g[1]:4d} {g[2]:4d} {g[5]:5d} {g[3]:>28s} {g[4]:10.2f} {'Y' if alarm else 'N'}")
chk("投毒全部报警", nalarm == len(poisons), f"{nalarm}/{len(poisons)}")
print("""
  ⟹ 门有分辨力：§C 的 5/8 与 §E 的变体判决都**不是恒真式**。
     P4（μ* 端点反序）是最容易犯、最难自查的转录错，门抓得到。
     P3（互换两行 ω̄₂）只靠逐点重叠抓得到，靠 max|rel| 倍数抓不到
     ⟹ 顺带量出一件事：**Tc^AD 列对 ω̄₂ 的敏感度本来就弱**
        （dln f₂/dln r ≈ 0.49），这正是反解精度只有 ±0.02 量级的根源。
""")

# ---------------------------------------------------------------- §I
print("=" * 78)
print("§I  自由度审计")
print("=" * 78)
print("""
  印出的数 ≠ 独立检验数。逐项数：

  | 断言 | 印出的量 | 有效自由度 | 说明 |
  |---|---|---|---|
  | Eq.(10) 自洽 8/8 | μ*_E, ω_ph, ω_max | **8**（两端点 × 4 材料） | 恒等式，只测转录；不是物理证据 |
  | Jensen r≥1 4/4 | ω̄₂, ω_log | **0** | 单位已印出 ⟹ 无需推断 ⟹ 该检验在本篇**空洞** |
  | (M1a) Mc 列 5/8 | λ, ω_log, μ* | **4 行 × 1**（两端点对同一乘性偏差高度相关） | 1 行结构性失效 + 1 点边缘 |
  | (M1b) AD 列 0/8 | + ω̄₂ | **4 行 × 1** 同上 | 判决一致 ⟹ 是系统性、非噪声 |
  | (b) 变体唯一 | 8 点 | **8 点 / 0 参数** | 每个候选都是零参数假设 ⟹ 比 max\\|rel\\| 合法 |
  | 天花板 8/8 | Tc, λ, ω̄₂ | AD 列 **≈0**（与天花板同源）· ME 列 **4** | ⟹ 只 ME 列有内容 |

  ★ 最重要的一条：**Tc(AD) 那 4 点对天花板的自由度近 0**
    ——AD 公式本身就是 λ、ω̄₂、ω_log、μ* 的函数，而天花板也是 λ、ω̄₂ 的函数
    ⟹ 「Tc(AD) 不超天花板」主要由 AD 公式的大 λ 渐近常数 0.1878 vs C_inf 0.1827 决定，
       不由材料信息决定（同 hydride CSV 的填充因子循环量事故）。
""")

print("=" * 78)
if FAIL:
    print(f"RESULT: {len(FAIL)} CHECK(S) FAILED")
    for f in FAIL:
        print("   -", f)
    sys.exit(1)
print("RESULT: ALL CHECKS PASSED")
print("=" * 78)
