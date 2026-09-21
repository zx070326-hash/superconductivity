# -*- coding: utf-8 -*-
"""
𝔣 分母的数字化误差预算：λ 对得上 ≠ √(λ⟨ω²⟩) 对得上。

被评估的计划（`_scratch_ir/_multiband_papers_circularity_reclassification_2026-09-20.md`）：
  「数字化 Golubov Fig.1 / Sanna / Ummarino 的 α²F 曲线 ⟹ 补 ω̄₂ ⟹ 净增 2–3 行非循环 𝔣 样本」，
  其隐含验收方案是「若重算出的 λ 复现论文印的 0.87，则数字化可信」。

本脚本查四件事：
  A. λ⟨ω²⟩ = 2∫α²F·ω dω 是否为恒等式（⟹ 分母只需一个积分，λ 被约掉）。
     ⚠ 该列在谐波/非谐间一致到 1e-6 **不是**非谐不变性的实证确认 ——
       `materials/deficit_budget_h3s_2026-09-17.md` §3 已拦下该假结论：
       SSCHA-on-DFPT 复用同一套电子-声子矩阵元，I₁=η/2M 在同一次计算内部形式上恒成立。
  B. 低频截断下 λ 与 √I₁ 各丢多少（⟹ I₁ 对数字化最不可靠的那一段免疫）。
  C. 分频段的权重分解（⟹ λ 与 I₁ 的权重分居谱的两端）。
  D. ★ 定向构造：能否造出「λ 精确不变、而 √I₁ 偏 1–4%」的数字化误差？
     按 CLAUDE.md「先写下使该命题恰好取等的那组参数」——δ′ 由 λ 不变条件**解析**确定，
     不是随机扫描（随机扫描打不到这个等号面）。

严格性：A 为 analytic（恒等式）+ 数值核对；B/C 为 numerical-evidence；D 为构造性。
配套笔记：`../materials/f_denominator_digitization_budget_2026-09-20.md`
"""
import os, sys
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

BASE = os.path.join(r"D:\physics\superconductivity", "theory", "data", "external",
                    "mishra_h3s_d3s", "Data-H3S-D3S-PH-EPW-Anharmonic")
tr = np.trapezoid if hasattr(np, "trapezoid") else np.trapz


def load(p):
    """EPW .a2f：第一列 ω(meV)、第二列 α²F。尾部有非数值行，逐行 try。"""
    rows = []
    for line in open(p, errors="replace"):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        try:
            v = [float(x) for x in s.split()]
        except ValueError:
            continue
        if len(v) >= 2:
            rows.append(v[:2])
    a = np.array(rows)
    return a[:, 0], a[:, 1]


def clean(w, a):
    m = w > 1e-9
    return w[m], np.clip(a[m], 0.0, None)


CASES = [("H3S-200 非谐", "H3S-200GPa/EPW_anharmonic/H3S.a2f"),
         ("H3S-200 谐波", "H3S-200GPa/EPW_harmonic/H3S.a2f"),
         ("D3S-200 非谐", "D3S-200GPa/EPW_anharmonic/D3S.a2f"),
         ("D3S-200 谐波", "D3S-200GPa/EPW_harmonic/D3S.a2f")]

print("=" * 86)
print("A. λ⟨ω²⟩ = 2∫α²F·ω dω 是恒等式吗？（若是，𝔣 的分母只需一个积分，λ 被约掉）")
print("=" * 86)
for tag, rel in CASES:
    w, a = clean(*load(os.path.join(BASE, rel.replace("/", os.sep))))
    lam = 2 * tr(a / w, w)
    w2 = np.sqrt((2 / lam) * tr(a * w, w))
    lhs, rhs = lam * w2 ** 2, 2 * tr(a * w, w)       # 两条路算同一个量
    print(f"  {tag}:  λ·ω̄₂²={lhs:14.6f}   2∫α²F·ω dω={rhs:14.6f}   相对差 {abs(lhs-rhs)/rhs:.1e}")
    print(f"      λ={lam:.6f}   ω̄₂={w2:.6f} meV   ω_log={np.exp((2/lam)*tr(a/w*np.log(w),w)):.6f} meV"
          f"   √I₁={np.sqrt(rhs):.6f} meV")
print("  ⚠ 第 3 列在谐波/非谐间一致 ≠ 非谐不变性的实证确认（见 docstring 与笔记 §1 的告警）。")

w, a = clean(*load(os.path.join(BASE, "H3S-200GPa", "EPW_anharmonic", "H3S.a2f")))
lam0, I10 = 2 * tr(a / w, w), 2 * tr(a * w, w)

print()
print("=" * 86)
print("B. 低频截断：λ 与 √I₁ 各丢多少？（H3S-200 非谐）")
print("=" * 86)
for w0 in (10, 20, 30, 40):
    k = w >= w0
    lam, I1 = 2 * tr(a[k] / w[k], w[k]), 2 * tr(a[k] * w[k], w[k])
    print(f"   截去 ω<{w0:2d} meV:  λ 丢 {(1-lam/lam0)*100:6.3f}%   "
          f"I₁ 丢 {(1-I1/I10)*100:6.3f}%   √I₁ 丢 {(1-np.sqrt(I1/I10))*100:6.3f}%")

print()
print("=" * 86)
print("C. 分频段权重：λ 与 I₁ 看的是谱的哪一段？")
print("=" * 86)
print("    频段(meV)     占 λ      占 I₁     [+10% 振幅 ⟹ λ 变 / √I₁ 变]")
for lo, hi in [(0, 40), (40, 80), (80, 120), (120, 160), (160, 200), (200, 270)]:
    k = (w >= lo) & (w < hi)
    if k.sum() < 2:
        continue
    b = a.copy(); b[k] *= 1.10
    lam, I1 = 2 * tr(b / w, w), 2 * tr(b * w, w)
    fl = 2 * tr(np.where(k, a, 0) / w, w) / lam0
    fi = 2 * tr(np.where(k, a, 0) * w, w) / I10
    print(f"    [{lo:3d},{hi:3d})  {fl*100:6.2f}%   {fi*100:6.2f}%      "
          f"{(lam/lam0-1)*100:+6.3f}%  /  {(np.sqrt(I1/I10)-1)*100:+6.3f}%")

print()
print("=" * 86)
print("D. ★ 定向构造：造一族使 λ **精确不变** 的数字化误差，看 √I₁ 偏多少")
print("=" * 86)
kl, kh = (w < 60), (w >= 140)
Al = 2 * tr(np.where(kl, a, 0) / w, w)      # 低频段对 λ 的贡献
Ah = 2 * tr(np.where(kh, a, 0) / w, w)      # 高频段对 λ 的贡献
for d in (0.05, 0.10, 0.20):
    dp = d * Al / Ah                        # 解析地选 δ′ 使 λ 净变化 = 0
    b = a.copy(); b[kl] *= (1 + d); b[kh] *= (1 - dp)
    lam, I1 = 2 * tr(b / w, w), 2 * tr(b * w, w)
    print(f"   低频<60 +{d*100:4.0f}% / 高频≥140 −{dp*100:5.2f}%  ⟹  "
          f"λ 变 {(lam/lam0-1)*100:+.4f}%   但 √I₁ 变 {(np.sqrt(I1/I10)-1)*100:+.3f}%")
print("   ⟹ 『数字化后复现论文印出的 λ』**不是** I₁ 的有效负对照。")
print("   （方向注记：δ<0 给同样大小的正偏，该族双向可走，非单边界。）")
print()
print("   消费点：𝔣 := Tc_exp/(C_∞·√I₁) ⟹ √I₁ 的相对误差 1:1 传进 𝔣。")
print("           ⚠ 初版在此处写『带宽 ±2% ⟹ 该路线不达标』——**已撤回**，那是 legacy")
print("             压强未匹配锚点下的值。真实带宽见下面 E 节（≈9.5%）。")
print("=" * 86)

# ---------------------------------------------------------------------------
# E. ★ 订正（2026-09-20）：本脚本初版把「𝔣 主张带宽」写成 ±2%，那是 legacy
#    （压强未匹配）锚点下的值；`deficit_factor_f_headline_2026-09-20.md` §2.3
#    已标明该紧带是**运气型抵消**，不是精度证据。这里从零重算两档带宽。
#    消费点：新行的数字化误差要低于**材料间散布**才有分辨力。
# ---------------------------------------------------------------------------
C_INF, KB = 0.1827262477, 0.08617333262          # 无量纲 / meV per K


def ceiling_K(rel):
    w, a = clean(*load(os.path.join(BASE, rel.replace("/", os.sep))))
    return C_INF * np.sqrt(2 * tr(a * w, w)) / KB


print()
print("=" * 86)
print("E. 𝔣 的真实带宽（判断 §D 的 1–4% 是大还是小）")
print("=" * 86)
H = ceiling_K("H3S-200GPa/EPW_anharmonic/H3S.a2f")
D = ceiling_K("D3S-200GPa/EPW_anharmonic/D3S.a2f")
print(f"   天花板：H3S-200anh = {H:.3f} K     D3S-200anh = {D:.3f} K")
for lab, tcH, tcD in (("legacy（压强未匹配）", 203.0, 150.0),
                      ("压强匹配（干净锚点）", 184.04, 145.5)):
    v = np.array([tcH / H, tcD / D])
    sp = (v.max() - v.min()) / v.mean()
    print(f"   {lab}: 𝔣 = {v[0]:.5f} / {v[1]:.5f}   散布 = {sp*100:.2f}%")
print("   ⟹ 干净锚点下带宽 ≈9.5%，§D 的 1–4% 是**小**的 ⟹ 数字化行有信息量。")
print("      初版据 ±2% 得出的『不达标』判决**撤回**（见笔记 §4 的两段式订正）。")
print("=" * 86)
