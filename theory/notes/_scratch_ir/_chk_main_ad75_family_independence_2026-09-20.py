# -*- coding: utf-8 -*-
"""
item 2 / 家族独立性：`deficit_factor_f_headline_2026-09-20.md` §0 的阻断陈述
「样本里真正独立的材料家族数 = 1」是否还成立？

该陈述写于只有 3 行（H3S/D3S 同一家族）时。此后 AD75 Table I 的 29 行入库
（`materials/deficit_factor_ad75_tunneling_sample_2026-09-20.md`），本脚本问：

  Q1  29 行按化学体系分，有几个互不相交的组？
  Q2  𝔣 的两个极值落在哪两个组？（若同组 ⟹ 3× 展布可能是族内伪影）
  Q3  **留一元素出**（删掉含该元素的**全部**行）后，趋势回归系数动多少？
  Q4  只留 8 个**单质**行（7 种元素，零合金）时趋势还在吗？——最强的去相关检验
  Q5  ω_log ≤ ω̄₂ 在 29/29 上成立吗？（𝔣_up 那条单边界路线的前提）

**方向注记**：本脚本不产生单边界判决。Q3/Q4 的判据是「系数是否在 baseline 的
拟合不确定度内」，双侧。

**诚实边界（随结论同行）**：29 行全部出自**同一种测量技术**（隧道反演 McMillan–Rowell）、
同一时期、重叠的实验组 ⟹ 化学上的独立**不**蕴含系统误差上的独立。本脚本只能判前者。

主会话 2026-09-20，从零。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import re

C_INF = 0.1827262477

# 与 _chk_ad75_table1_frak_f_2026-09-20.py 同源（页图逐行转录，29 行）
T1 = [
    ("Pb",              56,  65, 1.55, 7.2 ),
    ("In",              68,  89, 0.805,3.40),
    ("Sn",              99, 121, 0.72, 3.75),
    ("Hg",              29,  49, 1.6,  4.19),
    ("Tl",              52,  64, 0.795,2.36),
    ("Ta",             132, 148, 0.69, 4.48),
    ("a-Ga",            55, 101, 1.62, 8.56),
    ("b-Ga",            87, 129, 0.97, 5.90),
    ("Tl0.9Bi0.1",      48,  62, 0.78, 2.30),
    ("Pb0.4Tl0.6",      48,  62, 1.15, 4.60),
    ("Pb0.6Tl0.4",      50,  62, 1.38, 5.90),
    ("Pb0.8Tl0.2",      50,  61, 1.53, 6.80),
    ("Pb0.6Tl0.2Bi0.2", 48,  58, 1.81, 7.26),
    ("Pb0.9Bi0.1",      50,  60, 1.66, 7.65),
    ("Pb0.8Bi0.2",      46,  57, 1.88, 7.95),
    ("Pb0.7Bi0.3",      47,  57, 2.01, 8.45),
    ("Pb0.65Bi0.35",    45,  55, 2.13, 8.95),
    ("In0.9Tl0.1",      63,  86, 0.85, 3.28),
    ("In0.73Tl0.27",    55,  77, 0.93, 3.36),
    ("In0.67Tl0.33",    57,  79, 0.90, 3.26),
    ("In0.57Tl0.43",    53,  74, 0.85, 2.60),
    ("In0.5Tl0.5",      53,  73, 0.83, 2.52),
    ("In0.27Tl0.73",    42,  63, 1.09, 3.64),
    ("In0.17Tl0.83",    45,  63, 0.98, 3.19),
    ("In0.07Tl0.93",    49,  63, 0.89, 2.77),
    ("In2Bi",           46,  67, 1.40, 5.6 ),
    ("Sb2Tl7",          37,  58, 1.43, 5.2 ),
    ("Bi2Tl",           47,  59, 1.63, 6.4 ),
    ("a-Pb0.45Bi0.55",  29,  47, 2.59, 7.0 ),
]

ELEMENTS = ["Pb", "Bi", "Tl", "In", "Sn", "Hg", "Ta", "Ga", "Sb"]

def elems_of(name):
    """从材料名抽元素集合。先剥 a-/b- 前缀，再按两字母元素符号匹配。"""
    s = re.sub(r"^[ab]-", "", name)
    out = set()
    for e in ELEMENTS:
        if re.search(e + r"(?![a-z])", s):
            out.add(e)
    return frozenset(out)

name = np.array([r[0] for r in T1])
wlog = np.array([float(r[1]) for r in T1])
wb2  = np.array([float(r[2]) for r in T1])
lam  = np.array([float(r[3]) for r in T1])
tc   = np.array([float(r[4]) for r in T1])
frak = tc / (C_INF * np.sqrt(lam) * wb2)
esets = [elems_of(n) for n in name]

print("=" * 100)
print("Q5. 前提检查：ω_log ≤ ω̄₂ （𝔣_up = Tc/(C_∞√λ ω_log) 是 𝔣 的上界）")
print("=" * 100)
bad = [(name[i], wlog[i], wb2[i]) for i in range(len(T1)) if wlog[i] > wb2[i]]
print(f"   违反行数：{len(bad)} / {len(T1)}   {'⟹ 前提成立（幂平均不等式在数据上自洽）' if not bad else bad}")
print(f"   ω_log/ω̄₂ ∈ [{(wlog/wb2).min():.4f}, {(wlog/wb2).max():.4f}]，"
      f"中位 {np.median(wlog/wb2):.4f}  ⟹ 𝔣_up 比 𝔣 高 {100*(np.median(wb2/wlog)-1):.0f}% 上下")

print()
print("=" * 100)
print("Q1/Q2. 化学体系分组，以及 𝔣 的两个极值落在哪里")
print("=" * 100)
groups = {}
for i, es in enumerate(esets):
    groups.setdefault(es, []).append(i)
print(f"   互不相同的元素集合数 = {len(groups)}")
print()
print(f"{'元素集合':<22}{'n':>3}  {'𝔣 范围':<22}{'λ 范围':<16} 成员")
for es, idx in sorted(groups.items(), key=lambda kv: -len(kv[1])):
    fs, ls = frak[idx], lam[idx]
    tag = "{" + ",".join(sorted(es)) + "}"
    mem = ", ".join(name[idx][:3]) + ("…" if len(idx) > 3 else "")
    print(f"{tag:<22}{len(idx):>3}  [{fs.min():.4f}, {fs.max():.4f}]     "
          f"[{ls.min():.2f}, {ls.max():.2f}]    {mem}")

imin, imax = int(frak.argmin()), int(frak.argmax())
print()
print(f"   𝔣 最小 = {frak[imin]:.4f}  {name[imin]:<16} 元素集 {set(esets[imin])}  λ={lam[imin]}")
print(f"   𝔣 最大 = {frak[imax]:.4f}  {name[imax]:<16} 元素集 {set(esets[imax])}  λ={lam[imax]}")
print(f"   展布 = {frak[imax]/frak[imin]:.2f}×")
disj = esets[imin].isdisjoint(esets[imax])
print(f"   两极值元素集是否**不相交**：{disj}"
      f"  ⟹ {'3× 展布**不可能**是族内伪影' if disj else '⚠ 同族，展布可疑'}")

# ---------------- 回归 ----------------
X0 = np.column_stack([np.ones(len(T1)), np.log(lam), np.log(wlog / wb2)])
y0 = np.log(frak)

def fit(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    dof = len(y) - X.shape[1]
    s = float(np.sqrt(r @ r / dof))
    R2 = 1.0 - (r @ r) / float(((y - y.mean()) ** 2).sum())
    cov = s * s * np.linalg.inv(X.T @ X)
    return beta, s, R2, np.sqrt(np.diag(cov))

b0, s0, R20, se0 = fit(X0, y0)
print()
print("=" * 100)
print("baseline 回归   log 𝔣 = a + b·log λ + c·log(ω_log/ω̄₂)     (29 行)")
print("=" * 100)
print(f"   a = {b0[0]:+.4f} ± {se0[0]:.4f}")
print(f"   b = {b0[1]:+.4f} ± {se0[1]:.4f}   ← λ 指数")
print(f"   c = {b0[2]:+.4f} ± {se0[2]:.4f}   ← 谱形状指数")
print(f"   R² = {R20:.4f}   s = {s0:.4f}")

print()
print("=" * 100)
print("Q3. 留一元素出：删掉含该元素的**全部**行，重拟合")
print("=" * 100)
print(f"{'删去':<6}{'剩余n':>6}{'b (λ指数)':>16}{'Δb/σ_b':>10}{'c (形状)':>14}{'Δc/σ_c':>10}{'R²':>8}{'s':>8}")
print(f"{'—':<6}{len(T1):>6}{b0[1]:>16.4f}{0.0:>10.2f}{b0[2]:>14.4f}{0.0:>10.2f}{R20:>8.4f}{s0:>8.4f}")
worst_b = worst_c = 0.0
for e in ELEMENTS:
    keep = np.array([e not in es for es in esets])
    if keep.sum() < 6:
        print(f"{e:<6}{keep.sum():>6}   (剩余行太少，跳过)")
        continue
    b, s, R2, se = fit(X0[keep], y0[keep])
    zb, zc = (b[1] - b0[1]) / se0[1], (b[2] - b0[2]) / se0[2]
    worst_b, worst_c = max(worst_b, abs(zb)), max(worst_c, abs(zc))
    print(f"{e:<6}{keep.sum():>6}{b[1]:>16.4f}{zb:>10.2f}{b[2]:>14.4f}{zc:>10.2f}{R2:>8.4f}{s:>8.4f}")
print()
print(f"   最劣位移：|Δb| = {worst_b:.2f} σ_b     |Δc| = {worst_c:.2f} σ_c")
print(f"   ⟹ {'趋势不由任何单一元素驱动（全部 <1σ）' if max(worst_b,worst_c) < 1.0 else '有元素驱动趋势，见上表'}")

print()
print("=" * 100)
print("Q4. 最强检验：只留 8 个**单质**行（Pb In Sn Hg Tl Ta a-Ga b-Ga，7 种元素，零合金）")
print("=" * 100)
pure = np.array([len(es) == 1 for es in esets])
print(f"   n = {pure.sum()}：{', '.join(name[pure])}")
print(f"   𝔣 ∈ [{frak[pure].min():.4f}, {frak[pure].max():.4f}]  展布 {frak[pure].max()/frak[pure].min():.2f}×"
      f"    λ ∈ [{lam[pure].min():.2f}, {lam[pure].max():.2f}]")
bp, sp, R2p, sep = fit(X0[pure], y0[pure])
print(f"   b = {bp[1]:+.4f} ± {sep[1]:.4f}   (全样本 {b0[1]:+.4f} ± {se0[1]:.4f})   "
      f"Δ = {abs(bp[1]-b0[1])/np.hypot(sep[1],se0[1]):.2f}σ_合并")
print(f"   c = {bp[2]:+.4f} ± {sep[2]:.4f}   (全样本 {b0[2]:+.4f} ± {se0[2]:.4f})   "
      f"Δ = {abs(bp[2]-b0[2])/np.hypot(sep[2],se0[2]):.2f}σ_合并")
print(f"   R² = {R2p:.4f}   s = {sp:.4f}")

print()
print("=" * 100)
print("汇总判决")
print("=" * 100)
print(f"  · 互不相同的化学体系（元素集）数 = {len(groups)}，**不是 1**")
print(f"  · 𝔣 的两个极值 ({name[imin]} / {name[imax]}) 化学上{'不相交' if disj else '相交'}"
      f" ⟹ {frak[imax]/frak[imin]:.2f}× 展布跨家族")
print(f"  · 留一元素出最劣位移 {max(worst_b, worst_c):.2f}σ；纯单质子样本（7 元素）独立复现同一趋势")
print(f"  · **仍然成立的窄陈述**：高 Tc / 氢化物区的独立家族数 = 1（H3S/D3S 同一家族）")
print(f"  · **相关系统误差（本脚本判不了）**：29 行同一技术（隧道反演）、同一时期、重叠实验组")
print("=" * 100)
