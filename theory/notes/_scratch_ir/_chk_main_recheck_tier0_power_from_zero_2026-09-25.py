# -*- coding: utf-8 -*-
"""主会话从零复核 `tier0-power` 的 Q1/Q2，并补它未做的「四档 μ* × 氢化物端」。

写这个脚本的三个理由（都是复核发现，不是重跑）：
 1. 它的 Q1 自查门 `Tc_AD > 天花板` **是恒等式**：batch-q 的 Tc_AD 列 194 行空、55 行恰为 0.0，
    `0.0 > ceiling` 永假 ⟹ 自由度 0，门永不可能失败。且 194/249 = 78% 的行**从未进门**。
    （这条缺陷继承自主会话自己的 §4 Tier 0 表，不是交付方引入的。）
 2. 它的 Q1 头条「Tc≥77K 档代价恰好 0」是在**四档全部剪到 100%** 的区间量的 ⟹ 该区间四档不可分辨，
    自由度 0。要判「0.10 承不承重」，必须找到**四档给出不同答案**的目标 Tc 区间。
 3. 它的 Q2 对 YH6/YH9 用了 `C(λ,0.13)`，而这两个材料的 `mustar_used = 0.11`。
    C 对 μ* 单调降 ⟹ 用 0.13>0.11 给出的是**偏低**的天花板 = 上界的下估。
    方向后果分两半：假阴性门（PASS 型）**仍然成立且更强**；
    但**剪枝计数与 𝔣 被放大** ⟹「分辨力」被高估。必须用 μ*≤真值 的档重算。

全部结论 `numerical-evidence`（float64 + 有限数据集）。
常数来源与交付脚本同（`../ceiling/C_lambda_certified_curve_2026-09-24.md` §1 的 C_upper 列
+ `_chk_main_verify_certC_mu010_from_zero_2026-09-24.py` 的 μ*=0.10 六格），
但本文件**独立重打**这些表并与交付脚本逐格对撞，不 import 它。
"""
import os
import sys
import csv
import math

os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

KB = 0.086173332621452  # meV/K
C_INF = 0.1827262477

C_UP = {
    0.00: [(2.0, 0.151236028), (3.2953, 0.164718304), (4.0, 0.168304891), (6.0, 0.173820789)],
    0.10: [(1.5, 0.108973543), (2.0, 0.123102031), (2.5, 0.131332069), (3.0, 0.136767973),
           (3.2953, 0.139202029), (4.0, 0.143488274)],
    0.13: [(1.0, 0.085109595), (1.5, 0.101992146), (2.0, 0.116687205), (2.5, 0.125402058),
           (3.0, 0.130991149), (3.2953, 0.133361200), (4.0, 0.137762123), (5.0, 0.141888517),
           (6.0, 0.144677171), (8.0, 0.148030265), (10.0, 0.150564456)],
}
D = os.path.dirname(os.path.abspath(__file__))


def c_up(lam, key):
    """声音上界：取 λ 的上邻认证格（C 对 λ 单调升 ⟹ 用更大 λ 的格是上界）。表右端外回退 C_∞。"""
    for L, c in C_UP[key]:
        if lam <= L:
            return c
    return C_INF


def median(v):
    s = sorted(v)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


print("=" * 100)
print("A. Q1 自查门的自由度审计（我的头号发现）")
print("=" * 100)
p1 = os.path.join(D, "..", "..", "data", "external", "cerqueira6912_batchq_sample.csv")
rows = list(csv.DictReader(open(p1, encoding="utf-8")))
tcad_raw = [r["Tc_AD_mustar010_K"] for r in rows]
n_empty = sum(1 for t in tcad_raw if t == "")
vals = sorted({t for t in tcad_raw if t != ""})
print(f"batch-q 总行数            : {len(rows)}")
print(f"Tc_AD 列为空的行          : {n_empty}  ({100.0*n_empty/len(rows):.1f}%)  ⟹ **从未进过自查门**")
print(f"Tc_AD 列的全部相异非空取值: {vals}")
print(f"Tc_Eliashberg 列非空行数  : {sum(1 for r in rows if r['Tc_Eliashberg_mustar010_K'])}")
print("⟹ 判词：`Tc_AD > 天花板` 在本数据集上是 `0.0 > 正数`，**恒假**。")
print("   这不是「门通过」，是**门没有自由度**。0/55 不构成任何证据。")

lam = [float(r["lambda"]) for r in rows]
w2 = [float(r["w2av_meV"]) for r in rows]
M2 = [l * w * w for l, w in zip(lam, w2)]
print(f"\nλ ∈ [{min(lam):.4f}, {max(lam):.4f}]   ⟹ 全部落在 μ*=0.13 表最左格(λ=1.0)**之左**，")
print("   四档 c_up() 于是对每一行都返回同一个常数 ⟹ 本数据集上『分辨 λ』这件事自由度也是 0。")

print("\n" + "=" * 100)
print("B. Q1 重做：找出四档**真正可分辨**的目标 Tc 区间")
print("=" * 100)
takes = {"C_inf": [C_INF] * len(lam)}
for k in (0.00, 0.10, 0.13):
    takes[f"C(l,{k:.2f})"] = [c_up(l, k) for l in lam]
ceil = {k: [c * math.sqrt(m) / KB for c, m in zip(v, M2)] for k, v in takes.items()}
for k, v in ceil.items():
    print(f"  {k:<11} 常数 C={takes[k][0]:.9f}   天花板范围 [{min(v):.2f}, {max(v):.2f}] K")

print(f"\n{'目标Tc':>7}" + "".join(f"{k:>12}" for k in ceil) + "   四档是否可分辨")
prev_all_same = None
disc = []
for tgt in [1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 77, 100, 200, 300]:
    rates = [100.0 * sum(1 for c in ceil[k] if c < tgt) / len(lam) for k in ceil]
    same = (max(rates) - min(rates)) < 1e-12
    if not same:
        disc.append(tgt)
    print(f"{tgt:>6}K" + "".join(f"{r:>11.1f}%" for r in rates) + f"   {'—— 不可分辨（自由度 0）' if same else '可分辨'}")
print(f"\n⟹ 四档给出不同答案的目标 Tc = {disc} K")
print("   Tc≥77K 区间四档**全部 100%**，交付方正是在这里读出「代价 0pp」——")
print("   那不是『0.10 不承重』的证据，那是**该区间无自由度**。")

print("\n" + "=" * 100)
print("C. Q2 重做：四档 μ* × 真实氢化物（交付方只跑了 0.13 这一档）")
print("=" * 100)
p2 = os.path.join(D, "..", "..", "data", "hydrides_theory.csv")
hyd = list(csv.DictReader(open(p2, encoding="utf-8")))
TC_EXP = {("H3S", "H"): 203, ("D3S", "D"): 150, ("LaH10", "H"): 250, ("LaD10", "D"): 180,
          ("YH6", "H"): 224, ("YH9", "H"): 243, ("CaH6", "H"): 215}


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


mats = []
for d in hyd:
    L, W = fnum(d["lambda"]), fnum(d["omega_2_meV"])
    if L is None or W is None:
        continue
    tce = TC_EXP.get((d["material"], d["isotope"].split()[0]))
    if tce is None:
        continue
    mats.append({"mat": d["material"], "treat": d["treatment"], "P": d["pressure_GPa"],
                 "lam": L, "M2": L * W * W, "tce": tce,
                 "mustar_src": fnum(d["mustar_used"]), "tier": d["evidence_tier"]})
print(f"入选行 n={len(mats)}（有 λ 且有真实 ω̄₂ 且有实验 Tc 锚）")

print(f"\n{'material':<7}{'treat':<14}{'P':>5}{'λ':>7}{'μ*_src':>7}   " +
      "".join(f"{'ceil@'+k:>14}" for k in ("0.00", "0.10", "0.13")) + f"{'Tc_exp':>8}")
for m in mats:
    cs = {k: c_up(m["lam"], float(k)) * math.sqrt(m["M2"]) / KB for k in ("0.00", "0.10", "0.13")}
    m["ceil"] = cs
    flag = ""
    if m["mustar_src"] is not None and m["mustar_src"] < 0.13:
        flag = f"  ⚠ 0.13 档 > μ*_src={m['mustar_src']} ⟹ 天花板被**低估**"
    print(f"{m['mat']:<7}{m['treat'][:13]:<14}{m['P']:>5}{m['lam']:>7.3f}"
          f"{m['mustar_src'] if m['mustar_src'] is not None else float('nan'):>7.2f}   " +
          "".join(f"{cs[k]:>14.1f}" for k in ("0.00", "0.10", "0.13")) + f"{m['tce']:>8.0f}{flag}")

print("\n𝔣 = Tc_exp / 天花板 与假阴性门，逐档：")
print(f"{'档':>8}{'𝔣_min':>9}{'𝔣_中位':>9}{'𝔣_max':>9}{'假阴性(𝔣>1)':>13}{'首次开剪目标Tc':>16}{'325K剪掉':>10}")
for k in ("0.00", "0.10", "0.13"):
    fs = [m["tce"] / m["ceil"][k] for m in mats]
    fneg = sum(1 for f in fs if f > 1.0)
    flip = None
    for tgt in range(100, 501, 25):
        if any(m["ceil"][k] < tgt for m in mats):
            flip = tgt
            break
    k325 = sum(1 for m in mats if m["ceil"][k] < 325)
    print(f"{k:>8}{min(fs):>9.3f}{median(fs):>9.3f}{max(fs):>9.3f}{fneg:>13}{flip:>16}{k325:>8}/{len(mats)}")

print("\n⟹ 自由度审计：假阴性门 `∃tgt: tce≥tgt>ceil` ⟺ `ceil < tce` ⟺ `𝔣 > 1`。")
print("   17 个目标 Tc 的扫描**不是 17 个独立检验**，是同一个 `𝔣>1` 判据的 17 次重述。")
print(f"   真自由度 = 行数 = {len(mats)}；而这 {len(mats)} 行只来自 4 个物理材料"
      f"（{sorted({m['mat'] for m in mats})}）⟹ 独立材料数 = 4。")

print("\n" + "=" * 100)
print("D. 判词（只写我实际算出来的）")
print("=" * 100)
print("D1. 交付方 Q1 的自查门恒假、Q1 头条区间无自由度 ⟹ 『0.10 不承重』**在该证据上不成立**。")
print("    （该命题本身另有来源——§6.6 的 analytic 论证——但不是这张表给的。）")
print("D2. 四档真正可分辨的区间见 B 节；要判紧度代价只能在那里读。")
print("D3. Q2 的 0.13 档对 μ*_src=0.11 的 YH6/YH9 是**低估天花板** ⟹ 分辨力被高估；")
print("    改用 0.10 档（≤0.11，方向声音）后的剪枝与 𝔣 见 C 节表。")
print("D4. 假阴性 0 这一条在三档下都成立，且 0.13 档因天花板偏低而是**更强**的 PASS。")
