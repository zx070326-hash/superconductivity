# -*- coding: utf-8 -*-
"""Tier 0 / Tier 1 预算测量脚本（团队派活，2026-09-25）。

只量数，不证明任何新定理。全部结论 numerical-evidence（float64 + 有限数据集），
不写 certificate、不写「已证」。

Q1（对应 prescriptive_protocol §6.6）：μ*_min 从 0.10 回退到 0（或 C_∞）要付几个百分点的剪枝率？
Q2（对应 prescriptive_protocol §6.3）：强耦合真实氢化物上筛子的分辨力——目标 Tc 从 100→500 K 扫，
    每档剪掉几个/存活几个，筛子第一次翻转在哪，假阴性必须为 0。

方向规则（每一步代入前复述一遍）：kTc ≤ C(λ,μ*)·√M2，C 对 λ 单调升、对 μ* 单调降（仓库已认证）。
⟹ 声音上界要 λ 取上、μ* 取下。任何把 C 放大的取法只会**减少**剪枝（更保守），
   绝不会把一个真实 Tc≥目标 的材料误剪 ⟹ 这是本脚本两个自查门的理论根据。

数据/常数来源（读哪一行，写在这里，不改这些既有脚本）：
 - C_∞、C_UP(0.00)/C_UP(0.13) 四档表：逐字抄自
   `_chk_main_hydride_ceiling_budget_2026-09-24.py` 里的 C_UP 字典（该脚本转引
   `../ceiling/C_lambda_certified_curve_2026-09-24.md` §1 的 C_upper 列）。
 - C_UP(0.10) 六格：逐字抄自
   `_chk_main_verify_certC_mu010_from_zero_2026-09-24.py` 的 DEL 表第 4 列（C_upper），
   该脚本已完成 G1–G7 七门从零复核（PASS）。
 - Cerqueira batch-q：`theory/data/external/cerqueira6912_batchq_sample.csv`（n≈249，
   与 `_chk_main_prescriptive_design_budget_2026-09-24.py` 的 Q1 用同一份文件、同一列名）。
 - 氢化物：`theory/data/hydrides_theory.csv`，TC_EXP 锚表逐字抄自
   `_chk_main_hydride_ceiling_budget_2026-09-24.py`。
"""
import os
import sys
import csv
import math

os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

KB = 0.086173332621452  # meV/K
C_INF = 0.1827262477    # analytic, λ/μ*全盲

# 四档认证 C_upper 表（λ, C_upper），逐字转引（见文件头注释）
C_UP = {
    0.00: [(2.0, 0.151236028), (3.2953, 0.164718304), (4.0, 0.168304891), (6.0, 0.173820789)],
    0.10: [(1.5, 0.108973543), (2.0, 0.123102031), (2.5, 0.131332069), (3.0, 0.136767973),
           (3.2953, 0.139202029), (4.0, 0.143488274)],
    0.13: [(1.0, 0.085109595), (1.5, 0.101992146), (2.0, 0.116687205), (2.5, 0.125402058),
           (3.0, 0.130991149), (3.2953, 0.133361200), (4.0, 0.137762123), (5.0, 0.141888517),
           (6.0, 0.144677171), (8.0, 0.148030265), (10.0, 0.150564456)],
}
D = os.path.dirname(os.path.abspath(__file__))


def c_upper(lam, tbl_key):
    """声音上界：λ 取上邻认证格（C 对 λ 单调升，certificate）。超出表右端回退到 C_∞。"""
    for L, c in C_UP[tbl_key]:
        if lam <= L:
            return c
    return C_INF


print("=" * 112)
print("Q1 — μ*_min 回退代价：Cerqueira batch-q 上四种 C 取法 × 五个目标 Tc 的剪枝率")
print("=" * 112)

p1 = os.path.join(D, "..", "..", "data", "external", "cerqueira6912_batchq_sample.csv")
rows1 = [r for r in csv.DictReader(open(p1, encoding="utf-8")) if r["lambda"] and r["w2av_meV"]]
lam1 = np.array([float(r["lambda"]) for r in rows1])
w2_1 = np.array([float(r["w2av_meV"]) for r in rows1])
M2_1 = lam1 * w2_1 * w2_1
tc_ad1 = np.array([float(r["Tc_AD_mustar010_K"]) if r["Tc_AD_mustar010_K"] else np.nan for r in rows1])
print(f"n={len(rows1)}   λ∈[{lam1.min():.3f},{lam1.max():.3f}]   "
      f"M2∈[{M2_1.min():.0f},{M2_1.max():.0f}] meV²")

takes = {
    "C_inf":     np.full_like(lam1, C_INF),
    "C(l,0.00)": np.array([c_upper(l, 0.00) for l in lam1]),
    "C(l,0.10)": np.array([c_upper(l, 0.10) for l in lam1]),
    "C(l,0.13)": np.array([c_upper(l, 0.13) for l in lam1]),
}
ceilings = {k: v * np.sqrt(M2_1) / KB for k, v in takes.items()}

TARGETS_Q1 = (10, 30, 77, 200, 300)
q1_bad_total = 0
q1_table = {}
for name, ceil in ceilings.items():
    finite_mask = ~np.isnan(tc_ad1)
    bad = int((tc_ad1[finite_mask] > ceil[finite_mask]).sum())
    q1_bad_total += bad
    row = []
    for tgt in TARGETS_Q1:
        kill = ceil < tgt
        row.append((int(kill.sum()), kill.mean() * 100.0, int((~kill).sum())))
    q1_table[name] = (row, bad, int(finite_mask.sum()))

hdr = f"{'C 取法':<12}" + "".join(f"{'Tc='+str(t)+'K':>18}" for t in TARGETS_Q1) + f"{'自查(Tc_AD>ceil)':>20}"
print(hdr)
for name, (row, bad, ntot) in q1_table.items():
    cells = "".join(f"{k:>6}剪/{s:<4}活({r:4.1f}%)" for k, r, s in row)
    print(f"{name:<12}{cells}{f'{bad}/{ntot}':>20}")

print(f"\n自查门 Q1：所有取法下 Tc_AD > 天花板 的总行数 = {q1_bad_total}（应为 0；非 0 立即报警）")
if q1_bad_total != 0:
    print("  ⚠⚠⚠ 天花板被证伪！停止后续解读。")

# 剪枝率代价：以 0.10 为基准，看退到 0.00 / C_inf 各目标 Tc 下的剪枝率差
print("\n剪枝率代价表（相对 μ*=0.10 基准，单位：百分点，负号=剪枝率下降=更保守）：")
base = {t: q1_table["C(l,0.10)"][0][i][1] for i, t in enumerate(TARGETS_Q1)}
for name in ("C(l,0.00)", "C_inf"):
    diffs = [q1_table[name][0][i][1] - base[t] for i, t in enumerate(TARGETS_Q1)]
    print(f"  0.10 → {name:<10}: " + "  ".join(f"Tc={t}K: {d:+6.2f}pp" for t, d in zip(TARGETS_Q1, diffs)))

print("\n" + "=" * 112)
print("Q2 — 强耦合真实氢化物上的分辨力：目标 Tc 100→500 K 扫描（步长 25 K）")
print("=" * 112)

p2 = os.path.join(D, "..", "..", "data", "hydrides_theory.csv")
hyd = list(csv.DictReader(open(p2, encoding="utf-8")))


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


TC_EXP = {("H3S", "H"): 203, ("D3S", "D"): 150, ("LaH10", "H"): 250, ("LaD10", "D"): 180,
          ("YH6", "H"): 224, ("YH9", "H"): 243, ("CaH6", "H"): 215}

mats = []
for d in hyd:
    w2 = fnum(d["omega_2_meV"])
    lam = fnum(d["lambda"])
    if w2 is None or lam is None:
        continue  # 排除 not-reported / NA 行——只用真实 ω̄₂
    key = (d["material"], d["isotope"].split()[0])
    tce = TC_EXP.get(key)
    if tce is None:
        continue
    M2 = lam * w2 * w2
    c13 = c_upper(lam, 0.13)
    ceil = c13 * math.sqrt(M2) / KB
    f = tce / ceil
    mats.append({
        "mat": d["material"], "treat": d["treatment"], "P": d["pressure_GPa"],
        "lam": lam, "M2": M2, "ceil": ceil, "tce": tce, "f": f,
    })

print(f"逐材料表（n={len(mats)}，全部有真实 ω̄₂ 且有实验 Tc 锚，用 C(λ,0.13) 认证上界）：")
print(f"{'material':<8}{'treat':<12}{'P':>5}{'λ':>7}{'M2(meV²)':>10}{'天花板_C13':>10}{'Tc_exp':>8}{'𝔣':>7}{'松倍数':>8}")
for m in mats:
    print(f"{m['mat']:<8}{m['treat'][:11]:<12}{m['P']:>5}{m['lam']:>7.3f}{m['M2']:>10.0f}"
          f"{m['ceil']:>10.1f}{m['tce']:>8.0f}{m['f']:>7.3f}{1.0/m['f']:>8.3f}")

fs = np.array([m["f"] for m in mats])
loose = 1.0 / fs
print(f"\n𝔣 分布：min={fs.min():.3f}  中位={np.median(fs):.3f}  max={fs.max():.3f}")
print(f"松倍数(1/𝔣) 分布：min={loose.min():.3f}  中位={np.median(loose):.3f}  max={loose.max():.3f}")

print(f"\n{'目标Tc':>8}{'剪掉':>7}{'存活':>7}   剪掉材料（treat@P）")
TARGETS_Q2 = list(range(100, 501, 25))
false_neg_total = 0
flip_target = None
for tgt in TARGETS_Q2:
    killed = [m for m in mats if m["ceil"] < tgt]
    alive = [m for m in mats if m["ceil"] >= tgt]
    # 假阴性：真 Tc_exp >= 目标，却被剪掉（ceil < tgt 但 tce >= tgt）
    fneg = [m for m in killed if m["tce"] >= tgt]
    false_neg_total += len(fneg)
    if flip_target is None and len(killed) > 0:
        flip_target = tgt
    names = ", ".join(f"{m['mat']}/{m['treat'][:6]}@{m['P']}" for m in killed) if killed else "（无）"
    fneg_flag = f"  ⚠⚠ 假阴性 {len(fneg)} 个！" if fneg else ""
    print(f"{tgt:>6}K{len(killed):>7}{len(alive):>7}   {names}{fneg_flag}")

print(f"\n自查门 Q2：假阴性（真 Tc_exp≥目标 却被剪掉）总计数 = {false_neg_total}（应为 0；非 0 立即报警）")
if false_neg_total != 0:
    print("  ⚠⚠⚠ 天花板被证伪！停止后续解读。")
print(f"筛子第一次开始剪材料的目标 Tc = {flip_target} K（此前所有材料的天花板都 ≥ 该目标，筛子空转）")

print("\n" + "=" * 112)
print("★ 判词草稿（数字见上，不外推）")
print("=" * 112)
