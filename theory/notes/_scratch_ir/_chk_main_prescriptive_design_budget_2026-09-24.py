# -*- coding: utf-8 -*-
"""路线 3 预算 — 天花板能不能**反过来指导**计算与合成？三问，每问先算预算。

Q1【筛选力】`M2 = λ⟨ω²⟩ = 2∫α²F·ω dω = η/M`（Hopfield 求和规则）是**纯电子量 ÷ 质量**，
   原则上不需要声子计算。那么 `Tc ≤ C_∞√(M2)` 当预筛能剪掉多少候选？
   在 Cerqueira batch-q（249 个真实材料，独立于氢化物）上实测剪枝率。
   ⚠ 预算若 ≈0（剪不掉）或 ≈100%（把好的也剪了）都是废筛子，要看分布。

Q2【跨方法免疫】天花板对「声子那一半」是否免疫 ⟹ 能否用便宜的谐波声子替代 SSCHA？
   YH6 有三套**独立**声子处理（DFPT 谐波 / SSCHA 非谐 / SPIA 路径积分），同一电子结构。
   测 λ、ω̄₂、M2、天花板各自的跨方法跨度。

Q3【反设计】给定目标 Tc，天花板 + 实测填充因子 𝔣 反解出 `(λ, η/M)` 的**必要条件**曲线。
   这是唯一能直接交给实验/计算的产物：「要 300 K，η/M 至少多少」。
   ⚠ 必要条件不是充分条件——方向必须写清。

严格性：C_∞ 是 `analytic`；C(λ,μ*) 端点 `certificate`；𝔣 与全部材料数 `numerical-evidence`。
"""
import os
import sys
import csv
import math
os.environ["OMP_NUM_THREADS"] = "1"
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

KB = 0.086173332621452
C_INF = 0.1827262477
D = os.path.dirname(os.path.abspath(__file__))
C_UP13 = [(1.0, 0.085109595), (1.5, 0.101992146), (2.0, 0.116687205), (2.5, 0.125402058),
          (3.0, 0.130991149), (3.2953, 0.133361200), (4.0, 0.137762123), (5.0, 0.141888517),
          (6.0, 0.144677171), (8.0, 0.148030265), (10.0, 0.150564456)]
C_UP0 = [(2.0, 0.151236028), (3.2953, 0.164718304), (4.0, 0.168304891), (6.0, 0.173820789)]


def c_up(lam, tbl):
    for L, c in tbl:
        if lam <= L:
            return c
    return C_INF


# ═══════════════ Q1 筛选力 ═══════════════
p = os.path.join(D, "..", "..", "data", "external", "cerqueira6912_batchq_sample.csv")
rows = [r for r in csv.DictReader(open(p, encoding="utf-8")) if r["lambda"] and r["w2av_meV"]]
lam = np.array([float(r["lambda"]) for r in rows])
w2 = np.array([float(r["w2av_meV"]) for r in rows])
M2 = lam * w2 * w2
ceil_inf = C_INF * np.sqrt(M2) / KB
tc_ad = np.array([float(r["Tc_AD_mustar010_K"]) if r["Tc_AD_mustar010_K"] else np.nan for r in rows])

print("=" * 100)
print(f"Q1 筛选力 — Cerqueira batch-q，n={len(rows)}（λ∈[{lam.min():.3f},{lam.max():.3f}]）")
print(f"   M2 = λ⟨ω²⟩ ∈ [{M2.min():.0f},{M2.max():.0f}] meV²   天花板 ∈ [{ceil_inf.min():.1f},{ceil_inf.max():.1f}] K")
print(f"\n{'目标 Tc':>9}{'被剪枝':>9}{'剪枝率':>9}{'存活':>7}   | 存活者中 Tc_AD 最大值（应 ≥ 0，越小说明剪得越干净）")
for tgt in (10, 30, 77, 150, 200, 300):
    kill = ceil_inf < tgt
    surv = ~kill
    mx = np.nanmax(tc_ad[surv]) if surv.any() and not np.all(np.isnan(tc_ad[surv])) else float("nan")
    print(f"{tgt:>7} K{kill.sum():>9}{kill.mean()*100:>8.1f}%{surv.sum():>7}   | "
          f"max Tc_AD = {mx:6.2f} K" + ("  （存活 0 个）" if surv.sum() == 0 else ""))
print("\n★ 声音性（唯一要紧的一条）：天花板是**上界** ⟹ 「天花板 < 目标」蕴含「真 Tc < 目标」")
print("   ⟹ 剪枝**零假阴性**（不会丢掉好材料）。代价只是剪得不够狠（假阳性 = 白算）。")
bad = int((tc_ad[~np.isnan(tc_ad)] > ceil_inf[~np.isnan(tc_ad)]).sum())
print(f"   自查：Tc_AD > 天花板 的行数 = {bad} / {int((~np.isnan(tc_ad)).sum())}  （应为 0；非 0 即天花板被证伪）")
print(f"\n⚠ 限定：batch-q 全部 λ≤{lam.max():.2f}，是**弱耦合尾巴**，不代表氢化物总体。")
print("   这里量到的是『筛子在垃圾堆上剪枝率 100%』，**不是**『筛子能分辨好材料』。")

# ═══════════════ Q2 跨方法免疫 ═══════════════
hyd = list(csv.DictReader(open(os.path.join(D, "..", "..", "data", "hydrides_theory.csv"), encoding="utf-8")))
print("\n" + "=" * 100)
print("Q2 跨方法免疫 — 同一材料、同一电子结构、**不同声子处理**下各量的跨度")


def span(v):
    return (max(v) - min(v)) / min(v) * 100


CONSTRUCTED = []   # M2 跨度小到只能是冻结矩阵元的材料（不得计入免疫性证据）
CONFOUNDED = []    # 各套声子处理不在同一压力上 ⟹ 压力与方法混杂（不得计入）
for mat in ("YH6", "H3S", "D3S", "YH9", "CaH6"):
    sel = [d for d in hyd if d["material"] == mat and d["omega_2_meV"] not in ("not-reported", "NA", "")]
    if len(sel) < 2:
        sel2 = [d for d in hyd if d["material"] == mat and d["omega_log_meV"] not in ("not-reported", "NA", "")
                and d["pressure_GPa"] == sel[0]["pressure_GPa"]] if sel else []
        if len(sel2) < 2:
            print(f"  {mat:<6} 跳过（同压下少于 2 套声子处理，或缺 ω̄₂）")
            continue
        sel = sel2
    # ⚠⚠ 混杂门（2026-09-24 加）：Q2 问的是「**只**换声子处理，天花板动多少」。
    #   若两套不在同一压力上，量到的是「压力 + 方法」的合效应 ⟹ 该材料不能回答 Q2。
    #   （实测抓到 YH9：harmonic@260 GPa vs anharmonic+quantum@255 GPa。）
    pres = sorted({d["pressure_GPa"] for d in sel})
    if len(pres) > 1:
        CONFOUNDED.append((mat, pres))
        print(f"  {mat:<6} ⚠⚠ **剔除：压力不齐** {pres} GPa ⟹ 压力与声子方法混杂，不能回答 Q2")
        continue
    ls = [float(d["lambda"]) for d in sel]
    ws = [float(d["omega_2_meV"]) for d in sel]
    ms = [l * w * w for l, w in zip(ls, ws)]
    cs = [C_INF * math.sqrt(m) / KB for m in ms]
    c13 = [c_up(l, C_UP13) * math.sqrt(m) / KB for l, m in zip(ls, ms)]
    # ⚠ 构造性门（2026-09-24 加）：M2 跨度 < 0.1% 而 λ 跨度 > 10% ⟹ 不是物理免疫，
    #   是冻结矩阵元/本征矢的指纹（Σ_ν λ_ν ω_ν² 逐模守恒）。矩族邻居检验已在
    #   memory `moment-family-neighbour-test-for-constructed-invariance` 上判实。
    #   方向：该门只能**下调**正面主张，永不上调 ⟹ 单边安全。
    bad = span(ms) < 0.1 and span(ls) > 10.0
    if bad:
        CONSTRUCTED.append(mat)
    print(f"  {mat:<6} n={len(sel)} 套 @{pres[0]} GPa（{', '.join(d['treatment'] for d in sel)}）")
    print(f"         λ 跨 {span(ls):5.1f}%   ω̄₂ 跨 {span(ws):5.1f}%   |   "
          f"**M2 跨 {span(ms):4.1f}%   天花板_∞ 跨 {span(cs):4.1f}%   天花板_C13 跨 {span(c13):4.1f}%**"
          + ("   ⚠⚠ 构造性恒真，**不计入证据**" if bad else "   ✓ 百分量级，可计入"))
if CONSTRUCTED:
    print(f"\n⚠⚠ 构造性剔除 {CONSTRUCTED}：M2 跨度 <0.1% 而 λ 跨度 >10%，这是**冻结矩阵元**的指纹，")
    print("   不是物理免疫性。真做 SSCHA 时本征矢也变，一致度应是**百分量级**。")
    print("   （撤回本会话早先『H3S M2 跨 0.0% ⟹ 天花板完全免疫』的表述——那是自证，不是证据。）")
if CONFOUNDED:
    print(f"⚠⚠ 混杂剔除 {CONFOUNDED}：两套不在同一压力 ⟹ 量到的是压力+方法的合效应。")
print("\n★ 读法（已按上面两道门收窄）：在**未被剔除**的材料上，")
print("  M2/天花板的跨度 ≪ λ/ω̄₂ 的跨度 ⟹ 天花板对**声子计算的精度**不敏感。")
print("  ⚠⚠ 收窄后**干净样本只剩 YH6 一个材料**（n=3，同压 165 GPa）⟹ 这条结论的证据基础是")
print("     **单材料**，等级 `numerical-evidence`，**不得**写成『天花板对声子方法免疫』这种一般命题。")
print("  ⟹ **筛选阶段可以用便宜的谐波 DFPT，不必跑 SSCHA/PIMD**（后者贵 1–2 个量级）。")
print("  ⚠ 限定：三套共享同一电子结构（同泛函/结构/矩阵元）⟹ 检验的是『声子算法』，")
print("     **不是**『电子结构近似』。M2 对后者极敏感——见下面 H3S 两组数据的对照。")
mish = [d for d in hyd if d["material"] == "H3S" and d["source_author"] == "team-internal"]
err = [d for d in hyd if d["material"] == "H3S" and "Errea" in d["source_author"] and d["pressure_GPa"] == "200"]
if mish and err:
    m_m = float(mish[0]["lambda"]) * float(mish[0]["omega_2_meV"]) ** 2
    m_e = float(err[0]["lambda"]) * float(err[0]["omega_log_meV"]) ** 2
    print(f"     Mishra(EPW) M2={m_m:.0f} vs Errea(ω_log 代) M2≥{m_e:.0f} meV²，"
          f"即便按 H3S 实测 ω̄₂/ω_log=1.207 补偿，Errea→{m_e*1.207**2:.0f}，仍差 "
          f"{(m_m-m_e*1.207**2)/m_m*100:.0f}%。")
    print("     ⟹ **M2 是电子结构的指纹**：两套算的 M2 不一致 ⟹ 分歧在电子侧，不在声子侧。")
    print("        这条本身就是个可交付的诊断工具（免费，只需两组已发表的 λ 和 ω̄₂）。")

# ═══════════════ Q3 反设计 ═══════════════
print("\n" + "=" * 100)
print("Q3 反设计 — 给定目标 Tc，反解 (λ, η/M) 的**必要条件**")
print("   kTc ≤ C(λ,μ*)·√(η/M)  ⟹  η/M ≥ (kTc / C(λ,μ*))²        [必要，非充分]")
print("   加上实测填充因子 𝔣（真 Tc / 天花板）得**工程估计**： η/M ≳ (kTc / (𝔣·C))²")
F13 = 0.799       # 本轮实测中位（μ*≥0.13 档，7 个 ω̄₂ 真值行）
print(f"\n   用本轮实测 𝔣_C13 中位 = {F13:.3f}（n=7 氢化物，numerical-evidence）")
print(f"\n{'λ':>6}{'C(λ,.13)↑':>11} | " + "".join(f"{f'Tc={t}K':>22}" for t in (200, 300, 400)))
print(f"{'':>6}{'':>11} | " + "".join(f"{'必要 η/M':>11}{'工程 η/M':>11}" for _ in range(3)))
for L in (1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0):
    c = c_up(L, C_UP13)
    cells = ""
    for t in (200, 300, 400):
        nec = (KB * t / c) ** 2
        eng = (KB * t / (F13 * c)) ** 2
        cells += f"{nec:11.0f}{eng:11.0f}"
    print(f"{L:6.1f}{c:11.9f} | {cells}")
print(f"{'∞':>6}{C_INF:11.9f} | " + "".join(f"{(KB*t/C_INF)**2:11.0f}{(KB*t/(F13*C_INF))**2:11.0f}"
                                            for t in (200, 300, 400)))
print("\n   参照（本轮实测的真实材料 η/M，单位 meV²）：")
for d in hyd:
    if d["omega_2_meV"] not in ("not-reported", "NA", ""):
        print(f"      {d['material']:<6} {d['treatment'][:12]:<13} @{d['pressure_GPa']:>4} GPa  "
              f"λ={float(d['lambda']):.3f}  η/M = {float(d['lambda'])*float(d['omega_2_meV'])**2:7.0f}")
print("\n★ 方向（必须随表一起说）：这是**必要条件**。η/M 达标 ⟹ 天花板允许，**不保证**做得到；")
print("  η/M 不达标 ⟹ **确定性排除**，无论声子谱长什么样、无论用什么方法算。后者才是承重的一半。")
