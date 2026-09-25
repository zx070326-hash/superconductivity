#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
把氢化物的天花板余量换算到与 Tier-0 元素/A15 **同一信息类**，再比「谁最紧」。

## 为什么需要这个脚本（承重的约定审计）

`ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md` §6.3 报
氢化物 9 行的 `𝔣 = Tc_exp/天花板`：μ*=0 档 max **0.668**、μ*=0.13 档 max **0.866**。
`ceiling/tier0_eta_ceiling_32_elements_2026-09-25.md` 报 Nb3Ge `r = 0.626`
并称其为「本仓库迄今最紧的天花板数据点」。

**0.668 > 0.626，看起来直接反驳了那句话。但两个数不是同一个量**：
  - §6.3 的 𝔣 用 `C(λ, μ*)`（依赖 λ、依赖 μ* 档），
  - Tier-0 的 r 用 `C_∞`（λ 与 μ* 全盲的**封顶**常数）。
`C(λ,μ*) ≤ C_∞` 逐点成立（GATE M，certificate）=> 同一材料下
`天花板_C(λ) ≤ 天花板_C∞` => `𝔣 ≥ r`。**𝔣 天然偏大，不能与 r 直接比大小。**
（CLAUDE.md W1：两个数不一致时，先问「它们是同一个证书/同一个约定下的量吗」，再问「谁错了」。）

本脚本只做一件事：**用 C_∞ 重算氢化物**，得到可与 Nb3Ge 直接比较的 `r_C∞`。

## 单位
`M2 = λ⟨ω²⟩ = λ·ω̄₂²`（ω̄₂ 是二阶矩的平方根，**不是** ω_log）。
`kTc ≤ C_∞·√M2` => `Tc_ceiling[K] = C_INF * sqrt(lam) * w2[meV] * K_PER_MEV`。

严格性：C_∞ analytic；本比较 numerical-evidence（λ/ω̄₂ 是各文献的 DFT 值，Tc 是实验值）。
"""
import csv
import math
import os

C_INF = 0.1827262477
# 1 meV -> K，从 CODATA 独立算，不硬编码
EV = 1.602176634e-19
KB = 1.380649e-23
K_PER_MEV = 1e-3 * EV / KB

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "..", "..", "data", "hydrides_theory.csv")

# 实验 Tc 锚（theory/data/hydrides_experimental.csv 逐行读取；括号内是压力失配的幅度）
# 注意：理论行的压力与实验锚的压力**并不严格一致**，这是本比较的已知限度，见 §限度。
TC_EXP = {
    "H3S":  (203.0, "155 GPa, resistivity+susceptibility (理论行 @200 GPa)"),
    "D3S":  (150.0, "~150 GPa, resistivity (样品压力一致性未确认)"),
    "YH6":  (224.0, "166 GPa, resistivity"),
    "YH9":  (243.0, "201 GPa, resistivity (177 GPa 的 256 K 出自**已撤稿**论文，不用)"),
}


def ceiling_K(lam, w2_meV):
    """C_infty 档天花板，单位 K。"""
    return C_INF * math.sqrt(lam) * w2_meV * K_PER_MEV


def main():
    print("K_PER_MEV = %.6f K/meV   C_INF = %.10f" % (K_PER_MEV, C_INF))
    print("Tc_ceiling[K] = C_INF * sqrt(lam) * w2[meV] * K_PER_MEV")
    print()

    rows = []
    with open(CSV, newline="", encoding="utf-8") as fh:
        for rec in csv.DictReader(fh):
            w2 = rec["omega_2_meV"].strip()
            mat = rec["material"].strip()
            if w2 in ("", "not-reported", "NA") or mat not in TC_EXP:
                continue
            rows.append((mat, rec["treatment"].strip(), rec["pressure_GPa"].strip(),
                         float(rec["lambda"]), float(w2)))

    hdr = "%-6s %-20s %5s %8s %9s %10s %9s %8s" % (
        "mat", "treatment", "P", "lambda", "w2[meV]", "ceil_C8[K]", "Tc_exp", "r_C8")
    print(hdr)
    print("-" * len(hdr))
    worst = (0.0, "")
    n = viol = 0
    for mat, treat, p, lam, w2 in rows:
        c = ceiling_K(lam, w2)
        tce = TC_EXP[mat][0]
        r = tce / c
        n += 1
        viol += (r > 1.0)
        if r > worst[0]:
            worst = (r, "%s/%s" % (mat, treat))
        print("%-6s %-20s %5s %8.4f %9.2f %10.1f %9.1f %8.4f"
              % (mat, treat, p, lam, w2, c, tce, r))

    print()
    print("=== 判决（同信息类，全部用 C_infty）===")
    print("  氢化物 %d 行：%d 违反；最紧 %s  r_C8 = %.4f" % (n, viol, worst[1], worst[0]))
    print("  Tier-0 A15 最紧：Nb3Ge r = 0.6260")
    print("  Tier-0 元素最紧：In   r = 0.4460")
    if worst[0] < 0.6260:
        print("  => **Nb3Ge 仍是本仓库最紧的点**（氢化物在同信息类下只到 %.4f）。" % worst[0])
        print("     tier0 笔记那句话**成立**，但此前缺一步换算才能这么说。")
    else:
        print("  => **Nb3Ge 不再是最紧的点**，tier0 笔记那句话须订正。")

    print()
    print("=== 为什么不能直接引 §6.3 的 0.668 / 0.866 ===")
    print("  §6.3 的 𝔣 用 C(λ,μ*)，本表的 r 用 C_∞，且 C(λ,μ*) <= C_∞ 逐点成立")
    print("  => 同一材料 𝔣 >= r_C∞。两者比大小是**跨约定比较**，无效。")
    print("  实测放大倍数 𝔣/r 的量级：")
    print("    §6.3 μ*=0    档 max 𝔣 = 0.668  vs 本表 max r = %.4f  => 比值 %.3f"
          % (worst[0], 0.668 / worst[0]))
    print("    §6.3 μ*=0.13 档 max 𝔣 = 0.866  vs 本表 max r = %.4f  => 比值 %.3f"
          % (worst[0], 0.866 / worst[0]))
    print("  这两个比值就是「λ 与 μ* 信息值多少钱」的直接测量。")

    print()
    print("=== ⚠ 同材料 harmonic/anharmonic 的 ceiling 几乎相等：**构造性恒真，不是物理守恒** ===")
    print("  theory/data/hydrides_theory.csv 的 D3S 非谐行 notes 列逐字写着：")
    print("    「M2 与同材料谐波行相符到 4e-7 是**构造性恒真**（冻结矩阵元），不是物理守恒")
    print("      ……**不得用作跨方法免疫性的证据**」")
    print("  机制：两套 α²F 出自同一 EPW 数据集，电子矩阵元被冻结、只换声子；")
    print("        而 M2 = 2∫α²F·ω dω 只吃矩阵元 => 按构造相等。自由度 0。")
    print("  下面印出的离散是 CSV 三到五位有效数字的**舍入残差**，不是物理离散（全精度 4e-7）。")
    print("  => **不得**据此说「Tier 0 对非谐免疫 / 不需要 SSCHA」。")
    print("     该论断只能由求和规则本身（M2=η/M，η 纯电子量，analytic）支持，")
    print("     且仍只是近似：非谐改变平均原子位置 => 间接改变 η，本数据界定不了该二阶效应。")
    for mat in ("H3S", "D3S"):
        sub = [(t, ceiling_K(l, w)) for m, t, p, l, w in rows if m == mat]
        if len(sub) >= 2:
            lo = min(c for _, c in sub); hi = max(c for _, c in sub)
            print("    %-5s ceiling 跨 treatment 相对离散 = %.4f%%  [舍入残差]"
                  % (mat, 100 * (hi - lo) / lo))

    print()
    print("=== 限度（不许越级读）===")
    print("  1. 理论行压力与实验锚压力**不严格一致**（H3S 理论 200 GPa vs 实验 155 GPa）。")
    print("     天花板随压力变化，故这是量级比较，不是逐点认证。")
    print("  2. λ 与 ω̄₂ 是各文献的 DFT/SSCHA 值，非可认证量 => numerical-evidence。")
    print("  3. YH9 的 256 K 读数出自**已撤稿**论文，本脚本用 243 K（201 GPa）。")


if __name__ == "__main__":
    main()
