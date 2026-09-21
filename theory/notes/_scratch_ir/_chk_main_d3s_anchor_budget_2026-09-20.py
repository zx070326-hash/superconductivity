# -*- coding: utf-8 -*-
"""
D3S 锚点决策的**预算**：`h3s-omegalog-convention` 建议把 `Tc_exp` 从 150 K 改到 145.5 K
（Drozdov Fig.2c @200 GPa，**单路**数字化，该 agent 自己标注「承重前需双路复核」）。

按 CLAUDE.md「预算优先」：先问**这个数字改动会翻转哪个下游判决**，再决定要不要为
「把它升到双路」派活。若不翻转，就带限定条件采用，不派活。

下游消费点（唯一）：`deficit_factor_f_headline_2026-09-20.md` §2.1 残差表的 D3S 行，
判词是**残差符号**（「𝔣_pred 系统性高于 𝔣_obs」）。

方向：本脚本解一个**翻号阈值**（双侧求根），不产生单边界判决。
严格性：analytic（给定 `𝔣_pred` 与一个锚点下的 `𝔣_obs`，其余是除法）。
主会话 2026-09-20。
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# 来自 §2.1 表（仓库既有数字，本脚本不重解 Tc_P）
F_PRED   = 0.5203          # δ_FIN·δ_BRD·δ_COU
TC_A, F_OBS_A = 145.5, 0.4873     # 压强匹配锚点（单路数字化）
TC_B          = 150.0             # legacy 锚点（出处未追到一手）

# 𝔣_obs ∝ Tc_exp（天花板与 Tc_exp 无关）⟹ 可直接线性外推
def f_obs(tc):  return F_OBS_A*tc/TC_A
def resid(tc):  return F_PRED/f_obs(tc) - 1.0

print("="*82)
print("D3S 锚点：残差对 Tc_exp 的依赖，与翻号阈值")
print("="*82)
print(f"  {'Tc_exp [K]':>12}{'𝔣_obs':>10}{'残差':>12}   来源")
for tc, src in [(TC_A, "Drozdov Fig.2c @200GPa 单路数字化（建议采用）"),
                (147.75, "两锚中点"),
                (TC_B, "仓库 legacy，出处未追到一手"),
                (152.0, "Errea 自己的 ME 预测值（非实验，仅作标尺）")]:
    print(f"  {tc:>12.2f}{f_obs(tc):>10.4f}{100*resid(tc):>11.2f}%   {src}")

tc_flip = TC_A*F_PRED/F_OBS_A      # resid = 0 ⟺ 𝔣_obs = 𝔣_pred
print(f"\n  ★ 残差翻号阈值：Tc_exp = **{tc_flip:.2f} K**")
print(f"     争议区间 [145.5, 150] 的上端距阈值还有 {tc_flip-TC_B:.2f} K "
      f"（{100*(tc_flip/TC_B-1):.1f}%）")
print(f"     ⟹ 两个候选锚点给出**同号**残差（{100*resid(TC_A):+.2f}% / {100*resid(TC_B):+.2f}%）")

print("""
  ★ 预算判决：
    · 消费点 = §2.1 残差表 D3S 行，判词 = **残差符号**（𝔣_pred 系统性高于 𝔣_obs）。
    · 该判词对锚点在 [145.5, 150] 内的取值**不变**——翻号要 Tc_exp > 155.35 K，
      比任何已报值都高，且 Errea 自己的 ME 预测 152.0 K 都够不到。
    ⟹ **把这个数字升到「双路数字化」的预算 = ∞**（升了也不改任何判决）⟹ **不派活**。
    ⟹ 采用 145.5 K，**但必须随数字携带两条限定**：
       (a) 单路数字化，未双路复核；
       (b) 残差幅度（6.77% vs 3.56%）**确实**依赖锚点 ⟹ **幅度不得引用为精度声明**，
           只有符号可引用。
    ⚠ 反向提醒：CSV 第 9 行因 `ω̄₂` 缺失本就不能进 AD75 回归（该 agent 的判决 B），
       所以这条锚点**只**活在 §2.1 的望远镜记账里，不进头条回归。""")
print("="*82)
