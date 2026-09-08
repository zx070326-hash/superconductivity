# 跨材料 Migdal vs 顶角修正 Tc 数据检索

调研 agent：vertex-data-other（sonnet / researcher）　日期：2026-09-08。原始产出摘录。WebSearch 配额耗尽，用 arXiv export API + Semantic Scholar + PDF 落盘读全文。

## 找到的

**Mishra, Mori, Margine, npj Comput. Mater. 11, 342 (2025) = arXiv:2507.01897，Table I（p.12，全文核实）**，H3S 200 GPa / Pb，μ*_E = 0.16 / 0.10：

| 材料 | 声子 | λ | λ^V | Tc^FBW | Tc^FBW+μ | Tc^FSR | Tc^ver-FBW | Tc^ver-FBW+μ | Tc^ver-FSR | Tc^exp |
|---|---|---|---|---|---|---|---|---|---|---|
| H3S | har | 2.29 | 1.78 | 235 | 238 | 246 | 222 | 226 | 198 | 175–185 |
| H3S | anh | 1.72 | 0.71 | 189 | 195 | 192 | 178 | 183 | 157 | 175–185 |
| Pb | har | 1.29 | 0.20 | 7.2 | – | 7.2 | 7.2 | – | 7.2 | 7.2 |

（"+μ" = 化学势自洽更新的 FBW。）原文："This methodology ... is implemented in the EPW code and supports both the full-bandwidth (FBW) ... and Fermi-surface restricted (FSR) approaches"。

**Cappelluti, Ciuchi, Grimaldi, Pietronero, Strässler, PRL 88, 117003 (2002) = cond-mat/0105351，MgB2**：σ 带用 2D 抛物线空穴带模型（非全 DFT），E_F^σ ≈ 0.4–0.6 eV，ω_ph/E_F ≈ 0.1–0.2；Fig. 3(c) 给 **Tc(vertex)/Tc(ME) ≈ 1.0–2.5**（λ = 0.6–1.0，随空穴填充趋近带顶增大）——**顶角增强**（小 q、2D）。λ_LDA 0.7–0.9 而 ME 拟合实验需 1.4–1.7，非绝热顶角被提出来填这个缺口。无绝对 Tc 表。

**Schrödi, Cheenicode Kabeer, Aperis, Oppeneer, PRB 104, 094516 (2021) = arXiv:2103.16928，ThFeAsN**：Uppsala 组唯一真实材料；只有一次性非自洽顶角函数图（Fig. 9），无自洽 Tc_vertex；λ₀ = 0.40，ω_log = 17 meV，E_F^eff = 52 meV。

## 没找到
- UppSC 系列其余全是模型系统；PRB 104, 144511 存在性未证实
- Margine 组：Li2AgH6、CeH9 无顶角论文；CaH6/YH6 只有 FBW 无顶角（arXiv:2310.00056）
- Sadovskii 2019/2020 只有闭式无材料数字
- EPW 的顶角功能是否已进公开分支未核实

## 主会话备注
1. 2507.01897 Table I 的 FSR 数字（har 246/198，anh 192/157）与我们从 Figshare 能隙文件提取的（245/197，192/156）一致——P2 §1 的 δ_H 有文献发表值背书。
2. "+μ" 列：化学势更新使 H3S FBW Tc **升高** 3–6 K（235→238，189→195）。我们求解器的 `fix_mu="none"` 给 195.4、`"conserve"` 给 189.7——**方向相反**（none ≈ +μ 值，conserve ≈ 无 μ 值）。要么我们的 mu_shift 符号反了，要么两边"无 μ"的定义不同。留待求解器后续核查；不影响 P1/P2（四因子数字直接来自 EPW 文件）。
3. MgB2 是"动态区主导 → 增强"的真实材料例子，与 H3S"静态区主导 → 压制"形成对照，正好是式 (5.1) 适用边界的两端。
