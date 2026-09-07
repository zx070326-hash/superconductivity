# theory/data · 数据表说明

## ⚠ 2026-09-08 二轮更新：`isotope_curated.md` 取代本目录 CSV 中的同位素相关行

researcher 读本地 PDF 原文后发现一轮摘要级表有两处实质错误：
- **LaD10 Tc=180 K / α=0.46** 不在 Drozdov 2019 原文，是 Durajski 2020 理论值经两篇文献误转引；原文为 α=0.57（fcc）/0.62（P4/nmm），LaD10 Tc ~168 或 ~160 K
- **YD9 无实验数据**（Kong 2021 只做 YH6/YD6）；唯一来源是已撤稿论文
- H3S 作者给的 α_H ≈ 0.3（区间比较），不是两点比 0.44

同位素相关的量一律以 `isotope_curated.md` 为准。

## hydrides_experimental.csv

来源：researcher agent `hydride-exp-data`，2026-09-08，第 1/3 条回传。**原始产出，未手改。**

### 证据等级（agent 自述）
- 数据来自 WebSearch 摘要/二手转述，**未逐篇 WebFetch 原文核对**
- 数字置信度分"多源交叉"与"单源转述"，见各行 notes

### 主会话已补核（2026-09-08，WebFetch arXiv 摘要页）
- **LaH10 / LaD10**：Drozdov et al. arXiv:1812.01561（Nature 569, 528）摘要明确把"isotope effect"列为超导证据之一，所以 LaD10 数据的一手来源**就是这篇**，Durajski 2020 的转引有源头；具体 Tc(D) 与压力在正文图表，仍需读正文确认 180 K vs 169–170 K 哪个对。
- **YH9 / YD9**：Kong et al. arXiv:1909.10482（Nat. Commun. 12, 5075）摘要同样列出 isotope effect；YH9 243 K @201 GPa、YH6 227 K @237 GPa 确认；YD9 数字在正文。
- **YH6 / YD6**：Troyan et al. arXiv:1908.01534（Adv. Mater. 2021）摘要原文：Tc=224 K @166 GPa；**α=0.4**；Bc2(0)=116–158 T，比计算值大 2–2.5 倍；"calculations suggest unusually large impact of the Coulomb repulsion"——这句对 P1 的 δα_μ* 项是直接线索。
- **CaH6**：Ma et al. arXiv:2103.16282，215 K @172 GPa，WHH 外推 Hc2(0)=203 T；正式版 PRL 128, 167001 (2022)（DOI 按记忆填，未核）。

### 使用前必须处理的问题
1. **LaD10 的 α=0.46 溯源最弱**：数字出现在 Durajski et al. PRB 101, 214501 (2020) 理论论文的转引里，一手实验来源（Drozdov Nature 2019 正文是否含 LaD10）未确认。P1 用这个数之前要查原文。
2. **YH9 的 α=0.46–0.48 来自已撤稿论文**（Snider et al.），只能用 Kong et al. Nat. Commun. 2021 的 YH9/YD9（243 K @201 GPa / 172 K @205 GPa），α 需自己算：ln(243/172)/ln(9.0/... ) —— 注意 YH9 的质量比不是 2（Y 不变，只有 H→D），应按 α = −Δln Tc / Δln M_eff，M_eff 的定义本身是 P1 要处理的理论问题。
3. **CeD9 数据（67 K @111 GPa）与部分氘代阶梯（CeH8D 92 K、CeH6D3 82 K）出处未确认**，不建议定量使用，除非补核。
4. **H/D 样品压力不匹配**是所有体系的通病（H3S ~155 vs ~150；YH6 166/183 vs 173；YH9 201 vs 205）。P1 的 α 分解必须把 dTc/dP 的修正显式放进误差预算。
5. **三个同位素数据空白**：CaD6、ThD9/ThD10、LaBeD8 均无实验。P1 的检验集实际只有 H3S、LaH10（若溯源成立）、YH6、YH9、（CeH9 存疑）、PdH 五到六个体系。
6. H3S 130 GPa 处 α=2.37 是 Szczęśniak–Durajski 对已有 Tc(P) 的再分析，反映的是 R3m 相 + 结构相变区，不是单相内的同位素效应；应排除或单独处理。
7. Minkov et al. Nat. Commun. 2022 的 H3S/LaH10 磁化数据有公开可信度争议（arXiv:2409.12211）；本表未用其数据，但引用 Hc1/磁化时注意。

### 后续
- 第 2/3 条（数据质量说明 + 公式核对）、第 3/3 条（文献 JSON）到达后追加到本目录
- 补核清单：LaD10 一手来源；CeD9 出处；Ma et al. 2022 PRL 的 DOI；Chen et al. CeH9 正式版卷期
