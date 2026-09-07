# P1 支撑调研：氢化物同位素效应的数据质量与已有公式出处

调研 agent：hydride-exp-data（sonnet / researcher）　日期：2026-09-08
状态：完整（3 条回传：CSV → `theory/data/hydrides_experimental.csv`；正文 → 本文件；JSON → `library/index/P1.json`）。原始产出，未手改。

## 任务理解

主会话要把氢化物同位素系数 α = −d ln Tc / d ln M 分解为库仑 μ*、非谐、非绝热顶角、能量依赖 DOS 四项修正，需要（A）实验数据基线和（B）已有解析公式的出处核实。以下是B部分的公式核对结果，以及对A部分数据表的质量说明汇总。

## A部分数据质量问题汇总

**1. Tc定义不统一。** 不同论文对Tc用onset、midpoint、zero-resistance三种定义混用，且很少在摘要级信息中明确注明。LaH10两个独立发现组（Drozdov 2019 Nature 250K vs Somayazulu 2019 PRL 260K）读数不同，我怀疑主要是定义差异而非样品差异，但未核实原文方法节，这是推断非事实。

**2. H/D样品压力经常不严格匹配。** YH6(166GPa)/YD6(173GPa)、CeH9(116GPa)/CeD9(111GPa)等，同位素替代实验很少能在完全相同的压力点上完成两次独立高压合成。这意味着直接用CSV里的(Tc_H, Tc_D)算α，隐含了"用同一Tc(P)斜率做压力修正"这一假设，而多数原论文报告的α已经是作者做过修正后的值——但修正方法本身未必和理论工作的分解框架自洽，需要team-lead自行判断是否要重新用原始Tc(P)曲线插值再算α。

**3. 一手/二手来源混淆是本次检索最大的弱点。** 因为工具限制为WebSearch而非逐篇WebFetch全文，多条关键数据（LaD10的Tc=180K/α=0.46、CeD9的Tc=67K、部分氘化CeH8D/CeH6D3的数据）只在理论论文或综述转引中出现，我未能独立定位到"第一次报告这个数字"的原始实验论文。CSV中已用"来源存疑"字样标出，这些数字使用前建议至少让另一个agent或team-lead直接打开原文Table确认。

**4. 一例撤稿数据混入文献生态。** Snider et al. 关于YH9(报告为~262K钇超氢化物)的论文已被期刊撤稿，但其α=0.46-0.48数值仍在部分二手引用中流传，容易被误采。

**5. H3S/LaH10磁性测量数据的公开争议。** Minkov et al. Nature Comm 13:3194(2022)关于Hc1的数据可信度受到质疑（arXiv:2409.12211），本表未直接引用其数值，但若后续需要Hc1数据要单独核查。

## B部分：公式与理论出处核对

**B1. McMillan/Allen-Dynes框架下μ*对同位素系数的修正**

出处确认：W. L. McMillan, Phys. Rev. 167, 331 (1968)；P. B. Allen & R. C. Dynes, Phys. Rev. B 12, 905 (1975)；综述 J. P. Carbotte, Rev. Mod. Phys. 62, 1027 (1990)，题为"Properties of boson-exchange superconductors"。

标准Allen-Dynes公式：
Tc = (ω_log/1.2) f1 f2 exp[−1.04(1+λ)/(λ−μ*(1+0.62λ))]

对此求−d ln Tc/d ln M，在忽略ω_log之外M依赖（即λ、μ*、f1f2对M无显式依赖，只通过ω_log∝M^(−1/2)传递）的最简单近似下，得到教科书上常见的形式：
α = ½ [1 − (1.04μ*(1+λ)(1+0.62λ)) / (λ−μ*(1+0.62λ))² × ... ]（严格展开涉及对指数项中μ*(1+0.62λ)项关于ω_log隐式依赖的求导，若μ*本身通过ω_c=ε_F截断而与ω_log耦合，则多一项）

**注意（我的推断，非直接抄录某单一文献的现成公式）**：Carbotte 1990综述里给出的是α随λ、μ*变化的数值图表和麦克米兰原始的近似展开，而不是一个可以直接照抄的闭式代数表达式；不同教材（如Carbotte 1990 Eq附近段落 vs 一些讲义）给出的展开形式在处理μ*是否随ω_log变化（即是否用"固定μ*"还是"固定μ，让μ*=μ/[1+μln(E_F/ω_log)]随M变化"）上并不一致，这直接对应主会话要分解的"库仑修正对α的贡献"这一项——**这是需要team-lead明确选择哪种μ*处理约定的关键分歧点，我不建议我代为决定**。

**B2. 非谐性对同位素系数：Errea-Calandra-Mauri机制**

出处确认：I. Errea, M. Calandra, F. Mauri, Phys. Rev. Lett. 111, 177002 (2013)，arXiv:1305.7123。该文用SSCHA方法计算PdH/PdD，得到α=−0.38（理论），与实验反常值（α≈−0.25至−0.3）符号一致但量级略有差异。机制陈述（来源直接陈述）：非谐效应改变了氢原子的有效声子谱和电子-声子耦合矩阵元对同位素质量的依赖方式，使得原本谐振近似下预期的正常同位素效应被非谐修正翻转符号。

更早的非谐同位素效应理论：检索中未找到比Errea 2013更早、专门针对同位素系数给出解析公式的非谐理论论文（Hui-Allen等我检索未命中明确匹配，可能标题记忆有误，建议team-lead确认具体指的是哪篇，我可以再查）。**这是一个检索缺口，我没有验证到"Hui-Allen"这一具体提法的存在**。

**B3. 非绝热框架下的同位素效应：Grimaldi/Cappelluti/Pietronero**

出处确认（均已核对标题）：
- Pietronero, Strässler, Grimaldi, "Nonadiabatic superconductivity I. Vertex corrections for the electron-phonon interactions," Phys. Rev. B 52, 10516 (1995)
- Grimaldi, Pietronero, Strässler, "Nonadiabatic superconductivity II. Generalized Eliashberg equations beyond Migdal's theorem," Phys. Rev. B 52, 10530 (1995)
- 同作者, "Nonadiabatic superconductivity: electron-phonon interaction beyond Migdal's theorem," Phys. Rev. Lett. 75, 1158 (1995)
- Cappelluti & Pietronero, "Nonadiabatic superconductivity: the role of van Hove singularities," Phys. Rev. B 53, 932 (1996)；同组, Europhys. Lett. 36, 619 (1996)
- **同位素效应专论**：Grimaldi, Cappelluti, Pietronero, "Isotope effect on m* in high-Tc materials due to the breakdown of Migdal's theorem," Europhysics Letters 42, 667 (1998)。

这篇1998年EPL文章的核心论点（来源直接陈述）：非绝热顶点修正会给准粒子有效质量m*带来一个即使远未达到极化子极限也存在的同位素效应，这与极化子机制给出的m*同位素效应是两种不同的物理起源，论文将其用于解释YBCO/LSCO等体系实验观测到的m*同位素效应。

**适用前提（我的判断）**：这组理论的核心变量是ε_F/ω_ph比值（绝热性参数），只有当声子频率与费米能量可比时（如氢化物中的高频H振动模式相对较低的电子费米能标）才需要引入。这与主会话要分解的"非绝热顶角修正"直接对应，但**这些论文本身是针对高Tc铜氧化物/富勒烯体系写的，没有一篇是专门针对氢化物的**——检索中确认1995-1998这批原始论文都未提及氢化物（当时H3S尚未发现）。氢化物语境下的非绝热+非谐联合处理是后来的工作（见B5，2608.14949/2510.25902）。

**未检索到arXiv ID**：1995-1996这五篇论文均早于cond-mat大规模上传arXiv的年代，搜索未返回任何arXiv编号，判断为很可能没有对应预印本，JSON中将arxiv_id标为null。

**B4. 能量依赖DOS / van Hove奇点对同位素系数的影响**

出处确认：S. Villa-Cortés, O. De la Peña-Seaman, "Effect of van Hove singularity on the isotope effect and critical temperature of H3S hydride superconductor as a function of pressure," arXiv:2103.12055，正式发表于 J. Phys. Chem. Solids 161, 110451 (2022)。方法：DFT+线性化Migdal-Eliashberg方程+"functional derivative approach"，把费米面附近能量窗口内N(ε)的形状变化（而非常数N(0)近似）代入α的计算，覆盖H3S的Im-3m相（162-250 GPa）。**来源直接陈述的结论**：N(ε)的强能量依赖会同时重整化α中的电子-声子（正贡献）和电子-电子（负贡献）两部分，与恒定N(0)近似相比显著改善了对实验Tc(P)和α(P)的复现，尤其在高压端。

**一个未完全厘清的溯源问题（需要team-lead留意）**：我在检索"H3S isotope effect Solid State Communications"时，搜索引擎摘要曾把这篇van Hove论文的方法学描述和另一篇更早（2016/2017年, arXiv:1601.01531, Szczęśniak & Durajski）的"The isotope effect in H3S superconductor"论文的标题及α=0.31/2.37数值混在同一段摘要文字里呈现给我，我怀疑这是搜索工具对两篇标题相近论文的摘要做了错误拼接，而不是同一篇论文的两个版本。我没有WebFetch两篇原文来切割清楚，所以在CSV和这里都保留了两个独立条目并标注"归属存疑"，请勿假设α=0.31/2.37这组数字就是Villa-Cortés这篇论文给出的。

**B5. 已合并多项修正讨论氢化物同位素系数的论文**

- **Errea et al., Nature 532, 81 (2016)**（仅非谐/量子核效应，不含非绝热顶角、不含能量依赖DOS）：H→D替代使Im-3m对称化转变压力上移约12 GPa（0K），这是非谐修正对同位素效应的结构侧证据，但**该文本身未给出一个显式的α数值**（据我检索到的摘要信息，未见其报告α；这是我的检索局限，不是"论文没算"的断言，需要核实原文是否含此数值）。**arXiv ID未能确认**。
- **Szczęśniak & Durajski (arXiv:1601.01531, 2016/17)**：仅经典Migdal-Eliashberg（不含非谐、不含非绝热顶角），α=0.31在200GPa处与实验吻合，但无法解释130GPa处的反常高值2.37，作者把2.37归因于DOS效应（衔接到B4）。
- **Villa-Cortés & De la Peña-Seaman (arXiv:2103.12055, 2021/22)**：Migdal-Eliashberg + 能量依赖DOS（不含非谐、不含非绝热顶角），改善了对α(P)全曲线的复现。
- **Durajski, Szczęśniak, Li, Wang, Cho, PRB 101, 214501 (2020)/arXiv:1912.09663**：针对LaH10，经典Migdal-Eliashberg（不含非谐、不含非绝热），α理论=0.458-0.470，与"转引"实验值0.46吻合。
- **Mishra, Belli, Zurek, Margine, arXiv:2608.14949 (2026)**：**同时**含非谐重整化(A_λ)和非绝热顶角修正(R_V)两项，是我检索范围内唯一明确同时处理这两项修正的氢化物论文；发现LaBeH8受顶角修正压低Tc约4K，PdH/D/T强非谐但绝热（顶角修正可忽略）。**该文未纳入能量依赖DOS，也未纳入库仑μ*的独立处理**——即，目前没有任何一篇论文把题目要求的"四项修正"同时纳入同一个氢化物计算框架，这是主会话工作可能填补的空白，建议在报告中明确指出。
- **YH6/YH9的Ding, Chen, Shi, arXiv:2507.03383 (2025)**：非谐+μ*截止频率重整化两项联合（不含非绝热顶角、不含能量依赖DOS），结论是仅非谐仍高估Tc，必须同时调整μ*截止才能接近实验。

**B6. Tolmachev对数收缩论点**

**关键更正/澄清**：我在此前D2方向调研中记录的"Marsiglio 2022"这一出处，本轮反复检索（含直接搜索arXiv作者列表）**仍未能确认其存在**——所有搜索结果只把这个说法追溯到Semenok博士论文（arXiv:2307.13313）的转引，未找到一篇独立、可核实标题的2022年Marsiglio论文专门论证"氢化物中Tolmachev对数收缩导致μ*偏大"。

但检索找到了一个**可以确认、内容高度相关、可以替代或补充**的独立出处：**M. V. Sadovskii, "Antiadiabatic phonons, Coulomb pseudopotential and superconductivity in Eliashberg-McMillan theory," JETP Letters 109, 166 (2019), arXiv:1811.10184**，及其扩展版 arXiv:1908.00718（J. Supercond. Nov. Magn.特刊）。**来源直接陈述**：在离散声子频率模型中，若某个（或几个）声子模式进入"非绝热"区间（声子频率超过某个特征电子能标），则该模式对Tolmachev对数没有贡献，Tolmachev对数只由绝热声子模式的部分贡献决定；由于非绝热声子的存在会抑制总的Tolmachev对数，这会导致μ*偏大，进而显著压低Tc。**适用前提**：这是一个离散多模型（而非连续谱）的解析处理，物理图像与氢化物中高频H声子模式相对电子能标"过快"这一情形直接对应，但论文本身的应用背景是FeSe/STO单层，不是氢化物——**转用到氢化物语境是主会话或此前二手引用做的外推，不是Sadovskii原文的应用对象**。

## 矛盾与不确定点汇总

1. LaH10 Tc读数250K(Drozdov) vs 260K(Somayazulu)——定义差异还是样品差异未核实。
2. LaD10的Tc=180K/α=0.46一手来源未定位，仅见理论论文转引。（主会话补：Drozdov 2019 摘要列 isotope effect，一手就是该文，正文数字待读。）
3. CeD9的Tc=67K@111GPa来源未定位。
4. H3S"isotope effect in H3S superconductor"这一标题下可能存在两篇不同论文被搜索引擎摘要混淆。
5. Errea Nature 2016论文的arXiv ID仍未确认。
6. "Marsiglio 2022"关于Tolmachev对数的论文，多轮检索仍无法证实存在，建议放弃该引用，改用Sadovskii 2019。

## 建议的后续检索方向

- 请另一位agent对LaD10/CeD9/部分氘化CeH8D数据做WebFetch级别的原文定位。
- 若要采用B1的μ*-α解析式，需要team-lead先确定"μ*是否随M变化"这一建模选择。
- 若需要"四项修正同时出现在一篇论文"的更完整检索，建议扩大到2026年最新预印本。
