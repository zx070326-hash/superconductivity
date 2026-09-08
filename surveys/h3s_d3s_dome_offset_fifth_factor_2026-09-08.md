# H3S/D3S 圆顶错位与"第五因子"（结构侧量子核效应）的经验输入核验

调研 agent：researcher　日期：2026-09-08　服务对象：`theory/derivations/P1b_structural_fifth_factor.md` v0
检索手段：本地 PDF 全文逐页阅读（`library/pdf/P1`、`library/pdf/D3`）+ WebFetch 下载 arXiv PDF 后本地 Read + arXiv API + Crossref/Semantic Scholar 元数据。
**本轮 WebSearch 配额已用尽（200/200）**，全部检索通过 WebFetch/arXiv API 完成；这限制了"广撒网找竞品论文"的覆盖面，见 §6。

置信度标注约定：
- **[文本]** = 一手论文正文/图注/表格的文字直接给出的数字（最高）
- **[表]** = 一手论文表格中的数字
- **[读图]** = 我从 PDF 渲染图像上读的，附误差估计（最低，勿用于亚 10% 精度的拟合）
- **[二手]** = 论文 A 转述论文 B 的数字，且我未在 B 的原文中找到该数字
- **[自算]** = 我用上面的数字算的，非任何论文原文

---

## 0. 执行摘要：三条会改动 P1b 立论的发现

**(A) "低压 α ≈ 0.44"这个数字很可能是样品质量假象，不是结构污染。**
Drozdov 2015 **从未报过 0.44**——全文唯一的 α 陈述是 `α_H ~ 0.3`（P>170 GPa 区间拟合）。0.44 是我方在 155 GPa 用 (203 K, 150 K) 两点自算的。而 Minkov et al. 2020（同一个 Eremets 组，改用单质 S + 过量 H₂ 直接合成，避免贫氢副产物污染）明确指出**此前报道的 D3S ~150 K 是被污染样品的低估值，干净样品的 D3S 最高 Tc = 166 K @ 157 GPa** [文本，经 Crossref + Semantic Scholar 两路交叉核对摘要]。
用 (H3S 203 K @155 GPa, D3S 166 K @157 GPa) 重算：**α = −ln(166/203)/ln2 = 0.29** [自算]。
即：**低压端的 α 一旦换用干净 D3S 样品，就从 0.44 掉到 0.29，与高压端的 0.30 几乎相同**。P1b §0 的出发点"低压 0.44 vs 高压 0.30，差 0.14 待解释"可能根本不存在。
⚠️ 这条的唯一漏洞：我**没能拿到 Minkov 2020 全文**（Wiley 403、PMC 被 CAPTCHA 拦），所以不知道他们同批次的 **H3S** 峰值是多少。若他们的 H3S 也被"boost"到 203 K 以上，α 会回升。**这是第一优先级的后续动作**（§6.1）。

**(B) 对称化压力的同位素移动 ΔP_c，四个一手来源给出 12 / 16 / 17–18 / 30–35 GPa，没有一个是 30 GPa 的"共识值"。**
详见 §2 表。最常被引的 Errea 2016 给的是 **12 GPa**（103 → 115）。唯一给出 ~25–35 GPa 的是 Taureau 2024，且那是"局域偶极矩涨落"探针的值，同一篇论文的"密度分裂"探针只给 17–18 GPa。

**(C) 圆顶峰位错位 ~30 GPa 这个前提本身缺乏一手支撑。**
没有任何一手来源说 D3S 峰在 180 GPa。可溯源的峰位是：H3S ~150–155 GPa [文本]，D3S **157 GPa** [文本，Minkov 2020]，错位 ≈ 0–10 GPa。Bianco 2018 转述实验 kink 为 "150 GPa (H3S) / 160 GPa (D3S)"，错位 10 GPa [二手]。"180 GPa"来自我方上一轮对 Drozdov Fig. 2c 的读图，而该面板 D3S 顶部在 150–200 GPa 是**宽平台**（~140–152 K），峰位读图不确定度 ≥ ±20 GPa。

**对 P1b 的净影响**：分解的数学（式 (1)–(3)）没问题，但**要被解释的目标量（0.14 的低-高压 α 差）和主要输入（ΔP_c ≈ 30 GPa）两个都需要重新确定**。好消息是式 (3) 依然给出一个自洽且可证伪的预言，见 §5.4。

---

## 1. Drozdov et al. 2015：H3S/D3S 的 Tc(P) 与 α

**来源**：arXiv:1506.08190v1（本地 `library/pdf/P1/2015_A_P_Drozdov_1506.08190.pdf`，共 19 页含 Extended Data）；发表版 Nature **525**, 73–76 (2015)，DOI 10.1038/nature14964。

### 1.1 论文报出的 α（唯一一处）

正文 p.3 第 (2) 条，逐字：

> "A strong isotope effect: T_c shifts to lower temperatures for sulfur deuteride indicating phonon-assisted superconductivity (Fig. 2b,c). The BCS theory gives dependence of critical temperature on the atomic mass M as T_c~M^-α where α_H ~0.5. **Comparison T_c's in the pressure range P>170 GPa (Fig. 2c) gives α_H ~0.3.**"

- **α_H ~ 0.3**，明确限定在 **P > 170 GPa**，且是"区间比较"，**没有给出点对列表、没有给拟合细节、没有不确定度** [文本]。
- 正文 p.4 另引 Errea 2015：非谐把 α 从 0.5 降到 0.35，"both in agreement with the present work" [文本]。
- **全文没有 0.44、没有任何低压 α 值。** 我方 P1b §0 引的"低压 α ≈ 0.44"是自算量，**不是 Drozdov 报的**。这一点与已有 memory 记录一致，此处再次全文确认。

### 1.2 文字直接锚定的 (P, Tc) 点

| 体系 | P (GPa) | Tc (K) | 出处 | 置信度 |
|---|---|---|---|---|
| H3S | 145 | 170（退火前）→ ~195（退火后） | Fig. 2a 图注 + 正文 | [文本] |
| H3S | 141 | 195 | Fig. 2b 曲线标注 | [文本] |
| H3S | 155 | **203(1) onset**（ZFC 磁化率）；外推 M(H) 极小值给 203.5 K | Fig. 4a、Fig. 4e、正文 p.3(4) | [文本] |
| H3S | 195 | — （磁场测量用，Fig. 3b）| Fig. 3b | [文本] |
| D3S | 155 | **150** | Fig. 2b 曲线标注 | [文本] |
| D3S | （未注压力）| ~70（未退火，D₂S）| Ext. Data Fig. 2d | [文本] |

### 1.3 Fig. 2c（Tc–P 图，H3S 实心 + D3S 空心红）

图注逐字："**Dependence of T_c on pressure. The data on annealed samples are presented. Open color points refer to sulfur deuteride, and filled points to sulfur hydride. Magenta point was obtained in magnetic susceptibility measurements (Fig. 4a).**"

读图结果（面板尺寸小、刻度稀疏，**误差 ±8 GPa / ±10 K**）[读图]：

- H3S：~145 GPa 处黑星 ~195 K；~152 GPa 黑方 ~190 K；~155 GPa 品红方 203 K；170–225 GPa 一簇 ~165–190 K，缓慢下降；最高点在 ~150–160 GPa。
- D3S：~140 GPa 起 ~90–110 K 快速上升；**150–200 GPa 是一个宽平台，~140–152 K**；250–255 GPa 降到 ~125–130 K。
- **关键判断：D3S 顶部是平台不是尖峰，"峰在 180 GPa"无法从此图判定**，读图不确定度 ≥ ±20 GPa。上一轮 survey（`h3s_tc_of_p_alpha_of_p_2026-09-08.md` 第 14 行）记的"D3S ~180 GPa ~150 K（峰）"应降级为"平台内的一个点"，不应作为峰位使用。

### 1.4 Drozdov 自己对结构相变压力的说法（重要，与常引的 150 GPa 冲突）

正文 p.4 [文本]：
> "The phase transition between R3m and Im-3m structures is predicted at **~180 GPa**. This likely relates to the changes in Raman spectra observed in the **170-180 GPa** (see Extended Data Fig. 3)."

Extended Data Fig. 3 图注文字（PDF p.16）[文本]：
> "Raman spectra of sulfur hydride after the release of pressure from 208 GPa at room temperature. ... **There is an apparent phase transition at ≈180 GPa.**"
> "...the phase transformation are at different pressures: 160 GPa for sulfur and **180 GPa for sulfur hydride**."

即 **Drozdov 组自己的拉曼证据把 H3S 的结构相变放在 ≈180 GPa，不是 150 GPa**（拉曼谱在 126/143/155/169/177/184/208 GPa 测得，室温、卸压过程）。这与 Bianco 2018 转述的"150 GPa"直接冲突，见 §4.2。

### 1.5 未找到

- **未找到** Drozdov 2015 中任何 (P, Tc) 数值表；仅有图。
- **未找到** D3S 的磁化率测量（arXiv v1 的 Fig. 4 只有 H3S；Eremets 2022 综述 arXiv:2201.05137 提到发表版 Supplementary Fig. 6d 有 D3S 的 ZFC 台阶在 **136 K**，压力未注 [二手]）。
- **未找到** 逐点 α(P)。若 P1b 需要 α_app(P) 的实验曲线，只能从 Fig. 2c 读图重建，精度 ±0.1 量级——**建议用 WebPlotDigitizer 处理高分辨率原图，不要用本工具读数**。

---

## 2. Errea et al. 2016 及后续：对称化压力的同位素移动 ΔP_c

**arXiv ID 已核实：arXiv:1512.02933v1**（本地 `library/pdf/P1/2016_Ion_Errea_1512.02933.pdf`），
"Quantum Hydrogen-Bond Symmetrization and High-Temperature Superconductivity in Hydrogen Sulfide"，Nature **532**, 81–84 (2016)，DOI 10.1038/nature17175。
（此前 memory 里"arXiv ID 未核实、疑似 1607.02348"的疑点**已解决并更正**：正确 ID 是 1512.02933，1607.02348 是另一篇。）

### 2.1 Errea 2016 的数字

- 静态晶格（经典核）R3m→Im-3m：**175 GPa** [文本, p.2]
- 量子核（SSCHA, PBE, 0 K）：**H3S 103 GPa，D3S 115 GPa** [文本, p.4 + Fig. 3 图内标注]
  > "Our calculations show that this symmetry breaking occurs at **103 GPa in H3S and 115 GPa in D3S** (see Fig. 3). The higher transition pressure in D3S is due to weaker quantum effects."
- **量子核把 H3S 的对称化压力压低 72 GPa**（摘要逐字："the quantum nuclear motion lowers the symmetrization pressure by 72 GPa"）；D3S 对应 175−115 = 60 GPa [自算]
- **ΔP_c = 115 − 103 = 12 GPa** [自算，基于两个文本值]
- Fig. 3(b) 相图条带：Quantum H3S 在 90–103 GPa 为 R3m、103–200 GPa 为 Im-3m；Quantum D3S 在 90–115 为 R3m、115–200 为 Im-3m；Classical 在 90–175 为 R3m、175–200 为 Im-3m；另标注 <112 GPa H3S 进入 C2/c 相 [文本/图注]
- **Errea 自己算的同位素系数（同一相 Im-3m 内、同压比较）**：p.4 逐字
  > "We also find an isotope coefficient α = −[lnT_c(D3S) − lnT_c(H3S)]/ln2 for H→D substitution of **α = 0.35 at 210 GPa** and **α = 0.40 at 155 GPa** in good agreement with experiment (see Fig. 4(c))."
  **这一条对 P1b 极重要**：这是**纯配对** α（两个同位素都被强制放在 Im-3m 相、同压），而它**本身就随压力上升**（155 GPa 处 0.40 > 210 GPa 处 0.35）。也就是说"低压 α 偏大"至少有一部分是配对物理，不是结构污染。P1b 假设 α_pair 在立方相区弱 P 依赖（§5 第 3 条），与此有张力。
- 计算参数：μ* = 0.16 固定，512 个 Matsubara 频率，各向同性 Migdal-Eliashberg，SSCHA 用 3×3×3 超胞（H3S 与 D3S 均在 Im-3m 相）[文本, Methods]
- 作者自陈的不确定度 [文本, p.4]：
  > "the predicted transition pressure depends on the choice of the exchange correlation functional. ... we cannot exclude a small error on the transition pressure."
- 作者对实验的解读 [文本, p.4]：
  > "The sudden drop in T_c measured for D3S below 150 GPa, which is not present in H3S, could arise either from higher temperatures to anneal the D3S sample or could indicate the symmetry breaking that we predicted at 115 GPa."
  → 注意作者自己给了**两个**可能解释，其中第一个（退火温度/样品质量）正是 Minkov 2020 后来证实的那个。

### 2.2 Bianco et al. 2018（同一课题组的后续，方法更精，专门做 H 与 D 两种硫化物）

**arXiv:1802.07968v1**，Phys. Rev. B **97**, 214101 (2018)，DOI 10.1103/PhysRevB.97.214101。
（下载于 `C:\Users\charles\.claude-apexin\projects\D--physics/abaf4c9b-.../tool-results/webfetch-1788864617608-77419k.pdf`）
方法升级：用 SSCHA **自由能 Hessian**（Bianco 2017, PRB 96, 014111）算二阶相变临界压力，精度远高于 Errea 2016 的"自由能有限差分"法。

Fig. 3（PDF p.6）图内标注 + 正文，全部 [文本/图注]：

| 体系/泛函 | 经典（谐波，无 ZPE）P_c^cl | 量子（自由能 Hessian）P_c |
|---|---|---|
| H3S, PBE | **173 GPa** | **91 GPa** |
| D3S, PBE | **173 GPa** | **107 GPa** |
| H3S, BLYP | **256 GPa** | **114 GPa** |

正文逐字（p.6）：
> "In the classical case, the critical pressure is identical for H3S and D3S, with **P_c^cl = 173 GPa**. ... The critical pressures drastically reduce, with **P_c = 91 GPa and P_c = 107 GPa** for H3S and D3S, respectively (**around 10 GPa smaller than the ones presented in Ref. 10**). ... a consistent isotope effect appears with a **critical pressure 16 GPa larger for D3S**, which brings the result for the heavier isotope closer to the classical limit."

- **ΔP_c = 16 GPa** [文本，作者直接说出]
- 经典极限下 H3S 与 D3S 的 P_c **严格相同**（谐波动力学矩阵只差质量因子，无 ZPE 的经典 EoS 对两者一样）——这是"ΔP_c 纯粹是量子核效应"的干净理论证明 [文本]
- 泛函敏感性：谐波层面 PBE 与 BLYP 差 83 GPa，但**含量子非谐后只差 23 GPa**（91 vs 114）[文本]。即 ΔP_c 的绝对位置有 ~20 GPa 的泛函系统误差，但**同位素差** 16 GPa 相对稳健（两同位素用同一泛函，误差大部分抵消）。

### 2.3 Taureau et al. 2024（PIMD + QMC，给出最大的 ΔP_c）

**arXiv:2307.15684v2**，npj Computational Materials **10**, 56 (2024)。
（下载于 `...\tool-results\webfetch-1788864621571-q2zrr2.pdf`）

**Table I（PDF p.7）逐格转录** [表]："Transition pressures and volumes for the local moment formation yielded by PIMD according to the electronic description (DFT-BLYP or QMC), the probe used (density or local moment fluctuations), and the isotope considered (H3S or D3S)."

| | DFT-BLYP H3S Fluc. | BLYP H3S Dens. | BLYP D3S Fluc. | BLYP D3S Dens. | QMC H3S Fluc. | QMC H3S Dens. | QMC D3S Fluc. | QMC D3S Dens. |
|---|---|---|---|---|---|---|---|---|
| V_c [a₀³] | 104.2 | 99.6 | 97.6 | 96.6 | 100.9 | 96.3 | 95.6 | 93.6 |
| **P_c [GPa]** | **130** | **153** | **165** | **171** | **126** | **152** | **156** | **169** |

由此 [自算]：
- BLYP：ΔP_c(Fluc.) = 165−130 = **35 GPa**；ΔP_c(Dens.) = 171−153 = **18 GPa**
- QMC：ΔP_c(Fluc.) = 156−126 = **30 GPa**；ΔP_c(Dens.) = 169−152 = **17 GPa**

正文/摘要 [文本]：
- 摘要："Nuclear quantum effects play a major role in a significant reduction (**≈ 50 GPa**) of the classical ferroelectric transition pressure at 200 K and in a **large isotope shift (≈ 25 GPa)** upon hydrogen-to-deuterium substitution of the local dipole formation pressure, **in agreement with the corresponding change in the T_c maximum location**."
- 正文 p.7："we show that the hydrogen-to-deuterium substitution brings about an increase of the local polarization formation pressure that falls into the **[17-35] GPa** range."
- 正文 p.5："We estimate the peak splitting to take place in H3S at a volume of 99.6 a₀³ for DFT-BLYP, and of 96.3 a₀³ for QMC. ... these volumes correspond to pressures of **153 GPa** for DFT-BLYP and of **152 GPa** for QMC. These values ... **are in good agreement with the position of the maximum T_c measured in experiments**."

**⚠️ 关键概念区分（P1b 若引用此文必须写清）**：Taureau 定义了**两个不同的压力**，不能混用——
1. **局域偶极形成压力**（Table I，126–171 GPa）：局域、瞬时的对称破缺（预形成局域矩），**这个才与 Tc 极大值位置对应**；
2. **铁电相变压力**（全局 Im3m–R3m 长程偶极序）：正文 p.7 明确"**does not correspond to** the ferroelectric transition pressure ... which happens at **lower pressures**"，200 K 下 V_ferro = 117 a₀³，**P_ferro ≈ 82 GPa**。

也就是说：**Taureau 的"真正的 R3m→Im-3m 结构相变"在 82 GPa，远低于圆顶；与圆顶峰位重合的是"局域偶极形成"这个 crossover，不是热力学相变。** P1b §5 最后一条把 D(u) 的物理归因于"R3m 相 H 偏离对称位 → Tc 低"，在 Taureau 的图像下需要改写为"局域偶极/局域矩形成"。

### 2.4 Cherubini et al. 2026（最新，PIMD + MACE MLIP，只做 H3S）

**arXiv:2602.00833v2**（本地 `library/pdf/D3/2026_Marco_Cherubini_2602.00833.pdf`），标题 "Displacive quantum critical point in superconducting hydrides: The case of H3S"。

- **位移型量子临界点 p_QCP = 134 ± 2 GPa**（T→0 外推），属 **4D Ising** 普适类，ν ≈ 0.483 ± 0.014 [文本, p.1/p.4]
- T = 200 K 时位移临界压力 = **124 ± 2 GPa** [文本, p.3]
- 量子核效应相对经典 MD 把相变压力压低 **≈ 50 GPa**（200 K）[文本, p.3]
- **对 P1b 最有杀伤力的一句**（摘要 + p.1，逐字）：
  > "theoretical attempts to study this transition failed to locate it at a pressure compatible with the position of the T_c peak, the most reliable theories placing the transition pressure **more than 40 GPa below** [20–22]."
  > "the experimental T_c peak falls into a centrosymmetric region of large nuclear quantum fluctuations above the displacive QCP"
  即：**主流 SSCHA/PIMD 计算都把结构相变放在离 Tc 峰 40+ GPa 以下的地方**，Tc 峰**不在**相变线上，而在对称相一侧。
- **本文只算 H3S，没有 D3S**，因此不提供 ΔP_c [已通读 p.1–6 确认]

### 2.5 ΔP_c 汇总表（P1b 的核心输入）

| 来源 | 方法 | P_c(H3S) | P_c(D3S) | **ΔP_c** | 置信度 |
|---|---|---|---|---|---|
| Errea 2016 (1512.02933) | SSCHA 自由能差, PBE, 0 K | 103 | 115 | **12** | [文本]×2 |
| Bianco 2018 (1802.07968) | SSCHA 自由能 Hessian, PBE, 0 K | 91 | 107 | **16** | [文本]，作者亲述 |
| Taureau 2024 (2307.15684) | PIMD, 密度分裂探针, BLYP/QMC | 152–153 | 169–171 | **17–18** | [表] |
| Taureau 2024 | PIMD, 局域矩涨落探针, BLYP/QMC | 126–130 | 156–165 | **30–35** | [表] |
| Taureau 2024 摘要综述值 | — | — | — | **≈25（区间 17–35）** | [文本] |
| Bianco 2018 转述实验 kink | Tc(P) 拐点 | 150 | 160 | **10** | [二手]，见 §4.2 |
| 经典（谐波无 ZPE）理论极限 | — | 173 | 173 | **0** | [文本] |

**建议 P1b 采用**：中心值 **ΔP_c = 16 GPa**（Bianco 2018，方法最精、作者直接给出、与 Errea 2016 同组同泛函可比），敏感性区间 **10–35 GPa**。**不建议**直接用 30 GPa 作为中心值——只有 Taureau 的一个探针支持，且该探针测的是局域矩涨落而非结构相变。

---

## 3. Einaga et al. 2016（XRD）

**arXiv ID 已核实：arXiv:1509.03156v1**，"Crystal Structure of 200 K-Superconducting Phase of Sulfur Hydride System"，Nature Physics **12**, 835–838 (2016)，DOI 10.1038/nphys3760。
（下载于 `...\tool-results\webfetch-1788864697024-vgnf48.pdf`，8 页）

### 3.1 最重要的一条：XRD **无法**区分 R3m 与 Im-3m

正文 p.2 逐字 [文本]：
> "the obtained pattern and the calculated ... show the good agreement as follows: R_wp = 0.40%, S = 1.92 **in bcc lattice** and R_wp = 0.41%, S = 1.94 **in hexagonal lattice**. It is because that **both of them are of same periodicity in sulfur atoms** of H3S, which has bcc lattice."

即：两种结构模型的 Rietveld 拟合质量几乎相同（0.40% vs 0.41%），因为 **H 原子是弱散射体、S 亚晶格在两相中几乎相同**。
**结论：Einaga 2016 没有、也无法给出实验测定的 R3m→Im-3m 转变压力——无论对 H3S 还是 D3S。** 任何把"实验 XRD 测到 150 GPa 相变"归给 Einaga 的引用都是错的。

Bianco 2018 也独立确认了这一点，正文 p.2 [文本]：
> "Early x-ray diffraction (XRD) measurements³ confirmed the bcc arrangement of the S atoms at high pressure, **but could not distinguish between the R3m and Im3m phases, as H atoms are very weak scatters and their position cannot be determined.**"

### 3.2 D3S 的 XRD：有，但只在 173–190 GPa

- **有 D3S（sulfur deuteride）的 XRD** [文本]：Fig. 1b、Fig. 2b，压力点 **173、180、190 GPa**（室温）
- 正文 p.2–3："The XRD patterns of sulfur deuteride obtained on increase and decrease from 173 GPa to 190 GPa then 180 GPa are shown in Fig. 2b. **The pattern doesn't change by pressure except the shift of peak positions with pressure changes.**"
- **D3S 在 173 GPa 做了 13–300 K 的变温 XRD + 电阻同步测量**："the XRD patterns kept at all measured temperature region and the resistance measurement shows the drop indicates superconducting transition below 150 K. Therefore, the structure in the normal and SC states are the same - **no visible structural transition at SC transition**." [文本]
- **D3S 没有 173 GPa 以下的 XRD**——因此**无法**从 Einaga 得到 D3S 的低压结构信息。

### 3.3 H3S：92–150 GPa 卸压，未见结构变化

- Fig. 2a：150 → 133 → 123 → 111 → 103 → 92 GPa 卸压序列，"**No structural change was observed in main peaks.**" [文本]
- 正文 p.3 [文本]：
  > "The present data agree with the previous data² above 133 GPa. **T_c suddenly decreased from 123 GPa** and settled down. It seems that a phase transition from high-T_c phase to low-T_c phase occurs. In fact, the resistance at room temperature also increases suddenly at the pressure (Fig. 3a), but **no evidence of this structural phase transition is detected in XRD measurements.**"

  → Einaga 把 Tc 骤降放在 **123 GPa**（不是 150、不是 180），且明确说 XRD 看不到对应的结构变化。

### 3.4 晶格常数（可用于 EoS/体积-压力换算）

- H3S @150 GPa 室温：a_bcc = **3.089(1) Å**；hexagonal a = **4.379(10) Å**, c = **2.663(9) Å** [文本]
- D3S @173 GPa 室温：a_bcc = **3.0468(2) Å**；a = **3.351(3) Å**, c = **2.659(3) Å**（R_wp = 0.28%, S = 0.55）[文本]
- Birch EoS 拟合（Fig. 2c，用 3 个实验点：150 GPa 硫氢化物 + 173、190 GPa 硫氘化物）：**B₀ = 486(35) GPa, B₀' = 4（固定）** [文本]

### 3.5 Fig. 3c 的 Tc(P) 读图（低置信度）

图注：黑圆 = 硫氢化物，红方 = 硫氘化物；**实心 = 本工作（升压与降压），空心 = Drozdov 2015 的升压数据**。
读图 [读图，±5 GPa / ±8 K]：
- H3S 实心：~111/95、~121/117、~133/168、~150/189
- H3S 空心（Drozdov）：~123/195 起，~130–135 达 ~200–203，随后缓降至 ~220/170
- D3S 空心：~113/90、~120/78、~123/110、~133/125、~140/143、~150/152、~160/145、~190/141、~240/135
- D3S 实心：~172/145、~190/134

由 H3S 实心点算上升沿对数斜率 [自算，供 P1b 式 (3) 用]：
- 111→133 GPa：d lnTc/dP = ln(168/95)/22 = **0.026 /GPa**（陡沿）
- 133→150 GPa：ln(189/168)/17 = **0.0069 /GPa**（接近顶部）
- 峰顶（~150–157 GPa）：≈ **0**

---

## 4. 有没有人已经做过"结构 vs 配对"的分解？

### 4.1 结论：**没有找到任何论文做过这个分解**

用 arXiv API 做了 4 组检索（`abs:isotope AND abs:H3S`、`abs:"isotope coefficient" AND abs:pressure`、`abs:D3S OR "sulfur deuteride"`、`ti:"High-pressure phase diagram of hydrogen and deuterium sulfides"`），逐条查看返回条目的摘要：**没有任何一篇把测量的 α 拆成"结构相变贡献 + 配对贡献"，也没有任何一篇明确说"低压 α 被 R3m–Im-3m 相变的同位素移动污染"。**

⚠️ **这个"未找到"的强度有限**：WebSearch 配额用尽，只能用 arXiv API 的关键词匹配（不含全文检索、不含非 arXiv 期刊如 Angew. Chem./Nature Phys./PRB 的正文）。**这不是"证明不存在"**，只是"在可及范围内未命中"。新颖性主张可以提，但要写成"据我们所知"。

### 4.2 五个"最接近"的工作（必须在 P1b 里引用/区分）

**(a) Minkov et al. 2020 — 最接近的实验陈述**
Angew. Chem. Int. Ed. **59**, 18970–18974 (2020)，DOI 10.1002/anie.202007091（英文版）/ 10.1002/ange.202007091（德文版）；PMID 32633069；PMC7589447；**无 arXiv 预印本**（arXiv 全库搜 "sulfur deuteride" 返回 0 条）。
摘要要点（Crossref 与 Semantic Scholar 两路独立取回、互相印证）[二手，未拿到全文]：
- 从单质 S + 过量 H₂ 直接原位合成，避免"贫氢副产物污染"
- "**the Im3̄m phase of D3S samples exhibits a T_c significantly higher than previously reported values (ca. 150 K), reaching a maximum T_c of 166 K at 157 GPa**"
- **把 150 GPa 以下 Tc 的下降归因于晶格的"逐渐菱方畸变"（gradual rhombohedral distortion）**
- Cccm 相是含分子氢单元的绝缘体，非金属

→ 这是文献中**最接近"低压 Tc/α 被结构畸变污染"**的一手实验陈述，但它只讲 Tc，**没有讲 α，也没有做分解**。

**(b) Bianco et al. 2018 — 最接近的理论陈述，且直接质疑"kink = 相变"**
正文 p.6 逐字 [文本]：
> "In experiments the observed kink in the pressure dependence of T_c on pressure release, which could be ascribed to the Im3m → R3m transition, is estimated to happen at **150 GPa and 160 GPa for H3S and D3S**^{1,3}, thus at pressures **higher** than the critical pressure of the second-order Im3m→R3m transition obtained in our calculations"

正文 p.7 逐字 [文本]：
> "On the experimental side, one possible interpretation of our results is that **the observed kink in the T_c(P) curve is not related to a Im3m-R3m phase transition**. On lowering pressure the system in Im3m could undergo a very different, maybe even not simply displacive, phase transition, and different sulfur hydride stoichiometries could be involved. ... another, more intriguing, hypothesis is that non-hydrostatic experimental conditions have induced the rhombohedral transition at higher pressures."

⚠️ **溯源警告**：Bianco 的 "150 GPa and 160 GPa" 挂的引用是 ref 1 = Drozdov 2015、ref 3 = Einaga 2016。但我逐页读完这两篇后确认：**Drozdov 2015 自己的拉曼把 H3S 相变放在 ≈180 GPa（§1.4），Einaga 2016 把 Tc 骤降放在 123 GPa 且说 XRD 看不到结构变化（§3.3），两篇都没有出现 "160 GPa for D3S" 这个数字。** 因此 **"150/160 GPa" 是 Bianco 等人自己对那两篇 Tc(P) 图的读图判断，不是被引论文的原话** —— P1b 若用这个 10 GPa 的"实验 ΔP_c"，必须标注为 [二手/读图]，不能写成"实验测得"。

**(c) Taureau et al. 2024** — 明确把"局域偶极形成压力的同位素移动"与"Tc 极大值位置的相应变化"挂钩（摘要原话 "in agreement with the corresponding change in the T_c maximum location"）。**这是最接近 P1b 物理图像的一篇**，但他们做的是"结构量与 Tc 峰位的对应"，**没有写出 α 的分解式，也没有计算 α**。P1b 的新颖性主要就落在"写出 α_app = α_pair + α_struct 并给出式 (3) 的闭式预言"上。

**(d) Gor'kov & Kresin** — arXiv:1511.06926（Sci. Rep. **6**, 25608 (2016)）摘要 [文本]：
> "The isotope effect (change of T_c at the substitution of deuterium for hydrogen) owes its origin to the high frequency of phonons and **is different in different phases**."
> 并识别 "**P_cr = 123 GPa**" 处的一级相变（可能与 CDW 不稳定性相关）。
另有综述 arXiv:1802.02296（Rev. Mod. Phys. **90**, 011001 (2018)）指出 α 随压力变、且 150 GPa 附近有"一系列结构相变"，**但没有把两者联系起来** [二手，来自 arXiv API 摘要级]。
→ **"同位素效应在不同相里不同"这个定性观点已经被 Gor'kov-Kresin 提出过（2016）**，P1b 应引用并说明自己的增量是定量分解式。

**(e) Bianconi & Jarlborg — 竞争性解释，必须正面处理**
arXiv:1507.01093v7，摘要 [文本]：实验显示同位素系数随压力"**between 0.25 and 0.5**"变化，"in disagreement with Eliashberg theory"；他们归因于 **Lifshitz 相变（费米面拓扑变化）+ 零点晶格涨落**，并主张多能隙 BCS-BEC crossover / BPV 理论。
→ 这是对"α(P) 为何变化"的**另一套完整解释**，与结构相变无关。P1b 的证伪计划（§3）必须能把两者区分开——否则即使式 (3) 拟合上了，也不能排除 Lifshitz 机制。

**(f) 其他已排除的**：Szczęśniak arXiv:1601.01531（归因于正常态强重整化）、arXiv:1605.01259（多能带 + 氢键非谐）、Talantsev arXiv:2004.03155 与 2002.12859（Debye 温度/Tc 比值路线，结论是"起源仍待发现"）、arXiv:2607.18095（D3Se 的非谐压制 α，姊妹体系，可作旁证但非分解）。均**不构成**分解工作。

---

## 5. LaH10/LaD10：Fm-3m 稳定窗口的同位素差异

### 5.1 理论：Errea et al. 2020 只给出**单侧界**，没有给 ΔP_c

**arXiv:1907.11916v1**（本地 `library/pdf/D3/2020_Ion_Errea_1907.11916.pdf`），Nature **578**, 66–69 (2020)。

正文 p.4 逐字 [文本]：
> "While the Fm-3m phase of LaH10 **remains a minimum of E(R) as low as ~129 GPa**, the case of **LaD10 shows instabilities at 126 GPa**, implying that at this pressure the Fm-3m phase of LaD10 distorts to a new phase (as suggested by Drozdov et al. [3])."

Fig. 2 图注 [文本]："The case of deuterium develops an instability at low pressures (**126 GPa**) consistent with experimental evidence."

⚠️ **这不是一个 ΔP_c 测量**。Fig. 2 显示 LaH10 只算了 **129 / 163 / 214 / 264 GPa** 四个点，LaD10 只算了 **126 / 159 / 210 / 260 GPa** 四个点——**129 和 126 只是各自网格的最低压力点**。可推出的只有：
- P_c(LaH10) < 129 GPa（在 129 GPa 仍稳定）
- P_c(LaD10) > 126 GPa（在 126 GPa 已失稳）

**这两条甚至不能严格证明 ΔP_c > 0**（因为没有算 LaH10 @126 GPa）。定性方向与 H3S 一致（重同位素更早失稳），但**量级完全未定**。
**未找到**任何论文对 LaD10 做过 P_c 的压力扫描 / 二分搜索。

### 5.2 实验：Drozdov et al. 2019 自己承认结构未匹配

**arXiv:1812.01561**（本地 `library/pdf/P1/2019_A_P_Drozdov_1812.01561.pdf`），Nature **569**, 528 (2019)。正文 p.3 逐字 [文本]：

> "A complication arose from the existence of competing phases: in hydrides different phases with T_c =250 K (Fm3m structure) and 215 K (the structure to be determined) are formed under the same pressure-synthesis conditions. This assertion is valid also for the synthesized deuteride phases displaying T_c of ~168 K and ~140 K. ...
> **At present we cannot directly compare T_cs in hydrides and deuterides having the same crystallographic structures. We can naturally suppose that** the hydride with the maximum T_c (250 K) and the deuteride with maximum T_c (168 K) have the same crystal structure, that is the Fm3m phase ... Under this assumption, we obtained the isotope coefficient α of about **0.57** and **0.62** for the superconductivity in the **fcc** and **P4/nmm** phases, respectively."

**这是文献中最直白的"α 被结构问题污染"的作者自陈**，且是 P1b §4 需要的正是这一条。要点：
- α = 0.57 / 0.62 **建立在"结构相同"的假设上，作者明说无法直接比较**
- 两个配对**都不是同压比较**（250 K 与 168 K 各取自各自序列的最大值）
- 更糟的是**化学计量也不同**：Fig. 3 图注给出 Fm-3m 相精修为 **LaH9.6**（a = 5.1019(5) Å, V = 132.80(4) Å³），P4/nmm 相为 **LaD10.6**（a = 3.7258(6) Å, c = 5.0953(12) Å, V = 70.73(2) Å³）[文本]。即"氢化物"与"氘化物"样品的 H:La 比就不一样。
- LaH10 圆顶峰位 [文本]："after the initial increase and reaching the maximum T_c = **250-252 K at ~170 GPa**, T_c decreases abruptly at higher pressures (Fig. 1, inset)"。**LaD10 的峰位压力全文未给。**
- 作者对 α > 0.5 的解释 [文本]："The values of α>0.5 could be explained probably by a soft anharmonism of H/D vibrations (ω_H/ω_D < √2)"——即他们自己归因于非谐，不是结构。

### 5.3 补充：常被转引的 α = 0.46/0.47 依然溯源不明

Errea 2020 正文 p.4 [文本]："the isotope coefficient ... is close to 0.5 (**0.43** around 160 GPa) ... and it is **in agreement with the experimentally reported α of 0.46**."
但一手实验来源（Drozdov 2019 终版 + arXiv v1）报的是 **0.57/0.62**，无 0.46。此疑点在上一轮已记入 memory，本轮**未获新证据，维持"出处不明"**。

### 5.4 对 P1b §4 的判断

**LaH10 的 α 反常（0.57–0.62 > ½）确实很可能含结构成分**，但证据链与 H3S 不同：
- H3S：结构污染来自**同一相内的圆顶水平错位**（连续、可用式 (3) 处理）
- LaH10：污染来自**根本没在同一结构、同一化学计量、同一压力下比较**（离散、式 (3) 不适用）

因此 **P1b §4 目前的写法（"用式 (3) 解释 LaH10 反常"）在物理上不成立**，建议改为："LaH10 的 α 不满足式 (2) 的适用前提（同压、同相、同计量），应从数据集中剔除或单列，而不是用第五因子拟合。" 这与 P1 §4.2 已列的"LaD10 数据本身存疑"合流，且更强——因为是 Drozdov 自己写的。

---

## 6. 矛盾与不确定点（按对 P1b 的杀伤力排序）

**6.1 [最高] Minkov 2020 的 H3S 数据未知，直接决定 P1b 是否还有要解释的对象**
若 Minkov 2020 干净样品的 H3S 峰值仍是 ~203 K，则 α(157 GPa) = 0.29 ≈ α(P>170) = 0.30，**低压 α 异常消失，第五因子无用武之地**。若他们的 H3S 也被 boost，则异常可能仍在。**必须拿到全文。** Wiley `/doi/` 与 `/doi/pdfdirect/` 均 403，PMC 被 CAPTCHA 拦。可行路径见 §7。

**6.2 [高] H3S 结构相变压力，五个来源给出 82 / 91–103 / 123 / 134 / 180 GPa，无共识**

| 来源 | P_c(H3S) | 性质 |
|---|---|---|
| Taureau 2024 铁电（长程序）@200 K | ≈82 | PIMD 理论 |
| Bianco 2018 SSCHA-Hessian PBE | 91 | 理论 |
| Errea 2016 SSCHA PBE | 103 | 理论 |
| Bianco 2018 SSCHA-Hessian BLYP | 114 | 理论 |
| Einaga 2016 Tc 骤降点 | 123 | 实验（但 XRD 无对应信号）|
| Cherubini 2026 位移 QCP（T→0）| 134 ± 2 | 理论 |
| Taureau 2024 局域偶极形成 | 152–153 | 理论 |
| Bianco 2018 转述的"实验 kink" | 150 | 二手读图 |
| **Drozdov 2015 自己的拉曼** | **≈180** | **实验** |

**这直接决定式 (3) 里 d lnTc/dP 该在哪个压力取值**：若 P_c ≈ 91–103 GPa（Errea/Bianco），相变点远在圆顶上升沿以下，圆顶形状 D(u) 的"上升"根本不是相变驱动的，P1b §5 最后一条的物理归因失效（Cherubini 明确这么说）。若 P_c ≈ 152 GPa（Taureau 密度探针），相变正好在峰位，则 §1.3/§3.5 显示那里 d lnTc/dP ≈ 0，式 (3) 给 α_struct ≈ 0。**两种极端下第五因子都不产生 0.14 的贡献。**

**6.3 [高] 圆顶"错位 30 GPa"没有一手支撑**（§0-C）。可溯源的峰位错位 ≈ 0–10 GPa。

**6.4 [中] ΔP_c 的四个值分散在 12–35 GPa，且"哪个物理量才是对的"取决于机制**
若 Tc 圆顶由"局域偶极/局域矩形成"控制（Taureau/Cherubini 的图像），应取 ΔP_c ≈ 30–35 GPa；若由"热力学 R3m–Im3m 相变"控制，应取 12–18 GPa。**这两个选择差 2–3 倍，直接线性传到 α_struct。** P1b 必须把它作为显式敏感性参数，不能选一个值了事。

**6.5 [中] α_pair 本身在低压就更大，与"结构污染"竞争解释同一现象**
Errea 2016（同相同压、纯 Im-3m）给 α = 0.40 @155 GPa vs 0.35 @210 GPa [文本]。Villa-Cortés 2022 给的趋势相反（α_total 从 0.33@162 降到 0.26@250，也是低压更大）。**两个独立理论都说纯配对 α 在低压就更大**，这意味着 P1b §5 第 3 条"α_pair 在立方相区弱 P 依赖"的假设不成立，且**一部分（可能是全部）低压 α 上升已被配对物理解释掉**。

**6.6 [中] 与 Bianconi-Jarlborg 的 Lifshitz 解释无法区分**（§4.2e）。两套机制都预言 α 随 P 变、都能覆盖 0.25–0.5 的范围。P1b 需要一个**能区分二者**的判据（例如：Lifshitz 机制不预言 α 在两个同位素各自峰位相等，而错位机制预言相等）。

**6.7 [低] 术语陷阱**：Taureau 2024 的"transition pressure"在 Table I 里指**局域矩形成**，正文另有一个更低的**铁电相变压力** ≈82 GPa。转引时极易串。

---

## 7. 建议的后续检索方向

1. **[最高优先] 取 Minkov et al. 2020 全文**（DOI 10.1002/anie.202007091，CC BY-NC 开放获取）。已试失败的路径：`onlinelibrary.wiley.com/doi/...`（403）、`/doi/pdfdirect/...`（403）、`pmc.ncbi.nlm.nih.gov/articles/PMC7589447/`（CAPTCHA）。**建议由用户本地浏览器直接下载 PMC7589447 或 Wiley PDF 存入 `library/pdf/P1/`**，然后我用 Read 逐页读。需要抽取的是：他们的 **H3S 圆顶峰值与峰位**、D3S 完整 Tc(P) 点列、以及"gradual rhombohedral distortion"那段的原文与依据。这一篇决定 P1b 的存废。
2. **Bianco 2018 的附录 B（泛函/收敛测试）与 Fig. 3 的原始数据**：想要 ΔP_c 的不确定度而不只是中心值 16 GPa。本轮只读了 p.1–8，附录未读。
3. **Cherubini 2026 的 Supplemental Material**：是否有 D3S 的 PIMD？正文说只做 H3S，但 SM 可能有同位素测试。若有，会是 ΔP_c 的第五个独立值，且方法最新。
4. **Einaga 2016 发表版（Nat. Phys. 12, 835）与 arXiv v1 的差异**：arXiv v1 是 2015-09 投稿，Nature Physics 版 2016 年发表，图 3c 的数据点可能有增补。
5. **D3Se / H3Se 姊妹体系**（arXiv:2607.18095，2026）：若那里也有 ΔP_c 与 α，可作为"第五因子普适性"的独立检验，且不受 H3S 特有的样品污染争议影响。
6. **系统性反面检索**（本轮受配额限制未做）：用 Google Scholar / Web of Science 全文检索 `"isotope" AND "transition pressure" AND ("H3S" OR "sulfur hydride")` 在**非 arXiv 期刊**（Angew. Chem.、Nat. Phys.、PRB、J. Phys. Chem.）里的正文，以加强 §4.1 的"未找到"强度。

---

## 附录 A：本轮新增/更正的 arXiv ID（可直接入 `library/index.md`）

| 论文 | arXiv ID | 期刊 | 本地状态 |
|---|---|---|---|
| Errea et al., Quantum H-bond symmetrization | **1512.02933** | Nature 532, 81 (2016) | 已在 `library/pdf/P1/`（此前 ID 未核实，现已确认；**不是** 1607.02348）|
| Bianco, Errea, Calandra, Mauri, H & D sulfides phase diagram | **1802.07968** | PRB 97, 214101 (2018) | **未入库，建议下载** |
| Taureau, Cherubini, Morresi, Casula, Quantum symmetrization QMC+PIMD | **2307.15684** | npj Comput. Mater. 10, 56 (2024) | **未入库，建议下载** |
| Einaga et al., Crystal structure of 200 K SC phase | **1509.03156** | Nat. Phys. 12, 835 (2016) | **未入库，建议下载** |
| Cherubini, Raghav, Casula, Displacive QCP in H3S | 2602.00833 | (2026) | 已在 `library/pdf/D3/` |
| Minkov, Prakapenka, Greenberg, Eremets, D3S 166 K | **无 arXiv** | Angew. Chem. Int. Ed. 59, 18970 (2020) | 需人工下载（PMC7589447）|
| Gor'kov & Kresin, Pressure and high-Tc in sulfur hydrides | 1511.06926 | Sci. Rep. 6, 25608 (2016) | 未入库 |
| Bianconi & Jarlborg, Lifshitz transitions in sulfur hydride | 1507.01093 | (2015–2016) | 未入库 |

## 附录 B：明确"未找到"的条目（勿再重复检索）

- Drozdov 2015 的 (P, Tc) 数值表 —— **不存在**，只有图
- Drozdov 2015 的任何低压 α 值或 0.44 —— **不存在**
- Einaga 2016 实验测定的 R3m→Im-3m 转变压力（H3S 或 D3S）—— **不存在**（XRD 原理上分辨不了，作者自陈）
- Einaga 2016 中 173 GPa 以下的 D3S XRD —— **不存在**
- Cherubini 2026 的 D3S 数据 —— 正文中**不存在**（SM 未查，见 §7.3）
- LaD10 的 Fm-3m 临界压力扫描 —— **不存在**，Errea 2020 只给单侧界（126/129 GPa 是各自网格最低点）
- Drozdov 2019 中 LaD10 峰值 Tc 对应的压力 —— **不存在**
- Minkov 2020 的 arXiv 预印本 —— **不存在**（arXiv 全库搜 "sulfur deuteride" 返回 0 条）
- 任何把 α 分解为结构 + 配对两部分的论文 —— **在可及检索范围内未找到**（强度见 §4.1 的限定）
