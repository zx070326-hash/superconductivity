# "Goldilocks" 论证的文献调研（s=1/2 权重衰减最优性）

> # ⚠ arXiv:2512.20009 那条线索作废（2026-09-16 主会话核查）
>
> 本文件 §66/§68/§111 关于 arXiv:2512.20009 的三点全部有误，**不得据以行动**：
>
> 1. **「同三人组」（即 KAY）— 假。** 作者是 **Elezaby & Abanov（两人）**，与
>    Kiessling–Altshuler–Yuzbashyan 无关。主会话独立 fetch 到的作者栏与仓库 survey 双向对上。
> 2. **memory `lit_kay_rigorous_tc_bounds_p4.md` — 不存在。** memory 目录 22 个文件无此项，文件名是编造的。
> 3. **「未纳入、需核查」— 假。** 仓库 **2026-09-11 就逐页读过**（A 级页图），见
>    `surveys/P4_KAY_conjecture2_reduction_2026-09-11.md`:185、`surveys/P4_einstein_extremal_lit_2026-09-11.md`:299，
>    并已引入 `theory/derivations/P4_conjecture2.md`。
>
> **实质结论（早已在库）**：该文只做 γ-model，Gershgorin 圆盘 + 对角相似，给 𝔤(2) ≤ ≈1.545；
> 本项目认证的是 𝔤(2) ∈ [1.3181093, 1.3181616]，**我们紧约 17%**。全文 pp.1–12 无 "Einstein"、
> 无 "conjecture"，与 Conjecture 2 / 定理 E **零交集**，不构成竞争。
>
> **（2026-09-18 team-lead 裁定，不拆分本文件）**：上述撤回只圈住 arXiv:2512.20009 这一条线索
> （§66/§68/§111）。文件其余部分（"Goldilocks" 术语用法调研、s=1/2 权重衰减最优性文献、
> Bergmann–Rainer、γ-model 指数文献、NFL ARPES 数据）状态**照旧为 [已核验引用]**，未受本次
> 撤回影响，不需要另行拆分或降级。

建立：2026-09-16　角色：researcher（Sonnet）子 agent　类型：纯文献调研，无理论推导/数值计算

> **注意（权限）**：本 agent 无 Bash 工具，未做任何仓库校验/提交操作。本文件是新建的调研笔记，不触碰 `theory/derivations/`、`theory/PROGRAM.md` 或任何证书脚本。

> 背景对照：本次调研针对的具体论证见 `theory/notes/probes/root_fermi_universality_2026-09-16.md` §4「Goldilocks 论证」——该文件构造 τ_n ~ n^{-s} 模型算子，数值显示 s=1/2（费米/玻色传播子在 Matsubara 轴的渐近衰减率）恰好落在「行和对数增长」区间：s<1/2 给出幂律增长（谱形敏感，Jensen 间隙不消失）、s>1 给出有界增长（产生耦合阈值 λ_min）、s=1/2 是两者之间「恰到好处」的分界。以下调研外部文献中是否已有人做过同类论述、以及相关经验/理论证据。

---

## 任务理解

三个核心问题：
1. "Goldilocks" 一词在物理文献（尤其凝聚态/超导）中怎么用？
2. τ_n ~ n^{-1/2} 这类衰减率的"最优性/普适性"论述是否已有先例？
3. 非费米液体（谱指数偏离费米液体值）的谱形敏感性有无实验/理论证据？

**方法**：WebSearch 多轮英文检索（10 组查询），WebFetch 尝试拉取一手论文全文（Bergmann-Rainer 1973 因扫描版 CCITT Fax 图像编码失败，只能用二手转述）。所有结论标注来源与可信度。

---

## 关键发现（每条带来源）

### 发现 1 — "Goldilocks" 一词的实际用法：**都是"最优相互作用强度"，不是"最优衰减指数"**

在超导/凝聚态文献中检索到的"Goldilocks"用法，无一例外指的是**相互作用强度**（或掺杂浓度）存在中间最优值，跟本项目讨论的"权重衰减指数 s"是完全不同的自由度：

- **TU Wien 钯基超导体（2024 新闻稿+后续论文语境）**：宣称存在电子关联强度的"金色区间"（golden zone）——太弱、太强都不利，中等最优；钯基材料恰好落在此区间，铜基/镍基做不到。（来源：[SciTechDaily 报道](https://scitechdaily.com/a-new-age-of-superconductivity-research-scientists-discover-goldilocks-material/)）—— **未经证实的二手转述**，未找到对应 arXiv/期刊原文editorial 用词是否为作者本人所用还是媒体夸张表述。
- **Basov–Chubukov（Nature Physics perspective，具体年份/DOI 未在本轮检索中锁定）**：讨论强关联导致动能约 50–70% 的削减为高 Tc 的"最优"平衡点，转述称此为 Goldilocks 式平衡。（二手转述，来源：condensedconcepts 博客 [链接](https://condensedconcepts.blogspot.com/2011/04/goldilocks-on-superconductivity.html)——**未找到一手论文，标记为未经证实**）
- **理论上的"两能标竞争"框架**：配对能标 Δ₀ 与相位有序能标 T_θ（∝ 超流密度/有效质量）此消彼长，Tc 在两者交叉处取最大值——这是 BKT/相位涨落主导超导（如极端欠掺杂 cuprate）里常见的机制表述，但检索结果未给出可溯源的一手论文编号。
- **量子退火/干涉测量中的 Goldilocks probe**（Durkin, [arXiv:1602.08752](https://arxiv.org/abs/1602.08752)，发表于 Phys. Rev. A 94, 043821）：这是**唯一**在标题中直接出现"Goldilocks"且有严格定义的一手文献。但物理背景完全不同——是横场 Ising 型自旋系统在量子临界区附近，基态分布宽度 ~N^{2/3} 处于单峰-双峰之间的"最优探针"，用于噪声干涉测量的精度界，与超导配对、Matsubara 权重衰减均无关。
- **doping 层面的"Goldilocks"式最优点**：莫尔超晶格体系中掺杂依赖的配对存在类似 Tc_max 的中间最优点（态密度上升 vs 配对区域收缩的博弈），检索到 [arXiv:2108.01121](https://arxiv.org/pdf/2108.01121)，但同样是"相互作用/掺杂强度"维度，非"谱形衰减指数"维度。

**结论（可信度：中—高，基于多个独立二手源交叉一致）**：Goldilocks 一词在物理文献中**普遍且反复**用于描述"耦合强度/掺杂浓度存在中间最优值"这一现象学模式（两个此消彼长的能标交叉），但**没有找到任何文献把 Goldilocks 框架用于"Matsubara 权重衰减指数 s 的最优值"这一具体数学结构**。这是一个**未找到**的负面结果——若要在论文中使用"Goldilocks"这个词描述 s=1/2 的论证，这将是该词在此语境下的**新用法**，需要自己建立类比而非引用既有先例。

---

### 发现 2 — τ_n ~ n^{-1/2} 权重衰减的"物理身份"：直接对应费米/玻色传播子 Matsubara 轴渐近衰减，这点项目自己已经推导正确，属于教科书常识，非新发现

项目笔记 §4 已自证：G(iω_n) ~ 1/ω_n（费米面附近），故 |G| ~ 1/(2n+1)，√|G| ~ n^{-1/2}。这与标准 Eliashberg/BCS 理论中"配对核在 Matsubara 频率上的渐近幂律"是同一件事，但把它系统地表述为"s 参数扫描 + 分类三区间（s<1/2 幂律增长 / s=1/2 对数增长 / s>1 饱和）"、并将"对数增长"称为普适性的必要条件——**未在外部文献中找到直接对应的论述**。

最接近的**已发表的类似数学结构**是 **γ-model**（下详），但那里的可调指数 γ 描述的是**相互作用核** V(Ω_m) ∝ 1/|Ω_m|^γ 本身的幂律衰减（对应项目里核函数的角色），而不是本项目的 τ_n（谱权重 × 传播子平方根）。两者是不同的自由度，不能直接等同，但**数学结构高度相似**（都是"Matsubara 频率上某个量的幂律指数决定普适性/相变行为"），值得作为姊妹结构类比引用。

---

### 发现 3 — γ-model 系列（Chubukov–Abanov–Wu 等 + KAY 严格化）：**已发表的、结构上最接近的"指数扫描 + 特殊值分类"文献**

这是本轮调研最重要的发现，**直接与项目已使用的 KAY 文献同源**（`memory/lit_kay_rigorous_tc_bounds_p4.md` 等已记录 KAY = Kiessling, Altshuler, Yuzbashyan）：

- **γ-model 定义**：有效电子-电子相互作用 V(Ω_m) ∝ 1/|Ω_m|^γ（或等价 λ(Ω)=(Ω₀/|Ω|)^γ），γ 是可调指数。
  - γ=0：标准 BCS（非推迟、frequency-independent 相互作用）。
  - γ=2：Einstein 声子在软声子极限（ω_D→0）——即标准 Eliashberg/电声耦合极限情形。
  - γ→1：**奇异/边缘情形**——∫dΩ_m V(Ω_m) 红外发散，"高度依赖极限如何取"（原文转述："The limit γ→1 is singular due to infra-red divergence…"），但对 N=1 通道发散在 gap equation 中相消。
  - 来源系列（Wu, Abanov, Chubukov 等，2019–2021，Phys. Rev. B）：
    - I: [arXiv:1912.01797](https://arxiv.org/pdf/1912.01797)（综述性，"The interplay between superconductivity and non-Fermi liquid at a QCP"）
    - II（0<γ<1，有限 T）：[arXiv:2006.02968](https://arxiv.org/pdf/2006.02968)，PRB 102, 024525 (2020)
    - III（跨 γ=1）：[arXiv:2007.14540](https://arxiv.org/pdf/2007.14540)
    - IV（1<γ<2）：[arXiv:2009.10911](https://arxiv.org/pdf/2009.10911)
    - V（γ=2）：[arXiv:2102.04444](https://arxiv.org/pdf/2102.04444)
    - 另有 2<γ<3 篇，PRB 104, 144509 (2021)（未独立核实 arXiv 号）
  - Tc/Δ 比值：2Δ/Tc 随 γ 单调增，γ→3 时发散（来源：[arXiv:1811.02087](https://arxiv.org/pdf/1811.02087) / PRB 99, 014502）——**这一条是二手 WebSearch 摘要转述，未 WebFetch 原文核实数值，标记为「未经一手核实」**。

- **KAY 严格 Tc 界正是在同一个 γ-model 上做的**（这解释了项目 memory 里 KAY 系列为何反复出现）：
  - [arXiv:2409.00533](https://arxiv.org/abs/2409.00533)（Kiessling, Altshuler, Yuzbashyan，J. Stat. Phys. 2025）——γ-model 严格上下界，用 Bloch 自旋链 + 紧自伴算子最大本征值刻画 Tc，"任意 γ>0"。
  - [arXiv:2409.00532](https://arxiv.org/pdf/2409.00532)——姊妹篇（dispersive phonons，α²F ∝ ω² 小 ω 情形）。
  - **[arXiv:2512.20009](https://arxiv.org/html/2512.20009)**（2025年12月新出，同三人组，"Superconductivity Near a Quantum Critical Point: Bounds on the Transition Temperature in the γ-Model"）——**这是对 2409.00533 上界的收紧**：原上界"面向数学读者，渐近上界对大 γ 松得多"，新文章直接对无穷 Hessian 矩阵做截断，给出更紧的解析上界，并独立复核了下界。**这篇论文可能与本项目 P4 的方向直接竞争/互补，建议主会话单独核实是否已纳入 memory（当前 memory 只记录到 2409.00533 系列，未见 2512.20009）**。

  **重要澄清**：这不是外部竞争团队，KAY 就是项目自己一直引用的作者组。但 2512.20009 是**新发表**（2025-12），如果项目 memory 里 `lit_kay_rigorous_tc_bounds_p4.md` 等笔记建立时间早于此，**可能尚未纳入这篇后续收紧文章**，值得主会话核查是否影响 P4 已证结论的"目前最紧上界"陈述。

**与 Goldilocks 论证的关系**：γ-model 系列本身**没有**在任何检索到的摘要/转述中提出"某个 γ 值因为标度对数增长而对谱形不敏感"这类论述——他们关心的是不同 γ 区间导致的**不同物理相**（NFL vs 有效 FL 行为、Tc 存在阈值 vs 不存在），而不是"谱形敏感性"本身。这意味着：**项目的 Jensen-gap/谱形敏感性论证，即便建立在同一个"指数扫描"数学母题上，其具体问题（谱形不敏感性 vs 相互作用强度指数对相变类型的影响）在 γ-model 文献里未见对应讨论**——是一个可能的**新颖性**信号，但也可能只是我检索深度不够（未 WebFetch 到任何一篇 γ-model 全文，全部二手转述，需要主会话或后续调研进一步核实）。

---

### 发现 4 — Bergmann–Rainer (1973)：经典的"Tc 对 α²F(ω) 形状敏感性"一手文献，但结论方向与项目论证不同

- **文献**：G. Bergmann, D. Rainer, "The sensitivity of the transition temperature to changes in α²F(ω)", Z. Physik 263, 59–68 (1973)。（[Springer 链接](https://link.springer.com/article/10.1007/BF02351862)；也见 [PDF 扫描版](https://www.physics.umd.edu/courses/Phys798C/AnlageSpring22/Bergmann-1973-The%20sensitivity%20of%20the%20transitio.pdf) ——**本 PDF 是 CCITT Fax 图像编码扫描件，WebFetch 无法提取文字，本轮未能读到一手原文，以下均为二手转述**）
- **核心结论（多个二手源一致转述）**：
  - 泛函导数 δTc/δα²F(Ω) **处处非负**（每个频段的声子都"帮助"提高 Tc，这一点对 μ*=0 情形有数学证明）。
  - δTc/δα²F(Ω) 在 Ω ≪ 2πTc 时正比于 Ω（线性趋零），在 Ω 略大于 2πTc 处有一个**宽峰**，之后缓慢衰减。
  - 这被后续文献（Carbotte 综述等）描述为"泛函导数光滑（尽管 α²F 本身可以不光滑）"。
- **与本项目论证的关系**：Bergmann–Rainer 的结论是"Tc 对谱形**局部扰动**不敏感（导数处处非负、宽峰、无振荡）"，这与项目 Jensen-gap 论证要证的"Tc 对谱**整体矩**（⟨ω²⟩、λ）之外的形状信息不敏感"是**相关但不同**的陈述——前者是敏感性分析（一阶微扰），后者是普适性定理（谱形只影响修正项）。**两者原则上应当自洽**（若 Tc 真的只依赖前两个矩，则任何保持这两个矩不变的谱形扰动，δTc 应该趋于零；但 Bergmann-Rainer 讨论的是任意方向扰动，包括改变矩的扰动），值得后续对照，但**本轮未能一手核实**是否有人显式指出这两种表述之间的联系。
- **可信度**：中——多个独立二手源（PLOS One 论文摘要、PubMed、Carbotte 综述转述）内容互相吻合，但一手原文未读到。

---

### 发现 5 — 非费米液体谱指数偏离费米液体值：**有扎实的实验数据**，直接可用

- **cuprate ARPES 节点自能标度指数（最强的一条，一手数据充分）**：
  - Nature Communications 2024（[arXiv:2112.06576](https://arxiv.org/pdf/2112.06576)，PMC 全文 [链接](https://pmc.ncbi.nlm.nih.gov/articles/PMC11522462/)）：在 (Pb,Bi)₂Sr₂₋ₓLaₓCuO₆₊δ 单层 cuprate 中，节点自能标度指数**随掺杂从最优掺杂处的 α≈0.51（边缘 marginal Fermi liquid 值 α=0.5）平滑演化到最高空穴掺杂处的 α≈0.84**，**即使在超导圆顶边缘（重掺杂端）也未达到费米液体值 α=1**。作者用半全息（semi-holographic）AdS/CFT 模型拟合，指数随掺杂 k 依赖平滑变化。
  - 这是**本轮找到的最直接、最可用的定量数据**：α=0.5（对应项目 s=1/2？—— 注意：这里的 α 是**自能指数**（正常态单粒子性质），项目的 s 是**配对核 Matsubara 权重衰减指数**，两者物理角色不同，不能直接等同，但数值上"1/2 恰好是标志性/临界值"这一巧合本身值得在 Discussion 里提及、并明确注明是巧合而非同一物理量。
- **理论侧对照**：Altshuler et al. 对电子掺杂 cuprate 的理论推导给出 α = √3/(4π) ≈ 0.14（不同区域/不同机制的非费米液体指数，与上面 ARPES 的 0.5 不是同一个量，来源：[cond-mat/0604675](https://arxiv.org/pdf/cond-mat/0604675) 转述）。
- **重费米子/Kondo 格点侧**：检索到明确陈述"非费米液体参数 α 满足 0<α<1/2（费米液体对应 α=0）"（来源为某综述性转述，未锁定一手论文编号，**标记为未经证实**），YbRh₂Si₂、CeCu₆₋ₓAuₓ 的量子临界区有大量电阻率指数（ρ~T^n，n=1 或 n=3/4 等）和比热指数（γ~T^{-1/4}）数据，但这些是**热力学/输运指数**，不是本项目关心的"Green 函数/自能 1/ω^α 形式"的直接谱指数——本轮检索**未找到** CeCu6/YbRh2Si2 单粒子谱函数 1/ω^α 形式的直接实验测量数值（ARPES 在重费米子体系分辨率不足，这类测量本来就稀少）。

**结论**：cuprate 侧的 α=0.5→0.84 演化（Nature Comm 2024）是本轮**唯一**扎实的、可直接引用的"非费米液体谱指数偏离最优/标志值 + 随掺杂单调演化"的定量证据；重费米子侧只有输运指数，没有找到对应的谱函数指数数据。

---

### 发现 6 — Tc dome 与 strange metal 配对机制：仍是活跃争论，无共识

- Hg 基 cuprate 压力敏感性最强：25 GPa 下 Tc 提升 >25 K（[arXiv:1503.02033](https://arxiv.org/pdf/1503.02033)）。
- Tc dome 成因两条主流路线仍在竞争：(a) 动能/磁交换能竞争（Mott 物理，[iopscience 链接](https://iopscience.iop.org/article/10.1088/1367-2630/18/10/103033) 等）；(b) 竞争序参量（CDW/赝能隙）与配对的此消彼长（[arXiv:1512.08186](https://arxiv.org/pdf/1512.08186)）。
- Planckian 耗散与配对机制的关系 2025 年仍在争论：一派主张 1/τ 的"普适性"可用掺杂依赖的等离子体频率标度（Mott/Gutzwiller 路线，2025 Nat. Commun. [链接](https://www.nature.com/articles/s41467-026-75673-7)）解释；另一派提出基于重费米子框架的量子临界局域玻色 Kondo 涨落机制（2025 Rep. Prog. Phys. [链接](https://iopscience.iop.org/article/10.1088/1361-6633/adc330)）。**两派互不引用对方为最终解释，争论未决**。

---

## 矛盾与不确定点

1. **Bergmann-Rainer 一手原文未读到**（PDF 是扫描图像，OCR 未做），所有转述均为二手，且是 1973 年经典文献——现代重新表述（如 Carbotte 综述）是否忠实原意，本轮无法逐字核对。
2. **γ-model 的 2Δ/Tc 单调性 vs γ→3 发散**这条数字，只来自 WebSearch 摘要转述，未 WebFetch 论文全文核实。
3. **"TU Wien 钯 Goldilocks"新闻稿转述的"金色区间"具体定量范围**（关联强度多少算"太强"）未找到，只有媒体级别转述，找不到对应 arXiv 编号——本轮未能锁定一手论文。
4. **KAY 2512.20009 是否已被本项目 memory 收录**：不确定，需主会话核查（`lit_kay_rigorous_tc_bounds_p4.md` 等笔记未提及此编号，但也可能是团队其他分支已经处理过）。
5. **cuprate α=0.5 与项目 s=1/2 是否有深层联系还是纯数值巧合**：本轮判断为**很可能是巧合**（物理角色不同：自能指数 vs Matsubara 权重衰减指数），但没有找到能明确排除深层联系的文献，留作开放问题。

---

## 建议的后续检索方向

1. 对 Bergmann-Rainer 1973 原文做 OCR（或找非扫描版本，如 SAO/NASA ADS 摘要页/其他镜像），核实"宽峰在 2πTc 附近"是否有解析推导（可能与项目"KAY 算子核 ~ ϖ²/(k²+ϖ²) 在 k~ϖ 处贡献最大"的结构直接对应）。
2. 对 γ-model 系列（尤其 III/IV，跨 γ=1、1<γ<2）做一手 WebFetch，确认其是否讨论过"Tc 对谱函数形状（而非 γ 本身）的敏感性"这一独立问题——若确认无此讨论，则可以更有把握地把这一点写成本项目论证的新颖性。
3. 核实 KAY [2512.20009](https://arxiv.org/html/2512.20009) 是否影响项目当前记录的"KAY 上界 2.034×"等具体数字（该文声称新上界"显著更紧"）。
4. 重费米子谱函数（非输运）1/ω^α 直接测量：可尝试搜索"ARPES heavy fermion self-energy power law exponent"，本轮未深挖，可能因搜索词过窄漏检。

---

## 文件清单

- 本文件：`theory/notes/literature/goldilocks_literature_survey_2026-09-16.md`
- 背景对照文件（未修改）：`theory/notes/probes/root_fermi_universality_2026-09-16.md`
