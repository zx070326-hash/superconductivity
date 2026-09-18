# 自旋涨落媒介耦合矩阵：开源数据普查（不做代入计算）

> 2026-09-19 · worktree `spinfluc-quant`（独立分支，UNMERGED，等团队统一归档）
> 派活口径见 teammate 消息（team-lead）：为 `probes/spin_fluctuation_fese_material_probe_2026-09-19.md` 的
> 「10 对 patch 耦合只公开 5 对，另 5 对锁在 APS SM（403）」这一具体缺口，找**开源、免费、可直接取数值**的替代/补充数据集
> ——完整 patch/orbital 间自旋涨落耦合矩阵（或等价：RPA 磁化率 χ(q,ω) 全张量），覆盖任意材料。
> **本文性质：数据普查，不含任何代入计算。** 结论级别全部是「文献/仓库检索事实」，非 analytic/certificate/numerical-evidence。

---

## 0. 一句话结论

**没有找到任何"现成、完整、可直接下载数值"的自旋涨落耦合矩阵数据集**——无论是 Materials Project / AFLOW / NOMAD
这类通用材料数据库，还是论文自带的 Zenodo / figshare 补充数据。找到的最接近的东西是**三条开源代码工具链**，
它们能够**从头计算**出这样的矩阵，但仓库本身都不附带某个真实铁基/自旋涨落材料的完整现成输出——即"有工具，没有现成数据"。
此外还有一个**高通量 Wannier tight-binding Hamiltonian 库**，若覆盖目标材料，可作为上述代码的输入，
省去自己做 DFT+Wannierization 这一步（但仍需要跑 RPA/FLEX 计算本身）。

**对当前任务的直接含义**：`spin_fluctuation_fese_material_probe_2026-09-19.md` §4 指出的缺口（Essenberger FeSe
剩余 5/10 条边、APS SM 403）**无法靠"换一个已发表的现成数据集"绕过**——检索范围内不存在能替代它的公开数据。
唯一可行的路径变成"用开源代码自己算一遍"，这是一个新的、更大的工作量决定，超出本轮"数据普查"的授权范围。

---

## 1. 检索方向 1：开源材料数据库（Materials Project / AFLOW / NOMAD）

**结果：无对口条目。**

- Materials Project / NOMAD 均为通用 DFT 数据库（结构、能带、DOS、磁矩等静态量），**不存放** RPA 动态磁化率
  χ(q,ω) 或自旋涨落耦合矩阵这类"响应函数/微扰计算"级别的中间量。NOMAD 允许存原始计算归档，理论上有人可能上传过，
  但未检索到具体条目。
- 相邻的专门磁性数据库（NEMAD——东北大学磁性材料库；"Experimental Exchange Interaction Dataset"，*Sci. Data* 2025）
  存的是**实验测得的交换耦合常数**（中子散射拟合自旋波理论得到的 J_ij），属于**局域矩模型的有效交换积分**，
  不是巡游电子配对问题需要的 **k 空间/patch 分辨的配对核** ——概念上是两个不同层级的量，不能互相替代。
  （这条边界本身值得记：如果将来有人想"偷懒"用 INS 交换积分数据库冒充配对核数据，是**范畴错误**，不能这样代入。）

**判断：这个方向确认为空，不建议继续在通用数据库里搜。**

---

## 2. 检索方向 2：论文自带开源托管的补充材料（Zenodo/figshare/GitHub）

**结果：搜索到的 2023–2026 年自旋涨落配对论文（HgBa₂CuO₄ LDA+FLEX arXiv:2311.10030、ScAu₂Al ab-initio
arXiv:2311.06075、FeTe 五轨道 arXiv:2608.23467 等）均未确认有独立托管的补充数据文件**——WebSearch 未返回
Zenodo/figshare 记录，WebFetch 对 arXiv 摘要页/HTML 页也未发现补充数据链接（注：**只查了摘要页和搜索引擎索引，
未逐篇打开 arXiv ancillary files 列表页核实**，因此这条是"未查到"而非"确认不存在"，见 §5 待办）。
- Essenberger 系列本身（arXiv:1409.7968 / 1411.2121 / 1501.02200）均为 2014–2016 年论文，早于"补充数据上
  Zenodo/GitHub 而非期刊 SM"这个近年趋势变得普遍之前，与仓库已有笔记判定一致（APS SM 是唯一存放位置，403）。

**判断：这条趋势（新论文倾向开源托管补充数据）目前在"自旋涨落配对"这个具体子领域里没有找到实例。**
不代表以后不会有，但当前检索预算内一个都没找到。

---

## 3. 检索方向 3：专做自旋涨落理论的研究组公开代码

这是本次普查**唯一有实质产出**的方向。三个候选，均为**代码工具**而非数据集：

### 3.1 `Orth-Research/RPA_SrCo2As2`（GitHub，MIT 协议）

- 出处：Nedić, Christensen, Orth 论文（SrCo₂As₂ 的磁性相图）。
- 能力：**完整的 16 轨道 Hubbard-Hund 模型 RPA 自旋磁化率张量**代码，两个脚本
  `multiorbital_bare_susceptibility.py`（裸磁化率，可用 Wannier90 tight-binding 或对称化哈密顿量）
  与 `multiorbital_RPA_susceptibility.py`（RPA 磁化率张量 + 由最大本征值给出的相图）。
  形式上这**正是**任务要的"完整 patch/orbital 间耦合矩阵"层级的计算（q 依赖 + 频率依赖 + 轨道分辨）。
- **数据现状**：README 明确写"部分生成数据在 `results/` 目录，完整数据（画图用）因存储量大需按需索取"——
  ⟹ **仓库本身不含现成的完整 χ(q,ω) 数值文件**，需要跑代码或联系作者要数据。
- 材料是 SrCo₂As₂（阻挫磁性巡游磁体），不是铁基超导体，但代码架构与材料无关，理论上换一个
  Wannier90 `_hr.dat` 输入即可套用于 FeSe/铁基（前提：拿到 FeSe 的 Wannier90 哈密顿量，见 §3.4）。

### 3.2 `nikwitt/FLEX_IR`（GitHub + Zenodo DOI 10.5281/zenodo.10014333）

- 出处：Witt, van Loon, Nomoto, Arita, Wehling，*Phys. Rev. B* 103, 205148 (2021)。
- 能力：IR-basis 下的多轨道 FLEX（fluctuation exchange）近似，含 `multi_orbital/` 子目录，
  理论上可给出多轨道自旋磁化率 χ_s(q,iω) 全矩阵与配对本征值/本征函数——**FLEX 比 RPA 更进一步**，
  自洽包含自能修正，是文献里做"真实材料自旋涨落配对"更常用的框架（已发表应用于 Na_xCoO2·yH2O）。
- **数据现状**：GitHub 页面本身自陈"poorly commented/documented"，landing page 未暴露文件级细节；
  **未逐层打开子目录核实**是否含现成材料 example 输出（WebFetch 对该仓库两次尝试，一次网络中断、
  一次只拿到顶层目录名，未能深入 `multi_orbital/` 内部——这是本文唯一的检索缺口，见 §5）。

### 3.3 `johndeisz/flex-general`（GitHub，Fortran90）

- 多带 tight-binding 模型的 FLEX 近似实现。**仅搜到仓库存在，未做任何进一步核实**（时间预算约束下的取舍，
  优先深挖了 Python 生态的前两个候选）。

### 3.4 高通量 Wannier tight-binding Hamiltonian 库（可作为上述代码的输入）

- Garrity & Choudhary, *Scientific Data* 8, 106 (2021)，`github.com/usnistgov/jarvis`，
  1771 材料（1406 三维 + 365 二维）的 Wannier TB 哈密顿量数据库，配套网页应用 `ctcms.nist.gov/jarviswtb`。
- **⚠ 未确认是否覆盖 FeSe 或任何铁基超导体**——论文摘要明确该库的材料筛选标准是"强自旋轨道耦合、
  可能是拓扑绝缘体/半金属"，与"铁基超导体的自旋涨落配对"筛选标准不重合，**大概率不覆盖**，
  但未逐条查库确认（若团队后续想走这条路，第一步应是查该库的材料列表而非假设它没有）。

**判断**：三条代码路线里，`RPA_SrCo2As2` 的能力边界最清楚（明确"完整数据未附带，按需索取"），
是若决定"自己算"时**最值得优先尝试**的起点（代码开源、许可宽松、有已发表论文对照可验证复现）；
`FLEX_IR` 理论上更贴近真实材料计算但文档缺失，需要更多前期投入读代码本身。

**Essenberger 系列 TDDFT 自旋涨落配对理论本身没有找到开源实现**——检索未发现它被整合进 FLEUR/SPEX
（这两者本身开源，MIT 协议）等公开代码库，也没找到作者自己发布的独立实现。这意味着即便想"复现"
Essenberger FeSe 计算的完整 10×10 矩阵，也没有对应的开源代码可用；能用的开源代码（RPA/FLEX）
是**不同的理论框架**（静态/自洽 RPA、FLEX vs Essenberger 的 TDDFT 交换关联核），
不是同一套方法的开源版本，数值不会逐位对应。

---

## 4. 是否有产出可直接用于仓库 `_chk_spm_*.py` 框架？

**没有。** 现有 `_chk_spm_*.py` / `_chk_fese_essenberger_sf_bipartite_2026-09-18.py` 框架的输入是一个
**5×5（或更大）的标量耦合矩阵**（数值直接来自论文图说读数）。本次普查找到的三个代码工具链
**没有一个能直接产出这种"喂进去就能跑"的现成矩阵文件**——它们是需要自己配置输入（Wannier90 哈密顿量、
Hubbard U/J 参数）、跑计算（多轨道 RPA/FLEX 通常是数小时到数天量级的数值计算，需要 k 网格、Matsubara 求和收敛测试）
才能产出矩阵的**软件**，不是数据。因此按 team-lead 派活要求（"重心是数据普查，不做完整代入计算"），
**本文不尝试用这些代码跑出任何数字**，如实报告"找到工具、未产出数据"。

---

## 5. 未做 / 检索缺口（诚实边界）

1. **未逐篇打开候选论文（arXiv:2311.10030、2311.06075、2608.23467 等）的 ancillary files 页面**核实是否真的没有
   补充数据——只做了摘要页/HTML 全文检索，这是"未找到"而非"确认没有"的区别，按 CLAUDE.md 纪律不越级陈述。
2. **`nikwitt/FLEX_IR` 的 `multi_orbital/` 子目录内容未深入核实**（两次 WebFetch 一次网络中断，未重试第三次）——
   不确定其中是否含某个真实材料（如 Na_xCoO2）的现成 example 输出文件，这条若要继续追查是最短的下一步。
3. **`johndeisz/flex-general` 完全未核实**，只确认仓库存在。
4. **`usnistgov/jarvis` WTBH 库是否覆盖 FeSe/铁基超导体未查**——只读了论文摘要对材料筛选标准的描述，
   未实际查询该库的材料列表。
5. **未检查 Lischner et al. PRB 91,020502(R) (2015) 的标量 λ[s±]=0.30/0.28 数据来源是否可从 arXiv 获取**——
   这是 `spin_fluctuation_fese_material_probe` 已记录的"次要数据点"，本文未重新核实其可访问性，
   仍按原笔记标注为"已投影成标量、可用但丢失 patch 结构"处理。

---

## 6. 建议的下一步（若团队决定继续投入）

按工作量从小到大：

1. **（最小）** 深挖 `nikwitt/FLEX_IR` 的 `multi_orbital/` 目录 + 关联 PRB 103,205148 正文，
   确认它是否已含某个真实材料的可直接读取输出（哪怕不是 FeSe，只要是"有完整矩阵"的任意材料，
   即可满足 team-lead 原始派活里"不必是 FeSe，任意材料即可"这条放宽条件）。
2. **（中）** 查 `usnistgov/jarvis` WTBH 库材料列表是否含 FeSe 或其他铁基超导体的 Wannier90 哈密顿量；
   若有，配合 `RPA_SrCo2As2` 代码（换掉 `_hr.dat` 输入）可以**自己产出**一个 FeSe 的完整 RPA 磁化率张量——
   但这是"生成数据"而非"数据普查"，工作量与新颖性风险需要team-lead 重新评估。
3. **（最大）** 直接联系 Essenberger/Sanna/Gross 或后续跟进作者（PRB 94,014503 通讯作者）索取完整 10×10
   patch 矩阵原始数据——`spin_fluctuation_fese_material_probe_2026-09-19.md` §7 已建议此路，
   本文没有新增信息，只是重申它仍是"结构性缺口，非原则性缺口"里成本最低的一条。

---

相关笔记：
[`spin_fluctuation_fese_material_probe_2026-09-19.md`](spin_fluctuation_fese_material_probe_2026-09-19.md)（本文回应的缺口来源）·
[`spin_fluctuation_rigorous_tc_scoping_2026-09-18.md`](spin_fluctuation_rigorous_tc_scoping_2026-09-18.md)（立项级侦察，§1.1 右端物理量的循环性约束）·
[`../gamma/spin_fluctuation_structural_compatibility_2026-09-18.md`](../gamma/spin_fluctuation_structural_compatibility_2026-09-18.md)（结构相容性判定）·
[`../literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md`](../literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md)（Essenberger 数据 L1 转录来源）
