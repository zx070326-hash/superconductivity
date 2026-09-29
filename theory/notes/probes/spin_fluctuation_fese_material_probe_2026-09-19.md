# 自旋涨落配对代入 PF/Δ通道框架：FeSe 材料级数值探针

> ⚠ **2026-09-29 合入核查注**（`worktree-spinfluc-quant` 合入 master 时由核查 agent 复跑 + 主会话抽查；**正文未改，本注为唯一改动**；标「agent」者为核查 agent 报告、主会话未独立复算）：
> - ⛔ **头条结论「前提 1 在 FeSe 上不成立、spm_gap 定理不适用」（§0.1、§3 Part C、§0.2）站不住——勿引。** 主会话 2026-09-29 复核（exact 算术）：脚本 Part C 把 Λ^Δ 构造成「同族 −V、跨族 +V」再乘 D=diag(+Γ,−M)，**跨族边的相位符号被算了两次**。按 master `vertex/spm_gap_needs_no_exotic_physics_2026-09-18.md` 的 K1 约定 Λ^Δ=Λ^ph−Λ^sf，此处 Λ^ph=0，故 Λ^Δ=−V 对**全部 5 条边**成立。复算：空穴/电子 D 下非对角最小元 **−0.07**（仅 Γ1-Γ2、M1-M2 两条同族边违反，不是「全部 5 条」）；用 Part B 的二部划分 D={+,−,+,−,+} 则非对角最小元 **=0（Metzler 精确成立）**。λ_min=−λ_max 等价于二部，也就等价于 K1 下共轭后的 Metzler，所以 §0.2「二部机制不需要 Metzler 号约定」也不对。另：5 条边构成 4-圈加一个悬挂点，必然二部，Part B 几乎是平凡的；「二分不沿空穴/电子」只是因为数据里 Γ1-Γ2、M1-M2 两条同族边存在。（level：numerical-evidence，5×5 exact 算术。）
> - Part A/B/D 的数值可复现（agent）；标级仍是 numerical-evidence，没有越级。
> - 与 master `ceiling/zdom_iron_based_zero_dof_2026-09-21.md` 对 FeSe「2 带 vs 5 patch」的处理有张力：master 已订正为三通道共用同一个 2 带分区。

> 2026-09-19 · worktree `spinfluc-quant`（独立分支，UNMERGED，等团队统一归档）
> 脚本：`../_scratch_ir/_chk_fese_essenberger_sf_bipartite_2026-09-18.py`（已随本笔记入库）
> 前置笔记：[`../gamma/spin_fluctuation_structural_compatibility_2026-09-18.md`](../gamma/spin_fluctuation_structural_compatibility_2026-09-18.md)、
> [`spin_fluctuation_rigorous_tc_scoping_2026-09-18.md`](spin_fluctuation_rigorous_tc_scoping_2026-09-18.md)、
> [`../vertex/spm_gap_needs_no_exotic_physics_2026-09-18.md`](../vertex/spm_gap_needs_no_exotic_physics_2026-09-18.md)、
> [`../ceiling/lambda_ij_criterion_fails_on_FeSe_2026-09-18.md`](../ceiling/lambda_ij_criterion_fails_on_FeSe_2026-09-18.md)

---

## 0. 一句话

把 FeSe 的**真实、已发表**自旋涨落配对耦合数据（Essenberger et al., FeSe 五区 patch 分解）代入仓库现有的两套框架，得到两条具体、**numerical-evidence** 级的负向/定界结论：

1. **spm_gap 定理（PF/Metzler 机制）的前提 1 在这份真实数据上不成立**——FeSe 的自旋涨落耦合在同族（hole-hole、electron-electron）之间本身就是**排斥**的（Γ1-Γ2=0.01, M1-M2=0.07 > 0，且按 K1 号约定映射到 Λ^Δ 后仍为负），不存在任何"族内吸引"的声子项来补足该前提，因此该定理**不适用**于纯自旋涨落配对的 FeSe 模型。
2. 真正适用的是 `probes/spin_fluctuation_rigorous_tc_scoping` §3.3 的**另一套**机制（非负对称矩阵 V 的 λ_min，不需要 Metzler 号约定）：5 条已发表边构成的图恰好是**二部图**（`λ_min(V) = −λ_max(V)` 精确到浮点），但其二分方式是 `{Γ1,Γ3,M2} | {Γ2,M1}`——**不是**教科书期望的空穴/电子分裂 `{Γ1,Γ2,Γ3}|{M1,M2}`——且这个二分性质对被省略的另外 5 条边（未发表/APS SM 不可访问）的取值**敏感**，加入极小扰动即普遍破坏精确对称。

**没有产出**可与实验直接对照的"若 PF 成立则该材料 Tc/能隙对称性应为 X"陈述——数据不完整（只发表了 10 对 patch 耦合中的 5 对）是唯一障碍，已具体点名（§4）。

---

## 1. 现有框架对"媒介机制"的假设边界（任务步骤 1）

读 `spin_fluctuation_structural_compatibility`、`spin_fluctuation_rigorous_tc_scoping`、`spm_gap_needs_no_exotic_physics` 三篇后，可以把"机制无关的结构性质"与"需要具体机制才能算的量"分开：

**机制无关（对任何玻色媒介自动成立）：**
- 谱权重非负性（Källén–Lehmann：`χ''(q,ω)≥0`，任何玻色关联函数皆然，见 structural_compatibility §2.1）；
- `𝒪[f]` 算子在 f 非负、非增、ℓ¹ 时的 Perron/锥不变性机制（引理 M，不需要 f∈𝒫，只需要这三条解析性质，与"声子"这个标签无关）；
- spm_gap 定理本身（§3.2）：只要给定的多带耦合矩阵满足（1）族内 `Λ^ph_ii`（或任何扮演"吸引族内项"角色的量）≥0、（2）族间 `Λ^sf_ij` ≥0、（3）连通性，则不论这些系数来自声子还是自旋涨落还是任何其他玻色介质，PF 论证逐字成立——**这一步是纯代数，不依赖声子物理**。

**需要具体机制才能算（不能从结构论证白送）：**
- `f∈𝒫`（Bochner 正测度表示）本身，来自声子的谐振子谱表示，对一般自旋涨落核不自动成立（但也不是 spm_gap/引理 M 真正需要的——见上一条）；
- 单调性（关于 |n−m| 非增）——声子靠 `⟦k⟧` 的显式初等形式；自旋涨落配对核经过投影到符号变化通道（d 波/s± 形状因子）之后，**没有先验理由保持这条性质**（structural_compatibility §2.2(b)）；
- **最关键的一条，本文新确认的**：spm_gap 定理前提 1（"族内耦合非负/吸引"）本身是否成立，**取决于具体材料的真实耦合矩阵**，不是任何一般定理能保证的。声子物理里这一项通常由声子的族内吸引提供；纯自旋涨落配对（无声子项）不保证存在这样的族内吸引项——本文 §3 用 FeSe 的真实数据具体证实了这一点确实会失败。

---

## 2. 候选材料与文献源（任务步骤 2）

**材料：FeSe**（已由仓库 `ceiling/lambda_ij_criterion_fails_on_FeSe_2026-09-18.md` 独立确立为多带判据的材料侧反例来源，本文复用同一文献但换一套框架）。

- **一次文献**：Essenberger, Sanna, Buczek, Ernst, Sandratskii, Gross, *Ab-initio theory of iron based superconductors*, arXiv:1411.2121（= PRB 94, 014503 (2016)）——Fig.4(a) 插图 + Fig.5 图说给出 FeSe 自旋涨落配对耦合的**五区 patch 分解**（Fermi 面上 Γ1,Γ2,Γ3 空穴口袋、M1,M2 电子口袋）：

  | 边 | λ_ij | 图说标注 |
  |---|---|---|
  | Γ1-Γ2 | 0.01 | "同相位区域" → s± 通道内排斥 |
  | M1-M2 | 0.07 | 同上 |
  | Γ3-M1 | 0.19 | "连接 +/− 相位区域" → 有效吸引 |
  | Γ1-M1 | 0.18 | 同上 |
  | Γ2-M2 | 0.11 | 同上 |

  五个数**全为正**（幅值）；符号结构不在 `λ_ij` 矩阵本身，而由图说文字标注"哪两个 patch 的能隙相位相同/相反"给出——本笔记 `../literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md` §A.4/D.6 已 L1 转录。

- **微观定义与 DOS 平均约定的二次来源**：Essenberger, Sanna, Linscheid, Tandetzky, Profeta, Cudazzo, Gross, *Superconducting pairing mediated by spin-fluctuations from first principles*, arXiv:1409.7968（= PRB 90, 214504 (2014)）——Eq.(32) 给出 Λ^SF 的 TDDFT 定义，§VI.A 给出双指标平均公式（DOS 只除左指标，与声子侧同构），Eq.(37) 明确 `Im[Λ^SF_kk'(ω)]>0` 对所有 k,k'（自旋涨落通道**处处排斥**，符号变化必须由能隙函数本身携带，不由耦合矩阵携带——`iron_arsenide_lambda_ij_retrieval` 已记录为"符号不由 λ 矩阵携带"模式的第四次出现）。Eq.(38) 给出一个玩具模型，其中族内 SF 项**恰为零**（`I≠J` 条件）——一个已发表、物理动机明确、但仅限玩具模型（非材料计算）的 `λ11λ22=0<λ12λ21` 型反例。

- **次要/辅助数据点（未代入本文计算，仅记录）**：Lischner et al., PRB 91, 020502(R) (2015)——FeSe，给出对称投影后的**标量** s± 通道耦合 `λ[s±]=0.30`（裸能带）/`0.28`（含电荷涨落修正），对应 `T_c=5K`（下界，实验 8K），明确未考虑声子。这是一个**已经投影成标量**的自旋涨落耦合常数，可以直接喂给 γ-model/单通道 KAY 机制，但丢失了本文关心的 patch 分辨结构，故本文未使用。

---

## 3. 数值代入（任务步骤 3）

脚本 `../_scratch_ir/_chk_fese_essenberger_sf_bipartite_2026-09-18.py`（从零编写，仅用 `numpy`/`itertools`，无仓库依赖）。以 5 个 patch `{Γ1,Γ2,Γ3,M1,M2}` 为顶点、上表五条边为权重构造对称非负矩阵 `V`（未发表的另外 5 对——Γ1-Γ3、Γ2-Γ3、Γ1-M2、Γ2-M1、Γ3-M2——显式设为 0，这是本文唯一的近似，已在脚本文档字符串与下文 §4 明确标注）。

**Part A：直接谱分解。**
```
lambda_min(V) = -0.273570 ,  lambda_max(V) = +0.273570
lambda_min + lambda_max = 0  (精确，到浮点)
```
`lambda_min` 对应特征向量符号模式 `+-+-+`（按 `Γ1,Γ2,Γ3,M1,M2` 排序），**不匹配**教科书期望的空穴/电子分裂模式。

**Part B：二部图检验（穷举 `2^4=16` 种 2-染色，固定 Γ1=+1 去除全局符号冗余）。**
找到唯一二分：`{Γ1,Γ3,M2}` vs `{Γ2,M1}`——与 `probes/spin_fluctuation_rigorous_tc_scoping` §3.3 的"非负对称矩阵二部 ⟺ 谱关于 0 对称"机制精确对应，且与 Part A 的特征向量符号模式一致。**这不是**空穴 `{Γ1,Γ2,Γ3}` vs 电子 `{M1,M2}` 的朴素分裂——`Γ2` 和 `M1` 落在同一侧。判断：这几乎肯定是示意图只给出 10 对中 5 对的**产物**，不是物理预测（见 §4）。

**Part C：spm_gap.md §3.2 意义下的 Metzler/PF 检验。**
按仓库 K1 号约定构造 `Λ^Δ`：同族边（族内，Γ1-Γ2、M1-M2）取 `−V_ij`（排斥/破坏配对），跨族边（Γ-M）取 `+V_ij`（有效吸引）。取 `D=diag(+1 on Γ, −1 on M)`（朴素空穴/电子分裂，spm_gap.md 理论构造里预期的那个 D），计算 `D Λ^Δ D`：

**Metzler 条件在全部 5 条已知边上都失败**（10 个非对角项全为负，包括跨族边——因为 `D_iD_j=−1` 时符号被再翻转一次）。具体点名违反前提 1 的边：`Γ1-Γ2`（`D_ΓD_Γ·Λ^Δ = −0.01`）、`M1-M2`（`−0.07`）——这两条正是族内、且在这份数据里本身是排斥性质的，不存在声子项来补足 spm_gap 定理前提 1 要求的"族内 `Λ^ph_ii`≥0（吸引）"。

**结论（具体、可陈述）**：spm_gap 定理的三前提机制**不适用于** FeSe 的纯自旋涨落配对模型（在已发表的这 5 条边范围内）——原因不是定理错，是这份材料数据根本不满足定理要求的输入结构（缺一个族内吸引的"声子项"）。真正适用的是 Part A/B 里的 λ_min(V) 二部图机制，它不需要 Metzler 号约定，只需要 V 逐元非负——这与 `probes/spin_fluctuation_rigorous_tc_scoping` §3.3 判定 (b)（结构部分可搬，但需要额外前提）里"若不做单热点简化、保留真实多带结构，问题变成矩阵值算子谱问题、需要仓库 ceiling 方向的矩阵推广"这句预判，在本文里**得到一次具体材料实现的确认**。

**Part D：敏感性。** 对 5 条未发表边逐一开启幅值 `eps∈{1e-4,...,0.05}`（保持逐元非负）：`lambda_min+lambda_max` 从精确 0 单调偏离（`eps=1e-2` 时已到 `7e-3`，约 2.6% 相对偏移），且在 `eps=0.05` 时特征向量符号模式已从 `+-+-+` 变为 `+++-+`——即 Part B 找到的二分方式**对缺失数据敏感**，不是这 5 个数字之外还有鲁棒性的结构预测。

---

## 4. 缺什么（任务步骤 4，具体点名）

**唯一、具体的缺口**：Essenberger et al. 两篇论文均未在可访问的正文/图中给出完整的 patch 分辨 `Λ^SF_ij(ω)` 矩阵（10 对中只给了 5 对示意性数值，且是积分后的示意值不是逐频率函数）。完整数据的存放位置极可能是 PRB 94,014503 的 Supplemental Material——`../literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md` 已记录该 SM 链接在本环境下返回 **403**（APS harvest 路线此前已验证对全文 PDF 可用，但 SM 是另一条独立链接，未验证是否同样可绕过）。

若拿到完整 10×10（或更细的 patch 划分）矩阵，可以：
1. 直接重跑本文 Part B/C，检验二部图结构、Metzler 前提失败是否是"数据不全的假象"还是"FeSe 这份耦合矩阵本身结构性的"；
2. 若矩阵确认仍是二部但非空穴/电子分裂，这将是一个**可与角分辨光电子能谱（ARPES）能隙相位测量直接对照**的具体预测（"能隙相位分界不沿空穴/电子边界，而沿 `{Γ1,Γ3,M2}|{Γ2,M1}` 这条线"）——但**当前数据不足以支撑这个陈述**，只能停在"5 条边给出这个二分，对缺失数据敏感"这一步。

因此，任务要求的"若 PF 结构假设成立，该材料 Tc/能隙对称性应为 X"这类完全可与实验对照的陈述**未能产出**——障碍是具体的（缺 5/10 条边的发表数据，来源大概率是不可访问的 APS SM），不是原则性的。

---

## 5. 严格性标级（逐条，不越级）

| 陈述 | 级别 | 依据 |
|---|---|---|
| spm_gap 定理前提 1/2/3 的抽象结构（不含具体材料代入） | `analytic`（转引 `vertex/spm_gap_needs_no_exotic_physics_2026-09-18.md`，未重新证明） | 该文件 §3.2 |
| FeSe 五区数据 Γ1-Γ2=0.01 等五个数值 | `numerical-evidence`，输入为 **L1-读图**（转录自 `literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md` §A.4，该文件本身标注为读原文图说） | 见 §2 引用链 |
| `V` 的 5 边子图是二部图，二分为 `{Γ1,Γ3,M2}\|{Γ2,M1}` | `numerical-evidence`（本文脚本从零计算；穷举 16 种染色，非启发式） | 本文 §3 Part A/B，脚本已运行确认 |
| spm_gap 定理前提 1 在此数据上不满足（Γ1-Γ2、M1-M2 边违反 Metzler） | `numerical-evidence`（同一脚本 Part C，代数检验，非近似） | 本文 §3 Part C |
| 二分性质对缺失 5 条边敏感（Part D） | `numerical-evidence`（脚本 Part D，5 个 eps 取值下的定性观察，非穷举扫描） | 本文 §3 Part D |
| 完整 patch 矩阵的缺失定位到 APS SM / 403 | `numerical-evidence`（转引 `literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md`，未重新尝试该 SM 链接） | 见 §4 |

**未越级声明**：本文不含 certificate 或 rational-exact 级陈述（所有输入是浮点文献转录值，未做区间/有理算术）；没有产出任何可直接与实验对照的 Tc/能隙对称性预测（§4 已明确该缺口）；Part B/C 的"二分/Metzler 失败"结论建立在只有 5/10 条边数据的基础上，Part D 已明确标注该结论对缺失数据的敏感性，不宣称是 FeSe 真实耦合矩阵的鲁棒性质。

---

## 6. 与已有笔记的关系

- 不修改、不取代 `vertex/spm_gap_needs_no_exotic_physics_2026-09-18.md` 的任何定理陈述——该定理本身仍然是 `analytic` 且成立；本文只是找到一个**真实材料**使其**前提不满足**，这是对定理适用范围的具体材料级确认，不是对定理的挑战。
- 呼应 `gamma/spin_fluctuation_structural_compatibility_2026-09-18.md` §2.2(b) 的预判（"若保留真实多带结构，问题升级为矩阵值算子谱，需要仓库 ceiling 方向的矩阵推广，且那条线已知有非平凡反例"）——本文用 FeSe 五区数据给出这句预判的第一个具体材料实现。
- 与 `ceiling/lambda_ij_criterion_fails_on_FeSe_2026-09-18.md` 互补而非重复：后者检验的是**声子** `λ_ij` 矩阵在多带 ceiling 判据 `λ11λ22≥λ12λ21` 下的表现（FAIL，2.48×）；本文检验的是**自旋涨落** `Λ^SF_ij` 矩阵在 PF/Δ通道号结构判据下的表现（前提不满足）——两者用同一材料、同一篇文献的两个不同数据集（声子 Fig.3 vs 自旋涨落 Fig.4/5），回答不同问题，互不覆盖。

---

## 7. 待办 / 下一步（若决定继续投入）

- [ ] 尝试 APS SM 链接的 403 是否可用其他路线绕过（仓库已知 harvest.aps.org 全文 PDF 路线对正文有效，SM 未测试）。
- [ ] 若拿不到完整矩阵，可考虑直接联系作者索取数据（`literature/` 目录已有类似邮件草稿的先例，未检查是否已含此文）。
- [ ] 若要推进"可检验预测"，下一个最小步骤是把 Lischner et al. 的标量 `λ[s±]=0.30/0.28` 代入 γ-model/单通道 KAY 机制做一次独立交叉检验（不同文献源，不同框架），可能比等待 patch 分辨矩阵更快产出一个真正可对照实验的数字。
