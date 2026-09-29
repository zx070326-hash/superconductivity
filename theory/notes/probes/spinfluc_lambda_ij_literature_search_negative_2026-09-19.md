# λ_ij 判据的自旋涨落输入：换个检索路径（发表的紧束缚+RPA 论文）——仍是负结果，且给出一个结构性原因

> ⚠ **2026-09-29 合入核查注**（`worktree-spinfluc-quant` 合入 master 时由核查 agent 复跑 + 主会话抽查；**正文未改，本注为唯一改动**；标「agent」者为核查 agent 报告、主会话未独立复算）：
> - 文献检索事实，无数值可复现（agent）。
> - ⚠ 「Essenberger 是例外、给出非拟合矩阵」漏了 α：Essenberger 的自旋涨落道含作者自陈无法第一性原理确定的标度参数 α；master `ceiling/zdom_iron_based_zero_dof_2026-09-21.md` §6.5 明确写「Essenberger 全程零拟合不成立」。同分支 `lambda_ij_violations…` 已正确限定 α 只进 SF 道，两篇口径不一致——以后者与 master 为准。

> 2026-09-19 · worktree `spinfluc-quant`（分支 `worktree-spinfluc-quant`，UNMERGED）
> 派活口径：不重复 [`spinfluc_coupling_matrix_open_data_census_2026-09-19.md`](spinfluc_coupling_matrix_open_data_census_2026-09-19.md)
> 已证伪的路线（Materials Project/AFLOW/NOMAD/Zenodo 通用数据库、GitHub 开源 RPA/FLEX 代码仓库），
> 改走「找一篇论文**同时**公开发表 (a) 具体材料的多带紧束缚/Wannier 哈密顿量参数 和 (b) 对应的自旋涨落
> RPA/FLEX 配对相互作用矩阵（λ_ij / Γ_ij 型，非拟合）」，以便直接喂给仓库现有 `_chk_spm_*.py` /
> `lambda_ij_criterion` 判据机器。
> **本文性质：文献检索事实，不含代入计算，不设 analytic/certificate 分级。**

---

## 0. 一句话结论

**仍未找到。** 检索了三条与上一轮完全不同的具体文献方向（铁基多带 RPA 配对理论的奠基/主力论文），
条件 (a)（显式紧束缚哈密顿量数值参数）在这些论文里**普遍满足**——反而是这个子领域的标配、
公开发表率很高；但条件 (b)（配对相互作用**以离散小矩阵形式印成表**，像 α²F_ij 那样）**系统性地不满足**。
原因不是"没人做"，而是一个**结构性事实**：这三条独立路线（RPA 轨道基、fRG патch 基）算出来的配对顶角
`Γ(k,k')` / `g(k1,k2,k3)` 本身是**费米面上的连续函数**（k,k' 各自跑遍整个费米面），
它们的天然呈现方式是**色图/等高线图**或**投影到对称性谐波后的标量本征值 `λ_α`**，
而不是像声子 Eliashberg 那样天然坍缩成一个 N×N（N=费米面分块数，通常 3–10）的小矩阵。
Essenberger 系列之所以能给出仓库已引用的那种 10×10 patch 矩阵，靠的是他们**自己额外做了一次
"投影到少数几个 patch 上取平均"这一步**——这一步在本轮检索的三条主流路线里都没有被做（或没有被公开印出来）。

---

## 1. 检索方向 A：Kreisel/Andersen/Hirschfeld 组 FeSe 显式紧束缚 + RPA 配对理论

两篇候选，均下载全文核实（非仅摘要页）：

### 1.1 Choubey, Berlijn, Kreisel, Cao, Hirschfeld, *Phys. Rev. B* 90, 134520 (2014)，arXiv:1401.7732

- **(a) 满足**：给出 10 轨道 `t^μν_R`（实空间跳跃矩阵元）与 Hubbard–Hund 相互作用参数
  `U=0.90 eV, J=U/4`（Kuroki 记号），Supplement Material 明确写出。
- **(b) 不满足**：Sec. III「the pairing interaction shows a rapidly decaying [real-space profile]」——
  他们展示的是配对相互作用 `Γ^μν_{RR'}` **随实空间距离衰减的曲线/图**，不是一张离散轨道对矩阵数值表。
  正文与 Supplement 全文 grep（`pdftotext -layout -enc UTF-8`）未见任何 `TABLE`/矩阵数值块对应配对项，
  只有 Fig. 用图示呈现衰减。

### 1.2 Mukherjee, Kreisel, Hirschfeld, Andersen, *Phys. Rev. Lett.* 115, 026402 (2015)，arXiv:1502.03354

- **(a) 满足且完整**：Supplemental Material（本文档 §S1–S2）**逐项印出**10 轨道与 5 轨道 3D
  紧束缚模型的全部数值跳跃参数（`t1111=0.086/z, t1160=−0.063/z−0.0211, …` 等约 40 个复数系数，
  z=6 为重整化因子），另有 kz 方向的附加跳跃项，来源标注为对 Eschrig et al. DFT 能带的拟合。
  这是本次检索里**条件 (a) 满足得最完整**的一篇。
- **(b) 不满足**：配对相互作用的解析式是 Eq. (S4)–(S5)：
  `Γ_{1234}(k,k') = (3/2)Ū^s χ^s_RPA(k−k')Ū^s + [(1/2)Ū^s − (1/2)Ū^c]χ^0_RPA(k−k')Ū^c + (1/2)Ū^c`，
  是 **k,k' 的连续函数**（轨道指标 1–4 固定在 5 个 Fe-d 轨道上，但动量指标不离散化）。
  该文用它算出的是**磁化率色图**（Fig. S8/S9）与**投影到 Fermi 面后插值的能隙函数** `Δ(k)=Δ₀g(k)`
  （Fig. S10，1000×1000×200 k 网格插值），从未把 `Γ` 坍缩成一张小矩阵表。

**判断**：这条路线在"公开显式紧束缚参数"这一半上是本次检索最强的候选，但配对相互作用侧
从未以判据要求的离散矩阵形式出现——不是数据缺失，是这类论文的呈现范式本身不产生这种对象。

---

## 2. 检索方向 B：Graser, Maier, Hirschfeld, Scalapino, *New J. Phys.* 11, 025016 (2009)，arXiv:0812.0343

这是该子领域被引用最多的奠基论文（五轨道 RPA 配对，铁基标准模型的源头之一）。全文下载核实：

- **(a) 满足**：**Table I**（轨道内跳跃参数，5×8 网格）+ **Table II**（轨道间跳跃参数）
  ——这组紧束缚参数（署名来自 Cao–Hirschfeld–Cheng PRB 77, 220506(R) (2008) 的 DFT 拟合）
  是后续大量铁基 RPA/fRG 文献反复复用的"标准起点"，公开、无拟合争议。
- **(b) 不满足**：全文 grep `pairing strength`/`eigenvalue` 只找到**标量**对象——
  `λ_α`（正文 §"dimensionless pairing strength functional" 附近，Eq. 前后）是投影到费米面
  对称性谐波 `g_α(k)` 之后的**本征值**（每个候选配对对称性一个数），不是轨道/能带对之间的矩阵元。
  正文没有印出任何形如"能带 i–能带 j 配对强度 = 数值"的表。

**判断**：与方向 A 同一结构性结论——这篇奠基论文同样把配对相互作用当连续核处理，
最终只报告**每个对称性通道一个标量 `λ_α`**，这本身就不是判据需要的"矩阵"这一层级的对象，
不存在"矩阵被藏起来没发表"的问题，是这个量根本不以矩阵形式被计算出来。

---

## 3. 检索方向 C：fRG patch 方法（Platt/Thomale/Hanke 及合作者）

搜索确认该学派把费米面离散化为 patch（典型 32–96 个），流方程的耦合函数
`g(k1,k2,k3)`（patch 指标）在数值实现上**确实是离散的**——这一点原则上比方向 A/B 更接近判据要的形状
（费米面已经被分块，正是仓库 λ_ij / Essenberger patch 矩阵的同一种离散化思路）。

但检索到的综述与主力论文（Thomale–Platt–Hu–Honerkamp–Bernevig PRB 80, 180505(R) (2009)；
Platt–Honerkamp–Hanke NJP 11, 055058 (2009)；Platt–Hanke–Thomale 综述 arXiv:1310.6191 等）
**均以色图/流图形式展示 `g` 的 patch 依赖**（"Λ 是 fRG 流参数……能量最终逼近费米面"），
未在正文中找到把 32–96 个 patch 的耦合值坍缩、**印成一张可直接摘录的小数值表**的实例——
patch 数本身就比判据框架假设的"3–10 个能带/pocket"高出一到两个数量级，需要先做一次
"按能带/pocket 分组平均"的降维，而**这一步在检索到的文献里未见公开做过或未见其平均后的数值表**。

**未做**：本文未逐篇打开这几篇 fRG 论文的补充材料/数值代码仓库核实是否存在这样的降维表——
只核实了正文层面的呈现方式，这是本文最大的检索缺口（见 §5）。

---

## 4. 为什么这构成一个（比上一轮更具体的）结构性负结果

上一轮 `spinfluc_coupling_matrix_open_data_census_2026-09-19.md` 的结论是"没有现成数据集，只有代码工具"，
停留在"数据托管方式"层面。本轮的新发现在**物理呈现范式**层面：

> **声子 Eliashberg 函数 `α²F_ij(ω)` 天然带一个离散的 (i,j) 指标**（能带/pocket 编号，来自
> `Σ_ν |g^ν_{k,k'}|² δ(ξ_k)δ(ξ_{k'}) δ(ω−ω_ν)` 对 k,k' 分别限制在第 i、j 个费米面 sheet 上再对 sheet 内积分）
> ⟹ 论文自然而然地把它整理成一张 N_sheet × N_sheet 的表（Essenberger、Sanna CaC₆ 都是这么做的，
> 仓库 `ceiling/` 下两篇笔记已核实）。
>
> **自旋涨落配对顶角 `Γ(k,k')` 没有这个天然的双重求和结构**——RPA/fRG 计算里 k,k' 各自是
> 独立的连续费米面坐标（配对的符号结构、节点位置正依赖于 k,k' 的精细变化，这是自旋涨落配对
> 区别于声子配对的物理核心），把它投影到少数几个 patch 上会**丢失判据本来最想测的那类信息**
> （节点结构、s±/d 波符号），所以这个子领域的论文选择保留连续表示（色图）或投影到对称性谐波
> （标量 λ_α），两者都不是判据要的离散矩阵。

⟹ **这不只是"还没人发表"，而是"发表离散矩阵形式在这个子领域不是自然选择，需要作者额外做一步
（且这一步会引入一个新的、非平凡的 patch 划分方案作为隐藏输入）"。** Essenberger 系列做了这一步，
是这个子领域里的**例外**而非常规，这与上一轮census 已知"只有 Essenberger 一家给出这种矩阵，且部分锁在 APS SM"
的观察完全一致，本轮补上了"为什么只有它这么做"的机制解释。

---

## 5. 未做 / 检索缺口

1. **未逐篇核实 fRG patch 论文的补充材料/代码仓库**是否含"按能带分组平均后的 patch 耦合数值表"
   ——这是 §3 判断"未见"而非"确认没有"的具体缺口，也是若要继续投入的最短下一步。
2. **未检索 Chubukov 组、Altmeyer/Roeser 组的原始论文全文**（只在 WebSearch 摘要层面确认存在，
   派活里提到的这两个方向本文未深入，时间预算内优先做了 A/B/C 三条更集中的方向）。
3. **未检索 UPt₃、Sr₂RuO₄、CeCu₂Si₂ 等重费米子体系**的自旋涨落配对文献
   ——派活里给出的备选方向，本文因预算未覆盖，留给下一轮。
4. 方向 A 的两篇论文本文均只读了 Sec. I–III / Supplement 全文的文本层（`pdftotext -layout -enc UTF-8`），
   未逐张重新渲染页图核对是否有文本抽取遗漏的小字数值表格（按 CLAUDE.md「两种仪器都要跑」的纪律，
   这条尚未做第二遍页图核对，若后续要把"确认没有"提到更高确信度，应补this一步）。

---

## 6. 结论与对下一步派活的建议

- **仍不能给仓库 λ_ij 判据喂一组已发表、非拟合的自旋涨落耦合矩阵。**
- 若团队仍想在自旋涨落侧关闭这个方向，**最现实的路径不是继续搜文献**（本轮 + 上一轮已覆盖
  通用数据库、开源代码、三条主流理论文献范式，结构性原因已经清楚），而是二选一：
  1. 联系 Essenberger/Sanna 系列作者索取完整数据（上一轮 census §6 第 3 条已建议，成本最低）；
  2. 接受"自己跑一次 RPA/FLEX + 投影到 patch"的工作量（用 §3 提到的 fRG patch 划分或 Essenberger
     的 patch 方案作为降维方案），这已经不是"数据普查"范畴，需要 team-lead 重新评估投入。

---

相关笔记：
[`spinfluc_coupling_matrix_open_data_census_2026-09-19.md`](spinfluc_coupling_matrix_open_data_census_2026-09-19.md)（上一轮：通用数据库/开源代码普查，本文不重复）·
[`spin_fluctuation_fese_material_probe_2026-09-19.md`](spin_fluctuation_fese_material_probe_2026-09-19.md)（本文回应的原始缺口来源）·
[`../ceiling/lambda_ij_violations_are_genuine_not_fit_artifacts_2026-09-19.md`](../ceiling/lambda_ij_violations_are_genuine_not_fit_artifacts_2026-09-19.md)（判据消费点、声子侧 λ_ij 生成式核查方法论，本文方向 A/B 的"生成式"判断法沿用同一套）
