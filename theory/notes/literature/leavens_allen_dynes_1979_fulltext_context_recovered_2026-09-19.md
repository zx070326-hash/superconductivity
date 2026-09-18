# LAD79 优先权风险复查（第二轮）——通过现代引用上下文间接重建证明内容

**任务**：2026-09-18 笔记 `leavens_allen_dynes_1979_priority_check_2026-09-18.md` 判定"低-中风险偏低，
但未完全排除（正文始终未读到）"。本轮受 team-lead 委托，穷尽新渠道尝试获取
Leavens, Allen & Dynes, *Solid State Commun.* **30**(10), 595–596 (1979)（下称 LAD79）全文，
并重新评估该风险判定。

**结论先行：全文仍未获取（渠道已穷尽），但本轮找到一篇 2026 年开放获取论文，
在正文（非摘要）里逐字给出 LAD79 实际证明内容的现代转述，是本笔记链条里目前最强的一手证据。
基于这条新证据，风险判定从"低-中偏低、不能完全排除"下调为"低，且可以较高置信度排除
LAD79 与仓库定理 2 / Einstein 极值性构成实质重复"。**

---

## 1. 本轮新增检索渠道（均未获取原文全文）

| 渠道 | 结果 |
|---|---|
| NRC Publications Archive (`nrc-publications.canada.ca`) 站内搜索 "Leavens superconducting" | **零命中**——该仓库收录的是 NRC 技术报告，1970s 期刊论文（版权归 Elsevier）从未存档于此 |
| archive.org `advancedsearch.php` 查 `sim_solid-state-communications_1979*` | **`numFound: 0`**——该刊 1979 卷未被 archive.org 数字化（对比：1992/1997 卷有扫描件） |
| Google Books 检索该文标题 | 无命中（期刊短讯，非书籍收录） |
| Semantic Scholar `/citations` API（含 contexts 字段） | 拿到 4 篇引用该文的论文列表，其中两篇是 2025/2026 年 Thai 团队论文——**本轮的突破口** |

## 2. 突破：Tasee, Natkunlaphat & Pinsook (2026)，*J. Met. Mater. Miner.* **36**(2), e2393 —— 一手全文

**开放获取**（Chulalongkorn 大学期刊，`jmmm.material.chula.ac.th`，本轮亲自下载 PDF 并逐页 Read，
非 AI 摘要转述）。这是现代文献中**唯一**一篇在正文（不是摘要/参考文献列表）里明确说出
LAD79 具体证明了什么的来源。原文引言部分逐字（p.2，编号为文中 Eq.(11) 之后）：

> "In the limit of λ → ∞ (extreme case), or approximately λ ≥ 10,
> `T_c^AD → F(μ*)√(λω̄₂)` (11)
> where F(μ*) is a function of the Coulomb pseudopotential, μ*, due to the screened Coulomb
> interactions between electrons. … The expression in Equation (11) is regarded as the asymptotic
> behavior of Tc. **Nevertheless, it is more rigorous to derive Equation (11) under the condition
> ω_max/2πT_c >> 1, rather than relying solely on the parameter λ [7].**"

其中 `[7]` 就是 LAD79（参考文献列表逐字核对：`[7] C. R. Leavens, P. B. Allen, and R. C. Dynes,
"The superconducting transition temperature for very large electron-phonon coupling," Solid State
Communications, vol. 30, no. 10, pp. 595–596, 1979`）。

**这句话是本笔记链条里第一次读到、由现代作者在正文中明确归属给 LAD79 的具体技术内容**：
LAD79 的贡献被理解为——把"渐近公式 (11) 何时成立"这一判据，从"λ 需要多大"（一个不精确的、
依赖谱形状的说法）**改写为一个更严格的无量纲条件 `ρ_max ≡ ω_max/(2πT_c) ≫ 1`**（与谱截止频率
和 Tc 本身的比值相关，而非仅仅是耦合常数 λ）。

## 3. 这与仓库定理的关系：技术命题不同类

这条被现代作者复述的内容，与本仓库 `P4_g2_analytic.md` / `P4_einstein_extremal.md` 证的东西
是**两个不同结构的命题**：

- **LAD79（据此转述）**：对**固定形式**的渐近公式 `T_c → F(μ*)√(λ⟨ω²⟩)`（系数 F(μ*) 已经写定），
  刻画这个近似式**何时数值上准确**——用 `ρ_max = ω_max/2πT_c ≫ 1` 这个条件代替粗糙的
  "λ 很大"。这是一个**收敛域/近似有效性**命题，本质上关于渐近展开的下一阶修正何时可忽略。
- **仓库定理 2 / 定理 E（Einstein 极值性 + C_∞ 上界）**：对**任意**谱测度 P（不限定渐近区、
  不假设已经知道系数 F(μ*)），证明 `T_c ≤ C_∞√(λ⟨ω²⟩)` 对**所有** λ、**所有**谱形状成立，
  且刻画等号何时取到（Einstein 谱、μ*=0）。这是一个**变分/极值性**命题。

**这两类命题在逻辑上互不蕴含**：知道"公式 (11) 在 ρ_max≫1 时准确"不能推出"该公式的系数
是所有谱形状里的最大值"（后者需要证明其余谱形状不会超过这个系数，这是完全不同的论证）；
反过来，证明了 Einstein 谱是极大值也不会自动给出"这个渐近逼近何时数值准确"这个收敛速度问题
的答案。

## 4. 与第一轮判定的交叉验证

这条新证据**独立支持**了第一轮笔记 §2.1 纯靠摘要字面义做出的判读（"range of validity of the
asymptotic relation" 更像"数值收敛区间"而非"对任意谱测度的严格不等式"）——现在有了
**正文级**的现代转述，不再只是对摘要句子的语言学猜测。同时也与第一轮 §2.3 的负面证据一致：
Carbotte RMP62(1990) 与 arXiv:2504.17435（同一 Thai 团队更早一篇）都没有把 LAD79 当作
"极值性/上界"的来源来引用，而是用于别的场合。

本轮另找到一篇独立文献 **arXiv:2407.12922**（"Fundamental limits on the electron-phonon coupling
and superconducting Tc"，2024，全文 html 已读）作为第三个交叉验证点：该文在讨论渐近公式
`T_c → 0.18√(λ⟨ω²⟩)`（其 Eq.7）时，引用的是 **Leavens 1975（SSC 17,1499，单独作者）+
Allen-Dynes 1975 + KAY**，参考文献表**完全没有 LAD79 词条**。这是继 Carbotte(1990)、
KAY I/II/III 之后**第三个**独立的、专门处理这个常数历史渊源的现代文献都不引用 LAD79
作为该不等式来源的例子。

## 5. 更新后的风险判定

**判定：低风险（较高置信度），不构成对仓库"首次严格证明"措辞的实质威胁。**

支持这一判定的证据链（按强度排序）：
1. **（本轮新增，最强）** 现代同题材开放获取论文（Tasee et al. 2026）在正文中明确将 LAD79
   的贡献归属为"渐近公式的有效性条件（`ρ_max≫1`）"，而非"对任意谱形状的严格不等式"。
2. 三篇独立现代文献（Carbotte 1990、KAY I/II/III、Trachenko-Setty-Zaccone 组 2024）在专门
   讨论这同一个 0.18/0.1827 常数的历史来源时，均未引用 LAD79。
3. 标题聚焦"very large electron-phonon coupling"极限区间，该区间内 AD75 自己已指出谱形状
   效应趋于消失，极值性问题在此区间趋于平凡。
4. 参考文献列表（Crossref 一手-元数据）不含任何形状-变分证明工具类文献。

**仍未闭合的缺口（不影响上述判定方向，但严格意义上原文仍未读到）**：
- LAD79 原文（2 页）本身仍未获取——本笔记的判定基础是**第三方现代转述**，属于二手，
  按仓库纪律不能升级为 analytic/certificate 级别的"已排除"。
- 现代转述本身也可能简化或误读了 1979 年原文的确切数学内容——不能完全排除 LAD79 原文
  在这句"ρ_max≫1 条件"之外，还包含某种更接近极值性的旁支论证（可能性低，无证据支持，
  但不能证伪，因为没读到原文）。

## 6. 对仓库对外措辞的建议（更新）

维持第一轮笔记 §4 的建议方向，但可以把语气从"不确定、待核实"调整为"倾向性明确、但引用需注明间接来源"：

> "Leavens, Allen & Dynes (Solid State Commun. 30, 595, 1979) present 'a rigorous derivation of
> the range of validity' of the large-λ asymptotic relation `k_BT_c ≈ 0.1827[λ⟨ω²⟩]^{1/2}`. Based
> on a modern citing context (Tasee, Natkunlaphat & Pinsook, J. Met. Mater. Miner. 36, e2393, 2026,
> open access), this result concerns the convergence condition `ω_max/(2πT_c) ≫ 1` under which the
> fixed-coefficient asymptotic formula becomes accurate — a different technical statement from the
> shape-independent extremal inequality `T_c ≤ C_∞√(λ⟨ω²⟩)` for general spectral measures P proved
> here. We were unable to obtain the primary source directly (paywalled, no free full text located
> despite exhausting NRC Canada's institutional repository, archive.org, Google Books, and author
> self-archive checks); this characterization rests on third-party citation context, not a direct
> reading of the 1979 text."

**不建议**声称"已排除"或"已证实无重叠"——仍应标注"基于二手转述，非一手核实"。

## 7. 「读到的」vs「推断的」清单

**一手（本次亲自下载并逐页 Read 的 PDF 正文）**：Tasee, Natkunlaphat & Pinsook, *J. Met. Mater.
Miner.* 36(2), e2393 (2026) 全文 5 页，含 Eq.(11) 前后段落及参考文献列表 [7]。

**一手（本次 WebFetch 全文 html）**：arXiv:2407.12922v2 引言段落与完整参考文献列表 1–61
（未见 LAD79 词条）。

**一手-元数据（本次查询）**：NRC Publications Archive 站内搜索零命中；archive.org
advancedsearch API `numFound=0`；Semantic Scholar citations API 的 4 篇引用列表。

**二手（延续自第一轮笔记，本次未重新取件）**：Carbotte RMP62(1990) p.1042 相关段落；
Crossref 12 条参考文献列表。

**推断（本笔记新增的判断，非任何来源直接陈述）**：§3 的"两类命题互不蕴含"论证；
§5 的整体风险等级判定。

---

关联：[[leavens_allen_dynes_1979_priority_check_2026-09-18]]、
[[lit_coulomb-optimal-spectrum-variational-2026-09-18]]（researcher 跨会话 memory，非本仓库文件）、
`theory/derivations/P4_einstein_extremal.md`、`theory/derivations/P4_g2_analytic.md`
