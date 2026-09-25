# η 数值表第二轮语料搜索——0 个可验证命中，会话级 WebFetch 完全失能

**日期**：2026-09-25 · **等级**：见下 · **状态**：live（诚实的空结果 + 工具故障报告）
**派活**：team-lead 要求在**新语料**里找已发表的 McMillan–Hopfield 参数 η（= N(E_F)⟨I²⟩）数值表配原子质量
M，喂给 `kTc ≤ C_∞√(η/M)`（`C_∞ = 0.1827262477`，analytic，λ 无关）。先读了
[`eta_numerical_tables_not_published_2026-09-25.md`](eta_numerical_tables_not_published_2026-09-25.md)
（Shipley 2021 / Moussa–Cohen 2006 两篇 0 命中），本篇找的是**另外的**语料。

## 0. 本轮的硬约束：WebFetch 本会话完全失能

这**不是**已知的 `arxiv.org`/`www.nature.com` TLS 劫持（memory
`arxiv-and-nature-are-tls-intercepted-on-this-machine`——那个故障给的是证书错误，只挡两个域名）。
本轮遇到的是另一种错误，**任意域名**都给同一句：`Unable to verify if domain X is safe to fetch`。
测试过 5 个域名，逐一失败：`harvest.aps.org`、`en.wikipedia.org`、`github.com`、
`arxiv.org/abs/2508.21281`、以至 **`https://example.com`**（无害对照，同样失败）。
⟹ **工具级失能，不是网络/域名事实**，本仓 `library/pdf/` 111 个已有文件里也没有任何一篇匹配
Gaspari–Gyorffy / McMillan 1968 / Papaconstantopoulos / Klein–Boyer–Papaconstantopoulos / Savrasov 系列
（`Glob **/*aspari*` `**/*cMillan*` `**/*apaconstant*` `**/*Klein*` 全部 0 命中）。
**⟹ 本轮我读不到任何一篇候选论文的正文，一个仪器都没有，谈不上两台仪器交叉。**
已向主会话报告此故障（`SendMessage` to main，2026-09-25），未收到回复前不再重试 WebFetch。

按 CLAUDE.md「记录丢失时不得反向重建」的同一精神：**没有正文时，只能报告「查到了什么语料线索」，
不能把 WebSearch 摘要当作已验证的数值。** 下表全部标记 `未取得、未定`。

## 1. 候选语料清单（全部只有 WebSearch 摘要，未读正文，0 个可计入「命中」）

| 语料 | 年份 | 为什么可能有 η 表 | WebSearch 摘出的疑似数字（**未核实，禁止消费**） |
|---|---|---|---|
| Gaspari & Gyorffy, PRL **28**, 801 (1972) | 1972 | η 理论本身的原始论文，RMTA 公式来源 | 无 |
| McMillan, Phys. Rev. **167**, 331 (1968) | 1968 | η=N(0)⟨I²⟩ 概念提出处，含经验值表 | 摘要转述「N(0)⟨I²⟩ 只变化约 2 倍，N(0) 与 ⟨I²⟩ 各变化约 10 倍」，无表格数字 |
| Papaconstantopoulos, Klein 等, Phys. Rev. B **15**, 4221 (1977) | 1977 | WebSearch 自称"计算了 32 个 Z≤49 金属的超导性质"——**若属实，是迄今找到的最大候选表** | 无（未验证该篇是否真存在/真有此内容，仅 WebSearch 转述） |
| Klein, Boyer, Papaconstantopoulos, Mattheiss, Phys. Rev. B **18**, 6411 (1978) | 1978 | Nb₃Sn/V₃Si 的 APW 计算，据称是 A15 化合物 η 值的常被引用源头 | 无 |
| Klein, Boyer, Papaconstantopoulos, PRL **42**, 530 (1979) | 1979 | 10 个 A15 化合物（V₃X, Nb₃X）的 λ 与 Tc，理应基于 η 表 | 无 |
| Radevych, Shishidou, Weinert, Margine, Kolmogorov, Mazin, npj Comput. Mater. (**2026**)，DOI `10.1038/s41524-026-02141-7`（= arXiv:2508.21281） | **2026，本年内** | 明确以 RMTA 重算 η，逐角动量分解，**对 V/Mo/Nb 及其 A15/L1₂/B1 化合物与 FLAPW 基准比对**，5% 一致 | 摘要提及"5% 一致"字样，无具体数字；配套 GitHub `danylo-radevych/pp-rmta` 有 bcc Nb、A15 V₃Si 算例 |
| Savrasov & Savrasov, Phys. Rev. B **54**, 16487 (1996)（= cond-mat/9604124） | 1996 | Al/Cu/Mo/Nb/Pb/Pd/Ta/V 线性响应电声耦合；有非 arXiv/APS 直挂 PDF（`savrasov.physics.ucdavis.edu/Works/Publications/epiprb.pdf`）——若 WebFetch 恢复，这是**唯一不落在已知封锁域名上**的候选 | 摘要只给出 Pd 的 λ（0.35/0.37/0.69 三个不同口径），无 η 表数字 |
| arXiv:1910.08312（(TaNb)₀.₆₇(HfZrTi)₀.₃₃ 高熵合金） | ~2019 | WebSearch 摘要**声称**含 Ta/Nb/Hf/Zr/Ti 的 η 值（mRy/a_B²，"每原胞"） | **摘要给的数字：Ta 155, Nb 152, Hf 54, Zr 62, Ti 48**（单位 mRy/a_B²/primitive cell）——⚠ 这是 LLM 摘要转述，我一个字都没读到原文，**不得当命中使用** |
| arXiv:2308.03153（(ScZrNb)(RhPd) 高熵合金超导体） | 2023 | WebSearch 摘要提及 bcc Nb 的 η≈76 mRy/a_B²（"等效原胞"） | 同上，**未核实** |
| arXiv:cond-mat/0111385（高压硼的第一性原理超导性研究） | 2001 | WebSearch 摘要给出「η=7.627 eV/Å² 对应 Nb」「η=14.13 eV/Å² 对应 467 GPa 金属氢」 | 同上，**未核实** |

## 2. 判决

**本轮 0 个可用 (η, M) 数对**——但理由与
[`eta_numerical_tables_not_published_2026-09-25.md`](eta_numerical_tables_not_published_2026-09-25.md)
**完全不同、不可合并**：
那篇是**读了全文**之后确认两篇论文真的没有 η 表（干净否定，`numerical-evidence` 级）；
**本篇是工具故障导致连一篇候选论文的正文都没读到**，判词只能是**未取得、未定**，
不是「这些语料没有 η 表」。按仓库纪律
（memory `negative-existence-claim-must-grep-the-asserted-corpus`）：
没有全文时不得写「没有」，同理也不得写「有」——WebSearch 转述的具体数字**一律禁止写入任何承重表格
或消费点**，包括本篇自己列出的那几个（1910.08312 的五元素表、2308.03153 的 Nb 值、0111385 的 Nb/氢值）。

## 3. 优先级建议（供 WebFetch 恢复后直接执行，不必重新搜索）

1. **npj Comput. Mater. 2026 RMTA 论文**（`s41524-026-02141-7` / `arXiv:2508.21281`）——最高优先级：
   全新语料（今年发表）、明确是 η 表本身而非过路引用、且是首要候选中最不可能落空的一个
   （论文标题即声称对 V/Mo/Nb 及化合物给出 McMillan–Hopfield 参数）。GitHub `danylo-radevych/pp-rmta`
   若含 worked-example 输出文件，可能比论文正文更快拿到干净数字。
2. **Savrasov & Savrasov 1996**——唯一一个不落在已知封锁域名（arxiv/nature）上的直接 PDF 源
   （`savrasov.physics.ucdavis.edu`），值得优先测试 WebFetch 是否已恢复。
3. **Papaconstantopoulos et al. PRB 15, 4221 (1977)**——若 WebSearch 转述属实（32 个元素），
   量级上远超本仓已有的任何候选表，但**存在性本身未核实**，需先确认该文献真实存在、
   DOI/卷期正确，再决定是否值得走 harvest.aps.org 扫描件路线。
4. Gaspari–Gyorffy 1972 / McMillan 1968 原始论文——理论重要性最高但**当年论文极可能没有
   大规模数值表**（McMillan 1968 的摘要用词是"a number of metals"，量级未知，需实读才能判断
   够不够格当"数值表"）。

## 4. 与团队的分工请求

已 `SendMessage` 给 main 报告 WebFetch 故障、询问是否为本会话独有、
以及是否有 peer 持有可用的 WebFetch 可代为取件。本篇写作时尚未收到回复。

---

**未生成任何 `_scratch_ir/` 脚本或产物**（没有可供脚本处理的原始数据）。
**相关**：[`eta_numerical_tables_not_published_2026-09-25.md`](eta_numerical_tables_not_published_2026-09-25.md)（第一轮，读了全文的干净否定）、
`../ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md`（Tier 0/1 协议本体）。
