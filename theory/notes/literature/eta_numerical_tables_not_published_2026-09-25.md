# Tier 0 的数据荒：Shipley 2021 与 Moussa–Cohen 2006 都不印 η 数值表

**日期**：2026-09-25 · **等级**：见每条 · **状态**：live，但 **§5「外部效度」已被同日订正框推翻**
（§1–§4 全部有效；订正见 §5 末）
**起因**：团队 agent `eta-table-harvest` 受派去两篇论文里收 (η, M) 数对，喂给 Tier 0 判据
`kTc ≤ C_∞ √(η/M)`（`C_∞ = (1/2π)√𝔤(2) = 0.1827262477`，analytic）。
**结论**：**两篇合计 0 个可用 (η, M) 数对。**

> ⚠ 本笔记的每条逐字认定都是**我自己从零重跑**的（CLAUDE.md §1），不是转引 agent。
> agent 自陈它在 MC 一侧是**单仪器、未交叉**；而它以为自己持有的两份 Shipley 文本
> （`_shipley_2026-09-25.txt` 与 `_shipley_prb104.txt`）**md5 相同**（`613cd29e…`）
> ⟹ 那是**同一个文件的两个名字**，不构成两台仪器。见 §4。

---

## 1. Moussa–Cohen 2006（PRB **74**, 094520）—— 全文零表格

本机 `library/pdf/D1/2006_Jonathan_E_Moussa_cond-mat-0607832.pdf`，12 页。
我的抽取：`pdftotext -enc UTF-8` 全文 → `_scratch_ir/_mc2006_full_mine_utf8.txt`（57180 B）。

- `grep -ic "table"` = **5**，但逐条看，**5 条全部是 `stable` / `unstable` 里的子串**
  （`s·table`），**没有任何一条是表格标题**。⟹ **全文零 `TABLE` 环境。**
- η 只以**图轴标签**出现：`~η(E) (eV/Å2)` 共 8 处（Figs. 2/3/…/10 的纵轴），
  即 η(E) 是**曲线**，没有印成数字。
- 全文只有**三个 η 标量**，且都不可消费：
  - `The experimental value of η at E_F is 2 eV/A˚2`（fcc Al，FIG. 7 图注）；
  - `we get η_max values of 27 eV/Å2 for carbon, 11 eV/Å2 for silicon, and 1.6 eV/Å2 for aluminum`
    —— 紧接下一句是作者自己的免责：`While these values are not quantitatively correct, they
    reasonable reproduce the trend`，且它们出自一个**紧束缚玩具模型**（N↑(E) 从图里读、
    E_F 取占据价带带宽）。⟹ **作者自陈非定量 ⟹ 不得当数据消费。**

**等级**：`numerical-evidence`（「某论文印了什么」的事实认定）。
**仪器**：文本层（`-enc UTF-8`）**单台**。⚠ 未跑页图。
对**否定存在命题**（「全文没有 η 表」）而言，文本层是合适仪器——表格标题必进文本层——
但按 CLAUDE.md §2 这仍是单仪器，**若将来有人要把这条当承重用，须补页图**。
当前它只用来说明「我们拿不到数据」，**不支撑任何物理判词** ⟹ 单仪器可接受。

## 2. Shipley 2021（PRB **104**, 054501）—— 只有一张表，表里没有 η

全文经 `harvest.aps.org` 取得（此前 5 次 arXiv 失败实为**本机 TLS 劫持**，见
`../ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md` §6）。

`grep "TABLE"` ⟹ **全文只有一张表**。其标题逐字：

> `TABLE I. Allen-Dynes and Eliashberg Tc values for dynamically stable superconductors found in
> this work, along with calculated λ and stability measures.`

其**八列表头逐字**（我从文本层第 265 行读出）：

| # | 列名（逐字） |
|---|---|
| 1 | `Stoichiometry` |
| 2 | `Space group` |
| 3 | `Pressure (GPa)` |
| 4 | `Allen-Dynes Tc (K)` |
| 5 | `Eliashberg Tc (K)` |
| 6 | `λ` |
| 7 | `E_stoic (meV/unit)` |
| 8 | `E_struc (meV/unit)` |

⟹ **无 η、无 M、无 N(E_F)、也无 `ω_log` / `⟨ω²⟩`。**
η 在这篇里自始至终是 **GPR 模型的内部描述符**，不作为可发布的物理量列出
（逐字见 memory `shipley-2021-uses-eta-as-gpr-descriptor-not-a-filter`，本轮复核该句仍在第 105 行：
`we do not obtain direct estimates of λ or Tc from Gaspari-Gyorffy theory and instead use the
η j values themselves as descriptors in our Tc model`）。

**等级**：`numerical-evidence`（事实认定）。**仪器**：文本层单台 + 已有 p.5 页图
`_scratch_ir/_shipley_pages/page-05.png`（此前会话生成，本轮未重渲）。

## 3. ⚠ 两个必须分开写的副产物

### 3.1 正面副产物：Shipley Table I 是一个 36+ 行的 (λ, Tc_AD, Tc_Eliashberg, P) 数据集

这对**天花板检验没用**（缺 `⟨ω²⟩` ⟹ 算不出 `M2 = λ⟨ω²⟩`），
但对 **AD-vs-Eliashberg 亏损线**可能有用：它给出**同一结构上两种方法的配对 Tc**，
样本量远大于我们现有的氢化物 9 行。作者自己也做了这个比较
（`it systematically underestimates`，并拟合出 `a = 1.0061`、`b = 0.0663`）。
**状态：未消费、未核、仅登记**。要用之前须先算预算 + 定消费点（CLAUDE.md「预算优先」）。

### 3.2 ⚠ 陷阱：不要从 Tc_AD 反解 `ω_log` 来补上缺的 `⟨ω²⟩`

表里有 λ 和 Tc_AD，看起来可以反解 `ω_log` 把缺的频率补上。**不要做。**
这正是 `hydride-csv-audit-filling-factor-circular` 记下的那个形状：
**从 AD 公式反解出来的量，再拿去检验一个含 AD 结构的判据，自由度会塌掉**——
量到的是「AD 公式 vs 我们的天花板」，不是「真实材料 vs 天花板」。
且本例连反解都做不干净：AD 的强耦合修正 `f₁f₂` **同时依赖 `ω_log` 与 `ω̄₂`**，
一个方程两个未知数 ⟹ **欠定**；再加上表里**未声明 μ\***。
⟹ 判词：**这条路不通，且即使通了也是循环的。** 两条理由独立，各自成立。

## 4. 方法学：agent 的「两份文本」其实是同一个文件

`_shipley_2026-09-25.txt` 与 `_shipley_prb104.txt` 的 md5 都是 `613cd29ed0cd4d1b5ce4b148da24dabf`
（65663 B，同一天不同时刻复制）。agent 诚实地自陈了「不确定是否带 `-enc UTF-8`」，
但**没有查它们是不是同一份**。
⟹ 教训：**「我有两份文件」不等于「我有两台仪器」**——
两台仪器的定义是**两条独立的生成路径**（文本层 vs 页图），不是两个文件名。
与 `zero-difference-at-discontinuity-is-not-cross-validation`（两套「独立」算法给 0.00e+00
相对差反而可疑）同形：**先问两条路径是不是真的独立，再看它们一不一致。**

## 5. 对 Tier 0 的净影响

**声音性不受影响**：Tier 0 判据 `kTc ≤ C_∞√(η/M)` 的成立性与有没有数据无关
（`C_∞` 是 analytic、μ\*=0 的界对任何 μ\*≥0 逐字成立，见
`../../derivations/P4_einstein_extremal.md` §2 推论 μ\*）。
**外部效度受影响**：我们至今**没有一个来自文献的真实 (η, M) 数对**可以代入。
现有落点全部走的是 `M2 = λ⟨ω²⟩` 这条路（需要声子计算），
而 Tier 0 的卖点恰恰是 `M2 = η/M` 这条**零声子计算**的路 ⟹ **卖点尚未被任何真实数据演示过**。
⟹ 报告 Tier 0 时**必须写明这一条**，不得让读者以为已有 η 侧的实例。

> ### ★ 主会话订正（2026-09-25，同日）——上面这段已失效，不得再引
>
> **旧陈述错在哪**：上一段写「我们至今没有一个来自文献的真实 (η, M) 数对」、
> 「卖点尚未被任何真实数据演示过」。这两句在**写下时**对**本笔记扫的那两篇语料**
> （Moussa–Cohen 2006、Shipley 2021）成立，但被当成了**全仓状态**陈述。
> 它把「这两篇里没有」外推成了「哪里都没有」——正是
> `../../..` 记忆 `negative-existence-claim-must-grep-the-asserted-corpus` 的形状：
> **否定存在命题的作用域，只能是你实际扫过的那个语料。**
> §5 标题「对 Tier 0 的净影响」把作用域悄悄放大到了整个 Tier 0。
>
> **新陈述凭什么对**：同日的
> [`../ceiling/tier0_eta_ceiling_32_elements_2026-09-25.md`](../ceiling/tier0_eta_ceiling_32_elements_2026-09-25.md)
> 已经拿到并消费了**两个**真实 (η, M) 语料，我（主会话）从零复核过：
> - Papaconstantopoulos et al. 1977 TABLE III（p.4224），`η_sp/η_pd/η_df/η_tot(I)/η_tot(II)`，
>   单位 **eV/Å²**，**32 个元素** ⟹ `η_tot(I)` 与 `η_tot(II)` 两档各 **0/31 违反**，
>   最紧 In `r = 0.4460 / 0.4492`；
> - Klein–Boyer–Papaconstantopoulos, **PRL 42, 530 (1979)** TABLE I（p.531），
>   `η_A`、`η_B`（eV/Å²），**12 个 A15** ⟹ **0/10 违反**，最紧 **Nb₃Ge `r = 0.6260`**。
>
> 合计 **0/41**，全部走 `M2 = η/M` 这条零声子计算的路。
> ⟹ **「卖点尚未被真实数据演示」这句已被证伪；Tier 0 的卖点已在真实发表数据上跑通。**
>
> **本笔记剩下什么仍然有效**（作用域收回到它该在的地方）：
> §1 / §2 / §3 / §4 **逐字全部不动**——「MC 2006 全文零表格」「Shipley Table I 八列无 η」
> 「不要从 Tc_AD 反解 ω_log」「两个文件名不等于两台仪器」都只陈述那两篇本身，与新语料无冲突。
> §5 第一句（**声音性不受影响**）同样不动。
> 失效的只有「外部效度」那三句的**全仓外推**。
>
> **等级**：本订正框中「两个语料存在且判决为 0/41」是 `numerical-evidence`
> （两端都是他人发表的 DFT 值 + 实验 Tc；无新 analytic 命题、无证书）。
> **⚠ 不要把 0/41 读成 41 个独立检验**：它是 31 + 10，且 A15 那 10 行的 `M2_min`
> 取法是重数方向上**最难**的一档（单边界方向规则）；`η_A` 是否为子晶格合计这第三种读法
> **尚未关闭**，见该笔记 §A15。若它成立，Nb₃Ge 的 0.626 会上升。

**下一步的唯一诚实选项**：η 的数值表要到**别的语料**去找
（Gaspari–Gyorffy 系列、McMillan 1968 的原始 η 表、以及任何印 `η` 或 `N(E_F)⟨I²⟩` 的表）。
⚠ 按 `negative-existence-claim-must-grep-the-asserted-corpus`：
在拿到全文之前，对「某某论文有没有 η 表」只能写「**未取得、未定**」，不得写「没有」。

---

**脚本 / 产物**：
- `_scratch_ir/_mc2006_full_mine_utf8.txt`（我本轮生成，`pdftotext -enc UTF-8`，12 页全文）
- `_scratch_ir/_shipley_2026-09-25.txt`（既有；md5 见 §4）
- `_scratch_ir/_shipley_pages/page-05.png`（既有页图，本轮未重渲）

**相关**：`../ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md`（Tier 0/1 协议本体）、
`einstein_extremality_and_eta_criterion_priority_2026-09-23.md`（优先权）。
