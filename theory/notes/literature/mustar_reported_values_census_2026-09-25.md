# 氢化物 Eliashberg/SCDFT 计算中 μ\* 取值：有界文献普查（2026-09-25）

> 消费点：[`../ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md`](../ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md) §6 第 6 条——收紧 0.10/0.13 这档惯例区间的锐度，需要一条独立的 μ\* 下界，目前未有。本笔记只回答「文献里实际写了什么」，**不判定**该不该取 0.10。

## ① 结论（一句话 + 等级）

在本次读到的 9 篇一手文献里，氢化物超导 Eliashberg/SCDFT 计算实际使用（INPUT）或反解（FITTED）的 μ\* 覆盖 **0.085–0.35**；「0.10–0.13」是这批文献里 INPUT 类的常见下段，不是任何一篇给出的下界或上界。**等级：numerical-evidence（文献计数，N=9 篇一手 + 2 条二手转引未核实）；不是界。** 未见任何一篇文献给出 μ\* 的严格上/下界陈述（见 ③）；Errea 2020（LaH10）与 Ding et al. 2025（YH6/YH9）显示业内已注意到 μ\* 数值本身依赖 ω_c 截断约定，并非普适常数。

## ② 逐篇表

| 论文 | 材料 | 压力 | μ\* | INPUT / FITTED / 计算(无 μ\*) | 出处（表/式/页） | 仪器 |
|---|---|---|---|---|---|---|
| Klein, PRL 39, 574 (1977)（未在本轮重读，转引自 researcher agent memory `lit_pdh-mustar-klein-essenberger-source-hunt-2026-09-18`） | PdH | 常压 | **0.085** | INPUT（Bennemann–Garland 公式 `μ*(x)=0.26N(E_F;x)/[1+N(E_F;x)]` 算出，非拟合 Tc） | 原文式(见 memory 转述，未本轮重读原文，标"转引未核实原文") | 二手 |
| Errea et al., PRL 114, 157004 (2015)（H3S 首篇） | H3S / D3S | 155–250 GPa | **0.16** | INPUT（正文明写 "Tc=250K for μ\*=0.16"（谐波高估）、"account for the experimental Tc when μ\*=0.16 is used"（非谐）） | Table I 图注 + 正文，p.3 | 单仪器（页图，Read pages=，无独立文本层核对） |
| Errea et al., PRL 116/arXiv:1512.02933 (2016, "Quantum Hydrogen-Bond Symmetrization…") | H3S / D3S | 155 GPa 一带 | **0.16** | INPUT（Methods 明写 "the isotropic Migdal-Eliashberg equations were solved using 512 Matsubara frequencies and μ\*=0.16"） | 正文 Methods，p.4 | 单仪器（页图） |
| Szczesniak & Durajski, Solid State Commun. 249, 30 (2017)（H3S 同位素） | H3S | 130 GPa / 200 GPa | **0.239（130 GPa 谐波）/ 0.286（200 GPa 谐波）/ 0.146（非谐）** | **FITTED**（正文明写"μ\* corresponding to [Tc]exp were equal to…"，反解自实验 Tc） | 正文 p.2，ω_C=3Ω_max | 单仪器（页图） |
| Durajski, arXiv:1912.09663 (2020, LaH10 谐波同位素) | LaH10 | 250/300/350 GPa | 扫描 **0.10 / 0.15 / 0.20**；正文倾向 0.20（"better agreement with the measurements"） | INPUT 扫描（Table I 三值全给），**0.20 一档带准拟合倾向**（比较实验后择优，非严格反解单值） | Table I, 正文 p.2/p.4/p.5 | 双仪器（本轮读到文本层+页图） |
| Errea et al., Nature 578, 66 (2020)（LaH10） | LaH10 | ~150–250 GPa | **无 μ\***（ME 用 ab initio 静态 RPA 屏蔽库仑 V^C，SCDFT 同样"without an empirical μ\* parameter"）；仅在与 Allen-Dynes 对比时用 **μ\*=0.1** 作图示（非主结果） | **计算（无自由参数）**，Allen-Dynes 对比项为 INPUT | p.4/p.7/p.8，Extended Data Table I | 单仪器（页图，此前会话已读） |
| Villa-Cortés & De la Peña-Seaman, arXiv:2103.12055 (2022, H3S/D3S 范霍夫奇点) | H3S / D3S | 162–250 GPa | **FITTED 逐压力点**：Fig.7 显示 N(ε) 方案 H3S≈0.24→0.23、D3S 略低；N(0) 方案 H3S≈0.15→0.10、D3S≈0.11→~0.08（**图读，未标数值坐标刻度精度，视为近似**） | **FITTED**（正文明写"μ\* is fitted at each pressure of interest"，用已知 Tc(P₀) 反解，Sec. II/III） | Fig.7 + 正文 p.4 "the only adjustable parameter left is μ\*, which is fitted at each pressure of interest" | 单仪器（页图；数值系图读估计，未交叉） |
| Troyan et al., arXiv:1908.01534（YH6 实验+理论综述段） | YH6 | 120–166 GPa | 转引 Li et al. 2015：**μ\*=0.1–0.13**（Introduction 段引用，非本文自算值） | INPUT（引自他文，本文未独立算 μ\*） | 正文 Introduction p.2 | 单仪器（页图） |
| Ding, Chen & Shi, arXiv:2507.03383 (2025, YH6/YH9 非谐+μ\*截止重整化) | YH6 / YH9 | 165 GPa / 255 GPa | **0.11**（"Setting μ\*=0.11, as determined from GW methods[33,36]"） | INPUT（GW 方法算出，非拟合实验 Tc；但作者强调此值须经 Eq.(4) 的截止重整化 μ\*(N) 才能给出收敛 Tc） | Table I 图注 + 正文 p.4，"IV. C. Renormalization" | 单仪器（页图；本论文文本混排清晰，但仍未独立跑 pdftotext） |

**转引未核实的二手线索（不计入上表严格计数）**：Durajski–Szczesniak, Sci. Rep. 6, 38570 (2016)（PH3，报告 μ\*≈0.083–0.088 的拟合值，来自既往 agent memory 检索摘要，本轮未重新取得原文）；Akashi et al., Sci. Rep. (2018)/arXiv:1804.09760（H5S2，μ\*≥0.402，同样转引自既往 memory，未重读）。

## ③ INPUT vs FITTED 分类计数

- **INPUT**（取定值，非反解实验 Tc）：Klein 0.085；Errea 2015 0.16；Errea 2016 0.16；Durajski 2020 扫描 0.10/0.15/0.20；Errea 2020 对比项 0.1；Troyan 转引 0.1–0.13；Ding et al. 2025 0.11（GW 算出）。**N=7 篇/条目**，取值范围 0.085–0.20，中位数落在 0.11–0.16 一带。
- **FITTED**（反解自实验 Tc，含循环信息，不能支持独立下界）：Szczesniak-Durajski 2017（0.146–0.286）；Villa-Cortés 2022（图读 0.08–0.28，视方案而定）。**N=2 篇**，全部 ≥0.08，且都**高于**部分 INPUT 类的下段（0.085–0.10），说明 FITTED 并不系统性地把 μ\* 推向更低。
- **计算/无 μ\* 参数**：Errea 2020（LaH10 ab initio RPA 屏蔽库仑）。**N=1 篇**，本身不产出一个可比较的 μ\* 数字。

**在这 9 个条目里，最小值是 Klein 1977 的 0.085（PdH，转引未核实原文，且 PdH 是低-λ 常压体系，与高压强耦合氢化物不同族），其次是 Durajski/Ding/Errea 等主流高压氢化物 INPUT 类的 0.10–0.11。FITTED 类从未低于约 0.08（Villa-Cortés 图读下限，且该图读精度未交叉核实）。**

## ④ 语料范围与检索方式（可复现）

- **一手全文读取**：`library/pdf/P1/` 与 `library/pdf/D3/` 目录下本地 PDF，逐篇用 `Read(pages=...)` 抽取（本会话读了 Errea 1502.02832、1512.02933、Durajski 1912.09663、Villa-Cortés 2103.12055、Troyan 1908.01534、Ding 2507.03383；此前会话已读 Errea 1907.11916、Szczesniak-Durajski、Sadovskii 1811.10184[排除，非氢化物]）。**除 Durajski 1912.09663 外，全部是单仪器（仅页图，无独立 `pdftotext -enc UTF-8` 文本层交叉核对，本 agent 无 Bash 工具无法执行该指令）。**
- **未读的本地 PDF**（存在于 `library/pdf/P1/` 或 `D3/`，本轮受限于有界普查未打开）：`P1/1506.08190, 1812.01561, 1808.07695, 2103.12055`（已读）`, 1909.10482, 2101.01315, 1802.07968, 2307.15684, 1509.03156, 2510.25902, 2608.14949`；`D3/2406.13178, 2501.14420, 2212.09789, 2512.24849, 2602.00833, 2211.14474, 2508.19781`；`D1/2507.01897, 1911.12872`；`D2/2503.03559`；`D4/2307.10728, 2403.13496, 2502.18281`。这些文件**未取得内容，未判定**其 μ\* 用法。
- **纯转引/二手**（未本轮取得原文）：Klein 1977 PRL 39,574（依据既有 agent memory 转述）；Durajski–Szczesniak Sci. Rep. 6,38570 (2016)；Akashi et al. Sci. Rep. (2018)/arXiv:1804.09760。
- 未执行 WebSearch/WebFetch 补充检索（`arxiv.org`/`www.nature.com` 本机 TLS 劫持，见既有 memory `arxiv-and-nature-are-tls-intercepted-on-this-machine`；本轮未尝试 `harvest.aps.org` 或 `zenodo.org` 路线，属未做而非做了取不到）。

## ⑤ 未取得清单

- Klein 1977 (PRL 39, 574) 原文——未取得，转引自既有 agent memory，等级降为"二手、未核实"。
- Durajski–Szczesniak, Sci. Rep. 6, 38570 (2016)（PH3）原文——未取得。
- Akashi et al., Sci. Rep. (2018)/arXiv:1804.09760（H5S2）原文——未取得。
- `library/pdf/` 下列于 ④ 的约 22 篇未读 PDF——完全未取得，包括 Sanna 2508.19781、Mishra 2507.01897、Kogler/IsoME 2503.03559 等既有 memory 提示含 μ\* 相关讨论的文件。
- 关于「μ\* 下界」的任何直接陈述——**摘要层与全文层均未见**（基于本次读到的 9 篇一手全文；未读的 ~22 篇未纳入此否定判断，不得推广为"文献中没有"）。

---

**给主会话的判断材料，非结论**：若要收紧 §6 第 6 条的锐度，本普查显示的自然候选下界候选值是 **0.10**（多篇独立 INPUT 一致下限，且被 Ding et al. 2025 的 GW 计算值 0.11 与 Errea 2020 的对比线 0.1 双重印证），但这仍是**惯例区间的下段，不是任何作者主张的下界**；是否可将其提升为该方向的工作假设，由主会话裁定。
