# ω₂ / ⟨ω²⟩ 数据源调研（2026-09-17）

调研 agent: researcher（sonnet）。任务来自 team-lead 的两条独立线：(1) Cerqueira 数据集五问；(2) 四个氢化物体系 (λ, ω_log, ω₂) 的有/无核查。
**纪律**：每条结论标注证据等级。`[L1-PDF]` = 我直接读了 PDF 页面图像并看到数字/表格。`[L1-TXT]` = 我直接读了 pdftotext 全文并看到数字/表格。`[WebFetch]` = 工具对网页内容做了转述/摘要，我未见原始逐字文本，证据等级低于前两者。`[UNVERIFIED]` = 提到但未亲眼确认，或未检索。**禁止填空/估算**——找不到的一律写"未找到"。

---

## Part 1：Cerqueira et al. 数据集五问

主论文：Cerqueira, Sanna, Marques, "Sampling the Whole Materials Space for Conventional Superconducting Materials," arXiv:2307.10728, published as *Adv. Mater.* **36**, 2307085 (2024)。venue/卷期已通过companion paper引用列表核实 `[L1-TXT]`（见下）。

### Q1：数据集托管在哪？

`[WebFetch, API 直查]` Materials Cloud Archive：
- 记录页：**https://archive.materialscloud.org/records/3kbt5-r3n56**
- DOI：**https://doi.org/10.24435/materialscloud:qv-bq**
- 标题："Sampling the materials space for conventional superconducting compounds"
- 作者：Cerqueira (Coimbra), Sanna (MPI Halle), Marques (Ruhr Bochum)；发布日期 2023-10-30；License CC-BY-4.0。

我通过 `archive.materialscloud.org` 的 REST API（`/api/records?q=Cerqueira` 及 `/api/records/3kbt5-r3n56`）独立验证了这条记录存在、可访问，字段与仓库内既有调研笔记 `surveys/a2F_machine_readable_sources_2026-09-08.md` 第14行的记录完全一致（该笔记独立地、更早给出了同一 URL），构成**双重交叉验证**。**结论：SOLVED，高置信度。**

论文正文（arXiv v1）Section VII "Data Availability Statement" 原句 `[L1-TXT]`："All data used in or resulting from this work will be available in Materials Cloud and in NOMAD as soon as the manuscript is accepted for publication." —— 正文本身不给具体 URL；NOMAD 端我没有找到对应记录（未检索到，`[UNVERIFIED]`，只确认了 Materials Cloud 一侧）。

### Q2：数据集是否包含逐材料 ω₂（或 ⟨ω²⟩），还是只有 λ 和 ω_log？—— 全任务的承重点

**答案取决于问的是哪个数据集，必须分开回答：**

**(a) 主数据集（6912 材料，arXiv:2307.10728 / Materials Cloud 3kbt5-r3n56）：**
`[WebFetch]` 通过 `/api/records/3kbt5-r3n56/files/files_description.md/content`（经 302 重定向到 `rgw.cscs.ch` 签名 S3 URL 后取得内容，工具做了转述而非逐字引用，证据等级为 WebFetch 级）：每个 batch 内每个材料的文件包括 `a2F.dos*`（不同展宽的原始 α²F(ω) 谱函数）、`lambda`（不同展宽下的 λ）、`McMillan.dat`（McMillan/Allen-Dynes 公式给出的 Tc）、`Eliashberg.dat`（不同 μ* 下 Eliashberg 方程给出的 Tc）等。**没有一个文件被描述为直接存放提取出的 ω_log 或 ω₂ 数值。** ω_log/ω₂ 只是内部计算 McMillan.dat/Eliashberg.dat 所需的中间量，是否被单独写出为可检索标量**未确认**。

论文正文 Table I（`[L1-TXT]`，已读完整文本）只列 Formula/Spg/NSites/E_hull/**ω_log**/λ/Tc，**没有 ω₂ 列**。

**结论 (a)：主数据集"没有"直接给出逐材料 ω₂——只有 λ 和 ω_log 被明确存在/发表；ω₂ 若要拿到，需要下载原始 `a2F.dos` 文件后自己按论文 Eq.(1c) 积分。我没有下载任何一个 ~2.4-3.4 GB 的 batch 压缩包去验证 α²F 文件里是否连 ω_log 都不需重新积分——这是本次调研的诚实局限，未做，非"确认不存在"。**

**(b) 姊妹/后续数据集——常压氢化物专项（arXiv:2403.13496，见 Q5）：**
`[L1-TXT]`（已读全文）Table I 标题原句："Superconducting transition temperatures calculated with the Allen-Dynes formula (T_c^AD in K), logarithmic average of the phonon frequencies (ω_log in K), **average of the square of the phonon frequencies (⟨ω²⟩ in K)**, electron-phonon coupling constant (λ) and distance to the convex hull of stability (E_hull in meV/atom)."

表中约50个化合物逐条给出 Tc_AD、ω_log、⟨ω²⟩、λ、E_hull 数值，例如：Li₂CuH₆: 86.0, 491, 890, 1.96; Mg₂PtH₆: 78.3, 770, 922, 1.24; KInH₃: 72.9, 487, 698, 1.70; PdH: 34.9, 343, 429, 1.24; ZrH₃: 25.8, 363, 748, 0.97（顺序：Tc_AD, ω_log, ⟨ω²⟩, λ）。

**结论 (b)：这个姊妹数据集"有"⟨ω²⟩，且是逐材料发表在论文正文表格里的，我亲眼在 pdftotext 全文中看到了数字。但这是"论文里的表格数字"，不是一个独立的可下载结构化数据文件**——该论文 Data Availability Statement（`[L1-TXT]`）原句："All relevant data is present in the manuscript and in the Supplementary Information." 没有提到 Materials Cloud/NOMAD/Zenodo 等外部仓库。

**总结回答 Q2**：主 6912 材料数据集——只有 λ/ω_log 被确认可得，ω₂ 未确认存在于可下载文件中（需要从原始 α²F 积分）。约50材料的常压氢化物姊妹论文——ω₂ 逐材料存在，但只在论文表格里，不是独立下载文件。

### Q3：文件格式与体积

`[WebFetch]` 主数据集：raw Quantum ESPRESSO 逐材料输出（`.in/.out`、`a2F.dos*`、`lambda`、`McMillan.dat`、`Eliashberg.dat`、`geo_opt.cif`、`el-bs.pdf`、`ph-bs.pdf`、`qe.dyn*`）打包为 17 个 `batch-a.tar.bz2` … `batch-q.tar.bz2`，每包约 500 材料，单包 2.4–3.4 GiB，**总计约 48.1 GiB**。另有 `README.txt`（1.6 KiB）与 `files_description.md`（1.5 KiB）。

### Q4：许可与引用要求

`[WebFetch, API 直查]` License: **CC-BY-4.0**。引用应指向已发表版本：Cerqueira, Sanna, Marques, *Adv. Mater.* **36**, 2307085 (2024)（该卷期号是从姊妹论文 arXiv:2403.13496 的参考文献列表 `[L1-TXT]` 中读到的：ref [44] = "T. F. T. Cerqueira, A. Sanna, and M. A. L. Marques, Sampling the materials space for conventional superconducting compounds, Adv. Mater. 36, 2307085 (2024)"），以及 DOI 10.24435/materialscloud:qv-bq。已发表 Wiley 版本页面本身（`onlinelibrary.wiley.com/doi/10.1002/adma.202307085`）被我尝试 WebFetch，**返回 403 Forbidden（付费墙）**，未能直接核对其 Data Availability Statement 文字。

### Q5：同组后续/配套数据集

`[L1-TXT]`（全文已读）**arXiv:2403.13496** — Cerqueira, Fang, Errea, Sanna, Marques, "Searching Materials Space for Hydride Superconductors at Ambient Pressure" (2024)。这正是常压氢化物专项数据集，见 Q2(b)。**无独立下载 URL**，数据只在论文正文+SI 里。

`[WebFetch]` 另外两个通过 Materials Cloud API 搜索到、但**未深入核实**的同组相关记录（`[UNVERIFIED]`，仅看到标题/DOI，未确认内容是否含 EPC/ω 数据）：
- `hvq5r-dby55`："Large-scale machine-learning-assisted exploration of the whole materials space"，DOI `10.24435/materialscloud:m7-50`
- `52jw6-bh071`："A new dataset of 415k stable and metastable materials calculated with the PBEsol and SCAN functionals"，DOI `10.24435/materialscloud:5j-9m`（看标题像是底层结构数据库 Alexandria，不像是 EPC 专项数据）
- `2hg1p-63427`："Superconductivity in antiperovskites"，DOI `10.24435/materialscloud:6g-az`（同组 Hoffmann/Cerqueira/Marques 的另一个专项数据集，非氢化物）

---

## Part 2：四个氢化物体系的 (λ, ω_log, ω₂) 核查

### H3S / D3S
来源：Errea, Calandra, Pickard, Nelson, Needs, Li, Liu, Zhang, Ma, Mauri, "Hydrogen Sulfide at High Pressure: A Strongly-Anharmonic Phonon-Mediated Superconductor," PRL **114**, 157004 (2015), arXiv:1502.02832。`[L1-PDF]` 已读全文页面图像（1–5页，全文）。

Table I (p.3) 逐字数据：

| 体系 | 压力 | λ^har | ω_log^har (meV) | λ^anh | ω_log^anh (meV) | Tc^ME,har | Tc^ME,anh | Tc(Exp) |
|---|---|---|---|---|---|---|---|---|
| H3S | 200 GPa | 2.64 | 90.1 | 1.84 | 92.86 | 250 | 194.0 | 190 |
| H3S | 250 GPa | 1.96 | 109.1 | 1.71 | 101.3 | — | 226 | 190 |
| D3S | 200 GPa | 2.64 | 68.5 | 1.87 | 73.3 | 183 | 152.0 | 90 |

（μ*=0.16；频率单位 meV，Tc 单位 K）

- **λ**：有（谐波+非谐两套，200 GPa 与 250 GPa 两个压力点）。
- **ω_log**：有（同上）。
- **ω₂ / ⟨ω²⟩**：**表I中没有此列**。我通读了全文5页正文，没有发现任何 ω₂/⟨ω²⟩ 的数值。若需要，理论上可从文中提到的 α²F(ω) 计算结果反推积分得到，但论文正文本身不给出该数值——**结论：ω₂ 需从 α²F 积分得到，本文未直接给出**。
- **补充材料**：本文引用了 Supplemental Material [ref 24]，其中是否含 ω₂ 我**没有单独取得并核查**（`[UNVERIFIED]`）。
- **异号系数/其他**：D3S 非谐同位素系数 α≈0.35（200 GPa），偏离 BCS 的 0.5。
- **"2016 年"H3S 论文**：任务简述提到"Errea 2015/2016"，我检索到一个可能对应的后续文献 Errea et al., PRB **92**, 205125 (2016)（`http://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.92.205125/fulltext`），但**未打开阅读，未核实其内容**——标记为 `[UNVERIFIED]`，本次未能确认该文献是否补充了 ω₂ 数据。

### LaH10 / LaD10
来源：Errea, Belli, Monacelli, Sanna, Koretsune, Tadano, Bianco, Calandra, Arita, Mauri, Flores-Livas, "Quantum Crystal Structure in the 250 K Superconducting Lanthanum Hydride," *Nature* **578**, 66 (2020), arXiv:1907.11916。`[L1-PDF]` 已读全部 14 页正文+Methods+Extended Data 部分（pages 1–14 of 20；pages 15–20 为更多 Extended Data 图未读，见下方局限说明）。

**Extended Data Table I**（p.11，"Summary of calculated superconducting Tc"）逐字数据：

| System | Pressure (GPa) | λ | ω_log (meV) | Tc^Mc_μ*=0.1 (K) | Tc^AD_μ*=0.1 (K) | Tc^ME_ani (K) | Tc^SCDFT (K) |
|---|---|---|---|---|---|---|---|
| LaH₁₀ | 129 | 3.62 | 76.4 | 171.8 | 252.6 | 255.3 | 230 |
| LaH₁₀ | 163 | 2.67 | 96.4 | 197.1 | 247.0 | 242.8 | 225 |
| LaH₁₀ | 214 | 2.06 | 115.5 | 196.3 | 235.9 | 237.9 | 210 |
| LaH₁₀ | 264 | 1.73 | 126.6 | 189.5 | 219.2 | 216.9 | 201 |
| LaD₁₀ | 159 | 3.14 | 63.5 | 135.0 | 184.2 | 180.4 | 171 |
| LaD₁₀ | 210 | 2.21 | 81.7 | 145.5 | 176.5 | 172.9 | 158 |
| LaD₁₀ | 260 | 1.80 | 92.2 | 142.2 | 164.6 | 157.9 | 151 |

（μ*=0.1；此表用的声子频率来自 **SSCHA 非谐计算**，Tc 同时给出 McMillan/Allen-Dynes-McMillan/各向异性 Migdal-Eliashberg/SCDFT 四种方法对比。）

- **λ**：有（7个压力点，Fm-3m 相，非谐 SSCHA 声子）。
- **ω_log**：有（同表，单位 meV）。
- **ω₂ / ⟨ω²⟩**：**该表没有此列**。我逐页读了 Methods 部分及 Extended Data Figures 1–6（含 α²F(ω) 谱函数图 Extended Data Fig. 4、5，给出了 129/163/214/264 GPa 的 α²F(ω) 曲线图及累积 λ(ω) 曲线，但没有印出 ω₂ 数值）。**结论：ω₂ 需从图中 α²F(ω) 积分得到，本文未直接给出数值**。
- **压力与非谐**：129–264 GPa（LaH10），159–260 GPa（LaD10）；**非谐效应已包含**（全文核心方法即 SSCHA 量子离子效应）。
- **未读部分**：本 PDF 共 20 页，我读到第 14 页（Extended Data Fig. 6）。第 15–20 页（推测为更多 Extended Data 图与表，如结构相关的 Extended Data Table II）**未读**，标记 `[UNVERIFIED]`——但已经拿到的 Extended Data Table I 是最可能包含 λ/ω_log 汇总的地方，且已确认无 ω₂ 列，因此后续页面大概率不影响"ω₂ 未直接给出"这一结论，但我不能 100% 排除。

### YH6 / YH9
来源：Ding, Chen, Shi, "Anharmonicity and Coulomb pseudopotential effects on superconductivity in YH6 and YH9," arXiv:2507.03383v2 (2025)。`[L1-PDF]` 已读全文（7页正文全部）。

**Table I**（p.4）逐字数据：

| System | Method | λ(0) | λ(1) | λ(2) | ω̄₂ (K) | Tc (K) |
|---|---|---|---|---|---|---|
| YH6 | DFPT | 2.534 | 1.109 | 0.477 | 1399 | 270 (248) |
| YH6 | SSCHA | 1.779 | 0.948 | 0.446 | 1662 | 245 (218) |
| YH6 | SPIA | 1.712 | 0.910 | 0.444 | 1738 | 212 |
| YH9 | DFPT | 3.038 | 0.958 | 0.409 | 1308 | 252 (236) |
| YH9 | SPIA | 1.749 | 0.944 | 0.460 | 1905 | 240 |

（μ*=0.11；括号内 Tc 为用重整化 μ*(N) 后的值；DFPT/SSCHA 谱函数来自外部引用 [11]（YH6）与 [20]（YH9），本文作者在此基础上重新计算了各阶矩。）

> **★ 2026-09-25 订正（主会话一手核对参考文献表）——原文此处把 [11] 写成 Kong，错了。**
> **旧陈述错在哪**：原写 `[11]=Kong et al. Nat. Commun. 12, 5075 (2021)（YH6）`。
> **引用编号 [11] 是对的，挂在它上面的作者名是错的**——Kong 是 **[12]**。
> **新陈述凭什么对**（`_scratch_ir/_ding2507_txt.txt` 参考文献表逐字）：
> - 第 404 行：`[11] I. A. Troyan, D. V. Semenok, A. G. Kvashnin, ... and A. R. Oganov, Anomalous high-`
>   ⟹ **[11] = Troyan et al. 2021（YH6）**；该文作者含 Errea / Bianco / Calandra / Mauri / Monacelli
>   ——正是做 SSCHA 的那批人，与「SSCHA 谱来自 [11]」自洽。
> - 第 347–352 行：`[12] P. Kong, V. S. Minkov, ... and M. I. Eremets, Superconductivity up to 243 K in the`
>   `yttrium-hydrogen system under high pressure, Nature Communications 12, 5075 (2021)`
>   ⟹ **Kong 是 [12]**，且它是 **YH9 的 243 K 实验**文献，不是 YH6 的谱源。
> - 正文归属逐字：第 228 行 `parameters of YH6 calculated by SSCHA in Ref. [11]`；
>   第 260 行 `SSCHA result of YH6 from Ref. [11]`。
> - 第 384 行：`[20] M. Du, Z. Li, D. Duan, and T. Cui, Superconducting ...`
>   ⟹ **[20] 的作者是 Du–Li–Duan–Cui 四人**，原文缩写成「Du & Duan」漏了两位。
>
> ⟹ **YH6 的谱源是 Troyan 2021，不是 Kong 2021。**
> 下游影响：Tc_exp 的取值归属亦按此——**224 K 出自 Troyan=[11]（YH6）**、
> **243 K 出自 Kong=[12] 的标题（YH9）**。见 `../ceiling/hydride_ceiling_r_h3s_yh6_yh9_2026-09-25.md`。

- **λ**：有——但注意本文用的是"λ(m)"（Matsubara阶数 m 的推广定义），不是传统单一标量 λ；λ(0) 对应常规 λ。
- **ω_log**：**Table I 中没有 ω_log 列**。我通读了全文7页，正文中没有找到显式的 ω_log 数值报告——本文用的核心量是 ω̄₂（见下），文中未见 ω_log 的单独讨论或表格。**结论：ω_log 在本文中未找到（not found）**，若需要须自己从文中引用的外部 α²F(ω)（Kong 2021 / Du 2023）积分得到。
- **ω₂ / ⟨ω²⟩（本文记为 ω̄₂）**：**有**，Table I 直接给出（单位 K），定义为 ω̄₂(m→∞) = 2π/(ħβ)·√(m²λ(m)/λ(0)) 的渐近值，与标准 Allen-Dynes 二阶矩定义在物理含义上一致。
- **压力与非谐**：YH6 在 165 GPa；YH9 在 255 GPa（DFPT 结果为260 GPa，与正文255 GPa略有不同，原文脚注注明）。**非谐效应已包含**——通过非微扰 SPIA（stochastic path-integral approach）以及引用的 SSCHA 结果对比。
- **重要说明（证据等级降级）**：Table I 中标注"DFPT"和"SSCHA"两行的 ω̄₂/λ 数值，是 Ding et al. **依据 Kong 2021 / Du 2023 发表的 α²F(ω) 谱函数重新计算**得到的（本文自陈），而不一定是 Kong 2021 / Du 2023 原文自己印出的数字。Kong 2021 / Du 2023 原文本身是否直接报告 ω₂——**我没有检查这两篇原始文献，未核实（`[UNVERIFIED]`）**。

### CaH6
`[WebFetch]`（非直接 PDF 阅读，证据等级低于前三个体系）—— 本地库中**没有 CaH6 的 PDF**（已用 Glob 搜索 `library/pdf/**/*Ca*H6*.pdf`，无结果）。经 WebSearch 定位候选文献后，对 arXiv:2111.10797（Liang, Che, et al. 或相近作者，"Electron-phonon coupling and superconductivity in an alkaline earth hydride CaH6 at high pressures"）做了一次 WebFetch（读取 `arxiv.org/html/2111.10797` 的 HTML 全文渲染版，工具做了转述，非我逐字亲见原始 PDF）：

- **λ 与 ω_log**（按 WebFetch 转述的原文句子拼合）：
  - 150 GPa：λ = 2.75，ω_log = 684 cm⁻¹
  - 200 GPa：λ = 1.99，ω_log = 851 cm⁻¹
  - 250 GPa：λ = 1.71，ω_log = 933 cm⁻¹
  - 另有 120 GPa 处 λ 极大值 4.08（该点未给出 ω_log）
  - 注意：172 GPa（对应 Ma et al. 实验压力点）本文**没有**直接给出数值。
- **ω₂ / ⟨ω²⟩**：据 WebFetch 转述，全文检索"second moment"/"omega_2"/"average frequency"等词**均无匹配**，即**未找到**。
- **非谐效应**：据 WebFetch 转述，本文方法为标准 DFPT（谐波），**未提及** SSCHA/PIMD 等非谐处理。

**重要局限声明**：以上 CaH6 数据**全部经由 WebFetch 网页摘要获得，我没有亲自打开/核对 PDF 原文或逐字确认这些数字**，证据等级为本报告中最弱的一档（`[WebFetch]`，非 `[L1-PDF]`）。若该结论需要用于承重计算，**必须**后续直接读取 arXiv:2111.10797 的 PDF 页面加以核实。此外，可能还有更晚、更权威的 CaH6 文献（如 Darussalam & Koretsune 2024, *Sci. Rep.*，在另一篇 arXiv:2512.08126 的引用中提到——`[UNVERIFIED]`，我没有找到并打开该文献）未被本次调研覆盖。

---

## 汇总表

| 体系 | λ | ω_log | ω₂/⟨ω²⟩ | 压力 | 非谐 | 证据等级 |
|---|---|---|---|---|---|---|
| H3S/D3S | 有 | 有 | 需从α²F积分（本文未给） | 200/250 GPa | 是 | L1-PDF |
| LaH10/LaD10 | 有 | 有 | 需从α²F积分（本文未给，图给了曲线未给数值） | 129–264 GPa | 是 | L1-PDF（部分页未读，见上） |
| YH6/YH9 | 有 | **未找到** | 有（ω̄₂，Table I） | 165/255 GPa | 是 | L1-PDF（外部谱函数来源未核实） |
| CaH6 | 有 | 有 | 未找到 | 120–250 GPa | 否（未提及） | WebFetch（未直接核实，最弱） |

## Cerqueira 数据集五问汇总

| 问题 | 结论 | 证据等级 |
|---|---|---|
| Q1 托管处 | Materials Cloud, doi:10.24435/materialscloud:qv-bq，记录 3kbt5-r3n56 | WebFetch+API 双重交叉验证 |
| Q2 是否有逐材料ω₂ | 主6912材料集：**没有**（只有λ/ω_log，ω₂需自行从原始α²F积分，未下载验证）；姊妹常压氢化物专项集（arXiv:2403.13496,~50化合物）：**有**，但只在论文表格里非独立文件 | L1-TXT（正文表格）+ WebFetch（files_description.md） |
| Q3 格式/体积 | 17个.tar.bz2，共约48.1 GiB，内含raw QE输出 | WebFetch |
| Q4 许可 | CC-BY-4.0；引用 Adv. Mater. 36, 2307085 (2024) | WebFetch+L1-TXT(引用列表) |
| Q5 姊妹数据集 | arXiv:2403.13496（常压氢化物专项，有ω₂但数据只在论文里非独立下载）；另两个未核实的Materials Cloud记录见正文 | L1-TXT + WebFetch |

## 未完成/局限事项（诚实清单）

1. 没有下载任何一个 Materials Cloud batch 压缩包（~2.4–3.4 GB/个）去确认 `a2F.dos` 原始文件里 ω_log/ω₂ 是否已作为中间量被写出——只能确认 `files_description.md` 没有明确提到独立的 ω_log/ω₂ 文件。
2. LaH10 论文 PDF 第 15–20 页（部分 Extended Data）未读。
3. H3S "2016" PRB 92, 205125 后续论文完全未打开阅读。
4. YH6/YH9 论文引用的两篇原始谱函数来源（Kong 2021, Du 2023）未核实是否自己报告 ω₂。
5. CaH6 全部结论仅基于一次 WebFetch 网页摘要，未直接读取 PDF 原文核实；也未找到更权威/更完整的 CaH6 专项文献（如可能存在的 Darussalam & Koretsune 2024 Sci. Rep.）。
6. Materials Cloud 另外两个候选记录（`hvq5r-dby55`, `52jw6-bh071`）未核实内容与本任务的相关性。
