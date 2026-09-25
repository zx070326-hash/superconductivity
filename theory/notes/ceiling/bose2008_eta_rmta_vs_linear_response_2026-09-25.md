# Bose 2008 (hcp Sc)：**第一篇通过验收门的 RMTA-vs-非RMTA η 对照**

> 2026-09-25，主会话从零取件 + 从零复核。
> 复核脚本：[`_scratch_ir/_chk_main_bose2008_eta_rmta_vs_linear_response_2026-09-25.py`](../_scratch_ir/_chk_main_bose2008_eta_rmta_vs_linear_response_2026-09-25.py)（`ALL CHECKS PASSED`）
> 相关：[`radevych2025_external_witness_pap1977_eta_2026-09-25.md`](radevych2025_external_witness_pap1977_eta_2026-09-25.md)（缺口的来历与三起「实现一致≠方法正确」）、
> [`savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md`](savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md)、
> [`tier0_eta_ceiling_32_elements_2026-09-25.md`](tier0_eta_ceiling_32_elements_2026-09-25.md)

## 0. 这条为什么重要

Tier-0 天花板 `kTc ≤ C_∞√(η/M)` 的全部 η 来自 **RMTA**（Gaspari–Gyorffy 刚性 muffin-tin）。
这是整条链**唯一未经检验的输入**。此前该缺口的状态是：

- 唯一干净的方法偏差数只有 **sp 金属 Al** 的 `η_RMTA/η_true = 0.14/0.44 = 0.318`（Savrasov PRB 54,16487 同篇同代码）；
- **过渡金属侧被量化的跨方法 η 对照个数 = 0**——只有 `similar` / `close` / `complete agreement` 三句散文，
  且其中两句经查比的是 **RMTA vs 更早的 RMTA**（实现 vs 实现）。

**本文把这个计数从 0 变成 1。**

## 1. 文献与取件

> S. K. Bose, *Linear response results for phonons and electron-phonon coupling in hcp Sc —
> spin fluctuations and implications for superconductivity*,
> J. Phys.: Condens. Matter **20**, 045209 (2008)；预印本 **arXiv:0711.1543v1**。

取件：`curl https://arxiv.org/pdf/0711.1543` → HTTP 200、314 106 B、11 页。**全程未用 `-k`。**

⚠ **仪器记忆更新**：`arxiv.org` 本机现由 `Certainly / Certainly Intermediate R1`（arXiv 真实 CA，Fastly）签发，
python `ssl` 默认上下文**校验通过**。记忆 `arxiv-and-nature-are-tls-intercepted-on-this-machine` 所记的劫持态**已解除**。
判词是「该记忆已过期、须逐会话现测」，**不是**「旧观测为假」；**禁用 `-k` 的纪律不变**。

**两台仪器**（CLAUDE.md §2 强制）：

| 仪器 | 产物 |
|---|---|
| 文本层 | [`_scratch_ir/_bose0711_utf8.txt`](../_scratch_ir/_bose0711_utf8.txt)（`pdftotext -enc UTF-8 -layout`） |
| 页图 | `library/pdf/eta_hunt_2026-09-25/bose_p7-07.png`（`pdftoppm -r 200`，PDF 第 7 页） |

⚠ **两台仪器在本例确实分歧**：文本层把 `ω_ln` / `ω_m` / `ω̄` / `ω̃` 四行的**标签**相对数值**下移一行**。
η 的两行不受影响（逐位一致），但这正是「不许只信一种仪器」的又一次兑现——
若承重量落在那四行上，单靠文本层会读错。

## 2. TABLE III 逐字（两台仪器一致的部分）

图注逐字：

> `TABLE III: Hopfield parameters from the linear response calculation η and the rigid`
> `muffin-tin (atomic sphere) approximation η (RMT/RAS), the Fermi surface average of the`
> `square of the electron-phonon matrix element ⟨I²⟩ and ⟨I²⟩(RMT/RAS) from the linear`
> `response and the RMT results, respectively; ...`

| 行 | 单位 | a=6.255 | 6.15 | 6.05 | 5.85 | 5.75 |
|---|---|---|---|---|---|---|
| `a` | a_B | 6.255 | 6.15 | 6.05 | 5.85 | 5.75 |
| `η`（FP-LMTO 线性响应） | Ry/bohr² | 0.0381 | 0.0433 | 0.0519 | 0.0666 | 0.0771 |
| `⟨I²⟩` | (Ry/bohr)² | 0.0026 | 0.0031 | 0.0040 | 0.0056 | 0.0065 |
| **`η (RMT/RAS)`** | Ry/bohr² | **0.0439** | **0.0509** | **0.0591** | **0.0826** | **0.1077** |
| `⟨I²⟩(RMT/RAS)` | (Ry/bohr)² | 0.0030 | 0.0037 | 0.0045 | 0.0065 | 0.0081 |

（a 递减 = 压力递增。）

## 3. 验收门逐条

| 验收门项 | 判 | 依据 |
|---|---|---|
| 同一篇论文内 | PASS | TABLE III 一张表内并列两列 |
| 同一种过渡金属 | PASS | hcp Sc（最轻 3d 过渡金属），五个晶格常数 |
| 印的是 **η** 本身而非 λ | PASS | 行标 `η` 与 `η (RMT/RAS)` |
| 带单位 | PASS | `Ry/bohr²`，两行同单位 |
| 一端是 RMTA | PASS | `the RMT scheme ... known as the rigid atomic sphere (RAS) method` |
| 另一端非 RMTA | PASS | `the FP-LMTO linear response method` |

⟹ **六项全过。这是本课题至今第一篇通过验收门的文献。**

## 4. 逐列比值

| a [a_B] | η_LR | η_RMT/RAS | **比值** | ⟨I²⟩ 比值 | N(0)_LR / N(0)_RMT |
|---|---|---|---|---|---|
| 6.255 | 0.0381 | 0.0439 | 1.1522 | 1.1538 | 1.0014 |
| 6.150 | 0.0433 | 0.0509 | 1.1755 | 1.1935 | 1.0153 |
| 6.050 | 0.0519 | 0.0591 | **1.1387** | 1.1250 | 0.9879 |
| 5.850 | 0.0666 | 0.0826 | 1.2402 | 1.1607 | 0.9359 |
| 5.750 | 0.0771 | 0.1077 | **1.3969** | 1.2462 | 0.8921 |

**逐列全部 > 1 ⟹ RMTA 一致地高估 η。** 范围 `[1.1387, 1.3969]`。

### 4.1 一处对作者散文的精确化（我方新发现）

作者自陈：`with the level of the overestimation increasing with pressure`。
**逐列并不成立**：`1.1522 → 1.1755 → 1.1387 → 1.2402 → 1.3969`，在 a=6.05 处有真实回落。

舍入检验（印出值 3 位有效数字的最坏区间）：
`a=6.15 → [1.17301, 1.17803]`，`a=6.05 → [1.13667, 1.14079]`，**两区间不相交**
⟹ 回落**不是舍入伪影**。

⟹ 作者那句是**粗趋势**（端点 1.1522 → 1.3969 确实上升），**不得转述成「逐列单调」**。
此事不影响本轮承重量（消费点只用 `max`）。

### 4.2 一处必须说清的混杂（以及它为什么无害）

`η = N(0)⟨I²⟩`，但上表末列 ≠ 1 ⟹ **两端的 N(0) 也不同**
⟹ 这不是纯粹的「RMT vs 线性响应」，而是 **(RMT + LMTO-ASA) vs (线性响应 + FP-LMTO)** 的**复合**差。

**但对本预算这恰好是想要的量。** 消费点问的是：
「若把 Pap1977 的 RMTA η 换成一次现代非 RMTA 计算，数字会动多少？」——基组差**本来就该算在内**。

> 对照 [`savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md`](savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md) 的 `λ_Savr/λ_Pap`：
> 那里的 ⟨ω²⟩ 混杂是**有害**的（把误差变成**双边**，且 ⟨ω²⟩ 不可观测、无法归因）；
> 这里的基组混杂是**无害**的（复合量就是消费点要的量，单边、可直接代入）。
> **「有混杂」不是判词，「混杂朝哪边走、消费点要不要它」才是。**

## 5. ★ 方向审计：RMTA 的偏差方向**不统一**

`ceiling = 137.0959·√(η/M)`，`r = Tc_exp/ceiling ∝ 1/√η`。
RMTA **高估** η 达 κ ⟹ 真 η 小 ⟹ 真天花板低 ⟹ 真 `r = r_报 × √κ`（**上升，危险方向**）。

| 材料 | 类 | 方法偏差 | 方向 | 对 r |
|---|---|---|---|---|
| Al | sp | `η_RMTA/η_true = 0.318` | **低估 3.14×** | 下降（安全） |
| **Sc** | **早 3d** | **比值 ∈ [1.139, 1.397]** | **高估 ≤1.397×** | **上升（危险）** |
| Fe | 晚 3d | 作者自陈 `consistently UNDERestimate` | 低估 | 下降（安全） |

> `This is in contrast with the results for the late transition metal hcp Fe, where RAS was`
> `shown to consistently underestimate the Hopfield parameter, with the level of the`
> `underestimation increasing with pressure.` [ref 12 = Bose–Dolgov–Kortus–Jepsen–Andersen, PRB **67**, 214518 (2003)]

⟹ **不存在统一的 RMTA 偏差方向。**
「低估所以 PASS 保守」只能**逐材料**说，**不可推广**——这一条直接约束了
记忆 `rmta-underestimates-eta-for-sp-metals-by-3x` 后果 1 的适用域。
本条的价值不是给方向，而是给**危险方向上第一个实测量级**。

## 6. 代入头条

取 Sc 上实测的最大高估 `κ = 1.3969`（a=5.75 那列）：

> **`r(Nb₃Ge) ≤ 0.626 × √1.3969 = 0.7399 < 1`**，余量 **1.352×**。
> 翻判决需 `κ = (1/0.626)² = 2.5518×`，实测包封比它小 **1.827 倍** ⟹ **判决不翻**。

**外推限制（如实写）**：Sc 不是 Nb/V，更不是 A15 化合物。
本条把「过渡金属侧完全无数」升级为「有一个同篇 η 级实测点」，
**不**等于「已测得 Nb 的 RMTA 偏差」。仍属**指示性包封**。

### 6.1 与既有两个包封的关系（三者不可混用，来源不同）

| 来源 | κ | r(Nb₃Ge) 上界 | 性质 |
|---|---|---|---|
| 文献 RMTA λ(Nb) 区间 1.12–1.86 vs LR 1.26 | 1.4762 | 0.7606 | **实现扩散**，非方法偏差；且带 ⟨ω²⟩ 混杂 |
| **本文 Sc 的 η 级实测** | **1.3969** | **0.7399** | **方法偏差**，η 级，同篇 |
| 翻判决阈值 | 2.5518 | 1.0000 | — |

⟹ 头条现有**两条独立的、都远低于阈值的**包封。承重表述取更有方法学依据的一条：
**`r(Nb₃Ge) ≤ 0.7399`，理由是「即使按唯一实测到的过渡金属 RMTA 高估量级，仍 PASS」。**

> ⚠ **2026-09-25 同日订正**：「唯一实测到的过渡金属 RMTA 高估量级」已不再唯一——超快实测（Brorson PRL 64,2172）给出
> **V 上 κ=2.020**（1σ 2.176）、Klein 同源 **V₃Ga κ≤2.062**（1σ 2.461），均比本节 Sc 的 1.3969 大。承重包封改为
> Nb 锚 `r≤0.648` / V 锚 `r≤0.890–0.923` / 同源 A15 锚 1σ `r≤0.982`（2σ 翻）。本节 Sc 数值本身不动。
> 见 [`brorson1990_mansart2010_ultrafast_rows_2026-09-25.md`](brorson1990_mansart2010_ultrafast_rows_2026-09-25.md) §5。

## 7. 顺带：两份被 peer 判为「误抓取」的 PDF 实为本课题最核心的两篇

**旧陈述错在哪**：peer `eta-rmta-bias` 报 `library/pdf/eta_hunt_2026-09-25/PhysRevLett.28.801.pdf`
是「液晶声模」、`PhysRev.167.331.pdf` 是「TiO₂:V⁴⁺ 自旋回波」，判为**误抓取、与课题无关**。
成因：只读了 PDF **第 1 页开头几行**。APS 按**期刊页**取件，页首是**前一篇文章在同一页上的续尾**，
本篇标题在同页稍下方才开始。**页首 ≠ 本篇。**

**新陈述凭什么对**（全文 grep，非首行）：

- `PhysRevLett.28.801.pdf` 第 75–77 行逐字：
  `Electron-Phonon Interactions, d Resonances,` / `and Superconductivity in Transition Metals` /
  `G. D. Gaspari*f and B. L. Gyorffy`
  ⟹ = **Gaspari & Gyorffy, PRL 28, 801 (1972)**，即 **RMTA/η 的奠基论文**，
  Tier-0 语料里每一个 η 的方法出处。第 110 行并印出 `λ = n(e_F)⟨I²⟩/M⟨ω²⟩`。
- `PhysRev.167.331.pdf` 第 25 行逐字：`Transition Te1nperature of Strong-Coupled Superconductors`（OCR 的 m→1n）
  ⟹ = **McMillan, Phys. Rev. 167, 331 (1968)**，即 **McMillan 方程**本身。

⟹ 两份都**不是**误抓取，**不得删除**。教训同
[[audit-the-artifact-not-the-stale-build]]：判一份产物的身份要扫**全文**，首行不代表整份。

## 8. 引用链交叉印证

本文参考文献独立确认了仓库的一处书目钉死与两处方法归属：

| Bose ref | 内容 |
|---|---|
| [7] | S.Y. Savrasov, Phys. Rev. B**54, 16470** (1996) |
| [8] | S.Y. Savrasov and D.Y. Savrasov, Phys. Rev. B**54, 16487** (1996) |
| [12] | S.K. Bose, O.V. Dolgov, J. Kortus, O. Jepsen, O.K. Andersen, PRB **67**, 214518 (2003)（hcp Fe，反方向） |
| [31] | G.D. Gaspari and B.L. Gyorffy, PRL **28**, 801 (1972) |
| [32] | D. Glötzel, D. Rainer, H.R. Schober, Z. Phys. B **35**, 317 (1979)（RAS 方法） |

⟹ **第三方独立确认 16470 与 16487 是两篇不同论文**（我此前的钉死成立）。

## 9. 严格性与残余

- **全部数字 numerical-evidence**（两端都是他人发表的 DFT 计算值，非可认证量）。
- 单位常数 `137.0959 K` 从 `C_∞` + CODATA 独立重建、逐位一致。
- **残余未关闭**：
  1. Sc → Nb/V → A15 的外推无证据，只能当指示性包封；
  2. 未取一手的 hcp Fe 论文（PRB 67,214518），其「低估」目前只是 Bose 的转述；
  3. 仍无任何 **A15 化合物本身**的非 RMTA EPC 计算。
