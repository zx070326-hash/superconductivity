# cuprate 超快 `λ⟨ω²⟩` 行一手化：Gadermaier 2010 与 Perfetti 2007——二手表的 YBCO 行 G1 失败（r 虚高 1.5×），Bi2212 行余量 4.7× 且失效条件可陈述

**日期**：2026-09-25
**上游**：[`ceiling_as_phonon_mediation_discriminator_2026-09-25.md`](ceiling_as_phonon_mediation_discriminator_2026-09-25.md) §3–§5（二手表、三道门、采购 P2/P3）·
[`brorson1990_mansart2010_ultrafast_rows_2026-09-25.md`](brorson1990_mansart2010_ultrafast_rows_2026-09-25.md)（同日 P1 Mansart）
**严格性**：全篇 `numerical-evidence`（他人发表数值 + float64 点算）。「Perfetti 的 `λΩ₀²` 即 Einstein 谱的 `λ⟨ω²⟩`」为定义代入（`α²F=(λΩ₀/2)δ(Ω−Ω₀)` ⟹ `2∫α²F·Ω dΩ = λΩ₀²`），analytic 但平凡。
**脚本**：[`_scratch_ir/_chk_main_cuprate_ultrafast_rows_gadermaier_perfetti_2026-09-25.py`](../_scratch_ir/_chk_main_cuprate_ultrafast_rows_gadermaier_perfetti_2026-09-25.py)（三节，含正对照与投毒，退出码 0）

判别器读法与方向沿用上游：「声子介导的各向同性 Eliashberg 解释不了该 Tc」需 `X := λ⟨ω²⟩` 的**上**界 ⟹ **一律条件式**，条件写成「所需低估倍数 `r²`」并与**该文自身**已点名的模型因子比较。

---

## 1. 来源（G2）与仪器

检索由 `general-purpose` agent 完成（harvest 路线，两篇 APS 均为真文本层 PDF）；**全部承重字符由主会话读页图复核**：

| 文献 | 承重内容 | 页图 |
|---|---|---|
| C. Gadermaier *et al.*, **PRL 105, 257001 (2010)** | 摘要 `λ⟨ω²⟩ = 800 ± 200 meV²`（La₁.₈₅Sr₀.₁₅CuO₄）、`400 ± 100 meV²`（YBCO）；Eq.(1)(2)；样品 Tc | [`gadermaier_prl105_p1_abstract-1.png`](../../../library/pdf/eta_hunt_2026-09-25/gadermaier_prl105_p1_abstract-1.png)、[`…_p2_eq1eq2-2.png`](../../../library/pdf/eta_hunt_2026-09-25/gadermaier_prl105_p2_eq1eq2-2.png)、[`…_p1_samples_tc-1.png`](../../../library/pdf/eta_hunt_2026-09-25/gadermaier_prl105_p1_samples_tc-1.png) |
| L. Perfetti *et al.*, **PRL 99, 197001 (2007)** | 摘要 `λΩ₀² = 360 ± 30 meV²`；速率方程 (1)–(3)；Table I；样品 Tc | [`perfetti_prl99_p1_abstract-1.png`](../../../library/pdf/eta_hunt_2026-09-25/perfetti_prl99_p1_abstract-1.png)、[`…_p3_rate_eqs-3.png`](../../../library/pdf/eta_hunt_2026-09-25/perfetti_prl99_p3_rate_eqs-3.png)、[`…_p3_tableI-3.png`](../../../library/pdf/eta_hunt_2026-09-25/perfetti_prl99_p3_tableI-3.png)、[`…_p1_sample_tc-1.png`](../../../library/pdf/eta_hunt_2026-09-25/perfetti_prl99_p1_sample_tc-1.png) |

**两处派活侧缺陷（如实记）**：
1. **我的派活 prompt 把 arXiv:1406.6207 错标为 Giannetti 综述**。agent 查实 1406.6207 是另一篇（ar5iv 标题 *Ultrafast spectroscopy of quasiparticle dynamics in cuprate superconductors*），Giannetti 综述是 **arXiv:1601.07204**（Adv. Phys. 65, 58 (2016)），并改走后者。
   两个一次源（Gadermaier、Perfetti）已一手取得，所以不影响本篇结论；但上游笔记里「综述 arXiv:1406.6207 转述」那行的**出处**仍未核实，本篇不替它补。
2. **agent 把 Gadermaier Eq.(2) 抄成 `(2/3)k_BT_l/(ħτ)`，丢了一个 π**；页图是 **`λ⟨ω²⟩ = (2π/3)·k_BT_l/(ħτ_{e-ph})`**。
   不修的话 NEM/TTM 换算比会错 π 倍（脚本投毒项）。

## 2. G1：二手表的 YBCO 行用错了 Tc（主要发现）

Gadermaier 样品逐字（页图）："Single crystals of **YBa₂Cu₃O₆.₅ (YBCO, T_c = 60 K)** and La₁.₈₅Sr₀.₁₅CuO₄ (LSCO, T_c = 38 K)"，全部实验在室温下做。
上游二手表把 `400` 配给 "YBCO, Tc 90" 与 "YBa₂Cu₃O₆.₉, Tc 92"——**两行都是错配**：`400 ± 100` 属于**欠掺杂 O₆.₅、Tc = 60 K** 的样品。

| | Tc | r | 所需低估 `r²` |
|---|---|---|---|
| 二手（Tc=90，脚本复现上游 2.122 作正对照） | 90 | 2.1222 | 4.504× |
| **一手同样品** | **60** | **1.4148** | **2.0016×** |

⟹ 上游的 G1 门（「不得用体系公认 Tc 补位」）**在第一个一手化的 cuprate 行上就兑现了**：二手 r 虚高正好 90/60 = 1.5×，所需低估从 4.5× 掉到 2.0×。
LSCO（Tc 38，X 800）r = 0.6336 < 1，无判别力，不进 PASS 语料（前提存疑）。

## 3. G3：Gadermaier 用的不是 Allen 式

页图逐字：Eq.(1)（TTM）`λ⟨ω²⟩ = (π/3)·k_BT_e/(ħτ_{e-ph})` = Allen PRL 59,1460 Eq.(14) 首项；
Eq.(2)（NEM，该文 Ref.[14] 的非平衡模型解析解；Ref.[14] 未取）`λ⟨ω²⟩ = (2π/3)·k_BT_l/(ħτ_{e-ph})`，**报出的 400/800 用的是 Eq.(2)**（作者依据：τ 与通量无关，与 TTM 的 `τ ∝ T_e` 矛盾）。
同一 τ 下 `X_TTM/X_NEM = T_e/(2T_l)`；取文中 "T_e is in the range 400–800 K"（**一般性陈述，非本实验实测 T_e**）与室温 `T_l`：换算比 ∈ **[0.667, 1.333]**。
两式是同一个 Eliashberg 矩在两种动力学假设下的读法 ⟹ G3 记为「**同一量、模型依赖因子 ≤1.33×**」。NEM 本身相对真值的偏差方向**未知**。

| YBa₂Cu₃O₆.₅ 情形 | X | r |
|---|---|---|
| 中心 | 400 | 1.4148 |
| +1σ | 500 | 1.2654 |
| +1σ × TTM 上端 | 666.7 | 1.0959 |
| +2σ × TTM 上端 | 800 | **1.0004** |

⟹ **条件式判别器，不入账**：1σ 误差棒与模型换算叠在一起仍 r>1，2σ 叠加后正好落在阈值上。性质与 Mansart 行相同（数值上也巧合：60/√400 = 24/√64 = 3）。

## 4. Perfetti Bi2212：余量最大的一行，失效条件可以具体陈述

页图逐字：样品 "optimally doped Bi₂Sr₂CaCu₂O₈₊δ (transition temperature of 91 K)" ⟹ **G1 满足**。
模型：Einstein 谱 `F(Ω)=δ(Ω−Ω₀)`，`λ = 2∫Ω⁻¹α²F dΩ`，三浴速率方程 (1)–(3)，电子方程右端用**完整 Bose 因子** `3λΩ₀³(n_e−n_p)/(ħπk_BT_e)`
（高温极限退化为 Allen 首项形，脚本验到 0.9968）⟹ **没有 Allen Eq.(14) 那种被丢掉的负修正项**。
Table I：`λΩ₀²` = 300–380（T_l=30 K）/ 340–380（300 K），区间端点对应 Ω₀ = 40 / 70 meV 的选择；摘要 `360 ± 30`。
**低估来源只有一个，且是作者明写的**："Electron-phonon scattering with the (1 − f) lattice modes that are more weakly coupled barely contributes to the temporal evolution of T_e and has been therefore **neglected**"（f = 0.13–0.25）
⟹ 拟合出的是**强耦合子集**的 `λ⟨ω²⟩` ⟹ 全谱 X 被低估，方向对判别器不利。

| X（meV²） | r | 所需低估 | ⟺ 被忽略的 (1−f) 模须携带子集的 |
|---|---|---|---|
| 360（中心） | 2.2618 | 5.116× | ≥ 4.12 倍 |
| 380（Table I 上端） | 2.2015 | 4.847× | ≥ 3.85 倍 |
| 390（摘要 +1σ） | 2.1731 | 4.722× | ≥ 3.72 倍 |

⟹ **判词（条件式）**：「若被作者忽略的弱耦合模的 `λ⟨ω²⟩` 合计不超过强耦合子集的 **3.7 倍**，则 Bi2212 的 91 K 不能由各向同性声子 Eliashberg 解释」。
这比「所需低估 4.7×」多一层信息：它指出了**该由什么测量来关掉**——全谱 `λ⟨ω²⟩`（例如另一种不分浴的实验，或第一性原理的 `η/M`），不是再测一次同一个 τ。
作者自己的物理读法（这些模 "barely contribute"）与失效条件相悖，但那是**拟合敏感度**陈述，不是界 ⟹ 不计入。

## 5. 对上游二手表的净效果

| 上游二手行 | 一手化后 |
|---|---|
| YBCO 400 / Tc 90 / r 2.122 | **作废**：一手是 YBa₂Cu₃O₆.₅ Tc 60 ⟹ r 1.4148，条件式（§2–3） |
| YBa₂Cu₃O₆.₉ 400±100 / Tc 92 | **作废**：400±100 不属于 O₆.₉；O₆.₉ 的值（综述转述 Gadermaier PRX 4, 011056 (2014) 的 450）**未一手核实**，该 PRX 的逐材料 Tc 对位 agent 也未完成 |
| Bi2212 380 / Tc 90 / r 2.177 | **一手化**：Tc 91（同样品）、X 300–390 ⟹ r 2.17–2.48，条件式、余量 ≥4.7×（§4） |
| LSCO 800±200 / 38 / 0.634 | 数值一手确认；无判别力 |

⟹ cuprate 族的判别余量不再是上游写的「4.5–4.7× 相当稳固」：**YBCO 的那份余量来自 Tc 错配**；真正稳的只剩 Bi2212 一行，且它依赖一个可陈述、可检验的条件。

> ⚠ **同日订正（§7）**：「真正稳的只剩 Bi2212 一行」过度推广——PRX 4 给出同样品 Tc 的 YBa₂Cu₃O₆.₉（r²=3.87×）与 Hg1201（3.69×）两条一手行。上游那行的错配只作废**那一行**。

## 6. 未做 / 未取得

- 综述转述的其余值（Dal Conte Science 335,1600 的 1400±700；Carbone PNAS 105,20161；Chia NJP 15,103027 的 160–640）：**全部二手、未核实**。
- ~~PRX 4 的 Hg1201 580 / YBa₂Cu₃O₆.₉ 450 / Bi2212 720 未核实~~ ⟹ **已一手化，见 §7**。

## 7. 同日续：Gadermaier PRX 4, 011056 (2014)，六个 cuprate 样品，G1 全满足

检索与裁图由 sonnet `general-purpose` agent 完成；**承重字符（Table I、B 常数、Eq.(1) 上下文）主会话读页图复核，数值由主会话独立脚本从零重算**
（[`_scratch_ir/_chk_main_gadermaier_prx4_rows_2026-09-25.py`](../_scratch_ir/_chk_main_gadermaier_prx4_rows_2026-09-25.py)，退出码 0）。
页图：[`gadermaier_prx4_tableI-p2-2.png`](../../../library/pdf/eta_hunt_2026-09-25/gadermaier_prx4_tableI-p2-2.png)、
[`gadermaier_prx4_Bconst-p4-4.png`](../../../library/pdf/eta_hunt_2026-09-25/gadermaier_prx4_Bconst-p4-4.png)、
[`gadermaier_prx4_eq1-p3-3.png`](../../../library/pdf/eta_hunt_2026-09-25/gadermaier_prx4_eq1-p3-3.png)。

**该文不印 `λ⟨ω²⟩`**，只印 Table I 的（样品 Tc, `τ_e-ph`）与 Eq.(1) `1/τ = B·λ⟨ω²⟩/T`，图注 **`B = 3h/2πk_B`**（页图上 h **无横线**）。
**h 还是 ħ 是承重的**（差 2π）。由该文**自身**两处数字判定为 ħ（即 NEM 式 `λ⟨ω²⟩ = (2π/3)k_BT/(ħτ)`，与同组 PRL105 Eq.(2) 同构）：
① p.4「λ ∼ 0.3–0.5 for cuprates」（取 `√⟨ω²⟩ = 40 meV`）——ħ 读法给 0.29–0.53，落在一位有效数字舍入区间；字面 h 读法给 0.05–0.08；
② 同组 PRL105 的 LSCO `τ = 45 fs` ⟹ ħ 读法 792（印 800±200），字面 h 读法 126。
（脚本初版用拍的 ±0.02 / 1% 容差，被 0.530 与 792 判 FAIL；改为**印刷精度的舍入区间**后通过——判据换成印刷精度，不是为过而放宽。）

T 取正文 300 K（Fig.1 图注写 295 K，给更小 X、更大 r ⟹ 300 K 是保守侧）；`τ` 两位有效数字 ⟹ r 取舍入上端 X：

| 样品（Table I 逐字） | 该文 Tc | τ (fs) | X (meV²) | r（舍入上端 X） | 所需低估 r² |
|---|---|---|---|---|---|
| YBa₂Cu₃O₆.₅ | 63 | 72 | 495.0 | 1.3308 | 1.771× |
| **YBa₂Cu₃O₆.₉** | **90** | 77 | 462.8 | **1.9665** | **3.867×** |
| **HgBa₂CuO₄.₁** | **98** | 62 | 574.8 | **1.9199** | **3.686×** |
| Bi₂Sr₂CaCu₂O₈.₁₄ | 80 | 49 | 727.3 | 1.3918 | 1.937× |
| La₁.₉Sr₀.₁CuO₄ | 30 | 42 | 848.5 | 0.4828 | — |
| La₁.₈₅Sr₀.₁₅CuO₄ | 38 | 45 | 792.0 | 0.6333 | — |

二手三值（Hg1201 580 / O₆.₉ 450 / Bi2212 720）与反解值差 +0.9% / −2.8% / −1.0% ⟹ **二手值与这张表的反解一致，数值证实**（二手源是否正是这样反解的，未核）；
错的只是上游表把 PRL105 的 O₆.₅ 值 400 配给了 O₆.₉（§2）。

**模型换算方向与 §3 相反**：该文明写低通量下 `T_e ∼ T_L` ⟹ `X_TTM/X_NEM = T_e/(2T_L) ≈ 1/2` ⟹ Allen/TTM 读法只会**抬高** r（O₆.₉ 给 2.79）。
§3 的 [0.667, 1.333] 来自 PRL105 的高通量一般陈述，不适用于本文 ⟹ 这里模型换算**不是**风险项。
剩下的条件只有一条：**NEM 值相对真 `λ⟨ω²⟩` 的偏差方向未知**，判词仍是条件式。

**旧陈述错在哪**（§5 末句与上游框「稳的只剩 Bi2212 一行」）：那句只基于 PRL105 一个 YBCO 样品（O₆.₅，Tc 60），
把「这一个二手配对是错配」推广成了「YBCO 族没有大余量」——**错配作废的是那一行，不是那个材料**。
**新陈述凭什么对**：PRX 4 Table I 对 O₆.₉ 与 Hg1201 印了**同样品** Tc（G1 满足）、`τ` 两位有效数字、B 的读法由该文内证钉死；
所需低估 3.87× / 3.69×（舍入上端），模型换算只加大余量。
⟹ **cuprate 判别器现有三条 r²>3.5 的一手行**：Perfetti Bi2212（4.7–5.1×，失效条件已陈述）、YBa₂Cu₃O₆.₉（3.87×）、Hg1201（3.69×）；全部**条件式**。

**两处交叉（`numerical-evidence`，跨样品、只作指示）**：
- **Bi2212**：Gadermaier（Tc 80，NEM 全速率）X = 727 vs Perfetti（Tc 91，强耦合子集）360，比 2.02。
  若把 727 移植到 Perfetti 样品，r = 1.591、所需低估 2.53×——即 Perfetti 的失效条件（被忽略模 ≥3.7× 子集）**未被触发**（本数据对应 ≈1.0×），但**余量从 4.7× 显著缩小**。两样品掺杂不同、模型不同（G1/G3 不同构），不入账。
- **YBa₂Cu₃O₆.₅**：PRX 4（72 fs，Tc 63）vs PRL105（100±20 fs，Tc 60）⟹ X 495 vs 356–400，`τ` 差 >1σ。同名样品、同组、四年间的测量不一致——**同一组的 `τ` 复现性本身约 30%**，这是 NEM 行误差的一个实测量级。

**本节严格性**：全部 `numerical-evidence`；「B 是 ħ」为**文内证据判定**（两处独立内证 + 同组前文同构），非作者勘误。
- 上游二手表 "arXiv:1406.6207" 出处：未核实（见 §1）。
