# P3 · 氢化物之外（non-hydride）的声子介导超导图景 + α²F 数据源 + Semenok 系已发表非氢化物 λ_*

调研 agent：`p3-nonhydride`（文献核验专员）　日期：2026-09-09。

**证据等级**：**A** = 一手原文逐字/公式/表格/官方数据文件清单；**B** = 权威二手逐字转录、或对 A 级数字做的确定性算术；**C** = 普通二手/图注/推算/未核实推断。

**本次会话的工具限制（必读，影响 verified 判定）**：WebSearch 配额在本任务开始前已被本 session 用尽（200/200），全部检索改用 WebFetch 直取 arXiv / figshare API / Semantic Scholar API / 官方文档。**本 agent 无 Bash/下载能力**，故所有"数据文件"只能核到**页面/API 级别的文件清单与 license**，**不能核到文件内部是否真含目标材料的 α²F 数组**——下面表格的 `verified_exists` 严格按此口径填。Nature 系页面 303 重定向到 idp.nature.com（登录墙），PRL/PRB 摘要在 S2 API 里被出版社 elide，故若干条目只能核到"标题+卷页+作者"（标题本身含关键数字时仍算 A）。

---

## §0 一句话结论

**Semenok 自己的 SI 里就有非氢化物**：Table S1 给了 MgB₂/Nb/Nb₃Sn/非晶 Ga,Pb,Bi,PbBi 的 (Tc, λ, ξ, ω_log, ω_max)，而 Eq. 6 是 **λ_* = λ/ξ**（用 YH₉ 的 2.66/0.36 = 7.39 精确复现他正文亲口给的 λ_*=7.39，反演口径确证）——**因此可以直接读出已发表的非氢化物 λ_*：MgB₂ ≈ 5.44、Nb ≈ 4.1–5.25、Nb₃Sn ≈ 5.9–6.7、非晶 Pb ≈ 6.6、非晶 Ga ≈ 14.1**，这是 P3 §4 展开的**现成、材料无关性的交叉验证靶子**，且 **λ_* 全部 ≥ Debye 4.72（Nb 下端除外）、无一"近 Einstein 3.69"，与我们"λ_* 是谱宽 σ 的增函数、真实谱远离 δ 峰"的图像一致**。最值得做数值检验的非氢化物是 **MgB₂**（常压、λ 与 ω_log 与 Kong 2001 一手计算一致、α²F 有多条可得路径），但**逐点 α²F(ω) 数值文件本次未到手**（EPW 官方 tutorial 04 可生成 `mgb2.a2f`，JARVIS-SuperconDB 有 CC-BY-4.0 开放下载但内部是否含 MgB₂ 未核）。

---

## §1 非氢化物声子介导超导体的图景 + 一手 Tc

### 1.1 常压经典高 Tc 非氢化物

| 材料 | Tc | 压力 | 一手出处 | 等级 | 备注 |
|---|---|---|---|---|---|
| MgB₂ | **39 K** | 常压 | Nagamatsu, Nakagawa, Muranaka, Zenitani, Akimitsu, *Nature* **410**, 63–64 (2001)，标题逐字："Superconductivity at 39 K in magnesium diboride" | **A**（标题逐字，经 Semantic Scholar API `DOI:10.1038/35065039` 取回卷页作者；摘要被出版社 elide、未取得） | two-gap（σ/π）；**注意 Semenok Table S1 该行写 Tc=40 K**，与一手 39 K 差 1 K，引用时以 39 K 为准 |
| Cs₃C₆₀（A15，掺杂富勒烯） | **38 K** | 加压（"hydrostatic pressure"，**具体压力值本次未核到一手**；常见二手值 0.7 GPa/7 kbar 为 **C**） | Ganin, Takabayashi, Khimyak, Margadonna, Tamai, Rosseinsky, Prassides, *Nature Materials* **7**(5), 367–371 (2008)，标题逐字："Bulk superconductivity at 38 K in a molecular system" | **A**（标题+卷页作者，S2 API `DOI:10.1038/nmat2179`）；摘要仅取到片段"C(60)-based solids are archetypal molecular superconductors with transition temperatures (Tc) as high as 33 K"（A），且明确说加压在绝缘/超导之间调控而晶体结构不变 | 是否属"纯声子介导"本身有争议（Mott 邻近、Jahn–Teller）；对 P3 而言它**不是干净的声子介导检验对象** |
| Ba₀.₆K₀.₄BiO₃（铋酸盐） | **~30 K**（"near 30 K"） | 常压 | Cava, Batlogg, Krajewski, Farrow, Rupp, White, Short, Peck, Kometani, *Nature* **332**, 814–816 (1988)，标题逐字："Superconductivity near 30 K without copper: the Ba0.6K0.4BiO3 perovskite" | **A**（标题+卷页作者，S2 API `DOI:10.1038/332814a0`）；摘要 elide、未取得 | Semenok Table S2 另列 BaKBiO₃ 20 K / BaPbBiO₃ 13 K，且**把它归在 "Unconventional" 一栏**（A，见 §3）——即 Semenok 本人不把铋酸盐当声子介导 |

### 1.2 ★ 压力调控的**单质**超导体（P3 主题最贴的一类：有压力、有声子、无氢）

| 元素 | Tc | 压力 | 一手出处 | 等级 |
|---|---|---|---|---|
| **Sc（截至 2026 的单质纪录）** | **36 K** | **~260 GPa** | Ying, Liu, Lu, Wen, Gui, Zhang, Wang, Sun, Chen, **PRL 130, 256002 (2023)**，arXiv:2302.14378，标题逐字："Record-high-Tc elemental superconductivity in scandium"；摘要原句片段（A）："elemental scandium (Sc) can be increased to 36 K, which is a record-high Tc for superconducting elements" | **A**（arXiv 摘要原句 + 期刊卷页）。归因于 **Sc-V 相**，作者用第一性原理归结为 d 电子与中频声子的强耦合 |
| Sc（独立第二组，数值略低） | Tc^onset ~32 K（Tc^zero ~31 K） | ~283 GPa | He, Zhang, Li, ... Tse, Jin, **Chin. Phys. Lett. 40, 107403 (2023)**，arXiv:2303.01062，摘要原句片段（A）："Tconset of Sc increases from ~3 K at around 43 GPa to ~32 K at about 283 GPa (Tczero ~ 31 K)"；并称"no sign of saturation up to the maximum pressure achieved" | **A** |
| **Ti（Sc 之前的单质纪录）** | Tc^onset **>26 K**，Tc^zero 21 K | 峰值在 108–240 GPa 宽区间内，测到 310 GPa | Zhang, He, Liu, ..., Chen, Jin, **Nature Communications 13, 5411 (2022)**，arXiv:2112.12396，标题"Record High Tc Element Superconductivity Achieved In Titanium"；摘要原句片段（A）："maximum Tc^onset above 26 K and zero resistance Tc^zero of 21 K are record high values hitherto achieved among element superconductors" | **A** |
| **Ca** | **29 K** | **216 GPa** | Sakata, Nakamoto, Shimizu, Matsuoka, Ohishi, **PRB 83, 220512 (2011)**，标题逐字："Superconducting state of Ca-VII below a critical temperature of 29 K at a pressure of 216 GPa"（S2 API `DOI:10.1103/PhysRevB.83.220512` 取回标题+卷页+作者；摘要 elide、未取得） | **A**（标题即含全部关键数字） |
| S、Li、Y、La、黑磷/P、B 等 | — | — | **本次未核到一手**（WebSearch 配额已尽；arXiv 全文检索对"elemental sulfur"只返回 Struzhkin 等 cond-mat/0201520 的实验方法学文章，其中提到"observation of a Tc in sulfur to above 230 GPa"——**Tc 数值本身未取得**） | **未取得**，非"确认不存在" |

**单质纪录裁定（截至 2026-09-09，本次可核范围内）**：**Sc，36 K @ ~260 GPa，PRL 130, 256002 (2023)**。用 arXiv 全文检索 "record high Tc element superconductor pressure" 未见 2023 之后有超过 36 K 的单质报道；**但检索面有限（WebSearch 不可用），只能报"未见"，不能报"确认不存在"**。Talantsev, *Symmetry* **15**, 1632 (2023)（arXiv:2301.05148）作为独立二手也把 Sc/Ti 的 "up to 36 K" 认作单质纪录（**B**）。

**★ 对 P3 的意义**：Sc/Ti/Ca 这三条是**"高压 + 声子 + 无氢 + 有 Tc(P) 曲线"**的现成对象，正是 §4 λ_*(P) 展开想要的材料无关检验场。缺的只是它们的 α²F(ω;P) 数值——见 §2 的诚实结论（未到手）。

---

## §2 机器可读 / 可数字化的 α²F(ω) 数据源（本任务最高优先级）

```csv
system,pressure_GPa,source_type,url,format,license,verified_exists,notes
MgB2,0,tutorial+code(regenerable),https://docs.epw-code.org/doc/tutorials/tutorial_04/index.html,"mgb2.a2f (isotropic alpha2F(omega) + cumulative lambda(omega), 多个 smearing 列)",GPL(EPW/QE),yes(页面级),"★最实在的一条：EPW 官方 tutorial 04 Exercise 2 就是 MgB2 各向异性 Migdal-Eliashberg，输出文件明确含 mgb2.a2f；输入文件+赝势打包在 Google Drive 文件夹 1LZiO1c1LX4UVhvl7wWpoZzxwfiYiSw8l（gdown），赝势来自 pseudo-dojo.org。代价=要跑 EPW，不是即取即用的数值文件"
Pb,0,tutorial+code(regenerable),https://docs.epw-code.org/doc/tutorials/tutorial_04/index.html,"pb.a2f (col2=alpha2F, col3=cumulative lambda), pb.a2f_proj(按原子分解)",GPL(EPW/QE),yes(页面级),"同上 Exercise 1（isotropic, a2f_iso=.true.）。Pb 也在 Semenok Table S1 里（非晶 Pb 行），可做第二个非氢化物对照"
general(1058 materials incl. 3D),ambient,dataset,https://api.figshare.com/v2/articles/21370572,"jarvis_epc_data_figshare_1058.json.zip (8.83 MB), jarvis_epc_data_2d.json.zip (0.64 MB)",CC BY 4.0,partial,"★ figshare API 直取成功（上次 a2f-data-hunt 报的 403 是网页端，API 通）：文件名/大小/开放 download_url（ndownloader.figshare.com/files/38307921 与 /38950433）+ license 全部核实=A。与 NIST 页面 supercon_3d=1058 条一致（A）。未核：文件内部是否逐材料存 alpha2F 数组、MgB2 是否在内——本 agent 无下载能力。对应论文 Choudhary & Garrity, npj Comput. Mater.(2022), arXiv:2205.00060，摘要明确说把 Eliashberg function 当中间量预测（A），故存 alpha2F 的可能性高但未证"
general(818 materials),ambient,dataset+code,https://github.com/henniggroup/BETE-NET,"database.json + structures/(CIF) + indices/ + test_preds/",Apache-2.0,partial,"论文 arXiv:2401.16611 (Gibson et al., npj Comput. Mater. 2024) 摘要逐字：'compute alpha^2F(omega) for 818 dynamically stable materials'(A)，故 818 条 alpha2F 确实被算过；但 repo landing page 只说 database.json '含训练用数据库'，未给 schema，MgB2 是否在内未核"
MgB2,0,paper(figure/first-principles),https://arxiv.org/abs/cond-mat/0102499,PDF 图（无数值文件）,N/A,no(数值文件),"Kong, Dolgov, Jepsen, Andersen, PRB 64, 020501(R) (2001)，'Electron-phonon interaction in the normal and superconducting states of MgB2'（A：标题/作者/卷页/摘要经 arXiv abs 核实）。摘要给 mu*=0.10 复现 Tc=40 K、2Delta0/kBTc=3.9、sigma 带 + 590 cm^-1 bond-stretching 光学声子主导（A）。**只能数字化图**"
MgB2,0,tunneling-inversion,https://arxiv.org/abs/cond-mat/0305563,PDF 图,N/A,no,"Daghero, Gonnelli, Ummarino, Dolgov, Kortus, Golubov, Shulga, Physica C 408-410, 353 (2004)：从隧道谱反演 alpha2F(omega)（非第一性原理）。作为独立谱形对照有价值"
MgB2,0,anisotropic ab initio,https://arxiv.org/abs/cond-mat/0111182,PDF 图,N/A,no,"Choi, Roundy, Sun, Cohen, Louie, PRB 66, 020513 (2002)，各向异性 Eliashberg，摘要给 lambda=0.61、Tc=39 K（A）。注意其 lambda=0.61 与 Kong/Semenok 的 0.87 不同口径，选谱时须一致"
compressed S,~230,none_found,,,,no,"未取得。arXiv 检索只返回 Struzhkin et al. cond-mat/0201520（实验方法学，提到 Tc 测到 230 GPa 以上），无 alpha2F 数值源"
compressed Ca,216,none_found,,,,no,"未取得。arXiv 检索该关键词组只返回 Na-intercalated graphite 一条，无单质 Ca 的 alpha2F 源。**未核到 ≠ 确认不存在**"
compressed Sc,260,not_searched,,,,no,"WebSearch 配额已尽，未对 Sc 的第一性原理 alpha2F 做独立检索；PRL 130,256002 正文含第一性原理计算（摘要亲述），其 SI 可能有谱——**下一步优先项**"
```

**§2 结论**：**没有一条"即下即用的 MgB₂ α²F(ω) 数值文件"被本次核实到手。** 三条可行路径按代价排序：
1. **EPW tutorial 04 跑一遍 → `mgb2.a2f`**（最可靠、格式明确、含 α²F 与累积 λ(ω)，但要跑计算）；
2. **JARVIS-SuperconDB 下载 8.8 MB zip 核 schema 与 MgB₂**（CC BY 4.0、直链已核，只差一次下载）；
3. **数字化 Kong 2001 PRB 64,020501(R) 的 α²F 图**（谱形与 Semenok 用的 λ=0.87/ω_log=680 K 最可能同源，见 §3 注）。

---

## §3 ★ Semenok/Yuzbashyan 系是否已给非氢化物的 λ_*：**有，而且可以精确读出**

### 3.1 关键钥匙：λ_* = λ/ξ（口径已用作者自己的数字验证）

arXiv:2407.12922 HTML 全文（**A**）：
- Eq. 3 定义 ξ（**主定义是被极大化的积分，不是 λ/λ_***）：`ξ ≡ max_T {∫₀^∞ g(ω/2πT) · 2α²F(ω)/ω dω} < 1`
- Eq. S12：`λ_* = 1 / max_T {∫₀^∞ g(ω/2πT) P(ω) dω}`，`P(ω) = (1/λ)·2α²F(ω)/ω`
- Eq. 6：`λ_* = λ/ξ,  λ_c = λ ξ_c/ξ = λ_* ξ_c`
- 正文唯一显式给出的单材料 λ_*：**"we obtain λ_*=7.39 and (assuming ξ_c=0.5) λ_c=3.69 for YH9 at 205 GPa"**（逐字，A）
- 通用端点：Einstein `λ_* = 1/g_max = 3.6915`；Debye（α²F∝ω²）`λ_* = 4.72`（均逐字，A）

**验证**：Table S1 的 YH₉(205 GPa) 行是 λ=2.66、ξ=0.36 → **2.66/0.36 = 7.3889 ≈ 7.39** ✓。**故 λ_*=λ/ξ 的反演口径被作者自己的数字确证**，可放心套用到非氢化物行。

### 3.2 Semenok Table S1 的**非氢化物行**（A 级逐行）+ 反演出的 λ_*（B 级算术）

| 化合物 | Tc (K) | λ | ξ | ω_log (K) | ω_max (K) | **λ_* = λ/ξ**（本报告反演） | λ_*/3.6915 |
|---|---|---|---|---|---|---|---|
| **MgB₂** | 40 | **0.87** | **0.16** | 680 | 1624 | **5.44** | **1.473** |
| Nb | 9.2 | 0.82–1.05 | 0.2 | 163 | 325 | **4.10–5.25** | 1.11–1.42 |
| Nb₃Sn | 17.9 | 1.6–1.8 | 0.27 | 142 | 422 | **5.93–6.67** | 1.61–1.81 |
| Hg | 4.2 | 1.0–1.6 | 0.14 | 86 | 165 | **7.14–11.43** | 1.93–3.10 |
| Pb（非晶） | 7.2 | 1.91 | 0.29 | 35 | 155 | **6.59** | 1.78 |
| Ga（非晶） | 8.6 | 2.25 | 0.16 | 62 | 302 | **14.06** | 3.81 |
| Bi（非晶） | 6.1 | 1.84–2.46 | 0.19 | 42 | 163 | **9.68–12.95** | 2.62–3.51 |
| PbBi（非晶） | 7.0 | 3.0 | 0.26 | 33.3 | 150 | **11.54** | 3.13 |
| PbBi₃（非晶） | 6.8 | 2.78 | 0.23 | 33.4 | 154 | **12.09** | 3.28 |

（对照：同表氢化物 H₃S(157) λ=1.84/ξ=0.32→λ_*=5.75；LaH₁₀(214) 2.06/0.43→**4.79**；YH₆(165) 1.71/0.39→4.38；YH₉(205) 2.66/0.36→**7.39**✓；(La,Y)H₁₀(180) 3.87/0.47→8.23；CaH₆(172, DFT) 2.69/0.58→4.64；Li₂MgH₁₆(250, DFT) 3.30/0.49→6.73；金属氢(I4₁/amd, 500 GPa) 2.85/0.45→6.33。全部为 B 级算术。）

**精度警告**：ξ 只印到 2 位小数，ξ=0.16 的 ±0.005 舍入 → λ_* 有 **±3.1%** 不确定；ξ=0.14 时达 ±3.6%。故 MgB₂ 应报 **λ_* = 5.44 ± 0.17**（仅舍入误差，不含 λ 本身的口径不确定）。

### 3.3 ★ 与 P3 §4 展开的直接对照（这是本次最有价值的产出）

用 P3 §4 的 `λ_*/λ₀ = 1 + ½κσ² + c₃γ₃σ³`（λ₀=3.6915、½κ=2.802、c₃=0.8959）反解 σ（暂设 γ₃=0）：

| 材料 | λ_*（Semenok 反演） | λ_*/λ₀ | 反解 σ（γ₃=0） | 落点判读 |
|---|---|---|---|---|
| Nb（下端 λ=0.82） | 4.10 | 1.111 | 0.20 | **展开有效区**（对称谱 σ≲0.2 精度 ≤0.4%）——**最干净的检验点** |
| **MgB₂** | **5.44** | **1.473** | **0.41** | 已超 Debye（σ_Debye=0.354）→ **必须用 Eq. S12 直接积分、不能套展开** |
| Nb₃Sn | 5.93–6.67 | 1.61–1.81 | 0.47–0.54 | 远出展开区 |
| 非晶 Ga/Bi/PbBi | 9.7–14.1 | 2.6–3.8 | ≥0.75 | 极宽谱，展开完全失效（非晶，预期如此） |

**三条硬读数**：
1. **非氢化物的 λ_* 并不小、也不"近 Einstein"**：除 Nb 下端（4.10）外，全部 **≥ Debye 4.72**，非晶体系高达 14。这与 P3 §4 "λ_* 是谱宽 σ 的单调增函数、真实谱远离 δ 峰" 完全一致，且**把这个结论从氢化物推广到了非氢化物**——正是导师问的"Hx 之外还有没有"的直接答案：**有，而且 Semenok 自己算过。**
2. **MgB₂ λ_*≈5.44 比 LaH₁₀（4.79）还高**，即 MgB₂ 的 α²F 比 LaH₁₀ 的**更宽**（σ≈0.41 vs LaH₁₀ 的 Debye 尺度）。物理上说得通：MgB₂ 有低频 Mg/声学段 + 尖锐的 E₂g bond-stretching 峰（Kong 2001 亲述 590 cm⁻¹ 主导），是典型的"宽底+强峰"双结构、偏度大。**这使 MgB₂ 成为检验 c₃γ₃σ³ 偏度项的好靶子，但同时也意味着二阶展开在此处不够用。**
3. **★ 最干净的定量检验是 Nb，不是 MgB₂**：Nb 的 λ_*≈4.1–5.25 对应 σ≈0.20–0.35，下端正落在展开的高精度区；且 Nb 的 α²F 是全领域最标准、最多独立来源的谱（EPW tutorial 里 Pb 同理）。**建议数值检验的顺序改为 Nb/Pb（展开有效区，检验 κ）→ MgB₂（宽谱，检验 Eq. S12 直接积分 + c₃ 偏度项）。**

### 3.4 另两篇的核查结果（均为**否**）

- **Yuzbashyan–Altshuler–Patra, arXiv:2409.19562**（"Instability of metals with respect to strong electron-phonon interaction"）：摘要页无任何单材料数字、无 λ_* 表、无非氢化物条目（**A**，摘要页级；全文未取）。
- **Kiessling–Altshuler–Yuzbashyan, arXiv:2409.00532**：纯数学论文（紧自伴算子的本征值界），摘要明确给的是抽象的 `T_c(λ,P)` 上下界序列与大 λ 渐近 `T_c ~ C√⟨ω²⟩√λ`，并称其常数 C 比 sharp constant **大约 2.034 倍**（**A**，摘要级）——**无任何具体材料**（MgB₂/Nb/Pb 均未出现于摘要）。
  - **★ 顺带的一个有用交叉点**：这个"比 sharp constant 大 2.034 倍"的说法，与 P3 §7 讨论的"AD 拟合自带渐近 0.1878 vs 真 Eliashberg 渐近 0.1827（+2.8% 过冲）"**不是同一件事**（2.034 是数学上界相对 sharp 常数的倍数，不是 AD 拟合的定标差），**切勿混用**。全文核实该常数的确切定义是 P3 后续可选项。

**注（未证的溯源推断，C 级）**：Semenok 的 MgB₂ 行 λ=0.87 与 Kong–Dolgov–Jepsen–Andersen PRB 64, 020501(R) (2001) 的 λ 值一致（该文摘要用 μ*=0.10 复现 Tc=40 K，与 Table S1 的 Tc=40 K 也一致），**很可能 Semenok 的 MgB₂ ξ 就是从 Kong 的 α²F 算的**。若属实，则用 Kong 的 α²F 做 Eq. S12 直接积分应当**精确复现 ξ=0.16 / λ_*=5.44**——这将是 P3 §4 机器的一次闭环验证。**本次未核 Semenok 参考文献列表，此条仅为待验假设。**

---

## §4 2024–2026 非氢化物高压/常压新进展（简短，均 A 级摘要）

- **arXiv:2509.23137（2025）**：压缩下**钠插层石墨**超导 **22.3 K**。摘要同时给出前人对照："Tc limited to 11.5 K at ambient pressure and 15.1 K at 7.5 GPa in calcium-intercalated graphite"。→ 石墨插层化合物（GIC）是一条活跃的**无氢、有压力**路线，Tc 已进入 20 K 段。
- **arXiv:2509.20742（2025）**：**空穴掺杂 YB₂ 常压 Tc 22.83 K** 的第一性原理 + Migdal–Eliashberg 预测（未掺杂 0 GPa 时仅 2.14 K）。作者自陈这不是破纪录结果（未超 MgB₂）。→ 硼化物家族仍在被系统扫。
- **arXiv:2310.13945（2023，略早于窗口）**：**Sr₂B₅ 预测 Tc ~100 K @ 40 GPa**。→ 若成立，是"非氢化物 + 中等压力 + 三位数 Tc"的最激进预测，**对 P3 §4 极有价值**（可反过来问：这个预测的 α²F 对应 λ_* 是多少、λ 是否已越过稳定界）。**本次仅核到检索结果级摘要（B），未取原文。**
- **arXiv:2510.26605（2025）**：金刚石 NV 量子传感推到 240 GPa，并在**单质 Ti** 的超导转变上给出 180 GPa 的迈斯纳效应/俘获磁通证据。→ 单质高压超导的**表征手段**在升级（Meissner 而非仅电阻），意味着 Sc/Ti 一类数据的可信度在提高。
- **本次未见** 2023 之后有超过 Sc 36 K 的单质报道（检索面有限，见 §1.2 裁定）。

---

## §5 诚实边界（未取得 / 被墙 / 存疑）

**未取得（非"确认不存在"）**：
1. **MgB₂ 的可下载 α²F(ω) 数值文件**——三条路径都只核到"存在且可得"，无一核到手（§2）。
2. **压力下单质（S / Ca / Sc）的 α²F(ω;P)**——一条都没找到。Sc 的 PRL 130,256002 正文含第一性原理计算，其 SI 是**下一步最该开的一扇门**，本次因配额未开。
3. **S、Li、Y、La、黑磷/P、B 的高压 Tc 一手值**——全部未核。
4. **2409.19562 与 2409.00532 的全文**（只核到摘要页），故"无非氢化物 λ_*"是**摘要级判断**，不排除正文/附录里有。
5. **Semenok 2407.12922 的参考文献列表**——故 §3.4 末尾"MgB₂ 数据源自 Kong 2001"仅为待验假设（C）。
6. **Cs₃C₆₀ 38 K 对应的具体压力值**——一手未核到（常见二手 0.7 GPa/7 kbar 为 C 级）。

**被墙**：Nature / Nature Materials 全文页 303 重定向到 idp.nature.com 登录墙；PRL/PRB 摘要在 S2 API 中被出版社 elide（Nagamatsu 2001、Cava 1988、Sakata 2011 三条只核到标题+卷页+作者，好在这三条的**标题本身就含关键数字**）。figshare 网页端 403 但 **API 通**（这更正了 `a2F_machine_readable_sources_2026-09-08.md` 里 JARVIS 那行的 `partial/403` 判定）。

**存疑（建议下一步核）**：
1. **Semenok Table S1 的 MgB₂ ω_max = 1624 K ≈ 140 meV ≈ 1128 cm⁻¹**。MgB₂ 的声子谱通常认为截止在 ~100 meV（~800 cm⁻¹）附近，Kong 2001 强调的主导模在 590 cm⁻¹。**1624 K 看起来偏高**，可能是不同定义（如含 anharmonic 尾、或取某个计算的最高分支）。**凡用到 MgB₂ 的 Tc/ω_max 比值（§6 U7 那类分析）前必须先核这个数。**
2. **MgB₂ 的 λ 有两个流通口径**：0.87（Kong/Semenok）与 0.61（Choi 2002 各向异性）。λ_* 反演对 λ 线性敏感——用 0.61 会得 λ_*=3.81（几乎贴 Einstein 底），**结论会翻转**。做数值检验时 α²F 与 λ 必须同源。
3. **Semenok 把 BaKBiO₃/BaPbBiO₃ 归入 "Unconventional"**（Table S2，A 级）。故铋酸盐**不能**被我们当作"非氢化物声子介导"的支持证据来用，除非另立论证。
4. MgB₂ 的 Tc：一手 39 K vs Semenok 表 40 K（1 K 差，无害但引用须一致）。

---

## 附：对 `a2F_machine_readable_sources_2026-09-08.md` 的两处更正建议

1. JARVIS-SuperconDB 行：`verified_exists` 从 `partial(WebFetch 403)` 改为 `partial(API 已核：2 文件、8.83+0.64 MB、CC BY 4.0、直链可下；内部 schema 未核)`，并补 API 端点 `https://api.figshare.com/v2/articles/21370572`。**教训：figshare 网页端 403 时改用 `api.figshare.com/v2/articles/<id>` 通常可通**（此法上次也是这样拿下 H3S/D3S 的 31136635）。
2. 建议新增两行：EPW tutorial 04 的 `mgb2.a2f` 与 `pb.a2f`（GPL、页面级已核、需跑计算）。
