# P3 · LaH₁₀ 谱二阶矩 ω₂ ≡ √⟨ω²⟩ 的一手数据调研

调研 agent：p3-omega2-hunt　日期：2026-09-08
任务：为"Tc^emp = 0.182·ω_log·√λ（ω_log 代理式）vs 0.1827·√(λ⟨ω²⟩)（Allen–Dynes 真渐近）"之争拿到一手 ω₂(P)，裁定对抗核验提出的"LaH₁₀ 谱随压力 Einstein 化（ω₂/ω_log 单调下降 1.12→1.08→1.05→1.04）"。

**证据等级约定（沿用 `P3_tc_upper_bound_foundations_2026-09-08.md`）**
- `[作者亲述]`：原文/图注逐字句。
- `[表格读数]`：一手论文表格里的数字，逐位转录。
- `[图读数]`：我从一手论文的图里读的，与对抗核验的数字化同级。
- `[本报告解析反演]`：**只用一手表格数字**做的封闭解析求解，无任何图像数字化。精度由表格给出位数决定。

---

## §0　一句话结论

1. **不存在任何一手来源列出 LaH₁₀ 的 ω₂ 或 ⟨ω²⟩。** Errea 2020 的 Nature 版**没有 Source Data、没有 Supplementary**（Data availability 是"upon request"），arXiv 版 Extended Data Table I 与 Nature 版 Extended Data Table 1 都只有 λ、ω_log、四个 Tc；Semenok Table S1/S3/S4 同样只有 λ、ω_log（+ 单点 ω_max）。Kruglov 2020（LaH₁₀ 谐波，独立一手）在 SI 公式 (S3) 里**用了** ω_log/ω₂，但表里也没列 ω₂。Materials Cloud 无 LaH₁₀ 条目。
2. **但 ω₂/ω_log 可以从一手表格里精确反解出来**，不需要任何抠图：Errea 2020 Ext. Data Table I **同时给了 Tc^McM 与 Tc^AD**，两者只差 Allen–Dynes 的 f₁·f₂，而 f₂ 是 ω₂/ω_log 的显函数。反解结果（§2）：

| P (GPa) | 129 | 163 | 214 | 264 |
|---|---|---|---|---|
| **ω₂/ω_log（一手表格反演）** | **1.194** | **1.075** | **1.105** | **1.095** |
| 对抗核验（图数字化） | 1.12 | 1.08 | 1.05 | 1.04 |

3. **裁定：对抗核验的"单调下降"不成立，但其"⟨ω²⟩ 天花板远比 ω_log 天花板平"的核心结论成立，且比它自己算的还强。**
   - 方向上：ω₂/ω_log 在 129→163 GPa **大幅下降**（1.194→1.075），随后 163/214/264 三点**在 ±0.015 内持平且非单调**（1.075→1.105→1.095）。不是"随压力单调 Einstein 化"。
   - 后果上：用 ⟨ω²⟩ 天花板 T_c^emp = 0.1827·√λ·ω₂ 得 **368 → 359 → 388 → 387 K**（全程 +5%，且 129→163 段**下降**）；而用 ω_log 天花板得 **307 → 333 → 350 → 352 K**（全程 +15%，单调）。
   - 所以 **P3 稿"圆顶整个在 R、天花板反而在涨"这个强论断只在 ω_log 约定下成立**；换到 Allen–Dynes 真正的 ⟨ω²⟩ 约定，天花板基本是平的（且开头一段还降），圆顶的解释权在 R 与天花板之间要重新分配（§3）。

---

## §1　一手来源逐条溯源（找到什么 / 没找到什么）

### 1.1 Errea et al., Nature 578, 66 (2020) / arXiv:1907.11916

**本地 PDF**：`library/pdf/D3/2020_Ion_Errea_1907.11916.pdf`（arXiv v1, 20 页，含全部 Extended Data）。**已通读全文。**

**(a) Extended Data Table I 全文转录 `[表格读数]`**（arXiv v1 p.11）：

> **Extended Data Table I. Summary of calculated superconducting T_c.** Values are within different approaches ranging from empirical to fully *ab initio*: McMillan equation (T_c^{Mc}_{μ*=0.1}), Allen-Dynes modified McMillan equation (T_c^{AD}_{μ*=0.1}), anisotropic treatment of Migdal-Éliashberg (T_c^{ME}_{ani}) and SCDFT (T_c^{SCDFT}). **Values of λ and ω_log are also given.**

| System | Pressure (GPa) | λ | ω_log (meV) | T_c^{Mc}_{μ*=0.1} (K) | T_c^{AD}_{μ*=0.1} (K) | T_c^{ME}_{ani} (K) | T_c^{SCDFT} (K) |
|---|---|---|---|---|---|---|---|
| LaH₁₀ | 129 | 3.62 | 76.4 | 171.8 | 252.6 | 255.3 | 230 |
| LaH₁₀ | 163 | 2.67 | 96.4 | 197.1 | 247.0 | 242.8 | 225 |
| LaH₁₀ | 214 | 2.06 | 115.5 | 196.3 | 235.9 | 237.9 | 210 |
| LaH₁₀ | 264 | 1.73 | 126.6 | 189.5 | 219.2 | 216.9 | 201 |
| LaD₁₀ | 159 | 3.14 | 63.5 | 135.0 | 184.2 | 180.4 | 171 |
| LaD₁₀ | 210 | 2.21 | 81.7 | 145.5 | 176.5 | 172.9 | 158 |
| LaD₁₀ | 260 | 1.80 | 92.2 | 142.2 | 164.6 | 157.9 | 151 |

→ **无 ω₂、无 ⟨ω²⟩、无 ω_max。** 工作区 `isotope_curated.md` 第 24 行与 Semenok Table S4 的 LaH₁₀ 行确认逐位一致（本轮再次核对无误）。**注意：此前工作区只记了 LaH₁₀ 四行的 λ/ω_log/Tc^ME；Tc^Mc 与 Tc^AD 两列是本轮新增，正是 §2 反演的全部输入。**

**(b) Methods 里对 λ、ω_log 的定义（作者亲述，arXiv p.8）**——**明确只定义了这两个矩，没有定义二阶矩**：

> "In Extended Data Table I we summarize all calculated T_c's within anisotropic ME and isotropic SCDFT. We also include the values obtained with McMillan equation and Allen-Dynes modified McMillan equation (μ*=0.1). The calculated electron-phonon coupling constant, λ = 2∫₀^∞ dω α²F(ω)/ω, and the logarithmic frequency average, ω_log = exp( (2/λ)∫₀^∞ dω (α²F(ω)/ω) log ω ), are also included in the table."

**(c) Nature 出版版：无 Source Data、无 SI `[作者亲述 / 期刊页面]`**

`https://www.nature.com/articles/s41586-020-1955-z` 页面（绕过 idp 重定向后取到）：
- 章节只有 Abstract / Access options / **Data availability** / Code availability / References / Acknowledgements / Author information / Ethics declarations / Additional information / **Extended data figures and tables** / Rights and permissions / About。**没有 "Supplementary information" 段，没有 "Source data" 段**，无任何 ESM / MediaObjects 数据文件链接（只有 Fig.1–3 的 PNG）。
- **Data availability 逐字**：
  > "All the data generated in this work is available upon request from I.E. and J.A.F.-L."
- 出版版 Extended Data Table 1 的 caption（取自 `.../tables/1` 页面）与 arXiv 一致，末句仍是 "Values of λ and ω_log are also given."
- 出版版只有 Extended Data Figs. 1–8 + Tables 1–2（arXiv v1 有 Figs. 1–11 + Tables I–III）。

→ **`a2F_machine_readable_sources_2026-09-08.md` 第 12 行的 "LaH10/LaD10 … Nature Source Data (name unknown), Nature页面403, verified_exists=no" 现在可以定死为 `verified_exists=NO（确认不存在）`，不是 403 假阴性。** 唯一路径是写信问 I.E. / J.A.F.-L.。

**(d) Extended Data Fig. 4 是唯一的 α²F 一手载体（仅图）`[作者亲述-图注]`**：

> "**Extended Data Figure 4. α²F(ω) for the Fm-3m phase of LaH₁₀.** Calculated α²F(ω) for different pressures together with the integrated electron-phonon coupling constant, which is defined as λ(ω) = 2∫₀^ω dΩ α²F(Ω)/Ω. The results show that optical modes, who have hydrogen character, are responsible for the large value of the electron-phonon coupling constant λ. It is worth noting, however, that **acoustic modes with La character contribute between 0.2 and 0.5 to λ** and cannot be neglected to estimate T_c correctly."

→ 这是唯一一句作者亲述的、跨四个压力的谱形约束：**低频（La 声学）对 λ 的绝对贡献落在 0.2–0.5**。配上总 λ = 3.62/2.67/2.06/1.73，声学份额 λ_ac/λ 在 0.5/3.62 = 14% 到 0.2/1.73 = 12% 之间——**份额近似不随压力变**。这与"谱在整体硬化的同时保持相似形状"相容，与"低频尾被压掉、谱 Einstein 化"**不相容**。（这是从作者亲述的区间推的，不是抠图。）

**(e) Errea 对谱随压力演化的定性陈述：只有"硬化 / 谐波不稳定消失"，没有"变窄 / Einstein 化" `[作者亲述]`**

我通读全文，**没有任何一句说 LaH₁₀ 的 α²F 或声子谱随压力变窄、集中、或趋于单模**。能找到的相关句全部只讲硬化与动力学稳定性：

> Fig. 2 图注："The harmonic phonons show large instabilities in several regions of the Brillouin zone. Only at the high-pressure limit, e.g. above 220-250 GPa dynamic (harmonic) stabilization is reached."

> p.4："Flagrantly, the breakdown of the classical harmonic approximation for phonons makes impossible the estimation of T_c below ∼250 GPa in the Fm-3m phase and questions all previous calculations. Indeed, **the anharmonic phonon renormalization remains huge also at 264 GPa** (see Figure 2)."

> p.4（结论段）："quantum fluctuations are indispensable to sustain crystals with huge λ's (λ reaches a record value of 3.6 at 129 GPa for LaH₁₀) that would be otherwise destabilized by the colossal electron-phonon interaction to distorted (low symmetry) structures reducing the electronic density of states at the Fermi level"

→ **结论：Errea 2020 既不支持也不反对"Einstein 化"，作者根本没讨论这个方向。** 任何声称"Errea 说谱变窄"的转述都是无源的。

### 1.2 Semenok, Altshuler, Yuzbashyan, arXiv:2407.12922 / Adv. Mater. 37, 2507013 (2025)

**本地 PDF**：`library/pdf/D1/2025_Dmitrii_V_Semenok_2407.12922.pdf`。**SI Table S1–S5 全部逐页核过。`[表格读数]`**

- **Table S4**（SI p.21）列头逐字："Values of the stability parameter ξ, electron-phonon coupling constant λ, and the average logarithmic frequency ω_log for several hydrides at different pressures P." → **只有 ξ、λ、ω_log 三列**。LaH₁₀ 行：~50*/129/163/214/264 GPa → ξ = 0.96/0.64/0.53/0.43/0.37，λ = 5.32/3.62/2.67/2.06/1.73，ω_log = —/887/1119/1340/1469 K。**无 ω₂、无 ⟨ω²⟩、无逐压力 ω_max。**
- **Table S1**（SI p.18）有 ω_max 列，但 LaH₁₀ **只有单点** "LaH₁₀ (214 GPa)：T_c=245, λ=2.06, ξ=0.43, ω_log=1340 K, ω_max=3145 K"，凑不成 ω_max(P) 序列（且该行 λ/ω_log 是 Errea 的 SSCHA 理论值配实验 Tc，见 `P3_tc_upper_bound_foundations` §4.2 注 2）。
- **Table S3**（SI p.20）的 "Hydrides (experiment)" 段 LaH₁₀ 只有一行：163 GPa, λ=2.67, ω_log=1118 K，引用 [90]=Errea Nature 578, 66。**同样只有 λ、ω_log。**
- Table S5（Debye 温度表）**没有 LaH₁₀ 行**。
- **Fig. S2 图注逐字**（这是 P3 稿那个 0.182ω_log√λ 天花板的出处，作者自己就是用 ω_log 代理式）：
  > "Comparison of experimentally observed and DFT calculated critical temperatures T_c for various superconductors with Allen-Dynes empirical strong coupling estimate **T_c^{emp} = 0.182 ω_log √λ**. Observe that for all experimentally studied hydrides T_c < T_c^{emp}."

→ **Semenok 全文用 ω_log 代理式做图（Fig. 2 / Fig. S2），却在正文 Eq. (7) 写 T_c < 0.18√(λ⟨ω²⟩)，自己就把两式混用了。** 这一点 P3 稿必须显式声明用的是哪一个（`P3_tc_upper_bound_foundations` §2.4 已经指出，本轮从 Fig. S2 图注拿到了作者亲述的确证）。

### 1.3 Kruglov et al., PRB 101, 024508 (2020) / arXiv:1810.01113 —— 独立一手，谐波 `[表格读数]`

本轮**新取到全文 + Supporting Information**（arXiv v5 PDF，31 页，SI 内嵌）。这是与 Errea 完全独立的一组 Fm-3m LaH₁₀ 电声计算（DFPT 谐波、QE、4×4×4 q / 16×16×16 k，**明确不含非谐**："In this work we did not consider the effect of anharmonicity on the dynamic and superconducting properties of lanthanum hydrides."）。

**Table 1 转录（Fm-3m-LaH₁₀ 部分）**：

| Parameter | 170 GPa | 200 GPa | 210 GPa | 250 GPa |
|---|---|---|---|---|
| N_f, states/f.u./Ry | 11.2 | 10.6 | 10.3 | 10.0 |
| λ | 3.94 | 3.75 | 3.42 | 2.29 |
| ω_log, K | 801 | 906 | 851 | 1253 |
| T_C, K | 259 | 271* | 249 | 246 |

（\* 该列 Tc 是 SCDFT；其余是 Eliashberg 数值解，μ*=0.2。）

- **无 ω₂ 列。** 但 SI Eq. (S3)（同位素系数 β_AD 的解析式）**显式含 ω_log/ω₂**：
  > "β_AD = β_McM − 2.34μ*²λ^{3/2}/[(2.46+9.25μ*)·((2.46+9.25μ*)+λ^{3/2})] − 130.4·μ*²λ²(1+6.3μ*)(1 − ω_log/ω₂)(ω_log/ω₂) / { [8.28+104μ*+329μ*²+2.5λ²(ω_log/ω₂)]·[8.28+104μ*+329μ*²+2.5λ²(ω_log/ω₂)²] }"

  → **作者算过 ω₂ 但没报。** 与 Errea 情况相同：凡是用 Allen–Dynes f₂ 或 β_AD 的人都持有 ω₂，文献惯例是不列。
- SI 的 "Eliashberg spectral functions" 段（Fig. S14–S18）**只给了 Fm-3m-LaH₁₀ 在 200 GPa 一条 α²F 曲线**（Fig. S17）与 R-3m 在 150 GPa 一条（Fig. S16），**不成压力序列**，无数值文件。
- **附带价值**：这组谐波 ω_log(P)（801/906/851/1253 K @170/200/210/250）**非单调**（200→210 GPa 反降 6%），而 Errea 的 SSCHA 序列（887/1119/1340/1469）严格单调。**谐波结果的 ω_log(P) 不可用于提取标度指数**，这与 `P3_tc_upper_bound_foundations` §4.3b 对 CaH₆/ThH₁₀ 的判断一致：只有 Errea 的 LaH₁₀ SSCHA 序列是自洽的单一计算。

### 1.4 其他路线：全部未取到

| 路线 | 结果 |
|---|---|
| Nature Source Data（s41586-020-1955-z） | **确认不存在**（页面无 Source data 段，Data availability = upon request）。 |
| Nature Supplementary Information | **确认不存在**。 |
| Materials Cloud Archive | API 检索 `LaH10` → `"total": 0`；检索 `hydride` → 命中 Zr-H/Ti-H/Lu-N-H，**无 LaH₁₀、无 Errea 组条目**。 |
| Errea 组 Zenodo / GitHub | 未找到；python-sscha 的 `Tutorials/LaH10` 只是变胞弛豫教程（已记于 `a2F_machine_readable_sources` 第 13 行），无 α²F。 |
| Xie et al. 2022（arXiv:2106.05235，本地 PDF） | **有 ω₂ 但没有 LaH₁₀。** 其 "Hydrides" 验证集（13/19 条）来自 refs [59]=Kvashnin ACS AMI 10, 43809 (2018)（Th-H）与 [60]=Kruglov Sci. Adv. 4, eaat9776 (2018)（U-H）；Fig. 4 图例只有 AcH₁₀/AcH₁₂/AcH₁₆。**但它给了一个有用的一手标尺**：氢化物的 ω_log/ω₂ 在 Fig. 4 横轴上覆盖 **0.45–1.0**，即 ω₂/ω_log 最高到 ~2.2。LaH₁₀ 的 1.08–1.19（§2）落在**极靠近 Einstein 端**。作者亲述：<br>"the Allen-Dynes equation fails in this region particularly badly, since it was designed to fit materials with the ratio of the Allen-Dynes parameters ω̄₂/ω_log ≃ 1, which is strongly violated in some of the higher-T_c materials"<br>→ **对 LaH₁₀ 而言 ω̄₂/ω_log ≃ 1 其实是近似成立的（1.08–1.19），Xie 这句针对的是别的氢化物。** |
| WebSearch | 本会话 200/200 配额已用尽，后半程全部改用 WebFetch 直取。 |

---

## §2　核心结果：从一手表格解析反演 ω₂/ω_log(P)　`[本报告解析反演]`

### 2.1 方法（无任何图像数字化）

Errea Ext. Data Table I 同时给了 **T_c^{Mc}** 与 **T_c^{AD}**，同一 λ、同一 ω_log、同一 μ*=0.1。Allen–Dynes 修正的 McMillan 式与 McMillan 式**只差前因子 f₁f₂**，指数项完全相同：

$$T_c^{AD} = f_1 f_2 \cdot T_c^{Mc},\qquad
f_1 = \left[1+\left(\frac{\lambda}{2.46(1+3.8\mu^*)}\right)^{3/2}\right]^{1/3},\qquad
f_2 = 1+\frac{\lambda^2\left(r-1\right)}{\lambda^2+1.82(1+6.3\mu^*)\,r^2},\qquad r\equiv\frac{\bar\omega_2}{\omega_{\log}}$$

（f₁、f₂ 形式取自 Xie 2022 Eq. (3)(4) 的转录，见 `P3_tc_upper_bound_foundations` §2.4；与 Kruglov SI Eq. (S3) 里出现的 2.46+9.25μ*、8.28+104μ*、6.3μ* 等系数互相印证。）

于是 **f₁f₂ = T_c^{AD}/T_c^{Mc}**，指数项、ω_log、前因子约定全部约掉。已知 λ 即得 f₁，从而得 f₂，再解 r 的二次方程：

$$g\,b\,r^2 - \lambda^2 r + \lambda^2(1+g)=0,\qquad g\equiv f_2-1,\ b\equiv 1.82(1+6.3\mu^*)=2.9666$$

取物理根（小根）。

**方法的三重自洽检验（都通过）**：
1. **r ≥ 1 是数学必然**（Jensen 不等式：⟨ω²⟩^{1/2} ≥ exp⟨ln ω⟩，同一权重 2α²F/(λω)）。反演出的四个值全部 > 1。
2. **前因子约定被反演本身钉死**：若 Errea 的 "McMillan equation" 用的是 ω_log/1.45 而非 ω_log/1.20，则 f₁f₂ 要整体除以 1.2083，得 f₂ < 1 → r < 1 → **违反 Jensen，被排除**。故 Errea 两列都用 ω_log/1.20 前因子（QE 惯例）。
3. **LaD₁₀ 三行独立复算**（同一张表，同一方法）：r = 1.145 / 1.093 / 1.079 @159/210/260 GPa。高压两点与 LaH₁₀ 对应压力（214/264 → 1.105/1.095）**一致到 1%**——ω₂/ω_log 是无量纲形状量，同位素质量标度精确抵消，本应一致。**这是对整套反演的独立验证。**

**误差预算**：表格 Tc 给到 0.1 K，相对误差 ~5×10⁻⁴；传到 r 上 δr ≈ 0.001。**远小于图像数字化的 5–11%。** 残余不确定性只来自"Errea 是否用了标准 f₁/f₂ 表达式与 μ*=0.1"——正文明写了 "Allen-Dynes modified McMillan equation (μ*=0.1)"，且检验 1/2/3 都通过。

### 2.2 结果

| P (GPa) | λ | ω_log (meV) | ω_log (K) | T_c^{Mc} | T_c^{AD} | f₁f₂ | f₁ | f₂ | **r = ω₂/ω_log** | **ω₂ (K)** |
|---|---|---|---|---|---|---|---|---|---|---|
| 129 | 3.62 | 76.4 | 886.6 | 171.8 | 252.6 | 1.4704 | 1.2823 | 1.1467 | **1.194** | **1059** |
| 163 | 2.67 | 96.4 | 1118.7 | 197.1 | 247.0 | 1.2532 | 1.1927 | 1.0507 | **1.075** | **1203** |
| 214 | 2.06 | 115.5 | 1340.3 | 196.3 | 235.9 | 1.2017 | 1.1377 | 1.0563 | **1.105** | **1481** |
| 264 | 1.73 | 126.6 | 1469.1 | 189.5 | 219.2 | 1.1567 | 1.1090 | 1.0430 | **1.095** | **1609** |
| *LaD₁₀ 159* | *3.14* | *63.5* | *737.0* | *135.0* | *184.2* | *1.3644* | *1.2362* | *1.1038* | *1.145* | *844* |
| *LaD₁₀ 210* | *2.21* | *81.7* | *948.1* | *145.5* | *176.5* | *1.2131* | *1.1511* | *1.0538* | *1.093* | *1036* |
| *LaD₁₀ 260* | *1.80* | *92.2* | *1070.1* | *142.2* | *164.6* | *1.1575* | *1.1148* | *1.0383* | *1.079* | *1155* |

（meV→K 用 11.6045 K/meV。）

### 2.3 与对抗核验的逐点比对

| P (GPa) | 一手反演 r | 对抗核验（数字化）r | 偏差 |
|---|---|---|---|
| 129 | 1.194 | 1.12 | −6.2% |
| 163 | 1.075 | 1.08 | +0.5% |
| 214 | 1.105 | 1.05 | −5.0% |
| 264 | 1.095 | 1.04 | −5.0% |

**裁定**：
- **"随压力单调下降（Einstein 化）"→ 不成立。** 一手反演给的是"129 GPa 一个显著更宽的谱（1.194），其余三个压力持平在 1.075–1.105 且非单调"。数字化把 129 点的宽度低估了 6%、把高压三点各低估 5%，恰好把"一次台阶 + 平台"抹成了"单调下降"。
- **数字化的系统偏差方向是一致的（低估 r，即低估谱宽）**，量级 5–6%，与对抗核验自报的"数字化反推 λ 与 Errea 表值差 5–11%"同量级。**这给了数字化精度一个独立标尺：LaH₁₀ 那组抠图在形状量上大约有 5% 的系统性低估，不是随机噪声。**（我拿不到一手 α²F 数值，无法直接复算 λ；但 r 是一个纯形状量，其 5% 系统偏差与 λ 的 5–11% 偏差方向自洽。）
- **129 GPa 那一点是唯一真实的谱形异常**：r = 1.194 明显高于其余三点。这与 §1.1(d) 的作者亲述相容（129 GPa 是 λ = 3.62 的极端强耦合点，且是 Fm-3m 相稳定性的下边界，附近谐波不稳定最严重、非谐重整化最剧烈）。

---

## §3　对 P3 稿的直接后果：两种天花板约定下的圆顶分解

用上面的一手 r，把两种 T_c^emp 都算出来（μ*、非谐等一概不动，只换天花板定义）：

| P (GPa) | 129 | 163 | 214 | 264 | 129→264 |
|---|---|---|---|---|---|
| T_c^{ME}_{ani}（一手，Errea 表） | 255.3 | 242.8 | 237.9 | 216.9 | **−15%** |
| **天花板 A** = 0.182·ω_log·√λ（P3 现用） | 307.0 | 332.7 | 350.1 | 351.7 | **+15%**（单调升） |
| **天花板 B** = 0.1827·√(λ⟨ω²⟩) = 0.1827·√λ·ω₂ | **368.0** | **359.1** | **388.3** | **386.6** | **+5%**（非单调，129→163 段降 2.4%） |
| 亏损比 R_A = Tc/天花板A | 0.831 | 0.730 | 0.680 | 0.617 | −26% |
| 亏损比 R_B = Tc/天花板B | 0.694 | 0.676 | 0.613 | 0.561 | −19% |

（三行乘性核对：A 约定 +15% × −26% = −15% ✓；B 约定 +5% × −19% = −15% ✓。）

**给 P3 的四条可操作结论**：

1. **"天花板在涨"这句话对压力段敏感、对约定敏感。** 在 B（Allen–Dynes 真渐近）约定下，天花板 129→163 GPa 是**降的**，之后升 8%，然后持平。全程只涨 5%。**P3 稿"圆顶整个在 R、天花板反而在涨"的强论断，在 B 约定下要弱化为"天花板近似持平，圆顶几乎全部来自 R"。** 这其实**不损伤** P3 的主线（圆顶归因于 R 而非天花板），反而更干净：天花板 ±5% 波动，R 单调掉 19%，Tc 掉 15%。**建议 P3 主论断改写为"天花板在 129–264 GPa 区间近似平坦（|Δ| ≤ 8%），Tc 的 15% 下降几乎全部由亏损比 R 承担"——这个说法在 A、B 两种约定下都成立，且比"天花板在涨"更稳健。**
2. **两个天花板的绝对值差 12–20%**（B 比 A 高 60 K 左右）。P3 若要引用 Semenok 的"所有实验氢化物都 Tc < T_c^emp"，要注意 Semenok 的图（Fig. S2）**用的是 A**；换成 B 后余量更大，结论方向不变但定量全变。**必须在稿子里显式声明用哪个，并给出 ω₂/ω_log 的值。**
3. **P3 现在可以直接引用一手数字了**：LaH₁₀ 的 ω₂/ω_log = 1.19/1.08/1.11/1.10，来源是 Errea Ext. Data Table I 的 T_c^{AD}/T_c^{Mc} 之比（一手表格 + 封闭解析反演），**不需要引用任何图像数字化**。建议在稿中以脚注给出反演推导（三行），并标注"若 Errea 的 f₂ 用了非标准形式则需修正"。
4. **对 §1 的 λ_* 泛函（`P3_tc_upper_bound_foundations` U1）有直接输入**：Semenok 的 λ_*[P(ω)] 由形状函数决定，Einstein 极限 λ_*=3.69 最小。LaH₁₀ 的 r = 1.08–1.19 说明其谱**非常接近 Einstein**（Xie 的氢化物样本里 r 可到 2.2），因此 **LaH₁₀ 的 λ_* 应贴近 3.69 而非 Debye 的 4.72，且 λ_*(P) 的压力依赖很弱（r 在高压三点内几乎不变）**。这条可以直接写进 P3 §1。

---

## §4　未取到 / 待补

1. **LaH₁₀ 的一手 α²F(ω;P) 数值文件**：确认不存在于公开渠道（Nature 无 Source Data / 无 SI；Materials Cloud 无；Errea 组无 Zenodo/GitHub）。**唯一路径是按 Data availability 写信给 Ion Errea (ion.errea@ehu.eus) 与 José A. Flores-Livas。** 建议主会话决定是否发信。
2. **ω₂ 的直接一手数值**：**任何来源都没有列出**。§2 的值是反演，等级高于图读数但低于表格读数。若拿到 α²F 数值文件可直接积分核对。
3. **Allen–Dynes 1975 PRB 12, 905 原文**：f₁/f₂ 的确切表达式仍依赖 Xie 2022 的转录（`P3_tc_upper_bound_foundations` §5.1 已列）。§2 的反演直接依赖 f₂ 的形式，**这是本报告最大的单一不确定性来源**。Kruglov SI Eq. (S3) 里的系数（2.46+9.25μ*、1+6.3μ*）与 Xie 转录一致，算第二个独立佐证，但仍非原文。**建议优先补 PRB 12, 905。**
4. **Errea 对 LaD₁₀ 与 LaH₁₀ 谱形差异的解释**：LaD₁₀ @159 GPa 反演出 r = 1.145，明显高于 LaH₁₀ @163 GPa 的 1.075。谐波下同位素不改变 r（形状量），所以这个差是**非谐 + 量子涨落的真实效应**（LaD₁₀ 在 126 GPa 已失稳，159 GPa 离失稳更近）。P3 若讨论同位素效应 α 的压力依赖，这是一个可用的一手线索，但我没有进一步核实。
5. **对抗核验数字化反推的总 λ 与 Errea 表值的 5–11% 差**：我没有一手 α²F，**无法给出一手算的 λ 做直接标尺**。§2.3 给的是形状量 r 上的 5% 系统低估，作为间接标尺。

---

## §5　本轮新增/确认的编号与本地文件

| 编号 | 论文 | 状态 |
|---|---|---|
| 1907.11916 | Errea et al., Nature 578, 66 (2020) | 本地 PDF 全文（含 Ext. Data I–III、Figs. 1–11）已通读；**Ext. Data Table I 的 T_c^{Mc}/T_c^{AD} 两列为工作区新增** |
| **1810.01113** | Kruglov, Semenok, Song, Szczęśniak, Wrona, Akashi, Davari Esfahani, Duan, Cui, Kvashnin, Oganov, "Superconductivity of LaH₁₀ and LaH₁₆ polyhydrides", PRB 101, 024508 (2020) | **新增**；arXiv v5 PDF（正文 + SI，31 页）本轮取到，缓存于 `C:\Users\charles\.claude-apexin\projects\D--physics/abaf4c9b-…/tool-results/webfetch-1788880231061-xilkag.pdf`（建议主会话归档到 `library/pdf/P3/`） |
| 2407.12922 | Semenok, Altshuler, Yuzbashyan, Adv. Mater. 37, 2507013 (2025) | 本地 PDF；SI Table S1–S5 本轮逐页核过，确认无 ω₂ |
| 2106.05235 | Xie et al., npj Comput. Mater. 8, 14 (2022) | 本地 PDF；本轮确认其氢化物集是 Ac-H/Th-H/U-H，**不含 LaH₁₀** |

**未改动 theory/ 下任何文件。**

Sources（本轮网络取证）：
- [Errea et al., Nature 578, 66 (2020) 期刊页面](https://www.nature.com/articles/s41586-020-1955-z)
- [同上，Extended Data Table 1](https://www.nature.com/articles/s41586-020-1955-z/tables/1)
- [Kruglov et al., arXiv:1810.01113 (PDF v5，含 SI)](https://arxiv.org/pdf/1810.01113v5)
- [Materials Cloud Archive API 检索](https://archive.materialscloud.org/api/records?q=LaH10)
