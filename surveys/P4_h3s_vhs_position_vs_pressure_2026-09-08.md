# P4：H3S（Im-3m）vHs–ε_F 距离 δ(P) 与 dTc/dP 符号的输入数据

调研 agent：researcher（P4 任务）　日期：2026-09-08

**一句话结论**：在 150–250 GPa 区间，**没有任何一手 DFT 计算记录到 vHs 峰穿过 ε_F**。同一计算协议（Mishra EPW/Figshare 数值文件）下，DOS 峰顶始终在 ε_F **下方** 0.20 eV（160 GPa）→ 0.28 eV（200 GPa），即随加压**远离** ε_F；ε_F 被钉在峰的**高能侧肩部**（ε_F 上方 ~25–50 meV 内 DOS 塌到 N(ε_F) 的 85% 以下——这正是 Mishra E_F^eff ≈ 40 meV 的来历）。唯一声称穿越的是 Bianconi & Jarlborg 的 Lifshitz 场景（P_c1≈110、P_c2≈175–180 GPa），属该文自设模型假设，未被后续 EPW/DOS 工作确认。**跨代码的 vHs 绝对位置散差（±0.2–0.3 eV，同一压力）大于一切声称的压力移动量（~0.1 eV）**——"vHs 在哪个压力穿过 ε_F"这个问题以现有公开数据**不可判定**。

---

## 任务理解

为 P4「van Hove 奇点距 ε_F 的距离（压力可调）如何决定 dTc/dP 的符号」抽取 H3S 立方相第一性原理输入：δ(P)、ε_F 附近 DOS 形状、λ(P) 与 ω_log(P)、Tc(P) 圆顶斜率。同时核实 `theory/derivations/P2_vertex_sign_criterion.md` 转述 Mishra 的两句。

**对 P2 稿两句转述的核实结论**（先行给出，细节见 §1）：
1. **"Mishra 说 H3S 的 E_F^eff 由 ε_F 下方的 vHs 特征决定（≈40–50 meV @200 GPa）"——属实，但含义要澄清**：40–50 meV 是 vHs 窄宽度设定的**电子能标**（E_F^eff，Table S3 的 ±15% 窗口半宽定义，精值 ≈0.036 eV），**不是 DOS 峰顶到 ε_F 的距离**（后者是 0.2–0.3 eV）。出处：arXiv:2608.14949 p.3 与 arXiv:2510.25902 p.1（原句见下）。
2. **"vHs 随压力靠近 ε_F 使 m 增大"——Mishra 两篇原文均无此句，且与其数据方向相反**：本地 Mishra 数值 DOS 文件显示峰 160→200 GPa **远离** ε_F；"moves towards E_F as pressure increases" 只出现在 **Villa-Cortés 2022（LDA 泛函，162→250 GPa）**。建议 P2 稿改挂 Villa-Cortés 并注明与 Mishra 数据方向相反，或删。

---

## §1 vHs 到 ε_F 的距离 δ(P)

### 1.1 本地数值文件直接读数（最高置信度，[数值文件读数]）

数据源：Mishra–Margine Figshare 数据集（DOI 10.6084/m9.figshare.31136635，对应 arXiv:2510.25902），本地已解压于
`theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/`。

**H3S 200 GPa**（`H3S-200GPa/band/H3S.dos`，QE dos，5 meV 网格；文件头：`#  E (eV)   dos(E)     Int dos(E) EFermi =   17.619 eV`）：
- ε_F = 17.619 eV；N(ε_F) = 0.557 states/eV（17.620 处行值；同目录 `EPW_anharmonic/fsr_na/H3S.dos` 文件头给 dos_EFermi = 0.587，50 meV 网格，两者差异为网格/处理链不同）。
- **DOS 峰顶 0.667 states/eV @ E = 17.335 eV，即 ε = −0.284 eV（ε_F 下方）**。
- 峰区（>0.64）横跨 −0.32 ~ −0.18 eV，不对称：低能侧缓升（从 −1.0 eV 的 0.335 缓慢爬升），高能侧在 −0.18~−0.15 eV 处快速下折、并在 **ε_F 上方 +26–31 meV** 跌破 0.85×N(ε_F)，随后单调塌落。
- 归一化形状（N̄ = N/N(ε_F)，以 0.557 为基准）：N̄(−0.1 eV)=1.14，**N̄(+0.1 eV)=0.78**，N̄(−0.5)=0.83，N̄(+0.5)=0.62，N̄(+1.0)=0.56。
- 50 meV 网格的 fsr_na 文件交叉验证：峰 0.641 @ −0.15 eV（网格限定了分辨率），±15% 窗口在 ε_F 上方于 +25~50 meV 处闭合——与 Mishra 2608 Table S3 的 E_F^eff ≈ 0.04 eV 一致。

**H3S 160 GPa**（`H3S-160GPa/EPW_anharmonic/H3S.dos`，50 meV 网格；文件头：`EFermi = 1.6333784508E+01, dos_EFermi = 5.0767830591E-01`）：
- ε_F = 16.334 eV；N(ε_F) = 0.508 states/eV。
- **DOS 峰顶 0.654 states/eV @ E = 16.134 eV，即 ε = −0.20 eV（ε_F 下方）**。
- N̄(−0.1 eV)=1.27，**N̄(+0.1 eV)=0.89**，N̄(+0.95)=0.64。

**δ(P) 结论（同一协议内）**：峰顶 −0.20 eV（160 GPa）→ −0.28 eV（200 GPa）：**随加压远离 ε_F**；同时 ε_F 上方的 DOS 缺损加深（N̄(+0.1)：0.89→0.78），即峰的高能侧边沿向 ε_F 压近——**"峰顶远离、峰沿贴紧"同时发生**，这与 Mishra 2510.25902 的 η 上升（见 §3）以及 P2 已确认的"顶角效应 160→200 GPa 增强 3×"自洽。

### 1.2 Mishra 两篇的文字陈述（[作者亲述]）

- **arXiv:2510.25902 p.1**（Mishra & Margine, "Nonadiabatic and anharmonic effects in high-pressure H3S and D3S superconductors", Annalen der Physik 538, e00553 (2026)）：
  > "…a pronounced van Hove singularity (vHs) located near the Fermi level. **The narrow width of the vHs leads to a relatively small Fermi velocity (2.1−3.7×10⁵ m/s) [24], and hence a low Fermi energy scale (εF ≈ 40−50 meV [25]).** Together with the large phonon frequencies, this results in a high adiabatic ratio η = ωph/εF…"
- **arXiv:2510.25902 p.2**：
  > "The electronic density of states (DOS) at the Fermi level (εF) remains high at both pressures, characterized by a sharp vHs peak near εF (see Figure 1). **With the increase in pressure, the DOS at εF remains virtually unchanged**, while the free-electron-like band edge shifts to lower energies, reflecting a modest broadening of the electronic bandwidth."
- **arXiv:2608.14949 p.2**（Mishra, Belli, Zurek, Margine, "Anharmonicity and Nonadiabaticity in Hydride Superconductors"）：
  > "**In H3S, a pronounced DOS peak arising from strongly hybridized H-s/S-p bands lies just below εF [33].** This rapid variation near εF makes H3S sensitive to nonadiabatic e-ph vertex corrections [27, 35]…"
- **arXiv:2608.14949 p.3**：
  > "The adiabaticity ratio is governed by η = ω_log/E_F^eff where E_F^eff is the characteristic electronic energy scale relevant to pairing. **In H3S, the van Hove feature sets E_F^eff at approximately 40−50 meV [27, 32], well below ω_log^anh = 127 meV**…"
- **arXiv:2608.14949 Table S3 / SM S4**（经 `surveys/mishra_2608_tables_2026-09-08.md` 核实）：E_F^eff(H3S) = 0.04 eV（反算精值 ≈0.036），定义原句：
  > "we define the effective electronic scale E_F^eff entering η = ω_log^anh/E_F^eff as **the half-width of the energy window around εF over which the DOS remains within ±15% of N(εF)**."

**注意（η 定义陷阱）**：2510.25902 的 η 是 **ω_ph/ε_F**（正文给 H3S 从 160 GPa 的 4.5–5.6 升到 200 GPa 的 5.0–6.2），2608.14949 的 η 是 **ω_log^anh/E_F^eff**（H3S = 3.5）。两个 η 数值不可混用（`surveys/mishra_2608_tables_2026-09-08.md` 记录的 3.5 是后者）。

### 1.3 Quan & Pickett 2016（[ar5iv HTML 逐字提取，原文出处 arXiv:1508.04491 = PRB 93, 104526 (2016)，PRB 编号未独立重验]）

计算压力 210 GPa（a = 5.6 a.u.）。关键原句：
- "two van Hove singularities (vHs) separated by 300 meV very near the Fermi level"（Sec. I）
- "two vHs 0.25–0.30 eV apart in energy, juts upward. The Fermi energy… lies very near the upper vHs"（Sec. III.1）
- **代码间散差（同一压力！）**（Appendix A）：
  - FPLO："the vHs energies lie at **−0.43 eV and −0.11 eV**, i.e. the Fermi level lies 110 meV above the upper vHs"
  - Wien2k："In Wien2k, the vHs lie at **−0.20 eV and +0.05 eV; EF lies just below the upper vHs**"
- N(E_F)：0.44（FPLO）vs 0.63（Wien2k）states/eV-f.u.（Sec. III.1）；峰值 0.70 states/eV-f.u.、"the peak value of 0.70/eV-f.u. for N(E) dropping by 37%"（Sec. IV.2）。
- **该文没有给出 δ 随压力的系列**——只有 210 GPa 单点。

### 1.4 Sano, Koretsune, Tadano, Akashi, Arita（arXiv:1512.07365 = PRB 93, 094525 (2016)，PRB 编号据早前核对未重验；本地 PDF 全文已读）

计算压力 **250 GPa**（Im-3m-H3S）与 140 GPa（P1-H2S 对照）。原句：
- Fig. 2 图注：'In (b), there is a sharp peak around the Fermi level and **the peak width is comparable with the phonon energy scale**…'
- 正文 p.7：'there is a strong energy dependence of the DOS as a consequence of the vHs around the Fermi level for the Im3m structure (Fig. 2 (b)). **This narrow peak is a characteristic common feature of the DOS of the H3S phases** since one can also observe it in the R3m structure stabilized under pressures around 180 GPa.'
- [读图，低分辨率]：Fig. 2(b)（±1 eV 放大图）中峰顶紧贴 ε_F 虚线、略偏下方（δ ≈ 0 ~ −0.1 eV），峰高 ~0.2 states/eV/atom，背景 ~0.15。
- 该文只算 250 GPa 单点，无 δ(P) 系列；但给出 vHs 对 Tc 的影响：常数 DOS ME 高估 Im-3m-H3S @250 GPa 的 Tc 达 **225 K vs 168 K（−57 K，−34%）**（能量依赖 DOS one-shot），自洽后 193 K；Table I 全套：const DOS 225 / one-shot 168 / SC 193 / SC+ZPR 202 / SC+ZPR+AH 181 K。

### 1.5 Villa-Cortés & Rhone 2022（arXiv:2103.12055；本地 PDF p.5 已读，LDA 泛函、QHA）

原句（p.5，OCR 个别虚词可能有误差）：
> "As can be seen from the DOS, **there is a van Hove singularity close to N(0)**. In order to observe its behavior more clearly, we show in Fig. 3 the behavior of N(ε) near the Fermi level, where it can be seen how **the narrow peak moves towards E_F as pressure increases**. The inset shows how N(0) increases as pressure rises up on H3S."

[读图] Fig. 3：162 与 250 GPa 两条 N(ε)（±0.4 eV 窗），峰高均 ~0.63–0.64 states/eV·atom，162 GPa 峰位在 ε ≈ −0.05~−0.10 eV，250 GPa 峰位更贴近 0；插图 N(0) 从 ~0.52（160 GPa）升至 ~0.55（250 GPa）。
同页还有能带层面陈述：'The tops of the hole-like bands near the Γ point are above the Fermi level (E_F) and they move slowly far away from E_F as the applied pressure increases.'

### 1.6 Bianconi & Jarlborg 2015（arXiv:1507.01093，[该文自设模型场景，未被后续工作确认]）

原句（ar5iv 提取）：
- "The top of this band goes from **−0.2 eV below the chemical potential to +0.1 eV** when pressure increases"（Sec. V）
- "they move above the chemical potential **above 130 GPa** pressure as shown in Fig.4"（Sec. V）
- "**two topological Lifshitz transitions at two critical pressures Pc1=110 GPa, and Pc2=175 GPa**"（Sec. VII；同文他处写 Pc2=180 GPa）

这是唯一声称 ε_F 在 150–250 GPa 内穿越带顶/vHs 的一手来源，但其定性是半经验的（PBE + ZPE 能带边估计 + 掺杂式化学势分析），后续 EPW/DOS 数值工作（Mishra、Sano、Quan-Pickett）均未采用或确认该穿越。注意其 Pc2=175 GPa 恰等于 Errea 2016 的**经典**（无量子核效应）R3m→Im-3m 对称化压——量子化后该结构转变已移到 103 GPa，150–300 GPa 内 Im-3m 全程稳定（Errea 2016 p.1："quantum nuclear motion lowers the symmetrization pressure by 72 GPa. Consequently, we predict that the Im3m phase is stable over the whole pressure range within which a high Tc was measured."），故圆顶不是结构相变。

---

## §2 ε_F 附近 DOS 形状 N(ξ)（200 GPa 主、160 GPa 对比）

全部来自 §1.1 的本地数值文件（[数值文件读数]，置信度受网格限制：200 GPa 为 5 meV 网格、160 GPa 为 50 meV 网格）：

| 量 | 160 GPa | 200 GPa |
|---|---|---|
| N(ε_F)（states/eV） | 0.508 | 0.557（band 文件）/ 0.587（fsr 文件头） |
| 峰顶 N_peak（states/eV） | 0.654 | 0.667 |
| 峰位 δ（eV，相对 ε_F） | **−0.20** | **−0.284** |
| N_peak/N(ε_F) | 1.29 | 1.20 |
| 峰形对称性 | 不对称：低能侧缓坡 ~0.5 eV，高能侧陡降 | 同左，且高能侧在 ε_F 上方 +26–31 meV 跌破 85% N(ε_F) |
| N̄(−0.1 eV) | 1.27 | 1.14 |
| N̄(+0.1 eV) | 0.89 | 0.78 |

交叉参照（不同协议，量级一致）：
- Quan-Pickett @210 GPa：N(E_F) = 0.44（FPLO）/ 0.63（Wien2k）states/eV-f.u.，峰值 0.70；[作者亲述，ar5iv 提取]
- Villa-Cortés @162/250：峰高 ~0.63–0.64 states/eV·atom [读图]，N(0) 0.52→0.55 [读图]；
- Mishra 2608.14949 Fig. S3：N(E_F) = 0.55 states/eV/f.u.（另一独立计算，经 `surveys/mishra_2608_tables_2026-09-08.md` 转录）；
- Sano @250：峰高 ~0.2 states/eV/atom、峰宽与声子能标同量级 [读图 + Fig. 2 图注引文]。

注意 **峰顶位置在 200 GPa 上三个来源三种答案**：Mishra 数值文件 −0.284 eV、Quan-Pickett Wien2k 上鞍点 +0.05 eV / FPLO −0.11 eV、Sano 读图 ~0~−0.1 eV。DOS **峰顶**位置跨代码散差 ~0.3 eV；只有"ε_F 上方 ~40 meV 内 DOS 陡塌"（Mishra E_F^eff）与"峰在 ε_F 附近且偏下方"是跨来源共同的。

---

## §3 λ(P) 与 ω_log(P)

| 来源（协议） | P (GPa) | λ | ω_log (meV) | 备注 |
|---|---|---|---|---|
| Errea 2015（PRL 114, 157004；arXiv:1502.02832，SSCHA 非谐，μ*=0.16）Table I | 200 | 1.84（anh）/ 2.64（har） | 92.9（anh）/ 90.1（har） | 经 `lit_h3s-d3s-tcP-alphaP-data-extraction` 已核对的表值 |
| 同上 | 250 | 1.71（anh）/ 1.96（har） | 101.3（anh）/ 109.1（har） | |
| Mishra & Margine 2510.25902 Table I（PBE + ZG-A 非谐，μ*_E=0.16） | 160 | 1.88（anh） | 107 | 正文："resulting in a total λ that saturates at 1.88" |
| 同上 | 200 | 1.72（anh）/ 2.29（har） | 110（anh）/ 107（har） | 正文："integrated e-ph coupling drops by about 25%, from 2.29 to 1.72"；λ^V 从 1.78（har）被非谐压到 0.71（anh） |
| Villa-Cortés 2022（LDA + QHA，谐波谱 + LMEE）正文 p.5–6 | 162 | 2.38 | 81.18 | 原句："it decreases from 2.38 to 1.13 for H3S… from pressure 162 to 250 GPa"；"the calculated ωln value raises from 81.18 to 137.49 meV going from 162 to 250 GPa" |
| 同上 | 250 | 1.13 | 137.49 | |
| Sano 2016 @250（Im-3m，AHC/SCPH） | 250 | 1.83（har，一阶 Hermite-Gaussian）；2.06（含 ZPR）；1.83→SCPH 后 | ω_ln = 1521 K（1057 cm⁻¹，har）→ 1515 K（anh） | 原句 p.8–9、p.12 |

对数斜率（d ln x/dP，用于 P4 的"声子硬化背景项"）：

- **d ln λ/dP 全部为负**：Errea 2015 anh (200→250)：−1.5×10⁻³/GPa；Mishra 2510 anh (160→200)：−2.2×10⁻³/GPa；Villa-Cortés (162→250)：−8.5×10⁻³/GPa（LDA 谐波，幅度最大）。
- **d ln ω_log/dP 全部为正**：+1.7×10⁻³ / +0.7×10⁻³ / +10.5×10⁻³ /GPa（同上顺序）。
- 两项**部分抵消**——Errea 2016 p.4 原句（[作者亲述]）直接点明圆顶机制：
  > "This explains the smooth decrease of Tc with increasing pressure. **Between approximately 130 and 150 GPa the increase in λ is compensated by the decrease in the average phonon frequency and Tc saturates.**"
- 声子硬化的一手数字：2510.25902 p.2："a noticeable hardening of the high-frequency H-derived modes, which shift upward by approximately 22 meV"（160 GPa，Im-3m 稳定化）；200 GPa 处 H 模再硬化 ~25 meV、D3S 仅 ~14 meV（同页）。
- **来源冲突（不裁决）**：Villa-Cortés 的 250 GPa ω_ln = 137.5 meV 显著高于 Errea 2015 的 101.3 meV（同为非谐口径近似），且其 d ln λ/dP、d ln ω_log/dP 幅度比 Errea/Mishra 大 3–5 倍。三方 protocol（LDA/QHA/谐波谱 vs PBE/SSCHA vs PBE/ZG-EPW）不同，无法归一；引用斜率时建议以 Errea 2015 / Mishra 2510 两组 PBE 非谐值为下界、Villa-Cortés 为上界。

---

## §4 Tc(P) 圆顶的斜率与峰位

**峰位（[作者亲述]）**：arXiv:2510.25902 p.2：
> "**The superconducting Tc reaches its maximum near 150−160 GPa** and gradually decreases with further compression."

（与 `surveys/h3s_d3s_dome_offset_fifth_factor_2026-09-08.md` 的结论"H3S 可溯源峰位 ~150–155 GPa [文本]"一致。）

**实验点与斜率**（Drozdov 2015，arXiv:1506.08190，经 `surveys/h3s_tc_of_p_alpha_of_p_2026-09-08.md` 三轮核对）：
- 文本锚点：141 GPa / 195 K（Fig. 2b）；145 GPa / 170→195 K 退火后（Fig. 2a 文字）；155 GPa / 203(1) K 磁化率 onset（Fig. 4a）。
- Fig. 2c 读图：200–230 GPa → ~169–172 K [读图，±10 K]。
- 斜率：上升侧 (203−195)/(155−141) ≈ **+0.57 K/GPa**（两个文本锚点）；下降侧 (170−203)/(200−155) ≈ **−0.73 K/GPa**（分母中 200 GPa 点为读图值）。
- 更细的实验 Tc(P) 点（Mozaffari 2019、Einaga 2016、Osmond/Minkov 2022）：**本轮未取到**（不在本地库，网络检索预算已尽），见"未找到"清单。

**理论 Tc(P) 对照**：
- Errea 2016（arXiv:1512.02933 = Nature 538, 81 (2016)）Fig. 4(c) [读图]：非谐 Im-3m H3S 理论 ~215 K @140 → ~200 K @160–175 → ~190 K @225–250；理论斜率下降段 ≈ (190−205)/75 ≈ −0.2 K/GPa，比实验平缓。原句："The observed pressure-dependence of Tc is closely reproduced in our calculations for the Im3m phase, but not for the R3m phase."（p.1）
- Errea 2016 同表给 α：'α = 0.35 at 210 GPa and α = 0.40 at 155 GPa'（p.4）。
- Mishra 2510.25902 Table I：Tc^FSR(anh) 205（160 GPa）→ 192（200 GPa），即 −0.33 K/GPa（μ*=0.16）；ver-FSR 191→157，−0.85 K/GPa——**顶角修正使理论下降斜率变陡、向实验 −0.73 K/GPa 靠拢**（对 P4 有直接判别力：dTc/dP 的压力依赖里非绝热项贡献为负且随 P 增强，与 P2 的"顶角效应 160→200 增强 3×"一致）。

---

## 未找到一手来源的项

1. **任何来源的 δ(P) 连续曲线或表格**（3 个以上压力点的 vHs–ε_F 距离）——不存在；只有三点两两比较（Mishra 160/200、Villa-Cortés 162/250、Quan-Pickett 210、Sano 250），且趋势互相冲突。
2. **"vHs 穿过 ε_F 的压力"**：无一手计算声称；Bianconi-Jarlborg 的 P_c1=110 / P_c2=175–180 GPa 是其模型假设。
3. **Mozaffari et al. 2019**（Nat. Commun., H3S 相图）、**Einaga et al. 2016**（Nat. Phys. 12, 835，XRD+Tc）、**Osmond/Minkov 2022** 的细粒度实验 Tc(P) 点：本地无 PDF，本轮未取全文 [未取到全文]。
4. Quan & Pickett 的 PRB 卷页号（93, 104526）与 Sano 的 PRB 卷页号（93, 094525）：本轮据转引/记忆记录，arXiv 编号（1508.04491 / 1512.07365）已确证，投稿核对时以期刊页为准。
5. Quan & Pickett 原句经 ar5iv HTML 提取（逐字，但为二手工具转写）；本地无该 PDF。

## 最大不确定点

**跨代码/跨协议的 vHs 绝对能量位置散差（~0.2–0.3 eV）大于所有声称的压力移动（~0.1 eV）**：同一 210 GPa 下 Quan-Pickett 的两个代码给出上鞍点 −0.11 vs +0.05 eV；200 GPa 附近 Mishra 数值文件给峰顶 −0.28 eV、Sano（250）读图给 ~0~−0.1 eV、Villa-Cortés（162，LDA）给 −0.05~−0.10 eV。因此：
- "vHs 是否/在哪穿过 ε_F"在现有公开数据下**不可判定**；
- δ(P) 的**方向**依赖协议：PBE/EPW（Mishra）说峰顶 160→200 GPa 远离 ε_F，LDA/QHA（Villa-Cortés）说 162→250 GPa 靠近 ε_F——两者可能都对（不同泛函下带填充不同），也可能其一有系统误差；
- P4 若要用 δ(P) 解释 dTc/dP 符号，**不能依赖文献的 δ(P) 绝对值**，可行路径是固定单一协议（如直接用本地 Mishra 数值 DOS）做 ΔTc(δ) 的受控实验（把 N(ε) 平移 δ 后喂求解器），把"δ(P) 从哪来"降级为外部输入。
