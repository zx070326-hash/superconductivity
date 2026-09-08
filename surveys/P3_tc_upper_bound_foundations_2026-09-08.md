# P3「强耦合 Tc 上界的压力依赖」文献与数据地基

调研 agent：P3-foundations（researcher）　日期：2026-09-08
任务：核实 P3（`theory/PROGRAM.md` §2.3）所需的四组一手事实，不做推导。

**证据等级约定**
- **A（一手·原句）**：本地或抓取到的原文 PDF，附页码/公式号 + 逐字英文引文。
- **B（一手论文的编纂表）**：数字出现在某篇一手论文的自有表格里，但该表把数字归给了另一篇文献（即该论文对此数是"编纂者"）。可用，但引用时要连原始出处一起写。
- **C（二手）**：只在别人的转述里见到，原文未取到 → 标 `[二手，待核]`。
- `[读图]`：我从图里读的，非作者亲述。
- `[本报告自算]`：我从一手数字算出的导出量。

---

## §1　λ 的上界

### 1.1 结论摘要

工作区此前记的"Yuzbashyan–Altshuler 2025 λ≲4"**出处应归给 Semenok–Altshuler–Yuzbashyan（arXiv:2407.12922, Adv. Mater. 37, 2507013）而不是 PRL 那篇**。两篇分工是：

| 论文 | 给的是什么 | 数值 |
|---|---|---|
| Yuzbashyan, Altshuler, Patra, **PRL 135, 026503 (2025)**, arXiv:2409.19562 | 动理学不稳定性的**存在性证明**；Einstein 声子的临界 λ | λ_c = 1.27（Einstein，动理学）/ 1.4788（该文 Fig.2 谱）；C_el<0 阈值 λ > 3.69 |
| Semenok, Altshuler, Yuzbashyan, **Adv. Mater. 37, 2507013 (2025)**, arXiv:2407.12922 | 由 C_el>0 导出 λ 与 Tc 的**定量上界**，并与实验数据对照 | λ_* = 3.69（Einstein）/ 4.72（Debye）；"λ in actually existing metals is limited to λ < 4" |

### 1.2 界是多少、来自什么物理（等级 A）

**物理来源：电子比热 C_el < 0 → 电子-晶格热平衡的动理学不稳定 → 一阶结构相变。不是晶格软化，也不是 P/CDW，也不是极化子。**

arXiv:2409.19562（PRL）摘要，逐字：

> "We show that thermal equilibrium between conduction electrons and phonons becomes kinetically unstable when the renormalized electron-phonon coupling exceeds a certain threshold. We prove that negative electronic specific heat, C_el < 0, is sufficient to trigger the instability. Specifically, the instability sets in as soon as the quasiparticle weight becomes negative over a range of energies, even before C_el turns negative."

同文 p.1 右栏，明确与 P/CDW 划界（逐字）：

> "This instability is fundamentally distinct from the P/CDW instability also because it is a nonequilibrium phenomenon, does not involve phonon softening, and is not based on effective electron-phonon models."

同文 p.6，参考文献 [55] 处的注（逐字）——**这是两个阈值数字的一手出处**：

> "For example, for Einstein phonons, the kinetic equation becomes unstable for λ > 1.27, while C_el turns negative only for λ > 3.69 [33]."

同文 p.3 右栏，Einstein 声子的具体计算（逐字）：

> "Consider, for example, Einstein phonons, α²F(ω) = (λΩ/2)δ(ω − Ω), at T = 0.183g > T_c. Using Eq. (5), we find that Z_0(E) touches the horizontal axis at E_min = 1.44Ω and λ_c = 1.27."

同文 p.3，Fig. 2 说明（逐字）：

> "Results of the linear stability analysis: generalized eigenvalues γ vs the electron-phonon coupling λ at temperature T_0 = 0.20g. The metal is unstable at least when γ < 0. At this temperature, the instability takes place at λ_c = 1.4788."

### 1.3 在什么假设下成立（等级 A）

arXiv:2407.12922 §3、SI §S2 给出的稳定性判据（Eq. 3、Eq. 5、Eq. S12）：

$$\xi \equiv \max_T\left\{\int_0^\infty g\!\left(\frac{\omega}{2\pi T}\right)\frac{2\alpha^2F(\omega)}{\omega}\,d\omega\right\} < 1,\qquad
\lambda = \int_0^\infty \frac{2\alpha^2F(\omega)}{\omega}d\omega,\qquad
\lambda_* = \frac{1}{\max_T\left\{\int_0^\infty g\!\left(\frac{\omega}{2\pi T}\right)P(\omega)d\omega\right\}}$$

其中 g(x) = 6x² + 12x³ Im ψ'(ıx) + 6x⁴ Re ψ''(ıx)（Eq. 1），P(ω) ≡ (1/λ)·2α²F(ω)/ω 是归一化的**形状函数**（Eq. S10–S11）。

SI §S4 p.16 逐字：

> "The function g(x) reaches its single maximum g_max = 0.2709 at x_max = 0.3273. The upper bound λ_* is smallest, λ_* = 1/g_max = 3.6915, for Einstein phonons, i.e., for P(ω) = δ(ω − ω_E), because in this case the denominator of Equation (S12) assumes its maximum possible value g_max."

主文 §4 p.6 逐字（Debye 值 + λ<4 的论断）：

> "The value of λ_* varies with the shape of the Eliashberg function with a minimum λ_* = 3.69 attained for Einstein phonons (Section S4). For Debye phonons (α²F ∝ ω² for ω ≤ ω_D and zero otherwise), λ_* = 4.72."

> "In particular, we obtain λ_* = 7.39 and (assuming ξ_c = 0.5) λ_c = 3.69 for YH_9 at 205 GPa (see Table S1) and λ_c = 2.36 for Debye phonons. The values of λ_c for all experimentally realized hydrides are in the range from 2.19 for YH_6 to 4.12 for (La,Y)H_10. **This indicates that stability must be the reason why λ in actually existing metals is limited to λ < 4.**"

**关键假设（逐条，全部 A 等级原句）**：
1. E_F ≫ ω_max（"We define a metal as a good conductor with a Fermi energy E_F much larger than typical phonon frequencies. This is the only essential assumption that goes into the Eliashberg theory from which we derived Equation (9)." — 主文 p.7）。
2. 声子谱作为**外部输入**、用物理（非裸）声子，不在有效模型里重整化（SI §S3：Holstein/Fröhlich 的 λ₀≈0.5 晶格不稳定是双计数假象，"the instability seen in simplified electron-phonon models is essentially an artifact of double-counting (or 'overscreening')"）。
3. 界不依赖声子谱形状可被多模绕过（主文 p.4："this upper bound cannot be circumvented by invoking multiple phonon modes with comparable coupling strengths, as this does not alter the ME expression for C_el"）。
4. ξ_c = 0.5 是**猜测不是定理**（主文 p.6："Hydrides loose stability at ξ ≤ 0.5 ... suggesting that ξ_c = 0.5 for hydrides and ξ_c < 0.5 for conventional metals"；"our conjecture that ξ_c ≤ 0.5"）。ξ_* = 1 是严格上界，ξ_c 是一阶相变点，非普适。
5. 亚稳区：λ_c < λ < λ_* 内金属只能亚稳存在。

**λ 分层（务必区分，P3 若把"λ≲4"当成一个数会出错）**：

| 量 | Einstein | Debye | 物理含义 | 出处 |
|---|---|---|---|---|
| λ_c（动理学） | 1.27 | — | 动理学方程失稳（准经典，低估） | 2409.19562 p.6 |
| λ（δE>ΔE 量子不确定度可靠区） | ≳ 3.40 | — | 上一行判据只有 λ≳3.4 才自洽 | 2409.19562 p.3 |
| λ_*（C_el=0，绝对上界） | 3.69 | 4.72 | 严格上界，λ>λ_* 连亚稳都不行 | 2407.12922 Eq.(6)+SI S4 |
| λ_c = λ_* ξ_c（一阶相变，取 ξ_c=0.5） | 1.85 | 2.36 | 稳定/亚稳分界，材料相关 | 2407.12922 Eq.(6) |
| 实测氢化物 λ_c 范围 | 2.19（YH₆）– 4.12（(La,Y)H₁₀） | | | 2407.12922 p.6 |

**⚠ 对 P3 最重要的一点**：λ_* **不是常数**，它由 α²F 的**形状** P(ω) 决定（Eq. S12）。压力改变谱形（LaH₁₀ 从 129→264 GPa 谱明显硬化并变窄），因此 **λ_*(P) 本身有压力依赖**——这正是 P3 要推的东西之一，Semenok 等没有给 λ_*(P) 的显式形式。

### 1.4 竞争的上界论述（等级 A / B）

- **Esterlis–Kivelson–Scalapino 2018**：见 §2，那是 Tc 的界不是 λ 的界，但机制（λ* ~ 1 处 ME 崩溃 + 双极化子）与本节直接竞争。Semenok 等在 SI §S3 明确反驳其模型基础。
- **Migdal / Eliashberg 原始的 λ₀ ≈ 0.5 裸耦合晶格不稳定**（2409.19562 p.1）："it is well known that the metallic phase of the electron-phonon models, such as Holstein, Frölich, Su-Schrieffer-Heeger, etc., becomes unstable when the bare (unrenormalized) dimensionless electron-phonon coupling λ₀ > 0.5"，但作者强调"the physical electron-phonon coupling λ ... diverges as λ₀ → 0.5"，因此**不给 λ 的界**。
- 2407.12922 参考文献里列出但我未取全文的竞争上界（全部标 `[二手，待核]`）：
  - Hofmann, Chowdhury, Kivelson, Berg, "Heuristic bounds on superconductivity and how to exceed them", npj Quantum Materials 7, 83 (2022)。
  - Trachenko, Monserrat, Hutcheon, Pickard, arXiv:2406.08129，"Upper bounds on the highest phonon frequency and superconducting temperature from fundamental physical constants"（Gao 2025 用它的 ω_max = 3680 K，见 §3）。
  - Sadovskii, arXiv:2407.03602 / Phys. Usp. 65, 724 (2022)，"Upper limit for superconducting transition temperature in Eliashberg-McMillan theory"。
  - Kiessling, Altshuler, Yuzbashyan, "Bounds on Tc in the Eliashberg theory of Superconductivity. II: Dispersive phonons" (arXiv:2409.00532) / "III: Einstein phonons"（J. Stat. Phys. 192, 93 & 94, 2025）——**这两篇是 P3 最该补的**，是同一团队的严格数学版本。
  - Combescot, "Strong-coupling limit of Eliashberg theory", PRB 51, 11625 (1995)。
  - Leavens, "A least upper bound on the superconducting transition temperature", Solid State Commun. 17, 1499 (1974)。

---

## §2　强耦合 Tc 上界

### 2.1 Esterlis–Kivelson–Scalapino 2018（等级 A，全文已取）

**arXiv:1806.00488（v1, 2018-06-01）；npj Quantum Materials 3, 59 (2018)；标题是 "A bound on superconducting T_c's"**（工作区旧记录的年份/机构对，但要注意**没有** "in the Holstein model" 的标题）。

摘要 p.1，Eq. (1)(2) 逐字：

> "Specifically, drawing on both empirical results (shown in Figure 2 below) and recent results[1] of determinant quantum Monte Carlo (DQMC) studies of the paradigmatic Holstein model, we propose that
> **k_B T_c ≤ A_max ℏω̄**　(1)
> where A_max is a dimensionless number of order one that we estimate to be
> **A_max ≈ 1/10.**　(2)"

**ω̄ 是哪个频率？** 摘要里只说 "an appropriate characteristic phonon frequency"；**§2 明确把它当 Debye 温度**（p.4，逐字）：

> "In Fig. 2 we plot the superconducting T_c vs. the Debye temperature, Θ_D for various elemental superconductors and compounds for which data is available. ... It thus represents a specific, unambiguously defined (although somewhat crude) estimate of the characteristic phonon frequency ℏω̄ ∼ k_BΘ_D. Also shown in the figure is the proposed bound T_c ≤ Θ_D/10."

**所以是 ω_D，不是 ω_0 也不是 ω_log。** 但对 H₃S 他们改用 DFT 最高声子频率（见下）。

**Holstein 模型给的具体系数**（p.4，逐字）——注意是 0.08 和 0.12 两个数：

> "Combining the results from the ME theory (where valid) with the DQMC results, we obtained estimates of A_max for the Holstein model on the square lattice. For the range of parameters we have explored, the highest value of T_c we have inferred is **0.08 times the bare phonon frequency**, but because significant phonon softening occurs for λ ∼ λ*, this value of T_c is **0.12 times the maximal renormalized phonon frequency**."

**参数区间**（p.1 Eq.(3) 上下文 + p.2）：

> "k_B T_c = ℏω̄ A(w, λ, μ, ...)　(3) ... where w ≡ ℏω̄/E_F"

> "for a given small value of w we find that the ME approximation is extremely accurate for λ smaller than a characteristic value λ* ∼ 1, while for λ > λ* (even if λ ≪ w^{-1}) ME theory is both quantitatively and qualitatively incorrect."

界的来源是 **A(w,λ,...) 在 λ = λ* ≈ 1 处有极大值**（Max[A] ≡ A_max = A(0, λ*, 0, ...)），机制是向双极化子/CDW 的 crossover，**不是**晶格不稳定。DQMC 参数：ℏω₀/E_F = 0.1，ρ = 0.8，T ≈ E_F/30（Fig.1 说明）。

**经验支持（Fig. 2，等级 A）**：

> "Measured values of T_c and Θ_D for various putatively conventional crystalline superconductors. The shaded region is the proposed bound **T_c = (0.08 − 0.12) Θ_D**."

具体材料 T_c/Θ_D（p.4，作者亲述）：Pb 0.069、Nb 0.034、Hg 0.058；Nb₃Sn 0.066、Nb₃Ga 0.071、Nb₃Ge 0.072；MgB₂ 0.04（Θ_D = 884 K）；BKBO 0.09。

**H₃S 的自洽检验（p.5，逐字）——这条对 P3 直接有用**：

> "The highest T_c of all conventional superconductors is T_c = 203 K in H₃S at 155 GPa [12]. As far as we know, the Debye temperature has not been measured; however, if we identify ω̄ with the largest phonon frequencies found in DFT calculations of the phonon band structure we obtain the estimate ℏω̄ = 0.23 eV [13]. If we accept this theoretical value, then **k_B T_c/ℏω̄ = 0.08**, i.e. it comes very close to saturating our bound."

### 2.2 与 Semenok/Yuzbashyan 界的冲突（等级 A）

2407.12922 主文 p.7 明确说 0.1ω_D 这条界**被氢化物违反**（逐字）：

> "Note that both our upper bounds (9) and (10) are significantly more generous than the heuristic bound of 0.1ω_D proposed in [10]. However, observe that, for example, (La,Ce)H_{9−10} at 123 GPa with T_c = 190 K and ω_D = 1107 K, extracted from the temperature dependence of the electrical resistance, violates the latter bound, so it apparently does not hold for compressed hydrides. Indeed, to fall below the 0.1ω_D bound, (La,Ce)H_{9−10} must have a Debye frequency of about 2000 K, which exceeds all calculated and experimental values of ω_D ever obtained for synthesized superhydrides (Table S5)."

（(La,Ce)H₉ 的 190 K / ω_D=1107 K 见 2407.12922 Table S5；该行原始来源是 refs [43,44] = Chen et al. Nat. Commun. 14, 2660 (2023) 与 Bi et al. Nat. Commun. 13, 5952 (2022)，我未取原文 → 等级 B。190/1107 = 0.172。）

**这是 P3 最有价值的一处张力**：Esterlis 的 0.1Θ_D 在常规金属上几乎饱和（Pb/Nb₃Ge 0.07），在 H₃S 上用 ω_max 定义也几乎饱和（0.08），但用 ω_D 定义在 (La,Ce)H₉ 上被违反 1.7 倍。**分歧的全部内容在"ω̄ 取哪个频率"**——ω_D vs ω_max vs ω_log 在氢化物里差 2–3 倍（见 §4 表）。

### 2.3 Semenok–Altshuler–Yuzbashyan 的 Tc 界（等级 A）

主文 §4，Eq. (7)–(11) 逐字：

> "In the formal limit λ → ∞ of the ME theory, T_c asymptotically approaches 0.18√(λ⟨ω²⟩), while at finite λ it falls below this asymptote [8, 29, 53], i.e.,
> **T_c < 0.18√(λ⟨ω²⟩)**　(7)"
> "**⟨ω²⟩ = (2/λ)∫₀^∞ α²F(ω)ω dω**　(8)"

> "Therefore, the critical temperature for Einstein phonons (T_c^E) with maximal λ = λ_* = 3.69 and frequency ω_max provides a rigorous upper bound on T_c for stable and metastable superconducting metals,
> **T_c < 0.32 ω_max**　(9)
> where we used a standard algorithm [29, 46] to numerically compute T_c^E = 0.3175 ω_max."

> "Assuming ξ_c = 1/2 (see above), Equation (6) obtains λ_c = λ_*/2. Similarly calculating T_c^E = 0.1995 ω_max for Einstein phonons with λ = 3.69/2, we determine an upper bound on T_c for stable metals,
> **T_c < 0.20 ω_max**　(10)"

> "**T_c^stable < 600 K/√A,　T_c^metastable < 950 K/√A**　(11)
> where A is the atomic mass (in atomic units) of the lightest element in the material."

Eq.(11) 的输入（p.8，逐字）：

> "According to recent ab-initio calculations of the Eliashberg function for solid atomic hydrogen at 500 GPa [48], ω_max = 3000 cm⁻¹ = 4320 K for harmonic phonons ... It is therefore safe to replace ω_max in Equation. (9) and (10) with ω₀ = 3000 K"

**⚠ 压力依赖：作者自己承认 Eq.(11) 的压力处理是最弱的一环**（p.8，逐字）——**这是 P3 的直接入口**：

> "It is important to note that ω_max grows (roughly linearly) with pressure, see, e.g., [12] and Figure S1a. For the above bounds, we used ω_max at 500 GPa—the current upper limit for transport measurements. However, there is a negative correlation between ω_max and other characteristic phonon frequencies and the electron-phonon interaction strength (Figure 2 and Figure S1b), **which we did not take into account in our analysis**, such that increasing the pressure beyond 500 GPa will likely result in an overall decrease rather than increase of T_c."

### 2.4 Allen–Dynes 1975（等级 B——原文 PRB 12, 905 未取到，无 arXiv）

**⚠ 与任务书里写的形式有出入，请以下面为准。**

任务书写的强耦合渐近是 "Tc → (1/2π)√(λ⟨ω²⟩)"，即系数 0.1592。**三处独立的一手论文都给 0.18 而不是 0.159**：

| 出处 | 系数 | 逐字 |
|---|---|---|
| Semenok 2407.12922 Eq.(7) | **0.18** | "T_c asymptotically approaches 0.18√(λ⟨ω²⟩)" |
| Yuzbashyan 2409.19562 p.3 | **0.1827** | "T_c = 0.1827 in the limit λ → ∞ [20]"（[20] = Allen & Dynes PRB 12, 905） |
| Semenok 2407.12922 Fig.2/Fig.S2 说明 | **0.182**（但配 ω_log） | "T_c^emp ≡ 0.182 ω_log √λ = const, where T_c^emp is the empirical strong coupling estimate for T_c [29]" |

**注意 0.182 ω_log√λ 与 0.18√(λ⟨ω²⟩) 不是同一式**：前者用 ω_log，后者用 √⟨ω²⟩，只有 ω_log ≈ √⟨ω²⟩ 时才重合。氢化物里两者可差 20–30%（Xie 2022 强调 ω_log/ω̄₂ 在高 Tc 材料里"strongly violated"）。P3 若混用会引入系统偏差。

**Allen–Dynes 公式的完整形式（等级 B，经 Xie et al. 2022, arXiv:2106.05235 p.2, Eq. (2)–(4) 转录）**：

$$T_c = \frac{f_1 f_2 \omega_{\log}}{1.20}\exp\left(-\frac{1.04(1+\lambda)}{\lambda - \mu^*(1+0.62\lambda)}\right)$$
$$f_1 = \left[1 + \left(\frac{\lambda}{2.46(1+3.8\mu^*)}\right)^{3/2}\right]^{1/3},\qquad
f_2 = 1 + \frac{\lambda^2\left(\frac{\bar\omega_2}{\omega_{\log}} - 1\right)}{\lambda^2 + 1.82(1+6.3\mu^*)\left(\frac{\bar\omega_2}{\omega_{\log}}\right)^2}$$

Xie 2022 p.2 对适用区间的一手判断（逐字）：

> "Allen and Dynes showed that the true Eliashberg T_c did not obey such a bound in this limit but rather grew as √λ. They proposed an alternate approximate fit to Eliashberg theory based on data on a few low-T_c superconductors known in 1975"

> "the Allen-Dynes equation fails in this region particularly badly, since it was designed to fit materials with the ratio of the Allen-Dynes parameters ω̄₂/ω_log ≃ 1, which is strongly violated in some of the higher-T_c materials"

> "while the Allen-Dynes formula accurately predicts the Eliashberg T_c for λ values near 1.6 (the coupling constant for Hg and Pb), it nevertheless deviates from the Eliashberg T_c when λ is significantly larger or smaller than 1.6 and when the shape of α²F(ω) differs from the simple unimodal Einstein model."

且 Allen–Dynes 的 f₂ 本身是拟合而非推导（Xie p.3，逐字）：

> "In the Allen-Dynes formula, the 'arbitrarily chosen' shape-dependent factor f₂ is based on the numerical solutions using the spectral functions of Hg, Pb, and the Einstein model."

**未找到一手来源**：Allen–Dynes 1975 原文的 λ→∞ 渐近式（含前因子的确切值 0.1827 与其推导条件）我没有拿到原文。三处二手一致给 0.18x，但 P3 若要把这个系数写进解析式，建议主会话补 PRB 12, 905 原文核对（该刊无 arXiv 版，需机构订阅）。

---

## §3　常压 ω_log–λ 权衡（Gao et al. 2025）

**arXiv:2502.18281（v1, 2025-02-25），"The Maximum T_c of Conventional Superconductors at Ambient Pressure"，Kun Gao, Cerqueira, Sanna, Fang, Dangić, Errea, Wang, Botti, Marques。本地 PDF：`library/pdf/D4/2025_Kun_Gao_2502.18281.pdf`（D2 目录有同一份）。等级 A（全文已读）。**

### 3.1 它给的是什么（关键：**经验的散点趋势，不是解析上界**）

摘要逐字：

> "Through analysis of electron-phonon calculations for over 20 000 metals, we critically examine this question. We find that while hydride metals can exhibit maximum phonon frequencies of more than 5000 K, the crucial logarithmic average frequency ω_log rarely exceeds 1800 K. **Our data reveals an inherent trade-off between ω_log and the electron-phonon coupling constant λ, suggesting that the optimal Eliashberg function that maximizes T_c is unphysical.**"

Fig. 2 说明（逐字）：

> "Scatter plot of ω_log versus λ for all calculated systems, the contour lines are T_c from McMillan's formula with μ* = 0.1."

正文 p.3 对 trade-off 的说明（逐字）——**明确说明这是"两个矩"的关系，不是拟合出的公式**：

> "the parameters ω_log and λ are not entirely independent, as they are two different moments of the same α²F(ω) function. ... Obviously ω_log is favored by high frequencies, while λ by low frequencies, so the shape of the optimal α²F(ω) that maximizes T_c (and that only contains a *single*, very high-frequency, flat-band phonon mode that couples very strongly to the electrons) seems to be unattainable physically."

> "We see that, as expected, compounds with very large values of ω_log often have small λ and vice-versa. Furthermore, no material comes close to the optimal case, and compounds with the highest values of T_c are the ones that achieve a good compromise between reasonable high ω_log and λ."

**所以：经验的、散点形式的负相关；论文自己没有给权衡的解析式或拟合直线。** 唯一带公式的"上界"是用代表值代进 McMillan：

> "If, for the sake of the argument, we insert λ = 2, ω_log = 1800 K into the McMillan function, and assume a value of μ* = 0.1, we obtain T_c = 260 K. This is smaller than the value of Ref. 37, but still much larger than the current record of MgB₂."

> "It is most likely, however, that this limit is unattainable for any physical system at ambient pressure."

### 3.2 数据基于哪些材料（等级 A）

> "This is done by analyzing experimentally-known superconductors and our calculations of the electron-phonon interaction and superconducting properties for more than 20 000 metals [38, 39]."

ref [38] = Cerqueira, Fang, Errea, Sanna, Marques, Adv. Funct. Mater. 34, 2404043 (2024)（= arXiv:2403.13496，本地有 PDF）；ref [39] = Cerqueira, Sanna, Marques, Adv. Mater. 36, 2307085 (2024)。数据集是**常压凸包附近的金属**，DFPT 谐波。

### 3.3 结论数字（等级 A）

- ω_max 最高：AgTl₂H₂（假想）5396 K；"the overwhelming majority of our materials have ω_max below 3000 K, we find a few that extend to 5400 K"。
- ω_log 分布："Contrary to ω_max, ω_log has a single peak at very low frequency, and decays rapidly with frequency. At the end of the tail, at ω_log values in excess of 1800 K we find a few hydrides"。
- 最优实际候选：Li₂AgH₆（Eliashberg Tc = 108.8 K，SCDFT 83.0 K，λ ≈ 4，谐波），Li₂AuH₆（Eliashberg 116.1 K，SCDFT 91.0 K）；两者都热力学不稳定（凸包上方 0.319 / 0.172 eV/atom）。
- 实际上限："it is possible to design hypothetical compound with values of T_c reaching 100–120 K, much larger than the current experimental record, but still very far from room-temperature"；"achieving room-temperature conventional superconductivity at ambient pressure is extremely unlikely"。
- Fig. 4：凸包距离 vs Tc 的 Pareto 前沿——"all compounds with high T_c appear to be thermodynamically unstable"。这与 Semenok 的 ξ<1 稳定性论证**在结论上一致但机制不同**（一个是热力学凸包，一个是电子比热动理学）。

### 3.4 对 P3 的用途与限制

Gao 2025 **只覆盖常压**，明确排除了高压相（"A possible path through the high-pressure synthesis followed by quenching to ambient pressure is commonly suggested in the literature [58–60], however it is still unproved to this date, and remains highly speculative"）。P3 要的压力依赖，Gao 提供的是 **P=0 的锚点**：ω_log ≲ 1800 K、λ ≲ 4（Fig.1 下图 λ 直方图截到 4）、trade-off 存在但无解析形式。

**Semenok Fig. S1b 给了唯一一条 ω_log–λ 的显式经验线**（等级 A，2407.12922 SI 图说明逐字）：

> "(b) Correlation between ω_log in Kelvin and the electron-phonon coupling parameter λ. For about 95% of the compounds ω_log in Kelvin falls below the line **2500 − 500λ**. Red areas indicate unlikely combinations of ω_log, P and λ."

**这条 ω_log[K] < 2500 − 500λ 是 P3 目前能拿到的最接近"ω_log–λ 权衡解析式"的东西，且与 Gao 的散点定性一致（λ=2 → 1500 K；Gao 说 ω_log 极少超 1800 K）。两者独立数据集，互相印证。**

---

## §4　氢化物 Tc 纪录包络 + Grüneisen 标度

数据表见 `theory/data/hydride_tc_record_envelope.csv`。本节只记核实过程与关键原句。

### 4.1 每条 (材料, P, Tc) 的一手核实（等级 A，除标注外）

**H₃S — Drozdov, Eremets, Troyan, Ksenofontov, Shylin, Nature 525, 73 (2015), arXiv:1506.08190（本地 PDF）**
Fig. 4(a) 说明 p.11 逐字（**这是 203 K + 155 GPa 配对的一手出处**）：
> "Temperature dependence of the magnetization of sulfur hydride at a pressure of 155 GPa in a the zero-field cooled (ZFC) and 20 Oe field cooled (FC) modes. The onset temperature is T_onset = 203(1) K."

p.3 逐字（方法与同位素）：
> "The onset temperature of the superconducting state T_onset = 203(1) K, the width of the superconducting transition is nearly the same as in electrical measurements (Fig. 4a)."
> "Comparison T_c's in the pressure range P>170 GPa (Fig. 2c) gives α_H ~0.3."

Fig. 2(c) 说明 p.9："Dependence of T_c on pressure. The data on annealed samples are presented. Open color points refer to sulfur deuteride, and filled points to sulfur hydride." → **Tc(P) 圆顶的一手来源，175–225 GPa 区间 Tc 在 190 K 附近平台，250 GPa 附近降到 ~170**（[读图]，不作定量使用）。

**LaH₁₀ — Drozdov, Kong, Minkov et al., Nature 569, 528 (2019), arXiv:1812.01561（本地 PDF）**
摘要逐字：
> "Here, we report superconductivity with a record T_c ∼ 250 K within the Fm3̄m structure of LaH₁₀ at a pressure P ∼ 170 GPa. We proved the existence of superconductivity at 250 K through the observation of zero-resistance, isotope effect, and the decrease of T_c under an external magnetic field, which suggests an upper critical magnetic field of ∼120 T at zero-temperature. The pressure dependence of the transition temperatures T_c (P) has a maximum of 250-252 K at the pressure of about 170 GPa."

p.3 逐字（圆顶形状 + 与 Somayazulu 的冲突）：
> "It clearly has a 'dome'-like shape: after the initial increase and reaching the maximum T_c = 250-252 K at ∼170 GPa, T_c decreases abruptly at higher pressures (Fig. 1, inset). This is in a clear disagreement with the claims of superconductivity with T_c = 280 K and 290 K in Ref. 23."

p.2：H_c2(0) 外推 "yields values between 95 and 136 T"（与摘要的 ~120 T 是同一区间）。

**LaH₁₀（第二组）— Somayazulu et al., PRL 122, 027001 (2019), arXiv:1808.07695（本地 PDF）**
摘要逐字（**注意：是"电阻下降起点"，不是零电阻，作者自己用的是 "evidence"/"suggest"**）：
> "We detail the synthesis and characterization, including four-probe electrical transport measurements, of lanthanum superhydride samples that display a significant drop in resistivity on cooling beginning around 260 K and pressures of 190 GPa. Additional measurements on two additional samples synthesized the same way show resistance drops beginning as high as 280 K at these pressures."
> "We suggest that the resistance drop is a signature of the predicted superconductivity in LaH₁₀ near room temperature"

→ **250 K vs 260 K 的差异在工作区 CSV 里被记为"读数定义差异"，现在可以定死：两组压力不同（170 vs 190 GPa）、判据不同（零电阻+同位素+磁场 vs 仅电阻下降起点），且 Drozdov 明确不同意 280/290 K 的说法。** 纪录包络应采用 Drozdov 的 250–252 K @ 170 GPa。

**(La,Y)H₁₀ — Semenok et al., Mater. Today 48, 18 (2021), arXiv:2012.04787（未取全文，抓取到摘要，等级 A-摘要）**
> 报道 (La,Y)H₁₀ "a maximum T_C ~ 253 K" "at 183 GPa"；同文另有立方 (La,Y)H₆，T_c = 237±5 K。

**YH₆ — Troyan, Semenok, Kvashnin et al., Adv. Mater. 33, 2006832 (2021), arXiv:1908.01534（本地 PDF）**
摘要 + Highlights 逐字：
> "Here we report the synthesis of yttrium hexahydride Im3̄m-YH₆ that demonstrates the superconducting transition with T_C ∼ 224 K at 166 GPa, much lower than the theoretically predicted (>270 K). The measured upper critical magnetic field B_c2(0) of YH₆ was found to be 116–158 T, which is 2–2.5 times larger than the calculated value. A pronounced shift of T_C in yttrium deuteride YD₆ with the isotope coefficient 0.4 supports the phonon-assisted superconductivity."
> Highlights: "High-T_C superconductivity in Im3̄m-YH₆ at 224 K and 166 Gigapascals"

**YH₉ 与 YH₆（第二组）— Kong, Minkov, Kuzovnikov et al., Nat. Commun. 12, 5075 (2021), arXiv:1909.10482（本地 PDF）**
摘要逐字（**YH₉ 243 K @ 201 GPa 与 YH₆ 227 K @ 237 GPa 的一手出处**）：
> "the pressure dependence of T_c for YH₉ also displays a 'dome like shape' with the highest T_c of 243 K at 201 GPa. We also observed a T_c of 227 K at 237 GPa for the YH₆ phase. However, the measured T_cs are notably lower by ∼30 K than predicted. Evidence for superconductivity includes the observation of zero electrical resistance, a decrease of T_c under an external magnetic field and an isotope effect."

⚠ **YH₆ 有两个不相容的一手读数**：Troyan 224 K @166 GPa 与 Kong 227 K @237 GPa。压差 71 GPa 而 Tc 只差 3 K——这意味着 YH₆ 的 Tc(P) 在 166–237 GPa 之间非常平，或两组样品化学计量不同。**P3 画包络时这一点必须显式处理**（不能把两点连成 dTc/dP）。

**CaH₆ — Ma, Wang, Xie et al., PRL 128, 167001 (2022), arXiv:2103.16282（未取全文，抓取到摘要，等级 A-摘要）**
> "measured superconducting critical temperature T_c of 215 K at 172 GPa"，判据 "a sharp drop of resistivity to zero and a characteristic decrease of T_c under a magnetic field"，WHH 外推 "giant zero-temperature upper critical magnetic field of 203 T"。

（工作区 CSV 里另记有 207 K @181 GPa / 212 K @170 GPa 两组批次读数，来源为同一论文，我未取正文核实 → 保留 `[二手，待核]`。）

**CeH₉ / CeH₁₀ — Chen, Semenok, Huang et al., PRL 127, 117001 (2021), arXiv:2101.01315（本地 PDF）**
摘要逐字——**这条推翻了工作区 CSV 现有的 CeH₉ 行**：
> "SC-I has T_c=115 K at 95 GPa, showing expected decrease on further compression due to decrease of the electron-phonon coupling (EPC) coefficient λ (from 2.0 at 100 GPa to 0.8 at 200 GPa). SC-II has T_c = 57 K at 88 GPa, rapidly increasing to a maximum T_c ∼100 K at 130 GPa, and then decreasing on further compression. This maximum of T_c is due to a maximum of λ at the phase transition from P6₃/mmc-CeH₉ into a symmetry-broken modification C2/c-CeH₉."

其中 SC-I = Fm3̄m-CeH₁₀，SC-II = P6₃/mmc-CeH₉。
**⚠ 现有 `hydrides_experimental.csv` 第 22 行写 "CeH9, 88-116 GPa, 90-100 K"、第 23 行写 "CeH10, 95 GPa, 115 K"：CeH₁₀ 那行对，CeH₉ 那行应改为"最大 ~100 K @130 GPa"。**
**这条对 P3 特别有价值**：作者亲述 CeH₁₀ 的 λ 从 100 GPa 的 2.0 降到 200 GPa 的 0.8（λ ∝ P^{-0.72}[本报告自算]），且 CeH₉ 的 Tc 峰**不是**声子硬化/λ 竞争的结果而是**结构相变处 λ 的极大**。P3 若把 Tc(P) 圆顶一律归给 "ω↑ vs λ↓" 的连续竞争，CeH₉ 是反例。

**ThH₉ / ThH₁₀、(La,Ce)H₉₋₁₀ — 未取原文，等级 B/C**
Semenok 2407.12922 Table S1 给 ThH₉ 146 K @150 GPa、ThH₁₀ 161 K @170 GPa、(La,Ce)H₉ 190 K @123 GPa，原始来源分别是 [61]=Semenok et al. Mater. Today Phys. (2020)、[62]、[43,44]=Chen Nat.Commun.14,2660(2023) / Bi Nat.Commun.13,5952(2022)。**我未取这些原文** → CSV 中标 `二手(Semenok Table S1 编纂)`。

### 4.2 λ、ω_log、ω_max 的来源分层

**Semenok 2407.12922 Table S1 是目前唯一一张把 (Tc, λ, ξ, ω_log, ω_max) 放在一起的表**（SI p.18）。原样转录（实验材料段）：

| Compound | T_c, K | λ | ξ | ω_log, K | ω_max, K |
|---|---|---|---|---|---|
| H₃S (157 GPa) | 190 | 1.84 | 0.32 | 1078 | 2704 |
| LaH₁₀ (214 GPa) | 245 | 2.06 | 0.43 | 1340 | 3145 |
| YH₆ (165 GPa) | 224 | 1.71 | 0.39 | 1333 | 2450 |
| ThH₉ (150 GPa) | 146 | 1.73 | 0.32 | 957 | 3362 |
| ThH₁₀ (170 GPa) | 161 | 1.65 | 0.31 | 1116–1470 | 3122 |
| YH₉ (205 GPa) | 235 | 2.66 | 0.36 | 916–1054 | 3167 |
| (La,Y)H₁₀ (180 GPa) | 253 | 3.87 | 0.47 | 868 | 3119 |
| (La,Ce)H₉ (123 GPa) | 190 | 2.27 | 0.33 | 582–915 | 2894 |
| CaH₆ (172 GPa, **DFT 段**) | 215 | 2.69 | 0.58 | 950 | 2877 |

**⚠ 三处必须注意的问题（我核对出来的）：**

1. **H₃S 那行的压力标错了。** λ=1.84、ω_log=1078 K（=92.9 meV）与 Errea, Calandra, Pickard et al. PRL 114, 157004 (2015), arXiv:1502.02832, Table I 的 **200 GPa 非谐 SSCHA 值（λ_anh=1.84, ω_log=92.86 meV, Tc=194 K）逐位吻合**（见 `theory/data/isotope_curated.md` 第 20 行，该行经 researcher 读 PDF 核实）。Semenok 标 "157 GPa"、Tc 标 190 K。**这是 200 GPa 的理论值被贴上了 157 GPa 的标签。** P3 用这行时压力要按 200 GPa 记。

2. **LaH₁₀ 那行是理论值不是实验值。** Table S4 的 LaH₁₀ 序列（129/163/214/264 GPa → λ=3.62/2.67/2.06/1.73，ω_log=887/1119/1340/1469 K）与 Errea et al., Nature 578, 66 (2020), arXiv:1907.11916, Ext. Data Table I 的 SSCHA 值（76.4/96.4/115.5/126.6 meV = 886.6/1118.7/1340.4/1469.2 K）**逐位吻合**（见 `isotope_curated.md` 第 24 行）。所以 Table S1/S4 的 LaH₁₀ 行是 **Errea 2020 的 SSCHA 计算**，Table S1 却把它放在 "Experimental Materials" 段并配 Tc=245 K。**λ/ω_log 是理论、Tc 是实验，混在一行。**

3. **CaH₆ 在 Table S1 里明确在 DFT 段**，而 Table S3 又把 CaH₆ (170 GPa, λ=2.69, ω_log=950) 列在 "Hydrides (experiment)" 段。Semenok 自己在 p.6 说 "We attribute ξ = 0.58 in CaH₆ and MgH₆ to the fact that these are not actually existing hydrides at these pressures"——**但 CaH₆ 在 172 GPa 是被 Ma et al. 2022 合成并测到 215 K 的**。这处内部不一致我无法解释，标为关键不确定点。

**结论：Table S1/S3 里没有一个 λ 或 ω_log 是实验直接测的，全部是 DFT/SSCHA 计算值，只是 Tc 和压力来自实验。** CSV 里我按此分列。

### 4.3 ω(P) 的压力标度：能用的一手数据

**(a) LaH₁₀，最完整的一手 ω_log(P) 序列（等级 A，Errea 2020 Ext. Data Table I，经 workspace `isotope_curated.md` 核实；与 Semenok Table S4 独立吻合）**

| P (GPa) | 129 | 163 | 214 | 264 |
|---|---|---|---|---|
| ω_log^anh (meV) | 76.4 | 96.4 | 115.5 | 126.6 |
| λ^anh | 3.62 | 2.67 | 2.06 | 1.73 |
| Tc^aniso-ME (K), μ*=0.1 | 255.3 | 242.8 | 237.9 | 216.9 |

由此算出的标度指数（**全部 [本报告自算]**，定义 n ≡ d ln ω_log / d ln P、m ≡ d ln λ / d ln P，分段用相邻两点的对数斜率）：

| 区间 (GPa) | n = dlnω_log/dlnP | m = dlnλ/dlnP |
|---|---|---|
| 129 → 163 | **0.99** | −1.05 |
| 163 → 214 | **0.66** | −0.95 |
| 214 → 264 | **0.44** | −1.13 |
| 129 → 264（全程） | **0.70** | **−1.03** |

→ **λ ∝ P^{−1.0} 在整个 135 GPa 区间惊人地稳定；ω_log 的硬化指数从 ~1.0 单调降到 ~0.44（趋于饱和）。** 这两条是 P3 目前最硬的压力标度输入。

由同一组数算出的 Allen–Dynes 经验界 T_c^emp = 0.182 ω_log√λ（**[本报告自算]**，用 Semenok 的系数约定）：

| P (GPa) | 129 | 163 | 214 | 264 |
|---|---|---|---|---|
| T_c^emp = 0.182 ω_log√λ (K) | 307 | 333 | 350 | 352 |
| Tc^aniso-ME (K) | 255 | 243 | 238 | 217 |
| 比值 Tc/T_c^emp | 0.83 | 0.73 | 0.68 | 0.62 |

→ **界随压力单调升高并饱和（307→352 K），而实际 Tc 单调下降（255→217 K），两者的间隙从 17% 张到 38%。** P3 的核心问题——"Tc_max(P) 的形状"——在 LaH₁₀ 上表现为：**上界不解释圆顶，圆顶来自 Tc/T_c^emp 这个比值的压力依赖（即 μ*、非谐、谱形/f₂ 因子的压力依赖），而不是来自 λ⟨ω²⟩ 本身。** 这条我建议作为 P3 的第一个可证伪命题。

**(b) ThH₁₀、CaH₆、CeH₉ 的 ω_log(P) / ξ(P)（等级 B，Semenok Table S4，原始来源未追）**

| 材料 | P (GPa) | ξ | λ | ω_log (K) |
|---|---|---|---|---|
| ThH₁₀ | 100 / 200 / 300 | 0.46 / 0.32 / 0.25 | 2.57 / 1.58 / 1.33 | 1048 / 1184 / 1150 |
| CaH₆ | 150* / 160 / 170 / 180 / 190 | 0.58 / 0.41 / 0.40 / 0.39 / 0.39 | 2.97 / 1.96 / 1.89 / 1.83 / 1.84 | 964 / 1204 / 1419 / 1249 / 1333 |
| CeH₉ | 100 / 120 / 150 / 200 | 0.175 / 0.165 / 0.152 / 0.131 | （表中空） | （表中空） |
| LaH₁₀ | ~50** | 0.96 | 5.32 | — |
| H₃S | 135 / 157 | 0.42 / 0.39 | （空） | （空） |

\* 表注："Without accounting for the anharmonicity"　\*\* 表注："Extrapolation obtained by scaling the Eliashberg function at 129 GPa by a factor of 1.5. Extrapolated T_c = 353 K."
→ CaH₆ 的 ω_log(P) **非单调**（1419 @170 → 1249 @180 → 1333 @190），明显是不同来源拼接或收敛噪声，**不建议用来提取标度指数**。ThH₁₀ 的 ω_log 在 200→300 GPa 反而下降。**只有 LaH₁₀ 那一列是自洽的单一 SSCHA 计算。**

**(c) 全体系的经验 ω_log(P) 线（等级 A，Semenok Fig. S1a 说明逐字）**

> "(a) Correlation between the logarithmic average phonon frequency ω_log in Kelvin and pressure P in GPa. For about 95% of the compounds ω_log lies below the line **500 + 5P**."

即 **ω_log[K] ≲ 500 + 5·P[GPa]**（线性、非幂律）。在 170 GPa 处给 1350 K，与 LaH₁₀ @214 GPa 的 1340 K 量级一致。**这是 P3 唯一拿到的、覆盖全部氢化物的 ω(P) 经验包络，也是 Semenok 自己用来外推到 500 GPa 的依据。**
注意：它与 (a) 里 LaH₁₀ 的 n≈0.7 幂律不矛盾——线性 ω=500+5P 在 P≫100 GPa 时局部对数斜率 = 5P/(500+5P)，在 129/214/264 GPa 处分别为 0.56/0.68/0.73，与 LaH₁₀ 实测的 0.99/0.66/0.44 同量级但趋势相反（经验线的斜率随 P 上升，LaH₁₀ 的下降）。**这个矛盾没有解决，列为不确定点。**

**(d) E_F(P)：未找到一手数据。**
我在本工作区所有 PDF 与抓取到的原文里**没有找到任何氢化物的 E_F(P) 数据点或 Grüneisen 参数 γ = −dlnω/dlnV 的一手数值**。能找到的最接近的三条：
- Semenok 2407.12922 p.7 的定性陈述（等级 A）："We define a metal as a good conductor with a Fermi energy E_F much larger than typical phonon frequencies."——只给了不等式，无数值。
- 本工作区 `theory/derivations/P2_vertex_sign_criterion.md` §（P2 推导）记录：Mishra 2026 认为 H₃S 的有效费米能标由 ε_F 下方的 vHs 决定，≈40–50 meV @200 GPa，且"vHs 随压力靠近 ε_F 使 m 增大"。这是 P2 转录的 Mishra 结论，**我未回 Mishra 原文核对该数**，标 `[二手，待核]`。
- Mishra 数据集 `theory/data/external/mishra_h3s_d3s/band/H3S.dos` 里 E_F = 17.619 eV（从带底算的 QE 值，不是有效 E_F）。

**→ P3 若需要 E_F(P)，必须新做或新找：候选是 Sano et al. PRB 93, 094525 (2016)（arXiv:1512.07365，本地有 PDF，做了 H₃S 40 eV 范围的能量依赖 DOS 与压力依赖）——我本轮没读它，这是最明确的下一步。**

**(e) Grüneisen 参数 γ：未找到一手来源。**
氢化物文献里没有人直接报 γ。可替代路径（P3 需自行选）：
- 由 ω_log(P) + 状态方程 V(P) 合成 γ = −dlnω/dlnV = n·(dlnP/dlnV)⁻¹ ... 需要 V(P)。**Drozdov 2019 给了 LaH₁₀ 的一个体积点**（p.3，逐字）："face-centered cubic (fcc) lattice Fm3̄m with the refined lattice constant a = 5.1019(5) Å (V = 132.80(4) Å³)" @150 GPa；并引理论值 "131.9 Å³ at 172 GPa and 128.8 Å³ at 175 GPa"。**Troyan 2021 Table 1 给了 YH₆ 的完整 V(P)**（等级 A，本地 PDF p.3）：166/168/172/177/180 GPa → V = 45.82/45.91/45.53/45.34/45.07 Å³。由此 dlnV/dlnP = ln(45.07/45.82)/ln(180/166) = **−0.203**[本报告自算]。**这是本轮找到的唯一一组一手的氢化物 V(P)。** 但 YH₆ 没有配套的 ω_log(P) 序列，所以拼不出 γ。
- **结论：目前无法从一手数据给出任何氢化物的 Grüneisen 参数。**

---

## §5　未找到一手来源的项

1. **Allen–Dynes 1975, PRB 12, 905 原文**（无 arXiv，付费墙）。λ→∞ 渐近系数 0.1827 与 f₁/f₂ 的推导条件全部依赖二手转录（Xie 2022、Semenok 2025、Yuzbashyan 2025 三处一致）。
2. **Kiessling–Altshuler–Yuzbashyan, J. Stat. Phys. 192, 93 & 94 (2025)**（arXiv:2409.00532 及 III），Semenok 界的严格数学基础，未取。**P3 最该补的两篇。**
3. **Hofmann–Chowdhury–Kivelson–Berg, npj Quantum Materials 7, 83 (2022)**，"Heuristic bounds on superconductivity and how to exceed them"，未取。
4. **Sadovskii, arXiv:2407.03602 / Phys. Usp. 65, 724 (2022)**，Eliashberg-McMillan 框架内的 Tc 上限，未取。
5. **Trachenko–Monserrat–Hutcheon–Pickard, arXiv:2406.08129**，ω_max = E_electronic/√m = 3680 K 的推导，未取（Gao 2025 §转述）。
6. **ThH₉/ThH₁₀ 原文**（Semenok et al., Mater. Today Phys. 2020）、**(La,Ce)H₉₋₁₀ 原文**（Chen Nat.Commun.14, 2660 (2023) / Bi Nat.Commun.13, 5952 (2022)）：CSV 中这三行的 Tc/P 只有 Semenok Table S1 一处编纂来源。
7. **CaH₆ 的 207 K @181 GPa / 212 K @170 GPa 批次读数**（工作区旧 CSV），未回 Ma et al. 正文。
8. **CaH₆、ThH₁₀、CeH₉ 的 ω_log(P) 原始出处**（Semenok Table S4 未给逐行引用）。
9. **任何氢化物的 Grüneisen 参数 γ 与 E_F(P)**（见 §4.3d/e）。
10. **Semenok 2407.12922 Fig. S1a/S1b 的底层数据集**（图说明只给 "Data for the plot are from Refs. [52, 46, 48, 49, 50, 51, 47, 53]"，是编纂的编纂）。

---

## §6　关键不确定点

**U1（最大的一处）·「λ ≲ 4」这个数在 P3 里不能当常数用。**
它不是一个界，而是四个不同的量被压成了一句话：λ_*(Einstein)=3.69、λ_*(Debye)=4.72、λ_c=λ_*ξ_c（取 ξ_c=0.5 → 1.85–2.36）、以及"实测氢化物的 λ_c 落在 2.19–4.12"。其中 **ξ_c = 0.5 是 Semenok 等明说的 conjecture**，而 λ_* 由 α²F 的形状函数 P(ω) 决定 → **λ_*(P) 有压力依赖但文献没给**。P3 若要"λ 上限 + 求和规则 + Grüneisen"三合一，第一步必须是把 λ_*[P(ω)] 的泛函（Eq. S12）对压力展开——这是文献真空，也是 P3 最可能的原创点。

**U2 · ω̄ 取哪个频率决定了所有界之间是否矛盾。**
同一批氢化物：ω_log ≈ 900–1350 K、ω_D ≈ 1100–1350 K、ω_max ≈ 2700–3400 K（Semenok Table S1/S5）。Esterlis 的 0.1ω_D 在 (La,Ce)H₉ 上被违反 1.7 倍，但他自己对 H₃S 用 ω_max 时得 0.08 几乎饱和。Semenok 的 0.20/0.32 ω_max 与 Esterlis 的 0.1 ω_D 在数值上其实相容（ω_max/ω_D ≈ 2.5 → 0.1ω_D ≈ 0.04ω_max ≪ 0.20ω_max，Esterlis 更紧）。**P3 必须先固定一个频率约定并把三套界换算到同一约定下，否则"谁更紧"是无意义的比较。**

**U3 · Semenok Table S1 把 SSCHA 理论 λ/ω_log 与实验 Tc 混在同一行，且 H₃S 行的压力标签（157 GPa）与其数值来源（Errea 2015 的 200 GPa）不符。**
我通过与工作区 `isotope_curated.md` 的逐位比对确认了 H₃S 与 LaH₁₀ 两行的真实来源。**其余行（YH₆/YH₉/ThH₉/ThH₁₀/(La,Y)H₁₀/(La,Ce)H₉）我没有做同样的溯源**，它们的压力标签同样可能与数值来源不匹配。P3 用这张表前建议逐行做一次这种比对。

**U4 · LaH₁₀ 的 ω_log(P) 硬化指数（0.99→0.44，随 P 下降）与 Semenok 全局经验线 ω_log=500+5P（局部斜率 0.56→0.73，随 P 上升）趋势相反。**
两者一个是单材料 SSCHA 序列、一个是跨材料散点包络，不是同一个量（前者是沿压力路径，后者是跨材料的 P–ω_log 相关），但 P3 若要写 "ω(P) 的 Grüneisen 标度"就必须说清用哪一个。**我倾向用 LaH₁₀ 的单材料序列做 Grüneisen，用 500+5P 做包络约束，两者不可互换。**

**U5 · Tc(P) 圆顶的机制在不同材料里不同。**
LaH₁₀：上界升、实际 Tc 降，圆顶在"Tc/T_c^emp 比值"里（§4.3a）。CeH₉：作者亲述峰来自结构相变处 λ 的极大（P6₃/mmc → C2/c），不是连续竞争。H₃S：工作区 P1b 已确认峰在结构相变上方 40+ GPa 的对称相内（Cherubini 2026）。**"Tc_max(P) 包络"和"单材料 Tc(P) 圆顶"是两件事，P3 的立论必须先把它们分开。**

**U6 · CaH₆ 在 Semenok 论文内部自相矛盾**（Table S1 归 DFT / Table S3 归 experiment / p.6 说它"不是真实存在的氢化物"，而 Ma et al. 2022 已合成并测到 215 K @172 GPa）。ξ=0.58 > 猜测的 ξ_c=0.5，若 CaH₆ 确实稳定存在，则 ξ_c=0.5 的 conjecture 被证伪。**这是检验 Semenok 稳定性判据最锋利的一个数据点，值得单独追。**

---

## §7　新增/确认的 arXiv 编号

| 编号 | 论文 | 本轮状态 |
|---|---|---|
| **1806.00488** | Esterlis, Kivelson, Scalapino, "A bound on superconducting T_c's", npj Quantum Mater. 3, 59 (2018) | **新增**（工作区此前无编号）；全文已读 |
| 2409.19562 | Yuzbashyan, Altshuler, Patra, PRL 135, 026503 (2025) | 确认；补 DOI/卷期；全文已读 |
| 2407.12922 | Semenok, Altshuler, Yuzbashyan, Adv. Mater. 37, 2507013 (2025) | 确认；**"λ≲4"的真正出处**；全文+SI 已读 |
| 2502.18281 | Gao et al., "The Maximum T_c of Conventional Superconductors at Ambient Pressure" | 确认；全文已读 |
| **2012.04787** | Semenok et al., (La,Y)H₁₀ 253 K @183 GPa, Mater. Today 48, 18 (2021) | **新增**；摘要已核 |
| **2103.16282** | Ma et al., CaH₆ 215 K @172 GPa, PRL 128, 167001 (2022) | **新增**（工作区旧 CSV 记 "doi待查"）；摘要已核 |
| **2409.00532** | Kiessling, Altshuler, Yuzbashyan, "Bounds on T_c ... II: Dispersive phonons", J. Stat. Phys. 192, 94 (2025) | **新增**；未取，建议补 |
| **2406.08129** | Trachenko, Monserrat, Hutcheon, Pickard, "Upper bounds on the highest phonon frequency and superconducting temperature" | **新增**；未取 |
| **2407.03602** | Sadovskii, "Upper limit for superconducting T_c in Eliashberg-McMillan theory" | **新增**；未取 |
| 2106.05235 | Xie et al., npj Comput. Mater. 8, 14 (2022) | 已有；本轮用作 Allen–Dynes 公式的转录源 |
| 1506.08190 / 1812.01561 / 1808.07695 / 1908.01534 / 1909.10482 / 2101.01315 | H₃S、LaH₁₀(×2)、YH₆、YH₉、CeH₉ 实验 | 已有；本轮逐条回一手原句 |
| 1502.02832 / 1907.11916 | Errea SSCHA H₃S / LaH₁₀ | 已有；本轮用于揭穿 Semenok Table S1 的来源 |
