# 七项「拿不到」条目的取件调研（2026-09-18，researcher）

> **[2026-09-18 主会话复核附注]** 项 1 已由主会话**独立重取**：自行 fetch 同一 URL 得 267 KB 扫描 PDF
> （确认非付费墙、确认是 JBIG2 扫描件 ⟹ 纯文本抽取必失败，须读页图），读 p.576 页图，
> 上述 Table I 表头与「0.085 at x=1.0 … 0.10 at x=0.77」两处**逐字相符**，判决「0.085 无 ω_c 约定」成立。
> 项 2 / 5 / 6 / 7 **未逐项复核**；篇内自标「转引」「推断」的行（arXiv:2603.03492 全表、Bennemann–Garland 单位、
> Essenberger 部分读图）按仓库纪律**不得当已验结论引用**。
>
> **主会话在同一页图上另读到两条本文未提取、但对同位素主线承重的信息**（已另记，见
> `../isotope/` 与 auto-memory）：
> 1. Table I 的 `x=1.00` 行给出 Klein 自己的 **λ_tot(PdH)=0.54 / λ_tot(PdD)=0.62**（D 大于 H，
>    与 Errea SSCHA 的 0.40/0.46 **同序**）与 **Tc(calc) 7.9 / 9.6 K** ⟹ 反解 **α = −0.282**。
>    对照：Errea SSCHA −0.38、仓库亏损预算 −0.399、实验 −0.29/−0.32。
>    ⟹ 仓库预算与**当代非谐第一性原理**（Errea）相差 5%，而两者**同向**偏离实验；
>    1977 年的 Klein 反而正中实验区间，但走的是另一套机制（光学/声学 α² 比的 x 依赖），
>    **不构成对仓库链条的独立确认，也不构成反驳**——只说明残差与 SSCHA 共享，不像出在天花板+亏损那一段。
> 2. p.576 逐字：「our quantitative prediction for the inverse isotope effect **is independent of the exact
>    values of μ\* that we have used, provided only that μ\* is the same for both PdH_x and PdD_x**.」
>    ⟹ 与仓库「同位素库仑通道经 `μ̃=μ*/λ` 进场、且非可忽略」**表面相抵**（λ 随同位素变 ⟹ μ̃ 变）。
>    两者未必真矛盾（Klein 说的是其自身 Eliashberg 解对 μ\* 取值不敏感，未给量化），
>    但这是一条**该算的交叉检验**：把仓库 α_budget 的库仑分项单独关掉，看 α 动多少。**尚未做。**

> 指导原则（用户原话）：**「理论高价值的不会和一篇文章锁死。」**
> 本文对每项先写「真正需要的量」，再写取件结论。
> **可信度标注**：`亲读` = 我自己看到页图/正文原文；`转引` = 别人复述；`推断` = 我的算术或读图。
> 项 3、项 4 按派活说明跳过（算力授权 / 本地文件决定），本文不涉及。

**本轮打开的通用路线（重要）**：`https://harvest.aps.org/v2/journals/articles/10.1103/<PhysRevX.vol.page>/fulltext`
本轮 **4/4 命中**（PRB 94.014503、PRL 39.574、PRL 33.1297、PRB 12.905），全是免订阅整篇扫描 PDF。
「APS 付费墙 ⇒ 拿不到」这个前提在本环境**不成立**。（`journals.aps.org` 与 `link.aps.org/supplemental` 仍 403/302 到墙内。）

---

## 项 1 — Klein 1977 的 μ* = 0.085 配的是哪个 ω_c

**真正需要的量**：Errea 用的 `μ*=0.085` 是在哪个截断约定下定义的，以便做 Anderson–Morel 跑动换锚。

**结论：`拿到了`，但答案是「**该值根本没有配 ω_c**」——这是个否定结论，不是缺口。**

### 链条逐环（全部亲读）

1. **Errea, Calandra, Mauri, PRL 111, 177002 (2013) = arXiv:1305.7123v2**，p.4 左栏逐字：
   > "We estimate $T_c$ from the solution of the single-band Migdal-Eliashberg equations, using $\mu^* = 0.085$ **as calculated in Ref. [18]**."

   Ref. [18] = **B. M. Klein, E. N. Economou, and D. A. Papaconstantopoulos, Phys. Rev. Lett. 39, 574 (1977)**（亲读参考表 p.5）。
   ⚠ **不是 Phys. Rev. B，是 PRL**；派活单里的「PRB」是错的。
   ⚠ Errea 全文（5 页，无 arXiv 版 SM）**没有任何一处写出他自己 Matsubara 求和的截断频率**。

2. **Klein–Economou–Papaconstantopoulos, PRL 39, 574 (1977)**，标题 *Inverse Isotope Effect and the x Dependence of the
   Superconducting Transition Temperature in PdH_x and PdD_x*。
   来源：`harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.39.574/fulltext`（4 页扫描，pp.574–577，全部亲读）。

   **p.576 左栏逐字**：
   > "The Coulomb pseudopotential $\mu^*$ was estimated for both PdD$_x$ and PdH$_x$ by using the formula of
   > Bennemann and Garland.$^{31}$ With this formula, $\mu^*$ shows a fairly weak dependence on $x$, varying
   > from **0.085 at $x = 1.0$**, to 0.10 at $x = 0.77$."

   **Table I 表头逐字**：
   > "Here, $\mu^*(x) = 0.26 N(E_F; x)/[1 + N(E_F; x)]$."

   Ref. 31 = K. H. Bennemann and J. W. Garland, in *Superconductivity in d- and f-band Metals*,
   AIP Conf. Proc. No. 4, ed. D. H. Douglass (AIP, New York, 1972), p. 103.

### 判决

- **`0.085` 是一个经验 DOS 公式的输出，不是任何 ω_c 上的重整化值。** Klein 全文**没有**出现
  Morel–Anderson 对数、没有 `ω_c`、没有「取 N 倍最大声子频率」之类的陈述。
- ⇒ **「0.085 锚在哪个 ω_c」这个问题在一手文献里无解，不是我没查到。**
  任何把它按 Anderson–Morel 跑动搬到别的 ω_c 的做法，**起点是拍的**。这正是仓库 memory
  `mutilde-is-mustar-over-lambda-2026-09-18` 里那个「17% 锚点跨度」的根源，现在有了一手依据。

### 唯一能说的「隐含约定」：来自消费点，不是来自 Klein

Klein p.574 右栏逐字：
> "Our detailed calculations of $T_c$ are based on solutions to the linearized Eliashberg equations,$^{19}$ ...
> We calculated $T_c$ also using the McMillan-type$^{20}$ equation proposed by Allen and Dynes$^{21}$ with
> results which agree within 10% of the Eliashberg-equation solutions."

（Ref.19 = Bergmann–Rainer, Z. Phys. 263, 59 (1973) + Leavens, SSC 15, 1329 (1974)；Ref.20 = McMillan PR 167, 331；Ref.21 = Allen–Dynes PRB 12, 905。）

**Allen–Dynes, PRB 12, 905 (1975)**（`harvest.aps.org/.../PhysRevB.12.905/fulltext`，亲读 pp.905–907）把约定写死了：

- p.906 左栏 Eq. (4)：$\dfrac{1}{\mu^*} = \dfrac{1}{\mu} + \ln\!\left(\dfrac{\omega_{el}}{\omega_{ph}}\right)$
- p.907 左栏 Eq. (10)(11)：$\lambda(n) = 2\int_0^\infty d\omega\,\dfrac{\alpha^2F(\omega)\,\omega}{\omega^2 + (2n\pi T)^2}$，
  $\dfrac{1}{\mu^*(N)} = \dfrac{1}{\mu} + \ln\!\left(\dfrac{\omega_{el}}{\omega_N}\right)$
- p.907 左栏逐字（**这是本项最该抄进仓库的一句**）：
  > "A two-square-well analysis (such as McMillan's) provides an approximate $T_c$ equation, which depends on a
  > number $\mu^*$, which is scaled to a cutoff at some upper phonon frequency $\omega_{ph}$ as in Eq. (4).
  > **There is no known safe best choice for $\omega_{ph}$**, and we have used $\bar\omega_2$ in our numerical solutions:"
  > $\dfrac{1}{\mu^*} = \dfrac{1}{\mu} + \ln\!\left(\dfrac{\omega_{el}}{\bar\omega_2}\right)$  (12)
- p.907 Eq. (13)：$\dfrac{1}{\mu^*(N)} = \dfrac{1}{\mu^*} + \ln\!\left(\dfrac{\bar\omega_2}{\omega_N}\right)$

**ω_c 作用在哪个变量上（仓库纪律要求写清）**：
`ω_N = (2N+1)πT` 是**松原频率求和的截断**（Eq. (7)/(14) 的求和上限 N），不是 α²F 的谱截断。
而 AD **T_c 公式里那个 μ\*** 是重整化到 `ω̄₂ = ⟨ω²⟩^{1/2}` 的那一个（Eq. (12)），两者由 Eq. (13) 相联。

⇒ **最强可辩护的表述**：Klein 的 0.085 本身无约定；它被喂进 AD 公式时，AD 自己的约定是
`ω_c = ω̄₂ = ⟨ω²⟩^{1/2}`，且 AD 本人明写「没有已知的安全最佳选择」。

### 附带：Bennemann–Garland 公式的 N(E_F) 单位（**推断，未核实到原文**）

- 反解 Klein 自己的 Table I：`0.26N/(1+N) = 0.085` ⇒ **N = 0.4857**；`= 0.101` ⇒ N = 0.6352。（我的算术）
- Klein Fig. 1 纵轴是 "DENSITY OF STATES (STATES/RY–UNIT CELL)"，E_F ≈ 0.57 Ry 处读数约 6–7 states/Ry·cell
  = 0.44–0.51 states/(eV·cell)。（我的读图）
- ⇒ 与 N = 0.4857 相符 ⇒ **Klein 的 N 很可能是 states/(eV · 单胞)**。
- ⚠ **但文献在这一点上确实不自洽**：另有论文用 "0.46 states/Ry/spin ⇒ μ*≈0.082"（数值巧合、单位不同）。
  AIP Conf. Proc. 4 (1972) p.103 原文**未取到**，单位问题**标为未核实**。

---

## 项 2 — Errea 2013 的 PdH/PdD α²F(ω)

**真正需要的量**：PdH/PdD 的 λ、ω_log、⟨ω²⟩、声学带顶；理想情况是 α²F(ω) 数值曲线。

**结论：`拿到了`（λ/ω_log/Tc/α 全表 + α²F 图）；`拿不到`（机器可读的 α²F 数值文件；⟨ω²⟩ 该文根本没给）。**

来源：**arXiv:1305.7123v2**（页 1 斜纹戳 "arXiv:1305.7123v2 [cond-mat.supr-con] 16 Oct 2013"），5 页全部亲读。
= Errea, Calandra, Mauri, PRL **111**, 177002 (2013)。

### Table I（p.4 右栏，逐字转录；ω_log 单位 cm⁻¹，T_c 单位 K）

|        | 谐波 λ | 谐波 ω_log | 谐波 T_c | SSCHA λ | SSCHA ω_log | SSCHA T_c | Exp T_c [7] | Exp T_c [8] |
|--------|-------|-----------|---------|---------|-------------|-----------|-------------|-------------|
| PdH    | 1.55  | 205       | 47      | 0.40    | 405         | 5.0       | 9           | 8           |
| PdD    | 1.55  | 150       | 34      | 0.46    | 304         | 6.5       | 11          | 10          |
| PdT    | 1.55  | 125       | 30      | 0.48    | 257         | 6.9       | —           | —           |

| 同位素系数 | 谐波 | SSCHA | Exp [7] | Exp [8] |
|-----------|------|-------|---------|---------|
| α_PdH(D)  | 0.47 | −0.38 | −0.29   | −0.32   |
| α_PdH(T)  | 0.41 | −0.29 | —       | —       |

表注逐字：「$\omega_{\log}$ values are given in cm$^{-1}$ and $T_c$'s in K. ... The value of $\alpha$ between isotopes
PdA and PdB is calculated as $\alpha_{\mathrm{PdA(B)}} = -(\ln T_c(\mathrm{PdB}) - \ln T_c(\mathrm{PdA}))/(\ln M_B - \ln M_A)$.」

### 其它可用数字（全部亲读）

- 「In the harmonic approximation $\lambda$ is independent of the mass」⇒ 谐波三行 λ 恒为 1.55。
- 「At zone center the PdH optical modes are degenerate at **488 cm⁻¹**, in good agreement with inelastic
  neutron and Raman experiments [11–14] (around **450–472 cm⁻¹** at different temperatures and H concentrations).」(p.3)
- 「H-character optical modes have the largest contribution to $\lambda$, **between 75% and 79%** of the total
  depending on the isotope. The contribution of the low-energy Pd-character acoustic modes is similar for the
  three hydrides.」(p.4)
- 均方氢位移：0.55 a.u. (PdH) / 0.48 (PdD) / 0.44 (PdT) at 0 K，「does not scale according to the harmonic $M^{-0.5}$ relation」。
- λ 定义：$\lambda = 2\int_0^\infty d\omega\,\alpha^2F(\omega)/\omega$，$\omega_{\log} = \exp\!\big(\tfrac{2}{\lambda}\int_0^\infty d\omega\,\alpha^2F(\omega)\ln\omega/\omega\big)$（p.4 左栏）。
- **该文从未给出 ⟨ω²⟩ 或 ω̄₂。** 若需要，只能从 Fig. 3 的 α²F(ω) 数字化后自算。

### 图

- **Fig. 2**：PdH/PdD/PdT 的谐波 + SSCHA 声子色散（Γ–X–K–Γ–L）+ 右侧 Pd/H 投影 PDOS。
  纵轴范围：(a) PdH −100…900 cm⁻¹；(b) PdD −100…700；(c) PdT 0…500。
  含 PdD_0.63 (T=150 K, 295 K) 与 PdT_0.70 (T=80 K) 的实验点。
- **Fig. 3**：α²F(ω) 与 λ(ω)，谐波 + SSCHA，(a) PdH 0–900 cm⁻¹、(b) PdD 0–600、(c) PdT 0–500。
  **可数字化**，但 arXiv 版无数据文件。

### 机器可读数据：`拿不到`

`α²F` 的数值文件在 arXiv 版不存在；APS 的 SM 链接 `link.aps.org/supplemental/10.1103/PhysRevLett.111.177002`
**302 到 `journals.aps.org/prl/abstract/...#supplemental`（墙内，403）**，本环境取不到。
与已有 memory「PdH/PdD/PdT 完全找不到机器可读 α²F」一致，本轮**未推翻**。

**替代路线（未取全，供下一轮）**：`arXiv:2603.03492`（2026，"Enhanced superconductivity in palladium hydrides by
non-perturbative electron-phonon effects"）有 PdH/PdD 的 λ、ω_log 表（四种自能方案），并有 α²F^(2) 双声子分解：

| PdH 自能方案 | λ | ω_log (cm⁻¹) |   | PdD 自能方案 | λ | ω_log (cm⁻¹) |
|---|---|---|---|---|---|---|
| gGDg | 0.42 | 405.4 | | gGDg | 0.48 | 294.8 |
| ⟨g⟩GD⟨g⟩ | 0.27 | 353.5 | | ⟨g⟩GD⟨g⟩ | 0.33 | 272.4 |
| gGDg+g⁽²⁾GDDg⁽²⁾ | 1.10 | 497.0 | | 同 | 1.11 | 363.9 |
| ⟨g⟩GD⟨g⟩+⟨g⁽²⁾⟩GDD⟨g⁽²⁾⟩ | 0.64 | 373.9 | | 同 | 0.72 | 309.8 |

⚠ **可信度 = 转引**（经 WebFetch 摘要器读 arXiv HTML，我**没有**看页图核对）。该文用 μ*=0.13 + Allen–Dynes。
注意它第一行 `gGDg` 的 PdH λ=0.42 / ω_log=405.4 与 Errea SSCHA 的 0.40 / 405 几乎重合，是个交叉验证信号。

---

## 项 5 — J. Supercond. Nov. Magn. 22, 603 (2009) 的 Eliashberg 方程约定

**真正需要的量**：Λ^Z 与 Λ^Δ 的写法 —— 是 `λ^ph ± λ^sf` 的加减型（K1），还是单玻色子 × 形状因子型（K2）？

**目标文献识别**：G. A. Ummarino（单作者），*Iron-Based Layered Compounds: The Effect of Negative Interband
Coupling*, J. Supercond. Nov. Magn. **22**, 603–607 (2009), DOI `10.1007/s10948-009-0459-8`。

**结论：原文 `拿不到`；但 `有等价替代`，且替代源的证据强度足够做约定级判定。**

- 原文取件：Springer 登录墙；**OpenAlex 对该 DOI 给出 `is_oa: false` / `oa_status: "closed"` /
  `best_oa_location: null` / 无 `abstract_inverted_index`** ⇒ 可写进笔记的「**已确认无 OA、无摘要**」否定证据。
  arXiv 全库按标题检索（`ti:"Iron based layered compounds" OR ti:"effect of negative interband coupling"`）**零命中**。

### 等价替代（亲读，页图）：**arXiv:0904.1808v2 = Ummarino, Tortello, Daghero, Gonnelli, PRB 80, 172503 (2009)**

同一第一作者、**同一年（2009-04 提交）**、同一体系（铁基 s±）、同一形式体系。p.2 右栏 Eqs. (1)(2) 逐字：

$$\omega_n Z_i(i\omega_n) = \omega_n + \sum_j (\Gamma^N_{ij} + \Gamma^M_{ij}) N^Z_j(i\omega_n) + \pi T \sum_{m,j} \Lambda^Z_{ij}(i\omega_n, i\omega_m) N^Z_j(i\omega_m) \qquad (1)$$

$$Z_i(i\omega_n)\Delta_i(i\omega_n) = \sum_j (\Gamma^N_{ij} - \Gamma^M_{ij}) N^\Delta_j(i\omega_m) + \pi T \sum_{m,j}\big[\Lambda^\Delta_{ij}(i\omega_n,i\omega_m) - \mu^*_{ij}(\omega_c)\big]\theta(\omega_c - |\omega_m|) N^\Delta_j(i\omega_m) \qquad (2)$$

紧接其后逐字：
> "where $\Lambda^Z_{ij}(i\omega_n,i\omega_m) = \Lambda^{ph}_{ij}(i\omega_n,i\omega_m) + \Lambda^{sp}_{ij}(i\omega_n,i\omega_m)$,
> $\Lambda^\Delta_{ij}(i\omega_n,i\omega_m) = \Lambda^{ph}_{ij}(i\omega_n,i\omega_m) - \Lambda^{sp}_{ij}(i\omega_n,i\omega_m)$.
> $\theta$ is the Heaviside function and $\omega_c$ is a cut-off energy. In particular,
> $\Lambda^{ph,sp}_{ij}(i\omega_n,i\omega_m) = \int_0^{+\infty} d\Omega\, 2\alpha^2_{ij}F^{ph,sp}(\Omega)/[(\omega_n-\omega_m)^2+\Omega^2]$,
> where $ph$ means 'phonon' and $sp$ 'spin fluctuations'."

**⇒ 判定：(K1) 加减型，`Λ^Z = Λ^ph + Λ^sp`、`Λ^Δ = Λ^ph − Λ^sp`。** 这是 Ummarino 一贯的写法。

配套的定性陈述（p.2 左栏，逐字，对 (Z-dom) 普查直接有用）：
> "in the multiband Eliashberg equations (EE) the SF term in the intraband channel has **positive sign for the
> renormalization functions $Z_i$ and negative sign for the superconducting order parameters $\Delta_i$** thus
> leading to a strong reduction of $T_c$. However, if we consider negative interband contributions in the
> $\Delta_i$ equations, the final result can be an increase of the critical temperature [14]."

**ω_c 作用在哪个变量上**：Eq. (2) 里 `θ(ω_c − |ω_m|)`，即**松原频率求和的截断**。
数值取法（p.2 右栏逐字）：「The cut-off energy is $\omega_c = 12\cdot\Omega_0$ and the maximum quasi-particle
energy is $\omega_{max} = 16\cdot\Omega_0$.」，`Ω_0` 在 5–35 meV 之间，`Y_ij = 2 meV`。
谱函数形状：Lorentzian，`α²_ij F(Ω) = C_ij[L(Ω+Ω_ij, Y_ij) − L(Ω−Ω_ij, Y_ij)]`。

### 第二个独立替代（亲读）：Ummarino 的 Jülich 讲义

G. A. C. Ummarino, *Eliashberg Theory*, in E. Pavarini, E. Koch, U. Schollwöck (eds.),
*Emergent Phenomena in Correlated Matter*, Modeling and Simulation Vol. 3,
Forschungszentrum Jülich (2013), ISBN 978-3-89336-884-6, **Ch. 13, pp. 13.1–13.36**。
开放获取直链：`https://www.cond-mat.de/events/correl13/manuscripts/ummarino.pdf`

- §9 "Multi-band Eliashberg theory", p.13.21, **Eqs. (56)(57)(58)(59)**：只写**一个** `λ_ij`，`λ_ij`
  同时进 Z 与 Δ 通道（Δ 里多减 `μ*_ij(ω_c)`），并给出约束 `λ_ij/λ_ji = N_i(0)/N_j(0)` (Eq. 59)。
- §10.1, p.13.25：「intraband coupling could only be provided by phonons $ph$, and interband coupling only by
  antiferromagnetic spin fluctuations $sf$」；并设 `λ^ph_ii = λ^ph_ij = λ^sf_ii = 0`、`μ*_ij(ω_c) = 0`，
  耦合矩阵化为 Eq. (60) 的纯非对角形式。
- p.13.26 Table 2 给出四个化合物（BaFeCoAs / LaFeAsOF / BaKFeAs / SmFeAsOF）的
  `λ_tot`、`λ_12/21`、`λ_13/31`、`λ_23/32`、`Ω_0 (meV)`。**全部是唯象拟合值。**
- ⚠ **讲义的 ω_c 取法与 2009 论文不同**：讲义 p.13.26 写 `ω_c = 30Ω_0`、`ω_max = 40Ω_0`、`Ω_max = 10Ω_0`；
  2009 论文写 `ω_c = 12Ω_0`、`ω_max = 16Ω_0`。**同一作者在两处用了不同截断** —— 引用时必须注明是哪一篇。
- ⚠ 讲义的 Eq. (56)(57) 形式上**没有**显式写出 `Λ^Z / Λ^Δ` 的 ± 拆分（因为它已把 phonon 项设为 0），
  所以**要引「加减型」这个约定，应引 arXiv:0904.1808 p.2，不要引讲义**。

**残余风险（必须写明）**：JSNM 22, 603 本身我没读到。上面是**同作者 + 同年 + 同体系**的强推定，
不是对那一篇的直接核实。若这条进承重链，应在文中写成「Ummarino 体系的约定（据 PRB 80, 172503 (2009)）」，
不要写成「JSNM 22, 603 的约定」。

---

## 项 6 — Pd（及 PdH）的声学带顶 Ω_ac,max

**真正需要的量**：Pd 声子谱的最高频率（fcc Pd 一个原子/原胞 ⇒ 三支全是声学支，带顶 = 区边界最高值）；以及 PdH 里 Pd 特征声学带的带顶。

**结论：`拿到了`，两个独立来源一致。并且——**当前 ~190 cm⁻¹ 的占位值对纯 Pd 是错的（低约 23%），但对 PdH 大致是对的。**

### 来源 A（亲读页图，带数表）

Enaroseha, Osuhor, Oyibo, *Phonon Dispersion Relation of Lead (Pb) and Palladium (Pd)*, engrxiv preprint,
`https://engrxiv.org/preprint/download/3099/5661`，**p.9 Table 3.2.1**，全部单位 THz，
"EXP (a)" 行标注 a = **Miiller and Brockhouse (1971)**（即 Can. J. Phys. **49**, 704）：

| Pd | a (a.u.) | X_T | X_L | W_T | W_L | L_T | L_L |
|----|----------|-----|-----|-----|-----|-----|-----|
| EXP (Miiller–Brockhouse 1971) | 7.34 | 4.65 | **6.71** | 4.26 | 5.69 | 3.33 | **7.01** |
| LDA (this work) | 7.30 | 5.12 | 7.17 | 4.37 | 6.01 | 3.58 | 7.38 |
| GGA (this work) | 7.52 | 5.26 | 6.35 | 3.93 | 5.26 | 3.17 | 6.53 |

⚠ **首列 7.34 是晶格常数 a（原子单位），不是频率。** 一次网页检索摘要把它误读成 "X longitudinal 7.34 THz" —— 别抄。

⇒ **Ω_ac,max(Pd, 实验) = L_L = 7.01 THz**。
换算（我的算术）：`× 33.35641 cm⁻¹/THz = 233.8 cm⁻¹`；`× 4.135667 meV/THz = 28.99 meV`；`≈ 336 K`。
次高 X_L = 6.71 THz = 223.8 cm⁻¹ = 27.75 meV。

### 来源 B（亲读页图，图）

D. A. Stewart, *Ab-initio investigation of phonon dispersion and anomalies in palladium*,
**arXiv:cond-mat/0606767v2**（页 1 戳 "30 Oct 2007"）= New J. Phys. **10**, 043025 (2008)。
Fig. 1（PDF p.13）：纵轴 0–8 THz，Miiller–Brockhouse 120 K 实验点（红圆）最高约 **7.0 THz（在 L 点纵支）**，
DFPT 计算（黑三角）略高，约 7.2–7.3 THz。与来源 A 的 EXP L_L=7.01 / LDA 7.38 一致。
Stewart p.7 逐字解释了系统性偏高：「experimental studies were performed at 120 K, while first principle
calculations assume zero temperature ... Using LDA, the equilibrium lattice constant for fcc Pd was determined
to be 3.88 Å, while the measured lattice constant is slightly larger, 3.89 Å ... This trend of overestimating
the phonon frequency is greatest in the longitudinal acoustic branch.」

**两来源的实验数据同源**（都出自 Miiller–Brockhouse 1971），所以这是「同一实验、两处独立转录一致」，
不是两个独立实验。真正独立的第二条腿是 Stewart 的 DFPT 计算（7.2–7.3 THz）。

### 占位值 190 cm⁻¹ 的诊断（**推断**）

- `190 cm⁻¹ = 5.70 THz`，几乎正好等于表里的 **W_L = 5.69 THz**，也几乎正好等于 Pd 德拜温度换算
  `Θ_D(Pd) ≈ 274 K × 0.6950 cm⁻¹/K = 190 cm⁻¹`。
- ⇒ 占位值**很可能是 Θ_D 换算，不是带顶**。用作 Ω_ac,max(Pd) 时**偏低约 23%**。

### PdH / PdD 的声学带顶：`拿到了`（一手 + 逐字）

**J. M. Rowe, J. J. Rush, H. G. Smith, M. Mostoller, H. E. Flotow, *Lattice Dynamics of a Single Crystal of
PdD₀.₆₃*, Phys. Rev. Lett. 33, 1297 (1974)**
（`harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.33.1297/fulltext`，pp.1297–1299 全部亲读）。
**这正是 Klein 1977 Ref.17 用来做 Born–von Kármán 拟合的那份数据**，所以它同时给项 1 收口。

p.1299 左栏逐字：
> "The acoustic modes are **20–30% lower in frequency than the corresponding modes in pure Pd**, a result quite
> consistent with the observed lattice expansion."

摘要逐字：
> "The acoustic modes are considerably lower in frequency than the corresponding modes for pure Pd, in accord
> with the observed lattice expansion. The optic modes, in which D motions predominate, show considerable dispersion."

⇒ **Ω_ac,max(PdD₀.₆₃) ≈ (0.70–0.80) × 7.01 THz = 4.9–5.6 THz = 164–187 cm⁻¹ = 20.3–23.2 meV**（我的算术）。
Fig. 2（PDF p.3）纵轴 0–18 THz：声学支组在 [100]/[110] 区边界约 5.5–6.0 THz，虚线（纯 Pd, 300 K, Ref.12）约 6.9–7.0 THz——读图与 20–30% 的文字一致。
Fig. 3 给出由 Born–von Kármán 力常数导出的 **PdD 与 PdH 的频率分布（PDOS）**：声学带 0 至 ~6–7 THz，
然后有能隙，PdD 光学带约 8–11 THz、PdH 光学带约 11–16 THz。（读图，L2）

**交叉核对（读图，L2）**：Errea arXiv:1305.7123 Fig. 2(a) 的 PdH SSCHA 色散，低频（Pd 特征）支组带顶约 **200 cm⁻¹**，
与上面 164–187 cm⁻¹ 区间的上沿相符（Errea 是 SSCHA 非谐、化学计量 PdH；Rowe 是 PdD₀.₆₃ 实验，两者不该完全相等）。

⇒ **占位值 190 cm⁻¹ 对 PdH/PdD 的声学带顶大致是对的（落在 164–200 cm⁻¹ 区间内），对纯 Pd 是错的。**
**用之前必须先确认要的是哪一个材料。**

---

## 项 7 — Essenberger et al., PRB 94, 014503 (2016)：第一性原理多带 λ^ph_ij / λ^sf_ij

**真正需要的量**：**非唯象拟合**的、逐带分辨的自旋涨落与声子耦合矩阵（铁基）。

**结论：`拿到了`（整篇 6 页，亲读页图）。这是本轮最高价值产出。**

来源：`https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.94.014503/fulltext`（1.2 MB，pp. 014503-1 … -6）。
标题：**F. Essenberger, A. Sanna, P. Buczek, A. Ernst, L. Sandratskii, E. K. U. Gross, *Ab initio theory of
iron-based superconductors*, Phys. Rev. B 94, 014503 (2016)**，Received 16 Nov 2014; revised 13 Jun 2016; published 1 Jul 2016。

### 体系与近似（p.2，亲读）

- 计算了 **FeSe (11)、LiFeAs (111)、LaOFeAs (1111)** 的费米面与自旋磁化率；**耦合矩阵只给 FeSe**。
- 两带各向同性近似（Eq. (2)）：`+` = Γ 附近空穴带，`−` = M 附近电子带。
- 三个配对通道同时从头算：(1) 带分辨 Eliashberg 函数 `α²F_ij`；(2) RPA 屏蔽库仑矩阵元 `w_ij`（静态）；
  (3) paramagnon 介导的 SF 相互作用 `Λ^SF_ij`（引自 Ref.[46] = Essenberger et al., PRB **90**, 214504 (2014) = arXiv:1409.7968）。
- ⚠ **有一个非第一性原理的旋钮**：磁化率用 `χ_zz = χ^KS/(1 − γ f_xc χ^KS)` 按标度参数 `γ<1` 缩放（Eq. (1)），
  用来绕开 DFT 把 AFM 基态算得过稳的问题。临界值 `γ_c = 0.69 (FeSe) / 0.91 (LiFeAs) / 0.70 (LaOFeAs)`。
  第二条路线是调 Se 位置 `z_Se`。作者自陈：「the inability of standard DFT to describe the ground state of FeSe
  forced us to introduce the parameter $\gamma$ that we cannot determine from first principles」。
  ⇒ **它比唯象拟合强得多，但不是零参数。用它时必须写明这一条。**

### λ 的定义（Ref. [86]，p.6，逐字）

> "$\lambda = \frac{1}{N}\sum_{ij}\lambda_{ij}N_i$, where $N_i$ is the density of states at $\epsilon_F$ on band $i$,
> $N = N_+ + N_-$, and $\lambda_{ij}$ is defined as $\lambda_{ij} = 2\int \frac{|A_{ij}(\omega)|}{\omega}d\omega$,
> where $A_{ij}(\omega)$ is a spectral function ($\alpha^2F_{ij}(\omega)$ or $\Lambda^{SF}_{ij}(\omega)$)
> respectively for phononic or spin-fluctuation pairing). $\lambda^{\max}$ is the maximum eigenvalue of the
> $\lambda_{ij}$ matrix and whithin BCS determines $T_c$ in a multiband framework [98] acting as effective coupling."

⚠ **注意定义里有绝对值 `|A_ij|`** ⇒ `λ_ij` 全部非负，符号信息在谱函数本身，不在 λ 里。
⚠ `λ_ij ≠ λ_ji`（见下表 0.042 vs 0.048），符合 `λ_ij/λ_ji = N_i/N_j` 类的 DOS 加权。

### λ^ph_ij（Fig. 3(a)，在「theoretically optimized z_Se」处；图内标签逐字）

| | `+` (Γ 空穴带) | `−` (M 电子带) |
|---|---|---|
| `+` | λ_{++} = **0.028** | λ_{+−} = **0.042** |
| `−` | λ_{−+} = **0.048** | λ_{−−} = **0.029** |

图中 `++` 与 `−−` 曲线向上（正），`+−` 与 `−+` 曲线**向下画（负）**，频率窗 20–60 meV。
正文 p.3 右栏：「Calculations show [83] the el-ph coupling [Fig. 3(a)] to be very small ... With the integrated
coupling [86] λ being less that 0.1, $T_c$ would be exponentially small if this was the only pairing channel.」
脚注 [87] 逐字（**这条对多带天花板 (Z-dom) 判据直接有用**）：
> "The fact that phonons reduce the $T_c$ in this material is a direct consequence of the stronger interband than
> intraband coupling (**∼0.03 intraband and ∼0.05 interband**), effectively acting against the $s_\pm$ ordering."

### λ^sf（SF 通道）

- **Fig. 3(c)** 插图，随 γ/γ_c：`λ^max = 0.07`（γ/γ_c=0.9）→ `0.16` → `0.37` → `0.48`（γ→γ_c）。
  正文 p.3 右栏逐字：「Its effective coupling $\lambda^{\max}$ [86] rises from the negligible value of 0.07 at
  $\gamma/\gamma_c = 0.9$ up to 0.48 as $\gamma \to \gamma_c$.」
- **Fig. 4(a)** 图注逐字：「Spin fluctuation spectral function $\Lambda^{SF}$ computed at $z_{Se}=0.22$,
  leading to $\lambda^{\max} = 1.3$.」（纵轴 −0.2…0.3，`++` 与 `−−` 为负、`+−`/`−+` 为正，ω 到 ~400 meV）
- **Fig. 5(a)**（各向异性，BZ 分成五区：Γ₁,Γ₂,Γ₃ 空穴袋 + M₁,M₂ 电子袋），**区分辨 SF 耦合**：

  | 区对 | λ_ij | s± 通道中的作用 |
  |------|------|----------------|
  | Γ₁–Γ₂ | **0.01** | 排斥（连接同号能隙区域） |
  | M₁–M₂ | **0.07** | 排斥 |
  | Γ₃–M₁ | **0.19** | 吸引（连接 + 与 −） |
  | Γ₁–M₁ | **0.18** | 吸引 |
  | Γ₂–M₂ | **0.11** | 吸引 |

  图注逐字：「Interactions on the left act as repulsive in the $s_\pm$ channel (as they connect portions of the
  Fermi surface with the same gap phase), while those on the right act as attractive (by connecting regions with
  a + phase with those with a − phase).」
  ω 窗：左 0–800 meV，右 0–800 meV。

### 库仑通道

Fig. 3(b)：RPA 屏蔽库仑矩阵元 `w_ij(E,E')`，在 ±5 eV 的能量窗上按等能面平均，色标 0–0.6 eV。
正文：「The Coulomb pairing [Fig. 3(b)] is, as expected, diagonally dominated (++ and −− components).
Within the static approximation, this structure cannot lead to any pairing」。

### T_c 分解（p.3 §V + p.4，逐字/近逐字）

- 只算 SF：T_c 高达 **11 K**；加库仑 → **4 K**；再加声子 → **3 K**（估计上限）。实验 8 K。
- 第二条路线（z_Se = 0.22，λ^max=1.3）：三通道全开 **24 K**，只 SF **32 K**。
- 各向异性（§VI）：T_c 只升几个百分点，但能隙在 s± 内被显著调制（M₁, M₂ 处大振荡）。

### 与派活单里已知信息的对照

- Ref. [46] = Essenberger, Sanna, Linscheid, Tandetzky, Profeta, Cudazzo, Gross, **PRB 90, 214504 (2014)**
  = arXiv:1409.7968（两带电子气模型）——**确认**派活单说的「前作只有模型」。
  **2016 这篇正是把它搬到真实材料上的那一步。**
- 同组后续/姊妹工作（从参考表亲读，均未取正文，**供下一轮**）：
  - [51] A. Linscheid, A. Sanna, A. Floris, E. K. U. Gross, **PRL 115, 097002 (2015)**
  - [52] A. Linscheid, A. Sanna, F. Essenberger, E. K. U. Gross, **PRB 92, 024505 (2015)**
  - [53] A. Linscheid, A. Sanna, E. K. U. Gross, **PRB 92, 024506 (2015)**
  - [95] J. Lischner, T. Bazhirov, A. H. MacDonald, M. L. Cohen, S. G. Louie, **PRB 91, 020502 (2015)**
    （正文 p.4 说本文结果「in agreement with the work of Lischner and co-worker [95]」）
  - [78] Supplemental Material at `http://link.aps.org/supplemental/10.1103/PhysRevB.94.014503`
    （"a brief review of SCDFT and the density functional used in this work"）——**未尝试取**，
    但按项 2 的经验，`link.aps.org/supplemental` 会 302 到墙内。
  上述 4 条都可用 `harvest.aps.org/v2/journals/articles/10.1103/PhysRev<X>.<vol>.<page>/fulltext` 直接试。

### 缺口（诚实陈述）

- **只有 FeSe**（11 族硫族化合物），**没有磷化物/砷化物（LiFeAs, LaOFeAs, BaKFeAs…）的 λ_ij 矩阵**。
  若 A 类 s± 目标人群是 122/1111 铁砷，这篇**不直接覆盖**。
- **两带近似**下的 λ^ph_ij 是 2×2；五区 λ^sf 是另一套分区，**两者不在同一分区方案上**，
  不能直接拼成一个 5×5 的 (λ^ph, λ^sf) 配对矩阵。
- **γ / z_Se 是可调的**，且 T_c 在 3–24 K 之间随之变化（作者自陈"varies between 3 and 24 K"）。
  ⇒ 报「第一性原理 λ_ij」时必须连同 γ/z_Se 的取值一起报。

---

## 汇总表

| 项 | 结论 | 关键出处 | 最大的坑 |
|---|------|---------|---------|
| 1 | `拿到了`（否定结论） | Klein PRL **39**, 574 (1977) p.576 + Table I 表头；AD PRB **12**, 905 p.907 | **0.085 根本没有 ω_c**；Errea 引的是 PRL 不是 PRB |
| 2 | `拿到了`（表+图）/ `拿不到`（数值文件、⟨ω²⟩） | arXiv:1305.7123v2 Table I, Figs. 2–3 | 该文无 ⟨ω²⟩；APS SM 取不到 |
| 5 | `有等价替代` | arXiv:0904.1808v2 p.2 Eqs.(1)(2) | JSNM 原文确认 closed/no-OA；讲义与论文 ω_c 取法不同（30Ω₀ vs 12Ω₀） |
| 6 | `拿到了`（两来源） | engrxiv Table 3.2.1（Miiller–Brockhouse 1971）；cond-mat/0606767 Fig.1；PRL **33**, 1297 p.1299 | 190 cm⁻¹ 对纯 Pd 低 23%，对 PdH 大致对；表首列 7.34 是晶格常数不是频率 |
| 7 | `拿到了`（全文） | PRB **94**, 014503 (2016) Figs. 3(a), 3(c), 4(a), 5(a) + 脚注 [86][87] | 只有 FeSe；有 γ/z_Se 旋钮；两套分区不可拼 |

## 建议的后续检索方向

1. **AIP Conf. Proc. No. 4 (1972) p.103**（Bennemann–Garland 原文）——把 `N(E_F)` 的单位钉死。
   这是项 1 唯一剩下的未核实环节。老会议录，可试 OSTI 标题短语检索 / INIS。
2. **harvest.aps.org 批量扫**：[51][52][53][95] 四篇（Linscheid ×3 + Lischner），找铁**砷**化物的
   第一性原理 λ_ij。这是项 7 缺口的最直接补法，成本极低（每篇两次调用）。
3. **Errea 的 ⟨ω²⟩**：若承重，唯一路线是数字化 arXiv:1305.7123 Fig. 3 的 α²F(ω)。
   先确认 ⟨ω²⟩ 到底进不进链条（预算优先纪律）。
4. **arXiv:2603.03492 页图核对**：本轮该文全部数字都是摘要器转引，未看页图。它若要进承重链必须先补这一步。
