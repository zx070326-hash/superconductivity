# 铁基超导带分辨 λ_ij 取件：Linscheid ×3 + Lischner（2026-09-18）

> **性质**：纯取件，不做物理判定。派活单（team-lead，2026-09-18）：对四篇论文各取三件——
> (1) 带分辨 `λ_ij` 数值 + 带标签定义；(2) `λ_ij` 定义式逐字（**有没有绝对值**）；(3) 归一化约定（是否含 DOS 因子、`N_j` 取在哪）。
> 拿不到写「未找到」，不推断。
> **可信度标签**：`L1` = 本会话经 `harvest.aps.org` 落盘 PDF 后逐页读页图、逐字转录；`L1-读图` = 从图上读的近似值（非表格数字）；
> `转引` = 来自检索摘要/别篇引用、未见页图；`推断` = 我的综合。本仓四档严格性对本文不适用（无自算结论）。
> **上游**：`source_hunt_seven_items_2026-09-18.md` 项 7 的「下一轮补铁砷化物」建议——**四个目标是我自己上一轮提的**。

## 0. 一句话结果

**四篇全部取到全文（L1，harvest 4/4），但四篇都没有铁砷化物的带分辨 `λ_ij`。**
更糟的是：**Linscheid 的三篇根本不是铁基超导论文**（PRL 115 是 MgB₂/CaC₆/石墨烷实空间序参量；PRB 92 024505/024506 是磁场中的 SpinSCDFT 形式理论 + 自由电子气测试系统）。
上一轮我把它们列为「补铁砷化物最便宜路线」是**错的目标**——当时只凭作者组（Sanna/Gross 组 = Essenberger 同组）与年份推断，没有核对题目。这是我的失误，不是派活方的。
Lischner PRB 91 020502 是 FeSe（不是砷化物），给的是**对称性投影的单标量 `λ[g]`**，不是带分辨矩阵。
替代检索（§6）也**没有找到**任何第一性原理的铁砷化物带分辨 `λ_ij`；唯一取到的一手第一性原理砷化物耦合是 Boeri 2008 的**总** λ=0.21（明确不分带）。

> **这个阴性结果的估价（2026-09-18 补记，附录 A.7）**：派活方一度估价为「代价小」，理由是判据
> `λ11λ22 ≥ λ12λ21` 在任何第一性原理声子 `λ_ij` 上恒真（Gram/Bochner）。**该理由已被派活方自己撤回**
> （声子 `V_ij ∝ |g_ij|²` 是 `[g_ij]` 的 Hadamard 平方、不是 Gram 形，失效面是原子型）。
> 本笔记从未写入那条估价，故无旧句可删；当前估价是**阴性有代价**——那个测试有信息量。
> 附录 A.7 用本文已有的 FeSe 四个数做了一次算术代入，**判据在那张表上不成立**。

## 1. 总表

| 论文 | 取件 | 体系 | (1) 带分辨 λ_ij 数值 | (2) λ_ij 定义式 / 绝对值 | (3) 归一化 / DOS 因子 |
|---|---|---|---|---|---|
| Linscheid, Sanna, Floris, Gross, **PRL 115, 097002 (2015)** | L1 全 5 页 | MgB₂、CaC₆、C₂H₂（**无铁**） | **未找到**（全文无 λ） | **未找到** | **未找到** |
| Linscheid, Sanna, Essenberger, Gross, **PRB 92, 024505 (2015)** | L1 pp.1–4、16–19（pp.5–15 未读） | 形式理论 SpinSCDFT，**无任何材料** | **未找到** | **未找到**（无 λ_ij；有 Eliashberg 核 Eqs.(144)(145)，见 §3） | p.16 有 μ\*(ω_c) 约定句（见 §3） |
| Linscheid, Sanna, Gross, **PRB 92, 024506 (2015)** | L1 pp.1–3（pp.4–15 未读） | 自旋劈裂自由电子气测试系统，高斯 α²F | **未找到**（单带模型 λ=0.7，Eq.(2)） | 不适用 | 不适用 |
| Lischner, Bazhirov, MacDonald, Cohen, Louie, **PRB 91, 020502(R) (2015)** | L1 全 5 页 | **FeSe**（非砷化物） | **未找到带分辨矩阵**；只有对称性投影标量 `λ[s±]=0.30 / 0.28`（§4） | Eq.(1) 逐字（§4）：**无绝对值**，符号由 g(k) 与 V_SF 决定 | Eq.(1) 前因子 `−D_F`，`D_F` = **费米能处每自旋 DOS**，`g(k)` 归一化使 `D_F = Σ g²δ`（§4） |
| （替代）Boeri, Dolgov, Golubov, **PRL 101, 026403 (2008)** | L1 全 4 页 | LaFeAsO（砷化物，DFPT 电声） | **未找到带分辨**；总 λ=0.21（§6） | Eq.(1) 逐字：α²F 前因子 `1/N(0)`，`|g|²`（§6） | `N(0)` = 2.1 states/eV·spin（p.2 文字） |

## 2. Linscheid PRL 115, 097002 — 不是铁基论文（L1）

- 来源：`https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.115.097002/fulltext`（1.5 MB 扫描 PDF）。
- 题目逐字：「First-Principles Calculation of the Real-Space Order Parameter and Condensation Energy Density in Phonon-Mediated Superconductors」。
- 摘要逐字：「In a first application we present results for MgB₂, CaC₆ and hole-doped graphane.」
- 全文 5 页翻完：**没有 λ、没有 λ_ij、没有铁。** 唯一的方程是 SCDFT 能隙方程 Eq.(1)
  `Δ_s^{nk} = 𝒵_{nk}Δ_s^{nk} − Σ_{n'k'} 𝒦_{nk,n'k'} tanh(βE_{n'k'}/2)/(2E_{n'k'}) Δ_s^{n'k'}` 与实空间变换 Eqs.(3)(4)。
- 三件：**全部未找到**。

## 3. Linscheid PRB 92, 024505 — 形式理论，无材料（L1 部分）

- 来源：`.../PhysRevB.92.024505/fulltext`（637.7 KB，19 页）。读了 pp.1–4 与 pp.16–19（含 §V 总结、附录、参考文献）；**pp.5–15 未读**。
- 题目逐字：「Ab initio theory of superconductivity in a magnetic field. I. Spin density functional theory for superconductors and Eliashberg equations」。
- §V 总结逐字（p.17）：「In this work, we have developed fully ab initio methods to compute the SC phase of a material in a magnetic field Zeeman coupled to the spin magnetization.」— 通篇是形式推导，**没有任何材料计算**。
- (1)(2)：**未找到 λ_ij**。有的是 Eliashberg 核（p.16，L1 逐字）：
  - Eq.(144) `𝔎^σ_{n,n'}(J,J') = ∫dΩ 2Ω α²F^D_σ(0,J,0,J',Ω) / [(ω_n−ω_{n'})²+Ω²]`
  - Eq.(145) `𝔏_{n,n'}(J,J') = ∫dΩ 2Ω α²F(0,J,0,J',Ω) / [(ω_n−ω_{n'})²+Ω²] + C^stat(0,J,0,J')`
  - 这里没有出现绝对值；宗量是「能量中心 ε=0、劈裂 J」而不是带指标。
- (3) 与 μ\*/ω_c 相关的一句（p.16 右栏，L1 逐字，**顺带取到、与本仓「ω_c 作用在哪个变量」纪律直接相关**）：
  「Often, the effect of the Coulomb potential is mimicked by replacing C^stat with **μ\*θ(ω_c − |ω_n|)** where μ\* = C^stat/(1 + C^stat ln(ℰ/ω_c)) with ℰ, a parameter of the electronic band structure and **ω_c a phonon frequency cutoff** [59,60]. Usually, the so-called Morel-Anderson pseudopotential μ\* is fitted so that the calculated T_c matches the experimental one. μ\* usually ranges between 0.1 and 0.16 for conventional SC [28].」
  ⟹ 在这篇里 `ω_c` **作用在松原频率** `|ω_n|` 上（`θ(ω_c−|ω_n|)`），尽管作者称它 "a phonon frequency cutoff"。[59]=Morel–Anderson PR 125,1263；[60]=Scalapino–Schrieffer–Wilkins PR 148,263；[28]=Carbotte RMP 62,1027。
- 同页另一句（L1 逐字）：「The above equations imply that the coupling is isotropic in the sense that all states with equal center of energy and equal splitting share the same coupling matrix elements. Sometimes as in the well-known case of MgB₂ there are significant differences in the couplings and it is important to group states into bands for the isotropic approximation to hold. We refer to this case as the multiband approximation which simply means that all isotropic variables obtain another index for the band they correspond to.」— 多带只作为「加一个指标」提及，**没有给任何 λ_ij**。
- 参考文献里可能被误当线索的两条（L1 读参考文献页）：[50] Linscheid & Essenberger arXiv:1503.00970 = 「Hedin Equations for Superconductors」（arXiv abs 页核实，**与铁基无关**）；[54] A. Linscheid PhD thesis, MLU Halle-Wittenberg 2014（未取）。

## 4. Linscheid PRB 92, 024506 — 自由电子气测试系统（L1 部分）

- 来源：`.../PhysRevB.92.024506/fulltext`（3.8 MB，15 页）。读了 pp.1–3；**pp.4–15 未读**。
- 题目逐字：「Ab initio theory of superconductivity in a magnetic field. II. Numerical solution」。
- 摘要逐字：「As a test system, we employ a free-electron gas featuring an exchange splitting, a phononic pairing kernel, and a Coulomb repulsion.」
- §II 逐字（p.2）：Eq.(2) `α²F(ω) = λ (ω/2)(1/(ω_w√π)) e^{−½((ω−ω_0)/ω_w)²}`；「we fix these numbers to ω_0 = 2.2 mHa, ω_W = 0.5 mHa, and λ = 0.7 which lead to coupling properties that are loosely similar to those of MgB₂.」
- (1)(2)(3)：**不适用**——单带模型，无带分辨，无铁基。

## 5. Lischner PRB 91, 020502(R) — FeSe，对称性投影标量 λ[g]（L1 全文）

- 来源：`.../PhysRevB.91.020502/fulltext`（581 KB，5 页，全读）。arXiv 版：**arXiv:1501.02200**（转引：检索结果给出同题目同作者；页图未核）。
- 题目逐字：「First-principles theory of electron-spin fluctuation coupling and superconducting instabilities in iron selenide」。**体系是 FeSe，不是砷化物。**

### (2) 定义式（p.1 右栏 Eq.(1)，L1 逐字）

```
λ[g] = −D_F · Σ_{nk,n'k'} g(k) V_SF(nk,n'k') g(k') δ(ε_nk−ε_F) δ(ε_n'k'−ε_F)
             / [ Σ_{nk} g²(k) δ(ε_nk−ε_F) ]²                                   (1)
```

- **没有绝对值。** 符号由 `−D_F × V_SF × g g'` 决定；作者约定（p.3 逐字）「Positive coupling constants indicate that the system may undergo a transition to a superconducting phase at sufficiently low temperatures.」
- 这是**对称性泛函** `λ[g]`：给定试探间隙对称函数 g(k)，输出一个标量。**不是带分辨矩阵 λ_ij**；文中从未定义或打印任何 `λ_ij`。
- p.2 逐字：「We define V_SF(nk,n'k') as the matrix element of the spin-fluctuation-mediated interaction for the scattering of the spin-singlet Cooper pair (nk↑, n−k↓) to (n'k'↑, n'−k'↓) (with n denoting a band and k a k point in the Brillouin zone).」
- Eq.(2)（L1 逐字）：`V_SF(r,r',ω) = 3 ∫dr₁∫dr₂ I_xc(r,r₁) χ_S(r₁,r₂,ω) I_xc(r₂,r')`，`I_xc = δ²E_xc/[δm δm]`；Eq.(3) 是 TDDFT Dyson 方程 `χ_S = χ_0 + ∫∫ χ_0 I_xc χ_S`。

### (3) 归一化（p.2 左栏，L1 逐字）

- 「where ε_F denotes the Fermi energy and **D_F is the density of states per spin at the Fermi energy**.」
- 「We follow the standard convention [18] and normalize g(k) such that **D_F = Σ_{nk} g²(k) δ(ε_nk − ε_F)**.」 [18] = Maier, Graser, Scalapino, Hirschfeld, PRB 79, 224510 (2009)。
- ⟹ 前因子是**总的每自旋 DOS**（不是逐带 N_j）；分母是费米面上 g² 的平方。DOS 因子的带归属问题在这个标量泛函里不出现。
- 能带重整化（p.4 逐字）：「we use α = 3」（DFT-LDA 能量除以 3），「D'_F = αD_F」；Stoner 参数重标 「β = 4.17」；提取的 `I_xc = 0.43 eV (per Fe atom)`（p.3）。

### (1) 数值（L1）

- p.4 右栏逐字：「This yields **λ[s±] = 0.30**. Including the effect of charge fluctuations results in a slightly reduced value of **λ[s±] = 0.28**.」（重整化能带 + 重整化 Stoner 参数，实验 Se 高度）
- 同段：「T_c = ω_SF e^{−1/λ[s±]} (where we use ω_SF = 15 meV …) yields **T_c = 5 K** consistent with the experimental transition temperature of 8 K in FeSe」；「our calculation therefore provides a **lower bound** for the transition temperature.」
- p.4 逐字：「In the current Rapid Communication, we have not considered the effect of phonons」。⟹ **无 λ^ph，无 λ^ph_ij。**
- 对称函数（p.3 逐字）：s 波 `g=1`；s± `g ∝ cos(k_x)+cos(k_y)`；`d_{x²−y²} ∝ cos(k_x)−cos(k_y)`；`d_xy ∝ sin(k_x)sin(k_y)`。
- Fig. 3(a)（**L1-读图**，未重整化能带，λ[g] vs Se 高度 1.10–1.26 Å）：s 波约 −9.5 → −5（全负）；s± 约 −4.3 (1.10 Å)、−6 (1.15)、−1.2 (1.20)、−0.5 (1.23)、**+0.5 (1.26)**；d 波两条在 0 附近。Fig. 3(c)（L1-读图）：s± 的正贡献从 ≈+0.5 升到 ≈+4.3，负贡献 ≈ −4 ~ −6 大致平。
- 文字定性（p.4 逐字，可对上 s± 的「带间/带内」语言但**不是数字**）：「These spin fluctuations scatter electrons from the hole pockets at the Γ point to the electron pockets near the M point. In the s± scenario, the superconducting order parameter switches sign between the Γ and the M points, and therefore such scattering events are favorable … spin fluctuations with wave vectors in the vicinity of the Γ point scatter Cooper pairs only from electron pockets to other electron pockets or from hole pockets to other hole pockets. Their contribution to the s±-coupling constant is negative.」
- 「Note added」（p.4，L1 逐字）：「After submission of the present Rapid Communication, a similar calculation on spin-fluctuation-mediated superconductivity in iron selenide was published [40].」[40] = Essenberger, Sanna, Buczek, Ernst, Sandratskii, Gross, **arXiv:1411.2121**。
  ⟹ **订正上一轮的说法**：Essenberger PRB 94, 014503 (2016) **有 arXiv 版本 = arXiv:1411.2121 v1（2014-11-08，题「Ab-initio theory of Iron based superconductors」，同六位作者；arXiv abs 页核实）**。上一轮写「no arXiv version found」错了。⚠ v1 是 2014 年投稿版，与 2016 发表版可能有差异，逐带数字以 harvest 的发表版为准。

## 6. 替代检索：铁砷化物第一性原理带分辨 λ_ij（结果：未找到）

- 两次 WebSearch（关键词：ab initio band-resolved interband/intraband λ_ij BaFe₂As₂/LaFeAsO；first-principles spin-fluctuation coupling hole/electron pocket SCDFT paramagnon）均**没有给出任何第一性原理的砷化物带分辨 λ_ij**；检索引擎自己的结论（转引）：「none of the returned sources presents a fully ab initio, explicitly band-resolved λ_ij matrix for BaFe₂As₂ or LaFeAsO — the λ_ij matrices found are phenomenological Eliashberg fits」。
- 取了一篇作交叉：**Boeri, Dolgov, Golubov, PRL 101, 026403 (2008)**「Is LaFeAsO₁₋ₓFₓ an Electron-Phonon Superconductor?」（`.../PhysRevLett.101.026403/fulltext`，722.7 KB，4 页全读，L1）：
  - 摘要逐字：「For pure LaFeAsO, the calculated electron-phonon coupling constant **λ = 0.21** and logarithmic-averaged frequency **ω_ln = 206 K** give a maximum T_c of 0.8 K, using the standard Migdal-Eliashberg theory.」（p.3 正文写 ω_ln = 205 K，Allen–Dynes T_c=0.5 K 为上界，Eliashberg 数值解 0.8 K，μ\*=0。）
  - Eq.(1)（p.3，L1 逐字）：`α²F(ω) = (1/N(0)) Σ_{nmk} δ(ε_nk) δ(ε_mk+q) Σ_ν |g_{ν,nmk}|² …`（原式的 δ 在 ω−ω_{νq} 上；`λ(ω) = 2∫₀^ω dΩ α²F(Ω)/Ω`）。**有 |g|² 无 |α²F|；总量，不分带。**
  - N(0)（p.2 逐字）：「the DOS is high (2.1 states/eV spin)」。
  - 带分辨：**明确没有**。p.4 逐字：「In principle, multiband and/or anisotropic coupling could provide the missing factor 5 in the coupling, but this is very unlikely because it would require a very large anisotropy of the distribution [17].」只有定性分布（p.3）：「The largest couplings are concentrated around the Γ(Z) points, where the intraband nesting is large, and around the M point, where the interband nesting between the hole and electron cylinders takes place.」
- **未取、留作线索（全部未核实）**：Boeri, Calandra, Mazin, Dolgov, Mauri, PRB 82, 020506(R) (2010)（BaFe₂As₂ 磁性态电声，据记忆总 λ≈0.35，是否分带不知）；Essenberger et al. PRB 90, 214504 (2014) = arXiv:1409.7968（上一轮已判为两带电子气模型）；Linscheid, Maiti, Wang, Johnston, Hirschfeld, PRL 117, 077003 (2016)（RPA/模型，非第一性原理）。

## 7. 未核实清单

- Lischner 的 arXiv 号 1501.02200：仅检索结果标题匹配，未开页图。
- Essenberger arXiv:1411.2121 与 PRB 94, 014503 的**逐字一致性**：只核了 abs 页题目/作者/日期，未比对正文。
- PRB 92, 024505 的 pp.5–15、024506 的 pp.4–15 未读；「无 λ_ij」的判断基于题目、摘要、§V 总结与已读页。
- §6 的三条线索全部未取件。

## 8. 复核路线

四篇 + Boeri 的 PDF 均可用 `https://harvest.aps.org/v2/journals/articles/10.1103/<ID>/fulltext` 重取（本轮首次调用三篇 503，第二次全部命中；扫描件须落盘后 Read 页图）。

---

# 附录 A（2026-09-18 第三轮）：Essenberger `λ_ij` 定义式的绝对值——**结论：有，且不冗余**

> 派活（team-lead）：核 `arXiv:1411.2121` 的 `λ_ij` 定义式到底有没有绝对值，逐字抄公式 + `A_ij` 定义 + 求和/积分限 + 前因子 DOS 约定；两版若不同则**两版都抄并标版本**。三种结局中只有「没有 `|·|`」能推进判决，注意别把模棱两可读成想要的那种。

## A.1 两版对照：定义式只存在于发表版

| 版本 | 取件 | `λ_ij` 定义式 | 逐带数值 |
|---|---|---|---|
| **arXiv:1411.2121v1**（2014-11-08，`arxiv.org/pdf/1411.2121v1`，6 页全读，L1） | 有 | **全文没有任何 `λ_ij` 的定义式。** 正文 p.3 只写「The maximum eigenvalue[69] of the λ_ij matrix being λ^max=0.07」，而 **v1 的脚注 [69] 是一条纯文献引用**：「H. Suhl, B. T. Matthias, and L. R. Walker, Phys. Rev. Lett. **3**, 552 (1959)」——**不含定义**。 | Fig.3a 图标签 λ_++=0.028、λ_--=0.029、λ_+-=0.042、λ_-+=0.048（与发表版同） |
| **PRB 94, 014503 (2016)**（harvest 全文，p.014503-6，L1） | 有 | **脚注 [86] 给出定义，带绝对值**（逐字见 A.2） | 同上 |

⟹ **两版定义式「不同」的方式是：v1 根本没有定义式。** 发表版把 v1 的脚注 [69]（Suhl 引用）拆成了 [86]（定义 + 指向 [98]=同一条 Suhl 引用）。
v1 里对应发表版脚注 [87] 的是 **v1 脚注 [70]**，逐字：「The fact that phonons reduce the T_c in this material is a direct consequence of the stronger interband than the intraband coupling λ_ij (∼ 0.03 intraband and ∼ 0.05 interband), effectively acting against the s± ordering.」（与发表版 [87] 逐字相同。）

## A.2 发表版脚注 [86] 逐字（L1，PRB 94, 014503, p.014503-6 右栏）

> [86] λ = (1/N) Σ_ij λ_ij N_i, where N_i is the density of states at ε_F on band i,
> N = N_+ + N_−, and λ_ij is defined as **λ_ij = 2 ∫ |A_ij(ω)|/ω dω**, where
> A_ij(ω) is a spectral function (α²F_ij(ω) or Λ^SF_ij(ω) respectively for phononic
> or spin-fluctuation pairing). λ^max is the maximum eigenvalue of the λ_ij matrix
> and whithin BCS determines T_c in a multiband framework [98] acting as effective coupling.

（「whithin」是原文排印，照抄。[98] = Suhl–Matthias–Walker, PRL **3**, 552 (1959)。）

**逐项回答派活的四个点**：
1. **绝对值：有**，作用在被积函数 `A_ij(ω)` 上（`|A_ij(ω)|`，不是 `|∫…|`）。
2. **`A_ij` 定义**：`α²F_ij(ω)`（声子）或 `Λ^SF_ij(ω)`（自旋涨落）。**这两个都是带号量**，见 A.3。
3. **积分限**：脚注只写 `∫ … dω`，**没有写上下限**（无 `₀^∞`）。原文如此，不补。
4. **前因子 / DOS 约定**：`λ_ij` 的定义式里**没有任何 DOS 因子**——`λ_ij = 2∫|A_ij|/ω dω` 就是全部。
   DOS 只出现在把矩阵塌成标量的那一步：`λ = (1/N) Σ_ij λ_ij N_i`，`N_i` = 带 `i` 在 `ε_F` 的态密度，`N = N_+ + N_-`。
   ⚠ 注意求和里的权重是 **`N_i`（左指标）**，且 `N` 在前面做整体归一。这与「`λ_ij = N_j V_ij`」那类**把 DOS 吸进 λ_ij 自身**的约定**不是同一个约定**；是否等价取决于 `A_ij` 内部是否已含 DOS，脚注没说，**本文不判**。

## A.3 绝对值不是冗余的（关键，L1-读图）

`A_ij` 本身带号，图上直接可见：

- **Fig. 3(a)**（发表版 p.014503-3；v1 同图）**电声 α²F_ij(ω)**，纵轴 −0.2 … ≈0.5，横轴 20–60 meV：
  两个**正**峰被标注为 `λ_--=0.029`（≈34 meV）与 `λ_++=0.028`（≈47 meV）；
  两条**负**的迹被标注为 `λ_+-=0.042`（≈33–35 meV，谷深到约 −0.2）与 `λ_-+=0.048`（≈47–53 meV，负）。
  ⟹ **非对角 α²F_ij 画在零线以下，而引用的 λ_+-、λ_-+ 是正数。** 这正是 `|·|` 干的事。
- **Fig. 4(a)**（p.014503-4，z_Se=0.22 的两带 `Λ^SF(ω)`，纵轴 −0.2 … 0.3）：
  `+−`（红实线）与无标签黑虚线在正区；`++`（紫）约 −0.05、`--`（绿虚）约 −0.15，**均为负**。
- ⟹ **结局属于 team-lead 三分法的第一种：有 `|·|`，且 `A_ij` 确实带号 ⟹ 发表的 `λ_ij` 表已丢掉符号。**
  「`A_ij` 本身非负所以 `|·|` 冗余」这一支**被图否定**。

⚠ **方向自查**（按派活要求）：上面这条是「此路不通」的那一支，不是我想要的那一支，所以不存在把模棱两可往有利方向读的风险；反过来我也没有把它读得比证据更强——`|·|` 的存在是脚注逐字（L1），`A_ij` 带号是**读图**（L1-读图），后者若要承重应回到 Fig.3(a)/4(a) 的原始数据或作者数据。

## A.4 顺带取到的两条（未被要求，但同页）

- **五区分辨 SF 的 λ_ij 与符号规则**（Fig.4(a) 右上小图 + Fig.5 顶部图注，L1 逐字）：
  `Γ₁-Γ₂ λ_ij=0.01`、`M₁-M₂ 0.07`（左图）｜`Γ₃-M₁ 0.19`、`Γ₁-M₁ 0.18`、`Γ₂-M₂ 0.11`（右图）。
  图注逐字：「Interactions on the left act as repulsive in the s± channel (as they connect portions of the Fermi surface with the same gap phase), while those on the right act as attractive (by connecting regions with a + phase with those with a − phase).」
  ⟹ 在这个五区分解里，**吸引/排斥由「连接哪两块费米面 + 间隙相位」决定，不由 `λ_ij` 的符号携带**（五个数全是正的）。
- **近似方案三行**（p.014503-3，L1 逐字，标了 attractive/repulsive）：
  `α²F_ij(E,E',ω) ≈ α²F_ij(ε_F,ε_F,ω)` (attractive)；`w_ij(E,E',ω) ≈ w_ij(E,E',0)` (repulsive)；
  `Λ^SF_ij(E,E',ω) ≈ Λ^SF_ij(ε_F,ε_F,ω)` (repulsive)。

## A.5 一条独立第三方的「cut-off 名字不算证据」实例（应 team-lead 要求单列）

**Linscheid, Sanna, Essenberger, Gross, PRB 92, 024505 (2015), p.024505-16 右栏，逐字（L1）**：

> Often, the effect of the Coulomb potential is mimicked by replacing C^stat with **μ\*θ(ω_c − |ω_n|)**
> where μ\* = C^stat/(1 + C^stat ln(ℰ/ω_c)) with ℰ, a parameter of the electronic band structure
> and **ω_c a phonon frequency cutoff** [59,60].

同一句里：符号 `θ(ω_c − |ω_n|)` **作用在松原频率 `ω_n`** 上，作者却称它 "a phonon frequency cutoff"。
[59] = Morel & Anderson, PR **125**, 1263 (1962)；[60] = Scalapino, Schrieffer & Wilkins, PR **148**, 263 (1966)。
⟹ 与本仓「名字里带 cut-off 不算证据、必须由定义式确定它作用在哪个变量」是同一形状，且来自**不相干的第三方作者**。

## A.6 本附录的未核实项

- v1 与发表版的**正文**逐字一致性未做全面比对；只比对了 λ_ij 定义、脚注 [69]/[70]↔[86]/[87]、Fig.3a 的四个数值标签。
- A.3 的符号结构是**读图**，不是数表；曲线归属（哪条迹对应哪个标签）靠标注箭头位置判断。
- `A_ij` 内部是否已含 DOS 因子：脚注未说，**未查** Ref.[46]=PRB 90, 214504 (2014) 的 `Λ^SF` 原始定义。
  若 A.2 第 4 点的约定问题承重，这是下一个最短的取件目标（该文有 arXiv:1409.7968）。

## A.7 判据 `λ11λ22 ≥ λ12λ21` 在本文已取的 FeSe 声子 `λ_ij` 上的一次算术代入

> **来由**：派活方 2026-09-18 撤回了「声子 `λ_ij` 必 PSD ⟹ 判据恒真 ⟹ 砷化物阴性无所谓」的估价，
> 改判「判据是真判据、失效面原子型（对角耦合恰为 0）」。本节只做**算术代入**，不做物理判定；
> 判据本身与其物理含义都是派活方的，不是本笔记的。

用本文 A.3 / §5-上游已转录的 Essenberger PRB 94, 014503 (2016) Fig.3(a) 四个**电声** `λ_ij` 标签
（`λ_++=0.028`、`λ_--=0.029`、`λ_+-=0.042`、`λ_-+=0.048`；标签本身是 **L1**，见 A.1 两版一致）：

```
λ11·λ22 = 0.028 × 0.029 = 8.12e-4
λ12·λ21 = 0.042 × 0.048 = 2.016e-3          ⟹  比值 2.483
```

⟹ **`λ11λ22 < λ12λ21`，判据在这张第一性原理声子 `λ_ij` 表上不成立**（差 2.48 倍）。
对照派活方给的 MgB₂ 转引值（`λ_σσ=1.017, λ_ππ=0.448, λ_σπ=0.213, λ_πσ=0.155`）**通过、余量 13.8×**
⟹ 两个材料落在判据两侧，**判据在第一性原理声子 `λ_ij` 上确实会分**。

这与作者自己的物理陈述对得上（发表版脚注 [87] / v1 脚注 [70] 逐字，L1，见 A.1）：
「The fact that phonons reduce the T_c in this material is a direct consequence of the
**stronger interband than the intraband coupling** λ_ij (∼0.03 intraband and ∼0.05 interband),
effectively acting against the s± ordering.」
——「带间强于带内」正是判据失效的那个方向。

**四条限制（不可略）**：
1. **FeSe 不是砷化物**。派活单要的是铁砷化物，这条不替代 §0 的阴性结论。
2. 这四个数经过 `|·|`（A.2）。只有当 `A_ij` 在整个积分区间内符号不变时才有 `λ^abs_ij = |λ^signed_ij|`、
   从而乘积判据不受 `|·|` 影响。Fig.3(a) **读图**看是对角两条全正、非对角两条全负（**L1-读图**，非数表）
   ⟹ 这个前提**看起来**成立，但它是读图结论，若承重须回到原始数据。
3. **约定不匹配**：判据的推导设 `λ_ij = N_j V_ij`（DOS 吸在 λ 里），而脚注 [86] 明说 `λ_ij` 内**无 DOS 因子**
   （A.2 第 4 点）。在前一个约定下 `V` 对称 ⟹ `λ12/λ21 = N_2/N_1`；本表实测 `0.042/0.048 = 0.875`。
   这个比值是否真等于 `N_+/N_-`，**未核**（`A_ij` 内部是否含 DOS 未查，见 A.6）。
4. 该矩阵是 2 带（`±`），而同文的自旋涨落 `λ_ij` 是 5 区分解（A.4）且五个数全为正，
   **两套分区不能拼成一个矩阵**（§5 上游已记）。自旋涨落侧因此**没有**可直接代入的 2×2 表。

> **限制 3 已于第五轮解决（附录 B）**：约定不再是推断——`λ_ij = N_j V_ij` 由**定义式原文**确立（L1），
> 并在两张独立的 `λ_ij` 表上通过了细致平衡检验。限制 1/2/4 不变。

---

# 附录 B（2026-09-18 第五轮）：约定问题在源头解决 + 第二个失效材料（CaC₆）

> **来由**：A.7 限制 3（「脚注[86] 的 λ_ij 内无 DOS，与 `λ_ij=N_jV_ij` 不是同一约定」）被派活方判为承重。
> 派活方给的下一步是 Ref.**[46]**=arXiv:1409.7968，并随后用内部自洽性解析关掉了这个问题、通知不必取件。
> **本轮的起点是一处指针订正**：`α²F_ij` 的定义**不在 [46]**。

## B.1 指针订正：声子 `α²F_ij` 的定义在 Ref.[54]，不在 Ref.[46]（L1）

Essenberger PRB 94, 014503 **p.014503-2 右栏**逐字（本轮补读 p.2，上一轮未读）：

> Under this approximation we account for:
> (1) the phononic interaction by a band-resolved Eliashberg function α²F_ij (**detailed definitions
> can be found in Ref. [54]** and references therein);
> (2) the Coulomb repulsion by the average of the RPA screened matrix elements …;
> (3) the SF term containing an effective interaction mediated by paramagnon excitations Λ^SF
> **as introduced in Ref. [46]**, and averaged [78] as Λ_{nk,n'k'}(ω) ⟶ Λ_ij(E,E',ω) with i,j ∈ {+,−}.

参考文献页（p.014503-5/6，L1）：
- **[54] = A. Sanna, S. Pittalis, J. K. Dewhurst, M. Monni, S. Sharma, G. Ummarino, S. Massidda, E. K. U. Gross, Phys. Rev. B 85, 184514 (2012)**
- [46] = F. Essenberger, A. Sanna, A. Linscheid, F. Tandetzky, G. Profeta, P. Cudazzo, E. K. U. Gross, PRB **90**, 214504 (2014)（=arXiv:1409.7968，派活方的识别正确）
- [78] = 「See Supplemental Material at http://link.aps.org/supplemental/10.1103/PhysRevB.94.014503 for a brief review of SCDFT and the density functional used in this work.」 ⟹ **两指标平均的确切归一化在 SM 里，而 APS SM 本环境取不到**（已知 302→403）。
- [82] = P. B. Allen and B. Mitrović, *Theory of Superconducting Tc*, Solid State Physics Vol. 37 (1983), pp. 1–92。

⟹ 判掉「λ_ij 里 DOS 带哪个指标」要看的是 **[54]**（声子侧）。**该文在 harvest 上是 LaTeX 文本 PDF、8 页、非扫描件**。

## B.2 定义式原文：`λ_ij = N_j V_ij` 得证（L1，PRB 85, 184514, p.184514-2 右栏）

> α²F_{JJ'}(ω) = (1/N_J) Σ_{k∈J, k'∈J', ν} δ(ξ_k) δ(ξ_k') δ(ω − ω_{q,ν}) |g^ν_{k,k'}|²   (11)
>
> where **N_J is the density of states per spin at the Jth FS portion**.

配 Eq.(10)：`λ_{JJ'}(n,m) = ∫dω 2ω α²F_{JJ'}(ω) / [(ω_n−ω_m)² + ω²]`。

**归一化是不对称的，且方向明确**：`1/N_J` 只除**左**指标；右指标的 `k'∈J'` 求和**未归一**，
其标度为 `N_{J'}`。⟹ `α²F_{JJ'} ∝ N_{J'}·⟨|g|²⟩_{JJ'}`，`⟨|g|²⟩` 为费米面平均、**对称**
⟹ **`λ_{JJ'} = N_{J'} V_{JJ'}`，`V` 对称。这正是派活方推导所用的约定，现在是读来的、不是推的。**

⟹ 可检验的推论（细致平衡）：**`N_J λ_{JJ'} = N_{J'} λ_{J'J}`**，即 `λ_{JJ'}/λ_{J'J} = N_{J'}/N_J`。
这**不是恒真式**：若约定是 `λ_ij = N_i V_ij`，该比值会取倒数（CaC₆ 上是 0.25 vs 4.0，量级可分）。

> **⚠ 本节对判决**不承重**（派活方 2026-09-18 结果，[analytic]，见 `_scratch_ir/_chk_essenberger_lambda_ij_psd.py` 文件头 (b)）**：
> 对任意正的 `f(i), g(j)` 与对称 `V`，若 `λ_ij = f(i)g(j)V_ij`，则
> `λ11λ22 = f1f2g1g2·V11V22`、`λ12λ21 = f1f2g1g2·V12V21`，**前因子逐项相同**
> ⟹ 比值 `λ12λ21/(λ11λ22) = V12V21/(V11V22)` **精确无关于约定**。
> `f=g=1`、`λ_ij=N_jV_ij`、`λ_ij=N_iV_ij`、`λ_ij=√(N_iN_j)V_ij` 全部同判决。
> ⟹ **A.7 限制 3 其实从一开始就不影响 2.48×**（这是比「约定查清了」更强的关闭方式）。
> 本节仍有两个用处：**(i)** 它把约定从推断升为读来的；**(ii)** 它**解掉**派活方在同一文件头 (c) 里
> 标为「未解决的张力」的那条——脚注[86] 写 `λ_ij = 2∫|A_ij|/ω dω`「不含 DOS 因子」与
> 「`λ_ij` 实测不对称（0.042≠0.048）」**并不矛盾**：DOS 因子在 `A_ij = α²F_ij` **内部**
> （Eq.(11) 的 `1/N_J` + 未归一的 `k'∈J'` 求和），脚注只是没往下展开。
> 附带：派活方文件头 (a) 已自撤回「标量塌缩只在 `λ_ij=N_jV_ij` 下才说得通」这条推理；
> 本节的结论与它一致**纯属结论重合，论证不同**——这里走的是定义式。

## B.3 细致平衡在两张独立的 `λ_ij` 表上通过（算术，输入 L1）

同文 **Table I（p.184514-3）** 给出 CaC₆ 的 `λ_{JJ'}` **整表 + DOS 行**（states/eV/spin），
是本轮顺带取到的、**第二与第三个第一性原理带分辨声子 `λ_ij` 矩阵**：

三带（表注：「Isotropic DOS and λ are 0.787 and 0.870」）
```
         1       2       3
1     0.301   0.136   0.257
2     0.546   0.239   0.479
3     0.427   0.198   0.367
DOS   0.412   0.104   0.249
```
六带
```
        1a      1b      2a      2b      3a      3b
1a    0.163   0.126   0.099   0.033   0.201   0.046
1b    0.179   0.140   0.105   0.035   0.221   0.050
2a    0.331   0.245   0.151   0.084   0.384   0.096
2b    0.271   0.202   0.206   0.047   0.400   0.080
3a    0.252   0.194   0.145   0.061   0.309   0.073
3b    0.206   0.157   0.128   0.044   0.259   0.060
DOS   0.250   0.176   0.075   0.031   0.200   0.056
```

> **⚠ 转录订正（2026-09-18 第七轮，派活方从零重取同一 PDF 发现；我已自己重读 p.184514-3 页图逐行确认）**：
> 上表的六带块此前有**两处**错，均已改正——
> ① `(1b,2a)` 我印作 `0.055`，**原文就是 `0.105`**（我当时自标「疑似读错一位、应为 0.105」并隔离，猜对了、隔离方向也审对了）；
> ② `(3a,2a)` 我印作 `0.148`，**原文是 `0.145``——这处我没发现**，由派活方查出。
> 我自读该页 1b 行为 `0.179 0.140 0.105 0.035 0.221 0.050`、3a 行为 `0.252 0.194 0.145 0.061 0.309 0.073`。
> 三带表与 DOS 行两边一致。后果见 B.3 末与 B.4（六带违反 10 → **11**；`2a–3a` 1.218 → **1.193**；判决不变）。

**细致平衡检验（B.2 的推论，逐对算术；订正后的表，我逐对重算）**：
- 三带 **3/3 对通过，最劣 1.32%**（`(1,2)`：`N₁λ₁₂ = 0.412×0.136 = 0.056032` vs `N₂λ₂₁ = 0.104×0.546 = 0.056784`）；
  另两对 0.41%、1.03%。
- 六带 **15/15 对通过，最劣 1.80%**（`(1a,2b)`：`0.250×0.033 = 0.008250` vs `0.031×0.271 = 0.008401`）；
  次劣 1.63%（1b,2b）、1.61%（2b,3a），其余 ≤1.35%，多数 ≤0.7%。
- ⚠ **我上一版报的「最劣 2.8%、六带 14/15」两个数都错了**：`2.8%` 是我当时用了错的相对量（同一对我现在算得 1.32%），
  `14/15` 则是那格 `0.055` 的后果。订正后**全部 18 对通过**。

> **⚠ 方法论订正（我上一版的这句话不成立）**：我原写「这个检验**证的是** DOS 权重挂在哪个指标」。**不对。**
> 细致平衡只钉死**（排版读法 × 约定）的组合**：
> 「读法 A（行=第一指标）+ `λ_ij = N_j V_ij`」与「读法 B（列=第一指标，即**表注字面**）+ `λ_ij = N_i V_ij`」
> 是同一数组转置、约定也转置，**残差逐位相同**，该检验分不开二者。
> **真正钉死约定的只有 Eq.(11)**（`1/N_J` 挂第一指标）。**给定 Eq.(11) 之后**，细致平衡才反过来钉死
> **排版是读法 A** ⟹ 表注「The first index runs over the columns」不是让该表自洽的那个读法。
> （派活方 2026-09-18 指出；我复核同意——两条读法的残差确实同一组数。）
> **判据 `λ11λ22 ≥ λ12λ21` 对转置不变** ⟹ **B.4 的判决与这个歧义无关**，无需解决它。

**★ 一个比细致平衡更强的转录自检（派活方做的，我复核）**：按 Eq.(15) 后的标量塌缩 `λ = (1/N)Σ_ij N_i λ_ij`，
**两套分块各自独立**给 `ΣN = 0.765 / λ = 0.8685`（三带）与 `ΣN = 0.788 / λ = 0.8684`（六带），
表注印的是 `N = 0.787, λ = 0.870` ⟹ λ 差 **0.2%**、六带 DOS 差 **0.1%**。
它**同时**检验 λ 表与 DOS 行，且两套分块互为独立通道。
⚠ 但它**同样不定约定**——`(1/N)Σ_ij N_iλ_ij = Σ_i (N_i/N)λ_i` 在任何约定下都成立。**是转录检验，仅此而已。**
顺带解释了三带为何更脏：三带 `ΣN` 比印的低 2.8%，正与它那两对 3–8% 的违反同量级。

## B.4 判据在 CaC₆ 上也失效（算术；**物理判定与消费点仍归派活方**）

对每一对 `(i,j)` 代入 `λ_ii λ_jj` vs `λ_ij λ_ji`（比值 >1 即违反）：

| 分区 | 对数 | 违反 | 最劣比值 | 出处 |
|---|---|---|---|---|
| 三带 | 3 | **2** | **1.0322**（1–2 对）与 **1.0813**（2–3 对：0.087713 vs 0.094842） | (1,3) 0.9934 通过 |
| 六带 | 15 | **11** | **2.4382**（2a–2b：0.0070970 vs 0.0173040） | 次劣 (2b,3a) 1.6801、(2a,3b) 1.3563、(1a,2a) 1.3314、(2b,3b) 1.2482、(1b,2a) 1.2169、(2a,3a) 1.1933 |

（**订正**：上一版报「六带 10/15、`2a–3a` 1.218」。改正两处转录后：`(1b,2a)=0.105` 使该对从 0.6375 变 **1.2169**，
新增一个违反 ⟹ **11/15**；`(3a,2a)=0.145` 使 `2a–3a` 从 1.218 变 **1.1933**。**最劣值 2.4382 与三带的两个数不变。**
全部 18 对我已逐对重算，与派活方的独立脚本一致。）

⟹ **CaC₆ 的 `V` 不是 PSD**（任一 2×2 主子式为负即足够），**两种分区都失效**。

**这比 FeSe 那条多说了一句**：CaC₆ 是**常规电声超导体**，不是 s±、不是带间主导的非常规体系。
⟹ 判据的失效**不限于**「带间强于带内」的非常规人群。（这是我从两个数据点看到的**模式**，不是定理。）

**⛔ 必须随本节同行的限制**：
1. 这里的「带」是**费米面分块（FS portions）**——`1a/1b` 是同一条带的两块。判据用在分块分解上是否合法，
   是**物理判定，我不做**。三带划分（π 片 / Ca 球 / π 棱柱）更接近「带」。
2. **两张表不是两个独立数据点**（派活方查出，我复核）：`N_{1a}+N_{1b}=0.426` vs 三带 `N_1=0.412`（+3.4%）、
   `2a+2b=0.106` vs `0.104`（+1.9%）、`3a+3b=0.256` vs `0.249`（+2.8%）⟹ **同一费米面切两次**。
   **CaC₆ 算一个材料侧数据点，不是两个。**
3. 我**没有**钉死该判据在仓库哪条引理里被消费（派活方自陈亦未钉死）⟹ **不得传播成「天花板在 CaC₆ 上失效」**。
4. 违反幅度分布很宽（1.006–2.44），三带侧两对都在 3–8%，与表中数字自身的自洽噪声同量级
   （B.3：三带 `ΣN` 差 2.8%、细致平衡最劣 1.3–1.8%）；
   **只有六带的 1.2× 以上那几对超出噪声一个量级**。引用时不要只引 2.44×。

## B.5 顺带取到的两条（与本仓既有纪律直接相关，L1，PRB 85, 184514 p.184514-2）

- **`μ = W N` 这个约定在第三方论文里逐字出现**：「A conventional way to deal with this problem is to
  choose an energy cutoff of the order of the Fermi energy, and to assume that the product of W with the
  density of states (DOS) is constant: **μ = W N**, N being the **DOS per spin** at the Fermi energy.」
  （**订正**：我上一版漏了 "per spin" 三个字，已补；派活方查出。）
  ⟹ 与本仓 `μ_ij = N_j K_ij` 约定同形（**标量版**；矩阵版是否同样挂右指标，该句没说——但同文 p.184514-4
  的 `μ*_{JJ'} = μ*N_J/N` 半各向同性近似给了矩阵版的线索，见 B.8）。
  同段：「we choose μ*_{JJ'} in such a way to reproduce the gap structure obtained within SCDFT」
  ⟹ 该文的 **`μ*_{JJ'}` 是拟合量，不是第一性原理量**。
- **又一个 `ω_c` 作用在松原频率上的实例**：「the cutoff in the Matsubara frequency introduced for
  Coulomb terms can be uniformly applied to all the terms of the Eliashberg equations」
  （数值：Matsubara 截断 2500 meV、库仑截断 500 meV、最大声子频率 ~200 meV）。
  与 A.5 的 Linscheid 实例同形，作者与本仓无关。

## B.6 Essenberger 的带标签定义（L1，p.014503-2，补上派活单第 (1) 项）

> Δ^xc_{nk} ≈ Δ^xc_j(E) := Σ_{nk∈j} [δ(ε_nk − E)/N_j(E)] Δ^xc_{nk}   (2)
>
> … This way we can group together the **holes (Fig. 1) around the Γ point (labeled with j = +)**
> and the **electrons at the M point (labeled with j = −)**, **N_j being their density of states at the
> Fermi energy**. Within this approximation only s-wave pairing symmetries are possible, including the
> s± symmetry suggested in Ref. [4].

⟹ **`+` = Γ 点空穴口袋，`−` = M 点电子口袋**（此前我只有推断，现为逐字）。
⚠ Eq.(2) 是**单指标**量的平均（除以 `N_j(E)`）；**两指标核 `Λ_ij` 的平均只指向 SM [78]，本环境取不到**。
声子侧不受影响——B.2 的 Eq.(11) 就是声子的两指标定义式。

## B.7 本附录的未核实项（第七轮更新）

- ~~SM [78] 的两指标平均式取不到，`Λ^SF_ij` 的归一化仍未由定义式确认~~ **⟹ 已关，见附录 D**：
  该平均式的原始形式在 **Ref.[46] = arXiv:1409.7968 §VI A** 里（不在 SM 里），我已逐字取到。
- ~~CaC₆ Table I 的 `(1b,2a)` 格疑似读错一位~~ **⟹ 已关**：原文即 `0.105`；同页另有 `(3a,2a)=0.145` 一处我
  此前读错、由派活方查出（见 B.3 订正框）。两处我均已自读 p.184514-3 页图确认。
- ~~`λ_{JJ'}` 表中的值是否即 `2∫α²F_{JJ'}/ω dω`（推断）~~ **⟹ 已关，现为 L1**，见 B.8。
- PRB 85, 184514 只读了 pp.1–4；pp.5–8 未读（派活方判定 pp.4–8 为结果与讨论、对判据无增量；
  我只读了 p.4 用于核 Eq.(15) 的出处，见 B.8）。

## B.8 Eq.(15) 逐字 + 页码订正（L1，PRB 85, 184514, **p.184514-4**）

派活方报「Eq.(15) 就在 p.184514-3」。**页码差一页**：p.184514-3 上是 Eqs.(12)(13)(14) 与 §III/§IV 开头，
**Eq.(15) 在 p.184514-4 左栏**。我为核这一条读了 p.4（这是核对它引用的公式，不是去读结果与讨论）。
原文逐字，且**前一句直接把 Table I 绑到该式**：

> In Table I is presented the DOS at the Fermi level, the intraband and the interband electron-phonon couplings
>
> λ_{JJ'} = 2 ∫ α²F_{JJ'}(ω)/ω dω .   (15)

⟹ **B.7 的旧未核实项「表中的值是否即 `2∫α²F/ω dω`」现为 L1**，而且比「同文出现过该式」更强：
是**同一句话**把 Table I 指给 Eq.(15)。

**顺带（同页，L1 逐字，与本仓 `μ_ij` 约定直接相关）**：

> Interestingly, the simple *semi-isotropic* approximation μ*_{JJ'} = μ*N_J/N turns out to be sufficient to
> reproduce the SCDFT gap. … This choice of a semi-isotropic Coulombic pairing, which is often done in
> Eliashberg theory, can be validated by the present analysis. However, this can not be applied as a general
> rule; it has been shown, in cases such as in MgB₂, that a more detailed Coulombic structure is necessary
> to get the correct gap anisotropy.

⚠ **该下标在我的渲染分辨率下分不清是 `N_J` 还是 `N_{J'}`**，我不猜——照抄为 `N_J/N` 并标此不确定。
~~方向提示（**推断，不是读到的**）：与同文 Eq.(11) 的右指标约定、以及附录 D 的 Essenberger Eq.(41)
（`w_IJ ∝ N_J(E₂)`，右指标）一致的读法是 `N_{J'}`。~~
无论哪个读法，**乘积判据都不受影响**（派活方的 `f(i)g(j)` 不变性定理）。

> **第八轮订正（派活方 L1 文本层，我未能独立复核）**：派活方在同一 PDF 的**文本层**查到码点级证据——
> prime `U+2032` 出现在 `μ*` 的下标 `JJ′` 中，而 `N_J` 之后**无** prime、紧跟 `/`；并另渲染 7× 页图二次确认
> ⟹ **是 `N_J`（第一指标）**。即同一篇文章里 `λ_{JJ'}` 挂右指标、`μ*_{JJ'}` 的拟合 ansatz 挂左指标
> （不算自相矛盾：`μ*_{JJ'}` 是拟合 ansatz，不是由定义式导出的量）。
> **我上面那条「一致的读法是 `N_{J'}`」的推断方向错了**；不猜是对的，猜测方向是错的，照抄的 `N_J/N` 本身对。
> 我这边只有 Read 的页图渲染、没有文本层与高倍渲染手段，故此条**标为派活方 L1、我未独立复核**。
> 派活方顺带指出 [analytic]：该 ansatz 可分（`f(J)=N_J, g≡1, V≡μ*/N`）⟹ 按其 (S1) **恰好取等**
> `μ*_11μ*_22 = μ*_12μ*_21`，即半各向同性库仑 ansatz 正落在判据边界上。我核了这两行代数，成立。

---

# 附录 C（2026-09-18 第六轮）：铁砷化物的 `λ_ij` 矩阵——**唯象拟合侧有，且判据在它上面以原子型方式失效**

> **来由**：派活方列的两个仍缺项之一是「铁砷化物（非 FeSe）的 `λ_ij`」。§0/§6 的阴性只针对
> **第一性原理**带分辨 `λ_ij`。派活方更早提过的备选支路是「实验拟合的两带参数」——**那一支有货**，
> 而且货就在我第一轮已经取过的那篇里（`source_hunt_seven_items_2026-09-18.md` 项 3），
> 当时为回答「Λ^Z/Λ^Δ 是不是加减型」只读了 p.2 的两个方程，**没往下读到矩阵**。
> （这是 A 轮那条「取到一张表要整行抄完」教训的同形复发：**取到一篇文章要把承重对象找全**。）

## C.1 出处与体系（L1，arXiv:0904.1808v2 页图，pp.1–3）

**Ummarino, Tortello, Daghero, Gonnelli, "Three-band s± Eliashberg theory and the superconducting gaps
of iron pnictides", PRB 80, 172503 (2009) = arXiv:0904.1808v2。**
体系逐字：**Ba₀.₆K₀.₄Fe₂As₂（122，Tc=37 K）** 与 **SmFeAsO₀.₈F₀.₂（1111，Tc=52 K）**——**都是砷化物。**
Fig.1 图注逐字：「The two hole bands (1 and 2) are centered around the Γ point, while the equivalent
electron band (3) around the M point of the reduced Brillouin zone.」

⚠ **这是唯象/实验拟合，不是第一性原理**：耦合常数由拟合 Tc 与能隙定出（见 C.3）。
§0 的「第一性原理带分辨砷化物 `λ_ij` 未找到」**不因本附录改变**。

## C.2 耦合矩阵逐字（L1，p.2 右栏）

> It is obvious that a practical solution of these equations requires a drastic reduction in the number
> of free parameters of the model. On the other hand, from the work of Mazin et al. [8] we know that:
> i) λ^ph_ii >> λ^ph_ij ≈ 0 i.e. phonons mainly provide intraband coupling but the total electron-phonon
> coupling constant Σ_iλ^ph_ii should be very small [5], ii) λ^sp_ij >> λ^sp_ii ≈ 0, i.e. SF mainly provide
> interband coupling. We include these features in the most simple three-band model by **posing:
> λ^ph_ii = λ^ph_ij = 0, λ^sp_ii = 0 and μ*_ii(ω_c) = μ*_ij(ω_c) = 0.** …
>
> Under these approximations, the electron-boson coupling-constant matrix is then [8]:
> ```
>   (  0      0     λ31·ν1 )
>   (  0      0     λ32·ν2 )
>   ( λ31   λ32       0    )
> ```
> where **ν1 = N1(0)/N3(0), ν2 = N2(0)/N3(0)** and N_i(0) is the normal density of states at the Fermi
> level for the i-band.

⟹ **对角元被显式置零**：这正是派活方那条一行解析反例
`g = [[0, g₁₂],[g₁₂, 0]]`（`λ11λ22 = 0 < λ12λ21`）的**已发表实现**，且在铁砷化物上。
按本仓抽样纪律，这是**原子型**退化面（某量恰为 0）——不是连续型、不是零测搪塞得过去的。

**⚠ 但要害在于「posing」**：`λ_ii = 0` 是作者**按 Mazin 的定性序关系主动设的近似**，
不是算出来或测出来的。作者自己在同段说 `λ^ph_ii >> λ^ph_ij`（带内**大于**带间！）、
只是总电声耦合很小 ⟹ **真实的 `λ_ii` 是小而非零**。C.3 的 Table I 正是作者自己做的这个检验。

## C.3 数值与判据代入（数值 L1；代入为算术）

| 体系 | 出处 | λ31 | λ32 | ν1 | ν2 | λ_ii |
|---|---|---|---|---|---|---|
| Ba₀.₆K₀.₄Fe₂As₂（122） | p.3 左栏正文 | 4.267 | 0.569 | 1 | 2 | 0（模型设定） |
| SmFeAsO₀.₈F₀.₂（1111） | p.3 左栏正文 | 14.520 | 1.708 | 0.4 | 0.5 | 0（模型设定） |
| Ba₀.₆K₀.₄Fe₂As₂（122） | **Table I**（p.3 右下） | 3.866 | 0.471 | 1 | 2 | **0.4** |

（Table I 表注逐字：「The effect of a small contribution of intraband coupling, λ_ii = 0.4 for the case
of Ba₀.₆K₀.₄Fe₂As₂ at T = 2 K.」同表：`μ*_ij = 0`，`Δ1 = 10.30 meV`、`Δ2 = 5.62 meV`、`|Δ3| = 10.24 meV`。
ν1/ν2 的出处：122 取自 [8]，1111 取自 [18]，p.2 末–p.3 首逐字。）

判据 `λ_ii λ_jj ≥ λ_ij λ_ji` 逐对：

| 体系 / 情形 | 对 | `λ_iiλ_jj` | `λ_ijλ_ji` | 判决 |
|---|---|---|---|---|
| 两者，λ_ii=0 | (1,3)、(2,3) | **0** | >0 | **FAIL（原子型，任意大）** |
| 两者，λ_ii=0 | (1,2) | 0 | 0 | 取等（`λ12=λ21=0`） |
| 122，Table I（λ_ii=0.4） | (1,3) | 0.16 | 3.866×3.866 = **14.946** | **FAIL 93.4×** |
| 122，Table I（λ_ii=0.4） | (2,3) | 0.16 | 0.942×0.471 = **0.4437** | **FAIL 2.77×** |
| 122，Table I（λ_ii=0.4） | (1,2) | 0.16 | 0 | 通过 |

⟹ **即使把对角元从 0 放开到作者自己检验用的 0.4，判据仍以 93× 失效。**
这条比 C.2 的 `0` 更有说服力，因为它不依赖「恰好为零」这个建模选择。

## C.4 又一个不同的 DOS 约定（与派活方的不变性定理正相关）

本文的 `λ_13 = λ31·ν1 = λ31·N1/N3`、`λ_31 = λ31` ⟹ **`λ_ij/λ_ji = N_i/N_j`**，
即 DOS 挂在**左**指标（`λ_ij ∝ N_i V_ij`）——**与 B.2 的 Sanna 约定（右指标）正好相反。**

⟹ 两篇同领域论文用了两个相反的约定，**这恰好说明派活方那条不变性定理
（`λ_ij = f(i)g(j)V_ij` 下判决与 f,g 无关）在跨论文比较时是承重的**，
不是一条可有可无的观察。上表的四个判决在两种约定下逐位相同（`λ_ijλ_ji` 与对角元都不含 f,g 的差别）。

## C.5 符号在哪里：第三次「不在 λ 矩阵里」

矩阵元全为正（`λ31, λ32 > 0`）。s± 的符号反转不由 `λ_ij` 的符号携带，而由**核的 ± 结构**携带——
p.2 Eq.(2) 下逐字：`Λ^Z_ij = Λ^ph_ij + Λ^sp_ij`，**`Λ^Δ_ij = Λ^ph_ij − Λ^sp_ij`**
（与第一轮 `source_hunt_seven_items_2026-09-18.md` 项 3 记录一致）。
摘要逐字把机制称作「dominant **negative** interband-channel superconductivity (s±-wave symmetry)」。

⟹ **同一模式的第三次现身**：① Essenberger 脚注[86] 的 `|·|`（A.3）；② Essenberger Fig.5 图注把吸引/排斥
归给「连接哪两块费米面 + 间隙相位」而非 `λ_ij` 符号（A.4）；③ 本文把符号放进 `Λ^Δ` 的减号里。
**在这三处文献里，`λ_ij` 矩阵一律非负，符号信息一律在别处。**（这是我从三个数据点看到的**模式，不是定理**。）

## C.6 本附录的未核实项（第八轮更新）

- ~~只读了 pp.1–3（arXiv v2 共 4 页正文 + 参考文献）；p.4 未读。~~ **⟹ 已关，p.4 已读页图，见 C.7。**
- `ν1, ν2` 的原始出处 [8]（Mazin 等）与 [18]（**私人通信**，p.4 参考文献逐字「I.I. Mazin, private communication」）**未取**；
  本文只是引用它们的 DOS 比值。[18] 是私人通信 ⟹ **不可取**，1111 的 ν 只有此一来源。
- Table I 的 `λ_ii = 0.4` 是否三个对角元同值：表注写 `λ_ii = 0.4`（单一值），p.4 正文亦只写「we used λ_ii = 0.4」，
  **按字面理解为三个都取 0.4**，未见逐元列表。若实际只加在某一带上，(2,3) 那行的 2.77× 会变；**(1,3) 的 93× 不受影响**
  （它只需 `λ_11λ_33 ≤ 0.16`，而任何 ≤0.4 的对角元都给更小的积 ⟹ 方向保守）。
  **⟹ C.7 的救活阈把这一项变成无关**：阈的右端不含任何对角值。
- 该文是**唯象拟合**：这些 `λ_ij` 反映的是「能同时拟合 Tc 与三个能隙的参数」，
  **不是**第一性原理算出的耦合。用它说「判据在真实砷化物上失效」**要打这个折扣**。

## C.7（第八轮）：p.4 已读——Table II 补上、救活阈、计数；并对派活方「你漏了两组」一句作精确化

> **来由**：派活方从零复核 C.1–C.5（**本轮我无转录错，93.412×/2.773× 逐位复现**），并补了四点。
> 承重数字按纪律由我自读 p.3–p.4 页图 + 自算复核，不转引其数。

### C.7.1 p.4 逐字（L1，页图）——此前未读的一页

**Table II**（p.4 左栏顶）：

> `λ31 = 2.730 | λ32 = 0.758 | μ*_ij = 0.1, λ_ii = 0.4`
> `Δ1 = 7.49 meV | Δ2 = 5.72 meV | |Δ3| = 7.98 meV`
> TABLE II: The effect of the Coulomb interaction, μ*_ij for the case shown in Table I.

「for the case shown in Table I」⟹ 同一化合物（122）⟹ `ν1=1, ν2=2` 沿用。

正文三段，逐字摘承重句：

> In the case of Ba₀.₆K₀.₄Fe₂As₂ we used λ_ii = 0.4 since we know indeed that this coupling **cannot be very high** [5].
> It might be thought that this term can sensibly contribute to increase the gap values but, as can be seen in Table I,
> this is not the case: the gap values only show a slight increase (of the order of 1%).
>
> … We chose μ*_ij = 0.1 and, as expected, we found that the intraband Coulomb pseudopotential has a negligible effect
> while the interband one [14] strongly contributes to raise Tc and reduces in a considerable way the value of λ31:
> in this case, as shown in Table II, it is only possible to obtain the correct value of the small gap. As a consequence,
> **this result seems to exclude a strong interband Coulomb interaction in these compounds.**
>
> In conclusion, we have shown that the newly discovered iron pnictides very likely represent a case of **dominant
> negative interband-channel pairing superconductivity** … provided that the typical energies of the spectral functions
> are of the order of 10 meV and **the coupling constants are very high (λ31 > 4)**.

另：Fig. 4（1111，双峰谱 Ω₁=10/Ω₂=20 或 30 meV）的 `λ31, λ32` 随权重 `w_p` 变化**只有曲线、无数字**——
这才是「只在图里」的那组；p.3 正文的 F3a/F3b 两组是印出来的数字（见 C.3）。

### C.7.2 全文四组已发表 `(λ31, λ32)`，以及**免假设**的救活阈（算术，输入 L1）

派活方的改法：对 `(i,3)` 对，`λ_iiλ_33 ≥ λ_i3λ_3i = (λ3iν_i)(λ3i) = λ3i²ν_i` ⟺ **`√(λ_iiλ_33) ≥ λ3i√ν_i`**，
右端只含已发表数字 ⟹ 不必假设任何对角值，C.6 第三项的「λ_ii 是否同值」随之无关。我逐组重算：

| 集 | 化合物 | 出处（页/位置） | λ31 | λ32 | ν1 | ν2 | 阈(1,3)=λ31√ν1 | 阈(2,3)=λ32√ν2 |
|---|---|---|---|---|---|---|---|---|
| T1 | 122 | Table I, p.3（μ*=0, λ_ii=0.4） | 3.866 | 0.471 | 1 | 2 | **3.866** | **0.6661** |
| T2 | 122 | Table II, p.4（μ*=0.1, λ_ii=0.4） | 2.730 | 0.758 | 1 | 2 | **2.730** | **1.0720** |
| F3a | 122 | p.3 正文（Fig.3 参数，ω₀=10 meV） | 4.267 | 0.569 | 1 | 2 | **4.267** | **0.8047** |
| F3b | 1111 | p.3 正文（Fig.3 参数，ω₀=10 meV） | 14.520 | 1.708 | 0.4 | 0.5 | **9.1832** | **1.2078** |

八个阈**全部高于**作者自己背书的最大对角值 0.4（最小 0.6661 = 1.67×）；`(1,3)` 阈在 `ν1=1` 时**恰等于**带间 `λ31`
⟹ 救活需 `√(λ_11λ_33) ≥ λ31`，与 p.4「cannot be very high」和「dominant interband-channel pairing (λ31 > 4)」直接冲突。
**方向审计**：左端对对角单调增、右端固定 ⟹ 报阈无单边性问题（判据成立 ⟺ 真实对角几何平均 ≥ 阈）。

代入 `λ_ii = 0.4` 时的违反倍数（**只有 T1/T2 的 0.4 是发表值**；F3a/F3b 发表对角为 0，那里的 0.4 是保守替代）：

| 集 | (1,3) = λ31²ν1/0.16 | (2,3) = λ32²ν2/0.16 |
|---|---|---|
| T1 | 14.946/0.16 = **93.41×** | 0.4437/0.16 = **2.773×** |
| T2 | 7.4529/0.16 = **46.58×** | 1.1491/0.16 = **7.182×** |
| F3a | 18.207/0.16 = **113.8×** | 0.6475/0.16 = **4.047×** |
| F3b | 84.33/0.16 = **527.1×** | 1.4586/0.16 = **9.116×** |

与派活方脚本 `_scratch_ir/_chk_ummarino_arsenide_lambda_ij.py` 的输出逐位一致（我手算 4 位有效）。

### C.7.3 对派活方「你漏了两组已发表参数」的精确化

派活方写「你漏了两组已发表参数 … 这同时关掉了『1111 的 λ31/λ32 只在 Fig.2 曲线里』的旧判断」。**分两半说**：
- **C.3 第六轮就已列出 F3a（4.267/0.569）与 F3b（14.520/1.708）**，出处「p.3 左栏正文」——它们没有漏。
  「只在 Fig.2 曲线里」这句出自派活方**自己早先脚本**的注释（`1111's own lam31/lam32 are given only in Fig. 2`），不出自本笔记。
- **我真正漏的是 Table II（T2）**——它在 p.4，而 C.6 明写 p.4 未读；以及**没有在 F3a/F3b 上以 `λ_ii=0.4` 代入判据**
  （C.3 的判据表只对它们算了 `λ_ii=0` 那行）。这两项 C.7.1/C.7.2 已补。
⟹ 转录侧无缺；缺的是**一页没读 + 两组没算**。

### C.7.4 计数与折扣（采纳派活方，与 C.6 叠加）

- T1/T2/F3a 是**同一化合物**（122）三套假设下的拟合 ⟹ **一个**数据点、三条稳健性行；F3b（1111）是第二个。
  连同 FeSe（附录 A）与 CaC₆（附录 B，两张表算一个）⟹ **材料侧共 4 个**：FeSe / CaC₆ / 122 / 1111，
  其中**只有 FeSe 与 CaC₆ 来自第一性原理**；§0 的「第一性原理带分辨砷化物 `λ_ij` 未找到」不变。
- T2 是作者**自己否掉**的情形（C.7.1 第二段）⟹ 只当 μ* 方向的稳健性行。
- 零对角来自**整条声子通道被丢掉**（C.2：作者自陈 `λ^ph_ii >> λ^ph_ij`，声子是带内主导）⟹ 这里的带间主导只属**保留下来的 SF 通道**
  ⟹ 承重量必须是 C.7.2 的救活阈，不是 `λ_ii=0` 那行。
- ~~**消费点仍未钉死**~~ ⟹ **已由 `../ceiling/lambda_ij_criterion_mechanism_and_CaC6_2026-09-18.md` §10.1 钉死（2026-09-18 第十轮）**：判据只在库仑侧引理 μ-1 被消费，铁砷化物是**声子侧** ⟹ 仍不得传播成「天花板在铁砷化物上失效」。

### C.7.5 跨页句子（记录，我第六轮读对了）

p.2 末–p.3 首逐字：「In the **122** case (Tc = 37 K) we know that ν1 = 1 and ν2 = 2 [8] while in the **1111** case (Tc = 52 K)
we have ν1 = 0.4 and ν2 = 0.5 [18].」——C.3 第六轮已按此填表。派活方复核时只读到 p.3 那半句、险些把 `ν1=1,ν2=2` 判为转录错
（会把 (1,3) 93.4→37.4×、(2,3) 翻成 PASS），落笔前读了 p.2 末尾拦下。**可复用教训：跨页断开的句子，后半截单独读会给出自洽但错误的语义**——
这与我 B 轮「取到一张表要整行抄完」是同一形状（整句抄完）。

---

# 附录 D（2026-09-18 第七轮）：自旋涨落侧 `Λ^SF_ij` 的归一化——**在 [46] 正文里，不在取不到的 SM 里**

> **来由**：派活方列的优先项 ① 是「SF 侧 `Λ^SF_ij` 的归一化」。此前的阻塞点是：Essenberger PRB 94, 014503
> 把两指标平均外包给 **SM [78]**，而 APS Supplemental Material 在本环境取不到（302→403）。
> 派活方的建议是「试 arXiv:1409.7968（=PRB 90,214504）的正文，`Λ^SF` 的原始定义在那里，
> 两指标平均式可能在其附录」。
> **结果：建议成立，而且比预期好——平均式不在附录，就在正文 §VI A，且与声子侧 Eq.(11) 是同一形状。**

## D.1 出处（L1，`arxiv.org/pdf/1409.7968`，**13 页** LaTeX 文本 PDF（p.13 仅为参考文献溢出；第八轮订正，初记 12 页），正文 pp.1–12 全读）

Essenberger, Sanna, Linscheid, Tandetzky, Profeta, Cudazzo, Gross,
**"Superconducting pairing mediated by spin-fluctuations from first principles"**, arXiv:1409.7968v1（2014-09-28）
= PRB 90, 214504 (2014) = Essenberger PRB 94, 014503 的 **Ref.[46]**。
章节：I 引言 / II SCDFT 简述 / III 自能扩展（A 图集 B 粒子-空穴传播子 C 局域近似）/ IV 自能终形 /
V 泛函（Sham–Schlüter）/ VI **两带模型应用**（A 各向同性近似与两带模型 B SF 贡献讨论 C 库仑-SF-声子相互作用）/ VII 总结。

`Λ^SF` 的**微观**定义在 p.7 **Eq.(32)**（L1 逐字）：

> Λ^SF(x₁,x₂) := 2 ∬ dx dx′ f^xc_zz(x₁,x) χ_zz(x,x′) f^xc_zz(x′,x₂)   (32)

（`x = {r,τ}`；`f^xc_zz` 是 TDDFT 的 xc 核，`χ_zz` 是磁响应函数。）**这是两点函数，尚无带指标。**

## D.2 ★ 两指标带平均式逐字（L1，p.9 左栏，§VI A）——与声子侧同一形状

> Δ_nk ≈ Δ_n(E) := (1/N_n(E)) Σ_k δ(ε_nk − E) Δ_nk
> N_n(E) := Σ_k δ(ε_nk − E) .
>
> where N_n(E) is the density of states of band n. This simplification leads to an isotropic gap equation,
> where all interactions are replaced by energy and band dependent quantities. **As an example of how this
> averaging works, we consider the spin-fluctuation term:**
>
> Λ^SF_{nn'}(E,E′,ω) ≈ (1/N_n(E)) Σ_{kk′} Λ^SF_{nk n′k′}(ω) × δ(ε_nk − E) δ(ε_n′k′ − E′).

⟹ **`1/N_n(E)` 只除左指标；右指标的 `k′` 求和未归一，其标度为 `N_{n′}(E′)`。**
这与声子侧 Sanna Eq.(11)（B.2）**逐字同构** ⟹ 若 `Λ^SF_{nk,n′k′}` 在 `(nk)↔(n′k′)` 下对称，则

**`Λ^SF_ij = N_j V_ij`，`V` 对称——DOS 挂右指标，与声子侧同一约定。**

带标签也在同页逐字给出：「The portions (bands) of the Fermi surface nested by **q_c** will be labeled as
**n = + and n = −**.」——与 Essenberger 2016 的 `i,j ∈ {+,−}` 完全一致。

## D.3 模型形式独立佐证同一约定，且**对角元被显式置零**（L1，p.9 右栏）

> In such a situation the isotropic effective interaction is expected to be small for intra-band scattering
> (**Λ^SF_±± ≈ 0**) and peaked for inter-band scattering (Λ^SF_±∓).

> Λ^SF_{IJ}(E,E′,ω) = { **c₁ N_J(E′)** [1 − ((ω − ω̄ − c₂/2)/c₂)²]   if |ω − ω̄| ≤ c₂/2 **and I ≠ J** ;
>                        **0 elsewhere.**   (38)
>
> We fix the width c₂ to a value of 0.01 Ryd and the density of states times the peak height c₁ is determined
> by requiring a value for the effective coupling strength λ^SF. The effective coupling strength is given by
> the integral of Λ^SF_{IJ} with respect to ω:
>
> λ^SF_{IJ} := (3/2π) ∫₀^∞ dω  2Λ^SF_{IJ}(ε_F,ε_F,ω)/ω .

三件事：
1. **`c₁ N_J(E′)`：DOS 带的是右指标 `J`、且取在右能量 `E′`** ⟹ 与 D.2 从定义式推出的约定**独立吻合**。
2. **`I ≠ J` / `0 elsewhere` ⟹ 该模型的 `Λ^SF` 对角元恰为 0**——这是派活方那条一行反例
   （`λ11λ22 = 0 < λ12λ21`）在 **SF 侧**的又一个已发表实现，且理由是物理的（`q_c` 只嵌套 `+`↔`−`），
   不是为反例而构造。按本仓抽样纪律属**原子型**退化。
   ⚠ 但这是**两带电子气模型**，不是材料计算 ⟹ 与附录 C 的唯象折扣同级，**不算材料侧数据点**。
3. `λ^SF_IJ` 的定义带前因子 **`3/(2π)`**，见 D.5。

⚠ **一处我不修的转录**：Eq.(38) 抛物线的宗量按字面是 `(ω − ω̄ − c₂/2)/c₂`，与条件 `|ω − ω̄| ≤ c₂/2`
配起来不是关于 `ω̄` 对称的抛物线。**可能是我在该分辨率下把小分数读错**，也可能是原文如此。
照抄不改，并标此存疑——**它不承重**（承重的是 `c₁N_J(E′)` 与 `I≠J` 两点）。

## D.4 库仑侧同一约定（L1，p.10 右栏 Eq.(41)）——第三方对本仓 `μ_ij = N_j K_ij` 的逐字佐证

> w_{IJ}(E₁,E₂) = { **N_J(E₂)** (U₀ + U₁ e^{−κ(E₁²+E₂²)})   if I = J ;
>                   **N_J(E₂)** U₀                          if I ≠ J.   (41)

⟹ 库仑矩阵元也把 **DOS 挂在右指标、取在右能量**。
⟹ 本仓 `μ_ij = N_j K_ij` 的约定在一篇不相干第三方论文里有**逐字**同形（B.5 的 `μ = W N` 是标量版，
这里是**矩阵版**）。

## D.5 一处前因子不一致（记录，但**判决中性**）

- 本文 p.9：`λ^SF_IJ := (3/2π) ∫₀^∞ dω 2Λ^SF_IJ(ε_F,ε_F,ω)/ω`
- Essenberger PRB 94, 014503 脚注 [86]（A.2）：`λ_ij = 2∫|A_ij(ω)|/ω dω`，`A_ij ∈ {α²F_ij, Λ^SF_ij}`

二者差一个全局常数 `3/(2π)`（`3/2` 是自旋因子，见 Eq.(36) 的 `Im[w + (3/2)Λ^SF]`；`1/π` 来自谱表示）。
**对乘积判据无影响**：全局常数 `c` 使 `λ11λ22` 与 `λ12λ21` 同乘 `c²`，比值精确不变。
（这是派活方 `f(i)g(j)` 不变性定理的 `f=g=√c` 特例。）**记录在案，不据此改任何数。**

## D.6 顺带：符号「不在 λ 矩阵里」的第四次现身（L1，p.9 左栏 Eq.(37)）

> ℑm[Λ^SF_{kk′}(ω)] > 0 and ℑm[w_{kk′}(ω)] > 0.   (37)
>
> If these are the strongest terms in the gap equation 36, in order to have a non-trivial solution,
> **a sign change must occur in the gap function.**

⟹ SF 与库仑**都是排斥的**（谱函数定号为正），s± 的符号反转**由能隙函数承担**，不由 `λ_ij` 的符号承担。
与 A.3（`|·|`）、A.4（Fig.5 图注把吸引/排斥归给费米面块+相位）、C.5（`Λ^Δ = Λ^ph − Λ^sp` 的减号）
构成同一模式的第四例。（**模式，不是定理。**）

## D.7 本附录的未核实项

- ⚠ **这不是 SM [78] 本身。** Essenberger 2016 的原句是「Λ^SF **as introduced in Ref. [46]**, and
  **averaged [78]** as Λ_{nk,n′k′}(ω) ⟶ Λ_ij(E,E′,ω)」——即它把**平均**这一步引给 SM，把 `Λ^SF` 的**引入**
  引给 [46]。本附录取到的是 **[46] 正文里的平均式**。两者是同一式的证据是**记号与标签**层面的
  （同样的 `Λ_ij(E,E′,ω)` 三宗量形式、同样的 `i,j ∈ {+,−}` 带标签、同一作者组），**不是逐字比对**——
  SM 仍取不到。⟹ 归一化问题**实质上关闭**，但若要 100% 逐字闭合，缺口是 SM 的那一页。
- arXiv v1（2014-09-28）与 PRB 90, 214504 发表版的逐字一致性**未比对**（只取了 arXiv 版）。
- `χ_zz`、`f^xc_zz` 的具体近似（ALDA 等）与 `Λ^SF_{nk,n′k′}` 是否严格对称：**未逐字核**。
  ~~D.2 的「`V` 对称」依赖这个对称性（声子侧由 `|g_{kk′}|² = |g_{k′k}|²` 保证，SF 侧我**未见**同样直白的一句）。~~
  **⟹ 第八轮：此缺口不阻塞判据，已关**（派活方 §8.2，我核了代数）：设 `λ_ij = f(i)g(j)V_ij`，
  `λ_ijλ_ji = (f_if_jg_ig_j)V_ijV_ji`、`λ_iiλ_jj = (f_if_jg_ig_j)V_iiV_jj`，前因子相消**全程没用到 `V` 对称**
  ⟹ 乘积判据的比值对任意 `V` 都与 `f,g` 无关。[analytic] `V` 对称只在把判据改述成「`V` 的 2×2 主子式 ≥ 0 / PSD」时才需要；
  SF 侧 `Λ^SF_{nk,n′k′}` 的对称性本身仍未逐字核，只是不再承重。
- 本文 §VI 全部是**两带自由电子气模型**，无材料；Fig.1 的 `λ^SF` 是**输入参数**不是算出的耦合。
- 页数：我第七轮记「12 页」，**错**——派活方记 13 页，我读了第 13 页确认：是参考文献 [58]–[63] 的溢出页（同一 953.8 KB 文件）。
  正文与 §VI 均在 pp.1–12，无实质影响；D.1 已同步改。

---

# 附录 E（2026-09-18 第九轮）：文献里的 `λ_ij` 定义是否总落在 `f(i)g(j)V_ij` 形状里——**声子侧阴性；库仑侧找到一个非秩一映射**

> **来由**：派活方第九轮派活：「拟合侧 `λ_ij` 矩阵的约定是否总能写成 `f(i)g(j)V_ij` 的反例（有没有论文的 `λ_ij` 定义不落在这个形状里），
> 那会直接威胁约定无关定理」。本附录只取件：定义式逐字 + 算术；**是否「威胁」是派活方的判定**。
> 4 次检索 + 5 篇一手（其中 4 篇本轮读页图）。

## E.0 一句话结果

- **声子 `λ_ij` 定义：未找到任何非乘积形状。** 见到的全部是 `h(N_i,N_j) ∈ {N_j, N_i, √(N_iN_j), N_j/N_F, 1}`，每一种都是 `f(i)g(j)`。
  且有**两篇第三方把这个结构当一般事实写出来**（E.1 Mitrović：「`α²F_ij = α²f_ij N_j`, `α²f_ij = α²f_ji`」；E.2 Dolgov–Mazin–Parker–Golubov：「`B_12/B_21 = N_2/N_1`」）。
- **找到的唯一非秩一映射在库仑侧**：Morel–Anderson 截断重整化的**矩阵版** `μ̂*(ω_c) = (1̂ + μ̂ ln(E_F/ω_c))⁻¹ μ̂`（E.5，Mitrović Eq.(12)–(16) 逐字）。
  它不是 `f(i)g(j)×`：2×2 时 `det μ̂* = det μ̂ / D`（判决符号在 `D>0` 时保持），但比值 `μ*_12μ*_21/(μ*_11μ*_22)` 随 `ω_c` 变（算术，E.5）。
- **顺带关掉本仓一条「原文未读」**：MgB₂ 的 `λ_σσ=1.017, λ_ππ=0.448, λ_σπ=0.213, λ_πσ=0.155` 与 `N_σ=0.300, N_π=0.410` 在 Golubov JPCM 14,1353 = cond-mat/0111262v2 **p.4 逐字**（E.3）；
  同页另给 `μ*` 矩阵四个值。Mitrović 转引同一组 α²F 时印 `0.446/0.212`（第 3 位不同，E.1）。

## E.1 Mitrović, cond-mat/0402473v1（L1 页图，pp.1–3；「EPJ manuscript」，期刊卷页未核）

p.2 左栏 §2 逐字：

> In the case of several bands i = 1, 2, … with different partial densities of states N_i the Eliashberg equations for T_c do not
> have the form of a Hermitian eigenvalue problem [11]. That is because the interband electron-phonon coupling functions α²F_ij(Ω)
> and the corresponding Coulomb repulsion parameters μ_ij **are proportional to N_j** and are not symmetric under the exchange of the
> band indices i and j, i ≠ j. … one first takes the cutoff ω_c in the Matsubara frequency sums to be large enough so that
> μ*_ij(ω_c) = μ_ij ≡ V^c_ij N_j, where V^c_ij is the Fermi surface averaged screened Coulomb matrix element between the states in
> the bands i and j; clearly V^c_ji = V^c_ij. **The electron-phonon coupling functions can be written as α²F_ij(Ω) = α²f_ij(Ω)N_j
> with α²f_ij(Ω) = α²f_ji(Ω).**

⟹ 第三方把「对称核 × 右指标 DOS」当一般结构陈述（对声子与库仑两侧同时）。
p.2 右栏 Eq.(7)/(9)/(10)：`φ̄_i(n) = φ_i(n)√N_i/√(|ω_n|Z_i(n))`；`λ^s_ij(n−m) = √(N_i/N_j)∫α²F_ij(Ω)·2Ω/(Ω²+(ω_n−ω_m)²)dΩ`；`μ^s_ij = √(N_i/N_j) μ_ij`
⟹ 对称化版 `= √(N_iN_j)·(symmetric)`，仍是 `f(i)g(j)` 形（`f=g=√N`）。

p.3 右栏 §3 逐字：「We used the α²F's given in [11] with the coupling parameters (see Eq. (3)) **λ_σσ = 1.017, λ_ππ = 0.446, λ_σπ = 0.212**
and λ_πσ = 0.155. … μ_σσ : μ_σπ : μ_πσ : μ_ππ = 1.75 : 2.04 : 1.61 : 1.00, consistent with the density of states ratio N_π/N_σ with the
ratio λ_σπ/λ_πσ for the spectra that we used」。
⚠ 与 Golubov 原文（E.3）的 `0.448 / 0.213` 在第 3 位不同——Mitrović 说的是「用 [11] 的 α²F 自己算的 λ」，差异归属**未定**（可能是各自积分）；
本仓转引的 `0.448/0.213` 与原文 Golubov 一致，**不是转录错**。

## E.2 Dolgov, Mazin, Parker, Golubov, arXiv:0810.1476v1 = PRB 79, 060502(R) (2009)（L1 页图 pp.1–2；4 页）

⚠ 页图渲染丢失了希腊字母（`λ`、`Δ`、`ω` 显示为方块），以下 `λ` 为按上下文补入，其余逐字：

p.1 右栏 Eq.(2) 后：「where **λ_12 = N_2 V**, the dimensionless coupling constant, with a similar expression for λ_21. These equations
readily yield λ_eff = √(λ_12 λ_21) and Δ_2/Δ_1 = √(N_1/N_2).」
p.2 右栏 Eq.(8) 后：「This B_12 represents the electron-boson coupling function which supplants the pairing potential used in BCS theory,
and there is an exactly analogous equation for band 2. **Here B_12(ω)/B_21(ω) = N_2/N_1 = β.**」

⟹ 右指标 DOS，`f(i)g(j)` 形；且 `λ_eff = √(λ_12λ_21)` 直接用了乘积（约定无关的那个组合）。
这篇是 Ummarino 2009 的 Ref.[19]（附录 C）。

## E.3 Golubov et al., cond-mat/0111262v2 = J. Phys.: Condens. Matter 14, 1353 (2002)（L1 页图 pp.1–6）——**关掉本仓「原文未读」**

p.4 §2.2 逐字：「Their respective densities of states at the Fermi energy have values of **N_σ(0) = 0.300 states/cell·eV and
N_π(0) = 0.410 states/cell·eV**. … The coupling constants corresponding to the superconducting Eliashberg functions have been
calculated to be: **λ_σσ = 1.017, λ_ππ = 0.448, λ_σπ = 0.213, and λ_πσ = 0.155.** The small difference to the values given in
reference [30] may be attributed to the different first-principles methods used in the calculation of the Eliashberg functions.」
（Fig.1 图内标签同四值。）

p.5 逐字：「The σσ-, ππ- and σπ-values were in the ratio 2.23/2.48/1. … The μ*(ω_c) matrix elements determined by this procedure are
**μ*_σσ(ω_c)=0.210, μ*_σπ(ω_c)=0.095, μ*_πσ(ω_c)=0.069, and μ*_ππ(ω_c)=0.172.**」（Tc=39.4 K 拟合出唯一自由参数。）

算术（输入 L1）：`λ_σπ/λ_πσ = 1.3742`，`N_π/N_σ = 1.3667`（差 0.5%）⟹ 右指标 DOS，与 Mitrović 的一般陈述一致；
`μ*_σπ/μ*_πσ = 1.377` 亦同。判据代入：`λ`：0.4556 vs 0.03302（通过 13.8×，与本仓已有数一致）；
`μ*`：0.03612 vs 0.006555（通过 5.51×）。
**定义式**：本文未印 `λ_ij` 的显式定义（「More details can be found elsewhere [37]」）；Eq.(6) 的
`N_σ(0)(I_σσ+I_σπ)+N_π(0)(I_ππ+I_πσ)` 只是塌缩式，按 B.3 结论**不定约定**；定约定的是上面那对比值。

## E.4 EPW 文档（`docs.epw-code.org/Theory/Midgal-EliashbergTheory.html`，WebFetch 抽文本，非页图）

Eq.(60)：`α²F_{nk,mk+q}(ω) = N_F Σ_ν |g_{mnν}(k,q)|² δ(ω−ω_qν)`；FSR 方程带 `πk_BT/N_F Σ_{mk'} δ(ε_mk'−ε_F)` 权重。
文档**不定义**纯带指标的 `λ_ij`。
（算术，我的推导，非文档原句）：对 `k∈i` 取 `(1/N_i)Σ` 平均、`k'∈j` 按 FSR 权重求和 ⟹ `λ_ij = (1/N_i)Σ_{k∈i,k'∈j}δδ|g|² = N_j⟨|g|²⟩_ij`
⟹ 与 Sanna Eq.(11) 同形，`f(i)g(j)` 形。EPW 产出的带分辨表（如 Margine–Giustino MgB₂）**未读**，此推导未与其对照。

## E.5 ★ 找到的非秩一映射：Morel–Anderson 截断重整化的矩阵版（Mitrović p.2–3，L1 逐字）

Eq.(12)：`μ̂*(ω_c) = (1̂ + μ̂^s ln(E_F/ω_c))⁻¹ μ̂^s`，两带显式（Eqs.13–16）：

```
μ*_σσ(ω_c) = [μ^s_σσ + (μ^s_σσ μ^s_ππ − μ^s_σπ²) ln(E_F/ω_c)] / D
μ*_σπ(ω_c) =  μ^s_σπ / D
μ*_ππ(ω_c) = [μ^s_ππ + (μ^s_σσ μ^s_ππ − μ^s_σπ²) ln(E_F/ω_c)] / D
D = 1 + (μ^s_σσ + μ^s_ππ) ln(E_F/ω_c) + (μ^s_σσ μ^s_ππ − μ^s_σπ²) (ln(E_F/ω_c))²
```

（我重推了 `(1+AL)⁻¹A` 的 2×2 展开，与印出的四式逐项一致。）

**算术（不是判定）**：记 `L = ln(E_F/ω_c)`、`Δ = ac − b²`（`a=μ^s_σσ, b=μ^s_σπ, c=μ^s_ππ`）。
- `det μ̂* = Δ/D` ⟹ **2 带时判据 `μ*_11μ*_22 ≥ μ*_12μ*_21` 的符号与裸 `μ̂^s` 相同，当且仅当 `D>0`**；
- 但比值 `μ*_12μ*_21/(μ*_11μ*_22) = b²/((a+ΔL)(c+ΔL))` **随 `L` 变**（裸值是 `b²/(ac)`）⟹ **不是** `f(i)g(j)×` 重标（后者比值不变）；
- 3 带以上：`(1+AL)⁻¹A` 的各 2×2 主子式与 `A` 的**没有**这种简单关系（我未推，只标记）。
⟹ 这是文献里**同一物理量的两个约定（裸 `μ` 与 `μ*(ω_c)`）之间**唯一一处不落在 `f(i)g(j)` 形状里的映射。它在**库仑侧**，
不在派活方问的声子 `λ_ij` 上；对本仓 `μ11μ22 ≥ μ12μ21` 判据的含义由派活方判。

## E.6 同一篇内两种挂法（结构观察，重述 B.8/§8.3，不新增来源）

Sanna PRB 85,184514：`λ_JJ'` 挂右指标（Eq.(11)）、`μ*_JJ' = μ*N_J/N` 挂左指标（拟合 ansatz）。
⟹ 组合核 `λ_JJ' − μ*_JJ' = N_J'V_JJ' − μ*N_J/N` **不是单一 `f(i)g(j)` 对某个 `W` 的重标**（两项前因子不同）。
但派活方的判据只对 `λ` 单独用，且任何矩阵都平凡地满足 `f=g=1`；这条**只是结构观察**，是否重要归派活方。

## E.7 转引（未读原文）：Mazin–Antropov 多带 McMillan

> **⚠ 本节已于第十轮被附录 F 取代，且其中的形状判断被我自己订正。原文照录如下，勿再引用其结论。**

检索摘要给出 cond-mat/0212263 的 `Λ^eff_ij = [Λ_ij − μ*_ij(1 + 0.62Σ_nΛ_in)]/(1 + Σ_nΛ_in)`（**转引，未读页图**）。
形状：`f(i) = 1/(1+λ_i)` 乘以 `[Λ_ij − c(i)μ*_ij]`；若 `Λ, μ*` 同一挂法则整体仍是 `f(i)g(j)W_ij`（`W_ij = V_ij − c(i)K_ij`，不对称但按 §8.2 不需要对称）。
未核实，不承重。

> **订正（F.2）**：上面这句「整体仍是 `f(i)g(j)W_ij`」在**字面上**没错（任何矩阵都能这样写，取 `f=g=1`），
> 但它**答非所问**：要问的不是「能不能写成那个形状」，而是「**约定变换在这个对象上还是不是秩一重标**」。
> 对 `Λ^eff` **不是**——因为 `Σ_n Λ_in` 在 `Λ_ij → f_i g_j Λ_ij` 下**不协变**（见 F.2 的代数与 F.3 的数值演示）。
> 我当时把「写得成」当成了「协变」，这两件事不同。原文已取到（F.1），形状判断以 F.2 为准。

## E.8 本附录的未核实项

- Mitrović 的期刊卷页（EPJ B？）未核；只读 pp.1–3。
- Dolgov 2009 页图希腊字母丢失，`λ` 由上下文补入；`B_12/B_21 = N_2/N_1` 这一句字母完整。
- Golubov 的 `λ_ij` 显式定义在其 Ref.[37]，未取；「右指标」由 `λ_σπ/λ_πσ = N_π/N_σ` 数值反推（0.5%），不是读到的定义式。
- EPW 带分辨表（Margine–Giustino PRB 87,024505）未读；E.4 的推导未与其对照。
- Mazin–Antropov 只有检索摘要，未读。 → **第十轮已取全文，见附录 F。**
- 检索覆盖：4 次 WebSearch（两带/多带 Eliashberg 定义、`(N_i+N_j)` 对称化、总 DOS 归一化、`λ_ij/(1+λ_i)` 重整化）；
  「未找到非乘积形状」是**这四次检索 + 五篇一手的阴性**，不是「不存在」。

---

# 附录 F（2026-09-18 第十轮）：Mazin–Antropov 与 Mazin–Golubov 两篇原文取到——**声子 `λ_ij` 定义仍全在乘积形状里；但「约定变换 = 秩一重标」在 `Λ^eff` 上失效，且失效**只在 `μ*>0` 时发生

> **来由**：第九轮 E.7 只有检索摘要。本轮取到两篇原文页图，**并订正 E.7 自己的形状判断**。
> 本附录只取件 + 算术；**是否「威胁约定无关定理」、以及对本仓消费点的含义，是派活方的判定。**
> 本轮新增一手 2 篇（累计 7 篇），新增 WebSearch 2 次（累计 6 次）。

## F.0 一句话结果

1. **声子侧仍然阴性**，且证据更强：本轮两篇都**逐字印出** `Λ_ij = (对称核)_ij × N_j`（F.1.a、F.4.a），累计 7 篇一手无一例外。
2. **★ 真正的反例不在「定义」里，在「派生矩阵」里**：Mazin–Antropov Eq.(7) 的多带 McMillan 有效矩阵
   `Λ^eff_ij = [Λ_ij − μ*_ij(1 + 0.62 Σ_n Λ_in)]/(1 + Σ_n Λ_in)`——
   它含**行和** `Σ_n Λ_in`，而行和在 `Λ_ij → f_i g_j Λ_ij` 下**不协变**（`→ f_i Σ_n g_n Λ_in`，不是 `f_i g_i Σ_n Λ_in`）
   ⟹ 由 `Λ^eff` 构成的比值 `Λ^eff_ijΛ^eff_ji/(Λ^eff_iiΛ^eff_jj)` **随约定改变**（F.2 代数，F.3 数值演示 5.3%）。
   **分母 `1+Σ_nΛ_in` 本身无害**（它是纯 `f(i)`，在比值里精确相消）；**非协变全部来自乘在 `μ*` 上的 `1+0.62Σ_nΛ_in`**
   ⟹ **`μ* = 0` 时 `Λ^eff` 的比值与 `Λ` 的逐位相同，威胁只存在于 `μ*>0`。**
3. **顺带**：Mazin–Antropov Eq.(8) 是 E.5 那个非秩一映射（多带 Morel–Anderson）的**第二个独立一手来源**，且更强——
   它的 `log(W_n/ω_c)` 里的带宽 `W_n` **带指标**（F.5）。
4. **顺带**：Mazin–Golubov 1996 p.2 印出**全排斥耦合仍有超导解**的条件（F.4.b）；我独立算出的 2×2 精确条件恰是
   `Λ_12Λ_21 > Λ_11Λ_22`，即**判据的否定**——但与印出的那一行**不完全一致**，差异**未裁决**（F.4.b ⚠）。
5. **顺带**：该文 p.2 的教科书式「各向异性增益最大」例子 `Λ_11=Λ_22=Λ>0, Λ_12=Λ_21=−Λ` **恰好落在判据的等号上**（F.4.c）。

## F.1 Mazin & Antropov, cond-mat/0212263v1 = Physica C 385, 49 (2003)（L1 页图，14 页全读）

> 取件路线：`arxiv.org/pdf/cond-mat/0212263` → WebFetch 拒绝抽文本但落盘 → `Read(pages=)` 读页图。
> 标题 "ELECTRONIC STRUCTURE, ELECTRON-PHONON COUPLING, AND MULTIBAND EFFECTS IN MgB₂"，
> I.I. Mazin（Code 6391, NRL）与 V.P. Antropov（Ames Lab），(10/11/2002)。正文 pp.1–12，参考文献 pp.13–14（[1]–[68]）。

### F.1.a 定义式（p.6 右栏，§III A，L1 逐字）

> which can be broken into n × n matrix separating the interband pairing interaction from the intraband one
>
> α²F(ω)_ij = (1/N_i(E_F)) Σ_{q,ν} Σ_{k∈i, k+q∈j} |g_{k,k+q,ν}|² δ(ε_{k+q} − ε_k − ħω_{q,ν}) δ(ω − ω_{q,ν}),
>
> where i, j label different electronic bands or group of bands, e.g., i = σ, π.

⟹ `1/N_i` **只除左指标**，与 Sanna Eq.(11)（附录 B.2）**逐字同构** ⟹ `λ_ij = N_j V_ij`（DOS 挂右指标）。
⚠ 我照抄了印出的 δ 宗量 `δ(ε_{k+q} − ε_k − ħω_{q,ν})`（与 Sanna 的 `δ(ξ_k)δ(ξ_{k'})` 不同形）；
**这是页图上就是这么印的**，我不推测它是否为排印问题，也不据此做任何判定。

### F.1.b p.9 右栏 §IV A（L1 逐字）——同一事实的第二次印出，且给出对称化

> where i, j are the band indices and Λ is an *asymmetric* matrix related to the *symmetric* matrix of the pairing
> interaction, **Λ_ij = V_ij N_j**, where N_i is the contribution of the i-th band to the total DOS. … Note that although
> the matrix Λ is not symmetric, its eigenvalues are the same as those of the symmetric matrix **√N V √N**.

⟹ 两种写法（`V_ij N_j` 与 `√(N_iN_j) V_ij`）**都是 `f(i)g(j)`**；这是 E.1 Mitrović 对称化的第二个一手来源。

### F.1.c p.10 左栏（L1 逐字）——行和的物理含义，**这条使约定在该文里是被钉死的**

> The mass renormalization parameters for each band can be constructed from the matrix Λ: **λ_i = Σ_j Λ_ij**. These λ_i
> define, among other things, the de Haas-van Alphen thermal masses.

⟹ 行和 `Σ_j Λ_ij` 被赋予**可测量**（dHvA 热质量）⟹ 在这篇里 `Λ_ij = V_ij N_j` 不是自由约定。
⚠ 同句之后印有一条塌缩式（`λ̄ = …`），**该行在我的页图分辨率下指标读不准，我不逐字引用、也不使用**。

同栏另一句（L1 逐字，与本仓 §4 的 (S1)/(S2) 机制直接相关）：

> One can show that **λ ≤ λ_eff, the equality being achieved when and only when all elements of the V matrix are the same**
> (the relative magnitude of N_i is irrelevant).

### F.1.d ★ p.10 右栏 §III C（应为 IV C）Eq.(5)–(7)（L1 逐字）——**本轮的主对象**

> kT_c = (ℏω_log/1.2) exp[ −1.02(1+λ) / (λ − μ*(1+0.62λ)) ].   (5)
> … The multiband analog of the McMillan equation is, therefore,
> kT_c = (ℏω_log/1.2) exp[ −1 / (λ − μ*)_eff ],   (6)
> where (λ − μ*)_eff is defined as the maximum eigenvalue of the matrix
> **Λ^eff_ij = [Λ_ij − μ*_ij(1 + 0.62 Σ_n Λ_in)] / (1 + Σ_n Λ_in).**   (7)
> This expression gives the results very close to the full solution of the multiband Eliashberg equations.

紧接一句（p.10 右栏末 / p.11 左栏，L1 逐字）：

> The Coulomb pseudopotential matrix is not a constant, as it is sometimes believed. … Indeed, it is formally defined as
> **⟨⟨V_C⟩⟩_ij N_j** (where V_C is the screened Coulomb interaction, and the averaging is over the corresponding Fermi surfaces)

⟹ **`Λ` 与 `μ*` 在这篇里挂法相同（都是右指标 `N_j`）**——这一点是 F.2 代数的前提，且由原文逐字给出，不是我的假设。

### F.1.e p.9 右栏 TABLE II（L1 逐字，MgB₂ 三套第一性原理 2×2 EPC 矩阵）

> TABLE II. 2x2 EPC matrices in different calculations.
>
> | Ref. [10] | Ref. [55] | Ref. [30]ᵃ |
> |---|---|---|
> | 0.96  0.17 | 1.02  0.16 | 0.78  0.11 |
> | 0.23  0.29 | 0.21  0.45 | 0.15  0.21 |
>
> ᵃ Obtained by integrating λ(k,k′) distribution plots from Ref. [30].

判据代入（算术，输入 L1；判据转置不变故与排版歧义无关）：

| 来源 | `λ11λ22` | `λ12λ21` | 判决 |
|---|---|---|---|
| Ref. [10]（Liu–Mazin–Kortus） | 0.2784 | 0.0391 | pass **7.12×** |
| Ref. [55]（Golubov JPCM 14,1353） | 0.4590 | 0.0336 | pass **13.66×** |
| Ref. [30]（Choi et al.） | 0.1638 | 0.0165 | pass **9.93×** |

⚠ **一处与本文自身约定的张力，我只记录、不裁决**：Ref.[55] 就是 E.3 那篇，其原文印 `λ_σπ = 0.213 > λ_πσ = 0.155`；
而 Table II 的 (1,2) 格印 **0.16**、(2,1) 格印 **0.21**，即**较小者在 (1,2)**。若行列序为 (σ,π)（由对角 1.02/0.45 定），
则按本文自己的 `Λ_ij = V_ij N_j` 与 `N_π > N_σ` 应有 `Λ_σπ > Λ_πσ`，(1,2) 应是较大者。三列**内部彼此一致**（都把较小者放 (1,2)）。
可能是排版转置、也可能是我对行列序的读法错；**判据对转置不变 ⟹ 上表三行判决不受影响**，故不追。

## F.2 ★ 代数：`Λ^eff` 在约定变换下**不协变**——订正 E.7（算术，我从零推）

设约定变换为 `Λ_ij → f_i g_j Λ_ij`，且（按 F.1.d 的逐字）`μ*_ij → f_i g_j μ*_ij`（同一 `f,g`）。记

```
S_i := Σ_n Λ_in ,   A_i := 1 + 0.62 S_i ,   Λ^eff_ij = (Λ_ij − μ*_ij A_i)/(1 + S_i)
```

变换后：

```
S_i  →  Σ_n f_i g_n Λ_in = f_i · S̃_i ,   S̃_i := Σ_n g_n Λ_in           ← 一般地 ≠ f_i g_i S_i
Λ^eff_ij  →  f_i g_j · [Λ_ij − μ*_ij Ã_i] / (1 + f_i S̃_i) ,  Ã_i := 1 + 0.62 f_i S̃_i
```

前因子 `f_i g_j` 与分母 `1/(1+f_iS̃_i)`（纯 `f(i)`）在比值里**精确相消**，于是

```
R_ij := Λ^eff_ij Λ^eff_ji / (Λ^eff_ii Λ^eff_jj)
      = [(Λ_ij − μ*_ij A_i)(Λ_ji − μ*_ji A_j)] / [(Λ_ii − μ*_ii A_i)(Λ_jj − μ*_jj A_j)]
```

**`R_ij` 只通过 `A_i, A_j` 依赖约定，而 `A_i` 不协变** ⟹ `R_ij` **随约定改变**。三条随行结论：

- **(i) `μ* = 0` 时完全没问题**：`Λ^eff_ij = Λ_ij/(1+S_i)` 是纯 `f(i)` 重标 ⟹ `R_ij = Λ_ijΛ_ji/(Λ_iiΛ_jj)`，与 `Λ` 的**逐位相同**。
  **⟹ 非协变**只*在 `μ*>0` 时出现。
- **(ii) 质量重整化分母不是原因**，`0.62` 那一项才是：非协变的唯一入口是乘在 `μ*_ij` 上的 `A_i`。
- **(iii) 在这篇文章里约定并不自由**（F.1.c：行和 = dHvA 质量重整化钉死了 `Λ_ij = V_ij N_j`）
  ⟹ 这**不是**「已发表数字有歧义」，而是「**约定无关定理覆盖 `λ_ij`，不自动覆盖 `Λ^eff`**」。
  若某处消费点吃的是 `λ` 本身，定理照用；若吃的是 `λ − μ*` 型或 `Λ^eff` 型的组合核，**定理要重证**。这一句是取件方的陈述边界，判定归派活方。

## F.3 数值演示：同一组已发表数字，两种挂法给出**不同**的 `R`（算术，输入 L1；**这是我的构造，不是任何论文的计算**）

用 E.3 的 Golubov MgB₂ 两张已发表表（同一篇、同一体系）代入 Mazin–Antropov Eq.(7)：
`Λ = [[1.017, 0.213],[0.155, 0.448]]`，`μ* = [[0.210, 0.095],[0.069, 0.172]]`（行列序 σ,π）。

| | `S_σ` | `S_π` | `A_σ` | `A_π` | `R_σπ` | 判决 |
|---|---|---|---|---|---|---|
| 原挂法 | 1.230 | 0.603 | 1.7626 | 1.37386 | **0.02003** | pass 49.9× |
| 两表同时转置（= 相反挂法） | 1.172 | 0.661 | 1.72664 | 1.40982 | **0.02108** | pass 47.4× |

⟹ **差 5.3%**（比值本身，不是判决）。对照：**裸 `Λ` 的比值两种挂法逐位相同**，`0.033015/0.455616 = 0.072463`（= E.3 的 13.80×）
——这个不变性正是约定无关定理，它在 `Λ` 上成立、在 `Λ^eff` 上不成立。

⚠ 三条边界：① 两种挂法的**判决**在这个例子里都是 pass、没有翻转，我**不**主张能翻转；
② 把 Golubov 的 `λ` 与 `μ*` 一起喂进 Mazin 的 Eq.(7) 是**我**的组合（两表同篇同体系，但该文本身没这么算）；
③ 非协变的**幅度界**我没推，只给了一个点。手算 float，未上机（本会话无 Bash）。

## F.4 Mazin & Golubov, cond-mat/9602118v1（L1 页图 pp.1–4；1996 年 1 月，未查期刊出处）

> "Impurity scattering in highly anisotropic superconductors and inerband sign reversal of the order parameter."（标题即印作 "inerband"）
> I.I. Mazin（Carnegie Institution of Washington）、A. Golubov（ISI, FZ Jülich / Chernogolovka）。

### F.4.a 定义式（p.2，L1 逐字）——第 7 篇一手，仍是乘积形状

> The effective coupling constant λ_eff in this case is simply the maximal eigenvalue λ_max of the matrix
> **Λ_ij = g_ij N_j**, where N_j is the density of states at the Fermi level (per spin) in the j-th band. Λ_ij plays the
> role of the coupling constant λ in the one-band BCS theory. Note that conventional (isotropic) λ is also defined in
> terms of Λ_ij: **λ = Σ_ij Λ_ij N_i/N = Σ_i λ_i N_i/N**, where the mass renormalization for the i-th band, as measured,
> for instance, in de Haas-van Alphen experiments, is **λ_i = Σ_j Λ_ij**, and N = Σ_i N_i.

p.4 同一事实第三次（L1 逐字）：

> Coupling matrix Λ is defined in the same way as Allen's matrix λ_JJ′: **Λ_JJ′ = V^pairing_JJ′ N_J′(0)**.

> ⚠ **订正（第十二轮）**：上面这个 `N_J′(0)` **下标抄错了**，应是**两个指标** `N_{JJ′}(0)`。
> 派活方在真文本层做 span 级字号读（`N` 9.96 / `JJ` 6.97 / `′` 4.98 / `(0).` 9.96，与同行 `Λ_JJ′` 同模式）判为两指标；
> 我重读同一页的页图，**在这个分辨率下无法区分 `N_{J′}` 与 `N_{JJ′}`**（下标宽度与同行 `V^pairing_{J,J′}` 目视同宽），
> ⟹ **文本层 span 是更好的仪器，我采纳该订正、撤回我的单指标读法。**
> 语义见 **H.3**：`N_{JJ′}` 的两种可能读法**都**落在 `f(i)g(j)` 形状里 ⟹ 对约定无关定理无威胁（但语义本身仍未逐字定案）。

⚠ p.3 的库仑修正是**在 `g` 那一层**做减法（逐字）：

> If g's are electron-phonon pairing potentials, then Eq.1 should be corrected for a Coulomb repulsion, which can be
> readily done by substituting **g_ij ⟶ g_ij − U*_ij ≈ g_ij − U***, where the effective Coulomb repulsion U* is
> logarithmically renormalized in the same way as in one-band superconductivity theory (U* is assumed to be independent on i,j).

⟹ `Λ_ij = (g_ij − U*)N_j` **仍是 `f(i)g(j)` 形**（减法发生在对称核内部）。
**与 F.1.d 的 Eq.(7) 形成对照**：坏掉形状的不是「减 μ*」，而是「μ* 前面那个**带指标的行和因子**」。

### F.4.b p.2 逐字：全排斥耦合下仍有超导解的条件

> One can easily check that Eq.1 may have a superconducting solution even for all g_ij < 0, i.e., when no attractive
> interaction is present in the system. **The condition for that is |g_12| > (|g_1|N_1² + |g_2|N_2²)/2N_1N_2.**

⚠ `g_1`/`g_2` 的下标在我的页图分辨率下就是印作单指标；**我不猜它是否指 `g_11`/`g_22`**。

**我独立推的 2×2 精确条件（算术）**：设 `Λ = [[a,b],[c,d]]`，`a,d < 0`，`bc > 0`（全排斥），则
`λ_max = ((a+d) + √((a−d)² + 4bc))/2 > 0 ⟺ bc > ad`，即 **`Λ_12Λ_21 > Λ_11Λ_22`——恰是本仓判据的否定**。
以 `Λ_ij = g_ij N_j` 代回：`⟺ g_12² > g_11 g_22 ⟺ |g_12| > √(|g_11||g_22|)`。
而印出的右端按 AM–GM 满足 `(|g_1|N_1² + |g_2|N_2²)/(2N_1N_2) ≥ √(|g_1||g_2|)`，等号仅当 `|g_1|N_1² = |g_2|N_2²`
⟹ **印出的条件比我推的精确条件更强（充分不必要），二者不等价**。
**差异未裁决**：可能是印刷、可能是我把单指标 `g_1` 读成 `g_11` 读错了、也可能该句针对的是别的（如带上截断的完整解）。
**我只报「两者不一致」这一事实**；是否以及如何使用，归派活方。

### F.4.c p.2 逐字：教科书式「增益最大」例子**恰在判据等号上**

> An instructive example of the opposite case is the two-band model with **Λ_11 = Λ_22 = Λ > 0, Λ_12 = Λ_21 = −Λ.
> Then λ = 0, while λ_eff = 2Λ.**

算术：`Λ_11Λ_22 = Λ²`、`Λ_12Λ_21 = Λ²` ⟹ **恰好取等**。该矩阵 `= Λ·vvᵀ`（`v = (1,−1)`）是秩一 ⟹ 与本仓 §4 的 **(S1) 可分取等**同一现象。
（这是第三次现身：§4 的 `g_{kk'}=γ_k\bar γ_{k'}`、§8.3 的半各向同性库仑 ansatz、此处。）

## F.5 Eq.(8)：多带 Morel–Anderson 的第二个一手来源，且比 E.5 更强（p.11 左栏，L1 逐字）

> For a multiband case we have a matrix equation, which is a natural extension of the standard procedure [65]
>
> **μ*_ij = μ_ij − Σ_n μ_in log(W_n/ω_c) μ*_nj.**   (8)
>
> It is easy to show that renormalization enhances any nonuniformity in μ; indeed, assuming μ_σσ = μ_ππ = αμ_σπ, (α>1),
> and μ_σσ log(W_σ/ω_c) = L, we obtain **α* = α + (α − 1/α)L**. From the ratios of ⟨⟨V_C⟩⟩'s above, α ≈ 2.3, and L for
> MgB₂ is of the order of 0.5 - 1, so for μ*_ij it holds that **μ*_σσ = μ*_ππ ∼ 4μ*_σπ**.

- 这是 E.5（Mitrović `μ̂*(ω_c) = (1̂+μ̂^s L)⁻¹μ̂^s`）的**第二个独立一手来源**：Eq.(8) 就是它的线性方程写法。
- **更强的一点**：这里的对数是 `log(W_n/ω_c)`，**带宽 `W_n` 带指标**——被夹在 `μ_in` 与 `μ*_nj` 之间求和
  ⟹ 连「乘一个标量 `L`」都不是，更不可能是 `f(i)g(j)` 重标。
- 算术核对（输入 L1）：`α=2.3, L=1 ⟹ α* = 2.3 + (2.3 − 0.4348) = 4.165`，与印出的「∼4」自洽。
- 算术（不是判定）：若 `μ` 对称且 `μ_σσ=μ_ππ=αμ_σπ`，判据比值 `= 1/α²`
  ⟹ 裸 `α≈2.3` 给 **5.29×** 余量、重整化后 `α*≈4` 给 **16×** 余量 ⟹ **该重整化把库仑矩阵推离判据边界**（在这组参数上）。
  ⚠ 对照 E.3：Golubov 拟合出的 `μ*` 四值给 **5.51×**，更接近**裸** `α` 而非 `α*`——两篇的 `μ*` 定义与拟合过程不同，**不裁决**。

## F.6 本附录的未核实项

- Mazin–Golubov cond-mat/9602118 的期刊出处未查（只读 pp.1–4，全文 pp.1–? 未定；第 4 页仍在正文中）。
- F.1.a 的 δ 宗量、F.1.c 的塌缩式、F.4.b 的 `g_1/g_2` 下标：三处我按「印出什么抄什么 / 分辨率不足就不抄」处理，**均未推测**。
- F.1.e 的 Table II 排版转置张力未追（判据转置不变 ⟹ 与判决无关）。
- F.3 是手算 float64 级、单点、无脚本（本会话无 Bash）；非协变的幅度界未推。
- Ref.[58]（Allen & Mitrović, Solid State Phys. 37, 1）——Eq.(2)/(3) 与 Eq.(5)–(7) 的原始出处——**未取**。
- 检索覆盖累计：**6 次 WebSearch**、**7 篇一手**（Sanna / Essenberger / Ummarino / Mitrović / Dolgov / Golubov / Mazin–Antropov + Mazin–Golubov）。
  「声子 `λ_ij` 定义未找到非乘积形状」仍是**这些检索的阴性，不是「不存在」**。
  本轮另见两条**未取**的候选线索（检索摘要级，未读）：arXiv:1809.07000（Chevrel，据称有显式 `1/N_i(0)` 定义）、
  arXiv:cond-mat/0310284（Mg₁₋ₓAlₓB₂，据称有 `λ_ij/λ_ji` 关系）。

---

# 附录 G（第十一轮）：Golubov 的 `2.23/2.48/1` 有第三方对照值，两种读法被它指向相反的一侧

> **缘起**：主会话第九轮把 Golubov p.5 的 `2.23/2.48/1` 判为**逐行归一**（`μ*_σσ/μ*_σπ=2.2105` vs 2.23、
> `μ*_ππ/μ*_πσ=2.4928` vs 2.48，均 <1%；单一参照读法差 26.99%），并写进笔记 §9 / 脚本 `(F)` / 记忆第 25 条。
> 我手上有它当时没有的一手件（Mazin–Antropov 全文，第十轮取）。核对之后发现：
> **该篇印出了同一体系裸 `⟨⟨V_C⟩⟩` 的第三方比值，而那个比值指向另一种读法。**
> 本节只取件 + 算术，**不裁决**；判归主会话。

## G.1 [L1 逐字] Mazin–Antropov cond-mat/0212263 p.11 左栏

> …when different bands have different orbital character, the Coulomb matrix elements between these bands are
> suppressed compared to intraband matrix elements. **Jepsen and Andersen [64] estimated this effect, using the
> tight-binding LMTO method and found the ratio between ⟨⟨V_C⟩⟩_σσ, ⟨⟨V_C⟩⟩_ππ and ⟨⟨V_C⟩⟩_σπ to be ≈3:2.5:1.**
> Furthermore, any anisotropy in *bare* pseudopotential is further enhanced in the renormalized μ*_ij. …

`⟨⟨V_C⟩⟩` 是**裸、对称、无 DOS** 的费米面平均屏蔽库仑矩阵元（同页上文：`μ_ij` "is formally defined as
`⟨⟨V_C⟩⟩_ij N_j`"，见 F.1.d）。⟹ **3 : 2.5 : 1 是 `V` 层的比值，不是 `μ` 层的。**

**该篇 [64] 逐字（p.14）**：`I. I. Mazin, O. K. Andersen, O. Jepsen, A. A. Golubov, O. V. Dolgov, and J. Kortus, to be published.`

## G.2 [L1 逐字] Golubov cond-mat/0111262v2 印刷 p.5 的**完整**句子（比第九轮所引更长）

> Besides the spectral functions we need to know the the Coulomb matrix element **μ_ij**. With the help of the
> wavefunctions from our first-principles calculations we can approximately calculate **the ratios for the μ-matrix [37]**.
> The σσ-, ππ- and σπ-values were in the ratio **2.23/2.48/1**. This allows one to express **μ*_ij(ω_c)** by these ratios
> and one single free parameter which is fixed to get the experimental T_c of 39.4 K from the solution of the
> Eliashberg equations. The μ*(ω_c) matrix elements determined by this procedure are
> **μ*_σσ(ω_c)=0.210, μ*_σπ(ω_c)=0.095, μ*_πσ(ω_c)=0.069, and μ*_ππ(ω_c)=0.172.**

**该篇 [37] 逐字（p.9）**：`Kuz'menko A B et al 2001 Preprint cond-mat/0107092; Solid State Commun. in press 2002;
Dolgov O V et al 2002 to be published.`
⟹ 与 G.1 的 [64] **不是同一条著录**，但两条都是同一批人（Golubov/Dolgov/Andersen/Jepsen/Mazin/Kortus）的
「to be published」；**是否同一份计算，我判不了，不裁决。**

⚠ **转录分辨率**：`σσ-, ππ- and σπ-` 这三个下标的字号很小，我这次重读时无法在页图上单独确证顺序。
记录为「两次独立阅读一致」（我第九轮 E.3 + 主会话第九轮脚本 docstring 各自抄成同一顺序）。
另一种顺序（`σπ-, ππ-, σσ-`）会使 `σσ` 成为最小者，与印出的 `μ*_σσ=0.210` 为最大值直接冲突 ⟹ 与本节结论无关。

## G.3 [算术，输入 L1] 把三组数放到同一个 `V` 层上比

约定：`μ*_ij = c·N_j V_ij`，`V` 对称（该约定由印出数自证：`μ*_σπ/μ*_πσ = 0.095/0.069 = 1.37681` vs
`N_π/N_σ = 0.410/0.300 = 1.36667`，**0.74%**）。`N_σ=0.300, N_π=0.410`（同篇 p.4）。

| 来源 | `V_σσ : V_ππ : V_σπ` | 与 G.1 的 `≈3:2.5:1` 的差 |
|---|---|---|
| **印出的四个 `μ*`**（不依赖任何读法） | `0.210/0.300 : 0.172/0.410 : 0.095/0.410` = **3.0211 : 1.8105 : 1** | σσ **+0.70%**；ππ **−27.6%** |
| 句子按**逐行归一**读（主会话第九轮） | `2.23·N_π/N_σ : 2.48·N_σ/N_π : 1` = **3.0477 : 1.8146 : 1** | σσ +1.59%；ππ **−27.4%** |
| 句子按**单一参照**读（字面语法） | `2.23·N_π/N_σ : 2.48 : 1` = **3.0477 : 2.4800 : 1** | σσ +1.59%；ππ **−0.80%** |
| G.1 第三方（Jepsen–Andersen via [64]） | **3 : 2.5 : 1** | — |

**判别量只有一个**：`V_ππ/V_σπ`。两种读法在 `σσ` 上给**同一个数**（该比在两种读法下同式），
差别全部落在第二项：**2.48 ↔ 2.5（0.80%）vs 1.81 ↔ 2.5（27.6%）**。

**内部逻辑**（不是我的假设，是 G.2 那句话自己说的）：比值是**裸 `μ` 矩阵**的，`μ*` 由「这些比值 + 一个自由参数」
直接给出 ⟹ `μ* ∝ μ` 逐元 ⟹ **句子的读法与印出的四个 `μ*` 必须自洽**。于是只剩两支：

- **(i)** 句子是逐行归一、四个 `μ*` 全对 ⟹ 该篇的 `V_ππ/V_σπ = 1.81`，与 [64] 的 2.5 差 **27.6%**；
- **(ii)** 句子是单一参照 ⟹ `μ*_ππ` 应为 `2.48×0.095 = 0.2356`（印出 0.172，差 27%）⟹ 该篇的 `V` 与 [64] 差 1.6%/0.8%。

**我不裁决 (i)/(ii)。** 两支都要求某处有 27% 的不一致，只是位置不同（一个在跨文献、一个在同篇内）。

## G.4 [算术] 对判决与余量的影响

| 量 | 值 | 依赖读法？ |
|---|---|---|
| `μ*` 判据余量（印出四数直接算） | `0.210×0.172/(0.095×0.069) = **5.5103×**` | **否** |
| 主会话的「免转录形式」`2.23×2.48` | `5.5304` | **是**（只在逐行归一下成立） |
| 支 (ii) / [64] 比值下的余量 | `3.0477×2.48 = **7.558×**`（取整 `3×2.5 = 7.5×`） | 是 |

⟹ **判决不动**：三个数全部 **PASS**，MgB₂ 仍是通过的一例。受影响的只有**余量数值**（5.51× ↔ 7.5×，差 37%）
与两条**派生**陈述：主会话脚本 `(F)` 节印出的 `implied symmetric V^c`（`V_σσ/V_σπ=3.0473, V_ππ/V_σπ=1.8146`、
`V_σσ/V_ππ=1.679`）和记忆第 25 条的「歧义由数字定案」——后者在**同篇内部**成立（我复核了 0.87%/0.51% 与 26.99%），
但把 G.1 的第三方 `V` 比值放进来之后不再是单侧证据。

## G.5 顺带（均已在 F.5 记过，此处只补一处交叉）

Mazin–Antropov 同页由 `α≈2.3`（其 `μ_σσ=μ_ππ=αμ_σπ` 的单一参照式）出发。
把 [64] 的 `V = 3:2.5:1` 代进 `μ_ij = N_jV_ij`：`μ_σσ/μ_σπ = (N_σ/N_π)·3 = **2.195**`、`μ_ππ/μ_σπ = **2.500**`
——两者恰是 Golubov 印出的 **2.23 / 2.48** 的邻居（1.6% / 0.8%），且那篇把它们**都**近似成一个 `α≈2.3`。
这是支 (ii) 的算术来源，**不是裁决**：`α` 是两个不同数的共同近似，"≈3:2.5:1" 本身也带 `≈`。

## G.6 本附录的未核实项

- G.1 的 [64] 与 G.2 的 [37] 是否同一份未发表计算——**未查**（两条著录首作者不同，均为 "to be published"）。
- Jepsen–Andersen 的原始 LMTO 数值**未取**（只有 Mazin–Antropov 的转述，且原文自带 `≈`）。
- G.2 三个下标的印刷顺序：**分辨率不足，未独立确证**，按两次独立阅读一致处理（见 G.2 ⚠）。
- 支 (ii) 下 `μ*_ππ=0.2356` 能否仍给 `T_c=39.4 K`——**无法检验**（本会话无 Bash，且需重解 Eliashberg 方程）。
- 本轮零次 WebSearch：两份 PDF 均为第十轮已落盘的本地件（`theory/notes/_scratch_ir/_mazin_antropov0212263.pdf`、
  `_golubov_cond-mat0111262.pdf`），新读页：Mazin–Antropov pp.10–11 / pp.13–14、Golubov p.5 / pp.8–9。

---

# 附录 H（第十二轮）：读法之争由第三方矩阵式关闭（站 (A)）；Allen 原文取到，`N_{JJ′}` 的两种读法都不威胁定理

> **本轮做了三件事**：① 取到一个**印出整个 `μ*` 矩阵**的第三方来源，读法之争关闭；
> ② 取到派活方点名「在拿到之前谁都别写进承重链」的 **Allen PRB 13, 1416 (1976)**；
> ③ 采纳派活方的四条订正/裁决并记录随行 OPEN。
> 本附录仍只取件 + 算术；对本仓消费点的含义归派活方。

## H.0 先说一句：我第十一轮举的证据指向 (B)，本轮的新证据指向 (A)。我站 (A)。

第十一轮我拿 Mazin–Antropov 转述的 Jepsen–Andersen 裸比 `≈3:2.5:1` 指出它支持**单一参照**读法 (B)，并明确不裁决。
本轮取到的 **Ummarino et al. 2004 Eq.(10) 把整个 `μ*` 矩阵逐元印了出来**，形状是**逐行归一**（(A)）。
⟹ **我不再为 (B) 举证。** 派活方的第十一轮裁决站得住，且它给的「值域不可达」论证与本轮的直接证据独立同向。
JA 的 `≈2.5` 与 (A) 蕴含的 `1.815` 之间那 27.4% **没有被解释掉**，仍是 OPEN（H.4）。

## H.1 ★ [L1 逐字] Ummarino, Gonnelli, Massidda, Bianconi, cond-mat/0310284v3, p.3 Eq.(10)

正文引导句逐字：「The Coulomb pseudopotential, calculated, for the first time, in **ref. 15 and 16** is **[22]**」

```
μ*(p) = | μ*_σσ   μ*_σπ |
        | μ*_πσ   μ*_ππ |

                                  | 2.23/N^σ_N(E_F,p)     1/N^σ_N(E_F,p)  |
      = μ(ω_c) N^tot_N(E_F,p)  ×  |                                        |     (10)
                                  |  1/N^π_N(E_F,p)     2.48/N^π_N(E_F,p) |
```

p.4 左栏逐字：「where `μ(ω_c,x)` is a free parameter and `N^tot_N(ω=E_F,x)` is the total normal density of states at the
Fermi level. **The numbers 2.23 and 2.48 in the Coulomb matrix have been calculated in the MgB₂ case** … In conclusion
we put the cut-off energy **ω_c = 700 meV**. So we have now only one free parameter, `μ(ω_c,x)`…」

**著录（p.6）**：[15] = `Golubov, Kortus, Dolgov, Jepsen, Kong, Andersen, Gibson, Ahn, Kremer, J. Phys.: Condens. Matter 14 (2002) 1353`（= 本仓那篇）；
[16] = `Brinkman, Golubov, Rogalla, Dolgov, Kortus, Kong, Jepsen, Andersen, PRB 65 (2002) 180517(R); Golubov, Brinkman, Dolgov, Kortus, Jepsen, PRB 66, 054524 (2002)`；
**[22] = `O.V. Dolgov, talk at the 6th EUCAS, Sorrento, Naples, Italy, September 14–18, 2003`**（会议报告，**取不到**）；
另 [18] = `I. I. Mazin, V. P. Antropov, Physica C 385 (2003) 49` ⟹ **F.1 的期刊著录得到第三方确认**。

### H.1.a [算术] 它把读法之争关掉了

Eq.(10) 里**每一行共用同一个 DOS 分母**（第 σ 行都是 `1/N^σ`、第 π 行都是 `1/N^π`）⟹ 行内比值直接可读：

| 量 | Eq.(10) | Golubov 印出的四个数 | 差 |
|---|---|---|---|
| `μ*_σσ/μ*_σπ` | **2.23** | 0.210/0.095 = 2.2105 | **0.87%** |
| `μ*_ππ/μ*_πσ` | **2.48** | 0.172/0.069 = 2.4928 | **0.51%** |
| `μ*_σπ/μ*_πσ` | `N^π/N^σ` = 1.3667 | 0.095/0.069 = 1.3768 | 0.74% |

⟹ **2.23 与 2.48 各自绑在自己那一行的非对角元上 = 逐行归一 = 读法 (A)**，由一个**独立印出整个矩阵**的来源直说，
不再依赖对 Golubov 那句话的语法解读。判据余量 `= 2.23 × 2.48 = **5.5304**`（与派活方的免转录形式逐位一致）。

⚠ **这份证据是什么、不是什么**：它**不是**独立的第一性原理计算（2.23/2.48 明说是从 [15]/[16] 拿的），
而是**由第三方写下的、对同一组数的矩阵化重述**；且 Eq.(10) 的出处 [22] 是 **Dolgov 本人的报告**——
Dolgov 是 Golubov 那篇的**共同作者**。⟹ 作为「原作者一方怎么理解这两个数」的证据，这已接近上限；
作为「独立数据」，它**不算**。两种定性我都写下，不合并。

### H.1.b [L1 逐字，p.2 右栏] 又一处 `λ` 与 `μ*` 共用同一 DOS 约定

> …as shown in ref. 14, **λ_ij(x)/λ_ji(x) = μ*_ij(x)/μ*_ji(x) = N^j_N(ω=E_F,x)/N^i_N(ω=E_F,x)** where
> `N^j_N(ω=E_F,x)` is the density of states at the Fermi level in the j-band…

（[14] = `S. D. Adrian, S. A. Wolf, O. V. Dolgov, S. Shulga, V. Z. Kresin, PRB 56 (1997) 7878`。）
**形状普查加一条**：Eq.(10) 的放置是 `μ*_ij ∝ 1/N_i`（**倒数、挂左指标**），是我见到的第 5 种放置，
但 `1/N_i = f(i)·1` 仍是 `f(i)g(j)` ⟹ **判据比值不变**；且它与 `∝N_j` 给出**同一个** `μ*_ij/μ*_ji = N_j/N_i`，
所以 H.1.b 那句话对两种放置都成立、**不定约定**（与 B.3 撤回过的「细致平衡定约定」同一陷阱）。

## H.2 ★ [L1 逐字] Allen, PRB **13**, 1416 (1976), §V「Eliashberg equations」，p.1423

取件路线：`harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.13.1416/fulltext`（**JBIG2 扫描件 ⟹ 必须读页图**）。

> A similar pair of equations relate `Z(kn)` to `Z_Ln`. Let us write out explicitly the expansion of the Eliashberg
> interaction `V_ep` Eq. (31):
>
> **α²_{LL′}F(Ω) ≡ (1/N(0)) Σ_{kk′} δ(ε_k)δ(ε_k′) |M_{kk′}|² × B(k−k′, Ω) ψ_L(k) ψ_{L′}(k′)**   (33)
>
> **V_ep(kk′, nn′) = −(1/N(0)) Σ_{LL′} λ_{LL′}(n−n′) ψ_L(k) ψ_{L′}(k′)**   (34)
>
> **λ_{LL′}(n−n′) ≡ ∫₀^∞ dΩ α²_{LL′}F(Ω) · 2Ω/(Ω² + (ω_n−ω_{n′})²)**   (35)

以及 p.1416 的正交归一定义（Eq.(1)）与 p.1423 右栏「The isotropic components (`L=0, L′=0`) … `λ₀₀(0)` is the usual
coupling constant and mass enhancement parameter λ」。

> **★ 关键事实**：Allen 的 `λ_{LL′}` 里**唯一的态密度因子是标量 `1/N(0)`**（费米面总态密度），
> **全式没有任何带指标 / 谐波指标的 DOS**。带结构只经 `ψ_L(k)ψ_{L′}(k′)` 进场。

## H.3 [算术，我的推导，非 Allen 原句] ⟹ `N_{JJ′}(0)` 的两种读法都是 `f(i)g(j)`

Mazin–Golubov p.4 自陈「without intraband anisotropy `J` can be chosen to be the band indices
(so-called **disjoint representation**)」。在该表示下由 Eq.(1) 定出 `ψ_J(k) = √(N/N_J)·1_{k∈J}`
（代入 Eq.(1)：`c_J² N_J/N = 1`）。把它代进 Eq.(33)，记 `Σ_{k∈J,k′∈J′}δδ|M|²B = N_J N_{J′}⟨|M|²B⟩_{JJ′}`：

　`α²_{JJ′}F = (1/N)·√(N/N_J)√(N/N_{J′})·N_J N_{J′}⟨|M|²B⟩_{JJ′} = **√(N_J N_{J′})·⟨|M|²B⟩_{JJ′}**`

⟹ Allen 的 `λ_{JJ′}` 在 disjoint 表示下是 **`√(N_J N_{J′}) × 对称核`**。于是：

| `N_{JJ′}(0)` 读作 | 形状 | 威胁约定无关定理？ |
|---|---|---|
| `N_{J′}(0)`（我原来的单指标读法） | `1·g(j)` | 否 |
| `√(N_J(0)N_{J′}(0))`（与 Allen 的 disjoint 表示自洽） | `f(i)g(j)`，`f=g=√N` | 否 |
| 真正不可分的两指标 `N_{JJ′}` | Hadamard | **是**——但 **Allen 的形式里没有这种对象**（H.2） |

⟹ **对定理而言这条 OPEN 可以放下**：无论 `N_{JJ′}` 读成哪一种**与 Allen 自洽**的意思，判据比值都不变。
⚠ **语义本身我仍不裁决**：Allen 从不写 `N_{JJ′}(0)` 这个符号，`√(N_JN_{J′})` 是**我从 Eq.(1)+(33) 推出来的**，
不是任一篇的原句。另记一处**同篇内张力**（只标不裁决）：Mazin–Golubov p.2 印的是 `Λ_ij = g_ij N_j`（**非对称**），
p.4 引 Allen 的那条按 H.3 是**对称**的 `√(N_iN_j)` 形 —— 两者都在 `f(i)g(j)` 类里，但不是同一个约定。

## H.4 采纳的四条 + 随行 OPEN

1. **消费点** [采纳，未独立复核]：判据只在库仑侧引理 μ-1 被消费，声子侧吃逐元非负/(Z-dom)
   ⟹ 我们凑的材料侧 4 个失效（FeSe/CaC₆/122/1111）**全在声子侧，从未喂进天花板链条**，
   它们证的是「判据有分辨力」而非「天花板失效」。本文件 §0 / A.7 / B.4 / C 各处的「不得传播成天花板失效」
   现在有了**正面理由**（此前只是「不知道所以别说」）。⚠ 该判定是派活方的 [analytic]，**我没有能力复核**（无 Bash）。
2. **`N_{JJ′}(0)` 两指标** [采纳，见 F.4.a 订正横幅]：我的单指标读法撤回；页图分辨率不足，文本层 span 是更好的仪器。
3. **全排斥条件** [采纳]：我第十轮标「未裁决」的那条，锐条件恰是 `Λ12Λ21 > Λ11Λ22` = 本仓判据的否定；
   印出的那行按 AM–GM 是**充分不必要**。⟹ 「全排斥仍能超导」与「判据失效」是同一不等式的两侧。
4. **读法 (A)** [采纳并加证]：见 H.1。

**★ 随行 OPEN（派活方点名要我记的）**：按 (A)，Golubov 侧 `V_ππ/V_σπ = 1.815` 与 Mazin–Antropov 转述的
Jepsen–Andersen `≈2.5` **差 27.4%**。两条著录不同（MG [64] = `Mazin, Andersen, Jepsen, Golubov, Dolgov, Kortus,
to be published`；Golubov [37] = `Kuz'menko et al …; Dolgov et al, to be published`）、作者表有重叠、**是否同一份计算未查**；
JA 原始 LMTO 数值**未取**（本轮检索也没找到印出 `3:2.5:1` 的第二处）。
⟹ **任何下游引用该裸 `V` 比值时必须随行标这条分歧。** 判据余量 5.5103×（印出四数直算）不受影响。

## H.5 本附录的未核实项

- Ummarino Eq.(10) 的出处 [22] 是**会议报告，取不到** ⟹ 「2.23/2.48 怎么算出来的」仍无第一手推导。
- Eq.(10) 与 Golubov 印出四数的一致只到**比值**层；绝对值不同（Ummarino 自己拟合 `μ(ω_c)=0.0326`、`ω_c=700 meV`），
  我**没有**做绝对值复核（`N^tot` 只能从其 Fig.3 读图，不承重）。
- H.3 的 `√(N_JN_{J′})` 是**我的推导**，Allen 原文无此符号；disjoint 表示的归一化 `ψ_J=√(N/N_J)1_{k∈J}`
  也是我由 Eq.(1) 反解的，Allen 未逐字写出。
- Allen §V 我只读了 pp.1422–1423（Eqs.26–39 所在）+ pp.1424–1426；**Sec. IV（Boltzmann）与 Table I 未读**。
- 本轮检索：**1 次 WebSearch**（结果是二手综述式摘要，**其中的 Choi/Louie `1.75:2.04:1.61:1.00` 我未取源、不引用**）、
  2 次 WebFetch（两份均「raw binary」被拒但已落盘，再用 `Read(pages=)` 读页图，**2/2 成功，累计 7/7**）。




