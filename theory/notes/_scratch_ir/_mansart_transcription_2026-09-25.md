# 一手转录：Mansart et al., Phys. Rev. B 82, 024513 (2010)

**来源**：`https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.82.024513/fulltext`（PDF，5 页，264846 字节，非 JBIG2 扫描——文本层可直接抽取，与页图逐字核对一致，无歧义字符）。
**两仪器**：`pdftotext -enc UTF-8 -layout` 全文 + `pdftoppm -r 200 -png` 页图（Read 工具读了 page-1、page-3）。两者在本篇涉及的全部承重字符上**完全一致**（无需仲裁）。
**产物**：`theory/notes/_scratch_ir/_mansart_prb82_024513.pdf`、`_mansart_prb82_024513.txt`、`_mansart_page-{1..5}.png`（均未入库，仍在 `_scratch_ir/`）。

---

## 表：四项速览

| # | 结果 | 数值/状态 |
|---|---|---|
| A | 全文 `λ⟨ω²⟩` | 摘要与正文汇总值 **64 meV²**；Table I 给 4 个逐拟合值 **66.1 / 69.6 / 60.9 / 59.2 meV²**（无显式误差棒） |
| B | 定义式与 Allen (1987) 是否同构 | **不同构 / 未直接写出矩积分**——该文经由三温度模型的能量弛豫率关系 `g = 6ħγ/(πkB) · λ⟨ω²⟩`（引 Allen, PRL 59, 1460 (1987)）**反解** `λ⟨ω²⟩`，全文没有写出 `λ⟨ωⁿ⟩ = 2∫₀^∞ dΩ [α²F(Ω)/Ω] Ωⁿ` 这个积分形式 |
| C | 同样品 Tc | **有，本文自己声明**：x=0.08，`Tc = 24 K`（引第 18 号文献表征该批样品；本文未重新给出电阻/磁化曲线，是引用同组前作的表征值，但样品明确是"from the same batches"） |
| D | `λ` 是否独立拟合 | **是除出来的，非独立拟合**——`λ` 由 `λ⟨ω²⟩` 除以选定的声子能标 `⟨ħω⟩` 反解得到（两种除法，见下） |

---

## A. `λ⟨ω²⟩` 数值

**原文逐字（摘要，页 1）**：
> "Our analysis indicates that the fast relaxation should be attributed to preferential scattering of the electrons with only a subset of the lattice-vibration modes with a second moment of the Eliashberg function λ⟨ω²⟩ ≅ 64 meV²."

**原文逐字（正文，页 3，紧接 g 的关系式之后）**：
> "By using the relation g = 6ħγ/(πkB) · λ⟨ω²⟩ [Ref. 28], we obtain for the second moment of the Eliashberg function λ⟨ω²⟩ ≅ 64 meV², in good agreement with Refs. 15 and 16."

**原文逐字（结论段，页 4）**：
> "On the basis of our analysis, we can estimate a value of λ⟨ω²⟩ ≅ 64 meV²"

**Table I（页 3，标题 "Maximum electronic temperatures, fraction of strongly coupled phonons, and electron-phonon coupling constants obtained by three temperature model simulations."）逐行**（列：Ti(K), F(mJ/cm²), Te^max(K), α, λ⟨ω²⟩(meV²), λ）：

| Ti (K) | F (mJ/cm²) | Te^max (K) | α | λ⟨ω²⟩ (meV²) | λ |
|---|---|---|---|---|---|
| 25 | 1.27 | 1100 | 0.4 | 66.1 | 0.125 |
| 10 | 1.27 | 1100 | 0.4 | 69.6 | 0.132 |
| — | 1.91 | 1350 | 0.35 | 60.9 | 0.115 |
| — | 3.18 | 1740 | 0.35 | 59.2 | 0.112 |

正文取自页图与文本层双重核对，四个 `λ⟨ω²⟩` 值在 Fig. 3 的四张分图标题里也各自复现一次（`λ<ω²> = 66.1/69.6/60.9/59.2 meV²`），逐位一致。

**误差棒**：全文（摘要、正文、Table I、Fig. 3）**均未给出 `λ⟨ω²⟩` 的误差棒**——只印中心值（`≅`/无 `±`）。唯一给出的区间量是 e-ph 耦合常数 `g` 本身：**"we obtain an e-ph coupling 6.8 < g < 8.0 mJ/K/s/m³"**（页 3），`λ⟨ω²⟩` 是由这个 `g` 区间的（隐含）中点通过上式反解出的单点值 64 meV²，而非独立测得的区间。

**出现位置**：Abstract；正文第 II 节末尾（页 3，紧邻 Eq.(4) 之前，g-关系式段落）；结论段（页 4）；Table I（页 3）；Fig. 3 四张子图标题（页 3）。**没有专属公式编号**——`λ⟨ω²⟩` 是通过反解页 3 那行未编号的 `g = 6ħγ/(πkB)·λ⟨ω²⟩` 得到，该关系式本身也没有独立的 Eq.(n) 标号（它写在 Eq.(2)/(3) 之间的正文散句里，不在编号方程之列）。

---

## B. `λ⟨ω²⟩` 定义式（前提门）

**原文逐字**（页 3，紧接 g 的取值区间之后）：
> "As a result of our analysis, we obtain an e-ph coupling 6.8 < g < 8.0 mJ/K/s/m³. By using the relation
> g = 6ħγ/(πkB) · λ⟨ω²⟩ [ref 28],
> we obtain for the second moment of the Eliashberg function λ⟨ω²⟩ ≅ 64 meV², in good agreement with Refs. 15 and 16."

引用的 ref 28 就是 **P. B. Allen, Phys. Rev. Lett. 59, 1460 (1987)**——与 CLAUDE 记忆库中天花板路线的锚文献同一篇。

**我的判断**：该文**没有在自己文中写出** `λ⟨ωⁿ⟩ = 2∫₀^∞ dΩ [α²F(Ω)/Ω] Ωⁿ` 这个矩积分定义式本身；它只写出并使用了由该定义**导出**的电子-声子能量弛豫率关系式 `g = 6ħγλ⟨ω²⟩/(πkB)`（这正是 Allen 1987 那篇里从矩积分定义推出的二级结果，用于双温/三温模型的弛豫率）。

这与已入库先验（`ultrafast-lamW2-corpus`/`ultrafast-pump-probe-measures-the-ceiling-denominator`）里其它超快泵浦-探测文献（如 LaH₁₀ 那篇）的做法**同构**——都是通过这个 `g`-关系式反解 `λ⟨ω²⟩`，而不是直接测/算 `α²F` 再积分。所以从"文献族内部一致性"角度看，本文与该语料的其余成员用的是同一条换算链、同一篇锚文献（Allen 1987），**前提门可判为同构**（同一个 `g ↔ λ⟨ω²⟩` 恒等式），但**严格意义上"逐字定义式"本身未见于本文**——若要复核这条恒等式的推导本身，需要回到 Allen PRL 59,1460 (1987) 原文，本文不提供独立推导，只提供引用。

**未见**：本文全文没有出现 `α²F` 这个符号、也没有写出对 `Ω` 的积分号。`⟨ω²⟩`／`λ⟨ω²⟩` 全程作为一个通过 `g` 反解出的整体量出现，从未被拆成 `λ` 与 `⟨ω²⟩` 两个独立测得的因子（`λ` 是之后另用声子能标反解出来的，见 D 节）。

---

## C. 同样品 Tc（G1 门）

**原文逐字（页 1，引言末段）**：
> "In this work, we studied the iron-pnictide Ba(Fe₁₋ₓCoₓ)₂As₂, at almost optimal doping, x = 0.08. This material is superconductor below Tc = 24 K and does not present any magnetic transition.¹⁸ Single-crystal samples were grown by the self-flux method¹⁸ and fully characterized prior to our measurements."

**原文逐字（结论段末，页 4，二次确认同一 Tc）**：
> "By taking ⟨ħω⟩ = 23 meV and μ* = 0, we obtain Tc ≅ 0.03 K, which is far below the actual Tc of about 24 K: therefore, the electron pair formation cannot be explained only by electron-lattice interactions."
> …
> "This low value cannot explain a Tc of 24 K within a standard BCS framework"

**参考文献 18**（页 4 文献列表）：
> "¹⁸ F. Rullier-Albenque, D. Colson, A. Forget, and H. Alloul, Phys. Rev. Lett. 103, 057001 (2009)."

**我的判断**：这是**同一批样品**的 Tc——本文明确写 "Single-crystal samples were grown by the self-flux method¹⁸"（同一篇 Ref. 18），且高温比热测量段落（页 2）写 "we performed high-resolution, low-temperature heat-capacity measurements in samples coming from the same batches"，说明超快泵浦-探测所用晶体与该组自己表征 Tc 的晶体是**同一批次（same batches）**，同一实验组（Rullier-Albenque, Colson, Forget 均在作者列里）。Tc=24 K 本身的**测量方法**（电阻零点/磁化率/比热中点）**未在本文直接给出**——本文把 Tc 当作已知量引用自己组的前作 Ref. 18，未重复给出判据。

⟹ **G1 门判定：通过，但带一个附注**——Tc 数值是本文自己声明、样品可溯源为"同一批次"，满足"同样品"要求；但 Tc 的**测定判据**（哪种探针、哪个阈值）本文没有自己重述，只能追溯到 Ref. 18（Rullier-Albenque et al., PRL 103, 057001 (2009)），若要核实判据细节需另取该文。这不构成"未印"，只是判据的出处隔了一层引用。

---

## D. `λ` 值与是否独立拟合（循环性审计）

**原文逐字（页 3–4，紧接 λ⟨ω²⟩≅64 meV² 之后）**：
> "In the estimation of ⟨ω²⟩, we also have to take into account the fact that some vibrational modes are more efficiently coupled than others. The most natural choice is to take the frequency of the fully symmetric A1g mode, which is coherently excited by our photoexcitation and consequently efficiently coupled. …Then we used the A1g energy, 23 meV,²⁰ as ⟨ħω⟩, which gives λ ≅ 0.12. We notice that, by taking the average of the whole phonon spectrum as measured by neutron diffusion,³¹ we obtain an even smaller but very close value of λ ≅ 0.10; the coupling given here is therefore a higher limit estimation of λ."

**Table I 最后一列**逐值：0.125 / 0.132 / 0.115 / 0.112（页 3，见上表）。

**我的判断（循环性审计）**：`λ` **不是独立拟合出来的**，而是把上一步反解出的 `λ⟨ω²⟩`（=64 meV² 汇总值，或 Table I 各行逐值）**除以**一个另外选定的声子能标 `⟨ħω⟩` 得到的——文中明写了两种除法，给出两个不同的 `λ`：
1. 用 A1g 声干模能量 23 meV（引 Ref. 20，即同一作者组前作 PRB 80, 172504 (2009)）：`λ ≅ 0.12`
2. 用中子衍射测得的全声子谱平均：`λ ≅ 0.10`（本文称这是"更小但很接近"的值，作者明确说"the coupling given here is therefore a higher limit estimation of λ"，即 0.12 是取的**上限**估计）

Table I 逐行的 λ（0.125/0.132/0.115/0.112）是用 `λ⟨ω²⟩` 逐行值除以同一个 23 meV 得到（`λ = λ⟨ω²⟩ / ⟨ħω⟩²`，量纲核对：`66.1 meV² / (23 meV)² = 0.1250`，精确复现 Table I 第一行 0.125；`69.6/23² = 0.1316≈0.132`；`60.9/23²=0.1151≈0.115`；`59.2/23²=0.1119≈0.112`——四行全部用 `λ = λ⟨ω²⟩/529` 精确复现，**零剩余自由度**）。

⟹ 这与仓库记忆 `hydride-m2-hunt` / `lah10-omega2` 系列已识别的模式**同形**：`λ` 不是第二个独立观测量，而是 `λ⟨ω²⟩`（唯一真正测得的量）除以一个从别处（本文自己前作 Ref.20 的拉曼/相干声子频率，或 Ref.31 的中子谱）取来的能标平方。**若要把本文的 `λ` 列当独立数据点使用，需先确认 `⟨ħω⟩=23 meV` 本身不是从别的循环链条里来的**——本文明写它取自"the A1g energy, 23 meV"（引 Ref. 20，本组前作 PRB 80,172504(2009)，报的是相干声子振荡频率，属独立拉曼/泵浦探测测量，不是从本文的 λ⟨ω²⟩ 反推），所以这条链没有二次循环，但 `λ` 确实是**除出来的**、不是本文独立拟合的第二个自由参数。

---

## 门禁结论小结（供主会话判词参考，不代替判词）

- **B 门（定义同构）**：本文与 Allen 1987 的联系通过引用成立，但本文未展示矩积分定义式本身，只展示并使用 `g=6ħγλ⟨ω²⟩/(πkB)` 这条导出关系。
- **C 门（同样品 Tc）**：通过——Tc=24K 是本文自己写的，样品明确"from the same batches"，但 Tc 的测定判据本身要追溯到 Ref. 18。
- **D（循环性）**：`λ` 是除出来的（`λ⟨ω²⟩` 除以 `⟨ħω⟩²`，两种能标各给一个值，0.12 为上限估计），不是独立拟合值。`λ⟨ω²⟩` 本身也是通过 `g` 关系式反解、无直接误差棒。
