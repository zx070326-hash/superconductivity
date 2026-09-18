# 文献查新：γ-model + 静态库仑排斥（2026-09-16）

**任务**：判定三条命题的新颖性。立场是**对抗性的**——目标是找到先行工作，不是给安慰。
**检索范围**：arXiv（cond-mat.supr-con / str-el）、Google/Bing 网页检索、若干 PDF 全文直读。
**证据等级约定**：
- **L1** = 本轮亲自读到一手正文（PDF 逐页），逐字引用可信。
- **L2** = 只读到摘要 / 只有检索引擎转述，未见一手正文。
- **推断** = 我自己的综合，不是任何来源的直接陈述。

---

## 0. 结论速览

| 命题 | 判决 | 最近的先行工作 |
|---|---|---|
| 1. 库仑压制对 N_c **非单调**（有限 N_c 处取峰） | **NOT FOUND** | 无。最近的是 Elezaby–Abanov 的**单调**截断阶梯（不同对象） |
| 2. 秩一久期方程 → AM/BT 型屏蔽律 μ_eff = μ/(1+μ ln N_c/(2𝔤)) | **PARTIALLY KNOWN** | Pimenov–Chubukov 2112.06273 p.6 已写出 f/[1+λ f log Λ] → 1/(λ log Λ) |
| 3. μ 在配对 QCP 边缘无关；任意强排斥杀不死量子临界配对 | **PARTIALLY KNOWN（概念层面已发表）** | Raghu–Torroba–Wang 1507.06652 App. B p.27 逐字："there is no choice of λ₀ that can kill the superconducting instabilities" |

---

## 1. 命题 1：𝔤(γ) 的库仑压制对 N_c 非单调 — **NOT FOUND**

**没有找到**任何文献报告"配对本征值对库仑 Matsubara 截断的压制是非单调的"。

反复用不同措辞检索（"non-monotonic" + "pairing eigenvalue" + "Matsubara cutoff"；"optimal cutoff" + Tc + Morel-Anderson；"mu*(omega_c)" + non-monotonic）均未命中。

**最接近但不是同一件事的三处：**

1. **Elezaby & Abanov, arXiv:2512.20009**（PRL 2026，"Superconductivity Near a Quantum Critical Point: Bounds on the Transition Temperature in the γ-Model"）——**L2**（WebFetch 全文检索式问答，未逐页读 PDF）。他们的截断序列是**严格单调的**：`τ_c^(1) < τ_c^(2) < ⋯ < τ_c^(∞) = τ_c`，Fig. 3 给 N = 100/400/1000 从下方收敛。
   - **关键区分**：他们的 N 是**声子算子的变分截断**（= KAY Thm 4 的阶梯），不是库仑作用域 N_c。两个 N 的物理含义完全不同，不构成撞车。
   - 同一次全文检索确认：该文中 "Coulomb"/"pseudopotential"/"repulsion"/"repulsive" **零出现**。

2. **Pimenov & Chubukov, arXiv:2112.06273 p.6**（**L1**，逐页读）："Third, for any finite Λ, the repusive part of the interaction gets reduced, but cannot be completely eliminated." —— 这句只支持"随 Λ **单调**减弱但不归零"，**没有**非单调的意思。

3. 常规文献里 μ*(ω_c) = μ/[1+μ ln(ω_pl/ω_c)] 对 ω_c **单调递减**（**L2**，检索转述的标准结果，与本仓库 memory 中 Marsiglio eq.(14) 一致）。

**推断（我的，非来源陈述）**：命题 1 的非单调性只能来自"排斥算符与顶本征向量的重叠随 N_c 先增后被对数屏蔽压过"这一竞争，而现有文献只写出了屏蔽那一半。这解释了为什么它没被报道过。

⚠ **诚实边界**：我没有逐个 grep γ-model 系列 I–VI（2004.13220 / 2006.02968 / 2007.14540 / 2009.10911 / 2102.04444 / 2107.14340）的全文。本轮对"γ-model 系列不含库仑"的判断，证据来自 (a) 2026 年权威综述全篇无库仑、(b) 2512.20009 全文无库仑、(c) 多次针对性检索无命中。这是**强间接证据，不是逐篇核实**。

---

## 2. 命题 2：秩一久期方程 → Matsubara 窗口设定的对数 — **PARTIALLY KNOWN**

### 2.1 已知的部分（撞车点）

**Pimenov & Chubukov, arXiv:2112.06273v2 = npj Quantum Materials 7, 45 (2022)**，p.6，**L1 逐字**：

> "Second, optimizing ω₀ in the limit log(Λ) ≫ 1, one obtains that a repulsive f effectively gets reduced to f/[1 + λ f log (Λ)] → 1/λ log(Λ). This gives rise to the McMillan formula for T_c ∝ e^{−1/(λ−μ*)}, in which μ* ≈ 1/log Λ is the contribution from the repulsion. Third, for any finite Λ, the repusive part of the interaction gets reduced, but cannot be completely eliminated. As a result, superconductivity exists at f smaller than some critical f_c > 1."

同文 p.3 的引言版本（**L1 逐字**）：

> "For f > 1, V(Ω_m) > 0 for all frequencies, yet for 1 < f < f_c, superconductivity emerges below a finite T_c, which contains f in the combination f/[1 + const. × f log(E_F/Ω_1)]. For a given f and large enough log(E_F/Ω_1), the Coulomb repulsion becomes logarithmically small, and one recovers the McMillan formula for T_c."

**⟹ 撞车的具体内容**：
- 「μ_eff = μ/(1 + c·μ·log(截断)) 」这个**函数形式**已发表；
- 「对数前系数 c 由**耦合本身**给出（他们是 λ），不是裸的 1」——已发表；
- 「f → ∞ 时 μ_eff 饱和到与 f 无关的值 1/(λ log Λ)」——已发表。**这正是"T_c 渐近不依赖裸库仑强度"的机制**（见 §3）。

他们的模型（p.5, eq.(3)，**L1**）：`Ṽ(Ω) = λ(f − 1/(1+(Ω/Ω₁)²))`，Ω₁ ~ Debye 能标，Λ ~ E_F 的 UV 截断，作用于全部频率（无 θ 函数窗口）。

**溯源链**：Pimenov–Chubukov 把这条线归给 Refs. [8–12] = Rietschel & Sham PRB 28, 5100 (1983)；Ruhman & Lee PRB 94, 224515 (2016) 与 PRB 96, 235107 (2017)；Wölfle & Balatsky PRB 98, 104505 (2018)；Chubukov, Prokof'ev & Svistunov PRB 100, 064513 (2019)（arXiv:1905.03732，"Implicit renormalization approach to the problem of Cooper instability"）。这五篇我**只读到题录，未读正文（L2）**；若要写"我们的形式与前人的差别"，**最后那篇最该亲自核**——它的方法论（不用原始 gap 方程最大本征值、改用隐式重整化的替代本征值问题）与秩一久期方程最可能同型。

### 2.2 仍然未找到的部分（可主张为新）

- 把对数认定为 **Matsubara 窗口 ln N_c = ln(Λ/2πT)**，而不是带结构比 E_F/ω_D —— **NOT FOUND**。这在 γ-model 里才成立，因为那里根本没有玻色能标。
- 系数被认定为 **1/(2𝔤(γ))**，即由配对本征值本身给出 —— **NOT FOUND**。Pimenov–Chubukov 的 "const." 从未定值；p.6 版本的系数是 λ（他们模型的无量纲耦合）。
- 用**秩一算符的久期方程**（Matsubara 指标上的投影算符 |n|<N_c）推导 —— 在严格数学侧和物理侧都 **NOT FOUND**（见 §4）。

**推断（我的）**：referee 会说"形式是老的"。可主张的新颖性应集中在"**什么设定了对数**"与"**系数是谱量而非带结构量**"，而不是屏蔽律本身。

---

## 3. 命题 3：μ 在配对 QCP 边缘无关 — **PARTIALLY KNOWN（概念层面已发表）**

### 3.1 最强先行工作：Raghu, Torroba & Wang（2015）

**S. Raghu, G. Torroba, H. Wang, "Metallic quantum critical points with finite BCS couplings", arXiv:1507.06652, Phys. Rev. B 92, 205104 (2015)** —— **L1**，逐页读了 pp. 1–8, 15–24, 26–28。

符号约定（p.5，L1 逐字）："The sign in L_BCS is such that λ₀ > 0 corresponds to an attractive interaction." ⟹ **λ₀ < 0 就是排斥**。

BCS β 函数（eq. III.11，p.15，L1）：
```
μ dλ/dμ = −2π²α + 2γλ − λ²/(2π²N)
```

**★ 决定性的一句（Appendix B, p.27，L1 逐字）**：

> "We see that once we cross below αN = 4, the phase transition in λ₀ reduces to a step discontinuity across λ_disc = π²αN: **there is no choice of λ₀ that can kill the superconducting instabilities.**"

配套的 eq. (B.3)（p.26，L1）：`μ_sc^- = M exp[−2π/√(4α/N − α²)] · f(λ₀)`，其中 f(λ₀) 是有界前因子 ⟹ **超导标度的指数与 λ₀ 无关，λ₀ 只进前因子**。

**⟹ 命题 3 的物理内容（任意强排斥都杀不死量子临界配对；标度渐近不依赖裸排斥强度）在 RTW 2015 已经以 RG 语言发表。**

**与我们的差距（gap）**：
1. RTW 用的是 d = 3−ε 的 ε 展开 + 大 N + 一圈 β 函数，**不是** γ-model 的 Matsubara Eliashberg 方程，也不是有限温 T_c。
2. RTW 全文**从不写** "Coulomb"、"μ*"、"pseudopotential"、"Morel"、"Tolmachev"——他们的 λ₀ 是抽象 4-费米 BCS 耦合，没有屏蔽律，没有 μ_eff 公式。
3. RTW 的陈述只在 αN < 4 的参数窗内成立（αN > 4 时存在 λ_+ 不稳定不动点，SC 被 NFL 杀死），**不是无条件的**。
4. 他们**没有**说"no choice of λ₀"包括排斥这一层意思——原文是中性的"no choice"，覆盖了负 λ₀，但作者没有把它当作"库仑杀不死配对"来宣传。

### 3.2 次强：Anderson–Morel/Bogoliubov–Tolmachev 的经典内容

「窗口对数 → ∞ ⟹ 任意裸排斥被驯服」本身是教科书。两处文献确认：
- Pimenov–Chubukov p.3（L1，见 §2.1 引文）："For a given f and large enough log(E_F/Ω_1), the Coulomb repulsion becomes logarithmically small."
- RTW p.3（L1 逐字）："This is in sharp contrast to the behavior of Fermi liquids, where the BCS coupling flows to zero for repulsive forces, or grows indefinitely, leading to a BCS instability if the couplings are attractive." ——即"排斥 BCS 耦合边缘无关"在 FL 里是**常识**（Shankar RMP 66, 129 (1994)；**L2**，本轮未读）。

**⟹ "marginally irrelevant" 这个措辞本身没有新颖性可言。** 新颖性只能落在"**在 γ-model / 无 Debye 窗口的 QCP 上，窗口自动随 T→0 趋于无穷**"这一结构性理由。

### 3.3 ★ 反面证据：相邻模型里排斥**能**杀死超导

这是审稿人最可能拿来打的两枪，必须正面处理：

- **Pimenov & Chubukov 2112.06273**（L1）：固定有限 Λ 时存在**临界排斥 f_c**，f > f_c 超导消失；全文主题就是这个 T=0 量子相变。
- **Phan & Chubukov, arXiv:2112.03490 = PRB 105, 064518 (2022)**，"The Effect of Repulsion on Superconductivity at Low Density"（**L2**，只读摘要+检索转述）：低密度（μ ≪ ω_L，无对数窗口）时同样存在阈值 f*，且 T_c ∝ (f*−f) 的幂律而非指数。
- **Abanov & Chubukov 综述 2606.21731**（L1，读了 pp. 1–4, 17–19）：γ-model 自己的"配对可被杀死"旋钮是 **1/N**（pp 与 ph 通道耦合之比），存在阈值 N_cr："We showed that at a given γ, superconductivity develops only when N becomes smaller than the threshold value N_cr. At larger N, the ground state remains a non-FL."

**⟹ 命题 3 必须明确说清：我们说的"杀不死"是对**带 N_c 窗口的静态排斥**而言，与 γ-model 已知的 N_cr 阈值是**不同旋钮**，不矛盾。** 不做这个区分会被直接拒。

### 3.4 概念先例（低置信度，**未核一手**）

高密度 QCD 色超导里有结构同型的已知说法：长程未屏蔽磁胶子交换给 Δ ~ exp(−c/g)，而短程（Debye 屏蔽的电胶子、瞬子诱导的四费米项）**只进前因子、不改指数**（Son, PRD 59, 094019 (1999) 一系）。这与命题 3 的"排斥只改前因子"同构。
⚠ **本轮只有检索引擎的综合转述（L2/推断混合），未读任何一手正文。若要引用必须先核。**

---

## 4. 严格数学侧（Frank / Hainzl / Seiringer / Solovej）——秩一久期论证是否标准？

**NOT FOUND**。检索（**L2**，只读检索转述与题录）显示：
- Hainzl, Hamza, Seiringer & Solovej, "The BCS functional for general pair interactions", CMP 281, 349 (2008)：用 Birman–Schwinger 原理，把势拆成 V₋V₊ = 0 两部分，结论形如"排斥部分 V₊ **不太大**时 T_c > 0"。
- Frank, Hainzl, Naboko & Seiringer, J. Geom. Anal. 17, 559 (2007)；Frank–Hainzl–Seiringer–Solovej, JAMS 25, 667 (2012)、CMP 342, 189 (2016)；Hainzl–Seiringer, JMP 57, 021101 (2016)。

**没有找到**以"Matsubara 指标上的秩一排斥 + 久期方程"为工具的严格工作。该纲领**全是瞬时核**（与本仓库 memory 中 [[lit_kay_rigorous_tc_bounds_p4]] 的既有结论一致：推迟核严格 T_c 界是处女地）。

---

## 5. 拒绝引入库仑的脚注（逐字，**L1 本轮重新核实**）

**Kiessling, Altshuler & Yuzbashyan, "Bounds on T_c in the Eliashberg theory of Superconductivity. I: The γ-model", arXiv:2409.00533v2 = J. Stat. Phys. 192, 69 (2025)，p. 9，§1.3 "Mathematical critique and goals" 的脚注 4**：

> "⁴In Eq. (7) of [BR] the kernel λ(n − m) is shifted by the addition of an N-dependent "reduced Coulomb repulsion parameter" μ*(N). We emphasize that no such Coulomb term can be extracted from the Fröhlich hamiltonian (1), and so it has no place in Eliashberg theory as derived from (1). It is "smuggelt" into the truncated equations using some hand-waving physical reasoning, basically for the purpose of obtaining better fits to experimental data. Note that no non-zero constant can be added to λ(n − m) in Eliashberg theory because the summations over all Matsubara frequencies would diverge, then."

（[BR] = Bergmann & Rainer, Z. Phys. 263, 59 (1973)；原文拼写 "smuggelt" 与 "repusive" 等排印保留原样。）

**回复 / Comment 检索结果：NOT FOUND。** 未找到任何针对该脚注的 Comment/Reply。最接近的后续是 Elezaby & Abanov arXiv:2512.20009（收紧上界、独立复核下界），该文**同样一字不提库仑**。

**⟹ 对我们有利的一点**：脚注的第二句（"no non-zero constant can be added … the summations would diverge"）恰恰是**带 N_c 窗口的秩一构造绕开的**——加了截断就不发散。这是可以正面回应而非回避的切口。

---

## 6. 五篇最相关文献（每篇一行）

1. **Pimenov & Chubukov, arXiv:2112.06273 = npj Quantum Mater. 7, 45 (2022)** — 常数排斥 f + 动力学吸引，写出 f/[1+λf log Λ] → 1/(λ log Λ)，并给出 f > f_c 时超导消失的 T=0 量子相变；**命题 2 的主要撞车源**。
2. **Raghu, Torroba & Wang, arXiv:1507.06652 = PRB 92, 205104 (2015)** — QCP 上 BCS 耦合的 RG；Appendix B 明写 αN<4 时"no choice of λ₀ can kill the superconducting instabilities"；**命题 3 的主要撞车源**。
3. **Abanov & Chubukov, arXiv:2606.21731 = Annu. Rev. Condens. Matter Phys. 18 (2026)** — 2026 年 γ-model 权威综述；全篇**零库仑/μ***，Conclusions/Summary Points/Future Issues 也不列库仑为开放问题；**最强的"空白"锚点**。
4. **Kiessling, Altshuler & Yuzbashyan, arXiv:2409.00533 = J. Stat. Phys. 192, 69 (2025)** — γ-model 严格 T_c 界；脚注 4 明确拒绝库仑项（见 §5）。
5. **Phan & Chubukov, arXiv:2112.03490 = PRB 105, 064518 (2022)** — 低密度（无对数窗口）时排斥阈值 f* 与 T_c ∝ (f*−f) 幂律；**命题 3 的反面对照**。

（第六篇候选：**Elezaby & Abanov, arXiv:2512.20009, PRL (2026)** — γ-model 本征值上下界的最新收紧，单调截断阶梯，全文无库仑。）

---

## 7. 明确标注：哪些是我的推断而非读到的

- §1 关于"非单调性来自重叠 vs 屏蔽的竞争"的机制解释 —— **我的推断**，无来源。
- §3.1 "RTW 的 'no choice of λ₀' 覆盖了负（排斥）λ₀" —— **我的推断**（原文只写 "no choice of λ₀"，作者未点明排斥；但 p.5 的符号约定使这个推断是直接的）。
- §3.4 色超导先例 —— **检索引擎综合 + 我的类比**，未核一手。
- §2.2 "referee 会说形式是老的" —— **我的判断**，非文献陈述。
- "γ-model 系列 I–VI 不含库仑" —— **强间接证据下的推断**，非逐篇核实（见 §1 诚实边界）。

---

## 8. 建议的后续检索方向（按优先级）

1. **Chubukov, Prokof'ev & Svistunov, PRB 100, 064513 (2019) / arXiv:1905.03732** 全文 —— 最可能与"秩一久期方程"方法论同型；本轮只读题录。**最高优先级。**
2. **逐篇 grep γ-model 系列 I–VI 全文**（2004.13220 / 2006.02968 / 2007.14540 / 2009.10911 / 2102.04444 / 2107.14340）查 "Coulomb"/"repulsion"/"μ*"，把 §1 的间接证据升级为逐篇核实。
3. **Rietschel & Sham PRB 28, 5100 (1983)** —— Pimenov–Chubukov 归因链的源头，可能是 μ_eff 形式的最早出处；决定命题 2 的"已知"边界该画在哪一年。
4. **Ruhman & Lee PRB 94, 224515 (2016) / PRB 96, 235107 (2017)；Wölfle & Balatsky PRB 98, 104505 (2018)** —— 低 E_F / 铁电 QCP + 排斥，检查是否已有"窗口由别的东西设定"的陈述。
5. **Son PRD 59, 094019 (1999)** 一手 —— 若要用色超导作概念先例，必须核。
6. **Metlitski, Mross, Sachdev & Senthil, PRB 91, 115111 (2015)**（arXiv:1403.3694）全文 —— 本轮三次尝试取正文均失败（arXiv abs 页只有摘要；Harvard DASH 返回 HTTP 405）。他们的 RG 里若有排斥支路，会是命题 3 的第二个撞车源。**这是本轮最大的未闭合缺口。**

---

## 9. 检索限制（诚实声明）

- 本 agent **无 Bash**，未运行任何脚本，未改动仓库其它文件，**未 commit**。
- MMSS 1403.3694 正文三次取件失败（见 §8.6）。
- §2.1 溯源链的五篇（Rietschel–Sham / Ruhman–Lee ×2 / Wölfle–Balatsky / Chubukov–Prokof'ev–Svistunov）**均只读题录**。
- γ-model 系列 I–VI **未逐篇 grep**。
- Google Scholar 未直接使用（WebSearch 后端不支持 boolean 语法，多次查询被降级为关键词匹配）。
