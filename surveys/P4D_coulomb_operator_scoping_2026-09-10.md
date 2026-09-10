# P4D · 库仑赝势 μ* 进入 KAY 紧自伴算子框架：一手核实

调研 agent：researcher　日期：2026-09-10　worktree：`.claude/worktrees/p3-numerics`

**证据等级**：**A** = 原文逐字句/公式（本轮亲自读取 PDF 页面）；**B** = 我方组装/自算/推断（原文成分为 A 级但组合是我做的）；**C** = 二手转述。

**本轮读取的一手材料**（全部 WebFetch 落盘 PDF + `Read(pages=)` 逐页图像读取）：

| 文献 | 版本 | 读取页 |
|---|---|---|
| KAY I，arXiv:2409.00533v2 | 47pp | pp.4–9（§1.2–1.3 + 脚注 4）、pp.44–47（Appendix B + References） |
| KAY II，arXiv:2409.00532v3 | 42pp | pp.10–12（Thm 5/6, Cor 1/2）、pp.15–20（§4 全部算子定义 eqs.34–46）、pp.23–25（eq.51 + eqs.59–63） |
| KAY III，arXiv:2409.02121v2 | 36pp | pp.27–30（§4.1 Summary + §4.2 Outlook + Appendix A 开头） |
| Marsiglio & Carbotte, "Electron–Phonon Superconductivity"，cond-mat/0106143v1 | 123pp | pp.12–24、pp.27–35 |
| Marsiglio, "Eliashberg theory: a short review"，arXiv:1911.05065 (= Ann. Phys. 417, 168102 (2020) = KAY 的 [Ma]) | | pp.5–12 |

另外三篇 KAY 全文用 `arxiv.org/html/…` 做过 "Coulomb / pseudopotential / μ* / Morel / Anderson / repulsion / screened" 的全文关键词扫描（见 §1）。

---

## §0　一句话结论（先看这条）

**KAY 不是"为简化而设 μ*=0"，而是明确认为 μ* 项在 Eliashberg 理论里根本不合法。** KAY I 脚注 4（p.9）逐字批评 Bergmann–Rainer/Allen–Dynes 的 μ*(N) 是被 "smuggelt"（原文拼写）进截断方程的，并给出一条**恰好就是我们 Q3 技术难点的**理由："no non-zero constant can be added to λ(n−m) in Eliashberg theory because the summations over all Matsubara frequencies would diverge."

而这条"发散"警告在算子语言里可以被精确量化：库仑算子在 KAY 的 ℓ²(ℕ₀) 归一化下是 **𝔎_c = 2|u⟩⟨u| − D_{N_c}**（u_n = 1/√(2n+1)·1_{n<N_c}），它在每个固定 T>0 都是**有限秩（rank = N_c）因而紧、且 Hilbert–Schmidt**——紧性完全没有被破坏——但它的**算子范数 ‖𝔎_c‖ ≈ ln N_c + 1.9635 随截断对数发散，不一致有界**。团队负责人的两半猜测**都反了**：不是"有界但非紧"，而是**"紧（甚至有限秩）但不一致有界"**。

而且这个对数发散**不是病灶，它就是 Morel–Anderson/Tolmachev 对数本身**——它恰恰解释了为什么 μ*(ω_c) 必须随 ω_c 跑动。更巧的是：这个 2|u⟩⟨u| − (对角) 结构与 **KAY II eq.(51)（T→0 极限矩阵 −ℑ^(N) + 2Ξ*_N⊗Ξ*_N，ξ*_n = 1/√(2n+1)）逐位同型**，KAY 自己已经把这类算子的谱算干净了（Theorem 5, eq.24）。

---

## §1　Q1 · KAY 三篇里的库仑内容

### 1-1　全文扫描结果（A 级，负面结论）

- **KAY II (2409.00532v3)**：全文**零处**出现 Coulomb / pseudopotential / μ* / Morel / Anderson / repulsion / repulsive。
- **KAY III (2409.02121v2)**：全文**零处**。§4.2 Outlook（p.29，本轮逐页读取）也未提库仑，只说："There are further issues concerning the linearized Eliashberg gap equations that merit clarification or vindication, but these fall outside the thrust of our three papers about bounds on T_c for various realizations of Eliashberg theory, and will be addressed elsewhere." 下一步明写是**非线性方程 + Banach 空间**，不是库仑。
- **KAY I (2409.00533v2)**：**唯一一处**，在 §1.3 "Mathematical critique and goals" 的**脚注 4**（p.9）。

### 1-2　KAY I 脚注 4，逐字（A 级 ★核心证据★）

> ⁴In Eq. (7) of [BR] the kernel λ(n−m) is shifted by the addition of an N-dependent "reduced Coulomb repulsion parameter" μ*(N). We emphasize that no such Coulomb term can be extracted from the Fröhlich hamiltonian (1), and so it has no place in Eliashberg theory as derived from (1). It is "smuggelt" into the truncated equations using some hand-waving physical reasoning, basically for the purpose of obtaining better fits to experimental data. Note that no non-zero constant can be added to λ(n−m) in Eliashberg theory because the summations over all Matsubara frequencies would diverge, then.
>
> — KAY I, arXiv:2409.00533v2, p.9, footnote 4

同页正文（§1.3, p.9）配套引用 [BR] 原话（A 级）：

> ... in [BR] one also reads that "The frequency summation in the [linearized]³ gap equation (2) should be done up to a frequency of the order of the band width. For technical reasons one only does the summation up to a cut-off frequency somewhat larger than the maximum phonon frequency." This sounds very much like a BCS-type frequency cutoff Ω_cut having re-entered through a backdoor. Note that [AD] start directly with Eq. (7) of [BR], also their Eq. (7), a judiciously truncated version of the linearized Eliashberg gap equations.⁴

**[BR] 的完整出处（A 级，KAY I p.46 参考文献表逐字）**：
> [BR] Bergmann, G., and Rainer, D., *The sensitivity of the transition temperature to changes in α²F(ω)*, Z. Phys. **263**, 59–68 (1973).

（顺带，同表 p.46 逐字确认 **[Co] Combescot, R., *Strong-coupling limit of Eliashberg theory*, Phys. Rev. B **51**, 11625 (1995)** —— 此前只有卷页号，标题现为 A 级；[AD] = Phys. Rev. B **12**, 905–922 (1975)；[AM] = Allen & Mitrovic, Solid State Phys. **37**, 1–92 (1982)，注意是 **1982** 不是常引的 1983。）

### 1-3　直答"他们为什么设 μ*=0"、"有没有说若加入会如何"

- **"为什么"**：**不是**"为简化设为零"。KAY 的立场是**规范性的**：Fröhlich 哈密顿量（他们的 eq.(1)）里推不出库仑项，所以它 "has no place in Eliashberg theory as derived from (1)"。三篇里**没有任何一句** "we set μ*=0 for simplicity" 之类的软化表述。
- **"若加入会如何"**：**没有**。三篇全文找不到任何 "if one includes / were one to add …" 型的前瞻句。唯一接近的是脚注 4 末句那个**技术性否定**（加常数 ⇒ Matsubara 求和发散）。

> **⚠ 项目风险（我方判断，B 级）**：P4D 的整个立项前提（"把 KAY 上界推广到含 μ*"）与 KAY 本人对 μ* 的公开立场**直接对立**。这不影响数学正确性，但影响**表述方式**：论文里不能写成"沿 KAY 的思路自然推广"，而应写成"我们说明：即便接受超出 Fröhlich 模型的唯象 μ* 项（KAY [I, fn.4] 对此持保留意见），KAY 的紧算子/变分框架依然可以承载它，且 KAY 脚注 4 指出的发散正是 Morel–Anderson 对数"。把 KAY 的反对意见**正面引用并回应**，比回避它安全得多。

---

## §2　Q2 · 含库仑的标准各向同性线性化 Eliashberg 方程（Tc 处）

### 2-1　(a)+(b)　最干净的 canonical 形式：Marsiglio 2020（A 级，逐字）

来源：F. Marsiglio, "Eliashberg theory: a short review", arXiv:1911.05065（= Ann. Phys. **417**, 168102 (2020)；**正是 KAY I 参考表里的 [Ma]**），**p.5**，逐字：

> Z(iω_m) = 1 + (πT_c/ω_m) { λ + 2 Σ_{n=1}^{m−1} λ(iν_n) } .　　(22)
>
> Z(iω_m)Δ(iω_m) = πT_c Σ_{m′=−∞}^{+∞} [ λ(iω_m − iω_{m′}) − u*(ω_c) θ(ω_c − |ω_{m′}|) ] Δ(iω_{m′})/|ω_{m′}| .　　(23)

其非线性（T<T_c）母式是同文 **eq.(16)**（p.5），结构相同，只是分母换成 √(ω²_{m′}+Δ²(ω_{m′}))。同页定义（A 级）：

> where u*(ω_c) ≡ g(ε_F)U*(ω_c). One should note that U*(ω_c) < U, physically corresponding to the fact that retardation effects allow two electrons to exchange a phonon while not being at the same place at the same time. This means they do not feel the full direct Coulomb interaction with one another.

**⇒ 团队负责人的判断正确并已 A 级证实：μ* 只进能隙通道，Z 方程 eq.(22) 里完全没有 u\*。**

### 2-2　交叉核验：Marsiglio–Carbotte 综述（A 级，逐字）

cond-mat/0106143v1, **p.18**, eqs.(54)(56)（Fermi 面平均后的 Eliashberg 方程组）：

> Z_m = 1 + πT Σ_{m′} λ(iω_m − iω_{m′}) (ω_{m′}/ω_m)Z_{m′}/√(ω²_{m′}Z²_{m′}+φ²_{m′}) · A₀(m′)　　(54)
>
> φ_m = πT Σ_{m′} ( λ(iω_m − iω_{m′}) − N(μ)V_coul ) φ_{m′}/√(ω²_{m′}Z²_{m′}+φ²_{m′}) · A₀(m′)　　(56)

同文 **p.19** 第一段逐字（★这句决定截断的确切位置★）：

> If the electron density of states is assumed to be constant, then, with the additional approximation of infinite bandwidth, A₀(m′) ≡ 1 (actually a cutoff, θ(ω_c− |ω_{m′}|), is required in Eq. (56)), and A₁(m′) ≡ 0.

同文 **p.30**（§4 开篇）逐字，确认两个输入参数：

> In the conventional theory there are two input "parameters": a function of frequency, α²F(ν), … and μ*(ω_c), a number which summarizes the (reduced) Coulomb repulsion experienced by a Cooper electron pair.

**两源一致：库仑只出现在 φ（能隙/配对）通道，Z 通道无 μ\*。**

### 2-3　⚠ 一个必须由主会话拍板的分歧：截断是**单边**还是**双边**

这是本轮发现的、对"自伴算子"这件事有直接后果的技术分歧：

| 写法 | 出处 | 截断落在哪 | 算子是否对称 |
|---|---|---|---|
| θ(ω_c − \|ω_{m′}\|) **只在求和指标 m′ 上** | Marsiglio eq.(23)（A 级）；M-C eq.(56) 的 A₀(m′)（A 级） | 单边 | **否**（形如 \|v⟩⟨u\|，v 全支撑、u 截断） |
| 整个方程（含 n 指标）一起截到 N，且 μ* 随 N 变 | Bergmann–Rainer eq.(7) / Allen–Dynes eq.(7)，**由 KAY I 脚注 4 转述为 "N-dependent μ\*(N)"**（A 级转述，BR/AD 原文未取） | 双边 | **是**（对称秩一 + 对角） |

**我方建议（B 级）**：KAY 的整套机器（二次型 Q(Ξ)、变分原理 Λ=1/𝔨、Rayleigh 上界、Perron–Frobenius/Krein–Rutman）**要求自伴**，所以 P4D 必须采用**双边对称**版本 θ(ω_c−|ω_n|)θ(ω_c−|ω_m|)。这在物理上等价于"整条方程截在 ω_c"，正是 BR/AD 的做法，也正是 KAY 脚注 4 攻击的那个做法。文中必须**显式声明**采用双边约定，并说明与 Marsiglio eq.(23) 单边写法的差别只在 |ω_n|>ω_c 的分量（那里 Δ_n 已经很小），否则 reviewer 会抓这一点。

### 2-4　(c)　Tolmachev–Morel–Anderson 对数重整化

**A 级（本轮取到的唯一逐字源）**：Marsiglio arXiv:1911.05065 **eq.(14)**：

> U*(ω_c) ≡ U / ( 1 + U ln(W/2ω_c) )

即**对数宗量是"半带宽 / 截断频率"= W/2ω_c**，不是团队负责人问句里写的 E_F/ω_c。

- 两者在"半满平带、E_F ≈ W/2"约定下重合，但**引用时不能混写**。若 P4D 正文要写 μ*(ω_c) = μ_c/(1+μ_c ln(E_F/ω_c))，应加一句"（等价地，以半带宽标度，μ* = μ/(1+μ ln(W/2ω_c))，Marsiglio 2020 eq.(14)）"。
- 注意 Marsiglio 的 U 是**无量纲化前**的量；无量纲赝势是 u*(ω_c) ≡ g(ε_F)U*(ω_c)（同文 p.5）。
- **Morel–Anderson 原文（Phys. Rev. 125, 1263 (1962)）与 Allen–Mitrović SSP 37 (1982) 本轮均未取得**（前者无 arXiv，后者是书章）。**Allen–Dynes PRB 12, 905 原文本轮同样未取得**（APS 403，且本会话 WebSearch 配额已在任务开始前耗尽 200/200，无法再检索镜像）。故 (c) 目前是**单一 A 级源**，建议主会话把它标为"待第二源交叉"。

### 2-5　(d)　μ* 随 ω_c 跑动的显式陈述（A 级）

- Marsiglio eq.(14) 本身就是跑动关系；同文正文明写 "U changes to an effective U*(ω_c)"，并在 **Fig.1 图注**里给出一个具体数值约定（A 级）："The points marked with asterisks are again for the same broadened Lorentzian centred at ν₀=100 meV, but now with u*(ω_c = 1 eV) = 0.1. Note that ν_c is used as a practical cutoff for the phonon spectrum whereas ω_c is used for the Matsubara cutoff for the direct Coulomb repulsion."
- M-C p.34 补充了一个**同位素相关**的重要限定（A 级）："μ*(ω_c) is reduced from some larger value μ(ω_B) through the pseudopotential effect. … In analytical treatments this is often modelled by endowing a mass dependency to the Coulomb pseudopotential through the cutoff." — 即在同位素语境下 ω_c 与声子频率挂钩会引入质量依赖；P1/P2 的同位素工作若与 P4D 交叉，这条要注意。

---

## §3　Q3 · 算子语言下的库仑项（★本轮最大增量★）

### 3-1　KAY 的确切约定（A 级，逐页读取）

KAY II §4（pp.16–19），逐字关键片段：

- 变量代换 **eq.(37)**：`(𝔒^{1/2}Θ)_n = √(2n+1) θ_n =: ξ_n`（𝔒 的对角元是奇自然数）。
- 二次型 **eq.(38)**（p.18）逐字：
  > Q(Ξ) = Σ_n [ 1 + (2/(2n+1)) Σ_{k=1}^{n} λ(k) ] ξ_n²
  >   − Σ_{n≠m} ξ_n [ λ(n−m) / (√(2n+1)√(2m+1)) ] ξ_m
  >   − Σ_{n,m} ξ_n [ λ(n+m+1) / (√(2n+1)√(2m+1)) ] ξ_m
  >
  > where the contributions from the first line at r.h.s.(38) are positive, those from the second and third line negative.
- **eq.(39)**：`Q(Ξ) = ⟨Ξ, (ℑ − λ𝔎)Ξ⟩`，`𝔎 = −𝔎₁ + 𝔎₂ + 𝔎₃`（p.18）。
- 分量式 **eqs.(40)(41)(42)**（p.19）：ℌ₁ 对角、ℌ₂ 对角为零、ℌ₃ 全正。`𝔎_j = ∫₀^∞ ℌ_j(ϖ)P(dω)`，ϖ:=ω/2πT。
- **Proposition 2**（p.19）：三个 𝔎_j 都是 HS 算子，紧地映 ℓ²(ℕ₀)→ℓ²(ℕ₀)。
- **eq.(51)**（p.23）逐字：
  > lim_{T→0} 𝔎^(N)(P,T) = −ℑ^(N) + 2 Ξ*_N ⊗ Ξ*_N,
  >
  > with Ξ*_N = (ξ*_0, ξ*_1, …, ξ*_{N−1}) has components **ξ*_n = 1/√(2n+1)**. Since r.h.s.(50) is independent of ω, … The limiting matrix at r.h.s.(51) has an N−1-dimensional eigenspace for the eigenvalue −1 … The remaining eigenspace is {tΞ*_N : t∈ℝ}, associated with the eigenvalue **−1 + 2 Σ_{n=0}^{N−1} 1/(2n+1) ≥ 1**.
- **Theorem 5 / eq.(24)**（p.11）逐字：
  > lim_{T→0} 𝔩^(N)(P,T) = −1 + 2 Σ_{n=0}^{N−1} 1/(2n+1) =: 𝔩₀^(N).
  > … Moreover, 𝔩₀^(N) is strictly monotonically increasing with N, **diverging ∼ ln N** to ∞ as N→∞.

### 3-2　库仑项的确切算子结构（B 级组装，成分全 A 级）

对 eq.(38) 做**唯一合法的**替换：在**能隙通道**（第二、三行）把核 λ(j) → λ(j) − μ*，并把 n,m 都限制到 n,m < N_c（双边截断，§2-3）。第一行（Z 通道）不动。得

```
Q_μ(Ξ) = ⟨Ξ, ( ℑ − λ𝔎_ph + μ*·𝔎_c ) Ξ⟩ ,

𝔎_c := 2 |u⟩⟨u| − D_{N_c} ,

  u_n = 1/√(2n+1) · 1_{n<N_c} ,     D_{N_c} = diag( 1/(2n+1) · 1_{n<N_c} ) ,
  N_c := ω_c / 2πT   （费米 Matsubara 计数：|ω_m|<ω_c ⇔ m < N_c）
```

**推导要点**：`−Σ_{n≠m}(−μ*)ξ_nξ_m/(√√) − Σ_{n,m}(−μ*)ξ_nξ_m/(√√) = μ*[ 2⟨u,Ξ⟩² − Σ_n ξ_n²/(2n+1) ]`。因子 **2** 来自 λ(n−m) 与 λ(n+m+1) 两个通道各贡献一次；**−D** 来自 𝔎₂ 排除 n=m 而 𝔎₃ 不排除。

**自洽性检验（这一步把 B 级证据强度拉满）**：把上式的推导逻辑原样套到 KAY 自己的 T→0 极限（ϖ→∞ 时 λ(j)→λ 对所有固定 j，即"核变成常数"，与库仑项**完全同型**），应当复现 KAY 的 eq.(51)。检验：
- 非对角 (n≠m)：ℌ₂+ℌ₃ → 2/(√(2n+1)√(2m+1)) = 2ξ*_nξ*_m ✓ 与 eq.(51) 的 `2Ξ*⊗Ξ*` 一致。
- 对角：−ℌ₁ → −2n/(2n+1)，ℌ₃ 对角 → +1/(2n+1)，合计 −(2n−1)/(2n+1)；而 eq.(51) 的 `−1 + 2ξ*_n²` = −1 + 2/(2n+1) = −(2n−1)/(2n+1) ✓ **逐位相等**。

⇒ **库仑算子与 KAY 自己的 T→0 极限算子是同一个数学对象**（差别只在：库仑版没有 −𝔎₁ 的 −ℑ 那一份，因为 μ* 不进 Z 通道；且支撑被截断在 N_c）。

### 3-3　直答"秩一 vs 投影"（B 级）

**两个都不是。**
- **不是** −μ*·P_{N_c}（投影）。若是投影，谱只有 {0,−μ*}，那样物理上就没有 Tolmachev 对数了。
- **不是**纯秩一 −μ*|u⟩⟨u|。少了 −D 那一项，就与 KAY eq.(51) 的对角结构对不上。
- **正确形式**：`𝔎_c = 2|u⟩⟨u| − D_{N_c}`，**秩 = N_c**（D 在 span{e_0,…,e_{N_c−1}} 上可逆），是"对数发散的秩一 + 有界负对角（‖D‖≤1）"。
- 关键的一点：**u 不是截断指示向量 (1,1,…,1,0,…)，而是带 1/√(2n+1) 权重的向量**。这个权重是 KAY 的 eq.(37) 归一化带来的，也是全部对数行为的来源（Σ 1/(2n+1) ~ ½ln N_c）。

### 3-4　直答"一致有界？紧？"（B 级，★结论与团队负责人猜测相反★）

**范数（精确双边界，我自算）**：记 `S(N_c) := Σ_{n=0}^{N_c−1} 1/(2n+1) = H_{2N_c−1} − ½H_{N_c−1}`。
- 上界：`‖𝔎_c‖ ≤ ‖2|u⟩⟨u|‖ = 2‖u‖² = 2S(N_c)`。
- 下界（Rayleigh at ξ=u）：`‖𝔎_c‖ ≥ 2S(N_c) − (Σ_{n<N_c} 1/(2n+1)²)/S(N_c) ≥ 2S(N_c) − (π²/8)/S(N_c)`，用了 Σ_{n=0}^∞ 1/(2n+1)² = π²/8 = 1.2337。
- 渐近：`S(N_c) = ½ln N_c + ln2 + γ/2 + o(1) = ½ln N_c + 0.98176…`，故
  **‖𝔎_c‖ = ln N_c + 1.9635… + O(1/ln N_c)**，其中 N_c = ω_c/2πT。

于是：

| 性质 | 团队负责人猜测 | 一手核实结论 |
|---|---|---|
| 算子范数 | "= μ*，有界" | **错**。范数 = μ*·(ln N_c + 1.9635 + …)，**随截断对数发散**，不一致有界 |
| 秩 | "随 N_c 增长" | **对**（rank = N_c），但这不是坏消息 |
| 紧性 | "有界但非紧，对紧性论证要命" | **反了**。有限秩 ⇒ **紧、且 Hilbert–Schmidt**。**紧性完全没有被威胁** |

**要命的技术点因此不是紧性，而是 T-一致性**：固定 T>0 时 𝔎_c 是有限秩紧算子，KAY 的 Proposition 2（HS ⇒ 紧）无缝推广；出问题的是 T→0 时 ‖𝔎_c‖ → ∞。

**但这个发散并不新鲜、也不致命（B 级，★对项目最有利的一条★）**：KAY 自己的**声子**算子在 T→0 时同样对数发散——Theorem 5 逐字说 `𝔩₀^(N)` "diverging ∼ ln N"，Corollary 2（p.12）还给出 `𝔩(P,T) ≤ C/T²`。所以 μ* 项引入的是一个**与 KAY 已有行为同型、且量级更小（ln vs 1/T²）**的病态。P4D 不必为此发明新技术。

### 3-5　三条对推导链有直接后果的观察（B 级）

1. **不是 𝔎_ph − μ*𝔎_c，而是算子束（pencil）。** 因为 λ 乘 𝔎_ph 而 μ* **不**随 λ 缩放，正确的对象是
   `ℑ − λ𝔎_ph(P,T) + μ*𝔎_c(T) ⪰ 0`，
   临界条件是 `min spec( ℑ − λ𝔎_ph + μ*𝔎_c ) = 0`。**不能**简单写成 𝔎(P,T;μ*) = 𝔎_ph − μ*𝔎_c 再取最大本征值——那样会把 λ 和 μ* 的相对标度弄错。若坚持写成单算子，只能是 `𝔎_ph − (μ*/λ)𝔎_c`，代价是算子显式依赖 λ，变分原理 Λ=1/𝔨 的不动点结构随之改变。**这是本轮发现的最容易翻车的一处，建议主会话优先处理。**
2. **变分上界阶梯（KAY Thm 4）原样存活。** KAY 的 Λ < 1/𝔨^(N) 只依赖"把二次型限制到 N 维子空间给出 ≥ 的极小值"，与核的具体形式无关。加了 +μ*𝔎_c 项后 Q_μ 仍是二次型，同一论证逐字成立。**但**：𝔎_c 支撑在前 N_c 个坐标上，且它在 Q 里是**正**贡献（拆对），所以截断 N < N_c 时会**低估**库仑排斥；建议约束 N ≥ N_c，或显式估计被丢掉的那部分。
3. **𝔎_c 不是半正定。** 在 u^⊥ 上 ⟨ξ,𝔎_cξ⟩ = −Σξ_n²/(2n+1) < 0。谱结构：一个 ≈ ln N_c 的大正本征值（配对通道，就是 Tolmachev 对数），加上 N_c−1 个落在 [−1, 0) 的负本征值。这意味着 −μ*𝔎_c 在**顶通道**是压制 Tc（物理正确），但在若干次级通道上是**弱增强**——若 P4D 要做"μ* 单调降低 Tc"的严格陈述，这一点需要单独论证，不能想当然。

---

## §4　Q4 · 含库仑的严格 / 算子理论 Tc 界

### 4-1　(i) Combescot, PRB 51, 11625 (1995)

- **A 级新增**：标题现由 KAY I 参考表（p.46）逐字确认为 **"Strong-coupling limit of Eliashberg theory"**（此前工作区只有卷页号）。
- **全文仍未取得**：无 arXiv；`journals.aps.org/prb/abstract/10.1103/PhysRevB.51.11625` 返回 **HTTP 403**；本会话 WebSearch 配额已耗尽（200/200），无法再找镜像/引用链。
- **是否含库仑：未知，仍挂账。** 从标题（"strong-coupling limit"）推测它讨论的正是 λ→∞ 渐近，**因此它极可能是回答 §4-3 那个 sharp 问题的最直接来源**，优先级应当提高。需要人工访问 APS 或机构图书馆。
- KAY I 正文 pp.1–40 未找到引用 [Co] 的位置（前次调研结论，本轮未推翻）。

### 4-2　(ii) Hainzl–Seiringer 纲领是否碰库仑

- 沿用既有结论（记忆基线，C→B 级）：Hainzl–Seiringer 的严格 BCS 纲领是**纯瞬时核**，其综述全文无 "Eliashberg"/"retarded"。本轮**没有**新证据推翻或加强。
- **本轮可新增的一条 A 级情报**：KAY I 参考表（pp.46–47）列出的数学物理侧文献全集是 —— [HSa] Hainzl & Seiringer, PRB **77**, 184517 (2008)；[HSb] Hainzl & Seiringer, J. Math. Phys. **57**, 021101 (2016)；[HLR] Henheik, Lauritsen & Roos, Rev. Math. Phys. 2360005 (2023)；[LT] Langmann & Triola, PRB **108**, 104503 (2023)；[LHSB] Langmann, Hainzl, Seiringer & Balatsky, arXiv:2207.01825 (2022)；[FL] Freericks & Lieb, PRB **51**, 2812 (1995)。**KAY 自己在做文献综述时，找不出任何一篇"含库仑赝势的严格 Tc 界"来引用**——这本身是"该方向为空白"的间接但有力的证据（B 级推断）。
- KAY I p.5 脚注 1 逐字说明这批数学物理工作与声子超导的关系（A 级）："In recent years this model has received much attention by mathematical physicists, see e.g. [HSa], [HSb], [LT], [HLR], but these studies do not alleviate the criticism of the BCS theory of phonon-mediated superconductivity."

### 4-3　(iii) ★sharp 事实核实：强耦合渐近里 μ* 进入的阶★

这是本轮**唯一没有拿到 A 级直接答案**的问题，但我能给出一个界定清楚的结论。

**A 级（原文逐字，均无 μ*）**：
- M-C cond-mat/0106143 **p.32, eq.(94)**：
  > The correct asymptotic result for T_c is [140,193]: **T_c = 0.183 √λ ω_E**, as λ → ∞. In obtaining Eq. (94), an Einstein spectrum has been assumed (this is not required) …
  —— 式中**没有 μ\***，且**上下文完全没有提到 μ\* 的修正**。
- Marsiglio arXiv:1911.05065 **p.6**：
  > It is well known that Eliashberg theory predicts that T_c increases with both frequency and coupling strength as **T_c ≈ ν₀√λ in the asymptotic limit.**
  —— 同样无 μ*。

**但这两处都只是"没写 μ\*"，不等于原文断言"与 μ\* 无关"。** 我没有找到任何一句逐字说 "independent of μ*" 的原文。**AD 1975 原文（含其 eq.(22)）本轮未取得**（APS 付费墙 + 搜索配额耗尽）。

**B 级（我自算，可复核）**：把 Allen–Dynes **拟合公式**（f₁,f₂ 形式已由 `surveys/P3_allen_dynes_original_2026-09-09.md` 用 Xie 2022 eq.(2)(3)(4) 强 B 级锁定）取 λ→∞ 极限：
- f₁ → √( λ / (2.46(1+3.8μ*)) )
- f₂ → r ≡ ω̄₂/ω_log
- 指数 → exp( −1.04 / (1 − 0.62μ*) )

⇒ `T_c → C_AD(μ*) · √(λ⟨ω²⟩)`，其中 **C_AD(μ*) = (1/1.2)·[2.46(1+3.8μ*)]^{−1/2}·exp(−1.04/(1−0.62μ*))**（注意配对量自动变成 ω̄₂√λ = √(λ⟨ω²⟩)）。数值：

| μ* | C_AD(μ*) | 相对 μ*=0 |
|---|---|---|
| 0 | **0.18780** | — |
| 0.10 | 0.14927 | −20.5% |
| 0.13 | 0.14026 | −25.3% |

**⚠ 两个副产品**：(1) μ*=0 时得 **0.18780**，与工作区 `P3_tc_max_pressure.md` §7 自算的 **0.1878** 逐位吻合 —— 这独立确认了工作区那个数**不是代数错误，而是"AD 拟合公式的渐近"这个对象本身**（它与真渐近 0.1827 差 +2.8%，差异来自 AD 拟合 vs AD 真渐近，不是算错）。(2) **在 AD 拟合公式里，μ\* 确确实实进入领头 √λ 系数，而且是 20–25% 量级的大效应。**

**B 级（我的物理判断，明确标为推断）**：真渐近常数 C_∞ 应当**与 μ\* 无关**，理由是：μ* 项带一个**固定能标** ω_c，而 T_c ~ C√(λ⟨ω²⟩) 在 λ→∞ 时无界增长，故 N_c = ω_c/2πT_c → 0，库仑项支撑上的 Matsubara 点数趋于零，其贡献（§3-4 里 ~μ* ln N_c）退化为 O(μ*) 的有界移动，相对 √λ 是**次领头**。而 AD 拟合公式**不含 ω_c**，它把 μ* 的效应硬编码成与 λ 无关的乘性因子，所以它在 λ→∞ 外推时**人为地**保留了 μ* 依赖——上表那 20–25% 是**外推假象**。这也与 KAY 的 C_∞ = (1/2π)√𝔤(2) = 0.1827262477 推导里**根本没有出现任何截断能标**相容。

**⚠ 反向考虑（诚实边界）**：上述论证同时说明，在 λ→∞ 且 T_c ≳ ω_c 的区域，μ* 赝势构造本身（前提是 ω_ph ≪ ω_c ≪ E_F）已经失效。所以"μ* 不改领头常数"可能**为真但空洞**——它之所以不改，是因为在那个极限里 μ* 这个概念已经没有意义。P4D 若要用这条，应当写成"在 λ 有限但大的**物理相关**区间内，μ* 的效应是相对 √λ 领头项的 O(1) 修正"，而不是简单宣称"C_∞ 与 μ* 无关"。

**⇒ 直答团队负责人的二选一**：现有证据倾向 **"μ\* 不改领头常数、只给修正"**，但**证据等级只到 B**，且论证方式必须小心（见上一段）。**要把它升到 A，最短路径是拿到 Combescot PRB 51, 11625 (1995) 全文或 AD 1975 原文 eq.(22) 附近的正文**——这两件事现在合并成了同一个待办。

---

## §5　诚实边界（本轮**没有**做到的事）

1. **Allen–Dynes PRB 12, 905 (1975) 原文未取得**（APS 付费墙，无 arXiv）。故 AD eq.(22) 的确切形式、以及 μ* 在其中出现与否，**仍是 C/未知级**。本会话 WebSearch 配额在任务开始前已耗尽（200/200），无法检索镜像或代购路径。
2. **Combescot PRB 51, 11625 (1995) 全文未取得**（APS 403）。标题已 A 级确认，内容未知。
3. **Bergmann–Rainer Z. Phys. 263, 59 (1973) 原文未取得**。其 eq.(7) 的确切形式（尤其 μ*(N) 的 N 依赖是如何写的）目前**只有 KAY I 脚注 4 的转述**——这是 A 级的"KAY 怎么说"，但不是 A 级的"BR 原文怎么写"。§2-3 那个"双边截断"的判断因此带一点风险。
4. **Morel–Anderson PR 125, 1263 (1962) 原文、Allen–Mitrović SSP 37 (1982) 书章未取得**。§2-4 的 μ* 对数公式目前是**单一 A 级源**（Marsiglio 2020 eq.14）。
5. **"是否有人做过含库仑的严格/算子理论 Tc 界"没有做系统检索**（搜索配额耗尽）。§4-2 的结论是从 KAY 参考文献表的**缺席**推断的，属间接证据。
6. §3 的算子结构推导（3-2 到 3-5）**全部是我的组装**，成分（KAY eqs.37/38/39/40–42/51 与 Thm 5；Marsiglio eqs.22/23）逐字取自原文属 A 级，但**组装本身是 B 级，未经任何文献背书**。特别是 3-5 的第 1 条（算子束 vs 单算子）建议主会话独立复核。

---

## §6　建议的后续检索方向（按性价比排序）

1. **合并办：AD 1975 原文 + Combescot 1995 全文。** 两者共同决定 §4-3 那个 sharp 事实能否升 A。需人工访问 APS / 机构图书馆，或下个会话恢复搜索配额后找机构镜像。
2. **Bergmann–Rainer 1973 eq.(7) 原文**，用来锁死 §2-3 的截断约定（单边 vs 双边）与 μ*(N) 的确切 N 依赖形式。Z. Phys. 是 Springer，可能有 static-content 直链（参见既有经验：Springer static-content URL 曾直接命中）。
3. **Morel–Anderson 1962 原文** 作为 §2-4 的第二源。
4. **"rigorous bound Tc Coulomb pseudopotential" 型系统检索**（下会话配额恢复后），重点查 Frank / Hainzl / Seiringer / Solovej 系列里是否有"BCS + 库仑"的严格工作，以及 Henheik–Lauritsen–Roos 的后续。
5. **数值旁证（无需检索，主会话可自己做）**：直接对 `ℑ − λ𝔎_ph^(N) + μ*𝔎_c^(N_c)` 做数值谱计算，验证 §3-2 的算子形式能否复现 Marsiglio Fig.1 里 u*(ω_c=1 eV)=0.1 那组带星号的 Tc 数据点。这是检验整个 §3 组装是否正确的**最便宜的证伪实验**。
