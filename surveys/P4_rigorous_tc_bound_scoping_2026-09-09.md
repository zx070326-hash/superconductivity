# P4「严格的声子介导超导 Tc 上界」立项调研

调研 agent：researcher　日期：2026-09-09
任务：把 P3（λ_* 谱形展开，κ=5.6035 / c₃=0.8959）严格化升格所需的前沿地形钉死，并收敛出候选精确命题。

**证据等级**
- **A** = 一手 arXiv 原文逐字（我用 PDF→Read 或 HTML 直读到的原句/公式）
- **B** = 权威二手逐字，或我对 A 级数字做的确定性推算（标 `[本报告自算]`）
- **C** = 图注 / 推断 / 元数据级 / 未核实

**本轮取全文的一手材料**：arXiv:2409.00532v3（KAY II，PDF 页 1–18、27–33 逐页读）、arXiv:2409.19562v2（本地 PDF 页 1–3）。**HTML/摘要级**：2409.00533v2（KAY I）、2409.02121v2（KAY III）、1511.01995（Hainzl–Seiringer 综述，ar5iv）、2207.01825、2407.00796。**元数据级**：arXiv API 三次普查。

---

## §0　一句话结论

**"确定各向同性 Eliashberg 大 λ Tc 渐近的 sharp constant"这个提法本身已经不开放了**——KAY 2025 已经把它定义清楚、证明它是渐近常数、并算到 10 位有效数字：**C_∞ = 0.1827262477…**，而且**它就是 Allen–Dynes 1975 eq.(22) 那个常数**（KAY 原文明说）。**真正开放的是另外三件事**：(a) 这个 10 位数字只有 Maple/Mathematica 浮点，没有严格包围；(b) 已证的显式上界常数比它大 2.034 倍，这个缺口作者自陈 "not optimal"；(c) 把渐近上界升级为**全 λ** 上界是 KAY 白纸黑字的 Conjecture。这三条都良定义、可证、且与 P3 现有机器有接口。

---

## §1　第一块 · 算子本征值型的 Eliashberg Tc 界（KAY 系）

### 1.1 三部曲的完整定位（等级 A，arXiv API 普查 + 各文摘要）

工作区 P3 survey §5 把这条线记为"最该补的两篇"。实际是**三篇**，全部 2024-08-31 同日投稿、2025 全部发表在 J. Stat. Phys.：

| 编号 | 标题 | 期刊 | 本轮状态 |
|---|---|---|---|
| **2409.00533**v2 | Bounds on T_c … **I: The γ-model** | J. Stat. Phys. **192, art.69**, 35pp (2025) | HTML 已读 |
| **2409.00532**v3 | Bounds on T_c … **II: Dispersive phonons** | J. Stat. Phys. 92, art.94, 32pp (2025) | **PDF 全文逐页已读** |
| **2409.02121**v2 | Bounds on T_c … **III: Einstein phonons** | J. Stat. Phys. 92, art.93, 26pp (2025) | HTML 已读 |

作者三人全同：Michael K.-H. Kiessling（Rutgers 数学系）、B. L. Altshuler（Columbia 物理）、E. A. Yuzbashyan（Rutgers 物理）。II 的 arXiv comments 逐字："44 pages, 2 figures; several typos in the previous version have been corrected, bibliography has been updated. Appeared open access in J. Statist. Phys. vol.92, art.94 (2025)"。

⚠ **卷号在 arXiv 元数据里自相矛盾**：I 记 vol.192，II/III 记 vol.92。P3 survey §5 记的是 "J. Stat. Phys. 192, 93 & 94"。**主会话引用前请核 DOI**，我倾向 192 是对的（92 疑为 arXiv 录入脱字）。等级 C。

### 1.2 紧自伴算子与变分原理（等级 A，逐字）

**摘要（p.2）逐字**，这是全套定义的浓缩：

> "The standard Eliashberg theory of superconductivity is studied, in which the effective electron-electron interactions are modelled as mediated by generally dispersive phonons, with Eliashberg spectral function α²F(ω) ≥ 0 that is ∝ ω² for small ω > 0 and vanishes for large ω. The Eliashberg function also defines the electron-phonon coupling strength λ := 2∫_{ℝ₊} (α²F(ω)/ω)dω. Setting (2α²F(ω)/ω)dω =: λP(dω), formally defining a probability measure P(dω) with compact support, and assuming as usual that the phase transition between normal and superconductivity coincides with the linear stability boundary 𝒮_c of the normal region in the (λ,P,T) parameter space against perturbations toward the superconducting region, it is shown that this *critical hypersurface* 𝒮_c is a graph of a function Λ(P,T)."

**⚠ 这里有一个对 P3 极其重要的记号巧合**：KAY 的 **P(dω) 就是 P3 用的形状函数 P(ω)**。Semenok SI Eq.(S10–S11) 定义 P(ω) ≡ (1/λ)·2α²F(ω)/ω；KAY eq.(7) 定义 2α²F(ω)/ω =: λP'(ω)。**逐字同一个对象。** P3 的 σ、γ₃ 是这个 P 的相对宽度与偏度——可以原封不动搬进 KAY 框架。

**核心链条（全部 A 级逐字）**：

- **Theorem 1（p.7）**："The positive (λ,P,T) cone of the standard Eliashberg model consists of two simply connected regions. In one region the normal state is unstable against small superconducting perturbations, in the other region it is linearly stable. The boundary between the two regions, called the critical hypersurface 𝒮_c, is a graph over the set {(P,T)}, i.e. 𝒮_c = {(λ,P,T) : λ = Λ(P,T)}."

- **Theorem 3（p.8–9）**："The function Λ(P,T) is determined by the following variational principle, Λ(P,T) = 1/𝔨(P,T), (9) where 𝔨(P,T) > 0 is the largest eigenvalue of a compact self-adjoint operator 𝔎(P,T) on the Hilbert space of square-summable sequences over the non-negative integers, with 𝔎(P,T) explicitly constructed in section 4.1, see (39)–(42)."

  即 **ℓ²(ℕ₀)**，且 **𝔎 = −𝔎₁ + 𝔎₂ + 𝔎₃**（eq.39 后文），三块都是**关于 P(dω) 的平均**：𝔎^(4)(P,T) = ∫₀^∞ ℌ^(4)(ϖ)P(dω)，**ϖ := ω/2πT**。

- **Theorem 4（p.9，上界序列 → Tc 下界）**："For all N ∈ ℕ, Λ(P,T) < 1/𝔨^(N)(P,T), where 𝔨^(N)(P,T) is the largest eigenvalue of 𝔎^(N)(P,T), the restriction of 𝔎(P,T) to the first N components of ℓ²(ℕ₀). The eigenvalues 𝔨^(N)(P,T) can be explicitly computed for N ∈ {1,2,3,4}."

  最低阶（eq.10）**逐字**：𝔨^(1)(P,T) = ∫₀^∞ ω²/(ω² + (2πT)²) P(dω)。N=2 用二阶特征多项式（eq.11）、N=3 用三次方 arccos 解（eq.12）、N=4 用四次方 resolvent cubic（eq.15–22）。矩阵元（eq.23）：ℌ^(4) 的分量全由 **⟦n⟧(ϖ) := ϖ²/(n²+ϖ²) ≡ ω²/((2nπT)²+ω²)** 组合而成。

- **Theorem 6（p.11，下界 → Tc 上界）**："Let (P,T) be given. Then 𝔨(P,T) ≤ 𝔨\*(P,T), where 𝔨\*(P,T) = 𝔨^(1)(P,T) + 2√((2^{1+ε}−1)ζ(1+ε)ζ(5−ε))⟨ϖ²⟩ (26) with ε = 0.65, and where ⟨ϖ²⟩ := ∫₀^∞ ϖ²P(dω) with ϖ = ω/2πT."

- **Corollary 4（p.14，显式 Tc 上界，eq.32）逐字**：

  T_c^♯(λ,P) = (√⟨ω²⟩/2π)·√( ½( λ(1+b) − 1 + √((λ(1+b)−1)² + 4bλ) ) )，
  b := 2((2^{1+ε}−1)ζ(1+ε)ζ(5−ε))^{1/2}，ε = 0.65。

- **Corollary 4 续（eq.33，Tc 下界）**：对 λ > Ω̄²/⟨ω²⟩，T_c(λ,P) > T_c^♭(λ,P) = (1/2π)(λ⟨ω²⟩ − Ω̄²)^{1/2}(P)。

### 1.3 「2.034」的确切定义与出处（**这是本轮最核心的核实**）

**出处一，摘要 p.2 逐字**：

> "Also obtained is an upper bound on T_c(λ,P). Although not optimal, it agrees with the asymptotic form **T_c(λ,P) ~ C√⟨ω²⟩√λ** valid for λ ~ ∞, given P, though with a constant C that is a factor **≈ 2.034 larger than the sharp constant**; here, ⟨ω²⟩ := ∫_{ℝ₊} ω²P(dω)."

**出处二，正文 p.5 逐字（给出 sharp constant 的数值）**：

> "Also obtained is an upper bound on T_c(λ,P), which is not optimal yet agrees with the asymptotic behavior T_c(λ,P) ~ C√⟨ω²⟩√λ for large enough λ, given P, though with **C ≈ 2.034 C_∞, where C_∞ = 0.1827262477... is the sharp constant.** Here, ⟨ω²⟩ := ∫_{ℝ₊} ω²P(dω)."

**⟨ω²⟩ 的定义与 Semenok Eq.(8) 一致**（Semenok：⟨ω²⟩ = (2/λ)∫α²F(ω)ω dω；KAY：∫ω²P(dω)，而 P(dω) = 2α²F/(λω)dω ⟹ ∫ω²·2α²F/(λω)dω = (2/λ)∫α²F·ω dω）。**逐位同一。** 可直接互换。

**sharp constant 的来源，Remark 1（p.8）逐字**——这是 10 位数字的一手出处：

> "With the help of Maple and, independently, Mathematica, we found that when the rank is increased beyond **N = 200**, then 10 significant decimal places of the lower approximation to T_c(g, γ = 2) have stabilized, yielding **(1/g)T_c(g,2) = 0.1827262477....**"

**Conjecture 2（p.14，eq.31）逐字**——**这是全套的开放核心**：

> "**Conjecture 2**: There is a critical temperature T_c(λ,P) > 0 which for all P ∈ 𝒫 and λ > 0 is bounded above by T_c(λ,P) < T_c^~(λ,P), with **T_c^~(λ,P) := (1/2π)√(𝔤(2)⟨ω²⟩λ)** (31); here, 𝔤(2) is the spectral radius of 𝔊(γ) at γ = 2."

所以 **C_∞ = (1/2π)√𝔤(2)**。

### 1.4 **判定：C_∞ 与 Allen–Dynes 的 0.1827 是同一个常数——KAY 原文自己说的**

**Remark 3（p.14）逐字，这是决定性证据**：

> "The asymptotic behavior T_c(λ,P) ~ (1/2π)√(𝔤(2)⟨ω²⟩λ) was **first stated in [AD], see their eq.(22)**. Apparently they did not notice that this is presumably an upper bound on T_c for all λ > 0."

[AD] = Allen & Dynes, PRB 12, 905 (1975)。**KAY III 补充（HTML 级，等级 B）**：AD 报的 C = 0.182 是 **64-mode 截断**的数值结果，KAY 评其 "was not very accurate"；KAY 自己 N=4 的闭式已给 0.1825，四舍五入到 0.183。

**再一条独立佐证（等级 A，本地 PDF 2409.19562 p.2 eq.(6) + p.3）**：Yuzbashyan–Altshuler–Patra 定义 **g² = 2∫dω ω α²F(ω) = Σ_k λ_k ω_k²**，而 λ⟨ω²⟩ = (2/λ)∫α²F ω dω · λ = 2∫α²F ω dω = **g²**。所以 g = √(λ⟨ω²⟩)，PRL 里的 "T_c = 0.1827 in the limit λ → ∞" 与 Semenok Eq.(7) 的 "T_c → 0.18√(λ⟨ω²⟩)" 与 KAY 的 C_∞√⟨ω²⟩√λ **三者是同一个式子的三种写法**。

**⟹ 判定：可对接，且是同一个常数。P3 记的 0.1827 = KAY 的 C_∞ = AD eq.(22)。三处一致，无冲突。**

### 1.5 **但 0.1878 是另一个数——务必分开**（等级 B，`[本报告自算]`）

任务书提到的 "AD 拟合自带 0.1878"。我用 P3 survey §2.4 转录的 AD f₁f₂ 公式（经 Xie 2022 转录，B 级）自算其 λ→∞ 渐近：

取 μ\*=0：exp 因子 → exp(−1.04) = 0.353455；f₁ → (λ/2.46)^{1/2}；f₂ → r ≡ ω̄₂/ω_log。故

T_c → (0.353455/1.20)·(1/√2.46)·ω̄₂√λ = 0.294546 × 0.637578 × ω̄₂√λ = **0.187803 ω̄₂√λ**。

**四个常数必须分层（这是 P3 最容易翻车的一处）**：

| 记号 | 数值 | 是什么 | 地位 | 出处等级 |
|---|---|---|---|---|
| C^(1) | **0.159154943 = 1/2π** | N=1 截断的 Tc **下界**常数 | 已证（闭式） | A（KAY I，𝔤^(1)(γ)=1） |
| C^(4) | **0.1825137102** | N=4 截断的 Tc 下界常数 | 已证（闭式，"too unwieldy"） | A（KAY I） |
| **C_∞** | **0.1827262477…** | **真 Eliashberg 大 λ 渐近常数** | 渐近性已证；**数值只有浮点** | A（KAY I/II） |
| C_ADfit | **0.187803** | AD **拟合公式** f₁f₂ 的 λ→∞ 渐近（μ\*=0） | 拟合，非渐近理论 | B `[本报告自算]` |
| C^♯ | **≈ 0.371665** = 2.034·C_∞ | **已证的显式上界**常数 | 已证（Cor. 4） | A |

**AD 拟合公式比真常数高 +2.78%**（0.187803/0.1827262 = 1.0278）。P3 若用 AD 公式做"上界"会系统性高估 ~2.8%；若用 Semenok 的 0.182 则低估 0.4%。

⚠ **μ\* 敏感**：C_ADfit 只在 μ\*=0 时是 0.1878。μ\*=0.1 时我自算得 **0.1493**（f₁ 分母变 2.46×1.38、指数变 −1.04/0.938）。**KAY 全套框架 μ\* = 0**，所以与 C_∞ 可比的是 μ\*=0 那支。P3 用 Semenok Fig.2 的 "0.182 ω_log√λ" 时也要记住那是 μ\*=0 约定。

### 1.6 **2.034 缺口的解剖：它藏在哪三步**（这是 P4 可攻击性的关键，等级 A 原文 + B 自算）

**结构恒等式（`[本报告自算]`，B 级）**：把 Cor. 4 的 eq.(32) 展到 λ→∞，内根号 → λ(1+b)，故

C^♯ = (1/2π)√(1+b)，　而 C_∞ = (1/2π)√𝔤(2)　⟹　**2.034 = √((1+b)/𝔤(2))**。

数值核对：ε=0.65 ⟹ 2^{1.65}−1 = 2.13809，ζ(1.65) ≈ 2.163，ζ(4.35) ≈ 1.0617（我逐项求和到 n=10 + 尾）⟹ b ≈ **4.4318**，1+b ≈ 5.4318；𝔤(2) = (2π·0.1827262477)² = **1.318140**。比值 √(5.4318/1.31814) = **2.0300**，与原文 2.034 差 0.17%——**在我 ζ(1.65) 的精度内一致，恒等式确认**。

我另检查了 ε 是否优化过：ε=0.60 给 b=4.4356、ε=0.70 给 b=4.4349，**均高于 ε=0.65 的 4.4318**。⟹ **ε=0.65 是该族内的最优取值，缺口不能靠调 ε 缩小。** [B 级自算]

**⟹ 于是 λ→∞ 的"三明治"是：𝔨/⟨ϖ²⟩ ∈ [1, 5.432]，而真值 = 𝔤(2) = 1.318。缩窄缺口 ≡ 把上界 5.432 往 1.318 压。**

**证明链的损失点（PDF p.30–32 逐字读出）**：

1. **Lemma 1（eq.76）**：ρ(ℭ(η)) ≤ ρ((η𝔍+𝔎₁)^{-1})·ρ(𝔎₂+𝔎₃)。
   eq.(77) 显示 ρ(((1/λ)𝔍+𝔎₁)^{-1}) = max_{n≥0}[1/λ + (2/(2n+1))Σ_{k=1}^n ∫ϖ²/(ϖ²+k²)P]^{-1} = **λ**，取到最大值的是 **n=0**（求和为空）。而真正的优化本征模 Ξ^opt 铺在很多 n 上。**把对角因子在最有利的单一指标上取值、与 𝔎₂+𝔎₃ 脱耦——这是最大的一处松弛。**
2. **Lemma 2（eq.81）**：ρ(𝔎₂+𝔎₃) ≤ ρ(𝔎₂) + ρ(𝔎₃)。**丢掉了 𝔎₂、𝔎₃ 之间的一切抵消，也丢掉了 −𝔎₁ 的负贡献**（注意真的 𝔎 = −𝔎₁+𝔎₂+𝔎₃，𝔤(2)=1.318 已经含 −𝔊₁ 的压低；而界里 𝔎₁ 只以 (η𝔍+𝔎₁)^{-1} 的形式出现且被在 n=0 处取值 ⟹ 完全没吃到）。
3. **Proposition 7（eq.83）**：ρ(𝔎_j(P,T)) ≤ ρ(𝔊_j(2))∫ϖ²P(dω)，j∈{2,3}。证明逐字："Since **ϖ² + k² > k²**, the claim of the proposition follows." **注意：在 λ→∞（等价 ϖ→0）极限下这一步是渐近紧的**，所以它**不是**渐近缺口的来源。
4. ρ(𝔊₂(2))、ρ(𝔊₃(2)) 的上估来自 **KAY I 的 Propositions 5 与 6**（ζ(1+ε)ζ(5−ε) 就产在那里）。原文 p.32 逐字："closer inspection of the proof of Proposition 6 in [KAYa] reveals that the estimate ϖ²+1 > 1, used in the proof of Proposition 7 above, is not needed … avoiding the unnecessary estimate ϖ²+1>1 for this particular contribution to the upper bound on ρ(𝔎₃) gives the bound stated in Theorem 6."

**⟹ 结论：2.034 的绝大部分来自 Lemma 1 + Lemma 2 这两步"解耦 + 次可加"，而不是来自 ζ 估计的粗糙或 ε 的选择。攻击点明确。**

### 1.7 KAY 自陈的开放项清单（等级 A）

| # | 陈述 | 状态 |
|---|---|---|
| **Conj. 1**（II p.7 逐字） | "There is a critical temperature T_c > 0, depending on α²F(ω), such that for temperatures T ≥ T_c, the spin chain of the normal state **N** is the unique minimizer of H(**S**\|**N**), whereas at temperatures T < T_c a spin chain **S** ≠ **N** … Moreover, the phase transition at T_c from normal to superconductivity is continuous." | **未证**。"A confirmation of the continuous normal-to-superconductivity phase transition requires a study of the **nonlinear** Eliashberg gap equations, which we postpone to a later publication." |
| **Conj. 2**（II p.14, eq.31） | T_c(λ,P) < (1/2π)√(𝔤(2)⟨ω²⟩λ) 对**所有** λ>0 | **未证**（KAY III 里编号为其 Conjecture 1，仅有数值证据） |
| **Conj. 3**（II p.29, eq.74） | ⟨Θ_N,𝔊̂^(N)(γ)Θ_N⟩/⟨Θ_N,𝔇^(N)Θ_N⟩ ≥ (1/N²)Σ_{k=1}^{2N−1} min{k,2N−k}/k^γ，在非零、非负、递减序列集上，等号 iff Θ_N 为常序列 | **未证**（Remark 7：已证常序列是局部极小，"We suspect that these constant sequences in fact are global minimizers"） |
| 单调性（III） | ϖ ↦ 𝔥(ϖ) 在全 ℝ₊ 上严格递增 | **只证到 ϖ ≤ √2**（II 的 Lemma 3 逐字）；N=3,4 情形 KAY III 自陈 "vindicated through visual inspection of the plots" |
| 𝒞_c^(N) 是 [λ_N,∞) 上的图（III） | 一般 N | **只证了 N ∈ {1,2,3,4}** |
| 𝔤(2) 的值 | 0.1827262477… | **只有 Maple/Mathematica 浮点，无严格包围，无闭式** |

**⚠ 我做了后续检查（arXiv API，au:"Yuzbashyan" 按时间倒序）：截至 2026-09-09，Yuzbashyan 2024 年后共 10+ 篇，主题为 Kondo 精确解、Richardson–Gaudin、测量诱导相变、Trotter 混沌、LaSc₂H₂₄ 等，**没有任何一篇是 KAY IV 或"非线性 Eliashberg 能隙方程"的后续**。Kiessling 在 `cat:math-ph AND abs:"Eliashberg"` 普查里最新条目仍是 2409.02121。⟹ **上表所有开放项截至 2026-09 仍然开放。**（限制：普查基于 25/71 条分页 + 单一 cat 过滤，非穷尽。等级 B。）

### 1.8 收敛的 Tc 下界阶梯（等级 A，KAY I HTML；数值为 g=1 单位）

| N | T_c^(N)(2) | 相对 C_∞ 的误差 `[自算]` |
|---|---|---|
| 1 | 0.159154943（= 1/2π 精确） | −12.90% |
| 2 | 0.1796160944 | −1.702% |
| 3 | 0.1820383 | −0.377% |
| 4 | 0.1825137102 | −0.116% |
| ∞ | **0.1827262477…** | — |

局部收敛阶 `[自算]`：N=3→4 的误差比给 ln(3.237)/ln(4/3) ≈ **4.1**，即 ~N^{-4}。与"N>200 才稳定 10 位"一致（代数收敛）。**这条阶梯是 §5 候选命题 P4-B 可行性的直接依据。**

### 1.9 KAY III（Einstein 声子）的具体数字（等级 B，HTML 摘要级）

- Tc = Ω f(λ) 的**存在性**证到 **λ > 0.77**；Theorem 5 给 λ\* < **0.7670**。
- N=1 下界闭式：**T_c^(1)(λ,Ω) = (Ω/2π)√(λ−1)**，λ>1。
- λ₂ = 3/5，λ₃ = 15/31 = 0.48387…，λ₄ = 105/247 = 0.4251…；1/𝔥^(3)(2) = 1.0158，1/𝔥^(4)(2) = 0.7670。
- 与 McMillan **无**比较；与 AD 的比较只在"64-mode 不够准"那一处。

### 1.10 **"严格化税"：把 Semenok 的头条界换成已证的界会付出什么**（等级 B，`[本报告自算]`，A 级输入）

这是本轮对 P3 **最可直接引用**的一条推论。取 Semenok Eq.(9) 的场景（Einstein 谱、λ = λ_* = 3.6915、⟨ω²⟩ = ω_max²）：

| 版本 | T_c 上界 | 说明 |
|---|---|---|
| Semenok Eq.(9)（**数值 Eliashberg 求解，非定理**） | **0.3175 ω_max** | "we used a standard algorithm to numerically compute T_c^E = 0.3175 ω_max" |
| KAY **Conjecture 2** 若被证明 | C_∞√(3.6915)·ω_max = **0.3511 ω_max** | 0.1827262477 × 1.921327 |
| KAY **今天已证的** Cor. 4 | 2.034 × 0.3511 = **0.7141 ω_max** | 唯一有定理背书的那条 |

代入 Semenok 用的 ω₀ = 3000 K（固体原子氢 500 GPa 谐波谱）：

- Semenok Eq.(11) 头条：**T_c^metastable < 950 K/√A**（= 0.3175×3000）
- **今天可证的**：**< 2142 K/√A**
- Conjecture 2 一旦证明：**< 1053 K/√A**

**⟹ 缺口 2.034 在物理单位上值 ~1090 K/√A。** 而且注意 **0.3175 < 0.3511 ✓**——Semenok 的独立数值解落在 KAY Conjecture 2 的界内且差 9.6%，这是 Conjecture 2 在 λ=3.69 处的一个**独立数值佐证**（两组人、两套算法）。

---

## §2　Yuzbashyan–Altshuler–Patra 2409.19562：λ_* 稳定界的严格地位

**判定：不是定理，是物理层面的"证明"，且带一个自陈的近似。**

**证据（等级 A，本地 PDF 页 1–3 逐字）**：

- 摘要用词："We show that thermal equilibrium … becomes kinetically unstable … **We prove** that negative electronic specific heat, C_el < 0, is sufficient to trigger the instability."
- 正文 p.1 右栏："Here, we **prove rigorously** that thermal equilibrium … becomes *intrinsically unstable* whenever C_el < 0."
- **但全文没有 Theorem/Lemma/Proposition 环境**（这是 PRL）。论证结构是：ME 理论导出的动理学方程 eq.(1)–(3) → 自能 eq.(4) → 逆准粒子权重 Z_0^{-1}(E) eq.(5) → 判据 eq.(10) min_E Z_0^{-1}(E) < 0 → 线性稳定性分析 eq.(12)–(14)（广义本征值 γ）→ Fig.1 的**全非线性数值模拟**佐证。
- **作者自陈的近似（p.3 逐字，重要）**："It is important to mention that the stability condition (10) **underestimates λ_c due to the quasiclassical nature of the kinetic equation (1)**, which ignores the quantum uncertainty in the electron energy [40]. … the uncertainty ΔE ∼ Ω. However, **δE ≥ ΔE only for λ ≳ 3.40**."
- 数值：Einstein 声子 T=0.183g 时 Z_0(E) 触轴于 E_min = 1.44Ω、λ_c = 1.27；Fig.2 在 T_0=0.20g 给 λ_c = 1.4788。

**⟹ 分层结论（与 P3 survey §1.4 的记录一致并补强）**：
- **λ_* = 3.6915（Einstein）/ 4.72（Debye）本身不在这篇 PRL 里**，在 Semenok–Altshuler–Yuzbashyan 2407.12922（Adv. Mater.）。
- 该 PRL 提供的是"C_el<0 ⟹ 动理学失稳"这一**机制性论证**，其严格性是物理级（ME 理论内部 + 动理学方程的线性稳定性 + 数值模拟），**不是数学定理**。
- ξ_c = 0.5 是 Semenok 明说的 conjecture（P3 survey §1.3 已有 A 级原句）。
- **⟹ λ_*/λ ≲ 4 是"启发式 + 数值支撑的物理判据"，不是命题也不是定理。** 这一点对 P4 立项很关键：**P3 的 λ_* 谱形展开（κ、c₃）是在一个非严格的物理判据之上做的严格泛函分析**——数学是严格的，但被展开的对象不是定理。而 KAY 线是相反的：对象（Eliashberg 方程的线性稳定边界）有定理，但常数不紧。**P4 的升格价值正在于把 P3 的展开机器挪到有定理背书的对象上。**

---

## §3　别的组给过 Eliashberg 大 λ 渐近常数吗

**普查（等级 A/B）**：arXiv API，`cat:math-ph AND abs:"Eliashberg"`，按时间倒序 30 条 → **与超导相关的只有 4 篇**：KAY I/II/III（2024-08-31）+ **Langmann, Hainzl, Seiringer, Balatsky, arXiv:2207.01825**。其余 11 条全部是 **Yakov** Eliashberg（辛几何 / SFT / Gromov–Eliashberg 定理），同姓不同人。另有 Sigal–Tzaneteas 两篇是 Gorkov–Eliashberg–Schmid **时变 GL** 的 Abrikosov 涡旋格稳定性，不是 Tc。

| 来源 | 给的常数 | 严格性 | 等级 |
|---|---|---|---|
| **Allen–Dynes 1975, PRB 12, 905, eq.(22)** | C ≈ 0.182（64-mode 截断数值） | 非严格；KAY 评 "not very accurate" | B（经 KAY II Remark 3 / KAY III 转述） |
| **KAY I/II/III 2025** | C_∞ = 0.1827262477（N>200 浮点）；已证上界 2.034C_∞ | **唯一严格线** | A |
| **Esterlis–Kivelson–Scalapino 2018**（1806.00488） | T_c ≤ A_max ℏω̄，A_max ≈ 1/10 | 启发式 + DQMC；**且机制是 ME 崩溃/双极化子，不是 Eliashberg 内的渐近常数** | A（P3 survey §2.1 已有全文逐字） |
| **Combescot, PRB 51, 11625 (1995)**, "Strong-coupling limit of Eliashberg theory" | 未知 | 未取原文 | **C** — KAY III 把它列为先行工作，无 arXiv，本轮未取。**建议主会话补。** |
| Chubukov–Abanov γ-model（KAY I 记为 [MC]/[CAWW]） | γ 模型的物理来源 | 物理级 | C |
| Marsiglio | 未找到严格常数 | — | C（未穷尽检索） |

**⟹ 判定：除 KAY 外，没有第二个组给出过 Eliashberg 大 λ Tc 渐近的严格常数。** 限制：普查依赖 math-ph 交叉分类，未交叉分类的论文会漏；Combescot 1995 未核实。

---

## §4　第二块 · 数学 BCS 的严格 Tc 纲领（Hainzl–Seiringer 系）

### 4.1 **关键判别：这套纲领覆盖推迟核吗？——不覆盖。**

**证据（等级 A−，ar5iv 全文抓取，1511.01995 = Hainzl & Seiringer, J. Math. Phys. 57, 021101 (2016) 综述）**：

- 泛函 eq.(1.2)/(2.30) 的配对项是 **∬ V(x−y)|α(x,y)|²dxdy**——**瞬时**两体势。
- 综述原文只讨论"局域 vs 空间非局域"的推广："All our results can easily be extended to the case of non-local interactions as well"——**这是空间非局域，不是时间/频率推迟**。
- 哈密顿量 (2.4)–(2.5) 是标准等时两体相互作用。
- **"Eliashberg" 与 "retarded" 两个词在全文中一次都没出现。**

**⟹ 这是本轮最干净的一条结论：Hainzl–Seiringer 严格纲领 = 纯瞬时 BCS。声子推迟版的严格 Tc 界在这条线上是空白。**

### 4.2 这套纲领证到了什么（等级 A 综述 + A 元数据）

| 结果 | 出处 | 等级 |
|---|---|---|
| **线性判据定义 Tc**：T_c = inf{T : K_T + V ≥ 0}（eq. 3.28），等价于 K_{T_c}+V 最低本征值 = 0；"there is a unique critical temperature 0 ≤ T_c(V) < ∞ that separates two phases" | 综述 §3 | A |
| **弱耦合严格渐近**：lim_{λ→0} λ ln(μ/T_c(λV)) = −1/(√μ e_μ(V))，e_μ(V) = inf spec 𝒱_μ | Frank–Hainzl–Naboko–Seiringer, **J. Geom. Anal. 17, 559 (2007)**, arXiv:0704.3564 | A |
| Tc > 0 的充分条件：∫V(x)dx < 0 | 综述 §3 | A |
| **Tc 与能隙的普适比**（与势和化学势无关） | Hainzl–Seiringer, **PRB 77, 184517 (2008)**, arXiv:0801.4159 | B（元数据级） |
| **BCS → Ginzburg–Landau 微观导出** | Frank–Hainzl–Seiringer–Solovej, **JAMS 25, 667 (2012)**, arXiv:1102.4001 | A（元数据） |
| 低密度 → BEC / Gross–Pitaevskii | Hainzl–Seiringer, Lett. Math. Phys. 100, 119 (2012), arXiv:1105.1100 | A（元数据） |
| 弱均匀**磁场**下 Tc 线性下降（严格化 WHH） | Frank–Hainzl–Langmann, J. Spectr. Theory 9, 1005 (2019), arXiv:1706.05686 | A（元数据） |
| 一般外场下的 GL + Tc 移动 | Deuchert–Hainzl–Maier, **Calc. Var. PDE 62, 203 (2023)**, arXiv:2210.09356；Prob. Math. Phys. 4, 1 (2023), arXiv:2105.05623 | A（元数据） |
| **边界/角点增强 Tc**：1D Dirichlet 边界 Tc > 体值；象限 > 半空间 > 全空间 | Hainzl–Roos–Seiringer, J. Spectr. Theory 12, 1507 (2022), arXiv:2201.08090；Roos–Seiringer, **Forum Math. Sigma 13:e71 (2025)**, arXiv:2308.07115 | A（元数据） |
| **BCS Tc 的线性上界判据 T_u**（"a similar linear criterion which gives a temperature T_u above which no superconductivity occurs"），并给出 T_u > T_l 与 T_u = T_l 的例子 | **Barbara Roos, J. Math. Phys. 66, 011901 (2025)**, arXiv:2407.00796 | A（摘要逐字） |
| **多带超导 Tc 恒被带间耦合抬高** | Henheik–Langmann–Lauritsen, arXiv:2409.17297（无期刊 ref，可能仍是预印本） | A（元数据）/ C（发表状态） |
| **从 reduced BCS 哈密顿量精确导出零温 BCS 泛函**（体积无关误差） | Hainzl–Panza, **arXiv:2609.04060（2026-09-03，六天前）** | A（元数据） |
| BdG 临界点 ↔ GL 临界点双向对应（不止极小元） | Frank–Hainzl–Ou Yang, arXiv:2605.26008 | A（元数据） |

### 4.3 开放边界在哪（等级 A 综述 + B 推断）

- **强耦合：完全没有。** 综述可见正文只覆盖**弱耦合渐近**与**低密度极限**；无强耦合 Tc 公式，无一般的双侧 Tc 界。**这与 KAY 线正好互补——KAY 的全部力气在 λ→∞，数学 BCS 的全部力气在 λ→0。中间没有人。**
- **多体 ODLRO**：综述逐字判定 "is out of reach of present-day mathematics"（对 eq.(2.7) 与 eq.(2.32) 各说了一次）。
- **平移对称性破缺**：§4.6 只在"a certain positivity assumption (which can be proved for a class of interaction potentials V)"下确认平移不变极小元是真极小元 ⟹ **一般 V 未解**。
- **direct/exchange 项的舍弃**：只对"very short range"平移不变系统被严格辩护（引 [7]），一般情形是启发式（引 Leggett）。
- **多带 / 各向异性**：2015 综述里没有；多带在 2024 年由 Henheik–Langmann–Lauritsen 补上（仍是瞬时核）；**各向异性/d 波的严格 Tc 界未找到**。
- **推迟核：空白**（§4.1）。

### 4.4 1 ∩ 2 的处女地——**明确指出**

**两条线的唯一接触点是 arXiv:2207.01825，Langmann, Hainzl, Seiringer, Balatsky, "Obstructions to odd-frequency superconductivity in Eliashberg theory"**（等级 A 摘要）：

> "We present a necessary condition for odd-frequency (odd-f) superconductivity (SC) to occur in a large class of materials" … "We use this condition to prove a no-go theorem ruling out the occurrence of odd-f SC in standard one-band superconductors" with "pairing interactions mediated by phonon exchange"，并有自旋涨落版的对应 no-go 定理。

即：**数学 BCS 学派的人（Hainzl、Seiringer）确实进过 Eliashberg，但问的是 odd-frequency 的 no-go，不是 Tc 界**。该文摘要与元数据中**无任何 Tc 界**。arXiv 页无期刊 ref（v1, 2022-07-05）⟹ 发表状态未知，等级 C。是否引 KAY：无法从摘要页判定（KAY 晚两年）。

**⟹ 处女地的精确陈述**：

> **不存在这样的结果**：从多体哈密顿量（或一个变分泛函）出发、带**频率依赖/推迟**配对核、**严格**导出的 Tc 上界。
>
> KAY 从 Eliashberg 方程**出发**（方程本身是未受控的近似）做严格分析；Hainzl–Seiringer 从多体哈密顿量**严格导出** BCS 泛函（Hainzl–Panza 2026 刚把 T=0 情形做干净）但核是瞬时的。**两半从未接上。**

**同时必须诚实指出**：把这两半接上是"证明 Migdal–Eliashberg 近似的可控性"，那属于**极难**（相当于要控制 Ω/E_F 展开的所有阶）。**不建议作为 P4 的主攻方向**——见 §5 的难度评估。

### 4.5 格点模型：反射正性 / ODLRO（覆盖最弱的一块，等级 B/C）

我只跑了一次普查（`abs:"reflection positivity" AND (abs:"Holstein" OR abs:"superconductivity")`），命中 3 条：

- **Freericks, arXiv:2202.02271**（Lieb 90 岁纪念文集章节）："revisiting Lieb's 1989 spin-reflection positivity proofs for the Hubbard model's ground-state spin properties, their **1995 extension to the Holstein model**, and related open questions"。⟹ Freericks–Lieb 1995（PRB 51, 2812）把自旋反射正性推到电子-声子模型。等级 A（元数据）。
- **Macris–Piguet, cond-mat/9906008, PRB 60, 13484 (1999)**：在一个特意构造成 reflection positive 的 Holstein 变体上，**严格证明强耦合低温下的 Peierls / 半填充 CDW 相变**。⟹ **是 CDW，不是超导**。等级 A（元数据）。
- Tian–Tang, cond-mat/9808207, PRB 58, 12333 (1998)：小超导颗粒的宇称参数严格下界。与 Tc 界无关。

- **Kubo & Kishi, PRB 41, 4866 (1990), "Rigorous bounds on the susceptibilities of the Hubbard model"** — 经 Semantic Scholar API 确认**存在**（作者 K. Kubo, T. Kishi；1990；PRB；42 次引用），**但该记录的 abstract 字段被出版方抹去**，我**没有拿到任何内容**。⟹ 这是文献里最常被引作"Hubbard 模型超导磁化率的严格上界"的那篇，**等级 C，内容未核实，不要写进任何账目**。

**⟹ 诚实边界：这一小块我覆盖不足。** 没有找到任何 **Holstein / 吸引 Hubbard 模型的严格 Tc 上界**。找到的严格结果是基态自旋性质（Lieb/Freericks）与 CDW 相变（Macris–Piguet），**不是超导 Tc 界**。建议单独一轮检索（见 §7）。

---

## §5　第三块 · 候选精确命题

以下三条按"推荐度"排序。每条给 (i) 数学陈述 (ii) 已知/未知边界 (iii) 与 P3 机器的接口 (iv) 难度/可交付性/AI 契合度。

---

### 【P4-A】把 P3 的谱形展开搬到 KAY 的严格算子上　★ **首选**

**(i) 精确陈述.**
设 P(dω) 为 KAY II eq.(7) 的概率测度，ϖ̄ := ⟨ω⟩/2πT，ω = ⟨ω⟩(1+ξ)，⟨ξ⟩=0，⟨ξ²⟩=σ²，⟨ξ³⟩=γ₃σ³。设 𝔨^(N)(P,T) = ρ(∫ℌ^(N)(ϖ)P(dω)) 为 KAY Theorem 4 的第 N 个本征值，v_N(ϖ̄) 为 ℌ^(N)(ϖ̄) 的最大本征向量（由 Krein–Rutman 在正锥上单一）。**待证**：存在与 P 的谱形族无关的普适数 κ_N、c_{3,N}（只依赖 ϖ̄ 与 N），使

  𝔨^(N)(P,T) = 𝔥^(N)(ϖ̄)·[1 + ½κ_N(ϖ̄)σ² + c_{3,N}(ϖ̄)γ₃σ³ + O(σ⁴)]，
  κ_N(ϖ̄) = ϖ̄²⟨v_N|ℌ^(N)''(ϖ̄)|v_N⟩ / 𝔥^(N)(ϖ̄)，
  c_{3,N}(ϖ̄) = (1/6)ϖ̄³⟨v_N|ℌ^(N)'''(ϖ̄)|v_N⟩ / 𝔥^(N)(ϖ̄)，

并由此给出 T_c^(N)(λ,P) 的矩展开与**带误差界**的余项估计。**可证伪版**：若 κ_N 依赖于 P 的四阶以上矩（而非只依赖 σ、γ₃、ϖ̄、N），命题为假。

**(ii) 已知 / 未知.**
- **已知**：KAY 给出 𝔎^(N) = ∫ℌ^(N)(ϖ)P(dω) 的显式矩阵（eq.23，元素全由 ⟦n⟧(ϖ)=ϖ²/(n²+ϖ²) 组合）；ℌ^(N) 对称；最大本征值单一（Krein–Rutman）；ℌ^(N) 对 T² 解析（eq.86 用一阶微扰论算 ∂𝔨^(N)/∂T²，**KAY 自己已经在用一阶微扰论**）。
- **未知/空白**：**KAY Remark 2（p.10）明确点出这个洞但没填**——逐字："with the sole exception of the eigenvalue 𝔨^(1)(P,T) of the first approximation, **all eigenvalues 𝔨^(N)(P,T) for N > 1 are not simply averages w.r.t. P(dω) of the pertinent largest eigenvalues 𝔥^(N)(ϖ) of ℌ^(N)(ϖ)**. The same is true for 𝔨(P,T) in relation to 𝔥(ϖ)."
  ⟹ **"𝔨^(N) 与 𝔥^(N) 的平均差多少"这个问题被作者提出、未回答。P4-A 就是定量回答它：差在 O(σ²)，系数是 κ_N。**
- 文献里**没有人**做过 Eliashberg Tc 界的谱形（矩）展开。

**(iii) 与 P3 的接口——这是最强的一条.**
- **对象完全同构**：P3 展开的 Semenok 泛函 λ_*[P] = 1/max_T ∫g(ω/2πT)P(dω) 是**标量核的 P-平均**；KAY 的 𝔨^(1) = ∫[ϖ²/(1+ϖ²)]P(dω) **也是标量核的 P-平均**，形式一模一样。P3 的 envelope-theorem 论证（g'(x_max)=0 ⟹ 鞍点位移只在 O(σ⁴) 反馈）**逐字可移植**。
- N ≥ 2 的推广只需把 envelope 定理换成 **Hellmann–Feynman / 一阶本征值微扰**——而 KAY eq.(86) 已经在用这个工具，接口是现成的。
- **P(dω) 记号完全相同**（§1.2 已核）；⟨ω²⟩ 定义与 Semenok Eq.(8) 逐位相同（§1.3 已核）⟹ σ、γ₃ 的数值可**直接复用** P3 已有的 LaH₁₀ 反演结果（r = ω₂/ω_log = 1.194/1.075/1.105/1.095 @129/163/214/264 GPa，经 r−1 = σ²−⅓γ₃σ³ 换算）。
- **产出直接可比**：P3 现在给的是 λ_*(P) 的压力依赖天花板（4.0–4.85，非严格对象）；P4-A 给的是 **T_c^(N)(λ,P) 的压力依赖严格下界 + 谱形修正**。两者可并列成一张表。
- 求解器复用：`theory/solvers/lambda_star.py`（10/10 单测）的积分/求极大骨架可直接改造；新增的只是 4×4 以内矩阵的本征值与本征向量。

**(iv) 难度 / 可交付性 / AI 契合度.**
- 难度 **中低**。核心是有限维对称矩阵的解析微扰 + 矩展开，无需新思想。真正的工作量在**余项的严格化**（要一致控制 σ⁴ 项与 N→∞ 的换序）。
- 可交付性 **高**（≥85%）。即使严格余项没做完，形式展开 + 数值验证（照 P3 用多峰/非对称谱验普适性的同一套做法）已经是可发表的增量。
- AI 契合度 **很高**：符号微分（sympy 算 ℌ^(N)''、ℌ^(N)'''）+ 高精度数值（mpmath 验普适性）+ 与 KAY 闭式 N≤4 的交叉核对，全部是 AI 的强项，且**每一步都有独立可验证的对照物**（KAY 的 eq.10–15 闭式）。
- 风险：κ_N 可能依赖 ϖ̄ 而**不是**像 P3 的 κ=5.6035 那样是纯数（P3 的普适性来自 g'(x_max)=0 这个特殊点；KAY 这边没有取极大，是解 λ=1/𝔨）。**这不构成失败**——ϖ̄-依赖的 κ_N(ϖ̄) 仍是良定义的普适函数，只是形态不同。但报告里必须一开始就把这个预期讲清楚，别承诺"又一个纯数"。

---

### 【P4-B】把 C_∞ = 0.1827262477 从浮点升级为定理（严格包围 / 闭式）　★ 次选，最稳的可交付

**(i) 精确陈述.**
𝔤(2) := ρ(𝔊(2))，𝔊(γ) 为 KAY I 构造的 ℓ²(ℕ₀) 上紧自伴算子（𝔊 = −𝔊₁+𝔊₂+𝔊₃，核分别为对角部分幂和、1/|n−m|^γ、1/(n+m+1)^γ）。C_∞ = (1/2π)√𝔤(2)。**待证 (B1)**：给出计算机辅助的**严格区间包围** 𝔤(2) ∈ [a,b]，使 (1/2π)√· 的区间宽度 < 10^{-8}，方法为「有限秩截断的严格本征值界 + 尾部迹类估计」。**待证/待证伪 (B2)**：𝔤(2) 是否有闭式（是否为已知常数的初等/超越组合）。

**(ii) 已知 / 未知.**
- **已知**：𝔤^(N)(2) 单调**上升**收敛到 𝔤(2)（KAY I）⟹ **下界方向天然免费**，任何 N 的截断都给严格下界。N≤4 有闭式（0.159154943 / 0.1796160944 / 0.1820383 / 0.1825137102）。局部收敛 ~N^{-4}（§1.8 自算）。
- **未知**：(a) **上界方向没有**——需要尾部估计（𝔊 的核是 1/k^γ 型，γ=2，Hilbert–Schmidt 范数收敛，尾部可控，**技术上完全可做**）；(b) 现有 10 位数字来源逐字是 "With the help of Maple and, independently, Mathematica … 10 significant decimal places … have stabilized"——**这是浮点收敛观察，不是包围**；(c) 闭式：**KAY I/II/III 三篇里从未提出过闭式问题**，我自算 𝔤(2) = (2π×0.1827262477)² = **1.318140**，未匹配到明显的已知常数（≠ π²/6, 4/π, 2−ln2 …）。等级 B（自算）。

**(iii) 与 P3 的接口.**
中等。**直接效果**：P3 现在引 "0.1827"（转引自 Yuzbashyan PRL）与 Semenok 的 "0.18"，属 B 级；(B1) 完成后可升为 A 级带定理的常数，且把 P3 §2.4 那条"Allen–Dynes 原文没拿到、系数只有二手"的长期挂账**彻底关闭**（因为不再需要 AD 原文——常数有了独立的严格来源）。**不直接改进 P3 的压力依赖工作。**

**(iv) 难度 / 可交付性 / AI 契合度.**
- 难度 **低到中**（B1）/ **不可预测**（B2）。B1 是标准的计算机辅助分析：区间算术 + Weyl/Bauer–Fike + 迹类尾部界。
- 可交付性 **B1 极高（>90%）**，B2 是彩票。
- AI 契合度 **极高**：这是"良定义、可完全形式化"的典型题，甚至可以 Lean 形式化（与我记忆中量子信息"带证书"扫描记的那类选题同型）。B2 可用 inverse symbolic（ISC / RIES）搜索。
- **弱点：科学分量小。** 这是把一个已知到 10 位的数升级为定理，物理上什么也没改变。**适合作为 P4-A 的附属交付物，不适合当主线。**

---

### 【P4-C】收窄 2.034 缺口　★ 真开放问题，但风险最高

**(i) 精确陈述.**
在 KAY II 的框架（P ∈ 𝒫，λ>0）内，构造显式常数 C' 与显式函数 T_c^{♯♯}(λ,P)，使 (a) T_c(λ,P) < T_c^{♯♯}(λ,P) 对所有 λ>0、P∈𝒫 成立；(b) T_c^{♯♯} ~ C'√⟨ω²⟩√λ (λ→∞) 且 **C'/C_∞ < 2.034**。
**分级里程碑**：C'/C_∞ < 1.8（弱）/ < 1.4（中）/ → 1（即证明 **Conjecture 2**，强）。
**等价的算子陈述**：把 λ→∞ 下 𝔨/⟨ϖ²⟩ 的上界从 **1+b = 5.4318** 压向真值 **𝔤(2) = 1.3181**。

**(ii) 已知 / 未知——本轮把攻击面解剖清楚了（§1.6）.**
- **不能靠调 ε**：ε=0.65 已是该族最优（我验证 ε=0.60/0.70 都更差）。[B 自算]
- **不能靠改进 Proposition 7**：其 "ϖ²+k² > k²" 这一步在 ϖ→0（即 λ→∞）时**渐近紧**，不是渐近缺口的来源。[A 原文 + B 推断]
- **缺口的两个真来源**：
  1. **Lemma 1（eq.76）** 的 ρ(AB) ≤ ρ(A)ρ(B) 解耦 + eq.(77) 把对角因子在 **n=0**（求和为空、最有利）处取值，与优化本征模脱耦；
  2. **Lemma 2（eq.81）** 的次可加 ρ(𝔎₂+𝔎₃) ≤ ρ(𝔎₂)+ρ(𝔎₃)，**丢光了 −𝔎₁ 的负贡献与 𝔎₂/𝔎₃ 之间的抵消**（真的 𝔤(2)=1.318 是含 −𝔊₁ 压低后的值）。
- **一条具体路线**：把三个算子当成一个整体，对 −𝔊₁+𝔊₂+𝔊₃ 直接做**带试验向量的 Rayleigh 上界**（如用 KAY I 已知的近似 Ξ^opt 做 Temple/Kato–Temple 型上界），而不是逐项取谱半径再相加。这需要 𝔊 的次大本征值下界——是可算的。
- **另一条**：KAY 自己的 **Conjecture 3（eq.74）** 恰好就是关于同一族 Rayleigh 商的下界猜想（"我们怀疑常序列是全局极小元"）。**证明 Conjecture 3 很可能是收窄缺口的一块拼图**，且它是纯离散不等式（有限 N 的组合命题），**AI 契合度极高、可形式化**。

**(iii) 与 P3 的接口.**
**最高的物理回报**：§1.10 已把缺口折算成 **~1090 K/√A**。收窄到 C'/C_∞ = 1.4 会把"今天可证的" 2142 K/√A 压到 1474 K/√A；证明 Conjecture 2 则给 **1053 K/√A**，与 Semenok 的头条 950 K/√A 只差 11%，且**变成定理**。这正是 P3→P4 升格的全部意义。

**(iv) 难度 / 可交付性 / AI 契合度.**
- 难度 **高**（弱里程碑中等；Conjecture 2 极高）。
- 可交付性 **中偏低**（弱里程碑 ~50%；Conjecture 2 <10%）。需要真正的新估计，不是机械化的。
- AI 契合度 **中**：数值探路（试验向量搜索、Temple 界的数值优化）AI 很强；但最后的严格论证需要人的算子理论判断。**例外**：其子问题 **Conjecture 3** 与 **"ϖ↦𝔥(ϖ) 在全 ℝ₊ 单调"**（KAY 只证到 ϖ≤√2，N=3,4 靠 "visual inspection of the plots"）**都是完全良定义、可形式化的小定理，AI 契合度很高，且都是已发表论文里公开挂着的洞**。

---

### 【不推荐】把 Hainzl–Seiringer 纲领推到推迟核

处女地是真的（§4.4），但要做的实质是"控制 Migdal 定理"，即证明 Ω/E_F 展开的可控性。这是超导数学物理公认的极难问题，无中间里程碑，与 P3 现有机器无接口。**列为地形，不列为选题。**

---

## §6　矛盾与不确定点

**U1（须主会话裁决）· J. Stat. Phys. 卷号自相矛盾。** arXiv 元数据：KAY I 记 vol.**192** art.69；KAY II/III 记 vol.**92** art.94 / art.93。P3 survey §5 记 "192, 93 & 94"。**引用前核 DOI。** 等级 C。

**U2 · 2.034 的结构恒等式是我的推算，不是原文陈述。** 原文只说 "C ≈ 2.034 C_∞"，**没有**写出 2.034 = √((1+b)/𝔤(2))。我的数值复现给 **2.0300** vs 原文 **2.034**（差 0.17%），差异落在我手算 ζ(1.65)≈2.163 的精度内。**若主会话要引用这个分解，建议先用高精度 ζ 复算一遍确认。** 等级 B。

**U3 · Corollary 3 的 N→∞ 换序。** 原文逐字 "This result also holds when N → ∞ (with superscripts ^(N) purged)"——**这是一句断言，我没有读到它的证明**（应在我未读的页里）。所以"C_∞ 是真 Tc 的渐近常数"这一条的严格性依赖该换序。**这也解释了为什么 Conjecture 2（全 λ）还是猜想。** 等级 A（断言）/ C（其证明）。

**U4 · KAY 全套 μ\* = 0。** 框架里没有库仑赝势。与实验/Semenok 表对照时，Semenok 的 0.3175 ω_max 也是 μ\*=0 的数值解，可比；但**任何要落到真实材料的 Tc 界都缺 μ\* 这一层，而 P3 已知 μ\*(ω_c) 重整化是 D2 方向的开放问题**。这是 P4 的一条系统性限制，必须写进任何对外陈述。

**U5 · λ_* 与 Tc 界是逻辑独立的两件事，但被 Semenok 复合使用。** λ_*（动理学稳定性，非定理，§2）与 T_c(λ) 上界（KAY，有定理）没有共同的证明基础；Semenok Eq.(9) 的 "T_c < 0.32 ω_max" 是把两者串联（λ ≤ λ_* 代入 Tc(λ)）。**P4 若要给"严格的 Tc 上界"，λ 那一侧仍然不严格。** 诚实的产品形态是：**"给定 λ 的严格 Tc 上界" + "λ 的物理上界（非严格）"，两层分开写。** 这一点必须前置声明，否则会被审稿人一击致命。

**U6 · 一个未解决的悬案：Allen–Dynes 拟合公式 vs 真渐近。** AD 拟合公式的 λ→∞ 极限（μ\*=0）是 0.1878，比真常数 0.1827 高 2.78%；而 AD 自己在 eq.(22) 写的渐近常数是 0.182。**同一篇论文里的两个数不自洽（0.182 vs 0.1878）**，因为 f₁/f₂ 是拟合、eq.(22) 是渐近计算。P3 survey §2.4 已挂账"AD 原文未取"。**本轮没有解决——我仍未取到 PRB 12, 905 原文**，0.1878 是我从 Xie 2022 转录的公式自算（B 级），0.182 是 KAY 转述（B 级）。

**U7 · 格点模型那一小块覆盖不足（§4.5）。** Kubo–Kishi 1990 只确认了存在与标题，**内容零信息**。没有找到任何 Holstein / 吸引 Hubbard 的严格 **Tc 上界**。不要基于本报告对这一块下任何结论。

**U8 · 2207.01825 的技术细节未核。** 是否用完整 Eliashberg 方程、是否发表、是否引 KAY——三项全未确认（只取到 arXiv 摘要页）。等级 C。

**U9 · math-ph 普查的覆盖限制。** "除 KAY 外无人给严格 Eliashberg 常数"这一判断基于单一过滤 `cat:math-ph AND abs:"Eliashberg"`；未交叉分类到 math-ph 的论文会漏。**Combescot PRB 51, 11625 (1995) 未取原文**，是最可能藏着竞争性强耦合常数的一篇。

---

## §7　建议的后续检索方向（按优先级）

1. **Combescot, "Strong-coupling limit of Eliashberg theory", PRB 51, 11625 (1995)**——无 arXiv，需订阅。KAY III 把它列为先行工作。**是 U9 的唯一未堵住的口子。**
2. **KAY I（2409.00533）的 Propositions 5 与 6 全文**——ζ(1+ε)ζ(5−ε) 的产地，是 P4-C 的直接攻击面。本轮只取到 HTML 摘要级，**建议走 PDF→Read 路线**。
3. **KAY II 的 pp.19–26 与 34–44**——Corollary 3 的 N→∞ 换序证明（U3）、Theorem 7 的证明。本轮读了 pp.1–18、27–33。
4. **Allen–Dynes 1975 PRB 12, 905 原文的 eq.(22)**——长期挂账（P3 survey §5 第 1 项）。若 P4-B 成功，这条挂账可作废。
5. **格点模型严格 Tc 界**（U7）：Kubo–Kishi 1990 全文；检索 "attractive Hubbard" + "off-diagonal long-range order" + rigorous；Miyao 关于 Holstein–Hubbard 的系列。
6. **arXiv:2207.01825 全文**（U8）+ 查其发表状态。
7. **Henheik–Langmann–Lauritsen 2409.17297 是否已发表**——多带是 P3 各向异性讨论的潜在接口。

---

## §8　本轮新增/确认的编号

| 编号 | 论文 | 状态 |
|---|---|---|
| **2409.00532**v3 | KAY **II: Dispersive phonons**, J. Stat. Phys. 2025 | **新取；PDF pp.1–18, 27–33 逐页读**（P3 survey §5 挂账项，已关闭大半） |
| **2409.00533**v2 | KAY **I: The γ-model**, J. Stat. Phys. 192, art.69 | **新增**（P3 survey 未记此篇）；HTML 已读 |
| **2409.02121**v2 | KAY **III: Einstein phonons**, J. Stat. Phys. 2025 | **新增**（P3 survey 未记此篇）；HTML 已读 |
| **2206.07575** | Yuzbashyan–Kiessling–Altshuler, PRB 106, 064502 (2022)，γ 模型极端推迟区 | **新增**；未取 |
| **2207.01825** | Langmann–Hainzl–Seiringer–Balatsky, odd-frequency no-go in Eliashberg | **新增**；**1∩2 的唯一接触点**；摘要级 |
| **1511.01995** | Hainzl–Seiringer, J. Math. Phys. 57, 021101 (2016)，BCS 泛函综述 | **新增**；ar5iv 全文；**判定纲领为纯瞬时核** |
| **0704.3564** | Frank–Hainzl–Naboko–Seiringer, J. Geom. Anal. 17, 559 (2007)，弱耦合 Tc | **新增** |
| **0801.4159** | Hainzl–Seiringer, PRB 77, 184517 (2008)，Tc 与能隙 | **新增** |
| **1102.4001** | Frank–Hainzl–Seiringer–Solovej, JAMS 25, 667 (2012)，BCS→GL | **新增** |
| **2407.00796** | Roos, J. Math. Phys. 66, 011901 (2025)，**BCS Tc 的线性上界判据** | **新增**；摘要级 |
| **2409.17297** | Henheik–Langmann–Lauritsen，多带 Tc 增强 | **新增**；元数据级 |
| **2308.07115** | Roos–Seiringer, Forum Math. Sigma 13:e71 (2025)，角点增强 | **新增** |
| **2201.08090** | Hainzl–Roos–Seiringer, J. Spectr. Theory 12, 1507 (2022)，边界超导 | **新增** |
| **2609.04060** | Hainzl–Panza（2026-09-03），零温 BCS 泛函的精确导出 | **新增**；最新 |
| **2210.09356 / 2105.05623 / 1706.05686 / 2605.26008** | Deuchert–Hainzl–Maier ×2、Frank–Hainzl–Langmann、Frank–Hainzl–Ou Yang | **新增**；元数据级 |
| **2202.02271** | Freericks，Lieb 自旋反射正性（含 1995 Holstein 推广）综述章节 | **新增**；元数据级 |
| **cond-mat/9906008** | Macris–Piguet, PRB 60, 13484 (1999)，Holstein 变体的严格 Peierls/CDW | **新增**；元数据级；**是 CDW 不是超导** |
| PRB 41, 4866 (1990) | Kubo–Kishi，Hubbard 磁化率严格界 | 经 S2 API 确认存在；**内容零信息，等级 C** |
| 2409.19562 | Yuzbashyan–Altshuler–Patra, PRL 135, 026503 (2025) | 已有；本轮补读 pp.1–3，判定"非定理" |
