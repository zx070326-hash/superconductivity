# P4 补检三：KAY III（arXiv:2409.02121v2）全文 + KAY I Appendix A 的 HS 范数界

检索 agent：researcher　日期：2026-09-10
任务：派单 (b) = 读 KAY III；顺带 leftover (a) = 取 KAY I App.A 的 ‖𝔊_j(2)‖_HS 一手表达。

**一手材料**
- arXiv:**2409.02121v2**（"Revised version of April 12, 2025"，正文 29 页 + Appendix A 至 p.34 + References 至 p.36）**全文逐页读**（PDF 落盘 + Read）。
- arXiv:**2409.00533v2**（KAY I）**pp.41–45 逐页读**（Appendix A「Compactness of 𝔊」+ Appendix B 开头）。

**证据等级**：**A** = 一手逐字 / **B** = 对 A 级的确定性推算（标 `[自算]`）/ **C** = 推断或未取到

**三条结论速览**
| 目标 | 裁决 | 等级 |
|---|---|---|
| 1. KAY III 有无 Einstein 极限的**定量谱隙**？ | **无。全文零处。**§2.2a 量化部分维持「首个」定位。 | A（未见，全文核） |
| 2. KAY III 有无一般-P 的 𝔨^(N)↑𝔨 论证？ | **无**（它连 P=δ 的收敛都只是散文）。**且 KAY 自陈"don't know how fast"**——直接为增量 (b) 背书。 | A |
| 3. leftover (a)：‖𝔊_j(2)‖_HS | **𝔊₃ 有闭式 ζ(1+γ)（A 级）；𝔊₂ 可由原文成分组装（A+B）；𝔊₁ 原文无显式界（A 级未见）。** | A / A+B / A-未见 |

---

## §1 目标一【最重要】：KAY III 有无 Einstein 极限的定量谱隙？

### 1.1 裁决：**无。全文 36 页零处。** [A 级，逐页核]

KAY III 从未提及第二本征值、谱隙、本征值间距，也没有任何形式的 λ₁−λ₂ 定量陈述。**关于 ℌ(ϖ) 的谱，全文只给两类信息：顶本征值 𝔥(ϖ) 本身，以及顶本征值在 ϖ→0 / ϖ→∞ 两端的渐近展开。**

**⟹ 我方 §2.2a 的量化谱隙部分维持「首个」定位，不需降级为「Einstein 特例的推广」。这是本轮最重要的结论。**

### 1.2 逐条排除（A 级）

**(a) Proposition 2（p.20）——只给「顶本征值 = 谱半径」，不给单一性**
> **Proposition 2**: *Let ϖ > 0 be given. Then the largest eigenvalue 𝔥(ϖ) of ℌ(ϖ) is also the spectral radius ρ(ℌ(ϖ)).*
>
> 其证明一句话（p.20）："The proof of Proposition 2 is implied by the proof of the analogous statement about 𝔨(P,T) in [KAYb]."

**(b) p.20 的 Krein–Rutman 段——逐字复制 KAY II p.20，含同一瑕疵**

KAY III p.20 给出 eqs.(44)–(46)，与 KAY II eqs.(44)–(46) 结构完全相同：
> (44) (ℑ − L_E(ϖ)ℌ)Ξ^opt = 0；(45) ℭ(𝔥(ϖ))Ξ^opt = Ξ^opt；(46) **ℭ(η) := (ηℑ + ℌ₁)^{-1}(ℌ₂ + ℌ₃)**

然后逐字：
> "As in the proof of Theorem 1 in [KAYa] one shows that ℭ(η) for η > 0 is a compact operator that maps the positive cone ℓ²_{≥0}(ℕ₀) into itself, **in fact mapping any non-zero element of ℓ²_{≥0}(ℕ₀) into the interior of ℓ²_{≥0}(ℕ₀)**, and that the spectral radius of ℭ(𝔥) equals 1. Thus the Krein–Rutman theorem applies and guarantees that the nontrivial solution Ξ^opt of (45) is in the positive cone ℓ²_{≥0}(ℕ₀) (after at most choosing the overall sign) …"

**⟹ 重要副产物：我方在 v0.3 记录的「into the interior」瑕疵（ℓ²_{≥0} 在 ℓ²-范数下内部为空）在 KAY 三篇中全部出现，是同一段文字的三次复制，而非 KAY II 的孤立笔误。** 这提高了该观察的分量，但也意味着**三篇都不能作为单一性的可引用来源**。[A 级逐字 + B 级判断]

注意 KAY III 此处**连 "non-degenerate" 一词都没有**（KAY II 至少在 p.25 的 γ-模型段落里写了）。**⟹ Einstein 侧的单一性在 KAY III 中比 KAY II 更弱，不是更强。**

**(c) Theorem 6（p.11，eqs.22–24）——ϖ→∞ 的顶本征值渐近，不是谱隙**
> **Theorem 6**: *The eigenvalues 𝔥^(N)(ϖ) are analytic about ϖ = ∞, with*
> 𝔥^(N)(ϖ) = 𝔥^(N)(∞) − B_N·(1/ϖ²) + O(1/ϖ⁴)　(22)
> *where* 𝔥^(N)(∞) = −1 + 2Σ_{n=0}^{N−1} 1/(2n+1)　(23)
> *and* B_N = [ Σ_{n=0}^{N−1}Σ_{m=0}^{N−1} ((n−m)²+(n+m+1)²)/((2n+1)(2m+1)) − Σ_{ℓ=0}^{N−1} (2/(2ℓ+1)²)Σ_{k=1}^{ℓ}k² ] / [ Σ_{j=0}^{N−1} 1/(2j+1) ] ≥ 1　(24)

**只有顶本征值一条曲线的展开系数**，无任何次本征值信息。

**(d) eq.(52)（p.22）——秩一极限矩阵，与 KAY II eq.(51) 同型**
> lim_{ϖ→∞} ϖ²(ℌ^(N)(ϖ) + ℑ^(N) − 2Ξ*_N ⊗ Ξ*_N) = 𝔏^(N)　(52)

与 KAY II 一样，**这是 ϖ→∞（等价 T→0）极限矩阵的结构，不是有限 ϖ 的谱隙**。KAY III 甚至没有像 KAY II p.23 那样把极限矩阵的完整谱（−1 的 N−1 重 + 顶值）写出来。

**(e) ⚠ 全文唯一的 "gap" 一词是陷阱（p.10 逐字）**
> "For then **the gap between the T_c^(3) and T_c^(4) curves** is so small that the line width of the plotted curves fills it."

**这是 Fig.1 中两条相继下界曲线的视觉间距，与算子谱隙毫无关系。检索关键词时务必不要误采。**（另 p.10 "when λ is less than ≈ 0.5, then the gap between these two curves becomes clearly visible" 同理。）

### 1.3 一个可用的次级素材（A 级，但非谱隙）

KAY III **Appendix A（pp.30–34）给出了 ℌ^(2)、ℌ^(3)、ℌ^(4) 全部矩阵不变量的 ϖ-有理函数闭式**（eqs.66–151，含 tr ℌ^(2)、det ℌ^(2)、tr ℌ^(3)、tr adj ℌ^(3)、det ℌ^(3)、tr ℌ^(4)、tr(ℌ^(4))²、tr(ℌ^(4))³、det ℌ^(4)，系数如 P_j / Q_j / R_j / A_j / B_j / C_j / D_j 全部为显式整数）。例如：

> tr ℌ^(2) = (2/3)·ϖ²(5+ϖ²)/((1+ϖ²)(3²+ϖ²))　(66)
> det ℌ^(2) = −(1/3)·ϖ⁴(497+357ϖ²+81ϖ⁴+5ϖ⁶)/((1+ϖ²)²(2²+ϖ²)²(3²+ϖ²))　(67)
> tr ℌ^(3) = (1/15)·(4²+ϖ²)/∏_{j=1}^{5}(j²+ϖ²) · Σ_{j=0}^{3}P_j ϖ^{2j+2}，P₀=−1642, P₁=−1123, P₂=−56, P₃=1　(68)–(72)

**⟹ 若我方需要在 Einstein 特例上做 N≤4 的谱隙数值对照（λ₁ 与 λ₂ 都可由这些不变量解出），这是现成的一手闭式来源。** 但**这是我方的用法，KAY 自己没有用它算过任何谱隙**。[A 级（公式）+ B 级（用途判断）]

---

## §2 目标二：KAY III 是否含一般-P 的收敛论证？

### 2.1 裁决：**不含，且它连 Einstein 特例的收敛都只是散文。** [A 级]

**(a) Theorem 3（p.6）——与 KAY II Theorem 4 完全同型：只有不等式**
> **Theorem 3**: *For all N ∈ ℕ, L_E(ϖ) < 1/𝔥^(N)(ϖ), where 𝔥^(N)(ϖ) is the largest eigenvalue of ℌ^(N)(ϖ), the restriction of ℌ(ϖ) to the first N components of ℓ²(ℕ₀). The eigenvalues 𝔥^(N)(ϖ) can be explicitly computed for N ∈ {1,2,3,4}.*

无 "converge"、无 "N→∞"。**与 KAY II Theorem 4 逐字平行。**

**(b) §3.2.1（p.21）——KAY II p.21 那段散文的近乎逐字复制**
> "Since ℌ(ϖ) is compact, in principle one can get arbitrarily accurate upper approximations to L_E(ϖ) by restricting ℌ(ϖ) to suitably chosen finite-dimensional subspaces of ℓ²(ℕ₀). A sequence of decreasing rigorous upper bounds on L_E(ϖ) that converges to L_E(ϖ) is obtained by restricting the variational principle to a sequence of subspaces of ℓ²(ℕ₀) of vectors of the type **Ξ_N := (ξ₀, ξ₁, …, ξ_{N−1}, 0, 0, …)**, with ξ_j > 0 for j ∈ {0,…,N−1} and N ∈ ℕ. Evaluating (48) with Ξ_N in place of Ξ^opt yields a strictly monotonically decreasing sequence of upper bounds L_E^(N)(ϖ) on L_E(ϖ)"

**同样的散文、同样的「紧算子可被有限秩逼近」一句理由、同样没有证明。** §4.1 Summary（p.28）也重复一次同样的措辞。

**(c) 主动确认：KAY III 明说自己是 KAY II 的特例，不含一般 P**
> p.5（§2 开头）："Although most, though not all of our results in the present paper are special cases of the results we proved in [KAYb], we state them as theorems or propositions in their own right, rather than as corollaries."
> p.23（§3.3）："Turning to Theorem 4, it suffices to note that its proof is included as limiting case P(dω) → δ(ω−Ω)dω with Ω > 0 of the proof of Theorem 6 in [KAYb]."
> p.27（§4.1）："The results obtained in the present paper are **mostly special cases of the results of [KAYb]**."

**⟹ 隔离逻辑成立：KAY III 是 P=δ 单点测度特例，动不了我方推论(7) 的 (a) 一般-P 增量。** [A 级]

### 2.2 【本轮最有价值的发现】KAY 自陈「不知道收敛有多快」

**p.9 逐字（§2，Corollary 1 之后）**：
> "Besides knowing the best explicitly available 'sandwiching bounds' on T_c(λ,Ω), it is of interest to plot all these bounds in a single diagram to get some visual impression of the speed of convergence; recall that **we only know that the sequence L_E^(N)(ϖ) converges downward to L_E(ϖ) when N → ∞, but don't know how fast** — analogously for the lower critical-temperature bounds T_c^(N)(λ,Ω) when λ > max{λ_N, 0.7671}…"

**⟹ 这是 KAY 本人在最简单的 Einstein 特例上、以正面陈述的方式承认没有收敛速率。** 我方增量 (b)（定量速率 ρ_N = 2‖(I−P_N)ℌ‖_HS）由此从「未见即新」升级为**「作者自陈缺失」**——这是强得多的证据等级，建议在论文中直接引这一句。[A 级]

### 2.3 ⚠ 一个必须避开的速率陷阱（A 级，p.10）

> "On the other hand, when λ is less than ≈ 0.5, then the gap between these two curves becomes clearly visible. In fact, near the λ axis convergence is slow. **The sequence of upper bounds L_E^(N)(ϖ) to L_E(ϖ) meets the λ axis at explicitly computable locations λ_N that converge slowly to 0 like 1/ln N as N → ∞.**"

**这个 1/ln N 是 λ_N（第 N 个上界曲线与 λ 轴的交点）的收敛速率，不是 𝔥^(N)(ϖ) → 𝔥(ϖ) 在固定 ϖ 处的收敛速率。** 二者不可混。（其来源是 eq.(23) 的 𝔥^(N)(∞) ~ ln N 发散，即 §1.2(c)。）**若我方论文提到「KAY 唯一给出的 N-速率」，必须写清它是 λ_N 的、不是本征值的。**

### 2.4 KAY III 明确是系列终篇（A 级，p.29）

> "The present paper **completes our series of three papers** on the rigorous study of the Eliashberg gap equations as linearized about the normal state. … Our next goal is the study of the non-linear Eliashberg gap equations. … However, in general our Hilbert space analysis of papers I–III **will have to be replaced by an analysis of operators in a certain Banach space**."

**⟹ 线性化侧不会有 KAY IV，本轮之后不存在「被后续论文抢先」的剩余风险。** 记忆中「无 KAY IV」的旧判断由此获得一手确认。

---

## §3 leftover (a)：KAY I Appendix A 的 HS 范数界

KAY I（2409.00533v2）**Appendix A「Compactness of 𝔊」pp.42–43**，含 Proposition A / B / C 各配完整证明。逐条抄录：

### 3.1 Proposition A（𝔊₁）——**只证有限，无显式界** [A 级]

> **Proposition A**: *For any γ > 0, 𝔊₁ ∈ ℓ²(ℕ₀) × ℓ²(ℕ₀).*
> *Proof*: (𝔊₁)_{n,m} = [ (1/(2n+1))Σ_{k=1}^{n}(2/k^γ) ]δ_{n,m}　(108)
> Σ_nΣ_m (𝔊₁)²_{n,m} = Σ_{n∈ℕ₀} (1/(2n+1)²)[Σ_{k=1}^{n} 2/k^γ]²　(109)
> Σ_{k=1}^{n} 1/k^γ ≤ 1 + { (1/(1−γ))(n^{1−γ}−1) ; γ≠1　|　ln n ; γ=1 }　(110)
> **"Hence Σ_nΣ_m (𝔊₁)²_{n,m} < ∞. Q.E.D."**

**⚠ 证明到「< ∞」为止，不给任何数值或闭式上界。** [A 级]

**我方若需要，可由 (110) 自行组装** `[自算，B 级]`：γ=2 时 Σ_{k=1}^{n}1/k² ≤ 1 + (n^{−1}−1)/(−1) = 2 − 1/n < 2，故 [Σ2/k^γ]² < 16，得 ‖𝔊₁(2)‖²_HS < 16·Σ_{n≥1}1/(2n+1)² = 16(π²/8 − 1) ≈ 3.739，即 ‖𝔊₁(2)‖_HS < 1.934。（n=0 项为 0。）**此数为本报告自算，非 KAY 原文，对外引用须标注。**

### 3.2 Proposition B（𝔊₂）——**原文给出全部成分但未组装成一个数** [A 级成分 + B 级组装]

> **Proposition B**: *For any γ > 0, 𝔊₂ ∈ ℓ²(ℕ₀) × ℓ²(ℕ₀).*
> (𝔊₂)_{n,m} = (1/√(2n+1))·((1−δ_{n,m})/|n−m|^γ)·(1/√(2m+1))　(111)
> Σ_nΣ_m (𝔊₂)²_{n,m} = **2 Σ_{n≥0} Σ_{m>n} (1/(2n+1))·(1/(m−n)^{2γ})·(1/(2m+1))**　(113)
> Σ_{m>n} (1/(m−n)^{2γ})(1/(2m+1)) = Σ_{m>0} (1/m^{2γ})(1/(2m+2n+1))　(114)

**γ > 1/2 情形，原文逐字（p.43）**：
> "Now suppose first that γ > 1/2. We simply estimate 2m+2n+1 > 2n+1 for m > 0. **The sum over m in (114) then equals ζ(2γ), multiplied by 1/(2n+1). The resulting sum over n at r.h.s.(113) equals (1/8)π².**"

**⟹ 组装** `[自算，B 级]`：‖𝔊₂(γ)‖²_HS ≤ 2·ζ(2γ)·(π²/8) = **(π²/4)·ζ(2γ)**。
γ=2：ζ(4)=π⁴/90 ⟹ ‖𝔊₂(2)‖²_HS ≤ π⁶/360 ≈ 2.6705 ⟹ **‖𝔊₂(2)‖_HS ≤ 1.6342**。
（自洽性核对：survey 记 ρ(𝔊₂(2)) = 0.711 ≤ 1.634 ✓，HS ≥ 谱半径成立。）

（γ ≤ 1/2 情形见 eqs.115–118，结果为 (1/2^{1−γ})(1−2^{−(1+γ)})ζ²(1+γ)×2；γ=2 用不到。）

### 3.3 Proposition C（𝔊₃）——**有干净闭式，可直接引用** [A 级]

> **Proposition C**: *For any γ > 0, 𝔊₃ ∈ ℓ²(ℕ₀) × ℓ²(ℕ₀).*
> (𝔊₃)_{n,m} = (1/√(2n+1))·(1/(n+m+1)^γ)·(1/√(2m+1))　(119)
> Σ_nΣ_m (𝔊₃)²_{n,m} = Σ_nΣ_m (1/(2n+1))·(1/(n+m+1)^{2γ})·(1/(2m+1))　(120)
> 1/(n+m+1)^{2γ} = (1/(n+m+1)^γ)(1/(n+m+1)^γ) ≤ (1/(n+1)^γ)(1/(m+1)^γ)　(121)–(122)
> **Σ_nΣ_m (𝔊₃)²_{n,m} ≤ [ Σ_n 1/(n+1)^{1+γ} ]² = ζ²(1+γ)**　(123)

**⟹ ‖𝔊₃(γ)‖_HS ≤ ζ(1+γ)　——原文 eq.(123) 的直接改写，A 级，可对外引用。**
γ=2：**‖𝔊₃(2)‖_HS ≤ ζ(3) ≈ 1.2020569**。
（自洽性核对：ρ(𝔊₃(2)) = 1.025 ≤ 1.202 ✓，且相当紧。）

### 3.4 ⚠ 对 ρ_N 认证的实际可用性判断 [B 级]

**KAY I App.A 给的是三个算子各自的「全 HS 范数」有限性/上界，没有任何尾部（tail）估计。** 我方 ρ_N = 2‖(I−P_N)ℌ‖_HS 需要的是 **Σ_{(n,m)∉[0,N)²}** 的部分和，App.A 完全没有涉及。

**⟹ implementer 自算的 HS 尾仍是净新增，App.A 只能提供两件事**：
1. **可对外引用的一手锚点**：𝔊₃ 的 ζ(1+γ) 闭式（A 级），以及 𝔊₂ 的 (π²/4)ζ(2γ)（成分 A 级 + 组装 B 级）。
2. **一个反面事实**：**𝔊₁ 没有一手显式界**——若我方 ρ_N 的推导链条里出现 ‖𝔊₁‖_HS 或 ‖ℌ₁‖_HS，必须自证并标 B 级，不能挂 KAY I。

---

## §4 顺带确认与更正

### 4.1 U1（期刊卷号）——**KAY I 部分获一手解决** [A 级]

KAY III 的 References（p.35）逐字：
> [KAYa] Kiessling, M. K.-H., Altshuler, B. L., and Yuzbashyan, E. A., *Bounds on T_c in the Eliashberg theory of superconductivity. I: The γ model*, **J. Statist. Phys. 192, art.69, 35.pp (2025)**.
> [KAYb] Kiessling, M. K.-H., Altshuler, B. L., and Yuzbashyan, E. A., *Bounds on T_c in the Eliashberg theory of superconductivity. II: Dispersive phonons*, **44p., J. Statist. Phys. (in press, 2025)**.

**⟹ KAY I = J. Stat. Phys. 192, art.69 (2025)，由 KAY 本人的参考文献确认，A 级，U1 的 KAY I 部分关闭。**
**⟹ KAY II 在 KAY III 定稿（2025-04-12）时仍 "in press"，故 KAY III 无法确认其卷号。**

**但由此可作一个判断** `[B 级]`：KAY I 是 vol **192**；同刊同年的 KAY II 记为 "vol.92" 高度可疑，**大概率是 "192" 漏首位的转录滑落**（art. 号 94 与 KAY I 的 69 同量级，与同卷相容）。**引用前仍须核 DOI，但应优先假设 192 而非 92。** 我方记忆与 survey 的旧记录需据此标注。

### 4.2 KAY III 的 Remark 1（p.13）——关于 AD 常数的补充史料 [A 级]

> **Remark 1**: "We remark that the asymptotic behavior T_c(λ,Ω) ∼ CΩ√λ for large enough λ, given Ω, was anticipated in [AD], though based on nonrigorous arguments; see their eq.(22). Incidentally, in [AD] **the value C = 0.182 was stated to be computed with a 64 Matsubara mode approximation** to the linearized Eliashberg gap equation. Their computation was not very accurate, for **the rounded approximation to three significant digits reads C = 0.183**, and 0.183 is already the rounded value of the closed form approximation with merely four Matsubara modes obtained in our paper."

**⟹ 与我方记忆一致（AD 的 0.182 是 64-mode 数值），并新增一条：真值三位舍入为 0.183，而 KAY 的 N=4 闭式已达 0.183。** 可用于强化「AD 数值不精确」的叙述。

### 4.3 KAY III 的其它可复用常数 [A 级]

- **λ_* < 0.7670…**（Theorem 5, p.8；另 p.26 给 1/𝔥^(4)(√2) = 0.7670…）——Einstein 情形的存在性阈值，比 KAY II 的一般-P 版本 λ_*(P) ≤ (3/2)Ω̄²/⟨ω²⟩ 具体。
- **λ₂ = 3/5 = 0.6**（eq.64, p.25）、**λ₃ = 15/31 = 0.48387…**、**λ₄ = 105/247 = 0.4251…**（p.26）。
- **1/𝔥^(3)(√2) = 1.0158…**（p.26）。
- Corollary 1（p.9，eq.18）：**T_c^(1)(λ,Ω) = (Ω/2π)√(λ−1)**，λ>1。
- Proposition 1（p.14，eq.29）：lim_{λ→∞}(1/√λ)f(λ) = (1/g)T_c(g,2)，数值 0.1827262477…（"a 200 mode approximation yields"）。

---

## §5 对 P4 文档的净影响

| 我方部件 | 是否需改 | 理由 |
|---|---|---|
| **推论(7)（一般-P 收敛）** | **不需改** | KAY III 是 P=δ 特例且连该特例都未证；隔离逻辑成立 |
| **增量 (a)（推广到一般 P）** | **不需改，可加强** | KAY III 三处自陈「mostly special cases of [KAYb]」 |
| **增量 (b)（定量速率）** | **建议加强措辞** | 从「未见」升级为 KAY 自陈 "**but don't know how fast**"（KAY III p.9），可直接引用 |
| **§2.2a 量化谱隙「首个」定位** | **不需降级** | KAY III 全文无任何谱隙陈述 |
| **§2.2b（K–R 瑕疵标注）** | **建议加一句** | 该「into the interior」措辞在 KAY I/II/III **三篇复制出现**，非孤立笔误；且 KAY III 连 "non-degenerate" 都未写 |
| **ρ_N 的对外引用** | **可补一手锚点** | 𝔊₃：ζ(1+γ)（A 级，KAY I eq.123）；𝔊₂：(π²/4)ζ(2γ)（成分 A + 组装 B）；**𝔊₁：无，须自证并标 B 级** |
| **KAY I 卷号** | **可改为 A 级** | J. Stat. Phys. **192**, art.69 (2025)，KAY III 参考文献确认 |
| **KAY II 卷号** | **仍存疑但改默认假设** | "92" 大概率是 "192" 的转录滑落 `[B]`；仍须核 DOI |
| **「无 KAY IV」** | **可升为 A 级** | KAY III p.29 "completes our series of three papers" |

---

## §6 遗留

- **KAY II 的期刊卷号仍未一手证实**（三篇 PDF 内均无）。唯一出处仍是 arXiv comments 字段。若要 A 级，需查 DOI 或 J. Stat. Phys. 官网。
- **de Pagter / ideal-irreducible Krein–Rutman 的适用性**：本轮未涉及，维持「未验证提示」状态（团队 v0.3 已不依赖它）。
- **KAY III Appendix A 的 ℌ^(N)(ϖ) 不变量闭式（eqs.66–151）本轮已读但未逐字转录全部整数系数**（A_j / B_j / C_j / D_j 等共约 60 个大整数）。若我方要做 N≤4 的 Einstein 谱隙数值对照，需回 PDF pp.30–34 逐个抄；本报告已转录 eqs.(66)(67)(68) 与 P_j 作为样本，其余按需再取。
