# P4 补检二：KAY II（arXiv:2409.00532v3）关于「有限-N 截断 → N→∞ 真算子」的确切陈述

检索 agent：researcher　日期：2026-09-10
目的：确认把 Eliashberg 线性稳定算子的有限-N 截断展开抬到 N→∞ 真算子时，哪些环节可以直接引 KAY，哪些必须我方自证。

**一手材料**：arXiv:2409.00532**v3**（"Revised version of July 12, 2025"），PDF 全 42 页正文逐页读（webfetch 落盘 + Read）。
本报告所有页码 = **PDF 内印刷页码**（与 PDF 物理页一致，标题页 = p.1）。

**证据等级**
- **A** = 一手 PDF 逐字（原文语句/公式号/页码）
- **B** = 对 A 级做的确定性推理（标 `[本报告自算/自推]`）
- **C** = 推断 / 未取到

**先说结论（一句话）**：**一般 P(dω) 下的 𝔨^(N) ↑ 𝔨 在 KAY II 中被明确断言并被当作已知事实使用，但它不是任何编号定理/推论的组成部分，也没有独立证明——只有一句「紧算子可被有限秩逼近」的散文式理由。想严格引用，必须自证。**

---

## §0 KAY II 的结果清单（A 级，用于定位与消歧）

| 编号 | 页 | 内容一句话 |
|---|---|---|
| Conjecture 1 | 7 | 存在 T_c，相变连续（需非线性能隙方程） |
| Theorem 1 | 7 | (λ,P,T) 正锥分两个单连通区，临界超曲面 𝒮_c 是 Λ(P,T) 的图 |
| Theorem 2 | 8 | T ↦ Λ(P,T) 在 [T_*,∞) 严格单调增，T_*(P) ≤ Ω̄(P)/2√2π；eq.(8) |
| Remark 1 | 8 | 复述 KAY I 的 𝔤^(N)↑𝔤 与 N>200 数值 0.1827262477… |
| Theorem 3 | 8 | Λ(P,T) = 1/𝔨(P,T)，eq.(9) |
| **Theorem 4** | **9** | **∀N: Λ(P,T) < 1/𝔨^(N)(P,T)；𝔨^(N) 定义 + N=1..4 闭式 eqs.(10)–(23)** |
| Remark 2 | 10 | **警告：N>1 时 𝔨^(N)(P,T) ≠ 𝔥^(N)(ϖ) 的 P-平均** |
| Proposition 1 | 11 | 𝔨^(N)(P,T) 在 T∈[T_*,∞) 严格单调减 |
| Theorem 5 | 11 | T→0 极限 𝔨₀^(N) = −1+2Σ_{n=0}^{N−1}1/(2n+1)，eq.(24)，P-独立，~ln N |
| Corollary 1 | 11 | λ_*(P) < 1/𝔨^(4)(P,T_*(P))，eq.(25) |
| Theorem 6 | 11 | 𝔨 ≤ 𝔨*，eq.(26)（含 ε=0.65 的 ζ 因子） |
| Corollary 2 | 12 | Jensen ⟹ 更弱的 𝔨^♯，eq.(27) |
| Theorem 7 | 13 | 𝔨^(N)(P,T) 在 T=∞ 附近解析，eq.(28)–(29) |
| **Corollary 3** | **13** | **T_c^(N)(λ,P) 的大-λ 渐近闭式 eq.(30)；末句「N→∞ 亦成立」** |
| Conjecture 2 | 14 | T_c(λ,P) < (1/2π)√(𝔤(2)⟨ω²⟩λ)，eq.(31) |
| Remark 3 | 14 | 该渐近「first stated in [AD], see their eq.(22)」 |
| Corollary 4 | 14 | 显式上界 T_c^♯，eq.(32)，b := 2((2^{1+ε}−1)ζ(1+ε)ζ(5−ε))^{1/2}, ε=0.65；下界 T_c^♭ eq.(33) |
| Theorem 1⁺ | 18 | Theorem 1 的自旋链版；Q(Ξ) 在球面上取到最小；Ξ^opt ∈ ℝ_+^{ℕ₀} |
| **Proposition 2** | **19** | **𝔎_j ∈ ℓ²(ℕ₀×ℕ₀)，即 HS，故紧** |
| Proposition 3 | 20 | 最大本征值 𝔨(P,T) = 谱半径 ρ(𝔎(P,T)) |
| Proposition 4 | 24 | ⟨𝔊^(N)(γ')⟩_γ > 0，eq.(56) |
| Proposition 5 | 24 | n ↦ (1/√(2n+1))(Ξ_N^opt(γ))_n 正且递减 |
| Remark 4–6 | 28 | 重排不等式（Hardy–Littlewood–Pólya / Pruss / Riesz） |
| Remark 7 + Conjecture 3 | 29 | Rayleigh 商下界猜想 eq.(74) |
| Remark 8 | 30 | 线性稳定 ⟺ ℭ(1/λ) 是 ℓ²_{≥0} 上的压缩映射 |
| Lemma 1 | 30 | ρ(ℭ(η)) ≤ ρ((ηℑ+𝔎₁)^{-1}) ρ(𝔎₂+𝔎₃)，eq.(76) |
| Proposition 6 / Corollary 5 | 31 | Λ ≥ 1/ρ(𝔎₂+𝔎₃)，eq.(80) |
| Lemma 2 / Corollary 6 | 31 | 次可加 ρ(𝔎₂+𝔎₃) ≤ ρ(𝔎₂)+ρ(𝔎₃)，eqs.(81)(82) |
| Proposition 7 | 31 | ρ(𝔎_j(P,T)) ≤ ρ(𝔊_j(2))∫ϖ²P(dω)，j∈{2,3}，eq.(83) |
| Proposition 8 | 33 | Prop 5 的标准 Eliashberg 版 |
| Lemma 3 | 33 | G_ϖ(j)=j²/(ϖ²+j²)² 递减 ⟺ ϖ ≤ √2 |
| Corollary 7 | 33 | T ≥ Ω̄(P)/2√2π ⟹ r.h.s.(86) < 0 |
| Remark 9 | 34 | 猜测 eq.(91) |
| Proposition 9 | 37 | T ↦ Λ^(2)(P,T) 可逆；λ₂ = 3/5 = 0.6，eq.(102) |

**KAY II 只有 Lemma 1、2、3——没有 Lemma 6，也没有 KAY I 那个「紧算子有限秩逼近」引理的对应物。** (A 级)

---

## §1【最重要】一般 P 下 𝔨^(N) → 𝔨 的确切陈述状态

### 1.1 裁决

> **不存在任何编号定理/推论陈述「𝔨^(N)(P,T) 单调上收敛到 𝔨(P,T)」。**
> KAY II 的 **Theorem 4 只陈述了不等式**（每个 N 一个严格不等式），**不含收敛**。
> 收敛是在**散文中断言**的，出现三次，理由只有一句「紧算子可被有限秩逼近」。

### 1.2 Theorem 4 的确切陈述（A 级，p.9）

> **Theorem 4**: *For all N ∈ ℕ, Λ(P,T) < 1/𝔨^(N)(P,T), where 𝔨^(N)(P,T) is the largest eigenvalue of 𝔎^(N)(P,T), the restriction of 𝔎(P,T) to the first N components of ℓ²(ℕ₀). The eigenvalues 𝔨^(N)(P,T) can be explicitly computed for N ∈ {1,2,3,4}. They read:* …（eqs. 10–23）

**注意**：只有 `<`，**没有 "converges"、没有 "monotonically"、没有 "N→∞"**。[A 级，逐字]

### 1.3 收敛断言出现的三处（全部 A 级逐字）

**(i) §3 引入 Theorem 4 之前（p.9，散文）**：

> "Since compact operators on separable Hilbert spaces can be arbitrarily closely approximated by their finite-rank approximations obtained by a sequence of restriction to finite-dimensional subspaces that converge monotonically to the Hilbert space, our variational principle (9) furnishes a sequence of upper approximations to Λ(P,T) that converges monotonically downward to Λ(P,T). The first four of these can be computed in closed form, using only elementary functions."

**(ii) §5.1「Proof of Theorem 4」开头（p.21）**：

> "Since 𝔎(P,T) is compact, in principle one can get arbitrarily accurate upper approximations to Λ(P,T) by restricting 𝔎(P,T) to suitably chosen finite-dimensional subspaces of ℓ²(ℕ₀). A sequence of decreasing rigorous upper bounds on Λ(P,T) that converges to Λ(P,T) is obtained by restricting the variational principle to a sequence of subspaces of ℓ²(ℕ₀) consisting of vectors Ξ_N := (ξ₀,…,ξ_{N−1},0,…), with ξ_j > 0 for j ∈ {0,…,N−1} and N ∈ ℕ. Evaluating (48) with Ξ_N in place of Ξ^opt yields a strictly monotonically decreasing sequence of upper bounds Λ^(N)(P,T) on Λ(P,T), viz." （then eq. 49）

**⚠ 该段之后的 proof 全部用于「计算 N=1,2,3,4 的闭式」，没有一行在证明收敛。** [A 级]

**(iii) §9.1 Summary（p.38）**：

> "Since a compact operator on a separable Hilbert space can be arbitrarily closely approximated by truncating it to finite-dimensional subspaces, in this case spanned by the first N positive Matsubara frequencies, we obtained from our variational principle a strictly monotonically decreasing sequence of rigorous upper bounds on Λ(P,T), the first four of which we have computed explicitly in closed form … This is the content of our Theorem 4."

（作者自己把收敛算进 "the content of our Theorem 4"，但 Theorem 4 的印刷陈述并不含它。）

**另外，p.41 有一句总结性断言（A 级）**：
> "By contrast, we proved that our sequence of lower bounds T_c^(N)(λ,P) converges upward to T_c(λ,P), for each P ∈ 𝒫."

——「we proved」，但正文中找不到对应的证明段落。**[B 级判断：这是一处 overclaim，或至少是把 KAY I Theorem 2 的 γ-模型证明默认推广到了一般 P。]**

### 1.4 收敛在何处被**当作已知**使用（A 级，p.34）

**Proof of Theorem 2**（p.34）：

> "The strict monotonicity of T ↦ Λ^(N)(P,T), given P, for T ≥ Ω̄(P)/2√2π and all N ∈ ℕ, implies via the convergence Λ^(N)(P,T) → Λ(P,T) as N → ∞, also the monotonicity of T ↦ Λ(P,T) if T ≥ Ω̄(P)/2√2π."

**⟹ Theorem 2（KAY II 的主定理之一）在逻辑上依赖这条未编号、未证明的收敛断言。** [A 级 + B 级逻辑判断]

### 1.5 单调性的来源与「严格」的地位（B 级 `[本报告自推]`）

- **单调（非严格）**是免费的：Ξ_N 子空间嵌套 ⟹ max 在更大子空间不减 ⟹ 𝔨^(N) 不减。KAY 未明写这一步。
- **严格单调**（原文 "strictly monotonically decreasing"）**未给证明**。对 𝔎 需要类似「Ξ^opt 分量全正」（Prop 8, p.33）才能保证每加一维都真严格改进。Prop 8 确实给了正性，但 KAY 没有把它接到严格单调上。
- **收敛**需要的是：紧算子 + P_N → I 强收敛 ⟹ ρ(P_N𝔎P_N) → ρ(𝔎)。这对紧自伴算子确实成立，但 KAY 只给了「紧算子可被有限秩逼近」这句话，**而 P_N𝔎P_N 是一族特定的（坐标）压缩，不是抽象的「某个有限秩逼近」——两者之间的一步（‖𝔎−P_N𝔎P_N‖→0 需要紧性 + P_N 强收敛到 I）KAY II 没写。**

### 1.6 我方结论（对 P4 的操作含义）

> **不能把「一般 P 的 𝔨^(N)↑𝔨」当成「引 KAY II Theorem X」来引用。**
> 可以说：*"KAY II asserts this (p.9, p.21, p.38) and uses it in the proof of their Theorem 2 (p.34), but it is not the content of any numbered result and no proof is given beyond the remark that compact operators admit finite-rank approximations."*
> **我方仍需自证一般-P 收敛**（好消息：这是标准命题，紧自伴 + P_N↑I ⟹ ‖𝔎−P_N𝔎P_N‖→0 ⟹ 谱收敛，可完全严格地写出，篇幅约半页）。这反而是我方相对 KAY 的一个干净的增量。

**⚠ 增量边界（补记 2026-09-10，防 overclaim）**：**不可写「N→∞ 收敛首证」。** KAY I Theorem 2 + Lemma 6（KAY I p.29）**已完整证明 γ-model 特例**的 𝔤^(N)↑𝔤。我方增量精确地只有两条：
- **(a)** 从 γ-model 的**固定核**推广到**任意谱测度 P** 的 𝔎(P,T) = ∫ℌ(ϖ)P(dω)；
- **(b)** 补 KAY 两篇**都缺的定量速率**（ρ_N = 2‖(I−P_N)ℌ‖_HS 一类）。

**附带的风险隔离**：KAY III（2409.02121）是 Einstein 极限 P(dω) = δ(ω−Ω)，即**单点测度特例**。即使 KAY III 补证了收敛，也够不到一般-P 的定位。（但 KAY III 若讨论了 Einstein 极限下的**谱隙**，会影响量化部分的定位——见 §8。）

---

## §2 紧性 / Hilbert–Schmidt 的确切陈述

### 2.1 Proposition 2（A 级，p.19，逐字）

> **Proposition 2**: *For j ∈ {1,2,3}, each 𝔎_j ∈ ℓ²(ℕ₀ × ℕ₀) for all T > 0 and P ∈ 𝒫. Thus the operators 𝔎_j = 𝔎_j(P,T) for j ∈ {1,2,3} are Hilbert–Schmidt operators that map ℓ²(ℕ₀) compactly into ℓ²(ℕ₀).*

### 2.2 证明（A 级，p.19，逐字要点）

> *Proof of Proposition 2*: "We recall that lim_{ϖ→0} (1/ϖ²)ℌ_j(ϖ) = 𝔊_j(2) for j ∈ {1,2,3}, with 𝔊(γ) = −𝔊₁(γ) + 𝔊₂(γ) + 𝔊₃(γ) the interaction matrix of the linearized γ model. **In Appendix A of [KAYa] we proved that each 𝔊_j(γ), j ∈ {1,2,3} is a Hilbert–Schmidt operator for all γ > 0.** In particular, the 𝔊_j(2), j ∈ {1,2,3} are Hilbert–Schmidt operators. We let [ℌ_j(ϖ)]_{n,m} denote the expressions bracketed [ · · · ] in (40), (41), (42). Since the maps ϖ ↦ (1/ϖ²)[ℌ_j(ϖ)]_{n,m} are non-negative and strictly monotonically decreasing for each j ∈ {1,2,3}, it follows that also the operators (1/ϖ²)ℌ_j(ϖ) are Hilbert–Schmidt operators, hence so are the operators ℌ_j(ϖ) for all j ∈ {1,2,3}."
>
> "**A moment of reflection reveals that this now implies that the ℓ² norm of the operators 𝔎_j is bounded above by ∫₀^∞ ϖ² P(dω) times the ℓ² norm of the γ model operators 𝔊_j at γ = 2.** Q.E.D."

### 2.3 显式 HS 范数上界（A 级 + B 级）

**A 级**：`‖𝔎_j(P,T)‖_HS ≤ ⟨ϖ²⟩ · ‖𝔊_j(2)‖_HS`，其中 ⟨ϖ²⟩ := ∫₀^∞ ϖ² P(dω)，ϖ = ω/2πT。（这里 "ℓ² norm of the operator" = ℓ²(ℕ₀×ℕ₀) 范数 = HS 范数。）

**⚠ 但 KAY II 没有给 ‖𝔊_j(2)‖_HS 的数值或闭式**——那要回 KAY I Appendix A（pp.42–43）。**KAY II 全文不含任何 ζ 值形式的 HS 范数上界。** 文中所有 ζ 因子（Theorem 6 eq.26、Corollary 4 的 b）都是**谱半径**上界，不是 HS 范数上界。[A 级：已全文核查]

### 2.4 关于 ℌ(ϖ) 与 𝔎(P,T) 本身（而非分量）的紧性

- ℌ_j(ϖ) 各自 HS（上引证明末段，A 级）。**没有对整体 ℌ(ϖ) = −ℌ₁+ℌ₂+ℌ₃ 单独陈述 HS**（但三 HS 之和显然 HS，B 级）。
- 𝔎(P,T) 的紧自伴性在 **Theorem 3 / eq.(9) 的陈述里直接作为形容词给出**（p.8–9 逐字）："where 𝔨(P,T) > 0 is the largest eigenvalue of a compact self-adjoint operator 𝔎(P,T) on the Hilbert space of square-summable sequences over the non-negative integers, with 𝔎(P,T) explicitly constructed in section 4.1, see (39)–(42)."
- 分量定义（A 级，p.19，eqs. 40–42）：

  (ℌ₁(ϖ)Ξ)_n = [ (2/(2n+1)) Σ_{k=1}^{n} ϖ²/(ϖ²+k²) ] ξ_n　(40)
  (ℌ₂(ϖ)Ξ)_n = Σ_{m≠n} [ 1/√(2n+1) · ϖ²/(ϖ²+(n−m)²) · 1/√(2m+1) ] ξ_m　(41)
  (ℌ₃(ϖ)Ξ)_n = Σ_m [ 1/√(2n+1) · ϖ²/(ϖ²+(n+m+1)²) · 1/√(2m+1) ] ξ_m　(42)

  且 𝔎_j = ∫₀^∞ ℌ_j(ϖ) P(dω)，𝔎 = −𝔎₁+𝔎₂+𝔎₃（p.18–19，eq.39 附近逐字："Here, ℑ is the identity operator, and 𝔎 = −𝔎₁ + 𝔎₂ + 𝔎₃, where the 𝔎_j = 𝔎_j(P,T) for j ∈ {1,2,3} are given by the P averages ∫₀^∞ ℌ_j(ϖ)P(dω) of the operators ℌ_j(ϖ) of the non-dispersive model"）。

  **⟹ 我方笔记的记号 𝔎(P,T)=∫ℌ(ϖ)P(dω)、⟦n⟧(ϖ)=ϖ²/(n²+ϖ²) 与原文一致。** ⟦n⟧ 的定义在 p.10 eq.(23) 下方逐字："with ⟦n⟧(ϖ) := ϖ²/(n²+ϖ²) ≡ ω²/((2nπT)²+ω²) for n ∈ ℕ."

---

## §3 顶本征值的单一性（simplicity / non-degeneracy）

### 3.1 KAY 如何绕开「−𝔎₁ 破坏正性」——**这是本条最有价值的发现**（A 级，p.20）

KAY **不对 𝔎（或 𝔊）本身用 Perron–Frobenius**。他们把本征值问题改写成一个**正算子的不动点问题**，把 −𝔎₁ 塞进一个**对角正预解因子**里：

> eq.(44)：(ℑ − Λ(P,T)𝔎)Ξ^opt = 0
> eq.(45)：ℭ(𝔨(P,T)) Ξ^opt = Ξ^opt
> eq.(46)：**ℭ(η) := (ηℑ + 𝔎₁)^{-1}(𝔎₂ + 𝔎₃)**

然后逐字（p.20）：

> "As in the proof of Theorem 1 in [KAYa] one shows that ℭ(η) for η > 0 is a compact operator that maps the positive cone ℓ²_{≥0}(ℕ₀) into itself, in fact mapping any non-zero element of ℓ²_{≥0}(ℕ₀) into the interior of ℓ²_{≥0}(ℕ₀), and that the spectral radius of ℭ(𝔥) equals 1. **Thus the Krein–Rutman theorem applies** and guarantees that the non-trivial solution Ξ^opt of (45) is in the positive cone ℓ²_{≥0}(ℕ₀) (after at most choosing the overall sign), hence a superconducting perturbation of the normal state Ξ̲."

**关键机制（B 级 `[本报告自推]`，但推理是初等的）**：𝔎₁ 是对角非负 ⟹ (ηℑ+𝔎₁)^{-1} 是**对角正**算子；𝔎₂+𝔎₃ 是**逐元素正**（𝔎₂ 对角为 0 但非对角正，𝔎₃ 全正，见 p.19 逐字）⟹ 乘积 ℭ(η) **逐元素正**。**所以 −𝔎₁ 的负号被吸收进预解式，正性完好无损。这正是我方问题 3 里担心的那一步，KAY 已经处理。**

### 3.1a ⚠ 但 KAY 这一步有技术瑕疵（补记 2026-09-10，等级 A 逐字 + B 判断）

原文 p.20 与 p.25 都说 ℭ(η) 把 ℓ²_{≥0}(ℕ₀) 的任意非零元 **"into the interior of ℓ²_{≥0}(ℕ₀)"**。

**问题：ℓ²_{≥0}(ℕ₀) 在 ℓ²-范数下内部为空。** 取任意全正 x ∈ ℓ²_{≥0} 与任意 ε > 0：因 x_n → 0，可取 n 使 x_n < ε/2，则 x − (ε/2)e_n 在 x 的 ε-球内却有负分量，故 x 不是内点。**⟹ 要求锥有非空内部的经典强 Krein–Rutman 不能照搬。**

这不使结论错，但**单一性不是白拿的**：
- 标准补法是换用**不需要锥内部非空**的版本（de Pagter：不可约紧正算子若 r(T) > 0 则 r(T) 是单一本征值；或 K–R 的 ideal-irreducibility 版本）。ℭ(η) 的核逐元素严格正应当足以给不可约性——**但该适用性本报告未独立验证，记为开放。**
- 同一瑕疵**同时影响 γ 侧**：p.25 那句 "It is non-degenerate" 依赖完全相同的措辞，故 §3.2(a) 的 γ-模型单一性也不能直接当既成事实引用。

**团队处置（v0.3）**：主路线已改为不依赖任何抽象锥定理——单一性由认证谱隙直接给（μ > α ∧ λ₂ ≤ α ⟹ 超过 α 的本征值恰一个 ⟹ λ₁ 严格孤立、单一）。KAY 的 ℭ(η) 改写降为可选加分项并如实标注上述瑕疵。

且（B 级 `[自算]`，等价性验算）：
ℭ(η)Ξ = Ξ ⟺ (𝔎₂+𝔎₃)Ξ = (ηℑ+𝔎₁)Ξ ⟺ (−𝔎₁+𝔎₂+𝔎₃)Ξ = ηΞ ⟺ 𝔎Ξ = ηΞ。
**⟹ ℭ[η] 的本征值-1 本征空间 与 𝔎 的本征值-η 本征空间一一对应（同一个子空间）。故 ρ(ℭ)=1 的单一性 ⟺ 𝔎 顶本征值的单一性。**

### 3.2 「non-degenerate」一词的确切出处（A 级）

**(a) γ-模型侧，p.25（Proof of Proposition 5 内），逐字**：

> "As shown in [KAYa], and analogously in the proof of Theorem 1⁺ in this paper, the operator ℭ_γ[η] maps ℓ²_{≥0}(ℕ₀) compactly into itself. Any non-zero element of ℓ²_{≥0}(ℕ₀) is mapped into the interior of ℓ²_{≥0}(ℕ₀). **Iff η = 𝔤(γ) its spectral radius equals 1. It is non-degenerate and the associated normalized eigenmode Ξ^opt(γ) ∈ ℓ²_{≥0}(ℕ₀) has only positive components, by the Krein–Rutman theorem. The analogous statement holds for all its finite-rank, N-frequency approximations Ξ_N^opt(γ), by the Perron–Frobenius theorem.**"

  其中（p.25，eq.63）：**ℭ_γ[η] := (ηℑ + 𝔊₁(γ))^{-1}(𝔊₂(γ) + 𝔊₃(γ))**，eq.(62)：ℭ_γ^(N)[𝔤(γ)]Ξ_N = Ξ_N。

**(b) 有限-N 侧，p.22（Proof of Theorem 4，N=3 情形），逐字**：

> "In [KAYa] we determined the largest eigenvalue of 𝔊^(3)(2), which equals its spectral radius, then used a continuity argument, **combined with the Perron–Frobenius theorem that establishes the non-degeneracy of the spectral radius**, to extend the determination of the largest eigenvalue of 𝔊^(3)(2) to that of 𝔊^(3)(γ) for all γ > 0."

**(c) 标准 Eliashberg 侧（𝔎 而非 𝔊）**：**p.20 只说 Krein–Rutman 给出正锥中的解，"non-degenerate" 一词未出现。** Proposition 3（p.20）只陈述「最大本征值 = 谱半径」：

> **Proposition 3**: *Let T > 0 and P ∈ 𝒫 be given. Then the largest eigenvalue 𝔨(P,T) of 𝔎(P,T) is also the spectral radius ρ(𝔎(P,T)).*

  其证明（p.20）逐字："The proof is a straightforward adaptation from the proof of Proposition 1 in [KAYa], with the same replacements needed as stated in the proof of our Theorem 1⁺."

### 3.3 裁决

| 对象 | 单一性是否被明写 | 依据 | 等级 |
|---|---|---|---|
| ρ(ℭ_γ[𝔤(γ)]) = 1（γ 模型，∞ 维） | **是**，p.25 逐字 "It is non-degenerate" | Krein–Rutman | A |
| 𝔊^(N)(γ) 顶本征值（有限 N） | **是**，p.25 + p.22 逐字 | Perron–Frobenius | A |
| **𝔤(γ) = ρ(𝔊(γ)) 本身单一** | **未直接明写**，但由 §3.1 的等价性 + p.25 立得 | — | **B** `[自推]` |
| ρ(ℭ[𝔨(P,T)]) = 1（标准 Eliashberg，∞ 维） | **未写 non-degenerate**，只写 Krein–Rutman 给正锥解 | — | A（负面） |
| **𝔨(P,T) = ρ(𝔎(P,T)) 本身单一** | **未见任何陈述** | — | **A（未见）** |
| 𝔎^(N)(P,T) 顶本征值单一 | 未见；但 Prop 8（p.33）给了 Ξ_N^opt(P,T) 正且递减 | — | C |

**⟹ 对 P4 的操作含义**：γ-模型侧（𝔊、𝔤(2)、C_∞）的单一性可以引 KAY II p.25。**标准 Eliashberg 侧（𝔎、𝔨(P,T)）的单一性 KAY II 没有明写，我方若要用，需自己补一句**（补起来很容易：ℭ(η) 强正 + 紧 ⟹ Krein–Rutman 给单一，论证与 p.25 逐字平行，KAY 自己也说 "The proof is a straightforward adaptation"）。

---

## §4 截断 = 主子块压缩的确切定义

### 4.1 KAY II 的定义（A 级，多处互相印证）

**(a) Theorem 4，p.9 逐字**：
> "…where 𝔨^(N)(P,T) is the largest eigenvalue of 𝔎^(N)(P,T), **the restriction of 𝔎(P,T) to the first N components of ℓ²(ℕ₀)**."

**(b) §5.1，p.21 逐字**（把子空间写死）：
> "…by restricting the variational principle to a sequence of subspaces of ℓ²(ℕ₀) consisting of vectors **Ξ_N := (ξ₀, …, ξ_{N−1}, 0, …)**, with ξ_j > 0 for j ∈ {0,…,N−1} and N ∈ ℕ."
>
> eq.(49)：Λ^(N)(P,T) := 1 / max_{Ξ_N} [ ⟨Ξ_N, 𝔎(P,T)Ξ_N⟩ / ⟨Ξ_N, Ξ_N⟩ ]

**(c) eq.(11) 下方，p.9 逐字**：
> "where 𝔎^(2)(P,T) is **the upper leftmost 2 × 2 block** of the matrix 𝔎^(4)(P,T)"

  同样 eq.(12) 下方："where 𝔎^(3)(P,T) is **the upper leftmost 3 × 3 block** of the matrix 𝔎^(4)(P,T)"

**(d) γ-模型侧，p.24 eq.(58) 下方逐字**：
> "𝔊^(N)(γ) = −𝔊₁^(N)(γ) + 𝔊₂^(N)(γ) + 𝔊₃^(N)(γ) is **the projection of 𝔊(γ) = −𝔊₁(γ)+𝔊₂(γ)+𝔊₃(γ) onto the subspace spanned by the first N positive Matsubara frequencies**."
>
> 且 p.25，eqs.(59)–(61) 下方："**Analogously their N-Matsubara frequency truncations are defined by limiting n and m to the set {0,…,N−1}** in (59), (60), and (61)."

**(e) Proposition 8，p.33 逐字**：
> "It should suffice to recall that **𝔎^(N) = −𝔎₁^(N) + 𝔎₂^(N) + 𝔎₃^(N)** in complete analogy to 𝔊^(N) = −𝔊₁^(N)+𝔊₂^(N)+𝔊₃^(N)"

### 4.2 裁决

> **是。截断 = 取前 N 个坐标（n,m ∈ {0,…,N−1}）的左上主子块 = P_N 𝔎 P_N，𝔨^(N) = λ_max(P_N 𝔎 P_N)。KAY II 对此有 (a)(b)(c)(d)(e) 五处互证的逐字定义。** [A 级]
>
> **KAY I Theorem 2 的「on the subset of normalized N-angle sequences Ξ_N」与之等价**：KAY II p.21 把同一记号 Ξ_N 显式写成 (ξ₀,…,ξ_{N−1},0,…)，且 p.24 eq.(58) 用 "projection onto the subspace spanned by the first N positive Matsubara frequencies" 措辞。**⟹ 等价性确认。** [A 级（KAY II 侧逐字）+ B 级（跨篇同一记号的对齐）]
>
> 注意角度变量与 ξ 变量的换算是对角正矩阵 𝔇^{1/2}（KAY II eq.(37)：(𝔇^{1/2}Θ)_n = √(2n+1) θ_n），**对角 ⟹ 不混合坐标 ⟹ 前 N 个角 ↔ 前 N 个 ξ 一一对应，主子块身份不受影响**。[B 级 `[自推]`，但是初等的]

### 4.3 一个必须写进我方论证的警告（A 级，Remark 2，p.10 逐字）

> **Remark 2**: "We have stated our Theorem 4 in a way that makes it obvious that the matrices 𝔎^(N)(P,T) are averages w.r.t. P(dω) of matrices ℌ^(N)(ϖ) of non-dispersive models, with ϖ = ω/2πT. **However, with the sole exception of the eigenvalue 𝔨^(1)(P,T) of the first approximation, all eigenvalues 𝔨^(N)(P,T) for N > 1 are not simply averages w.r.t. P(dω) of the pertinent largest eigenvalues 𝔥^(N)(ϖ) of ℌ^(N)(ϖ). The same is true for 𝔨(P,T) in relation to 𝔥(ϖ).** Each 𝔨^(N)(P,T) with N ∈ ℕ is a function of averages w.r.t. P(dω) of the expressions ω²/(ω²+4n²π²T²) with n ∈ {1,…,2N−1}."

**⟹ 算子层面 𝔎 = ∫ℌ(ϖ)P(dω) 线性，但本征值层面 𝔨 ≠ ∫𝔥(ϖ)P(dω)。任何「先取顶本征值再对 P 平均」的换序都是错的。** 这条对我方「把有限-N 展开抬到 N→∞」的论证是必须显式声明的边界。

---

## §5 算子范数收敛 / 谱隙

### 5.1 算子范数收敛 ‖𝔎 − P_N𝔎P_N‖ → 0

> **全文未见。** [A 级：42 页正文逐页核查]

最接近的只有两处**定性散文**（已在 §1.3 引全）：
- p.9："compact operators on separable Hilbert spaces can be arbitrarily closely approximated by their finite-rank approximations obtained by a sequence of restriction to finite-dimensional subspaces that converge monotonically to the Hilbert space"
- p.21："Since 𝔎(P,T) is compact, in principle one can get arbitrarily accurate upper approximations to Λ(P,T) by restricting 𝔎(P,T) to suitably chosen finite-dimensional subspaces of ℓ²(ℕ₀)."

**没有不等式、没有范数记号、没有 N 的收敛率。** KAY II 全文对收敛速度唯一的量化陈述是 Remark 1（p.8）转述 KAY I 的**数值**观察（N>200 时前 10 位稳定），以及 Theorem 5（p.11）关于 **T→0 极限**下 𝔨₀^(N) ~ ln N 发散（这是关于 T→0 而非 N→∞ 的逼近质量）。

### 5.2 谱隙 λ₁ − λ₂

> **全文未见。** [A 级]
> KAY II 从未提及第二大本征值。全文出现的 "gap" 一律指物理的 gap equation / superconducting gap。

相关但**不是**谱隙的结果，避免误引：
- **Theorem 5 / eq.(51)（p.23）**：在 **T→0** 极限，lim_{T→0} 𝔎^(N)(P,T) = −ℑ^(N) + 2Ξ_N^* ⊗ Ξ_N^*，该极限矩阵谱为 {−1（N−1 重）, −1+2Σ_{n=0}^{N−1}1/(2n+1)}。**这是极限矩阵的完整谱，看起来像"隙"，但它只在 T→0 成立，且顶本征值 ~ln N 发散，不能当作有限 T 的谱隙证据。**
- **Remark 7 + Conjecture 3（p.29，eq.74）**：Rayleigh 商在「非零、非负、递减序列」集合上的**下界猜想**，不是本征值间隙。仍是 open。
- **Remark 9（p.34，eqs.91–93）**："we suspect that…"，也是猜想。

### 5.3 一个对我方极有利的观察（B 级 `[本报告自推]`）

**Corollary 3（p.13）末句逐字**：
> "**This result also holds when N → ∞ (with superscripts ^(N) purged).**"

Corollary 3 的公式 eq.(30) 里含 **⟨𝔊^(N)(4)⟩₂**，其定义（eq.29 / eq.55）为
⟨𝔊^(N)(4)⟩₂ := ⟨Ξ_N^opt(2), 𝔊^(N)(4) Ξ_N^opt(2)⟩ / ⟨Ξ_N^opt(2), Ξ_N^opt(2)⟩，
即**用 γ=2 的顶本征向量去测 γ=4 的算子**。

**⟹ 「N→∞ 时 Corollary 3 仍成立」隐含要求 Ξ_N^opt(2) → Ξ^opt(2)（本征向量收敛），而不只是本征值收敛。KAY II 对此一字未证。而本征向量收敛正是需要谱隙的地方（顶本征值单一 + 隙 > 0 ⟹ 谱投影收敛）。**

**这是 P4 「补谱隙」工作的最强的一手动机**：KAY 自己的 Corollary 3 在 N→∞ 处有一个未填的洞，而填它恰好就是我方要做的谱隙估计。（且 P4-C 已数值发现 𝔊(2) 的隙极大，λ₂ ≈ 5e-4 ≪ λ₁ = 1.318，见 `surveys/P4_KAY_I_and_combescot_2026-09-09.md` §7.2——**注意那是数值，尚非严格界**。）

---

## §6 U3 悬案裁决（终结）

**旧记录**（`surveys/P4_KAY_I_and_combescot_2026-09-09.md` §5.3）：
> "KAY I 中不存在 Corollary 3；N→∞ 收敛在 KAY I 是 Theorem 2；一般 P(dω) 版本应在 KAY II（未核）。"

**本轮裁决**：

1. **「KAY I 无 Corollary 3」——正确，维持。** [沿用 A 级]
2. **「KAY II 有 Corollary 3」——正确，但它与 N→∞ 换序无关。** KAY II Corollary 3（p.13）是 **T_c^(N)(λ,P) 的大-λ 渐近闭式 eq.(30)**，是 Theorem 7（T=∞ 附近解析性 + 一阶微扰论）的推论。[A 级]
3. **「一般 P(dω) 版本应在 KAY II」——部分错。** KAY II **确实断言**了一般 P 下的 Λ^(N)↓Λ（等价 𝔨^(N)↑𝔨），并在 Theorem 2 的证明里使用它；**但它不是 KAY II 任何编号结果的内容，也没有证明**（只有一句紧算子有限秩逼近的散文理由）。[A 级]
4. **我方笔记里那个叫「KAY II Corollary 3」的收敛结果——不存在。这是一次编号误挂。** 正确的引用应是：**KAY II, p.9（§3 引言）/ p.21（Proof of Theorem 4 首段）/ p.38（§9.1 Summary）的散文断言，以及 p.41 的 "we proved that our sequence of lower bounds T_c^(N)(λ,P) converges upward to T_c(λ,P)" 这句无对应证明的总结。** [A 级]

**⟹ U3 关闭。结论：一般-P 收敛「KAY 说了但没证」，我方必须自证。**

---

## §7 五条结论速查

| # | 问题 | 裁决 | 等级 | 出处 |
|---|---|---|---|---|
| 1 | 一般 P 下 𝔨^(N)↑𝔨 是否有编号定理？ | **否。Theorem 4 只有不等式。收敛仅散文断言（p.9/21/38），且被 Theorem 2 的证明依赖（p.34）。无 Lemma 6 类支撑。我方须自证。** | A | p.9, 21, 34, 38, 41 |
| 2 | 紧性 / HS 陈述与显式界？ | **有：Proposition 2（p.19）：𝔎_j ∈ ℓ²(ℕ₀×ℕ₀)，HS 且紧。显式界 ‖𝔎_j‖_HS ≤ ⟨ϖ²⟩·‖𝔊_j(2)‖_HS。但 ‖𝔊_j(2)‖_HS 的数值/ζ 表达式不在 KAY II（回 KAY I App. A）。ℌ_j(ϖ) 也是 HS。** | A | p.19, eqs.40–42 |
| 3 | 顶本征值单一性？如何绕过 −𝔎₁？ | **绕法：不对 𝔊/𝔎 用 PF，而对 ℭ(η)=(ηℑ+𝔎₁)^{-1}(𝔎₂+𝔎₃) 用 Krein–Rutman——负号被吸进对角正预解式（eq.46, p.20）。γ 模型侧明写「non-degenerate」（p.25）；标准 Eliashberg 侧从未明写。⚠ 但 KAY 这一步的 "into the interior of ℓ²_{≥0}" 措辞不成立（该锥在 ℓ²-范数下内部为空），经典强 K–R 不能照搬，γ 侧那句也因此不白拿；补法需 de Pagter / ideal-irreducible 版本，本报告未验证。见 §3.1a。** | A（逐字）/ B（瑕疵判断） | p.20, 22, 25 |
| 4 | 截断 = 左上 N×N 主子块？ | **是，五处逐字互证（Thm 4 "restriction to the first N components"；p.21 Ξ_N=(ξ₀,…,ξ_{N−1},0,…)；p.9 "upper leftmost 2×2 / 3×3 block"；p.24 eq.58 "projection onto the subspace spanned by the first N positive Matsubara frequencies"；p.25 "limiting n and m to {0,…,N−1}"）。与 KAY I Thm 2 的「N-angle sequences」等价（𝔇^{1/2} 对角）。⚠ 但 Remark 2（p.10）警告：本征值层面不可与 P-平均换序。** | A | p.9, 10, 21, 24, 25, 33 |
| 5 | 算子范数收敛？谱隙？ | **两者全文均未见。** 无 ‖𝔎−P_N𝔎P_N‖→0 的任何不等式，无收敛率，无第二本征值/谱隙的任何提及。**这正是我方要补的。** 加分观察：Corollary 3 末句「N→∞ 亦成立」隐含需要本征向量收敛（⟨𝔊^(N)(4)⟩₂ 用 Ξ_N^opt(2) 定义），而那正需谱隙——KAY 未证，是我方工作的直接动机。 | A（未见）+ B（动机推断） | 全文；p.13, 29, 34 |

---

## §8 遗留不确定点

- **U1（KAY II 期刊卷号）**：本轮读的是 arXiv v3 PDF，**正文/首页无期刊卷号与 DOI**（仅 arXiv DOI 10.48550/arXiv.2409.00532）。旧记录的 "J. Stat. Phys. vol.92, art.94" 来自 arXiv comments 字段，**本轮未能从 PDF 内证实**。引用前仍须独立核 DOI。[C 级，维持存疑]
- **p.41 的 "we proved that our sequence of lower bounds T_c^(N)(λ,P) converges upward to T_c(λ,P), for each P ∈ 𝒫"**：我在正文中找不到对应证明。可能作者视其为 KAY I Theorem 2 的直接移植（"the same replacements" 套路），但一般 P 的移植需要 §1.5 那一步。**建议我方在论文中以中性措辞描述此状态，不要断言 KAY 有错。** [B 级判断]
- **KAY I Appendix A 的 ‖𝔊_j(γ)‖_HS 显式上界数值**：本轮未取（要回 2409.00533v2 pp.42–43）。若我方要把 §2.3 的 ⟨ϖ²⟩·‖𝔊_j(2)‖_HS 变成可算的数，这是下一步。
- **KAY III（2409.02121，Einstein 声子）本轮完全未读**。若 KAY III 补证了收敛或给了谱隙，本报告的 §1/§5 结论需修订。**建议下一轮优先核 KAY III。**
