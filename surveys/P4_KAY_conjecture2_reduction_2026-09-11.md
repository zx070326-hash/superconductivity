# P4 · KAY III Conjecture 1 逐字核对 + 「一般 P ⇐ Einstein」化约是否已被 KAY 写出

检索 agent：researcher　日期：2026-09-11　派单：`conj2-kay-check`（见 `theory/derivations/P4_conjecture2.md` §4）

**一手材料（本轮全部 PDF 页图逐页读，非 HTML 摘要级）**
- KAY III = arXiv:**2409.02121v2**（页 1 斜纹戳 "arXiv:2409.02121v2 [cond-mat.supr-con] 7 Jun 2025"，扉页 "Revised version of April 12, 2025"）：**pp.1–30 全读**（正文至 p.29，App.A 起 p.30；App.A 其余 pp.31–34 与参考文献 pp.35–36 本轮未重读，2026-09-10 已读过）。
- KAY II = arXiv:**2409.00532v3**：**pp.3–42 全读**（pp.1–2 为扉页/摘要，未重读）。
- 竞争文献：arXiv:**2512.20009v2**（Elezaby–Abanov，"Dated: June 11, 2026"，PRB 10.1103/wpn8-8j2t）**pp.1–12 页图读**。

**证据等级**：**A** = 原文逐字（页码 = PDF 印刷页码）；**B** = 我方对 A 级材料的确定性推算/判断（标 `[自算]`）；**C** = 未取到/推断。

**速览**
| # | 问题 | 裁决 | 等级 |
|---|---|---|---|
| 1 | Conjecture 1 是否就是逐点 𝔥(ϖ)<𝔤(2)ϖ²？ | **不是逐字如此**：原文是 **T_c 形式、严格 "<"、且把「T_c 存在」也包进猜想**。Tc 无 sup/首次穿越定义，只有「临界面是图」时的反函数定义。**界的那一半与逐点不等式等价（B 级推导）；存在性那一半 = 单调性猜想，我方不证。** | A + B |
| 2 | KAY 是否写出「一般 P ⇐ Einstein」化约 / 特征值 Jensen？ | **未找到**（KAY II pp.3–42、KAY III pp.1–30 逐页）。最接近的四处见 §2，都不是。**KAY II Prop 7 (eq.83) 是「正部分」的同型步骤**，值得引用以示诚实。 | A（未见） |
| 3 | 单调性证到哪？ | 𝔥^(N) 对**所有 N** 在 **[0,ϖ_*)、ϖ_*≥√2** 严格增（Thm 5）；N=∞ 经「收敛」散文转移；全 ℝ₊ 明写 "have not succeeded"；N=3,4 "visual inspection"；**𝔥(ϖ)/ϖ² 的单调性全文零处**。 | A |
| 4 | 级数 Σ(−ϖ²)^j𝔊(2j+2) 是否出现？ | **只到二阶 "± ⋯"**（KAY II eq.53/54、KAY III eq.57）；**无一般项、无算子范数收敛半径、无 𝔊(6)**；‖𝔊(γ)‖ 的界只有 KAY I Prop 5/6 与 App.A 的 HS 界。 | A |
| 5 | Conjecture 1 的数值证据？ | 只有 **Fig.1（λ∈[0,5]）与 Fig.2（1/√λ∈[0,2]）里 N≤4 下界曲线落在虚线之下**；**无 𝔥(ϖ)/ϖ² 图、无 ϖ 区间、无 N=200 曲线**。小 ϖ：Thm 7 给二阶系数 ⟨𝔊^(N)(4)⟩₂>0（有限 N 定理），N=∞ 只有一句 "also holds"。 | A |
| 6 | 2024–2026 有人证了吗？ | **无**。唯一新的严格上界文献 2512.20009 只做 γ-model（Gershgorin，𝔤(2)≤≈1.545 `[自算]`），不碰 Einstein/Conjecture 1。 | A/B |

---

## §1 Conjecture 1 的精确陈述、严格性、T_c 定义

### 1.1 逐字（KAY III p.12，eq.28）[A]

> **Conjecture 1**: *There is a critical temperature T_c(λ,Ω) > 0 which for all Ω > 0 and λ > 0 is bounded above by T_c(λ,Ω) < T_c^~(λ,Ω), with*
> T_c^~(λ,Ω) := (Ω/2π)√(𝔤(2)λ);　(28)
> *here, 𝔤(2) is the spectral radius of the operator 𝔊(γ) at γ = 2.*

其前置句（p.12，紧接 Corollary 3 之后）[A]：
> "For large λ the r.h.s.(27) can be expanded to yield T_c^(N)(λ,Ω) ∼ (Ω/2π)√(𝔤^(N)(2)λ), with N ∈ ℕ. By a simple convexity estimate, r.h.s.(27) ≤ (Ω/2π)√(𝔤^(N)(2)λ). And so, for large enough λ the asymptotic expression (Ω/2π)√(𝔤^(N)(2)λ) is an upper bound on T_c^(N)(λ,Ω) that is asymptotically sharp as λ ∼ ∞. Moreover, in [KAYa] we showed that 𝔤^(N)(2) converges upward to 𝔤(2). Furthermore, as noted earlier, each T_c^(N)(λ,Ω) vanishes for λ ≤ λ_N, while √λ > 0 for all λ. All the above, plus the rapid convergence for λ > 0.77 of the T_c^(N)(λ,Ω) with N discernible in Fig. 1 now suggests"

KAY II 的对应物 Conjecture 2（p.14，eq.31）[A]（与 2026-09-09 survey §82 一致，本轮页图复核）：
> **Conjecture 2**: *There is a critical temperature T_c(λ,P) > 0 which for all P ∈ 𝒫 and λ > 0 is bounded above by T_c(λ,P) < T_c^~(λ,P), with* T_c^~(λ,P) := (1/2π)√(𝔤(2)⟨ω²⟩λ);　(31)

**三个逐字事实**：
1. **严格不等号 "<"**（两篇同）。
2. **不是逐点算子陈述**：全文没有 "𝔥(ϖ) ≤ 𝔤(2)ϖ²" 或 "L_E(ϖ) ≥ 1/(𝔤(2)ϖ²)" 这样的句子。
3. **猜想包含存在性**："There is a critical temperature T_c > 0 which for all … λ > 0 …"——即对**所有** λ>0 都存在 T_c。而 KAY III 已证的存在性只到 λ ≥ λ_*、λ_* < 0.7670（Thm 5）。

### 1.2 KAY 的 T_c 定义 [A]

- **临界面**（Thm 1，p.5）："The boundary between the two regions, called the critical surface 𝒮_c, is a graph over the positive (Ω,T)-quadrant, i.e. 𝒮_c = {(λ,Ω,T) ∈ ℝ³₊ : λ = Λ_E(Ω,T)}. … Λ_E(Ω,T) = L_E(ϖ)."
- **T_c 只在临界面对 (λ,Ω) 也是图时定义**（Thm 5，p.8，eq.17）："Thus the portion of the critical surface 𝒮_c over the region {λ ≥ λ_*} in the positive (λ,Ω)-quadrant is also a graph, yielding the critical temperature T_c(λ,Ω), viz. 𝒮_c|_{λ≥λ_*} = {(λ,Ω,T) ∈ ℝ³₊ : T = T_c(λ,Ω), λ ≥ λ_*}. Moreover, T_c(λ,Ω) = Ωf(λ), where f(λ) is continuous and strictly monotonically increasing for λ ≥ λ_*. Furthermore, λ_* < 0.7670…"
- **多次穿越的可能性被明写**（pp.7–8）："Our Theorems 1 – 4 do not rule out that some lines ℒ(λ,Ω) := {(λ,Ω,T) ∈ ℝ³₊ : λ = c₁ & Ω = c₂} could pierce 𝒮_c more than once, in which case the critical surface would not be a graph over the positive (λ,Ω)-quadrant … for a multiple piercing would mean that there is no unique critical temperature T_c(λ,Ω) for certain (λ,Ω)."
- **全文无 "sup"、无 "first crossing"、无 "largest solution" 式的备用定义。**（KAY II 同：Thm 2 p.8 eq.8 只在 λ ≥ λ_*(P) 定义。）
- 但**临界曲线对每个 λ>0 至少有一个点**（KAY III p.11）："the continuous critical curve 𝒞_c that divides the positive (λ,T/Ω)-quadrant into simply connected normal and superconducting regions, and which is a graph over the T/Ω axis, starts at (0,0)."

### 1.3 「界的那一半」与逐点不等式的等价 `[自算，B]`

设 Lemma B：𝔥(ϖ) < 𝔤(2)ϖ² 对所有 ϖ>0。由 Thm 2（L_E = 1/𝔥）：
- Lemma B ⟺ L_E(ϖ) > 1/(𝔤(2)ϖ²) ∀ϖ>0 ⟺ [(λ,Ω,T)∈𝒮_c ⟹ λ = L_E(ϖ) > 1/(𝔤(2)ϖ²) ⟹ (2πT)² < 𝔤(2)λΩ² ⟹ T < T_c^~(λ,Ω)]。
- 反向：若「𝒮_c 上每一点都满足 T < T_c^~」，取任意 ϖ>0、令 λ := L_E(ϖ)>0，则 (λ,Ω,Ω/2πϖ)∈𝒮_c，得 ϖ² > 1/(𝔤(2)λ) = 𝔥(ϖ)/𝔤(2)，即 Lemma B。
- 更强的表述（不依赖 T_c 定义）：Lemma B ⟹ 对 T ≥ T_c^~ 有 λ𝔥(ϖ) ≤ λ𝔤(2)ϖ² ≤ 1（T = T_c^~ 时 λ𝔥 < 1 因严格），故**正常态在 [T_c^~, ∞) 上全部线性稳定**；线 ℒ(λ,Ω) 与 𝒮_c 的**所有**交点（包括可能的多次穿越）都在 T < T_c^~。
- **Lemma B 不给「唯一穿越」**，故不证 Conjecture 1 的 "There is a critical temperature … for all λ>0" 那一半——那一半 = KAY III p.28 的单调性猜想（§3.3 引）。

**结论**：Conjecture 1 = (存在性猜想) ∧ (界)。**Lemma B ⟺ (界)，且给出比 (界) 更强的「整个不稳定集 ⊂ (0,T_c^~)」。**

---

## §2 「一般 P ⇐ Einstein」化约：KAY 未写出 [A，未见]

**逐页核对范围**：KAY II pp.3–42；KAY III pp.1–30。检索目标：任何形式的 (i) λ₁ 的凸性/Jensen 用于 𝔎 = ∫ℌ(ϖ)P(dω)；(ii) "it suffices to prove the Einstein case" / "reduces to the non-dispersive case" / "Conjecture 1 implies Conjecture 2"；(iii) 𝔨(P,T) ≤ ∫𝔥(ϖ)P(dω)。**三类均未出现。**

**最接近的五处（逐字，均不构成化约）**：

**(a) KAY II Remark 2（p.10）——只说「不是平均」，不说方向** [A]
> "We have stated our Theorem 4 in a way that makes it obvious that the matrices 𝔎^(N)(P,T) are averages w.r.t. P(dω) of matrices ℌ^(N)(ϖ) of non-dispersive models, with ϖ = ω/2πT. However, with the sole exception of the eigenvalue 𝔨^(1)(P,T) of the first approximation, all eigenvalues 𝔨^(N)(P,T) for N > 1 are not simply averages w.r.t. P(dω) of the pertinent largest eigenvalues 𝔥^(N)(ϖ) of ℌ^(N)(ϖ). The same is true for 𝔨(P,T) in relation to 𝔥(ϖ)."

⟹ KAY 明确意识到 𝔨 ≠ ∫𝔥P，但**没有写不等式 𝔨 ≤ ∫𝔥P**（这正是我方引理 A）。

**(b) KAY II Corollary 2（p.12）——用了 Jensen，但是标量凹函数、方向相反** [A]
> "**Corollary 2**: By Jensen's inequality, P(dω)-averaging the concave function ω² ↦ ω²/(ω²+4π²T²) yields the weaker upper bound 𝔨(P,T) ≤ 𝔨^♯(P,T), where 𝔨^♯(P,T) = ⟨ϖ²⟩/(⟨ϖ²⟩+1) + 2√((2^{1+ε}−1)ζ(1+ε)ζ(5−ε))⟨ϖ²⟩ (27)"

⟹ 这里的 Jensen 作用在 𝔨^(1) 这个**标量** P-平均上，把 Thm 6 的界放松成只含 ⟨ϖ²⟩ 的界。**不是特征值泛函 λ₁(·) 的凸性**。

**(c) KAY II p.13 / KAY III p.12 的 "simple convexity estimate"** [A]
> KAY II p.13："By a simple convexity estimate, r.h.s.(30) ≤ (1/2π)√(𝔤^(N)(2)⟨ω²⟩λ)."

⟹ 这是大-λ 渐近公式 (30)/(27) 里对 √ 的凸性估计，与 P 无关。

**(d) KAY III §3.3（p.23）——「Einstein 严格下界对一般 P 一致成立」，但走的是谱半径次可加路线** [A]
> "We add the remark that our lower bound L_E^*(ϖ) on L_E(ϖ) is a uniform lower bound on the analogous function in standard Eliashberg theory with dispersive phonons, i.e. on Λ(P,T), for all P(dω) ∈ 𝒫 supported on [0,Ω]. This follows by inspection of the proof of Theorem 6 in [KAYb]."

⟹ 这是**已证界**（Thm 4 的 𝔥*，含 2.034 因子）的 Einstein→一般 P 转移，机制是 KAY II §6 的 Lemma 1/2 + Prop 7，**不是**λ₁ 凸性，**不涉及** Conjecture 1/2。

**(e) KAY II Proposition 7（p.31，eq.83）——「正部分」的 P-平均 ≤ Einstein 系数 × ⟨ϖ²⟩** [A]
> "**Proposition 7**: Let T > 0 and P ∈ 𝒫 be given. Then, with ϖ = ω/2πT, for j ∈ {2,3} we have ρ(𝔎_j(P,T)) ≤ ρ(𝔊_j(2)) ∫₀^∞ ϖ² P(dω). (83)
> Proof: All kernels of the operators 𝔎_j(P,T), j ∈ {2,3}, are of the type ∫₀^∞ ϖ²/(ϖ²+k²) P(dω) > 0 for k ∈ ℕ, while all kernels of 𝔊_j(γ), j ∈ {2,3}, are of the type 1/k^γ > 0 for k ∈ ℕ. So the spectral radii of all these operators are achieved on the positive cone ℓ²_{≥0}(ℕ₀). Since ϖ² + k² > k², the claim of the proposition follows."

⟹ **这是 KAY 最接近我方化约的一步**：对**逐元素为正**的 𝔎₂、𝔎₃，用 Perron 型论证得到「P-平均算子的谱半径 ≤ 𝔊_j(2) 的谱半径 × ⟨ϖ²⟩」。它对含 −𝔎₁ 的完整 𝔎 **不适用**（KAY 因此走了次可加，丢掉 −𝔎₁，这正是 2.034 缺口的来源，见 2026-09-09 survey）。我方引理 A（λ₁ 是自伴算子上的凸泛函 ⟹ λ₁(∫ℌP) ≤ ∫λ₁(ℌ)P）**绕过了正性要求**，对完整 𝔎 成立。**论文里应引 Prop 7 并说明区别，不应写「KAY 未考虑过 P-平均与 Einstein 系数的关系」。**

**(f) 方向相反的陈述（KAY III p.4–5、p.14、p.23、p.27）** [A]：全部是「一般 P 的结果在 P→δ 极限下给出 Einstein 结果」，例如 p.23 "its proof is included as limiting case P(dω) → δ(ω−Ω)dω with Ω > 0 of the proof of Theorem 6 in [KAYb]"。**没有一处反向。** 且 KAY III p.12 把 Conjecture 1 当作独立猜想陈述，**未说** "Conjecture 1 would imply Conjecture 2 of [KAYb]"。

**裁决**：**「一般 P ⇐ Einstein」化约（引理 A + 引理 B ⟹ Conjecture 2 之界）在 KAY I–III 中不存在，是我方新写的**——但它是一行 Jensen，学术分量小；诚实的写法是「化约是初等的（λ₁ 的凸性），但据我们所知未被 [KAY II, III] 指出；其正部分版本见 [KAY II, Prop. 7]」。

---

## §3 单调性：已证范围与 "visual inspection" 的位置 [A]

### 3.1 定理级（KAY III Theorem 5，p.8）
> "**Theorem 5**: For all N ∈ ℕ the eigenvalues 𝔥^(N)(ϖ) increase strictly monotonically with ϖ ∈ [0,ϖ_*). Moreover, ϖ_* ≥ √2. As a consequence, the map ϖ ↦ L_E(ϖ) is strictly monotone decreasing on [0,ϖ_*], with ϖ_* ≥ √2. Thus the portion of the critical surface 𝒮_c over the region {λ ≥ λ_*} in the positive (λ,Ω)-quadrant is also a graph, yielding the critical temperature T_c(λ,Ω), viz. … Moreover, T_c(λ,Ω) = Ωf(λ), where f(λ) is continuous and strictly monotonically increasing for λ ≥ λ_*. Furthermore, λ_* < 0.7670…."

- **范围**：ϖ ∈ [0,√2]（至少）；**N**：定理陈述对**所有有限 N**；**N=∞** 的 "As a consequence" 一步靠 KAY II p.34 Proof of Theorem 2 的散文："implies via the convergence Λ^(N)(P,T) → Λ(P,T) as N → ∞, also the monotonicity of T ↦ Λ(P,T) if T ≥ Ω̄(P)/2√2π. Since T ↦ Λ(P,T) is analytic for each P ∈ 𝒫, it cannot have a constant piece, for this would violate its lower bound Λ*(P,T); hence T ↦ Λ(P,T) is strictly monotonic increasing, too."（严格性靠解析性 + 显式下界排除平台。）
- **机制**（KAY II p.33–34）：Prop 8（Ξ^opt_N/√(2n+1) 正且递减）+ **Lemma 3**："Let ϖ > 0. Then the map j ↦ G_ϖ(j) := j²/(ϖ²+j²)², j ∈ ℕ, is decreasing iff ϖ ≤ √2." + Hardy–Littlewood–Pólya 重排 ⟹ eq.(90) ⟨Ξ^opt_N,[d/dϖ² ℌ^(N)(ϖ)]Ξ^opt_N⟩ > 0。**√2 是核 j²/(ϖ²+j²)² 在 j∈ℕ 上单调的临界值，是证明技术的界限**（KAY II p.39："due to the limitations of our technique of proof and presumably not intrinsic to the model"）。
- **KAY III 的证明只是一句转引**（p.24）："the monotonicity of ϖ ↦ L_E(ϖ) for ϖ ∈ [0,√2] follows simply by the specification P(dω) = δ(ω−Ω)dω in the proof of monotonicity of T ↦ Λ(P,T) for T > T_*(P), given P, in [KAYb]."

### 3.2 全 ℝ₊ 未证 + "visual inspection"（KAY III p.8）
> "While we have not succeeded in showing that the map ϖ ↦ 𝔥(ϖ) is strictly monotonic increasing for all ϖ ∈ ℝ₊, our lower bounds to 𝔥(ϖ) stated explicitly in Theorem 3 for N ∈ {1,2,3,4} all are strictly monotonic increasing with ϖ. This is manifestly obvious only for 𝔥^(1)(ϖ). For 𝔥^(2)(ϖ) this is a consequence of Proposition 9 in [KAYb]. For 𝔥^(3)(ϖ) and 𝔥^(4)(ϖ) the monotonicity for λ < λ_* is vindicated through visual inspection of the plots (see below)."

p.26 §3.5.3（N=3）："With some extra (not too hard, but daunting) work, one should be able to rigorously prove that λ ↦ T_c^(3)(λ,Ω) is strictly monotonic increasing for all λ ≥ λ_3, but here we are content with pointing out that the plot of our parameter representation for T_c^(3)(λ,Ω) in Fig. 1 reveals that there is no sudden 'horizontal oscillation' in the critical curve 𝒞_c^(3) for λ ∈ (0.4838…, 1.0158…)." §3.5.4（N=4）同样措辞，区间 (0.4251…, 0.7670…)。

p.28 Summary："We also supplied compelling evidence for the conjecture that the map ϖ ↦ L_E(ϖ) is strictly monotone decreasing for all ϖ ∈ [0,∞), but to prove it would require a different strategy."

p.9 另一个未编号猜想："We conjecture that all 𝒞_c^~(N) are such graphs for general N ∈ ℕ."

### 3.3 ϖ→∞ 端（Theorem 6，p.11）
𝔥^(N)(ϖ) = 𝔥^(N)(∞) − B_N/ϖ² + O(1/ϖ⁴)，B_N ≥ 1 ⟹ "It establishes that the L_E^(N)(ϖ) are asymptotic to strictly monotonically decreasing functions when ϖ ∼ ∞." **只对有限 N**（𝔥^(N)(∞) ~ ln N 发散）；对 N=∞ 无大-ϖ 单调性陈述。

### 3.4 𝔥(ϖ)/ϖ² 的单调性：**全文零处** [A，未见]
KAY I–III 从未考虑比值 𝔥(ϖ)/ϖ²（或 ϖ²L_E(ϖ)）作为 ϖ 的函数。仅有两端信息：ϖ→0 时 → 𝔤(2)（Thm 7），有限 N 时 ϖ→∞ 趋 0。

⟹ **我方 §3「Hellmann–Feynman + 认证谱隙」把 N=∞ 的 𝔥 单调性从 [0,√2] 推到 [0.2,60] 类区间是净新增；比值 r(ϖ) 的严格单调下降更是净新增。** 但注意 KAY 已证的是 [0,√2] 含 0 端，我方证书网格不含 ϖ→0 端，两者互补而非包含。

---

## §4 级数展开 Σ(−ϖ²)^j 𝔊(2j+2) 与 ‖𝔊(γ)‖ 的界 [A]

### 4.1 KAY 写出的展开只到二阶
KAY II p.23（Proof of Theorem 7）：
> "Since all kernels of ℌ_j(ϖ), j ∈ {1,2,3}, are of the type ϖ²/(ϖ²+k²) > 0, with k ∈ ℕ, the operators ℌ_j(ϖ), j ∈ {1,2,3}, and their N-frequency truncations, are real analytic in ϖ, and in ϖ², about ϖ = 0 = ϖ². Maclaurin expansion in powers of ϖ² about ϖ² = 0 yields
> ϖ²/(ϖ²+k²) = ϖ²·(1/k²) − ϖ⁴·(1/k⁴) ± ⋯, (52)
> which implies that
> ℌ_j(ϖ) = ϖ²𝔊_j(2) − ϖ⁴𝔊_j(4) ± ⋯, j ∈ {1,2,3}, (53)
> and the analog holds for their N-frequency truncations. Averaging over P then yields the asymptotic large T expansion
> 𝔎_j(P,T) = 𝔊_j(2)⟨ω²⟩/(4π²T²) − 𝔊_j(4)⟨ω⁴⟩/(16π⁴T⁴) ± ⋯, j ∈ {1,2,3}, (54)
> … The claim (28) of the theorem now follows from first-order perturbation theory [K]."

KAY III p.23 eq.(57) 逐字同 (53)。另 KAY II p.22："We recall that lim_{ϖ→0} (1/ϖ²)ℌ^(N)(ϖ) = 𝔊^(N)(2)"；KAY II p.19 Prop 2 证明："lim_{ϖ→0} (1/ϖ²)ℌ_j(ϖ) = 𝔊_j(2)"。

**逐字事实**：
- **一般项 (−ϖ²)^j 𝔊_j(2j+2) 未写出**，只有 "± ⋯"；**𝔊(6) 及更高从未出现**。
- **"real analytic in ϖ²" 是断言**，无收敛半径、无算子范数意义下的收敛证明（逐元素显然 |ϖ|<1；作为 ℓ² 算子级数需 ‖𝔊(2j+2)‖ 一致有界，KAY 未讨论）。
- **一阶微扰论 [K]（Kato）用于 N 有限**；N=∞ 的 (28)/(26) 展开在 Corollary 3 后一句 "This result also holds when N → ∞ (with the superscripts ^(N) purged)" 断言，无证明（与 2026-09-10 报告 U3 一致）。

### 4.2 ‖𝔊(γ)‖（γ>2）的界：**KAY II/III 无新界** [A，未见]
KAY II 全文的 ζ 因子只在 Thm 6 eq.(26)、Cor 2 eq.(27)、Cor 4 的 b 中，全是 γ=2 谱半径界的转写；KAY III 的 eq.(16)/(19) 同。**γ>2 的 ‖𝔊_j(γ)‖ 界仍只有 KAY I Prop 5/6（谱半径，对所有 γ 成立）与 KAY I App.A Prop A/B/C（HS，见 2026-09-10 报告 §3）。** 我方段 I 需要的 Σ_{j≥2}ϖ^{2j}‖𝔊(2j+2)‖ 一致界须自证（B 级），可用的一手成分：‖𝔊₃(γ)‖_HS ≤ ζ(1+γ)（KAY I eq.123）、‖𝔊₂(γ)‖²_HS ≤ (π²/4)ζ(2γ)（成分 A，组装 B）、𝔊₁ 无一手界。

### 4.3 ⟨𝔊(4)⟩₂ > 0 的一手来源
KAY II Proposition 4（p.24，eq.56）= KAY III Proposition 5（p.23，eq.59）：
> "Let γ > 0 be given. Then for all γ′ > 0 and N ∈ ℕ₀, ⟨𝔊^(N)(γ′)⟩_γ := ⟨Ξ^opt_N(γ), 𝔊^(N)(γ′)Ξ^opt_N(γ)⟩ / ⟨Ξ^opt_N(γ), Ξ^opt_N(γ)⟩ > 0, with Ξ^opt_N(γ) any eigenvector of the top eigenvalue 𝔤^(N)(γ) of 𝔊^(N)(γ)."

证明（KAY II pp.24–29）：Prop 5 重排不等式 + eq.(72)–(73) Dirichlet 级数系数非负。**只对有限 N**；KAY II p.42 致谢 "the referee who spotted a mistake in our original proof of Proposition 4, now corrected"。**数值 ⟨𝔊(4)⟩₂ 在 KAY II/III 中未给**（我方 1.19317 是新数；KAY I 本轮未复查）。Remark 7 + Conjecture 3（p.29，eq.74）是 l.h.s.(70) 的下界猜想，与 ⟨𝔊(4)⟩₂ 的定量下界相关但**仍 open**。

---

## §5 Conjecture 1 的数值证据 [A]

**KAY III p.27 §3.6.2 逐字**：
> "The discussion in section 2 establishes that T_c^~(λ,Ω) in Conjecture 1 is an upper bound on T_c(λ,Ω) for large enough λ. This is a consequence of Corollary 3 to Theorem 7, which establishes Proposition 1.
> Recall that Conjecture 1 proposes that T_c^~(λ,Ω) is an upper bound on T_c(λ,Ω) for all λ > 0 and Ω > 0. Fig. 1 and Fig. 2 present numerical evidence for its veracity."

**证据内容**：
- **Fig.1（p.10）**：横轴 λ∈[0,5]，纵轴 T/Ω ∈ [0,0.5]；曲线 T_c^(N)/Ω（N=1..4）、T_c^*/Ω（红，已证上界）、T_c^~/Ω（黑虚线，猜想）。标题："All displayed curves are plots of λ as given by some explicitly computed elementary functions of T/Ω."
- **Fig.2（p.13）**：横轴 1/√λ ∈ [0,2]（即 λ ≥ 0.25），纵轴 T/(Ω√λ) ∈ [0,0.20]；N=1..4 下界 + 虚线 T_c^~。标题："The third and fourth approximates virtually agree with each other if 1/√λ < 0.5, and only minute discrepancies are visible when 0.5 < 1/√λ < 1, indicating the rapid convergence of our lower bounds to the exact result when λ > 1. All four lower bounds hit zero, the N-th one at 1/√λ_N ≥ 1. The conjectured global upper bound T_c^~(λ,Ω) (dashed line) is visibly seen to be asymptotically exact (as proved), yet not very accurate away from the asymptotic regime λ ∼ ∞. Our rigorous bound (not shown) would be an almost horizontal line twice as high as the dashed one."
- p.10："Further visual inspection of Figure 1 reveals that the sequence of lower bounds T_c^(N)(λ,Ω) on T_c(λ,Ω) appears to converge upward very rapidly to some limiting curve T_c(λ,Ω) = Ωf(λ) when λ > 0.7 (say)."
- p.12："it is a challenging task to improve it to the precision that is suggested by a small ϖ analysis of the operators ℌ^(N)(ϖ), and the apparent rapid convergence of the sequence of eigenvalues 𝔥^(N)(ϖ) for small ϖ."

**逐字事实**：
- **证据只是「N≤4 的下界曲线位于虚线之下」**。T_c^(N) ↑ T_c，所以这只说明 T_c^(4) < T_c^~，对 T_c 本身的支持依赖「N=4 已近收敛」的视觉判断。**无 N=200 曲线**（N=200 只用于 𝔤(2) 的 10 位数，Prop 1 p.14）。
- **无 𝔥(ϖ)/ϖ² 对 ϖ 的图**，无任何以 ϖ 为横轴的图；**无显式 ϖ 区间**。折算：Fig.2 的 λ≥0.25 与 T_c^(4) 曲线对应的 ϖ 范围未标。
- **小 ϖ 行为**：Theorem 7（p.12，eq.26）"The eigenvalues 𝔥^(N)(ϖ) are analytic about ϖ = 0, with 𝔥^(N)(ϖ) = 𝔤^(N)(2)ϖ² − ⟨𝔊^(N)(4)⟩₂ϖ⁴ + O(ϖ⁶)"，⟨𝔊^(N)(4)⟩₂ > 0 ⟹ **𝔥^(N)(ϖ)/ϖ² 在 ϖ→0 时从下方趋于 𝔤^(N)(2)，二阶系数 −⟨𝔊^(N)(4)⟩₂（有限 N 定理级）**。N=∞ 只靠 "also holds when N→∞"。**无数值**。
- **⟹ KAY 声称的「大 λ 已证」（§3.6.2 第一句）本身建立在 N=∞ 的展开断言之上**，而该断言无证明（无谱隙、无换序）。**我方段 I（显式余项 + 认证谱隙 + Davis–Kahan）如闭合，将是这一「已证」部分的首个完整证明**——这是一个可写进论文的附带贡献，措辞应为「we also close a gap in the large-λ argument of [KAY III, §3.6.2], which relied on the unproven N→∞ limit of Theorem 7」。

---

## §6 2024–2026 文献：无人声称证明 [A/B]

**检索**：WebSearch ×3（"Kiessling Altshuler Yuzbashyan Eliashberg upper bound conjecture proof"；论文标题精确检索；"Einstein phonons" conjecture rigorous bound 2026）；arXiv API `abs:"Eliashberg" AND abs:"upper bound" AND abs:"critical temperature"`（按时间降序，共 4 条：2505.02894、2409.02121、2409.00533、2409.00532）；Semantic Scholar `arXiv:2409.02121/citations`（3 条：2607.12213、2409.19562、2409.00532）。

| 候选 | 内容 | 是否碰 Conjecture 1/2 | 等级 |
|---|---|---|---|
| **Elezaby–Abanov, arXiv:2512.20009v2**（Dated June 11, 2026；PRB 10.1103/wpn8-8j2t）"Superconductivity Near a Quantum Critical Point: Bounds on the Transition Temperature in the γ-Model" | **只做 γ-model**。τ := 2πT/g（eq.17），τ_c^γ = 𝔤(γ)。主结果 eq.(44)：**τ_up^γ = Σ_{n≥0}(1/2)^{2n}ζ(γ+2n+1)**，Gershgorin 圆盘 + 对角相似变换 O = diag(1/(n+p))、p=1/2 取最优（App.A）。摘要自陈 "significantly tighter than those currently available in the literature"；Fig.1/4 与 KAY 上界（红线，γ=2 处 ≈5.3 于 τ² 单位）对比。下界部分 "confirm the lower bounds previously established in the literature [3]"（[3]=KAY I）。 | **否**：全文 pp.1–12 无 "Einstein"、无 "conjecture"、无 "dispersive"（仅参考文献列 KAY II/III 标题）。 | A（页图） |
| 同上，γ=2 数值 `[自算]` | τ_up² = ζ(3)+ζ(5)/4+ζ(7)/16+ζ(9)/64+… ≈ 1.2021+0.2592+0.0630+0.0157+0.0039+0.0010+0.0002+… ≈ **1.545**，即 𝔤(2) ≤ 1.545（真值 1.3181，比值 1.172；折 T_c 因子 √1.172 ≈ 1.083，对比 KAY 的 2.034）。**这是对 𝔤(2) 的严格上界，与我方 P4-C 的 [1.31714,1.31820] 包围相比宽得多，但它是纯解析闭式，值得在 P4-C 文献段引用。** | — | B |
| Gnezdilov, arXiv:2505.02894v3（2025-05）"Upper bound on T_c in a strongly coupled electron-boson superconductor" | Yukawa-SYK 模型，T_c ≈ 0.18ω_D√λ 交叉到 0.04ε_F；非严格。 | 否 | A（摘要级） |
| Escamilla–Wang, arXiv:2607.12213（2026-07，Ann. Phys. 10.1016/j.aop.2026.170637）"Analytical solution of the Eliashberg equations…" | Debye 模型近似解析解，非定理。引 KAY III。 | 否 | A（摘要级） |
| Sadovskii, arXiv:2506.19326（2025）"Upper limit … very strong coupling" | 物理稳定性论证，反驳 YAP；非严格。 | 否 | A（摘要级） |

**结论**：截至 2026-09-11，**无任何文献声称证明 KAY II Conjecture 2 / KAY III Conjecture 1，也无任何对标准 Eliashberg（Einstein 或色散声子）的全-λ 严格上界优于 KAY 的 2.034 因子**。唯一的新严格上界工作（2512.20009）限于 γ-model，且其 Gershgorin 技术在 γ=2 只到 1.083 倍。

---

## §7 对 `P4_conjecture2.md` 的具体修订建议

1. **§0 陈述**：补一句「Conjecture 2/1 均以 "There is a critical temperature … which for all … λ>0 is bounded above by" 起句，**含存在性**；本文只证界的那一半，存在性（= ϖ↦𝔥 全 ℝ₊ 单调）仍开放」。
2. **§1 化约**：保留；加注 KAY II Remark 2（p.10）与 Prop 7（p.31, eq.83）作为「KAY 最接近的两处」；措辞用「elementary but, to our knowledge, not stated in [KAY II, III]」。
3. **§1 推论**：把「Tc 由 𝔨(P,T_c)=1/λ（或 sup 定义）给出」改为不依赖 T_c 定义的强形式：「对所有 T ≥ T̃_c，λ𝔨(P,T) < 1，正常态线性稳定；故 𝒮_c 与线 ℒ(λ,P) 的所有交点满足 T < T̃_c」——KAY 本人没有 sup 定义，别替他们定义。
4. **§2 段 I**：级数一般项与算子范数收敛是**我方要证的**（KAY 只有 "± ⋯"）；⟨𝔊(4)⟩₂>0 可引 KAY II Prop 4（有限 N，A 级），N=∞ 的正性由我方 Davis–Kahan 给。
5. **§3 单调性**：改「KAY 只证 ϖ≤√2，N=3,4 靠看图」为「KAY III Thm 5：所有有限 N 在 [0,ϖ_*)、ϖ_*≥√2；N=∞ 经收敛散文；N=3,4 在 λ<λ_* 区间靠看图（p.8, p.26）；全 ℝ₊ 自陈 'would require a different strategy'（p.28）；𝔥/ϖ² 单调性 KAY 零处」。
6. **§5 诚实边界**：加「KAY III §3.6.2 的『大 λ 已证』依赖 Thm 7 的 N=∞ 断言，我方段 I 顺带补证」。

---

## §8 裁决（≤10 行）

1. **KAY III Conjecture 1 逐字是 T_c 形式（严格 "<"）且含「对所有 λ>0 存在 T_c」**；它**不是**逐点 𝔥(ϖ)≤𝔤(2)ϖ²，但界的那一半与逐点严格不等式**等价**（B 级，一行推导）。
2. **「一般 P ⇐ Einstein」的 Jensen 化约在 KAY I–III 中不存在**；KAY 只有 Remark 2「𝔨 不是 𝔥 的平均」与 Prop 7「正部分的 P-平均 ≤ 𝔊_j(2) 谱半径 × ⟨ϖ²⟩」。
3. **选项 (b)「KAY 已陈述化约」为假，不得采用。**
4. **选项 (a)「首证 Conjecture 2 含化约」过强**：Conjecture 2 的存在性半句（单调性）我方不证。
5. **建议措辞 (c)**："We prove the T_c bound asserted in Conjecture 2 of [KAY II] (equivalently Conjecture 1 of [KAY III]): for every P ∈ 𝒫 and λ > 0 the normal state is linearly stable for all T ≥ T̃_c(λ,P) := (1/2π)√(𝔤(2)⟨ω²⟩λ); in particular T_c(λ,P) < T̃_c(λ,P) whenever T_c is defined. The reduction to the Einstein case is an elementary convexity argument not stated in [KAY II, III]; the existence/uniqueness of T_c for all λ > 0 (monotonicity of ϖ ↦ 𝔥(ϖ) on ℝ₊) remains open."
6. 附带可宣称：**KAY III §3.6.2「大 λ 已证」依赖未证的 N=∞ 展开，段 I 补证**；**𝔥 与 𝔥/ϖ² 在认证区间上的单调性是 KAY 零处的净新增**。
7. **无竞争者**（2026-09-11）：2512.20009 只做 γ-model（𝔤(2)≤≈1.545）。
