# P4 · 定理 E（固定 λ、⟨ω²⟩ 时 Einstein 谱极大化 Tc）文献核查 + 弱耦合渐近 h(ϖ)=lnϖ+c₀ 的一手溯源

检索 agent：researcher　日期：2026-09-11　派单：`extremal-lit`（见 `theory/derivations/P4_einstein_extremal.md` §3）
上游：`P4_einstein_extremal.md` v0.1、`P4_conjecture2.md` v0.4、`surveys/P4_KAY_conjecture2_reduction_2026-09-11.md`

**证据等级**：**A** = 原文逐字（页图或 HTML 全文，附页码/式号）；**B** = 我方对 A 级材料的确定性推算（标 `[自算]`）；**C** = 仅摘要级/二手转述/未取到正文。

**本轮一手材料（全部本轮新读）**

| 文献 | 取得方式 | 读到的范围 |
|---|---|---|
| KAY III = arXiv:2409.02121v2 | arXiv PDF 落盘 + 页图 | pp.11–13（Thm 6 / Cor 2 / Thm 7 / Cor 3 / Conj 1 / Remark 1 / Fig.2） |
| Marsiglio, PRB **98**, 024523 (2018) = arXiv:1807.04907v2 | arXiv PDF 落盘 + 页图 | pp.1–7 全（eqs 13/14/26/37） |
| Chubukov–Abanov–Esterlis–Kivelson, arXiv:2004.01281 (Ann. Phys. 417, 168190) | arXiv PDF 落盘 + 页图 | pp.6–13（§III.A 弱耦合 eq.17；§III.B 强耦合 eq.19） |
| Bergmann & Rainer, Z. Physik **263**, 59–68 (1973) | **UMD 课程页扫描件**（见 §5 路线） | pp.59–62（摘要 + §1 + §2 + eqs 1–9 + App.A 结论） |
| Marsiglio & Carbotte, arXiv:cond-mat/0106143（综述） | arXiv PDF 落盘 + 页图 | pp.31–35（§4.1 / **§4.2 Maximum Tc, Asymptotic Limits, and Optimal Phonon Spectra** / §4.3） |
| Moussa & Cohen, PRB **74**, 094520 (2006) = cond-mat/0607832 | ar5iv HTML | 正文相关段 + 参考表 |
| Semenok–Altshuler–Yuzbashyan, arXiv:2407.12922v2（Adv. Mater. 2025, 2507013） | arXiv HTML | 正文 eqs 7/8/9 + **SI §S4（eq. S11/S13）** |
| Mityagin–Sobolev, arXiv:1106.0127 | arXiv HTML | Theorem 1 + 技术综述 |

**未取到正文（付费墙，多路线失败，记 C）**：Allen & Dynes PRB **12**, 905 (1975)；Leavens SSC **17**, 1499 (1975)；Leavens J. Phys. F **7**, 1911 (1977)；Combescot PRB **42**, 7810 (1990)；Carbotte RMP **62**, 1027 (1990)。

---

## §1（a）定理 E 的优先权与新颖性

### 1.1 四类不同的「最优谱」问题——必须先分清约束

文献里至少有四个**互不等价**的优化问题，混用会导致误判优先权：

| 类 | 固定量 | 代表工作 | 最优谱 | 严格性 |
|---|---|---|---|---|
| (1) | **面积** A=∫α²F dω（= λ⟨ω⟩/2 `[自算]`） | Leavens 1975/1977；Bergmann–Rainer 1973 泛函导数为其基础 | δ(ω−ω_opt)，ω_opt 由 ω̄_E≈1 定 | 数值 + 变分启发 |
| (2) | **λ 与 ⟨ω²⟩ 分别固定**（⟹ λ⟨ω²⟩ 固定） | **Moussa–Cohen 2006 明写「已知」并归功 Allen–Dynes 1975** | Einstein δ(ω−√⟨ω²⟩) | **无证明**（见 1.3） |
| (3) | **ω_log** | Combescot 1990（弱耦合精确解） | 零宽度（谱宽增加 ⟹ Tc 单调下降） | 解析（非严格），仅 λ→0 |
| (4) | **ω_max** | Semenok–Altshuler–Yuzbashyan 2024/25 SI §S4 | δ(ω−ω_max) | 变分导数，非严格 |

**定理 E 属于类 (2)。** 与类 (1)(3)(4) 都不是同一命题——尤其类 (1) 的最优频率是被优化出来的自由参数，而定理 E 的 ω_E=√⟨ω²⟩ 被约束钉死。

### 1.2 类 (2) 的命中：Moussa–Cohen 2006 [A]

Moussa & Cohen, *Two bounds on the maximum phonon-mediated superconducting transition temperature*, PRB **74**, 094520 (2006) = cond-mat/0607832，正文（ar5iv HTML，无印刷页码）：

> "0≤Tc≤0.183√(λ⟨ω²⟩)"

紧随其后描述达到上下界的谱：

> "the lower bound corresponding to ω_E→∞ and the upper bound corresponding to ω_E→0"

（并称这两个极限谱是 "unphysical pathologies"，对应最大声子软化与无限硬化。该处引 [allen_dynes] 与 [strong_bounds]。）

**关键句（本轮任务的核心命中）** [A]：

> "It is already known [allen_dynes] that the Einstein spectrum maximizes Tc for a given λ and ⟨ω²⟩."

其中 `[allen_dynes]` = 参考文献 (3)，文献表逐字：

> "P. B. Allen and R. C. Dynes, Phys. Rev. B 12, 905 (1975)."

上下文（该段其余内容，WebFetch 因单条引文 125 字符上限只能转述，记 C 级转述）：作者接着说即便同时固定 λ 与 ⟨ω²⟩ 并指定 μ\*，**下界仍为 0**，原因是一个病态双峰谱（一峰趋零频、一峰趋无穷频）；为得到非零下界他们进一步加上 ω_max 限制。

⟹ **定理 E 的命题内容在 2006 年已被当作「已知」印在 PRB 上**，且归功于 Allen–Dynes 1975。

### 1.3 但 Allen–Dynes 1975 原文未取到——归属链有一个未闭合环节

本轮四条路线全败：APS 付费墙；OSTI 仅元数据；ADS 无全文源；Allen 个人页 `felix.physics.sunysb.edu/~allen/` 连接被拒（socket closed）。**故「AD 1975 原文是否真的陈述了类 (2) 最优性」目前只有 Moussa–Cohen 的转述（A 级的「MC 怎么说」，不是 A 级的「AD 怎么写」）。**

**但 AD 拟合公式本身蕴含该结论** `[自算，B]`。AD 的形状因子（经 Marsiglio–Carbotte 综述 eq.93 与 repo `P3_allen_dynes_original_2026-09-09.md` 交叉确认的形式）：

  T_c = (f₁f₂ ω_log/1.20)·exp[−1.04(1+λ)/(λ−μ\*(1+0.62λ))]，
  f₁ = [1+(λ/Λ₁)^{3/2}]^{1/3}（只依赖 λ,μ\*），
  f₂ = 1 + λ²(r−1)/(λ²+Λ₂²)，r := ω₂/ω_log ≥ 1，Λ₂ = 1.82(1+6.3μ\*)·r。

固定 λ、μ\*、ω₂=√⟨ω²⟩，则 ω_log=ω₂/r，指数项与 r 无关，于是

  T_c ∝ ω_log·f₂ = ω₂·[f₂/r] = ω₂·[a(r) + (1−a(r))/r]，　a(r) := λ²/(λ²+Λ₂(r)²) ∈ (0,1)。

因 a<1 且 a(r) 随 r 递减（Λ₂ ∝ r），括号内对 r 严格递减 ⟹ **T_c 在 r=1 取最大，r=1 ⟺ ω_log=ω₂ ⟺ Einstein 谱**。∎`[自算]`

这与 Jensen 方向一致（ω₂ ≥ ω_log 恒成立，等号仅对 δ 谱），并与 λ→∞ 时 f₂ → r（完全补偿、渐近与形状无关）自洽。

⟹ **定理 E 的「数值内容」不新**：它是 AD 1975 拟合公式的一行代数推论，并且已被 MC 2006 当作既有事实引用。

### 1.4 类 (1)「optimum spectrum」的一手陈述 [A]

Marsiglio & Carbotte, cond-mat/0106143, **p.33**：

> "The reasoning above leads naturally to the concept of an optimum spectrum, first determined by Leavens [196], and elaborated upon in Refs. [197–201]. In an optimum spectrum calculation, one imagines having a fixed area of α²F(ν), and asks at what frequency it would best be situated in order to optimize some particular property."

> "An appropriate scaling of the linearized Eliashberg equations for an Einstein phonon spectrum with frequency ω_E leads to the result
> T_c/A = f(ω̄_E, μ\*),　(98)
> where A is the area, ω̄_E ≡ ω_E/A, and f is a universal function of ω̄_E, to be determined numerically for each choice of μ\* … The result is a curve with a maximum at ω̄_E ≈ 1; placing a spectral function at this frequency will yield the maximum T_c. This procedure yields a result,
> T_c ≤ A c(μ\*),　(99)
> where c(μ\*) is a function of μ\* shown in Fig. 19."

**三个逐字要点**：
1. 类 (1) 的优化是**在 δ 族内部对频率**做的（"one imagines having a fixed area … and asks at what frequency it would best be situated"）——**「δ 优于任意形状」这一步 MC 综述没有证明，也没有声称证明**。
2. **c(μ\*) 的数值 MC 正文未给**（只说见 Fig. 19）；网络二手来源给 c(0)≈0.25、经验关联 Tc≈0.148A（Leavens–Carbotte 1974, JLTP 14, 195），本轮**均未取到一手数值，记 C，不得入 key_numbers**。
3. 同页 **p.32** 的 λ→∞ 渐近（本轮再次逐字确认，与 repo 既有记录一致）：
> "T_c = 0.183√λ ω_E,　(94) … as λ → ∞. In obtaining Eq. (94), an Einstein spectrum has been assumed (this is not required), which in turn is characterized by two parameters: the weight, A ≡ λω_E/2, and the frequency, ω_E."

⟹ 由 A = λω_E/2 得 **A = λ⟨ω⟩/2**（Einstein 时 ⟨ω⟩=ω_E）`[自算]`，即类 (1) 固定的是**第一矩** λ⟨ω⟩，类 (2) 固定的是**第二矩** λ⟨ω²⟩。两者本质不同。

### 1.5 Bergmann–Rainer 1973：泛函导数——**不是**形状最优性 [A]

Z. Physik **263**, 59–68 (1973)，扫描件逐字。**p.59 摘要**：

> "For this purpose the functional derivative δT_c/δα²F(ω) is calculated. It is shown that all frequency regions of α²F(ω) yield a positive contribution to T_c and that the most effective range covers frequencies, slightly above 2πT_c."

**p.60**：

> "In all examples considered here the functional derivative is always positive. Furthermore we can prove this property rigorously for μ\*=0 (see Appendix A). Therefore we are convinced that in thermal equilibrium phonons never have a pair-breaking effect on the transition temperature."

**p.62**（低频行为，eq.9 之后）：

> "Thus the functional derivative goes to zero linearly at sufficiently low frequencies. … At T_c the superconductor cannot feel the dynamics of lattice oscillations with frequencies small compared to T_c. Such lattice modes act more like static lattice deformations which have no influence on T_c according to Andersons theorem."

> "In Appendix A we prove rigorously (for μ\*=0) that the functional derivative of ρ has always a positive sign."

**判读**：BR 给的是 **δT_c/δα²F ≥ 0**（单调性，μ\*=0 时严格证明）+ 峰位在 **"slightly above 2πT_c"**。这**不是**任何固定矩约束下的形状最优性陈述。二手材料常说的「峰在 ω ≈ 7k_BT_c」与 BR 原文 "slightly above 2πT_c"（2π≈6.28）**量级一致但不是原文措辞**——引用时应用 BR 原话或引 MC 综述 p.33 的 "a factor of order 10 times the critical temperature"（MC 自己的说法，与 BR 原话也不同，见 §4 矛盾表）。

BR 的 δρ/δα²F 闭式在 **p.62 eq.(8)**；低频渐近 eq.(9)；δT_c/δα²F 与 δρ/δα²F 的关系 **p.61 eq.(6)**；Hermitian 本征值问题 **p.61 eq.(7)** —— **eq.(7) 正是我方 ℌ 的 1973 年版本**（含 μ\* 与 pair-breaking 参数 ρ），对 P4-D 的 μ\* 算子序层也有用。

MC 综述 **p.32** 给了 BR 泛函导数的一个近似闭式（引 Ref.[195]）：

> "δT_c/δα²F(Ω) = 1/(1+λ) Σ_{n=1}^∞ 4Ω̄/(Ω̄²+4π²n²),　(96)
> where Ω̄ ≡ Ω/k_BT_c"

并在 **p.33** 说：

> "It illustrates that the optimal phonon frequency lies at some finite frequency (i.e. non-zero, and non-infinite), which is a factor of order 10 times the critical temperature."

### 1.6 类 (3) Combescot 1990：固定 ω_log 的宽度单调性 [C，摘要级]

R. Combescot, PRB **42**, 7810–7824 (1990), *Critical temperature of superconductors: Exact solution from Eliashberg equations on the weak-coupling side*。摘要级内容：λ, μ\*→0 极限下对一般 α²F(ω) **精确**求解 Tc；ω_log 是「最优特征频率」；**残余的形状依赖被详细研究，给出「固定 ω_log 时谱宽增大导致 Tc 一致下降」**；弱耦合下 Tc 只依赖谱的大尺度结构。

⟹ 这是**定理 E 在弱耦合端、换一个固定量（ω_log 而非 ⟨ω²⟩）的先例**，且是解析而非仅数值。**正文未取到（APS 付费墙），不得当 A 级引。** 注意 repo 既有记忆中的「Combescot PRB 51, 11625 = 弱耦合」是**错的**：51,11625 (1995) 是 *Strong-coupling limit of Eliashberg theory*；弱耦合那篇是 **42, 7810 (1990)**。本条**订正既有记忆**。

### 1.7 类 (4) Semenok–Altshuler–Yuzbashyan 2024/2025 [A]

arXiv:2407.12922v2（Adv. Mater. 2025, 2507013），正文：

> "Tc < 0.18√(λ⟨ω²⟩)"　(7)　（引 [8, 29, 53]）
> "⟨ω²⟩ = (2/λ)∫₀^∞ α²F(ω) ω dω"　(8)
> "Tc < 0.32 ω_max"　(9)
> "λ\* = 3.69"（"attained for Einstein phonons"）；"T_c^E = 0.3175 ω_max"（引 [29, 46]）

**SI §S4**（eq. S13 之后）：作者对 eq.(S13) 右端关于形状函数 P(x) 取变分导数，得到

> "we find that it is maximal for P(x) = δ(x−1)"

**但约束是 supp P ⊂ [0,1]（即 ω ≤ ω_max）、归一化 ∫P=1（eq. S11），被极大化的是 λ⟨ω²⟩ 那个比值（含对 τ 的极大化），目的是由 (7) 推出 (9)。** ⟹ 这是**固定 ω_max 的形状优化**，**不是**固定 λ 与 ⟨ω²⟩；且方法是变分导数（一阶条件），**非严格**。**不构成定理 E 的先例，但构成「Einstein 在某约束下极值」这一直觉在该团队内部的近期背书**——而 KAY 三部曲的两位作者（Altshuler、Yuzbashyan）正是该文作者，这一点在论文里值得点名（他们在非严格侧做过 δ-最优的变分论证，却在严格侧的 KAY II/III 没有写 Jensen 化约）。

### 1.8 KAY I–III：零处 [A，未见]

前轮（`surveys/P4_KAY_conjecture2_reduction_2026-09-11.md` §2）已对 KAY II pp.3–42、KAY III pp.1–30 逐页核对，**无任何「Einstein is extremal」/λ₁ 的 Jensen 化约/"suffices to prove the Einstein case"**。本轮复核 KAY III pp.11–13 未发现新增。

一个**结构性旁证**（`[自算]`，可写进论文）：KAY II Conjecture 2 的 T̃_c(λ,P) = (1/2π)√(𝔤(2)⟨ω²⟩λ) 与 KAY III Conjecture 1 的 T̃_c(λ,Ω) = (Ω/2π)√(𝔤(2)λ) **是同一函数在 ⟨ω²⟩=Ω² 处的取值**——即 KAY 的两个猜想天然与「界只经 ⟨ω²⟩ 进入、Einstein 是等号情形」相容，但 KAY 从未把它表述成极值原理，也从未指出 Conj.1 ⟹ Conj.2。

### 1.9 裁决（(a)）

- **作为「严格定理」：新。** 检索未发现任何对类 (2) 陈述的证明（严格或半严格）。MC 2006 的 "It is already known" 之后没有给证明，只给了归属；AD 1975 即便真的写了，也只能是数值/拟合层面（其全文本轮未取到）。
- **作为「陈述」：不新。** MC 2006 逐字印过，归 AD 1975。**论文里必须引 MC 2006 这一句**，不得写「据我们所知无人陈述」。
- **作为「数值内容」：不新。** AD 拟合公式的一行代数即给出同一结论（§1.3 `[自算]`），而 AD 公式本身是 200+ 次 Eliashberg 数值解的拟合。
- **真正的增量**：(i) 首个**严格**证明；(ii) 覆盖**全部 λ>0**（不只 AD 拟合的适用域，也不只 Combescot 的 λ→0）；(iii) 对**任意** P∈𝒫（只需 ⟨ω²⟩<∞，不需紧支、不需单峰）；(iv) 把 KAY Conj.1 ⟹ Conj.2 的化约显式化。
- **建议措辞**："That the Einstein spectrum should maximize Tc at fixed λ and ⟨ω²⟩ has been asserted in the physics literature — see [Moussa–Cohen 2006], who call it 'already known' and attribute it to [Allen–Dynes 1975] — and is implied by the Allen–Dynes interpolation formula; we are not aware of any proof. Theorem E supplies one."

---

## §2（b）弱耦合 / ϖ→∞ 端：h(ϖ) = ln ϖ + c₀ 的一手溯源

### 2.1 KAY III Theorem 6 到底证了什么 [A，逐字，p.11]

> **Theorem 6**: *The eigenvalues 𝔥^(N)(ϖ) are analytic about ϖ = ∞, with*
> 𝔥^(N)(ϖ) = 𝔥^(N)(∞) − B_N·(1/ϖ²) + 𝒪(1/ϖ⁴)　(22)
> *where*
> 𝔥^(N)(∞) = −1 + 2 Σ_{n=0}^{N−1} 1/(2n+1)　(23)
> *and*
> B_N = [ Σ_{n=0}^{N−1}Σ_{m=0}^{N−1} ((n−m)²+(n+m+1)²)/((2n+1)(2m+1)) − Σ_{ℓ=0}^{N−1} (2/(2ℓ+1)²)Σ_{k=1}^{ℓ}k² ] / [ Σ_{j=0}^{N−1} 1/(2j+1) ] ≥ 1.　(24)
> *Thus, as ϖ → ∞ the N-th upper approximation L_E^(N)(ϖ) to the critical curve given by ϖ ↦ L_E(ϖ) converges downward to λ_N = 1/𝔥^(N)(∞). Moreover, both 𝔥^(N)(∞) and B_N are strictly monotonically increasing with N, diverging to ∞ as N → ∞.*

紧接（p.11）：

> "Only the existence of the 𝔥^(N)(∞) follows as special case of the analogous result about the T → 0 limit of the eigenvalues 𝔨^(N)(P,T) of the Eliashberg model with dispersive phonons that we proved in [KAYb]. Theorem 6 in full therefore will be proved in this paper. It establishes that the L_E^(N)(ϖ) are asymptotic to strictly monotonically decreasing functions when ϖ ∼ ∞."

**Corollary 2**（p.11，eq.25）：as λ→λ_N⁺，T_c^(N)(λ,Ω) ∼ (Ω/2π)√( (1/B_N)(1/λ_N − 1/λ) )。

**三个必须写清的限定** [A]：
1. **Theorem 6 只对有限 N。** 它给的是「固定 N、ϖ→∞」的极限 𝔥^(N)(∞)（有限），**不是**「N=∞、ϖ→∞」。
2. **KAY III 全文没有 h(ϖ) ~ ln ϖ 这一句**，也没有任何 N=∞ 的大-ϖ 渐近。ϖ→∞ 与 N→∞ 的换序在 KAY 处**完全没碰**（𝔥^(N)(∞) 与 B_N 都 →∞）。
3. **「λ_N ~ 1/ln N」是上界曲线与 λ 轴交点的速率**（由 (23) 得 λ_N = 1/𝔥^(N)(∞)），**不是**固定 ϖ 处本征值的收敛速率——这是既有记忆里已经标记过的检索陷阱，本轮由 Thm 6 正文再次确认。

**𝔥^(N)(∞) 的显式渐近** `[自算，B]`：Σ_{n=0}^{N−1}1/(2n+1) = H_{2N−1} − ½H_{N−1} = ½ln N + ln 2 + γ/2 + O(1/N)，故

  **𝔥^(N)(∞) = ln N + (γ + ln 4 − 1) + O(1/N) = ln N + 0.9635101… + O(1/N)。**

（KAY III 只说 "diverging to ∞"，**未给这个常数**。）

### 2.2 已发表的弱耦合 Tc 前因子：1.13 是错的，正确的差一个 1/√e

**Marsiglio, PRB 98, 024523 (2018) = arXiv:1807.04907v2** [A]：

- **p.3, eq.(13)**（BCS）：
> "T_c = 1.13 ω_E exp(−1/λ)　(13)"
- **p.3, eq.(14)**（加入 Z 重整）：
> "The inclusion of the renormalization, Z, modifies this equation to read
> T_c = 1.13 ω_E exp(−(1+λ)/λ).　(14)"
- **p.3 右栏（关键警告）**：
> "One can immediately write this like Eq. (13) but with reduced pre-factor 1.13e^{−1}. This is *not* what is meant when we stated that the pre-factor in Eliashberg theory is actually modified from the BCS result — but rather an additional change occurs."
- **p.4, eq.(26)**（Z≡1 的「un-renormalized」情形）：
> "T_c = (1.13/√e) ω_E exp(−1/λ)　(26)
> in contrast to Eq. (13)."
- **p.6, eq.(37)（含 Z 重整——这才是与 KAY 的 ℌ 对口的那一个）**：
> "T_c = (1.13/√e) ω_E exp(−(1+λ)/λ).　(37)"
> "As previously mentioned, Eq. (37) can of course be written with a −1/λ in the exponential, along with a prefactor denominator of e^{3/2} instead of √e. However, the present form more explicitly shows the role of the 'normal-state' renormalization that gives rise to the usual 1+λ factor, *along* with the not-so-usual √e denominator in the pre-factor."
- 归属（p.3）：
> "McMillan and others adopted phenomenological pre-factors, whose justification is now more readily understood after Karakozov et al.²⁸ solved the gap equation on the real axis with an iterative method and obtained the result that T_c attains a pre-factor significantly different from that obtained with BCS theory."

**独立第二来源：Chubukov–Abanov–Esterlis–Kivelson, arXiv:2004.01281, p.12, eq.(17)** [A]：

> "The only substantial difference between the Eliashberg and BCS theories at weak coupling is that the latter requires a high-energy cutoff, which sets the pre-exponential factors for T_c and Δ, while in ETh the cutoff is effectively provided by the frequency dependence of V(Ω). As a consequence, both T_c and Δ are obtained within ETh with the exact prefactors, as has been discussed in several papers, using different computational tools¹⁸⁻²⁶. The result is
> T_c = 1.13 e^{−1/2} Ω₀ e^{−(1+λ)/λ} = 0.252 Ω₀ e^{−1/λ}　(17)"

同文 **p.12, eq.(19)** 的强耦合端（与 KAY 的 C_∞ 对照）："T_c ≈ 0.1827 α"，α=√(λ)Ω₀。

⟹ **两个独立 A 级来源一致**：含 Z 重整的 Einstein 弱耦合精确前因子是 **1.13 e^{−3/2}**（= 2e^γ/π · e^{−3/2}）。

### 2.3 由此定出 c₀ —— 对 implementer `extremal-asymptotics` 的可证伪预测 `[自算，B]`

设 h(ϖ) = ln ϖ + c₀ + o(1)。线性化失稳判据 λh(ϖ_c)=1 给 ϖ_c = e^{1/λ−c₀}，而 ϖ=Ω/2πT ⟹

  **T_c = (Ω/2π)·e^{c₀}·e^{−1/λ}。**

与 eq.(17)/eq.(37) 的 T_c = (2e^γ/π)e^{−3/2}Ω e^{−1/λ} 对比（1.13 = 2e^γ/π）：

  e^{c₀} = 2π·(2e^γ/π)·e^{−3/2} = 4e^{γ−3/2} ⟹ **c₀ = γ + ln 4 − 3/2 = 0.4635101…**

对照「朴素硬截断」直觉（把 ϖ 当成 Matsubara 截断数 N 代进 (23)）会给 **c₀ = γ + ln4 − 1 = 0.9635101…**，正差 **1/2**，即

  **h(ϖ) ≈ 𝔥^(N)(∞)|_{N = ϖ/√e}**　`[自算，B]`

——**√e 的物理意义就是「有效 Matsubara 截断是 ϖ/√e 而非 ϖ」**，这正是 Marsiglio p.3 那句 "an additional change occurs" 与 eq.(24)（I₀ = λI₀² − I₀/2）所刻画的软截断修正。

**给 implementer 的判据**：高精度拟合若得 c₀ ≈ 0.4635，则 (i) 数值管线与已发表 Eliashberg 弱耦合结果一致，(ii) 常数有闭式可猜、可作为段 C-II 的目标；若得 ≈ 0.9635，则说明实现里丢了软截断（大概率是 −𝔎₁ 或 n=0 行的处理），**是 bug 信号而非新物理**。两者差 0.5，在 ϖ~10³ 时对应 h 的相对差 ~7%，float64 完全可分辨。

### 2.4 是否存在任何**严格**（带误差界）的大-ϖ 渐近？ [A/C，未见]

- **KAY I–III：无。** Thm 6 是有限 N 的（且 KAY III p.9 自陈 "we only know that the sequence L_E^(N)(ϖ) converges downward to L_E(ϖ) when N→∞, but don't know how fast"，见既有记忆）。
- **Marsiglio 2018 / CAEK 2020 / Karakozov et al. 1975–76：非严格**（迭代解、digamma 渐近展开、忽略 O(λ²) 项；Marsiglio eq.(21) 原文自陈 "This equation is clearly an approximation"，p.4）。
- **Combescot 1990：「exact solution」是物理意义的精确渐近，非定理**（正文未取到，记 C）。
- **结论**：**该 Matsubara 矩阵族顶本征值的 ϖ→∞ 渐近，目前没有任何带误差界的严格结果。** 这是段 C-II 的净新增空间，也是本项目在 (b) 上的真实增量。

---

## §3（c）数学侧：顶本征值的凹性 / 标度极限 / 比较算子

### 3.1 一般理论给的是**凸**性，不是凹性——没有免费午餐

- λ₁ 是 Rayleigh 商的 sup，故对算子**线性**依赖时 λ₁ 是参数的**凸**函数（Rayleigh 极大原理，教科书级）。要得到凹性只能靠 x ↦ ℌ(√x) 的**算子凹性**（Loewner），而 repo 的 `loewner_probe.py` 已数值排除（λ_min(D)≈−0.47）。**⟹ 抽象路线关闭是可复核的，应在论文里明写这一句。**
- 非仿射参数依赖下凸性何时成立的判据：J.-B. Hiriart-Urruty 等，*When is the greatest eigenvalue of a parametrized symmetric matrix a convex function of the parameter?*, Linear Algebra Appl. **169** (1992) 129–133。**这是判凸的判据，不给凹性**，但可用来说明「凹性必须来自模型的特殊结构」。[C，摘要级]
- 「本征值和的凹/凸性」（Hermitian V ↦ Σ 负本征值 是凹的）：arXiv:math/0112279（Concavity of Eigenvalue Sums and the Spectral Shift Function）与 J. Stat. Phys. 的 *Convexity and concavity of eigenvalue sums*。**对象是本征值之和，不是单个 λ₁，与我方不对口。** [C]

⟹ **(c) 第一条裁决：文献里没有可直接搬运的「Matsubara 核顶本征值凹性」定理。** (C\*) 式那条「二阶微扰 + 谱隙」路线是自建的，必须自证。

### 3.2 标度极限 n = ϖs：确实把 ℌ(ϖ) 变成差核积分算子 `[自算，B]`

ℌ 的矩阵元由 ⟦k⟧(ϖ) = ϖ²/(k²+ϖ²) 与 (2n+1)^{−1/2} 型权组成（KAY II eqs.(40)–(42)）。令 n = ϖs、m = ϖt：

  ⟦n−m⟧(ϖ) = ϖ²/((n−m)²+ϖ²) = 1/(1+(s−t)²)，　Σ_n → ϖ∫ds，　(2n+1)^{−1/2} → (2ϖs)^{−1/2}。

⟹ ℌ(ϖ) 的标度极限是 L²(0,∞) 上核为 **[1/(1+(s−t)²)] · (st)^{−1/2}** 型的算子（加上 −ℌ₁ 的对角部分）。**分母 1+(s−t)² 是纯差核**，权 (st)^{−1/2} 是齐次的。

这与下面这支数学文献的核**同型**：

### 3.3 可直接对接的数学文献（本轮最有价值的 (c) 线索）

- **Mityagin & Sobolev, arXiv:1106.0127, *A family of anisotropic integral operators and behaviour of its maximal eigenvalue*** [A，HTML 全文]：L²(ℝ) 上紧自伴积分算子 **K_β**，核
  K_β(x,y) = (1/π)·1/(1+(x−y)² + β²Θ(x,y))，Θ 对称非负、**齐次 γ≥1 次**。
  **Theorem 1**（逐字）：记 M_β 为 K_β 的最大本征值，则
  > "lim_{β→0} β^{−2/(γ+1)}(1 − M_β) = λ₁"，
  > "the rescaled eigenfunctions … converge in norm to φ₁ as β → 0"，
  其中 λ₁ 是**模型算子 A = |D_x| + 2^{−1}θ(x)** 的最低本征值，θ(x)=Θ(x,x)。
  **技术**：Widom 的差核方法（原文 "Some ideas of this paper are used in the proof"）；**拟微分算子符号演算**（对去对称化后的核，符号近似为 1 − α a(x,ξ)，a(x,ξ)=|ξ| + 2^{−1}θ(x)）；正性改进半群（Reed–Simon IV）给顶本征值非退化 + 本征函数正；Schur 检验与 Cordes 型 L² 有界性判据（Prop. 15–16）。
  **物理来源**：Θ(x,y)=(x²+y²)² 一例 "was suggested by P. Krotkov and A. Chubukov … as a simplified model of high-temperature superconductivity"（把 Matsubara 求和丢掉后的 Eliashberg 型间隙方程）。
  **⟹ 这是我们能找到的、与段 C-II 最同构的严格工作**：同样是「小参数 → 差核 1/(1+(x−y)²) 的顶本征值渐近」，同样给出**模型算子 |D| + 势**，同样有**误差阶**与**本征函数收敛**（后者正是 (C\*) 里 ‖(I−P)ℌ′v‖ 需要的东西）。

- **Mityagin, arXiv:0803.3159, *An anisotropic integral operator in high temperature superconductivity*** [A，摘要级]：同一族的前身，严格证明 Krotkov–Chubukov 猜的 **τ(a) ~ a^{2/5}**（E₀(a,T)=1 的解 T(a)=1−τ(a)），见其 Theorem 8 与 Proposition 10；摘要自陈 **Landau–Pollack–Slepian 加性不确定性原理**在分析中起重要作用（**注意**：Mityagin–Sobolev 1106.0127 全文**未提** LPS，见 §4 矛盾表）。

- **经典工具**：
  - H. Widom, *Extreme eigenvalues of translation kernels*, Trans. Amer. Math. Soc. **100** (1961), 252–262。
  - H. Widom, *Extreme eigenvalues of N-dimensional convolution operators*, Trans. Amer. Math. Soc. **106** (1963), 391–414（**AMS 免费 PDF**：`ams.org/journals/tran/1963-106-03/S0002-9947-1963-0145294-7/S0002-9947-1963-0145294-7.pdf`）。框架：k∈L¹、k(−x)=k(x)、域 tΩ 膨胀、t→∞；主阶由符号（Fourier 变换）最大值 M 给出，μ_{m,t}→M，**论文主体是收敛速率**，速率由符号在极大点附近的光滑性/衰减决定（对数型修正正是从这里来的）。
  - H. Kesten, *On the extreme eigenvalues of translation kernels and Toeplitz matrices*, J. Analyse Math. **10** (1962), 117–138。
  - A. Sobolev, *Extreme eigenvalues of an integral operator*, arXiv:1609.02052（现代版，兼引 Slepian 的 prolate 函数）。

### 3.4 Hilbert / Carleman 型比较算子 [A/C，未见]

- **未找到**任何把 Hilbert 矩阵或 Carleman 算子用作线性化 Eliashberg 顶本征值上下界比较对象的文献。
- KAY 自己用的是 **ζ 型 Schur / Hilbert–Schmidt 界**（KAY I Prop 5/6 eqs.(89)/(98)、App.A Prop A/B/C），**不是** Hilbert 矩阵比较。
- 唯一的「新严格上界」竞争者 Elezaby–Abanov arXiv:2512.20009v2 用的是 **Gershgorin 圆盘 + 对角相似**（既有记忆已记），也不是 Hilbert 型。
- **⟹ Hilbert 型夹逼若可行，是空白地**；但注意 ℌ 的核是 1/(k²+ϖ²) 型（平方可和），Hilbert 矩阵 1/(n+m+1) 只平方可和到边缘（非 HS），**两者衰减阶不同，直接 Loewner 夹逼未必成立**，须先做逐元比较验证。

---

## §4 矛盾与不确定点（必须在论文/工作区如实保留）

| # | 矛盾/不确定 | 双方 | 处理建议 |
|---|---|---|---|
| M1 | **AD 1975 是否真的陈述了类 (2) 最优性** | Moussa–Cohen 2006 逐字说 "already known [AD]"；AD 原文本轮四路线全败未取到 | 引用写成「MC 2006 称其为已知并归功 AD 1975」，**不要写「AD 1975 证明了/陈述了」**；列为待补 |
| M2 | **BR 泛函导数峰位** | BR 1973 p.59 原文 "slightly above 2πT_c"；MC 综述 p.33 "a factor of order 10 times the critical temperature"；网络二手 "ω ≈ 7k_BT_c" | 三者量级相容但措辞互不相同；**只引 BR 原话**，其余标转述 |
| M3 | **Combescot 弱耦合论文的卷号** | 既有记忆记 PRB 51, 11625 (1995)；本轮检索确认 51,11625 是 *Strong-coupling limit*，弱耦合是 **PRB 42, 7810 (1990)** | **订正记忆**；两篇都未取到正文 |
| M4 | **LPS 不确定性原理的角色** | Mityagin 0803.3159 摘要称其 "plays important role"；Mityagin–Sobolev 1106.0127 全文**未出现** LPS | 两篇技术路线不同（前者 LPS，后者 Widom+PDO）；引用时分别标注，勿混 |
| M5 | **Leavens 界的常数 c(μ\*=0)** | 二手给 0.25（另有 Tc≈0.148A 的经验关联，来自 Leavens–Carbotte JLTP 14,195 (1974)）；MC 综述 p.33 只说「见 Fig.19」未给数 | **一手数值未取到，不得入 key_numbers**；若要用须取 Leavens 1977 或 MC Fig.19 |
| M6 | **c₀ 的闭式是「预测」不是「已证」** | 它依赖 Marsiglio eq.(37)/CAEK eq.(17) 的前因子为**精确**（两文都自称 exact prefactor），以及 h 的展开确为 lnϖ+c₀+o(1) | 标 B 级预测；数值一旦复现 0.4635 才升级 |

---

## §5 方法笔记（可复用）

- **Bergmann–Rainer 1973 正文的取件路线**：`physics.umd.edu/courses/Phys798C/AnlageSpring22/Bergmann-1973-The%20sensitivity%20of%20the%20transitio.pdf`（UMD Phys798C 课程页扫描件，CCITTFax 扫描 PDF）。**WebFetch 的文本抽取必然失败**（返回 "raw binary"），但落盘后用 `Read` 带 `pages` 读**页图**即可逐字。这再次证实既有记忆里的两条：① 大学课程页扫描件是取 pre-arXiv 老论文的首选路线；② PDF 一律「WebFetch 落盘 → Read pages」。
- **arXiv PDF 同法**：KAY III、Marsiglio 2018、CAEK 2020、Marsiglio–Carbotte 综述，本轮全部经此法逐字读到。
- **ar5iv HTML** 对 2006 年的 cond-mat 论文可用（Moussa–Cohen），但 WebFetch 的下游模型有**单条引文 125 字符上限**，长段落取不全；对策是分多次、按关键词定点取短句。
- **arXiv 新版 HTML**（`arxiv.org/html/…v2`）对 2024 年后的论文可用，且**能读到 SI**（Semenok 2407.12922 的 §S4 即由此取得）。

---

## §6 建议的后续检索方向

1. **Allen & Dynes PRB 12, 905 (1975) 全文**（闭合 M1）：机构订阅 APS / 馆际互借 / 检查 Allen–Mitrović, *Theory of Superconducting Tc*, Solid State Phys. **37**, 1–92 (1982)（同一作者的扩展综述，可能复述类 (2) 最优性且更易取得）。
2. **Combescot PRB 42, 7810 (1990) 全文**（闭合 M3/类 (3)）：其「固定 ω_log 下谱宽增大 ⟹ Tc 一致下降」如果是可写下的解析不等式，**可能直接给段 C-II 一个弱耦合端的解析对照**。
3. **Karakozov, Maksimov, Mashkov (ZhETF 1975/76)**：√e 的原始出处，Marsiglio ref.[28]。若其推导是对同一 Matsubara 本征值问题做的渐近，则是 c₀ 的最早来源。
4. **Leavens J. Phys. F 7, 1911 (1977) + Carbotte RMP 62, 1027 (1990) §optimum**：取 c(μ\*) 的一手数值（闭合 M5）。
5. **Widom 1961/1963 + Mityagin–Sobolev 1106.0127 全文精读**：把 Theorem 1 的证明骨架（PDO 符号 1−α(|ξ|+θ/2)、正性改进半群给非退化、Schur/Cordes 给余项）逐条对照我方标度极限，判断能否直接移植出 **h(ϖ) 的带误差界渐近 + 本征函数 φ(n/ϖ)/√ϖ 的收敛**——后者正是 (C\*) 里 ‖(I−P)ℌ′v‖ 与谱隙 g 在大 ϖ 处的来源。
6. **Krotkov–Chubukov 原始物理论文**（PRL/PRB，"Non-Fermi liquid and pairing in electron-doped cuprates"）：看他们如何从 Eliashberg 方程丢掉 Matsubara 求和得到那个积分算子——我方是**保留**求和做标度极限，路线不同但可借其记号。

---

## §7 裁决（≤12 行）

1. **定理 E 作为严格定理：新。** 检索未发现任何证明（严格或半严格）。
2. **定理 E 作为陈述：不新。** Moussa–Cohen PRB 74, 094520 (2006) 逐字写 "It is already known [Allen & Dynes 1975] that the Einstein spectrum maximizes Tc for a given λ and ⟨ω²⟩."，**论文必须引这句，不得声称首次陈述**。
3. **定理 E 的数值内容：不新。** Allen–Dynes 拟合公式的一行代数即给同一结论（f₂/r 对 r=ω₂/ω_log 严格递减，§1.3 `[自算]`）。**AD 1975 原文本轮未取到**，归属链最后一环仍是 MC 的转述（M1）。
4. 文献里的「optimum spectrum」（Leavens 1975/77，MC 综述 p.33 eqs.98/99）固定的是**面积 A=λ⟨ω⟩/2**（第一矩），Bergmann–Rainer 1973 给的是**泛函导数正性 + 峰在 "slightly above 2πT_c"**，Semenok–Altshuler–Yuzbashyan 2024/25 SI §S4 固定的是 **ω_max**——**三者都不是定理 E 的约束**。Combescot PRB 42, 7810 (1990) 的「固定 ω_log、谱宽↑ ⟹ Tc↓」是最接近的类比，但只在 λ→0 且正文未取到。
5. **KAY I–III 零处**（前轮逐页 + 本轮 pp.11–13 复核）。可点名：KAY 的两位作者在非严格侧（2407.12922 SI）做过 δ-最优的变分论证，却在严格侧没写 Jensen 化约。
6. **h(ϖ) 的最佳已发表渐近形式**：**没有任何文献直接写 h(ϖ)=lnϖ+c₀**。KAY III **Theorem 6 (p.11, eqs.22–24) 只对有限 N**，给 𝔥^(N)(∞)=−1+2Σ_{n<N}1/(2n+1)（⟹ = ln N + γ + ln4 − 1 + O(1/N) `[自算]`），λ_N=1/𝔥^(N)(∞)~1/ln N 是**交点速率不是收敛速率**。
7. **但 c₀ 被已发表的弱耦合前因子唯一确定** `[自算，B]`：由 Marsiglio PRB 98, 024523 **eq.(37)** T_c=(1.13/√e)ω_E e^{−(1+λ)/λ} 与 CAEK arXiv:2004.01281 **eq.(17)** T_c=1.13e^{−1/2}Ω₀e^{−(1+λ)/λ}=0.252Ω₀e^{−1/λ}（两独立 A 级来源一致），得
   **c₀ = γ + ln 4 − 3/2 = 0.4635101…**，等价于「有效 Matsubara 截断 N_eff = ϖ/√e」。
8. **给 `extremal-asymptotics` 的判据**：拟合出 0.4635 ⟹ 管线与已发表 Eliashberg 一致；拟合出 0.9635（= γ+ln4−1，硬截断值）⟹ **丢了软截断，是 bug 信号**。两者差 0.5。
9. **严格渐近（带误差界）：全域为零。** KAY 自陈 "don't know how fast"；Marsiglio eq.(21) 自陈 "clearly an approximation"。段 C-II 是真空地。
10. **大-ϖ 严格论证的具体线索（最有价值）**：标度 n=ϖs 把 ⟦n−m⟧ 变成**差核 1/(1+(s−t)²)** `[自算]`，与 **Mityagin–Sobolev arXiv:1106.0127** 的 K_β 族**同型**；其 **Theorem 1** 给出「顶本征值渐近 + 重标本征函数按范数收敛 + 模型算子 A=|D_x|+θ(x)/2」，技术是 **Widom 差核方法 + 拟微分符号演算 + 正性改进半群 + Schur/Cordes 余项界**。**本征函数收敛正是 (C\*) 里 ‖(I−P)ℌ′v‖ 与谱隙 g 在大 ϖ 端所缺的那一块。**
11. **抽象凹性路线确认关闭**：λ₁ 对算子恒**凸**（Rayleigh sup），凹性只能来自算子凹（Loewner），而 repo 探针已排除 ⟹ 论文可明写「no abstract route」，并引 Linear Algebra Appl. **169** (1992) 129 的凸性判据说明为何一般理论不给凹性。
12. **Hilbert / Carleman 型比较算子用于 Eliashberg 顶本征值界：文献零处**（KAY 用 ζ 型 Schur/HS；Elezaby–Abanov 用 Gershgorin）。空白地，但衰减阶不匹配（1/(k²+ϖ²) vs 1/(n+m+1)），须先做逐元比较验证再投入。
