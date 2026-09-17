# 定理 E（Einstein 极值性）对抗性查新 — 2026-09-16

**任务**：对「固定 λ 与 ⟨ω²⟩ 时 Einstein 谱极大化 T_c」做敌意查新。
**检索员**：researcher subagent。**未提交 git，未改动其他文件。**

**证据等级约定**
- **L1** = 本次会话逐页读了一手 PDF/HTML 页图，逐字转录
- **L1′** = 前序会话（2026-09-11 / 09-14）逐页读过一手源并记入我的持久记忆，本次未重读
- **L2** = 仅摘要 / 检索引擎综述转述
- **[自算]** = 我自己的推导或计算，不是任何来源的陈述

---

## 0. 裁决总表

| # | 项目 | 裁决 |
|---|------|------|
| 主命题 | 「固定 λ、μ*、⟨ω²⟩ 时 Einstein 谱极大化 T_c」这个**陈述** | **ALREADY KNOWN（逐字，Allen–Dynes 1975）** |
| 主命题 | 该陈述的**证明** | **NOT FOUND** — AD 只证了逐元引理，谱/T_c 那一步是无证断言 |
| 1 | Bergmann–Rainer 1973 是否蕴含或抵触极值性 | **既不蕴含也不抵触**（是不同的单调性） |
| 2a | Leavens / Carbotte 的「optimum spectrum」 | **ALREADY KNOWN 但约束不同**（固定**面积** A，不是固定 λ,⟨ω²⟩）；**且其证明同样只是局部变分论证**（§10.1，L1） |
| 2b | T_c ≤ 0.18√(λ⟨ω²⟩) 是否已证 | **NOT PROVEN** — 它正是 KAY II 的**公开 Conjecture 2** |
| 3 | 严格数学侧（KAY I–III / Frank–Hainzl–Seiringer 线） | **NOT FOUND** — 无人证，KAY 连猜想都没以此形式写 |
| 4 | 「高阶矩无关」的**定量**界（固定 ⟨ω²⟩ 与 ⟨ω⁴⟩ 后形状只买 ≤0.9%） | **NOT FOUND** |
| 5 | 反例文献（有谱在固定 λ,⟨ω²⟩ 下超过 Einstein） | **NOT FOUND**（并已排除一条看似成立的反例，见 §5.2） |

**一句话结论**：这**不是**「人人相信没人证」的民间传说，也**不是**已被证明——它是 **Allen & Dynes 1975 逐字写下、并把关键一步明确外包给「显然」的断言**。论文的新颖性只能落在「第一个证明」上，绝不能落在「第一个陈述」上。

---

## 1. 主命题：Allen & Dynes 1975 逐字已有（L1）

一手源：P. B. Allen and R. C. Dynes, *Transition temperature of strong-coupled superconductors reanalyzed*, Phys. Rev. B **12**, 905 (1975)。
取件路线：`harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.12.905/fulltext` → WebFetch 落盘 → Read(pages=) 读页图。本次逐页核对 p.908, 909, 911, 912, 918–921。

### 1.1 命题陈述（p.908 末 → p.909 左栏，L1 逐字）

> "As a final mathematical observation, it is interesting to ask what is the optimum shape for the function α²F **in the sense that T_c is maximized under the conditions of fixed λ, μ\*, and ⟨ω²⟩**. The asymptotic formula (22) indicates that at very large λ, the shape doesn't affect T_c. At smaller values of λ the shape does play a role. **It is clear from inspection of Eq. (15) that the eigenvalue ρ (and thus T_c) is maximized by making the parameters λ(n) as large as possible.** If we write λ(n) in terms of the normalized spectral function g(ω),
>
>     λ(n) = λ ∫₀^∞ dω g(ω)ω² / (ω² + (2πnT)²),     (29)
>
> we can then ask which shape of g(ω) (for fixed second moment ⟨ω²⟩) maximizes λ(n)/λ. **It is proved in Appendix A that the Einstein spectrum optimizes the integral (29) and thus is the preferred phonon spectrum.**"

### 1.2 再次陈述（p.911，L1 逐字）

左栏：
> "In order to discuss the role of the shape of α²F it is necessary to decide what parameters should remain invariant. **The natural choice is to fix λ and ⟨ω²⟩** … We are thus asking to know how T_c varies as the shape of α²F varies with λ and ⟨ω²⟩ fixed. **It has already been shown in Sec. II that the maximum T_c under these conditions occurs for an Einstein spectrum.**"

右栏：
> "The curves λ(n)/λ versus n normalized to n₀ = ⟨ω²⟩^{1/2}/2πT_c are shown in Fig. 6. **On this scale the Einstein spectrum is an upper bound for all n. This is proved in Appendix A and is the basis for the idea that the Einstein spectrum is optimal.**"

Fig. 6 图注（L1 逐字）：
> "In Appendix A it is shown that for an arbitrary phonon spectrum, λ(n)/λ never exceeds the value it takes for an Einstein spectrum."

### 1.3 Appendix A 真正证了什么（p.920–921，L1 逐字）

AD 定义 g(ω) ≡ (2/λω)α²F(ω)，⟨ωⁿ⟩ 为其矩，ω̄_n = ⟨ωⁿ⟩^{1/n}。要证的三条性质是 (A2)(A3)(A4)：

> (A4)   f(a) = ∫₀^∞ dω g(ω) ω²/(ω²+a²) ≤ f_E(a) = ⟨ω²⟩/(⟨ω²⟩+a²)

证法（p.921）：
> "The integrating factor x/(x+a²) has a negative second derivative for all positive values of x. **Thus (A4) is immediately proved using theorem (A10).**"

(A10) 是 **Wheeler–Gordon**（AD ref. 42）的矩逼近定理，前提逐字：
> "Theorem (A10) yields rigorous inequalities **provided** F(x) has a second derivative of fixed sign in the range of integration."

然后是 AD 把结论推广到一般矩的那段（p.921 右栏，**L1 逐字，这是最要命的一段**）：

> "be proved in a more general form, namely,
>
>     ∫dω g(ω) ω²/(ω²+a²) ≤ ⟨ωⁿ⟩^{2/n} / (⟨ωⁿ⟩^{2/n}+a²) ,     (A13)
>
> valid for n ≥ 2. **All of the inequalities mentioned here become exact equalities *only* if g(ω) is a δ function.** Thus the inequality (A13) would lead to a more general version of the statement made in the text; namely, **given λ, μ\*, and ω̄_n for an arbitrary value n ≥ 2, the maximum T_c consistent with these parameters is obtained when (and only when) α²F is an Einstein spectrum.** For simplicity we will give the proof only for n = 2 on the grounds that ⟨ω²⟩ is the most 'physical' moment of α²F."

### 1.4 精确的缺口（这是论文唯一能主张的东西）

AD 的逻辑链是两步：

| 步 | 内容 | AD 的状态 |
|----|------|-----------|
| (i) | 逐元：任意 P 的 λ(n)/λ ≤ Einstein 的，等号 iff δ | **已证**（Appendix A，(A4)，经 Wheeler–Gordon (A10)） |
| (ii) | 逐元 ⟹ 顶本征值 ρ ⟹ T_c 极大 | **无证断言**："clear from inspection of Eq. (15)"；"the basis for **the idea** that the Einstein spectrum is optimal" |

(ii) **不是平凡的**：λ(n) 同时进入配对核（增强）与 Z 质量重整化对角项（压制），两个方向相反，所以「λ(n) 变大 ⟹ ρ 变大」需要真正的证明，而不是「inspection」。[自算：这是我对为什么 (ii) 非平凡的论证，AD 没有讨论这一点。]

**另外三条必须在论文里说清的边界**：

1. **AD 的断言比定理 E 更宽**：AD 的 (A13) 段把结论写成「对任意 n ≥ 2 的 ω̄_n」，且**包含 μ\* ≠ 0**。定理 E 只覆盖 n=2 且 μ\*=0。也就是说，我们证的是 AD 所**声称**范围的一个特例。这必须显式承认，否则会被审稿人抓。
2. **AD 自己承认 μ\*≠0 时严格性丢失**（p.908 右栏，L1 逐字）："However, the convergence is no longer monotonic, and inequalities of the type (18) and (23) are no longer rigorously true." —— 这反而支持定理 E 把 μ\*=0 作为假设是**合理且与 AD 一致**的。
3. **等号条件 AD 已经写了**："become exact equalities *only* if g(ω) is a δ function"。所以「equality iff point mass」也不是新陈述（在逐元层面已证，在 T_c 层面仍是断言）。

### 1.5 二手确认：Moussa & Cohen 2006（L1 逐字，ar5iv）

R. Moussa & M. L. Cohen, *Two bounds on the maximum phonon-mediated superconducting transition temperature*, PRB **74**, 094520 (2006) = cond-mat/0607832。

> "We attempt to tighten the bounds by independently constraining λ and ⟨ω²⟩ as well as specifying μ\*. **It is already known [allen_dynes] that the Einstein spectrum maximizes T_c for a given λ and ⟨ω²⟩.** However, the lower bound on T_c is still zero even with these additional constraints."

其中 `[allen_dynes]` = AD 1975。同文另有：
> "0 ≤ T_c ≤ 0.183 √(λ⟨ω²⟩)"
> "Both bounds can be attained by an Einstein model phonon spectrum, α²F(ω)_E = λω_E δ(ω−ω_E)/2, the lower bound corresponding to ω_E → ∞ and the upper bound corresponding to ω_E → 0"
> "The spectra associated with the two bounds are both unphysical pathologies, corresponding to the cases of maximum phonon softening and infinite phonon hardening."

**关键**：MC **没有**独立核验这条，他们把它当既成事实引用。他们自己做的数值是**最小化**谱（固定 λ, ⟨ω²⟩, ω_max），脚注 [12] 逐字：
> "The dT_c/dα²F(ω) response for fixed λ and ⟨ω²⟩ is numerically calculated to be positive in the frequency interval from zero to ω_max."

（注：这条**受约束的**泛函导数正是本问题的正确对象，MC 只是数值算了它，没有证明。）

---

## 2. Bergmann & Rainer 1973：证的是别的单调性（L1）

一手源：G. Bergmann and D. Rainer, *The sensitivity of the transition temperature to changes in α²F(ω)*, Z. Physik **263**, 59–68 (1973)。
取件：`physics.umd.edu/courses/Phys798C/AnlageSpring22/Bergmann-1973-The%20sensitivity%20of%20the%20transitio.pdf`（UMD 课程页扫描件）→ Read(pages=)。

摘要（p.59，L1 逐字）：
> "It is shown that **all frequency regions of α²F(ω) yield a positive contribution to T_c** and that the most effective range covers frequencies, **slightly above 2πT_c**."

p.60（L1 逐字）：
> "In all examples considered here the functional derivative is always positive. **Furthermore we can prove this property rigorously for μ\*=0 (see Appendix A).** Therefore we are convinced that in thermal equilibrium phonons never have a pair-breaking effect on the transition temperature."

p.62（L1 逐字）：
> "In Appendix A we prove rigorously (for μ\*=0) that the functional derivative of ρ has always a positive sign. This implies that the functional derivative of T_c which differs only by the positive factor (−dρ/dT)^{−1}, is also a positive function."

**BR 证的是**：δT_c/δα²F(ω) > 0 对所有 ω，μ\*=0（严格，Appendix A）。即**加谱重必涨 T_c**。

**BR 没证的是**：形状极值性。定理 E 比较的两个谱有**相同**的 λ 与 ⟨ω²⟩，因此互相之间是**不定号扰动**，BR 的正定性对它一句话都说不上。[自算]

**对 team-lead 问的「峰在 ω≈7T_c 是否蕴含或抵触极值性」**：
- **数值订正**：BR 原文写的是 **"slightly above 2πT_c"**，不是 7T_c。2π ≈ 6.28，坊间的「7k_BT_c」是这个数的口口相传版本。引用时用原话。（Marsiglio–Carbotte 综述里又写成 "a factor of order 10 times the critical temperature"——三种措辞并存，只引原话。[L1′]）
- **既不蕴含也不抵触**。峰位置说的是「往哪儿加重最划算」，这是**固定面积**优化问题的答案（→ §3 Leavens/Carbotte）；极值性说的是**固定 λ 与 ⟨ω²⟩**下的形状比较，两个约束流形不同。BR 的结果与 Leavens–Carbotte 的 δ-最优是同一条线，与定理 E 不是。[自算]

---

## 3. Leavens / Carbotte / 「0.18√(λ⟨ω²⟩)」（L1′ + L2）

### 3.1 四类互不等价的「最优谱」问题——混用即误判

| 固定什么 | 结论 | 来源 | 状态 |
|---|---|---|---|
| **面积** A = ∫α²F dω | δ 峰在 ω\*(μ\*) 处最优；T_c ≤ A·c(μ\*) | Leavens SSC **17**, 1499 (1975)；Leavens–Carbotte JLTP **14**, 195 (1974)；Mitrović–Carbotte SSC **40**, 249 (1981)；Carbotte RMP **62**, 1027 (1990) | 变分/线性响应论证，**非严格** [L2] |
| **λ 与 ⟨ω²⟩** | Einstein 极大 | **AD 1975**（本报告 §1） | 逐元已证，谱步断言 [L1] |
| **ω_log** | 谱宽↑ ⟹ T_c 一致↓ | Combescot PRB **42**, 7810 (1990) eq.(56)，R ≥ 0 且 R=0 iff Einstein | 弱耦合展开级 [L1′] |
| **ω_max**（或固定 λ 求最大 ω_log） | δ(ω−ω_max) 最优 | Semenok–Altshuler–Yuzbashyan arXiv:2407.12922 SI §S4；Gao et al. arXiv:2502.18281 = Nat. Commun. **16**, 8253 (2025) | 变分导数，**非严格** [L1] |

**Carbotte RMP 62, 1027 (1990) 的 optimum spectrum 是固定面积那一类**，不是我们的约束。[**已升格 L1，见 §10** —— RMP 正文已逐页取到并核实，§II.H 逐字确认约束是固定 A，且其证明只是**局部**极大论证。]

### 3.2 「T_c ≤ 0.18√(λ⟨ω²⟩)」的真实身份（L1 + L1′）—— 三个数别混

AD 1975 原文（L1，p.908 右栏）：
- **eq.(22)**：`T_c → (λ⟨ω²⟩)^{1/2} f(μ*(N))`，**无数值常数**，"f is a function of μ\* only"。
- **eq.(24)**：`T_c ≥ T_c1 = 0.180(λ⟨ω²⟩)^{1/2}  (λ>10, μ*=0)` —— 这是 2×2 截断给的**下界**，方向是 ≥。
- **eq.(25)**：`T_c ≃ 0.182(λ⟨ω²⟩)^{1/2}  (λ≳10, μ*=0)` —— N+1=64 的**数值**，"exact behavior in this limit, as deduced from computer calculations"。
- AD 全文**没有** 0.1827 或 0.183。

⟹ **AD 从未证明过任何形如 0.18√(λ⟨ω²⟩) 的上界**。

**订正（§10.2，本次新取 Carbotte RMP 正文后）**：数字 **0.183 不是 Moussa–Cohen 首创**，它印在 Carbotte RMP 62, 1027 (1990) **eq.(2.29)**，且在那里是 Einstein 谱在 λ→∞ 的**等式**（不是对一般 P 的上界）。MC 的贡献是把这个**等式**配上 AD 的极值性**断言**，翻成对一般 P 的**上界**——因此 MC 的 0.183 上界继承了 §1.4 的缺口。**方向别混：Carbotte 是「=」，MC 是「≤」。**

**而真正的 0.1827 是 KAY 的常数**：C_∞ = (1/2π)√𝔤(2) = 0.1827262477…（KAY III p.28 印全十位）。[L1′]

**⟹ 这个上界目前是公开猜想，不是定理。** KAY II（arXiv:2409.00532 = J. Stat. Phys. 192, art. 94 (2025)）**Conjecture 2**：存在 T_c(λ,P) > 0，对**所有** P∈𝒫 与 λ>0 满足
`T_c(λ,P) < T̃_c := (1/2π)√(𝔤(2)⟨ω²⟩λ)`（严格不等号，含存在性断言）。[L1′ 页图核过；本次 HTML 复取在实质上一致（存在性 + 严格 <），仅页码标注 12 vs 14 不同——HTML/PDF 分页差异。]

KAY II 已证的上界常数是 **2.034 × C_∞**（约 11.4 倍松于猜想值的 …差 2.034 倍）。[L1′]

**⟹ 这是本次调研发现的最强卖点**：定理 E（μ\*=0）配上 T_c^E(λ,Ω)/(√λ Ω) 对 λ 单调趋于 C_∞，直接给出 KAY II Conjecture 2。**一个 2025 年发表的公开猜想被解决**，这个 framing 比「Einstein 极值性」强得多，且不与 AD 撞车。
（**⚠ 这一步是我的推断**：需要仓库侧独立确认 T_c^E/√λΩ 确实对 λ 单调上升到 C_∞——若只是渐近而非单调，蕴含关系不成立。记忆里仓库已记 "Conj. 2 COMPLETE (μ\*=0) + Theorem E"，但我没有核验那份证明。）

---

## 4. 严格数学侧：无人证，且 KAY 明确点名 AD 的无证断言（L1）

### 4.1 KAY I 的措辞（arXiv:2409.00533，§1.3，L1 逐字）

> "However, when the first N > 1 positive Matsubara frequencies are involved in the truncation the largest eigenvalue of the pertinent N×N matrix does not manifestly depend suitably monotonically on T, **claims to the contrary in [AD] notwithstanding**."

这**不是**极值性那条断言（它针对的是截断本征值对 T 的单调性），但它确立了一个极其有用的事实：**KAY 小组逐条审计过 AD 的无证断言，并公开点名了其中一条**——却**没有**采纳、也没有证明、甚至没有以猜想形式写下极值性那条。

### 4.2 KAY II 是否猜了极值性？否（L1，arXiv HTML）

本次对 2409.00532 的全文检索：Einstein 模型只在 eq.(3) 作为特例定义（`2α²F(ω)|_E := λΩδ(ω−Ω)`），并注明「will be evaluated in [KAYc]」；对 AD 的唯一实质引用是 **Remark 3**：
> "The asymptotic behavior T_c(λ,P) ∼ (1/2π)√(𝔤(2)⟨ω²⟩λ) was first stated in [AD], see their eq.(22)."
（**注意**：这句把 eq.(22) 说成 "asymptotic behavior … first stated"，是对的——AD eq.(22) 只有形式没有常数，见 §3.2。）

**没有任何猜想或定理声称点质量 P 在固定 λ 与 ⟨ω²⟩ 下极大化 T_c。**
⚠ 这次 WebFetch 对 Conjecture 1 的转述明显是 mis-parse（返回了一段关于 spin chain 的乱文），因此这条**否定结论我按 L1.5 记**：正面引文（Remark 3、eq.(3)）与我的记忆一致可信，但「全文无极值性猜想」这个**否定**结论建议 lead 亲自 grep 一遍 PDF 再落笔。

### 4.3 Frank–Hainzl–Seiringer 线（L2）

这条线是**瞬时核 BCS**（K_T(p) 算子、外场依赖、弱耦合/低密度渐近），**不覆盖推迟核 Eliashberg**，也不做 α²F 形状优化。与定理 E 无撞车。[与我记忆中「retarded-kernel rigorous T_c bounds are virgin territory」一致]

### 4.4 Elezaby & Abanov（L1，摘要）

arXiv:2512.20009v2（v1 2025-12-23，v2 2026-06-10），*Superconductivity Near a Quantum Critical Point: Bounds on the Transition Temperature in the γ-Model*。摘要逐字核过：全篇是 **γ-model** 的上下界（Hessian 线性代数 + Rayleigh–Ritz），**摘要零处提 Einstein 谱最优、零处引 Allen–Dynes**。不撞车。
（⚠ 提醒：Elezaby–Abanov 的单调阶梯 τ_c^(1) < τ_c^(2) < … 是**变分截断维数**的阶梯，与谱形状无关——上一轮已记录过这个易混点。）

### 4.5 Esterlis–Kivelson–Scalapino 那一线（L1 摘要）

arXiv:1806.00488 = npj Quantum Materials **3**, 59 (2018)。摘要逐字：
> "we propose that k_B T_c ≤ A_max ħω̄ where A_max is a dimensionless number of order one that we estimate to be A_max ≈ 1/10."

这是**基于 Holstein 模型 DQMC + 经验数据的近似界**（λ 有上限 ⟹ T_c/ω̄ 有上限），**不是** α²F 形状泛函的界，摘要零处提 α²F 形状、矩、Einstein 极值。后续 Hofmann–Chowdhury–Kivelson–Berg，npj Quantum Mater. **7**, 83 (2022) 明确论证这类界都不是基本界。**与定理 E 不同类，不撞车，但审稿人可能要求在引言里区分开。**

Sadovskii arXiv:2407.03602（L1 摘要逐字）：
> "We present simple qualitative estimates for the maximal superconducting transition temperature … in the limit of very strong coupling the upper limit for transition temperature is determined in fact by a combination of atomic constants and density of conduction electrons."
自称 "qualitative estimates"，不涉形状极值性。不撞车。

---

## 5. 「形状只买多少」的定量陈述 & 反例

### 5.1 定量的矩压缩界：NOT FOUND

检索了「固定 ⟨ω²⟩（甚至再固定 ⟨ω⁴⟩）后形状还能改变 T_c 多少」的**双边定量**陈述。**没找到任何来源给出这种界。** 文献里有的都是别的东西：

- **AD 1975 自己给的是反向信号**（p.912，L1 逐字）：
  > "There is a **surprisingly large sensitivity of T_c/ω̄₂ to the shape of the spectrum**, especially considering that for λ ≳ 10, this sensitivity must disappear."
  其 Fig. 8（T_c/⟨ω²⟩^{1/2} vs λ，μ\*=0.1，三个形状 Einstein / Pb / Hg）本质上**就是 η 的图**。目视读图 λ≈2 处 Einstein≈0.115、Pb≈0.105、Hg≈0.095 ⟹ η(Hg) ≈ 0.83。[自算读图，误差大，仅作量级参照]
  ⟹ **「η ≥ 0.86 对 LaH10/H3S/Nb/MgB2」不是新类型的观察**，是 AD Fig. 8 的现代材料版；但 AD 只有三个形状、只有 μ\*=0.1、且明说这是「surprisingly large」。
- **Xie et al., npj Comput. Mater. 8, 14 (2022) = arXiv:2106.05235**（L1，读了 p.1–4）量化的是 **AD 拟合公式 vs 数值 Eliashberg 的误差**，不是 Einstein 与一般谱的差。逐字：
  > p.1: "while the Allen-Dynes formula accurately predicts the Eliashberg T_c for λ values near 1.6 …, it nevertheless deviates from the Eliashberg T_c when λ is significantly larger or smaller than 1.6 **and when the shape of α²F(ω) differs from the simple unimodal Einstein model**."
  > p.4: "ratio ⟨ω_log⟩/⟨ω̄₂⟩ is 1, the shape of α²F is that of the unimodal Einstein model and the Allen-Dynes T_c accurately predicts the Eliashberg T_c regardless of the coupling strength."
  > p.3: "The root mean square error ratio in the Allen Dynes paper is around 5.6%"
  其 Fig. 4 显示 T_c(AD)/T_c(Eliashberg) 在双峰谱、ω_log/ω̄₂ ∈ [0.45, 1.0] 区间内跨越 **0.65 到 1.2**。
  ⟹ **这对论文的 corollary framing 是直接威胁**（见 §7）。
- **Combescot 1990 eq.(56)**：固定 ω_log 时形状因子 R ≥ 0，R=0 iff Einstein——这是**弱耦合展开级**的定量形状项，但约束是 ω_log 不是 ⟨ω²⟩。[L1′]

**结论**：定理 E 的「固定 ⟨ω²⟩ 与 ⟨ω⁴⟩ 后形状只买 ≤0.9%」若成立，**在文献中确实没有对应物**。这是本项目最干净的增量。工具上有讽刺意味：AD Appendix A 用的 Wheeler–Gordon 矩逼近定理正是做「固定前几阶矩的两点求积界」的标准机器——AD 把它用在 n=2 单矩上就停了。

### 5.2 反例：NOT FOUND，并驳回一条貌似成立的反例

未找到任何报告称某个谱在固定 λ 与 ⟨ω²⟩ 下 T_c 超过 Einstein。多峰/ML 优化文献（Xie 2106.05235；BETE-NET arXiv:2401.16611；Gao arXiv:2502.18281）做的都是**预测**或**固定 λ / 固定 ω_max** 下的优化，没人在我们的约束流形上找反例。

**但检索引擎给出了一条看起来很像反例的论证，必须当场驳回**（审稿人很可能提同一条）：

> 「AD 的形状因子 f₂ = 1 + λ²(ω̄₂/ω_log − 1)/[λ² + 1.82(1+6.3μ\*)(ω̄₂/ω_log)²] 在 ω̄₂ > ω_log 时 > 1，所以 AD 公式自己就预言多峰谱能在固定 λ、⟨ω²⟩ 下超过 Einstein。」

**这是错的。** [自算，本次从零推导]：AD 公式里 T_c ∝ f₁(λ,μ\*)·f₂·ω_log，而 f₁ 不含形状。固定 λ, μ\*, ω̄₂，令 r ≡ ω̄₂/ω_log ≥ 1，则 ω_log = ω̄₂/r，且记 c ≡ 1.82(1+6.3μ\*) > 0：

    T_c ∝ (ω̄₂/r)·[1 + λ²(r−1)/(λ²+cr²)]
        = ω̄₂ · (λ² + c r) / (λ² + c r²)  ≡ ω̄₂ · F(r)

F(1) = 1，且

    F′(r) = c[λ²(1−2r) − c r²] / (λ²+cr²)²

对 r ≥ 1 恒有 (1−2r) ≤ −1 < 0，故分子 < 0 ⟹ **F 在 r ≥ 1 上严格递减 ⟹ F(r) < 1 对所有 r > 1**。
⟹ **AD 拟合公式本身也给出 Einstein（r=1）是固定 λ、μ\*、⟨ω²⟩ 下的严格极大值。** 那条「反例」漏掉了 ω_log = ω̄₂/r 这个 1/r 因子。

**副产品（应写进论文）**：F 是 λ 的函数，λ→0 时 F(r)→1 对所有 r —— 即 **AD 公式在弱耦合端对形状完全失明**。而 Combescot 1990 证明弱耦合端形状确实压制 T_c（固定 ω_log）。二者不矛盾（f₂ 是强耦合修正因子，按构造 λ→0 时 →1），但说明 **AD 的「两参数规整性」在弱耦合端是平凡的，不是有内容的**。

---

## 6. 最相关的五篇

1. **Allen & Dynes, PRB 12, 905 (1975)** — 逐字陈述了完全相同的命题（固定 λ, μ\*, ⟨ω²⟩），Appendix A 证了逐元不等式 (A4)/(A13)（经 Wheeler–Gordon），谱/T_c 那一步是无证断言。**必引，必逐字引，必在引言前半页说明缺口。**
2. **Moussa & Cohen, PRB 74, 094520 (2006) = cond-mat/0607832** — 逐字 "It is already known [allen_dynes] that the Einstein spectrum maximizes T_c for a given λ and ⟨ω²⟩"，并由此得 0 ≤ T_c ≤ 0.183√(λ⟨ω²⟩)。**证明这条在 2006 年已是社区常识，不能主张陈述新颖性。**
3. **Kiessling–Altshuler–Yuzbashyan I & II, arXiv:2409.00533 / 2409.00532, J. Stat. Phys. 192 (2025)** — KAY II **Conjecture 2** 就是尖锐上界 T_c < (1/2π)√(𝔤(2)⟨ω²⟩λ)，公开未决；KAY I §1.3 "claims to the contrary in [AD] notwithstanding" 是 AD 断言不可采信的公开背书。**这是最有力的立论支点。**
4. **Bergmann & Rainer, Z. Physik 263, 59 (1973)** — 唯一一条严格证明的相关单调性（δT_c/δα²F > 0，μ\*=0，Appendix A），峰在 "slightly above 2πT_c"。**方法论上是「对 α²F 的严格结论是可以证的」的先例，但不覆盖定理 E。**
5. **Xie et al., npj Comput. Mater. 8, 14 (2022) = arXiv:2106.05235** — 现代大规模数值，量化 AD 公式在形状偏离单峰 Einstein 时误差可达 −35%/+20%。**是 corollary framing 的直接威胁，必须在论文里主动区分「精确上界定理」与「AD 拟合公式精度」。**

---

## 7. 对 framing 的具体建议（含一条必须改的）

**可以这样写**：
- 「AD 1975 陈述了这一极值性质并证明了其逐元前提；从逐元到本征值/T_c 的推断在原文中是无证断言（"clear from inspection"）。我们给出第一个完整证明，对全部耦合强度、任意谱测度、μ\*=0，无计算机辅助。」
- 「作为推论，我们解决了 Kiessling–Altshuler–Yuzbashyan 的 Conjecture 2（μ\*=0 情形）。」← **这是最强、最不撞车的卖点，建议提到摘要**。

**必须修掉的一条**：
> 「α²F compresses to (λ, ⟨ω²⟩) with a provably bounded residual — i.e. the Allen–Dynes two-parameter regularity is a theorem, not an empirical fit.」

三个问题：
1. 定理 E 是**单边上界**，不给双边压缩；「bounded residual」需要下界，定理 E 没有（Moussa–Cohen 逐字说固定 λ 与 ⟨ω²⟩ 后**下界仍是 0**）。
2. **AD 本人写的是反话**："surprisingly large sensitivity of T_c/ω̄₂ to the shape"。直接被原文打脸。
3. 「Allen–Dynes 两参数规整性」指的是 AD 的**拟合公式**（f₁,f₂），而 Xie et al. 2022 已发表证据显示它在形状变化下误差 −35%…+20%。把定理 E 说成「AD 拟合是定理」会被这篇直接反驳。

**建议改为**：「定理 E 表明 (λ, ⟨ω²⟩) 给出 T_c 的**尖锐上界**；再固定 ⟨ω⁴⟩ 后残余形状自由度在我们考察的谱族上只改变 T_c ≤0.9%（**numerical-evidence**）。这与 Allen–Dynes 拟合公式的精度是两回事，后者独立地被 Xie et al. 2022 测出可达数十个百分点的误差。」

---

## 8. 「读到的」vs「推断的」清单

**读到的（L1，本次会话逐页核过）**
- AD 1975 p.908 / 909 / 911 / 912 / 918–921 的全部逐字引文，含 (A4)(A13)(22)(24)(25)(29)(30)、Fig. 6/7/8 图注
- Moussa–Cohen cond-mat/0607832 的 "already known" 段、0.183 界、两端点谱、脚注 [12]
- Bergmann–Rainer 1973 摘要、p.60、p.62 的逐字段
- KAY I arXiv:2409.00533 §1.3 的 "notwithstanding" 整句；KAY II Remark 3 与 eq.(3)
- Xie et al. arXiv:2106.05235 p.1–4（含 eq.(2)(3)(4)(5)(6)(7)、Fig.4、Table I、5.6% RMS 句）
- EKS arXiv:1806.00488 摘要；Sadovskii arXiv:2407.03602 摘要；Elezaby–Abanov arXiv:2512.20009 摘要
- Semenok–Altshuler–Yuzbashyan arXiv:2407.12922 SI §S4 的 δ-最优句与 eq.(7)(8)(9)

**L1′（前序会话一手核过，本次未重读）**
- KAY II Conjecture 2 的 eq.(31) 精确形式与「严格 <、含存在性」
- KAY III p.28 的 C = 0.1827262477
- Combescot PRB 42, 7810 (1990) eq.(21)/(56) 与 R ≥ 0 iff Einstein
- Marsiglio–Carbotte cond-mat/0106143 p.33 eqs.(98)(99) 的 Leavens 归属

**L2（只看到摘要/检索综述，未取正文）**
- Carbotte RMP 62, 1027 (1990) 的 optimum-spectrum 小节内容（转述自 cond-mat/0101419 与 cond-mat/0106143）
- Leavens SSC 17, 1499 (1975) 与 Leavens–Carbotte JLTP 14, 195 (1974) 的具体陈述
- Gao et al. arXiv:2502.18281 / Nat. Commun. 16, 8253 (2025) 的「最优 α²F 不物理」论证细节
- Hofmann–Chowdhury–Kivelson–Berg npj QM 7, 83 (2022)

**我的推断（不是任何来源的陈述）**
- 「(ii) 非平凡，因为 λ(n) 同时进配对核与 Z 项」——AD 没讨论过
- 「BR 的正定性对定理 E 完全用不上，因为定理 E 的比较是不定号扰动」
- 「定理 E + T_c^E/√λΩ 对 λ 单调 ⟹ KAY Conjecture 2」——**单调性我没核验**
- §5.2 的 F(r) = (λ²+cr)/(λ²+cr²) 严格递减推导（本次从零算，驳回检索引擎提出的「反例」）
- AD Fig. 8 目视读图得 η(Hg) ≈ 0.83
- 「AD 的 f₂ 在 λ→0 时退化为 1，所以 AD 的两参数规整性在弱耦合端无内容」

---

## 9. 检索局限（诚实边界）

1. ~~**Carbotte RMP 62, 1027 (1990) 正文未取到**~~ —— **已闭合，见 §10**。RMP 全文（PDF 7.4MB，152 页）已取到并逐页核实 p.1039–1047, 1095–1096, 1120–1121。裁决未变（固定面积 A），并额外拿到三项对论文有用的一手材料。
2. **Leavens SSC 17, 1499 (1975) 与 Leavens–Carbotte JLTP 14, 195 (1974) 原文未取到**（1975 年 SSC 无 arXiv、无免费扫描）。
3. **KAY II「全文无极值性猜想」是否定结论**，本次 WebFetch 对 Conjecture 1 的转述已证明不可靠（返回乱文），建议 lead 亲自 grep PDF。
4. **Moussa–Cohen 只取了 ar5iv**，未核对印刷版页码/式号。
5. 未逐篇检索 1975–1990 之间可能存在的、以俄文或会议论文形式发表的极值性证明（AD eq.(25) 脚注 14 提到 Rainer 1973 Gatlinburg 会议 unpublished 与 Garland 私人通信有等价推导——**这两个来源原则上不可取**）。
6. WebSearch 的布尔/引号支持差，部分查询靠关键词组合；检索引擎给出的综述含至少一处实质性错误（§5.2 的伪反例），已当场驳回——**本报告中任何未标 L1/L1′ 的内容都不要直接落笔到论文。**

---

## 10. 补充轮：Carbotte RMP 62, 1027 (1990) 正文已取到（L1，最大缺口已闭合）

取件路线：`harvest.aps.org/v2/journals/articles/10.1103/RevModPhys.62.1027/fulltext` → WebFetch 落盘（PDF 7.4MB）→ Read(pages=) 读页图。**PDF 页码 = RMP 页码 − 1026。**
本次逐页核实：p.1027–1028（目录）、1039–1047（§II.B–J）、1095–1096（§VI.D–F）、1120–1121（§VIII.A–B）。

### 10.1 §II.H「The optimum spectrum for T_c」逐字（p.1045）—— §3 裁决**确认不变**，且发现同型缺口

> "We can prove that (2.48) represents a **local maximum** for any shape of α²F(ω) **with a fixed A**. In the lower part of Fig. 5, we show the functional derivative of T_c/A with respect to α²F(ω) for the optimum base delta-function spectrum
>
>     α²F(ω) = A δ(ω − ω\*_E(μ\*))     (2.49)
>
> with ω\*_E(μ\*) = 1.3 A for μ\* = 0.1. In this case A δ(T_c/A)/δα²F(ω) is negative definite and zero exactly at ω̄_E = 1.3 ≡ ω\*_E(μ\*). That is, the functional derivative is exactly zero at the frequency of the base delta function (2.49) and negative everywhere else. This means that adding weight to the spectral density at this frequency leaves T_c/A unchanged, and adding weight anywhere else reduces it. Thus we have maximized T_c/A with the optimum spectrum (2.49) for arbitrary A. It follows that the inequality
>
>     T_c ≤ A c(μ\*)     (2.50)
>
> holds. The equality applies to a delta function and the inequality to any other spectral shape. **This result was first established by Leavens (1975) using somewhat different arguments.**"

**三条结论**：

1. **约束确认是固定面积 A**（"with a fixed A"、"for a given spectral weight A"），**全篇没有一处用固定 λ 与 ⟨ω²⟩**。§3.1 表格的分类成立，现由 L2 升格 **L1**。
2. **Leavens 1975 的归属现有一手出处**（此前只有 Marsiglio–Carbotte 的转述）。
3. **Carbotte 的证明与 AD 有同型缺口**：他自陈只证了 "a **local** maximum"（一阶变分 + 泛函导数负定），然后把全局不等式 (2.50) 当作 "It follows that" 写下。**局部驻点 ⟹ 全局极大**这一步同样没有证明。

而且 Carbotte 在类比小节里**明确区分了局部与绝对**（p.1096，L1 逐字，说的是比热跳变而非 T_c，但方法完全相同）：
> "We can use functional derivatives to prove that any other shape will lower ΔC(T_c)/(γT_c). **Our experience with realistic spectra would lead us to believe that this is, in fact, an absolute maximum.** … They all fall below our theoretical local maximum, **indicating that** for physical systems, the solid line is, indeed, an absolute maximum."

⟹ **这是本轮对论文最有用的一条**：δ-谱最优性这一整类命题，在其**最权威的综述**里就是「局部已证 + 全局靠经验/数值」。这**不是**我们对 AD 的单点挑刺，而是该类结论 1974–1990 年的**共同状态**。可以据此写：「δ-函数极值性在文献中从 Leavens (1975) 到 Carbotte (1990) 一直是局部变分论证；全局命题（在任一约束流形上）此前没有证明。」

**补记**：§VI.D（比热跳变，p.1095）与 §VI.F（H_c(0)，p.1096）的最优谱结论同样是 "for a fixed spectral weight A" / "For a fixed spectral density A"。**RMP 全篇的最优谱机器都是固定面积的。**

### 10.2 §II.E eq.(2.29)：0.183 的真实出处（p.1042，L1 逐字）

> "T̄_c = 0.2584 or T_c = 0.258√(ω_E A) = **0.183√λ ω_E** , (2.29)
> a result first obtained by Allen and Dynes [(1975); **Owen and Scalapino (1971)**] which holds for μ\* = 0. We shall see later that a single-Matsubara-gap approximation would have given 1/(π√2) for the coefficient in (2.29) rather than 0.258."

**四点**：

1. **数字 0.183 印在 Carbotte 1990，不是 Moussa–Cohen 首创**（§3.2 已订正）。
2. **这里是等式不是上界**：它是 **Einstein 谱自身**在 λ→∞ 的渐近值。要变成对一般 P 的上界，必须外加 AD 的极值性断言——这正是 MC 干的事，也正是缺口继承的地方。
3. **[自算] Carbotte 的 0.2584 就是 KAY 的 C_∞**：对 α²F = Aδ(ω−ω_E) 有 λ = 2A/ω_E ⟹ A = λω_E/2 ⟹ √(ω_E A) = ω_E√(λ/2)，故系数换算为 0.2584/√2 = 0.182716。Carbotte 只印 4 位，其舍入区间 0.2584±0.00005 映到 [0.182681, 0.182751]，**包含** KAY III 的 C_∞ = 0.1827262477。⟹ 这个常数在 1990 年已有 4 位数值，KAY 的贡献是给它闭式 (1/2π)√𝔤(2) 并放进严格框架。论文里提 C_∞ 时应引 Carbotte (2.29) 作为数值先例。
4. **新线索：Owen & Scalapino (1971)** 被 Carbotte 列为与 AD 1975 并列的渐近式来源，**早于 AD 四年**。
   **本次已定位完整著录（L1.5）**：C. S. Owen & D. J. Scalapino, *S-state instabilities for retarded interactions*, **Physica 55, 691–695 (1971)**, DOI `10.1016/0031-8914(71)90320-X`。
   出处：Scalapino 本人在 UCSB 的著作列表（`web.physics.ucsb.edu/~djs/djspublist.php`）；标题经 Semantic Scholar API 独立复核一致（`api.semanticscholar.org/graph/v1/paper/DOI:10.1016/0031-8914(71)90320-X`，abstract 字段为 null）。
   **正文未取到**：doi.org → linkinghub.elsevier.com 只回重定向桩；Physica（1973 年拆分前）无免费扫描、无 arXiv。
   **[自算] 风险评估：它含有「固定 λ 与 ⟨ω²⟩ 下 Einstein 谱极大」这一陈述的概率很低**，四条理由：(a) 用 ⟨ω²⟩ / ω_ln 作为 α²F 的组织参数是 AD 1975 引入的——Carbotte RMP p.1040 逐字写 "ω_ln, which was first used by Allen and Dynes (1975)"；(b) Carbotte 只在**渐近式** (2.29) 的方括号里引它，且主语句是 "first obtained by Allen and Dynes"；(c) Moussa–Cohen 2006 把极值性**单独**归于 [allen_dynes]；(d) 标题是「推迟相互作用下的 s 态不稳定性」，指向配对不稳定性的存在性，不是谱形状优化。
   **但这仍是一条未读的一手文献**，若要在论文里写「AD 1975 是该陈述的最早出处」，应当降一格写成「据我们所知最早的明确陈述」，或先设法取到这 5 页。

### 10.3 §II.D 标度定理（p.1041，L1 逐字）—— 对定理 E 有正面用处

> "We start by establishing a scaling theorem (Coombes and Carbotte, 1986a) for T_c that **applies exactly when μ\* = 0** and is approximately true when μ\* is not an important player in T_c … Consider a base spectrum α²₀F(ω) and another denoted by α²_MF(ω), which is constructed from the base spectrum through scaling of the frequency axis by a constant factor of b; i.e.,
>
>     α²_MF(ω) = α²₀F(bω)     (2.20)
>
> … Reference to Eq. (2.16) gives that ω^M_ln = ω⁰_ln/b and reference to Eq. (2.5) yields λ_M = λ₀; i.e., the mass-renormalization function λ is unchanged. … we immediately conclude that
>
>     T^M_c = T⁰_c/b .     (2.23)
>
> This relationship implies that the characteristic strong-coupling ratio T_c/ω_ln is **completely unaffected** by a uniform stretch or compression of the frequency scale for the electron-boson spectral density."

⟹ 定理 E 把 ⟨ω²⟩ 归一到 1 用的正是这条齐次性。**它已发表**（Coombes & Carbotte 1986a，Carbotte RMP (2.20)–(2.23)），论文里应当直接引用而不是重证；同时它**独立佐证 μ\*=0 是这类结论的自然假设**（Carbotte 逐字："applies exactly when μ\*=0"）。

### 10.4 §II.E / §II.I：Eliashberg 理论对 T_c **无**上界（p.1042, 1046，L1 逐字）

> p.1042: "**There is, in fact, no limit to T_c in Eliashberg theory.** … Leavens (1975, 1977) and Leavens and Carbotte (1974, 1977) arrived at the conclusion that T_c is not bounded in Eliashberg theory by using very different arguments from those discussed so far."
> p.1046: "Equation (2.48), like the asymptotic form (2.29), implies that **there is no limit placed on T_c by the theory of superconductivity itself.** … We conclude, as previously described, that the criteria for maximum T_c need to be found outside Eliashberg theory."

⟹ 对 framing 有用：任何 T_c 界**必须**写成归一化形式（T_c/√(λ⟨ω²⟩) 或 T_c/A），定理 E 与 KAY Conjecture 2 正是这种形式。可以用这两句把定理 E 与 §4.5 的 EKS/Sadovskii/Hofmann 那类「绝对 T_c 上限」彻底切开——**Carbotte 1990 已经说清那类界不在 Eliashberg 理论内部**。

### 10.5 「7T_c」这个数的出处（p.1043–1044，L1 逐字）—— §2 的订正现在完全解释清楚

> p.1043（两方阱模型）: "the curve exhibits a single broad maximum at **Ω̄\* ≅ 10.0** with maximum value G(Ω̄\*) = 0.8 … This leads to the idea of the existence of an optimum frequency for T_c."
> p.1043 eq.(2.34)（经典图像估计）: "ω_op = πv_F/2ξ₀ **~ 7 k_B T_c**"
> p.1044（精确数值）: "in particular, the maximum in the exact curves occurs **around 7.0 k_B T_c** rather than 10.0 k_B T_c."

⟹ 三个并存措辞的来源现在齐了：**BR 1973 原文 "slightly above 2πT_c"（≈6.3，μ\*=0，严格）；Carbotte 精确数值 ≈7.0 k_BT_c（含 μ\*）；Carbotte 两方阱模型 ≈10 k_BT_c；Marsiglio–Carbotte 综述 "order 10"。** 是同一条曲线在不同近似下的读数，不是矛盾。论文若要提这个数，**引 BR 原文的 2πT_c**（最严格的那个），并可加一句 Carbotte 的 7.0。

### 10.6 §10 的「读到的 vs 推断的」

**读到的（L1，本次逐页页图）**：RMP p.1039 (2.13)(2.15)、p.1040 (2.16)(2.17)(2.18) + Table I、p.1041 (2.19)(2.20)(2.23)(2.24) + Fig.2、p.1042 (2.25)–(2.29) + §II.F 开头、p.1043 (2.31)(2.32)(2.34) + Fig.3、p.1044 (2.38)–(2.44) + Fig.4 + §II.G 开头、p.1045 (2.45)–(2.50) + Fig.5、p.1046 (2.51)(2.52)(2.53) + Fig.6、p.1047 §II.I–J + Fig.7/8、p.1095 §VI.D + Fig.70/71、p.1096 §VI.E–F、p.1120–1121 §VIII.A–B（确认 §VII/VIII 无新的最优谱约束）。

**我的推断（不是 Carbotte 的陈述）**：
- 「Carbotte (2.50) 的局部⟹全局这一步没证明」——Carbotte 写的是 "local maximum" 然后 "It follows that"，**「没证明」是我的判读**；他没有明说这一步缺证明（但 p.1096 的 "would lead us to believe" / "indicating that" 强烈支持这个判读）。
- 0.2584/√2 = 0.182716 与 C_∞ 的数值一致性（§10.2 第 3 点），及由此说「C_∞ 在 1990 年已有 4 位数值」。
- 「定理 E 可以直接引 (2.20)–(2.23) 而不必重证齐次性」——需仓库侧确认定理 E 用的齐次性与 (2.20) 同形。

**新的未闭合缺口（替代原 §9.1）**：**Owen & Scalapino, Physica 55, 691–695 (1971) 正文未取到**（著录已锁定，风险评估见 §10.2 第 4 点）。这是目前唯一可能影响优先权裁决的未读文献，且我判断风险低。

**⚠ 检索引擎误导（本轮第二次，须记）**：搜 Owen–Scalapino 时，检索引擎把 KAY II Remark 3 的句子「The asymptotic behavior T_c(λ,P) ∼ (1/2π)√(𝔤(2)⟨ω²⟩λ) was first stated in [1]」端出来，暗示 [1] = Owen–Scalapino。**错**——我已在 §4.2 逐字核过，KAY II 该句的 [AD] 就是 Allen–Dynes，指其 eq.(22)。同一轮检索还把 AD 的渐近系数报成 "0.15"（那是 μ\*=0.1 的值，不是 μ\*=0 的 0.182/0.183）。**本报告任何未标 L1/L1′ 的数字与归属都不要落笔。**
