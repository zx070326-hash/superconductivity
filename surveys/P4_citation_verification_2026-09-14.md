# P4 投稿前引用核实报告（2026-09-14）

**任务**：核实主稿 `theory/derivations/P4_main_theorems.md` §9 与 `theory/paper/P4_section3_existence_uniqueness.md` 中的二手引用，以及 `P4_weak_coupling_series.md` / `P4_b1_closed_form.md` 承重的 Marsiglio / Chubukov 公式号；独立复核 b₁ 首次主张；查证 Li₂ 恒等式。

**证据等级**：**L1** = 本轮亲读原文页图（PDF 落盘后 Read pages）；**L2** = 只读到摘要/二手转述；**自算** = 本报告推导。

**本轮新取到的全文**（此前仓库记录「从未拿到」）：
- **Allen & Dynes, PRB 12, 905 (1975)**，全文 18 页 L1。路线：`harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.12.905/fulltext`（APS 扫描件，不需订阅）。
- **Combescot, PRB 42, 7810 (1990)**，前 7 页 L1。同一路线。**仓库既有 survey 把这篇的标题写错了（见 §6）。**

---

## 0. 总裁决一览

| # | 引用 | 裁决 | 需改写？ |
|---|------|------|---------|
| 1 | [AD, eq. (22)] 给 0.1827 | **编号有误。** AD eq.(22) 是 $T_c\to(\lambda\langle\omega^2\rangle)^{1/2}f(\mu^*)$，**无数值常数**。数值 0.182 在 **eq.(25)**（64 模式数值），0.180 在 eq.(24)（2×2 闭式下界）。AD 原文任何地方都没有 0.1827 或 0.183。 | **是**（见 §1.5） |
| 2 | [KAY II, Remark 3] | **确认。** p.14 逐字："The asymptotic behavior … was first stated in [AD], see their eq.(22)." 但 KAY 自己转引的 eq.(22) 也是那个无常数的形式，"0.1827262477" 是 KAY 的数，不是 AD 的。 | 措辞微调 |
| 3 | [KAY II] p.14 eq.(31) = Conjecture 2 | **确认。** p.14, eq.(31): $T_c^\sim(\lambda,P):=\frac1{2\pi}\sqrt{\mathfrak g(2)\langle\omega^2\rangle\lambda}$。 | 否 |
| 4 | 𝔎 是 ℌ 的逐元 P-平均 | **确认，且出处应改为 KAY II p.10（eq.(22) 之后的一行 + Remark 2）**，不是 Remark 3。 | 出处改 |
| 5 | [KAY III, Thm. 5] 单调性 on [0,ϖ_*) | **确认。** p.8。 | 否 |
| 6 | "would require a different strategy" | **确认，p.28**（§4.1 Summary）。逐字见 §3.3。 | 否 |
| 7 | Marsiglio eq. 37 | **确认。** PRB 98, 024523 (2018) p.6 eq.(37): $T_c=\frac{1.13}{\sqrt e}\,\omega_E\exp(-(1+\lambda)/\lambda)$。 | 否 |
| 8 | Chubukov eq. 17 | **确认。** Ann. Phys. 417, 168190 (2020), arXiv:2004.01281 p.12 eq.(17): $T_c=1.13\,e^{-1/2}\Omega_0 e^{-\frac{1+\lambda}{\lambda}}=0.252\,\Omega_0 e^{-1/\lambda}$。**原文印 0.252，我方文件写 0.2532**（见 §4.3）。 | 数字改 |
| 9 | b₁ = 49π²/144 − 3 首次主张 | **维持「文献未报告」**，但须加两条限定：(a) Combescot 1990 **明确定义了这个展开并声明只算到 O(λ⁰)**，b₁ 是他留下的空位；(b) Marsiglio 2018 的 "first order in λ" 是能隙函数的一阶不是 T_c 的一阶。 | 措辞补引 |
| 10 | Li₂(−3)+2Li₂(−2) = −π²/3 − 2ln2·ln3 | **已知恒等式的直接推论**（Abel 五项关系在 x=y=1/2 处的特例 + Landen 变换），不新颖。不应写成「新恒等式」。 | **是**（见 §5） |

---

## 1. Allen & Dynes 1975 —— 0.1827 的出处

### 1.1 完整条目
P. B. Allen and R. C. Dynes, "Transition temperature of strong-coupled superconductors reanalyzed," *Phys. Rev. B* **12**, 905–922 (1 August 1975). DOI: 10.1103/PhysRevB.12.905. Received 20 January 1975.

### 1.2 原文逐字（L1，页图）

**p.906 §I 末**：「The asymptotic behavior of $T_c$ at large λ is rigorously shown to be $\sim(\lambda\langle\omega^2\rangle)^{1/2}$ rather than $\sim\langle\omega\rangle$ as in Eq. (2).」

**p.908 eq.(21)**（N=1 闭式下界，Einstein）：
$$T_c\ge T_{c0}=\frac{\omega_E}{2\pi}\Bigl(\frac{\lambda}{1+2\mu^*(N)}-1\Bigr)^{1/2}.\tag{21}$$

**p.908 eq.(22)**（渐近**形式**，无常数）：
$$T_c\to(\lambda\langle\omega^2\rangle)^{1/2}f(\mu^*(N)),\tag{22}$$
「where $f$ is a function of $\mu^*$ only. This follows because in the large-λ limit [using Eq. (16)], information about $T_c$, λ, and phonons only enters Eq. (15) in the combination $\lambda\langle\omega^2\rangle/T_c^2$.」

**p.908 eq.(24)**（2×2 截断下界，闭式极限）：
$$T_c\ge T_{c1}=0.180(\lambda\langle\omega^2\rangle)^{1/2}\quad(\lambda>10,\ \mu^*=0).\tag{24}$$

**p.908 eq.(25)**（64 模式数值，"exact"）：
$$T_c=0.182(\lambda\langle\omega^2\rangle)^{1/2}\quad(\lambda\gtrsim10,\ \mu^*=0).\tag{25}$$
「The exact behavior in this limit,¹⁴ as deduced from computer calculations with N+1=64, is … Thus the 2×2 approximation has converged to within 1% at large λ.」脚注 14：「A result equivalent to Eq. (25) has been derived by D. Rainer [Symposium on Superconductivity and Lattice Instabilities, Gatlinburg, Tenn., 1973 (unpublished)]. This and several other results of this paper were also derived independently by J. W. Garland (private communication).」

**p.910 Fig.3 图注**：「The asymptotic formulas for $T_c/\omega_E$ are $0.18\lambda^{1/2}$ and $0.15\lambda^{1/2}$ for μ*=0 and 0.1, respectively.」

**p.905 摘要**：「in the asymptotic regime of very large λ where $T_c$ is found to equal $0.15(\lambda\langle\omega^2\rangle)^{1/2}$ (assuming μ*=0.1).」

**p.918 Fig.14 图注 / 正文**：「The actual maximum is $T_c\sim0.15\times(\lambda\langle\omega^2\rangle)^{1/2}$ (assuming μ*=0.1), occurring when λ>10.」

### 1.3 AD 全文中的常数清单（逐页核过 pp.905–922）

| 数字 | 位置 | 性质 |
|------|------|------|
| 0.180 | eq.(24) | 2×2 截断下界的 λ→∞ 极限，μ*=0，闭式 |
| 0.182 | eq.(25) | N+1=64 数值，μ*=0，"exact" |
| 0.18 | Fig.3 图注、p.910 正文 | eq.(25) 的两位舍入 |
| 0.15 | 摘要、Fig.3、Fig.14、p.918 | μ*=0.1 |
| **0.1827 / 0.183** | **全文无** | — |

**AD 原文的 0.182 是数值不是定理**；KAY 的 0.1827262477 是 KAY 自己算的 $\frac1{2\pi}\sqrt{\mathfrak g(2)}$。KAY III Remark 1（我方旧记忆）也说 AD 的 C=0.182 是 64-mode 数值，与本轮 L1 一致。

### 1.4 附带核实：Einstein 极值性（定理 E 的归属链）

**这一条本轮闭合了归属链最后一环**（旧记录：Moussa–Cohen 2006 转引 AD，AD 原文未取）。

- **p.909 左栏**：「it is interesting to ask what is the optimum shape for the function $\alpha^2F$ in the sense that $T_c$ is maximized under the conditions of fixed λ, μ*, and $\langle\omega^2\rangle$. The asymptotic formula (22) indicates that at very large λ, the shape doesn't affect $T_c$. At smaller values of λ the shape does play a role. It is clear from inspection of Eq. (15) that the eigenvalue ρ (and thus $T_c$) is maximized by making the parameters λ(n) as large as possible. … we can then ask which shape of g(ω) (for fixed second moment $\langle\omega^2\rangle$) maximizes λ(n)/λ. It is proved in Appendix A that the Einstein spectrum optimizes the integral (29) and thus is the preferred phonon spectrum.」
- **p.911 右栏**：「On this scale the Einstein spectrum is an upper bound for all n. This is proved in Appendix A and is the basis for the idea that the Einstein spectrum is optimal.」
- **Appendix A, p.920–921, (A4)**：$f(a)\equiv\int_0^\infty d\omega\,g(\omega)\frac{\omega^2}{\omega^2+a^2}\le f_E(a)\equiv\frac{\langle\omega^2\rangle}{\langle\omega^2\rangle+a^2}$，用 Wheeler–Gordon（Padé）定理 (A10) 证；p.921 明言「Thus the inequality (A13) would lead to a more general theorem on $T_c$ than the one given in the text; namely, given λ, μ*, and $\langle\omega^n\rangle$ for an arbitrary value n≥2, the maximum $T_c$ consistent with these parameters is obtained when (and only when) $\alpha^2F$ is an Einstein spectrum.」

**AD 实际证明的**：核矩阵的每个元 λ(n) 被 Einstein 谱逐点最大化（等价于我方 (A4)/Jensen 型逐元比较）。**AD 没证的**：从「逐元 ≤」到「最大本征值 ≤」这一步（因 𝔎 含负号项 −𝔎₁，逐元优势不自动传到谱），AD 只说「It is clear from inspection of Eq. (15) that the eigenvalue ρ … is maximized by making the parameters λ(n) as large as possible」——**这是断言，不是证明**。这正是我方定理 2 补的那一步。**论文可以且应当写**：「AD 1975 (Sec. II end, Appendix A) proved the elementwise inequality and asserted the eigenvalue consequence; Theorem 2 supplies the missing step and the strictness.」

### 1.5 我方用法是否需要改写

**需要。** 当前 `P4_section3_existence_uniqueness.md` 第 969–970 行与 1367–1370 行写：

> This is precisely the constant in the strong-coupling asymptotics of Allen and Dynes, $T_c\simeq0.1827\sqrt{\lambda\langle\omega^2\rangle}$ as $\lambda\to\infty$ [AD, eq. (22)]

**三处不实**：(i) AD eq.(22) 无常数；(ii) AD 从未给 0.1827；(iii) 「the constant in the strong-coupling asymptotics of Allen and Dynes」暗示 AD 给了这个常数。

**建议改为**：

> Allen and Dynes [AD, eq. (22)] showed that $T_c\to(\lambda\langle\omega^2\rangle)^{1/2}f(\mu^*)$ as $\lambda\to\infty$, and reported the numerical value $f(0)\approx0.182$ from a 64-frequency computation [AD, eq. (25)]; the sharp value $f(0)=\frac1{2\pi}\sqrt{\mathfrak g(2)}=0.1827262477\ldots$ was identified in [KAY II, Remark 3] (see also [KAY I]). Our certified enclosure … agrees with it.

References 段第 1069 行「eq. (22) for the strong-coupling asymptotic constant」改为「eq. (22) for the $\sqrt{\lambda\langle\omega^2\rangle}$ asymptotic form; eq. (25) for the numerical constant 0.182」。

同理 `P4_main_theorems.md` §9 第 330 行「未决」可销账；`P4_einstein_extremal.md` §0.A「AD 原文从未拿到」可更新为已取，并加 §1.4 的归属细节。

---

## 2. KAY II —— Remark 3、eq.(31)、逐元 P-平均

### 2.1 完整条目
M. K.-H. Kiessling, B. L. Altshuler, E. A. Yuzbashyan, "Bounds on $T_c$ in the Eliashberg theory of Superconductivity. II: Dispersive phonons," arXiv:2409.00532 (v3, revised 12 July 2025). *J. Stat. Phys.* **192**, art. 94 (2025)【卷号：KAY III 参考文献页逐字 "J. Statist. Phys. 192, art.69" 给的是 KAY I；KAY II 在 KAY III 中是 "in press"。本轮 arXiv v3 首页仍无卷号，**KAY II 卷号 192 art.94 仍为 B 级（arXiv comments 转录），投稿前请由 lead 用 DOI 10.1007/s10955-025-03470-y 或 Springer 页最终核一次**】。

### 2.2 Remark 3（p.14，L1 逐字）
> **Remark 3**: *The asymptotic behavior $T_c(\lambda,P)\sim\frac1{2\pi}\sqrt{\mathfrak g(2)\langle\omega^2\rangle\lambda}$ was first stated in [AD], see their eq.(22). Apparently they did not notice that this is presumably an upper bound on $T_c$ for all λ>0.*

**注意**：KAY 把 $\frac1{2\pi}\sqrt{\mathfrak g(2)\cdots}$ 归给 AD eq.(22)，但 AD eq.(22) 只有 $f(\mu^*)$ 这个未定函数。KAY 的说法「first stated」指的是 $\sqrt{\lambda\langle\omega^2\rangle}$ 这个**形式**。我方引 Remark 3 时须把 "the constant" 与 "the form" 分开。

### 2.3 Conjecture 2 与 eq.(31)（p.14，L1 逐字）
> **Conjecture 2**: *There is a critical temperature $T_c(\lambda,P)>0$ which for all $P\in\mathcal P$ and λ>0 is bounded above by $T_c(\lambda,P)<T_c^\sim(\lambda,P)$, with*
> $$T_c^\sim(\lambda,P):=\tfrac1{2\pi}\sqrt{\mathfrak g(2)\langle\omega^2\rangle\lambda};\tag{31}$$
> *here, $\mathfrak g(2)$ is the spectral radius of $\mathfrak G(\gamma)$ at γ=2.*
> For the numerical approximation of $\frac1{2\pi}\sqrt{\mathfrak g(2)}$ to 10 significant decimal places, see Remark 1 after Theorem 2.

**确认**：编号、页码、形式（严格 <，含存在性）全部与我方引用一致。

### 2.4 「𝔎 是 ℌ 的逐元 P-平均」——出处修正

**正确出处是 p.10，不是 Remark 3。** 两处逐字：

- **p.10，eq.(22) 之后**：「with $\mathfrak K^{(4)}(P,T)=\int_0^\infty\mathfrak H^{(4)}(\varpi)P(d\omega)$, where we introduced the abbreviation $\varpi:=\omega/2\pi T$」；矩阵 (23) 的元是 $[\![n]\!](\varpi):=\frac{\varpi^2}{n^2+\varpi^2}$。
- **p.10，Remark 2**：「We have stated our Theorem 4 in a way that makes it obvious that the matrices $\mathfrak K^{(N)}(P,T)$ are averages w.r.t. $P(d\omega)$ of matrices $\mathfrak H^{(N)}(\varpi)$ of non-dispersive models, with $\varpi=\omega/2\pi T$. However, with the sole exception of the eigenvalue $\mathfrak k^{(1)}(P,T)$ of the first approximation, all eigenvalues $\mathfrak k^{(N)}(P,T)$ for N>1 are not simply averages w.r.t. $P(d\omega)$ of the pertinent largest eigenvalues $\mathfrak h^{(N)}(\varpi)$ of $\mathfrak H^{(N)}(\varpi)$. The same is true for $\mathfrak k(P,T)$ in relation to $\mathfrak h(\varpi)$.」

（N=∞ 的对应定义在 §4.1 eqs.(39)–(42)，旧记忆已记。）**引用应写 [KAY II, p. 10, eq. (23) and Remark 2]**。当前 References 段已有 "[KAY II, Remark 2, p. 10]"，只需把「逐元 P-平均」这一句挂到它而不是 Remark 3。

---

## 3. KAY III —— Thm. 5 与 "different strategy"

### 3.1 完整条目
M. K.-H. Kiessling, B. L. Altshuler, E. A. Yuzbashyan, "Bounds on $T_c$ in the Eliashberg theory of Superconductivity. III: Einstein phonons," arXiv:2409.02121 (v2). *J. Stat. Phys.* **192**, art. 93 (2025)【同上，卷号 B 级】。

### 3.2 Theorem 5（p.8，L1 逐字）
> **Theorem 5**: *For all N∈ℕ the eigenvalues $\mathfrak h^{(N)}(\varpi)$ increase strictly monotonically with $\varpi\in[0,\varpi_*)$. Moreover, $\varpi_*\ge\sqrt2$. As a consequence, the map $\varpi\mapsto L_E(\varpi)$ is strictly monotonic decreasing on $[0,\varpi_*]$, with $\varpi_*\ge\sqrt2$. Thus the portion of the critical surface $\mathscr S_c$ over the region $\{\lambda\ge\lambda_*\}$ in the positive (λ,Ω)-quadrant is also a graph, yielding the critical temperature $T_c(\lambda,\Omega)$, viz.*
> $$\mathscr S_c|_{\lambda\ge\lambda_*}=\{(\lambda,\Omega,T)\in\mathbb R^3_+:T=T_c(\lambda,\Omega),\lambda\ge\lambda_*\}.\tag{17}$$
> *Moreover, $T_c(\lambda,\Omega)=\Omega f(\lambda)$, where f(λ) is continuous and strictly monotonically increasing for λ≥λ_*. Furthermore, λ_*<0.7670….*

紧接着 p.8：「While we have not succeeded in showing that the map $\varpi\mapsto\mathfrak h(\varpi)$ is strictly monotonic increasing for all $\varpi\in\mathbb R_+$, our lower bounds to $\mathfrak h(\varpi)$ stated explicitly in Theorem 3 for N∈{1,2,3,4} all are strictly monotonic increasing with ϖ. This is manifestly obvious only for $\mathfrak h^{(1)}(\varpi)$. For $\mathfrak h^{(2)}(\varpi)$ this is a consequence of Proposition 9 in [KAYb]. For $\mathfrak h^{(3)}(\varpi)$ and $\mathfrak h^{(4)}(\varpi)$ the monotonicity for λ<λ_* is vindicated through visual inspection of the plots (see below).」

**确认**：Thm. 5 陈述的是 N 有限的 [0,ϖ_*) 单调 + N=∞ 的 [0,√2] 单调（via KAY II Prop 1 的收敛散文）。

### 3.3 "different strategy"（p.28，§4.1 Summary，L1 逐字）
> 「By taking the dispersionfree limit of our results in [KAYb], we showed that all our upper approximations to the map $\varpi\mapsto L_E(\varpi)$, and this map itself, are strictly monotone decreasing for $\varpi\in[0,\sqrt2]$. We also supplied compelling evidence for the conjecture that the map $\varpi\mapsto L_E(\varpi)$ is strictly monotone decreasing for all $\varpi\in[0,\infty)$, but to prove it would require a different strategy.」

**确认，页码 p.28**。我方 References 段「[KAY III, Thm. 5, p. 8] and [KAY III, p. 28]」正确。注意原文的对象是 $L_E(\varpi)=1/\mathfrak h(\varpi)$ 的**递减**，等价于 𝔥 递增；引用时按原文措辞。

同页还有可直接引的：「what we just wrote proves that a unique critical temperature $T_c(\lambda,\Omega)$ in the Eliashberg model with Einstein phonons is mathematically well-defined in terms of the untruncated linearized Eliashberg gap equations whenever λ>0.7670…. Also this λ value is not a sharp boundary but a consequence of our method of proof.」以及「$C=\frac1{2\pi}\sqrt{\mathfrak g(2)}=0.1827262477\ldots$」——**这是 0.1827262477 十位数在 KAY III 的原始位置（p.28），可作为该常数的引用点**。

---

## 4. Marsiglio eq.(37) 与 Chubukov eq.(17)

### 4.1 Marsiglio 2018
**条目**：F. Marsiglio, "Eliashberg theory in the weak-coupling limit," *Phys. Rev. B* **98**, 024523 (2018). arXiv:1807.04907. DOI: 10.1103/PhysRevB.98.024523.

**p.6 eq.(37)，L1 逐字**：
$$T_c=\frac{1.13}{\sqrt e}\,\omega_E\exp\bigl(-(1+\lambda)/\lambda\bigr).\tag{37}$$
前句：「Following with the same type of analysis as that leading to Eq. (26) and to Eq. (33) we find here that」。后句：「As previously mentioned, Eq. (37) can of course be written with a −1/λ in the exponential, along with a prefactor denominator of $e^{3/2}$ instead of $\sqrt e$. However, the present form more explicitly shows the role of the "normal-state" renormalization that gives rise to the usual 1+λ factor, along with the not-so-usual $\sqrt e$ denominator in the pre-factor.」

其它承重公式（L1）：eq.(13) $T_c=1.13\omega_E\exp(-1/\lambda)$；eq.(14) $T_c=1.13\omega_E\exp(-(1+\lambda)/\lambda)$；eq.(26) $T_c=\frac{1.13}{\sqrt e}\omega_E\exp(-1/\lambda)$（Z≡1）；eq.(24) $I_0=\lambda I_0^2-I_0/2$；eq.(25) $I_4\approx\tfrac12$。**p.3 右栏逐字**：「One can immediately write this like Eq. (13) but with reduced pre-factor $1.13e^{-1}$. This is *not* what is meant when we stated that the pre-factor in Eliashberg theory is actually modified from the BCS result — but rather an additional change occurs.」

**确认**：我方 c₀ = γ + ln4 − 3/2 的推导以 eq.(37) 为锚，公式号与内容对。**Fig.4 图注**：「This latter result becomes essentially exact for λ≲0.2」——这是 Marsiglio 对 (37) 有效范围的自陈，可引来对照 b₁λ 修正量（λ=0.2 时 +7%）。

### 4.2 Chubukov–Abanov–Esterlis–Kivelson 2020
**条目**：A. V. Chubukov, A. Abanov, I. Esterlis, S. A. Kivelson, "Eliashberg theory of phonon-mediated superconductivity — When it is valid and how it breaks down," *Ann. Phys.* **417**, 168190 (2020). arXiv:2004.01281. DOI: 10.1016/j.aop.2020.168190.

**p.12 eq.(17)，L1 逐字**：
$$T_c=1.13\,e^{-1/2}\,\Omega_0\,e^{-\frac{1+\lambda}{\lambda}}=0.252\,\Omega_0\,e^{-\frac1\lambda}.\tag{17}$$
前句：「As a consequence, both $T_c$ and Δ are obtained within ETh with the exact prefactors, as has been discussed in several papers, using different computational tools¹⁸⁻²⁶. The result is」。

**p.12 eq.(19)**：$T_c\approx0.1827\,\alpha$；**p.13**：「Eq. (19) was first obtained in Ref.⁹ [= AD]. These authors expressed the critical temperature as $T_c\sim\Omega_0\sqrt\lambda$…」——**Chubukov 也把 0.1827 归给 AD，与 AD 原文不符**（AD 是 0.182）。若我方论文要说「0.1827 在文献中的来源」，不能引 CAEK 当作 AD 给 0.1827 的证据。

### 4.3 数字不一致
`P4_weak_coupling_series.md` 第 170 行：「1.13/1.6487·ω_E e^{−1/λ}·e^{−1} ≈ 0.2532」。**自算**：$1.13/(\sqrt e\cdot e)=1.13/4.4817=0.25214$；用 $2e^\gamma/\pi=1.1339$ 得 $0.25301$。**0.2532 两条都对不上**（1.13/e^{3/2} 给 0.2521；精确前因子给 0.2530）。CAEK 印 0.252（对应 1.13）。建议正文统一用 $\frac{2e^{\gamma}}{\pi}e^{-3/2}=0.25301\ldots$，并注明 CAEK 的 0.252 是 1.13 舍入所致。

---

## 5. Li₂(−3) + 2Li₂(−2) = −π²/3 − 2 ln2·ln3 —— 已知

### 5.1 裁决
**不新颖**。它是 Abel 五项关系在 x=y=1/2 处的特例（给出 Li₂(¼)+2Li₂(⅓) 的归约）再经 Landen 变换换到负实轴的直接推论。目标组合本身检索未见显式印刷处，但推导只有两步、全部是 Lewin/Zagier 教科书内容，必须按「已知」处理。

### 5.2 推导（自算，全部数值复核到 7 位）

**步 1（Abel 五项关系，x=y=1/2）**。五项关系的一种标准形式：
$$\mathrm{Li}_2(x)+\mathrm{Li}_2(y)=\mathrm{Li}_2(xy)+\mathrm{Li}_2\!\Bigl(\frac{x(1-y)}{1-xy}\Bigr)+\mathrm{Li}_2\!\Bigl(\frac{y(1-x)}{1-xy}\Bigr)+\ln\!\frac{1-x}{1-xy}\,\ln\!\frac{1-y}{1-xy}.$$
取 x=y=1/2：xy=1/4，x(1−y)/(1−xy)=1/3，ln((1/2)/(3/4))=ln(2/3)，故
$$2\mathrm{Li}_2(\tfrac12)=\mathrm{Li}_2(\tfrac14)+2\mathrm{Li}_2(\tfrac13)+\ln^2(\tfrac23),$$
$$\Longrightarrow\quad \mathrm{Li}_2(\tfrac14)+2\mathrm{Li}_2(\tfrac13)=\tfrac{\pi^2}{6}-2\ln^2 2+2\ln 2\ln 3-\ln^2 3.\tag{★}$$
数值：左边 = 0.2676526 + 2×0.3662132 = 1.0000791；右边 = 1.6449341 − 0.9609060 + 1.5230000 − 1.2069490 = 1.0000791。✓

**步 2（Landen + 反射）**。Landen：$\mathrm{Li}_2(-z)=-\mathrm{Li}_2\!\bigl(\tfrac{z}{1+z}\bigr)-\tfrac12\ln^2(1+z)$；反射：$\mathrm{Li}_2(u)+\mathrm{Li}_2(1-u)=\tfrac{\pi^2}{6}-\ln u\ln(1-u)$。
- z=3：$\mathrm{Li}_2(-3)=-\mathrm{Li}_2(\tfrac34)-\tfrac12\ln^2 4$；
- z=2：$\mathrm{Li}_2(-2)=-\mathrm{Li}_2(\tfrac23)-\tfrac12\ln^2 3$；
- $\mathrm{Li}_2(\tfrac34)=\tfrac{\pi^2}{6}-\mathrm{Li}_2(\tfrac14)-\ln\tfrac34\ln\tfrac14$，$\mathrm{Li}_2(\tfrac23)=\tfrac{\pi^2}{6}-\mathrm{Li}_2(\tfrac13)-\ln\tfrac23\ln\tfrac13$。

代回并整理（记 a=ln2, b=ln3；对数项合计 (b−2a)(−2a) − 2a² + 2(a−b)(−b) − b² = 2a² + b² − 4ab）：
$$\mathrm{Li}_2(-3)+2\mathrm{Li}_2(-2)=-\tfrac{\pi^2}{2}+\bigl[\mathrm{Li}_2(\tfrac14)+2\mathrm{Li}_2(\tfrac13)\bigr]+2\ln^2 2+\ln^2 3-4\ln 2\ln 3.$$
代入 (★)：二次项全部相消，得
$$\boxed{\ \mathrm{Li}_2(-3)+2\mathrm{Li}_2(-2)=-\tfrac{\pi^2}{3}-2\ln 2\,\ln 3\ }$$
数值：左边 = −1.9393754 + 2×(−1.4367464) = −4.8128681；右边 = −3.2898681 − 1.5230000 = −4.8128681。✓（7 位一致）

### 5.3 出处与建议措辞
(★) 正是 Lewin 意义下的「dilogarithm ladder」型有理点关系（1/4 与 1/3 混合基）；(★)+Landen+反射 ⟹ 目标式，整条链是 Abel 五项关系的初等推论。**检索未找到目标式（−3, −2 组合）的显式印刷处**，但它的推导只有两步、全是教科书内容（Lewin 1981 Ch.1 的 (Landen) 与 five-term relation；Zagier 2007 §I.1 有五项关系与 Landen 的现成陈述），**必须按「已知」处理**。

建议措辞：「The combination $\mathrm{Li}_2(-3)+2\mathrm{Li}_2(-2)=-\pi^2/3-2\ln2\ln3$ is a two-step consequence of the Abel five-term relation and Landen's transformation [Lewin 1981, Ch. 1; Zagier 2007, §I.1]; it explains why $b_1$ is a rational multiple of $\pi^2$ although $\mathrm{Li}_2(-2)$ and $\mathrm{Li}_2(-3)$ are individually irreducible.」

**引用条目**：L. Lewin, *Polylogarithms and Associated Functions*, North-Holland, New York, 1981（ISBN 0-444-00550-1），Ch. 1；D. Zagier, "The Dilogarithm Function," in *Frontiers in Number Theory, Physics, and Geometry II*, Springer 2007, pp. 3–65（免费 PDF 在 MPIM 页）。**Lewin 书本身本轮未取到页图**（无免费全文），**不要引具体式号**；引 Zagier 2007 更稳妥。

---

## 6. b₁ 首次主张的独立复核

### 6.1 裁决
**维持「文献未报告 b₁ 的值」**，但既有 survey 有一处**实质性错误**必须改：

**`surveys/P4_weak_coupling_correction_lit_2026-09-13.md` §5 把 Combescot PRB 42, 7810 (1990) 的标题写成 "Analytical results for the longitudinal and transverse Tsuneto functions…" 并裁为「Wrong topic (unconventional SC)」。这是错的。** 该篇真实标题（L1 页图）：

> R. Combescot, "Critical temperature of superconductors: Exact solution from Eliashberg equations on the weak-coupling side," *Phys. Rev. B* **42**, 7810–7824 (1 November 1990). DOI: 10.1103/PhysRevB.42.7810. Received 16 April 1990.

**它正是 b₁ 问题的直接前身。**

### 6.2 Combescot 1990 说了什么（L1，pp.7810–7816）

**p.7811 左栏，逐字**：「More specifically, in this limit of small λ we have the following expansion in terms of the coupling constant: $\ln T_c=-1/\lambda+\ln A+O(\lambda)$. The first term is already present in the BCS result $T_c/\omega_c=1.13e^{-1/\lambda}$, but since the cutoff $\omega_c$ is ill defined, we have only an order of magnitude information on the constant A. Further work has not yet given a complete expression for A and its $\alpha^2F(\omega)$ dependence… The purpose of the present paper is to provide the exact result for A as a functional of $\alpha^2F(\omega)$ and to study it. More generally, we will show how one can generate systematically an expansion for $\ln T_c$ on the weak-coupling side.」

**p.7813 左栏，逐字**：「However, in this paper we want only to work out the expression of $T_c$ to lowest order in λ, i.e., neglecting all the terms which vanish for λ→0.」

**p.7813 eq.(21)**（Einstein，μ*=0）：
$$T_c=\frac2\pi e^{C-1/2}\,\omega_0\,e^{-(1+\lambda)/\lambda},\tag{21}$$
「where the coefficient $2e^C/\pi\simeq1.13$ is the one found in the BCS model. We have checked, down to λ=0.2, that the numerical solution extrapolates properly to a result which agrees with the prefactor of the exponential in Eq. (21).」

**p.7813 右栏**：「The next-order correction in the average (see Appendix A) would be $\langle\langle X_{i-1}X_i\rangle\rangle-\langle\langle X_i\rangle\rangle^2$. ln(4/r) and higher-order terms introduce higher-order correlations. This correction vanishes in the weak-coupling limit.」

**p.7816 eq.(56)**（一般谱 + μ*）：$T_c=\frac2\pi e^{C-1/2}\omega_{\log}\exp[-\frac{1+\lambda}{\lambda-\mu^*}-\frac{\lambda\mu^*}{(\lambda-\mu^*)^2}-\frac{\lambda^2}{(\lambda-\mu^*)^2}R]$，R≥0 形状因子。

**结论**：
1. Combescot 1990 **首次**给了 Einstein 谱 μ*=0 的精确前因子 $\frac2\pi e^{\gamma-1/2}=e^{-1/2}\times1.134$，即我方 c₀ = γ + ln4 − 3/2 的**原始出处**（早于 Marsiglio 2018 二十八年；Marsiglio 把 √e 归功 Karakozov et al. 1975，Combescot 引用的前身是 ref.16 "A short account … already been published"）。**我方文件与主稿目前只引 Marsiglio/Chubukov，应补引 Combescot 1990 eq.(21)。**
2. Combescot **明确定义了 $\ln T_c=-1/\lambda+\ln A+O(\lambda)$ 展开、声明方法可系统推高阶、并明确只做到 O(λ⁰)**。他的 Appendix A 累积展开给出了 O(λ) 项的**结构**（$\langle\langle X_{i-1}X_i\rangle\rangle-\langle\langle X_i\rangle\rangle^2$ 型关联），但**没有算出数值或闭式**。
3. **注意展开变量**：Combescot 的 O(λ) 是在 $\ln T_c$ 里，且指数写成 $-(1+\lambda)/\lambda=-1/\lambda-1$。我方的 $T_c=A\omega_E e^{-1/\lambda}(1+b_1\lambda)$ ⟺ $\ln T_c=-1/\lambda+\ln A+b_1\lambda+O(\lambda^2)$，**与 Combescot 的展开逐项同一**。所以 b₁ 就是 Combescot p.7811 那个 O(λ) 的系数。

### 6.3 其它来源（本轮再核）
- **Marsiglio 2018** 摘要句 "extend the analytical results to first order in the coupling constant" 指的是**能隙函数** Δ(ω_m) 的一阶（eq.(33)/(39)），不是 T_c 的一阶。T_c 侧只到 eq.(37) 前因子。Fig.4 图注自陈 (37) 「essentially exact for λ≲0.2」。**无 b₁。**
- **Mirabi–Boyack–Marsiglio PRB 101, 064506 (2020)**（L2 摘要）：1/√e 在实轴的确认，无 O(λ)。
- **CAEK 2020** eq.(17) 后：「A recipe for computing the weak coupling $T_c$ for an arbitrary non-critical bosonic propagator has been given in²³」（ref.23 = Wang–Abanov–Altshuler–Yuzbashyan–Chubukov PRL 117, 157001）。该 recipe 是前因子级，旧 survey 已裁「错误区域」。**无 b₁。**
- **Karakozov–Maksimov–Mashkov JETP 41, 971 (1975)**：旧 survey 只读到扫描件 6 页，无 Einstein λ-扫描。本轮未重跑（无新路线）。
- **Hainzl–Seiringer 系**：Frank–Hainzl–Naboko–Seiringer J. Geom. Anal. 17, 559 (2007) 与 Roos 2407.00796 均为瞬时 BCS 核；**本轮再搜确认无推迟核严格弱耦合结果**。
- **Freericks–Scalapino 1994**：旧 survey 已裁「Holstein QMC，形式太复杂」，本轮不重复。

### 6.4 建议措辞
> The leading weak-coupling asymptotics $T_c=\frac{2}{\pi}e^{\gamma-1/2}\,\omega_E\,e^{-(1+\lambda)/\lambda}$ for an Einstein spectrum with μ*=0 was obtained by Combescot [Co90, eq. (21)], who set up the expansion $\ln T_c=-1/\lambda+\ln A+O(\lambda)$ and computed $A$ exactly, stating explicitly that the $O(\lambda)$ term was left aside; the same prefactor was rederived on the Matsubara axis by Marsiglio [Ma18, eq. (37)] and is quoted in [CAEK20, eq. (17)]. To our knowledge the coefficient of the $O(\lambda)$ term, $b_1=49\pi^2/144-3$, has not previously been reported, either in closed form or numerically.

---

## 7. 仍拿不到 / 仍为 B 级的项

| 项 | 状态 | 说明 |
|----|------|------|
| KAY II / III 期刊卷号与 art. 号 | B 级 | arXiv v3/v2 首页无卷号；KAY III 参考文献确认 KAY I = 192 art.69，KAY II "in press"。请 lead 用 Springer DOI 页最终核。 |
| Lewin 1981 具体式号 | 放弃引式号 | 无免费全文。引书 + 章即可（Ch. 1），或改引 Zagier 2007（免费 PDF）。 |
| Karakozov 1975 是否含 O(λ) 系数 | 未解 | 扫描件不可读；但 Combescot 1990 p.7811 已断言「Further work has not yet given a complete expression for A」，间接说明 1975 年 KMM 未给完整 A，更不会有 b₁。 |
| Combescot 1990 pp.7817–7824（Appendix A–E） | 未读 | 前 7 页已足以裁决；Appendix A 的累积展开若 lead 想核 O(λ) 项结构可再取（同一 harvest URL，pages 8–15）。 |

---

## 8. 需要改写的文件清单

1. `theory/paper/P4_section3_existence_uniqueness.md` 第 969–970、1069–1070、1367–1370 行：按 §1.5 改。
2. `theory/derivations/P4_main_theorems.md` §9 第 330 行「未决」→ 销账，指向本报告。
3. `theory/derivations/P4_einstein_extremal.md` §0.A：AD 原文已取；补 §1.4 的「AD 证了逐元、断言了谱、我方补证谱」归属句。
4. `theory/derivations/P4_weak_coupling_series.md` 第 104、170 行：0.2532→0.2530；补引 Combescot 1990 eq.(21)。
5. `theory/derivations/P4_b1_closed_form.md` 第 190–194 行：Li₂ 恒等式改为「经典推论」措辞（§5.3）。
6. `surveys/P4_weak_coupling_correction_lit_2026-09-13.md` §5：**标题与裁决错误**，改为本报告 §6.1–6.2。
7. 主稿 References 加：[Co90] Combescot 1990；[Le81] Lewin 1981 或 [Za07] Zagier 2007。

---

## 9. 本轮方法记录

- **`harvest.aps.org/v2/journals/articles/10.1103/<DOI 尾>/fulltext`** 对 1975 与 1990 的 PRB 均直接返回全文扫描 PDF（无订阅、无 CAPTCHA）。WebFetch 落盘后 Read pages 读页图。这是此前四轮（APS 付费墙 / OSTI / ADS / 作者主页）都没试过的路线，**AD 1975 与 Combescot 1990 一次取到**。
- arXiv 的 `/abs/` 页只给元数据；`/pdf/` 落盘后 Read 可读。`/pdf/<id>v3` 对 2004.01281 返回 404（该文只有 v1/v2），不带版本号即可。
