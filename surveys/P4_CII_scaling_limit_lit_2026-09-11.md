# P4 · 段 C-II（x ≥ X 的解析段）文献核查：Mityagin–Sobolev / Widom 标度极限框架能否覆盖 ℌ(ϖ)

检索 agent：researcher　日期：2026-09-11　派单：`cii-lit`（主会话消息，branch p4-conj2）
上游：`theory/derivations/P4_einstein_extremal.md` §0 / §1 / §1.B(B1–B9) / §1.C / §1.A.2–1.A.5、`surveys/P4_einstein_extremal_lit_2026-09-11.md` §3

**证据等级**：**A** = 原文逐字（附式号/定理号）；**B** = 我方对 A 级材料的确定性推算，标 `[自算]`；**C** = 仅摘要级 / 二手 / 正文未取到。

**本轮一手材料**

| 文献 | 取得方式 | 读到的范围 | 等级 |
|---|---|---|---|
| Mityagin–Sobolev, arXiv:1106.0127 = J. Spectral Theory **1** (2012) 443–460 | ar5iv HTML（四次定点取） | 摘要、§1.1–1.3、Thm 1、Thm 3、Lemma 2、§6.1–6.3（Conj 14）、Prop 15/16、参考表 | A |
| Sobolev, arXiv:1609.02052 | ar5iv HTML | 摘要、Condition 1.1、Thm 1.2、Lemma 2.2/2.3/2.4、Lemma 3.1、对 Widom 1961/1963 的转述 | A |
| Frank–Hainzl–Naboko–Seiringer, arXiv:0704.3564 = JGA **17** (2007) 559–567 | ar5iv HTML | Thm 1(i)(ii)(iii)、eq.(7)、Lemma 2、特征值扰动估计 | A |
| Stubbe, arXiv:2605.24694 | arXiv HTML v1 | Thm 1.2 (1.17)–(1.19)、Thm 1.9 (1.21)(1.22)、Cor 4.3 | A |
| Widom TAMS **100** (1961) 252–262；TAMS **106** (1963) 391–414 | **ams.org 全站 403（见 §5）** | **正文未取到** | **C** |
| Mityagin, arXiv:0803.3159 = Israel J. Math. **181** (2011) 1–28 | 经 1106.0127 §1.2 + 参考表 | 仅转述 | C |
| Montgomery–Vaughan 加权 Hilbert 不等式最优常数（arXiv:2203.14950） | arXiv abs | 摘要 | C |

---

## §1（问 1）Mityagin–Sobolev arXiv:1106.0127 到底证了什么 [A]

**题名/出处**（逐字）：*A family of anisotropic integral operators and behaviour of its maximal eigenvalue*，B. S. Mityagin, A. V. Sobolev，**Journal of Spectral Theory 1 (2012), Issue 4, 443–460**（arXiv 元数据 16 页，MSC 45C05 / 47A75）。**注意**：既有 survey 把它写成 "…and behaviour of its maximal eigenvalue"，题名正确；期刊与年份此前未记录，本轮补齐。

### 1.1 算子族（逐字式号）

> (1)　"$\mathbf{K}_{\beta}u(x)=\int K_{\beta}(x,y)u(y)dy$"
> (2)　"$K_{\beta}(x,y)=\frac{1}{\pi}\frac{1}{1+(x-y)^{2}+\beta^{2}\Theta(x,y)}$"
> (3)　"$\Theta(tx,ty)=t^{\gamma}\Theta(x,y),\ \gamma>0$"
> (4)　"$c\leq\Theta(x,y)\leq C,\ |x|^{2}+|y|^{2}=1$" 与 "$\Theta(x,y)=\Theta(y,x)$"

域：**L²(ℝ)**（整条实轴，非半轴）。条件 (3)(4) "guarantee that the operator $\mathbf{K}_{\beta}$ is self-adjoint and compact"。

### 1.2 Theorem 1（逐字）

假设：γ≥1，Θ 满足 (3)(4) 及对角附近的 Lipschitz 条件 (8)
> "$|\Theta(t,1)-\Theta(1,1)|\leq C|t-1|,\ t\in(1-\epsilon,1+\epsilon)$"（另一支在 (−1,−1) 处对称）

结论两句：
> "$\lim_{\beta\to 0}\beta^{-\frac{2}{\gamma+1}}(1-{\sf{M}}_{\beta})=\lambda_{1}$"
> "the rescaled eigenfunctions $\alpha^{-\frac{1}{2}}\Psi_{\beta}(\alpha^{-1}\cdot)$ … converge in norm to $\phi_{1}$"，α = β^{2/(γ+1)}

§1.2 的等价写法给出**余项阶**：
> "$\mathsf{M}_{\beta}=1-\lambda_{1}\beta^{2/(\gamma+1)}+o(\beta^{2/(\gamma+1)}),\ \beta\to 0$"

**⟹ 余项只是 little-o，无任何定量速率；本征函数收敛只是「按范数收敛」，无速率。** 这两条是本轮对段 C-II 最重要的负面信息（见 §4）。

### 1.3 模型算子 A（逐字）

> (6)　"$(\mathbf{A}u)(x)=|D_{x}|u(x)+2^{-1}\theta(x)u(x)$"，θ(x)=Θ(x,x)
> (7)　符号 "$a(x,\xi)=|\xi|+2^{-1}\theta(x)$"
> 域 "$D(\mathbf{A})=\textup{H}^{1}(\mathbb{R})\cap\textup{L}^{2}(\mathbb{R},|x|^{2\gamma})$"

**它怎么来的**：§1.3 的**幺正重标** "$(U_{\alpha}f)(x)=\alpha^{-1/2}f(\alpha^{-1}x)$"，令 β²=α^{γ+1}，得 "$\mathbf{K}_{\beta}=U_{\alpha}^{*}\mathbf{B}_{\alpha}U_{\alpha}$"；**Theorem 3** 把 Thm 1 改写成 "$\lim_{\alpha\to 0}\alpha^{-1}(1-\mu_{\alpha})=\lambda_{1}$" + ψ_α→φ₁ 按范数。1−K_β 的符号在重标后近似为 α·a(x,ξ)，其中 |ξ| 来自差核 1/(1+(x−y)²) 的 Fourier 变换 e^{−|ξ|} 在 ξ=0 处的**尖点**，θ/2 来自 β²Θ。

### 1.4 其它编号结果（一行式）

Lemma 2（半群 e^{−tA} positivity improving ⟹ λ₁ 单重、φ₁>0）；Lemma 4（去对称化误差 ‖K_β−K_β^{(l)}‖ ≤ β 的幂）；Lemma 5/6/Cor 7（余项 e_α^{(1)}、E_α^{(2)} 的 O(α) 界）；Lemma 8（limsup 上界）；Lemma 9–13（ψ_α 的加权有界性、Fourier/空间局域化、弱收敛升级为范数收敛）；**Prop 15 = 经典 Schur 检验**；**Prop 16 = Cordes 判据**（PDO 的 L² 有界性）。

**§6 Miscellaneous 三条 [A]**：6.1 **Conjecture 14**："$\lim_{\beta\to 0}\beta^{-2/(\gamma+1)}(1-\mathsf{M}_{j,\beta})=\lambda_{j}$"（全谱版，仍未证；文中说其灵感来自 Widom 1961）；6.2 作者自陈**不知道** K_β 在 β>0 时是否有负本征值（"can't say whether or not"）；6.3 Θ 为偶函数时奇偶子空间不变，高阶本征值的大小关系是开问题。

### 1.5 物理来源与前身 [A]

参考表逐字：**[6]** P. Krotkov, A. Chubukov, PRL **96** (2006) 107002–107005；**[7]** 同作者 PRB **74** (2006) 014509–014524；**[9]** B. Mityagin, *An anisotropic integral operator in high temperature superconductivity*, **Israel J. Math. 181, No. 1 (2011), 1–28**；**[15]** H. Widom, *Extreme eigenvalues of translation kernels*, TAMS **100** (1961) 252–262。

§1.2 逐字：
> "The top eigenvalue of $\mathbf{K}_{\beta}$ was studied by B. Mityagin in [9]"
> "only the two-sided bound … $c\beta^{2/(2\sigma+1)}\leq 1-\mathsf{M}_{\beta}\leq C\beta^{2/(2\sigma+1)}$ … was proved."
> "J. Adduci found an approximate numerical value $\lambda_{1}=0.978\ldots$"
> "Similar eigenvalue asymptotics were investigated by H. Widom in [15]"
> "Some ideas of this paper are used in the proof of Theorem 1."

**⟹ 问 1 的「用了哪些工具」答案**：Widom 差核方法（**仅「some ideas」，非整段搬运**）+ PDO 符号演算（Cordes, Prop 16）+ Schur 检验（Prop 15）+ 正性改进半群（Lemma 2）。**既有 survey §3.3 的技术清单逐条属实**，本轮补上编号。

### 1.6 直接回答派单的两个具体问题 [A]

- **权 (st)^{−1/2}：没有。** 全文不出现任何 (xy)^{−1/2} 型权，K_β 是纯差核加 β²Θ 的扰动。
- **对角相减（我方 −L1）：没有对应物。** 文中唯一的「对角」是 θ(x)=Θ(x,x)，它是**去对称化**（Lemma 4 把 Θ(x,y) 换成 Θ(x,x)）后进入**模型算子的势项**，与我方从算子里减掉一个对角矩阵是两回事。

---

## §2（问 2）Widom 1961 / 1963 —— 正文未取到，只能经转述 [C]

**取件失败**：`ams.org/journals/tran/1961-100-02/…pdf`、`…/1963-106-03/…pdf` 及两个 landing page 全部 **HTTP 403**。既有 survey §3.3 标注的「AMS 免费 PDF」**本轮无法验证，应降级为未核实**（该标注此前也非 A 级）。

**经 Sobolev 1609.02052 的转述（A 级的「Sobolev 怎么说」，不是 A 级的「Widom 怎么写」）**：

> "H. Widom considered the function $V$ which was the indicator of an interval"
> "$a=a(\xi),\ \xi\in\mathbb{R}$, having one global maximum at $\xi=0$"
> "$\lambda_{\alpha}^{(n)}=A_{0}-\alpha^{\gamma}\Psi\mu^{(n)}+o(\alpha^{\gamma}),\ \alpha\to 0$"
> "$\mu^{(n)}$ … are eigenvalues of the fractional Dirichlet Laplacian $(-\Delta)^{\gamma/2}$ on $I$"

1963 年多维版：
> "A multi-dimensional analogue of this result was obtained by H. Widom"
> "We omit its formulation for the sake of brevity"

**结构判读** `[自算，B]`：Widom 框架里 λ_max(区间长 L) **从下方单调升向符号上确界 A₀**，修正阶 L^{−γ}。故在该框架内 λ_max(L) 自动是**增且凹**的（−cL^{−γ} 增且凹）——但这是**饱和型**渐近。**我方 λ(L)=h(e^L) 线性发散、不饱和**（`P4_einstein_extremal.md` (B1)(B2)），对应「符号无有限上确界」，**不在 Widom 的核类内**（他要求 k∈L¹ 从而符号有界）。**⟹ 问 2 的「凹性/单调性是否可搬」答案：形式上 Widom 框架确实给凹性，但适用前提我方不满足；不可搬。**

**Landau–Widom / Böttcher–Silbermann / Basor–Ehrhardt**：本轮检索在这一支只找到**行列式**渐近（Wiener–Hopf±Hankel det 的 Fisher–Hartwig 型结果），**没有任何关于顶本征值对区间长度的二阶导数或凹性的结果**。记为 **文献零处**（问 4 前半的裁决）。

---

## §3（问 3）Sobolev 1609.02052 与 Mityagin 0803.3159 各加了什么

### 3.1 Sobolev 1609.02052 [A]

算子族 **B_α = V A_α V**，"$B_{\alpha}=V\mathcal{F}^{*}a_{\alpha}\mathcal{F}V,\ a_{\alpha}(\xi)=a(\alpha\xi)$"。**Condition 1.1** 逐字要点：
> "$a$ and $V$ are real-valued $L^{\infty}$-functions"
> "attain their global maxima only at $\xi=0$ and $x=0$"
> "$-V_{0}+c\leq V(x)\leq V_{0}$, $x$ a.e., with a positive constant $c$"
> 局部行为 "$V(x)=V_{0}-\Phi_{\beta}(x)+o(|x|^{\beta})$"、"$a(\xi)=A_{0}-\Psi_{\gamma}(\xi)+o(|\xi|^{\gamma})$"，Φ_β、Ψ_γ 分别 β、γ 次齐次正函数

**Theorem 1.2**："$\lim_{\alpha\to 0}\alpha^{-\sigma}(A_{0}V_{0}^{2}-\lambda_{\alpha}^{(n)})=\mu^{(n)}$"，"$\frac{1}{\sigma}=\frac{1}{\beta}+\frac{1}{\gamma}$"；模型算子符号 "$t(x,\xi)=V_{0}^{2}\Psi_{\gamma}(\xi)+2A_{0}V_{0}\Phi_{\beta}(x)$"。

**对我方最重要的一条**：**这篇是唯一带乘性权 V 的版本**，形式上最接近 (st)^{−1/2}；但 **Condition 1.1 要求 V 有界、全局极大在单点、且 −V₀+c ≤ V**。我方 s^{−1/2} 在 s→0 **无界**，条件**硬失效** `[自算，B]`。

余项与本征函数：仍是 "$+o(\alpha^{\sigma})$"，**无速率**；本征函数只有 "the sequence $\psi_{\alpha_{k}}^{(n)}$ converges to $\psi^{(n)}$ in norm"（Lemma 3.1），局域化界 "$\liminf_{\alpha\to0}\|\hat\psi_{\alpha}^{(n)}\chi_{R}\|^{2}\geq 1-C\mu^{(n)}R^{-\gamma}$"。**同样无速率。** 全文**未出现** "Landau–Widom"。

### 3.2 Mityagin 0803.3159 [C，经 1106.0127 转述]

只证了**双侧界** cβ^{2/(2σ+1)} ≤ 1−M_β ≤ Cβ^{2/(2σ+1)}（不是渐近等式），**没有本征函数收敛**。既有记忆里的「τ(a)~a^{2/5}」与此一致（γ=4 ⟹ 2/(γ+1)=2/5）。**摘要自陈 Landau–Pollack–Slepian 起作用、而 1106.0127 全文不提 LPS**——既有 survey §4 的矛盾 M4 本轮复核**仍然成立**。

---

## §4（问 5）诚实评估：这套框架覆盖 ℌ(ϖ) 吗？

**裁决：不覆盖，且不是「补几个条件就能移植」的距离。** 四条互相独立的失配，任一条都足以否决：

1. **发散 vs 饱和（致命）** `[自算，B]`。MS/Sobolev/Widom 全族的小参数极限里顶本征值**收敛到有限值**（M_β→1、λ_α→A₀V₀²），研究对象是**第一修正项**。我方 h(ϖ)=lnϖ+c₀+O(1/lnϖ) **线性发散于 L=lnϖ**。这不是量的差别：他们的模型算子 A=|D|+θ/2 之所以存在，正因为重标后本征函数**集中在单一尺度** α^{−1}；而 repo (B4) 已证我方本征向量 v_n=(2n+1)^{−1/2}g(n/ϖ)/√𝒩、**质量在 1≲n≲ϖ 上按 log 均匀铺开**，𝒩=½lnϖ+0.77。**没有单一尺度可供重标。**
2. **朴素标度极限的算子根本不存在** `[自算，B]`。按 n=ϖs 取极限，L2 块给 L²(0,∞) 上核 **½(st)^{−1/2}/(1+(s−t)²)**；其「秩一部分」对应 u(s)=s^{−1/2}，而 ∫₀^∞ s^{−1}ds 在**两端都对数发散** ⟹ u∉L²(0,∞)，极限算子**无界、无顶本征值**。这与 §1.C 的 ℌ=2uu^T−I−W 中 u∉ℓ² 是同一件事。**既有 survey §3.2「标度极限与 K_β 同型」这句应订正为：核型相同，但极限算子不存在，同型是表面的。**
3. **权与对角项在 MS/Sobolev 里没有对应物**（§1.6、§3.1）。Sobolev 的 V 必须有界且 −V₀+c ≤ V ≤ V₀；(st)^{−1/2} 违反。我方的 −L1 对角相减与 L3 的 Hankel 块 f(n+m+1)/√(…) 在两篇里都零对应。
4. **离散性不是技术细节，而是常数的来源** `[自算，B]`。c₀=γ+ln4−3/2 里的 Euler γ 正是 Σ_{n<N}(2n+1)^{−1} 的离散求和常数（repo §0.A / KAY III eq.23）。n=O(1) 的**边缘区**决定 c₀，而连续极限把这个区抹掉。MS/Sobolev 全程在 ℝ 上连续。

**关于派单问的两项「段 C-II 所需」：**

- **(a) 残差 ‖(I−P)𝒪[ψ_x]v‖ 的一致控制：文献不提供。** MS Thm 1 与 Sobolev Thm 1.2 的本征函数结论都只是**按范数收敛、无速率**（§1.2、§3.1 逐字）。范数收敛推不出任何 ϖ-一致的定量残差界，更推不出 ϖ⁴-归一化下 ϖ⁴E 的上界。
- **(b) 谱隙：文献不提供。** MS 用 Lemma 2（正性改进半群）只得到 λ₁ **单重**，是定性的；Sobolev 同理。两篇都没有 ϖ-一致的隙下界。
- **结论**：即便把失配 1–4 全部克服，**这套框架交付的也只有主阶本征值渐近**，正是段 C-II **不需要**的那一半。

### 4.1 一个精确的等价重述（可直接写进 §1，`[自算，B]`，已与 (B1)(B7) 数值对表）

令 L:=lnϖ、λ(L):=h(e^L)。链式法则给**恒等式**

  **H″(x) = [λ″(L) − 2λ′(L)] / (4x²)　⟹　H 凹 ⟺ λ″(L) < 2λ′(L)。**

用 §1.A.2 的恒等式可把右端算穿。由 (I2)：**ℌ_L := ϖℌ′(ϖ) = 2ℌ − 2x²𝒪[ψ_x]**；再对 L 求导并用 φ_x = ψ_x − x/(k²+x)³：

  **ℌ_LL = 2ℌ_L − 8x²𝒪[φ_x]**　⟹　⟨v,ℌ_LL v⟩ = 2λ′(L) − 8x²F，

于是 λ″ = ⟨v,ℌ_LLv⟩ + 2⟨r_L,R r_L⟩ ≤ 2λ′ − 8x²F + 8x⁴E²/g，**λ″<2λ′ ⟺ x²E²/g < F**——与 §1.A.5 的判据 b⁴E_hi²/g_lo < F_lo **逐字相同**。

**两条可用读数**：(i) 这独立确认 **(C\*\*) 已经是「无相消」形式**，2λ′ 那一项恰好被 (I1) 的恒等式抵掉，没有进一步简化的余地；(ii) 在 L-变量下两端都是 **O(1)**（λ′→1、λ″→0），而 (C\*\*) 两端都 ~ϖ^{−4}。**这只是记账单位的改变，不减少工作量**，但它把「要证什么」讲成一句人话：**证 ϖh′(ϖ) 在 lnϖ 上的变化率小于 2ϖh′(ϖ)**。

### 4.2 两条本轮新识别的候选路线（均 `[自算，B]`，repo §1 的路线 (i)(ii)(iii) 都不是这个）

**路线 (iv)：解析性 + Cauchy 估计——把 ϖ⁴-相对精度换成 O(0.1) 的绝对精度。**
设 λ(L)=L+c+f(L)，f 在复盘 |L−L₀| ≤ w（L₀ 实、≥ lnX）上解析且 |f| ≤ ε。Cauchy 给 |f′| ≤ ε/w、|f″| ≤ 2ε/w²。需要 λ″=f″ < 2λ′=2(1+f′)，**充分条件 ε/w² + ε/w < 1**；取 w=0.5 得 **ε < 1/6 ≈ 0.167**。数值上 f(L)=b₁/L+… 在 ϖ=45（L=3.81）处只有 **0.090**，**已经在裕度内**。
⟹ 段 C-II 所需被换成：**(α) ℌ(x) 在 x 的复邻域内解析且顶本征值保持单重孤立（Kato + Rouché）；(β) 在该邻域上 |h(ϖ)−lnϖ−c| ≤ 0.16 的双侧界（c 任取，不必是真 c₀）。**（β）的量级是 **绝对 0.16**，而 (B5) 的试探向量 (ii) 在实轴上已做到 **6.0e−3**（ϖ=50）——低了 25 倍以上。这是本轮找到的**最大杠杆**。
**必须同时记下的代价**：复 ϖ 处 ℌ 非自伴，"h" 只能定义为解析延拓的孤立单重本征值，需要 Rouché/隙论证保持其孤立；且 Rayleigh 商在复域失效，上界要改走 §1.C/(B8) 的 secular 方程 2⟨u,B(λ)^{−1}u⟩=1 + Neumann 级数。
**为什么不能退回实变量**：对固定 ε 的实变量二阶差分永远推不出 H″<0——中点凹性必须对**所有**步长 a 成立，而 a→0 时 ½ln x 的凹性亏损 ~a²/(4x²)→0，ε 不随之趋零。**⟹ 解析性在这条路上是不可省的。**

**路线 (v)：λ′ 的一致双侧界给「粗尺度单调性」。**
H′(x)=λ′(L)/(2x)。若 0<c ≤ λ′(L) ≤ C 对 L ≥ L₀ 一致，则 **x₂/x₁ > C/c ⟹ H′(x₂) < H′(x₁)**。数值上 λ′=ϖh′ 从 0.9703(ϖ=10) 单调升到 0.99696(ϖ=2e4)，若能证 λ′∈[0.96,1]，则 H′ 在**任何 ϖ₂/ϖ₁ > 1.021 的跨度上自动下降**，只剩相对宽度 2% 的局部窗口要处理。这不闭合凹性，但**把段 C-II 的未知区间压成一族小窗口**，可与路线 (iv) 或证书机器组合。

### 4.3 真正的结构类比不在 Widom 支，在 BCS 支：FHNS arXiv:0704.3564 [A]

Frank–Hainzl–Naboko–Seiringer 处理的正是**对数发散的顶/底本征值**。其 Birman–Schwinger 算子分解（逐字）：

> (7)　"$B_{T}=\lambda\ln(1+\mu/2T)JX+\lambda Y_{T}$"

**Lemma 2**：
> "$Y_{T}$ … is Hilbert-Schmidt, and its Hilbert-Schmidt norm is bounded uniformly in $T$"
> "$\sup_{T>0}\mathrm{tr}[Y_{T}^{\dagger}Y_{T}]<\infty$"

扰动结论：
> "differs from the spectrum of $\alpha JX$ by at most $\Theta(\sqrt{\alpha\lambda})+\Theta(\lambda)$, uniformly in $T$"
> "the lowest eigenvalue of $B_{T}$ is bounded from above and below by $\lambda a_{\mu}(V)\ln(1+\mu/2T)+\Theta(\sqrt{\lambda})$"

**这是「对数发散系数 × 固定有限秩算子 + 一致有界余项」的模板，正是段 C-II 需要的形状** `[自算，B]`：我方 §1.C 的 ℌ=2uu^T−I−W(x) 形式上同构（2uu^T = (lnϖ+const)·P̂，P̂ 秩一投影）。**但关键差别必须写清**：FHNS 的 Y_T **一致 Hilbert–Schmidt**；我方 (B8) 的数值显示 λ_min(𝒪_N[m_ϖ]) 以 ~1/lnN 的速度朝地板 −(1+h) 漂移，**余项不一致有界**，两个对数发散互相抵消。**⟹ 段 C-II 必须先解决的新问题（本轮判定为 well-posed 且尚无文献）**：

  **是否存在 ϖ-无关常数 C 与秩一投影族 P_ϖ（如取归一化的 u1_{n<ϖ}），使 ‖ℌ(ϖ) − h̃(ϖ)P_ϖ‖ ≤ C 一致成立？**

若成立，(β) 的双侧界与 λ′、λ″ 的一致界都能顺着 FHNS 的扰动论证写下来。

### 4.4 问 4 的另外两半

- **顶本征值对尺度参数的二阶导数 / 凹性：文献零处**（Widom 支、Böttcher–Silbermann 支、Basor–Ehrhardt 支均只有行列式或主阶本征值）。最接近的通用机器是 **Stubbe arXiv:2605.24694**（2026-05-23），**Theorem 1.2** 的和规则 (1.17) 与 **Theorem 1.9** 的夹逼 (1.21)：
  > "$0\leq\sum\langle \ddot H u_{j},u_{j}\rangle-\sum\ddot\lambda_{j}\leq d_{m}^{-2}\sum\langle[\dot H,[H,\dot H]]u_{j},u_{j}\rangle$"
  **对我方无增益** `[自算，B]`：取单点谱集时 ⟨[Ḣ,[H,Ḣ]]v,v⟩ = 2⟨(λ₁−ℌ)r,r⟩（r=ℌ_Lv），于是 (1.21) 等价于 ⟨r,(λ₁−ℌ)⁺r⟩ ≤ d^{−2}⟨(λ₁−ℌ)r,r⟩；因我方 λ₁−ℌ 在 v^⊥ 上取值于 [g,g+1] 而 g~lnϖ≫1，该界给 (g+1)‖r‖²/g² > ‖r‖²/g，**比 repo 现用的 Davis–Kahan/隙界更松**。与 (B7) 的观察一致（C_bd/C 已降到 1.03，几乎无余地）。**不推荐。**
- **1/√(nm) 权的 Hilbert 型比较算子**：既有 survey §3.4 的「文献零处」本轮复核**仍成立**。补一条负面信息：**Montgomery–Vaughan 加权 Hilbert 不等式的最优常数本身仍是开问题**（arXiv:2203.14950 摘要：已有方法 "cannot achieve any value below 3.19497"，猜想值 π，"remains open"）[C]。⟹ 即便做逐元比较可行，也拿不到锐常数，与我方需要的 O(0.1) 绝对精度目标不匹配。

---

## §5 矛盾与不确定点

| # | 事项 | 状态 |
|---|---|---|
| N1 | **Widom 1961 / 1963 正文本轮未取到**（ams.org 全站 403，PDF 与 landing page 皆然） | 既有 survey 的「AMS 免费 PDF」**降级为未核实**；§2 全部内容是 Sobolev 的转述 |
| N2 | 既有 survey §3.2「标度极限与 K_β 同型」 | **订正**：核型相同但极限算子不存在（u∉L²(0,∞)），同型是表面的（§4 第 2 条） |
| N3 | M4（Mityagin 0803.3159 称 LPS 起作用 / 1106.0127 全文不提 LPS） | 本轮复核**仍然成立**，两文技术路线确实不同 |
| N4 | 路线 (iv) 的 ε<1/6 阈值 | 是 `[自算]` 的充分条件，依赖 (α) 复邻域解析 + 孤立单重；**(α) 本轮未验证**，是该路线的承重假设 |
| N5 | FHNS 类比中「我方余项不一致有界」 | 依据是 repo (B8) 的 **float64 数值**（λ_min 的 1/lnN 漂移），**不是证明**；有可能换一个 P_ϖ 就一致有界（§4.3 的新问题） |
| N6 | Mityagin 0803.3159 的双侧界指数 β^{2/(2σ+1)} vs 1106.0127 的 β^{2/(γ+1)} | 两处记号不同（σ vs γ），本轮未取 0803.3159 正文核对换算；引用时分开标注 |

---

## §6 建议的后续检索方向

1. **Widom TAMS 100 (1961) 正文**（闭合 N1）：ams.org 被封，改试 JSTOR 机构订阅 / Kesten J. Analyse Math. 10 (1962) 117–138（直接建立在其上，可能复述其定理）/ Böttcher–Silbermann *Introduction to Large Truncated Toeplitz Matrices* (1999) 的历史章。
2. **FHNS 全文 + Hainzl–Seiringer PRB 77, 184517 (2008) + Hainzl "Mathematical aspects of the BCS theory" 讲义**（`ceremade.dauphine.fr/~lewin/data/treilles-2019/Hainzl.pdf`）：抄 §4.3 那个「对数发散秩一 + 一致 HS 余项」论证的完整技术细节，判断能否对 ℌ(ϖ) 复制。**这是本轮认定最值得投入的下一步检索。**
3. **复域延拓的先例**：Kato *Perturbation Theory* II §1（解析族 (A) 型）+ 任何把 Rouché 用于保持 Matsubara 型算子顶本征值孤立性的文献——支撑路线 (iv) 的承重假设 (α)。
4. Mityagin 0803.3159 正文（闭合 N6）与 Adduci 的数值文献（λ₁=0.978…）——优先级低，只影响记号一致性。

---

## §7 裁决（≤12 行）

1. **Mityagin–Sobolev（1106.0127, J. Spectral Theory 1 (2012) 443–460）不覆盖 ℌ(ϖ)，也不是「补条件可移植」的距离。** 四条独立失配：(i) 他们的顶本征值**收敛到有限值**而我方**对数发散**；(ii) 朴素标度极限 n=ϖs 的算子**根本不存在**（u(s)=s^{−1/2}∉L²(0,∞) 两端对数发散）；(iii) 权 (st)^{−1/2} 与对角项 −L1 在 MS 与 Sobolev 1609.02052 里都零对应，且 Sobolev **Condition 1.1** 明确要求权有界（−V₀+c ≤ V ≤ V₀），我方硬失效；(iv) 离散性决定常数（c₀ 里的 Euler γ 来自 Σ(2n+1)^{−1}），连续极限把它抹掉。
2. **即便克服全部失配，它交付的也是段 C-II 不需要的那一半。** Thm 1 的余项只是 **o(β^{2/(γ+1)})**，本征函数只是**按范数收敛、无任何速率**（Sobolev Thm 1.2 同）。**⟹ 对 ‖(I−P)𝒪[ψ_x]v‖ 的 ϖ-一致界与隙下界，这支文献提供零。**
3. **Widom 1961/1963 正文本轮取不到（ams.org 全站 403）**；经 Sobolev 转述，其框架里 λ_max 从下方饱和到符号上确界 A₀、修正 O(L^{−γ})，**自动是增且凹**——但前提是核 ∈L¹、符号有界，**我方线性发散不满足，凹性不可搬**。
4. **顶本征值对尺度参数的二阶导数/凹性：Widom / Landau–Widom / Böttcher–Silbermann / Basor–Ehrhardt 全支零处。** 最接近的 Stubbe arXiv:2605.24694 Thm 1.9 (1.21) 经换算**比 repo 现用的隙界更松**（(g+1)/g² > 1/g），不推荐。1/√(nm) 权的 Hilbert 型比较仍是空白地，且 Montgomery–Vaughan 的最优常数本身未解，拿不到锐常数。
5. **最有价值的结构类比不在 Widom 支，而在 BCS 支：FHNS arXiv:0704.3564 eq.(7) + Lemma 2** 的「**对数发散系数 × 固定秩一算子 + 一致 Hilbert–Schmidt 余项**」，与我方 §1.C 的 ℌ=2uu^T−I−W 同构。**但我方余项按 (B8) 的数值不一致有界**（两个对数发散互相抵消）。⟹ 段 C-II 的**第一个必须新证的命题**：是否存在 ϖ-无关 C 与秩一投影族 P_ϖ 使 **‖ℌ(ϖ) − h̃(ϖ)P_ϖ‖ ≤ C** 一致。
6. **推荐路线（新，repo 的 (i)(ii)(iii) 之外）：解析性 + Cauchy 估计。** 由恒等式 **H″(x)=[λ″(L)−2λ′(L)]/(4x²)**（L=lnϖ），写 λ(L)=L+c+f(L)；若 f 在复盘 |L−L₀|≤w 上解析且 |f|≤ε，则 **ε/w²+ε/w<1 即给 H″<0**；取 w=0.5 得 **ε<1/6≈0.167**，而实轴上 f(L)≈b₁/L 在 ϖ=45 处仅 0.090、(B5) 的试探向量已做到 6e−3。**这把段 C-II 从「ϖ⁴-归一化下的相对精度」换成「绝对精度 0.16」，是本轮找到的最大杠杆。** 代价：需 (α) ℌ(x) 复邻域解析 + 顶本征值保持孤立单重（Kato+Rouché，**本轮未验证，是承重假设**），(β) 复域上的双侧界须改走 §1.C/(B8) 的 secular 方程而非 Rayleigh 商。**实变量退路不存在**：固定 ε 的二阶差分推不出 H″<0（中点凹性需对所有步长成立，而 a→0 时凹性亏损趋零）。
7. **辅助路线 (v)**：若证得 λ′(L)=ϖh′(ϖ)∈[c,C]，则 **x₂/x₁>C/c 时 H′ 自动下降**；以 [0.96,1] 计只剩相对宽度 2% 的局部窗口，可与 (iv) 或证书机器组合。
8. **一条交叉核**：L-变量重述给出的判据 **x²E²/g<F** 与 §1.A.5 的 b⁴E_hi²/g_lo<F_lo **逐字相同**，独立确认 **(C\*\*) 已是无相消形式**，2λ′ 恰被恒等式 (I1) 抵掉，该方向无进一步简化余地。
9. **须要新证的东西（清单）**：(a) §7.5 的一致秩一分解；(b) 复邻域内的解析性与孤立性；(c) 复域双侧界 |h−lnϖ−c|≤0.16；(d) 若走 (v)：λ′ 的一致双侧界。**四项都不依赖 Mityagin–Sobolev 的任何定理**，其技术来源是 FHNS + Kato + repo 自有的 secular 方程。
10. **对外措辞建议**：可写 "the scaled kernel is formally of the difference-kernel type studied by Mityagin–Sobolev [JST 1 (2012) 443] and, with a weight, by Sobolev [arXiv:1609.02052]; however the scaling limit here does not exist as a bounded operator, the top eigenvalue diverges logarithmically rather than saturating, and those results supply no rate for either the eigenvalue remainder or the eigenfunction convergence"——**不得写「我们推广了 Mityagin–Sobolev」。**
