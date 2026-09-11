# P4 · KAY Conjecture 2（全 λ 严格 Tc 上界）—— 化约为 Einstein 逐点不等式 + 三段证明计划

建立：2026-09-11　状态：v0.2（主会话立项 + 数值探针 + §1 严格化约 + researcher KAY 逐字核并入；证书未落地、未审）
上游：`P4_interchange_uniformity.md` v0.6（换序区间定理 [0.2,60]、认证谱隙、ϖ-导数包络）、`P4_rigorous_tc_bound.md` §P4-C（𝔤(2) 严格双侧包围 [1.31714,1.31820]、Cholesky-deflation α 证书、`prop5_tail`）、`p4_remainder_interchange.py`（二阶本征值余项 (*)）。
目标：把 KAY II Conjecture 2（= KAY III Conjecture 1 的一般-P 版本）由猜想升为定理（计算机辅助）。

## 0. 陈述（KAY II p.14 eq.31，逐字见 `surveys/P4_rigorous_tc_bound_scoping_2026-09-09.md` §82）

对所有 P∈𝒫、λ>0：T_c(λ,P) < T̃_c(λ,P) := (1/2π)·√(𝔤(2)⟨ω²⟩λ)，𝔤(2)=ρ(𝔊(2))，⟨ω²⟩=∫ω²P(dω)。

**★ researcher `conj2-kay-check` 逐字核（2026-09-11，`surveys/P4_KAY_conjecture2_reduction_2026-09-11.md`，A 级）**：KAY II Conj. 2 与 KAY III Conj. 1（p.12 eq.28，Einstein：T̃_c=(Ω/2π)√(𝔤(2)λ)）均以 "There is a critical temperature … which for all … λ>0 is bounded above by" 起句，**含存在性**（KAY 只证 λ≥λ_*<0.7670 时存在）。**本文件只证「界」的那一半**——对所有 T≥T̃_c 正常态线性稳定（λ𝔨(P,T)<1），故临界面与线 ℒ(λ,P) 的一切交点满足 T<T̃_c；存在性/唯一性（⟺ ϖ↦𝔥 在全 ℝ₊ 单调）仍开放。KAY 自身无 sup/首次穿越的 T_c 定义（只在临界面为图时由反函数定义，Thm 5 eq.17；pp.7–8 明写可能多次穿越），**不替他们定义**。KAY III 的 Conj. 1 是 T_c 形式而非逐点 𝔥(ϖ)≤𝔤(2)ϖ²，但界的那一半与逐点严格不等式等价（survey §1.3，B 级一行推导）。**「一般 P ⇐ Einstein」的 Jensen 化约在 KAY I–III 零处**（最接近：KAY II Remark 2 p.10「𝔨 不是 𝔥 的 P-平均」只说不等不说方向；Prop 7 p.31 eq.83 ρ(𝔎_j)≤ρ(𝔊_j(2))∫ϖ²P 只对正部分 j=2,3，对含 −𝔎₁ 的全 𝔎 不适用）。**竞争**：截至 2026-09-11 无人声称证明；唯一新严格上界 Elezaby–Abanov arXiv:2512.20009v2（PRB 2026）只做 γ-model，Gershgorin 闭式给 𝔤(2)≤≈1.545（T_c 因子 1.083），无 Einstein/猜想。

**对外措辞（researcher 裁决 (c)，采纳）**："We prove the T_c bound asserted in Conjecture 2 of [KAY II] (equivalently Conjecture 1 of [KAY III]): for every P∈𝒫 and λ>0 the normal state is linearly stable for all T ≥ T̃_c(λ,P); in particular T_c(λ,P) < T̃_c(λ,P) whenever T_c is defined. The reduction to the Einstein case is an elementary convexity argument not stated in [KAY II, III]; the existence/uniqueness of T_c for all λ>0 (monotonicity of ϖ↦𝔥(ϖ) on ℝ₊) remains open."

## 1. 化约（主会话严格写法，v0.1；待 researcher 核 KAY III 是否已写出）

**设定**。𝒫 = ℝ₊ 上的概率测度，⟨ω²⟩=∫ω²P(dω)<∞。记 ϖ=ω/2πT，KAY II 的算子按矩阵元定义为 P-平均（eq.39 后，survey §51 逐字）：

  𝔎(P,T)_{nm} = ∫ ℌ(ϖ)_{nm} P(dω)，　ℌ(ϖ) = −ℌ₁(ϖ)+ℌ₂(ϖ)+ℌ₃(ϖ)，矩阵元由 ⟦k⟧(ϖ)=ϖ²/(k²+ϖ²) 组成。

**引理 A（Jensen，逐元积分 ⟹ 二次型积分 ⟹ 顶本征值凸性）**。对任意有限支撑单位向量 u（有限支撑在 ℓ² 稠密，且 λ₁ 是 sup over 有限支撑的 Rayleigh 商，因 𝔎 有界自伴）：
  ⟨u,𝔎u⟩ = Σ_{nm} ū_n u_m ∫ℌ(ϖ)_{nm}P(dω) = ∫ ⟨u,ℌ(ϖ)u⟩ P(dω) ≤ ∫ h(ϖ) P(dω)，
有限和与积分交换无条件成立（有限和）；每个 ⟨u,ℌ(ϖ)u⟩ ≤ h(ϖ):=λ₁(ℌ(ϖ)) 是 Rayleigh 商定义。对 u 取 sup：
  𝔨(P,T) := λ₁(𝔎(P,T)) ≤ ∫ h(ϖ) P(dω)。　　(A)
可积性：由段 III 的解析界 h(ϖ) ≤ ‖ℌ(ϖ)‖ ≤ Cπϖ/2 及 ∫ϖP ≤ (∫ϖ²P)^{1/2} < ∞。

**引理 B（Einstein 逐点不等式，本文件主体）**。h(ϖ) < 𝔤(2)·ϖ² 对所有 ϖ>0（h(0)=0 时取等）。

**定理（Conjecture 2 之「界」，强形式、不依赖 T_c 定义）**。设引理 B 成立。则对所有 P∈𝒫（P≠δ₀）、λ>0 及所有 T ≥ T̃_c := (1/2π)√(𝔤(2)⟨ω²⟩λ)：λ𝔨(P,T) < 1，即正常态线性稳定；因此临界面 𝒮_c 与线 ℒ(λ,P) 的一切交点满足 T < T̃_c，特别地只要 T_c(λ,P) 有定义（任何合理定义：sup、首次穿越、反函数）就有 T_c < T̃_c。
**证**。由 (A) 与引理 B，对任意 T>0：
  𝔨(P,T) ≤ ∫h(ϖ)P < 𝔤(2)∫ϖ²P(dω) = 𝔤(2)⟨ω²⟩/(2πT)²，
严格性因 P 在 ω>0 上有正质量（P≠δ₀）且 h(ϖ)<𝔤(2)ϖ² 在 ϖ>0 逐点严格。于是当 T ≥ T̃_c 时 𝔨(P,T) < 𝔤(2)⟨ω²⟩/(2πT̃_c)² = 1/λ，即 [T̃_c,∞) 上 λ𝔨<1 严格。交点 T<T̃_c 由此直接得出；若还想要「T_c ≤ T̃_c−ε」形式的定量严格性，用 T↦𝔨(P,T) 的连续性（矩阵元关于 T 连续且一致有界 ⟹ 算子范数连续 ⟹ λ₁ 连续，Weyl）与 T̃_c 处的严格不等式。∎

**与 KAY 的关系（A 级）**：此化约 elementary but, to our knowledge, not stated in [KAY II, III]；KAY 最接近的两处是 Remark 2（p.10）与 Prop 7（p.31 eq.83，只对正部分 𝔎₂,𝔎₃ 用 Perron），论文应引以示区别。

**注**。(i) 证明**不需要** ϖ↦h(ϖ) 单调；单调性（§3）是独立的加分项，它对应 Conjecture 的**存在性/唯一性**半句（全 ℝ₊ 单调 ⟹ 每条线恰一次穿越），KAY III p.28 自陈 "would require a different strategy"。(ii) 不需要 P 紧支，只需 ⟨ω²⟩<∞。(iii) 引理 B 恰是 KAY III 对 Einstein 声子 P=δ_{ω_E} 的 Conjecture 1（待 researcher 逐字核对陈述与严格性）。**⟹ 一般-P 猜想被 Einstein 特例蕴含；本文件只需证引理 B。**

## 2. 引理 B 的三段证明计划

数值探针（`notes/_scratch_ir/conj2_probe.py`、`conj2_probe2.py`，float64 证据、非证书）：r(ϖ):=h(ϖ)/ϖ² 在 [0.01,3] 上严格单调下降；ϖ→0 时 r→𝔤(2)，r = 𝔤(2) − 1.193·ϖ² + O(ϖ⁴)（二阶系数 = ⟨v₀,𝔊(4)v₀⟩=1.19317，与网格逐位符）；r(0.3)=0.926𝔤(2)、r(0.7)=0.701𝔤(2)、r(1)=0.543𝔤(2)、r(2.4)=0.199𝔤(2)。极值只在 ϖ=0 从上方逼近。

**恒等式（级数；一般项与算子范数收敛是我方要证的——KAY II eq.52–54 / III eq.57 只写到 ϖ²𝔊(2)−ϖ⁴𝔊(4)±⋯，"real analytic in ϖ²" 是断言）**：⟦k⟧(ϖ)/ϖ² = 1/(k²+ϖ²) = Σ_{j≥0}(−1)^j ϖ^{2j}/k^{2j+2}（|ϖ|<k，即对全部 k≥1 当 ϖ<1），逐元代入 ℌ 的定义（与 𝔊(γ) 定义同型，`H_matrix` vs `G_matrix`）得

  ℌ(ϖ)/ϖ² = Σ_{j≥0} (−ϖ²)^j 𝔊(2j+2)，　ϖ<1，算子范数收敛（‖𝔊(γ)‖ 由 Prop5/6 型 ζ-界一致有界，或由 |𝔊(γ)| ≤ |𝔊(2)| 逐元 dominate 的 Schur 界）。

- **段 I（小 ϖ，解析微扰，ϖ∈(0,ϖ₀]，目标 ϖ₀≈0.5–0.7）**：A(ϖ)=𝔊(2) − ϖ²𝔊(4) + R(ϖ)，‖R‖ ≤ Σ_{j≥2}ϖ^{2j}‖𝔊(2j+2)‖。二阶本征值余项 (*)（`p4_remainder_interchange.py` 行 40）给 λ₁(𝔊(2)−ϖ²𝔊(4)) ≤ 𝔤(2) − ϖ²⟨v₀,𝔊(4)v₀⟩ + ϖ⁴‖𝔊(4)‖²/(g−2ϖ²‖𝔊(4)‖)，g = 𝔤(2)−λ₂(𝔊(2)) 由 P4-C 的 α 证书（λ₂≤2.0e-2）严格；⟨v₀,𝔊(4)v₀⟩ 的严格下界由数值试探向量 + Davis–Kahan 角度（认证谱隙 ≈1.30、残差小）给（有限 N 的正性 ⟨𝔊^(N)(4)⟩₂>0 是 KAY II Prop 4 p.24 eq.56 / III Prop 5 eq.59，A 级；N=∞ 的正性与定量值由我方给）。**附带贡献**：KAY III §3.6.2「大 λ 已证」依赖 Thm 7 的 N=∞ 展开断言（无谱隙、无换序证明），段 I 闭合即补证之。要证：−c₂ϖ² + Cϖ⁴ + ‖R‖ < 0 ⟺ ϖ² < c₂/(C+…)，粗算 ϖ₀≈0.7。**全部量须 mpmath、进界的无 float64**。
- **段 II（中 ϖ，区间证书，[ϖ₀, ϖ₁]，ϖ₁≈6）**：逐 cell 上界 h(ϖ) ≤ λ₁(A_head)+‖B‖+‖D‖ + L(b−a)/2（λ₁(A_head)≤c 由 Cholesky cI−A SPD，**不需 deflation**；‖B‖,‖D‖ 分段包络；L 加权第 0 行 `schur_full_seg_w0`），与 𝔤(2)_lo·a²（cell 左端）比较。余量 ≥6%（ϖ₀=0.7）到 80%（2.4），网格可粗，Na=70/140。
- **段 III（大 ϖ，纯解析，ϖ≥ϖ₁）**：‖ℌ(ϖ)‖ ≤ Schur 行和界 ≤ C·S1(ϖ)，S1=Σ_k ϖ²/(k²+ϖ²) ≤ πϖ/2 ⟹ 线性于 ϖ，而 𝔤(2)ϖ² 二次 ⟹ 对 ϖ ≥ ϖ₁ 闭合（C=2 时 ϖ₁=π/𝔤(2)=2.39；用现成 row_env 常数定 ϖ₁）。

三段接口：ϖ₀ 由段 I 的解析界定，段 II 从 ϖ₀ 起；段 II 上端 ≥ 段 III 的 ϖ₁。若段 I 只闭到 ϖ₀<0.5，段 II 向下延伸（ϖ=0.3 处余量仍 7.4%），但 ϖ→0 端 cell 余量 ~1.19ϖ⁴ 趋零、区间证书必失效——**段 I 的解析微扰是唯一真卡点**。

## 3. 顺带：𝔥(ϖ) 单调性（KAY III Thm 5 p.8：所有有限 N 在 [0,ϖ_*)、ϖ_*≥√2；N=∞ 经 KAY II p.34 收敛散文 + 解析性排平台；N=3,4 在 λ<λ_* 区间靠 "visual inspection" p.8/p.26；全 ℝ₊ 自陈 "would require a different strategy" p.28；𝔥/ϖ² 的单调性 KAY 零处）

h 在认证谱隙区间 [0.2,60] 上单一本征值 ⟹ Hellmann–Feynman h′(ϖ)=⟨v,ℌ′(ϖ)v⟩。逐 cell 严格下界：⟨v̂,ℌ′v̂⟩（mpmath 精确二次型，v̂ 有限支撑）− Davis–Kahan 角度项 − 尾包络 − Lipschitz（|d/dϖ⟨v,ℌ′v⟩| ≤ ‖ℌ″‖+2‖ℌ′‖²/gap，包络阶 1、2 齐备）。大 ϖ 端 h′ 小（20→30 约 0.04/单位）而 slack 不小，可能只闭到有限上端，如实报。

## 4. 分派（2026-09-11）

- researcher `conj2-kay-check`：KAY III 逐字——Conjecture 1 的精确陈述（是否即 h(ϖ)≤𝔤(2)ϖ²、严格性、Tc 定义）；KAY II/III 是否已写出 §1 的 Jensen 化约或任何"一般 P ⇐ Einstein"陈述；单调性已证范围；是否出现 𝔊(4)/级数展开；KAY III 对 Conjecture 1 的数值证据范围。产出 `surveys/P4_KAY_conjecture2_reduction_2026-09-11.md`。
- implementer `conj2-small-w`：段 I，脚本 `p4_conj2_small_w.py`（含 ‖𝔊(γ)‖ 严格界、二阶余项、Davis–Kahan、ϖ₀ 的认证值）。
- implementer `conj2-interval`：段 II + 段 III + §3 单调性，脚本 `p4_conj2_interval.py`。
- 主会话：§1 化约的严格写法、接口检查、采信复核。
- reviewer：三件合体后对抗核。

## 5. 诚实边界（前置）

- 结论是计算机辅助定理（mpmath 证书），非纯解析证明。
- ~~若 KAY III 已写出 Jensen 化约……~~ researcher 核毕：化约 KAY 零处，可写 "not stated in [KAY II, III]"；但**不写「首证 Conjecture 2」**——猜想含存在性半句（全 ℝ₊ 单调），我方不证；措辞按 §0 的 (c)。
- KAY III §3.6.2「大 λ 已证」实际依赖未证的 N=∞ 展开；段 I 顺带补证，措辞 "we also close a gap in the large-λ argument of [KAY III, §3.6.2]"。
- 引理 B 只在 ϖ>0 逐点，ϖ→0 端由解析微扰覆盖、非证书网格；段 I 的 ϖ₀ 若闭不到与段 II 相接则定理不成立、如实报。
