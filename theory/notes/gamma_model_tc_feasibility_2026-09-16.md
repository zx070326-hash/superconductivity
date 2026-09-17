# γ-model 量子临界 Eliashberg 严格 Tc 界：可行性报告

日期：2026-09-16
状态：**调研，未做任何新计算**——所有结论转引自 arXiv:2409.00533（KAY-I，最终版 2025-04-04，已发表 J. Stat. Phys.）与 arXiv:2512.20009（Elezaby–Abanov，v1 2025-12-23 / v2 2026-06-10，已发表 PRB，DOI 10.1103/wpn8-8j2t）。分级标注：`analytic`（论文自证的解析结果）/ `numerical`（论文自认为浮点数值，非证书级）/ `our-gap`（我们尚未验证、需独立核算的断言）。

---

## 0. 一句话结论

**CONDITIONAL-GO**：KAY-I 的 N≤4 下界是精确解析闭式（不是数值），N>4 的下界序列虽有单调收敛的解析证明，但**具体数值（N=400 等）是普通浮点 Maple 计算，不是区间算术 certificate**——这正是我们的技术优势可以嵌入的缺口。Elezaby–Abanov 的 zeta-级数上界是解析证明（Gershgorin + 符号论证），已经比 KAY 的上界紧得多，我们在"上界更紧"这条轴上很难超过他们。真正的新颖性窗口在：(a) 把 KAY 的 N>4 浮点下界升格为**声音的区间算术双侧界**（尤其是 Davis-Kahan 给出本征向量集中度，从而给出比单纯 Rayleigh 商更强的结构性陈述）；(b) **Coulomb μ\* 推广**——KAY-I 明确拒绝将 μ\* 纳入他们的"严格" Eliashberg 理论（称其为"smuggled in… hand-waving"），这是一条他们主动放弃的方向，也是我们 P4 项目已经有整套机器（P4-D Coulomb 全局 Tc certificate）可以直接迁移的方向。

---

## 1. γ-model 算子理论综述

模型定义（KAY-I eq.13）：
$$V_\gamma(\omega_n-\omega_m) := \frac{g^\gamma}{|\omega_n-\omega_m|^\gamma}, \quad n\neq m,\qquad V_\gamma(0):=0.$$

物理对应（KAY-I 文本，2D）：γ=1/3 nematic QCP；γ=1/2 磁性（SDW/AFM）QCP；γ≈7/10 cuprate 自旋液体模型；γ=1 无阻尼传播玻色子配对（这是最常引用的"γ-model"标准情况）；γ=2 恰好是 Holstein/Einstein 声子模型 λ→∞ 极限的特殊情形（对应我们 P4 的 Einstein 声子理论，是两个项目的自然衔接点）。

Eliashberg 线性化 gap 方程（KAY-I eq.4）：
$$\omega_n \Delta_n = \pi T \sum_{m\in\mathbb Z} \lambda_{n,m} \frac{\omega_n\Delta_m - \omega_m\Delta_n}{\sqrt{\omega_m^2+|\Delta_m|^2}}\Big|_{\Delta\to 0}.$$

Tc 由一个**显式构造的紧致自伴算子** 𝔊(γ)（作用在 ℓ² 上，KAY-I §4）的最大本征值 𝔤(γ) 定义：
$$T_c(g,\gamma) = \frac{g}{2\pi}\big[\mathfrak g(\gamma)\big]^{1/\gamma}.$$

结构关键点（与 P4 的对照，见第 4 节）：𝔊(γ) 的矩阵元（例如 KAY-I eq.31 给出的 𝔊^(4)）**目视上逐项非负**——这是 Perron-Frobenius/Krein-Rutman 正锥映射性质的来源（KAY-I p.12："largest eigenvalue 的本征模也是某个非对称算子把 ℓ² 正锥映入自身的谱半径本征模，输入 Krein-Rutman 定理"）。**这与我们 P4 的 Lemma M 路线不同**：我们的 Perron-Frobenius 靠 Bochner 正定性（f_hat(k)≥0，核在频域非负）得到；KAY 靠的是离散矩阵元本身逐项非负，不需要核的 Fourier 变换非负。这解释了为什么主会话的可行性脚本发现"Bochner 正定性对所有 γ∈(0,1] 失败"**不构成障碍**——KAY 根本没用 Bochner 路线，我们也不需要用它来做 γ-model 的秩一/Perron-Frobenius 论证。

---

## 2. KAY-I 与 Elezaby–Abanov 已有结果精确复述

### 2.1 KAY-I 下界（`analytic` for N≤4；`numerical` for N>4）

- N=1 下界是**精确闭式**（KAY-I eq.17）：$T_c^{(1)}=g/2\pi$。
- N=2,3,4 下界是**精确代数闭式**（eq.18–19，23），通过求解 2×2/3×3/4×4 特征多项式（三次/四次方程，用 resolvent cubic/quartic 显式公式，eq.20–30）得到——**这些不是数值解，是解析的代数根式表达式**。这点很重要：KAY 在低 N 端做到了我们做证书通常需要 mpmath.iv 才能达到的严格性，但用的是纯代数手段。
- N>4：$T_c^{(N)}(g,\gamma)$ 序列**单调递增**且收敛到真值 $T_c(g,\gamma)$（eq.15–16，这个单调收敛性本身是 analytic 证明——来自算子紧致性 + Rayleigh-Ritz 变分原理，标准且严格）。但**具体数值**（如 N=400 的图 1、图 2 数据点）是用 **Maple2021 LinearAlgebra[Eigenvalues]** 浮点计算的，KAY-I 原文明说"numerically computed data points"，交叉验证仅到 8–9 位小数一致（两种方法互检 + 与 N=4 闭式比对），**没有任何区间算术/证书级舍入控制**。这正是文中反复强调的方法论坑："多精度一致只测舍入，不测离散化/尾部精度"——KAY 的做法恰好落在这个坑的边缘：他们做了两种独立算法互检（Maple vs [YKA] 方法），但都是浮点，没有声音的误差包络。

- 上界（KAY-I eq.32）：$T_c^*(g,\gamma)$ 也是**闭式解析表达式**（含 Riemann zeta 函数），γ→∞ 时收敛到精确值 g/2π，γ→0 时按 (2/γ)^{1/γ} 发散——这个上界本身是 loose 的（图 1 显示在中等 γ 时上下界差距很大，"off the charts" 描述其在图 2 尺度下完全不在图内）。

### 2.2 Elezaby–Abanov（`analytic`，Gershgorin + 符号论证）

- Hessian 矩阵元显式公式：对角 $H_{nn}=(2n+1)\tau^\gamma+2\sum_{m=1}^n(1/m^\gamma-1/(2n+1)^\gamma)$，非对角 $H_{mn}=-(1/|m-n|^\gamma+1/(m+n+1)^\gamma)$。
- 上界推导链：Gershgorin 圆盘定理作用在相似变换矩阵 $h^{(N)}=O^{-1}H^{(N)}O$（O 是含参数 p 的对角矩阵，保行列式/本征值过零点条件）→ Gershgorin 半径显式求和，N→∞ 时级数对任意 γ>0 收敛 → **p=1/2 的最优性有附录 A 的闭式解析证明**（对所有 m 的符号论证，非数值搜索）→ 代入得
$$\tau_{\rm up}^\gamma=\sum_{n=0}^\infty 4^{-n}\zeta(\gamma+2n+1).$$
- 与 KAY 上界比较：KAY 的上界渐近 $\tau_c\sim(5.2991)^{1/\gamma}$，EA 称其"对有限甚至大 γ 定量上是松的"；EA 上界渐近 γ→∞ 时 →4/3=1.333…，接近数值真值渐近 ~1.1843。图 4（论文自述）展示 EA 界（蓝）vs KAY 界（红）vs 数值数据（黑叉），EA 显著更紧。
- **没有任何区间算术/浮点精度讨论**——EA 论文的严格性完全建立在算子理论（Weyl 定理、Sylvester 惯性定律、紧致性、Rayleigh-Ritz、Perron-Frobenius）上，是纯解析证明，不依赖数值截断的证书化。
- EA 独立复核确认了 KAY 下界（论文自述其独立计算"confirms the lower bounds previously established… by Kiessling and coauthors"）——说明双方在下界上互相印证，上界是 EA 的主要改进点。

### 2.3 两篇论文都没有 Coulomb μ\* 推广

KAY-I 在引言脚注 4 明确批评标准文献（[BR]）把 Coulomb 赝势 μ\*(N) 加到线性化 gap 方程核里的做法："We emphasize that no such Coulomb term can be extracted from the Fröhlich hamiltonian (1), and so it has no place in Eliashberg theory as derived from (1). It is 'smuggled' into the truncated equations using some hand-waving physical reasoning… Note that no non-zero constant can be added to λ(n−m) in Eliashberg theory because the summations over all Matsubara frequencies would diverge, then."

这段话值得逐字记住：**KAY 是刻意、公开地拒绝处理 μ\* 的**，理由是发散性（常数项加到 λ(n-m) 上对所有 Matsubara 频率求和会发散）。这不是他们没想到，是他们判断这条路走不通（至少在他们的框架下）。EA 论文（截至可提取内容）同样未讨论 μ\*。

---

## 3. 我们的 certificate 路线 vs 他们的方法

| 维度 | KAY-I | Elezaby–Abanov | 我们（P4 迁移路线） |
|---|---|---|---|
| 下界机制 | 显式紧算子 𝔊(γ) 最大本征值，N≤4 闭式代数根，N>4 浮点 Rayleigh-Ritz | 未做独立下界（复核 KAY 的） | 有限维矩阵 H_N 的**区间算术 Rayleigh 商**下界，任意 N |
| 上界机制 | zeta-级数闭式，但松 | Gershgorin + 相似变换，闭式且紧 | 若走这条路，可直接照抄 EA 的解析证明（不需要我们自己发明） |
| 数值截断的声音性 | 无（Maple 浮点，仅多方法互检） | 不适用（无截断，是解析级数） | mpmath.iv 区间算术，声音外向舍入 |
| 本征向量集中度 | 未讨论 | 未讨论 | Davis-Kahan sinΘ（我们 P4 c1² certificate 路线的核心工具），可以量化"截断本征向量离真本征向量多近"，而不只是本征值上下界 |
| Coulomb μ\* | 明确拒绝（发散论证） | 未处理 | 我们有 P4-D 全局 Tc(λ,μ\*) certificate 的完整机器（interval-Rayleigh+KT），需要重新核实"KAY 的发散论证"是否真的排除了所有正规化方案 |
| 正定性/单调性前提 | 矩阵元逐项非负（Krein-Rutman 直接用在离散矩阵上） | 承袭 KAY 的正定性结论 | 我们的 Lemma M 走 Bochner 正定性（核 Fourier 变换非负），**已确认对 γ∈(0,1] 全部失败**——这条路对 γ-model 走不通，但也不需要它，因为 KAY 已经证明了矩阵元层面的正定性，我们可以直接借用他们的 Krein-Rutman 论证而非重新发明 |

**核心判断**：我们**不能**在"下界怎么算"和"上界怎么算"这两个维度上原创性超过 KAY+EA——他们的方法（Gershgorin、resolvent cubic/quartic、Krein-Rutman）已经把解析闭式做到很紧。我们的比较优势窄且具体：**把 N>4 的浮点截断升格为声音区间算术界，并用 Davis-Kahan 给出他们完全没有的"截断有多接近真解"的定量控制（而不仅仅是本征值的松紧）**。这是增量贡献，不是突破性贡献——除非我们能在 Coulomb μ\* 或紧凑性边界（γ≤0.5）上做出他们没做的东西。

---

## 4. 技术迁移对照表（P4 → γ-model）

| P4 组件 | γ-model 对应 | 迁移难度 |
|---|---|---|
| Einstein 声子核 V(ω_n-ω_m)=g²/(Ω²+(ω_n-ω_m)²) | V_γ(ω_n-ω_m)=g^γ/\|ω_n-ω_m\|^γ | 核形式不同（洛伦兹 vs 幂律），但都是 n-m 的函数（Toeplitz 结构），矩阵元计算方式相同 |
| Lemma M（Perron-Frobenius via Bochner 正定性） | **不适用**（Bochner 对 γ∈(0,1] 失败）；改用 KAY 的矩阵元逐项非负 + Krein-Rutman | 需要独立核实 𝔊(γ) 矩阵元对所有 γ、所有 (n,m) 逐项非负——KAY 只展示了 N=4 的显式矩阵，未必对任意 N 有通用证明；这是我们要做的第一块地基工作 |
| Rayleigh 商区间算术下界（H_N 截断） | 直接照搬，替换核函数为 V_γ | 低难度，核是初等幂函数（无 arctan/atan2 那类 iv 陷阱），比 Einstein 核（含 atan2）更简单 |
| Davis-Kahan sinΘ → c1² 集中度 | 直接照搬 | 低难度，前提是先建立 γ-model 的秩一主导性（可行性脚本已报告 c1²~0.70-0.85，需独立复核） |
| 尾部/Schur 检验（行和列和衰减） | V_γ 的行和 ~ Σ 1/k^γ，γ≤1 时调和级数发散，需要不同的尾部控制手段（Toeplitz 但非绝对可和） | **中高难度**——这是可行性脚本报告"γ≤0.5 需要额外论证"的紧凑性问题的根源；HS 紧凑性要求 Σ 1/k^{2γ} <∞ ⟺ γ>0.5，γ≤0.5 时算子甚至不是 Hilbert-Schmidt，需要换紧凑性判据（例如紧算子的极限/紧化，或直接在原始 ℓ² 算子层面而非 HS 范数层面论证，正是 KAY 用紧致算子直接构造而非诉诸 HS 范数的原因） |
| Gate A/证书脚本约定（独立新文件扩展，探针不提交） | 直接照搬仓库机械约定 | 无迁移成本 |
| P4-D Coulomb 全局 Tc(λ,μ\*) certificate | γ-model + μ\*：需要先解决 KAY 声称的发散问题——μ\* 作为 n=m 处的一个常数项加进 λ(n-m) 会让对所有 Matsubara 频率的求和发散 | **高难度、高价值**——如果这个发散是可以通过某种正规化（例如截断到能带宽 + 重整化，类似 BCS 的 Ω_cut 处理）声音地控制的，这是唯一一个 KAY/EA 都没做、且新颖性明确的方向 |

---

## 5. 发表策略和 impact 评估

**竞争态势判断**：KAY-I 和 EA 都已经发表在同行评审期刊（J. Stat. Phys. / PRB），都自称"rigorous"，且这个"rigorous"的含义与我们仓库对"analytic"的定义相符（不是浮点数值证据，是解析或代数闭式证明），只是他们在有限 N 截断的具体数值（KAY 的 N=400 数据点）上退回到普通浮点。

- **选项 A（sound certificate 双侧界）**：价值有限。KAY 的 N≤4 已经是精确代数解，EA 的上界已经是解析闭式证明，二者组合起来在大多数 γ 范围内已经把"严格"这个词占满了。我们把 N>4 的浮点升格为区间算术，只是**堵上一个已经很窄的缝隙**（N=400 时 8-9 位小数一致，区间算术顶多把这个不确定性从"浮点误差量级不明"变成"证书宽度已知且极小"）——这不太可能构成 PRL 级别的新颖性，因为读者会问"KAY 的数值已经和两种独立方法吻合到 9 位，你证书化它有什么本质区别"。

- **选项 B（Coulomb μ\* 推广）**：**价值最高**。这是 KAY 明确声明放弃、EA 未触及的方向。如果能证明：(i) μ\* 项可以通过某种声音的正规化方案（保持 Eliashberg 理论的可推导性，不是"smuggled in"）引入 γ-model 而不发散；(ii) 给出 Tc(g,γ,μ\*) 的严格（哪怕只是 certificate 级）双侧界；这将是**第一个**处理 γ-model + Coulomb 排斥的严格结果，直接对标我们 P4-D 已完成的 Einstein 声子 + Coulomb 全局界工作，构成清晰的"论文二部曲"叙事（P4 声子 QCP-free vs P4-γ 量子临界 with Coulomb）。但**前提风险很大**：KAY 给出的发散论证（"summations over all Matsubara frequencies would diverge"）需要先被证明是可规避的，而不是真的不可规避——这是一个需要独立解析核实的开放问题，不能假设它一定有解。

- **选项 C（全 γ 相图，统一 SDW+nematic+FM QCP）**：中等价值，更像是一篇很扎实的后续/综述型论文而非突破。KAY-I 已经画出了 γ∈(0,3.2] 的完整 T_c(γ) 数值曲线（图1/图2），我们如果只是重新证书化这条曲线的不同区段，新颖性增量小。若能把 γ≤0.5（非 HS 紧凑区）也纳入严格框架，这是一块 KAY/EA 都语焉不详的空白（KAY 的下界公式在 γ→0 时如何定义 Tc、是否仍良定义未被两篇论文明确讨论紧凑性失效后的替代论证），有一定价值但物理相关性（真实材料对应的 γ 值都 ≥1/3）有限。

**总体判断**：单纯做"更紧的证书"（选项 A）**不足以发 PRL**——KAY/EA 已经把解析闭式做到很紧，我们的边际贡献主要是方法论上的（"我们的界是声音的区间算术，他们的数值截断不是"），这更适合作为方法论笔记或与 P4 论文合并的附录，而不是独立的 PRL。真正值得投入的是**选项 B（Coulomb 推广）**，但这需要先花 1-2 天做一个**纯解析的可行性判定**：KAY 的发散论证到底在多大范围内成立，是否存在正规化窗口。这个判定本身不需要大量证书基础设施，是纸笔/mpmath 探针级别的工作，应该作为 Phase 1 的第一步而不是最后一步。

---

## 6. 路线图

### Phase 0（先于 Phase 1，1-2 天，决策性）
- 独立核实 KAY 的 Coulomb 发散论证：写一个最小 mpmath 探针，检验 λ(n-m)+μ\*（常数）对 Matsubara 求和是否真的发散，还是仅在无穷截断极限下发散、有限截断/物理能带宽正规化下良定义（类比 BCS 的 Ω_cut，KAY 自己在引言就批评过这种"cutoff 走后门"的做法，所以要小心不要重蹈同样的批评）。
- 若发散不可规避 → 降级选项 B 为 NO-GO，转向选项 A/C 的降级目标（作为 P4 论文的方法论附录，而非独立论文）。
- 若存在声音的正规化方案 → 进入 Phase 1，全力做选项 B。

### Phase 1（1 周，若 Phase 0 判定 GO）
- 复核可行性脚本的原始数据（c1²~0.70-0.85 秩一主导性、γ=1 紧凑性、特征值收敛）——**从零独立重算**，按仓库纪律不采信已有"已完成"结论。
- 建立 𝔊(γ) 矩阵元逐项非负的通用（任意 N，任意 γ）证明或反例，取代 Lemma M（Bochner 路线对 γ∈(0,1] 已确认失败）。
- 实现 H_N 区间算术 Rayleigh 商下界脚本（复用 P4 证书框架，核函数替换为 V_γ+μ\* 正规化项）。
- milestone：γ=1（SDW QCP，最重要物理情况）的 certificate 级双侧界，μ\*=0 情况先做通，作为 sanity check against KAY 的 N=4 精确闭式（数值应严格落在 KAY 的 T_c^(4) 与 T_c^* 之间）。

### Phase 2（2-3 周）
- 引入 μ\*>0，仿照 P4-D 的 interval-Rayleigh+KT 路线，做出 Tc(g,γ,μ\*) 的全局 certificate 界。
- Davis-Kahan 集中度分析，量化截断本征向量偏差（KAY/EA 都没做的定量陈述）。
- γ≤0.5 紧凑性问题的处理（若时间允许；物理相关性有限，可降低优先级）。
- 论文草稿：聚焦"第一个 sound Coulomb-inclusive γ-model Tc certificate"，与 KAY-I 的发散论证正面对话（引用其脚注 4 的批评，说明我们的正规化方案为何不落入同样的"smuggled cutoff"陷阱）。

---

## 7. 结论

**CONDITIONAL-GO**，条件是 Phase 0 的 Coulomb 发散可规避性判定为正。理由：
1. KAY-I + Elezaby-Abanov 已经把"γ-model 严格 Tc 界"这个问题的核心（无 Coulomb、任意 γ 的解析上下界）做到很紧、很严格，我们在这条轴上很难产出 PRL 级别的增量贡献——不要低估他们的解析深度（N≤4 精确代数解、Gershgorin 闭式上界都是真正的 analytic 结果，不是可以轻易超越的浮点数值）。
2. 我们真正独特的技术资产（区间算术 certificate 流水线、Davis-Kahan 集中度、P4-D 的 Coulomb 全局界机器）价值最大的落点是 KAY 明确声明放弃的 Coulomb μ\* 方向，而不是重复他们已经做好的下界/上界。
3. 在投入 Phase 1/2 的证书基础设施之前，必须先用最小成本（mpmath 探针，不是完整证书）独立核实 KAY 的发散论证是否真的排除了所有正规化方案——这是决定 GO/NO-GO 的关键判定，且目前**尚未验证**，是本报告最大的 `our-gap`。
