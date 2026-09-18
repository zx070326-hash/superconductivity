# Impact ≥ 2 物理方向调研：技术栈迁移候选

日期：2026-09-16
严格性：调研报告（非证明）
状态：**未提交**，报回主会话

---

## 0. 我们的技术栈总结

| 工具 | 成熟度 | 标志性成果 |
|---|---|---|
| 紧自伴算子的谱不等式（变分原理、Weyl 不等式） | 核心 | 定理 1/2 全链、引理 M |
| 秩一 / 低秩扰动稳定性 | 核心 | 𝔎_c = 2\|u⟩⟨u\|、Coulomb 推广 |
| Davis–Kahan sin Θ 本征向量扰动 | 核心 | c₁² 集中界、D5 路线 |
| 区间算术 certificate（mpmath.iv / Fraction） | 核心 | 1668-check 换序定理、全局 Tc bracket |
| Schur 检验（行和 → 算子范数界） | 成熟 | λ₂ 界、rearrangement identity |
| Perron–Frobenius / Bochner 正定性 | 成熟 | 核单调性引理 M、锥不变性 |
| 非线性本征值问题的凹凸分析 | 成熟 | 引理 C（H 凹性）、q-单调化 |
| Bhatia–Davis 方差界 | 新 | 定理 F/G（Gap 1 闭合） |

---

## 1. 五个方向的深度调研

### 1a. 量子多体物理的严格界

**2D 面积律（Area Law）**
- 现状：1D 已证（Hastings 2007）。2D **无挫折**已证（Anshu-Arad-Gosset 2021, [STOC 2022](https://dl.acm.org/doi/10.1145/3519935.3519962)）。**一般 2D 有隙哈密顿量的面积律仍 OPEN** — 数学物理最重要的开放问题之一。
- 关键障碍：需要构造近似基态投影子（AGSP），本质上是张量网络 + 信息论工具，**不是**算子不等式。
- 技术匹配度：**Low**。我们的谱不等式工具无法直接构造 AGSP；入场成本 ≥ 6 个月文献。
- Impact：~1（Nature/Science），但 feasibility 极低。
- **结论：不适合。**

**MBL（多体局域化）**
- 现状：是否在热力学极限存在**仍有严重争议**（[2026 综述](https://iopscience.iop.org/article/10.1088/1361-6633/ad9756)）。数值表明临界系统尺寸 ~L=50，远超精确对角化能力。核心问题是雪崩机制（Griffiths 稀有区域是否摧毁局域化）。
- 技术匹配度：**Low**。MBL 的关键是多体本征态的统计性质，不是单个算子的谱结构。需要完全不同的工具（RSRG-X、量子太阳模型、精确对角化统计分析）。
- **结论：不适合。**

**谱隙 certificate（自旋系统）— ★★★ 重点**
- 现状：Rai-Kull-Emonts-Tura-Schuch-Baccari（[Quantum, 2026](https://quantum-journal.org/papers/q-2026-04-13-2065)）建立了**无挫折自旋系统谱隙的 SDP certificate 层级**：
  - 核心算子不等式：$H^2 - \delta H \succeq 0 \Leftrightarrow$ 谱隙 $\geq \delta$
  - 构造 n-体正算子的平移不变分解（LTI-SDP）
  - 已有界严格优于 Knabe 界和 Gosset-Mozgunov 界
  - **AKLT 链隙 certified 到 0.34976**（真值 ~0.350），极其接近
- **Haldane 猜想**（spin-1 Heisenberg 反铁磁链有隙）：**仍完全 OPEN**（[emergentmind 页面](https://www.emergentmind.com/open-problems/haldane-conjecture-integer-spin-heisenberg-chain-gap)）。AKLT 点已证（2D 也已证，[Pomata-Wei 2019](https://arxiv.org/html/1911.01410v2)，[Lemm-Sandvik-Wang PRL 2020](https://link.aps.org/doi/10.1103/PhysRevLett.124.177204)），但**物理 Heisenberg 点无任何严格结果**。
  - 困难：Heisenberg 链不是 frustration-free，所有已有数学工具（detectability lemma、AGSP）都要求 frustration-freeness。
- 2D AKLT 谱隙：已在蜂巢、4.8.8、3.12.12、4.6.12 等格上证明。**正方格 AKLT 隙仍 OPEN。**
- **技术匹配度：HIGH** — 见下方 Top 3 详细分析。

### 1b. 量子信息 / 量子计算

**量子纠错阈值**
- 现状：Google 2024 Nature 实验达到 surface code 阈值以下（[Nature 2024](https://www.nature.com/articles/s41586-024-08449-y)）。理论界主要关注解码器算法和噪声模型。
- 技术匹配度：**Low**。编码理论问题，非算子谱问题。
- **结论：不适合。**

**量子 Fisher 信息**
- 现状：QFI 与谱隙有已知联系，但关系已被较好理解。最近工作偏向机器学习预测（[arXiv:2606.02986](https://arxiv.org/html/2606.02986v1)）和多参数计量学（conic programming）。
- 技术匹配度：**Medium**。QFI 是密度矩阵的谱函数，我们的工具原则上可用，但领域已相当成熟，增量贡献空间小。
- **结论：不建议作为主方向。**

**Barren plateau**
- 现状：Barren plateau 已被广泛研究。最新进展是"provable absence of barren plateau ⟹ classical simulability"（[PMC 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12378457)）— 即 barren plateau 的*缺失*本身意味着可经典模拟，这对变分量子算法是坏消息。
- 技术匹配度：**Low-Medium**。核心是随机电路的集中不等式，不是确定性算子谱。
- **结论：不适合。**

### 1c. 冷原子 / 超流体

**Bertsch 参数严格界 — ★★ 值得考虑**
- 现状：幺正费米气体 Bertsch 参数 $\xi = E/E_{\rm FG} \approx 0.376$（QMC 高精度）。严格界方面：
  - **上界**：Seiringer 等人通过试探波函数给出（但远大于 QMC 值）
  - **下界**：Lieb-Thirring 型不等式给出（但远小于 QMC 值）
  - 严格界与真值之间有巨大 gap
- 技术匹配度：**High**。BCS 配对方程就是非线性本征值问题，与 Eliashberg 同构。秩一/低秩配对核结构、区间算术 certificate、变分谱界 — 全部直接可用。
- 但存在本质困难：Bertsch 参数是**多体**基态能量，不仅仅是 BCS 方程的最大本征值。从 BCS 到多体真值需要 beyond-mean-field 控制。
- **结论：有价值，但需要认真评估 mean-field → many-body 的 gap。**

**BCS 配对方程的严格结果**
- Seiringer 等近期工作：[Enhanced superconductivity at a corner](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/enhanced-superconductivity-at-a-corner-for-the-linear-bcs-equation/39F8C08749921E9C64E9EF8C5E5076A2) — 角几何增强超导的线性 BCS 方程严格结果。
- Frank-Hainzl-Seiringer 系列论文建立了 BCS 泛函的数学理论。
- 我们的 P4 工作与他们在**方法论上高度对标**，但聚焦不同（我们做 Eliashberg 强耦合，他们做 BCS 弱耦合几何效应）。
- **交叉机会：把我们的强耦合 certificate 技术应用到他们的框架。**

### 1d. 统计力学严格结果

**3D Ising 临界指数**
- 现状：conformal bootstrap 已给出极高精度数值界（Kos-Poland-Simmons-Duffin，[Caltech 2016](https://authors.library.caltech.edu/records/65v19-7hc54)）。但这些界不是传统意义上的数学证明 — 它们依赖于正定性条件的 SDP 截断。
- 技术匹配度：**Medium**。conformal bootstrap 的核心也是正定性（算子乘积展开的正性），但具体技术（SDP over 交叉方程）与我们的工具不同。入场成本高。
- **结论：技术 gap 太大，不建议。**

### 1e. 凝聚态新普适类

**Planckian 界 — ★★ 有趣但技术不匹配**
- 现状：两个重要新结果：
  1. **Cao (PRL 2026)**：[Planckian bound on quantum dynamical entropy](https://journals.aps.org/prl/abstract/10.1103/wd91-v58h) — 量子动力学熵增长率 $\leq cT$。证明用集簇分解 → emergent Gaussianity → 精确速率公式 → 数值优化界。
  2. **Chowdhury (arXiv 2512.09979)**：[Planckian bounds via spectral moments](https://arxiv.org/html/2512.09979v1) — 光电导谱矩比满足普适 Planckian 上界。证明核心是**热带中的解析性 + Schwarz-Pick 引理**。最紧常数 $C = 16/\pi$。
- 两者都是**复分析 / KMS 条件**驱动的界，而非算子不等式。
- 技术匹配度：**Low-Medium**。虽然谱矩方法表面上涉及算子，但核心技巧是复分析。我们要改行。
- **结论：方向极热（PRL 级），但不是我们的武器能打的仗。**

**SYK 模型严格界**
- 现状：两个 2026 年重大结果：
  1. [Sharp bounds on SYK ground state energy](https://arxiv.org/abs/2607.27185)（Basu-Kothari-Midha）— 用 Johnson scheme 代数方法。
  2. [Free energy and spectral edge of SYK](https://www.alphaxiv.org/zh/abs/2607.18998)（Yukun He）— Schwinger-Dyson 热力学严格推导。
- 技术匹配度：**Low**。SYK 的随机矩阵 / 组合代数工具与我们的确定性算子工具完全不同。
- **结论：不适合。**

**奇异金属 / 非费米液体**
- 没有适合严格算子界的切入点。输运理论的严格化需要 Kubo 公式 + 解析延拓，不是我们的专长。
- **结论：不适合。**

---

## 2. 量子临界 Eliashberg 的 Tc 界（与我们直接竞争的方向）

**必须讨论的背景**：在调研中发现了我们的**直接竞争者**：

| 组 | 论文 | 结果 | 发表 |
|---|---|---|---|
| Kiessling-Altshuler-Yuzbashyan | [Part I: γ-model](https://arxiv.org/abs/2409.00533) | 严格 Tc 变分下界序列 + 稳定性上界 | J. Stat. Phys. 192, 69 (2025) |
| Kiessling-Altshuler-Yuzbashyan | Part II: Dispersive phonons | 推广到色散声子 | [J. Stat. Phys. 2025](https://link.springer.com/article/10.1007/s10955-025-03468-z) |
| Kiessling-Altshuler-Yuzbashyan | [Part III: Einstein phonons](https://www.alphaxiv.org/abs/2409.02121v2) | 进一步推广 | J. Stat. Phys. 2025 |
| Elezaby-Abanov | [QCP Tc bounds](https://arxiv.org/abs/2512.20009) | γ-model 更紧上界（Hessian 线代方法） | PRB 114, 024512 (2026) |
| Chubukov et al. | [QC Eliashberg review](https://arxiv.org/html/2506.11952v1) | 综述：γ-model 是量子临界超导的标准框架 | Annual Reviews 2026 |

**我们相对于竞争者的独特优势**：
1. 我们有**全链解析**的定理（引理 M 消除计算机成分），他们仍有数值环节
2. 我们有**区间算术 certificate** 基础设施，可以给出 sound 双侧包围
3. 我们的 Coulomb 推广（μ* > 0）是 KAY 系列尚未触及的
4. 我们的矩量阈值定理（p ≥ 2 iff）和极值约束族是 KAY 系列没有的
5. 我们的弱耦合级数（b₁ 闭式）和强耦合展开给出的渐近结构更丰富

**但是**：把我们的 P4 推广到 γ-model 是 Eliashberg 框架内的工作，不是「换赛道」。如果目标是 impact ≥ 2 的**新方向**，这只是半换赛道。

---

## 3. Top 3 候选（按 impact × feasibility × urgency 排序）

### 候选 #1：量子自旋系统的谱隙 certificate（Haldane 猜想的计算机辅助方向）

**一句话**：用算子不等式 $H^2 - \delta H \succeq 0$ 加区间算术 certificate，为 spin-1 Heisenberg 链（或更广泛的量子自旋系统）提供严格谱隙下界。

**为什么 impact ≥ 2**：
- Haldane 猜想是凝聚态理论最著名的开放问题之一（1983 年至今未证）
- AKLT 模型的谱隙已证（1D 和多个 2D 格），但**物理 Heisenberg 点无任何严格结果**
- 任何对 Heisenberg 链谱隙的严格正下界都是 PRL 级别；若能在热力学极限证明则是 Nature 级别
- 2026 年 Rai et al. 的 LTI-SDP 层级刚刚发表，领域正热

**我们的具体切入点**：
1. **直接方法**：实现 LTI-SDP 层级，但用我们的**区间算术**替代浮点 SDP，使得 certificate 是 sound 的（不依赖 SDP solver 的数值精度）。Rai et al. 的 AKLT 结果是浮点的 — 如果我们用 mpmath.iv 重做，每个 certificate 就是严格定理。
2. **Davis-Kahan 扰动路线**：AKLT Hamiltonian 与 Heisenberg Hamiltonian 之间的关系：$H_{\rm Heis} = H_{\rm AKLT} + V$，其中 $V$ 是可以显式写出的算子。如果 AKLT 的谱隙 $\Delta_{\rm AKLT}$ 大于 $\|V\|$（在某种适当范数下），则 Heisenberg 也有隙。
   - AKLT 谱隙 ≈ 0.350 已知
   - 关键问题：$\|V\|$ 有多大？
   - 这正是我们 Davis-Kahan + 区间算术的战场
3. **finite-size + 外推**：用 LTI-SDP 在有限尺寸（n=10,12,14,...）给出严格下界序列，用 Richardson 外推到热力学极限 — 我们在 P4 中已有 Richardson-in-M 的成熟经验。

**技术匹配度：High**
- $H^2 - \delta H \succeq 0$ ←→ 我们日常做的算子不等式
- 正算子分解 ←→ Perron-Frobenius 正性工具
- 区间算术 certificate ←→ 核心基础设施（mpmath.iv，1668-check 级别）
- Davis-Kahan 扰动 ←→ 核心技术
- 有限 N 到 N→∞ 的换序 ←→ P4 换序定理的直接翻版

**Impact 级别**：
- AKLT 的严格区间算术 certificate → J. Stat. Phys. (impact ~3)
- Heisenberg 链 finite-size 严格下界 → **PRL (impact ~2)**
- Heisenberg 链热力学极限严格正隙 → **Nature/Science (impact ~1)**

**时间成本**：
- Phase 1（AKLT 区间 certificate）：2-3 周
- Phase 2（Heisenberg finite-size）：4-6 周
- Phase 3（热力学极限）：未知（可能需要全新 idea）

**数据/文献就绪度**：**今天就能开工**。
- Rai et al. 2026 论文有完整方法描述
- AKLT/Heisenberg Hamiltonian 完全显式
- 不需要任何实验数据

**竞争态势**：
- Rai/Schuch/Baccari 组：SDP 层级的发明者，但他们用浮点 SDP（不 sound），且他们的兴趣是 frustration-free 系统
- Lemm 组：AKLT 2D 隙的证明者，使用数值辅助但不同工具
- **我们的独特优势**：(1) 区间算术使 certificate sound，(2) Davis-Kahan 扰动理论可以从 AKLT 过桥到 Heisenberg，(3) 我们已有 Richardson 外推的成熟工程

**风险因子**：
- Davis-Kahan 路线：‖V‖ 可能太大（AKLT → Heisenberg 的扰动不小），使得隙不能直接传递
- LTI-SDP 的 SDP 规模：spin-1 的 Hilbert 空间维度 = 3^n，n=14 时矩阵 3^14 ≈ 4.8M 维 — 区间 SDP 可能太贵
- 热力学极限：可能需要本质新想法，不仅仅是 finite-size 外推

**第一步**：
1. 读 [Rai et al. arXiv:2411.03680](https://arxiv.org/html/2411.03680v1) 全文，理解 LTI-SDP 的精确构造
2. 写脚本构造 spin-1 AKLT 的 LTI-SDP 在 n=4,6,8（小尺寸），用 mpmath.iv 包围 SDP 的可行性值
3. 计算 $H_{\rm Heis} - H_{\rm AKLT}$ 的算子范数（精确，有限链）

---

### 候选 #2：量子临界 Eliashberg Tc 的严格双侧界（γ-model）

**一句话**：把 P4 的全链解析定理和区间算术 certificate 迁移到 γ-model（频率幂律配对核 $V \sim 1/|\omega|^\gamma$），给出超越 Kiessling-Yuzbashyan Part I 和 Elezaby-Abanov 的最紧严格 Tc 界。

**为什么 impact ≥ 2**：
- Chubukov 2026 Annual Reviews 综述确认 γ-model 是量子临界超导的标准框架
- Elezaby-Abanov PRB 2026 的上界虽 "significantly tighter"，但仍不是 sound certificate（是线性代数 + 浮点）
- 如果我们给出**第一个双侧 sound certificate**（带区间算术），直接跳到 PRL
- 量子临界超导是凝聚态最热方向之一（cuprates, Fe-pnictides, heavy fermions 的统一框架）

**我们的具体切入点**：
1. γ-model 的 Eliashberg 方程写成 $K_\gamma \phi = \lambda \phi$，其中 $K_\gamma$ 是 Matsubara 频率上的紧自伴算子 — **与我们 P4 的 K[P] 同构**
2. γ = 1 对应 SDW QCP，γ = 1/3 对应铁磁 QCP，γ = 2/3 对应 nematic QCP — 一次性覆盖所有热门 QCP
3. **引理 M 的推广**：核单调性在 γ-model 下是否成立？如果成立，Einstein 极值定理直接搬过来
4. **弱耦合级数**：我们的 b₁ 闭式技术可以给 γ-model 的弱耦合前因子闭式
5. **Coulomb 推广**：γ-model + μ* > 0 是文献中完全空白的——KAY 没做、Elezaby-Abanov 没做

**技术匹配度：Very High**
- 算子结构与 P4 同构
- 所有已有工具（谱不等式、certificate、引理 M、q-单调化）理论上可搬
- 主要工作量是**适配**而非发明

**Impact 级别**：**PRL (impact ~2)**
- 特别是如果能证明 γ-model 的 Einstein 极值并给出 sound 双侧 Tc bracket

**时间成本**：3-5 周
- 核搬迁（改核函数）：1 周
- 引理 M 推广检查：1 周
- 区间 certificate 生成 + 验证：1-2 周
- 写稿：1 周

**数据/文献就绪度**：**今天就能开工**。
- γ-model 是纯理论框架，无需实验数据
- KAY Parts I-III 和 Elezaby-Abanov 完全公开
- 我们自己的 P4 代码就是起点

**竞争态势**：
- KAY 组（Rutgers）：已发三篇 J.Stat.Phys.，但他们**没有**区间算术
- Elezaby-Abanov：更紧上界但非 sound；且只有上界，无下界 certificate
- Chubukov 组：做物理不做严格化
- **我们的优势**：唯一同时有 (a) 全链解析框架、(b) sound certificate、(c) Coulomb 推广 的组

**风险因子**：
- γ-model 的核 $1/|\omega_n - \omega_m|^\gamma$ 在 γ < 1 时不是 trace-class（甚至可能不紧凑）— 需要仔细检查算子理论前提
- 引理 M 的核单调性可能不成立（γ-model 的核形状与 Einstein/Lorentzian 不同）
- J. Stat. Phys. 级别的竞争者已经在前面，我们需要实质超越而非跟风

**第一步**：
1. 验证 γ-model 核的紧凑性 / trace-class 性质（不同 γ 区间）
2. 检查引理 M（核单调性）在 γ-model 下是否成立
3. 用现有 P4 代码框架替换核函数，先跑数值证据

---

### 候选 #3：BCS 配对方程的严格 Tc 界（超越 mean-field 的冷原子普适比例）

**一句话**：用我们的 Eliashberg 算子不等式技术对 BCS 配对方程建立严格 Tc/E_F 上界，为冷原子（BCS-BEC crossover）普适比提供 operator-level 证明。

**为什么 impact ≥ 2**：
- 冷原子实验在 BCS-BEC crossover 中精确测量 Tc/T_F（MIT, JILA, ENS 组长期竞争课题）
- BCS 配对方程的 Tc 上界 → 幺正极限的 Tc/T_F 严格上界 → 与 QMC 和实验直接对比
- 第一个从**算子不等式**给出的 BCS Tc/T_F 普适上界将是 PRL 级别
- Frank-Hainzl-Seiringer 系列的数学物理传统认可这类严格结果

**我们的具体切入点**：
1. BCS 线性化方程在 Tc 处：$K_V \Delta = \Delta$，其中 $K_V$ 是配对核 — 与 Eliashberg 方程同构
2. 接触势 $V = -g\delta(r)$ 的 BCS 方程有精确的 Matsubara 频率表示 — 就是 γ = 0（对数核）的特例
3. 我们的 Einstein 极值定理 → BCS 框架下的"最优声子谱"定理
4. 维度从 Matsubara 频率扩展到动量 — 需要处理 3D 问题，但算子结构不变

**技术匹配度：High**
- BCS 方程 ≡ 特殊 Eliashberg 方程
- 紧自伴算子谱理论直接适用
- 区间 certificate 可搬
- **新挑战**：需要同时处理频率和动量（不仅仅是 Matsubara 求和）

**Impact 级别**：**PRL / Nature Physics (impact ~2)**

**时间成本**：4-8 周
- BCS 配对核的算子分析：2 周
- 从频率到 (频率,动量) 的推广：2-3 周
- Certificate 生成：1-2 周
- 写稿：1 周

**数据/文献就绪度**：**基本就绪**。
- BCS 理论完全标准
- 需要读 Frank-Hainzl-Seiringer 系列论文（~1 周文献消化）
- 幺正极限有精确 QMC 数据可对比

**竞争态势**：
- Seiringer 组：BCS 泛函数学理论的主要建设者，但方法是变分 + 函数分析，不是 certificate
- 无人用区间算术做 BCS 的 Tc 界
- **我们的优势**：certificate 技术是全新的切入角

**风险因子**：
- BCS mean-field Tc 与真正多体 Tc 之间有 gap — 我们只能给 mean-field 框架内的严格界
- 如果这个 gap 太大，物理意义受限（但 BCS Tc 上界仍是 true 上界 — 严格性不受影响）
- 从 1D 频率到 3D (频率,动量) 空间的扩展可能有本质困难

**第一步**：
1. 读 Frank-Hainzl-Seiringer "The BCS Functional for General Pair Interactions"
2. 把 BCS 接触势方程写成 Matsubara 求和形式，验证与 P4 的 K[P] 的对应关系
3. 检验 Einstein 极值定理（引理 M + 定理 E）在 BCS 核下是否成立

---

## 4. 超级候选：Haldane 猜想的 Davis-Kahan 攻击（Impact = 1）

**一句话**：用 AKLT 模型已知谱隙 + Davis-Kahan sin Θ 扰动定理 + 区间算术 certificate，证明 spin-1 Heisenberg 反铁磁链在热力学极限有严格正谱隙。

**为什么 impact = 1**：
- Haldane 猜想（1983）是凝聚态理论/数学物理领域最著名的 40+ 年开放问题之一
- Haldane 因相关工作获 2016 诺贝尔奖
- 完整证明会是 Nature/Science 级别，可能是年度最重要数学物理结果

**为什么值得赌**：
- **没有人尝试过 Davis-Kahan 路线**。所有已有尝试都在 frustration-free 框架内工作（因为 detectability lemma 需要），而 Heisenberg 链 explicitly 不是 frustration-free。我们的 Davis-Kahan 扰动理论提供了**完全不同的逻辑**：不从 Heisenberg 直接出发，而是从已知有隙的 AKLT 出发，控制扰动。
- 具体来说：
  - 设 $H = H_{\rm AKLT} + \epsilon V$，其中 $\epsilon \in [0,1]$，$\epsilon = 0$ 是 AKLT，$\epsilon = 1$ 是 Heisenberg
  - AKLT 隙 $\Delta_0 \approx 0.350$
  - 如果我们能证明 $\|V\| < \Delta_0$（在适当意义下），则 $H(\epsilon)$ 对所有 $\epsilon$ 有隙
  - 即使 $\|V\| \geq \Delta_0$，如果我们有更细的信息（比如 $V$ 投影到低能子空间后的有效范数更小），Davis-Kahan 仍可能工作

**核心技术问题**（make-or-break）：
1. $H_{\rm AKLT}$ 和 $H_{\rm Heis}$ 之间的算子距离有多大？
   - $H_{\rm AKLT} = \sum_i P_i^{(2)}$（投影到 spin-2 扇区），$H_{\rm Heis} = \sum_i \vec{S}_i \cdot \vec{S}_{i+1}$
   - 它们之间的差 $V = H_{\rm Heis} - c \cdot H_{\rm AKLT}$ 是否"小"？
   - **这是可以计算的** — 这是我们第一步要做的
2. 如果 ‖V‖ 太大（很可能），能否用**逐步扰动**？
   - 定义插值族 $H(\theta) = \cos\theta \cdot H_{\rm AKLT} + \sin\theta \cdot (H_{\rm AKLT} + V')$
   - 在每一步用 finite-size certificate 确认谱隙 > 0
   - 然后用 Richardson 外推到热力学极限
3. 或者：用 **LTI-SDP + 区间算术**在 Heisenberg 点**直接**给出 finite-size 谱隙界，再外推？

**技术匹配度：Very High**（核心匹配），但问题的内在困难度极高。
- Davis-Kahan：核心工具 ✓
- 区间算术：核心基础设施 ✓
- 算子扰动：核心经验 ✓
- finite-size → thermodynamic limit：P4 换序定理经验 ✓
- 但：可能需要处理指数大 Hilbert 空间的 SDP，这是计算瓶颈

**时间成本**：
- Phase 1（可行性判定 — AKLT vs Heisenberg 的算子距离）：**1-2 周**
- Phase 2（如果可行，finite-size certificate）：6-12 周
- Phase 3（thermodynamic limit）：未知

**第一步**（今天就做）：
1. 显式计算 $V = H_{\rm Heis} - c \cdot H_{\rm AKLT}$ 在 2-site 和 3-site 子空间上的矩阵表示
2. 求 ‖V‖ / Δ_AKLT 的比值 — 如果 < 1，路线可行；如果 ≫ 1，直接 Davis-Kahan 路线死亡
3. 如果 ‖V‖ 太大，转向 LTI-SDP 直接方法（候选 #1 的 Phase 2）

**风险因子**：
- **高风险**。Haldane 猜想 open 了 43 年是有原因的。
- AKLT → Heisenberg 的扰动可能太大，Davis-Kahan 不够
- finite-size SDP 的计算瓶颈（spin-1 系统，3^n 维 Hilbert 空间）
- 即使 finite-size 有隙，外推到热力学极限可能需要本质新想法
- 但：**Phase 1 的成本只有 1-2 周**，即使失败我们只损失 1-2 周，而成功则是改变职业轨迹的结果

---

## 5. 方向比较矩阵

| 维度 | #1 自旋系统谱隙 | #2 γ-model Tc | #3 BCS Tc | 超级候选 Haldane |
|---|---|---|---|---|
| Impact (成功时) | 2-1 | 2 | 2 | **1** |
| 技术匹配 | High | Very High | High | Very High |
| 时间成本 (Phase 1) | 2-3 周 | 3-5 周 | 4-8 周 | **1-2 周** |
| 文献就绪 | 今天 | 今天 | 1 周文献 | 今天 |
| 竞争压力 | 中（Rai/Schuch 刚发） | 高（KAY 在前） | 低（Seiringer 不做 cert） | 低（无人试此路线） |
| 失败成本 | 低 | 低 | 中 | 低（Phase 1） |
| 独特优势 | 区间 cert + D-K | 唯一 sound cert | 全新角度 | 全新路线 |

---

## 6. 建议行动顺序

1. **立即**（本周内）：超级候选 Phase 1 — 计算 ‖H_Heis - c·H_AKLT‖ 的精确值。这是 1-2 天的脚本 + 区间运算工作。
   - 如果比值 < 1 → 全力投入 Haldane 方向
   - 如果比值 ≫ 1 → 转向候选 #1 的 LTI-SDP 方法
2. **并行**（本周起）：候选 #2 的 Phase 0 — 验证 γ-model 核的紧凑性和引理 M 是否成立（1 天数值探测）
3. **根据 1-2 的结果**：选定主方向，集中推进

---

## 文献引用索引

### 谱隙 / Haldane 猜想
- Rai, Kull, Emonts, Tura, Schuch, Baccari. "A Hierarchy of Spectral Gap Certificates for Frustration-Free Spin Systems." [Quantum (2026)](https://quantum-journal.org/papers/q-2026-04-13-2065), [arXiv:2411.03680](https://arxiv.org/html/2411.03680v1)
- Pomata, Wei. "Demonstrating the AKLT spectral gap on 2D degree-3 lattices." [arXiv:1911.01410](https://arxiv.org/html/1911.01410v2)
- Lemm, Sandvik, Wang. "Existence of a Spectral Gap in the AKLT Model on the Hexagonal Lattice." [PRL 124, 177204 (2020)](https://link.aps.org/doi/10.1103/PhysRevLett.124.177204)
- Haldane conjecture status: [emergentmind](https://www.emergentmind.com/open-problems/haldane-conjecture-integer-spin-heisenberg-chain-gap)

### Eliashberg Tc 界
- Kiessling, Altshuler, Yuzbashyan. "Bounds on Tc... I: The γ-model." [J. Stat. Phys. 192, 69 (2025)](https://arxiv.org/abs/2409.00533)
- Kiessling, Altshuler, Yuzbashyan. "Bounds on Tc... II: Dispersive Phonons." [J. Stat. Phys. (2025)](https://link.springer.com/article/10.1007/s10955-025-03468-z)
- Elezaby, Abanov. "Superconductivity near a QCP: Bounds..." [PRB 114, 024512 (2026)](https://arxiv.org/abs/2512.20009)
- Chubukov et al. "Quantum Critical Eliashberg Theory." [Annual Reviews (2026)](https://arxiv.org/html/2506.11952v1)

### BCS 严格结果
- Frank, Hainzl, Seiringer. "The BCS Functional for General Pair Interactions." [arXiv:math-ph/0703086](https://arxiv.org/html/math-ph/0703086v2)
- "Enhanced superconductivity at a corner for the linear BCS equation." [Forum Math. Sigma (2025)](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/enhanced-superconductivity-at-a-corner-for-the-linear-bcs-equation/39F8C08749921E9C64E9EF8C5E5076A2)
- Seiringer publication list 2024-2026: [ISTA](https://pub.ista.ac.at/~rseiring/Robert_Seiringer/Robert_Seiringers_Homepage_files/List%20of%20Publications.pdf)

### Planckian 界
- Cao. "Planckian Bound on Quantum Dynamical Entropy." [PRL 137, 030402 (2026)](https://journals.aps.org/prl/abstract/10.1103/wd91-v58h), [arXiv:2507.20914](https://arxiv.org/html/2507.20914v5)
- Chowdhury. "Planckian Bounds via Spectral Moments of Optical Conductivity." [arXiv:2512.09979](https://arxiv.org/html/2512.09979v1)

### SYK
- Basu, Kothari, Midha. "Sharp Bounds on Ground State Energy of the SYK Model." [arXiv:2607.27185](https://arxiv.org/abs/2607.27185)
- He. "Free energy and spectral edge of the SYK model." [arXiv:2607.18998](https://www.alphaxiv.org/zh/abs/2607.18998)

### 其他
- Anshu, Arad, Gosset. "An area law for 2D frustration-free spin systems." [STOC 2022](https://dl.acm.org/doi/10.1145/3519935.3519962)
- MBL review (2026): [Rep. Prog. Phys.](https://iopscience.iop.org/article/10.1088/1361-6633/ad9756)
- Google QEC below threshold: [Nature (2024)](https://www.nature.com/articles/s41586-024-08449-y)
- Hartnoll. "Planckian dissipation in metals." [RMP 94 (2022)](https://research-repository.st-andrews.ac.uk/bitstream/handle/10023/26867/Hartnoll_2022_RMP_Colloquium_Planckian_VoR.pdf?sequence=1)
