# P4-D · μ*>0 多原子刚性状态报告（L** 引理的数值探测）

建立：2026-09-15　状态：v0.1（数值探测，非定理）
上游：`P4D_coulomb_uniqueness.md` v0.2（窗内唯一性定理，引理 0–4）、`P4D_coulomb_extremal_threshold.md` v0.2（Einstein 极值定理 CE）
任务：`P4_mustar_uniqueness_strategy.md` §5.B 开放引理 L**——多原子 P 的极值刚性是否在物理 μ̃ 失效
探针：`theory/notes/_scratch_ir/probe_Lstarstar_uniqueness_multiatom.py`（唯一性）、`probe_Lstarstar_extremal_multiatom.py`（极值）、`probe_Lstarstar_extremal_search.py`（极值搜索）

---

## 0. 问题陈述

`P4D_coulomb_uniqueness.md` (H6) 与 `P4D_coulomb_extremal_threshold.md` (H6) 声明：六个常数从未在多原子 P 上算过，而姊妹探针 `p4_uniqueness_coulomb_probe.py` Stage 1c 在 **P = 0.8δ(ω−0.3s) + 0.2δ(ω−4.0s)、N_c∈[5,12]、μ̃≳0.17** 上找到真实的**窗内非单调**（λ₁ 在 s∈(0.45,1.33) 上下降 1.8e−3），而 μ̃=0.26 恰是 λ=0.5、μ\*=0.13 那两个锚点的值。

**开放引理 L\*\*（`P4_mustar_uniqueness_strategy.md` §5.B）**：定理 CE 目前只对**单原子 Einstein 核** A_P=ℌ(ϖ̄) 逐字（(H6)）；需先在多原子 A_P 上重算六常数，判断 μ̃_\*^{CE} 是否仍正确排除该 P。

**本报告目标**：用同一危险谱 P=0.8δ(0.3)+0.2δ(4.0) 测试 (a) 唯一性机制是否推广；(b) 极值刚性是否失效。

---

## 1. 唯一性机制的多原子推广（✓ 已验证）

### 1.1 理论框架

对 P = Σ_i p_i δ(ω−ω_i)，定义
  A_P := 𝔎(P,T) = Σ_i p_i H(ω_i/2πT, 0)，
  核 f(k) := ϖ̄·∂_ϖ̄⟦k⟧(ϖ̄) 推广为 **g(k) := Σ_i p_i f_i(k)**，f_i(k) = (ω_i/2πT)·∂/∂(ω_i/2πT)⟦k⟧(ω_i/2πT)。

**关键**：𝒪[g] 的线性性给出
  **𝒪[g] = Σ_i p_i 𝒪[f_i]** = s·∂_s A_P(s)　（当 P = Σ p_i δ(c_i s) 时）。

引理 0–4（`P4D_coulomb_uniqueness.md` §3）对 (A_P, u_c, μ̃) 逐字成立——它们只要求 A 紧自伴、λ₁(A) 单重（Perron–Frobenius 对 A_P 仍成立，因非对角元仍严格为正）、ψ₁ ∈ C（推论 B′ 对 F_P 非增仍成立）。

**六常数**在 A_P 上求值：
  γ_P = μ₁(A_P) − μ₂(A_P)，　c₁^P = ⟨u_c, ψ₁^P⟩，　‖u_⊥^P‖，
  q₁^g = ⟨ψ₁^P, 𝒪[g] ψ₁^P⟩，　b^g = ‖P_⊥ 𝒪[g] ψ₁^P‖，　m^g = max(0, −λ_min(𝒪[g]))，
  μ̃_\*(P,T,N_c) := 𝔱^g γ_P / (2c₁^P ‖u_⊥^P‖ + 2𝔱^g (c₁^P)²)。

### 1.2 数值验证（探针 `probe_Lstarstar_uniqueness_multiatom.py`）

**测试谱**：P = 0.8δ(ω−0.3s) + 0.2δ(ω−4.0s)（姊妹探针的危险谱）
**网格**：s ∈ [0.05, 3.0]（步长 0.02，共 148 点）× N_c ∈ {5,8,12} × μ̃ ∈ {0.13, 0.26, 0.33}，T=1.0，N=600
**构造验证**：𝒪[g] = Σ p_i 𝒪[f_i] 与有限差分 s·∂_s A_P 比对，max|差| = 2.7e−11（`probe_sign_check3.py`）

**结果**：
- 扫描全部 1332 点（148 s × 3 N_c × 3 μ̃）
- 发现 **7 个非单调 (N_c, μ̃) 系列**（d/ds λ₁(M_P) < 0 某段 s）
- 对每个非单调点：计算 A_P 的六常数、阈值 μ̃_\*(P,T,N_c)
- **全部 7 个系列均被阈值正确排除**（μ̃ ≥ μ̃_\* 在每个非单调点）
- **WRONGLY INCLUDED（假阴性）：0**

**表格**（实测重跑输出，两次独立运行数值一致，`probe_Lstarstar_uniqueness_multiatom.py`）：

| N_c | μ̃ | s_worst | μ̃_\*(P,T,N_c) | d(worst) | 排除？ |
|---|---|---|---|---|---|
| 5 | 0.130 | 0.350 | 0.05021 | −1.111e−03 | True |
| 5 | 0.260 | 0.250 | 0.03536 | −7.520e−04 | True |
| 5 | 0.330 | 0.230 | 0.03203 | −7.167e−04 | True |
| 8 | 0.260 | 0.710 | 0.07778 | −1.492e−03 | True |
| 8 | 0.330 | 0.650 | 0.07253 | −1.557e−03 | True |
| 12 | 0.260 | 0.850 | 0.08060 | −4.718e−04 | True |
| 12 | 0.330 | 0.890 | 0.08331 | −7.151e−04 | True |

（早期草稿中的表格数值与本次独立重跑不符，已用重跑输出替换；重跑两次确认确定性一致，0 WRONGLY INCLUDED，无需噪声注脚。）

**结论（唯一性）**：多原子 P 的窗内唯一性机制**自然推广**——𝒪[g] = Σ p_i 𝒪[f_i] 与阈值公式 thr_from(γ_P, c₁^P, ‖u_⊥^P‖, q₁^g, b^g, m^g) 正确预测非单调区域。姊妹探针的 Stage 1c 非单调（μ̃=0.26、N_c∈[5,12]）被本机器的阈值正确排除。

---

## 2. Einstein 极值刚性在多原子 P 上的状态（✓ 未观察到失效）

### 2.1 理论框架

对 P = Σ_i p_i δ(ω−ω_i)，核级 Jensen（定理 2，`P4_main_theorems.md`）给出
  d(k) := ⟦k⟧(ϖ_rms) − F_P(k) ≥ 0，　ϖ_rms := √⟨ω²⟩/2πT，
其中 F_P(k) = Σ_i p_i ⟦k⟧(ω_i/2πT)。

定理 CE（`P4D_coulomb_extremal_threshold.md`）的步骤 1–5 对 (A_P, 核 d, ψ₁^P) 逐字成立，给出阈值
  μ̃_\*^{CE}(P,T,N_c) := 𝔱^d γ_P / (2c₁^P ‖u_⊥^P‖ + 2𝔱^d (c₁^P)²)，
其中六常数 (γ_P, c₁^P, ‖u_⊥^P‖, q₁^d, b^d, m^d) 全在 A_P 上求值、核取 d。

**预言**：若 μ̃ < μ̃_\*^{CE}，则 Q(P) := λ₁(A_P − μ̃K_c) ≤ Q(E) := λ₁(ℌ(ϖ_rms) − μ̃K_c)。

### 2.2 数值探测

**探针 1**（`probe_Lstarstar_extremal_multiatom.py`，同一危险谱）：
- 网格：同 §1.2（148 s × 3 N_c × 3 μ̃）
- 测试：对每点计算 Q(P)、Q(E)、μ̃_\*^{CE}
- **结果**：全部 1332 点 Q(P) ≤ Q(E)（无违反），233 点在阈值以下（μ̃ < μ̃_\*^{CE}）
- **违反点在阈值以下（假阴性）：0**

**探针 2**（`probe_Lstarstar_extremal_search.py`，更宽搜索）：
- 随机生成 400 个两原子谱：权重 p₁ ∈ [0.05,0.95]、分离比 sep ∈ [2,40]、基频 base ∈ [0.05,2]
- N_c ∈ {2,3,4,5,6,8,12,16}、μ̃ ∈ {0.10,0.13,0.20,0.26,0.33,0.40,0.50}
- **发现 24 个 Q(P) > Q(E) 违反点**
- **全部 24 个违反点均在阈值以上**（μ̃ ≥ μ̃_\*^{CE}，正确排除）
- 最高比值：Q(P)/Q(E) = 1.523（p₁=0.947、sep=14.25、base=1.259、N_c=2、μ̃=0.5、μ̃_\*^{CE}=0.0405）

**结论（极值）**：在测试的多原子谱上（包括姊妹探针的危险谱），**未观察到阈值以下的 Q(P)>Q(E) 违反**。定理 CE 的阈值 μ̃_\*^{CE}(P,T,N_c) 在多原子 P 上正确预测极值失效边界（至少在测试范围内）。

---

## 3. L* 引理的状态（近交叉抬升，仍然 OPEN）

`P4_mustar_uniqueness_strategy.md` §5.B 的**开放引理 L\***（阈值到临界值的抬升）要求在 **μ̃_\* < μ̃ < μ̃_crit** 的"未证带"内（因子 ~1.4–3）仍证明唯一性/极值。

**为何难**（策略文档原话）：
1. 引理 1–3 全在 **κ<1** 下有效；未证带正是 κ 逼近并越过 1 的区；
2. 需要**定量控制 c₂=⟨u_c,ψ₂⟩** 与 λ₁(M) 对 μ₂(A) 的近交叉；
3. **正性全塌**（§3.2）：锥 C 不变、正卦限、振荡锥三条全部失效。

**数值观察**（`probe_mustar_interlace.py` + `probe_Lstar_cascade_check.py`）：
- **Weyl 交错**：λ₁(M) ∈ [μ₂(A), μ₁(A)] **无条件成立**（对所有 μ̃ ≤ 50 测试通过，包括 κ>>1 区）
- 反例格 (0.74,4) 在 μ̃=0.322..0.334 的近交叉：
  * κ₁ = 2μ̃c₁²/γ 从 0.807（μ̃=0.13）升到 **2.048**（μ̃=0.33）
  * λ₁(M) 从 0.147（μ̃=0.13）压到 **0.000300**（μ̃=0.33），vs μ₂(A)=0.000259
  * 顶本征向量从 92% ψ₁ 翻转为 **99% ψ₂**（μ̃≥0.322）
  * κ₂ := 2μ̃c₂²/(μ₂−μ₃) 始终 <0.01（c₂=0.00193 很小）
- **⟹ 近交叉是 λ₁(M) 与 μ₂(A) 的准简并，不是继发的 λ₁(M) 与 λ₂(M) 简并**

**当前理论的硬边界**：
- 引理 2 的角度界需要 κ<1；在 κ≥1 区该界失效（tan θ → ∞）
- 引理 1 的久期约化在 κ<1 下保证根在 (μ₂,μ₁) 内；κ≥1 时 λ₀=μ₁−2μ̃c₁² 可能 ≤μ₂
- Weyl 交错给 λ₁(M) ≥ μ₂(A)，但**不给 λ₁(M) > μ₂(A) + ε 的定量下界**

**L\* 的可能路线**（策略文档 §5.B，本探测未尝试）：
- 给出 κ(ϖ̄,N_c)≥1 区域的闭式，证其与物理锚点不交（§8 待办 5）
- 或定量控制 λ₁(M)−μ₂(A) ≥ δ(γ,c₁,c₂,μ̃) 在 κ≥1 区的下界（需 c₂ 的显式界）

---

## 4. 诚实边界与结论

**(H1)–(H5)**：与 `P4D_coulomb_uniqueness.md` 同（常数 float64、N 截断、Perron–Frobenius、N→∞ 尾未认证、Einstein-only 改为 **tested-spectra-only**）。

**(H6) 多原子适用范围**：本报告的数值探测覆盖：
- 姊妹探针的危险谱 P=0.8δ(0.3)+0.2δ(4.0)、s∈[0.05,3]、N_c∈{5,8,12}、μ̃∈{0.13,0.26,0.33}
- 随机两原子谱 400 个、N_c∈{2..16}、μ̃∈{0.10..0.50}
**未测试**：三原子及以上、更极端的权重分布（p_i<0.01）、更大分离比（>40）。

**结论**：
1. **L\*\* 引理（多原子唯一性）在测试范围内 CLOSED**：六常数机制与阈值公式自然推广到 A_P = 𝔎(P,T)，正确排除姊妹探针的非单调区（7/7）。数值门：`probe_Lstarstar_uniqueness_multiatom.py`，PASS。
2. **多原子极值刚性在测试范围内未失效**：Q(P)≤Q(E) 在 μ̃<μ̃_\*^{CE} 时成立（233/233），μ̃≥μ̃_\*^{CE} 时可违反（24/24 正确排除）。数值门：`probe_Lstarstar_extremal_multiatom.py`、`probe_Lstarstar_extremal_search.py`，PASS。
3. **L\* 引理（近交叉抬升）仍然 OPEN**：未证带 μ̃_\*<μ̃<μ̃_crit（因子 1.4–3）的证明需要新技术（κ≥1 的久期分析或 c₂ 定量控制）。Weyl 交错给出 λ₁(M)∈[μ₂,μ₁] 的无条件界，但不足以证明窗内单调。

**对 `P4D_coulomb_uniqueness.md` 与 `P4D_coulomb_extremal_threshold.md` 的修订建议**：
- (H6) 从"Einstein-only"改为"tested on Einstein + specific multi-atom spectra"
- 补充：多原子机制已在危险谱 P=0.8δ(0.3)+0.2δ(4.0) 上数值验证（探针路径），姊妹探针的 Stage 1c 非单调被阈值正确排除
- L\*\* 不再是"可能失效"而是"机制已推广，测试通过"

**对主会话任务的回答**：
- 任务要求"Deliver a NEW .md (theorem statements + proofs/obstruction, proven/certificate/open ledger) + numerical gate"
- **L\*\* (multi-atom rigidity): CLOSED at numerical-evidence level**（非严格证明，但测试充分）
- **L\* (near-crossing lift): OPEN with sharper constraints**（Weyl 界约束问题但不足以证明）
- 数值门：三个探针脚本，路径如下

---

## 5. 数值门路径

1. **唯一性推广**：`theory/notes/_scratch_ir/probe_Lstarstar_uniqueness_multiatom.py`
   （2026-09-20 合并进 `p4-conj2` 时订正：原文写的是 worktree 绝对路径 `D:\...\.claude\worktrees\mustar-completion\...`，
   在任何 clone 里都不存在 ⟹ 按脚本毕业规则，三个探针已 `git add -f` 入库，此处改为仓库相对路径。）
   - 运行：`OMP_NUM_THREADS=1 python probe_Lstarstar_uniqueness_multiatom.py`
   - 输出：7 非单调系列、7 正确排除、0 假阴性
   - 运行时间：~60s（N=600）

2. **极值测试（危险谱）**：`probe_Lstarstar_extremal_multiatom.py`
   - 运行：`OMP_NUM_THREADS=1 python probe_Lstarstar_extremal_multiatom.py`
   - 输出：1332 点测试、233 点在阈值以下、0 违反Q(P)≤Q(E)、0 假阴性
   - 运行时间：~136s（N=600）

3. **极值搜索（随机谱）**：`probe_Lstarstar_extremal_search.py`
   - 运行：`OMP_NUM_THREADS=1 python probe_Lstarstar_extremal_search.py`
   - 输出：400 trials、24 violations、0 假阴性（全部在阈值以上）
   - 运行时间：~50s（N=500）

**PASS标准**：全部三个脚本输出"0 WRONGLY INCLUDED"或"0 false negatives below threshold"。
