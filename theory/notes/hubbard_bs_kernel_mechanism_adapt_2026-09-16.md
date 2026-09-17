# Hubbard Bethe-Salpeter 核机制适配：定理 I（Δ 通道免疫）

**建立**：2026-09-16 · Track A 子会话（并行方向）
**状态**：机制对照 + 瓶颈诊断 · **未提交**

---

## 0. 执行摘要

**Phase 1 可行性判断**：**Medium**（4–6 周，2 个引理需新证明）

- **引理 1（Bochner）**：需适配 — 自旋涨落核 χ_spin(q) 不满足 Bochner 类条件（q 空间非正定）
- **引理 2（偶扇区压缩）**：需适配 — 动量空间对称性与 Matsubara 指标对称性不同构
- **引理 3（ψ≥0 不等式）**：直接套用 — 与核的具体形式无关
- **引理 4（秩一泄漏）**：**根本障碍** — d-波节点使 g(k) 变号，a_k≥0 条件失效

**核心机制匹配度**：
- ✓ 自伴正算子 + 秩一吸引扰动的一般结构**完全匹配**
- ✓ Δ/Z 通道分离的物理图像**直接对应**
- ✗ 具体核函数的解析性质（正定性、Toeplitz 结构）**不适用**
- ✗ d-波符号结构（节点）是**新障碍**

**建议优先级**：在完成 Phase 1 适配前，先平行启动 **Track C（实验数据直接匹配）** — Hubbard 的 QMC Tc(U,μ) 数据可以直接与定理 I 的推论（免疫预算 c₁²）做数值检验，不需要等解析适配完成。

---

## 1. 核结构对照

### 1.1 Einstein（定理 I 的基底）

**算子**：ℓ²(ℕ₀)，Matsubara 频率指标 n，
```
M = -L1[φ_Z] + P[φ_Δ]
L1_n = diag( 2Σ_{k=1}^n φ_Z(k) / (2n+1) )
P = L2 + L3,   L2_nm = (1-δ_nm)φ_Δ(|n-m|)u_n u_m,   L3_nm = φ_Δ(n+m+1)u_n u_m
u_n = (2n+1)^{-1/2}
```

**核类**：𝒫 = { f(k) = ∫₀^∞ ω²/(ω²+k²) dP(ω) }（Bochner 类，k∈ℕ₀ 上正定）
- φ_Z = f₊ + f₋（总谱权，≥0）
- φ_Δ = f₊ - f₋（配对通道，符号：f₊(0)<0 吸引，f₊(k)>0 for k≥1）

**关键性质**（支撑定理 I）：
1. f∈𝒫 ⟹ 延拓到 ℤ 上正定，f̂(θ)≥0
2. 反射延拓 ⟹ T_f + H_f ⪰ f̂(π)I（引理 2）
3. L1 对角增长 ∝ Σφ 逐点压住 P 的负部分（引理 4）

### 1.2 Hubbard Bethe-Salpeter（目标适配对象）

**算子**：动量空间 k∈BZ（第一 Brillouin 区），Matsubara 频率 ωₙ，
```
Δ_n(k) = (T/Z_n(k)) Σ_{k'} V_d(k,k') Δ_n(k')
```
其中 V_d(k,k') 是 d-波投影的配对核。

**核的物理形式**（从 [Scalapino et al., arXiv:cond-mat/0612579](https://arxiv.org/html/cond-mat/0612579v1)）：
```
V_d(k,k') = g(k) · K(k,k') · g(k')
g(k) = cos(k_x) - cos(k_y)    （d-wave form factor，**有节点**：g=0 沿对角线）
K(k,k') ≈ (3/2)Ū²(T) χ_spin(k-k')
```
χ_spin(q) 是自旋涨落磁化率，在 q=(π,π) 附近发散（反铁磁不稳定性）。

**分解到 Z/Δ 通道**（类比 Einstein 的 φ_Z/φ_Δ）：
- **Δ 通道**（配对）：K_Δ = (3/2)Φ_m（磁通道，S=1，吸引）
- **Z 通道**（质量重整化）：包含电荷通道 (1/2)Φ_d（S=0，排斥）+ 完全不可约顶点 Λ_irr
- 净效应：磁通道主导 ⟹ d-波配对，电荷通道部分抵消

**关键差异**：
- 核定义在**连续动量空间** k∈BZ（2D Torus），不是离散 ℕ₀
- χ_spin(q) 在 q 空间**非正定**（有符号振荡，只在 q≈(π,π) 处增强）
- d-波因子 g(k) **变号**（节点线），Einstein 的 u_n>0 恒正

---

## 2. 引理逐条对照

### 引理 1（Bochner / Poisson）

**原陈述**（Einstein）：f∈𝒫 ⟹ f 延拓到 ℤ 上正定，f̂(θ) = Σ_k f(k)e^{-ikθ} ≥ 0 ∀θ。

**Hubbard 适配**：**需适配**（不能直接套）

**瓶颈**：
1. χ_spin(q) 在动量空间**不是 Bochner 类**：
   - Einstein 的 ω²/(ω²+k²) 是 Laplace 变换（Matsubara 频率 → 虚时），Bochner 定理保证 k 空间正定
   - Hubbard 的 χ_spin(q) 是 RPA 或 DMFT 计算的物理响应函数，q 空间有**符号振荡**（短程反铁磁关联 vs 长程筛选）
   
2. 物理上 χ_spin(q) 在 q≈(π,π) 处峰化（commensurate 反铁磁），q≈0 处被库仑排斥压低 — 这是**非单调、有结构**的，不满足正定性。

**可能的适配路径**（未验证）：
- **放弃全局正定性**，只要求**局部正性**：在费米面附近的 k-patch，g(k)²χ_spin(k-k') 的主本征值远大于次本征值。
- **数值验证**：直接计算 Hubbard QMC 的 χ_spin(q)，做 Fourier 变换 χ̂(x)，检查实空间是否 χ̂(x)≥0（这等价于 q 空间正定）— 预期**不成立**，但可以定量评估"违反程度"。

**时间估计**：2 周（文献调研 + 数值检验 + 可能的局部化论证）。

---

### 引理 2（偶扇区压缩：T_f + H_f ⪰ f̂(π)I）

**原陈述**（Einstein）：对 ℤ 上正定 f，反射延拓到全整数 ⟹ Toeplitz + Hankel 矩阵和有下界。

**Hubbard 适配**：**需适配**（动量空间对称性不同）

**瓶颈**：
1. Einstein 的 T_f + H_f 来自 Matsubara 频率的**正负分支反射对称性**（ω_n ↔ ω_{-1-n}）
2. Hubbard 在**动量空间**，对称性是 BZ 的晶格对称群（C₄ᵥ，d-wave 的 B₁ᵍ 表示）
3. L2（Toeplitz 型，|n-m| 依赖）↔ 动量转移 q=k-k'；L3（Hankel 型，n+m+1）↔ **没有直接对应**（负频反射是 Matsubara 特有的）

**可能的适配路径**：
- **放弃 Hankel 项**：Hubbard 的核只有 K(k,k')=g(k)λ(k-k')g(k')（momentum-transfer 型），没有"n+m"型的配对。引理 2 的下界可能需要**换工具** — 用动量空间的 Schur 补引理或 Gershgorin 圆盘定理。
- **或者**：在 Matsubara 频率指标上仍有 L3 型项（Nambu 方程推导，见 `nambu_weight_audit` 第 61 行），但动量部分只有 L2。需要**分离频率与动量**的作用。

**时间估计**：1 周（明确 Hubbard 的 L3 是否存在 + 如果不存在，换下界工具）。

---

### 引理 3（ψ(ω) ≥ 0 初等不等式）

**原陈述**（Einstein）：ψ(ω) = 2ω²/(ω²+1) + πω/sinh(πω) - 1 ≥ 0。

**Hubbard 适配**：**直接套用**

**理由**：引理 3 只依赖 Matsubara 求和的 telescoping（Einstein 核 ∫dP(ω) 在 Σφ_Z - ρ 中的正性），与核的具体 k 空间形式无关。只要 Hubbard 的核仍是"玻色频率 ω 的 Lorentz 型"（χ_spin(q,Ω) ~ Ω²/(Ω²+ω_sf²)，自旋涨落特征频率），引理 3 的不等式逐字成立。

**无障碍**。

---

### 引理 4（秩一泄漏：a_n ≥ 0 for n≥1）

**原陈述**（Einstein）：𝒪'[f] = L1 + P ⪰ -ρ(f)e₀e₀ᵀ，证明路线：
```
diag(a),  a_n = [2Σ_{k=1}^n f(k) - f(0) + f̂(π)] / (2n+1)
a_0 = -ρ(f),  a_n ≥ 0 for n≥1  （引理 3 的 ψ≥0 保证）
```

**Hubbard 适配**：**根本障碍**

**瓶颈**：
1. Einstein 的论证依赖 **u_n = (2n+1)^{-1/2} > 0 恒正**，所以 diag(a) 的符号直接决定 𝒪' 的负部分。
2. Hubbard 的 d-wave 因子 **g(k) = cos(k_x) - cos(k_y) 变号**：
   - 节点线：k_x = k_y（对角线）
   - g(k) > 0 在 (0,π)×(π,0) 区域（吸引）
   - g(k) < 0 在 (π,0)×(0,π) 区域（排斥符号相反）
3. 类比 Einstein 的 u_n²，Hubbard 的权重应该是 |g(k)|²/DOS(k)，但 **g(k) 的符号在 BZ 上不恒定** ⟹ "diag(a_k) 的逐点非负性"论证**失效**。

**这不是技术细节，是结构性障碍**：
- Einstein 定理 I 的核心是"L1 的对角增长逐点压住 P 的负部分"，这要求**正权重 u_n>0**。
- d-wave **必须有节点**（配对对称性），节点处 g(k)=0 ⟹ 那些 k 态根本不参与配对 ⟹ 不存在"逐点压制"。

**可能的替代路径**（Phase 2 探索）：
1. **限制到 nodal pocket**：只看费米面上 g(k)>0 的区域，定义"有效维度降维"的算子 — 但这失去了全局陈述。
2. **换工具**：不用逐点论证，直接在 k 积分层面做**变分界** — 比如 Temple 不等式、Weyl 单调性，但这需要**更强的核信息**（比如 ‖P[f₋]‖ 的上界，而 Einstein 的优势恰恰是**只用单边界 ρ(f)**）。
3. **放弃免疫的一致性**：定理 I 的价值在于"λ₂ 界与 f₋ 的强度/形状**无关**"（推论 I.2 的免疫预算）。如果 Hubbard 的适配需要依赖 χ_spin 的具体形状，定理的**普适性**就丢了 — 虽然仍可能对特定参数给出界，但不再是"no-go 型"结论。

**时间估计**：如果走路径 1（nodal pocket），2 周；如果走路径 2（换工具），4–6 周（需要推导 Hubbard 的新引理）。

---

## 3. 推论适配

### 推论 I.1（静态库仑完全免疫）

**原陈述**：瞬时库仑排斥 P_C = 2μ̃ u_c u_cᵀ（秩一半正定）⟹ λ₂(M) ≤ μ₂(ℌ₊)。

**Hubbard 适配**：**可能直接套用**

**理由**：
- Hubbard 的库仑排斥 U 在 Bethe-Salpeter 方程里对应**瞬时相互作用**（Ω=0 的 Matsubara 频率），RPA/DMFT 处理后进入电荷通道 Φ_d。
- 如果电荷通道可以分离成**秩一或低秩**的形式（比如 q=0 的长程库仑 ∝ 1/q²，Fourier 到实空间是常数 ⟹ 秩一），推论 I.1 的 Cauchy 交错论证**逐字成立**。

**需验证**：Hubbard 的 Φ_d(q) 是否可以写成低秩形式。文献（[arXiv:cond-mat/0612579](https://arxiv.org/html/cond-mat/0612579v1) 图 6）显示电荷通道在 d-wave 投影下**弱排斥、几乎 q 独立** ⟹ 如果 Φ_d ≈ const，那就是秩一的。

**时间估计**：1 周（检查文献 + 数值验证秩结构）。

---

### 推论 I.2（免疫预算 c₁²）

**原陈述**：c₁² = |⟨ψ₁,e₀⟩|²（顶配对模在第一 Matsubara 频率上的权重）是免疫预算。

**Hubbard 适配**：**需重新定义 c₁**

**瓶颈**：
1. Einstein 的 e₀ 是"第一 Matsubara 基矢"（ω₁ = πT），c₁² 是**频率空间的浓度**。
2. Hubbard 的顶配对模是 Δ(k)，定义在**动量空间** — "第一 Matsubara"对应 ωₙ=πT，但**哪个 k**？

**可能的适配定义**：
- **费米面平均**：c₁² := ⟨|Δ(k_F,ω₁)|²⟩_{FS}，费米面上 ω₁ 权重的平均。
- **或者节点附近**：c₁² := ⟨|Δ(k,ω₁)|²⟩_{nodal pocket}，只在 g(k)≠0 的区域积分。
- **物理意义**：c₁² 控制"配对模在低频的集中度" — Hubbard 下这对应"准粒子峰 vs 高能 incoherent 谱重"的分配，与 Z 因子（质量增强 m*/m = 1/Z）相关。

**数值检验**（Track C 可先行）：
- QMC 数据通常报 Tc(U,μ)，可以直接与 Einstein 的 μ₁(ϖ) 表（定理 I §2）对比 — 如果 Hubbard 的 Tc 压制份额 s = 1-1/μ₁ 也落在某个"c₁²(U)"的范围内，**免疫预算的定性图像就成立**，即使解析定义还没对齐。

**时间估计**：2 周（定义对齐 + 数值检验）。

---

## 4. 核心机制诊断

### 4.1 哪些机制**直接匹配**？

| 机制 | Einstein | Hubbard | 匹配度 |
|---|---|---|---|
| 算子自伴性 | ✓（实对称 M） | ✓（Hermitian K） | **完全匹配** |
| 秩一吸引 + 排斥扰动 | ✓（f₊ 主导，f₋ 小扰动） | ✓（磁通道主导，电荷通道次要） | **完全匹配** |
| Δ/Z 通道分离 | ✓（配对 vs 质量重整化） | ✓（粒子-粒子 vs 粒子-空穴） | **完全匹配** |
| 平坦 DOS | ✓（u_n ∝ n^{-1/2}） | **?**（2D Hubbard 有 van Hove 奇点） | **需适配**（见 `nambu_weight_audit`） |
| 顶模孤立性（μ₁ ≫ μ₂） | ✓（Perron–Frobenius） | ✓（d-wave 主导 s-wave/p-wave） | **完全匹配**（实验支持） |

**结论**：**物理图像完全对应** — Hubbard 就是 Einstein 算子的"真实物理实现"。

---

### 4.2 哪些工具**不适用**？

| 工具 | Einstein 的角色 | Hubbard 的障碍 |
|---|---|---|
| Bochner 正定性（引理 1） | f∈𝒫 ⟹ f̂(θ)≥0 | χ_spin(q) 在 q 空间**非正定** |
| Toeplitz+Hankel 下界（引理 2） | 反射对称 ⟹ T+H ⪰ c·I | 动量空间无 Hankel（n+m）型配对 |
| 逐点压制（引理 4） | u_n>0 ⟹ diag(a)≥0 逐点 | g(k) 变号 ⟹ 无法逐点压制 |

**结论**：定理 I 的**具体证明路线**严重依赖 Matsubara 频率的离散性 + Toeplitz/Hankel 结构，这些在动量空间**不直接搬运**。

---

### 4.3 什么是**新障碍**？

1. **d-wave 节点**：
   - Einstein 的秩一主导来自"吸引部分 f₊ 在所有 n 上逐点正"（u_n>0）⟹ ψ₁(n)>0 全正（Perron–Frobenius）。
   - Hubbard 的 d-wave **必须有节点** g(k)=0，那些 k 态不参与配对 ⟹ 顶模 Δ(k) 在节点处也为零 ⟹ **不存在"全 k 正"的配对模**。
   - 这导致引理 4 的逐点论证**结构性失效**。

2. **van Hove 奇点**（2D square lattice）：
   - Einstein 假设平坦 DOS ⟹ u_n = (2n+1)^{-1/2}。
   - 2D Hubbard 在 E=0（van Hove）处 DOS 对数发散 N(E)~ln|E| ⟹ 费米面附近态密度**不平坦**。
   - `nambu_weight_audit` 已证明：van Hove 型 DOS 使引理 4 的充分条件失效（a_n<0 出现），虽然**实际谱隙仍安全**（gap/λ₁≈1.00）。
   - Hubbard 适配**必须**处理非平坦 DOS — 好消息是 `nambu_weight_audit` 给了精确 L1 形式（B_n 卷积），坏消息是**没有解析充分条件**（只有数值证据）。

3. **动量依赖的 Z 因子**：
   - Einstein 的 Z 通道是对角 L1，只依赖 Matsubara 频率 n。
   - Hubbard 的 Z_n(k) 是**动量依赖**的准粒子权重（Self-Energy Σ(k,ω) 的导数），在费米面附近变化剧烈（pseudogap 开口 ⟹ Z→0）。
   - 这使 L1 不再是"简单对角"，可能有**非对角泄漏** — 虽然 `nambu_weight_audit` 的精确 B_n 已包含这个效应（g_m 因子），但解析控制更难。

---

## 5. Phase 1 可行性判断

### 判断：**Medium**（4–6 周，2 个引理需新证明，1 个根本障碍需换工具）

**可以做的**（2 周内）：
- ✓ 引理 3 直接套用
- ✓ 推论 I.1（库仑免疫）验证秩结构，可能逐字成立
- ✓ 推论 I.2（c₁² 预算）重新定义后，**数值检验可先行**（Track C）

**需要新证明的**（2–4 周）：
- 引理 1（Bochner）：放弃全局正定性，改用局部正性 / 数值验证 χ̂(x) 的符号分布
- 引理 2（T+H 下界）：明确 Hubbard 是否有 Hankel 项，如果没有则换下界工具（Schur 补 / Gershgorin）

**根本障碍**（4–6 周，可能需要 Phase 2）：
- **引理 4（秩一泄漏）**：d-wave 节点使逐点论证失效，必须换工具：
  - 路径 A：限制到 nodal pocket（丧失全局性）
  - 路径 B：变分界（需要更强的核信息，丧失一致性）
  - 路径 C：接受"数值证据级结论"（`nambu_weight_audit` 已证明一切所测 DOS 下 gap/λ₁≥0.976），放弃解析证明

**推荐策略**：
1. **Phase 1A（2 周）**：完成引理 1/2/3 + 推论 I.1 的适配，给出"部分解析 + 数值支持"的中间结论。
2. **并行启动 Track C（实验数据匹配）**：用 Hubbard QMC 的 Tc(U,μ) 数据直接检验"免疫预算"的定性图像 — **这不需要等引理 4 解决**，可以立即做。
3. **Phase 1B（4 周）**：如果 Track C 的数值检验**强支持**免疫机制，再投入引理 4 的换工具（路径 B/C）；如果 Track C **不支持**，说明 Hubbard 有定理 I 未覆盖的逃逸口（比如各向异性 / 多带），Phase 1B 停止，转 Track B（多带推广）。

---

## 6. 与仓库已有结论的关系

### 6.1 与 `nambu_weight_audit` 的接口

**直接可用**：
- L1 的精确形式 B_n（卷积 - 负频反射）已推导，**analytic** 级别
- van Hove / 赝能隙等 5 种 DOS 形状下，数值验证 gap/λ₁≥0.976（**numerical-evidence**）
- ⟹ Hubbard 适配的"非平坦 DOS"部分**已有现成工具**，不需从零开始

**遗留问题**：
- `nambu_weight_audit` §6.1 报"没有解析充分条件" — 这正是引理 4 适配的瓶颈
- 如果 Hubbard 适配走路径 C（接受数值证据），`nambu_weight_audit` 的 gap 表可以直接引用

---

### 6.2 与定理 E（Einstein 极端性）的关系

**定理 E**（`p4-theorem-status-2026-09-12`）：Einstein 核在 𝒫 类中**最大化** λ₁（给定 ϖ）。

**Hubbard 推论**：
- 如果 Hubbard 的 χ_spin(q) 可以"嵌入"到 𝒫 类的某个等价类（比如"局部 Bochner"），定理 E 给出 Tc 的**上界**：Hubbard 的 Tc ≤ Einstein 的 Tc（同 ϖ）。
- 这可以作为 **Track C 的交叉检验**：QMC 的 Tc(U) 是否被 Einstein 的 Tc(ϖ_eff(U)) 封顶？

---

### 6.3 与 P4D（μ*>0 库仑）的类比

**P4D 的教训**（`p4d-coulomb-step-model-downjumps-2026-09-13`）：
- 静态库仑 μ̃ 是**秩一半正定扰动**，推论 I.1 给出 λ₂≤μ₂ 的完全免疫。
- Hubbard 的 U（on-site 排斥）对应 P4D 的 μ̃ ⟹ **库仑免疫在 Hubbard 下也应该成立**（推论 I.1 优先验证）。

**预期**：Hubbard 的"U 增大 → Tc 下降"曲线应该与 P4D 的"μ̃ 增大 → Tc 下降（但 gap/λ₁≈常数）"定性一致 — Track C 可检验。

---

## 7. 建议的三轨并行策略

### Track A（本文档）：解析适配

**目标**：把定理 I 的引理 1–4 逐条适配到 Hubbard，给出解析级结论。
**时间线**：4–6 周（Phase 1A 2周 + Phase 1B 4周）
**交付物**：`hubbard_theorem_I_adapted.md`（带引理编号 + 证明 / "数值证据"标注）

### Track B（多带推广）：推广到 t-J / 三带 Hubbard

**目标**：如果单带 Hubbard 的适配遇到**结构性障碍**（d-wave 节点），检查多带模型（CuO₂ 三带、t-J）是否有"更接近 Einstein"的算子形式。
**前提**：Track C 的数值检验显示 Hubbard **不满足**免疫预算 ⟹ 说明定理 I 的假设（各向同性 / 单带）在 Hubbard 下被违反。
**时间线**：6–8 周（需要新推导）

### Track C（实验数据直接匹配）：**优先启动**

**目标**：用 Hubbard QMC 的 Tc(U,μ,T) 数据**直接检验**推论 I.2 的免疫预算图像，不等解析适配完成。
**具体步骤**：
1. 提取 QMC 数据：Tc(U) 在不同 U（on-site 排斥）下的 Tc 压制
2. 拟合 Einstein 的等效参数：ϖ_eff(U)、μ₁(U) = 1/Z(U)
3. 检验不等式：s(U) := 1 - Tc(U)/Tc(U=0) 是否 ≤ c₁²(ϖ_eff) （Einstein 表，定理 I §2）
4. 如果**成立** ⟹ 免疫机制在 Hubbard 下定性正确，Track A 的解析适配有物理保证
5. 如果**违反** ⟹ Hubbard 有定理 I 未覆盖的逃逸口，Track A 停止，转 Track B

**时间线**：1–2 周（数据提取 + 拟合 + 画图）
**交付物**：`hubbard_qmc_immunity_budget_check.md` + 图（s vs c₁²，多个 U 点）

**为什么 Track C 优先**：
- **不依赖引理 4** — 只需要推论 I.2 的**最终结论**（c₁² 预算），不需要引理 1–4 的逐条证明
- **2 周出结果** — 如果数值检验强支持，给 Track A 的 4–6 周投入以**物理正当性**；如果不支持，节省 Track A 的时间成本
- **已有数据** — Hubbard QMC 的 Tc(U,μ) 是高引用的标准 benchmark（Scalapino、White、最新 normalizing-flow 突破），不需自己跑

---

## 8. 诚实边界

1. **本文档未做任何数值计算** — 全部是从源文件（`delta_channel_immunity_2026-09-16.md`、`nambu_weight_audit_2026-09-16.md`）+ 文献（[arXiv:cond-mat/0612579](https://arxiv.org/html/cond-mat/0612579v1)、[arXiv:0706.0241](https://arxiv.org/html/0706.0241v1)）的机制对照，**未验证任何数值断言**。

2. **引理 4 的"根本障碍"判断可能过于悲观** — 也许存在某种"d-wave 对称性适配的 Perron–Frobenius"（比如在 irreducible representation 内部做秩一论证），本文档没有探索。如果 Track C 强支持免疫，值得再深挖。

3. **χ_spin(q) 的 Bochner 性质未查文献** — 本文档根据"RPA 有符号振荡"判断非正定，但没有引用具体的 Hubbard χ_spin 计算。可能存在某些参数区（比如弱耦合 U≪t）χ_spin **近似正定**，那时引理 1 可以直接套。

4. **Track C 的"免疫预算检验"依赖 c₁²(ϖ) 的数值表**（定理 I §2），那张表是 **numerical-evidence** 级别（float64 `eigvalsh`），不是 certificate。如果 Track C 要升格到 certificate 级，需要先把 c₁²(ϖ) 用区间算术重算 — 但这不影响定性检验（数量级对齐）。

---

## 9. 脚本与文件

- 本文档：`D:\physics\superconductivity\theory\notes\hubbard_bs_kernel_mechanism_adapt_2026-09-16.md`（**未提交**）
- 源文件：
  - `theory/notes/delta_channel_immunity_2026-09-16.md`（定理 I，analytic）
  - `theory/notes/nambu_weight_audit_2026-09-16.md`（一般 DOS 的 L1 精确形式 + van Hove 数值，analytic + numerical-evidence）
- 相关 memory：
  - `p4-theorem-status-2026-09-12.md`（定理 E / 引理 C / M，Conj. 2 状态）
  - `p4d-coulomb-step-model-downjumps-2026-09-13.md`（P4D 的库仑唯一性，推论 I.1 的物理实例）
  - `delta-channel-immunity-theorem-2026-09-16.md`（定理 I 的 memory 索引）

**未提交文件**：本文档 + 定理 I 源文件（在 `git status` 的 `??` 列表中）

---

## 参考文献

- [Scalapino et al., Phys. Rev. Lett. 2007 (arXiv:cond-mat/0612579)](https://arxiv.org/html/cond-mat/0612579v1) — 自旋涨落磁化率表示的配对相互作用
- [Maier et al., Phys. Rev. B 2007 (arXiv:0706.0241)](https://arxiv.org/html/0706.0241v1) — 2D Hubbard 的 d-wave 配对核系统分析
- [Antiferromagnetic correlations in doped Hubbard (arXiv:2405.04277)](https://arxiv.org/html/2405.04277v3) — 自旋磁化率与反铁磁关联
- [Magnetic Susceptibility in 2D Hubbard (arXiv:1904.10782)](https://arxiv.org/pdf/1904.10782.pdf) — QMC 的自旋磁化率计算

---

**状态**：机制对照完成，交付主会话。推荐：**Track C 优先启动**（QMC 数据检验），2 周出结果后再决定 Track A 的投入深度。
