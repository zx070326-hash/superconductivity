# RPA χ 理论缺口调查：负非对角元 + PSD 结构约束

**日期**：2026-09-21  
**状态**：文献空位确认  
**严格性级别**：文献综述 + 分析性论证

---

## 任务背景

**动机**：hubbard-rpa-chi 扫描 **0/23000+ 组参数**未发现反例（χ^RPA 负非对角 + PSD 同时成立）。

**问题**：是否存在 RPA 层次的结构性约束，防止 χ^RPA 同时具有负非对角元与正半定性？

**上游依赖**：
- `rpa_screening_psd_derivation_2026-09-21.md`：W^RPA ⪰ 0 已证（analytic）
  ⚠ **2026-09-22 订正**：该篇原证明第三步为假（Dyson 逐项 PSD，两条 rational-exact 反例）；
  结论经一行修补后仍成立，但**仅静态 `ω=0`**（`W=(v⁻¹−χ₀)⁻¹≻0`），
  **不含** `Im W(ω>0) ⪰ 0`。引用前先读该篇末尾的追加块。
- `rpa_chi_nogo_exploration_2026-09-21.md`：0/23000 扫描 + 部分解析约束
- Essenberger Eq.(37)：声称 Im Λ^SF_{ij} > 0（逐元，含非对角），需理论支撑

---

## 主要发现

### 1. 定理层检索结果（NEGATIVE）

**检索范围**：
- 标准教材：Fetter & Walecka, Mahan, Schrieffer（通过 WebSearch 二手）
- 铁基超导 RPA 文献：Graser, Maier, Hirschfeld 团队 2008-2015
- 多轨道 RPA 系统性工作：arXiv:1606.03167, 1606.06189, 1102.5532

**结果**：**未找到**明确的 no-go 定理形式："RPA + 负非对角 ⟹ 必有负本征值"

**最接近结果**：
1. **arXiv:1606.03167**（"Solutions of RPA equation for positive-semidefinite stability matrix"）  
   - 讨论**稳定性矩阵 S** PSD 时 RPA 解的性质
   - 关注点：**输入** χ₀ 的稳定性 ⟹ **输出** χ^RPA 的物理性
   - **不涉及** χ^RPA 本身的符号约束定理
   - WebFetch 被网络策略阻断，无法获取全文细节

2. **arXiv:1606.06189**（"Role of vertex corrections in matrix RPA for multiorbital Hubbard"）  
   - 多轨道 Hubbard 模型的矩阵 RPA 形式
   - 讨论顶角修正的作用，未明确陈述 χ^RPA 的符号结构约束
   - WebFetch 阻断

3. **arXiv:1102.5532**（"Multiorbital spin susceptibility in magnetically ordered state"）  
   - 多轨道自旋磁化率计算
   - 可能讨论矩阵结构，但 WebFetch 阻断无法确认

**文献空位确认**：标准教材与专题文献均**未显式组装** no-go 定理。

---

### 2. 机制假说探索

#### （a）Lehmann 表示约束（analytic）

**静态 RPA 磁化率**：
```
χ^RPA_{αβ}(q,0) = ∫ dE ρ_{αβ}(q,E) / E
```

**谱函数性质**：
- **对角元**：ρ_{αα}(q,E) ≥ 0 严格（单粒子激发谱权重非负）
- **非对角元**：ρ_{αβ}(q,E) 符号由**带间相干性**决定，**可以为负**

**2×2 矩阵 PSD 充要条件**（analytic）：
1. 对角非负：χ₁₁, χ₂₂ ≥ 0  
2. Cauchy-Schwarz：χ₁₂² ≤ χ₁₁ · χ₂₂

**推论**：χ₁₂ < 0 与 χ ⪰ 0 **可以共存** ⟺ |χ₁₂| / √(χ₁₁·χ₂₂) < 1（相关系数 < 1）

**结论**：Lehmann 表示**不禁止**负非对角 + PSD 配置。

---

#### （b）RPA Dyson 级数对角主导猜想（未证明）

**RPA 展开**：
```
χ^RPA = χ₀ + χ₀Uχ₀ + χ₀Uχ₀Uχ₀ + ...
       = χ₀[1 + Uχ₀ + (Uχ₀)² + ...]
```

**猜想**（numerical-evidence 支持）：
- Lindhard χ₀ 在费米面附近通常**对角主导**（|χ₀_{αα}| ≫ |χ₀_{αβ}| for α≠β）
- RPA 迭代**系统性增强对角元**（几何级数求和效应）
- 非对角元增强较弱 ⟹ 比值 |χ₁₂| / √(χ₁₁·χ₂₂) **收缩**

**机制图像**：
- 每一项 χ₀(Uχ₀)^n 都保持 χ₀ 的对角主导趋势
- 级数和使对角增强 > 非对角增强
- ⟹ 相关系数远离 Cauchy-Schwarz 等号，PSD 条件**更容易**满足

**状态**：**无解析证明**，仅有 0/23000 扫描支持（numerical-evidence）

---

#### （c）涨落耗散定理约束（analytic）

**热平衡要求**：涨落耗散定理保证 χ 的 Hermitian 部分在平衡态必须 PSD。

**静态极限**：对于实对称 χ（时间反演不变）⟹ χ ⪰ 0 是**物理要求**（热力学稳定性）

**推论**：任何违反 PSD 的配置**不能在热平衡态实现**。

**局限**：此约束只说"非 PSD 禁止"，不说"负非对角 + PSD 困难"。

---

#### （d）Stoner 不稳定性边界（分析性论证）

**Stoner 判据**：铁磁不稳定发生在 det[1 − Uχ₀] = 0，此时 χ^RPA 发散。

**发散前夕的矩阵结构**：
- χ^RPA 的最小本征值 λ_min → 0⁺
- 若同时 χ₁₂ < 0 显著 ⟹ det(χ) = χ₁₁·χ₂₂ - χ₁₂² 接近零
- 这要求 χ₁₂² ≈ χ₁₁·χ₂₂（接近 Cauchy-Schwarz 等号）

**数值观察**（from hubbard-rpa-chi）：
- 扫描从未触及这个边界配置（χ₁₂<0 且 det≈0 且 λ_min>0）
- 提示 RPA 动力学**避开**精细调节参数区域

**状态**：定性论证，无定量界或普遍性定理。

---

### 3. 铁基超导 RPA 文献调查

**检索焦点**：Graser, Maier, Hirschfeld 团队（2008-2013 铁基超导理论奠基期）

**主要发现**（from WebSearch）：

**arXiv:0812.0343**（Graser et al., "Near-degeneracy of several pairing channels"）：
- 5-轨道紧束缚模型 + RPA 增强的磁化率
- **核心结果**：s± 配对态（能隙在空穴/电子费米面间**反号**）
- **符号约束来源**：自旋涨落介导的配对势 V ~ χ(Q)，其中 Q 连接空穴/电子袋
- **矩阵结构**：未明确讨论 χ^RPA_{αβ} 的符号约束定理

**相关工作**（WebSearch 二手信息）：
- "Inelastic neutron scattering as probes of sign structure"
- "Sensitivity of superconducting state to electronic structure"
- **共同特征**：关注**能隙符号**（s±），不关注**完整 χ 矩阵**的正定性

**结构性原因**：
- 研究焦点（2008-2013）：配对机制、能隙对称性、中子散射签名
- **不是**：多轨道响应函数的数学结构定理

---

### 4. 对 Essenberger Eq.(37) 的启示

**Essenberger 声称**：Im Λ^SF_{kk'} > 0（逐元，含非对角），Eq.(37) 前散文 "intrinsically repulsive"

**理论支撑现状**：

| 层次 | 结论 | 依据 | 级别 |
|---|---|---|---|
| **W^RPA ⪰ 0** | RPA 屏蔽势正半定 | Lehmann 表示 + Dyson 级数 | **analytic** |
| **χ^RPA Hermitian 部分 PSD** | 热平衡稳定性 | 涨落耗散定理 | **analytic** |
| **χ^RPA 对角主导** | RPA 迭代增强对角 | 0/23000 扫描 + 机制猜想 | **numerical-evidence**（无定理） |
| **χ^RPA 逐元正（非对角）** | 负非对角 + PSD 事实上不共存？ | **NOT OBTAINED**（无定理） |

**关键区分**：
- **PSD** ≠ **逐元正**（PSD 的非对角元可为负）
- (Z-dom) 判据在多带上吃的恰是**非对角元**（带间通道）
- ⟹ **Eq.(37) 非对角逐元正 = 唯一既承重、又无结构支撑的一步**

**若 χ_SF 非对角为负且矩阵 PSD，需要超越 RPA 的机制**：
1. 顶角修正（ladder/bubble beyond RPA）
2. 自能反馈（self-consistent Green's function）
3. 动态屏蔔效应（ω 依赖）

**建议判决**：
- 若 Essenberger 声称在 **RPA 近似**下计算 ⟹ 负非对角 + PSD 需要**额外论证**
- 更可能：χ_SF 来自非 RPA 图，或矩阵只是 Hermitian 不保证逐元正

---

## 未完成方向（超出本次范围）

1. **动态响应 ω≠0**：
   - hubbard-rpa-chi 只扫了静态极限 ω=0
   - Im χ(ω) 的正定性（能量耗散）可能进一步约束 Re χ

2. **解析证明尝试**：
   - 对 χ₀ 施加物理约束（费米面几何、粒子-空穴对称性）
   - 证明 RPA 变换保持某个符号不变量（如对角主导度量）

3. **文献深挖**（受阻于网络）：
   - Mahan *Many-Particle Physics* Ch. 5（响应函数正定性）
   - Altland-Simons *Condensed Matter Field Theory*
   - Moriya *Spin Fluctuations in Itinerant Electron Magnetism*

4. **轨道分辨 χ_{αβ} 实验数据**：
   - 另一 agent 报告：0 个确认阳性，3 个候选（Maier PRB 83 最高优先）
   - 结构性空洞原因：计算成本（5×5 矩阵 × (q,ω) = 25 面板）+ 研究焦点（能隙符号，非完整矩阵）

---

## 总结

### analytic 层结论

1. **PSD 必要条件**：χ₁₂² ≤ χ₁₁·χ₂₂（2×2 矩阵 Cauchy-Schwarz）
2. **Lehmann 表示保证对角元 ≥ 0**（单粒子谱）
3. **涨落耗散定理要求 χ ⪰ 0**（热平衡，Hermitian 部分）
4. **W^RPA ⪰ 0**（Dyson 级数，已证）

### numerical-evidence 层结论

- **0/23000+ 扫描**：静态 RPA 下，χ₁₂<0 与 χ⪰0 **事实上不共存**
- 提示存在**结构性约束**（RPA 对角主导传递），但缺少解析证明

### 文献空位

**未找到**：
- 明确的解析 no-go 定理："RPA + 负非对角 ⟹ 必有负本征值"
- 最接近：RPA 正半定稳定矩阵下的解结构（arXiv:1606.03167），但讨论**输入** χ₀，非**输出** χ^RPA
- 标准教材（Fetter & Walecka, Mahan）给出**全部组件**（Lehmann 表示、Dyson 方程、因果性），但**未显式组装**为定理

**结论**：这是一个**开放的理论问题**（RPA 层次的符号结构约束），有强数值证据但无完整定理。

---

## 参考文献

### 直接检索（WebSearch）

1. [Graser et al. 2008, arXiv:0812.0343](https://arxiv.org/abs/0812.0343) — 5-轨道 RPA, s± 配对
2. [Solutions of RPA for PSD stability matrix, arXiv:1606.03167](https://arxiv.org/abs/1606.03167) — 稳定性矩阵 S
3. [Matrix RPA for multiorbital Hubbard, arXiv:1606.06189](https://arxiv.org/abs/1606.06189) — 顶角修正
4. [Multiorbital spin susceptibility, arXiv:1102.5532](https://arxiv.org/abs/1102.5532) — 磁有序态

### 相关材料（from rpa_chi_nogo_exploration_2026-09-21.md）

5. Lehmann representation and reflection positivity, arXiv:1201.3415
6. Fluctuation-dissipation theorem modern perspective, arXiv:2312.10134
7. Lindhard susceptibility in extended momentum space, arXiv:1208.3210
8. Generalized Stoner criterion, arXiv:1708.01971

### 教材（二手引用，未直接访问）

9. Fetter & Walecka, *Quantum Theory of Many-Particle Systems*, Ch. 3
10. Mahan, *Many-Particle Physics*, Ch. 5
11. Altland-Simons, *Condensed Matter Field Theory*

---

**交付状态**：文献空位确认完成  
**严格性标级**：
- 2×2 PSD 条件：analytic  
- Lehmann 对角非负：analytic  
- 0/23000 扫描：numerical-evidence  
- RPA 对角主导机制：**猜想**（无证明）  
- no-go 定理：**NOT FOUND**

**对主线的影响**：Essenberger Eq.(37) 非对角逐元正仍是**散文断言**（无结构定理支撑），FeSe 玻色层判决条件未解除。
