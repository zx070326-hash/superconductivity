# RPA 磁化率 no-go 定理探索

**日期**：2026-09-21  
**状态**：探索完成，部分 analytic / 部分 numerical-evidence  
**动机**：hubbard-rpa-chi 扫描 23000+ 组参数，0 个反例（χ_12<0 且 χ⪰0 同时成立）

---

## 问题陈述

**目标**：证明或证伪 RPA 层次的结构性约束：

> 若 RPA 自旋磁化率 `χ^RPA = χ₀[1−Uχ₀]⁻¹` 满足非对角元 `χ_ij < 0`（某 i≠j），  
> 则 χ^RPA **必有**负本征值（不 PSD）。

**背景**：
- Essenberger scdft-functional Eq.(37) 的自旋涨落层 χ_SF 可能有负非对角元
- 若 χ_SF 在 RPA 近似下计算，no-go 定理将限制其 PSD 性
- hubbard-rpa-chi 的 0/23000+ 扫描提示存在系统性约束

---

## 主要发现

### 1. Lehmann 表示约束（analytic）

**静态磁化率的谱表示**[^1]：
```
χ^RPA_αβ(q,ω) = ∫ dE ρ_αβ(q,E) / (E - ω - iη)
```
在 ω=0 静态极限：
```
χ^RPA_αβ(q,0) = ∫ dE ρ_αβ(q,E) / E
```

**谱函数性质**（从反射正定性 + 幺正性）：
- 对角：`ρ_αα(q,E) ≥ 0` 严格（单粒子激发谱权重非负）
- 非对角：`ρ_αβ(q,E)` 符号由相干性决定，**可以为负**

**热平衡约束**[^2]：涨落耗散定理保证 χ 的 Hermitian 部分在平衡态必须 PSD。  
对于实对称 χ（时间反演不变）⟹ **χ ⪰ 0 是物理要求**。

**2×2 矩阵 PSD 条件**（analytic）：
1. 对角非负：`χ_11, χ_22 ≥ 0`
2. Cauchy-Schwarz：`χ_12² ≤ χ_11 · χ_22`

**推论**：χ_12 < 0 与 χ ⪰ 0 共存 ⟺ `|χ_12| / √(χ_11·χ_22) < 1`（相关系数 < 1）

---

### 2. RPA Dyson 级数结构（分析性论证）

**RPA 展开**：
```
χ^RPA = χ₀ + χ₀Uχ₀ + χ₀Uχ₀Uχ₀ + ...
       = χ₀ · (1 + Uχ₀ + (Uχ₀)² + ...)
```

**观察**（numerical-evidence 支持的猜想）：
- Lindhard χ₀ 在费米面附近通常**对角主导**[^3]（`|χ₀_αα| ≫ |χ₀_αβ|` for α≠β）
- RPA 迭代**系统性增强对角元**（几何级数求和）
- 非对角元增强较弱 ⟹ **比值** `|χ_12| / √(χ_11·χ_22)` **收缩**

**结果**：即使 χ₀_12 < 0，RPA 修正使得 PSD 条件**更容易满足**，不是更难。

---

### 3. Stoner 不稳定性边界分析[^4]

**Stoner 判据**：铁磁不稳定发生在 `det[1 − Uχ₀] = 0`，此时 χ^RPA 发散。

**发散前夕的矩阵结构**：
- χ^RPA 的最小本征值 λ_min → 0⁺
- 若同时 χ_12 < 0 显著 ⟹ **det(χ) = χ_11·χ_22 - χ_12² 接近零**
- 这要求 `χ_12² ≈ χ_11·χ_22`（接近 Cauchy-Schwarz 等号）

**数值观察**（from hubbard-rpa-chi）：
- 扫描从未触及这个边界配置（χ_12<0 且 det≈0 且 λ_min>0）
- 提示 RPA 动力学**避开**这个精细调节的参数区域

---

### 4. 扫描结果统计（numerical-evidence）

**hubbard-rpa-chi 两份脚本**：
- **1D 扫描**：两轨道 Hubbard，q=π，n∈[0.5,0.9]，U/t∈[1,4]，Δ/t∈[0.5,2]  
  总计 **~9×13×7 = 819** 组参数，**0 个反例**

- **2D 扫描**：方形晶格，各向异性 hopping t2/t1∈{0.5,0.75,1.0,1.5}，三个 q 点  
  总计 **~23000** 组参数，**0 个反例**

**结论（numerical-evidence 层）**：在 RPA 层次、静态极限 ω=0 下，  
**χ_12<0 且 χ⪰0 的配置在常规参数区域事实上不存在**。

---

### 5. 文献空位

**已检索**：
- Lehmann representation & reflection positivity[^1]
- Fluctuation-dissipation theorem[^2]
- RPA in Hubbard model[^5]
- Multiorbital susceptibility[^6]
- Stoner criterion & magnetic instabilities[^4]

**未找到**：
- 明确的解析 no-go 定理："RPA + 负非对角 ⟹ 必有负本征值"
- 最接近：RPA 正半定稳定矩阵下的解结构[^7]，但讨论的是**输入** χ₀ 的性质，非输出 χ^RPA

---

## 理论图景总结

### analytic 层结论：
1. **PSD 必要条件**：`χ_12² ≤ χ_11·χ_22`（2×2 矩阵）
2. **Lehmann 表示保证对角元 ≥ 0**（单粒子谱）
3. **涨落耗散定理要求 χ ⪰ 0**（热平衡）

### numerical-evidence 层结论：
- **0/23000+ 扫描**：静态 RPA 下，χ_12<0 与 χ⪰0 **事实上不共存**
- 提示存在**结构性约束**，但缺少解析证明

### 机制猜想（未严格证明）：
- RPA Dyson 级数保持 χ₀ 的对角主导性
- 几何级数求和使对角增强 > 非对角增强
- ⟹ 相关系数 `|χ_12|/√(χ_11·χ_22)` 系统性收缩，远离 Cauchy-Schwarz 等号

---

## 对 Essenberger Eq.(37) 的启示

**若 χ_SF 非对角为负且矩阵 PSD**：
1. **在 RPA 层次极难实现**（0/23000+ 数值证据）
2. 需要**超越 RPA** 的机制：
   - 顶角修正（ladder/bubble beyond RPA）
   - 自能反馈（self-consistent Green's function）
   - 动态屏蔽效应（ω 依赖）

**建议判决**：
- 若 Essenberger 的 χ_SF 声称在 **RPA 近似**下计算 ⟹ 负非对角 + PSD 需要**额外论证**
- 更可能情况：χ_SF 来自非 RPA 图，或矩阵只是 Hermitian 不保证 PSD

---

## 未完成方向（留作后续）

1. **动态响应 ω≠0**：
   - hubbard-rpa-chi 只扫了静态极限
   - Im χ(ω) 的正定性（能量耗散）可能进一步约束 Re χ

2. **解析证明尝试**：
   - 对 χ₀ 施加物理约束（费米面几何、粒子-空穴对称性）
   - 证明 RPA 变换保持某个符号不变量

3. **文献深挖**：
   - 凝聚态场论教材（Mahan, Altland-Simons）关于响应函数正定性
   - 磁性理论（Moriya, Villain）关于 RPA 磁化率的约束

---

## 参考文献

[^1]: [Lehmann representation and reflection positivity](https://arxiv.org/abs/1201.3415)  
[^2]: [Fluctuation-dissipation theorem modern perspective](https://arxiv.org/abs/2312.10134)  
[^3]: [Lindhard susceptibility in extended momentum space](https://arxiv.org/abs/1208.3210)  
[^4]: [Generalized Stoner criterion](https://arxiv.org/abs/1708.01971)  
[^5]: [Dyson equation RPA in Hubbard model](https://arxiv.org/abs/cond-mat/9804011)  
[^6]: [Multiorbital spin susceptibility](https://arxiv.org/abs/1102.5532)  
[^7]: [RPA solutions for positive-semidefinite stability matrix](https://arxiv.org/abs/1606.03167)

---

**严格性标级**：
- 2×2 PSD 条件：`analytic`
- Lehmann 表示对角非负：`analytic`（量子统计标准结果）
- 0/23000+ 扫描：`numerical-evidence`
- RPA 对角主导机制：**猜想**（无证明）

**交付物状态**：探索完成，无完整 no-go 定理，但有强数值证据 + 部分解析约束。
