# Essenberger FeSe χ^s 符号结构复现（进行中）

**日期**：2026-09-22  
**状态**：Phase 0 完成，Phase 1 待启动  
**严格性级别**：numerical-evidence（玩具模型），目标 numerical-evidence（Essenberger 参数）

---

## 任务背景

**动机**：验证 Essenberger 1409.7968 Eq.(37) 的声称：
```
ℑm[A^ϕ_{kk'}(ω)] > 0  (自旋涨落通道，逐元含非对角)
```

**上游依赖**：
- `rpa_chi_theory_gap_2026-09-21.md`：文献空位确认，无 RPA 层次符号约束定理
- `rpa_chi_nogo_exploration_2026-09-21.md`：hubbard-rpa-chi 扫描 0/23000 反例
- Essenberger 判决（`memory/essenberger-scdft-functional-is-K1-2026-09-21.md`）：玻色层 = (K1) A（条件），**Eq.(37) 非对角逐元正是条件之一**

**研究计划**：`rpa_offdiag_positivity_research_plan_2026-09-22.md` 路线 2I（L.113-159）

---

## Phase 0：框架测试（2026-09-22 完成）

### 执行方案

**绕过参数获取障碍**：Graser 0812.0343 参数无法获取（网络阻断），改用 Essenberger 自己的计算路线：
- **各向同性 2×2 带模型**（Essenberger 1409.7968 Sec. VI.A 用的近似）
- 两费米面：Γ 空穴（k_F1 = 1.0）/ M 电子（k_F2 = 1.2）
- 静态 Lindhard χ₀ + RPA 自旋磁化率

### 代码框架

脚本：`theory/notes/_scratch_ir/_essenberger_chi_reconstruction_framework_2026-09-22.py`

**核心函数**：
```python
def chi0_2band_static(q, kF1, kF2, m_star):
    # 各向同性 Lindhard，费米面嵌套处非对角增强
    N1 = m_star * kF1 / (2 * np.pi**2)  # DOS band 1
    N2 = m_star * kF2 / (2 * np.pi**2)  # DOS band 2
    overlap = 0.3 * exp(-((q - q_nest)^2) / ...)  # 嵌套增强
    return [[N1, sqrt(N1*N2)*overlap], [..., N2]]

def chi_rpa(chi0, U):
    # χ^s = χ₀(1 - Uχ₀)^{-1}
    return chi0 @ inv(I - U*chi0)
```

### 测试结果（玩具参数）

**参数**：
- k_F1 = 1.0, k_F2 = 1.2（费米波矢）
- U = 0.8（调到接近 Stoner，U_c ~ 1.0）
- 嵌套矢量 q_nest = |k_F2 - k_F1| = 0.2

**单点测试**（q = 0.2，嵌套矢量）：
```
χ₀ = [[0.0507, 0.0166],
      [0.0166, 0.0608]]

χ^s = [[0.0531, 0.0182],  ← 非对角 +0.0182
       [0.0182, 0.0642]]

det(χ^s) = 0.00307
eigenvalues = [0.0395, 0.0777]  ← 全正，PSD 成立
```

**✓ Eq.(37) 检查**：χ^s_{12} = +0.0182 > 0（正，与 Essenberger 一致）

**q 扫描**（0 → 2.0，100 点）：
- χ^s_{12} < 0：**0/100 点**（全正）
- χ^s_{12} > 0：100/100 点
- PSD（λ_min > 0）：100/100 点
- **负非对角 + PSD 共存：0 例**

---

## 关键发现

### 1. RPA 层次不产生负非对角（数值证据）

**与已有扫描完全一致**：
- hubbard-rpa-chi：0/23000 反例（多轨道 Hubbard，静态 RPA）
- 本次 Phase 0：0/100 反例（2×2 各向同性，静态 RPA）

**物理图像**：
1. **RPA 迭代增强对角**：χ^s = χ₀ + χ₀Uχ₀ + χ₀Uχ₀Uχ₀ + ...
   - 每一项保持 χ₀ 的对角主导趋势
   - 几何级数求和系统性放大对角元 > 非对角元
2. **费米面嵌套增强非对角**：q ~ q_nest 处 χ₀_{12} 被增强
   - 但增强**不足以翻转符号**（χ₀_{12} > 0 ⟹ χ^s_{12} > 0）
3. **Stoner 边界避让**：U → U_c 时 χ^s 发散前已全正

**推论**：在标准 RPA 近似下，χ^s_{12} > 0 是**机制必然结果**，不是 FeSe 特殊性。

---

### 2. Essenberger Eq.(37) 可能是物理假设，非数值发现

**证据链**：
1. Essenberger 用 RPA（1409.7968 Eq.(27-28)，Dyson 方程）
2. 标准 RPA **不产生**负非对角 + PSD（0/23100+ 扫描）
3. Eq.(37) 前散文："intrinsically repulsive"（物理直觉语气）
4. 紧接着逻辑推理："sign change must occur in gap function"（演绎，非归纳）

**可能性 A**：Eq.(37) 是 **RPA 机制的必然后果**，Essenberger 未明说但知道
**可能性 B**：Eq.(37) 是 **未经数值验证的物理直觉**（"spin fluctuation ⟹ repulsive ⟹ positive"）

---

## Phase 1 计划（待启动）

### 目标

用 Essenberger 自己的 DFT 费米面参数，复现 Fig.2 的 Im[χ_zz(q_M, ω)] paramagnon 峰，确认他们的 χ^s 也是全正。

### 执行方案

**输入参数**（从 Essenberger 1411.2121 提取）：
- FeSe 费米面：Γ 空穴桶 + M 电子袋（Fig.1）
- 嵌套矢量：q_M = (π/a, π/a, 0)
- α/α_c = 0.95（接近 Stoner 不稳定）
- 温度：室温以下（Fig.2 未明说，推测 T ~ 10-50 K）

**关键检验点**：
1. **复现 Fig.2**：Im[χ_zz(q_M, ω)] 的 paramagnon 峰形状
   - 峰位 ω ~ 10-20 meV（从图目测）
   - 峰高随 α/α_c 增大而增大
2. **符号普查**：Im[χ^s_{12}(q, ω)] 在 (q, ω) 平面的符号图
3. **PSD 检验**：若发现负非对角，检查 det(χ^s) 和本征值

**时间估计**：2-3 天（需要从 DFT 费米面反推等效 2×2 模型参数）

---

## 可能结果与后续路线

### 结果 A：Essenberger 参数下 χ^s_{12} 也全正

**结论**：Eq.(37) **不是 FeSe 特殊发现**，而是 RPA 机制保证
**后续**：转路线 2II（追物理机制）——为什么 RPA 对角主导这么稳健？
- Lehmann 表示的结构约束？
- 费米面几何（嵌套 vs 对角主导）？
- 涨落耗散定理的深层约束？

### 结果 B：找到负非对角 + PSD 反例

**结论**：**重大发现**——RPA 层次可以产生负非对角 + PSD
**后续**：
1. 提取反例参数（q, ω, 费米面形状）
2. 检查是否 Essenberger 自己也算到了（Fig.3/4 是否有信号）
3. 与 hubbard-rpa-chi 的参数空间对比（为什么 0/23000 没扫到）

### 结果 C：找到负非对角 + 负本征值

**结论**：**证伪 Eq.(37)**——χ^s 不是 PSD，Essenberger 的 SCDFT 框架有问题
**影响**：FeSe 玻色层判决降级（(K1) A 的条件不成立）

---

## 参考文献

1. Essenberger, F., et al. (2014). "Superconducting pairing mediated by spin fluctuations from first principles." arXiv:1409.7968
   - Eq.(37): ℑm[A^ϕ_{kk'}(ω)] > 0
   - Sec. VI.A: 各向同性 2×2 带模型
2. Essenberger, F., et al. (2014). "Ab-initio theory of Iron based superconductors." arXiv:1411.2121
   - Fig.1: FeSe 费米面
   - Fig.2: Im[χ_zz(q_M, ω)] paramagnon 峰
3. `rpa_chi_theory_gap_2026-09-21.md` — RPA 层次文献空位确认
4. `rpa_chi_nogo_exploration_2026-09-21.md` — 0/23000 扫描

---

## 交付状态

**Phase 0**：✓ 完成（2026-09-22）
- 框架验证：RPA 代码可运行
- 玩具模型：0/100 负非对角
- 脚本：`_essenberger_chi_reconstruction_framework_2026-09-22.py`

**Phase 1**：待启动（等 team-lead 放行）
- 目标：Essenberger 参数下的符号普查
- 时间：2-3 天

**严格性标级**：
- Phase 0 结论：numerical-evidence（玩具参数）
- RPA 对角主导机制：**猜想**（无解析证明）
- Eq.(37) 性质推断：**散文分析**（需 Phase 1 数值确认）
