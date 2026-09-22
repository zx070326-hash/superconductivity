# SCDFT 形式论框架（Essenberger 2014）

> Phase 1 交付物：理解 SCDFT 构造，特别是 f^xc 顶角修正如何让 Λ^SF 恢复 PSD。
> 基于 Essenberger et al., PRB 90, 214504 (2014) [arXiv:1409.7968]

## 核心问题

**已知**：裸 RPA 自旋磁化率 χ^s(q,ω) 在 FeSe 100% BZ 有负本征值（Team B 已验证）。

**断言**（Essenberger Eq.(37)）：SCDFT 有效相互作用 Im[Λ^SF_{ij}(ω)] > 0（逐元，含非对角）。

**目标**：理解 SCDFT 顶角修正 f^xc 如何"筛选掉"负本征值。

---

## 一、SCDFT 自能构造（Sec. II–IV）

### 1.1 多体微扰框架

SCDFT 从 Kohn-Sham 体系出发，自能 Σ_k 分解为：

```
Σ_k = GW贡献 + Σ^SF（自旋涨落）+ Σ^Ph（声子）
```

**关键简化**（decoupling approximation）：
```
Δ^xc_k ≈ δ_{k,k'} Δ^xc_k
```
即排除杂化效应，电子 KS 系统可解析对角化 → 配对势方程（SCDFT gap equation, Eq.(6)）：

```
Δ^xc_k = -Δ^xc_k Z^D_k - ∑_{k'} K^C_{kk'} tanh(βE_{k'}/2)/(2E_{k'}) Δ^xc_{k'}
```

其中 K^C_{kk'} 是**有效配对核**，包含三项贡献（库仑 w、自旋涨落 Λ^SF、声子 Λ^Ph）。

### 1.2 粒子-空穴传播子（Sec. III.B）

SCDFT 核心对象是**不可约粒子-空穴传播子** A^P（irreducible particle-hole propagator）：

**定义**（Eq.(13)）：
```
A_0(I,2,3,4) = δΣ^γ(I,2)/δG(3,4) = A^c_0 - A^d_0 =: A^{c-d}_0
```

其中：
- A^c：crossed（交叉图）
- A^d：direct（直接图）

**BSE 递推**（Eq.(14)–(15)）：
```
A^c = ∑_{n=0}^∞ A^c_{(n)} = ∑_{n=0}^∞ A^d_{(n)}
A^c_{(n+1)} = A^c_0 G G A^c_{(n)}
A^d_{(n+1)} = A^d_0 G G A^c_{(n)} + A^c_0 G G A^d_{(n)} - A^d_0 G G A^d_{(n)}
```

**封闭形式**（Eq.(16)）：
```
A^{c-d} = A_0 + A_0 G G A^{c-d}   with  A_0 = δΣ^γ/δG
```

### 1.3 符号约定的关键（Sec. III.C, p.5）

**反常部分的对称性**（collinear 磁系统）：

Eq.(19)–(20) 给出 c/d 贡献在不同自旋通道的符号：

```
A^{c-d}_{σ,-σ,-σ} = -A^d_{σσ,-σ,-σ} = -A^{c-d}_{σσ,-σ,-σ}     (19)
A^{c-d}_{σ,-σ,-σ} = A^c_{-σσ,-σ,-σ} = A^{c-d}_{-σσ,-σ,-σ}      (20)
```

**结论**：A^c 和 A^d 在自洽方程中**符号相反**（loop rule：crossed 无额外负号，direct 有）。

---

## 二、f^xc 核与夹心结构（Sec. IV–V）

### 2.1 局域近似（Sec. III.C）

**两点形式**（Eq.(25)–(26)）：
```
Σ^SF(I,2)_{ab} ≈ G(I,2)_{ab} w^SF(I,2)_{ab}
```

其中 w^SF 是有效相互作用，通过 f^xc 核构造：

```
w^SF(I,2)_{ab} = ∬ dx dx' {-G A^{c-d} F A^{c-d} / -F^† A^{c-d} G^† A^{c-d}}_{ab} / G(I,2)_{ab}
```

### 2.2 交换-关联核 f^xc（Sec. V）

**定义**（functional derivative）：
```
f^xc_{σ₁σ₂σ₅σ₆}(x₁,x₅) := δv^xc_{σ₁σ₂}(x₁) / δρ_{σ₅σ₆}(x₅)
```

**Sham-Schlüter 连接**：将 KS 密度映射到相互作用系统，f^xc 作为质量算符抵消 Hartree/xc 贡献的非对角部分。

**ALDA（adiabatic local density approximation）**：
- f^xc 取决于局域密度 ρ(r) 和自旋极化
- 在均匀电子气（r_s ~ 2–3，铁基密度）下可用 LDA 核（如 Perdew-Wang、von Barth-Hedin）

### 2.3 夹心形式（Eq.(32)）

**非磁系统的 Λ^SF**：

```
Λ^SF(x₁,x₂) := 2 ∬ dx dx' × f^xc_zz(x₁,x₂) χ_zz(x,x') f^xc_zz(x₁,x₂)
```

这是 **f^xc · χ_zz · f^xc** 的夹心（因子 2 来自自旋求和）。

**约化到 Nambu 空间**：
```
Λ^{(c-d)}_{σ₁σ₁σ₂σ₂}(x₁,x₂) = z_{σ₁} z_{σ₂} (1/2) Λ^SF(x₁,x₂)
Λ^{(c-d)}_{σ̄σσ̄σ}(x₁,x₂) = Λ^SF(x₁,x₂)
```

---

## 三、最终配对核与符号断言（Eq.(36)–(37)）

### 3.1 SCDFT gap equation（Eq.(36)）

```
Δ^xc_k = -Δ^xc_k Z^D_k - ∑_{k'} K^C_{kk'} tanh(βE_{k'}/2)/(2E_{k'}) Δ^xc_{k'}
```

**配对核**：
```
K^C_{kk'} = (2/π) ∫₀^∞ dω [Im[w_{kk'}(ω) + (3/2)Λ^SF_{kk'}(ω)] + Λ^Ph_{kk'}(ω)] / [tanh(βζ_k/2) tanh(βζ_{k'}/2)] J⁻(ζ_k,ζ_{k'},ω)
```

其中：
- w_{kk'}(ω)：库仑散射（裸 Coulomb 的屏蔽版本）
- Λ^SF_{kk'}(ω)：自旋涨落（本工作核心）
- Λ^Ph_{kk'}(ω)：声子（吸引）
- J⁻：Fermi-Dirac 因子组合

### 3.2 符号断言（Eq.(37)）

```
Im[Λ^SF_{kk'}(ω)] > 0  and  Im[w_{kk'}(ω)] > 0
```

**物理解释**（Sec. II, p.2 footnote）：
> "spin-fluctuations and Coulomb terms are **intrinsically repulsive** while phonons provide attraction (Λ^Ph < 0)"

**但问题在于**：
- 裸 RPA χ^s 有负本征值（Team B 已验证）
- Λ^SF 通过 f^xc 夹心构造自 χ^s
- **为什么 f^xc · χ_zz · f^xc 能"翻转"符号？**

---

## 四、机制假说（待 Phase 2/3 验证）

### 假说 α：f^xc 的非局域性平滑负本征值

**推理**：
- χ^s(q,ω) 的负本征值可能集中在**尖锐的 q 点**（如嵌套向量 Q_AFM）
- f^xc 是**空间非局域**的（虽然 ALDA 是局域近似，但 BSE 迭代引入非局域性）
- 夹心 ∬ f^xc(x₁,x) χ(x,x') f^xc(x',x₂) dx dx' 相当于**卷积平滑**

**检验**：Phase 2 逐 q 点追踪 χ^s(q) → Λ^SF(q) 的本征值谱变化。

### 假说 β：自洽迭代的 PSD 收敛

**推理**：
- SCDFT 自洽求解 Δ^xc，类似变分原理
- 物理解：系统只接受**降低自由能**的解 → PSD 相互作用
- BSE 迭代 A^{c-d} = A_0 + A_0 G G A^{c-d} 可能向 PSD 方向收敛

**检验**：Phase 2 监控自洽迭代中 Λ^SF 本征值的演化。

### 假说 γ：频率积分的因果律约束

**推理**：
- Eq.(36) 的频率积分 ∫₀^∞ dω [...] / tanh(...) 有特定权重
- 涨落-耗散定理：Im χ(ω) 的符号受 causality 约束
- 可能**负贡献在积分中相消**

**检验**：Phase 3 分解 Im Λ^SF(ω) 的频率依赖，看哪些 ω 段主导。

---

## 五、与 Team C 的整合点

**Team C 已建立**：涨落-耗散定理 + Q_AFM 嵌套机制 → χ^s 符号约束（文献已知）。

**Team B（本工作）需补充**：
1. **SCDFT 顶角如何在此框架下工作**（f^xc 的作用）
2. **幅度增强 + PSD 恢复**的定量机制
3. **合起来**：完整的 Essenberger Eq.(37) 机制

**交付物**：
- Phase 1（本文档）：理论框架 ✓
- Phase 2：SCDFT 代码 + FeSe 复现（2–3 周）
- Phase 3：机制提取 + 物理图像（1–2 周）

---

## 六、Phase 2 准备清单

### 需要的组件

1. **f_xc 核计算**：
   - libxc（ALDA/PGG 近似）
   - 或手动实现 von Barth-Hedin 核（r_s ~ 2–3 参数）

2. **自洽迭代**：
   - 解 SCDFT gap equation（Eq.(36)）
   - 监控 Δ^xc_k 和 Λ^SF_{kk'} 的收敛

3. **频率网格**：
   - Matsubara → 实轴解析延拓（Padé 或最大熵）
   - 或直接在实轴计算（需处理因果律）

4. **FeSe 参数**：
   - 5-轨道紧束缚模型（Team B 已有）
   - U（Hubbard）、J（Hund's）参数（文献：U ~ 1 eV, J ~ 0.25 eV）
   - 费米面拓扑（已知：电子/空穴口袋）

### 技术门槛与风险

**高风险**：
- SCDFT 代码环境搭建（1–2 周）
- f_xc 数值不稳定性（需高精度）

**降低风险**：
- 先用现有代码（octopus/exciting）复现 Essenberger 结果
- 再从零实现关键步骤（追踪符号结构）

---

## 七、文献优先权

**Essenberger 的新颖性**：
- **断言 Im[Λ^SF] > 0 是 SCDFT 结果**，不是裸 RPA
- **FeSe s± 配对的 SCDFT 解释**（非对角主导）

**前置工作**：
- Marques et al. PRB 72, 024546 (2005)：SCDFT 基础理论
- Oliveira et al. PRL 60, 2430 (1988)：Sham-Schlüter 连接

**竞争工作**：
- RPA-based spin-fluctuation（无 SCDFT 顶角）
- Eliashberg 理论（无 DFT 自洽）

---

## Status

- [x] Phase 1：理论框架学习（本文档）
- [ ] Phase 2：SCDFT 代码搭建 + FeSe 复现
- [ ] Phase 3：机制提取 + 与 Team C 整合

**下一步**：开始 Phase 2 — 选择代码路线（现有 vs 从零）+ 搭建环境。
