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

---

## 八、2026-09-22 追加：Eq.(37) 逐字核对——是断言，不是证明，也不是数值核验（对 team-lead 问题的回答）

**任务**：team-lead 要求核实其回忆"SF 片是 PSD 推断，非对角逐元正仍是散文"，需引用原文逐字核对。
两篇论文均已 `pdftotext -enc UTF-8 -layout` 重新提取（旧提取可能有乱码，本次未加 `-enc` 的旧文件
未使用，已用带编码参数的新提取 `_essenberger1409.7968.txt` / `essenberger1411.2121.txt`，均已 `git add -f`）。

### 8.1 Eq.(37) 在 PRB90/1409.7968（形式论文）里的原文（第 8 页，本笔记 §3 已转录框架，这里是逐字复核）

配对核 `K^C_kk'` 定义式（Eq. 36 上方）列出三项贡献后，原文写：

> "The two last terms originate both from the Coulomb interaction and are therefore intrinsically repulsive:
> Im[Λ^SF_kk'(ω)] > 0 and Im[w_kk'(ω)] > 0. (37)"

**这是全篇唯一一处给出 Eq.(37) 的地方，前后未见任何数学推导**（已对全文 `positiv|definite|repulsiv|sign of|
fluctuation-dissipation` 等关键词逐一检索，命中仅此一处及其后文重复"repulsive"措辞，无 PSD/半正定/
涨落-耗散定理等推导性词汇出现在 Λ^SF 相关段落）。即：

- **不是 team-lead 回忆的"PSD 推断"**——文中甚至没有给出"sandwich 结构 `Λ^SF=2∬f^xc χ_zz f^xc`（Eq.32）
  ⟹ 半正定 ⟹ 逐元正"这一步论证链；`Λ^SF` 的正定性从未被从 Eq.(32) 的结构推出。
- **是纯物理断言（散文），比"PSD 推断"还弱一档**：理由只有一句"两项都源自库仑相互作用，因此内禀排斥"，
  没有引用任何定理、没有给出符号证明，属于**物理直觉陈述**，不属于"assumed as a mathematical premise"
  或"proven"，也不是"observed"（下面 §8.2 说明它甚至不是从计算中观察到的）。

### 8.2 全篇（两文合计）从未数值核验 `Im Λ^SF` 的逐元符号

- **PRB90 第六节的"二能带各向同性模型"** (Eq. 38)：`Λ^SF_IJ(E,E',ω)` 被**手工构造**为
  `c1·N_J(E)·[1-((ω-ω̄)/c2)²]`（clip 到 0 以下截断）——`N_J(E)≥0` 且括号项在支撑区间内 `≥0`，
  该模型**按构造方式天然满足 Eq.(37)**，不是从 χ_zz 算出来再验证符号；这是一个**玩具谱函数**，
  不触及"从第一性原理 χ_zz 出发的夹心结构是否保号"这个问题。
- **1411.2121（FeSe 应用文）** 第 174-176 行同样是**近似方案**层面的记号：
  ```
  w_ij(E,E',ω) ≈ w_ij(E,E',0)     (repulsive)
  Λ^SF_ij(E,E',ω) ≈ Λ^SF_ij(ε_F,ε_F,ω)   (repulsive)
  ```
  "(repulsive)" 是**贴在近似方案旁的标签**，不是计算后核验的结果陈述。全文唯一给出的、与真实
  计算相关的定量描述是（第 197-203 行附近）："Coulomb pairing (Fig. 3b) is as expected,
  **diagonally dominated** (++ and −− components)" ——"diagonally dominated" 是比"逐元正"弱得多
  的陈述（对角占优不蕴含非对角同号，更不蕴含正号），且这里说的是 Coulomb 项 `w`，不是 `Λ^SF`。
  对 `Λ^SF` 本身，全文只报告了 `λ_ij` 矩阵的**最大本征值**（`λmax=0.07→0.48`，随 α/αc）和一句定性
  描述"SF are dominated by off diagonal components"（导致 s± 对称性），**从未报告 `Im Λ^SF_ij` 的
  逐元符号表或逐元符号是否处处一致**。Fig. 3c/4a 画的是 `Λ^SF` 作为 α/αc 或能量的函数（图像本身
  未 Read，仅核对文本描述；若需逐元符号的图面证据需另外 Read 页图，但文本层已确认叙述中不含
  "we verified elementwise positivity"一类陈述）。

### 8.3 结论（回答 team-lead 的问题）

**逐元正性（Eq.37）在两篇原文中都只是一句无推导的物理断言（"故内禀排斥"），既不是从 Eq.(32) 的
夹心结构证明出来的（没有 PSD 论证链），也不是从真实第一性原理 χ_zz 计算中逐元核验观察到的**
（应用文仅报告最大本征值和"对角占优"这类弱得多的统计量，从未报告逐元符号表）。比 team-lead
"PSD 推断，非对角逐元正是散文"的回忆更弱一档——**"PSD 推断"这一步论证本身也不存在，只有断言**。

**对本方向（动力学 `Im χ^s(q,ω>0)` 逐元符号扫描，见
[`../misc/dynamic_im_chis_elementwise_sign_scan_2026-09-22.md`](../misc/dynamic_im_chis_elementwise_sign_scan_2026-09-22.md)）
的意义**：我们的 320 点扫描（离对角 45.8% 负 / 54.2% 正，逐元混合）与 Essenberger Eq.(37) 的断言
**方向相反**，但这不构成对 Eq.(37) 的证伪——Eq.(37) 针对的是**裸 RPA 之外经 SCDFT 顶角 `f^xc` 夹心
构造的** `Λ^SF`，我们算的是裸 `χ^s`（无顶角修正）。可以确立的是：**Essenberger 本人从未在原文中
给出任何独立于顶角修正机制之外的证据表明逐元正性会自动成立**——若顶角修正确实是逐元正性的
唯一来源，那么我们的裸-RPA 逐元混合结果与 Eq.(37) 完全相容（顶角修正是必要的，不是装饰性的）；
若顶角修正被证明不足以翻转全部离对角符号，Eq.(37) 本身就会失去支持——但**这一步验证需要真正实现
`f^xc` 夹心结构**（PRB90 §V ALDA/local density 近似），不在本次任务预算内，超出团队当前 1 周
Direction 2I 预算（同 §4.2 已评估的"2-3 周"工作量）。

**严格性**：这是文本层核对（`analytic`，直接引用原文字符串，非数值计算）。引用行号来自
`pdftotext -enc UTF-8 -layout` 输出，已用两种仪器交叉：本节引文均为完整英文句子/公式编号，
未见明显 OCR/编码问题（无 `�`），未额外核对页图；若需要更高确信度可 Read 对应页图复核，
当前未做（预算未要求逐字符级确信度，文本层已足够回答"证明/假设/观察"这一判词问题）。
