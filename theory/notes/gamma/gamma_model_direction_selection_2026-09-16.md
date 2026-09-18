# 新方向选择：γ-model + Coulomb μ* 严格 Tc 界

2026-09-16，分支 `p4-conj2`，**未提交**。

---

## 0. 方向选择过程

从 P4（Einstein 声子模型严格 Tc 界）出发，寻找 impact ≥ 2（PRL 级）的新物理方向。
调研了 5 个候选，经并行可行性验证后选定唯一存活方向。

| 候选 | 判定 | 关键原因 |
|---|---|---|
| #1 量子自旋谱隙 (Haldane) | **NO-GO** | LTI-SDP 硬性要求 frustration-free；Heisenberg 链不满足（§1 详析） |
| **#2 γ-model Tc + Coulomb μ\*** | **GO** | KAY/EA 明确放弃 μ\*；P4-D Coulomb 机器直接搬 |
| #3 BCS 冷原子 | **NO-GO** | 接触势无算子谱结构；FHS 占领；无声子无 P4 |
| 超级：Haldane Davis-Kahan | **NO-GO** | 扰动/隙 = O(L)→∞ (§1) |
| Einstein 延续（KAY 论文） | 低优先 | 已由 paper draft agent 完成初稿 |

---

## 1. Haldane 方向为什么走不通（完整技术分析）

### 1.1 Davis-Kahan 路线

代数恒等式（已验证，err < 1e-14）：
$$H_{\rm Heis} = 2H_{\rm AKLT} - \sum_i P_i^{(0)} - LI$$

每键扰动 $V_i = -P_i^{(0)}$，$\|V_i\| = 1$。
但 $V = \sum_{i=1}^L V_i$ 是**广延量**：$\|V\| \propto L$。

AKLT 谱隙 $\Delta_{\rm AKLT} \approx 0.350$（O(1)）。
$$\frac{\|V\|}{\Delta_{\rm AKLT}} \sim \frac{L}{0.350} \to \infty$$

Davis-Kahan 需要 $\|V\| < \Delta$，这在**原理上**对热力学极限不可能成立。
不是数值不够——是方法论层面的不匹配。

优化 $c$ in $\|H_{\rm Heis} - c H_{\rm AKLT}\|$ 也无帮助：S=0 扇区的特征值 $-2$ 固定，
每键扰动范数 ≥ 2 对任何 $c$。

### 1.2 LTI-SDP / Knabe / Gosset-Mozgunov 路线

**关键定义域限制**：全部要求 frustration-free（基态是所有局域项的公共零向量）。

**frustration-free 系统**：$E_0 = 0$（精确），谱隙 $\Delta = E_1$。
只需证一个数 > 0。

**Heisenberg 链**：$E_0 \approx -1.401L$，$E_1 \approx -1.401L + 0.411$。
谱隙是两个 O(L) 量的差：$\Delta = E_1 - E_0$。
需要相对精度 $O(1/L)$——无穷链极限需要**非局域全局论证**。

所有已知谱隙证明工具（Knabe/GM/LTI-SDP/detectability lemma）都是
局域→全局论证，frustration-freeness 是粘合步骤的本质前提。

Rai et al. 2026（arXiv:2411.03680）的 LTI-SDP 论文只测了
AKLT/Z₃-clock/Glauber，**从未测过 Heisenberg**——不是遗漏，是方法不适用。

### 1.3 差了什么

一个**不依赖 frustration-freeness 的、能从有限尺寸精确信息推到热力学极限的严格框架**。
43 年没人做到。这不在我们的工具箱内。

### 1.4 插值路径数值

| L | gap(Heis) | gap(AKLT) |
|---|---|---|
| 4 | 1.000 | 0.736 |
| 6 | 0.721 | 0.602 |
| 8 | 0.594 | 0.499 |
| 10 | 0.525 | 0.448 |
| ∞ (DMRG) | ~0.4105 | ~0.350 |

插值 H(t) = (1-t)·2H_AKLT + t·H_Heis 的谱隙全程 > 0（L=4,6,8 验证）。
物理上隙确实存在——但**证明它**需要超出当前所有工具的方法。

---

## 2. γ-model Phase 0 结果

### 2.1 KAY 的 𝔊(γ) 算子

𝔊(γ) = −𝔊₁ + 𝔊₂ + 𝔊₃，关键对角项 𝔊₁_{nn} = [1/(2n+1)] Σ_{k=1}^n 2/k^γ
（求和到 n，不到 ∞——γ≤1 时有限）。

γ=1 (SDW QCP) 校准结果：
- 𝔤(1) ≈ 1.6075（N=600 收敛，单调递增下界序列）
- c₁² = 0.723（强秩一主导）
- λ₁ − λ₂ = 1.58（巨大谱隔离）
- N=4 闭式 𝔤⁽⁴⁾(1) = 1.5398（与 KAY eq.31 一致）

### 2.2 Coulomb 发散论证

KAY 脚注 4："summations over all Matsubara frequencies would diverge"。

验证：$\sum_{m=0}^{N_c-1} u_m = \sum 1/\sqrt{2m+1} \sim \sqrt{N_c} \to \infty$。
发散**是真的**——常数 μ\* 加到所有 Matsubara 项确实发散。

但**物理正规化（固定 N_c 截断）数学上良定义**：
$$\mathfrak{G}(\gamma, \mu^*, N_c) = \mathfrak{G}(\gamma) - \mu^* P_{N_c}$$
其中 $P_{N_c}$ 是投影到 n < N_c 的有限秩算子。
有限秩扰动不破坏紧凑性，算子良定义。

### 2.3 Coulomb 效应（γ=1, N=400）

| μ\* | N_c | Tc/Tc₀ | c₁² | 谱隙 |
|---|---|---|---|---|
| 0.00 | — | 1.000 | 0.723 | 1.58 |
| 0.10 | 50 | 0.880 | 0.750 | 1.41 |
| 0.13 | 50 | 0.848 | 0.756 | 1.36 |
| 0.15 | 50 | 0.827 | 0.760 | 1.32 |

12–17% Tc 压制，物理合理（BCS 理论中 μ\*=0.13 典型 ~10–20%）。

---

## 3. 竞争者定位

| 组 | 成果 | 我们的优势 |
|---|---|---|
| KAY (J. Stat. Phys. 2025) | N≤4 代数闭式下界；N>4 浮点 | Certificate + Coulomb |
| Elezaby-Abanov (PRB 2026) | zeta 级数解析上界（更紧） | Coulomb 推广 |
| 无人 | γ-model + μ\* 的严格结果 | **我们是第一个** |

---

## 4. Phase 1 计划

1. iv Rayleigh 商下界（γ=1, N=200,400）→ certificate 级 𝔤⁽ᴺ⁾
2. Schur 尾部 ‖C‖ 和 ‖𝔊^BB‖ 上界
3. Davis-Kahan sinΘ → c₁² certificate
4. Coulomb μ\*=0.10/0.13/0.15 的 Tc(g, γ, μ\*, N_c) certificate 界
5. 交叉验证：KAY N=4 闭式 vs 我们的 certificate

---

## 5. 严格性分级

| 陈述 | 级别 |
|---|---|
| 代数恒等式 H_Heis = 2H_AKLT - ΣP⁽⁰⁾ - LI | numerical-evidence (float64, err<1e-14) |
| AKLT→Heisenberg 扰动比值发散 | analytic（广延量 vs O(1) 隙） |
| LTI-SDP 要求 frustration-free | analytic（方法定义域） |
| 𝔊₁ 对角 = Σ_{k=1}^n 2/k^γ / (2n+1) | analytic（KAY 定义） |
| 𝔤(1) ≈ 1.6075 | numerical-evidence (float64, N=600) |
| c₁² ≈ 0.723 | numerical-evidence (float64, N=600) |
| Coulomb 有限秩扰动良定义 | analytic（紧算子 + 有限秩） |
| Tc/Tc₀ 表 | numerical-evidence (float64, N=400) |
