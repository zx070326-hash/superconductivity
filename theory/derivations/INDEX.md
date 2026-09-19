# theory/derivations · 推导稿索引（研究图谱）

> **这是什么**：`theory/derivations/` 下每篇推导稿的一句话地图 —— 证了什么 / 严格性等级 / 数值门 / 版本状态。
> 目的是让任何一个新的 Claude Code 会话（或人）在不通读全部文件的情况下，快速定位「哪条结论硬、硬到什么程度、门在哪、下一步在哪」。
> **本索引只记录已入库（tracked，在 `p4-conj2` 分支上）稿件的 committed 状态**；六方向深挖轮的更深版本**部分已并入 `p4-conj2` 集成分支、部分仍 held 在各自 `deep-*` 分支**，逐分支合并状态见文末 §5 与 auto-memory。（`p4-conj2` 本身仍**未合并到 `main`**。）
> 维护：改动某稿的严格性/门/结论时，同步改这里对应行。最后更新 2026-09-17。

## 图例

**严格性等级**（房规，不可越级上调）：
- **解析** = 已证明，零计算机辅助（或仅有理数精确算术）。最高等级。
- **认证 / certificate** = 外向舍入区间算术（mpmath.iv / Fraction）给出的 sound 包围。
- **数值证据 / numerical-evidence** = float64 点求值 / 有限网格，仅为证据、非证明。
- 一篇稿常是**混合**（不同命题不同等级）；下方标注的是「承重结论」的等级。

**状态**：全部 29 篇均在 `p4-conj2` 分支（**未合并到 `main`**，按用户指令 hold）。版本号取自各稿头部。

**深挖轮方向记号**：A=普适速度极限 · B=全局双侧 Tc(λ) · C=η/M「为何是氢」 · D4=μ\*>0 唯一性/极值 · D5=引理 C 凹性 · D6=各向异性。权威状态见 auto-memory `p4-six-direction-round-2026-09-15.md`。

---

## §0. 从哪读起

- **`P4_main_theorems.md`（v0.8，349 行）** —— **P4 的骨架与总账**。逐条列出每个陈述是「解析证明 / 计算机辅助证书 / 数值常数」，并记录合并沿革、被推翻的路线、记账缺陷的标注。**要理解 P4 全貌，先读这一篇。** 它本身无独立数值门（是整理稿）。
  - 主脊：KAY **Conjecture 2**（= 定理 1，全 λ 严格 Tc 上界）与 **Einstein 极值**（= 定理 2 / 定理 E）经**引理 M**（核单调性，已证无 CHALLENGED）双双降为**三行解析推论**；**𝔤(2) 解析化**（见 §4.1）消掉定理 1、2 里最后的计算机成分 ⟹ 定理 1、2 全链无计算机辅助。**[2026-09-18 标级订正，保留原句]**：定理 1、2 的陈述与证明是 `analytic`；𝔤(2) 的具体数值区间（CERT-A/B/D）是 `rational-exact`（CLAUDE.md §1 新档位，解析定理 + 有限次精确有理算术求值，须标格数与有理数个数），不应简称「解析」——「全链无计算机辅助」在「无浮点/无区间/无数值包围」的意义上仍成立。**唯一仍靠证书的是引理 C**（H 凹性）。⚠ **定理 1、2 的陈述均限于 μ\*=0**；μ\*>0 的推广**不是无条件成立的**——Einstein 极值性在弱耦合近抵消区 λ ≲ 2.5μ\*_rms 真实失效（§4.5，物理结果非证明缺口）。

---

## §1. P1 —— 同位素系数的主公式

- **`P1_isotope_master_formula.md`**（v0.7）—— 同位素系数 α 的主公式 + 四因子分解，H3S 用 Mishra Figshare 数据做数值分解，LaH10 顶角路径。严格性：**数值证据 + 一手数据锚**（解析四因子结构）。门：`p1_sympy_checks.py`、`p1_ad_sensitivity.py`（机械核对）。
- **`P1b_structural_fifth_factor.md`**（v0.1）—— **负结果**：结构相变的同位素移动是真实量子核效应（12–16 GPa），但**不污染**报道的 H3S α；把 P2 §6 的「第五因子」猜想降级。严格性：文献 + 解析论证。无脚本（`p1b_alpha_of_P.py` 明确不建）。

## §2. P2 —— 非绝热顶角修正的符号与标度

- **`P2_vertex_sign_criterion.md`**（v0.7）—— **机制已解决**：正的顶角核经 Z/Δ 双通道进入、净压低 Tc，符号与同位素比与 EPW FSR 相容。严格性：**数值证据**。门：`p2_schrodi_fsr.py`（Schrödi FSR 双通道）。
- **`P2_numeric_SPEC.md`** —— 一阶电声顶角函数 P(q,ω;N(ε)) 的数值积分**任务规格**（源公式 Grimaldi-Pietronero-Scattoni 1999 Eq.6）。配 `p2_vertex_numeric.py` + `tests/test_p2_vertex.py`（4 个解析极限 <2%）。
- **`P2_vertex_momentum_report.md`** —— 补充探索：动量结构能否给负号顶角。门：`p2_vertex_momentum.py` + 测试 40 项。

## §3. P3 —— 强耦合 Tc 上界的压力依赖

- **`P3_tc_max_pressure.md`**（v0.3.4）—— 把 Tc(P)=Tc^emp(λ,ω_log)·R(P) 分解，**圆顶下降沿完全由亏损比 R 承担、约定无关**；第二线：λ\_\* 稳定性上界的半解析压力展开 **κ=−x²g″/g=5.6035、c₃=0.8959**（闭式，从 Semenok Eq.S12 精确复现两端点）。严格性：**混合**——κ/c₃ 闭式**解析**、LaH10 (ω_log,λ,Tc)(P) 一手数据锚、「λ 随压力降 ⟹ 亏损加深」约定无关**稳健结论**。门：`p3_lambda_star_kappa.py`、`p3_lambdastar_universality.py`、`theory/solvers/lambda_star.py`（10/10）。

---

## §4. P4 —— KAY 严格 Tc 上界（主攻方向 / crown jewel）

### 4.1 主定理与骨架

- **`P4_conjecture2.md`**（v0.4）—— 把 KAY **Conjecture 2**（= 定理 1(c) 界半句）严格化约为 **Einstein 逐点不等式**，段 I 认证 ϖ∈(0,6.0]、段 II/III 解析、单调性 [√2,30] ⟹ **引理 B 在全 (0,∞) 闭合**。严格性：化约**解析** + 分段**认证**。门：区间机器（`p4_interchange_uniformity.py` 等）+ float64 探针。
- **`P4_g2_analytic.md`**（v1.1）—— **𝔤(2) 的纯解析双侧界**：1.3181093 ≤ 𝔤(2) ≤ 1.3181616（相对宽 3.96e-5）。**只用有理数精确算术，无浮点、无舍入余量**。这就是消掉定理 1、2 计算机成分的那一块。严格性：**纯解析**。门：`p4_g2_analytic_bounds.py`（389 检查 / 8s）。
- **`P4_einstein_extremal.md`**（v0.4，1276 行）—— **定理 E（Einstein 极值）**：固定 λ、⟨ω²⟩ 时 Einstein 谱使 Tc 最大。段 C-I 认证 ϖ∈(0,45]、段 C-II 未开；经引理 M 已降为解析推论，仅**引理 C（H 凹性）**仍需证书。严格性：定理**解析**（经引理 M）+ 凹性段**认证**。门：主会话复跑 463 检查 + float64 核。
- **`P4_rigorous_tc_bound.md`**（v0.1）—— P4′ 主线：把 P3 谱形展开搬到 KAY 严格算子上，有限-N 定理完成、N→∞ 换序如实分层为 OPEN（由下条解决）。门：`p4_kappa_N.py`（1257）、`p4_tc_table.py`（79）、`p4_remainder_interchange.py`（8089）。
- **`P4_interchange_uniformity.md`**（v0.6）—— 把 N→∞ 换序从「7 节点网格证书」**升为 [0.3,5] 区间定理**（inf gap ≥ 0.093225）。严格性：**认证区间定理**。门：`p4_interchange_uniformity.py`（**1668 检查 / 112.5s**，主会话独立复跑全过）。

### 4.2 μ\*=0 唯一性与核心引理

- **`P4_uniqueness_strategy.md`**（v0.3，已证已审）—— **Conjecture 2 剩余半句（μ\*=0 唯一性）**：T↦𝔨(P,T) 严格单调递减，锥不变性 + Abel 分解把问题化为两条与 P、T 无关的组合不等式。引理 A/B 已初等证明，reviewer PROVED 无 CHALLENGED。严格性：**解析**。门：`p4_uniqueness_lemmas.py`（37）。
- **`P4_concavity_q_route.md`**（102 检查）—— **引理 C 的 q-单调化路线**：(Q1)(Q2)(A) 解析证成；(B) 归约为两个明确解析缺口（ϖ∈[1.46,1.63] 谱隙下界 + 交换子 Frobenius 范数上界），真隙下 (B) 全程闭合、余量 4.4 倍。严格性：**部分解析 + 认证缺口**。门：`p4_concavity_q_route.py`（102 ck / 9.9s）。
- **`P4_gap_via_variance.md`**（54 检查）—— **Gap 1 闭合（方差路线，绕开交换子 Frobenius 范数）**：定理 F（bound1≤bound2）+ 引理 V（Bhatia–Davis）+ 定理 G（|C|<E 判据）；B NOT ⪰ 0（β_min<0），全程闭合（worst 0.29）。严格性：定理 F/引理 V **analytic**，β_max/β_min/定理 G **numerical-evidence**。门：`p4_gap_via_variance.py`（54 ck / 1.4s）。
- **`P4_concentration_bound.md`**（123 检查）—— **顶本征向量集中界**（引理 C 与 smooth-coulomb 的共享瓶颈）：结构恒等式 (★)(★★)、先验界全部解析；主廓线下界 v_n ≥ c·v₀·u_n 在 **ϖ≥5** 给显式正常数 c；小 ϖ 仅数值证据。严格性：**解析（ϖ≥5）+ 数值证据（小 ϖ）**。门：`p4_concentration_bound.py`（**123 ck**，2026-09-15 主会话在 318cfe2 复跑 0 fail）。

### 4.3 弱耦合 / 矩量塔展开

- **`P4_weak_coupling_series.md`**（v1.0）—— **定理 W（双侧包围）+ W′（三项渐近）**，弱耦合久期级数，b₁ 闭式。严格性：**解析**。门：数值门（见稿内）。
- **`P4_b1_closed_form.md`**（v1.0）—— **b₁ = 49π²/144 − 3**（闭式 + 存在性引理 + 三路交叉验证）。严格性：**解析**。门：`p4_b1_closed_form.py`（84 检查 / 187s）。
- **`P4_b3_tower_weight4.md`**（EVIDENCE-LEVEL）—— 矩量塔 weight-4 档 κ₃/b₃/b₃^Tc：数值全阶可算，**PSLQ 双精度落空 ⟹ 无严格闭式**（已修一次求积偏差 bug）。严格性：**数值证据**。门：`p4_b3_tower.py`（22 检查）。⚠ 见 auto-memory `quadrature-bias-not-caught-by-dps-agreement`。

### 4.4 推论与拓展

- **`P4_moment_threshold.md`**（v0.1）—— **定理 2′（p-矩阈值，p≥2 iff）+ 命题 2″（p<2 反例构造）**。严格性：**解析（零计算机辅助）**。
- **`P4_tc_sandwich.md`** —— Tc 双侧夹逼定理 + Einstein Tc 凸性。严格性：理论推导完整 + 数值验证。门：`p4_tc_sandwich.py`。
- **`P4_extremal_constraints.md`**（v0.1 初稿）—— 引理 M + Jensen 套路推广到其他约束（固定 ω_max/⟨ω⟩/ω_log/面积/各向异性压缩），5 命题 + 数值门。门：`p4_extremal_constraints.py`。
- **`P4_strong_coupling_expansion.md`**（v0.2）—— **定理 S（强耦合展开含谱形修正）**（整理稿，a₂/a₃ 已订正）。门：小数值门。
- **`P4_universal_tc_envelope.md`**（定理，单边）—— **【方向 A 基版】** 对所有 λ、所有谱 P 的单一标量上界 Tc < C·√(λ⟨ω²⟩)；Einstein bracket 在 λ∈[0.345,34.4] finite-N 认证、外部解析渐近。门：`p4_universal_tc_envelope.py`（**PASS 66**）。深挖版见 §5-A。
- **`P4_global_tc_certified.md`**（定理级）—— **【方向 B 基版，μ\*=0】** 把示意性全局 Tc(λ) 曲线升为**严格双侧 bracket** F_lo≤Tc/ω_E≤F_hi（Einstein 谱）。门：`p4_global_tc_certified.py`（N=2000）。深挖版（全 λ 深带、零标定常数）见 §5-B。

### 4.5 库仑 μ\*>0（把 μ\* 引入 KAY 严格算子）

- **`P4D_coulomb_operator.md`**（v0.6）—— 攻 U4：库仑赝势的严格算子形式 **𝔎_c = 2|u⟩⟨u|（秩一、半正定）**，首版 −D 项已证伪并订正；换序区间定理据 μ̃=μ\*/λ 单参族认证。严格性：**解析算子形式 + 认证轮**。门：`p4d_coulomb_falsify.py`（134）+ `verify_coulomb_fold.py`（30）。
- **`P4D_coulomb_uniqueness.md`**（v0.2）—— **μ\*>0 唯一性（Conjecture 2 库仑推广）**：秩一预解式阈值策略；**窗内定理已给解析论证**（常数 float64、有限 N 为高估），跨窗全局唯一性在 r=3 处真实失效、已量化。严格性：**窗内解析 + 数值证据边界**。门：`p4d_coulomb_uniqueness.py`（**432 检查 / 371s**，PN=400，覆盖 11 物理锚点）。
- **`P4D_coulomb_extremal_threshold.md`**（v0.2）—— **库仑阈值定理（Einstein 极值版）** Q(P)≤Q(E) 的显式条件（六常数全在 A_P 上求值，可证形式）。严格性：定理陈述 **analytic**；阈值数值 **numerical-evidence**（float64、N 有限、六常数未区间化）。**§6 定理 CE′（久期精确根）**：Davis–Kahan 整段可去，14/14 锚点认证（旧 13/14），最劣余量 1.377（旧 1.002），κ<1 前提解除。门：`p4d_coulomb_extremal_threshold.py`（G1–G4，PN=300，≤10min）。
- **`P4_coulomb_global_tc_certified.md`**（v0.2）—— **【方向 B 库仑】** 物理库仑（μ\*>0）**首个严格认证全局 Tc(λ) 双侧包围**。严格性：**认证**。门：`p4_coulomb_global_tc_certified.py`（自检 44 / 主管线 78）。
- **μ\*>0 Einstein 极值性失效边界**（笔记 `mustar_extremality_failure_map_2026-09-17.md`）—— 固定 (λ, ⟨ω²⟩) 下极值性失效分两类：**I（小 ϖ）**需 μ\*≥1.57、对 Tc 无影响；**II（弱耦合近抵消）λ ≲ (1.773+3.760μ\*)μ\*_rms** 确实失效。机制：Anderson–Morel 延迟屏蔽。严格性：领头阶判据 λ<2μ\* **analytic**；精确边界 **numerical-evidence**（N 收敛）。
- **μ\*>0 凹性秩一精确预解式**（笔记 `concavity_exact_resolvent_2026-09-17.md`）—— H″=−2Ψ+R，恒等式 R=2⟨y_⊥,(h−A)^{−1}y_⊥⟩ **analytic**；j=1 负项可保留（抵消 1–97%）。N_c=⌊rϖ⌋ 绑定后失效区 ϖ<ϖ_c(r)：r=20→≤0.6、r=10→0.8、r=3 在 ϖ≥0.7 时 μ̃_c≥0.16>物理 μ̃。严格性：恒等式 **analytic**；区域边界 **numerical-evidence**（FD 1.2e-6 验证、N=400/800/1600 四位不变）。

### 4.6 物理应用与负结果

- **`P4_vhs_dtc_dp_sign.md`**（v0.3）—— van Hove 峰到 E_F 的距离 δ 如何定 dTc/dP 符号：**符号判据 sign(∂E/∂δ)=−sign(δ) 严格成立（28 区间零违反）是唯一无条件可交付结论**；单洛伦兹拟合失效、δ 在 H3S 区间**非主旋钮** ⟹ 用 δ(P) 解释 H3S 圆顶的路线被**证伪**（干净负结果）。严格性：符号判据数值严格 + 负结果。门：`p4_vhs_enhancement.py`（25 测试）。
- **亏损预算（deficit budget）：无拟合氢化物 Tc 预测**（笔记 `deficit_budget_h3s_2026-09-17.md`）—— 定义填充因子 𝔣:=Tc/(C_∞√(λ⟨ω²⟩))，C_∞=0.1827262477 [analytic]。用真实 α²F（EPW anharmonic）+ μ\*=0.13 给出 **H3S-200: 199.0 K / D3S-200: 155.4 K**（实验 203/150 K，残差 −2.0%/+3.6%）。三项亏损分解（乘性）：δ_FIN（有限 λ）~20%、δ_BRD（谱展宽，定理 E 的不等式变成数）~13%、**δ_COU（库仑）27–30%，最大项**。反解 μ\*∈[0.12, 0.15]。严格性：C_∞ **analytic**；全部 Tc **numerical-evidence**（float64、N 有限、α²F 来自外部）。
- **负结果 NO-GO**（笔记多篇 2026-09-16）—— (a) **Haldane 量子自旋谱隙**：LTI-SDP 硬性要求 frustration-free，Heisenberg 链不满足。(b) **BCS 冷原子普适比**：接触势无算子谱结构，FHS 已占领，无声子无 P4 技术栈可迁移性。(c) **γ>2 核不在定理 E 类内**：𝒫 类尾部恰为 k⁻² 或更慢（二阶矩有限/无穷），γ>2 要求更快衰减 ⟹ γ>2 ∩ 𝒫 = ∅ [analytic]。

### 4.7 γ-model 量子临界配对 + Coulomb（新方向 2026-09-16）

- **方向选择与可行性**（笔记 `gamma_model_direction_selection_2026-09-16.md`）—— KAY/Elezaby–Abanov 明确放弃 μ\*；P4-D Coulomb 机器直接可搬。Phase 0 通过：𝔤(1)≈1.6075（KAY 浮点下界），c₁²=0.723（float64）。
- **定理 I（Δ 通道免疫）**（笔记 `delta_channel_immunity_2026-09-16.md`）—— **配对通道排斥无法摧毁秩一主导性**：M ⪯ ℌ₊ + ρ(f₋)·e₀e₀ᵀ，λ₂(M) ≤ (1−c₁²)μ₁ + c₁²μ₂，对排斥强度/谱形/秩一致。严格性：引理 1–4、定理 I、推论全 **analytic**；c₁²(ϖ) 数值为 **numerical-evidence**。
- **c₁² certificate 下界**（笔记 `gamma_phase2_multi_gamma_scan_2026-09-16.md` + `gamma_tail_bound_audit_v2_2026-09-16.md`）—— γ=1, N=400, μ\*=0: **c₁²_lo = 0.704** [certificate v3, 块对角 Davis–Kahan]（float64 真值 0.723，gap 2.6%）；γ=0.67: 0.543（从空真 0.028 恢复）。Phase 2 全扫描 5γ × 4μ\*：𝔤_lo 全行 [certificate]；γ≤0.5 非紧（尾部发散）。生产脚本 `gamma_dk_coulomb_schur_v2.py` + `gamma_tail_bound_v2.py`（未提交，Phase 1 v1 有三处单边缺陷、数字事后成立）。
- **Anderson–Morel 屏蔽**（笔记 `am_plateau_2026-09-17.md`；另见 `deficit_budget_h3s_2026-09-17.md` §5(b)）—— Δ𝔤 对 N_c **非单调**，μ_eff=μ/(1+μ ln N_c/(2𝔤))；H3S-200 anharmonic 在 r=10/20/40/80/160 上 Tc = 199.0/198.4/198.3/198.4/198.4 K（0.1% 变化跨 8 倍截断），落在平台。严格性：恒等式 **analytic**；全部数值 **numerical-evidence**（N=400）。独立 claim 候选（但 (1+0.62λ) 归属为 McMillan 1968 非 AD，且与 f₁ 简并分不开）。

---

## §5. 六方向深挖轮（部分已并入 `p4-conj2`，部分仍 held）

> 这些是在 `p4-conj2` 之上、各自 `deep-*` 分支上完成、**均已个人核验闭环**的更深结果。截至 2026-09-15，另一会话已把其中若干**并入 `p4-conj2` 集成分支**（见下表「合并」列），其余仍 held 在分支上；**`p4-conj2` 整体仍未合并到 `main`**。**权威状态、门数、订正历史见 auto-memory `p4-six-direction-round-2026-09-15.md`**（本表只给指针，避免与易变的分支状态重复而失同步）。合并状态用 `git merge-base --is-ancestor <branch> p4-conj2` 核验。
>
> 另有若干**基版 / 共享瓶颈分支**已并入 `p4-conj2`：`deep-universal-tc-envelope`（A 基版，§4.4）、`deep-coulomb-global-tc`（B 库仑版，§4.5）、`deep-concentration-bound`（D5/smooth-coulomb 共享瓶颈，§4.2）、`deep-b3-tower`（§4.3 矩量塔）。
>
> **合并后主会话在 318cfe2 上重跑门确认无回归（2026-09-15，`OMP_NUM_THREADS=1`）**：`p4_g2_analytic_bounds.py` 389/0 · `p4_universal_tc_envelope.py` 66/0 · `p4_global_tc_certified.py` 2018/0 · `p4_concentration_bound.py` 123/0 · `p4_b3_tower.py` 22/0 · `p4_coulomb_global_tc_certified.py` 78/0。合并为纯增量、无冲突、交付文件与源分支逐字一致。详见 auto-memory `p4-conj2-consolidation-reverify-2026-09-15.md`。

| 方向 | 分支 @tip | 合并 | 头条结论 | 严格性 | 门 |
|---|---|---|---|---|---|
| **A** 普适速度极限 | 已并入主线 @318cfe2 | **已并入** | Tc < C_∞√(λ⟨ω²⟩)，C_∞=(1/2π)√𝔤(2)=0.18273，**所有 P、所有 μ\*≥0**（天花板**不依赖**包络门） | **解析**（单边） | 承重 `p4_g2_analytic_bounds.py` **389/0** + `p4_universal_tc_envelope.py` **66/0**〔memory 旧记的 `..._coulomb.py` 159/159 为已删 worktree 里的**未提交探针**、非承重、不可复现〕 |
| **B** 全局双侧 Tc(λ) | `deep-global-tc-allcoupling` @6a5572a | 未合并（+11） | 所有 λ>0 双侧 bracket；**深带经 Prop.C 闭式 (C1) 变为完全解析、零标定常数** | 深带**解析** / 其余认证 | `p4_global_tc_allcoupling.py` **104/104**（全，150.2s）已 team-lead 权威复跑 |
| **C** η/M「为何是氢」 | 已并入主线（`etaM-corrected` @ed6de52，经 9e1c83d 合并） | **已并入** | 定理 A：Tc < 0.18273·√(η/M)，η=N(0)⟨I²⟩ Hopfield，**μ\* 无关** | 定理 A **解析** | `P4_eta_over_M.py` 261/261（116.3s） |
| **D4** μ\*>0 唯一性/极值 | `deep-mustar-completion` @36ba77a | 未合并（+1） | 多原子唯一性 + Einstein 极值机制逐字推广；L\*\*（数值探测）；L\*（κ≥1 近交叉）诚实 **OPEN** | **数值证据**（自标非定理） | `probe_Lstarstar_uniqueness_multiatom.py` 7/7 |
| **D5** 引理 C 凹性 | `deep-lemmaC-analytic` @2ff4d54 | 未合并（+2） | E(x)>0、C≤0、q′=E−|C| 解析；残余 **|C|<E 仅数值证据**（诚实边界，廓线上界为已识别缺失件，四条路线已排除） | 结构**解析** + 边界**数值证据** | 见 `P4_concentration_bound.py` |
| **D6** 各向异性增强 | `deep-anisotropy-enhancement` @5e5653d | 未合并（+1） | 定理 A：finite-N Bhatia–Davis λ_max 下界；定理 B：Tc 转换 | finite-N **解析** / N→∞ 数值 | 门 1280/1280（2.8s） |

---

## 附：配套目录

- **数值门脚本**：与各 `.md` 同名的 `.py`（如 `p4_g2_analytic_bounds.py`）；运行前 `export OMP_NUM_THREADS=1`。
- **求解器**：`theory/solvers/`（`eliashberg_iso.py` 全非线性求解器、`lambda_star.py` λ\_\* 泛函）。见该目录 README。
- **过程脚本 / 探针**：`theory/notes/_scratch_ir/`（git-ignored，非承重、含 rebuild 源）。
- **审查记录**：`theory/notes/*/*_review_*.md`（对抗核验审计轨迹；2026-09-18 起 notes 已按主题分子目录，路由见 `theory/notes/INDEX.md`）。
- **笔记语料**：`theory/notes/`（69 篇 2026-09-16/17 新增）—— 证据轨迹、探索过程、订正记录。**权威状态在 auto-memory**（`MEMORY.md` + `memory/*.md`），笔记是原始材料不是结论。
- **文献 / 一手数据**：`library/`、`theory/data/`（各有 README/MANIFEST）。
- **总纲**：`theory/PROGRAM.md`（研究纲领叙事）；`README.md`（工作区总览）。
