# theory/derivations · 推导稿索引（研究图谱）

> **这是什么**：`theory/derivations/` 下每篇推导稿的一句话地图 —— 证了什么 / 严格性等级 / 数值门 / 版本状态。
> 目的是让任何一个新的 Claude Code 会话（或人）在不通读全部文件的情况下，快速定位「哪条结论硬、硬到什么程度、门在哪、下一步在哪」。
> **本索引只记录已入库（tracked，在 `p4-conj2` 分支上）稿件的 committed 状态**；六方向深挖轮的更深版本**部分已并入 `p4-conj2` 集成分支、部分仍 held 在各自 `deep-*` 分支**，逐分支合并状态见文末 §5 与 auto-memory。（`p4-conj2` 本身仍**未合并到 `master`**——仓库当前无 `main` 分支，只有 `master`；`p4-conj2` 领先 `master` 数十个提交，2026-09-19 复核确认。）
> 维护：改动某稿的严格性/门/结论时，同步改这里对应行。最后更新 2026-09-29（接线合并）。

## 图例

> 严格性等级定义与维护规则已搬至本目录的 [`CLAUDE.md`](CLAUDE.md)（单一真源，自动加载）。下方标注的是各稿「承重结论」的等级。

**状态**：全部 **41 篇**（= 目录内 `.md` 43 个减去 `INDEX.md`/`CLAUDE.md`）均在 `p4-conj2` 分支（**未合并到 `master`**，按用户指令 hold）。版本号取自各稿头部。　*（2026-09-29 接线：+4 篇 = `P4_cone_coercivity` / `P4_lemmaC_var_bound` / `P4_mustar_uniqueness_strategy` / `P4_strat_pi2_tower`，已合入 master。）*

> **2026-09-20 记账订正**：此处原写「29 篇」，但**逐文件核对发现 7 篇从未入索引**
> （`P4_commutator_trace_bound` / `P4_profile_upper_bound` / `P4_eta_over_M` 三篇是**旧欠账**，
> 另四篇是本日并入的 `deep-*` 交付）。七篇已全部补行，计数改为**逐文件实测值 37**，
> 不再手工维护。自查：`cd theory/derivations && for f in *.md; do grep -q "\`$f\`" INDEX.md || echo "MISSING: $f"; done`。

**深挖轮方向记号**：A=普适速度极限 · B=全局双侧 Tc(λ) · C=η/M「为何是氢」 · D4=μ\*>0 唯一性/极值 · D5=引理 C 凹性 · D6=各向异性。权威状态见 auto-memory `p4-six-direction-round-2026-09-15.md`。

---

## §0. 从哪读起

- **`P4_main_theorems.md`（v0.8，**352 行**；2026-09-25 复核行数，原写 349 已漂移）** —— **P4 的骨架与总账**。逐条列出每个陈述是「解析证明 / 计算机辅助证书 / 数值常数」，并记录合并沿革、被推翻的路线、记账缺陷的标注。**要理解 P4 全貌，先读这一篇。** 它本身无独立数值门（是整理稿）。
  - 主脊：KAY **Conjecture 2**（= 定理 1，全 λ 严格 Tc 上界）与 **Einstein 极值**（= 定理 2 / 定理 E）经**引理 M**（核单调性，已证无 CHALLENGED）双双降为**三行解析推论**；**𝔤(2) 解析化**（见 §4.1）消掉定理 1、2 里最后的计算机成分 ⟹ 定理 1、2 全链无计算机辅助。**[2026-09-18 标级订正，保留原句]**：定理 1、2 的陈述与证明是 `analytic`；𝔤(2) 的具体数值区间（CERT-A/B/D）是 `rational-exact`（CLAUDE.md §1 新档位，解析定理 + 有限次精确有理算术求值，须标格数与有理数个数），不应简称「解析」——「全链无计算机辅助」在「无浮点/无区间/无数值包围」的意义上仍成立。**唯一仍靠证书的是引理 C**（H 凹性）。⚠ **定理 1、2 的陈述均限于 μ\*=0**；μ\*>0 的推广**不是无条件成立的**——Einstein 极值性在弱耦合近抵消区 λ ≲ 2.5μ\*_rms 真实失效（§4.5，物理结果非证明缺口）。

---

## §1. P1 —— 同位素系数的主公式

- **`P1_isotope_master_formula.md`**（v0.7）—— 同位素系数 α 的主公式 + 四因子分解，H3S 用 Mishra Figshare 数据做数值分解，LaH10 顶角路径。严格性：**数值证据 + 一手数据锚**（解析四因子结构）。门：`p1_sympy_checks.py`、`p1_ad_sensitivity.py`（机械核对）。⚠ **2026-09-25：`p1_ad_sensitivity.py:23` 的 `f₂` 缺陷 `(r²−1)→(r−1)` 已修，§3.1/§3.2a 两张 AD 表与 `s_λ≥½` 的 3.7 万格点扫描全部重解**（符号规则结论存活；旧扫描证据作废）——见 `../notes/ceiling/f2_defect_repair_2026-09-25.md`。
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
- **`P4_einstein_extremal.md`**（v0.4，**1291 行**；2026-09-25 复核行数，原写 1276 已漂移）—— **定理 E（Einstein 极值）**：固定 λ、⟨ω²⟩ 时 Einstein 谱使 Tc 最大。段 C-I 认证 ϖ∈(0,45]、段 C-II 未开；经引理 M 已降为解析推论，仅**引理 C（H 凹性）**仍需证书。严格性：定理**解析**（经引理 M）+ 凹性段**认证**。门：主会话复跑 463 检查 + float64 核。**§2 另含两条便宜引理，此前未在本图谱登记 ⟹ 2026-09-25 补登**：**引理 X**（存在性：∀P∈𝒫、λ>0，∃T>0 使 λ𝔨(P,T)≥1，**analytic**）与 **★推论 μ\***（第 283 行：`K_c ⪰ 0 ⟹ 𝔎−μ̃K_c ⪯ 𝔎` ⟹ **定理 E 与 Conjecture 2 的界半句对任何 μ\*≥0 逐字成立**，**analytic**，P4-D 算子序 Layer 1）。⚠ **作用域上限（第 285–288 行的边界说明，不许顺手扩大）**：它断言 `Q_μ(P) ≤ Q_0(E)`（用 μ\*=0 的 Einstein 值当上界），**不是** μ\*>0 下的极值性 `Q_μ(P) ≤ Q_μ(E)` —— 后者**已知为假**（弱耦合近抵消区反例，边界 `λ_c/μ*_rms = 1.773 + 3.760 μ*_rms`，numerical-evidence）。**消费点**：天花板 Tier 0 判据据此**不需要任何 μ\* 输入**（`../notes/ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md` §4/§5/§6.6 订正 D）。**补登的理由本身是一次事故**：2026-09-25 我把这条**早已并入**的 analytic 推论在另一页重新登记成了「未关闭的声音性缺口」，共因就是它不在本图谱上、grep 不到。
- **`P4_rigorous_tc_bound.md`**（v0.1）—— P4′ 主线：把 P3 谱形展开搬到 KAY 严格算子上，有限-N 定理完成、N→∞ 换序如实分层为 OPEN（由下条解决）。门：`p4_kappa_N.py`（1257）、`p4_tc_table.py`（79）、`p4_remainder_interchange.py`（8089）。
- **`P4_interchange_uniformity.md`**（v0.6）—— 把 N→∞ 换序从「7 节点网格证书」**升为 [0.3,5] 区间定理**（inf gap ≥ 0.093225）。严格性：**认证区间定理**。门：`p4_interchange_uniformity.py`（**1668 检查 / 112.5s**，主会话独立复跑全过）。

### 4.2 μ\*=0 唯一性与核心引理

- **`P4_uniqueness_strategy.md`**（v0.3，已证已审）—— **Conjecture 2 剩余半句（μ\*=0 唯一性）**：T↦𝔨(P,T) 严格单调递减，锥不变性 + Abel 分解把问题化为两条与 P、T 无关的组合不等式。引理 A/B 已初等证明，reviewer PROVED 无 CHALLENGED。严格性：**解析**。门：`p4_uniqueness_lemmas.py`（37）。
- **`P4_concavity_q_route.md`**（102 检查）—— **引理 C 的 q-单调化路线**：(Q1)(Q2)(A) 解析证成；(B) 归约为两个明确解析缺口（ϖ∈[1.46,1.63] 谱隙下界 + 交换子 Frobenius 范数上界），真隙下 (B) 全程闭合、余量 4.4 倍。严格性：**部分解析 + 认证缺口**。门：`p4_concavity_q_route.py`（102 ck / 9.9s）。 **2026-09-29 合入 Phase 2**（`worktree-agent-af21da55…`，含 bf7b1bc「Fix overclaim」）：`p4_concavity_gap2_gates.py` **165 ck / 0 fail（主会话 2026-09-29 复跑）**；谱隙下界 γ_lo>0 **仅 ϖ>1.3188 纯解析**（引理 7.2a，e₀ 单试探 + Θ\* + φ 单调性），ϖ∈[1.0,1.3188] 证书级、ϖ<1.0 含证书级输入；**并行分支 a8932ec 的「γ_lo>0 对一切 ϖ 成立」正是被订正的过度声称——冲突取订正侧，原版留在合并的第二父提交**。
- **`P4_cone_coercivity.md`**（分支 `g2-cone-coercivity`，2026-09-29 合入；分支自注「勿合并」已由用户放行覆盖）—— **锥-矫顽性引理**：把 `P4_concavity_q_route.md` §3 定理 A 的定性 `E(x)>0` 升级为定量 `E(x) ≥ κ(x) > 0`（三条由弱到强的界 + 残余障碍图，§7.4 缺口 #4）。严格性：**稿内自述 rigorous（analytic），本会话未逐条复核**。门：`p4_cone_coercivity_probe.py`——稿内写 191 ck，**主会话 2026-09-29 复跑为 196 ck / 0 fail**（计数漂移，未影响通过）。被 `P4_concentration_bound.md` 单向引用。
- **`P4_lemmaC_var_bound.md`**（分支 `worktree-agent-ae6c83eb…`，2026-09-29 合入）—— **引理 C 方差腿**：把 `Var_v(B)<γ·xE/2`（`P4_concavity_q_route.md` §4 终端不等式）归约为单个非负核算子 `M:=𝒪[f−f²]` 的标量谱不等式；三条精确恒等式用 mpmath dps=30 验证。严格性：**归约 + 精确恒等式为稿内自述证明/证书级；残余为已知瓶颈、数值证据**。门：`p4_lemmaC_var_gate.py` **217 ck / 0 fail（主会话 2026-09-29 复跑）**。
- **`P4_gap_via_variance.md`**（54 检查）—— **Gap 1 闭合（方差路线，绕开交换子 Frobenius 范数）**：定理 F（bound1≤bound2）+ 引理 V（Bhatia–Davis）+ 定理 G（|C|<E 判据）；B NOT ⪰ 0（β_min<0），全程闭合（worst 0.29）。严格性：定理 F/引理 V **analytic**，β_max/β_min/定理 G **numerical-evidence**。门：`p4_gap_via_variance.py`（54 ck / 1.4s）。
- **`P4_concentration_bound.md`**（123 检查）—— **顶本征向量集中界**（引理 C 与 smooth-coulomb 的共享瓶颈）：结构恒等式 (★)(★★)、先验界全部解析；主廓线下界 v_n ≥ c·v₀·u_n 在 **ϖ≥5** 给显式正常数 c；小 ϖ 仅数值证据。严格性：**解析（ϖ≥5）+ 数值证据（小 ϖ）**。门：`p4_concentration_bound.py`（**123 ck**，2026-09-15 主会话在 318cfe2 复跑 0 fail）。
- **`P4_lemmaC_final_status.md`**（235 行，2026-09-20 由 `deep-lemmaC-analytic` 并入）—— **引理 C 的终局盘点：答案是「证书级，不是完全解析」**。链条 `引理 C ⟸ H″<0 ⟺ q′>0 ⟸ E>0 且 |C|<E` 中 `E>0` 已解析；`|C|<E` 归约到单一标量不等式 `(★) Var_v(B) < γ(q−m₃)`，证书级恒真（余量 ≥4×），但换成已知解析下界（cone-coercivity）只捕获真值 **0.58–0.73**，需 ≳0.90。本稿新做**三条**闭合尝试**全部失败**。严格性：结构 **analytic** + 残差边界 **numerical-evidence**。门：`p4_lemmaC_residual_gate.py`（**8 ck**，2026-09-20 主会话合并后从零复跑 0 fail）。
- **`P4_commutator_trace_bound.md`**（255 行）—— 交换子 Frobenius 范数上界的**死路记录** + Gap 1 闭合：(i) `‖[A,B]‖_F ≤2×` 的六条直接路线全部 ≥3.2× 松，判定**已排除**；(ii) 命题 4（`bound1≤bound2`，解析恒等式）使 (i) 对引理 C 的闭合**不是必需的**。严格性：命题 4 **analytic**；系数 0.4462 那条 **numerical-evidence**（自标脆弱）。
- **`P4_profile_upper_bound.md`**（271 行）—— 秩一预解式廓线上界：**结构性死路**（诚实负结果）。构造本身比 Cauchy–Schwarz 路线质变更优（比值随 ϖ、n **改善** 1.79×→1.12×，与 Attempt #4 的发散方向相反），但接入 L1+L2 朴素恢复判据后 `max_n c_needed` 超过任何合法浓度常数的理论上限 `c_best`，**即便代入真值仍在 ϖ≥8 不闭合**。对捕获率目标的贡献：**0**。严格性：**numerical-evidence**（负结果）。

### 4.3 弱耦合 / 矩量塔展开

- **`P4_weak_coupling_series.md`**（**v1.1**，2026-09-29 合入 W″ 折入版）—— **定理 W（双侧包围）+ W′（三项渐近）+ W″（矩量塔：Λ 变量下 γ-free 全阶展开 + b₂ 首算 + π²-塔猜想收回）**，弱耦合久期级数，b₁ 闭式。严格性：W/W′ **解析**；W″/b₂ 部分见稿内分级（b₂ 判决为**证据级**，条件于开问题 γ∉ℚ[π²]）。门：`p4_moment_tower.py` **21 检查全过（主会话 2026-09-29 复跑，87 s）**。
- **`P4_strat_pi2_tower.md`**（strategy + lemma-skeleton，分支 `worktree-agent-a0987668…`，2026-09-29 合入）—— **π²-塔猜想对 h(ϖ)**：约化、b₂、对抗性判决；**原猜想在 n=2 处被强证据否定（非严格证伪，条件于开问题 γ∉ℚ[π²]）**。严格性：**证据级**。探针 `pi2_tower_*.py`（已入库，本会话未复跑）。
- **`P4_b1_closed_form.md`**（v1.0）—— **b₁ = 49π²/144 − 3**（闭式 + 存在性引理 + 三路交叉验证）。严格性：**解析**。门：`p4_b1_closed_form.py`（84 检查 / 187s）。
- **`P4_b3_tower_weight4.md`**（EVIDENCE-LEVEL）—— 矩量塔 weight-4 档 κ₃/b₃/b₃^Tc：数值全阶可算，**PSLQ 双精度落空 ⟹ 无严格闭式**（已修一次求积偏差 bug）。严格性：**数值证据**。门：`p4_b3_tower.py`（22 检查）。⚠ 见 auto-memory `quadrature-bias-not-caught-by-dps-agreement`。

### 4.4 推论与拓展

- **`P4_moment_threshold.md`**（v0.1）—— **定理 2′（p-矩阈值，p≥2 iff）+ 命题 2″（p<2 反例构造）**。严格性：**解析（零计算机辅助）**。
- **`P4_tc_sandwich.md`** —— Tc 双侧夹逼定理 + Einstein Tc 凸性。严格性：理论推导完整 + 数值验证。门：`p4_tc_sandwich.py`。
- **`P4_extremal_constraints.md`**（v0.1 初稿）—— 引理 M + Jensen 套路推广到其他约束（固定 ω_max/⟨ω⟩/ω_log/面积/各向异性压缩），5 命题 + 数值门。门：`p4_extremal_constraints.py`。
- **`P4_strong_coupling_expansion.md`**（v0.3）—— **定理 S（强耦合展开含谱形修正）**（整理稿，a₂/a₃ 已订正）。门：小数值门。⚠ **2026-09-25：§5 推论 S.3 的 AD `f₂` 母式整条抄错已订正**（旧式连 `f₂≥1` 都不满足、两个极限全错），小 σ 展开随之重推为 **analytic** 闭式 `β·λ²/(λ²+Λ₂(1)²)`；σ³ 那格由「无」降为「未定」——见 `../notes/ceiling/f2_defect_repair_2026-09-25.md`。
- **`P4_universal_tc_envelope.md`**（定理，单边）—— **【方向 A 基版】** 对所有 λ、所有谱 P 的单一标量上界 Tc < C·√(λ⟨ω²⟩)；Einstein bracket 在 λ∈[0.345,34.4] finite-N 认证、外部解析渐近。门：`p4_universal_tc_envelope.py`（**PASS 66**）。深挖版见 §5-A。
- **`P4_global_tc_certified.md`**（定理级）—— **【方向 B 基版，μ\*=0】** 把示意性全局 Tc(λ) 曲线升为**严格双侧 bracket** F_lo≤Tc/ω_E≤F_hi（Einstein 谱）。门：`p4_global_tc_certified.py`（N=2000）。深挖版（全 λ 深带、零标定常数）见 §5-B。
- **`P4_global_tc_allcoupling.md`**（524 行，2026-09-20 由 `deep-global-tc-allcoupling` 并入）—— **【方向 B 深挖版】弱耦合上侧 bracket 已闭合**：基版在 λ≲0.582 的上界是带经验 pad（REM=0.75）的**渐近式**、不是证书；本稿用命题 B（经命题 C 升级后**行和有精确闭式、零标定常数、完全证明**）+ 命题 U（finite-N **certificate**）取代它 ⟹ **每个 λ>0 都有两侧皆非渐近的 bracket**。逐带严格性分档见稿内表（ϖ≤0.14 一档为 **fully analytic**，Weyl 𝔤(2) 包络）。门：`p4_global_tc_allcoupling.py`（**104/104**，2026-09-20 主会话合并后从零复跑 146.1s，0 fail）。
- **`P4_eta_over_M.md`**（v0.2，505 行）—— **定理 A：Tc ≤ 0.1827262478·√(η/M)**，η=N(0)⟨I²⟩ 为 Hopfield 纯电子量 ⟹ **μ\*-无关、声子频率无关、谱形无关**，且对该界同位素指数恰为 ½（「为什么是氢」的一行回答）。⚠ **优先权**：公式形式与常数 **0.1827 都是 Leavens–Allen–Dynes SSC 30, 595 (1979) Eq.(10) 的**；本稿增量仅三点——**不等号方向（上界，LAD 给的是下界）/ 全 λ 全谱 / 常数被认证包围**。严格性：定理 A **analytic**（定理 B 不是新结果，严格性以 `P4_universal_tc_envelope.md` 为准）。门：`p4_eta_over_M.py`（261/261，116.3s）。
- **`P4_anisotropy_quantitative.md`**（v0.1，217 行，2026-09-20 由 `deep-anisotropy-enhancement` 并入）—— **命题 5 的定量化**：各向异性 Tc 增强的**显式下界** λ₁(𝒜) ≥ λ̄h + h²·Var_w(z)/(λ̄h−μ)，分子是**纯纤维量**（与 Matsubara 结构无关）、各向同性时恒为零；Z 通道经恒等式 `Mu=0` **精确**消失（非近似）。严格性：定理 A/B + 引理 1 **全解析**；N→∞ 数值。门：`p4_anisotropy_quantitative.py`（**1280/1280**，2026-09-20 主会话合并后从零复跑 3.1s，0 fail）。

### 4.5 库仑 μ\*>0（把 μ\* 引入 KAY 严格算子）

- **`P4D_coulomb_operator.md`**（v0.6）—— 攻 U4：库仑赝势的严格算子形式 **𝔎_c = 2|u⟩⟨u|（秩一、半正定）**，首版 −D 项已证伪并订正；换序区间定理据 μ̃=μ\*/λ 单参族认证。严格性：**解析算子形式 + 认证轮**。门：`p4d_coulomb_falsify.py`（134）+ `verify_coulomb_fold.py`（30）。
- **`P4D_coulomb_uniqueness.md`**（v0.2）—— **μ\*>0 唯一性（Conjecture 2 库仑推广）**：秩一预解式阈值策略；**窗内定理已给解析论证**（常数 float64、有限 N 为高估），跨窗全局唯一性在 r=3 处真实失效、已量化。严格性：**窗内解析 + 数值证据边界**。门：`p4d_coulomb_uniqueness.py`（**432 检查 / 371s**，PN=400，覆盖 11 物理锚点）。
- **`P4_mustar_uniqueness_strategy.md`**（v0.1，分支 `worktree-agent-aa561c3a…`，2026-09-29 合入）—— **μ\*>0 唯一性/极值：Perron–Frobenius 与不可约性能否破除简并**（策略稿，**自述非定理**）：归并两种猜想的精确陈述，定位简并/双根发生位置与「锥/单指标 k=1」机制失效原因。严格性：**策略 + 诚实评估，无承重结论**。探针 `p4_mustar_pf_probe.py`（已按毕业规则入库）。
- **`P4D_coulomb_extremal_threshold.md`**（v0.2）—— **库仑阈值定理（Einstein 极值版）** Q(P)≤Q(E) 的显式条件（六常数全在 A_P 上求值，可证形式）。严格性：定理陈述 **analytic**；阈值数值 **numerical-evidence**（float64、N 有限、六常数未区间化）。**§6 定理 CE′（久期精确根）**：Davis–Kahan 整段可去，14/14 锚点认证（旧 13/14），最劣余量 1.377（旧 1.002），κ<1 前提解除。门：`p4d_coulomb_extremal_threshold.py`（G1–G4，PN=300，≤10min）。 **2026-09-29**：§7 已追加局限性补注——「14/14 认证」只覆盖每锚点**一个**近简并两原子形状，放开到两原子全族后 Rayleigh 步骤在物理锚点上为假（计数「12/14」不自洽，见 `notes/coulomb/mustar_strong_coupling_extremality_2026-09-19.md` 核查注）；极值性结论本身未被推翻。
- **`P4_coulomb_global_tc_certified.md`**（v0.2）—— **【方向 B 库仑】** 物理库仑（μ\*>0）**首个严格认证全局 Tc(λ) 双侧包围**。严格性：**认证**。门：`p4_coulomb_global_tc_certified.py`（自检 44 / 主管线 78）。
- **`P4D_coulomb_multiatom_status.md`**（v0.1，177 行，2026-09-20 由 `deep-mustar-completion` 并入）—— **开放引理 L\*\* 的多原子数值探测（自标：探测报告，非定理）**：在姊妹探针找到的危险谱 `P=0.8δ(0.3s)+0.2δ(4.0s)` 上重算六常数，(a) 唯一性机制**逐字推广**、阈值 μ̃\*^{CE} 正确排除全部 7 个窗内非单调系列（**7/7**）；(b) 多原子极值刚性 Q(P)≤Q(E) 在 μ̃<μ̃\*^{CE} 时 233/233 成立、μ̃≥阈值时 24/24 正确排除。L\*（κ≥1 近交叉）诚实 **OPEN**。严格性：**numerical-evidence**。门：`probe_Lstarstar_uniqueness_multiatom.py`（**7/7**，2026-09-20 主会话合并后从零复跑 45.2s）+ `probe_Lstarstar_extremal_multiatom.py` + `probe_Lstarstar_extremal_search.py`（三个探针按毕业规则已 `git add -f` 入 `theory/notes/_scratch_ir/`；原稿写的 worktree 绝对路径已订正为仓库相对路径）。
- **μ\*>0 Einstein 极值性失效边界**（笔记 `mustar_extremality_failure_map_2026-09-17.md`）—— 固定 (λ, ⟨ω²⟩) 下极值性失效分两类：**I（小 ϖ）**需 μ\*≥1.57、对 Tc 无影响；**II（弱耦合近抵消）λ ≲ (1.773+3.760μ\*)μ\*_rms** 确实失效。机制：Anderson–Morel 延迟屏蔽。严格性：领头阶判据 λ<2μ\* **analytic**；精确边界 **numerical-evidence**（N 收敛）。 **2026-09-29 合入 `worktree-spinfluc-quant`**：新增 `notes/coulomb/mustar_strong_coupling_extremality_2026-09-19.md`（定理 TA 两原子约化 / 定理 SC 强耦合端 / (AB) 记账恒等式，analytic 部分 agent 逐行读过无漏洞）与 `mustar_stageB_budget_probe`（Rayleigh 两段式路线的预算探针）；据此，物理锚点上 CE/CE′ 的 Rayleigh 步骤为假是**存在性数值反例（float64、N 收敛）**，**不是 analytic 证明**——「已证的负结果」的说法**不采用**（核查判定，降为 numerical-evidence）；「12/14」计数不自洽（表内 11/14、Stage B 细网格 14/14），勿引。
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

### 4.8 多带/各向异性 Tc 天花板（新方向 2026-09-18）

- **规范来源（唯一现在真值入口）：`theory/notes/ceiling/multiband_ceiling_CANONICAL_2026-09-18.md`**——四篇更早的长笔记
  （`multiband_ceiling_route_Gprime_2026-09-18.md`、`multiband_ceiling_lift_2026-09-18.md`、
  `lemma_G_offdiag_GF_factor_2026-09-18.md`、`zdom_literature_census_2026-09-18.md`）已**自行退为证明存档**，
  不再单独引用，仅作证明细节备查。**定理（无条件版）**：设 `α²F^Z=α²F^ph+α²F^sf`、`α²F^Δ=α²F^ph−α²F^sf`
  两片各自逐元非负，`c_GF:=(π²/6)/𝔤(2)=1.2479204046`，则 `T_c ≤ C_∞√(ρ(M2″_Z))`，
  `M2″_Z := diag(M2^ph+c_GF·M2^sf) + c_GF·offdiag(M2^ph+M2^sf)`。严格性（文件自标）：链上每一步 **analytic**
  （(ML1)/(ML3-abs)/分块 Cauchy–Schwarz/引理 M/定理 GF-mix/Horn–Johnson 8.1.18 逐元谱半径单调）；
  唯二非 analytic 成分——① `𝔤(2)`/`C_∞` 数值 **rational-exact**（沿用 §4.1 的 CERT 档位）；
  ② 物理前提 **(Z-dom)**（非数学缺口，是对文献计算约定的核对，见下）。
  ⚠ **`𝒫` 对减法不封闭**是该定理形状的根因：`M2^ph−M2^sf`（或其绝对值）**不界任何东西**
  （`Λ^ph−Λ^sf` 可有 `k⁻⁴` 尾，掉出 `𝒫`；定向扫 `r:=M2^sf/M2^ph→1` 时 45/45 违反）——
  **早两版「(G′-abs) 覆盖 s± 全族」的主张（v1、v2）均已被证伪并撤回**，此定理是取代它们的唯一存活结论。
- **(Z-dom) 覆盖普查**（同一规范文件 §3）：两条 **analytic** 充分条件 (K1)（`Λ^Z=Λ^ph+Λ^sf` 型加减约定，自动满足）
  / (K2)（`Λ^Δ_ij=g·F_ij·λ_ij`，`g|F_ij|≤1` 充分）——普查文献 A(自动覆盖)/B(未覆盖)/C(需读数值表) = 5/0/1（83%），
  材料计 11/1（92%）；**分类单位是论文的约定不是材料**（同材料族换约定可从 A 掉到 C）。s± 族 (K2) 是"刀口"
  （`|F|≡1` 使条件退化成 `g≤1`、零余量）。μ\*>0 延拓（笔记 `ceiling/mustar_psd_multiband_ceiling_2026-09-18.md`，
  主会话从零重导）：引理 μ-1（Gram 形取代 Kronecker 形，放开逐带权重廓线）/ μ-2 / μ-3 三条全部 **analytic**，
  但判据须写成约定无关形式 `μ11μ22≥μ12μ21`（不能用 `μ12²`，平方形式在 `μ_ij=N_jK_ij` 约定下 3.5% 判错），
  且 Bochner 保证第一性原理正核构造的 `μ` 恒 PSD ⟹ 该判据在现有文献上**无信息量**（材料落地见 §4.12）。
- **非对角 GF 因子的材料分辨版**（同一规范文件 §4）：**analytic** 一步把 `sup_ϖ g(ϖ)`（无条件常数 1.248）
  换成实际谱的加权平均 `ḡ`；铁基 s± 族（Ummarino 系列，Ω₀=4.65k_BT_c 经验律）算出 `ḡ/𝔤(2)<1` 恒成立
  ⟹ **该族上非对角因子恰为 1、缺口买回 0%**（`numerical-evidence`，两套独立实现符合 4–5 位）；
  纯声子族（MgB₂/H₃S/LaH₁₀）`Λ^sf≡0` 结构性覆盖，同样不付代价。一处已撤回的错误方向核对：
  误把非对角截断读成谱截断（`τ_max` 一度错报 2.98，订正为 4.07——`θ` 的宗量是松原频率非能量截断，已逐字核对原文）。
- **标量引理 (G) 的部分结果**（笔记 `ceiling/lemma_G_partial_bounds_2026-09-18.md`，`lemmaG` 会话）——目标
  `|λmin(𝒪[f])|≤λ₁(𝒪[f])`（c=1）**未证**；已得：恒等式 (Z)、定理 G3（c=3）、定理 GF（`|λmin|≤max(f(1),2C_f)`）
  **analytic**；定理 GA `c=2501567/1174507≈2.1299` **rational-exact（7 格 / 13 个有理数）**；阈值 `ϖ≥1.0496`、
  barrier 1.296 为 numerical-evidence；no-go（**rational-exact**，CERT-G2+CERT-B）：「逐元绝对值 / Perron」整族证法
  给不出 `c<1.106366`。注意：无条件定理（本节首条）**不消费 (G)**——(G) 只服务已被证伪撤回的 (G′-abs) 路线
  （规范文件 §5 明示「本篇与两篇存档均不依赖」GA），故其未证不构成缺口。
- **未复核项（规范文件 §5 明示，不得引用）**：上条定理 GA 与 no-go 的 CERT-G1/G2 主会话尚未从零复核；
  MgB₂ `1.017/0.448/0.213/0.155` 四个数（转引未读原文）；Essenberger PRB 94,014503（付费墙，挂用户）。

### 4.9 零参数一行公式：C_∞√(η/M)（η/M 定理落地 + 独立材料复核，新方向 2026-09-18）

- **规范来源**：`theory/derivations/P4_eta_over_M.md`（定理 A，分支 `etaM-corrected`@`ed6de52`，已并入
  `p4-conj2`——见 §5 表 C 行，2026-09-19 复核订正）+ `P4_main_theorems.md` 定理 1(c)/定理 2 + `P4_g2_analytic.md`。
  一行公式 `T_c ≤ C_∞√(η/M) = C_∞√(λ⟨ω²⟩)`，`C_∞=(1/2π)√𝔤(2)=0.1827262477…`，两种写法是**同一条定理**，
  `λ⟨ω²⟩=η/M` 是**恒等式**（McMillan 1968 Eq.22 重排），不是近似。
- **独立复核件**（笔记 `ceiling/oneliner_zero_parameter_2026-09-18.md`）：逐项标级——定理陈述 **analytic**；
  `λ⟨ω²⟩=η/M`（固定 M，单组分）在 (A1)(A2)(A3) 下 **analytic**；`C_∞` 数值 **rational-exact，不得简称 analytic**
  （CERT-B：12 格/约 18 个有理数，宽 3.964e-5；CERT-D：200 格/同法，宽 4.634e-10，取上端 0.1827262477890）；
  材料表 `numerical-evidence`。**四条承重限定**（缺一都会读错）：①单边界方向（上界非估计）；
  ②μ\* 地位——定理对**一切 μ\*≥0** 成立，但"达到上界的谱是 Einstein 谱"只在 **μ\*=0** 严格成立
  （μ\*>0 下 Einstein 极值性在弱耦合近抵消区 λ≲2.5μ\* 已被证伪，两句不可合并成一句）；
  ③ (A2)(A3) 只在把结果读成"η 与 M 无关"（同位素指数 ½ 的解读）时才进场，非谐下会漂移（PdH 实测 30%）；
  ④ 单带/各向同性框架，套多带材料须换 §4.8 的定理并核 (Z-dom)。0/48 常压氢化物候选独立重算 0 违反
  （`numerical-evidence`，`_verify_48_oneliner_2026-09-18.py`）。
- **材料表两轮独立复核**（`materials/etaM_material_table_recheck_2026-09-18.md` + 主会话再复核
  `materials/etaM_table_mainsession_reverify_2026-09-18.md`）：H3S/LaH10(4 压力点)/YH6/PdH 七行 η/M、界A、界B、
  A/exp **全部独立重算相符**（各自独立手写算子，未共用代码），仅两处**未披露的转引**值得关注——
  YH6 的 `Tc_thy=241.5`/`μ\*=0.125` 是区间中点非原文单点值；LaH10@163GPa 理论值配了 150GPa 的实验值——
  两者不影响"上界"性质，只影响"紧度"类比较性陈述的精确解读。严格性：全部 **numerical-evidence**
  （float64 点算，η/M 材料输入来自外部 DFT，未区间化）。
- **独立第三方复核**（笔记 `theoremI/Cinf_oneliner_independent_recheck_2026-09-18.md`）——文件自标
  **analytic（链条）+ rational-exact（`C_∞` 包围）**：判决表判定链条走通、无缺陷；并指出 CERT-B 只稳到 5 位，
  凡引用超过 5 位的 `C_∞` 数字必须点名用的是 CERT-D，否则位数不可信。
- **优先权风险核查**（笔记 `literature/leavens_allen_dynes_1979_priority_check_2026-09-18.md`）——核查
  Leavens–Allen–Dynes, SSC 30, 595 (1979) 是否已给出与 `P4_g2_analytic.md`/`P4_einstein_extremal.md` 同等
  或更强的结果；三个独立数据库（ScienceDirect/OSTI/Semantic Scholar）均未能取得全文或摘要，Crossref 给出
  完整参考文献列表（一手-元数据，该笔记专用标注，非仓库四档），P.B. Allen 个人主页精选出版物列表未收录该文。
  ~~结论：全文未能获取，优先权风险未能证实也未能排除——留作待办~~
  **【2026-09-25 订正：以上为 09-18 定格，已被同一笔记的 §8 取代，且与本文件第 82 行自相矛盾】**
  **缺口 1 已于 2026-09-20 关闭**：LAD79 一手全文到手（**用户直接提供 PDF**，存 `../notes/_scratch_ir/LAD1979.pdf`，2 页）。
  ⚠ 该 PDF 是 TextBridge OCR 扫描件，**文本层完全不可用**（实测输出为字符乱序的粘连串）⟹ 全部逐字内容取自**页图**
  `_lad1979_p-1.png` / `_lad1979_p-2.png`——这是「两种仪器都要跑」的又一个实例。
  **判决（三段，不许压成一句）**：
  ① **常数与公式形式是 LAD 的**——`0.1827` 与一行式 `kTc = 0.1827√(η/M)` 都出自 LAD79 Eq.(10)（本文件第 82 行已如实标注）；
  ② **但方向相反**：LAD 的 Eq.(8) 代数解出的是 `k_BT_c ≥ (1/2π)√(λ⟨ω²⟩ − ω²_max)`，是 Tc 的**下界（地板）**，
  且来源是**数值解 gap 方程**、只刻画「极限成立的 range of validity」；他们**从未**声称 `0.1827√(λ⟨ω²⟩)` 是对一切 λ、一切 P 的**上界**；
  ③ ⟹ **全局天花板 `T_c < C_∞√(λ⟨ω²⟩)`（`𝔣 ≤ 1`，∀λ>0、∀P∈𝒫）仍然是我们的**，「首次严格证明」措辞不需要改——
  **但理由变了**：新颖性只剩「不等号方向（上界）+ 全 λ 全谱 + 常数被认证包围」三点，**不含常数本身**。
  副产物：LAD 的地板可与我们的天花板配成 `𝔣` 的双侧夹（该侧 **numerical-evidence**，且其前提 `λ(ν₁)≤1` 的证明在 AD75、**我们未独立复核**）。
  详见该笔记 §7–§8；脚本 `../notes/_scratch_ir/_chk_main_lad1979_floor_2026-09-20.py`、`_chk_main_lad1979_table1_2026-09-20.py`、`_chk_main_lad1979_invert_2026-09-20.py`。
  **★ 本行过期了五天没人发现，共因与 `P4_einstein_extremal.md` 那条「推论 μ\*」未登记相同**：
  同一文件里第 82 行（新）与第 167–171 行（旧）各说各话，而**谁也没去核对另一行** ⟹ 改某篇稿的结论时，必须 grep 本文件内**所有**提到它的行，不只改自己那一行。

### 4.10 λ_ij 判据：Hadamard 平方 ≠ Gram 形（新方向 2026-09-18）

- **机制**（笔记 `ceiling/lambda_ij_criterion_mechanism_and_CaC6_2026-09-18.md`）——**analytic**：μ\*_ij 的
  Bochner/PSD 论证（第一性原理正核构造 ⟹ 恒 PSD）**不能**外推到 λ_ij：声子耦合 `λ_ij ∝ |g_ij|²` 是耦合矩阵
  `[g_ij]` 的 **Hadamard 平方**，不是 Gram 形；纯带间反例 `g=[[0,g₁₂],[g₁₂,0]]` 一行给出 PSD 判据违反，
  失效面**原子型**；根因是动量守恒令 `g` 只依赖 `k−k'`（Toeplitz 结构），`|g(q)|²≥0` 不蕴含 (k,k') 上 PSD。
  判据的约定无关形式是 `λ11λ22 ≥ λ12λ21`。**消费点已钉死**：该判据只在库仑侧引理 μ-1 被消费，声子侧从未
  被天花板链条（§4.8）实际用到——材料侧的通过/违反因此**在链条外**，不影响 §4.8 定理本身。
- **材料证据**（同文件 + `ceiling/lambda_ij_criterion_fails_on_FeSe_2026-09-18.md`）——**numerical-evidence**
  （文献转引 λ_ij 数值 + 简单算术）：MgB₂ 通过、余量 13.8×（源头一手核对文本层）；**FeSe 违反 2.48×**；
  CaC₆（常规电声六带）**11/15 对违反、最劣 2.438×**。铁砷化物一手带分辨 `λ_ij` 取件
  （笔记 `literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md`，四篇 harvest 全文 L1）**全部拿不到**——
  三篇 Linscheid 不是铁基论文（当初以作者组推断为"最便宜路线"是误判，已在文中自认）、Lischner FeSe 只给对称性
  投影标量非带分辨矩阵；替代检索（Boeri 2008 LaFeAsO）只有未分带总 λ=0.21。**自由度纪律**（新规则，
  §12.2 同文件）：判据数值前须先数生成式自由度——如 MgB₂ 的四个拟合 μ\*_ij 实际只有三个自由度。
- ⚠ **多次"订正别人的错"里新陈述被免检放行**（CLAUDE.md W1 热区）在本文件内被作者自己计数发生**十余次**、
  逐次标注并订正——本索引不复述过程，仅记最终判决（如上两条）。

### 4.11 变号能隙 / Δ 通道顶角修正：Perron–Frobenius 定理（新方向 2026-09-18）

- **「变号能隙需要 𝒫 之外物理」推论证伪**（笔记 `vertex/spm_gap_needs_no_exotic_physics_2026-09-18.md`）——
  **analytic**：原推论断在第一步——定理 I（配对通道免疫）从未对动量/带指标的符号发言，
  `γ>2∩𝒫=∅` 谈的是 Matsubara 差指标 k 上的量，两者都不是"能隙变号"（带指标符号）这件事，
  是 CLAUDE.md「一个符号只许一个意思」纪律要防的对象混淆。**替代结果（正定理）**：带内声子 `Λ^ph≥0`
  逐元非负 + 带间自旋涨落 `Λ^sf≥0` 逐元非负 + 耦合图不可约（三条前提，第三条是本文件经三轮订正才补全的
  必要前提，缺失时有反例）⟹ Perron–Frobenius 结构性必然给出顶本征向量 band-sign-changing（s±）——
  不需要两带对称、不需要定理 I 的任何引理。该文件经历三轮自我订正（脚本 bug、前提遗漏两次），
  末版判决为独立于团队复核脚本的自证 analytic。
- **顶角侦察报告复核 + 结构性升级**（笔记 `vertex/vertex_sign_recon_verification_2026-09-18.md`）——
  逐条复核 `vertex_sign_theorem_recon_2026-09-17.md`：§5.6 通道表/McMillan 导数 **逐位复现**；
  报告"问题 2"（线性化后核仍是四指标张量）与"问题 4"（`γ^Δ` 非凹）**均为假**（sympy 精确求导：
  `∂²[Δ/√(ω²+Δ²)]/∂Δ²` 对一切 ω≠0、Δ>0 严格为负，`γ^Δ` 严格凹；报告把 `γ^ω` 的性质错安到 `γ^Δ` 头上，
  且凹性本该问的是声子测度变量 `ϖ²` 而非 Δ）；总判断"不在 P4 算子类内、先做版本 1"仍站得住。
  顶角矩阵 `K^V_Δ`（Schrödi FSR factorized 核线性化，对称矩阵 `-Λ∘(C+Cᵀ)+D`）**660/660 负定**
  （`numerical-evidence`，float64 `eigvalsh`）⟹ 一行变分给出 Δ 通道顶角修正压低 Tc（**analytic**，相似变换
  Weyl 单调性）。
- **Δ 通道 no-go 推到同位素指数（条件式预言）**（笔记 `vertex/nogo_to_isotope_conditional_2026-09-18.md`）——
  条件式 **analytic**：**若** `K^V_Δ≼0` 且质量依赖只经 `Ω_ph(M)∝M^{-1/2}` 进入（标准绝热假设，非本仓已证命题），
  **则** `Δα_Δ-vertex<0`；`δ_D/δ_H≈√(M_H/M_D)` 数值支持 **numerical-evidence**（本机独立重跑
  `theory/derivations/p2_schrodi_fsr.py` 复现，非转引）。**关键反直觉发现**：该负位移方向与 PdH 亏损预算的
  系统性过冲**同向叠加而非抵消**——加入顶角修正只会让过冲更大，不能用来解释它。量级本身承认不可信
  （factorized 核已知过估 3–4×，本文实测部分组态达 8×）；PdH 论证的"顶角效应预期很小"未核实外部数字，
  标注需用户授权外部检索。
- **支持性文献取件**（笔记 `literature/source_hunt_seven_items_2026-09-18.md`）——七项文献取件任务，
  含 Klein 1977 μ\*=0.085 的 ω_c 约定核查（**结论：该值无 ω_c 约定**，APS harvest 路线全文 4/4 命中）；
  其余各项标注**亲读/转引/推断**，未汇总单一头条结论。此文件同时是同位素判别检验相关背景材料的来源
  （详细结论在 `theory/notes/isotope/`，本次任务范围排除该目录，此处只记指针）。

### 4.12 材料侧亏损预算与数据集普查（新方向 2026-09-18）

- **f_pred vs f_obs 零参数化尝试**（笔记 `materials/deficit_budget_fpred_vs_fobs_2026-09-18.md`）——三轮自我订正
  （原始版本混淆"独立确认"与共享系统误差、混淆"交叉验证"与代数恒等式、传播公式反向）后判定 **PASS**：
  H3S/D3S 共享同一份 α²F 工作流，不构成两次独立检验，需正交分解为**共同偏移 +7.44%**（几何平均，样本量 1，
  ≈3.9σ，实验散布解释不掉但不能反推模型被证伪——系统误差来源未量化）+ **独立同位素比值检验 −1.24%**
  （`numerical-evidence`，落在噪声带内 0.33σ，不构成"最紧的数字"）。严格性：全部 **numerical-evidence**
  （传播估计经三轮口径统一，最终版本用 RMS 合并 4 对 D3S 数据独立复算 σ(ratio_obs)=3.81%）。
- **H3S 实验 Tc 压强赋值审计**（笔记 `materials/h3s_tc_exp_pressure_audit_2026-09-18.md`）——两篇一手文献
  （Drozdov Nature 2015、Errea PRL 2015）逐页读图核对，级别 **L1**（本文专用标注，非仓库四档）。判定 1
  （压强不匹配已被原笔记自陈，非隐藏缺陷）为**文本核对**（文献层面，不套四档）；发现 Errea Table II 的 190K 锚点
  源自被 Drozdov 正式版点名撤回的未退火预印本数据——本轮用 PyMuPDF 直接解析 Drozdov Fig.2(c) 内嵌矢量图元
  重新数字化出干净的退火后锚点，取代有争议的 190K，脚本已 `git add -f`。
- **Cerqueira 6912 材料数据集 batch-q 剪枝/人口扫描**（笔记 `materials/cerqueira6912_batchq_pruning_population_2026-09-18.md`）——
  只覆盖 17 个 batch 中的 1 个（`batch-q`，250 个材料，249 个含完整 `McMillan.dat`），选择理由是体积最小
  非材料属性挑选（非随机 λ 分布偏置需注意）。**关键更正**：`McMillan.dat` 直接打印 `w2av`（⟨ω²⟩，与
  一行公式 Eq.1c 定义一致）——此前 2026-09-17 笔记「ω̄₂ 未确认存在于可下载文件」的判断**不成立**，
  本文档下载解压后独立发现。落盘 `theory/data/external/cerqueira6912_batchq_sample.csv`（249 行）。
  严格性：**numerical-evidence**（外部数据集直接转录 + 简单解析）。
- **μ\*>0 弱耦合反转经验命中率**（笔记 `materials/dir2_population_overlap_2026-09-18.md`）——判据
  `Ξ:=λ/μ\*(ω_log)`（`Ξ<2.6` 反转候选 / `Ξ>3.1` 安全）领头阶 **analytic** + 精确算子两点下界
  `numerical-evidence`；本文全部人口统计数字 **numerical-evidence**（有限本地样本）；物理 μ\*=0.10–0.13 区间
  该判据"没有测量"，只有拟合外推，且实测阈值比外推线高 3.4–9.0%（方向：真实阈值不会更低）。
- **多带 μ\*_ij 矩阵收割（MgB₂ 快照）**（笔记 `materials/mustar_psd_materials_2026-09-18.md`）——**阴性头条结论**：
  已发表多带 Eliashberg 的 μ\* 只有三种形态（标量、正核构造、整体设零），全部落在 Gram 保护伞下 ⟹ μ\* 通道
  不构成 §4.8 天花板的例外；主流做法（Choi et al. PRB 66,020513）用**标量** μ\*，非矩阵；两组转引比值
  （Mazin/Moon 均未重整化裸 μ）与真正 μ\*_ij 之间存在未做的重整化步骤。严格性（主会话 2026-09-19 裁定，
  与 `theory/notes/INDEX.md` 该行一致）：材料数字 **numerical-evidence**（文献收割 + 2×2 主子式算术）；
  阴性头条所依赖的引理 μ-2/μ-3 为 **analytic**（转引 §4.8 的 `ceiling/mustar_psd_multiband_ceiling`，非本文自证）。
  ⚠ **partially superseded**：该文末尾的「衍生结果 (G′-abs)」（`Tc ≤ C_∞√ρ(|M2|)`）**已被证伪**（见 §4.8），
  只有阴性头条 live。
- **η/M 材料表数据集覆盖普查**（笔记 `materials/etaM_material_table_recheck_2026-09-18.md` /
  `materials/etaM_table_mainsession_reverify_2026-09-18.md`）——见 §4.9，两份独立复核已并入该节。

### 4.13 非声子机制 / λ-晶格自洽 / μ\*>0 深挖（H_μ Tolmachev 透镜）（新方向 2026-09-18）

- **自旋涨落与晶格结构相容性**（笔记 `gamma/spin_fluctuation_structural_compatibility_2026-09-18.md`）与
  **自旋涨落严格 Tc 定界侦察**（笔记 `probes/spin_fluctuation_rigorous_tc_scoping_2026-09-18.md`）——均为
  **侦察级**：多数结论 `numerical-evidence`/文献转引，含少量 `analytic` 恒等式（各文件内部自行标注，
  未汇总出统一头条结论，**本索引不代入单一标级**，具体见各文件自身的逐条标注）。
- **λ 稳定性界 + KAY Tc(λ) 联合天花板**（笔记 `misc/lambda_stability_tc_joint_ceiling_2026-09-18.md`，
  已由主会话独立复核 `misc/lambda_stability_tc_joint_ceiling_review_2026-09-18.md`，裁决 **CONFIRMED**）——
  "记号同构"前提（KAY 的 P(dω) 与 Semenok λ\_\* 链是同一测度上的 Jensen 型泛函）确认为仓库既有事实
  （`P4_rigorous_tc_bound.md:25`），属**对既有一手推导原文的阅读复核**（文献层面，不套四档）。核心数值结论：对称
  两点谱上联合自洽解给出 `d ln Tc_joint/dσ²|₀ = +1.88`（谱展宽**提升**联合天花板，与"Einstein 谱总最优"
  直觉相反）——两条独立方法（直接自洽求解+外推 vs 隐函数定理闭式）5–6 位一致，**numerical-evidence**。
  偏度项 (γ₃) 已测（原头号未测风险已解除）：量级约为 σ² 项的 1/5，扫描范围内（σ 到 0.20、|γ₃| 到 4）
  符号未翻转；N=1→N=2 截断稳健（同号，差 1.3%）。**限制（构造固有，未被解除）**：这是 σ→0 局部微扰结果，
  非"Einstein 谱全局最优"的变分定理；构造本身是探索性的人为耦合，不推翻 U5（λ\_\* 与 Tc(λ) 上界逻辑独立）。
- **H_μ Tolmachev 透镜**（笔记 `coulomb/hmu_tolmachev_lens_2026-09-18.md`，已由团队复核
  `coulomb/hmu_tolmachev_lens_review_2026-09-18.md`）——**mixed** 严格性（文件自标）：§1.2/§1.3/§2.6/§4.1/§4.4
  为 **analytic**，其余全部 **numerical-evidence**，零 certificate/零 rational-exact；复核文件判定核心链条
  全部 CONFIRMED，另有 2 项因预算未复核、未计入通过率（详见文件自身 §0 十四项表）。
- **μ\*>0 五处失效有没有共同机制**（笔记 `coulomb/mustar_failure_common_mechanism_2026-09-18.md`）——
  **假说被证伪**：五处失效不能归到同一个控制量，分四类、阈值彼此相差最多 6.5×10⁴ 倍；但存在一个共同
  *充分条件*（顶本征向量 v∈锥 C，经同一 Rayleigh 步）是其中三处（跨窗下跳、窗内单调性、Einstein 极值性）
  的共同**祖先**——库仑项在这三个论证里逐项精确抵消，`analytic`（§1 久期恒等式、§3 命题 P1/P2/P3、§5.2
  领头阶判据、闭式 (★A)）；其余全部数字 `numerical-evidence`。真正的结构性张力：锥条件在 μ̃~N_c⁻³ 就早失效，
  远早于任何结论真正失效（超前 3.6–6.5×10⁴ 倍）——"锥不是免疫预算，是过期得极早的通行证"。此结论是
  §4.5"μ\*>0 Einstein 极值性失效边界"条目的上游统一视角，未取代其数值结论。

---

## §5. 六方向深挖轮（部分已并入 `p4-conj2`，部分仍 held）

> 这些是在 `p4-conj2` 之上、各自 `deep-*` 分支上完成、**均已个人核验闭环**的更深结果。截至 2026-09-15，另一会话已把其中若干**并入 `p4-conj2` 集成分支**（见下表「合并」列），其余仍 held 在分支上；**`p4-conj2` 整体仍未合并到 `master`**（仓库无 `main` 分支）。**权威状态、门数、订正历史见 auto-memory `p4-six-direction-round-2026-09-15.md`**（本表只给指针，避免与易变的分支状态重复而失同步）。合并状态用 `git merge-base --is-ancestor <branch> p4-conj2` 核验。
>
> 另有若干**基版 / 共享瓶颈分支**已并入 `p4-conj2`：`deep-universal-tc-envelope`（A 基版，§4.4）、`deep-coulomb-global-tc`（B 库仑版，§4.5）、`deep-concentration-bound`（D5/smooth-coulomb 共享瓶颈，§4.2）、`deep-b3-tower`（§4.3 矩量塔）。
>
> **合并后主会话在 318cfe2 上重跑门确认无回归（2026-09-15，`OMP_NUM_THREADS=1`）**：`p4_g2_analytic_bounds.py` 389/0 · `p4_universal_tc_envelope.py` 66/0 · `p4_global_tc_certified.py` 2018/0 · `p4_concentration_bound.py` 123/0 · `p4_b3_tower.py` 22/0 · `p4_coulomb_global_tc_certified.py` 78/0。合并为纯增量、无冲突、交付文件与源分支逐字一致。详见 auto-memory `p4-conj2-consolidation-reverify-2026-09-15.md`。
>
> **2026-09-19 复核**：逐行重跑 `git merge-base --is-ancestor <branch> p4-conj2`——**C**（`etaM-corrected`）现为真，
> 已改标「已并入」（见下表，合并 commit `9e1c83d`）；**B**（`deep-global-tc-allcoupling`）、**D4**（`deep-mustar-completion`）、
> **D5**（`deep-lemmaC-analytic`）、**D6**（`deep-anisotropy-enhancement`）四支复核结果均为 NOT-ancestor，表中「未合并」状态维持不变。
>
> **★ 2026-09-20：B / D4 / D5 / D6 四支经用户放行，已全部并入 `p4-conj2`**
> （合并 commit 依次 `a3ebb5d` / `0e75b0f` / `d8d716d` / `c4ac93a`，四次 `--no-ff` 全部纯增量、零冲突）。
> ⟹ **六方向深挖轮的分支 held 状态至此清零**；`p4-conj2` **本身仍未合并到 `master`**（按用户指令继续 hold）。
> **四个门由主会话在合并后的工作树上从零复跑**（`OMP_NUM_THREADS=1`，**不采信表内记录的计数**）：
> `p4_global_tc_allcoupling.py` **104/0**（146.1s）· `probe_Lstarstar_uniqueness_multiatom.py` **7 非单调 / 7 正确排除 / 0 假阴性**（45.2s）·
> `p4_lemmaC_residual_gate.py` **8 ck / 0 fail** · `p4_anisotropy_quantitative.py` **1280/0**（3.1s）。四项均与源分支记录一致。
> 连带两处记账修补：`P4D_coulomb_multiatom_status.md` 的三个被引探针原以 **worktree 绝对路径**引用（任何 clone 里都不存在），
> 已按脚本毕业规则 `git add -f` 入 `theory/notes/_scratch_ir/` 并改为相对路径；四篇新并入的稿件已补 §4 索引行。

| 方向 | 分支 @tip | 合并 | 头条结论 | 严格性 | 门 |
|---|---|---|---|---|---|
| **A** 普适速度极限 | 已并入主线 @318cfe2 | **已并入** | Tc < C_∞√(λ⟨ω²⟩)，C_∞=(1/2π)√𝔤(2)=0.18273，**所有 P、所有 μ\*≥0**（天花板**不依赖**包络门） | **解析**（单边） | 承重 `p4_g2_analytic_bounds.py` **389/0** + `p4_universal_tc_envelope.py` **66/0**〔memory 旧记的 `..._coulomb.py` 159/159 为已删 worktree 里的**未提交探针**、非承重、不可复现〕 |
| **B** 全局双侧 Tc(λ) | `deep-global-tc-allcoupling` @6a5572a | **已并入**（`a3ebb5d`，2026-09-20；门 104/0 已从零复跑） | 所有 λ>0 双侧 bracket；**深带经 Prop.C 闭式 (C1) 变为完全解析、零标定常数** | 深带**解析** / 其余认证 | `p4_global_tc_allcoupling.py` **104/104**（全，150.2s）已 team-lead 权威复跑 |
| **C** η/M「为何是氢」 | `etaM-corrected` @ed6de52 | **已并入**（commit `9e1c83d`，2026-09-19 用 `git merge-base --is-ancestor etaM-corrected p4-conj2` 复核确认为真；§4.9 已按合并态更新引用方式） | 定理 A：Tc < 0.18273·√(η/M)，η=N(0)⟨I²⟩ Hopfield，**μ\* 无关** | 定理 A **解析** | `P4_eta_over_M.py` 261/261（116.3s） |
| **D4** μ\*>0 唯一性/极值 | `deep-mustar-completion` @36ba77a | **已并入**（`0e75b0f`，2026-09-20；门 7/7 已从零复跑） | 多原子唯一性 + Einstein 极值机制逐字推广；L\*\*（数值探测）；L\*（κ≥1 近交叉）诚实 **OPEN** | **数值证据**（自标非定理） | `probe_Lstarstar_uniqueness_multiatom.py` 7/7 |
| **D5** 引理 C 凹性 | `deep-lemmaC-analytic` @2ff4d54 | **已并入**（`d8d716d`，2026-09-20；门 8/0 已从零复跑） | E(x)>0、C≤0、q′=E−|C| 解析；残余 **|C|<E 仅数值证据**（诚实边界，廓线上界为已识别缺失件，四条路线已排除） | 结构**解析** + 边界**数值证据** | 见 `P4_concentration_bound.py` |
| **D6** 各向异性增强 | `deep-anisotropy-enhancement` @5e5653d | **已并入**（`c4ac93a`，2026-09-20；门 1280/0 已从零复跑） | 定理 A：finite-N Bhatia–Davis λ_max 下界；定理 B：Tc 转换 | finite-N **解析** / N→∞ 数值 | 门 1280/1280（2.8s） |

---

## 附：配套目录

- **数值门脚本**：与各 `.md` 同名的 `.py`（如 `p4_g2_analytic_bounds.py`）；运行前 `export OMP_NUM_THREADS=1`。
- **求解器**：`theory/solvers/`（`eliashberg_iso.py` 全非线性求解器、`lambda_star.py` λ\_\* 泛函）。见该目录 README。
- **过程脚本 / 探针**：`theory/notes/_scratch_ir/`（git-ignored，非承重、含 rebuild 源）。
- **审查记录**：`theory/notes/*/*_review_*.md`（对抗核验审计轨迹；2026-09-18 起 notes 已按主题分子目录，路由见 `theory/notes/INDEX.md`）。
- **笔记语料**：`theory/notes/`（69 篇 2026-09-16/17 新增 + 33 篇 2026-09-18 新增，见 §4.8–§4.13）—— 证据轨迹、探索过程、订正记录。**权威状态在 auto-memory**（`MEMORY.md` + `memory/*.md`），笔记是原始材料不是结论。
- **文献 / 一手数据**：`library/`、`theory/data/`（各有 README/MANIFEST）。
- **总纲**：`theory/PROGRAM.md`（研究纲领叙事）；`README.md`（工作区总览）。
