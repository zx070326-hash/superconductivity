# 弱耦合久期级数（定理 W / W′）对抗审查记录

日期：2026-09-14　审查者：reviewer `weakcoupling-review`（sonnet）　落盘：主会话（reviewer 角色无写权限，报告文件未能落盘，本文按其两次消息全文转录，主会话只加「主会话处置」段）
对象：`theory/derivations/P4_weak_coupling_series.md`（合入 `p4-conj2` 的版本，含主会话对 W(a)、W(d) 的订正）与 `p4_weak_coupling_series.py`（37 检查）。

## 总裁决

- **定理 W（双侧包围）：CONFIRMED**，数学严格。
- **定理 W′（三项渐近）：CHALLENGED**，证明不完整——b₁ 的极限存在性未证。
- 数值侧：37/37 复现；reviewer 自写代码独立重算一致；b₁ ≈ 0.3584 ± 0.0001（4 位可信），文献首次。
- 无 showstopper。

## 逐条

1. **W(a) Woodbury 与唯一性**：数学正确，但叙述顺序有歧义（`:46-50`）——读起来像在证「λ₁(ℌ) 满足久期方程」，实际是在证「久期方程的根唯一」。**建议措辞**：「给定 λ₁(ℌ) 由 Woodbury 是 F=1 的根，则由 F′<0 它是 (λ_max(M),∞) 上唯一的根。」（DOWNGRADE，措辞）
2. **W(b) 逐项正性**：R(s,σ) ≥ 0 对所有 s,σ ≥ 0 成立（分子 2s²σ²(5+s²+σ²)，分母四个正因子）。CONFIRMED。附注：收敛条件 h > ‖𝓡‖ 偏保守，真正需要的是 h > ‖𝓡‖ − min 𝒟_nn；应用中 h ~ ln ϖ ≫ Θ\*，无害。（`:52-54`）
3. **W(c) 尾项几何级数**：逐步核对，代数正确。Θ\* 是认证常数（`P4_einstein_extremal.md` §1.F，λ₂ ≤ ‖𝓡‖ ≤ Θ\* 对所有 ϖ 成立，p-加权 Schur）。CONFIRMED。
4. **W′ 步骤 4→5（`:77-106`）**：δ = −[2⟨τ,(I−𝒟)τ⟩ + 2⟨τ,𝓡τ⟩]/(2B) + O(1/ln²ϖ) 到 b₁/ln ϖ = δ 这一步**要求方括号内两项在 ϖ→∞ 收敛到常数**，稿子只在 §3.A 给了标度极限的启发式论证，没有控制格点和与积分之差。数值收敛证据强（ϖ=100/200/400：0.358369 / 0.358396 / 0.358404，差 2.7e-5 → 8e-6），但这是验证不是证明。**CHALLENGED**。修法：Euler–Maclaurin 型余项控制，或把 W′ 标为「数值定理」。
5. **b₁ 的定义**：b₁ 为 1/ln ϖ 项系数，分母用 2B 还是 ln ϖ 的差别只进 O(1/ln²ϖ)（改变 b₂ 不改变 b₁）。CONFIRMED。
6. **数值门**：(M1) 恒等式 4.44e-16；包围在 K=0,1,2,4、ϖ∈[3,100] 成立；ϖ=3 端 2B−1 = 1.5526 > Θ\* = 0.6349 算术核对无误；独立重算 ϖ=10、N=800 与文中数字一致。CONFIRMED。
7. **物理翻译**：(ω_E/2π)e^{c₀} = 0.2530ω_E 与 Marsiglio PRB 98, 024523 一致；(1+b₁λ) 一阶项 reviewer 核对无误。CONFIRMED。
8. **文献比对**：稿子未沿用 survey 里「Conj. 2 等价于 b₁ 展开」的错句。CONFIRMED。

## 未发现问题的检查项

算术精度（手算到 4 位）、边界情形（ϖ=3、ϖ→∞）、尾界代数、上游常数 Θ\*（两处独立来源一致）、物理前因子。

## reviewer 建议的对外措辞

"We establish bilateral bounds via a secular series (Theorem W), and provide rigorous numerical verification of the weak-coupling expansion with b₁ ≈ 0.358, the first reported coefficient in the literature."

## 主会话处置（2026-09-14）

- CHALLENGED 4 的补法**正在进行**：opus agent `b1-closed-form` 的任务就是把 2⟨τ,𝓡τ⟩ 与 2⟨τ,(I−𝒟)τ⟩ 分别写成「格点精确和 + 标度积分 + O(1/ϖ) 或 O(ln ϖ/ϖ) 修正」，这恰是极限存在性的构造性证明；若它交付带余项界的闭式，W′ 升为定理，否则按 reviewer 建议标「数值定理」。
- DOWNGRADE 1 的措辞在 b₁ 结果合并时一并改。
- 第 2 条的保守收敛条件不改（无害，且简洁）。

---

## 附：reviewer 完整报告（第二次消息全文，主会话转录，2026-09-14）

**执行摘要**：定理 W CONFIRMED；定理 W′ CHALLENGED on proof completeness。37/37 全过。1 个措辞问题、1 个证明完整性 gap。无 showstopper。

**Q1 (M1) 恒等式**：CONFIRMED。w=7、N=300：4.440892e-16，独立复现相同；R 核、𝒟 对角、τ 归一化逐一验证。位置 `:24, 29-50`。

**Q2 W(a) 唯一性逻辑**：CHALLENGED → DOWNGRADE（措辞）。`:26-28`（陈述）、`:46-50`（证明）。陈述读起来像在证「久期方程成立」，实际在证「久期根唯一」；Woodbury 只说「根存在时分母消失」，实际逻辑是：已知 λ₁(ℌ) 存在且满足特征方程，证明它是久期方程的唯一根。数学正确，叙述顺序有歧义。建议：陈述改「λ₁(ℌ) 是久期方程 F=1 在 (λ_max(M),∞) 上的唯一根」，证明开头加「λ₁(ℌ) 满足 F=1 由 Woodbury 保证，下证唯一性」。**主会话已按此改（同日）。**

**Q3 W(b) Neumann 收敛条件**：CONFIRMED，偏保守。真实收敛域 h > ‖𝓡‖ − min 𝒟_nn；应用中 h ≥ 1.099（ϖ≥3）≫ Θ\*=0.6349，无害。

**Q4 W(c) 尾界**：CONFIRMED。R_K ≤ 2B·(Θ\*/h)^{K+1}·(1/h)/(1−Θ\*/h) = 2BΘ\*^{K+1}/[h^{K+1}(h−Θ\*)]，逐步正确，与 gate 公式一致。

**Q5 W(d) ϖ₀=3**：CONFIRMED（v0.1 订正后）。独立算 w=3、N=120：2B=2.5526，2B−1=1.5526 ≥ 1.5135 = ln3+c₀−1/36−3/144；> Θ\*。订正注里「0.51 > 0.6349」是描述旧错误时的笔误（可选修：改为「原公式多减 1/2，下界 1.01 而非 1.51」）。

**Q6 W′ 的 b₁ 推导**：CHALLENGED，证明不完整。`:77-106` 步骤 4→5：从 δ = −[2⟨τ,(I−𝒟)τ⟩+2⟨τ,𝓡τ⟩]/(2B)+O(1/ln²ϖ) 到 b₁/ln ϖ = δ，要求方括号在 ϖ→∞ 的极限存在且有限；§3.A 只有标度极限的启发式草图，未证格点和与积分之差 → 0，未证 ℓ_G(n) 在标度极限下的行为。需补：(1) Σ_n τ_n²(1−𝒟_nn) 与标度积分之差为 O(1/ln ϖ)；或 (2) 直接证 lim Σ_n w_nG(s_n)²(1−𝒟_nn) 存在。数值（gate `:209-233`）：0.358369 / 0.358396 / 0.358404（ϖ=100/200/400，Δ 2.7e-5 → 8e-6），4 位可信。置信：数值 CONFIRMED，解析 PLAUSIBLE but INCOMPLETE。建议标「Numerical Theorem W′」或补 Euler–Maclaurin 控制。严重性中等。

**Q7 b₁ 闭式策略**：CONFIRMED 可行。𝓡 部分被积函数 8 次有理、无奇点、两端可积；(I−𝒟) 部分格点和 + 标度积分。文献无人报告 b₁。

**Q8 数值门覆盖**：CONFIRMED。(M1) 1、包围 24、正性 5、单调 4、b₁ 分解 2、前因子 1。未覆盖（非 critical）：宽度单调只 print 未 assert；F′<0 未直接验证；注释 `:196-198` 可更明确。

**上游依赖**：Θ\*（认证值，`p4_cii_m1.py:145` 一致到 7 位）、2B 渐近（Marsiglio eq.37 + Chubukov eq.17，c₀ 拟合差 2.1e-5，排除 hard cutoff）、(M1)（1e-15）——全部可靠。

**独立重算**：ϖ=3 有效性、(M1)、单调收敛 (h−L−c₀)·L：0.24739 / 0.26905 / 0.28712 / 0.29648 / 0.30363（ϖ=10…200）。

**最终判决**：数值 PASS；解析 CONDITIONAL PASS——定理 W：A（措辞小瑕疵）；W′ 数值：A−；W′ 证明：B（有 gap 但数值强支持）。论文口径："We establish bilateral bounds via secular series (Theorem W), and provide rigorous numerical verification of h = ln ϖ + c₀ + 0.358/ln ϖ + O(1/ln²ϖ), where b₁ ≈ 0.358 is the first reported coefficient in the literature (verified to 4 significant figures)."

**修改清单**：必修无；应修 (1) `:92-106` 补极限存在证明或 `:69-70` 标「Numerical Theorem W′」——**主会话已加状态框，待 b1-closed-form**；(2) `:26-28`、`:46` 措辞——**已改**；可选 (3) `:42` 订正注措辞——未改。
