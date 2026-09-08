# H3S 电声耦合的 q 分布 + Mishra 2507 SM（researcher 摘录）

调研 agent：H3S q 分布 + 2507 SM（sonnet / researcher）　日期：2026-09-08。
用途：判 P2 §5.5 的"偏大 q 权重"假设是否有文献支持；查 Mishra 2507 SM 是否给 Z/φ 通道 λ^V 分解。

## H3S 耦合在 q 上的分布

- **Errea et al., PRL 114, 157004 (2015)**（"High-Pressure Hydrogen Sulfide from First Principles"，arXiv:1502.02832），p.3：H 振动的声子线宽（∝ q 分辨的电声耦合）**"fairly uniform throughout the spectrum"**——耦合在整个声子谱/布里渊区上大致均匀，**不是**小 q 峰化、也**不是**大 q 主导。作者明确对比：与 MgB2 那种"单一模主导"（σ 带 E2g 声子）**相反**。
- 推论：H3S 的 α²F(q) 没有强 q 结构。P2 §5.5 里"偏大 q 权重给负 P̄_eff"的构造在 H3S **没有文献支持**；均匀 q ≈ flat 相空间权重，对应 §5.5 表里 flat 权重那一行——P̄_eff **为正**。
- 所以：各向同性一阶顶角函数 P 本身在 H3S 的物理 q 权重下是正的，给不出 EPW 的负号。负号必须来自顶角函数之外的结构 → 两通道（§5.6，`p2_schrodi_fsr.py` 已确认）。

## Mishra 2507 (arXiv:2507.01897) SM

- SM 5 页，已取得（Springer static-content CDN）。
- FSR 对峰化 DOS 的放大机制：**"N_F enters explicitly in the numerator, amplifying the impact of a peaked DOS at ε_F"**——常数 DOS 的 FSR 把费米面 DOS 当全局尺度，尖 vHs 被线性放大；FBW 里能量依赖 DOS 把这个放大摊平，是 FSR/FBW 顶角效应差异（19.5% vs 5%）的来源之一。
- **SM 没有** Z 通道 vs φ 通道的 λ^V 分解，也没有把顶角净效应拆成两通道贡献的表——正文 Eq. 50–54 的 P̃^ω/P̃^Δ 矩阵是唯一的通道结构信息。所以 §5.6 的"Δ 通道主导"来自我们自己的 `p2_schrodi_fsr.py`（factorized 核），不是 Mishra 给的。

## 对 P2 的净影响

两条都指向同一结论、互相独立：(a) H3S q 均匀 → 顶角函数 P 正 → 各向同性路线给不出负号；(b) 负号在两通道结构。与 §5.6 的数值裁决一致。
