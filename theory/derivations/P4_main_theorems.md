# P4 · 主定理整理稿（论文形态）—— KAY Conjecture 2 与 Einstein 极值定理的证明状态一览

建立：2026-09-12　状态：v0.1（主会话整理；所有内容均已在分散笔记中经 reviewer 对抗核，本文件只做**逻辑依赖图 + 统一陈述 + 状态标注**，不引入新论证）
用途：论文起草的骨架；逐条给出每个引理是「解析证明」「计算机辅助证书」还是「数值常数」，以及对应脚本与审查记录。
分散来源：`P4_conjecture2.md`（界半句、存在性）、`P4_uniqueness_strategy.md`（唯一性）、`P4_einstein_extremal.md`（极值定理）、`P4_interchange_uniformity.md`（谱隙/换序机器）、`P4_rigorous_tc_bound.md` §P4-C（𝔤(2) 包围）。审查：`notes/2026-09-1{0,1,2}_*review.md` 共 9 份。

---

## 0. 设定与记号（KAY II，逐字对齐）

- 谱测度 P∈𝒫（ℝ₊ 上概率测度，P((0,∞))=1，⟨ω²⟩:=∫ω²P<∞）；ϖ:=ω/2πT；u_n:=(2n+1)^{-1/2}；⟦k⟧(ϖ):=ϖ²/(k²+ϖ²)。
- 核-通用算子 𝒪[f] := −L1[f]+L2[f]+L3[f]，L1 对角 (2n+1)^{-1}Σ_{k=1}^{n}2f(k)，L2 非对角 f(|n−m|)u_nu_m，L3 f(n+m+1)u_nu_m。**线性于 f**。
- Einstein 算子 ℌ(ϖ):=𝒪[⟦·⟧(ϖ)]，h(ϖ):=λ₁(ℌ(ϖ))；一般算子 𝔎(P,T)=𝒪[F_T]，F_T(k):=∫⟦k⟧(ω/2πT)P(dω)（逐元 P-平均，KAY II Remark 2）；𝔨(P,T):=λ₁(𝔎(P,T))（最大本征值，KAY II Thm 3/eq.9）。
- Tc 判据：λ𝔨(P,T)≥1 ⟺ 正常态线性不稳定；𝔤(2):=ρ(𝔊(2))，𝔊(γ)=𝒪[k^{-γ}]。
- **结构恒等式**（本项目发现，逐元验证 + 代数证明）：(S1) 𝒪[1]=2uuᵀ−I（u∉ℓ²）；(S2) ℌ(ϖ)/ϖ²=𝔊(2)−ϖ²𝒪[k↦1/(k²(k²+ϖ²))]（所有 ϖ>0）；(S3) Lorentz 加法 ℌ(ϖ)=2ττᵀ+𝓡−𝒟，τ_n=u_nG((n+½)/ϖ)∈ℓ²，𝓡≥0 逐元，0≤𝒟≤1 对角。

## 1. 主定理

**定理 1（KAY II Conjecture 2，完整，μ\*=0）.** 对所有 P∈𝒫、λ>0：
(a) [存在] 存在 T>0 使 λ𝔨(P,T)≥1；(b) [唯一] T↦𝔨(P,T) 于 (0,∞) 严格递减、连续、T→∞ 时 →0；故临界面与线 ℒ(λ,P) **恰交一次**，T_c(λ,P) 良定义；(c) [界] T_c(λ,P) < T̃_c(λ,P):=(1/2π)√(𝔤(2)⟨ω²⟩λ)；更强地，对所有 T≥T̃_c 有 λ𝔨(P,T)<1。
状态：**(a) 解析（一页）；(b) 解析（两条初等组合引理 + 锥不变性）；(c) 半页解析化约 + 引理 B（计算机辅助证书 (0,6] ∪ 纯解析 ϖ>2.08）。** 三部分各经一轮对抗核，无 CHALLENGED。

**定理 2（Einstein 极值，全 λ）.** 对所有 P∈𝒫、T>0：𝔨(P,T) ≤ h(ϖ_rms)，ϖ_rms:=√⟨ω²⟩/2πT。等价：固定 λ 与 ⟨ω²⟩，Einstein 谱 δ_{√⟨ω²⟩} 是线性化 Eliashberg Tc 的极大者；T_c(λ,P) ≤ T_c^E(λ,√⟨ω²⟩)。定理 1(c) 是其 λ→∞ 推论。
状态：**两行解析化约 + 引理 C（凹性）：证书 (0,45]（服务器到 120）∪ 解析路线 R ϖ≥8.52**；解析段含 52 个显式常数，28 个区间包围、24+12 个 float64+自收敛余量（M3，部分严格化；M4 进行中）。三轮对抗核，无 CHALLENGED。文献定位：陈述为 Allen–Dynes 1975 folklore（Moussa–Cohen 2006 转述），严格证明新。

**推论.** (i) Einstein：h 于 (0,∞) 严格递增（KAY III Thm 5 的 [0,ϖ_*) 推广到全 ℝ₊）、H(x)=h(√x) 于 [0,∞) 凹；(ii) 弱耦合常数：h(ϖ) ≥ ln ϖ+γ+ln4−3/2−1/(12ϖ)−3/(16ϖ²)（首个严格下半边）；(iii) 𝔤(2)∈[1.318140211762671, 1.318140212136387]；(iv) μ\*≥0：定理 1(c)、定理 2 的上界对库仑束逐字成立（K_c⪰0 算子序）；(v) 一般-P 收敛 𝔨^(N)↑𝔨 带显式尾率（KAY 只有散文）。

## 2. 引理与依赖图

```
定理1(a) ← 引理X（存在）           ← (S1) + 压缩单调 + 控制收敛              [解析]
定理1(b) ← 引理A + 引理B + 锥论证    ← Abel 分解、Metzler 半群、谱投影强极限     [解析]
定理1(c) ← 引理J（Jensen）+ 引理B'   ← 段I证书 (0,6] + 段III解析 ϖ>2.08         [混合]
定理2    ← 引理J + 引理C（凹性）     ← 段C-I证书 (0,45] + 路线R解析 ϖ≥8.52      [混合]
```

| 引理 | 陈述 | 证明形态 | 脚本 / 检查 | 审查 |
|---|---|---|---|---|
| **X**（存在） | 𝔥^(N)(∞)=2S_N−1→∞ ⟹ 每 λ>0 有穿越 | 解析 | — | uniq-review 顺带核 |
| **A**（矩形部分和） | ⟨e_L,𝒪[δ_k]e_{L′}⟩≥0，取等 iff L+L′≤k | 解析（三组逐项非负） | `p4_uniqueness_lemmas.py` 37 | uniq-review PROVED |
| **B**（Kamke） | 𝒪[1_{k≤K}] 阶梯基 Metzler，取等 iff K=0 或 L≥L′+K+1 | 解析（恒等式 u_n²−u_{n+1}²=d_n(u_n+u_{n+1}) + 凸性） | 同上 | uniq-review PROVED |
| 锥论证 | 半群保锥 ⟹ 顶本征空间含锥向量 ⟹ Rayleigh 比较 | 解析（不需单重/Krein–Rutman/紧支/孤立） | — | uniq-review（两处文字订正已并入） |
| **J**（Jensen） | 𝔨(P,T)≤∫λ₁(ℌ(ϖ))P(dω) | 解析（逐元积分 → 二次型 → sup） | — | conj2-review CONFIRMED |
| **B′**（Einstein 逐点） | h(ϖ)<𝔤(2)ϖ²，ϖ>0 | 段 I：恒等式 (S2) + 锐化二阶余项（𝔤(2) 两边相消）证书 (0,6]；段 III：Schur 行和 ‖ℌ‖≤(⅓+2/√3+1/√15)S₁≤Cπϖ/2 解析 ϖ>2.08 | `p4_conj2_small_w.py` 540、`p4_conj2_interval.py` 1878(+1656) | conj2-review CONFIRMED |
| **C**（凹性） | H(x)=h(√x) 于 (0,∞) 严格凹 | 段 C-I：(C\*\*) H″≤−2F+2x²E²/g 无相消写法，证书 (0,45]（→120）；路线 R：(S3) + Θ\*=0.6349 隙界 + 显式试探向量 + DK，闭式判据 ϖ≥8.52 | `p4_extremal_concavity.py` 1940(+3253)、`p4_cii_m1.py` 219、`p4_cii_m2.py` 443、`p4_cii_m3_interval.py` 137 | concavity-/m1-/m2-review CONFIRMED |
| 𝔤(2) 包围 | 见推论 (iii) | Rayleigh 下 + Kato–Temple 上，mpmath | `p4_conj2_small_w.py` Stage 2 | conj2-review |
| 换序/谱隙机器 | δ(ϖ)>0、inf_N gap_N≥g_alt、[0.2,60] 区间 | 证书 | `p4_interchange_uniformity.py` 1668 + 扩展 2209/5561 | 三轮 review |

## 3. 残余假设（对外必须前置）

1. **计算机辅助证书的约定**：mpmath dps≥30、显式舍入垫（1e-12/1e-25）、Cholesky 成功当作正定证书；非区间算术。涉及：引理 B′ 段 I、引理 C 段 C-I、𝔤(2) 包围。
2. **定理 2 解析段的常数**：52 个显式收敛积分/sup，28 个已 mpmath.iv 区间包围，24 个嵌套双重积分 + 12 个 pad/sup 目前为 float64 值加自收敛余量（M3 §1.H；M4 进行中）。非承重（ϖ=45 处对 <200% 单常数误差稳健），但形式上未闭。
3. 全部结果在 **μ\*=0**（KAY 猜想原文亦无 μ\*）；μ\*>0 只有上界推论 (iv)，唯一性不覆盖（库仑束 Perron 向量出锥；且 N_c(T) 阶梯使 f(ϖ) 可有下降步，Tc 须取 sup）。
4. 线性化 Eliashberg 层面；λ 的物理上界（U5）、μ\* 的第一性（U7）、有界峭度 U6 等系统性限制照旧。

## 4. 与文献的关系（措辞已定）

- KAY I–III：Conjecture 2 = KAY II p.14 eq.31 = KAY III Conj. 1 一般-P 版；化约 "elementary but not stated in [KAY II, III]"（最近 Remark 2、Prop 7）；单调性 KAY 只证 [0,ϖ_*)、自陈 "would require a different strategy"；一般-P 收敛与速率 KAY 只散文；§3.6.2「大 λ 已证」依赖未证的 N=∞ 展开，段 I 补证。
- Einstein 极值：Allen–Dynes 1975 folklore（Moussa–Cohen PRB 74, 094520 逐字），我方首个严格证明；勿写「无人陈述」。
- 竞争：Elezaby–Abanov PRB 2026 只做 γ-model（𝔤(2)≤1.545）；无人证猜想。
- 弱耦合常数 c₀=γ+ln4−3/2 与 Marsiglio PRB 98, 024523 eq.37、Chubukov et al. 2004.01281 eq.17 一致。

## 5. 论文骨架建议

1. Introduction（KAY 框架、猜想、主定理 1/2、方法概述：结构恒等式 + 锥不变性 + 计算机辅助段）
2. Setting and structural identities（§0；(S1)–(S3) 的证明）
3. Theorem 1(a),(b)：存在与唯一（引理 X、A、B、锥论证——全解析，可独立成节）
4. Theorem 1(c) 与 Theorem 2 的化约（引理 J；Einstein 逐点不等式与凹性的陈述）
5. Einstein pointwise inequality（段 I/III）
6. Concavity（段 C-I 证书 + 路线 R 解析）
7. Computer-assisted certificates: conventions, reproducibility（脚本表、检查数、复现说明、残余假设）
8. Corollaries and discussion（推论 (i)–(v)、μ\*、与 Allen–Dynes/Marsiglio 常数的一致、开放问题：μ\*>0 唯一性、完全区间化、λ 侧）

**建议下一步**：先按本骨架把 §3（纯解析部分）写成可投稿的定理-证明形态，因为它不依赖任何证书、可先行严审；§5–§6 待 M4 收口后定稿。
