# Hubbard d-波配对 Track B：QMC 数据对接 + certificate 管线可行性（2026-09-16）

分支 `p4-conj2`，**未提交**。任务来源：main 会话指派（Track B）。
本笔记只做**调研 + 数据核验 + 管线可行性评估**，不产出新证书结论。

---

## 0. 结论先行（Phase 1 数据就绪度）

**判定：Blocked（核心 χ_d(q,ω) 数据）/ Partial（标量诊断量数据）。**

- 2026 年 normalizing-flow 突破论文（2603.18205、2601.18273）**没有公开数据或代码仓库**（无
  GitHub/Zenodo 链接），且论文本身聚焦小系统（8–18 site 六角格）验证方法正确性，**尚未报告
  大尺寸 χ_d(q,ω) 结果**——这正是"突破之后干净数据"这一前提目前不成立的地方。
- 找到唯一一个**真实可下载**的 2026 年数据集（Zenodo, Ghosh, DOI 10.5281/zenodo.19484945），
  但它是**传统 DQMC**（非 normalizing-flow，L=4,6,8，β=2–4，U/t=6，sign 0.85–0.999，
  不需要符号问题突破就能算），且只给**标量**诊断量（S_d 顶点、P_d 配对关联、χ_d、S(π,π)
  vs 掺杂），**没有动量分辨的 χ_d(q,ω)** 表格。
- 该数据集的 `processed_data/real_space/*_Pd_map.dat`（本应是 BS 核构造最需要的实空间/
  动量核信息）存在**数据完整性问题**（见 §2.3）：不能直接用于构造可信的 Bethe-Salpeter 矩阵。
- 因此：**没有找到任何公开的、动量/频率分辨的 χ_d(q,ω) QMC 数据表**可供 Track B 直接对接。
  Phase 1（1 周内上 certificate）不可行；即便是 Partial 路线（2–3 周清洗/外推）也缺少
  最基本的原始材料——不是清洗问题，是"数据还不存在"的问题。

---

## 1. 文献清单（2024–2026，normalizing-flow 突破前后）

| arXiv | 标题 | 日期 | 相关性 |
|---|---|---|---|
| [2603.18205](https://arxiv.org/abs/2603.18205) | Tackling the Sign Problem in the Doped Hubbard Model with Normalizing Flows | 2026-03-20 | ★核心突破论文。spin-basis + annealing NF，有限化学势下把符号问题的实际遍历性问题解决；benchmark 在 8/18-site 六角格，U=2,3，β=6,8。**无公开数据/代码**。 |
| [2601.18273](https://arxiv.org/abs/2601.18273) | Toward Scalable Normalizing Flows for the Hubbard Model | 2026-01 | Lattice2024 会议稿，讨论把 NF 方法推向更大格点/更低温的可扩展性路径，**未给出具体 U/t、格点尺寸的新结果**，无代码/数据仓库。 |
| [2604.01737](https://arxiv.org/abs/2604.01737) | Sign-Free Evidence for a d-Wave Superfluid Stiffness Dome in the Doped Hubbard Model (Wang & Lin) | 2026-04-03 | **非 normalizing-flow**，但同样是 2026 年"绕开符号问题"的技术（K_eff 矩阵对数平均法，符号扇区分解）。U/t=4, t'/t=−0.30，L=8,10,12（也测到 16，但 sector 分解在 L=16 失效），μ/t∈[−1.5,0]（掺杂 1–7%），β 最深到 20。给出 R_g（能隙比）、ρ_s（超流刚度）、S(π,π)，**无 χ_d(q,ω)**，**无公开数据仓库**。作者明确说明 vertex correction（配对通道的真正两体散射）仍未计入——ρ_s 是正常态单粒子谱的性质，不是配对本身。 |
| [2501.07371](https://arxiv.org/abs/2501.07371) | Simulating the Hubbard Model with Equivariant Normalizing Flows | 2025-01 | 更早的 NF-on-Hubbard 概念验证，1+1D 小格点（N_x=2），graphene 应用背景，非掺杂 2D 方格子超导问题，**仅方法学**。 |
| [Ghosh, EPJ B (2026)](https://link.springer.com/article/10.1140/epjb/s10051-026-01197-0) | Doping and temperature dependence of pairing correlations in the 2D Hubbard model: a DQMC study | 2026-06 | ★**唯一找到有真实可下载数据的论文**。传统 DQMC（QUEST 包），U/t=6，L=4,6,8，β=2–4，掺杂 x≈0–0.27。d-wave 通道识别为主导配对通道，vertex 随掺杂由负转正。详见 §2。 |
| [LeBlanc et al., PRB 107, 195138 (2023)](https://arxiv.org/pdf/2301.08176) | Pairing susceptibility of the 2D Hubbard model in the thermodynamic limit | 2023 | 非 2024–2026，但是**唯一给出热力学极限、动量分辨配对磁化率**的系统性工作（算法 Matsubara 积分，微扰展开）。若 2026 数据长期不可得，这是退而求其次的理论对照（非 QMC 原始数据，是微扰级数外推）。未在本轮深入提取，留作候选。 |

**符号问题突破的定位说明**：2603.18205 与 2604.01737 是两条**不同**技术路线，都发表于 2026
年、都号称"符号无关"：前者是训练一个生成模型直接采样有正权重的量（normalizing flow，
用户任务描述里指的技术）；后者是把 DQMC 传播子取矩阵对数再对符号扇区平均（K_eff 构造，
不训练神经网络）。两者都还没有报告 χ_d(q,ω) 的干净大尺寸数据表。

---

## 2. Ghosh (2026) Zenodo 数据集：提取与核验

DOI: `10.5281/zenodo.19484945`（真实下载验证，MD5 校验通过，296KB zip）。
本地路径：`theory/notes/_scratch_ir/hubbard_qmc_papers/ghosh_data/DQMC_Hubbard_Pairing_Ghosh_2026/`

### 2.1 数据结构

```
raw_data/L{4,6,8}/beta{2,3,4}/n{density}/L{L}U6beta{beta}_n{density}.out   — QUEST 风格原始输出
processed_data/dqmc_summary_full.csv                                      — 12 行汇总表
processed_data/real_space/{...}/L{L}U6beta{beta}_n{density}_Pd_map.dat    — 实空间 d-wave 关联图
plots/Fig1-6_*.pdf
supplementary_material/supplementary_material.pdf
```

参数：t=1, U/t=6, Δτ=0.1, L∈{4,6,8}, β∈{2,3,4}（T/t∈{0.25,0.33,0.5}），
掺杂 x=1−n ∈ [0, 0.27]，avg sign 0.33–1.0（β=4 时 sign 降到 0.33，符号问题已经在起作用，
但仍是传统 DQMC 可处理范围，**不需要 normalizing-flow 突破**）。

`dqmc_summary_full.csv` 列：L, U, beta, density, doping, avg_sign, S_xx, S_zz, S_total,
Sd（d-wave 顶点）, chi_d, s/s*/s**/d/d*-wave pair-field 关联。共 **12 个数据点**。

### 2.2 原始 .out 文件核验

抽查 `raw_data/L4/beta3/n0.777/L4U6beta3_n0.777.out`：格式与 QUEST DQMC 包标准输出高度吻合
（Number of sites=16=4², Time slice=30=β/Δτ=3/0.1，Pair-field correlation function 分 s/s*/s**/d/d*
五个对称通道，数值与 `dqmc_summary_full.csv` 对应行**完全一致**（0.51003172 vs d-wave=0.51003172,
误差 0.0032319264 vs err_d_wave 同值）——**汇总 CSV 与原始输出自洽**，这部分数据看起来是真实生成的。

### 2.3 数据完整性问题（发现，需标注）

检查 `processed_data/real_space/L6_beta3_n0.781/L6U6beta3_n0.781_Pd_map.dat`：

- 坐标范围 r_x, r_y ∈ [−30, 30]，**与 L=6 lattice（周期边界，物理上 r 应被约束在格点范围内，
  典型约 [−3,2]）完全不匹配**。
- 更严重：对比 `L8_beta3_n0.782_Pd_map.dat`，**坐标范围同样是 [−30,30]**——两个不同格点尺寸
  （L=6 vs L=8）给出了完全相同的实空间关联范围。这在周期性 L×L 格点上是**不可能**的（真实
  DQMC 实空间关联最多到 r=L/2 附近就应按周期性折叠或截断）。
- 3721 行中仅 55 行非零，非零值集中在 r_x∈[−13,0] 附近的窄带，衰减形态可疑（不像标准
  DQMC pairing correlation 应有的、以晶格矢量为周期的分布）。

**判断**：`Pd_map.dat` 系列文件**很可能是占位符/通用模板生成，而非从原始 DQMC 组态真实算出
的实空间关联**——这是 Track B 最需要的、用于构造截断 Bethe-Salpeter 核的原始材料，但恰恰是
这部分数据不可信。`dqmc_summary_full.csv` 里的标量（Sd, Pd, chi_d 等）与原始 `.out` 文件自洽，
可信度较高，但**这些都是 q=0（均匀）标量**，不含动量分辨信息，无法重建 λ(k−k′) 核。

**严格性标注**：以上是**numerical-evidence 级**的数据审计结论（人工抽查 + 一致性脚本核对，
未穷举全部 12×若干实空间文件），但坐标范围不随 L 变化这一点是**确定性矛盾**（两个文件
的字面内容直接比较），不依赖数值精度，可视为**analytic 级**的证伪。

---

## 3. 对照我们的 certificate 工具链

读过 `theory/notes/c1sq_schur_improvement_2026-09-16.md` 与
`theory/notes/_scratch_ir/c1sq_schur_bound.py` / `c1sq_lower_bound.py`。

**现状**：我们的证书管线是为**1D 理论模型**（ℓ²(ℕ₀) 上的 Einstein 核
f(k)=ϖ²/(ϖ²+k²)，权重 u_n=(2n+1)^{-1/2}）构造的，核心步骤：

1. 截断矩阵 H_N（N×N，iv 精确构造，闭式求和 a_n）。
2. float64 求 H_N 的近似特征向量 v1，Rayleigh 商给 μ1 的**certificate 下界**（无需尾控制）。
3. 区间 LDLᵀ + Sylvester 惯性计数，给 λ2(H_N) 的**certificate 上界**。
4. 尾部扰动 E=H−A 分块为 [[0,C],[Cᵀ,H^BB]]，‖C‖、‖H^BB‖ 各自解析/Schur 检验给 certificate 上界。
5. Davis-Kahan sinΘ 组合，给 c₁²=|ψ_1(0)|² 的 certificate 下界。

**可直接复用的部分**：
- **区间 LDLᵀ + Sylvester 惯性计数**（`inertia_count_below`）：与核的具体形式无关，只要求
  能把候选核写成 iv 精确矩阵——**可直接复用**于任何截断后的有限维对称/厄米矩阵，包括
  Bethe-Salpeter 截断矩阵。
- **Rayleigh 下界 + Davis-Kahan sinΘ 组合框架**（`rayleigh_lower`, `residual_norm`,
  `run_one_dk` 骨架）：给"截断矩阵最大特征值/特征向量分量的 certificate 下界"这一**通用模式**
  可以直接搬过来，只要能把 λ₁(核) 关联到某个物理量（比如 top pairing eigenvalue 对
  d-wave form factor 的投影）。
- **Schur 检验的方法论**（利用逐点尾部衰减而非全局算子范数夹逼）：**方法可以照搬**，但
  具体的 Tail(J) 解析闭式（arctan 有理下界、1/j² 伸缩）是**专门针对 f(k)=ϖ²/(ϖ²+k²) 这个
  Lorentzian 核**推出来的，换核要重新推导。

**需要重写的部分**（工作量大，是主要瓶颈）：
- **1D → 2D 动量空间求和**：Hubbard BS 核是在 2D 布里渊区 k=(k_x,k_y) 上定义的
  λ(k,k′) ∝ χ_d(k−k′) · g_d(k) g_d(k′)（d-wave form factor g_d(k)=cos k_x − cos k_y），
  离散化后是 N_k×N_k 矩阵（N_k = 布里渊区网格点数，不是我们现有代码里的"半轴序号 n"）。
  当前代码的 u_n、f(k) 全部是 1D 序号的函数，**没有可直接套用的 2D 网格模块**。
- **χ_d(q,ω) 核本身没有闭式**：Einstein 核 f(k)=ϖ²/(ϖ²+k²) 的闭式求和（fhat(0), fhat(π) 等）
  依赖它是一个具体解析函数；真实 χ_d(q,ω) 是**QMC 数值表**（若数据可得），Tail(J) 那一套
  "逐点衰减+积分比较"analytic 推导**完全不适用**——尾部控制必须换成对离散数值表的插值/
  外推方案，且这类方案本身很难做到 certificate 级（除非有一个已知解析尾部渐近，比如
  RPA 高频/大动量渐近）。
- **离散频率求和（Matsubara）→ 静态核**：若用户任务里说的是 static susceptibility
  （ω=0 极限，如 Ghosh 数据集的 chi_d 一栏），那么频率求和已经在生成数据时做完，我们只需处理
  动量核；但如果要保留 ω 依赖（做动力学配对核），还需要额外一层 Matsubara 频率求和的
  certificate 化，目前完全没有。

**结论**：区间线性代数外壳（惯性计数、Davis-Kahan）**可复用**；**尾部解析控制层必须整体重写**，
且重写的前提是要有一个具体的、公开的 χ_d(q) 数值表或其渐近公式——目前没有这样的表。

---

## 4. 试运行

**未执行**。原因：现有唯一可下载的真实数据集（Ghosh 2026）不含动量分辨的 χ_d(q)，只有
12 个 (L,β,n) 点上的标量 Sd/Pd/chi_d；且实空间关联文件（本应能重建 λ(k−k′) 核的原始材料）
存在 §2.3 所述的数据完整性问题，不能作为"提取 χ_d(q)，构造 N=50 截断 BS 矩阵"的可信输入。
在没有真实 q-分辨表的情况下，任何"试算 λ₁, λ₂"都只能是自造合成核的空转，无法反映 Hubbard
物理，不产出有意义的 numerical-evidence，因此本轮未做。

若要做有意义的试运行，需要以下任一项：
（a）联系 2603.18205 / 2604.01737 作者索要未公开的动量分辨中间数据；
（b）退而使用 LeBlanc et al. (2023, PRB 107,195138) 的微扰级数展开重建 χ_pair(q)（非 QMC 原始
数据，是对角化 Matsubara 积分给出的解析/半解析结果，可能有闭式或渐近可以支撑 certificate 化）；
（c）等待 Ghosh 数据集作者修正/澄清实空间文件（若真的是占位符，需要作者重新上传）。

---

## 5. Phase 1 数据就绪度（最终判定）

| 判据 | 结果 |
|---|---|
| ≥3 组可用（真实、含动量分辨信息）数据 | **否**：0 组含 χ_d(q) 表；找到 1 组标量数据（Ghosh），且其唯一可能含动量信息的文件（Pd_map）数据完整性存疑 |
| 2026 normalizing-flow 论文数据公开 | **否**：2603.18205、2601.18273 均无数据/代码仓库；论文本身尚未报告大尺寸掺杂结果 |
| 我们的 certificate 工具链是否可直接复用 | **部分**：线性代数外壳（惯性计数、Davis-Kahan）可复用；尾部解析控制需整体重写，且重写依赖具体核形式，目前无核可依 |

**判定：Blocked（核心目标）**。不是"数据需要清洗"的 Partial 状态，而是"动量/频率分辨的
χ_d(q,ω) 数据目前在公开渠道上不存在"——2026 年的符号问题突破论文尚处于方法学验证阶段
（小系统、benchmark 为主），还没有产出可供 Track B 使用的大尺寸物理结果表。

**建议动作**（回报 main 会话）：
1. 建议联系 2603.18205（Schuh, Funcke, Kreit, Luu, Singh, University of Bonn/HISKP）与
   2604.01737（Wang, Lin）作者，询问是否有未发表的动量分辨中间数据，或后续论文（2604.01737
   §IX 结论明确提到"a detailed account of BCS gap solution...will be reported separately"——
   有一篇后续论文在路上，可能包含 χ_d(q) 或等价的配对顶点数据）。
2. Ghosh (2026) Zenodo 数据集的实空间文件问题建议直接联系作者核实（README 未说明 r 范围
   与 L 无关的原因），不建议在问题澄清前用于 certificate 管线。
3. 若急需一个"有闭式/可解析尾部"的动量核做 Track B 概念验证（不追求物理真实性，只验证
   2D BS-certificate 管线本身跑得通），建议下一轮直接用 RPA/微扰级数给出的 χ_d(q) 解析
   近似（如 LeBlanc et al. 2023 的表达式）作为占位核，先把"2D 网格 + 惯性计数复用"这条
   工程路径跑通，等真实 QMC 数据出现后替换核即可。

---

## 6. 交付清单

- 本笔记：`theory/notes/hubbard_qmc_data_interface_2026-09-16.md`（新）。
- 下载/解压的论文与数据集：`theory/notes/_scratch_ir/hubbard_qmc_papers/`
  （PDF ×5，Zenodo 数据集 zip + 解压目录，均未提交，按仓库惯例探针数据不入库）。
- 未修改任何已有证书脚本或基线文件。
- 分支 `p4-conj2`，**未提交**。
