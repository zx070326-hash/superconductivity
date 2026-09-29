# 谱形因子层级：哪一个「新因子」能压紧天花板——三约束天花板 B(λ, ω̄₁, ω̄₂; μ\*)

> 主会话 2026-09-25。起因（用户）：「过去这些计算器是怎么做出来的？加入新的因子或公式能不能进一步提升这些算子？」
> 分支 `ceiling-shape-factor-hierarchy-2026-09-25`，**UNMERGED**。
> 规范陈述仍以 [`multiband_ceiling_CANONICAL_2026-09-18.md`](multiband_ceiling_CANONICAL_2026-09-18.md) 为准；本篇只**加**一层单带谱形约束，不改头条。

---

## 0. 结论（逐条标级）

| # | 结论 | 严格性 |
|---|---|---|
| F1 | **ω_log（AD75 的 f₂ 因子 ρ_log := ω_log/ω̄₂）不能压紧任何上界**：对任意 ρ_log，在 Einstein 上加 λ 份额 ε、频率 ω₀=exp(ln ρ_log/ε) 的软端原子即可把 ω_log 拉到指定值，而它对 G(k≥1) 的贡献 ~λεω₀²/k² → 0 ⟹ sup 恒等于 Einstein 值 | **analytic**（逃逸构造）；GATE L 数值核（ε=0.2→0.005 时 Tc/Tc_E = 0.9395→0.99830，ρ_log→0.80、M2 恒 1） |
| F2a | **ρ₁ := ω̄₁/ω̄₂ = ⟨ω⟩/√⟨ω²⟩（等价于 Leavens 面积 A = λω̄₁/2）真能压紧**：软端逃逸对 ∫dν 贡献 λεω₀/2→0，硬端逃逸 ≤ M2/(2Ω)→0 ⟹ A 在 𝒞̄ 上弱\*连续 ⟹ A<A_E 真把 Einstein 排除出松弛集 | **analytic** |
| F2b | 三约束松弛集 𝒞̄_ρ = {ν≥0 : λ(ν)≤λ, A(ν)≤λω̄₁/2, M2(ν)≤λω̄₂²} 随三个输入嵌套 ⟹ B 对每个输入单调不减 ⟹ **只需三个量的上界即声音** | **analytic** |
| F2c | 𝒞̄_ρ ⊂ (Leavens 面积集) ∩ (λ-分辨 Einstein 集) ⟹ **B ≤ min(Leavens, λ-分辨 Einstein)** | **analytic**（集合包含）；但 Leavens 常数 0.2309 本身只是数值局部极大（见 [`leavens_c_rational_exact_per_mustar_2026-09-24.md`](leavens_c_rational_exact_per_mustar_2026-09-24.md)、[`same_information_class_comparison_vs_leavens_2026-09-24.md`](same_information_class_comparison_vs_leavens_2026-09-24.md)），故「数值上 B ≤ Leavens 列」一栏是 numerical-evidence：0 违反 |
| F2d | **sup 由单原子取到**，原子落在极值曲线 a(w) = min(1, ρ₁/w, 1/w²) 上（λ 份额 a、频率 w，单位 ω̄₂≡1）；μ\*=0 闭式 `t₁(λ,ρ₁) = max_a min(ρ₁/a, 1/√a)·τ(aλ)`，τ(x)=t\*(x,1)，C = t/(2π√λ) | 曲线代数 **analytic**；「单原子即 sup」**numerical-evidence**（μ\*=0：55 271 个权重≥0.02 的非退化 2/3 原子随机成员 + 240 次 NM，max t/t₁ = 1+1e−12；μ\*=0.10：8 格 max 0.99924；非原子全谱 12/12 满足 Tc_true ≤ Tc_B） |
| F3 | 真实隧穿 α²F 42 行（AD75 Table I 29 行 + Kihlstrom 1984 Nb₃Ge 13 行）上：对仓库头条 C_∞ 紧化**中位 1.61×**（1.26–2.28×）；对 λ-分辨 Einstein 中位 1.105×；**对已有 min(Einstein, Leavens) 中位只 1.010×、最大 1.100×（Sn）** | numerical-evidence |
| F4 | 作计算器（Tc ≈ c·Tc_B）：LOO 误差 4.4–5.5%（随 μ\* 档），与 LOO 重标定的 AD75（4.0–5.4%）**统计上持平**；唯一的区别是它**同时是定理支撑的上界**（29/29 行 r_B ≤ 0.941） | numerical-evidence；**不得声称「比 AD 准」** |
| F5 | H₃S/D₃S 六条 EPW α²F：真实谱达到同 (λ, ω̄₁, ω̄₂) 下最优单原子谱的 **90–97%** ⟹ 给定三矩后「谱形工程」余量只有 3–10%；大杠杆仍是 λ 与 ω̄₂（η/M） | numerical-evidence |

**一句话**：过去的计算器（McMillan、AD75）是「在一族谱上拟合 Eliashberg 数值解」；天花板是「在所有正谱上取 sup」。在这个框架里，**能不能加一个因子，取决于该因子对应的泛函在逃逸序列下是否弱\*连续**——ω_log 不是（一分钱不买，已证），⟨ω⟩ 是（真压紧，但在已有 min(Leavens, Einstein) 之上只多 1–10%）。

---

## 1. 历史：旧计算器是怎么造的（背景，非本篇结论）

| 公式 | 构造方式 | 信息类 | 方向 |
|---|---|---|---|
| McMillan 1968 | Nb 声子谱形状上解 Eliashberg，拟合 Tc(λ, μ\*, Θ_D) | 单一谱形 | 无（拟合） |
| Allen–Dynes 1975 | 217 个 Eliashberg 数值解拟合 f₁(λ)、f₂(ω̄₂/ω_log) 两个修正因子（**f₂ 无开方**，见 [`f2_defect_repair_2026-09-25.md`](f2_defect_repair_2026-09-25.md)） | (λ, ω_log, ω̄₂, μ\*) | 无（拟合，λ→∞ 渐近来自 μ\*=0 解析） |
| LAD 1979 | μ\*=0 强耦合渐近 kTc → 0.1827√(λ⟨ω²⟩) | (λ⟨ω²⟩) | 渐近，非全 λ 上界 |
| Leavens 1975 | 面积 A 固定时 sup Tc，数值局部极大 kTc ≤ 0.2309A | (A) | 上界（数值） |
| Moussa–Cohen 2006 | (λ, λ⟨ω²⟩, ω_max) 下的 Tc 极值问题 | 多矩 | 上界/下界（数值） |
| 本仓库 | 认证 C(λ,μ\*) 曲线 + η/M 一行式 | (λ, λ⟨ω²⟩, μ\*) | 上界（certificate） |

「多用一个矩」这个**做法**属于 Moussa–Cohen 2006 与 Leavens（见 §7）。

---

## 2. 为什么 ω_log 不行而 ⟨ω⟩ 行（F1/F2a 的论证）

线性化 Eliashberg 算子只经 G(k) = ∫ g_k dν，g_k(ω) = 2ω/(ω²+(2πTk)²)（k≥1）依赖 ν；g_k ∈ C₀((0,∞))。
新约束 ∫h dν ≤ c 能排除 Einstein，**当且仅当**不存在一个对 G 贡献 → 0 的逃逸序列能单独满足它。

- **ω_log**：ln ω_log = (2/λ)∫ ln ω · dν/ω。软端原子（份额 ε，频率 ω₀）贡献 ε ln ω₀，取 ω₀ = exp(ln ρ_log/ε) 恰补上所需的对数偏移；
  它对 λ 的份额 ε→0，对 M2 贡献 ∝ εω₀²→0，对 G(k≥1) 贡献 ∝ εω₀²/k²→0。⟹ 约束被满足，Tc → Tc_E。**analytic**。
- **⟨ω⟩（面积）**：软端原子对 ∫dν 贡献 λεω₀/2 → 0；硬端原子（频率 Ω→∞、受 M2 约束）贡献 ≤ M2/(2Ω) → 0。两端都**不能**替约束「买单」
  ⟹ A 弱\*连续 ⟹ A < A_E 真正把 Einstein 挤出集合。**analytic**。

推论（**analytic**）：凡含 ln ω、1/ω 的高负幂等在软端发散/可被软端逃逸调节的因子，**不能**用作上界因子，只能作拟合因子。
GATE F2（`_out_shape_factor_budget_2026-09-25.txt`）40 格给出 C/C_E 随 ρ₁ 单调下降，例如 λ=2：ρ₁=0.80 时 0.8685。

---

## 3. 单原子归约（F2d）

在单位 ω̄₂≡1 下，λ 份额 a、频率 w 的原子满足三约束 ⟺ a≤1、aw≤ρ₁、aw²≤1，最大允许份额

- a = 1（w ≤ ρ₁）
- a = ρ₁/w（ρ₁ ≤ w ≤ 1/ρ₁）
- a = 1/w²（w ≥ 1/ρ₁）

（代数 **analytic**。）Einstein 标度给 t = w·τ(aλ)，故 t₁ = max_a min(ρ₁/a, 1/√a)·τ(aλ)。

**「sup 在单原子取到」目前只是 numerical-evidence**，不是定理。Richter–Rogosinski 只给「≤3 原子」；
单原子性由【S】（面枚举 40 格 max t_budget/t₁ = 1.000000000000）、【R】（55 271 非退化成员 + NM）、【R′】（μ\*=0.10）
和 H₃S 全谱门（§6）支撑。**未做**：区间证书。

μ\*>0 时库仑截断不随 w 标度，脚本中用 ω_c/w，并扫 a(w) 三支全长而非只扫端点；μ\*→0 极限与 t₁ 一致到 1e−8（门）。
根的取法沿用 inf/sup 规则（自上而下几何扫 260 点取最高交点，库仑饱和会造出赝根，判词必须用最高交点）。

---

## 4. 门禁表

| 门 | 内容 | 结果 | 级别 |
|---|---|---|---|
| GATE 0 | Einstein C_E 落在仓库认证括号内（λ=2, 3.2953, 4, 6） | 4/4 IN | numerical-evidence 对 certificate 括号 |
| GATE N | N=64/128/256 收敛 | t\* 差 ≤ 3.6e−7（N=64→256），≤ 2.1e−8（N=128→256）；截断是 λ₁ 下估 ⟹ fail-open，量级远低于任何判决余量 | numerical-evidence |
| GATE L | ω_log 软端逃逸 | Tc/Tc_E → 0.9983（ε=0.005） | 数值核 analytic 构造 |
| GATE F2 | 40 格压紧曲线 | 单调、ρ₁=1 恰回 Einstein | numerical-evidence |
| 【S】 | 闭式 vs 面枚举 | max ratio 1.000000000000 | numerical-evidence |
| 【R】 | 55 271 随机非单原子 + 240 NM（μ\*=0） | max 1+1e−12 | numerical-evidence |
| 【R′】 | μ\*=0.10，ω_c = 1.8/3 ω̄₂，8 格 | max 0.99924 | numerical-evidence |
| 【L】 | new ≤ min(Leavens·ρ₁, Einstein)（48 格） | 违反 0 | numerical-evidence |
| 实数据包含 | r_new ≥ max(r_E, r_L)（42 行） | 违反 0 | numerical-evidence |
| AD 复现 | 裸 AD75 在印刷 μ\* 档 | 3.85% / 12.62%，逐位复现仓库既有门 | numerical-evidence |
| H₃S 全谱门 | Tc_true ≤ Tc_B（6 谱 × μ\*∈{0, 0.13}） | 12/12 | numerical-evidence |

**过程中修掉的两个缺陷**（均已重验修正本身）：
① 首版 NM 频率跑飞产生 NaN/inf ⟹ log-w 夹到 [log 0.01, log 100] 并改为对松弛集的**全面枚举**；
② 公平比较首版误差度量用了 r/mean−1 而非 Tc_pred/Tc−1 ⟹ 裸 AD 给 3.90% 而非仓库门的 3.85%，逐行对比确认转录与实现一致后改正度量，3.85%/12.62% 逐位复现。

---

## 5. 真实数据（F3、F4）

### 5.1 紧度（μ\*=0 天花板，`_out_shape_factor_real_data_2026-09-25.txt`）

r = Tc_exp / 天花板。r_∞ = 仓库头条 C_∞；r_E = λ-分辨 Einstein；r_L = Leavens 面积；r_new = 三约束 B。

- r_new/max(r_E, r_L)：中位 **1.0103**，最大 **1.1000（Sn）**——弱耦合、ρ₁ 小的地方收益最大
- r_new/r_E：中位 1.1052，[1.0459, 1.2639]
- r_new/r_∞：中位 **1.6066**，[1.2601, 2.2756]
- max r_new = 0.8787（Kihlstrom Nb₃Ge #2）；max r_∞ = 0.6102

**读法**：相对头条的 1.6× 大部分是「λ 分辨 + 面积」两份**旧**信息的功劳；ρ₁ 在它们之上只多付 1–10%。
三约束天花板的价值是把两条旧界**统一成一个包封**（它恒 ≤ 两者之 min，F2c），而不是大幅压紧。

### 5.2 计算器（μ\*>0，ω_c = ω_ph，`_out_shape_factor_mustar_predictor_*.txt`、`_out_shape_factor_predictor_fair_*.txt`）

| μ\* 档 | CV r_B | CV r_AD | B·LOO 均/最大 | AD·LOO 均/最大 | AD 裸 均/最大 | max r_B |
|---|---|---|---|---|---|---|
| (a) 印刷 μ\* | 5.69% | 5.29% | 4.81% / 12.27% | 4.01% / 12.80% | 3.85% / 12.62% | 0.857 |
| (b) 0.10 | 6.73% | 6.36% | 5.51% / 16.60% | 5.44% / 14.01% | 5.25% / 13.71% | 0.864 |
| (c) 0.13 | 5.45% | 5.54% | 4.42% / 10.94% | 4.40% / 10.51% | 13.64% / 22.63% | 0.941 |

- ⚠ **μ\* 回流**：AD75 Table I 的印刷 μ\* 是 McMillan–Rowell 反演拟合出的，带 Tc 信息 ⟹ (a) 档对两者都偏乐观；公平比较看 (b)(c) 常数 μ\* 档。
- 裸 AD 在 μ\*=0.13 崩到 13.6%，是因为它的绝对标度对 μ\* 敏感；B 与 AD 的 LOO 重标定版都对 μ\* 档鲁棒。
- 加第四因子 ρ_log（LOO 线性）**不改善**：5.01/5.86/4.33% vs 4.81/5.51/4.42%——与 F1 一致（ω_log 对上界无信息；作预测因子也几乎无增量，corr(r_B, ρ_log)=+0.22–0.34）。
- **不得声称比 AD 准。** 可声称的只有：同等 LOO 精度下，它是**唯一同时是上界**的计算器（r_B ≤ 1 由 F2b+F2d 保证，F2d 是 numerical-evidence）。
- 29 行样本内拟合常数 c=0.7737（印刷 μ\*）；**未做**：语料外预测检验。

---

## 6. H₃S/D₃S 全谱（F5，`_out_shape_factor_h3s_2026-09-25.txt`）

Mishra EPW α²F（`theory/data/external/mishra_h3s_d3s/`，外部数据目录被 `.gitignore` 屏蔽，与仓库既有用法一致），梯形权重作离散测度；ω_c = 10ω̄₂。

| 谱 | μ\* | Tc_true | Tc_B | Tc_E | Leavens | true/B |
|---|---|---|---|---|---|---|
| H3S-200 anharm | 0 | 283.3 | 302.9 | 327.4 | 306.0 | 0.9355 |
| H3S-200 harm | 0 | 328.7 | 340.0 | 350.2 | 370.4 | 0.9668 |
| H3S-200 anharm | 0.13 | 199.0 | 219.9 | 246.8 | 219.6† | 0.9046 |
| D3S-200 anharm | 0.13 | 155.4 | 166.0 | 180.2 | 166.1† | 0.9360 |

（完整 12 行见输出文件。）μ\*=0：B/min(旧) ∈ [0.964, 0.990]。
† μ\*=0.13 的 Leavens 列用 c(0.13, r=10) = 0.165741，其 ω_c 约定是 10ω_E，与此处 10ω̄₂ 只在 Einstein 上一致 ⟹ **H3S-200 anharm 的 B/min = 1.0014 是约定赝物，不是包含关系被违反**；该列只作参照。**未做**：约定对齐重算。

**谱形工程余量**：给定 (λ, ω̄₁, ω̄₂) 后，把 H₃S 的谱「整形」成最优单原子形状最多只能涨 3–10%（μ\*=0.13 时 4.7–9.5%）。
这对材料发现的含义：**形状不是杠杆，λ 和 ω̄₂（即 η/M）才是**——与 [`ultrafast_lamW2_measures_the_ceiling_2026-09-25.md`](ultrafast_lamW2_measures_the_ceiling_2026-09-25.md) 一线「分母可测」的叙事一致。

---

## 7. 优先权与新颖性

- 两矩 (A, M2) 包封：本仓库 09-21/09-22 已做预算（`_chk_two_moment_envelope_budget_2026-09-21.py` 头部已记 MC2006 优先权判断：信息类不同、「做法是他们的」）。
  `_chk_main_two_moment_envelope_2026-09-22.py`（本轮因被引用而 `git add -f`）重跑：Leavens c₀ 重建 0.2308679（印 0.2309，差 0.014%）；在两约束 (A=1, M=m) 集上（λ 不约束），3 个 m 值的 2 原子 sup **恰等于单原子（Einstein）值**，并比 min(Leavens, C_∞√(2m)) 低 5.2–14.3%
  ——即 09-22 已在「无 λ」的两约束集上看到单原子性，本篇把它推到加了 λ 约束的三约束集（numerical-evidence，3 个点，不作承重）。
- Leavens 1975：面积界（单约束）。Moussa–Cohen 2006：(λ, λ⟨ω²⟩, ω_max) 多约束极值。LAD 1979：常数与 η/M 一行式。
- **本篇可主张的只有**：(λ, A, M2) 三约束集、ω_log 的逃逸 no-go、单原子归约与 a(w) 曲线，**在本仓库与已查文献中「未见反例」**。「多用一个矩」的做法不是我们的。

---

## 8. 未做清单

1. 单原子 sup 的 **certificate** 级区间证明（现为 numerical-evidence）。
2. μ\*>0 下三约束天花板的认证（现只有【R′】8 格 + 全谱 12 例）。
3. 计算器的**语料外**预测检验（现 LOO 在同一 29 行内）。
4. Leavens μ\*>0 列的 ω_c 约定对齐重算。
5. 预算：ρ₁ 的上界从哪来（隧穿 α²F / 超快泵浦-探测只给 λ⟨ω²⟩，不给 ⟨ω⟩）——**消费点**：F2b 的单调性要求 ω̄₁ 的**上界**；若只有 ω̄₂ 可测，取 ρ₁=1 退回 λ-分辨 Einstein，三约束不花钱也不买东西。

---

## 9. 脚本与输出（均在 `../_scratch_ir/`，已 `git add -f`）

| 脚本 | 输出 | 内容 |
|---|---|---|
| `_chk_main_shape_factor_budget_2026-09-25.py` | `_out_shape_factor_budget_2026-09-25.{txt,npy}` | GATE 0/N/L/F2，面枚举 |
| `_chk_main_shape_factor_single_atom_and_leavens_2026-09-25.py` | `_out_shape_factor_single_atom_2026-09-25.txt` | 【S】【R】【L】 |
| `_chk_main_shape_factor_real_data_tightness_2026-09-25.py` | `_out_shape_factor_real_data_2026-09-25.txt` | 42 行紧度（数据取自已入库的 `_chk_ad75_table1_frak_f_2026-09-20.py` 与 `_chk_main_kihlstrom1984_nb3ge_tunneling_kappa_2026-09-25.py`） |
| `_chk_main_shape_factor_mustar_predictor_2026-09-25.py` | `_out_shape_factor_mustar_predictor_2026-09-25.txt` | μ\*>0 求解器、【R′】、【P】 |
| `_chk_main_shape_factor_predictor_fair_mustar_2026-09-25.py` | `_out_shape_factor_predictor_fair_2026-09-25.txt` | 三档 μ\* 公平比较 |
| `_chk_main_shape_factor_h3s_spectra_2026-09-25.py` | `_out_shape_factor_h3s_2026-09-25.txt` | H₃S/D₃S 全谱门与谱形余量 |
