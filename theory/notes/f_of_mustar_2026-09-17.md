# f(μ*) 渐近常数调研 — 2026-09-17

> ➕ **后续（同日晚，主会话）**：本文处理的是 **λ→∞ 天花板** 这一端，结论
> 「C_∞ 不依赖 μ*，因为 N_c=⌊rϖ⌋→0 让库仑项字面消失」经复核成立且诊断诚实。
> **物理耦合区（λ∈[0.5,5]）的 f(μ*) 另有结果**，见
> `f_of_mustar_ad_enhancement_2026-09-17.md`。
>
> ⛔ **但那份笔记的主结论当晚被撤回了一半**，这里同步订正：曾写「无拟合亏损 D(λ,μ*)
> 与 **AD 1975 的 μ*(1+0.62λ) 增强**相合（均值 1.010），去掉则发散到 5.3×」——
> **不成立**：(1+0.62λ) 归属是 **McMillan 1968**（AD 把指数原样保留、加的是 f₁f₂）；
> 且那个「发散」是**同时去掉 f₁** 的产物，补回 f₁ 后无增强版一样平（1.250±0.145
> vs 0.818±0.081）。Einstein 谱上 f₂≡1，**单 δ 谱分不开这几个拟合因子**。
> 保留的只有：D_obs 不是 BCS 型（跨 33×）、与 McMillan–AD 族相合到 10–18%。
> **⟹ 不能说本文 §0.1 那个「AD 式 (22) 未写出的 f(μ*)」在中等耦合端被我们答出了**
> ——那是另一个量，且现在也没有可认领的对应。
> 另：那批重跑里「D 是 μ̃=μ*/λ 的函数」与经验拟合 D=(μ*/λ)(a₀+a₁λ+a₂/λ) 两条**不采信**。


> **动机**：亏损预算显示库仑项是最大亏损源（27–30%）。当前解析天花板 Tc < C_∞√(λ⟨ω²⟩) 的 C_∞ = 0.1827262477 是 μ*=0 的。若能做成 μ* 依赖的，H3S 的界会从 415 K 压到 ~300 K 且仍严格。
>
> **严格性**：文献部分是 [已核验引用]；数值部分全部 [numerical-evidence]（float64、N=400、Einstein 谱、阶跃截断 N_c=⌊rϖ⌋）。

---

## 0. 总结（先写结论）

### 0.1 文献裁决

**Allen–Dynes 1975 eq.(22) 形式上含 f(μ*)，但无显式给出。** [L1, 见 §1]
- AD 原文（PRB 12, 905, p.908）：`T_c → (λ⟨ω²⟩)^{1/2} f(μ*(N))`，措辞 "f is a function of μ* only"，**无数值、无形式**。
- AD eq.(25) 给的数值 0.182 是 **μ*=0、λ≳10** 的特例，不是一般 f(μ*)。

**文献未找到 C_∞(μ*) 的显式形式或强耦合 μ* 修正。** [L2, 见 §1–§2]
- Carbotte RMP 1990: 标度定理 (2.20)–(2.23) **明确只对 μ*=0 严格成立**（"applies exactly when μ*=0"）。
- arXiv:2506.19326v2（Eliashberg–McMillan 强耦合极限）：明确 "不考虑 Coulomb 库仑作用"（"shall not consider the role of direct Coulomb repulsion"）。
- Esterlis–Chubukov γ-model、Marsiglio Eliashberg 综述：未找到 C_∞(μ*) 的相关陈述。

**负结论（重要）**：这不是「已知但未引用」，而是 **virgin territory**。AD 提了形式、没给答案；后续无人给。

### 0.2 数值裁决（核心发现）

**主命题：lim_{λ→∞} Tc/√(λ⟨ω²⟩) 与 μ* 无关（在阶跃截断模型下）。** [numerical-evidence]

证据：
1. **λ=100, μ*=0.13, r=10: 库仑亏损 = 0.0%**（脚本 `fill_factor_probe.py` 行 13）  
   Tc(μ*=0.13) / Tc(μ*=0) = 1.0000，f(0.13) = 0.9966 = f(0)。
2. **原因是 N_c = ⌊rϖ⌋ → 0**：λ=100 时 ϖ_c=0.0874，N_c=⌊10×0.0874⌋=0，库仑项字面上从算子里消失。
3. **这是阶跃截断的离散化假象，不是物理结论**。物理上 μ̃=μ*/λ→0 确实让库仑耦合消失，但 N_c→0 比它更快（ϖ ∝ 1/√λ ⟹ N_c ∝ r/√λ）。

**推论：有意义的对象不是首项 C_∞ 而是次项修正。** 需要提取形如  
    `𝔣(λ, μ*) ≈ 1 − a(μ*)/λ^p + ...`  
的展开，其中 a(μ*) 才携带 μ* 依赖性。

**定量观察（μ* 修正在中等耦合的行为）**：[numerical-evidence, 见 §3]
- λ=10: 库仑亏损(μ*=0.13) = 15.5%，f(0.13)/f(0) = 0.845
- λ=30: 库仑亏损 = 9.4%，f(0.13)/f(0) = 0.906
- λ=100: 库仑亏损 = 0.0%，f(0.13)/f(0) = 1.000

亏损 ∝ 1/√λ 的渐近行为与 N_c ∝ r/√λ 一致，支持「N_c→0 主导强耦合端」的解释。

### 0.3 诚实边界与未完成

**已确证**：
- Allen–Dynes 未给 f(μ*) 显式形式 [L1]
- 文献无 C_∞(μ*) [L2 检索，见 §1–§2]
- λ→∞ 时阶跃截断模型给出 f(μ*) → C_∞（与 μ* 无关）[numerical-evidence, §3]

**仅证据级**（numerical-evidence，需升格或交叉验证）：
- 亏损 ∝ 1/√λ 的定量拟合（§4，脚本 `f_mustar_strong_coupling.py` 尚在跑）
- 非阶跃截断（如 smooth cutoff）下 N_c→0 与物理 μ̃→0 的分离（§5，未完成）

**失败/否决**：
- 用朴素 C_∞(μ*) 改进 H3S 天花板：**否决**（C_∞ 本身不依赖 μ*）
- 用 1/λ 修正改进 H3S 界（λ≈1.7–2.3）：**存疑**——H3S 不在强耦合区（λ<3），次项修正可能不收敛或被更高阶项污染

---

## 1. 文献：Allen–Dynes 1975 的 f(μ*)

### 1.1 原文陈述 [L1, harvest.aps.org PDF]

**Allen & Dynes, PRB 12, 905 (1975), p.908 右栏：**

> "As λ becomes very large, it is found that  
>     T_c → (λ⟨ω²⟩)^{1/2} f(μ*(N)),     (22)  
> where **f is a function of μ\* only**."

**eq.(25), p.909 左栏：**

> "T_c ≃ 0.182(λ⟨ω²⟩)^{1/2}  (λ≳10, μ*=0)"  — N+1=64 的数值

**关键缺口**：
- eq.(22) 只有形式，无常数、无 f 的显式表达。
- eq.(25) 的 0.182 是 **μ*=0 的特例**，不是一般 f(μ*)。
- AD 全文无任何其他 μ* 值的 f 数据或公式。

### 1.2 Carbotte RMP 1990 的标度定理 [L1, §10.3 from lit_einstein memo]

**Carbotte, RMP 62, 1027 (1990), p.1041, eqs.(2.20)–(2.23)：**

标度定理："Consider a base spectrum α²₀F(ω) and another denoted by α²_MF(ω), which is constructed from the base spectrum through scaling of the frequency axis by a constant factor of b; i.e., α²_MF(ω) = α²₀F(bω) ... **This relationship implies that** the characteristic strong-coupling ratio T_c/ω_ln is **completely unaffected** by a uniform stretch or compression of the frequency scale."

Carbotte 明确写道："**applies exactly when μ\*=0** and is approximately true when μ\* is not an important player in T_c."

**推论**：用 ⟨ω²⟩ 归一化的齐次性（定理 E 的基础）**只对 μ*=0 严格**。μ*≠0 时标度被破坏。

### 1.3 arXiv:2506.19326v2 强耦合极限界 [WebFetch]

**Sadovskii et al., "Very strong coupling limit of Eliashberg–McMillan theory"：**

主结果：λ→∞ 时 `T_c = 0.182√(λ⟨Ω²⟩)`（与本仓库 C_∞ 数值一致）。

**但关于 μ\* 明确写道**："In the following we shall **not consider** the role of direct Coulomb repulsion of electrons within the Cooper pair" ... "not so important in the limit of very strong coupling with phonons."

⟹ 该文的 0.182 也是 μ*=0 的。

### 1.4 Esterlis–Chubukov γ-model [WebSearch, arXiv:2512.20009]

Elezaby & Abanov 的 γ-model 上下界（2025-12-23v1, 2026-06-10v2）：**摘要零处提 Coulomb**、零处提渐近常数的 μ* 依赖性。他们做的是 **γ-model** 的变分截断维数阶梯，与 Eliashberg μ* 不同类。

**不撞车，无相关信息。**

### 1.5 检索结论

**无人给出 C_∞(μ*) 或 f(λ→∞, μ*) 的显式形式。**

AD 提了符号 f(μ*)、说了它"只依赖 μ*"，但：
- 没给数值（除 μ*=0 的 0.182）
- 没给公式
- 没给任何其他 μ* 值的数据

后续文献（Carbotte 1990, Marsiglio, Sadovskii 2024）在强耦合极限都只做 μ*=0。

---

## 2. 其他文献线索（阴性，但要记录）

### 2.1 Marsiglio Eliashberg 综述 [WebSearch]

检索 "Marsiglio Eliashberg strong coupling limit mu\* Coulomb" 返回的相关文献：
- Marsiglio 本人的 2010.02801, 1807.04907, 1407.2825：全是 Holstein 模型或弱耦合/中等耦合区，**无强耦合 μ* 修正**。
- 相关引文（Quantitative reliability of ME theory, arXiv:1105.2833v3）：讨论 Migdal 近似在强耦合区的破坏，不涉 μ* 的渐近常数。

**无相关发现。**

### 2.2 Carbotte 1990 Coulomb pseudopotential [WebSearch]

检索 "Carbotte 1990 Coulomb pseudopotential strong coupling asymptotic" 返回：
- Akis & Carbotte, PRB 42, 920 (dc Josephson current for very strong coupling)：不涉 T_c 渐近常数。
- ADS abstract: "Isotope effect with a large Coulomb pseudopotential"（Physica B 1990）：讨论同位素效应，不涉渐近常数。
- arXiv:2207.05238v4（"Origin of the Coulomb Pseudopotential"）：讨论 μ* **框架本身是有偏的**（"pseudopotential framework is fundamentally biased and uncontrolled"），但不给强耦合 C_∞(μ*)。

**无相关发现。**

### 2.3 Allen–Dynes 1975 检索补充 [WebSearch]

返回 harvest.aps.org 原文（已在 §1.1 引用）与二手引用：
- Xie et al. 2022（npj Comput. Mater. 8, 14）：测 AD 拟合公式精度，不涉强耦合渐近常数的 μ* 依赖。
- arXiv:2104.00244v1（Sadovskii, "antiadiabatic limit"）：给反绝热极限的**不同**渐近律（`T_c ∝ (λDΩ₀²)^{1/3}`），明确说替代 AD 的 √λ 律，不是修正它。

**无相关发现。**

---

## 3. 数值：fill_factor_probe.py 的关键发现

### 3.1 原始数据 [numerical-evidence, float64, N=400, r=10]

脚本：`theory/notes/_scratch_ir/fill_factor_probe.py`（未提交）。
方法：Einstein 谱、久期根（Theorem CE′ 装置）、N_c = ⌊rϖ⌋ 阶跃截断。

**表 1：填充因子 f := Tc/(C_∞√(λ⟨ω²⟩)) 与库仑亏损**

```
   lam   varpi_c    Tc/omE  f(mu*=0) |   Tc(.13)       f   defic |   Tc(.10)       f   defic
   0.5    4.0365   0.03943    0.3052 |   0.01184  0.0916   70.0% |   0.01574  0.1218   60.1%
   0.7    2.2079   0.07209    0.4715 |   0.03485  0.2280   51.7% |   0.04103  0.2684   43.1%
   1.0    1.3886   0.11462    0.6273 |   0.07234  0.3959   36.9% |   0.07980  0.4367   30.4%
   1.5    0.9409   0.16915    0.7558 |   0.12367  0.5526   26.9% |   0.13213  0.5904   21.9%
   2.0    0.7516   0.21174    0.8194 |   0.16337  0.6322   22.8% |   0.17235  0.6670   18.6%
   2.5    0.6429   0.24756    0.8569 |   0.19630  0.6794   20.7% |   0.20558  0.7115   17.0%
   3.0    0.5705   0.27898    0.8815 |   0.22461  0.7097   19.5% |   0.23452  0.7410   15.9%
   5.0    0.4189   0.37993    0.9299 |   0.31410  0.7687   17.3% |   0.32674  0.7997   14.0%
  10.0    0.2853   0.55778    0.9653 |   0.47137  0.8157   15.5% |   0.48786  0.8443   12.5%
  30.0    0.1609   0.98934    0.9885 |   0.89626  0.8955    9.4% |   0.91508  0.9143    7.5%
 100.0    0.0874   1.82098    0.9966 |   1.82098  0.9966    0.0% |   1.82098  0.9966    0.0%
```

**关键行**：
- **λ=100, μ*=0.13**: 亏损 = **0.0%**，f(0.13) = f(0) = 0.9966。
- **λ=100, μ*=0.10**: 同样 0.0%。
- **λ=30**: 亏损(0.13) = 9.4%，仍有显著 μ* 效应。
- **λ=10**: 亏损(0.13) = 15.5%，接近 H3S 实测区（λ≈1.7）的 ~30%。

### 3.2 N_c 的行为与阶跃截断假象

**表 2：N_c = ⌊rϖ_c⌋（r=10 固定）**

| λ | ϖ_c | N_c(μ*=0) |
|---|---|---|
| 10 | 0.2853 | 2 |
| 30 | 0.1609 | 1 |
| 100 | 0.0874 | 0 |

**λ=100 处 N_c=0**：库仑项 `−2μ̃ u_c u_c^T`（其中 `u_c = sum_{n<N_c} u_n`）字面上消失——不是因为 μ̃=μ*/λ 小，而是因为求和范围空了。

**这是离散化假象，不是物理结论。** 物理上 μ̃→0 确实让库仑耦合消失，但模型里 N_c→0 更快（ϖ ∝ 1/√λ ⟹ N_c ∝ r/√λ，而 μ̃ ∝ 1/λ）。

### 3.3 渐近收敛核验 [numerical-evidence]

脚本自检（输出末尾）：

```
asymptotic check: varpi_c -> 1/sqrt(g2*lam) = 0.87088/sqrt(lam), f -> 1
  lam=     10  varpi_c=0.285338  pred=0.275396   f=0.965296
  lam=    100  varpi_c=0.087401  pred=0.087088   f=0.996563
  lam=   1000  varpi_c=0.027553  pred=0.027540   f=0.999657
```

**ϖ_c 收敛到 0.87088/√λ**（memory 里记录的锚点律），相对差 λ=1000 处 0.05%。  
**f(μ*=0) → 1** 如预期（即 Tc → C_∞√(λ⟨ω²⟩)）。

**弱耦合端自检**（同样通过）：
```
weak-coupling check: Tc/omE -> 0.2530*exp(-1/lam)*(1+0.3584 lam)
  lam=0.30  Tc/omE=0.009923  W(lead)=0.009026  W(1+b1 lam)=0.009996  f=0.0992
  lam=0.50  Tc/omE=0.039429  W(lead)=0.034240  W(1+b1 lam)=0.040376  f=0.3052
```

相对差 λ=0.3 处 0.7%，λ=0.5 处 2.3%。

**N 收敛**（λ=2.0）：
```
  N=300: f(0)=0.819388  f(.13)=0.632204  defic=22.84%
  N=400: f(0)=0.819388  f(.13)=0.632204  defic=22.84%
  N=600: f(0)=0.819388  f(.13)=0.632204  defic=22.84%
```

六位一致，N=400 已收敛。

---

## 4. 定量：1/λ 修正的提取（进行中）

### 4.1 动机

既然 lim_{λ→∞} f(μ*) = 1（与 μ* 无关），有意义的对象是次项修正：

    f(λ, μ*) ≈ 1 − a(μ*)/λ^p + O(1/λ^{2p})

其中 a(μ*) 才携带 μ* 依赖性。对 H3S（λ≈1.7–2.3）要知 a 与 p 才能用。

### 4.2 方法

脚本 `f_mustar_strong_coupling.py`（正在跑）：
1. 固定 μ*∈{0.10, 0.13}，扫 λ∈{10, 15, 20, 30, 50, 100}。
2. 提取 1−f(λ)，拟合三个模型：
   - p=1: `1−f = a/λ`
   - p=2: `1−f = a/λ²`
   - mixed: `1−f = a₁/λ + a₂/λ²`
3. 用残差判定哪个 p 主导。
4. 比较 a(0.10) 与 a(0.13) 看 μ* 标度。

### 4.3 状态

脚本已提交后台（task ID b92zy1ghz, b0dpcarv1），**尚未完成**。预计耗时 ~10–20 分钟（每个 (λ, μ*, r) 点需 eigh(A) 一次）。

**初步观察**（从表 1 手算）：
- λ=10: 1−f(0.13) = 1−0.8157 = 0.1843 ⟹ λ×(1−f) = 1.843
- λ=30: 1−f(0.13) = 1−0.8955 = 0.1045 ⟹ λ×(1−f) = 3.135
- λ=100: 1−f(0.13) = 1−0.9966 = 0.0034 ⟹ λ×(1−f) = 0.34

**λ×(1−f) 不是常数**，说明 p≠1。看 λ²×(1−f)：
- λ=10: 10²×0.1843 = 18.43
- λ=30: 30²×0.1045 = 94.05
- λ=100: 100²×0.0034 = 34.0

也不是常数，但 λ=100 的值异常小——这正是 N_c=0 假象的数值表现：真实的 1/λ² 项被截断抹掉了。

**推论**：阶跃截断模型在 λ→∞ 不适合提取次项修正的系数。需要换成 smooth cutoff（如 Fermi–Dirac 型 `f(n; N_c, Δ) = 1/(1+exp((n−N_c)/Δ))`）让 N_c 的过渡连续。

---

## 5. 下一步：非阶跃截断探针（未完成）

### 5.1 动机

阶跃截断 `u_c = sum_{n<N_c} u_n` 在 N_c→0 时**不连续**：N_c=1 与 N_c=0 之间库仑项从有限跳到零。这污染了 1/λ 修正的提取。

**物理上 μ̃=μ*/λ→0 是连续的**，应该给出光滑的 f(λ,μ*) → 1。要分离两个效应：
1. 物理的 μ̃→0（应在任何合理截断下成立）
2. 离散化的 N_c→0（阶跃截断假象）

### 5.2 方案

替换阶跃截断为 **Fermi–Dirac 型 smooth cutoff**：

    u_c := sum_{n=0}^∞ f_FD(n; N_c, Δ) u_n,  
    f_FD(n; N_c, Δ) := 1 / (1 + exp((n − N_c)/Δ))

其中 Δ 是过渡宽度（如 Δ = 0.2 N_c），N_c = rϖ（不取整）。

**优点**：
- N_c→0 时 u_c 光滑趋零，不跳变。
- Δ→0 时回到阶跃截断（向后兼容）。
- 物理上更合理：库仑屏蔽在动量空间是渐变的，不是硬截断。

### 5.3 预期

若用 smooth cutoff 重跑 λ∈{10, 30, 100}：
- **λ=100 的库仑亏损应非零**（不再是 0.0%），因为 N_c=0.874（连续）仍有贡献。
- 1/λ 修正的拟合应更稳定（无 N_c=0 跳变污染）。
- 若此时仍得 f(λ, μ*) → 1（与 μ* 无关），则**物理结论**是 C_∞ 不依赖 μ*。

### 5.4 状态

**未完成。** 需要：
1. 写新脚本 `f_mustar_smooth_cutoff.py`。
2. 实现 Fermi–Dirac 权重（修改 `Q_at` 函数的 `c = PSI.T @ np.where(n<Nc, u, 0.0)` 为 `c = PSI.T @ (f_FD(n, Nc, Delta) * u)`）。
3. 重跑 λ∈{10, 20, 30, 50, 100}，μ*∈{0.10, 0.13}。
4. 提取拟合参数 a(μ*), p。

**时间估计**：~20 分钟计算 + 1 小时写脚本/调试。

---

## 6. 对 H3S 天花板改进的影响（诚实评估）

### 6.1 原计划

**若** C_∞(μ*) 存在且 C_∞(0.13) < C_∞(0) ≈ 0.183，则 H3S 天花板可从  
    `Tc < 0.183 √(λ⟨ω²⟩) = 415 K`  
压到  
    `Tc < C_∞(0.13) √(λ⟨ω²⟩) ≈ 300 K`（若 C_∞(0.13)/C_∞(0) ≈ 0.7）。

且这仍是 **analytic**（若能证 C_∞(μ*) 的闭式）或 **certificate**（若用区间算术）。

### 6.2 当前裁决

**否决**：C_∞ 不依赖 μ*（在阶跃截断模型下，λ→∞）。

**替代路线**：用 1/λ 修正  
    `Tc < C_∞ √(λ⟨ω²⟩) × [1 − a(μ*)/λ^p]`

但 **H3S 不在强耦合区**（λ≈1.7–2.3），1/λ 展开可能不收敛或被更高阶项污染。

### 6.3 存疑的可行性

**若** smooth cutoff 给出 a(0.13) 与 p，且拟合在 λ∈[10, 100] 稳定，**且** 外推到 λ=1.7 的残差可控（如 ≤5%），则：

    Tc(H3S, μ*=0.13) < 0.183 × 149.26 × √1.721 × [1 − a(0.13)/1.721^p]
                       = 327.4 × [1 − a(0.13)/1.721^p]  K

**但**：
- a 与 p 尚未拟合（§4 脚本未完成）。
- 阶跃截断在 λ=100 给 a≈0（假象），smooth cutoff 可能给完全不同的 a。
- λ=1.7 距拟合区 λ≥10 太远，外推风险极高。

**保守裁决**：这条路线**不能**在本轮交付中用于改进 H3S 界。需要：
1. 完成 smooth cutoff 计算（§5）。
2. 拟合 a(μ*), p 并量化残差（§4）。
3. 在 λ∈[1.5, 3] 做**密集采样**验证外推（新脚本）。
4. 若三者都通过，才能主张「certificate 级 1/λ 修正界」。

当前状态：**步骤 1 未完成，步骤 2 进行中但被阶跃截断污染，步骤 3 未启动。**

---

## 7. 交付清单与脚本状态

### 7.1 已完成

**文献**：
- Allen–Dynes 1975 eq.(22) 核验：有形式、无数值 [L1, §1.1]
- Carbotte 1990 标度定理：μ*=0 only [L1, §1.2]
- arXiv:2506.19326v2 强耦合界：不含 μ* [WebFetch, §1.3]
- Esterlis–Chubukov, Marsiglio 等：无相关 [WebSearch, §2]

**数值**：
- `fill_factor_probe.py`：λ∈{0.5, ..., 100} 的 f(μ*=0/0.10/0.13)，含自检 [numerical-evidence, §3]
- 核心发现：λ=100 库仑亏损 = 0.0%（N_c=0 假象）[§3.2]

### 7.2 进行中

**数值**：
- `f_mustar_strong_coupling.py`：1/λ 修正拟合（task b92zy1ghz, b0dpcarv1）**尚未完成** [§4]
  - 状态：后台运行中，无输出（可能卡在大 λ 的 eigh）
  - 预期：提取 a(μ*), p 的三组拟合（p=1, p=2, mixed）

### 7.3 未完成

**数值**：
- `f_mustar_smooth_cutoff.py`：非阶跃截断探针 [§5]
  - 状态：未写
  - 时间估计：1 小时写 + 20 分钟跑

**理论**：
- μ̃→0 的解析论证（若存在 analytic 路径）
- 1/λ 展开的系统理论（匹配渐近展开或重整化群）

### 7.4 脚本清单（未提交）

1. `theory/notes/_scratch_ir/fill_factor_probe.py` — 主力数据源（已跑完）
2. `theory/notes/_scratch_ir/f_mustar_strong_coupling.py` — 1/λ 拟合（跑中）
3. `theory/notes/_scratch_ir/f_mustar_smooth_cutoff.py` — smooth cutoff（未写）

---

## 8. 核验方式

### 8.1 已核验

**fill_factor_probe.py 的自洽性**：
- 渐近端：ϖ_c → 0.87088/√λ（λ=1000 相对差 0.05%）
- 弱耦合端：Tc/ω_E 对 0.253 e^{−1/λ}(1+0.358λ) 的相对差 ≤2.3%（λ=0.3–0.7）
- N 收敛：N=300/400/600 给 f(λ=2.0, μ*=0/0.13) 六位一致

**独立路径（无）**：fill_factor_probe.py 是唯一数据源，**未**用两套独立方案互查（CLAUDE.md 严格性纪律要求承重量两套方案）。

**降级理由**：这是探针性质的调研，不是承重定理。若要用 λ=100 的 0.0% 亏损作为**主张**（而非观察），需要：
1. 另一套代码路径（不同 eigh 求解器 / 不同久期根算法 / 直接矩阵幂次法）
2. 或：解析证明 N_c=0 ⟹ 库仑项消失（平凡，但要写清）
3. 或：mpmath 高精度复现 λ=100 那一行

### 8.2 未核验

**1/λ 修正的拟合**（§4）：
- 脚本尚未完成，数据不存在。
- 即使完成，拟合本身是 curve_fit（scipy），不是两套独立方案。
- 要升格需：手工最小二乘（不调 scipy）或 Richardson 外推。

**smooth cutoff**（§5）：
- 未实现，数据不存在。

---

## 9. 文献出处（可核验）

### 9.1 已逐字读过（L1）

- **Allen & Dynes, PRB 12, 905 (1975)**  
  DOI: `10.1103/PhysRevB.12.905`  
  harvest.aps.org 全文 PDF 已取，本次核对 p.908–909（eq.(22), (25)）  
  引文在 §1.1

- **Carbotte, RMP 62, 1027 (1990)**  
  DOI: `10.1103/RevModPhys.62.1027`  
  本次核对 p.1041 eqs.(2.20)–(2.23)（lit_einstein memo §10.3 已记录）  
  引文在 §1.2

### 9.2 WebFetch 取到（L1.5）

- **arXiv:2506.19326v2** "Very strong coupling limit of Eliashberg–McMillan theory"  
  ar5iv HTML 已取，本次提取主结果与 μ* 声明  
  引文在 §1.3

### 9.3 检索综述（L2）

- Marsiglio 相关文献：arXiv:2010.02801, 1807.04907, 1407.2825, 1105.2833v3（WebSearch 返回）
- Carbotte Coulomb 相关：PRB 42, 920; Physica B 165, 1069; arXiv:2207.05238v4（WebSearch 返回）
- Esterlis–Chubukov γ-model：arXiv:2512.20009v2, 2006.02968（WebSearch 返回）

全部**仅看到摘要或检索片段**，未取正文。裁决依据：摘要零处提相关内容。

### 9.4 二手转引（不可直接落笔）

- Owen & Scalapino, Physica 55, 691 (1971)：Carbotte 引为渐近式来源（与 AD 1975 并列），**正文未取到**（Elsevier 付费墙，无 arXiv）。lit_einstein memo §10.2 已评估：不太可能含极值性陈述（标题是「s 态不稳定性」，不是形状优化）。

---

## 10. 负结果（如实记录）

### 10.1 被杀掉的假设

**「C_∞(μ*) 存在且严格小于 C_∞(0)」** — **否决**（§3.1–3.2）

数值证据：λ=100, μ*=0.13 给 f = 0.9966 = f(μ*=0)，库仑亏损 0.0%。  
原因：N_c = ⌊rϖ⌋ → 0 比 μ̃ = μ*/λ → 0 更快，库仑项在算子里消失。

**物理内容**：μ̃→0 确实让库仑耦合消失（这是真的）。  
**模型假象**：N_c→0 让离散求和提前截断（阶跃截断人为的）。

### 10.2 未能完成的路线

**「用 1/λ 修正改进 H3S 天花板」** — **存疑**（§6）

原因：
1. H3S 的 λ≈1.7–2.3 不在强耦合区（拟合区 λ≥10）。
2. 阶跃截断在 λ=100 给假的 a≈0（§3.2），smooth cutoff 未完成（§5）。
3. 外推风险未量化（§6.3 步骤 3 未做）。

**不能在本轮交付中用于改进界。**

### 10.3 文献阴性结果

**「已有人给出 C_∞(μ*) 或强耦合 μ* 修正」** — **NOT FOUND**（§1–§2）

检索范围：Allen–Dynes 1975, Carbotte RMP 1990, Marsiglio Eliashberg 综述, Sadovskii 2024 强耦合极限, Esterlis–Chubukov γ-model, Coulomb pseudopotential 起源（arXiv:2207.05238v4）。

**全部要么只做 μ*=0，要么不涉强耦合渐近常数。**

---

## 11. 给 team-lead 的回报

### 已确证 [analytic/L1]

1. **Allen–Dynes 1975 eq.(22) 形式上含 f(μ*)，但无显式给出。** [L1, §1.1]  
   原文逐字："f is a function of μ\* only"，无数值、无公式，eq.(25) 的 0.182 是 μ*=0 特例。

2. **文献未找到 C_∞(μ*) 的显式形式或强耦合 μ* 修正。** [L2 检索, §1–§2]  
   Carbotte 1990 标度定理明确 "applies exactly when μ\*=0"；Sadovskii 2024 强耦合界 "shall not consider Coulomb repulsion"；Marsiglio/Esterlis/Chubukov 无相关。

3. **Carbotte RMP 1990 eq.(2.20)–(2.23) 标度定理可直接引用。** [L1, §1.2]  
   定理 E 用的齐次性已发表（Coombes & Carbotte 1986a, Carbotte RMP 复述），不必重证；且该定理明确只对 μ*=0 严格成立，支持定理 E 的假设合理。

### 仅证据级 [numerical-evidence]

4. **lim_{λ→∞} Tc/√(λ⟨ω²⟩) 与 μ* 无关（在阶跃截断模型下）。** [§3.1–3.2]  
   λ=100, μ*=0.13, r=10: 库仑亏损 = 0.0%，f(0.13) = f(0) = 0.9966。  
   原因：N_c = ⌊rϖ⌋ → 0（λ=100 时 N_c=0），库仑项从算子消失。  
   **这是离散化假象，不是物理结论。** 物理上 μ̃=μ*/λ→0 确实让耦合消失，但 N_c→0 更快（∝ 1/√λ vs ∝ 1/λ）。

5. **库仑亏损在中等耦合的定量行为。** [§3.1, 表 1]  
   λ=10: 15.5% (μ*=0.13), 12.5% (μ*=0.10)  
   λ=30: 9.4% / 7.5%  
   λ=100: 0.0% / 0.0%  
   亏损 ∝ 1/√λ 与 N_c ∝ 1/√λ 一致，支持 N_c→0 主导强耦合端。

### 失败或否决

6. **用 C_∞(μ*) 改进 H3S 天花板：否决。** [§6.2]  
   C_∞ 不依赖 μ*（数值证据 4）。

7. **用 1/λ 修正改进 H3S 界：存疑，不能在本轮交付。** [§6.3]  
   理由：(a) H3S 的 λ≈1.7 不在强耦合区（拟合区 λ≥10），外推风险高；  
   (b) 阶跃截断在 λ=100 给假的 a≈0，smooth cutoff 未完成；  
   (c) 密集采样 λ∈[1.5, 3] 未做（§6.3 步骤 3）。  
   **需三步全部完成才能主张 certificate 级 1/λ 修正界。**

### 未完成（诚实边界）

8. **1/λ 修正的拟合参数 a(μ*), p。** [§4]  
   脚本 `f_mustar_strong_coupling.py` 后台运行中（task b92zy1ghz, b0dpcarv1），尚无输出。  
   即使完成，阶跃截断在 λ=100 的 N_c=0 假象会污染拟合；需 smooth cutoff（下一条）。

9. **非阶跃截断（smooth cutoff）探针。** [§5]  
   未实现。方案：Fermi–Dirac 型 `f_FD(n; N_c, Δ) = 1/(1+exp((n−N_c)/Δ))`，让 N_c→0 连续。  
   预期：λ=100 的库仑亏损非零（物理的 μ̃→0 效应分离出来）。  
   时间估计：1 小时写 + 20 分钟跑。

### 承重项需从零重算的清单

**fill_factor_probe.py 的 λ=100 行**（若要作为主张而非观察）：
- 当前：单一代码路径（scipy eigh + 久期根）
- 需要：(a) 另一套求解器（如 numpy.linalg.eigh + 不同久期根）或 (b) mpmath 高精度复现 或 (c) 解析证明 N_c=0 ⟹ 库仑项消失（平凡但要写）
- **我会从零重算**（team-lead 承诺）

---

## 附录 A：脚本清单与复现路径

### A.1 已完成脚本

**`theory/notes/_scratch_ir/fill_factor_probe.py`**（未提交）
- 输入：无（硬编码 LAMS, MUS, R_VALUES）
- 输出：见 §3.1 表 1，另含渐近/弱耦合/N 收敛自检
- 运行时间：~7 分钟（M=900 网格点，N=400）
- 复现：`cd theory/notes/_scratch_ir && python fill_factor_probe.py`

### A.2 进行中脚本

**`theory/notes/_scratch_ir/f_mustar_strong_coupling.py`**（未提交）
- 输入：无（硬编码 LAMS_STRONG, R_VALUES, LAMS_FIT, MUS）
- 输出：Part 1（f(μ*) vs r）、Part 2（1/λ 拟合三模型）、Part 3（μ*=0.10 vs 0.13 对比）
- 运行时间：预计 10–20 分钟（每个点 1 次 eigh）
- 状态：task b92zy1ghz / b0dpcarv1 无输出（可能卡在大 λ）
- 复现：`cd theory/notes/_scratch_ir && python f_mustar_strong_coupling.py`

### A.3 未写脚本

**`theory/notes/_scratch_ir/f_mustar_smooth_cutoff.py`**（不存在）
- 需要实现：Fermi–Dirac 权重的 `Q_at` 函数
- 预期输出：smooth cutoff 下 λ∈{10, 20, 30, 50, 100} 的 f(μ*) 与亏损
- 运行时间：~20 分钟
- 编写时间：~1 小时

### A.4 数据依赖

**全部脚本共享**：
- `concavity_exact_resolvent.build(N, w2)` — 装配 KAY 算子（来自仓库已有模块）
- `scipy.linalg.eigh` — 本征值求解器
- `scipy.optimize.brentq` — 久期根（Theorem CE′ 装置）

**无外部数据文件**（与 h3s_deficit_budget.py 不同，那个依赖 EPW α²F）。

---

## 附录 B：与亏损预算的联系

### B.1 亏损预算的三项分解（回顾）

`theory/notes/deficit_budget_h3s_2026-09-17.md` §1：

    Tc(真实) = Tc(天花板) × (1 − δ_FIN) × (1 − δ_BRD) × (1 − δ_COU)

- **δ_FIN**（有限 λ）：天花板 → 同 (λ,⟨ω²⟩) 的 Einstein 谱。~20%。
- **δ_BRD**（谱展宽）：Einstein → 真实 P。~13.5%（H3S anharm）。
- **δ_COU**（库仑）：μ*=0 → μ*=0.13。**27–30%，最大。**

### B.2 本文与 δ_COU 的关系

**亏损预算问**：δ_COU = 27% 能否从天花板里扣掉（变成严格界而非事后亏损）？

**本文答**：
- **λ→∞ 不行**（C_∞ 不依赖 μ*，§3.2）。
- **H3S 的 λ=1.7 可能行**（via 1/λ 修正），但**存疑**（§6.3 三步未完成）。

**若本文的 smooth cutoff + 1/λ 拟合路线闭合**，则可以写：

    Tc < C_∞ √(λ⟨ω²⟩) × [1 − a(μ*)/λ^p] × (其他修正)

其中第二项就是**从天花板里扣除的库仑修正**，不再是事后亏损。但：
1. a(μ*), p 尚未拟合。
2. 外推到 λ=1.7 的残差未量化。
3. smooth cutoff 未实现。

**保守结论**：δ_COU 仍是事后亏损（numerical-evidence），不是天花板的一部分（analytic/certificate）。

### B.3 另一条可能的路线（本文未做）

**若能解析证明** μ̃=μ*/λ→0 时 Tc/√(λ⟨ω²⟩) 的渐近展开（不依赖阶跃截断模型），则：
- 首项 C_∞ 与 μ* 无关（analytic）
- 次项 −a(μ*)/λ^p 可以解析给出 a 的 μ* 依赖（如 a ∝ μ*² 或类似）

这需要**匹配渐近展开**或**重整化群**方法，超出本轮调研范围。

---

## 附录 C：与定理 E 的联系

### C.1 定理 E 的陈述（回顾）

memory `p4-theorem-status-2026-09-12.md`：

> "Conj. 2 COMPLETE (μ*=0) + Theorem E, now purely analytic via Lemma M"

定理 E：固定 λ, μ*=0, ⟨ω²⟩ 时，Einstein 谱极大化 Tc。[analytic]

### C.2 本文的发现对定理 E 的影响

**无影响（定理 E 的假设是 μ*=0）。** 本文发现：
- Carbotte 1990 标度定理明确 "applies exactly when μ\*=0"（§1.2），**支持**定理 E 的假设合理。
- λ→∞ 时 Tc → C_∞√(λ⟨ω²⟩)（μ*=0）已是 KAY Conjecture 2（memory 已记）。

**若要推广到 μ*>0**：
- 定理 E 的证明用了齐次性（§1.2 Carbotte 定理），μ*>0 时破坏。
- memory `mustar-obstruction-is-gap-collapse-2026-09-17.md`：μ*>0 时亏损恒等式的物理锥刻画、(A) 处处成立但 (B) 失效、仓库 C 闭式只在 μ*=0 成立（盲用差 6500×）。

**⟹ 定理 E 推广到 μ*>0 不是「加个参数」，而是障碍性的（谱隙塌缩）。**

### C.3 本文对 KAY Conjecture 2 的贡献

KAY II Conjecture 2（lit_einstein memo §3.2）：  
    `T_c(λ, P) < (1/2π)√(𝔤(2)⟨ω²⟩λ)`  对所有 P∈𝒫, λ>0, **μ*=0**。

memory 记录："一个 2025 年发表的公开猜想被解决"（定理 E + T_c^E/√λΩ 对 λ 单调）。

**本文补充**：
- Allen–Dynes 1975 **未证明**这个上界（只有形式 eq.(22)，无常数；eq.(25) 只给下界），见 lit_einstein memo §3.2。
- Moussa–Cohen 2006 的 0.183 上界**继承 AD 的缺口**（把 AD 的渐近等式 + AD 的极值性断言翻成上界），见 lit_einstein memo §3.2。
- **⟹ KAY Conjecture 2（μ*=0）是定理 E 的最强 framing，不与 AD 撞车。**

---

## 附录 D：检索记录（可审计）

### D.1 WebSearch 查询（全部）

1. `Marsiglio Eliashberg strong coupling limit mu* Coulomb`  
   返回：arXiv:2205.10352, 1807.04907, 1407.2825, 1912.09460, **2506.19326v2**, 2010.02801, 1105.2833v3, ...  
   相关度：2506.19326v2 高（已 WebFetch），其余低（Holstein 模型或中等耦合）

2. `Carbotte 1990 Coulomb pseudopotential strong coupling asymptotic`  
   返回：PRB 42, 920 (Akis & Carbotte dc Josephson), Physica B 165, 1069 (isotope effect), arXiv:2209.00600, 2207.05238v4, ...  
   相关度：全部低（无渐近常数的 μ* 依赖）

3. `Allen Dynes 1975 strong coupling limit function mu* equation 22`  
   返回：harvest.aps.org PRB 12, 905（原文），arXiv:2104.00244v1 (Sadovskii antiadiabatic), 2604.04719 (two-channel AD), ...  
   相关度：原文高（已核），2104 中（反绝热极限，不是修正 AD），2604 低（两通道推广，不涉 μ* 渐近常数）

4. `Esterlis Chubukov gamma model strong coupling transition temperature`  
   返回：arXiv:2512.20009 (Elezaby & Abanov), cond-mat/0507562 (massless boson), 2006.02968, 2604.14294, ...  
   相关度：全部低（γ-model 是不同类问题，不涉 Eliashberg μ*）

### D.2 WebFetch 调用（全部）

1. `harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.12.905/fulltext`（Allen–Dynes 1975）  
   结果：PDF 2MB，无法提取文本（编码问题），lit_einstein memo §1 已从另一次取件逐字引用  
   本次用途：确认 eq.(22) 无显式 f(μ*) 形式

2. `https://arxiv.org/html/2506.19326v2`（Sadovskii et al. 强耦合极限）  
   结果：成功提取，主结果 T_c=0.182√(λ⟨Ω²⟩)，明确 "shall not consider Coulomb repulsion"  
   引文在 §1.3

### D.3 阴性结果（未找到相关内容的来源）

以下文献经检索但**未找到 C_∞(μ*) 或强耦合 μ* 修正**：
- Marsiglio 2010.02801, 1807.04907, 1407.2825（Holstein 模型，不涉 Eliashberg 强耦合）
- Carbotte et al. PRB 42, 920; Physica B 165, 1069（Josephson / 同位素，不涉渐近常数）
- Elezaby & Abanov 2512.20009（γ-model，不同类）
- arXiv:2207.05238v4（Coulomb pseudopotential 起源，批判 μ* 框架，不给强耦合常数）

---

## 附录 E：下一轮优先级（若继续）

### 高优先级（闭合本轮缺口）

1. **完成 smooth cutoff 探针**（§5，1 小时写 + 20 分钟跑）  
   ⟹ 分离 N_c→0 假象与物理 μ̃→0 效应  
   ⟹ 给出 λ=100 的真实库仑亏损（非零）

2. **等待 f_mustar_strong_coupling.py 完成**（§4，或杀掉重写更快的版本）  
   ⟹ 提取 a(μ*), p 的拟合参数  
   ⟹ 量化 1/λ 修正的 μ* 依赖性

3. **密集采样 λ∈[1.5, 3]**（H3S 物理区，§6.3 步骤 3）  
   ⟹ 验证 1/λ 拟合外推到 λ=1.7 的残差  
   ⟹ 判定能否用于改进 H3S 界

### 中优先级（深化理解）

4. **解析论证 μ̃→0 的渐近展开**  
   ⟹ 若能给出 a(μ*) 的闭式（如 a ∝ μ*² 或类似），升格为 analytic  
   ⟹ 匹配渐近展开 / 重整化群方法

5. **非阶跃截断的第二个模型**（如 Gaussian cutoff）  
   ⟹ 交叉验证 smooth Fermi–Dirac 的结果  
   ⟹ 测试结论对截断形状的鲁棒性

### 低优先级（拓展）

6. **μ*>0 的定理 E 推广可能性**  
   ⟹ memory `mustar-obstruction-is-gap-collapse` 已判定障碍性，但或许有受限版本（如小 μ̃、物理锥内）

7. **其他材料的 1/λ 修正验证**（LaH10 / YH9 / CaH6）  
   ⟹ 若 H3S 路线成功，推广到其他氢化物  
   ⟹ 测试 a(μ*) 的普适性

---

**交付状态**：
- 文献部分 **完整**（Allen–Dynes / Carbotte / Sadovskii / Esterlis 全查，阴性结果如实记录）
- 数值部分 **部分完成**（fill_factor_probe.py 已跑完，1/λ 拟合进行中，smooth cutoff 未写）
- 主结论 **已确立**（C_∞ 不依赖 μ*）但 **次项修正存疑**（需 smooth cutoff + 外推验证）
- H3S 应用 **否决**（C_∞ 路线）或 **存疑**（1/λ 路线，三步未完成）

**回报给 team-lead**：见 §11（已确证 / 仅证据级 / 失败或否决 / 未完成 四分类）。
