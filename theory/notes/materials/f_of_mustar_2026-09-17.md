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

~~亏损 ∝ 1/√λ 的渐近行为与 N_c ∝ r/√λ 一致，支持「N_c→0 主导强耦合端」的解释。~~

> ⛔ **订正（2026-09-18，主会话从零重算，脚本 `_scratch_ir/_chk_fill_anchor_indep.py`，未提交）**：
> 上面这句**把取整阶梯读成了幂律**。第二路径（稠密 `eigvalsh` 直接对 M，不走久期根）
> 逐位复现了 15.5 / 9.4 / 0.0，但同时给出 **`N_c = ⌊rϖ⌋ ∈ {2, 1, 0}`**（ϖ_c = 0.28534 / 0.16087 / 0.08740）。
> 三个点坐落在三个相邻整数上，**建立不了任何标度律**；λ=100 的「亏损 0.0%」是
> `N_c=0`（窗口里一个 Matsubara 指标都没有，μ* 被取整函数整个关掉）的 **floor 赝物**，
> 就是本文 §178 待办选项 (c) 的那个平凡情形。
> 「N_c→0 主导强耦合端」这个**解释**仍成立；`∝1/√λ` 这个**标度**不成立。
> §178 的「team-lead 承诺从零重算」已兑现并关闭。

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

**原文陈述 [L1, harvest.aps.org PDF]**（AD, PRB 12, 905 (1975), p.908 右栏）：
> "As λ becomes very large, it is found that T_c → (λ⟨ω²⟩)^{1/2} f(μ*(N)), (22), where **f is a function of μ\* only**."
eq.(25), p.909：`T_c ≃ 0.182(λ⟨ω²⟩)^{1/2}` (λ≳10, **μ\*=0**, N+1=64) — the only numeric value given; no other μ* value, no formula for f anywhere in the paper.

**Carbotte RMP 1990 标度定理 [L1]**（RMP 62, 1027, p.1041, eqs.(2.20)–(2.23)）：scaling of α²F(ω)→α²F(bω) leaves T_c/ω_ln invariant — but Carbotte states this **"applies exactly when μ\*=0** and is approximately true when μ\* is not an important player." ⟹ the homogeneity underlying Theorem E is strict only at μ*=0.

**arXiv:2506.19326v2** (Sadovskii et al., strong-coupling ME limit) [WebFetch]: reproduces T_c=0.182√(λ⟨Ω²⟩) but explicitly states "shall not consider the role of direct Coulomb repulsion" — again μ*=0 only.

**Esterlis–Chubukov γ-model** (arXiv:2512.20009, Elezaby & Abanov) [WebSearch]: zero mentions of Coulomb or asymptotic-constant μ* dependence — orthogonal problem (variational truncation-dimension ladder), no collision.

**检索结论**：no one has given C_∞(μ*) or f(λ→∞, μ*) in closed form. AD names the symbol and asserts μ*-only dependence but supplies neither values (beyond μ*=0) nor a formula; Carbotte 1990 / Marsiglio / Sadovskii 2024 all restrict the strong-coupling limit to μ*=0. This is virgin territory, not an overlooked citation.

**Negative search trail** (all checked, no relevant hits — kept for audit): Marsiglio's own Holstein-model papers (arXiv 2010.02801, 1807.04907, 1407.2825, 1105.2833v3); Carbotte Coulomb-pseudopotential literature (PRB 42,920 Josephson; Physica B 165,1069 isotope effect; arXiv:2207.05238v4 "pseudopotential framework is fundamentally biased and uncontrolled" — critiques μ* itself but gives no C_∞(μ*)); AD secondary citations (Xie et al. 2022 npj Comput. Mater. 8,14 — AD-fit accuracy only; arXiv:2104.00244v1 Sadovskii antiadiabatic limit gives a **different** T_c∝(λDΩ₀²)^{1/3} scaling replacing, not correcting, AD's √λ law). Owen & Scalapino, Physica 55, 691 (1971) — Carbotte cites as asymptotic-formula co-source with AD 1975, but full text unreachable (Elsevier paywall, no arXiv); title ("s-state instability") makes an extremality statement unlikely (per lit_einstein memo §10.2, second-hand only, not citable directly).

---

## 2. 数值：fill_factor_probe.py 的关键发现

脚本：`theory/notes/_scratch_ir/fill_factor_probe.py`（未提交）。方法：Einstein 谱、久期根（Theorem CE′ 装置）、N_c = ⌊rϖ⌋ 阶跃截断，N=400, r=10。

**填充因子 f := Tc/(C_∞√(λ⟨ω²⟩)) 与库仑亏损**（worst/best rows; full 11-row table in script output）：

```
   lam   varpi_c    Tc/omE  f(mu*=0) |   Tc(.13)       f   defic |   Tc(.10)       f   defic
   0.5    4.0365   0.03943    0.3052 |   0.01184  0.0916   70.0% |   0.01574  0.1218   60.1%   [worst deficit]
  10.0    0.2853   0.55778    0.9653 |   0.47137  0.8157   15.5% |   0.48786  0.8443   12.5%
  30.0    0.1609   0.98934    0.9885 |   0.89626  0.8955    9.4% |   0.91508  0.9143    7.5%
 100.0    0.0874   1.82098    0.9966 |   1.82098  0.9966    0.0% |   1.82098  0.9966    0.0%   [key row]
```

**关键行**：λ=100, μ*=0.13/0.10: 亏损 = **0.0%**，f(μ*) = f(0) = 0.9966. λ=30: 9.4% still significant. λ=10: 15.5%, approaching H3S's physical-regime ~30%.

**N_c 阶跃截断假象**：N_c = ⌊rϖ_c⌋ → {2, 1, 0} at λ={10, 30, 100}. At λ=100, N_c=0: the Coulomb term `−2μ̃ u_c u_c^T` (u_c = Σ_{n<N_c} u_n) literally vanishes from the operator — not because μ̃=μ*/λ is small, but because the summation range is empty. **This is a discretization artifact, not a physical conclusion.** ϖ ∝ 1/√λ ⟹ N_c ∝ r/√λ, while μ̃ ∝ 1/λ — N_c hits zero faster.

**渐近/弱耦合/N-收敛自检**（脚本输出，全部通过）：
- ϖ_c → 0.87088/√λ（锚点律），λ=1000 相对差 0.05%；f(μ*=0)→1 as expected.
- Tc/ω_E → 0.2530·e^{−1/λ}(1+0.3584λ)：相对差 λ=0.3 处 0.7%，λ=0.5 处 2.3%。
- N=300/400/600 给 f(λ=2.0, μ*=0/0.13) 六位一致（0.819388/0.632204）— N=400 已收敛。

---

## 3. 次项 1/λ 修正的提取（进行中，未完成）

**动机**：既然 C_∞ 与 μ* 无关，有意义的对象是 f(λ,μ*) ≈ 1 − a(μ*)/λ^p，其中 a(μ*) 携带 μ* 依赖。H3S（λ≈1.7–2.3）需要 a、p 才能用。

脚本 `f_mustar_strong_coupling.py`（后台 task b92zy1ghz, b0dpcarv1，**尚未完成**，可能卡在大 λ 的 eigh；预计 10–20 分钟）：固定 μ*∈{0.10,0.13}，扫 λ∈{10,...,100}，拟合 p=1 / p=2 / mixed 三模型，比较 a(0.10) vs a(0.13)。

**手算初步观察**：λ×(1−f) 给 {1.843, 3.135, 0.34} at λ={10,30,100} — 不是常数，p≠1；λ²×(1−f) 给 {18.4, 94.1, 34.0} — 也不是常数，且 λ=100 异常小（N_c=0 假象把真实 1/λ² 项抹掉了）。**推论**：阶跃截断模型在 λ→∞ 不适合提取次项系数，需换 smooth (Fermi–Dirac) cutoff。

**下一步（未完成）**：新脚本 `f_mustar_smooth_cutoff.py` 实现 `f_FD(n;N_c,Δ)=1/(1+exp((n−N_c)/Δ))`（N_c=rϖ 不取整），重跑 λ∈{10,20,30,50,100}×μ*∈{0.10,0.13}，期望 λ=100 库仑亏损非零（连续过渡取代跳变），提取 a(μ*), p。未写脚本，估时 1 小时写 + 20 分钟跑。

---

## 4. 对 H3S 天花板改进的影响（诚实评估）

**原计划**：若 C_∞(0.13) < C_∞(0)，H3S 天花板可从 415 K 压到 ~300 K，且仍是 analytic/certificate。

**当前裁决：否决**——C_∞ 不依赖 μ*（λ→∞ 极限，§2）。

**替代路线（存疑）**：`Tc < C_∞√(λ⟨ω²⟩)×[1−a(μ*)/λ^p]`，但 H3S 的 λ≈1.7–2.3 不在强耦合拟合区（λ≥10），展开可能不收敛/被高阶项污染。若走通，需：(1) 完成 smooth cutoff（§3）；(2) 拟合 a(μ*),p 并量化残差；(3) λ∈[1.5,3] 密集采样验证外推残差≤5%。**三步全未完成**（步骤1未做，步骤2跑中但被阶跃截断污染，步骤3未启动）——**不能在本轮交付中用于改进 H3S 界**。

---

## 5. 核验方式

**fill_factor_probe.py 自洽性已核验**（渐近端/弱耦合端/N 收敛三项自检，见 §2）。**独立路径缺失**：fill_factor_probe.py 是唯一数据源，未用两套独立方案互查（CLAUDE.md 严格性纪律要求承重量两套方案）。若要把 λ=100 的 0.0% 亏损从「观察」升格为「主张」，需要：(a) 另一套求解器路径，或 (b) mpmath 高精度复现，或 (c) 解析证明 N_c=0⟹库仑项消失（平凡但要写清）。

**未核验**：1/λ 修正拟合（§3）——数据不存在，即使完成也只是单一 scipy curve_fit，非两套独立方案。smooth cutoff（§3）——未实现。

---

## 6. 负结果（如实记录）

1. **「C_∞(μ*) 存在且严格小于 C_∞(0)」— 否决**（§2）。数值：λ=100, μ*=0.13 给 f=0.9966=f(0)，亏损 0.0%。原因：N_c→0 比 μ̃→0 更快（离散化假象，非物理）。
2. **「用 1/λ 修正改进 H3S 天花板」— 存疑，不能在本轮交付**（§4）。H3S 的 λ 不在拟合区、smooth cutoff 未完成、外推残差未量化。
3. **「已有人给出 C_∞(μ*) 或强耦合 μ* 修正」— NOT FOUND**（§1）。检索覆盖 AD 1975 / Carbotte RMP 1990 / Marsiglio 综述 / Sadovskii 2024 / Esterlis–Chubukov / pseudopotential 起源论文，全部要么只做 μ*=0，要么不涉渐近常数。

---

## 7. 给 team-lead 的回报

**已确证 [analytic/L1]**
1. AD 1975 eq.(22) 形式上含 f(μ*) 但无显式给出（§1）。
2. 文献未找到 C_∞(μ*) 或强耦合 μ* 修正（§1）。
3. Carbotte RMP 1990 eqs.(2.20)–(2.23) 标度定理可直接引用（定理 E 的齐次性已发表，不必重证；明确只对 μ*=0 严格，支持定理 E 假设合理）。

**仅证据级 [numerical-evidence]**
4. lim_{λ→∞} Tc/√(λ⟨ω²⟩) 与 μ* 无关（阶跃截断模型下）；原因是 N_c→0 比 μ̃→0 更快，离散化假象非物理结论（§2）。
5. 库仑亏损定量行为：λ=10: 15.5%/12.5%；λ=30: 9.4%/7.5%；λ=100: 0.0%/0.0%（μ*=0.13/0.10）；∝1/√λ，与 N_c∝1/√λ 一致（§2）。

**失败或否决**
6. 用 C_∞(μ*) 改进 H3S 天花板：否决（§4）。
7. 用 1/λ 修正改进 H3S 界：存疑，不能在本轮交付（§4）。

**未完成（诚实边界）**
8. 1/λ 修正拟合参数 a(μ*), p — 脚本后台运行中，尚无输出（§3）。
9. smooth cutoff 探针 — 未实现，方案已定（§3）。

**承重项需从零重算的清单**：fill_factor_probe.py 的 λ=100 行（若要作为主张），需独立求解器路径 / mpmath 复现 / 解析证明三选一。**我会从零重算**（team-lead 承诺）。

---

## 附录：脚本清单与复现路径

1. `theory/notes/_scratch_ir/fill_factor_probe.py`（未提交，已跑完，~7分钟, M=900, N=400）— 主力数据源。复现：`cd theory/notes/_scratch_ir && python fill_factor_probe.py`
2. `theory/notes/_scratch_ir/f_mustar_strong_coupling.py`（未提交，跑中，task b92zy1ghz/b0dpcarv1，预计10–20分钟）— 1/λ 拟合三模型。
3. `theory/notes/_scratch_ir/f_mustar_smooth_cutoff.py`（不存在，未写，估计1小时写+20分钟跑）— Fermi–Dirac smooth cutoff。

共享依赖：`concavity_exact_resolvent.build(N,w2)`（仓库已有模块装配 KAY 算子）、`scipy.linalg.eigh`、`scipy.optimize.brentq`（久期根，Theorem CE′ 装置）。无外部数据文件。

**与亏损预算的联系**（`deficit_budget_h3s_2026-09-17.md` §1，三项分解 `Tc(真实)=Tc(天花板)×(1−δ_FIN)(1−δ_BRD)(1−δ_COU)`，δ_COU=27–30% 最大）：本文问「δ_COU 能否从天花板里扣掉（变严格界而非事后亏损）」——λ→∞ 不行（C_∞ 不依赖 μ*）；H3S 的 λ=1.7 可能行（via 1/λ 修正）但存疑（§4 三步未完成）。**保守结论**：δ_COU 仍是事后亏损（numerical-evidence），非天花板的一部分。

**与定理 E / KAY Conjecture 2 的联系**：定理 E（μ*=0, analytic）不受本文发现影响，Carbotte 标度定理支持其假设合理。μ*>0 推广不是「加参数」而是障碍性的（谱隙塌缩，见 `mustar-obstruction-is-gap-collapse-2026-09-17.md`）。AD 1975 未证明 KAY Conjecture 2 的上界（只有形式，无常数）；Moussa–Cohen 2006 的 0.183 上界继承 AD 缺口；⟹ KAY Conjecture 2（μ*=0）是定理 E 最强 framing，不与 AD 撞车。

---

**交付状态**：文献部分完整（阴性结果如实记录）；数值部分部分完成（fill_factor_probe.py 已跑完，1/λ 拟合进行中，smooth cutoff 未写）；主结论已确立（C_∞ 不依赖 μ*）但次项修正存疑；H3S 应用否决（C_∞ 路线）或存疑（1/λ 路线）。
