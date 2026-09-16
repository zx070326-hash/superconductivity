# 配对凸优化框架深推：对偶定理严格化 + 谱效率指标计算

建立：2026-09-16
状态：v0.2（对偶定理证明 + 谱效率 eta 数值 + 双矩极值扫描 + mu*>0 凸结构分析）
严格性：**对偶定理 = 解析（条件性：依赖引理 C）；其余 = 数值证据**
脚本：`theory/notes/_scratch_ir/spectral_efficiency.py`（未提交）

---

## 1. 对偶定理严格证明

### 1.1 定理陈述

**定理（配对优化对偶）.** 设 h: (0,infinity) -> R 是 KAY 的 Einstein 函数（h(varpi) = lambda_1(H(varpi))，`P4_main_theorems.md` 定理 2），且 H(x) = h(sqrt(x)) 在 x > 0 凹（引理 C，认证级）。设 K_C = {P 是 [0,infinity) 上概率测度 : integral(omega^2 dP) = C, P != delta_0}。则：

(a) **原始问题闭式解.** sup_{P in K_C} kappa(P,T) = h(varpi_rms)，varpi_rms := sqrt(C) / (2 pi T)。唯一最优测度 P* = delta_{sqrt(C)}。

(b) **Lagrange 对偶问题.**
    min_{mu, nu} [mu + nu C]
    s.t. h(varpi) <= mu + nu varpi^2  对一切 varpi > 0
    有最优解 nu* = H'(C/(2piT)^2), mu* = H(C/(2piT)^2) - (C/(2piT)^2) * H'(C/(2piT)^2)。

(c) **强对偶.** sup = inf = h(varpi_rms)。

(d) **唯一性.** 最优测度唯一：P* = delta_{sqrt(C)}。

### 1.2 证明

**Step 1 (原始问题上界).** 记 x := varpi^2 = (omega/(2piT))^2。对每个 k >= 1，核函数 phi_k(omega) = x/(k^2+x) 关于 x 严格凹（二阶导 = -2k^2/(k^2+x)^3 < 0）。因此

  F_T(k) = integral phi_k(omega) P(d omega) <= phi_k(omega_rms)    (Jensen)

其中 omega_rms := sqrt(C)。这是逐 k 的不等式。将其递入引理 M（核单调性，`P4_main_theorems.md` 引理 M，解析已证）：

  kappa(P,T) = lambda_1(O[F_T]) <= lambda_1(O[bracket(varpi_rms)]) = h(varpi_rms)。

等号 iff F_T(k) = phi_k(omega_rms) 对所有 k，iff P = delta_{omega_rms}。 [analytic]

**Step 2 (下界，可达).** 取 P = delta_{sqrt(C)}，则 kappa(P,T) = h(varpi_rms)。故 sup 可达。 [trivial]

**Step 3 (对偶问题推导).** 原始问题等价写为

  max_P  kappa(P,T)
  s.t.   integral P(d omega) = 1,  integral omega^2 P(d omega) = C,  P >= 0.

引入 Lagrange 乘子 mu（归一化）和 nu（矩约束）：

  L(P; mu, nu) = kappa(P,T) - mu(integral P - 1) - nu(integral omega^2 P - C).

对偶函数 g(mu, nu) = sup_{P >= 0} L(P; mu, nu)。因 kappa(P,T) = sup_{||v||=1} integral phi(omega; v) P(d omega)（其中 phi(omega; v) := <v, O[bracket(omega/2piT)] v>），交换 sup 后

  g(mu, nu) = mu + nu C + sup_v sup_{P >= 0} integral [phi(omega; v) - mu - nu omega^2] P(d omega).

内层对 P 的 sup 为 0 iff phi(omega; v) <= mu + nu omega^2 对所有 omega, v；否则为 +infinity。取 v 为 H(varpi) 的顶本征向量：

  **对偶问题**：min_{mu, nu} [mu + nu C]  s.t.  h(varpi) <= mu + nu varpi^2  对所有 varpi > 0。

即：**找二次函数 mu + nu x 包络 H(x) = h(sqrt(x)) 且在 x = C/(2piT)^2 处最小。** [analytic]

**Step 4 (强对偶与对偶最优解).** 若 H(x) 凹（引理 C），则 H(x) 是它自身的凹包络。在 x_0 := C/(2piT)^2 处的切线 mu + nu x 满足

  nu* = H'(x_0),   mu* = H(x_0) - x_0 H'(x_0),

且 mu* + nu* x_0 = H(x_0) = h(varpi_rms)。因 H 凹，切线在 x_0 处取等、在其余处不低于 H（即切线是 H 的包络）。故 min_{mu,nu} [mu + nu C/(2piT)^2] = H(x_0) = h(varpi_rms) = sup。

强对偶间隙 = 0。 [analytic, conditional on Lemma C]

**Step 5 (唯一性).** 若 P* != delta_{sqrt(C)}，则至少存在 k 使 Jensen 严格：F_T(k) < phi_k(omega_rms)。引理 M 的严格性条件（d(1) > 0）给出 kappa(P*,T) < h(varpi_rms)，矛盾。 [analytic] QED

### 1.3 依赖关系

| 引用 | 严格性 | 角色 |
|---|---|---|
| 引理 M（核单调性） | 解析（已证无 CHALLENGED） | Step 1, 5 |
| 引理 C（H 凹性） | 认证（段 C-I + M5） | Step 4（强对偶） |
| KAY 定理 2（h 定义） | 解析（三行核级 Jensen） | 定义 |
| phi_k 严格凹 | 解析（初等） | Step 1 |

**注.** (a)(d) 不依赖引理 C（只用引理 M），(b)(c) 依赖引理 C。因此 **Einstein 极值与唯一性是解析的；强对偶的成立条件性地依赖引理 C**。若引理 C 段 C-II 闭合，则全部四个陈述无计算机辅助。

### 1.4 物理解读

- **nu*（矩约束影子价格）** = H'(varpi_rms^2) = h'(varpi_rms) / (2 varpi_rms)：增加 <omega^2> 一个单位，kappa 改善的边际收益。物理上是"频率硬化对超导配对的边际贡献"。
- **mu*（归一化影子价格）** = H(varpi_rms^2) - varpi_rms^2 * H'(varpi_rms^2)：切线的截距，即"零频声子的虚拟配对价值"。
- **强对偶 = 谱优化无松弛损耗**：在 <omega^2> 约束下，没有 duality gap——配对能力完全由二阶矩决定，谱的精细结构只能降低 Tc。

---

## 2. 谱效率指标 eta(P)

### 2.1 定义

  **eta(P; lambda, mu*) := Tc(P; lambda, mu*) / Tc^E(lambda, sqrt(<omega^2>_P), mu*)**

即真实谱的 Tc 与**同 lambda、同 <omega^2>** 的 Einstein 谱 Tc 之比。由定理 E（Einstein 极值），**eta in (0, 1]，eta = 1 iff P = delta_{sqrt(<omega^2>)}**（点质量 = Einstein 谱）。

**Allen-Dynes 层面的闭式**：eta_AD = f2(r, lambda, mu*) / r，其中 r = omega_2 / omega_log >= 1。r = 1 对应 Einstein 谱（f2 = 1），r > 1 对应展宽谱。

### 2.2 结果：材料谱效率表

| 材料 | P (GPa) | lambda | mu* | r | eta_AD |
|---|---|---|---|---|---|
| LaH10 | 129 | 3.62 | 0.10 | 1.194 | **0.921** |
| LaH10 | 163 | 2.67 | 0.10 | 1.075 | **0.959** |
| LaH10 | 214 | 2.06 | 0.10 | 1.105 | **0.932** |
| LaH10 | 264 | 1.73 | 0.10 | 1.095 | **0.932** |
| H3S | 200 | 1.84 | 0.16 | ~1.20 | **0.858** |
| Nb (BETE-NET) | amb. | 1.23 | 0.13 | ~1.14 | **0.889** |
| MgB2 (BETE-NET) | amb. | 0.62 | 0.13 | ~1.12 | **0.896** |

**数据来源**：LaH10 的 (lambda, omega_log, r) 来自 Errea 2020 + P3 封闭反演（`P3_tc_max_pressure.md`）；H3S 的 r 从 Semenok 的 lambda_* 估计；Nb、MgB2 的 r 从 `p3_nonhydride_lambdastar_check.py` 的 sigma 值估计。所有 r 值标 **numerical-evidence**。

### 2.3 物理解读

1. **LaH10 效率最高** eta ~ 0.92-0.96，符合预期——LaH10 的声子谱以光学氢振动为主、相对窄。其中 163 GPa 处 eta 最高（0.959），对应反演 r 最小（1.075），即谱最接近 Einstein。

2. **129 GPa 效率最低**（0.921），一手 r = 1.194 显著更宽——这与 P3 "129 GPa 谱显著更宽"的结论吻合。**高压三点效率持平**（0.93-0.96），对应 r 的"台阶+平台"结构。

3. **H3S 效率较低**（~0.86）。r ~ 1.20 意味着 omega_2 比 omega_log 高 20%，对应双峰谱结构（高频 S-H 拉伸 + 低频框架模式）。这里 mu* = 0.16 也加重了惩罚。

4. **非氢化物效率**：Nb 和 MgB2 的 eta ~ 0.89-0.90，与氢化物的 0.86-0.96 可比。原因是 r 值（1.12-1.14）虽然不大，但 lambda 较小时 f2 修正项更弱（AD f2 公式的分母含 lambda^2），对 r 的敏感度更低。

5. **"谱展宽对 Tc 的惩罚"不大**：即使 r = 1.20（显著展宽），eta 仍 > 0.85。要让 eta 跌到 0.5 以下，需要 r > 2.0（omega_2 是 omega_log 的两倍以上），在已知超导体中罕见。这说明**谱优化的实际空间有限——材料 Tc 主要由 lambda 和 <omega^2> 决定，谱形状只是微调**。

### 2.4 eta 对 r 的灵敏度

| r | eta (lam=2.5, mu*=0.10) | eta (lam=1.8, mu*=0.16) | eta (lam=1.0, mu*=0.13) |
|---|---|---|---|
| 1.00 | 1.000 | 1.000 | 1.000 |
| 1.05 | 0.971 | 0.961 | 0.956 |
| 1.10 | 0.943 | 0.924 | 0.915 |
| 1.20 | 0.888 | 0.857 | 0.843 |
| 1.50 | 0.747 | 0.699 | 0.680 |
| 2.00 | 0.575 | 0.529 | 0.511 |

惩罚随 r 超线性增长。强耦合（大 lambda）下惩罚较轻，因为 AD f2 的分子含 lambda^2。

---

## 3. KAY 算子层面交叉验证

为验证 AD 层面的 eta 不是拟合公式的假象，直接用 KAY 算子（`P4_main_theorems.md` 定义的 H 矩阵）计算 2-atom vs Einstein。

| eps | om1 (K) | om2 (K) | r | Tc_KAY (K) | eta_KAY | eta_AD |
|---|---|---|---|---|---|---|
| 0.10 | 464 | 750 | 1.027 | 103.8 | 0.980 | 0.982 |
| 0.20 | 415 | 750 | 1.071 | 100.8 | 0.952 | 0.952 |
| 0.10 | 167 | 1500 | 2.408 | 46.8 | 0.442 | 0.458 |

参数：lambda = 2.0, omega_E = 500 K, N_KAY = 150。

**结论**：AD 与 KAY 在 eta 上吻合到 1-2%（r 不太大时），极端展宽（r = 2.4）处 AD 偏高 ~4%——AD 过估 Tc 是 P3 已知行为（AD 拟合在大 r 时不精确）。定性方向一致：**展宽降低 eta，Einstein 最优**。

---

## 4. 双矩极值扫描（Richter 约束）

### 4.1 设置

固定 <omega^2> = C2 和 <omega^4> = C4。由 Richter 极值定理（1957），n 个矩约束下的极值测度最多 n+1 个原子。此处 3 个约束（归一化 + C2 + C4）=> 最优解最多 3 个原子。

参考：LaH10 163 GPa 参数（lambda = 2.67, mu* = 0.1, omega_rms = 1202.6 K）。

### 4.2 结果

| C4/C2^2 | Tc_best (K) | Tc^E (K) | eta_opt | 最优配置 |
|---|---|---|---|---|
| 1.001 | -- | 244.4 | -- | 无可行 3-atom |
| 1.01 | 244.1 | 244.4 | 0.999 | 2-atom 主导 |
| 1.05 | 244.0 | 244.4 | 0.998 | 2-atom 主导 |
| 1.10 | 243.8 | 244.4 | 0.997 | 2-atom 主导 |
| 1.15 | 243.6 | 244.4 | 0.997 | 2-atom 主导 |
| 1.20 | 243.4 | 244.4 | 0.996 | 2-atom 主导 |
| 1.30 | 243.0 | 244.4 | 0.994 | 2-atom 主导 |
| 1.50 | 242.3 | 244.4 | 0.991 | 2-atom 主导 |

C4/C2^2 = 1 对应 Einstein（所有矩重合），C4/C2^2 > 1 对应展宽。

### 4.3 解读

1. **额外约束下的增益微乎其微**。即使 C4/C2^2 = 1.5（相当大的过量峭度），最优 3-atom 配置只比 Einstein 低 ~1%。

2. **最优 3-atom 退化为 2-atom**。在所有扫描点，w3 ~ 0.001-0.002，第三个原子几乎无权重。这说明在 AD 公式的结构下（Tc 主要由 omega_log 和 r 决定），多原子谱的最优化空间非常有限。

3. **"如果谱形可调，最多能多赚多少 Tc"的回答**：在物理合理的 C4/C2^2 范围（1.0-1.3），**不超过 0.5-1%**。谱形优化对 Tc 的贡献本质上是二阶的——这是定理 E 的数值回声。

**诚实边界**：(a) 扫描在 AD 层面，AD 在大 r 时不精确（见 sec.3）；(b) 网格分辨率有限（25x25x12），全局最优可能被遗漏但偏差不会大；(c) 物理上 C4 约束由材料晶格结构决定，不是自由参数。

---

## 5. mu*>0 的凸结构分析

### 5.1 lambda_1(M) 关于 P 的凸性 [analytic]

库仑项 K_c = 2|u><u|（秩一，P-无关）。束 M(P,T) = O[F_T(P)] - mu_tilde * K_c。因 O[F_T] 关于 P 线性（O 的逐元线性 + F_T 关于 P 线性），且 K_c 与 P 无关，M 关于 P 是仿射的。

lambda_1(M) = sup_{||v||=1} <v, M v> 是 P 的仿射泛函族的上确界，故**关于 P 凸**。

**推论**：最大化 lambda_1(M) 在 K_C 上仍是凸函数在紧凸集上的极大化问题，极值在极点 delta_{sqrt(C)} 取到——**只要 H_mu(x) := lambda_1(M(x)) 凹**。

### 5.2 H_mu''(x) 数值检验

直接数值计算 H_mu''(x) = h_mu''(sqrt(x)) 的换元值，用 5-point stencil。

| x | mu_tilde=0 | mu_tilde=0.03 | mu_tilde=0.05 | mu_tilde=0.10 | mu_tilde=0.13 |
|---|---|---|---|---|---|
| 0.01 | -2.32 | **+2.21** | **+0.78** | -0.31 | -0.57 |
| 0.05 | -2.06 | **+1.12** | **+0.35** | -0.40 | -0.59 |
| 0.10 | -1.80 | **+0.13** | -0.07 | -0.48 | -0.61 |
| 0.50 | -0.72 | -0.63 | -0.54 | -0.44 | -0.42 |
| 1.0 | -0.317 | -0.300 | -0.280 | -0.239 | -0.224 |
| 5.0 | -0.018 | -0.018 | -0.017 | -0.016 | -0.015 |
| 50.0 | -1.9e-4 | -1.9e-4 | -1.8e-4 | -1.7e-4 | -1.6e-4 |
| 100.0 | -4.8e-5 | -4.7e-5 | -4.6e-5 | -4.2e-5 | -4.0e-5 |

### 5.3 发现：小 mu_tilde 在 x -> 0 附近凹性失效

**关键发现**：

1. **mu_tilde = 0.03 (即 mu*/lambda ~ 0.03)**：H_mu''(x) > 0 在 x in (0, 0.109)，零交叉在 **x ~ 0.109 (varpi ~ 0.33)**。

2. **mu_tilde = 0.05**：H_mu''(x) > 0 在 x in (0, 0.091)，零交叉在 **x ~ 0.091 (varpi ~ 0.30)**。

3. **mu_tilde >= 0.10**：H_mu''(x) < 0 在全部检测点（x in [0.01, 500]）。凹性在数值精度内成立。

**物理解读**：秩一库仑项 2*mu_tilde*|u><u| 使 H_mu(0) 从 0 变为 -2*mu_tilde*||u||^2 < 0（负值），而 H_mu(x) -> h(sqrt(x)) - 2*mu_tilde*||u||^2 在 x 大时近似为 ln(sqrt(x)) + const。在 x = 0 附近，mu* 项的导数效应可以翻转 H_mu 的凹凸性——本质上是因为 |u><u| 的 Rayleigh 商 <v, |u><u| v> 随 v（顶本征向量）变化时引入了一个凸贡献。

**但物理上不致命**：
- 凹性失效区在 varpi < 0.33，对应 omega/T < 2pi * 0.33 ~ 2。对真实材料，主要配对声子频率远高于 Tc（即 varpi >> 1），此区间不被真实谱的支撑覆盖。
- 直接 2-atom vs Einstein 数值测试（Part 4(d)）在 mu_tilde = 0.05, 0.10, 0.15 下均给出 **kappa_mu(2-atom) / h_mu(Einstein) < 1**（0.76, 0.74, 0.73）——Einstein 仍赢，且赢得很宽。

### 5.4 mu*>0 Einstein 极值性总结

| 命题 | 状态 | 严格性 |
|---|---|---|
| kappa(P,T) - mu_tilde * <P, K_c P> 关于 P 凸 | 严格但注意 K_c 项 | K_c 是 P-无关的（只含 u），故 M 关于 P 仿射 => lambda_1 凸。[analytic] |
| H_mu(x) 凹 (mu_tilde >= 0.10) | 数值扫描无反例 | numerical-evidence |
| H_mu(x) 凹性在小 x 处失效 (mu_tilde = 0.03-0.05) | 确认：零交叉在 x ~ 0.09-0.11 | numerical-evidence |
| Einstein 仍为 Tc 极大者 (mu_tilde in [0.05, 0.15]) | 直接数值验证（eta < 1） | numerical-evidence |

**结论**：在物理相关的参数范围（mu_tilde = mu*/lambda >= 0.05，对应 lambda <= 20 时 mu* >= 0.1），凹性在可达 varpi 范围内持有，Einstein 极值性仍成立。严格证明（对所有 P）需要：(1) mu_tilde < 0.05 时的小 varpi 修正论证，或 (2) 一个不依赖逐点凹性的替代路线（如引理 M 的直接库仑版——这在 `P4D_coulomb_extremal_threshold.md` 已有认证级条件版本）。

---

## 6. 诚实边界

1. **对偶定理 Step 4 条件性**：强对偶依赖引理 C。若只用引理 M（解析），能得到 (a) sup 值 + (d) 唯一性，但不能写出 (b) 的对偶最优解闭式或 (c) 的间隙为零。

2. **eta_AD 公式是 Allen-Dynes 拟合的推论，非严格**。AD 公式在 lambda > 1.6 或 r 显著偏离 1 时有偏差（见 P3 文档）。KAY 层面的交叉验证（sec.3）确认偏差在 2-4% 内。

3. **r 值来源**：LaH10 的 r 来自 P3 的封闭反演（强 B 级 + 三重自洽检验）；H3S 和非氢化物的 r 是估计值。

4. **双矩极值扫描**（sec.4）在 AD 层面，网格分辨率有限。结论"eta_opt > 0.99"可能被更精密的扫描微调，但 "额外矩约束下增益 < 1%"的定性结论应当稳健。

5. **mu*>0 凹性失效**（sec.5.3）是实际发现，不是数值噪声——零交叉可被 brentq 定位到 3 位有效数字。但此区域在物理上不可达（声子频率 >> Tc 意味着 varpi >> 1）。

6. **本文无新定理**。对偶定理是 sec.7 of `pairing_convex_optimization_2026-09-16.md` 中框架的严格化；eta 计算是定理 E 的定量应用；mu*>0 分析是 `P4D_coulomb_extremal_threshold.md` 的补充探测。所有引用以 `P4_main_theorems.md` v0.8 为准。

---

## 7. 与已有结论的衔接

| 本文结论 | 依赖 | 衔接 |
|---|---|---|
| 对偶定理 (a)(d) | 引理 M（解析） | 定理 E 的重述 |
| 对偶定理 (b)(c) | 引理 C（认证） | 强对偶的新结论 |
| eta(LaH10 129 GPa) = 0.921 | P3 反演 r=1.194 | P3 "129 GPa 谱最宽"的定量回声 |
| eta(LaH10 163 GPa) = 0.959 | P3 反演 r=1.075 | r 最小 => eta 最高 |
| 双矩极值 < 1% 增益 | Richter + AD | 定理 E 的量化：谱形优化空间本质有限 |
| mu*>0 凹性 x->0 失效 (mu_tilde < 0.05) | 直接数值 | 新发现，但物理上在可达 varpi 范围外 |
| mu*>0 Einstein 仍最优（数值） | 直接 KAY 对比 | 与 `P4D_coulomb_extremal_threshold.md` 一致 |
