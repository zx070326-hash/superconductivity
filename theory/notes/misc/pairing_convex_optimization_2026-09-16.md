# 配对问题的凸优化框架探索

建立：2026-09-16
状态：v0.1（框架搭建与初步分析，未经对抗核）

---

## 0. 动机与背景

Einstein 极值定理（定理 E / 定理 2，`P4_einstein_extremal.md`，全部经对抗核）：
> 固定 lambda 和 ⟨omega^2⟩ 时，Einstein 谱（单频 delta_{omega_E}）使 Tc 最大。

证明路线经过两个阶段：
- **旧路线**（引理 C + Jensen）：若 H(x) = lambda_1(script-H(sqrt{x})) 凹，则由 Jensen 得结论。需要证书段 (0,45] + 解析段 [8.52, infty) 的混合证明，52 个常数。
- **新路线**（引理 M + 核级 Jensen）：x -> x/(k^2+x) 逐 k 严格凹，由引理 M（核单调性，锥不变性 + Abel 分解）一步给出。三行解析证明，零计算机辅助。

核心问题：这个极值结构能否形式化为**凸优化问题**？如果能，则：
1. 理论上，给出更一般的极值刻画（不只 Einstein，也不只固定 ⟨omega^2⟩）。
2. 应用上，SDP / 内点法可做声子谱的逆向设计（给定 Tc 目标反推谱）。
3. 对偶理论能提供新的物理洞察（Lagrange 乘子的物理意义）。

---

## 1. 形式化为优化问题

### 1.1 变量空间

变量：谱测度 P in P = {P 是 R_+ 上的概率测度，⟨omega^2⟩ := int omega^2 P(d omega) < infty，P != delta_0}。
P 是 C_0(R_+)^* 的闭凸子集（弱-* 拓扑下紧化需要在 R_+ 的一点紧化上加条件，但 ⟨omega^2⟩ < infty 给出了充分的尾部控制）。

### 1.2 目标泛函

**目标**：max Tc(lambda, P)。

线性化 Eliashberg 框架（KAY）中 Tc 的判据是 lambda * kappa(P,T) = 1，kappa(P,T) := lambda_1(K(P,T))。由定理 1(b)（T -> kappa(P,T) 严格递减，已证），Tc 良定义：

  Tc(lambda, P) = 唯一的 T 使 lambda * kappa(P,T) = 1

等价地，固定 T，最大化 kappa(P,T) 等价于最大化 Tc。

### 1.3 约束

- **归一化**：int P(d omega) = 1
- **矩约束**：int omega^2 P(d omega) = C（固定二阶矩）
- **非负性**：P >= 0（测度）
- **支撑**：P((0,infty)) = 1（排除 delta_0）

### 1.4 可行域的几何

可行域 K_C := {P in P : ⟨omega^2⟩ = C} 是弱-* 紧的（因为在 {⟨1⟩=1, ⟨omega^2⟩=C} 约束下由 Prokhorov 定理紧）。极点集 = {delta_omega : omega > 0, omega^2 = C} = {delta_{sqrt{C}}} 是单点——这在 Choquet 意义下就是 Einstein 谱。

**注意**：这不是标准 Choquet 定理的直接推论。Choquet 定理要求的是"紧凸集上的连续仿射函数在极点取极值"，而 kappa(P,T) 既非 P 的仿射函数（它是 lambda_1，凸的 sup），也非在弱-* 拓扑下连续（只是上半连续）。所以 Einstein 极值性不从 Choquet 免费得来；引理 M 给出的是更强的结论。

---

## 2. 凸性分析

### 2.1 kappa(P,T) 关于 P 的凸性

由 P4_main_theorems.md §1 注（Remark 3.17'）：

  kappa(P,T) = lambda_1(O[F_T]) = sup_{||v||=1} ⟨v, O[F_T] v⟩

其中 F_T(k) = int [k](omega/2piT) P(d omega) 关于 P 线性。因此 kappa 是 P 的**线性泛函族的上确界**，故**关于 P 凸**。

**推论**：max_{P in K_C} kappa(P,T) 是一个**凸函数在紧凸集上的最大化**问题。凸函数的极值在极点取到——这与 Choquet 方向一致。但如 §1.4 所述，K_C 的唯一极点就是 delta_{sqrt{C}}，所以 Einstein 最优性在这个层面是一致的。

### 2.2 kappa(P,T) 不是 P 的凹泛函

这一点至关重要：如果 kappa 是凹的，则 Tc 的最大化（等价于最大化凹泛函）会是标准凸优化。但 kappa 是**凸**的，最大化凸函数是 NP 难的一般问题。

幸运的是，我们不需要一般性：约束集 K_C 只有一个极点（固定 ⟨omega^2⟩ 的约束太强），所以极值的唯一性来自集合的几何，而不是目标的凹性。

### 2.3 Frechet 导数

设 P_t = (1-t)P_0 + t P_1（P_0, P_1 in K_C），则 F_T^t = (1-t)F_T^0 + t F_T^1，故

  d/dt|_{t=0} kappa(P_t, T) = ⟨v_0, O[F_T^1 - F_T^0] v_0⟩

其中 v_0 是 K(P_0,T) 的顶本征向量。即 Frechet 导数是线性的、由顶本征向量的二次型给出。

**二阶**：
  d^2/dt^2 kappa(P_t, T) = 2 ⟨r, (lambda_1 - K)^+ r⟩ >= 0

其中 r = (I - P_0) O[F_T^1 - F_T^0] v_0，(lambda_1 - K)^+ 是约化预解式。**严格非负**（除非 r=0），确认凸性。

### 2.4 引理 C 视角下的"拟凹性"

虽然 kappa 关于 P 是凸的，但 H(x) = h(sqrt{x}) 关于 x 是凹的（引理 C，已证）。这意味着：

  int H(varpi^2) P(d omega) <= H(int varpi^2 P(d omega)) = H(varpi_{rms}^2) = h(varpi_{rms})

即 kappa(P,T) <= int h(varpi) P <= h(varpi_{rms})。**第一个不等式是凸性给的，第二个是凹性给的**。引理 M 绕开了第二步，直接给出更强的 kappa <= h(varpi_{rms})。

---

## 3. 对偶问题

### 3.1 Lagrange 对偶

原始问题：
  maximize_{P in K_C} kappa(P, T)
  subject to: int P = 1, int omega^2 P = C, P >= 0

引入 Lagrange 乘子 mu（归一化）和 nu（矩约束）：

  L(P; mu, nu) = kappa(P,T) - mu(int P - 1) - nu(int omega^2 P - C)

对偶函数：
  g(mu, nu) = sup_{P >= 0} L(P; mu, nu)
            = mu + nu C + sup_{P >= 0} [kappa(P,T) - int (mu + nu omega^2) P(d omega)]

**关键观察**：kappa(P,T) = sup_{||v||=1} ⟨v, O[F_T] v⟩ = sup_v int ⟨v, H(varpi) v⟩ P(d omega)（逐元线性 + Rayleigh 商）。设 q(omega; v) := ⟨v, H(omega/2piT) v⟩，则

  sup_{P >= 0} [int q(omega; v) P - int (mu + nu omega^2) P] = sup_{P >= 0} int [q(omega; v) - mu - nu omega^2] P

内层 sup 对 P 的结果：若对某个 omega_0 有 q(omega_0; v) > mu + nu omega_0^2，则取 P = c delta_{omega_0}（c -> infty）使目标 -> +infty。因此对偶有限 iff 对所有 omega, v 有

  **q(omega; v) <= mu + nu omega^2**，即 ⟨v, H(varpi) v⟩ <= mu + nu varpi^2 对所有 varpi, ||v||=1。

取 v 为 H(varpi) 的顶本征向量：**h(varpi) <= mu + nu varpi^2 对所有 varpi > 0**。

### 3.2 对偶的物理意义

对偶问题：
  minimize_{mu, nu} mu + nu C
  subject to: h(varpi) <= mu + nu varpi^2 对所有 varpi > 0

即：**找一个二次函数 mu + nu x（x = varpi^2）使它包络（majorize）Einstein 函数 h(varpi)，且在约束点 x = C/(2piT)^2 处最小**。

这正是 h(varpi) 的**凸包络**（convex hull / concave envelope）！如果 H(x) = h(sqrt{x}) 是凹的（引理 C），则凸包络就是 H 本身，对偶最优值 = H(⟨varpi^2⟩) = h(varpi_{rms})，且取等条件给出 nu = H'(⟨varpi^2⟩)，mu = H(⟨varpi^2⟩) - ⟨varpi^2⟩ H'(⟨varpi^2⟩)——这些是 H 在 ⟨varpi^2⟩ 处的**切线参数**。

**物理解读**：
- **nu**（矩约束的 Lagrange 乘子）= H 在 ⟨varpi^2⟩ 处的导数 = 谱矩的"影子价格"：增加 ⟨omega^2⟩ 一个单位，Tc 判据改善 nu 个单位。
- **mu**（归一化的 Lagrange 乘子）= 切线截距 = 在 omega=0 处的"虚拟价值"。
- **强对偶**：由引理 C（H 凹），对偶间隙为零。即 max_P kappa(P,T) = min_{mu,nu} [mu + nu C] = h(varpi_{rms})。这就是定理 E 的对偶表述。

### 3.3 Slater 条件

Slater 条件要求可行域有内点。K_C 在弱-* 拓扑下确实有内点（任何绝对连续的 P 且 ⟨omega^2⟩ = C 都是内点）。但严格来说，此处对偶性的成立不需要 Slater——强对偶由 H 的凹性 + Fenchel 对偶直接给出。

---

## 4. 与 Choquet 定理的关系

### 4.1 极点定理

**命题**：K_C = {P in P : ⟨omega^2⟩ = C} 的极点恰为 {delta_{sqrt{C}}}。

证明：设 P = alpha P_1 + (1-alpha) P_2（alpha in (0,1)，P_1, P_2 in K_C）。若 P = delta_{omega_0}（omega_0 = sqrt{C}），则 P_1 和 P_2 都必须集中在 {omega_0}（因为 alpha P_1 <= P = delta_{omega_0}），故 P_1 = P_2 = delta_{omega_0}。反过来，若 P 有不止一个支撑点，则可以分裂为两个不同的测度仍满足约束——因此 P 非极点。

**推论**：固定 ⟨omega^2⟩ 的约束把极点坍缩成单点，极值只能在那里取到。这**不是**因为 kappa 是凸的或凹的，而是因为**约束集的极端结构**。

### 4.2 多约束下的极点

如果放松到固定 ⟨omega^2⟩ 和 ⟨omega^4⟩（即固定前两个矩），则极点集变成：
- 单频 delta_omega（一个参数）
- 双频 a delta_{omega_1} + (1-a) delta_{omega_2}（三个参数 a, omega_1, omega_2，两个矩约束留一个自由度）

此时极值**不一定**在 Einstein 谱取到——可能在某些双频谱处取更大值。这与 P3 的 lambda_star 展开（kappa = g(2) 的 P-平均展开，谱形修正 sigma^2, gamma_3 sigma^3）吻合：宽谱可以增大或减小 Tc，取决于 varpi 区域与 kappa(varpi) 函数的精细形状。

### 4.3 Bergmann-Rainer 泛函导数

文献定位（researcher `extremal-lit`，已并入 `P4_einstein_extremal.md` §0.A）：Bergmann-Rainer Z. Physik 263, 59 (1973) 给出了泛函导数 delta Tc / delta alpha^2F(omega) > 0 与峰位的变分结果。他们的框架是**固定面积** A = lambda ⟨omega⟩/2（Leavens-Carbotte 约束），与本文的**固定 ⟨omega^2⟩** 不同——两者回答的是不同的优化问题。

Semenok-Altshuler-Yuzbashyan 2407.12922 SI §S4 固定 omega_max，变分非严格。Combescot PRB 42, 7810 (1990) 固定 omega_log（仅 lambda -> 0）。这些都是不同约束下的优化，与本文的 ⟨omega^2⟩ 约束互补。

---

## 5. 算法可行性

### 5.1 有限维离散化

将谱离散为 N 个频率 {omega_1, ..., omega_N}，P = sum_i p_i delta_{omega_i}：
- 变量：p = (p_1, ..., p_N) in R^N_+
- 约束：sum p_i = 1, sum p_i omega_i^2 = C
- 目标：lambda_1(sum_i p_i H(omega_i/2piT))

这是一个 **N 维有限集上的凸最大化问题**——在凸体的极点（= 顶点）取极值。离散版的顶点就是 e_j（只选一个频率的基向量），约束 omega_j^2 = C 选定 j。

但如果 omega_j 不精确等于 sqrt{C}，则需要两个频率的混合——最优解在 K_C 的边上，是一个二元 delta 的混合。这对应 SDP 的一个角。

### 5.2 SDP 形式（当前不可行的方向）

如果把 lambda_1(A(p)) >= t 写成 SDP 约束 A(p) >= t I（LMI），则：
  maximize t
  subject to: sum_i p_i H_i >= t I, sum p_i = 1, sum p_i omega_i^2 = C, p >= 0

这是一个标准的 SDP（半定规划），可用 MOSEK / SCS / SDPA 求解。

**但有两个实际障碍**：
1. **H_i = H(omega_i/2piT)** 是无穷维矩阵，必须截断到 N_trunc x N_trunc。截断引入的误差由 P4 的换序一致性定理控制（`P4_interchange_uniformity.md`：kappa^(N) -> kappa 带显式尾率），但 N_trunc 必须足够大（~10 varpi 以上，对弱耦合 varpi ~ 100-1000 则需 N_trunc ~ 10^3-10^4）。
2. **温度 T 不固定**——Tc 是 lambda kappa(P,T) = 1 的隐定义，优化时 T 也是变量，但 T 以高度非线性的方式进入 H。

**当前评估**：SDP 在 N_trunc ~ 100 时可行（秒级），但物理相关的参数需要 N_trunc ~ 1000+，这对 SDP solver 是边界情况。更重要的是，定理 2 已经给出了闭式答案（Einstein 谱最优），SDP 的计算价值在于**验证**和**多约束扩展**，而非发现新极值。

### 5.3 与 DFT + EPW 工作流的接口

逆向设计的真实场景：给定 lambda, omega_max, ⟨omega^2⟩ 的窗口，找使 Tc 最大的 alpha^2F(omega)。这需要：
1. 用本文的框架（或 SDP 有限维版本）确定最优谱形。
2. 用 DFT + EPW 的 alpha^2F(omega) 计算检验该谱形是否可达。
3. 若不可达，在可达集内投影（最近点）。

这一链条**目前不可行**：alpha^2F(omega) 是材料性质而非设计变量，可达集取决于晶格结构 + 化学组成 + 压力，远非简单约束。但它提供了一个**评估方向**：给定 alpha^2F(omega)，其 Tc 占理论上限（定理 E 给的 Einstein 极值）的百分比是一个有意义的效率指标。

---

## 6. 与经典优化文献的联系

### 6.1 谱测度优化

本问题属于 **infinite-dimensional linear programming with spectral measure variables**，与以下领域相关：
- **矩问题**（Hausdorff, Hamburger）：固定 ⟨omega^k⟩ 的约束正是截断矩问题。Richter (1957) 的极值定理说 n 个矩约束下的极值测度最多 n+1 个原子——这给出 §4.2 的结论。
- **D-optimal design**（统计学）：选择实验设计点以最大化 Fisher 信息矩阵的行列式，结构完全类似（目标是矩阵函数，变量是离散测度）。Kiefer-Wolfowitz 定理是那边的 Choquet 类比。
- **最优控制中的 relaxed controls**：Young 测度 = 参数化测度，优化目标是凸关于测度的——同类问题。

### 6.2 矩阵凸性 / Loewner 理论

P4_main_theorems.md §1 注明确指出：lambda_1 关于**算子**是凸的（Rayleigh 商定义），关于**核**也是凸的（O 的线性性）。引理 C 要求的 H(x) 的凹性与 lambda_1 的凸性**方向相反**——这正是引理 C 困难的结构性原因。

Loewner 矩阵单调/凸性理论（f: R -> R 是算子凸 iff f 的 Loewner 矩阵非负定）**不直接适用**，因为 varpi -> H(varpi) 不是算子凸的（`P4_einstein_extremal.md` §0 探针确认 lambda_min(D) ~ -0.47）。

---

## 7. 可行性评估与下一步

### 7.1 已建立的（硬事实）

1. **kappa(P,T) 是 P 的凸泛函**——`P4_main_theorems.md` Remark 3.17'，解析。
2. **H(x)=h(sqrt{x}) 在 (0,infty) 凹**——引理 C，计算机辅助定理。
3. **Einstein 谱是 ⟨omega^2⟩ 固定下的 Tc 极大者**——定理 2，解析。
4. **对偶间隙为零**（强对偶成立）——由引理 C + Fenchel 直接。
5. **Lagrange 乘子 nu = H'(⟨varpi^2⟩)** 有明确的物理意义。

### 7.2 值得推进的方向

1. **多约束极值**：固定 ⟨omega^2⟩ 和 ⟨omega^4⟩（或 omega_max）时，极值谱是**双频**的。引理 M 不直接覆盖（它只给出与 Einstein 的比较），但 Richter 极值定理给出结构。数值上可以用 SDP（小规模）或参数扫描确定最优双频谱。
2. **效率指标** eta(P) := Tc(lambda, P) / Tc^E(lambda, sqrt{⟨omega^2⟩})——定理 E 保证 eta <= 1，它量化了谱展宽对 Tc 的损耗。对具体材料（H3S, LaH10）算 eta 可以回答"谱形还有多少优化空间"。
3. **逆问题的松弛**：给定目标 Tc*，问是否存在满足矩约束的 P 使 Tc >= Tc*。这是对偶问题 nu C + mu >= kappa*（kappa* = 1/lambda 在 Tc*）的可行性判定，可用对偶 LP 在矩约束内求解。
4. **mu* > 0 的优化**：库仑赝势使问题**非凸**（K_c 的秩一结构破坏了 O[F_T] 的线性性），但算子序给出上界。mu* > 0 下的极值谱是否仍是 Einstein？定理 E 的推论 (v) 说上界成立，但极值的取到（下界）是开放的。

### 7.3 不值得推进的方向

1. **通用 SDP solver 暴力求解**——定理 2 已给出闭式答案，SDP 能做的只是数值验证。
2. **Loewner 型论证**——已被 §0 探针排除（H 不是算子凸或算子凹）。
3. **写软件做"材料逆向设计"**——alpha^2F 的可达集取决于材料自由度，远非本框架能控制。

### 7.4 时间评估

- 效率指标 eta 的计算：1-2 天（用已有的 Eliashberg 求解器 + 文献 alpha^2F）。
- 双频极值的参数扫描：2-3 天（数值，无需新证明）。
- 多约束极值的严格理论：数周级（需要新的核级 Jensen 类论证或引理 M 的扩展）。

---

## 8. 诚实边界

1. 本文是**探索性框架搭建**，没有新定理。所有引用的定理状态以 `P4_main_theorems.md` v0.8 与 `P4_einstein_extremal.md` 为准。
2. 对偶问题的强对偶性**依赖引理 C**（H 凹性），后者目前只有计算机辅助证明。但引理 M 给出的 kappa <= h(varpi_{rms}) 已经是解析的，所以**定理 E 本身不需要引理 C 或对偶理论**——对偶提供的是理论框架，不是更强的结论。
3. mu* > 0 下的凸优化框架**严格更难**：K(P,T) - mu_tilde K_c 不再关于 P 线性（K_c 的 N_c(T) 阶梯破坏线性），引理 M 的单调性论证不直接搬到束上。
4. "SDP 可做材料逆向设计"的陈述**极度乐观**：物理约束（晶格结构、稳定性、声子色散）未被任何凸约束捕获。

---

## 9. 本文与已有结果的关系总表

| 已有结论 | 严格性 | 本文对应 |
|---|---|---|
| kappa 关于 P 凸 | 解析（§1 注 / Remark 3.17'） | 凸优化的目标结构 §2.1 |
| 定理 E：Einstein 极值 | 解析（定理 2，核级 Jensen） | 原始问题的闭式解 §1 |
| 引理 C：H 凹 | 认证（证书 + 解析） | 对偶间隙为零的充分条件 §3.2 |
| 引理 M：核单调 | 解析（锥 + Abel） | 极值性的直接证明（绕开对偶） §2.4 |
| Richter 极值 | 经典概率（1957） | 多约束极点结构 §4.2 |
| KAY 换序 | 认证（区间定理 [0.2,60]） | SDP 截断误差控制 §5.2 |
| Bergmann-Rainer 变分 | 文献（1973） | 不同约束下的优化 §4.3 |
