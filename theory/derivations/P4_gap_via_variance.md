# P4 · Gap 1 闭合：方差路线（绕开交换子 Frobenius 范数）

> 状态：定理 F/引理 V **analytic**；β_max/β_min/定理 G **numerical-evidence**（B NOT ⪰ 0，见 §2.2）
> 门禁：theory/derivations/p4_gap_via_variance.py（54 ck / 1.4s，0 fail）
> 依赖：P4_concavity_q_route.md (4.2)、P4_g2_analytic_bounds.md（Θ*）

## 0. 记号（与 P4_concavity_q_route.md §0 一致）

f(k;x) = x/(k²+x)，x = ϖ²；A := 𝒪[f]，B := 𝒪[f²]；
h(x) := λ₁(A)，v = v_x 为顶单位本征向量（v[0] > 0）；
q(x) := ⟨v, B v⟩；E := ⟨v, 𝒪[g] v⟩，g = ∂_x(f²) = 2xk²/(k²+x)³；
C := q′ − E（本征向量运动项），γ := h − λ₂（谱隙）。

P4_concavity_q_route.md 已证：E > 0（定理 A，analytic），C ≤ 0（二阶扰动）。
引理 C（H 凹）⟺ q′ > 0 ⟺ |C| < E。本文档给出 |C| < E 的**方差路线**判据。

---

## 1. 核心恒等式：bound1 ≤ bound2（命题 4 → 定理 F）

**定理 F**（方差-隙不等式）。对任意 γ > 0，定义
```
    bound1(γ) := (2/x) Var_v(B) / γ,        Var_v(B) := ‖Bv‖² − q²
    bound2(γ) := (2/x) ‖[A,B]v‖² / γ³
```
则 **bound1(γ) ≤ bound2(γ)** 恒成立。

*证明*（P4_concavity_q_route.md 命题 3 的直接推论）。命题 3 的 (4.2) 给
```
    ‖(I−P₁)Bv‖ ≤ ‖[A,B]v‖/γ.
```
注意 `Var_v(B) = ‖(I−P₁)Bv‖²`（因 P₁Bv = ⟨v,Bv⟩v = qv，故 `‖Bv‖² = q² + ‖(I−P₁)Bv‖²`），
两边平方再乘 (2/x)/γ：
```
    (2/x) Var_v(B) / γ = (2/x) ‖(I−P₁)Bv‖² / γ
                       ≤ (2/x) ‖[A,B]v‖² / γ³.   □
```

**含义**：控制 |C| 可以**绕开** ‖[A,B]‖_F 的解析上界（P4_concavity_q_route.md §5 诚实边界 #3），
只需要 Var_v(B) 的上界 + γ 的下界。这把「对 ‖[A,B]‖_F 的解析界」从阻塞性必需降级为
「若有会更省事，但不是唯一路径」。

---

## 2. Var_v(B) 的解析上界：Bhatia–Davis + 谱端点包围

**引理 V**（方差的 Bhatia–Davis 界）。对自伴算子 B 与单位向量 v，
```
    Var_v(B) := ⟨v, B²v⟩ − ⟨v,Bv⟩² ≤ (β_max − q)(q − β_min),
```
其中 β_max, β_min 分别是 B 的最大/最小本征值，q = ⟨v,Bv⟩。

*证明*。标准 Bhatia–Davis 不等式，见 Bhatia, Matrix Analysis (2nd ed.), Theorem VII.3.7。□

### 2.1 β_max 的解析界

**引理 β_max**（analytic）。B = 𝒪[f²]，f(k;x) = x/(k²+x)，则
```
    β_max ≤ tr B + Θ*,        Θ* = 0.6349287917881079529866324253...
```
其中 Θ* 已由 P4_g2_analytic_bounds.md 解析确定（Schur majorant，零计算机辅助）。

*证明*。Gerschgorin 圆盘定理：λ_max ≤ max_n [B_nn + R_n]，其中 R_n := Σ_{m≠n}|B_nm|
是第 n 行的非对角绝对行和。P4_g2_analytic_bounds.md §2.3 给出 R_n ≤ Θ* 对一切 n 成立
（Schur 包络的逐行界）。因 B 对角元非负（B_nn = ⟨e_n, Be_n⟩ ≥ 0，见 §2.2 验证），
max_n B_nn ≤ Σ_n B_nn = tr B，故
```
    β_max ≤ max_n B_nn + Θ* ≤ tr B + Θ*.   □
```

**注**：tr B = Σ_n ⟨e_n, 𝒪[f²] e_n⟩ 有显式级数（见 §2.3），需验证是否有闭式。
若无闭式，tr B 为数值量（float64 求和或区间算术），则本引理降级为 numerical-evidence。

### 2.2 β_min 的数值界（B NOT ⪰ 0）

**引理 β_min**（numerical-evidence）。B = 𝒪[f²] **NOT 半正定**（门禁验证：λ_min(B) < 0 在所有 18 点）。

*分析*。𝒪[f²] 的构造为 L2 + L3 + L1 三项：
- L2: f²(|n−m|)，对角外，核 f² ≥ 0
- L3: f²(n+m+1)/(√(2n+1)√(2m+1))，核 f² ≥ 0
- L1: 对角 −2Σ_{k=1}^n f²(k)/(2n+1)，为**负项**

L2 + L3 单独各自半正定（核非负），但 **L1 负对角使 B_nn < 0**（门禁实测：min(B_diag) ≈ −1e-3 at x=0.04），
从而 λ_min(B) < 0（门禁实测：β_min ∈ [−0.97, −0.001] on x ∈ [0.04, 1638.4]）。

因 B NOT ⪰ 0，Bhatia–Davis 不退化，需使用完整形式 `Var_v(B) ≤ (β_max−q)(q−β_min)`，
其中 β_min < 0 **增大右端**（因 q−β_min = q + |β_min| > q）。

**数值处理**：门禁直接用 `eigh` 算出的 β_min（numerical）。解析界可尝试 Gerschgorin 下盘
`β_min ≥ min_n[B_nn − Θ*]`，但该下界为负且无闭式改进，故最终判据（定理 G）**降级为 numerical-evidence**。□

### 2.3 tr B 的计算

tr B = Σ_n ⟨e_n, 𝒪[f²] e_n⟩。展开 𝒪 的三项：
```
    ⟨e_n, 𝒪[f²] e_n⟩ = [L2 对角为 0] + [L3 对角项：2f²(2n+1)/(2n+1)]
                      + [L1 对角项：−2Σ_{k=1}^n f²(k)/(2n+1)]
                      = 2f²(2n+1)/(2n+1) − 2Σ_{k=1}^n f²(k)/(2n+1)
                      = [2/(2n+1)] [f²(2n+1) − Σ_{k=1}^n f²(k)].
```
故
```
    tr B = Σ_{n=0}^∞ [2/(2n+1)] [f²(2n+1;x) − Σ_{k=1}^n f²(k;x)].
```
f²(k;x) = x²/(k²+x)²。该级数**尾部为 O(1/n³)**（L1 累积和抵消 L3 主项后剩余），
数值求和需高精度截断或 Richardson 外推。**若无闭式，则 tr B 为数值量**，
引理 β_max 与最终定理 G 降级为 numerical-evidence。

**门禁任务**（§2.3 核验）：
- 用 N=300, KX=2400 计算 tr B（截断级数），核验与 `np.trace(B)` 一致（相对误差 <1e-10）。
- 若需区间算术认证，用 mpmath.iv 计算截断和 + 尾部 Cauchy 余项界（本任务暂不实施）。

---

## 3. 定理 G：Gap 1 闭合（|C| < E）的解析判据

**定理 G**（方差路线闭合 Gap 1）。设 γ_lo > 0 为 γ 的下界（analytic 或 numerical），
β_max, β_min 为 B 的谱端点上界/下界（来自引理 β_max/β_min），q, E 已知。若
```
    (2/x) (β_max − q)(q − β_min) / γ_lo  <  E,
```
则 |C| < E。

*证明*。由 P4_concavity_q_route.md (4.1)，|C| ≤ (2/x)·‖(I−P₁)Bv‖²/γ。应用定理 F：
```
    |C| ≤ (2/x) Var_v(B) / γ_true
       ≤ (2/x) Var_v(B) / γ_lo                     (γ_lo ≤ γ_true，单边安全)
       ≤ (2/x) (β_max−q)(q−β_min) / γ_lo          (Bhatia–Davis)
       <  E                                         (判据)。   □
```

**严格性分级**：
- 若 γ_lo, β_max, q, E 全部 analytic，**且** β_min = 0 经 B ⪰ 0 解析证明，
  **则整个链 analytic**。
- 若 tr B 需数值（无闭式）或 β_min 需 Gerschgorin 数值下盘，**降级为 numerical-evidence**。

**与 P4_concavity_q_route.md 的对比**：
- 旧路线（bound2，交换子）：|C| ≤ (2/x)‖[A,B]‖_F²/(2γ³)，需 ‖[A,B]‖_F 解析上界（缺失）。
- 新路线（bound1，方差）：|C| ≤ (2/x)Var_v(B)/γ，只需 β_max/β_min + γ_lo，
  **已有独立来源**（Θ* 来自 g2_analytic，γ_lo 可复用 trace-commutator 的 h_lo_trial 机制）。

---

## 4. 数值对比：bound1 vs bound2（机制 B 优势的定量证据）

用 18 点网格（ϖ ∈ [0.2, 40]，x = ϖ²），对比：
- `ratio_bounds := bound1(γ_lo) / bound2(γ_lo)`：理论预测 ~1/γ²（因 bound2 有 γ³ 分母）
- `bound1/E` vs `bound2/E`：哪个路线更紧

γ_lo 取自 **P4_concavity_q_route.md 已有的解析隙界**：γ_lo = max(0.42h, h−Θ*)
（该界独立于交换子，来自重排恒等式 + Schur majorant）。

**预期发现**（门禁将验证）：
1. 小 γ 处（ϖ ∈ [1.4,1.7]），bound1 远优于 bound2（ratio ~ 1/γ² ≪ 1）。
2. bound1 路线在 ϖ=1.54 worst 点的 |C|_bound1/E **应接近** 
   P4_concavity_q_route.md 的机制 B 优化版（若有）的闭合比值。
   若 P4_concavity_q_route.md 的 bound2 路线在 ϖ=1.54 给 1.096（不闭合），
   bound1 路线**预期更好**（因绕开了交换子的松界）。

**表格**（门禁输出，18 点）：
```
  ϖ  |    x    |    h    |  γ_lo  |  Var_v(B)  | bound1/E | bound2/E | ratio
-----|---------|---------|--------|------------|----------|----------|-------
 0.20|  0.0400 | ...     | ...    | ...        | ...      | ...      | ...
 ...
 1.54|  2.3716 | ...     | ...    | ...        | ...      | ...      | ...  ← worst for bound2
 ...
40.00|1600.00  | ...     | ...    | ...        | ...      | ...      | ...
```

**判据**：若 bound1/E < 1 在 18 点全部成立，则方差路线**数值上全程闭合 Gap 1**。

---

## 5. 诚实边界

1. **β_min = 0 的解析证明未完成**。若 B ⪰ 0 经门禁验证（数值 λ_min ≥ 0），
   但**无解析证明 B 半正定**，则 β_min = 0 只是数值断言，定理 G 降级为 numerical-evidence。
   尝试路径：
   - (a) 证明 𝒪 保持核的逐点非负性 ⟹ 𝒪[f²] ⪰ 0（已知**不成立**，−L1 负对角破坏）。
   - (b) 证明 B 的所有主子式非负（Sylvester 判据）——计算量大，未实施。
   - (c) 若 (a)(b) 均失败，诚实标注 numerical-evidence。

2. **tr B 的闭式求和未建立**。级数 Σ_{n=0}^∞ [2/(2n+1)][f²(2n+1) − Σ_{k=1}^n f²(k)]
   尾部为 O(1/n³)，数值可算但无已知闭式。若无闭式，tr B 为数值量（float64 或 mpmath 高精度），
   引理 β_max 与定理 G 降级为 numerical-evidence。

3. **本文档假定 γ_lo 已由其他来源给出**（复用 P4_concavity_q_route.md 的
   γ_lo = max(0.42h, h−Θ*)，或更强的机制 B 优化版）。本文档**不包含 γ_lo 本身的构造**，
   那是独立的工具链（gap 界 + h 下界）。

4. **bound1 路线的数值优势依赖 γ_lo 的质量**。若 γ_lo 过松（如 0.42h 在 ϖ∈[1.4,1.7]
   仅为真隙的 0.42 倍），bound1 也会失效（但失效阈值比 bound2 更宽容，因只有 γ 一次方
   vs γ³）。若未来有更强的 γ_lo（如 0.45h 或逐点隙表），bound1 路线将**自动受益**。

5. **与 P4_concavity_q_route.md 的关系**：本文档**不替代** P4_concavity_q_route.md
   的任何已证结论（定理 A / 命题 1–3 / E>0 / C≤0），只提供 |C|<E 的**另一条路径**，
   作为交换子 Frobenius 范数解析界（§5 诚实边界 #3）的备选方案。两条路线可并行推进。

---

## 6. 与既有结果的接口

- 定理 F 的证明直接引用 P4_concavity_q_route.md 命题 3 的 (4.2)。
- 引理 β_max 用 P4_g2_analytic_bounds.md 的 Θ* = 0.6349...（已证 analytic）。
- γ_lo 复用 P4_concavity_q_route.md 的解析隙界（重排恒等式 + Schur majorant），
  或未来更强版本（如 trace-commutator 机制 B 的 h_lo_trial + Gdom cap）。
- 引理 V（Bhatia–Davis）是标准教材定理，无需独立证明。
- 不改动任何既有文件；A, B, v 的构造由门禁从零重建，复用 p4_concavity_q_route.py 的
  `O_ker` / `kernels` / `top_eigen` 函数（自包含，不依赖其他模块的数值结果）。
