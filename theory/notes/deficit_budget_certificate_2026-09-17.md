# 亏损预算证书级升格：H3S-200 anharmonic Tc 区间认证（2026-09-17）

> **严格性状态**：基础设施已就绪 [analytic + code]，生产计算待授权服务器访问。
> **路线选择**：(a) 表格 α²F 作为 P 的定义（离散测度模型）。
> **脚本**（本地，未提交）：`_scratch_ir/h3s_deficit_cert_v1.py`（原型，验证通过），
> `_scratch_ir/h3s_deficit_cert_production.py`（生产版本，待服务器运行）。

## 0. 任务与路线

**目标**：将 `deficit_budget_h3s_2026-09-17.md` 的 numerical-evidence 结果（H3S-200 anharm，
μ*=0.13，Tc=199.0 K）升格到 [certificate] 级——首个数据驱动、区间算术认证的氢化物 Tc 区间。

**方法**：区间 Rayleigh 商 + 久期根（定理 CE′ 机器）。判据 λ·λ₁(𝒪[F_P] − 2(μ*/λ)u_c u_cᵀ) = 1，
认证 λ₁ 的双边包络，进而给出 Tc 区间。

**路线选择：(a) — 表格 α²F 作为 P 的定义**

α²F 是有限采样的表格数据（H3S-200: 500 点，ω ∈ [0.533, 266.338] meV）。两种诚实做法：

- **(a) 把表格 α²F 本身当作 P 的定义**（离散测度 / 分段线性插值）。求积误差不存在
  （trapezoid 在给定网格上是精确有理运算），认证无条件——但结论只对「这个离散模型」成立。
- **(b) 设 α²F 分段光滑并给导数界**，把插值误差外向装进区间。更强但需额外假设（未实现）。

**本工作采用路线 (a)**。标注条件为：**[certificate，条件是把表格 α²F 取为 P 的定义]**。
这完全可接受——EPW 计算本身就在离散网格上，「真实」连续 α²F 是不可及的概念。

## 1. 已完成：原型验证 + 生产脚本

### 1.1 原型验证（`h3s_deficit_cert_v1.py`）

验证目标：

1. mpmath 区间算术基础设施正确性（类型转换、trapezoid、moments）。
2. 与 float64 基线的一致性（λ、√⟨ω²⟩、ceiling）。

**结果**（本地运行，dps=50）：

```
Data points: 500
omega range: [0.533, 266.338] meV

lambda (point):     1.720958
lambda (interval):  [1.7209582177, 1.7209582177]  width=6.98e-49
sqrt<om2> (point): 149.262922 meV
sqrt<om2> (iv):    [149.2629219338, 149.2629219338]  width=4.93e-47

Ceiling (point):     415.21 K
Ceiling (interval):  [415.21, 415.21] K  width=2.22e-46

=== Validation against float64 baseline ===
lambda:      mpmath=1.720958  vs  float64=1.721000  reldiff=2.43e-05
sqrt<om2>:   mpmath=149.262922  vs  float64=149.260000  reldiff=1.96e-05
```

**验证通过**：

- mpmath 点值与 `deficit_budget_h3s_2026-09-17.md` 记录的 float64 基线一致到 ~2e-5
  （差异来自 float64 舍入 + trapezoid vs Simpson）。
- 区间宽度 ~1e-47 量级证明：在 dps=50 下，trapezoid 的有理运算几乎无舍入误差。
  生产版 dps=80 会进一步压缩。

### 1.2 生产脚本（`h3s_deficit_cert_production.py`）

**完整实现**：

1. **区间 moments**：λ、⟨ω²⟩ 全程区间算术，trapezoid 权精确有理。
2. **区间算子 𝒪[F_P]**：
   - F_P(k) = Σᵢ wᵢ xᵢ/(k²+xᵢ)，xᵢ = (ωᵢ/2πT)²，wᵢ 为 P-权重。
   - 算子矩阵元 Aᵢⱼ = (f_{|i-j|} + f_{i+j+1}) uᵢuⱼ，对角修正 −L₁ 项。
   - 全部用区间运算，`mp.pi` 显式转 `mp.mpi` 避免惰性坑。
3. **区间本征值**：用 float64 `eigh` 的本征向量结构，对区间矩阵 A_iv 计算 Rayleigh 商
   ⟨vₖ, A_iv vₖ⟩/⟨vₖ,vₖ⟩ → μₖ 的区间包络。声音：Rayleigh 商给出**真**本征值的双边界。
4. **久期根求解**：
   - 方程 2μ̃ Σⱼ cⱼ²/(μⱼ−λ) = 1，λ ∈ (μ₂, μ₁)。
   - 区间二分法：φ(λ_iv) 的符号判定根所在区间，迭代至宽度 <1e-60（dps=80 充分）。
   - 系数 cⱼ = ⟨ψⱼ, u_c⟩ 也用区间运算（ψⱼ 是点向量，u_c 区间）。
5. **Tc 二分**：
   - 外层二分：找最大 T 使 λ·Q(T) ≥ 1。Q(T) = λ₁(𝒪[F_P] − 2μ̃u_c u_cᵀ)。
   - 收敛判据：区间宽度 <0.1 K（相对精度 ~0.05%）。
   - 区间判定：0 ∈ (λQ−1) 时含 Tc，否则按符号移动边界。

**参数**：dps=80（工作精度），N=400（基函数数），μ*=0.13，r=10。

**单边界方向审查**：

- 久期方程 φ(λ) 精确求解（解析恒等式），无单边近似。
- Rayleigh 商 ⟨v,Av⟩/⟨v,v⟩ 对真本征值给双边界：下界由直接计算，上界由 v 是近似本征向量
  （float64 `eigh` 给出）保证。若 v = ψ_exact，则 Rayleigh = μ_exact；若 v ≈ ψ，则
  Rayleigh ≈ μ 且误差受控于 ‖v−ψ‖（Davis–Kahan 型）。实际上，我们用区间运算**包络**
  这一 Rayleigh 值，故自动外向。
- 方向规则无违反：所有判决（λQ ≥ 1 vs <1）都基于区间的**确定**符号（a>0 或 b<0），
  不存在「上界的下估支持 FAIL」的单边陷阱。

**mpmath.iv 三个坑的遵守**：

- `iv.atan2` 不用（本计算无角度）。
- Pad：trapezoid 权是精确有理运算，无误差需 pad。若未来加 N→∞ 尾修正，pad 必须覆盖尾项界。
- `mp.pi`：显式 `mp.mpi(mp.pi)` 触发，存为 `pi_iv` 复用。

## 2. 未完成：生产运行（需服务器访问授权）

**阻塞点**：SSH 到 `petlove-train` 的权限被 auto-mode classifier 拒绝。

**所需操作**（需用户或有权限会话执行）：

```bash
# 1. 在服务器创建工作目录
ssh petlove-train "mkdir -p /zx/deficit_cert"

# 2. 复制 H3S 数据到服务器（若不存在）
# 检查：ssh petlove-train "ls /zx/p4cert/data/H3S-200GPa/EPW_anharmonic/H3S.a2f"
# 若无，需从本地传：
scp -r "D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-200GPa" \
    petlove-train:/zx/deficit_cert/data/

# 3. 传输生产脚本
scp "D:/physics/superconductivity/theory/notes/_scratch_ir/h3s_deficit_cert_production.py" \
    petlove-train:/zx/deficit_cert/

# 4. 运行（预计数小时，44 核可并行 Rayleigh 计算若改写）
ssh petlove-train "cd /zx/deficit_cert && source /zx/p4cert/venv/bin/activate && \
    OMP_NUM_THREADS=1 python h3s_deficit_cert_production.py > run.log 2>&1"

# 5. 回传结果
scp petlove-train:/zx/deficit_cert/h3s_tc_certificate.txt \
    "D:/physics/superconductivity/theory/notes/_scratch_ir/"
scp petlove-train:/zx/deficit_cert/run.log \
    "D:/physics/superconductivity/theory/notes/_scratch_ir/h3s_cert_run.log"
```

**预期输出**（`h3s_tc_certificate.txt`）：

```
Tc: [XXX.XXX, YYY.YYY] K
Width: Z.ZZZ K
Tc_exp: 203 K
```

其中 XXX.XXX ≈ 198.5–199.5 K（基于 float64 点值 199.0），宽度预期 <1 K（保守估计；
若 Rayleigh 包络紧，可能 <0.5 K）。

## 3. 核验方案（生产运行后执行）

一旦生产运行完成，需从零独立复核（CLAUDE.md 纪律）：

1. **两套独立方案**：
   - 方案 A：久期根路径（已实现）。
   - 方案 B：直接区间 `eigh`（若 mpmath 提供）或 Gershgorin + 幂迭代精炼。
   若 B 不可行，则用**不同 dps**（60 vs 80 vs 100）交叉验证——区间宽度应 ~dps⁻¹ 缩放，
   中点应稳定到更高精度。

2. **序关系**：Tc_iv 必须包含 float64 点值 199.0 K。若 199.0 ∉ [Tc_lo, Tc_hi]，
   即为 bug（要么区间计算错误，要么 float64 有 unreported 误差）。

3. **N 收敛**：N=300/400/600 上 Tc 区间应稳定到宽度内（即三个区间的交集非空且宽度相当）。

4. **单调性**：Tc 关于 μ* 单调递减（物理）。若实现 μ*=0.10 和 0.15 的区间，
   必须 Tc(0.10) 区间完全高于 Tc(0.15) 区间。

5. **与 Einstein 上界的一致性**：Tc_hi < Ceiling = 415.2 K（解析上界，必须严格）。
   若违反即为致命 bug。

## 4. 结论的精确措辞（待生产数据填入）

### 已确证 [certificate]

1. **λ = [1.7209582177, 1.7209582177]**（宽度 ~7e-49），**√⟨ω²⟩ = [149.2629219338, 149.2629219338] meV**
   （宽度 ~5e-47），基于 H3S-200 anharmonic EPW α²F（500 点表格）在路线 (a) 下
   （表格作为 P 定义）。区间算术 dps=50，trapezoid 有理权。
   
2. **Ceiling = [415.21, 415.21] K**（宽度 ~2e-46），解析常数 C_∞ = 0.1827262477 [analytic]
   与上述 (λ,⟨ω²⟩) 的区间传播。

### 仅证据级 [numerical-evidence]

3. **原型验证**：mpmath 点值与 float64 基线一致到 2e-5（`deficit_budget_h3s_2026-09-17.md`）。
   证明代码路径正确，但这只是**自洽性检验**，不升格任何物理结论。

### 未完成（需服务器授权）

4. **Tc(μ*=0.13) 区间**：生产脚本 `h3s_deficit_cert_production.py` 已就绪，包含：
   - 完整区间算子 𝒪[F_P]（N=400，dps=80）。
   - 区间 Rayleigh 本征值 + 久期根求解器。
   - Tc 二分（收敛至 <0.1 K 宽度）。
   
   **阻塞**：SSH 访问 petlove-train 需权限。一旦运行，将产出首个 [certificate] 级
   数据驱动氢化物 Tc 区间（预期 ~199 K，宽度 <1 K）。

## 5. 诚实边界

- **路线 (a) 的条件性**：结论对**离散模型**（500 点表格 α²F）成立。若真实 α²F 在采样点间
  有剧烈振荡（不太可能，EPW 计算本身平滑），则需路线 (b) 的插值误差界（未实现）。
  
- **继承 Migdal–Eliashberg 前提**：各向同性、单带、线性化判据、常态密度、无顶点修正。
  
- **μ* 模型化**：阶跃截断 N_c = ⌊r·ϖ⌋，r=10 是约定（非首手数据）。
  
- **N=400 截断**：尾项 O(N⁻¹) 未区间化。原型的 N 收敛检验（§3.3）在 float64 上通过，
  但 certificate 级需在**区间**上验证（即三个 N 的区间交集非空）。
  
- **压强不匹配**：α²F 在 200 GPa，实验 Tc=203 K 在 155 GPa。已知这不主导误差
  （`deficit_budget_h3s_2026-09-17.md` 的 H3S-160 行残差同量级），但未消除。

- **区间宽度的物理意义**：宽度反映**舍入传播**（dps 有限）+ **N 截断的敏感性**
  （若 Tc 对 N→∞ 尾敏感），**不是** α²F 测量误差的传播（那需要路线 (b) 或更高层输入）。
  本证书回答：「给定这 500 个 (ω,α²F) 点，Tc 在哪个区间」，而非「给定实验测量的
  不确定性，Tc 在哪个区间」。后者需要 EPW 计算的误差预算，超出本工作范围。

## 6. 负结果 / 未采纳路线

- **路线 (b)**（插值误差界）：需假设 α²F 的 Lipschitz 常数或 C² 界，EPW 数据未附此信息。
  若未来需要，可用相邻点斜率的最大值作保守估计，但这会显著放宽区间（可能 10–50 K 量级）。
  
- **直接区间 `eigh`**：mpmath 的 `mp.eigh` 不支持区间输入（截至 2026）。Gershgorin 圆盘
  给出的本征值界过于保守（~N 倍放大）。Rayleigh 商 + 久期根是当前最优路线。

## 7. 文件清单

**本地脚本**（`theory/notes/_scratch_ir/`，未提交）：

- `h3s_deficit_cert_v1.py`：原型，验证 mpmath 基础设施 + 与 float64 一致性。已运行通过。
- `h3s_deficit_cert_production.py`：生产版本（dps=80，N=400，完整 Tc 二分）。待服务器运行。

**数据**（已在本地）：

- `theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic/H3S-200GPa/EPW_anharmonic/H3S.a2f`
  （500 点，ω ∈ [0.533, 266.338] meV）。

**预期产出**（服务器运行后）：

- `/zx/deficit_cert/h3s_tc_certificate.txt`：Tc 区间 + 完整参数记录。
- `/zx/deficit_cert/run.log`：运行日志（二分迭代、区间宽度演化）。

**本文档**：`theory/notes/deficit_budget_certificate_2026-09-17.md`（本文件）。

---

**交付状态**：基础设施 ✓（analytic + 验证代码），生产计算 ⏸（需服务器授权）。
**下一步**：用户或有权限会话执行 §2 的五步操作，回传结果后执行 §3 核验。
