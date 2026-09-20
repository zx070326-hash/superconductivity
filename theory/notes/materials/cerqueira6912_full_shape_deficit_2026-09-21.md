# Cerqueira 6912 数据管线：谱形状亏损分布与定理 2 紧检验

> 2026-09-21  
> 状态：**BLOCKED — 关键前提被证伪**  
> 分支：`p4-conj2`（worktree agent-a251c94abbc4f0ae2，UNMERGED）

## § 0  执行摘要（前提验证结果）

### 0.1  关键发现：Eliashberg 数据在高 λ 材料上几乎全部缺失

**前提假设（主会话手上 0/249 行证据）：** Cerqueira 数据集在高 λ 材料上会有 Eliashberg 全谱数值解 `Tc(μ*)`，从而可检验天花板与定理 2 的 λ 分辨界。

**验证结果（batch-p，500 材料，λ ∈ [0.107, 0.213]，中位 0.168）：**
- **Eliashberg 有数：2/500 材料（0.4%）**
- 有数的两个：`H6Ca2NiPd_agm004930255`（λ=0.191），`H6Sr2RhIr_agm004930587`（λ=0.195）
- **μ* 覆盖：仅 μ*=0.00 一个点，其余 19 个 μ* 网格点全为 `-----`**
- batch-q（249 材料，λ ≤ 0.21，**全部弱耦合**）：0/249 有 Eliashberg 数据

**推论：**
1. **「定理 2 紧检验」（推广点 2）的前提为假**——数据集不提供 `Tc_Eliashberg(μ*=0)` 在足够多材料上的覆盖，无法做 `𝔣/𝔣_E(λ)` 的分布统计或违反扫描。
2. **仅剩天花板检验**——`Tc ≤ C_∞√λ·ω̄₂` 可用 McMillan/AD 列作为管线自检（它们是 λ、ω_log、ω̄₂ 的确定性函数，零新信息），但**不构成独立确证**。
3. **唯一非退化检验路线：SuperCon 交叉匹配**——把 Cerqueira 的 (λ_DFT, ω̄₂) 与 SuperCon 实验 Tc 交叉，得到「非隧道 (λ, ω̄₂, Tc_expt) 家族」，这是推广点 1 的唯一外部效度来源。

---

## § 1  数据集概况

### 1.1  Cerqueira et al. 2023 数据集
- **来源：** Materials Cloud 3kbt5-r3n56（17 批次，6912 材料）
- **文件结构（每材料）：**
  - `McMillan.dat`：表头给 λ、ω_log[K,meV]、ω̄₂[K,meV]；数据行给 `(μ*, Tc_McMillan, Tc_AD)` 对 μ* 网格
  - `Eliashberg.dat`：表头 `# mustr - TcEliashberg`；数据行 `(μ*, Tc_Eliashberg)` 或 `-----`（求解失败）
  - `a2F.dos1..10`：不同展宽的 α²F(ω)
- **本地已有：** batch-q（249 材料）
- **服务器（petlove-train）：** `/zx/cerqueira6912/`，16 批已下载并解压（batch-p 先到，500 材料）

### 1.2  解析器验证（batch-q 负对照）
- **工具：** `tools/cerqueira_collect.py`（输出 CSV：batch, material_dir, formula, mcid, lambda, wlog_K/meV, w2av_K/meV, r, Tc_ME_mu000/010/013, Tc_McM_mu010, Tc_AD_mu010, n_ME_numeric）
- **验证：** 对 batch-q 249 材料，λ、ω_log、ω̄₂、r 与参考 CSV 逐位一致（PASS）
- **batch-q 特征：**
  - λ ∈ [0.024, 0.209]，**全部弱耦合**（Tc_McM/AD 在 μ*=0.10 给 0 或空）
  - Eliashberg 列：0/249 材料有数（`n_ME_numeric=0`）
  - ⟹ 负对照符合预期：弱耦合材料 Eliashberg 求解器不跑或全 `-----`

### 1.3  batch-p 验证结果（§ 0.1 的详细数据）
- **材料数：** 500
- **λ 分布：** 最小 0.107，中位 0.168，p90 0.204，最大 0.213
- **Eliashberg 有数：** 2/500（0.4%）
  - `H6Ca2NiPd_agm004930255`：λ=0.191，`Tc(μ*=0)=0.998 K`，其余 19 个 μ* 点全 `-----`
  - `H6Sr2RhIr_agm004930587`：λ=0.195，`Tc(μ*=0)=0.858 K`，同样只有 μ*=0.00 一点
- **μ* 网格：** 0.00, 0.02, 0.04, …, 0.38（步长 0.02，20 点）
- **结论：** Eliashberg 求解器在 λ>0.1 区间**系统性失败**，数据集不提供承重的 `Tc(μ*=0)` 覆盖。

---

## § 2  阻断分析

### 2.1  为何 Eliashberg 数据缺失阻断定理 2 检验

**定理 2（λ 分辨天花板，μ*=0）：** `𝔣 := Tc/(C_∞√λ ω̄₂) ≤ 𝔣_E(λ)`，其中 `𝔣_E(λ)` 是 Einstein 谱给出的亏损因子上界（λ=0.5→0.305，λ=1.0→0.627，λ=10→0.965）。

**检验需要：**
1. **分子 Tc(μ*=0)** — 各向同性 Eliashberg 全谱数值解，无拟合
2. **分母 C_∞√λ ω̄₂** — McMillan.dat 表头给出，解析

**当前状态：**
- McMillan.dat 的 Tc_McM/Tc_AD 列**不是** Eliashberg 解，它们是 McMillan/Allen–Dynes **半经验公式**，本质是 λ、ω_log、ω̄₂ 的确定性函数（见 § 2.2）。
- 用它们「检验」天花板 = 用 `f(λ,r)` 检验 `f(λ,r) ≤ C`，恒真或恒假，零信息。
- **Eliashberg 列才是独立数值解**，但 batch-q+p 共 749 材料中只有 2 个有 μ*=0 的值。

### 2.2  McMillan/AD 公式为何是循环量

Allen–Dynes (1975) 公式（简化版）：
```
Tc_AD = (ω_log/1.20) exp[−1.04(1+λ) / (λ − μ*(1+0.62λ))]
```
输入：λ、ω_log、μ*；输出：Tc。  
⟹ `Tc_AD(λ, ω_log, μ*)` 与 `C_∞√λ ω̄₂` 的比值 **= F(λ, r=ω̄₂/ω_log, μ*)**，是三个已知量的闭式函数。

**推论：** 把 AD 的 Tc 代入天花板或定理 2，检验的不是「数据是否满足理论」，而是「AD 公式是否满足理论」。  
- AD 公式本身**不**满足天花板（它在某些 (λ,r) 区间会超 Einstein 极值），这已在仓库笔记中记录。
- 用它做检验会混淆「公式缺陷」与「物理违反」，不可取。

### 2.3  剩余可行路线

**路线 A（本任务原计划，现不可行）：** 用 Cerqueira Eliashberg 数据做定理 2 的 λ 分辨检验。  
⟹ **前提为假，阻断。**

**路线 B（退化为天花板自检，零新信息）：** 用 McMillan/AD 列做管线负对照（乘 1.5 应报 FAIL），确认解析脚本正确。  
⟹ 可做，但不产出物理结论。

**路线 C（唯一非退化路线）：** SuperCon 交叉匹配——把 Cerqueira 的 (λ_DFT, ω̄₂) 与 SuperCon2/3DSC 的实验 Tc 交叉，得到「非隧道 (λ, ω̄₂, Tc_expt) 家族」。  
- **为何「非隧道」：** AD75 隧道样本（29 条）已在 `ad75_family_independence_and_lambda_curvature_2026-09-20.md` 中用于趋势回归，Cerqueira×SuperCon 交叉的材料**必不与 AD75 重叠**（时间跨度、测量手段）。
- **可匹配数量级估计：** Cerqueira 6912 材料 × SuperCon2 约 3 万条记录，化学式规范化后预期 100–500 条匹配（见 § 3）。
- **用途：** 推广点 1 的外部效度检验——在**非隧道**家族上，天花板与 λ 曲率是否成立。

---

## § 3  SuperCon 交叉匹配设计（仅设计，不执行）

### 3.1  匹配流程

1. **化学式规范化：**
   - Cerqueira：`material_dir` 给 `<formula>_<mcid>`，如 `Al2CoNi_agm003191160` → `Al2CoNi`
   - SuperCon2：化学式字段（如 `La2-xSrxCuO4`）需去除掺杂标记 `x`、空格、下标格式统一
   - 规范化目标：`AlNbRu2` vs `AlRu2Nb` → 字典序排序元素符号，统一为 `AlNbRu2`

2. **匹配策略：**
   - **一级匹配：** 规范化化学式**完全相等**
   - **二级匹配（如需要）：** 元素组成相同但化学计量数不同（如 `MgB2` vs `Mg2B4`）→ 归一化计量数到最小整数比
   - **排除：** 多相材料（`/`）、压力依赖记录（同一化学式多个 Tc）需人工裁决或取最高 Tc

3. **歧义处理：**
   - **多个 SuperCon 记录对应一个 Cerqueira 材料：** 报告所有匹配，附带 Tc 范围与测量方法
   - **多个 Cerqueira 材料对应一个 SuperCon 记录：** 报告 (λ, ω̄₂) 散布（DFT 计算的结构依赖）

4. **匹配失败模式：**
   - Cerqueira 材料不在 SuperCon 中（预期 90%+，DFT 筛选出的候选多数未实验合成）
   - SuperCon 记录无 DFT 计算（预期 95%+）
   - 化学式无法规范化（特殊字符、非化学计量）

### 3.2  可匹配数量级估计

- **Cerqueira：** 6912 材料，化学式去重后约 6500（少量重复结构）
- **SuperCon2：** 约 30000 条记录，去重后约 15000 种材料
- **预期匹配：** 10%（Cerqueira 是 DFT 筛选候选，SuperCon 是实验已知超导体）→ 约 **600 条**
- **可用于天花板检验：** 筛除 Tc 未报告、λ 缺失、ω̄₂ 缺失 → 预期 **100–300 条**

### 3.3  交叉匹配的输出

**CSV 列：** formula_normalized, cerqueira_mcid, cerqueira_lambda, cerqueira_wlog_K, cerqueira_w2av_K, cerqueira_r, supercon_Tc_K, supercon_Tc_method, supercon_year, match_quality（exact/stoich_normalized/ambiguous）

**用途：**
- 计算 `Tc_expt / (C_∞√λ ω̄₂)` 分布，检验天花板违反率
- 按 λ 分箱，检验 `𝔣(λ)` 的单调性与 Einstein 极值性
- **与 AD75 隧道样本对比**：两个独立家族的趋势是否一致

---

## § 4  已完成工作

### 4.1  解析器 `tools/cerqueira_collect.py`
- **输入：** `<root_dir>`（包含材料子目录），`<batch_name>`
- **输出：** CSV（16 列）
- **验证：** batch-q 249 材料与参考 CSV 逐位一致（PASS）
- **诊断输出（stderr）：** 材料数、Eliashberg 有数材料数、μ* 网格范围

### 4.2  batch-q 与 batch-p 解析结果
- **batch-q：** 249 材料，0/249 Eliashberg 有数（负对照 PASS）
- **batch-p：** 500 材料，2/500 Eliashberg 有数（前提验证，§ 0.1）
- **已生成：** `batch-q-parsed.csv`（250 行含表头），`batch-p.csv`（服务器 `/tmp/`，506 行含表头+诊断）

### 4.3  脚本入库状态
- `tools/cerqueira_collect.py`：worktree 中已创建，待 commit
- 分析脚本 `_chk_cerqueira_shape_deficit_2026-09-21.py`：**未完成**（前提为假，无数据可分析）

---

## § 5  交付物与后续

### 5.1  当前交付（本 worktree 分支）
1. **本笔记**：`theory/notes/materials/cerqueira6912_full_shape_deficit_2026-09-21.md`
2. **解析器**：`tools/cerqueira_collect.py`（已验证）
3. **验证脚本**：`verify_batchq_parse.py`（本地，不入库）
4. **结论：** Eliashberg 数据缺失阻断定理 2 检验；唯一可行路线是 SuperCon 交叉匹配（推广点 1）

### 5.2  未完成项（阻断，不派活）
- ~~分析脚本 `_chk_cerqueira_shape_deficit_2026-09-21.py`~~（无数据）
- ~~𝔣/𝔣_E(λ) 分布统计~~（无数据）
- ~~定理 2 违反扫描~~（无数据）
- ~~与氢化物谱形状亏损对比~~（无数据）

### 5.3  建议后续（需用户决策）
1. **立即：** 决定是否执行 SuperCon 交叉匹配（§ 3）——这是唯一非退化路线，但需额外工作量（SuperCon 数据清洗、化学式规范化、歧义裁决）
2. **备选：** 放弃 Cerqueira 数据集，回到氢化物 + AD75 隧道样本的已有结论
3. **低优先级：** 继续解析其余 15 批（batch-o, n, …, a）确认 Eliashberg 数据缺失是系统性的（预期结果相同）

---

## 附录 A：batch-p 的两个有数材料详情

### A.1  H6Ca2NiPd (agm004930255)
```
λ = 0.191
ω_log = 786.7 K (67.79 meV)
ω̄₂ = 1088.6 K (93.81 meV)
r = 1.384

Tc_Eliashberg(μ*=0.00) = 0.998 K
Tc_Eliashberg(μ*≥0.02) = ----- (全部求解失败)

Tc_McMillan(μ*=0.00) = 0.999 K
Tc_AD(μ*=0.00) = 1.008 K
```

### A.2  H6Sr2RhIr (agm004930587)
```
λ = 0.195
ω_log = 790.5 K (68.12 meV)
ω̄₂ = 1084.8 K (93.48 meV)
r = 1.372

Tc_Eliashberg(μ*=0.00) = 0.858 K
Tc_Eliashberg(μ*≥0.02) = ----- (全部求解失败)
```

**共同特征：**
- 都是含氢复杂化合物（H₆ 高氢含量）
- λ 都在 0.19 附近（batch-p 的 p90 分位）
- Eliashberg 求解器**只在 μ*=0.00 成功**，μ*>0 时全部失败
- 两个样本不足以做任何统计

---

## 附录 B：μ* 网格与 Eliashberg 求解失败的可能原因

### B.1  观察到的模式
- batch-q（λ ≤ 0.21，弱耦合）：0/249 Eliashberg 有数
- batch-p（λ ∈ [0.11, 0.21]）：2/500 Eliashberg 有数，且仅 μ*=0.00 一点
- 求解失败标记：`-----`（非数值错误，是求解器显式放弃）

### B.2  可能原因（推测，需查 Cerqueira et al. 2023 论文方法部分）
1. **收敛失败：** Eliashberg 迭代求解在某些 (λ, μ*, α²F 形状) 组合下不收敛
2. **计算成本筛选：** 作者可能只对「预期 Tc > 阈值」的材料跑 Eliashberg（弱耦合预期 Tc 低，跳过）
3. **μ*>0 的库仑项数值不稳定：** 两个有数材料都只在 μ*=0.00 成功，μ*>0 时失败，提示库仑项积分可能有数值问题
4. **α²F 质量问题：** DFT 声子计算的 α²F 可能在某些材料上有虚频、负值、或其他病态特征，导致 Eliashberg 求解器拒绝

### B.3  对本任务的影响
- **无法通过「跑更多批」解决**——batch-p 已是高 λ 区间（中位 0.17），若这里只有 0.4% 有数，其余批预期更差
- **无法通过「改求解器」解决**——我们没有原始 α²F 数据与 Eliashberg 求解器代码
- **结论：** 数据集设计时未以「提供 μ*=0 的 Tc 分布用于天花板检验」为目标，无法用于定理 2 验证

---

**提交状态：** 本笔记与 `tools/cerqueira_collect.py` 在 worktree 分支 `agent-a251c94abbc4f0ae2` 中，待用户审阅后 commit。  
**索引更新：** 待补 `theory/notes/INDEX.md` 行。
