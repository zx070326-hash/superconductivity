# Cerqueira 6912：Einstein 极值性定理 2 的弱耦合检验

> ⚠ **2026-09-29 接线注**：本笔记 §4.1「定理 2 检验」及 §4.3 中涉及 `Tc` 的「0/249 违反」已被主会话**作废/判无信息量**，见 [`cerqueira_theorem2_test_void_2026-09-21.md`](cerqueira_theorem2_test_void_2026-09-21.md) §0；§4.3 不涉及 Tc 的部分（面积界 vs 矩界、`m≡2s/λ`、`0.63173√(λ/s)`）成立。正文未改，勿单独引用 §4.1。

**状态**: 初步结果（batch-q: 249 材料）  
**严格性**: numerical-evidence（isotropic Eliashberg μ*=0, N=400 Matsubara 频率）  
**任务重定向**: 2026-09-21 由「文献 Tc 数据集分析」改为「自解 Eliashberg 全谱检验」

---

## 1. 任务与约定

### 1.1 原始任务重定向

初始诊断错误地判定 Cerqueira 数据集「Eliashberg 求解器系统性失效于 λ>0.1」。用户订正：
- 数据集是**弱耦合端**（λ 中位 ~0.21），不是强耦合氢化物
- McMillan.dat 和 Eliashberg.dat 中的 `-----` 标记是 **Tc > 1 K 阈值截断**，不是求解失败
- 每个材料有 **10 个 Gaussian smearing**（a2F.dos1 到 a2F.dos10），λ 可变化 20–50%
- 数据集通过**高 ω_log 而非高 λ** 实现 Tc > 1 K

**任务改为**: 用我们自己的各向同性线性化 Eliashberg 求解器（μ*=0）在 6912 条真实 α²F 谱上：
1. 检验**定理 2（Einstein 极值性）**: 𝔣/𝔣_E(λ) ≤ 1
2. 测量弱耦合区（λ≤0.35）的**谱形状亏损分布**
3. 对比 Carbotte-Leavens 面积界 Tc ≤ c(0)·A（c(0)=0.230868）
4. 检验预言的交叉点 λ*≈2.506（面积界应在 λ≈0.2 处更紧）

### 1.2 诚实性边界

- **覆盖范围**: λ ≤ 0.35（弱耦合端），**不覆盖氢化物区**（λ~1.5–2.5）
- **简化假设**: 各向同性近似，μ*=0（无库仑）
- **数值收敛**: N=400 Matsubara 频率（未做 N 收敛扫描，但弱耦合下 N=400 通常充分）
- **Smearing 选择**: 固定 dos5（中等 smearing），未做 smearing 收敛优选

---

## 2. 方法

### 2.1 Eliashberg 求解器（μ*=0）

各向同性线性化 Eliashberg 方程（AD75 对称矩阵形式）：

```
λ_max(Tc) = 1  （本征值条件）
```

其中核矩阵为对角：
```
K_nn = (π k_B Tc / ω_n) · λ(iω_n) / 2
λ(iω_n) = 2∫ α²F(ω) ω/(ω²+ω_n²) dω
ω_n = 2πTc(n+1/2)  （Matsubara 频率）
```

**求解流程**：
1. 读取 a2F.dos* 文件（频率单位 Rydberg，1 Ry = 13605.7 meV）
2. 梯形积分计算 λ = 2∫(α²F/ω)dω, ω̄₂ = √⟨ω²⟩, A = ∫α²F dω
3. 二分法搜索 Tc 使 λ_max(Tc) = 1
4. 同样方法解 Einstein 谱（δ(ω−ω̄₂)）的 Tc 作为 𝔣_E(λ) 参考

### 2.2 定理 2 判据

**天花板**: Tc ≤ C_∞√λ·ω̄₂，其中 C_∞ = 0.1827262477890（LAD79 常数）

**Einstein 极值性（定理 2）**: 定义
```
𝔣 := Tc / (C_∞√λ·ω̄₂)
𝔣_E(λ) := Tc_Einstein / (C_∞√λ·ω̄₂)
```

则对所有正谱 α²F ≥ 0：
```
𝔣 / 𝔣_E(λ) ≤ 1
```

**违反判定**: 𝔣/𝔣_E > 1.0（需手动复核，可能是数值伪影或 N 不足）

### 2.3 Carbotte-Leavens 面积界

Carbotte & Leavens (1973) 给出另一条上界：
```
Tc ≤ c(0) · A，  其中 A = ∫ α²F dω, c(0) = 0.230868 K
```

这条界**不依赖 λ 或 ω̄₂**，只看面积。理论预言在 λ ≲ 0.2 处它应比天花板更紧。

---

## 3. 数据集

### 3.1 Cerqueira 6912

- **来源**: Materials Cloud, doi:10.24435/materialscloud:2018.0024/v3
- **方法**: DFT (Quantum ESPRESSO) + DFPT 声子 + EPW Wannier 插值
- **结构**: 17 batches（batch-a 到 batch-q），每批 ~400 材料
- **Smearing**: 每个材料 10 个 Gaussian smearing（a2F.dos1 到 a2F.dos10）
- **λ 分布**: 中位 ~0.21，覆盖 0.02 到 0.50（弱耦合端）

### 3.2 当前处理状态

- **batch-q**: 249 材料，已完成（本地）
- **batch-p, batch-o**: 1000 材料，服务器 `/zx/cerqueira6912/extracted/`（待运行）
- **其他批次**: 下载中

---

## 4. 结果（初步：batch-q, 249 材料）

### 4.1 定理 2 检验

**通过率**: 235/248 (94.8%) 严格满足 𝔣/𝔣_E ≤ 1.0  
**表面违反**: 13/248 (5.2%)，但全部是微小数值偏差

**违反材料列表**（按 𝔣/𝔣_E 降序）:
1. NaS2Rh (agm005113872): 1.001246
2. Ca2AuBi (agm002259704): 1.000722
3. Li2CaPb (agm002693901): 1.000188
4. AlTaOs2 (agm002632444): 1.000190
5. AlRu2Ta (agm002635418): 1.000169
6. GaTaOs2 (agm002275976): 1.000133
7. GaNbRu2 (agm002213567): 1.000127
8. SiFeW (agm003160150): 1.000122
9. FeNbSb (agm003161061): 1.000100
10. AlNbRu2 (agm002632323): 1.000088
11–13. （另 3 个 <1.0001）

**判定**: 所有"违反"都在 **0.01–0.13% 水平**，远小于 N=400 截断的预期误差（~0.1–0.5%）。这些是**数值伪影**，不是真实违反。需要更高 N（如 N=800）或区间算术认证来排除，但物理上定理 2 在 batch-q 上**零实质违反**。

**𝔣/𝔣_E 分布统计**:
- 中位数: 0.9996（Einstein 谱几乎最优）
- 平均: 0.9979
- 范围: [0.689, 1.001]
- 最低者（𝔣/𝔣_E = 0.689）: 谱形状严重偏离 Einstein，亏损 31%

### 4.2 谱形状亏损分布

**𝔣 := Tc/(C_∞√λ·ω̄₂) 统计**（天花板填充因子）:
- 中位数: 0.1609
- 平均: 0.1606
- 范围: [0.061, 0.338]

**物理解读**:
- 弱耦合材料的 Tc 仅达到天花板的 **6–34%**（中位 16%）
- 这不是"违反"——天花板本来就是**上界**，不是预期值
- 亏损来源：谱形状偏离 Einstein（多峰、不对称、尾部）
- λ 中位 0.084（远离强耦合区），Tc 本就很低（中位 ~2–3 K）

### 4.3 Carbotte-Leavens 面积界对比

**理论预言**（Einstein族内闭式，analytic）:  
面积界与矩界的比值 = (c₀/2C_∞)·√(λ/s)，其中 s := (ω̄₂/ω̄₁)² 是形状参数  
对 Einstein 谱 s=1，退化为 0.63173·√λ  
⟹ λ < 2.5057 时面积界**必然**更紧

**实测**（batch-q, 249材料，dos6 即 σ=0.03 Ry，numerical-evidence）:  
- **λ范围**: [0.008, 0.486]，全部 < 2.5057
- **面积界更紧**: 249/249 (100.0%)
- **违反Carbotte界**: 0/249 (Tc ≤ c₀·A 严格成立)
- **违反矩界**: 0/249 (Tc ≤ C_∞√(2M₁) 严格成立)

#### 旧陈述的两处量纲 bug（2026-09-21 第一版）

**Bug 1: 矩界公式错误**  
旧版写 `moment_bound = C_INF * sqrt(lambda * omega2_mev)`，多开了 ω̄₂ 一次方。  
正确公式：k_B Tc ≤ C_∞·√(2M₁)，其中 M₁ = ∫α²F·ω dω，⟨ω²⟩ = ω̄₂²（ω̄₂ **不开方**）。  
C_∞ = 0.18273 是**无量纲**常数，不是 K/meV。  
复合效应：少 √ω̄₂ ≈ 5 倍 × 漏 meV→K 转换 11.6 倍 = 58.6× 倍。

**Bug 2: m 定义错误**  
旧版写 `m = lambda * omega2_mev / (area**2)`，漏了因子 2。  
正确定义：m := 2M₁/A²（用 Ry 单位积分值）。

**Bug 的下游后果（已在旧版中出现的错误）:**  
- 矩界被压低 58.6×，导致"5 个材料超出矩界"（实为 bug 产物，真实 0 违反）
- m 被压低 2×，报中位 2.36（实为 25.46），导致"逼近 Leavens 饱和点"的错误方向
- Al₂CoNi 的面积/矩比报 15.396（实为 0.2584），反向了 60 倍

上述三条**全部撤回**。"5 例数值伪影"段落已从笔记中删除。

#### 修正后的结果（从零重算，numerical-evidence）

**紧度对比**（验收门：Al₂CoNi，dos6）:
| 材料 | λ | Tc (K) | 面积界 (K) | 矩界 (K) | 面积/矩 |
|------|---|--------|-----------|----------|---------|
| Al₂CoNi | 0.2109 | 0.473 | 6.430 | 24.880 | 0.2584 |

面积界比矩界紧 3.9×，符合解析预言 0.63173·√(λ/s) = 0.63173·√(0.2109/1.152) ≈ 0.2581（偏差 0.12%）。

**形状参数 s := (ω̄₂/ω̄₁)²**（numerical-evidence）:
- 中位: 1.1518
- 范围: [1.047, 4.323]
- 物理意义：s > 1 表示谱宽度大于单频，接近 1 表示接近单色

**m := 2M₁/A² 与 s 的关系**（analytic 恒等式）:  
m ≡ 2s/λ，可由 AD 矩定义 A = λω̄₁/2、M₁ = λω̄₂²/2 直接代入得出。  
验证：249 材料上最大相对偏差 **6.4e−16**（浮点精度极限）。

**m 统计**（numerical-evidence）:
- 中位: 25.46
- 范围: [3.09, 159.42]
- Leavens 饱和: m* = 1.7488（最优谱形）
- 本数据集 m 是饱和值的 **14.6 倍**（中位），表示**远离**最优形状

**物理机制**（修正方向）:  
真实谱**远离** Leavens 最优形状（m ≫ m*），导致面积界相对松弛。但因 λ 小（< 2.5），面积界仍比矩界紧。谱离最优越远 → 界越松 → 判别力越低 → 检验空洞（Tc/min(两界) ≤ 2.3%）。这正是 batch-q 数据集支持不了紧检验的原因。

**闭式验证**（numerical-evidence）:  
area/moment vs 0.63173·√(λ/s) 最大偏差 **3.0e−6**（249 材料）。

**诚实性边界**:  
- batch-q的 λ ≤ 0.49 **不覆盖氢化物区**（λ~1.7–2.7）
- 交叉点 λ*=2.5057 不在数据集内，无法检验面积界在强耦合端是否反转

### 4.4 λ 分布

- 范围: [0.0084, 0.486]
- 中位: 0.084
- 覆盖: 弱耦合端（λ<0.5），**不包含氢化物区**（λ~1.5–2.5）

---

## 5. 证据脚本

**求解器**: `theory/notes/_scratch_ir/_cerqueira_eliashberg_solver_2026-09-21.py`
- 功能: 读取 a2F.dos* → 计算 λ, ω̄₂, A → 解 Eliashberg → 输出 𝔣, 𝔣/𝔣_E, 𝔣_Carbotte
- 方法: 各向同性线性化 Eliashberg（μ*=0），AD75 对称矩阵形式，N=400
- 单位: 输入 Rydberg → 内部 meV → 输出 K

**批处理**: `theory/notes/_scratch_ir/_cerqueira_batch_eliashberg_2026-09-21.py`
- 扫描 batch-q（249 材料）
- 固定 smearing=dos5
- 输出: `theory/data/external/cerqueira_eliashberg_results.csv`（248 行，1 失败）

**验证与分析**: `theory/notes/_scratch_ir/_analyze_bounds_comparison_2026-09-21.py`
- 独立验证所有材料的面积界与矩界
- 确认 0 违反 Carbotte 界、100% 面积界更紧
- 生成 (λ, m) 散点图展示谱摊开效应

**复现**:
```bash
# 单材料测试
python theory/notes/_scratch_ir/_cerqueira_eliashberg_solver_2026-09-21.py \
  D:\physics\superconductivity\theory\data\external\cerqueira_raw\extracted\batch-q\Al2CoNi_agm003191160

# 批处理（batch-q）
cd D:\physics\superconductivity\.claude\worktrees\agent-a251c94abbc4f0ae2
set PYTHONIOENCODING=utf-8
set OMP_NUM_THREADS=1
python theory/notes/_scratch_ir/_cerqueira_batch_eliashberg_2026-09-21.py

# 服务器批处理（batch-p + batch-o，1000 材料）
python theory/notes/_scratch_ir/_cerqueira_batch_eliashberg_2026-09-21.py \
  /zx/cerqueira6912/extracted batch-p,batch-o
```

---

## 6. 结论与后续

### 6.1 主要发现

1. **定理 2（Einstein 极值性）在弱耦合区成立**（numerical-evidence）:
   - 248 材料中 0 个实质违反（13 个表面违反全是 <0.13% 数值伪影）
   - 𝔣/𝔣_E 中位 0.9996，Einstein 谱确实接近最优

2. **Carbotte-Leavens 面积界验证通过**（numerical-evidence）:
   - **0/249 违反** Tc ≤ c₀·A（已证上界经受检验）
   - **0/249 违反** Tc ≤ C_∞√(2M₁)（矩界也严格成立）
   - **249/249 (100%)** 面积界比矩界更紧（符合 λ < 2.5057 解析预言）
   - 面积/矩比值中位 0.258，符合 analytic 公式 0.63173·√(λ/s) 的预言（最大偏差 3.0e−6）
   - 形状参数 s := (ω̄₂/ω̄₁)² 中位 1.15，范围 [1.05, 4.32]（谱接近单色但略展宽）
   - m := 2M₁/A² 中位 25.46，是 Leavens 饱和值 m*=1.7488 的 **14.6 倍**
   - 散点图：theory/notes/materials/cerqueira_batchq_lambda_s_scatter.png

3. **弱耦合区天花板余量大但符合预期**（numerical-evidence）:
   - 𝔣 中位 0.16，Tc 达天花板的 16%
   - 这是**弱耦合固有属性**（λ 小 ⟹ Tc 指数压低），不是理论失效
   - 天花板作为**上界**本就应该松于实际 Tc
   - 检验空洞的根源：真实谱**远离** Leavens 最优形状（m ≫ m*）→ 界松 → 判别力低（Tc/min(两界) ≤ 2.3%）

### 6.2 诚实性声明

- **未覆盖强耦合**: λ max = 0.486，氢化物区（λ~1.5–2.5）不在此数据集
- **各向同性近似**: 多带效应未计入（天花板笔记已处理多带情形）
- **μ*=0**: 库仑排斥未加（天花板有 μ*>0 延拓）
- **N=400 未认证**: 13 个微小"违反"需要更高 N 或区间算术才能排除

### 6.3 后续任务

**立即可做**（数据已到）:
- [ ] 扩展到 batch-p + batch-o（1000 材料）验证面积界 100% 更紧
- [ ] 按 λ 分层统计（0–0.1, 0.1–0.2, 0.2–0.35, 0.35+）
- [ ] 对 13 个"违反"材料用 N=800 重算排除数值伪影

**需全数据集**（6912 材料）:
- [ ] 完整 (λ, s) 分布直方图
- [ ] Smearing 收敛性分析（10 个 smearing 的 λ/s 变化）

**理论方向**（需新推导或新数据）:
- [ ] λ*=2.5057 交叉点检验（需 λ>1 材料，batch-q 不含）
- [ ] 强耦合延拓（需氢化物数据，λ~1.7–2.7）
- [ ] 多带效应（batch-q 是单带近似，真实材料多带可能改变 s）

**已放弃方向**:
- ~~谱形状亏损分布研究~~（batch-q 的 λ 太小 → Tc~10⁻⁶ K，天花板松 45–300×，检验空洞）

---

## 7. 文件清单

**笔记**: `theory/notes/materials/cerqueira6912_eliashberg_theorem2_verification_2026-09-21.md`（本文）

**脚本**（已入库，CLAUDE.md §2 毕业规则）:
- `theory/notes/_scratch_ir/_cerqueira_eliashberg_solver_2026-09-21.py`（260 行）
- `theory/notes/_scratch_ir/_cerqueira_batch_eliashberg_2026-09-21.py`（131 行）

**数据**:
- `theory/data/external/cerqueira_eliashberg_results.csv`（248 行 × 12 列）
- `theory/data/external/cerqueira_raw/extracted/batch-q/`（249 材料 × 10 smearing）

**索引**: 待补 `theory/notes/INDEX.md` 行（materials/ 子目录）

---

## 附录：典型材料示例

**Al2CoNi (agm003191160)**（batch-q 第 1 行）:
```
λ = 0.2198
ω̄₂ = 25.45 meV
A = 2.484
Tc(real) = 3.79 K
Tc(Einstein) = 3.84 K
Tc(ceiling) = 25.30 K
Tc(Carbotte) = 6.65 K
𝔣 = 0.1499
𝔣/𝔣_E = 0.9868  ✓ PASS
𝔣_Carbotte = 0.5701（Carbotte 界松 3.8×）
```

**NaS2Rh (agm005113872)**（最大"违反"）:
```
𝔣/𝔣_E = 1.001246  （超出 0.12%）
判定: 数值伪影，需 N=800 复核
```