# D3S Tc=150K 压强配对核查 (2026-09-20)

> **任务**：核查 `theory/data/hydrides_theory.csv`:9 的 D3S Tc=150K 是否配 200 GPa。
> **驱动**：team-lead 在 `h3s_on_ad75_trend_2026-09-20.md` §11 中指出，H3S 已栽过压强错配（203K 实为 155 GPa，200 GPa 处是 184.04K），D3S 的 150K 极可能同形。

---

## §1  CSV 记录与引用链

`theory/data/hydrides_theory.csv`:9（D3S-200GPa anharmonic）：

| 字段 | 值 | 备注 |
|---|---|---|
| `material` | D3S | |
| `pressure_GPa` | 200 | **声明压强** |
| `lambda` | 1.87 | 来自 Errea PRL 2015 Table II |
| `omega_log_meV` | 73.3 | 来自 Errea PRL 2015 Table II |
| `Tc_predicted_K` | 152.0 | Errea 自己用 ME 解的 |
| notes 列 | "Tc(Exp)=90 K" | **引用 Errea Table II 最后一列** |

**问题 1**：notes 列的 `Tc(Exp)=90 K` 与 CSV 其他地方使用的 150 K 不一致。  
**问题 2**：Errea Table II 的 90 K 已被 team-lead 明确标记为「出自被撤回的未退火预印本，不可用」。

---

## §2  Errea PRL 2015 Table II 原文核验

**来源**：`theory/notes/_scratch_ir/errea2015.pdf` p.3, Table II（页图 + pdftotext 双路核验）

| Compound | λ^anh | ω^anh_log | T^ME,anh_c | **T_c(expt.)** |
|---|---|---|---|---|
| D₃S (200 GPa) | 1.87 | 73.3 | 152.0 | **90** |

**确认**：Errea Table II 确实写 `Tc(expt.) = 90 K`。

**但 Errea 自己在正文中警告**（未在 Table 里标注）：这个 90 K 出自 arXiv:1412.0460（Drozdov 预印本），  
且 Drozdov Nature 2015 正式版 Fig. 2(b) caption 明确说：  
> "Note here that **the only point for sulfur deuteride presented in arXiv:1412.0460 was determined without sample annealing**, and Tc would increase after annealing at room temperature."

**结论 1**：Errea Table II 的 `Tc(expt.)=90 K` 是**未退火值**，被 Drozdov 正式版点名撤回，不可用作锚点。

**严格性**：`事实核验（Drozdov Nature 2015 Fig. 2(b) caption, 页图原文）`

---

## §3  Drozdov Nature 2015 正式版 Fig. 2(c) 数字化

**来源**：`theory/notes/_scratch_ir/drozdov2015.pdf` p.9, Fig. 2(c)

Fig. 2(c) caption 明确说：
> "Open color points refer to **sulfur deuteride** (D3S), and filled points to sulfur hydride (H3S).  
> Data on **annealed samples** are presented."

### 本次目测读数（承重前需精确数字化）

从 Fig. 2(c) 目测，D3S（空心圆）在 200 GPa 附近：
- P ≈ 195 GPa → Tc ≈ 138 K
- P ≈ 205 GPa → Tc ≈ 140 K

**线性插值**：200 GPa → Tc ≈ **139 K**

**严格性**：`目测估计（±5 K 不确定度）`，承重前需矢量数字化。

---

### 主会话已完成的精确数字化（2026-09-18）

**来源**：`theory/notes/materials/h3s_tc_exp_pressure_audit_2026-09-18.md`

主会话用 `fitz` 解析 Drozdov Fig. 2c 内嵌矢量坐标，得出：

> **"D3S（红色空心圆，annealed）在 P≈200 GPa 处：200.0 GPa → 145/146 K（两个相邻标记）"**

并明确说：
> "这与仓库 `hydrides_experimental.csv` 记的 ~150 K **一致**（差 3%），  
> 与 Errea Table II 的 `T_c(expt.)=90 K` **不一致**（差 60 K，62%）——  
> **独立确认了 team-lead 的判断：Errea 那一列的 90 K 就是被 Drozdov 正式版点名撤回的未退火赝值，  
> `hydrides_experimental.csv` 的 ~150 K 才是干净值**。"

**严格性**：`矢量数字化（主会话已完成，fitz 解析 PDF 内嵌坐标）`

---

## §4  仓库 150 K 的出处追溯

通过 grep 搜索，发现：
1. `deficit_budget_h3s_2026-09-17.md` 使用 D3S Tc=150 K 作为实验锚点
2. `h3s_tc_exp_pressure_audit_2026-09-18.md` 已追溯到 Drozdov Fig. 2c，并矢量数字化出 **145-146 K @200 GPa**

**结论 2**：仓库的 150 K **不是**来自 Errea Table II 的 90 K（那个是未退火赝值），  
而是来自 **Drozdov Nature 2015 Fig. 2(c) 的 annealed 数据**。

主会话的数字化给出 **145-146 K @200 GPa**，与 150 K 的差异为：
- 150 K vs 145.5 K（取中值）= **+4.5 K = +3.1%**

---

## §5  ~~与 Mishra α²F 的压强匹配性~~ → **订正（2026-09-20 晚）：CSV 第 9 行的 ω̄₂ 缺失，不是混用**

> ⚠️ **本节 §5–§6 判决 B 的「混用发现」已撤回**（2026-09-20 team-lead 指正）。  
> **旧陈述错在哪**：我把 team-lead 消息里给的 Mishra ω̄₂=105.43 meV 当成了 CSV 的内容，犯了「转引当读原件」的错误。  
> **实际情况**（已读 CSV 原件核验）：`hydrides_theory.csv`:9 的 `omega_2_meV` 列是 **`not-reported`**，不存在混用。  
> **新陈述凭什么对**：CSV 第 9 行**不能喂 AD75 回归**的真实理由是 **ω̄₂ 缺失**（`not-reported`），不是「来源不一致」。这两件事的修法完全不同：前者只要补一个同源的 ω̄₂，后者要推翻整行。

---

### Mishra vs Errea 的 ω_log 差异（仍然成立，但约束的是**未来补数动作**）

team-lead 在消息中给出（2026-09-20）：
> "直接积分给 `ω_log(D3S) = 85.535 meV`（我自写解析器+自写求积）。"

**对比**：
- Mishra α²F → ω_log = **85.535 meV**（team-lead 计算）
- Errea Table II → ω_log = **73.3 meV**（差 16.7%）

这与 H3S 的情况**同形**：Mishra 与 Errea 是**两个独立 SSCHA 计算**，给出不同的 α²F（H3S 侧差 18.5%，D3S 侧差 16.7%）。

---

### AD75 回归的三元组约束（约束未来补数，不是描述现存缺陷）

team-lead 在 H3S §11 中已确认：
> "把 Errea 的 `ω_log=92.86` 代进一个 `λ`、`ω̄₂` 都来自 Mishra α²F 的三元组，  
> **不是合法的参数选择**——AD75 回归要求三个量出自同一条谱。"

**同理适用于 D3S**（但这是**未来补数时的约束**，不是当前状态的缺陷描述）：
- 若未来补 Mishra 路径：需 (λ, ω_log=85.535, ω̄₂) 全部从 Mishra α²F 重算
- 若未来补 Errea 路径：需 Errea 自己的完整谱（但 Errea Table II 无 ω̄₂ 列 ⟹ 该路径阻塞）

**CSV 第 9 行当前状态**（订正后）：
- `lambda=1.87, omega_log=73.3`（Errea Table II）
- `omega_2_meV = not-reported`（**缺失**，不是混用）
- ⟹ **不能喂进 AD75 回归**，因为缺第三个量

---

## §6  200 GPa 压强配对的最终判决

### 判决 A：Tc_exp 数值的压强匹配性 ✅ **基本匹配**

| 来源 | 压强 | Tc_exp | 与 200 GPa 的偏差 | 限定 |
|---|---|---|---|---|
| Drozdov Fig. 2(c) 数字化（主会话） | 200 GPa | **145-146 K** | 0 GPa（精确） | **单路**矢量数字化（:114 限定） |
| 仓库 CSV / 笔记普遍使用 | **未明确标注** | 150 K | +4.5 K = +3.1% | |

**结论 A**：仓库的 150 K 与 200 GPa **基本配对**，偏差 +3.1% 在实验不确定度内。  
但**应改用主会话数字化的 145.5 K**（取 145-146 K 中值）以消除该 +3.1% 的系统偏移。

**限定**（team-lead 2026-09-20 指正）：主会话 `h3s_tc_exp_pressure_audit_2026-09-18.md:65` 的 145/146 K 是 Fig. 2c 顺带产出，  
该轮只对 **H3S 的 184.04 K 做了双路复核**（:73 第二次独立数字化），**D3S 的 145/146 K 目前是单路**，不是双路。

**严格性**：`矢量数字化（主会话完成，单路，承重前需第二路复核）`

---

### 判决 B：~~α²F 的来源一致性~~ → **ω̄₂ 缺失，不可评估** ❌

~~CSV 第 9 行当前状态：~~  
~~- `lambda=1.87, omega_log=73.3`（Errea Table II）~~  
~~- `omega_2=105.43`（应为 Mishra，需确认）~~  
~~- ⟹ **混用了两个独立计算的中间量**~~

**订正后的判决 B**（2026-09-20 晚）：

CSV 第 9 行当前状态（已读原件核验）：
- `lambda=1.87, omega_log_meV=73.3`（Errea Table II）
- `omega_2_meV = not-reported`（**缺失**）
- ⟹ **不能用于 AD75 回归**，因为 AD75 回归式需要三个量 (λ, ω_log, ω̄₂)

**未来补数的两条路径**（二选一）：

#### 路径 1：全用 Mishra α²F（推荐，因 Errea 路径阻塞）
- λ = ?（需从 `D3S.a2f` 重算）
- ω_log = 85.535 meV（team-lead 已算）
- ω̄₂ = ?（需从 `D3S.a2f` 重算，team-lead 消息中给的 105.43 meV 需确认来源）
- Tc_exp = 145.5 K（Drozdov Fig. 2c @200 GPa，单路数字化）

#### 路径 2：全用 Errea α²F（**阻塞**）
- λ = 1.87
- ω_log = 73.3 meV
- ω̄₂ = **未报**（Errea Table II 无此列）
- ⟹ **路径 2 阻塞**，无法评估

**判决 B（订正后）**：当前 CSV 第 9 行**不能用于 AD75 回归判决**，因为 **ω̄₂ 缺失**（`not-reported`）。  
未来补数时，必须整行改用 Mishra α²F（λ, ω_log, ω̄₂ 三个量同源），不能只补第三个量。

---

### ~~判决 C：余量计算的有效性~~ → **撤回，未提供消费点**

~~team-lead 报告：~~  
~~> | D3S，Tc=150 K | **−4.48%** | PASS | **+5.43 K** |~~

~~但该余量计算基于：~~  
~~- Tc_exp = 150 K（应改为 145.5 K）~~  
~~- 混用的 α²F 三元组（不合法）~~

~~**若改用 Tc_exp=145.5 K**（下调 4.5 K）：~~  
~~- 原残差 −4.48% 对应余量 +5.43 K~~  
~~- 下调 4.5 K 后，余量从 +5.43 K 降到 **+0.93 K**（几乎触及 FAIL 边界）~~

~~**判决 C**：下压 **4.5 K** 即从 PASS 变 FAIL，team-lead 的警告精确命中。~~

**订正（2026-09-20 team-lead 指正）**：  
我未给出该余量的**消费点**（哪张表、哪个不等式里算出 +5.43 K），  
且「混用 α²F」的前提已撤回 ⟹ **判决 C 整段撤回**。

按 CLAUDE.md：「目标常数必须随消费点一起交付」，  
若要报「接近 FAIL」，需补上：该余量在哪个判据、哪一行里被消费。  
当前未提供 ⟹ 该判决不成立。

---

## §7  完整证据链汇总

### 证据 1：Errea Table II 的 90 K 不可用 ✅

- **来源**：Errea PRL 114, 157004 (2015) Table II, p.3
- **事实**：`Tc(expt.) = 90 K`
- **撤回依据**：Drozdov Nature 2015 Fig. 2(b) caption 明确说该值出自**未退火**预印本
- **文件**：`theory/notes/_scratch_ir/errea2015.pdf` p.3, `drozdov2015.pdf` p.9
- **严格性**：`事实核验（页图 + 文本层双路）`

---

### 证据 2：Drozdov Fig. 2(c) 给出 D3S @200 GPa = 145-146 K ✅

- **来源**：Drozdov Nature 2015 Fig. 2(c), p.9
- **数字化**：主会话 `h3s_tc_exp_pressure_audit_2026-09-18.md` 用 `fitz` 矢量解析
- **结果**：200.0 GPa → 145/146 K（annealed）
- **文件**：`theory/notes/_scratch_ir/drozdov2015.pdf` p.9
- **严格性**：`矢量数字化（主会话完成）`

---

### 证据 3：仓库 150 K 与 145.5 K 差 +3.1% ✅

- **计算**：(150 - 145.5) / 145.5 = +3.1%
- **解释**：150 K 可能是早期粗读或四舍五入；主会话数字化更精确
- **严格性**：`算术（精确）`

---

### 证据 4：Mishra α²F 的 ω_log = 85.535 meV（vs Errea 73.3 meV，差 16.7%） ✅

- **来源**：team-lead 消息（2026-09-20）
- **方法**：从 `D3S-200GPa/EPW_anharmonic/D3S.a2f` 直接积分
- **对比**：与 Errea Table II 的 73.3 meV 差 **16.7%**
- **结论**：两个独立 SSCHA 计算，与 H3S 的 18.5% 差异同形
- **严格性**：~~`numerical-evidence（team-lead 计算，本 agent 未独立复现）`~~
  → **已升级**：`numerical-evidence（主会话从零重算，三重交叉）`，见本文 §9 项 3 的结案块。

---

### ~~证据 5：CSV 第 9 行混用了两个 α²F 的中间量~~ **已撤回** ❌

> **主会话 2026-09-20 补：本条与本文 §8「发现 3」直接矛盾，且 §8 才是正确的那一版。**
> §8 已自行撤回「混用」说（CSV:9 的 `omega_2_meV` 是 `not-reported`，**缺失**而非混用），
> 但 §7 这条漏改，仍挂着 ✅。读者只扫 §7 会取到已撤回的陈述 ⟹ 此处补标撤回。
> **正确陈述见 §8「发现 3」。** 不删除原文（保留撤回痕迹）。

~~- **事实**：λ, ω_log 来自 Errea；ω̄₂ 应来自 Mishra（需确认）~~
~~- **违反**：AD75 回归要求三元组 (λ, ω_log, ω̄₂) 出自同一条谱~~
~~- **后果**：该行**不能用于 AD75 回归判决**~~
~~- **严格性**：`逻辑推理（基于 team-lead 在 H3S 案例中确立的规则）`~~

**存活的部分**：「该行不能用于 AD75 回归」这个**结论**仍然成立，但**理由换了**——
是 ω̄₂ **缺失**，不是来源混用。两者修法不同（补一个同源数 vs 推翻整行）。

---

## §8  给 team-lead 的交付（订正版，2026-09-20 晚）

### 一句话回答

**D3S Tc=150K 与 200 GPa 基本配对（主会话单路数字化给 145-146K，偏差 +3.1%），应改用 145.5K 以消除系统偏移；CSV 第 9 行因 ω̄₂ 缺失（`not-reported`）不能用于 AD75 回归；未来补数时必须整行改用 Mishra α²F（Mishra ω_log=85.535 vs Errea 73.3 差 16.7%，与 H3S 同形），不能只补第三个量。**

---

### 两段式订正：我的「严重性 A 混用发现」是假的

#### 旧陈述错在哪

我在初版交付中声称：
> "CSV 第 9 行混用了 Errea 的 (λ, ω_log) + Mishra 的 ω̄₂=105.43 meV，  
> 该行不能喂进 AD75 回归，必须整行改用同源 α²F。"

并把它列为「严重性 A，阻塞判决」。

**错误根源**：我把 team-lead 消息里给的 Mishra D3S ω̄₂=105.43 meV **当成了 CSV 的内容**，  
犯了「转引当读原件」的错误，违反了 CLAUDE.md「承重项自己从零重算」的要求。

**实际情况**（已读 `hydrides_theory.csv`:9 核验）：  
该行的 `omega_2_meV` 列是 **`not-reported`**，不存在混用。

---

#### 新陈述凭什么对

**CSV 第 9 行当前状态**（已读原件核验）：
- `lambda=1.87, omega_log_meV=73.3`（Errea Table II）
- `omega_2_meV = not-reported`（**缺失**，不是混用）
- ⟹ **不能用于 AD75 回归**的真实理由是 **ω̄₂ 缺失**，不是「来源不一致」

**两种缺陷的修法完全不同**：
- **ω̄₂ 缺失** → 只要补一个同源的 ω̄₂
- **来源混用** → 要推翻整行

---

#### 仍然成立的约束（改为「未来补数时的约束」，不是「现存缺陷」）

Mishra vs Errea 的 ω_log 差异（team-lead 2026-09-20 计算）：
- Mishra D3S: ω_log = **85.535 meV**
- Errea Table II: ω_log = **73.3 meV**
- **差 16.7%**（与 H3S 侧 18.5% 同形）

⟹ **任何人往 CSV 第 9 行补 ω̄₂ 时，必须整行改用 Mishra 的 (λ, ω_log, ω̄₂)，不能只补第三个量。**

Team-lead 认可这条顾虑，但要求改措辞：  
写成「现存缺陷」是假的，写成「补数时的约束」是真的。

---

### 详细发现（订正后，三项）

#### 发现 1：Tc_exp 数值 ✅ 基本匹配，但需修正（附限定）

| | 压强 | Tc_exp | 来源 | 限定 |
|---|---|---|---|---|
| 主会话数字化 | 200 GPa | **145-146 K** | Drozdov Fig. 2(c), annealed | **单路**（未双路复核） |
| 仓库 CSV / 笔记 | **未明确** | 150 K | 出处未追到一手 | 与数字化一致 ±3% |

**建议**：改用 **145.5 K**（取 145-146 K 中值）。

**限定**（team-lead 2026-09-20 指正）：  
主会话 `h3s_tc_exp_pressure_audit_2026-09-18.md:65` 的 145/146 K 是 Fig. 2c 顺带产出，  
该轮只对 **H3S 的 184.04 K 做了双路复核**（:73 第二次独立数字化），  
**D3S 的 145/146 K 目前是单路**，承重前需第二路复核。

~~**余量影响**：从 +5.43 K 降到 **+0.93 K**（下压 4.5 K 即 FAIL）。~~  
**余量影响**（撤回）：我未提供消费点（哪张表、哪个不等式），team-lead 不采纳该余量计算。

---

#### 发现 2：Errea Table II 的 90 K 是未退火赝值 ✅ 已知（非新发现）

- Drozdov Nature 2015 Fig. 2(b) caption 明确点名
- 与 annealed 值 145 K 差 **60 K（62%）**
- CSV 第 9 行 notes 列的 "Tc(Exp)=90 K" 应删除或标记"撤回"

**Team-lead 指正**：这条早已在仓库定案（`h3s_tc_exp_pressure_audit_2026-09-18.md` §D5 + memory），  
我的复核是**第三次确认**，不是新发现——不必再当作「发现 3」计入交付。

---

#### 发现 3：CSV 第 9 行 ω̄₂ 缺失 ❌ 不可评估（订正后）

~~CSV 第 9 行混用 Errea (λ, ω_log) + Mishra (ω̄₂)~~（**假发现，已撤回**）

**实际情况**（已读 CSV 原件核验）：
- `omega_2_meV = not-reported`（**缺失**）
- ⟹ **不能用于 AD75 回归**，因为回归式需要三个量 (λ, ω_log, ω̄₂)

**未来补数的约束**：
- Mishra ω_log=85.535 vs Errea 73.3 差 **16.7%**（与 H3S 的 18.5% 同形）
- ⟹ 补数时必须整行改用 Mishra α²F，不能只补 ω̄₂ 而保留 Errea 的 λ/ω_log

---

### 文件清单（绝对路径）

#### 本次新建

1. `D:\physics\superconductivity\theory\notes\materials\d3s_tc_pressure_pairing_resolved_2026-09-20.md`  
   （本文件，含订正）

2. `D:\physics\superconductivity\theory\notes\_scratch_ir\_digitize_drozdov_fig2c_d3s_2026-09-20.py`  
   （目测数字化脚本，承重前需精确方法）

#### 已有证据

3. `D:\physics\superconductivity\theory\notes\_scratch_ir\errea2015.pdf` p.3 Table II  
   （Errea 的 90 K，已核验）

4. `D:\physics\superconductivity\theory\notes\_scratch_ir\drozdov2015.pdf` p.9 Fig. 2(b)(c)  
   （annealed D3S 数据点 + 未退火警告）

5. `D:\physics\superconductivity\theory\notes\materials\h3s_tc_exp_pressure_audit_2026-09-18.md`  
   （主会话已完成的 145-146 K @200 GPa 单路数字化，:65）

6. `D:\physics\superconductivity\theory\data\hydrides_theory.csv`:9  
   （D3S 行，ω̄₂ 缺失）

---

### 建议行动（订正后，两项）

#### 优先级 1：修正 Tc_exp 为 145.5 K ⚠️

- **当前**：150 K（来源模糊，但与单路数字化一致 ±3%）
- **应改为**：145.5 K（主会话单路数字化，取 145-146 K 中值）
- **限定**：承重前需对 D3S 做第二路复核（主会话只对 H3S 做了双路）
- ~~**影响**：余量从 +5.43 K 降到 +0.93 K~~（撤回，未提供消费点）

---

#### 优先级 2：CSV :9 notes 列删除 "Tc(Exp)=90 K" 📝

- 该值已被 Drozdov 正式版点名撤回（未退火）
- 建议改为："Tc(Exp,annealed)=145.5 K @200 GPa (Drozdov Nature 2015 Fig. 2c, 主会话单路数字化，承重前需双路复核)"

---

~~#### 优先级 3：统一 D3S 的 α²F 来源~~（撤回，这不是「修现存缺陷」而是「补新数据」）

~~**当前状态**：CSV :9 混用 Errea 的 (λ, ω_log) + 疑似 Mishra 的 ω̄₂~~

**订正后**：  
CSV :9 的 ω̄₂ 是 `not-reported`（缺失），不存在混用。  
未来若要补数以评估 AD75 回归，需整行改用 Mishra α²F（因 Errea 路径阻塞：Table II 无 ω̄₂ 列）。

~~**我可以做第 1 项**（重算 Mishra D3S 的 λ），若 team-lead 需要。~~  
**Team-lead 说**：不要现在算，先完成订正。

---

## §9  遗留问题（不影响当前判决）

1. **D3S 的 145/146 K 需要第二路复核**：  
   主会话只对 H3S 做了双路矢量数字化，D3S 目前是单路（顺带产出）。

2. **仓库 150 K 的第一手出处**：  
   虽然与单路数字化一致 ±3%，但未追到「谁在哪篇笔记里第一次写下 150 K」。

3. ~~**未来补数若选 Mishra 路径，需算**~~ → **已结案，见下方结案块**：  
   - λ = 2∫(α²F/ω)dω（从 `D3S.a2f`）
   - ω̄₂ 来源确认（team-lead 给的 105.43 meV 需核对文件）

---

**交付完成时间**：2026-09-20（初版）/ 2026-09-20 晚（订正版）  
**执行 agent**：h3s-omegalog-convention (sonnet-5, teammate context)  
**订正原因**：team-lead 指正「严重性 A 混用发现」是转引错误  
**复核状态**：awaiting team-lead decision on Tc_exp correction & 单路数字化的承重性

---

## §10  §9 项 3 的结案块（主会话 2026-09-20，从零重算，不转引）

脚本：`../_scratch_ir/_chk_main_mishra_moments_2026-09-20.py`；
配套笔记：`mishra_moments_mainsession_verification_2026-09-20.md`。
直接读 `D3S-200GPa/EPW_anharmonic/D3S.a2f`，**三重交叉**（梯形 · Richardson ·
文件自报的 `Integrated el-ph coupling` 运行积分端点——一个免费的第三方锚点）：

| 量 | 值 | 交叉核对 |
|---|---|---|
| λ | **1.783761** | 文件自报 1.7838036，相对差 **2.4e−5** |
| ω_log | **85.535326 meV** | 两求积方案一致 ≤1.3e−7 |
| ω̄₂ | **105.430104 meV** | 同上 ⟹ **105.43 已核对文件，为真数** |

⟹ **路径 1（全用 Mishra α²F）所需的三元组已全部备齐、同源。**
配套 H3S 非谐同源三元组：λ=1.720958 / ω_log=109.995163 / ω̄₂=149.262922（文件自报 1.7209949）。

**两条限定随结论同行**：

1. 填 CSV:9 时必须**整行**替换成 Mishra 的 (λ, ω_log, ω̄₂) = (1.783761, 85.535326, 105.430104)，
   **不得**保留 Errea 的 λ=1.87 / ω_log=73.3 只补第三个量（本文 §8 已确立的约束；
   Mishra vs Errea 的 ω_log 差 16.7%，与 H3S 侧 18.5% 同形）。
2. **本次未改 CSV**——数据文件的改动等用户放行。本块只把「需算」的两项算完并备案。

**严格性**：`numerical-evidence`（float64 求积，两方案 + 文件自报锚点三重交叉）。

> **连带**：`mishra_moments_mainsession_verification_2026-09-20.md` 同时**推翻**了姊妹笔记
> `h3s_omega_log_discrepancy_resolved_2026-09-20.md` §3 的物理归因
> （「ω_log 偏离 1/√2 因为非谐」——谐波对照组已有 +4.99%/+9.97%，且其解析前提自相矛盾）。
> 真实机制是**部分替换**（S 原子不伸缩）。这条对本文的 Tc/压强配对判决**无影响**，
> 但凡引用「D/H 频率按 1/√2 均匀伸缩」的地方都受影响。
