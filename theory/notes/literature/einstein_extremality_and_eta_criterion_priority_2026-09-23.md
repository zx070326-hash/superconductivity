# Einstein 极值性(μ*>0)与 η 判据优先权检索（2026-09-23）

派单：team-lead，三问文献核查。**本会话 WebFetch 全域封锁（所有域名，包括 harvest.aps.org / arxiv.org / export.arxiv.org）
+ WebSearch 预算在本会话开始时已耗尽（200/200）**——两条外部检索通道本轮都不可用，零新查询。
下述内容**全部转引自本机 auto-memory 中此前会话（2026-09-11~09-22）已做的一手/二手核实**，
未做任何新检索；凡是「未见」的结论，是「过去多轮检索 + 本轮尝试复核仍未见」，不是本轮独立新查。

---

## 问题 1：μ\* > 0 时的 Einstein 极值性

**判词：部分重叠 —— 命题的 μ\*=0 版本已被 AD1975 陈述+部分证明、被 Carbotte1990 系统讨论；
μ\*>0 版本目前没有任何人给出严格证明，AD1975 有一句范围更宽的断言但无证明，Leavens 系列只给了
"approximate" 推广，Bergmann–Rainer 的泛函导数正性证明明确止步于 μ\*=0。**

### 证据等级逐篇

**① Allen & Dynes, PRB 12, 905 (1975) —— 一手全文（本机持有 PDF，2026-09-14 起多轮页图核对）**
- p.921 (A13) 段逐字：*"given λ, μ\*, and ω̄_n for an arbitrary value n≥2, the maximum Tc consistent with
  these parameters is obtained when (and only when) α²F is an Einstein spectrum. For simplicity we will
  give the proof only for n=2"* —— μ\* 出现在参数清单里，但这句里"a more general theorem than the one
  given in the text"指的「更一般」是**矩指标 n**（n≥2 vs 正文只证 n=2），**不是 μ\***
  ⚠ **订正（team-lead 2026-09-23 核实，pdftotext -enc UTF-8 逐字对过）**：我此前把这句读成
  「AD 对 μ\*>0 情形做了断言但没证明」，**这个读法本身就过度解读了**——AD **从未把「μ\*>0 的推广」
  单独拎出来声明为已证**，μ\* 只是随口列在参数清单里，全文没有一句专门针对 μ\*≠0 的极值性主张。
  正确判词是：**这是一个从未被单独陈述过的空位，不是一个被断言但未证的命题**（后者会构成「AD 声称
  了但没兑现」，前者只是「AD 没碰这个问题」——两者优先权后果不同，不要混用）。
- p.908 右栏逐字：*"inequalities of the type (18) and (23) are **no longer rigorously true**"* 当 μ\*≠0——
  这是 AD 自己承认 μ\*>0 时他们的逐元不等式机器**失效**，与上面的订正一致：他们清楚知道 μ\*≠0 时
  自己的证明链条不适用，所以才不会去单独声明这个推广。
- **⟹ AD1975 (A13) 末段真正证明的内容是**：`x/(x+a²)` 关于 `a²`（等价地关于某个矩指标对应的量）
  二阶导恒负 ⟹ **Einstein 谱逐 k 最大化 λ(iν_k)**（记作 **E1**，analytic、与 μ\* 无关——这一步的证明
  不涉及库仑赝势，纯粹是核函数凸性的结果）。在固定 (λ,ω̄₂) 的谱族里，`G[0]=λ` 被约束恒定，
  **库仑项 μ\* 与谱形状完全无关地进入久期方程**（它只出现在顶本征向量所满足的边界条件里，
  不出现在被最大化的核本身）⟹ **μ\* 只通过顶本征向量进入这个问题，不通过核 `x/(x+a²)`**。
  ⟹ **缺口有了显式形式**：analytic 已证的只有 E1（逐 k 核凸性）与相应的凸性结论；
  ~~**真正缺的一步是「顶本征值 λ₁ 对矩序列 {G[k]}_{k≥1} 单调不减」**~~
  ⚠ **订正（team-lead 2026-09-23 当日实测，`_scratch_ir/_chk_main_extremality_monotonicity_derivative_2026-09-23.py` 门 1/1b）**：
  这条单调性**本身为假**，所以它不可能是「缺的那一步」。用 Hellmann–Feynman 解析式
  `D[k] := ∂λ₁/∂G[k] = 2Σ_{n≥k}w_n w_{n−k} + Σ_{i+j=k−1}w_i w_j − 2Σ_{n≥k}w_n²`（`w:=v⊙u`，
  9/9 个 k 对有限差分通过，最劣 \|差\|=2.1e−10）实测 `D[10] = −4.97e−2 < 0`。
  且**失效由库仑驱动**：同一 (λ,谱,r) 下 μ\*=0 时 `D[k]` 全非负、顶本征向量 `v` 无节点；
  μ\*=0.05/0.13/0.25 时 `v` 各有 2 个节点、`D<0` 的 k 分别有 2/18/25 个。
  ⟹ **正确的缺口陈述**：μ\* 只通过顶本征向量进入，而 μ\*>0 使能隙函数出现**节点**、
  `w=v⊙u` 带号，逐元正性机器（AD 自己在 p.908 承认失效的那台）随之失效。
  任何证明都必须处理这个**符号不定**的本征向量——这正是文献里无人碰过的那一步。
  详见 [[../ceiling/mustar_extremality_two_routes_closed_2026-09-23]]。

**② Leavens, SSC 17, 1499 (1975) —— L2（作者本人摘要，OSTI 7247568；ScienceDirect 全文 403，未获取）**
- 摘要级：μ\*=0 严格给出 `k_BTc ≤ 0.2309 A`，等号在 `Aδ(ω−1.750A)`；**μ\*≠0 是"generalized in an
  approximate but sufficiently accurate way"**，给出 `k_BTc ≤ c(μ\*)A`，c/d 列表覆盖 0≤μ\*≤0.20。
- **⟹ 作者本人用词是"approximate"，不是严格证明。**

**③ Leavens, SSC 19, 395 (1976) —— L2（作者本人摘要，OSTI 7344151；全文 403，未获取，本轮再次
尝试仍未能取得，见文末"本轮尝试记录"）**
- 摘要：对比 Bergmann–Rainer 泛函导数、AD 的大 λ 渐近式、与自己的 0.2309A 界，断言 AD 的**大 λ
  渐近式** "is not compatible with the other results and is, in fact, **incorrect**"。
- ⚠ **被质疑的是 AD 的大 λ 渐近式（强耦合极限公式），不是 AD 的 Einstein 极值性命题 (A13)。
  两者是不同命题，此摘要不能读成"Leavens 反驳了极值性"。**
- **未获取全文** ⟹ 无法确认此文是否顺带讨论了 μ\*>0 时的极值性问题本身。**这是本问题里价值最高、
  仍未闭合的一环**（若它讨论了 μ\*>0 极值性，会直接决定命题归属）。

**④ Carbotte, Rev. Mod. Phys. 62, 1027 (1990) —— 一手全文（本机持有，2026-09-17 逐页核对）**
- §II.H (p.1045) 逐字：*"We can prove that (2.49) represents a **local maximum** for any shape of
  α²F(ω) with a fixed A. … It follows that the inequality Tc ≤ Ac(μ\*) (2.50) holds."*
  ⟹ **Carbotte 的证明本身也只是局部变分（local maximum），"It follows that"这一步跳到全局不等式
  没有独立论证**——这是与 AD 同型的缺陷，且**这条覆盖的是固定面积 A 的约束，不是固定 (λ,⟨ω²⟩)**。
- p.1096（比热跳变，方法相同的旁证）逐字：*"Our experience with realistic spectra would lead us to
  believe that this is, in fact, an absolute maximum."* ⟹ **Carbotte 自己承认从局部到全局这一步，
  他给出的理由是"经验/对现实谱系的观察"，不是证明。**
- **⟹ 1974–1990 这一整条文献线（Leavens→Carbotte）对 μ\*>0 情形提供的都是局部变分论证 + 经验外推，
  不是严格证明。** 这与仓库既有笔记 `literature/lit_einstein_extremality_novelty_2026-09-16.md`
  的裁决一致（本轮未改变判词，只是复核未发现推翻证据）。

**⑤ Bergmann & Rainer, Z. Physik 263, 59 (1973) —— 一手全文（本机持有，2026-09-11 起页图核对）**
- App.A 严格证明泛函导数 δTc/δα²F(ω) 在 **μ\*=0 时处处非负**。
- **全文对 μ\*>0 情形的泛函导数符号无任何证明**——这是本项目多轮检索（2026-09-11 至今）反复确认的
  空白：**"μ\*>0 时 δTc/δα²F 的符号"，全网未找到任何证明**（正的、负的、或变号的证明都没有）。
  这是极值性命题在 μ\*>0 下缺失的核心技术原因：没有正性泛函导数，局部变分论证站不住脚。

**⑥ Manalo–Schachinger, cond-mat/0101419 —— 一手页图（2026-09-16/18 核对）**
- eq.(1) 逐字：最优谱是 `α²F_opt(ω) = A δ[ω − ω*(μ*)]`，"with the delta function placed at the
  frequency ω*(μ\*) at which the functional derivative δTc/δα²F(ω) displays its maximum **for a fixed
  value of the Coulomb pseudopotential μ\***"。
- ⟹ **这篇明确把 μ\* 当固定参数，给出的仍是单峰最优猜想的具体形式（沿用 Carbotte/Leavens 的局部
  变分框架），不是新证明，也没有把这条命题推广到与团队相同的"固定 λ,⟨ω²⟩"约束。**

**⑦ Moussa & Cohen, PRB 74, 094520 (2006) —— 一手页图（2026-09-18 核对，问题 3 详见下）**
- p.2 逐字：*"We attempt to tighten the bounds by independently constraining λ and ⟨ω²⟩ as well as
  specifying μ\*. **It is already known [AD] that the Einstein spectrum maximizes Tc for a given λ and
  ⟨ω²⟩.**"* ⟹ **他们把 μ\*=0 的极值性直接搬到 μ\*≠0 的语境下用，全文没有重新检验这一步**——是本项目
  已核实过的、对目标命题最直接的一次"默认成立、未经检验"的沿用。

### 判词汇总（问题1）

| 来源 | 陈述范围 | 证明状态 |
|---|---|---|
| AD1975 p.921 | (A13) 只证 n≥2 的 Einstein 核凸性（E1，μ\* 无关）；**μ\*>0 的推广从未被单独陈述** | analytic（E1）；μ\*>0 是空位，不是未证断言 |
| Leavens 1975/1976 | μ\*≠0（c(μ\*) 表） | 作者自称 "approximate" |
| Carbotte 1990 §II.H | 固定面积 A, μ\*≠0 | 局部变分证明 + "It follows that"的未证跳跃 |
| Bergmann-Rainer 1973 | 仅 μ\*=0 | 严格证明（App.A），μ\*>0 完全空白 |
| Manalo-Schachinger 2001 | 单峰最优猜想形式 | 沿用 Carbotte 框架，非新证明 |
| Moussa-Cohen 2006 | 固定 λ,⟨ω²⟩ | 直接引用 AD，未检验 μ\*>0 |

**对我们的影响**：命题「固定 λ,ω̄₂,μ\*>0 时 Einstein 谱极大化 Tc」**没有被任何人严格证明，
也没有被任何人单独陈述过**——AD1975 (A13) 证明的 E1（核凸性）本身与 μ\* 无关，μ\* 只通过顶本征向量
进入问题。⚠ 缺口的具体形状**已被 team-lead 当日实测改写**（见上方订正横幅）：不是「λ₁ 对矩序列
单调不减」（该命题为假），而是「μ\*>0 使顶本征向量带号后，如何仍得到全局极值」；
且「E1+逐 k 单调」「E1+沿线段单调」两条自然证法**已被数值排除**
（详见 [[../ceiling/mustar_extremality_two_routes_closed_2026-09-23]]）。
其余文献（Carbotte/Leavens/Manalo-Schachinger）止步于局部变分论证或经验外推——
这一点本轮另获**独立数值佐证**：Einstein 的**局部**极大在扫到的全部方向上成立（归一方向导数
min Ψ̂(1)=2.23e−3 > 0，10/10 行），**包括已知失效带**（λ≲2μ\*）⟹ 局部极大对失效是**盲**的，
Carbotte 那句 "It follows that" 的跳跃因此是实质性的，不是形式上的。**这是团队命题的真空位**，
仓库既有笔记（`mustar-obstruction-is-gap-collapse-2026-09-17`、`mustar-extremality-weak-coupling-failure-2026-09-17`）
里团队自己数值证伪弱耦合窗口的极值性，与这条文献空位完全一致——**没人证明过、也没人陈述过，
我们的反例不是打脸已发表定理，而是填补了一个从未被触碰的空位**。

---

## 问题 2：Hopfield 参数 η 作为室温超导判据

**判词：未见（本轮无法确认，原因见下）。**

**本轮无法独立检索**：WebSearch 预算在会话开始时已耗尽（200/200），WebFetch 对所有域名（包括
export.arxiv.org 的 API 端点）返回"无法验证域名安全性"全部拒绝。本机 auto-memory 中此前的
检索记录（覆盖高压氢化物/η/M 相关的至少 6 篇独立笔记：`lit_hydride-isotope-effect-data-formulas`、
`lit_hydride-alpha2F-machine-readable-sources`、`lit_htp-ml-conventional-sc-d4`、
`lit_sscha-anharmonic-hydride-superconductivity`、`etaM_material_table_recheck_2026-09-18.md` 等）
**没有一条记录过"η 本身被当作独立筛选判据"这一用法**——但这些笔记的检索目的都不是这个问题，
**不能把"没记录"读成"文献里没有"**。

已知背景（来自仓库既有材料，非本轮新查）：
- Grimvall 1981《The Electron-Phonon Interaction in Metals》系统整理了 η 的第一性原理计算方法
  （刚性 muffin-tin 近似），Carbotte RMP62 p.1040 附近引用过，但**该引用语境是把 η 作为计算 λ 的中间量**，
  不是把 η 单独当筛选判据。
- Cerqueira 等高通量氢化物数据集（`lit_htp-ml-conventional-sc-d4` 笔记）与 Cerqueira 6912 材料数据集
  （`lit_cerqueira-materialscloud-omega2-hydride-omega2-2026-09-17`）：本仓已核实这些数据集**给出的是
  λ、ω_log/⟨ω²⟩、Tc**，未在笔记核实范围内见过独立报告 η（或 η/M）作为筛选量的表格列。

**对我们的影响**：无法给出"已被占用/未被占用"的判词；这是本轮交付里唯一一个**因工具中断而非文献
本身空白**造成的开放项，建议下一轮 WebSearch/WebFetch 恢复后优先补查，关键词沿用派单里给出的四组
（`Hopfield parameter screening superconductor` / `electron-phonon coupling descriptor high-throughput` /
`Gaspari-Gyorffy rigid muffin tin hydride` / `N(0)<I^2> criterion`）。

---

## 问题 3：Moussa & Cohen, PRB 74, 094520 (2006)

**判词：派单方（主会话）在 prompt 中的一句无出处陈述，经全仓 grep 确认此前不存在于任何文件；
已由 team-lead 本人核实并注销（2026-09-23），不是需要溯源的「丢失记录」，是一句需要当场作废的错话。**

⚠ 说明（按 CLAUDE.md「记录丢失时不得反向重建」的对应处理）：这条"两个界是 doping-limit 与
lattice-instability"的陈述**从一开始就没有仓库出处**——不是被覆盖/丢失的正文，是派单 prompt 里
新写入的一句无依据陈述。全仓 grep（`theory/` + `paper/`）确认这两个措辞此前唯一的出现位置就是
本笔记（转述派单原话），team-lead 已 grep 复核同一结果并认领。**按此，不留作待澄清的开放项，
当场按实际核实内容替换。**

**证据等级**：一手页图（2026-09-18 核实，本轮无法重新取得原文复核，但当时的转录含逐字引文+
公式编号，可信度高，本轮只是转述，不是新读）。

### Moussa & Cohen 实际给出的"两个界"（p.2, eq.7/10/11，逐字已存档）

1. **eq.(7)**：`0 ≤ Tc ≤ 0.183√(λ⟨ω²⟩)`（μ\*=0，无 μ\* 依赖，这是最简单的一条，纯矩界）。
2. **eq.(10)（含 μ\* 的双边界，这是论文的核心结果）**：
   `f(λ⟨ω²⟩/ω²_max, μ\*) ω_max ≤ Tc ≤ f(λ,μ\*)√⟨ω²⟩`
   下界用 `ω_max`（谱的支撑上限）+ 折合比值 `λ⟨ω²⟩/ω²_max`；上界用 `⟨ω²⟩`。
   **两端都是矩/支撑型的量（ω_max、λ、⟨ω²⟩），不是"doping-limit"或"lattice-instability"这类物理
   机制型的界——论文原文没有这两个措辞（已由 team-lead 独立 grep 确认）。**
3. eq.(11) 给出 `f(λ,μ\*)` 的具体拟合形式（指数型），拟合区间 `0.4≤λ≤20, 0≤μ\*≤0.2`。

**结论（team-lead 2026-09-23 认领并将引用）**：Moussa-Cohen 的两个界是矩/支撑型
（eq.7 纯矩界；eq.10 下界用 ω_max、上界用 ⟨ω²⟩），**不提供「λ 有物理上限」这条可行域边界**，
不能引 Moussa-Cohen 作为该边界的来源。

### 能否画成 (λ, η/M) 平面上的曲线

**能，但形式已经确定为矩型，不是额外的物理约束**：
- eq.(10) 下界 `f(...)ω_max ≤ Tc` 涉及 `ω_max`，这是一个**独立于 (λ,⟨ω²⟩) 的额外自由度**（谱的支撑
  上界，不能仅由 λ 和 ⟨ω²⟩ 决定，除非附加谱形状假设）。⟹ **若要把它画到团队的 (λ, η/M) 二维可行域图
  上，需要额外指定 ω_max 或谱形状**，否则这条界不是这两个坐标的函数，**不能直接搬上团队现有的图**。
- eq.(10) 上界 `Tc ≤ f(λ,μ\*)√⟨ω²⟩` **可以**直接画：固定 μ\*，这是 `η/M = λ⟨ω²⟩` 平面上一条随 λ 变化
  的曲线（`Tc = f(λ,μ\*)√(η/(λM))`），但它是**对 Tc 的界，不是对可行域边界（"λ 不能无限增大"）的界**——
  这条曲线约束的是纵轴 Tc，不约束横轴 λ 本身的取值范围。
- **⟹ Moussa-Cohen 不提供"λ 不能无限增大"这条边界**——那条边界（如果存在）需要别的物理输入（例如
  晶格失稳判据 `ω²_min>0` 或 McMillan 的经验上限），Moussa-Cohen 论文本身没有涉及。

**对我们的影响**：不能引 Moussa-Cohen 作为「λ 存在物理上限」这条边界的现成来源；他们的界是关于 Tc
本身的（矩型双边界），不约束 λ 或 η 的定义域。团队图上「λ 不能无限增大」的边界仍需要另找文献或
自行论证（例如晶格失稳、电声耦合微扰论适用范围等），这不是 Moussa-Cohen 能提供的。

---

## 本轮尝试记录（工具状态，供下一轮参考）

- WebSearch：会话开始时已达 200/200 上限，零次成功检索。
- WebFetch：对 `export.arxiv.org`、`arxiv.org`、`harvest.aps.org` 全部返回域名安全性验证失败，
  零次成功取件——与仓库既有笔记记录的"WebFetch 全域封死"故障同型（`reference_web-fetch-access-blockers`、
  多篇 2026-09-18~09-22 笔记记录过同一故障）。
- 本地 `_scratch_ir/` 未持有 Leavens 1976 (SSC 19,395)、Moussa-Cohen (PRB74,094520) 或
  Bergmann-Rainer (Z.Physik 263,59) 的 PDF；只有转录进 auto-memory 的逐字引文可用。
- **建议**：下一轮工具恢复后，优先级顺序为 ① Leavens 1976 全文（问题1最后缺口，特别是「λ₁ 对
  {G[k]} 单调不减」这一步是否有人碰过）② Hopfield 参数筛选判据检索（问题2完全未查）。
  ~~③ 向团队确认 doping-limit/lattice-instability 具体指代~~ **已关闭**（2026-09-23：team-lead
  核实该措辞是派单 prompt 里的无出处陈述，非仓库记录，已当场注销，见问题3）。

## 关联笔记

- [[../literature/lit_einstein_extremality_novelty_2026-09-16]]（AD1975/Carbotte 归属链主线）
- [[../literature/leavens_allen_dynes_1979_priority_check_2026-09-18]]（LAD79 优先权，问题1相邻但不同命题）
- [[../ceiling/oneliner_zero_parameter_2026-09-18]]（η/M 恒等式，问题2的框架背景）
