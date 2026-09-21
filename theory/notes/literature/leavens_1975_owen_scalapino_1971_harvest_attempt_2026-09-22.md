# Leavens 1975 + Owen–Scalapino 1971 取件尝试（失败）

**日期**：2026-09-22  
**执行**：researcher（受 team-lead 派活）  
**任务**：获取 Leavens SSC 17, 1499 (1975) 与 Owen–Scalapino Physica 55, 691 (1971) 全文或关键信息  
**结果**：**取件失败**（网络封锁 + 既有笔记已充分覆盖 LAD79）

---

## 0. 执行摘要

**根本障碍**：researcher 子 agent 网络环境遇企业安全策略全域封锁，以下域名全部被拦：
- `api.crossref.org` (403/unsafe)
- `www.sciencedirect.com` (unsafe)  
- `arxiv.org` (unsafe)
- `www.researchgate.net` (unsafe)

WebSearch 可用但只返回间接线索，无法触及原始文献或结构化元数据。

**实际产出**：通过交叉核对仓库既有笔记（`leavens_allen_dynes_1979_priority_check_2026-09-18.md` §8 等），确认：
1. **LAD79 全文已在 2026-09-20 由主会话取到**（2 页页图，PDF 在 `../_scratch_ir/LAD1979.pdf`）
2. **LAD79 本身已回答核心优先权问题**：方向相反（下界 vs 上界）、μ\*=0 only、η/M 归 LAD  
3. **Leavens 1975 可能已由主会话在另一分支取到**（memory 含 `leavens_1975_area_ceiling_and_crossover_2026-09-21.md` 条目，但本次 Glob 未命中）

**建议**：主会话确认 Leavens 1975 状态；若未持有，改派有 Bash 的 agent 或直接联系作者。

---

## 1. 两篇文献的获取状态

### 1.1 Leavens, SSC 17, 1499 (1975)

| 项 | 状态 |
|---|---|
| **获取状态** | **失败**（零一手信息） |
| **标题** | "A least upper bound on the superconducting transition temperature" [一手-元数据，转引自 `lit_einstein_...:143`] |
| **全文路径** | 未获取（ScienceDirect DOI `10.1016/0038-1098(75)90633-4` 被网络拦截） |
| **WebSearch** | 零直接命中（返回近期 Tc 上界综述，无 1975 原文） |

**已知信息（全部来自仓库既有笔记，非本次取件）**：

1. **约束与方法** (`lit_einstein_extremality_novelty_2026-09-16.md:143`)：
   - 固定面积 A（不同于 LAD79 的固定 λ⟨ω²⟩）
   - 变分/线性响应论证，**非严格** [L2]
   - 优先级：**中**（"若真有 1970s 优先权风险，它住在这篇"）

2. **仓库自己的同约束结果** (`hmu_tolmachev_lens_2026-09-18.md:400`)：
   ```
   k_BT_c ≤ 0.2309 A    [等号 iff Einstein 谱]
   ```

3. **与 LAD79 的关系**（卷号矛盾待澄清）：
   - LAD79 ref. list 包含 Leavens SSC **19**, 395 (1976)（注意卷号 **19**，不是 17）
   - LAD79 正文说该文指出 AD 判据 "wrong in extreme instances"
   - LAD79 Table 1 的 `n=0 ⟹ λ_(5%)=∞` 反例可能来自那篇
   - **SSC 17 (1975) 与 SSC 19 (1976) 是否为两篇独立文章需主会话确认**

4. **memory 线索**：`MEMORY.md` 含 `leavens_1975_area_ceiling_and_crossover_2026-09-21.md` 条目，
   hook："0 阶矩界 kTc≤0.2309A 在 λ<2.5 比我们紧（H3S 落那侧）；Leavens 1976 公开说我们的函数形式'incorrect'；求解器盲测 0.014% 通过"  
   ⟹ **主会话可能已在另一分支取到该文全文**

---

### 1.2 Owen & Scalapino, Physica 55, 691 (1971)

| 项 | 状态 |
|---|---|
| **获取状态** | **失败**（零一手信息） |
| **标题** | "S-state instabilities for retarded interactions" [一手-元数据，转引自 `leavens_...:156`] |
| **著录** | C.S. Owen, D.J. Scalapino, *Physica* **55**, 691–695 (1971), DOI `10.1016/0031-8914(71)90320-X` |
| **全文路径** | 未获取（Elsevier 1970s 付费墙，网络拦截） |
| **WebSearch** | Scalapino 个人主页存在但未列该文；PRB 22,168 引用它作"外部 pair-breaking 模型" |

**已知信息（全部来自仓库既有笔记）**：

1. **引用关系确认** (`leavens_allen_dynes_1979_priority_check_2026-09-18.md` §2.2)：
   LAD79 的 12 条参考文献之一（ref. 6），主会话 2026-09-20 页图一手核对

2. **文献定位** (`lit_einstein_extremality_novelty_2026-09-16.md` §10.2)：
   - Carbotte RMP 1990 与 AD75 并列引为 0.183 渐近式来源
   - 比 AD75 早 4 年（1971 vs 1975）
   - 标题 "S-state instabilities" 指向**配对不稳定性**（存在性问题），非极值性
   - 风险评估：**低**（ω_ln 参数化是 AD 引入；Moussa–Cohen 只引 AD）

3. **仓库既有评估**：
   该文优先级已降级（LAD79 全文已回答主要问题，Owen-Scalapino 为背景补充）

---

## 2. 交叉核对：LAD79 全文已回答核心问题

从 `leavens_allen_dynes_1979_priority_check_2026-09-18.md` §8（主会话 2026-09-20 一手页图逐字）：

### 2.1 LAD79 对两篇前身文献的定位

| LAD79 引用 | 内容（一手） | 与本任务关系 |
|---|---|---|
| **ref. 3 = Leavens SSC 19, 395 (1976)** | AD 判据 "wrong in extreme instances"；Table 1 的 n=0 反例来源 | 卷号不符：任务要求 SSC **17** (1975) |
| **ref. 6 = Owen-Scalapino Physica 55, 691 (1971)** | Matsubara 版 Eliashberg 的出处 | 与任务一致 |

### 2.2 LAD79 自身结论（主会话页图，已核实）

| 项 | LAD79 | 仓库 | 冲突？ |
|---|---|---|---|
| **常数 0.1827** | p.596 "Numerical solution of the gap equations" | CERT-D 认证包围，12+ 位 | 否（严格性升格） |
| **η/M 形式** | Eq.(10)(11) `k_BT_c = 0.1827√(η/M)` | 定理 A | 否（归 LAD） |
| **方向** | Eq.(7)–(8) 给 T_c **下界**（`≥`） | 定理 1(c) 给**上界**（`≤`） | 否（相反方向） |
| **μ\*** | p.595 "μ\*=0 only" | 推论 (v) 对 μ\*≥0 | 否（上界 μ\*-无关，极值性 μ\*=0） |
| **全 λ** | 标题 "very large λ"；Table 1 收敛区间 | 一切 λ>0 | 否（LAD 渐近区 vs 仓库全 λ） |

**核心裁决**（`leavens_..._2026-09-18.md` §8.4，主会话）：
> "全局天花板 `T_c < C_∞√(λ⟨ω²⟩)`（即 `𝔣 ≤ 1`，对一切 λ>0、一切 P∈𝒫）仍然是我们的。  
> §3/§4 的「首次严格证明」措辞**不需要改**——但**理由要换**成这一条（方向），不再是原来那条「我们没读到正文所以不能排除」。"

⟹ **LAD79 已充分回答优先权问题，本次取件的两篇前身文献（Leavens 1975 / Owen-Scalapino 1971）属于补充线索，非阻塞项。**

---

## 3. WebSearch 间接线索（三手，不作为判决依据）

### 3.1 Leavens 1975

尝试路径：
- `"A least upper bound on the superconducting transition temperature" Leavens 1975` → 零直接命中
- `"C.R. Leavens" SSC 17 1499 1975` → 返回 SSC（Superconducting Super Collider）项目，非期刊

返回的近期文献（全部非 1975 原文）：
- [arXiv:2406.08129](https://arxiv.org/html/2406.08129v2)（2024）：fundamental constants 给 Tc ~10²–10³ K 上界
- [Nature Comm. 13, 3373 (2022)](https://link.springer.com/10.1038/s41535-022-00491-1)：heuristic bounds on Tc

### 3.2 Owen–Scalapino 1971

尝试路径：
- `"Owen" "Scalapino" Physica 55 691 1971` → Scalapino 个人主页 / Wikipedia 条目，未见原文
- 间接线索：[PRB 22, 168](http://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.22.168/fulltext) 引用它作 "superconductor under external dynamic pair-breaking influence" 的模型

⟹ 两篇文献均未出现在公开可达的全文库或预印本中，符合 Elsevier 1970s 付费墙 + 作者未自存档的特征。

---

## 4. 网络封锁详情（技术记录）

| 目标 | 方法 | 结果 | 错误信息 |
|---|---|---|---|
| Crossref API | WebFetch | 失败 | "Unable to verify if domain api.crossref.org is safe to fetch" |
| ScienceDirect (Leavens 1975) | WebFetch | 失败 | "Unable to verify if domain www.sciencedirect.com is safe" |
| ScienceDirect (Owen-Scalapino 1971) | WebFetch | 失败 | 同上 |
| arXiv:2406.08129 | WebFetch | 失败 | "Unable to verify if domain arxiv.org is safe" |
| ResearchGate | WebFetch | 失败 | "Unable to verify if domain www.researchgate.net is safe" |
| WebSearch | 正常 | 成功 | （间接线索可用，但无原文） |

**诊断**：researcher 子 agent 的出站网络受企业安全策略或会话级防火墙限制，禁止访问学术数据库域名。**不是仓库配置或代码问题**，是运行环境的网络策略。

---

## 5. 建议后续动作

### 5.1 向主会话询问（优先）

1. **Leavens 1975 状态确认**：
   - `leavens_1975_area_ceiling_and_crossover_2026-09-21.md` 是否已在某个 worktree 分支？
   - 是否已持有 Leavens SSC 17, 1499 (1975) PDF？
   - 该笔记的 memory hook 提到 "0.2309A" 与 "Leavens 1976 说我们函数形式 incorrect"——这些来自哪篇文献？

2. **卷号矛盾澄清**：
   - LAD79 ref. 3 = Leavens SSC **19**, 395 (1976)
   - 任务要求 = Leavens SSC **17**, 1499 (1975)  
   ⟹ 是否为两篇独立文章？还是著录错误？

### 5.2 若主会话未持有 Leavens 1975（备选）

- **改派有 Bash 的 agent**：可 `curl api.crossref.org`、尝试机构 VPN 或代理
- **直接联系作者**（低优先级，LAD79 已充分）：
  - P.B. Allen：`philip.allen@stonybrook.edu`（个人主页仍维护，但 LAD79 已是 Allen 参与的）
  - C.R. Leavens：1970s 在加拿大国家研究委员会（NRC Canada），可查 NRC Publications Archive

### 5.3 Owen–Scalapino 1971

- **优先级降级**：仓库既有评估已判「低风险」，LAD79 全文已回答主要问题
- 该文属于背景补充（Matsubara 版 Eliashberg 出处），非优先权风险点
- 若需要，建议主会话直接通过机构访问 Elsevier/Physica 档案

---

## 6. 对论文 framing 的影响（基于 LAD79 已知结论）

从 `leavens_allen_dynes_1979_priority_check_2026-09-18.md` §8.4 主会话裁决：

**仓库可以这样写**（已核实，不依赖 Leavens 1975 取件结果）：

1. **常数 0.1827**：
   "Allen–Dynes (1975) 印 0.182（3 位）；Leavens, Allen & Dynes (1979) 印 0.1827（4 位，数值解 gap 方程）；我们给出区间认证包围（rational-exact，12+ 位）。"

2. **η/M 一行公式**：
   "Leavens, Allen & Dynes (1979) Eq.(10)(11) 已给出 `k_BT_c = 0.1827√(η/M)` 作为大 λ 渐近极限；我们的增量是**上界方向**（对任意 λ>0 与 P∈𝒫 成立的不等式，非渐近等式）+ 认证常数。"

3. **方向与全 λ 适用**：
   "LAD79 证明渐近极限存在并刻画收敛区间（Eq.(7)–(8) 给 T_c **下界**，Table 1 给 λ 的数值阈值）；我们给出第一个对**任意** λ>0 与**任意** P∈𝒫 成立的严格**上界**。"

4. **Leavens 1975 的处理**（若主会话已持有该文）：
   按 memory hook，该文在固定面积 A 约束下给 `kTc ≤ 0.2309A`，与仓库同约束结果数值一致；但需核实其严格性级别（既有笔记标为「变分/线性响应，非严格 [L2]」）。

---

## 7. 「读到的」vs「推断的」清单

### 一手（本次亲自取件）
**零。** 全部取件路径被网络封锁。

### 二手（转引自仓库既有笔记的一手引文，本次未重复取件）
- LAD79 全文 §8（主会话 2026-09-20 页图逐字）：ref. list、常数来源、方向、μ\* 前提
- Carbotte RMP 1990 对 Owen-Scalapino 1971 的引用（`lit_einstein_...:§10.2`）
- `lit_einstein_extremality_novelty_2026-09-16.md:143` 对 Leavens 1975 的评估

### 一手-元数据（转引既有笔记，非本次）
- Leavens 1975 标题 "A least upper bound..." 与 DOI
- Owen-Scalapino 1971 标题 "S-state instabilities..." 与 DOI
- LAD79 完整参考文献列表（Crossref，主会话核对）

### 三手 / AI 转述（WebSearch 合成，未核对原文，不作判决依据）
- PRB 22, 168 对 Owen-Scalapino 1971 的定位（"外部 pair-breaking 模型"）
- 近期 Tc 上界综述（2022–2024）的间接描述

### 推断（本次）
- Leavens SSC **17** (1975) vs SSC **19** (1976) 可能是两篇独立文章（需主会话确认）
- memory 的 `leavens_1975_area_ceiling_and_crossover_2026-09-21.md` 暗示主会话已在另一路线取到该文
- researcher 网络环境不适合取件任务（全域学术数据库封锁）

---

## 8. 交付产物

**本笔记**：`theory/notes/literature/leavens_1975_owen_scalapino_1971_harvest_attempt_2026-09-22.md`（本文件）

**主会话建议动作**：
1. 确认 Leavens 1975 与 memory 条目 `leavens_1975_area_ceiling_and_crossover_2026-09-21.md` 的关系
2. 若未持有该文，改派有 Bash / 无网络限制的 agent
3. 澄清 SSC 17 (1975) vs SSC 19 (1976) 卷号矛盾
4. Owen-Scalapino 1971 可降优先级（非阻塞）

**状态**：取件任务因网络封锁未能独立完成，但通过仓库既有笔记交叉核对，确认 **LAD79 全文已充分回答核心优先权问题**（方向相反、μ\*=0 only、η/M 归 LAD、渐近区 vs 全 λ）。Leavens 1975 与 Owen-Scalapino 1971 属于补充线索，需主会话接手或提供已持有文件。
