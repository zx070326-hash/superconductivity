# 超导理论 · autoresearch 工作区

建立日期：2026-09-07
状态：**方向调研中**（7 个方向并行调研，等待选定聚焦方向）

---

## 0. 这个工作区做什么

从理论出发推进超导 Tc 的第一性原理预测与配对机制研究。第一阶段不写代码、不算材料，只做一件事：把领域的前沿方向摸清楚，选出 1–2 个 1–2 人、6–12 个月能做出实质结果的聚焦点。

---

## 1. 理论层级（作为所有方向的坐标系）

电声超导的 Tc 计算按近似层级排，每上一层去掉一个假设、多付一次算力：

| 层 | 方法 | 去掉的假设 | 典型误差来源 | 高通量可行性 |
|---|---|---|---|---|
| L0 | BCS：k_B Tc = 1.13 ħω_D exp[−1/N(0)V] | — | 单参数，无频率信息 | 只能定性 |
| L1 | McMillan / Allen-Dynes（λ, ω_log, μ*） | 单参数 V → α²F 的两个矩 | λ>2 系统偏离；μ* 经验值 | ✅ 已有 10³–10⁴ 规模数据库 |
| L2 | 各向同性 Migdal-Eliashberg 数值解 | 拟合公式 | μ* 与截止频率绑定；各向异性丢失 | 勉强（EPW 等） |
| L3 | 各向异性 ME / SCDFT / ab initio 库仑 | 费米面平均；经验 μ* | DFT 泛函；谐波声子 | ❌ 单材料级 |
| L4 | + SSCHA 非谐 / + 顶角修正 / + 能量依赖 DOS | 谐波近似；Migdal 定理；N(ε)≈N(0) | 计算量 10–100× | ❌ |
| — | 非常规（Hubbard / 自旋涨落） | 声子介导本身 | 模型级别就无定论 | ❌ |

高通量能到 L1，候选进入 L2–L4 才逐级获得定量意义。**筛选阶段的 Tc 只能当排序变量用，不能当绝对预测。**

---

## 2. 七个候选方向

| 编号 | 方向 | 核心问题 | 属于哪一层 |
|---|---|---|---|
| D1 | 超越 Migdal-Eliashberg | 顶角修正、能量依赖 DOS、各向异性、Allen-Dynes 强耦合失效 | L2→L4 |
| D2 | 库仑 / μ* 从头算 | SCDFT、ab initio 屏蔽库仑 Eliashberg、等离激元 | L3 |
| D3 | 非谐与量子核效应 | SSCHA、氢化物稳定性、MLIP 加速 | L4 |
| D4 | 高通量 + 机器学习 | 数据库、ML 代理、描述符漏斗、逆向设计 | L1 及其可靠性 |
| D5 | Hubbard / 非常规配对数值 | 2D Hubbard 基态之争、DCA/DΓA/DMRG/AFQMC、ab initio 下折 | 非常规 |
| D6 | 镍酸盐理论 | La3Ni2O7 双层两轨道模型、s± vs d、DFT+DMFT | 非常规（热点） |
| D7 | 平带 / 量子几何 / 莫尔 | 量子度规与超流刚度、魔角石墨烯配对、kagome | 新范式 |

D1–D4 是「第一性原理 Tc」主线的四个断点；D5–D7 是主线之外的前沿。

---

## 3. 目录约定

```
superconductivity/
├── README.md                  ← 本文件：总览、层级、方向表、进度
├── directions/                ← 每个方向一个文件夹
│   └── Dn_xxx/
│       ├── README.md          ← 方向卡：范围、核心问题、现状摘要、切入点、关键文献
│       └── notes/             ← 读文献/推导的工作笔记，一篇一文件
├── surveys/                   ← 调研 agent 的原始报告（只读，不改）
│   └── Dn_xxx_YYYY-MM-DD.md
├── library/                   ← 统一文献库
│   ├── index/Dn.json          ← 各方向文献元数据（agent 产出的 JSON）
│   ├── index.md               ← 自动生成的总索引
│   └── pdf/Dn/                ← PDF，命名 {year}_{first_author}_{arxiv_id}.pdf
├── decisions/                 ← 方向选择等决策记录，一事一文件，带日期
└── tools/                     ← fetch_arxiv.py 等辅助脚本
```

规则：
- `surveys/` 是 agent 原始产出，不手改；修正写到 `directions/Dn/README.md` 或 `notes/`。
- 文献只存一份，在 `library/pdf/Dn/`；方向卡用相对路径引用。
- 决策记录写清「选了什么、放弃了什么、为什么、什么条件下回头」。

---

## 4. 进度

- [x] 2026-09-07 建立目录结构，7 个方向并行调研启动（sonnet researcher agents）
- [x] 2026-09-07/08 7 份调研报告归档到 `surveys/`；文献 JSON 归档到 `library/index/`（D1–D7 全部到齐）
- [x] 2026-09-08 下载 core/important 级 PDF：86 篇，378 MB（16 条无 arXiv 版）
- [x] 2026-09-08 横向对比与聚焦建议：`SYNTHESIS_2026-09-08.md`
- [x] 2026-09-08 **方向已定**：高压超导的强耦合修正理论（理论为主 + 实验落实 + 适度计算验证），见 `decisions/2026-09-08_focus_theory_high_pressure.md` 与 `theory/PROGRAM.md`
- [x] 2026-09-08 启动氢化物实验数据表 + 同位素公式核对（researcher agent）
- [ ] P1（同位素系数四因子分解）推导初稿；求解器 v0（需先 `git init` 供 implementer 用 worktree）
- [ ] 各方向 README 填入现状摘要（D1–D4 优先，D5–D7 降级为参考）
- [ ] 补核 16 条无 arXiv ID 的元数据（多为 Nature/PRL 正式版，无预印本）

## 5. 工作区补充

```
theory/
├── PROGRAM.md        ← 研究纲领：核心论题、P1–P6 候选问题、闭环结构、里程碑
├── derivations/      ← sympy 脚本与推导稿
├── solvers/          ← 各向同性 Eliashberg 求解器（含库仑核、能量依赖 DOS、质量缩放、模型顶角核）
├── data/             ← 实验数据表 hydrides_experimental.csv、理论值表 hydrides_theory.csv
└── notes/            ← 每轮闭环一篇，固定结构见 PROGRAM.md §3
```
