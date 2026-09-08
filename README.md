# 超导理论 · autoresearch 工作区

建立日期：2026-09-07
状态：**理论推进中**（P1 同位素四因子分解 v0.7、P2 顶角判据 v0.4；方向已定，见 `decisions/`）

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
- [x] 2026-09-08 启动氢化物实验数据表 + 同位素公式核对（researcher agent）→ `theory/data/hydrides_experimental.csv`、`surveys/P1_isotope_data_formulas_2026-09-08.md`、`library/index/P1.json`（17 PDF）
- [x] 2026-09-08 P1 推导 v0 → reviewer 发现 3 个致命问题（McMillan 饱和伪象、主公式结构、部分同位素缺推导）→ **v0.2** 修正并新增"谱重分布通道"（Errea 2015 表 I）与 H3S 第一轮账
- [x] 2026-09-08 求解器 v0 到位（32 测试全绿）；P1 v0.3（二轮审稿 N1–N4）→ v0.4（数值：闭式高估库仑项 2–4×）
- [x] 2026-09-08 **第一篇闭环记录** `theory/notes/2026-09-08_P1_round1.md`
- [x] 2026-09-08 二轮/三轮原文抽取：H3S α_exp=0.3（非 0.44）、LaH10 α_exp=0.57/0.62（文献 0.46 无一手来源）、YH9 无数据、LaD10 结构未确认 → `notes/2026-09-08_P1_round2_data.md`
- [x] 2026-09-08 **Mishra Figshare 数据集**（H3S/D3S × 谐波/非谐 × fsr/fbw × adia/na）→ H3S 四因子数值分解闭合；DOS×顶角交叉项不可忽略 → `notes/2026-09-08_P1_round3_mishra_dataset.md`；推导稿 v0.6
- [x] 2026-09-08 P2 公式抄录到位（PGS 原始论文不在 arXiv；EPL 42,667 = cond-mat/9710159 等替代源；Schrödi 2020 挑战 Q_c 图像）→ `surveys/P2_vertex_formulas_2026-09-08.md`
- [x] 2026-09-08 **P2 启动稿** `theory/derivations/P2_vertex_sign_criterion.md`：从 Mishra 数据反推顶角修正的同位素标度 δ_D/δ_H ≈ 0.66–0.72 ≈ M^{-1/2}；顶角效应 160→200 GPa 增强 3×；LaH10 顶角路径基本关闭
- [x] 2026-09-08 reviewer 核验 P1 v0.6 §4.1 + P2 v0：无致命；s_λ 更正（0.95→1.27）、交叉项加误差棒（四点同号为主证据）、160 GPa 论证改绝对量、数据集身份（Figshare = 2510.25902）→ P2 v0.1
- [x] 2026-09-08 求解器第二阶段（DOS + 顶角钩子，40 测试）；常 DOS 路径与 EPW 标定一致到 1–2 K；DOS 路径 α 移动对但 Tc 差 6–8 K
- [x] 2026-09-08 求解器 `fix_mu="conserve"` 三次迭代：参考态修正后非谐 Tc 与 EPW fbw 差 <1 K，谐波差 5–7 K；全带表因求积混叠需 n_eps≥1600（400 s/次）。**定稿：DOS 路径用于 α 差分，不再深挖绝对 Tc**
- [x] 2026-09-08 **P2 数值**：GPS99 一阶顶角在裸 E_F=24 eV 下符号错、量级差 10–20×（四个解析极限复现 <1%）→ P2 v0.2 重定位为"有效 E_F 由什么决定"（与 P4 合并）
- [x] 2026-09-08 **P2 v0.3：DOS 缺损泛函** P_s/λ = (ω₀/2)∫N̄/(|ε|+ω₀)² − 1 —— 真实 H3S DOS 给 −0.103 vs 第一性原理 −0.095，压力比 0.37 vs 0.37，E_F^eff = 0.96 eV；P2 判据 = 费米面附近 ~ω₀ 内的 DOS 缺损（与 P4 合并）
- [x] 2026-09-08 researcher：2608.14949 无氘化体系顶角数据；第二家族 H/D 顶角需向作者索取（用户决定）
- [x] 2026-09-08 reviewer 核验 P2 §5：推导正确；压力比吻合依赖 EPW 窗口、高能空带贡献 15–30%、判据降级为静态区部分答案 → P2 v0.4
- [x] 2026-09-08 LaH10 DOS 形状（文献）：与 H3S 同在 vHs 下坡但更平坦 → 顶角压制约 H3S 一半、同号；LaH10 α>½ 顶角路径关闭
- [x] 2026-09-08 **P2 第一轮闭环记录** `theory/notes/2026-09-08_P2_round1.md`
- [x] 2026-09-08 其他材料 adia/vertex Tc 数据：Mishra 2507 Table I 确认我们的 FSR 提取；MgB2 σ 带（Cappelluti 2002）为动态区增强的反例；H3S/D3S 圆顶错开 ~30 GPa（结构侧第五因子）
- [ ] P2 E_F 扫描（implementer）；(Q, ω_q) 动态版泛函（implementer）；求解器 muchem 方向核查（implementer）
- [ ] P2 下一步：式 (5.1) 的 (Q, ω_q) 动态修正版；求解器顶角钩子接入
- [ ] 用户决定：是否向 Mishra/Margine 发数据请求（`theory/notes/draft_email_mishra_margine.md`）
- [x] 用 Mishra 的 H3S.dos 标定求解器 DOS 实现（非谐 <1 K，谐波 5–7 K）
- [ ] 用求解器扫 DOS×顶角交叉项
- [x] 2026-09-08 D1–D4 方向卡填入现状摘要与切入点（commit 15fd936）；D5–D7 保留为参考
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
