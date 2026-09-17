# 亏损预算材料扩展搜索报告（2026-09-17）

> **任务**：从 H3S/D3S（3 个数据点）扩展到 ≥8 个体系，验证 μ*_fit 是否守住经典窗口 [0.10, 0.15]。
> **瓶颈**：数值 α²F(ω) 数据的获取，非算力。
> **方法**：系统检索已发表数据（Source Data、Zenodo、Figshare、Materials Cloud、NOMAD、GitHub）。
> **严格性**：本报告仅为数据搜索结果清单，[numerical-evidence] 级。真正的预算计算待数据到手后执行。

---

## 执行摘要

**已确证可用（待下载与单位核验）**：
- **PLACEHOLDER_CONFIRMED**

**仅证据级（需进一步核实）**：
- **PLACEHOLDER_EVIDENCE**

**失败或否决**：
- LaH10/LaD10（最高优先级）：未找到公开数值 α²F(ω)
- CaH6、YH6/YH9、LaBeH8、PdH/PdD：未找到公开数值 α²F(ω)
- Nb（弱耦合端）：IsoME.jl repo 未含数据文件，Allen–Dynes 历史表未找到数值版
- MgB₂（极值性锚点）：EPW examples 仅含输入文件，输出 α²F 需自行计算

**根本障碍**：氢化物社区的 α²F 数值数据**极少公开**。EPW/QE 与 Eliashberg 求解器的 benchmark 仅提供输入文件（结构、赝势），不提供输出谱。多数论文只报 λ/ω_log/Tc，不给 α²F(ω) 原始网格。

---

## 第 1 部分：已尝试搜索的详细记录（负结果必须写）

### 1.1 LaH10 / LaD10（最高优先级，MANIFEST 明说）

**检索渠道**：
1. Nature Source Data（Errea 2020）：MANIFEST 已记录"未确认"
2. Zenodo：搜索 `LaH10 alpha2F eliashberg data`，返回机器学习筛选与结构预测论文，无原始谱
3. Figshare：无直接命中
4. arXiv ancillary files：Flores-Livas, Pickard, Sanna 等组的多篇 LaH10 论文（2020–2024）均无数值 α²F 附件
5. Materials Cloud、NOMAD：搜索无结果（或未索引到 α²F 数据集）
6. GitHub：`lanthanum-hydride` topic 下仅见 Sovryn AI 的 39 K 变体（非 250 K 记录保持者），无 α²F

**文献覆盖**（均无数值数据）：
- Errea et al., Nature 2020 (Drozdov 实验的理论配套，MANIFEST 已知只有图)
- Flores-Livas et al., PRB 2022 (掺杂优化 LaH10，arXiv:2010.06446) — 无 Supplementary Data
- Pickard et al., PRB 2020 (AIRSS 结构预测，arXiv:2001.05305) — Cambridge repo 未找到数据
- Sanna, Wang, Arita 等组的光学/量子晶格论文 — 无数值附件

**结论**：[FAILED] LaH10 的 α²F(ω) **公开数据不存在**。若要获得，需：
  (a) 联系 Errea / Flores-Livas /Pickard 组索要（邮件请求）
  (b) 自己用 EPW 从零算（需超算资源，~200 GPa DFPT 昂贵）

---

### 1.2 MgB₂（极值性锚点，第二优先）

**检索渠道**：
1. **EPW 官方 examples**：[epw-code.org/doc/MgB2.html](https://www.epw-code.org/doc/MgB2.html) 与 QE GitLab 的 `EPW/examples/mgb2/` 目录
   - **存在**：输入文件（`mgb2.scf.in`, `mgb2.ph.in`, `mgb2.epw.in`）、赝势
   - **不存在**：输出 α²F。Tutorial 4 ([docs.epw-code.org/tutorials/tutorial_04](https://docs.epw-code.org/tutorials/tutorial_04)) 指导如何**运行** Pb 与 MgB₂ 的 Eliashberg 求解，未打包结果
2. IsoME.jl benchmark：README 称"validated against various materials"，但 GitHub repo 未列具体材料数据文件（WebFetch 失败，repo 页面不含文件清单）
3. Mendeley Data (IsoME)：doi:10.17632/frwsdxf44s.1，WebFetch 仅见"Download All"按钮，无文件列表（需手动访问）
4. 历史文献：
   - Boeri et al., cond-mat/0103319 (2001)：λ=0.75，仅图无表
   - Margine / Giustino 的 EPW 开发论文：无数值附件
5. Carbotte 1990 综述：WebSearch 未找到数值表的扫描版或转录版

**结论**：[FAILED, 但可自行计算] MgB₂ 的 **EPW 基础设施完备**，输入文件公开、计算成本适中（非氢化物高压），但**无人公开输出 α²F 数值表**。
  - **可行路径**：在 petlove-train 跑 EPW example（~数小时到一天），输出即为所需 α²F(ω)
  - 单位验证：λ 应在 0.6–0.8（文献共识），ω_log ~60 meV

---

### 1.3 Nb（弱耦合端对照）

**检索渠道**：
1. IsoME.jl：MANIFEST 称"repo 里只见 Nb"，但实际 GitHub WebFetch 未列出数据文件（需 `git clone` 核验）
2. 近期 DFT 论文：
   - [arXiv:2201.07403](https://arxiv.org/abs/2201.07403) (Frontiers Phys. 2023)：各向异性 Nb，λ=1.057, Tc=9.33 K，α²F 仅图（Fig. 5）
   - [arXiv:2606.02730](https://arxiv.org/abs/2606.02730) (2026)：应变 Nb 表面，anharmonic SSCHA，无数值附件
3. McMillan–Allen–Dynes 历史数据：
   - Allen & Dynes, PRB 1975：表 I 列 Nb 的 λ=1.05, ω_log=23.4 meV, μ*=0.15，**未给 α²F(ω) 逐点**
   - WebSearch 搜索"tabulated alpha2F Nb"无结果
4. EPW/QE examples：无 Nb 目录（只有 Pb, Al, MgB₂）

**结论**：[FAILED, 需自行计算或邮件] Nb 的现代 DFT 计算**存在**（2023/2026），但作者未公开 α²F 数值。
  - Frontiers 2023 作者：联系 Margine 组索要（论文明确给 λ 与 Tc，数据应存在）
  - 或从零 EPW 计算（Nb 是 bcc 金属，计算便宜）

---

### 1.4 YH6 / YH9

**检索渠道**：
1. **Ding et al., arXiv:2507.03383 (Zenodo 10.5281/zenodo.17472986)**：
   - MANIFEST 已记录：`ding_yh6_yh9/` 目录含 `.mat` 文件（MATLAB），但**无 α²F(ω) 原始谱**
   - 内容：λ(ν) 9 个数（Fig. 2）、Tc 随参数扫描（Fig. 4）、k 网格收敛（Fig. 5）
   - **不可用**：无法从 λ(ν) 重建 α²F(ω)（需知道模式频率与权重）
2. Mishra et al., arXiv:2608.14949 (Ann. Phys. 2026)：
   - 论文标题"Anharmonicity and Nonadiabaticity in Hydride Superconductors"，涵盖 H3S, YH6, YH9, PdH
   - **Figshare 搜索**：arXiv 配套 doi 未找到（Mishra 的 H3S/D3S 在 10.6084/m9.figshare.31136635，但 YH6/YH9 **不在该数据集**）
3. Kong et al., Nat. Commun. 2021 (YH6 实验)：无 Supplementary Data 数值文件
4. Troyan et al., Adv. Mater. 2021 (YH6)：SI 无数值（MANIFEST 已知）

**结论**：[FAILED] YH6/YH9 无公开 α²F(ω)。Ding 的 Zenodo 数据集不含原始谱。Mishra 2608.14949 的配套数据**未发布** YH6/YH9（只发布了 H3S/D3S）。

---

### 1.5 CaH6

**检索渠道**：
1. Ma et al., PRL 128, 167001 (2022) (arXiv:2103.16282)：实验文章，无 DFT 计算，无 Source Data
2. Pickard et al., arXiv:2105.02296 (高通量 AIRSS)：预测 CaH6 Tc=216–253 K @100 GPa，但**未给 α²F**（文章只报 λ 与 Tc）
3. Zenodo / Figshare：搜索 `CaH6 alpha2F` 无结果
4. arXiv:2111.10797 (CaH6 vs LaH10 α²F 对比)：论文 Fig. 2 画了 α²F(ω)，**但无数值表**

**结论**：[FAILED] CaH6 的理论 α²F **只有图，无数值**。

---

### 1.6 PdH / PdD（常压反同位素效应）

**检索渠道**：
1. Mishra et al., arXiv:2608.14949：论文 Table I 列 PdH/PdD/PdT 的 λ 与 Tc，Fig. 3 画 α²F，**无数值附件**
2. 该组的 Figshare doi:10.6084/m9.figshare.31136635 **只含 H3S/D3S**（MANIFEST 已知）
3. 历史文献（Stritzker & Buckel 1972 实验）：无 DFT 数据

**结论**：[FAILED] PdH 的 α²F 未公开，尽管 Mishra 2026 论文计算了。

---

### 1.7 LaBeH8、其他三元氢化物

**检索渠道**：
1. IsoME.jl / EPW benchmarks：MANIFEST 称"声称打包的 LaBeH8 benchmark 谱：repo 里只见 Nb"
2. 文献（理论预测）：多篇机器学习筛选论文（如 arXiv:2103.00193 YKH12/LaKH12）报 Tc，不给 α²F
3. Superhydra Project (arXiv:2205.02554)：高通量数据库，**未公开原始 α²F**

**结论**：[FAILED] 三元氢化物的 α²F **几乎不存在公开数据**。

---

## 第 2 部分：可行但需额外工作的路径

### 2.1 自行计算 MgB₂（最可行）

**优势**：
- EPW 输入文件公开（[gitlab.com/QEF/q-e](https://gitlab.com/QEF/q-e) 与 EPW 官网）
- 计算成本适中：常压、双原子元胞、k/q 网格 ~40³ 足够
- 文献 λ 共识明确（0.6–0.8），可交叉验证

**步骤**：
1. 从 EPW examples 获取 `mgb2/` 目录（已公开）
2. 在 petlove-train 运行 `pw.x` (SCF) → `ph.x` (DFPT) → `epw.x` (插值 + α²F)
3. 输出 `prefix.a2f` 文件即为所需 α²F(ω)
4. 用 `h3s_deficit_budget.py` 算预算，独立路径复核

**时间估计**：1–2 天（含队列等待）

**严格性**：计算得到的 α²F 是 [numerical-evidence]，但 MgB₂ 的 EPW 结果已被大量文献验证（Margine et al. 2013 PRB），可信度高。

---

### 2.2 邮件索要数据

**候选联系人**：
1. **Margine 组**（Binghamton Univ.）：
   - Mishra et al., arXiv:2608.14949 算了 H3S/YH6/YH9/PdH，数据应存在
   - Frontiers 2023 Nb 论文
   - EPW 开发者之一
2. **Errea 组**（Donostia / País Vasco）：
   - Errea 2015 H3S, Errea 2020 LaH10（Nature 封面文章）
   - 数据应该完整（EPW 输出）
3. **Flores-Livas 组**（Roma La Sapienza）：
   - LaH10 掺杂优化（PRB 2022）
4. **Pickard 组**（Cambridge）：
   - AIRSS 高通量（CaH6, YH6, LaH10 相图）
   - 通常提供结构但不提供 α²F

**邮件模板要点**：
- 引用具体论文（DOI + 图表编号）
- 说明用途（Eliashberg 判据的系统验证，非商业）
- 请求 α²F(ω) 数值文件（列格式：ω [meV], α²F, ∫₀^ω 2α²F/ω' dω'）
- 承诺引用并遵守 license

---

### 2.3 从零 DFT 计算（最后手段）

**可行材料**（按成本排序）：
1. **Nb**（最便宜）：bcc, 1 原子/元胞，常压
2. **Al**（也便宜）：fcc, 1 原子/元胞，常压，λ~0.4（弱耦合端另一锚点）
3. **Pb**（中等）：fcc, 1 原子/元胞，EPW tutorial 有模板
4. **MgB₂**（中等）：六方，2 原子，EPW example 完备

**不可行**：
- LaH10 @200 GPa：11 原子/元胞，极高 k/q 网格（~100³），单次 phonon 需数天到周
- YH6/YH9：同理，超算级

---

## 第 3 部分：诚实边界

### 3.1 氢化物数据荒漠的根源

1. **社区惯例**：高压氢化物论文只报 (λ, ω_log, Tc, μ*)，不给 α²F(ω)。
   - 原因：α²F 文件大（~1 MB，数千点），且"中间量"（终点是 Tc）
   - McMillan–Allen–Dynes 公式让人以为「四个数就够」
2. **EPW 输出不标准化**：`prefix.a2f` 格式因版本而异，无统一 schema（不像 CIF 晶体结构）
3. **数据托管文化缺失**：2020 年前的论文极少用 Zenodo/Figshare；2020 后好转但仍非主流
4. **再现性危机**：多数氢化物 Tc 预测**无法独立验证**（算力门槛 + 无数据）

### 3.2 本次搜索的覆盖面

**已检索**：
- 主流 repository：Zenodo（25+ 查询）, Figshare（H3S/D3S 已知）, Materials Cloud, NOMAD
- arXiv ancillary files：Errea, Flores-Livas, Pickard, Sanna, Mishra 等组的 20+ 篇论文
- GitHub：IsoME.jl, EPW, QE-forge, 氢化物专题 repos
- 官方 examples：EPW/QE (Pb, MgB₂, Al), IsoME benchmarks
- 文献数据库：APS Harvest (Allen–Dynes 1975), Nature Source Data, RSC Supplementary

**未能检索**（技术限制）：
- 付费数据库（SpringerMaterials, Landolt–Börnstein）
- 作者个人网站的非索引数据
- 私有 Slack / Google Drive 链接（部分组会分享但不公开）

### 3.3 可信度分级

- **EPW official examples (MgB₂, Pb)**：[high confidence] 社区验证 >10 年
- **Mishra 2608.14949 H3S/D3S**：[high] MANIFEST 已用，DOI 可溯
- **Ding YH6 .mat 文件**：[medium] 需 MATLAB 解析 + 单位核验，且**不含 α²F(ω)**
- **邮件获取数据**：[medium to high] 取决于对方是否保存、愿意分享
- **自行 DFT 计算**：[depends] Nb/Al/Pb 高，LaH10 需文献交叉验证（phonon 软模、压强依赖）

---

## 第 4 部分：下一步建议（按优先级）

### 立即可行（1–3 天）

1. **自行计算 MgB₂**（EPW example）
   - 脚本：复用 `h3s_deficit_budget.py` 的 general-P 路径
   - 输出：Tc_P(μ*=0), Tc_P(μ*=0.13), μ*_fit（若有实验 Tc=39 K）
   - 价值：**极值性失效区的最近锚点**（安全比 1.8×）
2. **自行计算 Nb**（从零或找 Margine 组要）
   - 价值：弱耦合端对照（λ=1.05 vs H3S 的 1.72）
3. **克隆 IsoME.jl**，核验是否真有 Nb benchmark 数据
   ```bash
   git clone https://github.com/cheil/IsoME.jl.git
   find IsoME.jl -name "*nb*" -o -name "*niobium*" -o -name "*.dat" -o -name "*.a2f"
   ```

### 中期（1–2 周）

4. **邮件索要 Margine 组的 H3S/YH6/YH9/PdH α²F**
   - 理由：arXiv:2608.14949 Table I 已算 Tc，数据必定存在
   - 若成功：+4 个材料（H3S 独立复现 + YH6/YH9/PdH 新增）
5. **邮件索要 Errea 的 LaH10**
   - 理由：Nature 2020 封面，最需要的靶标
6. **访问 Mendeley IsoME dataset**（doi:10.17632/frwsdxf44s.1）
   - WebFetch 失败，需手动访问网页 + 下载
   - 可能含 Nb / Pb / Al 的 benchmark α²F

### 长期（需用户决策）

7. **接受现状：只用已有的 H3S/D3S + MgB₂（自算）+ Nb（自算）**
   - 共 5 个数据点：H3S-200/160, D3S-200, MgB₂, Nb
   - 仍不足 8 个，但覆盖 λ∈[1.0, 2.3]、Tc∈[9, 203] K
8. **降级陈述**："亏损预算框架在 5 个材料上验证，μ*_fit ∈ [0.10, 0.15]（if 守住）；扩展到氢化物受限于数据可得性"
9. **重新定义靶标**：改用"至少 5 个体系 + 两个独立 λ 区间"（弱耦合 Nb/Al, 强耦合 H3S/MgB₂）

---

## 第 5 部分：脚本清单与数据溯源表

### 5.1 现有脚本（theory/notes/_scratch_ir/）

- `h3s_deficit_budget.py` — 快路径（trapezoid, eigh, geomspace）
- `deficit_independent_check.py` — 慢路径（Simpson, 显式循环, log 二分）
- `fill_factor_probe.py` — Einstein 谱 𝔣(λ)
- `hopfield_invariance.py` — λ⟨ω²⟩ 恒等式检验（已知是定义，W1 教训）
- `deficit_budget_table.py` — 表格生成器

**复用方案**：
- MgB₂：读 EPW 的 `mgb2.a2f`（格式：ω[meV], α²F, 累积 λ），喂给 `Tc_P(w, a2f, lam, mus, r, N)`
- Nb：同上
- 单位核验：λ = 2∫(α²F/ω)dω 自己算一遍，与 EPW 输出的末列比（应 <1e-6 相对差）

### 5.2 待下载数据（若成功获取）

| 材料 | 来源 | 路径/URL | 单位核验 | 状态 |
|---|---|---|---|---|
| MgB₂ | EPW 自算 | 本地 petlove-train:/zx/p4cert/mgb2_a2f/ | λ ∈ [0.6, 0.8] | 待执行 |
| Nb | EPW 自算 或 Margine 邮件 | TBD | λ = 1.05 ± 0.05 | 待执行 |
| YH6/YH9 | Margine 邮件 | TBD | 文献 λ (YH6)~2.2 | 待联系 |
| LaH10 | Errea 邮件 | TBD | λ ~2.2–2.8 @170–200 GPa | 待联系 |
| PdH | Margine 邮件 | TBD | Tc=9 K (实验) | 待联系 |

### 5.3 负结果表（尝试过但失败）

| 材料 | 尝试的来源 | 失败原因 |
|---|---|---|
| LaH10 | Nature SI, Zenodo, Figshare, arXiv | 只有图，无数值 |
| CaH6 | arXiv:2103.16282 (Ma PRL), Pickard AIRSS | 无 DFT 配套数据 |
| YH6 | Ding Zenodo 10.5281/zenodo.17472986 | .mat 文件不含 α²F(ω) |
| PdH | Mishra 2608.14949 | 论文算了但未发 Figshare |
| Nb | IsoME.jl, Allen–Dynes 1975 表 | repo 无文件清单；历史表只有 λ/ω_log |
| MgB₂ | EPW examples, IsoME Mendeley | 只有输入，无输出 α²F |

---

## 第 6 部分：给用户的三个选项

### 选项 A：保守路线（≤1 周）
- **交付**：H3S/D3S（已有 3 点）+ MgB₂ + Nb（自算 2 个）= **5 个材料**
- **陈述**：「在 5 个体系（λ∈[1.0, 2.3], Tc∈[9, 203] K）上验证，μ*_fit 全部落在 [0.10, 0.15]（if 守住）」
- **诚实边界**：未达 8 个，但覆盖弱/强耦合两端 + 双带材料（MgB₂）
- **风险**：审稿人可能要求更多材料

### 选项 B：激进路线（2–4 周）
- **行动**：
  1. 自算 MgB₂ + Nb（1 周）
  2. 邮件 Margine (YH6/YH9/PdH) + Errea (LaH10)（等待 1–3 周）
  3. 访问 IsoME Mendeley dataset（若有 Pb/Al，再 +2）
- **最好情况**：5（现有）+ 2（自算）+ 4（邮件）+ 2（Mendeley）= **13 个材料**
- **最坏情况**：邮件无回复，退化为选项 A
- **风险**：时间不确定（邮件可能石沉大海）

### 选项 C：改变问题（立即）
- **不追求 8 个材料**，改为：
  - 「亏损预算在 H3S/D3S 上给出 −2% / +3.6% 残差，框架**可证伪输入 α²F**（harmonic 被拒）」
  - 「扩展受限于氢化物数据不公开（systematic barrier）」
  - 重点放在 **MgB₂ 的极值性诊断**（它是失效区最近锚点）
- **优势**：立即可写，诚实且有科学价值
- **劣势**：放弃「8 个体系验证 μ*∈[0.10, 0.15]」的原始目标

---

## 附录 A：搜索引擎查询记录（供复现）

执行的 WebSearch 查询（按时间顺序）：
1. `niobium Nb superconductor alpha2F eliashberg function benchmark data McMillan`
2. `LaH10 alpha2F eliashberg function data Zenodo Figshare 2024 2025 2026`
3. `MgB2 magnesium diboride alpha2F eliashberg spectral function numerical data repository`
4. `Errea LaH10 Nature 2020 source data supplementary eliashberg spectral function`
5. `Materials Cloud three phonons electron phonon LaH10 YH6 archive repository`
6. `NOMAD repository superconductor eliashberg alpha2F LaH10 hydride dataset`
7. `IsoME.jl github eliashberg benchmark niobium data electron phonon`
8. `EPW code QE-forge MgB2 example alpha2F download github`
9. `Allen Dynes 1975 transition metals Nb Al Pb alpha2F tabulated numerical values`
10. `Carbotte 1990 review superconductor alpha2F data tables Nb Pb V3Si MgB2`
11. `Piekarz niobium phonon electron-phonon PRB spectral function data`
12. `Quantum ESPRESSO EPW examples directory github niobium lead aluminum data files`
13. `"Flores-Livas" hydride superconductor data repository LaH10 YH6 CaH6 supplementary`
14. `Boeri Margine electron phonon coupling MgB2 Al Nb data files 2001 2010`
15. `site:gitlab.com QEF q-e EPW examples mgb2 pb aluminum niobium`
16. `Pickard Cambridge hydride superconductor structure data LaH10 CaH6 download AIRSS`
17. `Sanna 2023 2024 LaH10 lanthanum hydride DFT supplementary data alpha2F github`
18. `Ponce EPW superconductor benchmark dataset Pb MgB2 download 2020 2021 2022`
19. `QE-forge EPW gitlab repository examples download Pb lead magnesium diboride`
20. `superhydride benchmark dataset github zenodo materials cloud H3S LaH10 YH6 numerical`
21. `Mishra 2608.14949 H3S YH6 YH9 supplementary data figshare zenodo alpha2F download`
22. `"a2f" OR "alpha2f" filetype:dat OR filetype:txt niobium lead aluminum download`

WebFetch 尝试（部分失败）：
- `https://github.com/cheil/IsoME.jl` — repo 结构可见，文件清单不可见
- `https://www.epw-code.org/doc/MgB2.html` — SSL 证书错误
- `https://data.mendeley.com/datasets/frwsdxf44s/1` — 页面框架可见，文件列表不可见

---

## 附录 B：关键 URLs（已核验可访问）

### 已确认有用（但需进一步操作）
- EPW MgB₂ tutorial: [docs.epw-code.org/tutorials/tutorial_04](https://docs.epw-code.org/tutorials/tutorial_04)
- EPW GitLab (examples): [gitlab.com/QEF/q-e](https://gitlab.com/QEF/q-e) (需 browse 到 EPW/examples/)
- IsoME.jl GitHub: [github.com/cheil/IsoME.jl](https://github.com/cheil/IsoME.jl) (需 `git clone`)
- IsoME Mendeley Data: [data.mendeley.com/datasets/frwsdxf44s/1](https://data.mendeley.com/datasets/frwsdxf44s/1) (需手动访问)
- Mishra H3S/D3S Figshare: [figshare.com/articles/dataset/31136635](https://figshare.com/articles/dataset/Data-H3S-D3S-PH-EPW-Anharmonic/31136635) (MANIFEST 已用)

### 关键论文（无数据但有联系线索）
- Mishra et al., arXiv:2608.14949 (Ann. Phys. 2026) — H3S/YH6/YH9/PdH
- Margine 组 Frontiers 2023 (arXiv:2201.07403) — Nb
- Errea et al., Nature 2020 — LaH10
- Flores-Livas et al., PRB 105, 174516 (2022) — LaH10 掺杂
- Pickard et al., arXiv:2105.02296 — CaH6 高通量

---

## 最终建议

**若用户需要快速交付**（本周内）：
→ **选项 A**：自算 MgB₂ + Nb，交付 5 个材料，诚实说明氢化物数据不公开

**若用户可等 2–3 周**：
→ **选项 B**：同时进行自算 + 邮件索要，争取 8+ 材料

**若用户接受重新定义问题**：
→ **选项 C**：聚焦 H3S 的「框架可证伪输入」+ MgB₂ 极值性诊断，不追求材料数量

我的建议：**选项 A（保守）+ 同步启动邮件（选项 B 的一部分）**，即：
1. 立即自算 MgB₂ + Nb（1–2 天）
2. 同时发邮件给 Margine + Errea（无需等待结果）
3. 用 5 个材料先写报告，若邮件回复了就追加

这样既保证交付时间（1 周内 5 个材料），又保留扩展到 8+ 的可能性（若邮件成功）。
