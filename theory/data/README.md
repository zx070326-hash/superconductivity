# theory/data · 数据表说明

## ⚠ 2026-09-08 二轮更新：`isotope_curated.md` 取代本目录 CSV 中的同位素相关行

researcher 读本地 PDF 原文后发现一轮摘要级表有两处实质错误：
- **LaD10 Tc=180 K / α=0.46** 不在 Drozdov 2019 原文，是 Durajski 2020 理论值经两篇文献误转引；原文为 α=0.57（fcc）/0.62（P4/nmm），LaD10 Tc ~168 或 ~160 K
- **YD9 无实验数据**（Kong 2021 只做 YH6/YD6）；唯一来源是已撤稿论文
- H3S 作者给的 α_H ≈ 0.3（区间比较），不是两点比 0.44

同位素相关的量一律以 `isotope_curated.md` 为准。

## hydrides_theory.csv

来源：researcher agent，2026-09-17。收集氢化物理论计算参数：电声耦合常数 λ、特征频率 ω_log 和 ω̄₂（Allen-Dynes 1975 第二矩）。

### 列说明

- **material, isotope, pressure_GPa**：材料标识，与 `hydrides_experimental.csv` 对应
- **treatment**：计算方法，`harmonic`（谐波 DFT 声子）/ `anharmonic`（SSCHA/TDEP 等）/ `anharmonic+nonadiabatic`（含顶角修正）
- **lambda, lambda_uncertainty**：电声耦合常数 λ = 2∫[α²F(ω)/ω]dω
- **omega_log_meV, omega_log_uncertainty_meV**：对数平均频率 ω_log = exp⟨ln ω⟩，单位 meV
- **omega_2_meV, omega_2_uncertainty_meV**：第二矩频率 ω̄₂ = √⟨ω²⟩，其中 ω̄₂² = (2/λ)∫α²F(ω)·ω dω（Allen-Dynes 1975 定义），单位 meV
- **mustar_used**：计算中使用的库仑赝势 μ*
- **Tc_predicted_K, Tc_uncertainty_K**：理论预测 Tc
- **evidence_tier**：证据等级（2026-09-17 重建时新增列），见下方「标记约定」
- **source_author, source_year, source_venue, source_identifier, source_table_figure**：数据来源
- **notes**：说明

### 为什么 ω̄₂ 重要

本仓库的严格 Tc 天花板（**方向 A「普适速度极限」，见 `theory/derivations/P4_universal_tc_envelope.md`；`theory/derivations/INDEX.md` §4.4 / §5-A。不是 KAY 2016 的结果，也不是定理 E——定理 E 是 Einstein 极值性**）为：

**Tc ≤ C_∞ · √λ · ω̄₂**，  C_∞ = (1/2π)√𝔤(2) = 0.1827262477 [analytic]

只要有 (λ, ω̄₂) 两个数，就能给任何材料算天花板，**不需要完整 α²F(ω)**。
（注意方向：ω_log ≤ ω̄₂ 恒成立，所以拿 ω_log 顶替 ω̄₂ 会**低估**天花板 ⟹ 是上界的下估 ⟹ 只能支持 PASS、永不支持 FAIL。见 auto-memory `p4-onesided-bound-direction-rule`。）

ω̄₂ 恰好是 Allen-Dynes 强耦合修正因子 f₂ 要用的量，因此大量论文本来就计算它（常记作 ω̄₂、ω_2、ω₂、⟨ω²⟩^{1/2}）。

### 校准点

前两行（H3S 200 GPa）是本仓库已有完整 α²F 计算的校准点：
- **Harmonic**：λ=2.2929, ω_log=107.165 meV, ω̄₂=129.32 meV
- **Anharmonic**：λ=1.7210, ω_log=109.995 meV, ω̄₂=149.26 meV

从文献收集的 H3S 数据应与此一致。

### 数据完整性

> ✅ **2026-09-17 已重建**（脚本 `theory/notes/_scratch_ir/_build_hydrides_theory.py`，单位换算由脚本计算）。
> 旧版几乎整列 `not-reported`，因为生成它的 agent 只检索了**论文摘要**
> （表内多条 notes 原文即「未在摘要中直接给出」）。ω̄₂ / ω_log 从来不写在摘要里，
> 它们在表格和方法节里。该 agent 同时忽略了 `library/pdf/` 下的 111 篇 PDF
> （其中就有 Errea 2015 H3S、Errea 2020 LaH10、Ding 2025 YH6/YH9）和
> `surveys/a2F_machine_readable_sources_2026-09-08.md`。
>
> **重建后 23 行**：ω_log 18/23 有值，**ω̄₂ 只有 7/23**。ω̄₂ 稀缺是**真实的**，
> 但原因与旧版所说的不同——不是「文献不报道」，而是**高压旗舰体系的论文报 ω_log 不报 ω̄₂**
> （H3S/LaH10 都要自己从 α²F 积分）。唯一直接印出 ω̄₂ 的是 Ding 2025 的 YH6/YH9，
> 它反过来**不报 ω_log**。另两行 ω̄₂ 来自本仓库自己的 α²F 积分。
>
> **反证**：Cerqueira–Sanna–Marques, arXiv:2307.10728（库内 `2024_Tiago_F_T_Cerqueira_2307.10728.pdf`）
> 的 Eq.(1c) 就是 ⟨ω²⟩ = (2/λ)∫α²F(ω)ω dω。Xie 2022（`2106.05235`）以
> **ω̄₂/ω_log = 1.1 和 1.6** 为界画阴影带。⟹ ω̄₂ 是这一行的标准报道量。
>
> **眼下可直接用的 (λ, ω̄₂) 样本**：姊妹论文 **arXiv:2403.13496**（库内已有 PDF）
> TABLE I + II 逐材料印出 (Tc_AD, ω_log[K], ω₂[K], λ, E_hull)，约 48 个常压氢化物候选，
> μ*=0.10。主会话已逐行转录并抽样复核。
> （**订正**：arXiv:2307.10728 那 6912 个材料的**公开数据集只含 λ 和 ω_log**，
> 不含逐材料 ω₂；ω₂ 需自行从 Materials Cloud 的 raw `a2F.dos`（~48 GiB）积分，
> 属于需要用户批准的资源动作。本条曾被写成「6912 个逐个给出 ω₂」，过强。）
>
> 源头调查见 `theory/notes/materials/omega2_data_sources_2026-09-17.md`；
> 详见 auto-memory `omega2-not-missing-cerqueira-6912-2026-09-17`。
>
> 检索提示：`pdftotext` 会把希腊字母/上下标糊掉（ω→`w` 或乱码，ω̄₂→`w?2`，ω_log→`wlog`），
> **用 Unicode 符号搜 PDF 文本会静默返回零命中**。搜 PDF 一律用 ASCII 词或物理术语。

**重建后的覆盖（2026-09-17）**：

| 体系 | 行数 | λ | ω_log | ω̄₂ | 档 |
|---|---|---|---|---|---|
| H3S 200 GPa（本仓库 α²F） | 2 | ✓ | ✓ | **✓** | repo |
| H3S/D3S（Errea PRL 2015） | 6 | ✓ | ✓ | ✗ | L1-self |
| LaH10/LaD10（Errea Nature 2020） | 7 | ✓ | ✓ | ✗ | L1-self |
| YH6/YH9（Ding 2025） | 5 | ✓ | ✗ | **✓** | L1-self |
| CaH6（Liang 2021） | 3 | ✓ | ✓ | ✗ | WebFetch |

**主会话对 teammate 转录的两处订正**（均据本地 PDF 全文）：
1. H3S 200 GPa 谐波 ω_log = **90.4** meV，不是 90.1（原文两处独立出现 90.4）。
2. H3S 250 GPa 的 226 K 归属为**谐波**、非非谐——摘要称谐波 Tc 随压力下降、非谐 Tc
   几乎与压力无关（200 GPa 时非谐为 194.0）。该行 pdftotext 缺一个数，**仍需原图复核**。

**注意：同一个「H3S 200 GPa harmonic」有两套不相容的数**——本仓库 α²F 给 λ=2.2929 / ω_log=107.165 meV，
Errea PRL 2015 给 λ=2.64 / ω_log=90.4 meV。不同的计算，不是转录错误；**引用时必须指明哪一套**。

**标记约定**：
- `evidence_tier`：`repo`（本仓库自有 α²F 积分）/ `L1-self`（主会话亲自读原文文本核对）/
  `WebFetch`（仅网页摘要转述，**最弱，承重使用前必须补核**）
- `not-reported`：**该文献确实没有报道这个量**（已读过表格/正文，非检索失败）
- `NA`：不适用或无不确定度信息

### 证据等级

所有文献数据标记为 **numerical-evidence**（按 `CLAUDE.md` §1 四级分类）。这是别人的 DFT 计算结果，不是 analytic 或 certificate 级别。

### 与实验表的关系

`hydrides_experimental.csv` 提供实验 Tc 和同位素系数 α；本表提供理论输入参数。两表材料清单应一致，以便交叉检验理论预测。

## ambient_hydride_candidates.csv

来源：Cerqueira, Fang, Errea, Sanna, Marques, "Searching Materials Space for Hydride
Superconductors at Ambient Pressure", **arXiv:2403.13496** TABLE I + TABLE II（库内 PDF）。
48 行由主会话从 pdftotext 全文逐行转录并抽样复核（Li2CuH6 = 86.0 / 491 / 890 / 1.96）。
生成脚本 `theory/notes/_scratch_ir/_build_ambient_candidates.py`（K→meV 换算由脚本计算）。

用途：在**一整个材料族**上检验方向 A 的天花板 Tc ≤ C_∞·√λ·ω̄₂。

**结果**：𝔣 = Tc/(C_∞√λ·ω̄₂) ∈ [0.1334, 0.4384]，均值 0.2475，**违反 0/48**。
ω̄₂/ω_log ∈ [1.197, 3.358]，均值 1.831（比 Xie 2022 画的 1.1–1.6 示意带宽得多）。

> ⚠ **这批 Tc 是 Allen–Dynes 公式值（μ*=0.10），不是实验值。**
> 在这批数据上 𝔣 只检验「AD 公式与 P4 界相容」，**不检验自然**。
> 要拿实验 Tc 检验的点见 `hydrides_theory.csv` 的 YH6/YH9 行
> （𝔣_exp = 0.528–0.583）与 H3S（𝔣_exp ≈ 0.49）。**两件事不能混着报。**

> ⚠ `column_assignment = UNCERTAIN` 的 6 行（Perovskites ABH3）：pdftotext 把
> (ω₂, λ) 两列拆到了另一行，逐行配对是按顺序 + 正文旁证（"λ … reach 1.70 for KInH3"）复原的。
> 对该块另算了**最坏配对上界**（每个 Tc 配该块最小的 √λ·ω₂）：**𝔣 ≤ 0.9629**，
> 故 0/48 的判决对该块的任意列置换都成立。按单边界方向规则，这是 𝔣 的**上估**，
> **只支持 PASS、永不支持 FAIL**。

严格性：C_∞ = (1/2π)√𝔤(2) [analytic]；全部 𝔣 值 [numerical-evidence]（外部 DFT + AD 公式 + float64）。

## hydrides_experimental.csv

来源：researcher agent `hydride-exp-data`，2026-09-08，第 1/3 条回传。**原始产出，未手改。**

### 证据等级（agent 自述）
- 数据来自 WebSearch 摘要/二手转述，**未逐篇 WebFetch 原文核对**
- 数字置信度分"多源交叉"与"单源转述"，见各行 notes

### 主会话已补核（2026-09-08，WebFetch arXiv 摘要页）
- **LaH10 / LaD10**：Drozdov et al. arXiv:1812.01561（Nature 569, 528）摘要明确把"isotope effect"列为超导证据之一，所以 LaD10 数据的一手来源**就是这篇**，Durajski 2020 的转引有源头；具体 Tc(D) 与压力在正文图表，仍需读正文确认 180 K vs 169–170 K 哪个对。
- **YH9 / YD9**：Kong et al. arXiv:1909.10482（Nat. Commun. 12, 5075）摘要同样列出 isotope effect；YH9 243 K @201 GPa、YH6 227 K @237 GPa 确认；YD9 数字在正文。
- **YH6 / YD6**：Troyan et al. arXiv:1908.01534（Adv. Mater. 2021）摘要原文：Tc=224 K @166 GPa；**α=0.4**；Bc2(0)=116–158 T，比计算值大 2–2.5 倍；"calculations suggest unusually large impact of the Coulomb repulsion"——这句对 P1 的 δα_μ* 项是直接线索。
- **CaH6**：Ma et al. arXiv:2103.16282，215 K @172 GPa，WHH 外推 Hc2(0)=203 T；正式版 PRL 128, 167001 (2022)（DOI 按记忆填，未核）。

### 使用前必须处理的问题
1. **LaD10 的 α=0.46 溯源最弱**：数字出现在 Durajski et al. PRB 101, 214501 (2020) 理论论文的转引里，一手实验来源（Drozdov Nature 2019 正文是否含 LaD10）未确认。P1 用这个数之前要查原文。
2. **YH9 的 α=0.46–0.48 来自已撤稿论文**（Snider et al.），只能用 Kong et al. Nat. Commun. 2021 的 YH9/YD9（243 K @201 GPa / 172 K @205 GPa），α 需自己算：ln(243/172)/ln(9.0/... ) —— 注意 YH9 的质量比不是 2（Y 不变，只有 H→D），应按 α = −Δln Tc / Δln M_eff，M_eff 的定义本身是 P1 要处理的理论问题。
3. **CeD9 数据（67 K @111 GPa）与部分氘代阶梯（CeH8D 92 K、CeH6D3 82 K）出处未确认**，不建议定量使用，除非补核。
4. **H/D 样品压力不匹配**是所有体系的通病（H3S ~155 vs ~150；YH6 166/183 vs 173；YH9 201 vs 205）。P1 的 α 分解必须把 dTc/dP 的修正显式放进误差预算。
5. **三个同位素数据空白**：CaD6、ThD9/ThD10、LaBeD8 均无实验。P1 的检验集实际只有 H3S、LaH10（若溯源成立）、YH6、YH9、（CeH9 存疑）、PdH 五到六个体系。
6. H3S 130 GPa 处 α=2.37 是 Szczęśniak–Durajski 对已有 Tc(P) 的再分析，反映的是 R3m 相 + 结构相变区，不是单相内的同位素效应；应排除或单独处理。
7. Minkov et al. Nat. Commun. 2022 的 H3S/LaH10 磁化数据有公开可信度争议（arXiv:2409.12211）；本表未用其数据，但引用 Hc1/磁化时注意。

### 后续
- 第 2/3 条（数据质量说明 + 公式核对）、第 3/3 条（文献 JSON）到达后追加到本目录
- 补核清单：LaD10 一手来源；CeD9 出处；Ma et al. 2022 PRL 的 DOI；Chen et al. CeH9 正式版卷期
