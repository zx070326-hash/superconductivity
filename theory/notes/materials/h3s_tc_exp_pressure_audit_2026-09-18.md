# H3S 实验 Tc 压强赋值审计（2026-09-18）

> 任务来源：team-lead 派活，核「亏损预算」headline 对照（H3S 199.0 K vs 「实验 203 K」）背后的压强赋值是否自洽。
> 方法：全仓 grep `203`（.md/.csv），逐处溯源；**两篇一手文献均已下载并用 Read 工具逐页读图**（非 WebFetch/pdftotext）：
> - Drozdov et al., Nature 525, 73 (2015)，`arxiv.org/pdf/1506.08190` → `_scratch_ir/drozdov2015.pdf`（6 页，全文）。
> - Errea et al., PRL 114, 157004 (2015)，APS harvest 全文 PDF → `_scratch_ir/errea2015.pdf`（5 页，全文，JBIG2 扫描件路线）。
> 严格性标注：文献转录标 **L1**（本会话亲自读页图/全文逐字）；仓库内部转引标**转引**；无来源支持的数字标**推断**。
> 两个 PDF 是数据文件不是脚本，未 `git add -f`（毕业规则只管 `_scratch_ir/*.py` 等被笔记点名当证据的**代码**；本文档直接点名的是可复核的 arXiv/APS URL，不依赖这两个本地缓存文件）。
> **2026-09-18 二轮更新**：team-lead 复核后指出判定 3 第 0 版（190 K 锚）依赖的 Errea Table II `T_c(expt.)` 列源自被 Drozdov 正式版点名撤回的未退火预印本数据（arXiv:1412.0460）。本轮用 `fitz`（PyMuPDF）直接解析 Drozdov (2015) Nature 论文 Fig. 2(c) 内嵌的矢量绘图坐标（散点标记本身是矢量图形，非位图），重新数字化出干净的退火后压强匹配锚点，取代有争议的 190 K。脚本 `_scratch_ir/digitize_fig2c.py` 已 `git add -f`，`sh tools/check_cited_scripts_tracked.sh` PASS。判定 3 与 §D 已按两段式订正改写，旧版「+4.7%」结论已废弃，见判定 3 正文。

## A. 赋值表：每一处「H3S 实验 Tc」出现的地方

| 文件:行号 | 声称压强 | Tc 值 | 溯源 | 直读/推断 |
|---|---|---|---|---|
| `theory/data/hydrides_experimental.csv:2` | 155 GPa | 203 K | Drozdov 2015, Nature 525:73 | **L1 直读**（见下方 §C.1） |
| `theory/data/hydride_tc_record_envelope.csv:2` | 155 GPa | 203(1) K onset | Drozdov 2015 Fig.4(a) 原句转录 | **L1 直读**（笔记自带原句，本次复核逐字确认） |
| `theory/notes/isotope/gap_ratio_from_eigenvector_2026-09-17.md:66,152` | 155 GPa（§References 处写明） | 203 K, 2Δ=60 meV | Du/Drozdov/Minkov et al. 2025 Nature（隧道谱论文，doi:10.1038/s41586-025-08895-2） | 转引（本次未独立核对 2025 隧道谱原文，但其压强/Tc 与 2015 磁化率一手值一致，无内部冲突） |
| `theory/notes/materials/deficit_budget_h3s_2026-09-17.md:36`（H3S-200 anharm 行） | 理论侧 200 GPa，Tc_exp 列写 203 | 199.0 vs 203，残差 −2.0% | **文件自身已标注压强不匹配**（第 172 行原文） | 见 §B 判定 1 |
| `theory/notes/materials/deficit_budget_h3s_2026-09-17.md:40`（H3S-160 anharm 行） | 理论侧 160 GPa，Tc_exp 列**同样写 203** | 211.2 vs 203，残差 +4.0% | 同上 | 见 §B 判定 2（这是 team-lead 记忆中「203 同时赋给两行」的出处） |
| `theory/notes/materials/deficit_budget_certificate_2026-09-17.md:136,201` | 未标注理论侧压强的行内文本 Tc_exp=203，但 §5 诚实边界（201行）已写明「实验 Tc=203 K 在 155 GPa」vs「α²F 在 200 GPa」 | 203 K | 同 §C.1 | 转引（该文件本身已履行 CLAUDE.md 纪律，主动标注不匹配） |
| `theory/notes/materials/eliashberg_crosscheck_2026-09-17.md:20` | 未显式标压强，行文「对应实验 203/150 K」 | 203 / 150 K | 转引自 hydrides_experimental.csv | 转引 |
| `theory/notes/materials/meta/session_resume_verification_2026-09-17.md`（`theory/notes/meta/`） | 未显式标压强 | 203/150 | 转引 | 转引 |
| `theory/notes/materials/materials_expansion_2026-09-17.md:166,214` | 只作为 Tc 值域端点罗列「Tc∈[9,203]K」 | 203 K | 转引 hydrides_experimental.csv | 转引，非配对压强的判据用途，不受本审计影响 |
| `theory/data/hydrides_theory.csv`（H3S-250 harmonic 行） | 250 GPa | Tc=226 K | Errea PRL 2015 Table II | **本次升格为 L1 直读**（原标「主会话推断，非直读」，见 §C.2） |
| `theory/data/hydrides_theory.csv`（H3S-250 anharmonic 行） | 250 GPa | 原写 `not-reported` | Errea PRL 2015 Table II / 正文 | **发现错误，见 §B 判定 3** |

## B. 判定

### 判定 1（确认，非新发现）：`deficit_budget_h3s_2026-09-17.md` 与 `deficit_budget_certificate_2026-09-17.md` 均已自陈压强不匹配

**旧陈述**（team-lead 派活时的疑虑）：亏损预算笔记拿「实验 203 K」去对照理论侧 200 GPa 的计算，二者压强不一致，这个对照可能在悄悄比较不同压强的量而未声明。

**核验结果**：这条疑虑**部分属实但已被原作者自己声明**——两份笔记都各自写了一行「压强不匹配」的诚实边界（`deficit_budget_h3s_2026-09-17.md:172` 与 `deficit_budget_certificate_2026-09-17.md:201`），并不是隐藏缺陷。**级别：analytic**（这是文本核对，不是数值结论）。

**但**：这两条声明本身只说「不匹配，但同量级、未消除」，**没有触及本审计发现的更严重问题**——见判定 3。

### 判定 2：203 K 被同时用作 200 GPa 行和 160 GPa 行的 Tc_exp，是同一实验值被两次复用，不是「两个不同的 203」冲突

`deficit_budget_h3s_2026-09-17.md` 第 36 行（H3S-200 anharm）和第 40 行（H3S-160 anharm）的 Tc_exp 列都写 203。核对后这**不是数据输入错误**：全仓只有一个「一手、磁化率确认」的 H3S 实验 Tc 值（155 GPa / 203(1) K），该文件把它同时拿来对照 200 GPa 和 160 GPa 两套理论计算，是因为**没有其他压强点的一手实验 Tc 可用**。文件第 137 行也写明理由（「H3S-160 行压强更接近」）。

**判定**：不是内部冲突，是「唯一实验锚点被迫复用」——但这恰恰说明**没有一行理论计算的压强真正等于 155 GPa**，见判定 3。

### 判定 3（两轮订正——本节按 CLAUDE.md「订正分两段写」纪律记录完整过程）

**第 0 版陈述**（本审计初稿）：203 K 不是 200 GPa 的合适对照值；Errea (2015) 自己在 200 GPa 用的实验对照是 190 K；据此把 H3S-200anharm 的残差从 −2.0%（对 203 K）改判为 +4.7%（对 190 K）。

**第 0 版错在哪**（由 team-lead 指出，本会话独立复核后确认属实）：190 K 出自 Errea 论文 Table II 的 `T_c(expt.)` 列，其来源引用是 **arXiv:1412.0460**（Drozdov/Eremets/Troyan 更早预印本）。Drozdov et al. (2015) **正式发表版** Fig. 2a 图注原句（page 9，本次重读逐字确认）：

> 「(Note here that the only point for sulfur deuteride presented in **arXiv:1412.0460** was determined **without sample annealing**, and Tc would increase after annealing at room temperature.)」

同一图注给出退火前后的量级对照：145 GPa 处未退火 Tc=170 K，室温退火后升至 195 K（+14.7%）。这说明 arXiv:1412.0460 里的数据点**系统性偏低**（至少 D3S 那个点如此），Errea 引用的「190 K」继承自这份有偏数据源，**不能当作干净的压强匹配锚**——它同时存在压强匹配（好）与制备方法有偏（坏）两个问题，与 203 K 的「制备方法干净、压强不匹配」构成两个独立缺陷的对照，**残差符号在两者之间被夹住、不能简单地判成 +4.7%**。第 0 版把「压强匹配」当成唯一判据，漏检了「样品退火状态」这第二个判据，是真实的分析缺陷。

**第 1 版陈述凭什么对**（本会话独立验证，用**矢量位图数字化**而非肉眼读图——见方法说明）：

Drozdov (2015) Nature 论文本身的 **Fig. 2(c)**（图注原句：「The data on annealed samples are presented.」——这批数据本来就是退火后的，不依赖有争议的 arXiv:1412.0460）给出了压强匹配 **且** 退火干净的第三个锚点。本次没有用肉眼估读，而是用 `fitz`（PyMuPDF）直接解析 PDF 内嵌的矢量绘图坐标（散点标记是矢量矩形/圆，而非位图像素），把标记中心坐标按坐标轴刻度线（同样是矢量线段，量出的刻度间距）转换成 (P, Tc)。**校准的内部自检**：数字化流程独立给出的洋红色方块（Fig. 4(a) 磁化率点在 Fig. 2c 上的标注）坐标为 **(155.0 GPa, 203.0 K)**，与图注原文「155 GPa, T_onset=203(1) K」**精确重合**（零误差），证明该数字化流程本身的系统误差可忽略（脚本：`_scratch_ir/digitize_fig2c.py`，已 `git add -f`）。

**H3S（黑色实心方块，annealed，电阻法）在 P≈200 GPa 附近的读数**（全部标记，按压强排序）：

| P (GPa) | 145.0 | 150.9 | 175.0 | 176.0 | 176.9 | 183.0 | 183.0 | 190.0 | **195.0** | **200.0** | 206.0 | 220.0 | 226.0 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Tc (K) | 195.4 | 190.0 | 190.0 | 179.0 | 185.0 | 176.0 | 185.0 | 176.4 | **190.0** | **184.0** | 171.4 | 171.4 | 170.0 |

**P=200 GPa 处直接命中一个标记：Tc = 184.0 K**（不需要插值）。邻近点（195 GPa→190 K，206 GPa→171.4 K）给出局部散布 ~±10–15 K，反映样品间/批次间的实验波动（Drozdov 原文多次强调结果对制样条件敏感），不是数字化误差。

**D3S（红色空心圆，annealed）在 P≈200 GPa 处**：200.0 GPa → 145/146 K（两个相邻标记）。这与仓库 `hydrides_experimental.csv` 记的 ~150 K **一致**（差 3%），与 Errea Table II 的 `T_c(expt.)=90 K` **不一致**（差 60 K，62%）——**独立确认了 team-lead 的判断：Errea 那一列的 90 K 就是被 Drozdov 正式版点名撤回的未退火赝值，`hydrides_experimental.csv` 的 ~150 K 才是干净值**。

**结论（重新表述，符号已可判定）**：

- 干净的压强匹配锚：**Tc(H3S, 200 GPa, annealed, resistivity) = 184.0 K**（Fig. 2c 矢量数字化，本次校准自检零误差；局部散布 ±10–15 K）。
- H3S-200anharm 理论预测 199.0 K vs 184.0 K ⟹ **残差 = +8.2%（高估）**，若计入局部散布（171.4–190 K 范围），残差落在 **+4.7% 到 +16.1% 之间，符号稳定为高估**（不像第 0 版那样只有一个点估计，也不像 team-lead 指出问题后那样「符号未定」——三个可能锚点 171.4/184.0/190 全部给出正残差）。
- 原仓库用 203 K 算出的 −2.0%（低估）**确实是压强不匹配导致的假象**：203 K 是 155 GPa 处的穹顶峰值，200 GPa 处（无论退火与否）的干净读数都明显低于它。

**独立复核（team-lead，2026-09-18，自零重取，未读/未引用本会话的 `digitize_fig2c.py`）**：team-lead 自写 `_scratch_ir/_chk_drozdov_fig2c_independent.py`，改用刻度线段（非文字标签 bbox）定标，次刻度自检 5/5 通过，洋红已知点复现 (154.98 GPa, 203.01 K)（误差 −0.02 GPa/+0.01 K），200 GPa 处直接命中同类矩形图元 **184.04 K** ⟹ 残差 **+8.13%**，与本节 184.0 K/+8.2% 在两位小数内一致，方法完全独立（不同定标基准、不同脚本）。稳健区间口径有一处分歧并已查明原因：team-lead 排除了多段线（星形/三角）图元的承重资格（其 bbox 中心不等于数据点，只有矩形/贝塞尔闭合图元的 bbox 中心才精确等于数据点），得到承重区间 **+4.74%~+17.03%**（略宽于本节 +4.7%~+16.1%，因保留了 225.97 GPa 处 170.04 K 这一精确圆心点）；两个区间的**下端一致、符号结论一致**（194–227 GPa 窗内三个精确图元 189.99/184.04/170.04 K 全部 <199 K，不翻号）。三点边界随结论同行：①199.0 K 是仓库既有输出，两次核验均未重算；②该锚点"干净"的理由是 Fig. 2c 图注逐字声明"annealed samples"，不是数值上更合意；③CSV 与 `deficit_budget_h3s_2026-09-17.md` 均未改动，等用户发话。

**级别**：190/226 K（Table II 数字）与「样品未退火」的文字证据是 **L1（原文逐字/逐格）**；Fig. 2c 的 184.0 K 是 **L1（矢量坐标直接解析，非位图估读，校准自检零误差——比通常的「读图」更强，接近「rational-exact 式」的坐标读出，但因刻度线本身仍是有限精度的印刷坐标、且标记半径引入 ±1–2 GPa/±2 K 的量化不确定度，仍标 L1 而非升格）**；「残差 +8.2%，正负号在扰动下不变」是 **numerical-evidence**（199.0 K 未重新独立算过，继承自仓库既有 float64 结果；除法本身是 analytic）。

### 判定 4（新发现，纠正 `hydrides_theory.csv`）：H3S-250 GPa anharmonic 行的 Tc 不是 `not-reported`，而是 190 K；226 K 确属 harmonic（已确认无误）

**旧陈述**（`theory/data/hydrides_theory.csv` 现状 + `deficit_budget_certificate_2026-09-17.md` 记录）：
- 「226 K 的归属为主会话推断（非直读）……需原图复核」；且怀疑曾被 teammate 错误归给非谐。
- H3S-250 GPa anharmonic 行 Tc 标 `not-reported`（notes 写「Tc^ME,anh 在 pdftotext 中缺失」）。

**新陈述凭什么对**（L1 直读 Table II，见上表转录）：
1. **226 K 确实属于 harmonic**（T_c^ME,har 列，250 GPa 行）——原「推断」判定现已由本次读图**证实**，可升格为 **L1**，无需再打「推断」标签。之前 pdftotext 缺数是路线问题（扫描 PDF 文本层不可靠），不是数字本身有误。
2. **250 GPa anharmonic 的 Tc 并非缺失，是 190 K**——Table II 该行第七列（T_c^ME,anh）清楚印着 190。CSV 现状 `not-reported` 是错的，应改为 190。这与判定 3 引用的正文原句「anharmonic 时 190K 几乎不随压强变化（200–250 GPa 区间）」互相印证，双重确认。

**行动建议**（未执行，因不 commit/改文件超出被授权范围——见下方「未做的事」）：`theory/data/hydrides_theory.csv` 对应行的 `Tc_predicted_K` 应从 `not-reported` 改为 `190`，`notes` 列中「Tc^ME,anh 在 pdftotext 中缺失」应改为「Tc^ME,anh=190 K，已 L1 读图确认（`h3s_tc_exp_pressure_audit_2026-09-18.md`）」。

## C. 250 GPa / 226 K 原图复核（细节记录）

### C.1 Drozdov 2015 原句核对

Fig. 4 图注（page 11，L1 读图）：「Temperature dependence of the magnetization of sulfur hydride **at a pressure of 155 GPa** in the zero-field cooled (ZFC) and 20 Oe field cooled (FC) modes. The onset temperature is **T_onset = 203(1) K**.」——与仓库两处 csv 记录的「155 GPa / 203(1) K」逐字一致，**确认无误，可从「转引」升格为「L1 直读」**（本次核验之前，仓库内是转引 team-lead 之前一轮工作的读图结果；本次是主会话独立重读同一份 PDF，逐字复现)。

摘要另称「confirmed by magnetic susceptibility measurements with Tc=203 K」（page 1）——未在摘要级出现压强，压强信息只在 Fig.4 图注给出，**必须读图注才能确定 155 GPa**，仅读摘要/正文会漏掉压强（这正是本审计要防的错误类型）。

### C.2 Errea PRL 2015 Table II 逐格转录（page 3，L1）

已完整转录于判定 3。补充核对：
- Table I（harmonic-only, McMillan 方程估计，page 2）另有三行（C2/c HS2、C2/m HS2、C2/m HS，均非 H3S/Im-3m 相），与 H3S Table II 无关，仓库未引用，不受影响。
- 正文（page 3）：「the anharmonic level ... λ=1.84 ... 30% smaller than the harmonic result」——与 Table II 的 λ^har=2.64/λ^anh=1.84 一致。
- 正文（page 4）：「using the smaller value of μ*=0.1, we obtain Tc=222 K」——这是 H3S 200 GPa anharmonic 在 μ*=0.1（而非表中默认 μ*=0.16）下的另一个数，仓库当前未收录此行，**若未来要用需新增一行并标明 μ*=0.1**。

## D. 受影响的下游结论

1. **`deficit_budget_h3s_2026-09-17.md` 的 headline「μ*=0.13 给出 199.0 K，对应实验 203 K，残差 −2.0%」**——该表述本身没有算错（199.0 vs 203 确实是 −2.0%），但**「203」不是与该理论计算压强匹配的实验值**这件事，比该文件已承认的「不主导误差」更严重。判定 3 的最终结论（Fig. 2c 矢量数字化，非 Errea 转引的有偏 190 K）：换成压强匹配 **且** 退火干净的锚点 184.0 K，残差从「低估 2.0%」变为**稳定的「高估 8.2%」（局部散布下 +4.7% 到 +16.1%，符号在此范围内不变号）**。**这直接影响「μ*_fit=0.122 落在经典窗口」这一结论的解释**——μ*_fit 是反解自 Tc_exp=203 K 得到的，若改用 184.0 K 反解，μ*_fit 会系统性变小（未重新计算反解值，仅指出该链条现在依赖一个压强不匹配的输入，且新锚点的方向是「需要更小的 μ*_fit 才能压低理论 Tc 到位」）。
2. **同位素指数** §4（`deficit_budget_h3s_2026-09-17.md:91-109`）：「实验两点比 203/150 → α=0.4365」——150 K 是 D3S 在 ~150 GPa（同一篇 Drozdov 2015 报道），203 K 是 H3S 在 155 GPa；这两个数已经是同一篇文献内部对照，压强大致匹配（155 vs ~150 GPa），**不受判定 3 影响**（判定 3 只影响「203 K vs 200 GPa 理论」这一条链，不影响「203 K vs 150 K 两点实验对比」这一条链）。补充：本轮 Fig. 2c 数字化顺带给出 D3S(200 GPa, annealed)=145–146 K，与仓库 ~150 K 一致（差 3%），进一步支持仓库这条链本身没有问题；但独立证实了 **Errea Table II 的 D3S T_c(expt)=90 K 不能用作任何锚点**（同源于被撤回的未退火数据，偏差达 60 K/62%，比 H3S 侧的偏差大得多——D3S 点正是 Drozdov 图注点名"未退火"的那一个）。
3. **`deficit_budget_certificate_2026-09-17.md`**（证书化基础设施笔记）：生产运行尚未执行（阻塞于服务器权限），**Tc_exp: 203 K** 只是预期输出模板里的占位对照值。判定 3 的结论应在生产运行前一并采纳，即认证目标压强应与 184.0 K（200 GPa、annealed、Fig. 2c 矢量数字化）而非 203 K（155 GPa 穹顶峰值）核对，或者干脆放弃把「实验残差」当作认证目标，只认证 Tc 本身的区间。
4. **`hydrides_theory.csv` 数据表本身**：H3S-250GPa anharmonic 行需要改（判定 4：`not-reported`→190），226 K 的「推断/需复核」标记应解除（升格为 L1）。
5. **`hydrides_experimental.csv` 的 D3S 行**（team-lead 要求追加）：应补一条 notes，明确指出「不要把 Errea Table II 的 D3S T_c(expt)=90 K 当作实验锚点——该数据源自 arXiv:1412.0460 的未退火测量，已被 Drozdov 正式发表版 Fig. 2a 图注点名撤回；仓库当前记录的 ~150 K 已被本次 Fig. 2c 矢量数字化独立复核（145–146 K，annealed，200 GPa）」。本审计**未直接编辑该 CSV**（权限边界，见下节），此处只记录建议文字供 team-lead/用户采纳。

## E. 未做的事（诚实边界）

- ~~未独立核实 D3S ~150 GPa 那条实验数据的图注压强~~ **本轮已部分补上**：Fig. 2c 矢量数字化顺带给出 D3S(200 GPa, annealed)=145–146 K，与仓库 ~150 K 一致，但这是**读图交叉核验**、不是重读该值原始出处的图注文字（未再单独确认 ~150 K 这个数字最初来自哪一张图/哪一句话，只确认了它与本次独立读出的坐标一致）。
- 未核实 2025 年 Du/Drozdov/Minkov 隧道谱论文（doi:10.1038/s41586-025-08895-2）原文中 203 K 的压强图注，只确认了它与 2015 论文的数字一致、无内部冲突；若要升格为 L1 需单独下载核对。
- 未修改 `theory/data/hydrides_theory.csv`、`theory/data/hydrides_experimental.csv`、`deficit_budget_h3s_2026-09-17.md`、`deficit_budget_certificate_2026-09-17.md` 本身——按纪律「不要 commit/改配置」，且这些文件的修改权限不在本次派活范围内，留给 team-lead/用户决定是否采纳判定 3/4 与 §D5 的建议文字并编辑。
- 未重新反解 μ*_fit（若改用 184.0 K 作为 H3S-200 的实验基准）——这需要跑 `deficit_independent_check.py` 或等价代码，本审计只是数据溯源，未获派算这一步。
- Fig. 2c 数字化的刻度标定依赖两组矢量刻度线坐标（各 2 点线性插值），未对全部刻度线做冗余交叉核验（只用了自检点 155/203 GPa/K 精确重合这一个证据支撑标定无系统偏差）；局部散布 ±10–15 K 的成因（批次/样品差异 vs 标记读数量化误差）未拆分。

---

**总结（给 team-lead 的一句话，已按两段式订正更新）**：203 K 的溯源本身没有错（155 GPa、Drozdov 2015 Fig.4a，L1 确认）；但第 0 版审计把 Errea 转引的 190 K 当作「干净的压强匹配锚」是错的——那个数继承自被 Drozdov 正式发表版点名撤回的未退火预印本（team-lead 抓出的问题，本会话独立复核确认属实）。本轮改用 `fitz` 直接解析 Drozdov Fig. 2c 的内嵌矢量坐标（非肉眼读图，校准自检零误差），读出干净的退火后锚点 **184.0 K @200 GPa**，并顺带证实 D3S(200GPa,annealed)=145–146K 与仓库 ~150K 一致、反证了 Errea 的 D3S 90K 同样不可用。结论：headline 残差从 −2.0%（对 203K，低估）变为**稳定的 +8.2%（高估，局部散布下 +4.7%~+16.1%，符号不翻转）**，需要 team-lead 决定是否用 184.0 K 重算亏损预算表，以及是否在 `hydrides_experimental.csv` 的 D3S 行加上 §D5 建议的撤回说明。
