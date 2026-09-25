# Radevych 2025 给出 Pap1977 η 的外部见证（7/7），但**关不掉** RMTA 方法偏差

**日期**：2026-09-25 · **等级**：见每条 · **状态**：live
**一句话**：新取得的 arXiv:2508.21281 TABLE I 逐元素复现了我们从 Pap1977 页图读出的
`η_tot(I)`（7/7，`η_tot(II)` 0/7），并由**作者自己的话**给出方向——旧 η **偏低** ⟹
我们的 `0/31` 是保守 PASS；但该文**全程只比较 RMTA 的不同实现**，
⟹ 「过渡金属 RMTA 相对真值的方法偏差」**仍然未定**。
**§5 追加（同日）**：Savrasov PRL 72, 372 也已从零读完——**过渡金属侧零个数，只有散文**；
唯一干净的方法偏差数是 sp 金属 Al 的「RMTA 低估 3.14×」。两篇合计仍关不掉该缺口。

复核脚本：[`../_scratch_ir/_chk_main_radevych2025_vs_pap1977_eta_2026-09-25.py`](../_scratch_ir/_chk_main_radevych2025_vs_pap1977_eta_2026-09-25.py)
（`ALL CHECKS PASSED`；不 import 任何既有脚本）

---

## 0. 取件与仪器

> D. Radevych, T. Shishidou, M. Weinert, E. R. Margine, A. N. Kolmogorov, I. I. Mazin,
> *Rigid muffin-tin approximation in plane-wave codes*, arXiv:**2508.21281**
> （= npj Comput. Mater. `s41524-026-02141-7`）。

- 取得：`curl -s -L https://arxiv.org/pdf/2508.21281` → HTTP 200，15 页，
  本地 `library/pdf/eta_hunt_2026-09-25/arxiv_2508.21281.pdf`。**未用 `-k`。**
- **两台仪器都跑了**：
  文本层 `pdftotext -enc UTF-8 -layout` → [`../_scratch_ir/_radevych2508_utf8.txt`](../_scratch_ir/_radevych2508_utf8.txt)；
  页图 `pdftoppm -png -r 230 -f 8 -l 8` → `library/pdf/eta_hunt_2026-09-25/radevych_p8-08.png`。
  两者在 TABLE I 的 7 行 × 6 列上**逐字一致**，无 OCR 分歧。

> ### ⚠ 顺带订正一条仪器记忆（与本笔记的物理判词无关，但承重）
> **旧陈述**：记忆 `arxiv-and-nature-are-tls-intercepted-on-this-machine` 记 arxiv.org
> 被 TLS 劫持（证书 subject 为空、对该主机名无效）。
> **本次实测**：`curl -s -o /dev/null -w '%{http_code} %{ssl_verify_result}' https://arxiv.org/`
> 给 **`200 0`**——`ssl_verify_result=0` 即**校验通过**，PDF 正常到手。
> **两者都不是错的**：旧观测在当时为真，本次状态已变（或此前走的是另一条路径）。
> ⟹ 判词应为「**该记忆已过期，须逐会话现测**」，不是「旧观测是假的」。
> **禁令不变**：任何情况下不得用 `-k` / `--insecure` 绕过校验取承重原文。
> 同轮实测：`nature.com` 与 `raw.githubusercontent.com` **连不上**（HTTP 000，
> 但 `ssl_verify_result=0` ⟹ **是连接失败，不是证书错**，两者不可混为一谈）。

---

## 1. TABLE I（p.8）逐字

表注逐字：

> `TABLE I. Estimated semi-phenomenological EPC constant, λ, based on calculated
> McMillan-Hopfield parameters, η, at specified MT-radii, r_MT, and experimental Debye
> temperature, θ_D, for simple metals. θ_D values are taken from corresponding references.
> MH parameters η_ref, calculated with the **APW method and MJW potentials [14]**, are
> provided for the reference.`

参考文献 [14] 逐字 = `D. A. Papaconstantopoulos, L. L. Boyer, B. M. Klein, A. R. Williams,
V. L. Morruzzi, J. F. Janak` —— **正是我们 32 元素语料的来源** PRB 15, 4221 (1977)。

| Structure | r_MT (bohr) | η (eV/Å²) | η_ref [14] (eV/Å²) | θ_D (K) | λ |
|---|---|---|---|---|---|
| V (BCC) | 2.48 | 8.01 | 6.89 | 385 | 1.19 |
| Cu (FCC) | 2.41 | 0.57 | 0.46 | 343 | 0.09 |
| Nb (BCC) | 2.70 | 9.07 | 7.63 | 275 | 1.45 |
| Mo (BCC) | 2.57 | 5.56 | 5.80 | 430 | 0.35 |
| Rh (FCC) | 2.56 | 4.95 | 4.72 | 480 | 0.24 |
| Pd (FCC) | 2.66 | 2.65 | 1.99 | 299 | 0.31 |
| Ag (FCC) | 2.73 | 0.28 | 0.29 | 225 | 0.06 |

**单位与我们一致**（eV/Å²，页图上是 `eV/Å²` 的分式排版，文本层也给出 `eV` / `Å2` 两行）。

---

## 2. ★ 收获一：对我们抽取 Pap1977 的**外部见证**，7/7

我们此前的风险登记（见 [`tier0_eta_ceiling_32_elements_2026-09-25.md`](tier0_eta_ceiling_32_elements_2026-09-25.md)）：
Pap1977 的 OCR 把 `η_tot(I)` 列标签印成 `&tot(II)`，**列归属只能靠页图定死**。
这是一个**单点、无第二见证**的读数——正是最该被外部核的地方。

| 元素 | 我们读的 `η_tot(I)` | 我们读的 `η_tot(II)` | Radevych `η_ref[14]` | \|I−ref\| | \|II−ref\| |
|---|---|---|---|---|---|
| V | 6.894 | 4.370 | 6.89 | 0.004 | 2.520 |
| Cu | 0.456 | 0.380 | 0.46 | 0.004 | 0.080 |
| Nb | 7.627 | 5.410 | 7.63 | 0.003 | 2.220 |
| Mo | 5.797 | 3.440 | 5.80 | 0.003 | 2.360 |
| Rh | 4.722 | 2.490 | 4.72 | 0.002 | 2.230 |
| Pd | 1.992 | 1.050 | 1.99 | 0.002 | 0.940 |
| Ag | 0.288 | 0.250 | 0.29 | 0.002 | 0.040 |

容差取 0.005（他们印 3 位有效数字 ⟹ 末位半个单位）。

**`η_tot(I)` 命中 7/7；`η_tot(II)` 命中 0/7。**

⟹ **我们读的确实是 (I) 列，OCR 的标签错位没有污染数值。**
这是一条**真见证**而非恒等式：两端是**两条独立的生成路径**（我们的页图 OCR vs
他们 2025 年重新排印的表），且判别器本身有分辨力——`η_tot(II)` 被 0/7 排除。

> ⚠ **这 7 行不是 7 个独立的物理检验**，是**同一次抽取动作的 7 个采样点**。
> 它证的命题只有一个：「我们读对了列」。不得计入天花板的通过率。

**等级**：`numerical-evidence`（两端都是他人发表的表格数值）。

---

## 3. ★ 收获二：作者自陈的方向 ⟹ 我们的 `0/31` 是**保守** PASS

正文 Table I 讨论段逐字：

> `It should be noted that the MH parameters from Ref. [14] were obtained with the augmented
> plane wave (APW) method and Moruzzi-Janak-Williams (MJW) LDA potentials [72]. The noticeable
> differences from the present results are likely due to the evolution of potentials and the use
> of different exchange-correlation functionals. **While the prior MH parameters tend to be lower
> than those obtained in this work**, the overall trend across elemental metals remains similar.`

**方向审计**：`Tc_ceiling = 137.0959·√(η/M)`，`r = Tc_exp/Tc_ceiling ∝ 1/√η`。
旧 η **偏低** ⟹ 天花板**偏低** ⟹ 判据**更难**通过 ⟹ 由单边界方向规则，
**我们在偏低的 η 上取得 PASS，蕴含在更高的新 η 上同样 PASS。**

这是作者自己的话，不是我们的推断——比我们自己论证「RMTA 低估 ⟹ 方向安全」多一层外部支撑。

逐元素的漂移（**仅供定位，见 §4 的限定**）：

| 元素 | η_new/η_ref | r 的乘子 `1/√比值` | 方向 |
|---|---|---|---|
| Pd | 1.332 | 0.867 | r 降 |
| Cu | 1.239 | 0.898 | r 降 |
| **Nb** | **1.189** | **0.917** | r 降 |
| **V** | **1.163** | **0.927** | r 降 |
| Rh | 1.049 | 0.976 | r 降 |
| Ag | 0.966 | 1.018 | ⚠ r 升 |
| Mo | 0.959 | 1.021 | ⚠ r 升 |

**7 个里 5 降 2 升，最不利是 Mo（r 最多 +2.1%）。**
包封：元素栏最紧点 In `r=0.4460` × 最坏乘子 `1.0214` = **0.4555 < 1** ⟹ 仍 PASS。

A15 的主元素 **Nb 与 V 的漂移方向都使 r 下降**（乘子 0.917 / 0.927）。
若把元素 Nb 的漂移**指示性地**搬到 Nb₃Ge：`0.626 × 0.917 = 0.574`。
> ⚠ **这只是指示，不是新判决。** 元素 Nb 的势漂移不能严格外推到化合物 Nb₃Ge
> （不同晶体结构、不同 MT 半径、不同 semicore 处理）。它给的是方向，不是数。
> **不得把 0.574 写进任何表。** 承重的仍是 `0.626`。

---

## 4. ⚠ 收获三（否定的）：该文**关不掉** RMTA 的方法偏差——我派活时的定位是错的

**旧陈述错在哪**：我（主会话）派 `eta-second-method` 时，依据摘要转述把这篇定位成
「RMTA vs FLAPW 的方法误差测量」，并据此认为它能关掉「过渡金属 RMTA 是否高估 η」的缺口。
**这个定位是错的**，且错在一个从名字推语义的地方——`FLAPW` 是**基组/势的实现**，
不是「非 RMTA 的正确计算」。

**新陈述凭什么对**（逐字证据三条，全部在手）：
1. TABLE I 的 `η` 与 `η_ref` **两列都是 RMTA**：`η` 是他们的 PS+PW(赝势平面波) RMTA，
   `η_ref` 是 Pap1977 的 APW RMTA。表注已逐字引在 §1。
2. FIG. 6 / FIG. 7 的 `PS+PW` vs `FLAPW` 对照，**FLAPW 侧跑的同样是 RMTA**——
   FIG. 7 图注逐字：`Comparison of partial, η_{l,l+1}, and total, η, type-resolved MH
   parameters between PS+PW and FLAPW methods`。比的是 **MH 参数**，即 RMTA 的产物。
3. 作者结语逐字（补充材料 §III 末）：
   `All results are correct within the frameworks they are derived from.`

⟹ **全文从未把 RMTA 与一个非 RMTA 的 EPC 计算对照。**
摘要里的 `excellent agreement` 指的是**实现一致性**（同一近似的两种数值实现），
**不是方法正确性**。

⟹ **「过渡金属 RMTA 相对真值的方法偏差」这个缺口，本轮未关闭，仍是「未定」。**
现有唯一的干净方法误差数仍是 Al 的 `η_RMTA/η_true = 0.14/0.44 = 0.318`
（见 [`savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md`](savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md) §7），
那是 **sp 金属**；过渡金属侧至今只有 Savrasov 正文的一句散文（`complete agreement` / `close`）。

**预算（留给下一轮）**：要把 Nb₃Ge 的 `r=0.626` 推过 1，需 RMTA 对 Nb 基**高估** η 达
`(1/0.626)² = 2.55×`。现有证据都指向该因子接近 1（且 §3 的方向是反的）。
⟹ 这个缺口**大概率不翻判决**，但它目前只由散文支撑 ⟹ 值得一次定向取数，不值得大规模投入。

---

## 5. Savrasov PRL 72, 372 (1994) 已从零读完 —— 过渡金属侧**零个数**

teammate `eta-second-method` 报了一张 §5 表，称由该文的 λ 得出 RMTA 的 η 方法偏差。
**我自己把 PDF 和页图都读了**，三处缺陷，与它的结论存活与否分开写。

复核脚本：[`../_scratch_ir/_chk_main_savrasov_prl72_lambda_direction_2026-09-25.py`](../_scratch_ir/_chk_main_savrasov_prl72_lambda_direction_2026-09-25.py)
（`ALL CHECKS PASSED`；不 import 任何既有脚本）
页图：`library/pdf/eta_hunt_2026-09-25/savrasov_prl72_tableI-3.png`（p.374，`pdftoppm -r 250`）
文本层：[`../_scratch_ir/_savrasov_prl72_utf8.txt`](../_scratch_ir/_savrasov_prl72_utf8.txt)

### 5.1 TABLE I 逐字（列头**只有三个**）

文本层是 OCR 产物，把 `λ` 印成 `A`/`p`（行标签 `pcalc`/`pexp`/`pctarlc` 实为
`λ^calc`/`λ^exp`/`λ_tr^calc`）⟹ **数值与列归属以页图为准**，两仪器在数值上逐字一致。

| 行 | Al | Nb | Mo |
|---|---|---|---|
| λ^calc | 0.44 | 1.26 (1.44, 1.37, 1.08) | 0.42 |
| λ^exp | 0.42 | 1.33 (1.71, 1.43, 1.10), 1.04 | 0.44 |
| λ_tr^calc | 0.36 | 1.17 | 0.35 |

### 5.2 三处缺陷（旧陈述错在哪）

1. ~~**`λ_Savr(V) = 1.19` 不存在于该文**——表里没有 V 行。1.19 恰等于 §1 中
   Radevych TABLE I 里 **V 的 λ** ⟹ 串源。该行作废。~~
   **⟹ 这条判词是我（主会话）的假订正，同日已撤回，见 §5.5。**
   活着的那一半只有：`PRL 72,372` 的表**确实没有 V 行**（Nb 1.26 / Mo 0.42 是对的）。
2. **方向反了**（单边界方向规则的靶心）：`r ∝ 1/√η`，η 变大 ⟹ r 变**小**。
   比值 `λ_LR/λ_RMTA > 1` ⟹ 正确操作是 **除以** `√ratio`：
   Nb `0.626/√1.0311 = 0.6165`（它写 `×` 得 0.636）、Mo `0.626/√1.25 = 0.5599`（它写 0.700）。
3. **λ 比值混了 ⟨ω²⟩**：`λ = η/(M⟨ω²⟩)` ⟹
   `λ_Savr/λ_Pap = (η_Savr/η_Pap)·(⟨ω²⟩_Pap/⟨ω²⟩_Savr)`。
   Pap1977 取 `⟨ω²⟩ = Θ_D²/2`（德拜代餐），Savrasov 用算出的声子谱。
   我们只观测到**乘积**；把它整个归给 η 使误差变成**双边**
   ⟹ 单边界豁免**不适用**（同 `sandwich-kernel-budget-must-be-quadratic-in-fxc` 形状）。

### 5.3 它的 PASS 存活（新陈述凭什么对）

判据是 `r ≤ 1`。方向反的结果恰好**高估** r（0.636 > 0.6165；0.700 > 0.5599），
而上界量的**高估只支持 PASS、永不支持 FAIL** ⟹ **PASS 结论存活，垮掉的是理由**，
且正确方向下余量比它报的更大。
⚠ 但因缺陷 3，这两个数**本来就不该**被当作 η 方法偏差的证据消费。

### 5.4 ★ 核心否定发现：过渡金属侧该文给出的 RMTA 数值个数 = **0**

唯一干净的方法偏差数是 **Al**，正文 p.374 左栏逐字：

> `In order to check previous conclusions [7,22] about the inapplicability of the RMTA
> for sp metals, we also performed such a calculation and indeed found λ = 0.14.`

对上同表 `λ^calc(Al) = 0.44` ⟹ `λ_RMTA/λ_LR = 0.318`，RMTA **低估 3.14×**。
它干净是因为**同一篇、同一套势、同一个 ⟨ω²⟩ 约定** ⟹ 缺陷 3 的混杂不出现。
但它是 **sp 金属**，且作者明说这正是 RMTA 的已知失效域。

过渡金属侧只有两句**散文，无数**：

> Nb（p.374 左栏）：`We obtain similar results with the RMTA, in full agreement with
> earlier calculations [17,24]`
> Mo（p.374 右栏）：`Our linear-response calculations of α²F(ω) are found to be close
> to our RMTA calculations and to earlier ones [30]`

⟹ **「过渡金属 RMTA 方法偏差」判词仍是「未定」**，本轮两篇（Radevych §4 + 本篇）都关不掉。
**等级**：`numerical-evidence`（「某论文印了什么」的事实认定）。

**预算复算**：要把 Nb₃Ge 的 `r = 0.626` 推过 1，需 RMTA 对 Nb 基**高估** η 达 `2.552×`。
已知因子：Al 是**低估** 3.14×（方向相反）、Radevych 代际漂移最坏 `1.021×`（Nb/V 使 r 下降）。
⟹ 该缺口**大概率不翻判决**，但只由散文支撑 ⟹ 不得写「已关闭」。
**要关它，需要的是同一篇论文里同时印出 `η_RMTA` 与 `η_全势` 的表**
（或 λ 比值 + 两端 `⟨ω²⟩` 约定同时给出）——λ 比值这条路因缺陷 3 已堵死。

### 5.5 ⚠ 我自己的一次假订正：「V 行串源」已撤回

**旧陈述错在哪**（§5.2 第 1 条，同日我自己写下）：
我写「`λ_Savr(V)=1.19` 恰等于 Radevych TABLE I 里 V 的 λ ⟹ 串源、该行作废」。
**这条是假的。** 我拿「两个数相等」当了「同一个来源」——
而 CLAUDE.md W1 热区那条恰恰是：**两个数一致时，先问「它们是同一个证书/同一个约定下的量吗」，
再问「谁错了」**。我把它用反了：用巧合去证出处，且**只查了一个候选出处就下判词**，
没有先 grep 仓库里 1.19 还有没有别的来源。
这是 W1 的原型形状——**注意力全在「对方错在哪」，新陈述被免检放行**（本轮第 N 起，同一天里第二起）。

**新陈述凭什么对**：真出处是
**Savrasov & Savrasov, PRB 54, 16487 (1996) TABLE II**（期刊页 **16 494**），
`λ_calc` 行逐字（我本轮**亲自读页图** `library/pdf/eta_hunt_2026-09-25/savr_rmta_p8-08.png`）：

> `λ_calc:  Al 0.44 | Pb 1.68 | V 1.19 | Nb 1.26 | Ta 0.86 | Mo 0.42 | Cu 0.14 | Pd 0.35`

表注逐字：`TABLE II. Comparison between the calculated electron-phonon coupling constants
λ_calc and the values of λ_tun deduced from the tunneling experiments. Also listed are the
values of λ_{s−h} extracted from the measured specific-heat coefficient γ with the use of our
calculated density of states N(ε_F).`
⟹ **V 行不作废**，它是已发表值。

⚠ **一个必须钉死的编号区分**：`PRB 54, 16470` 与 `PRB 54, 16487` 是**两篇不同的论文**
（同卷相邻）。16470 = `Linear-response theory and lattice dynamics: A muffin-tin-orbital approach`，
只有 Si + NbC 的声子频率两张表；16487 = 金属的电声耦合那篇，TABLE II 在此。
我此前把它们当同一篇的不同页，是这次假订正的诱因之一。

**V 行真正的缺陷**（这才是该说的那条）：交付方把**三篇不同论文**的数字揉进一张标题写
「本轮三篇 PDF」的表，**却没有逐行标出处** ⟹ 读者无法验证归属。
交付方已自陈并接受删除该行的处置。另外 §5.2 第 3 条的 `⟨ω²⟩` 混杂对 V 行**同样适用**
⟹ 它在任何情况下都不承重。

**⟹ 净变化**：§5.4 的核心否定发现（过渡金属侧 RMTA 数值个数 = 0）与 §5.3 的 PASS 存活判决
**都不受影响**——V 行从来不是 RMTA-vs-非RMTA 的对照，它是 λ_calc(线性响应) 对 λ_Pap(RMTA+德拜)，
仍然带 `⟨ω²⟩` 混杂。改变的只有「那个数有没有出处」。

### 5.6 ★ 第三起「实现一致 ≠ 方法正确」——并换来一个更强的头条理由

复核脚本：[`../_scratch_ir/_chk_main_rmta_nb_reliability_retraction_2026-09-25.py`](../_scratch_ir/_chk_main_rmta_nb_reliability_retraction_2026-09-25.py)
（`ALL CHECKS PASSED`）

**旧陈述错在哪**：仓库/记忆此前记「RMTA 不是全局失效，作者明言 Nb 上
`performed our own RMTA calculations and obtained complete agreement`
⟹ 失效是 sp 金属专属 ⟹ A15 头条落在 RMTA 可靠侧」。
**那句的两端不是 RMTA vs 真值。** 上下文逐字（页图 `savr_rmta_p8-08.png`，期刊页 16 494）：

> `The discrepancy found by us has already been reported in the past **RMTA-based**
> calculations.[41,42] To check the consistency of our results with **the earlier ones**,
> we have performed our own RMTA calculations and obtained complete agreement between them.`

`the earlier ones` 在前一句被明确定义为 `the past RMTA-based calculations`
⟹ 比较集合 = {他们自己的 RMTA, 更早的 RMTA} ⟹ **实现 vs 实现**。
PRL 72,372 是同一句式：`in full agreement with **earlier calculations** [17,24]`。
⟹ **Nb 侧不存在 RMTA-vs-线性响应的一致性陈述。**
（Mo 那句 `Our **linear-response** calculations … close to our **RMTA** calculations` 两端**是**对的，
但散文无数 ⟹ 过渡金属侧被**量化**的跨方法对照个数 = **0**，与 §5.4 一致。）

⚠ 这是**同一轮里第三起**同形事故：① Radevych 的 FLAPW-RMTA vs PS+PW-RMTA（§4）；
② PRL 72,372 的 `full agreement with earlier calculations`；③ 本条。
**`agreement` 这个词在 RMTA 文献里默认指实现一致**，读到它必须先定位两端。

**新陈述凭什么对**：同一段里有一个此前被漏掉的定量把手，**两篇独立论文各自给出同一区间**：

> PRB 54,16487（文本层第 576–577 行 + 页图）：
> `the RMTA calculations give the values of λ varying from **1.12 to 1.86**`
> PRL 72,372（文本层第 252 行 + 页图）：
> `a large variety of λ (**1.12–1.86**) have been obtained`

线性响应（同作者同代码）`λ_calc(Nb) = 1.26` ⟹ 最大**高估**包封 `1.86/1.26 = 1.4762×`。

**方向审计（与 Al 那条相反，必须分开写）**：`r ∝ 1/√η`。
RMTA **高估** η 达 κ ⟹ 真 η 更小 ⟹ 真天花板更低 ⟹ 真 `r = r_报 × √κ`（**上升** = 危险方向）。
Al 那条是**低估**（r 下降 = 安全）。两条方向相反，**不得共用一个「保守 PASS」判词**。

> **`r(Nb₃Ge) ≤ 0.626 × √1.4762 = 0.7606 < 1`**，余量 `1.315×`。
> 翻判决需 `κ = (1/0.626)² = 2.5518×`，包封比它小 `1.729` 倍。

⟹ **头条结论存活，理由被换成更强的一条**：从「**假定** RMTA 在 Nb 上准」
变为「**即使**它差到文献区间上端也仍 PASS」——更弱的前提、更强的结论。
**等级**：`numerical-evidence`。
⚠ 限定：`1.12–1.86` 是**区间宽度**（跨不同实现/不同势），不是系统偏差；
且 λ→η 仍带 `⟨ω²⟩` 混杂 ⟹ 本条是**指示性包封**，不是认证界。

---

## 6. 未取得 / 未定

| 项 | 状态 | 说明 |
|---|---|---|
| 过渡金属 RMTA 方法偏差（V/Nb/Mo 的 η_RMTA/η_true） | **未定** | Radevych 关不掉（§4）、Savrasov 只有散文（§5.4） |
| Savrasov PRB 54, 16470 有无 η 表 | **已读，无 η 表** | 交付方逐字：正文标题 `Linear-response theory and lattice dynamics: A muffin-tin-orbital approach`、`IV. APPLICATIONS` 子节 `A. Si` / `B. NbC`、TABLE I/II 均为 `phonon frequencies ... (THz)`；全文仅两张表 ⟹ 无金属超导数据 |
| Savrasov PRB 54, **16487** 有无 η 表 | **已读，无 η 表** | 与 16470 是**两篇不同论文**；全文 4 表、无 `Hopfield` 一词；TABLE II 给 λ_calc/λ_tun/λ_{s−h}/N(ε_F)/γ（V 的 1.19 出自此处，见 §5.5）|
| Savrasov PRL 72, 372 有无 η 表 | **已读，无 η 表** | 全文只有 λ/ρ/w，无 η 数值（§5） |
| 补充材料 [1] / [2] 的身份 | **已定**（逐字到手） | 交付方从 `_radevych2508_utf8.txt` 第 1103–1109 行独立重读：`[1] Papaconstantopoulos, Boyer, Klein, Williams, Morruzzi, Janak, Calculations of the superconducting properties of 32 metals with Z ≤ 49, Phys. Rev. B 15, 4221 (1977)`；`[2] W. E. Pickett, Transferability and the electron-phonon interaction: A reinterpretation of the rigid-muffin-tin approximation, Phys. Rev. B 25, 745 (1982)` ⟹ [2] 标题自陈是 **RMTA 的再解释**，**不是**非 RMTA 对照 |
| 补充材料 §III FIG. 2 的 A15 数值 | **未取得（判定为不可取）** | 交付方逐字：Y 轴仅刻度 `η (eV/Å2)` 0.0/2.5/5.0/7.5/10.0/12.5，X 轴类别 `V3Sn V3Si V Mo Nb Nb3Si NbN`，图例 `PS+PW / FLAPW / Papaconstantopoulos et al. / Pickett`，**柱子无数字标注** ⟹ 不目测填数 |
| Pickett PRB 25, 745 的 A15 η | **未取得** | 标题已表明是 RMTA 的再解释；若与我们 Klein 1979 的 A15 冲突，需单独审 |

---

**相关**：
[`tier0_eta_ceiling_32_elements_2026-09-25.md`](tier0_eta_ceiling_32_elements_2026-09-25.md)（被本笔记加固的那条判决）、
[`savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md`](savrasov_vs_pap1977_crossmethod_eta_2026-09-25.md)（Al 的 3.14× 方法偏差）、
[`../literature/eta_numerical_tables_not_published_2026-09-25.md`](../literature/eta_numerical_tables_not_published_2026-09-25.md)（§5 已被同日订正）。
