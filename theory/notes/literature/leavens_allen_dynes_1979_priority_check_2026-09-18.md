# Leavens–Allen–Dynes, SSC 30, 595 (1979) — 优先权风险核查

**任务**：核查 C.R. Leavens, P.B. Allen, R.C. Dynes, *Solid State Communications* **30**, 595 (1979)
（下称 LAD79）是否已经给出与本仓库 `C_∞=0.1827262477…` 上界（`P4_g2_analytic.md`）/
Einstein 谱全 λ 极值性（`P4_einstein_extremal.md`）同等或更强的严格结果，从而冲击仓库
「首次严格证明」的对外措辞。同时顺路核查两篇背景文章
（Leavens SSC 17,1499(1975)；Leavens SSC 19,395(1976)）与 Owen–Scalapino Physica 55,691(1971)
的可得性。

**证据等级约定（本笔记专用，不是仓库四档）**：
- **一手** = 亲自读到该文原文（全文或摘要）的逐字内容
- **一手-元数据** = 未读到摘要/正文，但独立数据库（Crossref/OSTI/Semantic Scholar）给出的
  官方著录字段（标题、作者、参考文献列表）——这些是出版商登记的结构化数据，可信度高于检索引擎转述，
  但**不是**正文内容
- **二手** = 另一篇已读原文的论文，逐字引用或转述 LAD79 的具体陈述
- **三手/AI 转述** = WebSearch 工具自己合成的摘要段落，**未经我本人核对原文**，可信度最低，
  本笔记中一律显式标注，且不作为风险判决的独立依据

---

## 1. 检索过程记录

| 渠道 | 结果 |
|---|---|
| ScienceDirect（`sciencedirect.com/.../pii/0038109879901030`） | 403 Forbidden，连摘要都拿不到（WebFetch 直接被拒） |
| OSTI/ETDEWEB (`osti.gov/etdeweb/biblio/5841608`) | 完整著录字段（作者/期刊/页码/DOI），**无摘要正文、无全文 PDF** |
| NASA ADS (`ui.adsabs.harvard.edu/abs/1979SSCom..30..595L`) | 405 Method Not Allowed（WebFetch 打不进去） |
| Semantic Scholar API (`api.semanticscholar.org/graph/v1/paper/DOI:...`) | 返回记录但 **`abstract: null`**，`openAccessPdf.status: "CLOSED"`，
  disclaimer 明示「出版商已隐去 abstract 字段」 |
| Crossref API (`api.crossref.org/works/10.1016/0038-1098(79)90103-0`) | **拿到完整参考文献列表（12 条）**——这是本次最有价值的一手-元数据，见 §2.1 |
| P.B. Allen 个人主页 (`felix.physics.sunysb.edu/~allen/`，及其 `plist.html`) | curl 直接拿到完整页面（40KB「Favorite Publications」+ 36KB 主页）；
  **`grep -i "leavens\|1979"` 零命中**——LAD79 既不在他的「精选出版物」列表，也不在专门收录
  「难以获取论文」PDF 重印本的 `Pdffiles/` 区块（该区块列出的都是会议文集/书章节，
  不含这篇期刊短讯） |
| Google Scholar「所有版本」 | WebSearch 未能直接抓到该页面结构化数据，多次搜索均落回
  ScienceDirect/OSTI 同一批链接，无独立免费全文线索 |
| scholar.archive.org (fatcat) | 被反爬虫验证码墙挡住（ALTCHA challenge），未能检索 |
| ResearchGate | 未见该文（1979 年 Elsevier 短讯，作者当时均无 ResearchGate 习惯性自存档） |

**结论：全文未能获取。** 三个独立数据库（ScienceDirect 本身、OSTI、Semantic Scholar）
都确认没有免费全文；Allen 本人的个人主页——按用户建议专门查过——**没有**把这篇收进他的
「精选」或「难以获取论文」重印本列表，说明连作者本人的自存档习惯都没有覆盖到它。

---

## 2. 核实结果逐条

### 2.1 摘要/标题（一手-元数据，多源交叉确认一致）

标题：**"The superconducting transition temperature for very large electron-phonon coupling"**
（注意：不是任务背景里写的书名式全称，Crossref/OSTI/ScienceDirect 三处著录一致）。
摘要（OSTI/ScienceDirect 元数据一致）：
> "A rigorous derivation of the range of validity of the asymptotic relation
> k_B T_c ≈ 0.1827[λ⟨ω²⟩]^(1/2) is presented."

这句摘要本身**没有**出现在我读到的任何原始 PDF/HTML 页面里（我只在 WebSearch 综述和 OSTI/SD
元数据摘录中见到它），所以严格说这是**转引自出版商索引摘要字段**，不是我本人读到 PDF 里的
逐字段落；但它在 OSTI（美国能源部 INIS 提交渠道）与多个独立检索引擎摘要中反复一致出现，
可信度足以当作标题+摘要文本使用（列为「一手-元数据」而非「一手」）。

**读法上的关键分歧点**："range of validity of the asymptotic relation"——这句话字面上讨论的
是一个**已经写定的渐近关系式**（固定数值系数 0.1827）"在什么范围内成立"，这在数学结构上
**不同于**证明 `T_c ≤ 0.1827√(λ⟨ω²⟩)` 是对**任意**谱测度 P、**任意** λ 都成立的不等式
（后者才是本仓库 `P4_g2_analytic.md` + `P4_einstein_extremal.md` 证的东西，也是 KAY 论文
公开未决的 Conjecture 2）。"range of validity" 更自然的读法是：**这条渐近式在多大的 λ 范围内
数值上准确**，是对 Allen–Dynes 1975 自己 eq.(24)(25)（λ>10, μ*=0，由 N+1=64 阶截断数值
矩阵对角化"deduced"出系数 0.182）的**延伸/收紧**，回答"这条渐近式从多大的 λ 开始生效"，
而不是"证明它对所有形状 P 都是上界"。——**这只是读法，不是核实结论**，因为我没有拿到正文；
下面 §2.2–2.5 用间接证据检验这个读法。

### 2.2 参考文献列表给出的内容线索（一手-元数据，Crossref）

LAD79 的 12 条参考文献：McMillan PR167,331(1968)；Allen PRB12,905(1975)（即 Allen–Dynes 1975）；
Leavens SSC19,395(1976)；Louie SSC22,1(1977)；Wu 1977 *Scientia Sinica*；
**Owen & Scalapino, Physica 55,691(1971)**；Bergmann Z.Phys263,59(1973)；Scalapino 1969（书章节）；
Hopfield PR186,443(1969)；Poon PRB12,4816(1975)；Butler PRB14,3823(1976)；Leavens J.Phys.F7,1911(1977)。

这份列表全部是电子-声子强耦合超导 Tc 理论文献（McMillan/AD/Leavens 自己的前作、Owen-Scalapino
的配对不稳定性、Bergmann-Rainer 的敏感性分析），**没有**任何看起来专门处理"谱形状变分/极值性
证明工具"的条目（对比：AD 1975 自己在 Appendix A 用 Wheeler–Gordon 矩逼近定理证明逐元不等式，
若 LAD79 也做类似的形状极值证明，通常会引用相应的分析工具文献，这里没有出现）。
这与"这是一篇纯粹做大 λ 渐近展开/数值收敛区间讨论的短讯"的读法一致，**不构成决定性证据**，
只是弱支持。

### 2.3 二手引用：谁怎么用这篇文章

- **Allen & Mitrović, *Theory of Superconducting Tc*, Solid State Physics 37, 1 (1982/83)**——
  Allen 本人合著的权威综述。Crossref 参考文献列表核实**确实引用了 LAD79**（DOI 匹配，条目
  `bib key ≈ 90` 区间，同段落还引用了 Leavens SSC17,1499(1975)）。**但没有拿到正文**，无法读到
  Allen 自己如何描述 LAD79 的结果——这是本次调研最大的未闭合缺口（Allen 本人对自己论文的
  转述，本应是判定这篇文章意义的最佳二手来源）。
  - 仓库已有笔记 `theory/notes/literature/lit_einstein_extremality_novelty_2026-09-16.md` §10
    逐页读过 Carbotte RMP 62,1027(1990) 正文（另一位与 Allen 同时代的权威作者），其中 p.1042
    逐字（**该笔记内一手引文，本次未重复取件**）："There is, in fact, no limit to T_c in
    Eliashberg theory. … **Leavens (1975, 1977) and Leavens and Carbotte (1974, 1977)** arrived
    at the conclusion that T_c is not bounded in Eliashberg theory by using very different
    arguments from those discussed so far." ——注意这里点名的是 **Leavens 1975/1977 单独**
    的论文，**不含** 1979 年这篇三人合著的 LAD79，且讨论的是"T_c 无上界"这个**不同**的命题
    （与"渐近关系式的有效范围"不是一回事）。Carbotte RMP 同段 (2.29) 式把 **0.183√λω_E**
    这个渐近值的首次给出归功于 "Allen and Dynes (1975); Owen and Scalapino (1971)"——
    **依然不含 LAD79**。这是一个值得记录的负面证据：**同一位权威综述作者在专门讨论这个
    精确渐近系数的段落里，两次点名相关文献都没有点 LAD79**，倾向于支持"LAD79 不是该系数
    最初/最主要的严格来源"这一判断，但不能排除 Carbotte 只是漏引或认为不需要重复引用（毕竟
    LAD79 本身晚于 AD75 四年，作为"补充说明"类短讯常被综述省略）。

- **Natkunlaphat, Tasee & Pinsook, arXiv:2504.17435（2025，投给 Physica Scripta）**——
  一手全文（PDF 已下载、pdftotext 逐行核对）。正文第 134 行明确把 LAD79 与 AD75 并列引用
  （`[4, 6]`）：
  > "Another characteristic phonon frequency, ω̄₂, naturally arises in the asymptotic behavior
  > of T_c [4, 6]: T_c → C(μ*)√(λ⟨ω²⟩), (5) where C(μ*) is a function of the Coulomb
  > pseudopotential, μ*..."
  这是**一手**读到的现代二手转述（我亲自 pdftotext 核对了行号与上下文，非 AI 摘要）。
  它把 LAD79 定位为"渐近式里 μ*-函数系数 C(μ*) 的来源之一"，**间接确认 LAD79 确实处理了
  μ*>0 的情形**（否则不会被拉来和 AD75 一起支撑一个显式写着 `C(μ*)` 的公式）。
  但这仍是 2025 年一篇非权威来源（投稿中、非 Eliashberg 理论核心圈子作者）对 45 年前文献的
  笼统引用，**不能当作 LAD79 具体做了什么证明的权威转述**，只能当作"LAD79 涉及 μ*"这一点的
  弱佐证。

- **Kiessling–Altshuler–Yuzbashyan (KAY) I/II/III**（arXiv:2409.00533/00532/02121，
  2024–2025，J. Stat. Phys. 已发表）——**全文 grep 零命中 "Leavens"**（我下载了 PDF 并
  `pdftotext` 全文搜索 II 和 III，两篇都没有）。这三篇是目前对"C_∞=0.1827262477…这个常数
  的严格数学地位"讨论最细致、最挑剔的文献（KAY I §1.3 甚至专门点名 AD75 的某条无证断言
  "claims to the contrary in [AD] notwithstanding"，说明他们对这个领域历史文献做过审计式
  梳理）。**这三篇论文完全没有提到 LAD79**，这是一个有分量的负面证据：如果 LAD79 真的已经
  给出了 C_∞ 的严格证明或 μ*>0 下的严格推广，以 KAY 这种审计级别的严谨作者，大概率会引用
  并讨论（哪怕是为了指出它的缺口，就像他们对 AD75 做的那样）。他们没有引用，更符合
  "LAD79 是一篇处理数值渐近收敛区间的短讯，不构成严格不等式定理"这一判断。

  **附带发现（超出本任务范围，但直接相关，需单独跟进）**：KAY III 明确写道 C_∞ 常数
  "determined rigorously in the first one, (I)"——即 KAY 自称在他们 2024 年的 Part I
  论文里**已经严格确定**了 g(2)=C_∞ 这个谱半径常数。这与本仓库 `C_∞` 常数的"首次严格证明"
  措辞存在**潜在的独立优先权问题**，但这是**另一个**风险点（KAY vs 本仓库，不是 LAD79 vs
  本仓库），仓库 memory 里已有大量"KAY/EA"相关的历史核对记录（`kay-me-bias-is-cutoff-
  convention-2026-09-17.md` 等），说明这条线索**很可能已被主会话追踪过**；本笔记只是提醒
  这个交叉点仍然存在，不在此展开核实。

### 2.4 是否处理 μ*>0

**间接证据倾向于"是"**（§2.3 的 arXiv:2504.17435 引用），但**没有一手确认**具体做到什么
严格程度。任务背景关心的"仓库 μ*>0 深挖工作是否被 LAD79 抢先"——目前证据不足以下判断，
需要正文。

### 2.5 是否证明 Einstein 谱极值性

**没有任何证据（一手或二手）表明 LAD79 讨论谱形状极值性**。标题聚焦"very large
electron-phonon coupling"（大 λ 渐近区），这正是 AD75 自己已经指出"当 λ 很大时谱形状不再
影响 T_c"的区域（仓库既有笔记 `lit_einstein_extremality_novelty_2026-09-16.md` §1.1 逐字
引用 AD75 p.908："The asymptotic formula (22) indicates that **at very large λ, the shape
doesn't affect T_c**."）——也就是说，在 LAD79 讨论的这个极限区间里，Einstein 极值性问题
本身趋于**平凡**（一切谱形状渐近同值），LAD79 大概率**不是**在做极值性证明，而是在给这个
渐近合并点找一个更精确的数值/解析刻画。这与标题、摘要、参考文献列表三条独立信号都自洽。

### 2.6 Owen–Scalapino 1971（顺路核查，未新增进展）

`C.S. Owen and D.J. Scalapino, "S-state instabilities for retarded interactions", Physica 55,
691–695 (1971)`，DOI `10.1016/0031-8914(71)90320-X`——本次通过 Crossref 独立确认了 LAD79
自己也引用了这篇文献（是它 12 条参考文献之一），与仓库既有笔记 §10.2 第 4 点的定位一致。
**依然没有拿到全文**（同样的 Elsevier 1970 年代付费墙问题，且 *Physica* 在 1973 年拆分之前
的老刊没有免费扫描/无 arXiv）。这条"唯一未闭合的链接"本次**未能推进**，与既有笔记的判断
（"风险评估：低，但仍是一条未读的一手文献"）保持一致，不重复展开。

---

## 3. 风险等级判定

**判定：低到中风险，偏低，但不构成完全闭合（因为正文始终未读到）。**

支持"低风险"的证据链（全部间接，无一手正文）：
1. 摘要字面义（"range of validity of the asymptotic relation"）在数学结构上更像"数值收敛
   区间刻画"，不像"对任意谱测度 P 的严格不等式"——这是两种不同陈述，本仓库结论属于后者。
2. 标题聚焦的大 λ 渐近区正是 AD75 自己指出谱形状效应消失的区域，极值性问题在此区间趋于平凡，
   LAD79 不太可能在此证明极值性。
3. 现代最挑剔的同题材严格数学论文（KAY I/II/III）完全未引用 LAD79，而它们对 AD75 的无证
   断言做过审计式点名——如果 LAD79 已给出严格不等式或极值性证明，被 KAY 遗漏的概率较低。
4. Allen 本人的个人主页收录了不少 1970s "难以获取"的会议论文重印本，唯独没有把 LAD79
   收进去，暗示他本人也不认为这是一篇需要特别抢救留存的关键结果（弱证据，纯粹是选择偏好，
   不能排除遗漏）。

支持"不能完全排除中风险"的理由：
1. **全文始终未读到**——按仓库 CLAUDE.md 的纪律（"记录丢失时不得反向重建""承重结论必须能
   从零复核"），本笔记的判定**只能算 numerical-evidence 级别的间接推断**，不能当作
   analytic/certificate 级别的"已排除"结论。
2. LAD79 确实处理 μ*>0（间接证据），这与仓库正在进行的 μ*>0 深挖方向有实质重叠，具体重叠到
   什么程度取决于正文，目前无法判定。
3. Allen 本人是 LAD79 的共同作者，若他认为这篇文章确立了某种严格结果，理论上会在 1982/83
   的权威综述里体现，但**该综述正文本次同样未能读到**，仅确认了引用关系。

**不建议现在就大改仓库对外措辞**——现有证据不支持"LAD79 抢先"这个假设，风险判定为低到中、
偏低；但也不建议在没读到正文之前把这条风险彻底销号。

---

## 4. 对仓库现有措辞的具体修改建议

1. **暂不需要修改 `P4_g2_analytic.md` / `P4_einstein_extremal.md` 的"首次严格证明"核心措辞**——
   目前没有一手证据表明 LAD79 给出了同等或更强的结果。
2. **建议在两份定理文档或对应论文草稿的"相关工作"/历史脚注段补一句**（若尚未有），
   明确 LAD79 的存在并说明其未决状态，例如：
   > "Leavens, Allen & Dynes (Solid State Commun. 30, 595, 1979) present, per their abstract,
   > 'a rigorous derivation of the range of validity' of the large-λ asymptotic relation
   > k_BT_c ≈ 0.1827[λ⟨ω²⟩]^{1/2}; we were unable to obtain the full text (paywalled, no
   > author self-archive found). Based on the title, reference list, and the absence of any
   > citation in the most recent rigorous treatment of this constant (Kiessling–Altshuler–
   > Yuzbashyan 2024–2025), we believe this refers to a numerical convergence-range statement
   > for the fixed-coefficient asymptotic formula rather than a proof of the shape-independent
   > inequality T_c ≤ C_∞√(λ⟨ω²⟩) for general spectral measures P — but this has not been
   > verified against the primary source and should be revisited if the text becomes available."
3. **不建议**现在就把这条历史脚注升级为"经核实，LAD79 与本仓库结果不同"这种确定性措辞——
   证据链是间接的，只到"倾向性判断"，写成确定性陈述会越级。
4. **建议的后续动作（如果用户希望彻底闭合）**：
   - 通过馆际互借（ILL）或直接联系作者/机构申请扫描件——P.B. Allen（Stony Brook，
     `philip.allen@stonybrook.edu`，主页仍在维护）本人可能愿意提供私人重印本；
     C.R. Leavens 当时在加拿大国家研究委员会（NRC Canada）；R.C. Dynes 后在 UCSD。
   - 若拿到全文，优先核对：(a) 摘要里"range of validity"具体指的是 λ 的哪个下界数值，
     是否给出误差估计的严格上下界；(b) 是否有形如 Appendix 的严格证明段落，用了什么工具；
     (c) μ*>0 情形处理到什么程度（是否给出 C(μ*) 的显式界或只是数值拟合）。

---

## 5. 「读到的」vs「推断的」清单（诚实边界）

**一手-元数据（本次亲自查询数据库拿到）**：LAD79 标题/作者/期刊/页码/DOI（Crossref+OSTI+SD）；
LAD79 完整 12 条参考文献列表（Crossref）；Owen-Scalapino 1971 标题/DOI（Crossref）；
Leavens 1975 SSC17,1499 标题="A least upper bound on the superconducting transition
temperature"（Crossref）；Leavens 1976 SSC19,395 标题="On some recent results concerning the
superconducting transition temperature"及其 5 条参考文献（Crossref）；Semantic Scholar 确认
LAD79 `openAccessPdf.status="CLOSED"`、`abstract=null`；Allen 个人主页 `plist.html`+主页全文
（curl 拿到，逐行 grep 确认零命中 "leavens"/"1979"）；Allen-Mitrović 1982/83 综述的 Crossref
参考文献列表中存在指向 LAD79 的 DOI 匹配条目。

**一手（本次亲自读到正文/PDF）**：arXiv:2504.17435 全文 pdftotext，第 134 行 `[4,6]` 引用段；
KAY II/III 全文 pdftotext，grep 确认零 "Leavens" 命中，及 C_∞ 常数出现处的上下文。

**二手（转引自仓库既有笔记里已读过的一手引文，本次未重新取件）**：Carbotte RMP62,1027(1990)
p.1042 关于 Leavens 1975/1977 与 Owen-Scalapino(2.29) 的逐字段落——转引自
`theory/notes/literature/lit_einstein_extremality_novelty_2026-09-16.md` §10.2/§10.4。

**三手/AI 转述（WebSearch 工具自己合成，未经核对，不作为判决依据，仅供背景参照）**：
关于 Allen-Mitrović 正文如何具体描述大 λ 渐近系数、"Eliashberg 理论对 Tc 无严格上界"等
WebSearch 摘要段落——本笔记正文中所有此类内容均已明确标注为不可信来源，未计入 §3 判决。

**推断（不是任何来源的直接陈述，是我本人的推理）**：§2.1 关于"range of validity"字面义
与本仓库定理陈述结构不同的判读；§2.5 关于大 λ 渐近区极值性问题平凡化、故 LAD79 不太可能
证极值性的推理；§3 的整体风险等级判定。

---

## 6. 未闭合的缺口

1. LAD79 全文（2 页）始终未获取，是本笔记最大的缺口。
2. Allen & Mitrović 1982/83 综述正文（本应是判定 LAD79 意义的最佳二手来源）本次同样未获取，
   只确认了引用关系存在。
3. Owen–Scalapino 1971 全文仍未获取（与仓库既有笔记状态一致，未新增进展）。
