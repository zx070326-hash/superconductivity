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

---

## 7. 主会话追加（2026-09-20）：三条新证据 + 一次判决重构

> 本节所有数字由主会话本人取件，不转引上文。分「内容命题」与「影响力命题」两段写
> （`CLAUDE.md` 订正纪律：两者是独立命题，不许互相冒充）。

### 7.1 【内容命题】"range of validity" 是这门文献的**行话**，有一手样板 [一手]

上文 §2.1 把「range of validity ≠ 对所有 P 的不等式」列为**我的读法**（§5 明列为「推断」）。
**现在它不再是读法**：仓库自己持有的 Carbotte RMP 62, 1027 (1990) 全文抽取
（`../_scratch_ir/CARB.txt`，12487 行）里，p.1046 右栏有一段**同体裁样板**，
讨论的恰是一条 `T_c = const·√(λ⟨ω²⟩)` 型渐近式的有效范围（`CARB.txt:2026–2046` 逐字）：

> `T_c = 0.148 θ`　(2.53)
> …"In the range 0.6 < ω_E/θ < 1.5, the curves for T_c/θ are reasonably flat… For μ*=0.1,
> we should take the constant somewhat below the maximum value of 0.175, and 0.148 is reasonable.
> **The range defined above, over which Eq. (2.53) should be valid, is roughly 1.3 < λ < 3.3.
> The more rigorous approach on which (2.53) is based gave the slightly smaller range quoted above.**"

其中「the more rigorous approach」= **Leavens–Carbotte (1974)**，「quoted above」= `CARB.txt:2019`
的 **1.2 < λ < 2.4**（"For λ in the range 1.2 < λ < 2.4, Leavens and Carbotte (1974) noticed that
the function on the right-hand side of (2.52) was approximately constant"）。

⟹ **在这门文献自己的用语里，「渐近关系式的 range of validity」= 一个 λ 区间，在其中固定系数的
近似式数值上够准**，且「more rigorous」在同一句里就是用来形容这种 λ 区间估计的。
这是一个**限定 λ 区间上的双边近似陈述**；本仓库的 `T_c ≤ C_∞√(λ⟨ω²⟩)` 是**对所有 P、所有 λ 的
单边不等式**。两者数学结构**不相交**（互不蕴含）。LAD79 摘要用的正是这个词组。
**级别：这条把 §2.1 的判读从「推断」升为「有一手同体裁样板支持的推断」，仍不是「已读正文」。**

### 7.2 【内容命题】篇幅硬约束 [一手-元数据，Crossref，主会话本人拉取]

`api.crossref.org/works/10.1016/0038-1098(79)90103-0` 返回：
`page = 595-596`（**两页**）、`reference-count = 12`、`volume 30 / issue 10 / 1979-06`。

⟹ 两页、12 条参考的 SSC 短讯，**不可能同时**容纳 (i) 无穷矩阵顶本征值的严格控制
与 (ii) 对所有谱测度 P 的极值性证明——而这两件事**都**是本仓库 `C_∞` 结论的必需件。
（AD75 自己的 0.182 是 N+1=64 阶截断**数值**对角化得来的。）
配合 §2.5 的论证（LAD79 标题限定的大 λ 区，正是 AD75 p.908 自陈「shape doesn't affect T_c」
的区域 ⟹ 极值性在该区**平凡**），(ii) 的缺席有两条独立理由。

### 7.3 【影响力命题——与上两条**不同**的命题】Carbotte RMP 1990 整本书目零条目 [一手]

上文 §2.3 只说「Carbotte 正文两处点名相关文献都不含 LAD79」。**可以升一档**：
`CARB.txt:12000–12030` 是该综述书目的 Leavens 区块，逐条为
1974 SSC15,1329 / 1975 SSC17,1499 / 1977 J.Phys.F7,1911 / &Carbotte 1971,1972,1974,1977 /
&Fenton 1980 / &Ritchie 1984,1985。按 RMP 字母序，「Leavens, C. R., P. B. Allen, and R. C. Dynes, 1979」
应落在 `1977 J.Phys.F` 与 `&Carbotte 1971` 之间（第二作者 Allen < Carbotte）——**该位置为空**。
⟹ **这本 ~1000 条参考的领域权威综述完全没有引用 LAD79**，而作者 Carbotte 是 Leavens 的长期合作者，
且该综述**专门讨论了同一条渐近式的有效范围**（§7.1）却引了 Leavens–Carbotte 1974。

Crossref 同时给出 **`is-referenced-by-count = 4`**（47 年被引 4 次）。

> ⚠ **方向审计（必须写，否则是 W1 形状的越级）**：§7.3 界的是**影响力 / referee 会不会知道它**，
> **不界内容**。一篇被遗忘的论文照样可以持有优先权。承担「内容」判决的只有 §7.1/§7.2/§2.5，
> 承担「风险暴露」判决的才是 §7.3。**两者不许互相冒充。**

### 7.4 取件路线的重新评估

- **Allen & Mitrović SSP 37**（§6 缺口 2）：本轮确认 **无免费全文**（ScienceDirect
  `S0081194708606657` 付费墙；Semantic Scholar 仅元数据）。且由 §7.3 可下调其预期收益——
  同代、同题材、同圈子的 Carbotte 综述连引都没引。
- **NRC Publications Archive**（Leavens 雇主机构库）：本轮仅经搜索引擎间接试探，**未直接站内检索**，
  仍是一条**未走过**的路。
- **★ 探查对象应当换人**：LAD79 的前身 **Leavens, SSC 17, 1499 (1975)，标题
  "A least upper bound on the superconducting transition temperature"** ——
  **这个标题离本仓库的主张比 LAD79 近得多**，且仓库在**同一个约束**（固定面积 A）下自己有
  `k_BT_c ≤ 0.2309 A`、等号 iff Einstein 谱（`../coulomb/hmu_tolmachev_lens_2026-09-18.md:400`）。
  仓库现有评估（`lit_einstein_extremality_novelty_2026-09-16.md:143`）把它列为「变分/线性响应论证，
  **非严格** [L2]」——**即二手评估，正文未读**。⟹ **若真有 1970s 优先权风险，它住在这篇，不住在 LAD79。**
- 待核线索（未验证）：搜索中出现的 **Leavens, C. R., SSC 31, 837 (1979)** 同年独著，未经 Crossref 确认。

### 7.5 顺带：两篇仓库视野外的近期文献 [一手-元数据/摘要]

- `arXiv:2505.02894` = **Gnezdilov & Boyack, PRB 112, L180504 (2025)**,
  "Upper bound on T_c in a strongly coupled electron-boson superconductor"。**仓库全库零命中。**
  内容是二维 Yukawa-SYK 模型**超出 Migdal** 的 `λ_E = λω_D/ε_F` 大值区，`T_c → 0.04 ε_F` 饱和；
  小 `λ_E` 时复现 `T_c ≈ 0.18 ω_D√λ`。**不是本仓库定理的竞争者**（不同约束、不同理论边界、无 μ*），
  但属于「T_c 上界」相关工作，应进文献表。
- `arXiv:2506.19326`（Sadovskii）仓库已有记录（`../materials/f_of_mustar_2026-09-17.md`：μ*=0、
  "elementary estimate"），无新增风险。

---

## 8. 一手全文到手（2026-09-20）：缺口 1 关闭，判决重写

> **取件路线**：用户直接提供 PDF（微信文件），已复制到 `../_scratch_ir/LAD1979.pdf`。
> 2 页，TextBridge OCR 扫描件，**文本层糊得不可用**
> （实测：`"ArerliagtoiroonuskBdTecr~iOva.tIi8o2n7o(fA<twh2e>]ria~sngepreosfenvtaeldi.dity"`）。
> **本节全部逐字内容取自页图** `../_scratch_ir/_lad1979_p-1.png` / `_lad1979_p-2.png`，
> 由主会话本人阅读，**不采信文本层任何字符**（`CLAUDE.md`：文本层与页图互为对照）。
>
> 按订正纪律，「旧陈述错在哪」与「新陈述凭什么对」分两段写（§8.2 / §8.3）。

### 8.1 全文骨架（一手，逐字）

**标题**：*The superconducting transition temperature for very large electron-phonon coupling*
（C.R. Leavens, NRC Canada；P.B. Allen, SUNY Stony Brook；R.C. Dynes, Bell Labs；
Received 8 January 1979 by M.F. Collins）。

**摘要（页图逐字，与 §2.1 的元数据摘要一致）**：

> "A rigorous derivation of the range of validity of the asymptotic relation
> k_BT_c ≈ 0.1827[λ⟨ω²⟩]^½ is presented."

**全文的问题设定**（p.595 左栏）：AD75 对「λ 充分大」的判据
*"was too loosely investigated in ref. 2, and is **wrong in extreme instances**"*
（ref. 3 = Leavens SSC 19, 395 (1976)）。
本文把 AD 的 `λ ≳ 10` 换成 **`λ⟨ω²⟩ ≫ ω²_max`**，并证明后者充分。
右栏逐字：*"**The condition λ≫1 is thus necessary but not sufficient.**"*

**逐式（页图）**：

| 式 | 内容 | 地位 |
|---|---|---|
| (1) | `Limit_{λ⟨ω²⟩≫ω²_max}  k_BT_c ≈ 0.1827[λ⟨ω²⟩]^½` | **渐近极限**，注意原文是 `≈` |
| (2)(3) | `λ⟨ω²⟩ ≡ 2∫dω α²F(ω)ω`；`λ ≡ 2∫dω α²F(ω)/ω` | 定义，与仓库 (1.3)(1.4) 逐字同 |
| (4) | `λ(n) ≡ 2∫₀^{ω_max} dω α²F(ω)ω/(ω²+(2πnk_BT_c)²)` | Matsubara 核 = AD Eq.(10) |
| (5) | `λ(n) ≈ λ⟨ω²⟩/(2πnk_BT_c)² + …` | 大 T_c 展开 |
| (6) | `2πk_BT_c ≫ ω_max` | 判据（"also a sufficient condition"） |
| (7) | `λ⟨ω²⟩/(ω²_max+(2πk_BT_c)²) ≤ 2∫₀^{ω_max}dω α²F(ω)ω/(ω²+(2πk_BT_c)²) ≤ 1` | 左端 "obvious"，**右端 "is proved in ref. 2"（= AD75）** |
| (8) | `λ⟨ω²⟩/ω²_max ≤ 1 + (2πk_BT_c/ω_max)²` | (7) 两端相连 |
| (9) | `λ⟨ω²⟩ ≫ ω²_max` 充分保证 (6) | 本文主结果 |
| (10)(11) | `Limit_{x≫ω_max} k_BT_c = 0.2584 x`，`x ≡ [∫₀^{ω_max}dω α²F(ω)ω]^½` | η/M 形式 |

**常数的地位**（p.596 左栏，逐字）：

> "…must take the form `T_c/√(λ⟨ω²⟩) = c(μ*)` where the function c depends only on μ*.
> **Numerical solution of the gap equations gives c(0) = 0.1827.**"

**μ\***（p.595 左栏，逐字）：

> "**For simplicity, we consider explicitly only the case in which the Coulomb
> pseudopotential, μ\*, is zero.**"

**η/M**（p.596 右栏，逐字）：

> "(For elements `x² = η/2M` where η is the McMillan-Hopfield parameter^{1,9} and M is the
> ionic mass.) This suggests that in the limit of very strong electron-phonon coupling
> **it is the first moment of α²F, rather than λ, which is the single most important
> parameter** characterizing α²F.^{2,10,11,12}"

**Table 1**（幂律谱 `α²F(ω) = (nλ/2)(ω/ω_max)^n θ(ω_max−ω)`；
`λ_(5%)` = 使 Eq.(1) 精确到 5% 所需的 λ 下界）：

| n | λ_(5%) | (⟨ω²⟩λ/ω²_max)_(5%) | (x/ω_max)_(5%) |
|---|---|---|---|
| 0 | ∞ | 3.4 | 1.3 |
| ½ | 18.8 | 3.8 | 1.4 |
| 1 | 12.3 | 4.1 | 1.4 |
| 2 | 9.2 | 4.6 | 1.5 |
| ∞ | 7.0 | 7.0 | 1.9 |

表注："To achieve accuracy to within 2%, λ must be larger than λ_(5%) by a factor of ~2.5."
收尾："**λ is not a very useful parameter while the first moment of α²F is.**"

### 8.2 §2.4 错在哪（撤回）

§2.4 写「LAD79 是否处理 μ\*>0：**间接证据倾向于"是"**」，依据是 §2.3 里 arXiv:2504.17435
引 LAD79 谈 `C(μ*)`。

**这条推断为假，现撤回。** 原文那句
*"For simplicity, we consider explicitly only the case in which the Coulomb pseudopotential,
μ\*, is zero"* 逐字否掉它。LAD 只是**写下** `T_c/√(λ⟨ω²⟩) = c(μ*)` 这个函数形式
（即「渐近律的谱形无关性对任意 μ\* 成立」），然后**只算 c(0)**。
arXiv:2504.17435 引它作 `C(μ*)` 的出处是**宽引**——LAD79 不提供 `C(μ*)`。

连带：§3「支持不能完全排除中风险」第 2 条（"LAD79 确实处理 μ\*>0…与仓库 μ\*>0 深挖方向
有实质重叠"）**前提消失，该条作废**。仓库的 μ\*>0 工作与 LAD79 **零重叠**。

### 8.3 §7 的两条推断被一手确认（不是撤回，是升级）

§7.1（"range of validity" 是行话，指**定系数渐近式在数值上准确的 λ 区间**，不是对所有 P 的
不等式）与 §7.2（2 页 12 引的篇幅装不下极值性证明）——**两条都由正文逐字确认**：

- Table 1 的说明栏字面就是 *"λ_(5%) is the numerical value which λ must exceed for eq. (1)
  to be accurate to within 5%"* ⟹ **"range of validity" = λ 的数值收敛区间**，逐字命中 §7.1。
- 文中**没有任何**谱形极值性论证；全文的 "rigorous derivation" 指的是 (7)→(8)→(9) 这条
  **充分条件**链，不是常数 0.1827 的推导——后者明说是 "Numerical solution of the gap equations"。
  §2.5（"LAD79 大概率不是在做极值性证明"）同样确认。

### 8.4 三分法：什么是 LAD 的、什么是我们的

这是本节唯一承重的输出。**三件事必须分开，不许互相冒充：**

**(a) 常数 0.1827 —— LAD 有，且比 AD75 多一位。**
AD75 Eq.(25) 印 `0.182`；LAD79 Eq.(1) 印 `0.1827`，与我们认证的
`C_∞ = (1/2π)√𝔤(2) = 0.1827262477…` **逐位相符**（相对差 **+0.0144%**，落在 LAD 的末位内）。
⟹ **仓库不得声称「首次得到这个常数的数值」。** 我们的增量是它的**认证包围**
（`P4_g2_analytic.md` CERT-D，12+ 位，有理算术），而 LAD 的来源是**数值解 gap 方程**。
这是严格性升格，不是数值发现。

**(b) 不等号方向 —— LAD 的 Eq.(8) 与我们的天花板方向相反。**
把 (7) 两端相连解出 T_c（纯代数，方向逐步审；脚本
`../_scratch_ir/_chk_main_lad1979_floor_2026-09-20.py` F1）：

```
λ⟨ω²⟩ ≤ ω²_max + (2πk_BT_c)²   ⟹   k_BT_c ≥ (1/2π)√(λ⟨ω²⟩ − ω²_max)     [右端为负时空洞]
```

**这是 T_c 的下界。** LAD79 **从未**声称 `0.1827√(λ⟨ω²⟩)` 是对一切 λ、一切谱的**上界**——
他们证的是极限存在、并刻画极限成立的区间。
⟹ **全局天花板 `T_c < C_∞√(λ⟨ω²⟩)`（即 `𝔣 ≤ 1`，对一切 λ>0、一切 P∈𝒫）仍然是我们的。**
§3/§4 的「首次严格证明」措辞**不需要改**——但**理由要换**成这一条（方向），
不再是原来那条「我们没读到正文所以不能排除」。

**(c) η/M 一行公式 —— 是 LAD 的，仓库必须补引。**
`x² = η/2M` 加 Eq.(10) `k_BT_c = 0.2584x` 合起来就是
`k_BT_c = 0.2584√(η/2M) = 0.1827√(η/M)`，**即 `P4_eta_over_M.md` 定理 A 的数值形式**。
而且 Eq.(10) 与 Eq.(1) **不是两条独立信息**：`λ⟨ω²⟩ = 2x²` ⟹
`√2 × 0.1827 = 0.258377` vs LAD 印的 `0.2584`（差 0.0086%，即末位）。
⟹ **`P4_eta_over_M.md` §8 必须加 LAD79 条目**，且 §6「本文新结果」里「定理 A」一项的措辞要改：
新的是**「对一切 λ 成立的上界」**这一步，**不是**「把界写成 η/M 的形式」——
后者 LAD 1979 已经做了，连 McMillan–Hopfield 的归属（refs 1,9 = McMillan 1968 / Hopfield 1969）
都写在同一句括号里。"first moment rather than λ is the single most important parameter"
这句**物理论断**同样是 LAD 的。

### 8.5 另一条探查方案：Table 1 是可跑的第三方检验，且跑过了

用户问「有没有别的探查方案」——**有，Table 1 本身就是**。它给了 5 个数值，
每个都是「在这条谱上，λ 要多大 Eq.(1) 才准到 5%」，
**这完全可以用我们的求解器重算**，从而对 `(求解器, C_∞)` 做一次 1979 年的独立检验。

脚本：`../_scratch_ir/_chk_main_lad1979_table1_2026-09-20.py`（A–D 段）
+ `../_scratch_ir/_chk_main_lad1979_invert_2026-09-20.py`（E 段）。
求解器 `theory/solvers/eliashberg_iso.py`，μ\*=0，N=512（截断误差 1.3e−6，C0 段实测），
谱按 LAD 页图的**原式** `α²F = (nλ/2)(ω/ω_max)^n θ(ω_max−ω)` 构造（该式给 `λ_check = λ`，已验）。

> **仪器自检先于结论**（C1 段）：`n=½` 的 `α²F/ω ~ ω^{−1/2}` 在 0 附近有可积奇点，
> 均匀 Simpson 网格只有 ~1% 精度，第一版因此给出 `λ(𝔣=.95) = 21.54`，**看起来像 LAD 偏了 14.6%**。
> 换成聚点网格 `ω = ω_max·t³`（ngrid=20001）后 λ 与 ⟨ω²⟩ 的解析值都对到 **1e−6**，
> 该行变成 18.850。**这是我的仪器，不是 LAD 的数**——
> 正是 `CLAUDE.md` 点名的 W1 热区（订正别人时新陈述被免检放行），落笔前抓住。
> 另：`frak_f` 一律用**求解器实际那条谱**的 `λ_total` 与 `ω₂` 归一，不用标称 λ，
> 否则网格误差会藏进比值里、自己抵消掉。

**先数自由度**（`CLAUDE.md`：把一组数字当 N 个独立检验之前先数自由度）。
幂律谱下 `⟨ω²⟩/ω²_max = n/(n+2)` 是**纯形状常数**，于是

```
col3 = col2 × n/(n+2) :  3.760 / 4.100 / 4.600 / 7.000   vs LAD 印的 3.8 / 4.1 / 4.6 / 7.0
col4 = √(col3/2)      :  1.378 / 1.432 / 1.517 / 1.871   vs LAD 印的 1.4 / 1.4 / 1.5 / 1.9
```

⟹ **Table 1 的 15 个数只有 5 个自由度（每行 1 个）**；第 3/4 列是第 2 列的函数，不是独立检验。

**检验结果**（4 个有限谱行 + 1 条表注）：

| n | LAD λ_(5%) | 我们的 λ(𝔣=0.95) | 比值 | 我们在 LAD 的 λ 处的 𝔣 |
|---|---|---|---|---|
| ½ | 18.8 | 18.850 | 1.003 | 0.9499 |
| 1 | 12.3 | 12.291 | 0.999 | 0.9500 |
| 2 | 9.2 | 9.171 | 0.997 | 0.9502 |
| ∞ | 7.0 | 6.973 | 0.996 | 0.9502 |
| 表注 "~2.5" | — | λ(𝔣=.98)/λ(𝔣=.95) = 2.52 / 2.51 / 2.50 / 2.48 | — | — |

**「5%」的读法由数据定死，不是我挑的**：若分母取 `T_c` 本身（`𝔣 = 1/1.05`），
四行全部偏高 4.5–5.3%；取渐近式作分母（`𝔣 = 0.95`）四行全部落在 0.4% 内（E3 段）。

**反解 + 先算灵敏度**（`CLAUDE.md`：判据若对该量不敏感，吻合就是空的）：

```
dln𝔣/dlnλ = 0.0387 (n=1) / 0.0216 (n=∞)   ⟹  c 错 1% 会把 λ_(5%) 挪 29% / 59%
LAD 印 3 位有效数字（±0.4% in λ）          ⟹  本检验对 c 的分辨力 = ±0.0159% = ±2.9e−5
```

分辨力算清楚之后再看吻合：

```
c_implied = C_∞ · 𝔣_ours(λ_LAD) / 0.95  =  0.1827011 / 0.1827335 / 0.1827570 / 0.1827645
          = 0.1827390 ± 0.0000285  (4 行，1σ)      vs   C_∞ = 0.1827262477  (+0.0070%)
```

**4 行散布 2.85e−5 恰好等于分辨力 2.9e−5** ⟹ 残余散布**全部**是 LAD 的印数舍入，
没有可归因于物理的偏差。

**这个检验测什么、不测什么（逐箭头）**：
- ✔ 测到：`(我们的 Eliashberg 求解器) × (C_∞ 的值)` 这个**组合**，与 1979 年一次独立的
  gap 方程数值解在 `±0.016%` 内一致。
- ✘ **分不开这两者**：若求解器偏高 x% 而 C_∞ 偏低 x%，本检验全盲。
  要分开必须另找只依赖其中一个的观测量。
- ✘ **不测天花板**。D 段实测四条谱上 `𝔣 < 1` 且单调升（max 0.998857），
  这是天花板的**推论没被证伪**，**不是**天花板的证据（LAD 只证极限，没证上界）。
- 级别：**numerical-evidence**（float64 求解器 + 有限 N + LAD 的 3 位印数）。

### 8.6 附带收获：LAD Eq.(7)(8) 给出 𝔣 的**地板**

由 §8.4(b) 的代数，两边除以 `C_∞√(λ⟨ω²⟩)`：

> **`𝔣 ≥ (1/√𝔤(2))·√(1 − ω²_max/(λ⟨ω²⟩))`**，　`1/√𝔤(2) = 0.8710020867`
> （CERT-D 包围两端同值到 10 位），　μ\*=0，`λ⟨ω²⟩ > ω²_max` 时非空。

这是**第三方来源的 `𝔣` 下界**，正好补在我们天花板 `𝔣 ≤ 1` 的另一侧
——对「𝔣 做成头条」那条线：天花板是我们的，地板可以引 LAD79 + AD75。

**直测**（`_chk_main_lad1979_floor_2026-09-20.py`，四条谱 × 六个 λ）：
19 个非空行**零违反**，最小余量 **+0.130**。
**前提单独直测**（唯一非平凡输入 = (7) 右端 `λ(ν₁) ≤ 1`，由 AD75 证）：
16 个点上 `max λ(ν₁) = 0.7581 ≤ 1`，成立。

**三条诚实标注**：
1. **严格性继承 AD75**：(7) 右端的证明在 AD75，我们**未独立复核**那一步。
   ⟹ 这条地板目前只能标 **numerical-evidence**（我们这侧）+ 引 AD75（他们那侧）。
2. **渐近松 13%**：`λ→∞` 时地板 → 0.8710 而真值 → 1。
3. **`λ(ν₁) → 0.758` 不是 𝔤(2) 的独立确认**——在渐近区 `λ(ν₁) = 1/(𝔤(2)𝔣²)`，
   故 `𝔣→1 ⟺ λ(ν₁)→1/𝔤(2) = 0.758644` 是**恒等式**。
   （λ=100 Einstein 实测 0.758094，与恒等式预测 0.763884 差 0.76%，恰好等于 (5) 里被略去的
   `ω²/(2πk_BT_c)²` 项；对得上，但这只是自洽，**不是证据**。写下来是因为它差一点被我当成
   「一个完全不同的观测量独立确认了 𝔤(2)」——它不是。）

### 8.7 §6 缺口状态更新

1. ~~LAD79 全文未获取~~ ——**关闭**（用户提供 PDF，主会话已读两页页图）。
2. Allen & Mitrović 1982/83 综述正文 —— **仍未获取**。
   但其判定价值已大幅下降：LAD79 正文已直接回答本笔记要问的全部问题。
3. Owen–Scalapino 1971 全文 —— **仍未获取**（LAD79 的 ref. 6，Matsubara 版 Eliashberg 的出处）。

**新缺口**：LAD79 ref. 3 = Leavens SSC **19**, 395 (1976)，正文说它指出 AD 判据
"wrong in extreme instances"。那篇是 Table 1 里 `n=0 ⟹ λ_(5%)=∞` 这个反例的来源，
与我们「高度弥散谱」的边界讨论同题材，**值得取件**（优先级：低——结论已由 LAD79 正文转述）。
