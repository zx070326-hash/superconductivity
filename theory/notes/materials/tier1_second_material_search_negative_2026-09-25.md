# Tier 1 第二材料检索：**干净否定** —— 仍是 YH₆ 单材料

**日期**：2026-09-25　**分支**：`ceiling-prescriptive-protocol-2026-09-24`（UNMERGED）
**上游**：`../ceiling/prescriptive_protocol_and_moussa_cohen_priority_2026-09-24.md` §0-C（Q2 干净样本塌到 YH₆ 单材料）

---

## 0. 这篇在回答什么

规定性协议的 **Tier 1** 主张是「用便宜谐波 DFPT 代替 SSCHA」。支持它的干净证据目前**只有 YH₆ 一个材料**
（n=3 @165 GPa：λ 跨 48.0% / ω̄₂ 跨 24.2%，而 `M2` 跨 **5.2%**、天花板_∞ 跨 **2.6%**），
且已明令**不得写成一般命题**。本轮派 `tier1-second-material`（researcher）去找**第二个**满足三条门的材料：

1. 同一材料、**同一压力**（同压门，见 §0-C 的自我订正 C）；
2. 至少两套声子方法（谐波 vs SSCHA/非谐）；
3. **两套都印 `⟨ω²⟩` 或 `ω̄₂`** —— 因为天花板吃的是 `M2 = λ⟨ω²⟩`，只有 `ω_log` **换算不出** `M2`。

## 1. 结论

**没找到第二个材料。Tier 1 的干净证据仍然是 YH₆ 单材料。**

这是**干净否定**，不是「检索失败」：卡住的位置每一处都定位到了具体的表和具体缺哪一列。

## 2. 逐条落点

| 候选 | 出处 | 卡在第几条门 | 逐字依据 |
|---|---|---|---|
| **PdH / PdD / PdT** | Errea et al. 2013，arXiv:1305.7123，**Table I** | **门 3**：只印 `λ, ω_log, Tc`，**无 `⟨ω²⟩`/`ω̄₂`** | 见 §3（主会话双仪器复核） |
| **H₃S / D₃S** | Mishra & Margine 2026（agent 报告） | **门 2 的独立性**：与 CSV 现有行**同构造**，不是独立第二套 | ⚠ 未经主会话复核 |
| **5 材料集** | arXiv:2608.14949 | **取不到**：Figshare DOI 自陈 `temporary, to be finalized upon publication` | ⚠ 未经主会话复核 |

## 3. 主会话从零复核：Errea 2013 Table I（**双仪器已齐**）

`tier1-second-material` 明确自陈**只用了页图、未跑 Bash** ⟹ 单仪器，按仓库规矩不得作为承重依据。
我用第二台仪器补上（我有 Bash，agent 没有）：

```
pdftotext -enc UTF-8 <library>/2013_Ion_Errea_1305.7123.pdf \
          theory/notes/_scratch_ir/_errea2013_recheck_utf8_2026-09-25.txt
```

文本层 23613 字符 / 670 行，正则计数：**`ω̄` 命中 0 次、`ω2` 命中 0 次、`omega.?2` 命中 0 次**；
`TABLE` 仅 1 处，其表头逐字：

> `TABLE I. Calculated λ, ωlog and Tc values for PH. ωlog values are given in cm−1 and Tc's in K.`

表体逐字（谐波三列 / SSCHA 三列 / 实验两列）：

```
        Harmonic          SSCHA            Exp.
        λ   ωlog  Tc      λ    ωlog  Tc    Tc[7] Tc[8]
PdH   1.55  205   47     0.40  405   5.0     9     8
PdD   1.55  150   34     0.46  304   6.5    11    10
PdT   1.55  125   30     0.48  257   6.9
```

⟹ **两台仪器一致**：Errea 2013 **确实不印 `⟨ω²⟩`**，PdH 族**进不了天花板检验**。
agent 的单仪器警示**至此解除**；判词由「未交叉」升为**已交叉、成立**。
**等级**：`numerical-evidence`（这是一条关于某篇论文印了什么的事实认定，非数值结论）。

### 3.1 顺带记两条（都**不**承重，仅备查）

- **λ 在谐波行三个同位素上恰好相同（1.55/1.55/1.55）** —— 与仓库已有结论
  「Errea 自陈 λ 与质量无关」（memory `dt-residual-is-spectral-shape-not-ill-conditioning-2026-09-20`）**同向**。
  ⚠ 但这里印的是**保留两位小数的表值**，「恰好相同」可能只是四舍五入的产物，
  **不得**当作「λ 严格质量无关」的证据（自由度审计：三个数可能只有一个自由度）。
- **不得用 `ω_log` 代替 `ω̄₂` 硬凑 `M2`。** 方向未审且必然错向：Jensen 只给 `ω_log ≤ ω̄₂`，
  用 `ω_log` 会**低估** `M2` ⟹ 低估天花板 ⟹ 产生**假 FAIL**（单边界方向规则：低估上界只支持 PASS）。

## 4. 未复核项（明确标出，不得当结论用）

以下两条**全部来自 agent 报告，主会话未取一手**，按「否定存在命题必须扫被断言的那个语料」的纪律，
在拿到全文前判词只能是 **未取得、未定**：

- Mishra & Margine 2026 与 CSV 现有 H₃S/D₃S 行「同构造」；
- arXiv:2608.14949 的 5 材料集、及其 Figshare DOI 的临时状态。

⚠ 本机 `arxiv.org` / `www.nature.com` 走 TLS 被劫持（证书 subject 为空，memory
`arxiv-and-nature-are-tls-intercepted-on-this-machine`）⟹ **不得 `-k` 绕过去取承重原文**。
这两条要闭合，需用户侧提供 PDF，或走 `harvest.aps.org` / `zenodo.org`（已实测干净）。

## 5. 对 Tier 1 主张的影响（**判词不变**）

Tier 1 的**声音性**不依赖材料数——它依赖的是上界方向，与紧度无关（`../../PROGRAM.md` §7 P4-E 条）。
本轮否定影响的只是**外部效度**：

- ~~「谐波 DFPT 足以定 `M2`」是跨材料一般命题~~ —— **仍未获得任何跨材料支持**；
- 当前可写的最强表述仍是 §0-C 已冻结的那句：**YH₆ 单材料、`numerical-evidence`、三套共享同一电子结构
  ⟹ 测的是声子算法、不是电子结构近似**。

⟹ **本轮不改任何判词**，只把「第二材料仍缺」这件事从口头状态变成可复核的落档。

---

## 引用产物

- `../_scratch_ir/_errea2013_recheck_utf8_2026-09-25.txt` —— 我这次跑的 UTF-8 文本层（第二仪器），已 `git add -f`。
- 源 PDF 在 `library/pdf/`（受 `*.pdf` 版权 gitignore，不入库）；复现命令见 §3。
