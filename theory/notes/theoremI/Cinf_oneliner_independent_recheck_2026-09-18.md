# C_∞ 与「一行公式」链条的主会话独立复算（2026-09-18）

> 起因：会话意外中断后接续。`ceiling/oneliner_zero_parameter_2026-09-18.md`（sonnet subagent 产出）
> 把 **`T_c ≤ C_∞√(η/M)`** 当作方向 1 的交付件，而它的规范来源 `theory/derivations/P4_eta_over_M.md`
> **不在本分支的工作树里**（在未合并分支 `etaM-corrected`，HEAD `ed6de52`）。
> 按 CLAUDE.md §1「承重结论从零独立重算」，本篇对该链条的**数学部分**逐步重算，不转引任何一方。
>
> 复现：`python theory/notes/_scratch_ir/_chk_Cinf_from_g2_rational.py`（**8/8 通过**，含投毒负对照）。
> 快照：`theory/notes/_scratch_ir/_etaM_snapshot.md` = `git show ed6de52:theory/derivations/P4_eta_over_M.md`。

---

## 1. 判决表

| 环节 | 我的判决 | 级别 |
|---|---|---|
| (1.2) `λF_T(k) = λ(n)`（本文核 ⟺ AD Eq.(10) 的松原核） | **成立**，逐 k 恒等，代入 `P=2α²F/(λω)` 直接得 | `analytic` |
| (1.3) `⟨ω²⟩_𝒫 = (2/λ)∫ωα²F = ω̄₂²`（AD Eq.(5), n=2） | **成立** | `analytic` |
| (1.4) `λ⟨ω²⟩ = 2∫ωα²F dω` | **成立**（离散 Einstein 混合上有理算术验，偏差恒为 0） | `analytic` |
| (1.5) `λ⟨ω²⟩ = η/M`（McMillan Eq.(22) 代入） | 代数**成立**；一手性依赖对方所述的 McMillan 原文核对，**我未取原文**（见 §3） | `analytic`（代数）＋ **未复核**（文献） |
| 定理 A ⟸ 定理 1(c) | **成立**：`P4_main_theorems.md:56` 的 (c) 逐字给 `T_c < (1/2π)√(𝔤(2)⟨ω²⟩λ)`，方向对（`𝔨` 严格递减、判据 `λ𝔨=1`） | `analytic` |
| `C_∞ = √𝔤(2)/2π` 的数值 | **成立**，我的外向包围与对方逐位相符（见 §2） | `rational-exact` |

**结论：这一行公式的数学链条我从零走通，没有发现缺陷。**
唯一我没有亲自落地的是 McMillan Eq.(19)–(24) 的一手转录（§3）。

---

## 2. `C_∞` 的独立外向包围（纯有理算术，零浮点）

输入只有两样：`𝔤(2)` 的 CERT-D 有理端点（`P4_g2_analytic.md:13`）与 `π` 的保守有理包围。
`√` 用 `math.isqrt` 做整数开方取 40 位标度的上下端，除法按**方向**分配（下端 = 小分子/大分母）：

```
𝔤(2) ∈ [1.318140211762671, 1.318140212373478]        (CERT-D)
⟹ C_∞ = √𝔤(2)/2π ∈ [0.182726247746701, 0.182726247789037]      宽 4.234e−11
```

- 与 `P4_eta_over_M.md` §2.A 写的 `[0.182726247746, 0.182726247790]` **逐位相符**，且对方那一对端点是我这个区间的**正确外舍入**（未凭空变锐）。
- 与 `P4_main_theorems.md:62` 的注（同一对端点 + 真值 `[0.18272624774670140, 0.18272624778901]`）一致。
- 报告常数 **`0.1827262478` 是合法上界**（≥ 包围上端 `0.1827262477890`）——方向对：上界常数必须取包围的**上端**。
- KAY 的 `f(0)=0.1827262477…` 落在包围内。
- **投毒负对照**：把 `𝔤(2)` 换成 `[1.33,1.34]`，上端升到 `0.184235…`，`0.1827262478` 立刻失去上界资格 ⟹ 该断言确实被数据驱动，不是恒真。

### ★ 一条引用纪律（新，本次复算发现）

**CERT-B 撑不住这 10 位。** 用 CERT-B（`P4_main_theorems.md:118` 指定的**论文正文口径**，12 个整数 / N₁=40）：

```
C_∞ ∈ [0.182724107209718, 0.182727728710722]      ⟹ 只稳得住 0.18272（5 位）
```

`0.1827262478` **不是** CERT-B 下的合法上界。两份文档本身是自洽的（`:118` 已写明「CERT-D 列在附录 A 给需要 10 位的读者」），
但 `:62` 的正文句子直接写 `T̃_c = 0.1827262477…·√(λ⟨ω²⟩)`。
⟹ **凡引用 10 位常数处必须点名 CERT-D**；正文若只陈述 CERT-B，那一句得带附录 A 的交叉引用。
（CERT-C 同样不够：`[0.182726196, 0.182726282]`，只稳 7 位。）

---

## 3. 我**没有**做的（诚实边界）

1. **McMillan 1968 Eq.(19)–(24) 的一手核对**。`P4_eta_over_M.md` §1.C 自陈「已取到原文 PDF 并逐式核对」，
   我**没有**重做这一步，因此 (1.5) 的**文献一手性**在本篇里仍是转引。
   代数重排本身我验了：`λ⟨ω²⟩ = 2∫ωα²F = 2·[N(0)ħ⟨𝒢²⟩/(2M)] = η/M`（ħ=1）。
2. **`P4_eta_over_M.md` §4 的材料数值表**（LaH10 / H3S 的 `η/M` 与界/SSCHA 比值）未复算。
   注意它与 `materials/` 的两条已知坑相邻：`hydride_csv_audit`（`Tc_AD` 是循环量）与
   `omega2-not-missing`（`ω̄₂` 就是 Eq.1c 的 `⟨ω²⟩`）。**要引用那张表必须先重算。**
3. **定理 B / §2.C′ 的认证包围**：对方自陈严格性以已合并的 `P4_universal_tc_envelope.md` 为准，本篇未验。
4. **(A5) 多组分分解**：对方自标「未独立于 McMillan 框架重证、`η_a` 二手转引」，我同意这个标注，未验。

---

## 4. 分支/交付状态（这是本次接续要先说的事）

- `P4_eta_over_M.md` 在 **`etaM-corrected`（`ed6de52`）**，`p4-conj2` 工作树里没有它。
  `theory/derivations/INDEX.md:114` 把它列为「方向 C · 未合并（+2）」——一致。
- 后果：`ceiling/oneliner_zero_parameter_2026-09-18.md` 现在**指向一个本分支不存在的文件**，
  与 CLAUDE.md §2 的「脚本毕业规则」同形（引用的证据必须在同一个 clone 里可见）。
  **不自行合并**（CLAUDE.md §3：合并要用户放行）；本篇 + 快照文件先把可复核性补上。

---

相关笔记：[`../ceiling/oneliner_zero_parameter_2026-09-18.md`](../ceiling/oneliner_zero_parameter_2026-09-18.md)、
[`../isotope/isotope_pdh_discriminating_test_2026-09-18.md`](../isotope/isotope_pdh_discriminating_test_2026-09-18.md)（界的 α=½ vs 真实 α）、
[`../materials/hydride_csv_independent_audit_2026-09-17.md`](../materials/hydride_csv_independent_audit_2026-09-17.md)、
[`../materials/omega2_data_sources_2026-09-17.md`](../materials/omega2_data_sources_2026-09-17.md)。
相关 memory：`[[p4-weak-coupling-prefactor-0p2530]]`（同形：继承来的数字必须先复核）、
`[[isotope-exponent-is-the-sharper-second-observable-2026-09-18]]`。
