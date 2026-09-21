# theory/derivations/ 目录约定

> 本文件是该目录的单一真源；`INDEX.md` 头部只保留一行指针。CLAUDE.md §2 指向此处。

## 严格性等级（房规，不可越级上调）

- **解析 / analytic** = 已证明，零计算机辅助（或仅有理数精确算术）。最高等级。
- **认证 / certificate** = 外向舍入区间算术（mpmath.iv / Fraction）给出的 sound 包围。
- **数值证据 / numerical-evidence** = float64 点求值 / 有限网格，仅为证据、非证明。
- 一篇稿常是**混合**（不同命题不同等级）；`INDEX.md` 标注的是「承重结论」的等级。

## 维护规则

改动某稿的严格性/门/结论时，同步改 `INDEX.md` 对应行。
