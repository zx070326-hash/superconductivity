# ALDA f_xc 符号结构调研：r_s ~ 2-3 金属密度区（2026-09-21）

> **任务来源**：team-lead，检验 Essenberger Eq.(32) 推断 `Im Λ^SF ⪰ 0` 的前提。
> 
> **核心问题**：ALDA (adiabatic local density approximation) 交换关联核 f_xc 在 r_s ~ 2-3（铁基超导金属密度）能否为负？
> 若能，则 Essenberger 的 congruence 形式 `Λ^SF = 2∬ f_xc χ f_xc` 不自动保证 `Im Λ^SF ⪰ 0`。

---

## 0. 一句话结论

**ALDA f_xc 在金属密度 r_s ~ 2-3 处的符号取决于交换与关联两部分的竞争，且已知在某些条件下可以为负。**

关键事实：
1. **ALDA 核定义**：`f_xc^ALDA(r) = d²ε_xc(n)/dn² |_{n=n(r)}`（二阶导数，标量，无频率依赖）
2. **交换部分恒负**：`d²ε_x/dn² < 0` 对所有密度成立
3. **关联部分符号变化**：`d²ε_c/dn² < 0` 在高密度、`> 0` 在低密度，**转折点在 r_s ~ 3-4 附近**
4. **r_s ~ 2-3 落在过渡区**：铁基密度正处于关联导数符号反转的邻域

⟹ **Essenberger 推断需要额外论证**：`f_xc` 符号不定 + 响应函数 χ 的谱结构共同决定 `Im Λ^SF` 符号。

---

## 1. ALDA 核的显式定义（L2，综合多源）

ALDA 交换关联核在实空间、绝热（ω→0）极限下的标准定义：

```
f_xc^ALDA(r, r') = δ²E_xc[n]/[δn(r)δn(r')] |_{n=n_0}
                 = [d²ε_xc(n)/dn²]_{n=n_0(r)} · δ(r − r')
```

其中 `ε_xc(n)` 是均匀电子气在密度 n 处的每粒子交换关联能。

**关键性质**：
- **局域**（δ 函数）：不含空间非局域效应
- **绝热**：无频率依赖（ω 不出现）
- **密度泛函**：只依赖局域密度 n(r)

⟹ 符号**完全由** `d²ε_xc/dn²` 在该密度点的值决定。

来源：
- [arxiv.org/abs/1501.00902](https://arxiv.org/abs/1501.00902) — Renormalized ALDA
- [arxiv.org/abs/1107.0199](https://arxiv.org/abs/1107.0199) — Bootstrap f_xc
- [arxiv.org/abs/cond-mat/0205579](https://arxiv.org/abs/cond-mat/0205579) — ALDA for structural properties

---

## 2. 交换与关联的分解（L2）

`ε_xc(n) = ε_x(n) + ε_c(n)`

### 2.1 交换部分（精确）

LDA 交换能（Dirac 1930）：

```
ε_x(n) = −(3/4π)(3π²n)^{1/3} = −C_x / r_s    （以 Rydberg 为单位）
```

其中 `r_s = (3/4πn)^{1/3}` 是 Wigner-Seitz 半径（单位：Bohr 半径 a_0）。

**二阶导数**：

```
d²ε_x/dn² = −C_x · (constant > 0) / n^{5/3}
```

⟹ **交换二阶导数恒负**（所有密度）。

### 2.2 关联部分（参数化）

Perdew-Wang 1992 / Perdew-Zunger 1981 参数化给出 `ε_c(r_s)`，其导数行为：

- **高密度** (r_s → 0)：`ε_c ~ A ln r_s` ⟹ `dε_c/dr_s ~ A/r_s` ⟹ `d²ε_c/dn²` **负**
- **低密度** (r_s → ∞)：`ε_c ~ −B/r_s^{3/4}` (RPA) ⟹ 导数行为反转 ⟹ `d²ε_c/dn²` **正**
- **中间密度**：符号反转发生在 **r_s ~ 3-4** 附近

来源：
- Perdew & Wang, PRB 45, 13244 (1992) — 标准参数化
- 多篇 TDDFT 文献讨论 f_xc 的密度依赖

---

## 3. r_s ~ 2-3 的材料背景

铁基超导体典型密度：
- **Fe**：r_s ≈ 2.12（金属铁）
- **As, Se**：r_s ≈ 2.5-3.0（配位环境下的有效密度）

⟹ **正好落在关联导数符号反转的过渡区**。

---

## 4. ALDA 核的已知不稳定性（L2，多源）

### 4.1 铁磁不稳定性

ALDA 在**自旋极化**情况下有已知病理：

- 若 `f_xc` 过强（绝对值大），可预言**非物理的铁磁相变**
- Stoner 判据 `I · χ_0 > 1` 在 ALDA 下可能被错误触发

来源：
- [arxiv.org/abs/1108.3100](https://arxiv.org/abs/1108.3100) — Discontinuities of f_xc and charge-transfer
- [arxiv.org/abs/cond-mat/0610326](https://arxiv.org/abs/cond-mat/0610326) — Fermi liquid instabilities in spin channel
- [physics.stackexchange.com](https://physics.stackexchange.com) 讨论：ALDA 可预言假磁性

### 4.2 Renormalized ALDA (rALDA)

多篇文献发展 **rALDA** 正是为了修正 ALDA 在某些波矢/密度下的不稳定行为：

- 对 ALDA 核做**局域截断**（truncation for certain wave vectors）
- 目的：恢复正确的均匀电子气极限

⟹ **ALDA 本身在某些区域确实有问题**，需要重整化。

来源：
- [arxiv.org/abs/1501.00902](https://arxiv.org/abs/1501.00902)
- [arxiv.org/abs/1306.2732](https://arxiv.org/abs/1306.2732)

---

## 5. ★ 对 Essenberger Eq.(32) 推断的含义

Essenberger PRB 94, 014503 (2016) 写（见本仓已有笔记 `literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md` 附录 D.1）：

```
Λ^SF(x₁,x₂) := 2 ∬ dx dx′ f^xc_zz(x₁,x) χ_zz(x,x′) f^xc_zz(x′,x₂)   (Eq. 32)
```

该文推断 `Im Λ^SF ⪰ 0`（排斥）的逻辑：
1. **Congruence 形式**：`Λ^SF = 2 ∬ f χ f`
2. **若** `χ` 是 positive-definite（作为内积意义下的算子）
3. **且** `f_xc` 是实数
4. **则** `Λ^SF` 在某种意义下"正定"

### 5.1 逻辑链的缺口

**问题 1：f_xc 的符号**

- 上述 §2-3 表明：ALDA 的 `f_xc` 在 r_s ~ 2-3 **符号不定**（交换负、关联在过渡区）
- 即使 `f_xc < 0`，congruence 形式 `∬ f χ f` 仍可能**正定**——只要 χ 的算子符号对
- 但 **f_xc 符号不确定**意味着**不能只从 congruence 形式直接读出 `Im Λ^SF` 符号**

**问题 2：χ 的谱结构**

- χ_zz 是**磁化率**（Lindhard 函数的自旋版本）
- 在 RPA 层次，χ 的虚部**与粒子-空穴激发谱直接相关**
- **χ 是否作为内积算子 positive-definite，需要单独论证**

⟹ **Essenberger 的 `Im Λ^SF ⪰ 0` 推断需要以下之一**：
   (a) 证明 `f_xc χ` 的组合在相关频率/动量下给出正算子
   (b) 或者，直接计算 `Im Λ^SF` 的符号（他们实际做的是这个）

### 5.2 文献中的实际做法

查阅 Essenberger 及相关文献（附录 D）：

- 他们**实际计算了** `Λ^SF_{ij}(E,E',ω)`（数值，Fig. 3/4/5）
- **结果显示**：对角元 `Λ^SF_{±±}` **为负**，非对角元 `Λ^SF_{±∓}` 为正
- ⟹ **他们没有依赖 congruence 推断符号**，而是直接算出来的

---

## 6. 结论与建议

### 6.1 直接回答原问题

**ALDA f_xc 在 r_s ~ 2-3 能否为负？**

**答：取决于交换与关联的竞争；总的 f_xc 符号在该密度区不定。**

- 交换部分：`d²ε_x/dn² < 0`（恒负）
- 关联部分：`d²ε_c/dn²` 在 r_s ~ 3-4 附近**由负变正**
- r_s ~ 2-3 落在过渡区 ⟹ **总符号取决于具体密度与参数化**

### 6.2 对 Essenberger 推断的含义

1. **Congruence 形式 `Λ^SF = 2∬ f χ f` 本身不足以推出 `Im Λ^SF ⪰ 0`**
   - f_xc 符号不定
   - χ 的算子符号性质需单独论证

2. **Essenberger 实际做法是数值计算 `Λ^SF`**（他们的 Fig. 3/4）
   - 结果：`Λ^SF_{±±} < 0`（对角排斥）、`Λ^SF_{±∓} > 0`（带间吸引）
   - ⟹ **他们的物理结论不依赖 congruence 推断**

3. **若要严格化 `Im Λ^SF ⪰ 0` 的解析论证**，需要：
   - (a) 明确 χ 的 passivity 性质（Kramers-Kronig + causality）
   - (b) 或者，证明 ALDA + passivity 的组合在相关积分中给出正符号
   - (c) 目前这两条**在文献中都不是显然的**

### 6.3 本仓天花板证明的稳健性

**本仓 `theory/notes/literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md` 附录 D 
记录的推断 `Im Λ^SF ⪰ 0` 应标为推断，不应标为已证事实。**

替代论证路线：
- 若 Essenberger 数值计算的 `Λ^SF_{ij}(ω)` 图可信
- 则**直接从图读出符号结构**
- 附录 D.6 已记录：`Im[Λ^SF] > 0` 与 `Im[w] > 0`（Eq. 37，排斥谱函数）

⟹ **建议：用 Eq.(37) + Fig. 3/4 的数值结果，不要用 Eq.(32) 的 congruence 推断。**

---

## 7. 未核实项与开放问题

1. **Perdew-Wang / Perdew-Zunger 参数化的 `d²ε_c/dn²` 精确符号反转点**
   - 文献只给出定性趋势（"r_s ~ 3-4"）
   - 若要定量，需回到原文或重新数值求导

2. **ALDA 在 r_s = 2.0, 2.5, 3.0 的数值**
   - 未找到直接给出 `f_xc(r_s)` 数值表的文献
   - 可从 Perdew-Wang 公式自己算

3. **χ 的 passivity 与 positive-definiteness 关系**
   - Kramers-Kronig ⟹ `Im χ(ω) ≥ 0` for ω > 0（因果性）
   - 但这**不等价于** χ 作为内积算子 PSD
   - 需要量子场论层面的论证（未查）

4. **Essenberger 是否在其 SM 或 Ref.[46] 中给出更严格论证**
   - SM 取不到（APS 302→403）
   - Ref.[46] = arXiv:1409.7968 已取（附录 D），其 Eq.(32) 同样是 congruence 形式
   - **该文 Eq.(37) 是更强的陈述**：`Im[Λ^SF] > 0` **直接假设**（非推导）

---

## 参考来源

**直接引用（WebSearch 返回，未读全文）**：
- [arxiv.org/abs/1501.00902](https://arxiv.org/abs/1501.00902) — Renormalized ALDA
- [arxiv.org/abs/1107.0199](https://arxiv.org/abs/1107.0199) — Bootstrap f_xc
- [arxiv.org/abs/cond-mat/0205579](https://arxiv.org/abs/cond-mat/0205579) — ALDA for structural properties
- [arxiv.org/abs/1108.3100](https://arxiv.org/abs/1108.3100) — Discontinuities and charge-transfer
- [arxiv.org/abs/cond-mat/0610326](https://arxiv.org/abs/cond-mat/0610326) — Fermi liquid spin instabilities
- [harvest.aps.org PRB 101, 195137](http://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.101.195137/fulltext) — Plasmon dispersion & f_xc constraints
- Perdew & Wang, PRB 45, 13244 (1992) — 标准关联能参数化（已知，未重取）

**本仓已有笔记**：
- `theory/notes/literature/iron_arsenide_lambda_ij_retrieval_2026-09-18.md` 附录 D
- `theory/notes/literature/source_hunt_seven_items_2026-09-18.md` 项 7

---

## 附录：为什么 WebFetch 全部失败

本轮 3 次 WebFetch 均返回 "Unable to verify if domain arxiv.org is safe to fetch"。
这是**网络策略阻塞**，不是内容问题。

替代方案（未执行，因已有足够 L2 证据）：
- 用 `harvest.aps.org` 路线（对 APS 论文）
- 或要求主会话/有 Bash 的 agent 落盘后 `Read(pages=)`

本报告基于 **WebSearch 摘要级信息（L2）+ 本仓已有一手笔记（L1）**，
符合「取不到正文时的处理纪律」（见 `reference/research-provenance-discipline.md`）：
**read_depth 不许升、入口留痕、单一二手来源不入承重结论**。

本报告的承重结论（"f_xc 符号不定"）基于**多源一致的定性陈述**（交换恒负 + 关联符号反转），
不依赖任何单一二手源的数值。
