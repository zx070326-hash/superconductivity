# `P4_eta_over_M.md` §4 材料表 —— 主会话对 matcheck 复核的独立再复核（2026-09-18）

> CLAUDE.md §1：subagent 的「已完成」是证据线索不是结论。本篇不转引
> [`etaM_material_table_recheck_2026-09-18.md`](etaM_material_table_recheck_2026-09-18.md)（matcheck 交付）
> 的任何数字，也不读它的脚本；算子与算术全部按 `P4_main_theorems.md` §0 与 `P4_eta_over_M.md` §1.A 的
> 定义重写。
>
> 复现：`python theory/notes/_scratch_ir/_chk_etaM_table_mainsession_2026-09-18.py`（**68/71**，三个 FAIL
> 是我的阈值卡得比文档印出精度还紧，见 §1.2）、
> `python theory/notes/_scratch_ir/_chk_r_budget_jensen_floor_2026-09-18.py`（§3 的预算）。

---

## 1. 判决：matcheck 的「零错误」我确认，但覆盖面小于它自述的范围

### 1.1 我独立重算并确认的

| 对象 | 我的做法 | 结果 |
|---|---|---|
| §4.B 七行 η/M = λω₂² | 纯算术，CODATA 常数自行写入 | **7/7 相符** |
| §4.B 七行 η_H^ub（eV/Å²） | 自写单位链 `η = M c²E²/(ħc)²`，未看仓库换算代码 | **7/7 相符**（≤3e−3） |
| §4.B 七行 界 A = C_∞√(η/M) | C_∞ 取我自己 rational-exact 包围的上端 0.1827262477890 | **7/7 相符**（≤0.06 K） |
| §4.B 三行 A/exp | — | **3/3 相符** |
| §2.C 表 ϖ_c / G / G/C_∞（6 行抽样） | **自写 𝒪[f] = −L1+L2+L3 与 ℌ(ϖ)**，N=800，`brentq` 解 h(ϖ_c)=1/λ | **18/18 相符**（ϖ_c 相对差 ≤4e−7，λ=100 处 4e−6 是 N=800 vs 1200 截断） |
| §4.B 七行 界 B = G(λ)√(η/M) | 用我自己算的 G | **7/7 相符**（≤0.06 K） |
| §5.3 四档门槛 | — | **4/4 相符** |
| 方向哨兵 界A > Tc^SSCHA、界B ≤ 界A | 逐行 | **14/14 成立** |

**算子自洽哨兵**（这一条是我加的，用来证明我手写的 ℌ 不是抄的也不是错的）：
文档 §2.B 引理第 1 步断言 `ℌ(ϖ)_{00} ≡ ϖ²/(1+ϖ²)` 是**恒等式**（与 N 无关）。
我的实现在 ϖ∈{0.3,1,5} 上给出 0.082568807339 / 0.500000000000 / 0.961538461538，
与闭式差 <1e−14 ⟹ **L2 不含对角、L3 含对角** 这个读法被独立钉死（读反了这一格就不对）。

### 1.2 三个 FAIL 不是缺陷

`PdH η/M`：我 1008.42，文档印 `1008`（整数位）；`PdH η_H^ub`：我 0.2432，文档印 `0.24`。
两者都是文档的打印舍入，我的阈值（5e−5 / 3e−3）比印出精度紧。`λ=100 的 ϖ_c` 相对差 4.25e−6 是
我用 N=800、文档用 N=1200 的截断差，方向与量级都与文档 Stage 5 自述的 N 收敛一致。

### 1.3 我**没有**做的

- §2.C′ 的认证包围 `F_lo^E/F_hi^E` 本身（区间算术管线）——与 matcheck 同，只做了方向核对：
  我算的 float64 G(λ)√λ 落在文档区间内、且 `F_hi^E ≥` 点值（方向对）。
- AD 1975 Table III 的 Pb/a-Pb 锚（2.4 / 2.0 eV/Å²）原文未取。
- McMillan 1968 Eq.(19)–(24) 一手转录（见 [`../theoremI/Cinf_oneliner_independent_recheck_2026-09-18.md`](../theoremI/Cinf_oneliner_independent_recheck_2026-09-18.md) §3）。

---

## 2. matcheck 的两处发现：我逐条溯源确认

两条我都回到 `theory/data/isotope_curated.md` 亲读，**不是转引**：

1. **YH6 的 `Tc_thy=241.5` / `μ*=0.125` 是区间中点。** 源表第 29 行印的是 `236–247` 与
   `μ*(6Ry)=0.10–0.15`；(236+247)/2 = 241.5、(0.10+0.15)/2 = 0.125。**确认。**
   中点化本身在 `p4_tc_sandwich.py:245` 的注释里（`# Tc 取 236-247 中点`），但
   `P4_eta_over_M.md` 正文与表脚注都没有这一句。
2. **LaH10@163 GPa 的 `Tc_exp=250` 是 150 GPa 的实验点。** 源表第 8 行明写 `150（H，XRD）` 配 `250`，
   而理论行 `Tc_thy=242.8` 是 163 GPa 的 SSCHA+aniso-ME。**确认，是跨压力对照且未标注。**
   与 [`h3s_tc_exp_pressure_audit_2026-09-18.md`](h3s_tc_exp_pressure_audit_2026-09-18.md) 同形。

两条都**不影响上界成立**（界 A 对 YH6 是 318.5 K > 区间上端 247 K；对 LaH10@163 是 359.1 K，
对 250 K 有 44% 余量），只影响「紧度 / A/exp」这类比较性措辞。

---

## 3. ★ matcheck 没问的那一条：`r = ω₂/ω_log` 的来历，以及它的预算

`P4_eta_over_M.md` §4.A 只写「r 只有 LaH10 四个压力点有（`p4_tc_table.py LAH10_POINTS`，**P3 一手反演**）」。
matcheck 把这条核成「一手（P3 反演），逐位相符」——**逐位相符是真的，但「一手」二字掩盖了它反演自什么。**
我追到 `P3_tc_max_pressure.md` §2 / 可信度分层第 2 条，那里自己写得很清楚：

> 一手 ω₂ 列表经核实**确认不存在**（Errea 无 Source Data、Semenok SI 只有 λ/ω_log），
> 但从 Errea Ext. Data Table I **同给的 Tc^Mc 与 Tc^AD 两列封闭反演**（两式只差 f₁f₂）得 r。

⟹ **r 不是 ω₂ 的测量值，是用 Allen–Dynes `f₁f₂` 的确切函数形式 + 前因子 1.20 从两列 Tc 解出来的。**
它继承三个前提：(a) Errea 的 Tc^AD 用的就是这一支 f₂；(b) 前因子是 1.20（P3 用 `r>1` 反过来钉死它）；
(c) Tc^AD 与 Tc^Mc 出自同一份 α²F。P3 对这三条有披露（并自标「整套 r 反演依赖 AD f₂ 的确切形式」），
**但 `P4_eta_over_M.md` 没把这个依赖带过来。**

**方向审计（为什么这条比 §2 的两条更要紧）**：r 若被**高估** ⟹ `η/M = λ(rω_log)²` 高估 ⟹ 界 A 高估
⟹ 判决「界A > Tc^ME」是 **fail-open**。而这四行恰恰是表里**唯一不带 `*`**、被当作合法上界引用的四行。

**预算（CLAUDE.md「预算优先」）**：把 r 朝判据变难的方向退到它的**解析下确界** r ≥ 1
（Jensen：⟨ω²⟩ ≥ exp⟨ln ω²⟩，`analytic`，与任何 AD 约定无关）：

| 行 | r_反演 | 界A(r_反演) | 界A(r=1) | Tc^ME | 紧度(r_反演) | 紧度(r=1) | 判决@r=1 |
|---|---:|---:|---:|---:|---:|---:|---|
| LaH10 129 | 1.194 | 368.0 | 308.2 | 255.3 | 1.442 | 1.207 | **PASS** |
| LaH10 163 | 1.075 | 359.1 | 334.0 | 242.8 | 1.479 | 1.376 | **PASS** |
| LaH10 214 | 1.105 | 388.4 | 351.5 | 237.9 | 1.633 | 1.478 | **PASS** |
| LaH10 264 | 1.095 | 386.6 | 353.1 | 216.9 | 1.783 | 1.628 | **PASS** |

**结论（`analytic` 的部分 + `numerical-evidence` 的部分分开说）**：

- **判决对 r 的预算是 ∞**：把反演来的 r 整个删掉、只用 Jensen 的 r ≥ 1，四行**全部仍 PASS**。
  ⟹ 「界 A 是 LaH10 四点的合法上界」这个结论**不承 r 的重**，AD-f₂ 那三个前提垮掉也不影响它。
  这一步的逻辑是 `analytic`（r≥1 由 Jensen）+ `numerical-evidence`（四个数值比较）。
- **紧度数字承 r 的重**：`1.44–1.78×` 退到 `1.21–1.63×`，最大变动 19.4%（129 GPa 行）。
  ⟹ **凡引用「紧度 1.44–1.78×」处，必须同时写出 r 是从 Errea 的 Tc^AD/Tc^Mc 两列经 AD-f₂ 反演的。**
  只说「P3 一手反演」不够——读者会理解成有人测了 ω₂。

**建议的最小改动**（未擅自改 `P4_eta_over_M.md`，留给用户放行）：§4.A 那句加半行出处，
§4.B 的紧度句后加「（r 为 AD-f₂ 反演值；退到 Jensen 下确界 r=1 时紧度为 1.21–1.63×，判决不变）」。

---

## 4. 自由度复核（CLAUDE.md「先数自由度」）

matcheck 的 §4 已正确指出 η/M 与 η_H^ub 不独立、界A/界B 不独立。我另外确认一条它没算的：

**r 是不是循环量？不是。** Tc^Mc 在给定 (λ, ω_log, μ\*) 下**完全被 McMillan 公式决定**，零新信息；
但 Tc^AD = f₁f₂·Tc^Mc 中的 f₂ 依赖 r，故「给定 (λ,ω_log,μ\*) 后的 Tc^AD」与 r **信息等价**——
每个压力点确实带进一个新数。⟹ r 反演**不是**从 Tc 反推 Tc 的循环（与
[`hydride_csv_independent_audit_2026-09-17.md`](hydride_csv_independent_audit_2026-09-17.md) 的
`filling_factor_f` 不同形）。**但它是「另一个公式的输出」，不是测量**，这才是 §3 要说的事。

**与 Tc^ME 的独立性**：界 A 用 (λ, ω_log, r)，被界量 Tc^ME 是 Errea 的各向异性 Eliashberg 数值解，
与 AD 公式那两列是不同的计算 ⟹ 「界 > Tc^ME」**是一个真检验**，不是恒等式自证。

---

相关：[`etaM_material_table_recheck_2026-09-18.md`](etaM_material_table_recheck_2026-09-18.md)（被复核对象）、
[`../theoremI/Cinf_oneliner_independent_recheck_2026-09-18.md`](../theoremI/Cinf_oneliner_independent_recheck_2026-09-18.md)（同一链条的数学侧）、
[`h3s_tc_exp_pressure_audit_2026-09-18.md`](h3s_tc_exp_pressure_audit_2026-09-18.md)（同形的跨压力对照问题）、
`theory/derivations/P3_tc_max_pressure.md` §2（r 反演的原始出处）。
相关 memory：`[[p4-onesided-bound-direction-rule]]`、`[[hydride-csv-audit-filling-factor-circular-2026-09-17]]`。
