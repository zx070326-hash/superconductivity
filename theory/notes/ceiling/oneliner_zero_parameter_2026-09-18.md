# 零参数一行公式（方向 1 交付件）

建立：2026-09-18　执行：sonnet subagent（受 team-lead 派活）

---

## 0. 那一行公式

$$
\boxed{\;T_c \;\le\; C_\infty\sqrt{\eta/M}\;=\;C_\infty\sqrt{\lambda\langle\omega^2\rangle}\;},\qquad
C_\infty=\frac{1}{2\pi}\sqrt{\mathfrak g(2)}=0.1827262477\ldots
$$

**这是同一条定理的两种写法，不是两个定理**——第二个等号 `λ⟨ω²⟩ = η/M` 是一条**恒等式**（§1），
不是近似或第二个证明。规范来源是 `theory/derivations/P4_eta_over_M.md`（**定理 A**）§2.A，
它上游依赖已合并主线的 `theory/derivations/P4_main_theorems.md` 定理 1(c) + 定理 2 +
`theory/derivations/P4_g2_analytic.md`（𝔤(2) 的解析包围）。

**分支状态（2026-09-19 订正，原⚠已过期）**：`P4_eta_over_M.md` 所在分支 `etaM-corrected`（HEAD `ed6de52`）
已于 commit `9e1c83d`（「合并 etaM-corrected：η/M 定理（一行公式的规范来源）进入 p4-conj2」）**合并进 `p4-conj2`**；
`git merge-base --is-ancestor etaM-corrected p4-conj2` 核实为真，文件当前**就在本分支工作树里**
（`theory/derivations/P4_eta_over_M.md`），可直接 `Read` 复核，不再需要切分支或 `git show`。
下方引用「规范来源」处均按此更新状态理解；§7 末尾的流程问题（原「建议先合并该分支」）已解决，不再适用。

符号定义：

| 符号 | 定义 | 出处 |
|---|---|---|
| $\mathfrak g(2)$ | $\lambda_1(\mathcal O[1/k^2])$——算子 $\mathcal O[f]$（`P4_main_theorems.md` §0 定义）在检验函数 $f(k)=1/k^2$ 处的顶本征值。**这不是任何测度的矩**，是一个纯算子常数。 | `P4_g2_analytic.md` |
| $C_\infty$ | $(1/2\pi)\sqrt{\mathfrak g(2)}$ | 同上 |
| $\langle\omega^2\rangle$ | **测度 $P$ 的二阶矩** $\int\omega^2 P(d\omega)$，$P\in\mathcal P$ 是 $\mathbb R_+$ 上概率测度。等于 Allen–Dynes 的 $\bar\omega_2^2=(2/\lambda)\int\alpha^2F(\omega)\,\omega\,d\omega$（AD 1975 Eq. 5，$n=2$） | `P4_universal_tc_envelope.md` §1；`P4_eta_over_M.md` (1.3) |
| $\eta$ | $N(0)\langle I^2\rangle$（Hopfield 参数；McMillan–Hopfield，$I$ 是形变势矩阵元） | McMillan 1968 Eq. (22)；`P4_eta_over_M.md` §1.C |
| $M$ | 离子（或含权多组分求和后的有效）质量 | 同上 |

⚠ **CLAUDE.md §1「一个符号在链条里只许有一个意思」的核对**：本仓的 `M2`/「矩」有两个不等价读法——
`lim k²f(k)`（尾系数）与 `∫ϖ²P`（正测度的二阶矩），只在 $f\in\mathcal P$ 上重合。**本一行公式里的
$\langle\omega^2\rangle$ 走的是后一种读法**（测度 $P$ 的矩，见上表），不是任何核函数的尾系数；$\mathfrak g(2)$
则完全是另一个对象（算子顶本征值，不是矩）。三者没有混用风险，但读者若把 $\langle\omega^2\rangle$ 当成
「某个核的尾部系数」就会读错——本定理链条里 $P$ 全程是概率测度，从不对核做加减，所以不会撞上
「$\mathcal P$ 对减法不封闭」那条坑（那条坑只出现在多带 ceiling 的 $M2^{ph}-M2^{sf}$ 里，见 §4 (iv)）。

---

## 1. λ⟨ω²⟩ = η/M：恒等式的地位（不是近似）

在框架内（各向同性 Migdal–Eliashberg，$P(d\omega):=2\alpha^2F(\omega)d\omega/(\lambda\omega)$，
$\lambda:=2\int\alpha^2F(\omega)/\omega\,d\omega$，AD Eq. 3/5），

$$
\lambda\langle\omega^2\rangle = 2\int_0^\infty \omega\,\alpha^2F(\omega)\,d\omega
$$

是纯算术恒等式（无假设）。McMillan 1968 Eq. (22) 进一步给出这个一阶矩的微观形式
$\int\omega\alpha^2F\,d\omega = N(0)\hbar\langle\mathcal G^2\rangle/(2M)$（$\mathcal G$ 是电子形变势矩阵元，
McMillan 原话："the first moment of $\alpha^2F(\omega)$ is independent of the phonon frequencies"）,
于是 $\lambda\langle\omega^2\rangle=\eta/M$，$\eta:=N(0)\langle\mathcal G^2\rangle\equiv N(0)\langle I^2\rangle$。
`P4_eta_over_M.md` §1 把这条链拆成「定义 / 恒等式 / 近似」三层，我复核了这个拆分，判定合理：

- $P=2\alpha^2F/(\lambda\omega)$，$\langle\omega^2\rangle=\bar\omega_2^2$ —— **定义**
- $\lambda\langle\omega^2\rangle=2\int\omega\alpha^2F\,d\omega$ —— **恒等式**（算术）
- $\int\omega\alpha^2F\,d\omega=N(0)\langle I^2\rangle/(2M)$ —— **恒等式**（McMillan Eq. 22 的重排；固定 $M$ 下精确，依赖 Born–Oppenheimer 因子化）
- 「$\eta$ 与 $M$ 无关」——**近似级物理命题**，需要 (A2) Born–Oppenheimer 解耦 + (A3) 谐波/准谐波声子；非谐下 $\eta$ 随同位素替换漂移（PdH 实测 30%）

**这一层区分不影响 $T_c\le C_\infty\sqrt{\lambda\langle\omega^2\rangle}$ 本身**——那条不等式只用到 $\lambda,\langle\omega^2\rangle$
的定义，不依赖 (A2)(A3)。(A2)(A3) 只在把它**改写**成 $\eta/M$ 并赋予「$\eta$ 不随 $M$ 变」这句物理解读时才进场。

---

## 2. 严格性逐项标级

| 结论 | 级别 | 依据 |
|---|---|---|
| 定理陈述 $T_c(P;\lambda) < C_\infty\sqrt{\lambda\langle\omega^2\rangle}$，对一切 $\lambda>0$、一切 $P\in\mathcal P$ | **analytic** | 定理 1(c)（`P4_main_theorems.md` 引理 M 取 $G(k)=\langle\varpi^2\rangle/k^2$）；无计算机辅助的算子不等式链 |
| $\lambda\langle\omega^2\rangle=\eta/M$（固定 $M$ 的单组分） | **analytic**（在 (A1)(A2)(A3) 下为恒等式） | McMillan 1968 Eq. (22)–(23)，一手核对 |
| $C_\infty$ 的数值区间 | **rational-exact**——**不得简称 analytic** | 见下 |
| 材料表上的 $\eta/M$、天花板值、$\mathfrak f$ 比值 | **numerical-evidence** | 外部 DFT（$\lambda,\bar\omega_2$）+ float64 算术，无区间化 |

**$C_\infty$ 的格数与有理数个数（CLAUDE.md §1 要求）**：`P4_g2_analytic.md` §8 给出四档证书，我独立重算了闭式端点：

- **CERT-B**（"手可列表"版本）：试探向量 $M=12$ 格（12 个正整数分量：`10000, 2286, 704, 275, 129, 70, 41, 27, 18, 13, 9, 7`）+ 尾部分割点 $N_1=40$；推导链上显式出现的有理数除 12 个分量外还有 $\|v\|^2$、$\rho$（26 位分子的精确有理数）、$W$、$T(40)$、$r^2$、$\rho-13/25$ 共 6 个中间有理量，合计约 **18 个有理数**。给出
  $$1.3181093293757127 \le \mathfrak g(2) \le 1.3181615784372800\quad(\text{相对宽度}\ 3.964\times10^{-5})$$
- **CERT-D**（脚本生成、本文引用的数值端点）：$M=200$ 格、$N_1=500$，同法生成 200 个整数分量 + 同一组中间有理量，宽度收窄到 $4.634\times10^{-10}$：
  $$1.318140211762671 \le \mathfrak g(2) \le 1.318140212373478$$

我用这两组端点各自算了一遍 $C_\infty=\sqrt{\mathfrak g(2)}/2\pi$（Python，`math.sqrt`，仅用于本次核验，不承重）：

```
CERT-D:  C_inf ∈ [0.1827262477467, 0.1827262477890]
CERT-B:  C_inf ∈ [0.1827241072097, 0.1827277287107]
```

CERT-D 上端 `0.1827262477890` 与仓库到处引用的 `0.1827262477` / `0.1827262478`（十三位口径）一致，**独立复现**。
**界的常数用上端**（CERT-D 上端 0.182726247790）才是合法的上界常数；下端只用于陈述常数本身的位置。

---

## 3. 限定条件（本任务重点）

从规范来源（`P4_eta_over_M.md` §1.D 假设清单、§2.A/§6 诚实边界、`P4_universal_tc_envelope.md` §3、
`ceiling/multiband_ceiling_CANONICAL_2026-09-18.md`）汇总全部限定条件，然后判定承重与否。

**全部限定条件清单**（不筛选，先列全）：

1. 单边界方向：这是**上界**，不是 Tc 的估计，也不保证可达
2. 适用核类 $\mathcal P$：概率测度，有限二阶矩；$\gamma>2\cap\mathcal P=\emptyset$（核尾恰为 $k^{-2}$ 或更慢）
3. $\mu^*$ 的地位：待厘清（见下，任务书要求的核实点）
4. 框架 (A1)：各向异性/多带效应不计入——单带各向同性 Migdal–Eliashberg
5. (A2)(A3)：Born–Oppenheimer 解耦 + 谐波/准谐波声子——只在把 $\lambda\langle\omega^2\rangle$ 读成「$\eta$ 与 $M$ 无关」时才需要
6. 多带扩展需要额外结构条件（(Z-dom)⟺(K1)/(K2)、GF 非对角因子 1.248）
7. $\eta$ 是反推量，非独立测量（材料表专属，不影响定理本身）
8. $C_\infty$ 只在 $\lambda\to\infty$ 渐近锐；有限 $\lambda$ 下比定理 B（$G(\lambda)\sqrt{\eta/M}$）松 10–90%

**判定：承重条款共 4 条**（不是 3 条——理由见每条末尾；删掉任何一条都会让读者把这行公式用错）：

### （承重 1）单边界方向

上界，不是 $T_c$ 的估计，也不是充要条件。$T_c$ 达到 $C_\infty\sqrt{\eta/M}$ 需要 $\lambda\to\infty$ **且** Einstein 谱
（$\delta$ 峰）；真实材料的谱有展宽、$\mu^*>0$，所以真实 $T_c$ 总是显著低于这个天花板（§4：LaH10 上界比理论值高 1.44–1.78×）。
反向不成立：$\eta/M$ 大不蕴含 $T_c$ 高；唯一严格的逆命题方向是 $\eta/M$ 小 $\Rightarrow T_c$ 低（(A1)–(A4) 框架内）。
**判定必要**：不写这条，读者会把 267–388 K 的天花板数字当成"预测"，而不是"必要条件"。

### （承重 2）μ* 的地位——两句话必须分开，否则会读成矛盾

任务书要求核实这一点，**我核实后判定：两个来源不矛盾，但极易被合并成一句假话**。

- `P4_universal_tc_envelope.md` 标题写明 `(μ*=0)`，`P4_main_theorems.md:141`："全部结果在 μ*=0"。
- 但 `P4_main_theorems.md:76` 推论 (v)："μ*≥0：定理 1(c) 与定理 2 的上界对库仑束逐字成立（$K_c\succeq0$ 算子序）"——
  即用 **μ*=0 的 Einstein 值**去界**任意 μ*≥0 的一般谱**，这一步是 analytic，对全体 μ*≥0 成立。
  `P4_eta_over_M.md` (A6) 正是引这条给出「定理 A μ\*-无关」的陈述。
- 该推论同一处**明确排除**了另一件事：它**不蕴含** μ*>0 下的 Einstein 极值性 $Q_\mu(P)\le Q_\mu(E)$（即"在同一个 μ*>0
  的世界里，Einstein 谱仍是最大化 $T_c$ 的谱"）——后者在弱耦合区（$\lambda\lesssim(1.77+3.76\mu^*)\mu^*$）**已被精确算子反例证伪**
  （`theory/notes/coulomb/mustar_extremality_weak_coupling_failure_2026-09-17.md`，本机 memory 同名条目复核过）。

**结论（不是来源间的真矛盾，是同一份来源里两个范围不同的陈述，合并转述会出错）**：
一行公式 $T_c\le C_\infty\sqrt{\eta/M}$ 对**一切 $\mu^*\ge0$** 成立（因为它只用到推论 (v) 那个更弱的"逐字成立"），
**但**"达到这个上界的谱是 Einstein 谱"这句只在 $\mu^*=0$ 严格成立；$\mu^*>0$ 时不能再说"最坏情形就是 Einstein 谱"。
**判定必要**：这正是派活消息要求核对的疑点，且真实存在混淆两句话的历史事故（同型混淆已在仓库别处发生过）。

### （承重 3）(A2)(A3)：η 与 M 无关只在谐波 + BO 解耦下成立

若只用原始形式 $T_c\le C_\infty\sqrt{\lambda\langle\omega^2\rangle}$，不涉及 (A2)(A3)——那条本身是 analytic、条件只是 §2 的框架。
但**一旦改写成 $\eta/M$ 并用来谈"为什么是氢""同位素指数 ½"**，就隐含用到 (A2) Born–Oppenheimer 解耦
（形变势与核质量无关）+ (A3) 谐波声子（$\omega_j(M)\propto1/\sqrt M$）。非谐下 (A3) 失效到"$\eta$ 漂移"的水平：
PdH 的 SSCHA 数据显示 $\eta(H\to D)$ 漂移 $\times1.296$（不是 $\times1$）。**这不破坏上界本身**（上界只需要 $\langle\omega^2\rangle$
的实际值），**破坏的是"同位素指数恰为 ½"这句解读**，而这句解读正是"为什么是氢"论证链条的一部分。
**判定必要**：不写这条，读者会把 PdH 的 $\alpha_{\rm true}=-0.38$（真实 Tc 同位素指数为负）误判为对本定理的证伪——
实际上两者谈的是不同对象（天花板的指数 vs 真实 $T_c$ 的指数，§3.C 的 (3.2) 分解），且 PdH 处天花板本身仍未被违反
（67.3 K vs 5.0 K，天花板高 13 倍，完全有空间容纳两者反号）。

### （承重 4）单带/各向同性框架——多带需要额外条件

一行公式的证明路径（定理 1(c)/2）是在**单一各向同性 Eliashberg 通道**（框架 (A1)）内建立的。
多带/各向异性推广（`ceiling/multiband_ceiling_CANONICAL_2026-09-18.md`）给出的是**另一条定理**
（$T_c\le C_\infty\sqrt{\rho(M2''_Z)}$），需要额外的结构性前提：
- $\alpha^2F^Z=\alpha^2F^{ph}+\alpha^2F^{sf}$、$\alpha^2F^\Delta=\alpha^2F^{ph}-\alpha^2F^{sf}$ 两片各自逐元非负；
- 物理前提 **(Z-dom)**（配对通道排斥不"倒灌"进 Z 通道——不是数学缺口，是对文献计算约定的核对：(K1) 加减型
  或 (K2) 单玻色子×形状因子 $g|F|\le1$，两条 analytic 的约定级充分条件；去掉它，天花板可被违反到 115×）；
- 非对角块额外乘常数 $c_{GF}=\pi^2/(6\mathfrak g(2))=1.2479$（引理 GF，analytic，MgB₂ 上代价 1.2%、LiFeAs 式 10.1%）。

**判定必要**：氢化物族（H3S、LaH10）本身在物理上是纯声子/单通道的场景（$\alpha^2F^{sf}\equiv0$，(Z-dom) 自动取等），
所以本一行公式**对本节 §4 的材料表直接适用、不需要多带条件**；但若把它套到多带材料（铁基超导等）上，
必须换成多带版本并核对 (Z-dom)，直接套单带公式**不是**已证的结论。

**未列为承重、但值得一提的次要限定**：
- $\gamma>2\cap\mathcal P=\emptyset$（核尾恰为 $k^{-2}$ 或更慢）——这是定理 2/推论 (ii) 证明机制内部的一个不变量，
  在本一行公式的**应用**层面不构成额外限制（$P\in\mathcal P$ 本身就是全部要求，读者不需要另外核对尾部指数）；
  它只在**修改**证明链条（例如尝试对核做加减）时才会咬人，本任务不涉及这类操作，故不升为承重项。
- "$\eta$ 是反推量、非独立测量"——只影响材料表的解读方式，不影响定理成立性，判定为文档级注记而非承重限定。

---

## 4. 材料表复核（0/48 独立重算）

数据源：`theory/data/ambient_hydride_candidates.csv`（Cerqueira et al. arXiv:2403.13496 Table I+II，
48 个常压氢化物候选，$\mu^*=0.10$ 的 Allen–Dynes 公式值，非实验值）。我写了一个新脚本从零重算
（不转引 `_build_ambient_candidates.py` 的打印结果，只读它写出的 csv 并独立复算 $\mathfrak f$）：

```
theory/notes/_scratch_ir/_verify_48_oneliner_2026-09-18.py
```

输出（本次独立运行，$C_\infty=0.1827262477890$，CERT-D 上端）：

```
n = 48
f: min 0.1334  max 0.4384  mean 0.2475
violations f>1: 0
omega2/omegalog: min 1.197  max 3.358  mean 1.831
```

**与 `theory/data/README.md` 记载的数字逐位一致**（$\mathfrak f\in[0.1334,0.4384]$，均值 0.2475，违反 0/48）——
独立复现成立。⚠ 这批 $T_c$ 是 **Allen–Dynes 公式值**（$\mu^*=0.10$），不是实验值；这批数据只检验
「AD 公式给出的 $T_c$ 与本一行公式相容」，不检验自然界。有实验值可比的点（H3S、YH6）在 `hydrides_theory.csv`：
$\mathfrak f_{\rm exp}\approx0.49$（H3S）、$0.528$–$0.583$（YH6），同样 $<1$，**未独立重算**（超出本任务的一天预算，
留给下一轮；标注为**尚待复核**）。

**紧度**：LaH10 四个有真实 $\bar\omega_2$ 的压力点上，界 A（$C_\infty$ 版）/ $T_c^{SSCHA}$ = 1.44–1.78×
（转引 `P4_eta_over_M.md` §5.2，本任务未独立重跑该表，因为四个数值点无法在本预算内重建完整 SSCHA 链；
标记为**转引未独立复核**，与 0/48 的独立复核区分开）。

---

## 5. 一句话：为什么是氢

$\lambda\langle\omega^2\rangle=\eta/M$ 把天花板写成电子量 $\eta=N(0)\langle I^2\rangle$ 除以核质量 $M$：在
Born–Oppenheimer 解耦下 $\eta$ 与 $M$ 无关（(A2)），所以天花板 $\propto M^{-1/2}$——**同样的电子结构**下，
把最重的原子换成氢（$M_{\rm La}/M_{\rm H}=138$）能把天花板抬高 $\sqrt{138}\approx11.7\times$；更一般地，
多组分分解 $\lambda\langle\omega^2\rangle=\sum_a\eta_a/M_a$（(A5)）里，重原子那一项被 $1/M_a$ 压制，
所以含氢化合物里**几乎必然是氢通道主导** $\eta/M$——不是因为氢的 $\eta_a$ 特别大，而是除以的 $M_a$ 特别小。
这是"why hydrogen"在本框架内唯一不依赖可调参数的形式。（"靠声子软化拿高 $T_c$"是假问题：软化会抬高 $\lambda$，
但 $\lambda\langle\omega^2\rangle$ 不变，天花板不动——这正是定理 2 的内容。）

---

## 6. 交付清单自检（CLAUDE.md §4）

1. 每条结论标级——见 §2 表。
2. 承重项自己从零重算——CERT-D→$C_\infty$ 的换算、48 候选的 $\mathfrak f$ 统计均独立重跑（§2、§4）；
   LaH10 四点紧度未独立重跑，已标注"转引未独立复核"（§4）。
3. 本笔记不是"订正"，是首次落盘，不适用该条。
4. 单边界行——§3 承重 1 已审方向；(v) 的方向（用 μ*=0 界任意 μ*≥0）已核。
5. 分支 UNMERGED（本笔记本身待用户放行）；工作树只新增本文件 + `_verify_48_oneliner_2026-09-18.py` +
   `INDEX.md` 一行 + 脚本入库（`git add -f`，见下）。
6. 本任务不是"收紧某个量"，不适用预算条款。

---

## 7. 给 team-lead 的回报摘要

- **最终一行公式**：$T_c\le C_\infty\sqrt{\eta/M}=C_\infty\sqrt{\lambda\langle\omega^2\rangle}$，$C_\infty=0.1827262477890$（CERT-D 上端）。
  两个写法是**同一条定理**（`P4_eta_over_M.md` 定理 A），$\lambda\langle\omega^2\rangle=\eta/M$ 是恒等式不是近似。
- **标级**：定理陈述 analytic；$C_\infty$ 数值 rational-exact（CERT-B：12 格/约 18 个有理数，宽 4.0e−5；
  CERT-D：200 格/同法，宽 4.6e−10，本文用它）；材料表 numerical-evidence。
- **承重限定共 4 条**（非 3 条）：①单边界方向；②μ* 地位（对全 μ*≥0 成立，但 Einstein 极值性只在 μ*=0 证过，
  两句不可合并）；③(A2)(A3) 只保证上界本身，"η 与 M 无关"这句解读非谐下会漂移（PdH 30%）；④单带/各向同性框架，
  套多带材料需换公式并核 (Z-dom)。
- **0/48 独立复核**：$\mathfrak f\in[0.1334,0.4384]$，均值 0.2475，违反 0/48，与仓库记载逐位一致（脚本见 §4）。
- **来源间矛盾**：**没有发现真矛盾**，但发现一处高危混淆点（承重 2）——"上界对 μ*≥0 成立"与"Einstein 是
  μ*>0 下的极值谱"是两句不同范围的话，仓库措辞容易把它们并成一句假话，已在 §3 承重 2 拆开写清。
- **一处流程问题（2026-09-29 复核：已解决——`etaM-corrected`（`ed6de52`）已并入 master，`P4_eta_over_M.md` 在 master 存在，`9e1c83d`、`ed6de52` 均为 master 祖先；以下为历史记录）**：规范来源 `P4_eta_over_M.md` 所在分支 `etaM-corrected`（`ed6de52`）**未合并进 p4-conj2**，
  当前工作树没有这个文件，只能通过 `git show` 读取——建议 team-lead 决定是否需要先合并该分支，
  否则这份一行公式笔记引用的"定理 A"在本分支不可直接复核。
