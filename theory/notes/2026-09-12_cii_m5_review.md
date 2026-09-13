# M5 对抗核验记录（segment C-II 最后 18 个对象的端到端区间化）

- 被审对象：`theory/derivations/p4_cii_m5_interval.py`
- 审查形式：reviewer 子 agent（opus）一轮对抗核验 + 本 implementer 的变异测试、区间类型清扫与守卫复验
- 分支：`worktree-agent-a8687e419d54121f8`（**未** merge 进 `p4-conj2`）
- 日期：2026-09-12/13

本文件记录的是**被挑战的东西和它们怎么被关掉的**，不是脚本说明。脚本本身的方法写在 `P4_einstein_extremal.md` §1.J 与模块头注释里。

---

## 0. 一句话结论

reviewer 提 5 项 MUST-FIX、11 项 CHALLENGED。**5 项必修全部关闭并逐条用运行验证**；C1–C9、C11 关闭，C10 改为在守卫自身的输出里如实打印。此外本会话自查又发现两个**同类**缺陷（(R11)，见 §6），其中一个是 **fail-open** 的，两个都不是任何检查发现的、只能由审计发现。

---

## 1. 五项 MUST-FIX

| # | 问题 | 为什么严重 | 关闭方式 |
|---|---|---|---|
| M1 | `chk("PI_IV encloses pi", _contains0(PI_IV - 4*iv_atan_nn(iv.mpf(1))))` | **空检查**：`PI_IV` 本身就由 `4*iv_atan_nn(iv.mpf(1))` 定义，左右两边是同一次调用的两个副本，恒真 | 换成对 **精确有理 Machin 括号** 的包含检查（`_pi_bracket`，纯整数/`Fraction`），另加 §4 的绝对锚点 |
| M2 | 交错级数余项 pad 无任何检查覆盖 | 删掉 pad、或把双侧 pad 改成单侧，**没有一个检查会失败**（N=20 时截断误差 ~2⁻¹⁶⁴，在 prec 103 下是亚 ulp，所有既有检查都看不见） | 新增 (A3)：在 `nterms=1..5` 共 25 个用例上，用**精确有理**部分和 `[S_n, S_{n+1}]` 验证 pad 是双侧的且支配截断。`iv_atan_nn` 加 `nterms` 形参专门给检查用 |
| M3 | 毒化守卫被 M3 的 `from ... import` 绕过 | `check_transcendental_unreachable` 替换 `M3.mlog` 等属性，但如果 M3 内部用 `from mpmath import log` 取了本地引用，替换不生效，守卫**静默变弱** | BAN/WATCH 改为 `(obj, ctx, nm)` 三元组，显式加入 `(M3,'M3','mlog')`、`(M3,'M3','digamma')`、`(M3,'M3','polygamma')` |
| M4 | BAN 条目缺失被静默跳过 | 名字拼错或上游改名 ⇒ 该原语根本没被毒化，守卫仍报「全部不可达」 | 新增 `missing` 列表并用一条 chk 报告；缺失即失败，而不是少毒化一个 |
| M5 | `check_constants_vs_m2` 只比宽度不比包含 | 常数可以整体偏移而检查照过 | 改为非对称包含 `lo <= ref*(1+1e-6) and hi >= ref*(1-2e-3)`，宽度另设一条 sanity floor（`width < 1.0`；见 §8 对这条容差的定级） |

---

## 2. 十一项 CHALLENGED

| # | 内容 | 处置 |
|---|---|---|
| C1 | σ-尾图表（`eta_tail_bj`、MEXP、σ⁵ 括号）**没有任何检查** | 新增 `check_tail_chart()`：两条精确图表恒等式各 108 个点，做区间相交判定。实测最差相对间隙 **4.54e-30 / 9.62e-28** |
| C2 | `1.0/B`、`1.0/y0` 用 float 倒数，可能向下舍入，在 B 上方留下无人覆盖的缝 | 改为向外舍入：`(iv.mpf(1)/iv.mpf(B)).b` |
| C3 | Peano 中点法则未计入**舍入后的中点**与精确中点之差 | `_cell_view` 显式包围位移 `E = m' - (a+b)/2`，余项加上 `w·|E|·sup|h'| + (ΔQ/2)·E²` |
| C4 | `iv_atan_nn` 对退化输入的断言是个脚坑 | 去掉断言，改为**端点分裂递归**：宽输入按两个端点分别算再取包络 |
| C5 | `theta_star_certified` 的尾 majorant 需要 `s0 >= 2`，但未强制 | 加 `assert s0 >= 2.0` |
| C6 | `_isect` 发现两个包围**不相交**时只是返回较窄的一个 | 加 `ISECT_HIT` 计数并在 `main()` 里 chk 为 0 |
| C7 | 未被执行的检查与通过的检查混在一起计数 | 加 `skip(name, why)`，摘要分别打印 `passed / failed / not exercised` |
| C8 | 源码守卫的行区间可能是空的或搞错了 | 加「区间非空且合理」（`hi-lo > 500`）的 chk，并把 teeth 探针移到区间**内部** |
| C9 | `find_wA2_wide` 的判据用了区间比较 | 改为 `_lohi(...)[1] > 1`（精确 `mp.mpf`） |
| C10 | 源码守卫只覆盖 M5 自己的 Sec. 1–7b，不覆盖 M3/M4 | **未修**，如实记入守卫自身的打印输出与模块头（见 §7） |
| C11 | 多处用 `ivmpf` 做比较 | 见 §6：这条后来被证明比 reviewer 提出时更严重 |

---

## 3. 两个守卫，以及必须原样保留的限度说明

守卫被拆成两个有标签的半边，因为它们的**声明强度不同**：

**(a) 源码守卫（`check_no_transcendental_source`）**——`tokenize` 扫描 M5 自己的 Sec. 1–7b（1329 行），断言其中没有 `mp.` / `iv.` 的超越函数调用。这是**纯文本缺席**的声明，比 (b) 严格更弱。它**不覆盖 M3/M4**（C10）。teeth 探针（一个故意放在区间内的 `iv.make_mpf@1130`）证明扫描范围确实包含该区间，而不是扫了个空集。

**(b) 可达性守卫（`check_transcendental_unreachable`）**——毒化 16 个原语后重算判据。它的限度必须原样陈述，不得软化：

> **it demonstrates unreachability along the path exercised with the inputs used, not a static proof for all inputs.**

不要把这句改写成 "verified unreachable"。被行使的路径确实包含产出最终 ϖ_A″ 的那次 `bounds3()` 调用，但这仍然是**动态**证据。

两个无法毒化的例外，守卫自己打印出来而不是假装没有：`iv.pi`（是常量不是函数，无法包装）、`M3.Z3H_IV`（import 时就由 `mp.zeta` 算好并套了 1e-13 的 pad）。

---

## 4. 空守卫规则（M1 的一般化）

M1 那条检查之所以是空的，不是因为写错了，而是因为它的**参照物和被测物出自同一条代码路径**。用协调者的原话记录这条规则：

> **a guard whose reference is produced by the same code path, or the same arithmetic type, as the thing under test is not a guard — it is a consistency check wearing a guard's name.**

同一条规则也解释了 M3 原来的 `theta_star_iv` 为什么自欺（§9）：它用**二分宽度**当 pad，而真正的误差源是 dps-40 → float64 的转换。**pad 必须覆盖它所替代的那个误差源。**

---

## 5. (R9)：arctan 验证集的同质性盲点

### 5.1 发现

M5 的 arctan 门原本靠三条恒等式（余角、倍角、Machin）+ 一次 700 位参照。同侪 agent 2-b 指出：**三条恒等式在 atan 下都是一次齐次的**，而 `PI_IV` 自己也由 `iv_atan_nn` 构造 ⇒ 整个恒等式集在 `atan → C·atan` 下不变，**结构性地看不见**回乘因子 `2**k` 的错误。M5 的 `iv_atan_nn` 最后一行正是 `return s * (2 ** k)`。

### 5.2 变异实验（每个变异体一个 fork 子进程，`PI_IV` 从变异体重算）

`PI_IV` 必须从变异体重算：真实 bug 里 π 也会由带 bug 的代码算出来；只替换 `iv_atan_nn` 会留下一个由原版构造的 `PI_IV`，让 Machin 检查**因为错误的理由**抓到变异体。协调者对这点的定级：**Recomputing `PI_IV` from the mutant was necessary, not thorough.**

修复前（`m5_atanmut.py`）：

| 变异体 | 被谁抓到 |
|---|---|
| `scaleback_x2` (`2**(k+1)`) | **没有任何恒等式**；只有两条 mp.atan 参照检查 |
| `scaleback_half` (`2**(k-1)`) | 同上 |
| `sqrt_1pu` | 恒等式 |
| `den_2n+2` | 恒等式 |
| `n4_norem` | 恒等式 |
| `n20_norem` | **正确地通过**（误差 ~2⁻¹⁶⁴，prec 103 下亚 ulp） |
| `red_4.0` | 宽度检查 |

修复后（`_m5_mut2.py`，覆盖修好的门）：

| 变异体 | 失败检查数 | 备注 |
|---|---|---|
| `baseline` | **0**（12 通过） | 应当通过、且确实通过的变异体是**承重证据** |
| `scaleback_x2` | 6 | 含绝对锚点、其负控制、精确有理 Machin |
| `scaleback_half` | 6 | 同上 |
| `rem_deleted` | 1 | 新增的 pad 双侧性检查——修复前**无人能抓** |
| `rem_onesided` | 1 | 同上 |
| `den_2n2` | 8 | |
| `sqrt_1pu` | 7 | |
| `n8` | 1（仅宽度） | 可接受：仍是合法包围，只是更宽 |

### 5.3 协调者对配置的修正（原话）

协调者先前说恒等式验证「strictly better than a reference test」，随后**撤回**：

> Appealing to no external value is what makes a homogeneous identity set unable to see `atan → C·atan`. The right configuration is not either/or: **identities plus a rational absolute anchor.** ... Keep both; don't drop the 700-bit reference check either, now that you know it was carrying the scale claim alone.

一般形式：

> **first ask what group your validation set is invariant under, then ask what it can detect.**

以及：**The anchor set has to span the parameter the defect lives in.**

### 5.4 锚点的构造与三重验证

锚点用的是交错级数的相邻部分和括号：对 0 < u < 1，

    u − u³/3  <  atan u  <  u − u³/3 + u⁵/5

两端都是有理数，**不引用任何超越函数、库原语或高精度参照**。比较用 `_cmp_ratio` 做纯整数运算（`_dyadic` 读 `x._mpf_`），所以比较本身是精确的。锚点取 **dyadic** 的 3/4, 1/2, 1/4, 1/8——都严格大于 `_ATAN_RED = 1/16`，因此**减半循环、级数与回乘三者被一起钉住**；用 dyadic 而不是 2-b 建议的 9/10，是为了让 u 本身在二进制下精确可表示。

三重验证（`m5_anchor.py`）：(i) 括号在四个点上都包含 700 位真值（宽度 1.18e-01 降到 6.10e-06）；(ii) baseline 通过；(iii) 两个回乘变异体**全部触发**（8/8 拒绝）。全套锚点扫描耗时 0.9 ms。

### 5.5 一个顺带的定量事实

arctan 的输出宽度在 **N=12 起饱和于 1.8e-30**（prec 103 的 ulp 量级）。N=20 是安全余量，不是必需；这也解释了 `n20_norem` 为什么应当通过。

---

## 6. (R11)：偏序当全序用——本会话自查发现的两个缺陷

### 6.1 规则

> 区间携带的是**偏序**；`max`、`min`、`sorted`、`bisect` 以及 `>` 的布尔用法都假定**全序**。`ivmpf > mpf` 在两者重叠时返回 **`None`**，`None` 是 falsy，于是 `if a > b` 静默走 else 分支，`max(...)` 可能返回一个**非最大元**。

### 6.2 两个实例

1. **`fam_sups`（fail-open）**：`out.append(max(max(val[ab] for ab in boxes), tl[k]))`，其中 `val[ab]` 是**非退化** `ivmpf`。返回值是族 sup，在 `tbar_tail_bound` / `prof_tail_bounds` 里当**上界**用 ⇒ 返回非最大元就是**上界被低估**，属于 fail-open。
2. **`tbar_bb`（仅松紧）**：`top <= lo_best * (1 + rel)`，右边在向外舍入下会宽出一个 ulp ⇒ 退出测试可能变 `None`-falsy。因为该循环从不剪枝（见 §7），只影响松紧不影响可靠性。

### 6.3 修法与理由

两处都是**恢复全序**——循环内部一律用 `mp.make_mpf(z._mpi_[0/1])` 取精确 `mp.mpf` 端点——而**不是加一个容差**。加容差会把一个类型错误糊过去。`fam_sups` 出口再用 `_iv_pt` 转回退化区间，因为 `tbar_tail_bound` 里的 `M[0]*g[1]/2` 必须保持**向外**舍入，而不是 `mp` 上下文的就近舍入。

### 6.4 检测方式，以及为什么它必须是审计

**两个都不是任何检查发现的，而且必然如此**：它们的误差方向让界看起来**更好**。审计是句法的：走遍产界路径上的每一个 `max(` / `min(` / `sorted(` / `sort(` / 裸比较，逐个问「参数重叠时这个返回什么」。

**清扫结果：用 AST 枚举全文件 181 个排序/比较站点（`_scratch_ir/m5_r11_sweep.py`），逐个按操作数类型定级；两个缺陷已修，第三处 `_isect` 的宽度比较（`(x.b-x.a) <= (y.b-y.a)`，两个退化区间相减仍向外舍入、可能宽出一个 ulp）也一并改为精确 `mp.mpf`。清扫 181 站点，2 个缺陷，剩余 0。** 顺带把 `_abs_range` 的 docstring 订正（它声称返回 exact mpf，实际返回退化 `ivmpf`——正是这个说明把我引向 `fam_sups` 那处误用），并把它跨零分支的 `mp.mpf(0)` 统一成 `iv.mpf(0)`。

---

## 7. (R10)：`tbar_bb` 的剪枝对比与作用域

协调者接受了关于 `tbar_bb` 可靠性的论证，并要求记录下面这句：

> **`theta_star_certified` prunes, so its soundness genuinely depends on `best` being a verified lower bound. Same branch-and-bound skeleton; whether it prunes decides whether an untested exit path is a tightness risk or a soundness risk.**

`tbar_bb` 与 `fam_sups` 都做 `boxes = keep + new`，**从不丢弃盒子**，所以覆盖不变量恒成立、`max hi` 无论哪条退出路径触发都是合法上界；`theta_star_certified` 则真的做 `boxes = [ab for ab in boxes if val[ab] > best]`，是真剪枝。

关于 `tbar_bb` 探针覆盖范围，必须原样保留的句子：

> *the convergence exit is exercised on 2 of 6 keys at loose settings, both in the easiest (G) family; the cap/stagnation exit — the path that actually changed — is not.*

探针实测（`srv_m5_tb.py` / `tb_probe.log`，4 进程，`rel=3e-2, rel_prof=3e-3, s0=512`）：`TBAR[G]` 291 盒 / 565 次求值 / 523.8 s，`TBAR[UG]` 208 盒 / 399 次 / 477.1 s，两者都由**收敛测试**退出（上限 900 未触及），结果 OK 0 bad of 2。

**后记（smoke 运行，协调者指出）**：上面那句在 smoke 之后需要更新——**cap 退出路径已被行使 2 次**：`TBAR[kap]` 与 `TBAR[phi]` 在 smoke 设置（`s0=64, tbrel=2e-1, tbrelprof=2e-2, --procs=20`）下都以 `boxes=900` 退出，分别得 `[0.45343049, 2.50756953]`（1790 次求值，1476.9 s）与 `[0.47233251, 5.13818246]`（1790 次求值，2660.5 s）——**无用地宽、但合法**的区间，与覆盖不变量论证的预测一致。这是论证被观察确认而不只是被阅读确认。限度：smoke 容差下的 cap 退出**不能**说明生产容差（`tbrel=1e-2, tbrelprof=1e-2, s0=4096`，见下一段）下 cap 不会触发；生产运行的日志逐 key 记录 `boxes=`，§11 按实际情况填。

**盒预算不是时间预算**：上面两个 key **花掉同样的 900 个盒子、同样的 1790 次求值，墙钟却差 1.4 倍**（1476.9 s vs 2660.5 s）。`max_boxes` 限的是盒数，而单盒成本随族（G ≪ κ < φ）和随盒内 profile 求值难度变化，所以调高 `max_boxes` 的墙钟代价**不能**按盒数线性外推，必须按最贵的族算。另需分清两个旋钮：`--maxc` 是 **CV 外层格子上限**（`cv_quad` 的 `max_cells`，只影响 12 个 pad 常数），`tbar_bb` 的 `max_boxes=900` 是函数默认值（本次加了 `--maxbox` 开关，生产运行保持 900）——即 smoke 里触发 cap 的那个上限与 `--maxc` 无关。

**生产容差的一次反转（协调者，发车 12 分钟后叫停并重发）**：`tbrelprof` 由 1e-3 改回 **1e-2**。理由是上面的数字：预期 6 个 key 都会骑在 900 盒的 cap 上，而**骑 cap 的 key 的宽度由 cap 决定、不由容差决定**——容差从未达到，这正是 cap 的含义。于是更紧的内层容差只是花 20 倍单次求值成本买来更窄的逐盒区间，然后被 cap 整个扔掉。协调者的原话：*Running a tight inner tolerance into a binding outer cap is the worst of both: you pay for the precision and the cap discards it.* 再加上上限的度量（§10b）：6 个 sup 合计只占 0.116 预算的 17%，即整个 TBAR 半边最多争得 ~0.015 的 ϖ_A″。**两个定价放在一起看**（不是对谁的修正）：按 smoke 实测的单次求值成本外推，`tbrelprof=1e-3` 时 TBAR 半边单独就要 15–30 h，压倒一切；`1e-2` 时它约 2.4 h，低于 CV 半边的 ~2.8 h（后者由 `--maxc=1500` 下的 φ 族一项主导，~1.7 h）。这个从「压倒一切」到「不再主导」的落差，正是「把紧的内层容差灌进会绑定的外层 cap」这条论证所预言的——论证成立；只是容差一旦不再相对于 cap 荒谬，它就不再起作用。被叫停的那次运行（`m5_final_tbrelprof1e-3_aborted.log`）尚在单进程的准备段（534 s，M4 常数刚算完），没有任何 M5 对象开始计算，无损失。重发命令：`--procs=32 --rel=5e-2 --tbrel=1e-2 --tbrelprof=1e-2 --maxc=1500 --dump=m5_final.json`。`--procs` 未按建议提到 40：机器不受内存（334 GB，可用 311 GB，每 worker RSS ~110 MB）或磁盘（worker 不写盘，dump 只由主进程写、KB 级）限制，留在 32 的唯一理由是用户的原始约束「44 核别开满，留一半」——32 已经越过它，不再往上。协调者的记录里没有这条、而有一条更晚的「允许多用几个核」；两条用户指令的出入由协调者提交用户裁定，裁定前不动。

---

## 7c. 交付表的第三列：退出原因在退出处记录，不从计数倒推（(R10)/(R13) 落到交付物上）

协调者从上面的估算里推出一条对交付物的新要求：κ 与 φ 这两个占掉整个预算的族预期在 `--maxc=1500` 处骑 cap，于是**它们的宽度是 cap 定的，不是容差定的**——把一个 cap 定的宽度拿去和「要求的容差」并排比较，是拿一个没人选过的数去比一个我们指定过的数。所以逐对象表要三列而不是两列：**实际宽度、要求的容差、哪条退出触发**；退出原因要**在触发处记录**，不许事后从格数/盒数倒推。

实现（`p4_cii_m5_interval.py`，md5 `da8dfb5e`；`533fa042` 是同一逻辑、少一段注释的前一版，只用于本机牙）：
- `cv_quad` 的 `why` 在触发的那个分支里赋值——`tolerance met` / `cell cap N` / `no cell cut` / `no cell split` / `round cap N`；容差与 cap 同时成立时写 `tolerance met (cells=N also >= cap N)`（宽度是容差定的，但 cap 对下一轮已绑定，两件事都说）；即时打印 `cv exit: ... after round R, cells=N`；存入 `QuadInfo.exit` 与模块级 `EXITS['cv:<族>']`。
- `tbar_bb` 原本就即时打印三条退出（`converged` / `box cap N` / `round cap 40`），现在同样写入 `EXITS['tb:<key>']`。
- `dump_constants` 把 `EXITS` 与 `SETTINGS`（实际使用的 rel / tb_rel / 两个 cap）与端点一起写盘；`load_constants` 读回，缺失时在日志里说明。
- 新的 `report_objects()`：18 行表，列为包围 / 实际相对宽度 / 要求容差 / 退出，末尾一行 `cap-set: k of n recorded exits (...)`。放在 `report_ratio()` **之前**打印——表是构建产物，不依赖装配，`bounds3` 抛异常不能把表一起带走。dump 里没有记录时每行印 `n/a (dump predates exit logging)`，**不猜**。
- 牙（`_scratch_ir/m5_exit_chk.py`，走真实函数）：松容差 → `tolerance met`；cap 低于基础格数 → `cell cap 3`；松 `tbar_bb` → `converged`；盒 cap 低于基础盒数 → `box cap 5`；dump→清空→load 后 `EXITS`/`SETTINGS` 逐项相等；表 18 行且各行带记录的字符串；去掉两个键的旧格式 dump → 18 个 n/a。RESULT OK。打桩的 `build_all` 管道（`m5_tinybuild.py`）仍通过（桩不记录，打印 `exit: ?`，正确）。

一个顺带的小证明，按协调者要求写进了 `cv_quad` 的 `no cell cut` 分支（(R4) 类陈述：不是「未测试故可疑」，而是「在陈述的条件下可证不可达」，条件随证明一起给出）：**容差未达时 `no cell cut` 不可达**。条件：`done` 为假，即对某个 kind 有 `sl > rel·E`（E = `max(est, 1e-300)`，与 `thr` 用的是**同一个**加了地板的量）。于是 `sl` 是 c ≥ 1 个正 pw 之和（c = 0 给 sl = 0，与 sl > rel·E ≥ 0 矛盾），`n_i = max(8, c) ≥ c`，最大 pw ≥ sl/c ≥ sl/n_i > rel·E/n_i，而阈值是 0.3·rel·E/n_i；1/0.3 = 3.3 的裕量吸收浮点求和误差（≤ 30000 项、eps 1e-16）。V 情形同理。分支保留并在触发时记日志：证明是关于**这个** `thr` 的，下一次改 `thr` 可能使它失效。所以实际可达的退出是容差、cap、`no cell split`（格子到浮点分辨率）、轮数上限四种。退出日志使这个证明对**报告**不再必要；它对**代码**仍然是一条有价值的陈述。

**关于本次运行的 CV 三行本可如何重建（记录，不采用）**：cap ⇔ 末轮 `cells ≥ 1500`；否则若末轮打印的四个松弛都明显 < rel ⇒ 容差退出；否则 ⇒ `no cell split`（`no cell cut` 已证不可达）。除打印舍入恰好落在 rel 上的刀锋情形外是精确的——但精确的重建仍是重建，交付表不用它。

**对正在跑的生产进程的影响与裁定**：它是从 md5 `82cca55f` 的脚本 import 的，6 个 sup 的退出会即时打印（那段代码早就在），但 3 个 CV 族只有逐轮 `cells=`/松弛行、dump 里没有 `EXITS`。两个选项（A：杀掉重发；B：让它跑完，CV 三行以「从轮日志重建」标注、dump 永远缺记录、交付脚本与产出 artifact 的脚本 md5 分裂）提交协调者，裁定 **A**，但决定因素不是那 40 分钟：**B 交付的脚本不是产出 artifact 的脚本**——这与 layer-two-b 曾不得不用字节级 `SEC.` 块 diff（`42a58b6`→`adb4c7f`）事后补救、并在其 §4.2 留下永久 caveat 的是同一种出处缺陷；在 §8 要求「每个认证常数出自一条命令」的里程碑里**明知故犯**站不住。次要代价同样真实：没有 `EXITS` 的 dump 在任何未来的 `--load` 上都答不出「这一行是哪条退出产生的」，这个洞不会愈合。

**换脚本前的门**（协调者，零墙钟成本——用生产作业没占的核）：新脚本先在最小设置下**端到端**跑完（`--procs=6 --rel=5e-1 --top=20 --maxc=50 --s0=8 --tbrel=5e-1 --tbrelprof=5e-1 --dump=tiny.json`），再 `--load=tiny.json` 一次，必须走到 `dump_constants → report_objects → report_ratio` 的末尾。理由是 §9b 那次崩溃本身：`--dump` 从未被行使、首次使用即失败、发生在末尾、代价 3.1 小时；`report_objects` 与 `EXITS`/`SETTINGS` 两个 dump 字段此刻同样从未被真实的 `build_all` 走到末尾行使过——打桩管道不是同一件事。**一条唯一职责是保存昂贵工作的路径，必须在昂贵工作存在之前被行使。** 门绿则杀 `82cca55f`、以 `da8dfb5e` 重发；门红则无损失，旧作业仍在跑。门的结果见下。

---

## 7d. 门的结果，(A) 的撤销，以及「一个误报自身参数的检查」

**门：绿。** `tiny_e2e.log`（pid 2936527，md5 `da8dfb5e`，`--maxbox=40`，2882.0 s）走完全程：`dumped 12 constants + 6 sups -> tiny.json` → 检查电池 → 认证原语安装与重建 → 毒化守卫 → **18 行对象表，末行 `cap-set: 9 of 9 recorded exits`** → `rigorous wA'' = 1520.136655` → 四个 ϖ 的 ratio 表（含 ϖ=1000，**未抛 `ComplexResult`**，sign clamp 报 6 次而不是静默）→ SUMMARY。新增路径无一崩溃。`--load=tiny.json`（第二半）单进程另跑，同样走到末尾。

**17 个 FAIL 逐条核对，无一是内部不变量**：7× `relative width < 1 (sanity floor)`（GU_G、GU_kap、GU_phi、V_W_phi、W_G、W_kap、W_phi）、6× `Sec.7 majorant ... does not bind`、1× `12 pad constants >= 0 ... sign clamp applied 6x`、3× `ratio.hi < 1 at w=10/45/1000`。没有分割断言、覆盖不变量断言、类型断言，没有 traceback。全部是在故意荒谬的宽度下不可能成立的**精度**断言（W_phi = [0, 1.76e15] ⇒ ratio = [−inf, +inf]）。

**(A) 被协调者撤销，理由是价格变了 3 倍。** 我给出 (A)（杀 `82cca55f` 重发）的报价时只有 G 族算完；到 07:09 生产日志已在 ~8100 s，**G（1164 s）与 κ（4487 s）都已落盘**在增量 `m5_final.json` 里，φ 在第 5 轮。重发要扔掉约 4500 s 的算术且**无法续跑**——`load_constants` 按我自己的设计拒绝部分 dump。于是 (A) 的真实价格是约 2 小时而不是 40 分钟，换来的是一个日志列。协调者的原话：「**一个以 40 分钟为输入的决定，在价格翻三倍时要重做。**」这一条记在这里不是作为对我的纠正，而是因为它是一条**可复用的规则**：批准是对某个价格的批准，价格变了批准就不再自动有效。

替代方案（现行）：让生产在 `82cca55f` 下跑完；交付表由**第二条命令** `--load=m5_final.json --exitlog=m5_final.log` 在冻结脚本下产出。dump→report 的切分正是 dump 的用途，§8 禁的是**手抄端点**，不是第二个进程读一个机器可读的文件。代价如实写在 §1.J：**artifact 的 md5 与交付脚本的 md5 不同**，并为此做两次出处核销（下）。

**一个真实缺陷，而且只有 tiny 跑才使它可见。** 6 行 `Sec.7 majorant for s>4096 does not bind` 打印的是模块常量 `TBAR_S0`，而该次运行用的是 `--s0=8`，紧挨着的数据行正确地写着 `(tail s>8: ...)`。检查的**消息里硬编码了一个它并不读取的参数**。这不是排版问题：它会把设置引起的失败归档成结构性失败、或反过来；而在生产设置下 s0 *恰好就是* 4096，这个缺陷将**永远不可观测**。协调者的定性，照录：**一个误报自身参数的检查是贴错标签的警报，不是无关紧要的字符串。** 修法：标签改读 `SETTINGS['s0']`（随 dump 一起走）；dump 里没有该键时印 `?(s0 not recorded in dump)` 而**不**代入默认值——**未知的参数与已知的参数不能印成同一个样子**。

同一条规则顺手修了第二处：`--load` 下 builder 根本没跑，`sign clamp applied 0x` 读起来像「从未触发」，实际是「未被行使」。现在印 `NOT EXERCISED here (--load replays post-clamp endpoints; the clamp lives in the builder) -- take the count from the build log`。这对交付物是实质性的：交付表现在**出自一次 `--load`**，协调者要的「生产宽度下的 clamp 触发次数」只能从生产日志读，表不该再引人往错的地方读。(R12) 的反面：**一个从未触发的检查与一个不可能触发的检查，在日志里必须长得不一样。**

### 7d.1 日志导出的退出：规则必须是全的，证据等级必须标出来

协调者对替代方案的要求：退出列从生产日志导出，但「不要写成 `cells > maxc ⇒ cap`」——要对 `cv_quad` 的四条与 `tbar_bb` 的三条各写出在打印行里的**判别签名**，并证明签名互斥且穷尽；两条退出共享签名就**标为未决**，不许猜；列要标 **log-derived**，不能和 dump-recorded 长成一个样。

实现为代码而不是散文：`--exitlog=FILE` → `exits_from_log()`。规则与拒绝点：

- **cv_quad**，取该族 summary 行之前最后一条 `round n: cells=N I-slack ... V-slack ...`，rel/max_cells/rounds 从两条 header 行读。松弛以 `%.1e` 打印，故打印值 s 只把真值框在 s·(1±0.051)（末位尾数的一半）。*确定达标* = 所有 s·1.051 ≤ rel；*确定未达标* = 某个 s/1.051 > rel；两者之间 → **未决**，绝不把打印舍入四舍五入成判决。其后：容差达标 / `cell cap`（cells ≥ max_cells）/ `round cap`（n == rounds）/ **stagnation**——这一桶同时点名 `no cell cut` 与 `no cell split`：两者都在下一轮打印之前 break，留下**完全相同**的最后一行。Sec. 4 证明了容差未达时 `no cell cut` 不可达，但那是**证明，不是日志证据**，而这个函数报告的是证据。
- **tbar_bb**，只有 `TBAR[k] = [lo,hi] (tail s>S: T) boxes=B evals=E` 一行——生产运行对 bb 的逐轮**不是** verbose，那些 hi/lo 行根本不存在。converged ⟺ hi ≤ lo·(1+rel)（两边相差在 1e-6 以内即印刷噪声 → 未决）。另两处拒绝：(a) 打印的 hi 是 `max(bb top, Sec.7 tail)`，当印出来的是尾界时 bb 自己的 top 只被它从上面框住，**非收敛**的判决因此未决（收敛的判决仍可靠：top ≤ hi）；尾界以 `%.3e` 对 hi 的 `%.8f`，故「尾界就是 hi」按 1e-3 相对判定——取两者中较粗的印刷精度，而不是一个永不成立的精确相等。(b) **`box cap` 与 `round cap 40` 共享签名**：轮数退出也可能带着 B == max_boxes 离开（切分被 cap 截断后 for-else 触发）。报为一桶并同时点名两者。它们只差在是哪个 cap，而表真正消费的那一栏——**容差定的 vs cap 定的**——依然是判定了的。

**签名穷尽性**：`cv_quad` 的循环只有这五个出口，且每轮恰打印一行；三个判别式（done / cells vs cap / n vs rounds）按固定次序作用在互斥谓词上，故互斥。`tbar_bb` 同理，三个出口两个签名，已如实标注。

**牙：`compare_exits()` 按「矛盾」而非「相等」计分。** 规则被有意做得比记录更粗，且**允许拒绝**——未决是正确输出，不是漏判；真正的缺陷是一个**给出了判决**却与退出分支自己写下的类别不同的标签。所以断言是「零矛盾」，精确命中数与拒绝数**并排另印**：「同意了」与「拒绝判断」是两种结果，不能加成一个数。对 `tiny.json` 的实测（先把日志里显式的 `cv exit:` / `bb exit:` 行**剥掉**，因为 `82cca55f` 的日志里根本没有这些行——不剥就成了同义反复）：**3/9 精确、4/9 同类更粗、2/9 拒绝、0 矛盾**。两处拒绝是 tb:G 与 tb:UG，即 `--s0=8` 下尾界恰是打印上端的那两行；s0=4096 时尾界远低于 hi，这种拒绝不应出现。另有 15 个变异用例走真实解析器，每个分支一个，全过（`_scratch_ir/m5_exitlog_chk.py`）。

### 7d.2 出处核销两次，不是一次

(i) **AST 层面的消去比较**（比协调者要的字节级 `SEC.` 块 diff 更强）。行级 diff 只能说「这些行变了，看起来像日志」。改为：解析两个文件，取 Sec. 1–7b 范围内每个顶层 def（各 36 个），删掉**恰好那些不可能影响取值的语句形态**——docstring；`log(...)`/`print(...)` 表达式语句；对 `why`、`EXITS[...]`、`SETTINGS[...]`、`self.exit` 的赋值；消去后变空的 `if` 包装；`QuadInfo` 的 `exit` 形参连同其默认值与调用点多出的那个实参——然后比较 `ast.dump`。结果：**36 个 def 逐个相同**，无单边缺失、无差异。负控制：把 `cv_quad` 细化阈值里的常数 0.3 改成 0.5，立刻被抓出。所以这次消去是有牙的，不是把什么都删光了（`_scratch_ir/m5_provenance_ast.py`）。
(ii) **正向的同一性测试**，diff 不是它。以完全相同的命令行（`--procs=6 --rel=5e-1 --top=20 --maxc=50 --s0=8 --tbrel=5e-1 --tbrelprof=5e-1 --maxbox=40`）在冻结脚本与 `82cca55f` 下各跑一次 tiny 端到端，比较 dump 里的端点元组是否**逐位相同**。协调者的理由：**diff 是对源码的一次阅读，逐位相同是来自算术的证据。** 结果见 §11。

**冻结**：协调者的 (b)——「对一个中间 md5 做的同一性测试是看起来像证据的陈旧证据，比没有更糟」——所以先停手、冻结、记 md5，再跑。第一次冻结 md5 **`019bc17db5d21a5e918687d5b5cb2e4a`**（`da8dfb5e` → 加 s0 标签修复、clamp 计数修复、`--exitlog`）。artifact 的产出脚本 md5 **`82cca55f`**。

**这次冻结随后被打破了，如实记。** §7e.3 的生产日志格式发现（82cca55f 打印 `bb exit:`；refusal (b) 的理由是错的）与 §7e.1 的 (R15) 拆分，都必须改在交付脚本里，md5 因此从 `019bc17d` 动到了后续值。按协调者自己的规矩，**在飞的那两个同一性进程测的已经是一个中间 md5**，它们的结果可以记录、但**不构成核销**。处理方式不是把它说成「差异反正只是报告层」就算了，而是：(1) 用同一个 `--all` AST 比较**机械地**给出 `019bc17d` → 当前文件的差异集——Sec. 1–7b **36 个 def 逐个相同**，全部差异只有 `exits_from_log`、`compare_exits`、`report_objects`、`main` 四个报告/CLI 阶段的 def 加两个新辅助函数 `_exit_class`、`rounds_bb_last`，`build_all` 不在其中；(2) 等所有修改落定后**重新冻结**，并对**最终提交的那个文件**重跑同一性测试。(1) 说明在飞的进程没有白跑（算术路径没被碰过），但**核销以 (2) 为准**——否则就正是那条规矩要禁的东西。

**中间 md5 的同一性结果（记录，非核销）**：`019bc17d` 与 `82cca55f` 的 tiny 端到端各跑 2300 s 后，`tiny_frozen.json`（complete=True，EXITS=9，SETTINGS=8）与 `tiny_82.json`（complete=True，EXITS=0，SETTINGS=0——82cca55f 按构造没有这两节，不是发现）的 12 个 CV + 6 个 TB + 6 个 TAILS 共 **24 个端点对，原始 mpf 元组逐位相同，0 处差异**（`_scratch_ir/m5_identity_cmp.py`，拒绝比较任何 partial dump）。这说明「消去论证删掉的东西确实对端点只写」这一点在算术上成立于 `019bc17d`；对最终文件的判决仍等 (2)。

### 7d.3 (2) 已完成：对最终冻结文件的正式核销

**最终冻结 md5 `3fdcfc7e30b4f8a5a0a97147c62d6986`**（本机文件与服务器上的 `p4_cii_m5_final.py` 两侧各算一次，相同）。**行尾歧义先说清**：这个值是 **LF 行尾**那串字节的 md5，也是 `git cat-file blob` 给出的值；仓库 `core.autocrlf=true` 且无 `.gitattributes`，同一文件在 Windows 工作树展成 CRLF 后是 `eb8945747d570cf97b99ddd9659b1617`。跑生产、跑两次同一性测试的都是 LF 那份。对它重跑同一条 tiny 命令行（两次运行的日志头逐项相同：`rel=0.5 top=20 max_cells=50 procs=6`、`s0=8 rel_prof=0.5 maxc=4000 max_boxes=40 rho0=1.35 -> 31 boxes`），与 `tiny_82.json` 比较：

```
A = tiny_82.json     complete=True  EXITS=0  SETTINGS=0
B = tiny_final.json  complete=True  EXITS=9  SETTINGS=8
24 entries compared, 0 differ
RESULT: IDENTICAL -- every endpoint tuple is bit-for-bit equal
```

**这一次才是核销**；中间那次 24/24 仍只作记录。两次都是 24/24 不能用来省掉这一次——「反正结果会一样」正是让陈旧证据活下来的那句话。

AST 侧同步重跑于**最终**文件：Sec. 1–7b **36 个 def 逐个相同、CLEAN**；`--all` 下差异仍只有 `check_tbar`、`dump_constants`、`main` 三个 def 加 5 个冻结脚本独有的报告阶段 def（`_exit_class`、`compare_exits`、`exits_from_log`、`report_objects`、`rounds_bb_last`），**`build_all` 不在其中**。读审计（`m5_read_audit.py`）CLEAN。

**顺带：两份 tiny 日志的检查区并排读，给了三条正面证据。** 两个脚本**失败同一组 17 条**（7 条相对宽 < 1 的 sanity floor、6 条尾 majorant 不绑定、3 条 `ratio.hi < 1`、1 条 clamp）——tiny 的宽度本就是 100%，这些是质量检查不是可靠性检查，逐条相同本身是「差异只在报告层」的旁证。三处差异恰好是宣称修掉的那三处，且**可观测**：

1. **s₀ 标签**：`82cca55f` 印 `Sec.7 majorant for s>4096`，而这次 tiny 跑的是 `--s0=8`；冻结脚本印 `s>8`。这正是 §7e 记的那个「检查误报自身参数」的缺陷，现在有了现场对照，不只有源码阅读。
2. **clamp 检查的合取拆分**：`82cca55f` 印成一条 `FAIL: 12 pad constants >= 0 (...) negative: -; sign clamp applied 6x`——**通过的那半**（端点无负）与**失败的那半**（clamp 触发 6 次）挤在一行共享一个判决；冻结脚本拆成两条，前者 PASS、后者 `FAIL: sign clamp on the 2 I-constants never had to fire  applied 6x`。
3. **检查总数 68 → 69**，增量正是这次拆分。


---

## 7e. (R15) 与 (R16)：拆分阶段之后欠的两笔账，以及我把它们记岔的地方

### 7e.1 (R15)：每条检查要申明哪个阶段能证伪它

协调者原话：*when you split a pipeline into build and report, every check must declare which stage can falsify it, and a check evaluated in a stage that cannot falsify it must print "not exercised" — never its default-valued result.*

本里程碑踩了一次，而且是在我刚刚「修好」clamp 计数之后还留着。符号 clamp 的检查是一个合取式 `not neg and SIGN_CLAMP[0] == 0`，两个合取项的证伪者不同：

- 「12 个端点 ≥ 0」——**两个阶段都能证伪**，端点正是 `--load` 回放的东西；
- 「clamp 从未触发」——**只有构建阶段能证伪**。clamp 活在构建器里，`--load` 下 `SIGN_CLAMP[0]` 还停在初值 0，于是这一项**因为不可能失败而通过**。

我先前的「修复」只是把说明文字改成「NOT EXERCISED here」，却仍把这一项折进 pass/fail 判定——**解释挂在旁边，判决还是错的**。现已拆成两条：端点照常 `chk`，「clamp 从未触发」在 `--load` 下走 `skip`。附带收益是另一个方向也正了：构建阶段的 0 现在是**挣来的真零**并按真零打印，于是「从未触发的检查」与「不可能触发的检查」不再长得一样（协调者：*a check that has never fired and a check that cannot fire must not look the same in the log.*）。

协调者另记一笔，值得抄进 §1.J：§8 的「一键复跑」条款**既是把我们推向拆分的力，也是造出这个隐患的源头**；该条款禁的是手抄端点，不是第二个进程去读机器可读的 artifact（主文档已在 §8 item 7 修订）。

并且**拆分不是「选择」的，是被迫的**：`report_objects` 在 82cca55f 里根本不存在（它只有 `report_ratio`），产出 artifact 的脚本**印不出**这张交付表。所以两个 md5 从一开始就不可能合成一个，(R15) 的分档记账是不可回避的，不是可选的谨慎。

### 7e.2 (R16)：消去论证的**前提**本身要审

负控制只证明「消去不是把什么都删光」。它**不**证明被删掉的四个名字对任何返回的界是**只写**的——而那恰好是「这只是日志」式论证会在比较前先删掉的东西。协调者要的读审计（`_scratch_ir/m5_read_audit.py`）结果**干净**：

- frozen 在 Sec. 1–7b 内对 `SETTINGS` 的 LOAD 次数为 **0**。被点名的 `SETTINGS['s0']` 只在 Sec. 8 `check_tbar` 出现一次，且只喂给 f-string 标签。
- 真正的取值源是 `TBAR_S0`，只有三个承载取值的读取点——`tbar_bb` 与 `tbar_tail_bound` 的默认实参、`build_all` 的 `s0 = (512.0 if quick else TBAR_S0) if s0 is None else s0`——**两个脚本逐行相同**，两者算出的尾不可能不同。
- `why` 的 4 次 in-range 读取全部流向 `EXITS[...]`、`log(...)`、`QuadInfo(...)`；`QuadInfo.exit` 全文件只被读一次，在 `if verbose:` 下的 f-string 里；`cv_quad` 只有一个调用者，取值路径是 `CV5.update(out)`，`info` 只在 `if verbose:` 内被碰。

**读审计顺带暴露了比较范围本身的一个洞，这才是这一节真正的产出。** `build_all` 计算着两个尾函数消费的 s0，却位于 Sec. 1–7b **之外**——一个承载取值的驱动函数从未进入过比较。把比较放宽到全部顶层 def（`m5_provenance_ast.py --all`）后：`build_all` **消去后完全相同**；全部差异只有 `check_tbar`、`dump_constants`、`main` 三个 def，加上 frozen 独有的 `exits_from_log`、`compare_exits`、`report_objects`，六个全在报告/序列化/CLI 阶段、全在常数产生之后运行。一般形式（协调者已收入主文档）：**一个只覆盖「取值在哪里算出」的比较，仍然必须覆盖「谁算出了这段范围的输入」。** Sec. 1–7b 是「什么产出界」的正确边界，却是「比较必须覆盖多少」的错误边界。

### 7e.3 解析器被拿去测了它不会面对的东西

§7d.1 里那个 3/9 精确、4/9 同类、2/9 拒绝、0 矛盾的记分，是**对 tiny 日志**量的；要解析的却是生产日志，而两者的 bb verbosity 不同。这与 s0 标签同类——**检查被拿去测了它不会面对的东西**——记分作废重做。重做时发现两件事：

1. **82cca55f 打印 `bb exit:` 的全部三个分支，只是没有 `cv exit:`。** 所以 9 个退出是**两个证据档次**：6 个 sup 是 **log-printed**（触发分支自己打的字符串，机器逐字读回），3 个 pad 族才是 **log-derived**（真正的重建，也是唯一可能返回「未决」的一档）。我先前把好的那一半也标成了坏的。`EXIT_SRC` 现在逐 key 记档次，表里多一列 `exit.src`。
2. **refusal (b) 不成立，而且我给的理由是错的。** 有了逐轮 `boxes`/`hi`/`lo`：循环在**每轮开头、打印之后**测 cap，所以盒上限退出发生在**其自身打印的 B 已达 max_boxes** 的那一轮；轮上限退出则要求最后一条打印行（round 39）的 B **仍小于** max_boxes、之后才分裂越界。汇总行报的是**分裂后**的盒数，这才是它分不开而轮次行分得开的原因。我原来的理由是「它们只差在哪个 cap，反正都是 cap-set」——**一个正确的拒绝若靠错误的理由支撑，比不拒绝更糟，因为理由才是读者会拿去复用的那部分**（协调者原话）。尾遮蔽同理消失：轮次行的 `hi` 是 bb 自己的 top，只有汇总行被 `max(top, tail)` 遮住。汇总行规则保留为兜底，两条拒绝只在那里为真。

牙改成三套共 31 项，全过（`_scratch_ir/m5_exitlog_chk.py`）：**A** 只有汇总行（确认两条拒绝在那里为真）；**B** 有轮次行（确认两条拒绝取消**且取消到正确的判决**——6 个 tb 行与 dump 记录逐字符相同——外加 5 个分支变异，与一条「轮次行压过矛盾的汇总行」的优先级测试）；**C** 直接对生产日志跑。C 的输出是这一节的落点：

```
cv:G   [log-derived] 'tolerance met'     tb:G    [log-printed] 'converged'
cv:kap [log-derived] 'cell cap 1500'     tb:UG   [log-printed] 'converged'
cv:phi [log-derived] 'cell cap 1500'     tb:kap  [log-printed] 'box cap 900'
```

把生产日志的退出行蒙掉后重建，仍复现全部判决，`box cap 900` 在内。最后这一条是决定性的：`kap w=0 bb exit: box cap 900 after round 6, boxes=900` 是在重测期间落盘的，**在旧规则下它会落进「两个 cap 都可能」的模糊桶**。取消 refusal (b) 因此在真实 artifact 上承重，不是假设性的。

**一个自伤的 bug，由蒙眼本身暴露。** 第一版蒙眼按「行里含 `exit:` 就整行删掉」做，而冻结脚本在 `TBAR[...]` 汇总行**末尾**追加了 `  exit: <why>`——于是蒙眼删掉了承载数字的那一行，套件 A 报「规则拒绝了」，真相是**测试毁掉了自己的输入**。现改为只截掉后缀，并断言 TBAR 行数与轮次行数在蒙眼前后相等。

### 7e.4 一条横幅，不是三处脚注

协调者：`m5_final.json` 的键只有 `CV5`/`TB5`/`TAILS`/`complete`，**没有 SETTINGS、没有 EXITS**，所以 s0 标签、9 个退出、clamp 计数**三样都是日志来的**，应当放在**一条**横幅后面，而不是散成三处脚注。`report_objects()` 现在印：三列是日志来的（点名）、**区间是 dump 来的**（原始 mpf 端点元组，无十进制往返）、以及 9 个退出里 log-printed 与 log-derived 各几个。

### 7e.5 tiny 对拍的两个副产品：clamp 的正面负控制，与一个硬编码的计数

两个 tiny 端到端（`019bc17d` 与 `82cca55f`，同命令行）跑完后，除了 §7d.2 记的 24/24 逐位相同，还给出两件事。

**(1) 符号 clamp 在 tiny 宽度下触发 6 次——production 的 0 因此是有意义的零。** §1.J.7 (vii) 说「生产宽度下 clamp 触发 0 次，且这个 0 是行使过的」，此前的依据只是 (R15) 的阶段论证（构建阶段能证伪它）。tiny 给了**正面**证据：同一份代码、同一条路径、只把容差放粗（`--rel=5e-1 --maxc=50 --maxbox=40`），`sign clamp on I-constants (header (vii)): 6x`，交付表里 `W_kap`、`W_phi`、`GU_kap`、`GU_phi` 等的下端正是被夹到 `0` 的。**一条从未触发的检查与一条不可能触发的检查不能长得一样**——tiny 把「能触发」演出来了，于是 production 的 0x 是「买下来的」而不是「构造上不会发生」。这也复证了 §9b 的机制：粗格上 Peano 余项确实把这些下端压过零。

**(2) 一个硬编码的计数，正是本里程碑在批评的那类记账缺陷。** 摘要行印的是 `of the 52 CV constants: Tier-1 46, Tier-2 0`——46 + 0 ≠ 52。查主文档：权威划分是 **46 个求积常数 + 6 个网格 sup = 52 个对象**（`:259`；`:267(c)` 把早期的「~30 个常数」改正为「46+6」）。所以 52 是**对象**数，被当成 **CV 常数**数硬编码进了 f-string，而 `TIER` 里实际只有 46 项。没有任何检查会抓到它：分子由表算出、分母是字面量，二者从不相互核对。已修为**分母一律由表计数**（`len(M3.TIER)`、`len(M3.TBAR3)`），并把既非 Tier-1 也非 Tier-2 的项显式印成 `OTHER n` 而不是让它消失在总数里；结尾那句 `OK: all 52 CV constants ...` 同样改为计数。**这条一并订正了我自己的 §1.J.0 草稿**，它原本写「52 个求积常数 + 6 个网格 sup」（即 58），现为「46 个求积常数 + 6 个网格 sup（合 52 个对象）」。修改落在 `main()` 的打印段，属 AST 消去比较里的报告阶段，不触及产界区 Sec. 1–7b；因文件尚未重新冻结，此时改的代价为零。


---

## 7b. (R10) 对 `cv_quad` 的追问：它剪枝吗？——不剪，而且不变量比 `tbar_bb` 的**更强**

协调者在定 `--maxc` 之前挡下发车问了这一条，理由是对的：smoke 日志里 φ 的 `cells=556` 撞上了 `--maxc=400`（外层格上限在**用来选 cap 的那次运行里就已经触发了**）。若 `cv_quad` 像 `theta_star_certified` 那样对着一个运行中的 `best` 丢格子，则 cap 退出是**可靠性**风险，9.769080 这个基线连同据此得出的一切结论都要作废。

**答：`cv_quad` 不剪枝。** 函数里没有 `best`、没有任何与运行界的比较、没有任何丢弃。细化循环对每个格子只做两件事之一——原样留下，或换成它的两个半格 `(a,b) → (a,m),(m,b)`（`a < m < b`，两处用的是**同一个** float `m`）：

```python
m = _split_pt(a, b) if i in cut else None
if m is None or not (a < m < b):
    nc.append((a, b))
else:
    nc.append((a, m)); nc.append((m, b))
cells = nc
```

四条退出路径（容差达标 / 格数上限 `max_cells` / 无格被选中 `cut` 为空 / 无格可分 `fresh` 为空）全都在这个替换之后离开，故**每条退出都只影响宽度，不影响可靠性**——与 `tbar_bb` 的结论相同。

**但它不是同一个不变量，而且这个区别是承重的**：`tbar_bb` 对盒子取 **max**，max 容忍重叠，所以它只需要**覆盖**；`cv_quad` 对格子取**和**（`I_acc = Σ (wᵢ·hᵢ + peᵢ)`），并且把定号段沿 `cells[i][0] .. cells[j][1]` 望远镜式地合并，所以它需要的是**分割**：一个缝隙会**静默丢掉质量**——积分、进而所报的界会**偏小**，是 fail-open——一个重叠会重复计数。

> **更强的不变量要配一条检查，而不是一条注释。**

协调者把这一条提为总文档的 **(R13)**，原话：*the aggregation operator decides what invariant you owe — max needs a cover, sum needs a partition.* 它是「剪枝 / 不剪枝」之外的第三根轴：两个骨架看上去一样（分治、细化、从不丢弃），fail-closed 与 fail-open 的差别只在最后一步是 `max` 还是 `Σ`。

所以装配前加了 fail-closed 断言：`cells[0][0] == 0.0`、`cells[-1][1] == B`、逐 i `cells[i][0] < cells[i][1]` 且 `cells[i][1] == cells[i+1][0]`（浮点**精确**相等——半格两侧共用同一个 `m`，所以精确相等是可以要求的，不需要容差；参见 (R11)：加容差就是把类型/结构错误糊过去）。

**这条断言有牙**（`_scratch_ir/m5_partition_chk.py`，通过 patch `outer_edges` 让**真实**的 `cv_quad` 吃到坏分割，而不是复制一份循环来测——见 §4 的空守卫规则）：

| 变异 | 结果 |
|---|---|
| 真实边（基线） | 通过，8 格，`W_G` 正常返回（46 s） |
| 边不从 0 开始 | **AssertionError**（43 s） |
| 退化格 `a == b` | **AssertionError**（45 s） |
| 接缝处有缝隙 | **AssertionError**（49 s） |

第三个变异需要说明：平铺的边表**无法**表达缝隙（格子由 `[(e[i], e[i+1]) ...]` 构造，相邻格共享同一个下标）。所以缝隙用一个 `list` 子类制造——下标 2 第二次被读时返回略大的值，那正好是「格 1 的右端 / 格 2 的左端」这条接缝。这也顺带说明：**在现有代码里缝隙只能由细化循环产生，而细化循环两侧共用同一个 `m`，所以产生不了缝隙**——断言目前恒真。

**恒真的断言守的是下一次编辑，不是这一次。** 这必须写明，因为一条今天不可能触发的断言，在不知道它为何在那里的人眼里，长得和一个空守卫（§4）一模一样，而区别正是全部要点。空守卫的定义是：参照物由**同一条代码路径或同一种算术类型**产出，所以它检测不到任何东西——它的牙在结构上被拔掉了。这条分割断言的参照物（首边 0、末边 B、接缝精确相等）**独立于**被测的细化循环，牙是有的：上表证明只要生产分割的一方坏了（patch 掉 `outer_edges`），它就咬。它今天咬不到，是因为**当前的生产者**恰好正确，不是因为它咬不动——「当前不可被本生产者证伪」与「空」是两回事。这也是它该是 **assert** 而不是 clamp 的理由（§9b 的分工）：clamp 是给「有证明、预期会触发」的场合的，把一个预期永不触发的性质做成 clamp，就是把下一次有人改坏细化循环时本该响的警报静音成一次静默的收窄。

---

## 8. 与 M2 的对照该怎么定级

`check_constants_vs_m2` 和 `check_tbar` 把 M5 的区间与 M2 的 float64 值比。协调者要求的定级，原样记录：

> a float64 grid max on a **different** code path is a real cross-check ... but with 1e-3/5e-3 one-sided slop it is **evidence that the answer is in the right place, not a bound**.

即：它能抓住「算错了一个数量级」「符号反了」「常数名接错了」，**不能**当成包围性证明。真正的包围性来自区间算术本身。

---

## 9. Θ* 的缺陷、中心形，以及它顺带清掉的一处引用

M3 出厂的 `theta_star_iv` 在 dps 40 算 Θ*、转成 float64、再按**二分宽度**（~1e-45）加 pad。这个 pad 覆盖的是错误的误差源，而且**自我掩盖**：根找得越好、pad 越小、包围越差。出厂区间**漏掉真实 Θ\* 达 4.834991e-17**，并且是**向下**漏——而 Θ\* 在 `bounds3` 里是当**上界**进去的，所以漏的方向是**不安全**的那一侧。

M5 的替代是对 Θ 做区间分支限界。这里需要**中心（均值）形** `Θ(box) ⊆ Θ(m) + Θ'(box)·(box − m)`，`Θ'(s) = N(s)/(1+s²)²`，`N(s) = 2s + (1−s²)atan s − s·log(1+s²)`。要记录的是**为什么**需要它：

> **是区间依赖误差、而不是极大值的几何形状，决定了分辨率。**

直接形在极大值附近的宽度由 `s` 在分子分母里重复出现的依赖误差主导，二分收敛得极慢；中心形把一阶依赖消掉。

副产品：中心形**清掉了 M3 对一处未记录的「N 单调下降、单根」证明的引用**——M5 不需要该性质，分支限界对 Θ 的极大值不做任何单调性假设。

数值上这两处原语的修复**不动任何打印位**（M4 基线上 ϖ_A″ 9.711925 → 9.711925，δ = +1.6e-07），但方法声明需要这个修复。

---

## 9b. `W_φ` 的负下端（smoke 运行，协调者指出；(R12)）

smoke 在 `--maxc=400` 下给出 `W_phi = [-0.0229346119588, 5.21463976575]`。问的是：W_φ 能否合法为负；若不能，下游是否对符号敏感。

**答**：(a) **不能为负**——见下面两条非负性证明。负下端纯粹是 Peano 双侧余项在粗格上把下端压过零。(b) **代码不是符号无关的**——`M3.bounds3` 把 `W_f + (D/2)·V_W_f` 送进 `iv.sqrt`，而 mpmath 的 `iv.sqrt` 对负下端**抛 `ComplexResult`**（本机验证）。也就是说：跨零的 W_φ 不会让界变宽，而是让装配**异常退出**。

**预测（从源码，服务器时钟约 02:20，smoke 起于 00:30:24、当时在 TBAR 3/6，运行秒数 ≈ 6400 s）**：按 smoke 的数值，`report_ratio` 在 ϖ=1000 处（D/2=5e-4，3.281·5e-4 = 1.6e-3 < 0.0229）会抛 `ComplexResult`；`find_wA2_wide` 的括号 [0.5, 64]（64 处 3.281/128 = 0.0256 > 0.0229）、payload 的 ϖ=45（0.0365 > 0.0229）与 dump 都在崩之前完成。运行中的进程已 import 旧模块，磁盘上的修改不影响它；按协调者要求**不杀不重启，让它崩**，回溯写进 `smoke.log`。

**观察（服务器时钟 03:37，运行秒数 11247.9 s 之后）**：**没有到达预测的那一行。** 进程在 `dump_constants`（源码第 2536 行，`'TAILS': {k: list(v._mpf_) ...}`）以 `AttributeError: 'ivmpf' object has no attribute '_mpf_'` 崩溃——在全部 18 个对象算完之后、在任何字节写入磁盘之前。原因是我自己的类型错误：`tbar_tail_bound` 返回 `(...).b`，是**退化 `ivmpf`**，而 dump 假定 `mp.mpf`。`--dump` 路径此前**从未被行使过**（run_a 没活到那里，tb 探针不 dump）。所以：(1) ComplexResult 的预测仍未被观察，`report_ratio` 从未被执行；(2) smoke 的 dump **没有幸存**，3.1 小时的 18 个区间只剩日志里的 12/8 位十进制打印。这是 smoke 作为管道测试**本该抓到、也确实抓到**的那类缺陷，只是代价比预想的大。

**关闭方式**：`dump_constants` 改为按 `_mpi_`/`_mpf_` 任一类型取精确端点对；`build_all` 在**每个族、每个 sup 之后**重写 dump 文件（先写 `.tmp` 再 `os.replace`，不留半截文件），晚崩不丢；`load_constants` 拒绝不完整的文件；往返在真实类型上做了单元测试（`_scratch_ir/m5_dump_chk.py`：18 个值逐位往返、部分文件被拒、加载后全为 `ivmpf`）；`build_all → 增量 dump → load` 的整条管道（含新的 `--maxbox` 到 6 个 sup 的传递、`complete` 标志、逐位重载）以打桩积分器跑通（`m5_tinybuild.py`，RESULT OK）。

**但这是同一个病因的第三例，这才是发现**：`_abs_range` 一类函数返回**退化 `ivmpf`** 而不是 `mp.mpf`，同一处歧义已经产出 (1) `fam_sups` 的 `max()` 误用（fail-open，§6.2）、(2) `tbar_bb` 退出测试的 `None`-falsy（松紧，§6.2）、(3) 本次 `dump_constants` 的 `_mpf_`（崩溃）。**三处补丁加一个类型无关的取端点器不是修复，只是把第四次推迟。** 真正的修复是在**产生处**消除歧义：给退化区间一个显式的构造/析构边界——要么让所有「标量」出口一律返回 `mp.mpf`（`_abs_range`、`tbar_tail_bound`、`.a`/`.b` 的所有包装）并在需要区间时显式 `_iv_pt()` 回升，要么反过来一律返回 `ivmpf` 并禁止裸 `.a`/`.b` 外泄，两者都配一个类型断言在边界上把关。**本里程碑没有做这件事**（它要动 M3/M4 的共享工具函数，超出本次范围），记在 §12。

**预测的补测（服务器，`--load=smoke_transcribed.json`，单进程，约 20 分钟）**：把日志里的 18 个区间转录为 dump 格式（`m5_transcribe.py`，`iv.mpf(str)` 向外舍入，**只作 steering 输入、不是精确 dump**）后回放。`--load` 不经过 clamp，W_φ 的下端保持 −0.0229。**结果：预测被确认，位置与操作数都对。** 日志打印到 `w=45.0000 ratio in [0.427985, 0.481477]`、`ok: ratio.hi < 1 at w=45.0`，随即在 `report_ratio` 的下一个 ϖ（列表是 `{wA2, 10, 45, 1000}`，即 **ϖ=1000**）抛出：

```
  File "/zx/p4cert/m5/p4_cii_m5_interval.py", line 2898, in main
    wA2 = report_ratio()
  File "/zx/p4cert/m5/p4_cii_m3_interval.py", line 941, in bounds3
    + q * iv.sqrt(CV3['W_' + nm] + (D / 2.0) * CV3['V_W_' + nm]))
  mpmath.libmp.libmpf.ComplexResult: square root of a negative number
```

即 §9b 开头写的那条表达式（`dZ` 项，nm='phi'，D/2 = 5e-4，V_W_phi 下端 3.281，5e-4·3.281 = 1.6e-3 < 0.0229），不是 `NR_hi` 里的 `W_G` 那条。**第三次独立复现**发生在本机：灵敏度跟进脚本把 `V_GU_phi` 展宽 1e6 倍作负控制时，`find_wA2_wide` 的括号被推到更大的 ϖ，撞上同一行同一异常——说明这不是某个 ϖ 的巧合，而是「W_φ 下端为负 + `iv.sqrt`」这一对组合在 D 足够小时必然发生。

**这同时说明 clamp 不是可选的**：没有 clamp，一个**合法但粗**的包围会让装配在大 ϖ 处**异常退出**，而不是给出一个更宽的界。

**关闭方式与理由（只允许这一种写法）**：在 `cv_quad` 出口对两个 I-常数做 `I ∩ [0,∞)`。**它可靠，是因为真值落在交集里**：真值在区间算术的包围里，也在 [0,∞) 里。**不是**因为「只动下端所以上界不变」——那句话一般不成立（下端经任何减法、取负、除法都能到达上界，这里 W_f 的下端正坐在一个大表达式内部的 sqrt 里），此处碰巧无害，而「碰巧无害」正是它危险的地方。clamp 所依赖的两条假设，任一条错则 clamp 不可靠，已写进模块头诚实边界 (vii)：
  (a) h_W(s) = s³Û(s)²/2 ≥ 0（s ≥ 0）：Û = T_f(s)/s² = TV_σ(ψ) 是全变差（绝对值之和，Sec. 4），平方 ≥ 0；
  (b) h_GU(s) = G(s)·s·Û(s)/2 ≥ 0（s ≥ 0）：**G(s) = 1/(1+s²) > 0 是定义**（g₁，M3 Sec. 1），外层格 s ≥ 0，Û ≥ 0 同 (a)。
于是 W_f = ∫h_W ≥ 0、GU_f = ∫h_GU ≥ 0 精确成立。

**clamp 与 assert 的分工**：有证明且预期可能触发 ⇒ clamp（两个 I-常数，计数 `SIGN_CLAMP`，`main()` 一条 chk「12 个 pad 常数 ≥ 0 且 clamp 0 次」，`--load` 下只能测存储端点、标签如实注明）；有证明且预期永不触发 ⇒ **fail-closed `assert`**（V-常数：定号 run 逐段 clamp 到 0、歧义格 `[inf|h′|, sup|h′|]`、尾 `[0,t]`；TBAR 下端：`lo_best` 从 0 起单调增；`cv_quad` 的分割不变量，§7b——当前生产者结构上不可能违反它，它守的是下一次编辑）。一条永远不会失败的检查与一条从未失败过的检查在日志里**不能长得一样**。

**clamp 的第二个理由，与可靠性论证并列、不替代它**：§9b 前文的 ComplexResult 三次复现说明，没有 clamp，一条**合法但粗**的路径（Peano 余项在粗格上把下端压过零，包围本身正确）在 D 足够小时会让装配**抛异常**而不是**返回一个更宽的界**。clamp 的可靠性来自「真值在交集里」；clamp 的**必要性**来自「它让一个合法但宽的包围以合法但宽的形式呈现，而不是以失败的形式呈现」。两条理由回答的是两个不同的问题（能不能做 / 为什么非做不可），都要写，写在一起。

**披露**：`load_constants` 原样回放存储端点，clamp 在 `cv_quad` 里、`--load` 不经过它，所以灵敏度扫描的基线 W_φ 下端仍是 dump 里的 −0.0229；扫描脚本打印这一点，并**直接测量** clamp 前后的 ϖ_A″ 差（见 §11），而不是论证它。

---

## 10. 受信任基座（不是空的，如实列出）

替换掉 `atan` 之后，装配路径上仍在调用：

| 原语 | 状态 | 证据 |
|---|---|---|
| `iv.sqrt` | SOUND | 0 / 40000 未包围 |
| `iv.log` | SOUND | 0 / 120000 未包围（含 M4 用的 `log(1+s²)` 形状） |
| `iv.pi` | SOUND | 端点与 `libmp.mpf_pi(prec, floor/ceiling)` 逐位相同 |
| `iv.atan2` | **UNSOUND，已替换** | 2502/40000 随机 + 1429/20000 系统操作数未包围，最差 1.41e-32（~0.2 ulp） |
| `mp.sqrt`（M3 `V_FPP_iv` 内） | 带 pad 的 float 路径 | **27 次调用**，dps 40 计算坐在一个显式 **1e-12 绝对 pad** 之下，比自身舍入误差高 28 个数量级 |

最后一行是本次唯一一处「pad 顶着 float 路径」的残留。按 §4 的 pad 规则检查过：pad 覆盖的正是 dps-40 的舍入误差源。**这次的缺陷是 sizing error，不是 design error**——pad 的位置对，只是当初选大小时没算过。

`check_transcendental_unreachable` **报告**调用计数而不是断言它们，这样一个新出现的原语会显式冒出来，而不是悄悄加入基座。

---

## 10b. 18 个对象的逐个灵敏度，以及其中 5 个不承重

**为什么先量再发车**：协调者的规矩是「measure first, then launch」——`--maxc` 该定多少，不靠论证，靠测。做法是在 smoke 的 18 个区间上，每次只动一个对象，重算 ϖ_A″：**(W)** 把上端再加一个自身宽度（「cap 再松一档最多亏多少」的上界），**(C)** 塌缩到下端（「完美包围最多赚多少」）。

**基线 ϖ_A″ = 9.769080**（smoke 宽度：`--maxc=400, rel=2e-1, top=200, s0=64`）。

| # | 对象 | 相对宽 | d(W) | d(C) |
|---|---|---|---|---|
| 1 | V_GU_phi | 0.899 | +0.033887 | −0.034320 |
| 2 | V_W_phi | 0.985 | +0.011755 | −0.017167 |
| 3 | V_W_kap | 0.926 | +0.010565 | −0.013245 |
| 4 | **TBAR[kap]** | 0.819 | +0.008823 | −0.008849 |
| 5 | W_phi | 1.000 | +0.006129 | −0.007214 |
| 6 | W_kap | 0.415 | +0.005798 | −0.006499 |
| 7 | GU_phi | 0.087 | +0.004042 | −0.004045 |
| 8 | **TBAR[phi]** | 0.908 | +0.003678 | −0.003684 |
| 9 | **TBAR[G]** | 0.134 | +0.002651 | −0.002653 |
| 10 | V_GU_kap | 0.776 | +0.000683 | −0.000684 |
| 11 | V_W_G | 0.185 | +0.000526 | −0.000533 |
| 12 | W_G | 0.018 | +0.000361 | −0.000364 |
| 13 | GU_kap | 0.032 | +0.000079 | −0.000079 |
| 14–18 | GU_G, V_GU_G, TBAR[UG], TBAR[Ukap], TBAR[Uphi] | — | **+0.000000** | **−0.000000** |

（`m5_sens.py` / `m5_sens.log`。d(W) 之和 +0.088977，d(C) 之和 −0.099335。）

**哪个 cap 在卡**：6 个 sup 合计只占总余量的 17%，12 个 CV 常数占 83%；真正撞上 900 盒上限的两个 key（kap、phi）合计 +0.0125——即**把它们做到完美**也只赚 0.012，而在这两族上每盒 1.5–3 s。结论：**`max_boxes` 保持 900，`--maxc` 才是该调的旋钮**。同时把 `max_boxes` 提成 CLI 开关 `--maxbox=N`（穿到全部 6 个 sup，日志头打印实际值），理由是这么承重的数不该躺在函数默认值里没人看得见。

**联合塌缩**（`m5_sens2.py`）：18 个全部塌缩到下端得 **9.653067**（d = −0.116013）。注意它**大于**逐个 d(C) 之和（−0.099335）——装配是非线性的，逐项之和不是联合下限，不能当界用。所以 M5 的全部宽度在当前 M3/M4 基线上只值 **0.116**。

**那 5 个 0 是结构性的，不是舍入**。读源码：`bounds3` 只在 908 行读 `TBAR3['G']`、948 行读 `TBAR3['phi'|'kap']`；`CV3['GU_'+nm]`/`['V_GU_'+nm]` 的唯一消费者是 `padX_iv`，而 `bounds3` 只对 nm ∈ {phi, kap} 调它。于是 `TBAR3['UG'|'Ukap'|'Uphi']` 与 `GU_G`、`V_GU_G` **一次都没被读**。M2 的 `bounds` 同理。三个数值佐证：

- 把这 5 个的宽度**乘 1e6**，ϖ_A″ = 9.769080，**d = +0.000000**；
- 只塌缩 13 个「活」对象、5 个保持原样，得 **9.653067**，与 18 个全塌缩**完全相同**；
- 负控制：`m5_sens.py` 表里把 V_GU_phi 只展宽一倍就动了 +0.033887（`m5_sens2.py` 里乘 1e6 的那次没给出数——界坏到把 `find_wA2_wide` 的括号推出可行域，撞上 §9b 那条 `iv.sqrt`；作为负控制它仍然有效：动「活」对象会把计算改到面目全非，动这 5 个**一个比特都不动**）。

**「代码里不承重」不等于「本该不承重」**。协调者的定级必须原样保留：

> A structurally unread object is one of two things: a harmless leftover from an earlier assembly, or **a term the derivation says should be there and the code dropped**. Those look identical from inside `bounds3`, identical in the widen-by-1e6 stress test, and identical in the dump — and only one of them is fine. The second is fail-open and is exactly the shape of defect this whole exercise exists to catch.

所以在 §1.J 写「18 个里 13 个承重」之前，必须**正面**从 `P4_einstein_extremal.md` 的推导确认：被装配的那条不等式**本来就没有**这几项。

**核查结论：(a)，是遗留物，不是被丢掉的项。** 两条论据都是推导里的原话，不是「观察到它没被读」。

**(a1) G 族没有 X pad，是因为那个二次型被按符号整个丢掉了。** `⟨τ,𝓡_G τ⟩` 正是「其格点→连续 pad 会是 q·I[G·U_G] = q·GU_G」的那个对象。推导（`P4_einstein_extremal.md:475–477`）把它丢了，并且**明说**丢完之后不需要 pad：

> **h ≥ ρ = 2B − 1 + (c_Δ + ⟨τ,𝓡τ⟩)/B ≥ 2B − 1 ≥ …**
> 两个被丢项都 ≥ 0（𝒟 ⪯ I；𝓡 逐元 ≥0 且 τ ≥ 0），**无任何离散化 pad**。

方向是对的：这是 h 的**下**界，丢掉一个 ≥ 0 的项仍是合法下界。紧接着的 `:497–499` 说明 φ/κ **不能**照此办理，这正是 `padX`（从而 `GU_φ`、`GU_κ`）只对这两族存在的原因：

> （注意 𝓡_φ 与 𝓡_κ **不**逐元 ≥0，min ≈ −0.32，故须保留其二次型常数，不能丢）

即：**pad 只在核不逐元非负时才需要**；𝓡_G 逐元非负，于是丢弃取代了 pad。

**(a2) 装配里的两个 sup 都是对「不加权」行和取的，所以要的是 T̄ 不是 Ū。** G-加权的内和误差在推导里始终是**逐点**消费、再用**积分**聚合，从不取 sup。`p4_cii_m2.py:126–128`（装配的记录式推导，SEC. 4）：

> `sqrt(N_R) = sqrt( I[Ups_G^2] + q V2[Ups_G^2] ) + q sqrt( I[U_G^2] + (D/2)V1[U_G^2] )`
> *(Minkowski in the w-weighted l^2: the inner-sum error of Ups_G^lat(n) is **<= q U_G(s_n) pointwise**, by (L1) in sigma)*

这是 `W_G = I[U_G²]` 的来源，也是 `U_G` 唯一合法的用处——逐点。而两处取 sup 的地方 `:131` `g >= 2B - 1 - Theta* - q Tbar_G` 与 `:142` `||O[f]|| <= sup_s Phi_f(s) + q Tbar_f + D^2/3` 用的都是**无 G 权**的 T̄（行和权只有 w_m = D/(2s_m)）。§1.G.3 的 pad 目录（`:603–619`）逐行按**用处**列出 11 个 pad，**没有一行含 `I[G·U_G]`，也没有一行含 U 的 sup**；其仅有的两行 sup（`:607`、`:617`）都是 T̄。M2 对抗核验里 reviewer 独立重写装配、逐位对上（`2026-09-12_cii_m2_review.md:13`），其装配路径清单也只有 `qT̄_G`、`padX`、`qT̄_f`。

**遗留物的来路**：`sup_pad_profile()` 在 `p4_cii_m2.py:604–606` 按 3 族 × 2 权的**统一循环**调用，产出 6 个 sup 而装配只消费 3 个；`build_constants`（`:491–495`）同样按三族统一循环建 `W_`/`GU_`，于是产出了 `GU_G`。**46+6 的清单是按「脚本算了什么」枚举的，不是按「装配读了什么」枚举的**——这就是 5 个从不被读的对象被一路认证、检查、加余量带到今天的原因。

**因此这不是 fail-open，而是一处账目/文档缺陷**：§1.G.4(3)（`:630`）与 `p4_cii_m2.py:116–117` 称 U-sup 是「为界」严格化的，而界从不用它们；§1.H.2（`:691`）、§1.H.6（`:723`）、§1.I.6（`:811`）、§1.I.7（`:829`）的 52 计数把承重集**多算了 5 个**。已交回协调者（§1.H 是 M3 的章节，本里程碑不改）。

**顺带**：`P4_einstein_extremal.md` §1.H（第 829 行）称这 6 个 sup 是「当前 ϖ_A″=9.71 而非逼近 M2 的 8.48 的主因」。上面的测量与之矛盾——其中 3 个贡献恒零，6 个加起来还不到 V_GU_phi + V_W_phi。**本里程碑不改 §1.H**（那是 M3 的章节），已交回协调者处理。

## 10c. 五个对象的逐个推导侧归类（专项纯阅读，2026-09-13）

§10b 的 (a1)/(a2) 是一次合并论证；协调者要求**逐对象**落到三种结局之一，各给一句可引用的论证，并明确这是「唯一一件不能由代码侧解决的」事——灵敏度为 0 对「无害残留」与「推导要求有、代码漏掉」读数完全相同。于是派了一个只读推导、不碰代码、不跑东西的 reviewer 子 agent，任务书里把 §10b 的结论明确标为**待检验的假设**。它的结论如下，其中每一处引文我（implementer）都回到原文核对过。

三种结局：(a) 推导里根本没有此项，须点名生成模式；(b) 推导有、在具名步骤按具名理由合法丢弃；(c) 推导有、不可丢弃、代码未读 = fail-open。

**先厘清一个命名陷阱**（写附录时会踩）：`phi_GU`/`kap_GU`（M2 `build_constants` :498）= I[G·Υ_f]，是**主**常数 X_f = ⟨τ,𝓡_fτ⟩ 的连续部分（M4 起 Tier 1）；`GU_phi`/`GU_kap`/`GU_G`（:495）= I[G·U_f]，是 **pad** 伴随量，U_f(s) = ∫|∂_σ²(G(σ)R_f(s,σ)/2σ)|dσ 是 G-加权内和的 (L1) σ-误差剖面。本节五个对象全是第二种及其 sup。

**论证走正面可数的事实，不走「找遍了没找到」。** 缺席是最弱的一类证据，也正是 fail-open 缺陷最擅长伪装成的那类。§1.G.3 的 pad 目录（`:607–617`）按**用处**逐行列出全部 **11 行** pad，读者扫一眼就能自己数：

- **T̄ 出现 2 行、共 3 个实例**：`:607` `q·T̄_G`（g_lo 里 sup_n 的 𝓡 行和）与 `:617` `q·T̄_f + D²/3`（一行覆盖两族，f = κ 与 φ₀ 的 Schur 范数）。**没有一行带 T̄_{U_f}**——三个 `TBAR[U*]` 的 (a) 由此直接可数。
- **I[G·U_f] 出现 2 行**：`:611`（φ）与 `:613`（κ，「同上（f=κ）」）。G 族对应的那行是 `:610`，用的是 `q·√(I[U_G²]+(D/2)V₁[U_G²])`——**ℓ²/Minkowski 路线，不是 G-加权**。**没有 `I[G·U_G]` 行**——`GU_G`/`V_GU_G` 的 (b) 由此直接可数。

配套的正面事实是**每个 U_f 实际在哪里被消费**（不是「没在哪里」）：格点内和只有两种。(i) **无权**行和 Σ_m w_m R_f(s_n,s_m) = Θ_f^lat(s_n)，误差 ≤ q·T_f(s_n)（`:576–578`），只在 **p-加权 Schur 检验**里出现、即只在 sup_n 之下（λ₂ ≤ ‖𝓡_G‖ ≤ sup_n Σ_m R w_m，`:486–488`；‖𝒪[f]‖ ≤ sup_s Φ_f，`:505`）——要的恰是 sup_n T_f ≤ T̄_f，对应 `bounds3` :908/:948 读 `TBAR3['G'|'phi'|'kap']`。(ii) **G-加权**内和 Σ_m w_m G(s_m) R_f(s_n,s_m) = Υ_f^lat(s_n)，误差 ≤ q·U_f(s_n) **逐 n**，只从 𝓡_fτ（τ_m = p_m G(s_m)）产生，恰在三处被消费，三处都由 **w-加权求和**聚合：C_r 里 ‖𝓡_Gτ‖ ≤ √S[Υ_G²] + q·√S[U_G²]（Minkowski，`p4_cii_m2.py:126–128`）；X_f = ⟨τ,𝓡_fτ⟩ ∈ S[GΥ_f] ∓ q·S[G·U_f]，f = φ,κ（`:133`）；e_f 里 ‖ΔZ‖ ≤ √ε2 + q·√S[U_f²]（`:140`）。再由 (L2) 转积分得 W_f + (D/2)V_W_f、GU_f + (D/2)V_GU_f。**三处都是求和，一处 sup 也没有。**

**本段最强的一句：对两处 U-站点，sup 路线不是「未选」，而是「不可用」。** Σ_n w_n = Σ 1/(2n+1) **发散**，故若用 Ū_f·√(Σ_n w_n) 去界 √(Σ_n w_n U_f(s_n)²)，得到的是 **+∞**——一个无用的界。三处 G-加权站点中只有 X_f 一处原则上可用 Ū_f（那里的权是 w_n G(s_n)，Σ w_n G(s_n) = A < ∞），而推导在那里也用了积分（`:133`）。所以「推导取 T̄ 而不取 Ū」不是风格选择，是两处**别无选择**、第三处**另有选择而未取**。

| # | 对象 | 结局 | 可引用论证 | 位置 |
|---|---|---|---|---|
| 1 | `GU_G` = I[G·U_G] | **(b)** | G-加权内和误差剖面 U_G 只进入 ⟨τ,𝓡_Gτ⟩ 的界，而该二次型只出现在下界 h ≥ ρ = 2B−1 + (c_Δ+⟨τ,𝓡_Gτ⟩)/B 中，并在那里作为非负项被丢掉（𝓡_G 由 (M1) 逐元 ≥ 0，τ > 0）；⟨τ,𝓡_Gτ⟩ ≥ 0 在格点上**精确**成立，故不需要任何格点→连续 pad，也就不需要 I[G·U_G]。 | `:446`（R ≥ 0）、`:450`（𝓡_nm ≥ 0）、`:473–477`（丢弃，原话「两个被丢项都 ≥ 0 … 无任何离散化 pad」）；对照 `:497–499`（𝓡_φ、𝓡_κ 不 ≥ 0，「不能丢」，这就是 padX 只对 φ/κ 存在的原因）。**丢弃的合法性另需 §10c.1 的方向性论证。** 来路：`build_constants` 对 `('G','phi','kap')` 循环（`p4_cii_m2.py:491–495`）每族同产 `W_` 与 `GU_`，`W_G` 是 C_r 要的，`GU_G` 是副产品。 |
| 2 | `V_GU_G` = V₁[G·U_G] | **(b)** | 它是 `GU_G` 在 (L2) 下的伴随项，只以 q(GU_f + (D/2)V_GU_f) 的形式出现（`p4_cii_m2.py:133, 632`），GU_f 不需要则它也不需要。 | 同上；同一生成循环（`:495`）。 |
| 3 | `TBAR[UG]` = 1.5·sup_s U_G | **(a)** | G-加权剖面 U_G 逐点消费、由 w-加权 ℓ² 范数（Minkowski）聚合，得 √(I[U_G²] + (D/2)V₁[U_G²]) = √(W_G + …)；对 𝓡_G 行和取的唯一 sup 是 p-加权 Schur 检验，其行和**无权**、pad 是 T̄_G，故 sup_s U_G 从不进入任何界。 | `:481–483`、`:486–488`、pad 目录 `:607`（g_lo 用 q·T̄_G）与 `:610`（C_r 内和用 q·√(I[U_G²]+…)）；`p4_cii_m2.py:126–128`。生成模式：`build_sups`（`p4_cii_m2.py:604–606`）对 3 族 × 2 tag 一律调 `sup_pad_profile(nm, tag)`——U-tag 的**剖面**确实需要（逐点，供 W_f/GU_f），其 **sup** 不需要，循环对两个 tag 一律取了 sup。 |
| 4 | `TBAR[Ukap]` = 1.5·sup_s U_κ | **(a)** | U_κ 进入 ϖ⁴E 的界有两处——\|q_κ\| 里 ⟨τ,𝓡_κτ⟩ 的双侧 pad，与 ‖ΔZ‖ = ‖𝓡_κτ − p·Υ_κ‖——两处逐点误差 q·U_κ(s_n) 分别以 w_n G(s_n) 加权求和、以 w_n 加权平方求和，得 I[G·U_κ] 与 I[U_κ²]；κ 族出现的 sup 是 ‖𝒦‖ 的 Schur 界里**无权**剖面 T_κ 的 sup，故 sup_s U_κ 不是任何界的项。 | `:499–503`（Z、‖Z‖²、⟨τ,Z⟩）、`:505`（Schur 的 sup_s，无权行和）、pad 目录 `:613, :616, :617`；`p4_cii_m2.py:133, 140, 142`；装配 `p4_cii_m3_interval.py:919–924, 938–943, 946–948`。同一生成模式。 |
| 5 | `TBAR[Uphi]` = 1.5·sup_s U_φ | **(a)** | 与 #4 同构，f = φ₀：两处是 ϖ⁴F(t) 的下界（`:497`，pad 行 `:611`）与 e_φ 里的 ‖ΔZ‖（`:616`）；φ 族的 sup 是 ‖𝒪[φ₀]‖ 里的 T̄_φ（`:617`）。 | 同上。 |

**对「从不需要 sup」做的对抗检查**（子 agent 在原文核过，我复核了引文）：有没有以 τ 为权的 Schur 检验（那会产生 G-加权行和、从而需要 Ū_f）？没有——三处 Schur（`:486`、`:505`、`p4_cii_m2.py:557`）都是 p-加权。h（或 ρ）是否在任何地方作为**上界**被需要（那会禁止丢弃 ⟨τ,𝓡_Gτ⟩）？没有——判据 x²E²/g < F 里 g 只在分母，DK 用 ρ−λ₂ ≥ g_lo（`:508–514`）。无权剖面是否在任何地方被积分消费（那需要一个不存在的 I[T_f²] 型常数）？没有——Θ_f 只在 sup_s Φ_f 之内。§1.H/§1.I 是否把这五个之一算作被消费？没有——它们按「脚本算了什么」列（§1.H.2 `:691`、§1.I.7 `:829`），不按 `bounds3` 读了什么。

**唯一一处推断而非阅读**：「统一循环」的来路是从 M2 代码结构（`:371–373`、`:491–495`、`:604–606`）推断的，推导里没有一句「顺带也算了这些」。**这条是关于这五个常数「为何存在」的推测，不承重**：论文需要的命题只是「推导从不消费它们」，那一条已由上面两组可数事实（pad 目录里没有 T̄_U 行、没有 I[G·U_G] 行）与 §10c.1 的方向性论证证成，与它们**为何被算出来**无关。放在这里只为解释现象，写进 §1.J 时须同样标明其等级，不得让它挂在已证成的论证末尾拉低整段等级。

### 10c.1 丢弃的方向性：这条论证唯一可能塌的地方

（协调者补的一步，我已回原文核对。）丢掉一个非负项只有在 ρ **自始至终只作下界**时才合法；若 ρ 在任何地方被当作上界用，丢掉 ⟨τ,𝓡_Gτ⟩ ≥ 0 就会把界推向不安全的方向，(b) 随之塌掉。原文两行把方向坐实：

> `:514` **判据：E_hi²/(g_lo·F_lo) < 1，g_lo = 2B−1−Θ\*，s = C_r/g_lo。**
> `:508` sinθ := ‖(I−P)t̂‖ ≤ ‖r‖/(ρ−λ₂) ≤ **C_r/g_lo =: s**

g_lo 变小的后果有三条，方向一致：(1) 它在判据 `:514` 里**只出现在分母**，g_lo ↓ ⟹ E²/(g_lo·F_lo) ↑；(2) 经 `:508` 有 s = C_r/g_lo ↑，而 `:510` 的 F ≥ (1−s²)·ϖ⁴F(t) − 2s·e_φ − s²‖𝒪[φ₀]‖ 对 s **单调下降** ⟹ F_lo ↓；(3) 同样经 s，`:510` 的 E ≤ e_κ + s(|q_κ| + ‖𝒦‖) 对 s **单调上升** ⟹ E_hi ↑。三条都让判据更难满足，即**丢弃是 fail-closed 的**。

补充的两条封口（代码侧，方向性之外）：ρ 只经 g = ρ − λ₂ 进装配（`p4_cii_m2.py:641–642`、`p4_cii_m3_interval.py:907–908`）；残差 r（`:479`）不需要 ρ 的值；E、F、q_κ 用 𝒪[κ]、𝒪[φ₀] 而非 𝓡_G——没有第二处能让 ⟨τ,𝓡_Gτ⟩ 以禁止丢弃的符号再进来。**这一步成立，(b) 才成立。**

**判决：五个全部落在 (a)/(b)，没有 (c)。**「18 个对象里只有 13 个承重」可以写。二分法干净且由推导自身结构给出，并且是**可数的**而非「找遍了没找到」：pad 目录 11 行里，T̄ 占 2 行/3 个实例（G、κ、φ₀）、I[G·U] 占 2 行（φ、κ）、G 族走 ℓ² 路线占 1 行（`:610`），没有任何一行是 T̄_U 或 I[G·U_G]。无权内和只在 sup_n（p-加权 Schur）下消费、取 T̄_f；G-加权内和只在 w-加权求和下消费、取 I[U_f²]、I[G·U_f]，且其中两处 sup 路线因 Σ w_n 发散而**不可用**；G 族唯一的二次型 ⟨τ,𝓡_Gτ⟩ 在 `:475–477` 按符号从一个下界里精确丢掉，方向性由 §10c.1 的三条 fail-closed 论证支撑。§10b 把 #1/#2 标为 (a) 是在更粗的分类下；按本节的三分法它们是 (b)——**项在推导里，pad 不在**，丢弃发生在 pad 产生之前。9.769080 **不需要**重评。

两条文档后果（与 §10b 一致）：§1.G.4(3)（`:630`）与 `p4_cii_m2.py:116–117` 称 U-sup 是「为界」严格化的，而没有界用它们；§1.I.6(2)/§1.I.7（`:812`、`:829`）称六个 T̄ sup 是 ϖ_A″ = 9.71 对 8.48 的主因，对三个 U-sup 不可能成立——它们按推导贡献恒为零。

---

## 10d. 在**交付宽度**上重量一次灵敏度（2026-09-13）

§10b 那张表是在 smoke 区间上量的（基线 9.769080）。它做出的决定（`--maxc` 提到 1500）已经执行完了，于是它的输入不再存在：生产区间与 smoke 区间不仅宽度不同，**下端也不同**（`W_phi` 下端 −0.0229 → +2.1395）。继续引用 smoke 的 −0.099335 / −0.116013 / 9.653067 来描述交付品，就是「拿一次检查去测它不会面对的东西」——与已作废的 tiny 退出记分、与 s0 标签同类。所以重量一次。

做法：`_scratch_ir/srv_m5_sens_prod.py`，载入生产 dump `m5_final.json`，走**最终冻结脚本 `3fdcfc7e`** 的装配路径（不是工作副本），每次只动一个对象，(W) 上端再加一个自身宽度、(C) 塌缩到下端。基线与两次还原都精确复现 **9.716960**（`match=True`），所以那些 `+0.000000` 不是还原坏掉的假象。dump 里**没有**负下端（`stored lower endpoints < 0: none`）——这与 §9b 的 smoke 情形相反，生产不需要 clamp。

| 对象 | 相对宽 | d(W) | d(C) |
|---|---|---|---|
| V_GU_phi | 7.83e-01 | +0.014608 | **−0.014689** |
| **T̄[κ]** | 8.19e-01 | +0.008987 | **−0.009014** |
| V_W_phi | 9.19e-01 | +0.005072 | −0.006425 |
| V_W_kap | 7.89e-01 | +0.004756 | −0.005342 |
| **T̄[φ]** | 9.08e-01 | +0.003752 | −0.003758 |
| W_kap | 1.32e-01 | +0.002048 | −0.002148 |
| W_phi | 2.96e-01 | +0.001831 | −0.001977 |
| GU_phi | 1.65e-02 | +0.000751 | −0.000751 |
| **T̄[G]** | 8.52e-03 | +0.000150 | −0.000150 |
| V_GU_kap | 2.67e-01 | +0.000136 | −0.000136 |
| W_G | 5.51e-03 | +0.000111 | −0.000111 |
| V_W_G | 3.62e-02 | +0.000095 | −0.000095 |
| GU_kap | 6.28e-03 | +0.000016 | −0.000016 |
| GU_G, V_GU_G, T̄[U_G], T̄[U_κ], T̄[U_φ] | — | **+0.000000** | **−0.000000** |

**三条结论。**

1. **五个惰性对象在交付宽度上仍然精确为零。** §10c 的推导侧归类因此不只在 smoke 上成立；这是那条论证在**实际发运的数**上的一次独立确认（仍然只是佐证——灵敏度为 0 分不清 (a)/(b)/(c)，(c) 的排除靠的是 §10c 的可数事实，不是这张表）。
2. **排序变了，而且是被 §10b 自己的决定改变的。** `T̄[κ]` 从第 4 升到第 2：它的 d(C) 从 −0.008849 到 −0.009014 几乎没动，而其它项普遍减半——因为 `--maxc` 加了四倍而它是 **box-cap-set** 的，那个旋钮碰不到它。**下一档该调的是 `--maxbox`，不再是 `--maxc`。** 一次在 smoke 上量出的旋钮选择只对它量的那一轮有效。
3. **(R14) 在生产宽度上仍然成立，方向也一样**：逐项 d(C) 之和 **−0.044611**，联合塌缩 **−0.046032**（ϖ_A″ = **9.670928**）——逐项相加**低估**联合效应，这一列是排序不是预算。

**最有用的一个数：M5 的全部剩余宽度只值 0.046。** 把 18 个对象做到完美包围，ϖ_A″ 也只从 9.716960 到 **9.670928**。对 M2 的 float64 真值 8.4761，差距 1.241——**M5 的宽度占 3.7%**。由区间算术的包含单调性（把任一常数换成自身的退化子区间只能收窄输出），任何**子集**的联合效应都被 0.046 控制，**六个 T̄ sup 也在内**。所以 §1.I.6(2)/§1.I.7 的「六个 T̄ sup 是 ϖ_A″=9.71 而非逼近 8.48 的主因」在数量上就不成立：三个 U-sup 按 §10c 恒为零，另外三个连同全部 12 个 pad 常数一起最多值 3.7%。**剩下的 96% 在 M3/M4 自身的保守性与判据装配里**，那是下一个里程碑该看的地方。

---

## 11. 最终数字

**头条：严格 ϖ_A″ = 9.716960。** 46 个求积常数 + 6 个网格 sup（合 52 个对象）全部 Tier 1，无任何 float64 余量记账。与段 C-I 的重叠 **4.631×**（认证到 45）、**12.350×**（认证到 120）。

沿革：M2 float64 8.4761 → M3 9.831910 → M4 9.711925 → **M5 9.716960**。**M5 比 M4 退了 0.005035**：M4 里这 18 个对象仍带 float64 自收敛余量，余量比真正的向外舍入包围窄，所以 M4 的数偏乐观。退这一步是把乐观换成界。**这四个数不是同类物，不能当收敛曲线读**：8.4761 是 float64 真值、根本不是界；9.831910 / 9.711925 是界，但分别有 36 / 18 个对象的内容来自余量记账；9.716960 是第一个每一项都由包围支撑的数。

两条命令、两份日志、两套记分：

| | 构建 82cca55f | 报告 3fdcfc7e |
|---|---|---|
| checks | 68 passed / 0 failed / **0 not exercised** | 66 / 0 / **3 not exercised** |
| 墙钟 | 16317.2 s（32 进程） | 1515.3 s（`--load`） |
| 剖面求值 | 92 013 914 jets / 23 777 896 σ-格 | 0 / 0 |
| sign clamp | **0x**（挣来的真零） | not exercised（`--load`） |

3 条 not exercised 恰好是 `--load` 阶段**不能证伪**的三条（解析 cap 不相交、`_isect` 不相交、clamp 计数），按 (R15) 印 not exercised 而不是印默认值结果。两份日志检查数不同不是不一致，是分档记账在起作用。

判据表：

```
w=    9.7170   ratio in [0.764429, 1.000000]   width=0.235571   ok at w=9.716959756753118
w=   10.0000   ratio in [0.748197, 0.971564]   width=0.223367
w=   45.0000   ratio in [0.428009, 0.481449]   width=0.053440
w= 1000.0000   ratio in [0.320613, 0.339820]   width=0.019206
```

两个 Θ\* 的 25 位并排（`--load` 报告原样）：

```
Theta*  M3 as shipped = [0.6349287917881079046367176, 0.6349287917881079046367176]  width 0.00e+00  <- does NOT enclose
Theta*  M5 certified  = [0.6349287917881079529866324, 0.6349287917881079529866324]  width 4.50e-29
wA'' with M3's iv.atan2 + float64 Theta* = 9.716960
wA'' with the certified arctan + Theta*  = 9.716960   (delta +0.000e+00)
```

δ 为 0 是**数值上的巧合，不是可靠性论证**：M3 出厂的 Θ\* 区间宽度为零且漏掉真值 4.834991e-17、方向向下，而它在 `bounds3` 里当上界。替换的理由是这条，不是判据动没动。

18 个对象的交付表（包围为 dump-derived；9 条退出中 6 条 log-printed、3 条 log-derived）：

| 对象 | 包围 | 相对宽 | 要求 | exit.src | 退出 |
|---|---|---|---|---|---|
| W_G | [1.151331608, 1.157714903] | 5.51e-03 | 0.05 | derived | tolerance met |
| V_W_G | [2.833053003, 2.939526939] | 3.62e-02 | 0.05 | derived | tolerance met |
| GU_G | [0.7688478014, 0.7698585162] | 1.31e-03 | 0.05 | derived | tolerance met |
| V_GU_G | [1.915850586, 1.921235073] | 2.80e-03 | 0.05 | derived | tolerance met |
| W_kap | [1.553803067, 1.791000512] | 1.32e-01 | 0.05 | derived | cell cap 1500 |
| V_W_kap | [2.953197326, 13.98702232] | 7.89e-01 | 0.05 | derived | cell cap 1500 |
| GU_kap | [1.054020017, 1.060679567] | 6.28e-03 | 0.05 | derived | cell cap 1500 |
| V_GU_kap | [3.077318622, 4.199943452] | 2.67e-01 | 0.05 | derived | cell cap 1500 |
| W_phi | [2.139498617, 3.037942395] | 2.96e-01 | 0.05 | derived | cell cap 1500 |
| V_W_phi | [4.541005954, 55.97724536] | 9.19e-01 | 0.05 | derived | cell cap 1500 |
| GU_phi | [1.377281751, 1.400457852] | 1.65e-02 | 0.05 | derived | cell cap 1500 |
| V_GU_phi | [2.431596078, 11.21679419] | 7.83e-01 | 0.05 | derived | cell cap 1500 |
| TBAR[G] | [1.431161363, 1.443465245] | 8.52e-03 | 0.01 | printed | converged (442 盒) |
| TBAR[UG] | [1.409631847, 1.418835459] | 6.49e-03 | 0.01 | printed | converged (329 盒) |
| TBAR[kap] | [0.4535260849, 2.507085948] | 8.19e-01 | 0.01 | printed | **box cap 900** |
| TBAR[Ukap] | [1.794956157, 1.807715813] | 7.06e-03 | 0.01 | printed | converged (710 盒) |
| TBAR[phi] | [0.4723879218, 5.138182464] | 9.08e-01 | 0.01 | printed | **box cap 900** |
| TBAR[Uphi] | [2.214284233, 2.231277434] | 7.62e-03 | 0.01 | printed | converged |

**9 条退出里 4 条 cap-set、5 条容差达标**（cap-set：cv:kap、cv:phi、tb:kap、tb:phi）。

**§10b 的 measure-first 预测对了一半，另一半的细账见 §10d。** 它说「`max_boxes` 保持 900、该调的是 `--maxc`」；生产把 `--maxc` 从 400 提到 1500（3.75×），ϖ_A″ 从 smoke 的 9.769080 到 9.716960（赚 0.0521）。**旋钮选对了，但「调完哪条宽度会掉」全反了**：按 smoke 灵敏度排前四的四个对象（`V_GU_phi`、`V_W_phi`、`V_W_kap`、`T̄[κ]`）恰好是几乎没动的四个（−13% / −7% / −15% / **0%**），而改善 66–94% 的全是排在后面的。κ/φ 两族 8 个 CV 常数**仍骑在 `cell cap 1500` 上**。**把这两件事记成同一次成功是不诚实的。**

**正式同一性核销（§7d.2 / §1.J.8 第 2 条）**：`82cca55f` 与最终冻结的 `3fdcfc7e` 在同一条 tiny 命令行下各跑一次（两次日志头逐项相同），**24 个条目逐位相同、0 处不同**。中途对中间 md5 `019bc17d` 的那次 24/24 只作记录、不构成核销。

**退出重建核验**：`m5_exitlog_chk.py` 四组 **35 项全过**（前置 4 / A 10 / B 13 / C 8）。SUITE C 现在跑的是**完整**生产日志，tb 侧是 4 converged / 2 box cap 的**真混合**，所以「box cap」那条分支是在有活的对立选项下被判别的。

---

## 12. 还没关上的

1. **C10**：源码守卫不覆盖 M3/M4，只覆盖 M5 自己的 Sec. 1–7b。已在守卫输出与模块头如实打印。
2. **可达性守卫是动态的**（§3(b) 原话）。
3. ~~**`tbar_bb` 的 cap 退出路径在生产容差下未被行使**~~ —— **这条已被生产日志证伪，留痕不删**。生产（s0=4096, tbrel=1e-2, max_boxes=900）下 `TBAR[kap]` 走的正是盒上限退出（`box cap 900 after round 6`），给出 `[0.45352608, 2.50708595]`（smoke 预测 `[0.453, 2.508]`）；同一次运行里 `TBAR[G]`、`TBAR[UG]` 在 442 / 329 盒收敛退出。所以 tb 侧是**真混合**，不是一律骑 cap。按 (R10)，cap 退出是**松紧**风险而非可靠性风险，因为 `tbar_bb` 不剪枝、从不丢盒；行使与否逐 key 记在交付表的退出列里。
4. **`M3.V_FPP_iv` 的 27 次 `mp.sqrt`** 仍在 1e-12 pad 之下（§10）。
5. **M3 Tier-1 内部 8 个 V-配对仍是 200% 相对宽度**（M3 §1.H.3 已证不承重，本里程碑未动）。
6. **退化 `ivmpf` 与 `mp.mpf` 的类型歧义未在产生处消除**（§9b）。已产出三个缺陷（fail-open 一个、松紧一个、崩溃一个），现在每处都单独修好了，但病因还在；要在产生处关掉需要改 M3/M4 共享的工具函数，超出本里程碑范围。
7. **5 个对象在装配里不承重**（§10b 量出、§10c 从推导侧定性）：`GU_G`、`V_GU_G`、`TBAR[UG]`、`TBAR[Ukap]`、`TBAR[Uphi]` 被 `bounds3` 完全不读；推导侧按 §1.G.3 pad 目录的**可数事实**（11 行里 T̄ 占 2 行 3 实例、无一行带 T̄_{U_f}；`I[G·U_f]` 占 2 行、G 族走 ℓ²/Minkowski 而无 `I[G·U_G]` 行）与 Σw_n 发散使 sup 路线**不可用**，定为 (a)/(a)/(a)/(b)/(b)，无一为 (c)。它们仍被认证、检查、记账——**这本身不是缺陷，但 46+6 的计数把承重集多算了 5 个**，对应的文字在 M3/M2 的章节里，已交回协调者。
