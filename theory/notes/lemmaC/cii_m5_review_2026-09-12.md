# M5 对抗核验记录（segment C-II 最后 18 个对象的端到端区间化）

- 被审对象：`theory/derivations/p4_cii_m5_interval.py`
- 审查形式：reviewer 子 agent（opus）一轮对抗核验 + 本 implementer 的变异测试、区间类型清扫与守卫复验
- 分支：`worktree-agent-a8687e419d54121f8`（**未** merge 进 `p4-conj2`）
- 日期：2026-09-12/13

本文件记录**被挑战的东西和它们怎么被关掉的**，不是脚本说明（方法见 `P4_einstein_extremal.md` §1.J 与模块头注释）。

---

## 0. 一句话结论

reviewer 提 5 项 MUST-FIX、11 项 CHALLENGED。**5 项必修全部关闭并逐条用运行验证**；C1–C9、C11 关闭，C10 改为在守卫自身输出里如实打印。此外本会话自查又发现两个**同类**缺陷（(R11)，见 §6），其中一个是 **fail-open**，两个都不是任何既有检查发现的、只能由审计发现。

---

## 1. 五项 MUST-FIX

| # | 问题 | 为什么严重 | 关闭方式 |
|---|---|---|---|
| M1 | `chk("PI_IV encloses pi", _contains0(PI_IV - 4*iv_atan_nn(iv.mpf(1))))` | **空检查**：`PI_IV` 本身就由同一次调用定义，恒真 | 换成对**精确有理 Machin 括号**的包含检查（`_pi_bracket`，纯整数/`Fraction`），另加 §4 的绝对锚点 |
| M2 | 交错级数余项 pad 无任何检查覆盖 | 删 pad 或改单侧，**没有一个检查会失败**（N=20 截断误差~2⁻¹⁶⁴，prec 103 下亚 ulp，既有检查看不见） | 新增 (A3)：`nterms=1..5` 共 25 个用例，用**精确有理**部分和 `[S_n,S_{n+1}]` 验证 pad 双侧且支配截断 |
| M3 | 毒化守卫被 `from...import` 绕过 | 若 M3 内部用 `from mpmath import log` 取本地引用，替换 `M3.mlog` 不生效，守卫**静默变弱** | BAN/WATCH 改为 `(obj, ctx, nm)` 三元组，显式加入 `(M3,'M3','mlog'/'digamma'/'polygamma')` |
| M4 | BAN 条目缺失被静默跳过 | 名字拼错/改名 ⇒ 原语未被毒化，守卫仍报「全部不可达」 | 新增 `missing` 列表，缺失即失败 |
| M5 | `check_constants_vs_m2` 只比宽度不比包含 | 常数可整体偏移而检查照过 | 改为非对称包含 `lo<=ref*(1+1e-6) and hi>=ref*(1-2e-3)`，宽度另设 sanity floor `width<1.0`（见 §8 定级） |

## 2. 十一项 CHALLENGED

| # | 内容 | 处置 |
|---|---|---|
| C1 | σ-尾图表无检查 | 新增 `check_tail_chart()`：两条精确图表恒等式各 108 点区间相交判定；最差相对间隙 4.54e-30 / 9.62e-28 |
| C2 | `1.0/B` float 倒数可能向下舍入 | 改向外舍入 `(iv.mpf(1)/iv.mpf(B)).b` |
| C3 | Peano 中点法则未计入舍入后中点与精确中点之差 | `_cell_view` 显式包围位移 `E=m'-(a+b)/2`，余项加 `w·|E|·sup|h'| + (ΔQ/2)·E²` |
| C4 | `iv_atan_nn` 对退化输入断言是脚坑 | 去断言，改端点分裂递归 |
| C5 | `theta_star_certified` 的尾 majorant 需 `s0>=2` 未强制 | 加 `assert s0>=2.0` |
| C6 | `_isect` 发现不相交时只返回较窄的一个 | 加 `ISECT_HIT` 计数、`main()` 里 chk 为 0 |
| C7 | 未执行检查与通过检查混在一起计数 | 加 `skip(name,why)`，摘要分打印 passed/failed/not exercised |
| C8 | 源码守卫行区间可能空/错 | 加「非空且合理」chk，teeth 探针移到区间内部 |
| C9 | `find_wA2_wide` 判据用区间比较 | 改为 `_lohi(...)[1]>1`（精确 `mp.mpf`） |
| C10 | 源码守卫只覆盖 M5 自己 Sec.1–7b，不覆盖 M3/M4 | **未修**，如实记入守卫输出与模块头（见 §7） |
| C11 | 多处用 `ivmpf` 做比较 | 见 §6：比 reviewer 提出时更严重 |

---

## 3. 两个守卫，及必须原样保留的限度说明

**(a) 源码守卫**（`check_no_transcendental_source`）——`tokenize` 扫 M5 自己 Sec.1–7b（1329 行），断言无 `mp.`/`iv.` 超越函数调用；纯文本缺席声明，**不覆盖 M3/M4**（C10）。teeth 探针（`iv.make_mpf@1130`）证明扫描范围确实覆盖该区间。

**(b) 可达性守卫**（`check_transcendental_unreachable`）——毒化 16 个原语后重算判据。限度必须原样陈述，不得软化：

> **it demonstrates unreachability along the path exercised with the inputs used, not a static proof for all inputs.**

不要改写成 "verified unreachable"。行使的路径含产出最终 ϖ_A″ 的那次 `bounds3()` 调用，仍是**动态**证据。

两个无法毒化的例外，守卫自己打印：`iv.pi`（常量非函数）、`M3.Z3H_IV`（import 时已由 `mp.zeta` 算好并套 1e-13 pad）。

---

## 4. 空守卫规则（M1 的一般化）

M1 之所以是空的，是因为**参照物和被测物出自同一条代码路径**：

> **a guard whose reference is produced by the same code path, or the same arithmetic type, as the thing under test is not a guard — it is a consistency check wearing a guard's name.**

同一条规则解释了 M3 原来的 `theta_star_iv` 为什么自欺（§9）：它用**二分宽度**当 pad，真正误差源是 dps-40→float64 转换。**pad 必须覆盖它所替代的那个误差源。**

---

## 5. (R9)：arctan 验证集的同质性盲点

**发现**：M5 arctan 门原靠三条恒等式（余角/倍角/Machin）+ 一次 700 位参照。三条恒等式对 atan 都一次齐次，`PI_IV` 自己也由 `iv_atan_nn` 构造 ⇒ 整个恒等式集在 `atan → C·atan` 下不变，**结构性看不见**回乘因子 `2**k` 的错误（`iv_atan_nn` 末行正是 `return s * (2**k)`）。

**变异实验**（`PI_IV` 须从变异体重算，只替换 `iv_atan_nn` 会让 Machin 检查因错误理由抓到变异体；协调者定级：*Recomputing PI_IV from the mutant was necessary, not thorough*）：

修复前（`m5_atanmut.py`）：`scaleback_x2`/`scaleback_half`（`2**(k±1)`）**没有任何恒等式抓到**，只靠两条 mp.atan 参照；`n20_norem` 正确通过（误差~2⁻¹⁶⁴亚ulp）。

修复后（`_m5_mut2.py`）：baseline 0 失败（12 通过，承重证据）；`scaleback_x2`/`_half` 各 6 失败（含绝对锚点+负控制+精确有理 Machin）；`rem_deleted`/`rem_onesided` 各 1 失败（新 pad 双侧性检查，修复前无人能抓）；`den_2n2` 8、`sqrt_1pu` 7；`n8` 仅宽度检查失败（可接受，仍合法只是更宽）。

**配置修正**（协调者撤回先前「恒等式严格优于参照测试」的说法）：

> Appealing to no external value is what makes a homogeneous identity set unable to see `atan → C·atan`. The right configuration is not either/or: **identities plus a rational absolute anchor.**

一般形式：**first ask what group your validation set is invariant under, then ask what it can detect.** 且 **the anchor set has to span the parameter the defect lives in.**

**锚点构造**：对 0<u<1，交错级数相邻部分和括号 `u−u³/3 < atan u < u−u³/3+u⁵/5`，两端有理数，不引用任何超越函数/库原语/高精度参照；比较用 `_cmp_ratio` 纯整数运算。锚点取 dyadic 3/4,1/2,1/4,1/8（严格 > `_ATAN_RED=1/16`，钉住减半循环+级数+回乘三者；用 dyadic 而非 2-b 建议的 9/10 以保证二进制精确可表示）。三重验证（`m5_anchor.py`）：700 位真值包含（宽度 1.18e-01→6.10e-06）、baseline 通过、两个回乘变异体 8/8 全部触发拒绝，耗时 0.9ms。

**顺带**：arctan 输出宽度 N=12 起饱和于 1.8e-30（prec 103 ulp 量级），N=20 是余量非必需。

---

## 6. (R11)：偏序当全序用——本会话自查发现的两个缺陷

**规则**：

> 区间携带的是**偏序**；`max`、`min`、`sorted`、`bisect` 以及 `>` 的布尔用法都假定**全序**。`ivmpf > mpf` 在两者重叠时返回 **`None`**，`None` 是 falsy，于是 `if a > b` 静默走 else 分支，`max(...)` 可能返回一个**非最大元**。

**两个实例**：(1) `fam_sups`（**fail-open**）：`max(max(val[ab] for ab in boxes), tl[k])`，`val[ab]` 是非退化 `ivmpf`；返回值当**上界**用 ⇒ 非最大元 ⇒ **上界被低估**。(2) `tbar_bb`（仅松紧）：`top <= lo_best*(1+rel)`，右边向外舍入下可能宽出一个 ulp ⇒ 退出测试变 `None`-falsy；因该循环从不剪枝（见 §7）只影响松紧不影响可靠性。

**修法**：两处都**恢复全序**——循环内部一律用 `mp.make_mpf(z._mpi_[0/1])` 取精确 `mp.mpf` 端点，**不是加容差**（加容差会把类型错误糊过去）。`fam_sups` 出口再转回退化区间（`_iv_pt`），因 `tbar_tail_bound` 里的 `M[0]*g[1]/2` 须保持**向外**舍入而非就近舍入。

**检测方式**：两个都不是任何检查发现的，**必然如此**——误差方向让界看起来更好。审计是句法的：走遍产界路径每一个 `max(`/`min(`/`sorted(`/`sort(`/裸比较，逐个问「参数重叠时返回什么」。

**清扫结果**：AST 枚举全文件 **181 个排序/比较站点**（`_scratch_ir/m5_r11_sweep.py`），逐个按操作数类型定级；**2 个缺陷已修，第三处 `_isect` 的宽度比较**（两个退化区间相减仍可能向外舍入宽出一个 ulp）也改为精确 `mp.mpf`。**清扫 181 站点，2 个缺陷，剩余 0。** 顺带订正 `_abs_range` docstring（它声称返回 exact mpf 实际返回退化 `ivmpf`——正是这个说明引向 `fam_sups` 误用），并把跨零分支的 `mp.mpf(0)` 统一成 `iv.mpf(0)`。

---

## 7. (R10)：`tbar_bb` 的剪枝对比与作用域

协调者要求记录：

> **`theta_star_certified` prunes, so its soundness genuinely depends on `best` being a verified lower bound. Same branch-and-bound skeleton; whether it prunes decides whether an untested exit path is a tightness risk or a soundness risk.**

`tbar_bb`/`fam_sups` 都 `boxes=keep+new`，**从不丢弃盒子**，故覆盖不变量恒成立、`max hi` 无论哪条退出都是合法上界；`theta_star_certified` 真的剪枝（`[ab for ab in boxes if val[ab]>best]`）。

必须原样保留：

> *the convergence exit is exercised on 2 of 6 keys at loose settings, both in the easiest (G) family; the cap/stagnation exit — the path that actually changed — is not.*

探针实测（`srv_m5_tb.py`/`tb_probe.log`，4 进程，`rel=3e-2,rel_prof=3e-3,s0=512`）：`TBAR[G]` 291盒/565次/523.8s，`TBAR[UG]` 208盒/399次/477.1s，均由收敛测试退出，OK 0 bad of 2。

**smoke 后记**（协调者指出上句需更新）：cap 退出路径**已被行使 2 次**——`TBAR[kap]`/`TBAR[phi]`（`s0=64,tbrel=2e-1,tbrelprof=2e-2,--procs=20`）都以 `boxes=900` 退出，得 `[0.45343049,2.50756953]`（1790次/1476.9s）与 `[0.47233251,5.13818246]`（1790次/2660.5s）——无用地宽但合法，与覆盖不变量论证一致。限度：smoke 容差下的 cap 触发**不能**说明生产容差下 cap 是否触发。

**盒预算不是时间预算**：上面两个 key 花掉**同样**900盒/1790次求值，墙钟却差 1.4 倍——单盒成本随族（G≪κ<φ）和盒内难度变化，`max_boxes` 加大的墙钟代价不能按盒数线性外推，必须按最贵族算。两个旋钮：`--maxc`（CV 外层格子上限，只影响 12 个 pad 常数）与 `tbar_bb` 的 `max_boxes=900`（新增 `--maxbox`，生产保持 900）——两者无关。

**生产容差一次反转**（协调者发车 12 分钟后叫停重发）：`tbrelprof` 由 1e-3 改回 **1e-2**。理由：预期 6 个 key 都会骑 900 盒 cap，**骑 cap 的宽度由 cap 决定不由容差决定**——更紧内层容差只是花 20 倍成本换更窄的逐盒区间，再被 cap 整个扔掉。协调者原话：*Running a tight inner tolerance into a binding outer cap is the worst of both: you pay for the precision and the cap discards it.* 定量对比：`tbrelprof=1e-3` 外推 TBAR 半边单独要 15–30h（压倒一切）；`1e-2` 时约 2.4h，低于 CV 半边~2.8h（φ 族在 `--maxc=1500` 下主导，~1.7h）。被叫停运行（`m5_final_tbrelprof1e-3_aborted.log`）尚在单进程准备段（534s，M4 常数刚算完），无损失。重发命令：`--procs=32 --rel=5e-2 --tbrel=1e-2 --tbrelprof=1e-2 --maxc=1500 --dump=m5_final.json`。`--procs` 未提到 40 的理由：机器不受内存/磁盘限制，留在 32 唯一理由是用户原始约束「44核别开满留一半」，32 已越过它不再往上；两条用户指令的出入提交用户裁定，裁定前不动。

---

## 7c. 交付表第三列：退出原因在退出处记录，不从计数倒推（(R10)/(R13) 落到交付物）

κ、φ 两族预期在 `--maxc=1500` 骑 cap，**宽度是 cap 定的，不是容差定的**——拿一个没人选过的数比一个我们指定的数不对。所以逐对象表要三列：**实际宽度、要求的容差、哪条退出触发**；退出原因要**在触发处记录**，不许事后从格数/盒数倒推。

实现（`p4_cii_m5_interval.py`，md5 `da8dfb5e`；`533fa042` 为前一版少一段注释，仅本机牙用）：`cv_quad` 的 `why` 在触发分支赋值（`tolerance met`/`cell cap N`/`no cell cut`/`no cell split`/`round cap N`，容差与 cap 同时成立写 `tolerance met (cells=N also >= cap N)`）；即时打印+存入 `QuadInfo.exit`/`EXITS['cv:<族>']`；`tbar_bb` 同样即时打印三条退出（`converged`/`box cap N`/`round cap 40`）并写入 `EXITS['tb:<key>']`；`dump_constants`/`load_constants` 把 `EXITS`/`SETTINGS`（实际 rel/tb_rel/两个 cap）与端点一起写盘/读回，缺失时说明；新的 `report_objects()` 打印 18 行表（包围/相对宽度/要求容差/退出/`cap-set: k of n`），放在 `report_ratio()` **之前**——表是构建产物，不依赖装配。牙（`_scratch_ir/m5_exit_chk.py`，走真实函数）全过。

**分支保留的解析证明**（(R4) 类陈述，写进 `cv_quad` 的 `no cell cut` 分支）：容差未达时 `no cell cut` **不可达**——`sl` 是 c≥1 个正 pw 之和，`n_i=max(8,c)>=c`，最大 pw ≥ sl/n_i > rel·E/n_i，而阈值是 0.3·rel·E/n_i，1/0.3=3.3 裕量吸收浮点求和误差。证明关于**当前** `thr`，改 `thr` 可能使其失效；实际可达退出是容差/cap/`no cell split`/轮数上限四种。

**对正在跑的生产进程的裁定**（决定因素不是 40 分钟，而是脚本出处）：它由 md5 `82cca55f` import，dump 里没有 `EXITS`。选项 A（杀掉重发）vs B（跑完但脚本 md5 与产出脚本分裂）交协调者，**裁定 A**——理由：**B 交付的脚本不是产出 artifact 的脚本**，与 layer-two-b 曾需字节级 diff 事后补救、留永久 caveat 是同一种出处缺陷；「每个认证常数出自一条命令」的里程碑要求下明知故犯站不住。

**换脚本前的门**（协调者，零墙钟成本）：新脚本先端到端跑最小设置（`--procs=6 --rel=5e-1 --top=20 --maxc=50 --s0=8 --tbrel=5e-1 --tbrelprof=5e-1 --dump=tiny.json`），再 `--load=tiny.json` 一次，须走到 `dump_constants→report_objects→report_ratio` 末尾。理由：`--dump` 此前从未被行使、首次使用即失败在末尾（§9b，代价 3.1h）；`report_objects`/`EXITS`/`SETTINGS` 同样从未被真实 `build_all` 走到末尾行使过。**一条唯一职责是保存昂贵工作的路径，必须在昂贵工作存在之前被行使。**

---

## 7d. 门的结果、(A) 的撤销、一个误报自身参数的检查

**门：绿。** `tiny_e2e.log`（pid 2936527，md5 `da8dfb5e`，`--maxbox=40`，2882.0s）走完全程：dump→检查电池→毒化守卫→**18行对象表，末行`cap-set: 9 of 9`**→`rigorous wA''=1520.136655`→ratio 表（ϖ=1000 未抛 `ComplexResult`，sign clamp 报 6 次而非静默）→SUMMARY。`--load=tiny.json` 单进程同样走到末尾。

**17 个 FAIL 逐条核对，无一是内部不变量**：7×`relative width<1 (sanity floor)`、6×`Sec.7 majorant...does not bind`、1×`sign clamp applied 6x`、3×`ratio.hi<1 at w=10/45/1000`——全是在故意荒谬宽度下不可能成立的**精度**断言（W_phi=[0,1.76e15]⇒ratio=[−inf,+inf]）。

**(A) 被协调者撤销，理由是价格变了 3 倍。** 我报价 (A) 时只有 G 族算完；到重估时生产日志（G 1164s、κ 4487s 已落盘，φ 在第5轮）已 ~8100s，重发要扔掉~4500s 且**无法续跑**（`load_constants` 拒绝部分 dump）。真实价格约 2h 而非 40 分钟。协调者原话：**「一个以 40 分钟为输入的决定，在价格翻三倍时要重做。」**（可复用规则：批准是对某个价格的批准，价格变了批准不再自动有效。）

**替代方案（现行）**：生产在 `82cca55f` 下跑完；交付表由**第二条命令**`--load=m5_final.json --exitlog=m5_final.log` 在冻结脚本下产出。代价如实写在 §1.J：**artifact 的 md5 与交付脚本的 md5 不同**，做两次出处核销（下）。

**一个真实缺陷，只有 tiny 跑才使它可见**：6 行 `Sec.7 majorant for s>4096 does not bind` 打印的是模块常量 `TBAR_S0`，该次运行用 `--s0=8`。**检查的消息里硬编码了一个它并不读取的参数**——生产设置下 s0 恰好就是 4096，此缺陷将**永远不可观测**。协调者定性：**一个误报自身参数的检查是贴错标签的警报，不是无关紧要的字符串。** 修法：标签改读 `SETTINGS['s0']`，dump 无该键印 `?(s0 not recorded in dump)` 而不代默认值。同一规则修第二处：`--load` 下 clamp 计数印 `NOT EXERCISED here`（不是「0x=从未触发」）。**(R12) 的反面：一个从未触发的检查与一个不可能触发的检查，在日志里必须长得不一样。**

### 7d.1 日志导出的退出：规则必须是全的，证据等级必须标出来

要求：不许写成 `cells>maxc⇒cap`，要对 `cv_quad` 四条与 `tbar_bb` 三条各写判别签名，证明互斥穷尽；两条退出共享签名标为**未决**，不许猜；列标 **log-derived** 区别于 dump-recorded。

实现为代码（`--exitlog=FILE`→`exits_from_log()`）：**cv_quad** 取该族 summary 行前最后一条 `round n: cells=N I-slack...V-slack...`；`%.1e` 打印值 s 框在 `s·(1±0.051)`；确定达标/未达标之间→**未决**。**tbar_bb** 只有一行 `TBAR[k]=[lo,hi]...boxes=B evals=E`，converged ⟺ `hi<=lo·(1+rel)`（1e-6内→未决）；另两处拒绝：(a) 打印 hi 是 `max(bb top, Sec.7 tail)`，尾界即 hi 时非收敛判决未决；(b) `box cap` 与 `round cap 40` **共享签名**（切分被截断可能带 B==max_boxes 离开），报为一桶同时点名。

**签名穷尽性**：`cv_quad` 循环只有五个出口且每轮恰打印一行，三个判别式互斥；`tbar_bb` 同理三个出口两个签名。**牙**：`compare_exits()` 按「矛盾」而非「相等」计分（未决是正确输出）。对 `tiny.json` 实测（剥掉冻结脚本才有的 `cv exit:`/`bb exit:` 行）：**3/9 精确、4/9 同类更粗、2/9 拒绝、0 矛盾**（两处拒绝是 `--s0=8` 特例）；另 15 个变异用例全过（`_scratch_ir/m5_exitlog_chk.py`）。

### 7d.2 出处核销两次，不是一次

**(i) AST 层面消去比较**（比字节级 `SEC.` diff 更强）：解析两文件 Sec.1–7b 内每个顶层 def（各36个），消去 docstring/log/print/`why`赋值/`EXITS`/`SETTINGS`赋值/`QuadInfo.exit`形参等**不可能影响取值**的语句形态，比较 `ast.dump`。结果：**36 个 def 逐个相同，无差异**。负控制：改 `cv_quad` 阈值常数 0.3→0.5 立刻被抓（`_scratch_ir/m5_provenance_ast.py`）。

**(ii) 正向同一性测试**（diff 不是它）：完全相同命令行在冻结脚本与 `82cca55f` 下各跑一次 tiny，比较 dump 端点元组是否**逐位相同**。协调者理由：**diff 是对源码的一次阅读，逐位相同是来自算术的证据。**

**冻结**：协调者要求先停手冻结记 md5 再跑（对中间 md5 的同一性测试是「看起来像证据的陈旧证据，比没有更糟」）。第一次冻结 md5 `019bc17db5d21a5e918687d5b5cb2e4a`。artifact 产出脚本 md5 `82cca55f`。

**这次冻结随后被打破了，如实记**：§7e.3 的日志格式发现与 §7e.1 的 (R15) 拆分都要改进交付脚本，md5 因此变动。处理：(1) 机械 diff 出 `019bc17d`→当前文件的差异集——**36 个 def 逐个相同**，全部差异只在 `exits_from_log`/`compare_exits`/`report_objects`/`main` 四个报告/CLI 阶段 def + 2 个新辅助函数，`build_all` 不在其中；(2) 待修改落定后**重新冻结**并对**最终文件**重跑同一性测试。**核销以 (2) 为准**。中间同一性结果（记录非核销）：`019bc17d` vs `82cca55f`，24 个端点对逐位相同、0 处差异（`_scratch_ir/m5_identity_cmp.py`）。

### 7d.3 (2) 已完成：对最终冻结文件的正式核销

**最终冻结 md5 `3fdcfc7e30b4f8a5a0a97147c62d6986`**（本机与服务器 `p4_cii_m5_final.py` 相同；这是 **LF 行尾**字节的 md5，也是 `git cat-file blob` 值；Windows CRLF 工作树展开后是 `eb8945747d570cf97b99ddd9659b1617`，跑生产/同一性测试的都是 LF 那份）。对最终文件重跑 tiny 命令行，与 `tiny_82.json` 比较：

```
A = tiny_82.json     complete=True  EXITS=0  SETTINGS=0
B = tiny_final.json  complete=True  EXITS=9  SETTINGS=8
24 entries compared, 0 differ
RESULT: IDENTICAL -- every endpoint tuple is bit-for-bit equal
```

**这一次才是核销**；中间那次 24/24 仍只作记录——「反正结果会一样」正是让陈旧证据活下来的那句话。AST 侧同步重跑：Sec.1–7b **36 个 def 逐个相同、CLEAN**；`--all` 差异仅报告阶段 def，**`build_all` 不在其中**。

**两份 tiny 日志检查区并排，给了三条正面证据**：两脚本**失败同一组 17 条**（差异只在报告层的旁证）。三处差异恰是宣称修掉的三处，且可观测：(1) s₀标签：`82cca55f` 印 `s>4096`，冻结脚本印 `s>8`；(2) clamp 检查合取拆分：`82cca55f` 一行混判两个含义，冻结脚本拆两条；(3) 检查总数 68→69，增量正是这次拆分。

---

## 7e. (R15) 与 (R16)：拆分阶段之后欠的两笔账

### 7e.1 (R15)：每条检查要申明哪个阶段能证伪它

协调者原话：*when you split a pipeline into build and report, every check must declare which stage can falsify it, and a check evaluated in a stage that cannot falsify it must print "not exercised" — never its default-valued result.*

符号 clamp 检查是合取式：「端点≥0」两阶段都能证伪；「clamp 从未触发」**只有构建阶段能证伪**（`--load` 下 `SIGN_CLAMP[0]` 停在初值 0，因不可能失败而通过）。先前「修复」只改了说明文字却仍折进 pass/fail 判定——**解释挂在旁边，判决还是错的**。已拆成两条：端点照常 `chk`，「clamp 从未触发」在 `--load` 下走 `skip`。协调者：*a check that has never fired and a check that cannot fire must not look the same in the log.*

另记：§8「一键复跑」条款既是推向拆分的力也是隐患源头，该条款禁的是手抄端点不是第二个进程读机器可读 artifact。**拆分不是选择的，是被迫的**——`report_objects` 在 82cca55f 里根本不存在，两个 md5 从一开始就不可能合成一个。

### 7e.2 (R16)：消去论证的**前提**本身要审

负控制只证明「消去不是把什么都删光」，不证明被删名字对返回的界是**只写**的。读审计（`_scratch_ir/m5_read_audit.py`）**干净**：frozen 在 Sec.1–7b 内对 `SETTINGS` 的 LOAD 次数为 0；真正取值源是 `TBAR_S0`，两脚本三个承载取值的读取点逐行相同；`why` 的 4 次读取全流向 `EXITS`/`log`/`QuadInfo`。

**读审计暴露了比较范围本身的一个洞**：`build_all` 计算着两个尾函数消费的 s0，却位于 Sec.1–7b **之外**——从未进入过比较。放宽到全部顶层 def 后：`build_all` **消去后完全相同**；全部差异仅报告/序列化/CLI 阶段的 6 个 def。一般形式：**一个只覆盖「取值在哪里算出」的比较，仍然必须覆盖「谁算出了这段范围的输入」。**

### 7e.3 解析器被拿去测了它不会面对的东西

§7d.1 的 3/9-4/9-2/9-0 记分是**对 tiny 日志**量的，要解析的却是生产日志（bb verbosity 不同）——记分作废重做，发现两件事：

1. **82cca55f 打印 `bb exit:` 全部三个分支，只是没有 `cv exit:`。** 9 个退出是两个证据档次：6 个 sup 是 **log-printed**，3 个 pad 族是 **log-derived**（唯一可能未决的一档）。先前把好的一半也标错。
2. **refusal (b) 不成立，理由也是错的。** 有了逐轮数据：cap 退出发生在**其自身打印的 B 已达 max_boxes** 的那一轮；轮上限退出则最后打印行 B 仍小于 max_boxes、之后才分裂越界——汇总行报的是分裂后盒数，这才是它分不开而轮次行分得开的原因。协调者：**一个正确的拒绝若靠错误的理由支撑，比不拒绝更糟，因为理由才是读者会拿去复用的那部分。**

牙改成三套共 **35 项全过**（`_scratch_ir/m5_exitlog_chk.py`）。SUITE C 输出（生产日志）：`cv:G [derived] tolerance met`、`cv:kap/phi [derived] cell cap 1500`、`tb:G/UG [printed] converged`、`tb:kap [printed] box cap 900`。**决定性**：`kap w=0 bb exit: box cap 900 after round 6` 是重测期间落盘的真实 artifact，旧规则会误判为模糊桶。

**一个自伤的 bug**：第一版蒙眼按「行含 `exit:` 整行删掉」，而冻结脚本在 `TBAR[...]` 汇总行**末尾**追加 `exit: <why>`——蒙眼删掉了承载数字的那一行，套件报「规则拒绝了」，真相是**测试毁掉了自己的输入**。改为只截后缀，并断言行数前后相等。

### 7e.4 一条横幅，不是三处脚注

`m5_final.json` 只有 `CV5`/`TB5`/`TAILS`/`complete`，无 `SETTINGS`/`EXITS`——s0 标签、9个退出、clamp 计数**三样都是日志来的**，合并为一条横幅而非三处脚注。

### 7e.5 tiny 对拍的两个副产品

**(1) 符号 clamp 在 tiny 宽度下触发 6 次——production 的 0 因此是有意义的零。** tiny 给出**正面**证据（同代码同路径只放粗容差）：`sign clamp on I-constants: 6x`，W_kap/W_phi/GU_kap/GU_phi 下端被夹到 0。**一条从未触发的检查与一条不可能触发的检查不能长得一样**——tiny 演出「能触发」，production 的 0x 是「买下来的」。

**(2) 一个硬编码的计数，正是本里程碑在批评的那类记账缺陷。** 摘要行印 `of the 52 CV constants: Tier-1 46, Tier-2 0`——46+0≠52。权威划分是**46 个求积常数 + 6 个网格 sup = 52 个对象**，52 是对象数被当成 CV 常数数硬编码进 f-string。已修为分母一律由表计数（`len(M3.TIER)`等），非 Tier-1/2 项显式印 `OTHER n`。**这条一并订正了 §1.J.0 草稿**（原「52 个求积常数+6 个网格 sup」即 58，现为「46+6 合 52 个对象」）。

---

## 7b. (R10) 对 `cv_quad` 的追问：它剪枝吗？——不剪，而且不变量比 `tbar_bb` 更强

协调者定 `--maxc` 前挡下发车问了这条：smoke 日志 φ 的 `cells=556` 撞上 `--maxc=400`（cap 在选 cap 的那次运行里就已触发）——若 `cv_quad` 像 `theta_star_certified` 那样对运行中 `best` 丢格子，cap 退出是**可靠性**风险，基线连同一切结论要作废。

**答：`cv_quad` 不剪枝。** 函数里没有 `best`、没有与运行界的比较、没有丢弃；细化循环只做原样留下或替换成两半格（`a<m<b`，两处用同一个 float `m`）。四条退出路径全在替换之后离开，故每条退出只影响宽度不影响可靠性。

**但不是同一个不变量，区别承重**：`tbar_bb` 对盒子取 **max**（容忍重叠，只需**覆盖**）；`cv_quad` 对格子取**和**（`I_acc=Σ(wᵢ·hᵢ+peᵢ)`）+望远镜合并，需要**分割**——一个缝隙会**静默丢掉质量**（fail-open），重叠则重复计数。

> **更强的不变量要配一条检查，而不是一条注释。** ——协调者提为总文档 **(R13)**：*the aggregation operator decides what invariant you owe — max needs a cover, sum needs a partition.*

装配前加了 fail-closed 断言（首边0/末边B/逐i精确相等），**有牙**（`_scratch_ir/m5_partition_chk.py`，patch `outer_edges` 让真实函数吃坏分割）：

| 变异 | 结果 |
|---|---|
| 真实边（基线） | 通过，8格，46s |
| 边不从0开始 | **AssertionError**（43s） |
| 退化格 a==b | **AssertionError**（45s） |
| 接缝处有缝隙 | **AssertionError**（49s，需 list 子类构造缝隙） |

**恒真的断言守的是下一次编辑，不是这一次**——须写明区分它与空守卫（§4）：空守卫参照物由同代码路径/同算术类型产出，结构上拔掉了牙；这条断言参照物独立于被测循环，牙是有的（上表证明坏分割会咬），今天咬不到只因当前生产者恰好正确。该性质做 **assert** 而非 clamp——clamp 是给「有证明、预期会触发」的场合，把预期永不触发的性质做成 clamp 就会把下次警报静音成静默收窄。

---

## 8. 与 M2 的对照该怎么定级

`check_constants_vs_m2`/`check_tbar` 把 M5 区间与 M2 float64 值比。协调者要求的定级原样保留：

> a float64 grid max on a **different** code path is a real cross-check ... but with 1e-3/5e-3 one-sided slop it is **evidence that the answer is in the right place, not a bound**.

即：能抓「数量级错」「符号反」「常数名接错」，**不能**当包围性证明。真正包围性来自区间算术本身。

---

## 9. Θ* 的缺陷、中心形，及顺带清掉的一处引用

M3 出厂的 `theta_star_iv` 在 dps 40 算 Θ*、转 float64、按**二分宽度**（~1e-45）加 pad——覆盖错误误差源，且**自我掩盖**（根找得越好、pad越小、包围越差）。出厂区间**漏掉真实 Θ\* 达 4.834991e-17**，方向**向下**，而 Θ\* 在 `bounds3` 里当**上界**——漏的方向不安全。

M5 替代：对 Θ 做区间分支限界，需**中心（均值）形** `Θ(box)⊆Θ(m)+Θ'(box)·(box−m)`。为什么需要：

> **是区间依赖误差、而不是极大值的几何形状，决定了分辨率。**

副产品：中心形清掉了 M3 对一处未记录的「N 单调下降、单根」证明的引用——M5 分支限界对 Θ 极大值不需要该性质。数值上修复不动任何打印位（M4 基线上 ϖ_A″ 9.711925→9.711925，δ=+1.6e-07），但方法声明需要它。

---

## 9b. `W_φ` 的负下端（smoke 运行；(R12)）

smoke（`--maxc=400`）给出 `W_phi=[-0.0229346119588, 5.21463976575]`。

**答**：(a) **不能为负**——见 §9b 末非负性证明；负下端是 Peano 双侧余项在粗格上把下端压过零。(b) **代码不是符号无关的**——`M3.bounds3` 把 `W_f+(D/2)V_W_f` 送进 `iv.sqrt`，mpmath 对负下端**抛 `ComplexResult`**（实测），即跨零 W_φ 不会让界变宽而是让装配**异常退出**。

**预测→观察确认**：源码预测 `report_ratio` 会在 ϖ=1000 处抛 `ComplexResult`；实际进程在 `dump_constants` 因**类型错误**（`tbar_tail_bound` 返回退化 `ivmpf` 而非 `mp.mpf`）先崩溃于 `AttributeError`，`--dump` 路径此前从未被行使过。**smoke 的 dump 没有幸存**（3.1h 算术只剩日志打印）。补测（转录日志区间为 dump 格式，`--load` 重放）**确认了预测**：`report_ratio` 在 ϖ=1000 处准确抛出同一行 `ComplexResult`（`p4_cii_m3_interval.py:941`）。第三次独立复现发生在灵敏度扫描展宽 `V_GU_phi` 时，撞上同一行同一异常——不是巧合，是「W_φ 下端为负+`iv.sqrt`」在 D 足够小时必然发生。**这同时说明 clamp 不是可选的**：没有 clamp，合法但粗的包围会让装配在大 ϖ 处异常退出，而非返回更宽的界。

**关闭方式与理由（只允许这一种写法）**：`cv_quad` 出口对两个 I-常数做 `I∩[0,∞)`。**它可靠是因为真值落在交集里**（不是「只动下端所以上界不变」——下端经减法/取负/除法能到达上界，此处碰巧无害，「碰巧无害」正是危险所在）。非负性证明（模块头诚实边界 (vii)）：(a) `h_W(s)=s³Û(s)²/2>=0`（Û=T_f(s)/s²=TV_σ(ψ) 全变差平方）；(b) `h_GU(s)=G(s)·s·Û(s)/2>=0`（G(s)=1/(1+s²)>0 是定义）。于是 W_f、GU_f 精确非负。

**clamp 与 assert 的分工**：有证明且预期可能触发⇒clamp（两个 I-常数，计数 `SIGN_CLAMP`）；有证明且预期永不触发⇒fail-closed `assert`（V-常数、TBAR 下端、`cv_quad` 分割不变量）。**clamp 的第二个理由，与可靠性论证并列不替代**：三次复现说明没有 clamp 会让合法但粗的路径抛异常而非返回更宽的界；clamp 可靠性来自「真值在交集里」，clamp 必要性来自「让宽包围以合法形式呈现」。

**披露**：`--load` 不经过 clamp，灵敏度扫描基线 W_φ 下端仍是 dump 里的 −0.0229，扫描脚本**直接测量** clamp 前后 ϖ_A″ 差（见 §11），不论证它。

**这是同一病因的第三例**：`_abs_range` 一类函数返回退化 `ivmpf` 而非 `mp.mpf`——同一歧义已产出 (1) `fam_sups` 的 fail-open（§6.2）、(2) `tbar_bb` 的 松紧问题（§6.2）、(3) 本次 `dump_constants` 崩溃。**三处补丁不是修复，只是推迟第四次**：真正修复要在**产生处**消除歧义（所有标量出口统一返回 `mp.mpf` 或统一返回 `ivmpf` 并加类型断言），本里程碑未做（要动 M3/M4 共享工具，记在 §12）。

---

## 10. 受信任基座（不是空的，如实列出）

| 原语 | 状态 | 证据 |
|---|---|---|
| `iv.sqrt` | SOUND | 0/40000 未包围 |
| `iv.log` | SOUND | 0/120000 未包围（含 M4 用的 `log(1+s²)`） |
| `iv.pi` | SOUND | 端点与 `libmp.mpf_pi` 逐位相同 |
| `iv.atan2` | **UNSOUND，已替换** | 2502/40000随机+1429/20000系统操作数未包围，最差 1.41e-32 |
| `mp.sqrt`（M3 `V_FPP_iv` 内） | 带 pad 的 float 路径 | 27 次调用，dps 40 坐在 **1e-12 绝对 pad** 下，比自身舍入误差高 28 个数量级——**sizing error 不是 design error**（pad 位置对，只是选大小时没算过） |

`check_transcendental_unreachable` **报告**调用计数而非断言，新原语会显式冒出来而非悄悄加入基座。

---

## 10b/10c/10d. 18 个对象的灵敏度、推导侧归类、生产宽度重量

**基线 ϖ_A″（smoke）= 9.769080**（`--maxc=400,rel=2e-1,top=200,s0=64`）。灵敏度做法：每次只动一个对象，(W) 上端再加自身宽度、(C) 塌缩到下端（`m5_sens.py`/`.log`）。

**smoke 表（worst/best 行）**：

| # | 对象 | 相对宽 | d(W) | d(C) |
|---|---|---|---|---|
| 1 | V_GU_phi | 0.899 | +0.033887 | −0.034320 |
| 4 | TBAR[kap] | 0.819 | +0.008823 | −0.008849 |
| 12 | W_G | 0.018 | +0.000361 | −0.000364 |
| 14–18 | GU_G, V_GU_G, TBAR[UG/Ukap/Uphi] | — | **+0.000000** | **−0.000000** |

d(W)之和+0.088977，d(C)之和−0.099335；**联合塌缩 9.653067（d=−0.116013）> 逐项d(C)之和**——装配非线性，逐项和不是联合下限，不能当界用。**哪个 cap 在卡**：6个sup合计仅占余量17%，12个CV常数占83%；kap/phi 撞900盒上限合计仅值+0.0125⇒**`max_boxes` 保持900，`--maxc` 才是该调的旋钮**（新增 CLI `--maxbox=N`）。

**那 5 个 0 是结构性的，不是舍入**：`bounds3` 只在特定行读 `TBAR3['G'|'phi'|'kap']` 与 `CV3['GU_'+nm]`（nm∈{phi,kap}），`TBAR3['UG'|'Ukap'|'Uphi']`、`GU_G`、`V_GU_G` **一次都没被读**。三处数值佐证：展宽1e6倍 d=+0.000000；只塌缩13个活对象与18全塌缩完全相同；负控制 V_GU_phi 展宽1倍即动 +0.033887。**协调者定级必须原样保留**：

> A structurally unread object is one of two things: a harmless leftover from an earlier assembly, or **a term the derivation says should be there and the code dropped**. Those look identical from inside `bounds3`... and only one of them is fine. The second is fail-open and is exactly the shape of defect this whole exercise exists to catch.

**§10c 核查结论：(a)/(a)/(a)/(b)/(b)，无一为 (c)。** 只读推导给出正面可数事实（`P4_einstein_extremal.md` §1.G.3 pad 目录 11 行，逐行用处）：T̄ 占2行3实例（G/κ/φ₀，无一行 T̄_{U_f}）；`I[G·U_f]` 占2行（φ/κ，G族走ℓ²/Minkowski无 `I[G·U_G]` 行）；G族唯一二次型 `⟨τ,𝓡_Gτ⟩` 在 `:475–477` 因 `𝓡_G` 逐元≥0（(M1)）按符号从下界精确丢掉（原文：「两个被丢项都≥0…无任何离散化pad」），而 `:497–499` 说明 φ/κ 的 𝓡 不逐元≥0（min≈−0.32）「不能丢」——这正是 padX 只对φ/κ存在的原因。三个 U-sup 因 `Σw_n` **发散**使 sup 路线不可用（Ū_f·√(Σw_n)=+∞），无关风格选择。

**丢弃方向性论证**（§10c.1，唯一可能塌的地方）：g_lo=2B−1−Θ* 只出现在判据分母，g_lo↓⟹判据更难满足；经 s=C_r/g_lo↑，F_lo↓、E_hi↑——三条方向一致，**丢弃是 fail-closed 的**。这一步成立，(b) 才成立。

**判决**：五个全部落在(a)/(b)，无(c)。「18个对象里只有13个承重」可以写；**9.769080/9.716960 不需要重评**。文档后果：§1.G.4(3)、`p4_cii_m2.py:116–117` 与 §1.H/§1.I 的52计数把承重集多算了5个（已交回协调者，M3章节非本里程碑改）。§1.H(:829) 称6个sup是ϖ_A″=9.71而非8.48主因——与测量矛盾（3个恒零，6个加起来不到 V_GU_phi+V_W_phi），**本里程碑不改**，已交回。

**§10d 生产宽度重量**（`_scratch_ir/srv_m5_sens_prod.py`，载生产 dump `m5_final.json`，走最终冻结脚本 `3fdcfc7e`）：基线与两次还原精确复现 **9.716960**，dump 里**无负下端**（生产不需要 clamp）。

| 对象 | 相对宽 | d(W) | d(C) |
|---|---|---|---|
| V_GU_phi | 0.783 | +0.014608 | −0.014689 |
| **T̄[κ]** | 0.819 | +0.008987 | **−0.009014** |
| GU_kap | 0.00628 | +0.000016 | −0.000016 |
| GU_G,V_GU_G,T̄[U_G/κ/φ] | — | **+0.000000** | **−0.000000** |

三条结论：(1) 五个惰性对象在**交付宽度**上仍精确为零——§10c 推导侧论证的独立确认（仍只是佐证，(c)排除靠可数事实非此表）。(2) **排序变了**：`T̄[κ]` 升到第2——`--maxc` 加四倍不影响 box-cap-set 的它，**下一档该调 `--maxbox` 不再是 `--maxc`**。(3) (R14) 仍成立：逐项d(C)和−0.044611 vs 联合塌缩−0.046032（ϖ_A″=**9.670928**）——逐项和低估联合效应。

**最有用的一个数：M5 全部剩余宽度只值 0.046。** 对 M2 float64 真值 8.4761，差距 1.241——**M5 宽度占 3.7%**，由区间算术包含单调性，任何子集（含六个T̄ sup）都被 0.046 控制。**⟹ §1.I.6(2)/§1.I.7「六个T̄ sup是主因」数量上不成立**：三个U-sup恒零，另三个连同12个pad常数一起最多值3.7%，剩下96%在M3/M4自身保守性与判据装配里。

---

## 11. 最终数字

**头条：严格 ϖ_A″ = 9.716960。** 46个求积常数+6个网格sup（合52对象）全部Tier 1，无float64余量记账。与段C-I重叠4.631×（认证到45）、12.350×（认证到120）。

沿革：M2 float64 8.4761（非界）→ M3 9.831910（界，36对象带余量）→ M4 9.711925（界，18对象带余量）→ **M5 9.716960（第一个每项都由包围支撑）**。M5比M4退0.005035——退这一步是把乐观换成界，四数不是同类物不能当收敛曲线读。

| | 构建 82cca55f | 报告 3fdcfc7e |
|---|---|---|
| checks | 68 passed/0 failed/**0 not exercised** | 66/0/**3 not exercised** |
| 墙钟 | 16317.2s（32进程） | 1515.3s（`--load`） |
| 剖面求值 | 92,013,914 jets/23,777,896 σ格 | 0/0 |
| sign clamp | **0x**（挣来的真零） | not exercised |

3条 not exercised 恰是 `--load` 阶段不能证伪的三条（cap不相交、`_isect`不相交、clamp计数）。

判据表：`w=9.7170 ratio[0.764429,1.000000] ok at 9.716959756753118`；`w=10 [0.748197,0.971564]`；`w=45 [0.428009,0.481449]`；`w=1000 [0.320613,0.339820]`。

Θ*对照：M3出厂 `[0.6349287917881079046367176]` width 0（**不包含**）；M5认证 `[...529866324]` width 4.50e-29；两者装配出的 wA'' 均=9.716960（δ=+0.000e+00，**数值巧合非可靠性论证**——M3出厂宽度为零且漏真值4.834991e-17、方向向下，替换的理由是这条，不是判据动没动）。

18对象交付表（worst 行摘录；包围为dump-derived，9条退出中6条log-printed、3条log-derived）：

| 对象 | 相对宽 | 要求 | exit.src | 退出 |
|---|---|---|---|---|
| W_G | 5.51e-03 | 0.05 | derived | tolerance met |
| V_W_phi | 9.19e-01 | 0.05 | derived | cell cap 1500 |
| TBAR[kap] | 8.19e-01 | 0.01 | printed | **box cap 900** |
| TBAR[phi] | 9.08e-01 | 0.01 | printed | **box cap 900** |

**9条退出里4条cap-set、5条容差达标**（cv:kap、cv:phi、tb:kap、tb:phi）。

**§10b预测对了一半**：`max_boxes`保持900、该调`--maxc`——旋钮选对（400→1500，赚0.0521）；但按smoke灵敏度排前四的对象（V_GU_phi/V_W_phi/V_W_kap/T̄[κ]）恰是几乎没动的四个（−13%/−7%/−15%/0%），改善66–94%的全在后面。κ/φ两族8个CV常数**仍骑cap**。**把这两件事记成同一次成功是不诚实的。**

**正式同一性核销**：`82cca55f`与`3fdcfc7e`同tiny命令行各跑一次，**24条目逐位相同、0处不同**（`019bc17d`那次只作记录）。

**退出重建核验**：`m5_exitlog_chk.py` 四组**35项全过**；SUITE C 跑完整生产日志，tb侧4converged/2box cap真混合。

---

## 12. 还没关上的

1. **C10**：源码守卫不覆盖M3/M4；已在守卫输出与模块头如实打印。
2. **可达性守卫是动态的**（§3(b)原话）。
3. ~~`tbar_bb`的cap退出路径在生产容差下未被行使~~ ——**已被生产日志证伪，留痕不删**：`TBAR[kap]`走盒上限退出 `[0.45352608,2.50708595]`（与smoke预测`[0.453,2.508]`一致），同一次运行`TBAR[G]`/`TBAR[UG]`在442/329盒收敛退出——tb侧是真混合。按(R10)cap退出是松紧风险非可靠性风险。
4. **`M3.V_FPP_iv`的27次`mp.sqrt`**仍在1e-12 pad下（§10）。
5. **M3 Tier-1内部8个V-配对仍是200%相对宽度**（M3 §1.H.3已证不承重，本里程碑未动）。
6. **退化`ivmpf`与`mp.mpf`的类型歧义未在产生处消除**（§9b）——已产出3个缺陷（fail-open一个/松紧一个/崩溃一个），每处单独修好但病因还在；要在产生处关掉需改M3/M4共享工具函数，超出本里程碑范围。
7. **5个对象在装配里不承重**（§10b量出、§10c推导侧定性）：`GU_G`、`V_GU_G`、`TBAR[UG]`、`TBAR[Ukap]`、`TBAR[Uphi]`被`bounds3`完全不读；定为(a)/(a)/(a)/(b)/(b)，无一为(c)。它们仍被认证/检查/记账——**这本身不是缺陷，但46+6的计数把承重集多算了5个**，对应文字在M3/M2章节，已交回协调者。
