# P4 · 段 C-I 区间化 —— 对抗性核查记录

日期：2026-09-12　worktree `agent-a64367a0e7fbab99b`
被审对象：`theory/derivations/p4_cii_ci_interval.py`（段 C-I 的区间算术证书）
基线（只读对照，零改动）：`theory/derivations/p4_extremal_concavity.py`、`p4_interchange_uniformity.py`
计划与断言清单：`theory/notes/2026-09-12_ci_intervalization_plan.md`
记账：`theory/derivations/P4_einstein_extremal.md` §2.A

立场：**假设这份区间化有洞**。分两路——(甲) implementer 自己的**变异测试**（故意把已知的 bug 类重新种回去，看门是否报警）；(乙) reviewer 子 agent 的独立审查。

---

## 甲. 变异测试（implementer，`notes/_scratch_ir/mutate.py`，未提交）

方法：只读 import 被审模块，在**内存中**替换某个函数/常数为「带已知 bug 的版本」，跑 `run_gate(full=True)`（或 `run_L`），记录门是否报警、以及**报警的那条检查是否语义正确**（报警但报在无关的检查上，说明门只是碰巧撞上，不算真覆盖）。

| # | 变异 | 语义 | 是否被捕获 | 捕获它的检查 |
|---|---|---|---|---|
| M1 | `PSI_P` → `{0:1,1:2,2:3,3:4,4:5,5:6}` | **历史 bug #1**：包络指数写错 ⟹ 包络**不**支配导数 | ✅ | `psi_j dominates \|d_w^j [[k]]\| pointwise`（j=1, w=0.05, k=2） |
| M2 | `_L1_hint` 用 p=2 | **历史 bug #2**：选格启发式指数写错 ⟹ 粗步长 10× 过宽 | ✅ | `_L1_hint tracks the certified rsum_seg to within a factor of 1.5`（w=1.0，hint 0.898 vs 认证 2.0） |
| M3 | `PAD = 1e-12` | 把旧的 `ROUND_PAD` 种回去 | ✅ | `PAD is identically zero` |
| M4 | `sub_lo(x,y) → x - y`（裸 mpf，**就近**舍入） | 单边性被破坏 1 ulp | ✅（**修复后**） | `sub_lo rounds OUTWARD: 1 - 2^-200 lands strictly below 1` |
| M5 | `chol_pd` → 「没抛异常即通过」 | 把旧的正定性论证种回去 | ✅ | `chol_pd REJECTS a singular (psd, not pd) matrix` |
| M6 | `LB` → `UB`（全局） | 单边塌缩方向翻转 | ✅ | `x**2 of a sign-straddling interval is [0,4]` |
| M7 | `cell_table(f_b,f_a)` → `cell_table(f_a,f_b)` | 细胞一致包夹方向翻转 | ✅ | `endpoints must be properly ordered`（`ivlh` 的构造断言） |
| M8 | `tail_int → 0` | **有限截断不补尾** | ✅ | `tail_int(A=0,p=2,x=5) = 1/375`（闭式对拍） |
| M9 | `subdivide` 把最后一个节点缩 2^-16 | 铺砌留缝 | ✅ | `tiling: consecutive fine cells share a bit-identical endpoint`（0.3749847412109375 vs 0.375） |
| M10 | `subdivide` 把全部内节点**整体**平移 2^-16 | 「另一种合法铺砌」 | ➖ 不报警 | **正确行为**：平移后细胞仍逐位相接、仍覆盖 [A,B]，是一个合法的不同铺砌，不该报警 |

### M4 是这一轮唯一的真实发现（已修，`c18c0f1`）

旧的四条检查写成

```python
chk("sub_lo rounds OUTWARD (never above the exact difference)",
    sub_lo(mp.mpf(1)/3, mp.mpf(1)/7) <= mp.mpf(1)/3 - mp.mpf(1)/7)
```

右端 `mp.mpf(1)/3 - mp.mpf(1)/7` 本身就是**就近舍入**的 mpf 值，于是「朴素实现 `x - y`」会使左右两端**逐位相等**，`<=` 平凡成立——**这四条检查对它们要防的那个 bug 完全没有约束力**。（注意：被审的实现本身是**对的**（`sub_lo = LB(ivf(x) - ivf(y))`），这是**门的盲区**，不是证书的洞；但既然全节的论点就是「端到端外向舍入」，门必须能证伪它。）

修法：把四条检查改成**判别性**的——操作数选得让精确结果在 prec 103 下**不可表示**、且**就近舍入朝错误方向走**，参考值在 `mp.extraprec` 下计算（mpf 之间的比较跨精度是精确的），并**同时断言朴素实现会失败**：

| 组合子 | 操作数 | 精确结果 | 就近舍入给出 | 判别式 |
|---|---|---|---|---|
| `sub_lo` | 1, 2^-200 | 1 − 2^-200 | **1**（偏大，违反下界） | `sub_lo ≤ exact` 且 `(1−2^-200)_nearest > exact` |
| `add_hi` | 1, 2^-200 | 1 + 2^-200 | **1**（偏小，违反上界） | `add_hi ≥ exact` 且 naive `< exact` |
| `mul_hi` | (1+2^-100)² | 1+2^-99+2^-200（201 bit） | 1+2^-99（偏小） | `mul_hi ≥ exact` 且 naive `< exact` |
| `div_hi` | 1/3 | 1/3（不可表示，prec 103 处**向下**舍入） | 3·(1/3) < 1 | `3·div_hi(1,3) > 1` 且 `3·(1/3)_nearest < 1` |

（`div_hi` 这一条自身还踩了一个精度上下文坑：若把 `3*(1/3)` 放在 `mp.extraprec(20)` 内计算，`1/3` 会在 prec 123 下重算、乘 3 后又被舍回 123 位、恰好等于 1，断言假失败。已改为在基准精度下先取 `_third = 1/3`，再在 `extraprec(60)` 内做 `3*_third`（≤105 位，精确）。）

修复后重跑：M4 被 `sub_lo rounds OUTWARD` 精确捕获，其余 9 条结论不变。本机 `(0,3]` 全档 **1093 检查 / 71.7 s**（Stage A′ 635 + Stage S 23 + Stage L 435），16 粗节点 / 39 细胞，最小认证枢轴 5.19445e-8，承重细胞 [1.02246, 1.11523] pos/F = 0.5536——与修复前逐位相同（改的只是门，不是数学）。

---

## 乙. reviewer 子 agent 的独立审查

立场：「假设这份区间化有洞」。reviewer 独立跑通 `GATE`（当时 635 检查 / 32 s）与默认档 `(0,3]`（1093 检查 / 76.5 s，4 进程），另写 13 个探测脚本 `theory/notes/_scratch_ir/pr1.py … pr13.py`（未提交）逐项证伪。

结论：**发现 1 条必修项**（一个区间基元在 mpmath 1.3.0 里并非外向舍入，且落在承重路径上），其余 6 条为「应修」（门的盲区与措辞），另列出 7 项「门未覆盖但实测通过」的盲区与一份「可接受的非区间成分」清单。**数值结论本身没有被动摇**：违规量级 1e-30 相对，而结构余量 38%。

### 乙.1 必修项 M-1（已修，见下）

**`iv.atan2` 不是外向舍入的** —— `tail_int` 的 `A>0` 分支因此不是「可证包含」，而 Stage L 的每个细胞都经它承重。

- 位置：原 `p4_cii_ci_interval.py` 中 `at = iv.atan2(A, x)`；下游 `S1_of` → `tail_sum` → `row_env`/`col_env` → `rsum`/`rsum_seg` → `schur_coupling`/`coupling_restricted` ⟹ `alpha_cell_iv`、`L1`、`N_ub`、`R3_ub`、`rtail`——即 Stage L 的 α、g_lo、E_hi、F_lo **全部**间接依赖。Stage S 不受影响（`A=0` 分支是初等的 π/2 常数，但那时用的是 `iv.pi`，同样被一并替换）。
- 原因（mpmath 源码 `libmp/libelefun.py`）：`mpf_atan2` 先算 `mpf_atan(mpf_div(y, x, prec+4), prec+4)`（**就近舍入**），只在最后一步 `mpf_pos(tquo, prec, rnd)` 施加方向舍入。方向舍入作用在一个已带 ±2^-(prec+4) 误差的中间值上，**不保证包含**。
- reviewer 的复现：20 万对随机点区间操作数与 600 bit 参考比较，**9.3–9.4% 的结果不包含真值**，最坏违规 0.235 ulp@prec103；生产操作数类（A 为 [0.3125,45] 内 2^-16 倍数、x 为 1..3000 整数）60000 次中 5735 次不包含；**传到 `tail_int`** 后 20000 次中 14 次不包含，最坏相对 6.2e-30。
- 为什么算洞：模块 docstring (a) 与 §2.A.1 第 1 条的核心论点就是「每个初等运算返回一个**可证**包含精确实数结果的区间」。这条对 `atan2` 不成立；它恰恰是本节声称消灭的「未量化舍入假设」的同类。

**修法（已实施）**：不用「高精度算完再加垫」（那正是本文件要去掉的东西），而是从零写一个只用 `+ − × ÷ √` 与整数幂的**认证反正切** `iv_atan_nn`：

1. 反复用 `atan(t) = 2·atan(t/(1+√(1+t²)))` 把参数减半（`ATAN_HALVINGS = 80` 上限）直到 `UB(u) ≤ 2^-4`；
2. 对交错级数 `Σ (−1)^j u^(2j+1)/(2j+1)` 取 `ATAN_TERMS = 20` 项，余项用**交错级数余项定理**给出对称区间 `±UB(u^41/41)`；
3. 乘回 `2^m`。`iv_pi_half() = 2·iv_atan_nn(1)`。

`tail_int` 的 `x=0` 分支改用 `iv_pi_half()`、`A>0` 分支改用 `iv_atan_nn(A/x)`（原 `iv.atan2(A,x)`，同值）。验证：3000 个操作数对 600 bit 参考 **0 次不包含**（同一测试下 `iv.atan2` 为 200/3000），相对宽度 6.3e-30，π/2 相对宽度 4.5e-30；2000 次 `tail_int` 调用 0.35 s（性能无退化）。

门里同时加了 **4 条严格自洽恒等式 + 2 条证据检查**：`atan t + atan(1/t) = π/2`（6 个数量级尺度）、倍角 `atan(2t/(1−t²)) = 2 atan t`（3 个）、Machin `π = 16 atan(1/5) − 4 atan(1/239)`、区间宽度 ≤ 1e-27 相对、以及 240 个操作数对 600 bit 参考的包含性。前四类是**严格**的（只用 `iv_atan_nn` 自己，不引用任何外部高精度值），足以钉死减半次数、项数与余项常数。

### 乙.2 应修项 S-1…S-6

| # | reviewer 的意见 | 处置 |
|---|---|---|
| S-1 | `check_tiling` 的「ends at reported W」是同义反复（`Wdone = cells[-1][1]` 后再 `check_tiling(cells, W0, Wdone)`）；若某粗细胞 6 次二分仍不闭，`stopped` 被设置但 `main` 不断言，RESULT 会**静默报更小的 W** | **已修**：`run_L` 中新增 `chk("the ladder reached the REQUESTED W_max …", stopped is None and Wdone == float(WMAX))`，并在 `check_tiling` 的调用处注明那两条是**内部一致性**检查。`check_tiling` 本身的「逐位相接 / 端点为 2^-16 倍数 / 非退化」三条仍是实质检查 |
| S-2 | 门的 (I1) 算子级检查注释说「must CONTAIN 0」，实现却是 `bad < 1e-25` | **已修**：改为对每个矩阵元累积 `_contains0(d)` 的真包含判定 |
| S-3 | 需要一条基元级检查：`iv.atan2` 端点 vs 600 bit 参考，要求 100% 包含（当时会立刻失败） | **以更强的形式实施**：`iv.atan2` 已被彻底移除，改为 (i) `iv_atan_nn` 的 240 操作数包含性证据检查，(ii) `iv.sqrt`（现在**唯一**留在承重路径上的库基元）的 300 操作数包含性证据检查，(iii) 一条**源码级回归闸**（见乙.4） |
| S-4 | docstring/§2.A 的「every elementary operation … PROVABLY contains」措辞过强 | **已修**：docstring (a) 现在点名基元集合 `+ − × ÷ √`、整数幂、min/max，并**明确写出** `iv.atan2` 的失效机理与实测数据、以及 `iv_atan_nn` 是自己造的；§2.A.1/§2.A.5 同步 |
| S-5 | `rsum_seg` docstring 说「additionally intersect the bulk with the rsum bulk bound」，代码没做 | **已修**：删去失实句 |
| S-6 | Stage A′ 没有 `cell_S_iv` vs `C0.cell_S` 的对拍（reviewer 实测 [0,0.3125]、[0,0.25] 十个量一致到 ≤3.4e-9 相对） | **已加**：Stage A′ (g2) 段，两个细胞 × 12 个量逐项对拍，外加一条「iv 的结论不得**优于** mpf 的结论」的保守性检查。实测 hpp_ub：mpf −4.36902975e-01 / iv −4.36902973e-01，mpf −1.50004212e+00 / iv −1.50004212e+00 |

### 乙.3 可接受的非区间成分（reviewer 逐条核过，与 §2.A.5 一致）

| 入口 | 用途 | 为什么不是界 |
|---|---|---|
| `np.linalg.eigh` / `trial_vec_np` | 试探向量 v̂、紧缩向量 u | 任意固定向量给合法 Rayleigh 下界与 Davis–Kahan 角；`λ₂(A) ≤ λmax(PAP)` 对任意 u 成立（x⊥u ⟹ Px=x）。分量以精确二进制进入 `ivf(float(t))`，归一化在 iv 内 |
| `c0=float(ev[-2])`、`CHOL_C0` | 区间 Cholesky 位移 c 的**种子** | c 本身是被区间 Cholesky 认证的那个精确 mpf；种子错只会触发升级或失败。reviewer 用精确奇异 PSD 整数 Gram 矩阵 ±1e-25/±1e-40 位移测 `chol_pd`：全部拒绝，+1e-15 接受 |
| `_L1_hint`、`support_size`、`na_band`、`coarse_grid` 的步长 | 选格 / 选维数 | 只决定**认证哪些**不等式，不决定它们是否成立；端点经 `int(d*1024)/1024` 落到 2^-10 格，`is_dyadic` 逐个验 |
| `repr(float(a))` 缓存键 | 缓存 | a,b 本来就是 float 端点 |
| `info` 字典里的 `float(...)`、`ratio` | 报表 / 细化判据 | `ratio > TARGET` 只决定是否再细分；`hpp_ub < 0` 用 mpf 判定 |
| 细胞中点 `p = 0.5*(a+b)` | 展开点 | 断言 `p−a == hw == b−p` 精确成立 |
| Stage A′ 的全部 `np.` / `mp.quad` / `mp.zeta` | 证据 | 不进入任何界（并被乙.4 的源码级闸挡在 Secs 0–6 之外） |

reviewer 另把 `head_iv` 与 float64 `op_np` 对拍（Na=6,30；w=0.3125,3,44）：区间中点与 f64 差 ≤2.2e-16，**确认头块确实在 iv 内从核表构造，不是 f64 拷贝**——这正是 (b) 条要消灭的第二个洞。

### 乙.4 门的盲区（reviewer 列了 7 条；处置）

reviewer 对这 7 条都做了独立探测，**当时的实现全部通过**，但门里没有对应检查：

**最终处置**：1、2、4、5、6、7 已加检查（门从 635 条涨到 854 条）；**只有第 3 条（L₁ 的细胞内部对拍）仍然敞开**，理由写在该条下。

1. **`eph_hi/eps_hi` 的导数包络常数 4b、6b**：`|∂_w ψ_w(k)| ≤ 4b/(k²+a²)³`、`|∂_w φ_w(k)| ≤ 6b/(k²+a²)³` 完全没被门检查——正是 `PSI_P` 那类「常数写错 ⟹ 包络不支配导数」的位置。实测最大比值 0.99999987 / 0.99999973：合法但**极紧**，把 6 改成 4 门不会报警。→ **已加** 56 条逐点支配检查（4 个代表细胞 × 7 个 k × 2 个包络），与 (e0) 同型。
2. `coupling_restricted`、`tail_sq`、`row_env`/`col_env` 无 brute-force 对拍。实测 `‖B_restricted‖` 界/真值 ≈ 5.8×、`tail_sq` ≈ 17×、`row_env(N)` 覆盖 n∈[N,N+50) 行和，全部支配。→ **已加**：取包络**恰好取等**的核 f(k)=C/(k²+A²)^p（三组 (C,A,p) = (1,0,2)、(2,1,2)、(6,0.5,3)），在 float64 里造 420×420 的 𝒪[f]，验 `row_env0`/`row_env(n)`/`col_env(m,N)`/`schur_tailblock`/`schur_coupling`/`coupling_restricted`/`tail_sq` 支配**截断后的**行和 / 列尾和 / 块谱范数 / 尾平方和（截断只丢非负项，故支配截断是必要条件）。实测比值：

   | 量 | 界/截断真值 |
   |---|---|
   | `row_env(n)` | 1.18 – 8.9 |
   | `col_env(m,N)` | 1.10 – 35.8 |
   | `schur_tailblock(Na)` = ‖D‖ | 1.30 – 1.49 |
   | `schur_coupling(Na)` = ‖B‖ | 4.2 – 13.5 |
   | `coupling_restricted(Ns,Nr)` | 5.0 – 9.8 |
   | `tail_sq(spec,M)` | 5.8 – 55 |

   前三个紧到 1.1–1.2，**确有牙齿**：把 `row_env` 或 `col_env` 整体乘 0.5，门立刻在这两条新检查上报警（已验证）。后三个松，只是量级卫兵：乘 0.5 不由它们报警，而是由「同输入逐量对拍」与「Stage S 保守性」捕获（也已验证）。这一保留写进了 §2.A.5 第 6 条。
3. L1 只在端点 b、N=400 中心差分处对拍；reviewer 核了细胞内部与 N=1200：L1/‖H′‖ ≥ 1.8×，四个细胞都支配。→ **未加**（余量 1.8×；且 `rsum_seg` 的组成部分 `row_env`/`col_env`/`tail_int` 现已被第 2 条的暴力对拍覆盖）。
4. `chol_pd` 缺「浮点 Cholesky 会误放行」的对抗用例。→ **已有**：门里 `chol_pd REJECTS a singular (psd, not pd) matrix`（变异 M5 即由它捕获）。
5. `cell_L_iv` 对**输入**的安全单调性（α↑、L1↑、μ₂↓ 应使 `hpp_ub` 不降）无检查。→ **已加**：Stage A′ 在 [1.0,1.125] 上三向各扰动 5%，断言 `hpp_ub` 不得改善。
6. `alpha_cell_iv` vs 细胞内 f64 λ₂（N=1200）、`mu2_cell_iv` vs 细胞内 λ₁(H⁽²⁾(w)) 无内部点对拍；实测全部支配。→ **已加**（Stage A′ 新增 34 条）：

   - `mu2_cell_iv`：4 个细胞 × 5 个 w（t = 0, ¼, ½, ¾, 1）× 「μ₂ ≤ λ₁ 的 N=400 f64 压缩」。压缩只会**降低** λ₁，故支配压缩是必要条件。实测 μ₂/λ₁^f64 = **0.38–0.96**，最紧处（[0.3125,0.375] 的 w=a）只剩 4% 余量——**有牙齿**（把 μ₂ 换成在右移一格的细胞上取值，门立刻在这条上失败，已验证）。但**如实说明**：μ₂ 的段论证依赖 k₁,k₂,k₃ ≥ 0（门里单独断言），在该断言成立时约束点就是 w = a，内部点是冗余的；加它是为了直接检验复合命题而不是它的证明。
   - `alpha_cell_iv`：3 个细胞 × 4 个 w × 「α ≥ λ₂ 的 N=600 f64 压缩」。实测 α/λ₂^f64 = **220–1000 倍**，只是量级卫兵。**这不是疏漏而是 α 的固有性质**：α = max(c + L₁·hw, ‖D‖) + ‖B‖ 由 Schur 包络项主导，而真 λ₂ 比它小 2–3 个数量级（ϖ=1 处 λ₂ ≈ 4e−4 而 α ≈ 0.197），没有任何 f64 对拍能把它压紧。
   - 另加 3 × 2 条：`alpha_cell_iv` 里**唯一未被其它检查触及的那一步**是 Weyl/Lipschitz 覆盖 `λ₂(head(w)) ≤ max(c_A,c_B) + L₁·(B−A)/2`。现在在细胞中点再跑一次头块区间 Cholesky（Na=70/140，0.2–1 s），直接验证该不等式。它对 L₁·hw 项的余量约 10³ 倍（L₁ 界 ‖∂_ϖℌ‖ ~ O(1)，λ₂ 随 ϖ 变化极慢），但**恰好能捕获它要防的那个 bug**——把 Lipschitz 项整个删掉：实测 c(中点) = 4.327671e−5 **>** c(端点) = 3.600594e−5，于是不等式立刻破。
7. 超越基元的包含性 —— 「本次唯一真正失守的一项」。→ **已修**（M-1）+ **源码级回归闸**：Stage A′ 用 `tokenize` 解析本文件自身（因此注释与 docstring 里的散文既不会误报也不会掩盖），断言 **Secs 0–6（一切产生界的代码）里不出现任何 `mp.`/`iv.` 超越函数调用**——白名单只有 `iv.sqrt`（mpmath 的 `mpi_sqrt` 对两个端点分别调用 `mpf_sqrt(·, prec, round_floor/round_ceiling)`，方向舍入直接作用在开方本身，这与 `mpf_atan2` 正相反），且只经 `iv_sqrt_nn` 到达；`mp.sqrt`（就近舍入）照样禁止。这条闸的**牙齿已验证**：把 `_ = iv.pi + iv.atan2(...)` 注入 SEC.3，门立刻以 `iv.pi@684; iv.atan2@684` 失败。

### 乙.5 reviewer 的一句话结论与现状

reviewer 原话：「**可以**称为『段 C-I 的区间算术证书』，但**必须附带限定语**：『除 `tail_int` 中的 `iv.atan2`（…实测约 9% 操作数不包含真值，相对违规 ≤ 1e-29）外，端到端外向舍入』——**或者按 M-1 的修法用严格 atan 上界替换后，去掉限定语**。」

M-1 已按第二条路径关闭（且做得更强：不是只给单调方向的上界，而是给了一个双侧认证的 `iv_atan_nn`），因此该限定语**已可去掉**。reviewer 明确核过而未发现问题的其余各项：单边记账方向、区间 Cholesky 实现与定理陈述、格的精确铺砌、Stage S↔L 在 ϖ=0.3125 的对接、(V1)–(V3) 每个迭代值的独立合法性、`c_k` 跨零时的 `cell_table` 包夹。

---

## 丙. 关闭后的复跑

| 档 | 机器 | 检查数 | 墙钟 | 关键数字 |
|---|---|---|---|---|
| `GATE`（Stage A′ 全档） | 本机 1 进程 | **854** | 50.1 s | 三条 mpf↔iv 对拍与历次逐位相同 |
| `(0,3]` | 本机 4 进程 | **1308**（854+18+436） | 87.5 s | 16 粗节点 / 39 细胞，最小枢轴 5.19445e-8，承重 [1.02246,1.11523] pos/F 0.5536 |
| `(0,45]` | petlove-train 22 进程 | 见 §2.A.4 | 见 §2.A.4 | 与修前逐位相同 |
| `(0,120]` | petlove-train 22 进程 | 见 §2.A.4b | 见 §2.A.4b | 同上 |

顺带修掉一个**计数**（非正确性）缺陷：`_do_node`/`_do_cell` 在串行路径下把 `chk` 记在自己的 `PASS[0]` 上，父进程又加一遍 `nch`，于是 `-j 1` 的 Stage L 计数是 `-j 4` 的近两倍。现在 worker 返回前把 `PASS[0]` 复位到进入时的值，串行与并行报出**完全相同**的检查数（实测 `1.0 FAST`：`-j 1` 与 `-j 3` 均为 118 / 760）。

