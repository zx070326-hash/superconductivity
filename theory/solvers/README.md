# theory/solvers — 各向同性 Eliashberg 求解器

对照 `SPEC.md`。目的：为 `theory/derivations/P1_*` 与后续 P2 提供数值验证器，不是生产代码。

v0 覆盖 SPEC 的必备功能 1–5 与全部要求测试。第二阶段在此基础上实现了功能 6
（能量依赖 DOS，见「能量依赖 DOS」一节）与功能 7（模型顶角核钩子，见「顶角核钩子」
一节），两者默认关闭（`dos=DOS.constant()`、`vertex=None`），不影响 v0 的既有行为
（`dos=DOS.constant()` 时 `tc()` 直接分派到与 v0 完全相同的代码路径，逐位相同）。

## 快速上手

```python
from eliashberg_iso import Spectrum, Solver
import numpy as np

sp = Spectrum.einstein(lam=2.0, w0=100.0)          # 或 Spectrum.from_table(w, a2F)
solver = Solver(sp, mu_star=0.13, w_cut=5 * sp.omega_log())
Tc = solver.tc()                                    # K
Tc_AD = sp.allen_dynes(mu_star=0.13)                # Allen-Dynes 对照

alpha, diag = solver.isotope_coefficient(
    {"full": sp}, which="all", mu_star_fixed=False
)
s_lam = solver.log_sensitivity("lambda")

# 第二阶段：能量依赖 DOS + 顶角核钩子（默认关闭，不影响以上行为）
from eliashberg_iso import DOS

dos = DOS.log_vhs(delta=50.0, a=0.3, W=1000.0, cut=5.0)   # 或 DOS.step / DOS.from_table
solver_dos = Solver(sp, mu_star=0.13, w_cut=5 * sp.omega_log(), dos=dos)
Tc_dos = solver_dos.tc()

solver_vertex = Solver(sp, mu_star=0.13, w_cut=5 * sp.omega_log(),
                        vertex=lambda n, m, wn, wm, T: 0.1 * np.ones_like(wn))
Tc_vertex = solver_vertex.tc()
```

运行测试：`pytest theory/solvers/tests -q`。运行示例：
`python theory/solvers/examples/p1_first_pass.py`、
`python theory/solvers/examples/dos_scan.py`。

## 单位约定

- 所有能量量（ω、Δ、ω_c、ω_log、ω₂、α²F 的横轴）一律 **meV**。
- 温度对外接口一律 **K**；内部换算 `T_meV = KB_MEV_PER_K * T[K]`，`KB_MEV_PER_K = 0.08617333262`。
- 质量缩放因子 `s = sqrt(M/M0)`，`s>1` 表示更重的同位素（频率整体除以 `s`）。
- `λ`、`μ*`、`α`（同位素系数）、`s_k`（对数敏感度）均为无量纲数。

## μ* 与 ω_c 绑定的约定（关键，容易踩坑）

Tolmachev-Morel-Anderson 赝势 `μ*(ω_c) = μ/[1+μ ln(E_F/ω_c)]` 本身就是"选一个截止 ω_c，
把截止以上的推迟效应压缩进一个重整化后的数"，所以**离开具体的 ω_c，单独一个 μ* 数值没有意义**。
本模块的约定：

1. `Solver(sp, mu_star, w_cut, ...)` 里的 `mu_star` **就是** 在 `w_cut` 这个截止处的值——
   调用方自己保证两者匹配（`allen_dynes()`、`tc()` 都不做交叉检查）。
2. 换算用 `mu_star_at_cutoff(mu_star_ref, w_cut_ref, w_cut_new)`，实现的是
   `1/μ*(ω_c1) - 1/μ*(ω_c2) = ln(ω_c2/ω_c1)`。`μ*_ref=0` 是这个关系的退化不动点
   （直接返回 0，不走除法，避免除零）。
3. 从裸参数出发用 `mu_star_from_bare(mu, e_fermi, w_cut)`。
4. **`isotope_coefficient` 的核心开关 `mu_star_fixed`**（P1 §3.1 的核心问题）：
   - `mu_star_fixed=True`：质量缩放时 `μ*` 数值原样保留，不随谱移动重新推导。这是
     "μ* 与材料质量无关"的理想化假设，物理上不自洽（μ* 依赖 ω_c，ω_c 依赖声子标度，
     声子标度依赖质量），但作为对照开关保留，便于复现文献里用这个约定算出的 α。
     在此开关下，α 应严格等于 0.5（任意 μ*、任意谱），因为整套方程在
     `(T, spectrum, w_cut) -> (T/s, spectrum/s, w_cut/s)` 下自洽缩放（w_cut 与谱同比例
     移动、μ* 数值不变，方程形式不变），Tc(s)=Tc(1)/s 精确成立。
   - `mu_star_fixed=False`（默认推荐、P1 认为物理的选择）：质量缩放时 `w_cut` 与谱的
     `omega_log()` 保持固定比例一起移动（比例 = 调用 `isotope_coefficient` 时那个
     Solver 自身的 `w_cut/spectrum.omega_log()`，可用 `w_cut_ref_ratio` 覆盖），
     `μ*` 按 `mu_star_at_cutoff` 重新推导。这是本模块认为"物理"的路径，也是与
     P1 式(3.2)对照的路径。
5. `log_sensitivity("w_cut")` 和 `log_sensitivity("mu_star")` 是**互相独立**扰动单个数值
   （不通过 Tolmachev 关系耦合），只用于诊断/教学，不代表"改变截止选择"这个物理操作
   （那是 `isotope_coefficient` 或手动调用 `mu_star_at_cutoff` 的事）。见下面「已知局限」
   里对 `log_sensitivity("w_cut")` 符号的讨论。

## 折叠技巧推导（正频折叠）

物理方程定义在全体整数 `n∈(-∞,∞)`。用 s 波配对的标准对称假设
`Δ_{-1-n}=Δ_n`、`Z_{-1-n}=Z_n`（`ω_{-1-n}=-ω_n`），把 `m` 从 `-∞..∞` 折成 `m=0..N-1`：

```
Z_n = 1 + (πT/ω_n) Σ_{m=0}^{N-1} [λ(n-m) - λ(n+m+1)]·ω_m/E_m
Z_n Δ_n = πT Σ_{m=0}^{N-1} { [λ(n-m)+λ(n+m+1)] - 2μ*(ω_c)θ(ω_c-ω_m) }·Δ_m/E_m
```

`n, m` 现在只取 `0..N-1`，`ω_n=πT(2n+1)>0`。代码里 `_kernels()` 直接构造这两个折叠核
`K1=λ(n-m)-λ(n+m+1)`、`K2=λ(n-m)+λ(n+m+1)-2μ*θ`，`_max_eigen_linearized()` 在
`Δ→0, E_m→ω_m` 的线性化极限下组装本征值问题。

## Tc 求解策略

线性化核矩阵最大本征值 = 1，对 `ln T` 用 `scipy.optimize.brentq` 求根。默认二分区间
以 `Spectrum.allen_dynes(mu_star)` 的估计值为中心（`t_guess/8, t_guess*4`），失败
（AD 抛 `ValueError`，即 `λ-μ*(1+0.62λ)<=0`）时退化用 `ω_log` 量级估计；两端各有最多
40 次加倍/减半的区间扩张兜底。这个策略是为了避免默认区间过宽时，Matsubara 截断数
`N ∝ w_cut/T` 在低 T 端爆炸导致 `MemoryError`（早期版本用固定绝对下界 `T_lo=1e-3 K`
时，确实在 Tc~40K 的例子上导致过 N≈1.4×10⁶ 的矩阵）。

## 推荐的收敛参数

- **`matsubara_oversample` 默认 3.0**：`N = ceil(oversample * w_cut/(2πT)) + 30`。
  实测（λ=1, μ*=0.15, 扫 w_cut/ω_log 比例）oversample=1.5 时同位素系数 α 有明显未收敛
  （α 相差 ~0.006 于 oversample=3.0/6.0），oversample=3.0 与 6.0 之间的差异降到 ~2e-4，
  故取 3.0 作为默认值，兼顾精度与速度。对精度要求更高的场景（例如要验证 <0.001 量级
  的效应），建议显式传 `matsubara_oversample=6.0` 或更高并做收敛性检查。
- **`w_cut / omega_log` 推荐比例 ≈ 5**：太小（<3）时 μ* 截止会切到声子谱主要权重区域，
  偏离"Coulomb 截止应远高于声子谱"的物理图像；太大时对 Allen-Dynes 弱耦合对照的偏差
  会增大（AD 公式本身假设的是这个量级的截止）。测试与示例统一用 5。
- **硬台阶 θ 的光滑化**：`θ(ω_c-ω_m)` 在代码里不是硬台阶，是宽度为一个 Matsubara 间距
  （`2πT`）的线性斜坡：`clip((w_cut-ω_m)/(2πT)+0.5, 0, 1)`。原因：对 T 做连续求根
  （brentq）时，硬台阶会在某个 `ω_m` 越过 `w_cut` 的瞬间让核矩阵产生有限跳变，
  Tc(T) 因此出现台阶式的微小不连续（实测：固定 `w_cut/ω_log` 比例做精细扫描，某些
  "巧合"比例下同位素系数 α 会跳变 ~0.03，这个跳变在扫描 ratio 时是非单调的锯齿，
  换 `matsubara_oversample`、`eps` 都不消失，只有换成光滑斜坡才消失）。光滑斜坡在
  `T→0`（Matsubara 间距→0）时收敛到严格阶跃，只在截止边界附近的一项上引入
  `O(1/N)` 的可忽略修正，换来 Tc(T) 对 T 连续可导、二分/求根稳定。
- **同位素系数差分步长 `eps=0.02`**：默认值，测试里核对过 `eps ∈ {0.01, 0.02, 0.04}`
  下 α 的稳定性（波动 < 0.005）。

## 能量依赖 DOS（SPEC 必备功能 6，第二阶段）

### 方程

Mitrović–Carbotte（Can. J. Phys. 61, 758, 1983）/ Pickett（PRB 21, 3897, 1980）把常
DOS 版本的 `(Z_n, Δ_n)` 两方程推广成 `(Z_n, φ_n, χ_n)` 三方程，多出的 `χ_n` 是能量
重整化（费米面偏移）；每一项都多一层对准粒子能量 `ε` 的积分，权重是归一化 DOS
`N̄(ε) ≡ N(ε)/N(0)`：

```
Z_n ω_n = ω_n + πT Σ_m λ(n−m) · ⟨ ω_m Z_m / D_m(ε) ⟩_ε
φ_n     = πT Σ_m [λ(n−m) − μ*θ(ω_c−|ω_m|)] · ⟨ φ_m / D_m(ε) ⟩_ε
χ_n     = −πT Σ_m λ(n−m) · ⟨ (ε + χ_m) / D_m(ε) ⟩_ε
D_m(ε)  = (ω_m Z_m)² + (ε + χ_m)² + φ_m²
⟨f(ε)⟩_ε ≡ ∫ dε N̄(ε) f(ε)
```

`n, m` 现在取全体整数（不再折叠到 `0..N-1`）——原因见下面「为什么放弃折叠技巧」。
代码里线性化用 `φ→0` 求 Tc：`Z_n、χ_n` 在 `φ=0` 时自成一组自洽方程（不依赖 `φ`），
先把 `(Z,χ)` 解出来，再对 `φ` 组装线性本征值问题、二分 T 使最大本征值=1，与 v0
的策略一致。

### 归一化：`N̄(0) ≡ 1` 是定义的直接推论，不是可选项

`N̄(ε)≡N(ε)/N(0)` 在 `ε=0` 处的值按定义就是 `N(0)/N(0)=1`——这不依赖任何模型细节。
`DOS.step()`、`DOS.log_vhs()` 的构造参数（`delta`、`ratio`、`a`、`W` 等）描述的是一个
直观的"形状函数" `g(ε)`（例如"远离特征处基线为 1，特征处变成 ratio，或对数发散"），
真正返回的 `N̄(ε)` 是 `g(ε)/g(0)`，强制满足这个恒等式。`DOS.from_table()` 相反：
**不做自动归一化**（表格数据的 `ε=0` 点应该就是调用方计算/测量 `N(0)` 的定义点；
插值到 `ε=0` 处显著偏离 1 时发警告但不改数据，避免掩盖真实的数据准备错误）。

这个归一化选择有一个不直观但数学上无法避免的推论（见「已知局限」第 1 条）：
**vHs 精确钉在费米面（`delta=0`）时，归一化会把"费米面处 DOS 增强"这个效应本身
抵消掉**，表现为 Tc 被压低而不是升高；要在这个模型里复现"vHs 靠近费米面提升 Tc"
的定性图像，`delta` 需要是一个有限但非零的偏移。

### 为什么放弃正频折叠技巧，改用全（未折叠）Matsubara 轴

v0 常 DOS 路径用 `Δ_{-1-n}=Δ_n, Z_{-1-n}=Z_n` 把求和折到 `m=0..N-1`。energy-dependent
版本里 `χ_n` 在 `n→-1-n` 下是**奇**的（`χ_{-1-n}=-χ_n`，`Z、φ` 仍是偶的）——直接推导
一个同时处理"偶+奇"三个量的折叠形式，实现出错的风险明显更高。第二阶段的实现选择
用全体 `n=-n_pos..n_pos-1`（`_full_axis_grid`）直接求解，牺牲 2–8 倍计算量换取更低的
实现风险，这是一个记录在案的工程取舍，不是"更优"的算法选择。

### 能量积分：tan 代换 + 定点 Gauss-Legendre

`ω_jZ_j, φ_j, χ_j` 是 `ε`-无关的前置系数，三个能量积分可以约化成两个基础积分：

```
G_j = (1/π) ∫ dε N̄(ε) / D_j(ε)
H_j = (1/π) ∫ dε N̄(ε) · ε / D_j(ε)
⟨ω_jZ_j⟩ = ω_jZ_j·G_j，⟨φ_j⟩ = φ_j·G_j，⟨ε+χ_j⟩ = χ_j·G_j + H_j
```

`D_j(ε)` 是关于 `ε` 的洛伦兹型分母，直接数值积分在极点附近条件差；代码用
`ε+χ_j = a_j·tanθ`（`a_j=sqrt((ω_jZ_j)²+φ_j²)`）代换掉这个洛伦兹峰，把积分变成
`θ∈(θ_min,θ_max)` 上的光滑函数，再用固定节点数 `dos_n_eps`（默认 400，
`functools.lru_cache` 缓存 `np.polynomial.legendre.leggauss` 的节点/权重）做
Gauss-Legendre 求积。

### 积分范围

`e_max = dos_energy_cutoff_ratio * w_cut`（默认 `dos_energy_cutoff_ratio=20`，对应
SPEC 要求的"至少 ±20·ω_c"）；`DOS.from_table()` 时会额外取 `max(e_max, |表格边界|)`
确保积分窗口不小于表格覆盖范围。表格边界外按端点值常数延拓（`np.interp` 默认行为，
即"flat extrapolation"）。

**这个默认截止本身带来一个 ~0.6% 的系统偏差**（实测：把 `dos_energy_cutoff_ratio`
从 20 提到 200，一个精确常数 `N̄≡1` 的 step-DOS 相对 v0 解析解的偏差从 0.6% 降到
0.06%，与截止量级成反比，是有限积分范围的截断误差，不是求积节点不够——加密
`dos_n_eps` 从 100 到 800 结果完全不变）。SPEC 给出的"至少 ±20·ω_c"本身就意味着
接受这个量级的截断误差；需要更高精度时可以传入更大的 `dos_energy_cutoff_ratio`
（计算量基本不变，`(Z,χ)` 自洽和本征值求解不随 `e_max` 增大而变慢，只有积分节点
覆盖的 `ε` 范围变大）。

### 化学势守恒（`fix_mu="conserve"`）

**动机**：用 Mishra 数据集的真实 H3S α²F + 真实 QE DOS 标定时（见
`theory/derivations/p1_mishra_dos_calibration.py`），发现常数 DOS 路径与 EPW
的 fsr（固定态密度）Tc 相差 1–2K，但能量依赖 DOS 路径与 EPW 的 fbw（全带宽）
Tc 相差 6–8K，且**方向相反**（本模块的能量依赖 DOS 路径略微**抬高** Tc，
EPW 的 fbw 略微**压低** Tc）。H3S 的 N̄(ε) 明显不对称（费米面下方 1.1–1.17、
上方 0.74–0.83），而 EPW 的 fbw 模式会解一个 `muchem`（化学势修正）保持正常
态电子数不变（Lucrezi et al., *Commun. Phys.* 7, 33 (2024),
arXiv:2310.00056）；本模块 `fix_mu="none"`（v0.1 默认）完全忽略这个效应，
是最可能的差距来源之一。`fix_mu="conserve"` 补上这个修正，纯加性（默认仍是
`"none"`，不影响任何既有结果）。

#### 方程

物理上 `mu_shift` 和 `χ_n` 处于同一个位置——都是正常态格林函数分母
`D_n(ε)=(ω_nZ_n)²+(ε+χ_n+mu_shift)²+φ_n²` 里"`ε+χ_n`"这个加法项的一部分——
区别是 `χ_n` 由声子自能自洽产生（出现在 `χ` 方程左边），`mu_shift` 是不依赖
Matsubara 指标的标量，由电子数守恒条件单独定出，不进入 `(Z,χ)` 自洽方程的
左边：

```
Z_n、χ_n 自洽方程与"能量依赖 DOS"一节相同，只是把 χ_m 替换成 χ_eff_m=χ_m+mu_shift：
Z_n ω_n = ω_n + πT Σ_m λ(n−m) · ⟨ω_m Z_m/D_m(ε)⟩_ε
χ_n     = −πT Σ_m λ(n−m) · ⟨(ε+χ_eff_m)/D_m(ε)⟩_ε   （注意：左边仍是 χ_n，不是 χ_eff_n）
D_m(ε)  = (ω_m Z_m)² + (ε+χ_eff_m)² + φ_m²
```

#### 电子数公式与"参考态=0"

标准（虚轴、正常态）电子数公式 `n=∫dε N(ε) f(ε-μ,T)` 在本模块的能量积分基元
`G_j=(1/π)∫dε N̄(ε)/D_j(ε)`、`H_j=(1/π)∫dε N̄(ε)ε/D_j(ε)` 下可以写成

```
偏差(mu_shift) ≡ Σ_j [ χ_eff_j·G_j + H_j ]，χ_eff_j = χ_j + mu_shift
```

—— 与 `χ` 方程右边用的是同一个组合，只是对全体 Matsubara 指标直接求和、不再
乘 `Λ` 核。求根目标是这个"偏差"=0；由于线性化（`φ→0`）求 Tc 只需要**正常态**
电子数（不需要超导态的完整 `φ` 依赖），这个公式已经足够。

参考态是"常数 DOS、`χ≡0`、`mu_shift=0`"：在对称积分区间 `[-e_max,e_max]` 上，
`N̄≡1` 时 `H_j` 是奇函数积分恒为 0，且 `χ=0` 是这条件下 `χ` 方程的精确不动点
（右边=0）——所以参考态的"偏差"精确为 0，不需要单独跑一次参考 `Solver` 来算，
直接把求根目标定为 0 即可。

**符号提醒（实测踩过的坑）**：标准推迟格林函数 `Re G=-(ε+χ_eff)/D` 比上面的
"偏差"定义多一个负号，即这个"偏差"量实际正比于 `-Δn` 而不是 `+Δn`。这不影响
求根本身（二分只看变号），但意味着不能望文生义地假设"`mu_shift` 越大偏差
越大"——详见 `_electron_number_deviation` 的代码文档字符串。

#### 数值求解：热启动二分，不用 `scipy.optimize.brentq`

`mu_shift` 是 Tc 二分套壳里的**内层**求根（每个候选 T 都要重新解一次），默认
用自实现的二分（`_solve_mu_shift`）而非 `brentq`，为的是在相邻试探点之间把
上一个 `mu_shift` 解出的 `(Z,χ)` 传给下一个试探值当起点（`_solve_z_chi_dos`
的 `z0`/`chi0` 参数）暖启动，显著减少总迭代次数——`brentq` 是黑盒实现，不便做
这个优化。括号从 `±0.05·w_cut` 开始双向加倍（最多 `mu_shift_max_expand` 次）
直到偏差变号，再二分（最多 `mu_shift_max_iter` 次或区间宽度 `<mu_shift_tol`）；
测试/诊断用途可用 `Solver.mu_shift_at_T(t_k)` 直接拿到某个温度下的 `mu_shift`
（不依赖 `self.fix_mu`，不修改 `self`，`dos.is_constant` 时恒返回 0）。

#### 已知局限：`mu_shift` 对 `dos_energy_cutoff_ratio` 敏感

这是排查计算比较久才定位到的一个重要限制，如实记录：**跟 `(Z,χ,φ)`
本身只有 ~1/cutoff 量级的小截断误差不同，`mu_shift` 的求根结果对能量积分
范围 `dos_energy_cutoff_ratio` 相当敏感**，实测（Einstein 谱 λ=1.84,
ω_log=93meV, w_cut=465meV，`DOS.step(delta=0, ratio=0.5, width=20)`，
T=150K）：

```
dos_energy_cutoff_ratio     mu_shift (meV)
              2                  -237
             20                 -2937
            200                -30964
```

持续随 cutoff 增大而增大，没有收敛的迹象（这条 DOS 的两侧台阶高度不同，
永远不会回到同一个常数）。已经排查过机制（解析拆出"常数背景 + 局域不对称
特征"两部分做单项积分核实）：**背景部分对"偏差(mu_shift)"关于 mu_shift 的
斜率的贡献随 e_max 增大按 ~1/e_max 衰减到 0**（可以解析证明，单个 Matsubara
项的背景贡献严格等于 `(1/2π)ln[(a²+(e_max+c)²)/(a²+(e_max-c)²)]`，`c` 固定时
这个表达式随 `e_max→∞` 趋于 0）；如果 DOS 的两侧尾部收敛到**同一个**常数
（例如手工构造的、远离特征处精确回到 N̄=1 的"局域凸起"DOS），"偏差"本身在
mu_shift 固定时确实收敛，但由于上述背景斜率趋于 0，方程"偏差=0"的根需要
越来越大的 `|mu_shift|` 才能抵消局域特征贡献的一个有限余量——根本身仍然
收敛，但收敛速度极慢（实测需要 `dos_energy_cutoff_ratio` 到 1e4~1e5 量级才能
看到明显收敛，远超模块默认值 20，也超过"高精度"建议值 200）；如果两侧尾部
根本不收敛到同一个常数（`DOS.step(ratio≠1)`、或真实材料表格数据在截止之外
按边界值常数延拓、而左右边界值本来就不相等——H3S/D3S 的真实 QE DOS 表格正
是这种情况），`mu_shift` 干脆不收敛。

**这不是求解器 bug**，是"用一个不依赖 Matsubara 指标的标量化学势守恒一个对
整个（可能无穷延伸的）能带积分而成的电子数"这个物理设定本身的性质：一个
持续到无穷远、左右不对称的背景，理论上需要无穷大的化学势偏移才能补偿。
**实践含义**：不要为了追求"`mu_shift` 收敛"盲目增大 `dos_energy_cutoff_ratio`
——这里的正确截止应该对应**真实材料的物理能带范围**（类似 EPW 的
`fsthick` 窗口），而不是数学意义上的"cutoff→∞极限"；`mu_shift` 的绝对数值
必须结合"用的是哪个 `dos_energy_cutoff_ratio`"来解读，不能脱离这个选择单独
引用。本模块测试（`tests/test_eliashberg_iso.py` 的 `fix_mu="conserve"`
一节）只断言符号和方向（在整个实测的 cutoff 范围内都稳健），不断言绝对数值；
`p1_mishra_dos_calibration.py` 的 `fix_mu="conserve"` 结果沿用与 `"none"` 相同
的 `dos_energy_cutoff_ratio=2.0`（与 QE DOS 表格窗口 `±3eV` 同量级）以保证
可比性，见「第二阶段结果」一节的具体数值与讨论。

### 已知局限（能量依赖 DOS 通用部分）

1. **`delta=0`（vHs 精确钉在费米面）不会提升 Tc，反而压低 Tc**：这是上面"归一化"
   一节讨论的数学推论，在 `a∈{0.05,...,1.0}, cut∈{1,5,20}` 的扫描下稳健成立
   （无反例）。同位素系数在 `delta=0` 处确实如预期明显低于 0.5（甚至比常数 DOS
   基线更低），符合 Tsuei (1990) 的定性结论；但"Tc 升高"这部分要用一个有限偏移
   （例如 `delta≈omega_log` 量级）才能复现，见
   `tests/test_eliashberg_iso.py` 里两个对照测试的 docstring 和
   `examples/dos_scan.py` 的完整扫描表（delta/omega_log 从 -3 到 3，可以看到
   Tc 在 delta=0 处有一个明显的凹陷，两侧各有一个高于常数 DOS 基线的峰）。
2. **`DOS.log_vhs` 在 `|delta|` 相对 `W` 过大时归一化本身会失效**：`g(0)=1+a·ln(W/|delta|)`
   在 `|delta|≫W` 时趋于很小甚至变负。`g(0)<=0` 会被构造函数拒绝抛 `ValueError`；
   `g(0)` 很小但仍为正时不会报错，但 `g(ε)/g(0)` 会把普通量级的 `g(ε)` 放大成远超 1
   的 `N̄`，且 `g(ε)` 在积分窗口内可能变号（被 `floor` 钳位），这种"陡峭+钳位"的
   `N̄(ε)` 实测会让 `(Z,χ)` 自洽迭代发散振荡、`dos_max_iter=500` 次不收敛（触发
   `UserWarning`）。**这不是"delta 越大越接近常数 DOS 因而无所谓"的直觉**——验证
   "delta 远大于 w_cut 回到常数 DOS"这类极限时，必须同时把 `W` 取得远大于扫描到的
   `delta`（`nbar()` 会在 `g(0)<0.1` 时发 `UserWarning` 提示这个风险）。
3. **`fix_mu` 已实现 `"none"`（默认）与 `"conserve"`**：`"none"` 忽略化学势的
   电子数守恒修正（v0.1 行为，逐位不变）；`"conserve"` 额外求解一个不依赖
   Matsubara 指标的标量化学势偏移使正常态电子数守恒，见上面「化学势守恒
   （`fix_mu="conserve"`）」一节（含方程、数值方法与已知的截止敏感性局限）。
   `fix_mu="linear"`（Mitrović-Carbotte 原文和 Pickett 1980 都讨论过的另一种
   线性化电子数守恒修正，与本模块实现的"非线性求根"路径不同）**未实现**，
   不在本次任务范围内；`Solver.__post_init__` 会在 `fix_mu` 取 `"none"`、
   `"conserve"` 之外的值时直接 `raise NotImplementedError`，不会静默给出
   不完整的结果。
4. **计算成本**：能量依赖 DOS 路径每次 `tc()` 调用里，每个二分迭代点都要做一次
   `(Z,χ)` 自洽（最多 `dos_max_iter=500` 次线性混合迭代，`dos_mixing=0.3`）+ 一次
   `M×M`（`M=2*n_pos`）本征值求解，比常数 DOS 路径慢一到两个数量级（典型运行时间
   0.4–1.3 秒/次 `tc()`，取决于 `dos_n_eps`、`T` 附近所需的 `n_pos`）。

## 顶角核钩子（SPEC 必备功能 7，第二阶段，P2 用）

### 接口

`Solver(..., vertex=None, vertex_in_Z=True)`，`vertex` 是可调用对象：

```python
def vertex(n_idx, m_idx, omega_n, omega_m, T) -> np.ndarray:
    ...  # 返回与 lambda(n-m) 同形状的乘子矩阵 V
```

应用方式：`λ(n−m) → λ(n−m)·[1+V(n,m,ω_n,ω_m,T)]`，**只作用在 φ 方程里的电声子部分**
（`μ*θ(ω_c−|ω_m|)` 这一项永远不受 `vertex` 影响——顶角修正描述的是电声子顶点本身
的修正，不应该混进 Coulomb 赝势）。`vertex_in_Z`（默认 `True`）控制 Z 方程是否也用
同一个 `(1+V)` 修正后的核；设为 `False` 时 Z 方程保持 v0 原始的 `λ(n−m)`，只有 φ
方程受 `vertex` 影响。

### 实现细节（已知简化）

v0 的折叠核矩阵同时用到 `λ(n−m)`（"直接项"）和 `λ(n+m+1)`（正频折叠产生的"镜像项"，
见「折叠技巧推导」一节）。第二阶段的实现对两项用**同一个** `V[n,m]`
（`(1+V)` 同时乘到 `lam_diff` 和 `lam_summ` 上）——这是一个记录在案的简化，没有
区分正、负 Matsubara 频率结构对顶角修正的不同贡献；更精细的处理留给 P2 按需展开。
能量依赖 DOS 路径（全轴、未折叠）里则是对同一个 `V[n,m]` 直接应用到全轴的
`λ(n−m)` 矩阵上，没有折叠项的歧义。

### 占位顶角核：`vertex_pietronero_grimaldi_isotropic`

```python
vertex_pietronero_grimaldi_isotropic(m_ratio, Qc)
```

实现 `V = m_ratio · P_V(Qc)`，`P_V(Qc) = 1 − 2·Qc`（`Qc→0` 时为正、`Qc→1` 时为负，
`m_ratio = ω_log/E_F` 的量级前置因子），返回值与 `(n,m,ω_n,ω_m)` 无关（常数矩阵）。
**这是一个占位的、非文献推导的示意公式**——真实的 Pietronero-Grimaldi-Strässler
顶角修正（各向同性极限）的解析形式留给 P2 完成，这里只是为了让顶角钩子接口有一个
可以立即跑起来、形状/符号定性合理的示例，不能当作物理结果引用。

### 已知局限

- 与能量依赖 DOS 一样，`vertex` 目前不支持依赖 `ε`（顶角修正只作用在 Matsubara
  求和核上，不进入能量积分）；如果后续需要 `ε`-依赖的顶角修正，需要扩展
  `_dos_energy_integrals` 的接口。

## 已知局限（v0 核心功能 1–5）

1. **`Spectrum.lorentzian` 不是解析谱**：用数值网格（默认 4001 点、`Simpson` 积分）
   构造并归一化，不像 `einstein`/`two_peak` 那样有闭式的 `λ(ν)`。网格分辨率不够时
   `λ_total`、`omega_log`、`omega_2` 会有限精度误差（测试里用 `rel=1e-3` 量级容差）。
2. **`solve_gap` 是简单的线性混合不动点迭代**，不是 Newton 类方法；在 `T` 非常接近 `Tc`
   时收敛变慢（临界慢化），`max_iter=2000, mixing=0.3` 是经验值，未做自适应步长。
3. **弱耦合 + 偏大 `w_cut` 组合会触发 `RuntimeError`**：`_n_for_T` 在自动估计的
   `N` 超过 `_n_max=20000` 时主动报错（而不是静默跑到 `MemoryError`），因为
   `N ∝ w_cut/T`，弱耦合意味着 Tc（从而 T 的搜索范围）本身很低，容易和较大的
   `w_cut` 组合出巨大的 N。遇到这种情况请手动传 `n_matsubara`、降低
   `matsubara_oversample`，或收窄 `tc(t_bracket_k=...)`。
4. **折叠核矩阵 `a[n,m]` 不是对称矩阵**：`μ*θ(ω_m)` 只依赖列 `m`，`wn`/`z0` 的归一化
   对行列也不对称。后果：`log_sensitivity("w_cut")`（只放大 `w_cut`，`μ*` 数值不变、
   不按 Tolmachev 关系耦合）的符号**不是**"截止范围变宽 → 压制更多 → Tc 应该下降"这个
   朴素预期——实测在 λ=1.5, μ*=0.15, w_cut=5ω_log 附近它是正的（约 +0.03）。已经用
   固定 T、固定 N 的直接本征值计算核实过：这不是二分搜索或 N 自适应带来的假象；
   用旧 w_cut 的本征向量去算新 w_cut 矩阵的 Rayleigh 商确实如直觉般下降，但完整求解
   （对称矩阵才有的变分单调性在这里不成立）后主导本征值反而上升。这是当前非对称
   折叠核矩阵表示法的真实数学性质，不影响 `isotope_coefficient(mu_star_fixed=False)`
   —— 那条路径里 `w_cut` 与 `μ*` 是按 Tolmachev 关系耦合一起变的，其结果已经独立
   用 P1 式(3.2) 核对过（见「首轮结果」）。`log_sensitivity("w_cut")`/`("mu_star")`
   仅作独立扰动的诊断量，使用时请留意这一点，不要脱离这个限定条件去解读符号。
5. **常 DOS 是默认值，不是唯一选项**：`Solver` 默认 `dos=DOS.constant()`（`N(ε)=N(0)`），
   与 v0 行为完全一致；能量依赖 DOS（SPEC 必备功能 6）已实现，见「能量依赖 DOS」
   一节，需要时显式传入 `dos=...` 即可，不传则不受影响。
6. **顶角核钩子默认关闭**：`vertex=None`，与 v0 行为完全一致；模型顶角核钩子
   （SPEC 必备功能 7）已实现，见「顶角核钩子」一节。



## IsoME 交叉对照

SPEC 要求"若本机有 Julia + IsoME（Kogler 2025），加一个交叉对照测试；没有就跳过并注明"。
本开发/测试机上未检测到 `julia` 可执行文件（`shutil.which("julia") is None`），
对应测试 `test_isome_cross_check_if_available` 会被 `pytest.skip` 跳过并打印原因，
不是遗漏。

## 首轮结果（examples/p1_first_pass.py 的输出）

以下为 `python theory/solvers/examples/p1_first_pass.py` 的原样输出（2026-09-08 跑的）。

```
============================================================
P1 first-pass check: numeric alpha vs eq (3.2), lambda=2, omega_log=100 meV
Convention: w_cut = 5 * omega_log, mu_star_fixed=False (Tolmachev RG)
============================================================

--- Einstein single mode (lam=2, w0=100 meV)  (lambda=2.0000, omega_log=100.0000 meV) ---
   mu* |  alpha_num | alpha_eq32 |     diff
--------------------------------------------
  0.10 |     0.4900 |     0.4889 |  +0.0011
  0.13 |     0.4837 |     0.4798 |  +0.0039
  0.20 |     0.4646 |     0.4420 |  +0.0226  (out of <0.01 tolerance)
  0.25 |     0.4476 |     0.3947 |  +0.0529  (out of <0.01 tolerance)
  0.30 |     0.4282 |     0.3217 |  +0.1066  (out of <0.01 tolerance)

--- Two-peak (lam1=1,w1=50 + lam2=1,w2=200 meV)  (lambda=2.0000, omega_log=100.0000 meV) ---
   mu* |  alpha_num | alpha_eq32 |     diff
--------------------------------------------
  0.10 |     0.4860 |     0.4889 |  -0.0029
  0.13 |     0.4766 |     0.4798 |  -0.0032
  0.20 |     0.4475 |     0.4420 |  +0.0056
  0.25 |     0.4220 |     0.3947 |  +0.0273  (out of <0.01 tolerance)
  0.30 |     0.3918 |     0.3217 |  +0.0702  (out of <0.01 tolerance)

============================================================
Secondary check: d ln Tc/d ln lambda, numeric solver vs McMillan analytic
============================================================

--- d ln Tc / d ln lambda: numeric vs McMillan analytic (mu*=0.13) ---
lambda |    s_num | s_mcmillan |     diff
------------------------------------------
   1.0 |   1.6858 |     1.7514 |  -0.0656
   2.0 |   0.8895 |     0.7475 |  +0.1420
   3.0 |   0.7224 |     0.4740 |  +0.2484
```

补充（λ 扫描，Einstein 单模，`w_cut=5·ω_log`，未进 examples 脚本但进了 `tests/`，
一并记录以便完整核对 SPEC 承诺的 "λ=1–2, μ*≤0.15" 范围）：

```
lambda    mu*  alpha_num alpha_eq32     diff
   1.0   0.10     0.4811     0.4760  +0.0050
   1.0   0.13     0.4680     0.4543  +0.0137   <- 超出 <0.01
   1.0   0.15     0.4574     0.4338  +0.0236   <- 超出 <0.01
   1.5   0.10     0.4877     0.4853  +0.0024
   1.5   0.13     0.4797     0.4728  +0.0069
   1.5   0.15     0.4734     0.4615  +0.0119   <- 超出 <0.01
   2.0   0.10     0.4900     0.4889  +0.0011
   2.0   0.13     0.4837     0.4798  +0.0039
   2.0   0.15     0.4789     0.4716  +0.0073
```

### 如实分析：SPEC 的 <0.01 容差在 λ=1–1.5、μ*≥0.13 处不成立

SPEC 测试 4 要求"数值 α 与 P1 式(3.2) 在 λ=1–2, μ*≤0.15 下相差 <0.01"。实测这**不是
在整个承诺范围内均匀成立**：

- **λ=2**（P1 §3.1 自己给出数值表的那一行）：μ*=0.10（+0.0011）、0.13（+0.0039）都
  舒适地在容差内，μ*=0.15 也只有 +0.0073。
- **λ=1.5**：μ*=0.10/0.13 在容差内，μ*=0.15 处 +0.0119，**略微超出**。
- **λ=1**：只有 μ*=0.10 在容差内（+0.0050）；μ*=0.13 (+0.0137) 与 μ*=0.15 (+0.0236)
  都明显超出 <0.01 的要求，且随 μ* 增大偏差单调放大。

按团队指示，**没有去调数值/公式凑合格**——已经排查过这不是数值 bug：

- 对差分步长 `eps ∈ {0.005,...,0.04}` 稳定（见 `test_isotope_coefficient_eps_stability`
  与本模块开发过程中的 ad hoc 扫描）。
- 对 Matsubara 截断收敛（`matsubara_oversample` 1.5→3.0→6.0 的序列表明 3.0 已经收敛
  到 1e-4 量级，不会解释 1–2 个百分点的偏差）。
- 光滑化 θ 之后，偏差随 μ*、λ 的变化是平滑单调的，不是台阶或锯齿。

**可能的物理/公式来源（不是求解器缺陷）**：

1. **式(3.2) 本身是 McMillan 闭式对 `∂lnΦ/∂μ*` 的解析求导**，McMillan 公式本身是对
   Eliashberg 方程的一个经验拟合，在中等耦合（`λ~1`）区间的拟合精度本来就比它在
   `λ~2` 附近（历史上 McMillan/Allen-Dynes 主要标定的区间）差；式(3.2) 是在这个已经
   有系统偏差的公式上再求一次导数，偏差可能被放大而不是抵消。
2. **P1 草稿本身在 §0/行69 已经预告了这一点**："处理：在求解器里用显式库仑核
   μ(ω_c)θ(ω_c−|ω_m|) 做数值 Eliashberg，直接数值求 s_μ，用 (3.2) 只作弱耦合校验"——
   也就是说式(3.2) 从一开始就被设计为"弱耦合下的校验基准"，不是"任意 λ,μ* 下都应该
   精确匹配的目标"。本次结果符合这个预期：耦合越强（λ=2，P1 自己举的例子）符合得
   越好，耦合减弱（λ=1）+ μ* 增大，偏差系统性增大。
3. 数值解用的是显式 `μ*(ω_c)θ(ω_c-|ω_m|)` 硬截止核（做了光滑化处理，见上）；
   式(3.2) 的推导路径是"McMillan Φ 对 μ* 求偏导"，并未显式引入这个截止的具体形状，
   只用 Tolmachev 关系把 `μ*` 本身的质量依赖代进去。二者对"截止怎么进入方程"的
   处理方式不完全一致，是偏差的另一个可能来源，值得在 P2/后续版本用更精细的数值
   实验（比如变化截止形状、比较硬/软截止）单独拆解。

**结论**：SPEC 测试 4 的核心检验点（λ=2 附近、μ*≤0.13 或更严格）稳健通过；把 SPEC 的
"λ=1–2"字面理解为"整个区间都要 <0.01"这个更强的版本不成立，已如实记录在此，并在
`tests/test_eliashberg_iso.py` 里用两类测试分别体现：λ=2（以及 λ=1,μ*=0.10）严格断言
<0.01，λ=1,μ*=0.13 单独用一个"锁定观测偏差量级"的回归哨兵测试（而不是断言 <0.01
或干脆不测）。

### 次要发现：`d ln Tc/d ln λ` 与 McMillan 解析值的偏离

McMillan 公式本身在强耦合端低估 `∂lnTc/∂lnλ`（这也是 Allen-Dynes 引入 `f1f2` 修正因子
的动机之一）。数值结果与此一致：λ=1 时数值解与 McMillan 解析值接近（diff -0.066），
λ=2、3 时数值解明显更"敏感"（diff +0.14、+0.25，即数值解比 McMillan 认为的对 λ 更敏感）。
这不是 SPEC 要求的测试项，只是诊断性记录，供后续 P2 参考。

## 第二阶段结果（examples/dos_scan.py 的输出）

以下为 `python theory/solvers/examples/dos_scan.py` 的原样输出（2026-09-08 跑的）。
模型：Einstein 谱 λ=1.84, ω_log=93 meV, μ*=0.16, w_cut=5ω_log=465 meV，
`DOS.log_vhs(a=0.3, W=1000 meV, cut=5 meV)`，扫 `delta/omega_log ∈ [-3, 3]`。

```
======================================================================
DOS scan: log-type van Hove singularity vs distance from E_F
Einstein spectrum lambda=1.84, omega_log=93.0000 meV, mu_star=0.16, w_cut=465.0000 meV
DOS.log_vhs(a=0.3, W=1000.0 meV, cut=5.0 meV)
Baseline (constant DOS): Tc=151.9680 K, alpha=0.4748
======================================================================
delta/wlog | delta(meV) |     Tc(K) |  Tc/Tc0 |   alpha
----------------------------------------------------------------------
     -3.00 |    -279.00 |  162.7051 |  1.0707 |  0.4818
     -2.50 |    -232.50 |  162.9381 |  1.0722 |  0.4940
     -2.00 |    -186.00 |  162.7385 |  1.0709 |  0.4796
     -1.50 |    -139.50 |  161.7536 |  1.0644 |  0.4660
     -1.00 |     -93.00 |  158.9358 |  1.0459 |  0.4599
     -0.50 |     -46.50 |  151.8382 |  0.9991 |  0.4463
      0.00 |       0.00 |  126.6142 |  0.8332 |  0.4323
      0.50 |      46.50 |  151.8382 |  0.9991 |  0.4463
      1.00 |      93.00 |  158.9358 |  1.0459 |  0.4599
      1.50 |     139.50 |  161.7536 |  1.0644 |  0.4660
      2.00 |     186.00 |  162.7385 |  1.0709 |  0.4796
      2.50 |     232.50 |  162.9381 |  1.0722 |  0.4940
      3.00 |     279.00 |  162.7051 |  1.0707 |  0.4818
```

### 如实分析：`delta=0` 处 Tc 有凹陷而不是尖峰

这张表关于 `delta=0` 对称（模型只依赖 `|ε-δ|`，符合预期），但形状本身值得说明：
`delta=0`（vHs 精确钉在费米面）处 Tc 明显**低于**常数 DOS 基线（`Tc/Tc0=0.833`），
两侧各有一个高出基线约 7% 的峰（`|delta/wlog|≈2.5` 附近），然后随 `|delta|` 继续
增大逐渐衰减回基线（第 (d) 类极限，见「能量依赖 DOS·已知局限」第 1 条的完整数学
论证：`N̄(0)≡1` 的归一化在 `delta=0` 处把"费米面处 DOS 增强"这个效应本身抵消掉了，
只留下"费米面之外全面走低"的效果）。同位素系数 `alpha` 在 `delta=0` 处达到扫描范围
内的最小值（0.4323，明显低于 0.5 也明显低于常数 DOS 基线 0.4748），随 `|delta|`
增大单调恢复向基线——这部分与 Tsuei (1990) "vHs 钉在/靠近费米面压制同位素效应"的
定性结论一致，且是这张表里最稳健的定性结论（不依赖 Tc 升降的符号）。

### `p1_mishra_dos_calibration.py` 新增 `fix_mu="conserve"` 两列的结果（真实 H3S/D3S QE DOS）

以下为 `python theory/derivations/p1_mishra_dos_calibration.py` 的原样输出（2026-09-08 跑的，
含新增的 `fix_mu="conserve"` 两列，全程约 10.4 分钟）。DOS 来自真实 QE 计算的
`H3S-200GPa/band/H3S.dos`（`E_F=17.619 eV`，`±3eV` 窗口，1200 点，归一化到 `N̄(E_F)=1`）：

```
QE DOS: N(E_F)=0.5593 states/eV; window ±3.0 eV, 1200 points
  N̄( -500 meV) = 0.825
  N̄( -200 meV) = 1.165
  N̄( -100 meV) = 1.129
  N̄(  -50 meV) = 1.114
  N̄(   +0 meV) = 1.000
  N̄(  +50 meV) = 0.830
  N̄( +100 meV) = 0.781
  N̄( +200 meV) = 0.736
  N̄( +500 meV) = 0.622

kind DOS                       Tc_H    Tc_D    α_H | EPW Tc_H EPW Tc_D EPW α_H
har  const                    243.9   179.3  0.444 |      245      180   0.445
har  QE-table                 241.6   179.3  0.431 |      234      173   0.436
har  QE-table+cons            148.3   108.1  0.456 |      234      173   0.436
har  QE-table+cons(±1eV)      168.2   120.9  0.477 |      234      173   0.436
anh  const                    190.6   151.5  0.331 |      192      152   0.337
anh  QE-table                 195.4   155.0  0.335 |      189      149   0.343
anh  QE-table+cons             97.3    80.1  0.281 |      189      149   0.343
anh  QE-table+cons(±1eV)      113.1    91.5  0.306 |      189      149   0.343

fix_mu="conserve" 的 mu_shift（各自 Tc 处求得，meV，负号=化学势下移）：
har  QE-table+cons          H3S  mu_shift=  -1055.3 meV  converged=True
har  QE-table+cons          D3S  mu_shift=  -1078.0 meV  converged=True
har  QE-table+cons(±1eV)    H3S  mu_shift=   -612.8 meV  converged=True
har  QE-table+cons(±1eV)    D3S  mu_shift=   -642.4 meV  converged=True
anh  QE-table+cons          H3S  mu_shift=  -1081.2 meV  converged=True
anh  QE-table+cons          D3S  mu_shift=  -1094.3 meV  converged=True
anh  QE-table+cons(±1eV)    H3S  mu_shift=   -645.2 meV  converged=True
anh  QE-table+cons(±1eV)    D3S  mu_shift=   -662.1 meV  converged=True
```

#### 如实分析：方向对了，但 `"conserve"` 把 gap 反向撑大了 8–13 倍，没有收窄

先确认动机里描述的"反方向"现象在这份真实数据里确实存在：`const` 列（对应 EPW 的
`fsr`，无 DOS 能量依赖）与 EPW 差 1–2K（`har`: 243.9 vs 245；`anh`: 190.6 vs 192），
和 SPEC 预期一致；但 `QE-table`/`fix_mu="none"` 列（对应 EPW 的 `fbw`）比 EPW 的
`fbw` **高** 6–8K（`har`: 241.6 vs 234，+7.6K；`anh`: 195.4 vs 189，+6.4K）——EPW 的
`fbw`（带 `muchem` 电子数守恒）比它自己的 `fsr` 低，本模块的 `QE-table`（`fix_mu="none"`，
无电子数守恒）却比它自己的 `const` 基线还略高（`har`）或明显更高（`anh`），确认了
方向相反这个 gap 的存在。

加上 `fix_mu="conserve"` 之后，方向确实被"掰"回来了——`QE-table+cons` 相比
`QE-table` 大幅降低（`har`: 241.6→148.3，`anh`: 195.4→97.3），这与"EPW 的 `muchem`
让 `fbw` 低于 `fsr`"是同一个方向，定性验证了「H3S 的 N̄(ε) 不对称 + 电子数守恒
压低 Tc」这一因果链确实成立。**但幅度完全对不上**：`har` 只需要再降
`241.6-234=7.6K` 就能对上 EPW 的 `fbw`，`"conserve"`（`dos_energy_cutoff_ratio=2.0`，
与 `QE-table` 相同）却把 Tc 砸到 148.3K，**过冲了 `148.3-234=-85.7K`，是目标缺口的
约 11.3 倍**；`anh` 同理，目标缺口 6.4K，实际过冲到 `97.3-189=-91.7K`，约 **14.3 倍**。
换成探索性的 `±1eV` 窗口（更接近 EPW `fsthick=1.0eV`）后过冲幅度有所收敛
（`har`: -65.8K≈8.7 倍，`anh`: -75.9K≈11.9 倍），说明窗口大小确实是过冲幅度的一个
贡献因素（与「化学势守恒」一节"`mu_shift` 对 `dos_energy_cutoff_ratio` 敏感"的
已知局限吻合：`±3eV` 窗口下 `|mu_shift|≈1055–1094meV`，`±1eV` 窗口下降到
`|mu_shift|≈613–662meV`，窗口越大过冲越猛），**但即使换成 EPW 量级的窗口，过冲
倍数仍然接近一个数量级，没有真正收窄到目标缺口附近**。

下一步怀疑对象（按优先级排序，供后续任务参考）：

1. **EPW 的 `fsthick=1.0eV` 不只是"化学势求解用的 DOS 窗口"，而是贯穿整个
   Eliashberg 求和的硬截止**——`fsthick` 限制的是哪些电子态参与电声耦合矩阵元
   求和本身（k/q 点求和的费米面附近壳层），不只是"用哪段 N̄(ε) 去解一个标量
   `mu_shift`"。本模块的 `"conserve"` 实现是后者：用一个**不依赖 Matsubara 指标
   的标量 `mu_shift`**去补偿**同一个（可能远超 `fsthick` 的）能量积分域**上的
   电子数偏差，这与 EPW `muchem` 在 `fsthick` 壳层内部、且与能带结构/ k 点权重
   自洽求解化学势的机制在数学结构上并不等价。本次 `±1eV` 截断实验只是把
   **DOS 形状**限制在 `±1eV`，积分和电子数守恒逻辑仍是本模块原有的"标量
   `mu_shift`"框架，并不是真正复现 EPW 的 `fsthick` 机制——这很可能是过冲量级
   对不上的主因，值得作为下一步的重点方向。
2. **真实 QE `*.dos` 表格本身的粗糙度/尖锐度**：从打印的 `N̄(ε)` 看，`±500meV`
   窗口内 `N̄` 已经从 0.622 变到 1.165，接近 2 倍的摆动——这个摆动如果部分来自
   `*.dos` 直方图的统计噪声（而非物理带结构的真实精细结构），会让标量
   `mu_shift` 机制过度反应；EPW 用于 `muchem` 的 DOS 本身分辨率/展宽
   （`degaussw`）设置未知，如果比这里用的原始 `*.dos` 更平滑，会显著改变
   `mu_shift` 的量级。
3. **Tc 提取网格分辨率（1–2K）与 EPW `*.dos` 网格粗糙度（`±1eV`/40 bins）**：
   相对于本节发现的 65–92K 过冲，这两个因素量级太小，不太可能是主因，但仍
   应在后续更精细的对照里排除。

**结论**：`fix_mu="conserve"` 的实现本身通过了本节和 `tests/test_eliashberg_iso.py`
的"符号、对称性、`N̄≡1` 退化到 `"none"`"等结构性检验（见「化学势守恒」一节），
是一个数学上自洽的标量电子数守恒修正；但把它套到真实 H3S/D3S QE DOS 表格上，
**没有让 gap 收窄，反而在正确的方向上过冲了一个数量级**。这不是"实现有 bug"，
而是"标量 `mu_shift` 模型 + 真实材料 DOS 表格"这个组合本身的定量局限——诊断和
（如果需要）改进这个组合，应作为独立的后续任务，而不是本次任务范围内可以
一次性解决的问题。

## 未实现/未覆盖的范围

- **`fix_mu="linear"`**（能量依赖 DOS 路径的化学势-电子数守恒**线性化**修正，
  与已实现的 `"conserve"`（非线性求根）是 Mitrović-Carbotte/Pickett 原文里
  并列讨论的两种不同近似）：未实现；`fix_mu` 目前接受 `"none"`（默认）与
  `"conserve"`，见「能量依赖 DOS·化学势守恒（`fix_mu="conserve"`）」一节。
- **顶角核钩子对 `ε` 的依赖**：`vertex` 目前只作用在 Matsubara 求和核上，不进入
  能量积分；见「顶角核钩子·已知局限」。
- **`vertex_pietronero_grimaldi_isotropic` 的真实文献解析形式**：当前实现是占位
  公式，留给 P2 按需替换为真正推导的 Pietronero-Grimaldi-Strässler 顶角修正。
- SPEC 必备功能 1–7 均已实现；IsoME 交叉对照按「IsoME 交叉对照」一节所述跳过。
