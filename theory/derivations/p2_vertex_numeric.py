"""
theory/derivations/p2_vertex_numeric.py

一阶电声顶角修正 P(k+q,k) 的数值积分（自由电子传播子、Matsubara 求和已做的
形式），来源：Grimaldi, Pietronero, Scattoni, Eur. Phys. J. B 10, 247 (1999),
arXiv:cond-mat/9806012（下称 GPS99）。本文件按 theory/derivations/P2_numeric_SPEC.md
§1-§4 实现；所有公式已对照本地 PDF（library/pdf/P2/1999_C_Grimaldi_cond-mat-9806012.pdf）
第 3-6 页原文核对，Eq. 编号均指该文。

============================================================================
0. 与 SPEC.md 转录不一致之处（已核对原文，以原文/自洽推导为准，用的就是这里）
============================================================================

(a) **v_F q 与 Q 的关系**：SPEC.md 写 "Q≡q/(2k_F)，v_F q = 2E_F·Q，v_F=2E_F/k_F"，
    三者不自洽：由 v_F=2E_F/k_F（标准自由电子关系 E_F=k_F²/2m, v_F=k_F/m 的
    等价形式）与 Q=q/(2k_F) 直接代入，v_F q = (2E_F/k_F)(2k_F Q) = **4E_F·Q**，
    而非 2E_F·Q。本文件用 v_F q = 4 E_F Q（`_vFq`），并在线性化色散
    ε_{k'+q} ≈ ε' + v_F q·cosθ 中使用这个正确关系。用错的 2E_F·Q 会让 Q 的
    有效标度差一倍，四个解析极限的数值核对会直接暴露这个问题（已用两个版本
    交叉测试过，4E_F·Q 版本通过，2E_F·Q 版本静态/动态极限对不上）。

(b) **λ 的归一化**：SPEC.md 猜测 "λ=2N(0)g²/ω₀（原文 g²=2γ²/ω₀）"。这是标准
    McMillan-Hopfield 公式的形式，但不适用于 GPS99 的记号：GPS99 的 Eq.(6)/(8)/(10)
    里 ω₀ 已经作为**显式前因子**提出来（(ω₀/2)Σ...），耦合本身只以 g²_{k-k'}
    的形式出现，不含 1/ω₀。直接从 Eq.(11) 反推：取结构无关 g²、常数 DOS N(0)，
    对 P_pol(q=0,ω_q=0)（Eq. 10，单电子情形，k' 求和跑遍能带 0→E，无 Fermi
    函数限制）做出解析积分，要求结果等于 Eq.(11) 的 λ(E/2)/(ω₀+E)，得到
        **λ = N(0)·g²**（无 2、无 1/ω₀ 因子）。
    独立交叉验证：GPS99 Eq.(20) 给出 lim_{q→0,ω_q→0} P_ex = Σ_k' g² df/dε_k' ≈ -λ；
    T→0 时 df/dε=-δ(ε)，Σ_k'→∫dε N(ε)，故 Σ_k' g² df/dε = -g²N(0) = -λ，
    **同一关系**，两处独立公式互相印证。本文件采用 λ=N(0)g²；报告里给出
    Eq.(11) 数值核对作为证据。
    实现上不显式出现 g²：把传入的 `dos` 统一约定为**已按 N(0)（或 N(E_F)）
    归一化**的 N̄(ε)=N(ε)/N(0)（N̄(0)=1），并令 g²≡1。这样任何时候
        P_actual = λ · vertex_P(...)
    即 `vertex_P`/`vertex_PV` 的返回值本身就是 P/λ，与 SPEC §2 "最终 P 以 λ
    为单位输出" 的要求一致，不需要再单独除以 λ。

(c) **静态/动态极限的"先后次序"表述**：本文件早期版本曾根据对原文第 3 页倒数
    第二段英文的字面解读，认为 SPEC.md 的"静态：先 ω_q→0 再 q→0；动态：先 q→0
    再 ω_q→0"次序与原文相反，并按反过来的次序实现——**这是错的，已改正**。
    独立数值验证（对 Eq.(6) 直接积分，不依赖任何一侧的转录）明确显示（数值见
    下面"数值做法"给出的具体路径）：
        - "Q 固定为一个不太小也不太大的非零值、令 ω_q→0（内层）"这条路径的
          极限值是 **P_s = -ω₀/(ω₀+E_F)**；
        - "ω_q 固定为一个不太小也不太大的非零值、令 Q→0（内层）"这条路径的
          极限值是 **P_d = +E_F/(ω₀+E_F)**。
    这与 SPEC.md 的文字描述（先 ω_q→0 再 q→0 得 P_s；先 q→0 再 ω_q→0 得 P_d）
    **完全一致**——SPEC.md 是对的，本文件早期的"更正"才是错的。
    但如果把原文"P_s ≡ lim_{ω_q→0} P(q=0,ω_q)"里的记号"q=0"读成**真正的代入**
    （先把 Q 严格置零，再让 ω_q→0），数值上给出的却是 **P_d 的值**，不是 P_s！
    这不是巧合，是嵌套极限的定义本身决定的：把 Q 严格置零就是"最先、最彻底地"
    做了 Q→0 这一步，之后再让 ω_q→0，按定义就等于
    lim_{ω_q→0}[lim_{Q→0} P(Q,ω_q)]——这恰好是 P_d 的定义，不是 P_s 的定义。
    也就是说，原文自己给 P_s 写的公式"P(q=0,ω_q)"，如果按字面代入执行，算出来
    的是原文自己定义的 P_d；这是原文记号本身的自相矛盾（"q=0"在物理上想表达
    "取 q→0 的外层极限"，但写成公式字面意义变成了"先代入"，两者对嵌套极限而
    言不是一回事），不是本文件的转录错误。能分辨该按哪条路径走的，只有直接
    数值核实（下面给出）——这也解释了为什么本文件早期会把次序搞反：早期正是
    按"字面代入 q=0"的读法实现的。
    嵌套极限的准确写法（现在的实现和测试都按这个）：
        P_s = lim_{Q→0} [ lim_{ω_q→0} P(Q,ω_q) ]  （内层 ω_q→0，外层 Q→0）
        P_d = lim_{ω_q→0} [ lim_{Q→0} P(Q,ω_q) ]  （内层 Q→0，外层 ω_q→0）
    数值做法（tests/test_p2_vertex.py 用同样的取值）：
        - 静态 P_s：固定 Q=min(2·ω₀/E_F, 0.1)（"够小、非零"的一个折中——见下）
          、ω_q 用极小正则化 floor→0，看是否收敛到 -ω₀/(ω₀+E_F)；
        - 动态 P_d：固定 Q→0（用极小值，如 1e-6），ω_q=0.05ω₀（"够小、非零"），
          看是否收敛到 +E_F/(ω₀+E_F)。
    实测（EF=1, ω₀=0.1）：Q=0.1,ω_q→floor 给 -0.0915（P_s_ana=-0.0909，
    相对误差 0.7%）；Q→1e-6,ω_q=0.05ω₀ 给 +0.9083（P_d_ana=+0.9091，
    相对误差 0.09%）；而"严格 Q≡0，ω_q→floor"给 +0.9091（=P_d，印证上一段）。
    "Q=min(2ω₀/E_F,0.1)"不是随便调出来的：Q 必须远大于正则化 floor（否则数值
    上还没脱离 floor 的量级，见 §1 奇异性减除的讨论），又必须远小于"能带边缘"
    效应开始起作用的尺度（数值上发现约在 Q~0.3-0.5，与 E_F 同量级，同 ω₀/E_F
    无关）；两者中间有一个宽的安全窗口，min(2ω₀/E_F,0.1) 落在窗口正中，对三个
    ω₀/E_F∈{0.01,0.1,0.5} 都给出 <1% 的相对误差（见 tests 与报告的极限核对表）。

(c') **Fig.2 定性符号判据（v_Fq/ω_q 与正负号的对应）有同一种记号歧义，方向
    相反于原文字面表述**——根源与 (c) 完全相同，故记在这里而不是重开一条。
    原文与 SPEC.md 都写"v_Fq/ω_q>1 时 P 为正，<1 时为负"。用 (c) 已验证的嵌套
    极限去看这句话，会发现它与 P_s、P_d 的取值互相矛盾：
        - P_s 对应"固定非零 Q，令 ω_q→0（内层）"这条路径，这条路径上 v_Fq 不
          变、ω_q→0，所以 v_Fq/ω_q → **+∞**（远大于 1）；P_s 是**负**的。
        - P_d 对应"固定非零 ω_q，令 Q→0（内层，v_Fq→0）"这条路径，
          v_Fq/ω_q → **0**（远小于 1）；P_d 是**正**的。
    即按 (c) 验证过的、真正对应 P_s/P_d 定义的极限路径，v_Fq/ω_q≫1 的一端是
    **负**、≪1 的一端是**正**——与原文 Fig.2 那句话的方向刚好相反。
    直接数值验证（固定 Q=0.2、E_F=1、w0=0.1，扫 ω_q 从 5w0 到 0.001w0，
    v_Fq=4E_F Q=0.8 不变）：ω_q 越小（v_Fq/ω_q 越大）P 越负（v_Fq/ω_q=8000 时
    P≈-0.097，趋于 P_s）；ω_q 越大（v_Fq/ω_q 越小）P 越正（v_Fq/ω_q≈1.6 时
    P≈+0.083，趋于 P_d 方向）；穿越零点发生在 v_Fq/ω_q≈10-16（不是恰好 1，
    这只是本例参数下的具体数字，不是普适 claim）。整条曲线单调、连续，不是只
    在两个极限点上矛盾——不是极限点的孤立巧合，而是 v_Fq/ω_q 这个比值在原文
    里被反着标了方向。结论：原文 Fig.2 一句话的符号判据，用 (c) 已验证的嵌套
    极限约定去读，方向是**反的**；正确方向是 v_Fq/ω_q≫1 时 P<0（静态型），
    v_Fq/ω_q≪1 时 P>0（动态型）。这与 (c) 是同一件事的两个侧面：原文
    "P_s≡lim_{ω_q→0}P(q=0,ω_q)"和"Fig.2：v_Fq/ω_q>1→正"两处表述，都用了同一
    个不严谨的"q=0 即代入"读法，而这个读法本身对嵌套极限而言是错的（(c) 已
    独立数值证明）；一旦按验证过的嵌套极限改正，两处不一致同时自动改正，不
    需要分别调参数或分别改公式。SPEC.md 转录的是原文这句话，所以 SPEC.md 在
    这一点上也和原文一样方向反了；已在报告里向 team-lead 说明。

(d) P_s、P_d、P_pol 公式在原文里都没有显式写 λ 或 g²，但 P_pol 的 Eq.(11) 明
    写 "=λ(E/2)/(ω₀+E)"，三者应是同一意义下的量（都是"以 λ 为单位"）——
    否则量纲/线性关系对不上（P 对 g² 严格线性，Ps、Pd 不可能不含耦合强度）。
    本文件按"P_s、P_d 也是 P_s/λ、P_d/λ"处理，与 (b) 的 g²≡1、λ≡1 约定完全
    对应，不需要额外改动数值。

============================================================================
1. 公式与积分约定
============================================================================

核心公式（GPS99 Eq. 6，Matsubara 求和已做，T→0，外线电子频率 iω_k 换成小量
正则化 η_k>0 代表 "ω_k=πT, T→0"，见下）：

    P(Q,ω_q) = (ω₀/2) Σ_k' g²_{k-k'} · 1/(ε_{k'} - ε_{k'+q} + iω_q)
      × [ (f(ε_{k'}) + n(-ω₀)) / (ε_{k'} + ω₀ - iη_k)
        − (f(ε_{k'}) + n(ω₀))  / (ε_{k'} − ω₀ − iη_k)
        − (f(ε_{k'+q}) + n(-ω₀)) / (ε_{k'+q} + ω₀ − i(η_k+ω_q))
        + (f(ε_{k'+q}) + n(ω₀))  / (ε_{k'+q} − ω₀ − i(η_k+ω_q)) ]

k' 求和 → 各向同性能量-角度积分：

    Σ_k' A(k') → ∫ dε' N̄(ε') · (1/2)∫_{-1}^{1} d(cosθ) A(ε', cosθ)

N̄ 已按 N(E_F)=1 归一化（见 §0(b)）；ε' 从 band_bottom 到 band_top（相对 E_F）。

色散：默认线性化 ε_{k'+q} = ε' + v_F q·cosθ，v_F q = 4 E_F Q（§0(a)）。
可选精确自由电子色散（`exact_dispersion=True`）：
    ε_{k'+q} = ε' + 4E_F·Q·[Q + cosθ·sqrt(1+ε'/E_F)]
（由 E=k²/2m, q=2k_F Q 代入 Δε=(q²+2k'q cosθ)/2m 得到；要求 ε'>-E_F，
只在常数-DOS 自由电子模型下使用，H3S 真实能带远深于任何"等效自由电子 E_F"，
不适用，§2 已在报告里说明）。

T=0：f(ε)=θ(-ε)，n(ω₀)=0，n(-ω₀)=-1（玻色分布 n(-ω₀)=-1-n(ω₀) 的 T→0 特例）。
T>0：用真实 Fermi/Bose 分布（供以后扩展用，本任务的四个解析极限都在 T=0 下核对）。

正则化：
    η_k（外线电子频率，"ω_k=0 或 πT"）：默认 η_k = REG_FRAC · ω₀。
    ω_q floor（网格里 ω_q=0 这一行本身在 Eq.(6) 的第一个分母
      1/(ε_{k'}-ε_{k'+q}+iω_q) 处于 (Q=0,ω_q=0) 才真正奇异；ω_q=0 但 Q≠0
      时该分母是 -v_Fq·cosθ，仅在 cosθ=0 处主值奇异）：ω_q=0 时用
      ω_q_eff = max(ω_q, REG_FRAC·ω₀) 代替，即用一个小但非零的 ω_q 代表
      "ω_q→0" 这一极限点。
    REG_FRAC 默认 0.03（即 η_k = ω_q_floor = 0.03 ω₀），用于 A-E 的主网格
    （物理上代表有限小 T 的合理量级，网格里 ω_q 本身多数取非零值，不需要
    push 到极限）。四个解析极限核对（§0(c)）需要更干净地分离"ω_q→0"与
    "Q→0"两个尺度，用更小的 REG_FRAC_LIMIT=1e-4（即 η_k=ω_q_floor=1e-4 ω₀）
    ，详见 tests/test_p2_vertex.py；收敛性（REG_FRAC_LIMIT 取 1e-3～1e-5
    结果几乎不变）见该文件的 test_regularization_convergence。

数值积分中的奇异性减除（重要，direct 影响精度）：cosθ 积分的被积函数在
denom0=-v_Fq·cosθ+iω_q_eff 里有一个位于 cosθ=iω_q_eff/v_Fq 附近的近极点；
当 ω_q（或 floor）远小于 v_Fq 时，这个近极点在实轴 cosθ∈[-1,1] 上投影出的
"宽度"~ω_q_eff/v_Fq 可以远小于固定节点数 Gauss-Legendre 求积能分辨的尺度，
此时增加节点数也收敛得很慢——直接后果是：把 Eq.(6) 在 Q→0、ω_q→0 附近数值
积分时，结果会系统性地"看起来收敛了，但收敛到错误的数"（这正是本文件早期
版本里 P_s、P_d 解析极限对不上的根源，见下面的调试记录）。修复：把
cosθ=0 处的被积函数值（bracket_0）解析地从 denom0 的奇异部分中分离出来——
    bracket(cosθ)/denom0 = bracket_0/denom0 + [bracket(cosθ)-bracket_0]/denom0
第一项对 cosθ 的积分是初等函数（对 1/(-v_Fq x+iω_q_eff) 积分给出 log 的
闭式），不依赖离散化；第二项在 cosθ=0 处可去奇点（分子分母同阶趋零），
用普通固定 Gauss-Legendre（`N_COSTHETA` 个节点）求积即可稳定收敛。
`vertex_P` 内部按这个方式实现，不需要用户额外设置。

**调试记录**（供复核者参考，不是当前代码状态）：本文件最初直接对
bracket/denom0 做固定 257 点 Gauss-Legendre 求积，四个解析极限里 P_s、P_d
对不上（数量级、符号都不对）。逐步排查：(1) 重新核对 v_Fq=4E_F·Q、
λ=N(0)g² 两处转录更正，Ppol 极限已经能精确匹配（验证了这两处更正是对的）；
(2) 用与本模块完全独立、从零手写的脚本重新数值积分 Eq.(6)，发现 P_s、P_d
的"先后次序"数值上与 SPEC.md/原文标注刚好对应（见 §0(c) 的更正历史）；
(3) 排查 P_s 在小 Q 处随正则化尺度、网格分辨率的收敛行为时，发现结果强烈
依赖 `N_COSTHETA`（从 257 提到 2049 明显改善），指向 cosθ 求积本身有系统
误差而非物理极限的问题；(4) 定位到 denom0 的近极点，实现上面的奇异性减除，
四个极限的相对误差降到 <1%（见 tests/test_p2_vertex.py 与报告的极限核对
表），且不再需要非常大的 `N_COSTHETA`。这是本任务里遇到的唯一实质性数值
bug，未改动任何物理公式或 λ/v_Fq 约定。

数值积分：cosθ 用固定 Gauss-Legendre 求积（`N_COSTHETA` 个节点，默认 257，
奇数保证不含端点 ±1）；ε' 用非均匀复合网格——以 ±FINE_WIDTH_W0·ω₀ 为界，内部
用 `N_EPS_FINE` 个均匀细点覆盖 Eq.(6) 分母在 ε'≈0, ±ω₀ 附近的窄峰，外部用
`N_EPS_COARSE` 个均匀粗点覆盖到能带边界（远离 ω₀ 尺度的区域对积分贡献已被
分母 ~1/ε'² 压低，粗网格足够）；对 ε' 用梯形法（np.trapz，天然支持非均匀
网格）。全部用 numpy 向量化（(n_eps, n_costheta) 矩阵），不逐点调用 scipy.quad，
换取速度（一次 A 网格 ~2700 个 (Q,ω_q,ω₀/E_F) 组合在几十秒内跑完）。
收敛性：test_p2_vertex.py 里对网格分辨率加倍做了交叉检查。

============================================================================
2. H3S DOS 与"等效自由电子 E_F"（SPEC §2 的选择，见函数 load_h3s_dos）
============================================================================

真实 H3S 能带 N(ε) 用于 N̄(ε)（能量-角度积分里的"哪些态参与"），但线性化色散
里的 v_F q=4E_F·Q 这个"Q 的物理刻度"需要一个"等效自由电子 E_F"。本文件用
SPEC §2 给出的更简单选项："取 E_F^free 使 N_free(E_F)=N(E_F)"，展开为：
自由电子气 N_free(E)=(3n)/(2E)（n=电子数/胞），用 H3S.dos 文件里 E_F 处的
积分 DOS 读出 n≈9.00（每个 H3S 胞，与 SPEC 提示的"≈9 电子/胞"一致），
E_F 处的 DOS 列读出 N(E_F)，则
    E_F^free = 1.5 n / N(E_F)
本文件算得 E_F^free ≈ 24.14 eV（band bottom 到 E_F 的真实宽度 ≈25.24 eV，
数量级吻合，是该定义自洽性的一个弱检验）。不需要 scf.in 的晶胞体积。

============================================================================
3. 主要函数
============================================================================

vertex_P(Q, wq, w0, EF, dos, band_bottom=None, band_top=None, T=0.0, wk=None,
         *, eta_k=None, wq_floor=None, exact_dispersion=False,
         n_costheta=N_COSTHETA, n_eps_fine=N_EPS_FINE, n_eps_coarse=N_EPS_COARSE,
         fine_width_w0=FINE_WIDTH_W0, return_imag=False)
    返回 P/λ（float，实部；虚部按理论应 ~0，`return_imag=True` 时额外返回
    (real, imag) 供自检）。band_bottom/band_top 缺省时取半填充对称带
    (-EF, +EF)（对应 SPEC §1 四个解析极限的假设）。

vertex_PV(wq, Qc, w0, EF, dos, *, weight="flat", n_Q=N_Q_AVERAGE, **kwargs)
    P_V(ω_q;Q_c) = Q_c 内的 Q 平均。weight="flat" 是 (1/Qc)∫₀^{Qc}dQ·P(Q,ωq)
    （PGS "耦合在 Q<Qc 内均匀"字面意思）；weight="QdQ" 是
    (2/Qc²)∫₀^{Qc} Q dQ·P(Q,ωq)（三维角度积分自然出现的相空间权重，
    SPEC §3B 允许两种都给）。

命令行入口：`python p2_vertex_numeric.py` 依次算 SPEC §3 的 A-E 并把表/图
写到 theory/derivations/p2_out/。
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Tuple

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass

# --------------------------------------------------------------------------
# 数值default常数（见 docstring §1）
# --------------------------------------------------------------------------
REG_FRAC = 0.03          # eta_k, wq_floor 默认取 REG_FRAC * w0（主 A-E 网格用）
REG_FRAC_LIMIT = 1e-4    # 四个解析极限核对专用的更小正则化尺度（见 docstring §0(c),§1）
N_COSTHETA = 257         # cos(theta) 高斯-勒让德节点数（奇数）
N_EPS_FINE = 2400        # +-FINE_WIDTH_W0*w0 内的细网格点数
N_EPS_COARSE = 800       # 细网格外、到能带边界的粗网格点数（两侧合计）
FINE_WIDTH_W0 = 8.0      # 细网格半宽度，单位 w0
N_Q_AVERAGE = 121        # vertex_PV 里 Q 积分的节点数

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE / "p2_out"

# numpy>=2.0 把 trapz 改名成 trapezoid；本文件统一用这个别名兼容两边。
_trapz = getattr(np, "trapezoid", None) or np.trapz


# ==========================================================================
# 基础统计分布（T=0 解析特例 + 有限 T 的通用形式）
# ==========================================================================

def fermi(eps: np.ndarray, T: float) -> np.ndarray:
    """费米函数。T<=0 时用 T=0 的阶跃 f=theta(-eps)（eps=0 处取 0.5，测度为零不影响积分）。"""
    eps = np.asarray(eps, dtype=float)
    if T is None or T <= 0.0:
        return np.where(eps < 0.0, 1.0, np.where(eps > 0.0, 0.0, 0.5))
    return 1.0 / (np.exp(eps / T) + 1.0)


def bose_pair(w0: float, T: float) -> Tuple[float, float]:
    """返回 (n(w0), n(-w0))。T<=0: n(w0)=0, n(-w0)=-1（=-1-n(w0) 的 T->0 极限）。"""
    if T is None or T <= 0.0:
        return 0.0, -1.0
    n_w0 = 1.0 / (np.exp(w0 / T) - 1.0)
    return n_w0, -1.0 - n_w0


# ==========================================================================
# ε' 积分网格
# ==========================================================================

def _eps_grid(w0: float, band_bottom: float, band_top: float,
              n_fine: int = N_EPS_FINE, n_coarse: int = N_EPS_COARSE,
              fine_width_w0: float = FINE_WIDTH_W0) -> np.ndarray:
    """非均匀复合网格：|eps'|<fine_width_w0*w0 内密，外面稀疏（见模块 docstring §1）。"""
    fine_hi = min(fine_width_w0 * w0, band_top)
    fine_lo = max(-fine_width_w0 * w0, band_bottom)
    if fine_hi <= fine_lo:
        # 能带比 fine window 还窄（理论上不该发生在本任务参数区），退化为单一细网格
        return np.linspace(band_bottom, band_top, n_fine)
    fine = np.linspace(fine_lo, fine_hi, n_fine)
    pieces = [fine]
    if fine_lo > band_bottom:
        pieces.insert(0, np.linspace(band_bottom, fine_lo, max(n_coarse // 2, 2), endpoint=False))
    if fine_hi < band_top:
        pieces.append(np.linspace(fine_hi, band_top, max(n_coarse // 2, 2))[1:])
    eps = np.concatenate(pieces)
    eps = np.unique(eps)
    return eps


_GL_CACHE = {}


def _gauss_legendre(n: int):
    if n not in _GL_CACHE:
        _GL_CACHE[n] = np.polynomial.legendre.leggauss(n)
    return _GL_CACHE[n]


# ==========================================================================
# 核心：一阶顶角函数 P/lambda 的数值积分（GPS99 Eq. 6）
# ==========================================================================

def _vFq(Q: float, EF: float) -> float:
    """v_F q = 4 E_F Q，见模块 docstring §0(a)（更正 SPEC.md 的 2 E_F Q）。"""
    return 4.0 * EF * Q


def vertex_P(
    Q: float,
    wq: float,
    w0: float,
    EF: float,
    dos: Callable[[np.ndarray], np.ndarray],
    band_bottom: Optional[float] = None,
    band_top: Optional[float] = None,
    T: float = 0.0,
    wk: Optional[float] = None,
    *,
    eta_k: Optional[float] = None,
    wq_floor: Optional[float] = None,
    exact_dispersion: bool = False,
    n_costheta: int = N_COSTHETA,
    n_eps_fine: int = N_EPS_FINE,
    n_eps_coarse: int = N_EPS_COARSE,
    fine_width_w0: float = FINE_WIDTH_W0,
    return_imag: bool = False,
):
    """GPS99 Eq.(6) 的数值积分，返回 P/lambda（见模块 docstring §0(b)）。

    Q, wq, w0, EF 单位一致即可（本任务里统一用 eV 或统一用 w0 的倍数都可以，
    函数本身不关心具体单位，只要求 wq, w0, EF 同一单位、Q 无量纲）。

    dos: callable(eps_array) -> N̄(eps)，约定 N̄(0)=1（对应 N(E_F) 归一化，
    §0(b)）。band_bottom/band_top 缺省时取半填充对称带 (-EF, EF)。

    wk: 外线电子频率（若不给，用 eta_k）。T, wk 目前只在 T<=0 且 wk 取 eta_k
    正则化小量的情形下用于本任务（四个解析极限都是 T=0 检验）；T>0/wk 显式
    给出时按通用公式计算，供以后扩展，未做专门收敛性测试。
    """
    if EF <= 0 or w0 <= 0:
        raise ValueError("EF, w0 必须为正")
    eta = eta_k if eta_k is not None else REG_FRAC * w0
    floor = wq_floor if wq_floor is not None else REG_FRAC * w0
    if wk is None:
        wk_eff = eta if T <= 0.0 else max(np.pi * T, 1e-12)
    else:
        wk_eff = wk if wk > 0 else eta
    wq_eff = wq if wq > 0 else floor

    if band_bottom is None:
        band_bottom = -EF
    if band_top is None:
        band_top = EF

    eps = _eps_grid(w0, band_bottom, band_top, n_eps_fine, n_eps_coarse, fine_width_w0)
    x, wts = _gauss_legendre(n_costheta)  # cos(theta) in (-1,1), 权重

    vFq = _vFq(Q, EF)

    eps_col = eps[:, None]              # (n_eps, 1)
    x_row = x[None, :]                  # (1, n_costheta)

    if exact_dispersion:
        if np.any(1.0 + eps / EF < 0):
            raise ValueError("exact_dispersion 要求 eps' > -EF（自由电子模型），"
                              "H3S 真实能带不满足，不要在 H3S DOS 上用这个选项")
        sqrt_term = np.sqrt(1.0 + eps_col / EF)
        eps_q = eps_col + vFq * (Q + x_row * sqrt_term)
        eps_q0 = eps_col + vFq * Q                      # eps_q at costheta=0
    else:
        eps_q = eps_col + vFq * x_row    # (n_eps, n_costheta)
        eps_q0 = eps_col * np.ones_like(eps_col)         # eps_q at costheta=0 (=eps_col)

    f_e = fermi(eps_col, T)              # (n_eps, 1) 广播
    f_eq = fermi(eps_q, T)               # (n_eps, n_costheta)
    f_eq0 = fermi(eps_q0, T)             # (n_eps, 1)，costheta=0 处的值，用于奇异性减除
    n_w0, n_mw0 = bose_pair(w0, T)

    denom0 = -vFq * x_row + 1j * wq_eff  # (1, n_costheta) 广播, ε'-ε'_q + i wq
    term1 = (f_e + n_mw0) / (eps_col + w0 - 1j * wk_eff)
    term2 = (f_e + n_w0) / (eps_col - w0 - 1j * wk_eff)
    term3 = (f_eq + n_mw0) / (eps_q + w0 - 1j * (wk_eff + wq_eff))
    term4 = (f_eq + n_w0) / (eps_q - w0 - 1j * (wk_eff + wq_eff))
    bracket = term1 - term2 - term3 + term4   # (n_eps, n_costheta)

    # --- costheta=0 处的 bracket 值（term1、term2 与 costheta 无关，不必重算）---
    term3_0 = (f_eq0 + n_mw0) / (eps_q0 + w0 - 1j * (wk_eff + wq_eff))
    term4_0 = (f_eq0 + n_w0) / (eps_q0 - w0 - 1j * (wk_eff + wq_eff))
    bracket_0 = term1 - term2 - term3_0 + term4_0   # (n_eps, 1)

    # --- 奇异性减除（见模块 docstring §1 末尾）：denom0 在 costheta=wq_eff/vFq*i 附近
    # 几乎贴着实轴，wq_eff 越小、vFq 越大，这个近极点越窄，固定节点数的 Gauss-Legendre
    # 求积会系统性地积分不准（不是"没收敛"，是"收敛到错误的数"）。把
    # bracket = bracket_0 + (bracket-bracket_0) 拆开：奇异部分 bracket_0/denom0 的
    # costheta 积分有闭式解（对 1/(-vFq*x+i*wq_eff) 积分），正则部分
    # (bracket-bracket_0)/denom0 在 costheta=0 处可去奇点，改用普通 Gauss-Legendre
    # 求积即可，不再受近极点宽度限制。
    if vFq == 0.0:
        I_sing = 2.0 / (1j * wq_eff)
    else:
        a = wq_eff / vFq
        I_sing = (1.0 / vFq) * (np.log(1.0 + 1j * a) - np.log(1j * a - 1.0))
    inner_sing = 0.5 * bracket_0[:, 0] * I_sing                       # (n_eps,)
    regular_integrand = (bracket - bracket_0) / denom0                # (n_eps, n_costheta)
    inner_reg = 0.5 * (regular_integrand * wts[None, :]).sum(axis=1)  # (n_eps,)
    inner = inner_sing + inner_reg

    dos_vals = np.asarray(dos(eps), dtype=float)
    integrand_eps = dos_vals * inner
    outer = _trapz(integrand_eps, eps)

    P_over_lambda = (w0 / 2.0) * outer
    if return_imag:
        return P_over_lambda.real, P_over_lambda.imag
    return float(P_over_lambda.real)


def vertex_PV(
    wq: float,
    Qc: float,
    w0: float,
    EF: float,
    dos: Callable[[np.ndarray], np.ndarray],
    *,
    weight: str = "flat",
    n_Q: int = N_Q_AVERAGE,
    Q_min: float = 0.0,
    **kwargs,
) -> float:
    """Q_c 内的 Q 平均顶角核 P_V(wq;Qc)，SPEC §3B。

    weight="flat": (1/Qc) * integral_0^Qc dQ P(Q,wq)         （PGS 字面意思）
    weight="QdQ" : (2/Qc^2) * integral_0^Qc Q dQ P(Q,wq)      （相空间权重）
    """
    if Qc <= 0:
        raise ValueError("Qc 必须为正")
    Qs = np.linspace(Q_min, Qc, n_Q)
    Ps = np.array([vertex_P(Q, wq, w0, EF, dos, **kwargs) for Q in Qs])
    if weight == "flat":
        avg = _trapz(Ps, Qs) / Qc
    elif weight == "QdQ":
        avg = _trapz(Ps * Qs, Qs) * 2.0 / Qc**2
    else:
        raise ValueError("weight must be 'flat' or 'QdQ'")
    return float(avg)


# ==========================================================================
# 解析极限（用于测试 / 报告核对，SPEC §1）
# ==========================================================================

def analytic_Ps(w0: float, EF: float) -> float:
    """静态极限 P_s/lambda = -w0/(w0+EF)（GPS99 正文，第 3 页）。"""
    return -w0 / (w0 + EF)


def analytic_Pd(w0: float, EF: float) -> float:
    """动态极限 P_d/lambda = +EF/(w0+EF)（GPS99 正文，第 3 页）。"""
    return EF / (w0 + EF)


def analytic_Ppol(w0: float, E: float) -> float:
    """单电子极限 P_pol(q=0,wq=0)/lambda = (E/2)/(w0+E)（GPS99 Eq. 11）。"""
    return (E / 2.0) / (w0 + E)


def vertex_Ppol_direct(w0: float, E: float, n_eps: int = 20000) -> float:
    """直接数值积分 GPS99 Eq.(10) 在 q=0,wq=0,wk=0 的单电子极限，
    与 analytic_Ppol 独立核对 lambda=N(0)g^2 的归一化（模块 docstring §0(b)）。
    能带取 [0,E]（从带底量起，单电子情形没有 E_F 概念，见 docstring）。
    """
    eps = np.linspace(0.0, E, n_eps)
    integrand = 1.0 / (eps + w0) ** 2   # N̄=1（常数 DOS）, g^2=1
    integral = _trapz(integrand, eps)
    return (w0 / 2.0) * integral


# ==========================================================================
# H3S DOS 加载（SPEC §2）
# ==========================================================================

@dataclass
class H3SDos:
    EF_abs: float
    band_bottom: float   # 相对 EF
    band_top: float      # 相对 EF
    N_EF: float           # DOS(EF), states/eV
    n_cell: float          # 积分 DOS 到 EF 的电子数/胞
    EF_free: float         # 等效自由电子 EF（用于 vFq 标度），见 docstring §2
    dos_bar: Callable[[np.ndarray], np.ndarray]   # N(eps)/N(EF)，eps 相对 EF


def load_h3s_dos(path: Path, EF_abs: float = 17.619) -> H3SDos:
    from scipy.interpolate import interp1d

    data = np.loadtxt(path, skiprows=1)
    E_abs, dos_col, idos_col = data[:, 0], data[:, 1], data[:, 2]
    f_dos = interp1d(E_abs, dos_col, kind="linear", bounds_error=False, fill_value=0.0)
    f_idos = interp1d(E_abs, idos_col, kind="linear", bounds_error=False,
                       fill_value=(0.0, idos_col[-1]))
    N_EF = float(f_dos(EF_abs))
    n_cell = float(f_idos(EF_abs))
    EF_free = 1.5 * n_cell / N_EF

    eps_rel = E_abs - EF_abs
    band_bottom = float(eps_rel[0])
    band_top = float(eps_rel[-1])

    def dos_bar(eps: np.ndarray) -> np.ndarray:
        eps = np.asarray(eps, dtype=float)
        return f_dos(eps + EF_abs) / N_EF

    return H3SDos(EF_abs=EF_abs, band_bottom=band_bottom, band_top=band_top,
                  N_EF=N_EF, n_cell=n_cell, EF_free=EF_free, dos_bar=dos_bar)


def const_dos(eps: np.ndarray) -> np.ndarray:
    """常数 DOS，已归一化 N̄=1。"""
    return np.ones_like(np.asarray(eps, dtype=float))


# ==========================================================================
# 声子传播子权重（SPEC §3C 的 P-bar 定义）
# ==========================================================================

def einstein_weight(wq: np.ndarray, w0: float) -> np.ndarray:
    """D(wq) = w0^2/(wq^2+w0^2)，Einstein 声子传播子权重（未归一化）。"""
    return w0 ** 2 / (wq ** 2 + w0 ** 2)


def weighted_Pbar(wq_grid: np.ndarray, PV_grid: np.ndarray, w0: float) -> float:
    """P_V(wq) 在 wq in [0, wq_grid[-1]] 上按 D(wq) 加权平均。"""
    D = einstein_weight(wq_grid, w0)
    num = _trapz(PV_grid * D, wq_grid)
    den = _trapz(D, wq_grid)
    return float(num / den)


# ==========================================================================
# CLI：SPEC §3 A-E
# ==========================================================================

def _task_A(out_dir: Path) -> dict:
    """常数 DOS：P/lambda 在 (Q, wq) 网格上，+ 四个解析极限核对。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    EF = 1.0  # 用 EF 做单位；w0 按 w0/EF 比例给
    ratios = [0.01, 0.1, 0.5]
    Q_grid = np.concatenate([[0.01], np.linspace(0.05, 1.5, 20)])
    wq_grid = np.linspace(0.0, 5.0, 26)   # 单位 w0，稍后乘 w0

    limit_rows = []
    for ratio in ratios:
        w0 = ratio * EF
        wq_abs = wq_grid * w0
        P = np.array([[vertex_P(Q, wq, w0, EF, const_dos) for wq in wq_abs] for Q in Q_grid])

        np.savetxt(out_dir / f"A_grid_w0_over_EF_{ratio:.2f}.csv", P, delimiter=",",
                   header=f"rows=Q{list(Q_grid)}; cols=wq/w0{list(wq_grid)}", comments="# ")

        fig, ax = plt.subplots(figsize=(6, 4.5))
        pc = ax.pcolormesh(wq_grid, Q_grid, P, shading="auto", cmap="RdBu_r",
                            vmin=-np.abs(P).max(), vmax=np.abs(P).max())
        fig.colorbar(pc, ax=ax, label="P/lambda")
        ax.set_xlabel("wq / w0")
        ax.set_ylabel("Q = q/(2kF)")
        ax.set_title(f"P/lambda,  w0/EF={ratio}")
        fig.tight_layout()
        fig.savefig(out_dir / f"A_grid_w0_over_EF_{ratio:.2f}.png", dpi=140)
        plt.close(fig)

        # --- 四个解析极限核对（次序按 GPS99 原文/SPEC.md，见 docstring §0(c)；
        #     REG_FRAC_LIMIT 是比主网格默认 REG_FRAC 更小的正则化尺度，专门
        #     用来干净地分离 Q→0 与 ω_q→0 两个尺度，见 docstring §1） ---
        floor_limit = REG_FRAC_LIMIT * w0

        # P_s = lim_{Q->0}[lim_{wq->0} P(Q,wq)]：内层 wq->0（用 floor_limit
        # 正则化），外层取 Q 落在"远大于 floor、远小于能带边缘效应"的窗口中央。
        Q_static = min(2.0 * ratio, 0.1)
        Ps_num = vertex_P(Q_static, 0.0, w0, EF, const_dos,
                           eta_k=floor_limit, wq_floor=floor_limit)
        Ps_ana = analytic_Ps(w0, EF)

        # P_d = lim_{wq->0}[lim_{Q->0} P(Q,wq)]：内层 Q->0（用极小值代表），
        # 外层取 wq=0.05*w0（远小于 w0、远大于正则化尺度，见 docstring §0(c)）。
        Q_dynamic = 1e-6
        Pd_num = vertex_P(Q_dynamic, 0.05 * w0, w0, EF, const_dos,
                           eta_k=floor_limit, wq_floor=floor_limit)
        Pd_ana = analytic_Pd(w0, EF)

        E_half = 2.0 * EF  # 半填充
        Ppol_num = vertex_Ppol_direct(w0, E_half)
        Ppol_ana = analytic_Ppol(w0, E_half)

        limit_rows.append(dict(
            ratio=ratio, w0=w0,
            Ps_num=Ps_num, Ps_ana=Ps_ana, Ps_relerr=abs(Ps_num - Ps_ana) / abs(Ps_ana),
            Pd_num=Pd_num, Pd_ana=Pd_ana, Pd_relerr=abs(Pd_num - Pd_ana) / abs(Pd_ana),
            Ppol_num=Ppol_num, Ppol_ana=Ppol_ana,
            Ppol_relerr=abs(Ppol_num - Ppol_ana) / abs(Ppol_ana),
        ))


    return {"limits": limit_rows, "Q_grid": Q_grid, "wq_grid": wq_grid, "ratios": ratios}


# Task B/C/D/E 需要对每个 (ratio 或 H3S 参数, Qc, wq) 组合都做一次 Q 平均
# （vertex_PV 内部循环调用 vertex_P 共 N_Q_AVERAGE=121 次），组合数比 Task A
# 的单点网格多一个量级；若照搬 vertex_P 的默认高精度网格（N_EPS_FINE=2400,
# N_EPS_COARSE=800, N_COSTHETA=257，专为四个解析极限在小 Q/小 wq 端的严苛精度
# 设计），CLI 一次跑完 B-E 需要数小时。B/C/D/E 不涉及 Q→0 或 wq→0 的严格极限
# （vertex_PV 的 Q 积分范围是 [0,Qc]，Qc 最小 0.1，且 wq 网格从 0 到 3-5*w0
# 均匀取点，不专门逼近奇异点），不需要四个解析极限那么苛刻的分辨率。用下面
# 这组降低的网格设置（做过收敛性核对，见报告）把单次 vertex_P 调用提速约
# 8-9 倍、vertex_PV 的 Q 积分节点数从 121 降到 61 再提速 2 倍，合计约 17 倍，
# 使 B-E 可以在十几分钟内跑完，而不改动 vertex_P/vertex_PV 本身的默认值
# （Task A 的网格核对、四个解析极限、tests/ 仍然用完整默认精度）。
# 收敛性核对（EF=1,w0=0.1,wq=0.2*w0,Qc=0.5,weight=flat）：
#   完整精度 (N_EPS_FINE=2400,N_EPS_COARSE=800,N_COSTHETA=257,n_Q=121):
#       P_V = -0.024294
#   本设置 (n_eps_fine=600,n_eps_coarse=200,n_costheta=129,n_Q=61):
#       P_V = -0.024104   （相对误差 0.78%，远小于本任务"领头阶量级判断"
#                           所需的精度，不影响 C-E 的定性/半定量结论）
_FAST_GRID_KWARGS = dict(n_eps_fine=600, n_eps_coarse=200, n_costheta=129)
_FAST_N_Q = 61


def _task_B(out_dir: Path) -> dict:
    """Q 平均的顶角核 P_V(wq;Qc)，常数 DOS。（网格设置见 _FAST_GRID_KWARGS 注释）"""
    EF = 1.0
    ratios = [0.01, 0.1, 0.5]
    Qcs = [0.1, 0.3, 0.5, 1.0]
    wq_grid = np.linspace(0.0, 5.0, 21)

    rows = []
    for ratio in ratios:
        w0 = ratio * EF
        wq_abs = wq_grid * w0
        for Qc in Qcs:
            PV_flat = np.array([vertex_PV(wq, Qc, w0, EF, const_dos, weight="flat",
                                           n_Q=_FAST_N_Q, **_FAST_GRID_KWARGS) for wq in wq_abs])
            PV_QdQ = np.array([vertex_PV(wq, Qc, w0, EF, const_dos, weight="QdQ",
                                          n_Q=_FAST_N_Q, **_FAST_GRID_KWARGS) for wq in wq_abs])
            np.savetxt(out_dir / f"B_PV_w0EF{ratio:.2f}_Qc{Qc:.1f}.csv",
                       np.column_stack([wq_grid, PV_flat, PV_QdQ]), delimiter=",",
                       header="wq/w0,PV_flat,PV_QdQ", comments="")
            rows.append(dict(ratio=ratio, w0=w0, Qc=Qc,
                              PV_flat_wq0=PV_flat[0], PV_QdQ_wq0=PV_QdQ[0]))
    return {"rows": rows}


def _h3s_case(out_dir: Path, h3s: H3SDos, w0: float, label: str, Qcs) -> dict:
    """（网格设置见 _FAST_GRID_KWARGS 注释）"""
    wq_grid_over_w0 = np.linspace(0.0, 3.0, 25)
    wq_abs = wq_grid_over_w0 * w0
    rows = []
    for Qc in Qcs:
        PV_flat = np.array([vertex_PV(wq, Qc, w0, h3s.EF_free, h3s.dos_bar, weight="flat",
                                       band_bottom=h3s.band_bottom, band_top=h3s.band_top,
                                       n_Q=_FAST_N_Q, **_FAST_GRID_KWARGS)
                             for wq in wq_abs])
        Pbar_flat = weighted_Pbar(wq_abs, PV_flat, w0)
        PV_QdQ = np.array([vertex_PV(wq, Qc, w0, h3s.EF_free, h3s.dos_bar, weight="QdQ",
                                      band_bottom=h3s.band_bottom, band_top=h3s.band_top,
                                      n_Q=_FAST_N_Q, **_FAST_GRID_KWARGS)
                            for wq in wq_abs])
        Pbar_QdQ = weighted_Pbar(wq_abs, PV_QdQ, w0)
        np.savetxt(out_dir / f"{label}_PV_Qc{Qc:.1f}.csv",
                   np.column_stack([wq_grid_over_w0, PV_flat, PV_QdQ]), delimiter=",",
                   header="wq/w0,PV_flat,PV_QdQ", comments="")
        rows.append(dict(Qc=Qc, Pbar_flat=Pbar_flat, Pbar_QdQ=Pbar_QdQ))
    return {"rows": rows, "w0": w0, "EF_free": h3s.EF_free}


def _task_C_D_E(out_dir: Path, h3s_path: Path) -> dict:
    h3s = load_h3s_dos(h3s_path)
    Qcs = [0.1, 0.3, 0.5, 1.0]
    w0_H = 0.110    # eV, Mishra anharmonic omega_log
    w0_D = 0.110 / np.sqrt(2.0)

    C_const = _h3s_case(out_dir, H3SDos(h3s.EF_abs, -h3s.EF_free, h3s.EF_free, 1.0, 1.0, h3s.EF_free, const_dos),
                         w0_H, "C_const", Qcs)
    C_h3s = _h3s_case(out_dir, h3s, w0_H, "C_h3s", Qcs)
    D_const = _h3s_case(out_dir, H3SDos(h3s.EF_abs, -h3s.EF_free, h3s.EF_free, 1.0, 1.0, h3s.EF_free, const_dos),
                         w0_D, "D_const", Qcs)
    D_h3s = _h3s_case(out_dir, h3s, w0_D, "D_h3s", Qcs)

    summary = []
    for i, Qc in enumerate(Qcs):
        row = dict(
            Qc=Qc,
            Pbar_H_const=C_const["rows"][i]["Pbar_flat"],
            Pbar_H_h3s=C_h3s["rows"][i]["Pbar_flat"],
            Pbar_D_const=D_const["rows"][i]["Pbar_flat"],
            Pbar_D_h3s=D_h3s["rows"][i]["Pbar_flat"],
            Pbar_H_const_QdQ=C_const["rows"][i]["Pbar_QdQ"],
            Pbar_H_h3s_QdQ=C_h3s["rows"][i]["Pbar_QdQ"],
            Pbar_D_const_QdQ=D_const["rows"][i]["Pbar_QdQ"],
            Pbar_D_h3s_QdQ=D_h3s["rows"][i]["Pbar_QdQ"],
        )
        row["ratio_D_over_H_const"] = row["Pbar_D_const"] / row["Pbar_H_const"]
        row["ratio_D_over_H_h3s"] = row["Pbar_D_h3s"] / row["Pbar_H_h3s"]
        row["dos_suppression_H"] = row["Pbar_H_h3s"] / row["Pbar_H_const"]
        row["dos_suppression_D"] = row["Pbar_D_h3s"] / row["Pbar_D_const"]
        row["ratio_D_over_H_const_QdQ"] = row["Pbar_D_const_QdQ"] / row["Pbar_H_const_QdQ"]
        row["ratio_D_over_H_h3s_QdQ"] = row["Pbar_D_h3s_QdQ"] / row["Pbar_H_h3s_QdQ"]
        row["dos_suppression_H_QdQ"] = row["Pbar_H_h3s_QdQ"] / row["Pbar_H_const_QdQ"]
        row["dos_suppression_D_QdQ"] = row["Pbar_D_h3s_QdQ"] / row["Pbar_D_const_QdQ"]
        summary.append(row)

    return {"EF_free": h3s.EF_free, "N_EF": h3s.N_EF, "n_cell": h3s.n_cell,
            "band_bottom": h3s.band_bottom, "band_top": h3s.band_top,
            "w0_H": w0_H, "w0_D": w0_D, "summary": summary}


def _fmt_table(rows, cols, fmt="{:.4g}"):
    lines = [",".join(cols)]
    for r in rows:
        lines.append(",".join(fmt.format(r[c]) if isinstance(r[c], float) else str(r[c]) for c in cols))
    return "\n".join(lines)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    h3s_path = HERE.parent / "data" / "external" / "mishra_h3s_d3s" / \
        "Data-H3S-D3S-PH-EPW-Anharmonic" / "H3S-200GPa" / "band" / "H3S.dos"

    print("=== Task A: 常数 DOS 网格 + 四个解析极限核对 ===")
    A = _task_A(OUT_DIR)
    cols = ["ratio", "w0", "Ps_num", "Ps_ana", "Ps_relerr",
            "Pd_num", "Pd_ana", "Pd_relerr", "Ppol_num", "Ppol_ana", "Ppol_relerr"]
    table = _fmt_table(A["limits"], cols)
    print(table)
    (OUT_DIR / "A_limits_check.csv").write_text(table, encoding="utf-8")

    print("\n=== Task B: Q 平均顶角核 P_V(wq;Qc)，常数 DOS ===")
    B = _task_B(OUT_DIR)
    cols_b = ["ratio", "w0", "Qc", "PV_flat_wq0", "PV_QdQ_wq0"]
    table_b = _fmt_table(B["rows"], cols_b)
    print(table_b)
    (OUT_DIR / "B_summary.csv").write_text(table_b, encoding="utf-8")

    print("\n=== Task C/D/E: H3S 参数、同位素、DOS 效应 ===")
    CDE = _task_C_D_E(OUT_DIR, h3s_path)
    print(f"EF_free = {CDE['EF_free']:.4f} eV, N(EF) = {CDE['N_EF']:.4f} states/eV, "
          f"n_cell = {CDE['n_cell']:.4f}")
    print(f"band (rel EF): [{CDE['band_bottom']:.3f}, {CDE['band_top']:.3f}] eV")
    print(f"w0_H = {CDE['w0_H']} eV, w0_D = {CDE['w0_D']:.5f} eV")
    cols_c = ["Qc", "Pbar_H_const", "Pbar_H_h3s", "Pbar_D_const", "Pbar_D_h3s",
              "ratio_D_over_H_const", "ratio_D_over_H_h3s",
              "dos_suppression_H", "dos_suppression_D"]
    table_c = _fmt_table(CDE["summary"], cols_c)
    print(table_c)
    (OUT_DIR / "CDE_summary.csv").write_text(table_c, encoding="utf-8")

    # SPEC.md §3B 允许两种 Q 权重都给；flat 与 QdQ 权重下 P-bar 的符号/量级
    # 差别很大（见报告），单独存一份 QdQ 版本的表，不覆盖上面 flat 版本。
    cols_c_qdq = ["Qc", "Pbar_H_const_QdQ", "Pbar_H_h3s_QdQ", "Pbar_D_const_QdQ", "Pbar_D_h3s_QdQ",
                  "ratio_D_over_H_const_QdQ", "ratio_D_over_H_h3s_QdQ",
                  "dos_suppression_H_QdQ", "dos_suppression_D_QdQ"]
    table_c_qdq = _fmt_table(CDE["summary"], cols_c_qdq)
    print("\n--- 同一组 C/D/E 量，改用 QdQ 权重（SPEC §3B 的另一种选项） ---")
    print(table_c_qdq)
    (OUT_DIR / "CDE_summary_QdQ.csv").write_text(table_c_qdq, encoding="utf-8")

    print(f"\n输出目录: {OUT_DIR}")


if __name__ == "__main__":
    main()
