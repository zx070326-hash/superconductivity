"""
theory/derivations/p2_vertex_momentum.py

P2 补充任务（team-lead 2026-09-08 指派，v0.5 §5.6 之后的下一步）：
P2_vertex_sign_criterion.md §5.6 用 (Q,wq) 全结构的 Eliashberg 核加权检验推翻了
静态 DOS 缺损泛函的机制解释——各向同性、自由电子色散、常数/真实 H3S DOS 的一阶
顶角理论，在 flat/QdQ 两种"Q<Qc 内均匀耦合"权重下，Eliashberg 加权顶角 P_eff
都是小正数或 ~0，复现不出 EPW 的 -0.095。§5.6 结尾明确指出"负号必须来自各向
同性模型没有的东西：电声耦合 g(q) 的动量结构、真实能带（扩展鞍点/嵌套）"。

本文件做两组独立的探索，回答"缺什么才能给出负号"：

  A. 把 Q 权重从 flat / QdQ 换成物理上更具体的 g²(Q) 模型（小 q 峰、大 q 峰、
     区界主导、以及用 Mishra 数据集 D3S.a2f_proj 的模投影信息构造的
     "更接近 H3S" 的代理权重），看是否存在任何 g²(Q) 使 P_eff < 0。

  B. 把线性化色散 eps_{k'+q} = eps' + v_F q cosθ 换成"部分态 v_F -> 0"的
     鞍点(vHs)色散模型，看 P_eff 随慢电子权重 f_vHs 的变化与符号。

============================================================================
0. 依赖与约定（与 p2_vertex_dynamic.py 完全一致，不重新推导）
============================================================================

本文件**只 import** theory/derivations/p2_vertex_numeric.py 与
theory/derivations/p2_vertex_dynamic.py 的公开接口（无下划线前缀的名字），
不修改这两个文件（另一 agent 正在用 p2_vertex_dynamic.py 的产出被 team-lead
审阅）。所有物理公式、v_F q = 4 E_F Q 的约定、lambda = N(0) g^2 的归一化、
静态/动态极限的嵌套次序、"静态区 <-> v_F q/w_q>1 <-> P<0"、"动态区 <->
v_F q/w_q<1 <-> P>0" 的方向，均直接沿用 p2_vertex_numeric.py 模块 docstring
§0(c)/(c') 与 p2_vertex_dynamic.py 已验证的结论，不重新推导或重新验证
（p2_vertex_dynamic.py 的测试套件已经覆盖）。

Part A 全程复用 p2_vertex_dynamic.build_variant / p_grid（对同一个 (Q,wq)
网格只调用一次 vertex_P 计算 P(Q,wq)，权重是后处理——g²(Q) 权重不改变
P(Q,wq) 本身，只改变对 Q 的加权平均方式，所以可以在同一份 P 网格上套用
任意多个 g²(Q) 而不用重新积分）。

Part B 的"部分态 v_F->0"在 vertex_P 的接口里无法表达（vertex_P 的色散对
所有 k' 态用同一个 v_F q），所以本文件新写 `vertex_P_vfscale`——按
p2_vertex_numeric.py 模块 docstring §1 的积分结构（GPS99 Eq.6，含 §1 末尾
的奇异性减除）独立重新实现，允许 v_F q 被一个 `vf_scale` 因子（标量或
eps 的函数）重新标度。**不从 p2_vertex_numeric.py import 任何下划线前缀
的私有函数**（_eps_grid、_gauss_legendre、fermi、bose_pair 里 fermi/
bose_pair 虽无下划线但本文件仍选择独立重写，避免对被测模块产生隐式耦合，
便于对拍测试的独立性）。`vf_scale ≡ 1.0` 必须与 `vertex_P` 逐点符合
（<1%，见下面 `_sanity_check_vfscale_matches_vertex_P` 与
tests/test_p2_vertex_momentum.py）。

============================================================================
1. Part A：g²(Q) 权重模型
============================================================================

四类权重（不做绝对归一化——P_eff 是权重的比值，常数因子在分子分母中约去）：

  (i)   小 q 峰（前向散射主导）：g²(Q) ∝ exp(-(Q/Q_c)^2)，Q_c ∈ {0.05,0.1,0.2}
  (ii)  大 q 峰（背散射/声学支典型形状）：g²(Q) ∝ Q^2 与 g²(Q) ∝ Q^4 两个变体
  (iii) 区界主导：g²(Q) ∝ 1_{Q∈[0.8,1.0]}（指示函数，其余处为 0）
  (iv)  H3S 数据启发的代理权重：见下面 §1b

以及两个基线（flat、QdQ）作对照，直接与 P2_vertex_sign_criterion.md §5.6 表
的数字延续对比（H3S EF_free real DOS: flat +0.0051, QdQ -0.0002）。

对每个权重，用与 p2_vertex_dynamic.weighted_split 完全相同的解析静态/动态
切分算法（crossover Q*(wq)=wq/(4EF) 用线性插值精确处理，P_eff=
contrib_static+contrib_dynamic 恒等式精确成立），只是把 weight 参数从
字符串（"flat"/"QdQ"）泛化成任意 callable(Q_array)->weight_array
（`weighted_split_general`，本文件新写，用与 p2_vertex_dynamic.weighted_split
在 flat/QdQ 上数值一致来做自洽性检验，见 tests）。

§1b：H3S 数据启发的代理权重（(iv)，reviewer 请注意本节的假设边界）
--------------------------------------------------------------------------
Mishra 数据集里没有真正 q 分辨的电声耦合（`*.a2f_proj` 是**模/支投影**的
alpha^2F(w)——按声子支分解，不是按 q 分解）。本文件用它做一件力所能及的事：
用 D3S-200GPa/EPW_harmonic/D3S.a2f_proj 的 12 个模投影列（D3S 4 原子/胞，
3x4=12 支），按每支的 alpha^2F 谱权重算出"谱权重加权平均频率" mean_freq_i，
把频率最高的 3 支（H3S/D3S 文献里公认的、色散最弱、贯穿全 BZ 的 S-H(D) 键
伸缩光学支——例如 Errea et al. 2016 PRL Fig.1、Flores-Livas et al. Phys.
Rep. 2020 Fig. 12 的声子色散图都显示最上面 2-3 支在 BZ 内几乎平直）标记为
"高频支"，其余 9 支标记为"低频支"，得到 alpha^2F 谱权重意义下的
frac_high/frac_low 分配（这是从数据算出来的数字，不是手调阈值）。

**关键的、无法用这份数据集验证的假设**（必须显式声明，不能含糊过去）：
"高频支色散弱 -> 耦合在整个 BZ（含大 q/区界）近似均匀分布"、"低频支（声学
+ 低光学）色散强 -> 耦合更偏向区界（H3S 的 S 笼弯曲模在 R 点附近有软化，
文献认为这是结构失稳的来源）"——这两条 Q 空间局域化的陈述来自对文献声子
色散图和已发表 Eliashberg 函数分解讨论的**定性复述**，D3S.a2f_proj 本身
不含 q 信息，无法在本文件里独立证实或证伪。因此 (iv) 只是一个"用得上手头
数据、方向上不是瞎猜"的代理，不是从头算的 g(q)，报告里会明确标注这一点。

代理权重构造：g²_proxy(Q) = frac_low * zone_boundary_norm(Q)
                           + frac_high * flat_norm(Q)
两个子形状各自在与 flat/QdQ 相同的 Q 网格上归一化到相同的"总面积"再按
frac_low/frac_high 线性叠加（`make_h3s_proxy_weight`），避免子形状本身的
任意归一化常数污染 frac_low/frac_high 的相对权重。

============================================================================
2. Part B：vHs（鞍点）色散的两个模型
============================================================================

`vertex_P_vfscale(Q, wq, w0, EF, dos, vf_scale, ...)`：把 vertex_P 里的
v_F q = 4 E_F Q 换成 4 E_F Q * vf_scale(eps')（vf_scale 可以是标量，也可以
是 eps' 的函数），其余全部照抄 GPS99 Eq.(6) 的积分结构（含 §1 的奇异性
减除，因为 v_F q 现在逐 eps' 行不同，奇异性减除的解析部分也要逐行算，见
函数实现里的说明）。vf_scale≡1.0 时应与 vertex_P 逐点一致。

模型 1（两态混合，字面对应 team-lead 的表述"一部分态权重 f_vHs 具有
v_F->0，其余态保持自由电子"）：
    P_mix(Q,wq; f) = f * P(Q,wq; vf_scale=0) + (1-f) * P(Q,wq; vf_scale=1)
**这个线性叠加不是近似**：vertex_P 的 eps'-cosθ 积分对被积函数是线性的，
"一部分态"在本文件里精确理解成"总 DOS 权重 N̄(eps') 被分成两个独立通道
（各自复制同一条 N̄(eps')，各自贡献可加），一个通道 v_F=0（真正的鞍点态，
q 无法把它散射到不同能量——粒子-空穴激发能恒为 wq，与 Q 无关）、一个通道
是普通自由电子"——在这个两通道假设下，P_mix 的线性叠加公式是精确恒等式，
不是线性化近似。这是本文件对"部分态"最简单、最可核验的具体化，报告里会
说明这不是唯一可能的具体化（见模型 2）。

关键数值事实：vf_scale=0 时 v_F q ≡ 0，与 Q 完全无关（P(Q,wq;vf_scale=0)
是 Q 的常数函数）——只需在 wq 网格上算一次，不需要在整个 (Q,wq) 网格上算，
本文件按这个事实省了一次 2D 网格计算（用 tests 里的
test_vfscale_zero_is_Q_independent 钉死这个事实，防止未来改动破坏它）。

模型 2（连续 v_F 抑制窗口，附加的稳健性检查，只测 1 组参数）：
    vf_scale(eps') = 1 - f_vHs * exp(-(eps'/w_vHs)^2)
这不是"两态混合"，而是"每个态的 v_F 在费米面附近被连续压低"（更接近
"扩展鞍点在 eps_F 附近展宽"的图像，而不是"两个分立通道"）。w_vHs 取 0.05
eV（Mishra 2608 表 S3 对 H3S 报告 E_F^eff≈0.04 eV 的 vHs 特征宽度量级，见
P2_vertex_sign_criterion.md §5.4），只测 f_vHs=0.5 一组作为对模型 1 结论的
交叉检验，不做完整扫描（时间预算，见下面运行时间说明）。

============================================================================
3. 输出与运行时间
============================================================================

CSV 写到 theory/derivations/p2_out/mom_*.csv；单测见
theory/derivations/tests/test_p2_vertex_momentum.py；报告见
theory/derivations/P2_vertex_momentum_report.md。

网格精度沿用 p2_vertex_dynamic.FAST_KW（n_eps_fine=600, n_eps_coarse=200,
n_costheta=129，该文件已记录收敛性核对：相对误差 0.78%，足够本任务的
"符号 + 半定量量级"精度要求）；Q 网格线性 [0,1] 81 点、wq 网格线性
[0,4*w0] 21 点（比 p2_vertex_dynamic 的 121x25 略降，因为 Part A 的多个
g²(Q) 权重都是对同一份 P(Q,wq) 网格做后处理，网格本身只需算 2-3 次，
不是瓶颈；主要为 Part B 的 vf_scale 网格控制耗时）。全文件（含 CLI 主流程
+ pytest）设计在 15 分钟内跑完，见 main() 末尾打印的耗时。
"""
from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Callable, Optional, Union

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from p2_vertex_numeric import (   # noqa: E402  (只用公开接口，不改该文件)
    vertex_P,
    load_h3s_dos,
    const_dos,
    einstein_weight,
    H3SDos,
    REG_FRAC,
)
from p2_vertex_dynamic import (   # noqa: E402  (只用公开接口，不改该文件)
    build_variant,
    p_grid,
    weighted_split,
    FAST_KW,
    H3S_DOS_PATH,
)

OUT_DIR = HERE / "p2_out"
D3S_A2FPROJ_PATH = (HERE.parent / "data" / "external" / "mishra_h3s_d3s" /
                     "Data-H3S-D3S-PH-EPW-Anharmonic" / "D3S-200GPa" /
                     "EPW_harmonic" / "D3S.a2f_proj")

_trapz = getattr(np, "trapezoid", None) or np.trapz

# --------------------------------------------------------------------------
# 常量
# --------------------------------------------------------------------------
W0 = {"H3S": 0.110, "D3S": 0.0778}
EF_VALUES = {"free": None, "control_1eV": 1.0}   # "free" 在 main() 里填 h3s.EF_free
N_Q = 81
N_WQ = 21
Q_GRID = np.linspace(0.0, 1.0, N_Q)   # 线性、从 0 起步（见 p2_vertex_dynamic 的
                                       # "数值发现"：log 网格从 0.005 起步会系统性
                                       # 漏掉 EF~24eV 时 Q<0.005 的动态区贡献）


# ==========================================================================
# 0. vf_scale 广义顶角积分（Part B 用；Part A 不需要，Part A 全部复用
#    p2_vertex_dynamic.p_grid -> p2_vertex_numeric.vertex_P）
# ==========================================================================

_GL_CACHE: dict = {}


def _gauss_legendre(n: int):
    if n not in _GL_CACHE:
        _GL_CACHE[n] = np.polynomial.legendre.leggauss(n)
    return _GL_CACHE[n]


def _eps_grid(w0: float, band_bottom: float, band_top: float,
              n_fine: int, n_coarse: int, fine_width_w0: float) -> np.ndarray:
    """与 p2_vertex_numeric._eps_grid 相同的非均匀复合网格构造（本文件独立
    重写，不 import 私有函数，见模块 docstring §0）。"""
    fine_hi = min(fine_width_w0 * w0, band_top)
    fine_lo = max(-fine_width_w0 * w0, band_bottom)
    if fine_hi <= fine_lo:
        return np.linspace(band_bottom, band_top, n_fine)
    fine = np.linspace(fine_lo, fine_hi, n_fine)
    pieces = [fine]
    if fine_lo > band_bottom:
        pieces.insert(0, np.linspace(band_bottom, fine_lo, max(n_coarse // 2, 2), endpoint=False))
    if fine_hi < band_top:
        pieces.append(np.linspace(fine_hi, band_top, max(n_coarse // 2, 2))[1:])
    eps = np.concatenate(pieces)
    return np.unique(eps)


def _fermi(eps: np.ndarray, T: float) -> np.ndarray:
    eps = np.asarray(eps, dtype=float)
    if T is None or T <= 0.0:
        return np.where(eps < 0.0, 1.0, np.where(eps > 0.0, 0.0, 0.5))
    return 1.0 / (np.exp(eps / T) + 1.0)


def _bose_pair(w0: float, T: float):
    if T is None or T <= 0.0:
        return 0.0, -1.0
    n_w0 = 1.0 / (np.exp(w0 / T) - 1.0)
    return n_w0, -1.0 - n_w0


def vertex_P_vfscale(
    Q: float,
    wq: float,
    w0: float,
    EF: float,
    dos: Callable[[np.ndarray], np.ndarray],
    vf_scale: Union[float, Callable[[np.ndarray], np.ndarray]] = 1.0,
    band_bottom: Optional[float] = None,
    band_top: Optional[float] = None,
    T: float = 0.0,
    wk: Optional[float] = None,
    *,
    eta_k: Optional[float] = None,
    wq_floor: Optional[float] = None,
    n_costheta: int = 129,
    n_eps_fine: int = 600,
    n_eps_coarse: int = 200,
    fine_width_w0: float = 8.0,
    return_imag: bool = False,
):
    """GPS99 Eq.(6) 的数值积分，v_F q = 4 E_F Q 被逐 eps' 重新标度为
    4 E_F Q * vf_scale(eps')（vf_scale 可以是标量，广播到所有 eps'；也可以
    是 callable(eps_array)->array，代表"不同能量的态有不同的有效 v_F"）。
    vf_scale≡1.0 时应与 p2_vertex_numeric.vertex_P 逐点一致（<1%，见
    tests/test_p2_vertex_momentum.py::test_vfscale_one_matches_vertex_P）。

    实现与 p2_vertex_numeric.vertex_P 的唯一实质区别：原函数假设 vFq 是
    对所有 eps' 相同的标量，奇异性减除的闭式 log 积分只需算一次；这里 vFq
    可以逐 eps' 行不同，所以把 I_sing 的闭式积分向量化到每一行（vFq=0 的
    行单独处理，避免除以零，见下面 zero_mask 分支——vFq=0 时 Eq.(6) 第一个
    分母退化成纯 i*wq_eff，其 costheta 积分是平凡的 2/(i*wq_eff)，见模块
    docstring §2 "关键数值事实"）。
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
    x, wts = _gauss_legendre(n_costheta)

    eps_col = eps[:, None]
    x_row = x[None, :]

    vf_arr = vf_scale(eps) if callable(vf_scale) else np.full_like(eps, float(vf_scale))
    vFq_col = (4.0 * EF * Q) * vf_arr[:, None]      # (n_eps,1)，逐 eps' 行的有效 vFq

    eps_q = eps_col + vFq_col * x_row                # (n_eps, n_costheta)
    eps_q0 = eps_col * np.ones_like(eps_col)          # costheta=0 处 eps_q（与 vFq 无关）

    f_e = _fermi(eps_col, T)
    f_eq = _fermi(eps_q, T)
    f_eq0 = _fermi(eps_q0, T)
    n_w0, n_mw0 = _bose_pair(w0, T)

    denom0 = -vFq_col * x_row + 1j * wq_eff          # (n_eps, n_costheta)
    term1 = (f_e + n_mw0) / (eps_col + w0 - 1j * wk_eff)
    term2 = (f_e + n_w0) / (eps_col - w0 - 1j * wk_eff)
    term3 = (f_eq + n_mw0) / (eps_q + w0 - 1j * (wk_eff + wq_eff))
    term4 = (f_eq + n_w0) / (eps_q - w0 - 1j * (wk_eff + wq_eff))
    bracket = term1 - term2 - term3 + term4

    term3_0 = (f_eq0 + n_mw0) / (eps_q0 + w0 - 1j * (wk_eff + wq_eff))
    term4_0 = (f_eq0 + n_w0) / (eps_q0 - w0 - 1j * (wk_eff + wq_eff))
    bracket_0 = term1 - term2 - term3_0 + term4_0     # (n_eps,1)

    vFq_flat = vFq_col[:, 0]                          # (n_eps,)
    I_sing = np.empty(len(eps), dtype=complex)
    zero_mask = (vFq_flat == 0.0)
    if np.any(zero_mask):
        I_sing[zero_mask] = 2.0 / (1j * wq_eff)
    nz = ~zero_mask
    if np.any(nz):
        a = wq_eff / vFq_flat[nz]
        I_sing[nz] = (1.0 / vFq_flat[nz]) * (np.log(1.0 + 1j * a) - np.log(1j * a - 1.0))

    inner_sing = 0.5 * bracket_0[:, 0] * I_sing
    regular_integrand = np.where(denom0 != 0, (bracket - bracket_0) / denom0, 0.0)
    inner_reg = 0.5 * (regular_integrand * wts[None, :]).sum(axis=1)
    inner = inner_sing + inner_reg

    dos_vals = np.asarray(dos(eps), dtype=float)
    integrand_eps = dos_vals * inner
    outer = _trapz(integrand_eps, eps)

    P_over_lambda = (w0 / 2.0) * outer
    if return_imag:
        return P_over_lambda.real, P_over_lambda.imag
    return float(P_over_lambda.real)


def _sanity_check_vfscale_matches_vertex_P(cases) -> list:
    """vf_scale≡1.0 与 vertex_P 逐点核对（同一批 (Q,wq,w0,EF,dos) 组合，
    同一组降精度网格设置），返回每个 case 的相对误差列表。"""
    rows = []
    for label, Q, wq, w0, EF, dos, bb, bt in cases:
        a = vertex_P(Q, wq, w0, EF, dos, band_bottom=bb, band_top=bt, **FAST_KW)
        b = vertex_P_vfscale(Q, wq, w0, EF, dos, vf_scale=1.0, band_bottom=bb, band_top=bt, **FAST_KW)
        rel = abs(a - b) / max(abs(a), 1e-10)
        rows.append(dict(label=label, vertex_P=a, vertex_P_vfscale=b, relerr=rel))
    return rows


# ==========================================================================
# Part A: g²(Q) 权重
# ==========================================================================

def weight_flat():
    return lambda Q: np.ones_like(np.asarray(Q, dtype=float))


def weight_QdQ():
    return lambda Q: np.asarray(Q, dtype=float)


def weight_small_q(Qc: float):
    """(i) 小 q 峰：g²(Q) ∝ exp(-(Q/Qc)^2)。"""
    return lambda Q: np.exp(-(np.asarray(Q, dtype=float) / Qc) ** 2)


def weight_large_q(power: float):
    """(ii) 大 q 峰：g²(Q) ∝ Q^power（power=2 或 4）。"""
    return lambda Q: np.asarray(Q, dtype=float) ** power


def weight_zone_boundary(lo: float = 0.8, hi: float = 1.0):
    """(iii) 区界主导：g²(Q) ∝ 1_{Q∈[lo,hi]}（指示函数）。"""
    def w(Q):
        Q = np.asarray(Q, dtype=float)
        return np.where((Q >= lo) & (Q <= hi), 1.0, 0.0)
    return w


def analyze_a2f_proj_branch_split(path: Path, n_branch: int = 12, n_high: int = 3) -> dict:
    """读 D3S.a2f_proj 的模投影列，按"谱权重加权平均频率"给 12 支排序，
    取频率最高的 n_high 支为"高频支"（见模块 docstring §1b 的文献依据与
    假设边界），返回 frac_low/frac_high（alpha^2F 谱权重意义下的占比）。
    """
    # 文件末尾有一行 "lambda_int = ...  lambda_sum = ..." 汇总行（列数与数据行
    # 不同，np.loadtxt 会因列数变化报错），跳过表头(1行)与文件最后一行。
    all_lines = Path(path).read_text(encoding="utf-8").splitlines()
    data_lines = [ln for ln in all_lines[1:] if ln.strip() and not ln.strip().startswith("lambda")]
    from io import StringIO
    data = np.loadtxt(StringIO("\n".join(data_lines)))
    w_meV = data[:, 0]
    branch_cols = data[:, 2:2 + n_branch]
    integrals = np.array([_trapz(branch_cols[:, i], w_meV) for i in range(n_branch)])
    mean_freq = np.array([
        _trapz(branch_cols[:, i] * w_meV, w_meV) / integrals[i] if integrals[i] > 1e-10 else np.nan
        for i in range(n_branch)
    ])
    order = np.argsort(mean_freq)
    high_idx = order[-n_high:]
    low_idx = order[:-n_high]
    total = integrals.sum()
    frac_high = float(integrals[high_idx].sum() / total)
    frac_low = 1.0 - frac_high
    return dict(w_meV=w_meV, integrals=integrals, mean_freq=mean_freq,
                high_idx=high_idx, low_idx=low_idx,
                frac_high=frac_high, frac_low=frac_low)


def make_h3s_proxy_weight(Q_grid: np.ndarray, frac_low: float, frac_high: float,
                           zb_lo: float = 0.8, zb_hi: float = 1.0):
    """(iv) H3S 数据启发的代理权重（模块 docstring §1b）：
    g²_proxy(Q) = frac_low * zone_boundary_norm(Q) + frac_high * flat_norm(Q)，
    两个子形状先在 Q_grid 上归一化到相同"总面积"再按 frac_low/frac_high
    线性叠加。归一化常数在构造时（用完整 Q_grid）算好并固定，使返回的
    callable 可以安全地在单点（如 weighted_split_general 处理 crossover 点
    时传入的长度为 1 的数组）上求值而不改变相对权重。
    """
    zb_grid = np.where((Q_grid >= zb_lo) & (Q_grid <= zb_hi), 1.0, 0.0)
    zb_area = max(_trapz(zb_grid, Q_grid), 1e-12)
    flat_area = _trapz(np.ones_like(Q_grid), Q_grid)

    def w(Q):
        Q = np.asarray(Q, dtype=float)
        zb = np.where((Q >= zb_lo) & (Q <= zb_hi), 1.0, 0.0) / zb_area
        flat = np.ones_like(Q) / flat_area
        return frac_low * zb + frac_high * flat
    return w


def weighted_split_general(w0: float, variant: H3SDos, Q_grid: np.ndarray,
                            wq_grid: np.ndarray, P: np.ndarray,
                            weight_func: Callable[[np.ndarray], np.ndarray]) -> dict:
    """p2_vertex_dynamic.weighted_split 的泛化版：把 weight 参数从字符串
    ("flat"/"QdQ") 换成任意 callable(Q_array)->weight_array，其余算法
    （crossover Q*(wq)=wq/(4EF) 线性插值精确切分、P_eff=contrib_static+
    contrib_dynamic 恒等式）逐字照抄 p2_vertex_dynamic.weighted_split（本
    文件不 import 它，因为它的 weight 参数签名是字符串，不能直接复用；
    一致性用 tests 里对 flat/QdQ 两个特例的回归比较来保证，见
    test_weighted_split_general_matches_dynamic_module）。
    """
    EF = variant.EF_free
    D = einstein_weight(wq_grid, w0)
    Qc = Q_grid[-1]
    Qmin = Q_grid[0]
    w_full = weight_func(Q_grid)

    dyn_row = np.empty(len(wq_grid))
    stat_row = np.empty(len(wq_grid))

    for j, wq in enumerate(wq_grid):
        Pv = P[:, j]
        Qstar = np.inf if wq <= 0 else wq / (4.0 * EF)
        if Qstar <= Qmin:
            dyn_row[j] = 0.0
            stat_row[j] = _trapz(w_full * Pv, Q_grid)
        elif Qstar >= Qc:
            dyn_row[j] = _trapz(w_full * Pv, Q_grid)
            stat_row[j] = 0.0
        else:
            Pstar = np.interp(Qstar, Q_grid, Pv)
            wstar = float(weight_func(np.array([Qstar]))[0])
            mask_lo = Q_grid < Qstar
            mask_hi = Q_grid > Qstar
            Q_dyn = np.append(Q_grid[mask_lo], Qstar)
            P_dyn = np.append(Pv[mask_lo], Pstar)
            w_dyn = np.append(w_full[mask_lo], wstar)
            Q_stat = np.insert(Q_grid[mask_hi], 0, Qstar)
            P_stat = np.insert(Pv[mask_hi], 0, Pstar)
            w_stat = np.insert(w_full[mask_hi], 0, wstar)
            dyn_row[j] = _trapz(w_dyn * P_dyn, Q_dyn)
            stat_row[j] = _trapz(w_stat * P_stat, Q_stat)

    num_dyn = _trapz(dyn_row * D, wq_grid)
    num_stat = _trapz(stat_row * D, wq_grid)
    denom_Q = _trapz(w_full, Q_grid)
    denom = denom_Q * _trapz(D, wq_grid)

    contrib_static = num_stat / denom
    contrib_dynamic = num_dyn / denom
    P_eff = contrib_static + contrib_dynamic
    frac_dynamic = contrib_dynamic / P_eff if abs(P_eff) > 1e-9 else float("nan")
    return dict(P_eff=P_eff, contrib_static=contrib_static,
                contrib_dynamic=contrib_dynamic, frac_dynamic=frac_dynamic)


def build_qweight_dict(a2f_split: dict) -> dict:
    """Part A 用到的全部权重（含两个基线 flat/QdQ），返回
    {name: callable(Q_array)->weight_array}。"""
    weights = {
        "flat": weight_flat(),
        "QdQ": weight_QdQ(),
        "small_q_Qc0.05": weight_small_q(0.05),
        "small_q_Qc0.10": weight_small_q(0.10),
        "small_q_Qc0.20": weight_small_q(0.20),
        "large_q_Q2": weight_large_q(2.0),
        "large_q_Q4": weight_large_q(4.0),
        "zone_boundary_0.8_1.0": weight_zone_boundary(0.8, 1.0),
        "h3s_proxy": make_h3s_proxy_weight(Q_GRID, a2f_split["frac_low"], a2f_split["frac_high"]),
    }
    return weights


# ==========================================================================
# Part B: vHs 色散
# ==========================================================================

def vHs_population_mix(Q_grid: np.ndarray, wq_grid: np.ndarray, w0: float,
                        variant: H3SDos, f_vHs: float,
                        P_normal: np.ndarray, P_vHs_wq: np.ndarray) -> np.ndarray:
    """模型 1（两态混合，模块 docstring §2）：
        P_mix(Q,wq) = f_vHs * P_vHs(wq) + (1-f_vHs) * P_normal(Q,wq)
    P_vHs_wq 是 vf_scale=0 时的 P(wq)（与 Q 无关，形状 (n_wq,)），
    P_normal 是标准 vertex_P 的 (Q,wq) 网格。"""
    return f_vHs * P_vHs_wq[None, :] + (1.0 - f_vHs) * P_normal


def vHs_localized_window_grid(Q_grid: np.ndarray, wq_grid: np.ndarray, w0: float,
                               variant: H3SDos, f_vHs: float, w_vHs: float,
                               **kwargs) -> np.ndarray:
    """模型 2（连续 v_F 抑制窗口，模块 docstring §2）：
    vf_scale(eps') = 1 - f_vHs * exp(-(eps'/w_vHs)^2)，对整张 (Q,wq) 网格
    重新算 vertex_P_vfscale（Q 依赖没有像模型 1 那样退化，必须算全网格）。
    """
    def vf_scale(eps):
        return 1.0 - f_vHs * np.exp(-(np.asarray(eps, dtype=float) / w_vHs) ** 2)

    P = np.empty((len(Q_grid), len(wq_grid)))
    for i, Q in enumerate(Q_grid):
        for j, wq in enumerate(wq_grid):
            P[i, j] = vertex_P_vfscale(Q, wq, w0, variant.EF_free, variant.dos_bar,
                                        vf_scale=vf_scale,
                                        band_bottom=variant.band_bottom,
                                        band_top=variant.band_top, **kwargs)
    return P


# ==========================================================================
# main
# ==========================================================================

def _write_csv(path: Path, rows: list, cols: list):
    lines = [",".join(cols)]
    for r in rows:
        lines.append(",".join(f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c]) for c in cols))
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    h3s = load_h3s_dos(H3S_DOS_PATH)
    EF_VALUES["free"] = h3s.EF_free
    w0 = W0["H3S"]
    wq_grid = np.linspace(0.0, 4.0 * w0, N_WQ)
    print(f"H3S DOS: EF_free={h3s.EF_free:.3f} eV, band=[{h3s.band_bottom:.2f},{h3s.band_top:.2f}] eV (rel EF)")

    # ---- 0. vf_scale=1 与 vertex_P 的对拍 ----
    print("\n=== 0. vertex_P_vfscale(vf_scale=1.0) 与 vertex_P 对拍 ===")
    variant_free = build_variant(h3s, h3s.EF_free, "real")
    variant_1eV = build_variant(h3s, 1.0, "real")
    variant_const_free = build_variant(h3s, h3s.EF_free, "const")
    cases = [
        ("H3S real EF_free Q=0.1 wq=0.3w0", 0.1, 0.3 * w0, w0, h3s.EF_free, variant_free.dos_bar,
         variant_free.band_bottom, variant_free.band_top),
        ("H3S real EF_free Q=0.02 wq=1e-3", 0.02, 1e-3, w0, h3s.EF_free, variant_free.dos_bar,
         variant_free.band_bottom, variant_free.band_top),
        ("H3S real EF=1eV Q=0.3 wq=2w0", 0.3, 2.0 * w0, w0, 1.0, variant_1eV.dos_bar,
         variant_1eV.band_bottom, variant_1eV.band_top),
        ("const EF_free Q=0.5 wq=w0", 0.5, w0, w0, h3s.EF_free, variant_const_free.dos_bar,
         variant_const_free.band_bottom, variant_const_free.band_top),
    ]
    check_rows = _sanity_check_vfscale_matches_vertex_P(cases)
    for r in check_rows:
        print(f"  {r['label']}: vertex_P={r['vertex_P']:+.6f}  vertex_P_vfscale(1.0)={r['vertex_P_vfscale']:+.6f}  relerr={r['relerr']:.2e}")
        assert r["relerr"] < 0.01, f"vf_scale=1 与 vertex_P 不一致: {r}"

    # ---- 1. Part A: 计算一次 P(Q,wq) 网格（H3S real DOS, 两个 EF），套用多种 g^2(Q) ----
    print("\n=== 1. Part A: g^2(Q) 权重扫描（同一份 P(Q,wq) 网格，后处理不同权重） ===")
    a2f_split = analyze_a2f_proj_branch_split(D3S_A2FPROJ_PATH)
    print(f"  D3S.a2f_proj 12 支按谱权重平均频率排序，高频 3 支 frac_high={a2f_split['frac_high']:.3f}, "
          f"低频 9 支 frac_low={a2f_split['frac_low']:.3f}")
    print(f"  mean_freq (meV，按 index 排序): {np.round(a2f_split['mean_freq'], 1)}")
    print(f"  高频支 index(0-based): {sorted(a2f_split['high_idx'].tolist())}")

    weights = build_qweight_dict(a2f_split)

    qgrid_rows = []
    P_by_ef = {}
    for ef_label, EF in EF_VALUES.items():
        variant = build_variant(h3s, EF, "real")
        print(f"\n  -- EF={ef_label} ({EF:.4f} eV) --")
        P = p_grid(w0, variant, Q_GRID, wq_grid, **FAST_KW)
        P_by_ef[ef_label] = (variant, P)
        for wname, wfun in weights.items():
            res = weighted_split_general(w0, variant, Q_GRID, wq_grid, P, wfun)
            sign = "neg" if res["P_eff"] < 0 else ("zero" if abs(res["P_eff"]) < 1e-6 else "pos")
            qgrid_rows.append(dict(weight=wname, ef_label=ef_label, EF=EF,
                                    P_eff=res["P_eff"], contrib_static=res["contrib_static"],
                                    contrib_dynamic=res["contrib_dynamic"], sign=sign))
            print(f"    {wname:24s}: P_eff={res['P_eff']:+.5f}  stat={res['contrib_static']:+.5f}  "
                  f"dyn={res['contrib_dynamic']:+.5f}  sign={sign}")

    cols_a = ["weight", "ef_label", "EF", "P_eff", "contrib_static", "contrib_dynamic", "sign"]
    _write_csv(OUT_DIR / "mom_qweights_summary.csv", qgrid_rows, cols_a)

    a2f_rows = [dict(branch=i, integral=float(a2f_split["integrals"][i]),
                      mean_freq_meV=float(a2f_split["mean_freq"][i]),
                      classification=("high" if i in a2f_split["high_idx"] else "low"))
                for i in range(len(a2f_split["integrals"]))]
    _write_csv(OUT_DIR / "mom_a2f_proj_branch_split.csv", a2f_rows,
               ["branch", "integral", "mean_freq_meV", "classification"])

    any_negative = any(r["P_eff"] < 0 for r in qgrid_rows)
    print(f"\n  Part A 结论：存在使 P_eff<0 的 g^2(Q) 权重？ {any_negative}")

    # ---- 2. Part B: vHs 色散 ----
    print("\n=== 2. Part B: vHs (v_F->0) 色散扫描 ===")
    # vf_scale=0 只是 wq 的函数（与 Q 无关，见模块 docstring §2 的"关键数值事实"）
    P_vHs_wq = {}
    for ef_label in EF_VALUES:
        variant, _ = P_by_ef[ef_label]
        vals = np.array([
            vertex_P_vfscale(0.37, wq, w0, variant.EF_free, variant.dos_bar, vf_scale=0.0,
                              band_bottom=variant.band_bottom, band_top=variant.band_top, **FAST_KW)
            for wq in wq_grid
        ])
        # 用另一个 Q 值复核 Q 无关性（钉死"关键数值事实"，tests 里也有正式版本）
        vals_check = np.array([
            vertex_P_vfscale(0.83, wq, w0, variant.EF_free, variant.dos_bar, vf_scale=0.0,
                              band_bottom=variant.band_bottom, band_top=variant.band_top, **FAST_KW)
            for wq in wq_grid
        ])
        maxdiff = float(np.max(np.abs(vals - vals_check)))
        print(f"  vf_scale=0 的 Q 无关性核对 ({ef_label}): max|P(Q=0.37)-P(Q=0.83)| = {maxdiff:.2e}")
        assert maxdiff < 1e-6, "vf_scale=0 时 P 应与 Q 无关"
        P_vHs_wq[ef_label] = vals

    vhs_rows = []
    for ef_label, EF in EF_VALUES.items():
        variant, P_normal = P_by_ef[ef_label]
        for f_vHs in [0.0, 0.3, 0.5]:
            P_mix = vHs_population_mix(Q_GRID, wq_grid, w0, variant, f_vHs,
                                        P_normal, P_vHs_wq[ef_label])
            for wname in ["flat", "QdQ"]:
                res = weighted_split_general(w0, variant, Q_GRID, wq_grid, P_mix, weights[wname])
                sign = "neg" if res["P_eff"] < 0 else ("zero" if abs(res["P_eff"]) < 1e-6 else "pos")
                vhs_rows.append(dict(model="pop_mix", ef_label=ef_label, EF=EF, f_vHs=f_vHs,
                                      weight=wname, P_eff=res["P_eff"],
                                      contrib_static=res["contrib_static"],
                                      contrib_dynamic=res["contrib_dynamic"], sign=sign))
                print(f"  [pop_mix] EF={ef_label:12s} f_vHs={f_vHs:.1f} weight={wname:4s}: "
                      f"P_eff={res['P_eff']:+.5f}  sign={sign}")

    # 模型 2：连续抑制窗口，只测 1 组（EF=free, f_vHs=0.5, w_vHs=0.05eV），时间预算
    print("\n  -- 模型 2（连续 v_F 抑制窗口）：只测 EF=free, f_vHs=0.5, w_vHs=0.05eV --")
    variant, _ = P_by_ef["free"]
    P_localized = vHs_localized_window_grid(Q_GRID, wq_grid, w0, variant,
                                             f_vHs=0.5, w_vHs=0.05, **FAST_KW)
    for wname in ["flat", "QdQ"]:
        res = weighted_split_general(w0, variant, Q_GRID, wq_grid, P_localized, weights[wname])
        sign = "neg" if res["P_eff"] < 0 else ("zero" if abs(res["P_eff"]) < 1e-6 else "pos")
        vhs_rows.append(dict(model="localized_window", ef_label="free", EF=h3s.EF_free, f_vHs=0.5,
                              weight=wname, P_eff=res["P_eff"], contrib_static=res["contrib_static"],
                              contrib_dynamic=res["contrib_dynamic"], sign=sign))
        print(f"  [localized_window] EF=free f_vHs=0.5 w_vHs=0.05eV weight={wname:4s}: "
              f"P_eff={res['P_eff']:+.5f}  sign={sign}")

    cols_b = ["model", "ef_label", "EF", "f_vHs", "weight", "P_eff",
              "contrib_static", "contrib_dynamic", "sign"]
    _write_csv(OUT_DIR / "mom_vHs_summary.csv", vhs_rows, cols_b)

    any_negative_b = any(r["P_eff"] < 0 for r in vhs_rows if r["f_vHs"] > 0)
    print(f"\n  Part B 结论：f_vHs>0 时存在使 P_eff<0 的情形？ {any_negative_b}")

    print(f"\n输出目录: {OUT_DIR}  (耗时 {time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
