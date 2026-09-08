"""
theory/derivations/p2_vertex_dynamic.py

P2 补充任务（team-lead 2026-09-08 指派）：把 P2_vertex_sign_criterion.md
§5 式 (5.1) 的静态极限 DOS 缺损泛函

    P_s/lambda = (w0/2) * integral d(eps) Nbar(eps) / (|eps|+w0)^2  -  1

升级为保留 (Q, wq) 结构、用真实 H3S DOS、按 Eliashberg 核加权的"动态版"，
检验 (5.1) 在多大程度上代表全 (Q, wq) 平面的平均。

============================================================================
0. 依赖与约定
============================================================================

本文件**只 import** theory/derivations/p2_vertex_numeric.py 的公开接口
（vertex_P, vertex_PV, load_h3s_dos, const_dos, einstein_weight,
weighted_Pbar, H3SDos），不修改该文件——另一个 agent 正在用它跑 E_F 扫描
（_task_F）。所有物理公式、v_F q = 4 E_F Q 的约定、lambda = N(0) g^2 的归一化、
静态/动态极限的嵌套次序（Q->0 outer, wq->0 inner 给 P_s；反过来给 P_d，见
p2_vertex_numeric.py 模块 docstring §0(c)/(c')）均直接沿用该模块，不重新推导。

`vertex_P(Q, wq, w0, EF, dos, band_bottom=..., band_top=..., **kwargs)` 的
`dos` 参数本来就接受任意 callable（N(eps)/N(EF) 归一化的表格函数），所以本
文件**不需要**另写一个"vertex_P_dos"——真实 H3S DOS 就是通过
`load_h3s_dos(...).dos_bar` 这个 callable 直接喂给 `vertex_P`/`vertex_PV`
的。常数 DOS 与真实 DOS 用的是同一套积分代码路径，"对拍"退化为同一个函数
在两种 dos 输入下的直接比较（下面 `_sanity_check_dos_callable` 做了这个
形式对照，不是独立实现）。

============================================================================
1. (Q, wq) 网格与两个 DOS 变体
============================================================================

Q ∈ [0.005, 1.0]，对数网格 30 点（team-lead 指定，用于描述性输出/热图）；
wq ∈ [0, 4*w0]，线性网格 25 点（wq=0 这一行由 vertex_P 内部自动 floor 到
REG_FRAC*w0，代表 wq->0 但非零，与主模块 Task A/B 的网格约定一致）。

**数值发现（本文件独立验证，超出原始指派但影响定量结果是否可信）**：
log 网格从 Q=0.005 起步，在 EF~24 eV 时会系统性漏掉 Q<0.005 里的动态区
贡献——crossover 尺度 Q*(wq)=wq/(4EF) 在这个 EF 下常常本身就 <0.005（例如
wq=w0=0.11eV, EF=24.1eV 给 Q*=0.00114）。直接验证：用该 log 网格算
weighted_split 得到的 P_eff 比 vertex_PV+weighted_Pbar 参考值（§4 Task C
本身记录的数字）小 ~4.75 倍。因此本文件另开一个线性网格 Q_GRID_INT=[0,1.0]
121 点（与 vertex_PV 默认 n_Q=121,Q_min=0.0 一致），**只用它做 P_eff/
静态-动态分区这些定量计算**；log 网格 Q_GRID 仅保留用于 CSV/热图这类
展示"P 在多数量级 Q 上如何变化"的描述性输出，不参与任何定量结论。

两种 DOS：
  - "real"：真实 H3S 能带边界 (band_bottom, band_top) 与真实 Nbar(eps) 表
    保持不变，只有传给 vertex_P 的 EF（决定 v_F q = 4*EF*Q 这个运动学标度）
    替换成扫描值——与 p2_vertex_numeric.py `_task_F` 的做法完全一致（该函数
    docstring 明确说"EF^eff 与 DOS 积分范围是两件独立的事"）。
  - "const"：半填充常数 DOS，带宽 ±EF（与 `_task_F` 的 const_scan 用同一
    约定），作为"裸费米能量级"的对照组。

两个 EF：24.1 eV（≈ h3s.EF_free = 24.14 eV，自由电子等效值，与 P2 §4 一致）
与 1.0 eV（对照，§2/§5 反推的 E_F^eff 量级）。

D3S 与 H3S 共用同一份 H3S DOS 表（体系电子结构不随同位素变化，只有 w0 变，
与 §4/§5 的既有约定相同）；w0_H=110 meV，w0_D=77.8 meV（team-lead 指定）。

============================================================================
2. Eliashberg 核加权的有效顶角与静态/动态分区
============================================================================

  P_eff = [ integral dwq D(wq) integral dQ w(Q) P(Q,wq) ]
          / [ integral dwq D(wq) integral dQ w(Q) ]

D(wq) = w0^2/(wq^2+w0^2)（Einstein 传播子权重）；w(Q)：flat -> 1（PGS"Q<Qc
内均匀耦合"字面意思），QdQ -> Q（三维相空间权重）；Qc 取网格上界 1.0（常数
因子在分子分母中约去，不影响 P_eff，只影响归一化，故 flat/QdQ 都不显式写
1/Qc 或 2/Qc^2）。

分区：GPS 静态区对应 v_F q/wq > 1（Q 相对 wq"大"，极限就是 P_s，负号）、动态
区对应 v_F q/wq < 1（Q 相对 wq"小"，极限就是 P_d，正号）——方向已由
p2_vertex_numeric.py 模块 docstring §0(c') 的数值验证钉死，这里直接采用。
边界 Q*(wq) = wq/(4*EF)（由 v_F q = 4 EF Q = wq 解出）在每个 wq 行上用线性
插值精确切分 Q 积分（而不是按最近网格点近似），使
    P_eff = contrib_static + contrib_dynamic   （恒等式，不是近似）
对同一个 (P_eff, contrib_static, contrib_dynamic) 三元组严格成立（同一个
分母、同一套梯形积分）。

============================================================================
3. 静态极限的独立复核（对 (5.1) = -0.103 / -0.085 的验证）
============================================================================

用与 p2_vertex_numeric.py `_task_A` 相同的嵌套极限数值路径（外层 Q->0：取
Q_static = min(2*w0/EF, 0.1)"远大于正则化 floor、远小于带边效应尺度"；内层
wq->0：用 REG_FRAC_LIMIT*w0 的 floor 正则化），但把 dos 换成真实 H3S
Nbar(eps)、EF 换成 h3s.EF_free，直接复核 P2_vertex_sign_criterion.md §5.2
reviewer 交叉核对结论（"收敛到 -0.09~-0.094，与闭式 (5.1) 的 -0.103 同号、
10% 内一致"）。目标值 -0.103 (H3S) / -0.085 (D3S) 直接抄自该文档式 (5.1) 的
QE 全带结果（不 import p2_dos_deficit.py，只引用其已发表的数字作为比较基准）。

============================================================================
4. 输出
============================================================================

CSV/PNG 写到 theory/derivations/p2_out/dyn_*；单测见
theory/derivations/tests/test_p2_vertex_dynamic.py。
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from p2_vertex_numeric import (   # noqa: E402  (只用公开接口，不改该文件)
    vertex_P,
    vertex_PV,
    load_h3s_dos,
    const_dos,
    einstein_weight,
    weighted_Pbar,
    H3SDos,
    REG_FRAC_LIMIT,
)

OUT_DIR = HERE / "p2_out"
H3S_DOS_PATH = (HERE.parent / "data" / "external" / "mishra_h3s_d3s" /
                 "Data-H3S-D3S-PH-EPW-Anharmonic" / "H3S-200GPa" / "band" / "H3S.dos")

_trapz = getattr(np, "trapezoid", None) or np.trapz

# --------------------------------------------------------------------------
# 常量（团队指定的物理参数）
# --------------------------------------------------------------------------
W0 = {"H3S": 0.110, "D3S": 0.0778}          # eV，team-lead 指定（D3S=77.8 meV）
EF_VALUES = {"free": None, "control_1eV": 1.0}   # "free" 在 main() 里填 h3s.EF_free
QC = 1.0
N_Q = 30
N_WQ = 25
# Q ∈ [0.005, 1.0] 对数网格，30 点（team-lead 指定的范围/点数）：只用于"描述性"
# 的 P(Q,wq) 输出表/热图（dyn_Pgrid_*.csv、dyn_*_grid.png），展示 P 在多个数量级
# 的 Q 上如何变化。
Q_GRID = np.logspace(np.log10(0.005), np.log10(1.0), N_Q)
# Q ∈ [0, 1.0] 线性网格，121 点（与 p2_vertex_numeric.py `vertex_PV` 的默认
# n_Q=121、Q_min=0.0 完全一致）：用于 P_eff / 静态-动态分区的**定量**积分。
# 见下面"数值发现"：对数网格从 0.005 起步会系统性漏掉 Q<0.005 里的动态区
# 贡献（在 EF~24 eV 时，动态区宽度 ~ wq/(4EF) 常常本身就 <0.005），必须用
# 从 Q=0 起步的线性网格才能与 vertex_PV 的参考值一致（见 cross_check_vs_vertex_PV
# 与下面 main() 里的核对打印）。
N_Q_INT = 121
Q_GRID_INT = np.linspace(0.0, 1.0, N_Q_INT)

# 降采样网格设置：复刻 p2_vertex_numeric.py `_FAST_GRID_KWARGS` 的数值
# （该模块的注释里给出收敛性核对：相对误差 0.78%，足够本任务"半定量"精度要求）。
# 本文件独立定义同样的数字，不 import 该模块的私有变量（避免耦合到另一个
# agent 可能改动的内部实现细节）。
FAST_KW = dict(n_eps_fine=600, n_eps_coarse=200, n_costheta=129)

# P2_vertex_sign_criterion.md §5.2 式 (5.1)，QE 全带 DOS 的目标值（直接抄自
# 文档，不 import p2_dos_deficit.py）。
STATIC_TARGET = {"H3S": -0.103, "D3S": -0.085}
STATIC_TARGET_REVIEWER_RANGE = (-0.094, -0.09)   # reviewer 用完整数值代码复核给出的区间（H3S）


# ==========================================================================
# DOS 变体构造（沿用 p2_vertex_numeric.py `_task_F` 的约定）
# ==========================================================================

def build_variant(h3s: H3SDos, EF: float, kind: str) -> H3SDos:
    """kind="real"：真实能带边界+真实 Nbar(eps)，只替换 vFq 标度用的 EF。
    kind="const"：半填充常数 DOS，带宽 ±EF。两者都与 _task_F 的
    h3s_scan / const_scan 完全一致的构造方式。
    """
    if kind == "real":
        return H3SDos(h3s.EF_abs, h3s.band_bottom, h3s.band_top,
                       h3s.N_EF, h3s.n_cell, EF, h3s.dos_bar)
    if kind == "const":
        return H3SDos(h3s.EF_abs, -EF, EF, 1.0, 1.0, EF, const_dos)
    raise ValueError(f"unknown kind {kind!r}")


def _sanity_check_dos_callable() -> float:
    """§0：确认 vertex_P 的 `dos` 参数按 callable 使用（不是常数假设写死）。
    对同一个 (Q,wq,w0,EF) 分别传 const_dos 与"手写返回 1.0 的 lambda"，两者
    应逐位相同（这不是物理验证，只是确认 dos 参数确实被当作 callable 调用，
    没有被内部 shortcut 成常数）。
    """
    same_const = lambda eps: np.ones_like(np.asarray(eps, dtype=float))
    a = vertex_P(0.1, 0.05, 0.11, 1.0, const_dos, **FAST_KW)
    b = vertex_P(0.1, 0.05, 0.11, 1.0, same_const, **FAST_KW)
    rel = abs(a - b) / max(abs(a), 1e-12)
    assert rel < 1e-8, f"dos callable 未按预期生效: {a} vs {b}"
    return rel


# ==========================================================================
# 2D 网格：P(Q, wq)
# ==========================================================================

def p_grid(w0: float, variant: H3SDos, Q_grid: np.ndarray, wq_grid: np.ndarray,
            **kwargs) -> np.ndarray:
    """P(Q,wq)/lambda 的 2D 网格，形状 (len(Q_grid), len(wq_grid))。"""
    P = np.empty((len(Q_grid), len(wq_grid)))
    for i, Q in enumerate(Q_grid):
        for j, wq in enumerate(wq_grid):
            P[i, j] = vertex_P(Q, wq, w0, variant.EF_free, variant.dos_bar,
                                band_bottom=variant.band_bottom,
                                band_top=variant.band_top, **kwargs)
    return P


# ==========================================================================
# Eliashberg 核加权平均 + 静态/动态分区（§2 的解析切分）
# ==========================================================================

def weighted_split(w0: float, variant: H3SDos, Q_grid: np.ndarray,
                    wq_grid: np.ndarray, P: np.ndarray, weight: str = "flat"):
    """返回 dict(P_eff, contrib_static, contrib_dynamic, frac_dynamic)。

    P_eff = contrib_static + contrib_dynamic 精确成立（同一分母）。
    frac_dynamic = contrib_dynamic / P_eff（P_eff 过零附近数值不稳定，
    调用方应结合符号判断，见报告）。
    """
    EF = variant.EF_free
    D = einstein_weight(wq_grid, w0)
    if weight == "flat":
        wfun = lambda Q: np.ones_like(Q)
    elif weight == "QdQ":
        wfun = lambda Q: Q
    else:
        raise ValueError("weight must be 'flat' or 'QdQ'")

    Qc = Q_grid[-1]
    Qmin = Q_grid[0]
    dyn_row = np.empty(len(wq_grid))
    stat_row = np.empty(len(wq_grid))

    for j, wq in enumerate(wq_grid):
        Pv = P[:, j]
        Qstar = np.inf if wq <= 0 else wq / (4.0 * EF)
        if Qstar <= Qmin:
            dyn_row[j] = 0.0
            stat_row[j] = _trapz(wfun(Q_grid) * Pv, Q_grid)
        elif Qstar >= Qc:
            dyn_row[j] = _trapz(wfun(Q_grid) * Pv, Q_grid)
            stat_row[j] = 0.0
        else:
            Pstar = np.interp(Qstar, Q_grid, Pv)
            mask_lo = Q_grid < Qstar
            mask_hi = Q_grid > Qstar
            Q_dyn = np.append(Q_grid[mask_lo], Qstar)
            P_dyn = np.append(Pv[mask_lo], Pstar)
            Q_stat = np.insert(Q_grid[mask_hi], 0, Qstar)
            P_stat = np.insert(Pv[mask_hi], 0, Pstar)
            dyn_row[j] = _trapz(wfun(Q_dyn) * P_dyn, Q_dyn)
            stat_row[j] = _trapz(wfun(Q_stat) * P_stat, Q_stat)

    num_dyn = _trapz(dyn_row * D, wq_grid)
    num_stat = _trapz(stat_row * D, wq_grid)
    denom_Q = _trapz(wfun(Q_grid), Q_grid)
    denom = denom_Q * _trapz(D, wq_grid)

    contrib_static = num_stat / denom
    contrib_dynamic = num_dyn / denom
    P_eff = contrib_static + contrib_dynamic
    frac_dynamic = contrib_dynamic / P_eff if abs(P_eff) > 1e-9 else float("nan")
    return dict(P_eff=P_eff, contrib_static=contrib_static,
                contrib_dynamic=contrib_dynamic, frac_dynamic=frac_dynamic)


# ==========================================================================
# 静态极限独立复核（§3）
# ==========================================================================

def static_limit_scan(w0: float, variant: H3SDos,
                       floor_frac: float = REG_FRAC_LIMIT,
                       Q_candidates: np.ndarray = None,
                       hi_precision: bool = False) -> dict:
    """Q->0 外层极限的收敛扫描，返回 dict(Qs, vals)。

    **数值发现（本文件独立做的验证，超出 team-lead 原始指派范围，但是
    §3 静态极限验证能否达到 10% 目标的关键）**：p2_vertex_numeric.py
    `_task_A` 用于常数 DOS 的启发式 Q_static = min(2*w0/EF, 0.1)，在真实
    （能量依赖）DOS、EF~24 eV（H3S 自由电子等效值）的情形下**不适用**——
    对 H3S: Q_static=min(2*0.11/24.14, 0.1)=0.0091，数值给出 P_s=-0.056，
    与式(5.1)目标 -0.103 相差 46%，远超 10% 容差。逐点扫描发现：Q 越小，
    数值反而越接近目标（Q=0.0091→-0.056，Q=0.001→-0.085，Q=2e-4→-0.089~
    -0.093，与 reviewer 的 -0.09~-0.094 吻合），直到 Q 小到 ~1e-4 才开始被
    floor 正则化污染而重新变差。

    物理解释：真实 DOS 有自己的能量尺度 w0（式(5.1)的 1/(|eps|+w0)^2 权重
    在 |eps|~w0 内取样），静态极限数值路径要求 Q 固定非零、看 wq->0 的内层
    极限——这本身没问题，但要让这条数值路径真正代表"P_pol 在 cosθ 积分里只
    感受到 DOS 在 eps'~0 附近的局域结构"，还需要额外满足 v_F q = 4EF·Q
    << w0（否则线性化色散 eps'+v_Fq·cosθ 会把 cosθ 积分展宽到超出 DOS 结构
    尺度，相当于在做一个"模糊过的"静态极限，偏离真正的 eps'-only 依赖）。
    对常数 DOS，没有这个额外尺度，"2*ratio" 启发式只需要避开 floor 和带边
    两个尺度即可（_task_A 就是这么验证的，见下面 const 对照）；对真实 DOS，
    还需要 Q << w0/(4EF)（H3S: w0/(4EF)=0.11/(4*24.14)=0.00114），这正是
    Q_static=0.0091 失败、Q~2e-4 成功的原因。用同一批候选 Q 对常数 DOS 做
    对照扫描，结果确认：常数 DOS 在 Q<~0.001 反而开始变差（无额外尺度约束，
    floor 噪声主导），方向与真实 DOS 相反——见 main() 打印的对照表。

    因此本文件对真实 DOS 用 Q_static ~ 2e-4（而不是模块默认的启发式），并
    用更高精度设置（n_costheta=1025, n_eps_fine=4000, n_eps_coarse=1200，
    只有 2-4 次调用，代价可忽略）抑制在这个更小 Q 下的离散化噪声。
    """
    if Q_candidates is None:
        Q_candidates = np.geomspace(0.02, 5e-5, 12)
    EF = variant.EF_free
    floor = floor_frac * w0
    kw = dict(n_costheta=1025, n_eps_fine=4000, n_eps_coarse=1200) if hi_precision else {}
    vals = np.array([
        vertex_P(Q, 0.0, w0, EF, variant.dos_bar,
                 band_bottom=variant.band_bottom, band_top=variant.band_top,
                 eta_k=floor, wq_floor=floor, **kw)
        for Q in Q_candidates
    ])
    return dict(Qs=Q_candidates, vals=vals)


def static_limit_numeric(w0: float, variant: H3SDos, dos_kind: str = "real",
                          floor_frac: float = REG_FRAC_LIMIT) -> float:
    """外层 Q->0、内层 wq->0（floor 正则化）的静态极限数值。

    dos_kind="real"（真实、能量依赖 DOS）：用 static_limit_scan 发现的
    Q_static~2e-4（见该函数 docstring 的"数值发现"）+ 高精度网格。
    dos_kind="const"（常数 DOS）：沿用 p2_vertex_numeric.py `_task_A` 的
    启发式 Q_static=min(2*w0/EF,0.1)，默认精度（该启发式已被模块自身的
    <1% 误差测试覆盖，不需要重新调）。
    """
    EF = variant.EF_free
    floor = floor_frac * w0
    if dos_kind == "real":
        Q_static = 2e-4
        kw = dict(n_costheta=1025, n_eps_fine=4000, n_eps_coarse=1200)
    else:
        Q_static = min(2.0 * w0 / EF, 0.1)
        kw = {}
    return vertex_P(Q_static, 0.0, w0, EF, variant.dos_bar,
                     band_bottom=variant.band_bottom, band_top=variant.band_top,
                     eta_k=floor, wq_floor=floor, **kw)


# ==========================================================================
# 常数 DOS 交叉核对（与静态极限解析式 -w0/(w0+EF) 对拍 <1%，§0 要求）
# ==========================================================================

def const_dos_cross_check(w0: float, EF: float) -> dict:
    """常数 DOS、EF 任取，静态极限数值 vs 解析式 -w0/(w0+EF)（§0 的对拍要求）。"""
    variant = H3SDos(0.0, -EF, EF, 1.0, 1.0, EF, const_dos)
    num = static_limit_numeric(w0, variant, dos_kind="const")
    ana = -w0 / (w0 + EF)
    return dict(num=num, ana=ana, relerr=abs(num - ana) / abs(ana))


# ==========================================================================
# 与 p2_vertex_numeric.py `vertex_PV`+`weighted_Pbar`（§4 Task C 用的公开
# 接口）的一致性交叉核对：weighted_split 是本文件独立写的 2D 网格+插值切分
# 实现，理应在忽略 (a) Q 网格下限（本文件 0.005 vs vertex_PV 默认 0.0）、
# (b) wq 上限（本文件 4*w0 vs §4 Task C 用的 3*w0）这两处约定差异后，与
# vertex_PV+weighted_Pbar 给出的 P_eff 一致（<5%，两个网格分辨率也不同）。
# ==========================================================================

def cross_check_vs_vertex_PV(w0: float, variant: H3SDos, weight: str = "flat",
                              wq_max_ratio: float = 4.0) -> dict:
    wq_grid = np.linspace(0.0, wq_max_ratio * w0, 25)
    PV = np.array([vertex_PV(wq, QC, w0, variant.EF_free, variant.dos_bar, weight=weight,
                              band_bottom=variant.band_bottom, band_top=variant.band_top,
                              n_Q=121, **FAST_KW) for wq in wq_grid])
    P_eff_ref = weighted_Pbar(wq_grid, PV, w0)
    return dict(P_eff_ref=P_eff_ref)


# ==========================================================================
# main：跑全部组合，写表/图
# ==========================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    print("=== 0. dos callable 形式核对 + 常数 DOS 静态极限对拍 ===")
    rel = _sanity_check_dos_callable()
    print(f"  const_dos vs 手写 callable(=1)：逐位相对误差 {rel:.2e}（应 ~0，确认 dos 按 callable 使用）")
    for EF_test in [1.0, 24.1]:
        cc = const_dos_cross_check(0.110, EF_test)
        print(f"  常数DOS EF={EF_test} eV: 数值 P_s={cc['num']:+.5f}  解析 -w0/(w0+EF)={cc['ana']:+.5f}  相对误差={cc['relerr']:.2%}")

    h3s = load_h3s_dos(H3S_DOS_PATH)
    EF_VALUES["free"] = h3s.EF_free
    print(f"\nH3S DOS: EF_free={h3s.EF_free:.3f} eV, N(EF)={h3s.N_EF:.4f} states/eV, "
          f"n_cell={h3s.n_cell:.3f}, band=[{h3s.band_bottom:.2f},{h3s.band_top:.2f}] eV (rel EF)")

    print("\n=== 0b. weighted_split 与 vertex_PV+weighted_Pbar 的一致性核对 ===")
    print("  （用 Q_GRID_INT：与 vertex_PV 默认 n_Q=121, Q_min=0.0 同一网格，"
          "见模块顶部'数值发现'——log 网格 Q_GRID 会系统性漏掉 Q<0.005 的动态区贡献）")
    v_const = build_variant(h3s, h3s.EF_free, "const")
    wq_probe = np.linspace(0.0, 4.0 * W0["H3S"], N_WQ)
    P_probe_int = p_grid(W0["H3S"], v_const, Q_GRID_INT, wq_probe, **FAST_KW)
    ws = weighted_split(W0["H3S"], v_const, Q_GRID_INT, wq_probe, P_probe_int, weight="flat")
    ref = cross_check_vs_vertex_PV(W0["H3S"], v_const, weight="flat")
    print(f"  H3S 常数DOS EF_free flat Qc=1.0: weighted_split(Q_GRID_INT) P_eff={ws['P_eff']:+.4f}  "
          f"vertex_PV+weighted_Pbar 参考值={ref['P_eff_ref']:+.4f}  "
          f"(§4 Task C 表记录值 +0.0076，Qc=1.0,flat,常数DOS,H)")
    # 同一份数据在 log 网格 Q_GRID 上的对照打印，量化"漏掉多少"（诊断信息，
    # 不用于任何下游定量结果）。
    P_probe_log = p_grid(W0["H3S"], v_const, Q_GRID, wq_probe, **FAST_KW)
    ws_log = weighted_split(W0["H3S"], v_const, Q_GRID, wq_probe, P_probe_log, weight="flat")
    print(f"  （对照：同一计算若误用 Q_GRID(log,[0.005,1]) 会得到 P_eff={ws_log['P_eff']:+.4f}，"
          f"偏差 {abs(ws_log['P_eff']-ref['P_eff_ref'])/abs(ref['P_eff_ref']):.1%} —— "
          f"这就是本文件改用 Q_GRID_INT 做定量计算的原因）")


    print("\n=== 1. 静态极限独立复核（对照 P2_vertex_sign_criterion.md §5.2 式(5.1)） ===")
    print("  （见 static_limit_scan 文档字符串的数值发现：真实DOS+大EF 需要 "
          "Q_static << w0/(4EF)，用模块默认启发式会差 ~46%；下面先打印 H3S 的收敛扫描表，"
          "再用发现的 Q_static=2e-4+高精度网格给出正式复核值）")
    scan = static_limit_scan(W0["H3S"], build_variant(h3s, h3s.EF_free, "real"), hi_precision=True)
    print("  H3S 真实DOS 收敛扫描 (Q_static -> P_s数值，目标 -0.103):")
    for Q_s, v in zip(scan["Qs"], scan["vals"]):
        print(f"    Q={Q_s:.5g}: P_s={v:+.4f}")
    static_rows = []
    for mat in ["H3S", "D3S"]:
        w0 = W0[mat]
        variant = build_variant(h3s, h3s.EF_free, "real")
        num = static_limit_numeric(w0, variant, dos_kind="real")
        target = STATIC_TARGET[mat]
        relerr = abs(num - target) / abs(target)
        static_rows.append(dict(mat=mat, w0=w0, P_s_num=num, P_s_target=target, relerr=relerr))
        print(f"  {mat}: 数值 P_s(真实DOS,EF={h3s.EF_free:.2f}eV,Q_static=2e-4,高精度) = {num:+.4f}   "
              f"式(5.1)目标 = {target:+.4f}   相对误差 = {relerr:.1%}")
    print(f"  （reviewer 用完整代码复核 H3S 得到的区间：{STATIC_TARGET_REVIEWER_RANGE}）")

    print("\n=== 2. (Q,wq) 2D 网格 + Eliashberg 核加权 + 静态/动态分区 ===")
    print("  （P_eff/分区用 Q_GRID_INT 定量计算；dyn_Pgrid_*.csv 用 Q_GRID(log,30点) 做描述性输出）")
    combos = []
    for mat in ["H3S", "D3S"]:
        w0 = W0[mat]
        wq_grid = np.linspace(0.0, 4.0 * w0, N_WQ)
        for ef_label, EF in EF_VALUES.items():
            for dos_kind in ["real", "const"]:
                variant = build_variant(h3s, EF, dos_kind)
                # 描述性输出：log 网格，用于 CSV 展示 P 在多个数量级 Q 上的形状
                P_desc = p_grid(w0, variant, Q_GRID, wq_grid, **FAST_KW)
                np.savetxt(OUT_DIR / f"dyn_Pgrid_{mat}_{ef_label}_{dos_kind}.csv", P_desc,
                           delimiter=",",
                           header=f"rows=Q{list(Q_GRID)}; cols=wq[eV]{list(wq_grid)}",
                           comments="# ")
                # 定量计算：线性网格 Q_GRID_INT（与 vertex_PV 默认网格一致），
                # 避免 log 网格漏掉 Q<0.005 的动态区贡献（见模块顶部"数值发现"）
                P_int = p_grid(w0, variant, Q_GRID_INT, wq_grid, **FAST_KW)
                row = dict(mat=mat, ef_label=ef_label, EF=EF, dos_kind=dos_kind)
                for weight in ["flat", "QdQ"]:
                    res = weighted_split(w0, variant, Q_GRID_INT, wq_grid, P_int, weight=weight)
                    row[f"P_eff_{weight}"] = res["P_eff"]
                    row[f"stat_{weight}"] = res["contrib_static"]
                    row[f"dyn_{weight}"] = res["contrib_dynamic"]
                    row[f"fracdyn_{weight}"] = res["frac_dynamic"]
                combos.append(row)

    cols = ["mat", "ef_label", "EF", "dos_kind",
            "P_eff_flat", "stat_flat", "dyn_flat", "fracdyn_flat",
            "P_eff_QdQ", "stat_QdQ", "dyn_QdQ", "fracdyn_QdQ"]
    lines = [",".join(cols)]
    for r in combos:
        lines.append(",".join(
            f"{r[c]:.5f}" if isinstance(r[c], float) else str(r[c]) for c in cols))
    table = "\n".join(lines)
    print(table)
    (OUT_DIR / "dyn_Pbar_eff_summary.csv").write_text(table, encoding="utf-8")

    # ---- 衍生量：同位素比、DOS 削弱因子 ----
    print("\n=== 3. 衍生量：P_D/P_H、DOS 削弱因子（真实/常数） ===")
    idx = {(r["mat"], r["ef_label"], r["dos_kind"]): r for r in combos}
    deriv_rows = []
    for ef_label in EF_VALUES:
        for dos_kind in ["real", "const"]:
            for weight in ["flat", "QdQ"]:
                H = idx[("H3S", ef_label, dos_kind)][f"P_eff_{weight}"]
                D = idx[("D3S", ef_label, dos_kind)][f"P_eff_{weight}"]
                ratio_DH = D / H if abs(H) > 1e-9 else float("nan")
                deriv_rows.append(dict(ef_label=ef_label, dos_kind=dos_kind, weight=weight,
                                        P_H=H, P_D=D, ratio_D_over_H=ratio_DH))
        for weight in ["flat", "QdQ"]:
            for mat in ["H3S", "D3S"]:
                real = idx[(mat, ef_label, "real")][f"P_eff_{weight}"]
                const = idx[(mat, ef_label, "const")][f"P_eff_{weight}"]
                supp = real / const if abs(const) > 1e-9 else float("nan")
                deriv_rows.append(dict(ef_label=ef_label, dos_kind=f"suppression_{mat}",
                                        weight=weight, P_H=real, P_D=const, ratio_D_over_H=supp))
    cols_d = ["ef_label", "dos_kind", "weight", "P_H", "P_D", "ratio_D_over_H"]
    lines_d = [",".join(cols_d)]
    for r in deriv_rows:
        lines_d.append(",".join(f"{r[c]:.5f}" if isinstance(r[c], float) else str(r[c]) for c in cols_d))
    table_d = "\n".join(lines_d)
    print(table_d)
    (OUT_DIR / "dyn_derived_ratios.csv").write_text(table_d, encoding="utf-8")

    # ---- 图：H3S 真实 DOS、EF=free 的 (Q,wq) 热图 ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        w0 = W0["H3S"]
        wq_grid = np.linspace(0.0, 4.0 * w0, N_WQ)
        variant = build_variant(h3s, h3s.EF_free, "real")
        P = p_grid(w0, variant, Q_GRID, wq_grid, **FAST_KW)
        fig, ax = plt.subplots(figsize=(6, 4.5))
        vmax = np.abs(P).max()
        pc = ax.pcolormesh(wq_grid * 1e3, Q_GRID, P, shading="auto", cmap="RdBu_r",
                            vmin=-vmax, vmax=vmax)
        ax.set_yscale("log")
        fig.colorbar(pc, ax=ax, label="P/lambda")
        ax.set_xlabel("wq (meV)")
        ax.set_ylabel("Q = q/(2kF)  (log)")
        ax.set_title(f"H3S real DOS, EF={h3s.EF_free:.1f} eV, w0={w0*1e3:.0f} meV")
        fig.tight_layout()
        fig.savefig(OUT_DIR / "dyn_H3S_real_EFfree_grid.png", dpi=140)
        plt.close(fig)
    except Exception as exc:  # pragma: no cover
        print(f"  (跳过画图: {exc})")

    print(f"\n输出目录: {OUT_DIR}  (耗时 {time.time()-t0:.1f}s)")


if __name__ == "__main__":
    main()
