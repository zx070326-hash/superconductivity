"""P4 §5.1-5.2：数值检验 van Hove 窗口平均增强因子 E(delta)，并标定常数 c。

解析预言（`theory/derivations/P4_vhs_dtc_dp_sign.md` §2）：

  N(xi)/N_b = 1 + a*Gam^2/((xi+dlt)^2+Gam^2)                        (1)
  Ntilde(W) = (W/pi) int dxi (N/N_b)/(xi^2+W^2)                      (2)
            = 1 + a*Gam*(Gam+W)/(dlt^2+(Gam+W)^2)                    (3)
  E(dlt)    = 1 + a*Gam*W_eff/(dlt^2+W_eff^2),  W_eff = Gam + c*w_log (4)
  lam_eff(dlt) = lam_b * E(dlt),   sign(dE/ddlt) = -sign(dlt)        (6)

本脚本用 `theory/solvers/eliashberg_iso.py` 的能量依赖 DOS 路径（Mitrovic-Carbotte /
Pickett 三方程 Z, chi, phi）扫 dlt，反解 E_num(dlt)，拟合 W_eff，标定 c。

------------------------------------------------------------------------------
关键实现决策（读结果前必须知道）
------------------------------------------------------------------------------
1) **N_b 归一化 vs N(0) 归一化：lam 输入必须乘 n0。**
   求解器的 DOS 接口用的是 Nbar(eps) = N(eps)/N(0)（费米面归一化，`DOS` 类的
   核心不变量），而 P4 (1) 式用的是 N/N_b（背景归一化）。两者相差常数
       n0 = N(0)/N_b = 1 + a*Gam^2/(dlt^2+Gam^2)
   电声耦合谱 alpha^2F ~ N(0)|g|^2，所以"同一材料、峰相对 E_F 移动、矩阵元与
   声子不变"这个物理设定对应
       lam_in = lam_b * n0,      Nbar = (N/N_b)/n0
   于是求解器核里 lam_in * Ntilde_solver = lam_b*n0 * Ntilde/n0 = lam_b*Ntilde，
   正是 P4 (4) 式。**如果忘掉 n0 因子（固定 lam_in = lam_b 只换 Nbar），
   N(0) 归一化会把"峰落在 E_F 抬高 DOS"这个效应本身除掉**——这正是
   `theory/solvers/README.md`「已知局限」第 1 条记录的"delta=0 反而压低 Tc"，
   那是归一化约定的产物，不是物理。

2) **窗口因子同时乘 lam 和 mu*，所以 (4) 式"lam_eff = lam_b*E"只在 mu*=0 时严格。**
   求解器（和 P4 §2）的 phi 方程核是 [lam(n-m) - mu*] * Ntilde(|wt_m|)，窗口平均
   是对**整个括号**作用的。两种 mu* 约定：
     mu_mode="scaled": mu*_in = mu*_b * n0（裸 mu ~ N(0)，与 lam 同步缩放）
        => 核 = (lam_b*shape - mu*_b) * Ntilde，逐字复现 P4 §2 的方程。
     mu_mode="fixed" : mu*_in = mu*_b（P4 v0 的"mu* 折进背景"）
   相应地有两种反解：
     E_kernel: 常数 DOS 参考解里把 (lam, mu*) 一起乘 s，找 s 使 Tc 相符
               —— 与 mu_mode="scaled" 严格匹配，是 (3)/(4) 式的干净检验。
     E_lambda: 常数 DOS 参考解里只缩放 lam、mu* 固定，找 lam_eff/lam_b
               —— 就是字面意义的"有效耦合"，但把 mu* 的窗口增强也吸进去了。
   mu*=0 时两者恒等。主结论用 E_kernel + mu_mode="scaled"，并报告 mu*=0 的对照。

3) **参考态用同一条数值路径的平坦 DOS（a=0），不用折叠常数 DOS 路径。**
   能量依赖 DOS 路径有 e_max = 20*w_cut 的有限积分范围截断（~0.6% 的 Tc 系统偏差，
   见 solvers/README.md「积分范围」）。E_num 是比值，用 a=0 的 `PeakDOS` 作参考
   可让这个偏差在分子分母间抵消（实测：flat table Tc=156.00 K vs 折叠路径
   155.10 K，w_cut=5*w_log、lam=1.5、mu*=0.13）。

4) 求积精度：tan 代换 + Gauss-Legendre 对本文件的有理型 Nbar 是谱精度的，
   `dos_n_eps=400` 时 G_j 与高精度 quad 参考的相对误差在 |dlt| <= 500 meV
   （本脚本生产用的 dlt 范围）内是 ~1e-14。**但它随 dlt 变差**（远处窄峰在
   theta 变量里被压成极窄的一段）：|dlt|=1000 meV 时 1.3e-7，2000 meV 时
   2.4e-5，4000 meV 时 2.2e-4。要扫到 |dlt| > 1 eV 必须同步加大 dos_n_eps
   （4000 meV 处 n_eps=6400 恢复到 5e-15）。

------------------------------------------------------------------------------
结论摘要（2026-09-08 全量运行，见 p4_out/*.csv 与 §5 的完整表）
------------------------------------------------------------------------------
* **成立**：dE/ddlt 的符号 (6) 式、峰位、偶性。E(dlt) 对对称峰形**严格**是偶
  函数（实测 max|E(+dlt)-E(-dlt)| = 1.4e-14，是 xi->-xi/dlt->-dlt/chi->-chi
  精确对称性的推论），极大精确在 dlt=0；28 个有限差分区间无一违反
  sign(dE/ddlt) = -sign(dlt)，且 E 沿 |dlt| 严格单调。
* **成立且是精确的**：(2) 式的"窗口平均"图像。求解器的 G_m 与闭式
  Ntilde(a_m; dlt-chi_m)/(n0 a_m) 在 2e-5 上一致（a_m = |w_m|Z_m）。
  **但窗口中心被 chi_m 挪走**：忽略 chi 用 dlt 代替 dlt-chi_m 会差 3%-6%
  （dlt=300 meV 处 max|chi_m| = 62 meV，是 dlt 的 21%）。P4 §2 没有 chi。
* **不成立**：(4) 式"E(dlt) 是单个洛伦兹、宽度 W_eff = Gam + c*w_log"。
  E(dlt) 是**一族宽度 a_m = |w_m|Z_m 各不相同的洛伦兹的加权和**，尾部被最宽的
  窗主导，比单洛伦兹重。量化（基线，dlt in [-500,500]）：
    - 振幅锁定 a*Gam 的 1-参数拟合 R2 = 0.65（远低于 0.99）；
    - 由峰值定的宽度 W_amp = 204 meV vs 由曲率定的 W_curv = 413 meV
      （比值 2.02）vs 半高宽 W_half = 319 meV（1.56）——三者本应相等；
    - 核心/尾部分别拟合给 W_core/W_tail = 0.73；
    - 放开振幅的 3-参数洛伦兹 R2 = 0.973，但要 1.52 倍的振幅。
  这个失效在 23 组参数上是系统性的（W_curv/W_amp 中位 1.86，范围 1.42-2.15）。
* **c 是 O(1) 但不是常数**：单洛伦兹拟合给 c = 0.95 +- 0.47（23 组），主要
  依赖 lam_b（lam_b=1.0/1.5/2.2/3.0 -> c = 0.48/0.92/1.61/2.48）——因为窗宽
  a_m ~ pi*T_c*Z 本身随 lam 增大。零拟合的 c_mean = <a_m>_w/w_log 同样随
  lam_b 从 0.79 涨到 4.10。c 对 w_cut（+-60%）与整体能标缩放不敏感
  （w_log 55/110/220 meV 给出逐位相同的 c，标度不变性的严格检验）。
* **mu* 约定影响很大**：mu*_in 固定（不随 n0 缩放）时拟合出 c = -0.05、
  R2 = 0.50——因为窗口因子同时乘 mu*，固定 mu* 的反解把这部分也吸了进去。
* **P4 (5) 式的前因子偏大 ~2 倍**：字面的 Lam*lam_b（Lam = dlnTc/dlam_eff，
  固定 mu*）= 0.52，而正确的 dlnTc/dE（E 同时乘 lam 与 mu*）= 0.24。

运行：python theory/derivations/p4_vhs_enhancement.py [--quick]
输出：theory/derivations/p4_out/{E_of_delta,c_calibration,
      sensitivity_delta_a_gamma,p4_diagnostics}.csv
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy import interpolate, optimize

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parents[0] / "solvers"))

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from eliashberg_iso import KB_MEV_PER_K, Solver, Spectrum, mu_star_at_cutoff  # noqa: E402

OUT_DIR = _HERE / "p4_out"

# ---------------------------------------------------------------------------
# 基线参数（H3S 量级，meV；见 P4 §0 与 p2_dos_deficit.py 的 w_log=110 meV）
# ---------------------------------------------------------------------------
W_LOG_0 = 110.0     # Einstein 模频率 = w_log = w_2（Einstein 谱下 w_2/w_log=1，f2=1）
LAM_B_0 = 1.5
MU_B_0 = 0.13
A_0 = 2.0
GAM_0 = 100.0
W_CUT_RATIO_0 = 5.0  # w_cut = 5*w_log，沿用 solvers/tests 与 examples/dos_scan.py 的约定
N_EPS = 400


# ===========================================================================
# 峰形与模型 DOS
# ===========================================================================

def shape_lorentz(u: np.ndarray, gamma: float, **kw) -> np.ndarray:
    """洛伦兹峰形 s(u) = Gam^2/(u^2+Gam^2)，s(0)=1，HWHM=Gam。P4 (1) 式。"""
    return gamma ** 2 / (u ** 2 + gamma ** 2)


def shape_cusp(u: np.ndarray, gamma: float, **kw) -> np.ndarray:
    """指数尖点峰形 s(u) = exp(-|u|/Gam)：对称但峰顶有 cusp、尾巴比洛伦兹轻。

    用于 P4 §7「峰形敏感性」：真实 3D 鞍点是单边 sqrt 型、2D 是对数型，都在峰顶
    非光滑。这里取一个可积、对称、有 cusp 的替身——它的 dlt-依赖**不是**洛伦兹
    （指数与洛伦兹窗的卷积），所以能分离"窗口平均图像对不对"与"E 在 dlt 上是不是
    洛伦兹"这两个不同的命题。
    """
    return np.exp(-np.abs(u) / gamma)


def shape_asym(u: np.ndarray, gamma: float, asym_ratio: float = 4.0, **kw) -> np.ndarray:
    """左右宽度不等的洛伦兹：u<0（峰的低能侧）宽度 Gam，u>0 宽度 Gam/asym_ratio。

    s(0)=1 连续，极大严格在 u=0。用于检验 P4 §7 的预言"单边/不对称峰会让 E(dlt)
    的极大偏离 dlt=0"。
    """
    g = np.where(u < 0.0, gamma, gamma / asym_ratio)
    return g ** 2 / (u ** 2 + g ** 2)


SHAPES: Dict[str, Callable[..., np.ndarray]] = {
    "lorentz": shape_lorentz,
    "cusp": shape_cusp,
    "asym": shape_asym,
}


class PeakDOS:
    """Nbar(eps) = N(eps)/N(0)，峰在 eps = -delta。duck-typed 兼容 `Solver.dos`。

    为什么不用 `DOS.from_table`：(i) 解析求值没有插值误差（峰宽可以低到 ~10 meV，
    而积分范围要覆盖 +-20*w_cut = +-11 eV，表格要同时满足两端很浪费）；(ii)
    `from_table` 不做归一化且会对 Nbar(0)!=1 报警，而我们要的恰恰是显式的 n0
    归一化 + lam_in 补偿（见模块 docstring 决策 1）。

    求解器只用到 dos 的 `is_constant`、`kind`（仅 kind=="table" 时读 `params`）与
    `nbar(eps)`，故 duck typing 是安全的、无需改动求解器。
    """

    kind = "p4_peak"
    params: dict = {}

    def __init__(self, a: float, gamma: float, delta: float,
                 shape: str = "lorentz", **shape_kw):
        if gamma <= 0:
            raise ValueError("gamma 必须 > 0")
        if a < 0:
            raise ValueError("a 必须 >= 0（a=0 即平坦 DOS 参考态）")
        self.a = float(a)
        self.gamma = float(gamma)
        self.delta = float(delta)
        self.shape = shape
        self.shape_kw = shape_kw
        self._s = SHAPES[shape]
        # n0 = N(0)/N_b：xi=0 处 u = 0+delta = delta
        self.n0 = 1.0 + self.a * float(self._s(np.array(float(delta)), self.gamma, **shape_kw))

    @property
    def is_constant(self) -> bool:
        return False

    def ratio_over_nb(self, xi) -> np.ndarray:
        """N(xi)/N_b，P4 (1) 式（未做 N(0) 归一化）。"""
        xi = np.asarray(xi, dtype=float)
        return 1.0 + self.a * self._s(xi + self.delta, self.gamma, **self.shape_kw)

    def nbar(self, eps) -> np.ndarray:
        """N(eps)/N(0) = [N/N_b] / n0，求解器要求的归一化。"""
        return self.ratio_over_nb(eps) / self.n0


# ===========================================================================
# 解析窗口平均
# ===========================================================================

def ntilde_lorentz(w: float, a: float, gamma: float, delta: float) -> float:
    """P4 (3) 式闭式：Ntilde(W) = 1 + a*Gam*(Gam+W)/(dlt^2+(Gam+W)^2)。"""
    weff = gamma + w
    return 1.0 + a * gamma * weff / (delta ** 2 + weff ** 2)


def e_analytic(delta, a: float, gamma: float, w_eff: float):
    """P4 (4) 式：E(dlt) = 1 + a*Gam*W_eff/(dlt^2+W_eff^2)（W_eff 已含 Gam）。"""
    delta = np.asarray(delta, dtype=float)
    return 1.0 + a * gamma * w_eff / (delta ** 2 + w_eff ** 2)


def de_ddelta_analytic(delta, a: float, gamma: float, w_eff: float):
    """P4 (6) 式：dE/ddlt = -2*dlt*(E-1)/(dlt^2+W_eff^2)。"""
    delta = np.asarray(delta, dtype=float)
    e_m1 = a * gamma * w_eff / (delta ** 2 + w_eff ** 2)
    return -2.0 * delta * e_m1 / (delta ** 2 + w_eff ** 2)


_NT_NODES: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}


def ntilde_numeric(w: float, dos_shape: str, a: float, gamma: float, delta,
                   n_nodes: int = 8001, **shape_kw):
    """任意峰形的窗口平均 (2) 式，用 xi = W*tan(theta) 代换后 Gauss-Legendre 求积。

    Ntilde(W) = 1 + (a/pi) int_{-pi/2}^{pi/2} dtheta s(W*tan(theta) + dlt)
    （代换把洛伦兹窗的极点消掉，被积函数光滑；积分范围是完整的 xi in (-inf, inf)，
    与解析式一致，不含求解器的 e_max 截断——见模块 docstring 决策 3）。
    """
    if n_nodes not in _NT_NODES:
        _NT_NODES[n_nodes] = np.polynomial.legendre.leggauss(n_nodes)
    u, wq = _NT_NODES[n_nodes]
    half = 0.5 * (np.pi - 1e-12)
    theta = half * u
    xi = w * np.tan(theta)
    delta_arr = np.atleast_1d(np.asarray(delta, dtype=float))
    s = SHAPES[dos_shape](xi[None, :] + delta_arr[:, None], gamma, **shape_kw)
    out = 1.0 + (a / np.pi) * half * np.sum(wq[None, :] * s, axis=1)
    return out[0] if np.isscalar(delta) or np.ndim(delta) == 0 else out


# ===========================================================================
# Tc 计算与反解
# ===========================================================================

@dataclass
class Setup:
    """一组固定的"背景"参数（不含 dlt）。"""
    lam_b: float = LAM_B_0
    mu_b: float = MU_B_0
    w_log: float = W_LOG_0
    w_cut_ratio: float = W_CUT_RATIO_0
    a: float = A_0
    gamma: float = GAM_0
    shape: str = "lorentz"
    mu_mode: str = "scaled"     # "scaled" | "fixed"
    n_eps: int = N_EPS
    oversample: float = 3.0
    shape_kw: Optional[dict] = None

    @property
    def w_cut(self) -> float:
        return self.w_cut_ratio * self.w_log

    def key(self) -> tuple:
        """决定平坦参考映射的参数组合（不含 a/gamma/dlt/峰形）。"""
        return (self.lam_b, self.mu_b, self.w_log, self.w_cut_ratio, self.n_eps, self.oversample)


def _solver(lam_in: float, mu_in: float, setup: Setup, dos) -> Solver:
    return Solver(
        Spectrum.einstein(lam=lam_in, w0=setup.w_log),
        mu_star=mu_in, w_cut=setup.w_cut, dos=dos,
        dos_n_eps=setup.n_eps, matsubara_oversample=setup.oversample,
    )


def tc_flat(setup: Setup, scale: float, mode: str) -> float:
    """平坦 DOS（a=0，走同一条能量依赖 DOS 数值路径）参考解的 Tc。

    mode="kernel": lam = lam_b*scale, mu* = mu_b*scale（整核缩放）
    mode="lambda": lam = lam_b*scale, mu* = mu_b（只缩放 lam）
    """
    dos = PeakDOS(0.0, setup.gamma, 0.0, shape=setup.shape)
    mu_in = setup.mu_b * scale if mode == "kernel" else setup.mu_b
    return _solver(setup.lam_b * scale, mu_in, setup, dos).tc()


def tc_peak(setup: Setup, delta: float) -> Tuple[float, float, bool]:
    """峰 DOS 的 Tc。返回 (Tc_K, n0, no_warning)。"""
    kw = setup.shape_kw or {}
    dos = PeakDOS(setup.a, setup.gamma, delta, shape=setup.shape, **kw)
    mu_in = setup.mu_b * dos.n0 if setup.mu_mode == "scaled" else setup.mu_b
    solver = _solver(setup.lam_b * dos.n0, mu_in, setup, dos)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        tc = solver.tc()
    return tc, dos.n0, len(caught) == 0


class TcMap:
    """平坦 DOS 参考解的 Tc(scale) 单调映射及其逆（用于反解 E_num）。

    在 (ln scale, ln Tc) 上建立 cubic spline。Tc(scale) 光滑单调，用 spline 逆
    比每个 dlt 都做一次 brentq 求根便宜一个量级（每次 tc() ~1 s），精度在
    `check_accuracy()` 里用真求根抽查（实测相对误差 <1e-5）。
    """

    def __init__(self, setup: Setup, mode: str, n_grid: int = 22,
                 s_lo: float = 0.8, s_hi: float = 4.0, verbose: bool = False):
        self.setup = setup
        self.mode = mode
        self.s_grid = np.geomspace(s_lo, s_hi, n_grid)
        self.tc_grid = np.empty_like(self.s_grid)
        for i, s in enumerate(self.s_grid):
            self.tc_grid[i] = tc_flat(setup, float(s), mode)
            if verbose:
                print(f"      map[{mode}] scale={s:.4f} -> Tc={self.tc_grid[i]:.4f} K", flush=True)
        if np.any(np.diff(self.tc_grid) <= 0):
            raise RuntimeError(f"平坦参考 Tc(scale) 非单调，无法反解：{self.tc_grid}")
        self._inv = interpolate.CubicSpline(np.log(self.tc_grid), np.log(self.s_grid))
        self._fwd = interpolate.CubicSpline(np.log(self.s_grid), np.log(self.tc_grid))

    @property
    def tc_range(self) -> Tuple[float, float]:
        return float(self.tc_grid[0]), float(self.tc_grid[-1])

    def invert(self, tc: float) -> float:
        lo, hi = self.tc_range
        if not (lo <= tc <= hi):
            raise ValueError(
                f"Tc={tc:.4f} K 超出平坦参考映射范围 [{lo:.4f}, {hi:.4f}] K"
                f"（scale in [{self.s_grid[0]:.3g}, {self.s_grid[-1]:.3g}]），请放宽 s_lo/s_hi"
            )
        return float(np.exp(self._inv(np.log(tc))))

    def dln_tc_dln_scale(self, scale: float) -> float:
        """d ln Tc / d ln scale，直接从求解器建的 spline 取导数（不用任何闭式）。"""
        return float(self._fwd.derivative()(np.log(scale)))

    def dln_tc_dscale(self, scale: float) -> float:
        """d ln Tc / d scale。mode="kernel" 时这就是链式法则 (5) 里应该用的
        d lnTc/dE——因为 mu_mode="scaled" 下窗口因子 E 同时乘 lam 与 mu*。
        """
        return self.dln_tc_dln_scale(scale) / scale

    def big_lambda(self, scale: float) -> float:
        """Lam = d ln Tc / d lam_eff（mode="lambda" 下有意义），lam_eff = lam_b*scale。

        reviewer 要求：Lam 一律由求解器解出的 Tc(lam) 数值求导得到，不套 McMillan /
        Allen-Dynes 闭式（后者在 lam~1.8-2.6 低估 Lam 达 27%-61%）。
        """
        if self.mode != "lambda":
            raise ValueError("big_lambda 只对 mode='lambda' 的映射有意义")
        return self.dln_tc_dln_scale(scale) / (self.setup.lam_b * scale)

    def check_accuracy(self, scales: Sequence[float]) -> List[float]:
        """抽查：对给定 scale 真算 Tc，再反解回 scale，返回相对误差。"""
        errs = []
        for s in scales:
            tc = tc_flat(self.setup, float(s), self.mode)
            errs.append(abs(self.invert(tc) / s - 1.0))
        return errs


# ===========================================================================
# 扫描与拟合
# ===========================================================================

def scan_delta(setup: Setup, deltas: Sequence[float],
               maps: Optional[Dict[str, TcMap]] = None,
               verbose: bool = True, tag: str = "") -> dict:
    """扫 dlt，返回 Tc 与两种约定下的 E_num。"""
    if maps is None:
        maps = {"kernel": TcMap(setup, "kernel"), "lambda": TcMap(setup, "lambda")}
    n = len(deltas)
    out = dict(delta=np.asarray(deltas, dtype=float), tc=np.empty(n), n0=np.empty(n),
               e_kernel=np.empty(n), e_lambda=np.empty(n), clean=np.empty(n, dtype=bool))
    for i, d in enumerate(deltas):
        t0 = time.time()
        tc, n0, clean = tc_peak(setup, float(d))
        out["tc"][i], out["n0"][i], out["clean"][i] = tc, n0, clean
        out["e_kernel"][i] = maps["kernel"].invert(tc)
        out["e_lambda"][i] = maps["lambda"].invert(tc)
        if verbose:
            print(f"    [{tag}] dlt={d:+8.1f} meV  Tc={tc:9.4f} K  n0={n0:.4f}  "
                  f"E_ker={out['e_kernel'][i]:.5f}  E_lam={out['e_lambda'][i]:.5f}"
                  f"{'' if clean else '  (WARN)'}   ({time.time()-t0:.1f}s)", flush=True)
    return out


def _r2(y: np.ndarray, yfit: np.ndarray) -> float:
    ss_res = float(np.sum((y - yfit) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def fit_weff_1p(delta: np.ndarray, e_num: np.ndarray, a: float, gamma: float,
                w_guess: Optional[float] = None) -> dict:
    """按 P4 (4) 式拟合，唯一自由参数 W_eff（振幅锁定为 a*Gam）。"""
    y = e_num - 1.0
    if w_guess is None:
        w_guess = gamma + 3.0 * W_LOG_0

    def resid(p):
        return a * gamma * p[0] / (delta ** 2 + p[0] ** 2) - y

    res = optimize.least_squares(resid, [w_guess], bounds=([1e-3], [1e6]))
    w_eff = float(res.x[0])
    fit = a * gamma * w_eff / (delta ** 2 + w_eff ** 2)
    return dict(w_eff=w_eff, r2=_r2(y, fit), max_abs_res=float(np.max(np.abs(y - fit))),
                rms_res=float(np.sqrt(np.mean((y - fit) ** 2))))


def fit_lorentz_3p(delta: np.ndarray, e_num: np.ndarray, a: float, gamma: float,
                   w_guess: Optional[float] = None) -> dict:
    """3 参数洛伦兹拟合 E-1 = A*W/((dlt-dlt0)^2+W^2)：检验振幅 A/(a*Gam)=1、峰位 dlt0=0。"""
    y = e_num - 1.0
    if w_guess is None:
        w_guess = gamma + 3.0 * W_LOG_0

    def resid(p):
        amp, w, d0 = p
        return amp * w / ((delta - d0) ** 2 + w ** 2) - y

    res = optimize.least_squares(
        resid, [a * gamma, w_guess, 0.0],
        bounds=([0.0, 1e-3, -1e5], [1e8, 1e6, 1e5]))
    amp, w, d0 = (float(v) for v in res.x)
    fit = amp * w / ((delta - d0) ** 2 + w ** 2)
    return dict(amp=amp, amp_ratio=amp / (a * gamma), w_eff=w, delta0=d0,
                r2=_r2(y, fit), max_abs_res=float(np.max(np.abs(y - fit))))


def fit_window_width(delta: np.ndarray, e_num: np.ndarray, setup: Setup,
                     w_guess: Optional[float] = None) -> dict:
    """通用窗口平均拟合：E_num(dlt) =? Ntilde(W; dlt)（(2) 式，任意峰形），自由参数 W。

    这是比"E 在 dlt 上是洛伦兹"更基本的命题：DOS 只经过一个宽度 W ~ c*w_log 的
    洛伦兹窗平均进入 Tc。对洛伦兹峰形二者等价（W_eff = Gam + W）；对 cusp/asym
    峰形只有本式成立、(4) 式不成立。
    """
    y = e_num - 1.0
    kw = setup.shape_kw or {}
    if w_guess is None:
        w_guess = 3.0 * setup.w_log

    def resid(p):
        nt = ntilde_numeric(max(p[0], 1e-6), setup.shape, setup.a, setup.gamma, delta, **kw)
        return (nt - 1.0) - y

    res = optimize.least_squares(resid, [w_guess], bounds=([1e-3], [1e6]))
    w = float(res.x[0])
    fit = ntilde_numeric(w, setup.shape, setup.a, setup.gamma, delta, **kw) - 1.0
    return dict(w=w, c=w / setup.w_log, r2=_r2(y, fit),
                max_abs_res=float(np.max(np.abs(y - fit))))


def fit_core_vs_tail(delta: np.ndarray, e_num: np.ndarray, a: float, gamma: float,
                     w_eff: float) -> dict:
    """把 dlt 网格切成"核心"(|dlt| <= W_eff) 与"尾部"(|dlt| >= W_eff) 分别做 1-参数拟合。

    reviewer 指出：E(dlt) 严格来说是**多个宽度不同的 Matsubara 窗之和**，不是单个
    洛伦兹；等效宽度在 dlt=0 与尾部会不同（玩具模型 ~18%）。这个函数直接量化
    "单窗口近似"的误差：W_core / W_tail 偏离 1 的幅度。
    """
    m_core = np.abs(delta) <= w_eff
    m_tail = np.abs(delta) >= w_eff
    out = dict(n_core=int(m_core.sum()), n_tail=int(m_tail.sum()),
               w_core=float("nan"), w_tail=float("nan"), ratio=float("nan"))
    if m_core.sum() >= 2:
        out["w_core"] = fit_weff_1p(delta[m_core], e_num[m_core], a, gamma, w_eff)["w_eff"]
    if m_tail.sum() >= 2:
        out["w_tail"] = fit_weff_1p(delta[m_tail], e_num[m_tail], a, gamma, w_eff)["w_eff"]
    if np.isfinite(out["w_core"]) and np.isfinite(out["w_tail"]) and out["w_tail"] > 0:
        out["ratio"] = out["w_core"] / out["w_tail"]
    return out


def log_sensitivities(setup: Setup, delta0: float, maps: Dict[str, TcMap],
                      w_eff: float, rel_eps: float = 0.06) -> dict:
    """在 (dlt0, a, Gam) 处中心差分测 d lnTc/d(dlt), d lnTc/d a, d lnTc/d Gam。

    **所有 Tc 都由求解器直接解全能量依赖 DOS 的 Eliashberg 方程得到**，不套用任何
    McMillan / Allen-Dynes 闭式（reviewer 要求 2）。Lam = dlnTc/dlam_eff 同样取自
    求解器建的平坦 DOS Tc(lam) 映射的数值导数。

    返回值同时给出解析链式法则 (5) 的对照：dlnTc/dX =? Lam*lam_b*dE/dX，
    dE/dX 由 (4) 式对 a / Gam / dlt 求偏导（W_eff = Gam + c*w_log，dW_eff/dGam = 1）。
    """
    h_d = max(0.03 * w_eff, 2.0)
    h_a = rel_eps * setup.a
    h_g = rel_eps * setup.gamma

    def ln_tc(st: Setup, d: float) -> float:
        return float(np.log(tc_peak(st, d)[0]))

    from dataclasses import replace as _replace
    ln_tc0 = ln_tc(setup, delta0)
    d_ln_d = (ln_tc(setup, delta0 + h_d) - ln_tc(setup, delta0 - h_d)) / (2 * h_d)
    d_ln_a = (ln_tc(_replace(setup, a=setup.a + h_a), delta0)
              - ln_tc(_replace(setup, a=setup.a - h_a), delta0)) / (2 * h_a)
    d_ln_g = (ln_tc(_replace(setup, gamma=setup.gamma + h_g), delta0)
              - ln_tc(_replace(setup, gamma=setup.gamma - h_g), delta0)) / (2 * h_g)

    # 解析对照：dE/dX（(4) 式）
    a_, g_, w_, d_ = setup.a, setup.gamma, w_eff, delta0
    den = d_ ** 2 + w_ ** 2
    e_m1 = a_ * g_ * w_ / den
    dE_dd = -2.0 * d_ * e_m1 / den
    dE_da = g_ * w_ / den
    dE_dg = a_ * (w_ / den + g_ * (d_ ** 2 - w_ ** 2) / den ** 2)  # dW_eff/dGam = 1
    tc0 = float(np.exp(ln_tc0))
    e_lam = maps["lambda"].invert(tc0)
    big_lam = maps["lambda"].big_lambda(e_lam)
    # 链式法则 (5) 里应该用哪个"dlnTc/dE"：mu_mode="scaled" 下 E 同时乘 lam 与 mu*，
    # 所以要用 kernel 映射的斜率；mu_mode="fixed" 下才是 Lam*lam_b（只缩放 lam）。
    key = "kernel" if setup.mu_mode == "scaled" else "lambda"
    e_ref = maps[key].invert(tc0)
    pre = maps[key].dln_tc_dscale(e_ref)
    pre_naive = big_lam * setup.lam_b   # P4 (5) 字面写法：Lam * lam_b
    return dict(
        delta0=delta0, tc=tc0, E_lambda=e_lam, Lambda=big_lam,
        dlnTc_dE=pre, dlnTc_dE_naive=pre_naive, E_ref=e_ref, ref_map=key,
        h_delta=h_d, h_a=h_a, h_gamma=h_g,
        dlnTc_ddelta=d_ln_d, dlnTc_da=d_ln_a, dlnTc_dgamma=d_ln_g,
        dlnTc_dlndelta=d_ln_d * delta0, dlnTc_dlna=d_ln_a * setup.a,
        dlnTc_dlngamma=d_ln_g * setup.gamma,
        pred_ddelta=pre * dE_dd, pred_da=pre * dE_da, pred_dgamma=pre * dE_dg,
        pred_naive_ddelta=pre_naive * dE_dd,
        step_delta_1pct=0.01 / abs(d_ln_d) if d_ln_d != 0 else float("inf"),
        step_a_1pct=0.01 / abs(d_ln_a) if d_ln_a != 0 else float("inf"),
        step_gamma_1pct=0.01 / abs(d_ln_g) if d_ln_g != 0 else float("inf"),
    )


def ntilde_multi(a_m: np.ndarray, w_m: np.ndarray, setup: Setup, delta) -> np.ndarray:
    """**多窗口**预言：E(dlt) = 1 + sum_m w_m [Ntilde(a_m; dlt) - 1]。

    P4 (4) 式把"每个 Matsubara 项 m 各自带一个宽度 |wt_m|=|w_m|Z_m 的窗"这件事
    近似成"一个有效窗 W_eff = Gam + c*w_log"。真实的 Tc 方程是一族宽度从
    pi*T*Z 一直排到 Matsubara 截断的窗的**加权和**——加权和的尾巴由最宽的窗主导，
    因此比单洛伦兹重。本函数按线性响应（一阶微扰，权重 w_m 见 flat_window_weights）
    给出这个加权和，用来判断 (4) 式的残差是不是就来自这个多宽度展开。
    """
    delta = np.atleast_1d(np.asarray(delta, dtype=float))
    if setup.shape == "lorentz":
        weff = setup.gamma + a_m[None, :]
        nt_m1 = setup.a * setup.gamma * weff / (delta[:, None] ** 2 + weff ** 2)
    else:
        kw = setup.shape_kw or {}
        nt_m1 = np.stack([
            ntilde_numeric(float(am), setup.shape, setup.a, setup.gamma, delta,
                           n_nodes=2001, **kw) - 1.0
            for am in a_m], axis=1)
    return 1.0 + nt_m1 @ w_m


def flat_window_weights(setup: Setup, scale: float = 1.0) -> dict:
    """平坦 DOS 参考解在其 Tc 处的线性化 phi 核的主本征模，给出各 Matsubara 列的
    权重 w_m（sum_m w_m = 1）与对应窗宽 a_m = |w_m| Z_m。

    推导：Tc 由 rho[A] = 1 定出，A[n,m] = piT(lam(n-m)-mu*th)*G_m。能量依赖 DOS 把
    第 m 列乘上 Ntilde(a_m)（见模块 docstring 决策 1），故一阶微扰
        drho = sum_{n,m} u_n A[n,m] v_m (Ntilde_m - 1) = sum_m w_m (Ntilde_m - 1)
    （u、v 为左右主本征向量，u.v=1）。整体标量缩放 E 给 drho = E-1，于是
    E ~= sum_m w_m Ntilde(a_m)，即 ntilde_multi。
    """
    dos = PeakDOS(0.0, setup.gamma, 0.0)
    solver = _solver(setup.lam_b * scale, setup.mu_b * scale, setup, dos)
    tc = solver.tc()
    t_mev = KB_MEV_PER_K * tc
    n_pos = solver._n_for_T(t_mev)
    n_arr, wn = solver._full_axis_grid(t_mev, n_pos)
    lam_full = solver._lambda_toeplitz_full(t_mev, n_pos)
    lam_z, lam_phi = solver._apply_vertex_full(lam_full, n_arr, wn, t_mev)
    e_max = solver._dos_e_max()
    z, chi, _ = solver._solve_z_chi_dos(t_mev, lam_z, lam_full, wn, e_max)
    g_int, _ = solver._dos_energy_integrals(z, chi, np.zeros_like(wn), wn, e_max)
    theta = solver._theta_mask_full(wn, t_mev)
    kernel = (np.pi * t_mev) * (lam_phi - solver.mu_star * theta[None, :]) * g_int[None, :]
    vals, vecs = np.linalg.eig(kernel)
    i = int(np.argmax(vals.real))
    rho = float(vals[i].real)
    v = vecs[:, i].real
    lvals, lvecs = np.linalg.eig(kernel.T)
    j = int(np.argmax(lvals.real))
    u = lvecs[:, j].real
    u = u / (u @ v)
    w_m = ((u[:, None] * kernel * v[None, :]).sum(axis=0)) / rho
    a_m = np.abs(wn) * z
    wm_pos = np.clip(w_m, 0.0, None)
    mean_a = float(np.sum(wm_pos * a_m) / np.sum(wm_pos))
    return dict(tc=tc, rho=rho, a_m=a_m, w_m=w_m, z=z, wn=wn, t_mev=t_mev,
                mean_width=mean_a, c_mean=mean_a / setup.w_log,
                w_sum=float(w_m.sum()), w_min=float(w_m.min()))


def window_identity_check(setup: Setup, delta: float) -> dict:
    """检验 (2) 式的窗口平均恒等式在求解器里是不是**逐项精确**成立。

    求解器的能量积分基元是 G_m = (1/pi) int deps Nbar(eps)/[(w_m Z_m)^2 + (eps+chi_m)^2]。
    对洛伦兹 DOS 这个积分有闭式（两个洛伦兹的卷积）：

        G_m = Ntilde(a_m; dlt - chi_m) / (n0 * a_m),   a_m = |w_m| Z_m

    注意 **chi_m 把窗口中心从 eps=0 挪到 eps=-chi_m**，等价于第 m 项看到的有效
    vHs 距离是 dlt - chi_m，而不是 dlt——P4 §2 的推导里没有 chi。本函数同时给出
    忽略 chi（用 dlt）与包含 chi（用 dlt-chi_m）两种预言与求解器实际 G_m 的偏差，
    以判断 (i) 窗口平均图像是否精确、(ii) chi 的位移有多大。
    """
    dos = PeakDOS(setup.a, setup.gamma, delta, shape=setup.shape, **(setup.shape_kw or {}))
    mu_in = setup.mu_b * dos.n0 if setup.mu_mode == "scaled" else setup.mu_b
    solver = _solver(setup.lam_b * dos.n0, mu_in, setup, dos)
    tc = solver.tc()
    t_mev = KB_MEV_PER_K * tc
    n_pos = solver._n_for_T(t_mev)
    n_arr, wn = solver._full_axis_grid(t_mev, n_pos)
    lam_full = solver._lambda_toeplitz_full(t_mev, n_pos)
    lam_z, _ = solver._apply_vertex_full(lam_full, n_arr, wn, t_mev)
    e_max = solver._dos_e_max()
    z, chi, conv = solver._solve_z_chi_dos(t_mev, lam_z, lam_full, wn, e_max)
    g_int, _ = solver._dos_energy_integrals(z, chi, np.zeros_like(wn), wn, e_max)
    a_m = np.abs(wn) * z
    w_eff_m = setup.gamma + a_m

    def _g(d_eff):
        """闭式 G_m。**必须带 e_max 截断修正**：求解器的 eps 积分只到 +-e_max
        (=20*w_cut)，而最高的几个 Matsubara 频率满足 a_m ~ e_max，此时"无穷积分"
        的闭式会高估 30%+（背景项 (1/pi)int_{-e}^{e} deps/(a^2+eps^2) 只有
        (1/pi a)[arctan((e-chi)/a)+arctan((e+chi)/a)]，不是 1/a）。峰项局域在
        |eps| ~ |dlt|+几倍 Gam << e_max，其截断可忽略。
        """
        bg = (np.arctan((e_max - chi) / a_m) + np.arctan((e_max + chi) / a_m)) / (np.pi * a_m)
        peak = (setup.a * setup.gamma / a_m) * w_eff_m / (d_eff ** 2 + w_eff_m ** 2)
        return (bg + peak) / dos.n0

    g_chi0 = _g(np.full_like(a_m, float(delta)))
    g_chi = _g(delta - chi)
    return dict(tc=tc, converged=conv, max_abs_chi=float(np.max(np.abs(chi))),
                chi_over_delta=float(np.max(np.abs(chi)) / abs(delta)) if delta else float("nan"),
                rel_err_no_chi=float(np.max(np.abs(g_chi0 / g_int - 1.0))),
                rel_err_with_chi=float(np.max(np.abs(g_chi / g_int - 1.0))),
                a_m=a_m, chi=chi, g_int=g_int)


def width_diagnostics(delta: np.ndarray, e_num: np.ndarray, a: float, gamma: float) -> dict:
    """从 E_num(dlt) 用三种互相独立的方式读出"洛伦兹宽度"。

    若 (4) 式 E-1 = a*Gam*W/(dlt^2+W^2) 严格成立，三者必须相等：
      W_amp  : 由峰值定 —— W = a*Gam/(E(0)-1)（(4) 式在 dlt=0）
      W_curv : 由 dlt=0 附近的曲率定 —— W = dlt/sqrt[(E(0)-1)/(E(dlt)-1) - 1]
               （取最小的非零 |dlt|）
      W_half : 由半高宽定 —— (E-1) 降到峰值一半处的 |dlt|（线性插值）
    三者的离散程度是"单洛伦兹形式失效"的最直观、与拟合权重无关的量化。
    """
    d = np.asarray(delta, dtype=float)
    y = np.asarray(e_num, dtype=float) - 1.0
    i0 = int(np.argmin(np.abs(d)))
    y0 = float(y[i0])
    out = dict(y0=y0, w_amp=a * gamma / y0 if y0 > 0 else np.nan,
               w_curv=np.nan, w_half=np.nan)
    nz = np.argsort(np.abs(d))
    for i in nz:
        if abs(d[i]) > 1e-9 and y[i] > 0 and y0 > y[i]:
            out["w_curv"] = float(abs(d[i]) / np.sqrt(y0 / y[i] - 1.0))
            out["w_curv_at_delta"] = float(abs(d[i]))
            break
    pos = np.argsort(np.abs(d))
    dd = np.abs(d[pos])
    yy = y[pos]
    for k in range(1, len(dd)):
        if yy[k] <= 0.5 * y0 <= yy[k - 1] and yy[k - 1] > yy[k]:
            t = (yy[k - 1] - 0.5 * y0) / (yy[k - 1] - yy[k])
            out["w_half"] = float(dd[k - 1] + t * (dd[k] - dd[k - 1]))
            break
    if np.isfinite(out["w_amp"]) and np.isfinite(out["w_curv"]):
        out["curv_over_amp"] = out["w_curv"] / out["w_amp"]
    else:
        out["curv_over_amp"] = np.nan
    if np.isfinite(out["w_amp"]) and np.isfinite(out["w_half"]):
        out["half_over_amp"] = out["w_half"] / out["w_amp"]
    else:
        out["half_over_amp"] = np.nan
    return out


def calibrate_c(setup: Setup, maps: Optional[Dict[str, TcMap]] = None,
                n_delta: int = 7, verbose: bool = True, tag: str = "") -> dict:
    """标定一组参数下的 c。

    对称峰形（lorentz / cusp）只在 dlt >= 0 上真算 Tc，另一半按严格偶性镜像补上
    （偶性已由主扫描数值确认到 ~1e-14，见 main() §2；它是 xi -> -xi, dlt -> -dlt,
    chi -> -chi 这个精确对称性的推论）。这样 3-参数拟合的峰位 dlt0 才有意义——
    单边数据拟合 3 个参数时 dlt0 是欠约束的。asym 峰形不对称，两侧都真算。
    """
    w_guess = setup.gamma + 3.0 * setup.w_log
    symmetric = setup.shape in ("lorentz", "cusp")
    if symmetric:
        frac = np.array([0.0, 0.35, 0.7, 1.1, 1.6, 2.3, 3.2])[:n_delta]
    else:
        frac = np.array([-2.5, -1.5, -0.75, -0.3, 0.0, 0.3, 0.75, 1.5, 2.5])
    deltas = frac * w_guess
    sc = scan_delta(setup, deltas, maps=maps, verbose=verbose, tag=tag)
    if symmetric:
        m = sc["delta"] > 0
        fit_d = np.concatenate([-sc["delta"][m][::-1], sc["delta"]])
        fit_e = np.concatenate([sc["e_kernel"][m][::-1], sc["e_kernel"]])
    else:
        fit_d, fit_e = sc["delta"], sc["e_kernel"]
    f1 = fit_weff_1p(fit_d, fit_e, setup.a, setup.gamma, w_guess)
    f3 = fit_lorentz_3p(fit_d, fit_e, setup.a, setup.gamma, w_guess)
    fw = fit_window_width(fit_d, fit_e, setup, 3.0 * setup.w_log)
    c = (f1["w_eff"] - setup.gamma) / setup.w_log
    core = np.abs(fit_d) <= f1["w_eff"]
    f1_core = (fit_weff_1p(fit_d[core], fit_e[core], setup.a, setup.gamma, f1["w_eff"])
               if core.sum() >= 2 else dict(w_eff=np.nan, r2=np.nan))
    ct = fit_core_vs_tail(fit_d, fit_e, setup.a, setup.gamma, f1["w_eff"])
    # 单点标定（只用 dlt=0）：W_eff = a*Gam/(E(0)-1)
    i0 = int(np.argmin(np.abs(sc["delta"])))
    w_1pt = setup.a * setup.gamma / (sc["e_kernel"][i0] - 1.0) if sc["delta"][i0] == 0 else np.nan
    # 多窗口（零自由参数）预言：展开点取该扫描的 E 中位数
    ww = flat_window_weights(setup, scale=float(np.median(sc["e_kernel"])))
    e_multi = ntilde_multi(ww["a_m"], ww["w_m"], setup, fit_d)
    r2_multi = _r2(fit_e - 1.0, e_multi - 1.0)
    wd = width_diagnostics(fit_d, fit_e, setup.a, setup.gamma)
    return dict(scan=sc, fit1=f1, fit3=f3, fitw=fw, c=c, core=ct, fit1_core=f1_core,
                weights=ww, r2_multi=r2_multi, widths=wd,
                c_1pt=(w_1pt - setup.gamma) / setup.w_log if np.isfinite(w_1pt) else np.nan,
                setup=setup, tag=tag)


# ===========================================================================
# CSV 输出
# ===========================================================================

def write_csv(path: Path, rows: List[dict], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(fieldnames))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"  -> {path}  ({len(rows)} 行)")


E_FIELDS = ["tag", "shape", "mu_mode", "a", "gamma_meV", "omega_log_meV", "lam_b",
            "mu_star_b", "w_cut_meV", "delta_meV", "Tc_K", "n0", "lam_in",
            "E_num", "E_num_lambda", "E_analytic", "E_window_num", "E_multiwindow",
            "converged"]


def e_rows(sc: dict, setup: Setup, tag: str, w_eff: float, w_window: float,
           e_multi: Optional[np.ndarray] = None) -> List[dict]:
    ea = e_analytic(sc["delta"], setup.a, setup.gamma, w_eff)
    kw = setup.shape_kw or {}
    ew = ntilde_numeric(w_window, setup.shape, setup.a, setup.gamma, sc["delta"], **kw)
    em = np.full(len(sc["delta"]), np.nan) if e_multi is None else np.asarray(e_multi)
    rows = []
    for i, d in enumerate(sc["delta"]):
        rows.append(dict(
            tag=tag, shape=setup.shape, mu_mode=setup.mu_mode, a=setup.a,
            gamma_meV=setup.gamma, omega_log_meV=setup.w_log, lam_b=setup.lam_b,
            mu_star_b=setup.mu_b, w_cut_meV=setup.w_cut,
            delta_meV=f"{d:.6g}", Tc_K=f"{sc['tc'][i]:.6f}", n0=f"{sc['n0'][i]:.6f}",
            lam_in=f"{setup.lam_b * sc['n0'][i]:.6f}",
            E_num=f"{sc['e_kernel'][i]:.6f}", E_num_lambda=f"{sc['e_lambda'][i]:.6f}",
            E_analytic=f"{float(np.atleast_1d(ea)[i]):.6f}",
            E_window_num=f"{float(np.atleast_1d(ew)[i]):.6f}",
            E_multiwindow=f"{float(em[i]):.6f}",
            converged=int(bool(sc["clean"][i])),
        ))
    return rows


C_FIELDS = ["tag", "shape", "mu_mode", "gamma_meV", "omega_log_meV", "lam_b", "a",
            "mu_star_b", "w_cut_meV", "gamma_over_wlog", "W_eff_meV", "c", "c_1pt",
            "R2_1p", "rms_res_1p", "W_core_meV", "W_tail_meV", "W_core_over_W_tail",
            "R2_1p_core", "W_amp_meV", "W_curv_meV", "W_half_meV", "W_curv_over_W_amp",
            "W_half_over_W_amp", "W_eff_3p_meV", "delta0_3p_meV", "amp_ratio_3p",
            "R2_3p", "W_window_meV", "c_window", "R2_window",
            "c_mean_matsubara", "R2_multiwindow", "n_delta"]


def c_row(res: dict) -> dict:
    s: Setup = res["setup"]
    ct = res.get("core") or {}
    fc = res.get("fit1_core") or {}
    ww = res.get("weights") or {}
    wd = res.get("widths") or {}
    nan = float("nan")
    return dict(
        tag=res["tag"], shape=s.shape, mu_mode=s.mu_mode, gamma_meV=s.gamma,
        omega_log_meV=s.w_log, lam_b=s.lam_b, a=s.a, mu_star_b=s.mu_b, w_cut_meV=s.w_cut,
        gamma_over_wlog=f"{s.gamma / s.w_log:.4f}",
        W_eff_meV=f"{res['fit1']['w_eff']:.4f}", c=f"{res['c']:.4f}",
        c_1pt=f"{res['c_1pt']:.4f}", R2_1p=f"{res['fit1']['r2']:.8f}",
        rms_res_1p=f"{res['fit1']['rms_res']:.3e}",
        W_core_meV=f"{ct.get('w_core', nan):.4f}", W_tail_meV=f"{ct.get('w_tail', nan):.4f}",
        W_core_over_W_tail=f"{ct.get('ratio', nan):.5f}",
        R2_1p_core=f"{fc.get('r2', nan):.8f}",
        W_amp_meV=f"{wd.get('w_amp', nan):.4f}", W_curv_meV=f"{wd.get('w_curv', nan):.4f}",
        W_half_meV=f"{wd.get('w_half', nan):.4f}",
        W_curv_over_W_amp=f"{wd.get('curv_over_amp', nan):.5f}",
        W_half_over_W_amp=f"{wd.get('half_over_amp', nan):.5f}",
        W_eff_3p_meV=f"{res['fit3']['w_eff']:.4f}", delta0_3p_meV=f"{res['fit3']['delta0']:.4f}",
        amp_ratio_3p=f"{res['fit3']['amp_ratio']:.5f}", R2_3p=f"{res['fit3']['r2']:.8f}",
        W_window_meV=f"{res['fitw']['w']:.4f}", c_window=f"{res['fitw']['c']:.4f}",
        R2_window=f"{res['fitw']['r2']:.8f}",
        c_mean_matsubara=f"{ww.get('c_mean', nan):.4f}",
        R2_multiwindow=f"{res.get('r2_multi', nan):.8f}",
        n_delta=len(res["scan"]["delta"]),
    )


SENS_FIELDS = ["delta0_meV", "Tc_K", "lam_eff", "Lambda", "dlnTc_dE", "dlnTc_dE_naive",
               "dlnTc_ddelta_per_meV", "dlnTc_da", "dlnTc_dgamma_per_meV",
               "pred_ddelta_per_meV", "pred_da", "pred_dgamma_per_meV",
               "dlnTc_dln_delta", "dlnTc_dln_a", "dlnTc_dln_gamma",
               "delta_for_1pct_meV", "a_for_1pct", "gamma_for_1pct_meV",
               "h_delta_meV", "h_a", "h_gamma_meV",
               "a", "gamma_meV", "omega_log_meV", "lam_b", "mu_star_b", "shape", "mu_mode"]


# ===========================================================================
# main
# ===========================================================================

def main(quick: bool = False) -> int:
    t_start = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    e_all: List[dict] = []
    c_all: List[dict] = []
    diag: List[dict] = []

    base = Setup()
    print("=" * 78)
    print("P4 §5：van Hove 窗口平均增强因子 E(delta) 的数值检验")
    print("=" * 78)
    print(f"基线：w_log={base.w_log} meV, lam_b={base.lam_b}, mu*_b={base.mu_b}, "
          f"a={base.a}, Gam={base.gamma} meV, w_cut={base.w_cut} meV, mu_mode={base.mu_mode}")

    # ---------------------------------------------------------------- 0. 诊断
    print("\n--- 0. 数值诊断 ---")
    map_ker = TcMap(base, "kernel")
    map_lam = TcMap(base, "lambda")
    maps0 = {"kernel": map_ker, "lambda": map_lam}
    print(f"  平坦参考映射范围: kernel Tc in [{map_ker.tc_range[0]:.3f}, {map_ker.tc_range[1]:.3f}] K, "
          f"lambda Tc in [{map_lam.tc_range[0]:.3f}, {map_lam.tc_range[1]:.3f}] K")
    errs = map_ker.check_accuracy([1.13, 1.77, 2.55])
    print(f"  spline 反解抽查相对误差 (kernel): {', '.join(f'{e:.2e}' for e in errs)}")
    diag.append(dict(check="tcmap_invert_relerr_kernel", value=f"{max(errs):.3e}",
                     note="scale=1.13/1.77/2.55 真算 Tc 再反解"))

    # Matsubara 截断收敛性
    tc_a, _, _ = tc_peak(base, 0.0)
    s2 = Setup(oversample=6.0)
    tc_b, _, _ = tc_peak(s2, 0.0)
    print(f"  Matsubara 过采样 3->6 (dlt=0): Tc {tc_a:.4f} -> {tc_b:.4f} K "
          f"(相对差 {abs(tc_b/tc_a-1):.2e})")
    diag.append(dict(check="matsubara_oversample_3_to_6_reldiff", value=f"{abs(tc_b/tc_a-1):.3e}",
                     note=f"dlt=0, Tc {tc_a:.4f} vs {tc_b:.4f} K"))

    # 窗口平均数值求积 vs 闭式 (3)
    err_nt = max(abs(ntilde_numeric(w, "lorentz", base.a, base.gamma, d)
                     / ntilde_lorentz(w, base.a, base.gamma, d) - 1.0)
                 for w in (10.0, 110.0, 800.0) for d in (0.0, 137.0, 600.0))
    print(f"  (2) 式数值求积 vs (3) 式闭式最大相对误差: {err_nt:.2e}")
    diag.append(dict(check="ntilde_numeric_vs_closed_form", value=f"{err_nt:.3e}", note=""))

    # ------------------------------------------------- 1. 主 dlt 扫描（细网格）
    print("\n--- 1. 主 dlt 扫描 (Lorentz, a=2, Gam=100 meV, mu*=0.13, mu_mode=scaled) ---")
    if quick:
        d_main = np.array([-500., -300., -150., -60., 0., 60., 150., 300., 500.])
    else:
        d_main = np.concatenate([
            np.arange(-500.0, -200.0, 50.0), np.arange(-200.0, 201.0, 25.0),
            np.arange(250.0, 501.0, 50.0)])
    sc_main = scan_delta(base, d_main, maps=maps0, tag="main")
    f1 = fit_weff_1p(sc_main["delta"], sc_main["e_kernel"], base.a, base.gamma)
    f3 = fit_lorentz_3p(sc_main["delta"], sc_main["e_kernel"], base.a, base.gamma)
    fw = fit_window_width(sc_main["delta"], sc_main["e_kernel"], base)
    f1_lam = fit_weff_1p(sc_main["delta"], sc_main["e_lambda"], base.a, base.gamma)
    c_main = (f1["w_eff"] - base.gamma) / base.w_log
    print(f"\n  1-参数拟合 (4) 式: W_eff = {f1['w_eff']:.2f} meV -> c = {c_main:.3f}, "
          f"R2 = {f1['r2']:.6f}, rms残差 = {f1['rms_res']:.2e}, max|res| = {f1['max_abs_res']:.2e}")
    print(f"  3-参数拟合:        W = {f3['w_eff']:.2f} meV, dlt0 = {f3['delta0']:+.3f} meV, "
          f"A/(a*Gam) = {f3['amp_ratio']:.5f}, R2 = {f3['r2']:.6f}")
    print(f"  通用窗口拟合 (2):  W = {fw['w']:.2f} meV -> c_window = {fw['c']:.3f}, R2 = {fw['r2']:.6f}")
    print(f"  [对照] E_lambda(固定 mu*) 的 1-参数拟合: W_eff = {f1_lam['w_eff']:.2f} meV -> "
          f"c = {(f1_lam['w_eff']-base.gamma)/base.w_log:.3f}, R2 = {f1_lam['r2']:.6f}")
    ct = fit_core_vs_tail(sc_main["delta"], sc_main["e_kernel"], base.a, base.gamma, f1["w_eff"])
    print(f"  单窗口近似误差（reviewer 3）：核心 |dlt|<=W_eff 拟合 W = {ct['w_core']:.2f} meV "
          f"({ct['n_core']} 点)，尾部 |dlt|>=W_eff 拟合 W = {ct['w_tail']:.2f} meV "
          f"({ct['n_tail']} 点)，W_core/W_tail = {ct['ratio']:.4f}"
          f"（偏离 1 即单洛伦兹不足，{abs(ct['ratio']-1)*100:.1f}%）")
    wd = width_diagnostics(sc_main["delta"], sc_main["e_kernel"], base.a, base.gamma)
    print(f"  三种独立读法的宽度（(4) 式成立则三者相等）：")
    print(f"    W_amp  (由峰值 E(0)-1={wd['y0']:.4f} 定) = {wd['w_amp']:.2f} meV")
    print(f"    W_curv (由 dlt=+-{wd.get('w_curv_at_delta', float('nan')):.0f} meV 处的曲率定) "
          f"= {wd['w_curv']:.2f} meV   -> W_curv/W_amp = {wd['curv_over_amp']:.3f}")
    print(f"    W_half (由半高宽定) = {wd['w_half']:.2f} meV   "
          f"-> W_half/W_amp = {wd['half_over_amp']:.3f}")
    diag.append(dict(check="W_core_over_W_tail", value=f"{ct['ratio']:.5f}",
                     note=f"W_core={ct['w_core']:.2f} meV ({ct['n_core']}点), "
                          f"W_tail={ct['w_tail']:.2f} meV ({ct['n_tail']}点)"))
    diag.append(dict(check="width_amp_vs_curv_vs_half",
                     value=f"W_amp={wd['w_amp']:.2f},W_curv={wd['w_curv']:.2f},"
                           f"W_half={wd['w_half']:.2f}",
                     note=f"W_curv/W_amp={wd['curv_over_amp']:.4f}, "
                          f"W_half/W_amp={wd['half_over_amp']:.4f}；(4) 式成立则都=1"))
    main_res = dict(scan=sc_main, fit1=f1, fit3=f3, fitw=fw, c=c_main, core=ct, widths=wd,
                    fit1_core=fit_weff_1p(
                        sc_main["delta"][np.abs(sc_main["delta"]) <= f1["w_eff"]],
                        sc_main["e_kernel"][np.abs(sc_main["delta"]) <= f1["w_eff"]],
                        base.a, base.gamma, f1["w_eff"]),
                    c_1pt=(base.a * base.gamma / (sc_main["e_kernel"][
                        int(np.argmin(np.abs(sc_main["delta"])))] - 1.0)
                        - base.gamma) / base.w_log,
                    setup=base, tag="main")

    # 偶性 / 峰位 / dE/ddlt 符号
    print("\n--- 2. 峰位、偶性与 dE/ddlt 符号 ---")
    d, e = sc_main["delta"], sc_main["e_kernel"]
    pairs = [(i, int(np.argmin(np.abs(d + d[i])))) for i in range(len(d)) if d[i] > 0]
    pairs = [(i, j) for i, j in pairs if abs(d[j] + d[i]) < 1e-9]
    even_res = max(abs(e[i] - e[j]) for i, j in pairs) if pairs else float("nan")
    i_max = int(np.argmax(e))
    print(f"  E_num 最大值在 dlt = {d[i_max]:+.1f} meV（网格分辨率 {np.min(np.diff(d)):.0f} meV）")
    print(f"  偶性检验 max|E(+dlt)-E(-dlt)| = {even_res:.3e}（洛伦兹峰形下应为 0，"
          f"这是 chi 的奇对称性的严格推论）")
    print(f"  3-参数拟合峰位 dlt0 = {f3['delta0']:+.3f} meV")
    de = np.diff(e) / np.diff(d)
    d_mid = 0.5 * (d[1:] + d[:-1])
    bad = [(float(dm), float(s)) for dm, s in zip(d_mid, de)
           if abs(dm) > 1.0 and np.sign(s) == np.sign(dm)]
    print(f"  (6) 式符号 sign(dE/ddlt) = -sign(dlt)：{len(de)} 个有限差分区间中 "
          f"{len(bad)} 个违反" + (f"  -> {bad}" if bad else ""))
    # 严格单调性（比"洛伦兹形式"更稳健的命题）：dlt>0 侧 E 必须逐点递减，dlt<0 侧递增
    pos = d > 0
    neg = d < 0
    mono_pos = int(np.sum(np.diff(e[pos]) >= 0))
    mono_neg = int(np.sum(np.diff(e[neg]) <= 0))
    print(f"  严格单调性：dlt>0 侧 {mono_pos} 处非递减 / {int(pos.sum())-1} 个区间；"
          f"dlt<0 侧 {mono_neg} 处非递增 / {int(neg.sum())-1} 个区间")
    print(f"  （E_num 的数值噪声量级可由偶性残差估计：{even_res:.1e}）")
    diag.append(dict(check="E_even_symmetry_max_dev", value=f"{even_res:.3e}", note="Lorentz 峰形"))
    diag.append(dict(check="peak_position_delta0_meV", value=f"{f3['delta0']:.4f}", note="3-参数拟合"))
    diag.append(dict(check="dE_ddelta_sign_violations", value=str(len(bad)),
                     note=f"共 {len(de)} 个区间"))
    diag.append(dict(check="monotonicity_violations",
                     value=f"pos={mono_pos},neg={mono_neg}",
                     note=f"dlt>0 侧 {int(pos.sum())-1} 区间, dlt<0 侧 {int(neg.sum())-1} 区间"))

    # ---------------------------- 2b. Lam = dlnTc/dlam_eff（求解器数值导数）
    print("\n--- 2b. Lam = d lnTc / d lam_eff（求解器映射的数值导数，不用 McMillan/AD）---")
    print("   lam_eff   dlnTc/dlnlam    Lam=dlnTc/dlam    Lam_AD(f1f2)  偏差    Lam_McM(裸)  偏差")
    lam_rows = []
    for s_scale in (1.0, 1.2, 1.5, 1.8):
        lam_eff = base.lam_b * s_scale
        lam_solver = map_lam.big_lambda(s_scale)
        h = 0.02 * lam_eff

        def _ad(x):
            return Spectrum.einstein(lam=x, w0=base.w_log).allen_dynes(base.mu_b)

        def _mcm(x):
            den = x - base.mu_b * (1 + 0.62 * x)
            return (base.w_log / 1.2) * np.exp(-1.04 * (1 + x) / den) / KB_MEV_PER_K

        lam_ad = (np.log(_ad(lam_eff + h)) - np.log(_ad(lam_eff - h))) / (2 * h)
        lam_mcm = (np.log(_mcm(lam_eff + h)) - np.log(_mcm(lam_eff - h))) / (2 * h)
        print(f"  {lam_eff:8.3f}   {map_lam.dln_tc_dln_scale(s_scale):12.4f}   "
              f"{lam_solver:14.5f}   {lam_ad:12.5f} {lam_ad/lam_solver-1:+6.1%}   "
              f"{lam_mcm:10.5f} {lam_mcm/lam_solver-1:+6.1%}")
        lam_rows.append(dict(check=f"Lambda_at_lam_eff_{lam_eff:.2f}",
                             value=f"{lam_solver:.5f}",
                             note=f"AD(f1f2) {lam_ad:.5f} ({lam_ad/lam_solver-1:+.1%}); "
                                  f"裸 McMillan {lam_mcm:.5f} ({lam_mcm/lam_solver-1:+.1%})"))
    diag += lam_rows

    # ------------------- 2c. 多窗口预言（解释 (4) 式的残差来自哪里）
    print("\n--- 2c. 多窗口预言：E(dlt) = sum_m w_m * Ntilde(a_m; dlt)（a_m = |w_m|Z_m）---")
    ww_lo = flat_window_weights(base, scale=1.0)
    e0_num = float(sc_main["e_kernel"][int(np.argmin(np.abs(sc_main["delta"])))])
    ww_hi = flat_window_weights(base, scale=e0_num)
    for name, ww in (("lam=lam_b", ww_lo), (f"lam=lam_b*E(0)={base.lam_b*e0_num:.3f}", ww_hi)):
        a_m, w_m = ww["a_m"], ww["w_m"]
        pos = w_m > 0
        aa = a_m[pos]
        wwt = w_m[pos] / w_m[pos].sum()
        order = np.argsort(aa)
        cum = np.cumsum(wwt[order])
        med = float(aa[order][int(np.searchsorted(cum, 0.5))])
        p90 = float(aa[order][int(np.searchsorted(cum, 0.9))])
        print(f"  [{name}] Tc_flat={ww['tc']:.3f} K, sum w_m={ww['w_sum']:.4f}, "
              f"min w_m={ww['w_min']:+.4f}")
        print(f"           窗宽 a_m 的权重分布: 均值 {ww['mean_width']:.1f} meV "
              f"(= {ww['c_mean']:.3f}*w_log), 中位 {med:.1f} meV (={med/base.w_log:.3f}*w_log), "
              f"90% 分位 {p90:.1f} meV (={p90/base.w_log:.3f}*w_log)")
        diag.append(dict(check=f"window_width_dist[{name}]",
                         value=f"mean={ww['mean_width']:.2f},median={med:.2f},p90={p90:.2f}",
                         note=f"以 w_log={base.w_log} 为单位: {ww['c_mean']:.3f}/"
                              f"{med/base.w_log:.3f}/{p90/base.w_log:.3f}"))
    e_multi_lo = ntilde_multi(ww_lo["a_m"], ww_lo["w_m"], base, sc_main["delta"])
    e_multi_hi = ntilde_multi(ww_hi["a_m"], ww_hi["w_m"], base, sc_main["delta"])
    e_sing = e_analytic(sc_main["delta"], base.a, base.gamma, f1["w_eff"])
    print("\n     dlt(meV)   E_num     多窗口(lam_b)  多窗口(lam_eff)  单洛伦兹(4)")
    step = max(1, len(sc_main["delta"]) // 12)
    for i, dv in enumerate(sc_main["delta"]):
        if i % step and abs(dv) > 1e-9:
            continue
        print(f"    {dv:+9.1f}  {sc_main['e_kernel'][i]:8.5f}  {e_multi_lo[i]:12.5f}  "
              f"{e_multi_hi[i]:14.5f}  {float(np.atleast_1d(e_sing)[i]):12.5f}")
    r2_multi_lo = _r2(sc_main["e_kernel"] - 1.0, e_multi_lo - 1.0)
    r2_multi_hi = _r2(sc_main["e_kernel"] - 1.0, e_multi_hi - 1.0)
    print(f"\n  多窗口预言（零自由参数）R2: lam_b 展开点 {r2_multi_lo:.5f}, "
          f"lam_eff 展开点 {r2_multi_hi:.5f}   vs  单洛伦兹(4) 1-参数拟合 R2 {f1['r2']:.5f}")
    diag.append(dict(check="R2_multiwindow_vs_single_lorentz",
                     value=f"multi_lamb={r2_multi_lo:.5f},multi_lameff={r2_multi_hi:.5f},"
                           f"single_fit={f1['r2']:.5f}",
                     note="多窗口是零自由参数预言，单洛伦兹是 1 参数拟合"))
    main_res["weights"] = ww_hi
    main_res["r2_multi"] = r2_multi_hi
    c_all.append(c_row(main_res))
    e_all += e_rows(sc_main, base, "main", f1["w_eff"], fw["w"], e_multi=e_multi_hi)

    # ------------------ 2d. (2) 式窗口平均恒等式在求解器里是否逐项精确
    print("\n--- 2d. 窗口平均恒等式 (2) 的逐项检验（G_m =? Ntilde(a_m; dlt-chi_m)/(n0 a_m)）---")
    print("   dlt(meV)   max|chi_m|   忽略 chi 的最大相对偏差   含 chi 的最大相对偏差")
    for dv in (0.0, 150.0, 300.0):
        wi = window_identity_check(base, dv)
        print(f"  {dv:9.1f}   {wi['max_abs_chi']:9.3f} meV   "
              f"{wi['rel_err_no_chi']:22.3e}   {wi['rel_err_with_chi']:20.3e}")
        diag.append(dict(check=f"window_identity_delta{dv:.0f}",
                         value=f"no_chi={wi['rel_err_no_chi']:.3e},"
                               f"with_chi={wi['rel_err_with_chi']:.3e}",
                         note=f"max|chi_m|={wi['max_abs_chi']:.3f} meV; 剩余偏差主要来自 "
                              f"e_max={base.w_cut*20:.0f} meV 的积分截断"))


    # -------------------------------------------------- 3. mu*=0 干净对照
    print("\n--- 3. mu*=0 对照（(4) 式在 mu*=0 时才严格：窗口因子不再乘到 mu* 上）---")
    s_mu0 = Setup(mu_b=0.0)
    maps_mu0 = {"kernel": TcMap(s_mu0, "kernel"), "lambda": TcMap(s_mu0, "lambda")}
    d_mu0 = d_main if quick else np.array([-500., -350., -250., -175., -110., -55., 0.,
                                           55., 110., 175., 250., 350., 500.])
    sc_mu0 = scan_delta(s_mu0, d_mu0, maps=maps_mu0, tag="mu0")
    f1_0 = fit_weff_1p(sc_mu0["delta"], sc_mu0["e_kernel"], s_mu0.a, s_mu0.gamma)
    f3_0 = fit_lorentz_3p(sc_mu0["delta"], sc_mu0["e_kernel"], s_mu0.a, s_mu0.gamma)
    fw_0 = fit_window_width(sc_mu0["delta"], sc_mu0["e_kernel"], s_mu0)
    c_mu0 = (f1_0["w_eff"] - s_mu0.gamma) / s_mu0.w_log
    wd0 = width_diagnostics(sc_mu0["delta"], sc_mu0["e_kernel"], s_mu0.a, s_mu0.gamma)
    ct0 = fit_core_vs_tail(sc_mu0["delta"], sc_mu0["e_kernel"], s_mu0.a, s_mu0.gamma,
                           f1_0["w_eff"])
    print(f"  W_eff = {f1_0['w_eff']:.2f} meV -> c = {c_mu0:.3f}, R2 = {f1_0['r2']:.6f}, "
          f"max|res| = {f1_0['max_abs_res']:.2e}")
    print(f"  3-参数: dlt0 = {f3_0['delta0']:+.3f} meV, A/(a*Gam) = {f3_0['amp_ratio']:.5f}, "
          f"R2 = {f3_0['r2']:.6f}")
    print(f"  三种宽度: W_amp = {wd0['w_amp']:.2f}, W_curv = {wd0['w_curv']:.2f} "
          f"(比值 {wd0['curv_over_amp']:.3f}), W_half = {wd0['w_half']:.2f} "
          f"(比值 {wd0['half_over_amp']:.3f}) meV")
    diag.append(dict(check="width_amp_vs_curv_vs_half_mu0",
                     value=f"W_amp={wd0['w_amp']:.2f},W_curv={wd0['w_curv']:.2f},"
                           f"W_half={wd0['w_half']:.2f}",
                     note=f"mu*=0；W_curv/W_amp={wd0['curv_over_amp']:.4f}, "
                          f"W_half/W_amp={wd0['half_over_amp']:.4f}"))
    e_all += e_rows(sc_mu0, s_mu0, "mu0", f1_0["w_eff"], fw_0["w"])
    c_all.append(c_row(dict(scan=sc_mu0, fit1=f1_0, fit3=f3_0, fitw=fw_0, c=c_mu0,
                            core=ct0, widths=wd0,
                            c_1pt=(s_mu0.a * s_mu0.gamma / (sc_mu0["e_kernel"][
                                int(np.argmin(np.abs(sc_mu0["delta"])))] - 1.0)
                                   - s_mu0.gamma) / s_mu0.w_log,
                            setup=s_mu0, tag="mu0")))

    # ------------------------------------------- 4. c 的标定与稳定性（OAT 设计）
    print("\n--- 4. c 的标定与稳定性 ---")
    combos: List[Tuple[str, Setup]] = []
    if quick:
        gam_list = [27.5, 110.0, 440.0]
        lam_list = [1.0, 2.2]
        a_list = [1.0]
        mu_list = [0.16]
        wcr_list = []
        scale_list = []
    else:
        gam_list = [11.0, 27.5, 55.0, 110.0, 220.0, 440.0, 880.0]
        lam_list = [1.0, 1.5, 2.2, 3.0]
        a_list = [0.5, 1.0, 2.0, 3.0]
        mu_list = [0.10, 0.16]
        wcr_list = [4.0, 8.0]
        scale_list = [55.0, 220.0]
    for g in gam_list:
        combos.append((f"gam{g:g}", Setup(gamma=g)))
    for lb in lam_list:
        combos.append((f"lamb{lb:g}", Setup(lam_b=lb)))
    for aa in a_list:
        combos.append((f"a{aa:g}", Setup(a=aa)))
    for mu in mu_list:
        combos.append((f"mu{mu:g}", Setup(mu_b=mu)))
    for wcr in wcr_list:
        # w_cut 变了要按 Tolmachev 重整化 mu*，否则比较的是不同的物理
        mu_new = mu_star_at_cutoff(MU_B_0, W_CUT_RATIO_0 * W_LOG_0, wcr * W_LOG_0)
        combos.append((f"wcut{wcr:g}", Setup(w_cut_ratio=wcr, mu_b=mu_new)))
    for wl in scale_list:
        # 整体能标缩放：Gam/w_log 固定，检验 c 的标度不变性
        combos.append((f"scale_wlog{wl:g}", Setup(w_log=wl, gamma=GAM_0 * wl / W_LOG_0)))

    map_cache: Dict[tuple, Dict[str, TcMap]] = {base.key(): maps0}
    for tag, st in combos:
        k = st.key()
        if k not in map_cache:
            print(f"  [{tag}] 建平坦参考映射 ...", flush=True)
            map_cache[k] = {"kernel": TcMap(st, "kernel"), "lambda": TcMap(st, "lambda")}
        res = calibrate_c(st, maps=map_cache[k], verbose=False, tag=tag)
        c_all.append(c_row(res))
        em = ntilde_multi(res["weights"]["a_m"], res["weights"]["w_m"], st, res["scan"]["delta"])
        e_all += e_rows(res["scan"], st, tag, res["fit1"]["w_eff"], res["fitw"]["w"], e_multi=em)
        print(f"  [{tag:16s}] Gam={st.gamma:6.1f} w_log={st.w_log:5.1f} lam_b={st.lam_b:.2f} "
              f"a={st.a:.1f} mu*={st.mu_b:.3f} w_cut={st.w_cut:6.1f} | "
              f"W_eff={res['fit1']['w_eff']:8.2f}  c={res['c']:6.3f}  "
              f"c_1pt={res['c_1pt']:6.3f}  R2={res['fit1']['r2']:.5f}  "
              f"W_core/W_tail={res['core']['ratio']:.3f}  "
              f"W_curv/W_amp={res['widths']['curv_over_amp']:.3f}  "
              f"c_mean={res['weights']['c_mean']:6.3f}  R2_multi={res['r2_multi']:.5f}  "
              f"dlt0={res['fit3']['delta0']:+7.2f}  A/(aGam)={res['fit3']['amp_ratio']:.4f}",
              flush=True)

    # ------------------------------- 5. 窄峰 vs 宽峰极限（P4 §4 的两行）
    print("\n--- 5. 窄峰/宽峰极限 (P4 §4) ---")
    print("  宽峰极限 Gam>>w_log：W_eff->Gam（推迟无关），E(0)-1 -> a")
    print("  窄峰极限 Gam<<w_log：W_eff->c*w_log（与 Gam 无关），E(0)-1 -> a*Gam/(c*w_log) ~ Gam")
    print("\n  Gam/w_log   W_eff/Gam   (W_eff-Gam)/w_log=c   E(0)-1   [E(0)-1]/a(宽峰->1)   "
          "[E(0)-1]*w_log/(a*Gam)(窄峰->1/c)")
    lim_rows = []
    for r in c_all:
        if not (str(r["tag"]).startswith("gam") or r["tag"] == "main"):
            continue
        if r["shape"] != "lorentz" or r["mu_mode"] != "scaled":
            continue
        gam, wl, aa = float(r["gamma_meV"]), float(r["omega_log_meV"]), float(r["a"])
        cc, weff = float(r["c"]), float(r["W_eff_meV"])
        e0 = None
        for er in e_all:
            if er["tag"] == r["tag"] and abs(float(er["delta_meV"])) < 1e-9:
                e0 = float(er["E_num"]) - 1.0
        if e0 is None:
            continue
        print(f"  {gam/wl:9.3f}   {weff/gam:9.3f}   {cc:18.3f}   {e0:7.4f}   "
              f"{e0/aa:19.4f}   {e0*wl/(aa*gam):32.4f}")
        lim_rows.append(dict(check=f"limit_gam_over_wlog_{gam/wl:.3f}",
                             value=f"c={cc:.4f},E0m1={e0:.5f},E0m1_over_a={e0/aa:.4f},"
                                   f"E0m1_wlog_over_aGam={e0*wl/(aa*gam):.4f}",
                             note=f"Gam={gam} w_log={wl} W_eff/Gam={weff/gam:.4f}"))
    diag += lim_rows

    # ---------------------------------------- 6. 峰形敏感性（cusp / 不对称）
    print("\n--- 6. 峰形敏感性：对称 cusp 与左右不等宽（P4 §7）---")
    d_shape = (np.array([-400., -200., -80., 0., 80., 200., 400.]) if quick
               else np.array([-500., -350., -220., -140., -80., -40., 0., 40., 80.,
                              140., 220., 350., 500.]))
    for tag, st in [("cusp", Setup(shape="cusp")),
                    ("asym4", Setup(shape="asym", shape_kw=dict(asym_ratio=4.0)))]:
        sc = scan_delta(st, d_shape, maps=maps0, verbose=False, tag=tag)
        f1s = fit_weff_1p(sc["delta"], sc["e_kernel"], st.a, st.gamma)
        f3s = fit_lorentz_3p(sc["delta"], sc["e_kernel"], st.a, st.gamma)
        fws = fit_window_width(sc["delta"], sc["e_kernel"], st)
        # 通用窗口预言的峰位（把拟合到的 W 代回 (2) 式，在细网格上找极大）
        dd = np.linspace(-500.0, 500.0, 4001)
        nt = ntilde_numeric(fws["w"], st.shape, st.a, st.gamma, dd, **(st.shape_kw or {}))
        d_pred = float(dd[int(np.argmax(nt))])
        i_num = int(np.argmax(sc["e_kernel"]))
        print(f"  [{tag}] 洛伦兹式(4) 拟合 R2 = {f1s['r2']:.5f} (3参数 R2 = {f3s['r2']:.5f}, "
              f"dlt0 = {f3s['delta0']:+.2f} meV)")
        print(f"         通用窗口(2) 拟合 R2 = {fws['r2']:.5f}, W = {fws['w']:.2f} meV -> "
              f"c_window = {fws['c']:.3f}")
        print(f"         E_num 极大在网格点 dlt = {sc['delta'][i_num]:+.1f} meV；"
              f"(2) 式用同一 W 预言极大在 dlt = {d_pred:+.1f} meV")
        e_all += e_rows(sc, st, tag, f1s["w_eff"], fws["w"])
        c_all.append(c_row(dict(scan=sc, fit1=f1s, fit3=f3s, fitw=fws,
                                c=(f1s["w_eff"] - st.gamma) / st.w_log, c_1pt=np.nan,
                                setup=st, tag=tag)))
        diag.append(dict(check=f"shape_{tag}_peak_delta_meV",
                         value=f"num_grid={sc['delta'][i_num]:.1f},fit3={f3s['delta0']:.2f},"
                               f"window_pred={d_pred:.1f}",
                         note=f"R2: lorentz(4)={f1s['r2']:.5f}, window(2)={fws['r2']:.5f}"))

    # -------------------------------------------------- 7. mu_mode=fixed 对照
    print("\n--- 7. mu_mode 对照（mu*_in 固定 vs 随 n0 缩放）---")
    st_fx = Setup(mu_mode="fixed")
    d_fx = d_shape
    sc_fx = scan_delta(st_fx, d_fx, maps=maps0, verbose=False, tag="mufixed")
    f1_fx = fit_weff_1p(sc_fx["delta"], sc_fx["e_kernel"], st_fx.a, st_fx.gamma)
    fw_fx = fit_window_width(sc_fx["delta"], sc_fx["e_kernel"], st_fx)
    print(f"  mu_mode=fixed: W_eff = {f1_fx['w_eff']:.2f} meV -> "
          f"c = {(f1_fx['w_eff']-st_fx.gamma)/st_fx.w_log:.3f}, R2 = {f1_fx['r2']:.6f}")
    print(f"  （对照 mu_mode=scaled 同一 dlt 网格：c = {c_main:.3f}）")
    e_all += e_rows(sc_fx, st_fx, "mufixed", f1_fx["w_eff"], fw_fx["w"])
    c_all.append(c_row(dict(scan=sc_fx, fit1=f1_fx,
                            fit3=fit_lorentz_3p(sc_fx["delta"], sc_fx["e_kernel"],
                                                st_fx.a, st_fx.gamma),
                            fitw=fw_fx, c=(f1_fx["w_eff"] - st_fx.gamma) / st_fx.w_log,
                            c_1pt=np.nan, setup=st_fx, tag="mufixed")))

    # ------------------------------------------------------------ 输出
    print("\n--- 8. dlt / a / Gam 三个旋钮的敏感度（检验『dlt 是主旋钮』的假设）---")
    print("  所有 Tc 由求解器直接解能量依赖 DOS 的 Eliashberg 得到；Lam 取自求解器映射的数值导数。")
    sens_rows: List[dict] = []
    d0_list = [45.0, 110.0, 250.0] if not quick else [45.0]
    for d0 in d0_list:
        s = log_sensitivities(base, d0, maps0, f1["w_eff"])
        print(f"\n  基线点 dlt0 = {d0:.0f} meV (Tc = {s['tc']:.3f} K, lam_eff = "
              f"{base.lam_b*s['E_lambda']:.4f}, Lam = dlnTc/dlam_eff = {s['Lambda']:.4f})")
        print(f"    链式法则前因子 dlnTc/dE = {s['dlnTc_dE']:.4f}（用 {s['ref_map']} 映射，"
              f"E_ref = {s['E_ref']:.4f}）；P4 (5) 字面的 Lam*lam_b = {s['dlnTc_dE_naive']:.4f}，"
              f"高估 {s['dlnTc_dE_naive']/s['dlnTc_dE']-1:+.0%}")
        print(f"    d lnTc/d dlt   = {s['dlnTc_ddelta']:+.6e} /meV   "
              f"[(5)+(4) 预言 {s['pred_ddelta']:+.6e}, 比值 "
              f"{s['pred_ddelta']/s['dlnTc_ddelta']:.3f}]")
        print(f"    d lnTc/d a     = {s['dlnTc_da']:+.6e}          "
              f"[(5)+(4) 预言 {s['pred_da']:+.6e}, 比值 {s['pred_da']/s['dlnTc_da']:.3f}]")
        print(f"    d lnTc/d Gam   = {s['dlnTc_dgamma']:+.6e} /meV   "
              f"[(5)+(4) 预言 {s['pred_dgamma']:+.6e}, 比值 "
              f"{s['pred_dgamma']/s['dlnTc_dgamma']:.3f}]")
        print(f"    无量纲：dlnTc/dln(dlt) = {s['dlnTc_dlndelta']:+.4f}, "
              f"dlnTc/dln(a) = {s['dlnTc_dlna']:+.4f}, "
              f"dlnTc/dln(Gam) = {s['dlnTc_dlngamma']:+.4f}")
        print(f"    使 Tc 变 1% 所需的改变量：|Ddlt| = {s['step_delta_1pct']:.2f} meV, "
              f"|Da| = {s['step_a_1pct']:.4f} (相对 a={base.a}: "
              f"{s['step_a_1pct']/base.a:.1%}), |DGam| = {s['step_gamma_1pct']:.2f} meV "
              f"(相对 Gam={base.gamma}: {s['step_gamma_1pct']/base.gamma:.1%})")
        sens_rows.append(dict(
            delta0_meV=f"{d0:.4g}", Tc_K=f"{s['tc']:.6f}",
            lam_eff=f"{base.lam_b*s['E_lambda']:.6f}", Lambda=f"{s['Lambda']:.6f}",
            dlnTc_dE=f"{s['dlnTc_dE']:.6f}", dlnTc_dE_naive=f"{s['dlnTc_dE_naive']:.6f}",
            dlnTc_ddelta_per_meV=f"{s['dlnTc_ddelta']:.6e}",
            dlnTc_da=f"{s['dlnTc_da']:.6e}",
            dlnTc_dgamma_per_meV=f"{s['dlnTc_dgamma']:.6e}",
            pred_ddelta_per_meV=f"{s['pred_ddelta']:.6e}",
            pred_da=f"{s['pred_da']:.6e}", pred_dgamma_per_meV=f"{s['pred_dgamma']:.6e}",
            dlnTc_dln_delta=f"{s['dlnTc_dlndelta']:.6f}",
            dlnTc_dln_a=f"{s['dlnTc_dlna']:.6f}",
            dlnTc_dln_gamma=f"{s['dlnTc_dlngamma']:.6f}",
            delta_for_1pct_meV=f"{s['step_delta_1pct']:.4f}",
            a_for_1pct=f"{s['step_a_1pct']:.6f}",
            gamma_for_1pct_meV=f"{s['step_gamma_1pct']:.4f}",
            h_delta_meV=f"{s['h_delta']:.4g}", h_a=f"{s['h_a']:.4g}",
            h_gamma_meV=f"{s['h_gamma']:.4g}",
            a=base.a, gamma_meV=base.gamma, omega_log_meV=base.w_log,
            lam_b=base.lam_b, mu_star_b=base.mu_b, shape=base.shape, mu_mode=base.mu_mode))
    print("\n  说明：把上表乘以各自的 dX/dP 才是对 dTc/dP 的贡献。p4-vhs-data 尚未回传")
    print("  da/dP、dGam/dP，因此这里只给『每个旋钮走多远等价于 1% Tc』，供拿到数据后直接换算。")

    # ------------------------------------------------------------ 输出
    print("\n--- 输出 ---")
    write_csv(OUT_DIR / "E_of_delta.csv", e_all, E_FIELDS)
    write_csv(OUT_DIR / "c_calibration.csv", c_all, C_FIELDS)
    write_csv(OUT_DIR / "sensitivity_delta_a_gamma.csv", sens_rows, SENS_FIELDS)
    write_csv(OUT_DIR / "p4_diagnostics.csv", diag, ["check", "value", "note"])

    # ------------------------------------------------------------ 汇总
    print("\n" + "=" * 78)
    print("汇总")
    print("=" * 78)
    cs = [(r["tag"], float(r["c"])) for r in c_all
          if r["shape"] == "lorentz" and r["mu_mode"] == "scaled"]
    cvals = np.array([v for _, v in cs])
    print(f"  Lorentz + mu_mode=scaled 的 {len(cvals)} 组标定：c(单洛伦兹拟合) = "
          f"{cvals.mean():.3f} +- {cvals.std():.3f}（min {cvals.min():.3f}, max {cvals.max():.3f}）")
    cm = np.array([float(r["c_mean_matsubara"]) for r in c_all
                   if r["shape"] == "lorentz" and r["mu_mode"] == "scaled"
                   and np.isfinite(float(r["c_mean_matsubara"]))])
    if len(cm):
        print(f"  同一批的 c_mean = <a_m>_w / w_log（Matsubara 窗宽的权重均值，零拟合）："
              f"{cm.mean():.3f} +- {cm.std():.3f}（min {cm.min():.3f}, max {cm.max():.3f}）")
    r2s = np.array([float(r["R2_1p"]) for r in c_all if r["shape"] == "lorentz"])
    r2m = np.array([float(r["R2_multiwindow"]) for r in c_all
                    if r["shape"] == "lorentz" and np.isfinite(float(r["R2_multiwindow"]))])
    print(f"  单洛伦兹 (4) 式 1-参数拟合 R2（Lorentz 峰形，{len(r2s)} 组）：min {r2s.min():.6f}, "
          f"中位 {np.median(r2s):.6f}, max {r2s.max():.6f}")
    if len(r2m):
        print(f"  多窗口零参数预言 R2（同批 {len(r2m)} 组）：min {r2m.min():.6f}, "
              f"中位 {np.median(r2m):.6f}, max {r2m.max():.6f}")
    ratios = np.array([float(r["W_core_over_W_tail"]) for r in c_all
                       if np.isfinite(float(r["W_core_over_W_tail"]))])
    if len(ratios):
        print(f"  W_core/W_tail（单窗口近似的自洽性，1 = 单洛伦兹成立）：中位 "
              f"{np.median(ratios):.3f}, 范围 [{ratios.min():.3f}, {ratios.max():.3f}]")
    cv = np.array([float(r["W_curv_over_W_amp"]) for r in c_all
                   if np.isfinite(float(r["W_curv_over_W_amp"]))])
    if len(cv):
        print(f"  W_curv/W_amp（峰值定的宽度 vs 曲率定的宽度，1 = (4) 式成立）：中位 "
              f"{np.median(cv):.3f}, 范围 [{cv.min():.3f}, {cv.max():.3f}]")
    n_warn = sum(1 for r in e_all if int(r["converged"]) == 0)
    print(f"  非收敛警告点：{n_warn} / {len(e_all)}")
    print(f"  总耗时 {time.time()-t_start:.1f} s")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="粗网格快跑（自测用）")
    args = ap.parse_args()
    raise SystemExit(main(quick=args.quick))
