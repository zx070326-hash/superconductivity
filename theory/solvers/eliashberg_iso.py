"""
各向同性 Migdal-Eliashberg 求解器（虚轴 Matsubara，常 DOS），v0。

目的：为 theory/derivations/P1_* 与后续 P2 提供数值验证器。不是生产代码，
是可信、可测、小而透明的研究工具。规格见 theory/solvers/SPEC.md。

单位约定（全文件统一，务必遵守）：
    - 能量（频率、耦合谱 α²F 的横轴、Δ、ω_c、ω_log、ω₂ ...）一律用 meV。
    - 温度 T 对外接口一律用 K；内部方程用 T_meV = k_B[meV/K] * T[K]。
    - 质量缩放因子 s ≡ sqrt(M / M0)（M0 为参考质量），s>1 表示更重的同位素。
    - λ、μ*、s_k 等对数敏感度、同位素系数 α 均为无量纲数。

物理方程（各向同性、常 DOS，SPEC.md 给出）：
    ω_n = πT(2n+1)                                   （n 取遍全体整数）
    Z_n = 1 + (πT/ω_n) Σ_m λ(n-m)·ω_m/√(ω_m²+Δ_m²)
    Z_n Δ_n = πT Σ_m [λ(n-m) - μ*(ω_c)θ(ω_c-|ω_m|)]·Δ_m/√(ω_m²+Δ_m²)
    λ(n-m) = 2∫₀^∞ dω ω α²F(ω) / (ω² + (ω_n-ω_m)²)

数值实现用标准的"正频折叠"技巧：假设物理解满足 Δ_{-1-n}=Δ_n、Z_{-1-n}=Z_n
（s 波奇宇称配对的标准假设），把 m 从 -∞..∞ 的求和折成 m=0..N-1：

    Z_n = 1 + (πT/ω_n) Σ_{m=0}^{N-1} [λ(n-m) - λ(n+m+1)]·ω_m/E_m
    Z_n Δ_n = πT Σ_{m=0}^{N-1} { [λ(n-m)+λ(n+m+1)] - 2μ*(ω_c)θ(ω_c-ω_m) }·Δ_m/E_m

（推导见 README.md「折叠技巧推导」一节；n, m 现在只取 0..N-1，ω_n=πT(2n+1)>0。）

Tc：线性化（Δ→0，E_m→ω_m）后是标准线性本征值问题
    Z_n^(0) Δ_n = πT Σ_m K_Δ(n,m) Δ_m/ω_m
对 T 用对数二分（scipy.optimize.brentq）找"最大本征值=1"处。
"""
from __future__ import annotations

import dataclasses
import warnings
from dataclasses import dataclass, field
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
from scipy import integrate, optimize

# --------------------------------------------------------------------------
# 物理常数
# --------------------------------------------------------------------------
KB_MEV_PER_K = 0.08617333262  # 玻尔兹曼常数, meV/K


# ==========================================================================
# Spectrum: α²F(ω) 的表示、矩、模型谱生成、质量缩放
# ==========================================================================

@dataclass
class _Component:
    """α²F(ω) 的一个分量：Einstein（δ 函数，解析）或 table（数值网格）。

    kind == "einstein": α²F_i(ω) = (lam/2)·w0·δ(ω-w0)，使得 2∫α²F_i/ω dω = lam。
    kind == "table":    由 (w, a2f) 网格给出，w 严格递增且 w[0]>0。
    """
    kind: str
    lam: float = 0.0
    w0: float = 0.0
    w: Optional[np.ndarray] = None
    a2f: Optional[np.ndarray] = None
    name: str = ""

    def scaled(self, s: float) -> "_Component":
        """质量缩放：新频率 = 旧频率 / s，s = sqrt(M/M0)。本分量 λ_i 不变（见模块顶部约定）。

        table 分量的缩放是精确的（不引入额外插值误差）：定义 α²F_s(ω) ≡ α²F(ω·s)，
        若原表在网格点 w_i 处取值 a2f_i = α²F(w_i)，则新分量在网格点 w_i/s 处取值
        α²F_s(w_i/s) = α²F(w_i) = a2f_i，恰好等于旧值——只需把网格整体除以 s。
        """
        if self.kind == "einstein":
            return _Component(kind="einstein", lam=self.lam, w0=self.w0 / s, name=self.name)
        return _Component(kind="table", w=self.w / s, a2f=self.a2f, name=self.name)

    def lambda_total(self) -> float:
        if self.kind == "einstein":
            return self.lam
        return float(2.0 * integrate.simpson(self.a2f / self.w, x=self.w))

    def moment_I2(self) -> float:
        """2∫ (α²F/ω) ln ω dω，用于组装 ω_log = exp(I2_total/λ_total)。"""
        if self.kind == "einstein":
            return self.lam * np.log(self.w0)
        return float(2.0 * integrate.simpson(self.a2f / self.w * np.log(self.w), x=self.w))

    def moment_2wA(self) -> float:
        """2∫ ω α²F dω，用于组装 ω₂ = sqrt(m2_total/λ_total)。"""
        if self.kind == "einstein":
            return self.lam * self.w0 ** 2
        return float(2.0 * integrate.simpson(self.w * self.a2f, x=self.w))

    def lambda_of_nu(self, nu: np.ndarray) -> np.ndarray:
        """λ_i(ν) = 2∫ dω ω α²F_i(ω)/(ω²+ν²)，ν 为一维数组（meV），返回同长度数组。"""
        if self.kind == "einstein":
            return self.lam * self.w0 ** 2 / (self.w0 ** 2 + nu ** 2)
        integrand = (self.w * self.a2f)[None, :] / (self.w[None, :] ** 2 + nu[:, None] ** 2)
        return 2.0 * integrate.simpson(integrand, x=self.w, axis=-1)

    def w_max(self) -> float:
        return self.w0 if self.kind == "einstein" else float(self.w[-1])


@dataclass
class Spectrum:
    """α²F(ω) = Σ_i α²F_i(ω)，由若干分量（Einstein 或 table）组成。

    用 Spectrum.from_table / Spectrum.einstein / Spectrum.two_peak / Spectrum.lorentzian
    构造；用 .combine() 拼接多分量谱（用于部分同位素替换：H 投影谱 + 宿主谱）。
    """
    components: List[_Component] = field(default_factory=list)

    # ---- 构造器 ----
    @classmethod
    def from_table(cls, w: Sequence[float], a2f: Sequence[float], name: str = "table") -> "Spectrum":
        """表格谱：w（meV，严格递增，w[0]>0）+ α²F(w)。内部用梯形/Simpson 数值积分。"""
        w = np.asarray(w, dtype=float)
        a2f = np.asarray(a2f, dtype=float)
        if w.ndim != 1 or a2f.shape != w.shape:
            raise ValueError("w, a2f 必须是等长一维数组")
        if w[0] <= 0:
            raise ValueError("网格必须从 w>0 开始，避免 α²F/ω 的 1/w 奇点")
        if np.any(np.diff(w) <= 0):
            raise ValueError("w 必须严格递增")
        return cls(components=[_Component(kind="table", w=w, a2f=a2f, name=name)])

    @classmethod
    def einstein(cls, lam: float, w0: float, name: str = "einstein") -> "Spectrum":
        """单一 Einstein 模：α²F(ω) = (lam·w0/2) δ(ω-w0)。λ(ν)、各矩全部解析。"""
        if lam <= 0 or w0 <= 0:
            raise ValueError("lam, w0 必须 > 0")
        return cls(components=[_Component(kind="einstein", lam=lam, w0=w0, name=name)])

    @classmethod
    def two_peak(cls, lam1: float, w1: float, lam2: float, w2: float,
                 names: Tuple[str, str] = ("mode1", "mode2")) -> "Spectrum":
        """双峰谱 = 两个 Einstein 分量之和（λ(ν) 精确解析，不需要网格积分）。"""
        return cls(components=[
            _Component(kind="einstein", lam=lam1, w0=w1, name=names[0]),
            _Component(kind="einstein", lam=lam2, w0=w2, name=names[1]),
        ])

    @classmethod
    def lorentzian(cls, lam: float, w0: float, gamma: float,
                   n_grid: int = 4001, w_max_factor: float = 30.0,
                   name: str = "lorentzian") -> "Spectrum":
        """Lorentzian 展宽的 Einstein 模型谱（数值网格，非解析——见 README「已知局限」）。

        形状：α²F(ω) ∝ ω·[ γ/((ω-w0)²+γ²) - γ/((ω+w0)²+γ²) ]，ω>0。
        减去镜像项是为了让 α²F(ω)/ω 在 ω→0 处正则（否则纯 Lorentzian 尾巴在 0 附近不为零，
        会造成非物理的低频拖尾），最后按 2∫α²F/ω dω = lam 归一化。
        网格范围 [1e-3·w0, w0 + w_max_factor·γ]，共 n_grid 点，Simpson 数值积分/插值。
        """
        w_hi = w0 + w_max_factor * gamma
        w = np.linspace(1e-3 * w0, w_hi, n_grid)
        shape = gamma / ((w - w0) ** 2 + gamma ** 2) - gamma / ((w + w0) ** 2 + gamma ** 2)
        shape = np.clip(shape, 0.0, None) * w  # ω·(...)，保证非负
        raw_lam = 2.0 * integrate.simpson(shape / w, x=w)
        a2f = shape * (lam / raw_lam)
        return cls(components=[_Component(kind="table", w=w, a2f=a2f, name=name)])

    # ---- 组合 / 缩放 ----
    def combine(self, other: "Spectrum") -> "Spectrum":
        return Spectrum(components=list(self.components) + list(other.components))

    def scale_mass(self, s: float, which: Optional[Iterable[str]] = None) -> "Spectrum":
        """质量缩放 α²F_s(ω)=α²F(ω·s)，s=sqrt(M/M0)。

        which=None：缩放全部分量（全同位素替换，保持 λ_total 不变）。
        which=[name,...]：只缩放指定 name 的分量（部分同位素替换，H 投影谱场景）。
        """
        names = set(which) if which is not None else None
        out = []
        for c in self.components:
            if names is None or c.name in names:
                out.append(c.scaled(s))
            else:
                out.append(c)
        return Spectrum(components=out)

    # ---- 谱的矩 ----
    @property
    def w_max(self) -> float:
        return max(c.w_max() for c in self.components)

    def lambda_total(self) -> float:
        return sum(c.lambda_total() for c in self.components)

    def omega_log(self) -> float:
        """ω_log = exp[(2/λ)∫(dω/ω)α²F(ω) ln ω]（McMillan/Allen-Dynes 对数矩），单位 meV。"""
        lam = self.lambda_total()
        i2 = sum(c.moment_I2() for c in self.components)
        return float(np.exp(i2 / lam))

    def omega_2(self) -> float:
        """ω₂ = sqrt[(2/λ)∫dω ω α²F(ω)]（二阶矩），单位 meV。"""
        lam = self.lambda_total()
        m2 = sum(c.moment_2wA() for c in self.components)
        return float(np.sqrt(m2 / lam))

    def lambda_of_nu(self, nu: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """λ(ν) = Σ_i λ_i(ν)，ν 可为标量或数组（meV），返回同形状。"""
        scalar_in = np.isscalar(nu)
        nu_arr = np.atleast_1d(np.asarray(nu, dtype=float))
        total = np.zeros_like(nu_arr)
        for c in self.components:
            total = total + c.lambda_of_nu(nu_arr)
        return float(total[0]) if scalar_in else total

    def allen_dynes(self, mu_star: float) -> float:
        """Allen-Dynes (1975) 强耦合修正公式（含 f1·f2），返回 Tc（K）。

        Tc = (ω_log/1.2)·exp[-1.04(1+λ)/(λ-μ*(1+0.62λ))]·f1·f2
        Λ1 = 2.46(1+3.8μ*)，Λ2 = 1.82(1+6.3μ*)(ω₂/ω_log)
        f1 = [1+(λ/Λ1)^1.5]^(1/3)，f2 = 1 + [(ω₂/ω_log-1)λ²]/(λ²+Λ2²)

        注意：这里的 μ* 是"直接使用的数值"，调用方需自行保证它与所设想的 ω_c 约定一致
        （本函数不做 μ*(ω_c) 的一致性检查，见 README「μ* 与 ω_c 绑定的约定」）。
        """
        lam = self.lambda_total()
        wlog = self.omega_log()
        w2 = self.omega_2()
        denom = lam - mu_star * (1 + 0.62 * lam)
        if denom <= 0:
            raise ValueError(f"λ-μ*(1+0.62λ) = {denom:.4g} <= 0，McMillan/AD 公式在此参数下无意义")
        lam1 = 2.46 * (1 + 3.8 * mu_star)
        lam2 = 1.82 * (1 + 6.3 * mu_star) * (w2 / wlog)
        f1 = (1 + (lam / lam1) ** 1.5) ** (1.0 / 3.0)
        f2 = 1 + ((w2 / wlog - 1) * lam ** 2) / (lam ** 2 + lam2 ** 2)
        tc_mev = (wlog / 1.2) * np.exp(-1.04 * (1 + lam) / denom) * f1 * f2
        return float(tc_mev / KB_MEV_PER_K)


# ==========================================================================
# μ*(ω_c) 一致重整化（Tolmachev-Morel-Anderson）
# ==========================================================================

def mu_star_at_cutoff(mu_star_ref: float, w_cut_ref: float, w_cut_new: float) -> float:
    """把参考截止 w_cut_ref 处的 μ* 换算到新截止 w_cut_new。

    Tolmachev 赝势重整化群方程：1/μ*(ω_c1) - 1/μ*(ω_c2) = ln(ω_c2/ω_c1)。
    等价于 μ*(ω_c) = μ/[1+μ ln(E_F/ω_c)]（μ、E_F 与 ω_c 无关）的直接推论。

    μ*_ref = 0（无库仑排斥）是这个关系的退化极限：裸 μ=0 在任何截止下都重整化为 0，
    直接返回 0（否则 1/μ* 发散会触发除零）。
    """
    if mu_star_ref == 0.0:
        return 0.0
    inv_new = 1.0 / mu_star_ref - np.log(w_cut_new / w_cut_ref)
    if inv_new <= 0:
        raise ValueError("换算后 1/μ* <= 0（μ* 发散或变负），w_cut_new 超出该近似的适用范围")
    return float(1.0 / inv_new)


def mu_star_from_bare(mu: float, e_fermi: float, w_cut: float) -> float:
    """从裸库仑赝势 μ 和费米能 E_F 生成给定截止 ω_c 处的 μ*：μ* = μ/[1+μ ln(E_F/ω_c)]。"""
    denom = 1.0 + mu * np.log(e_fermi / w_cut)
    if denom <= 0:
        raise ValueError("1+μ ln(E_F/ω_c) <= 0，该 (μ, E_F, ω_c) 组合下 μ* 未良定义")
    return float(mu / denom)


# ==========================================================================
# Solver: Matsubara Eliashberg 方程数值解
# ==========================================================================

@dataclass
class Solver:
    """各向同性 Matsubara Eliashberg 求解器。

    参数
    ----
    spectrum : Spectrum
    mu_star  : 库仑赝势数值 μ*，与 w_cut 绑定（见 README）。
    w_cut    : Matsubara 截止 ω_c（meV），μ* 只作用于 |ω_m| <= w_cut 的 m。
    n_matsubara : 正 Matsubara 频率个数 N；None 时按 T、w_cut 自动估计（见 _n_for_T）。
    matsubara_oversample : 自动估计 N 时的安全系数（>1 表示比"刚好覆盖 w_cut"多留裕量）。
    """
    spectrum: Spectrum
    mu_star: float
    w_cut: float
    n_matsubara: Optional[int] = None
    matsubara_oversample: float = 3.0
    _n_margin: int = 30
    _n_max: int = 20000

    # ---------------- 内部：Matsubara 网格与核矩阵 ----------------

    def _n_for_T(self, t_mev: float) -> int:
        if self.n_matsubara is not None:
            return self.n_matsubara
        n = int(np.ceil(self.matsubara_oversample * self.w_cut / (2 * np.pi * t_mev))) + self._n_margin
        n = max(n, 8)
        if n > self._n_max:
            raise RuntimeError(
                f"自动估计的 Matsubara 截断数 N={n} 过大（T={t_mev / KB_MEV_PER_K:.4g} K 相对 "
                f"w_cut={self.w_cut:.4g} meV 过低，即 w_cut/T 过大，稠密本征值分解 O(N^3) 会耗尽内存/时间）。"
                "常见于很弱耦合（Tc 本身很低）却仍用较大 w_cut 的组合。请手动设置 n_matsubara 上限、"
                "降低 matsubara_oversample，或提供更靠谱的 t_bracket_k 避免二分区间过宽。"
            )
        return n

    def _lambda_table(self, t_mev: float, n: int) -> np.ndarray:
        """λ(k)，k=0..2N-1，一次向量化调用覆盖折叠核矩阵所需的全部值。"""
        k = np.arange(0, 2 * n)
        nu = 2 * np.pi * t_mev * k
        return self.spectrum.lambda_of_nu(nu)

    def _kernels(self, t_mev: float, n: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """构造折叠后的 K1（Z 方程核）、K2（Δ 方程核）、正 Matsubara 频率数组 ω_m。"""
        lam_k = self._lambda_table(t_mev, n)
        idx = np.arange(n)
        diff = np.abs(idx[:, None] - idx[None, :])
        summ = idx[:, None] + idx[None, :] + 1
        wn = np.pi * t_mev * (2 * idx + 1)
        # θ(ω_c-ω_m) 用宽度为一个 Matsubara 间距 (2πT) 的线性斜坡代替硬台阶。
        # 原因：对 T 做连续求根（brentq）时，硬台阶会在某个 ω_m 越过 ω_c 的瞬间
        # 让核矩阵产生一个有限跳变，导致 Tc(T) 出现台阶式的微小不连续；这个不连续
        # 本身量级很小，但会被同位素系数的有限差分放大成明显噪声（实测：固定比例
        # w_cut/ω_log 时，某些"巧合"比例下 α 会跳变 ~0.03）。线性斜坡在 T→0（即
        # Matsubara 间距→0）时收敛到严格阶跃，只在硬截止边界附近的一项上引入
        # O(1/N) 的修正，物理上可忽略，但换来了 Tc(T) 对 T 连续可导，二分/求根稳定。
        theta = np.clip((self.w_cut - wn) / (2 * np.pi * t_mev) + 0.5, 0.0, 1.0)
        k1 = lam_k[diff] - lam_k[summ]
        k2 = lam_k[diff] + lam_k[summ] - 2.0 * self.mu_star * theta[None, :]
        return k1, k2, wn

    def _max_eigen_linearized(self, t_mev: float, n: int) -> float:
        """线性化（Δ→0）本征值问题的最大实本征值。"""
        k1, k2, wn = self._kernels(t_mev, n)
        z0 = 1.0 + (np.pi * t_mev / wn) * k1.sum(axis=1)
        a = (np.pi * t_mev) * k2 / wn[None, :]
        a = a / z0[:, None]
        eigvals = np.linalg.eigvals(a)
        i = np.argmax(eigvals.real)
        if abs(eigvals.imag[i]) > 1e-6 * max(abs(eigvals.real[i]), 1e-12):
            warnings.warn(f"主导本征值有明显虚部 ({eigvals[i]!r})，请检查参数是否物理")
        return float(eigvals.real[i])

    # ---------------- 公开接口：Tc ----------------

    def tc(self, t_bracket_k: Optional[Tuple[float, float]] = None, xtol: float = 1e-10) -> float:
        """求 Tc（K）：线性化核矩阵最大本征值=1 处，对 ln T 用 brentq 求根。"""

        def f(ln_t_k: float) -> float:
            t_k = np.exp(ln_t_k)
            t_mev = KB_MEV_PER_K * t_k
            n = self._n_for_T(t_mev)
            return self._max_eigen_linearized(t_mev, n) - 1.0

        if t_bracket_k is None:
            # 用 Allen-Dynes 估计值作为初始区间中心，避免二分区间过宽导致
            # N 自动估计（∝ 1/T）在低 T 端爆炸。AD 估计通常与数值解同量级
            # （见 SPEC 测试 1），widening 循环兜底处理偏差较大的情形。
            try:
                t_guess = self.spectrum.allen_dynes(self.mu_star)
                t_lo, t_hi = t_guess / 8.0, t_guess * 4.0
            except ValueError:
                wlog = self.spectrum.omega_log()
                t_lo, t_hi = 1e-2 * wlog / KB_MEV_PER_K, 0.5 * wlog / KB_MEV_PER_K
        else:
            t_lo, t_hi = t_bracket_k
        lo, hi = np.log(t_lo), np.log(t_hi)
        f_lo, f_hi = f(lo), f(hi)
        tries = 0
        while f_hi > 0 and tries < 40:
            hi += np.log(2.0)
            f_hi = f(hi)
            tries += 1
        tries = 0
        while f_lo < 0 and tries < 40:
            lo -= np.log(2.0)
            f_lo = f(lo)
            tries += 1
        if f_lo < 0 or f_hi > 0:
            raise RuntimeError(
                "无法为 Tc 定位合适的二分区间（f_lo={:.4g}, f_hi={:.4g}）；"
                "请检查 λ、μ*、w_cut 是否给出了超导解".format(f_lo, f_hi)
            )
        ln_tc = optimize.brentq(f, lo, hi, xtol=xtol)
        return float(np.exp(ln_tc))

    # ---------------- 公开接口：非线性 Δ(T) ----------------

    def solve_gap(self, t_k: float, delta0: Optional[np.ndarray] = None,
                  tol: float = 1e-10, max_iter: int = 2000, mixing: float = 0.3
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
        """T<Tc 时的非线性自洽解（线性混合不动点迭代，v0 用简单稳健的方案）。

        返回 (wn, Delta, Z, converged)：wn、Delta 为 meV，长度 N 的正 Matsubara 频率数组。
        """
        t_mev = KB_MEV_PER_K * t_k
        n = self._n_for_T(t_mev)
        k1, k2, wn = self._kernels(t_mev, n)

        if delta0 is None:
            a = (np.pi * t_mev) * k2 / wn[None, :]
            z0 = 1.0 + (np.pi * t_mev / wn) * k1.sum(axis=1)
            a = a / z0[:, None]
            vals, vecs = np.linalg.eig(a)
            i = np.argmax(vals.real)
            v = vecs[:, i].real
            if v[0] < 0:
                v = -v
            delta = 0.05 * wn[0] * v / (np.max(np.abs(v)) + 1e-30)
        else:
            delta = np.array(delta0, dtype=float)

        converged = False
        for _ in range(max_iter):
            e = np.sqrt(wn ** 2 + delta ** 2)
            z = 1.0 + (np.pi * t_mev / wn) * (k1 @ (wn / e))
            rhs = (np.pi * t_mev) * (k2 @ (delta / e))
            delta_new = rhs / z
            denom = np.max(np.abs(delta_new)) + 1e-30
            diff = np.max(np.abs(delta_new - delta)) / denom
            delta = mixing * delta_new + (1 - mixing) * delta
            if diff < tol:
                converged = True
                break
        e = np.sqrt(wn ** 2 + delta ** 2)
        z = 1.0 + (np.pi * t_mev / wn) * (k1 @ (wn / e))
        return wn, delta, z, converged

    # ---------------- 公开接口：同位素系数 ----------------

    def isotope_coefficient(
        self,
        parts: Dict[str, Spectrum],
        which: Union[str, Sequence[str]],
        mu_star_fixed: bool,
        eps: float = 0.02,
        w_cut_ref_ratio: Optional[float] = None,
    ) -> Tuple[float, Dict[str, float]]:
        """同位素系数 α = -Δln Tc/Δln M（中心差分，M -> M(1±eps)）。

        parts: {分量名: Spectrum}，组合谱 = Σ parts.values()（自身的 self.spectrum 应等于
               这个组合，调用方负责保持一致；这样设计是为了支持"只缩放其中一个分量"）。
        which: "all"（全同位素替换，缩放全部分量）或分量名 / 分量名列表（部分同位素替换）。
        mu_star_fixed: P1 §3.1 的核心开关。
            True  -> μ* 数值原样保留，不随 ω_c 移动重新推导（"μ* 与 M 无关"这一理想化假设）。
            False -> ω_c 随谱整体移动（比例 = self.w_cut/self.spectrum.omega_log()，
                     可用 w_cut_ref_ratio 覆盖），μ*(ω_c) 按 Tolmachev 关系重新推导。
        eps: 中心差分步长；SPEC 建议 0.02，并检查对 eps 的稳定性（见测试）。

        返回 (alpha, diagnostics)，diagnostics 含 Tc_plus/minus、mu_star_plus/minus、
        w_cut_plus/minus，供检查用了哪个约定。
        """
        names = list(parts.keys())
        which_list = names if which == "all" else ([which] if isinstance(which, str) else list(which))
        for w in which_list:
            if w not in parts:
                raise ValueError(f"which 中的分量名 {w!r} 不在 parts 里")

        base_spectrum = Spectrum(components=[])
        for sp_i in parts.values():
            base_spectrum = base_spectrum.combine(sp_i)
        ratio = w_cut_ref_ratio if w_cut_ref_ratio is not None else self.w_cut / base_spectrum.omega_log()

        def build(s: float) -> Tuple[float, float, float]:
            spec = Spectrum(components=[])
            for name, sp_i in parts.items():
                spec = spec.combine(sp_i.scale_mass(s) if name in which_list else sp_i)
            w_cut_new = ratio * spec.omega_log()
            if mu_star_fixed:
                mu_star_new = self.mu_star
            else:
                mu_star_new = mu_star_at_cutoff(self.mu_star, self.w_cut, w_cut_new)
            solver = Solver(spec, mu_star=mu_star_new, w_cut=w_cut_new,
                             n_matsubara=self.n_matsubara,
                             matsubara_oversample=self.matsubara_oversample)
            return solver.tc(), mu_star_new, w_cut_new

        s_plus = np.sqrt(1 + eps)
        s_minus = np.sqrt(1 - eps)
        tc_plus, mu_plus, wc_plus = build(s_plus)
        tc_minus, mu_minus, wc_minus = build(s_minus)
        alpha = -(np.log(tc_plus) - np.log(tc_minus)) / (np.log(1 + eps) - np.log(1 - eps))
        diag = dict(
            Tc_plus=tc_plus, Tc_minus=tc_minus,
            mu_star_plus=mu_plus, mu_star_minus=mu_minus,
            w_cut_plus=wc_plus, w_cut_minus=wc_minus,
        )
        return float(alpha), diag

    # ---------------- 公开接口：对数敏感度 ----------------

    def log_sensitivity(self, param: str, eps: float = 0.02) -> float:
        """数值计算 ∂ln Tc/∂ln(param)（中心差分）。

        param 支持：
            "lambda"    : α²F 幅度整体缩放 (1±eps)，ω_log、ω₂、谱形状不变（纯 λ 敏感度，
                          对应 McMillan/AD 里 ∂lnΦ/∂lnλ 的定义）。
            "omega_log" : 全谱频率轴整体缩放（等价于全同位素质量缩放，但 μ*、w_cut 的
                          绑定比例保持不变——不等同于 isotope_coefficient(mu_star_fixed=False)，
                          见 README）。
            "mu_star"   : μ* 数值直接缩放，w_cut、谱不变。
            "w_cut"     : w_cut 直接缩放，μ*、谱不变（μ* 不随 w_cut 移动重新推导）。
        """
        spec = self.spectrum

        if param == "lambda":
            def make(delta: float) -> float:
                comps = []
                for c in spec.components:
                    if c.kind == "einstein":
                        comps.append(dataclasses.replace(c, lam=c.lam * (1 + delta)))
                    else:
                        comps.append(dataclasses.replace(c, a2f=c.a2f * (1 + delta)))
                new_spec = Spectrum(components=comps)
                return Solver(new_spec, self.mu_star, self.w_cut,
                              self.n_matsubara, self.matsubara_oversample).tc()

        elif param == "omega_log":
            def make(delta: float) -> float:
                s = 1.0 / (1 + delta)
                new_spec = spec.scale_mass(s)
                new_w_cut = self.w_cut * (1 + delta)
                return Solver(new_spec, self.mu_star, new_w_cut,
                              self.n_matsubara, self.matsubara_oversample).tc()

        elif param == "mu_star":
            def make(delta: float) -> float:
                return Solver(spec, self.mu_star * (1 + delta), self.w_cut,
                              self.n_matsubara, self.matsubara_oversample).tc()

        elif param == "w_cut":
            def make(delta: float) -> float:
                return Solver(spec, self.mu_star, self.w_cut * (1 + delta),
                              self.n_matsubara, self.matsubara_oversample).tc()

        else:
            raise ValueError(f"不支持的 param: {param!r}")

        tc_p, tc_m = make(eps), make(-eps)
        return float((np.log(tc_p) - np.log(tc_m)) / (np.log(1 + eps) - np.log(1 - eps)))
