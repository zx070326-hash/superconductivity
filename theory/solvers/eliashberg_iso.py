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
import functools
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
# DOS: 能量依赖态密度 N̄(ε) ≡ N(ε)/N(0)（SPEC 必备功能 6，第二阶段）
# ==========================================================================

@dataclass
class DOS:
    """归一化态密度 N̄(ε) ≡ N(ε)/N(0) 的模型，供 Mitrović–Carbotte/Pickett 推广用。

    核心不变量：**N̄(0) ≡ 1 恒成立**——这不是可选的归一化，是定义本身的直接推论
    （N(0) 就是 ε=0 处的值，比值在 ε=0 处必为 1）。各构造器的处理方式：

    - ``constant()``：N̄≡1，显然满足。
    - ``step()``/``log_vhs()``：内部先写出一个直观的"形状函数" g(ε)（形如
      "远离特征处基线为 1，特征处/另一侧变成 ratio，或 log 发散"），然后返回
      g(ε)/g(0) 强制 N̄(0)=1。这样做的物理理由：ε 是相对费米面测量的准粒子能量，
      "费米面处的 DOS" 天然就是归一化基准，不应该因为 δ（特征相对 E_F 的位置）
      取到 0 附近就让基准本身漂移。副作用：当 delta 明显偏离 0（特征远离 E_F，
      正常使用场景）时，g(0)≈基线值≈1，归一化近似是恒等变换，不影响直觉理解；
      只有 delta→0（特征恰好或接近落在 E_F）时归一化才会显著改变数值——而这恰好
      是 Tsuei (1990) 式 vHs-钉扎-费米面 场景需要精确处理的地方。
    - ``from_table()``：**不做自动归一化**。表格数据的 ε=0 点应该就是调用方计算/
      测量 N(0) 的定义点；如果表里插值到 ε=0 处不是 1，多半是数据准备错误，不应该
      被静默"修正"掉（会掩盖真实的 bug）。``from_table`` 会在明显偏离时发出警告，
      但不强制改动数据。
    """
    kind: str
    params: dict = field(default_factory=dict)

    @property
    def is_constant(self) -> bool:
        return self.kind == "constant"

    def nbar(self, eps: Union[float, np.ndarray]) -> np.ndarray:
        """N̄(ε)，ε（meV）可为标量或任意形状数组，返回同形状数组（float64）。"""
        eps = np.asarray(eps, dtype=float)
        if self.kind == "constant":
            return np.ones_like(eps)

        if self.kind == "step":
            delta, ratio, width = (self.params[k] for k in ("delta", "ratio", "width"))

            def g(e: np.ndarray) -> np.ndarray:
                return 1.0 + (ratio - 1.0) * 0.5 * (1.0 + np.tanh((e - delta) / width))

            g0 = float(g(np.array(0.0)))
            return g(eps) / g0

        if self.kind == "log_vhs":
            delta, a, w_band, cut, floor = (
                self.params[k] for k in ("delta", "a", "W", "cut", "floor")
            )

            def g(e: np.ndarray) -> np.ndarray:
                return 1.0 + a * np.log(w_band / np.maximum(np.abs(e - delta), cut))

            g0 = float(g(np.array(0.0)))
            if g0 <= 0:
                raise ValueError(
                    f"log_vhs 在 ε=0 处的形状函数 g(0)={g0:.4g} <= 0（delta/a/W/cut 组合"
                    "导致归一化分母非正），请调整参数"
                )
            if g0 < 0.1:
                warnings.warn(
                    f"log_vhs 归一化分母 g(0)={g0:.4g} 接近 0（|delta| 相对 W 过大，"
                    "即 |ln(W/|0-delta|)| 接近 1/a）：g(ε)/g(0) 会把 g(ε) 的普通量级"
                    "（~1）放大成远大于 1 的 N̄，且 g(ε) 在积分窗口内可能变号（被 floor"
                    "强行钳到下限），二者共同作用容易让自洽迭代不收敛或给出不可靠结果。"
                    "建议增大 W（使其显著大于 |delta|）或缩小 |delta|，见 README「能量"
                    "依赖 DOS」已知局限",
                    stacklevel=2,
                )
            return np.clip(g(eps) / g0, floor, None)

        if self.kind == "table":
            eps_tab, nbar_tab = self.params["eps"], self.params["nbar"]
            return np.interp(eps, eps_tab, nbar_tab)  # 表格边界外按端点值延拓（np.interp 默认行为）

        raise ValueError(f"未知 DOS kind: {self.kind!r}")

    # ---- 构造器 ----
    @classmethod
    def constant(cls) -> "DOS":
        """N̄(ε) ≡ 1（v0 的常 DOS 假设）。"""
        return cls(kind="constant")

    @classmethod
    def step(cls, delta: float, ratio: float, width: float) -> "DOS":
        """台阶型：形状函数在 ε≪δ 侧为 1、ε≫δ 侧为 ratio，用 tanh 做宽度 width 的平滑
        过渡；返回值再除以形状函数在 ε=0 处的值（见类 docstring）。

        典型用法假设 |δ| ≳ 几倍 width（台阶明显偏离费米面），此时归一化近似是恒等
        变换，"ε≫δ 处 N̄≈ratio" 的直观描述近似成立；δ 落在 0 附近时 N̄(0) 强制为 1，
        两侧的相对比例仍是 ratio，但绝对基线会相应偏移，这是设计使然，不是 bug。
        """
        if width <= 0:
            raise ValueError("width 必须 > 0")
        if ratio <= 0:
            raise ValueError("ratio 必须 > 0（态密度非负）")
        return cls(kind="step", params=dict(delta=float(delta), ratio=float(ratio), width=float(width)))

    @classmethod
    def log_vhs(cls, delta: float, a: float, W: float, cut: float, floor: float = 1e-3) -> "DOS":
        """对数型 van Hove 奇点：形状函数 g(ε) = 1 + a·ln(W/max(|ε-δ|, cut))，
        cut 正则化 ε=δ 处的对数发散；最终返回 g(ε)/g(0)（见类 docstring），
        并用 floor 兜底裁剪，避免远离 δ 处（|ε-δ|≫W）g 可能变号/趋于 -∞ 导致
        N̄ 出现非物理的负值或跳变（真实材料的 DOS 不会是负的；这个模型只在
        |ε-δ| ≲ W 量级内有意义，floor 只是数值安全网，不是物理内容，用前须知）。

        delta：vHs 相对费米面的位置（meV，可正可负，可为 0）。
        a：vHs 强度（无量纲，a>0 表示 δ 处是增强峰）。
        W：背景能量尺度（meV，通常取能带宽量级，例如 ~1 eV）。
        cut：正则化尺度（meV，通常取远小于 W、量级与 k_BTc 或 ω 的最小 Matsubara
             间距相当，避免正则化本身影响感兴趣的物理）。

        参数选取的重要约束（实测踩过的坑，见 README「能量依赖 DOS·已知局限」）：
        **|delta| 应明显小于 W**（|delta| ≲ few×W 的量级，不要取 |delta|~10W 以上）。
        归一化用的是 g(ε)/g(0)，g(0)=1+a·ln(W/|delta|)——delta 一旦显著超过 W，
        ln(W/|delta|) 变成较大的负数，g(0) 可能接近 0 甚至变负（变负会被上面的检查
        直接拒绝，接近 0 会触发 UserWarning）。g(0)→0 时哪怕 g(ε) 本身量级正常
        （~1），除出来的 N̄=g(ε)/g(0) 也会被放大到远超 1 的量级，且 g(ε) 在感兴趣的
        ε 窗口内可能变号（被 floor 钳位），这种"陡峭+钳位"的 N̄(ε) 会让能量积分的
        自洽迭代（`_solve_z_chi_dos`）难以收敛，实测在 delta=60·w_cut、W=1000 meV
        时会看到 Z、χ 发散式振荡、500 次迭代不收敛。这不是"delta 越大越接近常数 DOS
        因而无所谓"的直觉——恰恰相反，delta 相对 W 过大时模型本身的归一化先失效。
        正确做法：验证"δ≫ω_c 回到常数 DOS"这类极限时，应同时把 W 取得足够大
        （W ≳ 几倍 delta），而不是固定 W 只增大 delta。
        """
        if W <= 0 or cut <= 0:
            raise ValueError("W, cut 必须 > 0")
        return cls(kind="log_vhs", params=dict(
            delta=float(delta), a=float(a), W=float(W), cut=float(cut), floor=float(floor)
        ))

    @classmethod
    def from_table(cls, eps_meV: Sequence[float], nbar: Sequence[float]) -> "DOS":
        """表格 N̄(ε)：eps_meV 严格递增（meV，可含负值，ε=0 应在范围内）；
        nbar 是对应的 N(ε)/N(0) 数值。线性插值，边界外按端点值延拓。

        不做自动归一化（见类 docstring）；若插值到 ε=0 处显著偏离 1（|差|>1e-3），
        发出警告但不修改数据——这通常意味着表格没有以 N(0) 归一化，请检查数据来源。
        """
        eps_meV = np.asarray(eps_meV, dtype=float)
        nbar = np.asarray(nbar, dtype=float)
        if eps_meV.ndim != 1 or nbar.shape != eps_meV.shape:
            raise ValueError("eps_meV, nbar 必须是等长一维数组")
        if np.any(np.diff(eps_meV) <= 0):
            raise ValueError("eps_meV 必须严格递增")
        if eps_meV[0] > 0 or eps_meV[-1] < 0:
            raise ValueError("eps_meV 范围必须包含 ε=0（费米面）")
        n0 = float(np.interp(0.0, eps_meV, nbar))
        if abs(n0 - 1.0) > 1e-3:
            warnings.warn(
                f"DOS.from_table：插值到 ε=0 处 N̄={n0:.4g}，偏离 1 超过 1e-3；"
                "N̄(ε)≡N(ε)/N(0) 在 ε=0 处按定义应为 1，请检查表格是否已经用 N(0) 归一化"
                "（本构造器不会自动修正，见类 docstring）。"
            )
        return cls(kind="table", params=dict(eps=eps_meV, nbar=nbar))


# ==========================================================================
# 模型顶角核钩子（SPEC 必备功能 7，第二阶段）
# ==========================================================================

def vertex_pietronero_grimaldi_isotropic(m_ratio: float, Qc: float) -> Callable:
    """P2 占位顶点核：形式上对应 Pietronero-Grimaldi-Strässler 一类"超越 Migdal"
    顶点修正的**领头阶、常数近似**——不是从任何具体文献公式精确抄录的结果。

    物理动机（占位，真正的动量/频率依赖形式留给 P2 推导）：非绝热顶点修正的量级
    应正比于绝热比 m_ratio ≡ ω_log/E_F（Migdal 展开的小参数），符号/相对大小由
    声子特征动量转移与费米动量之比 Qc（例如 q_c/k_F）调制。这里用一个不依赖任何
    文献细节的简单解析插值：

        P_V(Qc) = 1 - 2·Qc        （Qc=0 时 P_V=+1，Qc=1 时 P_V=-1，Qc=0.5 时为 0）
        V ≡ m_ratio · P_V(Qc)      （常数，不依赖 n、m、ω_n、ω_m、T）

    即 Qc→0（前向散射主导，声子动量转移小）时顶点修正为正（增强配对核），
    Qc→1（大动量转移显著）时符号翻转、抑制配对，与文献里定性报道的趋势一致；
    但这只是一个占位插值，量级和函数形式都未经过 P2 推导核实，**不要用于任何
    定量结论**，只用于验证"顶点钩子"这条数值管线本身能正确工作。

    返回的 callable 签名与 ``Solver(vertex=...)`` 要求一致：
        vertex(n_idx, m_idx, omega_n, omega_m, T) -> ndarray（与输入同形状的 V）
    """
    def vertex(n_idx, m_idx, omega_n, omega_m, T):
        p_v = 1.0 - 2.0 * Qc
        val = m_ratio * p_v
        shape = np.broadcast(np.asarray(n_idx), np.asarray(m_idx),
                              np.asarray(omega_n), np.asarray(omega_m)).shape
        return np.full(shape, val, dtype=float)

    return vertex


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
    dos : 能量依赖 DOS（SPEC 功能 6，第二阶段）。默认 DOS.constant()，此时求解器
          走与 v0 完全相同的折叠代数路径（结果逐位一致，见 README「能量依赖 DOS」）。
    fix_mu : 化学势固定策略。v0.1 只实现 "none"（不做电子数守恒修正，见 README
             已知局限）；传其他值会在构造时报 NotImplementedError。
    vertex : 模型顶角核钩子（SPEC 功能 7，第二阶段）。None 时无影响（结果与 v0
             逐位一致）；否则是可调用对象
             ``vertex(n_idx, m_idx, omega_n, omega_m, T) -> ndarray``，
             返回与 λ(n−m) 同形状的乘子矩阵 V，作用方式见 README。
    vertex_in_Z : vertex 是否也作用于 Z 方程的电声部分（默认 True）；恒作用于
                  φ/Δ 方程的电声部分，从不作用于 μ* 项。
    dos_n_eps, dos_energy_cutoff_ratio, dos_max_iter, dos_mixing, dos_tol :
        非常 DOS 路径（Mitrović–Carbotte 推广）的数值参数，见 README。
    """
    spectrum: Spectrum
    mu_star: float
    w_cut: float
    n_matsubara: Optional[int] = None
    matsubara_oversample: float = 3.0
    dos: DOS = field(default_factory=DOS.constant)
    fix_mu: str = "none"
    vertex: Optional[Callable[..., np.ndarray]] = None
    vertex_in_Z: bool = True
    dos_n_eps: int = 400
    dos_energy_cutoff_ratio: float = 20.0
    dos_max_iter: int = 500
    dos_mixing: float = 0.3
    dos_tol: float = 1e-9
    _n_margin: int = 30
    _n_max: int = 20000

    def __post_init__(self) -> None:
        if self.fix_mu not in ("none",):
            raise NotImplementedError(
                f"fix_mu={self.fix_mu!r} 未实现：v0.1 只支持 'none'（忽略化学势的电子数"
                "守恒修正），见 README「能量依赖 DOS」已知局限"
            )

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
        """构造折叠后的 K1（Z 方程核）、K2（Δ 方程核）、正 Matsubara 频率数组 ω_m。

        vertex 钩子（SPEC 功能 7）：若 self.vertex 不为 None，把折叠前的两个
        λ 值 λ(n-m)、λ(n+m+1) 都乘以 (1+V(n,m,ω_n,ω_m,T))——n、m 是折叠后的
        正频指标 (0..N-1)，ω_n、ω_m 是对应的正频率；即 V 同时施加于"直接项"
        λ(n-m) 和"折叠镜像项"λ(n+m+1)，不区分二者对应的原始负频率细节
        （这是 P2 占位钩子的简化，若后续需要精确到 Nambu/正负频结构，需要重新
        推导，见 README「模型顶角核钩子」）。μ* 项永远不乘 (1+V)；vertex_in_Z
        控制 K1（Z 方程）是否也乘，K2（φ/Δ 方程电声部分）恒乘。
        vertex=None 时代码路径与 v0 完全一致（不引入任何额外乘法），逐位相同。
        """
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
        lam_diff = lam_k[diff]
        lam_summ = lam_k[summ]
        if self.vertex is not None:
            n_grid = idx[:, None] * np.ones((1, n), dtype=int)
            m_grid = idx[None, :] * np.ones((n, 1), dtype=int)
            wn_n_grid = wn[:, None] * np.ones((1, n))
            wn_m_grid = wn[None, :] * np.ones((n, 1))
            v = np.asarray(
                self.vertex(n_grid, m_grid, wn_n_grid, wn_m_grid, t_mev), dtype=float
            )
            lam_diff_phi = lam_diff * (1.0 + v)
            lam_summ_phi = lam_summ * (1.0 + v)
            if self.vertex_in_Z:
                lam_diff_z, lam_summ_z = lam_diff_phi, lam_summ_phi
            else:
                lam_diff_z, lam_summ_z = lam_diff, lam_summ
        else:
            lam_diff_phi, lam_summ_phi = lam_diff, lam_summ
            lam_diff_z, lam_summ_z = lam_diff, lam_summ
        k1 = lam_diff_z - lam_summ_z
        k2 = lam_diff_phi + lam_summ_phi - 2.0 * self.mu_star * theta[None, :]
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
        """求 Tc（K）。dos.is_constant 时走 v0 折叠代数路径（结果逐位一致于 v0，
        vertex 钩子已接入该路径，见 _kernels）；否则走能量依赖 DOS 的
        Mitrović–Carbotte 推广路径（见 README「能量依赖 DOS」）。
        """
        if self.dos.is_constant:
            return self._tc_constant_dos(t_bracket_k, xtol)
        return self._tc_energy_dependent_dos(t_bracket_k, xtol)

    def _bisect_tc(self, max_eig_fn: Callable[[float], float],
                    t_bracket_k: Optional[Tuple[float, float]], xtol: float) -> float:
        """对 ln T 用 brentq 找 max_eig_fn(T_K)=1 的根；两条 Tc 路径共用的二分脚手架
        （bracket 猜测、区间扩张兜底、报错信息），从原 tc() 里提取出来，逻辑不变。
        """

        def f(ln_t_k: float) -> float:
            return max_eig_fn(np.exp(ln_t_k)) - 1.0

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

    def _tc_constant_dos(self, t_bracket_k: Optional[Tuple[float, float]] = None,
                          xtol: float = 1e-10) -> float:
        """v0 路径：线性化核矩阵最大本征值=1 处，对 ln T 用 brentq 求根。"""

        def max_eig(t_k: float) -> float:
            t_mev = KB_MEV_PER_K * t_k
            n = self._n_for_T(t_mev)
            return self._max_eigen_linearized(t_mev, n)

        return self._bisect_tc(max_eig, t_bracket_k, xtol)

    def _tc_energy_dependent_dos(self, t_bracket_k: Optional[Tuple[float, float]] = None,
                                  xtol: float = 1e-10) -> float:
        """能量依赖 DOS 路径：在 φ→0 线性化极限下，(Z,χ) 先自洽，再用 φ 的核矩阵
        求最大本征值=1。方程与折叠约定见 README「能量依赖 DOS」。
        """

        def max_eig(t_k: float) -> float:
            t_mev = KB_MEV_PER_K * t_k
            n_pos = self._n_for_T(t_mev)
            n_arr, wn = self._full_axis_grid(t_mev, n_pos)
            lam_full = self._lambda_toeplitz_full(t_mev, n_pos)
            lam_z, lam_phi = self._apply_vertex_full(lam_full, n_arr, wn, t_mev)
            e_max = self._dos_e_max()
            z, chi, converged = self._solve_z_chi_dos(t_mev, lam_z, lam_full, wn, e_max)
            if not converged:
                warnings.warn(
                    f"能量依赖 DOS 路径 (Z,χ) 自洽在 T={t_k:.6g} K 未在 "
                    f"{self.dos_max_iter} 次迭代内收敛到 tol={self.dos_tol:.2g}"
                )
            g_int, _ = self._dos_energy_integrals(z, chi, np.zeros_like(wn), wn, e_max)
            theta = self._theta_mask_full(wn, t_mev)
            kernel = (np.pi * t_mev) * (lam_phi - self.mu_star * theta[None, :]) * g_int[None, :]
            eigvals = np.linalg.eigvals(kernel)
            i = np.argmax(eigvals.real)
            if abs(eigvals.imag[i]) > 1e-6 * max(abs(eigvals.real[i]), 1e-12):
                warnings.warn(f"能量依赖 DOS 路径主导本征值有明显虚部 ({eigvals[i]!r})")
            return float(eigvals.real[i])

        return self._bisect_tc(max_eig, t_bracket_k, xtol)

    # ---------------- 能量依赖 DOS：全轴（未折叠）辅助方法 ----------------

    @staticmethod
    @functools.lru_cache(maxsize=32)
    def _gauss_legendre_nodes(n_eps: int) -> Tuple[np.ndarray, np.ndarray]:
        """Gauss-Legendre 节点/权重（[-1,1] 上），缓存避免重复计算。"""
        return np.polynomial.legendre.leggauss(n_eps)

    def _full_axis_grid(self, t_mev: float, n_pos: int) -> Tuple[np.ndarray, np.ndarray]:
        """全体 Matsubara 指标 n=-n_pos..n_pos-1（共 M=2·n_pos 点，对称覆盖正负频率，
        不做 v0 式的正频折叠——见 README「能量依赖 DOS」关于为什么在这条路径上放弃
        折叠：χ_n 在 n→-1-n 下是奇函数（χ_{-1-n}=-χ_n），若要折叠须引入"镜像 DOS"
        N̄(-ε) 的额外积分，直接在全轴上求解更简单、风险更低，代价是矩阵边长翻倍。
        """
        n_arr = np.arange(-n_pos, n_pos)
        wn = np.pi * t_mev * (2 * n_arr + 1)
        return n_arr, wn

    def _lambda_toeplitz_full(self, t_mev: float, n_pos: int) -> np.ndarray:
        """Λ[i,j] = λ(n_i-n_j)，i,j 覆盖全轴 M=2·n_pos 个指标。由于全轴等间距枚举，
        n_i-n_j = i-j（数组下标之差等于 Matsubara 指标之差，位移抵消），
        且 λ(ν) 是偶函数，故只需 λ_table[0..M-1] 按 |i-j| 取值。
        """
        m = 2 * n_pos
        lam_k = self.spectrum.lambda_of_nu(2 * np.pi * t_mev * np.arange(m))
        idx = np.arange(m)
        diff = np.abs(idx[:, None] - idx[None, :])
        return lam_k[diff]

    def _apply_vertex_full(self, lam_full: np.ndarray, n_arr: np.ndarray, wn: np.ndarray,
                            t_mev: float) -> Tuple[np.ndarray, np.ndarray]:
        """vertex 钩子在全轴路径上的应用：返回 (Λ_Z, Λ_φ)。vertex=None 时两者都
        原样返回 lam_full（逐位一致）。χ 方程的 Λ 恒不受 vertex 影响（SPEC 只要求
        vertex 作用于 φ 方程电声部分，Z 方程由 vertex_in_Z 控制，未提及 χ 方程）。
        """
        if self.vertex is None:
            return lam_full, lam_full
        m = len(n_arr)
        n_grid = n_arr[:, None] * np.ones((1, m), dtype=int)
        m_grid = n_arr[None, :] * np.ones((m, 1), dtype=int)
        wn_n_grid = wn[:, None] * np.ones((1, m))
        wn_m_grid = wn[None, :] * np.ones((m, 1))
        v = np.asarray(self.vertex(n_grid, m_grid, wn_n_grid, wn_m_grid, t_mev), dtype=float)
        lam_phi = lam_full * (1.0 + v)
        lam_z = lam_phi if self.vertex_in_Z else lam_full
        return lam_z, lam_phi

    def _dos_e_max(self) -> float:
        """ε 积分范围的半宽：至少 dos_energy_cutoff_ratio·w_cut，table 型 DOS 再取
        表格边界（取绝对值较大者）与之的较大值，边界外由 DOS.nbar 自身延拓。
        """
        e_max = self.dos_energy_cutoff_ratio * self.w_cut
        if self.dos.kind == "table":
            eps_tab = self.dos.params["eps"]
            e_max = max(e_max, abs(float(eps_tab[0])), abs(float(eps_tab[-1])))
        return e_max

    def _theta_mask_full(self, wn: np.ndarray, t_mev: float) -> np.ndarray:
        """θ(ω_c-|ω_m|) 的平滑版本（同 _kernels 里的处理，理由见那里的注释），
        在全轴上直接用 |ω_m|（ω_m 可正可负）。
        """
        return np.clip((self.w_cut - np.abs(wn)) / (2 * np.pi * t_mev) + 0.5, 0.0, 1.0)

    def _dos_energy_integrals(self, z: np.ndarray, chi: np.ndarray, phi: np.ndarray,
                               wn: np.ndarray, e_max: float) -> Tuple[np.ndarray, np.ndarray]:
        """能量积分核心：对每个 Matsubara 指标 j，返回
            G_j = (1/π) ∫dε N̄(ε)/D_j(ε)，
            H_j = (1/π) ∫dε N̄(ε)·ε/D_j(ε)，
        D_j(ε) = (ω_j Z_j)² + (ε+χ_j)² + φ_j²。

        由于 (ω_jZ_j)²、φ_j² 都不依赖 ε，三个方程里出现的 ε 积分
        ⟨ω_jZ_j⟩=ω_jZ_j·G_j、⟨φ_j⟩=φ_j·G_j、⟨ε+χ_j⟩=χ_j·G_j+H_j 都能从这两个
        基础积分组装出来（见 README 推导），不需要对每个方程分别数值积分。

        数值方法：代换 ε+χ_j = a_j·tanθ（a_j=sqrt((ω_jZ_j)²+φ_j²)），把
        D_j(ε)=a_j²sec²θ 的 Lorentzian 型分母代换掉，积分变成 θ 上的光滑函数
        （无极点），再用固定的 Gauss-Legendre 节点数值积分，θ 范围裁剪到对应
        ε∈[-e_max, e_max] 的区间。这个代换是能量依赖 DOS 数值 Eliashberg 文献
        （Pickett 1980 附近）的标准做法：既能精确处理低阶 Matsubara 频率处很窄
        （宽度~πT）的分母峰，又能覆盖远达 e_max 的 DOS 变化，不需要对"窄峰区域"
        和"远场区域"分别处理。

        归一化核对（N̄≡1 时应回到 v0 的常 DOS 形式）：N̄≡1 时
            (1/π)∫_{-∞}^{∞} dε/(a²+ε²) = (1/π)·(π/a) = 1/a，
        即 G_j → 1/a_j = 1/sqrt((ω_jZ_j)²+φ_j²)，H_j → 0（对称区间上奇函数积分为
        零）。代入 Z_iω_i=ω_i+πTΣ_jΛ[i,j]·ω_jZ_j/a_j，在线性化极限 φ_j→0 时
        a_j→|ω_j|Z_j，恰好是 v0 折叠前（未折叠）方程的形式，见模块顶部文档字符串。
        """
        u, w = self._gauss_legendre_nodes(self.dos_n_eps)
        a = np.sqrt((wn * z) ** 2 + phi ** 2)
        theta_max = np.arctan((e_max - chi) / a)
        theta_min = np.arctan((-e_max - chi) / a)
        half = 0.5 * (theta_max - theta_min)
        mid = 0.5 * (theta_max + theta_min)
        theta = half[:, None] * u[None, :] + mid[:, None]
        eps = a[:, None] * np.tan(theta) - chi[:, None]
        nbar_vals = self.dos.nbar(eps)
        weight = w[None, :] * half[:, None]
        g_int = (1.0 / np.pi) * np.sum(weight * nbar_vals / a[:, None], axis=1)
        h_int = (1.0 / np.pi) * np.sum(weight * nbar_vals * eps / a[:, None], axis=1)
        return g_int, h_int

    def _solve_z_chi_dos(self, t_mev: float, lam_z: np.ndarray, lam_chi: np.ndarray,
                          wn: np.ndarray, e_max: float) -> Tuple[np.ndarray, np.ndarray, bool]:
        """Tc 处 (Z,χ) 的自洽解（φ→0，线性混合不动点迭代，与 v0 solve_gap 同思路）。

        Z_i = 1 + (πT/ω_i)·Σ_j Λ_Z[i,j]·ω_j Z_j G_j
        χ_i = -πT·Σ_j Λ_χ[i,j]·(χ_j G_j + H_j)
        """
        z = np.ones_like(wn)
        chi = np.zeros_like(wn)
        zero_phi = np.zeros_like(wn)
        converged = False
        for _ in range(self.dos_max_iter):
            g_int, h_int = self._dos_energy_integrals(z, chi, zero_phi, wn, e_max)
            z_new = 1.0 + (np.pi * t_mev / wn) * (lam_z @ (wn * z * g_int))
            chi_new = -np.pi * t_mev * (lam_chi @ (chi * g_int + h_int))
            denom = np.max(np.abs(z_new)) + np.max(np.abs(chi_new)) + 1e-30
            diff = (np.max(np.abs(z_new - z)) + np.max(np.abs(chi_new - chi))) / denom
            z = self.dos_mixing * z_new + (1 - self.dos_mixing) * z
            chi = self.dos_mixing * chi_new + (1 - self.dos_mixing) * chi
            if diff < self.dos_tol:
                converged = True
                break
        return z, chi, converged

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
            # dataclasses.replace 而非手写 Solver(...)：自动带上 dos/vertex/vertex_in_Z
            # 等 v0 之后新增的字段（它们是电子结构/顶角模型参数，不随同位素质量变化，
            # 应该原样传给子 Solver），也避免每次加新字段都要回来改这里。
            solver = dataclasses.replace(self, spectrum=spec, mu_star=mu_star_new, w_cut=w_cut_new)
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

        # 下面四个 make() 都用 dataclasses.replace(self, ...) 而非手写 Solver(...)：
        # 自动带上 dos/vertex/vertex_in_Z 等字段，新增字段不需要回来同步这里。
        if param == "lambda":
            def make(delta: float) -> float:
                comps = []
                for c in spec.components:
                    if c.kind == "einstein":
                        comps.append(dataclasses.replace(c, lam=c.lam * (1 + delta)))
                    else:
                        comps.append(dataclasses.replace(c, a2f=c.a2f * (1 + delta)))
                new_spec = Spectrum(components=comps)
                return dataclasses.replace(self, spectrum=new_spec).tc()

        elif param == "omega_log":
            def make(delta: float) -> float:
                s = 1.0 / (1 + delta)
                new_spec = spec.scale_mass(s)
                new_w_cut = self.w_cut * (1 + delta)
                return dataclasses.replace(self, spectrum=new_spec, w_cut=new_w_cut).tc()

        elif param == "mu_star":
            def make(delta: float) -> float:
                return dataclasses.replace(self, mu_star=self.mu_star * (1 + delta)).tc()

        elif param == "w_cut":
            def make(delta: float) -> float:
                return dataclasses.replace(self, w_cut=self.w_cut * (1 + delta)).tc()

        else:
            raise ValueError(f"不支持的 param: {param!r}")

        tc_p, tc_m = make(eps), make(-eps)
        return float((np.log(tc_p) - np.log(tc_m)) / (np.log(1 + eps) - np.log(1 - eps)))
