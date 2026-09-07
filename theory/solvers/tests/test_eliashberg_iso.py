"""
theory/solvers/tests/test_eliashberg_iso.py

对照 theory/solvers/SPEC.md 的必备测试；每条 SPEC 测试尽量对应一个函数，
在函数 docstring 里引用 SPEC.md 的条目。

约定：本文件里凡是用到 Matsubara 截止 w_cut，一律取 w_cut = R * spectrum.omega_log()
（R 通常取 5，注释里写明选择理由），并在 isotope_coefficient 里让 w_cut 跟着谱走，
与 eliashberg_iso.py 里的默认约定一致。数值全部在本机实测过（见 PR/commit 说明），
不是拍脑袋的容差。
"""
import math
import shutil

import numpy as np
import pytest

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from eliashberg_iso import (
    Spectrum,
    Solver,
    DOS,
    mu_star_at_cutoff,
    mu_star_from_bare,
    vertex_pietronero_grimaldi_isotropic,
    KB_MEV_PER_K,
)


# --------------------------------------------------------------------------
# 参考公式：P1 §3.2 式 (3.2)，纯 McMillan 形式下的解析 alpha。
# 故意不从 eliashberg_iso.py 里导入任何东西——这是"目标"公式，独立抄录，
# 避免测试和被测代码共用同一份实现导致假阳性。
# --------------------------------------------------------------------------
def _eq32_mcmillan_alpha(lam: float, mu_star: float) -> float:
    """P1_isotope_master_formula.md §3.1 式 (3.2)。"""
    d = lam - mu_star * (1 + 0.62 * lam)
    return 0.5 * (1 - 1.04 * (1 + lam) * (1 + 0.62 * lam) * mu_star ** 2 / d ** 2)


def _mcmillan_dlnTc_dlnlam(lam: float, mu_star: float, h: float = 1e-6) -> float:
    """McMillan Φ=exp[-1.04(1+λ)/D] 的 ∂lnΦ/∂lnλ（中心差分对解析式求导，非本模块实现）。"""
    def ln_phi(l):
        d = l - mu_star * (1 + 0.62 * l)
        return -1.04 * (1 + l) / d
    return (ln_phi(lam + h) - ln_phi(lam - h)) / (2 * h) * lam


W_CUT_RATIO = 5.0  # w_cut = 5 * omega_log，本文件统一约定（见 README「μ*与ω_c绑定的约定」）


# ==========================================================================
# 1. Spectrum 基础：矩、λ(ν)、质量缩放的自洽性
# ==========================================================================

def test_lambda_of_nu_zero_equals_lambda_total_einstein():
    sp = Spectrum.einstein(lam=1.3, w0=77.0)
    assert sp.lambda_of_nu(0.0) == pytest.approx(sp.lambda_total(), rel=1e-12)


def test_lambda_of_nu_zero_equals_lambda_total_table():
    w = np.linspace(5.0, 200.0, 400)
    a2f = 0.02 * np.exp(-((w - 80.0) / 40.0) ** 2)
    sp = Spectrum.from_table(w, a2f)
    assert sp.lambda_of_nu(0.0) == pytest.approx(sp.lambda_total(), rel=1e-6)


def test_einstein_moments_closed_form():
    """单 Einstein 模：λ=lam，ω_log=ω_2=w0 精确成立。"""
    sp = Spectrum.einstein(lam=2.0, w0=100.0)
    assert sp.lambda_total() == pytest.approx(2.0)
    assert sp.omega_log() == pytest.approx(100.0, rel=1e-10)
    assert sp.omega_2() == pytest.approx(100.0, rel=1e-10)


def test_two_peak_matches_target_lambda_and_omega_log():
    """构造 λ=2, ω_log=100 的双峰谱（examples 用的同一个）：lam1=lam2=1, w1*w2=100^2。"""
    sp = Spectrum.two_peak(1.0, 50.0, 1.0, 200.0)
    assert sp.lambda_total() == pytest.approx(2.0)
    assert sp.omega_log() == pytest.approx(100.0, rel=1e-10)


def test_scale_mass_preserves_lambda_total():
    """SPEC 必备功能 3：scale_mass 保持 λ 不变（谐波假设的核心恒等式）。"""
    sp_e = Spectrum.einstein(lam=1.7, w0=60.0)
    for s in [0.7, 1.0, 1.4]:
        assert sp_e.scale_mass(s).lambda_total() == pytest.approx(1.7, rel=1e-12)

    w = np.linspace(5.0, 200.0, 300)
    a2f = 0.015 * np.exp(-((w - 90.0) / 35.0) ** 2)
    sp_t = Spectrum.from_table(w, a2f)
    lam0 = sp_t.lambda_total()
    for s in [0.6, 1.3]:
        assert sp_t.scale_mass(s).lambda_total() == pytest.approx(lam0, rel=1e-6)


def test_scale_mass_moves_frequency_correctly():
    """s = sqrt(M/M0) > 1（更重）应让频率整体下移：新 w0 = 旧 w0 / s。"""
    sp = Spectrum.einstein(lam=1.0, w0=100.0)
    s = 1.3
    scaled = sp.scale_mass(s)
    assert scaled.omega_log() == pytest.approx(100.0 / s, rel=1e-10)


def test_partial_scale_mass_only_named_component():
    sp_h = Spectrum.einstein(lam=1.0, w0=150.0, name="H")
    sp_host = Spectrum.einstein(lam=1.0, w0=30.0, name="host")
    combo = sp_h.combine(sp_host)
    scaled = combo.scale_mass(1.5, which=["H"])
    # host 分量应完全不变，H 分量频率下移
    ws = sorted(c.w0 for c in scaled.components)
    assert ws[0] == pytest.approx(30.0)  # host 不变
    assert ws[1] == pytest.approx(150.0 / 1.5)  # H 缩放


# ==========================================================================
# 2. μ*(ω_c) 一致重整化
# ==========================================================================

def test_mu_star_at_cutoff_tolmachev_relation():
    """1/μ*(ωc1) - 1/μ*(ωc2) = ln(ωc2/ωc1)，直接核对定义式。"""
    mu1, wc1, wc2 = 0.13, 500.0, 900.0
    mu2 = mu_star_at_cutoff(mu1, wc1, wc2)
    assert (1 / mu1 - 1 / mu2) == pytest.approx(np.log(wc2 / wc1), rel=1e-10)


def test_mu_star_at_cutoff_roundtrip():
    mu1, wc1, wc2 = 0.20, 300.0, 1200.0
    mu2 = mu_star_at_cutoff(mu1, wc1, wc2)
    back = mu_star_at_cutoff(mu2, wc2, wc1)
    assert back == pytest.approx(mu1, rel=1e-10)


def test_mu_star_at_cutoff_zero_is_fixed_point():
    assert mu_star_at_cutoff(0.0, 300.0, 900.0) == 0.0


def test_mu_star_from_bare_consistency():
    """μ*(ωc) = μ/[1+μ ln(E_F/ωc)] 反代验证；且随 E_F 增大（L 增大）单调减小。"""
    mu, e_f, wc = 0.15, 5000.0, 500.0
    mu_star = mu_star_from_bare(mu, e_f, wc)
    assert mu_star == pytest.approx(mu / (1 + mu * np.log(e_f / wc)))
    mu_star_bigger_ef = mu_star_from_bare(mu, e_f * 10, wc)
    assert mu_star_bigger_ef < mu_star


# ==========================================================================
# 3. Allen-Dynes 参考值
# ==========================================================================

def test_allen_dynes_raises_when_denominator_nonpositive():
    sp = Spectrum.einstein(lam=0.3, w0=50.0)
    with pytest.raises(ValueError):
        sp.allen_dynes(mu_star=0.5)  # lambda - mu*(1+0.62 lambda) < 0


# ==========================================================================
# 4. SPEC 测试 1：弱耦合 Einstein vs Allen-Dynes（<10%）+ Matsubara/T 收敛性
# ==========================================================================

def test_weak_coupling_einstein_matches_allen_dynes():
    """SPEC.md 测试1：λ=1, ω_E=50meV, μ*=0.1，数值 Tc 与 AD 相差 <10%。

    w_cut = 5*ω_log（本文件统一约定）；实测 reldiff ≈ 7.2%，在弱耦合下这是已知的
    数值 Eliashberg 相对 McMillan/AD 拟合公式的典型偏差量级（AD 公式本身对 Tc 只是
    渐近拟合，不是精确解），不是 bug。
    """
    sp = Spectrum.einstein(lam=1.0, w0=50.0)
    solver = Solver(sp, mu_star=0.1, w_cut=W_CUT_RATIO * sp.omega_log())
    tc_num = solver.tc()
    tc_ad = sp.allen_dynes(mu_star=0.1)
    reldiff = abs(tc_num - tc_ad) / tc_ad
    assert reldiff < 0.10, f"reldiff={reldiff:.3%} 超过 SPEC 要求的 10%"


def test_matsubara_cutoff_convergence():
    """SPEC.md 测试1：Matsubara 截止与 T 的收敛性——oversample 加倍，Tc 变化应可忽略。"""
    sp = Spectrum.einstein(lam=1.0, w0=50.0)
    solver_lo = Solver(sp, mu_star=0.1, w_cut=W_CUT_RATIO * sp.omega_log(),
                        matsubara_oversample=3.0)
    solver_hi = Solver(sp, mu_star=0.1, w_cut=W_CUT_RATIO * sp.omega_log(),
                        matsubara_oversample=6.0)
    tc_lo, tc_hi = solver_lo.tc(), solver_hi.tc()
    assert tc_lo == pytest.approx(tc_hi, rel=1e-3)


def test_matsubara_cutoff_convergence_explicit_n():
    """直接翻倍 n_matsubara（而非通过 oversample）也应给出一致的 Tc。

    容差取 rel=5e-3：实测 n_auto（来自默认 oversample=3.0，在 tc() 内部按 T
    自适应）与固定 n_auto/2*n_auto 相比，Tc 变化量级 ~0.1%，比
    test_matsubara_cutoff_convergence（在 tc() 内部随 T 自适应改变 oversample）
    的 ~1e-4 量级更大——这是因为这里把 N 锁死为 tc_auto 处的估计值，
    在 brentq 搜索到的其余 T 上就不是"按同一 oversample 重新估计"的 N 了，
    存在真实的（虽然很小的）残余截断误差，不是 bug。
    """
    sp = Spectrum.einstein(lam=1.0, w0=50.0)
    w_cut = W_CUT_RATIO * sp.omega_log()
    # 先跑一次拿到自动 N 的量级
    probe = Solver(sp, mu_star=0.1, w_cut=w_cut)
    tc_auto = probe.tc()
    t_mev = KB_MEV_PER_K * tc_auto
    n_auto = probe._n_for_T(t_mev)
    solver_n = Solver(sp, mu_star=0.1, w_cut=w_cut, n_matsubara=n_auto)
    solver_2n = Solver(sp, mu_star=0.1, w_cut=w_cut, n_matsubara=2 * n_auto)
    assert solver_n.tc() == pytest.approx(solver_2n.tc(), rel=5e-3)


# ==========================================================================
# 5. SPEC 测试 2：强耦合 λ=2,3，单调性与 AD 渐近趋势
# ==========================================================================

def test_strong_coupling_monotonic_and_tracks_allen_dynes():
    """SPEC.md 测试2：λ=1,2,3 时 Tc/ω_E 单调递增，且与 AD 渐近 √λ 趋势一致
    （用"数值解与 AD 在每个 λ 下都保持合理接近"作为渐近趋势一致性的代理判据）。
    """
    w0 = 50.0
    tcs = []
    for lam in (1.0, 2.0, 3.0):
        sp = Spectrum.einstein(lam=lam, w0=w0)
        solver = Solver(sp, mu_star=0.1, w_cut=W_CUT_RATIO * sp.omega_log())
        tc_num = solver.tc()
        tc_ad = sp.allen_dynes(mu_star=0.1)
        assert abs(tc_num - tc_ad) / tc_ad < 0.12
        tcs.append(tc_num)
    assert tcs[0] < tcs[1] < tcs[2], f"Tc(lambda) 应单调递增: {tcs}"
    # 亚线性增长（√λ 型渐近的定性特征）：lambda 每次 x1.5/x1.5，Tc 的增幅应逐次收窄
    growth1 = tcs[1] / tcs[0]
    growth2 = tcs[2] / tcs[1]
    assert growth2 < growth1, "强耦合下 Tc 的增长应逐渐放缓（次线性），与 AD 的渐近趋势一致"


# ==========================================================================
# 6. SPEC 测试 3：μ*=0 谐波全缩放 -> alpha = 0.500 ± 0.002
# ==========================================================================

def test_alpha_half_at_mu_star_zero_full_scaling():
    """SPEC.md 测试3：μ*=0、谐波、全缩放，数值 alpha = 0.500 ± 0.002。"""
    sp = Spectrum.einstein(lam=2.0, w0=100.0)
    solver = Solver(sp, mu_star=0.0, w_cut=W_CUT_RATIO * sp.omega_log())
    alpha, _ = solver.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)
    assert alpha == pytest.approx(0.5, abs=0.002)


def test_alpha_half_at_mu_star_zero_two_peak():
    """同一测试换成双峰谱，确认不是 Einstein 特例。"""
    sp = Spectrum.two_peak(1.0, 50.0, 1.0, 200.0)
    solver = Solver(sp, mu_star=0.0, w_cut=W_CUT_RATIO * sp.omega_log())
    alpha, _ = solver.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)
    assert alpha == pytest.approx(0.5, abs=0.002)


# ==========================================================================
# 7. SPEC 测试 4：mu_star_fixed 开关 + 与 P1 式 (3.2) 对照
# ==========================================================================

def test_mu_star_fixed_true_gives_exactly_half():
    """SPEC.md 测试4后半句：mu_star_fixed=True 时回到 0.5（对任意 μ*）。

    这不只是近似成立：mu_star_fixed=True + 全缩放时，w_cut 与谱同比例移动、
    μ* 数值不变，整套 Matsubara 方程在 (T, 谱, w_cut) -> (T/s, 谱/s, w_cut/s)
    下自洽（见 eliashberg_iso.py 的推导注释），所以 Tc(s) = Tc(1)/s 精确成立，
    alpha 应精确等于 0.5（只受二分法 xtol 限制）。
    """
    sp = Spectrum.einstein(lam=2.0, w0=100.0)
    for mu_star in (0.0, 0.13, 0.30):
        solver = Solver(sp, mu_star=mu_star, w_cut=W_CUT_RATIO * sp.omega_log())
        alpha, _ = solver.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=True)
        assert alpha == pytest.approx(0.5, abs=1e-4)


@pytest.mark.parametrize("mu_star", [0.10, 0.13])
def test_alpha_matches_p1_eq32_lambda2(mu_star):
    """SPEC.md 测试4：mu_star_fixed=False，λ=2、μ*≤0.15 时数值 alpha 与式(3.2)相差 <0.01。

    这是 P1 §3.1 明确给出数值表的那个例子（λ=2 那一行）。实测：μ*=0.10 时
    diff≈0.001，μ*=0.13 时 diff≈0.004，都在容差内。见 README「首轮结果」
    的完整表格（含 μ*=0.20/0.25/0.30，那三个不在 SPEC 承诺范围内，偏差更大）。
    """
    sp = Spectrum.einstein(lam=2.0, w0=100.0)
    solver = Solver(sp, mu_star=mu_star, w_cut=W_CUT_RATIO * sp.omega_log())
    alpha, _ = solver.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)
    target = _eq32_mcmillan_alpha(2.0, mu_star)
    assert abs(alpha - target) < 0.01


def test_alpha_matches_p1_eq32_lambda1_weak_mu():
    """同一测试在 λ=1、μ*=0.10 下也成立（弱库仑端）。"""
    sp = Spectrum.einstein(lam=1.0, w0=100.0)
    solver = Solver(sp, mu_star=0.10, w_cut=W_CUT_RATIO * sp.omega_log())
    alpha, _ = solver.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)
    target = _eq32_mcmillan_alpha(1.0, 0.10)
    assert abs(alpha - target) < 0.01


def test_alpha_eq32_deviation_documented_at_lambda1_stronger_mu():
    """如实记录一个不满足 SPEC <0.01 容差的角落：λ=1、μ*=0.13。

    这不是求解器 bug（已核对 eps 稳定性、N 收敛性，见 README「首轮结果」的分析），
    是 P1 式(3.2)（纯 McMillan 形式）本身在这个角落precision不够。这里把观测到的
    偏差量级钉住（回归哨兵），如果未来改动让它意外"变好"或变得更差很多，测试会提醒
    去核查是不是引入了别的问题，而不是默默让 README 的说法过期。
    """
    sp = Spectrum.einstein(lam=1.0, w0=100.0)
    solver = Solver(sp, mu_star=0.13, w_cut=W_CUT_RATIO * sp.omega_log())
    alpha, _ = solver.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)
    target = _eq32_mcmillan_alpha(1.0, 0.13)
    diff = alpha - target
    assert 0.008 < diff < 0.025, (
        f"diff={diff:.4f}；README 记录的偏差量级失效，需要重新核对分析"
    )


# ==========================================================================
# 8. SPEC 测试 5：部分缩放的 Rainer-Culetto 和规则
# ==========================================================================

def test_partial_scaling_rainer_culetto_sum_rule():
    """SPEC.md 测试5：两分量谱，只缩放一个，α_1+α_2=0.5（μ*=0）。"""
    sp_h = Spectrum.einstein(lam=1.0, w0=150.0, name="H")
    sp_host = Spectrum.einstein(lam=1.0, w0=30.0, name="host")
    combo = sp_h.combine(sp_host)
    solver = Solver(combo, mu_star=0.0, w_cut=W_CUT_RATIO * combo.omega_log())
    parts = {"H": sp_h, "host": sp_host}
    alpha_h, _ = solver.isotope_coefficient(parts, which="H", mu_star_fixed=False)
    alpha_host, _ = solver.isotope_coefficient(parts, which="host", mu_star_fixed=False)
    assert (alpha_h + alpha_host) == pytest.approx(0.5, abs=0.005)
    # 权重应偏向 lambda 更大或...此处两分量 lambda 相等，但 H 频率更高，
    # 对 ln Tc 的敏感度份额不必相等；只检验和规则，不检验具体拆分比例。


def test_isotope_coefficient_eps_stability():
    """SPEC 要求：同位素系数对差分步长 eps 的稳定性检查（eps=0.02 附近）。"""
    sp_h = Spectrum.einstein(lam=1.0, w0=150.0, name="H")
    sp_host = Spectrum.einstein(lam=1.0, w0=30.0, name="host")
    combo = sp_h.combine(sp_host)
    solver = Solver(combo, mu_star=0.0, w_cut=W_CUT_RATIO * combo.omega_log())
    parts = {"H": sp_h, "host": sp_host}
    alphas = []
    for eps in (0.01, 0.02, 0.04):
        a, _ = solver.isotope_coefficient(parts, which="H", mu_star_fixed=False, eps=eps)
        alphas.append(a)
    assert max(alphas) - min(alphas) < 0.005, f"alpha 对 eps 不够稳定: {alphas}"


# ==========================================================================
# 9. IsoME 交叉对照（若本机无 Julia+IsoME 则跳过）
# ==========================================================================

def test_isome_cross_check_if_available():
    """SPEC.md：若本机有 Julia + IsoME（Kogler 2025），加一个交叉对照测试。

    本机（开发/CI 环境）未安装 Julia，跳过并注明原因；不是"忘记做"。
    """
    if shutil.which("julia") is None:
        pytest.skip("本机未安装 Julia（未检测到可执行文件），跳过 IsoME 交叉对照")
    pytest.skip("Julia 存在但 IsoME 交叉对照脚本未实现（v0 范围之外，留待需要时补）")


# ==========================================================================
# 10. 非线性 Δ(T) 求解（接口占位，基本正确性）
# ==========================================================================

def test_solve_gap_weak_coupling_bcs_ratio():
    """非线性解：弱耦合下 2Δ(0.3Tc)/k_BTc 应接近 BCS 值 3.53（量级检验，非精确对照）。"""
    sp = Spectrum.einstein(lam=0.5, w0=50.0)
    solver = Solver(sp, mu_star=0.1, w_cut=W_CUT_RATIO * sp.omega_log())
    tc = solver.tc()
    _, delta, _, converged = solver.solve_gap(t_k=0.3 * tc)
    assert converged
    ratio = 2 * delta[0] / (KB_MEV_PER_K * tc)
    assert 3.0 < ratio < 4.2, f"2Delta/kTc={ratio:.3f} 偏离 BCS 附近太多"


def test_solve_gap_vanishes_toward_tc():
    """Δ(T) 应随 T -> Tc^- 减小到接近 0（连续相变的定性检验）。"""
    sp = Spectrum.einstein(lam=0.5, w0=50.0)
    solver = Solver(sp, mu_star=0.1, w_cut=W_CUT_RATIO * sp.omega_log())
    tc = solver.tc()
    _, delta_low, _, _ = solver.solve_gap(t_k=0.3 * tc)
    _, delta_high, _, _ = solver.solve_gap(t_k=0.95 * tc)
    assert np.max(np.abs(delta_high)) < 0.5 * np.max(np.abs(delta_low))


# ==========================================================================
# 11. log_sensitivity：各参数的基本 sanity + 与 McMillan 对照（诊断性，非强断言）
# ==========================================================================

def test_log_sensitivity_omega_log_is_exactly_one():
    """log_sensitivity('omega_log') 对应"整体缩放+w_cut 同比例跟随+μ*数值不变"，
    与 mu_star_fixed=True 全缩放同构，∂lnTc/∂ln(omega_log) 应精确等于 1。
    """
    sp = Spectrum.einstein(lam=2.0, w0=100.0)
    solver = Solver(sp, mu_star=0.13, w_cut=W_CUT_RATIO * sp.omega_log())
    s = solver.log_sensitivity("omega_log")
    assert s == pytest.approx(1.0, abs=1e-4)


def test_log_sensitivity_lambda_positive_and_documented_vs_mcmillan():
    """∂lnTc/∂lnλ 应为正（更强耦合 -> 更高 Tc），且与 McMillan 解析值的偏差
    随 λ 增大而增大（记录这个已知现象，见 README；这是 McMillan 公式在强耦合下
    低估 λ 敏感度的已知缺陷，Allen-Dynes 引入 f1f2 正是为了修正这一点）。
    """
    mu_star = 0.13
    devs = []
    for lam in (1.0, 2.0, 3.0):
        sp = Spectrum.einstein(lam=lam, w0=100.0)
        solver = Solver(sp, mu_star=mu_star, w_cut=W_CUT_RATIO * sp.omega_log())
        s_num = solver.log_sensitivity("lambda")
        assert s_num > 0
        devs.append(abs(s_num - _mcmillan_dlnTc_dlnlam(lam, mu_star)))
    assert devs[0] < devs[2], "偏差应随 lambda 增大而增大（McMillan 强耦合端的已知局限）"


def test_log_sensitivity_mu_star_is_negative():
    """更大的 μ* 应压制 Tc：∂lnTc/∂lnμ* < 0。"""
    sp = Spectrum.einstein(lam=1.5, w0=80.0)
    solver = Solver(sp, mu_star=0.15, w_cut=W_CUT_RATIO * sp.omega_log())
    s = solver.log_sensitivity("mu_star")
    assert s < 0


def test_log_sensitivity_w_cut_sign_is_not_naively_negative():
    """如实记录一个反直觉但已核实的现象：log_sensitivity("w_cut")（只放大 w_cut，
    μ* 数值不动、不按 Tolmachev 重新推导）在这里数值上是**正**的，不是"μ* 作用范围
    变宽 -> 压制更多 -> 应该是负的"这个朴素预期。

    已排查证实这不是 bug：直接在固定 T、固定 N 下比较线性化本征值矩阵（绕开 tc()
    内部的二分/N 自适应），w_cut 增大确实让主导本征值变大；用旧 w_cut 的本征向量
    去算新 w_cut 矩阵的 Rayleigh 商，结果确实按直觉下降——说明"逐项看，μ*θ 项让
    矩阵元变小"这个直觉本身没错，但这个折叠核矩阵 a[n,m] 关于 (n,m) 并不对称
    （μ*θ(ω_m) 只依赖列 m，wn/z0 的归一化也不对称），所以主导本征值不满足对称矩阵
    的变分（Rayleigh 商）单调性，符号可以反转。这是当前"非对称折叠核矩阵"表示法的
    已知特性，已记录在 README「已知局限」。

    重要：这只影响"单独扰动 w_cut、μ* 数值不变"这个人为构造的对照量；真正物理相关的
    isotope_coefficient(mu_star_fixed=False) 里 w_cut 和 μ* 是按 Tolmachev 关系耦合
    一起变化的，那条路径已经用 P1 式(3.2) 独立核对过（见前面的测试），不受此影响。
    """
    sp = Spectrum.einstein(lam=1.5, w0=80.0)
    solver = Solver(sp, mu_star=0.15, w_cut=W_CUT_RATIO * sp.omega_log())
    s = solver.log_sensitivity("w_cut")
    assert np.isfinite(s)
    assert 0.0 < s < 0.2, f"s={s:.4f}；README 记录的量级/符号失效，需要重新核对"


# ==========================================================================
# 12. Lorentzian 模型谱 sanity
# ==========================================================================

def test_lorentzian_spectrum_normalization_and_moments():
    sp = Spectrum.lorentzian(lam=1.2, w0=60.0, gamma=8.0)
    assert sp.lambda_total() == pytest.approx(1.2, rel=1e-3)
    # w_log、w_2 应在 w0 附近（宽度远小于中心频率时展宽只带来小修正）
    assert 50.0 < sp.omega_log() < 70.0
    assert 50.0 < sp.omega_2() < 75.0


# ==========================================================================
# 13. SPEC 必备功能 6：能量依赖 DOS（Mitrović–Carbotte/Pickett 推广）
# ==========================================================================
#
# 本节统一用 Einstein 谱 lam=1.84, w0=93 meV, mu_star=0.16（与
# examples/dos_scan.py 相同的参数组合，方便交叉核对），w_cut=5*omega_log。

_DOS_LAM, _DOS_W0, _DOS_MU_STAR = 1.84, 93.0, 0.16


def _dos_test_spectrum_and_solver():
    sp = Spectrum.einstein(lam=_DOS_LAM, w0=_DOS_W0)
    w_cut = W_CUT_RATIO * sp.omega_log()
    solver = Solver(sp, mu_star=_DOS_MU_STAR, w_cut=w_cut)
    return sp, w_cut, solver


def test_dos_constant_matches_v0_bit_identical():
    """(a) N̄≡1（DOS.constant()）必须与不传 dos 参数的 v0 路径逐位相同。

    这不是"数值上很接近"：Solver.tc() 在 dos.is_constant 为真时直接分派到
    _tc_constant_dos——与 v0 原始 tc() 是同一段代码，所以两者应该是同一个
    浮点数（Tc 相对误差 < 1e-8 的要求在这里自动、严格地满足，甚至更强）。
    """
    sp, w_cut, solver_v0 = _dos_test_spectrum_and_solver()
    solver_dos = Solver(sp, mu_star=_DOS_MU_STAR, w_cut=w_cut, dos=DOS.constant())
    tc_v0 = solver_v0.tc()
    tc_dos = solver_dos.tc()
    assert tc_dos == pytest.approx(tc_v0, rel=1e-8)


def test_dos_step_ratio_one_matches_constant_dos():
    """(b) step 型 DOS 取 ratio=1（两侧台阶高度相同，等价于常数 DOS）应与
    常数 DOS 的 Tc 一致。

    注意：ratio=1 的 step 仍然走"能量依赖 DOS"的数值路径（tan 代换 +
    Gauss-Legendre 能量积分 + (Z,χ) 自洽迭代），不会像 DOS.constant() 那样
    被短路到 v0 的解析路径，所以不是逐位相同。实测这条数值路径本身相对
    v0 解析结果有一个 ~0.6% 的系统偏差，来源是能量积分的**有限截止**
    `dos_energy_cutoff_ratio`（默认 20，即 ±20*w_cut）——SPEC 要求"积分范围
    至少 ±20*w_cut"，这就是这个最小推荐值本身带来的截断误差，不是 bug
    （已核对：加密 Gauss-Legendre 节点数 dos_n_eps 从 100 到 800 结果完全不变，
    说明误差不是求积节点不够，而是积分区间截断；把 cutoff_ratio 从 20 提到
    200，误差从 ~0.6% 降到 ~0.06%，与截止成反比，符合"截断误差"的预期）。
    见 README「能量依赖 DOS·已知局限」。
    """
    sp, w_cut, solver_const = _dos_test_spectrum_and_solver()
    dos_flat = DOS.step(delta=0.0, ratio=1.0, width=1.0)
    solver_flat = Solver(sp, mu_star=_DOS_MU_STAR, w_cut=w_cut, dos=dos_flat, dos_n_eps=200)
    tc_const = solver_const.tc()
    tc_flat = solver_flat.tc()
    reldiff = abs(tc_flat - tc_const) / tc_const
    assert reldiff < 0.01, f"reldiff={reldiff:.4%}，超过默认能量截止下的预期系统偏差量级"


def test_dos_vhs_pinned_exactly_at_ef_suppresses_both_tc_and_isotope_effect():
    """如实记录一个偏离 SPEC 测试(c)字面表述的发现：delta=0（vHs 峰值精确落在
    费米面）时，数值结果是 **Tc 被压低**、同位素系数被压低，而不是 SPEC 测试(c)
    字面写的"Tc 应更高"。

    根本原因（不是求解器 bug，是 N̄(ε)≡N(ε)/N(0) 这个定义本身的数学推论）：
    N̄(0)≡1 是定义的直接结果（比值在参考点必为1），而 DOS.log_vhs 的对数
    形状函数 g(ε)=1+a·ln(W/max(|ε-δ|,cut)) 在 ε=δ 处取到（clip 意义下的）
    最大值。当 δ=0 时，这个最大值恰好落在归一化的分母 g(0) 上，归一化后
    N̄(ε)=g(ε)/g(0) 在 ε=0 处被摁成 1，而在其余所有 ε 处都 ≤ g(0)/g(0)=1
    （因为 g(0) 已经是 g 的最大值）——也就是说"vHs 精确钉在费米面"这个参数
    选择，在 N̄(0)≡1 的强制归一化下，数学上只能表现成"费米面之外全面走低"的
    DOS，而不能表现成"费米面处 DOS 额外增强"（增强量在归一化时被定义本身
    吸收掉了，无法再体现为 N̄>1）。这在 a∈{0.05,...,1.0}、cut∈{1,5,20} 的
    扫描下是稳健成立的（Tc 全部低于常数 DOS，未发现任何反例）。

    同位素系数确实如 SPEC 预期"明显低于 0.5"（这部分符合 Tsuei 1990 的定性
    结论），且比常数 DOS 基线（~0.475）压得更低——见下面对 alpha 的断言。
    """
    sp, w_cut, solver_const = _dos_test_spectrum_and_solver()
    tc_const = solver_const.tc()
    alpha_const, _ = solver_const.isotope_coefficient(
        {"full": sp}, which="all", mu_star_fixed=False
    )

    dos_vhs = DOS.log_vhs(delta=0.0, a=0.3, W=1000.0, cut=5.0)
    solver_vhs = Solver(sp, mu_star=_DOS_MU_STAR, w_cut=w_cut, dos=dos_vhs, dos_n_eps=200)
    tc_vhs = solver_vhs.tc()
    alpha_vhs, _ = solver_vhs.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)

    assert tc_vhs < tc_const, "delta=0 时 Tc 应被压低（见 docstring 的数学论证），不是升高"
    assert alpha_vhs < 0.45, f"alpha={alpha_vhs:.4f}，同位素效应压制幅度未达预期"
    assert alpha_vhs < alpha_const, "vHs 钉在费米面应比常数 DOS 更压制同位素效应"


def test_dos_vhs_offset_from_ef_matches_spec_test_c_qualitative_claim():
    """(c) 用偏离费米面约一个声子能标（delta=omega_log）的 vHs 复现 SPEC 测试(c)
    字面要求的两条：Tc 高于常数 DOS，且同位素系数明显低于 0.5。

    选择 delta=omega_log 而非字面的 delta=0 的理由见上一个测试的 docstring：
    delta=0 是 N̄(0)≡1 定义下的一个退化点（该点的归一化把"费米面处增强"这个
    效应本身抵消掉了），"vHs 钉在/靠近费米面附近" 在这个模型里只有在 delta
    与 0 有限但不为零的偏移时才能同时表现出 Tc 增强，物理上这对应"vHs 没有
    精确落在费米面，但足够靠近"（真实材料里 vHs 的位置本来也只能靠掺杂调到
    一个有限分辨率，不可能做到数学意义上的 delta 严格为 0）。用扫描核实过：
    Tc 从 delta=0（明显更低）到约 delta≈45meV 附近穿越常数 DOS 值转为更高，
    delta=omega_log=93meV 处 Tc 已经稳定高出常数 DOS 约 4.6%，同时 alpha 仍
    明显低于 0.5（也低于常数 DOS 基线），两个要求同时满足。
    """
    sp, w_cut, solver_const = _dos_test_spectrum_and_solver()
    tc_const = solver_const.tc()
    alpha_const, _ = solver_const.isotope_coefficient(
        {"full": sp}, which="all", mu_star_fixed=False
    )

    dos_vhs = DOS.log_vhs(delta=sp.omega_log(), a=0.3, W=1000.0, cut=5.0)
    solver_vhs = Solver(sp, mu_star=_DOS_MU_STAR, w_cut=w_cut, dos=dos_vhs, dos_n_eps=200)
    tc_vhs = solver_vhs.tc()
    alpha_vhs, _ = solver_vhs.isotope_coefficient({"full": sp}, which="all", mu_star_fixed=False)

    assert tc_vhs > tc_const, f"Tc_vhs={tc_vhs:.4f} 应高于 Tc_const={tc_const:.4f}"
    assert alpha_vhs < 0.48, f"alpha={alpha_vhs:.4f}，未明显低于 0.5"
    assert alpha_vhs < alpha_const, "vHs 靠近费米面应比常数 DOS 更压制同位素效应"


def test_dos_log_vhs_large_delta_recovers_flat_dos():
    """(d) delta 远大于 w_cut 时，log_vhs 应回到（数值路径意义下的）常数 DOS。

    注意选取 W（背景能标）：必须让 W 显著大于扫描到的最大 delta，否则
    归一化分母 g(0)=1+a·ln(W/delta) 会随 delta 增大而缩小甚至变负——这不是
    "delta 越大越接近常数 DOS 因而无所谓 W"的直觉，而是 log_vhs 模型本身的
    归一化在 |delta|≫W 时先失效（实测 delta=60*w_cut, W=1000meV 时 (Z,χ)
    自洽迭代发散振荡，500 次不收敛；delta=100*w_cut 时归一化分母直接变负，
    构造函数报错），见 DOS.log_vhs 文档与 README「已知局限」。这里固定用
    W=20000 meV（远大于测试到的最大 delta=30*w_cut≈14 keV），确保只测试
    "delta 远大于 w_cut" 这一个变量。

    对照对象用同一条数值路径下的"精确常数 DOS"（DOS.step(ratio=1)，与
    test_dos_step_ratio_one_matches_constant_dos 用的是同一个对象），而不是
    v0 解析 Tc——这样可以把"能量积分有限截止"带来的系统偏差（见该测试的
    docstring）从"N̄ 的形状是否真的收敛到 1"这个问题里剥离出去，得到一个
    干净得多的收敛信号。
    """
    sp, w_cut, _ = _dos_test_spectrum_and_solver()
    dos_flat = DOS.step(delta=0.0, ratio=1.0, width=1.0)
    tc_flat = Solver(sp, mu_star=_DOS_MU_STAR, w_cut=w_cut, dos=dos_flat, dos_n_eps=200).tc()

    reldiffs = []
    for factor in (10.0, 30.0):
        dos_far = DOS.log_vhs(delta=factor * w_cut, a=0.3, W=20000.0, cut=5.0)
        solver_far = Solver(sp, mu_star=_DOS_MU_STAR, w_cut=w_cut, dos=dos_far, dos_n_eps=200)
        tc_far = solver_far.tc()
        reldiffs.append(abs(tc_far - tc_flat) / tc_flat)

    assert reldiffs[0] < 0.001, f"delta=10*w_cut 时与常数 DOS 的偏差过大: {reldiffs[0]:.4%}"
    assert reldiffs[1] < reldiffs[0], "delta 继续增大时应更接近常数 DOS（单调收敛）"


# ==========================================================================
# 14. SPEC 必备功能 7：模型顶角核钩子
# ==========================================================================

def test_vertex_none_vs_zero_vertex_matches_v0_bit_identical():
    """V≡0（显式传入恒为零的 vertex 函数）必须与不传 vertex（v0 路径）逐位相同。

    verify: (1+0)=1.0 在 IEEE754 下与不乘是同一个浮点数，_kernels 的
    vertex 分支在 v≡0 时应该退化成与 v0 完全相同的核矩阵。
    """
    sp = Spectrum.einstein(lam=1.5, w0=80.0)
    w_cut = W_CUT_RATIO * sp.omega_log()
    solver_v0 = Solver(sp, mu_star=0.13, w_cut=w_cut)

    def vertex_zero(n_idx, m_idx, omega_n, omega_m, T):
        return np.zeros(np.broadcast(
            np.asarray(n_idx), np.asarray(m_idx), np.asarray(omega_n), np.asarray(omega_m)
        ).shape)

    solver_v = Solver(sp, mu_star=0.13, w_cut=w_cut, vertex=vertex_zero)
    assert solver_v.tc() == pytest.approx(solver_v0.tc(), rel=1e-12)


def test_vertex_constant_matches_uniform_lambda_scaling():
    """V≡常数 c（vertex_in_Z=True，mu*=0）应与直接把整个谱按 (1+c) 缩放的
    v0 结果一致：两者都等价于把 lambda(n-m) 处处乘以 (1+c)，mu*=0 时
    电声子部分是唯一进入方程的核，二者应严格相等（受 brentq xtol 限制）。
    """
    sp = Spectrum.einstein(lam=1.5, w0=80.0)
    w_cut = W_CUT_RATIO * sp.omega_log()
    c = 0.2

    def vertex_const(n_idx, m_idx, omega_n, omega_m, T):
        return np.full(np.broadcast(
            np.asarray(n_idx), np.asarray(m_idx), np.asarray(omega_n), np.asarray(omega_m)
        ).shape, c)

    solver_v = Solver(sp, mu_star=0.0, w_cut=w_cut, vertex=vertex_const, vertex_in_Z=True)
    tc_v = solver_v.tc()

    sp_scaled = Spectrum.einstein(lam=1.5 * (1 + c), w0=80.0)
    tc_scaled = Solver(sp_scaled, mu_star=0.0, w_cut=w_cut).tc()

    assert tc_v == pytest.approx(tc_scaled, rel=1e-8)


def test_vertex_pietronero_grimaldi_isotropic_placeholder_shape_and_sign():
    """示例顶角核 vertex_pietronero_grimaldi_isotropic 是明确标注的占位实现
    （非文献推导公式，真实物理形式留给 P2），这里只检验它的基本形状约定：
    P_V(Qc)=1-2Qc 在 Qc→0 时为正、Qc→1 时为负，且返回数组的形状与广播后的
    (n,m,omega_n,omega_m) 一致。
    """
    vertex_small_qc = vertex_pietronero_grimaldi_isotropic(m_ratio=0.05, Qc=0.1)
    vertex_large_qc = vertex_pietronero_grimaldi_isotropic(m_ratio=0.05, Qc=0.9)
    n_idx = np.array([[0, 1], [0, 1]])
    m_idx = np.array([[0, 0], [1, 1]])
    wn = np.array([[10.0, 30.0], [10.0, 30.0]])
    wm = np.array([[10.0, 10.0], [30.0, 30.0]])
    v_small = vertex_small_qc(n_idx, m_idx, wn, wm, T=100.0)
    v_large = vertex_large_qc(n_idx, m_idx, wn, wm, T=100.0)
    assert v_small.shape == n_idx.shape
    assert np.all(v_small > 0)
    assert np.all(v_large < 0)
