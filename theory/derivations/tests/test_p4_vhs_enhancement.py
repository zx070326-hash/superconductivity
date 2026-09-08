"""
theory/derivations/tests/test_p4_vhs_enhancement.py

对 p4_vhs_enhancement.py（P4 §5 数值验证计划）的单元测试。风格与
tests/test_p2_vertex*.py 一致：中文 docstring 引用 P4_vhs_dtc_dp_sign.md 的
相应条目；解析公式在测试文件里**独立抄写**，不从被测模块 import，避免测试
与实现共用同一份公式导致假阳性。

覆盖范围：
  A. 模型 DOS 与归一化（P4 (1) 式 + 模块决策 1 的 n0 因子）
  B. 窗口平均闭式 (3) 与数值求积 (2) 的一致性、以及 §4 的四个极限
  C. 拟合机器本身的正确性（喂精确洛伦兹进去必须还原 W、R2=1）
  D. 求解器侧：a=0 时 E=1；dlt->inf 时 E->1；E(dlt) 严格偶；
     sign(dE/ddlt) = -sign(dlt)（P4 (6) 式）
  E. **窗口平均恒等式 (2) 在求解器里逐项精确成立**（含 chi 位移），
     以及不含 chi 时的偏差明显更大
  F. 单洛伦兹形式 (4) 的**失效**幅度（W_amp / W_curv / W_half 三种读法不一致）

关于任务书里"洛伦兹拟合 R2>0.99 的断言"：**这一条按实测不成立**，见
test_single_lorentzian_form_fails_quantitatively 的 docstring。测试如实断言
实测到的失效幅度，而不是把容差放宽到能通过——把一个不成立的预言测成"通过"
才是真正的 bug。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from p4_vhs_enhancement import (  # noqa: E402
    PeakDOS,
    Setup,
    TcMap,
    e_analytic,
    fit_lorentz_3p,
    fit_weff_1p,
    flat_window_weights,
    ntilde_lorentz,
    ntilde_multi,
    ntilde_numeric,
    scan_delta,
    tc_flat,
    tc_peak,
    width_diagnostics,
    window_identity_check,
)

# 测试用的轻量基线（比脚本默认更小的 dos_n_eps / 映射网格，为了跑得快；
# dos_n_eps 从 400 降到 200 对本文件的有理型 Nbar 影响 <1e-10，见脚本 docstring 决策 4）
BASE = Setup(lam_b=1.5, mu_b=0.13, w_log=110.0, a=2.0, gamma=100.0, n_eps=200)


@pytest.fixture(scope="module")
def maps():
    """平坦 DOS 参考映射（module 级缓存：每次构造要解 ~2x12 次 Tc）。"""
    return {"kernel": TcMap(BASE, "kernel", n_grid=12, s_lo=0.85, s_hi=2.6),
            "lambda": TcMap(BASE, "lambda", n_grid=12, s_lo=0.85, s_hi=2.6)}


# ==========================================================================
# A. 模型 DOS 与 n0 归一化
# ==========================================================================

def test_peak_dos_matches_eq1_and_normalizes_at_ef():
    """PeakDOS 的 N/N_b 就是 P4 (1) 式；Nbar = (N/N_b)/n0 且 Nbar(0)=1。"""
    a, gam, dlt = 2.0, 100.0, 137.0
    dos = PeakDOS(a, gam, dlt)
    xi = np.array([-800.0, -dlt, -50.0, 0.0, 50.0, 800.0])
    ref = 1.0 + a * gam ** 2 / ((xi + dlt) ** 2 + gam ** 2)   # (1) 式，独立抄写
    assert dos.ratio_over_nb(xi) == pytest.approx(ref, rel=1e-14)
    n0_ref = 1.0 + a * gam ** 2 / (dlt ** 2 + gam ** 2)
    assert dos.n0 == pytest.approx(n0_ref, rel=1e-14)
    assert float(dos.nbar(np.array(0.0))) == pytest.approx(1.0, abs=1e-14)
    # 峰顶（xi = -dlt）处 N/N_b = 1+a
    assert float(dos.ratio_over_nb(np.array(-dlt))) == pytest.approx(1.0 + a, rel=1e-14)


def test_peak_dos_a_zero_is_exactly_flat():
    """a=0（平坦参考态）时 Nbar 恒等于 1、n0=1——参考态与峰态走同一条数值路径。"""
    dos = PeakDOS(0.0, 100.0, 250.0)
    assert dos.n0 == 1.0
    assert np.allclose(dos.nbar(np.linspace(-5000, 5000, 101)), 1.0, atol=0.0)


# ==========================================================================
# B. 窗口平均 (2)/(3) 与 P4 §4 的极限
# ==========================================================================

@pytest.mark.parametrize("w", [1.0, 30.0, 110.0, 900.0])
@pytest.mark.parametrize("dlt", [0.0, 77.0, 400.0])
def test_ntilde_closed_form_matches_quadrature(w, dlt):
    """(3) 式闭式 = (2) 式数值求积（两个洛伦兹的卷积）。"""
    a, gam = 1.7, 60.0
    closed = 1.0 + a * gam * (gam + w) / (dlt ** 2 + (gam + w) ** 2)   # (3)，独立抄写
    num = float(ntilde_numeric(w, "lorentz", a, gam, dlt))
    assert num == pytest.approx(closed, rel=1e-10)
    assert ntilde_lorentz(w, a, gam, dlt) == pytest.approx(closed, rel=1e-14)


def test_ntilde_limits_of_p4_section4():
    """P4 §4 的极限：a->0 给 1；W->0 给 N(eps_F)/N_b；dlt->inf 给 1；W 大给 1+a*Gam/W。"""
    a, gam = 2.0, 100.0
    assert ntilde_lorentz(300.0, 0.0, gam, 50.0) == pytest.approx(1.0)
    # W->0：回到 BCS-with-DOS，就是费米面处的 DOS 值
    for dlt in (0.0, 150.0):
        n_ef = 1.0 + a * gam ** 2 / (dlt ** 2 + gam ** 2)
        assert ntilde_lorentz(1e-9, a, gam, dlt) == pytest.approx(n_ef, rel=1e-8)
    # dlt->inf
    assert ntilde_lorentz(200.0, a, gam, 1e7) == pytest.approx(1.0, abs=1e-8)
    # W 很大：增强 ~ 峰面积/窗宽
    big = 1e6
    assert ntilde_lorentz(big, a, gam, 0.0) - 1.0 == pytest.approx(a * gam / big, rel=1e-4)


def test_ntilde_wide_and_narrow_peak_limits():
    """P4 §4 两行：Gam>>W 时 W_eff->Gam（推迟无关）；Gam<<W 时增强 ~ a*Gam/W。"""
    a = 2.0
    w = 110.0
    # 宽峰 Gam >> W：Ntilde -> 1 + a*Gam^2/(dlt^2+Gam^2)（就是费米面 DOS 值）
    gam = 100.0 * w
    for dlt in (0.0, 0.5 * gam, 2.0 * gam):
        local = 1.0 + a * gam ** 2 / (dlt ** 2 + gam ** 2)
        assert ntilde_lorentz(w, a, gam, dlt) == pytest.approx(local, rel=2e-2)
    # 窄峰 Gam << W：增强 ~ a*Gam*W/(dlt^2+W^2)，正比于 Gam
    gam_small = 0.001 * w
    e1 = ntilde_lorentz(w, a, gam_small, 0.0) - 1.0
    e2 = ntilde_lorentz(w, a, 2 * gam_small, 0.0) - 1.0
    assert e2 / e1 == pytest.approx(2.0, rel=1e-2)


# ==========================================================================
# C. 拟合机器本身（喂精确模型进去必须完美还原）
# ==========================================================================

def test_fit_machinery_is_exact_on_synthetic_lorentzian():
    """把一个**精确**的 (4) 式曲线喂给拟合器：必须还原 W_eff、R2=1、dlt0=0、A/(a*Gam)=1。

    这条测试的意义：它保证后面测出来的 R2 明显低于 1 是**物理结论**（E(dlt) 不是
    单洛伦兹），不是拟合器的 bug。
    """
    a, gam, w_true = 2.0, 100.0, 260.0
    d = np.linspace(-600.0, 600.0, 41)
    y = 1.0 + a * gam * w_true / (d ** 2 + w_true ** 2)
    f1 = fit_weff_1p(d, y, a, gam)
    assert f1["w_eff"] == pytest.approx(w_true, rel=1e-8)
    assert f1["r2"] == pytest.approx(1.0, abs=1e-12)
    f3 = fit_lorentz_3p(d, y, a, gam)
    assert f3["amp_ratio"] == pytest.approx(1.0, rel=1e-6)
    assert f3["delta0"] == pytest.approx(0.0, abs=1e-6)
    assert f3["r2"] == pytest.approx(1.0, abs=1e-12)
    wd = width_diagnostics(d, y, a, gam)
    for key in ("w_amp", "w_curv"):
        assert wd[key] == pytest.approx(w_true, rel=1e-6), key
    # W_half 用相邻两点线性插值定半高点，精度受 dlt 网格间距限制（这里 30 meV -> ~0.2%）
    assert wd["w_half"] == pytest.approx(w_true, rel=5e-3)
    assert e_analytic(d, a, gam, w_true) == pytest.approx(y, rel=1e-14)


# ==========================================================================
# D. 求解器侧：E 的基本性质
# ==========================================================================

def test_flat_dos_gives_E_exactly_one(maps):
    """常数（a=0）DOS 时 E=1：峰态与参考态走同一条数值路径，比值必须为 1。

    残差只来自 TcMap 的 spline 反解（本文件用 n_grid=12 的粗映射跑得快，实测
    ~4e-5；主脚本用 n_grid=22，抽查相对误差 ~1e-5，见 §0 诊断）。
    """
    st = Setup(lam_b=BASE.lam_b, mu_b=BASE.mu_b, w_log=BASE.w_log, a=0.0,
               gamma=BASE.gamma, n_eps=BASE.n_eps)
    tc, n0, clean = tc_peak(st, 137.0)
    assert clean and n0 == 1.0
    assert maps["kernel"].invert(tc) == pytest.approx(1.0, abs=1e-4)
    assert maps["lambda"].invert(tc) == pytest.approx(1.0, abs=1e-4)
    # 同时核对 tc_flat(scale=1) 就是这个 Tc（这一步与 spline 无关，必须精确）
    assert tc_flat(st, 1.0, "kernel") == pytest.approx(tc, rel=1e-9)


def test_large_delta_recovers_flat_dos(maps):
    """P4 §4 最后一行：dlt->inf（峰远离费米面）时 E->1。

    定量：(4) 式给 E-1 ~ a*Gam*W/dlt^2，dlt=4000 meV、W~250 meV 时约 3e-3。
    实测 E-1 在 dlt=4000 meV 处是 **-1e-4 量级（略低于 1，不是略高）**：窗口因子
    Ntilde(a_m) 是 m-依赖的（高 m 的相对增强更大），它对 Z 方程与 phi 方程的加权
    不同，所以"一个标量 E"无法精确代表 m-依赖的窗口因子；在增强本身已经衰减到
    ~1e-3 的远场，这个残差就成了主项。已核对不是数值假象：把
    dos_energy_cutoff_ratio 从 20 加到 200，Tc(peak)/Tc(flat)-1 从 -1.00e-4
    收敛到 -1.26e-4；dos_n_eps 从 400 加到 6400 也收敛。
    """
    e = [maps["kernel"].invert(tc_peak(BASE, d)[0]) for d in (1000.0, 2000.0, 4000.0)]
    assert e[0] > 1.01, f"dlt=1000 meV 时仍应有 ~1% 的正增强，实测 {e[0]}"
    assert abs(e[1] - 1.0) < 2e-3, f"dlt=2000 meV: E={e[1]}"
    assert abs(e[2] - 1.0) < 2e-3, f"dlt=4000 meV: E={e[2]}"
    assert abs(e[1] - 1.0) < abs(e[0] - 1.0), f"|E-1| 应先快速衰减: {e}"


def test_E_is_exactly_even_in_delta(maps):
    """对称峰形下 E(dlt)=E(-dlt) 严格成立（xi->-xi, dlt->-dlt, chi->-chi 的精确对称性），
    因此 E 的极大必然恰在 dlt=0，不需要拟合就能断言峰位。
    """
    for dlt in (60.0, 250.0):
        tp, _, _ = tc_peak(BASE, dlt)
        tm, _, _ = tc_peak(BASE, -dlt)
        assert tp == pytest.approx(tm, rel=1e-10)
        assert maps["kernel"].invert(tp) == pytest.approx(maps["kernel"].invert(tm), rel=1e-9)


def test_dE_ddelta_sign_and_peak_at_zero(maps):
    """P4 (6) 式 sign(dE/ddlt) = -sign(dlt)：E 沿 |dlt| 严格单调下降，极大在 dlt=0。"""
    deltas = [0.0, 80.0, 200.0, 400.0, 700.0]
    e = [maps["kernel"].invert(tc_peak(BASE, d)[0]) for d in deltas]
    assert e[0] == max(e)
    assert all(e[i] > e[i + 1] for i in range(len(e) - 1)), f"E 非单调: {e}"
    assert e[0] > 1.5, "a=2, Gam=100 meV 时 dlt=0 处应有可观增强"


# ==========================================================================
# E. 窗口平均恒等式 (2) 在求解器里是否逐项精确
# ==========================================================================

def test_window_average_identity_is_exact_once_chi_included():
    """**P4 §2 的窗口平均图像是精确的**，但窗口中心被 chi_m 挪走了。

    求解器的能量积分基元 G_m 与闭式
        G_m = Ntilde(a_m; dlt - chi_m) / (n0 * a_m),  a_m = |w_m| Z_m
    在 <1e-4 相对精度上一致（残差是 e_max 有限积分范围的截断 + 求积节点数）。
    **用 dlt 而不是 dlt - chi_m 会差 3%-6%**——P4 §2 的推导里没有 chi_m，
    这是 v0.1 要补的一项。这里用 dos_n_eps=400（脚本生产值），n_eps=200 时
    残差退化到 ~1.6e-4。
    """
    st = Setup(lam_b=BASE.lam_b, mu_b=BASE.mu_b, w_log=BASE.w_log, a=BASE.a,
               gamma=BASE.gamma, n_eps=400)
    r0 = window_identity_check(st, 0.0)
    assert r0["max_abs_chi"] == pytest.approx(0.0, abs=1e-9), "dlt=0 时 DOS 对称，chi 必为 0"
    assert r0["rel_err_with_chi"] < 1e-4

    r = window_identity_check(st, 300.0)
    assert r["rel_err_with_chi"] < 1e-4, "含 chi 的窗口平均恒等式应精确成立"
    assert r["max_abs_chi"] > 10.0, "dlt=300 meV 时 chi 应是几十 meV 量级，不可忽略"
    assert r["rel_err_no_chi"] > 20 * r["rel_err_with_chi"], (
        "忽略 chi 应显著变差；实测 no_chi ~6e-2 vs with_chi ~2e-5")


def test_flat_window_weights_are_normalized():
    """一阶微扰权重 w_m 必须归一（sum_m w_m = 1），否则 ntilde_multi 的绝对值无意义。"""
    ww = flat_window_weights(BASE, scale=1.0)
    assert ww["w_sum"] == pytest.approx(1.0, rel=1e-8)
    assert np.all(ww["a_m"] > 0)
    # 权重均值给出的"c" = <a_m>_w / w_log 应是 O(1) 的几
    assert 0.3 < ww["c_mean"] < 5.0, ww["c_mean"]
    # 多窗口预言在 dlt->inf 时回到 1
    assert float(ntilde_multi(ww["a_m"], ww["w_m"], BASE, np.array([1e7]))[0]) == \
        pytest.approx(1.0, abs=1e-6)


# ==========================================================================
# F. 单洛伦兹形式 (4) 的失效幅度
# ==========================================================================

def test_single_lorentzian_form_fails_quantitatively(maps):
    """任务书要求断言"洛伦兹拟合 R2>0.99"——**实测不成立，本测试如实锁定失效幅度**。

    物理原因（reviewer 也独立指出过）：Tc 方程里每个 Matsubara 项 m 各自带一个宽度
    a_m = |w_m|Z_m 的窗，E(dlt) 是**一族宽度不同的洛伦兹的加权和**，尾部被最宽的窗
    主导，因此比单个洛伦兹重。表现为：
      - 由峰值定的宽度 W_amp = a*Gam/(E(0)-1) 明显小于由曲率定的 W_curv；
      - 振幅锁定为 a*Gam 的 1-参数拟合 R2 只有 ~0.8-0.9（远低于 0.99）；
      - 放开振幅的 3-参数拟合要用 ~1.5 倍的振幅才能同时拟合宽度。
    这些是 P4 v0.1 必须记录的定量结论，不是数值噪声（E_num 的噪声由偶性残差
    定出在 ~1e-14 量级，见主脚本 §2）。
    """
    deltas = np.array([-600., -400., -250., -150., -80., 0., 80., 150., 250., 400., 600.])
    sc = scan_delta(BASE, deltas, maps=maps, verbose=False)
    e = sc["e_kernel"]
    assert np.all(sc["clean"]), "本组参数下 (Z,chi) 自洽不应报不收敛"

    f1 = fit_weff_1p(sc["delta"], e, BASE.a, BASE.gamma)
    f3 = fit_lorentz_3p(sc["delta"], e, BASE.a, BASE.gamma)
    wd = width_diagnostics(sc["delta"], e, BASE.a, BASE.gamma)

    # (i) 单洛伦兹 1-参数拟合达不到 R2>0.99
    assert f1["r2"] < 0.99, (
        f"R2={f1['r2']:.4f}：若这条断言开始失败，说明 E(dlt) 变成了单洛伦兹，"
        "应重新检查 (4) 式的适用范围并更新 P4 §5 的结论")
    assert 0.6 < f1["r2"] < 0.98
    # (ii) 三种宽度读法互相矛盾，比值明显偏离 1
    assert wd["curv_over_amp"] > 1.3, wd
    assert wd["half_over_amp"] > 1.2, wd
    # (iii) 放开振幅后能拟合得多，但振幅要显著超过 (4) 式预言的 a*Gam
    assert f3["r2"] > f1["r2"] + 0.05
    assert f3["amp_ratio"] > 1.15, f3
    # (iv) 峰位仍严格在 dlt=0（对称峰形的精确对称性，与洛伦兹形式无关）
    assert abs(f3["delta0"]) < 1.0, f3


def test_c_is_order_one_but_not_a_constant(maps):
    """标定的 c = (W_eff - Gam)/w_log 是 O(1) 的量，但随 lam_b 明显变化——
    因为窗宽 a_m = |w_m| Z_m 本身正比于 T_c 与 Z~1+lam，不是一个材料无关的常数。
    """
    deltas = np.array([0.0, 150.0, 350.0, 700.0])
    e = np.array([maps["kernel"].invert(tc_peak(BASE, d)[0]) for d in deltas])
    f1 = fit_weff_1p(deltas, e, BASE.a, BASE.gamma)
    c = (f1["w_eff"] - BASE.gamma) / BASE.w_log
    assert 0.2 < c < 4.0, f"c={c:.3f} 应是 O(1)"
    # 单点标定（只用 dlt=0，(4) 式在峰值处）应与拟合值同量级
    c_1pt = (BASE.a * BASE.gamma / (e[0] - 1.0) - BASE.gamma) / BASE.w_log
    assert 0.2 < c_1pt < 4.0
    assert abs(c - c_1pt) < 1.0, (c, c_1pt)
