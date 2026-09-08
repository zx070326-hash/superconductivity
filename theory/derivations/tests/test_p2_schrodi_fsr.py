"""
theory/derivations/tests/test_p2_schrodi_fsr.py

对 p2_schrodi_fsr.py 的三类检验：
  1. P=0（或两个 vertex 开关都 False）时，SchrodiFSRSolver.tc() 必须与
     eliashberg_iso.Solver.tc()（折叠正频算法）逐位一致，容差 <0.1 K——
     两者是对同一组标准各向同性 FSR Migdal 线性化方程的两种独立数值实现
     （全轴 vs 折叠正频轴），一致性是最基本的正确性检验。
  2. 顶角核的对称性/量纲自洽检验（独立重新推导，不从被测模块 import 参考值）：
       - Z_m^(0)（含顶角修正）在折叠反射 m -> -1-m 下不变（这是折叠假设
         Z_{-1-n}=Z_n 成立所需要的不变量；E(m) 本身是反对称的，不是对称的，
         见 p2_schrodi_fsr.py 模块 docstring §1a 末尾的更正说明）。
       - K^V_Δ[m,k] 在同时反射 (m,k)->(M-1-m,M-1-k) 下不变。
  3. Einstein 谱下，P 很小时 δ=ln(Tc^V/Tc^0) 与 P 成正比（线性响应区），
     用有限差分交叉验证斜率的收敛性（不依赖 P_3D 的具体取值是否精确，只
     验证线性化实现本身在小 P 下行为正确）。

按本目录 test_p2_vertex.py 的约定：sys.path 插入被测模块所在目录；参考公式
在测试文件内独立写一遍，不从被测模块 import，避免共享同一份（可能有错的）
公式导致假阳性。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "solvers"))

from p2_schrodi_fsr import (
    SchrodiFSRSolver,
    vertex_kernels,
    kv_delta,
    _toeplitz,
    _sign_array,
    p3d_dimensional_estimate,
)
from eliashberg_iso import Spectrum, Solver


# ==========================================================================
# 1. P=0 基准：与 eliashberg_iso.Solver.tc() 逐位一致（<0.1 K）
# ==========================================================================

@pytest.mark.parametrize("lam,w0,mu_star,w_cut", [
    (0.8, 40.0, 0.12, 500.0),
    (1.5, 80.0, 0.16, 800.0),
])
def test_p0_matches_reference_solver_einstein(lam, w0, mu_star, w_cut):
    sp = Spectrum.einstein(lam=lam, w0=w0, name="e")
    tc_mine = SchrodiFSRSolver(sp, mu_star, w_cut, P=0.0).tc()
    tc_ref = Solver(sp, mu_star, w_cut).tc()
    assert abs(tc_mine - tc_ref) < 0.1, f"mine={tc_mine}, ref={tc_ref}"


def test_p0_matches_reference_solver_vertex_flags_off():
    """P!=0 但两个 vertex_in_* 都 False 时也应精确退化为标准 Migdal 方程。"""
    sp = Spectrum.einstein(lam=1.0, w0=60.0, name="e")
    tc_mine = SchrodiFSRSolver(
        sp, 0.14, 600.0, P=1.0, vertex_in_Z=False, vertex_in_Delta=False
    ).tc()
    tc_ref = Solver(sp, 0.14, 600.0).tc()
    assert abs(tc_mine - tc_ref) < 0.1


# ==========================================================================
# 2. 顶角核的对称性/自洽检验（独立重新推导）
# ==========================================================================

def _build_kernels(lam_of_nu, n_pos, t_mev):
    """独立（不依赖 SchrodiFSRSolver 内部方法）搭建 idx/wn/s/lam/C/D/E，供本文件的
    检验函数直接验证 vertex_kernels() 与 kv_delta() 的输出。"""
    m = 2 * n_pos
    idx = np.arange(-n_pos, n_pos)
    wn = np.pi * t_mev * (2 * idx + 1)
    lam_pos = lam_of_nu(2 * np.pi * t_mev * np.arange(m))
    lam, C, D, E = vertex_kernels(lam_pos, n_pos)
    return idx, wn, lam, C, D, E


def test_E_is_antisymmetric_under_fold_reflection():
    """独立推导（见测试文件顶部说明）：E(-1-m)=-E(m)，不是 E(-1-m)=E(m)。

    换元 m'=-1-a, m''=-1-b，用 λ(x)=λ(-x)（lambda_of_nu 是偶函数）与
    s_{-1-x}=-s_x，逐项代入 E(m)=Σλ(m-m')λ(m'-m'')s_{m'}s_{m''}s_{m''-m'+m}
    可得 E(-1-m) = -E(m)（详细代数见 p2_schrodi_fsr.py 模块 docstring §1a）。
    """
    sp = Spectrum.einstein(lam=0.9, w0=50.0, name="e")
    n_pos = 30
    _, _, _, _, _, E = _build_kernels(sp.lambda_of_nu, n_pos, t_mev=7.0)
    E_expected_antisym = -E[::-1]
    assert np.allclose(E, E_expected_antisym, atol=1e-9), "E(m) 应满足 E(-1-m)=-E(m)"
    # 顺带确认它不是（错误地）对称的，避免两个测试互相掩盖对方的失败
    assert not np.allclose(E, E[::-1], atol=1e-6)


def test_Z0_symmetric_under_fold_reflection():
    """物理不变量：Z_m^(0)（含顶角修正）在 m->-1-m 下不变，
    是折叠假设 Z_{-1-n}=Z_n 成立所需要的组合（标准项 S(m)/ω_m 与顶角项
    E(m)/ω_m 各自反对称，但除以同样反对称的 ω_m 后都变回对称，见
    p2_schrodi_fsr.py 模块 docstring §1a 的更正说明）。
    """
    sp = Spectrum.einstein(lam=0.9, w0=50.0, name="e")
    n_pos = 30
    t_mev = 7.0
    P = 3e-4
    idx, wn, lam, C, D, E = _build_kernels(sp.lambda_of_nu, n_pos, t_mev)
    s = _sign_array(n_pos)
    z0 = 1.0 + (np.pi * t_mev / wn) * (lam @ s) - P * (np.pi ** 3 * t_mev ** 2 / wn) * E
    assert np.allclose(z0, z0[::-1], atol=1e-9), "Z_m^(0) 应在 m->-1-m 反射下不变"


def test_kv_delta_symmetric_under_joint_fold_reflection():
    """K^V_Δ[m,k] 在同时反射 (m,k)->(M-1-m,M-1-k) 下不变。"""
    sp = Spectrum.einstein(lam=0.9, w0=50.0, name="e")
    n_pos = 30
    _, _, lam, C, D, _ = _build_kernels(sp.lambda_of_nu, n_pos, t_mev=7.0)
    KV = kv_delta(lam, C, D)
    assert np.allclose(KV, KV[::-1, ::-1], atol=1e-9)


def test_toeplitz_matches_direct_lambda_lookup():
    """_toeplitz 的位置差=值差 简化本身的正确性（与直接按定义构造的矩阵比较）。"""
    sp = Spectrum.einstein(lam=0.9, w0=50.0, name="e")
    n_pos = 12
    m = 2 * n_pos
    t_mev = 5.0
    lam_pos = sp.lambda_of_nu(2 * np.pi * t_mev * np.arange(m))
    lam = _toeplitz(lam_pos, m)
    idx = np.arange(-n_pos, n_pos)
    direct = np.array([[sp.lambda_of_nu(2 * np.pi * t_mev * abs(idx[i] - idx[j]))
                         for j in range(m)] for i in range(m)])
    assert np.allclose(lam, direct, rtol=1e-10)


# ==========================================================================
# 3. 小 P 线性响应检验（Einstein 谱，有限差分斜率收敛）
# ==========================================================================

@pytest.mark.parametrize("vz,vd", [(True, False), (False, True), (True, True)])
def test_small_P_linear_response(vz, vd):
    """δ=ln(Tc^V/Tc^0) 在 P->0 时应与 P 成正比：两个不同小步长算出的差分斜率
    应在 <1% 相对差内收敛（不依赖 P 的具体物理取值，只检验线性化实现本身）。
    """
    sp = Spectrum.einstein(lam=0.8, w0=40.0, name="e")
    mu_star, w_cut = 0.12, 500.0
    tc0 = SchrodiFSRSolver(sp, mu_star, w_cut, P=0.0).tc()

    def eig_at_tc0(P):
        return SchrodiFSRSolver(
            sp, mu_star, w_cut, P=P, vertex_in_Z=vz, vertex_in_Delta=vd, n_pos_cap=60
        ).max_eigen(tc0)

    e0 = eig_at_tc0(0.0)
    eps1, eps2 = 1e-6, 2e-6
    slope1 = (eig_at_tc0(eps1) - e0) / eps1
    slope2 = (eig_at_tc0(eps2) - e0) / eps2
    rel_diff = abs(slope1 - slope2) / max(abs(slope1), 1e-12)
    assert rel_diff < 0.01, f"斜率不收敛: slope(eps)={slope1}, slope(2eps)={slope2}"


def test_both_equals_sum_of_channels_in_linear_regime():
    """线性区内，'both' 配置的 δ 应约等于 Z_only 与 Delta_only 两个 δ 之和
    （线性响应的可加性：两个顶角开关独立作用于矩阵的不同部分，P 一次项下
    互不耦合）。

    n_pos_cap 必须充分覆盖该 (w_cut, Tc) 组合下的自然 Matsubara 窗口
    （由 _n_pos_for_T 给出，这里约 150），否则三个配置各自被 bisection
    找到的不同 T 会截到不同程度，截断偏差不会在三者之间抵消，人为放大
    "不可加性"（调试记录：cap=60 时相对偏差 89%，cap=200 时降到 4e-5，
    是截断而非物理，见模块 docstring 关于 n_pos_cap 的说明）。
    """
    sp = Spectrum.einstein(lam=0.8, w0=40.0, name="e")
    mu_star, w_cut = 0.12, 500.0
    tc0 = SchrodiFSRSolver(sp, mu_star, w_cut, P=0.0).tc()
    P = 1e-6

    def delta_for(vz, vd):
        tc = SchrodiFSRSolver(
            sp, mu_star, w_cut, P=P, vertex_in_Z=vz, vertex_in_Delta=vd, n_pos_cap=220
        ).tc(t_bracket_k=(0.8 * tc0, 1.3 * tc0))
        return float(np.log(tc / tc0))

    d_z = delta_for(True, False)
    d_d = delta_for(False, True)
    d_both = delta_for(True, True)
    assert abs(d_both - (d_z + d_d)) < 0.05 * max(abs(d_both), 1e-8)


def test_p3d_dimensional_estimate_positive_and_scales_with_inverse_EF():
    assert p3d_dimensional_estimate(1000.0) == pytest.approx(1.0 / 4000.0)
    assert p3d_dimensional_estimate(500.0) == 2.0 * p3d_dimensional_estimate(1000.0)
    with pytest.raises(ValueError):
        p3d_dimensional_estimate(-1.0)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
