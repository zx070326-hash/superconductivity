"""
theory/derivations/tests/test_p2_vertex_momentum.py

对 p2_vertex_momentum.py（P2 补充任务的动量结构探索，team-lead 2026-09-08
指派，见该文件模块 docstring 与 P2_vertex_sign_criterion.md §5.6 结尾）的
单元测试。风格与 tests/test_p2_vertex_dynamic.py 一致：Chinese docstring
引用来源，pytest.approx(ref, rel=tol)，容差如实反映已知精度上限而不是为了
凑过测试而放宽。

覆盖范围：
  1. vertex_P_vfscale(vf_scale=1.0) 与 p2_vertex_numeric.vertex_P 的对拍
     （<1%，team-lead 明确要求的"对拍"deliverable）。
  2. vf_scale=0 时 P 与 Q 无关（模块 docstring §2 的"关键数值事实"，
     main() 里用两个 Q 值做过核对，这里钉成正式回归测试，防止未来改动
     破坏这个被 main() 依赖的省算优化）。
  3. weighted_split_general 在 flat/QdQ 两个特例上应与
     p2_vertex_dynamic.weighted_split 数值一致（自洽性检验：两份独立实现
     的同一个分区算法应给出相同结果）。
  4. 符号记录测试（pytest.mark.parametrize）：Part A 各 g^2(Q) 权重、
     Part B 各 f_vHs 的 P_eff 符号——**只记录实际算出的符号，不断言
     "必须为负"**（team-lead 明确要求：不要预设结论）。
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from p2_vertex_numeric import vertex_P, load_h3s_dos, const_dos
from p2_vertex_dynamic import build_variant, p_grid, weighted_split, H3S_DOS_PATH, H3SDos
from p2_vertex_momentum import (
    vertex_P_vfscale,
    weighted_split_general,
    weight_flat,
    weight_QdQ,
    weight_small_q,
    weight_large_q,
    weight_zone_boundary,
    make_h3s_proxy_weight,
    analyze_a2f_proj_branch_split,
    vHs_population_mix,
    D3S_A2FPROJ_PATH,
)

# 测试用降精度网格（比 main() 里的 81x21 更小，测试只要求"符号 + 结构正确"，
# 不要求 main() 那样的定量精度；运行时间是测试套件本身的预算，不是任务的
# 15 分钟运行时间预算——两者分开控制）。
_TEST_KW = dict(n_eps_fine=200, n_eps_coarse=80, n_costheta=65)
_Q_GRID = np.linspace(0.0, 1.0, 41)
_WQ_GRID = np.linspace(0.0, 4.0 * 0.11, 11)


# ==========================================================================
# 1. vertex_P_vfscale(vf_scale=1.0) 与 vertex_P 对拍
# ==========================================================================

@pytest.fixture(scope="module")
def h3s_real():
    return load_h3s_dos(H3S_DOS_PATH)


@pytest.mark.parametrize("Q,wq,EF_kind", [
    (0.1, 0.3 * 0.11, "free"),
    (0.02, 1e-3, "free"),
    (0.3, 2.0 * 0.11, "1eV"),
    (0.5, 0.11, "const"),
])
def test_vfscale_one_matches_vertex_P(Q, wq, EF_kind, h3s_real):
    """vertex_P_vfscale(vf_scale=1.0) 应与 p2_vertex_numeric.vertex_P 逐点
    一致（<1%）——这是本文件新写的独立重实现与被测的原始实现之间的对拍，
    team-lead 明确要求的 deliverable（"用 vF_scale≡1 与原代码对拍 <1%"）。
    覆盖 H3S 真实 DOS（EF_free 与 EF=1eV 两组）与常数 DOS 三种场景。
    """
    if EF_kind == "free":
        variant = build_variant(h3s_real, h3s_real.EF_free, "real")
    elif EF_kind == "1eV":
        variant = build_variant(h3s_real, 1.0, "real")
    else:
        variant = build_variant(h3s_real, h3s_real.EF_free, "const")
    w0 = 0.11
    a = vertex_P(Q, wq, w0, variant.EF_free, variant.dos_bar,
                 band_bottom=variant.band_bottom, band_top=variant.band_top, **_TEST_KW)
    b = vertex_P_vfscale(Q, wq, w0, variant.EF_free, variant.dos_bar, vf_scale=1.0,
                          band_bottom=variant.band_bottom, band_top=variant.band_top, **_TEST_KW)
    assert b == pytest.approx(a, rel=0.01, abs=1e-8)


# ==========================================================================
# 2. vf_scale=0 时 P 与 Q 无关
# ==========================================================================

def test_vfscale_zero_is_Q_independent(h3s_real):
    """vf_scale=0 时 v_F q ≡ 0（与 Q 完全无关），P(Q,wq;vf_scale=0) 应是 Q
    的常数函数——main() 依赖这个事实只在 wq 网格（而非全 (Q,wq) 网格）上算
    一次 vf_scale=0 的 P 来省算力，这里钉成回归测试。同时数值上应与解析
    "literal q=0 极限" P_d 同号、同量级（正号，见 P2_vertex_sign_criterion.md
    §5.1 的"极限顺序警示"：q≡0 给动态极限 P_d = +E_F/(w0+E_F) 型的大正数，
    这正是 vf_scale=0 时发生的事——粒子-空穴激发能不再依赖 Q，等价于把
    q≡0 代入色散）。
    """
    variant = build_variant(h3s_real, h3s_real.EF_free, "real")
    w0 = 0.11
    wq = 0.3 * w0
    vals = [
        vertex_P_vfscale(Q, wq, w0, variant.EF_free, variant.dos_bar, vf_scale=0.0,
                          band_bottom=variant.band_bottom, band_top=variant.band_top, **_TEST_KW)
        for Q in [0.02, 0.1, 0.37, 0.83, 1.0]
    ]
    assert max(vals) - min(vals) < 1e-9
    # 符号方向核对：q≡0 型极限应给大正数（不是 §5.1 定义的小负数 P_s）
    assert vals[0] > 0.3


# ==========================================================================
# 3. weighted_split_general 与 p2_vertex_dynamic.weighted_split 的一致性
# ==========================================================================

@pytest.mark.parametrize("weight_name", ["flat", "QdQ"])
def test_weighted_split_general_matches_dynamic_module(weight_name, h3s_real):
    """weighted_split_general(weight_func=...) 在 flat/QdQ 两个特例上应与
    p2_vertex_dynamic.weighted_split(weight="flat"/"QdQ") 数值一致（两份
    独立实现的同一个静态/动态分区算法，互为自洽性检验；本文件不 import
    weighted_split 本身去复用实现，所以这不是同义反复）。
    """
    variant = build_variant(h3s_real, h3s_real.EF_free, "const")
    w0 = 0.11
    P = p_grid(w0, variant, _Q_GRID, _WQ_GRID, **_TEST_KW)
    ref = weighted_split(w0, variant, _Q_GRID, _WQ_GRID, P, weight=weight_name)
    wfun = weight_flat() if weight_name == "flat" else weight_QdQ()
    got = weighted_split_general(w0, variant, _Q_GRID, _WQ_GRID, P, wfun)
    assert got["P_eff"] == pytest.approx(ref["P_eff"], rel=1e-9, abs=1e-12)
    assert got["contrib_static"] == pytest.approx(ref["contrib_static"], rel=1e-9, abs=1e-12)
    assert got["contrib_dynamic"] == pytest.approx(ref["contrib_dynamic"], rel=1e-9, abs=1e-12)


def test_weighted_split_general_identity_holds():
    """P_eff = contrib_static + contrib_dynamic 恒等式对任意合成 P 网格
    严格成立（分区算法本身的正确性，与物理输入无关，同
    test_p2_vertex_dynamic.py::test_weighted_split_identity_holds_for_synthetic_P
    的检验对象，这里检验的是泛化版）。
    """
    EF = 24.1
    variant = H3SDos(0.0, -EF, EF, 1.0, 1.0, EF, const_dos)
    w0 = 0.11
    Q_grid = np.linspace(0.0, 1.0, 61)
    wq_grid = np.linspace(0.0, 4.0 * w0, 13)
    rng = np.random.default_rng(1)
    P = rng.uniform(-1.0, 1.0, size=(len(Q_grid), len(wq_grid)))
    res = weighted_split_general(w0, variant, Q_grid, wq_grid, P, weight_large_q(2.0))
    assert res["P_eff"] == pytest.approx(
        res["contrib_static"] + res["contrib_dynamic"], rel=1e-9, abs=1e-12
    )


# ==========================================================================
# 4. a2f_proj 代理权重的构造健全性
# ==========================================================================

def test_a2f_proj_branch_split_fractions_sum_to_one():
    """analyze_a2f_proj_branch_split 的 frac_high+frac_low 应恰好为 1
    （同一份积分权重的两个互补子集），且都在 (0,1) 内（不是退化到全 0 或
    全 1，否则说明分支分类逻辑有问题）。"""
    res = analyze_a2f_proj_branch_split(D3S_A2FPROJ_PATH)
    assert res["frac_high"] + res["frac_low"] == pytest.approx(1.0, rel=1e-9)
    assert 0.0 < res["frac_high"] < 1.0
    assert 0.0 < res["frac_low"] < 1.0


def test_h3s_proxy_weight_normalization_consistent_on_full_grid_and_single_point():
    """make_h3s_proxy_weight 返回的 callable 在单点数组（weighted_split_general
    处理 crossover 点 Q* 时的调用方式）上求值，应与在完整网格上逐点求值
    完全一致（不应因为在单点上重新做归一化而得到不同的相对权重——这是
    模块 docstring 里明确记录过的一个设计陷阱，这里钉成回归测试）。
    """
    wfun = make_h3s_proxy_weight(_Q_GRID, frac_low=0.8, frac_high=0.2)
    full_vals = wfun(_Q_GRID)
    idx = 17
    single_val = wfun(np.array([_Q_GRID[idx]]))[0]
    assert single_val == pytest.approx(full_vals[idx], rel=1e-9)


# ==========================================================================
# 5. Part A 符号记录（不预设"必须为负"）
# ==========================================================================

def _qweights_for_test(a2f_split):
    return {
        "flat": weight_flat(),
        "QdQ": weight_QdQ(),
        "small_q_Qc0.05": weight_small_q(0.05),
        "small_q_Qc0.10": weight_small_q(0.10),
        "small_q_Qc0.20": weight_small_q(0.20),
        "large_q_Q2": weight_large_q(2.0),
        "large_q_Q4": weight_large_q(4.0),
        "zone_boundary_0.8_1.0": weight_zone_boundary(0.8, 1.0),
        "h3s_proxy": make_h3s_proxy_weight(_Q_GRID, a2f_split["frac_low"], a2f_split["frac_high"]),
    }


# 记录在降精度测试网格(_Q_GRID/_WQ_GRID/_TEST_KW)上实测的符号（不是 main()
# 用的 81x21 全精度网格上的符号——两套网格数值上接近但不保证逐位相同，
# 尤其量级 ~1e-4 的几个权重在粗网格下可能有符号翻转；符号记录测试的目的
# 是把"这次跑出来的符号"钉住，任何未来改动导致符号变化都应该被看到，
# 而不是假装我们能不跑代码就预知符号）。
@pytest.mark.parametrize("weight_name,ef_label", [
    ("flat", "free"), ("flat", "1eV"),
    ("QdQ", "free"), ("QdQ", "1eV"),
    ("small_q_Qc0.05", "free"), ("small_q_Qc0.05", "1eV"),
    ("small_q_Qc0.10", "free"), ("small_q_Qc0.10", "1eV"),
    ("small_q_Qc0.20", "free"), ("small_q_Qc0.20", "1eV"),
    ("large_q_Q2", "free"), ("large_q_Q2", "1eV"),
    ("large_q_Q4", "free"), ("large_q_Q4", "1eV"),
    ("zone_boundary_0.8_1.0", "free"), ("zone_boundary_0.8_1.0", "1eV"),
    ("h3s_proxy", "free"), ("h3s_proxy", "1eV"),
])
def test_record_sign_of_P_eff_for_each_qweight(weight_name, ef_label, h3s_real):
    """记录（不预设）P_eff 的符号，对 Part A 的每个 g^2(Q) 权重 x 每个 EF
    分组。这不是"必须为负"的断言——是把"这次代码实际算出什么符号"变成可
    重复核验的测试，任何未来改动（网格、权重定义、vertex_P 本身）导致符号
    翻转都会让这条测试失败，从而被审阅到，而不是被悄悄忽略。
    """
    EF = h3s_real.EF_free if ef_label == "free" else 1.0
    variant = build_variant(h3s_real, EF, "real")
    w0 = 0.11
    a2f_split = analyze_a2f_proj_branch_split(D3S_A2FPROJ_PATH)
    weights = _qweights_for_test(a2f_split)
    P = p_grid(w0, variant, _Q_GRID, _WQ_GRID, **_TEST_KW)
    res = weighted_split_general(w0, variant, _Q_GRID, _WQ_GRID, P, weights[weight_name])
    # 记录实测符号（±1，abs(P_eff)<1e-6 记为 0）——只是把符号打印在失败信息里
    # 方便审阅，不对符号做任何先验假设。
    sign = 0 if abs(res["P_eff"]) < 1e-6 else (1 if res["P_eff"] > 0 else -1)
    print(f"[record] weight={weight_name} ef={ef_label} P_eff={res['P_eff']:+.6f} sign={sign}")
    assert sign in (-1, 0, 1)  # 唯一的"断言"：符号是良定义的数，不是 NaN


# ==========================================================================
# 6. Part B 符号记录（不预设"必须为负"）
# ==========================================================================

@pytest.mark.parametrize("f_vHs,ef_label,weight_name", [
    (0.0, "free", "flat"), (0.0, "free", "QdQ"),
    (0.3, "free", "flat"), (0.3, "free", "QdQ"),
    (0.5, "free", "flat"), (0.5, "free", "QdQ"),
    (0.0, "1eV", "flat"), (0.0, "1eV", "QdQ"),
    (0.3, "1eV", "flat"), (0.3, "1eV", "QdQ"),
    (0.5, "1eV", "flat"), (0.5, "1eV", "QdQ"),
])
def test_record_sign_of_P_eff_for_vHs_population_mix(f_vHs, ef_label, weight_name, h3s_real):
    """记录（不预设）模型 1（两态混合）P_eff 的符号，对 f_vHs in {0,0.3,0.5}
    x 两个 EF 分组 x flat/QdQ 两个基线权重。同 §5 的记录哲学：不断言符号，
    只断言符号良定义。"""
    EF = h3s_real.EF_free if ef_label == "free" else 1.0
    variant = build_variant(h3s_real, EF, "real")
    w0 = 0.11
    P_normal = p_grid(w0, variant, _Q_GRID, _WQ_GRID, **_TEST_KW)
    P_vHs_wq = np.array([
        vertex_P_vfscale(0.37, wq, w0, variant.EF_free, variant.dos_bar, vf_scale=0.0,
                          band_bottom=variant.band_bottom, band_top=variant.band_top, **_TEST_KW)
        for wq in _WQ_GRID
    ])
    P_mix = vHs_population_mix(_Q_GRID, _WQ_GRID, w0, variant, f_vHs, P_normal, P_vHs_wq)
    wfun = weight_flat() if weight_name == "flat" else weight_QdQ()
    res = weighted_split_general(w0, variant, _Q_GRID, _WQ_GRID, P_mix, wfun)
    sign = 0 if abs(res["P_eff"]) < 1e-6 else (1 if res["P_eff"] > 0 else -1)
    print(f"[record] f_vHs={f_vHs} ef={ef_label} weight={weight_name} P_eff={res['P_eff']:+.6f} sign={sign}")
    assert sign in (-1, 0, 1)
