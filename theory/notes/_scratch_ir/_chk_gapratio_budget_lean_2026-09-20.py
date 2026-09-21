"""
任务 0 预算（轻量版）：2Delta/kTc 判别性比值。

策略：用 Allen-Dynes 快速估计 Tc，仅在基准点做一次完整 N 收敛检查，
扰动点用固定的较小 N（N=1600）快速评估梯度。严格性明确标注为
numerical-evidence（非证书级），但足以判断「值不值得推进任务1/2」。

判别性比值 R = |dO/d(lambda)*Dlambda| / |dO/d(mu*)*Dmu*|
R<1 即被干扰量压住。
"""
import sys
sys.path.insert(0, r"D:\physics\superconductivity")
import os
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K


def gap_ratio_fast(lam, mu_star, w0_mev, w_cut_mev, t_low_frac=0.03):
    """用 Allen-Dynes 快速估 Tc，固定 N=1600 求 gap（H3S 量级下已够用）"""
    spec = Spectrum.einstein(lam=lam, w0=w0_mev)
    tc = spec.allen_dynes(mu_star)  # 快速闭式，跳过 eigensolve
    t_low = t_low_frac * tc
    solver = Solver(spectrum=spec, mu_star=mu_star, w_cut=w_cut_mev, n_matsubara=1600)
    wn, delta, z, conv = solver.solve_gap(t_low, mixing=0.25, max_iter=3000)
    ratio = 2.0 * delta[0] / (KB_MEV_PER_K * tc)
    return tc, ratio, conv


if __name__ == "__main__":
    lam0 = 1.876
    w0 = 106.9
    w_cut = 500.0
    mu0 = 0.13

    print("=== 基准点 N 收敛快速核验（只算三点：800/1600/3200）===")
    spec0 = Spectrum.einstein(lam=lam0, w0=w0)
    tc_base = spec0.allen_dynes(mu0)
    print(f"Tc(AD 估计, 基准)={tc_base:.3f} K")

    for N in [800, 1600, 3200]:
        solver = Solver(spectrum=spec0, mu_star=mu0, w_cut=w_cut, n_matsubara=N)
        t_low = 0.03 * tc_base
        wn, delta, z, conv = solver.solve_gap(t_low, mixing=0.25, max_iter=3000)
        r = 2.0 * delta[0] / (KB_MEV_PER_K * tc_base)
        print(f"  N={N:4d}: 2D(iw0)/kTc={r:.4f}  [收敛: {'是' if conv else '否'}]")

    print("\n采用 N=1600 用于全部扰动点（基于上述收敛性判断）")

    print("\n=== mu* 扰动 (Dmu*=0.05) ===")
    dmu = 0.05
    tc_mu_p, r_mu_p, conv_p = gap_ratio_fast(lam0, mu0 + dmu, w0, w_cut)
    tc_mu_m, r_mu_m, conv_m = gap_ratio_fast(lam0, mu0 - dmu, w0, w_cut)
    print(f"  mu*={mu0+dmu:.2f}: Tc={tc_mu_p:.3f} K, O={r_mu_p:.4f}  [收敛: {'是' if conv_p else '否'}]")
    print(f"  mu*={mu0-dmu:.2f}: Tc={tc_mu_m:.3f} K, O={r_mu_m:.4f}  [收敛: {'是' if conv_m else '否'}]")
    dO_dmu = 0.5 * (r_mu_p - r_mu_m)
    print(f"  中心差分 dO*dmu = {dO_dmu:.5f}")

    print("\n=== lambda 扰动 (Dlam/lam=0.15) ===")
    dlam_frac = 0.15
    lam_p = lam0 * (1 + dlam_frac)
    lam_m = lam0 * (1 - dlam_frac)
    tc_lam_p, r_lam_p, conv_lp = gap_ratio_fast(lam_p, mu0, w0, w_cut)
    tc_lam_m, r_lam_m, conv_lm = gap_ratio_fast(lam_m, mu0, w0, w_cut)
    print(f"  lam={lam_p:.4f}: Tc={tc_lam_p:.3f} K, O={r_lam_p:.4f}  [收敛: {'是' if conv_lp else '否'}]")
    print(f"  lam={lam_m:.4f}: Tc={tc_lam_m:.3f} K, O={r_lam_m:.4f}  [收敛: {'是' if conv_lm else '否'}]")
    dO_dlam = 0.5 * (r_lam_p - r_lam_m)
    print(f"  中心差分 dO*dlam = {dO_dlam:.5f}")

    print("\n=== 判别性比值 ===")
    R = abs(dO_dlam) / abs(dO_dmu)
    print(f"  R = |dO/dlam*Dlam| / |dO/dmu*Dmu| = {R:.4f}")
    if R < 1:
        verdict = "FAIL: 被 mu* 干扰压住（同位素指数式下场）"
    elif R < 1.5:
        verdict = "MARGINAL: 信号与干扰同量级，需任务1实轴修正后重判"
    else:
        verdict = "PASS: 机制轴信号显著强于 mu* 干扰，值得推进"
    print(f"  判决: {verdict}")

    print("\n=== Carbotte 公式交叉核对 ===")
    def carbotte_O(tc_k, wlog_mev):
        t = KB_MEV_PER_K * tc_k / wlog_mev
        return 3.53 * (1 + 12.5 * t**2 * np.log(wlog_mev / (2 * KB_MEV_PER_K * tc_k)))

    O_mu_p_c = carbotte_O(tc_mu_p, w0)
    O_mu_m_c = carbotte_O(tc_mu_m, w0)
    O_lam_p_c = carbotte_O(tc_lam_p, w0)
    O_lam_m_c = carbotte_O(tc_lam_m, w0)
    dO_dmu_c = 0.5 * (O_mu_p_c - O_mu_m_c)
    dO_dlam_c = 0.5 * (O_lam_p_c - O_lam_m_c)
    R_c = abs(dO_dlam_c) / abs(dO_dmu_c)
    print(f"  Carbotte: dO*dmu={dO_dmu_c:.5f}, dO*dlam={dO_dlam_c:.5f}, R={R_c:.4f}")
    print(f"  一致性检查: Eliashberg vs Carbotte 比值差异 = {abs(R - R_c):.3f}")
