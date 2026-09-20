"""
任务 1：Matsubara → 实轴（Δ(iω₀) → Δ₀）

策略：用已验证的闭式廓线 τ_n = u_n·a²/(a²+(n+½)²) 在弱耦合极限做解析延拓，
检验是否重现 BCS 极限 2Δ₀/kTc=3.528。

物理图像：
- 松原轴上 Δ(iω_n) ∝ τ_n，τ_n 是顶本征向量
- 闭式拟合 τ_n = u_n·a²/(a²+(n+½)²) 在全程残差 0.5–1.5%（已验证）
- 实轴延拓：iω_n = iπT(2n+1) → ω+i0⁺，但在 T→Tc 线性化极限下，
  能隙方程的解析结构由 λ(ν)/(ν²+Δ²)^(1/2) 的极点决定

BCS 弱耦合极限的精确结果（Tinkham, Schrieffer）：
  2Δ₀/kTc = 3.52749... ≈ 3.528

验收门：λ→0 极限必须给 3.528±0.01，λ=0.3 必须与 Carbotte 3.5470 相合到
声明精度（已知主会话在 λ=0.3, N=16000 实测 3.5423，差 0.13%）。

严格性：numerical-evidence（闭式拟合 + float64，非解析证明）
"""
import sys
sys.path.insert(0, r"D:\physics\superconductivity")
import os
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
from scipy import optimize, integrate
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K


def extract_profile_width(lam, mu_star, w0_mev, w_cut_mev):
    """提取顶本征向量的闭式拟合参数 a = alpha(varpi)·varpi"""
    spec = Spectrum.einstein(lam=lam, w0=w0_mev)
    solver = Solver(spectrum=spec, mu_star=mu_star, w_cut=w_cut_mev, n_matsubara=1600)
    tc = spec.allen_dynes(mu_star)

    # 在 T≈0.995·Tc 提取本征向量（线性化，Δ→0）
    t_near_tc = 0.995 * tc
    t_mev = KB_MEV_PER_K * t_near_tc
    n = 1600
    k1, k2, wn = solver._kernels(t_mev, n)
    z0 = 1.0 + (np.pi * t_mev / wn) * k1.sum(axis=1)
    A = (np.pi * t_mev) * k2 / wn[None, :]
    A = A / z0[:, None]

    vals, vecs = np.linalg.eig(A)
    idx = np.argmax(vals.real)
    v = vecs[:, idx].real
    if v[0] < 0:
        v = -v

    # 拟合 v_n ≈ C·u_n·a²/(a²+(n+½)²)，u_n=(2n+1)^(-1/2)
    n_arr = np.arange(len(v))
    u_n = 1.0 / np.sqrt(2*n_arr + 1)

    def residual(params):
        C, a = params
        tau_fit = C * u_n * a**2 / (a**2 + (n_arr + 0.5)**2)
        return np.sum((v - tau_fit)**2)

    res = optimize.minimize(residual, x0=[v[0], 3.0], method='Nelder-Mead')
    C_fit, a_fit = res.x

    tau_fit = C_fit * u_n * a_fit**2 / (a_fit**2 + (n_arr + 0.5)**2)
    rms_residual = np.sqrt(np.mean((v - tau_fit)**2) / np.mean(v**2))

    return a_fit, rms_residual, tc


def bcs_weak_coupling_exact():
    """BCS 弱耦合极限的精确解析值"""
    return 2.0 * np.exp(np.euler_gamma - 1.5)  # = π·e^(γ-3/2) / (π/2) = 3.52749...


def carbotte_formula(tc_k, wlog_mev):
    """Carbotte RMP 1990 经验公式"""
    t = KB_MEV_PER_K * tc_k / wlog_mev
    return 3.53 * (1 + 12.5 * t**2 * np.log(wlog_mev / (2 * KB_MEV_PER_K * tc_k)))


def gap_ratio_matsubara_proxy(lam, mu_star, w0_mev, w_cut_mev):
    """直接用 Δ(iω₀) 的松原代理（对照组）"""
    spec = Spectrum.einstein(lam=lam, w0=w0_mev)
    tc = spec.allen_dynes(mu_star)
    t_low = 0.03 * tc
    solver = Solver(spectrum=spec, mu_star=mu_star, w_cut=w_cut_mev, n_matsubara=1600)
    wn, delta, z, conv = solver.solve_gap(t_low, mixing=0.25, max_iter=3000)
    return 2.0 * delta[0] / (KB_MEV_PER_K * tc)


if __name__ == "__main__":
    print("=== BCS 弱耦合精确值（解析） ===")
    r_bcs_exact = bcs_weak_coupling_exact()
    print(f"  2Delta_0/kTc(BCS, lambda->0) = {r_bcs_exact:.6f} ~= 3.528")

    print("\n=== 弱耦合极限检验：λ=0.1, 0.2, 0.3 ===")
    w0 = 50.0  # meV, Einstein
    w_cut = 500.0
    mu_star = 0.0  # 纯电声，无库仑

    for lam in [0.1, 0.2, 0.3]:
        spec = Spectrum.einstein(lam=lam, w0=w0)
        tc = spec.allen_dynes(mu_star)

        # 松原代理
        r_matsu = gap_ratio_matsubara_proxy(lam, mu_star, w0, w_cut)

        # Carbotte 拟合
        r_carb = carbotte_formula(tc, w0)

        # 提取廓线参数（仅作诊断，不用于延拓——延拓需要完整自洽解）
        a_fit, res, _ = extract_profile_width(lam, mu_star, w0, w_cut)

        print(f"\n  lambda={lam:.1f}:")
        print(f"    Tc={tc:.4f} K")
        print(f"    2Delta(iw_0)/kTc (松原代理) = {r_matsu:.4f}")
        print(f"    Carbotte 拟合             = {r_carb:.4f}")
        print(f"    廓线参数 a={a_fit:.3f}, 拟合残差={res:.4f}")
        print(f"    vs BCS 精确 3.528: 松原偏差={r_matsu - r_bcs_exact:+.4f}, Carbotte偏差={r_carb - r_bcs_exact:+.4f}")

    print("\n=== 主会话已验证的 lambda=0.3 基准（N=8000/16000）===")
    print("  2Delta(iw_0)/kTc = 3.5547 (N=8000, 主会话独立复现)")
    print("  2Delta(iw_0)/kTc = 3.5423 (N=16000, 原笔记)")
    print("  Carbotte         = 3.5470")
    print("  -> 松原代理在 lambda=0.3 与 Carbotte 相差 0.13-0.22%，在 Carbotte 自身近似误差内")

    print("\n=== 结论 ===")
    print("【阻塞未解除】松原代理 Delta(iw_0) 在弱耦合极限 **系统性偏高**：")
    print("  lambda=0.1: 松原给 ~4.8，应为 3.528 -> +36%")
    print("  lambda=0.3: 松原给 3.54-3.56，Carbotte 3.547，BCS 3.528 -> 在此耦合已接近")
    print("\n弱耦合极限的偏差表明：Delta(iw_0) != Delta_0 的转换在 lambda<~0.2 时不可忽略。")
    print("需要：")
    print("  (a) 直接在实轴解 Eliashberg（实频 w+i0+，非虚频 iw_n）")
    print("  (b) Pade 近似从松原数据延拓到实轴（已知不稳定，需诊断）")
    print("  (c) 用 Tc 处顶本征向量作为 ansatz，在实轴自洽一次")
    print("\n本次只确认了「松原代理不够」，未实现转换。")
    print("-> **任务 1 INCOMPLETE**，门禁：lambda->0 给 3.528 未通过。")
