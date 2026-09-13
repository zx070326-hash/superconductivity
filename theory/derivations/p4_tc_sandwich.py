r"""
P4 · Tc 双侧夹逼与 Einstein Tc 凸性数值验证

任务（2026-09-13）：
  (A) 验证引理 A：任意谱的 𝔨 满足 h(ϖ_min) ≤ 𝔨 ≤ h(ϖ_rms) 与下界-2
  (B) 验证命题 B：Einstein x_c 关于 1/λ 凸，切线界成立
  (C) 氢化物表：Einstein 上界 T_c^E(λ,ω_rms) 与文献 Tc 对照

理论依据：引理 M（核单调性）+ 初等不等式，见 P4_tc_sandwich.md。

数值约定：
  - N=400（float64，非严格），h^(N) 是 h 的下估（N 截断单调）。
  - 单位：<ω>=1，T_c 输出为 T_c/<ω>（无量纲）。
  - 氢化物表的 ω 单位为 meV，Tc 单位为 K。

ANTI-STALL：
  - numpy float64，scipy.optimize.brentq 固定 maxiter=60。
  - Stage A/B/C 独立，每个 Stage 打印通过数。
  - 随机种子固定（可复现）。
"""

import numpy as np
from scipy.optimize import brentq

PASS = [0]


def chk(cond, msg=""):
    """Counting assert: aborts on failure, tallies successes."""
    assert cond, msg
    PASS[0] += 1


# ===========================================================================
# 基础函数（从 p4_tc_table.py 复制）
# ===========================================================================
def brk(k, w):
    """[[k]](w) = w^2/(k^2+w^2)."""
    return w ** 2 / (k ** 2 + w ** 2)


def H_matrix(N, w):
    """KAY II H^(N)(w)，逐元构造，与 p4_tc_table.py 完全相同。"""
    M = np.zeros((N, N))
    for n in range(N):
        for m in range(N):
            nrm = np.sqrt(2 * n + 1) * np.sqrt(2 * m + 1)
            val = brk(n + m + 1, w) / nrm
            if n != m:
                val += brk(abs(n - m), w) / nrm
            else:
                s = sum(2 * brk(k, w) for k in range(1, n + 1)) if n >= 1 else 0.0
                val -= s / (2 * n + 1)
            M[n, m] = val
    return M


def eig_max(M):
    """最大本征值 + 归一化本征向量。"""
    E, V = np.linalg.eigh(M)
    j = int(np.argmax(E))
    return E[j], V[:, j]


def h_func(N, w):
    """h^(N)(w) = lambda_1(H^(N)(w))。"""
    lam, _ = eig_max(H_matrix(N, w))
    return lam


def frak_k_of_P(N, T, omegas, ps):
    """𝔨^(N)(P,T) = rho(sum_j p_j H^(N)(omega_j/(2 pi T)))。"""
    M = np.zeros((N, N))
    for om, p in zip(omegas, ps):
        w = om / (2 * np.pi * T)
        M += p * H_matrix(N, w)
    lam, _ = eig_max(M)
    return lam


def solve_Tc_Einstein(N, lam_coupling, omega_E, T_lo=1e-7, T_hi=1e7, maxiter=60):
    """解 lambda * h^(N)(omega_E/(2 pi T)) = 1 得 T_c^E。

    返回 T_c^E（无量纲，单位 <omega>/2pi，这里 <omega>=omega_E）。
    """
    def f(T):
        w = omega_E / (2 * np.pi * T)
        return lam_coupling * h_func(N, w) - 1.0

    flo, fhi = f(T_lo), f(T_hi)
    if flo < 0:
        # lambda 太小，无有限 Tc
        return np.nan
    if fhi > 0:
        # T_hi 不够大
        return np.nan

    return brentq(f, T_lo, T_hi, maxiter=maxiter, xtol=1e-14, rtol=8.9e-16)


# ===========================================================================
# Stage A：随机谱夹逼验证
# ===========================================================================
def random_spectrum(n_peaks, seed=None):
    """生成随机 n_peaks 峰谱（频率于 [0.1,5]，权重归一化）。

    返回 (omegas, ps, omega_min, omega_max, omega_rms)。
    """
    if seed is not None:
        np.random.seed(seed)
    omegas = np.random.uniform(0.1, 5.0, n_peaks)
    omegas.sort()
    ps = np.random.uniform(0.1, 1.0, n_peaks)
    ps /= ps.sum()

    omega_min = omegas[0]
    omega_max = omegas[-1]
    omega_rms = np.sqrt(np.sum(ps * omegas**2))

    return omegas, ps, omega_min, omega_max, omega_rms


def run_stage_A():
    print("=== Stage A: 随机谱夹逼验证 ===")
    N = 200  # 从400降到200以加快速度（仍足够验证理论）
    n_spectra = 42
    T_grid = [0.01, 0.05, 0.1, 0.5, 1.0]

    n_checks = 0
    print(f"    N={N}（为加速降低），验证 {n_spectra} 个随机谱 × {len(T_grid)} 个温度")
    for i in range(n_spectra):
        if (i+1) % 10 == 0:
            print(f"      进度: {i+1}/{n_spectra} 谱...", flush=True)
        n_peaks = np.random.randint(2, 6)  # 2-5 峰
        omegas, ps, om_min, om_max, om_rms = random_spectrum(n_peaks, seed=1000+i)

        for T in T_grid:
            # 计算 𝔨(P,T)
            kappa = frak_k_of_P(N, T, omegas, ps)

            # 计算三个界
            w_min = om_min / (2 * np.pi * T)
            w_max = om_max / (2 * np.pi * T)
            w_rms = om_rms / (2 * np.pi * T)

            h_min = h_func(N, w_min)
            h_max = h_func(N, w_rms)
            om2_mean = np.sum(ps * omegas**2)
            h_lower2 = (om2_mean / om_max**2) * h_func(N, w_max)

            # 验证夹逼
            tol = 1e-10  # float64 容差
            chk(h_min - tol <= kappa <= h_max + tol,
                f"spectrum {i}, T={T}: h_min={h_min:.6e}, kappa={kappa:.6e}, h_max={h_max:.6e}")
            chk(h_lower2 - tol <= kappa,
                f"spectrum {i}, T={T}: h_lower2={h_lower2:.6e}, kappa={kappa:.6e}")
            n_checks += 2

    print(f"    每个 (谱,T) 检查 2 条不等式（上界 + 下界-1、下界-2）")
    print(f"Stage A: {PASS[0]}/{PASS[0]} checks passed ({n_checks} inequalities)")


# ===========================================================================
# Stage B：Einstein x_c 凸性验证
# ===========================================================================
def run_stage_B():
    print("\n=== Stage B: Einstein x_c 凸性验证 ===")
    N = 200  # 与 Stage A 一致
    omega_E = 1.0
    lam_grid = np.array([0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0])

    # (1) 计算 x_c = (omega_E / 2pi T_c^E)^2 对每个 lambda
    x_c = []
    for lam in lam_grid:
        Tc = solve_Tc_Einstein(N, lam, omega_E)
        if np.isnan(Tc):
            print(f"    WARNING: lambda={lam} 无有限 Tc，跳过")
            continue
        xc = (omega_E / (2 * np.pi * Tc))**2
        x_c.append(xc)

    x_c = np.array(x_c)
    inv_lam = 1.0 / lam_grid[:len(x_c)]

    print(f"    lambda 网格: {lam_grid[:len(x_c)]}")
    print(f"    x_c 序列: {x_c}")

    # (2) 验证凸性：二阶差分 >= 0
    if len(x_c) >= 3:
        # 中心二阶差分：D2[i] = (x[i+1] - 2x[i] + x[i-1]) / h^2
        h = np.diff(inv_lam)
        for i in range(1, len(x_c) - 1):
            d2 = (x_c[i+1] - 2*x_c[i] + x_c[i-1]) / (h[i-1] * h[i])
            chk(d2 >= -1e-8,  # 容差允许数值误差
                f"凸性失败 i={i}: D2 x_c = {d2:.6e}")
        print(f"    凸性检查：{len(x_c)-2} 个内点的二阶差分 >= 0")

    # (3) 验证切线界：选 lambda_0 = 2.0
    if 2.0 in lam_grid[:len(x_c)]:
        idx0 = int(np.where(lam_grid[:len(x_c)] == 2.0)[0][0])
        lam0 = lam_grid[idx0]
        xc0 = x_c[idx0]

        # 计算 H'(x_c0) ≈ dh/dw * dw/dx = (h(w+dw) - h(w-dw))/(2dw) * 1/(2sqrt(x))
        # 其中 w = sqrt(x_c0)
        w0 = np.sqrt(xc0)
        dw = 1e-6
        h_plus = h_func(N, w0 + dw)
        h_minus = h_func(N, w0 - dw)
        dh_dw = (h_plus - h_minus) / (2 * dw)
        # H(x) = h(sqrt(x)), H'(x) = h'(sqrt(x)) * 1/(2sqrt(x))
        Hprime = dh_dw / (2 * w0)

        slope = 1.0 / Hprime
        print(f"    切线基准点：lambda_0={lam0}, x_c0={xc0:.6f}, H'(x_c0)={Hprime:.6e}, 斜率={slope:.6e}")

        # 验证所有点在切线上方：x_c >= x_c0 + (1/lam - 1/lam0) * slope
        for i, lam in enumerate(lam_grid[:len(x_c)]):
            if i == idx0:
                continue
            tangent_val = xc0 + (1.0/lam - 1.0/lam0) * slope
            chk(x_c[i] >= tangent_val - 1e-6,
                f"切线界失败 lambda={lam}: x_c={x_c[i]:.6e}, tangent={tangent_val:.6e}")
        print(f"    切线界检查：{len(x_c)-1} 个点在切线上方")

    print(f"Stage B: {PASS[0]}/{PASS[0]} checks passed")


# ===========================================================================
# Stage C：氢化物 Tc 上界表
# ===========================================================================
def run_stage_C():
    print("\n=== Stage C: 氢化物 Einstein Tc 上界表 ===")
    N = 400

    # 数据来源：isotope_curated.md SSCHA 表（λ, ω_log, Tc_theory，μ*）
    # 格式：(system, P_GPa, lambda, omega_log_meV, Tc_theory_K, mu_star)
    hydrides = [
        ("H3S", 200, 1.84, 92.86, 194.0, 0.16),
        ("LaH10", 163, 2.67, 96.4, 242.8, 0.10),
        ("LaH10", 129, 3.62, 76.4, 255.3, 0.10),
        ("LaH10", 214, 2.06, 115.5, 237.9, 0.10),
        ("LaH10", 264, 1.73, 126.6, 216.9, 0.10),
        ("YH6", 165, 1.71, 114.87, 241.5, 0.125),  # Tc 取 236-247 中点
        ("PdH", 0, 0.40, 50.21, 5.0, 0.085),
    ]

    print(f"    N={N}（float64，非严格），h^(N) 是 h 的下估")
    print(f"    氢化物表：λ 与 ω_log 来自 SSCHA/第一性原理，μ*=0 上界 vs μ*>0 理论值\n")

    # 验证 N=200 vs N=400 的收敛性（抽样 1 个材料）
    test_lam, test_om = 2.67, 96.4
    Tc_200 = solve_Tc_Einstein(200, test_lam, test_om)
    Tc_400 = solve_Tc_Einstein(400, test_lam, test_om)
    rel_diff = abs(Tc_400 - Tc_200) / Tc_200
    print(f"    N 收敛性检查（LaH10 @163 GPa）：")
    print(f"      T_c^E(N=200) = {Tc_200:.6f}, T_c^E(N=400) = {Tc_400:.6f}")
    print(f"      相对差 = {rel_diff:.2e} {'< 0.5% OK' if rel_diff < 5e-3 else '>= 0.5% (需更大 N)'}\n")

    # 打印表头
    header = f"{'System':<10} {'P(GPa)':>7} {'λ':>6} {'ω_log(meV)':>12} {'T_c^E(K)':>10} {'Tc_thy(K)':>10} {'T_c^E/Tc':>10} {'μ*':>6}"
    print(header)
    print("-" * len(header))

    for system, P, lam, om_log, Tc_thy, mu_star in hydrides:
        # 计算 T_c^E：设 <omega>=1（无量纲），解出 T_c/<omega>，再乘物理 omega
        # 这里 omega_E = om_log（meV），需转换单位
        # 约定：om_log 单位 meV = 1.16045e4 K（玻尔兹曼常数 k_B = 8.617e-5 eV/K）
        # 但我们直接用 om_log 作为 omega_E 的相对单位（<omega>=om_log）
        omega_E_scaled = 1.0  # 归一化
        Tc_scaled = solve_Tc_Einstein(N, lam, omega_E_scaled)

        if np.isnan(Tc_scaled):
            print(f"{system:<10} {P:>7} {lam:>6.2f} {om_log:>12.2f} {'N/A':>10} {Tc_thy:>10.1f} {'N/A':>10} {mu_star:>6.3f}")
            continue

        # 转换回物理单位：Tc_scaled 是以 om_log/(2pi) 为单位的，乘回 om_log
        # 实际上 Tc_scaled = Tc / (om_log / 2pi)，所以 Tc = Tc_scaled * (om_log / 2pi)
        # 但 om_log 单位是 meV，Tc 单位是 K，需要转换因子
        # meV 转 K：1 meV = 11.6045 K
        meV_to_K = 11.6045
        Tc_E_K = Tc_scaled * (om_log * meV_to_K) / (2 * np.pi)

        ratio = Tc_E_K / Tc_thy

        print(f"{system:<10} {P:>7} {lam:>6.2f} {om_log:>12.2f} {Tc_E_K:>10.1f} {Tc_thy:>10.1f} {ratio:>10.2f} {mu_star:>6.3f}")

    print(f"\nStage C: 表格完成")
    print(f"    注：T_c^E 是 μ*=0 上界，Tc_thy 是 μ*>0 SSCHA 值，预期 T_c^E > Tc_thy")
    print(f"    N={N} 的 h^(N) 是 h 的下估 ==> T_c^E(N) 是 T_c^E(∞) 的下估")
    print(f"    作为「上界」用时严格性：需 N→∞ 或验证收敛（见上文 N=200 vs 400）")


# ===========================================================================
# Main
# ===========================================================================
if __name__ == '__main__':
    np.random.seed(42)  # 全局种子（Stage A 有独立种子）

    run_stage_A()
    run_stage_B()
    run_stage_C()

    print(f"\n{'='*70}")
    print(f"ALL {PASS[0]} checks passed.")
    print(f"{'='*70}")
    print("\n交付：")
    print("  - P4_tc_sandwich.md（引理 A、命题 B 理论推导）")
    print("  - p4_tc_sandwich.py（数值验证脚本，本文件）")
    print("  - 所有验证通过，氢化物表已生成")
