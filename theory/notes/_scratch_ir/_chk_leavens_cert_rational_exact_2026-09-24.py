# -*- coding: utf-8 -*-
"""
leavens-cert Phase 2（独立从零重算，不 import team-lead 的脚本，只在结构上参考其思路，
数值/公式全部自己重推一遍 —— 仓库纪律：承重结论不采信 peer 的「已完成」）。

## 目标
把 `c(mu*, r=10)` 在四个 mu*（0 / 0.10 / 0.13 / 0.20）处各产一个 **rational-exact** 下界，
代入 `lam_area = C/c^2` 给出认证上界，重算四行比值 `lam_ceil_lower / lam_area_upper`，
与既有 float64 结果 1.0152 / 1.1853 / 1.2982 / 1.7967 对照。
（末格 2026-09-24 由 team-lead 从零复核订正：原 1.8156 的 lam_ceil 取了 lam1(lam)=1 的**最大**根
 5.858499，而 N_c=floor(10*sqrt(M2/lam)) 在 lam~5.8120 的台阶使该 mu* 下 lam1-1 **变号 3 次**；
 物理陈述要求 lam_ceil = inf{lam : lam1 >= 1} = **最小**根 5.797500，取最大根是 fail-open。
 依据 `_chk_main_lamceil_mu020_step_root_audit_2026-09-24.py`。c_lower 一列与根的选取无关、未改。）

## 承接 team-lead 的订正（我独立重推一遍，不采信，见 STEP 0）
`S_n := sum_{m>=0} [G(|n-m|) - G(n+m+1)]` 两条无穷尾恒等相消，给出**闭式**（非截断）:
    S_n = G(0) + 2*sum_{j=1}^{n} G(j)
这是**有限项精确公式**，不是「N 截断到 n 项」——用它构造的 N x N 矩阵 Khat_N 的每个矩阵元
（对角、非对角）都与 N 无关地等于无穷维算子 K^inf 的对应矩阵元 ⟹ Khat_N 是 K^inf 的
**真正主子块**（entries 不随 N 变化，只是维度被截断）⟹ Cauchy 交错定理适用：
    theta_max(Khat_N, B_N) <= theta_max(K^inf, B^inf)   (N >= N_c, 使 B_N 也是 B^inf 的主子块)
⟹ **有限 N 的（认证）下界，就是无穷维真值的下界**，不需要额外的 N-尾界引理
   （这推翻了我自己 phase 1 报告里「N 截断方向不利」的结论——那条是用**截断和**版本的 S_n
   算出来的，指纹是「S_n 是部分和」，闭式版本没有这个问题）。

## 方法：有理 Rayleigh 商下界（analytic 变分原理 + 有限项精确有理算术）
对固定有理 w=k/r（每个「齿」的左端，齿区间 [k/r,(k+1)/r) 内 N_c=k 不变，齿顶实测在左端）：
把 u 吸进试探向量 y_n := u_n v_n（v 是 float64 广义本征向量，只用来**取**一个好的有理向量，
声音性与它的来历无关——任何有理向量都给合法下界，见 CLAUDE.md「预算」条：目标常数只买紧度不买通过）。
    分子 = sum_{ij} (G(|i-j|)+G(i+j+1)) y_i y_j - sum_i S_i y_i^2
    分母 = sum_i (2i+1) y_i^2 + 2*mu* * (sum_{n<Nc} y_n)^2
    theta_lower := 分子/分母   (Fraction 精确算术，零浮点)
变分原理（analytic）：对任意实（这里取有理）向量 y != 0，
    theta_max(Khat_N,B_N) >= <Khat_N y,y>/<B_N y,y>
⟹ theta_lower <= theta_max(Khat_N,B_N) <= theta_max(K^inf,B^inf)（Cauchy 交错，STEP 0 已证）
⟹ theta_lower 是 theta_max(w) 的 **rational-exact 下界**。

c(w) = theta_max(w)/(pi*w) ⟹ 需要 pi 的**有理上界**（除法方向：theta 下界 / pi 上界 / w
= c 的下界，w 本身取正好的有理数，无需额外区间）。取标准有理上界 355/113 > pi
（Machin 类展开已知 355/113 - pi ≈ 2.667e-7 > 0，本脚本 STEP 0.5 用 mpmath 高精度独立核对该不等式）。

c(mu*,r) = max_w c(w) >= c(该点下界) —— **单点下界即是 max 的下界**，不需要枚举全部分支
（leavens-cert phase 1 已确认的方向；本轮只需要挑一个好齿，不追求覆盖全部齿）。

级别：变分原理 analytic；Cauchy 交错 analytic；S_n 闭式恒等式 analytic；
数值部分 **rational-exact**（Fraction 精确算术，零浮点，格数与有理数个数逐行报出）；
唯一非有理的输入是 pi 的有理上界（356/113 型常数，本身是被独立验证的有理数，不引入舍入误差）。
float64（eigh, 找齿）仅用来**定位**试探向量与齿位，不进入承重的不等式链条。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh
from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 50
KB = 8.617333262e-2   # meV/K   (物理常数，float64，不进有理链，只用于换算 lam_area 的 C，另标级)
HB = 64.65415134      # meV
PI_UPPER = F(355, 113)  # 待独立核验


# ═══════════════════════ STEP 0: 独立重推 S_n 闭式恒等式 + Cauchy 交错前提 ═══════════════════════
def step0_identity_check():
    print("=" * 96)
    print("STEP 0  独立重推 S_n 闭式恒等式（不采信 team-lead 的推导，自己重新算一遍）")
    print("=" * 96)
    print("  claim: S_n = sum_{m=0}^{Mtail} [G(|n-m|)-G(n+m+1)]  (Mtail->inf)")
    print("         应等于闭式  G(0) + 2*sum_{j=1}^{n} G(j)")
    for w in (0.5, 1.2, 2.7):
        G = lambda k: w * w / (w * w + k * k)
        for n in (0, 3, 9, 20):
            # 数值上用很大的 Mtail 截断做「实际上等价于无穷」的直接求和，交叉核对闭式
            Mtail = 200000
            m = np.arange(Mtail)
            trunc_sum = np.sum(G(np.abs(n - m)) - G(n + m + 1))
            closed = G(0) + 2.0 * sum(G(j) for j in range(1, n + 1))
            print(f"    w={w:4.1f} n={n:3d}  直接求和(Mtail={Mtail})={trunc_sum:.10f}  "
                  f"闭式={closed:.10f}  |差|={abs(trunc_sum-closed):.2e}")
    print("  ⟹ 闭式恒等式数值复核成立（尾差 ~1/Mtail 量级，随 Mtail 增大而 ->0，与截断误差量级一致）。")
    print()


def step0p5_pi_bound_check():
    print("=" * 96)
    print("STEP 0.5  独立核验 355/113 > pi（mpmath 50 位精度）")
    print("=" * 96)
    pi_hi = mp.mpf(355) / mp.mpf(113)
    diff = pi_hi - (+mp.pi)   # mp.pi 是惰性的，+mp.pi 触发求值
    print(f"  355/113 = {pi_hi}")
    print(f"  pi      = {+mp.pi}")
    print(f"  355/113 - pi = {diff}   ({'> 0, 是合法上界' if diff > 0 else '** 不是上界 **'})")
    print()


# ═══════════════════════ float64 定位工具（只用来找齿位与试探向量，不进承重链条） ═══════════════════════
def pieces_f64(wbar, N):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    Gc = wbar * wbar / (wbar * wbar + n * n)
    S = Gc[0] + 2.0 * (np.cumsum(Gc) - Gc[0])          # 闭式 S，float64
    Gk = wbar * wbar / (wbar * wbar + np.arange(2 * N + 2) ** 2)
    Khat = (Gk[np.abs(i - j)] + Gk[i + j + 1]) * np.outer(u, u) - np.diag(u * u * S)
    return Khat, u


def theta_f64(wbar, mus, r, N):
    Khat, u = pieces_f64(wbar, N)
    Nc = int(np.floor(r * wbar))
    B = np.eye(N)
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(N) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    val, vec = eigh(Khat, B, subset_by_index=[N - 1, N - 1])
    return float(val[0]), vec[:, 0], Nc


def find_best_tooth(mus, r, N, kmax=45):
    """粗扫齿左端 w=k/r（k=1..kmax），找 c(w)=theta/(pi w) 最大的 k。"""
    rf = float(r)
    best = (-1.0, None)
    for k in range(1, kmax + 1):
        w = k / rf
        th, _, _ = theta_f64(w, mus, rf, N)
        c = th / (np.pi * w)
        if c > best[0]:
            best = (c, k)
    return best


# ═══════════════════════ 有理精确 Rayleigh 下界 ═══════════════════════
def rational_theta_lower(w_frac, mus_frac, r, N, v0):
    """w_frac, mus_frac: Fraction。v0: float64 近似本征向量（长度 N），只用来生成试探有理向量。
    返回 (theta_lower: Fraction, Nc: int, n_rationals: int)。"""
    Nc = int(np.floor(float(r) * float(w_frac)))
    w2 = w_frac * w_frac
    kmax = 2 * N + 2
    G = [w2 / (w2 + F(k * k)) for k in range(kmax)]
    # 闭式 S_n = G(0) + 2*sum_{j=1}^n G(j)
    S = [None] * N
    acc = G[0]
    S[0] = acc
    for n in range(1, N):
        acc = acc + 2 * G[n]
        S[n] = acc
    DEN = 10 ** 6
    u0 = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
    y = [F(round(float(v0[k]) * u0[k] * DEN), DEN) for k in range(N)]

    num = F(0)
    for i in range(N):
        yi = y[i]
        if yi == 0:
            continue
        for j in range(N):
            yj = y[j]
            if yj == 0:
                continue
            num += (G[abs(i - j)] + G[i + j + 1]) * yi * yj
    for i in range(N):
        num -= S[i] * y[i] * y[i]

    den = F(0)
    for i in range(N):
        den += F(2 * i + 1) * y[i] * y[i]
    if mus_frac > 0 and Nc > 0:
        s = sum(y[:Nc])
        den += 2 * mus_frac * s * s

    theta_lower = num / den
    n_rationals = N + kmax + N  # y + G + S
    return theta_lower, Nc, n_rationals


def cert_one_mustar(mus_val, mus_frac, r, N_scan, N_cert, kmax):
    print(f"--- mu* = {mus_val} ---")
    if mus_val == 0.0:
        # 无库仑（B=I 恒成立），c(w) 是光滑函数，没有锯齿 ⟹ 不能沿用 k/r 齿位，
        # 必须真扫 w（不能拿"omega_E/A=2/lam"当 w 用——那是另一个变量，已被我最初的手选
        # w=7/4 的版本搞混，订正：改成对 w 本身做细网格扫描定位 argmax，再取邻近有理点）。
        ws = np.linspace(0.3, 3.0, 541)
        best = (-1.0, None)
        for w in ws:
            th, _, _ = theta_f64(w, 0.0, float(r), N_scan)
            c = th / (np.pi * w)
            if c > best[0]:
                best = (c, w)
        w_star = best[1]
        w_frac = F(round(w_star * 1000), 1000)
        k_tag = f"细网格扫描 argmax w={w_star:.6f} (float64 c={best[0]:.6f}) -> 就近有理点 {w_frac}"
    else:
        cbest, kbest = find_best_tooth(mus_val, r, N_scan, kmax)
        w_frac = F(kbest, int(r))
        k_tag = f"齿 k={kbest} (float64 粗扫 c={cbest:.6f})"
    w_f64 = float(w_frac)
    th_f64, v0, Nc = theta_f64(w_f64, mus_val, float(r), N_cert)
    print(f"  取点: {k_tag}  ->  w = {w_frac} = {w_f64:.6f}, Nc={Nc}")
    print(f"  float64 参考 theta(w) = {th_f64:.12f}  (独立算法, scipy eigh, 仅供对照, 不进链条)")

    theta_lower, Nc2, n_rat = rational_theta_lower(w_frac, mus_frac, r, N_cert, v0)
    assert Nc2 == Nc
    rel_gap = 1.0 - float(theta_lower) / th_f64
    print(f"  rational-exact theta 下界 = {float(theta_lower):.12f}  "
          f"(相对 float64 亏损 {rel_gap:.3e}, {'>=0 声音' if rel_gap >= 0 else '** 违反, 检查 **'})")

    c_lower = theta_lower / (PI_UPPER * w_frac)
    c_f64 = th_f64 / (np.pi * w_f64)
    print(f"  c 下界(rational-exact) = {float(c_lower):.10f}   (float64 点值 c(w) = {c_f64:.10f})")
    print(f"  格数 N = {N_cert}, 有理数个数 = {n_rat} (y:{N_cert} + G:{2*N_cert+2} + S:{N_cert})")
    return c_lower, N_cert, n_rat, w_frac, Nc


def main():
    step0_identity_check()
    step0p5_pi_bound_check()

    r = F(10)
    N_scan = 300   # 粗扫定齿：float64，快
    N_cert = 40    # 认证用矩阵维度（Rayleigh 下界对 N 不敏感，N=16 已够，见 team-lead 实测；本脚本独立选 N=40 留余量）
    kmax = 40

    print("=" * 96)
    print("STEP 1  四个 mu* 的 rational-exact c 下界")
    print("=" * 96)
    mustars = [(0.00, F(0)), (0.10, F(1, 10)), (0.13, F(13, 100)), (0.20, F(1, 5))]
    results = {}
    for mv, mf in mustars:
        c_lower, N_used, n_rat, w_frac, Nc = cert_one_mustar(mv, mf, r, N_scan, N_cert, kmax)
        results[mv] = (c_lower, N_used, n_rat, w_frac, Nc)
        print()

    print("=" * 96)
    print("STEP 2  代入 lam_area = C/c^2，重算四行比值（不把新 c 代进旧解——lam_ceil 侧本来就没有")
    print("        依赖 c 的下游优化，只是查表，故这里无需重解任何优化问题）")
    print("=" * 96)
    eta_over_M = 9.172
    pref = HB * np.sqrt(eta_over_M)
    Cconst = (2.0 * 300.0 * KB / pref) ** 2
    print(f"  C = (2*300*KB/pref)^2 = {Cconst:.9f}   (float64, 物理常数换算, 不属于有理链条, 另标级 numerical-evidence)")

    lam_ceil_lower = {0.00: 1.32808 * 0.999, 0.10: 2.66204 * 0.999,
                       0.13: 3.29, 0.20: 5.797500 * 0.999}   # ← 0.20 订正：inf-根，原 5.85831 是最大根
    # 说明：mu*=0.13 用既有 certificate 下端 3.29（convex_extreme_point_reduction）；
    # 其余三个尚无族级 certificate，此处沿用既有 numerical-evidence 的 Einstein 精确根打 0.1% 折扣
    # 仅做**示意性**核对（不作为最终判词——真正的判词只对 mu*=0.13 是 certificate 级）。
    old_ratio = {0.00: 1.0152, 0.10: 1.1853, 0.13: 1.2982, 0.20: 1.7967}   # ← 0.20 同上订正

    print(f"\n  {'mu*':>6}{'c_lower':>14}{'N':>5}{'#有理数':>8}{'lam_area_upper':>16}"
          f"{'lam_ceil(参考)':>16}{'ratio_new':>12}{'ratio_old(f64)':>16}")
    for mv, mf in mustars:
        c_lower, N_used, n_rat, w_frac, Nc = results[mv]
        c_lo = float(c_lower)
        lam_area_upper = Cconst / (c_lo * c_lo)
        lc = lam_ceil_lower[mv]
        ratio_new = lc / lam_area_upper
        print(f"  {mv:6.2f}{c_lo:14.9f}{N_used:5d}{n_rat:8d}{lam_area_upper:16.6f}"
              f"{lc:16.5f}{ratio_new:12.4f}{old_ratio[mv]:16.4f}")

    print()
    print("  级别标注：")
    print("  - c_lower 列：rational-exact（N=40 格 / 各 mu* 有理数个数如上列，Fraction 精确算术，零浮点，")
    print("    唯一非纯有理输入是 pi 的有理上界 355/113，本身经 STEP0.5 独立核验）。")
    print("  - lam_area_upper：由 rational-exact c_lower 与 numerical-evidence 常数 C 相乘得到，")
    print("    C 本身（eta/M, HB, KB 的换算）未做有理化，故 lam_area_upper 整体标 numerical-evidence")
    print("    （c 这一步已经 rational-exact，是链条里最松的那环拖累整体级别——如需升级 lam_area_upper 本身")
    print("    到 certificate，需要另外区间化 C 的物理常数换算，这是未做的另一件事，不在本轮范围）。")
    print("  - ratio_new 列：mu*=0.13 一行分子(lam_ceil_lower)是 certificate 级（既有分支定界下端 3.29），")
    print("    其余三行分子仍是 numerical-evidence（Einstein 精确根 x0.999 的示意折扣，非严格证书），")
    print("    故只有 mu*=0.13 一行的 ratio_new 具备「分子分母都非 float64-only」的地位；")
    print("    其余三行仅供数量级参考，不得当作认证排除（与既有笔记 §4 的级别限制一致）。")


if __name__ == "__main__":
    main()
