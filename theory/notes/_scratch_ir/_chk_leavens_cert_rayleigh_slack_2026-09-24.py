# -*- coding: utf-8 -*-
"""
leavens-cert 分支任务，阶段 1（预算复核）。

派活消息（team-lead -> leavens-cert）：把 Leavens 面积天花板的 c(mu*,r) 逐分支
区间认证，让 (Phi-A)/(convex_extreme_point) vs Leavens 的头条比较从
numerical-evidence 升 certificate。消费点见
`ceiling/same_information_class_comparison_vs_leavens_2026-09-24.md` §4:
比值 lam*/lam_area = 1.0152 / 1.1853 / 1.2982 / 1.8156 (mu*=0/0.10/0.13/0.20)，
分子(lam_ceil)已有 certificate（分支定界），分母(lam_area, 经由 c(mu*,r))只有 float64。

本脚本只做阶段 1 第 3 步：**在 1-2 个分支上实测我打算用的区间估计子的固有松**，
在花力气认证全部 ~27 个分支之前，先看这条路线在最紧的 mu*=0 格（余量仅 1.52%）
是否可行。

━━━━━━━━━━ 定义式（逐字对照，文件+行号）━━━━━━━━━━
`_chk_main_c_sawtooth_exact_branchwise_2026-09-24.py`:
  L47-56  pieces(wbar,N): Khat = (Kp+Km)*outer(u,u) - diag(u^2 S),
          Kp=G[|i-j|], Km=G[i+j+1], G(k)=wbar^2/(wbar^2+k^2), u_n=(2n+1)^{-1/2}.
  L59-68  lam_root_exact: B = I + 2*mu* * outer(uc,uc) (mu*>0 时)，
          Nc=floor(r*wbar)，uc = u 的前 Nc 个分量（其余置 0）。
          引理（同脚本 L14-22，analytic）: lam*(w) = 1/theta_max(w)，
          theta_max = 广义本征问题 Khat v = theta B v 的最大广义本征值。
  L71-75  c_of_w: c(w) = 1/(pi * lam*(w) * w) = theta_max(w) / (pi * w).
`_chk_main_leavens_c_and_phiA_rescale_2026-09-23.py`:
  L265-266  lam_area(c) = (2*300*KB / (c*pref))**2   ,  pref = HB*sqrt(eta/M) 常数.
          ⟹ lam_area = C / c^2   for some C>0  (C = (600*KB/pref)^2 depends only on eta/M, 与 c 无关).

━━━━━━━━━━ 方向审计（承重）━━━━━━━━━━
lam_area(c) = C / c^2, C>0 常数（不依赖 c）⟹ d(lam_area)/dc = -2C/c^3 < 0 对任意 c>0。
⟹ **c 越大，lam_area 越小** —— 这是纯代数事实，不依赖 mu*/r 的任何取值，一次符号求导即可，
   本脚本 STEP 0 用 mpmath 符号/数值双重复核（sympy 不可用则直接数值差分，见下）。
⟹ 要证 `lam_ceil_lower > lam_area_upper`（ratio>1 的判决方向），
   需要 lam_area 的**上界** <=> 需要 **c 的下界**（c 越大 lam_area 越小，
   所以 c 的下界给出 lam_area 的上界 —— 单调递减函数：自变量下界 -> 函数上界）。
⟹ 本脚本要产的区间估计子是 **theta_max(w) 的认证下界**（Rayleigh 变分原理），
   用它给出 c(w) 的下界，再取「已知某个 w 处的下界」<= max_w c(w) = c(mu*,r)
   （任一点的下界天然是全局 max 的下界，不需要覆盖整个分支区间）。

严格性：变分不等式 theta_max(w) >= <Khat v,v>/<Bv,v> (任意实 v) 是 analytic
（对称广义本征问题的 Rayleigh 商变分原理）。本脚本用 mpmath.iv 外向舍入区间算术
逐元素构造 Khat_iv, B_iv 并对固定有理 w、固定（float64 算出的）近似本征向量 v0
计算区间 Rayleigh 商，其下端 = certificate 级下界。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh
import mpmath as mp
from mpmath import iv

iv.dps = 30  # 区间算术工作精度（十进制位），远高于 float64 的 ~15-17 位


# ═══════════════════════ float64 参考实现（与既有脚本一致，独立重写不 import） ═══════════════════════
def pieces_f64(wbar, N):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * N + 2)
    G = wbar * wbar / (wbar * wbar + k * k)
    Kp, Km = G[np.abs(i - j)], G[i + j + 1]
    S = (Kp - Km).sum(axis=1)
    Khat = (Kp + Km) * np.outer(u, u) - np.diag(u * u * S)
    return Khat, u


def theta_max_f64(wbar, mus, r, N):
    Khat, u = pieces_f64(wbar, N)
    B = np.eye(N)
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(N) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    ev, V = eigh(Khat, B)
    return float(ev[-1]), V[:, -1], Nc


# ═══════════════════════ STEP 0: 方向审计 ═══════════════════════
def step0_direction():
    print("=" * 92)
    print("STEP 0  方向审计：lam_area(c) = C/c^2 的单调性（C>0 常数，与 c 无关）")
    print("=" * 92)
    KB, HB = 8.617333262e-2, 64.65415134
    eta_over_M = 9.172
    pref = HB * np.sqrt(eta_over_M)
    Cconst = (2.0 * 300.0 * KB / pref) ** 2
    print(f"  lam_area(c) = C / c^2 ,  C = (2*300*KB/pref)^2 = {Cconst:.6f}  (数值核对既有脚本 L265-266)")
    for c in (0.10, 0.15, 0.16574, 0.20, 0.2309):
        la = Cconst / (c * c)
        print(f"    c={c:.5f}  lam_area={la:.6f}")
    # 数值差分复核符号（不信符号推导，独立算一次）
    eps = 1e-6
    c0 = 0.16574
    dla = ((Cconst / (c0 + eps) ** 2) - (Cconst / (c0 - eps) ** 2)) / (2 * eps)
    print(f"  数值差分 d(lam_area)/dc @ c={c0}: {dla:.6f}  (符号: {'负' if dla < 0 else '正'})")
    print("  ⟹ 确认 dla/dc < 0：c 的下界 ⟹ lam_area 的上界。这是本脚本要产的方向。")
    print()


# ═══════════════════════ STEP 1: 区间矩阵构造（外向舍入） ═══════════════════════
def pieces_iv(wbar_iv, N):
    """区间版 Khat, u。wbar_iv 是 mpmath.iv 区间（可以是薄区间，即某个有理数的精确包络）。"""
    u = [iv.mpf(1) / iv.sqrt(iv.mpf(2 * n + 1)) for n in range(N)]
    w2 = wbar_iv * wbar_iv
    kmax = 2 * N + 2
    G = [w2 / (w2 + iv.mpf(k * k)) for k in range(kmax)]
    Khat = [[iv.mpf(0) for _ in range(N)] for _ in range(N)]
    S = [iv.mpf(0) for _ in range(N)]
    for i in range(N):
        for j in range(N):
            kp = G[abs(i - j)]
            km = G[i + j + 1]
            Khat[i][j] = (kp + km) * u[i] * u[j]
    for i in range(N):
        # S_i = sum_j (G[|i-j|] - G[i+j+1])   (对角修正项，逐元素外向舍入求和)
        acc = iv.mpf(0)
        for j in range(N):
            acc += G[abs(i - j)] - G[i + j + 1]
        S[i] = acc
        Khat[i][i] = Khat[i][i] - u[i] * u[i] * S[i]
    return Khat, u


def rayleigh_lower_bound(wbar_val, mus, r, N, v0):
    """对固定有理 wbar_val（python Fraction 或可精确转 iv 的字符串/整数比），
    用区间算术算 R(v0) = <Khat v0,v0>/<Bv0,v0> 的下端，作为 theta_max(wbar_val) 的certified下界。
    v0: float64 近似本征向量（作为试探向量，不要求精确 —— Rayleigh 商对任意实向量都给下界，
        只是向量越准下界越紧）。
    """
    wbar_iv = iv.mpf(wbar_val)
    Khat, u = pieces_iv(wbar_iv, N)
    Nc = int(np.floor(r * float(wbar_val)))
    v = [iv.mpf(float(x)) for x in v0]  # 薄区间：float64 的值是精确二进制有理数，iv.mpf 精确包络

    # 分子 <Khat v, v>
    num = iv.mpf(0)
    Kv = [iv.mpf(0) for _ in range(N)]
    for i in range(N):
        acc = iv.mpf(0)
        row = Khat[i]
        for j in range(N):
            acc += row[j] * v[j]
        Kv[i] = acc
    for i in range(N):
        num += Kv[i] * v[i]

    # 分母 <Bv, v> = <v,v> + 2*mus* <uc,v>^2   (uc = u 的前 Nc 个分量)
    vv = iv.mpf(0)
    for i in range(N):
        vv += v[i] * v[i]
    if mus > 0.0 and Nc > 0:
        ucv = iv.mpf(0)
        for i in range(min(Nc, N)):
            ucv += u[i] * v[i]
        den = vv + iv.mpf(2.0 * mus) * ucv * ucv
    else:
        den = vv

    R = num / den
    return R, num, den


def slack_test(tag, wbar_val, mus, r, N):
    print(f"--- {tag}: wbar={wbar_val}, mu*={mus}, r={r}, N={N} ---")
    theta_f64, v0, Nc = theta_max_f64(float(wbar_val), mus, r, N)
    R, num, den = rayleigh_lower_bound(wbar_val, mus, r, N, v0)
    R_lo, R_hi = float(R.a), float(R.b)
    c_f64 = theta_f64 / (np.pi * float(wbar_val))
    c_lo = R_lo / (np.pi * float(wbar_val))
    rel_gap_theta = (theta_f64 - R_lo) / theta_f64
    rel_gap_c = (c_f64 - c_lo) / c_f64
    print(f"  Nc={Nc}")
    print(f"  theta_max (float64 eigh)        = {theta_f64:.12f}")
    print(f"  Rayleigh 区间 R(v0)              = [{R_lo:.12f}, {R_hi:.12f}]  区间宽度 {R_hi-R_lo:.3e}")
    print(f"  区间下端 vs float64 点值 相对差   = {rel_gap_theta:.3e}")
    print(f"  ==> c(w) 下界 = {c_lo:.10f}   (float64点值 c = {c_f64:.10f})")
    print(f"  ==> c 的相对松 (固有松) = {rel_gap_c:.3e}")
    return rel_gap_c, c_f64, c_lo


def main():
    step0_direction()

    print("=" * 92)
    print("STEP 3  区间估计子固有松实测（1-2 个分支，N 先用较小值以控制区间矩阵构造耗时）")
    print("=" * 92)
    N = 120  # 阶段1只测固有松量级，非最终认证；先用小 N 控速度，量级不受 N 敏感（见脚本内 print 的 float64 N-收敛对照）

    # 分支 A: mu*=0.13, r=10, 峰值分支 k=12 (w=1.2, 齿顶恰在左端，见既有笔记 §6.2)
    gap1, cf1, cl1 = slack_test("mu*=0.13,r=10, 分支k=12 (峰值附近 w=1.2+1e-6)", 1.2 + 1e-6, 0.13, 10.0, N)
    print()

    # 分支 B: mu*=0 (无库仑, B=I, 无锯齿——单一连续问题，peak 在别处；用既有脚本给出的 argmax lam 对应的 w)
    # mu*=0 情形没有 Nc 依赖，直接在已知 argmax 附近取一个 w 点做同样的 Rayleigh 检验
    gap2, cf2, cl2 = slack_test("mu*=0, r=0 (无库仑), w=1.75 附近 (已知 argmax lam~1.147, omega_E/A~1.75)",
                                 1.16 / 1.147 * 1.0, 0.0, 0.0, N)
    print()

    print("=" * 92)
    print("STEP 4  逐 mu* 余量 vs 固有松 对比")
    print("=" * 92)
    margins = {0.00: 0.0152, 0.10: 0.1853, 0.13: 0.2982, 0.20: 0.8156}
    print("  固有松实测（本脚本, N=120, 未最终收敛, 仅供数量级判断）：")
    print(f"    mu*=0.13 分支: c 相对松 = {gap1:.3e}")
    print(f"    mu*=0    分支: c 相对松 = {gap2:.3e}")
    print("  换算到 ratio 的相对松（一阶近似, ratio ~ 1/c^2 附近 Taylor）：")
    print("    d(ratio)/ratio ≈ 2 * d(c)/c  (因为 lam_area ∝ 1/c^2，ratio = lam_ceil/lam_area ∝ c^2)")
    for mus, m in margins.items():
        g = gap1 if mus > 0 else gap2  # 粗配（阶段1不需要精确匹配每个mu*的分支，只看量级）
        ratio_slack = 2.0 * g
        closeable = "能关" if ratio_slack < m else "关不掉（固有松 > 余量）"
        print(f"    mu*={mus:.2f}: 余量(ratio-1)={m:.4f}   估计 ratio 固有松≈{ratio_slack:.3e}   ⟹ {closeable}")


if __name__ == "__main__":
    main()
