# -*- coding: utf-8 -*-
"""
主会话从零独立复核 mustar-scan 的两张表（2026-09-24）。**不 import 它的脚本、不看它的日志**，
只用本会话自己的 `lam_root_exact`（广义 eigh 精确最小根，见
`_chk_main_c_sawtooth_exact_branchwise_2026-09-24.py` 里的引理）。

两张表是两条**不同的轴**，必须分开核（混起来是上一轮已记的事故形状）：

(I) 主表（**M2 固定** = 1.453151466 ⟺ eta/M=9.172 @ 300K）：sweep mu*，求临界 lam。
    M2 固定 ⟹ w = sqrt(M2/lam) 被 lam 绑死 ⟹ 解一维方程
        lam_root_exact(w; mu*) = M2 / w^2
    左边随 w 增而降（w 大=声子硬=需要的耦合小），右边随 w 增而降得更快 ⟹ 单交点。
    解出 w* 后 lam* = M2/w*^2。

(II) C(mu*) 轴（**lam 固定** = 3.2953）：sweep mu*，求 w*，C = 1/(2*pi*w*sqrt(lam))。
    推导：kTc = C*sqrt(lam<w^2>) = C*sqrt(M2)；松原单位下 kTc -> 1/(2pi)、sqrt(M2)=sqrt(lam)*w
    ⟹ C = 1/(2 pi w sqrt(lam))。对照已并入的 analytic 常数 C_inf = 0.1827262477890，
    一致性门：C(mu*) < C_inf 必须处处成立（mu*>0 只会压低 Tc）。

方向审计：本脚本给的是 float64 点值定位（numerical-evidence），**不是**认证界，
因此只能用来与 mustar-scan 的 float64 列对照、以及做一致性门；
它**不**复核 mustar-scan 的 certificate 级 PASS/FAIL 夹逼（那要另跑分支定界）。

严格性：全部 numerical-evidence（float64, 广义 eigh, 二分）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

M2_BASE = 1.453151466
C_INF = 0.1827262477890


def pieces(wbar, N):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    k = np.arange(2 * N + 2)
    G = wbar * wbar / (wbar * wbar + k * k)
    Kp, Km = G[np.abs(i - j)], G[i + j + 1]
    S = (Kp - Km).sum(axis=1)
    return (Kp + Km) * np.outer(u, u) - np.diag(u * u * S), u


def lam_root_exact(wbar, mus, r, N):
    """最小 lam 使 lam1(lam*Khat - 2mu* u_c u_c^T) = 1，等于 1/theta_max。"""
    Khat, u = pieces(wbar, N)
    B = np.eye(N)
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(N) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    th = float(eigh(Khat, B, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])
    return (1.0 / th) if th > 0 else np.inf, Nc


def solve_fixed_M2(mus, r, N, M2=M2_BASE, lo=0.10, hi=4.0, nbis=64):
    """(I) 解 lam_root_exact(w) = M2/w^2。返回 (lam*, w*, Nc)。"""
    f = lambda w: lam_root_exact(w, mus, r, N)[0] - M2 / (w * w)
    a, b = lo, hi
    fa = f(a)
    if not np.isfinite(fa):
        return None
    for _ in range(nbis):
        m = 0.5 * (a + b)
        fm = f(m)
        if not np.isfinite(fm) or (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    w = 0.5 * (a + b)
    return M2 / (w * w), w, int(np.floor(r * w))


def solve_fixed_lam(lam, mus, r, N, lo=0.05, hi=12.0, nbis=64):
    """(II) 解 lam_root_exact(w) = lam。lam_root_exact 在 w 上递减 ⟹ 直接二分。"""
    a, b = lo, hi
    for _ in range(nbis):
        m = 0.5 * (a + b)
        if lam_root_exact(m, mus, r, N)[0] > lam:
            a = m
        else:
            b = m
    w = 0.5 * (a + b)
    return w, int(np.floor(r * w))


def main():
    r, N = 10.0, 400
    MUS = (0.00, 0.08, 0.10, 0.13, 0.16, 0.20)
    # mustar-scan 报的数（只用于**对照**，不作为我的计算输入）
    THEIRS_I = {0.00: 1.328106, 0.08: 2.323594, 0.10: 2.662098,
                0.13: 3.295300, 0.16: 4.149504,
                0.20: 5.858499}   # ⚠ 0.20 这个是 mustar-scan 报的**最大**根；物理 λ_ceil 是
                                 # inf-根 5.797500（N_c 台阶使 λ₁−1 变号 3 次，2026-09-24 订正）。
                                 # 本脚本只把它当作「复现 mustar-scan 的输出」的对照，不作 λ_ceil 用。
    THEIRS_NC_I = {0.00: 10, 0.08: 7, 0.10: 7, 0.13: 6, 0.16: 5, 0.20: 4}
    THEIRS_II = {0.00: (0.5376450, 0.1630711), 0.08: (0.6172714, 0.1420354),
                 0.10: (0.6361974, 0.1378100), 0.13: (0.6640610, 0.1320276),
                 0.16: (0.6914153, 0.1268042), 0.20: (0.7244636, 0.1210197)}

    print("=" * 84)
    print(f"(I) 主表复核：M2 固定 = {M2_BASE}, r={r}, N={N}  —— 我的 lam* vs mustar-scan")
    print("=" * 84)
    print("   mu*      我的 lam*      我的 w*    我Nc  他Nc      他的 lam*     相对差")
    for mus in MUS:
        res = solve_fixed_M2(mus, r, N)
        lam, w, Nc = res
        th = THEIRS_I[mus]
        print(f"  {mus:4.2f}   {lam:12.6f}  {w:10.7f}  {Nc:4d}  {THEIRS_NC_I[mus]:4d}"
              f"   {th:12.6f}   {abs(lam-th)/th:.2e}")

    print()
    print("=" * 84)
    print(f"(II) C(mu*) 轴复核：lam 固定 = 3.2953, r={r}, N={N}")
    print("=" * 84)
    print("   mu*      我的 w*      Nc      我的 M2        我的 C        他的 C      C/C_inf")
    lam = 3.2953
    bad = 0
    for mus in MUS:
        w, Nc = solve_fixed_lam(lam, mus, r, N)
        M2 = lam * w * w
        C = 1.0 / (2.0 * np.pi * w * np.sqrt(lam))
        wt, Ct = THEIRS_II[mus]
        ratio = C / C_INF
        flag = ""
        if ratio >= 1.0:
            flag = "  <<< CONSISTENCY ALARM"
            bad += 1
        print(f"  {mus:4.2f}   {w:10.7f}  {Nc:4d}   {M2:11.7f}   {C:11.7f}   "
              f"{Ct:11.7f}   {ratio:8.5f}{flag}")
    print(f"\n  一致性门（C < C_inf = {C_INF}）：{'PASS' if bad==0 else 'FAIL'}  "
          f"（{len(MUS)-bad}/{len(MUS)}）")

    print()
    print("=" * 84)
    print("(III) 锐化量（他报 C(0)-C(0.13)=+0.0310435 = 1.235x；C(0)-C(0.20)=+0.0420514 = 1.348x）")
    print("=" * 84)
    w0, _ = solve_fixed_lam(lam, 0.00, r, N)
    C0 = 1.0 / (2.0 * np.pi * w0 * np.sqrt(lam))
    for mus in (0.13, 0.20):
        w, _ = solve_fixed_lam(lam, mus, r, N)
        C = 1.0 / (2.0 * np.pi * w * np.sqrt(lam))
        print(f"  C(0)-C({mus}) = {C0-C:+.7f}   C(0)/C({mus}) = {C0/C:.4f}x")

    print()
    print("=" * 84)
    print("(IV) 自由度审计：mu*=0 时 Nc 是否进判据？")
    print("=" * 84)
    print("  mu*=0 ⟹ 秩一项系数恰为 0 ⟹ B = I ⟹ Nc **在代数上**根本不出现。")
    print("  ⟹ 「实测 Nc in {0,1,6,...} 给逐位相同的 lam1」自由度为 0，是恒等式不是检验。")
    w_ = 0.5376450
    for Nc_try in (0, 1, 6, 50, 400):
        Khat, u = pieces(w_, N)
        uc = np.where(np.arange(N) < Nc_try, u, 0.0)
        B = np.eye(N) + 2.0 * 0.0 * np.outer(uc, uc)
        th = float(eigh(Khat, B, subset_by_index=[N-1, N-1], eigvals_only=True)[0])
        print(f"    Nc={Nc_try:4d}  theta_max={th:.15f}  (mu*=0)")
    print("  ⟹ 印出来的确是逐位相同，但这只重复了 `0 * X = 0`。")
    print("     真正有内容的是 mu*>0 的 Nc 台阶普查（主表 Nc 列 [10,7,7,6,5,4]）。")


if __name__ == "__main__":
    main()
