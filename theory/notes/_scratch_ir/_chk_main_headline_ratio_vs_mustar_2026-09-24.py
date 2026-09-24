# -*- coding: utf-8 -*-
"""
把头条判词（我们 vs Leavens 面积天花板）从 mu*=0.13 一格扩到 6 格（2026-09-24）。

这是 `ceiling/same_information_class_comparison_vs_leavens_2026-09-24.md` §8 列的第一条
「未做」：mu* 依赖未扫。现在可以做，因为两个缺口同时补上了：
  (1) `c(mu*,r)` 的**锯齿**结构已被显式分支化（见
      `_chk_main_c_sawtooth_exact_branchwise_2026-09-24.py`）：`Nc=floor(r*w)` 的台阶
      只在 w=k/r 处跳，每个分支内 c(w) 光滑、**齿顶恰在分支左端 w=k/r**
      ⟹ 逐 k 枚举即穷尽，不再是「扫到的最高齿」。
  (2) 固定 M2 轴上的临界 lam*(mu*) 已由本会话独立重算
      （`_chk_main_mustar_axes_independent_recheck_2026-09-24.py`，与 mustar-scan 相对差 ~2e-5）。

两侧都用**同一个** mu*、同一个 r、同一个 M2 ⟹ 比值 lam*/lam_area 是同信息类、同约定的判词。

━━━ 方向审计（承重）━━━
* `lam_area = 1/(pi^2 c^2 M2)`，`d lam_area/dc < 0`。
  c 取「逐齿枚举的 max」⟹ 若还有未枚举的更高齿，c 只会更大 ⟹ lam_area 只会更**小**
  ⟹ 比值只会更**大** ⟹ 现报的比值是**下估**，对「我们更强」这个判词**保守**。
* `lam*` 这一列是 float64 点值（不是认证下端）⟹ 本表的比值是 numerical-evidence 级；
  只有 mu*=0.13 一格另有 certificate 级族证书 [3.29,3.31]，其比值单独列出。

严格性：全部 numerical-evidence（float64, 广义 eigh）。分支结构的「齿顶在左端」是实测，
不是证明——故仍标 numerical-evidence，不上升到 analytic。
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
    Khat, u = pieces(wbar, N)
    B = np.eye(N)
    Nc = int(np.floor(r * wbar))
    if mus > 0.0 and Nc > 0:
        uc = np.where(np.arange(N) < Nc, u, 0.0)
        B = B + 2.0 * mus * np.outer(uc, uc)
    th = float(eigh(Khat, B, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])
    return (1.0 / th) if th > 0 else np.inf, Nc


def c_max_branchwise(mus, r, N, kmax=40, npts=40):
    """逐齿枚举 max_w 1/(pi*lam*(w)*w)。mu*=0 时无台阶，退化为一次光滑扫描。"""
    best = (-1.0, None, None, None)
    if mus == 0.0:
        ws = np.linspace(0.05, 8.0, 1200)
        for w in ws:
            lam, Nc = lam_root_exact(w, mus, r, N)
            if not np.isfinite(lam) or lam <= 0:
                continue
            c = 1.0 / (np.pi * lam * w)
            if c > best[0]:
                best = (c, w, lam, Nc)
        h = 8.0 / 1199
        for w in np.linspace(best[1] - h, best[1] + h, 61):
            lam, Nc = lam_root_exact(w, mus, r, N)
            c = 1.0 / (np.pi * lam * w)
            if c > best[0]:
                best = (c, w, lam, Nc)
        return best
    for k in range(1, kmax + 1):
        lo, hi = k / r + 1e-9, (k + 1) / r - 1e-9
        for w in np.linspace(lo, hi, npts):
            lam, Nc = lam_root_exact(w, mus, r, N)
            if not np.isfinite(lam) or lam <= 0:
                continue
            c = 1.0 / (np.pi * lam * w)
            if c > best[0]:
                best = (c, w, lam, Nc)
    # 齿顶实测在分支左端 ⟹ 在最优 k 的左端加密
    kbest = int(np.floor(r * best[1]))
    for w in np.linspace(kbest / r + 1e-12, kbest / r + 0.5 / r, 81):
        lam, Nc = lam_root_exact(w, mus, r, N)
        c = 1.0 / (np.pi * lam * w)
        if c > best[0]:
            best = (c, w, lam, Nc)
    return best


def solve_fixed_M2(mus, r, N, M2=M2_BASE, lo=0.10, hi=4.0, nbis=64):
    f = lambda w: lam_root_exact(w, mus, r, N)[0] - M2 / (w * w)
    a, b = lo, hi
    fa = f(a)
    for _ in range(nbis):
        m = 0.5 * (a + b)
        fm = f(m)
        if not np.isfinite(fm) or (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    w = 0.5 * (a + b)
    return M2 / (w * w), w, int(np.floor(r * w))


def main():
    r, N = 10.0, 400
    MUS = (0.00, 0.08, 0.10, 0.13, 0.16, 0.20)
    print("=" * 96)
    print(f"头条表 x mu*：同信息类 (lam, M2) 下 我们 vs Leavens 面积天花板"
          f"  (M2={M2_BASE}, r={r}, N={N})")
    print("=" * 96)
    print("  mu*     c(mu*,r)      齿顶lam  齿Nc |  lam*(我们)  Nc |   lam_area   "
          "| lam*/lam_area")
    print("  " + "-" * 92)
    out = {}
    for mus in MUS:
        c, wc, lc, ncc = c_max_branchwise(mus, r, N)
        lam_s, w_s, nc_s = solve_fixed_M2(mus, r, N)
        la = 1.0 / (np.pi ** 2 * c * c * M2_BASE)
        out[mus] = (c, lam_s, la, lam_s / la)
        print(f"  {mus:4.2f}   {c:.9f}  {lc:9.5f}  {ncc:3d} |  {lam_s:9.5f}  {nc_s:3d} |"
              f"  {la:9.6f}  |   {lam_s/la:8.4f}")
    print()
    print("  ⟹ 比值在 mu* in [0, 0.20] 上的范围： "
          f"{min(v[3] for v in out.values()):.4f} – {max(v[3] for v in out.values()):.4f}")

    print()
    print("=" * 96)
    print("GATE M1  mu*=0 的结构检验：无库仑 ⟹ 无台阶 ⟹ c 应是光滑函数的 max")
    print("=" * 96)
    c0, w0, l0, _ = c_max_branchwise(0.0, r, N)
    print(f"  c(0) = {c0:.9f} @ (lam={l0:.6f}, w={w0:.6f})")
    print("  邻域三点（应无跳变）：")
    for dw in (-0.02, -0.01, 0.0, 0.01, 0.02):
        lam, _ = lam_root_exact(w0 + dw, 0.0, r, N)
        print(f"    w={w0+dw:.6f}  lam*={lam:.8f}  c={1.0/(np.pi*lam*(w0+dw)):.9f}")

    print()
    print("=" * 96)
    print("GATE M2  mu*=0.13 的 certificate 级比值（唯一有族证书的一格）")
    print("=" * 96)
    c13, _, _, _ = c_max_branchwise(0.13, r, N)
    la13 = 1.0 / (np.pi ** 2 * c13 * c13 * M2_BASE)
    for tag, lo in (("认证下端 3.29", 3.29), ("Einstein 根(float)", out[0.13][1]),
                    ("认证上端 3.31", 3.31)):
        print(f"  {tag:22s} / lam_area={la13:.6f}  =  {lo/la13:.4f}")
    print("  ⟹ certificate 级可报的是**下端** 3.29 那一行（比值的认证下估）。")

    print()
    print("=" * 96)
    print("GATE M3  N 收敛（mu*=0.13 与 mu*=0.20 两格的比值）")
    print("=" * 96)
    for NN in (300, 400, 600):
        row = []
        for mus in (0.13, 0.20):
            c, _, _, _ = c_max_branchwise(mus, r, NN)
            lam_s, _, _ = solve_fixed_M2(mus, r, NN)
            la = 1.0 / (np.pi ** 2 * c * c * M2_BASE)
            row.append(f"mu*={mus}: ratio={lam_s/la:.5f}")
        print(f"  N={NN:4d}   " + "   ".join(row))


if __name__ == "__main__":
    main()
