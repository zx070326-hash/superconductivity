# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）：凸性 + 极点归约的四道门（承重量 = **精确 lam1**）。

为什么要重跑：`_chk_main_phiA_extreme_point_reduction_2026-09-23.py` 的 GATE X/Z 测的是
**(Phi-A) 的 ceiling 与 mu1**。本轮承重量已换成**精确 lam1**
（CLAUDE.md：「改了链条上任何一个对象的性质，必须把链条上每个用到该性质的箭头重走一遍」）
⟹ 凸性与极点归约必须在 lam1 上重测，旧的通过率不可转引。

待验的解析链（每一步都标出本脚本测哪一条）：
  (A1) 固定 (T, N_c) 时 M[nu] = K[nu] - 2 mu* u_c u_c^T 对测度 nu **仿射**   -> GATE C4（凸性）
  (A2) lam_max(仿射族) 对 nu **凸**                                          -> GATE C4
  (A3) 约束集 {nu>=0 : ∫w^{-1}dnu = lam/2, ∫w dnu = lam wbar2^2/2} 紧凸       -> 解析（本脚本不测）
  (A4) Bauer 极大原理 ⟹ max 在极点取到                                       -> 间接由 C1/C2 负对照
  (A5) Richter–Rogosinski：两个线性约束 ⟹ 极点支撑 <= 2 原子                 -> 间接由 C1/C2 负对照
  (A6) 两原子族上的 sup 落在退化成员 rho=1（= Einstein）                      -> 已由
        `_chk_main_exact_lam1_extreme_family_2026-09-23.py` 在 20 组 (mu*,lam) 上测过（<=1.8e-15）

GATE C1  三原子随机成员（固定 (lam, wbar2)）能否超过 Einstein 的 lam1？——按归约不应能。
GATE C2  四原子 + **连续谱**（截断高斯 / 方窗，各自重标度固定 wbar2）能否超过？——同上。
         连续谱是比多原子更强的负对照：它不在极点集里，却在约束集里。
GATE C3  **投毒负对照**：把库仑截断换成 w_c = r*w_max（逐成员不同的 N_c）⟹ 破坏仿射性。
         若此时 Einstein 仍恒为 sup，说明本组门**没有分辨力**（凸性不是结论的来源）；
         若此时能被超过，说明 w_c = r*wbar2 这个约定是**承重**的，门有牙。
GATE C4  精确 lam1 对 nu 的凸性（随机混合，固定 T 与 N_c）——(A1)+(A2) 的数值负对照。

严格性：(A1)-(A5) analytic（本脚本只做负对照）；一切数值为 numerical-evidence（float64, N=300）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh
from itertools import product

import importlib.util as _il
_d = os.path.dirname(os.path.abspath(__file__))
_s = _il.spec_from_file_location(
    "sup", os.path.join(_d, "_chk_main_phiA_general_spectrum_sup_2026-09-23.py"))
_m = _il.module_from_spec(_s); _s.loader.exec_module(_m)
Kfull, bundle, members, einstein_root = _m.Kfull, _m.bundle, _m.members, _m.einstein_root


def lam1(lams, wbars, mus, Nc, N):
    K, u = Kfull(np.asarray(lams, float), np.asarray(wbars, float), N)
    uc = np.where(np.arange(N) < Nc, u, 0.0)
    M = K - 2.0 * mus * np.outer(uc, uc)
    return float(eigh(M, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])


def rescale_to(lams, wbars, lam, wbarE):
    """把任意 (lams, wbars) 重标度到固定 (lam, wbar2)：先归一权重，再整体缩频率。"""
    a = np.asarray(lams, float)
    a = a / a.sum()
    w = np.asarray(wbars, float)
    s = wbarE / np.sqrt(float(a @ (w * w)))
    return lam * a, s * w


def main():
    N, r = 300, 10.0
    rng = np.random.default_rng(20260923)
    CASES = [(0.0, 2.0), (0.13, 2.0), (0.13, 3.2953), (0.13, 4.5), (0.20, 6.0)]

    # ───────────────── GATE C4: 精确 lam1 对 nu 的凸性 ─────────────────
    print("=" * 100)
    print("GATE C4  精确 lam1 对 nu 凸？  查 lam1(t*a+(1-t)*b) <= t*lam1(a)+(1-t)*lam1(b)")
    print("         固定 T 与 N_c（否则 M 对 nu 不仿射）。analytic：M 仿射 + lam_max 凸。")
    for mus, Nc in ((0.0, 8), (0.13, 8), (0.13, 20), (0.30, 5)):
        worst, nbad, n = -1e9, 0, 0
        for _ in range(400):
            na, nb = int(rng.integers(1, 5)), int(rng.integers(1, 5))
            la, wa = rng.uniform(0.05, 3.0, na), np.exp(rng.uniform(np.log(0.05), np.log(20.0), na))
            lb, wb = rng.uniform(0.05, 3.0, nb), np.exp(rng.uniform(np.log(0.05), np.log(20.0), nb))
            t = rng.uniform(0.05, 0.95)
            fa, fb = lam1(la, wa, mus, Nc, N), lam1(lb, wb, mus, Nc, N)
            fm = lam1(np.concatenate([t * la, (1 - t) * lb]), np.concatenate([wa, wb]), mus, Nc, N)
            d = fm - (t * fa + (1 - t) * fb)
            worst = max(worst, d); nbad += d > 1e-9; n += 1
        print(f"         mu*={mus:<5} N_c={Nc:<3d}  违反 {nbad}/{n}  最劣越界 {worst:+.3e}"
              f"   {'(凸)' if nbad == 0 else '** 有违反 **'}")
    print()

    # ───────────────── GATE C1 / C2 / C3 ─────────────────
    for mus, lam in CASES:
        wE = einstein_root(lam, mus, r, N)
        Nc = int(np.floor(r * wE))
        l1E = lam1([lam], [wE], mus, Nc, N)
        print("=" * 100)
        print(f"mu*={mus}  lam={lam:.4f}   Einstein wbarE*={wE:.6f}  N_c={Nc}  lam1_E={l1E:.12f}")

        # C1: 三原子
        best, bat, n = -1e9, None, 0
        for _ in range(800):
            a = rng.dirichlet([1.0, 1.0, 1.0])
            w = np.exp(rng.uniform(np.log(0.05), np.log(8.0), 3))
            L, W = rescale_to(a, w, lam, wE)
            v = lam1(L, W, mus, Nc, N); n += 1
            if v > best:
                best, bat = v, (np.round(a, 3), np.round(W, 3))
        print(f"  GATE C1  三原子随机 {n} 组：max lam1 = {best:.12f}  超出 Einstein {best-l1E:+.3e}"
              f"   {'OK（未超）' if best <= l1E + 1e-9 else '** 超了 **'}")

        # C2: 四原子 + 连续谱
        best2, tag, n2 = -1e9, None, 0
        for _ in range(400):
            a = rng.dirichlet(np.ones(4))
            w = np.exp(rng.uniform(np.log(0.05), np.log(8.0), 4))
            L, W = rescale_to(a, w, lam, wE)
            v = lam1(L, W, mus, Nc, N); n2 += 1
            if v > best2:
                best2, tag = v, "4-atom"
        for _ in range(200):   # 截断高斯（用 60 点离散逼近连续谱）
            c0 = np.exp(rng.uniform(np.log(0.2), np.log(4.0)))
            sg = c0 * rng.uniform(0.05, 1.2)
            g = np.linspace(max(1e-3, c0 - 4 * sg), c0 + 4 * sg, 60)
            a = np.exp(-0.5 * ((g - c0) / sg) ** 2)
            L, W = rescale_to(a, g, lam, wE)
            v = lam1(L, W, mus, Nc, N); n2 += 1
            if v > best2:
                best2, tag = v, "gauss"
        for _ in range(200):   # 方窗
            lo = np.exp(rng.uniform(np.log(0.03), np.log(3.0)))
            hi = lo * np.exp(rng.uniform(np.log(1.05), np.log(30.0)))
            g = np.linspace(lo, hi, 60)
            L, W = rescale_to(np.ones(60), g, lam, wE)
            v = lam1(L, W, mus, Nc, N); n2 += 1
            if v > best2:
                best2, tag = v, "box"
        print(f"  GATE C2  4原子+连续谱 {n2} 组：max lam1 = {best2:.12f}  超出 {best2-l1E:+.3e}"
              f"   (argmax 类型 {tag})   {'OK（未超）' if best2 <= l1E + 1e-9 else '** 超了 **'}")

        # C3: 投毒——破坏仿射性（w_c = r*w_max ⟹ 逐成员不同的 N_c）
        best3, bat3, n3 = -1e9, None, 0
        xs = np.linspace(0.02, 0.98, 25)
        rhos = np.concatenate([[1.0], np.geomspace(1.02, 40.0, 24)])
        for x, rho in product(xs, rhos):
            L, W = members(lam, wE, x, rho)
            Nc_bad = int(np.floor(r * W.max()))      # 逐成员不同 ⟹ M 不再对 nu 仿射
            v = lam1(L, W, mus, Nc_bad, N); n3 += 1
            if v > best3:
                best3, bat3 = v, (x, rho)
        l1E_bad = lam1([lam], [wE], mus, int(np.floor(r * wE)), N)
        vio = best3 > l1E_bad + 1e-9
        print(f"  GATE C3  投毒（w_c=r*w_max，N_c 逐成员变）{n3} 成员：max lam1 = {best3:.12f}"
              f"  超出 {best3-l1E_bad:+.3e} at (x,rho)=({bat3[0]:.2f},{bat3[1]:.2f})")
        print(f"           {'** 被超过 ⟹ 门有分辨力：w_c=r*wbar2 的仿射约定是承重的 **' if vio else '未被超过 ⟹ 本门无分辨力（警告：凸性可能不是结论来源）'}")
        print()


if __name__ == "__main__":
    main()
