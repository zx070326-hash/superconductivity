# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）：`_chk_main_phiA_general_spectrum_sup_2026-09-23.py` 的三道后续门。

上一个探针的意外发现（**不是我预期的形状**）：固定 (lam, wbar2) 的两原子族里，
  · ceiling = mu1 - Delta 的族内 sup **落在 Einstein（退化成员）上**；
  · 整族衰减**几乎全部来自 mu1**（1.3795 -> 0.655），**Delta 几乎不动**（D/D_E ~ 1.00）。
若属实，(Phi-A) 的一般谱化就**不需要** Delta 的一致下界那种精细工作，只需要
  sup mu1 = mu1^E  （AD1975 mu*=0 极值性）  +  Delta 的一个**粗**下界。

本脚本查三件事：

GATE X（**解析归约的数值前提**）：mu1 是测度 nu := alpha^2F dw 的**凸**泛函？
  K[nu] 对 nu 线性（G(k)=∫ 2 w/(w^2+(2 pi T k)^2) dnu），lam_max 是对称矩阵的凸函数
  ⟹ mu1 = lam_max(K[nu]) 对 nu 凸。**这是 analytic**，本门只做数值负对照（找不到违反）。
  意义：约束集 {nu>=0 : ∫w^{-1}dnu = lam/2, ∫w dnu = lam wbar2^2/2} 是**两个线性约束**
  ⟹ 其**极点支撑 <= 2 个原子**（Richter–Rogosinski / Carathéodory）；
  凸泛函在凸紧集上的 **max 在极点取到** ⟹ sup mu1 只需在**两原子族**上找。
  ⟹ 上一个探针扫的那个族**恰好就是极点集**，不是抽样。

GATE Y：min Delta/Delta_E 在整个两原子族上是多少？（预算要求 >= 0.76）

GATE Z（负对照）：三原子成员能否超过两原子族的 sup？按上面的归约**不应该**能。
  若能，说明归约错了（或我实现错了）。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh
from itertools import product

import importlib.util as _il
_s = _il.spec_from_file_location(
    "sup", os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "_chk_main_phiA_general_spectrum_sup_2026-09-23.py"))
_m = _il.module_from_spec(_s); _s.loader.exec_module(_m)
Kfull, bundle, members, einstein_root = _m.Kfull, _m.bundle, _m.members, _m.einstein_root


def mu1_of(lams, wbars, N):
    K, _ = Kfull(lams, wbars, N)
    return float(eigh(K, subset_by_index=[N - 1, N - 1], eigvals_only=True)[0])


def main():
    N, mus, r = 300, 0.13, 10.0
    rng = np.random.default_rng(20260923)

    # ───────────── GATE X: mu1 对测度的凸性（数值负对照）─────────────
    print("=" * 92)
    print("GATE X  mu1 = lam_max(K[nu]) 对 nu 凸？  查 mu1(t*a+(1-t)*b) <= t*mu1(a)+(1-t)*mu1(b)")
    print("        （analytic：K 对 nu 线性 + lam_max 凸。本门只找反例。）")
    worst = -1e9
    nbad = 0
    for _ in range(300):
        na, nb = rng.integers(1, 4), rng.integers(1, 4)
        la, wa = rng.uniform(0.1, 2.0, na), rng.uniform(0.2, 8.0, na)
        lb, wb = rng.uniform(0.1, 2.0, nb), rng.uniform(0.2, 8.0, nb)
        t = rng.uniform(0.05, 0.95)
        m_a, m_b = mu1_of(la, wa, N), mu1_of(lb, wb, N)
        m_mix = mu1_of(np.concatenate([t * la, (1 - t) * lb]),
                       np.concatenate([wa, wb]), N)
        d = m_mix - (t * m_a + (1 - t) * m_b)
        worst = max(worst, d)
        nbad += d > 1e-9
    print(f"        300 组随机混合：违反 {nbad}/300，最劣越界 = {worst:.3e}"
          f"   {'(凸性成立)' if nbad == 0 else '** 有违反 **'}")
    print()

    # ───────────── GATE Y / Z ─────────────
    for lam in (2.0, 3.2953, 4.5):
        wE = einstein_root(lam, mus, r, N)
        Nc = int(np.floor(r * wE))
        l1E, cE, dE, gE, m1E, c1E, _ = bundle(np.array([lam]), np.array([wE]), mus, Nc, N)
        print("=" * 92)
        print(f"lam = {lam:.4f}   wbarE* = {wE:.6f}   N_c = {Nc}   "
              f"Einstein: ceiling={cE:.6f} Delta={dE:.6f} mu1={m1E:.6f}")

        # GATE Y: 两原子族细扫 min D/D_E, max ceiling
        dmin, dmin_at = 1e9, None
        cmax, cmax_at = -1e9, None
        m1max = -1e9
        npa = n = 0
        xs = np.concatenate([np.linspace(0.02, 0.98, 25)])
        rhos = np.concatenate([[1.0], np.geomspace(1.02, 40.0, 24)])
        for x, rho in product(xs, rhos):
            lams, wbars = members(lam, wE, x, rho)
            _, c, d, g, m1, _, ok = bundle(lams, wbars, mus, Nc, N)
            n += 1
            if not ok:
                npa += 1
                continue
            if d / dE < dmin:
                dmin, dmin_at = d / dE, (x, rho)
            if c > cmax:
                cmax, cmax_at = c, (x, rho)
            m1max = max(m1max, m1)
        print(f"  GATE Y  两原子族 {n} 成员（含 rho=1 定向命中退化集）：(PA) 失效 {npa}/{n}")
        print(f"          min Delta/Delta_E = {dmin:.4f} at (x,rho)={dmin_at[0]:.3f},{dmin_at[1]:.2f}"
              f"   预算门槛 0.76  ->  {'PASS' if dmin >= 0.76 else '**FAIL**'}")
        print(f"          max ceiling = {cmax:.9f} at (x,rho)={cmax_at[0]:.3f},{cmax_at[1]:.2f}"
              f"   Einstein {cE:.9f}   超出 {cmax-cE:+.3e}")
        print(f"          max mu1     = {m1max:.9f}   Einstein {m1E:.9f}   超出 {m1max-m1E:+.3e}")

        # GATE Z: 三原子负对照（按归约不应超过两原子 sup）
        zmax, zat = -1e9, None
        zm1 = -1e9
        nz = 0
        for _ in range(600):
            a = rng.dirichlet([1.0, 1.0, 1.0])
            rr = np.exp(rng.uniform(np.log(0.05), np.log(6.0), 3))
            # 缩放频率使 sum a_j w_j^2 = wE^2 （固定 wbar2）
            s = wE / np.sqrt(float(a @ (rr * rr)))
            lams, wbars = lam * a, s * rr
            _, c, d, g, m1, _, ok = bundle(lams, wbars, mus, Nc, N)
            nz += 1
            if not ok:
                continue
            if c > zmax:
                zmax, zat = c, (np.round(a, 3), np.round(s * rr, 3))
            zm1 = max(zm1, m1)
        print(f"  GATE Z  三原子随机 {nz} 组：max ceiling = {zmax:.9f}"
              f"   vs 两原子 sup {cmax:.9f}   超出 {zmax-cmax:+.3e}"
              f"   {'OK（未超，与极点归约相容）' if zmax <= cmax + 1e-9 else '** 超了，归约有问题 **'}")
        print(f"          max mu1 = {zm1:.9f}  vs Einstein {m1E:.9f}  超出 {zm1-m1E:+.3e}")
        print()


if __name__ == "__main__":
    main()
