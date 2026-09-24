# -*- coding: utf-8 -*-
"""
主会话（2026-09-23）：把承重量从 (Phi-A) 的 ceiling 换回**精确 lam1**，在极点族上细扫。

理由（本轮的关键转向）：`_chk_main_phiA_extreme_point_reduction_2026-09-23.py` 的 GATE Y
证伪了我自己在订正笔记 §4 里写的路线（「须给 Delta 的一致下界，>=76%」）——
两原子族上 **min Delta/Delta_E = 0.042**，根本没有有用的一致下界。
但同一批扫描里 **ceiling 的族内 sup 仍精确落在 Einstein**（1e-15）⟹ 拆成
mu1 与 Delta 分别界**不是**正确的分解：Delta 塌的地方 mu1 塌得更快。

正确的分解根本不需要 (Phi-A)：固定 (T, N_c) 时
        M[nu] = K[nu] - 2 mu* u_c u_c^T
对测度 nu **仿射**（K 线性、库仑项与 nu 无关，因为 w_c = r*wbar2 对全族相同）
⟹ lam1(M[nu]) = lam_max(仿射) 是 nu 的**凸**泛函
⟹ 在两个线性矩约束下的凸紧集上，max 在**极点**取到，而极点支撑 <= 2 原子
   （Bauer 极大原理 + Richter–Rogosinski）。
⟹ **无穷维 sup 精确塌成 2 参数问题**，且承重量是**精确 lam1**、不是任何上界
   ⟹ (PA) 前提、Weyl 松弛、Delta 的一致下界**全部不需要**。

本脚本在极点族（两原子，含定向命中的退化成员 rho=1）上细扫 **lam1**，查
  max_{x,rho} lam1  是否 = lam1(Einstein)，以及超出多少。
这是新路线第 4 步（2 参数问题）的 float64 先验；升级到 certificate 需 2D 区间细分。

严格性：凸性与极点归约 analytic（本脚本只做数值负对照）；一切数值 numerical-evidence。
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from itertools import product

import importlib.util as _il
_d = os.path.dirname(os.path.abspath(__file__))
_s = _il.spec_from_file_location(
    "sup", os.path.join(_d, "_chk_main_phiA_general_spectrum_sup_2026-09-23.py"))
_m = _il.module_from_spec(_s); _s.loader.exec_module(_m)
bundle, members, einstein_root = _m.bundle, _m.members, _m.einstein_root


def main():
    N, r = 300, 10.0
    print("=" * 96)
    print("极点族（两原子，固定 lam 与 wbar2）上 **精确 lam1** 的族内 sup")
    print("约定 w_c = r*wbar2 ==> N_c 对全族相同 ==> M[nu] 对 nu 仿射（凸性归约的前提）")
    print()
    for mus in (0.0, 0.10, 0.13, 0.20):
        for lam in (1.5, 2.0, 3.2953, 4.5, 6.0):
            wE = einstein_root(lam, mus, r, N)
            if wE is None:
                print(f"mu*={mus} lam={lam}: 无根，跳过")
                continue
            Nc = int(np.floor(r * wE))
            l1E, *_ = bundle(np.array([lam]), np.array([wE]), mus, Nc, N)
            best, bat = -1e9, None
            n = 0
            xs = np.linspace(0.01, 0.99, 33)
            rhos = np.concatenate([[1.0], np.geomspace(1.005, 60.0, 32)])
            for x, rho in product(xs, rhos):
                lams, wbars = members(lam, wE, x, rho)
                l1, *_ = bundle(lams, wbars, mus, Nc, N)
                n += 1
                if l1 > best:
                    best, bat = l1, (x, rho)
            flag = "Einstein 是 sup" if best <= l1E + 1e-12 else "** Einstein 不是 sup **"
            print(f"mu*={mus:<5} lam={lam:<7.4f} wbarE*={wE:.5f} Nc={Nc:<4d} "
                  f"lam1_E={l1E:.12f}  max lam1={best:.12f} at (x={bat[0]:.3f},rho={bat[1]:.3f})"
                  f"  超出={best-l1E:+.2e}  [{n} 成员]  {flag}")
        print()


if __name__ == "__main__":
    main()
