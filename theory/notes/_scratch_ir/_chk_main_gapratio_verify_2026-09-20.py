# -*- coding: utf-8 -*-
"""
主会话从零复核 gap-ratio-observable 的两条承重陈述（2026-09-20）。

争议 1：它报「弱耦合门禁未过，lambda=0.5 给 ~3.63、lambda=0.1 预期 ~4.8（+36%）」。
        但这两个数都是**外推**，且方向（lambda 下降 ratio 上升）与 BCS 极限反向。
        仓库 memory 已记：lambda=0.3 给够 N 收敛到 3.5547@N=8000，对上 Carbotte 3.5470。
        ==> 检验：把 N 给够，ratio 是否单调趋向 BCS 3.528？

争议 2：它的预算脚本用 `spec.allen_dynes(mu_star)` 当分母 Tc，而分子 delta[0] 来自
        Eliashberg 自洽解。**两个不同理论**。==> 检验：换成自洽的 solver.tc()，R 变多少？

严格性：numerical-evidence（float64）。
"""
import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, r"D:\physics\superconductivity")
import numpy as np
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

BCS = 2 * np.pi / np.exp(np.euler_gamma)
print(f"BCS 2D/kTc = {BCS:.6f}")
print()

# =====================================================================
# 争议 1：弱耦合极限的 N 收敛（mu*=0，Einstein w_E=50 meV，与 09-17 脚本同约定）
# =====================================================================
print("=" * 88)
print("CHECK 1: weak-coupling limit, N convergence  (mu*=0, w_E=50meV, w_cut=500, t=0.02Tc)")
print("=" * 88)
print(f"{'lambda':>7} {'N':>7} {'Tc(K)':>11} {'w_max_cov(meV)':>15} {'2D(iw0)/kTc':>13} {'conv':>5}")
for lam in (0.5, 0.3, 0.2):
    for N in (400, 1600, 4000, 8000):
        spec = Spectrum.einstein(lam, 50.0)
        solver = Solver(spectrum=spec, mu_star=0.0, w_cut=500.0, n_matsubara=N)
        tc = solver.tc()
        t_low = 0.02 * tc
        wmax = 2 * np.pi * KB_MEV_PER_K * t_low * N
        wn, delta, z, conv = solver.solve_gap(t_low, max_iter=4000, mixing=0.2)
        r = 2 * delta[0] / (KB_MEV_PER_K * tc) if conv else float("nan")
        print(f"{lam:7.2f} {N:7d} {tc:11.5f} {wmax:15.1f} {r:13.4f} {str(conv):>5}", flush=True)
    print()

# =====================================================================
# 争议 2：AD 公式 Tc vs 自洽本征值 Tc，对判别性比值 R 的影响
# =====================================================================
print("=" * 88)
print("CHECK 2: R with AD-formula Tc  vs  R with self-consistent solver.tc()")
print("=" * 88)
LAM0, W0, WCUT, MU0 = 1.876, 106.9, 500.0, 0.13
DMU, DLAM_FRAC = 0.05, 0.15
N_USE = 1600


def obs(lam, mu, use_ad):
    spec = Spectrum.einstein(lam=lam, w0=W0)
    solver = Solver(spectrum=spec, mu_star=mu, w_cut=WCUT, n_matsubara=N_USE)
    tc = spec.allen_dynes(mu) if use_ad else solver.tc()
    wn, delta, z, conv = solver.solve_gap(0.03 * tc, mixing=0.25, max_iter=3000)
    return tc, 2.0 * delta[0] / (KB_MEV_PER_K * tc), conv


for use_ad in (True, False):
    tag = "AD-formula Tc (agent's choice)" if use_ad else "solver.tc() self-consistent"
    tcp, op, _ = obs(LAM0, MU0 + DMU, use_ad)
    tcm, om, _ = obs(LAM0, MU0 - DMU, use_ad)
    dmu = 0.5 * (op - om)
    tlp, olp, _ = obs(LAM0 * (1 + DLAM_FRAC), MU0, use_ad)
    tlm, olm, _ = obs(LAM0 * (1 - DLAM_FRAC), MU0, use_ad)
    dlam = 0.5 * (olp - olm)
    print(f"\n--- {tag} ---")
    print(f"  Tc(mu+)={tcp:8.3f}  Tc(mu-)={tcm:8.3f}   O(mu+)={op:.4f} O(mu-)={om:.4f}  dO_mu ={dmu:+.5f}")
    print(f"  Tc(l+) ={tlp:8.3f}  Tc(l-) ={tlm:8.3f}   O(l+) ={olp:.4f} O(l-) ={olm:.4f}  dO_lam={dlam:+.5f}")
    print(f"  R = {abs(dlam)/abs(dmu):.4f}")

print()
print("=" * 88)
print("CHECK 3: R 对 Dlambda 的线性依赖 —— R 只有在 Dlambda 是真实谐/非谐差时才有意义")
print("=" * 88)
for frac in (0.05, 0.15, 0.30):
    tcp, op, _ = obs(LAM0, MU0 + DMU, False)
    tcm, om, _ = obs(LAM0, MU0 - DMU, False)
    dmu = 0.5 * (op - om)
    tlp, olp, _ = obs(LAM0 * (1 + frac), MU0, False)
    tlm, olm, _ = obs(LAM0 * (1 - frac), MU0, False)
    dlam = 0.5 * (olp - olm)
    print(f"  Dlam/lam = {frac:.2f}  ==>  R = {abs(dlam)/abs(dmu):.4f}")
