# -*- coding: utf-8 -*-
"""独立方案交叉验证：不是加大 N，而是减小 w_cut（两种方案的失效方式不同）。
若 (w_cut=500,N=4000) 与 (w_cut=200,N=4000) 一致，则 N 已收敛。
另附 Carbotte 强耦合闭式作第三方对照。主会话 2026-09-20。"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, r"D:\physics\superconductivity")
import numpy as np
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

BCS = 2*np.pi/np.exp(np.euler_gamma)
WE = 50.0
def carbotte(tc_k, wlog):
    t = KB_MEV_PER_K*tc_k/wlog
    return 3.53*(1 + 12.5*t*t*np.log(wlog/(2*KB_MEV_PER_K*tc_k)))

print(f"BCS = {BCS:.6f}   (Einstein w_E = {WE} meV, mu*=0, t_low=0.02Tc)")
print(f"{'lam':>5} {'w_cut':>6} {'N':>5} {'Tc(K)':>10} {'cover':>7} {'ratio':>8} {'Carbotte':>9} {'dev':>7}")
for lam in (0.5, 0.3):
    for wcut, N in ((500.0, 4000), (200.0, 2000), (200.0, 4000)):
        t0 = time.time()
        spec = Spectrum.einstein(lam, WE)
        s = Solver(spectrum=spec, mu_star=0.0, w_cut=wcut, n_matsubara=N)
        tc = s.tc(); tl = 0.02*tc
        cov = 2*np.pi*KB_MEV_PER_K*tl*N/wcut
        wn, d, z, c = s.solve_gap(tl, max_iter=4000, mixing=0.2)
        r = 2*d[0]/(KB_MEV_PER_K*tc); cb = carbotte(tc, WE)
        print(f"{lam:5.2f} {wcut:6.0f} {N:5d} {tc:10.5f} {cov:6.2f}x {r:8.4f} {cb:9.4f} {(r/cb-1)*100:+6.2f}%"
              f"  conv={c} [{time.time()-t0:.0f}s]", flush=True)
