# -*- coding: utf-8 -*-
"""CHECK 2: 分母 Tc 用 AD 闭式 vs 用自洽 solver.tc()，对 R 的影响。主会话 2026-09-20。"""
import os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
sys.path.insert(0, r"D:\physics\superconductivity")
import numpy as np
from theory.solvers.eliashberg_iso import Spectrum, Solver, KB_MEV_PER_K

LAM0, W0, WCUT, MU0, N = 1.876, 106.9, 500.0, 0.13, 1600
DMU, DLAM = 0.05, 0.15

def obs(lam, mu, use_ad):
    spec = Spectrum.einstein(lam=lam, w0=W0)
    s = Solver(spectrum=spec, mu_star=mu, w_cut=WCUT, n_matsubara=N)
    tc = spec.allen_dynes(mu) if use_ad else s.tc()
    wn, d, z, c = s.solve_gap(0.03*tc, mixing=0.25, max_iter=3000)
    return tc, 2.0*d[0]/(KB_MEV_PER_K*tc)

for use_ad in (True, False):
    tag = "AD-formula Tc (agent)" if use_ad else "solver.tc() self-consistent"
    t0 = time.time()
    tcp, op = obs(LAM0, MU0+DMU, use_ad); tcm, om = obs(LAM0, MU0-DMU, use_ad)
    tlp, olp = obs(LAM0*(1+DLAM), MU0, use_ad); tlm, olm = obs(LAM0*(1-DLAM), MU0, use_ad)
    dmu, dlam = 0.5*(op-om), 0.5*(olp-olm)
    print(f"\n--- {tag} ---", flush=True)
    print(f"  Tc(mu+)={tcp:8.3f} Tc(mu-)={tcm:8.3f}  O={op:.4f}/{om:.4f}  dO_mu ={dmu:+.5f}", flush=True)
    print(f"  Tc(l+) ={tlp:8.3f} Tc(l-) ={tlm:8.3f}  O={olp:.4f}/{olm:.4f}  dO_lam={dlam:+.5f}", flush=True)
    print(f"  R = {abs(dlam)/abs(dmu):.4f}      [{time.time()-t0:.0f}s]", flush=True)
