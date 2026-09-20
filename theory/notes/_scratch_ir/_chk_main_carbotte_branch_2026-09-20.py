# -*- coding: utf-8 -*-
"""agent 的 Carbotte 交叉核对喂的是 AD 公式 Tc（= 它自己那处缺陷的分母）。
换成自洽 solver.tc() 的 Tc 后，Carbotte 的 R 是多少？符号是否与 Eliashberg 一致？
闭式，零求解器。主会话 2026-09-20。"""
import numpy as np
KB = 0.08617333262; WLOG = 106.9
def carb(tc):
    t = KB*tc/WLOG
    return 3.53*(1 + 12.5*t*t*np.log(1.0/(2*t)))

# Tc 四元组，两支（来自 _chk_main_gapratio_tcsrc_2026-09-20.py 实测）
BR = {"AD-formula Tc (agent)":      dict(mup=151.061, mum=201.905, lap=198.971, lam=149.176,
                                         eli_mu=+0.10614, eli_lam=+0.21371),
      "solver.tc() self-consistent":dict(mup=170.807, mum=209.333, lap=213.481, lam=160.660,
                                         eli_mu=-0.10342, eli_lam=+0.20767)}
for tag, d in BR.items():
    dmu  = 0.5*(carb(d["mup"]) - carb(d["mum"]))
    dlam = 0.5*(carb(d["lap"]) - carb(d["lam"]))
    R_c, R_e = abs(dlam)/abs(dmu), abs(d["eli_lam"])/abs(d["eli_mu"])
    agree = "AGREE" if np.sign(dmu) == np.sign(d["eli_mu"]) else "**SIGN CLASH**"
    print(f"--- {tag} ---")
    print(f"  Carbotte   dO_mu={dmu:+.5f}  dO_lam={dlam:+.5f}  R={R_c:.4f}")
    print(f"  Eliashberg dO_mu={d['eli_mu']:+.5f}  dO_lam={d['eli_lam']:+.5f}  R={R_e:.4f}")
    print(f"  sign(dO_mu): Carbotte vs Eliashberg -> {agree}")
    print(f"  R discrepancy Eliashberg/Carbotte = {R_e/R_c:.2f}x\n")
