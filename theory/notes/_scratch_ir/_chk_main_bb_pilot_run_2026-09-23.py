# -*- coding: utf-8 -*-
"""试跑：先在余量最大的 λ 上跑一次分支定界，量盒子数与耗时，再决定要不要上更紧的 λ。"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import importlib.util
import sys
import numpy as np

_d = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_d, name + ".py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


BB = _load("_chk_main_box_branch_bound_certificate_2026-09-23")
KT = BB.KT

r, mus, lam_star = 10.0, 0.13, 3.2953
M2 = lam_star*KT.einstein_root(lam_star, mus, r, 600)**2
lam = float(sys.argv[1]) if len(sys.argv) > 1 else 3.20
maxbox = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
wE = float(np.sqrt(M2/lam)); Nc = int(np.floor(r*wE))
lE = KT.lam_k([lam], [wE], mus, Nc, 600)
print(f"λ={lam}  ω̄₂={wE:.9f}  N_c={Nc}  λ₁(Einstein,float)={lE:.9f}  余量 {1-lE:.3e}", flush=True)
R = BB.certify(lam, wE, mus, Nc, 1.0, maxbox=maxbox, log_every=50)
print(f"盒子 {R['nleaf']}（R₀ {R['n_R0']} / KT {R['n_KT']}）  分裂 {R['nsplit']}  "
      f"对称跳过 {R['nskip']}  深度 {R['maxdepth']}  最劣上界 {R['worst']:.9f}  "
      f"γ_max {R['gamma_max']:.4f}  {R['secs']:.0f}s  失败 {len(R['fails'])}")
for f in R["fails"][:10]:
    print("  FAIL", f)
