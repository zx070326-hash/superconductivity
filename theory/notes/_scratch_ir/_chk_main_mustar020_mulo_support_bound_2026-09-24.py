"""mu_lo 若是支撑在 [0,Nv) 的 Rayleigh 商，必 <= lam1(P_Nv M P_Nv) <= lam1(P_600 M P_600)。
测这条链。Nv 取基线参数 40，另列 80/160/600 作梯子。"""
import os, importlib.util; os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
_d = os.path.join("theory","notes","_scratch_ir")
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(_d, n+".py"))
    m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
KT = _load("_chk_main_kato_temple_two_atom_prototype_2026-09-23")
print(f"{'lam':>8}{'w_lo':>10}{'Nc':>4} | " + "".join(f"{'P_'+str(n):>15}" for n in (40,80,160,600))
      + f"{'mu_lo':>15}{'判':>26}")
for lam, w, ml in [(2.0,1.100033,1.021190594),(3.2953,0.731708,1.020081705),
                   (4.0,0.640491,1.021883712),(6.0,0.495741,1.025180308)]:
    Nc = int(np.floor(10*w))
    vs = [KT.lam_k([lam],[w],0.20,Nc,n) for n in (40,80,160,600)]
    ok = ml <= vs[0] + 1e-12
    verdict = "OK (<=P_40)" if ok else ("** > P_40 且 > P_600 **" if ml > vs[-1] else "** > P_40, <=P_600 **")
    print(f"{lam:8}{w:10.6f}{Nc:4d} | " + "".join(f"{v:15.9f}" for v in vs)
          + f"{ml:15.9f}{verdict:>26}")
