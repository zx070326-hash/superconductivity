import os, sys, math, glob
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
sys.path.insert(0, "/d/physics/superconductivity/theory/notes/_scratch_ir")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_budget import Tc_P
from h3s_fill_factor import moments, MEV_K, ROOT

WSCUT_MEV = 2500.0          # authors' wscut = 2.5 eV
AUTH = {"H3S-200GPa":208.2, "D3S-200GPa":163.0, "H3S-160GPa":220.1, "D3S-160GPa":168.8}
print("our KAY at the AUTHORS' muc=0.12 AND the authors' cutoff wscut=2.5 eV")
print(f"{'dataset':>14} {'sq<om2>':>8} {'r_auth':>7} | {'Tc(r=10)':>9} {'Tc(r_auth)':>11} "
      f"| {'author ME':>9} {'dev@r=10':>9} {'dev@r_auth':>11}")
for tag in ("H3S-200GPa","D3S-200GPa","H3S-160GPa","D3S-160GPa"):
    g=[x for x in glob.glob(f"{ROOT}/{tag}/EPW_anh*/*.a2f") if not x.endswith(".a2fv")]
    if not g: print(f"{tag:>14}  (no a2f)"); continue
    lam,om2,_,w,a2f = moments(g[0]); sq=math.sqrt(om2)
    r_auth = WSCUT_MEV/sq
    t10 = Tc_P(w,a2f,lam,0.12,10.0,400,sq)*MEV_K
    tra = Tc_P(w,a2f,lam,0.12,r_auth,400,sq)*MEV_K
    A=AUTH[tag]
    print(f"{tag:>14} {sq:>8.2f} {r_auth:>7.2f} | {t10:>9.2f} {tra:>11.2f} "
          f"| {A:>9.1f} {100*(t10/A-1):>8.2f}% {100*(tra/A-1):>10.2f}%", flush=True)
