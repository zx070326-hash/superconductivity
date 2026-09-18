"""Independent re-run of the wscut cross-check (main session, own implementation).

Question: the -2.38% constant bias of our KAY linearized criterion vs the authors' full
Migdal-Eliashberg solution -- is it the criterion, or is it our r=10 cutoff convention?
The authors used a fixed Matsubara cutoff wscut = 2.5 eV and muc = 0.12.  Our r is defined
by omega_c = r * sqrt(<om^2>), so their cutoff corresponds to r_auth = 2500 meV / sqrt(<om^2>),
which is 16.7-25.7 depending on the dataset -- NOT 10.

Author ME Tc values are my own sqrt-extrapolation of their finite-Matsubara data
(verified earlier, m=2/3/5 windows agreeing to 0.4 K): 208.2 / 163.0 / 220.1 / 168.8 K.

[numerical-evidence] float64, N=400, own operator implementation (_main_indep_baseline).
"""
import os, sys, glob, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_indep_baseline import read_a2f, P_measure, Tc, MEV_K

ROOT = ("D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/"
        "Data-H3S-D3S-PH-EPW-Anharmonic")
WSCUT_MEV = 2500.0
AUTH_ME = {"H3S-200GPa": 208.2, "D3S-200GPa": 163.0,
           "H3S-160GPa": 220.1, "D3S-160GPa": 168.8}

print("=== independent wscut check: our KAY at authors' muc=0.12 ===")
print(f"{'dataset':>12} {'sq<om2>':>8} {'r_auth':>7} {'Tc(r=10)':>9} {'Tc(r_auth)':>11} "
      f"{'authorME':>9} {'dev@10':>8} {'dev@auth':>9}")
devs10, devsA = [], []
for tag in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa", "D3S-160GPa"):
    g = [x for x in glob.glob(f"{ROOT}/{tag}/EPW_anh*/*.a2f") if not x.endswith(".a2fv")]
    if not g:
        print(f"{tag:>12}   (no a2f)"); continue
    w, a2f, cum = read_a2f(g[0])
    lam, p, om2 = P_measure(w, a2f)
    sq = math.sqrt(om2)
    r_auth = WSCUT_MEV / sq
    t10 = Tc(w, p, lam, om2, 0.12, 10.0, 400) * MEV_K
    tra = Tc(w, p, lam, om2, 0.12, r_auth, 400) * MEV_K
    A = AUTH_ME[tag]
    d10, dA = 100 * (t10 / A - 1), 100 * (tra / A - 1)
    devs10.append(d10); devsA.append(dA)
    print(f"{tag:>12} {sq:>8.2f} {r_auth:>7.2f} {t10:>9.2f} {tra:>11.2f} {A:>9.1f} "
          f"{d10:>7.2f}% {dA:>8.2f}%", flush=True)

import statistics as st
print()
print(f"r=10 convention : mean dev {st.mean(devs10):+.3f}%,  std {st.pstdev(devs10):.3f}%")
print(f"authors' cutoff : mean dev {st.mean(devsA):+.3f}%,  std {st.pstdev(devsA):.3f}%")
