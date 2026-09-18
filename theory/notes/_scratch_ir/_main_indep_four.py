"""Independent four-dataset Tc table (main session). Extends _main_indep_baseline.py
by a SEPARATE file per repo convention (never edit the baseline script).

Checks the headline of deficit_budget_h3s_2026-09-17.md: parameter-free (bar mu*)
Tc for H3S/D3S at 160/200 GPa vs experiment.  [numerical-evidence] float64, N=400.
"""
import os, sys, glob, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_indep_baseline import read_a2f, P_measure, Tc, MEV_K

ROOT = ("D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/"
        "Data-H3S-D3S-PH-EPW-Anharmonic")
EXPT = {"H3S-200GPa": 203.0, "D3S-200GPa": 150.0, "H3S-160GPa": 203.0, "D3S-160GPa": 150.0}

print("=== independent 4-dataset Tc, mu*=0.13, r=10, N=400 [numerical-evidence] ===")
print(f"{'dataset':>12} {'lambda':>9} {'sq<om2>':>9} {'Tc(0)':>9} {'Tc(.13)':>9} "
      f"{'deficit':>8} {'expt':>7} {'dev':>8}")
for tag in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa", "D3S-160GPa"):
    g = [x for x in glob.glob(f"{ROOT}/{tag}/EPW_anh*/*.a2f") if not x.endswith(".a2fv")]
    if not g:
        print(f"{tag:>12}   (no a2f file found)")
        continue
    w, a2f, cum = read_a2f(g[0])
    lam, p, om2 = P_measure(w, a2f)
    t0 = Tc(w, p, lam, om2, 0.0, 10.0, 400) * MEV_K
    t13 = Tc(w, p, lam, om2, 0.13, 10.0, 400) * MEV_K
    e = EXPT[tag]
    print(f"{tag:>12} {lam:>9.5f} {math.sqrt(om2):>9.3f} {t0:>9.2f} {t13:>9.2f} "
          f"{100*(1-t13/t0):>7.2f}% {e:>7.1f} {100*(t13/e-1):>7.2f}%", flush=True)
