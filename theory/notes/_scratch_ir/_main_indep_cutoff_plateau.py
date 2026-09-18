"""Downstream of the wscut finding: does the deficit-budget headline move when the
cutoff convention is made consistent with the authors', and does the Anderson-Morel
plateau claim survive on real alpha^2F data?

Three columns per dataset (all mu* referenced to the r=10 convention value 0.13):
  (a) r=10,     mu*=0.13                     -- the headline convention
  (b) r=r_auth, mu*=0.13                     -- naive cutoff swap (mu* NOT re-run)
  (c) r=r_auth, mu*=mu_AM(r_auth)            -- cutoff swap WITH Anderson-Morel running
          1/mu*(r) = 1/mu*(10) - ln(r/10)

If the AM plateau statement (am_plateau_2026-09-17.md) is right, (c) should sit close to
(a) -- that is the whole content of "r is a renormalization scale, not a second knob".
Column (b) measures how much of the move is just un-run mu*.

[numerical-evidence] float64, N=400, own operator implementation.
"""
import os, sys, glob, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_indep_baseline import read_a2f, P_measure, Tc, MEV_K

ROOT = ("D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/"
        "Data-H3S-D3S-PH-EPW-Anharmonic")
EXPT = {"H3S-200GPa": 203.0, "D3S-200GPa": 150.0, "H3S-160GPa": 203.0, "D3S-160GPa": 150.0}
WSCUT_MEV = 2500.0
MUS0, R0 = 0.13, 10.0

print("=== cutoff convention vs deficit-budget headline (mu*_ref=0.13 at r=10) ===")
print(f"{'dataset':>12} {'r_auth':>7} {'mu*_AM':>7} | {'(a) r=10':>9} {'(b) rA,mu.13':>12} "
      f"{'(c) rA,muAM':>11} | {'(c)-(a)':>8} {'expt':>6} {'dev(a)':>8} {'dev(c)':>8}")
for tag in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa", "D3S-160GPa"):
    g = [x for x in glob.glob(f"{ROOT}/{tag}/EPW_anh*/*.a2f") if not x.endswith(".a2fv")]
    if not g:
        print(f"{tag:>12}  (no a2f)"); continue
    w, a2f, cum = read_a2f(g[0])
    lam, p, om2 = P_measure(w, a2f)
    r_auth = WSCUT_MEV / math.sqrt(om2)
    mu_am = 1.0 / (1.0 / MUS0 - math.log(r_auth / R0))

    ta = Tc(w, p, lam, om2, MUS0, R0, 400) * MEV_K
    tb = Tc(w, p, lam, om2, MUS0, r_auth, 400) * MEV_K
    tc = Tc(w, p, lam, om2, mu_am, r_auth, 400) * MEV_K
    e = EXPT[tag]
    print(f"{tag:>12} {r_auth:>7.2f} {mu_am:>7.4f} | {ta:>9.2f} {tb:>12.2f} {tc:>11.2f} "
          f"| {tc-ta:>+8.2f} {e:>6.0f} {100*(ta/e-1):>+7.2f}% {100*(tc/e-1):>+7.2f}%",
          flush=True)
