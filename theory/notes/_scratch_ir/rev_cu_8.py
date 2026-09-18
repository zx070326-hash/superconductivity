"""rev_cu_8: (a) r=10 on (6,8] (the docstring claims [0.34,8], the code scans [0.34,6));
(b) how much headroom does the promised Tier-2 route for gamma have?
    mut_* is EXACTLY proportional to gamma, so a certified lower bound gam_cert
    gives mut_*^cert = mut_* * gam_cert/gam.  Compare the existing certified gap
    bounds (delta / g_alt of p4d_coulomb_certify, mu*=0 node) to the true gap."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
import numpy as np
from rev_cu_lib import RCell, uc, myH

NN = 400
def lam1(w, Nc, mt):
    A = myH(NN, w); u = uc(NN, Nc)
    return float(np.linalg.eigvalsh(A - 2*mt*np.outer(u, u))[-1])

print("(a) r=10 stepped path on [6,8.2] (the .py docstring H2 claims [0.34,8])")
for mt in (0.065, 0.13, 0.26):
    ws = np.arange(5.9, 8.2, 0.004)
    vals = np.array([lam1(w, int(np.floor(10*w)), mt) for w in ws])
    run = np.maximum.accumulate(vals)
    bad = np.where(vals < run - 1e-12)[0]
    print(f"   mut={mt}: {len(bad)} depressed points on [5.9,8.2]"
          + (f"  worst {(vals-run)[bad].min():.2e} at wbar={ws[bad[0]]:.3f}" if len(bad) else " (all jumps upward)"))

print("\n(b) headroom for a certified gamma")
print(f"   {'wbar':>7} {'true gam':>10} {'certified (log)':>18} {'ratio':>7} "
      f"{'needed ratio at anchor':>24}")
cert = {1.0: ("delta=0.657276 / g_alt=0.610122 (p4d_coulomb_certify, mu*=0)", 0.657276),
        5.0: ("delta=1.803759 / g_alt=1.134466 (probe_hi.log node w=5, mu*=0)", 1.803759)}
need = {13.5: (135, 0.26), 15.6: (46, 0.26), 51.2: (153, 0.325), 4.56: (45, 0.13/0.7),
        2.2: (22, 0.13), 1.0: (10, None), 5.0: (50, None)}
for w in (1.0, 2.2, 4.56, 5.0, 13.5, 15.6, 51.2):
    Nc, mut = need[w]
    C = RCell(w, Nc, NN)
    s = ""
    if w in cert:
        s = f"{cert[w][1]:.6f}  ratio {cert[w][1]/C.gam:.3f}"
    nr = f"{mut/C.thr:.3f}" if mut else "-"
    print(f"   {w:>7} {C.gam:>10.5f} {s:>40}   need gam_cert/gam >= {nr}")
print("\n   (mut_*^cert = mut_* * gam_cert/gam exactly, since d ln mut_*/d ln gam = 1.)")
