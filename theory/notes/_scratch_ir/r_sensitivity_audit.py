"""SELF-AUDIT of my own claim "zero free parameters except mu*".
r = omega_c/sqrt<om^2> was CHOSEN = 10. If the headline 199 K moves much with r,
then r is a second fit knob and the claim must be weakened.
Also: how much of [0.10,0.15] does mu*_fit sweep as r ranges over defensible values?
numerical-evidence (float64)."""
import os, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np, sys
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_budget import Tc_P
from h3s_fill_factor import moments, MEV_K, ROOT

ROWS = [("H3S-200GPa", "EPW_anharmonic", 203.0), ("D3S-200GPa", "EPW_anharmonic", 150.0),
        ("H3S-160GPa", "EPW_anharmonic", 203.0)]
RS = (3.0, 5.0, 10.0, 20.0, 40.0, 80.0)

print("Tc_P(mu*=0.13) [K] vs r      (r=10 is what I used for the headline)")
print(f"{'dataset':>26} " + " ".join(f"{'r='+str(int(r)):>8}" for r in RS) + "   spread")
tab = {}
for tag, sub, texp in ROWS:
    gl = [x for x in glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f") if not x.endswith(".a2fv")]
    lam, om2, olog, w, a2f = moments(gl[0])
    sq = math.sqrt(om2)
    vals = []
    for r in RS:
        vals.append(Tc_P(w, a2f, lam, 0.13, r, 400, sq)*MEV_K)
    tab[tag] = (lam, om2, sq, w, a2f, vals)
    sp = (max(vals)-min(vals))/vals[RS.index(10.0)]*100
    print(f"{tag+'/'+sub:>26} " + " ".join(f"{v:>8.1f}" for v in vals) + f"   {sp:>5.1f}%", flush=True)

print("\nNc = floor(r*varpi) at the solution (integer -> step artifacts possible)")
for tag, sub, texp in ROWS:
    lam, om2, sq, w, a2f, vals = tab[tag]
    line = f"{tag:>26} "
    for r, v in zip(RS, vals):
        vp = sq/(2*math.pi*(v/MEV_K))
        line += f" Nc={int(r*vp):<5}"
    print(line + f"  (varpi={sq/(2*math.pi*(vals[2]/MEV_K)):.3f})")

print("\nmu*_fit needed to hit Tc_exp, vs r   -- does it stay in [0.10,0.15]?")
print(f"{'dataset':>26} " + " ".join(f"{'r='+str(int(r)):>8}" for r in RS))
for tag, sub, texp in ROWS:
    lam, om2, sq, w, a2f, _ = tab[tag]
    line = f"{tag+'/'+sub:>26} "
    for r in RS:
        f = lambda m: Tc_P(w, a2f, lam, m, r, 400, sq)*MEV_K - texp
        try:
            line += f" {brentq(f, 1e-4, 0.8, xtol=1e-6):>7.4f}"
        except ValueError:
            line += f" {'n/a':>7}"
    print(line, flush=True)

print("\nisotope exponent alpha_H vs r (H3S-200 / D3S-200, mu*=0.13)")
for i, r in enumerate(RS):
    a = math.log(tab["H3S-200GPa"][5][i]/tab["D3S-200GPa"][5][i])/math.log(2)
    print(f"  r={r:>5.0f}: alpha_H = {a:.4f}")

# ---- consistency: mu* RUNS with omega_c (Anderson-Morel), so scanning r at FIXED
# mu* is inconsistent. AM: 1/mu*(wc) = 1/mu*(wc0) - ln(wc/wc0)  =>  anchored at r=10.
print("\n=== r-scan WITH Anderson-Morel running of mu* (anchored mu*(r=10)=0.13) ===")
def mu_run(r, mu0=0.13, r0=10.0):
    inv = 1.0/mu0 - math.log(r/r0)
    return 1.0/inv if inv > 0 else float("nan")
print(f"{'r':>5} {'mu*(r)':>8} | " + " ".join(f"{t[:8]:>9}" for t,_,_ in ROWS))
for r in RS:
    m = mu_run(r)
    line = f"{r:>5.0f} {m:>8.4f} | "
    for tag, sub, texp in ROWS:
        lam, om2, sq, w, a2f, _ = tab[tag]
        line += f" {Tc_P(w, a2f, lam, m, r, 400, sq)*MEV_K:>8.1f}"
    print(line, flush=True)
