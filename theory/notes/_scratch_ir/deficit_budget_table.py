"""Deficit budget + inverse mu*: given Tc_exp, what mu* does the exact operator demand?
numerical-evidence (float64). Uses the validated fast path; cross-checked against
deficit_independent_check.py at mu*=0 and 0.13."""
import os, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import brentq
import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_budget import Tc_P
from h3s_fill_factor import moments, MEV_K, ROOT
from fill_factor_probe import CINF, fill

ROWS = [("H3S-200GPa", "EPW_anharmonic", 203.0), ("H3S-200GPa", "EPW_harmonic", 203.0),
        ("D3S-200GPa", "EPW_anharmonic", 150.0), ("D3S-200GPa", "EPW_harmonic", 150.0),
        ("H3S-160GPa", "EPW_anharmonic", 203.0)]

print("deficit budget: ceiling -> Einstein -> actual shape -> Coulomb -> experiment")
print(f"{'dataset':>26} {'ceil':>7} {'Tc_E':>7} {'dFIN':>6} {'Tc_P':>7} {'dBRD':>6} "
      f"{'Tc.13':>7} {'dCOU':>6} | {'Tc_exp':>7} {'resid':>7} {'mu*_fit':>8}")
out = {}
for tag, sub, texp in ROWS:
    gl = [x for x in glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f") if not x.endswith(".a2fv")]
    lam, om2, olog, w, a2f = moments(gl[0])
    sq = math.sqrt(om2)
    ceil_K = CINF*math.sqrt(lam*om2)*MEV_K
    tE = Tc_P(np.array([sq]), np.array([lam*sq/2.0]), lam, 0.0, 10.0, 400, sq)*MEV_K
    tP0 = Tc_P(w, a2f, lam, 0.0, 10.0, 400, sq)*MEV_K
    tP1 = Tc_P(w, a2f, lam, 0.13, 10.0, 400, sq)*MEV_K
    f = lambda m: Tc_P(w, a2f, lam, m, 10.0, 400, sq)*MEV_K - texp
    try:
        mfit = brentq(f, 1e-4, 0.60, xtol=1e-6)
    except ValueError:
        mfit = float("nan")
    print(f"{tag+'/'+sub:>26} {ceil_K:>7.1f} {tE:>7.1f} {(1-tE/ceil_K)*100:>5.1f}% "
          f"{tP0:>7.1f} {(1-tP0/tE)*100:>5.1f}% {tP1:>7.1f} {(1-tP1/tP0)*100:>5.1f}% | "
          f"{texp:>7.1f} {(tP1/texp-1)*100:>+6.1f}% {mfit:>8.4f}", flush=True)
    out[(tag, sub)] = (tP1, mfit)

a_pred = math.log(out[("H3S-200GPa","EPW_anharmonic")][0]/out[("D3S-200GPa","EPW_anharmonic")][0])/math.log(2)
a_harm = math.log(out[("H3S-200GPa","EPW_harmonic")][0]/out[("D3S-200GPa","EPW_harmonic")][0])/math.log(2)
print(f"\nisotope exponent alpha_H = ln(Tc_H/Tc_D)/ln 2, mu*=0.13 both:")
print(f"  anharmonic a2F : {a_pred:.4f}")
print(f"  harmonic   a2F : {a_harm:.4f}")
print(f"  experiment 2-pt: {math.log(203.0/150.0)/math.log(2):.4f}   (authors quote ~0.3, interval fit)")
print("\nstep-cutoff artifact guard: Nc = floor(r*varpi) hits 0 for large lambda")
for L in (10.0, 30.0, 100.0):
    vp = 1.0/(2*math.pi*fill(L,0.0)[0]); print(f"  lam={L:>5.0f} varpi_c={vp:.4f} Nc(r=10)={int(10*vp)}"
          + ("   <-- Coulomb term VANISHES, deficit 0% is an artifact" if int(10*vp)==0 else ""))
