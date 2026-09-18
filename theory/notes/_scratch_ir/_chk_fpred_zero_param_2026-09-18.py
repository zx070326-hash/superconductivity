"""Zero-free-parameter check: compute f_pred = Tc_P(mu*=0.13)/ceiling WITHOUT
ever touching Tc_exp, then separately compare to f_obs = Tc_exp/ceiling under
two anchors (203 K legacy vs 184.04 K corrected, per
theory/notes/materials/h3s_tc_exp_pressure_audit_2026-09-18.md).

Reuses the verified independent numerics from deficit_independent_check.py
(Simpson P-weights, explicit double loop, full eigh(M), secant-in-log-T root;
already cross-checked to 6 digits against the fast path, see that file's
header). This script does NOT re-derive the numerics -- it adds the Einstein
comparison stage (Tc_E(0): same lambda, same <omega^2>, single delta-function
mode) needed for delta_FIN, and keeps a strict separation between the
"prediction" stage (no Tc_exp anywhere) and the "compare to data" stage.

Severity: numerical-evidence (float64, N finite, alpha2F from external files).
"""
import os, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

import sys
sys.path.insert(0, os.path.dirname(__file__))
from deficit_independent_check import read_a2f, simpson_w, spec, lam1, Tc, MEV_K, ROOT

CINF = 0.1827262477890   # CERT-D upper end, rational-exact (see P4_g2_analytic.md / oneliner note)


def einstein_Tc(lam, om2, mus, r, N):
    """Tc for a single delta-function mode at omega_E = sqrt(om2), same lambda.
    Einstein Tc can exceed the [0.5,30] K bracket used by the generic Tc()
    (it's the hottest stage in the ladder, ceiling > Tc_E > Tc_P), so use a
    local secant with a wider bracket instead of importing Tc()'s fixed one."""
    sqE = math.sqrt(om2)
    W = np.array([sqE])
    p = np.array([1.0])
    g = lambda s: math.log(lam * lam1(W, p, math.exp(s), lam, mus, r, N, sqE))
    s0, s1 = math.log(0.5), math.log(600.0)
    g0, g1 = g(s0), g(s1)
    assert g0 > 0 > g1, (g0, g1)
    for _ in range(60):
        sm = 0.5*(s0+s1); gm = g(sm)
        if gm > 0: s0, g0 = sm, gm
        else:      s1, g1 = sm, gm
        if s1 - s0 < 1e-13:
            break
    return math.exp(0.5*(s0+s1))


ROWS = [
    ("H3S-200GPa", "EPW_anharmonic", "H3S-200 anharm"),
    ("H3S-200GPa", "EPW_harmonic",   "H3S-200 harm"),
    ("D3S-200GPa", "EPW_anharmonic", "D3S-200 anharm"),
    ("D3S-200GPa", "EPW_harmonic",   "D3S-200 harm"),
    ("H3S-160GPa", "EPW_anharmonic", "H3S-160 anharm"),
]

N = 250
R = 10.0

print("=== STAGE 1: PREDICTION (no Tc_exp used anywhere below this line until STAGE 2) ===")
print(f"{'row':>16} {'lam':>7} {'sq<om2>':>8} {'ceil[K]':>8} "
      f"{'Tc_E(0)':>8} {'dFIN':>6} {'Tc_P(0)':>8} {'dBRD':>6} {'Tc_P(.13)':>10} {'dCOU':>6} {'f_pred':>7}")

results = {}
for tag, sub, name in ROWS:
    gl = [x for x in glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f") if not x.endswith(".a2fv")]
    assert gl, (tag, sub)
    lam, om2, W, p, psum = spec(gl[0])
    assert abs(psum - 1.0) < 1e-9
    sq = math.sqrt(om2)
    ceil_K = CINF * math.sqrt(lam * om2) * MEV_K

    tE0 = einstein_Tc(lam, om2, 0.0, R, N) * MEV_K
    tP0 = Tc(W, p, lam, 0.0, R, N, sq) * MEV_K
    tP1 = Tc(W, p, lam, 0.13, R, N, sq) * MEV_K

    dFIN = (ceil_K - tE0) / ceil_K
    dBRD = (tE0 - tP0) / tE0
    dCOU = (tP0 - tP1) / tP0
    f_pred = tP1 / ceil_K
    f_pred_telescope = (1 - dFIN) * (1 - dBRD) * (1 - dCOU)
    assert abs(f_pred - f_pred_telescope) < 1e-9, (f_pred, f_pred_telescope)

    results[name] = dict(lam=lam, om2=om2, ceil=ceil_K, tE0=tE0, tP0=tP0, tP1=tP1,
                          dFIN=dFIN, dBRD=dBRD, dCOU=dCOU, f_pred=f_pred)
    print(f"{name:>16} {lam:7.4f} {sq:8.3f} {ceil_K:8.2f} "
          f"{tE0:8.2f} {100*dFIN:5.1f}% {tP0:8.2f} {100*dBRD:5.1f}% {tP1:10.3f} {100*dCOU:5.1f}% {f_pred:7.4f}")

print("\n=== STAGE 2: COMPARE f_pred to f_obs = Tc_exp/ceiling (Tc_exp now introduced) ===")
# Two anchors for H3S-200 GPa anharmonic specifically (per pressure audit note):
#   legacy: 203 K (155 GPa dome peak, pressure-mismatched with the 200 GPa calc)
#   corrected: 184.04 K (200 GPa, annealed, Fig.2c vector digitization, pressure-matched)
anchors = {
    "H3S-200 anharm": [("legacy 203K", 203.0), ("corrected 184.04K", 184.04)],
    "H3S-200 harm":   [("legacy 203K", 203.0), ("corrected 184.04K", 184.04)],
    "D3S-200 anharm": [("legacy ~150K", 150.0), ("corrected 145.5K (Fig.2c digitized @200GPa)", 145.5)],
    "D3S-200 harm":   [("legacy ~150K", 150.0)],
    "H3S-160 anharm": [("legacy 203K (pressure-mismatched)", 203.0)],
}
for name, anchs in anchors.items():
    r = results[name]
    for label, tc_exp in anchs:
        f_obs = tc_exp / r["ceil"]
        resid_pred_vs_obs = 100 * (r["f_pred"] - f_obs) / f_obs
        print(f"{name:>16} | anchor {label:>28} Tc_exp={tc_exp:6.1f} f_obs={f_obs:.4f} "
              f"f_pred={r['f_pred']:.4f}  (f_pred-f_obs)/f_obs={resid_pred_vs_obs:+.2f}%")

print("\n=== STAGE 3: orthogonal decomposition (H3S vs D3S are isotopes of ONE material, not two) ===")
# team-lead's point: H3S-200anharm and D3S-200anharm share the same crystal/electronic
# structure/alpha2F workflow, only the phonon mass scale differs -- their residuals are
# correlated. Decompose into (a) a common multiplicative offset shared by both, and
# (b) the isotope RATIO f(D)/f(H), theory vs experiment -- the latter is the genuinely
# independent test (systematic errors common to both isotopes cancel in the ratio).
fH, fD = results["H3S-200 anharm"]["f_pred"], results["D3S-200 anharm"]["f_pred"]
oH, oD = 184.04 / results["H3S-200 anharm"]["ceil"], 145.5 / results["D3S-200 anharm"]["ceil"]
ratio_pred = fD / fH
ratio_obs = oD / oH
common_offset = math.sqrt((fH/oH) * (fD/oD))  # geometric mean of the two (f_pred/f_obs) ratios
isotope_resid = ratio_pred / ratio_obs - 1.0
print(f"f_pred: H3S={fH:.4f} D3S={fD:.4f}  ratio_pred=f(D)/f(H)={ratio_pred:.4f}")
print(f"f_obs (corrected anchors): H3S={oH:.4f} D3S={oD:.4f}  ratio_obs=f(D)/f(H)={ratio_obs:.4f}")
print(f"common multiplicative offset (geometric mean of f_pred/f_obs) = {100*(common_offset-1):+.2f}%")
print(f"isotope-ratio residual (ratio_pred/ratio_obs - 1)            = {100*isotope_resid:+.2f}%")
# cross-check: (1+r_H)/(1+r_D) should equal ratio_obs/ratio_pred
rH = fH/oH - 1.0; rD = fD/oD - 1.0
print(f"cross-check: (1+rH)/(1+rD)={((1+rH)/(1+rD)):.5f}  vs  ratio_obs/ratio_pred={ratio_obs/ratio_pred:.5f}")

