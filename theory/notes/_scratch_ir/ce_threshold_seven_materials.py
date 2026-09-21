"""
CE threshold evaluation at 7 material anchor points (2026-09-17)

Task revision: Use the 7 actual material points (LaH10 4×, H3S, Nb, MgB2)
with ϖ_rms proxy for unmeasured points.

ϖ_rms proxy rule:
- LaH10 163 GPa: measured ω_rms=1202.6 K, Tc^E=244.4 K ⟹ ϖ_rms = 0.783
- Others: analytic lower bound ϖ_rms ≥ 0.8709/√λ (from Tc ≤ 0.1827√(λ⟨ω²⟩))
  μ*>0 suppresses Tc ⟹ raises ϖ ⟹ bound is a fortiori valid

Material parameters from pairing_duality_theorem §4.1:
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
import numpy as np

sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
from p4d_coulomb_extremal_threshold import CellExtremal

# Material anchor points: (name, λ, μ*, r_values)
MATERIALS = [
    ("LaH10 129 GPa", 3.62, 0.10),
    ("LaH10 163 GPa", 2.67, 0.10),
    ("LaH10 214 GPa", 2.06, 0.10),
    ("LaH10 264 GPa", 1.73, 0.10),
    ("H3S 200 GPa",   1.84, 0.16),
    ("Nb",            1.23, 0.13),
    ("MgB2",          0.62, 0.13),
]

# r values to test
R_VALUES = [3, 10]

N = int(os.environ.get("PN", "400"))

print(f"CE threshold at 7 material anchors (N={N})")
print("=" * 100)
print()
print("varpi_rms proxy: LaH10 163 GPa measured (0.783), others use floor 0.8709/sqrt(lambda)")
print("P = Einstein 2-atom equal-weight proxy")
print()

# Table header
header = (
    f"{'Material':<16} {'lam':>6} {'mu*':>6} {'r':>3} {'varpi':>7} {'Nc':>4} "
    f"{'mutilde':>7} {'mu*_CE':>9} {'margin':>7} {'verdict':>7} {'note':<10}"
)
print(header)
print("-" * 100)

results = []

for name, lam, mu_star in MATERIALS:
    # Compute ϖ_rms proxy
    if name == "LaH10 163 GPa":
        w_rms = 0.783  # measured
        note = "measured"
    else:
        w_rms = 0.8709 / np.sqrt(lam)  # analytic floor
        note = "floor"

    mu_tilde = mu_star / lam

    for r in R_VALUES:
        Nc = max(1, int(np.floor(r * w_rms)))

        if Nc > N // 2:
            print(f"{name:<16} {lam:6.2f} {mu_star:6.2f} {r:3d} {w_rms:7.4f} {Nc:4d} "
                  f"{mu_tilde:7.4f} {'SKIP':>9} {'N/A':>7} {'N/A':>7} Nc>N/2")
            continue

        # Einstein 2-atom proxy: equal weights at ω±
        # Use spread = 1.05 (near case, most conservative)
        T_fixed = 1.0
        omega_center = 2 * np.pi * w_rms * T_fixed
        spread = 1.05

        weights = np.array([0.5, 0.5])
        omegas = np.array([omega_center / spread, omega_center * spread])

        try:
            C = CellExtremal(weights, omegas, T_fixed, Nc, NN=N)
            mu_tilde_CE = C.thr

            margin = mu_tilde_CE / mu_tilde
            verdict = "PASS" if mu_tilde < mu_tilde_CE else "FAIL"

            results.append({
                'name': name,
                'lam': lam,
                'mu_star': mu_star,
                'r': r,
                'w_rms': w_rms,
                'Nc': Nc,
                'mu_tilde': mu_tilde,
                'mu_tilde_CE': mu_tilde_CE,
                'margin': margin,
                'verdict': verdict,
                'note': note,
            })

            print(f"{name:<16} {lam:6.2f} {mu_star:6.2f} {r:3d} {w_rms:7.4f} {Nc:4d} "
                  f"{mu_tilde:7.4f} {mu_tilde_CE:9.6f} {margin:7.3f} {verdict:>7} {note:<10}")

        except Exception as e:
            print(f"{name:<16} {lam:6.2f} {mu_star:6.2f} {r:3d} {w_rms:7.4f} {Nc:4d} "
                  f"{mu_tilde:7.4f} {'ERROR':>9} {'N/A':>7} {'N/A':>7} {str(e)[:20]}")

print()
print("=" * 100)
print("Summary:")
print(f"Total points: {len(MATERIALS) * len(R_VALUES)}")
print(f"Computed: {len(results)}")

if results:
    n_pass = sum(1 for r in results if r['verdict'] == 'PASS')
    n_fail = sum(1 for r in results if r['verdict'] == 'FAIL')
    print(f"PASS: {n_pass}, FAIL: {n_fail}")

    if n_pass > 0:
        margins_pass = [r['margin'] for r in results if r['verdict'] == 'PASS']
        print(f"Margin range (PASS): {min(margins_pass):.3f} – {max(margins_pass):.3f}")

    if n_fail > 0:
        print("\nFAIL details:")
        for r in results:
            if r['verdict'] == 'FAIL':
                print(f"  {r['name']}, r={r['r']}: mutilde={r['mu_tilde']:.4f} >= mu*_CE={r['mu_tilde_CE']:.6f}")

print()
print("Rigorous classification: numerical-evidence")
print("(float64 + N=400, varpi_rms proxy for 6/7 points)")
