"""
MgB2 r=10 FAIL investigation (2026-09-17)

Critical finding: MgB2 at r=10 gives FAIL (μ̃=0.2097 ≥ μ̃*_CE=0.154356, margin=0.736)
This is the ONLY failure among 14 material×r points.

Question: Is this a genuine failure, or a one-sided direction artifact?

One-sided direction audit:
- μ̃*_CE is a THRESHOLD (not a bound)
- If six constants are UNDERESTIMATED → μ̃*_CE underestimated → PASS is conservative (safe)
- If six constants are OVERESTIMATED → μ̃*_CE overestimated → FAIL could be spurious

For MgB2: varpi_rms = 1.106 is the FLOOR (0.8709/sqrt(0.62))
- True varpi_rms could be HIGHER (Tc suppression from mu*>0)
- Higher varpi → higher Nc → potentially different six constants

Need to check:
1. Does team-lead's version also give FAIL at this point?
2. What happens if we increase N to check convergence direction?
3. Is the floor proxy too conservative for the low-lambda corner?
"""

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_v] = "1"

import sys
import numpy as np

sys.path.insert(0, "D:/physics/superconductivity/theory/derivations")
from p4d_coulomb_extremal_threshold import CellExtremal

# Focus on MgB2 r=10
lam = 0.62
mu_star = 0.13
mu_tilde = mu_star / lam  # 0.2097
r = 10

# Floor proxy
w_rms_floor = 0.8709 / np.sqrt(lam)
Nc_floor = max(1, int(np.floor(r * w_rms_floor)))

print(f"MgB2 r=10 investigation")
print("=" * 80)
print(f"lambda = {lam}, mu* = {mu_star}, mu_tilde = {mu_tilde:.4f}")
print(f"varpi_rms floor = {w_rms_floor:.4f}, Nc = {Nc_floor}")
print()

# Test N convergence: does threshold go up or down as N increases?
N_values = [200, 400, 800]

print("N convergence test:")
print(f"{'N':>6} {'mu*_CE':>10} {'margin':>8} {'delta':>10} {'verdict':>8}")
print("-" * 60)

T_fixed = 1.0
omega_center = 2 * np.pi * w_rms_floor * T_fixed
spread = 1.05
weights = np.array([0.5, 0.5])
omegas = np.array([omega_center / spread, omega_center * spread])

prev_thr = None
for N in N_values:
    C = CellExtremal(weights, omegas, T_fixed, Nc_floor, NN=N)
    mu_tilde_CE = C.thr
    margin = mu_tilde_CE / mu_tilde
    verdict = "PASS" if mu_tilde < mu_tilde_CE else "FAIL"

    if prev_thr is not None:
        delta = mu_tilde_CE - prev_thr
        delta_str = f"{delta:+.6f}"
    else:
        delta_str = "---"

    print(f"{N:6d} {mu_tilde_CE:10.6f} {margin:8.4f} {delta_str:>10} {verdict:>8}")
    prev_thr = mu_tilde_CE

print()
print("Interpretation:")
print("- If threshold INCREASES with N: N=400 gives conservative FAIL (real threshold higher)")
print("- If threshold DECREASES with N: N=400 gives optimistic FAIL (real threshold lower)")
print()

# Also test varpi sensitivity: what if true varpi is 20% higher than floor?
print("varpi_rms sensitivity test (N=400, varpi = floor × factor):")
print(f"{'factor':>8} {'varpi':>8} {'Nc':>4} {'mu*_CE':>10} {'margin':>8} {'verdict':>8}")
print("-" * 70)

for factor in [1.0, 1.1, 1.2, 1.3, 1.5]:
    w_test = w_rms_floor * factor
    Nc_test = max(1, int(np.floor(r * w_test)))

    if Nc_test > 200:  # N=400, so Nc > N/2 would fail
        print(f"{factor:8.2f} {w_test:8.4f} {Nc_test:4d} {'SKIP':>10} {'N/A':>8} {'Nc>N/2':>8}")
        continue

    omega_center_test = 2 * np.pi * w_test * T_fixed
    omegas_test = np.array([omega_center_test / spread, omega_center_test * spread])

    C_test = CellExtremal(weights, omegas_test, T_fixed, Nc_test, NN=400)
    mu_tilde_CE_test = C_test.thr
    margin_test = mu_tilde_CE_test / mu_tilde
    verdict_test = "PASS" if mu_tilde < mu_tilde_CE_test else "FAIL"

    print(f"{factor:8.2f} {w_test:8.4f} {Nc_test:4d} {mu_tilde_CE_test:10.6f} {margin_test:8.4f} {verdict_test:>8}")

print()
print("Conclusion:")
print("MgB2 r=10 is the ONLY failure among 14 points.")
print("If N convergence shows threshold rising, FAIL is conservative (safe).")
print("If varpi floor is too tight, true material (higher varpi) could PASS.")
