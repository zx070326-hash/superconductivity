"""Derive f(mu*): explicit functional form for Coulomb correction to Tc.

Direct computation without expensive fits. Focus on:
1. mu_tilde scaling collapse
2. Weak coupling analytical check
3. Strong coupling decay law
4. N convergence at key points
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import math
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fill_factor_probe import fill, varpi_c

print("# f(mu*) derivation - Einstein spectrum, KAY framework")
print("# [numerical-evidence] float64, N=400")
print()

# Compact grid
lam_vals = [0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
mus_vals = [0.05, 0.10, 0.13, 0.15, 0.20]
r = 10.0
N = 400

print("# === Data Table ===")
print("# lambda   mu*   mu_tilde   Tc(0)      Tc(mu*)    D(%)    f=Tc/Tc0  varpi   Nc")

data = []
for lam in lam_vals:
    tc0, _ = fill(lam, 0.0, r, N)
    if tc0 is None:
        continue

    for mus in mus_vals:
        tc_mu, _ = fill(lam, mus, r, N)
        if tc_mu is None:
            continue

        vp_mu = varpi_c(lam, mus, r, N)
        Nc = int(math.floor(r * vp_mu))

        # Check validity: 0 < Nc < N
        if Nc <= 0 or Nc >= N:
            continue

        D = 1.0 - tc_mu / tc0
        f = tc_mu / tc0
        mu_tilde = mus / lam

        print(f"{lam:7.3f} {mus:6.3f} {mu_tilde:9.5f} {tc0:10.6f} {tc_mu:10.6f} "
              f"{D*100:7.3f} {f:9.6f} {vp_mu:7.4f} {Nc:4d}")

        data.append((lam, mus, mu_tilde, tc0, tc_mu, D, f, vp_mu, Nc))

print(f"\n# Valid points: {len(data)}")

# mu_tilde scaling collapse
print("\n# === mu_tilde Scaling Collapse ===")
print("# If D = F(mu*/lambda), then points with same mu_tilde should have same D")
print("# Sorted by mu_tilde:")
print("# mu_tilde   lambda   mu*      D(%)     f")

sorted_data = sorted(data, key=lambda x: x[2])
for lam, mus, mt, tc0, tc_mu, D, f, vp, Nc in sorted_data:
    print(f"{mt:9.5f} {lam:7.3f} {mus:6.3f} {D*100:8.4f} {f:9.6f}")

# Check weak coupling: compare to BCS/Morel-Anderson
print("\n# === Weak Coupling: Comparison to BCS/Morel-Anderson ===")
print("# BCS: f_BCS = exp(-mu*/(lam*(lam-mu*)))")
print("# Leading order: D ~ mu*/lam^2")
print("# lambda   mu*   D_obs(%)  D_BCS(%)  D_lead(%)  |obs-BCS|  |obs-lead|")

for lam, mus, mt, tc0, tc_mu, D, f, vp, Nc in data:
    if lam <= 1.0:  # Weak coupling
        f_BCS = math.exp(-mus / (lam * (lam - mus)))
        D_BCS = 1.0 - f_BCS
        D_lead = mus / (lam * lam)

        diff_BCS = abs(D - D_BCS) * 100
        diff_lead = abs(D - D_lead) * 100

        print(f"{lam:7.3f} {mus:6.3f} {D*100:8.4f}  {D_BCS*100:8.4f}  {D_lead*100:8.4f}  "
              f"{diff_BCS:8.4f}  {diff_lead:8.4f}")

# Strong coupling decay
print("\n# === Strong Coupling: Decay Law ===")
print("# Large lambda: mu_tilde -> 0, D -> 0")
print("# Check power law: log(D) vs log(mu_tilde)")
print("# If D ~ (mu*/lam)^alpha, then log(D) ~ alpha*log(mu_tilde)")
print("# lambda   mu*   mu_tilde   D(%)    log(D)    log(mu_t)  alpha_local")

for lam, mus, mt, tc0, tc_mu, D, f, vp, Nc in data:
    if lam >= 2.0 and D > 1e-8:  # Strong coupling
        log_D = math.log(D)
        log_mt = math.log(mt)

        # Local slope estimate (if we have another point nearby)
        alpha_local = None
        for lam2, mus2, mt2, tc02, tc_mu2, D2, f2, vp2, Nc2 in data:
            if abs(lam2 - lam) < 0.01 and abs(mus2 - mus) > 0.01 and D2 > 1e-8:
                log_D2 = math.log(D2)
                log_mt2 = math.log(mt2)
                if abs(log_mt2 - log_mt) > 0.01:
                    alpha_local = (log_D2 - log_D) / (log_mt2 - log_mt)
                    break

        alpha_str = f"{alpha_local:7.3f}" if alpha_local is not None else "    ---"
        print(f"{lam:7.3f} {mus:6.3f} {mt:9.5f} {D*100:7.4f} {log_D:9.4f} {log_mt:9.4f} {alpha_str}")

# Simple empirical formula
print("\n# === Empirical Formula Fit ===")
print("# Try: D = (mu*/lam) * (a0 + a1*lam + a2/lam)")
print("# [numerical-evidence]")

lams = np.array([x[0] for x in data])
muss = np.array([x[1] for x in data])
Ds = np.array([x[5] for x in data])
mu_tildes = muss / lams

# Fit
A = np.column_stack([mu_tildes, mu_tildes * lams, mu_tildes / lams])
coeffs, residuals, rank, s = np.linalg.lstsq(A, Ds, rcond=None)
D_fit = A @ coeffs
res = Ds - D_fit
rmse = np.sqrt(np.mean(res**2))

print(f"\nD = (mu*/lam) * ({coeffs[0]:.6f} + {coeffs[1]:.6f}*lam + {coeffs[2]:.6f}/lam)")
print(f"RMSE = {rmse*100:.4f}%")
print(f"Max |residual| = {np.max(np.abs(res))*100:.4f}%")

print("\n# Residuals:")
print("# lambda   mu*   D_obs(%)  D_fit(%)  residual(%)")
for i, (lam, mus, mt, tc0, tc_mu, D, f, vp, Nc) in enumerate(data):
    print(f"{lam:7.3f} {mus:6.3f} {D*100:8.4f}  {D_fit[i]*100:8.4f}  {res[i]*100:8.4f}")

# N convergence checks
print("\n# === N Convergence Checks ===")

test_points = [(1.5, 0.13), (3.0, 0.13), (0.7, 0.10)]

for lam_test, mus_test in test_points:
    print(f"\n# (lambda={lam_test}, mu*={mus_test})")
    print("# N      Tc(0)        Tc(mu*)      D(%)      Nc")

    for N_test in [200, 400, 800]:
        tc0, _ = fill(lam_test, 0.0, r, N_test)
        tc_mu, _ = fill(lam_test, mus_test, r, N_test)
        if tc0 is None or tc_mu is None:
            print(f"{N_test:5d}  FAIL")
            continue

        vp_mu = varpi_c(lam_test, mus_test, r, N_test)
        Nc = int(math.floor(r * vp_mu))
        D = 1.0 - tc_mu / tc0

        print(f"{N_test:5d}  {tc0:10.7f}  {tc_mu:10.7f}  {D*100:8.4f}  {Nc:4d}")

print("\n# === Summary ===")
print("# 1. mu_tilde scaling: check if D values cluster by mu*/lambda")
print("# 2. Weak coupling: compare to BCS exact and leading-order")
print("# 3. Strong coupling: characterize decay as lambda -> infinity")
print("# 4. Empirical formula: [numerical-evidence] fit for practical use")
