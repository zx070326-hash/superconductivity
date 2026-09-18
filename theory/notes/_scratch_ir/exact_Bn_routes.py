#!/usr/bin/env python3
"""
Explore three routes for analytic gap condition under exact B_n.
Route 1: Perron-Frobenius / nonnegative matrix theory
Route 2: Gershgorin + Cassini
Route 3: Piecewise (small vs large ϖ)
"""
import numpy as np
from mpmath import mp
mp.dps = 50

def make_kernel_lorentz(varpi, N):
    """f(k) = ϖ²/(ϖ²+k²), exact Lorentz."""
    return np.array([varpi**2 / (varpi**2 + k**2) for k in range(N)])

def make_g_shapes(N):
    """Generate several DOS weight profiles."""
    shapes = {}

    # Flat
    shapes['flat'] = np.ones(N)

    # van Hove step: dip at low n
    g_vH = np.ones(N)
    g_vH[:5] = 0.05
    shapes['vH_step'] = g_vH

    # Pseudogap step: enhancement at low n
    g_pg = np.ones(N)
    g_pg[:5] = 1.5
    shapes['pg_step'] = g_pg

    # Smooth van Hove (Gaussian dip)
    g_smooth_vH = 1 - 0.5 * np.exp(-((np.arange(N) - 3)**2) / 8.0)
    shapes['vH_smooth'] = g_smooth_vH

    return shapes

def compute_exact_Bn(g, lam, n):
    """
    Exact B_n from audit:
    B_n = Σ_{m≥0, m≠n} λ(|n−m|)·g_m  −  Σ_{i≥0} λ(n+1+i)·g_i
    """
    N = len(g)

    # First sum: m!=n
    sum1 = sum(lam[abs(n - m)] * g[m] for m in range(N) if m != n)

    # Second sum: negative frequency reflection
    sum2 = sum(lam[n + 1 + i] * g[i] for i in range(N) if (n + 1 + i) < len(lam))

    return sum1 - sum2

def build_L1_exact(g, lam, N):
    """Build L1[g,f] = diag(B_n/(2n+1))."""
    diag = np.array([compute_exact_Bn(g, lam, n) / (2*n + 1) for n in range(N)])
    return np.diag(diag)

def build_P(lam, N):
    """
    P[f] = L2 + L3, with u_n = 1/sqrt(2n+1).
    L2_nm = (1-δ_nm) f(|n-m|) u_n u_m
    L3_nm = f(n+m+1) u_n u_m
    """
    u = 1.0 / np.sqrt(2 * np.arange(N) + 1)
    P = np.zeros((N, N))

    for n in range(N):
        for m in range(N):
            if n != m:
                P[n, m] += lam[abs(n - m)] * u[n] * u[m]  # L2
            P[n, m] += lam[n + m + 1] * u[n] * u[m]  # L3

    return P

def compute_rho(lam, N_extend=1000):
    """
    ρ(f) = f(0) - f̂(π)
    f̂(π) = Σ_k f(k) (-1)^k
    """
    f_0 = lam[0]
    f_hat_pi = sum(lam[k] * ((-1)**k) for k in range(min(len(lam), N_extend)))
    return f_0 - f_hat_pi

print("=" * 80)
print("Route exploration: analytic gap condition under exact B_n")
print("=" * 80)

# Setup
N = 200
M = max(2 * N + 10, 500)  # Extended for lambda tail
varpi_test = [0.5, 1.0, 3.0, 5.0]

for varpi in varpi_test:
    print(f"\n{'=' * 80}")
    print(f"varpi = {varpi}")
    print("=" * 80)

    lam = make_kernel_lorentz(varpi, M)
    rho = compute_rho(lam, M)

    g_shapes = make_g_shapes(N)

    for g_name, g in g_shapes.items():
        # Build M = -L1 + P
        L1 = build_L1_exact(g, lam, N)
        P = build_P(lam, N)
        M_op = -L1 + P

        # Compute spectrum
        eigs = np.linalg.eigvalsh(M_op)
        eigs = np.sort(eigs)[::-1]  # Descending
        lam1, lam2 = eigs[0], eigs[1]
        gap_ratio = (lam1 - lam2) / lam1 if lam1 > 1e-10 else 0.0

        print(f"\n  {g_name:12s}: lam1={lam1:8.5f}, lam2={lam2:8.5f}, gap/lam1={gap_ratio:.5f}")

        # Route 1: Check if M + L1 is nonnegative (Perron structure)
        # M' = L1 + P should be ≥ -ρ·e0·e0ᵀ for flat case
        # For general g, check diagonal dominance
        M_prime = L1 + P

        # Diagonal dominance check
        diag_entries = np.diag(M_prime)
        offdiag_rowsum = np.array([np.sum(np.abs(M_prime[i, :])) - abs(M_prime[i, i])
                                   for i in range(N)])

        # Gershgorin radius
        gersh_rad = offdiag_rowsum
        gersh_center = diag_entries

        # Check if n=0 circle separates from rest
        circle_0_min = gersh_center[0] - gersh_rad[0]
        circle_0_max = gersh_center[0] + gersh_rad[0]

        rest_max = max(gersh_center[i] + gersh_rad[i] for i in range(1, min(20, N)))
        rest_min = min(gersh_center[i] - gersh_rad[i] for i in range(1, min(20, N)))

        separated = (circle_0_min > rest_max) or (circle_0_max < rest_min)

        print(f"    Route 1 (diagonal dominance): diag[0]={diag_entries[0]:.5f}, row_sum[0]={offdiag_rowsum[0]:.5f}")
        print(f"    Route 2 (Gershgorin): circle_0=[{circle_0_min:.5f}, {circle_0_max:.5f}], rest_max={rest_max:.5f}")
        print(f"              Separated: {separated}")

        # Check B_n structure
        B_n = np.array([compute_exact_Bn(g, lam, n) for n in range(min(10, N))])
        a_n = (B_n - rho * g[:len(B_n)]) / (2 * np.arange(len(B_n)) + 1)

        min_a = np.min(a_n)
        neg_count = np.sum(a_n < -1e-10)

        print(f"    a_n range: [{np.min(a_n):.6f}, {np.max(a_n):.6f}], neg_count={neg_count}")

print("\n" + "=" * 80)
print("Route 3: Small varpi analysis")
print("=" * 80)

# For small varpi, kernel decays slowly, L1 weights are small
for varpi in [0.3, 0.5, 1.0]:
    lam = make_kernel_lorentz(varpi, M)
    rho_val = compute_rho(lam, M)

    # Kernel tail behavior
    tail_10 = lam[10]
    tail_50 = lam[50] if len(lam) > 50 else 0

    print(f"\nvarpi={varpi}: rho(f)={rho_val:.6f}, f(10)={tail_10:.6f}, f(50)={tail_50:.6f}")
    print(f"  rho/f(0) = {rho_val/lam[0]:.6f}")

print("\nDone.")
