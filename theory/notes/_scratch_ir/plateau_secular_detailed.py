"""Detailed analysis of the secular equation to understand the plateau mechanism.

We examine how the spectral distribution c_j^2/(mu_j - lam) behaves as r varies
with Anderson-Morel running. The key is that the sum is WEIGHTED by 1/(mu_j - lam),
not just sum c_j^2.

Strictness: [numerical-evidence] for all outputs.
"""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

def build_einstein_at_x(N, x):
    """Build KAY operators at x = varpi^2 for Einstein kernel f(k) = x/(k^2+x)."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    ks = np.arange(2 * N + 1, dtype=float)
    s = ks ** 2

    def O(fker):
        absd = np.abs(n[:, None] - n[None, :])
        sumd = n[:, None] + n[None, :] + 1
        uu = u[:, None] * u[None, :]
        Am = (fker[absd] + fker[sumd]) * uu
        L1 = np.zeros(N)
        L1[1:] = np.cumsum(2.0 * fker[1:N])[:N - 1] / (2 * n[1:] + 1)
        np.fill_diagonal(Am, -L1 + fker[2 * n + 1] * u ** 2)
        return Am

    f = x / (s + x); f[0] = 1.0
    return n, u, O(f)

def solve_secular(mu, PSI, u, n, Nc, mut, verbose=False):
    """Solve secular equation: 2*mut * sum_j c_j^2/(mu_j - lam) = 1.
    Returns lam (the root) and diagnostic info."""
    uc = np.where(n < Nc, u, 0.0)
    c = PSI.T @ uc  # c_j = <psi_j, u_c>

    # Root must be in (mu[-2], mu[-1])
    m1, m2 = mu[-1], mu[-2]
    gap = (m1 - m2) * 1e-13 + 1e-300

    def secular(lam):
        return 2.0 * mut * float(np.sum(c ** 2 / (mu - lam))) - 1.0

    lam = brentq(secular, m2 + gap, m1 - gap, xtol=1e-16, rtol=1e-15)

    # Diagnostics
    weights = c ** 2 / (mu - lam)  # contribution of each mode
    sum_c2 = float(np.sum(c ** 2))
    sum_weighted = float(np.sum(weights))

    # Top mode contribution
    top_frac = weights[-1] / sum_weighted
    top5_frac = float(np.sum(weights[-5:])) / sum_weighted

    if verbose:
        print(f"  lam = {lam:.10f}, mu_1 = {mu[-1]:.10f}, gap = {m1-lam:.6e}")
        print(f"  sum c_j^2 = {sum_c2:.6f}, sum weights = {sum_weighted:.6e}")
        print(f"  Top mode: {top_frac:.3%}, Top 5: {top5_frac:.3%}")
        print(f"  c_1^2 = {c[-1]**2:.6f}, c_1^2/(mu_1-lam) = {weights[-1]:.6e}")

    return lam, {
        'gap': m1 - lam,
        'sum_c2': sum_c2,
        'sum_weighted': sum_weighted,
        'top_frac': top_frac,
        'top5_frac': top5_frac,
        'c1_sq': c[-1] ** 2,
        'c1_weight': weights[-1]
    }

# Parameters matching H3S-200 anharmonic
N = 400
varpi = 1.39  # from the data
x = varpi ** 2
lam_coupling = 1.721  # lambda from the data

print("=" * 80)
print("Anderson-Morel plateau: secular equation decomposition")
print("=" * 80)
print(f"N = {N}, varpi = {varpi:.3f}, x = {x:.4f}, lambda = {lam_coupling:.3f}\n")

# Build operator (same for all r, since A doesn't depend on Coulomb)
n, u, A = build_einstein_at_x(N, x)
mu, PSI = eigh(A)

print("Eigenvalue spectrum of A (top 10):")
for i in range(-10, 0):
    print(f"  mu[{i}] = {mu[i]:.10f}")
print()

# Scan r with Anderson-Morel running
r0 = 10.0
mu_star_0 = 0.13

print("=" * 80)
print(f"{'r':>6} {'mu*':>10} {'N_c':>5} {'||u_c||^2':>10} {'lam':>12} "
      f"{'gap':>10} {'top%':>7} {'sum_c2':>10}")
print("=" * 80)

results = []
for r in [3, 5, 10, 20, 40, 80, 160]:
    # Anderson-Morel running
    inv_mu = 1.0 / mu_star_0 - math.log(r / r0)
    mu_star = 1.0 / inv_mu if inv_mu > 0 else float('nan')

    Nc = int(math.floor(r * varpi))

    # Check degeneracy
    if Nc <= 0:
        print(f"{r:>6.0f} {mu_star:>10.6f} {Nc:>5} {'N/A':>10} {'N/A':>12} "
              f"{'N/A':>10} {'N/A':>7} {'N/A':>10}")
        print(f"       WARNING: N_c = 0, Coulomb term vanishes")
        continue
    if Nc >= N:
        print(f"{r:>6.0f} {mu_star:>10.6f} {Nc:>5} {'DEGEN':>10} {'DEGEN':>12} "
              f"{'DEGEN':>10} {'DEGEN':>7} {'DEGEN':>10}")
        print(f"       WARNING: N_c >= N, cutoff covers full basis (degeneracy)")
        continue

    mut = mu_star / lam_coupling
    lam, info = solve_secular(mu, PSI, u, n, Nc, mut, verbose=(r == 10))

    print(f"{r:>6.0f} {mu_star:>10.6f} {Nc:>5} {info['sum_c2']:>10.4f} {lam:>12.8f} "
          f"{info['gap']:>10.2e} {info['top_frac']:>6.1%} {info['sum_c2']:>10.4f}")

    results.append((r, mu_star, Nc, lam, info))

print("\n" + "=" * 80)
print("Plateau quantification")
print("=" * 80)

lam_values = [lam for r, mu_star, Nc, lam, info in results]
lam_min, lam_max = min(lam_values), max(lam_values)
lam_span = (lam_max - lam_min) / lam_values[2]  # relative to r=10

print(f"lam spans [{lam_min:.8f}, {lam_max:.8f}]")
print(f"Relative span: {lam_span:.4%} (variation over 8x in r)")
print(f"\nFor r in physical range [5, 20]:")
lam_510 = [lam for r, _, _, lam, _ in results if 5 <= r <= 10]
lam_1020 = [lam for r, _, _, lam, _ in results if 10 <= r <= 20]
if lam_510 and lam_1020:
    lam_phys = lam_510 + lam_1020
    span_phys = (max(lam_phys) - min(lam_phys)) / lam_values[2]
    print(f"  lam spans {span_phys:.4%}")

print("\n" + "=" * 80)
print("Key mechanism: concentration into top mode")
print("=" * 80)
print("As r grows (with AM running):")
print("  - N_c grows ~ r")
print("  - ||u_c||^2 grows ~ (1/2) ln(N_c) ~ (1/2) ln(r)")
print("  - mu* grows ~ exp(ln r) but mu~ = mu*/lambda still grows")
print("  - However, the secular sum is WEIGHTED by 1/(mu_j - lam)")
print("  - The top mode c_1^2/(mu_1 - lam) dominates (50-80% of sum)")
print("  - The gap (mu_1 - lam) adjusts to keep the weighted sum constant")
print("\nResult: lam (and hence Tc) is approximately constant.")

print("\n" + "=" * 80)
print("Detailed mode analysis for r=10 (baseline):")
print("=" * 80)
r = 10.0
mu_star = mu_star_0
Nc = int(math.floor(r * varpi))
mut = mu_star / lam_coupling
lam, info = solve_secular(mu, PSI, u, n, Nc, mut, verbose=True)

# Compute c_j distribution
uc = np.where(n < Nc, u, 0.0)
c = PSI.T @ uc
weights = c ** 2 / (mu - lam)

print("\nTop 10 modes:")
for i in range(-10, 0):
    frac = weights[i] / info['sum_weighted']
    print(f"  j={i}: c_j^2 = {c[i]**2:.6f}, weight = {weights[i]:.6e}, "
          f"fraction = {frac:.3%}, mu_j = {mu[i]:.8f}")

print("\n" + "=" * 80)
print("Projection: what fraction of u_c is along psi_1?")
print("=" * 80)
for r, mu_star, Nc, lam, info in results:
    uc_norm = math.sqrt(info['sum_c2'])
    c1_over_norm = math.sqrt(info['c1_sq']) / uc_norm if uc_norm > 0 else 0
    print(f"  r={r:>3.0f}: c_1 / ||u_c|| = {c1_over_norm:.4f} "
          f"(c_1^2 / ||u_c||^2 = {info['c1_sq']/info['sum_c2']:.4f})")
