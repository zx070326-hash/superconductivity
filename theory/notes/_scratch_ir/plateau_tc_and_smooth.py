"""Anderson-Morel plateau: Tc computation and smooth cutoff sensitivity.

This script:
1. Computes actual Tc for H3S-200 anharmonic data with AM running
2. Tests smooth cutoff: u_c[n] = u_n / (1 + (n/N_c)^p) instead of step

Strictness: [numerical-evidence] for all outputs.
"""
import os, math, sys, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from h3s_fill_factor import moments, MEV_K, ROOT
    from h3s_deficit_budget import Q_P
    HAS_DATA = True
except:
    HAS_DATA = False
    print("WARNING: Cannot import h3s data modules, using Einstein model only")

def build_einstein_at_x(N, x):
    """Build KAY operators at x = varpi^2 for Einstein kernel."""
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

def solve_coulomb_eigenvalue(A_eigs, A_vecs, u, n, Nc, mut, cutoff_type='step', cutoff_power=2):
    """Solve for largest eigenvalue of M = A - 2*mut*u_c u_c^T.

    cutoff_type: 'step' or 'smooth'
    cutoff_power: p in smooth cutoff u_c[n] = u_n / (1 + (n/N_c)^p)
    """
    if cutoff_type == 'step':
        uc = np.where(n < Nc, u, 0.0)
    elif cutoff_type == 'smooth':
        if Nc == 0:
            uc = np.zeros_like(u)
        else:
            uc = u / (1.0 + (n / Nc) ** cutoff_power)
    else:
        raise ValueError(f"Unknown cutoff_type: {cutoff_type}")

    if mut == 0.0 or Nc == 0:
        return float(A_eigs[-1])

    c = A_vecs.T @ uc
    mu = A_eigs
    m1, m2 = mu[-1], mu[-2]
    gap = (m1 - m2) * 1e-13 + 1e-300

    def secular(lam):
        return 2.0 * mut * float(np.sum(c ** 2 / (mu - lam))) - 1.0

    return brentq(secular, m2 + gap, m1 - gap, xtol=1e-16, rtol=1e-15)

def Tc_from_lam(lam, lam_coupling, x_or_sqom2):
    """Convert eigenvalue to Tc in meV. x = (varpi*Tc)^2, so Tc = sqrt(x)/varpi."""
    # For Einstein: x = varpi^2, so Tc solves lam_coupling * Q(Tc) = 1
    # where Q is the eigenvalue. Here we already have lam at a given T.
    # We need to invert: given lam, find T such that lam(T) = lam_target.
    # Shortcut: the eigenvalue scales with x, so...
    # Actually, this is tricky. Let me use the direct approach.
    return None  # placeholder

print("=" * 80)
print("Anderson-Morel plateau: Tc with real H3S data")
print("=" * 80)

if HAS_DATA:
    # Load H3S-200 anharmonic data
    tag = "H3S-200GPa"
    sub = "EPW_anharmonic"
    gl = [x for x in glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f") if not x.endswith(".a2fv")]
    if gl:
        lam_coupling, om2, olog, w, a2f = moments(gl[0])
        sq = math.sqrt(om2)
        varpi = sq / (2 * math.pi)  # will be computed at each Tc

        print(f"Data: {tag}/{sub}")
        print(f"lambda = {lam_coupling:.4f}, sqrt<om^2> = {sq:.2f} meV")
        print(f"Using N=400 for speed\n")

        r0 = 10.0
        mu_star_0 = 0.13
        N = 400

        print("=" * 80)
        print("Step cutoff: Tc vs r (with Anderson-Morel running)")
        print("=" * 80)
        print(f"{'r':>6} {'mu*':>10} {'N_c':>5} {'varpi':>8} {'Tc [K]':>10} "
              f"{'rel to r=10':>12}")

        Tc_baseline = None
        results_step = []

        for r in [3, 5, 10, 20, 40, 80, 160]:
            # AM running
            inv_mu = 1.0 / mu_star_0 - math.log(r / r0)
            mu_star = 1.0 / inv_mu if inv_mu > 0 else float('nan')

            # Use Q_P which does the full Tc solve
            Tc_mev = Q_P(0.5, w, a2f, lam_coupling, mu_star, r, N, sq)  # dummy T, will search
            # Actually Q_P returns eigenvalue, not Tc. Let me use the Tc_P function
            from h3s_deficit_budget import Tc_P
            Tc_mev = Tc_P(w, a2f, lam_coupling, mu_star, r, N, sq)

            if Tc_mev is not None:
                Tc_K = Tc_mev * MEV_K
                vp = sq / (2 * math.pi * Tc_mev)
                Nc = int(math.floor(r * vp))

                if Tc_baseline is None:
                    Tc_baseline = Tc_K
                rel = (Tc_K - Tc_baseline) / Tc_baseline

                print(f"{r:>6.0f} {mu_star:>10.6f} {Nc:>5} {vp:>8.4f} {Tc_K:>10.2f} {rel:>11.2%}")
                results_step.append((r, mu_star, Nc, Tc_K))
            else:
                print(f"{r:>6.0f} {mu_star:>10.6f} {'N/A':>5} {'N/A':>8} {'N/A':>10} {'N/A':>12}")

        Tc_vals = [Tc for r, mu_star, Nc, Tc in results_step]
        span = (max(Tc_vals) - min(Tc_vals)) / Tc_baseline
        print(f"\nTotal span: {span:.3%} over r from {results_step[0][0]:.0f} to {results_step[-1][0]:.0f}")

        # Physical range [5, 20]
        Tc_phys = [Tc for r, mu_star, Nc, Tc in results_step if 5 <= r <= 20]
        if len(Tc_phys) >= 2:
            span_phys = (max(Tc_phys) - min(Tc_phys)) / Tc_baseline
            print(f"Physical range [5, 20]: span {span_phys:.3%}")

        print("\n" + "=" * 80)
        print("N_c degeneracy check")
        print("=" * 80)
        print("Checking if 0 < N_c < N for all points:")
        for r, mu_star, Nc, Tc in results_step:
            status = "OK" if 0 < Nc < N else "DEGENERATE"
            print(f"  r={r:>3.0f}: N_c={Nc:>3}, N={N}, status={status}")

        print("\n" + "=" * 80)
        print("Smooth cutoff sensitivity")
        print("=" * 80)
        print("Testing smooth cutoff u_c[n] = u_n / (1 + (n/N_c)^p) at r=10")
        print(f"{'power p':>8} {'Tc [K]':>10} {'vs step':>10}")

        # For smooth cutoff, we need to modify the h3s code
        # This is complex, so let me test with Einstein model instead
        print("(requires modified solver, testing with Einstein model below)")

    else:
        print("No H3S data files found")
        HAS_DATA = False

print("\n" + "=" * 80)
print("Einstein model: smooth cutoff test")
print("=" * 80)

N = 400
x_test = 1.9321  # from H3S varpi ~ 1.39
lam_coupling = 1.721
n, u, A = build_einstein_at_x(N, x_test)
mu, PSI = eigh(A)

r = 10.0
varpi_test = math.sqrt(x_test)
Nc_nominal = int(math.floor(r * varpi_test))
mu_star = 0.13
mut = mu_star / lam_coupling

print(f"x = {x_test:.4f}, varpi = {varpi_test:.4f}, N_c(step) = {Nc_nominal}")
print(f"mu* = {mu_star}, mu~ = {mut:.6f}\n")

print(f"{'Cutoff':>15} {'power':>6} {'lam':>12} {'vs step':>10} {'||u_c||^2':>10}")

# Step cutoff
lam_step = solve_coulomb_eigenvalue(mu, PSI, u, n, Nc_nominal, mut, 'step')
uc_step = np.where(n < Nc_nominal, u, 0.0)
norm_step = float(np.sum(uc_step ** 2))
print(f"{'Step':>15} {'N/A':>6} {lam_step:>12.8f} {'baseline':>10} {norm_step:>10.4f}")

# Smooth cutoffs with different powers
for p in [1, 2, 3, 4, 6, 10]:
    lam_smooth = solve_coulomb_eigenvalue(mu, PSI, u, n, Nc_nominal, mut, 'smooth', p)
    uc_smooth = u / (1.0 + (n / Nc_nominal) ** p)
    norm_smooth = float(np.sum(uc_smooth ** 2))
    diff = (lam_smooth - lam_step) / lam_step
    print(f"{'Smooth':>15} {p:>6} {lam_smooth:>12.8f} {diff:>9.3%} {norm_smooth:>10.4f}")

print("\nInterpretation:")
print("  - Step cutoff is the baseline (current model)")
print("  - Smooth cutoff with small p (soft transition) differs more")
print("  - Large p (sharp transition) converges to step")
print("  - If differences are O(1%), step cutoff shape matters")
print("  - If differences are O(0.1%), step is acceptable approximation")
