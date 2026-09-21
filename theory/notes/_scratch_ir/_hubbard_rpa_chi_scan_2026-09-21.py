#!/usr/bin/env python3
"""
Two-orbital Hubbard RPA spin susceptibility negative off-diagonal scan.

Goal: Find parameters where χ^RPA_12 < 0 but χ^RPA ⪰ 0 (PSD).

Model:
  H = Σ_k ε_α(k) c†_kα c_kα + U Σ_i,α n_iα↑ n_iα↓

RPA: χ^RPA = χ₀ [1 − U χ₀]^{-1}

Dispersion: ε₁(k) = −2t cos(k), ε₂(k) = −2t cos(k) + Δ
"""
import numpy as np
from numpy.linalg import eigvalsh, inv

def fermi(E, T=1e-3):
    """Fermi function at small T (mimics T→0)."""
    beta = 1.0 / T
    x = beta * E
    # Clip to avoid overflow
    x = np.clip(x, -50, 50)
    return 1.0 / (1.0 + np.exp(x))

def dispersion(k, orbital, t=1.0, Delta=0.0):
    """1D tight-binding dispersion."""
    if orbital == 0:
        return -2.0 * t * np.cos(k)
    else:
        return -2.0 * t * np.cos(k) + Delta

def lindhard_element(q, omega, orbital_a, orbital_b, t, Delta, mu, Nk=200, T=1e-3):
    """
    Lindhard function χ₀_{αβ}(q,ω) at T≈0.

    χ₀_{αβ}(q,ω) = (1/Nk) Σ_k [f(ε_α(k)) − f(ε_β(k+q))] / [ε_α(k) − ε_β(k+q) + ω + iη]

    For static (ω=0) and real part, we use principal value approximation.
    """
    k_grid = np.linspace(-np.pi, np.pi, Nk, endpoint=False)

    eps_a_k = np.array([dispersion(k, orbital_a, t, Delta) for k in k_grid])
    eps_b_kq = np.array([dispersion(k + q, orbital_b, t, Delta) for k in k_grid])

    f_a = fermi(eps_a_k - mu, T)
    f_b = fermi(eps_b_kq - mu, T)

    denom = eps_a_k - eps_b_kq + omega
    eta = 1e-4  # Small imaginary part for principal value

    # Real part (static limit ω=0)
    chi = np.sum((f_a - f_b) / (denom + 1j * eta))

    return chi.real / Nk

def compute_chi0_matrix(q, omega, t, Delta, mu, Nk=200, T=1e-3):
    """Compute 2×2 Lindhard matrix χ₀."""
    chi0 = np.zeros((2, 2))
    for a in range(2):
        for b in range(2):
            chi0[a, b] = lindhard_element(q, omega, a, b, t, Delta, mu, Nk, T)
    return chi0

def compute_rpa_chi(chi0, U):
    """
    Compute RPA susceptibility χ^RPA = χ₀ [1 − U χ₀]^{-1}.
    """
    I = np.eye(2)
    M = I - U * chi0
    try:
        chi_rpa = chi0 @ inv(M)
        return chi_rpa
    except np.linalg.LinAlgError:
        return None

def find_mu_for_filling(target_n, t, Delta, Nk=200, T=1e-3):
    """Find chemical potential for target filling n per orbital."""
    from scipy.optimize import brentq

    def filling_error(mu):
        k_grid = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
        n_total = 0.0
        for orb in range(2):
            eps = np.array([dispersion(k, orb, t, Delta) for k in k_grid])
            n_total += np.sum(fermi(eps - mu, T))
        avg_filling = n_total / (2 * Nk)  # Average over 2 orbitals
        return avg_filling - target_n

    # Search range for mu
    mu = brentq(filling_error, -10.0 * t, 10.0 * t)
    return mu

def scan_parameters():
    """
    Scan parameter space for negative χ_12 with PSD χ^RPA.
    """
    t = 1.0  # Energy unit

    # Parameter grids
    n_vals = np.linspace(0.5, 0.9, 9)
    U_vals = np.linspace(1.0, 4.0, 13) * t
    Delta_vals = np.linspace(0.5, 2.0, 7) * t

    q = np.pi  # Momentum transfer
    omega = 0.0  # Static limit
    Nk = 300

    results = []
    count = 0
    total = len(n_vals) * len(U_vals) * len(Delta_vals)

    print(f"Scanning {total} parameter combinations...")
    print(f"{'n':<6} {'U/t':<6} {'Δ/t':<6} {'χ_11':<10} {'χ_12':<10} {'χ_22':<10} {'λ_min':<10} {'λ_max':<10}")
    print("-" * 70)

    for n in n_vals:
        for U in U_vals:
            for Delta in Delta_vals:
                count += 1
                if count % 50 == 0:
                    print(f"Progress: {count}/{total}", flush=True)

                # Find chemical potential
                try:
                    mu = find_mu_for_filling(n, t, Delta, Nk)
                except:
                    continue

                # Compute χ₀
                chi0 = compute_chi0_matrix(q, omega, t, Delta, mu, Nk)

                # Compute RPA χ
                chi_rpa = compute_rpa_chi(chi0, U)
                if chi_rpa is None:
                    continue

                # Check condition
                chi_12 = chi_rpa[0, 1]
                eigs = eigvalsh(chi_rpa)
                lambda_min, lambda_max = eigs[0], eigs[1]

                # Store if χ_12 < 0 and matrix is PSD
                if chi_12 < 0 and lambda_min > -1e-8:  # Small tolerance for numerical noise
                    results.append({
                        'n': n,
                        'U_over_t': U / t,
                        'Delta_over_t': Delta / t,
                        'chi_11': chi_rpa[0, 0],
                        'chi_12': chi_12,
                        'chi_22': chi_rpa[1, 1],
                        'lambda_min': lambda_min,
                        'lambda_max': lambda_max,
                        'chi_rpa': chi_rpa.copy()
                    })

                    print(f"{n:<6.2f} {U/t:<6.2f} {Delta/t:<6.2f} {chi_rpa[0,0]:<10.4f} {chi_12:<10.4f} "
                          f"{chi_rpa[1,1]:<10.4f} {lambda_min:<10.6f} {lambda_max:<10.4f}")

    return results

def main():
    results = scan_parameters()

    print("\n" + "=" * 70)
    print(f"SCAN COMPLETE: Found {len(results)} counterexamples")
    print("=" * 70)

    if results:
        print("\nTop 10 counterexamples (most negative χ_12):")
        results_sorted = sorted(results, key=lambda x: x['chi_12'])

        for i, r in enumerate(results_sorted[:10], 1):
            print(f"\n--- Example {i} ---")
            print(f"Parameters: n={r['n']:.3f}, U/t={r['U_over_t']:.3f}, Δ/t={r['Delta_over_t']:.3f}")
            print(f"χ^RPA matrix:")
            print(f"  [[{r['chi_11']:9.5f}, {r['chi_12']:9.5f}]")
            print(f"   [{r['chi_12']:9.5f}, {r['chi_22']:9.5f}]]")
            print(f"Eigenvalues: λ_min={r['lambda_min']:.6f}, λ_max={r['lambda_max']:.6f}")
            print(f"χ_12 / sqrt(χ_11 χ_22) = {r['chi_12'] / np.sqrt(r['chi_11'] * r['chi_22']):.4f}")
    else:
        print("\nNo counterexamples found in the scanned parameter range.")
        print("\nSuggestions:")
        print("1. Try different dispersions (e.g., square lattice nearest-neighbor)")
        print("2. Extend parameter ranges (larger U, different Δ)")
        print("3. Consider frequency-dependent χ(q, ω≠0)")
        print("4. Check if Mermin-Wagner / fluctuation-dissipation constraints forbid this")

if __name__ == '__main__':
    main()
