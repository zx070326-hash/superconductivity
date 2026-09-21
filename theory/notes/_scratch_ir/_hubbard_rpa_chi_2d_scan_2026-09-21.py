#!/usr/bin/env python3
"""
Two-orbital Hubbard RPA on 2D square lattice.

Try anisotropic hopping to break symmetry and induce negative χ_12.
"""
import numpy as np
from numpy.linalg import eigvalsh, inv

def fermi(E, T=1e-3):
    beta = 1.0 / T
    x = np.clip(beta * E, -50, 50)
    return 1.0 / (1.0 + np.exp(x))

def dispersion_2d(kx, ky, orbital, t1, t2, Delta):
    """
    2D square lattice with anisotropic hopping.

    Orbital 0: ε₁ = -2t1(cos kx + cos ky)
    Orbital 1: ε₂ = -2t2(cos kx + cos ky) + Δ
    """
    if orbital == 0:
        return -2.0 * t1 * (np.cos(kx) + np.cos(ky))
    else:
        return -2.0 * t2 * (np.cos(kx) + np.cos(ky)) + Delta

def lindhard_2d(qx, qy, omega, orb_a, orb_b, t1, t2, Delta, mu, Nk=40, T=1e-3):
    """Lindhard function in 2D."""
    k_grid = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
    kx_grid, ky_grid = np.meshgrid(k_grid, k_grid)

    eps_a = dispersion_2d(kx_grid, ky_grid, orb_a, t1, t2, Delta)
    eps_b = dispersion_2d(kx_grid + qx, ky_grid + qy, orb_b, t1, t2, Delta)

    f_a = fermi(eps_a - mu, T)
    f_b = fermi(eps_b - mu, T)

    denom = eps_a - eps_b + omega
    eta = 1e-4

    chi = np.sum((f_a - f_b) / (denom + 1j * eta))
    return chi.real / (Nk * Nk)

def compute_chi0_2d(qx, qy, omega, t1, t2, Delta, mu, Nk=40, T=1e-3):
    chi0 = np.zeros((2, 2))
    for a in range(2):
        for b in range(2):
            chi0[a, b] = lindhard_2d(qx, qy, omega, a, b, t1, t2, Delta, mu, Nk, T)
    return chi0

def compute_rpa_chi(chi0, U):
    I = np.eye(2)
    M = I - U * chi0
    try:
        chi_rpa = chi0 @ inv(M)
        return chi_rpa
    except np.linalg.LinAlgError:
        return None

def find_mu_2d(target_n, t1, t2, Delta, Nk=40, T=1e-3):
    from scipy.optimize import brentq

    def filling_error(mu):
        k_grid = np.linspace(-np.pi, np.pi, Nk, endpoint=False)
        kx_grid, ky_grid = np.meshgrid(k_grid, k_grid)

        n_total = 0.0
        for orb in range(2):
            eps = dispersion_2d(kx_grid, ky_grid, orb, t1, t2, Delta)
            n_total += np.sum(fermi(eps - mu, T))

        return n_total / (2 * Nk * Nk) - target_n

    mu = brentq(filling_error, -20.0, 20.0)
    return mu

def scan_2d():
    """Scan 2D parameter space with anisotropy."""

    # Try anisotropic hopping
    t1_vals = [1.0]  # Orbital 1 hopping
    t2_vals = [0.5, 0.75, 1.0, 1.5]  # Orbital 2 hopping (anisotropy)

    n_vals = np.linspace(0.5, 0.9, 9)
    U_vals = np.linspace(1.0, 5.0, 17)
    Delta_vals = np.linspace(0.3, 2.5, 12)

    # Try different q points
    q_points = [
        (np.pi, np.pi),    # AFM
        (np.pi, 0),        # Stripe
        (0.5*np.pi, 0.5*np.pi),  # Intermediate
    ]

    Nk = 50
    results = []

    print(f"2D scan with anisotropic hopping...")
    print(f"{'t2/t1':<6} {'q':<15} {'n':<6} {'U':<6} {'Δ':<6} {'χ_12':<10} {'λ_min':<10}")
    print("-" * 75)

    count = 0
    for t1 in t1_vals:
        for t2 in t2_vals:
            for (qx, qy) in q_points:
                q_label = f"({qx/np.pi:.2f}π,{qy/np.pi:.2f}π)"

                for n in n_vals:
                    for U in U_vals:
                        for Delta in Delta_vals:
                            count += 1
                            if count % 100 == 0:
                                print(f"Progress: {count}", flush=True)

                            try:
                                mu = find_mu_2d(n, t1, t2, Delta, Nk)
                            except:
                                continue

                            chi0 = compute_chi0_2d(qx, qy, 0.0, t1, t2, Delta, mu, Nk)
                            chi_rpa = compute_rpa_chi(chi0, U)

                            if chi_rpa is None:
                                continue

                            chi_12 = chi_rpa[0, 1]
                            eigs = eigvalsh(chi_rpa)
                            lambda_min = eigs[0]

                            if chi_12 < 0 and lambda_min > -1e-8:
                                results.append({
                                    't2_over_t1': t2/t1,
                                    'q': (qx, qy),
                                    'n': n,
                                    'U': U,
                                    'Delta': Delta,
                                    'chi_rpa': chi_rpa.copy(),
                                    'chi_12': chi_12,
                                    'lambda_min': lambda_min
                                })

                                print(f"{t2/t1:<6.2f} {q_label:<15} {n:<6.2f} {U:<6.2f} {Delta:<6.2f} "
                                      f"{chi_12:<10.4f} {lambda_min:<10.6f}")

    return results

def main():
    results = scan_2d()

    print("\n" + "=" * 75)
    print(f"2D SCAN COMPLETE: Found {len(results)} counterexamples")
    print("=" * 75)

    if results:
        print("\nTop 5 examples (most negative χ_12):")
        results_sorted = sorted(results, key=lambda x: x['chi_12'])

        for i, r in enumerate(results_sorted[:5], 1):
            print(f"\n--- Example {i} ---")
            print(f"t2/t1={r['t2_over_t1']:.2f}, q={r['q']}, n={r['n']:.3f}, U={r['U']:.3f}, Δ={r['Delta']:.3f}")
            chi = r['chi_rpa']
            print(f"χ^RPA = [[{chi[0,0]:9.5f}, {chi[0,1]:9.5f}]")
            print(f"         [{chi[1,0]:9.5f}, {chi[1,1]:9.5f}]]")
            print(f"Eigenvalues: {eigvalsh(chi)}")
    else:
        print("\nNo counterexamples found.")

if __name__ == '__main__':
    main()
