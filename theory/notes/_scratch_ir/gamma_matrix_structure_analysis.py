"""
Analytical structure of KAY's G(gamma) matrix in the γ→∞ limit.

For large γ, the power-law kernels 1/k^γ decay extremely rapidly:
  - G2_{nm} = u_n u_m / |n-m|^γ  →  only |n-m|=1 survives (tri-diagonal)
  - G3_{nm} = u_n u_m / (n+m+1)^γ  →  diagonal-dominated

This script computes the effective sparsity pattern and dominant contributions
for γ=1,10,50 to verify the tri-diagonal + diagonal-correction structure.
"""
import numpy as np
from mpmath import mp

mp.dps = 20


def build_G_components(gamma_val, N):
    """Build G1, G2, G3 separately as numpy arrays, returning their contributions."""
    gamma = gamma_val
    dim = N + 1

    def u(n):
        return 1.0 / np.sqrt(2*n + 1)

    u_arr = np.array([u(n) for n in range(dim)])

    # G1: diagonal only, a[n] = (1/(2n+1)) * sum_{k=1}^n 2/k^gamma
    G1 = np.zeros((dim, dim))
    running = 0.0
    for n in range(1, dim):
        running += 2.0 / n**gamma
        G1[n, n] = running / (2*n + 1)

    # G2: Toeplitz, off-diagonal only
    G2 = np.zeros((dim, dim))
    for n in range(dim):
        for m in range(dim):
            if n != m:
                G2[n, m] = u_arr[n] * u_arr[m] / abs(n - m)**gamma

    # G3: Hankel, all entries
    G3 = np.zeros((dim, dim))
    for n in range(dim):
        for m in range(dim):
            G3[n, m] = u_arr[n] * u_arr[m] / (n + m + 1)**gamma

    G = -G1 + G2 + G3

    return G, G1, G2, G3, u_arr


def sparsity_analysis(G, threshold=1e-6):
    """Count entries above threshold (relative to max |G_ij|)."""
    max_abs = np.max(np.abs(G))
    significant = np.abs(G) >= threshold * max_abs
    return np.sum(significant), significant


def tridiagonal_approximation(G2, G3, u_arr, gamma_val):
    """
    Build the γ→∞ leading-order approximation:
      G2_tri: only |n-m|=1 terms
      G3_diag: only n=m terms
    """
    dim = len(u_arr)
    G2_tri = np.zeros((dim, dim))
    for n in range(dim):
        if n > 0:
            G2_tri[n, n-1] = u_arr[n] * u_arr[n-1] / 1**gamma_val
        if n < dim - 1:
            G2_tri[n, n+1] = u_arr[n] * u_arr[n+1] / 1**gamma_val

    G3_diag = np.zeros((dim, dim))
    for n in range(dim):
        G3_diag[n, n] = u_arr[n]**2 / (2*n + 1)**gamma_val

    return G2_tri, G3_diag


if __name__ == '__main__':
    N = 50  # Small N for readability
    gammas = [1.0, 5.0, 10.0, 50.0]

    print("γ-model matrix structure analysis (N={})".format(N))
    print("=" * 70)

    for gamma_val in gammas:
        G, G1, G2, G3, u_arr = build_G_components(gamma_val, N)

        # Sparsity at 1e-6 threshold
        n_sig_G2, _ = sparsity_analysis(G2, threshold=1e-6)
        n_sig_G3, _ = sparsity_analysis(G3, threshold=1e-6)

        # Tri-diagonal approximation
        G2_tri, G3_diag = tridiagonal_approximation(G2, G3, u_arr, gamma_val)

        # Error of tri-diagonal + diagonal approximation
        G_approx = -G1 + G2_tri + G3_diag
        err_fro = np.linalg.norm(G - G_approx, 'fro')
        err_2 = np.linalg.norm(G - G_approx, 2)

        # Eigenvalues
        evals_full = np.linalg.eigvalsh(G)
        evals_approx = np.linalg.eigvalsh(G_approx)
        lam1_full = evals_full[-1]
        lam1_approx = evals_approx[-1]

        print(f"\nγ = {gamma_val}")
        print(f"  G2 significant entries (>1e-6 * max): {n_sig_G2} / {(N+1)**2}")
        print(f"  G3 significant entries (>1e-6 * max): {n_sig_G3} / {(N+1)**2}")
        print(f"  Tri-diagonal + diagonal approximation:")
        print(f"    ||G - G_approx||_F = {err_fro:.6e}")
        print(f"    ||G - G_approx||_2 = {err_2:.6e}")
        print(f"    λ₁(G_full) = {lam1_full:.8f}")
        print(f"    λ₁(G_approx) = {lam1_approx:.8f}")
        print(f"    relative error = {abs(lam1_full - lam1_approx)/lam1_full * 100:.3f}%")

    print("\n" + "=" * 70)
    print("CONCLUSION:")
    print("  For γ≳10, G2 is effectively tri-diagonal (|n-m|=1 dominates)")
    print("  and G3 is diagonal-dominated (n=m dominates).")
    print("  This suggests an analytical route via tri-diagonal eigenvalue theory.")
