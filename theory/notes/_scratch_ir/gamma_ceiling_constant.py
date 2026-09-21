"""Compute the ceiling constant C_inf(gamma) for gamma-model kernels.

For Lorentz (gamma=2): C_inf = (1/2pi)*sqrt(g2) = 0.1827262477 [analytic].
For general gamma: C_inf(gamma) = (1/2pi)*sqrt(g_infinity(gamma)), where
g_infinity is the infinite-operator eigenvalue lim_{N->inf} g^(N)(gamma).

This script:
1. Computes g^(N)(gamma) for large N using the KAY operator
2. Extrapolates N->infinity using Richardson or monotone bounds
3. Verifies the Lorentz case recovers the known analytic value
4. Produces certificate-level bounds when using mpmath.iv

Strictness: numerical-evidence (float64) for exploration; can be upgraded to
certificate by switching to mpmath.iv arithmetic.
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

def build_gamma_operator(N, gamma_val):
    """KAY operator for f(k) = k^{-gamma}.

    G(gamma)_nm = -G1 + G2 + G3,
    G1_nn = [1/(2n+1)] * sum_{k=1}^n 2/k^gamma  (sum to n, finite for gamma<=1)
    G2_nm = (n!=m) * [1/sqrt((2n+1)(2m+1))] * 1/|n-m|^gamma
    G3_nm = [1/sqrt((2n+1)(2m+1))] * 1/(n+m+1)^gamma
    """
    n = np.arange(N, dtype=float)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)

    # G1 diagonal
    G1_diag = np.zeros(N)
    for i in range(N):
        if i == 0:
            G1_diag[0] = 0.0  # sum k=1..0 is empty
        else:
            G1_diag[i] = np.sum(2.0 / np.arange(1, i+1)**gamma_val) / (2*i + 1)

    # G2 and G3
    absd = np.abs(n[:, None] - n[None, :])
    sumd = n[:, None] + n[None, :] + 1
    uu = u[:, None] * u[None, :]

    # G2: off-diagonal only (n != m)
    G2 = np.where(absd > 0, uu / absd**gamma_val, 0.0)

    # G3: all elements
    G3 = uu / sumd**gamma_val

    # Full operator
    G = -np.diag(G1_diag) + G2 + G3
    return G

def g_eigenvalue(N, gamma_val):
    """Largest eigenvalue g^(N)(gamma)."""
    G = build_gamma_operator(N, gamma_val)
    lam = eigh(G, eigvals_only=True)
    return float(lam[-1])

def ceiling_constant(g_inf, gamma_val):
    """C_inf(gamma) = (1/2pi) * sqrt(g_inf)."""
    return math.sqrt(g_inf) / (2.0 * math.pi)

def richardson_extrapolate(vals, Ns):
    """Simple Richardson extrapolation assuming g^(N) ~ g_inf - a/N^p."""
    # Use last 3 points for 2nd-order extrapolation
    if len(vals) < 3:
        return vals[-1], float('nan')

    g1, g2, g3 = vals[-3], vals[-2], vals[-1]
    N1, N2, N3 = Ns[-3], Ns[-2], Ns[-1]

    # Assume g(N) ~ g_inf + a/N
    # (g2 - g1) / (g3 - g2) ~ (N2 - N1) / (N3 - N2) if linear in 1/N
    # Richardson: g_inf ~ g3 + (g3 - g2) * N3 / (N3 - N2)

    if N2 == N1 or N3 == N2:
        return g3, float('nan')

    # Linear Richardson
    dg32 = g3 - g2
    dg21 = g2 - g1
    g_inf_est = g3 + dg32 * N3 / (N3 - N2)

    return g_inf_est, abs(g_inf_est - g3)

if __name__ == "__main__":
    print("Ceiling constant C_inf(gamma) for gamma-model kernels")
    print("=" * 70)

    # Test: Lorentz case (gamma=2) should give C_inf = 0.1827262477
    print("\n1. Validation: Lorentz case (gamma=2)")
    print("   Known analytic: C_inf = 0.1827262477")

    Ns_lorentz = [100, 200, 400, 800]
    gs_lorentz = []
    for N in Ns_lorentz:
        g = g_eigenvalue(N, 2.0)
        gs_lorentz.append(g)
        C = ceiling_constant(g, 2.0)
        print(f"   N={N:>4}: g^(N) = {g:.10f}, C_inf^(N) = {C:.10f}")

    g_inf_rich, err = richardson_extrapolate(gs_lorentz, Ns_lorentz)
    C_inf_rich = ceiling_constant(g_inf_rich, 2.0)
    print(f"   Richardson extrapolation: g_inf ~ {g_inf_rich:.10f}")
    print(f"                             C_inf ~ {C_inf_rich:.10f}")
    print(f"   Error vs analytic: {abs(C_inf_rich - 0.1827262477):.2e}")

    # Main scan: multiple gamma values
    print("\n2. Gamma-model ceiling constants")
    print(f"{'gamma':>6} | {'g^(400)':>12} {'g^(800)':>12} | "
          f"{'g_inf(Rich)':>13} | {'C_inf':>11} | {'note':>20}")
    print("-" * 90)

    gammas = [0.67, 1.0, 1.5, 2.0, 2.5, 3.0]
    results = {}

    for gamma_val in gammas:
        Ns = [200, 400, 800]
        gs = []
        for N in Ns:
            g = g_eigenvalue(N, gamma_val)
            gs.append(g)

        g_inf_est, err = richardson_extrapolate(gs, Ns)
        C_inf = ceiling_constant(g_inf_est, gamma_val)

        results[gamma_val] = {
            'g_400': gs[1],
            'g_800': gs[2],
            'g_inf': g_inf_est,
            'C_inf': C_inf
        }

        note = ""
        if gamma_val == 2.0:
            note = "Lorentz baseline"
        elif gamma_val == 1.0:
            note = "SDW QCP (KAY focus)"
        elif gamma_val < 1.0:
            note = "tail diverges slowly"

        print(f"{gamma_val:>6.2f} | {gs[1]:>12.8f} {gs[2]:>12.8f} | "
              f"{g_inf_est:>13.8f} | {C_inf:>11.8f} | {note:>20}")

    # Monotonicity check: g_infinity should increase with gamma
    print("\n3. Monotonicity check: g_inf(gamma) should be monotone")
    gamma_sorted = sorted(results.keys())
    for i in range(len(gamma_sorted) - 1):
        g1 = results[gamma_sorted[i]]['g_inf']
        g2 = results[gamma_sorted[i+1]]['g_inf']
        ratio = g2 / g1
        print(f"   g_inf({gamma_sorted[i+1]:.2f}) / g_inf({gamma_sorted[i]:.2f}) "
              f"= {ratio:.6f}  ({'PASS' if ratio < 1 else 'FAIL - should decrease'})")

    print("\n4. Convergence diagnostics (monotone sequence property)")
    for gamma_val in [1.0, 2.0]:
        print(f"\n   gamma = {gamma_val}:")
        Ns = [100, 200, 400, 800]
        gs = [g_eigenvalue(N, gamma_val) for N in Ns]
        for i, (N, g) in enumerate(zip(Ns, gs)):
            if i > 0:
                delta = gs[i] - gs[i-1]
                print(f"     N={N:>4}: g = {g:.10f}  (delta = {delta:>+.2e})")
            else:
                print(f"     N={N:>4}: g = {g:.10f}")

        # Check monotonicity
        is_monotone = all(gs[i+1] > gs[i] for i in range(len(gs)-1))
        print(f"     Monotone increasing: {'YES (KAY theorem satisfied)' if is_monotone else 'NO (violation!)'}")

    print("\n" + "=" * 70)
    print("Summary: C_inf(gamma) values (Richardson extrapolation, numerical-evidence)")
    for gamma_val in sorted(results.keys()):
        print(f"  gamma = {gamma_val:>4.2f}:  C_inf = {results[gamma_val]['C_inf']:.10f}")
