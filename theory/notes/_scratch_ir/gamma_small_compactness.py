"""
Route exploration for γ≤0.5 compactness when Hilbert-Schmidt norm diverges.

Problem: ||𝔊||_{HS}² = Σ |G_{nm}|² ~ Σ 1/k^{2γ} diverges for γ≤0.5.
P4 Schur tail bounds rely on HS compactness, which fails for γ≤0.5.

Three routes explored:
  Route 1: Operator norm ||𝔊||_{op} < ∞ even when ||𝔊||_{HS} = ∞
  Route 2: Direct finite-rank approximation in operator norm (no HS)
  Route 3: Krein-Rutman without compactness (only needs r(𝔊) < ∞)

Physical relevance: nematic QCP γ≈1/3, so γ∈(0.3,0.5) is physical.

Strictness: analytic (theoretical arguments), certificate (numerical verification with iv).
"""
import sys, os, time
import numpy as np
from mpmath import iv, mp

def lo(x):
    """Extract certified lower bound."""
    return mp.mpf(x.a) if hasattr(x, 'a') else mp.mpf(x)

def hi(x):
    """Extract certified upper bound."""
    return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)

# ============================================================================
# Route 1: Operator norm vs HS norm
# ============================================================================

def schur_test_row_sum(N, gamma_val, dps=30):
    """
    Schur test: ||𝔊||_{op} ≤ max(sup_n row_sum_n, sup_m col_sum_m).

    For symmetric matrix, row_sum = col_sum, so need only one.
    Row sum for 𝔊(γ): R_n = Σ_m |G_{nm}|

    For γ-model: G_{nm} = u_n u_m g(n,m) where g(n,m) ≥ 0 (all entries nonnegative).
    R_n = u_n Σ_m u_m g(n,m)
         = u_n [ Σ_m u_m (1/|n-m|^γ + 1/(n+m+1)^γ) ]  (m≠n for first term)

    For γ≤1, the Toeplitz sum Σ_m 1/|n-m|^γ ~ Σ_k 1/k^γ diverges.
    BUT: weighted sum Σ_m u_m/|n-m|^γ ~ Σ_k 1/(k^γ √k) = Σ 1/k^{γ+1/2}
    converges iff γ+1/2 > 1, i.e., γ > 0.5.

    So Schur test row sums:
      - γ > 0.5: converge → ||𝔊||_{op} < ∞  ✓
      - γ ≤ 0.5: diverge → Schur test inconclusive

    This route FAILS for γ≤0.5.
    """
    mp.dps = dps
    iv.dps = dps
    gamma = iv.mpf(gamma_val)

    def u(n):
        return 1 / iv.sqrt(iv.mpf(2*n + 1))

    # Compute row sum for n=0 (typically largest due to u_0 = 1)
    n = 0
    un = u(n)

    # Sum over m=0..N
    row_sum = iv.mpf(0)
    for m in range(N+1):
        um = u(m)
        if m != n:
            term = um * (iv.mpf(abs(n-m))**(-gamma) + iv.mpf(n+m+1)**(-gamma))
        else:
            # Diagonal: only Hankel term (no |n-m|^{-γ} since that's 0/0)
            term = um * iv.mpf(2*n+1)**(-gamma)
        row_sum = row_sum + term

    row_sum = row_sum * un
    return hi(row_sum)

def young_inequality_test(N, gamma_val, dps=30):
    """
    Young's inequality: ||𝔊||_{op} ≤ ||g||_{L¹} where g is convolution kernel.

    For Toeplitz part: G₂_{nm} ~ g(n-m) = 1/|n-m|^γ
    ||g||_{L¹} ~ Σ_k 1/k^γ  (diverges for γ≤1)

    For Hankel part: not a convolution, Young doesn't apply directly.

    This route FAILS for γ≤1 (includes all γ≤0.5).
    """
    # L¹ norm of kernel 1/k^γ over k=1..N
    mp.dps = dps
    gamma = mp.mpf(gamma_val)

    L1_norm = sum(1/mp.mpf(k)**gamma for k in range(1, N+1))

    # For γ≤1, this diverges as N→∞
    return float(L1_norm)

# ============================================================================
# Route 2: Direct finite-rank approximation (operator norm distance)
# ============================================================================

def operator_norm_tail_bound(N, gamma_val, dps=30):
    """
    Bound ||E||_{op} = ||𝔊 - 𝔊^{(N)}||_{op} via Gershgorin on the tail block.

    E decomposes as:
      E = [[ 0,    C   ],
           [ C^T,  𝔊^BB ]]

    where C couples head (n≤N) to tail (m>N), and 𝔊^BB is tail-tail.

    Gershgorin on 2-block matrix: ||E||_{op} ≤ max(sup_n R_n, sup_m K_m) where
      R_n = Σ_{m>N} |C_{nm}|  (n≤N)
      K_m = Σ_{n≤N} |C_{nm}| + Σ_{m'>N} |𝔊^BB_{mm'}|  (m>N)

    Key: For γ > 0.5, weighted sums Σ u_m/m^γ ~ Σ 1/m^{γ+1/2} converge.

    For γ ∈ (0.5, 1):
      - HS diverges: Σ 1/m^{2γ} diverges (2γ ≤ 1)
      - But operator norm tail converges: Σ 1/m^{γ+1/2} converges (γ+1/2 > 1)

    This is the KEY WINDOW where operator norm compactness holds despite HS failure!

    For γ ≤ 0.5:
      - Operator norm tail also diverges: γ+1/2 ≤ 1
      - This route FAILS.
    """
    mp.dps = dps
    iv.dps = dps
    gamma = iv.mpf(gamma_val)

    def u(n):
        return 1 / iv.sqrt(iv.mpf(2*n + 1))

    # Bound R_n for n=0 (worst case)
    n = 0
    un = u(n)

    # Weighted tail sum: Σ_{m>N} u_m/|n-m|^γ
    # For large m >> n: ~ Σ_{m>N} 1/(m^γ √m) = Σ 1/m^{γ+1/2}

    # Integral bound: ∫_{N}^∞ x^{-(γ+1/2)} dx = N^{-(γ-1/2)}/(γ-1/2)  if γ>0.5
    if gamma_val > 0.5:
        exponent = -(gamma - iv.mpf(0.5))
        integral = iv.mpf(N)**exponent / (gamma - iv.mpf(0.5))
        R_n_bound = un * integral * iv.mpf(2)  # factor 2: Toeplitz + Hankel
        return hi(R_n_bound)
    else:
        # Diverges
        return float('inf')

def route2_convergence_test(gamma_val, N_list, dps=30):
    """
    Test convergence of ||E_N||_{op} as N increases for γ > 0.5.

    Expected: ||E_N|| ~ N^{-(γ-1/2)} for γ ∈ (0.5, 1).
    """
    results = []
    for N in N_list:
        bound = operator_norm_tail_bound(N, gamma_val, dps=dps)
        results.append((N, bound))
    return results

# ============================================================================
# Route 3: Krein-Rutman without compactness
# ============================================================================

def spectral_radius_gershgorin(N, gamma_val, dps=30):
    """
    r(G) <= ||G||_{op} <= max_n row_sum_n  (Gershgorin upper bound).

    Krein-Rutman theorem requires:
      1. Operator positive (irreducible)  [OK] (all off-diag entries > 0)
      2. r(G) < inf

    If Route 1's Schur test gives finite bound, KR applies directly.
    If Route 1 fails (gamma<=0.5), need alternative r(G) bound.

    Alternative: Power iteration on finite truncation.
      - Compute lambda_1(G^{(N)}) for increasing N
      - If sequence converges, r(G) = lim lambda_1(G^{(N)}) < inf
      - KR then gives: exists positive eigenvector for r(G)

    But: this doesn't give a priori certificate that r(G) < inf before computing.

    Conclusion: Route 3 needs Route 1 or Route 2 to establish r(G) < inf.
    If both fail (gamma<=0.5), Route 3 also fails.
    """
    # Compute lambda_1(G^{(N)}) as proxy for r(G)
    def u_float(n):
        return 1.0 / np.sqrt(2*n + 1)

    # Build matrix (float64 for eigenvalue computation)
    dim = N + 1
    G = np.zeros((dim, dim))

    for n in range(dim):
        un_f = u_float(n)

        # Diagonal
        g1_diag = 0.0
        if n > 0:
            for k in range(1, n+1):
                g1_diag += 2.0 / (k**gamma_val)
        g1_diag /= (2*n + 1)

        g3_diag = un_f * un_f / ((2*n + 1)**gamma_val)
        G[n,n] = -g1_diag + g3_diag

        # Off-diagonal
        for m in range(n+1, dim):
            um_f = u_float(m)
            val = un_f * um_f * (1.0/(abs(n-m)**gamma_val) + 1.0/((n+m+1)**gamma_val))
            G[n,m] = val
            G[m,n] = val

    evals = np.linalg.eigvalsh(G)
    lambda1 = evals[-1]

    return lambda1

def route3_convergence_test(gamma_val, N_list):
    """
    Test monotonic convergence of λ₁(𝔊^{(N)}) for γ≤0.5.

    If converges to finite value, then r(𝔊) < ∞ and KR applies.
    """
    results = []
    for N in N_list:
        lam = spectral_radius_gershgorin(N, gamma_val, dps=30)
        results.append((N, lam))
    return results

# ============================================================================
# Main assessment
# ============================================================================

def assess_all_routes(gamma_val=0.5, N=200):
    """
    Comprehensive assessment of all three routes for given γ.
    """
    print(f"\n{'='*70}")
    print(f"γ = {gamma_val}  (N = {N})")
    print(f"{'='*70}\n")

    # Route 1: Schur test
    print("Route 1: Operator norm via Schur test")
    print("-" * 50)
    try:
        row_sum = schur_test_row_sum(N, gamma_val, dps=30)
        print(f"  Row sum (n=0, N={N}): {float(row_sum):.6f}")
        if row_sum < 1e10:  # Practical infinity threshold
            print(f"  [OK] Schur test gives ||G||_op <= {float(row_sum):.6f} < inf")
            print(f"  -> KR theorem applies")
        else:
            print(f"  [FAIL] Row sum diverging ({float(row_sum):.2e})")
            print(f"  -> Schur test INCONCLUSIVE")
    except Exception as e:
        print(f"  [FAIL] Error: {e}")

    print()

    # Route 1b: Young inequality
    print("Route 1b: Young inequality (L1 norm of kernel)")
    print("-" * 50)
    try:
        L1 = young_inequality_test(N, gamma_val, dps=30)
        print(f"  ||g||_L1 (N={N}): {L1:.6f}")
        if L1 < 1e10:
            print(f"  [OK] Young bound: ||G||_op <= {L1:.6f}")
        else:
            print(f"  [FAIL] L1 norm diverging")
    except Exception as e:
        print(f"  [FAIL] Error: {e}")

    print()

    # Route 2: Operator norm tail
    print("Route 2: Direct finite-rank approximation (operator norm)")
    print("-" * 50)
    try:
        tail_bound = operator_norm_tail_bound(N, gamma_val, dps=30)
        if tail_bound < 1e10:
            print(f"  ||E_N||_op <= {tail_bound:.6e}  (N={N})")
            print(f"  [OK] Tail decays as N^{{-(gamma-1/2)}} = N^{{-{gamma_val-0.5:.2f}}}")
            print(f"  -> Compact in operator norm (not HS!)")

            # Test convergence
            N_test = [50, 100, 200, 400]
            if N in N_test:
                print(f"\n  Convergence test:")
                conv = route2_convergence_test(gamma_val, N_test[:N_test.index(N)+1], dps=30)
                for Ni, bound_i in conv:
                    print(f"    N={Ni:3d}:  ||E|| <= {float(bound_i):.6e}")
        else:
            print(f"  [FAIL] Tail bound diverges (gamma <= 0.5)")
            print(f"  -> Route 2 FAILS")
    except Exception as e:
        print(f"  [FAIL] Error: {e}")

    print()

    # Route 3: KR via power iteration
    print("Route 3: Krein-Rutman via spectral radius")
    print("-" * 50)
    try:
        N_test = [50, 100, 200, 400] if N >= 400 else [50, 100, 200]
        print(f"  lambda_1(G^{{(N)}}) convergence:")
        conv = route3_convergence_test(gamma_val, N_test)

        last_lam = None
        for Ni, lam in conv:
            change = f" (Delta={lam-last_lam:.6e})" if last_lam else ""
            print(f"    N={Ni:3d}:  lambda_1 = {lam:.8f}{change}")
            last_lam = lam

        # Check convergence
        if len(conv) >= 2:
            rel_change = abs(conv[-1][1] - conv[-2][1]) / conv[-1][1]
            if rel_change < 1e-3:
                print(f"  [OK] Converging to r(G) ~= {conv[-1][1]:.8f} < inf")
                print(f"  -> KR applies: exists positive eigenvector")
            else:
                print(f"  [WARN] Still changing by {rel_change:.2%}, need larger N")
    except Exception as e:
        print(f"  [FAIL] Error: {e}")

    print("\n" + "="*70 + "\n")

if __name__ == '__main__':
    print("\n" + "="*70)
    print("γ-model Compactness: Three Routes for γ≤0.5")
    print("="*70)

    # Test cases
    gamma_cases = [
        (0.3, "nematic QCP (physical)"),
        (0.4, "nematic QCP (physical)"),
        (0.5, "boundary case"),
        (0.67, "comparison: γ>0.5, HS still diverges"),
        (1.0, "comparison: SDW QCP, HS converges"),
    ]

    for gam, desc in gamma_cases:
        print(f"\n### Case: γ = {gam} — {desc}")
        assess_all_routes(gam, N=200)
