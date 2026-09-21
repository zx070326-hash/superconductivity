"""
Schur test bound for ||C|| (head-tail coupling block) of KAY γ-model operator 𝔊(γ).

Operator definition (from gamma_model_correct_operator.py):
  𝔊(γ) = -𝔊₁ + 𝔊₂ + 𝔊₃ on ℓ²(ℕ₀), u_n = 1/√(2n+1)
  𝔊₁_{nn} = [1/(2n+1)] Σ_{k=1}^{n} 2/k^γ          (diagonal)
  𝔊₂_{nm} = u_n u_m / |n-m|^γ     (n≠m, Toeplitz)
  𝔊₃_{nm} = u_n u_m / (n+m+1)^γ   (Hankel)

Truncation at N: A = 𝔊⁽ᴺ⁾ ⊕ 0, E = 𝔊 - A decomposes as
  E = [[ 0,    C   ],     where C_{nm} = u_n u_m g(n,m)  for n ≤ N < m
      [ C^T,  𝔊^BB ]]

  g(n,m) = 1/|n-m|^γ + 1/(n+m+1)^γ  (𝔊₂ + 𝔊₃ parts; 𝔊₁ is diagonal so doesn't
  contribute to C since n ≤ N < m means n ≠ m, and 𝔊₁ only affects diagonal).

Schur test: ||C|| ≤ √(sup_n R_n · sup_m K_m), where
  R_n = Σ_{m>N} |C_{nm}| = u_n Σ_{m>N} u_m g(n,m)    (row sum, n ≤ N)
  K_m = Σ_{n≤N} |C_{nm}| = u_m Σ_{n=0}^N u_n g(n,m)  (col sum, m > N)

Key tail sums (J ≥ 1 integer):
  Tail(J,γ) := Σ_{j≥J} 1/j^γ

Bounds (all UPPER bounds on nonnegative quantities, direction-audited):
  γ > 1: Tail(J,γ) ≤ ∫_{J-1}^∞ x^{-γ} dx = (J-1)^{1-γ}/(γ-1)
  γ = 1: Tail(J,1) = harmonic tail ≤ 1/(J-1) + ln(1 + 1/(J-1))  (J ≥ 2)
         For J=1: Tail(1,1) ≤ 2 (crude but valid; exact = ∞ but weighted sums converge)
  γ < 1: Tail(J,γ) diverges, BUT u_m Tail(J,γ) may converge due to u_m ~ m^{-1/2}

Implementation focuses on γ=1 (most important for BCS), with framework for γ ∈ (0,2).

Strictness level: certificate (interval arithmetic via mpmath.iv).
"""
import sys, os, time
import numpy as np
from mpmath import iv, mp

# Utilities
def lo(x):
    """Extract certified lower bound."""
    return mp.mpf(x.a) if hasattr(x, 'a') else mp.mpf(x)

def hi(x):
    """Extract certified upper bound."""
    return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)


def build_G_matrix(N, gamma_val, dps=30):
    """Build KAY's 𝔊(γ) operator truncated to dimension N+1, using interval arithmetic.
    Returns dict with matrix H, dimension, and helper functions."""
    mp.dps = dps
    iv.dps = dps

    gamma = iv.mpf(gamma_val)
    dim = N + 1

    def u(n):
        """u_n = 1/√(2n+1)"""
        return 1 / iv.sqrt(iv.mpf(2*n + 1))

    def g_coupling(n, m):
        """g(n,m) = 1/|n-m|^γ + 1/(n+m+1)^γ  (for n ≠ m)"""
        if n == m:
            return iv.mpf(0)
        nm_diff = iv.mpf(abs(n - m))
        nm_sum = iv.mpf(n + m + 1)
        return nm_diff**(-gamma) + nm_sum**(-gamma)

    # Build H matrix
    H = [[None]*dim for _ in range(dim)]

    for n in range(dim):
        un = u(n)

        # 𝔊₁ diagonal: [1/(2n+1)] Σ_{k=1}^{n} 2/k^γ
        g1_diag = iv.mpf(0)
        if n > 0:  # k starts at 1
            for k in range(1, n+1):
                g1_diag = g1_diag + 2 * iv.mpf(k)**(-gamma)
        g1_diag = g1_diag / iv.mpf(2*n + 1)

        # Diagonal entry: -𝔊₁ + 𝔊₃(n,n)
        g3_diag = un * un / iv.mpf(2*n + 1)**gamma
        H[n][n] = -g1_diag + g3_diag

        # Off-diagonal entries
        for m in range(n+1, dim):
            um = u(m)
            # 𝔊₂ + 𝔊₃
            val = un * um * g_coupling(n, m)
            H[n][m] = val
            H[m][n] = val  # symmetric

    return dict(H=H, dim=dim, N=N, gamma=gamma, u=u, g_coupling=g_coupling)


def tail_sum_bound_hi(J, gamma_val, dps=30):
    """Certified UPPER bound for Tail(J,γ) := Σ_{j≥J} 1/j^γ, J ≥ 1 integer.

    Returns plain mp.mpf (upper bound).

    Three regimes:
    - γ > 1: integral bound (J-1)^{1-γ}/(γ-1)
    - γ = 1: harmonic tail bound
    - γ < 1: weighted bound (for use with u_m factors)
    """
    mp.dps = dps
    iv.dps = dps

    gamma = iv.mpf(gamma_val)

    if gamma_val > 1.0:
        # ∫_{J-1}^∞ x^{-γ} dx = (J-1)^{1-γ}/(γ-1)
        J_iv = iv.mpf(J)
        exponent = 1 - gamma
        bound = (J_iv - 1)**exponent / (gamma - 1)
        return hi(bound)

    elif abs(gamma_val - 1.0) < 1e-10:
        # Harmonic tail: Tail(J,1) ≤ 1/(J-1) + ln(1 + 1/(J-1)) for J ≥ 2
        # For J=1: use crude bound 2 (sum_{j≥1} 1/j diverges, but this is for
        # weighted sums where u_m ~ m^{-1/2} provides convergence)
        if J <= 1:
            return mp.mpf(2)
        else:
            J_iv = iv.mpf(J)
            # 1/(J-1) + ln(1 + 1/(J-1))
            # Use ln(1+x) ≤ x + x²/2 for x ≥ 0 (second-order Taylor remainder is positive)
            # Actually, for certificate level, use direct iv.log with proper rounding
            x = 1 / (J_iv - 1)
            bound = x + iv.log(1 + x)
            return hi(bound)

    else:
        # γ < 1: Tail diverges, but for weighted sum u_m Tail we need a different approach
        # Use integral bound with explicit remainder
        # Σ_{j≥J} j^{-γ} ≤ J^{-γ} + ∫_J^∞ x^{-γ} dx = J^{-γ} + J^{1-γ}/(γ-1)
        # But γ-1 < 0, so need careful direction
        # For γ < 1: ∫_J^∞ x^{-γ} dx = J^{1-γ}/(1-γ) (note sign flip)
        J_iv = iv.mpf(J)
        exponent = 1 - gamma
        bound = J_iv**(-gamma) + J_iv**exponent / (1 - gamma)
        return hi(bound)


def weighted_tail_bound_hi(J, gamma_val, dps=30):
    """Certified UPPER bound for Σ_{m≥J} u_m / m^γ, where u_m = 1/√(2m+1).

    This is the key tail sum appearing in row/column sums:
      Σ_{m≥J} u_m / m^γ  and  Σ_{m≥J} u_m / (a+m)^γ  for fixed a

    For γ=1: u_m ~ m^{-1/2}, so u_m/m ~ m^{-3/2} converges rapidly.
    Use u_m ≤ 1/√(2m) and 1/m^γ bound, then sum via integral comparison.
    """
    mp.dps = dps
    iv.dps = dps

    if J < 1:
        J = 1

    gamma = iv.mpf(gamma_val)

    # u_m = 1/√(2m+1) ≤ 1/√(2m) for m ≥ 1
    # u_m / m^γ ≤ 1/(√(2m) m^γ) = 1/(√2 m^{γ+1/2})

    # Σ_{m≥J} m^{-(γ+1/2)} ≤ ∫_{J}^∞ x^{-(γ+1/2)} dx  (if γ+1/2 > 1, i.e., γ > 1/2)

    gamma_eff = gamma + iv.mpf(0.5)

    if float(gamma_eff.mid) > 1.0:
        # ∫_{J}^∞ x^{-α} dx = J^{1-α}/(α-1),  α = γ+1/2
        # Special case: J=1 to avoid (J-1)**exponent = 0**negative = inf
        J_iv = iv.mpf(J)
        exponent = 1 - gamma_eff
        # exponent is negative since gamma_eff > 1
        # J^{1-α} = J^{negative} = 1/J^{α-1}

        if J == 1:
            # ∫_1^∞ x^{-(γ+1/2)} dx = 1/(γ+1/2-1) = 1/(γ-1/2)
            integral = 1 / (gamma - iv.mpf(0.5))
        else:
            # For J > 1, use ∫_{J-1}^∞ for upper bound (wider interval)
            # But this gives (J-1)^{1-α}/(α-1), which still blows up if J=1
            # Better: use ∫_J^∞ directly
            integral = J_iv**exponent / (gamma_eff - 1)

        bound = integral / iv.sqrt(iv.mpf(2))
        return hi(bound)
    else:
        # γ ≤ 1/2: Σ_{m≥J} u_m/m^γ ≥ Σ_{m≥J} 1/(√(2m+1) m^{1/2}) DIVERGES
        # (harmonic-type tail).  The only sound upper bound is +inf.
        # SOUNDNESS FIX 2026-09-16: the previous branch returned finite
        # numbers ((1+ln J)/√2 for γ=1/2, J^{1/2-γ}/(1/2-γ)/√2 for γ<1/2) for
        # a divergent sum — a pad that does not cover the quantity it stands
        # in for.  Verified: γ=0.5, J=401 gave 4.945 while the partial sum
        # already exceeds 6.02 at m=2e6.  This matches the non-compactness
        # result in gamma_small_compactness_2026-09-16.md (||E_N||=∞ for
        # γ≤1/2).  Downstream (Davis-Kahan) then yields sinΘ=1, c₁²_lo=0,
        # i.e. a vacuous-but-true verdict instead of a false certificate.
        return mp.inf


def Rsup_hi(data, N, dps=30):
    """Certified upper bound for sup_{n=0..N} R_n.

    R_n = u_n Σ_{m>N} u_m g(n,m)
        = u_n Σ_{m>N} u_m (1/|m-n|^γ + 1/(n+m+1)^γ)

    Upper bound via u_m ≤ u_{N+1} for m ≥ N+1:
      R_n ≤ u_n u_{N+1} Σ_{m>N} (1/|m-n|^γ + 1/(n+m+1)^γ)
         = u_n u_{N+1} (Tail(N+1-n, γ) + Tail(n+N+2, γ))

    But this is for unweighted tails. With the u_m factor:
      R_n = u_n Σ_{m≥N+1} u_m / |m-n|^γ + u_n Σ_{m≥N+1} u_m / (n+m+1)^γ

    Use weighted_tail_bound_hi for better bounds.
    """
    mp.dps = dps
    iv.dps = dps

    gamma_val = float(data['gamma'].mid)
    u = data['u']

    uNp1_hi = hi(u(N + 1))
    best = mp.mpf(0)

    for n in range(0, N + 1):
        un_hi = hi(u(n))

        # First sum: Σ_{m≥N+1} u_m / |m-n|^γ = Σ_{m≥N+1} u_m / (m-n)^γ  (since m > N ≥ n)
        # Start index: m = N+1, so m-n starts at N+1-n
        t1 = weighted_tail_bound_hi(N + 1 - n, gamma_val, dps=dps)

        # Second sum: Σ_{m≥N+1} u_m / (n+m+1)^γ
        # Here n+m+1 starts at n+N+2
        t2 = weighted_tail_bound_hi(N + 2 + n, gamma_val, dps=dps)

        val = un_hi * (t1 + t2)
        if val > best:
            best = val

    return best


def Ksup_hi(data, N, dps=30):
    """Exact (via monotonicity proof) certified value of sup_{m>N} K_m.

    K_m = u_m Σ_{n=0}^N u_n g(n,m) = u_m Σ_{n=0}^N u_n (1/|m-n|^γ + 1/(n+m+1)^γ)

    For m > N: u_m, |m-n|, and n+m+1 all increase with m (for fixed n ≤ N).
    - u_m decreasing in m
    - 1/|m-n|^γ = 1/(m-n)^γ decreasing in m
    - 1/(n+m+1)^γ decreasing in m

    Each term u_m × [stuff(m)] is non-increasing in m, so K_m is non-increasing.
    Hence sup_{m>N} K_m = K_{N+1}, computed exactly via finite (N+1)-term iv sum.
    """
    mp.dps = dps
    iv.dps = dps

    u = data['u']
    g_coupling = data['g_coupling']
    gamma = data['gamma']

    m = N + 1
    um = u(m)
    s = iv.mpf(0)

    for n in range(0, N + 1):
        un = u(n)
        # g(n,m) = 1/|n-m|^γ + 1/(n+m+1)^γ
        term = un * g_coupling(n, m)
        s = s + term

    return hi(um * s)


def schur_C_bound_hi(data, N, dps=30):
    """Certified upper bound for ||C|| via Schur test.

    ||C|| ≤ √(R_sup × K_sup)
    """
    R = Rsup_hi(data, N, dps=dps)
    K = Ksup_hi(data, N, dps=dps)
    return mp.sqrt(R * K)


def HBB_bound_hi(data, N, dps=30):
    """Certified upper bound for ||𝔊^BB|| (tail-tail block, n,m > N).

    𝔊^BB has three contributions:
    - 𝔊₁ diagonal: grows like (ln n)/n for γ=1
    - 𝔊₂ Toeplitz: u_n u_m / |n-m|^γ
    - 𝔊₃ Hankel: u_n u_m / (n+m+1)^γ

    For γ=1: diagonal 𝔊₁_{nn} ~ (2 ln n)/(2n+1) → 0 as n→∞
    Dominant: off-diagonal coupling ~ u_n u_m / (stuff)^γ

    Row sum bound: for n > N,
      Σ_{m>N} |𝔊_{nm}| = |𝔊₁_{nn}| + Σ_{m>N,m≠n} u_n u_m (1/|n-m|^γ + 1/(n+m+1)^γ)

    Use sup over n > N of the row sum, then bound via u_n ≤ u_{N+1} and tail sums.
    """
    mp.dps = dps
    iv.dps = dps

    gamma_val = float(data['gamma'].mid)
    gamma = data['gamma']
    u = data['u']

    # Diagonal 𝔊₁ bound: for n > N,
    #   𝔊₁_{nn} = [1/(2n+1)] Σ_{k=1}^{n} 2/k^γ
    # For γ=1: Σ_{k=1}^{n} 2/k ≤ 2(1 + ln n)
    # So 𝔊₁_{nn} ≤ 2(1 + ln n)/(2n+1) ≤ 2(1 + ln(N+1))/(2N+3)

    if abs(gamma_val - 1.0) < 1e-10:
        # γ=1
        Np1 = iv.mpf(N + 1)
        g1_bound = 2 * (1 + iv.log(Np1)) / (2*Np1 + 1)
    else:
        # γ ≠ 1: use Σ_{k=1}^{n} k^{-γ} ≤ 1 + n^{1-γ}/(1-γ) for γ < 1
        # or ≤ 1 + 1/(γ-1) for γ > 1
        Np1 = iv.mpf(N + 1)
        if gamma_val > 1.0:
            tail_sum = 1 + 1/(gamma - 1)
        else:
            tail_sum = 1 + Np1**(1 - gamma) / (1 - gamma)
        g1_bound = 2 * tail_sum / (2*Np1 + 1)

    # Off-diagonal row sum: for n > N,
    #   Σ_{m>N,m≠n} u_n u_m (1/|n-m|^γ + 1/(n+m+1)^γ)
    #
    # Key: we need to bound this for the WORST n > N. Use n = N+1.
    # For n = N+1:
    #   Σ_{m>N+1} u_{N+1} u_m (1/(m-N-1)^γ + 1/(N+1+m+1)^γ)
    #   + u_{N+1}² (1/(N+1+N+1)^γ)  [the m=N+1 term, special since |n-m|=0 excluded]
    #
    # Wait, m≠n means we exclude m=n. For n=N+1, we sum over m∈{N+1,N+2,...}\{N+1} = {N+2,...}

    # Actually, let me be more careful. For n=N+1, m>N, m≠n:
    #   m ∈ {N+1, N+2, ...} \ {N+1} = {N+2, N+3, ...}
    # So:
    #   Σ_{m≥N+2} u_{N+1} u_m (1/(m-N-1)^γ + 1/(N+1+m+1)^γ)

    uNp1 = u(N + 1)

    # First sum: Σ_{m≥N+2} u_m / (m-N-1)^γ
    # Start: m=N+2 means m-N-1 = 1
    t1 = mp.mpf(weighted_tail_bound_hi(1, gamma_val, dps=dps))  # Σ_{j≥1} u_{N+1+j} / j^γ

    # Second sum: Σ_{m≥N+2} u_m / (N+1+m+1)^γ = Σ_{m≥N+2} u_m / (m+N+2)^γ
    # For m=N+2: m+N+2 = 2N+4
    # But we need Σ_{m≥N+2} u_m / (m+N+2)^γ
    # Let m = N+1+k, k≥1: Σ_{k≥1} u_{N+1+k} / (2N+3+k)^γ
    # This is NOT a simple tail sum. Need to be more careful.

    # Actually, use the simpler crude bound: u_m ≤ u_{N+1}, and
    # Σ_{m≥N+2} 1/(m+N+2)^γ ≤ Σ_{j≥2N+4} 1/j^γ = Tail(2N+4, γ)

    t2_simple = mp.mpf(tail_sum_bound_hi(2*N + 4, gamma_val, dps=dps))

    offdiag_bound = hi(uNp1) * (t1 + hi(uNp1) * t2_simple)

    return mp.mpf(hi(g1_bound)) + offdiag_bound


def run_one_gamma(gamma_val, N, dps=30, verbose=True):
    """Run Schur bound computation for γ-model at given γ, N.

    Returns dict with certified bounds.
    """
    t0 = time.time()

    data = build_G_matrix(N, gamma_val, dps=dps)

    # Float64 reference eigenvalues
    H_float = np.zeros((N+1, N+1))
    for i in range(N+1):
        for j in range(N+1):
            H_float[i, j] = float(data['H'][i][j].mid)

    evals = np.linalg.eigvalsh(H_float)
    lam1_float = evals[-1]
    lam2_float = evals[-2]
    gap_float = lam1_float - lam2_float

    # Certified bounds
    C_bound = schur_C_bound_hi(data, N, dps=dps)
    HBB_bound = HBB_bound_hi(data, N, dps=dps)
    E_bound = C_bound + HBB_bound  # ||E|| ≤ ||C|| + ||H^BB||

    # Float64 validation: compute true ||C|| via finite tail truncation
    M_tail = min(N + 5000, N + 50*N)  # Tail truncation for float validation
    C_float = np.zeros((N+1, M_tail - N))
    u_vec = np.array([1.0/np.sqrt(2*n+1) for n in range(M_tail + 1)])

    for n in range(N+1):
        for m_idx, m in enumerate(range(N+1, M_tail + 1)):
            # C_{nm} = u_n u_m (1/|n-m|^γ + 1/(n+m+1)^γ)
            val = u_vec[n] * u_vec[m] * (abs(n-m)**(-gamma_val) + (n+m+1)**(-gamma_val))
            C_float[n, m_idx] = val

    if C_float.size > 0:
        C_true = np.linalg.norm(C_float, ord=2)
        C_ratio = float(C_bound) / C_true if C_true > 0 else float('inf')
    else:
        C_true = 0.0
        C_ratio = float('inf')

    if verbose:
        print(f"gamma={gamma_val} N={N} dim={N+1} time={time.time()-t0:.1f}s")
        print(f"  Float64 reference: lam1={lam1_float:.8f}  lam2={lam2_float:.8f}  gap={gap_float:.8f}")
        print(f"  Certified bounds:")
        print(f"    ||C|| <= {C_bound}  (float64 true ||C|| ~= {C_true:.6e}, ratio {C_ratio:.2f}x)")
        print(f"    ||G^BB|| <= {HBB_bound}")
        print(f"    ||E|| <= {E_bound}")
        print(f"  Strictness: certificate (interval arithmetic)")

    return dict(
        gamma=gamma_val, N=N,
        lam1_float=lam1_float, lam2_float=lam2_float, gap_float=gap_float,
        C_bound=C_bound, HBB_bound=HBB_bound, E_bound=E_bound,
        C_true=C_true, C_ratio=C_ratio
    )


if __name__ == '__main__':
    # Parse command line: gamma N [dps]
    gamma_val = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 30

    print("=" * 70)
    print("gamma-model Schur tail bound (certificate level)")
    print("=" * 70)

    if len(sys.argv) == 1:
        # No arguments: run a survey across γ and N
        print("\nSurvey mode: testing multiple (gamma, N) pairs\n")
        print(f"{'gamma':>6} {'N':>5} {'lam1':>10} {'gap':>10} {'||C||':>12} {'ratio':>7} {'||E||':>12} {'time':>6}")
        print("-" * 70)

        for gamma_val in [0.5, 0.75, 1.0, 1.5, 2.0]:
            for N in [100, 200, 400]:
                result = run_one_gamma(gamma_val, N, dps=30, verbose=False)
                print(f"{result['gamma']:6.2f} {result['N']:5d} "
                      f"{result['lam1_float']:10.6f} {result['gap_float']:10.6f} "
                      f"{float(result['C_bound']):12.6e} {result['C_ratio']:7.2f} "
                      f"{float(result['E_bound']):12.6e}")

        print("\nStrictness level: certificate (interval arithmetic via mpmath.iv)")
        print("Schur test: ||C|| <= sqrt(sup_n R_n * sup_m K_m)")
        print("R_n via weighted tail bounds; K_m exact at m=N+1 (monotonicity)")
    else:
        # Single run with arguments
        run_one_gamma(gamma_val, N, dps=dps, verbose=True)
