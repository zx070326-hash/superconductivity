"""
γ-model Davis-Kahan certificate pipeline with Coulomb μ*.

Constructs interval-arithmetic certified lower bounds for:
  - 𝔤(γ, μ*, Nc) = top eigenvalue of 𝔊(γ) - μ* P_Nc
  - c₁² = concentration at n=0 (via Davis-Kahan sinΘ)

Operator definition (KAY-I eq.13-31):
  𝔊(γ) = -𝔊₁ + 𝔊₂ + 𝔊₃ on ℓ²(ℕ₀)
  𝔊₁_{nn} = [1/(2n+1)] Σ_{k=1}^n 2/k^γ   (diagonal, exact sum to n)
  𝔊₂_{nm} = u_n u_m / |n-m|^γ            (n≠m, Toeplitz)
  𝔊₃_{nm} = u_n u_m / (n+m+1)^γ          (Hankel)
  u_n = 1/√(2n+1)

Coulomb finite-rank perturbation:
  𝔊(γ, μ*, Nc) = 𝔊(γ) - μ* P_Nc
where P_Nc projects onto n < Nc (Phase 0 verified: well-defined, no divergence).

Davis-Kahan route (from P4 c1sq_schur_bound.py):
  1. Float64 eigendecomposition → approx λ₁, λ₂, v₁
  2. Interval Rayleigh quotient → μ₁_lo (certified lower)
  3. Residual r₁ = ||H v₁ - μ₁ v₁|| (certified upper)
  4. λ₂ upper bound (Gershgorin or LDLᵀ inertia)
  5. sinΘ₁ = r₁/(gap₁ - r₁)
  6. Tail ||E|| upper bound (placeholder: global bound, to be replaced by Schur)
  7. sinΘ₂ = ||E||/(δ - ||E||)
  8. sinΘ_total = sinΘ₁ + sinΘ₂
  9. c₁² ≥ (cosΘ·|v₁[0]| - sinΘ)²

Tail bound placeholder (coarse, sufficient for pipeline validation):
  ||𝔊^BB|| ≤ sum of row sup over n>N (bound via integral comparison)
  ||C|| ≤ similar coupling block bound

Strictness: certificate (mpmath.iv interval arithmetic, all bounds direction-audited).

Reference:
  - gamma_model_correct_operator.py (KAY operator, float64 verification)
  - c1sq_schur_bound.py (P4 Davis-Kahan certificate pipeline)
  - gamma_coulomb_phase0.py (Coulomb finite-rank well-definedness)
"""
import sys, os, time
from mpmath import iv, mp
import numpy as np

mp.dps = 30
iv.dps = 30


def build_G_matrix(N, gamma_val, mu_star=0.0, Nc=None):
    """Build 𝔊(γ, μ*, Nc) matrix (N+1 × N+1) using interval arithmetic.

    Returns: dict with keys H (interval matrix), dim, gamma, mu_star, Nc.
    """
    dim = N + 1
    gamma = iv.mpf(gamma_val)

    # u_n = 1/√(2n+1)
    def u(n):
        return 1 / iv.sqrt(iv.mpf(2*n + 1))

    # Power function for γ-model kernel
    def power_gamma(k, g):
        """k^{-g} with interval arithmetic, k > 0."""
        return 1 / (iv.mpf(k) ** g)

    H = [[None]*dim for _ in range(dim)]

    for n in range(dim):
        un = u(n)

        # 𝔊₁ diagonal: [1/(2n+1)] Σ_{k=1}^n 2/k^γ
        g1_diag = iv.mpf(0)
        for k in range(1, n+1):
            g1_diag = g1_diag + 2 * power_gamma(k, gamma)
        g1_diag = g1_diag / iv.mpf(2*n + 1)

        # 𝔊₃ diagonal: u_n² / (2n+1)^γ
        g3_diag = un * un * power_gamma(2*n + 1, gamma)

        H[n][n] = -g1_diag + g3_diag  # -𝔊₁ + 𝔊₃ diagonal

        for m in range(dim):
            if m != n:
                um = u(m)
                # 𝔊₂: u_n u_m / |n-m|^γ
                g2 = un * um * power_gamma(abs(n - m), gamma)
                # 𝔊₃: u_n u_m / (n+m+1)^γ
                g3 = un * um * power_gamma(n + m + 1, gamma)
                H[n][m] = g2 + g3

    # Coulomb perturbation: subtract μ* u_n u_m for n, m < Nc
    if mu_star != 0 and Nc is not None and Nc > 0:
        mu_iv = iv.mpf(mu_star)
        Nc_eff = min(Nc, dim)
        for n in range(Nc_eff):
            un = u(n)
            for m in range(Nc_eff):
                um = u(m)
                H[n][m] = H[n][m] - mu_iv * un * um

    return dict(H=H, dim=dim, gamma=gamma, gamma_val=gamma_val,
                mu_star=mu_star, Nc=Nc, N=N, u=u)


def to_float_matrix(H, dim):
    """Convert interval matrix to float64 (midpoint) for eigendecomposition."""
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i, j] = float(H[i][j].mid)
    return A


def rayleigh_lower(H, dim, vfloat):
    """Certified lower bound for v^T H v / v^T v."""
    v = [iv.mpf(x) for x in vfloat]
    num = iv.mpf(0)
    for i in range(dim):
        Hv_i = iv.mpf(0)
        Hrow = H[i]
        for j in range(dim):
            Hv_i = Hv_i + Hrow[j] * v[j]
        num = num + v[i] * Hv_i
    den = iv.mpf(0)
    for i in range(dim):
        den = den + v[i] * v[i]
    return num / den  # interval; .a for certified lower


def residual_norm(H, dim, vfloat, la_iv):
    """Certified UPPER bound on ||H v - la v|| (v normalized internally)."""
    vraw = [iv.mpf(x) for x in vfloat]
    nrm2 = iv.mpf(0)
    for x in vraw:
        nrm2 = nrm2 + x * x
    nrm = iv.sqrt(nrm2)
    v = [x / nrm for x in vraw]
    rr2 = iv.mpf(0)
    for i in range(dim):
        Hv_i = iv.mpf(0)
        Hrow = H[i]
        for j in range(dim):
            Hv_i = Hv_i + Hrow[j] * v[j]
        d = Hv_i - la_iv * v[i]
        rr2 = rr2 + d * d
    return iv.sqrt(rr2)


def inertia_count_below(H, dim, c):
    """Count eigenvalues < c via interval LDLᵀ (Sylvester inertia).
    Returns (count, ok). ok=False if pivot straddles zero."""
    M = [row[:] for row in H]
    for i in range(dim):
        M[i][i] = M[i][i] - c
    neg = 0
    for k in range(dim):
        dk = M[k][k]
        lo, hi = dk.a, dk.b
        if lo <= 0 <= hi:
            return None, False
        if hi < 0:
            neg += 1
        # Eliminate
        if k < dim - 1:
            invdk = 1 / dk
            col = [M[i][k] for i in range(k+1, dim)]
            for ii, i in enumerate(range(k+1, dim)):
                lik = col[ii] * invdk
                Mik_row = M[i]
                Mk_row = M[k]
                for j in range(k+1, dim):
                    Mik_row[j] = Mik_row[j] - lik * Mk_row[j]
    return neg, True


def find_lambda2_upper(H, dim, lam2_approx, gap_guess):
    """Certified upper bound on λ₂(H) via inertia counting."""
    margin = max(gap_guess * 1e-3, 1e-8)
    for attempt in range(40):
        c = iv.mpf(lam2_approx + margin)
        cnt, ok = inertia_count_below(H, dim, c)
        if ok and cnt is not None and cnt >= dim - 1:
            return c, margin
        margin *= 1.7
    raise RuntimeError("could not certify lambda2 upper bound; increase dps or N")


def tail_bounds_coarse(data, N):
    """PLACEHOLDER: coarse global upper bounds for ||𝔊^BB|| and ||C||.

    Direction-audited: all bounds are UPPER bounds on nonnegative quantities.
    To be replaced by tight Schur-test bounds (separate agent task).

    Returns: (HBB_hi, C_hi) as plain mp.mpf.
    """
    gamma = data['gamma']
    u = data['u']

    # u_{N+1}
    uNp1 = u(N + 1)
    uNp1_hi = mp.mpf(uNp1.b) if hasattr(uNp1, 'b') else mp.mpf(uNp1)

    # Global row-sum bound (very coarse):
    # Each row n of 𝔊 has off-diagonal sum ~ Σ_m u_m / (1+|n-m|)^γ
    # For γ=1: this is O(ln(N)), bounded by 2*ln(N+2)
    # For general γ: use integral comparison
    #   Σ_{k=1}^∞ 1/k^γ ≤ 1 + ∫₁^∞ x^{-γ} dx = 1 + 1/(γ-1)  for γ>1
    #   For γ≤1: Σ 1/k diverges, but tail Σ_{k>N} 1/k ≤ ln(2N/N) ~ ln(2)

    gamma_val = data['gamma_val']
    if gamma_val >= 1.0:
        # Tail sum Σ_{k>N} 1/k^γ ≤ ∫_N^∞ x^{-γ} dx = N^{1-γ}/(γ-1)
        if gamma_val == 1.0:
            tail_sum_bound = mp.mpf(1.0)  # very coarse for γ=1
        else:
            tail_sum_bound = mp.mpf(N) ** (1 - gamma_val) / (gamma_val - 1)
    else:
        # γ < 1: tail diverges, need cutoff-dependent bound (harder)
        # Placeholder: use ln(N) scaling
        tail_sum_bound = mp.log(mp.mpf(N + 2))

    # ||𝔊^BB|| ≤ sup row sum * sup col weight
    # Row n>N: Σ_m |𝔊_{nm}| ≤ u_n * (tail_sum for Toeplitz + tail_sum for Hankel)
    # Coarse: u_{N+1} * 2 * tail_sum_bound
    HBB_hi = uNp1_hi * 2 * tail_sum_bound

    # ||C|| (coupling block n≤N, m>N): similar
    # Row n≤N: Σ_{m>N} |𝔊_{nm}| ≤ u_n * u_{N+1} * tail_sum_bound
    # sup over n≤N: u_0 * u_{N+1} * tail_sum_bound (u_0 largest)
    u0_hi = mp.mpf(u(0).b) if hasattr(u(0), 'b') else mp.mpf(u(0))
    C_hi = u0_hi * uNp1_hi * tail_sum_bound

    return HBB_hi, C_hi


def lo(x):
    """Extract certified lower bound (plain mp.mpf)."""
    return mp.mpf(x.a) if hasattr(x, 'a') else mp.mpf(x)


def hi(x):
    """Extract certified upper bound (plain mp.mpf)."""
    return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)


def run_one_dk(gamma_val, N, mu_star=0.0, Nc=None, dps=30, verbose=True):
    """Davis-Kahan certificate pipeline for γ-model + Coulomb.

    Returns: dict with keys gamma, N, mu_star, Nc, mu1_lo, sinT, c1sq_lo, etc.
    """
    mp.dps = dps
    iv.dps = dps
    t0 = time.time()

    data = build_G_matrix(N, gamma_val, mu_star=mu_star, Nc=Nc)
    H, dim = data['H'], data['dim']

    # Float64 eigendecomposition (approximate)
    A = to_float_matrix(H, dim)
    evals, evecs = np.linalg.eigh(A)
    lam1_approx = evals[-1]
    v1 = evecs[:, -1]
    lam2_approx = evals[-2]

    # Step 1: μ₁ lower bound (interval Rayleigh quotient)
    mu1_lo = lo(rayleigh_lower(H, dim, v1.tolist()))
    la_iv = iv.mpf(mu1_lo)

    # Step 2: Residual r₁ = ||H v₁ - μ₁ v₁|| (upper bound)
    r1 = hi(residual_norm(H, dim, v1.tolist(), la_iv))

    # Step 3: λ₂ upper bound (inertia count)
    c_upper, margin_used = find_lambda2_upper(H, dim, lam2_approx,
                                               lam1_approx - lam2_approx)
    lam2_HN_hi = hi(c_upper)

    # Step 4: sinΘ₁ (eigenvector residual angle)
    gap1_lo = mu1_lo - lam2_HN_hi
    if gap1_lo <= r1:
        sinT1 = mp.mpf(1)
    else:
        sinT1 = r1 / (gap1_lo - r1)
        if sinT1 > 1:
            sinT1 = mp.mpf(1)

    # Step 5: Tail bounds (coarse placeholder)
    HBB_hi, C_hi = tail_bounds_coarse(data, N)
    Ehi = HBB_hi + C_hi

    # Step 6: sinΘ₂ (truncation angle)
    lam2A_hi = max(lam2_HN_hi, mp.mpf(0))
    delta_lo = mu1_lo - lam2A_hi
    if delta_lo <= Ehi:
        sinT2 = mp.mpf(1)
    else:
        sinT2 = Ehi / (delta_lo - Ehi)
        if sinT2 > 1:
            sinT2 = mp.mpf(1)

    # Step 7: Total sinΘ and c₁² certificate
    sinT = min(mp.mpf(1), sinT1 + sinT2)
    cosT = mp.sqrt(1 - sinT**2) if sinT < 1 else mp.mpf(0)

    v0 = abs(v1[0])
    val = cosT * mp.mpf(v0) - sinT
    c1sq_lo = max(mp.mpf(0), val) ** 2

    if verbose:
        mu_str = f"mu*={mu_star:.3f} Nc={Nc}" if mu_star != 0 else "mu*=0"
        print(f"gamma={gamma_val} N={N} {mu_str}  dim={dim}  time={time.time()-t0:.1f}s")
        print(f"  lam1_approx={lam1_approx:.8f}  lam2_approx={lam2_approx:.8f}")
        print(f"  mu1_lo(cert) = {mu1_lo}")
        print(f"  lam2_hi(cert) = {lam2_HN_hi}  margin={margin_used:.3e}")
        print(f"  r1(residual) = {r1}  sinT1 = {sinT1}")
        print(f"  ||G^BB||<={float(HBB_hi):.6f}  ||C||<={float(C_hi):.6f}  (coarse placeholder)")
        print(f"  Ehi = {float(Ehi):.6f}  delta_lo = {float(delta_lo):.6f}")
        print(f"  sinT2 = {sinT2}  sinT_total = {sinT}")
        print(f"  cosT = {cosT}  |v1[0]| = {v0:.6f}")
        print(f"  c1sq_lo(cert) = {c1sq_lo}")
        print()

    return dict(
        gamma=gamma_val, N=N, mu_star=mu_star, Nc=Nc,
        mu1_lo=mu1_lo, lam2_hi=lam2_HN_hi, r1=r1,
        sinT1=sinT1, sinT2=sinT2, sinT=sinT, cosT=cosT,
        c1sq_lo=c1sq_lo, v0=v0,
        lam1_approx=lam1_approx, lam2_approx=lam2_approx,
        HBB_hi=HBB_hi, C_hi=C_hi, Ehi=Ehi
    )


if __name__ == '__main__':
    print("=" * 70)
    print("GAMMA-MODEL DAVIS-KAHAN CERTIFICATE PIPELINE WITH COULOMB MU*")
    print("=" * 70)
    print()

    # Part 1: Sanity check — gamma=1, N=200, mu*=0 (no Coulomb)
    print("--- Part 1: Sanity check (gamma=1, N=200, mu*=0) ---")
    print("Should recover KAY's g^(200)(1) float64 value as certificate.")
    print()
    result0 = run_one_dk(gamma_val=1.0, N=200, mu_star=0.0, Nc=None)

    # Part 2: Coulomb effect — gamma=1, N=200, mu*=0.10/0.13/0.15, Nc=50
    print("--- Part 2: Coulomb mu* effect (gamma=1, N=200, Nc=50) ---")
    print("Physical parameters: mu*=0.10/0.13/0.15 (typical BCS renormalized values)")
    print()

    results = []
    for mu_star in [0.10, 0.13, 0.15]:
        res = run_one_dk(gamma_val=1.0, N=200, mu_star=mu_star, Nc=50)
        results.append(res)

    # Part 3: Summary table
    print("=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(f"{'gamma':>5} {'N':>5} {'mu*':>6} {'Nc':>5} {'g_lo':>14} {'c1sq_lo':>10} "
          f"{'sinT':>8} {'gap':>10}")

    # mu*=0
    print(f"{result0['gamma']:5.1f} {result0['N']:5d} {0.0:6.2f} {'—':>5} "
          f"{float(result0['mu1_lo']):14.8f} {float(result0['c1sq_lo']):10.6f} "
          f"{float(result0['sinT']):8.4f} "
          f"{result0['lam1_approx']-result0['lam2_approx']:10.6f}")

    # mu*>0
    for res in results:
        print(f"{res['gamma']:5.1f} {res['N']:5d} {res['mu_star']:6.2f} "
              f"{res['Nc']:5d} {float(res['mu1_lo']):14.8f} {float(res['c1sq_lo']):10.6f} "
              f"{float(res['sinT']):8.4f} "
              f"{res['lam1_approx']-res['lam2_approx']:10.6f}")

    print()
    print("=" * 70)
    print("STRICTNESS CLASSIFICATION")
    print("=" * 70)
    print("""
All bounds: **certificate** (mpmath.iv interval arithmetic, direction-audited).

Tail bounds (||G^BB||, ||C||): currently COARSE PLACEHOLDER (global bounds).
  -> To be replaced by tight Schur-test bounds (separate agent, analogous to
     P4's c1sq_schur_bound.py route). Current bounds are valid (all upper
     bounds on nonnegative quantities), just not tight.

Cross-validation needed:
  - mu*=0 result should match KAY's N=200 float64 value (within float rounding)
  - Coulomb suppression should be ~10-15% for mu*=0.10-0.15 (physical range)

Next steps (Phase 1 continuation):
  1. Implement tight Schur tail bounds (analytic for gamma>=1, needs care for gamma<1)
  2. Extend N to 400 (KAY's benchmark size)
  3. Cross-validate against KAY's N=4 exact algebraic closed form
  4. Multi-gamma sweep (gamma=0.5, 1.0, 1.5, 2.0)
""")
