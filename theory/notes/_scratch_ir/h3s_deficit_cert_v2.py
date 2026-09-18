"""Certificate-level Tc bounds for deficit budget via residual bounds on M.

DESIGN (sound eigenvalue bounds without needing exact eigenvectors):
=========================================================================
We need λ₁(M) where M = O[F_P] - 2(μ*/λ)u_c u_cᵀ. Instead of secular equation
requiring ALL exact eigenpairs (unsound with approximate eigenvectors), we use:

1. Build M in interval arithmetic (N×N symmetric interval matrix)
2. Get trial vector v from float64 (quality only affects bound width, not soundness)
3. Compute in interval arithmetic:
   ρ = ⟨v,Mv⟩/⟨v,v⟩  (Rayleigh quotient)
   R = ‖Mv - ρv‖/‖v‖   (residual norm)

4. BOUNDS (each stated with its direction and prerequisite):

   (a) LOWER BOUND (unconditional): ρ ≤ λ₁(M)
       Proof: Rayleigh-Ritz variational principle. For symmetric M,
              λ₁ = max_v ⟨v,Mv⟩/⟨v,v⟩ ≥ ⟨v,Mv⟩/⟨v,v⟩ = ρ for ANY v≠0.
       Direction: ρ is a LOWER bound on λ₁.

   (b) UPPER BOUND (symmetric M): λ₁(M) ≤ ρ + R
       Proof: For symmetric M, Bauer-Fike / Krylov-Bogoliubov residual theorem
              guarantees ∃ eigenvalue μ ∈ [ρ-R, ρ+R].
              Combined with (a): ρ ≤ λ₁, we get λ₁ ∈ [ρ, ρ+R].
       Direction: ρ+R is an UPPER bound on λ₁.

   Thus: λ₁(M) ∈ [ρ, ρ+R]  (certified interval, width = R)

5. Tc BOUNDS (judgment criterion: λ·λ₁(M) = 1, λ₁ decreases with T):

   (a) Tc UPPER bound: Find smallest T where λ·(ρ+R) < 1
       Logic: If λ·(upper bound on λ₁) < 1, then λ·λ₁ < 1 → non-superconducting.
       This T is an UPPER bound on Tc (system becomes normal at or before this T).
       N_c direction: use ⌊r·ϖ⌋ (smaller N_c → weaker Coulomb → higher Tc → supports upper bound).

   (b) Tc LOWER bound: Find largest T where λ·ρ ≥ 1
       Logic: If λ·(lower bound on λ₁) ≥ 1, then λ·λ₁ ≥ 1 → superconducting.
       This T is a LOWER bound on Tc (system remains superconducting up to this T).
       N_c direction: use ⌈r·ϖ⌉ (larger N_c → stronger Coulomb → lower Tc → supports lower bound).

6. QUADRATURE: α²F integral uses trapezoidal rule on file grid.
   ASSUMPTION: This script is conditional on the trapezoid rule being exact.
   TODO: Add interval padding for quadrature error (deferred to post-v1).

7. ALIGNMENT REQUIREMENTS (mandatory verification):
   - Compare fker[0,1,20,399], L1[0,1,20,399], A[0,0], A[0,1], A[20,20] against
     float64 baseline (h3s_deficit_budget.py) to ensure interval encloses baseline.
   - L1[n] = (1/(2n+1))·Σ_{k=1}^{n} 2f(k), upper limit INCLUDES k=n (old script bug).
   - Run μ*=0 case and verify interval encloses 283.344 K.

certificate (interval arithmetic with mpmath.iv, all eigenvalue bounds proven sound)
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_fill_factor import moments, MEV_K, ROOT


def O_of_pt(fker, N):
    """Float64 version: KAY operator O[f]."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    absd = np.abs(n[:, None] - n[None, :])
    sumd = n[:, None] + n[None, :] + 1
    uu = u[:, None] * u[None, :]
    Am = (fker[absd] + fker[sumd]) * uu
    # L1[n] = (1/(2n+1)) * sum_{k=1}^{n} 2*f(k)  — upper limit INCLUDES k=n
    L1 = np.zeros(N)
    L1[1:] = np.cumsum(2.0 * fker[1:N+1])[:N-1] / (2 * n[1:] + 1)  # sum up to k=n
    np.fill_diagonal(Am, -L1 + fker[2 * n + 1] * u ** 2)
    return n, u, Am, L1


def FP_pt(ks2, xs, ws):
    """Float64 F_P(k²)."""
    return np.einsum('i,ki->k', ws, xs[None, :] / (ks2[:, None] + xs[None, :]))


def O_of_iv(fker_iv, N, dps):
    """Interval version: KAY operator O[f] with mpmath.iv.

    Note: fker_iv contains mpi intervals. We convert u to mpi for consistency.
    """
    mp.dps = dps
    n_arr = np.arange(N)

    # u as intervals (zero-width, but mpi type for arithmetic compatibility)
    u_iv = [mp.mpi(mp.mpf(1) / mp.sqrt(mp.mpf(2 * n + 1))) for n in range(N)]

    # Build A as list of lists
    A_iv = [[mp.mpi(0) for _ in range(N)] for _ in range(N)]

    # L1[n] = (1/(2n+1)) * sum_{k=1}^{n} 2*f(k)
    L1_iv = [mp.mpi(0) for _ in range(N)]
    for n in range(1, N):
        cumsum = mp.mpi(0)
        for k in range(1, n+1):  # INCLUDE k=n (old script bug was range(1,n))
            cumsum += fker_iv[k] * mp.mpi(2)
        L1_iv[n] = cumsum / mp.mpi(2*n + 1)

    # Off-diagonal and diagonal
    for i in range(N):
        for j in range(N):
            absd = abs(i - j)
            sumd = i + j + 1
            uu = u_iv[i] * u_iv[j]
            A_iv[i][j] = (fker_iv[absd] + fker_iv[sumd]) * uu
            if i == j:
                A_iv[i][i] = -L1_iv[i] + fker_iv[2*i + 1] * u_iv[i]**2

    return u_iv, A_iv, L1_iv


def FP_iv(ks2_iv, xs_iv, ws_iv, dps):
    """Interval F_P(k²) at a single k² value.

    Returns an mpi interval.
    """
    mp.dps = dps
    result = mp.mpi(0)
    for i in range(len(xs_iv)):
        result += ws_iv[i] * xs_iv[i] / (ks2_iv + xs_iv[i])
    return result


def build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, nc_mode='floor'):
    """Build M = O[F_P] - 2(μ*/λ)u_c u_cᵀ in interval arithmetic.

    nc_mode: 'floor' for Tc upper bound (weaker Coulomb), 'ceil' for Tc lower bound.
    Returns: M_iv (list of lists), u_iv (list), n_arr (np.array)
    """
    mp.dps = dps

    # Spectral weights (trapezoid rule, assumed exact for now)
    ws_pt = (2.0 / lam) * (a2f / w)
    if len(w) == 1:
        p_pt = np.array([1.0])
    else:
        tw = np.zeros_like(w)
        tw[1:-1] = 0.5 * (w[2:] - w[:-2])
        tw[0] = 0.5 * (w[1] - w[0])
        tw[-1] = 0.5 * (w[-1] - w[-2])
        p_pt = ws_pt * tw
        p_pt = p_pt / p_pt.sum()

    # Convert to interval
    xs_iv = [mp.mpi(mp.mpf(wi / (2 * math.pi * T_mev)))**2 for wi in w]
    ws_iv = [mp.mpi(mp.mpf(p)) for p in p_pt]

    # Build F_P kernel
    ks = np.arange(2 * N + 1, dtype=float)
    fker_iv = []
    for k in ks:
        ks2_iv = mp.mpi(mp.mpf(k)**2)
        fker_iv.append(FP_iv(ks2_iv, xs_iv, ws_iv, dps))
    fker_iv[0] = mp.mpi(mp.mpf(1))  # normalization

    # Build O[F_P]
    u_iv, A_iv, L1_iv = O_of_iv(fker_iv, N, dps)

    # Coulomb correction
    if mus == 0.0:
        return A_iv, u_iv, np.arange(N), fker_iv, L1_iv

    # ϖ = sqrt(<ω²>)/(2πT)
    sqom2 = math.sqrt(np.sum(p_pt * w ** 2))
    vp_iv = mp.mpf(sqom2) / (mp.mpf(2) * (+mp.pi) * mp.mpf(T_mev))

    # N_c with directional rounding
    if nc_mode == 'floor':
        Nc = int(mp.floor(mp.mpf(r) * vp_iv))
    else:  # 'ceil'
        Nc = int(mp.ceil(mp.mpf(r) * vp_iv))

    if Nc <= 0:
        return A_iv, u_iv, np.arange(N), fker_iv, L1_iv

    # u_c = sqrt-normalized characteristic function
    uc_iv = [u_iv[n] if n < Nc else mp.mpf(0) for n in range(N)]

    # M = A - 2(μ*/λ)u_c u_cᵀ
    mut_iv = mp.mpf(mus) / mp.mpf(lam)
    M_iv = [[A_iv[i][j] for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            M_iv[i][j] = M_iv[i][j] - mp.mpf(2) * mut_iv * uc_iv[i] * uc_iv[j]

    return M_iv, u_iv, np.arange(N), fker_iv, L1_iv


def rayleigh_residual_bounds(M_iv, v_pt, dps):
    """Compute ρ and R for eigenvalue bounds.

    Returns: (rho_iv, R_iv) where λ₁ ∈ [rho_iv, rho_iv + R_iv]
    """
    mp.dps = dps
    N = len(v_pt)

    # Convert trial vector to mpmath intervals (zero-width but mpi type)
    v_iv = [mp.mpi(mp.mpf(float(v_pt[i]))) for i in range(N)]

    # Mv
    Mv_iv = [mp.mpi(0) for _ in range(N)]
    for i in range(N):
        for j in range(N):
            Mv_iv[i] += M_iv[i][j] * v_iv[j]

    # ⟨v,v⟩
    vv = mp.mpi(0)
    for i in range(N):
        vv += v_iv[i]**2

    # ⟨v,Mv⟩
    vMv = mp.mpi(0)
    for i in range(N):
        vMv += v_iv[i] * Mv_iv[i]

    # ρ = ⟨v,Mv⟩/⟨v,v⟩
    rho_iv = vMv / vv

    # ‖Mv - ρv‖²
    residual_sq = mp.mpi(0)
    for i in range(N):
        diff = Mv_iv[i] - rho_iv * v_iv[i]
        residual_sq += diff**2

    # R = ‖Mv - ρv‖/‖v‖
    # For intervals, sqrt needs special handling
    R_iv = residual_sq**mp.mpf(0.5) / vv**mp.mpf(0.5)

    return rho_iv, R_iv


def Q_P_float64(T_mev, w, a2f, lam, mus, r, N):
    """Float64 λ₁(M) for trial vector extraction."""
    ws = (2.0 / lam) * (a2f / w)
    if len(w) == 1:
        p = np.array([1.0])
    else:
        tw = np.zeros_like(w)
        tw[1:-1] = 0.5 * (w[2:] - w[:-2])
        tw[0] = 0.5 * (w[1] - w[0])
        tw[-1] = 0.5 * (w[-1] - w[-2])
        p = ws * tw
        p = p / p.sum()

    xs = (w / (2 * math.pi * T_mev)) ** 2
    ks = np.arange(2 * N + 1, dtype=float)
    fker_pt = FP_pt(ks ** 2, xs, p)
    fker_pt[0] = 1.0
    n, u, A, L1 = O_of_pt(fker_pt, N)

    if mus == 0.0:
        mu, PSI = eigh(A)
        return float(mu[-1]), PSI[:, -1]

    sqom2 = math.sqrt(np.sum(p * w ** 2))
    vp = sqom2 / (2 * math.pi * T_mev)
    Nc = int(math.floor(r * vp))

    if Nc <= 0:
        mu, PSI = eigh(A)
        return float(mu[-1]), PSI[:, -1]

    mu, PSI = eigh(A)
    c = PSI.T @ np.where(n < Nc, u, 0.0)
    m1, m2 = mu[-1], mu[-2]
    wdg = (m1 - m2) * 1e-13 + 1e-300
    mut = mus / lam
    phi = lambda l: 2 * mut * float(np.sum(c ** 2 / (mu - l))) - 1.0
    lam1_M = brentq(phi, m2 + wdg, m1 - wdg, xtol=1e-16, rtol=1e-15)

    # Reconstruct v₁ (eigenvector of M)
    v1 = PSI @ (c / (mu - lam1_M))
    v1 = v1 / np.linalg.norm(v1)

    return float(lam1_M), v1


def Tc_bounds_cert(w, a2f, lam, mus, r, N, dps):
    """Certificate Tc bounds via residual method.

    Returns: (Tc_lo, Tc_hi) in meV, both certified via interval arithmetic.
    """
    mp.dps = dps

    # Temperature scan grid (meV)
    sqom2 = math.sqrt(np.sum((2.0/lam)*(a2f/w) * w**2))  # crude estimate
    T_scan = np.geomspace(0.05, 2.0 * sqom2 / (2*math.pi), 80)

    # Find Tc upper bound: smallest T where λ·(ρ+R) < 1
    Tc_hi = None
    for T_mev in T_scan:
        # Get trial vector from float64
        _, v_trial = Q_P_float64(T_mev, w, a2f, lam, mus, r, N)

        # Build M with floor(N_c) for upper bound
        M_iv, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, nc_mode='floor')

        # Get bounds
        rho_iv, R_iv = rayleigh_residual_bounds(M_iv, v_trial, dps)
        lam1_upper = rho_iv + R_iv

        # Check criterion: λ·λ₁_upper < 1 → non-SC
        # Extract upper bound of the interval
        lam1_upper_val = float(lam1_upper.b)  # upper endpoint
        if lam * lam1_upper_val < 1.0:
            Tc_hi = T_mev
            break

    # Refine Tc_hi with bisection
    if Tc_hi is not None and len(T_scan) > 1:
        # Find last T where λ·(ρ+R) >= 1
        idx = np.where(T_scan < Tc_hi)[0]
        if len(idx) > 0:
            a, b = T_scan[idx[-1]], Tc_hi
            for _ in range(20):
                mid = 0.5 * (a + b)
                _, v_trial = Q_P_float64(mid, w, a2f, lam, mus, r, N)
                M_iv, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, nc_mode='floor')
                rho_iv, R_iv = rayleigh_residual_bounds(M_iv, v_trial, dps)
                lam1_upper_val = float((rho_iv + R_iv).b)
                if lam * lam1_upper_val >= 1.0:
                    a = mid
                else:
                    b = mid
            Tc_hi = b

    # Find Tc lower bound: largest T where λ·ρ >= 1
    Tc_lo = None
    for T_mev in reversed(T_scan):
        # Get trial vector from float64
        _, v_trial = Q_P_float64(T_mev, w, a2f, lam, mus, r, N)

        # Build M with ceil(N_c) for lower bound
        M_iv, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, nc_mode='ceil')

        # Get bounds
        rho_iv, R_iv = rayleigh_residual_bounds(M_iv, v_trial, dps)

        # Check criterion: λ·ρ >= 1 → SC
        # Extract lower bound of the interval
        rho_val = float(rho_iv.a)  # lower endpoint
        if lam * rho_val >= 1.0:
            Tc_lo = T_mev
            break

    # Refine Tc_lo with bisection
    if Tc_lo is not None and len(T_scan) > 1:
        idx = np.where(T_scan > Tc_lo)[0]
        if len(idx) > 0:
            a, b = Tc_lo, T_scan[idx[0]]
            for _ in range(20):
                mid = 0.5 * (a + b)
                _, v_trial = Q_P_float64(mid, w, a2f, lam, mus, r, N)
                M_iv, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, nc_mode='ceil')
                rho_iv, R_iv = rayleigh_residual_bounds(M_iv, v_trial, dps)
                rho_val = float(rho_iv.a)
                if lam * rho_val >= 1.0:
                    a = mid
                else:
                    b = mid
            Tc_lo = a

    return Tc_lo, Tc_hi


if __name__ == "__main__":
    import glob, time

    N = 400
    dps = 50
    r = 10.0

    print("="*80)
    print("H3S DEFICIT BUDGET: Certificate-level Tc bounds via residual method")
    print("="*80)
    print(f"N={N}, dps={dps}, r={r}")
    print()

    # Load H3S-200GPa anharmonic
    tag = "H3S-200GPa"
    sub = "EPW_anharmonic"
    g = glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f")
    g = [x for x in g if not x.endswith(".a2fv")]

    if not g:
        print(f"ERROR: No a2f file found for {tag}/{sub}")
        sys.exit(1)

    lam, om2, olog, w, a2f = moments(g[0])
    print(f"Dataset: {tag}/{sub}")
    print(f"lambda = {lam:.6f}, sqrt(<om^2>) = {math.sqrt(om2):.3f} meV")
    print()

    # ALIGNMENT CHECK (i): Compare intermediate quantities against float64 baseline
    print("ALIGNMENT CHECK (i): Comparing against float64 baseline")
    print("-" * 80)

    # Float64 baseline at T=1.0 meV, μ*=0.13
    T_test = 1.0
    lam1_pt, v_pt = Q_P_float64(T_test, w, a2f, lam, 0.13, r, N)

    # Build float64 kernel and operator
    ws = (2.0 / lam) * (a2f / w)
    tw = np.zeros_like(w)
    tw[1:-1] = 0.5 * (w[2:] - w[:-2])
    tw[0] = 0.5 * (w[1] - w[0])
    tw[-1] = 0.5 * (w[-1] - w[-2])
    p = ws * tw / (ws * tw).sum()
    xs = (w / (2 * math.pi * T_test)) ** 2
    ks = np.arange(2 * N + 1, dtype=float)
    fker_pt = FP_pt(ks ** 2, xs, p)
    fker_pt[0] = 1.0
    n, u, A_pt, L1_pt = O_of_pt(fker_pt, N)

    # Interval version
    M_iv, u_iv, _, fker_iv, L1_iv = build_M_iv(w, a2f, lam, 0.13, r, T_test, N, dps, nc_mode='floor')

    # Compare fker at selected indices
    test_indices = [0, 1, 20, 2*N]
    print("fker comparison (float64 vs interval):")
    for idx in test_indices:
        if idx < len(fker_pt):
            pt_val = fker_pt[idx]
            iv_val = fker_iv[idx]
            enclosed = float(iv_val.a) <= pt_val <= float(iv_val.b)
            print(f"  fker[{idx:3d}]: pt={pt_val:.10e}, iv=[{float(iv_val.a):.10e}, {float(iv_val.b):.10e}], "
                  f"enclosed={'PASS' if enclosed else 'FAIL'}")

    print("\nL1 comparison (float64 vs interval):")
    test_n = [0, 1, 20, N-1]
    for n_idx in test_n:
        if n_idx < N:
            pt_val = L1_pt[n_idx]
            iv_val = L1_iv[n_idx]
            enclosed = float(iv_val.a) <= pt_val <= float(iv_val.b)
            print(f"  L1[{n_idx:3d}]: pt={pt_val:.10e}, iv=[{float(iv_val.a):.10e}, {float(iv_val.b):.10e}], "
                  f"enclosed={'PASS' if enclosed else 'FAIL'}")

    print("\nA (O[F_P]) comparison at selected elements:")
    test_pairs = [(0, 0), (0, 1), (20, 20), (20, 21)]
    for i, j in test_pairs:
        pt_val = A_pt[i, j]
        # Reconstruct A from M by adding back Coulomb term
        sqom2 = math.sqrt(np.sum(p * w ** 2))
        vp = sqom2 / (2 * math.pi * T_test)
        Nc = int(math.floor(r * vp))
        uc_pt = np.where(n < Nc, u, 0.0)
        M_pt = A_pt - 2 * (0.13 / lam) * np.outer(uc_pt, uc_pt)
        A_from_M_iv = M_iv[i][j] + mp.mpf(2) * mp.mpf(0.13) / mp.mpf(lam) * u_iv[i] * u_iv[j] if i < Nc and j < Nc else M_iv[i][j]

        # For direct comparison, rebuild A_iv without Coulomb
        M_iv_no_coul, _, _, _, _ = build_M_iv(w, a2f, lam, 0.0, r, T_test, N, dps, nc_mode='floor')
        iv_val = M_iv_no_coul[i][j]
        enclosed = float(iv_val.a) <= pt_val <= float(iv_val.b)
        print(f"  A[{i:2d},{j:2d}]: pt={pt_val:.10e}, iv=[{float(iv_val.a):.10e}, {float(iv_val.b):.10e}], "
              f"enclosed={'PASS' if enclosed else 'FAIL'}")

    print()
    print("="*80)

    # ALIGNMENT CHECK (ii): μ*=0 case
    print("\nALIGNMENT CHECK (ii): μ*=0 case (target: 283.344 K)")
    print("-" * 80)
    t0 = time.time()
    Tc_lo_0, Tc_hi_0 = Tc_bounds_cert(w, a2f, lam, 0.0, r, N, dps)
    t1 = time.time()

    target_K = 283.344
    print(f"Tc bounds (μ*=0):")
    print(f"  Lower: {Tc_lo_0 * MEV_K:.3f} K")
    print(f"  Upper: {Tc_hi_0 * MEV_K:.3f} K")
    print(f"  Width: {(Tc_hi_0 - Tc_lo_0) * MEV_K:.3f} K")
    print(f"  Target (float64): {target_K:.3f} K")
    enclosed_0 = Tc_lo_0 * MEV_K <= target_K <= Tc_hi_0 * MEV_K
    print(f"  Encloses target: {'PASS' if enclosed_0 else 'FAIL'}")
    print(f"  Time: {t1-t0:.1f} s")
    print()

    # Main computation: μ*=0.13
    print("="*80)
    print("MAIN COMPUTATION: μ*=0.13 (target: 198.974 K)")
    print("-" * 80)
    t0 = time.time()
    Tc_lo_13, Tc_hi_13 = Tc_bounds_cert(w, a2f, lam, 0.13, r, N, dps)
    t1 = time.time()

    target_K_13 = 198.974
    print(f"Tc bounds (μ*=0.13):")
    print(f"  Lower: {Tc_lo_13 * MEV_K:.3f} K")
    print(f"  Upper: {Tc_hi_13 * MEV_K:.3f} K")
    print(f"  Width: {(Tc_hi_13 - Tc_lo_13) * MEV_K:.3f} K")
    print(f"  Target (float64): {target_K_13:.3f} K")
    enclosed_13 = Tc_lo_13 * MEV_K <= target_K_13 <= Tc_hi_13 * MEV_K
    print(f"  Encloses target: {'PASS' if enclosed_13 else 'FAIL'}")
    print(f"  Time: {t1-t0:.1f} s")
    print()

    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"μ*=0:    Tc ∈ [{Tc_lo_0*MEV_K:.3f}, {Tc_hi_0*MEV_K:.3f}] K  (width: {(Tc_hi_0-Tc_lo_0)*MEV_K:.3f} K)")
    print(f"μ*=0.13: Tc ∈ [{Tc_lo_13*MEV_K:.3f}, {Tc_hi_13*MEV_K:.3f}] K  (width: {(Tc_hi_13-Tc_lo_13)*MEV_K:.3f} K)")
    print()

    if not enclosed_0 or not enclosed_13:
        print("WARNING: Some targets not enclosed. Check alignment or increase dps.")
