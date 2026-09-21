"""Certificate-level Tc bounds for deficit budget via Temple inequality.

DESIGN (sound eigenvalue bounds via Temple + Frobenius):
=========================================================================
We need λ₁(M) where M = O[F_P] - 2(μ*/λ)u_c u_cᵀ.

LOWER BOUND (unconditional, sound):
  ρ ≤ λ₁(M)  where ρ = ⟨v,Mv⟩/⟨v,v⟩ for any trial vector v.
  Proof: Rayleigh-Ritz variational principle.

UPPER BOUND (Temple inequality chain):
  Step 1: λ₂(M) ≤ λ₂(A) via Weyl (M = A - rank-one PSD)
  Step 2: λ₂(A) ≤ λ₂(G) where G = A + diag(L1) (requires L1 ≥ 0, checked at runtime)
  Step 3: λ₂(G) ≤ β where β = sqrt(‖G‖_F² - ρ_G²) (requires ρ_G ≥ 0, checked)
          Proof: ρ_G ≥ 0 and ρ_G ≤ λ₁(G) ⟹ λ₁(G)² ≥ ρ_G²
                 λ₂(G)² ≤ Σ_{i≥2}λᵢ² = ‖G‖_F² - λ₁(G)² ≤ ‖G‖_F² - ρ_G²
          NOTE: Does NOT require G ⪰ 0 (G has min eig ≈ -0.238 in practice)
  Step 4: Temple: If ρ > β, then λ₁(M) ≤ ρ + R²/(ρ - β)
          where R = ‖Mv - ρv‖/‖v‖
  Step 5: Fallback when ρ ≤ β: λ₁(M) ≤ ‖M‖_F (trivial bound, always sound)

Tc BOUNDS (judgment: λ·λ₁(M) = 1, λ₁ decreases with T):
  - Tc UPPER: Find smallest T where λ·(upper bound on λ₁) < 1
    N_c direction: use ⌊r·ϖ⌋ (weaker Coulomb → higher Tc → supports upper bound)
  - Tc LOWER: Find largest T where λ·ρ ≥ 1
    N_c direction: use ⌈r·ϖ⌉ (stronger Coulomb → lower Tc → supports lower bound)

ALIGNMENT: Compare against EXTERNAL baseline (_main_indep_baseline_N.py) at matching N
  at T=1.0 meV: F[0,1,20,2N or 800], L1[0,1,20,N-1], A[0,0], A[0,1], A[20,20], A[20,21]

NEGATIVE CONTROL (--poison mode): Break L1 sum to k=n-1 (old bug), must FAIL L1[20], A[20,20]

certificate (interval arithmetic with mpmath.iv, all bounds proven sound)
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_fill_factor import moments, MEV_K, ROOT


# Baseline values by N (from _main_indep_baseline_N.py, T=1.0 meV)
BASELINES = {
    40: {
        'F': {0: 1.0000000000000000e+00, 1: 9.8926738719898688e-01,
              20: 4.7083503166965779e-01, 80: 7.7154780791224328e-02},
        'L1': {0: 0.0000000000000000e+00, 1: 6.5951159146599125e-01,
               20: 6.9152796470690037e-01, 39: 5.1792707943366445e-01},
        'A': {(0,0): 9.8926738719898688e-01, (0,1): 1.1290937977616782e+00,
              (20,20): -6.8615029568296637e-01, (20,21): 2.8643164353339162e-02}
    },
    100: {
        'F': {0: 1.0000000000000000e+00, 1: 9.8926738719898688e-01,
              20: 4.7083503166965779e-01, 200: 1.3788506251693295e-02},
        'L1': {0: 0.0000000000000000e+00, 1: 6.5951159146599125e-01,
               20: 6.9152796470690037e-01, 99: 2.7311602015739078e-01},
        'A': {(0,0): 9.8926738719898688e-01, (0,1): 1.1290937977616782e+00,
              (20,20): -6.8615029568296637e-01, (20,21): 2.8643164353339162e-02}
    },
    400: {
        'F': {0: 1.0000000000000000e+00, 1: 9.8926738719898688e-01,
              20: 4.7083503166965779e-01, 800: 8.8050658178847501e-04},
        'L1': {0: 0.0000000000000000e+00, 1: 6.5951159146599125e-01,
               20: 6.9152796470690037e-01, 399: 7.8271859984520367e-02},
        'A': {(0,0): 9.8926738719898688e-01, (0,1): 1.1290937977616782e+00,
              (20,20): -6.8615029568296637e-01, (20,21): 2.8643164353339162e-02}
    }
}


def O_of_pt(fker, N, poison=False):
    """Float64 version: KAY operator O[f].

    poison: if True, use buggy L1 sum k=1..n-1 (old bug for negative control)
    """
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    absd = np.abs(n[:, None] - n[None, :])
    sumd = n[:, None] + n[None, :] + 1
    uu = u[:, None] * u[None, :]
    Am = (fker[absd] + fker[sumd]) * uu

    # L1[n] = (1/(2n+1)) * sum_{k=1}^{n} 2*f(k)
    # CORRECT: upper limit INCLUDES k=n
    # POISON: upper limit is k=n-1 (old bug)
    L1 = np.zeros(N)
    if poison:
        # Bug: sum_{k=1}^{n-1} instead of sum_{k=1}^{n}
        for nn in range(1, N):
            L1[nn] = np.sum(2.0 * fker[1:nn]) / (2 * nn + 1)  # stops at nn-1
    else:
        # Correct: includes k=n
        L1[1:] = np.cumsum(2.0 * fker[1:N+1])[:N-1] / (2 * n[1:] + 1)

    np.fill_diagonal(Am, -L1 + fker[2 * n + 1] * u ** 2)
    return n, u, Am, L1


def FP_pt(ks2, xs, ws):
    """Float64 F_P(k²)."""
    return np.einsum('i,ki->k', ws, xs[None, :] / (ks2[:, None] + xs[None, :]))


def O_of_iv(fker_iv, N, dps, poison=False):
    """Interval version: KAY operator O[f] with mpmath.iv.

    poison: if True, use buggy L1 sum k=1..n-1 (old bug for negative control)
    """
    mp.dps = dps
    u_iv = [mp.mpi(mp.mpf(1) / mp.sqrt(mp.mpf(2 * n + 1))) for n in range(N)]
    A_iv = [[mp.mpi(0) for _ in range(N)] for _ in range(N)]

    # L1[n] = (1/(2n+1)) * sum_{k=1}^{n} 2*f(k)
    # CORRECT: INCLUDES k=n
    # POISON: stops at k=n-1
    L1_iv = [mp.mpi(0) for _ in range(N)]
    for n in range(1, N):
        cumsum = mp.mpi(0)
        upper = n if not poison else (n - 1)  # poison: exclude k=n
        if upper >= 1:
            for k in range(1, upper + 1):
                cumsum += fker_iv[k] * mp.mpi(2)
        L1_iv[n] = cumsum / mp.mpi(2*n + 1)

    # Build A
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
    """Interval F_P(k²) at a single k² value."""
    mp.dps = dps
    result = mp.mpi(0)
    for i in range(len(xs_iv)):
        result += ws_iv[i] * xs_iv[i] / (ks2_iv + xs_iv[i])
    return result


def rel_width(iv):
    """Relative width (b-a)/|mid| of an interval, or 0 if mid=0."""
    a, b = float(iv.a), float(iv.b)
    mid = 0.5 * (a + b)
    if abs(mid) < 1e-100:
        return 0.0
    return (b - a) / abs(mid)


def build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, nc_mode='floor', poison=False):
    """Build M = O[F_P] - 2(μ*/λ)u_c u_cᵀ in interval arithmetic.

    nc_mode: 'floor' for Tc upper bound (weaker Coulomb), 'ceil' for Tc lower bound.
    poison: if True, use buggy L1 for negative control.
    Returns: M_iv, u_iv, fker_iv, L1_iv
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
    u_iv, A_iv, L1_iv = O_of_iv(fker_iv, N, dps, poison=poison)

    # Coulomb correction
    if mus == 0.0:
        return A_iv, u_iv, fker_iv, L1_iv

    # ϖ = sqrt(<ω²>)/(2πT)
    sqom2 = math.sqrt(np.sum(p_pt * w ** 2))
    vp_iv = mp.mpf(sqom2) / (mp.mpf(2) * (+mp.pi) * mp.mpf(T_mev))

    # N_c with directional rounding
    if nc_mode == 'floor':
        Nc = int(mp.floor(mp.mpf(r) * vp_iv))
    else:  # 'ceil'
        Nc = int(mp.ceil(mp.mpf(r) * vp_iv))

    if Nc <= 0:
        return A_iv, u_iv, fker_iv, L1_iv

    # u_c = sqrt-normalized characteristic function
    uc_iv = [u_iv[n] if n < Nc else mp.mpf(0) for n in range(N)]

    # M = A - 2(μ*/λ)u_c u_cᵀ
    mut_iv = mp.mpf(mus) / mp.mpf(lam)
    M_iv = [[A_iv[i][j] for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            M_iv[i][j] = M_iv[i][j] - mp.mpf(2) * mut_iv * uc_iv[i] * uc_iv[j]

    return M_iv, u_iv, fker_iv, L1_iv


def rayleigh_residual_iv(M_iv, v_pt, dps):
    """Compute ρ and R in interval arithmetic.

    Returns: (rho_iv, R_iv)
    """
    mp.dps = dps
    N = len(v_pt)
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

    # ρ = ⟨v,Mv⟩/⟨v,v⟩  (DIVISION)
    rho_iv = vMv / vv

    # ‖Mv - ρv‖²
    residual_sq = mp.mpi(0)
    for i in range(N):
        diff = Mv_iv[i] - rho_iv * v_iv[i]
        residual_sq += diff**2

    # R = ‖Mv - ρv‖/‖v‖  (TWO SQRTS, ONE DIVISION)
    R_iv = residual_sq**mp.mpf(0.5) / vv**mp.mpf(0.5)

    return rho_iv, R_iv


def frobenius_norm_sq_iv(M_iv, dps):
    """Compute ‖M‖_F² in interval arithmetic."""
    mp.dps = dps
    result = mp.mpi(0)
    N = len(M_iv)
    for i in range(N):
        for j in range(N):
            result += M_iv[i][j]**2
    return result


def temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_pt, dps):
    """Compute Temple upper bound for λ₁(M) via the chain:

    Step 1: λ₂(M) ≤ λ₂(A) (Weyl, M = A - rank-one PSD)
    Step 2: λ₂(A) ≤ λ₂(G) where G = A + diag(L1) (requires L1 ≥ 0)
    Step 3: λ₂(G) ≤ β where β = sqrt(‖G‖_F² - ρ_G²) (requires ρ_G ≥ 0)
    Step 4: Temple: if ρ > β, λ₁(M) ≤ ρ + R²/(ρ - β)
    Step 5: Fallback: λ₁(M) ≤ ‖M‖_F

    Returns: (upper_bound_iv, status_dict, width_warnings)
    """
    mp.dps = dps
    N = len(M_iv)
    width_warnings = []

    # Check L1 ≥ 0 (prerequisite for Step 2)
    L1_nonneg = all(float(L1_iv[i].a) >= 0 for i in range(N))

    if not L1_nonneg:
        # Skip Temple chain, use fallback
        normF_sq = frobenius_norm_sq_iv(M_iv, dps)
        upper = normF_sq**mp.mpf(0.5)  # SQRT
        rw = rel_width(upper)
        if rw == 0.0:
            width_warnings.append("‖M‖_F has zero rel width after sqrt")
        return upper, {'L1_nonneg': False, 'method': 'fallback_normF'}, width_warnings

    # Build G = A + diag(L1)
    G_iv = [[A_iv[i][j] for j in range(N)] for i in range(N)]
    for i in range(N):
        G_iv[i][i] = G_iv[i][i] + L1_iv[i]

    # Compute ρ_G with trial vector (same v used for M)
    rho_G_iv, _ = rayleigh_residual_iv(G_iv, v_pt, dps)  # CONTAINS DIVISION
    rho_G_nonneg = float(rho_G_iv.a) >= 0

    rw = rel_width(rho_G_iv)
    if rw == 0.0:
        width_warnings.append("ρ_G has zero rel width after division")

    if not rho_G_nonneg:
        # Skip Temple, use fallback
        normF_sq = frobenius_norm_sq_iv(M_iv, dps)
        upper = normF_sq**mp.mpf(0.5)  # SQRT
        rw = rel_width(upper)
        if rw == 0.0:
            width_warnings.append("‖M‖_F has zero rel width after sqrt")
        return upper, {'L1_nonneg': True, 'rho_G_nonneg': False, 'method': 'fallback_normF'}, width_warnings

    # Compute β = sqrt(‖G‖_F² - ρ_G²)
    normG_sq = frobenius_norm_sq_iv(G_iv, dps)
    beta_sq = normG_sq - rho_G_iv**2
    beta_iv = beta_sq**mp.mpf(0.5)  # SQRT

    rw = rel_width(beta_iv)
    if rw == 0.0:
        width_warnings.append("β has zero rel width after sqrt")

    # Compute ρ and R for M
    rho_iv, R_iv = rayleigh_residual_iv(M_iv, v_pt, dps)  # CONTAINS DIVISION + SQRT

    rw_rho = rel_width(rho_iv)
    rw_R = rel_width(R_iv)
    if rw_rho == 0.0:
        width_warnings.append("ρ has zero rel width after division")
    if rw_R == 0.0:
        width_warnings.append("R has zero rel width after division+sqrt")

    # Check ρ > β (prerequisite for Temple)
    rho_gt_beta = float(rho_iv.a) > float(beta_iv.b)

    if rho_gt_beta:
        # Temple: λ₁(M) ≤ ρ + R²/(ρ - β)  (DIVISION)
        temple_upper = rho_iv + R_iv**2 / (rho_iv - beta_iv)
        rw = rel_width(temple_upper)
        if rw == 0.0:
            width_warnings.append("Temple upper has zero rel width after division")
        return temple_upper, {
            'L1_nonneg': True, 'rho_G_nonneg': True, 'rho_gt_beta': True,
            'method': 'temple', 'beta': float(beta_iv.b), 'rho': float(rho_iv.a)
        }, width_warnings
    else:
        # Fallback: ‖M‖_F
        normF_sq = frobenius_norm_sq_iv(M_iv, dps)
        fallback_upper = normF_sq**mp.mpf(0.5)  # SQRT
        rw = rel_width(fallback_upper)
        if rw == 0.0:
            width_warnings.append("‖M‖_F has zero rel width after sqrt (fallback path)")

        return fallback_upper, {
            'L1_nonneg': True, 'rho_G_nonneg': True, 'rho_gt_beta': False,
            'method': 'fallback_normF', 'beta': float(beta_iv.b), 'rho': float(rho_iv.a)
        }, width_warnings


def Q_P_float64(T_mev, w, a2f, lam, mus, r, N, poison=False):
    """Float64 λ₁(M) for trial vector extraction.

    poison: if True, use buggy L1 for negative control.
    """
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
    n, u, A, L1 = O_of_pt(fker_pt, N, poison=poison)

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

    # Reconstruct v₁
    v1 = PSI @ (c / (mu - lam1_M))
    v1 = v1 / np.linalg.norm(v1)

    return float(lam1_M), v1


def Tc_bounds_cert(w, a2f, lam, mus, r, N, dps, poison=False):
    """Certificate Tc bounds via Temple + Rayleigh lower bound.

    poison: if True, use buggy L1 for negative control.
    Returns: (Tc_lo, Tc_hi, diagnostics)
    """
    mp.dps = dps

    # Temperature scan grid (meV)
    sqom2 = math.sqrt(np.sum((2.0/lam)*(a2f/w) * w**2))
    T_scan = np.geomspace(0.05, 2.0 * sqom2 / (2*math.pi), 80)

    diagnostics = {'upper_checks': [], 'lower_checks': [], 'width_warnings': []}

    # Find Tc upper bound: smallest T where λ·(upper bound on λ₁) < 1
    Tc_hi = None
    for T_mev in T_scan:
        _, v_trial = Q_P_float64(T_mev, w, a2f, lam, mus, r, N, poison=poison)
        M_iv, u_iv, fker_iv, L1_iv = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, nc_mode='floor', poison=poison)

        # Rebuild A_iv for Temple chain
        A_iv, _, _, _ = build_M_iv(w, a2f, lam, 0.0, r, T_mev, N, dps, nc_mode='floor', poison=poison)

        lam1_upper, status, warns = temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_trial, dps)
        lam1_upper_val = float(lam1_upper.b)

        if warns:
            diagnostics['width_warnings'].extend([(T_mev, 'upper', w) for w in warns])

        diagnostics['upper_checks'].append({
            'T_mev': T_mev, 'lam1_upper': lam1_upper_val,
            'criterion': lam * lam1_upper_val, 'status': status
        })

        if lam * lam1_upper_val < 1.0:
            Tc_hi = T_mev
            break

    # Refine Tc_hi with bisection
    if Tc_hi is not None and len(T_scan) > 1:
        idx = np.where(T_scan < Tc_hi)[0]
        if len(idx) > 0:
            a, b = T_scan[idx[-1]], Tc_hi
            for _ in range(20):
                mid = 0.5 * (a + b)
                _, v_trial = Q_P_float64(mid, w, a2f, lam, mus, r, N, poison=poison)
                M_iv, u_iv, fker_iv, L1_iv = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, nc_mode='floor', poison=poison)
                A_iv, _, _, _ = build_M_iv(w, a2f, lam, 0.0, r, mid, N, dps, nc_mode='floor', poison=poison)
                lam1_upper, _, warns = temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_trial, dps)
                lam1_upper_val = float(lam1_upper.b)
                if warns:
                    diagnostics['width_warnings'].extend([(mid, 'upper_refine', w) for w in warns])
                if lam * lam1_upper_val >= 1.0:
                    a = mid
                else:
                    b = mid
            Tc_hi = b

    # Find Tc lower bound: largest T where λ·ρ ≥ 1
    Tc_lo = None
    for T_mev in reversed(T_scan):
        _, v_trial = Q_P_float64(T_mev, w, a2f, lam, mus, r, N, poison=poison)
        M_iv, _, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, nc_mode='ceil', poison=poison)

        rho_iv, R_iv = rayleigh_residual_iv(M_iv, v_trial, dps)
        rho_val = float(rho_iv.a)

        rw_rho = rel_width(rho_iv)
        rw_R = rel_width(R_iv)
        warns = []
        if rw_rho == 0.0:
            warns.append("ρ has zero rel width after division (lower bound)")
        if rw_R == 0.0:
            warns.append("R has zero rel width after division+sqrt (lower bound)")
        if warns:
            diagnostics['width_warnings'].extend([(T_mev, 'lower', w) for w in warns])

        diagnostics['lower_checks'].append({
            'T_mev': T_mev, 'rho': rho_val, 'criterion': lam * rho_val
        })

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
                _, v_trial = Q_P_float64(mid, w, a2f, lam, mus, r, N, poison=poison)
                M_iv, _, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, nc_mode='ceil', poison=poison)
                rho_iv, R_iv = rayleigh_residual_iv(M_iv, v_trial, dps)
                rho_val = float(rho_iv.a)
                if lam * rho_val >= 1.0:
                    a = mid
                else:
                    b = mid
            Tc_lo = a

    return Tc_lo, Tc_hi, diagnostics


if __name__ == "__main__":
    import glob, time

    # Parse command-line arguments
    poison_mode = '--poison' in sys.argv
    args = [a for a in sys.argv[1:] if a != '--poison']
    N = int(args[0]) if len(args) > 0 else 400
    dps = 50
    r = 10.0

    mode_label = "POISON MODE (L1 sum k=1..n-1, must FAIL)" if poison_mode else "NORMAL MODE"

    print("="*80)
    print(f"H3S DEFICIT BUDGET v3: Temple inequality + Frobenius β")
    print(f"{mode_label}")
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

    # ALIGNMENT CHECK: Compare against external baseline (N-matched)
    print("ALIGNMENT CHECK: Comparing against external baseline")
    print(f"  Script N={N}, baseline lookup...")
    print("-" * 80)

    if N not in BASELINES:
        print(f"  N={N} not in baseline table → SKIP all alignment checks")
        print(f"  (Available: N={list(BASELINES.keys())})")
        print()
    else:
        print(f"  Baseline N={N} found → proceeding with element comparison")
        print()

        T_test = 1.0
        M_iv, u_iv, fker_iv, L1_iv = build_M_iv(w, a2f, lam, 0.0, r, T_test, N, dps, nc_mode='floor', poison=poison_mode)

        baseline = BASELINES[N]

        # Compare F kernel
        print("F kernel comparison (baseline vs interval):")
        for idx in baseline['F'].keys():
            if idx < len(fker_iv):
                baseline_val = baseline['F'][idx]
                iv_val = fker_iv[idx]
                enclosed = float(iv_val.a) <= baseline_val <= float(iv_val.b)
                rw = rel_width(iv_val)
                print(f"  F[{idx:3d}]: baseline={baseline_val:.16e}, "
                      f"iv=[{float(iv_val.a):.16e}, {float(iv_val.b):.16e}], "
                      f"relwidth={rw:.2e}, enclosed={'PASS' if enclosed else 'FAIL'}")

        # Compare L1
        print("\nL1 comparison (baseline vs interval):")
        for idx in baseline['L1'].keys():
            if idx < len(L1_iv):
                baseline_val = baseline['L1'][idx]
                iv_val = L1_iv[idx]
                enclosed = float(iv_val.a) <= baseline_val <= float(iv_val.b)
                rw = rel_width(iv_val)
                print(f"  L1[{idx:3d}]: baseline={baseline_val:.16e}, "
                      f"iv=[{float(iv_val.a):.16e}, {float(iv_val.b):.16e}], "
                      f"relwidth={rw:.2e}, enclosed={'PASS' if enclosed else 'FAIL'}")

        # Compare A (which is M_iv when mus=0)
        print("\nA (O[F_P]) comparison at selected elements:")
        for (i, j) in baseline['A'].keys():
            baseline_val = baseline['A'][(i,j)]
            iv_val = M_iv[i][j]
            enclosed = float(iv_val.a) <= baseline_val <= float(iv_val.b)
            rw = rel_width(iv_val)
            print(f"  A[{i:2d},{j:2d}]: baseline={baseline_val:.16e}, "
                  f"iv=[{float(iv_val.a):.16e}, {float(iv_val.b):.16e}], "
                  f"relwidth={rw:.2e}, enclosed={'PASS' if enclosed else 'FAIL'}")

        print()
        print("="*80)

    # μ*=0 case
    print("\nCOMPUTATION: μ*=0")
    print("-" * 80)
    t0 = time.time()
    Tc_lo_0, Tc_hi_0, diag_0 = Tc_bounds_cert(w, a2f, lam, 0.0, r, N, dps, poison=poison_mode)
    t1 = time.time()

    print(f"Tc bounds (μ*=0):")
    print(f"  Lower: {Tc_lo_0 * MEV_K:.3f} K  ({Tc_lo_0:.6f} meV)")
    print(f"  Upper: {Tc_hi_0 * MEV_K:.3f} K  ({Tc_hi_0:.6f} meV)")
    print(f"  Width: {(Tc_hi_0 - Tc_lo_0) * MEV_K:.3f} K")
    print(f"  Time: {t1-t0:.1f} s")

    # Prerequisite checks summary
    upper_methods = [c['status']['method'] for c in diag_0['upper_checks']]
    temple_count = sum(1 for m in upper_methods if m == 'temple')
    print(f"  Upper bound: {temple_count}/{len(upper_methods)} used Temple, "
          f"{len(upper_methods)-temple_count} used fallback")

    # Width warnings summary
    if diag_0['width_warnings']:
        print(f"\n  WIDTH WARNINGS ({len(diag_0['width_warnings'])} occurrences):")
        # Group by warning message
        warn_groups = {}
        for T, stage, msg in diag_0['width_warnings']:
            warn_groups.setdefault(msg, []).append((T, stage))
        for msg, occurrences in warn_groups.items():
            print(f"    {msg}: {len(occurrences)} times")
            if len(occurrences) <= 3:
                for T, stage in occurrences:
                    print(f"      T={T:.2f} meV ({stage})")
    print()

    # μ*=0.13 case
    print("="*80)
    print("COMPUTATION: μ*=0.13")
    print("-" * 80)
    t0 = time.time()
    Tc_lo_13, Tc_hi_13, diag_13 = Tc_bounds_cert(w, a2f, lam, 0.13, r, N, dps, poison=poison_mode)
    t1 = time.time()

    print(f"Tc bounds (μ*=0.13):")
    print(f"  Lower: {Tc_lo_13 * MEV_K:.3f} K  ({Tc_lo_13:.6f} meV)")
    print(f"  Upper: {Tc_hi_13 * MEV_K:.3f} K  ({Tc_hi_13:.6f} meV)")
    print(f"  Width: {(Tc_hi_13 - Tc_lo_13) * MEV_K:.3f} K")
    print(f"  Deficit: {100 * (1 - Tc_lo_13 / Tc_lo_0):.2f}% to {100 * (1 - Tc_hi_13 / Tc_hi_0):.2f}%")
    print(f"  Time: {t1-t0:.1f} s")

    # Prerequisite checks summary
    upper_methods = [c['status']['method'] for c in diag_13['upper_checks']]
    temple_count = sum(1 for m in upper_methods if m == 'temple')
    print(f"  Upper bound: {temple_count}/{len(upper_methods)} used Temple, "
          f"{len(upper_methods)-temple_count} used fallback")

    # Show where Temple failed (if any)
    temple_failures = [c for c in diag_13['upper_checks'] if c['status']['method'] != 'temple']
    if temple_failures:
        print(f"\n  Temple prerequisite failures:")
        for c in temple_failures[:5]:  # Show first 5
            st = c['status']
            print(f"    T={c['T_mev']:.2f} meV: ", end='')
            if not st.get('L1_nonneg', True):
                print("L1 not nonneg")
            elif not st.get('rho_G_nonneg', True):
                print("ρ_G < 0")
            elif not st.get('rho_gt_beta', True):
                print(f"ρ ≤ β (ρ={st.get('rho',0):.4f}, β={st.get('beta',0):.4f})")

    # Width warnings summary
    if diag_13['width_warnings']:
        print(f"\n  WIDTH WARNINGS ({len(diag_13['width_warnings'])} occurrences):")
        warn_groups = {}
        for T, stage, msg in diag_13['width_warnings']:
            warn_groups.setdefault(msg, []).append((T, stage))
        for msg, occurrences in warn_groups.items():
            print(f"    {msg}: {len(occurrences)} times")
            if len(occurrences) <= 3:
                for T, stage in occurrences:
                    print(f"      T={T:.2f} meV ({stage})")
    print()

    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"μ*=0:    Tc ∈ [{Tc_lo_0*MEV_K:.3f}, {Tc_hi_0*MEV_K:.3f}] K  (width: {(Tc_hi_0-Tc_lo_0)*MEV_K:.3f} K)")
    print(f"μ*=0.13: Tc ∈ [{Tc_lo_13*MEV_K:.3f}, {Tc_hi_13*MEV_K:.3f}] K  (width: {(Tc_hi_13-Tc_lo_13)*MEV_K:.3f} K)")
    print(f"Deficit: {100*(1-Tc_lo_13/Tc_lo_0):.2f}% to {100*(1-Tc_hi_13/Tc_hi_0):.2f}%")
    print()
