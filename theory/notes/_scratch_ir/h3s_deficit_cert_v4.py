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
  Step 4: Temple: If ρ > β STRICTLY, then λ₁(M) ≤ ρ + R²/(ρ - β)
          where R = ‖Mv - ρv‖/‖v‖
  Step 5: Fallback when ρ ≤ β: λ₁(M) ≤ ‖M‖_F (trivial bound, always sound)
          NO SLACK, NO MIN — use Temple alone when premise holds, else fallback alone.

N_c TREATMENT (F4 fix):
  N_c = floor(r·ϖ) is the MODEL DEFINITION, not an uncertain quantity to bracket.
  Compute Nc_lo = floor(r·ϖ_iv.a), Nc_hi = floor(r·ϖ_iv.b).
  If Nc_lo == Nc_hi: use that single N_c on BOTH legs (normal case).
  Else: evaluate at every N_c in [Nc_lo, Nc_hi], take min/max over endpoints.
  Print Nc_lo, Nc_hi, and which branch on every leg.

Tc BOUNDS (judgment: λ·λ₁(M) = 1, λ₁ decreases with T):
  - Tc UPPER: Find smallest T where λ·(upper bound on λ₁) < 1
  - Tc LOWER: Find largest T where λ·ρ ≥ 1

ALIGNMENT: Compare against EXTERNAL baseline (_main_indep_baseline_N.py) at matching N
  at T=1.0 meV: F[0,1,20,2N or 800], L1[0,1,20,N-1], A[0,0], A[0,1], A[20,20], A[20,21]

NEGATIVE CONTROL (--poison mode): Break L1 sum to k=n-1 (old bug), must FAIL L1[20], A[20,20]

FIXES FROM V3:
  F1: Strict rho > beta check, no slack, no min (v3 line 288-291)
  F2: Set mp.mp.dps AND mp.iv.dps (all dps assignments)
  F3: All arithmetic in interval form: u_n, xs, vp, mut, lam, p_pt
  F4: N_c is model definition, not bracketing parameter; handle straddle explicitly

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
    # F2: set BOTH contexts
    mp.mp.dps = dps
    mp.iv.dps = dps

    # F3: u_n in interval arithmetic (no point mpf injection)
    u_iv = [mp.mpi(1) / mp.sqrt(mp.mpi(2 * n + 1)) for n in range(N)]
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
    mp.mp.dps = dps
    mp.iv.dps = dps
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


def build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=None, poison=False):
    """Build M = O[F_P] - 2(μ*/λ)u_c u_cᵀ in interval arithmetic.

    Nc_fixed: if provided, use this N_c value (for straddle case); otherwise compute from vp_iv.
    poison: if True, use buggy L1 for negative control.
    Returns: M_iv, u_iv, fker_iv, L1_iv, Nc_used, vp_iv
    """
    mp.mp.dps = dps
    mp.iv.dps = dps

    # F3: Trapezoid rule and lambda in interval arithmetic
    # Build weights as intervals (file data w, a2f treated as exact points)
    ws_iv_raw = [(mp.mpi(2) / mp.mpi(lam)) * (mp.mpi(a2f[i]) / mp.mpi(w[i])) for i in range(len(w))]

    if len(w) == 1:
        p_iv = [mp.mpi(1)]
    else:
        # Trapezoid weights in interval arithmetic
        tw_iv = []
        for i in range(len(w)):
            if i == 0:
                tw_iv.append((mp.mpi(w[1]) - mp.mpi(w[0])) / mp.mpi(2))
            elif i == len(w) - 1:
                tw_iv.append((mp.mpi(w[-1]) - mp.mpi(w[-2])) / mp.mpi(2))
            else:
                tw_iv.append((mp.mpi(w[i+1]) - mp.mpi(w[i-1])) / mp.mpi(2))

        # p = ws * tw / sum(ws * tw)
        wstw_iv = [ws_iv_raw[i] * tw_iv[i] for i in range(len(w))]
        sum_wstw = mp.mpi(0)
        for x in wstw_iv:
            sum_wstw += x
        p_iv = [x / sum_wstw for x in wstw_iv]

    # F3: xs in interval arithmetic (no float64 division)
    xs_iv = [(mp.mpi(w[i]) / (mp.mpi(2) * mp.iv.pi * mp.mpi(T_mev)))**2 for i in range(len(w))]

    # Build F_P kernel
    ks = np.arange(2 * N + 1, dtype=float)
    fker_iv = []
    for k in ks:
        ks2_iv = mp.mpi(k)**2
        fker_iv.append(FP_iv(ks2_iv, xs_iv, p_iv, dps))
    fker_iv[0] = mp.mpi(1)  # normalization

    # Build O[F_P]
    u_iv, A_iv, L1_iv = O_of_iv(fker_iv, N, dps, poison=poison)

    # Coulomb correction
    if mus == 0.0:
        return A_iv, u_iv, fker_iv, L1_iv, 0, None

    # F3: ϖ = sqrt(<ω²>)/(2πT) in interval arithmetic
    # Compute <ω²> from p_iv and w
    om2_iv = mp.mpi(0)
    for i in range(len(w)):
        om2_iv += p_iv[i] * mp.mpi(w[i])**2
    sqom2_iv = om2_iv**(mp.mpi(1)/mp.mpi(2))
    vp_iv = sqom2_iv / (mp.mpi(2) * mp.iv.pi * mp.mpi(T_mev))

    # F4: N_c is model definition, not bracketing parameter
    if Nc_fixed is not None:
        Nc = Nc_fixed
    else:
        Nc = int(mp.floor(mp.mpf(r) * vp_iv.mid))  # use midpoint for single-value case

    if Nc <= 0:
        return A_iv, u_iv, fker_iv, L1_iv, 0, vp_iv

    # u_c = sqrt-normalized characteristic function
    uc_iv = [u_iv[n] if n < Nc else mp.mpi(0) for n in range(N)]

    # F3: mut in interval arithmetic
    mut_iv = mp.mpi(mus) / mp.mpi(lam)
    M_iv = [[A_iv[i][j] for j in range(N)] for i in range(N)]
    for i in range(N):
        for j in range(N):
            M_iv[i][j] = M_iv[i][j] - mp.mpi(2) * mut_iv * uc_iv[i] * uc_iv[j]

    return M_iv, u_iv, fker_iv, L1_iv, Nc, vp_iv


def rayleigh_residual_iv(M_iv, v_pt, dps):
    """Compute ρ and R in interval arithmetic.

    Returns: (rho_iv, R_iv)
    """
    mp.mp.dps = dps
    mp.iv.dps = dps
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
    R_iv = residual_sq**(mp.mpi(1)/mp.mpi(2)) / vv**(mp.mpi(1)/mp.mpi(2))

    return rho_iv, R_iv


def frobenius_norm_sq_iv(M_iv, dps):
    """Compute ‖M‖_F² in interval arithmetic."""
    mp.mp.dps = dps
    mp.iv.dps = dps
    result = mp.mpi(0)
    N = len(M_iv)
    for i in range(N):
        for j in range(N):
            result += M_iv[i][j]**2
    return result


def temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_pt, dps):
    """Compute Temple upper bound for λ₁(M) via the chain.

    F1 FIX: No slack, no min. Strict rho > beta → Temple alone, else fallback alone.

    Returns: (upper_bound_iv, status_dict, width_report)
    """
    mp.mp.dps = dps
    mp.iv.dps = dps
    N = len(M_iv)
    width_report = {}

    # Check L1 ≥ 0 (prerequisite for Step 2)
    L1_nonneg = all(float(L1_iv[i].a) >= 0 for i in range(N))

    if not L1_nonneg:
        normF_sq = frobenius_norm_sq_iv(M_iv, dps)
        upper = normF_sq**(mp.mpi(1)/mp.mpi(2))
        width_report['normF'] = rel_width(upper)
        return upper, {'L1_nonneg': False, 'method': 'fallback_normF'}, width_report

    # Build G = A + diag(L1)
    G_iv = [[A_iv[i][j] for j in range(N)] for i in range(N)]
    for i in range(N):
        G_iv[i][i] = G_iv[i][i] + L1_iv[i]

    # Compute ρ_G
    rho_G_iv, _ = rayleigh_residual_iv(G_iv, v_pt, dps)
    rho_G_nonneg = float(rho_G_iv.a) >= 0

    width_report['rho_G'] = rel_width(rho_G_iv)

    if not rho_G_nonneg:
        normF_sq = frobenius_norm_sq_iv(M_iv, dps)
        upper = normF_sq**(mp.mpi(1)/mp.mpi(2))
        width_report['normF'] = rel_width(upper)
        return upper, {'L1_nonneg': True, 'rho_G_nonneg': False, 'method': 'fallback_normF'}, width_report

    # Compute β = sqrt(‖G‖_F² - ρ_G²)
    normG_sq = frobenius_norm_sq_iv(G_iv, dps)
    beta_sq = normG_sq - rho_G_iv**2
    beta_iv = beta_sq**(mp.mpi(1)/mp.mpi(2))

    width_report['beta'] = rel_width(beta_iv)

    # Compute ρ and R for M
    rho_iv, R_iv = rayleigh_residual_iv(M_iv, v_pt, dps)

    width_report['rho'] = rel_width(rho_iv)
    width_report['R'] = rel_width(R_iv)

    # F1 FIX: STRICT check ρ > β (no slack)
    rho_gt_beta = float(rho_iv.a) > float(beta_iv.b)

    if rho_gt_beta:
        # Temple: λ₁(M) ≤ ρ + R²/(ρ - β)  — use ALONE
        temple_upper = rho_iv + R_iv**2 / (rho_iv - beta_iv)
        width_report['temple_upper'] = rel_width(temple_upper)
        return temple_upper, {
            'L1_nonneg': True, 'rho_G_nonneg': True, 'rho_gt_beta': True,
            'method': 'temple', 'beta': float(beta_iv.b), 'rho': float(rho_iv.a)
        }, width_report
    else:
        # Fallback: ‖M‖_F — use ALONE
        normF_sq = frobenius_norm_sq_iv(M_iv, dps)
        fallback_upper = normF_sq**(mp.mpi(1)/mp.mpi(2))
        width_report['fallback_normF'] = rel_width(fallback_upper)
        return fallback_upper, {
            'L1_nonneg': True, 'rho_G_nonneg': True, 'rho_gt_beta': False,
            'method': 'fallback_normF', 'beta': float(beta_iv.b), 'rho': float(rho_iv.a)
        }, width_report


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

    F4: N_c is model definition; handle straddle case explicitly.

    poison: if True, use buggy L1 for negative control.
    Returns: (Tc_lo, Tc_hi, diagnostics)
    """
    mp.mp.dps = dps
    mp.iv.dps = dps

    # Temperature scan grid (meV)
    sqom2 = math.sqrt(np.sum((2.0/lam)*(a2f/w) * w**2))
    T_scan = np.geomspace(0.05, 2.0 * sqom2 / (2*math.pi), 80)

    diagnostics = {
        'upper_checks': [], 'lower_checks': [],
        'width_samples': [],
        'premise_failures': [],
        'nc_info': []  # Track N_c on each leg
    }

    # Find Tc upper bound: smallest T where λ·(upper bound on λ₁) < 1
    Tc_hi = None
    for T_mev in T_scan:
        _, v_trial = Q_P_float64(T_mev, w, a2f, lam, mus, r, N, poison=poison)

        # F4: Check N_c straddle (skip for mus=0, no Coulomb term)
        if mus == 0.0:
            M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=None, poison=poison)
            A_iv = M_iv  # A = M when mus=0
            lam1_upper, status, widths = temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_trial, dps)
            lam1_upper_val = float(lam1_upper.b)
            diagnostics['nc_info'].append({'T_mev': T_mev, 'leg': 'upper', 'Nc_lo': 0, 'Nc_hi': 0, 'straddle': False, 'Nc_used': 0})
        else:
            M_iv_test, _, _, _, _, vp_iv = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=None, poison=poison)
            Nc_lo = int(mp.floor(mp.mpf(r) * vp_iv.a))
            Nc_hi = int(mp.floor(mp.mpf(r) * vp_iv.b))

            if Nc_lo == Nc_hi:
                # Normal case: single N_c
                Nc_used = Nc_lo
                M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=Nc_used, poison=poison)
                A_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, 0.0, r, T_mev, N, dps, Nc_fixed=Nc_used, poison=poison)
                lam1_upper, status, widths = temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_trial, dps)
                lam1_upper_val = float(lam1_upper.b)
                diagnostics['nc_info'].append({'T_mev': T_mev, 'leg': 'upper', 'Nc_lo': Nc_lo, 'Nc_hi': Nc_hi, 'straddle': False, 'Nc_used': Nc_used})
            else:
                # Straddle case: evaluate at all N_c in [Nc_lo, Nc_hi], take max
                lam1_upper_val = -float('inf')
                for Nc_try in range(Nc_lo, Nc_hi + 1):
                    M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=Nc_try, poison=poison)
                    A_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, 0.0, r, T_mev, N, dps, Nc_fixed=Nc_try, poison=poison)
                    lam1_upper_try, status_try, widths_try = temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_trial, dps)
                    lam1_upper_val = max(lam1_upper_val, float(lam1_upper_try.b))
                    if Nc_try == Nc_lo:
                        status, widths = status_try, widths_try  # use first for reporting
                diagnostics['nc_info'].append({'T_mev': T_mev, 'leg': 'upper', 'Nc_lo': Nc_lo, 'Nc_hi': Nc_hi, 'straddle': True, 'Nc_range': f'[{Nc_lo}, {Nc_hi}]'})

        diagnostics['upper_checks'].append({
            'T_mev': T_mev, 'lam1_upper': lam1_upper_val,
            'criterion': lam * lam1_upper_val, 'status': status
        })

        if not status.get('rho_gt_beta', False) and status.get('L1_nonneg', False) and status.get('rho_G_nonneg', False):
            diagnostics['premise_failures'].append({
                'T_mev': T_mev, 'rho': status.get('rho', 0), 'beta': status.get('beta', 0),
                'method': status['method']
            })

        if len(diagnostics['width_samples']) < 3:
            diagnostics['width_samples'].append({'T_mev': T_mev, 'stage': 'upper', 'widths': widths})

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
                if mus == 0.0:
                    M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, Nc_fixed=None, poison=poison)
                    A_iv = M_iv
                    lam1_upper, _, _ = temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_trial, dps)
                    lam1_upper_val = float(lam1_upper.b)
                else:
                    M_iv_test, _, _, _, _, vp_iv = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, Nc_fixed=None, poison=poison)
                    Nc_lo = int(mp.floor(mp.mpf(r) * vp_iv.a))
                    Nc_hi_val = int(mp.floor(mp.mpf(r) * vp_iv.b))
                    if Nc_lo == Nc_hi_val:
                        M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, Nc_fixed=Nc_lo, poison=poison)
                        A_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, 0.0, r, mid, N, dps, Nc_fixed=Nc_lo, poison=poison)
                        lam1_upper, _, _ = temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_trial, dps)
                        lam1_upper_val = float(lam1_upper.b)
                    else:
                        lam1_upper_val = -float('inf')
                        for Nc_try in range(Nc_lo, Nc_hi_val + 1):
                            M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, Nc_fixed=Nc_try, poison=poison)
                            A_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, 0.0, r, mid, N, dps, Nc_fixed=Nc_try, poison=poison)
                            lam1_upper_try, _, _ = temple_upper_bound_iv(M_iv, A_iv, L1_iv, v_trial, dps)
                            lam1_upper_val = max(lam1_upper_val, float(lam1_upper_try.b))
                if lam * lam1_upper_val >= 1.0:
                    a = mid
                else:
                    b = mid
            Tc_hi = b

    # Find Tc lower bound: largest T where λ·ρ ≥ 1
    Tc_lo = None
    for T_mev in reversed(T_scan):
        _, v_trial = Q_P_float64(T_mev, w, a2f, lam, mus, r, N, poison=poison)

        # F4: Check N_c straddle (skip for mus=0, no Coulomb term)
        if mus == 0.0:
            M_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=None, poison=poison)
            rho_iv, R_iv = rayleigh_residual_iv(M_iv, v_trial, dps)
            rho_val = float(rho_iv.a)
            diagnostics['nc_info'].append({'T_mev': T_mev, 'leg': 'lower', 'Nc_lo': 0, 'Nc_hi': 0, 'straddle': False, 'Nc_used': 0})
        else:
            M_iv_test, _, _, _, _, vp_iv = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=None, poison=poison)
            Nc_lo = int(mp.floor(mp.mpf(r) * vp_iv.a))
            Nc_hi = int(mp.floor(mp.mpf(r) * vp_iv.b))

            if Nc_lo == Nc_hi:
                Nc_used = Nc_lo
                M_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=Nc_used, poison=poison)
                rho_iv, R_iv = rayleigh_residual_iv(M_iv, v_trial, dps)
                rho_val = float(rho_iv.a)
                diagnostics['nc_info'].append({'T_mev': T_mev, 'leg': 'lower', 'Nc_lo': Nc_lo, 'Nc_hi': Nc_hi, 'straddle': False, 'Nc_used': Nc_used})
            else:
                rho_val = float('inf')
                for Nc_try in range(Nc_lo, Nc_hi + 1):
                    M_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, T_mev, N, dps, Nc_fixed=Nc_try, poison=poison)
                    rho_iv_try, R_iv_try = rayleigh_residual_iv(M_iv, v_trial, dps)
                    rho_val = min(rho_val, float(rho_iv_try.a))
                    if Nc_try == Nc_lo:
                        rho_iv, R_iv = rho_iv_try, R_iv_try
                diagnostics['nc_info'].append({'T_mev': T_mev, 'leg': 'lower', 'Nc_lo': Nc_lo, 'Nc_hi': Nc_hi, 'straddle': True, 'Nc_range': f'[{Nc_lo}, {Nc_hi}]'})

        rw_rho = rel_width(rho_iv)
        rw_R = rel_width(R_iv)

        diagnostics['lower_checks'].append({
            'T_mev': T_mev, 'rho': rho_val, 'criterion': lam * rho_val
        })

        if len([s for s in diagnostics['width_samples'] if s['stage'] == 'lower']) < 3:
            diagnostics['width_samples'].append({
                'T_mev': T_mev, 'stage': 'lower',
                'widths': {'rho': rw_rho, 'R': rw_R}
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
                if mus == 0.0:
                    M_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, Nc_fixed=None, poison=poison)
                    rho_iv, R_iv = rayleigh_residual_iv(M_iv, v_trial, dps)
                    rho_val = float(rho_iv.a)
                else:
                    M_iv_test, _, _, _, _, vp_iv = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, Nc_fixed=None, poison=poison)
                    Nc_lo_val = int(mp.floor(mp.mpf(r) * vp_iv.a))
                    Nc_hi_val = int(mp.floor(mp.mpf(r) * vp_iv.b))
                    if Nc_lo_val == Nc_hi_val:
                        M_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, Nc_fixed=Nc_lo_val, poison=poison)
                        rho_iv, R_iv = rayleigh_residual_iv(M_iv, v_trial, dps)
                        rho_val = float(rho_iv.a)
                    else:
                        rho_val = float('inf')
                        for Nc_try in range(Nc_lo_val, Nc_hi_val + 1):
                            M_iv, _, _, _, _, _ = build_M_iv(w, a2f, lam, mus, r, mid, N, dps, Nc_fixed=Nc_try, poison=poison)
                            rho_iv_try, _ = rayleigh_residual_iv(M_iv, v_trial, dps)
                            rho_val = min(rho_val, float(rho_iv_try.a))
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
    print(f"H3S DEFICIT BUDGET v4: Temple + Frobenius (F1/F2/F3/F4 fixes)")
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

    # ALIGNMENT CHECK: Compare against external baseline (always run)
    print("ALIGNMENT CHECK: Comparing against external baseline")
    print(f"  Script N={N}")
    print("-" * 80)

    if N not in BASELINES:
        print(f"  N={N} not in baseline table → SKIP all alignment checks")
        print(f"  (Available: N={list(BASELINES.keys())})")
        print()
    else:
        print(f"  Baseline N={N} found → proceeding with element comparison")
        print()

        T_test = 1.0
        M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(w, a2f, lam, 0.0, r, T_test, N, dps, Nc_fixed=None, poison=poison_mode)

        baseline = BASELINES[N]

        # Compare F kernel
        print("F kernel comparison (baseline vs interval):")
        for idx in baseline['F'].keys():
            if idx < len(fker_iv):
                baseline_val = baseline['F'][idx]
                iv_val = fker_iv[idx]
                a_val, b_val = float(iv_val.a), float(iv_val.b)
                enclosed = a_val <= baseline_val <= b_val
                rw = rel_width(iv_val)
                print(f"  F[{idx:3d}]: baseline={baseline_val:.16e}, "
                      f"iv=[{a_val:.16e}, {b_val:.16e}], "
                      f"relwidth={rw:.2e}, enclosed={'PASS' if enclosed else 'FAIL'}")

        # Compare L1
        print("\nL1 comparison (baseline vs interval):")
        for idx in baseline['L1'].keys():
            if idx < len(L1_iv):
                baseline_val = baseline['L1'][idx]
                iv_val = L1_iv[idx]
                a_val, b_val = float(iv_val.a), float(iv_val.b)
                enclosed = a_val <= baseline_val <= b_val
                rw = rel_width(iv_val)
                print(f"  L1[{idx:3d}]: baseline={baseline_val:.16e}, "
                      f"iv=[{a_val:.16e}, {b_val:.16e}], "
                      f"relwidth={rw:.2e}, enclosed={'PASS' if enclosed else 'FAIL'}")

        # Compare A
        print("\nA (O[F_P]) comparison at selected elements:")
        for (i, j) in baseline['A'].keys():
            baseline_val = baseline['A'][(i,j)]
            iv_val = M_iv[i][j]
            a_val, b_val = float(iv_val.a), float(iv_val.b)
            enclosed = a_val <= baseline_val <= b_val
            rw = rel_width(iv_val)
            print(f"  A[{i:2d},{j:2d}]: baseline={baseline_val:.16e}, "
                  f"iv=[{a_val:.16e}, {b_val:.16e}], "
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

    upper_methods = [c['status']['method'] for c in diag_0['upper_checks']]
    temple_count = sum(1 for m in upper_methods if m == 'temple')
    print(f"  Upper bound: {temple_count}/{len(upper_methods)} used Temple, "
          f"{len(upper_methods)-temple_count} used fallback")

    if diag_0['width_samples']:
        print(f"\n  WIDTH SAMPLES (first few temperatures):")
        for sample in diag_0['width_samples'][:3]:
            print(f"    T={sample['T_mev']:.2f} meV ({sample['stage']}):")
            for key, val in sample['widths'].items():
                print(f"      {key}: {val:.2e}")

    if diag_0['nc_info']:
        print(f"\n  N_c INFO (first few):")
        for info in diag_0['nc_info'][:3]:
            if info['straddle']:
                print(f"    T={info['T_mev']:.2f} meV ({info['leg']}): Nc_lo={info['Nc_lo']}, Nc_hi={info['Nc_hi']}, STRADDLE, range={info['Nc_range']}")
            else:
                print(f"    T={info['T_mev']:.2f} meV ({info['leg']}): Nc_lo={info['Nc_lo']}, Nc_hi={info['Nc_hi']}, NO STRADDLE, Nc_used={info['Nc_used']}")

    if diag_0['premise_failures']:
        print(f"\n  PREMISE FAILURES (ρ ≤ β, fallback used):")
        for fail in diag_0['premise_failures'][:5]:
            print(f"    T={fail['T_mev']:.2f} meV: ρ={fail['rho']:.4f}, β={fail['beta']:.4f}, method={fail['method']}")
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

    upper_methods = [c['status']['method'] for c in diag_13['upper_checks']]
    temple_count = sum(1 for m in upper_methods if m == 'temple')
    print(f"  Upper bound: {temple_count}/{len(upper_methods)} used Temple, "
          f"{len(upper_methods)-temple_count} used fallback")

    if diag_13['width_samples']:
        print(f"\n  WIDTH SAMPLES (first few temperatures):")
        for sample in diag_13['width_samples'][:3]:
            print(f"    T={sample['T_mev']:.2f} meV ({sample['stage']}):")
            for key, val in sample['widths'].items():
                print(f"      {key}: {val:.2e}")

    if diag_13['nc_info']:
        print(f"\n  N_c INFO (first few):")
        for info in diag_13['nc_info'][:3]:
            if info['straddle']:
                print(f"    T={info['T_mev']:.2f} meV ({info['leg']}): Nc_lo={info['Nc_lo']}, Nc_hi={info['Nc_hi']}, STRADDLE, range={info['Nc_range']}")
            else:
                print(f"    T={info['T_mev']:.2f} meV ({info['leg']}): Nc_lo={info['Nc_lo']}, Nc_hi={info['Nc_hi']}, NO STRADDLE, Nc_used={info['Nc_used']}")

    if diag_13['premise_failures']:
        print(f"\n  PREMISE FAILURES (ρ ≤ β, fallback used):")
        for fail in diag_13['premise_failures'][:5]:
            print(f"    T={fail['T_mev']:.2f} meV: ρ={fail['rho']:.4f}, β={fail['beta']:.4f}, method={fail['method']}")
    print()

    print("="*80)
    print("SUMMARY")
    print("="*80)
    print(f"μ*=0:    Tc ∈ [{Tc_lo_0*MEV_K:.3f}, {Tc_hi_0*MEV_K:.3f}] K  (width: {(Tc_hi_0-Tc_lo_0)*MEV_K:.3f} K)")
    print(f"μ*=0.13: Tc ∈ [{Tc_lo_13*MEV_K:.3f}, {Tc_hi_13*MEV_K:.3f}] K  (width: {(Tc_hi_13-Tc_lo_13)*MEV_K:.3f} K)")
    print(f"Deficit: {100*(1-Tc_lo_13/Tc_lo_0):.2f}% to {100*(1-Tc_hi_13/Tc_hi_0):.2f}%")
    print()
