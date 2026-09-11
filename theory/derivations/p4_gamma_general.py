r"""
P4-C generalization: geometric-vector Kato-Temple bracket on g(gamma), swept over gamma
=========================================================================================

Extends the VERIFIED gamma=2 result (project P4-C, `p4_gap_conjecture2.py`, main worktree
`p3-numerics`) to a grid of gamma. Target operator (Kiessling-Altshuler-Yuzbashyan I,
arXiv:2409.00533v2), the gamma-model operator on l^2(N_0) (index n from 0):

    G(gamma) = -G1 + G2 + G3        (compact, self-adjoint for gamma>0)
    (G1)_{n,m} = [ (1/(2n+1)) sum_{k=1}^{n} 2/k^gamma ] delta_{n,m}
    (G2)_{n,m} = (1-delta_{n,m}) / ( sqrt(2n+1) |n-m|^gamma sqrt(2m+1) )
    (G3)_{n,m} = 1 / ( sqrt(2n+1) (n+m+1)^gamma sqrt(2m+1) )

g(gamma) := rho(G(gamma)) (largest eigenvalue).  T_c(gamma) = (1/2pi) g(gamma)^{1/gamma}.

KAY's PROVEN analytic upper-bound chain (survey P4_KAY_I_and_combescot):
    rho(G) <= rho(G2+G3) <= rho(G2)+rho(G3) <= 2*Prop5(gamma)+1
    Prop5(gamma) = sqrt( (2^{1+eps}-1) zeta(1+2gamma-eps) zeta(1+eps) ), eps=min{gamma,0.65}.
At gamma=2 this chain factor is ~2.03x C_inf; the VERIFIED P4-C result showed this is an
ARTIFACT of the chain's subadditivity+Cauchy-Schwarz steps, not of the operator: a
one-parameter geometric trial vector v_n=t^n (t*~0.394) + Kato-Temple gives a rigorous
two-sided bracket g(2) in [1.31714, 1.31820] (Tc factor 1.00002x, not 2.03x).

THIS SCRIPT asks: does the geometric-vector trick generalize to gamma != 2, and where does
it degrade?  Method per gamma:
  (a) Rayleigh LOWER bound mu(t) = <v|G|v>/<v|v>, v_n=t^n, t on a fixed grid -- EXACT finite
      sum (v has genuine finite support Nv < NG, so no truncation error in mu itself).
  (b) Kato-Temple UPPER bound g <= mu + eta^2/(mu-alpha), eta^2=<v|G^2|v>/<v|v>-mu^2 (with a
      rigorous tail correction for the NG-truncation of the matvec), alpha a rigorous upper
      bound on the SECOND eigenvalue lambda_2(G), built via Cholesky-deflation on a finite
      head block A (size Na) + a GENERIC (gamma-parametrized) tail bound on the operator
      norm of the remainder Y=[[0,B],[B^T,D]].

The gamma=2 tail-bound machinery (Gershgorin row sums, hardcoded exponent 2) does NOT carry
over verbatim -- exponents are gamma-specific and (for the D-block, tail-tail) row sums only
converge for gamma>1/2 in the naive sense.  We instead ADAPT KAY's own Cauchy-Schwarz proof
of Propositions 5 & 6 (which bounds rho(G2), rho(G3) via a Cauchy-Schwarz split that produces
a factor sum_{n>=0}(2n+1)^{-(1+eps)} = (1-2^{-(1+eps)})zeta(1+eps)) to the TAIL sub-block
(n,m>=N instead of n,m>=0).  The identical C-S argument (it never uses n=0 specifically) gives
    rho(G2|_{n,m>=N}) <= sqrt( zeta(1+2gamma-eps) * HurwitzZeta(1+eps, N+1/2) )
using sum_{n>=N}(2n+1)^{-(1+eps)} = 2^{-(1+eps)} HurwitzZeta(1+eps,N+1/2) exactly (mpmath).
SELF-CONSISTENCY CHECK (Stage 0): setting N=0 must reproduce KAY's own Prop5(gamma) exactly,
since HurwitzZeta(s,1/2)=(2^s-1)zeta(s) is a standard identity -- this is checked numerically
below and is the main correctness gate on the new (non-copy-pasted) tail-bound derivation.
The head-tail coupling block B's operator-norm bound (beta=rho(BB^T) via Collatz-Wielandt) DOES
generalize cleanly (same style as the gamma=2 script, just gamma-parametrized exponents) and is
valid for ALL gamma>0 (finite head, so no convergence threshold issue there).

Numerics discipline (per project convention): numpy float64 + eigvalsh for all heavy linear
algebra; mpmath dps=30 ONLY for final one-matvec-per-gamma certificates and Cholesky calls.
N truncation size capped at 400. t-optimization on a fixed grid, never open search. Fixed,
small iteration caps for the alpha-certificate's head-size escalation (Na in {150,250,350}).

HONESTY MANDATE: the alpha<mu Kato-Temple hypothesis and the tightness of the bracket are
NOT forced to succeed. Where they fail or loosen for small gamma, this is reported as the
actual finding (search for gamma*, the threshold below which the SINGLE geometric vector
degrades), not papered over.

Run:  python theory/derivations/p4_gamma_general.py
"""
import time
import numpy as np
import mpmath as mp

mp.mp.dps = 30

CHECKS = []


def _check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "FAIL"
    print(f"  [{tag}] {name}" + (f"  {detail}" if detail else ""))
    assert ok, f"{name} FAILED: {detail}"


def _soft_check(name, ok, detail=""):
    """Like _check but does NOT assert -- used for gamma-dependent behavior we expect
    may legitimately fail (e.g. alpha<mu at small gamma). Still logged into CHECKS."""
    CHECKS.append((name, bool(ok)))
    tag = "PASS" if ok else "INFO (expected possible failure)"
    print(f"  [{tag}] {name}" + (f"  {detail}" if detail else ""))


def _stage_summary(k, note=""):
    n_pass = sum(1 for _, ok in CHECKS if ok)
    n_total = len(CHECKS)
    print(f"\nStage {k}: {n_pass}/{n_total} checks passed (cumulative){note}\n")


# ============================================================================
# Operator construction (inlined & self-contained; structure copied from the
# VERIFIED gamma=2 builder p4_gap_conjecture2.py -- already gamma-generic there,
# NO import, no change needed beyond re-declaring here for self-containedness).
# ============================================================================

def build_blocks_np(N, gamma=2.0):
    """float64 blocks (-G1 diag, G2, G3) and G = -G1+G2+G3 as N x N arrays."""
    gamma = float(gamma)
    n = np.arange(N)
    sq = np.sqrt(2.0 * n + 1.0)
    g1d = np.zeros(N)
    running = 0.0
    for k in range(N):
        if k >= 1:
            running += 2.0 / (k ** gamma)
        g1d[k] = running / (2 * k + 1)
    idx = np.arange(N)
    Dm = np.abs(idx[:, None] - idx[None, :])
    G2 = np.where(Dm > 0,
                  1.0 / (sq[:, None] * np.where(Dm > 0, Dm, 1) ** gamma * sq[None, :]),
                  0.0)
    Sm = (idx[:, None] + idx[None, :] + 1.0)
    G3 = 1.0 / (sq[:, None] * Sm ** gamma * sq[None, :])
    G1 = np.diag(g1d)
    return G1, G2, G3, (-G1 + G2 + G3)


def build_G_mp(N, gamma=2):
    """Exact mpmath G = -G1 + G2 + G3 (dps=30). O(N^2), seconds for N<=400."""
    gamma = mp.mpf(gamma)
    sq = [mp.sqrt(2 * n + 1) for n in range(N)]
    G1neg = []
    run = mp.mpf(0)
    for n in range(N):
        if n >= 1:
            run += 2 / mp.mpf(n) ** gamma
        G1neg.append(-run / (2 * n + 1))
    G = mp.matrix(N)
    for n in range(N):
        sqn = sq[n]
        for m in range(N):
            if n == m:
                G[n, m] = G1neg[n] + 1 / (sqn * mp.mpf(n + m + 1) ** gamma * sq[m])
            else:
                G[n, m] = (1 / (sqn * mp.mpf(abs(n - m)) ** gamma * sq[m])
                           + 1 / (sqn * mp.mpf(n + m + 1) ** gamma * sq[m]))
    return G


def prop5_rhs(gamma):
    """Prop 5 RHS = sqrt((2^{1+eps}-1) zeta(1+2gamma-eps) zeta(1+eps)), eps=min{gamma,0.65}."""
    gamma = mp.mpf(gamma)
    eps = min(gamma, mp.mpf("0.65"))
    return mp.sqrt((mp.mpf(2) ** (1 + eps) - 1) * mp.zeta(1 + 2 * gamma - eps) * mp.zeta(1 + eps)), eps


# ============================================================================
# NEW generic tail-bound machinery for the Kato-Temple alpha certificate.
# ============================================================================

def prop5_tail(N, gamma):
    """Rigorous upper bound on rho(G2 restricted to indices n,m>=N), ADAPTED from KAY's
    own Cauchy-Schwarz proof of Prop 5 (identical steps, just n ranges over n>=N instead
    of n>=0): rho(G2|_{>=N}) <= sqrt( zeta(1+2gamma-eps) * HurwitzZeta(1+eps, N+1/2) ).
    At N=0, HurwitzZeta(1+eps,1/2)=(2^{1+eps}-1)zeta(1+eps) recovers prop5_rhs EXACTLY
    (checked numerically in Stage 0) -- this is the correctness gate on this adaptation."""
    gamma = mp.mpf(gamma)
    eps = min(gamma, mp.mpf("0.65"))
    return mp.sqrt(mp.zeta(1 + 2 * gamma - eps) * mp.zeta(1 + eps, mp.mpf(N) + mp.mpf("0.5")))


def _harmonic_power_sum_bound(n, gamma):
    """Rigorous bound: sum_{k=1}^{n} k^{-gamma} <= 1 + integral_1^n x^{-gamma} dx
    (integral test for the decreasing tail k=2..n), closed form, valid all gamma>0, n>=1."""
    n = mp.mpf(n)
    gamma = mp.mpf(gamma)
    if abs(gamma - 1) < mp.mpf("1e-12"):
        return 1 + mp.log(n)
    return 1 + (n ** (1 - gamma) - 1) / (1 - gamma)


def g1_diag_tail_bound(N, gamma, grid_points=300, span_mult=200):
    """Rigorous bound on max_{n>=N} |G1_nn|, gamma-generic.
    G1_nn = (2/(2n+1)) sum_{k=1}^n k^-gamma <= f(n,gamma) := 2*H(n,gamma)/(2n+1), a closed
    form (via _harmonic_power_sum_bound) that ->0 as n->inf. f is NOT hand-proven monotone
    on n>=N for all gamma (no closed monotonicity proof attempted here -- would require a
    derivative sign argument depending on gamma), so we instead take the max of f over a
    fine geometrically-spaced grid n in [N, span_mult*N], which upper-bounds the true sup
    given f's provable eventual decay (f ~ n^{-gamma}/(2(1-gamma)) or ~ln(n)/(2n) for large
    n -- both ->0). Reported HONESTLY as a numerically-maximized closed-form pointwise
    bound, not a fully hand-proven supremum."""
    ns = np.unique(np.round(np.geomspace(N, N * span_mult, grid_points)).astype(np.int64))
    best = mp.mpf(0)
    for n in ns:
        n = int(n)
        if n < 1:
            continue
        H = _harmonic_power_sum_bound(n, gamma)
        f = 2 * H / (2 * mp.mpf(n) + 1)
        if f > best:
            best = f
    return best


def g3_diag_tail_bound(N, gamma):
    """max_{n>=N} G3_nn = 1/(2N+1)^{1+gamma} EXACTLY (G3_nn=1/(2n+1)^{1+gamma}, manifestly
    monotone decreasing in n -- no ambiguity, unlike the G1 case)."""
    return 1 / (2 * mp.mpf(N) + 1) ** (1 + mp.mpf(gamma))


def d_tail_bound(N, gamma):
    """Rigorous bound: rho(D_tail) <= ||G1_tail||_op + ||G2_tail||_op + ||G3_tail||_op
    (triangle inequality, D_tail = -G1|_{>=N} + G2|_{>=N} + G3|_{>=N}), using
    g1_diag_tail_bound (diagonal op norm), prop5_tail (G2 off-diag, C-S adaptation),
    and g3_diag_tail_bound + prop5_tail (G3: diag + off-diag, SAME C-S structure as G2's
    since KAY's Prop 6 proof uses the identical Cauchy-Schwarz off-diagonal estimate)."""
    return g1_diag_tail_bound(N, gamma) + g3_diag_tail_bound(N, gamma) + 2 * prop5_tail(N, gamma)


def beta_bbT_bound(N, gamma, w, M_cut):
    """Rigorous beta >= rho(B B^T), B = head(N) x tail coupling, GENERIC gamma (generalizes
    the gamma=2 script's function: only the exponents and the K(gamma) tail constant change).

    Collatz-Wielandt on nonnegative S=BB^T with positive w: rho(S) <= max_n (Sw)_n/w_n.
    Partial float64 sums m in [N, M_cut) (all positive; 1e-9 inflated) plus a closed-form
    remainder for m >= M_cut >= 4N, using: for n<N<=M_cut/4<=m, m-n>3m/4 and m+n+1>m, so
        B_nm < K(gamma)/(sqrt(2n+1) sqrt(2m+1) m^gamma) < K(gamma)/(sqrt(2n+1) m^{gamma+1/2}),
        K(gamma) := (4/3)^gamma + 1,
    giving  sum_{m>=M_cut} m^{-(2gamma+1)} <= (M_cut-1)^{-2gamma}/(2gamma)  (integral bound,
    valid for ALL gamma>0 since 2gamma+1>1 always)."""
    assert M_cut >= 4 * N, "need M_cut >= 4N for the tail bound"
    gamma_f = float(gamma)
    inv_sqrt = 1.0 / np.sqrt(2.0 * np.arange(float(N)) + 1.0)
    w_np = np.array([float(x) for x in w])
    n_arr = np.arange(float(N))
    Sw = np.zeros(N)
    m = N
    while m < M_cut:
        sqm = np.sqrt(2.0 * m + 1.0)
        mdm = m - n_arr
        Bcol = (1.0 / mdm ** gamma_f + 1.0 / (m + n_arr + 1.0) ** gamma_f) * inv_sqrt / sqm
        Sw += Bcol * float(np.dot(Bcol, w_np))
        m += 1
    M = mp.mpf(M_cut)
    gamma_mp = mp.mpf(gamma)
    K = (mp.mpf(4) / 3) ** gamma_mp + 1
    W1 = sum(w)
    tail_max = K ** 2 * W1 * (M - 1) ** (-2 * gamma_mp) / (2 * gamma_mp)
    beta = mp.mpf(0)
    for n in range(N):
        Swn = mp.mpf(float(Sw[n])) * (1 + mp.mpf("1e-9"))
        cand = (Swn + tail_max / mp.sqrt(2 * n + 1)) / w[n]
        if cand > beta:
            beta = cand
    return beta * (1 + mp.mpf("1e-9"))


def alpha_certificate(Na, gamma):
    """Rigorous upper bound on lambda_2(G(gamma)) via Cholesky-deflation on head block A
    (size Na) + generic tail bound on the operator norm of the off-head remainder Y.
    Returns (alpha, chol_ok, lam2A_up, beta, d, elapsed_s). gamma-generic; structure
    identical to the VERIFIED gamma=2 script's alpha-certificate, just build_G_mp/blocks
    and the tail bounds (beta_bbT_bound, d_tail_bound) are now gamma-parametrized."""
    t0 = time.time()
    _, _, _, A_np = build_blocks_np(Na, gamma)
    wA, vA = np.linalg.eigh(A_np)
    u = vA[:, -1]
    # Anchor c to the ACTUAL float64 top eigenvalue of the deflated matrix PAP, not the raw
    # wA[-2] estimate of lambda_2(A): u is only a float64-precision approximate eigenvector,
    # so PAP = A - u(Au)^T - (Au)u^T + (u.Au)uu^T retains a small numerical residual (~1e-15,
    # NOT exactly 0) in the deflated direction. When lambda_2(A) itself is tiny (or negative,
    # as happens at large gamma where the spectrum below the top eigenvalue is compressed
    # near 0), lambda2A_num+1e-7 can be SMALLER than this residual, making cI-PAP fail to be
    # SPD even though the deflation is doing its job -- anchoring on PAP's own top eigenvalue
    # fixes this uniformly for all gamma.
    Au_f = A_np @ u
    uAu_f = float(u @ Au_f)
    PAP_f = A_np - np.outer(u, Au_f) - np.outer(Au_f, u) + uAu_f * np.outer(u, u)
    lam_top_PAP_num = float(np.linalg.eigvalsh(PAP_f)[-1])
    c = mp.mpf(lam_top_PAP_num) + mp.mpf("1e-7")
    G_mp_a = build_G_mp(Na, gamma)
    umb = mp.matrix([mp.mpf(float(x)) for x in u])
    Au = G_mp_a * umb
    uAu = sum(umb[i] * Au[i] for i in range(Na))
    M = mp.matrix(Na)
    for i in range(Na):
        for jj in range(Na):
            pap = G_mp_a[i, jj] - umb[i] * Au[jj] - Au[i] * umb[jj] + umb[i] * umb[jj] * uAu
            M[i, jj] = (c if i == jj else mp.mpf(0)) - pap
    chol_ok = True
    try:
        mp.cholesky(M)
    except Exception:  # noqa
        chol_ok = False
    if not chol_ok:
        return None, False, None, None, None, time.time() - t0
    lam2A_up = c
    d = d_tail_bound(Na, gamma)
    beta = beta_bbT_bound(Na, gamma, [mp.mpf(1)] * Na, 4 * Na)
    alpha = lam2A_up + mp.sqrt(beta) + d
    return alpha, True, lam2A_up, beta, d, time.time() - t0


# ============================================================================
# Semi-analytic bonus: closed-form pieces of the Rayleigh sums for v_n=t^n.
# ============================================================================

def vv_closed(t, Nv):
    """<v|v> = sum_{n=0}^{Nv-1} t^{2n} = (1-t^{2Nv})/(1-t^2). Exact finite geometric series."""
    t = mp.mpf(t)
    return (1 - t ** (2 * Nv)) / (1 - t ** 2)


def g1_rayleigh_lerch(t, Nv, gamma):
    """<v|G1|v> = sum_{n=0}^{Nv-1} t^{2n} G1_nn, v_n=t^n (n<Nv), CLOSED FORM via Lerch Phi:
        <v|G1|v> = 2 sum_{k=1}^{Nv-1} k^-gamma *
                   [ (t^{2k}/2) LerchPhi(t^2,1,k+1/2) - (t^{2Nv}/2) LerchPhi(t^2,1,Nv+1/2) ]
    Derivation: swap summation order (k<=n<Nv) then use the standard finite-tail identity
    sum_{n=k}^{Nv-1} t^{2n}/(2n+1) = (t^{2k}/2)LerchPhi(t^2,1,k+1/2) - (t^{2Nv}/2)LerchPhi(t^2,1,Nv+1/2)
    (difference of two infinite LerchPhi tails, exact since |t|<1). Verified numerically
    against the direct O(Nv) running sum in Stage 0 (exact match to dps=30).
    NOTE: this reduces the (n,k) double sum to a single k-sum of LerchPhi calls -- a genuine
    closed form per term, though NOT asymptotically cheaper than the direct O(Nv) running-sum
    build (both are O(Nv)). G2 and G3's off-diagonal double sums over (n,m) do NOT reduce this
    way: the cross term 1/(sqrt(2n+1)sqrt(2m+1)) couples n and m non-separably with the
    |n-m|^-gamma / (n+m+1)^-gamma kernel, so no closed polylog/Lerch form was found for them;
    they are kept as exact truncated finite sums (rigorous since v has genuine finite support)."""
    t = mp.mpf(t)
    gamma = mp.mpf(gamma)
    total = mp.mpf(0)
    tail_Nv = (t ** (2 * Nv) / 2) * mp.lerchphi(t ** 2, 1, mp.mpf(Nv) + mp.mpf("0.5"))
    for k in range(1, Nv):
        head_k = (t ** (2 * k) / 2) * mp.lerchphi(t ** 2, 1, k + mp.mpf("0.5"))
        total += (1 / mp.mpf(k) ** gamma) * (head_k - tail_Nv)
    return 2 * total


# ============================================================================
# Reference constants (gamma=2, VERIFIED upstream P4-C bracket).
# ============================================================================
G_TRUE_REF_2 = mp.mpf("1.3181402121")   # corrected 2026-09-11 (P4-B: (2piC_inf)^2 = 1.31814021; conj2 Segment I enclosure width 3.7e-10); reference constant for evidence asserts only
C_INF_REF_2 = mp.mpf("0.1827262477")
KT_LOWER_REF_2 = mp.mpf("1.31714")
KT_UPPER_REF_2 = mp.mpf("1.31820")

GAMMA_GRID = [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]
NG = 400          # matvec-certificate truncation size (cap 400 per project discipline)
NV = 90           # trial-vector support (4*NV=360<=NG=400, safety margin)
NA_CANDIDATES = [150, 250, 350]   # fixed escalation ladder for the alpha certificate
T_GRID = np.linspace(0.03, 0.75, 145)   # fixed grid, no open search


def rayleigh_mu_eta2_float(G_np, t, Nv):
    """float64 mu(t), eta2(t) (fast, for grid search only -- NOT the final certificate)."""
    N = G_np.shape[0]
    v = np.zeros(N)
    v[:Nv] = t ** np.arange(Nv)
    Gv = G_np @ v
    vv = v @ v
    mu = (v @ Gv) / vv
    eta2 = (Gv @ Gv) / vv - mu * mu
    return mu, eta2


def gevp_2x2_top(H, S):
    """Largest generalized eigenvalue of 2x2 symmetric (H,S), S SPD: solve det(H-lam S)=0."""
    a = S[0, 0] * S[1, 1] - S[0, 1] ** 2
    b = -(H[0, 0] * S[1, 1] + H[1, 1] * S[0, 0] - 2 * H[0, 1] * S[0, 1])
    c = H[0, 0] * H[1, 1] - H[0, 1] ** 2
    disc = b * b - 4 * a * c
    if disc < 0 or a == 0:
        return -np.inf
    lam1 = (-b + np.sqrt(disc)) / (2 * a)
    lam2 = (-b - np.sqrt(disc)) / (2 * a)
    return max(lam1, lam2)


def two_param_bonus_mu2(G_np, t, Nv):
    """Rigorous improved LOWER bound via v = a*t^n + b*n*t^n (2-parameter family), the
    generalized Rayleigh quotient over span{t^n, n t^n} (n<Nv), a genuine 2D subspace of
    l^2 -- max Rayleigh quotient over ANY finite-dim subspace is <= true top eigenvalue, so
    mu2>=mu is guaranteed valid whenever numerically well-conditioned (S SPD)."""
    N = G_np.shape[0]
    n_arr = np.arange(Nv)
    v1 = np.zeros(N)
    v1[:Nv] = t ** n_arr
    v2 = np.zeros(N)
    v2[:Nv] = n_arr * (t ** n_arr)
    Gv1 = G_np @ v1
    Gv2 = G_np @ v2
    H = np.array([[v1 @ Gv1, v1 @ Gv2], [v1 @ Gv2, v2 @ Gv2]])
    S = np.array([[v1 @ v1, v1 @ v2], [v1 @ v2, v2 @ v2]])
    return gevp_2x2_top(H, S)


def run_gamma(gamma, tag=""):
    """Full pipeline for one gamma: returns a dict of results (table row)."""
    print(f"\n  --- gamma = {gamma} {tag} ---")
    _, _, _, G_NG_np = build_blocks_np(NG, gamma)
    g_ref = float(np.linalg.eigvalsh(G_NG_np)[-1])  # REFERENCE ONLY, not part of the bound

    # (a) float64 mu(t) grid search
    mus, eta2s = [], []
    for t in T_GRID:
        mu_f, eta2_f = rayleigh_mu_eta2_float(G_NG_np, t, NV)
        mus.append(mu_f)
        eta2s.append(eta2_f)
    mus = np.array(mus)
    i_mu = int(np.argmax(mus))
    t_mu_best, mu_float_best = T_GRID[i_mu], mus[i_mu]
    boundary_hit = i_mu == 0 or i_mu == len(T_GRID) - 1

    # (b) 2-parameter bonus lower bound (grid search, float64)
    mu2s = np.array([two_param_bonus_mu2(G_NG_np, t, NV) for t in T_GRID])
    i_mu2 = int(np.argmax(mu2s))
    t_mu2_best, mu2_float_best = T_GRID[i_mu2], mu2s[i_mu2]

    # (c) alpha certificate: escalate Na until alpha < mu_float_best (with margin) or exhausted
    alpha = None
    Na_used = None
    alpha_ok = False
    for Na in NA_CANDIDATES:
        a, chol_ok, lam2A_up, beta, d, elapsed = alpha_certificate(Na, gamma)
        print(f"    [alpha cert] Na={Na}: chol_ok={chol_ok}, alpha={float(a) if a is not None else None}"
              f"  ({elapsed:.1f}s)")
        if chol_ok and a < mu_float_best * (1 - 1e-6):
            alpha, Na_used, alpha_ok = a, Na, True
            break
        elif chol_ok:
            alpha, Na_used = a, Na  # keep the last (largest) attempt for honest reporting

    # (d) KT-bound t-optimization (only meaningful if alpha_ok)
    if alpha_ok:
        kts = []
        for i, t in enumerate(T_GRID):
            mu_f, eta2_f = mus[i], eta2s[i]
            if mu_f > float(alpha):
                kts.append(mu_f + eta2_f / (mu_f - float(alpha)))
            else:
                kts.append(np.inf)
        i_kt = int(np.argmin(kts))
        t_star = T_GRID[i_kt]
    else:
        t_star = t_mu_best

    # (e) rigorous mpmath certificate at t_star
    t0 = time.time()
    G_mp = build_G_mp(NG, gamma)
    tt = mp.mpf(float(t_star))
    v = mp.matrix(NG, 1)
    for n in range(NV):
        v[n, 0] = tt ** n
    vv = sum(v[i, 0] ** 2 for i in range(NG))
    Gv = G_mp * v
    xGx = sum(v[i, 0] * Gv[i, 0] for i in range(NG))
    mu = xGx / vv  # rigorous LOWER bound (exact finite sum, v has genuine finite support)
    normGv2 = sum(Gv[i, 0] ** 2 for i in range(NG))
    gamma_mp = mp.mpf(gamma)
    K = (mp.mpf(4) / 3) ** gamma_mp + 1
    S_geo = 1 / (1 - tt)  # rigorous bound on sum_n t^n/sqrt(2n+1) <= sum_n t^n
    tail_eta = K ** 2 * S_geo ** 2 * (mp.mpf(NG) - 1) ** (-2 * gamma_mp) / (2 * gamma_mp)
    eta2_up = (normGv2 + tail_eta) / vv - mu ** 2
    mp_elapsed = time.time() - t0

    kt_upper = None
    if alpha_ok and mu > alpha:
        kt_upper = mu + eta2_up / (mu - alpha)

    # Tc(gamma) = (1/2pi) g^{1/gamma}; C_inf(gamma) proxy = (1/2pi) g_ref^{1/gamma}
    tc_up_factor = None
    if kt_upper is not None:
        tc_up_factor = (kt_upper / mp.mpf(g_ref)) ** (mp.mpf(1) / gamma_mp)
    prop5, eps = prop5_rhs(gamma)
    chain_val = 2 * prop5 + 1
    kay_chain_tc_factor = (chain_val / mp.mpf(g_ref)) ** (mp.mpf(1) / gamma_mp)

    deficit = g_ref - float(mu)
    mu2_deficit = g_ref - mu2_float_best
    bracket_width = float(kt_upper - mu) if kt_upper is not None else None

    row = dict(
        gamma=gamma, t_star=float(t_star), mu=float(mu), g_ref=g_ref,
        kt_upper=(float(kt_upper) if kt_upper is not None else None),
        bracket_width=bracket_width, alpha_ok=alpha_ok, Na_used=Na_used,
        alpha=(float(alpha) if alpha is not None else None),
        deficit=deficit, mu2=mu2_float_best, mu2_deficit=mu2_deficit, t_mu2=float(t_mu2_best),
        tc_up_factor=(float(tc_up_factor) if tc_up_factor is not None else None),
        kay_chain_tc_factor=float(kay_chain_tc_factor),
        boundary_hit=boundary_hit, mp_elapsed=mp_elapsed,
    )
    print(f"    t*={t_star:.3f}  mu={float(mu):.6f}  g_ref={g_ref:.6f}  deficit={deficit:.3e}"
          f"  alpha_ok={alpha_ok}(Na={Na_used})"
          + (f"  kt_upper={float(kt_upper):.6f}  width={bracket_width:.3e}"
             f"  Tc_factor={float(tc_up_factor):.5f}" if kt_upper is not None else "  kt_upper=N/A"))
    print(f"    2-param bonus: mu2={mu2_float_best:.6f} (deficit {mu2_deficit:.3e} vs single-vector"
          f" {deficit:.3e}, improvement {'%.1fx' % (deficit/max(mu2_deficit,1e-300)) if mu2_deficit>0 else 'n/a'})")
    print(f"    KAY chain Tc factor = {float(kay_chain_tc_factor):.5f}")
    return row


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    T_START = time.time()

    # ------------------------------------------------------------------
    # Stage 0: operator + new tail-bound machinery self-tests
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Stage 0: operator construction + tail-bound-derivation self-tests")
    print("=" * 70)

    # g^(1)(gamma)=1 exactly for ANY gamma (G1^(1)=0 empty sum, G3_{00}=1/(0+0+1)^gamma=1)
    for gm in [0.5, 1.0, 2.0, 5.0]:
        _, _, _, G1x1 = build_blocks_np(1, gm)
        _check(f"g^(1)(gamma={gm}) = 1 exactly (float64)", abs(G1x1[0, 0] - 1.0) < 1e-14,
               f"got {G1x1[0,0]:.15f}")

    # 4x4 explicit entries at gamma=2 (known closed forms, survey Sec.1.4) -- validates builder
    _, _, _, G4 = build_blocks_np(4, 2.0)
    _check("G[0,0]=1 (4x4, gamma=2)", abs(G4[0, 0] - 1.0) < 1e-14, f"got {G4[0,0]}")
    _check("G[0,1]=(1/4+1)/sqrt3 (gamma=2)", abs(G4[0, 1] - (1.25) / np.sqrt(3)) < 1e-14,
           f"got {G4[0,1]:.12f}")
    for gm in [0.5, 1.0, 1.5, 3.0]:
        _, _, _, G4g = build_blocks_np(4, gm)
        _check(f"G symmetric (4x4, gamma={gm})", np.allclose(G4g, G4g.T), "")

    # Theorem 2: g^(N) monotone increasing with N, for gamma>0 (holds for all gamma per KAY)
    for gm in [0.5, 1.0, 2.0, 5.0]:
        g_prev = -1.0
        for Nt in [50, 100, 200, 400]:
            _, _, _, GN = build_blocks_np(Nt, gm)
            g_here = float(np.linalg.eigvalsh(GN)[-1])
            # tolerance 1e-9 (not 1e-15): at large gamma, g^(N) saturates to float64 noise
            # floor well before N=400, and LAPACK eigvalsh noise at that scale exceeds 1e-15.
            _check(f"monotone-up g^({Nt})(gamma={gm}) > g^(prev)", g_here > g_prev - 1e-9,
                   f"g^({Nt})={g_here:.6f}")
            g_prev = g_here

    # gamma=2 400-truncation matches the known reference (5e-7)
    _, _, _, G400_2 = build_blocks_np(400, 2.0)
    g400_2 = float(np.linalg.eigvalsh(G400_2)[-1])
    _check("g^(400)(2) matches known reference 1.3181402 (5e-7)",
           abs(g400_2 - float(G_TRUE_REF_2)) < 5e-7, f"got {g400_2:.9f}")

    # CRITICAL correctness gate: prop5_tail(N=0,gamma) == prop5_rhs(gamma) EXACTLY
    # (Hurwitz zeta identity zeta(s,1/2)=(2^s-1)zeta(s) -- validates the tail-bound adaptation)
    for gm in [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0]:
        p5_orig, _ = prop5_rhs(gm)
        p5_tail0 = prop5_tail(0, gm)
        diff = abs(p5_orig - p5_tail0)
        _check(f"prop5_tail(N=0,gamma={gm}) == prop5_rhs(gamma) (Hurwitz-zeta identity)",
               diff < mp.mpf("1e-25"), f"diff={mp.nstr(diff,4)}")

    # d_tail_bound / prop5_tail must DECAY as N grows (sanity: tail shrinks with more head)
    for gm in [0.75, 2.0]:
        d_small = d_tail_bound(50, gm)
        d_large = d_tail_bound(200, gm)
        _check(f"d_tail_bound decreasing in N (gamma={gm})", d_large < d_small,
               f"d(50)={mp.nstr(d_small,4)}, d(200)={mp.nstr(d_large,4)}")

    # semi-analytic bonus: Lerch-Phi closed form for <v|G1|v> matches direct finite sum
    for (tt, Nvv, gm) in [(0.4, 20, 1.3), (0.3, 30, 0.75), (0.35, 25, 3.0)]:
        direct = mp.mpf(0)
        run = mp.mpf(0)
        for nn in range(Nvv):
            if nn >= 1:
                run += 2 / mp.mpf(nn) ** mp.mpf(gm)
            direct += mp.mpf(tt) ** (2 * nn) * run / (2 * nn + 1)
        lerch = g1_rayleigh_lerch(tt, Nvv, gm)
        _check(f"Lerch-Phi closed form <v|G1|v> matches direct sum (t={tt},Nv={Nvv},gamma={gm})",
               abs(direct - lerch) < mp.mpf("1e-25"), f"diff={mp.nstr(direct-lerch,4)}")

    # <v|v> closed geometric form matches direct sum
    for (tt, Nvv) in [(0.4, 50), (0.6, 30)]:
        direct = sum(mp.mpf(tt) ** (2 * n) for n in range(Nvv))
        closed = vv_closed(tt, Nvv)
        _check(f"<v|v> closed geometric form matches direct sum (t={tt},Nv={Nvv})",
               abs(direct - closed) < mp.mpf("1e-25"), f"diff={mp.nstr(direct-closed,4)}")

    _stage_summary(0)

    # ------------------------------------------------------------------
    # Stage 1: gamma=2 SELF-CHECK GATE -- must reproduce known [1.31714,1.31820]
    #          using the NEW generic pipeline (run_gamma), not the old script.
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Stage 1: gamma=2 reproduction gate (new generic pipeline vs known bracket)")
    print("=" * 70)

    row2 = run_gamma(2.0, tag="[SELF-CHECK GATE]")

    _check("gamma=2: alpha<mu (Kato-Temple hypothesis) holds", row2["alpha_ok"],
           f"alpha={row2['alpha']}")
    _check("gamma=2: mu is a valid LOWER bound on known g(2) (mu <= 1.3181402+1e-6)",
           row2["mu"] <= float(G_TRUE_REF_2) + 1e-6, f"mu={row2['mu']:.8f}")
    _check("gamma=2: mu within 2e-3 of KNOWN mu=1.31714 (new pipeline reproduces old script)",
           abs(row2["mu"] - float(KT_LOWER_REF_2)) < 2e-3, f"mu={row2['mu']:.8f}")
    _check("gamma=2: kt_upper within 2e-3 of KNOWN kt_upper=1.31820",
           abs(row2["kt_upper"] - float(KT_UPPER_REF_2)) < 2e-3, f"kt_upper={row2['kt_upper']:.8f}")
    _check("gamma=2: kt_upper >= known g(2) reference (valid upper bound)",
           row2["kt_upper"] >= float(G_TRUE_REF_2) - 1e-6, f"kt_upper={row2['kt_upper']:.8f}")
    _check("gamma=2: Tc upper factor < 1.001 (matches VERIFIED P4-C ~1.00002x)",
           row2["tc_up_factor"] < 1.001, f"factor={row2['tc_up_factor']:.6f}")
    _check("gamma=2: Tc upper factor beats KAY chain factor by >50x margin",
           (row2["kay_chain_tc_factor"] - 1) / (row2["tc_up_factor"] - 1) > 50,
           f"ours={row2['tc_up_factor']:.5f}, chain={row2['kay_chain_tc_factor']:.5f}")
    print("\n  *** GATE PASSED: new generic (gamma-parametrized) pipeline reproduces the")
    print("      VERIFIED gamma=2 bracket [1.31714, 1.31820] from p4_gap_conjecture2.py. ***")
    _stage_summary(1)

    # ------------------------------------------------------------------
    # Stage 2: sweep gamma grid
    # ------------------------------------------------------------------
    print("=" * 70)
    print(f"Stage 2: sweep gamma in {GAMMA_GRID}")
    print("=" * 70)

    rows = {2.0: row2}
    for gm in GAMMA_GRID:
        if gm in rows:
            continue
        rows[gm] = run_gamma(gm)

    # basic structural checks on the sweep (soft where gamma-dependent behavior is expected)
    for gm in GAMMA_GRID:
        r = rows[gm]
        _check(f"gamma={gm}: mu <= g_ref (valid lower bound)", r["mu"] <= r["g_ref"] + 1e-9,
               f"mu={r['mu']:.6f}, g_ref={r['g_ref']:.6f}")
        _check(f"gamma={gm}: mu2 (2-param) >= mu (single-vector)", r["mu2"] >= r["mu"] - 1e-9,
               f"mu2={r['mu2']:.6f}, mu={r['mu']:.6f}")
        if r["kt_upper"] is not None:
            _check(f"gamma={gm}: kt_upper >= g_ref (valid upper bound vs reference)",
                   r["kt_upper"] >= r["g_ref"] - 1e-6, f"kt_upper={r['kt_upper']:.6f}")
        _soft_check(f"gamma={gm}: alpha<mu (Kato-Temple hypothesis) held",
                    r["alpha_ok"], f"Na_used={r['Na_used']}")

    _stage_summary(2)

    # ------------------------------------------------------------------
    # Stage 3: table + honest gamma* verdict
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Stage 3: table + effective-gamma-range verdict")
    print("=" * 70)

    header = (f"{'gamma':>6} {'t*':>6} {'mu':>10} {'g_ref':>10} {'width':>10} "
              f"{'a<mu':>5} {'TcFac':>9} {'KAYchain':>9} {'mu2def/mudef':>13}")
    print(header)
    print("-" * len(header))
    for gm in GAMMA_GRID:
        r = rows[gm]
        width_s = f"{r['bracket_width']:.2e}" if r['bracket_width'] is not None else "N/A"
        tcf_s = f"{r['tc_up_factor']:.5f}" if r['tc_up_factor'] is not None else "N/A"
        ratio = (r["deficit"] / r["mu2_deficit"]) if r["mu2_deficit"] > 0 else float("inf")
        print(f"{gm:6.2f} {r['t_star']:6.3f} {r['mu']:10.6f} {r['g_ref']:10.6f} {width_s:>10} "
              f"{'Y' if r['alpha_ok'] else 'N':>5} {tcf_s:>9} {r['kay_chain_tc_factor']:9.4f} "
              f"{ratio:13.2f}")

    # find gamma*: smallest gamma in grid (scanning from large to small) where factor<1.01 holds
    tight = [gm for gm in GAMMA_GRID if rows[gm]["tc_up_factor"] is not None
             and rows[gm]["tc_up_factor"] < 1.01]
    loose_or_failed = [gm for gm in GAMMA_GRID if gm not in tight]
    gamma_star = min(tight) if tight else None

    print()
    print("  VERDICT:")
    if gamma_star is not None:
        print(f"    Geometric single-vector Kato-Temple bound is TIGHT (Tc factor < 1.01x) for"
              f" gamma >= gamma* = {gamma_star} (grid-resolution).")
    else:
        print("    Geometric single-vector Kato-Temple bound did NOT reach <1.01x tightness"
              " anywhere on this grid.")
    if loose_or_failed:
        print(f"    Degrades / fails for gamma in {sorted(loose_or_failed)}:")
        for gm in sorted(loose_or_failed):
            r = rows[gm]
            reason = []
            if not r["alpha_ok"]:
                reason.append("alpha<mu (Kato-Temple hypothesis) FAILED even at Na=350")
            if r["tc_up_factor"] is not None and r["tc_up_factor"] >= 1.01:
                reason.append(f"Tc factor {r['tc_up_factor']:.4f} >= 1.01 (bracket too wide)")
            print(f"      gamma={gm}: " + "; ".join(reason) if reason else f"      gamma={gm}: (see table)")
    _check("at least gamma=2,3,5 (large-gamma end) are in the tight set",
           all(gm in tight for gm in [2.0, 3.0, 5.0] if gm in GAMMA_GRID),
           f"tight set={tight}")
    _stage_summary(3)

    # ------------------------------------------------------------------
    # Final synthesis
    # ------------------------------------------------------------------
    print("=" * 70)
    print("SYNTHESIS")
    print("=" * 70)
    print(f"  Stage 1 gate: gamma=2 reproduces known bracket [1.31714,1.31820] via the NEW")
    print(f"    generic pipeline -> mu={row2['mu']:.6f}, kt_upper={row2['kt_upper']:.6f},")
    print(f"    Tc factor={row2['tc_up_factor']:.6f} (vs KAY chain {row2['kay_chain_tc_factor']:.4f}).")
    print(f"  Stage 2/3 sweep: tight-bracket set (Tc factor<1.01x) = {sorted(tight)}")
    print(f"                   degraded/failed set = {sorted(loose_or_failed)}")
    print("  Semi-analytic bonus: <v|v> has an exact closed geometric form; <v|G1|v> has an")
    print("    exact closed form via Lerch Phi (verified). <v|G2|v>,<v|G3|v> off-diagonal")
    print("    double sums resist closed form (sqrt(2n+1)sqrt(2m+1) couples n,m non-")
    print("    separably with the power-law kernel) -- kept as exact truncated finite sums.")
    print("  New tail-bound derivation (prop5_tail via Hurwitz zeta, adapting KAY's own")
    print("    Cauchy-Schwarz proof of Prop 5/6 to the n,m>=N sub-block) is validated by an")
    print("    exact N=0 reduction to KAY's published Prop5(gamma) (Stage 0 gate).")
    n_pass = sum(1 for _, ok in CHECKS if ok)
    n_total = len(CHECKS)
    print("=" * 70)
    print(f"FINAL: {n_pass}/{n_total} checks passed   (total {time.time()-T_START:.1f}s)")
    print("=" * 70)
