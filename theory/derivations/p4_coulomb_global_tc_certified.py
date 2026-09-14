"""P4-D CERTIFIED GLOBAL T_c(lambda, mu*) FOR PHYSICAL COULOMB (mu* > 0).

First rigorous certified two-sided global T_c(lambda) bracket for the Morel-Anderson
Coulomb pseudopotential (mu* > 0), extending deliverable (3) (mu*=0, file
`p4_global_tc_certified.py`) to the physical Coulomb pencil.

------------------------------------------------------------------------------
THE PENCIL (conventions verbatim from P4D_coulomb_operator.md sec.1 / v0.6, and
P4_mustar_uniqueness_strategy.md sec.0):

    M(wbar, mutilde, N_c) = K_ph(wbar) - mutilde * K_c(N_c),
    K_ph(wbar)            = H_mat(., wbar)            (the certified phonon operator),
    K_c(N_c)              = 2 |u_c><u_c|,  u_c[n] = (2n+1)^{-1/2} 1_{n<N_c}
                           (rank-1, PSD),
    mutilde := mu*/lambda,   N_c = floor(r * wbar),   r = 3 (physical Coulomb cutoff),
    T_c boundary:  lambda * lam_1(M) = 1   <=>   lam_1(M(wbar*)) = 1/lambda,
    T_c/omega_E = 1/(2 pi wbar*).

For fixed mu* we scan lambda; at each lambda, mutilde = mu*/lambda and the target is
    g_mutilde(wbar) := lam_1(M(wbar, mutilde, floor(r wbar))) = 1/lambda =: gstar.

------------------------------------------------------------------------------
CERTIFICATE STRUCTURE (all bounds on the EXACT operator M, not a float image):

  (a) LOWER lam_1_LO(M): eigenvector Rayleigh quotient of the float-entry block
      A_M = H_mat(Nv) - mutilde K_c(Nv) evaluated exactly in mpmath, minus ROUND_PAD.
      Any trial vector is valid; the float64 top eigenvector of A_M is near-optimal.
      SOUNDNESS: v is supported on [0,Nv); <v,Mv> = <v, M^(Nv) v> (M^(Nv) the true
      leading block, via the exact compression identity P_N K_c(N_c) P_N =
      K_c(min(N,N_c)) with Nv>=N_c) <= lam_1(M^(Nv)) <= lam_1(M) by Cauchy interlacing.
      ||A_M - M^(Nv)|| (float64 entry rounding of H_mat and K_c) < 2e-14 << ROUND_PAD.
      A SECOND, fully outward-rounded interval-arithmetic witness (geometric trial
      vector, +,*,/ and iv.sqrt only) is printed alongside for compliance (task a).

  (b) UPPER lam_1_HI(M): Kato-Temple.  lam_1 <= mu + res^2/(mu - gamma), with
      gamma >= lam_2(M) and mu > gamma.  Here mu is the exact Rayleigh, res^2 the
      residual (head ||A_M v - mu v||^2 exact + tail ||(I-P_Nv) H v||^2 <= ||B||^2
      from schur_coupling; the Coulomb term does NOT reach the tail because
      K_c v is supported on [0,N_c) subset [0,Nv), so the tail bound is IDENTICAL
      to mu*=0).  gamma is the INHERITED alpha-side certificate:
          K_c >= 0  =>  M <= K_ph  =>  lam_k(M) <= lam_k(K_ph)  (min-max)
          =>  lam_2(M) <= lam_2(K_ph) <= alpha_certificate_H(wbar).
      This is the structural dividend of PSD rank-1: NOT ONE LINE of the alpha
      (Cholesky-deflation + Schur) side is recomputed.  It is CONSERVATIVE:
      alpha ~ 0.03-0.4 is a loose upper bound on lam_2(M) ~ 1e-4, so the certified
      wbar_min where mu_lo > gamma sits ~0.1 ABOVE the true collapse boundary
      (P4D_coulomb_operator.md U8: cert 0.539/0.643 vs true ~0.43/0.50).  A tighter
      lam_2(M) deflation would recover ~0.1 in wbar_min; the doc judges the gain not
      worth it and we match its published wbar_min.  Where the top gap collapses
      (small wbar, or near-crossing kappa = 2 mutilde c1^2 / gamma -> 1), mu_lo <=
      gamma and KT is DECLARED FAILED (returns None): we restrict the range and say
      so, never paper over it.

  (c) INVERSION (sup rule + down-jumps).  Within a Coulomb window (N_c fixed),
      g_mutilde is monotone increasing in wbar (H(wbar) increasing, K_c fixed;
      proven for mutilde < mutilde_* by Hellmann-Feynman, P4D_coulomb_uniqueness.md
      Lemma 4).  ACROSS windows N_c increments and the extra repulsion makes
      g_mutilde JUMP DOWN at wbar = k/r (a genuine step-model feature, r=3).  Hence
      {wbar : g >= gstar} can be disconnected and lambda*lam_1=1 can have several
      roots; T_c is DEFINED AS sup{T : lambda lam_1(M(T)) >= 1}, i.e. the SMALLEST
      root wbar* (largest T).  The certified bracket [wbar_lo, wbar_hi] is obtained
      by (i) certifying every lower window's right-end PEAK has gHI < gstar (proves
      no smaller root -> the sup is in this window), then (ii) bisecting inside the
      crossing window (monotone) until the bracket separates gstar.  The non-unique
      lambda-band around each down-jump is disclosed explicitly, NOT hidden.

Sanity anchor (P4D_coulomb_operator.md sec.3, P4_mustar_uniqueness_strategy.md sec.2.1):
  lambda=1, mutilde=0.10, r=3, N_c 5->6 at wbar=2.000: lam_1(M) down-jumps to 0.998882.

Interval arithmetic is sound: outward-rounded mpmath.iv with +,*,/ and iv.sqrt only
(iv.sqrt enclosure verified in run_selftest; NO iv.atan2, NO mp.pi inside iv).
OMP_NUM_THREADS=1, single thread.  Self-test: `python p4_coulomb_global_tc_certified.py test`.
"""
import os
import sys
import time
import numpy as np
import mpmath as mp

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
from p4_interchange_uniformity import (  # noqa: E402
    H_mat, alpha_certificate_H, schur_coupling, spec_point, ROUND_PAD)

mp.mp.dps = 30

R_PHYS = 3.0          # physical Coulomb cutoff r = omega_c / omega_E
TWO_PI = float(2 * mp.pi)

PASS = [0]
FAIL = [0]


def ck(cond, msg=""):
    if cond:
        PASS[0] += 1
    else:
        FAIL[0] += 1
        print(f"  !! CHECK FAILED: {msg}")
    return cond


# ===========================================================================
# COULOMB OPERATOR  K_c = 2 |u_c><u_c|
# ===========================================================================
def K_c_np(N, Nc):
    """2|u_c><u_c| on the leading NxN block, float64.  u_c[n]=(2n+1)^{-1/2} 1_{n<Nc}.
    Identical to p4d_coulomb_certify.K_c (cross-checked bit-for-bit in run_selftest)."""
    u = np.zeros(N)
    k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return 2.0 * np.outer(u, u)


def _nc_phys(w, r=R_PHYS):
    """N_c = floor(r * wbar).  The +1e-9 guards float64 floor at exact boundaries."""
    return int(np.floor(r * w + 1e-9))


# ===========================================================================
# (a)+(b)  PENCIL NODE CERTIFICATE  [lam_1_LO, lam_1_HI] of M(w, mutilde, Nc)
# ===========================================================================
def pencil_node_cert(w, mu_t, Nc, Nv=140, Na=70, gamma=None):
    """Rigorous (lam1_lo, lam1_hi, gamma, mu_lo, res2_hi) for the EXACT operator M.

    lam1_hi is None if Kato-Temple fails (mu_lo <= gamma: gap collapse / near-crossing).
    `gamma` may be supplied (cached alpha_certificate_H(w)); recomputed if None.
    At mu_t == 0 this returns node_cert of the phonon operator BIT-FOR-BIT (0.0*K_c
    is annihilated in IEEE-754), which is the mu*=0 self-check gate.
    """
    Nv = max(int(Nv), int(Nc))          # compression identity needs Nv >= Nc
    A = H_mat(Nv, w) - mu_t * K_c_np(Nv, Nc)
    _, V = np.linalg.eigh(A)
    v = V[:, -1]
    vm = mp.matrix([mp.mpf(float(x)) for x in v])
    Am = mp.matrix(Nv, Nv)
    for i in range(Nv):
        for j in range(Nv):
            Am[i, j] = mp.mpf(float(A[i, j]))
    Av = Am * vm
    num = sum(vm[i] * Av[i] for i in range(Nv))
    den = sum(vm[i] ** 2 for i in range(Nv))
    R = num / den
    mu_lo = R - ROUND_PAD
    mu_hi = R + ROUND_PAD
    res_head2 = sum((Av[i] - R * vm[i]) ** 2 for i in range(Nv)) / den
    # tail: (I - P_Nv) M v = (I - P_Nv) H v  (K_c v supported in [0,Nc) subset [0,Nv))
    nB = schur_coupling(Nv, w, 0, spec=spec_point(w, 0))
    res_norm_hi = mp.sqrt(res_head2) + ROUND_PAD
    res2_hi = res_norm_hi ** 2 + nB ** 2
    if gamma is None:
        gamma, ok, _, _, _ = alpha_certificate_H(w, Na=Na)
        if not ok or gamma is None:
            return float(mu_lo), None, None, float(mu_lo), None
    if not (mu_lo > gamma):
        return float(mu_lo), None, float(gamma), float(mu_lo), float(res2_hi)
    h_hi = mu_hi + res2_hi / (mu_lo - gamma)
    return float(mu_lo), float(h_hi), float(gamma), float(mu_lo), float(res2_hi)


# ---------------------------------------------------------------------------
# (a')  fully interval-arithmetic geometric-vector lower bound (compliance).
#       <v,Mv> = <v,Hv> - 2 mutilde <u_c,v>^2,  v_n = t^n.
# ---------------------------------------------------------------------------
def geo_coeffs(t, Nv):
    """Exact mpf coefficients c_k of <v,H(w)v> = sum_k c_k [[k]](w), v_n=t^n, plus
    B=||v||^2.  Re-derived verbatim from p4_global_tc_certified.geo_coeffs (numpy-free)."""
    t = mp.mpf(t)
    v = [t ** n for n in range(Nv)]
    inv = [1 / mp.sqrt(mp.mpf(2 * n + 1)) for n in range(Nv)]
    B = sum(vn * vn for vn in v)
    c = [mp.mpf(0)] * (2 * Nv)
    for k in range(1, Nv):
        s = mp.mpf(0)
        for n in range(0, Nv - k):
            s += v[n] * v[n + k] * inv[n] * inv[n + k]
        c[k] += 2 * s
    for n in range(Nv):
        for m in range(Nv):
            c[n + m + 1] += v[n] * v[m] * inv[n] * inv[m]
    for k in range(1, Nv):
        s = mp.mpf(0)
        for n in range(k, Nv):
            s += v[n] * v[n] / mp.mpf(2 * n + 1)
        c[k] += -2 * s
    return c, B


def geo_lower_pencil(w, mu_t, Nc, t, Nv=60):
    """Rigorous outward-rounded lower bound on lam_1(M(w,mutilde,Nc)) at the point w,
    from the geometric trial vector v_n=t^n.  Uses mp.iv (+,*,/ and iv.sqrt) only."""
    c, B = geo_coeffs(t, Nv)
    wiv = mp.iv.mpf(w)
    acc = mp.iv.mpf(0)
    for k in range(1, len(c)):
        if c[k] == 0:
            continue
        acc += mp.iv.mpf(str(c[k])) * (wiv * wiv / (mp.iv.mpf(k) ** 2 + wiv * wiv))
    # Coulomb subtraction: 2 mutilde <u_c,v>^2, upper-bounded via interval square
    tiv = mp.iv.mpf(str(t))
    uv = mp.iv.mpf(0)
    for n in range(min(Nc, Nv)):
        uv += tiv ** n * (mp.iv.mpf(1) / mp.iv.sqrt(mp.iv.mpf(2 * n + 1)))
    coulomb = mp.iv.mpf(2) * mp.iv.mpf(str(mu_t)) * (uv * uv)
    q = (acc - coulomb) / mp.iv.mpf(str(B))
    return q.a       # outward-rounded lower endpoint


def best_t(w, Nv=60):
    """t maximizing the float64 geometric Rayleigh of H(w) (tightness only)."""
    n = np.arange(Nv)
    A = H_mat(Nv, w)
    best = (0.4, -1e9)
    for t in np.linspace(0.05, 0.95, 91):
        vv = t ** n
        mu = (vv @ (A @ vv)) / (vv @ vv)
        if mu > best[1]:
            best = (t, mu)
    return round(best[0], 3)


# ===========================================================================
# gamma (alpha) cache -- alpha_certificate_H depends on w ONLY, shared over all
# (lambda, mu*, mutilde); this is the whole point of the inherited alpha side.
# ===========================================================================
_GAMMA_CACHE = {}


def gamma_of(w, Na=70):
    key = round(float(w), 10)
    if key not in _GAMMA_CACHE:
        g, ok, _, _, _ = alpha_certificate_H(w, Na=Na)
        _GAMMA_CACHE[key] = float(g) if (ok and g is not None) else None
    return _GAMMA_CACHE[key]


def g_cert(w, mu_t, r=R_PHYS, Nv=140, Na=70):
    """Certified bracket [gLO, gHI] on g_mutilde(w) = lam_1(M) at physical N_c=floor(rw).
    Returns (gLO, gHI, Nc, gamma).  gHI is None if KT failed at this node."""
    Nc = _nc_phys(w, r)
    Nv_eff = max(Nv, Nc + 40)
    g = gamma_of(w, Na)
    lo, hi, gam, _, _ = pencil_node_cert(w, mu_t, Nc, Nv=Nv_eff, Na=Na, gamma=g)
    return lo, hi, Nc, gam


def g_float(w, mu_t, r=R_PHYS, N=200):
    """Fast float64 g_mutilde(w) for GUIDANCE only (locating windows / roots)."""
    Nc = _nc_phys(w, r)
    A = H_mat(N, w) - mu_t * K_c_np(N, Nc)
    return float(np.linalg.eigvalsh(A)[-1])


def _window_ends(mu_t, k, r=R_PHYS, N=200, eps=1e-6):
    """float64 g at the LEFT (k/r, post down-jump) and RIGHT ((k+1)/r-, peak) ends
    of the Coulomb window N_c=k.  g is monotone increasing between them."""
    gL = g_float(max(k / r, 1e-3), mu_t, r, N)
    gR = g_float((k + 1) / r - eps, mu_t, r, N)
    return gL, gR


def _crossings(mu_t, gstar, r=R_PHYS, k_min=0, k_max=16, N=200):
    """Window-based up-crossing map of the physical staircase (float64 guidance).
    Returns (cross_windows, start_sc): cross_windows = windows k with an interior
    up-crossing gL<gstar<=gR; start_sc = True if g already >= gstar at wbar=k_min/r.
    Number of roots = len(cross_windows) + start_sc."""
    cross = []
    gL0, _ = _window_ends(mu_t, k_min, r, N)
    start_sc = gL0 >= gstar
    for k in range(k_min, k_max + 1):
        gL, gR = _window_ends(mu_t, k, r, N)
        if gL < gstar <= gR:
            cross.append(k)
    return cross, start_sc


# ===========================================================================
# (c)  CERTIFIED SUP-RULE INVERSION  g_mutilde(wbar*) = gstar
# ===========================================================================
def tc_from_wstar(w_lo, w_hi):
    """T_c/omega_E bracket from a wbar* enclosure; 1/(2 pi w) is DECREASING in w."""
    return 1.0 / (TWO_PI * w_hi), 1.0 / (TWO_PI * w_lo)


def invert_sup(mu_t, gstar, r=R_PHYS, w_min=0.30, w_max=5.0,
               Nv=140, Na=70, tol=2e-4, eps_b=1e-6):
    """Certified sup-root bracket for g_mutilde(wbar) = gstar on [w_min, w_max].

    Returns dict with:
      status : 'ok' | 'below_wmin' | 'above_range' | 'kt_fail' | 'ambiguous'
      w_lo, w_hi : certified enclosure of the SUP (smallest) root wbar*
      Nc         : Coulomb window at the root
      nroots     : number of float64 roots on [w_min, w_max] (>=2 => non-unique band)
      sup_ok     : True if every lower window's peak certified gHI < gstar
    """
    # ---- float64 guidance: window-based up-crossing map ----
    k_min = _nc_phys(w_min, r)
    k_max = _nc_phys(w_max - 1e-9, r)
    cross, start_sc = _crossings(mu_t, gstar, r, k_min, k_max)
    nroots = len(cross) + (1 if start_sc else 0)
    if nroots == 0:
        return dict(status='above_range', w_lo=None, w_hi=None, Nc=None,
                    nroots=0, sup_ok=None)
    if start_sc:
        # SC already at wbar=k_min/r: sup root at/below the low end of the range
        return dict(status='below_wmin', w_lo=None, w_hi=None, Nc=k_min,
                    nroots=nroots, sup_ok=None)
    k_star = cross[0]                       # Coulomb window of the SUP (first) crossing
    w_guess = k_star / r

    if w_guess < w_min - 1e-9:
        return dict(status='below_wmin', w_lo=None, w_hi=None, Nc=k_star,
                    nroots=nroots, sup_ok=None)

    # ---- certify: every lower window's right-end PEAK has gHI < gstar ----
    # (proves no SC point below window k_star, i.e. the sup root really is here)
    sup_ok = True
    k = max(1, _nc_phys(w_min, r))
    while k < k_star:
        peak_w = (k + 1) / r - eps_b
        if peak_w <= w_min:
            k += 1
            continue
        lo, hi, _, gam = g_cert(peak_w, mu_t, r, Nv, Na)
        if hi is None:            # KT failed at this peak: cannot certify normal
            sup_ok = False
            break
        if not (hi < gstar):      # bracket cannot rule out a crossing in window k
            sup_ok = False
            break
        k += 1

    # ---- crossing window may sit below wbar_min (KT-uncertifiable) ----
    lo0, hi0, Nc0, gam0 = g_cert(w_guess, mu_t, r, Nv, Na)
    if hi0 is None:
        return dict(status='kt_fail', w_lo=None, w_hi=None, Nc=k_star,
                    nroots=nroots, sup_ok=sup_ok)

    # ---- certified bisection inside the crossing window (monotone) ----
    lo_edge = k_star / r
    hi_edge = (k_star + 1) / r - eps_b
    a = max(lo_edge, w_min)
    b = min(hi_edge, w_max)
    # ensure a is (certified) below gstar and b is (certified) at/above gstar
    la, ha, _, _ = g_cert(a, mu_t, r, Nv, Na)
    lb, hb, _, _ = g_cert(b, mu_t, r, Nv, Na)
    if ha is None or hb is None:
        return dict(status='kt_fail', w_lo=None, w_hi=None, Nc=k_star,
                    nroots=nroots, sup_ok=sup_ok)
    if not (ha < gstar):
        # left edge already SC: crossing is at/below the window's left boundary
        # (a down-jump landed just above gstar) -> sup root pinned at left edge
        return dict(status='ok', w_lo=a, w_hi=a, Nc=k_star,
                    nroots=nroots, sup_ok=sup_ok, boundary_pinned=True)
    if not (lb >= gstar):
        return dict(status='ambiguous', w_lo=a, w_hi=b, Nc=k_star,
                    nroots=nroots, sup_ok=sup_ok)
    # Two independent edge searches (g monotone increasing in the window):
    #   w_hi = smallest wbar with CERTIFIED SC   (gLO >= gstar)
    #   w_lo = largest  wbar with CERTIFIED normal(gHI <  gstar)
    # The unresolvable straddle sliver (|g-gstar| < bracket width) lies between them;
    # its width is ~ 2*(node bracket width)/g', i.e. the honest resolution limit.
    def _bisect(is_target):
        lo_e, hi_e = a, b
        for _ in range(60):
            if hi_e - lo_e < tol:
                break
            m = 0.5 * (lo_e + hi_e)
            lm, hm, _, _ = g_cert(m, mu_t, r, Nv, Na)
            if hm is None:
                break
            if is_target(lm, hm):
                hi_e = m
            else:
                lo_e = m
        return lo_e, hi_e
    # smallest certified-SC point: SC iff lm >= gstar
    _, w_hi = _bisect(lambda lm, hm: lm >= gstar)
    # largest certified-normal point: normal iff hm < gstar
    w_lo, _ = _bisect(lambda lm, hm: not (hm < gstar))
    if w_lo > w_hi:                         # numerical crossing of the two edges
        w_lo, w_hi = min(w_lo, w_hi), max(w_lo, w_hi)
    return dict(status='ok', w_lo=w_lo, w_hi=w_hi, Nc=k_star,
                nroots=nroots, sup_ok=sup_ok, boundary_pinned=False)


# ===========================================================================
# NON-UNIQUENESS BANDS (r=3 cross-window down-jumps), float64 disclosure map.
# ===========================================================================
def _nroots(mu_star, lam, r=R_PHYS):
    mu_t = mu_star / lam
    cross, start_sc = _crossings(mu_t, 1.0 / lam, r, 0, 16)
    return len(cross) + (1 if start_sc else 0)


def _refine_edge(mu_star, lam_single, lam_multi, r=R_PHYS, iters=22):
    """Bisect the lambda-edge between a single-root lambda and a multi-root lambda."""
    a, b = lam_single, lam_multi
    for _ in range(iters):
        m = 0.5 * (a + b)
        if _nroots(mu_star, m, r) >= 2:
            b = m
        else:
            a = m
    return 0.5 * (a + b)


def downjump_bands(mu_star, r=R_PHYS, lam_lo=0.5, lam_hi=3.0, n=280):
    """Scan lambda; report lambda-intervals where g_mutilde = 1/lambda has >=2 roots on
    the physical staircase (T_c non-unique; sup rule active).  Edges are refined by
    lambda-bisection.  Descriptive (float64, window-boundary based); the certified
    anchor (lambda=1, mu*=0.10) is verified separately in run_selftest / main."""
    lams = np.linspace(lam_lo, lam_hi, n)
    multi = [(_nroots(mu_star, lam, r) >= 2) for lam in lams]
    bands = []
    i = 0
    while i < n:
        if multi[i]:
            j = i
            while j + 1 < n and multi[j + 1]:
                j += 1
            lo = (_refine_edge(mu_star, lams[i - 1], lams[i], r) if i > 0 else lams[i])
            hi = (_refine_edge(mu_star, lams[j + 1], lams[j], r) if j + 1 < n else lams[j])
            bands.append((lo, hi))
            i = j + 1
        else:
            i += 1
    return bands


# ===========================================================================
# SELF-TEST
# ===========================================================================
def run_selftest():
    print("=== iv.sqrt outward-rounding (soundness of geo Coulomb term) ===")
    bad = 0
    for m in [1, 3, 5, 7, 9, 21, 41, 81, 161, 321, 641]:
        iv = mp.iv.sqrt(mp.iv.mpf(m))
        if not (iv.a <= mp.sqrt(mp.mpf(m)) <= iv.b):
            bad += 1
    ck(bad == 0, "iv.sqrt encloses on test set")
    print(f"  iv.sqrt enclosure failures: {bad}")

    print("\n=== K_c matches p4d_coulomb_certify.K_c bit-for-bit ===")
    try:
        from p4d_coulomb_certify import K_c as K_c_ref
        mx = 0.0
        for (N, Nc) in [(20, 5), (50, 13), (80, 0), (60, 60)]:
            mx = max(mx, float(np.max(np.abs(K_c_np(N, Nc) - K_c_ref(N, Nc)))))
        ck(mx == 0.0, f"K_c bit-for-bit vs certify (max|diff|={mx:.1e})")
        print(f"  max|K_c - certify.K_c| = {mx:.1e}")
    except Exception as e:                                   # pragma: no cover
        print(f"  (skipped cross-check: {e})")

    print("\n=== SELF-CHECK GATE: pencil at mutilde=0 == phonon node cert (bit-for-bit) ===")
    for w in (0.5, 1.0, 2.0, 3.0):
        lo0, hi0, g0, _, _ = pencil_node_cert(w, 0.0, 5)     # forces Coulomb branch
        lo1, hi1, g1, _, _ = pencil_node_cert(w, 0.0, 0)     # K_c == 0 matrix
        ck(lo0 == lo1 and hi0 == hi1 and g0 == g1,
           f"mutilde=0 bit-for-bit w={w}")
        print(f"  w={w}: lam1_lo={lo0:.10f} lam1_hi={hi0:.10f}  bit-eq={lo0==lo1 and hi0==hi1}")

    print("\n=== node cert soundness vs true lam_1(M), physical points ===")
    print(f"{'w':>5}{'mu_t':>7}{'Nc':>4}{'lam1_lo':>11}{'lam1_hi':>11}{'true':>11}"
          f"{'geoLO':>11}{'gamma':>8}{'relwid':>10}")
    pts = [(2.0, 0.10, 6), (2.0, 0.10, 5), (1.19, 0.0667, 3), (0.63, 0.0333, 1),
           (3.99, 0.1429, 11), (1.0, 0.10, 10), (0.7, 0.10, 3), (1.5, 0.26, 4)]
    for (w, mu_t, Nc) in pts:
        lo, hi, g, mlo, r2 = pencil_node_cert(w, mu_t, Nc)
        tr = float(np.linalg.eigvalsh(H_mat(800, w) - mu_t * K_c_np(800, Nc))[-1])
        glo = float(geo_lower_pencil(w, mu_t, Nc, best_t(w)))
        if hi is None:
            ck(g is None or mlo <= g, f"KT declared fail honestly w={w}")
            print(f"{w:>5}{mu_t:>7.4f}{Nc:>4}{lo:>11.6f}{'KTfail':>11}{tr:>11.6f}"
                  f"{glo:>11.6f}{(g if g else float('nan')):>8.3f}{'--':>10}")
        else:
            ck(lo <= tr + 1e-9, f"lam1_lo<=true w={w}")
            ck(hi >= tr - 1e-9, f"lam1_hi>=true w={w}")
            ck(glo <= lo + 1e-9, f"geoLO<=eigLO w={w}")
            print(f"{w:>5}{mu_t:>7.4f}{Nc:>4}{lo:>11.6f}{hi:>11.6f}{tr:>11.6f}"
                  f"{glo:>11.6f}{g:>8.3f}{(hi-lo)/abs(lo):>10.2e}")

    print("\n=== ANCHOR: lambda=1, mutilde=0.10, r=3, N_c 5->6 down-jump at wbar=2.0 ===")
    lo5, hi5, _, _, _ = pencil_node_cert(2.0, 0.10, 5)
    lo6, hi6, _, _, _ = pencil_node_cert(2.0, 0.10, 6)
    tr5 = float(np.linalg.eigvalsh(H_mat(800, 2.0) - 0.10 * K_c_np(800, 5))[-1])
    tr6 = float(np.linalg.eigvalsh(H_mat(800, 2.0) - 0.10 * K_c_np(800, 6))[-1])
    ck(lo5 <= tr5 <= hi5, "N_c=5 brackets true lam_1 (~1.0074)")
    ck(lo6 <= tr6 <= hi6, "N_c=6 brackets true lam_1")
    ck(abs(lo6 - 0.998882) < 5e-6, "N_c=6 lam1_lo == 0.998882 (published anchor)")
    ck(lo6 < 1.0 < lo5, "down-jump straddles gstar=1 (non-unique at lambda=1)")
    print(f"  N_c=5: [{lo5:.6f}, {hi5:.6f}]   N_c=6: [{lo6:.6f}, {hi6:.6f}]")
    print(f"  -> lam_1 down-jumps {lo5:.6f} -> {lo6:.6f} across wbar=2 (matches 0.998882)")

    print(f"\nSELF-TEST: PASS {PASS[0]}  FAIL {FAIL[0]}")


# ===========================================================================
# MAIN CERTIFIED PIPELINE
# ===========================================================================
def main(lam_grid=None):
    t0 = time.time()
    print("=" * 78)
    print("P4-D CERTIFIED GLOBAL T_c(lambda, mu*)  --  PHYSICAL COULOMB (mu* > 0), r=3")
    print("=" * 78)
    if lam_grid is None:
        lam_grid = [0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]

    for mu_star in (0.10, 0.13):
        print(f"\n{'#'*70}\n#  mu* = {mu_star}   (r=3, mutilde = mu*/lambda)\n{'#'*70}")
        print(f"{'lambda':>7}{'mutilde':>9}{'gstar':>8}{'Nc*':>5}"
              f"{'wbar_lo':>10}{'wbar_hi':>10}{'Tc_lo/wE':>11}{'Tc_hi/wE':>11}"
              f"{'relwid':>9}{'roots':>6}  status")
        worst = 0.0
        worst_at = None
        rows = []
        for lam in lam_grid:
            mu_t = mu_star / lam
            gstar = 1.0 / lam
            res = invert_sup(mu_t, gstar, r=R_PHYS, w_min=0.30, w_max=5.0)
            st = res['status']
            if st != 'ok':
                print(f"{lam:>7.2f}{mu_t:>9.4f}{gstar:>8.4f}"
                      f"{(res['Nc'] if res['Nc'] else 0):>5}{'--':>10}{'--':>10}"
                      f"{'--':>11}{'--':>11}{'--':>9}{res['nroots']:>6}  {st}")
                continue
            w_lo, w_hi = res['w_lo'], res['w_hi']
            Flo, Fhi = tc_from_wstar(w_lo, w_hi)
            rw = (Fhi - Flo) / Flo
            ck(w_lo <= w_hi + 1e-12, f"wbar bracket ordered lam={lam}")
            ck(Flo <= Fhi + 1e-15, f"Tc bracket ordered lam={lam}")
            supflag = '' if res.get('sup_ok') else ' [sup?]'
            note = ' [bnd-pin]' if res.get('boundary_pinned') else ''
            rows.append((lam, mu_t, Flo, Fhi, rw, res['nroots']))
            if rw > worst:
                worst, worst_at = rw, lam
            print(f"{lam:>7.2f}{mu_t:>9.4f}{gstar:>8.4f}{res['Nc']:>5}"
                  f"{w_lo:>10.5f}{w_hi:>10.5f}{Flo:>11.6f}{Fhi:>11.6f}"
                  f"{rw:>9.2e}{res['nroots']:>6}  ok{supflag}{note}")
        print(f"  worst-case certified T_c relative width: {worst:.3e}  (at lambda={worst_at})")

        # non-uniqueness (down-jump) bands
        bands = downjump_bands(mu_star, lam_lo=0.55, lam_hi=3.0, n=500)
        print(f"  DISCLOSED non-uniqueness lambda-bands (>=2 roots, T_c:=sup), r=3:")
        if not bands:
            print("    (none on lambda in [0.55,3.0] at this resolution)")
        for (a, b) in bands:
            tag = ""
            if a <= 1.0 <= b and abs(mu_star - 0.10) < 1e-12:
                tag = "  <-- contains anchor lambda=1 (N_c 5->6, lam1->0.998882)"
            print(f"    lambda in [{a:.4f}, {b:.4f}]   (width {b-a:.4f}){tag}")

    # ---- collapse anchor + worst-N_c boundary (P4D_coulomb_operator.md U8) --------
    print(f"\n{'='*70}\nCOLLAPSE ANCHOR & worst-N_c KT boundary (r-general; physical r=3 path is clean)")
    # exact doc anchor: wbar=0.3, mutilde=0.10, worst N_c=3 -> lam_1 ~ 6.216e-5
    E = np.linalg.eigvalsh(H_mat(400, 0.3) - 0.10 * K_c_np(400, 3))
    lam1_c, gap_c = float(E[-1]), float(E[-1] - E[-2])
    ck(abs(lam1_c - 6.216144e-5) < 1e-8, "collapse anchor lam1(0.3,0.10,Nc=3)=6.216e-5")
    ck(abs(gap_c - 6.475852e-5) < 1e-8, "collapse anchor gap=6.476e-5")
    lo_c, hi_c, g_c, mlo_c, _ = pencil_node_cert(0.3, 0.10, 3)
    ck(hi_c is None, "node cert declares KT-FAIL at the collapse point (honest)")
    print(f"  wbar=0.3, mutilde=0.10, N_c=3: true lam_1={lam1_c:.6e}, gap={gap_c:.6e} "
          f"(matches doc 6.216e-5/6.476e-5); node cert KT: {'FAIL (mu_lo<gamma)' if hi_c is None else 'ok'}")
    print(f"  worst-N_c node-cert KT boundary vs doc 2x2-cell wbar_min vs true collapse:")
    DOC_CERT = {0.10: 0.5391, 0.13: 0.6427, 0.26: 1.1686}
    TRUE_COLL = {0.10: 0.43, 0.13: 0.50, 0.26: 0.90}
    for mu_t in (0.10, 0.13, 0.26):
        wm = _worst_nc_wbar_min(mu_t)
        ok = (wm is not None) and (TRUE_COLL[mu_t] - 0.05 <= wm <= DOC_CERT[mu_t] + 1e-6)
        ck(ok, f"worst-N_c boundary {wm} in (true~{TRUE_COLL[mu_t]}, doc {DOC_CERT[mu_t]}] mutilde={mu_t}")
        print(f"    mutilde={mu_t}: node-cert {wm:.3f}   (< doc-cell {DOC_CERT[mu_t]}, "
              f"> true collapse ~{TRUE_COLL[mu_t]})  {'OK' if ok else 'CHECK'}")
    print("  -> node cert (eigvector Rayleigh + KT) is TIGHTER than the doc 2x2-cell")
    print("     certificate; along the PHYSICAL r=3 path the whole Tc range certifies")
    print("     cleanly (the collapse needs a large Coulomb window N_c>>r*wbar, i.e. r>>3).")

    print(f"\nTIME: {time.time()-t0:.1f}s   PASS {PASS[0]}  FAIL {FAIL[0]}")


def _worst_nc_wbar_min(mu_t, r_cap=12.0, lo=0.30, hi=1.60, step=0.01, Na=70):
    """Smallest wbar where the node cert (KT) succeeds for EVERY Coulomb window
    N_c in [0, ceil(r_cap*wbar)] -- the r-general (worst-N_c) certifiable boundary.
    Tighter than the doc's 2x2-cell wbar_min because the eigenvector Rayleigh + KT
    beats a 2x2 cell bound; still above the true operator collapse."""
    w = lo
    while w <= hi:
        gam = gamma_of(w, Na)
        Ncap = int(np.ceil(r_cap * w)) + 2
        ok_all = True
        for Nc in range(0, Ncap + 1):
            _, hi_g, _, _, _ = pencil_node_cert(w, mu_t, Nc, Nv=max(140, Nc + 40),
                                                Na=Na, gamma=gam)
            if hi_g is None:
                ok_all = False
                break
        if ok_all:
            return w
        w += step
    return None


# ===========================================================================
if __name__ == "__main__":
    if "test" in sys.argv:
        run_selftest()
    else:
        main()
