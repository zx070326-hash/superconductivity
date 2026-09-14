"""P4 CERTIFIED CROSSOVER -- certified two-sided envelope h_LO(w) <= h(w) <= h_HI(w)
on the crossover window w in [0.15, 3], and its inversion to a certified global
T_c(lambda)/omega_E bracket for the Einstein spectrum (mu*=0).

h(w) := lam_1(H(w)),  H(w) = O[[[.]](w)],  [[k]](w) = w^2/(k^2+w^2),  w = w_E/2pi T.
T_c criterion (mu*=0): lambda*h(w*) = 1, T_c/w_E = 1/(2pi w*).

CERTIFICATE STRUCTURE (all bounds are on the EXACT limit operator, not a float image):
  (a) LOWER h_LO(w): the larger of two rigorous Rayleigh lower bounds
        * GEOMETRIC v_n = t^n (t rational, fixed per node): <v,H(w)v>/<v,v> computed
          by EXACT coefficient decomposition <v,H(w)v> = sum_k c_k [[k]](w), evaluated
          with OUTWARD-rounded interval arithmetic (mpmath.iv).  [task (a)]
        * EIGENVECTOR (rayleigh_mu): float64 top eigenvector as a finite-support trial
          vector, exact mpmath Rayleigh quotient minus ROUND_PAD (any vector is valid;
          this one is near-optimal).  Reused verbatim from p4_interchange_uniformity.
      Both are <= h; their max is <= h.
  (b) UPPER h_HI(w): Kato-Temple.  With a certified gamma >= lam_2(H(w)) and the
      exact Rayleigh mu, h <= mu + res^2/(mu - gamma), res = ||H v - mu v||_l2.
      gamma from alpha_certificate_H (Cholesky-deflation + Schur tail).  The residual
      tail ||(I-P_Nv) H v|| is bounded by the Schur coupling norm.  The enormous
      spectral gap (lam_2 ~ 1e-3..0.2 << mu ~ 0.03..1.7) makes KT near-lossless.
  (c) CONTINUOUS w-coverage: h is PROVEN strictly increasing (P4_main_theorems.md
      Cor (i), via Lemma M).  Hence on a cell [a,b]:  h_LO(a) <= h(w) <= h_HI(b)
      for every w in [a,b] -- monotone staircase, ZERO slack.  (The Schur/Lipschitz
      w-envelope of p4_interchange_uniformity is the fallback if monotonicity were
      unavailable; here it is owned, so we use it and pay nothing.)

Cross-checked against Theorem W (K=8 Neumann bracket) on the overlap [2,3].
No scratch imports from /tmp.  Single-thread, OMP_NUM_THREADS=1, < 10 min.
"""
import os
import sys
import time
import numpy as np
import mpmath as mp

ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
from p4_interchange_uniformity import (  # noqa: E402
    H_mat, rayleigh_mu, alpha_certificate_H, schur_coupling,
    spec_point, ROUND_PAD)
# NOTE: p4_weak_coupling_series runs its full numerical GATE at module level
# (lines ~116-275 execute on import, ~1-2 min up to w=200), so we do NOT import
# it.  pieces_M1 / F_partial are re-derived here verbatim from that module's
# algorithm (they only use numpy, no H_mat) for the Theorem W [2,3] cross-check.

mp.mp.dps = 30

# ---------------------------------------------------------------------------
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
# (a) GEOMETRIC TRIAL VECTOR:  <v,H(w)v> = sum_{k>=1} c_k [[k]](w),  v_n = t^n.
# ===========================================================================
def geo_coeffs(t, Nv):
    """Exact coefficients c_k (mpf) of <v,H(w)v> = sum_{k>=1} c_k [[k]](w) for the
    finite-support geometric trial vector v_n = t^n, n=0..Nv-1.  Also returns B=||v||^2.

    Derivation (H(w) entries, 0-indexed):
      off-diag n!=m : (v_n v_m)(1/sqrt((2n+1)(2m+1))) [[ |n-m| ]] + [[ n+m+1 ]]
      diag n==m     : v_n^2 ( [[2n+1]]/(2n+1) - (2/(2n+1)) sum_{j=1}^n [[j]] )
    Collecting the coefficient of [[k]]:
      c_k = 2 sum_{n>=0} v_n v_{n+k} / sqrt((2n+1)(2(n+k)+1))            (|n-m|=k)
          +   sum_{n+m=k-1}  v_n v_m / sqrt((2n+1)(2m+1))               (n+m+1=k)
          - 2 sum_{n>=k} v_n^2 / (2n+1)                                  (diag correction)
    """
    t = mp.mpf(t)
    v = [t**n for n in range(Nv)]
    inv = [1 / mp.sqrt(mp.mpf(2 * n + 1)) for n in range(Nv)]
    B = sum(vn * vn for vn in v)
    kmax = 2 * Nv - 1
    c = [mp.mpf(0)] * (kmax + 1)
    # off-diagonal |n-m|=k term (k>=1)
    for k in range(1, Nv):
        s = mp.mpf(0)
        for n in range(0, Nv - k):
            s += v[n] * v[n + k] * inv[n] * inv[n + k]
        c[k] += 2 * s
    # n+m+1=k term  (k = n+m+1, n,m in [0,Nv-1])
    for n in range(Nv):
        for m in range(Nv):
            k = n + m + 1
            c[k] += v[n] * v[m] * inv[n] * inv[m]
    # diagonal correction: -2 sum_{n>=k} v_n^2/(2n+1) as coeff of [[k]], k=1..Nv-1... actually k up to Nv-1? j runs 1..n, n<=Nv-1 => k<=Nv-1
    # coefficient of [[k]] is -2 sum_{n>=k, n<=Nv-1} v_n^2/(2n+1)
    for k in range(1, Nv):
        s = mp.mpf(0)
        for n in range(k, Nv):
            s += v[n] * v[n] / mp.mpf(2 * n + 1)
        c[k] += -2 * s
    return c, B


def brk_iv(k, wiv):
    """[[k]](w) = w^2/(k^2+w^2) as an outward-rounded interval for w in wiv (iv.mpf).
    Monotone increasing in w, and iv arithmetic here is +,*,/ only (no atan2/pi),
    so this is a sound enclosure."""
    k2 = mp.iv.mpf(k) ** 2
    return wiv * wiv / (k2 + wiv * wiv)


def geo_lower_point(w, t, Nv):
    """Rigorous lower bound on h(w): outward-rounded <v,H(w)v>/<v,v> at the point w."""
    c, B = geo_coeffs(t, Nv)
    wiv = mp.iv.mpf(w)
    acc = mp.iv.mpf(0)
    for k in range(1, len(c)):
        if c[k] == 0:
            continue
        ck_iv = mp.iv.mpf(str(c[k]))
        acc += ck_iv * brk_iv(k, wiv)
    q = acc / mp.iv.mpf(str(B))
    return q.a  # lower endpoint (outward-rounded down)


def geo_lower_cell(a, b, t, Nv):
    """EXACT rigorous lower bound on h(w) valid for EVERY w in [a,b]:
    each [[k]] is increasing in w and appears once, so min over [a,b] of the
    affine-in-{[[k]]} combination is attained termwise.  Outward-rounded."""
    c, B = geo_coeffs(t, Nv)
    wa = mp.iv.mpf(a)
    wb = mp.iv.mpf(b)
    acc = mp.iv.mpf(0)
    for k in range(1, len(c)):
        if c[k] == 0:
            continue
        ck_iv = mp.iv.mpf(str(c[k]))
        # minimize c_k [[k]](w) over [a,b]:  [[k]] increasing, so use a if c_k>0 else b
        if c[k] > 0:
            acc += ck_iv * brk_iv(k, wa)
        else:
            acc += ck_iv * brk_iv(k, wb)
    q = acc / mp.iv.mpf(str(B))
    return q.a


def best_t(w, Nv):
    """Pick t maximizing the (float64) geometric Rayleigh quotient at w (tightness
    only; ANY t gives a valid bound).  Returns t rounded to 3 decimals (rational)."""
    n = np.arange(Nv)
    A = H_mat(Nv, w)
    best = (0.4, -1e9)
    for t in np.linspace(0.05, 0.95, 91):
        v = t ** n
        mu = (v @ (A @ v)) / (v @ v)
        if mu > best[1]:
            best = (t, mu)
    t0 = best[0]
    for t in np.linspace(max(0.02, t0 - 0.02), min(0.98, t0 + 0.02), 41):
        v = t ** n
        mu = (v @ (A @ v)) / (v @ v)
        if mu > best[1]:
            best = (t, mu)
    return round(best[0], 3)


# ===========================================================================
# (b) EIGENVECTOR NODE CERTIFICATE: rigorous [h_lo, h_hi] via Rayleigh + Kato-Temple.
# ===========================================================================
def node_cert(w, Nv=140, Na=70):
    """Return (h_lo_eig, h_hi, gamma, mu_lo, res2_hi) rigorous for the EXACT operator."""
    A = H_mat(Nv, w)
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
    R = num / den                         # Rayleigh of the float-entry matrix
    mu_lo = R - ROUND_PAD                 # lower bound on true Rayleigh <= h
    mu_hi = R + ROUND_PAD                 # upper bound on true Rayleigh
    # residual head: ||A v - R v||^2 / ||v||^2 (exact, float entries)
    res_head2 = sum((Av[i] - R * vm[i]) ** 2 for i in range(Nv)) / den
    # residual tail: ||(I-P_Nv) H v||/||v|| <= ||coupling block|| = schur_coupling
    nB = schur_coupling(Nv, w, 0, spec=spec_point(w, 0))
    res_norm_hi = mp.sqrt(res_head2) + ROUND_PAD
    res2_hi = res_norm_hi ** 2 + nB ** 2
    # gamma >= lam_2 (certified)
    gamma, ok, c, _, _ = alpha_certificate_H(w, Na=Na)
    if not ok or gamma is None:
        return float(mu_lo), None, None, float(mu_lo), None
    if not (mu_lo > gamma):
        return float(mu_lo), None, float(gamma), float(mu_lo), float(res2_hi)
    h_hi = mu_hi + res2_hi / (mu_lo - gamma)
    return float(mu_lo), float(h_hi), float(gamma), float(mu_lo), float(res2_hi)


# ===========================================================================
# THEOREM W (K) bracket -- reused for the [2,3] cross-check.
# ===========================================================================
THETA_STAR = 0.6349287917881079


def pieces_M1(N, w):
    """tau, calR, D from the (M1) Lorentz-addition identity.  Re-derived verbatim
    from p4_weak_coupling_series.pieces_M1 (numpy-only, no H_mat)."""
    n = np.arange(N)
    wn = 1.0 / (2 * n + 1.0)
    p = np.sqrt(wn)
    s = (n + 0.5) / w
    G = 1.0 / (1.0 + s ** 2)
    tau = p * G
    S, Sg = s[:, None], s[None, :]
    R = 2 * S**2 * Sg**2 * (5 + S**2 + Sg**2)
    R /= (1 + (S - Sg)**2) * (1 + (S + Sg)**2) * (1 + S**2) * (1 + Sg**2)
    calR = np.outer(p, p) * R
    k = np.arange(1, N)
    Gk = 1.0 / (1.0 + (k / w) ** 2)
    lG = np.concatenate(([0.0], np.cumsum(Gk))) * (2.0 / (2 * n + 1.0))
    D = wn + lG
    return tau, calR, D


def F_partial(h, tau, calR, D, K_max):
    """Partial Neumann sum S_K(h) and its tail bound (Theorem W).  Re-derived
    verbatim from p4_weak_coupling_series.F_partial."""
    B = np.dot(tau, tau)
    A_inv_diag = 1.0 / (h + D)
    S_K = 2 * np.dot(tau**2, A_inv_diag)
    if K_max >= 1:
        v = A_inv_diag * tau
        for _ in range(K_max):
            v = A_inv_diag * (calR @ v)
            S_K += 2 * np.dot(tau, v)
    tail = 2 * B * THETA_STAR**(K_max + 1) / (h**(K_max + 1) * (h - THETA_STAR))
    return S_K, tail


def secular_bounds_K(w, N, K):
    """h_lo^(K) <= h(w) <= h_hi^(K) from the Neumann secular bracket (Theorem W)."""
    tau, calR, D = pieces_M1(N, w)
    B = float(np.dot(tau, tau))

    def f_lo(h):
        S_K, _ = F_partial(h, tau, calR, D, K)
        return S_K - 1.0

    def f_hi(h):
        S_K, tail = F_partial(h, tau, calR, D, K)
        return S_K + tail - 1.0

    def bisect(f, a, b):
        for _ in range(80):
            m = 0.5 * (a + b)
            if f(m) > 0:
                a = m
            else:
                b = m
        return 0.5 * (a + b)

    hi_bracket = 2 * B - 1 + 2.0
    h_lo = bisect(f_lo, THETA_STAR + 1e-6, hi_bracket)
    h_hi = bisect(f_hi, THETA_STAR + 1e-6, hi_bracket)
    return h_lo, h_hi


# ===========================================================================
# SMALL-GRID VALIDATION
# ===========================================================================
def run_selftest():
    print("=== SELF-TEST: geo coefficient decomposition vs v@H@v ===")
    for w in (0.15, 0.5, 1.0, 2.0, 3.0):
        Nv = 40
        t = best_t(w, Nv)
        c, B = geo_coeffs(t, Nv)
        # reconstruct <v,H(w)v> from coeffs
        acc = mp.mpf(0)
        for k in range(1, len(c)):
            acc += c[k] * (mp.mpf(w) ** 2 / (mp.mpf(k) ** 2 + mp.mpf(w) ** 2))
        n = np.arange(Nv)
        v = t ** n
        A = H_mat(Nv, w)
        direct = float(v @ (A @ v))
        rel = abs(float(acc) - direct) / abs(direct)
        ck(rel < 1e-10, f"geo coeff decomposition w={w}: rel={rel:.2e}")
        print(f"  w={w}: t={t} <v,Hv> coeff={float(acc):.8f} direct={direct:.8f} rel={rel:.2e}")

    print("\n=== SELF-TEST: node cert + geo lower, small grid ===")
    print(f"{'w':>6} {'h_true':>9} {'geoLO':>9} {'eigLO':>9} {'h_HI':>9} "
          f"{'gamma':>9} {'relwid':>9}")
    for w in (0.15, 0.5, 1.0, 2.0, 3.0):
        ht = float(np.linalg.eigvalsh(H_mat(600, w))[-1])
        l2 = float(np.linalg.eigvalsh(H_mat(600, w))[-2])
        Nv = 60
        t = best_t(w, Nv)
        glo = geo_lower_point(w, t, Nv)
        elo, hhi, g, mlo, r2 = node_cert(w)
        ck(float(glo) <= ht + 1e-12, f"geoLO<=h_true w={w}")
        ck(elo <= ht + 1e-12, f"eigLO<=h_true w={w}")
        ck(hhi is not None and hhi >= ht - 1e-12, f"h_HI>=h_true w={w}")
        ck(g >= l2 - 1e-9, f"gamma>=lam2 w={w}")
        relwid = (hhi - max(float(glo), elo)) / max(float(glo), elo)
        print(f"{w:>6.2f} {ht:>9.5f} {float(glo):>9.5f} {elo:>9.5f} {hhi:>9.5f} "
              f"{g:>9.4f} {relwid:>9.2e}")

    print("\n=== SELF-TEST: Theorem W K=8 on [2,3] ===")
    for w in (2.0, 2.5, 3.0):
        N = 40 * int(np.ceil(w))
        hlo8, hhi8 = secular_bounds_K(w, N, 8)
        ht = float(np.linalg.eigvalsh(H_mat(max(400, N), w))[-1])
        ck(hlo8 <= ht <= hhi8, f"W K=8 brackets h_true w={w}")
        print(f"  w={w}: W8=[{hlo8:.7f},{hhi8:.7f}] width={hhi8-hlo8:.2e} h_true={ht:.7f}")

    print(f"\nSELF-TEST: PASS {PASS[0]}  FAIL {FAIL[0]}")


# ===========================================================================
# STRONG-COUPLING END (Weyl envelope, PROVEN analytic; P4_main_theorems.md Cor (ii))
#   g(2) w^2 - 4 zeta(4) w^4  <=  h(w)  <=  g(2) w^2 ,  valid for all w
#   (upper strict via Lemma M; lower via ||O[1/k^2 - 1/(k^2+w^2)]|| <= 4 zeta(4) w^2).
# We use the certified g(2) bracket [G2_LO, G2_HI] (CERT-B).
# ===========================================================================
G2_LO = mp.mpf("1.3181093293757127")   # certified lower on g(2)=lam_1(O[1/k^2]) (CERT-B)
G2_HI = mp.mpf("1.3181615784372800")   # certified upper on g(2) (CERT-B)
Z4 = mp.zeta(4)             # pi^4/90, exact


def strong_h_bracket(w):
    """[h_lo, h_hi] on h(w) from the proven Weyl/strong envelope (mpf)."""
    w = mp.mpf(w)
    h_hi = G2_HI * w**2
    h_lo = G2_LO * w**2 - 4 * Z4 * w**4
    return h_lo, h_hi


def strong_invert(u):
    """Given u = 1/lambda, return [w_lo, w_hi] enclosing w* with h(w*)=u, using the
    strong envelope.  w* >= w where h_hi(w) <= u  (w_lo = sqrt(u/G2_HI));
    w* <= w where h_lo(w) >= u  (smallest positive root of G2_LO w^2 - 4Z4 w^4 = u).
    Returns None if u outside the reliable strong range."""
    u = mp.mpf(u)
    w_lo = mp.sqrt(u / G2_HI)
    # h_lo(w) = G2_LO x - 4 Z4 x^2 with x=w^2; solve = u -> 4Z4 x^2 - G2_LO x + u = 0
    disc = G2_LO**2 - 16 * Z4 * u
    if disc <= 0:
        return None  # w^4 term no longer controllable; outside strong regime
    x_up = (G2_LO - mp.sqrt(disc)) / (8 * Z4)   # smaller root (physical, small w)
    if x_up <= 0:
        return None
    w_hi = mp.sqrt(x_up)
    return w_lo, w_hi


# ===========================================================================
# WEAK-COUPLING END (Theorem W finite-N Neumann bracket, valid w >= 3).
# ===========================================================================
def weak_h_bracket(w, K=8):
    """[h_lo^(K), h_hi^(K)] on h(w) from Theorem W (w>=2.5). Matrix size grows with w
    (capped for cost; the Neumann tail bound, not N, controls the certified width)."""
    N = min(1200, max(400, 100 * int(np.ceil(w))))
    return secular_bounds_K(w, N, K)


# ===========================================================================
# MONOTONE-STAIRCASE INVERSION of a node envelope to a T_c bracket.
# ===========================================================================
def invert_envelope(ws, hlo, hhi, u):
    """Given monotone-increasing node arrays and u=1/lambda, return [w_lo, w_hi]
    enclosing w* (h(w*)=u) via the PROVEN monotonicity of h:
       w_lo = max{ ws[i] : hhi[i] <= u }   (h(ws[i]) <= hhi[i] <= u => w* >= ws[i])
       w_hi = min{ ws[i] : hlo[i] >= u }   (h(ws[i]) >= hlo[i] >= u => w* <= ws[i])
    Returns None if u is outside [hhi[0], hlo[-1]] (not covered by this grid)."""
    ws = np.asarray(ws, float)
    hlo = np.asarray(hlo, float)
    hhi = np.asarray(hhi, float)
    lo_ok = np.where(hhi <= u)[0]
    hi_ok = np.where(hlo >= u)[0]
    if len(lo_ok) == 0 or len(hi_ok) == 0:
        return None
    w_lo = ws[lo_ok[-1]]
    w_hi = ws[hi_ok[0]]
    if w_hi < w_lo:
        return None
    return w_lo, w_hi


TWO_PI = float(2 * mp.pi)


def tc_from_w(w_lo, w_hi):
    """T_c/omega_E bracket from a w* enclosure: 1/(2 pi w) is DECREASING in w."""
    return 1.0 / (TWO_PI * w_hi), 1.0 / (TWO_PI * w_lo)   # (F_lo, F_hi)


# ===========================================================================
# MAIN CERTIFIED PIPELINE
# ===========================================================================
def main(N_cross=1000, Nv=140, Na=70):
    t_start = time.time()
    print("=" * 78)
    print("P4 CERTIFIED GLOBAL T_c(lambda)  (Einstein spectrum, mu*=0)")
    print("=" * 78)

    # ---- crossover h-envelope on [0.15, 3] ------------------------------
    print(f"\n[1] Crossover h-envelope: {N_cross} geometric nodes on w in [0.15,3]")
    print(f"    node cert = eigenvector Rayleigh (lower) + Kato-Temple (upper), "
          f"Nv={Nv}, Na={Na}")
    ws = np.geomspace(0.15, 3.0, N_cross)
    hlo = np.zeros(N_cross)
    hhi = np.zeros(N_cross)
    gam = np.zeros(N_cross)
    bad = 0
    for i, w in enumerate(ws):
        elo, hh, g, mlo, r2 = node_cert(float(w), Nv=Nv, Na=Na)
        if hh is None:
            bad += 1
            hhi[i] = np.inf
            hlo[i] = elo
            gam[i] = g if g is not None else np.nan
            continue
        hlo[i] = elo
        hhi[i] = hh
        gam[i] = g
        # soundness at each node: lower<upper, gamma below the Rayleigh
        ck(hlo[i] < hhi[i], f"node {i} w={w:.4f}: hlo<hhi")
        ck(gam[i] < mlo, f"node {i} w={w:.4f}: gamma<mu")
    ck(bad == 0, f"all crossover nodes produced a finite Kato-Temple upper (bad={bad})")
    # monotonicity of the certified arrays (h proven increasing => must hold)
    ck(np.all(np.diff(hlo) > 0), "hlo strictly increasing along grid")
    ck(np.all(np.diff(hhi) > 0), "hhi strictly increasing along grid")
    node_relwid = (hhi - hlo) / hlo
    print(f"    node bracket relative width: max={node_relwid.max():.2e}  "
          f"min={node_relwid.min():.2e}  at w in "
          f"[{ws[node_relwid.argmax()]:.3f}, {ws[node_relwid.argmin()]:.3f}]")
    print(f"    h range: [{hlo[0]:.5f}, {hlo[-1]:.5f}]  =>  crossover covers "
          f"lambda in [{1/hlo[-1]:.4f}, {1/hhi[0]:.4f}]")

    # ---- cross-check against Theorem W (K=8) on the overlap [2,3] --------
    print("\n[2] Cross-check vs Theorem W (K=8) on the overlap w in [2,3]")
    for wc in (2.0, 2.5, 3.0):
        # evaluate BOTH brackets at the SAME grid node w = ws[j]
        j = int(np.argmin(np.abs(ws - wc)))
        wn = float(ws[j])
        w8lo, w8hi = weak_h_bracket(wn, K=8)
        my_lo, my_hi = hlo[j], hhi[j]
        # consistency: the two brackets must overlap (both contain h(wn))
        overlap = (my_lo <= w8hi) and (w8lo <= my_hi)
        ck(overlap, f"crossover & W8 brackets overlap at w={wn:.4f}")
        print(f"    w={wn:.4f}: crossover=[{my_lo:.6f},{my_hi:.6f}]  "
              f"W8=[{w8lo:.6f},{w8hi:.6f}]  overlap={overlap}")

    # ---- task-(a) compliance: geometric interval-Rayleigh at sample nodes -
    print("\n[3] Task-(a) geometric interval-Rayleigh lower bound (compliance check)")
    print(f"    {'w':>6} {'geoLO(iv)':>12} {'eigLO':>12} {'h_true':>12} {'geoLO<=eig':>11}")
    for wc in (0.15, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0):
        Nvg = 60
        t = best_t(wc, Nvg)
        glo = geo_lower_point(wc, t, Nvg)          # rigorous lower at EXACTLY wc
        elo, hh, g, mlo, r2 = node_cert(wc, Nv=Nv, Na=Na)  # eig lower at EXACTLY wc
        ht = float(np.linalg.eigvalsh(H_mat(600, wc))[-1])
        okc = float(glo) <= elo + 1e-9 and elo <= ht + 1e-12 and float(glo) <= ht + 1e-12
        ck(okc, f"geoLO<=eigLO<=h_true at w={wc}")
        print(f"    {wc:>6.2f} {float(glo):>12.6f} {elo:>12.6f} {ht:>12.6f} {str(okc):>11}")
    print("    -> geometric interval-Rayleigh is a rigorous but looser lower bound;")
    print("       the eigenvector Rayleigh (eigLO) dominates and is used as h_LO.")

    # ---- INVERSION to a certified T_c(lambda) table ---------------------
    print("\n[3.5] Validation: weak-W (finite-N) & asymptotic brackets vs h_true")
    for wv in (3.0, 5.0, 8.0):
        a, b = weak_h_bracket(wv, K=6)
        Nbig = min(3000, max(2000, 250 * int(np.ceil(wv))))
        ht = float(np.linalg.eigvalsh(H_mat(Nbig, wv))[-1])
        ck(a <= ht <= b, f"weak-W brackets h_true at w={wv}")
        print(f"    weak-W  w={wv}: [{a:.6f},{b:.6f}] contains h_true={ht:.6f}  "
              f"{a <= ht <= b}")
    for wv in (8.0, 20.0):
        L = np.log(wv)
        alo = _h_lo_proven(wv)                       # PROVEN (Cor iii)
        ahi = L + _C0 + _B1_HI / L + _REM / L**2      # asymptotic upper (labelled)
        Nbig = min(4000, max(2000, 200 * int(np.ceil(wv))))
        ht = float(np.linalg.eigvalsh(H_mat(Nbig, wv))[-1])
        ck(alo <= ht <= ahi, f"deep-tail [Cor(iii),asympt] brackets h_true at w={wv}")
        print(f"    deep    w={wv}: [{alo:.6f},{ahi:.6f}] contains h_true={ht:.6f}  "
              f"{alo <= ht <= ahi}")

    print("\n[4] Certified inversion  lambda -> [F_lo, F_hi] = T_c/omega_E bracket")
    lam_grid = sorted(set(
        [0.05, 0.07, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.582,
         0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 34.0,
         40.0, 50.0]))
    print(f"    {'lambda':>8} {'F_lo':>13} {'F_hi':>13} {'rel_wid':>9} {'method':>12}")
    rows = []
    cross_max = 0.0
    strong_max = 0.0
    weak_max = 0.0        # certified weak-W (finite-N) worst case
    for lam in lam_grid:
        u = 1.0 / lam
        method = None
        Flo = Fhi = None
        res = invert_envelope(ws, hlo, hhi, u)
        if res is not None:
            method = "crossover"
            w_lo, w_hi = res
            Flo, Fhi = tc_from_w(w_lo, w_hi)
        elif u < hhi[0]:        # h needed < h(0.15): strong / small-w
            sres = strong_invert(u)
            if sres is not None:
                method = "strong-Weyl"
                Flo, Fhi = tc_from_w(float(sres[0]), float(sres[1]))
        else:                   # u >= hlo[-1]: weak / large-w
            w_lo, w_hi = _weak_invert(u)
            if w_lo is not None:
                method = "weak-W"
                Flo, Fhi = tc_from_w(w_lo, w_hi)
            else:               # beyond finite-N weak grid: asymptotic (NOT hard-cert)
                w_lo, w_hi = _weak_asympt_invert(u)
                if w_lo is not None:
                    method = "weak-asympt*"
                    Flo, Fhi = tc_from_w(w_lo, w_hi)
        if method is None or Flo is None:
            print(f"    {lam:>8.3f}   (seam / not covered)")
            continue
        rw = (Fhi - Flo) / Flo
        rows.append((lam, Flo, Fhi, rw, method))
        if method == "crossover":
            cross_max = max(cross_max, rw)
        elif method == "strong-Weyl":
            strong_max = max(strong_max, rw)
        elif method == "weak-W":
            weak_max = max(weak_max, rw)
        print(f"    {lam:>8.3f} {Flo:>13.6e} {Fhi:>13.6e} {rw:>9.2e} {method:>12}")

    print("    (*) weak-asympt*: F_lo (lower T_c) side is PROVEN (Cor iii lower on h);")
    print("        F_hi (upper T_c) side is Theorem W' asymptotic form (not finite-N).")

    # dense worst-case scan across the crossover lambda range (honest max width)
    lam_lo_cross = 1.0 / hlo[-1]     # smallest lambda still crossover-covered
    lam_hi_cross = 1.0 / hhi[0]      # largest lambda still crossover-covered
    dense = np.geomspace(lam_lo_cross * 1.001, lam_hi_cross * 0.999, 4000)
    cross_scan_max = 0.0
    cross_scan_at = None
    for lam in dense:
        r = invert_envelope(ws, hlo, hhi, 1.0 / lam)
        if r is None:
            continue
        Flo, Fhi = tc_from_w(r[0], r[1])
        rw = (Fhi - Flo) / Flo
        if rw > cross_scan_max:
            cross_scan_max = rw
            cross_scan_at = lam
    cross_max = max(cross_max, cross_scan_max)

    print("\n[5] CERTIFIED WORST-CASE RELATIVE WIDTHS  (F_hi - F_lo)/F_lo")
    print(f"    crossover   (lambda in [{lam_lo_cross:.3f},{lam_hi_cross:.2f}]) : "
          f"{cross_max:.3e}  <-- HEADLINE  (dense-scan max at lambda={cross_scan_at:.3f})")
    print(f"    strong-Weyl (lambda in (34,50])    : {strong_max:.3e}")
    print(f"    weak-W      (finite-N, lambda<0.58): {weak_max:.3e}")

    # ---- mnemonic F2c sanity (NOT certified-tight) ----------------------
    print("\n[6] Mnemonic F2c(lambda)=sqrt(0.06401+0.03339 lambda) exp(-1/lambda)")
    print("    (2-param closed form; check it lands INSIDE the certified crossover)")
    inside = 0
    tot = 0
    devs = []
    worst = (0.0, None)
    for (lam, Flo, Fhi, rw, method) in rows:
        if method != "crossover":
            continue
        f2c = np.sqrt(0.06401 + 0.03339 * lam) * np.exp(-1.0 / lam)
        tot += 1
        mid = 0.5 * (Flo + Fhi)
        dev = abs(f2c - mid) / mid
        devs.append(dev)
        if dev > worst[0]:
            worst = (dev, lam)
        if Flo <= f2c <= Fhi:
            inside += 1
    med = float(np.median(devs)) if devs else 0.0
    print(f"    F2c inside certified bracket at {inside}/{tot} crossover lambda")
    print(f"    F2c vs certified midline: median {med*100:.2f}%, "
          f"worst {worst[0]*100:.2f}% (at lambda={worst[1]})  -> mnemonic only")

    print(f"\nTIME: {time.time()-t_start:.1f}s   PASS {PASS[0]}  FAIL {FAIL[0]}")
    return dict(ws=ws, hlo=hlo, hhi=hhi, rows=rows,
                cross_max=cross_max, strong_max=strong_max, weak_max=weak_max)


# --- weak-end inversion helpers (Theorem W finite-N, then asymptotic) --------
_WEAK_WS = None
_WEAK_LO = None
_WEAK_HI = None


def _build_weak_grid(w_max=10.0, n=18, K=6):
    """Precompute a monotone Theorem-W (finite-N) envelope on w in [2.5, w_max].
    Starts at 2.5 so it OVERLAPS the crossover (whose top node is w=3) and closes
    the seam at lambda ~ 0.58.  Matrix size is capped for cost."""
    global _WEAK_WS, _WEAK_LO, _WEAK_HI
    ws = np.geomspace(2.5, w_max, n)
    lo = np.zeros(n)
    hi = np.zeros(n)
    for i, w in enumerate(ws):
        a, b = weak_h_bracket(float(w), K=K)
        lo[i], hi[i] = a, b
    _WEAK_WS, _WEAK_LO, _WEAK_HI = ws, lo, hi


def _weak_invert(u):
    if _WEAK_WS is None:
        _build_weak_grid()
    res = invert_envelope(_WEAK_WS, _WEAK_LO, _WEAK_HI, u)
    if res is None:
        return None, None
    return res


# Weak-coupling deep tail.
#   LOWER on h: PROVEN (P4_main_theorems.md Cor (iii), "first rigorous lower half"):
#       h(w) >= ln w + c0 - 1/(12 w) - 3/(16 w^2),   c0 = EG + ln4 - 3/2.
#   UPPER on h: Theorem W' asymptotic FORM h ~ L + c0 + b1/L + O(1/L^2) with an
#       empirically-calibrated remainder pad (NOT a finite-N certificate).
# => The F_lo (lower T_c) side of the deep-tail bracket is PROVEN; the F_hi side
#    is asymptotic.  Labelled accordingly.
_EG = 0.5772156649015329
_C0 = _EG + np.log(4.0) - 1.5
_B1_HI = 0.360
_REM = 0.75   # upper-side remainder pad; validated numerically for w in [5,100]


def _h_lo_proven(w):
    """Proven lower bound on h(w) (Cor (iii))."""
    return np.log(w) + _C0 - 1.0 / (12 * w) - 3.0 / (16 * w * w)


def _weak_asympt_invert(u):
    """Enclose w* with h(w*)=u.  hlo = PROVEN Cor(iii); hhi = asymptotic (labelled)."""
    from scipy.optimize import brentq

    def hlo(w):
        return _h_lo_proven(w)

    def hhi(w):
        L = np.log(w)
        return L + _C0 + _B1_HI / L + _REM / L**2

    # w* >= w where hhi(w)<=u ; w* <= w where hlo(w)>=u.  Search on a wide bracket;
    # both hlo,hhi are increasing in w, so if the target sign is absent the root is
    # outside [wa,wb] and we clamp to the endpoint (monotone => safe one-sided bound).
    wa, wb = 2.5, 1e13

    def solve(f):
        fa, fb = f(wa), f(wb)
        if fa == 0:
            return wa
        if fb == 0:
            return wb
        if fa * fb > 0:
            return wa if abs(fa) < abs(fb) else wb
        return brentq(f, wa, wb)

    w_lo = solve(lambda w: hhi(w) - u)
    w_hi = solve(lambda w: hlo(w) - u)
    if w_hi < w_lo:
        w_lo, w_hi = w_hi, w_lo
    return w_lo, w_hi


# ===========================================================================
# ENTRY POINT
# ===========================================================================
if __name__ == "__main__":
    if "test" in sys.argv:
        run_selftest()
    else:
        N = 1000
        for a in sys.argv[1:]:
            if a.startswith("N="):
                N = int(a[2:])
        main(N_cross=N)
