"""Phi_E(lam, mu*, r) := Tc_E(lam,mu*,r)/Tc_E(lam,0)  --  SHARED float64 engine.

r := omega_c / omega_E.  Repo pencil convention (P4D_coulomb_operator.md sec.1):
    M(wbar, mutilde, Nc) = H_mat(., wbar) - mutilde * 2|u_c><u_c|,
    u_c[n] = (2n+1)^{-1/2} 1_{n<Nc},   Nc = floor(r*wbar),   mutilde = mu*/lam,
    T_c boundary: lam * lam_1(M) = 1,  T_c/omega_E = 1/(2 pi wbar*),
    T_c := sup{T : lam lam_1(M(T)) >= 1}  ==  SMALLEST root wbar*.

WINDOW-AWARE root finder (this is load-bearing): the staircase Nc=floor(r wbar)
makes g_mutilde(wbar) piecewise-continuous with DOWN-JUMPS at wbar=k/r.  A naive
uniform grid with step > 1/r can step over a whole window and bisect across a jump.
Here every window k is scanned on its OWN sub-grid and bisection never crosses a
window edge.

LEVEL: numerical-evidence (float64, finite N).
"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DERIV = os.path.normpath(os.path.join(HERE, "..", "..", "derivations"))
if DERIV not in sys.path:
    sys.path.insert(0, DERIV)

from p4_interchange_uniformity import H_mat  # noqa: E402

TWO_PI = 2.0 * np.pi
C_INF = 0.1827262477890
MEV_PER_K = 8.617333262e-2

_HCACHE = {}


def H_cached(N, w):
    key = (int(N), round(float(w), 12))
    A = _HCACHE.get(key)
    if A is None:
        A = H_mat(N, w)
        if len(_HCACHE) > 4000:
            _HCACHE.clear()
        _HCACHE[key] = A
    return A


def K_c_np(N, Nc):
    u = np.zeros(N)
    k = max(0, min(int(Nc), N))
    u[:k] = 1.0 / np.sqrt(2.0 * np.arange(k) + 1.0)
    return 2.0 * np.outer(u, u)


def nc_of(w, r):
    return int(np.floor(r * w + 1e-9))


def g_float(w, mu_t, r, N=220):
    Nc = nc_of(w, r)
    N_eff = max(N, Nc + 40)
    A = H_cached(N_eff, w) - mu_t * K_c_np(N_eff, Nc)
    return float(np.linalg.eigvalsh(A)[-1])


def g_fixedNc(w, mu_t, Nc, N=220):
    """g at FIXED Nc (no staircase) -- used for the within-window bisection."""
    N_eff = max(N, Nc + 40)
    A = H_cached(N_eff, w) - mu_t * K_c_np(N_eff, Nc)
    return float(np.linalg.eigvalsh(A)[-1])


def wstar_sup(lam, mu_star, r, N=220, w_max=40.0, sub=12, eps=1e-9,
              w_free=None, nbis=70):
    """SMALLEST root wbar* of g_mutilde(wbar)=1/lam (=> largest Tc = the sup rule).

    Returns (wstar, Nc_at_root, n_windows_scanned) or (None, None, k) if no root
    below w_max.  Window-by-window: within window k, Nc is CONSTANT, so g is the
    continuous function g_fixedNc(., k) and bisection is safe.

    w_free : optional PROVEN lower cut on the search.  Because K_c >= 0 we have
        M <= K_ph, hence g_mutilde(w) <= h(w) with h strictly increasing; so if
        h(w_free) = gstar (the mu*=0 root) then NO root of g_mutilde can lie below
        w_free and every window k < floor(r*w_free) is root-free.  Skipping those
        windows is therefore sound, not a heuristic.  (Direction audit: skipping a
        window that DID contain a root would return a LARGER wbar*, i.e. a SMALLER
        Tc -- the wrong direction for a ceiling -- which is exactly why the cut must
        be the proven free bound and nothing looser.)
    """
    mu_t = 0.0 if mu_star == 0.0 else mu_star / lam
    gstar = 1.0 / lam
    kmax = int(np.ceil(r * w_max))
    w_lo_first = 0.02
    k0 = 0 if w_free is None else max(0, int(np.floor(r * w_free)) - 1)
    for k in range(k0, kmax + 1):
        a = max(k / r, w_lo_first)
        b = (k + 1) / r
        if b <= w_lo_first:
            continue
        if a >= w_max:
            break
        b = min(b, w_max)
        ws = np.linspace(a + eps, b - eps, sub)
        gs = [g_fixedNc(w, mu_t, k, N) for w in ws]
        if gs[0] >= gstar:
            # Already >= gstar at the LEFT edge of this window: the up-crossing sits
            # exactly ON the window boundary.  This is NOT confined to k=0: the
            # boundary jump of g is DOWN only while <u_c,v>^2 grows with Nc, and
            # <u_c,v>^2 is NON-monotone in Nc (the top eigenvector of M changes
            # sign), so boundary jumps turn UP at large Nc.  Measured instance
            # (N-converged to 10 digits, N=220/400/800):
            #   lam=1, mu*=0.13, r=10: g(Nc=21, w=2.2-) = 0.9983828740
            #                          g(Nc=22, w=2.2 ) = 1.0001528407  (jump +1.77e-3)
            # so gstar=1 is crossed AT wbar=2.2.  Any root finder that looks only for
            # INTERIOR up-crossings misses this root (see the note for the knock-on
            # effect on p4_coulomb_global_tc_certified.invert_sup's nroots/status).
            return float(ws[0]), k, k
        for i in range(1, sub):
            if gs[i - 1] < gstar <= gs[i]:
                lo, hi = ws[i - 1], ws[i]
                for _ in range(nbis):
                    m = 0.5 * (lo + hi)
                    if g_fixedNc(m, mu_t, k, N) < gstar:
                        lo = m
                    else:
                        hi = m
                return float(0.5 * (lo + hi)), k, k
    return None, None, kmax


def tc_over_wE(lam, mu_star, r, N=220, **kw):
    w, Nc, _ = wstar_sup(lam, mu_star, r, N=N, **kw)
    if w is None:
        return None, None, None
    return 1.0 / (TWO_PI * w), w, Nc


def phi_E(lam, mu_star, r, N=220, **kw):
    t0, w0, _ = tc_over_wE(lam, 0.0, 1.0, N=N, **kw)   # mu*=0: r irrelevant
    t1, w1, nc = tc_over_wE(lam, mu_star, r, N=N, **kw)
    if t0 is None or t1 is None:
        return None, None, None, None
    return t1 / t0, t0, t1, nc


def f_E_from_solver(lam, N=220):
    t0, _, _ = tc_over_wE(lam, 0.0, 1.0, N=N)
    return t0 / (C_INF * np.sqrt(lam))


if __name__ == "__main__":
    # smoke test: reproduce the certified mu*=0 anchors of P4_global_tc_certified.md
    print("mu*=0 anchors vs certified table (P4_global_tc_certified.md sec.4):")
    ref = {0.7: (0.07202285, 0.07213087), 1.0: (0.1146122, 0.1147841),
           2.0: (0.2115556, 0.2118729), 5.0: (0.3795354, 0.3801047)}
    for lam, (lo, hi) in ref.items():
        t, w, _ = tc_over_wE(lam, 0.0, 1.0)
        ok = lo <= t <= hi
        print(f"  lam={lam:4.1f}  Tc/wE={t:.8f}  cert[{lo:.7f},{hi:.7f}]  in={ok}")
    print()
    print("mu*>0 anchors vs certified table (P4_coulomb_global_tc_certified.md sec.3, r=3):")
    ref2 = {(1.0, 0.10): (0.080230, 0.080290), (2.0, 0.10): (0.178667, 0.178734),
            (3.0, 0.10): (0.253132, 0.253264), (1.0, 0.13): (0.071121, 0.071189),
            (2.0, 0.13): (0.170535, 0.170595), (1.5, 0.13): (0.126174, 0.126240)}
    for (lam, mus), (lo, hi) in ref2.items():
        t, w, nc = tc_over_wE(lam, mus, 3.0)
        ok = lo <= t <= hi
        print(f"  lam={lam:4.1f} mu*={mus:.2f}  Tc/wE={t:.8f} Nc={nc}  "
              f"cert[{lo:.6f},{hi:.6f}]  in={ok}")
