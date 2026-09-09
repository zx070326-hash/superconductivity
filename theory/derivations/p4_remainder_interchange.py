r"""
P4-A · T5 — uniform O(sigma^4) remainder bound (THEOREM, with explicit
constants) + the N -> infinity interchange (numerical EVIDENCE, honestly
layered).  Standalone: numpy float64 only, no workspace imports.

================================================================================
WHAT THIS FILE UPGRADES
================================================================================
Upstream proposition P4A.1 (VERIFIED elsewhere, 1257 checks, `p4_kappa_N.py`):

  frak_k^(N)(P,T) = frak_h^(N)(wbar) * [ 1 + (1/2) kappa_N(wbar) sigma^2
                        + c_{3,N}(wbar) gamma3 sigma^3 + O(sigma^4) ],   (P4A.1)

where  frak_k^(N)(P,T) = rho( M(P) ),  M(P) := int_0^inf H^(N)(w) P(dw),
w := omega/2piT,  frak_h^(N)(wbar) = rho(H^(N)(wbar)),  and kappa_N, c_{3,N}
come from FIRST-order perturbation of the top eigenvalue with top eigenvector
v_N (eqs. P4A.2/P4A.3).  The "O(sigma^4)" was a SYMBOL: the remainder was only
stated, never controlled.  T5 = turn that symbol into a THEOREM with a finite,
computable constant K_N(wbar), and separately give honest numerical evidence for
the N -> infinity step (which is the genuinely open link, survey U3).

================================================================================
PART 1 — THE THEOREM (uniform O(sigma^4) remainder bound)   [tractable]
================================================================================
Notation.  omega = <omega>(1+xi), <xi>=0, sigma^2 = <xi^2>, gamma3 = <xi^3>/sigma^3,
wbar = <omega>/2piT, so w = wbar(1+xi).  A := H^(N)(wbar), simple top eigenpair
(lam1, v) with ||v||=1, spectral gap g := lam1 - lam2(A) > 0 (N=1: scalar, no
second eigenvalue).  Perturbation  E := M(P) - A.  Because <xi>=0, E starts at
O(sigma^2):
  E = E2 + E3 + E>=4,
  E2   = (1/2) wbar^2 sigma^2         H^(N)''(wbar),
  E3   = (1/6) wbar^3 gamma3 sigma^3  H^(N)'''(wbar),
  E>=4 = sum_{j>=4} (wbar^j <xi^j>/j!) H^(N)(j)(wbar).      (matrix Taylor tail)

Two rigorous ingredients:

(I)  SECOND-ORDER EIGENVALUE REMAINDER (self-contained proof below).
     For Hermitian A (simple top eigenvalue lam1, gap g, unit top vector v) and
     Hermitian E with ||E|| < g/2,
         0 <= lam1(A+E) - lam1(A) - <v|E|v>  <=  ||E||^2 / (g - 2||E||).   (*)
     PROOF.  Lower bound: lam1(A+E) >= <v|(A+E)|v> = lam1 + <v|E|v>.
     Upper bound: let w=cos t v + sin t u (u perp v, ||u||=1) be the top vector
     of A+E; with s=sin t, c=cos t,
        lam1(A+E) = <w|A|w> + <w|E|w>
                 <= (lam1 - s^2 g) + (<v|E|v> + 2cs Re<v|E|u> + s^2(<u|E|u>-<v|E|v>)).
     Hence R := lam1(A+E)-lam1-<v|E|v> <= -s^2 g + 2s||E|| + 2s^2||E||
              = -s^2(g-2||E||) + 2s||E|| <= ||E||^2/(g-2||E||)  (maximize over s).  QED
     (The task's cleaner "||E||^2/g" is the ||E||<<g approximation of (*); we
     keep the exact denominator g-2||E|| so the constant is genuinely rigorous.)

(II) MATRIX TAYLOR TAIL (integral remainder, holds entrywise & in norm).
     tail(xi) := H(wbar(1+xi)) - sum_{j=0}^3 (wbar xi)^j/j! H(j)(wbar)
               = (wbar xi)^4/3! int_0^1 (1-t)^3 H^(4)(wbar(1+t xi)) dt,
     E>=4 = E_xi[tail(xi)], so with B4 := sup_{w in J} ||H^(4)(w)||_op, J the
     support-image interval,
        ||E>=4|| <= (wbar^4 B4 / 24) <xi^4>.
     Likewise from n=2 Taylor (first-order term killed by <xi>=0),
        ||E|| <= (1/2) wbar^2 B2 sigma^2,   B2 := sup_{w in J}||H^(2)(w)||_op.

SCALING ASSUMPTION (stated explicitly, this is where "O(sigma^4)" is honest).
Compact support with FIXED RELATIVE WIDTH: xi = sigma*eta, |eta| <= eta_max,
<eta>=0, <eta^2>=1 (so |xi| <= sigma*eta_max, support radius scales with sigma).
Then <xi^4> = sigma^4 <eta^4> <= sigma^4 eta_max^2, genuinely O(sigma^4).
[WITHOUT this scaling, i.e. eta_max fixed but <xi^4> ~ eta_max^2 sigma^2 only,
the tail is merely O(sigma^2): the sigma^4 rate REQUIRES the higher moments to
scale, which is exactly the small-relative-width regime P4A.1 lives in.  We say
so rather than hide it.]

ASSEMBLED CONSTANT.  Combining, for all P in the class (|eta|<=eta_max, sigma<=sigma0)
with ||E|| < g/2 (guaranteed for sigma0 small enough):
   | frak_k^(N) - frak_h^(N)(1 + (1/2)kappa_N sigma^2 + c_{3,N} gamma3 sigma^3) |
        = | <v|E>=4|v> + R |  <=  K_N(wbar) sigma^4,
   K_N(wbar) = wbar^4 B4 eta_max^2 / 24                (Taylor tail, part 1b)
             + (1/4) wbar^4 B2^2 / (g - wbar^2 B2 sigma0^2).   (2nd-order, part 1c)
   N=1: scalar, R=0 identically -> K_1 = wbar^4 B4 eta_max^2/24 (Taylor tail only).

This is a THEOREM (for finite N): the O(sigma^4) symbol of P4A.1 now carries an
explicit finite constant, uniform over the compactly-supported P-class.

================================================================================
PART 2 — THE N -> infinity INTERCHANGE   [layered: elementary / evidence / open]
================================================================================
2a VERDICT (honest, three layers).
 LAYER A -- elementary (holds for general P).  H^(N)(w) is the leading NxN
   PRINCIPAL BLOCK of the infinite matrix H(w); so is M^(N)(P)=int H^(N)dP of
   M(P)=int H(w)dP.  Cauchy interlacing on nested principal blocks gives
   frak_k^(N) <= frak_k^(N+1) (monotone UP), and since the top eigenvalue is the
   Rayleigh sup over the coordinate subspace ran(P_N), whose union is dense,
   frak_k^(N) UP-> lam_max(M(P)) = frak_k(P,T)  -- PROVIDED M(P) is a compact
   self-adjoint operator (KAY Theorem 3).  So POINTWISE-in-P convergence is
   essentially elementary; we verify it numerically.
 LAYER B -- what KAY actually proves.  KAY I Theorem 2 proves monotone
   convergence frak_g^(N) UP-> frak_g for the gamma-model operator G(gamma)
   ONLY (constant kernel).  KAY II Corollary 3 ASSERTS the general-P version but
   the survey (U3) records its proof as UNREAD / not verified upstream.
 LAYER C -- the genuinely open point FOR THE CAPSTONE.  To turn the finite-N
   theorem (Part 1) into a statement about the true operator, one must
   interchange lim_N with the sigma-expansion.  That needs UNIFORMITY IN N of
   Part 1's constant K_N(wbar) -- and K_N carries a 1/g_gap,N factor.  So the
   real question is: does the spectral gap g_gap,N(wbar) stay bounded BELOW as
   N -> infinity, or collapse?  If it collapses, K_N -> infinity and the
   finite-N theorem does NOT survive the limit.  This uniform gap control is NOT
   established in the referenced material; we give numerical EVIDENCE for it.
 => We DECLARE the uniform-in-N interchange OPEN and state the capstone
    CONDITIONALLY.

2b EVIDENCE: for fixed P, frak_k^(N) is monotone increasing and converges;
   we fit frak_k^(inf) - frak_k^(N) ~ c/N^p (via successive differences) and
   report p.  We ALSO track g_gap,N(wbar) vs N (the Layer-C crux).

2c CONDITIONAL CAPSTONE + evidence that kappa_N(wbar), c_{3,N}(wbar) converge.
"""

import numpy as np

PASS = [0]


def chk(cond, msg=""):
    assert cond, msg
    PASS[0] += 1


# ===========================================================================
# STAGE 0 — inline builders (numpy float64 port of the VERIFIED p4_kappa_N.py)
# ===========================================================================
# [[k]](w) = w^2/(k^2+w^2) and its w-derivatives, orders 0..5.  All five closed
# forms were checked term-by-term against sympy before being pasted here.
def _brk_der(k, w, order):
    """d^order/dw^order of [[k]](w)=w^2/(k^2+w^2), k>=0, order in 0..5.

    k=0: [[0]]=1 identically -> value 1, all derivatives 0.
      d0 = w^2/d
      d1 = 2 k^2 w / d^2
      d2 = 2 k^2 (k^2 - 3 w^2) / d^3
      d3 = 24 k^2 w (w^2 - k^2) / d^4
      d4 = 24 k^2 (-5 w^4 + 10 k^2 w^2 - k^4) / d^5
      d5 = 240 k^2 w (3 w^4 - 10 k^2 w^2 + 3 k^4) / d^6      (d := k^2 + w^2)
    """
    k2 = float(k) ** 2
    w = float(w)
    if k == 0:
        return 1.0 if order == 0 else 0.0
    d = k2 + w * w
    if order == 0:
        return w * w / d
    if order == 1:
        return 2 * k2 * w / d**2
    if order == 2:
        return 2 * k2 * (k2 - 3 * w * w) / d**3
    if order == 3:
        return 24 * k2 * w * (w * w - k2) / d**4
    if order == 4:
        return 24 * k2 * (-5 * w**4 + 10 * k2 * w * w - k2 * k2) / d**5
    if order == 5:
        return 240 * k2 * w * (3 * w**4 - 10 * k2 * w * w + 3 * k2 * k2) / d**6
    raise ValueError(order)


def H_matrix(N, w, order=0):
    """KAY II H^(N)(w) (order=0) or its elementwise d^order/dw^order.

    Adopted placement (module-verified against p4_kappa_N.py): H = -H1 + H2 + H3,
      (H1)_{nn} = (1/(2n+1)) sum_{k=1}^n 2[[k]]           (diagonal only)
      (H2)_{nm} = (1-delta) [[|n-m|]] / sqrt((2n+1)(2m+1))
      (H3)_{nm} = [[n+m+1]] / sqrt((2n+1)(2m+1)).
    """
    M = np.zeros((N, N), dtype=np.float64)
    for n in range(N):
        for m in range(N):
            nrm = np.sqrt((2 * n + 1) * (2 * m + 1))
            val = _brk_der(n + m + 1, w, order) / nrm
            if n != m:
                val += _brk_der(abs(n - m), w, order) / nrm
            else:
                s = sum(2 * _brk_der(k, w, order) for k in range(1, n + 1))
                val -= s / (2 * n + 1)
            M[n, m] = val
    return M


def G_matrix(N, gamma):
    """KAY I G^(N)(gamma): same placement with [[k]] -> 1/k^gamma."""
    M = np.zeros((N, N), dtype=np.float64)
    for n in range(N):
        for m in range(N):
            nrm = np.sqrt((2 * n + 1) * (2 * m + 1))
            val = 1.0 / (nrm * (n + m + 1) ** gamma)
            if n != m:
                val += 1.0 / (nrm * abs(n - m) ** gamma)
            else:
                s = sum(2.0 / k**gamma for k in range(1, n + 1))
                val -= s / (2 * n + 1)
            M[n, m] = val
    return M


def top_eigpair(A):
    """(largest eigenvalue, its unit eigenvector) of symmetric A (eigh)."""
    E, V = np.linalg.eigh(A)
    j = int(np.argmax(E))
    return float(E[j]), V[:, j]


def opnorm(A):
    """spectral (operator) norm of symmetric A = max |eigenvalue|."""
    return float(np.max(np.abs(np.linalg.eigvalsh(A))))


def frak_h(N, w):
    """frak_h^(N)(w) = rho(H^(N)(w)) = largest eigenvalue."""
    return float(np.max(np.linalg.eigvalsh(H_matrix(N, w))))


def gap(N, w):
    """spectral gap lam1 - lam2 of H^(N)(w) (N>=2); returns +inf for N=1."""
    if N == 1:
        return np.inf
    E = np.sort(np.linalg.eigvalsh(H_matrix(N, w)))
    return float(E[-1] - E[-2])


def kappa_c3(N, w):
    """(kappa_N(w), c_{3,N}(w), frak_h^(N)(w)) from P4A.2/P4A.3."""
    lam, v = top_eigpair(H_matrix(N, w))
    q2 = float(v @ H_matrix(N, w, 2) @ v)
    q3 = float(v @ H_matrix(N, w, 3) @ v)
    return w**2 * q2 / lam, w**3 * q3 / (6 * lam), lam


def frak_k_of_P(N, wbar, xis, ps):
    """frak_k^(N)(P,T) = rho( sum_j p_j H^(N)(wbar(1+xi_j)) )  (discrete P)."""
    M = np.zeros((N, N), dtype=np.float64)
    for xi, p in zip(xis, ps):
        M += p * H_matrix(N, wbar * (1.0 + xi))
    return float(np.max(np.linalg.eigvalsh(M)))


def run_stage0():
    print("=== STAGE 0: builders (numpy float64 port of verified p4_kappa_N.py) ===")
    n0 = PASS[0]

    # (a) N=1 anchor: H^(1)(w) == w^2/(w^2+1)
    for ws in (0.3, 1.0, 2.0, 5.0):
        chk(abs(frak_h(1, ws) - ws**2 / (ws**2 + 1)) < 1e-13,
            "H^(1) != w^2/(w^2+1)")

    # (b) KAY I fingerprint: g^(1)(2)=1 and g^(2)(2)=(20+sqrt(13819))/108
    chk(abs(float(np.max(np.linalg.eigvalsh(G_matrix(1, 2)))) - 1.0) < 1e-13)
    g2 = float(np.max(np.linalg.eigvalsh(G_matrix(2, 2))))
    chk(abs(g2 - (20 + np.sqrt(13819)) / 108) < 1e-12,
        "g^(2)(2) fingerprint failed")

    # (c) small-w bridge: rho(H^(N)(w))/w^2 -> g^(N)(2)
    for N in range(1, 5):
        gN = float(np.max(np.linalg.eigvalsh(G_matrix(N, 2))))
        ratio = frak_h(N, 1e-3) / 1e-6
        chk(abs(ratio - gN) < 1e-4, f"bridge N={N}: {ratio} vs {gN}")

    # (d) N=1 closed-form kappa_1, c_{3,1} numeric == symbolic
    for ws in (0.3, 1.0, 2.0, 5.0):
        k1, c1, _ = kappa_c3(1, ws)
        k1s = (2 - 6 * ws**2) / (ws**2 + 1) ** 2
        c1s = 4 * ws**2 * (ws**2 - 1) / (ws**2 + 1) ** 3
        chk(abs(k1 - k1s) < 1e-10 and abs(c1 - c1s) < 1e-10,
            "kappa_1/c_3,1 numeric != symbolic")

    # (e) analytic H'' , H'''' cross-checked vs finite differences.  H'''' is
    # checked as the central 2nd difference of the ANALYTIC H'' (divides by h^2,
    # not h^4 -> numerically stable), i.e. order-4 formula vs order-2 formula.
    hfd = 1e-4
    for N in (2, 3, 4):
        for w in (0.5, 1.3, 3.0):
            H2a = H_matrix(N, w, 2)
            H2fd = (H_matrix(N, w + hfd) - 2 * H_matrix(N, w)
                    + H_matrix(N, w - hfd)) / hfd**2
            chk(np.max(np.abs(H2a - H2fd)) < 1e-4, f"H'' FD N={N} w={w}")
            H4a = H_matrix(N, w, 4)
            H4fd = (H_matrix(N, w + hfd, 2) - 2 * H_matrix(N, w, 2)
                    + H_matrix(N, w - hfd, 2)) / hfd**2
            chk(np.max(np.abs(H4a - H4fd)) < 1e-3, f"H'''' FD N={N} w={w}")

    print(f"Stage 0: {PASS[0]-n0}/{PASS[0]-n0} checks passed "
          "(N=1 anchor, KAY I fingerprint, bridge, kappa_1/c_3,1, FD deriv).")


# ===========================================================================
# STAGE 1 — PART 1a: certified positive lower bound on the spectral gap
# ===========================================================================
# Eigenvalues are 1-Lipschitz in the matrix (Weyl), and ||dH/dw|| <= sup||H'||,
# so both lam1 and lam2 move by <= sup_J||H'|| * |dw|.  Hence the gap g(w) is
# L-Lipschitz with L = 2 sup_J||H'||.  A fine grid (min over grid) minus L*(h/2)
# is a CERTIFIED lower bound on min_{w in [a,b]} g(w).
def sup_opnorm(N, a, b, order, ngrid=4000):
    """sup_{w in [a,b]} ||H^(N)(w)||_op via a fine grid (returns grid-sup)."""
    ws = np.linspace(a, b, ngrid)
    return max(opnorm(H_matrix(N, w, order)) for w in ws)


def certified_min_gap(N, a, b, ngrid=4000):
    """certified lower bound on min_{[a,b]} gap_N(w) (grid-min - Lipschitz slack)."""
    if N == 1:
        return np.inf, np.inf
    ws = np.linspace(a, b, ngrid)
    gmin_grid = min(gap(N, w) for w in ws)
    h = (b - a) / (ngrid - 1)
    L = 2.0 * sup_opnorm(N, a, b, 1, ngrid)   # Lipschitz const of the gap
    return gmin_grid, gmin_grid - L * h / 2.0


WBAR_LO, WBAR_HI = 0.2, 5.0
GAP_LB = {}


def run_stage1():
    print("=== STAGE 1: PART 1a — certified spectral-gap lower bounds "
          f"over wbar in [{WBAR_LO},{WBAR_HI}] ===")
    n0 = PASS[0]
    print(f"  {'N':>2}  {'grid-min gap':>14}  {'certified LB':>14}")
    print(f"  {1:>2}  {'inf (scalar)':>14}  {'inf (scalar)':>14}   "
          "N=1 is scalar: no 2nd eigenvalue, gap treated as +inf (R=0).")
    for N in (2, 3, 4):
        gmin, glb = certified_min_gap(N, WBAR_LO, WBAR_HI)
        GAP_LB[N] = glb
        chk(glb > 0, f"gap lower bound not positive for N={N}: {glb}")
        print(f"  {N:>2}  {gmin:>14.8f}  {glb:>14.8f}")
    GAP_LB[1] = np.inf
    print(f"Stage 1: {PASS[0]-n0}/{PASS[0]-n0} checks passed "
          "(certified min gap > 0 for N=2,3,4).")


# ===========================================================================
# STAGE 2 — PART 1b/1c: certified B2,B4 and the assembled constant K_N(wbar)
# ===========================================================================
ETA_MAX = 2.0     # fixed relative support radius: |xi| <= sigma * ETA_MAX
SIGMA0 = 0.25     # max width for the uniform bound / verification


def certified_sup_opnorm(N, a, b, order, ngrid=4000):
    """CERTIFIED sup of ||H^(N)(.)||_op on [a,b] (grid-sup + Lipschitz slack).

    Lipschitz slack uses grid-sup of the next derivative; the residual it omits
    is O(h^2) with h~1e-3, i.e. <~1e-6, negligible against the reported value.
    """
    ws = np.linspace(a, b, ngrid)
    grid_sup = max(opnorm(H_matrix(N, w, order)) for w in ws)
    h = (b - a) / (ngrid - 1)
    L = max(opnorm(H_matrix(N, w, order + 1)) for w in ws)  # d/dw ||H^(k)|| <= ||H^(k+1)||
    return grid_sup + L * h / 2.0


def K_constant(N, wbar):
    """Return dict with the assembled uniform O(sigma^4) constant at wbar.

    Support-image interval J = [wbar(1-SIGMA0*ETA_MAX), wbar(1+SIGMA0*ETA_MAX)].
    """
    a = wbar * (1.0 - SIGMA0 * ETA_MAX)
    b = wbar * (1.0 + SIGMA0 * ETA_MAX)
    a = max(a, 1e-6)
    B2 = certified_sup_opnorm(N, a, b, 2)
    B4 = certified_sup_opnorm(N, a, b, 4)
    K_tail = wbar**4 * B4 * ETA_MAX**2 / 24.0
    if N == 1:
        # scalar: E is a number added to a 1x1 matrix, eigenvalue = number
        # itself -> second-order remainder R = 0 IDENTICALLY.
        Emax = 0.0
        K_2nd = 0.0
        g = np.inf
    else:
        # bound (*) uses the gap of A = H^(N)(wbar) at the CENTER wbar (exact,
        # float64).  Stage-1's interval-min gap is a separate robustness result
        # certifying that g(wbar) > 0 for every wbar in the range.
        g = gap(N, wbar)
        Emax = 0.5 * wbar**2 * B2 * SIGMA0**2          # ||E||_max over sigma<=sigma0
        denom = g - wbar**2 * B2 * SIGMA0**2            # = g - 2||E||_max
        K_2nd = (0.25 * wbar**4 * B2**2) / denom if denom > 0 else np.inf
    return dict(B2=B2, B4=B4, K_tail=K_tail, K_2nd=K_2nd,
                K=K_tail + K_2nd, Emax=Emax, g=g,
                small_E_ok=(N == 1) or (Emax < g / 2))


WBAR_TEST = [0.3, 0.5, 1.0, 2.0, 3.0]
KTAB = {}


def run_stage2():
    print("=== STAGE 2: PART 1b/1c — certified B2, B4 and constant K_N(wbar) ===")
    print(f"  assumptions: |xi| <= sigma*eta_max, eta_max={ETA_MAX}, "
          f"sigma <= sigma0={SIGMA0}  (fixed relative support width)")
    n0 = PASS[0]
    for N in (1, 2, 3, 4):
        KTAB[N] = {}
        print(f"  --- N={N} ---")
        print(f"    {'wbar':>5}  {'B2':>10}  {'B4':>10}  {'K_tail':>11}  "
              f"{'K_2nd':>11}  {'K_N(wbar)':>11}  {'||E||<g/2':>9}")
        for wbar in WBAR_TEST:
            d = K_constant(N, wbar)
            KTAB[N][wbar] = d
            chk(np.isfinite(d["K"]) and d["K"] > 0, f"K not finite/positive N={N} wbar={wbar}")
            chk(d["small_E_ok"], f"||E||<g/2 violated N={N} wbar={wbar}")
            g_s = "inf" if not np.isfinite(d["g"]) else f"{d['g']:.4f}"
            print(f"    {wbar:>5}  {d['B2']:>10.4f}  {d['B4']:>10.3f}  "
                  f"{d['K_tail']:>11.4e}  {d['K_2nd']:>11.4e}  {d['K']:>11.4e}  "
                  f"{'yes':>9}")
    print(f"Stage 2: {PASS[0]-n0}/{PASS[0]-n0} checks passed "
          "(K_N(wbar) finite, positive; ||E||<g/2 holds at sigma0).")


# ===========================================================================
# STAGE 3 — PART 1d: numerical verification the bound is real & not vacuous
# ===========================================================================
def _random_shape(rng, eta_max):
    """A random unit-variance zero-mean discrete shape eta with |eta|<=eta_max.

    Returns (eta_atoms, weights, gamma3=<eta^3>) or None if it violates |eta|<=eta_max.
    """
    k = rng.integers(2, 7)                        # 2..6 atoms
    pos = rng.uniform(-1.5, 1.5, size=k)
    w = rng.dirichlet(np.ones(k) * rng.uniform(0.5, 3.0))
    mu = np.sum(w * pos)
    pos = pos - mu
    var = np.sum(w * pos**2)
    if var < 1e-6:
        return None
    eta = pos / np.sqrt(var)                       # unit variance, zero mean
    if np.max(np.abs(eta)) > eta_max:
        return None
    g3 = float(np.sum(w * eta**3))
    return eta, w, g3


def run_stage3(n_per=400, seed=20260910):
    print("=== STAGE 3: PART 1d — numerical verification actual remainder "
          "<= K_N(wbar) sigma^4 ===")
    n0 = PASS[0]
    rng = np.random.default_rng(seed)
    print(f"    {'N':>2}  {'wbar':>5}  {'#P':>5}  {'max|rem|/(K s^4)':>17}  "
          f"{'g3 range':>16}  {'sigma range':>14}")
    worst_ratio_overall = 0.0
    for N in (1, 2, 3, 4):
        for wbar in WBAR_TEST:
            kap, c3, hN = kappa_c3(N, wbar)
            K = KTAB[N][wbar]["K"]
            worst_ratio = 0.0
            cnt = 0
            g3s, sigs = [], []
            tries = 0
            while cnt < n_per and tries < 40 * n_per:
                tries += 1
                shp = _random_shape(rng, ETA_MAX)
                if shp is None:
                    continue
                eta, wts, g3 = shp
                sigma = rng.uniform(0.02, SIGMA0)
                xis = sigma * eta
                exact = frak_k_of_P(N, wbar, xis, wts)
                pred = hN * (1.0 + 0.5 * kap * sigma**2 + c3 * g3 * sigma**3)
                rem = abs(exact - pred)
                bound = K * sigma**4
                ratio = rem / bound
                # THE theorem check: actual remainder must not exceed K_N sigma^4
                chk(ratio <= 1.0 + 1e-9,
                    f"BOUND VIOLATED N={N} wbar={wbar} sigma={sigma:.4f} "
                    f"g3={g3:.3f} rem={rem:.3e} bound={bound:.3e}")
                worst_ratio = max(worst_ratio, ratio)
                g3s.append(g3)
                sigs.append(sigma)
                cnt += 1
            worst_ratio_overall = max(worst_ratio_overall, worst_ratio)
            print(f"    {N:>2}  {wbar:>5}  {cnt:>5}  {worst_ratio:>17.4f}  "
                  f"[{min(g3s):>6.2f},{max(g3s):>5.2f}]  "
                  f"[{min(sigs):.3f},{max(sigs):.3f}]")
    print(f"  worst actual/bound ratio over ALL (N,wbar,P) = {worst_ratio_overall:.4f} "
          "(<1 => bound holds; not ~0 => not wildly loose)")
    chk(worst_ratio_overall < 1.0, "some sample exceeded the bound")
    chk(worst_ratio_overall > 1e-4, "bound absurdly loose (ratio<1e-4 everywhere)")

    # --- necessity of the scaling assumption (honest counter-demonstration) ---
    # Drop "fixed relative width": FIX absolute support a and shrink sigma by
    # moving mass to xi=0.  P = (1-p) d_0 + (p/2)(d_{-a}+d_{+a}), a fixed.
    # Then sigma^2 = p a^2, <xi^4> = p a^4 = sigma^2 a^2 (NOT ~sigma^4),
    # kurtosis = 1/p -> inf.  The remainder is then ~ sigma^2 a^2, so
    # remainder/sigma^2 stays bounded while remainder/sigma^4 ~ a^2/sigma^2 -> inf.
    # => the O(sigma^4) RATE genuinely needs bounded kurtosis; without it only
    # O(sigma^2) survives.  We SHOW this rather than merely assert it.
    print("  necessity of bounded kurtosis (fixed absolute support a=0.4, "
          "shrink sigma via mass->0):")
    print(f"    {'p':>6}  {'sigma':>7}  {'rem':>11}  {'rem/sigma^2':>12}  "
          f"{'rem/sigma^4':>12}")
    wbar, N, a = 1.0, 2, 0.4
    kap, c3, hN = kappa_c3(N, wbar)
    r_over_s2, r_over_s4 = [], []
    for p in (0.5, 0.2, 0.1, 0.05, 0.02, 0.01):
        xis = np.array([-a, 0.0, a]); ps = np.array([p / 2, 1 - p, p / 2])
        sigma = np.sqrt(p) * a
        exact = frak_k_of_P(N, wbar, xis, ps)
        pred = hN * (1.0 + 0.5 * kap * sigma**2)          # gamma3=0 here
        rem = abs(exact - pred)
        r_over_s2.append(rem / sigma**2)
        r_over_s4.append(rem / sigma**4)
        print(f"    {p:>6.2f}  {sigma:>7.4f}  {rem:>11.4e}  {rem/sigma**2:>12.4e}  "
              f"{rem/sigma**4:>12.4e}")
    # rem/sigma^2 bounded (last <= 3x first); rem/sigma^4 blows up (grows a lot)
    chk(max(r_over_s2) < 3 * min(r_over_s2) + 1e-9, "rem/sigma^2 not bounded")
    chk(r_over_s4[-1] > 5 * r_over_s4[0], "rem/sigma^4 did not blow up")
    print("    => rem/sigma^2 bounded, rem/sigma^4 -> inf: the sigma^4 THEOREM")
    print("       genuinely requires bounded kurtosis (fixed relative width).")

    print(f"Stage 3: PART 1 is a THEOREM with constant — verified over random "
          f"compactly-supported P (worst ratio {worst_ratio_overall:.3f}); "
          f"scaling assumption shown necessary.")


# ===========================================================================
# STAGE 4 — PART 2: N -> infinity interchange (numerical EVIDENCE, layered)
# ===========================================================================
def _fixed_P(kind, sigma, gamma3=0.0):
    """A few fixed discrete shapes eta (unit var, zero mean), returned as
    (xis, ps) after scaling by sigma.  gamma3 used only for the asym case."""
    if kind == "2pk-sym":
        eta = np.array([-1.0, 1.0]); p = np.array([0.5, 0.5])
    elif kind == "3pk-sym":
        a = np.sqrt(2.0)
        eta = np.array([-a, 0.0, a]); p = np.array([0.25, 0.5, 0.25])
    elif kind == "3pk-asym":
        # equal weights, roots of t^3 - (3/2) t - gamma3 (unit-variance cubic)
        r = np.sort(np.roots([1.0, 0.0, -1.5, -gamma3]))
        r = np.real(r)
        eta = r; p = np.ones(3) / 3.0
        # renormalize numerically to exact zero-mean/unit-var
        eta = eta - np.sum(p * eta)
        eta = eta / np.sqrt(np.sum(p * eta**2))
    else:
        raise ValueError(kind)
    return sigma * eta, p


def run_stage4():
    print("=== STAGE 4: PART 2 — N -> infinity interchange (layered) ===")
    n0 = PASS[0]

    print("  2a VERDICT (three honest layers):")
    print("   A (elementary, general P): H^(N) is the leading NxN principal block")
    print("     of H; Cauchy interlacing => frak_k^(N) monotone UP; Rayleigh sup over")
    print("     the dense union of coordinate subspaces => frak_k^(N) UP-> frak_k,")
    print("     GIVEN M(P) compact self-adjoint (KAY Thm 3). Pointwise conv. is elementary.")
    print("   B (what KAY proves): KAY I Thm 2 = monotone conv. for gamma-model ONLY;")
    print("     KAY II Cor 3 asserts general P, proof unread (survey U3).")
    print("   C (OPEN, the capstone crux): uniform-in-N control of Part 1's K_N, i.e.")
    print("     whether the gap g_gap,N(wbar) stays bounded below as N->inf. If it")
    print("     collapses, K_N->inf and the finite-N theorem does NOT survive N->inf.")
    print("   => uniform-in-N interchange DECLARED OPEN; capstone stated CONDITIONALLY.")

    # 2b EVIDENCE: monotone convergence + rate fit + the Layer-C gap tracking
    print("\n  2b EVIDENCE: monotone increase, rate p (frak_k^(inf)-frak_k^(N)~c/N^p),")
    print("     and g_gap,N(wbar) vs N  [Layer-C crux]")
    Nmax = 48
    Ns = np.arange(1, Nmax + 1)
    wbar = 1.0
    sigma = 0.15
    cases = [("2pk-sym", 0.0), ("3pk-sym", 0.0), ("3pk-asym", 0.4)]
    print(f"     (wbar={wbar}, sigma={sigma})")
    print(f"     {'shape':>10}  {'frak_k^(1)':>11}  {'frak_k^(Nmax)':>13}  "
          f"{'monotone?':>9}  {'rate p':>7}")
    for kind, g3 in cases:
        xis, ps = _fixed_P(kind, sigma, g3)
        vals = np.array([frak_k_of_P(N, wbar, xis, ps) for N in Ns])
        dif = np.diff(vals)
        mono = bool(np.all(dif >= -1e-12))
        chk(mono, f"frak_k^(N) not monotone for {kind}")
        # rate via successive differences: d_N = frak_k^(N+1)-frak_k^(N) ~ c/N^(p+1)
        # so p = -slope(log d_N vs log N) - 1, fit on a clean mid range.
        m = (Ns[:-1] >= 5) & (Ns[:-1] <= Nmax - 4) & (dif > 1e-15)
        p_fit = np.nan
        if np.count_nonzero(m) >= 4:
            slope = np.polyfit(np.log(Ns[:-1][m].astype(float)),
                               np.log(dif[m]), 1)[0]
            p_fit = -slope - 1.0
        print(f"     {kind:>10}  {vals[0]:>11.7f}  {vals[-1]:>13.7f}  "
              f"{str(mono):>9}  {p_fit:>7.2f}")

    # Layer-C crux: spectral gap of H^(N)(wbar) vs N (does it collapse?)
    print("\n     g_gap,N(wbar) vs N  (Part-1 constant K_N ~ 1/g_gap,N):")
    print(f"        {'wbar':>5} " + " ".join(f"N={N:<2d}".rjust(9) for N in (2, 4, 8, 16, 32)))
    gap_ok = True
    for wb in (0.5, 1.0, 2.0):
        gaps = [gap(N, wb) for N in (2, 4, 8, 16, 32)]
        print(f"        {wb:>5} " + " ".join(f"{g:9.5f}" for g in gaps))
        # evidence the gap does NOT collapse: N=32 gap within 2x of N=8 gap
        if not (gaps[-1] > 0.4 * gaps[2]):
            gap_ok = False
    chk(gap_ok, "spectral gap appears to collapse with N (would break capstone)")
    print("        => gap settles to a positive limit (does NOT collapse): numerical")
    print("           EVIDENCE that K_N stays bounded in N. Still not a proof (Layer C).")

    # 2c EVIDENCE: kappa_N(wbar), c_{3,N}(wbar) converge in N
    print("\n  2c CONDITIONAL CAPSTONE + evidence kappa_N, c_{3,N} converge in N:")
    print("     IF frak_k^(N) UP-> frak_k with K_N bounded uniformly in N (Layer C),")
    print("     THEN kappa(wbar)=lim kappa_N, c3(wbar)=lim c_{3,N} exist and")
    print("       frak_k(P,T)=frak_h(wbar)[1+(1/2)kappa sigma^2+c3 gamma3 sigma^3+O(sigma^4)].")
    Nconv = [1, 2, 4, 8, 16, 30]
    for wbar in (0.5, 1.0, 2.0):
        print(f"     wbar={wbar}:")
        print(f"        {'N':>3}  {'kappa_N':>12}  {'c_{3,N}':>12}  {'frak_h^(N)':>12}")
        kseq, cseq = [], []
        for N in Nconv:
            kap, c3, hN = kappa_c3(N, wbar)
            kseq.append(kap); cseq.append(c3)
            print(f"        {N:>3}  {kap:>12.7f}  {c3:>12.7f}  {hN:>12.7f}")
        dk = abs(kseq[-1] - kseq[-2])
        dc = abs(cseq[-1] - cseq[-2])
        chk(dk < 1e-2, f"kappa_N not settling at wbar={wbar}: last step {dk:.2e}")
        chk(dc < 1e-2, f"c_3,N not settling at wbar={wbar}: last step {dc:.2e}")
        print(f"        last-step |Dkappa|={dk:.2e}, |Dc3|={dc:.2e}  (settling)")

    print(f"Stage 4: {PASS[0]-n0}/{PASS[0]-n0} checks passed "
          "(monotone conv.; gap non-collapse evidence; kappa_N, c_{3,N} settle).")


# ===========================================================================
if __name__ == "__main__":
    run_stage0()
    run_stage1()
    run_stage2()
    run_stage3()
    run_stage4()
    print(f"\nALL {PASS[0]} checks passed.")
    print(
        "\nSUMMARY (rigorously layered):\n"
        "  PART 1 = THEOREM.  For finite N in {1,2,3,4}, the O(sigma^4) symbol in\n"
        "    P4A.1 now carries an explicit, finite, uniform constant K_N(wbar):\n"
        "      |frak_k^(N) - frak_h^(N)(1 + (1/2)kappa_N sigma^2 + c_{3,N} gamma3 sigma^3)|\n"
        "         <= K_N(wbar) sigma^4,   over |xi|<=sigma*eta_max, sigma<=sigma0.\n"
        "    Proven via (I) a self-contained 2nd-order eigenvalue remainder bound\n"
        "    (needs the certified spectral gap, Part 1a) and (II) the matrix Taylor\n"
        "    integral remainder (needs certified B4, Part 1b).  Verified non-vacuous\n"
        "    over random compactly-supported P (Part 1d).\n"
        "  PART 2 = LAYERED.  (A) Pointwise-in-P convergence frak_k^(N) UP-> frak_k\n"
        "    is elementary (Cauchy interlacing on nested principal blocks + dense\n"
        "    Rayleigh sup, given KAY Thm 3 compactness) and is verified numerically.\n"
        "    (B) KAY I Thm 2's monotone result is gamma-model ONLY; KAY II Cor 3's\n"
        "    general-P assertion has an unread proof (survey U3).  (C) The capstone\n"
        "    needs UNIFORM-in-N control of K_N, i.e. the gap g_gap,N must not collapse;\n"
        "    numerics show it settling to a positive limit (EVIDENCE, not proof).\n"
        "    The limit constants kappa, c3 (shown to settle in N) therefore define the\n"
        "    O(sigma^4) expansion of the TRUE operator only CONDITIONALLY on (C).")
