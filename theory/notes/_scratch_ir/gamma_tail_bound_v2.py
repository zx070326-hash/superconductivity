"""
gamma_tail_bound_v2.py — SOUND, direction-audited tail bounds for the γ-model.

New file (does NOT modify gamma_tail_bound.py), per repo convention.

WHY A V2
--------
Independent audit (team-lead, 2026-09-16) of gamma_tail_bound.py found three
one-sided defects in the γ>1/2 branch (in addition to the γ<=1/2 divergent
branch already patched to +inf):

  (D1) weighted_tail_bound_hi used  Σ_{k>=J} k^{-α} <= ∫_J^∞ x^{-α} dx.
       WRONG direction: for decreasing k^{-α} the left Riemann sum dominates
       the integral, so ∫_J^∞ is a LOWER estimate.  Sound:  J^{-α} + ∫_J^∞.
       (J=1: the whole k=1 term, i.e. ~50% of the sum, was dropped; that J=1
       call is used inside HBB_bound_hi.)
  (D2) Rsup_hi's Hankel term called weighted_tail_bound_hi(N+2+n), i.e. it
       replaced Σ_{m>N} u_m (n+m+1)^{-γ} by Σ_{k>=n+N+2} u_k k^{-γ}.  That
       substitution k=n+m+1 carries u_{k-n-1} -> u_k, and u is decreasing,
       so u_k <= u_{k-n-1}: an UNDER-estimate (factor up to √2 at n≈N).
  (D3) HBB_bound_hi took only the row n=N+1 as sup_{n>N} of the row sums.
       For n deep in the tail the Toeplitz sum is TWO-sided (m<n as well as
       m>n); row N+1 has no left part at all, so it is not the sup.

Consequence: every ||E||, sinΘ, c₁² previously printed as `certificate`
(Phase 1 and the Phase-2 multi-γ scan) rests on a derivation that is not a
chain of inequalities in the claimed direction.  This file (i) gives a sound
derivation and (ii) brute-force checks every sup against float64 row sums on
a long tail truncation.

METHOD: weighted Schur test (p = q = u)
---------------------------------------
For A with nonnegative entries and positive weights p,q:
    ||A|| <= sqrt( sup_n (1/p_n) Σ_m A_nm q_m  ·  sup_m (1/q_m) Σ_n A_nm p_n ).
Take p = q = u, u_m = (2m+1)^{-1/2}.  Off the diagonal G_nm = u_n u_m g(n,m),
so the weighted row sums collapse to  Σ_m g(n,m)/(2m+1)  — decaying like
m^{-1-γ}, convergent for EVERY γ>0 and far tighter than the unweighted sums
(which decay only like m^{-γ-1/2} and need γ>1/2).  This is both sound and
strictly better than v1.

Notation
    H̄_γ(n) := 1 + ∫_1^n x^{-γ}dx          >= Σ_{k=1}^n k^{-γ}
    T̄_β(J) := J^{-β} + J^{1-β}/(β-1)      >= Σ_{k>=J} k^{-β}     (β>1)

Strictness: certificate (mpmath.iv), except brute_force_check() which is
numerical-evidence and only audits the DIRECTION of the certified bounds.
"""
import sys, os, time
import numpy as np
from mpmath import iv, mp

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)


def hi(x):
    return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)


def _setdps(dps):
    mp.dps = dps
    iv.dps = dps


# ---------------------------------------------------------------- primitives
def tail_hi(J, beta):
    """UPPER bound on Σ_{k>=J} k^{-β}, β>1, integer J>=1 (iv in, iv out).

    Σ_{k>=J} k^{-β} = J^{-β} + Σ_{k>=J+1} k^{-β}
                   <= J^{-β} + ∫_J^∞ x^{-β} dx          [k^{-β} <= ∫_{k-1}^k]
                    = J^{-β} + J^{1-β}/(β-1).
    """
    J = iv.mpf(J)
    return J**(-beta) + J**(1 - beta) / (beta - 1)


def partial_hi(n, gamma):
    """UPPER bound on Σ_{k=1}^n k^{-γ}, γ>0, integer n>=1.
    = 1 + Σ_{k=2}^n k^{-γ} <= 1 + ∫_1^n x^{-γ} dx."""
    n = iv.mpf(n)
    g = gamma
    if abs(float(g.mid) - 1.0) < 1e-12:
        return 1 + iv.log(n)
    return 1 + (n**(1 - g) - 1) / (1 - g)


def weighted_tail_hi(J, gamma):
    """UPPER bound on Σ_{k>=J} u_k k^{-γ} with u_k=(2k+1)^{-1/2}, γ>1/2.
    u_k <= (2k)^{-1/2} => <= 2^{-1/2} Σ_{k>=J} k^{-α}, α=γ+1/2 > 1.  [D1 fixed]"""
    alpha = gamma + iv.mpf(1) / 2
    if float(alpha.mid) <= 1.0:
        return mp.inf
    return tail_hi(J, alpha) / iv.sqrt(2)


# ---------------------------------------------- unweighted Schur (v1 repaired)
def C_bound_unweighted_hi(N, gamma_val, dps=30):
    """||C|| <= sqrt(R_sup K_sup), the v1 route with D1+D2 repaired.
    R_n = u_n Σ_{m>=N+1} u_m [(m-n)^{-γ} + (n+m+1)^{-γ}]
      Toeplitz: k=m-n>=N+1-n, u_m = u_{k+n} <= u_k  => <= W(N+1-n)
      Hankel  : (n+m+1)^{-γ} <= m^{-γ}              => <= W(N+1)   [D2 fix]
    K_m non-increasing in m => exact iv sum at m=N+1.
    Kept only for comparison; the weighted route below is much tighter."""
    _setdps(dps)
    g = iv.mpf(gamma_val)
    u = lambda n: 1 / iv.sqrt(iv.mpf(2 * n + 1))
    W_N1 = weighted_tail_hi(N + 1, g)
    if W_N1 == mp.inf:
        return mp.inf, mp.inf, mp.inf
    R = mp.mpf(0)
    for n in range(N + 1):
        R = max(R, hi(u(n) * (weighted_tail_hi(N + 1 - n, g) + W_N1)))
    m = N + 1
    K = iv.mpf(0)
    for n in range(N + 1):
        K += u(n) * (iv.mpf(m - n)**(-g) + iv.mpf(n + m + 1)**(-g))
    K = hi(u(m) * K)
    return mp.sqrt(R * K), R, K


# --------------------------------------------------- weighted Schur for ||C||
def C_bound_weighted_hi(N, gamma_val, dps=30):
    """||C|| <= sqrt(A_sup * B_sup), weights p=q=u.

    A_n = Σ_{m>N} g(n,m)/(2m+1),  n<=N:
      Hankel : (n+m+1)^{-γ}/(2m+1) <= m^{-γ}/(2m)         [n>=0, 2m+1>2m]
               => Σ <= T̄_{1+γ}(N+1)/2                     (n-independent)
      Toeplitz: j=m-n>=max(1,N+1-n); enlarge range to j>=1, m=n+j>=N+1:
               j<=N+1 : 1/(2m+1) <= 1/(2N+3), Σ_{j<=N+1} j^{-γ} <= H̄_γ(N+1)
               j> N+1 : 1/(2m+1) <= 1/(2j)   , Σ_{j>N+1} j^{-1-γ}/2
                                              <= T̄_{1+γ}(N+2)/2
      => A_sup <= H̄_γ(N+1)/(2N+3) + T̄_{1+γ}(N+2)/2 + T̄_{1+γ}(N+1)/2

    B_m = Σ_{n<=N} g(n,m)/(2n+1),  m>N: every term non-increasing in m,
      so B_sup = B_{N+1}, an exact finite iv sum.
    """
    _setdps(dps)
    g = iv.mpf(gamma_val)
    N1 = iv.mpf(N + 1)
    A = (partial_hi(N + 1, g) / (2 * N1 + 1)
         + tail_hi(N + 2, 1 + g) / 2
         + tail_hi(N + 1, 1 + g) / 2)
    m = N + 1
    B = iv.mpf(0)
    for n in range(N + 1):
        B += (iv.mpf(m - n)**(-g) + iv.mpf(n + m + 1)**(-g)) / iv.mpf(2 * n + 1)
    A, B = hi(A), hi(B)
    return mp.sqrt(A * B), A, B


# ----------------------------------------------- weighted Schur for ||G^BB||
def HBB_bound_weighted_hi(N, gamma_val, dps=30):
    """||G^BB|| <= sup_{n>N} Σ_{m>N} |G_nm| u_m/u_n   (symmetric => p=q=u).

    Row n>N splits into diagonal + Hankel + Toeplitz-right + Toeplitz-left.
    Each piece is bounded by a quantity that is NON-INCREASING in n FOR n>=3
    (or is n-independent), so evaluating at n=N+1 IS a sup whenever N>=2.
    [D3 fixed: the left Toeplitz part, absent from row N+1, is included for
    all n.]  Caveat on the n>=3 restriction: pieces (a) and (c) are decreasing
    for all n>=1, but d1(n)=H̄_γ(n)/(n+1) RISES from n=1 to n=2 for γ<2.  The
    decreasing condition for d1 is n^{-γ}(n+1) <= H̄_γ(n); for γ∈(1/2,1) this
    reads n^{1-γ}[1/(1-γ) - 1 - 1/n] >= γ/(1-γ), true for n>=3 since
    1/(1-γ)>2; for γ∈[1,2] it is immediate for n>=2.  We always call with
    n = N+1 >= 101.

    (a) diagonal    |G_nn| <= 2H̄_γ(n)/(2n+1) + (2n+1)^{-1-γ}
        decreasing: d/dn of H̄(n)/(2n+1) <= 0  <=>  n^{-γ}(2n+1) <= 2H̄_γ(n),
        true for all n>=1, γ∈(1/2,2]  (γ<1: 2n^{1-γ} vs 2n^{1-γ}/(1-γ);
        γ=1: 2 vs 2(1+ln n); γ>1: ->0 vs ->2(1+1/(γ-1))).
    (b) Hankel      Σ_{m>N} (n+m+1)^{-γ}/(2m+1) <= T̄_{1+γ}(N+1)/2   (n-indep)
    (c) Toep right  Σ_{j>=1} j^{-γ}/(2(n+j)+1)
                    <= H̄_γ(n)/(2n+1) + T̄_{1+γ}(n+1)/2   (both decreasing)
    (d) Toep left   Σ_{j=1}^{n-N-1} j^{-γ}/(2(n-j)+1), split at j = n/2:
        j<=n/2: 2(n-j)+1 >= n+1 => <= H̄_γ(n)/(n+1)       (decreasing in n)
        j> n/2: j^{-γ} <= (n/2)^{-γ}; m=n-j runs in [N+1, n/2], and
                Σ_{m=N+1}^{n/2} 1/(2m+1) <= 1/(2N+3) + (1/2)ln((n+1)/(2N+3)).
                Nonzero only for n >= 2N+3.  Put x=(n+1)/(2N+3) >= 1; then
                n >= (x(2N+3))/2 so (2/n)^γ <= (4/(2N+3))^γ x^{-γ}, and
                sup_{x>=1} x^{-γ}(a + (1/2)ln x) = a + 1/(2γe)  (a=1/(2N+3)).
                => <= (4/(2N+3))^γ [1/(2N+3) + 1/(2γe)]
    """
    _setdps(dps)
    g = iv.mpf(gamma_val)
    n = N + 1
    n_iv = iv.mpf(n)
    N3 = iv.mpf(2 * N + 3)
    Hn = partial_hi(n, g)
    a = 2 * Hn / (2 * n_iv + 1) + (2 * n_iv + 1)**(-1 - g)
    b = tail_hi(N + 1, 1 + g) / 2
    c = Hn / (2 * n_iv + 1) + tail_hi(n + 1, 1 + g) / 2
    d1 = Hn / (n_iv + 1)
    d2 = (4 / N3)**g * (1 / N3 + 1 / (2 * g * iv.exp(1)))
    tot = a + b + c + d1 + d2
    return hi(tot), dict(diag=hi(a), hankel=hi(b), toep_right=hi(c),
                         toep_left_near=hi(d1), toep_left_far=hi(d2))


def lamB_bound_hi(N, gamma_val, dps=30):
    """UPPER bound on  λ_max(G^BB)  — NOT on ||G^BB||.  (certificate)

    Weighted Gershgorin (D=diag(u); D^{-1}G^BB D has the same spectrum, discs
    centred at G_nn with radius rad(n) = (1/u_n) Σ_{m≠n} |G_nm| u_m):
        λ_max(G^BB) <= sup_{n>N} [ G_nn + rad(n) ]
                    <= sup_{n>N} G_nn  +  sup_{n>N} rad(n)
                    <= 0 + rad_hi,
    because G_nn = -𝔊₁_nn + (2n+1)^{-1-γ} < 0 for every n>=1:  𝔊₁_nn >= 2/(2n+1)
    (the k=1 term alone) and 2/(2n+1) > (2n+1)^{-1-γ}  <=>  2 > (2n+1)^{-γ}. ✓

    rad_hi = (b) + (c) + (d1) + (d2) from HBB_bound_weighted_hi, i.e. that
    bound MINUS its diagonal piece (a).

    WHY THIS MATTERS: G^BB is dominated by its NEGATIVE diagonal -𝔊₁, so
    ||G^BB|| is huge compared with λ_max(G^BB) (factor ~80 at γ=1, N=400).
    Davis-Kahan against the block-DIAGONAL reference head ⊕ G^BB only needs
    λ_max(G^BB) (it enters the separation δ) plus ||C|| (the whole coupling),
    never ||G^BB||.
    """
    H, parts = HBB_bound_weighted_hi(N, gamma_val, dps)
    rad = parts['hankel'] + parts['toep_right'] + parts['toep_left_near'] + parts['toep_left_far']
    return rad, parts


def E_bound_hi(N, gamma_val, dps=30):
    """||E|| <= ||C|| + ||G^BB||, both weighted-Schur certified (certificate)."""
    C, A, B = C_bound_weighted_hi(N, gamma_val, dps)
    H, parts = HBB_bound_weighted_hi(N, gamma_val, dps)
    return C + H, C, H, parts, A, B


# ------------------------------------------------------- brute-force checks
def brute_force_check(N, gamma_val, M, verbose=True):
    """float64 tail truncation at M >> N.  Direction audit only.

    Checks that the certified sups A_sup, B_sup, H_sup really dominate the
    TRUE weighted row/column sums of every row, and reports how loose the
    resulting ||C||, ||G^BB|| bounds are against the true spectral norms.
    """
    gam = gamma_val
    idx = np.arange(M + 1)
    u = 1.0 / np.sqrt(2 * idx + 1)
    ks = np.arange(1, M + 1, dtype=float)
    g1 = np.concatenate([[0.0], np.cumsum(2.0 / ks**gam)]) / (2 * idx + 1.0)

    # ---- C block rows (n<=N): A_n = Σ_{m>N} g(n,m)/(2m+1)
    mtail = idx[N + 1:]
    w = 1.0 / (2 * mtail + 1.0)
    A_rows = np.empty(N + 1)
    for n in range(N + 1):
        A_rows[n] = (((mtail - n)**(-gam) + (n + mtail + 1.0)**(-gam)) * w).sum()
    # ---- C block cols (m>N): B_m = Σ_{n<=N} g(n,m)/(2n+1)
    head = idx[:N + 1]
    wh = 1.0 / (2 * head + 1.0)
    B_cols = np.empty(M - N)
    for i, m in enumerate(mtail):
        B_cols[i] = (((m - head)**(-gam) + (head + m + 1.0)**(-gam)) * wh).sum()
    # ---- BB rows (n>N): |G_nn| + Σ_{m>N,m!=n} g(n,m)/(2m+1)
    BB_rows = np.empty(M - N)
    for i, n in enumerate(mtail):
        d = np.abs(mtail - n).astype(float)
        with np.errstate(divide='ignore'):
            toe = np.where(d > 0, d**(-gam), 0.0)
        s = ((toe + (n + mtail + 1.0)**(-gam)) * w).sum()
        s -= (2 * n + 1.0)**(-gam) * w[i]          # remove m=n Hankel self-term
        BB_rows[i] = s + abs(-g1[n] + (2 * n + 1.0)**(-1.0 - gam))

    # ---- true spectral norms (on the truncation)
    C = (np.outer(u[:N + 1], u[N + 1:])
         * ((mtail[None, :] - head[:, None])**(-gam)
            + (head[:, None] + mtail[None, :] + 1.0)**(-gam)))
    C_true = np.linalg.norm(C, 2)
    nb = min(M - N, 3000)
    sub = mtail[:nb]
    d = np.abs(sub[:, None] - sub[None, :]).astype(float)
    with np.errstate(divide='ignore'):
        toe = np.where(d > 0, d**(-gam), 0.0)
    BB = np.outer(u[sub], u[sub]) * (toe + (sub[:, None] + sub[None, :] + 1.0)**(-gam))
    BB[np.arange(nb), np.arange(nb)] = -g1[sub] + (2.0 * sub + 1.0)**(-gam) * u[sub]**2
    BB_true = float(np.abs(np.linalg.eigvalsh(BB)).max())

    Ehi, Chi, Hhi, parts, A, B = E_bound_hi(N, gam)
    okA = A_rows.max() <= float(A)
    okB = B_cols.max() <= float(B)
    okH = BB_rows.max() <= float(Hhi)
    if verbose:
        print(f"[direction audit] gamma={gam}  N={N}  M={M}")
        print(f"  A_n   true max {A_rows.max():.8f}  cert {float(A):.8f}  "
              f"{'OK' if okA else '*** VIOLATION ***'}  ({float(A)/A_rows.max():.2f}x)")
        print(f"  B_m   true max {B_cols.max():.8f}  cert {float(B):.8f}  "
              f"{'OK' if okB else '*** VIOLATION ***'}  ({float(B)/B_cols.max():.2f}x)")
        print(f"  BBrow true max {BB_rows.max():.8f} @n={mtail[BB_rows.argmax()]}  "
              f"cert {float(Hhi):.8f}  {'OK' if okH else '*** VIOLATION ***'}  "
              f"({float(Hhi)/BB_rows.max():.2f}x)")
        print(f"  ||C||2   true {C_true:.8f}  cert {float(Chi):.8f}  ({float(Chi)/C_true:.2f}x)")
        print(f"  ||G^BB|| true {BB_true:.8f} (sub {nb})  cert {float(Hhi):.8f}  "
              f"({float(Hhi)/BB_true:.2f}x)")
        print(f"  ||E|| cert {float(Ehi):.8f}")
        print("  HBB parts: " + ", ".join(f"{k}={float(v):.6f}" for k, v in parts.items()))
    return dict(ok=okA and okB and okH, C_true=C_true, BB_true=BB_true,
                E_hi=Ehi, C_hi=Chi, H_hi=Hhi,
                A_true=A_rows.max(), B_true=B_cols.max(), BB_row_true=BB_rows.max())


if __name__ == '__main__':
    gam = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    M = int(sys.argv[3]) if len(sys.argv) > 3 else N + 20 * N
    t0 = time.time()
    r = brute_force_check(N, gam, M)
    # v1 comparison
    try:
        from gamma_tail_bound import build_G_matrix, schur_C_bound_hi, HBB_bound_hi
        data = build_G_matrix(N, gam, dps=30)
        C_old = schur_C_bound_hi(data, N, dps=30)
        H_old = HBB_bound_hi(data, N, dps=30)
        print(f"\n  v1 (DEFECTIVE) ||C||<={float(C_old):.8f}  ||G^BB||<={float(H_old):.8f}"
              f"  ||E||<={float(C_old + H_old):.8f}")
    except Exception as e:
        print(f"  [v1 comparison skipped: {e}]")
    Cu, R, K = C_bound_unweighted_hi(N, gam)
    print(f"  v2 unweighted-repaired ||C||<={float(Cu):.8f}")
    print(f"  v2 weighted            ||C||<={float(r['C_hi']):.8f}"
          f"  ||G^BB||<={float(r['H_hi']):.8f}  ||E||<={float(r['E_hi']):.8f}")
    print(f"  [{time.time() - t0:.1f}s]")
