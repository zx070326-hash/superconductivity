"""
How far does a fully EIGENVECTOR-FREE analytic chain reach from x=0?

  F_lo(x)  via the rearrangement identity (general nonneg kernel, g2-branch note
           2026-09-13 Thm 1.1):  for p_n = u_n v_n,
             <v,O[phi]v> = sum p_n^2 T_phi(n) + <p,Hank_phi p> - (1/2) D_phi[p]
           and D_phi[p] <= c D_f[p] with c := max_k phi(k)/f(k), while
           (1/2) D_f[p] = sum p^2 T_f + <p,Hank_f p> - h.  Hence
             F >= c h - sum p_n^2 [c T_f(n) - T_phi(n)] - <p,(c Hank_f - Hank_phi) p>
               >= c h_lo - max_n [c T_f - T_phi](n) * S - lam_max(c Hank_f - Hank_phi) * S,
           S := sum p_n^2 <= 1  (u_n <= 1, |v|=1).                       [analytic]
  E_hi(x)  <= ||O[psi_x]||  (variance <= norm^2), norm by u-weighted Schur.  [analytic]
  g_lo(x)  = 0.4 h_lo (w<=2, trace-bound route) or h_lo - Theta* (w>=2).   [analytic-with-series]
  h_lo     : here the true h (any real proof uses a Rayleigh quotient; loss small).
Criterion: 2 x^2 E_hi^2 / g_lo < 2 F_lo.
float64, numerical-evidence for the *reach* of this chain.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from lemmaC_compact_probe import O_of, top, THETA_STAR

KX = 20000


def tails(phi, K=KX):
    k = np.arange(1, K + 1, dtype=float)
    pk = phi(k)
    T = np.concatenate([[pk.sum()], pk.sum() - np.cumsum(pk)])  # T(n) = sum_{k>n}, n=0..K
    return T


def schur_u(phi, N=2000):
    """u-weighted Schur bound for ||O[phi]|| using |entries|: max_n sum_m u_m^2 |M_nm|."""
    M = np.abs(O_of(phi, N))
    n = np.arange(N)
    w = 1.0 / (2 * n + 1)
    # entries of O are u_n u_m [..] ; weighted row sum sum_m |O_nm| u_m/u_n = sum_m u_m^2 |[..]| + |L1_n|
    R = (M * (w[None, :] / w[:, None]) ** 0.5).sum(axis=1)
    return R.max()


def lam_max_hankel_diff(c, f, phi, N=2000):
    n = np.arange(N)
    S = n[:, None] + n[None, :] + 1
    H = (c * f(S) - phi(S)) * np.outer(1 / np.sqrt(2 * n + 1), 1 / np.sqrt(2 * n + 1))
    return np.linalg.eigvalsh(H)[-1]


if __name__ == "__main__":
    print(f"{'w':>5} {'x':>7} {'h':>7} {'c':>9} {'F_true':>8} {'F_lo':>8} {'capt':>6} {'E_true':>8} {'E_hi':>8} {'g_lo':>7} {'ratio':>8}")
    for w in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 2.0]:
        x = w * w
        f = lambda k: x / (k * k + x)
        phi = lambda k: k * k / (k * k + x) ** 3
        psi = lambda k: 1.0 / (k * k + x) ** 2
        N = 600
        A = O_of(f, N)
        h, v, l2, _, _ = top(A)
        F_true = v @ O_of(phi, N) @ v
        r = O_of(psi, N) @ v; r -= v * (v @ r); E_true = np.linalg.norm(r)
        c = 1.0 / (x * (1 + x) ** 2)          # max_k phi/f attained at k=1
        Tf, Tp = tails(f), tails(phi)
        m1 = np.max(c * Tf - Tp)
        m2 = lam_max_hankel_diff(c, f, phi)
        F_lo = c * h - m1 - m2
        E_hi = schur_u(psi)
        g_lo = 0.4 * h if w <= 2 else h - THETA_STAR
        ratio = x * x * E_hi ** 2 / (g_lo * F_lo) if F_lo > 0 else np.inf
        print(f"{w:5.2f} {x:7.4f} {h:7.4f} {c:9.3f} {F_true:8.4f} {F_lo:8.4f} {F_lo/F_true:6.3f} {E_true:8.4f} {E_hi:8.4f} {g_lo:7.4f} {ratio:8.4f}")
