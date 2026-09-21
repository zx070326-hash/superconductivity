"""CORRECTED diagnosis of Lemma C's two conditions at mu* > 0.

WHY THE PREVIOUS VERSION OF THIS SCRIPT WAS WRONG
-------------------------------------------------
The repo (P4_concavity_q_route.md) uses, at mu*=0,
    A' = (A - B)/x,   <w_j, A v> = h <w_j,v> = 0  (j != 1),
hence  <w_j, A' v> = -<w_j, B v>/x  and
    C := q' - E = -(2/x) sum_{j>1} <w_j,B v>^2/(h - lam_j)  <= 0,   H'' = -(E+C)/x.
At mu*>0 the top vector obeys  M v = h v  with  M = A - Cm,  so
    <w_j, A v> = <w_j, Cm v>  !=  0,
and the sign-definite closed form above is NO LONGER VALID.  Using it anyway
(as I first did) reports |C| ~ 1e-3 and a false "(B) OK" everywhere.

WHAT IS ACTUALLY TRUE AT ANY mu*
--------------------------------
M(x) = O[f_x] - Cm  with Cm x-independent, so plain 2nd-order perturbation
theory for the top eigenvalue gives, with no structural assumption,
    H''(x) = <v, M'' v> + R,      R := 2 sum_{j!=1} <w_j, M' v>^2/(h-lam_j) >= 0
    M' = O[d_x f] = (A - B)/x,    M'' = -2 O[psi_0],  psi_0 = k^2/(k^2+x)^3
    Psi := <v, O[psi_0] v>   =>   H'' = -2 Psi + R.
In repo units  E = 2 x Psi  and  C = -x R, so
    (A)  E > 0        <=>  Psi > 0        [2nd-order EXTREMALITY of Einstein]
    (B)  |C| < E      <=>  R < 2 Psi      [needed on top of (A) for CONCAVITY]
This form is verified below against an independent 5-point finite difference
of h(x), and the cross-term formula for C is verified against -xR.
Strictness: numerical-evidence (float64, finite basis).
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

N = 300
_n = np.arange(N)
_u = 1.0 / np.sqrt(2.0 * _n + 1.0)
_ks = np.arange(2 * N + 1, dtype=float)
_absd = np.abs(_n[:, None] - _n[None, :])
_sumd = _n[:, None] + _n[None, :] + 1
_uu = _u[:, None] * _u[None, :]


def O(f):
    A = (f[_absd] + f[_sumd]) * _uu
    L1 = np.zeros(N)
    L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * _n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * _n + 1] * _u ** 2)
    return A


def f_ein(x):
    f = x / (_ks ** 2 + x)
    f[0] = 1.0
    return f


def coulomb(mu_t, Nc):
    uc = np.where(_n < (Nc if Nc else N), _u, 0.0)
    return 2.0 * mu_t * np.outer(uc, uc)


print("Lemma C at mu*>0:  (A) Psi>0 [extremality]   (B) R<2Psi [concavity]")
print("naive|C| = repo's mu*=0 closed form, shown to demonstrate it is invalid")
print()
print(f"{'x':>8} {'varpi':>6} {'mu_t':>5} {'Nc':>4} {'gap':>8} {'E=2xPsi':>10} "
      f"{'|C|=xR':>10} {'ratio':>7} {'naive|C|':>10} {'H\"(pert)':>10} "
      f"{'H\"(FD)':>10} {'A':>3} {'B':>4}")
for mu_t, Nc in ((0.0, None), (0.05, 6), (0.13, 6), (0.13, 45), (0.26, 6)):
    Cm = coulomb(mu_t, Nc) if mu_t else np.zeros((N, N))
    for x in (0.01, 0.09, 0.1089, 0.16, 0.25, 1.0, 9.0, 100.0):
        s = _ks ** 2
        f = f_ein(x)
        A_ = O(f)
        B_ = O(f ** 2)
        M = A_ - Cm
        lam, W = eigh(M)
        v = W[:, -1]
        gam = lam[-1] - lam[-2]
        Mp = (A_ - B_) / x                              # dM/dx
        Psi = float(v @ (O(s / (s + x) ** 3) @ v))
        cj = W[:, :-1].T @ (Mp @ v)
        R = float(2.0 * np.sum(cj ** 2 / (lam[-1] - lam[:-1])))
        # repo's mu*=0 closed form, applied blindly (should disagree at mu*>0)
        bj = W[:, :-1].T @ (B_ @ v)
        naive = abs(-(2.0 / x) * float(np.sum(bj ** 2 / (lam[-1] - lam[:-1]))))
        Hpert = -2 * Psi + R
        dx = x / 50.0
        vals = [float(eigh(O(f_ein(max(x + i * dx, 1e-12))) - Cm,
                           eigvals_only=True)[-1]) for i in (-2, -1, 0, 1, 2)]
        Hfd = (-vals[0] + 16 * vals[1] - 30 * vals[2]
               + 16 * vals[3] - vals[4]) / (12 * dx ** 2)
        E, absC = 2 * x * Psi, x * R
        print(f"{x:>8.4f} {math.sqrt(x):>6.3f} {mu_t:>5.2f} {str(Nc):>4} "
              f"{gam:>8.4f} {E:>10.3e} {absC:>10.3e} {absC/E:>7.3f} "
              f"{naive:>10.3e} {Hpert:>10.3e} {Hfd:>10.3e} "
              f"{'ok' if Psi>0 else 'NO':>3} "
              f"{'ok' if R<2*Psi else 'FAIL':>4}")
    print()
