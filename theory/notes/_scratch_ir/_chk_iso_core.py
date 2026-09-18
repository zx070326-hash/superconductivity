"""Independent KAY-operator / Tc machinery, written from P4_main_theorems.md sec.0 only.

O[f] := -L1[f] + L2[f] + L3[f]   on basis index n = 0..N-1, u_n = (2n+1)^{-1/2}
  L1 diagonal   : (2n+1)^{-1} sum_{k=1}^{n} 2 f(k)      (zero for n=0)
  L2 OFF-diagonal: f(|n-m|) u_n u_m                      (n != m only)
  L3 full        : f(n+m+1) u_n u_m

Einstein spectrum P = delta(w - wE):  F_T(k) = vp^2/(k^2+vp^2), vp = wE/(2 pi T).
Coulomb (step cutoff): M = O[F] - 2*(mus/lam) * uc uc^T,  uc_n = u_n * 1{n < Nc},
  Nc = floor(r*vp).
Tc criterion: lam * lam1(M) = 1.

Nothing is imported from any other repo script.  Two independent assemblers
(explicit double loop / vectorised) are cross-checked in selftest().
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np
from numpy.linalg import eigvalsh

C_INF = 0.1827262477          # (1/2pi) sqrt(g(2)); re-derived in selftest()


# ---------------------------------------------------------------- assemblers
def assemble_loop(fvals, N):
    """Explicit double loop.  fvals[k] = f(k) for k >= 1; fvals[0] unused."""
    u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
    M = np.zeros((N, N))
    for n in range(N):
        for m in range(N):
            x = fvals[n + m + 1] * u[n] * u[m]          # L3
            if n != m:
                x += fvals[abs(n - m)] * u[n] * u[m]    # L2, off-diagonal only
            M[n, m] = x
        s = 0.0
        for k in range(1, n + 1):
            s += 2.0 * fvals[k]
        M[n, n] -= s / (2.0 * n + 1.0)                  # -L1
    return M


def assemble_vec(fvals, N):
    """Vectorised twin of assemble_loop."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    M = (fvals[np.abs(i - j)] + fvals[i + j + 1]) * np.outer(u, u)
    M[n, n] -= fvals[0] * u * u                          # kill L2's diagonal
    csum = np.concatenate(([0.0], np.cumsum(2.0 * fvals[1:N])))
    M[n, n] -= csum[n] / (2.0 * n + 1.0)                 # -L1
    return M


def f_einstein(vp, N):
    k = np.arange(2 * N + 1, dtype=float)
    return vp * vp / (k * k + vp * vp)


def lam1_einstein(lam, mus, T, N=400, r=10.0, loop=False):
    vp = 1.0 / (2.0 * np.pi * T)
    fv = f_einstein(vp, N)
    M = assemble_loop(fv, N) if loop else assemble_vec(fv, N)
    if mus > 0.0:
        Nc = int(np.floor(r * vp))
        if Nc > 0:
            u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
            uc = np.where(np.arange(N) < Nc, u, 0.0)
            M = M - 2.0 * (mus / lam) * np.outer(uc, uc)
    return eigvalsh(M)[-1]


def Tc_einstein(lam, mus, N=400, r=10.0, tol=1e-13):
    """Tc in units of wE.  lam*lam1 is strictly decreasing in T (Theorem 1b)."""
    g = lambda T: lam * lam1_einstein(lam, mus, T, N, r) - 1.0
    lo, hi = 1e-7, 10.0
    assert g(lo) > 0.0 > g(hi), (g(lo), g(hi))
    while hi / lo > 1.0 + tol:                 # bisection in log T
        mid = np.sqrt(lo * hi)
        if g(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return np.sqrt(lo * hi)


def frak_f(lam, mus, N=400, r=10.0):
    """Filling factor f := Tc / (C_INF sqrt(lam <w^2>)), Einstein spectrum."""
    return Tc_einstein(lam, mus, N, r) / (C_INF * np.sqrt(lam))


# ---------------------------------------------------------------- self-tests
def selftest():
    print("=" * 74)
    print("SELFTEST")
    print("=" * 74)

    # (1) two assemblers agree
    N = 60
    fv = f_einstein(1.3, N)
    A, B = assemble_loop(fv, N), assemble_vec(fv, N)
    print(f"[1] loop vs vectorised assembler, max|diff| = {np.abs(A-B).max():.3e}")
    print(f"    symmetry max|M-M^T|               = {np.abs(B-B.T).max():.3e}")

    # (2) g(2) = lam1(O[1/k^2]) -> C_INF, and O[1] = 2uu^T - I  (identity S1)
    for NN in (400, 1600, 6400):
        k = np.arange(2 * NN + 1, dtype=float)
        f2 = np.zeros_like(k)
        f2[1:] = 1.0 / k[1:] ** 2
        g2 = eigvalsh(assemble_vec(f2, NN))[-1]
        print(f"[2] N={NN:5d}  g(2)={g2:.9f}  (1/2pi)sqrt(g2)={np.sqrt(g2)/(2*np.pi):.10f}")
    one = np.ones(2 * 200 + 1)
    u = 1.0 / np.sqrt(2.0 * np.arange(200) + 1.0)
    d = assemble_vec(one, 200) - (2.0 * np.outer(u, u) - np.eye(200))
    print(f"[3] structural identity S1: O[1]=2uu^T-I, max|diff| = {np.abs(d).max():.3e}")

    # (4) weak coupling: Tc/wE -> 0.2530 e^{-1/lam}(1+0.3584 lam)
    pref = (2.0 / np.pi) * np.exp(np.euler_gamma - 1.5)
    print(f"[4] weak-coupling prefactor (2/pi)e^(gamma-3/2) = {pref:.10f}")
    for lam in (0.3, 0.4, 0.5):
        t = Tc_einstein(lam, 0.0)
        w = pref * np.exp(-1.0 / lam) * (1.0 + 0.3584 * lam)
        print(f"    lam={lam:.2f}  Tc/wE={t:.7f}  series={w:.7f}  rel={100*(t/w-1):+.2f}%")

    # (5) strong coupling: vp_c -> 0.87088/sqrt(lam), f -> 1
    for lam in (30.0, 100.0, 1000.0):
        t = Tc_einstein(lam, 0.0)
        vp = 1.0 / (2 * np.pi * t)
        print(f"    lam={lam:7.1f}  vp_c={vp:.6f}  pred={0.87088/np.sqrt(lam):.6f}"
              f"  f={t/(C_INF*np.sqrt(lam)):.6f}")

    # (6) N-convergence at the small-lambda end that PdH actually lives in
    print("[6] N-convergence, lam=0.40, mus=0.085, r=10:")
    for NN in (200, 400, 800, 1600):
        print(f"    N={NN:5d}  Tc/wE={Tc_einstein(0.40, 0.085, N=NN):.9f}")
    print("    lam=0.40, mus=0 :")
    for NN in (200, 400, 800, 1600):
        print(f"    N={NN:5d}  Tc/wE={Tc_einstein(0.40, 0.0, N=NN):.9f}")

    # (7) reproduce the published f-table (independent check of an inherited number)
    print("[7] f(lam, mus=0) vs note deficit_budget_h3s sec.0:")
    ref = {0.5: 0.305, 1.0: 0.627, 1.5: 0.756, 2.0: 0.819,
           3.0: 0.882, 5.0: 0.930, 10.0: 0.965, 30.0: 0.989}
    for lam, rv in ref.items():
        print(f"    lam={lam:5.1f}  mine={frak_f(lam,0.0):.4f}  note={rv:.3f}")


if __name__ == "__main__":
    selftest()
