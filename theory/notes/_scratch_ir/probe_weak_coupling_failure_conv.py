"""N-convergence of the weak-coupling Tc-ordering violation found in
probe_weak_coupling_failure.py (lambda=0.50, mu=0.30 at omega_c=5 omega_rms):
at N = 10 varpi the kernel tail f(N) = varpi^2/(N^2+varpi^2) ~ 1%, so the
truncation is NOT obviously negligible at the 0.5% margin observed.  Repeat at
N = 10, 20, 40 varpi, re-solving Tc(E) at each N, and also the neighbouring
'no' row (lambda=0.55) and a finer lambda ladder to locate the boundary.
Strictness: numerical-evidence (float64).
"""
import os, math, time
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

_c = {}


def setup(N):
    if N in _c:
        return _c[N]
    n = np.arange(N); u = 1.0 / np.sqrt(2.0 * n + 1.0)
    ks = np.arange(2 * N + 1, dtype=float)
    _c[N] = dict(N=N, n=n, u=u, s=ks ** 2, absd=np.abs(n[:, None] - n[None, :]),
                 sumd=n[:, None] + n[None, :] + 1, uu=u[:, None] * u[None, :])
    return _c[N]


def top(M):
    N = M.shape[0]
    return float(eigh(M, eigvals_only=True, subset_by_index=[N - 1, N - 1])[0])


def Q(atoms, mut, vp, r, Nfac):
    Nc = max(1, int(round(r * vp)))
    N = int(max(200, 1.7 * Nc, Nfac * vp))
    D = setup(N)
    s = D["s"]; F = np.zeros_like(s)
    for w, x in atoms:
        F += w * x / (s + x)
    F[0] = sum(w for w, _ in atoms)
    n, u = D["n"], D["u"]
    A = (F[D["absd"]] + F[D["sumd"]]) * D["uu"]
    L1 = np.zeros(N); L1[1:] = np.cumsum(2.0 * F[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(A, -L1 + F[2 * n + 1] * u ** 2)
    uc = np.where(n < Nc, u, 0.0)
    return top(A - 2.0 * mut * np.outer(uc, uc)), Nc, N


def pair(xbar, w, sp):
    x1 = xbar / (w + (1 - w) * sp)
    return [(w, x1), (1 - w, sp * x1)]


r = 5.0
mu = 0.30
mus = mu / (1 + mu * math.log(r))
print(f"mu={mu}, omega_c={r:g} omega_rms, mu*_rms={mus:.4f}")
print(f"{'lambda':>7} {'2mu*/lam':>9} {'Nfac':>5} {'varpi_E':>9} {'Nc':>5} {'N':>6} "
      f"{'Q(P)-1/lam':>12} {'rel':>10} {'viol':>5} {'t[s]':>6}")
for lam in (0.55, 0.52, 0.50, 0.48):
    mut = mu / lam
    for Nfac in (10, 20, 40):
        t0 = time.time()
        def g(vp):
            return Q([(1.0, vp * vp)], mut, vp, r, Nfac)[0] - 1.0 / lam
        lo, hi = 2.0, 8.0
        while g(hi) < 0:
            lo, hi = hi, hi * 1.5
        vpE = brentq(g, lo, hi, xtol=1e-3)
        qE, Nc, N = Q([(1.0, vpE * vpE)], mut, vpE, r, Nfac)
        best = max(Q(pair(vpE * vpE, 0.99, sp), mut, vpE, r, Nfac)[0]
                   for sp in (10.0, 1e2, 1e3))
        d = best - 1.0 / lam
        print(f"{lam:>7.2f} {2*mus/lam:>9.3f} {Nfac:>5} {vpE:>9.2f} {Nc:>5} {N:>6} "
              f"{d:>12.4e} {d*lam:>10.3e} {'YES' if d > 0 else 'no':>5} {time.time()-t0:>6.0f}")
        _c.clear()
    print()
