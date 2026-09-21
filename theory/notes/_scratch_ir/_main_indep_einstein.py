"""Independent Einstein-spectrum Tc(lam, mu*), main session, from the derivation only.

P(dw) = delta(w - w_E)  =>  F_T(k) = varpi^2/(k^2+varpi^2),  varpi = w_E/(2 pi T).
Units: w_E = 1, so T is in units of w_E.  N_c = floor(r*varpi), r = 10.
Reuses NOTHING but numpy/scipy; operator written from P4_main_theorems.md:17.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from scipy.linalg import eigh

def lam1(lam, mus, T, N=400, r=10.0):
    vp = 1.0/(2.0*np.pi*T)
    k  = np.arange(2*N+1, dtype=float)
    F  = vp**2/(k**2+vp**2)          # F[0] = 1 automatically
    n  = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    L1 = np.zeros(N)
    acc = 0.0
    for nn in range(1, N):
        acc += 2.0*F[nn]             # sum_{k=1}^{nn}, INCLUDES k=nn
        L1[nn] = acc/(2.0*nn+1.0)
    i, j = np.meshgrid(n, n, indexing='ij')
    M = (F[np.abs(i-j)] + F[i+j+1])*u[i]*u[j] - np.diag(L1)
    # L2 is OFF-DIAGONAL only: remove the spurious F[0]*u_n^2 the vectorised
    # expression puts on the diagonal (matches P4_main_theorems.md:17).
    M[n, n] -= F[0]*u**2
    Nc = int(np.floor(r*vp))
    if mus > 0 and Nc > 0:
        uc = np.where(n < Nc, u, 0.0)
        M -= 2.0*(mus/lam)*np.outer(uc, uc)
    return eigh(M, eigvals_only=True)[-1]

def Tc(lam, mus, N=400):
    f = lambda T: lam*lam1(lam, mus, T, N) - 1.0
    lo, hi = 1e-6, 5.0
    assert f(lo) > 0 > f(hi), (f(lo), f(hi))
    for _ in range(120):
        mid = np.sqrt(lo*hi)
        if f(mid) > 0: lo = mid
        else: hi = mid
    return np.sqrt(lo*hi)

print(f"{'lam':>5}{'mus':>7}{'Tc(0)':>11}{'Tc(mus)':>11}{'D% (mine)':>11}"
      f"{'D% (agent)':>12}{'diff':>9}")
for lam, mus, dref in [(0.5,0.13,69.972),(1.0,0.13,36.884),(1.5,0.13,26.889),
                       (2.0,0.13,22.844),(3.0,0.13,19.488),(5.0,0.13,17.327),
                       (1.0,0.05,17.036),(1.0,0.20,49.296),(5.0,0.20,23.727)]:
    t0, t1 = Tc(lam,0.0), Tc(lam,mus)
    D = 100*(1-t1/t0)
    print(f"{lam:5.1f}{mus:7.3f}{t0:11.7f}{t1:11.7f}{D:11.3f}{dref:12.3f}{D-dref:+9.3f}")
