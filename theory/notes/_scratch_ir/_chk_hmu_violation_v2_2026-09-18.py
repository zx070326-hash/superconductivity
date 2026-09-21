import numpy as np
from scipy.sparse.linalg import eigsh
import time

def build_O(fvals, u):
    N = len(u)
    idx = np.arange(N)
    diffmat = np.abs(idx[:, None] - idx[None, :])
    summat = idx[:, None] + idx[None, :] + 1
    M = np.outer(u, u) * (fvals[diffmat] + fvals[summat])
    diagcorr = (u**2) * fvals[0]
    np.fill_diagonal(M, np.diag(M) - diagcorr)
    cumf = np.cumsum(fvals[1:N+1])
    diagextra = np.zeros(N)
    diagextra[1:] = 2 * (u[1:]**2) * cumf[:N-1]
    np.fill_diagonal(M, np.diag(M) - diagextra)
    return M

def make_u(N):
    n = np.arange(N)
    return 1.0/np.sqrt(2*n+1)

def f_x(kmax, x):
    k = np.arange(kmax+1)
    return x / (k**2 + x)

def top_eig(Mx, w0=None):
    N = Mx.shape[0]
    v0 = w0
    w = eigsh(Mx, k=1, which='LA', return_eigenvectors=False, maxiter=5000)
    return w[0]

def H_mu(x, N, Nc, mu_tilde, u_full):
    fv = f_x(2*N+2, x)
    Mx = build_O(fv, u_full)
    ucv = u_full.copy()
    ucv[Nc:] = 0.0
    Mx = Mx - 2*mu_tilde*np.outer(ucv, ucv)
    return top_eig(Mx)

def kappa_twoatom(x1, x2, w, N, Nc, mu_tilde, u_full):
    """lambda1 of the TRUE mixed kernel M(P) for P = w*delta_x1 + (1-w)*delta_x2."""
    kmax = 2*N+2
    Fv = w*f_x(kmax, x1) + (1-w)*f_x(kmax, x2)
    Mx = build_O(Fv, u_full)
    ucv = u_full.copy()
    ucv[Nc:] = 0.0
    Mx = Mx - 2*mu_tilde*np.outer(ucv, ucv)
    return top_eig(Mx)

mu_tilde_phys = 0.60
Nc_phys = 774
lam = 0.50
varpiE = 154.800
xbar = varpiE**2

cases = [(0.99, 3.0), (0.99, 10.0), (0.99, 30.0), (0.995, 10.0), (0.995, 30.0)]

for N_phys in [3000, 6000, 10000]:
    u_phys = make_u(N_phys)
    h_xbar = H_mu(xbar, N_phys, Nc_phys, mu_tilde_phys, u_phys)
    print(f"\n=== N={N_phys}: H_mu(xbar)={h_xbar:.6f}  (lambda*H should ~1: {lam*h_xbar:.6f}) ===")
    for w, sp in cases:
        x1 = xbar/(w + (1-w)*sp**2)
        x2 = sp**2 * x1
        t0 = time.time()
        kP = kappa_twoatom(x1, x2, w, N_phys, Nc_phys, mu_tilde_phys, u_phys)
        excess = lam*kP - 1.0
        dt = time.time()-t0
        print(f"  w={w:.3f} sp={sp:5.1f} x1={x1:.3f}(vp1={x1**0.5:.2f}) "
              f"x2={x2:.3f}(vp2={x2**0.5:.2f}) kappa(P)={kP:.6f} "
              f"lambda*kappa-1={excess:+.4e}  [{dt:.1f}s]")
