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

def H_mu(x, N, Nc, mu_tilde, u_full):
    fv = f_x(2*N+2, x)
    Mx = build_O(fv, u_full)
    ucv = u_full.copy()
    ucv[Nc:] = 0.0
    Mx = Mx - 2*mu_tilde*np.outer(ucv, ucv)
    w = eigsh(Mx, k=1, which='LA', return_eigenvectors=False, maxiter=5000)
    return w[0]

mu_tilde_phys = 0.60
Nc_phys = 774
varpiE = 154.800
xbar = varpiE**2
w = 0.99
sp = 10.0
x1 = xbar/(w + (1-w)*sp**2)
x2 = sp**2 * x1

for N_phys in [3000, 5000, 8000]:
    t0 = time.time()
    u_phys = make_u(N_phys)
    h_x1 = H_mu(x1, N_phys, Nc_phys, mu_tilde_phys, u_phys)
    h_x2 = H_mu(x2, N_phys, Nc_phys, mu_tilde_phys, u_phys)
    h_xbar = H_mu(xbar, N_phys, Nc_phys, mu_tilde_phys, u_phys)
    lamQ = w*h_x1 + (1-w)*h_x2
    excess = lamQ/h_xbar - 1.0
    dt = time.time()-t0
    print(f"N={N_phys}: H(x1)={h_x1:.6f} H(x2)={h_x2:.6f} H(xbar)={h_xbar:.6f} "
          f"excess={excess:+.4e}  [{dt:.1f}s]")
