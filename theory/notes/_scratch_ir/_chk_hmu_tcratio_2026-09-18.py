import numpy as np
from scipy.sparse.linalg import eigsh

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

def top_eig(Mx):
    w = eigsh(Mx, k=1, which='LA', return_eigenvectors=False, maxiter=5000)
    return w[0]

def kappa_mixed(xs, ws, N, Nc, mu_tilde, u_full):
    kmax = 2*N+2
    Fv = np.zeros(kmax+1)
    for xi, wi in zip(xs, ws):
        Fv += wi * f_x(kmax, xi)
    Mx = build_O(Fv, u_full)
    ucv = u_full.copy()
    ucv[Nc:] = 0.0
    Mx = Mx - 2*mu_tilde*np.outer(ucv, ucv)
    return top_eig(Mx)

# reference Einstein solution: lambda=0.50, mu*=0.30 => mu_tilde=0.60, r=5
lam = 0.50
mu_tilde = 0.60
NcE = 774
varpiE = 154.800
xE = varpiE**2
w_ = 0.99
sp = 10.0
x1_ref = xE/(w_ + (1-w_)*sp**2)
x2_ref = sp**2 * x1_ref
# these (x1_ref, x2_ref, NcE) were the peer's fixed-Nc snapshot AT rho=1 (T=Tc(E)).
# For self-consistent Tc(P), let rho = T/Tc(E). Then Nc(rho) = NcE/rho, x_i(rho) = x_i_ref/rho^2.

N_phys = 8000
u_phys = make_u(N_phys)

def lam_kappa_at_rho(rho):
    Nc_r = max(1, int(round(NcE/rho)))
    x1_r = x1_ref/rho**2
    x2_r = x2_ref/rho**2
    kP = kappa_mixed([x1_r, x2_r], [w_, 1-w_], N_phys, Nc_r, mu_tilde, u_phys)
    return lam*kP

# check rho=1 reproduces the earlier snapshot number
print("sanity @ rho=1:", lam_kappa_at_rho(1.0), " (should be ~1.005751)")

# bisect for rho such that lam*kappa(rho) = 1  (Tc(P) boundary)
lo, hi = 0.5, 1.5
flo, fhi = lam_kappa_at_rho(lo)-1, lam_kappa_at_rho(hi)-1
print(f"bracket: f({lo})={flo:+.4e}  f({hi})={fhi:+.4e}")
for _ in range(40):
    mid = 0.5*(lo+hi)
    fm = lam_kappa_at_rho(mid)-1
    if (fm > 0) == (flo > 0):
        lo, flo = mid, fm
    else:
        hi, fhi = mid, fm
rho_star = 0.5*(lo+hi)
print(f"rho* = Tc(P)/Tc(E) = {rho_star:.6f}   (peer claims 1.0279)")
