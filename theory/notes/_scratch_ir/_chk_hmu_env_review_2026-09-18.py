"""
Independent from-zero verification of (E-env) claim in
worktree-agent-abc83d1b9a2991916:theory/notes/coulomb/hmu_tolmachev_lens_2026-09-18.md

Rebuilds the operator O[f], the Coulomb rank-one term K_c = 2 u_c u_c^T, and
H_mu(x) = lambda_1( O[f_x] - 2*mu_tilde*u_c u_c^T ) completely from scratch, without
importing ANY script from that worktree or from the repo's certificate machinery.

Checks performed:
  (1) self-check: O[1] == 2 u u^T - I  (matrix identity, analytic)
  (2) reproduce h(10), h(20) against P4_einstein_extremal.md (B1) table values
      2.873537558690 and 3.549054516265 (mu_tilde=0, i.e. pure Einstein operator)
  (3) (CVX) sanity: 12 random 2-3 atom measures, mu_tilde>0, check
      lambda1(M(P)) <= sum w_i H_mu(x_i)
  (4) direct reproduction of the claimed violation point:
      lambda=0.50, mu*(w_c)=0.30 => mu_tilde=0.60, r=5, ScE=154.800, N_c=774,
      (w,sp)=(0.99,10) claimed lambdaQ(P)-1 = +5.751e-3
  (5) upper concave envelope (least concave majorant) sandwich check at
      N_c=200, several mu_tilde values: H_mu <= Hhat_mu <= H0 on a grid
  (6) N-shape / Tolmachev lens qualitative check at small-ish N_c (fast)
"""
import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh

def build_O(fvals_diff, fvals_sum, u):
    """
    O[f]_{nm} = u_n u_m [ f(|n-m|) (1-delta_nm) + f(n+m+1) ] - delta_nm * 2 u_n^2 * sum_{k=1}^n f(k)
    fvals_diff, fvals_sum: precomputed via a callable f(k) for k=0..2N
    u: array of u_n
    """
    N = len(u)
    idx = np.arange(N)
    diffmat = np.abs(idx[:, None] - idx[None, :])
    summat = idx[:, None] + idx[None, :] + 1
    M = np.outer(u, u) * (fvals_diff[diffmat] + fvals_sum[summat])
    # remove the erroneous double count of f(|n-m|) on the diagonal (diffmat=0 there)
    # since (1-delta_nm) should zero it; subtract u_n^2*f(0) from diagonal
    diagcorr = (u**2) * fvals_diff[0]
    np.fill_diagonal(M, np.diag(M) - diagcorr)
    # subtract the -delta_nm 2 u_n^2 sum_{k=1}^n f(k) term
    cumf = np.cumsum(fvals_diff[1:N+1])  # cumf[n-1] = sum_{k=1}^n f(k), careful indexing
    # cumf[j] = sum_{k=1}^{j+1} f(k)
    diagextra = np.zeros(N)
    diagextra[1:] = 2 * (u[1:]**2) * cumf[:N-1]
    np.fill_diagonal(M, np.diag(M) - diagextra)
    return M

def f_const1(kmax):
    return np.ones(kmax+1)

def f_x(kmax, x):
    k = np.arange(kmax+1)
    return x / (k**2 + x)

def make_u(N):
    n = np.arange(N)
    return 1.0/np.sqrt(2*n+1)

# ---------- (1) self-check O[1] = 2uu^T - I ----------
N = 60
u = make_u(N)
f1 = f_const1(2*N+2)
O1 = build_O(f1, f1, u)
ref = 2*np.outer(u, u) - np.eye(N)
err1 = np.max(np.abs(O1 - ref))
print(f"(1) self-check O[1]=2uu^T-I : max|diff| = {err1:.3e}")

# ---------- (2) reproduce h(10), h(20) ----------
def h_lambda1(varpi, N, u=None):
    if u is None:
        u = make_u(N)
    x = varpi**2
    fv = f_x(2*N+2, x)
    Mx = build_O(fv, fv, u)
    # top eigenvalue via dense eigh (N moderate) or eigsh for large N
    if N <= 1500:
        w = eigh(Mx, eigvals_only=True, subset_by_index=[N-1, N-1])
        return w[0]
    else:
        w = eigsh(Mx, k=1, which='LA', return_eigenvectors=False)
        return w[0]

for varpi, target, Ns in [(10, 2.873537558690, [400, 800]), (20, 3.549054516265, [800, 1600])]:
    vals = []
    for Nn in Ns:
        vals.append(h_lambda1(varpi, Nn))
    # simple Richardson extrapolation like repo's N=80*w,160*w convention isn't used here;
    # just report raw values and extrapolated estimate assuming O(1/N) error via Richardson-2x
    v1, v2 = vals
    rich = v2 + (v2 - v1)  # crude 2x Richardson assuming N-doubling halves error ~1/N
    print(f"(2) h({varpi}): N={Ns[0]} -> {v1:.10f}, N={Ns[1]} -> {v2:.10f}, "
          f"Richardson~{rich:.10f}, target={target:.10f}, diff(Richardson)={rich-target:.2e}")

# ---------- (3) (CVX) sanity with mu_tilde>0 ----------
rng = np.random.default_rng(12345)
Nc_test = 40
N_test = 800
u_test = make_u(N_test)
uc_test = u_test.copy()
uc_test[Nc_test:] = 0.0
mu_tilde_test = 0.30

def H_mu(x, N, Nc, mu_tilde, u_full=None):
    if u_full is None:
        u_full = make_u(N)
    fv = f_x(2*N+2, x)
    Mx = build_O(fv, fv, u_full)
    ucv = u_full.copy()
    ucv[Nc:] = 0.0
    Mx = Mx - 2*mu_tilde*np.outer(ucv, ucv)
    if N <= 1500:
        w = eigh(Mx, eigvals_only=True, subset_by_index=[N-1, N-1])
        return w[0]
    else:
        w = eigsh(Mx, k=1, which='LA', return_eigenvectors=False)
        return w[0]

def kappa_of_P(xs, ws, N, Nc, mu_tilde, u_full=None):
    """lambda1(M(P)) for P = sum w_i delta_{x_i}, using the linear kernel F_P(k) = sum w_i f_{x_i}(k)."""
    if u_full is None:
        u_full = make_u(N)
    kmax = 2*N+2
    Fv = np.zeros(kmax+1)
    for xi, wi in zip(xs, ws):
        Fv += wi * f_x(kmax, xi)
    Mx = build_O(Fv, Fv, u_full)
    ucv = u_full.copy()
    ucv[Nc:] = 0.0
    Mx = Mx - 2*mu_tilde*np.outer(ucv, ucv)
    if N <= 1500:
        w = eigh(Mx, eigvals_only=True, subset_by_index=[N-1, N-1])
        return w[0]
    else:
        w = eigsh(Mx, k=1, which='LA', return_eigenvectors=False)
        return w[0]

print("\n(3) (CVX) sanity checks, mu_tilde=0.30, Nc=40, N=800:")
worst_slack = None
for i in range(12):
    natoms = rng.integers(2, 4)
    xs = rng.uniform(0.5, 200, size=natoms)
    ws = rng.dirichlet(np.ones(natoms))
    xbar = np.sum(ws*xs)
    kP = kappa_of_P(xs, ws, N_test, Nc_test, mu_tilde_test, u_test)
    rhs = sum(wi*H_mu(xi, N_test, Nc_test, mu_tilde_test, u_test) for xi, wi in zip(xs, ws))
    slack = rhs - kP
    if worst_slack is None or slack < worst_slack:
        worst_slack = slack
    tag = "OK" if slack >= -1e-8 else "VIOLATION"
    print(f"  trial {i}: xbar={xbar:8.2f} kP={kP:10.6f} sum(w_i H_mu(x_i))={rhs:10.6f} slack={slack:+.3e} {tag}")
print(f"  worst slack over 12 trials: {worst_slack:+.3e}  (must be >=0 for CVX to hold)")

# ---------- (4) direct reproduction of the claimed reflex point ----------
print("\n(4) Reproducing claimed violation: lambda=0.50, mu*(wc)=0.30 => mu_tilde=0.60, r=5,")
print("    ScE=154.800, N_c=774, (w,sp)=(0.99,10): claimed lambdaQ(P)-1=+5.751e-3")
mu_tilde_phys = 0.60
Nc_phys = 774
varpiE = 154.800
xbar = varpiE**2
N_phys = 2200  # >> Nc and >> sqrt-scaled support of top eigenvector at these x's

u_phys = make_u(N_phys)

w = 0.99
sp = 10.0
# solve for varpi1, varpi2 s.t. w*x1+(1-w)*x2=xbar and x2=sp^2*x1 (matching sp=varpi2/varpi1)
# w*x1 + (1-w)*sp^2*x1 = xbar => x1 = xbar/(w+(1-w)*sp^2)
x1 = xbar/(w + (1-w)*sp**2)
x2 = sp**2 * x1
varpi1 = np.sqrt(x1)
varpi2 = np.sqrt(x2)
print(f"    solved: x1={x1:.4f} (varpi1={varpi1:.4f}), x2={x2:.4f} (varpi2={varpi2:.4f})")

h_x1 = H_mu(x1, N_phys, Nc_phys, mu_tilde_phys, u_phys)
h_x2 = H_mu(x2, N_phys, Nc_phys, mu_tilde_phys, u_phys)
h_xbar = H_mu(xbar, N_phys, Nc_phys, mu_tilde_phys, u_phys)
lamQ = w*h_x1 + (1-w)*h_x2
excess = lamQ/h_xbar - 1.0
print(f"    H_mu(x1)={h_x1:.6f}  H_mu(x2)={h_x2:.6f}  H_mu(xbar)={h_xbar:.6f}")
print(f"    lambdaQ(P) = w*H(x1)+(1-w)*H(x2) = {lamQ:.6f}")
print(f"    lambdaQ(P)/H_mu(xbar) - 1 = {excess:+.4e}   (claimed +5.751e-3)")

# ---------- (5) upper concave envelope sandwich at Nc=200 ----------
print("\n(5) Upper concave envelope sandwich check, N_c=200:")
Nc5 = 200
N5 = 1200
u5 = make_u(N5)

def upper_concave_envelope(xgrid, ygrid):
    """Return the least concave majorant values on xgrid via upper convex hull (in x,y)."""
    pts = list(zip(xgrid, ygrid))
    # standard upper hull (monotone chain), keep right turns for concave envelope from above
    hull = []
    for p in pts:
        while len(hull) >= 2:
            x1_, y1_ = hull[-2]; x2_, y2_ = hull[-1]
            x3_, y3_ = p
            # cross product to check if hull[-1] is "below" chord hull[-2]-p (i.e. convex from below)
            cross = (x2_-x1_)*(y3_-y1_) - (y2_-y1_)*(x3_-x1_)
            if cross >= 0:
                hull.pop()
            else:
                break
        hull.append(p)
    return hull

for mu_tilde5 in [0.20, 0.60, 1.00]:
    xs_grid = np.geomspace(1.0, 3.0e5, 14)
    Hvals = [H_mu(xv, N5, Nc5, mu_tilde5, u5) for xv in xs_grid]
    H0vals = [H_mu(xv, N5, 1, 0.0, u5) for xv in xs_grid]  # mu_tilde=0 => H0, Nc irrelevant
    hull = upper_concave_envelope(list(xs_grid), Hvals)
    hx = np.array([p[0] for p in hull]); hy = np.array([p[1] for p in hull])
    ok_sandwich = True
    for xv, Hv, H0v in zip(xs_grid, Hvals, H0vals):
        # piecewise-linear interpolation of hull to get envelope value at xv
        Hhat = np.interp(xv, hx, hy)
        if Hv - Hhat > 1e-6 or Hhat - H0v > 1e-6:
            ok_sandwich = False
            print(f"    mu_tilde={mu_tilde5}: SANDWICH VIOLATION at x={xv:.3g}: "
                  f"H_mu={Hv:.6f} Hhat={Hhat:.6f} H0={H0v:.6f}")
    print(f"    mu_tilde={mu_tilde5}: sandwich H_mu<=Hhat<=H0 holds on grid: {ok_sandwich}")

# ---------- (6) N-shape / lens qualitative check ----------
print("\n(6) N-shape check, N_c=200, mu_tilde=0.60 (smaller than the paper's Nc=774 for speed):")
Nc6 = 200
N6 = 1200
u6 = make_u(N6)
mu_tilde6 = 0.60
varpis = np.array([1.0, 5.0, 15.0, 30.0, 45.0, 60.0, 90.0, 130.0, 200.0, 400.0, 900.0])
Hvs = [H_mu(v**2, N6, Nc6, mu_tilde6, u6) for v in varpis]
for v, Hv in zip(varpis, Hvs):
    print(f"    varpi={v:7.1f}  H_mu={Hv:.6f}")
d2 = np.diff(Hvs, 2)
print(f"    2nd differences sign pattern: {['+' if d>0 else '-' for d in d2]}")
