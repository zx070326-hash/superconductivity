"""INDEPENDENT recomputation of the H3S deficit budget. Nothing imported from
h3s_deficit_budget.py: own file parser, own P-weights (Simpson, not trapezoid),
own kernel assembly, own operator, own Tc root (secant in log T, not bisection
on a geomspace).  Purpose: catch a bug in the headline 199 K.
numerical-evidence (float64)."""
import os, math, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

MEV_K = 11.604518
ROOT = ("D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/"
        "Data-H3S-D3S-PH-EPW-Anharmonic")


def read_a2f(path):
    W, A = [], []
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) < 2:
                continue
            try:
                x = float(t[0]); y = float(t[1])
            except ValueError:
                continue
            W.append(x); A.append(y)
    W = np.asarray(W); A = np.asarray(A)
    k = W > 0.0
    return W[k], A[k]


def simpson_w(x):
    """composite-Simpson weights on a (possibly non-uniform) grid, pairwise."""
    n = len(x); w = np.zeros(n)
    i = 0
    while i + 2 < n:
        a, b, c = x[i], x[i+1], x[i+2]
        h1, h2 = b - a, c - b
        # exact quadratic quadrature on [a,c] through the three nodes
        H = h1 + h2
        w[i]   += H * (2*h1 - h2) / (6*h1)
        w[i+1] += H**3 / (6*h1*h2)
        w[i+2] += H * (2*h2 - h1) / (6*h2)
        i += 2
    if i + 1 < n:                              # odd leftover -> trapezoid
        w[i] += 0.5*(x[i+1]-x[i]); w[i+1] += 0.5*(x[i+1]-x[i])
    return w


def spec(path):
    W, A = read_a2f(path)
    sw = simpson_w(W)
    lam = 2.0 * float(sw @ (A / W))
    p = (2.0 / lam) * (A / W) * sw             # P weights, should already sum to 1
    om2 = float(p @ W**2)
    return lam, om2, W, p, float(p.sum())


def lam1(W, p, T, lam, mus, r, N, sqom2, want_gap=False):
    x = (W / (2*math.pi*T))**2
    k2 = np.arange(2*N+1, dtype=float)**2
    # f(k) = sum_i p_i x_i/(k^2+x_i), computed by an explicit loop over modes
    f = np.zeros(2*N+1)
    for xi, pi in zip(x, p):
        f += pi * xi / (k2 + xi)
    f[0] = float(p.sum())                      # = 1 up to quadrature
    n = np.arange(N); u = 1.0/np.sqrt(2.0*n+1.0)
    Am = np.empty((N, N))
    for i in range(N):                         # explicit double loop, no fancy index
        for j in range(N):
            Am[i, j] = (f[abs(i-j)] + f[i+j+1]) * u[i] * u[j]
    for i in range(N):
        L1 = 0.0
        for kk in range(1, i+1):
            L1 += 2.0*f[kk]
        Am[i, i] = -L1/(2*i+1) + f[2*i+1]*u[i]**2
    if mus == 0.0:
        return float(eigh(Am, eigvals_only=True)[-1])
    Nc = int(math.floor(r * sqom2 / (2*math.pi*T)))
    if Nc <= 0:
        return float(eigh(Am, eigvals_only=True)[-1])
    uc = np.where(n < Nc, u, 0.0)
    M = Am - 2.0*(mus/lam)*np.outer(uc, uc)    # FULL eigh of M, no secular root
    return float(eigh(M, eigvals_only=True)[-1])


def Tc(W, p, lam, mus, r, N, sqom2):
    """secant in s = log T on g(s) = log(lam*lam1) ; root = Tc."""
    g = lambda s: math.log(lam * lam1(W, p, math.exp(s), lam, mus, r, N, sqom2))
    s0, s1 = math.log(0.5), math.log(30.0)
    g0, g1 = g(s0), g(s1)
    assert g0 > 0 > g1, (g0, g1)
    for _ in range(60):
        sm = 0.5*(s0+s1); gm = g(sm)
        if gm > 0: s0, g0 = sm, gm
        else:      s1, g1 = sm, gm
        if s1 - s0 < 1e-13:
            break
    return math.exp(0.5*(s0+s1))


if __name__ == "__main__":
    CINF = 0.1827262477
    print("independent path: Simpson P-weights, explicit loops, full eigh(M), secant root")
    for tag in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa"):
        for sub in ("EPW_anharmonic", "EPW_harmonic"):
            gl = [x for x in glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f") if not x.endswith(".a2fv")]
            if not gl: continue
            lam, om2, W, p, psum = spec(gl[0])
            sq = math.sqrt(om2)
            ceil_K = CINF*math.sqrt(lam*om2)*MEV_K
            N = 250
            tP0 = Tc(W, p, lam, 0.0,  10.0, N, sq)
            tP1 = Tc(W, p, lam, 0.13, 10.0, N, sq)
            print(f"{tag}/{sub:>15} lam={lam:.4f} sq<om2>={sq:.3f} sum(p)={psum:.10f} | "
                  f"ceil={ceil_K:.2f} Tc_P(0)={tP0*MEV_K:.3f} Tc_P(.13)={tP1*MEV_K:.3f} "
                  f"f(.13)={tP1*MEV_K/ceil_K:.4f}", flush=True)
