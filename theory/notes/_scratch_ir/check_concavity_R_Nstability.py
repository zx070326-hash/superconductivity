"""Independent N-stability check of the concavity-mechanism agent's finding 1:
"R is dominated by the HIGHEST j (j ~ 296-299 of 300)" at (mu~=0.05, N_c=6,
x=0.09).  The bottom of the finite spectrum is where truncation lives
(diagonal -L1_n ~ -2S/(2n+1) at n ~ N), so a bottom-dominated R could be an
N-artifact.  Recompute R, Psi, ratio = R/(2 Psi) and the j-distribution of R
at N = 300, 600, 1200, and cross-check H'' = -2Psi + R by finite differences.
Strictness: numerical-evidence.
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh


def make(N):
    n = np.arange(N); u = 1.0 / np.sqrt(2.0 * n + 1.0)
    ks = np.arange(2 * N + 1, dtype=float)
    return dict(N=N, n=n, u=u, s=ks ** 2, absd=np.abs(n[:, None] - n[None, :]),
                sumd=n[:, None] + n[None, :] + 1, uu=u[:, None] * u[None, :])


def O(f, D):
    N, n, u = D["N"], D["n"], D["u"]
    A = (f[D["absd"]] + f[D["sumd"]]) * D["uu"]
    L1 = np.zeros(N); L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * n + 1] * u ** 2)
    return A


def M_of(x, mut, Nc, D):
    s = D["s"]; f = x / (s + x); f[0] = 1.0
    uc = np.where(D["n"] < Nc, D["u"], 0.0)
    return O(f, D) - 2.0 * mut * np.outer(uc, uc)


def analyse(x, mut, Nc, N):
    D = make(N); s = D["s"]
    M = M_of(x, mut, Nc, D)
    lam, W = eigh(M)
    h = lam[-1]; v = W[:, -1]
    f1 = s / (s + x) ** 2; f1[0] = 0.0          # d/dx of x/(k^2+x); k=0: f=1 const
    psi0 = s / (s + x) ** 3; psi0[0] = 0.0
    Mp = O(f1, D)
    Psi = float(v @ (O(psi0, D) @ v))
    c = W[:, :-1].T @ (Mp @ v)                  # <w_j, M' v>, j != top
    terms = 2.0 * c ** 2 / (h - lam[:-1])
    R = float(terms.sum())
    # finite-difference H''
    def H(xx):
        return float(eigh(M_of(xx, mut, Nc, D), eigvals_only=True)[-1])
    dx = 1e-3 * x
    Hpp = (H(x + dx) - 2 * H(x) + H(x - dx)) / dx ** 2
    # where does R come from?  j indexed so that j=0 is the most negative eigenvalue
    order = np.argsort(terms)[::-1]
    top5 = order[:5]
    frac_bottom_half = float(terms[: (N - 1) // 2].sum() / R)
    frac_top10 = float(terms[-10:].sum() / R)    # the 10 eigenvalues just below h
    gap = h - lam[-2]
    return dict(h=h, gap=gap, Psi=Psi, R=R, ratio=R / (2 * Psi), Hpp=Hpp,
                Hpp_formula=-2 * Psi + R, top5=[(int(N - 1 - j), float(terms[j] / R)) for j in top5],
                frac_bottom_half=frac_bottom_half, frac_top10=frac_top10)


for (x, mut, Nc) in [(0.09, 0.05, 6), (0.04, 0.05, 6), (0.09, 0.13, 6), (0.25, 0.13, 45)]:
    print("=" * 96)
    print(f"x={x} (varpi={math.sqrt(x):.2f})  mu~={mut}  N_c={Nc}")
    print(f"  {'N':>5} {'gap':>9} {'Psi':>10} {'R':>10} {'ratio':>7} {'H'' FD':>10} {'-2Psi+R':>10} "
          f"{'R frac: 10 nearest j':>21} {'bottom half':>12}   top-5 (j from top, frac)")
    for N in (300, 600, 1200):
        a = analyse(x, mut, Nc, N)
        t5 = " ".join(f"({j},{fr:.2f})" for j, fr in a["top5"])
        print(f"  {N:>5} {a['gap']:>9.2e} {a['Psi']:>10.3e} {a['R']:>10.3e} {a['ratio']:>7.3f} "
              f"{a['Hpp']:>10.3e} {a['Hpp_formula']:>10.3e} {a['frac_top10']:>21.3f} "
              f"{a['frac_bottom_half']:>12.3f}   {t5}")
