"""Exact Laplacian identity (analytic) and the Dirichlet relaxation D(f).

With y on Z even about -1/2 (y_{-1-n}=y_n), omega_n=|2n+1|, nu(n)=(omega_n-1)/2,
S(m)=sum_{k<=m} f(k), R(m)=sum_{k>m} f(k), ||f||_1 = S(inf):

  -<w,O[f]w>  = (1/2) sum_{n!=m in Z} f(|n-m|)(y_n-y_m)^2  -  2 sum_n R(nu(n)) y_n^2
  ||w||^2     = (1/2) sum_n omega_n y_n^2          [w_n = y_n/u_n, n>=0]

WAIT: factor. Check numerically.

Then r(f) := -lmin(O[f]) = sup_even  [Dir - 2Neg]/(1/2)<y,Omega y>
and dropping -2Neg gives the upper bound  D(f) := sup_even Dir / ((1/2)<y,Omega y>).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _lemG_reduce import build, ev, step


def even_basis(N):
    E = np.zeros((2 * N, N))
    for n in range(N):
        E[N + n, n] = 1.0
        E[N - 1 - n, n] = 1.0
    return E


def pieces(f, N):
    idx = np.arange(-N, N)
    om = np.abs(2.0 * idx + 1.0)
    nu = (om - 1) / 2
    fk = np.array([0.0] + [f(k) for k in range(1, 4 * N + 2)])  # fk[k]=f(k), fk[0]=0
    T = fk[np.abs(idx[:, None] - idx[None, :])]
    np.fill_diagonal(T, 0.0)
    f1 = fk[1:].sum()
    S = np.cumsum(fk)            # S[m] = sum_{k<=m} f(k)
    Rr = f1 - S[nu.astype(int)]  # R(nu(n))
    return om, T, Rr, f1


def rvals(f, N):
    om, T, Rr, f1 = pieces(f, N)
    E = even_basis(N)
    Om = np.diag(om)
    Lap = 2.0 * f1 * np.eye(2 * N) - T    # PSD on Z (row sums 0); restriction PSD too
    Neg = np.diag(2.0 * Rr)
    B = 0.5 * (E.T @ Om @ E)
    bi = 1.0 / np.sqrt(np.diag(B))

    def gen(M):
        A = E.T @ M @ E
        return float(np.linalg.eigvalsh(A * bi[:, None] * bi[None, :])[-1])

    # exact identity operator:  Lap - Neg  should equal  Phi*I - T  with Phi=2S(nu)
    Phi = np.diag(2.0 * (f1 - Rr))
    assert np.allclose(np.diag(Lap) - np.diag(Neg), np.diag(Phi)), "diag mismatch"
    return gen(Lap - Neg), gen(Lap), gen(Phi - T)


if __name__ == "__main__":
    print("K    r_direct     r_identity   D_K(relax)   check(Phi-T)")
    for K in [1, 2, 3, 5, 8, 16, 32, 64, 128]:
        N = max(24 * K + 900, 1400)
        B = build(step(K), N)
        r_direct = -ev(B["O"])[1]
        r_id, D_K, chk = rvals(step(K), N)
        print(f"{K:4d} {r_direct:12.8f} {r_id:12.8f} {D_K:12.8f} {chk:12.8f}", flush=True)

    print()
    print("other kernels")
    for name, f in [("1/k^2", lambda k: 1.0 / k ** 2),
                    ("1/k^3", lambda k: 1.0 / k ** 3),
                    ("E w=1", lambda k: 1.0 / (k * k + 1.0)),
                    ("E w=7", lambda k: 49.0 / (k * k + 49.0))]:
        N = 1600
        B = build(f, N)
        r_direct = -ev(B["O"])[1]
        l1 = ev(B["O"])[0]
        r_id, D_K, chk = rvals(f, N)
        print(f"{name:8s} r={r_direct:.8f} r_id={r_id:.8f} D={D_K:.8f} "
              f"l1={l1:.6f} f(1)={f(1):.6f}  D/f(1)={D_K/f(1):.6f}", flush=True)
