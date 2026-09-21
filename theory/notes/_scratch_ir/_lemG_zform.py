"""Verify the exact reduction of O[1_{k<=K}] to a Dirichlet-kernel form on Z.

With x_n = u_n w_n (n>=0) and even reflection y_n = x_n (n>=0), y_{-1-n} = x_n,
omega_n := |2n+1| on Z, D_K(j) := 1_{|j|<=K}:

  -<w, A_K w>  =  (1/2) [ <y, min(Omega,2K+1) y> - <y, D_K * y> ]
  ||w||^2      =  (1/2)  <y, Omega y>

so r_K = sup_even-y  <y,(min(Omega,2K+1) - D_K*)y> / <y,Omega y>,
and relaxation  rt_K := sup <y,(Omega - D_K*)y>/<y,Omega y>  >= r_K.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _lemG_reduce import build, ev, step


def ZZ(K, N):
    """Build the Z-side operators on indices n = -N..N-1 (even about -1/2)."""
    idx = np.arange(-N, N)
    om = np.abs(2 * idx + 1.0)
    D = (np.abs(idx[:, None] - idx[None, :]) <= K).astype(float)
    return idx, om, D


def basis_even(N):
    """Orthonormal-ish basis of even (about -1/2) vectors on -N..N-1: columns."""
    E = np.zeros((2 * N, N))
    for n in range(N):
        E[N + n, n] = 1.0        # index n
        E[N - 1 - n, n] = 1.0    # index -1-n
    return E


def rK_Z(K, N, relax=False):
    idx, om, D = ZZ(K, N)
    E = basis_even(N)
    M = np.diag(om if relax else np.minimum(om, 2 * K + 1.0)) - D
    A = E.T @ M @ E
    B = E.T @ np.diag(om) @ E
    # generalized symmetric eigenproblem A z = lam B z ; B is diagonal here
    dB = np.diag(B)
    Bi = 1.0 / np.sqrt(dB)
    C = A * Bi[:, None] * Bi[None, :]
    return float(np.linalg.eigvalsh(C)[-1])


if __name__ == "__main__":
    print("check identity  r_K(direct)  vs  r_K(Z-form)  vs relaxation rt_K")
    for K in [1, 2, 3, 5, 10, 20, 40, 80]:
        N = max(30 * K + 800, 1200)
        B = build(step(K), N)
        r_direct = -ev(B["O"])[1]
        r_z = rK_Z(K, N)
        rt = rK_Z(K, N, relax=True)
        print(f"K={K:4d} N={N:6d}  r_direct={r_direct:.9f}  r_Z={r_z:.9f}  "
              f"diff={abs(r_direct-r_z):.2e}   rt_K={rt:.6f}")
