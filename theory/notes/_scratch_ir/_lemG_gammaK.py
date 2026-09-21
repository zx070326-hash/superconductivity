"""gamma_K = -lambda_min(G_K),  G_K[n,m] = 1_{|n-m|<=K}/sqrt(|2n+1||2m+1|)  on Z.
   r_K <= 1 + gamma_K   (relaxation min(Omega,2K+1) -> Omega).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def GK(K, N):
    """indices n = -N .. N-1 ; banded."""
    idx = np.arange(-N, N)
    w = 1.0 / np.sqrt(np.abs(2.0 * idx + 1.0))
    rows, cols, vals = [], [], []
    a = np.arange(2 * N)
    rows.extend(a); cols.extend(a); vals.extend(w * w)
    for j in range(1, min(K, 2 * N - 1) + 1):
        b = np.arange(2 * N - j)
        v = w[b] * w[b + j]
        rows.extend(b); cols.extend(b + j); vals.extend(v)
        rows.extend(b + j); cols.extend(b); vals.extend(v)
    return sp.coo_matrix((vals, (rows, cols)), shape=(2 * N, 2 * N)).tocsr()


if __name__ == "__main__":
    for K in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]:
        N = min(20 * K + 400, 60000)
        A = GK(K, N)
        lo = spl.eigsh(A, k=1, which="SA", tol=1e-9, maxiter=50000,
                       return_eigenvectors=False)
        print(f"K={K:6d} N={N:7d}  gamma_K={-float(lo[0]):.6f}  "
              f"(r_K <= {1-float(lo[0]):.6f})", flush=True)
