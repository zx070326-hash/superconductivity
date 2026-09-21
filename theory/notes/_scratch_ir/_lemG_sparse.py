"""Sparse r_K = -lambda_min(O[1_{k<=K}]) to large K, plus trial-vector scheme scan."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def AK_sparse(K, N):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    rows, cols, vals = [], [], []
    # band |n-m| <= K, n != m
    for j in range(1, K + 1):
        if j >= N:
            break
        a = np.arange(N - j)
        v = u[a] * u[a + j]
        rows.extend(a); cols.extend(a + j); vals.extend(v)
        rows.extend(a + j); cols.extend(a); vals.extend(v)
    # L3: n+m+1 <= K  (triangle, includes diagonal contributions handled separately)
    for a in range(min(K, N)):
        bmax = min(K - a - 1, N - 1)
        if bmax < 0:
            continue
        b = np.arange(bmax + 1)
        v = u[a] * u[b]
        rows.extend([a] * len(b)); cols.extend(b); vals.extend(v)
    d = -2.0 * np.minimum(n, K) / (2 * n + 1.0)
    rows.extend(n); cols.extend(n); vals.extend(d)
    A = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()
    return A


def rK(K, N, tol=0):
    A = AK_sparse(K, N)
    w = spl.eigsh(A, k=1, which="SA", tol=tol, maxiter=20000, return_eigenvectors=False)
    return -float(w[0])


def l1K(K, N, tol=0):
    A = AK_sparse(K, N)
    w = spl.eigsh(A, k=1, which="LA", tol=tol, maxiter=20000, return_eigenvectors=False)
    return float(w[0])


if __name__ == "__main__":
    import sys
    Ks = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024]
    for K in Ks:
        for mult in (30, 90):
            N = min(mult * K + 1000, 200000)
            r = rK(K, N)
            print(f"K={K:5d} N={N:7d}  r_K={r:.8f}", flush=True)
