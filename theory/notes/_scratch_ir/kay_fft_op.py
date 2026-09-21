"""Matrix-free (FFT) top-eigenvalue solver for the KAY-Coulomb operator, for the
mu*=0.13-0.15 physical-line check of the weak-coupling extremality failure.

M = O[F] - 2 mu~ u_c u_c^T on l^2(0..N-1),  u_n=(2n+1)^{-1/2}
  off-diag  (F[|n-m|] + F[n+m+1]) u_n u_m
  diag      -L1_n + F[2n+1] u_n^2,   L1_n = (1/(2n+1)) sum_{k=1}^n 2F(k)

The dense build is O(N^2) memory; the physical line needs N ~ 10^4-10^5, so the
matvec is done by circulant embedding (Toeplitz + Hankel), O(N log N), and the
top eigenvalue by Lanczos (eigsh, which='LA').

This file cross-checks FFT-matvec vs the dense operator at small N before any
production run (two independent schemes, per repo discipline).
Strictness: numerical-evidence (float64).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.linalg import LinearOperator, eigsh


def _toep_fft(col, row, L):
    """FFT of the circulant embedding of the Toeplitz matrix with given first
    column/row (col[0] must equal row[0]); L >= 2N."""
    N = col.size
    c = np.zeros(L)
    c[:N] = col
    c[L - N + 1:] = row[:0:-1]
    return np.fft.rfft(c)


class KAY:
    """M = O[F] - 2 mu~ u_c u_c^T as a LinearOperator."""

    def __init__(self, F, N, mut, Nc):
        # F must be defined on 0..2N-1 at least
        self.N = N
        n = np.arange(N)
        self.u = 1.0 / np.sqrt(2.0 * n + 1.0)
        L = 1
        while L < 2 * N:
            L *= 2
        self.L = L
        h = F[1:2 * N]                      # h[j] = F[j+1], j=0..2N-2
        self.fT = _toep_fft(F[:N].copy(), F[:N].copy(), L)          # T[n,m]=F[|n-m|]
        self.fH = _toep_fft(h[N - 1:2 * N - 1].copy(), h[N - 1::-1].copy(), L)
        L1 = np.zeros(N)
        if N > 1:
            L1[1:] = np.cumsum(2.0 * F[1:N])[:N - 1] / (2 * n[1:] + 1)
        # diagonal correction: replace u_n^2 (F[0]+F[2n+1]) by -L1_n + F[2n+1] u_n^2
        self.dcorr = -L1 - self.u ** 2 * F[0]
        self.uc = np.where(n < Nc, self.u, 0.0)
        self.mut = mut

    def _toep(self, fc, x):
        xp = np.zeros(self.L)
        xp[:self.N] = x
        return np.fft.irfft(fc * np.fft.rfft(xp), self.L)[:self.N]

    def matvec(self, x):
        x = np.asarray(x, float).ravel()
        w = self.u * x
        y = self._toep(self.fT, w) + self._toep(self.fH, w[::-1])
        y = self.u * y + self.dcorr * x
        return y - 2.0 * self.mut * self.uc * float(self.uc @ x)

    def op(self):
        return LinearOperator((self.N, self.N), matvec=self.matvec, dtype=float)

    def top(self, tol=1e-10, k=1):
        v0 = self.u.copy()
        return float(eigsh(self.op(), k=k, which="LA", tol=tol, v0=v0,
                           return_eigenvectors=False)[-1])


def dense(F, N, mut, Nc):
    n = np.arange(N); u = 1.0 / np.sqrt(2.0 * n + 1.0)
    A = (F[np.abs(n[:, None] - n[None, :])] + F[n[:, None] + n[None, :] + 1]) * (u[:, None] * u[None, :])
    L1 = np.zeros(N); L1[1:] = np.cumsum(2.0 * F[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(A, -L1 + F[2 * n + 1] * u ** 2)
    uc = np.where(n < Nc, u, 0.0)
    return A - 2.0 * mut * np.outer(uc, uc)


def Fspec(atoms, N):
    s = np.arange(2 * N + 1, dtype=float) ** 2
    F = np.zeros_like(s)
    for w, x in atoms:
        F += w * x / (s + x)
    F[0] = sum(w for w, _ in atoms)
    return F


def pair(xbar, w, sp):
    x1 = xbar / (w + (1 - w) * sp)
    return [(w, x1), (1 - w, sp * x1)]


if __name__ == "__main__":
    print("cross-check FFT matvec vs dense operator")
    rng = np.random.default_rng(0)
    for (N, vp, mut, Nc, atoms) in [
        (200, 3.0, 0.26, 15, [(1.0, 9.0)]),
        (400, 7.0, 0.26, 35, pair(49.0, 0.99, 100.0)),
        (777, 20.0, 0.40, 100, pair(400.0, 0.9, 1e3)),
    ]:
        F = Fspec(atoms, N)
        K = KAY(F, N, mut, Nc)
        A = dense(F, N, mut, Nc)
        x = rng.standard_normal(N)
        e_mv = np.abs(K.matvec(x) - A @ x).max() / np.abs(A @ x).max()
        t_d = float(eigh(A, eigvals_only=True, subset_by_index=[N - 1, N - 1])[0])
        t_f = K.top()
        print(f"  N={N:4d} Nc={Nc:4d} mu~={mut}  matvec rel err {e_mv:.2e}   "
              f"top dense {t_d:.12f}  fft {t_f:.12f}  diff {abs(t_d-t_f):.2e}")
