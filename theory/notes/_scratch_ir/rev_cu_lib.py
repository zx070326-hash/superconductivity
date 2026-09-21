"""Independent re-implementation for the P4-D Coulomb uniqueness review.
Built from the ENTRY FORMULAS (KAY II), not by importing the module under review.
Only p4_interchange_uniformity.H_mat is imported, and only as a cross-check."""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
import numpy as np


def brk(K, w):
    """[[k]](w) = w^2/(k^2+w^2); [[0]] = 1."""
    K = np.asarray(K, dtype=float)
    return np.where(K == 0, 1.0, w * w / (K * K + w * w))


def brk_wder(K, w):
    """w * d/dw [[k]](w) = 2 k^2 w^2/(k^2+w^2)^2 ; 0 at k=0."""
    K = np.asarray(K, dtype=float)
    d = K * K + w * w
    return np.where(K == 0, 0.0, 2.0 * K * K * w * w / (d * d))


def Ocal(N, fvals):
    """O[f] from the elementwise definition, fvals[k] = f(k) for k=0..2N."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    A = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    off = fvals[A].astype(float).copy()
    np.fill_diagonal(off, 0.0)
    M = (off + fvals[S]) * np.outer(u, u)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * fvals[1:N])))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M


def myH(N, w):
    return Ocal(N, brk(np.arange(2 * N + 2), w))


def myOf(N, w):
    return Ocal(N, brk_wder(np.arange(2 * N + 2), w))


def uc(N, Nc):
    u = 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)
    u[min(Nc, N):] = 0.0
    return u


class RCell:
    """Threshold ingredients, recomputed independently."""

    def __init__(self, w, Nc, N):
        self.w, self.Nc, self.N = w, Nc, N
        A = myH(N, w)
        Of = myOf(N, w)
        e, V = np.linalg.eigh(A)
        self.e, self.V = e, V
        self.G = V.T @ Of @ V
        self.c = V.T @ uc(N, Nc)
        self.mu1, self.mu2, self.mu3 = float(e[-1]), float(e[-2]), float(e[-3])
        self.gam = self.mu1 - self.mu2
        self.q1 = float(self.G[-1, -1])
        self.b = float(np.linalg.norm(self.G[:-1, -1]))
        self.m = max(0.0, -float(np.linalg.eigvalsh(Of)[0]))
        self.tmax = (np.inf if self.m == 0 else
                     (-self.b + np.sqrt(self.b ** 2 + self.m * self.q1)) / self.m)
        self.c1 = abs(float(self.c[-1]))
        self.c2 = abs(float(self.c[-2]))
        self.up = float(np.linalg.norm(self.c[:-1]))
        if np.isinf(self.tmax):
            self.thr = self.gam / (2 * self.c1 ** 2)
        else:
            self.thr = (self.tmax * self.gam /
                        (2 * self.c1 * self.up + 2 * self.tmax * self.c1 ** 2))
        self.A, self.Of = A, Of

    def lam_dense(self, mut):
        u = uc(self.N, self.Nc)
        return float(np.linalg.eigvalsh(self.A - 2.0 * mut * np.outer(u, u))[-1])

    def vhat_dense(self, mut):
        u = uc(self.N, self.Nc)
        e, V = np.linalg.eigh(self.A - 2.0 * mut * np.outer(u, u))
        return float(e[-1]), float(e[-2]), V[:, -1]

    def q_dense(self, mut):
        """<vhat, O[f] vhat> via a DENSE eigh of the pencil (no secular engine)."""
        _, _, v = self.vhat_dense(mut)
        return float(v @ (self.Of @ v))

    def tan_theta_dense(self, mut):
        _, _, v = self.vhat_dense(mut)
        c = abs(float(self.V[:, -1] @ v))
        return np.sqrt(max(0.0, 1.0 - c * c)) / max(c, 1e-300)

    def kappa(self, mut):
        return 2 * mut * self.c1 ** 2 / self.gam
