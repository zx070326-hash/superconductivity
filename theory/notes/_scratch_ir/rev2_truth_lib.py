"""rev2: INDEPENDENT float64 H builder written from the entry formulas in the
docstring of p4_interchange_uniformity.py; nothing imported from the certificate
modules."""
import numpy as np

def brk(k, w, order=0):
    k = np.asarray(k, float); k2 = k * k; d = k2 + w * w
    if order == 0: return w * w / d
    if order == 1: return 2 * k2 * w / d**2
    if order == 2: return 2 * k2 * (k2 - 3 * w * w) / d**3
    if order == 4: return 24 * k2 * (-5 * w**4 + 10 * k2 * w * w - k2 * k2) / d**5
    raise ValueError(order)

def H(N, w, order=0):
    n = np.arange(N)
    s = 1.0 / np.sqrt(2 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :]); S = n[:, None] + n[None, :] + 1
    off = np.where(D == 0, 0.0, brk(np.maximum(D, 1), w, order))
    M = (off + brk(S, w, order)) * s[:, None] * s[None, :]
    diag = np.concatenate(([0.0], np.cumsum(2 * brk(np.arange(1, N), w, order)))) / (2 * n + 1.0)
    M[n, n] -= diag
    return M

def lam12(N, w):
    E = np.linalg.eigvalsh(H(N, w)); return E[-1], E[-2]

def lam1_2x2(w):
    return np.linalg.eigvalsh(H(2, w))[-1]
