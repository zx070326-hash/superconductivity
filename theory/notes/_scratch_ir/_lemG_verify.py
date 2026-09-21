"""INDEPENDENT end-to-end verification of the Lemma (G) partial theorems.

Rebuilt from scratch (own operator builder, own Abel decomposition) -- does NOT
reuse _lemG_cert.py logic.  Checks, on a wide kernel family:

  (V1)  lmax(L1[f])  <= sum_K a_K theta_K                   [Abel majorant]
  (V2)  lam1(O[f])   >= sum_K a_K q_K(v)                    [trial vector, v=(9,5,2,1)]
  (V3)  |lmin(O[f])| <= lam1 + 2 lmax(L1)                   [Perron / abs-value]
  (V4)  |lmin(O[f])| <= (1+2A) lam1,  A = 663530/1174507    [main uniform bound]
  (V5)  |lmin(O[f])| <= max(f(1), 2 C_f),  C_f = -inf_th sum_k f(k) cos(k th)
  (V6)  |lmin(O[f])| <= 3 lam1                              [crude analytic]

Direction note: at finite N both |lmin_N| <= |lmin| and lam1_N <= lam1 (compression),
so these runs are CONSISTENCY checks, not proofs; the proof is the analytic chain.
A poisoned control (A halved) must FAIL, otherwise the test is vacuous.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from fractions import Fraction as Fr

A_CERT = Fr(663530, 1174507)
C_UNIF = float(1 + 2 * A_CERT)
TRIAL_T = [9.0, 5.0, 2.0, 1.0]


def op(fk, N):
    """fk: array with fk[k]=f(k), k=1..; build O[f] head of size N (own code)."""
    n = np.arange(N)
    u = (2.0 * n + 1.0) ** -0.5
    def f(k):
        k = np.asarray(k)
        out = np.zeros(k.shape, float)
        m = k < len(fk)
        out[m] = fk[k[m]]
        return out
    I, J = np.meshgrid(n, n, indexing="ij")
    O = (f(np.abs(I - J)) * (I != J) + f(I + J + 1)) * np.outer(u, u)
    S = np.concatenate([[0.0], np.cumsum(fk[1:N])])
    L1 = 2.0 * S / (2 * n + 1.0)
    O[np.diag_indices(N)] -= L1
    return O, L1, u


def abel(fk):
    """a_K = f(K)-f(K+1) >= 0 for nonincreasing f."""
    g = np.append(fk[1:], 0.0)
    return np.maximum(g[:-1] - g[1:], 0.0), np.arange(1, len(g))


def qK_float(t, K):
    p = len(t)
    n = np.arange(p)
    u = (2.0 * n + 1.0) ** -0.5
    v = u * np.array(t)
    fk = np.zeros(2 * p + 3)
    fk[1:K + 1] = 1.0
    O, _, _ = op(fk, p)
    return float(v @ O @ v) / float(v @ v)


def Cf(fk, M=1 << 16):
    """C_f = -inf_theta sum_{k>=1} f(k) cos(k theta), evaluated on an M-point grid.
    A grid UNDERestimates the sup, so V5 becomes a harder (conservative) test."""
    g = np.zeros(M)
    L = min(len(fk) - 1, M - 1)
    g[1:L + 1] = fk[1:L + 1]
    return float((-np.fft.rfft(g).real).max())


def check(name, fk, N=800, A=A_CERT, verbose=True):
    O, L1, u = op(fk, N)
    w = np.linalg.eigvalsh(O)
    l1, r = w[-1], -w[0]
    lL1 = L1.max()
    a, Ks = abel(fk)
    theta = 2.0 * Ks / (2.0 * Ks + 1.0)
    maj = float(a @ theta)
    qs = np.array([qK_float(TRIAL_T, int(K)) for K in Ks])
    low = float(a @ qs)
    C = Cf(fk)
    tol = 1e-9
    V = [lL1 <= maj + tol, l1 >= low - tol, r <= l1 + 2 * lL1 + tol,
         r <= (1 + 2 * float(A)) * l1 + tol, r <= max(fk[1], 2 * C) + tol,
         r <= 3 * l1 + tol]
    if verbose:
        print(f"{name:20s} r={r:9.6f} l1={l1:9.6f} r/l1={r/l1:7.5f} | "
              f"L1={lL1:8.6f}<=maj={maj:8.6f} | low={low:8.6f}<=l1 | "
              f"2C={2*C:9.6f} f1={fk[1]:.4f} | V1..V6={''.join('T' if x else 'F' for x in V)}")
    return V, r, l1


if __name__ == "__main__":
    fams = []
    for K in [1, 2, 3, 5, 10, 30, 100, 300]:
        fk = np.zeros(1400); fk[1:K + 1] = 1.0
        fams.append((f"step K={K}", fk))
    k = np.arange(0, 1400, dtype=float); k[0] = 1
    for nm, arr in [("1/k^2", 1 / k ** 2), ("1/k^3", 1 / k ** 3),
                    ("1/k^1.1", 1 / k ** 1.1), ("1/k", 1 / k),
                    ("E w=1", 1 / (k * k + 1)), ("E w=7", 49 / (k * k + 49)),
                    ("E w=100", 1e4 / (k * k + 1e4)), ("exp(-k)", np.exp(-k)),
                    ("exp(-k/20)", np.exp(-k / 20))]:
        a = arr.copy(); a[0] = 0.0
        fams.append((nm, a))
    allV = []
    for nm, fk in fams:
        V, _, _ = check(nm, fk)
        allV.append(V)

    print("\n--- 300 random nonneg nonincreasing kernels (V1..V6) ---")
    rng = np.random.default_rng(918)
    fails = np.zeros(6, int)
    worst = 0.0
    for _ in range(300):
        p = int(rng.integers(1, 30))
        fk = np.zeros(1400)
        for _ in range(p):
            K = int(rng.integers(1, 400))
            fk[1:K + 1] += rng.random() ** rng.uniform(0.3, 4)
        V, r, l1 = check("", fk, N=600, verbose=False)
        fails += ~np.array(V)
        worst = max(worst, r / l1)
    print(f"failures per check V1..V6: {fails.tolist()}   max r/l1 seen: {worst:.6f}")

    print("\n--- POISON control: A halved (must FAIL) ---")
    fk = np.zeros(1400); fk[1] = 1.0
    V, r, l1 = check("step K=1 POISON", fk, A=A_CERT / 2)
    print(f"V4 under poisoned A: {V[3]}  (expected False)")
