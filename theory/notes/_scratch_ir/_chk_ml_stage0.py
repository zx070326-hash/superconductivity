"""Stage 0: validate the multiband builder against the scalar operator + known constants.
Level: numerical-evidence."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _chk_ml_build import (u_vec, brk, scalar_O, multiband_O, lam_factorized,
                           lam1, lammin, h_of, g2)

FAIL = 0


def chk(name, cond, extra=""):
    global FAIL
    tag = "PASS" if cond else "FAIL"
    if not cond:
        FAIL += 1
    print(f"[{tag}] {name} {extra}")


N = 300

# --- 1. B=1 reduces to the scalar operator, entry by entry
for vp in (0.5, 1.0, 3.0, 12.0):
    A = multiband_O(lam_factorized([[1.0]], vp, N), N)
    ks = np.arange(2 * N + 2)
    Bm = scalar_O(brk(ks, vp), N)
    d = np.max(np.abs(A - Bm))
    chk(f"B=1 matches scalar O  varpi={vp}", d < 1e-13, f"maxdiff={d:.3e}")

# --- 2. g(2) reproduction (rational-exact enclosure in repo: [1.318140211762671, 1.318140212373478])
for Nn in (200, 400, 800, 1600):
    print(f"      g2(N={Nn}) = {g2(Nn):.9f}")
chk("g2 approaches 1.3181402", abs(g2(1600) - 1.3181402) < 3e-3)

# --- 3. scalar lam1 >= |lammin| for the Einstein kernel (needed later)
print("      varpi :   h(varpi)      lammin       |lmin|/lmax")
for vp in (0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0, 100.0):
    ks = np.arange(2 * N + 2)
    M = scalar_O(brk(ks, vp), N)
    a, b = lam1(M), lammin(M)
    print(f"      {vp:6.2f}: {a: .6f}  {b: .6f}   {abs(b)/a:.4f}")

# --- 4. structural decomposition:  OO = sum_{k>=1} O[delta_k] (x) Lt(k)  -  sum_{k>=0} Z_k (x) L(k)
def decomp(lam, N, w):
    B = lam.shape[1]
    u = u_vec(N)
    n = np.arange(N)
    R = np.einsum("kij,j->ki", lam, w) / w[None, :]
    out = np.zeros((N, B, N, B))
    ks = np.arange(2 * N + 2)
    for k in range(1, 2 * N + 2):
        dk = np.zeros(2 * N + 2)
        dk[k] = 1.0
        Ok = scalar_O(dk, N)
        out += Ok[:, None, :, None] * lam[k][None, :, None, :]
    # Z_0 = diag(u^2), Z_k = diag(2 u^2 1_{k<=n})
    for k in range(0, 2 * N + 2):
        Lk = np.diag(R[k]) - lam[k]
        if k == 0:
            z = u ** 2
        else:
            z = 2.0 * u ** 2 * (n >= k)
        if not np.any(z):
            continue
        out[n, :, n, :] -= z[:, None, None] * Lk[None, :, :]
    return out.reshape(N * B, N * B)


Ns = 60
for w in ([1.0, 1.0], [1.0, 2.0]):
    for Lt in ([[1.0, 0.3], [0.3, 0.7]], [[0.4, -0.25], [-0.25, 1.3]]):
        lam = lam_factorized(Lt, 2.0, Ns)
        A = multiband_O(lam, Ns, w=np.array(w))
        Bm = decomp(lam, Ns, np.array(w))
        d = np.max(np.abs(A - Bm))
        chk(f"decomposition identity w={w} Lt={Lt}", d < 1e-11, f"maxdiff={d:.3e}")

# --- 5. Laplacian PSD when Lt elementwise >= 0 (any positive w); and its Dirichlet form
rng = np.random.default_rng(0)
worst = 1e9
worst_neg = 0.0
for _ in range(4000):
    B = rng.integers(2, 6)
    Lt = rng.random((B, B))
    Lt = (Lt + Lt.T) / 2
    w = rng.random(B) + 0.05
    R = (Lt @ w) / w
    L = np.diag(R) - Lt
    ev = np.linalg.eigvalsh((L + L.T) / 2)
    worst = min(worst, ev[0])
    # L is NOT symmetric in general?  check
    worst_neg = max(worst_neg, np.max(np.abs(L - L.T)))
chk("L = diag(R)-Lt is symmetric", worst_neg < 1e-12, f"max asym={worst_neg:.3e}")
chk("L PSD for elementwise-nonneg Lt (4000 random)", worst > -1e-10, f"min eig={worst:.3e}")

print("\nFAILURES:", FAIL)
