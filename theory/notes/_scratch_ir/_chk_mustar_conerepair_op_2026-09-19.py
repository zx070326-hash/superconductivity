"""From-scratch kernel-generic operator O[f] + Coulomb pencil, for the mu*>0
Einstein-extremality "cone repair" (Theorem CR) work of 2026-09-19.

Gate 0 (non-negotiable): the Einstein operator built here must reproduce the
certified baseline `p4_interchange_uniformity.H_mat` bit-for-bit, and the
Coulomb operator must equal 2 u_c u_c^T elementwise.

No production file is imported except the baseline H_mat (read-only, gate only).
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "theory", "derivations"))


# ---------------------------------------------------------------- operator ---
def O_of_kernel(N, fvals):
    """O[f] on C^N.  fvals[k] = f(k) for k = 0 .. 2N-1 (f(0) never used).

    O[f] = -L1[f] + L2[f] + L3[f]  (P4_main_theorems.md sec.0):
      L2 off-diagonal  u_n u_m f(|n-m|)
      L3              u_n u_m f(n+m+1)
      L1 diagonal     (2n+1)^-1 sum_{k=1}^{n} 2 f(k)
    """
    n = np.arange(N)
    inv = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(inv, inv)
    A = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    off = fvals[A].astype(np.float64).copy()
    np.fill_diagonal(off, 0.0)
    M = (off + fvals[S]) * nrm
    cum = np.concatenate(([0.0], np.cumsum(2.0 * fvals[1:N])))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M


def brk(k, w):
    """[[k]](w) = w^2 / (k^2 + w^2)."""
    k = np.asarray(k, dtype=np.float64)
    return w * w / (k * k + w * w)


def kernel_einstein(N, w):
    return brk(np.arange(2 * N + 2), w)


def kernel_measure(N, ws, probs):
    """F_P(k) = sum_i p_i [[k]](w_i)."""
    k = np.arange(2 * N + 2)
    out = np.zeros(len(k))
    for p, w in zip(probs, ws):
        out += p * brk(k, w)
    return out


def u_vec(N):
    return 1.0 / np.sqrt(2.0 * np.arange(N) + 1.0)


def u_c_vec(N, Nc):
    u = u_vec(N)
    u[Nc:] = 0.0
    return u


def pencil(A, N, Nc, mut):
    uc = u_c_vec(N, Nc)
    return A - 2.0 * mut * np.outer(uc, uc)


def top_pair(M):
    vals, vecs = np.linalg.eigh(M)
    v = vecs[:, -1]
    if v.sum() < 0:
        v = -v
    return vals[-1], v, vals[0]


# ------------------------------------------------------------- Abel / cone ---
def abel(v):
    """a_L = v_{L-1} - v_L for L = 1..N  (v_N := 0).  v = sum_L a_L e_L."""
    N = len(v)
    a = np.empty(N)
    a[: N - 1] = v[: N - 1] - v[1:]
    a[N - 1] = v[N - 1]
    return a  # a[L-1] is a_L


def from_abel(a):
    """Inverse of `abel`: v_n = sum_{L > n} a_L."""
    return np.cumsum(a[::-1])[::-1]


# ---------------------------------------------------------------- Gate 0 -----
def gate0(verbose=True):
    import p4_interchange_uniformity as base

    checks, fails = 0, 0

    def chk(name, cond):
        nonlocal checks, fails
        checks += 1
        if not cond:
            fails += 1
            print("  FAIL", name)
        elif verbose:
            print("  ok  ", name)

    for N in (5, 17, 60, 200):
        for w in (0.3, 0.7, 1.0, 2.0, 5.0, 13.0):
            mine = O_of_kernel(N, kernel_einstein(N, w))
            theirs = base.H_mat(N, w)
            chk(f"H_mat bitwise N={N} w={w}", np.array_equal(mine, theirs))

    # (S1)  O[1] = 2 u u^T - I
    for N in (7, 40):
        ones = np.ones(2 * N + 2)
        mine = O_of_kernel(N, ones)
        u = u_vec(N)
        chk(f"(S1) O[1]=2uu^T-I  N={N}", np.allclose(mine, 2 * np.outer(u, u) - np.eye(N), atol=1e-14))

    # K_c = P_Nc (O[1] + I) P_Nc = 2 u_c u_c^T
    for N, Nc in ((40, 1), (40, 6), (40, 25)):
        ones = np.ones(2 * N + 2)
        P = np.zeros((N, N))
        P[np.diag_indices(N)] = [1.0 if i < Nc else 0.0 for i in range(N)]
        lhs = P @ (O_of_kernel(N, ones) + np.eye(N)) @ P
        uc = u_c_vec(N, Nc)
        chk(f"K_c=2u_cu_c^T N={N} Nc={Nc}", np.allclose(lhs, 2 * np.outer(uc, uc), atol=1e-14))

    # linearity of O in the kernel: O[F_P] = sum_i p_i O[[[.]](w_i)]
    N = 60
    ws, ps = (0.5, 3.0), (0.4, 0.6)
    lhs = O_of_kernel(N, kernel_measure(N, ws, ps))
    rhs = sum(p * O_of_kernel(N, kernel_einstein(N, w)) for p, w in zip(ps, ws))
    chk("O linear in kernel", np.allclose(lhs, rhs, atol=1e-14))

    # Abel round-trip
    rng = np.random.default_rng(0)
    for _ in range(5):
        v = rng.normal(size=50)
        chk("Abel round-trip", np.allclose(from_abel(abel(v)), v, atol=1e-14))

    # e_L representation: v = sum a_L e_L
    v = rng.normal(size=30)
    a = abel(v)
    rec = np.zeros(30)
    for L in range(1, 31):
        e = np.zeros(30)
        e[:L] = 1.0
        rec += a[L - 1] * e
    chk("v = sum a_L e_L", np.allclose(rec, v, atol=1e-13))

    print(f"Gate 0: {checks - fails}/{checks} passed")
    return fails == 0


if __name__ == "__main__":
    ok = gate0(verbose="-v" in sys.argv)
    sys.exit(0 if ok else 1)
