"""Multiband (matrix-valued) lift of the scalar operator O[f].

Level: numerical-evidence (float64 point evaluation, finite N truncation).

Scalar convention (repo, P4_main_theorems.md Sec.0 + P4_uniqueness_strategy.md Sec.2.A):
  u_n = (2n+1)^{-1/2}
  O[f]_{nm} = u_n u_m [ f(|n-m|) 1_{n!=m} + f(n+m+1) ]  -  delta_{nm} 2 u_n^2 sum_{k=1}^n f(k)
  (the |n-m|=0 entry is absent: f is defined on k>=1 only)

Multiband: derived from scratch from the linearized multiband Eliashberg equations
  Z_i(n) = 1 + u_n^2 ( R_i(0) + 2 sum_{k=1}^n R_i(k) ),    R_i(k) = sum_j lam_ij(k)
  Z_i(n) Delta_i(n) = sum_{j,m} [ lam_ij(|n-m|) + lam_ij(n+m+1) ] Delta_j(m) / (2m+1)
with psi_i(n) = u_n Delta_i(n).  Eigenvalue-1 criterion.

  OO_{(n,i),(m,j)} = u_n u_m [ lam_ij(|n-m|) 1_{n!=m} + lam_ij(n+m+1) ]
                     + delta_{nm} [ u_n^2 lam_ij(0) - delta_ij u_n^2 (R_i(0) + 2 sum_{k=1}^n R_i(k)) ]

NOTE the k=0 term does NOT cancel when the band matrix is non-diagonal: Z carries the
ROW SUM R_i, the pairing carries the MATRIX lam_ij.  That mismatch is the whole story.

Symmetrization: lam_ij = N_j V_ij with V symmetric (detailed balance N_i lam_ij = N_j lam_ji).
We parametrize by the SYMMETRIC Lt_ij = sqrt(N_i N_j) V_ij and weights w_i = sqrt(N_i);
then lam_ij = sqrt(N_j/N_i) Lt_ij and R_i = (Lt w)_i / w_i.
Conjugating by diag(w) makes the pairing block symmetric and leaves Z untouched
(band-diagonal), so OO is similar to a symmetric matrix; we build that symmetric form.
"""
import numpy as np

OMEGA_KEY = "w"


def u_vec(N):
    n = np.arange(N)
    return 1.0 / np.sqrt(2.0 * n + 1.0)


def brk(k, varpi):
    """[[k]](varpi) = varpi^2/(k^2+varpi^2).  k may be an array."""
    k = np.asarray(k, dtype=float)
    return varpi ** 2 / (k ** 2 + varpi ** 2)


def scalar_O(fvals, N):
    """fvals[k] for k=0..2N (fvals[0] is IGNORED, matching the repo convention)."""
    u = u_vec(N)
    n = np.arange(N)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    M = np.outer(u, u) * (np.where(D == 0, 0.0, fvals[D]) + fvals[S])
    csum = np.concatenate([[0.0], np.cumsum(fvals[1:N])])  # csum[n] = sum_{k=1}^{n} f(k)
    M[n, n] -= 2.0 * u ** 2 * csum
    return M


def multiband_O(lam, N, w=None):
    """lam[k] is a (B,B) SYMMETRIC array Lt(k) (the symmetrized band matrix), k=0..2N.

    w: sqrt(DOS) weights, default all-ones.  Returns the symmetric (N*B, N*B) matrix
    in the ordering index = n*B + i.
    """
    lam = np.asarray(lam, dtype=float)
    B = lam.shape[1]
    if w is None:
        w = np.ones(B)
    w = np.asarray(w, dtype=float)
    u = u_vec(N)
    n = np.arange(N)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1

    # R_i(k) = (Lt(k) w)_i / w_i   -- the UNsymmetrized row sums driving Z
    R = np.einsum("kij,j->ki", lam, w) / w[None, :]  # (K, B)

    uu = np.outer(u, u)
    # pairing part: sum over band indices kept as blocks
    OO = np.zeros((N, B, N, B))
    pair = np.where(D == 0, 0.0, 1.0)[:, :, None, None] * lam[D] + lam[S]
    OO += uu[:, None, :, None] * np.transpose(pair, (0, 2, 1, 3))
    # Z part (band-diagonal), and the k=0 pairing entry is already included above? No:
    # lam[D] with D==0 was zeroed out, so re-add the k=0 pairing entry explicitly.
    for i in range(B):
        for j in range(B):
            OO[n, i, n, j] += u ** 2 * lam[0][i, j]
    csumR = np.concatenate([np.zeros((1, B)), np.cumsum(R[1:N], axis=0)])  # (N,B)
    for i in range(B):
        OO[n, i, n, i] -= u ** 2 * (R[0][i] + 2.0 * csumR[:, i])
    return OO.reshape(N * B, N * B)


def lam_factorized(Lt, varpi, N):
    """lam(k) = Lt * [[k]](varpi), k=0..2N."""
    ks = np.arange(2 * N + 2)
    f = brk(ks, varpi)
    return f[:, None, None] * np.asarray(Lt, dtype=float)[None, :, :]


def lam_from_measures(P, N):
    """P: list of (weight_matrix (B,B) elementwise>=0, varpi) pairs.
    lam(k) = sum_a W_a [[k]](varpi_a).  This is the general matrix-valued measure.
    """
    B = np.asarray(P[0][0]).shape[0]
    ks = np.arange(2 * N + 2)
    out = np.zeros((len(ks), B, B))
    for W, vp in P:
        out += brk(ks, vp)[:, None, None] * np.asarray(W, dtype=float)[None, :, :]
    return out


def lam1(M):
    return float(np.linalg.eigvalsh(M)[-1])


def lammin(M):
    return float(np.linalg.eigvalsh(M)[0])


def h_of(varpi, N):
    ks = np.arange(2 * N + 2)
    return lam1(scalar_O(brk(ks, varpi), N))


def g2(N):
    ks = np.arange(2 * N + 2).astype(float)
    f = np.zeros_like(ks)
    f[1:] = 1.0 / ks[1:] ** 2
    return lam1(scalar_O(f, N))
