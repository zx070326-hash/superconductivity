"""
Probe: reuse the P4 Cholesky-deflation certificate machine
(theory/derivations/p4_interchange_uniformity.py :: alpha_certificate_H)
to produce a SOUND lower bound on the spectral gap of a small
transverse-field Ising model (TFIM) Hamiltonian, as a feasibility
check for "certified spectral gap bounds for quantum algorithms".

Convention: H(s) = -(1-s) sum_i X_i - s sum_i Z_i Z_{i+1 mod N}   (periodic TFIM, N qubits)
Ground-state gap Delta(s) = E_1(s) - E_0(s).

Let A := -H(s).  Then lambda_1(A) = -E_0 (top eigenvalue), lambda_2(A) = -E_1
(second-largest).  Delta(s) = lambda_1(A) - lambda_2(A).

Bound directions (must each be checked -- one-sided-bound discipline):
  - lambda_1(A) >= Rayleigh(A, v)      for ANY real vector v  (variational; sound
      in exact arithmetic, numerical-evidence here since float64 quotient, not
      interval-hardened -- flagged explicitly, this is a probe not a certificate)
  - lambda_2(A) <= c   via the SAME Cholesky-deflation test as alpha_certificate_H:
      deflate the float64 top eigenvector u of A, build P A P (P = I - u u^T),
      certify c*I - PAP is SPD via mpmath dps=30 Cholesky (rank-1 exact deflation
      here, no tail/off-diagonal Schur terms needed since A is the WHOLE matrix,
      not a truncation -- so this step *is* certificate-level, modulo the
      residual-of-inexact-eigenvector argument copied verbatim from P4-C).
  => Delta(s) >= Rayleigh(A,v) - c   is a sound-modulo-Rayleigh-hardening lower
     bound on the true gap.
"""
import numpy as np
import mpmath as mp

mp.mp.dps = 30

X = np.array([[0.0, 1.0], [1.0, 0.0]])
Z = np.array([[1.0, 0.0], [0.0, -1.0]])
I2 = np.eye(2)


def op(pauli, site, N):
    mats = [I2] * N
    mats[site] = pauli
    M = mats[0]
    for m in mats[1:]:
        M = np.kron(M, m)
    return M


def build_H(N, s, periodic=True):
    dim = 2 ** N
    H = np.zeros((dim, dim))
    for i in range(N):
        H -= (1 - s) * op(X, i, N)
    bonds = range(N) if periodic else range(N - 1)
    for i in bonds:
        j = (i + 1) % N
        ZZ = op(Z, i, N) @ op(Z, j, N)
        H -= s * ZZ
    return H


def cholesky_deflation_upper_bound(A_np, Na_check_full=True):
    """Certified upper bound on lambda_2(A) via rank-1 Cholesky deflation.
    A_np: full symmetric matrix (small enough for exact eigh, dim<=few hundred).
    Returns (c, chol_ok, lam2_float_reference)."""
    dim = A_np.shape[0]
    w, V = np.linalg.eigh(A_np)
    u = V[:, -1]
    lam2_float = w[-2]

    Au_f = A_np @ u
    uAu_f = float(u @ Au_f)
    PAP_f = A_np - np.outer(u, Au_f) - np.outer(Au_f, u) + uAu_f * np.outer(u, u)
    c = mp.mpf(float(np.linalg.eigvalsh(PAP_f)[-1])) + mp.mpf("1e-9")

    Am = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            Am[i, j] = mp.mpf(float(A_np[i, j]))
    um = mp.matrix([mp.mpf(float(x)) for x in u])
    Au = Am * um
    uAu = sum(um[i] * Au[i] for i in range(dim))
    M = mp.matrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            pap = Am[i, j] - um[i] * Au[j] - Au[i] * um[j] + um[i] * um[j] * uAu
            M[i, j] = (c if i == j else mp.mpf(0)) - pap
    try:
        mp.cholesky(M)
        chol_ok = True
    except Exception:
        chol_ok = False
    return c, chol_ok, lam2_float, u, w[-1]


def rayleigh_lower_bound(A_np, v):
    """Rayleigh quotient (float64) -- numerical-evidence direction check only."""
    Av = A_np @ v
    return float(v @ Av) / float(v @ v)


def main():
    for N in (6, 8):
        print(f"\n=== TFIM, N={N} qubits, dim={2**N} ===")
        for s in (0.3, 0.5, 0.5, 0.7):
            H = build_H(N, s, periodic=True)
            A = -H
            c, chol_ok, lam2_float, u, lam1_float = cholesky_deflation_upper_bound(A)
            rq = rayleigh_lower_bound(A, u)
            gap_true = lam1_float - lam2_float          # plain float64 reference
            gap_cert_lb = float(rq) - float(c)           # our sound-direction bound
            print(f"  s={s:4.2f}  E0(float)={-lam1_float:9.5f}  gap(float)={gap_true:9.6f}  "
                  f"chol_ok={chol_ok}  c>=lam2>={float(c):9.6f}  "
                  f"gap_cert_lb={gap_cert_lb:9.6f}  slack={gap_true-gap_cert_lb:.2e}")


if __name__ == "__main__":
    main()
