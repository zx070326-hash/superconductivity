"""Case 2: Heisenberg XXZ chain, same Cholesky-deflation gap certificate.
H = sum_i [Sx_i Sx_{i+1} + Sy_i Sy_{i+1} + Delta Sz_i Sz_{i+1}]  (periodic, N sites)
Also time the N=8, N=10 runs to see where the mpmath dps=30 Cholesky wall is
(this is the Na^3 bottleneck flagged in the P4 breakpoint-4 audit)."""
import time
import numpy as np
import mpmath as mp
from probe_qc_gap_cert import op, cholesky_deflation_upper_bound, rayleigh_lower_bound

mp.mp.dps = 30

X = np.array([[0.0, 1.0], [1.0, 0.0]])
Y = np.array([[0.0, -1.0j], [1.0j, 0.0]]).real  # placeholder, build XX/YY via real Sx,Sy below
Z = np.array([[1.0, 0.0], [0.0, -1.0]])
I2 = np.eye(2)


def build_xxz(N, Delta, periodic=True):
    dim = 2 ** N
    H = np.zeros((dim, dim), dtype=complex)
    Sx = 0.5 * np.array([[0, 1], [1, 0]], dtype=complex)
    Sy = 0.5 * np.array([[0, -1j], [1j, 0]], dtype=complex)
    Sz = 0.5 * np.array([[1, 0], [0, -1]], dtype=complex)
    bonds = range(N) if periodic else range(N - 1)
    for i in bonds:
        j = (i + 1) % N
        H += op(Sx, i, N) @ op(Sx, j, N) + op(Sy, i, N) @ op(Sy, j, N) \
             + Delta * (op(Sz, i, N) @ op(Sz, j, N))
    assert np.max(np.abs(H.imag)) < 1e-10
    return H.real


def main():
    for N in (6, 8, 10):
        Delta = 0.5
        H = build_xxz(N, Delta)
        A = -H
        t0 = time.time()
        c, chol_ok, lam2_float, u, lam1_float = cholesky_deflation_upper_bound(A)
        dt = time.time() - t0
        rq = rayleigh_lower_bound(A, u)
        gap_true = lam1_float - lam2_float
        gap_cert_lb = float(rq) - float(c)
        print(f"N={N:2d} dim={2**N:4d}  Delta_xxz={Delta}  gap(float)={gap_true:9.6f}  "
              f"chol_ok={chol_ok}  gap_cert_lb={gap_cert_lb:9.6f}  slack={gap_true-gap_cert_lb:.2e}  "
              f"wall={dt:6.2f}s")


if __name__ == "__main__":
    main()
