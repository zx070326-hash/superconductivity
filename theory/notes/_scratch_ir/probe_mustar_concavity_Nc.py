"""Independent from-zero recheck (float64, numerical-evidence) of the claimed
convexity window "varpi <~ 0.30-0.33 for mu_tilde <~ 0.05" in
theory/notes/pairing_duality_theorem_2026-09-16.md Sec 5.2-5.3.

Three things are being tested, each an alleged defect in the source:

  (A) spectral_efficiency.py:358 builds  u_n = (2n+1)^{-1/2}  over the WHOLE
      basis (Nuse=200), i.e. the Coulomb window is the numerical truncation,
      not the physical N_c = omega_c/(2 pi T) ~ 6-45.  Redo with a mask.

  (B) the source's "2-atom vs Einstein" check uses varpi = 1 and sqrt(17),
      which is a GLOBAL comparison and says nothing about local concavity.
      Redo as a merging test: two atoms at x(1 -+ d) with d -> 0, where the
      Jensen gap is ~ (1/2) H''(x) x^2 d^2 and its SIGN reads off H''.

  (C) H''(x) was taken from a 5-point stencil on float64 eigenvalues with
      dx = max(5e-3 x, 1e-3).  At small x the signal may be below the
      cancellation floor.  Redo at several dx and check stability.

Operator (rebuilt here from scratch, not imported):
    f(k)   = x/(k^2 + x),  x = varpi^2,  f(0) = 1
    u_n    = (2n+1)^{-1/2}
    H_nm   = (f(|n-m|) + f(n+m+1)) u_n u_m            (n != m)
    H_nn   = -(1/(2n+1)) sum_{k=1}^{n} 2 f(k) + f(2n+1) u_n^2
    M      = H - 2 mu_tilde |u_c><u_c|,  u_c = u * 1[n < N_c]
The factor 2 on the Coulomb term is the +/- frequency image (see note).
"""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigvalsh


def build_op(weights_x, N):
    """weights_x: [(w, x)] with x = varpi^2.  Returns the N x N operator."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    ks = np.arange(2 * N + 1, dtype=float)
    f = np.zeros(2 * N + 1)
    for w, x in weights_x:
        f += w * x / (ks ** 2 + x)
    f[0] = sum(w for w, _ in weights_x)
    A = (f[np.abs(n[:, None] - n[None, :])] + f[n[:, None] + n[None, :] + 1]) \
        * (u[:, None] * u[None, :])
    L1 = np.zeros(N)
    L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * n + 1] * u ** 2)
    return A, u


def kappa(weights_x, mu_t, Nc, N):
    A, u = build_op(weights_x, N)
    if mu_t:
        uc = np.where(np.arange(N) < (Nc if Nc else N), u, 0.0)
        A = A - 2.0 * mu_t * np.outer(uc, uc)
    return float(eigvalsh(A)[-1])


def H_pp(x, mu_t, Nc, N, dx):
    v = [kappa([(1.0, max(x + i * dx, 1e-9))], mu_t, Nc, N) for i in range(-2, 3)]
    return (-v[0] + 16 * v[1] - 30 * v[2] + 16 * v[3] - v[4]) / (12 * dx ** 2)


if __name__ == "__main__":
    N = 200

    # --- sanity: reproduce build_H_fast2 at mu=0 (same object?) ---
    for vp in (0.3, 1.0, 3.0):
        print(f"[sanity] varpi={vp}: h_mine={kappa([(1.0, vp**2)], 0, None, N):.12f}")

    print("\n(C) step-size stability of H''(x) at mu_tilde=0.03, Nc=None (source setup)")
    print(f"{'x':>7} " + " ".join(f"{'dx='+s:>13}" for s in ['x/200', 'x/50', 'x/10']))
    for x in (0.01, 0.05, 0.09, 0.1089, 0.15, 0.25, 0.5, 1.0):
        row = [H_pp(x, 0.03, None, N, x / d) for d in (200.0, 50.0, 10.0)]
        print(f"{x:>7.4f} " + " ".join(f"{v:>13.4e}" for v in row))

    print("\n(A) H''(x) with the PHYSICAL Coulomb window N_c (dx = x/50)")
    print(f"{'x':>7} {'varpi':>7} " + " ".join(f"{'Nc='+str(c):>12}"
          for c in (6, 15, 45, 200)))
    for x in (0.01, 0.04, 0.09, 0.1089, 0.16, 0.25, 0.5, 1.0, 4.0):
        row = [H_pp(x, 0.05, c, N, x / 50.0) for c in (6, 15, 45, 200)]
        print(f"{x:>7.4f} {math.sqrt(x):>7.4f} " + " ".join(f"{v:>12.4e}" for v in row))

    print("\n(B) merging 2-atom test: x_+- = x(1 -+ d), weights 1/2, <x> = x exactly")
    print("    gap := kappa(2-atom) - kappa(Einstein);  gap>0 => Einstein NOT extremal")
    for (mu_t, Nc) in ((0.05, 6), (0.05, 45), (0.05, None), (0.13, 6), (0.0, None)):
        print(f"  mu_tilde={mu_t}  Nc={Nc}")
        for x in (0.09, 0.1089, 0.25, 1.0, 9.0):
            ke = kappa([(1.0, x)], mu_t, Nc, N)
            out = []
            for d in (0.4, 0.2, 0.1, 0.05):
                k2 = kappa([(0.5, x * (1 - d)), (0.5, x * (1 + d))], mu_t, Nc, N)
                out.append((k2 - ke) / d ** 2)
            print(f"    x={x:>7.4f} kE={ke:.8f}  gap/d^2 @d=.4,.2,.1,.05: "
                  + " ".join(f"{v:>11.3e}" for v in out))
