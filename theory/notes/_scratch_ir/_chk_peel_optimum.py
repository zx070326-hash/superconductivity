"""[numerical-evidence] How much room is left in the peeled lambda_2 bound?

Chain of inequalities (each step analytic):

    lambda_2(A)  <=  lambda_1(A^{(0)})            [Cauchy interlacing, exact]
                 <=  sup_n [ (B~ p)_n / p_n - Delta_n ]   for any p > 0
                                                   [Schur/AM-GM, B~ >= 0]

where A^{(0)} = A with row/column 0 deleted, B~ = B likewise.  Collatz-Wielandt
says the second step is TIGHT at the optimal p (the Perron vector of A^{(0)}),
so lambda_1(A^{(0)}) is exactly the best any weight can do -- it is the target.

Columns:
  lam2      true second eigenvalue of A
  lam1_del  lambda_1(A^{(0)})   = the optimum of the weighted route
  P_1       the p = 1 row-sum bound (what _chk_laplacian_peel.py computed)
  P_pow     best p_n = (2n+1)^{-beta/2} over a beta grid  (explicit analytic p)

Run:  OMP_NUM_THREADS=1 python _chk_peel_optimum.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

MT = 40000


def O_of(phi, N):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    K = np.arange(0, 2 * N + 2, dtype=float)
    pk = phi(K); pk[0] = 0.0
    M = (pk[D] + pk[S]) * np.outer(u, u)
    M[np.diag_indices(N)] -= 2.0 * np.cumsum(pk[:N]) / (2.0 * n + 1.0)
    return M


def rowbound(w, beta, nrow=200):
    """sup_{1<=n<nrow} [ sum_{m>=1} B_nm p_m/p_n - Delta_nn ],  p_n=(2n+1)^{-b/2}."""
    x = w * w
    f = lambda t: x / (t * t + x)
    m = np.arange(1, MT, dtype=float)
    um = 1.0 / np.sqrt(2.0 * m + 1.0)
    pm = (2.0 * m + 1.0) ** (-beta / 2.0)
    kk = np.arange(0, nrow + 1, dtype=float); fk = f(kk); fk[0] = 0.0
    cum = np.cumsum(fk)
    best = -np.inf
    for lo in range(1, nrow, 25):
        hi = min(lo + 25, nrow)
        ii = np.arange(lo, hi, dtype=float)[:, None]
        t = f(np.abs(ii - m[None, :])) + f(ii + m[None, :] + 1.0)
        for r, i in enumerate(range(lo, hi)):
            t[r, i - 1] = f(float(2 * i + 1))
        un = 1.0 / np.sqrt(2.0 * ii[:, 0] + 1.0)
        pn = (2.0 * ii[:, 0] + 1.0) ** (-beta / 2.0)
        rs = un * (t * (um * pm)[None, :]).sum(axis=1) / pn
        idx = np.arange(lo, hi)
        best = max(best, (rs - 2.0 * cum[idx] / (2.0 * idx + 1.0)).max())
    return best


if __name__ == "__main__":
    BET = np.linspace(-1.0, 3.0, 17)
    hdr = ("  varpi        h       lam2     lam1_del  |      P_1     P_pow  "
           "beta* |  P_pow/h  lam1del/h")
    print("Room left in the peeled route   [numerical-evidence]")
    print(hdr); print("-" * len(hdr))
    for w in (0.30, 0.60, 1.00, 1.70, 2.40, 3.40, 4.80, 8.50, 20.0):
        x = w * w
        A = O_of(lambda k: x / (k * k + x), 1600)
        ev = np.linalg.eigvalsh(A)
        h, l2 = ev[-1], ev[-2]
        l1d = np.linalg.eigvalsh(A[1:, 1:])[-1]
        p1 = rowbound(w, 0.0)
        vals = [(rowbound(w, b), b) for b in BET]
        pp, bs = min(vals)
        flag = "" if (l1d >= l2 - 1e-12 and pp >= l1d - 1e-9) else "  <<< ORDER VIOLATED"
        print(f"{w:7.2f} {h:9.5f} {l2:10.3e} {l1d:11.3e} | {p1:9.5f} {pp:9.5f} "
              f"{bs:5.2f} | {pp/h:8.4f} {l1d/h:10.3e}{flag}")
