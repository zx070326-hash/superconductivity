"""
[numerical-evidence] Laplacian-defect route to an analytic lambda_2 bound.

Structure claim (analytic, to be tested numerically here):

    A = O[f]  =  B - Delta ,
        B_{nm} = u_n u_m [ f(|n-m|) 1_{n!=m} + f(n+m+1) ]   >= 0  entrywise,
        Delta  = diag( 2 u_n^2 sum_{k=1}^n f(k) )           >= 0.

Split B = B_off + diag(B_nn) and introduce the weighted graph Laplacian
    L := diag( rowsum(B_off) ) - B_off   ,   L >= 0  (B_off >= 0 entrywise).
Then
    A = -L + diag(delta),       delta_n := sum_m B_{nm} - Delta_nn .

Since L >= 0 we have lambda_1(-L) = 0 and lambda_k(-L) <= 0 for all k.  Weyl:
    lambda_1(A) <= lambda_1(-L) + lambda_1(diag delta) =  max_n delta_n
    lambda_2(A) <= lambda_1(-L) + lambda_2(diag delta) =  2nd largest delta_n
(the second uses lambda_2(X+Y) <= lambda_1(X) + lambda_2(Y), valid for
 symmetric X,Y.)

More generally, if delta_n <= c for every n outside a set S with |S| = k, then
    lambda_{k+1}(A) <= c.

This script measures the profile of delta_n and asks whether the 2nd-largest
(or the "all but one" sup) is small -- i.e. whether this gives a usable
analytic lambda_2 bound for varpi in (0,2].

Run:  OMP_NUM_THREADS=1 python _chk_laplacian_defect.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

THETA_STAR = 0.6349287917881079


def O_of(phi, N):
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    D = np.abs(n[:, None] - n[None, :])
    S = n[:, None] + n[None, :] + 1
    K = np.arange(0, 2 * N + 2, dtype=float)
    pk = phi(K)
    pk[0] = 0.0
    M = (pk[D] + pk[S]) * np.outer(u, u)
    cum = np.cumsum(pk[:N])
    M[np.diag_indices(N)] -= 2.0 * cum / (2.0 * n + 1.0)
    return M


def delta_profile(w, N, Mtail):
    """delta_n for n = 0..N-1, with the m-sum carried out to m = Mtail-1.

    Mtail >> N so the row sum is essentially complete (tail ~ m^{-5/2}).
    """
    x = w * w
    f = lambda k: x / (k * k + x)
    n = np.arange(N)
    m = np.arange(Mtail)
    un = 1.0 / np.sqrt(2.0 * n + 1.0)
    um = 1.0 / np.sqrt(2.0 * m + 1.0)
    # row sums of B
    rs = np.empty(N)
    for i in range(N):
        d = np.abs(i - m).astype(float)
        s = (i + m + 1).astype(float)
        t = f(d)
        t[i] = 0.0                      # f(0) never evaluated, n != m
        t = t + f(s)
        rs[i] = un[i] * np.dot(t, um)
    # Delta_nn = 2 u_n^2 sum_{k=1}^n f(k)
    kk = np.arange(0, N, dtype=float)
    fk = f(kk)
    fk[0] = 0.0
    cum = np.cumsum(fk)
    Dnn = 2.0 * cum / (2.0 * n + 1.0)
    return rs - Dnn


def true_eigs(w, N):
    x = w * w
    f = lambda k: x / (k * k + x)
    A = O_of(f, N)
    ev = np.linalg.eigvalsh(A)
    return ev[-1], ev[-2]


if __name__ == "__main__":
    print("delta_n = (full row sum of B)_n - Delta_nn   [numerical-evidence]")
    print()
    hdr = ("  varpi        h       lam2   |  max d    2nd d    3rd d   "
           "sup_{n>=1} d |  2nd/h    2nd/lam2")
    print(hdr)
    print("-" * len(hdr))
    for w in (0.30, 0.42, 0.60, 0.85, 1.00, 1.20, 1.70, 2.00, 3.40, 8.50, 20.0):
        N = 400
        h, l2 = true_eigs(w, 1600)
        d = delta_profile(w, N, 200000)
        srt = np.sort(d)[::-1]
        print(f"{w:7.2f} {h:9.5f} {l2:10.3e} | {srt[0]:8.5f} {srt[1]:8.5f} "
              f"{srt[2]:8.5f} {d[1:].max():12.5f} | {srt[1]/h:7.4f} "
              f"{srt[1]/l2:10.2e}")
    print()
    print("argmax location and the shape of delta_n (varpi = 1.0):")
    d = delta_profile(1.0, 40, 200000)
    for i in range(12):
        print(f"   n={i:3d}  delta={d[i]: .6f}")
    print("   argmax =", int(np.argmax(delta_profile(1.0, 400, 200000))))
