"""
[numerical-evidence] Weighted-Laplacian / index-peeling lambda_k bound.

Theorem shape (analytic, proof is three lines -- see the note):

  A = B - Delta,  B >= 0 entrywise symmetric,  Delta diagonal >= 0.
  For ANY positive weight p, the generalized Laplacian
      L_p := diag( sum_{m != n} B_{nm} p_m / p_n ) - B_off
  satisfies L_p >= 0  (AM-GM on each unordered pair).  Hence
      A = -L_p + diag(delta_p),   delta_p,n = sum_m B_{nm} p_m/p_n - Delta_nn,
  and Weyl's inequality lambda_{k+1}(X+Y) <= lambda_1(X) + lambda_{k+1}(Y)
  with lambda_1(-L_p) = 0 gives

      lambda_{k+1}(A)  <=  (k+1)-th largest entry of delta_p .        (*)

  Taking p_m -> 0 for m in a finite set S (|S| = k) sends delta_p,m -> +infty
  for m in S and delta_p,n -> sum_{m not in S} B_{nm} p_m/p_n - Delta_nn for
  n not in S.  So (*) becomes the PEELED bound

      lambda_{k+1}(A)  <=  sup_{n not in S} [ sum_{m not in S} B_{nm} p_m/p_n
                                              - Delta_nn ] .          (**)

  We want k = 1, S = {0}:  lambda_2(A) <= sup_{n>=1} [ sum_{m>=1} B_{nm} - Delta_nn ].

Run:  OMP_NUM_THREADS=1 python _chk_laplacian_peel.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np


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


def peeled_delta(w, k, N, Mtail, p=None):
    """delta^{(S)}_n for n = k..N-1 with S = {0,...,k-1}, weight p on m >= k.

    p is a callable m -> p_m (defaults to 1).  Sum over m from k to Mtail-1.
    """
    x = w * w
    f = lambda t: x / (t * t + x)
    m = np.arange(k, Mtail)
    um = 1.0 / np.sqrt(2.0 * m + 1.0)
    pm = np.ones_like(um) if p is None else p(m.astype(float))
    out = np.empty(N - k)
    kk = np.arange(0, N, dtype=float)
    fk = f(kk)
    fk[0] = 0.0
    cum = np.cumsum(fk)
    for j, i in enumerate(range(k, N)):
        un = 1.0 / np.sqrt(2.0 * i + 1.0)
        pn = 1.0 if p is None else float(p(np.array([float(i)]))[0])
        d = np.abs(i - m).astype(float)
        s = (i + m + 1).astype(float)
        t = f(d)
        t[i - k] = 0.0                       # n != m; f(0) never evaluated
        t = t + f(s)
        rs = un * np.dot(t * um, pm) / pn
        out[j] = rs - 2.0 * cum[i] / (2.0 * i + 1.0)
    return out


def true_eigs(w, N):
    x = w * w
    f = lambda t: x / (t * t + x)
    ev = np.linalg.eigvalsh(O_of(f, N))
    return ev[-1], ev[-2]


if __name__ == "__main__":
    print("PEELED bound (**) with S={0}, weight p = 1   [numerical-evidence]")
    print()
    hdr = ("  varpi        h       lam2_true |   bound_k1    bound/h   "
           "bound/lam2 |  target .4h   PASS?")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for w in (0.30, 0.42, 0.60, 0.85, 1.00, 1.20, 1.70, 2.00, 2.40,
              3.40, 4.80, 8.50, 20.0):
        h, l2 = true_eigs(w, 1600)
        d = peeled_delta(w, 1, 400, 200000)
        b = d.max()
        tgt = 0.4 * h
        ok = "yes" if b <= tgt else "NO"
        sound = "" if b >= l2 else "   <<< BOUND BELOW TRUTH -- BUG"
        print(f"{w:7.2f} {h:9.5f} {l2:11.4e} | {b:11.6f} {b/h:10.4f} "
              f"{b/l2:10.2f} | {tgt:11.5f}   {ok}{sound}")
        rows.append((w, h, l2, b))
    print()
    print("argmax of the peeled profile (which n saturates), varpi=1.0:")
    d = peeled_delta(1.0, 1, 30, 200000)
    for j in range(10):
        print(f"   n={j+1:3d}  delta^(0)={d[j]: .6f}")
