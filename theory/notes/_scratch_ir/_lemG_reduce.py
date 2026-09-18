"""Lemma (G) navigation: build O[f] from scratch, test reductions.

O[f]_{nn} = -2/(2n+1) * sum_{k=1}^n f(k) + f(2n+1) u_n^2
O[f]_{nm} = [f(|n-m|) + f(n+m+1)] u_n u_m   (n != m)
u_n = (2n+1)^{-1/2}
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np


def build(f, N):
    """f: callable on k>=1 (and we need f(0)? no: |n-m|>=1 off-diag, n+m+1>=1)."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    fk = np.array([f(k) for k in range(1, 2 * N + 2)])  # fk[k-1] = f(k)

    def F(k):  # k >= 1 array
        return fk[np.asarray(k) - 1]

    S = np.concatenate([[0.0], np.cumsum(fk[:N - 1])])  # S[n] = sum_{k<=n} f(k)
    I, J = np.meshgrid(n, n, indexing="ij")
    absd = np.abs(I - J)
    Moff = np.where(absd == 0, 0.0, F(np.maximum(absd, 1))) * np.outer(u, u)
    Mall = Moff + F(I + J + 1) * np.outer(u, u)   # L2 (offdiag) + L3 (all)
    L1 = 2.0 * S / (2 * n + 1.0)
    d = -L1 + F(2 * n + 1) * u ** 2
    O = Mall.copy()
    np.fill_diagonal(O, d)
    Moff_full = Mall.copy()
    np.fill_diagonal(Moff_full, F(2 * n + 1) * u ** 2)  # M = L2+L3 elementwise >=0
    # M_offdiag only (zero diagonal)
    Mo = Mall.copy()
    np.fill_diagonal(Mo, 0.0)
    return dict(O=O, d=d, Mo=Mo, M=Moff_full, L1=np.diag(L1), Mdiag=F(2 * n + 1) * u ** 2, u=u)


def ev(A):
    w = np.linalg.eigvalsh(A)
    return w[-1], w[0]


def report(name, f, N=600):
    B = build(f, N)
    O, d, Mo = B["O"], B["d"], B["Mo"]
    l1, lmin = ev(O)
    Ot = Mo + np.diag(-d)                 # tilde O : flip diagonal sign
    lt = ev(Ot)[0]
    Q = Mo + np.diag(np.maximum(-d, 0) + B["Mdiag"] * 0) + B["L1"] - np.diag(B["Mdiag"]) * 0
    # Q := M_off + L1 (drop -Mdiag which only lowers)  -- elementwise nonneg
    Q = Mo + B["L1"]
    np.fill_diagonal(Q, np.diag(Q) + 0.0)
    lq = ev(Q)[0]
    # Schur / Collatz-Wielandt with weight u on tilde O
    u = B["u"]
    cw = np.max((Ot @ u) / u)
    print(f"{name:28s} l1={l1:9.6f} |lmin|={abs(lmin):9.6f} "
          f"rho={abs(lmin)/l1:8.6f} rho#={lt/l1:8.6f} rhoQ={lq/l1:8.6f} cw/l1={cw/l1:8.4f}")
    return dict(l1=l1, lmin=lmin, lt=lt, lq=lq, cw=cw)


def step(K):
    return lambda k: 1.0 if k <= K else 0.0


if __name__ == "__main__":
    for K in [1, 2, 3, 5, 10, 100]:
        report(f"step 1_[k<={K}]", step(K))
    report("1/k^2", lambda k: 1.0 / k ** 2)
    report("1/k^3", lambda k: 1.0 / k ** 3)
    for w in [0.1, 1.0, 7.0, 1000.0]:
        report(f"Einstein w={w}", lambda k, w=w: w * w / (k * k + w * w))
    report("1/k^1.2", lambda k: 1.0 / k ** 1.2)
    report("exp(-k)", lambda k: np.exp(-float(k)))
