"""[numerical-evidence] From-zero check of the M1 identity + Lemma-Lambda applied
to (R - D) rather than to H itself.

Repo, P4_einstein_extremal.md 1.F.2 / 1.F.3(c):
    H(w) = 2 tau tau^T + R - D,
    tau_n = p_n G(s_n),  p_n = (2n+1)^{-1/2},  s_n = (n+1/2)/w,  G(y)=1/(1+y^2)
    R_{nm} = p_n p_m Rk(s_n,s_m) >= 0 entrywise
    Rk(s,sg) = 2 s^2 sg^2 (5+s^2+sg^2) / [(1+(s-sg)^2)(1+(s+sg)^2)(1+s^2)(1+sg^2)]
    D = diag(w_n + l_G(n)),  w_n=(2n+1)^{-1},  l_G(n)=(2/(2n+1)) sum_{k=1}^n G(k/w)

and the already-written min-max  lam2(H) <= lam1(R - D) <= ||R|| <= Theta*.

NEW here: R >= 0 entrywise and D diagonal, so Lemma Lambda applies to (R - D):

    lam_{k+1}(R - D) <= sup_{n not in S} [ sum_{m not in S} R_{nm} q_m/q_n - D_nn ]

With q_n = p_n (2n+1)^{-b/2} the summand is

    w_m * Rk(s_n,s_m) * ((2n+1)/(2m+1))^{b/2}

so b=0, S={} is EXACTLY the repo's Schur row sum (bounded by Theta*), except that
it keeps the -D_nn that the repo discards.  b>0 and peeling are free extras.

Tail pad (upward, valid for b>=0 and m > n):
    Rk(s,sg) <= G(s-sg)+G(s+sg) <= 2/(1+(sg-s)^2)
    => sum_{m>=M} w_m Rk (..)^{b/2} <= 2 w^2 / ((2M+1)(M-n-1))

Run: OMP_NUM_THREADS=1 python _chk_Rop_identity_peel.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

MT, NROW = 40000, 200
BETAS = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
PEELS = (0, 1, 2, 3)


def O(phi, N):
    """repo operator convention; phi(0) never used, Hankel keeps the diagonal."""
    n = np.arange(N); I = n[:, None]; J = n[None, :]
    pv = phi(np.arange(0, 2 * N + 2, dtype=float))
    T = pv[np.abs(I - J)].copy(); np.fill_diagonal(T, 0.0)
    u = 1 / np.sqrt(2 * n + 1.)
    M = (T + pv[I + J + 1]) * u[:, None] * u[None, :]
    M[n, n] -= 2.0 * np.concatenate(([0.], np.cumsum(pv[1:N]))) / (2 * n + 1.)
    return M


def G(y):
    return 1.0 / (1.0 + y * y)


def Rk(s, sg):
    num = 2.0 * s * s * sg * sg * (5.0 + s * s + sg * sg)
    den = ((1.0 + (s - sg) ** 2) * (1.0 + (s + sg) ** 2)
           * (1.0 + s * s) * (1.0 + sg * sg))
    return num / den


def pieces(w, N):
    n = np.arange(N, dtype=float)
    s = (n + 0.5) / w
    p = 1.0 / np.sqrt(2.0 * n + 1.0)
    tau = p * G(s)
    Rm = p[:, None] * p[None, :] * Rk(s[:, None], s[None, :])
    lG = 2.0 * np.concatenate(([0.0], np.cumsum(G(np.arange(1, N) / w)))) / (2 * n + 1.0)
    D = 1.0 / (2.0 * n + 1.0) + lG
    return tau, Rm, D


def peel_bound(w, beta, npeel):
    """Lemma Lambda on (R - D), S = {0..npeel-1}, q_n = p_n (2n+1)^{-beta/2}."""
    x = w * w
    lo0 = npeel
    m = np.arange(lo0, MT, dtype=float)
    sm = (m + 0.5) / w
    wm = 1.0 / (2.0 * m + 1.0)
    fac_m = wm * (2.0 * m + 1.0) ** (-beta / 2.0)
    # cumG[n] = sum_{k=1}^{n} G(k/w)   (cumG[0] = 0, matching the repo's l_G)
    cumG = np.concatenate(([0.0], np.cumsum(G(np.arange(1, NROW + 1) / w))))
    best = -np.inf
    for a in range(lo0, NROW, 20):
        b = min(a + 20, NROW)
        nn = np.arange(a, b, dtype=float)
        sn = (nn + 0.5) / w
        row = Rk(sn[:, None], sm[None, :]) @ fac_m
        row *= (2.0 * nn + 1.0) ** (beta / 2.0)
        pad = 2.0 * x / ((2.0 * MT + 1.0) * (MT - nn - 1.0))
        Dnn = (1.0 + 2.0 * cumG[a:b]) / (2.0 * nn + 1.0)
        best = max(best, float((row + pad - Dnn).max()))
    return best


if __name__ == "__main__":
    print("== (1) identity  H = 2 tau tau^T + R - D  (max abs residual) ==")
    for w in (0.3, 1.0, 3.4, 20.0):
        N = 400
        H = O(lambda k: w * w / (k * k + w * w), N)
        tau, Rm, D = pieces(w, N)
        res = np.abs(H - (2.0 * np.outer(tau, tau) + Rm - np.diag(D))).max()
        print(f"   varpi={w:6.2f}   max|H - (2tt+R-D)| = {res:.3e}")

    print("\n== (2) lam2(H) vs lam1(R-D) vs lam1(R) vs Theta*=0.63493 ==")
    print(f"{'varpi':>7} {'h':>9} {'lam2(H)':>11} {'lam1(R-D)':>11} {'lam1(R)':>10}"
          f" {'Schur b=0':>10}")
    for w in (0.3, 0.6, 1.0, 2.0, 3.4, 8.5, 20.0):
        N = 1200
        H = O(lambda k: w * w / (k * k + w * w), N)
        tau, Rm, D = pieces(w, N)
        ev = np.linalg.eigvalsh(H)
        l1RD = np.linalg.eigvalsh(Rm - np.diag(D))[-1]
        l1R = np.linalg.eigvalsh(Rm)[-1]
        schur = peel_bound(w, 0.0, 0)
        print(f"{w:7.2f} {ev[-1]:9.5f} {ev[-2]:11.4e} {l1RD:11.4e} {l1R:10.5f}"
              f" {schur:10.5f}")

    print("\n== (3) Lemma Lambda on (R-D): best over beta x peel ==")
    print(f"{'varpi':>7} {'lam1(R-D)':>11} {'best bnd':>10} {'beta*':>6} {'peel':>5}"
          f" {'ratio':>8} {'Theta*':>8}")
    for w in (0.3, 0.42, 0.6, 0.85, 1.2, 1.7, 2.4, 3.4, 4.8, 6.8, 9.6, 13.6, 20.0):
        N = 1200
        tau, Rm, D = pieces(w, N)
        l1RD = np.linalg.eigvalsh(Rm - np.diag(D))[-1]
        cand = [(peel_bound(w, b, k), b, k) for b in BETAS for k in PEELS]
        v, b, k = min(cand)
        print(f"{w:7.2f} {l1RD:11.4e} {v:10.5f} {b:6.2f} {k:5d}"
              f" {v / l1RD:8.1f} {0.63493:8.5f}")
