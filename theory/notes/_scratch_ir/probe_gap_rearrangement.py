"""Scratch: general rearrangement identity for ANY kernel f>=0 (not just 1/k^2):
  <v,O[f]v> = sum_n p_n^2 T_f(n) + <p, Hank_f p> - 1/2 D_f[p],  p_n = u_n v_n,
  T_f(n) = sum_{k>n} f(k),  Hank_f = f(n+m+1),  D_f[p] = sum_{n!=m} f(|n-m|)(p_n-p_m)^2 >= 0.
In v-coordinates:  H(w) <= B_w := diag(T(n)/(2n+1)) + L3[f]   (Loewner).
Test: (a) identity, (b) is B_w PSD / what is lam_min(L3)?, (c) trace bounds on lam_2 via B_w."""
import os, sys
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations")
sys.path.insert(0, os.path.abspath(ROOT))
from p4_interchange_uniformity import H_mat  # noqa: E402

N = 1500
n = np.arange(N)
u = 1.0 / np.sqrt(2 * n + 1.0)
print(f"{'w':>5} {'h':>8} {'lam2':>9} {'target.6h':>9} | {'minD':>9} {'lam_min(B)':>10} {'trB-l1B':>9} {'sqrt(trB2-l1^2)':>15} {'best/true':>9}")
for w in (0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0):
    f = lambda k: w * w / (k * k + w * w)
    K = 200000
    ks = np.arange(1, K + 1)
    fk = f(ks)
    tail = np.cumsum(fk[::-1])[::-1]              # tail[i] = sum_{k>=i+1} f(k)
    T = np.concatenate((tail[1:N + 1],))           # T(n) = sum_{k>n} f(k) = tail[n]  (index n -> tail[n])
    T = tail[n]                                    # tail[n] = sum_{k>=n+1}
    T = T + w * w * (1.0 / (K + 0.5))              # add ~integral tail beyond K (upper estimate, tiny)
    S = n[:, None] + n[None, :] + 1
    L3 = f(S) * np.outer(u, u)
    B = np.diag(T * u * u) + L3
    H = H_mat(N, w, 0)
    Dhalf = B - H
    E_H = np.linalg.eigvalsh(H)
    h, l2 = E_H[-1], E_H[-2]
    E_B = np.linalg.eigvalsh(B)
    E_D = np.linalg.eigvalsh(Dhalf)
    trB = np.trace(B)
    l1B = E_B[-1]
    b1 = trB - l1B
    b2 = np.sqrt(max(np.sum(E_B ** 2) - l1B ** 2, 0.0))   # valid if B >= 0 (else use |eigs|)
    b2abs = np.sqrt(max(np.trace(B @ B) - l1B ** 2, 0.0))
    best = min(b1 if E_B[0] >= -1e-12 else np.inf, b2abs)
    print(f"{w:>5.1f} {h:>8.4f} {l2:>9.2e} {0.6*h:>9.4f} | {E_D[0]:>9.2e} {E_B[0]:>10.2e} {b1:>9.4f} {b2abs:>15.4f} {best/l2:>9.1f}")
print("\nminD = lam_min(B-H) should be >= 0 (identity: B-H = 1/2 D in v-coords).")
print("lam_min(B): if >=0 the first-order trace bound lam2 <= trB - lam1(B) is valid.")
