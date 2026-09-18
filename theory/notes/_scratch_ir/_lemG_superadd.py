"""Superadditivity scheme:
  f = sum_K a_K 1_{k<=K},  a_K >= 0   (Abel, f nonincreasing >=0)
  -lmin is SUBadditive:  -lmin(O[f]) <= sum_K a_K r_K,  r_K := -lmin(O[1_{k<=K}])
  ell(f) := <v,O[f]v>/||v||^2 is LINEAR in f and <= lambda_1(O[f]).
  => if r_K <= c * q_K(v) for all K, then |lmin(O[f])| <= c lambda_1(O[f]) for ALL f.
Here we scan v to minimize sup_K r_K/q_K(v).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _lemG_reduce import build, ev, step


def AK_block(K, p):
    """exact top-left p x p block of O[1_{k<=K}]."""
    n = np.arange(p)
    u = 1.0 / np.sqrt(2 * n + 1.0)
    I, J = np.meshgrid(n, n, indexing="ij")
    ad = np.abs(I - J)
    off = ((ad <= K) & (ad >= 1)).astype(float) + (I + J + 1 <= K).astype(float)
    A = off * np.outer(u, u)
    d = (((2 * n + 1) <= K).astype(float) - 2.0 * np.minimum(n, K)) / (2 * n + 1.0)
    np.fill_diagonal(A, d)
    return A


def qK(v, K):
    p = len(v)
    A = AK_block(K, p)
    return float(v @ A @ v) / float(v @ v)


def rK_table(Ks, mult=30, base=800):
    out = {}
    for K in Ks:
        N = min(mult * K + base, 4000)
        B = build(step(K), N)
        l1, lmin = ev(B["O"])
        out[K] = (-lmin, l1)
    return out


if __name__ == "__main__":
    Ks = [1, 2, 3, 4, 5, 6, 8, 10, 14, 20, 30, 45, 70, 110]
    tab = rK_table(Ks)
    for K in Ks:
        print(f"K={K:4d}  r_K={tab[K][0]:.6f}  l1={tab[K][1]:.6f}  rho={tab[K][0]/tab[K][1]:.6f}")

    print()
    cands = {
        "e1": np.array([1.0]),
        "(1,1/sqrt3)": np.array([1.0, 1 / np.sqrt(3)]),
        "topvec(A_1) p=8": None,
    }
    w, V = np.linalg.eigh(AK_block(1, 12))
    tv = np.abs(V[:, -1])
    cands["topvec(A_1) p=12"] = tv
    cands.pop("topvec(A_1) p=8")
    for name, v in cands.items():
        vals = [tab[K][0] / qK(v, K) for K in Ks]
        print(f"{name:22s} sup_K r/q = {max(vals):.6f} at K={Ks[int(np.argmax(vals))]}   "
              f"q_1={qK(v,1):.5f} q_110={qK(v,110):.5f}")
        print("     ", " ".join(f"{K}:{val:.4f}" for K, val in zip(Ks, vals)))
