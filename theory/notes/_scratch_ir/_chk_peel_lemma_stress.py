"""Adversarial test of the LAPLACIAN-PEELING LEMMA itself (not of O[f]).

LEMMA [analytic].  Let A = B - Delta with B = B^T entrywise >= 0 and Delta
diagonal (Delta >= 0 not even needed).  For a finite index set S with |S| = k,

    lambda_{k+1}(A)  <=  sup_{n not in S} [ sum_{m not in S} B_{nm} p_m/p_n
                                            - Delta_nn ]         (**)

for every positive weight p on the complement of S.  (k = 0, p = 1 is the
classical row-sum/Gershgorin-type bound; k = 1 is what we need for lambda_2.)

If (**) is false, random search finds a violation fast.  This script hunts for
one over random sign patterns, sizes, weights, and peel sets, and reports the
tightest margin found.  A single negative margin refutes the lemma.

Run:  OMP_NUM_THREADS=1 python _chk_peel_lemma_stress.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

rng = np.random.default_rng(20260917)


def bound(B, Delta, S, p):
    keep = np.array([i for i in range(len(Delta)) if i not in S])
    Bk = B[np.ix_(keep, keep)]
    pk = p[keep]
    rs = (Bk * pk[None, :]).sum(axis=1) / pk
    return (rs - Delta[keep]).max()


worst = np.inf
worst_cfg = None
for trial in range(40000):
    N = int(rng.integers(2, 9))
    k = int(rng.integers(0, min(N - 1, 4)))
    # entrywise-nonnegative symmetric B, various sparsity / scale regimes
    style = trial % 4
    if style == 0:
        M = rng.random((N, N))
    elif style == 1:
        M = rng.random((N, N)) * (rng.random((N, N)) < 0.35)
    elif style == 2:
        M = rng.exponential(1.0, (N, N)) ** 3
    else:
        M = rng.random((N, N)) * 10.0 ** rng.integers(-4, 4, (N, N))
    B = M + M.T
    Delta = rng.normal(0, 2, N) if trial % 7 == 0 else rng.random(N) * 5
    A = B - np.diag(Delta)
    ev = np.sort(np.linalg.eigvalsh(A))[::-1]
    S = set(rng.choice(N, size=k, replace=False).tolist())
    p = 10.0 ** rng.normal(0, 1.2, N) if trial % 3 else np.ones(N)
    b = bound(B, Delta, S, p)
    margin = b - ev[k]               # (**) says margin >= 0
    if margin < worst:
        worst = margin
        worst_cfg = (N, k, sorted(S), ev[k], b)

print("Lemma (**) stress test: margin = bound - lambda_{k+1};  "
      "a NEGATIVE margin refutes it.")
print(f"  trials              : 40000")
print(f"  worst margin found  : {worst:.6e}")
print(f"  at (N,k,S,lam,bound): {worst_cfg}")
print()
print("VERDICT:", "no counterexample (lemma survives)" if worst >= -1e-12
      else "*** REFUTED ***")

# --- tightness: is the bound ever attained?  (k=0, B=J, Delta=0 -> equality) ---
for N in (3, 5, 8):
    B = np.ones((N, N)); Delta = np.zeros(N)
    ev = np.sort(np.linalg.eigvalsh(B))[::-1]
    print(f"  tightness check N={N}: lambda_1={ev[0]:.6f}  "
          f"bound(k=0,p=1)={bound(B, Delta, set(), np.ones(N)):.6f}")
