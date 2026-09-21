"""The negative control v3 was asked for and did not implement, run by the main session.

A cross-check that cannot fail is worth nothing (the discarded production script had
pt and iv wrong in the SAME way, so its self-check was all green). So:

  (1) POSITIVE: v3's interval L1/A must ENCLOSE the externally computed correct values
      (_main_indep_baseline.py, written from the derivation text).
  (2) NEGATIVE (poison): v3's interval L1/A must NOT enclose the values produced by the
      OLD off-by-one operator (L1 summed only to k=n-1).  If it does enclose them, the
      comparison has no discriminating power and its PASSes mean nothing.

N=40, T=1.0 meV, mu*=0 (so v3 returns A itself).  [numerical-evidence] + interval.
"""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from mpmath import mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _main_indep_baseline import read_a2f, P_measure, F_T, O_of_f, A2F
from h3s_deficit_cert_v3 import build_M_iv

N, T, DPS = 40, 1.0, 50


def O_poisoned(F, N):
    """The OLD buggy operator: L1 summed only to k = n-1 (off-by-one).

    Written as an explicit sum so the upper limit is unmistakable; a cumulative
    loop is exactly how the off-by-one hides (and how my first attempt at this
    poison silently reproduced the CORRECT formula).
    """
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    L1 = np.zeros(N)
    for nn in range(N):
        s = sum(2.0 * F[k] for k in range(1, nn))     # k = 1 .. nn-1  <-- the bug
        L1[nn] = s / (2.0 * nn + 1.0)
    M = np.empty((N, N))
    for i in range(N):
        for j in range(N):
            val = F[abs(i - j)] * u[i] * u[j] if i != j else 0.0
            val += F[i + j + 1] * u[i] * u[j]
            if i == j:
                val -= L1[i]
            M[i, j] = val
    return n, u, L1, M


w, a2f, cum = read_a2f(A2F)
lam, p, om2 = P_measure(w, a2f)
F = F_T(2 * N, w, p, T)
n, u, L1_ok, A_ok = O_of_f(F, N)
_, _, L1_bad, A_bad = O_poisoned(F, N)

# Sanity check on the CONTROL ITSELF: the poison must actually change something.
dmax = float(np.max(np.abs(L1_bad - L1_ok)))
print(f"poison self-test: max |L1_bad - L1_ok| = {dmax:.6e}  "
      f"({'ok, poison is live' if dmax > 1e-6 else '*** POISON IS INERT ***'})")
assert dmax > 1e-6, "the poison reproduces the correct formula; fix the control"

print(f"v3 interval build: N={N}, T={T} meV, mu*=0, dps={DPS} ... ", end="", flush=True)
A_iv, u_iv, fker_iv, L1_iv = build_M_iv(w, a2f, lam, 0.0, 10.0, T, N, DPS, nc_mode='floor')
print("done")
print()

fails = []
print(f"{'quantity':>12} {'correct baseline':>24} {'encloses?':>10} | "
      f"{'poisoned (old bug)':>24} {'encloses?':>10} {'verdict':>12}")
for idx in (1, 5, 20, N - 1):
    iv = L1_iv[idx]
    lo, hi = float(iv.a), float(iv.b)
    e_ok = lo <= L1_ok[idx] <= hi
    e_bad = lo <= L1_bad[idx] <= hi
    ok = e_ok and not e_bad
    fails.append(ok)
    print(f"{'L1[%d]' % idx:>12} {L1_ok[idx]:>24.16e} {str(e_ok):>10} | "
          f"{L1_bad[idx]:>24.16e} {str(e_bad):>10} {'OK' if ok else 'PROBLEM':>12}")

for (i, j) in ((5, 5), (20, 20), (N - 1, N - 1), (20, 21)):
    iv = A_iv[i][j]
    lo, hi = float(iv.a), float(iv.b)
    e_ok = lo <= A_ok[i, j] <= hi
    e_bad = lo <= A_bad[i, j] <= hi
    # off-diagonal elements are identical under the poison (L1 is diagonal-only):
    identical = abs(A_ok[i, j] - A_bad[i, j]) < 1e-300
    ok = e_ok and (not e_bad or identical)
    fails.append(ok)
    tag = "OK (off-diag unaffected)" if identical and ok else ("OK" if ok else "PROBLEM")
    print(f"{'A[%d,%d]' % (i, j):>12} {A_ok[i, j]:>24.16e} {str(e_ok):>10} | "
          f"{A_bad[i, j]:>24.16e} {str(e_bad):>10} {tag:>12}")

print()
print("relative interval widths (a width of exactly 0 after divisions/sqrt is a RED FLAG,")
print("it usually means the interval collapsed to a point rather than 'dps is ample'):")
for idx in (1, 20, N - 1):
    iv = L1_iv[idx]
    mid = (float(iv.a) + float(iv.b)) / 2
    relw = (float(iv.b) - float(iv.a)) / abs(mid) if mid else float('nan')
    print(f"  L1[{idx:3d}]: rel width = {relw:.3e}")

print()
print("VERDICT:", "negative control has discriminating power"
      if all(fails) else "*** CHECK IS NOT DISCRIMINATING -- fix the check first ***")
