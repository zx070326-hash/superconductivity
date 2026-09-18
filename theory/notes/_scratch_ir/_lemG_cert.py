"""Rational-exact certificates for Lemma (G) partial results.

CERT-G1:  A := sup_K theta_K/q_K(v),  theta_K = 2K/(2K+1),
          v_n = u_n t_n with t = (9,5,2,1).  Uses q_K = q_inf for K >= 2p-1.
          => |lmin(O[f])| <= (1+2A) lam1(O[f])  for every nonincreasing f >= 0 in l^1.

CERT-G2 (no-go): rational Rayleigh lower bound for lmax(|O[1_{k<=1}]|), to be compared
          with lam1(O[1_{k<=1}]) <= g(2) <= 1.3181615784372800 (repo CERT-B, Lemma M
          with G(k)=1/k^2).  If the Rayleigh value exceeds that, the elementwise-
          absolute-value (Perron) family provably cannot deliver c <= 1.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
from fractions import Fraction as F
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _lemG_exact_gamma import qK_exact, qinf_exact

# ---------------- CERT-G1 ----------------
T = [F(9), F(5), F(2), F(1)]
p = len(T)
print("CERT-G1   t =", [str(x) for x in T], " p =", p, " 2p-1 =", 2 * p - 1)

# v nonneg nonincreasing?  v_n = u_n t_n,  check v_n^2 >= v_{n+1}^2
vsq = [F(T[n] ** 2, 2 * n + 1) for n in range(p)]
print("  v_n^2 =", [str(x) for x in vsq], " nonincreasing:",
      all(vsq[n] >= vsq[n + 1] for n in range(p - 1)), "(v in cone C)")

qs = {}
for K in range(1, 2 * p):
    qs[K] = qK_exact(T, K)
qinf = qinf_exact(T)
print(f"  q_inf (via (S1): 2<u,v>^2/||v||^2 - 1) = {qinf} = {float(qinf):.12f}")
print(f"  q_{2*p-1} (direct)                     = {qs[2*p-1]} = {float(qs[2*p-1]):.12f}")
assert qs[2 * p - 1] == qinf, "saturation identity q_K = q_inf for K >= 2p-1 FAILED"
print("  saturation identity VERIFIED (exact rational equality)")

rows = []
for K in range(1, 2 * p - 1):
    th = F(2 * K, 2 * K + 1)
    rows.append((K, th, qs[K], th / qs[K]))
tail = F(1) / qinf
A = max([r[3] for r in rows] + [tail])
print("   K   theta_K        q_K                       theta_K/q_K")
for K, th, q, r in rows:
    print(f"  {K:2d}  {str(th):8s}  {str(q):26s}  {str(r):22s} = {float(r):.9f}")
print(f"  tail K>={2*p-1}: theta_K/q_K < 1/q_inf = {tail} = {float(tail):.9f}")
print(f"  A = {A} = {float(A):.9f}")
print(f"  ==> c = 1 + 2A = {1 + 2*A} = {float(1 + 2*A):.9f}")
print(f"  rational count: {p} integers t_n + {2*p-2} values q_K + q_inf + A + c "
      f"= {p + 2*p - 2 + 3} rationals; grid = {2*p-2} cells (K=1..{2*p-2}) + 1 tail cell")

# ---------------- CERT-G2 ----------------
print()
print("CERT-G2  (no-go for the elementwise-absolute-value / Perron family)")


def rayleigh_absO1(s):
    """<z,|O[1_{k<=1}]|z>/||z||^2 with z_n = u_n s_n, s rational, finite support.
    |O[1_{<=1}]| : diag (1, 2/3, 2/5, ...), offdiag u_n u_{n+1}."""
    P = len(s)
    us = [F(1, 2 * n + 1) for n in range(P)]
    nrm = sum(us[n] * s[n] ** 2 for n in range(P))
    dg = [F(1)] + [F(2, 2 * n + 1) for n in range(1, P)]
    tot = sum(dg[n] * us[n] * s[n] ** 2 for n in range(P))
    tot += 2 * sum(us[n] * us[n + 1] * s[n] * s[n + 1] for n in range(P - 1))
    return tot / nrm


import numpy as np
best = None
rng = np.random.default_rng(7)
# seed from the float top eigenvector of |O[1_{<=1}]| (float only picks the vector;
# the certified value is the exact rational Rayleigh quotient below)
P = 14
d = np.array([1.0] + [2.0 / (2 * n + 1) for n in range(1, P)])
u = 1.0 / np.sqrt(2 * np.arange(P) + 1.0)
Mx = np.diag(d)
for n in range(P - 1):
    Mx[n, n + 1] = Mx[n + 1, n] = u[n] * u[n + 1]
w, V = np.linalg.eigh(Mx)
vtop = np.abs(V[:, -1])
print(f"  float lmax(|O[1_<=1]|) on {P}x{P} head = {w[-1]:.9f} (navigation only)")
stop = vtop / u                       # s_n = z_n/u_n
stop = stop / stop[0] * 100
for scale in (10, 20, 50, 100, 200):
    s = [F(int(round(x * scale / 100))) for x in stop]
    while s and s[-1] == 0:
        s.pop()
    if not s:
        continue
    val = rayleigh_absO1(s)
    if best is None or val > best[0]:
        best = (val, s)
val, s = best
G2_UB = F(13181615784372800, 10000000000000000)   # repo CERT-B upper bound for g(2)
print(f"  s = {[str(x) for x in s]}")
print(f"  lmax(|O[1_<=1]|) >= {val} = {float(val):.12f}")
print(f"  lam1(O[1_<=1])   <= g(2) <= {float(G2_UB):.12f}  (Lemma M, G(k)=1/k^2; CERT-B)")
print(f"  ratio >= {float(val/G2_UB):.9f}   ->  no-go holds: {val > G2_UB}")
