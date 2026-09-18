"""PHYSICAL FILTER on the counterexamples: which violations Q(P)>Q(E) actually
break the Tc ordering  Tc(P) <= Tc(E)  at a physical (lambda, mu*)?

The eigenvalue comparison is made at a fixed temperature T (fixed varpi, N_c,
mu~ = mu*/lambda).  Tc is where Q(.;T) = 1/lambda.  A violation point breaks
the Tc ordering for exactly those mu* with
      mu~/Q(P)  <=  mu*  <  mu~/Q(E)                        (window W)
[left: T is at or below Tc(P), i.e. Q(P) >= 1/lambda = mu~/mu*;
 right: T is above Tc(E), i.e. Q(E) < 1/lambda].  The window is non-empty iff
Q(P) > Q(E).  So the physically decisive number for each violation is
      mu*_min := mu~ / Q(P),        with  lambda = 1/Q(P)  there.
If mu*_min >> 0.2 for EVERY violation, then no physical (lambda, mu*) has its
Tc ordering broken by any of them: the violations live at temperatures where
Coulomb has already killed both eigenvalues (Q ~ 1e-4), i.e. far above Tc of
both P and E for any physical coupling.

Direction: this script only CLASSIFIES existing violations; a violation with
small mu*_min would be a genuine physical counterexample (no bounds anywhere).
A large min over a FINITE grid proves nothing; it is numerical-evidence.
Strictness: numerical-evidence (float64, N=200).
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

N = 200
n = np.arange(N); u = 1.0 / np.sqrt(2.0 * n + 1.0)
ks = np.arange(2 * N + 1, dtype=float); s = ks ** 2
absd = np.abs(n[:, None] - n[None, :]); sumd = n[:, None] + n[None, :] + 1
uu = u[:, None] * u[None, :]


def O(f):
    A = (f[absd] + f[sumd]) * uu
    L1 = np.zeros(N); L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * n + 1] * u ** 2)
    return A


def F_of(atoms):
    F = np.zeros_like(s)
    for w, x in atoms:
        F += w * x / (s + x)
    F[0] = sum(w for w, _ in atoms)
    return F


def pair(xbar, w, sp):
    x1 = xbar / (w + (1 - w) * sp)
    return [(w, x1), (1 - w, sp * x1)]


XB = (0.04, 0.09, 0.25, 0.5, 1.0, 2.25, 4.0, 9.0, 25.0)
NC = (1, 2, 3, 4, 6, 10, 15, 25, 40)
MU = (0.03, 0.05, 0.08, 0.13, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0, 3.0, 5.0)
W = (0.5, 0.75, 0.9, 0.97, 0.99)
SP = (1.05, 2.0, 10.0, 1e2, 1e3, 1e5)

hits = []
tot = 0
for xbar in XB:
    AE = O(F_of([(1.0, xbar)]))
    APs = {(w, sp): O(F_of(pair(xbar, w, sp))) for w in W for sp in SP}
    for Nc in NC:
        uc = np.where(n < Nc, u, 0.0); K = 2.0 * np.outer(uc, uc)
        for mut in MU:
            QE = float(eigh(AE - mut * K, eigvals_only=True)[-1])
            for key, AP in APs.items():
                tot += 1
                QP = float(eigh(AP - mut * K, eigvals_only=True)[-1])
                if QP > QE:
                    hits.append((mut / QP, mut, xbar, Nc, key[0], key[1], QP, QE))

print("=" * 100)
print(f"scanned {tot} points; violations Q(P)>Q(E): {len(hits)}")
print("=" * 100)
hits.sort()
print("  violations ranked by mu*_min = mu~/Q(P)  (smallest = physically most dangerous)")
print(f"  {'mu*_min':>9} {'lambda':>8} {'mu~':>6} {'varpi':>7} {'Nc':>4} {'r=Nc/varpi':>10} "
      f"{'w':>5} {'spread':>8} {'Q(P)':>10} {'Q(E)':>10} {'Q(P)/Q(E)':>10}")
for h in hits[:25]:
    ms, mut, xbar, Nc, w, sp, QP, QE = h
    vp = math.sqrt(xbar)
    print(f"  {ms:>9.3f} {1/QP:>8.3f} {mut:>6.2f} {vp:>7.3f} {Nc:>4} {Nc/vp:>10.2f} "
          f"{w:>5.2f} {sp:>8.0e} {QP:>10.3e} {QE:>10.3e} {QP/QE:>10.3f}")
print()
print("  distribution of mu*_min over all violations:")
ms = np.array([h[0] for h in hits])
for q in (0, 1, 5, 25, 50):
    print(f"    percentile {q:>2}: mu*_min = {np.percentile(ms, q):.3f}")
print(f"    count with mu*_min < 0.25 : {int(np.sum(ms < 0.25))}")
print(f"    count with mu*_min < 1.0  : {int(np.sum(ms < 1.0))}")
print(f"    count with lambda=1/Q(P) < 5 : {int(np.sum(np.array([1/h[6] for h in hits]) < 5))}")
