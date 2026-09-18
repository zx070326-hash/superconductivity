"""Zoom on the LOW-mu~ violations: are the LaH10-regime points (mu~ ~ 0.03-0.06)
actually in danger, or are the violations confined to measure shapes that a real
alpha^2 F cannot have?

The broad scan (probe_mustar_counterexample.py) reported violations down to
mu~ = 0.05.  mu~ = mu*/lambda for LaH10 is 0.028-0.058, so this must be resolved
before any statement about physical anchors is made.  Three questions:

  Q1. At mu~ <= 0.13, WHICH (xbar, N_c, w, spread) violate, and by how much?
  Q2. Do those violations survive N-refinement (real) or vanish (truncation)?
  Q3. What distinguishes a violating P from a physical alpha^2 F?  In particular
      the violations seen so far have extreme weight asymmetry (w=0.9) and
      spread x2/x1 = 200, i.e. a tiny high-frequency atom carrying almost no
      weight -- check whether the violation needs spread -> large.

Strictness: numerical-evidence (float64, finite N).  Everything is a direct
top-eigenvalue computation of an explicitly assembled matrix: no bounds, no
one-sided substitutions, so a violation is a real finding at this N.
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

_c = {}


def setup(N):
    if N in _c:
        return _c[N]
    n = np.arange(N); u = 1.0 / np.sqrt(2.0 * n + 1.0)
    ks = np.arange(2 * N + 1, dtype=float)
    _c[N] = dict(N=N, n=n, u=u, s=ks ** 2,
                 absd=np.abs(n[:, None] - n[None, :]),
                 sumd=n[:, None] + n[None, :] + 1,
                 uu=u[:, None] * u[None, :])
    return _c[N]


def O(f, D):
    N, n, u = D["N"], D["n"], D["u"]
    A = (f[D["absd"]] + f[D["sumd"]]) * D["uu"]
    L1 = np.zeros(N); L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * n + 1] * u ** 2)
    return A


def Q(atoms, mut, Nc, D):
    s = D["s"]; F = np.zeros_like(s)
    for w, x in atoms:
        F += w * x / (s + x)
    F[0] = sum(w for w, _ in atoms)
    M = O(F, D)
    if mut:
        uc = np.where(D["n"] < Nc, D["u"], 0.0)
        M = M - 2.0 * mut * np.outer(uc, uc)
    return float(eigh(M, eigvals_only=True)[-1])


def pair(xbar, w, sp):
    x1 = xbar / (w + (1 - w) * sp)
    return [(w, x1), (1 - w, sp * x1)]


# -------------------------------------------------- Q1
print("=" * 104)
print("Q1: violations at LOW mu~ (the LaH10 regime is mu~ = 0.028-0.058)")
print("=" * 104)
D = setup(300)
low = []
for xbar in (0.04, 0.09, 0.25, 0.5, 1.0, 4.0, 9.0):
    for Nc in (1, 2, 3, 4, 6, 10, 20):
        for mut in (0.03, 0.05, 0.08, 0.13):
            QE = Q([(1.0, xbar)], mut, Nc, D)
            for w in (0.1, 0.25, 0.5, 0.75, 0.9, 0.97):
                for sp in (1.05, 1.3, 2.0, 5.0, 20.0, 200.0, 2000.0):
                    QP = Q(pair(xbar, w, sp), mut, Nc, D)
                    if QP > QE:
                        low.append((QP - QE, (QP - QE) / abs(QE), xbar, Nc,
                                    mut, w, sp, QP, QE))
print(f"scanned {7*7*4*6*7} points at mu~ in (0.03,0.05,0.08,0.13)"
      f"  ->  {len(low)} violations")
if low:
    low.sort(reverse=True)
    print()
    print(f"  {'Q(P)-Q(E)':>12} {'rel':>10} {'xbar':>7} {'varpi':>7} {'Nc':>4} "
          f"{'mu~':>6} {'w':>6} {'spread':>9} {'x1':>10} {'x2':>11}")
    for h in low[:20]:
        g, rel, xbar, Nc, mut, w, sp, QP, QE = h
        x1 = xbar / (w + (1 - w) * sp)
        print(f"  {g:>12.4e} {rel:>10.3e} {xbar:>7.3g} {math.sqrt(xbar):>7.3f} "
              f"{Nc:>4} {mut:>6.2f} {w:>6.2f} {sp:>9.4g} {x1:>10.4g} {sp*x1:>11.4g}")
    print()
    print("  spread values that appear among violations:",
          sorted({h[6] for h in low}))
    print("  weights  values that appear among violations:",
          sorted({h[5] for h in low}))
    print("  N_c      values that appear among violations:",
          sorted({h[3] for h in low}))
    print("  min mu~ with a violation:", min(h[4] for h in low))

    # ---------------------------------------------- Q2
    print()
    print("=" * 104)
    print("Q2: N-CONVERGENCE of the 4 strongest low-mu~ violations")
    print("=" * 104)
    for h in low[:4]:
        g, rel, xbar, Nc, mut, w, sp, QP, QE = h
        print(f"  xbar={xbar:g} N_c={Nc} mu~={mut:g} w={w:g} spread={sp:g}")
        print(f"    {'N':>6} {'Q(P)':>15} {'Q(E)':>15} {'Q(P)-Q(E)':>13} {'viol?':>7}")
        for N in (200, 300, 600, 1200, 2400):
            DD = setup(N)
            qp = Q(pair(xbar, w, sp), mut, Nc, DD)
            qe = Q([(1.0, xbar)], mut, Nc, DD)
            print(f"    {N:>6} {qp:>15.8e} {qe:>15.8e} {qp-qe:>13.4e} "
                  f"{'YES' if qp > qe else 'no':>7}")
        print()

    # ---------------------------------------------- Q3
    print("=" * 104)
    print("Q3: is a LARGE spread necessary?  sweep spread at a violating point")
    print("=" * 104)
    g, rel, xbar, Nc, mut, w, sp, QP, QE = low[0]
    D = setup(600)
    QE6 = Q([(1.0, xbar)], mut, Nc, D)
    print(f"  fixed xbar={xbar:g}, N_c={Nc}, mu~={mut:g}, w={w:g}, N=600")
    print(f"  {'spread':>10} {'x1':>11} {'x2':>12} {'Q(P)-Q(E)':>13} {'viol?':>7}")
    for s2 in (1.05, 1.5, 2, 5, 10, 50, 200, 1e3, 1e4, 1e5, 1e6):
        x1 = xbar / (w + (1 - w) * s2)
        qp = Q(pair(xbar, w, s2), mut, Nc, D)
        print(f"  {s2:>10.4g} {x1:>11.5g} {s2*x1:>12.5g} {qp-QE6:>13.4e} "
              f"{'YES' if qp > QE6 else 'no':>7}")
    print()
    print("  NOTE x2 = spread*x1 is in units of (2 pi T)^2, so x2 = 1e4 means")
    print("  omega = 100*2*pi*T ~ 100x the temperature scale: check whether the")
    print("  violating atom sits at a physically attainable phonon frequency.")
else:
    print("  NO violations at mu~ <= 0.13 on this grid.")
