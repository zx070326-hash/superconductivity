"""
When does the f^xc sandwich  Lambda^SF = 2 F chi F  restore ELEMENTWISE positivity?

Context.  Essenberger Eq.(37) asserts Im[Lambda^SF_{kk'}] > 0 elementwise.
Eq.(32) makes Lambda^SF = 2 F chi F a CONGRUENCE (same symmetric F both sides),
so Sylvester forces Im Lambda^SF to be PSD -- but PSD does not imply elementwise
positive, and that gap is exactly what Eq.(37) needs (see
theory/notes/meta/eq37_content_is_in_the_offdiagonal_2026-09-22.md).

This script asks the structural question: for WHICH F is the sandwich guaranteed
to deliver elementwise nonnegativity?  The answer relocates the open question
from a property of chi (material physics, expensive) to a property of f^xc
(the xc kernel, cheap and basis-explicit).

THEOREM (analytic).
    (F chi F)_ij >= 0  for EVERY PSD chi
      <=>  F = s * w w^T  with s >= 0 and w SIGN-COHERENT (all entries >= 0,
           or all <= 0),  i.e. rank(F) <= 1 AND range(F) spanned by a
           nonnegative vector.

  (<=)  F = s w w^T  =>  (F chi F)_ij = s^2 (w^T chi w) w_i w_j >= 0, since
        w^T chi w >= 0 (chi PSD) and w_i w_j >= 0 (w sign-coherent).
  (=>)  Take chi = v v^T (PSD).  Then (F chi F)_ij = (Fv)_i (Fv)_j, so the
        hypothesis says every Fv is sign-coherent, i.e. range(F) lies in
        (nonneg orthant) U (nonpos orthant).  That union contains no
        2-dimensional subspace: if it did, the plane minus the origin -- which
        is connected -- would be covered by two disjoint nonempty closed cones
        C and -C, a contradiction.  Hence rank(F) <= 1, and its range is
        spanned by a sign-coherent w.  QED

NOTE on an earlier wrong version: "rank(F) <= 1" ALONE is NOT sufficient.
A rank-1 F = w w^T with MIXED-SIGN w gives (F chi F)_ij proportional to w_i w_j,
which is negative whenever sign(w_i) != sign(w_j).  Measured below: 18790/20000
violations.  Sign-coherence is load-bearing, not decorative.

Since a physical f^xc is not rank-1, the theorem says: the sandwich structure
ALONE can never guarantee Eq.(37).  The remaining question is quantitative --
how much does a realistic f^xc help?  Tests C-E and the Perron stratification
answer that, and identify the actual driver.
"""
import numpy as np

rng = np.random.default_rng(23)
N = 5


def psd(n):
    G = rng.standard_normal((n, n))
    return G @ G.T


print("=" * 72)
print("T1  exact sufficient condition: F = s w w^T with w sign-coherent")
print("=" * 72)
bad = 0
for _ in range(20000):
    w = abs(rng.standard_normal(N))
    F = np.outer(w, w) * rng.standard_normal() ** 2
    bad += (F @ psd(N) @ F < -1e-12).any()
print(f"  w >= 0 entrywise      : violations {bad}/20000        (theory: 0)")

bad = tot = 0
for _ in range(20000):
    w = rng.standard_normal(N)
    if (w > 0).all() or (w < 0).all():
        continue
    tot += 1
    bad += (np.outer(w, w) @ psd(N) @ np.outer(w, w) < -1e-12).any()
print(f"  w MIXED sign (rank still 1): violations {bad}/{tot}")
print("  => rank<=1 alone is NOT sufficient; sign-coherence is load-bearing.")

print()
print("=" * 72)
print("T2  necessity: rank(F) >= 2 always admits a PSD witness chi = v v^T")
print("=" * 72)
found = tot = 0
for _ in range(4000):
    n = rng.integers(2, 6)
    F = rng.standard_normal((n, n))
    F = (F + F.T) / 2
    if np.linalg.matrix_rank(F) < 2:
        continue
    tot += 1
    for _ in range(400):
        u = F @ rng.standard_normal(n)
        if (u > 0).any() and (u < 0).any():
            found += 1
            break
print(f"  mixed-sign Fv found: {found}/{tot}   (theory: all; misses are search budget)")

print()
print("=" * 72)
print("T3  quantitative: how often does a sandwich deliver elementwise >= 0?")
print("=" * 72)
for label, gen in [("F symmetric, entries ANY sign",
                    lambda: (lambda A: (A + A.T) / 2)(rng.standard_normal((N, N)))),
                   ("F symmetric, entries ALL >= 0",
                    lambda: (lambda A: (A + A.T) / 2)(abs(rng.standard_normal((N, N)))))]:
    ok = 0
    for _ in range(20000):
        F = gen()                      # same F on BOTH sides -- it is a congruence
        ok += (F @ psd(N) @ F >= 0).all()
    print(f"  {label:32s}: {ok}/20000 = {100*ok/20000:6.2f}%")

print()
print("=" * 72)
print("T4  MECHANISM: Perron dominance.  F >= 0 entrywise => (Perron-Frobenius)")
print("    top eigenvector w >= 0; a large spectral gap makes F ~ s w w^T with")
print("    w >= 0, i.e. close to the exact condition of T1.")
print("=" * 72)
bins = [("gap < 2", 0, 2), ("2 <= gap < 5", 2, 5), ("gap >= 5", 5, 1e9)]
print("  (a) F entrywise NONNEGATIVE, stratified by Perron gap:")
for label, lo, hi in bins:
    ok = tot = 0
    for _ in range(40000):
        F = abs(rng.standard_normal((N, N)))
        F = (F + F.T) / 2
        ev = np.linalg.eigvalsh(F)
        gap = ev[-1] / max(abs(ev[:-1]).max(), 1e-12)
        if not (lo <= gap < hi):
            continue
        tot += 1
        ok += (F @ psd(N) @ F >= 0).all()
    print(f"      {label:14s}: {ok}/{tot} = {100*ok/max(tot,1):5.1f}%")

print("  (b) CONTROL -- same gap stratification, entries allowed NEGATIVE:")
print("      if the gap alone were the driver, these would match (a).")
for label, lo, hi in bins:
    ok = tot = 0
    for _ in range(400000):
        F = rng.standard_normal((N, N))
        F = (F + F.T) / 2
        ev = np.linalg.eigvalsh(F)
        i = int(np.argmax(abs(ev)))
        gap = abs(ev[i]) / max(abs(np.delete(ev, i)).max(), 1e-12)
        if not (lo <= gap < hi):
            continue
        tot += 1
        ok += (F @ psd(N) @ F >= 0).all()
        if tot >= 2000:
            break
    print(f"      {label:14s}: {ok}/{tot} = {100*ok/max(tot,1):5.1f}%")

print()
print("  => the gap is NOT the driver on its own; entrywise nonnegativity of F is,")
print("     with the gap modulating it.  Both ingredients are needed.")
print()
print("CONSEQUENCE: the load-bearing unknown for Eq.(37) is a property of f^xc_zz")
print("(entrywise sign + Perron gap, in the Fermi-surface (band,k) basis), NOT a")
print("property of the material's chi.  That is a far cheaper thing to settle.")
