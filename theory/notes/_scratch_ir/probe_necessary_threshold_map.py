"""NECESSARY-threshold map: mu~_viol(xbar, N_c) := the smallest mu~ at which some
2-atom P with the same <omega^2> has Q(P) > Q(E).

This is the COMPLEMENT of Theorem CE.  CE gives a SUFFICIENT threshold
mu~_*^CE below which extremality is certified; mu~_viol gives a NECESSARY one
above which extremality is FALSE.  The true breakdown threshold lies in
[mu~_*^CE, mu~_viol].  If mu~_viol is far above every physical mu~ = mu*/lambda
at physically realizable (varpi_rms, N_c), extremality holds at the anchors as
a matter of fact (numerical-evidence), whatever the sufficient side certifies.

Physical realizability: varpi_rms >= 0.8709/sqrt(lambda) (analytic floor from
Tc <= 0.1827 sqrt(lambda <w^2>) at mu*=0; mu*>0 lowers Tc, raising varpi, so
the floor holds a fortiori) and N_c = floor(omega_c/2 pi T) >= floor(r varpi_rms)
with r = omega_c/omega_rms >= ~3 (the cutoff sits above the phonon spectrum).

Search over P: w in W, spread in SP -- the extreme shapes (w->1, spread->inf,
i.e. a tiny weight pushed to very high frequency) were seen to be the worst.
mu~_viol is found by bisection on the indicator 'any violation over (w,sp)'.
Because the search over P is finite, mu~_viol here is an UPPER estimate of the
true necessary threshold: a finer P search could only lower it.  Direction:
"mu~ < mu~_viol" therefore proves NOTHING; "mu~ >= mu~_viol" is a genuine
counterexample.  The physical-safety conclusion rests on the anchors' mu~ being
below the SUFFICIENT threshold or on direct ground truth, not on this map alone.

Strictness: numerical-evidence (float64, N=200; convergence in N verified in
probe_lowmu_violations.py to <1e-8 relative for these quantities).
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


W = (0.5, 0.75, 0.9, 0.97, 0.99)
SP = (1.05, 2.0, 10.0, 1e2, 1e3, 1e5)


def any_violation(xbar, Nc, mut, AE, APs):
    uc = np.where(n < Nc, u, 0.0)
    K = 2.0 * np.outer(uc, uc)
    QE = float(eigh(AE - mut * K, eigvals_only=True)[-1])
    best = -np.inf; arg = None
    for key, AP in APs.items():
        QP = float(eigh(AP - mut * K, eigvals_only=True)[-1])
        if QP - QE > best:
            best = QP - QE; arg = key
    return best > 0, best, arg


def mu_viol(xbar, Nc, hi=2.0, tol=2e-3):
    AE = O(F_of([(1.0, xbar)]))
    APs = {(w, sp): O(F_of(pair(xbar, w, sp))) for w in W for sp in SP}
    ok_hi, g_hi, arg_hi = any_violation(xbar, Nc, hi, AE, APs)
    if not ok_hi:
        return None, None
    lo = 0.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        ok, g, arg = any_violation(xbar, Nc, mid, AE, APs)
        if ok:
            hi = mid; arg_hi = arg
        else:
            lo = mid
    return hi, arg_hi


XB = (0.04, 0.09, 0.16, 0.21, 0.25, 0.36, 0.5, 0.75, 1.0, 2.25)
NC = (1, 2, 3, 4, 5, 6, 8)

print("=" * 100)
print("NECESSARY threshold mu~_viol(varpi, N_c): smallest mu~ with a 2-atom violation")
print(f"  P-search: w in {W}, spread in {SP};  bisection tol 2e-3;  N={N}")
print("  '--' = no violation found even at mu~ = 2.0")
print("=" * 100)
print(f"{'varpi':>7} {'xbar':>6} {'lam_min':>8} |" +
      "".join(f"{'Nc='+str(Nc):>10}" for Nc in NC))
print(f"{'':>7} {'':>6} {'(floor)':>8} |")
table = {}
for xbar in XB:
    vp = math.sqrt(xbar)
    lam_floor = (0.8709 / vp) ** 2
    row = []
    for Nc in NC:
        mv, arg = mu_viol(xbar, Nc)
        table[(xbar, Nc)] = (mv, arg)
        row.append(f"{mv:>10.3f}" if mv is not None else f"{'--':>10}")
    print(f"{vp:>7.3f} {xbar:>6.3g} {lam_floor:>8.1f} |" + "".join(row))

print()
print("  lam_min(floor) = smallest lambda for which varpi_rms can be this small")
print("  (varpi_rms >= 0.8709/sqrt(lambda)).  Anchors: lambda <= 3.62 => varpi >= 0.458")
print()
print("  worst P-shape (w, spread) at the threshold, per cell:")
for xbar in XB:
    parts = []
    for Nc in NC:
        mv, arg = table[(xbar, Nc)]
        parts.append(f"({arg[0]:.2f},{arg[1]:.0e})" if arg else "   --   ")
    print(f"  varpi={math.sqrt(xbar):.3f}: " + "  ".join(parts))

# --- physical anchors against the map: interpolate nothing; use the nearest
#     harder cell (smaller varpi, N_c = floor(r*varpi) with r=3 -> smallest N_c)
print()
print("=" * 100)
print("ANCHORS vs the necessary threshold (nearest HARDER grid cell: varpi rounded DOWN)")
print("=" * 100)
ANCH = [("LaH10 129GPa", 3.62, 0.10, None), ("LaH10 163GPa", 2.67, 0.10, 0.783),
        ("LaH10 214GPa", 2.06, 0.10, None), ("LaH10 264GPa", 1.73, 0.10, None),
        ("H3S 200GPa", 1.84, 0.16, None), ("Nb", 1.23, 0.13, None),
        ("MgB2", 0.62, 0.13, None)]
print(f"{'anchor':>13} {'mu~':>7} {'varpi_rms':>10} {'cell varpi':>11} {'Nc(r=3)':>8} "
      f"{'mu~_viol':>9} {'ratio':>7}")
for label, lam, mus, vr in ANCH:
    vr = vr if vr else 0.8709 / math.sqrt(lam)
    mut = mus / lam
    cell = max(x for x in XB if math.sqrt(x) <= vr)
    Nc = max(1, int(3 * vr))
    Nc_cell = max(c for c in NC if c <= Nc)
    mv, arg = table[(cell, Nc_cell)]
    print(f"{label:>13} {mut:>7.4f} {vr:>10.3f} {math.sqrt(cell):>11.3f} {Nc_cell:>8} "
          f"{(mv if mv else float('nan')):>9.3f} "
          f"{(mv/mut if mv else float('inf')):>7.1f}")
print()
print("  ratio = mu~_viol / mu~_anchor : how far the anchor sits below the first")
print("  known counterexample.  (>1 is consistent with the 42/42 ground-truth PASS;")
print("  it is NOT a proof -- see direction note in the docstring.)")
