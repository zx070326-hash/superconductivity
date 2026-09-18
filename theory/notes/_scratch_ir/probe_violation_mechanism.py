"""Two follow-ups on the necessary-threshold map.

(A) LOCAL breakdown: mu~_loc(varpi, N_c) := smallest mu~ at which the NEAR-
    DEGENERATE pair (spread 1.05, w=0.5) already violates.  A violation by an
    infinitesimal perturbation of the Einstein point means the Einstein point
    is not even a LOCAL maximiser, i.e. the second-order coefficient Psi_mu
    (the mu*>0 analogue of Psi = <v,O[psi_0]v>) has turned NEGATIVE.  This
    decides whether 'Psi>0 at mu*>0' can be a theorem: if mu~_loc is finite it
    can only be a theorem BELOW a threshold.

(B) MECHANISM at a violating point: compare the top vectors of M_P and M_E.
    Does v^ leave the cone C = {v>=0, non-increasing}?  Where is the sign
    change relative to N_c?  How degenerate is the top of spec(M_P)?  This is
    the retardation (Anderson-Morel) picture: a Coulomb window of N_c Matsubara
    frequencies lets the gap function go negative just outside the window,
    which lowers the Coulomb cost; P with a soft dominant atom can exploit
    this better than the Einstein comparator.

Strictness: numerical-evidence (float64, N=200/400).
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


def K_of(Nc):
    uc = np.where(n < Nc, u, 0.0)
    return 2.0 * np.outer(uc, uc)


# ---------------------------------------------------------------- (A)
print("=" * 96)
print("(A) LOCAL breakdown mu~_loc(varpi,N_c): near-degenerate pair (w=0.5, spread=1.05)")
print("    finite entry => Einstein is NOT a local max above it (Psi_mu < 0)")
print("=" * 96)
XB = (0.04, 0.09, 0.16, 0.21, 0.25, 0.36, 0.5, 0.75, 1.0, 2.25)
NC = (1, 2, 3, 4, 5, 6, 8)
print(f"{'varpi':>7} |" + "".join(f"{'Nc='+str(Nc):>10}" for Nc in NC))
for xbar in XB:
    AE = O(F_of([(1.0, xbar)])); AP = O(F_of(pair(xbar, 0.5, 1.05)))
    row = []
    for Nc in NC:
        K = K_of(Nc)
        def viol(m):
            return (eigh(AP - m * K, eigvals_only=True)[-1]
                    > eigh(AE - m * K, eigvals_only=True)[-1])
        if not viol(2.0):
            row.append(f"{'--':>10}"); continue
        lo, hi = 0.0, 2.0
        while hi - lo > 2e-3:
            mid = 0.5 * (lo + hi)
            if viol(mid): hi = mid
            else: lo = mid
        row.append(f"{hi:>10.3f}")
    print(f"{math.sqrt(xbar):>7.3f} |" + "".join(row))
print()
print("  compare with the GLOBAL map (worst shape w=0.99, spread=1e5) in")
print("  probe_necessary_threshold_map.py: local thresholds are higher, so")
print("  the first counterexamples are FAR from Einstein, not local.")

# ---------------------------------------------------------------- (B)
print()
print("=" * 96)
print("(B) MECHANISM at violating points: top vectors, cone membership, gap")
print("=" * 96)


def cone(v):
    v = v * np.sign(v[np.argmax(np.abs(v))])
    return bool(np.all(v >= -1e-10) and np.all(np.diff(v) <= 1e-10)), v


def show(xbar, Nc, mut, w, sp, tag):
    AE = O(F_of([(1.0, xbar)])); AP = O(F_of(pair(xbar, w, sp))); K = K_of(Nc)
    lE, WE = eigh(AE - mut * K); lP, WP = eigh(AP - mut * K)
    l0E = eigh(AE, eigvals_only=True); l0P = eigh(AP, eigvals_only=True)
    inE, vE = cone(WE[:, -1]); inP, vP = cone(WP[:, -1])
    x1 = xbar / (w + (1 - w) * sp)
    print(f"  [{tag}] varpi={math.sqrt(xbar):.3f} N_c={Nc} mu~={mut:g} "
          f"P=({w:g} at varpi1={math.sqrt(x1):.4f}, {1-w:g} at varpi2={math.sqrt(sp*x1):.3f})")
    print(f"      Q(P)={lP[-1]:.5e}  Q(E)={lE[-1]:.5e}  "
          f"{'VIOLATION' if lP[-1] > lE[-1] else 'ok'}")
    print(f"      mu~=0 gaps: gamma_P={l0P[-1]-l0P[-2]:.4e}  gamma_E={l0E[-1]-l0E[-2]:.4e}"
          f"   (Q0(P)={l0P[-1]:.4e} <= Q0(E)={l0E[-1]:.4e}: {l0P[-1] <= l0E[-1]})")
    print(f"      perturbed gaps: P {lP[-1]-lP[-2]:.4e}   E {lE[-1]-lE[-2]:.4e}")
    for name, v, inc in (("v^_P", vP, inP), ("v^_E", vE, inE)):
        neg = np.nonzero(v < -1e-10)[0]
        first_neg = int(neg[0]) if neg.size else None
        print(f"      {name}: in cone={inc}  first negative index={first_neg}  "
              f"v[0..7]={np.array2string(v[:8], precision=3, suppress_small=True)}")
    # Coulomb overlap c1 of the top vectors (mu~=0 and perturbed)
    uc = np.where(n < Nc, u, 0.0)
    v0P = eigh(AP)[1][:, -1]; v0E = eigh(AE)[1][:, -1]
    print(f"      Coulomb overlap c1^2: mu~=0  P {float(uc@v0P)**2:.4f}  E {float(uc@v0E)**2:.4f}"
          f"  |  perturbed  P {float(uc@vP)**2:.4f}  E {float(uc@vE)**2:.4f}")
    print(f"      Coulomb energy cost 2mu~c1^2 (perturbed): P {2*mut*float(uc@vP)**2:.4e} "
          f"E {2*mut*float(uc@vE)**2:.4e}")
    print()


show(0.09, 2, 0.13, 0.97, 2000.0, "strong, low mu~")
show(0.04, 1, 0.13, 0.5, 1.05, "near-degenerate, local")
show(0.25, 3, 0.20, 0.99, 1e5, "at anchor varpi, above threshold")
show(0.25, 3, 0.05, 0.99, 1e5, "at anchor varpi, LaH10-like mu~ (should be ok)")
show(0.25, 10, 0.50, 0.99, 1e5, "at anchor varpi, N_c=10, large mu~ (should be ok)")
