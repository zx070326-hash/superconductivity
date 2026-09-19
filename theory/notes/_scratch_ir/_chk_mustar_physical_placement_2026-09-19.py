"""Physical placement of the mu*>0 extremality thresholds (2026-09-19).

Two questions, both answered from scratch (no number is taken from a note):

(A) Physical filter.  The triple (varpi, N_c, mutilde) is NOT free: on the Tc
    surface  lambda * lam1(M) = 1  and  mutilde = mu*/lambda, so
        lambda = 1/Q,   mu* = mutilde * lambda = mutilde / Q,   Q := lam1(M_E).
    A threshold point is physical only if the mu* it implies is a plausible
    Coulomb pseudopotential.  (Repo's failure map reaches the same conclusion by
    a different route; this is an independent recomputation.)

(B) Anchor placement.  For the repo's (lambda, mu*, r) anchors, solve the Tc
    condition for varpi, get N_c = floor(r*varpi) and mutilde = mu*/lambda, then
    test at that point: (i) is v_P in the cone C (Proposition P3), (ii) is the
    Rayleigh step <v_P,O[d]v_P> >= 0, (iii) does extremality hold -- each over a
    two-atom shape family (legitimate by Theorem TA).

float64, finite N, finite shape grid ==> numerical-evidence.
ONE-SIDED: "violation found" is robust; "no violation on this grid" is not a proof.
"""

import importlib.util
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "cr_op", os.path.join(_HERE, "_chk_mustar_conerepair_op_2026-09-19.py")
)
op = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(op)
_spec2 = importlib.util.spec_from_file_location(
    "cr_scan", os.path.join(_HERE, "_chk_mustar_conerepair_scan_2026-09-19.py")
)
scan = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(scan)


def lam1_E(w, Nc, mut, N=None):
    N = N or int(max(400, 26 * w, 5 * Nc + 40))
    H = op.O_of_kernel(N, op.kernel_einstein(N, w))
    return np.linalg.eigvalsh(op.pencil(H, N, Nc, mut))[-1]


# ------------------------------------------------------------ (A) filter -----
def physical_filter():
    print("\n### (A) physical filter at the measured thresholds")
    print("    lam = 1/Q(E),  mu* = mutilde/Q(E).")
    print("    ONE-SIDED WARNING: our mut_extr comes from a COARSE shape grid, so it is")
    print("    an UPPER bound on the true threshold; hence the implied mu* is also an")
    print("    upper bound.  This table therefore CANNOT establish a lower bound on")
    print("    min mu*; it can only be checked for CONSISTENCY with the repo's denser")
    print("    scan (mustar_extremality_failure_map sec.3 reports min mu* = 1.573 for")
    print("    varpi<=5 and 0.828 over an extended range).")
    print(f"{'wrms':>5} {'Nc':>4} {'r=Nc/w':>7} {'mut_extr':>10} {'Q(E)':>12} "
          f"{'lambda':>10} {'mu* implied':>12} {'>= repo 0.828?':>15}")
    rows = [(0.4, 1, 4.713e-1), (0.4, 3, 9.651e-2), (0.6, 1, 1.041e0),
            (0.6, 3, 1.904e-1), (1.0, 1, 2.886e0), (1.0, 3, 4.208e-1),
            (2.0, 3, 1.463e0), (2.0, 6, 5.911e-1), (2.0, 10, 1.306e0)]
    for (w, Nc, mut) in rows:
        Q = lam1_E(w, Nc, mut)
        lam = 1.0 / Q
        mus = mut / Q
        print(f"{w:5.1f} {Nc:4d} {Nc/w:7.2f} {mut:10.3e} {Q:12.5e} {lam:10.2f} "
              f"{mus:12.4f} {'consistent' if mus >= 0.828 else 'CONFLICT':>15}")


# ------------------------------------------------------------ (B) anchors ----
def solve_varpi(lam, mu, r, lo=0.08, hi=8.0, npts=90, iters=45):
    """Largest varpi on the bracket with lam*lam1(M_E(varpi)) = 1,
    N_c = floor(r*varpi), mutilde = mu/lam.  (Tc := sup, so take the largest root.)"""
    mut = mu / lam
    def g(w):
        return lam * lam1_E(w, int(np.floor(r * w)), mut) - 1.0
    ws = np.geomspace(lo, hi, npts)
    vals = np.array([g(w) for w in ws])
    sign_change = np.flatnonzero((vals[:-1] < 0) & (vals[1:] >= 0))
    if len(sign_change) == 0:
        return None
    i = sign_change[-1]
    a, b = ws[i], ws[i + 1]
    for _ in range(iters):
        m = 0.5 * (a + b)
        if g(m) < 0:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


MU_STAR = 0.13
# mutilde anchor set TRANSCRIBED from
#   theory/notes/coulomb/mustar_extremality_failure_map_2026-09-17.md  sec.6
# lambda := MU_STAR / mutilde is OUR arithmetic, not a literature lambda for any
# named material.  No material names are attached here on purpose: the material
# provenance of these seven numbers has not been re-verified in this session.
ANCHOR_MUTILDE = (0.0276, 0.0375, 0.0485, 0.0578, 0.0870, 0.1057, 0.2097)
ANCHORS = [(f"mut={m:.4f}", MU_STAR / m, MU_STAR) for m in ANCHOR_MUTILDE]


def anchor_table():
    print("\n### (B) anchor placement and criterion status (two-atom family, "
          f"{len(scan.SHAPES)} shapes)")
    print(f"{'anchor':>13} {'lam':>5} {'mu*':>5} {'r':>4} {'varpi':>8} {'Nc':>4} "
          f"{'mutilde':>9} {'cone':>6} {'rayl':>6} {'extr':>6} {'worst rayl':>12}")
    for (name, lam, mu) in ANCHORS:
        for r in (3.0, 10.0):
            w = solve_varpi(lam, mu, r)
            if w is None:
                print(f"{name:>13} {lam:5.2f} {mu:5.2f} {r:4.0f}   (no Tc root on bracket)")
                continue
            Nc = int(np.floor(r * w))
            mut = mu / lam
            cone_ok = rayl_ok = extr_ok = True
            worst = np.inf
            for (weight, spread) in scan.SHAPES:
                ws, ps = scan.two_atom(w, weight, spread)
                N = scan.nsize(w, Nc, ws)
                d = scan.diagnose(w, ws, ps, Nc, mut, N)
                cone_ok &= d["in_cone"]
                rayl_ok &= d["exact_ok"]
                extr_ok &= d["extr_ok"]
                worst = min(worst, d["exact"])
            print(f"{name:>13} {lam:5.2f} {mu:5.2f} {r:4.0f} {w:8.4f} {Nc:4d} "
                  f"{mut:9.4f} {'Y' if cone_ok else 'n':>6} {'Y' if rayl_ok else 'n':>6} "
                  f"{'Y' if extr_ok else 'n':>6} {worst:12.4e}")


if __name__ == "__main__":
    physical_filter()
    anchor_table()
