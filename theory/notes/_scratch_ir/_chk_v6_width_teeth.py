"""Teeth test for v6's NEW width assertion (rw > 1e-40 -> violation).

A check that cannot fail is worthless (the F5 lesson).  v6 added
    if abs(baseline_val) > 1e-100 and rw > WIDTH_THRESHOLD: width_violations.append(...)
and reports "All intervals < 1e-40 (PASS)".  Question: can it EVER say no?

Method: rebuild the same quantities at a range of dps and ask at which dps the
threshold fires.  If it fires when precision is degraded (the F2 regression this
assertion exists to catch) the check has teeth; if it never fires, it is decoration.

Also records the blind spot: the guard keys on baseline_val, not on rw, so an
interval that COLLAPSED to zero width still passes.
"""
import glob
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from h3s_deficit_cert_v6 import ROOT, moments, build_M_iv, rel_width

WIDTH_THRESHOLD = 1e-40
N = 40
T_test = 1.0
r = 10.0

g = glob.glob(f"{ROOT}/H3S-200GPa/EPW_anharmonic/*.a2f")
g = [x for x in g if not x.endswith(".a2fv")]
lam, om2, olog, w, a2f = moments(g[0])

# the same twelve quantities the alignment check looks at
F_IDX = [0, 1, 20, 80]
L1_IDX = [0, 1, 20, 39]
A_IDX = [(0, 0), (0, 1), (20, 20), (20, 21)]

print("=" * 78)
print("v6 width assertion: does it ever fire?   threshold rw > %.0e" % WIDTH_THRESHOLD)
print("=" * 78)
print()
print("%6s | %11s | %11s | %11s | %s" % ("dps", "max rw", "min rw", "A[20,20] rw", "verdict"))
print("-" * 78)

for dps in (12, 15, 20, 30, 50, 80):
    M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(
        w, a2f, lam, 0.0, r, T_test, N, dps, Nc_fixed=None, poison=False
    )
    rws = []
    a2020 = None
    for i in F_IDX:
        rws.append(("F[%d]" % i, rel_width(fker_iv[i], dps)))
    for i in L1_IDX:
        rws.append(("L1[%d]" % i, rel_width(L1_iv[i], dps)))
    for (i, j) in A_IDX:
        v = rel_width(M_iv[i][j], dps)
        rws.append(("A[%d,%d]" % (i, j), v))
        if (i, j) == (20, 20):
            a2020 = v

    # exactly v6's predicate, minus the baseline-magnitude guard (F[0], L1[0] are
    # exact by definition so they read 0 and are excluded there too)
    nz = [v for (nm, v) in rws if v > 0.0]
    fired = [nm for (nm, v) in rws if v > WIDTH_THRESHOLD]
    print("%6d | %11.3e | %11.3e | %11.3e | %s" % (
        dps, max(v for _, v in rws), min(nz) if nz else 0.0, a2020,
        ("FIRES on %d/%d: %s" % (len(fired), len(rws), ", ".join(fired[:4]) + ("..." if len(fired) > 4 else "")))
        if fired else "silent (PASS)"))

print()
print("=" * 78)
print("BLIND SPOT CHECK: the guard is `abs(baseline_val) > 1e-100`, keyed on the")
print("BASELINE, not on rw.  So a gauge that collapses to rw == 0 -- the exact F5")
print("failure this assertion was added to backstop -- still reports PASS.")
print("=" * 78)
print("  a one-sided 'too wide' assertion cannot detect 'too narrow / reads zero'.")
print("  v6 prints the per-element rw, so a human can see it; nothing automatic fires.")
