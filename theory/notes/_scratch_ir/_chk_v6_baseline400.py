"""Verify v6's hardcoded N=400 BASELINES entries against an independent dps=120 truth.

The N=400 production run's only self-check is the alignment block, which compares
against BASELINES[400].  If those hardcoded numbers are wrong, the check is worse
than useless at production size -- it would either fire spuriously or bless a bug.

F[1], F[20], L1[1], L1[20] and all four A entries are N-INDEPENDENT (already
verified at N=40).  The genuinely new N=400 entries are F[800] and L1[399].
"""
import glob
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_fill_factor import moments, ROOT

DPS_HI, T = 120, 1.0

g = [x for x in glob.glob(f"{ROOT}/H3S-200GPa/EPW_anharmonic/*.a2f") if not x.endswith(".a2fv")]
lam_f, om2_f, olog_f, w, a2f = moments(g[0])

mp.mp.dps = DPS_HI
W = [mp.mpf(float(x)) for x in w]          # EXACT binary float64, not repr()
Y = [mp.mpf(float(x)) for x in a2f]
n_pts = len(W)
tw = [None] * n_pts
tw[0] = (W[1] - W[0]) / 2
tw[-1] = (W[-1] - W[-2]) / 2
for i in range(1, n_pts - 1):
    tw[i] = (W[i + 1] - W[i - 1]) / 2
raw = [(Y[i] / W[i]) * tw[i] for i in range(n_pts)]
tot = mp.fsum(raw)
p = [raw[i] / tot for i in range(n_pts)]
xs = [(W[i] / (2 * mp.pi * mp.mpf(float(T)))) ** 2 for i in range(n_pts)]


def F_hi(k):
    if k == 0:
        return mp.mpf(1)
    k2 = mp.mpf(k) ** 2
    return mp.fsum([p[i] * xs[i] / (k2 + xs[i]) for i in range(n_pts)])


N = 400
Fh = [F_hi(k) for k in range(2 * N + 2)]
L1_hi = [mp.mpf(0)] * N
acc = mp.mpf(0)
for nn in range(1, N):
    acc += 2 * Fh[nn]
    L1_hi[nn] = acc / mp.mpf(2 * nn + 1)          # upper limit INCLUDES k=n
u_hi = [mp.mpf(1) / mp.sqrt(mp.mpf(2 * n + 1)) for n in range(N)]


def A_hi(i, j):
    if i == j:
        return -L1_hi[i] + Fh[2 * i + 1] * u_hi[i] ** 2
    return (Fh[abs(i - j)] + Fh[i + j + 1]) * u_hi[i] * u_hi[j]


# v6's hardcoded BASELINES[400]
CLAIMED = [
    ("F[0]",     1.0000000000000000e+00, Fh[0]),
    ("F[1]",     9.8926738719898688e-01, Fh[1]),
    ("F[20]",    4.7083503166965779e-01, Fh[20]),
    ("F[800]",   8.8050658178847501e-04, Fh[800]),   # N=400-specific
    ("L1[0]",    0.0000000000000000e+00, L1_hi[0]),
    ("L1[1]",    6.5951159146599125e-01, L1_hi[1]),
    ("L1[20]",   6.9152796470690037e-01, L1_hi[20]),
    ("L1[399]",  7.8271859984520367e-02, L1_hi[399]),  # N=400-specific
    ("A[0,0]",   9.8926738719898688e-01, A_hi(0, 0)),
    ("A[0,1]",   1.1290937977616782e+00, A_hi(0, 1)),
    ("A[20,20]", -6.8615029568296637e-01, A_hi(20, 20)),
    ("A[20,21]", 2.8643164353339162e-02, A_hi(20, 21)),
]

print("=" * 86)
print("v6 BASELINES[400] vs independent dps=120 truth    (ALIGN_TOL = 1e-14)")
print("=" * 86)
print("%-10s %24s %24s %12s  %s" % ("quantity", "v6 hardcoded", "dps=120 truth", "rel diff", "verdict"))
print("-" * 86)
worst = 0.0
nbad = 0
for nm, claimed, truth in CLAIMED:
    c = mp.mpf(float(claimed))
    if truth == 0:
        rel = abs(c - truth)
    else:
        rel = abs(c - truth) / abs(truth)
    relf = float(rel)
    worst = max(worst, relf)
    # must sit well inside the 1e-14 alignment tolerance to be usable as a baseline
    ok = relf < 1e-14
    nbad += (not ok)
    print("%-10s %24.16e %24.16e %12.2e  %s" % (
        nm, claimed, float(truth), relf, "ok" if ok else "*** OFF ***"))
print("-" * 86)
print("worst relative deviation: %.2e   (float64 datum is ~2e-16; ALIGN_TOL is 1e-14)" % worst)
print("baseline entries outside tolerance: %d/12" % nbad)

print()
print("SELF-TEST (the check must be able to say no): same comparison against a")
print("baseline poisoned the way the original L1 bug poisoned it (sum to k=n-1):")
L1_bad = [mp.mpf(0)] * N
acc = mp.mpf(0)
for nn in range(1, N):
    L1_bad[nn] = acc / mp.mpf(2 * nn + 1)      # excludes k=nn  <- the old bug
    acc += 2 * Fh[nn]
for nm, idx in (("L1[1]", 1), ("L1[20]", 20), ("L1[399]", 399)):
    claimed = dict(zip([c[0] for c in CLAIMED], [c[1] for c in CLAIMED]))[nm]
    # denominator is the BASELINE (nonzero): that is what the alignment check
    # scales its pad by.  L1_bad[1] is legitimately 0 under the poisoned formula.
    base = mp.mpf(float(claimed))
    rel = float(abs(base - L1_bad[idx]) / abs(base))
    print("  %-8s vs poisoned truth: rel diff %.2e  -> %s" % (
        nm, rel, "DETECTED" if rel > 1e-14 else "*** MISSED ***"))
