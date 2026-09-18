"""EXTERNAL high-precision check of h3s_deficit_cert_v4's intervals.

The v4 alignment check compares its intervals against a FLOAT64 baseline.  Now that
F3 is fixed the intervals are ~1e-51 wide, so a float64 baseline (own error ~1e-16)
CANNOT sit inside them -- the reported FAILs may be an artifact of that, not unsound
arithmetic.  This script settles it by recomputing the same quantities in POINT
arithmetic at dps=120, from the derivation definitions (same text as
_main_indep_baseline.py), and asking whether v4's intervals bracket THAT.

Definitions (P4_main_theorems.md:17):
    p_i = (2/lam)(a2F_i/w_i) * trapz_i,  normalised to sum 1
    x_i = (w_i/(2 pi T))^2 ,  F(k) = sum_i p_i x_i/(k^2+x_i) ,  F(0):=1
    u_n = (2n+1)^{-1/2} ,  L1[n] = (2n+1)^{-1} sum_{k=1}^{n} 2F(k)   (INCLUDES k=n)
    A[n,m] = (F(|n-m|)+F(n+m+1)) u_n u_m  for n!=m ;  A[n,n] = -L1[n] + F(2n+1) u_n^2
"""
import os, sys, glob
os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_cert_v4 import build_M_iv
from h3s_fill_factor import moments, ROOT

DPS_HI, DPS_V4, N, T = 120, 50, 40, 1.0

g = [x for x in glob.glob(f"{ROOT}/H3S-200GPa/EPW_anharmonic/*.a2f")
     if not x.endswith(".a2fv")]
lam_f, om2_f, olog_f, w, a2f = moments(g[0])

# ---------------- independent high-precision point computation -------------
mp.mp.dps = DPS_HI
# NB mp.mpf(float(x)) converts the EXACT binary float64; mp.mpf(repr(x)) would
# instead take the shortest round-trip DECIMAL string, shifting each datum by
# ~1e-17 relative -- enough to push a 1e-49-wide certificate interval off.
W = [mp.mpf(float(x)) for x in w]                # data points taken as exact
Y = [mp.mpf(float(x)) for x in a2f]
n_pts = len(W)
tw = [None] * n_pts
tw[0] = (W[1] - W[0]) / 2
tw[-1] = (W[-1] - W[-2]) / 2
for i in range(1, n_pts - 1):
    tw[i] = (W[i + 1] - W[i - 1]) / 2
raw = [(Y[i] / W[i]) * tw[i] for i in range(n_pts)]
tot = mp.fsum(raw)
p = [raw[i] / tot for i in range(n_pts)]          # sums to 1
xs = [(W[i] / (2 * mp.pi * mp.mpf(float(T)))) ** 2 for i in range(n_pts)]

def F_hi(k):
    if k == 0:
        return mp.mpf(1)
    k2 = mp.mpf(k) ** 2
    return mp.fsum([p[i] * xs[i] / (k2 + xs[i]) for i in range(n_pts)])

Fh = [F_hi(k) for k in range(2 * N + 2)]
u_hi = [mp.mpf(1) / mp.sqrt(mp.mpf(2 * n + 1)) for n in range(N)]
L1_hi = [mp.mpf(0)] * N
acc = mp.mpf(0)
for nn in range(1, N):
    acc += 2 * Fh[nn]
    L1_hi[nn] = acc / mp.mpf(2 * nn + 1)

def A_hi(i, j):
    if i == j:
        return -L1_hi[i] + Fh[2 * i + 1] * u_hi[i] ** 2
    return (Fh[abs(i - j)] + Fh[i + j + 1]) * u_hi[i] * u_hi[j]

# ---------------- v4's intervals ------------------------------------------
M_iv, u_iv, fker_iv, L1_iv, _, _ = build_M_iv(
    w, a2f, lam_f, 0.0, 10.0, T, N, DPS_V4, Nc_fixed=None, poison=False)

BASE64 = {                       # the float64 numbers v4 compares against
    ('F', 1): 9.8926738719898688e-01, ('F', 20): 4.7083503166965779e-01,
    ('F', 80): 7.7154780791224328e-02,
    ('L1', 1): 6.5951159146599125e-01, ('L1', 20): 6.9152796470690037e-01,
    ('L1', 39): 5.1792707943366445e-01,
    ('A', (0, 0)): 9.8926738719898688e-01, ('A', (0, 1)): 1.1290937977616782e+00,
    ('A', (20, 20)): -6.8615029568296637e-01,
    ('A', (20, 21)): 2.8643164353339162e-02,
}

def relw(iv):
    """correct relative width: subtract in mpf, THEN convert (v4:170 converts to
    float64 first, which collapses a 1e-51 width to exactly 0.0).
    NB iv.a / iv.b are themselves zero-width INTERVALS in mpmath, so they must be
    pulled into mpf before arithmetic or the subtraction re-widens."""
    a, b = mp.mpf(iv.a), mp.mpf(iv.b)
    d, m = b - a, (a + b) / 2
    return float(d / abs(m)) if m != 0 else float(d)

print(f"{'quantity':<12}{'relwidth':>11}{'encl TRUE':>11}{'encl f64':>10}"
      f"{'f64 rel-err vs true':>21}")
print("-" * 66)
rows = ([('F', k, fker_iv[k], Fh[k]) for k in (1, 20, 80)]
        + [('L1', k, L1_iv[k], L1_hi[k]) for k in (1, 20, 39)]
        + [('A', ij, None, A_hi(*ij)) for ij in ((0, 0), (0, 1), (20, 20), (20, 21))])
n_true_fail = 0
for nm, key, iv, truth in rows:
    if iv is None:                       # rebuild A element from v4's pieces
        i, j = key
        iv = (-L1_iv[i] + fker_iv[2 * i + 1] * u_iv[i] ** 2 if i == j else
              (fker_iv[abs(i - j)] + fker_iv[i + j + 1]) * u_iv[i] * u_iv[j])
    mp.mp.dps = DPS_HI                   # widen BEFORE pulling endpoints to mpf
    a, b = mp.mpf(iv.a), mp.mpf(iv.b)
    enc_true = a <= truth <= b
    n_true_fail += (not enc_true)
    f64 = BASE64[(nm, key)]
    enc_f64 = float(a) <= f64 <= float(b)
    err64 = float(abs(mp.mpf(f64) - truth) / abs(truth))
    print(f"{nm+str(key):<12}{relw(iv):11.2e}"
          f"{'YES' if enc_true else '*** NO ***':>11}"
          f"{'PASS' if enc_f64 else 'FAIL':>10}{err64:21.3e}")

print("-" * 66)
print(f"v4 intervals failing to enclose the TRUE (dps=120) value: {n_true_fail}")
print("=> if 0, every 'FAIL' in cert_v4_N40_normal.log is the float64 baseline's\n"
      "   own ~1e-16 error, not unsound interval arithmetic.")

# ---- does a tolerance-based alignment check still catch the poison? -------
_, _, fkerP, L1P, _, _ = build_M_iv(w, a2f, lam_f, 0.0, 10.0, T, N, DPS_V4,
                                    Nc_fixed=None, poison=True)
print("\npoison separation vs a 1e-14 relative tolerance:")
for k in (1, 20, 39):
    good, bad = mp.mpf(L1_iv[k].a), mp.mpf(L1P[k].a)
    sep = float(abs(bad - good) / abs(good))
    print(f"  L1[{k:2d}]: |poison-good|/good = {sep:.3e}  "
          f"{'DETECTED' if sep > 1e-14 else '*** MISSED ***'}")

# ---- SELF-TEST of this enclosure check: it must be able to say NO ---------
# v3 has F2 (dps inert, prec=53) and F3 (point-mpf injection) UNFIXED.  If the
# check above cannot fail, it proves nothing about v4.  Expect FAILs here.
from h3s_deficit_cert_v3 import build_M_iv as build_v3
_, u3, fker3, L13 = build_v3(w, a2f, lam_f, 0.0, 10.0, T, N, DPS_V4,
                             nc_mode='floor', poison=False)
mp.mp.dps = DPS_HI
n_v3_fail = 0
print("\nSELF-TEST -- same check applied to v3 (F2/F3 unfixed):")
for nm, key, iv, truth in ([('F', k, fker3[k], Fh[k]) for k in (1, 20, 80)]
                           + [('L1', k, L13[k], L1_hi[k]) for k in (1, 20, 39)]):
    mp.mp.dps = DPS_HI
    ok = mp.mpf(iv.a) <= truth <= mp.mpf(iv.b)
    n_v3_fail += (not ok)
    print(f"  {nm+str(key):<10} relwidth={relw(iv):9.2e}  "
          f"encloses TRUE = {'YES' if ok else '*** NO ***'}")
print(f"  v3 failures: {n_v3_fail}/6  "
      f"{'-> check HAS discriminating power' if n_v3_fail else '*** CHECK IS INERT ***'}")
