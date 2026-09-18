"""Does float(iv.a) / float(iv.b) round OUTWARD, or to-nearest?

v4 converts interval endpoints to float64 at ~20 sites and then uses each
one-sidedly (hi endpoints as upper bounds, lo endpoints as lower bounds,
both sides of the Temple premise rho > beta).  If float() rounds to nearest,
every one of those is a ~1e-16 one-sided hole in the certificate.
"""
import math, mpmath as mp
mp.mp.dps = 50; mp.iv.dps = 50

n_lo_bad = n_hi_bad = 0
worst_lo = worst_hi = 0.0
for m in range(3, 2000):
    x = mp.mpi(1) / mp.mpi(m)          # endpoints need >53 bits
    a, b = mp.mpf(x.a), mp.mpf(x.b)
    fa, fb = float(x.a), float(x.b)    # exactly what v4 does
    mp.mp.dps = 200
    if mp.mpf(fa) > a:                 # lower endpoint moved UP => hole
        n_lo_bad += 1
        worst_lo = max(worst_lo, float((mp.mpf(fa) - a) / a))
    if mp.mpf(fb) < b:                 # upper endpoint moved DOWN => hole
        n_hi_bad += 1
        worst_hi = max(worst_hi, float((b - mp.mpf(fb)) / b))
    mp.mp.dps = 50

print(f"m = 3..1999  ({1997} intervals)")
print(f"  float(iv.a) moved UP   (lower-bound hole): {n_lo_bad}  "
      f"worst rel {worst_lo:.2e}")
print(f"  float(iv.b) moved DOWN (upper-bound hole): {n_hi_bad}  "
      f"worst rel {worst_hi:.2e}")
print("  => float() rounds " + ("TO NEAREST -- one-sided handoffs are holes"
      if (n_lo_bad or n_hi_bad) else "OUTWARD -- handoffs are safe"))

# the correct outward handoff
def lo_f(iv):
    x = float(mp.mpf(iv.a))
    return x if mp.mpf(x) <= mp.mpf(iv.a) else math.nextafter(x, -math.inf)
def hi_f(iv):
    x = float(mp.mpf(iv.b))
    return x if mp.mpf(x) >= mp.mpf(iv.b) else math.nextafter(x, math.inf)

bad = 0
for m in range(3, 2000):
    x = mp.mpi(1) / mp.mpi(m)
    mp.mp.dps = 200
    if not (mp.mpf(lo_f(x)) <= mp.mpf(x.a) and mp.mpf(hi_f(x)) >= mp.mpf(x.b)):
        bad += 1
    mp.mp.dps = 50
print(f"  proposed lo_f/hi_f (nextafter): {bad} failures")
