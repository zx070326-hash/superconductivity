"""Is v4 line 234's  om2**(mpi(1)/mpi(2))  outward-rounded?  And line 212/235's
mp.iv.pi path?  A non-enclosing varpi could in principle move Nc_lo/Nc_hi.
"""
import mpmath as mp
DPS = 50
mp.mp.dps = DPS; mp.iv.dps = DPS
bad_pow = []
for m in range(1, 400):
    x = mp.mpi(m) + mp.mpi(1)/mp.mpi(7)          # a generic non-representable value
    s_pow = x**(mp.mpi(1)/mp.mpi(2))
    s_iv  = mp.iv.sqrt(x)
    mp.mp.dps = 200
    lo = mp.sqrt(mp.mpf(m) + mp.mpf(1)/mp.mpf(7))
    mp.mp.dps = DPS
    a, b = mp.mpf(s_pow.a), mp.mpf(s_pow.b)
    ai, bi = mp.mpf(s_iv.a), mp.mpf(s_iv.b)
    mp.mp.dps = 200
    # the truth for THIS x is an interval; test against x's own midpoint truth
    # is not meaningful -- instead demand s^2 encloses x (self-consistency)
    if not (a <= lo <= b):
        bad_pow.append(m)
    mp.mp.dps = DPS
print(f"interval-exponent sqrt, m=1..399: {len(bad_pow)} non-enclosing"
      f"{' -> ' + str(bad_pow[:8]) if bad_pow else ''}")

# squaring self-consistency: does s**2 contain x?  (soundness w/o external truth)
bad_sq = 0
for m in range(1, 400):
    x = mp.mpi(m) + mp.mpi(1)/mp.mpi(7)
    s = x**(mp.mpi(1)/mp.mpi(2))
    if not (s**2).a <= x.mid <= (s**2).b:
        bad_sq += 1
print(f"  (s^2 fails to contain x): {bad_sq}")
