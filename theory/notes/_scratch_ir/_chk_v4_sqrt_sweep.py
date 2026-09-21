"""Is v4 line 130's  mp.mpi(1)/mp.sqrt(mp.mpi(2n+1))  sound for EVERY n?

mp.sqrt applied to an interval returns a POINT mpf (verified in _chk_v4_widths),
so the square root's own rounding error is NOT outward-rounded.  The following
division by that point re-widens by ~2e-51, which may or may not cover the
missing sqrt rounding.  Sweep all n < 800 against dps=200 truth.
Also compare with the correct idiom mp.iv.sqrt.
"""
import mpmath as mp

DPS, NMAX = 50, 800
mp.mp.dps = DPS
mp.iv.dps = DPS

bad_v4, bad_iv, wid_v4 = [], [], []
for n in range(NMAX):
    m = 2 * n + 1
    u_v4 = mp.mpi(1) / mp.sqrt(mp.mpi(m))        # v4 line 130
    u_ok = mp.mpi(1) / mp.iv.sqrt(mp.mpi(m))     # correct idiom
    mp.mp.dps = 200
    truth = mp.mpf(1) / mp.sqrt(mp.mpf(m))
    mp.mp.dps = DPS
    a4, b4 = mp.mpf(u_v4.a), mp.mpf(u_v4.b)
    ai, bi = mp.mpf(u_ok.a), mp.mpf(u_ok.b)
    mp.mp.dps = 200
    if not (a4 <= truth <= b4):
        bad_v4.append((n, float(min(abs(truth - a4), abs(truth - b4)) / truth)))
    if not (ai <= truth <= bi):
        bad_iv.append(n)
    wid_v4.append(float((b4 - a4) / truth))
    mp.mp.dps = DPS

print(f"n = 0..{NMAX-1}, dps={DPS}")
print(f"  v4 idiom   mp.mpi(1)/mp.sqrt(mp.mpi(m)) : "
      f"{len(bad_v4)} non-enclosing")
if bad_v4:
    print("    first 10 misses (n, rel-miss):", bad_v4[:10])
print(f"  correct    mp.mpi(1)/mp.iv.sqrt(mp.mpi(m)) : "
      f"{len(bad_iv)} non-enclosing")
print(f"  v4 relwidth: min {min(wid_v4):.2e}  max {max(wid_v4):.2e}")
print("\nverdict: " + ("v4's u_n is sound on this whole range (the division's "
      "re-widening happens to\n  cover the un-rounded sqrt) -- but that is "
      "incidental, not guaranteed; still\n  should be mp.iv.sqrt."
      if not bad_v4 else "*** v4's u_n is NOT sound -- F3 residual is a real "
      "defect, not a nit ***"))
