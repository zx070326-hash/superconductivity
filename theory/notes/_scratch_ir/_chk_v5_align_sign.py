"""Does v5's F6 alignment predicate work for NEGATIVE baselines?

v5 (lines ~727, 750, 774):
    enclosed = (a <= base*(1+TOL)) and (base*(1-TOL) <= b)

For base < 0, base*(1+TOL) < base*(1-TOL): the band is INVERTED, so the
predicate demands a <= (more negative) and b >= (less negative) -- i.e. it
requires the interval to be WIDER than 2*TOL*|base|, the opposite of a tolerance.
Among the 12 aligned quantities exactly one is negative: A[20,20] = -6.86e-01,
and it is one of the three entries that carry the poison signal.
"""
import mpmath as mp
mp.mp.dps = 50; mp.iv.dps = 50
TOL = 1e-14

def v5_pred(a, b, base):
    return (a <= base * (1 + TOL)) and (base * (1 - TOL) <= b)

def fixed_pred(a, b, base):           # sign-safe additive band
    pad = TOL * abs(base)
    return (a <= base + pad) and (base - pad <= b)

# realistic case: certificate interval ~1e-49 wide sitting within 1e-16 of base
cases = [
    ("F[1]      (positive)",  9.8926738719898688e-01),
    ("L1[20]    (positive)",  6.9152796470690037e-01),
    ("A[20,20]  (NEGATIVE)", -6.8615029568296637e-01),
    ("A[20,21]  (positive)",  2.8643164353339162e-02),
]
print(f"TOL = {TOL:.0e};  interval = base*(1+3e-16) +/- 4e-49  (a sound, tiny certificate interval)")
print(f"{'quantity':<24}{'v5 predicate':>14}{'sign-safe fix':>16}")
print("-" * 54)
for name, base in cases:
    ctr = base * (1 + 3e-16)          # v4-scale offset from the float64 baseline
    a, b = ctr - 4e-49, ctr + 4e-49
    print(f"{name:<24}{'PASS' if v5_pred(a,b,base) else '*** FAIL ***':>14}"
          f"{'PASS' if fixed_pred(a,b,base) else '*** FAIL ***':>16}")

print("\nwhat v5's predicate actually demands for base < 0:")
base = -6.8615029568296637e-01
print(f"  base*(1+TOL) = {base*(1+TOL):.17e}   (a must be <= this)")
print(f"  base*(1-TOL) = {base*(1-TOL):.17e}   (b must be >= this)")
print(f"  => required interval width >= {base*(1-TOL) - base*(1+TOL):.3e}")
print(f"     v5's actual width       ~  8e-49   => IMPOSSIBLE, always FAIL")

# does the fix ALSO still reject a real defect (teeth)?
print("\nteeth of the sign-safe fix (poison-scale error, ~1e-2 relative):")
for name, base in cases:
    ctr = base * (1 + 3.3e-2)
    a, b = ctr - 4e-49, ctr + 4e-49
    print(f"  {name:<24}{'*** MISSED ***' if fixed_pred(a,b,base) else 'DETECTED':>16}")
