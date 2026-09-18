"""Why does v4 report relwidth == 0 everywhere?  Main session, direct probe.

Tests, in order of suspicion:
  (1) does mp.mp.dps / mp.iv.dps actually raise iv precision now?
  (2) does mp.sqrt(mp.mpi(x))  -- the mp-CONTEXT sqrt applied to an iv value,
      which is what v4 line 130 does -- return a widening interval, or does it
      silently fall back to point arithmetic?
  (3) does x**(mp.mpi(1)/mp.mpi(2))  -- v4's other sqrt idiom -- widen?
  (4) does mp.iv.pi enclose pi?
Ground truth for each: recompute at dps=200 in the POINT context and check
whether the claimed interval brackets it.
"""
import mpmath as mp

DPS = 50
mp.mp.dps = DPS
mp.iv.dps = DPS

def w(x):
    """relative width of an iv.mpf, as v4 would report it"""
    try:
        a, b = mp.mpf(x.a), mp.mpf(x.b)
    except AttributeError:
        return None
    if b == 0: return 0.0
    return float((b - a) / abs(b))

print(f"mp.mp.prec = {mp.mp.prec}   mp.iv.prec = {mp.iv.prec}"
      f"   (53 would mean dps is STILL inert)")
print(f"mp.mp.dps  = {mp.mp.dps}    mp.iv.dps  = {mp.iv.dps}\n")

# ---- (1) does a plain irrational interval widen at all? -------------------
third = mp.mpi(1) / mp.mpi(3)
print(f"(1) mpi(1)/mpi(3)          relwidth = {w(third):.3e}   "
      f"{'OK (widens)' if w(third) > 0 else '*** ZERO ***'}")

# ---- (2) v4 line 130 idiom:  mp.sqrt on an interval ----------------------
n = 20
s_ctx = mp.sqrt(mp.mpi(2 * n + 1))          # <-- exactly what v4 does
u_v4 = mp.mpi(1) / s_ctx
print(f"(2) mp.sqrt(mpi(41))       type = {type(s_ctx).__name__}, "
      f"relwidth = {w(s_ctx)}")
print(f"    u = mpi(1)/that        relwidth = {w(u_v4)}")

# ---- (3) the other idiom: interval exponent ------------------------------
s_pow = mp.mpi(2 * n + 1) ** (mp.mpi(1) / mp.mpi(2))
u_pow = mp.mpi(1) / s_pow
print(f"(3) mpi(41)**(mpi(1)/mpi(2))  relwidth = {w(s_pow)}")
print(f"    u = mpi(1)/that           relwidth = {w(u_pow)}")

# ---- (3b) the CORRECT idiom: iv context sqrt -----------------------------
s_iv = mp.iv.sqrt(mp.mpi(2 * n + 1))
u_iv = mp.mpi(1) / s_iv
print(f"(3b) mp.iv.sqrt(mpi(41))      relwidth = {w(s_iv)}")
print(f"     u = mpi(1)/that          relwidth = {w(u_iv)}")

# ---- ground truth: does each candidate ENCLOSE the true 1/sqrt(41)? ------
mp.mp.dps = 200
truth = mp.mpf(1) / mp.sqrt(mp.mpf(41))
mp.mp.dps = DPS
print(f"\n    true 1/sqrt(41) (dps=200) = {mp.nstr(truth, 30)}")
for nm, x in [("v4 mp.sqrt", u_v4), ("interval-exponent", u_pow),
              ("iv.sqrt", u_iv)]:
    try:
        ok = mp.mpf(x.a) <= truth <= mp.mpf(x.b)
        off = 0.0 if ok else float(min(abs(truth - mp.mpf(x.a)),
                                       abs(truth - mp.mpf(x.b))) / truth)
        print(f"    {nm:<20} encloses={'YES' if ok else '*** NO ***':<12}"
              f" rel-miss={off:.3e}  relwidth={w(x)}")
    except Exception as e:
        print(f"    {nm:<20} ERROR {e}")

# ---- (4) mp.iv.pi --------------------------------------------------------
p = mp.iv.pi
mp.mp.dps = 200
pi_true = +mp.pi
mp.mp.dps = DPS
try:
    ok = mp.mpf(p.a) <= pi_true <= mp.mpf(p.b)
    print(f"\n(4) mp.iv.pi  relwidth = {w(p)}  encloses true pi = "
          f"{'YES' if ok else '*** NO ***'}")
except AttributeError:
    print(f"\n(4) mp.iv.pi is {type(p).__name__}, not an interval: {p}")
