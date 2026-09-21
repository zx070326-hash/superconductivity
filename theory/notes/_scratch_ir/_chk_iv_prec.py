import mpmath as mp          # exactly as v3 imports it
print("module has own .dps?", 'dps' in mp.__dict__ or hasattr(type(mp), 'dps'))
for d in (15, 50, 200):
    mp.dps = d               # exactly as v3 does it
    x  = mp.mpi(1)/mp.mpi(3)
    y  = mp.mpf(1)/mp.mpf(3)
    print(f"  mp.dps={d:4d} -> mpmath.mp.prec={mp.mp.prec:5d}  iv.prec={mp.iv.prec:5d}"
          f"  iv width(1/3)={float(x.b)-float(x.a):.3e}"
          f"  mpf digits: {mp.nstr(y, 25)}")
print()
mp.mp.dps = 50; mp.iv.dps = 50
x = mp.mpi(1)/mp.mpi(3)
print(f"after setting mp.mp.dps=mp.iv.dps=50: iv width(1/3)={float(x.b)-float(x.a):.3e}"
      f"  mpf: {mp.nstr(mp.mpf(1)/mp.mpf(3), 25)}")
# u_n as v3 builds it vs interval-native
mp.mp.dps = 15; mp.iv.dps = 15
n = 20
u_v3   = mp.mpi(mp.mpf(1)/mp.sqrt(mp.mpf(2*n+1)))     # v3 line 114
u_true = mp.mpi(1)/mp.sqrt(mp.mpi(2*n+1))             # interval-native
mp.mp.dps = 60
ref = mp.mpf(1)/mp.sqrt(mp.mpf(2*n+1))
print(f"\nu_{n}: v3 width={float(u_v3.b)-float(u_v3.a):.3e}"
      f"  native width={float(u_true.b)-float(u_true.a):.3e}")
print(f"     v3 encloses true value? {mp.mpf(u_v3.a) <= ref <= mp.mpf(u_v3.b)}"
      f"   (miss = {float(ref)-float(mp.mpf(u_v3.a)):+.3e})")
