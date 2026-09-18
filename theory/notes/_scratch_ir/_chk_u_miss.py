import mpmath as mp
mp.mp.dps = 15; mp.iv.dps = 15
rows=[]
for n in (0,1,5,20,399):
    u_v3 = mp.mpi(mp.mpf(1)/mp.sqrt(mp.mpf(2*n+1)))   # v3 line 114: sqrt+div in POINT mpf
    a,b = mp.mpf(u_v3.a), mp.mpf(u_v3.b)
    mp.mp.dps = 60
    ref = mp.mpf(1)/mp.sqrt(mp.mpf(2*n+1))
    inside = (a <= ref <= b)
    rel = float(abs(ref-a)/ref)
    mp.mp.dps = 15
    rows.append((n, inside, rel))
print(f"{'n':>5} {'v3 interval encloses true u_n?':>32} {'|miss|/u_n':>14}")
for n,ins,rel in rows:
    print(f"{n:>5} {str(ins):>32} {rel:>14.3e}")
