import os, sys, math, glob, time
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_deficit_budget import Tc_P, Q_P
from h3s_fill_factor import moments, MEV_K, ROOT

lam,om2,_,w,a2f = moments(glob.glob(f"{ROOT}/H3S-200GPa/EPW_anharmonic/*.a2f")[0])
sq=math.sqrt(om2)
def mu_run(r, mu0=0.13, r0=10.0):
    inv=1.0/mu0-math.log(r/r0); return 1.0/inv if inv>0 else float("nan")

# fast local bisection seeded near the known plateau (avoids the 160-point geomspace)
def Tc_fast(ms, r, N, lo=14.0, hi=20.0):
    f=lambda T: lam*Q_P(T,w,a2f,lam,ms,float(r),N,sq)-1.0
    if f(lo)<0 or f(hi)>0: return None
    for _ in range(44):
        m=0.5*(lo+hi)
        if f(m)>=0: lo=m
        else: hi=m
    return lo

print(f"lam={lam:.4f}  sqrt<om2>={sq:.3f} meV")
print(f"{'r':>6} {'mu*(r)':>8} | " + " ".join(f"{'N='+str(N):>26}" for N in (400,800,1600)))
for r in (10,20,40,80,160):
    ms=mu_run(r); row=f"{r:>6} {ms:>8.4f} | "
    for N in (400,800,1600):
        t0=time.time(); t=Tc_fast(ms,r,N)
        if t is None: row+=f"{'out of bracket':>26} "; continue
        vp=sq/(2*math.pi*t); Nc=int(math.floor(r*vp))
        fr=Nc/N
        flag="DEGEN" if (Nc<=0 or Nc>=N) else ("Nc>N/2!" if fr>0.5 else "")
        row+=f"{t*MEV_K:>8.2f}K Nc={Nc:<4d}{fr:>5.2f}N {flag:>7} "
    print(row, flush=True)
