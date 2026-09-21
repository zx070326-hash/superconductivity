import os, sys, math, glob
os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import h3s_deficit_budget as H
from h3s_fill_factor import moments, MEV_K, ROOT

good = H.O_of
def bad(fker, N):                       # cert script's L1 (sum_{k=1}^{n-1})
    n=np.arange(N); u=1.0/np.sqrt(2.0*n+1.0)
    absd=np.abs(n[:,None]-n[None,:]); sumd=n[:,None]+n[None,:]+1
    uu=u[:,None]*u[None,:]
    Am=(fker[absd]+fker[sumd])*uu
    L1=np.zeros(N)
    for m in range(1,N): L1[m]=np.sum(2.0*fker[1:m])/(2*m+1)
    np.fill_diagonal(Am,-L1+fker[2*n+1]*u**2)
    return n,u,Am

lam,om2,olog,w,a2f = moments(glob.glob(f"{ROOT}/H3S-200GPa/EPW_anharmonic/*.a2f")[0])
sq=math.sqrt(om2)
for name,fn in (("baseline (KAY)",good),("cert script L1",bad)):
    H.O_of=fn
    t0=H.Tc_P(w,a2f,lam,0.0,10.0,400,sq); t1=H.Tc_P(w,a2f,lam,0.13,10.0,400,sq)
    print(f"{name:>16}:  Tc(mu*=0)={t0*MEV_K:8.3f} K   Tc(0.13)={t1*MEV_K:8.3f} K")
