"""Can any simple explicit family reach residual/gap <= 0.1 on w in [0.3,4]?  (numerical-evidence)"""
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from scipy.optimize import minimize
from lemmaC_compact_probe import O_of, top, THETA_STAR
N=600; n=np.arange(N); u=1/np.sqrt(2*n+1)
fams={
 'lor(a,b): u/(1+(n/a)^2)^b': lambda p: u/(1+(n/abs(p[0]))**2)**abs(p[1]),
 'lor2(a,b,c): u/(1+((n+c)/a)^2)^b': lambda p: u/(1+((n+p[2])/abs(p[0]))**2)**abs(p[1]),
 'exp(a,b): u*exp(-(n/a)^b)': lambda p: u*np.exp(-(n/abs(p[0]))**abs(p[1])),
 'e0+lor(a,b,c): c*e0 + u/(1+(n/a)^2)^b': lambda p: (np.eye(N)[0]*p[2] + u/(1+(n/abs(p[0]))**2)**abs(p[1])),
}
print(f"{'w':>5} {'h':>7} {'g_an':>7} | " + " | ".join(f"{k.split(':')[0]:>16} res  s_bnd" for k in fams))
for w in [0.3,0.5,0.7,1.0,1.5,2.0,3.0,4.0]:
    x=w*w; A=O_of(lambda k: x/(k*k+x),N); h,v,l2,_,_=top(A)
    g_an = 0.4*h if w<=2 else h-THETA_STAR
    out=[]
    for name,fam in fams.items():
        def res(p):
            t=fam(p); t=t/np.linalg.norm(t); r=A@t; return np.linalg.norm(r-(t@r)*t)
        best=np.inf
        for p0 in ([w,1,0.5],[2*w,2,0.0],[w,0.5,0.2],[0.5*w,1,0.0]):
            k=len(p0) if 'lor2' in name or 'e0' in name else 2
            try:
                o=minimize(res,p0[:k],method='Nelder-Mead',options={'xatol':1e-5,'fatol':1e-9,'maxiter':2000})
                best=min(best,o.fun)
            except Exception: pass
        out.append(f"{best:8.4f} {best/g_an:6.3f}")
    print(f"{w:5.2f} {h:7.4f} {g_an:7.4f} | " + " | ".join(f"{'':>16}{o}" for o in out))
