import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from scipy.optimize import minimize
from lemmaC_compact_probe import O_of, top, THETA_STAR
N=600; n=np.arange(N); u=1/np.sqrt(2*n+1)
def fam2(p): return u/(1+(n/abs(p[0]))**2)**abs(p[1])
def fam3(p): return u/(1+((n+p[2])/abs(p[0]))**2)**abs(p[1])
print(f"{'w':>5} | lor(a,b): a/w  b  res  s_true | lor2(a,b,c): a/w  b  c  res  s_true")
for w in [0.3,0.5,0.7,1.0,1.5,2.0,3.0,4.0,6.0,8.5]:
    x=w*w; A=O_of(lambda k: x/(k*k+x),N); h,v,l2,_,_=top(A)
    def res(fam,p):
        t=fam(p); t=t/np.linalg.norm(t); r=A@t; return np.linalg.norm(r-(t@r)*t)
    o2=min((minimize(lambda p: res(fam2,p),p0,method='Nelder-Mead',options={'xatol':1e-6,'fatol':1e-10,'maxiter':4000}) for p0 in ([w,1],[2*w,2],[0.7*w,0.7])),key=lambda o:o.fun)
    o3=min((minimize(lambda p: res(fam3,p),p0,method='Nelder-Mead',options={'xatol':1e-6,'fatol':1e-10,'maxiter':6000}) for p0 in ([w,1,0.5],[2*w,2,0.0],[w,0.8,0.3])),key=lambda o:o.fun)
    t2=fam2(o2.x); t2/=np.linalg.norm(t2); t3=fam3(o3.x); t3/=np.linalg.norm(t3)
    s2=np.sqrt(max(0,1-(t2@v)**2)); s3=np.sqrt(max(0,1-(t3@v)**2))
    print(f"{w:5.2f} | {abs(o2.x[0])/w:6.3f} {abs(o2.x[1]):6.3f} {o2.fun:8.5f} {s2:7.4f} | {abs(o3.x[0])/w:6.3f} {abs(o3.x[1]):6.3f} {o3.x[2]:7.3f} {o3.fun:8.5f} {s3:7.4f}")
