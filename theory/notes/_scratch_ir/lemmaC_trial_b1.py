import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from lemmaC_compact_probe import O_of, top, THETA_STAR
N=600; n=np.arange(N); u=1/np.sqrt(2*n+1)
def famac(p): a,c=abs(p[0]),p[1]; return u*a*a/(a*a+(n+c)**2)
def fama(a,c): return u*a*a/(a*a+(n+c)**2)
print(f"{'w':>5} {'h':>7} {'g_an':>7} | b=1,(a,c) free: a/w  c  res  s_true  s_bnd | b=1,c=1/2: a/w  res  s_bnd | b=1,c=1/2,a=w (M1 tau): res  s_bnd")
for w in [0.3,0.5,0.7,1.0,1.5,2.0,3.0,4.0,6.0,8.5,12.0,20.0]:
    x=w*w; A=O_of(lambda k: x/(k*k+x),N); h,v,l2,_,_=top(A)
    g_an = 0.4*h if w<=2 else h-THETA_STAR
    def res_t(t):
        t=t/np.linalg.norm(t); r=A@t; return np.linalg.norm(r-(t@r)*t)
    o=min((minimize(lambda p: res_t(famac(p)),p0,method='Nelder-Mead',options={'xatol':1e-7,'fatol':1e-11,'maxiter':4000}) for p0 in ([1.5*w,0.5],[2*w,0.2],[w,0.5])),key=lambda o:o.fun)
    t=famac(o.x); t/=np.linalg.norm(t); s=np.sqrt(max(0,1-(t@v)**2))
    oa=minimize_scalar(lambda a: res_t(fama(abs(a),0.5)),bounds=(0.3*w,8*w),method="bounded")
    r_m1=res_t(fama(w,0.5))
    print(f"{w:5.2f} {h:7.4f} {g_an:7.4f} | {abs(o.x[0])/w:6.3f} {o.x[1]:6.3f} {o.fun:8.5f} {s:7.4f} {o.fun/g_an:6.4f} | {abs(oa.x)/w:6.3f} {oa.fun:8.5f} {oa.fun/g_an:6.4f} | {r_m1:8.5f} {r_m1/g_an:6.3f}")
