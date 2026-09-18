"""r=4,mu=0.26,p=2 dip near w=0.5 is N-STABLE => REAL. Characterise cheaply (N=420)."""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
N=420; K=2*N+2; kk=np.arange(1,K+1,dtype=float); nn=np.arange(N)
sq=np.sqrt(2*nn+1.0)
def lam1(w,r,mu,p):
    x=w*w; A=O_of_kernel(N,x/(kk**2+x))
    t=(2*nn+1.0)/(2*r*w); u=(1.0/(1.0+t**(2*p)))/sq
    return np.linalg.eigvalsh(A-mu*2*np.outer(u,u)).max()
ws=np.arange(0.40,0.681,0.02)
v=np.array([lam1(w,4,0.26,2) for w in ws]); d=np.diff(v)
print("fine grid r=4 mu=0.26 p=2 (N=420)")
for i in range(len(d)):
    print(f"   w {ws[i]:.2f}->{ws[i+1]:.2f}  lam1={v[i]: .9f}  d={d[i]: .3e}{'  <== DOWN' if d[i]<0 else ''}")
print("\n(r,p) sweep, w in [0.35,0.85] step 0.05, mu=0.26")
print(f"   {'r':>3} {'p':>2} {'#down':>6} {'worst':>11}")
wg=np.arange(0.35,0.851,0.05)
for r in [3,4,6,10]:
    for p in [2,4]:
        vv=np.array([lam1(w,r,0.26,p) for w in wg]); dd=np.diff(vv); bb=np.where(dd<0)[0]
        print(f"   {r:3d} {p:2d} {len(bb):6d} {(dd[bb].min() if len(bb) else 0.0):11.2e}")
