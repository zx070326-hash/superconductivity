"""Independent check of peer claim 1: smooth chi kills the downjumps of lam1(M) in w.
My own operator build (not importing nx_core). Finer w-grid than peer used.
M = O[f(.;w^2)] - mu*2 u u^T,  u_n = chi(t_n)/sqrt(2n+1), t_n=(2n+1)/(2 r w).
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
def lam1(w,r,mu,mode,p=2,N=None):
    N=N or int(min(1400,max(320,5*max(w,r*w)))); K=2*N+2
    k=np.arange(1,K+1,dtype=float); x=w*w
    A=O_of_kernel(N,x/(k**2+x))
    n=np.arange(N); t=(2*n+1.0)/(2*r*w)
    chi=(t<1.0).astype(float) if mode=='sharp' else 1.0/(1.0+t**(2*p))
    u=chi/np.sqrt(2*n+1.0)
    return np.linalg.eigvalsh(A-mu*2*np.outer(u,u)).max()
ws=np.arange(0.5,8.001,0.05)
print(f"{'model':>10} {'r':>3} {'mu':>5} {'#down':>6} {'worst':>11}")
for r in [3,4]:
    for mu in [0.10,0.13,0.26]:
        for mode,p in [('sharp',0),('smooth',2),('smooth',4)]:
            v=np.array([lam1(w,r,mu,mode,p) for w in ws]); d=np.diff(v)
            b=np.where(d<0)[0]
            lab=mode if mode=='sharp' else f"smooth p{p}"
            print(f"{lab:>10} {r:3d} {mu:5.2f} {len(b):6d} {(d[b].min() if len(b) else 0.0):11.2e}")
