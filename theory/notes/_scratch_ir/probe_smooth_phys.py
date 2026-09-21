"""Is the smooth-cutoff dip confined to the unphysical regime?
lam1(M)=1/lambda at Tc, so lam1=3e-4 <=> lambda~3400. Physical: lambda in [0.3,4] <=> lam1 in [0.25,3.3].
Scan the whole w range and report lam1 together with implied lambda, marking dips.
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
N=420; K=2*N+2; kk=np.arange(1,K+1,dtype=float); nn=np.arange(N); sq=np.sqrt(2*nn+1.0)
def lam1(w,r,mu,p):
    x=w*w; A=O_of_kernel(N,x/(kk**2+x))
    t=(2*nn+1.0)/(2*r*w); u=(1.0/(1.0+t**(2*p)))/sq
    return np.linalg.eigvalsh(A-mu*2*np.outer(u,u)).max()
print("r=4 mu=0.26 p=2 : lam1 and implied lambda=1/lam1")
print(f"{'w':>5} {'lam1':>13} {'lambda':>11} {'phys?':>6}")
prev=None
for w in [0.4,0.5,0.6,0.8,1.0,1.2,1.5,2.0,3.0,5.0]:
    L=lam1(w,4,0.26,2); lam=1.0/L if L>0 else np.inf
    mark="YES" if 0.2<=lam<=6 else "no"
    d="" if prev is None else ("  DOWN" if L<prev else "")
    print(f"{w:5.2f} {L:13.6e} {lam:11.1f} {mark:>6}{d}")
    prev=L
print()
print("Where does lam1 enter the physical band (lambda<=6, i.e. lam1>=1/6)?")
for w in [0.9,0.95,1.0,1.05,1.1,1.2]:
    L=lam1(w,4,0.26,2); print(f"   w={w:.2f}  lam1={L:.6f}  lambda={1.0/L:8.2f}")
