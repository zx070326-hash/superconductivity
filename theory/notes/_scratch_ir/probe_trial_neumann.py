"""Why is tau's overlap with the top vector so low, and can one Neumann step fix it?
(M1): H = 2 tau tau^T + R - D.  Secular: v propto (h + D - R)^{-1} tau.
Candidates: t0=tau, t1=(h+D)^{-1}tau, t2=(h+D-R)^{-1}tau (exact, needs R), 
            t3 = two-term Neumann (h+D)^{-1}(I + R (h+D)^{-1}) tau.
Report overlaps.  Also commutator norm ||[O[f],O[f^2]]|| (route 3 of the plan).
"""
import numpy as np, sys
sys.path.insert(0,r"D:/physics/superconductivity/theory/derivations")
from p4_interchange_uniformity import H_mat
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
def ov(a,b): return abs(a@b)/(np.linalg.norm(a)*np.linalg.norm(b))
print(f"{'varpi':>6} {'tau':>7} {'(h+D)^-1tau':>11} {'Neu2':>7} {'exact':>7} {'|[Of,Of2]|':>10} {'|Of2|':>7}")
for vp in [0.2,0.5,1.0,2.0,3.0,5.0,10.0,20.0]:
    x=vp*vp; N=max(200,int(40*vp)+200); K=2*N+2
    k=np.arange(1,K+1,dtype=float); f=x/(k**2+x)
    A=O_of_kernel(N,f); w,V=np.linalg.eigh(A); i=np.argmax(w); h=w[i]; v=V[:,i]
    if v.sum()<0: v=-v
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); p=u/np.linalg.norm(u)
    s=(n+0.5)/vp; G=s**2/(1+s**2); tau=p*G
    # (M1): D = diag, R = A - 2 tau tau^T + D  => extract D from diagonal identity
    M = A - 2.0*np.outer(tau,tau)      # = R - D
    Dg = -np.diag(M).copy()            # D_nn = -(R-D)_nn since R has zero diagonal? check
    R  = M + np.diag(Dg)
    offdiag_ok = abs(np.diag(R)).max()
    t1=np.linalg.solve(np.diag(h+Dg), tau)
    t2=np.linalg.solve((h*np.eye(N)+np.diag(Dg)-R), tau)
    t3=np.linalg.solve(np.diag(h+Dg), tau + R@t1)
    Af2=O_of_kernel(N,f**2)
    comm=np.linalg.norm(A@Af2-Af2@A,2)
    print(f"{vp:6.2f} {ov(tau,v):7.4f} {ov(t1,v):11.4f} {ov(t3,v):7.4f} {ov(t2,v):7.4f} {comm:10.5f} {np.linalg.norm(Af2,2):7.4f}")
