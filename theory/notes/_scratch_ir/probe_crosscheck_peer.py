"""Cross-check peer's split  H'' = t1 + t2  against my  H'' = -q'/x.
Peer: t1 = <v,O[d2f/dx2]v> <= 0 ; t2 = 2<w,(hI-M)^{-1}w>, w=(I-P_v)(dO/dx)v, M=O[f]-2 tau tau^T.
Mine: H'' = -(E+C)/x with E=<v,O[g]v>, C=-(2/x)sum_{j>1}<v,O[f^2]w_j>^2/(h-lam_j).
NOTE peer's t2 uses M (rank-one REMOVED) whereas the standard 2nd-order formula uses
the reduced resolvent of O[f] itself on the orthogonal complement. Test BOTH:
   t2_std = 2*sum_{j>1} <v,(dO/dx)w_j>^2/(h-lam_j)
   t2_peer = 2*<w,(hI-M)^{-1}w>
and check which reproduces H''. This is the load-bearing question for peer's claim.
d2f/dx2 = -2k^2/(k^2+x)^3 ; df/dx = k^2/(k^2+x)^2
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
print(f"{'varpi':>6} {'H2_num':>11} {'t1':>11} {'t2_std':>10} {'t1+t2std':>11} {'errS':>8} {'t2_peer':>10} {'t1+t2peer':>11} {'errP':>8} {'t2/|t1|':>8}")
for vp in [0.5,1.0,2.0,5.0,20.0,100.0]:
    x=vp*vp; N=max(260,int(40*vp)+260); K=2*N+2
    k=np.arange(1,K+1,dtype=float); f=x/(k**2+x)
    A=O_of_kernel(N,f)
    w_,V=np.linalg.eigh(A); o=np.argsort(w_)[::-1]; w_=w_[o]; V=V[:,o]
    h=w_[0]; v=V[:,0]
    d1=O_of_kernel(N,k**2/(k**2+x)**2)          # dO/dx
    d2=O_of_kernel(N,-2.0*k**2/(k**2+x)**3)     # d2O/dx2
    t1=v@(d2@v)
    # standard second-order: 2 sum_{j>1} <v,d1 w_j>^2/(h-lam_j)
    c1=V[:,1:].T@(d1@v)
    t2s=2.0*np.sum(c1**2/(h-w_[1:]))
    # peer's: M = A - 2 tau tau^T with tau from (M1)
    n=np.arange(N); uu=1.0/np.sqrt(2*n+1.0); p=uu/np.linalg.norm(uu)
    s=(n+0.5)/vp; G=s**2/(1+s**2); tau=p*G
    M=A-2.0*np.outer(tau,tau)
    r=d1@v; wv=r-(v@r)*v
    try:
        t2p=2.0*(wv@np.linalg.solve(h*np.eye(N)-M,wv))
    except Exception:
        t2p=float('nan')
    dd=1e-4*x
    def H(y):
        Ay=O_of_kernel(N,y/(k**2+y)); return np.linalg.eigvalsh(Ay).max()
    H2=(H(x+dd)-2*h+H(x-dd))/dd**2
    print(f"{vp:6.2f} {H2:11.3e} {t1:11.3e} {t2s:10.3e} {t1+t2s:11.3e} {abs(H2-t1-t2s)/abs(H2):8.1e} {t2p:10.3e} {t1+t2p:11.3e} {abs(H2-t1-t2p)/abs(H2):8.1e} {t2s/abs(t1):8.4f}")
