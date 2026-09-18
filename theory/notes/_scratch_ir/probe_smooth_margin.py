"""Two questions peer left open:
(1) Is the p=2, r=4, mu=0.26 residual dip (-4e-7) numerical (N-convergence) or real?
(2) Where does the analytic margin mu*D_c/D_ph < 1 actually FAIL? Peer reports 0.985 at
    varpi=0.5 -- that is NOT safe. Scan small w to find the true breakdown point.
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
def build(w,r,mu,p,N):
    K=2*N+2; k=np.arange(1,K+1,dtype=float); x=w*w
    A=O_of_kernel(N,x/(k**2+x))
    n=np.arange(N); t=(2*n+1.0)/(2*r*w); chi=1.0/(1.0+t**(2*p))
    u=chi/np.sqrt(2*n+1.0)
    return A-mu*2*np.outer(u,u),A,u,t,n,k,x
print("(1) N-convergence of the worst residual dip, r=4 mu=0.26 p=2")
ws=np.arange(0.5,8.001,0.05)
for N in [400,700,1100,1600]:
    v=np.array([np.linalg.eigvalsh(build(w,4,0.26,2,N)[0]).max() for w in ws])
    d=np.diff(v); b=np.where(d<0)[0]
    print(f"   N={N:5d}  #down={len(b):2d}  worst={(d[b].min() if len(b) else 0.0):11.3e}  at w={[round(ws[i],2) for i in b]}")
print()
print("(2) margin ratio mu*D_c/D_ph, r=3 mu=0.26 p=2  (peer: 0.985 at w=0.5)")
print(f"   {'w':>5} {'Dph':>11} {'Dc':>11} {'ratio':>9} {'dlam1':>11}")
for w in [0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.5,2.0,3.0]:
    N=int(min(1400,max(320,5*max(w,3*w))))
    M,A,u,t,n,k,x=build(w,3,0.26,2,N)
    ev,V=np.linalg.eigh(M); v=V[:,-1]
    Ap=O_of_kernel(N,k**2/(k**2+x)**2)*(2*w)   # dO/dw = dO/dx * 2w
    Dph=v@(Ap@v)
    dchi=-4.0*t**3/(1.0+t**4)**2               # p=2
    up=-(t/w)*dchi/np.sqrt(2*n+1.0)
    Dc=4*(v@u)*(v@up)
    print(f"   {w:5.2f} {Dph:11.3e} {Dc:11.3e} {0.26*Dc/Dph:9.4f} {Dph-0.26*Dc:11.3e}")
