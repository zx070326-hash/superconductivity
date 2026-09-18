"""Same route but with ANALYTIC gap lower bounds only (no true lam_2), and
||[A,B]v|| replaced by ||[A,B]|| (Schur/Holmgren-boundable, no eigenvector).
   |C| <= (2/x) * ( ||[A,B]|| / g_lo )^2 / g_lo   [g_lo used twice: once in the
   residual bound, once in the C denominator]
   need  E_lo - |C|bnd > 0.   E also needs an analytic LOWER bound; test with the
   algebraic reduction E = 2(q - <v,O[f^3]v>)/x first (verify it), then bound.
Analytic gap: route R gives h - Theta* for varpi>=8; rearrangement gap note gives
lam_2/h <= 0.58 on [0.2,2].  Use g_lo = 0.42h on [0.2,2], h-0.6349 for varpi>=2.5.
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
def schur(M):
    r=np.abs(M).sum(1); return r.max()
THETA=0.6349
print(f"{'varpi':>6} {'E':>9} {'E_alg':>9} {'errE':>8} {'g_lo':>7} {'|[A,B]|':>8} {'|C|bnd':>10} {'E-|C|':>10} {'ok':>5}")
allok=True
for vp in [0.2,0.5,1.0,1.5,2.0,2.5,3.0,5.0,10.0,20.0]:
    x=vp*vp; N=max(240,int(40*vp)+240); K=2*N+2
    k=np.arange(1,K+1,dtype=float); f=x/(k**2+x)
    A=O_of_kernel(N,f); B=O_of_kernel(N,f**2); C3=O_of_kernel(N,f**3)
    w,V=np.linalg.eigh(A); o=np.argsort(w)[::-1]; w=w[o]; V=V[:,o]
    h=w[0]; v=V[:,0]
    E=v@(O_of_kernel(N,2.0*x*k**2/(k**2+x)**3)@v)
    E_alg=2.0*((v@(B@v))-(v@(C3@v)))/x
    g_lo = 0.42*h if vp<2.5 else max(h-THETA,1e-9)
    nc=schur(A@B-B@A)
    Cbnd=2.0*(nc/g_lo)**2/(x*g_lo)
    ok=(E-Cbnd)>0; allok&=ok
    print(f"{vp:6.2f} {E:9.6f} {E_alg:9.6f} {abs(E-E_alg):8.1e} {g_lo:7.4f} {nc:8.5f} {Cbnd:10.6f} {E-Cbnd:10.6f} {str(ok):>5}")
print("\nAll close with analytic gap + operator-norm commutator:", allok)
