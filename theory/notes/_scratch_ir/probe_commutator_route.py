"""Route 3: bound ||(I-P1)O[f^2]v|| via the commutator.
If [A, B] = 0 (A=O[f], B=O[f^2]) then B preserves eigenspaces => (I-P1)Bv = 0.
General bound: for v the top eigenvector of A with gap g,
   ||(I-P1)Bv|| = ||(I-P1)(Bv - <v,Bv>v)|| and
   (A-h)(I-P1)Bv = (I-P1)[A,B]v   =>   ||(I-P1)Bv|| <= ||[A,B]v|| / g.
That is EXACT (no trial vector, no Davis-Kahan) and needs only ||[A,B]v|| <= ||[A,B]||.
Test tightness, then feed into |C| <= (2/x)||(I-P1)O[f2]v||^2/g and check H'' < 0.
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
THETA=0.6349
print(f"{'varpi':>6} {'g_true':>7} {'rperp':>8} {'|[A,B]v|/g':>10} {'ratio':>6} {'E':>9} {'|C|bnd':>9} {'H\"bnd':>10} {'ok':>4}")
allok=True
for vp in [0.2,0.5,1.0,1.5,2.0,3.0,5.0,10.0,20.0]:
    x=vp*vp; N=max(240,int(40*vp)+240); K=2*N+2
    k=np.arange(1,K+1,dtype=float); f=x/(k**2+x)
    A=O_of_kernel(N,f); B=O_of_kernel(N,f**2)
    w,V=np.linalg.eigh(A); o=np.argsort(w)[::-1]; w=w[o]; V=V[:,o]
    h=w[0]; v=V[:,0]; g=h-w[1]
    Bv=B@v; rperp=np.linalg.norm(Bv-(v@Bv)*v)
    commv=np.linalg.norm(A@Bv-B@(A@v))
    bnd_r=commv/g
    E=v@(O_of_kernel(N,2.0*x*k**2/(k**2+x)**3)@v)
    Cbnd=2.0*bnd_r**2/(x*g)
    H2=-(E-Cbnd)/x
    ok=(E-Cbnd)>0; allok&=ok
    print(f"{vp:6.2f} {g:7.4f} {rperp:8.5f} {bnd_r:10.5f} {bnd_r/max(rperp,1e-15):6.2f} {E:9.6f} {Cbnd:9.6f} {H2:10.5f} {str(ok):>4}")
print("\nCommutator route closes with TRUE gap:", allok)
