"""Find an ANALYTIC upper bound for r_perp := ||(I-P1) O[fx] v||.

Key structural fact to exploit: fx(k) = k^2/(k^2+x)^2 = (1/x)*(f(k) - f(k)^2)
where f(k)=x/(k^2+x).  Since O is entrywise linear in the kernel,
   O[fx] = (1/x)( O[f] - O[f^2] ) = (1/x)( A - O[f^2] ).
Then O[fx] v = (1/x)( h v - O[f^2] v ), and the PERP part kills h v:
   (I-P1) O[fx] v = -(1/x) (I-P1) O[f^2] v.
So  r_perp = (1/x) * ||(I-P1) O[f^2] v|| <= (1/x) * ||O[f^2]||.
That is a CLOSED-FORM route: f^2(k) = x^2/(k^2+x)^2 in l^1 with explicit sum.
Test: (i) the identity O[fx]=(A-O[f^2])/x exactly; (ii) how tight (1/x)||O[f^2]||
is vs r_perp; (iii) whether term1 + 2((1/x)||O[f2]||)^2/g_anal < 0.
"""
import numpy as np, sys
def O_of_kernel(N, fvals):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j, fvals[np.clip(np.abs(i-j)-1,0,None)],0.0)+fvals[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fvals))); A[n,n]-=2.0*u**2*c[n]
    return A
def schur(A):
    r=np.abs(A).sum(1); return r.max()
THETA=0.6349
print(f"{'x':>7} {'id.err':>9} {'r_perp':>9} {'(1/x)|O[f2]|':>12} {'ratio':>7} {'term1':>10} {'H\"bnd':>10} {'closes':>7}")
ok=True
for x in [0.04,0.1,0.25,1.0,4.0,9.0,25.0,100.0,400.0]:
    N=max(200,int(40*np.sqrt(x))+200); K=2*N+2
    k=np.arange(1,K+1,dtype=float)
    f=x/(k**2+x); A=O_of_kernel(N,f)
    Ax=O_of_kernel(N,k**2/(k**2+x)**2)
    Af2=O_of_kernel(N,f**2)
    Axx=O_of_kernel(N,-2.0*k**2/(k**2+x)**3)
    iderr=np.abs(Ax-(A-Af2)/x).max()
    w,V=np.linalg.eigh(A); o=np.argsort(w)[::-1]; w=w[o]; V=V[:,o]
    h=w[0]; v=V[:,0]; g=h-w[1]
    r=Ax@v; rp=np.linalg.norm(r-(v@r)*v)
    bnd=schur(Af2)/x
    t1=v@(Axx@v)
    g_an = h*(1-0.58) if x<=4 else max(h-THETA,1e-9)
    H2=t1+2.0*bnd**2/g_an
    ok &= H2<0
    print(f"{x:7.2f} {iderr:9.1e} {rp:9.5f} {bnd:12.5f} {bnd/rp:7.2f} {t1:10.5f} {H2:10.5f} {str(H2<0):>7}")
print("\nAll close with closed-form (1/x)||O[f^2]|| + analytic gap:", ok)
