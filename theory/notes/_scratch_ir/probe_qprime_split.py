"""q'(x) = <v,O[d(f^2)/dx]v>  +  2*sum_{j>1} <v,O[f^2]w_j><w_j,O[fx]v>/(h-lam_j)
The FIRST term is the envelope part; d(f^2)/dx = 2x k^2/(k^2+x)^3 > 0 for every k.
Is O[.] of a POSITIVE kernel guaranteed to have positive quadratic form on v?
NO in general (the -L1 diagonal has a minus sign). So check term-by-term:
  E := <v, O[g]v>, g := 2x k^2/(k^2+x)^3 > 0
  C := q' - E   (the eigenvector-motion / cross term)
Numbers above show E>0 and E>q'>0, i.e. C<0 but |C| < E.  Quantify |C|/E.
If |C|/E is bounded away from 1 by an analytic argument, concavity follows
WITHOUT any gap certificate -> Lemma C becomes analytic.
Also: is E>0 provable? E = <v,O[g]v> with g>0. Use the SAME kernel-level Jensen /
Lemma M machinery: O[g] for g of the form  2x k^2/(k^2+x)^3  =  -d/dk-ish...
Test: is O[g] a POSITIVE OPERATOR (all eigenvalues >0)? or only positive on v?
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
print(f"{'x':>7} {'E=<v,O[g]v>':>12} {'q\'':>9} {'C=q\'-E':>10} {'|C|/E':>7} {'lam_min(O[g])':>13} {'lam_max':>8}")
for x in [0.05,0.2,0.8,1.6,3.2,6.4,25.6,102.4]:
    N=max(220,int(40*np.sqrt(x))+220); K=2*N+2
    k=np.arange(1,K+1,dtype=float); f=x/(k**2+x)
    A=O_of_kernel(N,f); w,V=np.linalg.eigh(A); i=np.argmax(w); v=V[:,i]; h=w[i]
    g=2.0*x*k**2/(k**2+x)**3
    Og=O_of_kernel(N,g)
    E=v@(Og@v)
    d=1e-4*x
    def q(y):
        Ay=O_of_kernel(N,y/(k**2+y)); wy,Vy=np.linalg.eigh(Ay); iy=np.argmax(wy)
        vy=Vy[:,iy]; return vy@(O_of_kernel(N,(y/(k**2+y))**2)@vy)
    qpr=(q(x+d)-q(x-d))/(2*d)
    eg=np.linalg.eigvalsh(Og)
    print(f"{x:7.2f} {E:12.6f} {qpr:9.6f} {qpr-E:10.6f} {abs(qpr-E)/E:7.4f} {eg.min():13.6f} {eg.max():8.5f}")
