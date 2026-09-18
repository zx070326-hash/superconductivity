"""The identity O[fx] = (O[f] - O[f^2])/x  suggests a route with NO gap at all.
Since H(x) = lam_1(O[f(.;x)]) and by envelope/Feynman-Hellmann H'(x) = <v,O[fx]v>
             = (h - <v,O[f^2]v>)/x.
So   x H'(x) = H(x) - q(x),   q(x) := <v_x, O[f^2(.;x)] v_x>.
Differentiate:  H' + x H'' = H' - q'  =>  H'' = -q'(x)/x.
CONCAVITY  <=>  q'(x) > 0  for all x>0.   (!!)
q(x) = <v,O[f^2]v> where f^2(k)=x^2/(k^2+x)^2 -- a kernel that is INCREASING in x
for each k, and O[.] is entrywise monotone in the kernel only if L1 sign allows.
This is a candidate ANALYTIC reformulation: concavity of H  <=>  monotonicity of q.
Test the identity and q' numerically; also test whether q is the top-eigenvalue
quadratic form of a MONOTONE family (which would give q' > 0 by envelope thm:
q'(x) >= <v, d/dx O[f^2] v> if v were the maximizer of q, but v maximizes H not q --
so check the sign of <v, O[d(f^2)/dx] v> separately).
"""
import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
def hq(x,N):
    K=2*N+2; k=np.arange(1,K+1,dtype=float); f=x/(k**2+x)
    A=O_of_kernel(N,f); w,V=np.linalg.eigh(A); i=np.argmax(w)
    v=V[:,i]; h=w[i]
    q=v@(O_of_kernel(N,f**2)@v)
    # d/dx of f^2 = 2x k^2/(k^2+x)^3
    dq_env=v@(O_of_kernel(N,2.0*x*k**2/(k**2+x)**3)@v)
    return h,q,dq_env
print(f"{'x':>8} {'H':>9} {'q':>9} {'xH\'':>9} {'H-q':>9} {'id.err':>9} {'q\'(num)':>9} {'<v,d(f2)v>':>10} {'-q\'/x':>10}")
prev=None
for x in [0.05,0.1,0.2,0.4,0.8,1.6,3.2,6.4,12.8,25.6,51.2]:
    N=max(220,int(40*np.sqrt(x))+220)
    h,q,dqe=hq(x,N)
    d=1e-4*x
    hp,qp,_=hq(x+d,N); hm,qm,_=hq(x-d,N)
    Hp=(hp-hm)/(2*d); qpr=(qp-qm)/(2*d); H2=(hp-2*h+hm)/d**2
    print(f"{x:8.3f} {h:9.5f} {q:9.5f} {x*Hp:9.5f} {h-q:9.5f} {abs(x*Hp-(h-q)):9.1e} {qpr:9.5f} {dqe:10.5f} {-qpr/x:10.5f}   H''num={H2:.5f}")
