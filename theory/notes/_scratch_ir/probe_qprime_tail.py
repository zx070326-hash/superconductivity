import numpy as np
def O_of_kernel(N,fv):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); i=n[:,None]; j=n[None,:]
    A=np.where(i!=j,fv[np.clip(np.abs(i-j)-1,0,None)],0.0)+fv[i+j]
    A=A*u[:,None]*u[None,:]
    c=np.concatenate(([0.0],np.cumsum(fv))); A[n,n]-=2.0*u**2*c[n]
    return A
print("large-x behaviour of |C|/E  (does it approach 1?)")
print(f"{'x':>9} {'varpi':>7} {'E':>10} {'C':>11} {'|C|/E':>7} {'N':>6}")
for x in [102.4,204.8,409.6,819.2,1638.4]:
    vp=np.sqrt(x); N=int(35*vp)+260; K=2*N+2
    k=np.arange(1,K+1,dtype=float); f=x/(k**2+x)
    A=O_of_kernel(N,f); w,V=np.linalg.eigh(A); i=np.argmax(w); v=V[:,i]
    E=v@(O_of_kernel(N,2.0*x*k**2/(k**2+x)**3)@v)
    d=1e-4*x
    def q(y):
        Ay=O_of_kernel(N,y/(k**2+y)); wy,Vy=np.linalg.eigh(Ay); iy=np.argmax(wy)
        vy=Vy[:,iy]; return vy@(O_of_kernel(N,(y/(k**2+y))**2)@vy)
    qpr=(q(x+d)-q(x-d))/(2*d); C=qpr-E
    print(f"{x:9.1f} {vp:7.2f} {E:10.6f} {C:11.6f} {abs(C)/E:7.4f} {N:6d}")
