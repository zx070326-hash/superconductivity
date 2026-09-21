import numpy as np, os
os.environ["OMP_NUM_THREADS"]="1"
def H_mat(N,w):
    u=1/np.sqrt(2*np.arange(N)+1); k=np.arange(N)
    f=lambda k: w*w/(k*k+w*w)
    L1=np.diag([sum(2*f(j) for j in range(1,n+1))/(2*n+1) for n in range(N)])
    n,m=np.meshgrid(k,k,indexing='ij')
    L2=np.where(n!=m,f(np.abs(n-m)),0.0)*np.outer(u,u)
    L3=f(n+m+1)*np.outer(u,u)
    return -L1+L2+L3
for N in (400,800,1600):
    w,Nc=0.74,4
    A=H_mat(N,w); uc=np.zeros(N); uc[:Nc]=1/np.sqrt(2*np.arange(Nc)+1)
    eA,VA=np.linalg.eigh(A); mu1,mu2=eA[-1],eA[-2]; psi1,psi2=VA[:,-1],VA[:,-2]
    c1,c2=uc@psi1,uc@psi2; gam=mu1-mu2
    print(f"N={N}: mu1={mu1:.6f} mu2={mu2:.6f} gamma={gam:.6f} c1={c1:.5f} c2={c2:.2e}")
    for mt in (0.30,0.322,0.334,0.40):
        M=A-2*mt*np.outer(uc,uc); e=np.linalg.eigvalsh(M); l1,l2=e[-1],e[-2]
        kap=2*mt*c1*c1/gam
        print(f"   mt={mt:.3f} kappa={kap:.3f}  true gap/l1={(l1-l2)/l1:.4f}   (l1-mu2(A))/l1={(l1-mu2)/l1:.4f}   l1={l1:.5f} l2={l2:.5f}")
