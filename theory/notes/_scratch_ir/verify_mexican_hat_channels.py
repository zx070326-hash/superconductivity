import numpy as np, os
os.environ["OMP_NUM_THREADS"]="1"
N=500
u=1/np.sqrt(2*np.arange(N)+1); k=np.arange(N); n,m=np.meshgrid(k,k,indexing='ij')
def parts(f):
    L1=np.diag([sum(2*f(j) for j in range(1,q+1))/(2*q+1) for q in range(N)])
    L2=np.where(n!=m,f(np.abs(n-m)),0.0)*np.outer(u,u); L3=f(n+m+1)*np.outer(u,u)
    return L1,L2,L3
fp=lambda k:1/(1+k**2); fm=lambda k:1/(1+(k/4)**2)
L1p,L2p,L3p=parts(fp); L1m,L2m,L3m=parts(fm)
A=-L1p+L2p+L3p
print("b    | full O[phi+]-b O[phi-]  | Delta-channel only (-b(L2+L3)) | Z-only (+b L1)")
for b in (0.0,0.3,0.5,0.8,0.95,1.2,1.5):
    def g(M): e=np.linalg.eigvalsh(M); return (e[-1]-e[-2])/e[-1] if e[-1]>0 else float('nan'), e[-1]
    gf,lf=g(A-b*(-L1m+L2m+L3m)); gd,ld=g(A-b*(L2m+L3m)); gz,lz=g(A+b*L1m)
    print(f"{b:4.2f} | gap/l1={gf:.4f} l1={lf:.4f}    | gap/l1={gd:.4f} l1={ld:.4f}         | gap/l1={gz:.4f} l1={lz:.4f}")
print("\nPhysical spin-fluctuation sign structure: Z gets +b*phi-, Delta gets -b*phi-  (Berk-Schrieffer)")
for b in (0.3,0.5,0.8,0.95,1.2):
    M=A-b*L1m-b*(L2m+L3m); e=np.linalg.eigvalsh(M)
    print(f"b={b:.2f}: l1={e[-1]:.4f} l2={e[-2]:.4f} gap/l1={(e[-1]-e[-2])/e[-1] if e[-1]>0 else float('nan'):.4f}")
print("\nOscillating kernel phi=cos(k/3)/(1+k^2)... full O vs Delta-only")
fo=lambda k: np.cos(k/3.0)/(1+(k/3.0)**2)
L1o,L2o,L3o=parts(fo)
for nm,M in (("full O[phi_osc]",-L1o+L2o+L3o),("Delta-only w/ Z from |phi|",-parts(lambda k:abs(fo(k)))[0]+L2o+L3o)):
    e=np.linalg.eigvalsh(M); print(f"{nm}: l1={e[-1]:.4f} l2={e[-2]:.4f} gap/l1={(e[-1]-e[-2])/e[-1]:.4f}")
