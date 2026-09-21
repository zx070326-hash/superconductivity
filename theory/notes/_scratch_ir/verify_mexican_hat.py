import numpy as np, os
os.environ["OMP_NUM_THREADS"]="1"
def O_mat(N,f):
    u=1/np.sqrt(2*np.arange(N)+1); k=np.arange(N)
    L1=np.diag([sum(2*f(j) for j in range(1,n+1))/(2*n+1) for n in range(N)])
    n,m=np.meshgrid(k,k,indexing='ij')
    L2=np.where(n!=m,f(np.abs(n-m)),0.0)*np.outer(u,u); L3=f(n+m+1)*np.outer(u,u)
    return -L1+L2+L3
for N in (500,1000):
    for b in (0.0,0.3,0.5,0.8,0.95):
        f=lambda k,b=b: 1/(1+k**2)-b/(1+(k/4)**2)
        e=np.linalg.eigvalsh(O_mat(N,f)); l1,l2=e[-1],e[-2]
        print(f"N={N} b={b:.2f}: l1={l1:.5f} l2={l2:.5f} gap/l1={(l1-l2)/l1:.4f}")
B=O_mat(500,lambda k:1/(1+(k/4)**2)); eB=np.linalg.eigvalsh(B)[::-1]
print("B=O[phi-] top eigenvalues:",np.round(eB[:6],4))
