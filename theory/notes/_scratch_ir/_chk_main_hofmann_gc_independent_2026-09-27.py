# Independent check (main session): strong-coupling limit g_k = g/k^2 (k>=1), mu*=0.
# Principal-submatrix kernel; find g_c with lambda1=1; C = 1/(2*pi*sqrt(g_c)); compare C_inf.
import numpy as np, os
os.environ.setdefault("OMP_NUM_THREADS","1")
def lam1(g,N):
    n=np.arange(N); u=1/np.sqrt(2*n+1)
    k=np.abs(n[:,None]-n[None,:]); s=n[:,None]+n[None,:]+1
    gk=np.where(k>0,g/np.maximum(k,1)**2,0.0)
    K=np.outer(u,u)*(gk+g/s**2)
    H=np.concatenate([[0],np.cumsum(g/np.arange(1,N)**2)])  # sum_{k=1}^{i} g_k
    K-=np.diag(u**2*2*H)
    return np.linalg.eigvalsh(K)[-1]
def gc(N):
    a,b=0.1,5.0
    for _ in range(60):
        m=(a+b)/2
        (a,b)=(a,m) if lam1(m,N)>1 else (m,b)
    return (a+b)/2
Cinf=0.1827262477
for N in [200,400,800,1600]:
    g=gc(N); print(N, g, 1/(2*np.pi*np.sqrt(g)), 1/(2*np.pi*np.sqrt(g))/Cinf)
