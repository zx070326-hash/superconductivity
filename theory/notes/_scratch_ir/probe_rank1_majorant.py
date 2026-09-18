import numpy as np, os
os.environ["OMP_NUM_THREADS"]="1"
def build(N,f):
    u=1/np.sqrt(2*np.arange(N)+1); k=np.arange(N); n,m=np.meshgrid(k,k,indexing='ij')
    L1=np.diag(np.concatenate(([0.0],np.cumsum([2*f(j) for j in range(1,N)])))/(2*k+1))
    L2=np.where(n!=m,f(np.abs(n-m)),0.0)*np.outer(u,u); L3=f(n+m+1)*np.outer(u,u)
    return L1,L2,L3,u
lor=lambda w:(lambda k:w*w/(k*k+w*w))

print("=== LEMMA CHECK:  O'[f] = L1+L2+L3  >= -f(0) e0 e0^T   when 2f(1)>=f(0) ===")
print("   (report lambda_min of O'[f] + f(0) e0e0^T ; must be >= 0)")
for N in (200,400):
  for w in (0.3,0.5,0.8,1.0,1.5,3.0,10.0,100.0):
    f=lor(w); L1,L2,L3,u=build(N,f); Op=L1+L2+L3
    E=np.zeros((N,N)); E[0,0]=f(0)
    lm_raw=np.linalg.eigvalsh(Op)[0]; lm_fix=np.linalg.eigvalsh(Op+E)[0]
    cond=2*f(1)-f(0)
    print(f"  N={N} w={w:6.1f}  2f(1)-f(0)={cond:+.4f}  lmin(O')={lm_raw:+.5f}   lmin(O'+f(0)e0e0^T)={lm_fix:+.3e}")
