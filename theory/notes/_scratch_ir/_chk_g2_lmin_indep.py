"""独立复核缺口 2 的两个数：g(2)=lam_max(O[1/k^2]) 与 lam_min(O[1/k^2])。
用我自己的 O() （不 import 任何生产构造器）。判据只需 |lam_min| < g(2)。"""
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
def O(phi,N):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0)
    pv=phi(np.arange(0,2*N+2,dtype=float))
    T=pv[np.abs(n[:,None]-n[None,:])].copy(); np.fill_diagonal(T,0.0)
    M=(T+pv[n[:,None]+n[None,:]+1])*np.outer(u,u)
    M[n,n]-=2.0*np.concatenate(([0.],np.cumsum(pv[1:N])))/(2*n+1.0)
    return M
phi=lambda k: np.divide(1.0,k*k,out=np.zeros_like(k),where=(k>0))
print(f"{'N':>6} {'g(2)=lam_max':>18} {'lam_min':>18} {'|lam_min|/g(2)':>16}")
for N in (200,400,800,1600,3200):
    e=np.linalg.eigvalsh(O(phi,N)); print(f"{N:>6} {e[-1]:>18.12f} {e[0]:>18.12f} {abs(e[0])/e[-1]:>16.9f}")
print("\n笔记/仓库报: g(2)=1.318140212 (CERT-D rational-exact 包围"
      " [1.318140211762671, 1.318140212373478]); lam_min=-0.989792220; 比 0.750900557")
