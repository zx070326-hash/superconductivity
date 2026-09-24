"""M3 下侧「反向」是不是 N 收敛噪声？闭式 S ⟹ lam1(N) 单调↑到 lam1^inf。
若 lam1(N) 升过 mu_lo，则 mu_lo <= lam1^inf 成立、下侧声音。"""
import os; os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
from scipy.linalg import eigh
R = 10.0
def lam1(w, lam, mus, nd):
    n = np.arange(nd); u = 1.0/np.sqrt(2.0*n+1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    G  = w*w/(w*w + np.arange(2*nd+2)**2)
    Gc = w*w/(w*w + np.arange(nd)**2)
    S  = Gc[0] + 2.0*(np.cumsum(Gc)-Gc[0])
    M = lam*((G[np.abs(i-j)]+G[i+j+1])*np.outer(u,u) - np.diag(u*u*S))
    Nc = int(np.floor(R*w))
    if mus>0 and Nc>0:
        uc = np.where(n<Nc,u,0.0); M = M - 2.0*mus*np.outer(uc,uc)
    return float(eigh(M, subset_by_index=[nd-1,nd-1], eigvals_only=True)[0])
CASES = [(2.0,1.100033,1.021190594),(3.2953,0.731708,1.020081705),
         (4.0,0.640491,1.021883712),(6.0,0.495741,1.025180308)]
for lam,wl,ml in CASES:
    vals=[(nd,lam1(wl,lam,0.20,nd)) for nd in (400,800,1600,2400)]
    s=" ".join(f"N={nd}:{v:.9f}" for nd,v in vals)
    inc=all(vals[k+1][1]>vals[k][1] for k in range(len(vals)-1))
    last=vals[-1][1]
    # Richardson 外推 (增量比 ~ 常数)
    d1,d2=vals[-2][1]-vals[-3][1], vals[-1][1]-vals[-2][1]
    rho=d2/d1 if d1 else 0.0
    ext=last + d2*rho/(1-rho) if 0<rho<1 else float('nan')
    print(f"lam={lam:<7} mu_lo={ml:.9f}")
    print(f"   {s}  单调↑={inc}")
    print(f"   N=2400 余量 {last-ml:+.2e};  外推 lam1^inf≈{ext:.9f} 余量 {ext-ml:+.2e}  "
          f"{'==> mu_lo <= lam1^inf 声音' if ext>ml else '** 仍反向, 需查 **'}")
