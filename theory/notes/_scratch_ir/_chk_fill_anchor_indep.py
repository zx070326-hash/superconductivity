"""独立复核 fill_factor_probe 的库仑亏损三点表 (mu*=0.13/0.10, r=10)：
15.5%/9.4%/0.0% @ lambda=10/30/100。
第二条路径：稠密 eigvalsh 直接对 M = A - 2 mu~ u_c u_c^T（不走久期根、不 import build）。
亏损 := 1 - Tc(mu*)/Tc(0) = 1 - varpi_0/varpi_mu  （Tc = w_E/(2 pi varpi)）"""
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
def Omat(w,N):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0); kk=np.arange(2*N+2,dtype=float)
    fv=w*w/(kk*kk+w*w)
    L1=np.zeros(N); L1[1:]=2*np.cumsum(fv[1:N])/(2*n[1:]+1.0)
    d=np.abs(n[:,None]-n[None,:]); F2=fv[d].copy(); np.fill_diagonal(F2,0.0)
    A=(F2+fv[n[:,None]+n[None,:]+1])*np.outer(u,u); A[np.diag_indices(N)]-=L1
    return A,u
def lam1(w,mut,r,N):
    A,u=Omat(w,N); Nc=int(np.floor(r*w))
    if Nc<=0 or mut==0.0: return np.linalg.eigvalsh(A)[-1], Nc
    assert Nc<N
    y=np.zeros(N); y[:Nc]=u[:Nc]
    return np.linalg.eigvalsh(A-2*mut*np.outer(y,y))[-1], Nc
def varpi_c(lam,mus,r,N,lo=1e-4,hi=None):
    hi=min(40.0,(N-2)/r) if hi is None else hi
    f=lambda w: lam*lam1(w,mus/lam if mus else 0.0,r,N)[0]-1.0
    assert f(hi)>0 and f(lo)<0, (f(lo),f(hi))
    for _ in range(60):
        m=0.5*(lo+hi)
        if f(m)<0: lo=m
        else: hi=m
    return 0.5*(lo+hi)
r=10.0
print(f"{'lam':>5} {'N':>5} | {'varpi(0)':>10} {'Nc(0)':>5} | {'varpi(.13)':>11} {'Nc':>4} {'亏损':>7}"
      f" | {'varpi(.10)':>11} {'Nc':>4} {'亏损':>7}")
for lam,N in ((10.,400),(30.,400),(100.,400),(10.,800),(100.,800)):
    w0=varpi_c(lam,0.0,r,N); _,nc0=lam1(w0,0.0,r,N)
    out=[]
    for mus in (0.13,0.10):
        w=varpi_c(lam,mus,r,N); _,nc=lam1(w,mus/lam,r,N)
        out.append((w,nc,1.0-w0/w))
    print(f"{lam:>5.0f} {N:>5} | {w0:>10.5f} {nc0:>5} | {out[0][0]:>11.5f} {out[0][1]:>4} "
          f"{out[0][2]*100:>6.1f}% | {out[1][0]:>11.5f} {out[1][1]:>4} {out[1][2]*100:>6.1f}%")
print("\n笔记报: lam=10 -> 15.5%/12.5% ; lam=30 -> 9.4%/7.5% ; lam=100 -> 0.0%/0.0%")
