# -*- coding: utf-8 -*-
"""(G) 的上确界搜索 + N 收敛方向审计。
方向审计：lam_min(P_N O P_N) >= lam_min(O)（变分）==> 有限 N 的 |lam_min| 是【下估】
==> 比值 |lmin|/lam1 的分子被低估、分母 lam1(P_N) <= lam1(O) 被低估。两边都单边，方向不干净
==> 必须看 N 单调趋势本身，不能只报一个 N。"""
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
from scipy.optimize import minimize
def O(fv,N):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0)
    T=fv[np.abs(n[:,None]-n[None,:])].copy(); np.fill_diagonal(T,0.0)
    M=(T+fv[n[:,None]+n[None,:]+1])*np.outer(u,u)
    M[n,n]-=2.0*np.concatenate(([0.],np.cumsum(fv[1:N])))/(2*n+1.0); return M
def rat(fv,N):
    e=np.linalg.eigvalsh(O(fv,N)); return abs(e[0])/e[-1], e[-1], e[0]

print("=== A. N 收敛（两端都是下估 ==> 看趋势）===")
def mk_step(K,N):
    f=np.zeros(2*N+2); f[1:K+1]=1.0; return f
print(f"{'核':>16} " + " ".join(f"{'N='+str(n):>22}" for n in (150,300,600,1200)))
for nm,mk in [("阶梯 K=1",lambda N:mk_step(1,N)),("阶梯 K=2",lambda N:mk_step(2,N)),
              ("1/k^2",lambda N:np.concatenate(([0.],1.0/np.arange(1,2*N+2)**2.0))),
              ("k^-3",lambda N:np.concatenate(([0.],1.0/np.arange(1,2*N+2)**3.0)))]:
    out=[]
    for N in (150,300,600,1200):
        r,l1,lm=rat(mk(N),N); out.append(f"{r:.6f}({l1:.3f},{lm:.3f})")
    print(f"{nm:>16} " + " ".join(f"{o:>22}" for o in out))

print("\n=== B. 上确界搜索：f = 非负非增，参数化为 M 个台阶高度的累加 ===")
N=300; Ks=np.array([1,2,3,4,6,9,14,21,32,50,80,130,210,299])
def f_from(z,N):
    h=np.exp(z); h=h/h.sum()                 # 各台阶权重 >0 归一
    fv=np.zeros(2*N+2)
    for w,K in zip(h,Ks): fv[1:K+1]+=w       # 非负非增（阶梯正组合）
    return fv
def neg(z):
    r,_,_=rat(f_from(z,N),N); return -r
best=(-1,None)
rng=np.random.default_rng(5)
for trial in range(14):
    z0=rng.normal(0,2.0,len(Ks)) if trial else np.array([6.]+[-6.]*(len(Ks)-1))
    res=minimize(neg,z0,method="Nelder-Mead",options=dict(maxiter=4000,xatol=1e-6,fatol=1e-10))
    if -res.fun>best[0]: best=(-res.fun,res.x)
h=np.exp(best[1]); h/=h.sum()
print(f"  sup 搜索最大比值 = {best[0]:.6f}")
print(f"  最优台阶权重（K:权重，只列 >1e-3）: "
      + ", ".join(f"{K}:{w:.4f}" for K,w in zip(Ks,h) if w>1e-3))
r,l1,lm=rat(f_from(best[1],N),N); print(f"  该点 lam1={l1:.6f} lmin={lm:.6f}")
for NN in (600,1200):
    rr,ll1,llm=rat(f_from(best[1],NN),NN); print(f"  同一 f 在 N={NN}: 比值={rr:.6f} lam1={ll1:.6f} lmin={llm:.6f}")
