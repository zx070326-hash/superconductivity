# -*- coding: utf-8 -*-
"""(1) 修正度量后重做缺口 1 定向攻击。 (2) 测猜想 (G): lam1(O[f]) >= |lam_min(O[f])| 对 f 非负非增。"""
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
def O(fv,N):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0)
    T=fv[np.abs(n[:,None]-n[None,:])].copy(); np.fill_diagonal(T,0.0)
    M=(T+fv[n[:,None]+n[None,:]+1])*np.outer(u,u)
    M[n,n]-=2.0*np.concatenate(([0.],np.cumsum(fv[1:N])))/(2*n+1.0); return M
def kern(c,w,N):
    k=np.arange(2*N+2,dtype=float); return c*w*w/(k*k+w*w)
def coneviol(v):
    """v 在 C={>=0, 非增} 内 <=> 返回值 <= 0（相对 max|v| 归一）。"""
    s=np.abs(v).max();  v = v if v[np.argmax(np.abs(v))]>0 else -v
    return max(-v.min(), np.diff(v).max())/s

N=300
print("=== 1. 缺口 1（度量已修正：违反 = max(-min v, max diff)/max|v|）===")
rng=np.random.default_rng(11); worst=(-1e9,None)
for t in range(1500):
    c11,c22=rng.uniform(0.05,3.0,2); c12=rng.uniform(0.0,3.0)
    w1,w2,wx=rng.uniform(0.2,25.0,3)
    M=np.zeros((2*N,2*N))
    M[:N,:N]=O(kern(c11,w1,N),N); M[N:,N:]=O(kern(c22,w2,N),N)
    X=O(kern(c12,wx,N),N); M[:N,N:]=X; M[N:,:N]=X
    p=np.linalg.eigh(M)[1][:,-1].reshape(2,N)
    i0,j0=np.unravel_index(np.argmax(np.abs(p)),p.shape)
    if p[i0,j0]<0: p=-p
    s=np.abs(p).max()
    vv=max(max(-p[b].min(), np.diff(p[b]).max()) for b in range(2))/s
    if vv>worst[0]: worst=(vv,(c11,c22,c12,w1,w2,wx))
print(f"  1500 组随机（c12 不受 sqrt(c11*c22) 约束）最劣违反 = {worst[0]:.3e}  参数={np.round(worst[1],3)}")
print("  ==> <=0 即 0 违反；agent 的 0/400 我用 1500 组更宽的采样复核")

print("\n=== 2. 猜想 (G): lam1 >= |lam_min| ? 报比值 |lam_min|/lam1（<1 即成立）===")
def ratio(fv,N):
    e=np.linalg.eigvalsh(O(fv,N)); return abs(e[0])/e[-1], e[-1], e[0]
kk=np.arange(2*N+2,dtype=float)
cases=[]
cases.append(("1/k^2 (缺口 2)", np.divide(1.0,kk*kk,out=np.zeros_like(kk),where=(kk>0))))
for w in (0.1,0.3,1.0,3.0,10.0,30.0,100.0,1000.0): cases.append((f"Einstein w={w}", kern(1.0,w,N)))
for K in (1,2,3,5,10,30,100,299): 
    f=np.zeros_like(kk); f[1:K+1]=1.0; cases.append((f"阶梯 1_[k<={K}]", f))
for a in (0.5,1.0,1.5,3.0):
    f=np.zeros_like(kk); f[1:]=kk[1:]**(-a); cases.append((f"幂 k^-{a}", f))
f=np.zeros_like(kk); f[1:]=np.exp(-0.3*kk[1:]); cases.append(("指数 e^-0.3k", f))
print(f"{'核':>22} {'|lmin|/lam1':>13} {'lam1':>12} {'lmin':>12}")
bad=0
for nm,fv in cases:
    r,l1,lm=ratio(fv,N); flag="  <== 违反!" if r>=1 else ""
    if r>=1: bad+=1
    print(f"{nm:>22} {r:>13.6f} {l1:>12.6f} {lm:>12.6f}{flag}")
rng=np.random.default_rng(3); wr=(0.0,None)
for t in range(600):   # 随机非负非增核
    inc=rng.exponential(1.0,60)*rng.binomial(1,0.4,60)
    f=np.zeros_like(kk); tail=np.cumsum(inc[::-1])[::-1]; f[1:61]=tail; f[61:]=0.0
    if f[1]<=0: continue
    r,_,_=ratio(f,N); bad+= (r>=1)
    if r>wr[0]: wr=(r,None)
print(f"\n  600 个随机非负非增核：最大比值 = {wr[0]:.6f}   总违反数 = {bad}")
