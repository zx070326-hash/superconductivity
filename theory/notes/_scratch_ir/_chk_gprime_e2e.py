# -*- coding: utf-8 -*-
"""(G') 端到端：lam1(O_pair) <= lam_max(N) <= g(2)*lam_max(M2)/kappa ?
kappa 归一：取 Lij(k)=c_ij*varpi_ij^2/(k^2+varpi_ij^2)，则 M2_ij/kappa = c_ij*varpi_ij^2
（因 sup_k k^2*Lij(k) = c_ij*varpi_ij^2，正是引理 M 里的 G(k)=<varpi^2>/k^2）。"""
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
def O(fv,N):
    n=np.arange(N); u=1.0/np.sqrt(2*n+1.0)
    T=fv[np.abs(n[:,None]-n[None,:])].copy(); np.fill_diagonal(T,0.0)
    M=(T+fv[n[:,None]+n[None,:]+1])*np.outer(u,u)
    M[n,n]-=2.0*np.concatenate(([0.],np.cumsum(fv[1:N])))/(2*n+1.0); return M
def kern(c,w,N):
    k=np.arange(2*N+2,dtype=float); return c*w*w/(k*k+w*w)
N=260; G2=1.318140212132
rng=np.random.default_rng(23); rows=[]
for B in (2,3,4):
    for t in range(40):
        c=rng.uniform(0.0,2.0,(B,B)); c=(c+c.T)/2
        w=rng.uniform(0.2,25.0,(B,B)); w=(w+w.T)/2
        blocks=[[O(kern(c[i,j],w[i,j],N),N) for j in range(B)] for i in range(B)]
        X=np.block(blocks); lam1=np.linalg.eigvalsh(X)[-1]
        M2=c*w*w                                  # = sup_k k^2 Lij(k)
        Nm=np.zeros((B,B))
        for i in range(B):
            for j in range(B):
                e=np.linalg.eigvalsh(blocks[i][j])
                Nm[i,j]= e[-1] if i==j else max(e[-1],-e[0])   # 对角用 lam1，非对角用范数
        lamN=np.linalg.eigvalsh(Nm)[-1]
        ceil=G2*np.linalg.eigvalsh(M2)[-1]
        rows.append((B,lam1,lamN,ceil, lam1<=lamN+1e-9, lamN<=ceil+1e-9, ceil/lam1))
bad1=[r for r in rows if not r[4]]; bad2=[r for r in rows if not r[5]]
print(f"总 {len(rows)} 组 (B=2,3,4 各 40)")
print(f"  步骤 1  lam1(O_pair) <= lam_max(N) : 违反 {len(bad1)}")
print(f"  步骤 2  lam_max(N)   <= g2*lmax(M2): 违反 {len(bad2)}")
sl=[r[6] for r in rows]
print(f"  端到端松紧 ceil/lam1: 中位 {np.median(sl):.4f}  最小 {min(sl):.4f}  最大 {max(sl):.4f}")
# 饱和检验：所有带共用同一核 ==> 应精确饱和
for B in (2,3):
    cc=rng.uniform(0.3,2.0,(B,B)); cc=(cc+cc.T)/2; ww=7.0
    bl=[[O(kern(cc[i,j],ww,N),N) for j in range(B)] for i in range(B)]
    lam1=np.linalg.eigvalsh(np.block(bl))[-1]
    M2=cc*ww*ww; ceil=G2*np.linalg.eigvalsh(M2)[-1]
    print(f"  [饱和检验 B={B}, 全带同 varpi={ww}] lam1={lam1:.6f} ceil={ceil:.6f} 比={ceil/lam1:.4f}"
          f"   (单带同 varpi 的 lam1(O)/varpi^2={np.linalg.eigvalsh(O(kern(1,ww,N),N))[-1]/ww**2:.6f} vs g2={G2:.6f})")
