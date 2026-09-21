# -*- coding: utf-8 -*-
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np
exec(open("theory/notes/_scratch_ir/_chk_ml_gap1_attack.py",encoding="utf-8").read().split("N, B = 300")[0])

N,B = 300,2
def build(c11,w1,c22,w2,c12,wx):
    return [[O(kern(c11,w1,N),N), O(kern(c12,wx,N),N)],
            [O(kern(c12,wx,N),N), O(kern(c22,w2,N),N)]]

for tag,(c11,w1,c22,w2,c12,wx) in [("B行1",(2.0,0.3,0.05,20.0,0.001,0.3)),
                                   ("A最劣",(0.2491,7.9005,1.5966,1.7533,0.1387,13.241))]:
    bl=build(c11,w1,c22,w2,c12,wx)
    M=np.zeros((B*N,B*N))
    for i in range(B):
        for j in range(B): M[i*N:(i+1)*N,j*N:(j+1)*N]=bl[i][j]
    e,V=np.linalg.eigh(M); p=V[:,-1].reshape(B,N)
    i0,j0=np.unravel_index(np.argmax(np.abs(p)),p.shape)
    if p[i0,j0]<0: p=-p
    print(f"\n--- {tag}  lam1={e[-1]:.6f}  次大={e[-2]:.6f}  间隙={e[-1]-e[-2]:.2e}")
    for b in range(B):
        d=np.diff(p[b]); k=int(np.argmin(d))
        print(f"  band{b}: ||.||={np.linalg.norm(p[b]):.4f}  头6={np.round(p[b,:6],4)}"
              f"  min={p[b].min():.3e}  最负差分 @n={k}: {d[k]:.3e}  (p[{k}]={p[b,k]:.4f}->p[{k+1}]={p[b,k+1]:.4f})")
    # 单带对照：各自 O 的顶本征向量是否在 C 内（仓库推论 B' 预言 yes）
    for b,(c,w) in enumerate([(c11,w1),(c22,w2)]):
        ev,Vv=np.linalg.eigh(O(kern(c,w,N),N)); v=Vv[:,-1]; v=v if v[np.argmax(abs(v))]>0 else -v
        print(f"  [标量对照 band{b}] lam1={ev[-1]:.6f} min={v.min():.3e} min diff={np.diff(v).min():.3e}")
