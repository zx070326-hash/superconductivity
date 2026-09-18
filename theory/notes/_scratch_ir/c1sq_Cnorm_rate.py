import numpy as np
# 头尾耦合块 C 的真实范数衰减率（决定"换更快的耦合块估计"值不值得投算力）
def Cnorm(N,w,tailmul=8):
    M2=int(N*tailmul); f=lambda k: w*w/(w*w+np.asarray(k,float)**2)
    n=np.arange(N+1); m=np.arange(N+1,M2+1)
    un=1/np.sqrt(2*n+1.0); um=1/np.sqrt(2*m+1.0)
    C=(f(np.abs(n[:,None]-m[None,:]))+f(n[:,None]+m[None,:]+1))*un[:,None]*um[None,:]
    return np.linalg.norm(C,2)
for w in (1.0,30.0):
    print(f"\nvarpi={w}")
    prev=None
    for N in (125,250,500,1000,2000):
        uN=1/np.sqrt(2*(N+1)+1.0); f0=np.pi*w/np.tanh(np.pi*w); fpi=np.pi*w/np.sinh(np.pi*w)
        cur=Cnorm(N,w); analytic=uN*(f0-fpi)/2
        rate="" if prev is None else f"  rate={np.log2(prev/cur):+.2f}"
        print(f"  N={N:5d}  ||C||_true={cur:.5e}  解析上界={analytic:.5e}  松弛={analytic/cur:7.1f}x{rate}")
        prev=cur
