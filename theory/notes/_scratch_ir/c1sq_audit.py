import numpy as np
# 主会话独立重算 c1sq_lower_bound 的表。从零构造，float64（区间层只影响末位）。
def H(N,w):
    n=np.arange(N+1); f=lambda k: w*w/(w*w+np.asarray(k,float)**2)
    u=1/np.sqrt(2*n+1.0)
    L1=np.zeros(N+1); L1[1:]=2*np.cumsum(f(np.arange(1,N+1)))/(2*n[1:]+1.0)
    T=f(np.abs(n[:,None]-n[None,:])); np.fill_diagonal(T,0.0)
    Hk=f(n[:,None]+n[None,:]+1)
    return -np.diag(L1)+(T+Hk)*u[:,None]*u[None,:]

def bound(N,w):
    f0=np.pi*w/np.tanh(np.pi*w); fpi=np.pi*w/np.sinh(np.pi*w)
    uN=1/np.sqrt(2*(N+1)+1.0)
    nBB=(f0-1)/(2*N+3)+uN**2*(f0+1); nC=uN*(f0-fpi)/2; E=nBB+nC
    ev,V=np.linalg.eigh(H(N,w)); l1,l2=ev[-1],ev[-2]
    d=l1-max(l2,0.0)
    if E>=d: return E,d,None,0.0
    s=E/(d-E)
    if s>=1: return E,d,s,0.0
    c=np.sqrt(1-s*s)*abs(V[0,-1])-s
    return E,d,s,max(0.0,c)**2

print(f"{'w':>5s} {'||E||':>8s} {'delta':>8s} {'sinTh':>8s} {'c1^2 LB':>9s} {'c1^2 true(N=4000)':>18s}")
for w in (0.5,1,2,3,5,10,30):
    E,d,s,lb=bound(500,w)
    ev,V=np.linalg.eigh(H(4000,w)); true=V[0,-1]**2
    print(f"{w:5.1f} {E:8.4f} {d:8.4f} {0 if s is None else s:8.4f} {lb:9.4f} {true:18.4f}")
E,d,s,lb=bound(1000,30); ev,V=np.linalg.eigh(H(4000,30))
print(f"\nw=30 N=1000: ||E||={E:.4f} delta={d:.4f} sinTh={s:.4f} c1^2 LB={lb:.4f}  true={V[0,-1]**2:.4f}")

# 附：若 ||C|| 用真值（理想 Schur 界的代理），varpi=30 在不加大 N 的情况下能到哪
def with_trueC(N,w,tailmul=8):
    f=lambda k: w*w/(w*w+np.asarray(k,float)**2); M2=int(N*tailmul)
    n=np.arange(N+1); m=np.arange(N+1,M2+1)
    un=1/np.sqrt(2*n+1.0); um=1/np.sqrt(2*m+1.0)
    nC=np.linalg.norm((f(np.abs(n[:,None]-m[None,:]))+f(n[:,None]+m[None,:]+1))*un[:,None]*um[None,:],2)
    f0=np.pi*w/np.tanh(np.pi*w); uN=1/np.sqrt(2*(N+1)+1.0)
    nBB=(f0-1)/(2*N+3)+uN**2*(f0+1); E=nBB+nC
    ev,V=np.linalg.eigh(H(N,w)); d=ev[-1]-max(ev[-2],0.0); s=E/(d-E)
    return nC,nBB,s,max(0.0,np.sqrt(1-s*s)*abs(V[0,-1])-s)**2
print()
for w in (10.0,30.0):
    nC,nBB,s,lb=with_trueC(500,w)
    print(f"w={w:4.0f} N=500: ||C||true={nC:.4f} ||H_BB||={nBB:.4f} sinTh={s:.4f} -> c1^2 LB={lb:.4f}")
