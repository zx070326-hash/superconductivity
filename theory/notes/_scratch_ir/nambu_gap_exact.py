import numpy as np
rng = np.random.default_rng(7)
M = 200000; mm = np.arange(M+1); ii = np.arange(M+1,dtype=float)
def spec(varpi,g,N=600):
    lam=lambda k: varpi**2/(np.asarray(k,dtype=float)**2+varpi**2)
    L=lam(np.arange(0,M+N+2)); n=np.arange(N)
    B=np.array([np.dot(L[np.abs(k-mm)],g)-np.dot(L[k+1+mm],g)-L[0]*g[k] for k in n])
    w=g[:N]/(2*n+1.0); u=np.sqrt(w)
    T=L[np.abs(n[:,None]-n[None,:])]; np.fill_diagonal(T,0.0); H=L[n[:,None]+n[None,:]+1]
    ev=np.linalg.eigvalsh(-np.diag(B/(2*n+1.0))+(T+H)*u[:,None]*u[None,:])   # M = -L1 + P
    return ev[-1], ev[-2]
shapes={'flat':np.ones(M+1),'vH_w3_d.02':np.where(ii<3,1.0,0.02),
        'vH_smooth':0.05+0.95/(1+ii/8),'pg_w5_d.05':np.where(ii<5,0.05,1.0),
        'wild':np.exp(rng.uniform(-2,2,M+1))}
print(f"{'varpi':>6s} {'shape':12s} {'lam1':>10s} {'lam2':>10s} {'gap/lam1':>9s}")
for varpi in (0.3,1.0,5.0):
    for nm,g in shapes.items():
        l1,l2=spec(varpi,g); print(f"{varpi:6.1f} {nm:12s} {l1:10.5f} {l2:10.5f} {(l1-l2)/l1:9.4f}")
