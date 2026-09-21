import numpy as np
# 直接测 Lemma 4 的算子结论（而非 a_n 这条充分路线）：O' = +L1 + P 是否 >= 秩一泄漏
rng = np.random.default_rng(7)
M = 200000; mm = np.arange(M+1); ii = np.arange(M+1,dtype=float)

def run(varpi, g, N=600):
    lam = lambda k: varpi**2/(np.asarray(k,dtype=float)**2+varpi**2)
    rho = lam(0) - np.pi*varpi/np.sinh(np.pi*varpi)
    L = lam(np.arange(0, M+N+2))
    n = np.arange(N)
    B = np.array([np.dot(L[np.abs(k-mm)],g)-np.dot(L[k+1+mm],g)-L[0]*g[k] for k in n])
    L1 = B/(2*n+1.0)
    w  = g[:N]/(2*n+1.0); u = np.sqrt(w)
    T  = L[np.abs(n[:,None]-n[None,:])]; np.fill_diagonal(T,0.0)
    H  = L[n[:,None]+n[None,:]+1]
    P  = (T+H)*u[:,None]*u[None,:]
    Op = np.diag(L1) + P                       # O' = +L1 + P
    a0 = L1[0] - rho*w[0]
    ev = np.linalg.eigvalsh(Op)
    Op2 = Op.copy(); Op2[0,0] -= a0            # 扣掉 n=0 那一格的秩一泄漏
    return a0, ev[0], np.linalg.eigvalsh(Op2)[0]

shapes = {
 'flat'      : np.ones(M+1),
 'vH_w3_d.02': np.where(ii<3, 1.0, 0.02),
 'vH_w10_d.1': np.where(ii<10,1.0, 0.10),
 'vH_smooth' : 0.05+0.95/(1+ii/8),
 'pg_w5_d.05': np.where(ii<5, 0.05, 1.0),
 'wild'      : np.exp(rng.uniform(-2,2,M+1)),
}
print(f"{'varpi':>6s} {'shape':12s} {'a_0':>10s} {'lmin(O\')':>11s} {'lmin(O\'-a0 e0e0)':>18s}")
for varpi in (0.3,1.0,5.0):
    for nm,g in shapes.items():
        a0,l0,l1 = run(varpi,g)
        flag = "  <-- 秩一漏网失效" if l1 < -1e-9 else ""
        print(f"{varpi:6.1f} {nm:12s} {a0:10.5f} {l0:11.3e} {l1:18.3e}{flag}")
