import numpy as np
# 用精确 L1（g 挂源指标 + 负频反射项）重算引理 4 的 a_n。从零构造。
varpi = 1.0
def lam(k): return varpi**2/(np.asarray(k,dtype=float)**2+varpi**2)
rho = lam(0) - np.pi*varpi/np.sinh(np.pi*varpi)      # f(0)-fhat(pi), 单 Einstein 模

M, N = 400000, 40
m = np.arange(M+1); i = np.arange(M+1, dtype=float)
shapes = {
  'flat'    : np.ones(M+1),
  'vH_decr' : 0.3+0.7/(1+i/20),
  'pg_incr' : 1-0.7/(1+i/20),
  'step_pg' : np.where(i<5, 0.05, 1.0),
  'step_vH' : np.where(i<5, 1.0, 0.05),
}
print(f"rho = {rho:.6f}\n")
print(f"{'shape':9s} {'a_0':>10s} {'min a_n (n>=1)':>15s} {'argmin':>7s} {'#neg n>=1':>10s}")
for name,g in shapes.items():
    a = np.empty(N+1)
    for n in range(N+1):
        S = np.sum(lam(np.abs(n-m))*g) - np.sum(lam(n+1+m)*g)
        B = S - lam(0)*g[n]
        a[n] = (B - rho*g[n])/(2*n+1)
    neg = np.sum(a[1:] < 0)
    print(f"{name:9s} {a[0]:10.5f} {a[1:].min():15.6f} {1+int(np.argmin(a[1:])):7d} {neg:10d}")
