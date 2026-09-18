import numpy as np
rng = np.random.default_rng(20260916)
M, N = 400000, 300
m = np.arange(M+1); i = np.arange(M+1, dtype=float)

def build(varpi):
    lam = lambda k: varpi**2/(np.asarray(k,dtype=float)**2+varpi**2)
    rho = lam(0) - np.pi*varpi/np.sinh(np.pi*varpi)
    return lam, rho

shapes = {}
for w in (3,10,40,120):                      # 赝能隙台阶宽度
    for d in (0.02,0.2):                     # 谷底深度
        shapes[f'pg_step_w{w}_d{d}']  = np.where(i<w, d, 1.0)
        shapes[f'vH_step_w{w}_d{d}']  = np.where(i<w, 1.0, d)
shapes['wild1'] = np.exp(rng.uniform(-2,2,M+1))      # 非单调，量级跨 e^4
shapes['wild2'] = np.exp(rng.uniform(-3,3,M+1))
shapes['osc']   = 1.0+0.9*np.sin(i/3.0)              # 强振荡
shapes['pg_deep']= np.where(i<8, 1e-3, 1.0)          # 近乎全空的赝能隙

worst = {}
for varpi in (0.3, 1.0, 5.0, 20.0):
    lam, rho = build(varpi)
    L = lam(np.arange(0, M+N+2))                     # 预算 lam 表
    for name,g in shapes.items():
        a = np.empty(N+1)
        for n in range(N+1):
            S = np.dot(L[np.abs(n-m)], g) - np.dot(L[n+1+m], g)
            a[n] = (S - L[0]*g[n] - rho*g[n])/(2*n+1)
        key=(varpi,name)
        worst[key]=(a[0], a[1:].min(), int(np.sum(a[1:]<0)))

bad=[k for k,v in worst.items() if v[1]<0]
print(f"总组合 {len(worst)}；min a_n(n>=1)<0 的组合数：{len(bad)}")
print(f"全局最小 min_n>=1 a_n = {min(v[1] for v in worst.values()):.6e}")
print(f"a_0 最大值（应 <=0 才保秩一漏网符号）= {max(v[0] for v in worst.values()):.6e}")
for k in bad[:10]: print("  NEG:", k, worst[k])
print("\n抽样：")
for k in [(1.0,'pg_step_w120_d0.02'),(1.0,'pg_deep'),(1.0,'wild2'),(0.3,'osc'),(20.0,'wild1')]:
    print(f"  {str(k):34s} a_0={worst[k][0]:+.5f}  min a_n={worst[k][1]:+.6e}")
