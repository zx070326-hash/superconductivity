import numpy as np
# 从零：不 import 仓库任何模块。审 nambu_weight_derivation 的 L1 索引收缩。
# lam(k)=f(k)=1/(1+k^2)  (Einstein, varpi=1)
def lam(k): return 1.0/(1.0+np.asarray(k,dtype=float)**2)

def g_shapes(M):
    i = np.arange(M+1, dtype=float)
    return {
      'flat'          : np.ones(M+1),
      'vH_decr'       : 0.3+0.7/(1+i/20),                 # 单调不增
      'pg_incr'       : 1-0.7/(1+i/20),                   # 单调不减
      'step_pg'       : np.where(i<5, 0.05, 1.0),         # 尖锐赝能隙
      'step_vH'       : np.where(i<5, 1.0, 0.05),         # 尖锐 van Hove
    }

M = 400000                      # 源指标截断（lam~k^-2，尾部可控）
N = 60                          # 场点 n
gs = g_shapes(M)
m = np.arange(M+1)

print(f"{'shape':10s} {'n':>3s} {'S_exact':>12s} {'note_C':>12s} {'ansatzA':>12s}")
for name,g in gs.items():
    for n in (0,1,3,10):
        # 精确：原始 Matsubara 求和
        S = np.sum(lam(np.abs(n-m))*g) - np.sum(lam(n+1+m)*g)
        B_exact = S - lam(0)*g[n]                      # 挖掉与 Delta 对角相消的那格
        k = np.arange(1,n+1)
        B_noteC = 2*np.sum(g[k]*lam(k)) if n>=1 else 0.0   # note 第 9/114 行
        B_ansA  = g[n]*2*np.sum(lam(k)) if n>=1 else 0.0   # dos-escape Ansatz A
        print(f"{name:10s} {n:3d} {B_exact:12.6f} {B_noteC:12.6f} {B_ansA:12.6f}")
