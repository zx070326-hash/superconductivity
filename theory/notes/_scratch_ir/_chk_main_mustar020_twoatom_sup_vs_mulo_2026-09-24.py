"""mu_lo 界的是 sup over 两原子族，不是 Einstein 单点。
自建两原子族（固定 lam 与 M2=lam*w^2），扫是否存在成员 lam1 > mu_lo。
若存在 ==> mu_lo 是该 sup 的合法下界，下侧声音；我的 M3 下侧判据对象错。"""
import os; os.environ["OMP_NUM_THREADS"]="1"
import numpy as np
from scipy.linalg import eigh
R, MUS, ND = 10.0, 0.20, 600

def lam1_member(lams, ws, Nc, nd=ND):
    """G(k)=sum_a lam_a w_a^2/(w_a^2+k^2); 对角用 G(0) 抵消后的写法。"""
    n = np.arange(nd); u = 1.0/np.sqrt(2.0*n+1.0)
    k = np.arange(1, 2*nd+2)
    g = np.zeros(2*nd+2)
    for la, wa in zip(lams, ws):
        g[1:] += la*wa*wa/(wa*wa + k*k)
    i, j = np.meshgrid(n, n, indexing="ij")
    off = g[np.abs(i-j)] + g[i+j+1]; np.fill_diagonal(off, 0.0)
    cum = np.concatenate([[0.0], np.cumsum(g[1:nd+1])])
    M = np.outer(u,u)*(off + np.diag(g[2*n+1] - 2.0*cum[:nd]))
    if MUS > 0 and Nc > 0:
        uc = np.where(n < Nc, u, 0.0); M = M - 2.0*MUS*np.outer(uc, uc)
    return float(eigh(M, subset_by_index=[nd-1,nd-1], eigvals_only=True)[0])

CASES = [(2.0,1.100033,1.021190594),(3.2953,0.731708,1.020081705),
         (4.0,0.640491,1.021883712),(6.0,0.495741,1.025180308)]
for lam, w, ml in CASES:
    Nc = int(np.floor(R*w)); M2 = lam*w*w
    ein = lam1_member([lam],[w],Nc)
    best = (ein, 1.0, w, None)
    for t in np.concatenate([np.linspace(0.02,0.98,25), [0.995,0.999]]):
        amax = w/np.sqrt(t)
        for a in np.linspace(1e-3, amax*0.9999, 40):
            w2sq = (w*w - t*a*a)/(1.0-t)
            if w2sq <= 0: continue
            v = lam1_member([lam*t, lam*(1-t)], [a, np.sqrt(w2sq)], Nc)
            if v > best[0]: best = (v, t, a, np.sqrt(w2sq))
    sup, t, a, b = best
    print(f"lam={lam:<7} w={w:.6f} Nc={Nc}  Einstein lam1={ein:.9f}   mu_lo={ml:.9f}")
    print(f"   两原子扫最优 lam1={sup:.9f}  (t={t:.3f}, w1={a:.4f}, w2={b if b is None else round(b,4)})")
    print(f"   sup >= mu_lo ? {'是 ==> 下侧声音, 我的 M3 对象错' if sup >= ml else '** 否, 需查 **'}"
          f"   (余量 {sup-ml:+.2e};  两原子比 Einstein 高 {sup-ein:+.2e})")
