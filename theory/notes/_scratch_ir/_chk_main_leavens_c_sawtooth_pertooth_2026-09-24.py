"""锯齿感知的 max_w c(w)：c(w) 在 N_c=floor(r*w) 跳点处不连续，
黄金分割/单峰假设在此失效。逐齿取最大，闭式 S。"""
import os; os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.linalg import eigh

def Khat(w, nd):
    n = np.arange(nd); u = 1.0/np.sqrt(2.0*n+1.0)
    i, j = np.meshgrid(n, n, indexing="ij")
    G = w*w/(w*w + np.arange(2*nd+2)**2)
    Gc = w*w/(w*w + np.arange(nd)**2)
    S = Gc[0] + 2.0*(np.cumsum(Gc) - Gc[0])
    return (G[np.abs(i-j)] + G[i+j+1])*np.outer(u,u) - np.diag(u*u*S), u

def c_of_w(w, mus, r, nd):
    K, u = Khat(w, nd); Nc = int(np.floor(r*w)); B = np.eye(nd)
    if mus > 0 and Nc > 0:
        uc = np.where(np.arange(nd) < Nc, u, 0.0); B = B + 2.0*mus*np.outer(uc,uc)
    th = float(eigh(K, B, subset_by_index=[nd-1,nd-1], eigvals_only=True)[0])
    return th/(np.pi*w)

def pertooth_max(mus, r, nd, wlo=0.6, whi=2.6, nper=13, eps=1e-9):
    best = (-1.0, None, None)
    ks = range(int(np.floor(r*wlo)), int(np.floor(r*whi))+1)
    rows = []
    for k in ks:                       # 齿 k: w in [k/r, (k+1)/r)
        a = max(wlo, k/r); b = min(whi, (k+1)/r - eps)
        if b <= a: continue
        ws = np.linspace(a, b, nper)
        cs = np.array([c_of_w(w, mus, r, nd) for w in ws])
        i = int(np.argmax(cs))
        # 齿内细化（齿内连续）
        lo2 = ws[max(i-1,0)]; hi2 = ws[min(i+1,nper-1)]
        ws2 = np.linspace(lo2, hi2, 25)
        cs2 = np.array([c_of_w(w, mus, r, nd) for w in ws2])
        i2 = int(np.argmax(cs2))
        rows.append((k, ws2[i2], cs2[i2], i2 in (0, 24)))
        if cs2[i2] > best[0]: best = (cs2[i2], ws2[i2], k)
    return best, rows

for mus, r in ((0.13, 10.0), (0.0, 10.0), (0.10, 10.0), (0.20, 10.0)):
    (cb, wb, kb), rows = pertooth_max(mus, r, 200)
    top = sorted(rows, key=lambda t: -t[2])[:5]
    print(f"mu*={mus:4.2f} r={r}:  c_max={cb:.9f}  at w*={wb:.6f}  齿 N_c={kb}")
    for k, w, c, edge in top:
        print(f"      齿 N_c={k:3d}  w={w:.6f}  c={c:.9f}  {'<-峰在齿端(跳点), 上确界不可达' if edge else ''}")
    print()
