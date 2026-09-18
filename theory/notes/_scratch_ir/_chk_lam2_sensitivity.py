"""现行 12 格 DK 判据比值，对「解析 lam2 界 eps」的敏感度。
g_an(eps) := (1-eps)*h。现行 eps=0.6（对 w<=2）。
ratio_dk = x^2 E_hi^2/(g_an F_lo)，E_hi/F_lo 也随 s=res/g_an 变。
这是可行性探针，不是证书：把 eps 当外生解析输入，问「若有人证到 lam2<=eps*h，最劣格变多少」。
真值 lam2/h ~ 7e-4，故 eps=0.05 这一列远未触及物理极限。
严格性：numerical-evidence (float64)。
"""
import os, io, contextlib, importlib
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

def O(phi, N):
    n = np.arange(N); I = n[:, None]; J = n[None, :]
    pv = phi(np.arange(0, 2*N+2, dtype=float))
    T = pv[np.abs(I-J)].copy(); np.fill_diagonal(T, 0.0)
    u = 1/np.sqrt(2*n+1.)
    M = (T + pv[I+J+1])*u[:, None]*u[None, :]
    M[n, n] -= 2.0*np.concatenate(([0.], np.cumsum(pv[1:N])))/(2*n+1.)
    return M

with contextlib.redirect_stdout(io.StringIO()):
    THETA_STAR = importlib.import_module("lemmaC_compact_probe").THETA_STAR

CELLS = [(0.30,0.42,5.0),(0.42,0.60,3.6),(0.60,0.85,2.7),(0.85,1.20,2.1),
         (1.20,1.70,1.8),(1.70,2.40,1.6),(2.40,3.40,1.5),(3.40,4.80,1.42),
         (4.80,6.80,1.35),(6.80,9.60,1.30),(9.60,13.6,1.27),(13.6,20.0,1.24)]
EPS = [0.6, 0.4, 0.2, 0.1, 0.05]

def ratios(w, al, N=600):
    x = w*w
    A = O(lambda k: x/(k*k+x), N)
    Ophi = O(lambda k: k*k/(k*k+x)**3, N); Opsi = O(lambda k: 1.0/(k*k+x)**2, N)
    h = np.linalg.eigvalsh(A)[-1]
    n = np.arange(N); a = al*w
    t = (1/np.sqrt(2*n+1.))*a*a/(a*a+(n+0.5)**2); t /= np.linalg.norm(t)
    res = np.linalg.norm(A@t-(t@(A@t))*t)
    Ft = t@(Ophi@t); e_phi = np.linalg.norm(Ophi@t-Ft*t)
    Kt = Opsi@t; qpsi = t@Kt; e_psi = np.linalg.norm(Kt-qpsi*t)
    nphi = np.abs(np.linalg.eigvalsh(Ophi)).max()
    npsi = np.abs(np.linalg.eigvalsh(Opsi)).max()
    out = []
    for eps in EPS:
        # w>2 仍用仓库现行的 h-THETA_STAR 与 (1-eps)h 取较大者（都是合法解析下界）
        g = max((1.0-eps)*h, h-THETA_STAR if w > 2 else 0.0)
        s = min(res/g, 1.0)
        F_lo = (1-s*s)*Ft - 2*s*e_phi - s*s*nphi
        E_hi = e_psi + s*(abs(qpsi)+npsi)
        out.append(x*x*E_hi**2/(g*F_lo) if (F_lo > 0 and g > 0) else np.inf)
    return out

print("现行 = eps 0.60 那一列（w<=2）。表内为每格 11 点扫描的最劣 ratio_dk。")
print(f"{'cell':>15} {'a':>5} | " + " ".join(f"{'eps='+str(e):>9}" for e in EPS))
worst = [0.0]*len(EPS)
for wa, wb, al in CELLS:
    rs = np.array([ratios(w, al) for w in np.geomspace(wa, wb, 11)])
    cw = rs.max(axis=0)
    worst = [max(a_, b_) for a_, b_ in zip(worst, cw)]
    print(f"[{wa:5.2f},{wb:5.2f}] {al:5.2f} | " + " ".join(f"{v:9.4f}" for v in cw))
print(f"{'12 格最劣':>15} {'':5} | " + " ".join(f"{v:9.4f}" for v in worst))
print(f"{'相对现行':>15} {'':5} | " + " ".join(f"{v/worst[0]:9.4f}" for v in worst))
