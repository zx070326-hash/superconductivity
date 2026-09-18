"""判决实验：在 s=0（真 v）极限下，routeRpM1 的 ||rho|| 三角放缩是否单独扼杀全部收益。

三个上界，全部用同一个真 v、同一个 y 搜索空间（m 维 Krylov）：
  UB_exact(y)    = -<y,Ty> + 2<y,r> + ||r-Ty||^2 / g        <- 主会话恒等式（精确第三项）
  UB_tri(y)      = -<y,Ty> + 2<y,r> + (||r||+||Ty||)^2 / g  <- routeRpM1 的三角放缩，s=0
  UB_tri_abs(y)  = -<y,Ty> + 2|<y,r>| + (||r||+||Ty||)^2/g  <- 再叠加它的 abs 缺陷

三者都是 R 的合法上界（单边：只支持 PASS）。看哪一个能降到 DK=E^2/g 以下。
严格性：numerical-evidence（float64）。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import minimize

def O(phi, N):
    n = np.arange(N); I = n[:, None]; J = n[None, :]
    k = np.arange(0, 2 * N + 2, dtype=float)
    pv = phi(k)
    T = pv[np.abs(I - J)].copy(); np.fill_diagonal(T, 0.0)
    Hk = pv[I + J + 1]
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    M = (T + Hk) * u[:, None] * u[None, :]
    cum = np.concatenate(([0.0], np.cumsum(pv[1:N])))
    M[n, n] -= 2.0 * cum / (2.0 * n + 1.0)
    return M

def setup(w, N=600):
    x = w * w
    A = O(lambda k: x / (k * k + x), N)
    Opsi = O(lambda k: 1.0 / (k * k + x) ** 2, N)
    ev, EV = np.linalg.eigh(A)
    h, v, g = ev[-1], EV[:, -1], ev[-1] - ev[-2]
    r = Opsi @ v; r = r - v * (v @ r)
    c = EV.T @ r
    R = float(np.sum(c[:-1] ** 2 / (h - ev[:-1])))
    return dict(w=w, A=A, h=h, v=v, g=g, r=r, E=np.linalg.norm(r), R=R, N=N)

def Top(s, y):                       # T = h I - A  作用
    return s["h"] * y - s["A"] @ y

def krylov(s, m):
    """y0 = r/||r||, 之后 y_{k+1} = T y_k 正交化（对 v 与前序）"""
    B = []
    z = s["r"].copy()
    for _ in range(m):
        z = z - s["v"] * (s["v"] @ z)
        for b in B:
            z = z - b * (b @ z)
        nz = np.linalg.norm(z)
        if nz < 1e-13:
            break
        b = z / nz; B.append(b); z = Top(s, b)
    return np.array(B)

def ubs(s, y):
    Ty = Top(s, y)
    yr = float(y @ s["r"])
    yTy = float(y @ Ty)
    rho = s["r"] - Ty
    tri = (s["E"] + np.linalg.norm(Ty)) ** 2 / s["g"]
    return (-yTy + 2 * yr + float(rho @ rho) / s["g"],
            -yTy + 2 * yr + tri,
            -yTy + 2 * abs(yr) + tri)

def best(s, B, idx):
    if len(B) == 0:
        return ubs(s, np.zeros(s["N"]))[idx]
    f = lambda cc: ubs(s, B.T @ cc, )[idx]
    b0 = ubs(s, np.zeros(s["N"]))[idx]
    out = b0
    for x0 in (np.zeros(len(B)), np.ones(len(B)) * 1e-3, np.ones(len(B)) * (1.0 / s["g"])):
        rr = minimize(f, x0, method="Nelder-Mead",
                      options=dict(xatol=1e-12, fatol=1e-16, maxiter=20000, maxfev=20000))
        out = min(out, float(rr.fun))
    return out

print("s=0（真 v）。列 = 上界/DK，DK = E^2/g。<1 才叫有收益。")
print(f"{'w':>6} {'R/DK':>8} | {'exact m=2':>10} {'exact m=5':>10} | "
      f"{'tri m=2':>9} {'tri m=5':>9} | {'tri+abs m=5':>11}")
for w in (0.3, 0.6, 1.0, 2.0, 3.0, 8.5, 20.0):
    s = setup(w)
    DK = s["E"] ** 2 / s["g"]
    B5 = krylov(s, 5); B2 = B5[:2]
    e2, e5 = best(s, B2, 0), best(s, B5, 0)
    t2, t5 = best(s, B2, 1), best(s, B5, 1)
    a5 = best(s, B5, 2)
    print(f"{w:6.2f} {s['R']/DK:8.4f} | {e2/DK:10.4f} {e5/DK:10.4f} | "
          f"{t2/DK:9.4f} {t5/DK:9.4f} | {a5/DK:11.4f}")
