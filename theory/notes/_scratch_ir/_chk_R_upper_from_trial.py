"""驳 routeRpM1 的「精确 R 需要全谱、试探向量算不出」。

R = <r,(h-A)^+ r> 只需解线性方程 (h-A)y = r（限制在 v^perp），不需全谱分解。
试探向量的变分原理 R = max_y [2<r,y> - <y,(h-A)y>] 给的是 R 的**下界**（方向错）。
但对**线性方程组**做试探可给上界：rho := r - (h-A)y，在 v^perp 上恒等式

    R = <y,(h-A)y> + 2<y,rho> + <rho,(h-A)^+ rho>
      <= <y,(h-A)y> + 2<y,rho> + ||rho||^2 / g          (UB)

y=0 时 UB 恰为 E^2/g = DK 界 => 只会更紧，不会更差。
方向审计：唯一的放缩是 <rho,(h-A)^+rho> <= ||rho||^2/g，用 g_lo <= g 即为上估；
这是 R 的**上估** => |C|=2x^3 R 的上估 => 判据「|C|<E」的保守侧 => 只支持 PASS。

严格性：本脚本 [numerical-evidence]（float64，精确 h/v）。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

def O(phi_vals, N):
    n = np.arange(N); I = n[:, None]; J = n[None, :]
    T = phi_vals[np.abs(I - J)].copy(); np.fill_diagonal(T, 0.0)
    Hk = phi_vals[I + J + 1]
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    M = (T + Hk) * u[:, None] * u[None, :]
    cum = np.concatenate(([0.0], np.cumsum(phi_vals[1:N])))
    M[n, n] -= 2.0 * cum / (2.0 * n + 1.0)
    return M

def setup(varpi, N):
    x = varpi ** 2
    k = np.arange(0, 2 * N, dtype=float); d = k * k + x
    A = O(x / d, N)
    w, V = np.linalg.eigh(A)
    h, v, lam2 = w[-1], V[:, -1], w[-2]
    g = h - lam2
    Opsi = O(1.0 / d ** 2, N)
    yv = Opsi @ v
    r = yv - v * (v @ yv)
    c = V.T @ r
    m = np.ones(len(w), bool); m[-1] = False
    R = float(np.sum(c[m] ** 2 / (h - w[m])))       # 真值（用全谱，仅作参照）
    return x, A, h, v, g, r, R

def ub_from_y(A, h, v, g, r, y):
    """UB = <y,(h-A)y> + 2<y,rho> + ||rho||^2/g，全部投影到 v^perp"""
    y = y - v * (v @ y)
    My = h * y - A @ y
    My = My - v * (v @ My)
    rho = r - My
    return float(y @ My + 2.0 * (y @ rho) + (rho @ rho) / g)

print(f"{'varpi':>6} {'R_true':>12} {'DK=E^2/g':>12} {'R/DK':>7} | "
      f"{'m=1':>7} {'m=2':>7} {'m=3':>7} {'m=5':>7} {'m=8':>7}   (UB/R_true)")
for varpi in (0.3, 1.0, 2.0, 3.0, 8.5, 20.0):
    N = 600
    x, A, h, v, g, r, R = setup(varpi, N)
    DK = float(r @ r) / g
    out = []
    # Krylov-Galerkin 求 y，再代入 UB
    for m in (1, 2, 3, 5, 8):
        K = np.zeros((len(r), m)); q = r.copy()
        for j in range(m):
            q = q - v * (v @ q)
            K[:, j] = q
            q = h * q - A @ q
        G = K.T @ (h * K - A @ K)
        b = K.T @ r
        coef = np.linalg.solve(G + 1e-300 * np.eye(m), b)
        out.append(ub_from_y(A, h, v, g, r, K @ coef) / R)
    print(f"{varpi:6.1f} {R:12.6e} {DK:12.6e} {R/DK:7.4f} | "
          + " ".join(f"{o:7.4f}" for o in out))
