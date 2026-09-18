"""独立复核 lemmaC 的核心断言：M1 试探向量 tau_n = u_n a^2/(a^2+(n+1/2)^2)
在 a = omega_bar 处残差 0.19-0.26，放开宽度 a = alpha*omega_bar 后降到 0.004-0.015。
算子从零自建（Einstein 谱），不 import 任何门禁脚本。"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

def build_A(w, N):
    """Einstein: F(k) = w^2/(k^2+w^2), F(0)=1.  仓库约定的 O[f] 矩阵。"""
    k = np.arange(0, 2*N+2, dtype=float)
    F = w*w/(k*k + w*w); F[0] = 1.0
    u = 1.0/np.sqrt(2*np.arange(N)+1.0)
    # L1 对角: (2n+1)^-1 * sum_{k=1}^{n} 2 F(k)   (上限含 k=n)
    L1 = np.zeros(N); acc = 0.0
    for n in range(1, N):
        acc += 2*F[n]; L1[n] = acc/(2*n+1)
    i = np.arange(N)[:, None]; j = np.arange(N)[None, :]
    A = (F[np.abs(i-j)] + F[i+j+1]) * u[i]*u[j]
    np.fill_diagonal(A, -L1 + F[2*np.arange(N)+1]*u**2)
    return A, u

def residual(A, tau):
    tau = tau/np.linalg.norm(tau)
    rho = tau @ (A @ tau)
    return np.linalg.norm(A @ tau - rho*tau), rho

print(f"{'omega':>7} {'lam1':>10} {'gap':>9} | {'a=omega':>9} {'best a/omega':>13} {'res_best':>10} {'ratio':>7}")
print("-"*78)
for w in (0.3, 0.5, 1.0, 2.0, 4.0, 8.48, 15.0, 20.0):
    N = max(400, int(60*w))
    A, u = build_A(w, N)
    ev = np.linalg.eigvalsh(A)
    lam1, lam2 = ev[-1], ev[-2]
    n_half = np.arange(N) + 0.5
    def res_at(alpha):
        a = alpha*w
        return residual(A, u * a*a/(a*a + n_half**2))[0]
    r_base = res_at(1.0)
    alphas = np.linspace(0.8, 12.0, 400)
    rs = np.array([res_at(al) for al in alphas])
    k = int(np.argmin(rs))
    print(f"{w:7.2f} {lam1:10.6f} {lam1-lam2:9.5f} | {r_base:9.4f} {alphas[k]:13.2f} {rs[k]:10.5f} {r_base/rs[k]:7.1f}x")
