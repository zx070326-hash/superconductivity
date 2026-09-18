"""必要条件检验：lam2(O[f_varpi]) <= Theta* = 0.6349287917881079 是否对一切 varpi 成立？
p4_concavity_q_route.py:343 用 max(0.42h, h-Theta*) 作解析谱隙下界，varpi 从 0.2 到 40。
若某处 lam2 > Theta*，那是已合并结果里的 fail-open。numerical-evidence (float64)。
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
TH = 0.6349287917881079

def O(phi, N):
    n = np.arange(N); I = n[:, None]; J = n[None, :]
    pv = phi(np.arange(0, 2*N+2, dtype=float))
    T = pv[np.abs(I-J)].copy(); np.fill_diagonal(T, 0.0)
    u = 1/np.sqrt(2*n+1.)
    M = (T + pv[I+J+1])*u[:, None]*u[None, :]
    M[n, n] -= 2.0*np.concatenate(([0.], np.cumsum(pv[1:N])))/(2*n+1.)
    return M

print(f"{'varpi':>7} {'N':>6} {'h':>9} {'lam2':>11} {'lam2/Theta*':>12} "
      f"{'0.42h':>8} {'h-Th*':>8} {'用哪个':>8}")
for w in (0.2, 0.5, 1.0, 1.55, 2.0, 3.4, 5.0, 10.0, 20.0, 40.0):
    x = w*w
    prev = None
    for N in (800, 1600, 3200):
        ev = np.linalg.eigvalsh(O(lambda k: x/(k*k+x), N))
        h, l2 = ev[-1], ev[-2]
        if prev is not None and abs(l2-prev) > 1e-6*max(abs(l2), 1e-12):
            note = f"  [N未收敛 {prev:.8f}->{l2:.8f}]"
        else:
            note = ""
        prev = l2
    a, b = 0.42*h, h-TH
    print(f"{w:7.2f} {3200:6d} {h:9.4f} {l2:11.4e} {l2/TH:12.4e} "
          f"{a:8.4f} {b:8.4f} {'h-Th*' if b>a else '0.42h':>8}{note}")
print(f"\n全部 lam2 是否 <= Theta*={TH}: 见上列 lam2/Theta* (须 < 1)")
