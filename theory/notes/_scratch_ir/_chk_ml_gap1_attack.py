# -*- coding: utf-8 -*-
"""缺口 1 定向攻击：O_pair 的顶本征向量是否真在锥 C^B 内。
C = {v: v>=0 且非增}。O_pair 分块 = [[O[L11],O[L12]],[O[L12],O[L22]]]，每个 Lij 非负非增。
结构提示：主带支配时 psi_2 ~ O[L12] psi_1 / mu，而 O[f] 不保锥 ==> 往 c12 小、两带 varpi 悬殊处打。
"""
import os; os.environ.setdefault("OMP_NUM_THREADS","1")
import numpy as np

def O(fv, N):
    """fv[k] for k=0..2N+1; fv[0] 不用。"""
    n = np.arange(N); u = 1.0/np.sqrt(2*n+1.0)
    T = fv[np.abs(n[:,None]-n[None,:])].copy(); np.fill_diagonal(T, 0.0)
    M = (T + fv[n[:,None]+n[None,:]+1])*np.outer(u,u)
    M[n,n] -= 2.0*np.concatenate(([0.], np.cumsum(fv[1:N])))/(2*n+1.0)
    return M

def kern(c, w, N):           # Einstein: c * w^2/(k^2+w^2)，非负非增
    k = np.arange(2*N+2, dtype=float); return c*w*w/(k*k+w*w)

def viol(psi, N, B):
    """返回 (最负分量, 最负一阶差分)；两者 >= -tol 才算在锥内。全局符号取 |.| 最大的那格为正。"""
    p = psi.reshape(B, N).copy()
    i, j = np.unravel_index(np.argmax(np.abs(p)), p.shape)
    if p[i, j] < 0: p = -p
    return p.min(), np.diff(p, axis=1).min(), p

def top(blocks, N, B):
    M = np.zeros((B*N, B*N))
    for i in range(B):
        for j in range(B): M[i*N:(i+1)*N, j*N:(j+1)*N] = blocks[i][j]
    e, V = np.linalg.eigh(M); return e[-1], V[:, -1]

N, B = 300, 2
rng = np.random.default_rng(7)
print("=== A. 复现 agent 的随机采样（同量级，400 组）===")
worst = (1e9, None)
for t in range(400):
    c11, c22 = rng.uniform(0.2, 3.0, 2); c12 = rng.uniform(0.0, np.sqrt(c11*c22))
    w11, w22, w12 = rng.uniform(0.3, 20.0, 3)
    bl = [[O(kern(c11,w11,N),N), O(kern(c12,w12,N),N)],
          [O(kern(c12,w12,N),N), O(kern(c22,w22,N),N)]]
    _, psi = top(bl, N, B); a, b, _ = viol(psi, N, B)
    if min(a,b) < worst[0]: worst = (min(a,b), (c11,c22,c12,w11,w22,w12))
print(f"  最劣违反量 = {worst[0]:.3e}   参数 {np.round(worst[1],4)}")

print("\n=== B. 定向攻击：c12 小 + 两带 varpi 悬殊 ===")
rows = []
for c12 in (1e-3, 1e-2, 0.05, 0.2):
    for w1, w2 in ((20.0, 0.3), (20.0, 1.0), (0.3, 20.0), (15.0, 0.5)):
        for wx in (0.3, 1.0, 5.0, 20.0):
            bl = [[O(kern(2.0,w1,N),N), O(kern(c12,wx,N),N)],
                  [O(kern(c12,wx,N),N), O(kern(0.05,w2,N),N)]]
            lam, psi = top(bl, N, B); a, b, p = viol(psi, N, B)
            rows.append((min(a,b), c12, w1, w2, wx, lam, a, b))
rows.sort()
print(f"{'viol':>12} {'c12':>7} {'w1':>6} {'w2':>6} {'wx':>6} {'lam1':>9} {'min psi':>12} {'min diff':>12}")
for r in rows[:8]: print(f"{r[0]:>12.3e} {r[1]:>7.3g} {r[2]:>6.1f} {r[3]:>6.1f} {r[4]:>6.1f} {r[5]:>9.4f} {r[6]:>12.3e} {r[7]:>12.3e}")
