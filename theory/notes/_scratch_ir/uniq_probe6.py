# Exhaustive check of the two combinatorial lemmas behind the uniqueness strategy (P- and w-independent):
#  (A) Q_{LL'}(delta_k) >= 0  for all k,L,L'  (=> <v,O[f]v> >= 0 for every f>=0 and every v in the cone C)
#  (B) Kamke: step-basis off-diagonal of O[1_{k<=K}] >= 0 for all K, L != L'  (=> e^{tO[F]} preserves C for F>=0 decreasing)
#  plus (B') the constant kernel 1 (K=infty): O[1] = 2uu^T - I.
import numpy as np
N = 500; n = np.arange(N); u = 1/np.sqrt(2*n+1)
def O_step(K):
    A = np.abs(n[:,None]-n[None,:]); S = n[:,None]+n[None,:]+1
    off = (A <= K)*1.0; np.fill_diagonal(off, 0.0); M = (off + (S <= K))*np.outer(u,u)
    cnt = np.minimum(n, K); M[np.diag_indices(N)] -= 2*cnt/(2*n+1); return M
worst = (np.inf, None)
for K in list(range(1, 120)) + [150, 200, 300, 499, 10**6]:
    M = O_step(K); Y = np.cumsum(M, axis=1); Wm = Y[:-1,:] - Y[1:,:]
    L = np.arange(1, N)[:,None]; Lp = np.arange(1, N+1)[None,:]
    off = np.where(L != Lp, Wm, np.inf); m = off.min()
    if m < worst[0]: worst = (m, (K,)+tuple(int(x)+1 for x in np.unravel_index(off.argmin(), off.shape)))
print(f"(B) step kernels 1_{{k<=K}}, N={N}: min Kamke off-diagonal = {worst[0]:+.3e} at (K,L,L')={worst[1]}")
# tight cases: count exact zeros for K=10
M = O_step(10); Y = np.cumsum(M, axis=1); Wm = Y[:-1,:] - Y[1:,:]; L = np.arange(1, N)[:,None]; Lp = np.arange(1, N+1)[None,:]
off = np.where(L != Lp, Wm, np.inf); print("    K=10: #off-diag entries with |W|<1e-15:", int(np.sum(np.abs(off) < 1e-15)), " (of", off.size, ")  examples:", [tuple(int(x)+1 for x in z) for z in np.argwhere(np.abs(off) < 1e-15)[:6]])
# (A) with weights c=1 exhaustively to N=500 (already done to 400): report min and the zero set structure
worstA = (np.inf, None)
for k in range(1, N):
    M = np.zeros((N, N)); i = np.arange(N-k); M[i, i+k] += u[i]*u[i+k]; M[i+k, i] += u[i]*u[i+k]
    for nn in range(0, min(k, N)):
        mm = k-1-nn
        if 0 <= mm < N: M[nn, mm] += u[nn]*u[mm]
    d = np.arange(k, N); M[d, d] -= 2/(2*d+1)
    Q = np.cumsum(np.cumsum(M, 0), 1); m = Q.min()
    if m < worstA[0]: worstA = (m, (k,)+tuple(int(x)+1 for x in np.unravel_index(Q.argmin(), Q.shape)))
print(f"(A) Q_LL'(delta_k) >= 0, k<N={N}: min = {worstA[0]:+.3e} at (k,L,L')={worstA[1]}")
