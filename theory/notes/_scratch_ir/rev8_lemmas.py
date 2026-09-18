"""rev8: independent brute-force + symbolic check of Lemma A / Lemma B (uniqueness strategy).
Own construction of O[f] from the entrywise definition; nothing imported from derivations."""
import numpy as np, sympy as sp, time, sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
t0 = time.time()
n = np.arange(NMAX + 2)
u = (2.0 * n + 1.0) ** -0.5

def O_of_f(f, N):
    """O[f] N x N from  O[f]_{nm} = -delta_{nm} 2u_n^2 sum_{k<=n} f(k) + (1-delta_{nm}) f(|n-m|) u_n u_m + f(n+m+1) u_n u_m.
    f: array with f[k] for k = 0..2N (f[0] unused)."""
    i = np.arange(N)
    A = np.abs(i[:, None] - i[None, :]); S = i[:, None] + i[None, :] + 1
    M = f[A] * np.outer(u[:N], u[:N]); np.fill_diagonal(M, 0.0)
    M += f[S] * np.outer(u[:N], u[:N])
    cf = np.cumsum(f[:N + 1]) - f[0]              # sum_{k=1}^{n} f(k)
    M[i, i] -= 2.0 * u[:N] ** 2 * cf[:N]
    return M

# ---------------- Lemma A: Q_{LL'}(delta_k) = sum_{n<L,m<L'} O[delta_k]_{nm} >= 0, = 0 iff L+L' <= k
N = NMAX
minQ = np.inf; tight_ok = True; minpos = np.inf; worst = None
for k in range(1, N + 1):
    f = np.zeros(2 * N + 2); f[k] = 1.0
    Q = np.cumsum(np.cumsum(O_of_f(f, N), 0), 1)          # Q[L-1, L'-1]
    L = np.arange(1, N + 1)
    tight = (L[:, None] + L[None, :] <= k)
    if not np.all(Q[tight] == 0.0): tight_ok = False; worst = ('nonzero on tight', k)
    if np.any(Q[~tight] <= 0.0): tight_ok = False; worst = ('nonpos off tight', k, Q[~tight].min())
    minQ = min(minQ, Q.min())
    if np.any(~tight): minpos = min(minpos, Q[~tight].min())
print(f"[A] k,L,L' <= {N}: min Q = {minQ:+.3e}; tight set == {{L+L'<=k}} exactly: {tight_ok} {worst or ''}; min Q off tight = {minpos:.4e}")
# closed form (three groups) vs brute force for L<=L'
maxdiff = 0.0
for k in [1, 2, 3, 5, 17, 100, 299]:
    f = np.zeros(2 * N + 2); f[k] = 1.0
    Q = np.cumsum(np.cumsum(O_of_f(f, N), 0), 1)
    for L in range(1, N + 1, 7):
        for Lp in range(L, N + 1, 11):
            g1 = 2 * sum(u[j] * (u[j - k] - u[j]) for j in range(k, L))
            g2 = sum(u[m] * u[m + k] for m in range(max(0, L - k), min(L, Lp - k)))
            g3 = sum(u[m] * u[k - 1 - m] for m in range(max(0, k - Lp), min(L, k)))
            maxdiff = max(maxdiff, abs(Q[L - 1, Lp - 1] - (g1 + g2 + g3)))
print(f"[A] three-group closed form vs brute force (L<=L', sampled): max diff {maxdiff:.2e}")
# Q_LL(delta_1) >= 1
f = np.zeros(2 * N + 2); f[1] = 1.0
Q1 = np.cumsum(np.cumsum(O_of_f(f, N), 0), 1)
print(f"[A] min_L Q_LL(delta_1) = {np.diag(Q1).min():.6f} (>= 1 claimed)")

# ---------------- Lemma B: W_{L,L'}(K) = (M e_{L'})_{L-1} - (M e_{L'})_L >= 0 for L != L', M = O[1_{k<=K}]
# also the last-row condition (M e_{L'})_{N-1} >= 0 for L' < N, and the delta_3 failure.
minW = np.inf; tightB_ok = True; minpos_lt = np.inf; at = None; lastrow_min = np.inf
Ls = np.arange(1, N + 1)
for K in range(0, N + 1):
    f = np.zeros(2 * N + 2); f[1:K + 1] = 1.0
    M = O_of_f(f, N + 1)                       # (N+1) x (N+1) so that (M e_{L'})_L exists for L = N
    Y = np.cumsum(M, axis=1)                   # Y[n, L'-1] = (M e_{L'})_n
    W = Y[:-1, :N] - Y[1:, :N]                 # W[L-1, L'-1], L = 1..N, L' = 1..N
    off = (Ls[:, None] != Ls[None, :])
    minW = min(minW, W[off].min())
    tight = (K == 0) | (Ls[:, None] >= Ls[None, :] + K + 1)
    if not np.all((W[off] == 0.0) == tight[off]): tightB_ok = False
    lt = (Ls[:, None] < Ls[None, :])
    if K >= 1:
        m = W[lt].min()
        if m < minpos_lt:
            minpos_lt = m; ii = np.argwhere(lt & (W == m))[0]; at = (K, ii[0] + 1, ii[1] + 1)
    # last-row condition for every truncation size Nt <= N: (M e_{L'})_{Nt-1} >= 0 for L' < Nt
    Yl = Y[:N, :N]  # Yl[n, L'-1]
    mask = (np.arange(N)[:, None] >= Ls[None, :])   # n >= L'  <=> L' < n+1 = Nt
    lastrow_min = min(lastrow_min, Yl[mask].min())
print(f"[B] K <= {N}, L != L' <= {N}: min W = {minW:+.3e}; tight set == {{K=0 or L>=L'+K+1}} exactly: {tightB_ok}")
print(f"[B] L<L', K>=1: min W = {minpos_lt:.4e} at (K,L,L') = {at}")
print(f"[B] last-row condition (M e_L')_{{n}} for n >= L': min = {lastrow_min:+.3e} (>= 0 needed)")
f = np.zeros(2 * N + 2); f[3] = 1.0
M3 = O_of_f(f, 6); Y3 = np.cumsum(M3, axis=1)
print(f"[B] control: delta_3, (L,L')=(2,1): W = {Y3[1,0]-Y3[2,0]:+.4f} (expected -0.4472 = -1/sqrt5)")

# ---------------- symbolic: II.b.1 and II.b.2 identities with FREE symbols u_m (no property of u used)
for K in range(0, 6):
    n_ = K + 3
    us = sp.symbols(f'u0:{n_+K+3}')
    d = lambda m: us[m] - us[m + 1]
    c = lambda m: d(n_) * us[m] + us[n_ + 1] * d(m)
    Sfun = lambda a, b: sum(us[m] for m in range(a, b))
    # II.b.1 : N' >= n+K+2, direct W from g(n)-g(n+1) with diagonal (2K+1) u^2
    g = lambda r: us[r] * Sfun(r - K, r + K + 1) - us[r] ** 2 * (1 + 2 * K)    # R(r) includes m=r, diag absorbs it: u_r R(r) - u_r^2(1+2K)
    Wdirect = g(n_) - g(n_ + 1)
    Wconv = sum(c(n_ - j) + c(n_ + j) - 2 * c(n_) for j in range(1, K + 1))
    # note: identity needs c_n = u_n^2 - u_{n+1}^2 which is d_n(u_n+u_{n+1}) -- true for free symbols
    ok1 = sp.expand(Wdirect - Wconv) == 0
    # II.b.2 : N' = n+2 ; F = sum_{m=n-K}^{n} c_m + u_n u_{n+1} - (2K+1) c_n
    if K >= 1:
        gc = lambda r, Np: us[r] * Sfun(r - K, min(Np, r + K + 1)) - us[r] ** 2 * (1 + 2 * K)
        Fdirect = gc(n_, n_ + 2) - gc(n_ + 1, n_ + 2)
        Fform = sum(c(m) for m in range(n_ - K, n_ + 1)) + us[n_] * us[n_ + 1] - (2 * K + 1) * c(n_)
        ok2 = sp.expand(Fdirect - Fform) == 0
        # monotone in N': W(N'+1) - W(N') = d_n u_{N'}
        inc = sp.expand(gc(n_, n_ + 3) - gc(n_ + 1, n_ + 3) - Fdirect - d(n_) * us[n_ + 2]) == 0
    else:
        ok2 = inc = True
    print(f"[sym] K={K}: II.b.1 convexity identity {ok1}; II.b.2 F identity {ok2}; N'-increment = d_n u_N' {inc}")
# key identity and convexity of d, c
x = sp.symbols('x', positive=True)
U = (2 * x + 1) ** sp.Rational(-1, 2)
print("[sym] (*) u_n^2-u_{n+1}^2 = 2u_n^2u_{n+1}^2:", sp.simplify(U**2 - U.subs(x, x+1)**2 - 2*U**2*U.subs(x, x+1)**2) == 0)
print("[sym] |u'| = (2x+1)^{-3/2}:", sp.simplify(-sp.diff(U, x) - (2*x+1)**sp.Rational(-3, 2)) == 0,
      "; |u'|'' > 0:", sp.simplify(sp.diff((2*x+1)**sp.Rational(-3, 2), x, 2)))
m = np.arange(1, NMAX)
dd = u[:-1] - u[1:]
print(f"[num] d_m decreasing: {np.all(dd[1:] < dd[:-1])}; d second difference min = {(dd[m-1]+dd[m+1]-2*dd[m]).min():.3e}; u second diff min = {(u[m-1]+u[m+1]-2*u[m]).min():.3e}")
print(f"[num] max_K 2K u_K u_{{K+1}} over K<=300 = {max(2*K*u[K]*u[K+1] for K in range(1,301)):.6f} (<1)")
print(f"total {time.time()-t0:.1f}s")
