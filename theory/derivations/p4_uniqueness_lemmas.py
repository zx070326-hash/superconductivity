r"""
P4 . Uniqueness half of KAY Conjecture 2 -- the two combinatorial lemmas of
`theory/derivations/P4_uniqueness_strategy.md` §2 (Lemma A, Lemma B).
==============================================================================
STATUS:  BOTH LEMMAS ARE PROVED ANALYTICALLY (elementary; no computer certificate is
needed for rigor).  This script (i) writes the proofs (below), (ii) verifies -- as a
regression gate against errors in the reduction "matrix definition -> closed form" --
that every closed form and every per-case lower bound used in the proofs agrees with
the brute-force matrix quantities for all indices <= N, and (iii) pins down the exact
sets where the inequalities are tight.

Notation.  u_n := (2n+1)^{-1/2} (n >= 0), d_n := u_n - u_{n+1} > 0.
    u is strictly decreasing and strictly convex;  d is strictly decreasing and convex
    (d_m = int_m^{m+1} (2x+1)^{-3/2} dx, and (2x+1)^{-3/2} is decreasing and convex).
    KEY IDENTITY:  u_n^2 - u_{n+1}^2 = 2 u_n^2 u_{n+1}^2 = d_n (u_n + u_{n+1}).        (*)
Kernel-generic operator (P4_conjecture2.md §6.2):  O[f] = -L1[f] + L2[f] + L3[f],
    L1[f]_{nn} = (2n+1)^{-1} sum_{k=1}^n 2 f(k),   L2[f]_{nm} = (1-delta_{nm}) f(|n-m|) u_n u_m,
    L3[f]_{nm} = f(n+m+1) u_n u_m   (L3 includes the diagonal n=m).
Partial sums  S(a,b) := sum_{m=a}^{b-1} u_m  (:= 0 if b <= a),  e_L := 1_{[0,L)}.

------------------------------------------------------------------------------
LEMMA A.   For all k >= 1, L, L' >= 1:   Q_{LL'}(delta_k) := <e_L, O[delta_k] e_{L'}> >= 0,
           with equality  IFF  L + L' <= k.
------------------------------------------------------------------------------
PROOF.  O[delta_k] is symmetric, so Q_{LL'} = Q_{L'L}; assume L <= L'.  Entrywise
    O[delta_k]_{nm} = u_n u_m [1_{|n-m|=k} + 1_{n+m+1=k}] - 2 u_n^2 1_{n=m, n>=k}.
Summing over n < L, m < L' (the diagonal entry (n,n) lies in the rectangle iff n < L):
    L2, m=n+k :  sum_{n=0}^{min(L, L'-k)-1} u_n u_{n+k}
    L2, m=n-k :  sum_{n=k}^{L-1} u_n u_{n-k}                 (n-k < L <= L' automatic)
    L3        :  sum_{n=max(0,k-L')}^{min(L,k)-1} u_n u_{k-1-n}   (all terms >= 0)
    L1        :  -2 sum_{n=k}^{L-1} u_n^2.
Since L <= L', min(L, L'-k) >= L-k, so the first sum contains sum_{n=0}^{L-k-1} u_n u_{n+k}
= sum_{j=k}^{L-1} u_{j-k} u_j.  Hence, EXACTLY,
    Q_{LL'}(delta_k) = 2 sum_{j=k}^{L-1} u_j (u_{j-k} - u_j)               [group 1]
                     + sum_{n=max(0,L-k)}^{min(L,L'-k)-1} u_n u_{n+k}      [group 2]
                     + sum_{n=max(0,k-L')}^{min(L,k)-1} u_n u_{k-1-n}      [group 3],
three sums of NON-NEGATIVE terms (group 1 by monotonicity of u).  So Q >= 0.
Tightness: Q = 0 iff all three index ranges are empty: L <= k, min(L,L'-k) <= 0 (i.e.
L' <= k) and k - L' >= L, i.e. L + L' <= k -- in which case the whole L x L' block of
O[delta_k] is identically zero (|n-m| <= max(L,L')-1 < k, n+m+1 <= L+L'-1 < k, n <= L-1 < k).
Conversely if L + L' > k then group 3 (if L <= k) or group 1 (if L > k) is non-empty and
strictly positive.                                                                   QED
Remark.  Corollary A' of the strategy note needs Q_{LL}(delta_1) >= u_0^2 = 1: here
group 3 is the single term n=0, u_0 u_0 = 1.  Remark 2: the explicit formula for Q_{LL'}
printed in §2 of the strategy note is valid only for L <= L' (the diagonal term needs
n < min(L,L') and the m=n-k term needs n-k < L'); for L > L' the "negative rows" it
mentions do not occur once the diagonal is restricted correctly -- and symmetry removes
the case anyway.

------------------------------------------------------------------------------
LEMMA B.   For all K >= 0 and L != L' (L, L' >= 1), with M := O[1_{k<=K}]:
           W_{L,L'}(K) := (M e_{L'})_{L-1} - (M e_{L'})_L >= 0,
           with equality  IFF  K = 0  or  L >= L' + K + 1.
------------------------------------------------------------------------------
PROOF.  Entrywise M_{nm} = u_n u_m [1_{1<=|n-m|<=K} + 1_{n+m+1<=K}] - 2 u_n^2 min(n,K) 1_{n=m}.
Write N' := L' and, for a row n,
    R(n) := S(max(0,n-K), min(N', n+K+1)) + S(0, min(N', K-n)),
    g(n) := (M e_{N'})_n = u_n R(n) - u_n^2 (1 + 2 min(n,K)) 1_{n<N'}
(the L2 window [n-K, n+K] cut to [0,N'), minus the m=n term; the L3 range m <= K-1-n cut
to [0,N'); the diagonal -2 min(n,K) u_n^2 present iff n < N').  W = g(n) - g(n+1), n := L-1.

CASE I:  L > L', i.e. n >= N'.  No diagonal in either row; g(n) = u_n R(n) with
    R(n) = S(max(0,n-K), N') + S(0, min(N', K-n)),  both pieces non-increasing in n.
    So g is non-increasing on n >= N' and W >= 0.  W = 0 iff R(n) = 0 iff n-K >= N',
    i.e. L >= L'+K+1 (these are the "many exact zeros" of the numerics, e.g. K=10,
    (L,L') = (12,1)).  Also g(n) >= 0 for n >= N' -- the last-row condition of Cor. B'.

CASE II: L < L', i.e. n+1 < N'.  Both rows carry the diagonal.
  II.a  n <= K-1 (so n+1 <= K):  u_n^2(1+2n) = 1, so g(n) = u_n R(n) - 1, same for n+1, and
        R(n+1) - R(n) = u_{n+K+1} 1_{n+K+1<N'} - u_{K-n-1} 1_{K-n-1<N'} <= 0
        (K-n-1 <= n+K+1 and u decreasing).  Hence W = u_n R(n) - u_{n+1} R(n+1)
        >= d_n R(n) >= d_n u_0 > 0.
  II.b  n >= K.  Now 1+2min(.,K) = 2K+1 in both rows, the L3 piece is empty, and
        R(n) = S(n-K, min(N', n+K+1)),  R(n+1) = S(n-K+1, min(N', n+K+2)).
        Define the row-n weights  c_m := d_n u_m + u_{n+1} d_m  (m >= 0): c is decreasing
        and convex in m, and by (*)  c_n = d_n(u_n+u_{n+1}) = u_n^2 - u_{n+1}^2, so the
        diagonal difference is  (2K+1)(u_n^2 - u_{n+1}^2) = (2K+1) c_n.
    II.b.1  N' >= n+K+2 (full window):
        W = u_n S(n-K,n+K+1) - u_{n+1} S(n-K+1,n+K+2) - (2K+1) c_n
          = d_n sum_{m=n-K}^{n+K} u_m + u_{n+1}(u_{n-K} - u_{n+K+1}) - (2K+1) c_n
          = sum_{m=n-K}^{n+K} c_m - (2K+1) c_n
          = sum_{j=1}^{K} [ c_{n-j} + c_{n+j} - 2 c_n ]  >= 0      (convexity of c),
        strictly > 0 for K >= 1 (u strictly convex, d_n > 0); = 0 for K = 0 (M = 0).
        Two orders cancel here:  W = 3K(K+1)(2K+1)/(2n+1)^4 + O((2n+1)^{-5})  (Stage 0).
    II.b.2  n+2 <= N' <= n+K+1 (window cut on the right; needs K >= 1):
        W(N') = d_n S(n-K, N') + u_{n+1} u_{n-K} - (2K+1) c_n  is non-decreasing in N'
        (increment d_n u_{N'} >= 0), so W(N') >= W(n+2) =: F, and with
        u_{n-K} = u_{n+1} + sum_{m=n-K}^{n} d_m,
        F = sum_{m=n-K}^{n} c_m + u_{n+1}(d_n + u_{n+1}) - (2K+1) c_n
          >= (K+1) c_n + u_n u_{n+1} - (2K+1) c_n            (c decreasing)
          = u_n u_{n+1} (1 - 2K u_n u_{n+1})  > 0,
        because n >= K gives 2K u_n u_{n+1} <= 2K/sqrt((2K+1)(2K+3)) < 1.
All cases covered (n <= K-1 / n >= K; N' >= n+K+2 / N' <= n+K+1).                    QED

Why the step structure is essential: for the non-monotone kernel delta_k the analogue
fails (k=3, L=2, L'=1 gives -0.447, reproduced in Stage 3); in the proof it is the
(2K+1)-fold multiplicity of the diagonal being EXACTLY the width of the L2 window that
makes the diagonal the "j=0 baseline" of the convexity sum in II.b.1.

------------------------------------------------------------------------------
WHAT THE SCRIPT CHECKS (float64 exhaustive to N = 500, mpmath dps=40 spot checks)
------------------------------------------------------------------------------
Stage 0  sympy: identity (*); u'' > 0 and (|u'|)'' > 0; asymptotic coefficient of
         W in case II.b.1 for K = 1..4 equals 3K(K+1)(2K+1).
Stage 1  Lemma A, all k < N, all L, L' <= N: matrix Q vs three-group closed form
         (max |diff|), symmetry, Q >= 0, tight set == {L+L' <= k}, min positive Q.
Stage 2  Lemma B, all K <= N, all L != L' <= N: matrix W vs closed form g(n)-g(n+1);
         per-case: I monotonicity of R and tight set; II.a  R(n+1) <= R(n), W >= d_n R(n);
         II.b.1 exact convexity-sum identity, all convexity terms >= 0, asymptotic ratio;
         II.b.2 W non-decreasing in N', W(n+2) = F identity, F >= u_n u_{n+1}(1-2K u_n u_{n+1}).
         Tight set == {K=0} u {L >= L'+K+1}; no zero for L < L', K >= 1.
Stage 3  Controls: delta_3 counterexample (-0.447); K = infinity gives 2uu^T - I.
Stage 4  mpmath (dps=40) recomputation of the smallest positive W's (cancellation-prone,
         ~1e-11) from the closed form: sign confirmed, float64 error << value.
Run:  python theory/derivations/p4_uniqueness_lemmas.py   (N=500, ~1 min)
"""
import sys, time
import numpy as np
import sympy as sp
import mpmath as mp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 500
T0 = time.time()
CHECKS = 0
def check(cond, msg):
    global CHECKS
    CHECKS += 1
    if not bool(cond):
        raise AssertionError("CHECK FAILED: " + msg)
    print("  ok  " + msg)

# ------------------------------------------------------------------ basic sequences
idx = np.arange(N + 1)
u = 1.0 / np.sqrt(2.0 * idx + 1.0)          # u_0 .. u_N
d = u[:-1] - u[1:]                           # d_0 .. d_{N-1}
cs = np.concatenate([[0.0], np.cumsum(u[:N])])   # cs[t] = sum_{m<t} u_m,  t = 0..N
cd = np.concatenate([[0.0], np.cumsum(d)])       # cd[t] = sum_{m<t} d_m,  t = 0..N
def Ssum(a, b):
    """vectorised S(a,b) = sum_{m=a}^{b-1} u_m (0 if b<=a); a,b integer arrays in [0,N]."""
    a = np.asarray(a); b = np.asarray(b)
    return cs[np.maximum(b, a)] - cs[a]
def Dsum(a, b):
    a = np.asarray(a); b = np.asarray(b)
    return cd[np.maximum(b, a)] - cd[a]

nn = np.arange(N)
def O_delta(k):
    """O[delta_k] on the N x N truncation (exact: entries do not depend on N)."""
    M = np.zeros((N, N)); i = np.arange(N - k)
    M[i, i + k] += u[i] * u[i + k]; M[i + k, i] += u[i] * u[i + k]
    n3 = np.arange(max(0, k - N), min(k, N)); M[n3, k - 1 - n3] += u[n3] * u[k - 1 - n3]
    dd = np.arange(k, N); M[dd, dd] -= 2.0 / (2 * dd + 1)
    return M
def O_step(K):
    """O[1_{k<=K}] on the N x N truncation (exact)."""
    A = np.abs(nn[:, None] - nn[None, :]); Sg = nn[:, None] + nn[None, :] + 1
    off = (A <= K).astype(float); np.fill_diagonal(off, 0.0)
    M = (off + (Sg <= K)) * np.outer(u[:N], u[:N])
    M[np.diag_indices(N)] -= 2.0 * np.minimum(nn, K) / (2 * nn + 1)
    return M

# ================================================================== Stage 0: sympy
print("=== Stage 0: symbolic identities")
x, eps = sp.symbols('x epsilon', positive=True)
U = (2 * x + 1) ** sp.Rational(-1, 2)
check(sp.simplify(U**2 - U.subs(x, x + 1)**2 - 2 * U**2 * U.subs(x, x + 1)**2) == 0,
      "(*)  u_n^2 - u_{n+1}^2 = 2 u_n^2 u_{n+1}^2")
check(sp.simplify(sp.diff(U, x, 2) - 3 * (2 * x + 1) ** sp.Rational(-5, 2)) == 0,
      "u''(x) = 3(2x+1)^{-5/2} > 0  (u convex);  u' < 0 (u decreasing)")
check(sp.simplify(sp.diff(-sp.diff(U, x), x, 2) - 15 * (2 * x + 1) ** sp.Rational(-7, 2)) == 0,
      "(|u'|)''(x) = 15(2x+1)^{-7/2} > 0  (d_m = int_m^{m+1}|u'| convex, decreasing)")
# asymptotics of case II.b.1:  a = 2n+1 = 1/eps,  W = sum_{m=n-K}^{n+K} c_m - (2K+1) c_n
for K in range(1, 5):
    ua = lambda j: (1 / eps + 2 * j) ** sp.Rational(-1, 2)     # u_{n+j}
    da = lambda j: ua(j) - ua(j + 1)
    c = lambda j: da(0) * ua(j) + ua(1) * da(j)
    Wsym = sum(c(j) for j in range(-K, K + 1)) - (2 * K + 1) * c(0)
    ser = sp.series(Wsym, eps, 0, 5).removeO()
    lead = sp.simplify(ser.coeff(eps, 4)); low = [sp.simplify(ser.coeff(eps, p)) for p in range(0, 4)]
    check(all(v == 0 for v in low) and lead == 3 * K * (K + 1) * (2 * K + 1),
          f"K={K}: W_II.b.1 = {lead}/(2n+1)^4 + O((2n+1)^-5)  [= 3K(K+1)(2K+1)]")

# ================================================================== Stage 1: Lemma A
print(f"=== Stage 1: Lemma A, k < {N}, L,L' <= {N}")
t1 = time.time()
Lg = np.arange(1, N + 1)[:, None]; Lpg = np.arange(1, N + 1)[None, :]
mask_le = (Lg <= Lpg)
maxdiff = 0.0; maxasym = 0.0; minQ = np.inf; minpos = np.inf; tight_ok = True; cells = 0
for k in range(1, N):
    M = O_delta(k); Q = np.cumsum(np.cumsum(M, 0), 1)           # Q[L-1, L'-1]
    maxasym = max(maxasym, np.abs(Q - Q.T).max())
    # group 1: 2 sum_{j=k}^{L-1} u_j (u_{j-k} - u_j)
    t = np.zeros(N); j = np.arange(k, N); t[j] = u[j] * (u[j - k] - u[j])
    G1 = 2.0 * np.concatenate([[0.0], np.cumsum(t)])[Lg]                  # index L (column)
    # group 2: prefix P[t] = sum_{n<t} u_n u_{n+k}
    P = np.concatenate([[0.0], np.cumsum(u[:N - k] * u[k:N])])
    hi = np.minimum(Lg, Lpg - k); lo = np.maximum(0, Lg - k)
    G2 = np.where(hi > lo, P[np.clip(hi, 0, N - k)] - P[np.clip(lo, 0, N - k)], 0.0)
    # group 3: w_n = u_n u_{k-1-n}, n in [max(0,k-N), k-1]
    w = np.zeros(k); n3 = np.arange(max(0, k - N), k); w[n3] = u[n3] * u[k - 1 - n3]
    Wp = np.concatenate([[0.0], np.cumsum(w)])
    hi3 = np.minimum(Lg, k); lo3 = np.maximum(0, k - Lpg)
    G3 = np.where(hi3 > lo3, Wp[np.clip(hi3, 0, k)] - Wp[np.clip(lo3, 0, k)], 0.0)
    Qc = G1 + G2 + G3
    maxdiff = max(maxdiff, np.abs((Qc - Q)[mask_le]).max())
    tight = (Lg + Lpg <= k)
    tight_ok &= bool(np.all((Q == 0.0) == tight))
    minQ = min(minQ, Q.min()); cells += Q.size
    if np.any(~tight): minpos = min(minpos, Q[~tight].min())
check(maxdiff < 1e-12, f"Q_LL'(delta_k) == group1+group2+group3 on L<=L' (max diff {maxdiff:.1e})")
check(maxasym < 1e-12, f"Q symmetric in (L,L') (max |Q-Q^T| {maxasym:.1e})")
check(minQ >= 0.0, f"min Q = {minQ:+.1e} >= 0 over {cells} cells")
check(tight_ok, "Q == 0 exactly IFF L+L' <= k")
check(minpos > 1e-6, f"min Q over L+L' > k is {minpos:.4e} > 0 (far above rounding)")
print(f"  Stage 1 time {time.time()-t1:.1f}s")

# ================================================================== Stage 2: Lemma B
print(f"=== Stage 2: Lemma B, K <= {N}, L != L' <= {N}")
t2 = time.time()
ng = np.arange(N - 1)[:, None]          # n = L-1, L = 1..N-1
Npg = np.arange(1, N + 1)[None, :]      # N' = L'
# global convexity of u and d on the integers (used in II.b.1): j in 1..(N-1)//2
conv_u = np.inf; conv_d = np.inf
for jj in range(1, (N - 1) // 2 + 1):
    m = np.arange(jj, N - jj)
    conv_u = min(conv_u, (u[m - jj] + u[m + jj] - 2 * u[m]).min())
    conv_d = min(conv_d, (d[m - jj] + d[m + jj] - 2 * d[m]).min())
check(conv_u > 0, f"u_{{n-j}}+u_{{n+j}}-2u_n > 0 for all 1<=j<=n<=N-1-j (min {conv_u:.2e})")
check(conv_d > 0, f"d_{{n-j}}+d_{{n+j}}-2d_n > 0 for all 1<=j<=n<=N-2-j (min {conv_d:.2e})")
maxdiff = 0.0; minW = np.inf; tight_ok = True; cells = 0
mono_I = True; tight_I = True; ok_IIa = True; id_IIb1 = 0.0; pos_IIb1 = np.inf; ratio_lo = np.inf; ratio_hi = 0.0
mono_IIb2 = True; id_IIb2 = 0.0; F_ok = True; minpos_lt = np.inf; minpos_at = None; zero_lt = 0
Ksweep = list(range(0, N + 1))
for K in Ksweep:
    M = O_step(K); Y = np.cumsum(M, axis=1)                    # Y[n, L'-1] = (M e_{L'})_n
    Wm = Y[:-1, :] - Y[1:, :]                                  # W[n, N'-1], n = L-1
    def R(m):   # m: integer array (broadcast with Npg)
        return Ssum(np.maximum(0, m - K), np.minimum(Npg, m + K + 1)) + Ssum(0 * m, np.minimum(Npg, K - m))
    def g(m):
        return u[m] * R(m) - u[m] ** 2 * (1 + 2 * np.minimum(m, K)) * (m < Npg)
    Wc = g(ng) - g(ng + 1)
    offd = (ng + 1 != Npg)
    maxdiff = max(maxdiff, np.abs((Wc - Wm)[offd]).max())
    minW = min(minW, Wm[offd].min()); cells += offd.sum()
    tight = (K == 0) | (ng >= Npg + K)
    tight_ok &= bool(np.all(((Wm == 0.0) == tight)[offd]))
    # ---- case I: n >= N'
    cI = (ng >= Npg)
    Rn, Rn1 = R(ng), R(ng + 1)
    mono_I &= bool(np.all((Rn1 <= Rn + 1e-15)[cI])) and bool(np.all((Rn >= -1e-15)[cI]))
    tight_I &= bool(np.all(((Rn == 0.0) == (ng >= Npg + K))[cI]))
    # ---- case II.a: n+1 < N', n <= K-1
    cIIa = (ng + 1 < Npg) & (ng <= K - 1)
    if cIIa.any():
        ok_IIa &= bool(np.all((Rn1 <= Rn + 1e-15)[cIIa])) and bool(np.all((Rn >= 1.0)[cIIa])) \
                  and bool(np.all((Wm >= d[ng] * Rn - 1e-13)[cIIa]))
    # ---- case II.b.1: n+1 < N', n >= K, N' >= n+K+2  (W independent of N' here)
    n1 = np.arange(K, N - 1 - K - 1)          # need n+K+2 <= N  ->  n <= N-K-2 ; and n+1<=N-1
    if n1.size:
        Wid = d[n1] * (Ssum(n1 - K, n1 + K + 1) - (2 * K + 1) * u[n1]) \
            + u[n1 + 1] * (Dsum(n1 - K, n1 + K + 1) - (2 * K + 1) * d[n1])
        Wmat = Wm[n1, N - 1]                  # N' = N >= n+K+2 for these n
        id_IIb1 = max(id_IIb1, np.abs(Wid - Wmat).max())
        if K >= 1:
            pos_IIb1 = min(pos_IIb1, Wid.min())
        if 1 <= K <= 5:                       # asymptotic regime n >> K^2 only (sanity, not used in proof)
            far = n1 >= N // 2
            asym = 3.0 * K * (K + 1) * (2 * K + 1) / (2.0 * n1[far] + 1) ** 4
            r = Wid[far] / asym; ratio_lo = min(ratio_lo, r.min()); ratio_hi = max(ratio_hi, r.max())
    # ---- case II.b.2: n+1 < N', n >= K, n+2 <= N' <= n+K+1  (K >= 1)
    if K >= 1:
        cb2 = (ng >= K) & (Npg >= ng + 2) & (Npg <= ng + K + 1)
        if cb2.any():
            # monotone in N': W(n,N') - W(n,N'-1) = d_n u_{N'-1} >= 0 within the case
            inc = (Wm[:, 1:] - Wm[:, :-1]); both = cb2[:, 1:] & cb2[:, :-1]
            mono_IIb2 &= bool(np.all((inc >= -1e-13)[both]))
            n2 = np.arange(K, N - 2)          # N' = n+2 <= N
            if n2.size:
                Fexact = d[n2] * Ssum(n2 - K, n2 + 2) + u[n2 + 1] * u[n2 - K] - (2 * K + 1) * (u[n2] ** 2 - u[n2 + 1] ** 2)
                Fmat = Wm[n2, n2 + 2 - 1]
                id_IIb2 = max(id_IIb2, np.abs(Fexact - Fmat).max())
                Flow = u[n2] * u[n2 + 1] * (1 - 2 * K * u[n2] * u[n2 + 1])
                F_ok &= bool(np.all(Fexact >= Flow - 1e-13)) and bool(np.all(Flow > 0))
    # ---- L < L' region: strictly positive for K >= 1, record min
    clt = (ng + 1 < Npg)
    if K >= 1:
        zero_lt += int((Wm[clt] <= 0.0).sum())
        mloc = Wm[clt].min()
        if mloc < minpos_lt:
            minpos_lt = mloc; ii = np.argwhere(clt & (Wm == mloc))[0]; minpos_at = (K, int(ii[0]) + 1, int(ii[1]) + 1)
check(maxdiff < 1e-12, f"W_LL'(K) == g(L-1)-g(L) closed form for all L != L' (max diff {maxdiff:.1e})")
check(minW >= 0.0, f"min W = {minW:+.1e} >= 0 over {cells} off-diagonal cells, K = 0..{N}")
check(tight_ok, "W == 0 exactly IFF (K == 0 or L >= L'+K+1)")
check(mono_I and tight_I, "case I (L>L'): R(n) >= 0 non-increasing in n; R(L-1)=0 iff L>=L'+K+1")
check(ok_IIa, "case II.a (L<L', L<=K): R(n+1) <= R(n), R(n) >= 1, W >= d_n R(n)")
check(id_IIb1 < 1e-12, f"case II.b.1: W == sum_j [c_(n-j)+c_(n+j)-2c_n] (convexity form; max diff {id_IIb1:.1e})")
check(pos_IIb1 > 0, f"case II.b.1: convexity form > 0 for K>=1 (min {pos_IIb1:.3e})")
check(0.6 < ratio_lo <= ratio_hi < 1.4, f"case II.b.1 (K<=5, n>={N//2}): W / [3K(K+1)(2K+1)/(2n+1)^4] in [{ratio_lo:.3f}, {ratio_hi:.3f}]")
check(mono_IIb2, "case II.b.2: W non-decreasing in N' on n+2 <= N' <= n+K+1")
check(id_IIb2 < 1e-12, f"case II.b.2: W(n, N'=n+2) == F closed form (max diff {id_IIb2:.1e})")
check(F_ok, "case II.b.2: F >= u_n u_{n+1} (1 - 2K u_n u_{n+1}) > 0")
check(zero_lt == 0, f"L < L', K >= 1: no zero (min W = {minpos_lt:.4e} at (K,L,L') = {minpos_at})")
print(f"  Stage 2 time {time.time()-t2:.1f}s")

# ================================================================== Stage 3: controls
print("=== Stage 3: controls")
M3 = O_delta(3); Y3 = np.cumsum(M3, axis=1); Wd = Y3[1, 0] - Y3[2, 0]        # L=2, L'=1
check(abs(Wd - (-0.4472)) < 1e-3, f"non-monotone kernel delta_3: W_(2,1) = {Wd:+.4f} < 0 (Lemma B needs the step)")
Minf = O_step(2 * N); Mref = 2 * np.outer(u[:N], u[:N]) - np.eye(N)
check(np.abs(Minf - Mref).max() < 1e-14, "K >= 2N-1 (K = infinity): O[1] = 2uu^T - I")

# ================================================================== Stage 4: mpmath
print("=== Stage 4: mpmath dps=40 recomputation of the smallest positive W (case II.b.1)")
mp.mp.dps = 40
def W_mp(K, n):
    um = lambda m: 1 / mp.sqrt(2 * m + 1)
    dm = lambda m: um(m) - um(m + 1)
    c = lambda m: dm(n) * um(m) + um(n + 1) * dm(m)
    return sum(c(m) for m in range(n - K, n + K + 1)) - (2 * K + 1) * c(n)
worst_rel = 0.0
for (K, n) in [(1, N - 3), (1, N - 4), (2, N - 4), (1, 1), (1, 5), (3, 3), (5, 200), (50, 300)]:
    if n < K or n + K + 2 > N: continue
    wmp = W_mp(K, n); wfl = float(np.cumsum(O_step(K), axis=1)[n, N - 1] - np.cumsum(O_step(K), axis=1)[n + 1, N - 1])
    asym = mp.mpf(3 * K * (K + 1) * (2 * K + 1)) / (2 * n + 1) ** 4
    rel = abs(wfl - float(wmp)) / float(wmp); worst_rel = max(worst_rel, rel)
    print(f"    K={K:3d} n={n:4d}: W_mp = {mp.nstr(wmp, 12)}  float64 = {wfl:.6e}  asym = {mp.nstr(asym, 6)}  rel.err {rel:.1e}")
    check(wmp > 0, f"    W_mp(K={K}, n={n}) > 0")
check(worst_rel < 1e-4, f"float64 vs mpmath relative error {worst_rel:.1e} << 1 (no sign ambiguity from rounding)")

print(f"\nALL {CHECKS} CHECKS PASSED   (N = {N}, total {time.time()-T0:.1f}s)")
print("Lemma A: PROVED (tight iff L+L' <= k).   Lemma B: PROVED (tight iff K=0 or L >= L'+K+1).")
