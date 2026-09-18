"""
gamma_frobenius_C.py — CERTIFIED Frobenius norm bound for the C block.

WHY THIS FILE
-------------
Section 3.2 of gamma_tail_bound_audit_v2_2026-09-16.md identifies ||C|| as the
remaining loose piece: the weighted-Schur certificate gives ||C|| <= 0.01744
at γ=1, N=400, while the true ||C||₂=0.0055 (3.2× loose).  For rank-one-dominant
blocks, the Frobenius norm ||C||_F may be tighter.

METHOD
------
||C||_F^2 = Σ_{n=0}^N Σ_{m=N+1}^∞ [g(n,m)^2 / ((2n+1)(2m+1))]

where g(n,m) = (m−n)^{−γ} + (n+m+1)^{−γ}, all entries positive.

Structure: exact iv double sum for N < m ≤ M, plus a CLOSED-FORM upper bound
for the tail m > M.  For the tail:

    g(n,m) <= (m−n)^{−γ} + (n+m+1)^{−γ}
           <= (m−n)^{−γ} + m^{−γ}           [n≥0 => n+m+1>m]
           <= (m−n)^{−γ} + (m−N)^{−γ}       [N<m => m−N<m]
           <= 2(m−N)^{−γ}                   [m>M≥N+1 => m−n≥m−N]

    g(n,m)^2 <= 4(m−N)^{−2γ}

So the tail contribution from row n is

    Σ_{m>M} g(n,m)^2 / ((2n+1)(2m+1)) <= (1/(2n+1)) · 4 · Σ_{m>M} (m−N)^{−2γ}/(2m+1)

For m>M≥N+1, (m−N)^{−2γ} is decreasing in m, and 1/(2m+1) <= 1/(2M+3), so

    Σ_{m>M} (m−N)^{−2γ}/(2m+1) <= 1/(2M+3) · Σ_{k≥M+1−N} k^{−2γ}

For a decreasing sequence k^{−β}, β>1, we have

    Σ_{k≥J} k^{−β} <= J^{−β} + ∫_J^∞ x^{−β}dx = J^{−β} + J^{1−β}/(β−1)

This needs β=2γ>1, i.e. γ>1/2.

Summing over n≤N:

    tail <= [Σ_{n≤N} 1/(2n+1)] · 4/(2M+3) · T̄_{2γ}(M+1−N)

where T̄_β(J) = J^{−β} + J^{1−β}/(β−1).

Strictness: certificate (mpmath.iv, outward rounded).
"""
import sys, os, time
import numpy as np
from mpmath import iv, mp

os.environ['OMP_NUM_THREADS'] = '1'

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)


def hi(x):
    return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)


def lo(x):
    return mp.mpf(x.a) if hasattr(x, 'a') else mp.mpf(x)


def _setdps(dps):
    mp.dps = dps
    iv.dps = dps


def tail_sum_upper(J, beta):
    """UPPER bound on Σ_{k>=J} k^{-β}, β>1, integer J>=1 (iv in, iv out).

    Uses  Σ_{k>=J} k^{-β} <= J^{-β} + ∫_J^∞ x^{-β} dx = J^{-β} + J^{1-β}/(β-1).
    Direction: for decreasing k^{-β}, the left Riemann sum is >= the integral,
    so we must include the J^{-β} term explicitly. NOT just the integral alone.
    """
    J = iv.mpf(J)
    return J**(-beta) + J**(1 - beta) / (beta - 1)


def frob_C_bound_hi(N, gamma_val, M=None, dps=30):
    """CERTIFIED upper bound on ||C||_F.

    Returns (C_F_hi, finite_part, tail_part, M_used).

    If M is None, uses the default guess M = N + max(10N, 200).  NOTE: this is a
    one-shot guess, NOT an adaptive choice — the tail/finite ratio is 16% at
    γ=0.67, 1.3% at γ=1, ~0 at γ≥1.5 (N=400).  The returned bound is sound
    regardless (the tail term can only increase it); do not claim "tail ≤ 1%".

    Requires γ > 1/2 for the tail bound to be valid.
    """
    _setdps(dps)
    g = iv.mpf(gamma_val)

    if float(g.mid) <= 0.5:
        raise ValueError(f"γ={gamma_val} <= 1/2: tail bound diverges (needs 2γ>1)")

    # Pre-compute sum of 1/(2n+1) for n<=N (appears in tail bound)
    sum_inv_2n1 = iv.mpf(0)
    for n in range(N + 1):
        sum_inv_2n1 += 1 / iv.mpf(2 * n + 1)

    # Choose M if not given
    if M is None:
        # Quick estimate: tail ~ sum_inv_2n1 * 4/(2M+3) * (M+1-N)^{1-2γ}/(2γ-1)
        # Want this <= 0.01 * finite (rough finite ~ 0.01 for γ=1, N=400)
        # Start with M = N + max(10N, 200) and adjust if needed
        M_guess = N + max(10 * N, 200)
        M = M_guess

    # Compute finite part: Σ_{n=0}^N Σ_{m=N+1}^M g(n,m)^2 / ((2n+1)(2m+1))
    print(f"[frob_C_bound_hi] Computing finite sum for N={N}, M={M}, γ={gamma_val}...")
    t0 = time.time()

    finite = iv.mpf(0)
    for n in range(N + 1):
        # Vectorize inner sum over m for speed
        row_sum = iv.mpf(0)
        for m in range(N + 1, M + 1):
            n_iv = iv.mpf(n)
            m_iv = iv.mpf(m)
            # g(n,m) = (m-n)^{-γ} + (n+m+1)^{-γ}
            g_nm = (m_iv - n_iv)**(-g) + (n_iv + m_iv + 1)**(-g)
            row_sum += g_nm**2 / ((2 * m_iv + 1))
        finite += row_sum / iv.mpf(2 * n + 1)

    finite_hi = hi(finite)
    t_finite = time.time() - t0
    print(f"  Finite part computed in {t_finite:.1f}s: {float(finite_hi):.8e}")

    # Compute tail part: m > M
    # tail <= sum_inv_2n1 * 4/(2M+3) * T̄_{2γ}(M+1−N)
    t0 = time.time()
    J = M + 1 - N
    if J < 1:
        raise ValueError(f"M={M} must be >= N+1={N+1} for tail bound")

    beta = 2 * g  # beta = 2γ > 1 by the guard above
    tail_sum = tail_sum_upper(J, beta)
    tail = sum_inv_2n1 * 4 / iv.mpf(2 * M + 3) * tail_sum
    tail_hi = hi(tail)
    t_tail = time.time() - t0
    print(f"  Tail bound (J={J}) computed in {t_tail:.1f}s: {float(tail_hi):.8e}")

    C_F_sq_hi = finite_hi + tail_hi
    C_F_hi = mp.sqrt(C_F_sq_hi)

    return C_F_hi, finite_hi, tail_hi, M


def verify_against_float64(N, gamma_val, M_long):
    """Independent float64 verification (numerical-evidence only).

    Computes ||C||_2 and ||C||_F on a long truncation M_long >> N.
    """
    gam = gamma_val
    idx = np.arange(M_long + 1)
    u = 1.0 / np.sqrt(2 * idx + 1)

    head = idx[:N + 1]
    tail = idx[N + 1:]

    # Build C block
    C = np.zeros((N + 1, M_long - N))
    for i, n in enumerate(head):
        for j, m in enumerate(tail):
            g_nm = (m - n)**(-gam) + (n + m + 1.0)**(-gam)
            C[i, j] = u[n] * u[m] * g_nm

    C_2 = np.linalg.norm(C, ord=2)
    C_F = np.linalg.norm(C, ord='fro')

    return C_2, C_F


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--gamma', type=float, default=1.0)
    parser.add_argument('--N', type=int, default=400)
    parser.add_argument('--M', type=int, default=None)
    parser.add_argument('--dps', type=int, default=30)
    args = parser.parse_args()

    print("=" * 78)
    print(f"Frobenius certificate for ||C||: γ={args.gamma}, N={args.N}")
    print("=" * 78)

    t_total = time.time()
    C_F_hi, finite, tail, M_used = frob_C_bound_hi(args.N, args.gamma, M=args.M, dps=args.dps)
    t_total = time.time() - t_total

    print(f"\nRESULT (certificate):")
    print(f"  ||C||_F <= {float(C_F_hi):.8f}")
    print(f"  Finite part: {float(finite):.8e}")
    print(f"  Tail part:   {float(tail):.8e}  ({100*float(tail/finite):.2f}% of finite)")
    print(f"  M used: {M_used}")
    print(f"  Total time: {t_total:.1f}s")

    # Compare against float64 reference
    print(f"\nFloat64 reference (numerical-evidence, M_long={args.N + 40*args.N}):")
    M_long = args.N + 40 * args.N
    C_2_f64, C_F_f64 = verify_against_float64(args.N, args.gamma, M_long)
    print(f"  ||C||_2   = {C_2_f64:.8f}")
    print(f"  ||C||_F   = {C_F_f64:.8f}")
    print(f"  cert/true_F = {float(C_F_hi)/C_F_f64:.3f}x")

    if float(C_F_hi) < C_F_f64:
        print("\n*** VIOLATION: certificate < float64 true value — DIRECTION BUG ***")

    # Compare against weighted Schur bound
    from gamma_tail_bound_v2 import C_bound_weighted_hi
    C_schur, A, B = C_bound_weighted_hi(args.N, args.gamma, dps=args.dps)
    print(f"\nWeighted Schur certificate: ||C|| <= {float(C_schur):.8f}")
    print(f"  Schur/Frob = {float(C_schur)/float(C_F_hi):.3f}x")
    print(f"  Schur/true_2 = {float(C_schur)/C_2_f64:.3f}x")
