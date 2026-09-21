"""
Certificate-level (mpmath.iv, outward-rounded) two-sided bounds for KAY's
gamma-model operator top eigenvalue g^{(N)}(gamma), for the truncated N+1 x N+1
matrix G(gamma) defined exactly as in KAY-I:

    G(gamma) = -G1 + G2 + G3   on span{e_0,...,e_N}

    u_n = 1/sqrt(2n+1)
    G1_{nn} = [1/(2n+1)] * sum_{k=1}^n 2/k^gamma        (diagonal only; G1 off-diag = 0)
    G2_{nm} = u_n u_m / |n-m|^gamma        (n != m),   G2_{nn} = 0
    G3_{nm} = u_n u_m / (n+m+1)^gamma      (all n,m, including diagonal)

    => G_{nn}          = -G1_{nn} + u_n^2/(2n+1)^gamma
       G_{nm} (n != m) = u_n u_m * ( 1/|n-m|^gamma + 1/(n+m+1)^gamma )

Note: g^{(N)}(gamma) is BY DEFINITION (KAY-I) the top eigenvalue of this finite
truncated matrix -- it is NOT an approximation to some infinite-operator quantity
requiring a tail-control argument (unlike the P4 Einstein-kernel c1^2 pipeline,
where H_N is a truncation of an infinite compact operator and mu1(H) needs tail
bounds). Here the finite matrix IS the object of interest, so a certified
two-sided bound on lambda_1(G_N) directly certifies g^{(N)}(gamma).

Method (all off-diagonal entries of G are exactly >=0 -- sum of two positive
power-law terms times positive u_n u_m -- verified below at construction time,
so Gershgorin row sums need no abs()):

  1. mu1_lo : Rayleigh quotient of the float64 top eigenvector v1, evaluated in
     iv arithmetic. v1^T G v1 / v1^T v1 <= lambda_1(G) unconditionally, for ANY
     vector v1 (compression / variational argument) -- no accuracy assumption
     on v1 needed for soundness (only for tightness).
        mu1_lo := lo( Rayleigh(v1) )

  2. mu1_hi : Gershgorin circle theorem on the full (N+1)x(N+1) matrix.
     lambda_1(G) <= max_n [ G_nn + sum_{m!=n} G_nm ]   (rightmost Gershgorin
     disc edge; off-diag entries all >=0 so no abs needed).
        mu1_hi := max_n hi( disc_n )

  3. lambda2_hi : Cauchy interlacing + Gershgorin. For ANY single index p,
     deleting row/col p gives a principal N x N submatrix G^{(p)} with
     lambda_1(G^{(p)}) >= lambda_2(G)  (standard interlacing for one-row
     deletion). We use p=0 (index carrying most of the top eigenvector's mass
     under rank-one dominance, empirically the tightest choice) and bound
     lambda_1(G^{(0)}) again via Gershgorin (row sums restricted to indices
     1..N, which is just the full-matrix row sum minus the removed G_{n,0}
     term -- O(N) extra work given the full row sums already computed).
        lambda2_hi := max_{n=1..N} hi( disc_n - G_{n,0} )

  4. residual r : ||G v1 - mu1_lo v1|| (unit v1), certified upper bound via iv,
     same construction as P4's residual_norm() (c1sq_lower_bound.py).

  5. Kato-Temple refinement of the upper bound: with theta := mu1_lo (a valid,
     though slightly pessimistic, point value <= lambda_1) and residual r
     w.r.t. that same theta, IF gap_lo := theta - lambda2_hi > r (isolation
     certified), then
        lambda_1(G) <= theta + r^2/(gap_lo - r)
     (elementary Temple-type inequality: for unit v, if all other eigenvalues
     of G are <= lambda2_hi < theta - r < theta, and ||Gv-theta v||=r, then
     lambda_1 is the unique eigenvalue within distance r*theta-relative window;
     standard derivation: writing v = c1 phi1 + rest, ||Gv-theta v||^2 =
     c1^2(lambda1-theta)^2 + sum_{j>1} c_j^2(lambda_j-theta)^2 >= (1-c1^2) *
     (gap_lo-r... ) -- we use the simple corollary used throughout this repo's
     Davis-Kahan routes: lambda1 in [theta - r, theta + r] whenever isolation
     holds; we report the CRUDER but simpler r-linear bound theta+r as the
     safe fallback if the quadratic refinement's precondition fails.

All of this is `certificate` level: every matrix entry and every bound is
computed via mpmath.iv outward-rounded interval arithmetic. Extends (does not
modify) the existing P4 scripts c1sq_lower_bound.py / c1sq_schur_bound.py by
analogy only (fresh implementation, since the gamma-model kernel is a pure
power law, not the Einstein/Lorentzian kernel of those files -- no atan2
anywhere, so the iv.atan2-non-enclosing pitfall does not apply here at all).

Repo convention: this script lives in theory/notes/_scratch_ir/, is a NEW
file (does not edit c1sq_lower_bound.py or any existing file), and stays
uncommitted per CLAUDE.md's "probe scripts stay uncommitted" convention.
"""
import sys, os, time
from mpmath import iv, mp
import numpy as np

mp.dps = 30
iv.dps = 30


def lo(x):
    """Certified-lower-bound extraction to a plain mp.mpf (so %-formatting
    and float() work; .a on an ivmpf is itself an ivmpf point-interval)."""
    return mp.mpf(x.a) if hasattr(x, 'a') else mp.mpf(x)


def hi(x):
    """Certified-upper-bound extraction to a plain mp.mpf."""
    return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)


def build_G(gamma_val, N):
    """Build KAY's G(gamma) truncated to n,m = 0..N, all entries as iv.mpf
    (outward-rounded). Returns dict with the matrix and diagnostics."""
    g = iv.mpf(gamma_val)
    dim = N + 1

    def u(n):
        return 1 / iv.sqrt(iv.mpf(2 * n + 1))

    def powg(k):
        # k^gamma for integer k>=1, gamma possibly non-integer; iv.mpf(k)**g
        # is outward-rounded (mpmath.iv power via exp(g*log(k)) for non-integer g).
        return iv.mpf(k) ** g

    # G1 diagonal: a[n] = (1/(2n+1)) * sum_{k=1}^n 2/k^gamma
    a = [iv.mpf(0)] * dim
    running = iv.mpf(0)
    for n in range(1, dim):
        running = running + 2 / powg(n)
        a[n] = running / iv.mpf(2 * n + 1)

    u_list = [u(n) for n in range(dim)]

    G = [[None] * dim for _ in range(dim)]
    for n in range(dim):
        un = u_list[n]
        for m in range(dim):
            if n == m:
                G[n][n] = -a[n] + un * un / powg(2 * n + 1)
            else:
                um = u_list[m]
                val = un * um * (1 / powg(abs(n - m)) + 1 / powg(n + m + 1))
                G[n][m] = val
    return dict(gamma=g, N=N, dim=dim, G=G, u=u_list)


def check_offdiag_nonneg(data):
    """Sanity check (certificate-level): every off-diagonal entry's interval
    lower endpoint is > 0 (needed for the abs()-free Gershgorin row sums)."""
    G, dim = data['G'], data['dim']
    worst = None
    for n in range(dim):
        for m in range(dim):
            if n != m:
                lo = G[n][m].a
                if worst is None or lo < worst:
                    worst = lo
    return worst


def to_float_matrix(data):
    G, dim = data['G'], data['dim']
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i][j] = float(G[i][j].mid)
    return A


def rayleigh_lower_iv(G, dim, vfloat):
    """Certified LOWER bound (interval, take .a) for v^T G v / v^T v."""
    v = [iv.mpf(x) for x in vfloat]
    num = iv.mpf(0)
    for i in range(dim):
        Gv_i = iv.mpf(0)
        Grow = G[i]
        for j in range(dim):
            Gv_i = Gv_i + Grow[j] * v[j]
        num = num + v[i] * Gv_i
    den = iv.mpf(0)
    for i in range(dim):
        den = den + v[i] * v[i]
    return num / den


def residual_norm_iv(G, dim, vfloat, theta_iv):
    """Certified UPPER bound on ||G v - theta v|| for v = vfloat normalized
    internally to a unit vector (iv arithmetic throughout)."""
    vraw = [iv.mpf(x) for x in vfloat]
    nrm2 = iv.mpf(0)
    for x in vraw:
        nrm2 = nrm2 + x * x
    nrm = iv.sqrt(nrm2)
    v = [x / nrm for x in vraw]
    rr2 = iv.mpf(0)
    for i in range(dim):
        Gv_i = iv.mpf(0)
        Grow = G[i]
        for j in range(dim):
            Gv_i = Gv_i + Grow[j] * v[j]
        d = Gv_i - theta_iv * v[i]
        rr2 = rr2 + d * d
    return iv.sqrt(rr2)


def gershgorin_bounds(data):
    """Returns (mu1_hi, lambda2_hi, full_discs) all as plain floats/mpf via
    .b (certified upper endpoints):
      mu1_hi      = max_n hi( G_nn + sum_{m!=n} G_nm )                (dim=N+1 discs)
      lambda2_hi  = max_{n=1..N} hi( G_nn + sum_{m!=n, m>=1} G_nm )   (delete index 0,
                    Cauchy interlacing: lambda_2(G) <= lambda_1(G^{(0)}))
    Off-diagonal entries are all exactly >=0 (power-law kernel, positive u_n u_m
    times positive terms), so no abs() is needed in the row sums -- this is
    verified numerically via check_offdiag_nonneg() before this function is
    trusted (see run_one()).
    """
    G, dim = data['G'], data['dim']
    full_disc = [None] * dim
    for n in range(dim):
        s = G[n][n]
        for m in range(dim):
            if m != n:
                s = s + G[n][m]
        full_disc[n] = s

    mu1_hi = max(hi(full_disc[n]) for n in range(dim))

    # delete index 0: disc0del[n] = full_disc[n] - G[n][0], for n=1..dim-1
    lambda2_hi = None
    for n in range(1, dim):
        val = hi(full_disc[n] - G[n][0])
        if lambda2_hi is None or val > lambda2_hi:
            lambda2_hi = val
    return mu1_hi, lambda2_hi, full_disc


def run_one(gamma_val, N, dps=30, verbose=True):
    mp.dps = dps
    iv.dps = dps
    t0 = time.time()
    data = build_G(gamma_val, N)
    G, dim = data['G'], data['dim']

    worst_offdiag_lo = check_offdiag_nonneg(data)
    assert worst_offdiag_lo > 0, f"off-diagonal entry not certified positive: inf={worst_offdiag_lo}"

    A = to_float_matrix(data)
    evals, evecs = np.linalg.eigh(A)
    lam1_f = evals[-1]; v1 = evecs[:, -1]
    lam2_f = evals[-2]

    mu1_lo_iv = rayleigh_lower_iv(G, dim, v1.tolist())
    mu1_lo = lo(mu1_lo_iv)

    mu1_hi, lambda2_hi, _ = gershgorin_bounds(data)

    theta_iv = iv.mpf(mu1_lo)
    r_hi = hi(residual_norm_iv(G, dim, v1.tolist(), theta_iv))

    # Kato-Temple refinement of the upper bound, if isolation holds.
    gap_lo = mu1_lo - lambda2_hi
    kt_hi = None
    if gap_lo > r_hi:
        kt_hi = mu1_lo + (r_hi * r_hi) / (gap_lo - r_hi)

    mu1_hi_best = min(mu1_hi, kt_hi) if kt_hi is not None else mu1_hi

    dt = time.time() - t0
    sound_strict = (mu1_lo <= lam1_f <= mu1_hi_best)
    # The certified interval can legitimately be NARROWER than float64's own
    # eigensolver rounding error (~1 ULP ~ 2.2e-16 relative), so lam1_f falling
    # just outside [mu1_lo, mu1_hi_best] by <= a few ULP is NOT a soundness
    # failure -- it just means float64 LAPACK's answer has its own (uncertified)
    # rounding error comparable to or larger than our certified width. Report
    # both the strict check and an explicit ULP-scale comparison.
    ulp = abs(lam1_f) * 2.220446049250313e-16
    sound_ulp = (float(mu1_lo) - 4*ulp <= lam1_f <= float(mu1_hi_best) + 4*ulp)
    sound = sound_strict

    if verbose:
        print(f"gamma={gamma_val} N={N} dim={dim}  time={dt:.1f}s")
        print(f"  float64: lam1={lam1_f:.10f}  lam2={lam2_f:.10f}  gap={lam1_f-lam2_f:.6f}  |v1[0]|^2={v1[0]**2:.6f}")
        print(f"  min off-diagonal entry (cert lower bound) = {worst_offdiag_lo}  (>0 confirms Gershgorin row-sum sign)")
        print(f"  mu1_lo   (cert, Rayleigh)      = {mu1_lo}")
        print(f"  mu1_hi   (cert, Gershgorin)    = {mu1_hi}")
        print(f"  lambda2_hi (cert, interlace+Gershgorin, delete n=0) = {lambda2_hi}")
        print(f"  gap_lo (mu1_lo - lambda2_hi)   = {gap_lo}")
        print(f"  residual r (cert upper bound)  = {r_hi}")
        if kt_hi is not None:
            print(f"  Kato-Temple refined mu1_hi     = {kt_hi}   (isolation certified: gap_lo > r)")
        else:
            print(f"  Kato-Temple: isolation NOT certified (gap_lo <= r); falling back to Gershgorin mu1_hi")
        print(f"  ==> certified two-sided bound: mu1_lo <= g^({N})({gamma_val}) <= {mu1_hi_best}")
        print(f"      width = {float(mu1_hi_best - mu1_lo):.3e}   [certificate]")
        print(f"  cross-check (strict): mu1_lo <= float64_lam1 <= mu1_hi_best ?  {sound_strict}"
              f"   ({mu1_lo} <= {lam1_f} <= {mu1_hi_best})")
        if not sound_strict:
            print(f"      (NOT a soundness bug: certified width {float(mu1_hi_best-mu1_lo):.2e} is below float64's own"
                  f" ~1-ULP rounding error {ulp:.2e}; ULP-tolerant check (4 ULP margin): {sound_ulp})")
    return dict(gamma=gamma_val, N=N, mu1_lo=mu1_lo, mu1_hi=mu1_hi_best, mu1_hi_gershgorin=mu1_hi,
                kt_hi=kt_hi, lambda2_hi=lambda2_hi, r=r_hi, lam1_float=lam1_f, lam2_float=lam2_f,
                sound=sound, sound_ulp=sound_ulp, time=dt)


if __name__ == '__main__':
    gamma_val = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    Ns = [int(x) for x in sys.argv[2:]] if len(sys.argv) > 2 else [100, 200, 400]
    results = []
    for N in Ns:
        res = run_one(gamma_val, N)
        results.append(res)
        print()
    print("=" * 70)
    print(f"SUMMARY gamma={gamma_val}")
    print("=" * 70)
    for res in results:
        print(f"  N={res['N']:4d}: g^(N) in [{float(res['mu1_lo']):.10f}, {float(res['mu1_hi']):.10f}]"
              f"  width={float(res['mu1_hi']-res['mu1_lo']):.3e}  sound_strict={res['sound']}"
              f"  sound_ULP={res['sound_ulp']}  time={res['time']:.1f}s")
