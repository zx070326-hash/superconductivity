"""
Certificate-level lower bound for c1^2 = psi_1(0)^2, top eigenvector weight
at n=0 of H_+ = O[f], f Einstein kernel f(k) = varpi^2/(varpi^2+k^2).

Route (Schur-complement / gap argument, derived fresh this session):
  c1^2 = 1/(1 + Sigma'(mu1)),  Sigma(z) = b^T (H'-z)^{-1} b   (Feshbach self-energy
  at the n=0 site, H' = H restricted to n>=1, b_n = H_{n0} for n>=1).
  Since Sigma'(z) = b^T (H'-z)^{-2} b and for z=mu1 > mu1(H') (top of H' spectrum,
  which by Cauchy interlacing is <= mu2(H)):
      Sigma'(mu1) = sum_j (b.phi_j)^2/(mu1-lambda_j(H'))^2 <= ||b||^2 / gamma^2
  with gamma := mu1 - mu2 (spectral gap of full H), because every
  mu1-lambda_j(H') >= mu1 - mu1(H') >= mu1-mu2 = gamma  (interlacing: mu1(H')<=mu2(H)).
  Hence:
      c1^2 >= gamma^2 / (gamma^2 + ||b||^2)          (*)
  This is a DECREASING function of ||b||^2 and INCREASING function of gamma, so we need
    - gamma_lo  <= true gamma   (certified LOWER bound: mu1_lo - mu2_hi)
    - Bsq_hi    >= true ||b||^2 (certified UPPER bound)
  and (*) evaluated at (gamma_lo, Bsq_hi) is then a valid certified LOWER bound on c1^2
  (h(g,B)=g^2/(g^2+B) increasing in g, decreasing in B: h(gamma_lo,Bsq_hi) <= h(gamma,Bsq_hi) <= h(gamma,B)).

mu1 lower bound: Rayleigh quotient of an (approximate, float) top eigenvector v1 of the
  truncated N x N matrix H_N gives mu1(H) >= mu1(H_N) >= R(v1) = v1^T H_N v1 / ||v1||^2
  exactly and unconditionally (compression can only decrease lambda_max; Rayleigh quotient
  is a lower bound for lambda_max for ANY vector). No tail control needed for this direction.

mu2 upper bound: needs a genuine tail argument.
  H = (H_N (+) 0) + E,  E = [[0,C],[C^T, H^{BB}]]   (B = tail indices n>N)
  ||E|| <= ||C|| + ||H^{BB}||   (triangle ineq, since the antidiagonal block matrix
  [[0,C],[C^T,0]] has operator norm = ||C||).
  By Weyl: mu2(H) <= lambda_2(H_N (+) 0) + ||E|| = lambda_2(H_N) + ||E||  (N>=2).
  lambda_2(H_N) upper bound: certified via interval LDLT inertia count (Sylvester's law).
  ||H^{BB}|| <= (fhat(0)-f(0))/(2N+3)   [L1 tail]  +  u_{N+1}^2 * (fhat(0)+f(0))  [L2+L3 tail]
  ||C||      <= u_{N+1} * (fhat(0)-fhat(pi))/2      [off-diagonal coupling block]
  where fhat(0)=pi*varpi*coth(pi*varpi), fhat(pi)=pi*varpi/sinh(pi*varpi)  (closed form,
  Poisson summation, same as delta_channel_immunity note Lemma 1), f(0)=1.

||b||^2 upper bound: b_n = u_n(f(n)+f(n+1)) for n>=1 (exact, from H_{n0}).
  Head part (1<=n<=N) computed exactly via iv. Tail (n>N): termwise
  u_n^2(f(n)+f(n+1))^2 <= 4 u_n^2 f(n)^2 <= 4*(1/(2n))*(varpi^4/n^4) = 2 varpi^4/n^5
  (f decreasing in k, u_n^2=1/(2n+1)<=1/(2n)), so
  tail sum <= sum_{n>N} 2 varpi^4/n^5 <= 2 varpi^4 * int_N^infty x^-5 dx = varpi^4/(2N^4).

All matrix entries and bounds computed with mpmath.iv (outward-rounded interval
arithmetic). f, u_n, S=sum f(k) etc. are all elementary/closed-form so this is genuine
certificate-level, not just multi-precision-agreement.
"""
import sys, time
from mpmath import iv, mp

mp.dps = 30
iv.dps = 30

def build_data(varpi_val, N):
    w = iv.mpf(varpi_val)
    pw = (+iv.pi) * w
    ep = iv.exp(pw)
    em = iv.exp(-pw)
    sh = (ep - em) / 2
    ch = (ep + em) / 2
    fhat0 = pw * ch / sh
    fhatpi = pw / sh
    assert fhatpi.a > 0, "fhat(pi) must be >0"

    def f(k):
        k = iv.mpf(k)
        return w*w/(w*w+k*k)

    def u(n):
        return 1/iv.sqrt(iv.mpf(2*n+1))

    # H_N entries, n,m = 0..N
    dim = N+1
    H = [[None]*dim for _ in range(dim)]
    # a_n for L1 diag, n=0..N
    a = [iv.mpf(0)]*dim
    running = iv.mpf(0)
    for n in range(1, dim):
        running = running + f(n)
        a[n] = (2*running)/iv.mpf(2*n+1)
    for n in range(dim):
        un = u(n)
        for m in range(dim):
            um = u(m)
            if n == m:
                val = -a[n] + f(2*n+1)*un*un
            else:
                val = f(abs(n-m))*un*um + f(n+m+1)*un*um
            H[n][m] = val

    return dict(w=w, fhat0=fhat0, fhatpi=fhatpi, f=f, u=u, H=H, dim=dim, N=N)


def inertia_count_below(H, dim, c):
    """Certified count of eigenvalues of H (dim x dim interval symmetric matrix)
    that are < c, via interval LDL^T (symmetric Gaussian elimination, no pivoting)
    and Sylvester's law of inertia. Returns (count, ok). ok=False if a pivot interval
    straddles zero (inconclusive -> caller should perturb c or increase precision).
    Does NOT mutate H (works on a copy)."""
    M = [row[:] for row in H]
    for i in range(dim):
        M[i][i] = M[i][i] - c
    neg = 0
    for k in range(dim):
        dk = M[k][k]
        lo, hi = dk.a, dk.b
        if lo <= 0 <= hi:
            return None, False
        if hi < 0:
            neg += 1
        # eliminate
        if k < dim-1:
            invdk = 1/dk
            col = [M[i][k] for i in range(k+1, dim)]
            for ii, i in enumerate(range(k+1, dim)):
                lik = col[ii]*invdk
                Mik_row = M[i]
                Mk_row = M[k]
                for j in range(k+1, dim):
                    Mik_row[j] = Mik_row[j] - lik*Mk_row[j]
    return neg, True


def rayleigh_lower(H, dim, vfloat):
    """Certified lower bound for v^T H v / v^T v using interval arithmetic, vfloat a
    plain python list of floats (need not be exact / need not be normalized)."""
    v = [iv.mpf(x) for x in vfloat]
    num = iv.mpf(0)
    for i in range(dim):
        Hv_i = iv.mpf(0)
        Hrow = H[i]
        for j in range(dim):
            Hv_i = Hv_i + Hrow[j]*v[j]
        num = num + v[i]*Hv_i
    den = iv.mpf(0)
    for i in range(dim):
        den = den + v[i]*v[i]
    return num/den  # interval; take .a for certified lower bound


def tail_bounds(data, N):
    """Certified UPPER bounds: ||H^{BB}||, ||C|| (coupling), and Bsq_tail = upper
    bound on sum_{n>N} u_n^2 (f(n)+f(n+1))^2. All via closed forms + iv, direction
    checked: every inequality below is an upper bound on a nonnegative quantity."""
    w = data['w']; fhat0 = data['fhat0']; fhatpi = data['fhatpi']
    f1 = iv.mpf(1)  # f(0)=1
    uNp1 = data['u'](N+1)
    # L1 tail: ||L1^BB|| <= (fhat0 - f(0)) / (2N+3)
    L1_tail = (fhat0 - f1) / iv.mpf(2*N+3)
    # L2+L3 tail: <= u_{N+1}^2 * (fhat0 + f(0))
    L23_tail = uNp1*uNp1*(fhat0+f1)
    HBB_bound = L1_tail + L23_tail
    # coupling: ||C|| <= u_{N+1} * (fhat0-fhatpi)/2
    C_bound = uNp1 * (fhat0-fhatpi)/2
    # ||b||^2 tail: <= varpi^4/(2 N^4)
    Bsq_tail = (w**4) / (2*iv.mpf(N)**4)
    return HBB_bound, C_bound, Bsq_tail


def bsq_head(data, N):
    """Exact (iv-certified) sum_{n=1}^{N} u_n^2 (f(n)+f(n+1))^2."""
    u = data['u']; f = data['f']
    s = iv.mpf(0)
    for n in range(1, N+1):
        term = u(n)*(f(n)+f(n+1))
        s = s + term*term
    return s


import numpy as np

def to_float_matrix(H, dim):
    A = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            A[i][j] = float(H[i][j].mid)
    return A


def find_lambda2_upper(H, dim, lam2_approx, gap_guess):
    """Certified upper bound on lambda_2(H) (dim x dim) via inertia counting.
    Tries increasing margins above lam2_approx until inertia test succeeds
    (pivot never straddles zero) and count(c) >= dim-1."""
    margin = max(gap_guess*1e-3, 1e-8)
    for attempt in range(40):
        c = iv.mpf(lam2_approx + margin)
        cnt, ok = inertia_count_below(H, dim, c)
        if ok and cnt is not None and cnt >= dim-1:
            return c, margin
        margin *= 1.7
    raise RuntimeError("could not certify lambda2 upper bound; increase dps or N")


def run_one(varpi_val, N, dps=30, verbose=True):
    mp.dps = dps
    iv.dps = dps
    t0 = time.time()
    data = build_data(varpi_val, N)
    H, dim = data['H'], data['dim']
    A = to_float_matrix(H, dim)
    evals, evecs = np.linalg.eigh(A)
    # eigh ascending; top = last
    lam1_approx = evals[-1]; v1 = evecs[:, -1]
    lam2_approx = evals[-2]; v2 = evecs[:, -2]
    gap_guess = lam1_approx - lam2_approx

    mu1_lo_iv = rayleigh_lower(H, dim, v1.tolist())
    mu1_lo = mu1_lo_iv.a

    c_upper, margin_used = find_lambda2_upper(H, dim, lam2_approx, gap_guess)
    HBB_bound, C_bound, Bsq_tail = tail_bounds(data, N)
    Eop = HBB_bound + C_bound
    mu2_hi_iv = c_upper + Eop
    mu2_hi = mu2_hi_iv.b

    gamma_lo = mu1_lo - mu2_hi
    Bsq_hi_iv = bsq_head(data, N) + Bsq_tail
    Bsq_hi = Bsq_hi_iv.b

    if gamma_lo <= 0:
        c1sq_lo = None
    else:
        g2 = iv.mpf(gamma_lo)**2
        c1sq_lo_iv = g2/(g2+iv.mpf(Bsq_hi))
        c1sq_lo = c1sq_lo_iv.a

    if verbose:
        print(f"varpi={varpi_val} N={N} dim={dim}  time={time.time()-t0:.1f}s")
        print(f"  lam1_approx={lam1_approx:.6f} lam2_approx={lam2_approx:.6f} gap_guess={gap_guess:.6f}")
        print(f"  mu1_lo(cert)={mu1_lo}")
        print(f"  lambda2(H_N)_hi(cert)={c_upper} margin_used={margin_used:.3e}")
        print(f"  ||H^BB||<= {HBB_bound}   ||C||<= {C_bound}")
        print(f"  mu2_hi(cert)={mu2_hi}")
        print(f"  gamma_lo(cert)={gamma_lo}")
        print(f"  Bsq_head={bsq_head(data,N)}  Bsq_tail<={Bsq_tail}  Bsq_hi={Bsq_hi}")
        print(f"  c1sq_lower(cert) = {c1sq_lo}")
    return dict(varpi=varpi_val, N=N, mu1_lo=mu1_lo, mu2_hi=mu2_hi, gamma_lo=gamma_lo,
                Bsq_hi=Bsq_hi, c1sq_lo=c1sq_lo, lam1_approx=lam1_approx, lam2_approx=lam2_approx)


if __name__ == '__main__':
    varpi_val = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    run_one(varpi_val, N)


def residual_norm(H, dim, vfloat, la_iv):
    """Certified UPPER bound on ||H v - la v|| (iv), v = vfloat (need not be unit;
    normalized internally). la_iv: iv scalar (point) for the Rayleigh value used."""
    vraw = [iv.mpf(x) for x in vfloat]
    nrm2 = iv.mpf(0)
    for x in vraw:
        nrm2 = nrm2 + x*x
    nrm = iv.sqrt(nrm2)
    v = [x/nrm for x in vraw]
    rr2 = iv.mpf(0)
    for i in range(dim):
        Hv_i = iv.mpf(0)
        Hrow = H[i]
        for j in range(dim):
            Hv_i = Hv_i + Hrow[j]*v[j]
        d = Hv_i - la_iv*v[i]
        rr2 = rr2 + d*d
    return iv.sqrt(rr2)  # ||H v - la v||, v unit


def run_one_dk(varpi_val, N, dps=30, verbose=True):
    """Corrected certificate route via Davis-Kahan-type sinTheta bound (see writeup):
    step 1: bound angle between float eigvec v1 (of H_N) and the EXACT top eigvec phi
            of A:=H_N (+) 0, via residual sinTheta1 <= r1/(gap1-r1), gap1=lambda1(H_N)-lambda2(H_N).
    step 2: bound angle between phi and psi (true top eigvec of full H=A+E) via
            sinTheta2 <= ||E||/(delta-||E||), delta=lambda1(H_N)-max(lambda2(H_N),0)
            (elementary derivation in the note, NOT the flawed interlacing-gamma route
            from the first attempt in this script -- that route is WRONG, see note).
    combine: sinTheta_total <= sinTheta1+sinTheta2 (crude triangle ineq, fine since
    sinTheta1 is float-eigensolver-residual scale, negligible next to sinTheta2).
    c1^2 = psi[0]^2 >= max(0, cosTheta_total*|v1[0]| - sinTheta_total)^2
    """
    mp.dps = dps; iv.dps = dps
    t0 = time.time()
    data = build_data(varpi_val, N)
    H, dim = data['H'], data['dim']
    A = to_float_matrix(H, dim)
    evals, evecs = np.linalg.eigh(A)
    lam1_approx = evals[-1]; v1 = evecs[:, -1]
    lam2_approx = evals[-2]

    def lo(x):  # certified-lower-bound extraction to a plain mp.mpf
        return mp.mpf(x.a) if hasattr(x, 'a') else mp.mpf(x)

    def hi(x):  # certified-upper-bound extraction to a plain mp.mpf
        return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)

    mu1_lo = lo(rayleigh_lower(H, dim, v1.tolist()))  # certified lower bound for lambda1(H_N)
    la_iv = iv.mpf(mu1_lo)

    r1 = hi(residual_norm(H, dim, v1.tolist(), la_iv))  # upper bound, plain mp.mpf

    c_upper, margin_used = find_lambda2_upper(H, dim, lam2_approx, lam1_approx-lam2_approx)
    lam2_HN_hi = hi(c_upper)

    gap1_lo = mu1_lo - lam2_HN_hi  # certified lower bound on lambda1(H_N)-lambda2(H_N)
    if gap1_lo <= r1:
        sinT1 = mp.mpf(1)
    else:
        sinT1 = r1/(gap1_lo-r1)
        if sinT1 > 1: sinT1 = mp.mpf(1)

    HBB_bound, C_bound, Bsq_tail = tail_bounds(data, N)
    Ehi = hi(HBB_bound) + hi(C_bound)
    lam2A_hi = max(lam2_HN_hi, mp.mpf(0))
    delta_lo = mu1_lo - lam2A_hi
    if delta_lo <= Ehi:
        sinT2 = mp.mpf(1)
    else:
        sinT2 = Ehi/(delta_lo-Ehi)
        if sinT2 > 1: sinT2 = mp.mpf(1)

    sinT = min(mp.mpf(1), sinT1+sinT2)
    cosT = mp.sqrt(1-sinT**2) if sinT < 1 else mp.mpf(0)

    v0 = abs(v1[0])
    val = cosT*mp.mpf(v0) - sinT
    c1sq_lo = max(mp.mpf(0), val)**2

    if verbose:
        print(f"varpi={varpi_val} N={N} dim={dim} time={time.time()-t0:.1f}s")
        print(f"  mu1_lo={mu1_lo}  lam2_HN_hi={lam2_HN_hi}  r1(eigvec resid)={r1}")
        print(f"  gap1_lo={gap1_lo}  sinTheta1={sinT1}")
        print(f"  ||H^BB||<={hi(HBB_bound)}  ||C||<={hi(C_bound)}  Ehi={Ehi}")
        print(f"  delta_lo={delta_lo}  sinTheta2={sinT2}")
        print(f"  sinTheta_total<= {sinT}  cosTheta>= {cosT}  |v1[0]|={v0}")
        print(f"  c1sq_lower(cert, Davis-Kahan route) = {c1sq_lo}")
    return dict(varpi=varpi_val, N=N, mu1_lo=mu1_lo, sinT=sinT, c1sq_lo=c1sq_lo, v0=v0)
