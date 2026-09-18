"""
Improved ||C|| (head-tail coupling block) upper bound via an UNWEIGHTED Schur test,
replacing the global-operator-norm-sandwich bound in c1sq_lower_bound.py's tail_bounds().

Does NOT modify c1sq_lower_bound.py (imported as a module, read-only). New file, per
convention (c1sq_ prefix, theory/notes/_scratch_ir/, uncommitted probe script).

Motivation (from team-lead's lead + this session's own float64 prototyping,
_scratch_ir/c1sq_schur_check.py, c1sq_schur_closed.py, c1sq_schur_analytic_proto*.py):
the old bound ||C|| <= u_{N+1}*(fhat(0)-fhat(pi))/2 sandwiches the WHOLE coupling
block against the global operator norm of T_f+H_f, discarding the pointwise decay
f(k)~varpi^2/k^2. A naive (unweighted) Schur test,

    C_{nm} = u_n u_m ( f(|n-m|) + f(n+m+1) ),   n<=N < m
    ||C|| <= sqrt( sup_n R_n * sup_m K_m ),   R_n = sum_m |C_nm|, K_m = sum_n |C_nm|

already recovers the correct N^{-1} decay rate (float64 measured; see
c1sq_Cnorm_rate.py from the previous round) instead of the old bound's N^{-1/2}.

Derivation, direction-audited term by term (all substitutions are UPPER bounds on
nonnegative quantities -- see NOTE lines):

1. Row sum R_n (n<=N, exact sum over infinite tail m>N -- needs real tail control):
     R_n = u_n * sum_{m>N} u_m (f(m-n)+f(n+m+1))
   NOTE(upper): u_m <= u_{N+1} for all m>N+1... m>=N+1 (u_n=(2n+1)^{-1/2} strictly
   decreasing in n, so m>=N+1 => u_m<=u_{N+1}). Substituting:
     R_n <= u_n u_{N+1} * sum_{m>N} (f(m-n)+f(n+m+1))
          = u_n u_{N+1} * ( Tail(N+1-n) + Tail(n+N+2) )
   where Tail(J) := sum_{j>=J} f(j) (J>=1 integer), bounded analytically below.
   sup_n R_n is NOT always at n=N (float64 check: for varpi/N large, e.g.
   varpi=30,100 at N=500, argmax is at n=0, not n=N -- team-lead's crude
   "worst at n~N" estimate does not hold uniformly). We therefore certify
   sup_n R_n by an O(N) loop over ALL n=0..N (finite, exact max of a finite
   set of individually-valid upper bounds -- trivially a valid upper bound
   for the true sup).

2. Column sup K_m (m>N): here the INNER sum over n=0..N is FINITE, so no tail
   control is needed for it -- it is computed EXACTLY via iv (not bounded).
     K_m = u_m * sum_{n=0}^N u_n (f(m-n)+f(n+m+1))
   Monotonicity in m (proved, not assumed): for FIXED n, u_m, f(m-n), f(n+m+1)
   are each non-increasing in m (u decreasing; f decreasing in its argument and
   both m-n, n+m+1 increase with m). A sum of non-increasing nonneg functions of
   m is non-increasing in m. Hence K_m is EXACTLY non-increasing in m>N, so
     sup_{m>N} K_m = K_{N+1}   (exact, not a bound -- proved via elementary monotonicity)
   computed as a single closed finite ((N+1)-term) iv sum -- certificate level,
   no analytic tail needed here at all.

3. Tail(J) := sum_{j>=J} f(j), f(j)=varpi^2/(varpi^2+j^2), J>=1 integer.
   Two INDEPENDENT closed-form upper bounds, take the min (min of two valid
   upper bounds is again a valid, tighter, upper bound):

   (a) Integral-comparison + arctan-avoidance (tight for J small vs varpi):
       NOTE(upper): f decreasing on x>=0 => f(j) <= integral_{j-1}^{j} f(x)dx
       (f(x)>=f(j) on [j-1,j]). Telescoping j=J..infty:
         Tail(J) <= integral_{J-1}^{infty} f(x)dx = varpi*(pi/2 - arctan((J-1)/varpi))
       NOTE(upper): arctan(y) >= y/(1+y^2) for y>=0 (standard: equality at y=0,
       d/dy[arctan(y)-y/(1+y^2)] = 2y^2/(1+y^2)^2 >= 0, so the difference is
       non-decreasing from 0). Hence pi/2-arctan(y) <= pi/2 - y/(1+y^2), giving
         Tail(J) <= varpi*( pi/2 - y/(1+y^2) ),   y=(J-1)/varpi
       -- purely rational in y plus the single constant pi/2 (mp.pi triggered via
       "+iv.pi" per repo convention); NO iv.atan / iv.atan2 call anywhere (known
       non-enclosing in this mpmath build -- repo pitfall list).

   (b) Termwise 1/j^2 comparison (tight for J large vs varpi):
       NOTE(upper): f(j) = varpi^2/(varpi^2+j^2) <= varpi^2/j^2 for j>=1 (denominator
       only larger on the LHS). NOTE(upper): sum_{j>=J} 1/j^2 <= 1/(J-1) for J>=2
       (telescoping 1/j^2 <= 1/(j(j-1)) = 1/(j-1)-1/j for j>=2). For J=1:
       sum_{j>=1}1/j^2 <= 1 + sum_{j>=2}1/j^2 <= 1+1 = 2 (crude but valid, purely
       rational, avoids needing pi^2/6). Hence
         Tail(J) <= varpi^2 * ( 2  if J<=1  else  1/(J-1) )

   Tail(J) := min( bound(a), bound(b) )  -- both individually valid upper bounds
   on the SAME nonnegative quantity sum_{j>=J} f(j), so the min is a valid
   (and strictly tighter) upper bound. Direction-audited: both (a) and (b) are
   monotone-substitution upper bounds throughout, no lower-bound step anywhere.

Float64 cross-validation (this session, c1sq_schur_analytic_proto2.py): this exact
bound (windowless R + exact K) vs true ||C|| (numpy SVD of the finite-tail-truncated
block, tail truncated at M=N+40000, converged) gives ratio 1.9x-5.6x across
varpi in [0.5,100], N in {500,1000} -- ALWAYS >= true value (no sign of the
"bound < true" bug pattern from the previous round's interlacing mistake), and a
factor 15-30x tighter than the old global-norm bound (which was 55-127x loose).
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import iv, mp
import numpy as np

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "c1sq_lower_bound", os.path.join(os.path.dirname(os.path.abspath(__file__)), "c1sq_lower_bound.py"))
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)

build_data = _base.build_data
rayleigh_lower = _base.rayleigh_lower
residual_norm = _base.residual_norm
find_lambda2_upper = _base.find_lambda2_upper
tail_bounds = _base.tail_bounds
to_float_matrix = _base.to_float_matrix


def lo(x):
    return mp.mpf(x.a) if hasattr(x, 'a') else mp.mpf(x)


def hi(x):
    return mp.mpf(x.b) if hasattr(x, 'b') else mp.mpf(x)


def tail_bound_hi(J, w_iv):
    """Certified UPPER bound (plain mp.mpf) for sum_{j>=J} f(j), J>=1 integer.
    min of two independent closed-form upper bounds (see module docstring §3)."""
    y = iv.mpf(J - 1) / w_iv
    boundA = w_iv * ((+iv.pi) / 2 - y / (1 + y * y))
    if J <= 1:
        boundB = 2 * w_iv * w_iv
    else:
        boundB = (w_iv * w_iv) / iv.mpf(J - 1)
    return min(hi(boundA), hi(boundB))


def Rsup_hi(data, N):
    """Certified upper bound for sup_{n=0..N} R_n (O(N) loop, each term individually
    a valid upper bound per the module docstring §1; max of valid upper bounds is
    a valid upper bound for the sup)."""
    w = data['w']
    u = data['u']
    uNp1_hi = hi(u(N + 1))
    best = mp.mpf(0)
    for n in range(0, N + 1):
        un_hi = hi(u(n))
        t1 = tail_bound_hi(N + 1 - n, w)
        t2 = tail_bound_hi(n + N + 2, w)
        val = un_hi * uNp1_hi * (t1 + t2)
        if val > best:
            best = val
    return best


def Ksup_hi(data, N):
    """Exact (not a bound -- proved sup at m=N+1, see module docstring §2) certified
    value of sup_{m>N} K_m, via a single finite (N+1)-term iv sum."""
    u = data['u']
    f = data['f']
    m = N + 1
    um = u(m)
    s = iv.mpf(0)
    for n in range(0, N + 1):
        term = u(n) * (f(m - n) + f(n + m + 1))
        s = s + term
    return hi(um * s)


def schur_C_bound_hi(data, N):
    R = Rsup_hi(data, N)
    K = Ksup_hi(data, N)
    return mp.sqrt(R * K)


def run_one_schur(varpi_val, N, dps=30, verbose=True):
    """Same Davis-Kahan sinTheta pipeline as c1sq_lower_bound.run_one_dk, but with
    ||C|| replaced by the Schur-test bound (schur_C_bound_hi) instead of the old
    global-operator-norm bound. ||H^BB|| unchanged (already O(N^-1), not the
    bottleneck -- team-lead's diagnosis, this round's scope is ||C|| only)."""
    mp.dps = dps
    iv.dps = dps
    t0 = time.time()
    data = build_data(varpi_val, N)
    H, dim = data['H'], data['dim']
    A = to_float_matrix(H, dim)
    evals, evecs = np.linalg.eigh(A)
    lam1_approx = evals[-1]; v1 = evecs[:, -1]
    lam2_approx = evals[-2]

    mu1_lo = lo(rayleigh_lower(H, dim, v1.tolist()))
    la_iv = iv.mpf(mu1_lo)

    r1 = hi(residual_norm(H, dim, v1.tolist(), la_iv))

    c_upper, margin_used = find_lambda2_upper(H, dim, lam2_approx, lam1_approx - lam2_approx)
    lam2_HN_hi = hi(c_upper)

    gap1_lo = mu1_lo - lam2_HN_hi
    if gap1_lo <= r1:
        sinT1 = mp.mpf(1)
    else:
        sinT1 = r1 / (gap1_lo - r1)
        if sinT1 > 1: sinT1 = mp.mpf(1)

    HBB_bound, C_bound_old, _ = tail_bounds(data, N)
    C_bound_schur = schur_C_bound_hi(data, N)
    Ehi = hi(HBB_bound) + C_bound_schur
    lam2A_hi = max(lam2_HN_hi, mp.mpf(0))
    delta_lo = mu1_lo - lam2A_hi
    if delta_lo <= Ehi:
        sinT2 = mp.mpf(1)
    else:
        sinT2 = Ehi / (delta_lo - Ehi)
        if sinT2 > 1: sinT2 = mp.mpf(1)

    sinT = min(mp.mpf(1), sinT1 + sinT2)
    cosT = mp.sqrt(1 - sinT ** 2) if sinT < 1 else mp.mpf(0)

    v0 = abs(v1[0])
    val = cosT * mp.mpf(v0) - sinT
    c1sq_lo = max(mp.mpf(0), val) ** 2

    if verbose:
        print(f"varpi={varpi_val} N={N} dim={dim} time={time.time()-t0:.1f}s")
        print(f"  mu1_lo={mu1_lo}  lam2_HN_hi={lam2_HN_hi}  r1={r1}  sinTheta1={sinT1}")
        print(f"  ||H^BB||<={hi(HBB_bound)}   ||C||_old<={hi(C_bound_old)}   ||C||_schur<={C_bound_schur}"
              f"   (improvement {float(hi(C_bound_old)/C_bound_schur):.2f}x)")
        print(f"  Ehi(schur)={Ehi}  delta_lo={delta_lo}  sinTheta2={sinT2}")
        print(f"  sinTheta_total<={sinT}  cosTheta>={cosT}  |v1[0]|={v0}")
        print(f"  c1sq_lower(cert, Schur route) = {c1sq_lo}")
    return dict(varpi=varpi_val, N=N, mu1_lo=mu1_lo, sinT=sinT, c1sq_lo=c1sq_lo, v0=v0,
                C_old=hi(C_bound_old), C_schur=C_bound_schur)


if __name__ == '__main__':
    varpi_val = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    run_one_schur(varpi_val, N)
