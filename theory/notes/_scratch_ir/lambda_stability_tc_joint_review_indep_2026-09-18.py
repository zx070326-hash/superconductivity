r"""
Independent re-derivation/verification (main-session review, not imported from
worktree-agent-a681a88291a430268's scripts) of the claimed result:

    d ln Tc_joint / d sigma^2 |_0  = +1.88   (symmetric two-point spectrum,
    N=1 KAY truncation + Semenok lambda_*[P] stability kernel, self-consistent)

and a follow-up test of the skewness (gamma_3) term that the original agent
flagged as untested.

All kernels re-typed from primary sources independently:
  - Semenok stability kernel g(x) = 6x^2 + 12x^3 Im psi'(ix) + 6x^4 Re psi''(ix)
    (theory/derivations/P3_tc_max_pressure.md Eq. S12, cross-checked there
    against Semenok arXiv:2407.12922 Eq. S12).
  - KAY N=1 truncation h1(w) = w^2/(w^2+1), w = omega/(2 pi T)
    (theory/derivations/P4_rigorous_tc_bound.md Sec.1, KAY Thm 4 eq.10).

lambda_*[P] = 1 / max_{u=2piT} int g(omega/u) P(domega)
Tc self-consistency (N=1 lower-bound operator): lambda * <h1>(P,Tc) = 1

Joint (stability-constrained) ceiling: at each spectrum shape, set
lambda := lambda_*[P] (saturate the stability bound), then solve self-consistently
for Tc via the N=1 truncated KAY operator. Track Tc(shape-parameter).
"""
import mpmath as mp
mp.mp.dps = 30


# ---------------------------------------------------------------------------
# 1. Semenok stability kernel g(x), independently re-typed from Eq. S12
# ---------------------------------------------------------------------------
def g_kernel(x):
    x = mp.mpf(x)
    ix = mp.mpc(0, 1) * x
    psi1 = mp.polygamma(1, ix)
    psi2 = mp.polygamma(2, ix)
    return 6*x**2 + 12*x**3*mp.im(psi1) + 6*x**4*mp.re(psi2)


def find_gmax_xmax():
    # locate the single peak by bracketing + findroot on g' via numerical diff
    x0 = mp.findroot(lambda x: mp.diff(g_kernel, x), mp.mpf('0.33'))
    return x0, g_kernel(x0)


xmax, gmax = find_gmax_xmax()
lam0 = 1/gmax
print("=== sanity: Einstein lambda_* from independent g(x) reimplementation ===")
print(f"x_max = {mp.nstr(xmax,10)}  (repo/lit: 0.3273214809)")
print(f"g_max = {mp.nstr(gmax,10)}  (repo/lit: 0.2708931513)")
print(f"lambda_0 = 1/g_max = {mp.nstr(lam0,10)}  (repo/lit Semenok Einstein: 3.6915)")

# kappa_S = -xmax^2 g''(xmax)/gmax via numerical second derivative (independent
# of the repo's hand-derived closed form -- cross-check, not reuse)
gpp_xmax = mp.diff(g_kernel, xmax, 2)
kappa_S = -xmax**2 * gpp_xmax / gmax
print(f"kappa_S (numeric 2nd deriv) = {mp.nstr(kappa_S,10)}  (repo rational-exact-ish value: 5.6035)")
print()


# ---------------------------------------------------------------------------
# 2. KAY N=1 kernel h1(w) = w^2/(w^2+1); kappa_1(wbar) closed form
#    independently re-derived: h1'(w)=2w/(w^2+1)^2, h1''(w)=(2-6w^2)/(w^2+1)^3
#    kappa_1(wbar) := wbar^2 h1''(wbar)/h1(wbar) = (2-6wbar^2)/(wbar^2+1)^2
# ---------------------------------------------------------------------------
def h1_kernel(w):
    w = mp.mpf(w)
    return w**2/(w**2+1)


def kappa1_closed(wbar):
    wbar = mp.mpf(wbar)
    return (2 - 6*wbar**2)/(wbar**2+1)**2


# independent finite-difference cross-check of kappa1 closed form
wbar_test = mp.mpf('0.6095')
h1pp_num = mp.diff(h1_kernel, wbar_test, 2)
kappa1_num = wbar_test**2 * h1pp_num / h1_kernel(wbar_test)
kappa1_cf = kappa1_closed(wbar_test)
print("=== sanity: kappa_1(wbar) closed form vs numeric 2nd derivative ===")
print(f"wbar={float(wbar_test)}  kappa1_closed={mp.nstr(kappa1_cf,10)}  kappa1_numeric={mp.nstr(kappa1_num,10)}")
print()


# ---------------------------------------------------------------------------
# 3. Joint self-consistent solve: symmetric two-point spectrum
#    P = 1/2 delta_{1-s} + 1/2 delta_{1+s},  <omega>=1, sigma=s, gamma_3=0
# ---------------------------------------------------------------------------
def lambda_star_two_point(s, dps_search=200):
    s = mp.mpf(s)
    w1, w2 = 1-s, 1+s

    def negF(u):
        if u <= 0:
            return mp.mpf('1e10')
        return -(g_kernel(w1/u) + g_kernel(w2/u))/2

    a, b = mp.mpf('1e-3'), mp.mpf(20)
    gr = (mp.sqrt(5)-1)/2
    c = b - gr*(b-a)
    d = a + gr*(b-a)
    for _ in range(dps_search):
        if negF(c) < negF(d):
            b = d
        else:
            a = c
        c = b - gr*(b-a)
        d = a + gr*(b-a)
    u_star = (a+b)/2
    return 1/(-negF(u_star)), u_star


def solve_tc_two_point(s, lam, iters=200):
    s = mp.mpf(s)
    w1, w2 = 1-s, 1+s

    def F(u):
        return (h1_kernel(w1/u) + h1_kernel(w2/u))/2 - 1/lam

    lo, hi = mp.mpf('1e-4'), mp.mpf(50)
    for _ in range(iters):
        mid = (lo+hi)/2
        if F(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


print("=== joint self-consistent solve: symmetric two-point spectrum ===")
s0 = mp.mpf('1e-6')
lam_s0, _ = lambda_star_two_point(s0)
u0 = solve_tc_two_point(s0, lam_s0)
print(f"baseline (s->0): lambda_*={mp.nstr(lam_s0,10)}  u0=2*pi*Tc0={mp.nstr(u0,12)}  wbar0=1/u0={mp.nstr(1/u0,10)}")
print()

svals = [mp.mpf(x) for x in ['0.002', '0.005', '0.01', '0.02', '0.05', '0.10', '0.15', '0.20', '0.30']]
print(f"{'s':>8}{'lambda_*(s)':>16}{'u(s)=2piTc':>16}{'Tc(s)/Tc(0)':>16}{'dlnTc/ds^2':>14}")
rows = []
for s in svals:
    lam_s, _ = lambda_star_two_point(s)
    u_s = solve_tc_two_point(s, lam_s)
    ratio = u_s/u0
    dlnTc_ds2 = mp.log(ratio)/s**2
    rows.append((s, lam_s, u_s, ratio, dlnTc_ds2))
    print(f"{float(s):>8.4f}{float(lam_s):>16.8f}{float(u_s):>16.8f}{float(ratio):>16.10f}{float(dlnTc_ds2):>14.6f}")

# Richardson-style extrapolation using the two smallest s values (independent
# of any particular s choice -- avoids relying on a single finite-diff point)
s_a, s_b = svals[0], svals[1]
d_a = rows[0][4]
d_b = rows[1][4]
# naive Richardson assuming O(s^2) leading correction to the ratio d(s^2)
# (linear extrapolation in s^2 of the finite-difference estimate itself)
extrap = d_a + (d_a - d_b) * (s_a**2) / (s_b**2 - s_a**2)
print(f"\nRichardson-style extrapolation (s={float(s_a)},{float(s_b)}): dlnTc/ds^2|0 ~ {mp.nstr(extrap,8)}")
print()


# ---------------------------------------------------------------------------
# 4. Implicit-function-theorem closed-form cross-check
#    G(u,sigma) := lambda(sigma) * <h1>(1/u, sigma) - 1 = 0
#    lambda(sigma) = lambda0*(1 + 0.5*kappa_S*sigma^2)
#    <h1>(wbar,sigma) = h1(wbar)*(1 + 0.5*kappa1(wbar)*sigma^2 + O(sigma^3))
# ---------------------------------------------------------------------------
wbar0 = 1/u0
kappa1_0 = kappa1_closed(wbar0)
dG_dsigma2 = mp.mpf('0.5')*(kappa_S + kappa1_0)

hp = lambda w: 2*w/(w**2+1)**2
dG_du = -lam0 * hp(wbar0) * wbar0**2  # d/du[h1(1/u)] = h1'(wbar)*(-1/u^2) = -h1'(wbar)*wbar^2

du_dsigma2 = -dG_dsigma2/dG_du
dlnu_dsigma2 = du_dsigma2/u0

print("=== implicit function theorem closed-form cross-check ===")
print(f"wbar0={mp.nstr(wbar0,10)}  kappa_S={mp.nstr(kappa_S,10)}  kappa1(wbar0)={mp.nstr(kappa1_0,10)}")
print(f"dG/dsigma^2|0 = {mp.nstr(dG_dsigma2,10)}")
print(f"dG/du|0       = {mp.nstr(dG_du,10)}")
print(f"du/dsigma^2   = {mp.nstr(du_dsigma2,10)}")
print(f"d ln Tc/dsigma^2|0 (closed form) = {mp.nstr(dlnu_dsigma2,10)}")
print()
print(f"COMPARISON: direct-solve extrapolation = {mp.nstr(extrap,8)}  vs  closed-form = {mp.nstr(dlnu_dsigma2,8)}")


# ---------------------------------------------------------------------------
# 5. Skewness (gamma_3) test: asymmetric three-point spectrum with small,
#    controlled sigma and a tunable gamma_3, gamma_3 != 0 while sigma is held
#    fixed-small so we can extract the gamma_3-dependence of the joint ceiling
#    at fixed sigma^2 (i.e. test the sign/magnitude of the sigma^3 term).
# ---------------------------------------------------------------------------
def make_three_point(sigma, skew_frac):
    """
    Three-point spectrum with mean 1, controlled sigma and skewness gamma_3.
    Parametrize as P = p_L delta_{1+a} + p_M delta_1 + p_R delta_{1+b}
    with a<0<b chosen so that mean=1 automatically and (sigma,gamma3) match
    target values, by solving a small algebraic system.

    Simpler explicit construction: use a two-parameter family
        omega_1 = 1 - x,  p_1 = q
        omega_2 = 1        p_2 = 1-q-r
        omega_3 = 1 + y,   p_3 = r
    with mean=1 automatically satisfied by  q*(-x) + r*y = 0  =>  x = r*y/q.
    Choose q=r=p (symmetric weights) but asymmetric offsets x != y to get
    nonzero skewness while controlling sigma via an overall scale.
    """
    # asymmetric offsets, equal weights p on each side, symmetric weight leftover in center
    p = mp.mpf('0.25')
    # skew_frac in (-1,1) controls x/y ratio asymmetry; sigma sets the overall scale
    # x*q = y*r (mean-preserving) with q=r=p => x=y actually forces symmetry.
    # Instead use unequal weights: q != r, with x=y=delta fixed ratio, tune q,r,delta
    # to hit (sigma, gamma3) approximately via direct solve below.
    raise NotImplementedError


def three_point_from_moments(sigma, gamma3, w_out=None):
    """
    Build an explicit 3-point distribution on xi = (omega-1) with mean(xi)=0,
    var(xi)=sigma^2, skew(xi)=gamma3, by placing point masses at
    xi in {a, 0, b} (a<0<b) with weights {p, 1-p-q, q} -- 3 unknowns (a,b,p)
    [taking q as a free 4th parameter we fix by symmetry choice] solving the
    3 moment equations mean=0, var=sigma^2, skew=gamma3.

    To keep this closed-form and robust we instead use a standard two-point-
    plus-atom construction: place equal masses p at +-c(1+eps) asymmetrically.
    Concretely: xi in {-c, d} with masses {q, 1-q}, mean=0 => q*(-c)+(1-q)*d=0.
    2nd moment: q c^2 + (1-q) d^2 = sigma^2.
    3rd moment: q(-c)^3 + (1-q) d^3 = gamma3 * sigma^3.
    Two-point (not three) but already asymmetric with a controllable skew:
    solve for (c,d,q) given (sigma,gamma3). Reduces to 1D root find in q.
    """
    sigma = mp.mpf(sigma)
    gamma3 = mp.mpf(gamma3)

    def solve_for_q(q):
        # mean=0: q*(-c) + (1-q)*d = 0  => d = q*c/(1-q)
        # Parametrize by c, solve var eq for c given q, then check skew.
        # var: q c^2 + (1-q) d^2 = sigma^2 ; d = q c/(1-q)
        # => q c^2 + (1-q) * (q c/(1-q))^2 = q c^2 (1 + q/(1-q)) = q c^2 /(1-q) = sigma^2
        # => c^2 = sigma^2 (1-q)/q  => c = sigma*sqrt((1-q)/q)
        c = sigma*mp.sqrt((1-q)/q)
        d = q*c/(1-q)
        skew_val = (q*(-c)**3 + (1-q)*d**3)/sigma**3
        return skew_val - gamma3, c, d

    # bisection on q in (0,1) excluding endpoints
    lo, hi = mp.mpf('0.01'), mp.mpf('0.99')
    flo, _, _ = solve_for_q(lo)
    fhi, _, _ = solve_for_q(hi)
    if flo*fhi > 0:
        # try to bracket differently: skew(q) is monotonic-ish; scan
        found = False
        qs = [mp.mpf(i)/1000 for i in range(10, 991)]
        prev_q, prev_f = qs[0], solve_for_q(qs[0])[0]
        for qv in qs[1:]:
            fv = solve_for_q(qv)[0]
            if prev_f*fv <= 0:
                lo, hi = prev_q, qv
                found = True
                break
            prev_q, prev_f = qv, fv
        if not found:
            raise RuntimeError(f"could not bracket skew root for gamma3={gamma3}")
    for _ in range(100):
        mid = (lo+hi)/2
        fmid, c, d = solve_for_q(mid)
        flo, _, _ = solve_for_q(lo)
        if flo*fmid <= 0:
            hi = mid
        else:
            lo = mid
    q = (lo+hi)/2
    _, c, d = solve_for_q(q)
    omega1, omega2 = 1-c, 1+d
    p1, p2 = q, 1-q
    return [(omega1, p1), (omega2, p2)]


def lambda_star_general(pts, dps_search=200):
    def negF(u):
        if u <= 0:
            return mp.mpf('1e10')
        return -sum(p*g_kernel(w/u) for w, p in pts)

    a, b = mp.mpf('1e-3'), mp.mpf(20)
    gr = (mp.sqrt(5)-1)/2
    c = b - gr*(b-a)
    d = a + gr*(b-a)
    for _ in range(dps_search):
        if negF(c) < negF(d):
            b = d
        else:
            a = c
        c = b - gr*(b-a)
        d = a + gr*(b-a)
    u_star = (a+b)/2
    return 1/(-negF(u_star)), u_star


def solve_tc_general(pts, lam, iters=200):
    def F(u):
        return sum(p*h1_kernel(w/u) for w, p in pts) - 1/lam

    lo, hi = mp.mpf('1e-4'), mp.mpf(50)
    for _ in range(iters):
        mid = (lo+hi)/2
        if F(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


print()
print("=== skewness (gamma_3) test: asymmetric two-point family at fixed small sigma ===")
sigma_fix = mp.mpf('0.05')
print(f"fixed sigma = {float(sigma_fix)}")
print(f"{'gamma3':>10}{'lambda_*':>16}{'u(2piTc)':>16}{'Tc/Tc0':>16}")
for gamma3v in ['-1.5', '-1.0', '-0.5', '0.0', '0.5', '1.0', '1.5']:
    g3 = mp.mpf(gamma3v)
    try:
        pts = three_point_from_moments(sigma_fix, g3)
    except RuntimeError as e:
        print(f"gamma3={gamma3v}: {e}")
        continue
    lam_s, _ = lambda_star_general(pts)
    u_s = solve_tc_general(pts, lam_s)
    ratio = u_s/u0
    print(f"{float(g3):>10.2f}{float(lam_s):>16.8f}{float(u_s):>16.8f}{float(ratio):>16.10f}")

print()
print("Comparing Tc-ratio at gamma3=+1.5 vs gamma3=-1.5 (same |sigma|) isolates the")
print("sign of the sigma^3*gamma3 contribution to the joint ceiling's sensitivity.")

# ---------------------------------------------------------------------------
# 6. Analytic (implicit-function-theorem) estimate of the gamma_3*sigma^3
#    coefficient, using the repo's independently-derived closed forms for
#    c_3 (P3, Semenok kernel) and c_{3,1}(wbar) (P4-A, N=1 KAY kernel).
#    c_3 = 0.8958993147 (P3_tc_max_pressure.md, sympy-cross-checked there)
#    c_{3,1}(wbar) = 4*wbar^2*(wbar^2-1)/(wbar^2+1)^3  (P4_rigorous_tc_bound.md eq. P4A.5)
# ---------------------------------------------------------------------------
print()
print("=== analytic estimate of d ln Tc / d(gamma3*sigma^3)|0 via IFT ===")
c3_P3 = mp.mpf('0.8958993147')  # repo closed-form value, independently re-derivable
c31 = 4*wbar0**2*(wbar0**2-1)/(wbar0**2+1)**3
print(f"c3 (Semenok/P3) = {mp.nstr(c3_P3,10)}")
print(f"c31(wbar0) (KAY N=1/P4-A) = {mp.nstr(c31,10)}")

# at eps := gamma3*sigma^3 -> 0, using lambda0*h(wbar0)=1 identity
coeff_gamma3 = -(c3_P3 + c31) / dG_du / u0
print(f"d ln Tc / d(gamma3 sigma^3)|0 (IFT, analytic-ish) = {mp.nstr(coeff_gamma3,8)}")
print(f"  compare to sigma^2 coefficient +1.8796: ratio = {mp.nstr(coeff_gamma3/dlnu_dsigma2,6)}")
print(f"  => net sign flips only if gamma3*sigma < -{mp.nstr(dlnu_dsigma2/coeff_gamma3,6)} (approx, leading order)")

# ---------------------------------------------------------------------------
# 7. Extended skewness test at larger sigma to see how large |gamma3| needs
#    to be, and whether it's achievable by a physically reasonable spectrum,
#    before the net joint-ceiling response flips sign.
# ---------------------------------------------------------------------------
print()
print("=== extended test: sigma=0.20, wide gamma3 range ===")
sigma_fix2 = mp.mpf('0.20')
lam0_s2, _ = lambda_star_two_point(mp.mpf('1e-6'))
print(f"{'gamma3':>10}{'lambda_*':>16}{'u':>16}{'Tc/Tc0':>16}{'ln(Tc/Tc0)':>16}")
for gamma3v in ['-4.0', '-3.0', '-2.0', '-1.0', '0.0', '1.0', '2.0', '3.0', '4.0']:
    g3 = mp.mpf(gamma3v)
    try:
        pts = three_point_from_moments(sigma_fix2, g3)
    except RuntimeError as e:
        print(f"gamma3={gamma3v}: unreachable ({e})")
        continue
    lam_s, _ = lambda_star_general(pts)
    u_s = solve_tc_general(pts, lam_s)
    ratio = u_s/u0
    print(f"{float(g3):>10.2f}{float(lam_s):>16.8f}{float(u_s):>16.8f}{float(ratio):>16.10f}{float(mp.log(ratio)):>16.8f}")

# ---------------------------------------------------------------------------
# 8. N=2 truncation-dependence check. Re-derives the KAY-II H^(N) matrix
#    construction independently (self-checked against KAY I truncation
#    constants T_c^(N)(2), following the "adopted placement" pinned down and
#    verified in theory/derivations/p4_kappa_N.py -- structure only reused,
#    not the numeric routines) and repeats the symmetric two-point joint
#    self-consistent solve with the N=2 truncated operator instead of N=1.
# ---------------------------------------------------------------------------
print()
print("=== N=2 truncation dependence check (independent re-implementation) ===")


def brk(k, w):
    w = mp.mpf(w)
    return w**2 / (mp.mpf(k)**2 + w**2)


def H2_matrix(w):
    w = mp.mpf(w)
    M = mp.matrix(2, 2)
    for n in range(2):
        for m in range(2):
            nrm = mp.sqrt(mp.mpf(2*n+1)) * mp.sqrt(mp.mpf(2*m+1))
            if n != m:
                val = (brk(n+m+1, w) + brk(abs(n-m), w)) / nrm
            else:
                s = sum(2*brk(k, w) for k in range(1, n+1)) if n >= 1 else mp.mpf(0)
                val = brk(n+m+1, w)/nrm - s/mp.mpf(2*n+1)
            M[n, m] = val
    return M


# self-check: small-w bridge to KAY I g^(2)(2) = (20+sqrt(13819))/108 = 1.2736503963
g2_target = (20 + mp.sqrt(13819))/108
w_small = mp.mpf('0.001')
E, _ = mp.eigsy(H2_matrix(w_small))
rho_over_w2 = max(E)/w_small**2
print(f"bridge check: rho(H2(w->0))/w^2 = {mp.nstr(rho_over_w2,10)}  vs g2(2) target = {mp.nstr(g2_target,10)}")


def h2_top(w):
    E, _ = mp.eigsy(H2_matrix(w))
    return max(E)


def lambda_star_two_point_N2_label(s):
    return lambda_star_two_point(s)  # lambda_* is defined by the stability kernel g, N-independent


def solve_tc_two_point_N2(s, lam, iters=150):
    s = mp.mpf(s)
    w1, w2 = 1-s, 1+s

    def F(u):
        Mavg = (H2_matrix(w1/u) + H2_matrix(w2/u))/2
        E, _ = mp.eigsy(Mavg)
        return max(E) - 1/lam

    lo, hi = mp.mpf('1e-4'), mp.mpf(50)
    for _ in range(iters):
        mid = (lo+hi)/2
        if F(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


s0 = mp.mpf('1e-6')
lam_s0, _ = lambda_star_two_point(s0)
u0_N2 = solve_tc_two_point_N2(s0, lam_s0)
print(f"N=2 baseline: u0={mp.nstr(u0_N2,10)}  (N=1 baseline was {mp.nstr(u0,10)})")

print(f"{'s':>8}{'u_N2(s)':>16}{'Tc_N2(s)/Tc_N2(0)':>20}{'dlnTc/ds^2 (N=2)':>18}")
for sv in ['0.01', '0.02', '0.05', '0.10']:
    s = mp.mpf(sv)
    lam_s, _ = lambda_star_two_point(s)
    u_s = solve_tc_two_point_N2(s, lam_s)
    ratio = u_s/u0_N2
    d = mp.log(ratio)/s**2
    print(f"{float(s):>8.4f}{float(u_s):>16.8f}{float(ratio):>20.10f}{float(d):>18.6f}")
