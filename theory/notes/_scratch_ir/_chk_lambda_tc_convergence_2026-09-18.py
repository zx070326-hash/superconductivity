import mpmath as mp
mp.mp.dps = 30


def g_kernel(x):
    x = mp.mpf(x)
    ix = mp.mpc(0, 1) * x
    psi1 = mp.polygamma(1, ix)
    psi2 = mp.polygamma(2, ix)
    return 6*x**2 + 12*x**3*mp.im(psi1) + 6*x**4*mp.re(psi2)


def h1_kernel(x):
    x = mp.mpf(x)
    return x**2/(x**2+1)


def lambda_star_two_point(s):
    w1, w2 = mp.mpf(1)-s, mp.mpf(1)+s

    def negF(u):
        if u <= 0:
            return mp.mpf('1e10')
        return -(g_kernel(w1/u) + g_kernel(w2/u))/2
    a, b = mp.mpf('1e-3'), mp.mpf(20)
    gr = (mp.sqrt(5)-1)/2
    c = b - gr*(b-a)
    d = a + gr*(b-a)
    for _ in range(300):
        if negF(c) < negF(d):
            b = d
        else:
            a = c
        c = b - gr*(b-a)
        d = a + gr*(b-a)
    u_star = (a+b)/2
    return mp.mpf(1)/(-negF(u_star)), u_star


def solve_tc_two_point(s, lam):
    w1, w2 = mp.mpf(1)-s, mp.mpf(1)+s

    def F(u):
        return (h1_kernel(w1/u) + h1_kernel(w2/u))/2 - 1/lam
    lo, hi = mp.mpf('1e-4'), mp.mpf(50)
    for _ in range(300):
        mid = (lo+hi)/2
        if F(mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo+hi)/2


lam0, u0 = lambda_star_two_point(mp.mpf('1e-6'))
u0b = solve_tc_two_point(mp.mpf('1e-6'), lam0)
print("baseline lam0,u0b:", lam0, u0b)

for sv in ['0.002', '0.005', '0.01', '0.02', '0.05', '0.1']:
    s = mp.mpf(sv)
    lam_s, _ = lambda_star_two_point(s)
    u_s = solve_tc_two_point(s, lam_s)
    ratio = u_s/u0b
    d = mp.log(ratio)/s**2
    print(f"s={float(s):<8} lam_s={float(lam_s):.8f}  Tc_ratio={float(ratio):.10f}  dlnTc/ds2={float(d):.6f}")
