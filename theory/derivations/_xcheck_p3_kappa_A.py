r"""Adversarial cross-check of P3 §4 kappa — part A (symbolic + x_max/kappa + perturbation structure).
Independent of p3_lambda_star_kappa.py. Run: python theory/derivations/_xcheck_p3_kappa_A.py"""
import time
import mpmath as mp
import sympy as sp
T0 = time.time()

# ---------------- Part 1: symbolic (sympy) ----------------
# For real x: conj(psi^(n)(ix)) = psi^(n)(-ix)  =>  Im/Re psi^(n)(ix) written as analytic combos, so sympy can differentiate.
x = sp.symbols('x', positive=True)
def IM(n): return (sp.polygamma(n, sp.I*x) - sp.polygamma(n, -sp.I*x))/(2*sp.I)
def RE(n): return (sp.polygamma(n, sp.I*x) + sp.polygamma(n, -sp.I*x))/2
A, B, C, D, E = IM(1), RE(2), IM(3), RE(4), IM(5)
g_s     = 6*x**2 + 12*x**3*A + 6*x**4*B
lead_g1 = 12*x + 36*x**2*A + 36*x**3*B - 6*x**4*C
lead_g2 = 12 + 72*x*A + 144*x**2*B - 60*x**3*C - 6*x**4*D
my_g3   = 72*A + 360*x*B - 324*x**2*C - 84*x**3*D + 6*x**4*E
def Z(e): return sp.simplify(sp.expand(e))
print("== Part 1: sympy symbolic identities (0 => confirmed)")
print("  A'-B     :", Z(sp.diff(A, x) - B))
print("  B'+C     :", Z(sp.diff(B, x) + C))
print("  C'-D     :", Z(sp.diff(C, x) - D))
print("  D'+E     :", Z(sp.diff(D, x) + E))
print("  g' -lead :", Z(sp.diff(g_s, x) - lead_g1))
print("  g''-lead :", Z(sp.diff(g_s, x, 2) - lead_g2))
print("  g'''-mine:", Z(sp.diff(g_s, x, 3) - my_g3))

# ---------------- Part 2: numeric (mpmath) ----------------
mp.mp.dps = 30
def Ipg(n, t): return mp.im(mp.polygamma(n, mp.mpc(0, t)))
def Rpg(n, t): return mp.re(mp.polygamma(n, mp.mpc(0, t)))
def g(t):  t = mp.mpf(t); return 6*t**2 + 12*t**3*Ipg(1, t) + 6*t**4*Rpg(2, t)
def g1(t): t = mp.mpf(t); return 12*t + 36*t**2*Ipg(1, t) + 36*t**3*Rpg(2, t) - 6*t**4*Ipg(3, t)
def g2(t): t = mp.mpf(t); return 12 + 72*t*Ipg(1, t) + 144*t**2*Rpg(2, t) - 60*t**3*Ipg(3, t) - 6*t**4*Rpg(4, t)
def g3(t): t = mp.mpf(t); return 72*Ipg(1, t) + 360*t*Rpg(2, t) - 324*t**2*Ipg(3, t) - 84*t**3*Rpg(4, t) + 6*t**4*Ipg(5, t)

print("\n== Part 2a: closed forms vs mp.diff ground truth (dps=30)")
for xs in ['0.1', '0.3273', '0.7', '2.0']:
    t = mp.mpf(xs)
    print("  x=%-7s |dg'|=%s |dg''|=%s |dg'''|=%s" % (xs, mp.nstr(abs(g1(t)-mp.diff(g, t)), 3),
          mp.nstr(abs(g2(t)-mp.diff(g, t, 2)), 3), mp.nstr(abs(g3(t)-mp.diff(g, t, 3)), 3)))

print("\n== Part 2b: shape of g")
grid = [mp.mpf(k)/200 for k in range(1, 401)] + [mp.mpf(k) for k in range(3, 21)]
sg = [mp.sign(g1(t)) for t in grid]
ch = [(mp.nstr(grid[i], 4), mp.nstr(grid[i+1], 4)) for i in range(len(grid)-1) if sg[i] != sg[i+1]]
print("  sign changes of g' on (0.005..20):", ch)
print("  g(0.05)=%s vs 6x^2=%s ; g(1)=%s g(3)=%s g(20)=%s g(200)=%s" % (
      mp.nstr(g('0.05'), 6), mp.nstr(6*mp.mpf('0.05')**2, 6), mp.nstr(g(1), 6), mp.nstr(g(3), 6), mp.nstr(g(20), 8), mp.nstr(g(200), 8)))

lo, hi = mp.mpf('0.2'), mp.mpf('0.5')
assert g1(lo) > 0 and g1(hi) < 0
for _ in range(110):
    mid = (lo+hi)/2
    if g1(mid) > 0: lo = mid
    else: hi = mid
xmax = (lo+hi)/2
xmax_num = mp.findroot(lambda s: mp.diff(g, s), mp.mpf('0.3'))
gmax, g2m, g3m = g(xmax), g2(xmax), g3(xmax)
g4m = mp.diff(g, xmax, 4)
kappa = -xmax**2*g2m/gmax
lam0 = 1/gmax
print("\n== Part 2c: x_max / g_max / kappa")
print("  x_max  bisection(closed g')=%s  findroot(numeric dg)=%s  lit 0.3273" % (mp.nstr(xmax, 12), mp.nstr(xmax_num, 12)))
print("  g_max=%s (lit 0.2709)   1/g_max=%s (lit 3.6915)" % (mp.nstr(gmax, 12), mp.nstr(lam0, 12)))
print("  g''(x_max)=%s  g'''(x_max)=%s  g''''(x_max)=%s" % (mp.nstr(g2m, 10), mp.nstr(g3m, 10), mp.nstr(g4m, 10)))
print("  kappa=%s   kappa/2=%s" % (mp.nstr(kappa, 10), mp.nstr(kappa/2, 10)))
print("  x^2 g''/g, x^3 g'''/g, x^4 g''''/g at x_max = %s, %s, %s" % (
      mp.nstr(xmax**2*g2m/gmax, 8), mp.nstr(xmax**3*g3m/gmax, 8), mp.nstr(xmax**4*g4m/gmax, 8)))

def maximize(Gf, s_guess, nscan=17, span=4, reltol=mp.mpf('1e-8')):
    ss = [s_guess*mp.power(span, mp.mpf(2*k)/(nscan-1) - 1) for k in range(nscan)]
    vs = [Gf(s) for s in ss]
    i = max(range(nscan), key=lambda k: vs[k])
    edge = i in (0, nscan-1)
    a, b = ss[max(i-1, 0)], ss[min(i+1, nscan-1)]
    phi = (mp.sqrt(5)-1)/2
    c, d = b - phi*(b-a), a + phi*(b-a)
    fc, fd = Gf(c), Gf(d)
    while (b-a) > reltol*(a+b)/2:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - phi*(b-a); fc = Gf(c)
        else:
            a, c, fc = c, d, fd
            d = a + phi*(b-a); fd = Gf(d)
    s = (a+b)/2
    return s, Gf(s), edge

print("\n== Part 3: symmetric 2-peak (sigma=eps exactly). ratio->1 confirms G_max-g_max = (1/2) sigma^2 x_max^2 g''(x_max); coef->kappa/2")
print("  %6s %14s %20s %18s" % ("eps", "ratio", "(lam/lam0-1)/eps^2", "(coef-k/2)/eps^2"))
for es in ['0.005', '0.01', '0.02', '0.04', '0.08']:
    eps = mp.mpf(es)
    Gf = lambda s, e=eps: (g(s*(1-e)) + g(s*(1+e)))/2
    s, G, _ = maximize(Gf, xmax)
    ratio = (G-gmax)/(eps**2*xmax**2*g2m/2)
    coef = (1/G/lam0 - 1)/eps**2
    print("  %6s %14s %20s %18s" % (es, mp.nstr(ratio, 8), mp.nstr(coef, 8), mp.nstr((coef-kappa/2)/eps**2, 6)))
print("  kappa/2 =", mp.nstr(kappa/2, 8))
print("\nelapsed %.1f s" % (time.time()-T0))
