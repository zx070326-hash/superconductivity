r"""
Check: what lambda_* do the first-hand LaH10 r=w2/wlog values actually imply?
Bridge (leading order): r = w2/wlog ~ 1 + sigma^2, sigma = std/mean of P(w).
But lambda_* needs the full shape. Bracket by fitting a shape family to the
target r, then compute EXACT lambda_*[P] = 1/max_s int g(s w) P(w) dw.
"""
import math
import mpmath as mp
mp.mp.dps = 25

A = lambda x: mp.im(mp.polygamma(1, 1j*x))
B = lambda x: mp.re(mp.polygamma(2, 1j*x))
def g(x):
    x = mp.mpf(x)
    return 6*x**2 + 12*x**3*A(x) + 6*x**4*B(x)
def gp(x):
    x = mp.mpf(x)
    return 12*x + 36*x**2*A(x) + 36*x**3*B(x) - 6*x**4*mp.im(mp.polygamma(3,1j*x))
xmax = mp.findroot(gp, mp.mpf('0.33')); gmax = g(xmax)
kappa = 5.603455916

def moments(ws, ps):
    pts = [(mp.mpf(w), mp.mpf(p)) for w, p in zip(ws, ps) if p > 0 and w > 0]
    tot = sum(p for _, p in pts); pts = [(w, p/tot) for w, p in pts]
    m1 = sum(p*w for w, p in pts)
    m2 = sum(p*w*w for w, p in pts)
    mln = sum(p*mp.log(w) for w, p in pts)
    wlog = mp.e**mln; w2 = mp.sqrt(m2)
    sigma = mp.sqrt(m2 - m1*m1)/m1
    return sigma, w2/wlog          # sigma, r=w2/wlog (both from SAME P weight)

def lam_star(ws, ps):
    pts = [(mp.mpf(w), mp.mpf(p)) for w, p in zip(ws, ps) if p > 0 and w > 0]
    tot = sum(p for _, p in pts); pts = [(w, p/tot) for w, p in pts]
    mean = sum(w*p for w, p in pts)
    G = lambda s: sum(p*g(mp.mpf(s)*w) for w, p in pts)
    s0 = mp.findroot(lambda s: mp.diff(G, s), mp.mpf(xmax)/mean)
    return 1/G(s0)

def twopeak(split): return [1-split, 1+split], [1.0, 1.0]
def gauss(rel, N=801, span=6):
    ws = [1 + rel*(-span + 2*span*k/(N-1)) for k in range(N)]
    ps = [math.exp(-0.5*((w-1)/rel)**2) if w > 0 else 0.0 for w in ws]
    return ws, ps

def fit_r(builder, target_r, lo, hi):
    f = lambda p: float(moments(*builder(p))[1]) - target_r
    for _ in range(60):
        mid = (lo+hi)/2
        if f(lo)*f(mid) <= 0: hi = mid
        else: lo = mid
    return (lo+hi)/2

LaH10 = [(129, 1.194), (163, 1.075), (214, 1.105), (264, 1.095)]
print(f"Einstein lambda_*=1/gmax={float(1/gmax):.4f}   Debye lambda_*~4.72")
print(f"{'P':>5}{'r(1st)':>8}{'sig=sqr(r-1)':>13}{'2pk:sig':>9}{'lam*':>8}{'gauss:sig':>10}{'lam*':>8}{'expand@sig':>11}")
for P, r in LaH10:
    sig_approx = math.sqrt(r-1)
    sp = fit_r(twopeak, r, 1e-4, 0.95)
    ws, ps = twopeak(sp); sg2, r2 = moments(ws, ps); L2 = lam_star(ws, ps)
    gr = fit_r(gauss, r, 1e-3, 0.6)
    ws, ps = gauss(gr); sgg, rg = moments(ws, ps); Lg = lam_star(ws, ps)
    exp2 = (1/gmax)*(1 + kappa/2*sg2**2)
    print(f"{P:>5}{r:>8.3f}{sig_approx:>13.3f}{float(sg2):>9.3f}{float(L2):>8.3f}{float(sgg):>10.3f}{float(Lg):>8.3f}{float(exp2):>11.3f}")
