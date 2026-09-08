r"""Adversarial cross-check of P3 §4 kappa — part B (exact lambda_*[P] for model spectra vs expansions).
Independent of p3_lambda_star_kappa.py. Run: python theory/derivations/_xcheck_p3_kappa_B.py"""
import time
import mpmath as mp
T0 = time.time()
mp.mp.dps = 15
def Ipg(n, t): return mp.im(mp.polygamma(n, mp.mpc(0, t)))
def Rpg(n, t): return mp.re(mp.polygamma(n, mp.mpc(0, t)))
def g(t):  t = mp.mpf(t); return 6*t**2 + 12*t**3*Ipg(1, t) + 6*t**4*Rpg(2, t)
def g1(t): t = mp.mpf(t); return 12*t + 36*t**2*Ipg(1, t) + 36*t**3*Rpg(2, t) - 6*t**4*Ipg(3, t)
def g2(t): t = mp.mpf(t); return 12 + 72*t*Ipg(1, t) + 144*t**2*Rpg(2, t) - 60*t**3*Ipg(3, t) - 6*t**4*Rpg(4, t)
def g3(t): t = mp.mpf(t); return 72*Ipg(1, t) + 360*t*Rpg(2, t) - 324*t**2*Ipg(3, t) - 84*t**3*Rpg(4, t) + 6*t**4*Ipg(5, t)
lo, hi = mp.mpf('0.2'), mp.mpf('0.5')
for _ in range(60):
    mid = (lo+hi)/2
    if g1(mid) > 0: lo = mid
    else: hi = mid
xmax = (lo+hi)/2
gmax, g2m, g3m = g(xmax), g2(xmax), g3(xmax)
kappa = -xmax**2*g2m/gmax
lam0 = 1/gmax
print("x_max=%s g_max=%s lam0=%s kappa=%s" % (mp.nstr(xmax, 8), mp.nstr(gmax, 8), mp.nstr(lam0, 8), mp.nstr(kappa, 8)))

def maximize(Gf, s_guess, nscan=17, span=4, reltol=mp.mpf('1e-6')):
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

def stats_from_moments(m1, m2, m3, m4, lnm):
    var = m2 - m1**2; sig = mp.sqrt(var)/m1
    mu3 = m3 - 3*m1*m2 + 2*m1**3
    mu4 = m4 - 4*m1*m3 + 6*m1**2*m2 - 3*m1**4
    skew = mu3/(m1*sig)**3 if sig > 0 else mp.mpf(0)
    kurt = mu4/(m1*sig)**4 if sig > 0 else mp.mpf(0)
    r = mp.sqrt(m2)/mp.exp(lnm)
    sig_a2F = mp.sqrt(m3*m1/m2**2 - 1)       # relative width if one (wrongly) used alpha2F weight (= w*P)
    return m1, sig, skew, kurt, r, sig_a2F

def discrete(pts):
    tot = sum(p for _, p in pts); pts = [(mp.mpf(w), mp.mpf(p)/tot) for w, p in pts]
    mom = [sum(p*w**k for w, p in pts) for k in (1, 2, 3, 4)]
    lnm = sum(p*mp.log(w) for w, p in pts)
    Gf = lambda s: sum(p*g(s*w) for w, p in pts)
    return stats_from_moments(*mom, lnm), Gf

def continuous(P, brk):
    Zn = mp.quad(P, brk)
    Pn = lambda w: P(w)/Zn
    mom = [mp.quad(lambda w, k=k: w**k*Pn(w), brk) for k in (1, 2, 3, 4)]
    lnm = mp.quad(lambda w: mp.log(w)*Pn(w), brk)
    Gf = lambda s: mp.quad(lambda w: Pn(w)*g(s*w), brk)
    return stats_from_moments(*mom, lnm), Gf

def gauss(rel, mean=1):
    sd = rel*mean
    a = max(mp.mpf(0), mean-7*sd); b = mean+7*sd
    return continuous(lambda w: mp.exp(-((w-mean)/sd)**2/2), [a, mean, b])

def power_law(n):   # P ~ w^n on (0,1]; n=1 is Debye (alpha2F ~ w^2 => P = 2 a2F/(lam w) ~ w)
    return continuous(lambda w: w**n, [0, 1])

def hydride(w_ac, w_cut=mp.mpf('0.35'), rel=mp.mpf('0.10')):
    # acoustic Debye-like tail P~w on (0,w_cut] carrying P-weight (=lambda-share) w_ac, + optical Gaussian at 1 (rel width 0.10)
    sd = rel
    tail = lambda w: 2*w/w_cut**2 if w <= w_cut else mp.mpf(0)
    Zg = mp.sqrt(2*mp.pi)*sd
    opt = lambda w: mp.exp(-((w-1)/sd)**2/2)/Zg if w > mp.mpf('0.4') else mp.mpf(0)
    P = lambda w: w_ac*tail(w) + (1-w_ac)*opt(w)
    return continuous(P, [0, w_cut, mp.mpf('0.4'), 1, 1+7*sd])

cases = [
 ("Einstein",                 discrete([(1, 1)])),
 ("2-peak eps=0.10",          discrete([(0.9, 1), (1.1, 1)])),
 ("2-peak eps=0.20",          discrete([(0.8, 1), (1.2, 1)])),
 ("2-peak eps=0.30",          discrete([(0.7, 1), (1.3, 1)])),
 ("3-peak (1,2,1) sig=0.2",   discrete([(1-mp.mpf('0.2')*mp.sqrt(2), 1), (1, 2), (1+mp.mpf('0.2')*mp.sqrt(2), 1)])),
 ("asym2 skew=+1.5 sig=0.2",  discrete([(1.4, 0.2), (0.9, 0.8)])),
 ("asym2 skew=-1.5 sig=0.2",  discrete([(0.6, 0.2), (1.1, 0.8)])),
 ("Gauss rel=0.10",           gauss(mp.mpf('0.10'))),
 ("Gauss rel=0.20",           gauss(mp.mpf('0.20'))),
 ("Gauss rel=0.30",           gauss(mp.mpf('0.30'))),
 ("P~w^4",                    power_law(4)),
 ("P~w^2",                    power_law(2)),
 ("Debye P~w (a2F~w^2)",      power_law(1)),
 ("P~w^0.84",                 power_law(mp.mpf('0.84'))),
 ("hydride w_ac=0.04",        hydride(mp.mpf('0.04'))),
 ("hydride w_ac=0.08",        hydride(mp.mpf('0.08'))),
 ("hydride w_ac=0.14",        hydride(mp.mpf('0.14'))),
 ("hydride w_ac=0.20",        hydride(mp.mpf('0.20'))),
]
print("\n== Part 4: exact lambda_*[P] = 1/max_s int g(s w) P(w) dw  vs expansions   [lam0=1/g_max=%s, kappa/2=%s]" % (mp.nstr(lam0, 7), mp.nstr(kappa/2, 6)))
hdr = ("spectrum", "sigma_P", "sig_a2F", "skew", "kurt", "r", "exact", "2nd", "err2", "3rd(+skew)", "err3", "2nd via r-1", "errR", "edge")
print("  %-26s%8s%8s%7s%6s%8s%9s%9s%8s%11s%8s%12s%8s%5s" % hdr)
for name, ((m1, sig, skew, kurt, r, siga), Gf) in cases:
    t1 = time.time()
    s, G, edge = maximize(Gf, xmax/m1)
    ex = 1/G
    lam2 = lam0*(1 + kappa/2*sig**2)
    G3 = gmax + sig**2*xmax**2*g2m/2 + skew*sig**3*xmax**3*g3m/6
    lam3 = 1/G3
    lamr = lam0*(1 + kappa/2*(r-1))
    f = float
    print("  %-26s%8.4f%8.4f%7.2f%6.2f%8.4f%9.5f%9.5f%7.2f%%%11.5f%7.2f%%%12.5f%7.2f%%%5s  (%.0fs)" % (
        name, f(sig), f(siga), f(skew), f(kurt), f(r), f(ex), f(lam2), 100*f((lam2-ex)/ex), f(lam3), 100*f((lam3-ex)/ex),
        f(lamr), 100*f((lamr-ex)/ex), "!" if edge else "", time.time()-t1), flush=True)

print("\n== Part 5: Allen-Dynes omega_log, omega_2 computed from alpha2F(w)=w^2 (Debye) directly vs P-moments")
a2F = lambda w: w**2
lam = 2*mp.quad(lambda w: a2F(w)/w, [0, 1])
wlog = mp.exp(2/lam*mp.quad(lambda w: mp.log(w)*a2F(w)/w, [0, 1]))
w2 = mp.sqrt(2/lam*mp.quad(lambda w: w*a2F(w), [0, 1]))
(m1, sig, skew, kurt, r, siga), _ = power_law(1)
print("  AD from a2F: lambda=%s omega_log=%s omega_2=%s r=%s" % (mp.nstr(lam, 8), mp.nstr(wlog, 8), mp.nstr(w2, 8), mp.nstr(w2/wlog, 8)))
print("  P-moments  : r=%s sigma_P=%s ; a2F-weighted sigma=%s ; analytic Debye: sigma_P=sqrt(1/8)=%s, r=e^(1/2)/sqrt2=%s, r-1=%s vs sigma^2=0.125" % (
      mp.nstr(r, 8), mp.nstr(sig, 8), mp.nstr(siga, 8), mp.nstr(mp.sqrt(mp.mpf(1)/8), 8), mp.nstr(mp.exp(mp.mpf('0.5'))/mp.sqrt(2), 8), mp.nstr(mp.exp(mp.mpf('0.5'))/mp.sqrt(2)-1, 6)))
print("\nelapsed %.1f s" % (time.time()-T0))
