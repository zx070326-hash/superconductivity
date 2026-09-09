r"""
P3 §4 — universality of the lambda_*(sigma, gamma3) expansion.

Claim to test:  lambda_*/lambda0 = 1 + (kappa/2) sigma^2 + c3 gamma3 sigma^3 + O(sigma^4),
with kappa = -x_max^2 g''(x_max)/g_max  and  c3 = -x_max^3 g'''(x_max)/(6 g_max)
being PURE NUMBERS independent of the spectral shape (peak count, kernel family).

Reason (moment/saddle expansion): with omega = mu(1+xi), t = s*mu,
    G(t) = \int g(t(1+xi)) P dxi
         = g(t) + (1/2) t^2 g''(t) sigma^2 + (1/6) t^3 g'''(t) gamma3 sigma^3 + O(sigma^4).
Maximizing over t: the saddle shift is O(sigma^2), and since g'(x_max)=0 the shift feeds
back only at O(sigma^4) (envelope theorem). Hence to O(sigma^3),
    G_max = g_max [ 1 + (x_max^2 g''/g_max)(sigma^2/2) + (x_max^3 g'''/g_max)(gamma3 sigma^3/6) ],
and lambda_* = 1/G_max gives the expansion above. The coefficients depend ONLY on the
kernel g at its own maximum and on the 2nd/3rd moments of P -- NOT on how many peaks P has.

So the test: build spectra with the SAME sigma (and same gamma3) but different shapes /
peak counts; the expansion prediction must track exact lambda_* equally well for all of them,
and the O(sigma^4) residual is the FIRST place shape (kurtosis, peak count) enters.
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
    return 12*x + 36*x**2*A(x) + 36*x**3*B(x) - 6*x**4*mp.im(mp.polygamma(3, 1j*x))

xmax = mp.findroot(gp, mp.mpf('0.33'))
gmax = g(xmax)
# derivatives at the max via robust numerical differentiation (cross-checked vs closed form elsewhere)
gpp_max  = mp.diff(g, xmax, 2)
gppp_max = mp.diff(g, xmax, 3)
kappa = -xmax**2 * gpp_max / gmax
c3    = -xmax**3 * gppp_max / (6*gmax)
lam0  = 1/gmax
print(f"x_max={mp.nstr(xmax,8)}  g_max={mp.nstr(gmax,8)}  lam0={mp.nstr(lam0,8)}")
print(f"kappa = {mp.nstr(kappa,8)}  (1/2 kappa={mp.nstr(kappa/2,7)})   c3 = {mp.nstr(c3,8)}")
print()

def moments(ws, ps):
    pts = [(mp.mpf(w), mp.mpf(p)) for w, p in zip(ws, ps) if p > 0 and w > 0]
    tot = sum(p for _, p in pts); pts = [(w, p/tot) for w, p in pts]
    mu = sum(w*p for w, p in pts)
    c2 = sum(p*(w-mu)**2 for w, p in pts)
    c3m = sum(p*(w-mu)**3 for w, p in pts)
    sig = mp.sqrt(c2)/mu
    gam3 = c3m/c2**mp.mpf('1.5')          # standard skewness (mu cancels)
    return sig, gam3

def lam_star(ws, ps):
    pts = [(mp.mpf(w), mp.mpf(p)) for w, p in zip(ws, ps) if p > 0 and w > 0]
    tot = sum(p for _, p in pts); pts = [(w, p/tot) for w, p in pts]
    mu = sum(w*p for w, p in pts)
    G = lambda s: sum(p*g(mp.mpf(s)*w) for w, p in pts)
    s0 = mp.findroot(lambda s: mp.diff(G, s), xmax/mu)
    return 1/G(s0)

def expand(sig, gam3):
    return lam0*(1 + kappa/2*sig**2 + c3*gam3*sig**3)

# ---- spectra: several shapes, some deliberately matched in sigma ----
def sympeaks(offsets, weights):          # symmetric/general discrete peaks about 1
    ws = [1+o for o in offsets]; return ws, list(weights)

cases = {
  # --- symmetric families at (near) sigma=0.20 : should agree through O(sigma^3), gamma3=0 ---
  "2-peak sym  s=.20":      sympeaks([-0.20, 0.20], [1, 1]),
  "3-peak sym  (eq wt)":    sympeaks([-0.2449, 0.0, 0.2449], [1, 1, 1]),   # tuned to sigma~0.20
  "3-peak sym  (heavy mid)":sympeaks([-0.2828, 0.0, 0.2828], [1, 2, 1]),   # tuned to sigma~0.20
  "4-peak sym":             sympeaks([-0.24, -0.08, 0.08, 0.24], [1,1,1,1]),
  # --- asymmetric (nonzero gamma3): tests the c3*gamma3*sigma^3 term across shapes ---
  "3-peak asym A":          sympeaks([-0.15, 0.0, 0.35], [2, 1, 1]),
  "3-peak asym B":          sympeaks([-0.30, 0.05, 0.20], [1, 1, 2]),
  "2-peak asym":            sympeaks([-0.10, 0.30], [3, 1]),
  # --- a wider symmetric case to see O(sigma^4) grow ---
  "3-peak sym  s=.30":      sympeaks([-0.3674, 0.0, 0.3674], [1, 1, 1]),
}

print(f"{'spectrum':24s}{'sigma':>8s}{'gamma3':>9s}{'exact':>10s}{'2nd':>10s}{'2nd+3rd':>10s}{'err2':>8s}{'err23':>8s}")
for name, (ws, ps) in cases.items():
    sig, gam3 = moments(ws, ps)
    ex = lam_star(ws, ps)
    e2 = lam0*(1 + kappa/2*sig**2)          # 2nd order only
    e23 = expand(sig, gam3)                  # + skewness sigma^3
    r2 = (e2-ex)/ex*100
    r23 = (e23-ex)/ex*100
    print(f"{name:24s}{float(sig):>8.4f}{float(gam3):>9.4f}{float(ex):>10.5f}"
          f"{float(e2):>10.5f}{float(e23):>10.5f}{float(r2):>7.2f}%{float(r23):>7.2f}%")

print()
print("KEY CHECKS:")
print("  * symmetric spectra (gamma3~0) at matched sigma=0.20: exact lambda_* must cluster,")
print("    2nd-order prediction identical for all (depends only on sigma) -> universality of kappa.")
print("  * asymmetric spectra: 2nd-order under/over-shoots; adding c3*gamma3*sigma^3 (SAME c3")
print("    for every shape) must pull all of them onto the exact value -> universality of c3.")
print("  * residual after 2nd+3rd is O(sigma^4): grows with sigma, differs by shape (kurtosis).")
