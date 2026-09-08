r"""
P3 §4 — analytic closed form for the curvature coefficient kappa in
    lambda_*(sigma) = (1/g_max) * [1 + (kappa/2) sigma^2 + O(sigma^4)]
of the Semenok-Altshuler-Yuzbashyan stability upper bound.

Semenok Eq. S12:  lambda_*[P] = 1 / max_T \int_0^inf g(w/2piT) P(w) dw
  g(x) = 6 x^2 + 12 x^3 Im psi'(ix) + 6 x^4 Re psi''(ix),   x = w/(2 pi T)
  single-peak maximum g_max @ x_max  (literature 0.2709 @ 0.3273).

Result (derived here): kappa = - x_max^2 g''(x_max) / g_max, with the closed
form of g''(x) obtained analytically via  d/dx psi^(n)(ix) = i psi^(n+1)(ix)
and [Im/Re, d/dx] = 0:
  A=Im psi'(ix)  -> A'=Re psi''(ix)=B
  B=Re psi''(ix) -> B'=-Im psi'''(ix)=-C
  C=Im psi'''(ix)-> C'=Re psi''''(ix)=D
  g'(x)  = 12x + 36x^2 A + 36x^3 B - 6x^4 C
  g''(x) = 12  + 72x A + 144x^2 B - 60x^3 C - 6x^4 D
Validity: 2nd-order in sigma (narrow-spectrum). Tested against exact lambda_*[P]
for model spectra below; good for near-Einstein (LaH10 regime), breaks for Debye.
"""
import math
import mpmath as mp
mp.mp.dps = 25

A = lambda x: mp.im(mp.polygamma(1, 1j*x))
B = lambda x: mp.re(mp.polygamma(2, 1j*x))
C = lambda x: mp.im(mp.polygamma(3, 1j*x))
D = lambda x: mp.re(mp.polygamma(4, 1j*x))

def g(x):
    x = mp.mpf(x)
    return 6*x**2 + 12*x**3*A(x) + 6*x**4*B(x)

def gp(x):   # g'(x), hand-derived closed form
    x = mp.mpf(x)
    return 12*x + 36*x**2*A(x) + 36*x**3*B(x) - 6*x**4*C(x)

def gpp(x):  # g''(x), hand-derived closed form
    x = mp.mpf(x)
    return 12 + 72*x*A(x) + 144*x**2*B(x) - 60*x**3*C(x) - 6*x**4*D(x)

# --- cross-check hand-derived closed forms vs numerical differentiation ---
print("cross-check analytic vs numeric derivatives:")
for s in ['0.20', '0.3273', '0.55']:
    x0 = mp.mpf(s)
    e1 = abs(gp(x0)  - mp.diff(g, x0))
    e2 = abs(gpp(x0) - mp.diff(g, x0, 2))
    print(f"  x={float(x0):.4f}  |dg'|={mp.nstr(e1,3)}   |dg''|={mp.nstr(e2,3)}")

# --- x_max, g_max, kappa ---
xmax = mp.findroot(gp, mp.mpf('0.33'))
gmax = g(xmax)
gpp_max = gpp(xmax)
kappa = -xmax**2 * gpp_max / gmax
print()
print(f"x_max      = {mp.nstr(xmax,10)}")
print(f"g_max      = {mp.nstr(gmax,10)}   (lit 0.2709)")
print(f"1/g_max    = {mp.nstr(1/gmax,10)}   = lambda_*(Einstein)  (lit 3.69)")
print(f"g''(x_max) = {mp.nstr(gpp_max,10)}")
print(f"kappa      = {mp.nstr(kappa,10)}   (prior numeric est ~2.8)")
print(f"=> lambda_*(sigma) = {mp.nstr(1/gmax,6)} [1 + {mp.nstr(kappa/2,6)} sigma^2 + O(sigma^4)]")

# --- exact lambda_*[P] for model spectra vs the 2nd-order expansion ---
def lam_star_exact(ws, ps):
    pts = [(mp.mpf(w), mp.mpf(p)) for w, p in zip(ws, ps) if p > 0 and w > 0]
    tot = sum(p for _, p in pts); pts = [(w, p/tot) for w, p in pts]
    mean = sum(w*p for w, p in pts)
    def G(s):
        s = mp.mpf(s)
        return sum(p*g(s*w) for w, p in pts)
    s0 = mp.findroot(lambda s: mp.diff(G, s), mp.mpf(xmax)/mean)
    return 1/G(s0)

def sigma_of(ws, ps):
    pts = [(mp.mpf(w), mp.mpf(p)) for w, p in zip(ws, ps) if p > 0 and w > 0]
    tot = sum(p for _, p in pts); pts = [(w, p/tot) for w, p in pts]
    m1 = sum(p*w for w, p in pts)
    m2 = sum(p*w*w for w, p in pts)
    return mp.sqrt(m2-m1*m1)/m1

def twopeak(mean, split):
    return [mean*(1-split), mean*(1+split)], [1.0, 1.0]

def gauss_disc(mean, rel, N=601, span=6):
    sd = rel*mean
    ws = [mean + sd*(-span + 2*span*k/(N-1)) for k in range(N)]
    ps = [math.exp(-0.5*((w-mean)/sd)**2) if w > 0 else 0.0 for w in ws]
    return ws, ps

def debye_P(N=1500):  # normalized shape P(w) = (2/lambda) a2F/w ; Debye a2F~w^2 => P~w on (0,wD]
    ws = [(k+0.5)/N for k in range(N)]
    ps = [w for w in ws]
    return ws, ps

cases = [
    ("Einstein (delta)",   ([1.0], [1.0])),
    ("2-peak split=0.05",  twopeak(1.0, 0.05)),
    ("2-peak split=0.10",  twopeak(1.0, 0.10)),
    ("2-peak split=0.20",  twopeak(1.0, 0.20)),
    ("gauss rel=0.10",     gauss_disc(1.0, 0.10)),
    ("gauss rel=0.20",     gauss_disc(1.0, 0.20)),
    ("gauss rel=0.30",     gauss_disc(1.0, 0.30)),
    ("Debye (P~w, a2F~w^2)", debye_P()),
]
print("\nmodel spectra: exact lambda_* vs 2nd-order expansion")
print(f"{'spectrum':22s}{'sigma':>9s}{'exact':>11s}{'expand':>11s}{'rel.err':>9s}")
for name, (ws, ps) in cases:
    sig = sigma_of(ws, ps)
    ex = lam_star_exact(ws, ps)
    exp = (1/gmax)*(1 + kappa/2*sig**2)
    rel = (exp-ex)/ex
    print(f"{name:22s}{float(sig):>9.4f}{float(ex):>11.5f}{float(exp):>11.5f}{float(rel)*100:>8.2f}%")
