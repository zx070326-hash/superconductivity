"""[numerical-evidence] Independent recompute of gbar / tau_max, plus the HIGH-frequency
cutoff budget that lemmaG's direction note (c) opens up.

Written from zero by team-lead; does not read lemmaG's _chk_gbar_material_tau.py.

gbar(T) := int 2*w*a2F(w) * g(w/2piT) dw  /  int 2*w*a2F(w) dw
factor   = max(1, gbar/g2),  g2 = 1.318140212132
"""
import os; os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp
mp.mp.dps = 40

G2 = mp.mpf('1.318140212132')
PI = +mp.pi

# ---- g(varpi), two independent representations -----------------------------
def g_closed(w):
    w = mp.mpf(w)
    if w < mp.mpf('1e-4'):                       # Taylor, avoids catastrophic cancellation
        x = PI*w
        # 1 - x/sinh x = x^2/6 - 7x^4/360 + 31 x^6/15120 - ...
        return (x**2/6 - 7*x**4/360 + 31*x**6/15120 - 127*x**8/604800)/w**2
    return (1 - PI*w/mp.sinh(PI*w))/w**2

def g_series(w):                                  # alternating partial-fraction form
    w = mp.mpf(w)
    return 2*mp.nsum(lambda j: (4*j-1)/((w**2+(2*j-1)**2)*(w**2+(2*j)**2)), [1, mp.inf])

for w in ('1e-6', '0.01', '0.458879', '0.74007', '3.0', '30.0'):
    a, b = g_closed(w), g_series(w)
    print(f"  g({w:>9}) closed={mp.nstr(a,12):>16}  series={mp.nstr(b,12):>16}"
          f"  rel diff={mp.nstr(abs(a-b)/a,3)}")

wstar = mp.findroot(lambda w: g_closed(w) - G2, mp.mpf('0.46'))
print(f"\n  varpi* = {mp.nstr(wstar,10)}   omega* = {mp.nstr(2*PI*wstar,10)} k_B T")
print(f"  g(4.65/2pi) = {mp.nstr(g_closed(mp.mpf('4.65')/(2*PI)),8)}   (pointwise factor 1 iff < {mp.nstr(G2,8)})\n")

# ---- Lorentzian-difference sf spectrum  ------------------------------------
# a2F(x) ∝ 1/((x-1)^2+y^2) - 1/((x+1)^2+y^2),  x := Omega/Omega_0, y := Y/Omega_0
def gbar(y, r, xmax=mp.inf):
    """r := Omega_0/(k_B T).  varpi = Omega/(2 pi T) = r*x/(2 pi)."""
    y = mp.mpf(y); r = mp.mpf(r)
    a2F = lambda x: 1/((x-1)**2+y**2) - 1/((x+1)**2+y**2)
    w   = lambda x: 2*x*a2F(x)                    # weight 2*omega*a2F, Omega_0 units
    pts = [0, 1, 2, 5] + ([10, xmax] if xmax == mp.inf else [xmax])
    pts = sorted(set([p for p in pts if p <= xmax])) + ([xmax] if xmax not in pts else [])
    num = mp.quad(lambda x: w(x)*g_closed(r*x/(2*PI)), pts)
    den = mp.quad(w, pts)
    return num/den

print("  [A] gbar/g2 and tau_max, Lorentz difference, Omega_0 = 4.65 k_B Tc, NO cutoff")
print(f"  {'Y/Om0':>7} {'gbar/g2 @Tc':>14} {'tau_max':>10}")
for y in ('0.1', '0.25', '0.5', '1.0', '2.0'):
    rat = gbar(y, mp.mpf('4.65'))/G2
    tau = mp.findroot(lambda t: gbar(y, mp.mpf('4.65')/t) - G2, mp.mpf('3'))
    print(f"  {y:>7} {mp.nstr(rat,6):>14} {mp.nstr(tau,6):>10}")

# ---- [B] the budget lemmaG's note (c) opens: how hard a HIGH cutoff to break it?
print("\n  [B] high-frequency cutoff budget at T = Tc (Y = Om0/2, Om0 = 4.65 k_B Tc)")
print(f"  {'xmax=Om_max/Om0':>17} {'gbar/g2':>10} {'factor':>8}")
for xm in ('20', '10', '5', '3', '2', '1.5', '1.0', '0.7', '0.5', '0.3'):
    rat = gbar('0.5', mp.mpf('4.65'), mp.mpf(xm))/G2
    print(f"  {xm:>17} {mp.nstr(rat,6):>10} {'1' if rat<=1 else mp.nstr(rat,6):>8}")

# ---- [C] critical high-frequency cutoff: where does the factor first exceed 1?
print("\n  [C] critical Om_max/Om_0 at which factor leaves 1  (T = Tc, Om_0 = 4.65 k_B Tc)")
print(f"  {'Y/Om0':>7} {'xmax_crit':>11}  {'(peak of weight 2*x*a2F at x =)':>32}")
for y in ('0.1', '0.25', '0.5', '1.0', '2.0'):
    try:
        xc = mp.findroot(lambda t: gbar(y, mp.mpf('4.65'), t) - G2, mp.mpf('0.8'))
    except Exception as e:
        xc = None
    yy = mp.mpf(y)
    wfun = lambda x: 2*x*(1/((x-1)**2+yy**2) - 1/((x+1)**2+yy**2))
    grid = [mp.mpf(k)/2000 for k in range(1, 12001)]          # coarse scan, no derivative
    xpk = max(grid, key=wfun)
    print(f"  {y:>7} {mp.nstr(xc,6) if xc else 'none':>11}  {mp.nstr(xpk,6):>32}")

# ---- [D] RETRACTED.  omega_c = 6.7568*Omega_0 is NOT a cutoff on alpha^2F.
#   Verified by me in the sources: _umm2301.txt:88 Eq.(2) carries Theta(omega_c - |omega_m|)
#   -- the argument is the MATSUBARA frequency omega_m, on the Sigma_{m,j} sum, paired with
#   mu*_ij(omega_c) in the same bracket (the standard pseudopotential-cutoff convention).
#   The kernel itself is 2\int_0^{+infty} dOmega ... (_umm2301.txt:100, upper limit +infty).
#   Same in _umm1304.txt:106/:115.  So M2 = \int_0^infty 2*Omega*alpha^2F takes the FULL
#   range and the "no cutoff" row is the correct one.  [D]/[E] below are kept only as a
#   hypothetical robustness scan; their "sound" label was WRONG.
print("\n  [D] RETRACTED (omega_c is the Matsubara-sum cutoff, not a spectral cutoff)")
for xm, tag in ((mp.mpf('6.7568'), 'omega_c = 6.7568 Om0  (#3, printed)'),
                (mp.mpf('10'),     'Om_max  = 10 Om0      (#2 model)'),
                (mp.mpf('18'),     'omega_c = 18 Om0      (#1, per lemmaG)'),
                (mp.inf,           'no cutoff             (fail-open side)')):
    rat = gbar('0.5', mp.mpf('4.65'), xm)/G2
    tau = mp.findroot(lambda t: gbar('0.5', mp.mpf('4.65')/t, xm) - G2, mp.mpf('3'))
    print(f"  {tag:<36} gbar/g2 = {mp.nstr(rat,6):>9}   tau_max = {mp.nstr(tau,6):>8}")

# ---- [E] the narrowest fixed-kernel rows, WITH the printed cutoff (sound direction)
print("\n  [E] #1 LiFeAs rows (Y = Om0/2), cutoff 18*Om0 per lemmaG transcription")
print(f"  {'Om0/kBTc':>9} {'gbar/g2 no-cut':>15} {'tau (no-cut)':>13} {'gbar/g2 @18Om0':>15} {'tau (sound)':>12}")
for r in ('4.99', '4.64', '4.62'):
    r = mp.mpf(r)
    r0 = gbar('0.5', r); t0 = mp.findroot(lambda t: gbar('0.5', r/t) - G2, mp.mpf('3'))
    r1 = gbar('0.5', r, mp.mpf('18'))
    t1 = mp.findroot(lambda t: gbar('0.5', r/t, mp.mpf('18')) - G2, mp.mpf('3'))
    print(f"  {mp.nstr(r,4):>9} {mp.nstr(r0/G2,6):>15} {mp.nstr(t0,6):>13}"
          f" {mp.nstr(r1/G2,6):>15} {mp.nstr(t1,6):>12}")
