"""[numerical-evidence] Material-resolved GF off-diagonal factor on the three (K1)
s+- papers, recomputed from zero (independent of team-lead's table).

factor(T0) = max(1, gbar(T0)/g2),   g2 := g(2) = 1.318140212 (CERT-D midpoint)
gbar(T0)   = int 2w a2F(w) g(w/2piT0) dw / int 2w a2F(w) dw      [Lemma g, GF-mix]
g(v)       = (1 - pi v / sinh(pi v)) / v^2   (strictly decreasing, sup pi^2/6)

Direction: g decreasing in v = w/(2 pi T0)  =>  gbar INCREASES with T0
=> factor is evaluated at an UPPER bound T0 >= Tc, and "factor == 1 at T0"
   is inherited by every T <= T0.  tau_max := largest T0/Tc with gbar <= g2.

Spectral shapes actually printed in the papers (x = w / Omega0):
  (S1) Lorentz difference  a2F ~ 1/((x-1)^2+y^2) - 1/((x+1)^2+y^2), y = Y/Omega0
       #1 arXiv:1304.3638  : Y = Omega0/2  (y=0.5), Omega0 = 8 meV, Tc_calc 18.6-20.1 K
       #3 arXiv:2301.04557 : Y = Omega0/2  (y=0.5), Omega0 = 4.65 kB Tc  (r = 4.65)
  (S2) "standard"  a2F ~ x theta(10-x)/(x^2+1)    (arXiv:1101.3473, usual-shape rows;
       the theta(Omega_max - w), Omega_max = 10 Omega0, is part of the printed model,
       so it is NOT an extra truncation)
       Omega0 / Tc*(no feedback) from Table I / Table II (read from the PDF page image).
mpmath quad at dps=30 so that g near v->0 is not hit by float64 cancellation.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import mpmath as mp

mp.mp.dps = 30
KB = mp.mpf("0.08617333262")          # meV / K
G2 = mp.mpf("1.318140212")            # CERT-D enclosure midpoint (rational-exact)
PI = +mp.pi


def g(v):
    v = mp.mpf(v)
    if v < mp.mpf("1e-6"):
        return PI**2 / 6 - 7 * PI**4 / 360 * v * v
    return (1 - PI * v / mp.sinh(PI * v)) / (v * v)


def shape_S1(y):
    y = mp.mpf(y)
    a2F = lambda x: 1 / ((x - 1)**2 + y * y) - 1 / ((x + 1)**2 + y * y)
    return (lambda x: 2 * x * a2F(x)), [0, 1, 10, mp.inf]


def shape_S2():
    return (lambda x: 2 * x * x / (x * x + 1)), [0, 1, 10]


def gbar(shape, v0):
    """v0 = Omega0 / (2 pi T0).  Returns weighted mean of g over the spectrum."""
    w, pts = shape
    num = mp.quad(lambda x: w(x) * g(x * v0), pts)
    den = mp.quad(w, pts)
    return num / den


def tau_max(shape, r):
    """r = Omega0/(kB Tc).  Largest tau = T0/Tc with gbar(Omega0/(2 pi tau Tc)) <= G2."""
    v_c = r / (2 * PI)                       # v0 at T0 = Tc
    f = lambda v: gbar(shape, v) - G2         # decreasing in v ; root v_crit
    lo, hi = mp.mpf("1e-4"), v_c
    if f(hi) > 0:
        return None, v_c, gbar(shape, v_c)    # already > 1 at Tc itself
    if f(lo) < 0:
        return mp.inf, v_c, gbar(shape, v_c)
    v_crit = mp.findroot(f, (lo, hi), solver="bisect", tol=1e-18, maxsteps=200)
    return v_c / v_crit, v_c, gbar(shape, v_c)


print("pointwise check   v0=4.65/(2pi) =", mp.nstr(mp.mpf("4.65") / (2 * PI), 8),
      "  g(v0) =", mp.nstr(g(mp.mpf("4.65") / (2 * PI)), 8), "  (team-lead: 0.98761)")
print("threshold v* = g^-1(G2) =", mp.nstr(mp.findroot(lambda v: g(v) - G2, 0.46), 8))
print()

print("(S1) Lorentz difference, r = 4.65  (team-lead's row, sanity)")
print("  y      gbar/G2    tau_max")
for y in ("0.1", "0.25", "0.5", "1.0", "2.0"):
    tm, vc, gb = tau_max(shape_S1(y), mp.mpf("4.65"))
    print(f"  {y:4s}   {mp.nstr(gb / G2, 4):>7s}   {mp.nstr(tm, 4) if tm else 'FAIL@Tc'}")
print()

print("=== per paper, printed Y and printed Omega0 / Tc ===")
rows = [
    # label, shape, Omega0 [meV], Tc used [K]  (model's own Tc, fixed kernel)
    ("#1 LiFeAs  Tc_calc=18.6 (row3)", shape_S1("0.5"), 8.0, 18.6),
    ("#1 LiFeAs  Tc_calc=20.0 (row4)", shape_S1("0.5"), 8.0, 20.0),
    ("#1 LiFeAs  Tc_calc=20.1 (row5)", shape_S1("0.5"), 8.0, 20.1),
    ("#1 LiFeAs  Tc_exp =18.0        ", shape_S1("0.5"), 8.0, 18.0),
    ("#2 BaFeCoAs usual  Tc*=33.00", shape_S2(), 9.04, 33.00),
    ("#2 LaFeAsOF usual  Tc*=37.22", shape_S2(), 11.44, 37.22),
    ("#2 BaKFeAs  usual  Tc*=55.26", shape_S2(), 14.80, 55.26),
    ("#2 SmFeAsOF usual  Tc*=74.13", shape_S2(), 20.80, 74.13),
]
print("  paper / row                      r=Om0/kTc   v0@Tc    gbar/G2@Tc   tau_max")
for lab, sh, om, tc in rows:
    r = mp.mpf(om) / (KB * mp.mpf(tc))
    tm, vc, gb = tau_max(sh, r)
    print(f"  {lab:32s} {mp.nstr(r, 5):>8s}  {mp.nstr(vc, 5):>7s}   {mp.nstr(gb / G2, 4):>8s}"
          f"   {mp.nstr(tm, 4) if tm else 'FAIL@Tc'}")

print()
print("#3 arXiv:2301.04557: r = 4.65 by construction, y = 0.5  ->", end=" ")
tm, vc, gb = tau_max(shape_S1("0.5"), mp.mpf("4.65"))
print("gbar/G2 =", mp.nstr(gb / G2, 5), " tau_max =", mp.nstr(tm, 5))

print()
print("#2 Lorentz rows: half-width Y NOT printed on the page I read -> scan y")
print("  y     BaFeCoAs(28.95) LaFeAsOF(31.81) BaKFeAs(46.18) SmFeAsOF(63.82)")
for y in ("0.25", "0.5", "1.0"):
    out = []
    for om, tc in ((9.04, 28.95), (11.44, 31.81), (14.80, 46.18), (20.80, 63.82)):
        r = mp.mpf(om) / (KB * mp.mpf(tc))
        tm, vc, gb = tau_max(shape_S1(y), r)
        out.append(mp.nstr(tm, 4) if tm else "FAIL")
    print(f"  {y:4s}  " + "  ".join(f"{o:>14s}" for o in out))
