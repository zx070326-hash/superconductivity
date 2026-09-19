"""
Third, fully independent confirmation path for the load-bearing constants of
../misc/hydride_margin_one_parameter_point_2026-09-19.md.

Deliberately shares NOTHING with hydride_margin_one_2026-09-19.py:
  * pure mpmath (dps=40), no numpy anywhere
  * ternary search instead of golden-section search for g_max
  * the spectral moments are rebuilt from the .a2f file's OWN cumulative column
    (the third column, = lambda(omega) = int 2 a2F/w), treated as a Stieltjes
    measure -- i.e. no quadrature of alpha^2 F by me at all.

Caveat this script is meant to expose: the Stieltjes/midpoint reconstruction is
only FIRST-order accurate for the second moment on the 0.477 meV grid, so its
absolute <w^2> differs from the (second-order) trapezoid/Simpson values by ~0.3%.
It therefore confirms the lambda<w^2> harmonic/anharmonic invariance only at the
~0.03% level, NOT at the 1.3e-6 level quoted in the note -- that tighter figure
rests on the two second-order quadratures, which agree with each other to 3e-10.
"""
import mpmath as mp

mp.mp.dps = 40

BASE = ("D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/"
        "Data-H3S-D3S-PH-EPW-Anharmonic/")
G2_HI = mp.mpf("1.318140212373478")   # CERT-D upper endpoint
KB = mp.mpf("0.08617333262")


def g(x):
    x = mp.mpf(x)
    ix = mp.mpc(0, 1) * x
    return 6*x**2 + 12*x**3*mp.im(mp.polygamma(1, ix)) + 6*x**4*mp.re(mp.polygamma(2, ix))


# ---- (1) g_max by ternary search (not golden section) -----------------------
a, b = mp.mpf('0.1'), mp.mpf('1')
for _ in range(200):
    m1 = a + (b - a)/3
    m2 = b - (b - a)/3
    if g(m1) < g(m2):
        a = m1
    else:
        b = m2
x_star = (a + b)/2
g_max = g(x_star)
lam_star_E = 1/g_max
print("(1) Semenok kernel maximum")
print("    x*                  =", mp.nstr(x_star, 15))
print("    g_max               =", mp.nstr(g_max, 15))
print("    lambda_*[Einstein]  =", mp.nstr(lam_star_E, 15),
      "   (repo P3 value 3.6915)")

# ---- (2) Einstein spectrum sitting at its own stability boundary -------------
xbar = mp.sqrt((1/lam_star_E)/(1 - 1/lam_star_E))   # from h1(xbar) = 1/lambda_*
sat_E = (1/xbar)/mp.sqrt(lam_star_E)
sq = mp.sqrt(G2_HI)
print("\n(2) Einstein at lambda = lambda_*  (closed form)")
print("    varpi_bar = w0/u    =", mp.nstr(xbar, 12), "  (toy note baseline 0.609542)")
print("    u/w0                =", mp.nstr(1/xbar, 12), "  (toy note baseline 1.640577)")
print("    sqrt(g(2))          =", mp.nstr(sq, 12))
print("    C_inf               =", mp.nstr(sq/(2*mp.pi), 13),
      "  (note/repo 0.1827262477890)")
print("    s^(1)               =", mp.nstr(sat_E, 12))
print("    m^(1)(lambda_*)     =", mp.nstr(sq/sat_E, 12))

# ---- (3) moments straight off the file's own cumulative column ---------------
def stieltjes(path):
    om, cum = [], []
    for line in open(path, errors='replace'):
        p = line.split()
        if len(p) != 3:
            continue
        try:
            v = [mp.mpf(t) for t in p]
        except (ValueError, TypeError):
            continue
        om.append(v[0]); cum.append(v[2])
    lam = cum[-1]
    w2 = mp.mpf(0)
    for i in range(1, len(om)):
        dP = (cum[i] - cum[i-1]) / lam
        wm = (om[i] + om[i-1]) / 2
        w2 += wm**2 * dP
    return lam, w2


print("\n(3) moments from the file's cumulative column (Stieltjes, midpoint rule)")
print(f"    {'case':17s} {'lambda':>10} {'w2av':>10} {'lambda<w^2>':>14} {'Tc_1L(K)':>10}")
vals = {}
for nm, rel in [("H3S-200GPa-har", "H3S-200GPa/EPW_harmonic/H3S.a2f"),
                ("H3S-200GPa-anh", "H3S-200GPa/EPW_anharmonic/H3S.a2f"),
                ("D3S-200GPa-har", "D3S-200GPa/EPW_harmonic/D3S.a2f"),
                ("D3S-200GPa-anh", "D3S-200GPa/EPW_anharmonic/D3S.a2f"),
                ("H3S-160GPa-anh", "H3S-160GPa/EPW_anharmonic/H3S.a2f"),
                ("D3S-160GPa-anh", "D3S-160GPa/EPW_anh/D3S.a2f")]:
    lam, w2 = stieltjes(BASE + rel)
    vals[nm] = lam * w2
    tc = (sq/(2*mp.pi)) * mp.sqrt(lam*w2) / KB
    print(f"    {nm:17s} {mp.nstr(lam,7):>10} {mp.nstr(mp.sqrt(w2),7):>10} "
          f"{mp.nstr(lam*w2,10):>14} {mp.nstr(tc,7):>10}")

print("\n(4) harmonic -> anharmonic invariance of lambda<w^2>, this coarse route")
for a_, b_ in [("H3S-200GPa-har", "H3S-200GPa-anh"), ("D3S-200GPa-har", "D3S-200GPa-anh")]:
    rd = abs(vals[a_] - vals[b_]) / vals[a_]
    print(f"    {a_} -> {b_}: rel diff = {mp.nstr(rd,4)}")
print("    (the note's 1.3e-6 / 4.0e-7 come from the two second-order quadratures;")
print("     this first-order route can only confirm the invariance to ~3e-4.)")
