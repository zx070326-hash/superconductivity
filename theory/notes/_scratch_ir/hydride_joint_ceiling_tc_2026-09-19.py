"""
numerical-evidence, non-interval.
Direct joint-ceiling test on real Mishra H3S/D3S spectra (not two-point toy):
for each real alpha^2F(omega) shape P, force lambda := lambda_*[P] (the Semenok
stability threshold for THAT shape) and solve the N=1 KAY critical equation
lambda_*[P] * integral h1(omega/u) P(domega) = 1 for u=2*pi*Tc_ceiling.
Then compare Tc_ceiling across spectra of different width r=w2av/wlog to test
the sign of the toy-model mechanism (+1.88: wider spectrum -> higher joint ceiling)
on genuine strong-coupling DFT spectra, without requiring the material's actual
lambda to sit at the boundary.
"""
import re
import numpy as np
import mpmath as mp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

mp.mp.dps = 25
KB_MEV_PER_K = 0.08617333262

def g_mp(x):
    x = mp.mpf(x)
    ix = mp.mpc(0, 1) * x
    psi1 = mp.polygamma(1, ix)
    psi2 = mp.polygamma(2, ix)
    return float(6*x**2 + 12*x**3*mp.im(psi1) + 6*x**4*mp.re(psi2))

print("building g(x) lookup table...")
xs_log = np.linspace(np.log(1e-5), np.log(200.0), 1400)
xs = np.exp(xs_log)
gs = np.array([g_mp(x) for x in xs])
g_spline = CubicSpline(xs_log, gs, extrapolate=True)
def g_fast(x):
    x = np.clip(np.asarray(x, dtype=float), 1e-5, 200.0)
    return g_spline(np.log(x))

def h1(x):
    x = np.asarray(x, dtype=float)
    return x*x/(x*x+1.0)

def parse_epw_a2f(path):
    om, a2f = [], []
    lam_reported = None
    with open(path, errors="replace") as f:
        lines = f.readlines()
    for line in lines[1:]:
        s = line.strip()
        m = re.match(r"^Summed el-ph coupling\s+([\-0-9.Ee+]+)", s)
        if m:
            lam_reported = float(m.group(1)); continue
        parts = s.split()
        if len(parts) < 2:
            continue
        try:
            w = float(parts[0]); a = float(parts[1])
        except ValueError:
            continue
        om.append(w); a2f.append(max(a, 0.0))
    return np.array(om), np.array(a2f), lam_reported

def analyze(om, a2f, u_lo_factor=0.02, u_hi_factor=50.0, n_u=500):
    mask = om > 1e-8
    om, a2f = om[mask], a2f[mask]
    lam = 2.0 * np.trapezoid(a2f/om, om)
    w2av = np.sqrt((2.0/lam) * np.trapezoid(a2f*om, om))
    wlog = np.exp((2.0/lam) * np.trapezoid((a2f/om)*np.log(om), om))
    r = w2av/wlog
    wmin, wmax = om.min(), om.max()
    u_grid = np.exp(np.linspace(np.log(wmin*u_lo_factor), np.log(wmax*u_hi_factor), n_u))
    def F_g(u):
        return (2.0/lam) * np.trapezoid((a2f/om) * g_fast(om/u), om)
    Fvals = np.array([F_g(u) for u in u_grid])
    imax = np.argmax(Fvals)
    u_fine = np.exp(np.linspace(np.log(u_grid[max(imax-1,0)]), np.log(u_grid[min(imax+1,n_u-1)]), 300))
    Ff = np.array([F_g(u) for u in u_fine])
    imf = np.argmax(Ff)
    Fmax = Ff[imf]
    lam_star = 1.0/Fmax
    def F_h1(u):
        return (2.0/lam) * np.trapezoid((a2f/om) * h1(om/u), om)
    lo, hi = wmin*1e-3, wmax*200.0
    def resid(u):
        return lam_star * F_h1(u) - 1.0
    r_lo, r_hi = resid(lo), resid(hi)
    if r_lo*r_hi > 0:
        for _ in range(20):
            lo *= 0.5; hi *= 2.0
            r_lo, r_hi = resid(lo), resid(hi)
            if r_lo*r_hi <= 0:
                break
    u_tc = brentq(resid, lo, hi, xtol=1e-8, rtol=1e-12)
    Tc_ceiling_K = (u_tc/(2*np.pi)) / KB_MEV_PER_K
    return dict(lam=lam, w2av=w2av, wlog=wlog, r=r, lam_star=lam_star, u_tc=u_tc, Tc_ceiling_K=Tc_ceiling_K)

BASE = "D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/Data-H3S-D3S-PH-EPW-Anharmonic"
FILES = {
    "H3S-200GPa-har": f"{BASE}/H3S-200GPa/EPW_harmonic/H3S.a2f",
    "H3S-200GPa-anh": f"{BASE}/H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S-200GPa-har": f"{BASE}/D3S-200GPa/EPW_harmonic/D3S.a2f",
    "D3S-200GPa-anh": f"{BASE}/D3S-200GPa/EPW_anharmonic/D3S.a2f",
    "H3S-160GPa-anh": f"{BASE}/H3S-160GPa/EPW_anharmonic/H3S.a2f",
    "D3S-160GPa-anh": f"{BASE}/D3S-160GPa/EPW_anh/D3S.a2f",
}

print(f"\n{'case':18s} {'lambda':>8} {'r':>7} {'lambda_*':>9} {'headroom':>9} {'Tc_ceiling(K)':>14}")
res = {}
for name, path in FILES.items():
    om, a2f, lam_rep = parse_epw_a2f(path)
    r = analyze(om, a2f)
    res[name] = r
    print(f"{name:18s} {r['lam']:8.4f} {r['r']:7.4f} {r['lam_star']:9.4f} {r['lam_star']/r['lam']:9.4f} {r['Tc_ceiling_K']:14.2f}")

print("\n=== paired comparisons (same material/pressure, width changes) ===")
pairs = [
    ("H3S-200GPa-har", "H3S-200GPa-anh", "H3S 200GPa: harmonic -> anharmonic (adds width)"),
    ("D3S-200GPa-har", "D3S-200GPa-anh", "D3S 200GPa: harmonic -> anharmonic (adds width)"),
    ("H3S-200GPa-anh", "H3S-160GPa-anh", "H3S anharmonic: 200GPa -> 160GPa"),
    ("D3S-200GPa-anh", "D3S-160GPa-anh", "D3S anharmonic: 200GPa -> 160GPa"),
]
for a, b, desc in pairs:
    ra, rb = res[a], res[b]
    dr = rb['r'] - ra['r']
    dlnTc = np.log(rb['Tc_ceiling_K']/ra['Tc_ceiling_K'])
    print(f"{desc:55s} dr={dr:+.4f}  dlnTc_ceiling={dlnTc:+.4f}  same_sign={ (dr>0)==(dlnTc>0) }")
