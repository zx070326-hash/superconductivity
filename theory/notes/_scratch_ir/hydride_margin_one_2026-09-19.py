"""
numerical-evidence level (float64/mpmath point evaluation, NOT interval arithmetic),
EXCEPT the two analytic inequalities that this script only *checks* numerically
(they are proved on paper in the accompanying note):
    (A1)  s^(1)(lambda,P) := u_Tc / sqrt(lambda*<w^2>)  <=  1       for all lambda>0, all P
    (A2)  lambda_*[P] >= lambda_*[Einstein] = 1/g_max               for all P

Independent re-derivation of hydride_lambda_star_joint_ceiling_2026-09-19.md, and
extension: locate the "margin = 1" parameter point for the two distinct ceilings.

INDEPENDENCE from the two prior scripts (CLAUDE.md: load-bearing numbers must be
recomputed from scratch, not inherited):
  * own .a2f parser (requires exactly 3 float fields; no shared code)
  * g(x) computed TWO independent ways: mpmath polygamma, and a directly summed
    Hurwitz-type series derived by hand here (see g_series); cross-checked
  * quadrature done TWO ways (trapezoid and Simpson) and compared
  * NO lookup-table spline for g on the load-bearing path (the prior scripts used a
    1400-point cubic spline in log x; here g is evaluated exactly at every grid point
    via mpmath for the load-bearing rows, spline only for the wide lambda scans)
"""
import os
import re
import numpy as np
import mpmath as mp

os.environ.setdefault("OMP_NUM_THREADS", "1")
mp.mp.dps = 30

# ---------------------------------------------------------------- constants
# C_inf = sqrt(g(2))/(2 pi), CERT-D upper endpoint (oneliner_zero_parameter_2026-09-18.md)
G2_LO = mp.mpf("1.318140211762671")
G2_HI = mp.mpf("1.318140212373478")
SQRT_G2_HI = float(mp.sqrt(G2_HI))
SQRT_G2_LO = float(mp.sqrt(G2_LO))
C_INF_HI = float(mp.sqrt(G2_HI) / (2 * mp.pi))
KB_MEV_PER_K = 0.08617333262


# ---------------------------------------------------------------- Semenok kernel g(x)
def g_polygamma(x):
    """g(x) = 6x^2 + 12 x^3 Im psi'(ix) + 6 x^4 Re psi''(ix)   (Semenok Eq. S12)."""
    x = mp.mpf(x)
    ix = mp.mpc(0, 1) * x
    return float(6 * x**2 + 12 * x**3 * mp.im(mp.polygamma(1, ix))
                 + 6 * x**4 * mp.re(mp.polygamma(2, ix)))


def g_series(x, nmax=2000000):
    """
    Independent evaluation #2, via the defining series, derived by hand here:
        psi'(ix)  = sum_{n>=0} 1/(n+ix)^2    => Im psi'(ix)  = -sum_{n>=1} 2 n x/(n^2+x^2)^2
        psi''(ix) = -2 sum_{n>=0} 1/(n+ix)^3 => Re psi''(ix) = -2 sum_{n>=1}(n^3-3n x^2)/(n^2+x^2)^3
    (the n=0 terms contribute nothing to Im psi' / Re psi'' for purely imaginary argument)
        g(x) = 6x^2 - 24 x^4 S1 - 12 x^4 S2,
        S1 = sum n/(n^2+x^2)^2 ,  S2 = sum (n^3-3 n x^2)/(n^2+x^2)^3
    Tail beyond nmax added by its leading integral (Euler-Maclaurin order 0).
    Slow; used ONLY for spot-checks, never on the load-bearing path.
    """
    x = float(x)
    n = np.arange(1, nmax + 1, dtype=float)
    d = n * n + x * x
    S1 = np.sum(n / d**2)
    S2 = np.sum((n**3 - 3 * n * x * x) / d**3)
    N = float(nmax)
    S1 += 0.5 / (N * N + x * x)
    S2 += 0.5 * (N * N - 2 * x * x) / (N * N + x * x) ** 2
    return 6 * x * x - 24 * x**4 * S1 - 12 * x**4 * S2


# --- evaluation #3: vectorised complex polygamma (recurrence + Bernoulli asymptotics).
#     This is the load-bearing path: fast, no lookup table, no spline, and it is
#     cross-checked against BOTH mpmath and the hand-derived series above.
_K_SHIFT = 40


def _psi1_psi2(z):
    """trigamma psi'(z) and tetragamma psi''(z) for complex z, Re z >= 0, vectorised."""
    z = np.asarray(z, dtype=complex)
    s1 = np.zeros_like(z)
    s2 = np.zeros_like(z)
    for k in range(_K_SHIFT):
        zk = z + k
        s1 += 1.0 / zk**2
        s2 += -2.0 / zk**3
    w = z + _K_SHIFT
    iw = 1.0 / w
    # psi'(w) ~ 1/w + 1/(2w^2) + 1/(6w^3) - 1/(30w^5) + 1/(42w^7) - 1/(30w^9) + 5/(66w^11)
    p1 = (iw + 0.5 * iw**2 + iw**3 / 6.0 - iw**5 / 30.0 + iw**7 / 42.0
          - iw**9 / 30.0 + 5.0 * iw**11 / 66.0)
    # psi''(w) = d/dw psi'(w)
    p2 = (-iw**2 - iw**3 - 0.5 * iw**4 + iw**6 / 6.0 - iw**8 / 6.0
          + 0.3 * iw**10 - 5.0 * iw**12 / 6.0)
    return s1 + p1, s2 + p2


def g_vec(x):
    """g(x) for an array of x > 0, vectorised."""
    x = np.asarray(x, dtype=float)
    z = 1j * x
    p1, p2 = _psi1_psi2(z)
    return 6 * x**2 + 12 * x**3 * p1.imag + 6 * x**4 * p2.real


def g_fast(x):
    """scalar wrapper"""
    return float(g_vec(np.atleast_1d(x))[0])


def h1(x):
    x = np.asarray(x, dtype=float)
    return x * x / (x * x + 1.0)


# ---------------------------------------------------------------- g_max  (Einstein lambda_*)
def golden_max(f, a, b, tol=1e-14, it=300):
    invphi = (np.sqrt(5.0) - 1.0) / 2.0
    c, d = b - invphi * (b - a), a + invphi * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(it):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - invphi * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + invphi * (b - a)
            fd = f(d)
        if abs(b - a) < tol * max(1.0, abs(a)):
            break
    return 0.5 * (a + b), f(0.5 * (a + b))


def scan_then_golden(f, lo, hi, n=4000):
    """Coarse log-grid scan to bracket the global max, then golden refine.
    (A bare golden search is unsafe here: g(x) is positive only on a small
    interval near x~0.33 and turns negative above it, so a wrong bracket
    silently returns an endpoint -- which is exactly what happened on the
    first run of this script.)"""
    xs = np.exp(np.linspace(np.log(lo), np.log(hi), n))
    vs = np.array([f(x) for x in xs])
    i = int(np.argmax(vs))
    if i in (0, n - 1):
        raise RuntimeError(f"max at scan edge: i={i}, x={xs[i]}")
    return golden_max(f, xs[i - 1], xs[i + 1], tol=1e-15, it=300)
print("=" * 78)
print(f"{'x':>10} {'g_polygamma(mpmath)':>22} {'g_series(hand)':>18} {'g_vec(load-bearing)':>21}"
      f" {'max rel diff':>13}")
for xt in [0.05, 0.1, 0.3, 0.60954, 1.0, 2.0, 5.0, 20.0]:
    a, b, c = g_polygamma(xt), g_series(xt), g_fast(xt)
    rd = max(abs(a - b), abs(a - c)) / abs(a)
    print(f"{xt:10.5f} {a:22.14f} {b:18.12f} {c:21.14f} {rd:13.3e}")

x_gmax, g_max = scan_then_golden(g_fast, 1e-3, 5.0)
lam_star_einstein = 1.0 / g_max
_, g_max_mp = scan_then_golden(g_polygamma, 0.05, 2.0, n=400)
print(f"\n  argmax x*          = {x_gmax:.12f}")
print(f"  g_max (g_vec)      = {g_max:.12f}  -> lambda_*[Einstein] = 1/g_max = "
      f"{lam_star_einstein:.6f}")
print(f"  g_max (mpmath)     = {g_max_mp:.12f}  -> 1/g_max = {1.0/g_max_mp:.6f}")
print(f"  repo value (P3_tc_max_pressure.md; toy note s->0 gave 3.691492):")
print(f"    mine {lam_star_einstein:.6f}   rel dev from 3.691492 = "
      f"{abs(lam_star_einstein-3.691492)/3.691492:.2e}")


# ---------------------------------------------------------------- data
def parse_a2f(path):
    """Independent parser: keep only lines with exactly 3 parseable floats."""
    om, a2f, cum, lam_sum, lam_int = [], [], [], None, None
    prev_int_tag = False
    with open(path, errors="replace") as fh:
        for line in fh:
            m = re.search(r"Summed el-ph coupling\s+([-\d.Ee+]+)", line)
            if m:
                lam_sum = float(m.group(1)); prev_int_tag = False; continue
            if "Integrated el-ph coupling" in line:
                prev_int_tag = True; continue
            p = line.split()
            if prev_int_tag and len(p) == 2 and p[0] == "#":
                lam_int = float(p[1]); prev_int_tag = False; continue
            prev_int_tag = False
            if len(p) != 3:
                continue
            try:
                vals = [float(t) for t in p]
            except ValueError:
                continue
            om.append(vals[0]); a2f.append(vals[1]); cum.append(vals[2])
    om = np.asarray(om); a2f = np.asarray(a2f); cum = np.asarray(cum)
    order = np.argsort(om)
    return om[order], np.maximum(a2f[order], 0.0), cum[order], lam_sum, lam_int


def simpson_nonuniform(y, x):
    """Composite Simpson for UNEQUAL intervals (Cartwright): parabola through each
    consecutive triple.  Independent of np.trapezoid.  Needed because the
    H3S-200GPa harmonic .a2f carries one irregular trailing panel (a 12 meV jump
    after a uniform 0.4769 meV grid) -- a uniform-h Simpson is simply invalid there."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    n = len(x) - 1
    total = 0.0
    i = 0
    while i + 2 <= n:
        h0 = x[i+1] - x[i]; h1 = x[i+2] - x[i+1]
        total += (h0 + h1) / 6.0 * ((2 - h1/h0) * y[i]
                                    + (h0 + h1)**2 / (h0*h1) * y[i+1]
                                    + (2 - h0/h1) * y[i+2])
        i += 2
    if i < n:  # odd number of panels: trapezoid the last one
        total += 0.5 * (x[n] - x[n-1]) * (y[n] + y[n-1])
    return total


BASE = ("D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/"
        "Data-H3S-D3S-PH-EPW-Anharmonic")
FILES = {
    "H3S-200GPa-har": f"{BASE}/H3S-200GPa/EPW_harmonic/H3S.a2f",
    "H3S-200GPa-anh": f"{BASE}/H3S-200GPa/EPW_anharmonic/H3S.a2f",
    "D3S-200GPa-har": f"{BASE}/D3S-200GPa/EPW_harmonic/D3S.a2f",
    "D3S-200GPa-anh": f"{BASE}/D3S-200GPa/EPW_anharmonic/D3S.a2f",
    "H3S-160GPa-anh": f"{BASE}/H3S-160GPa/EPW_anharmonic/H3S.a2f",
    "D3S-160GPa-anh": f"{BASE}/D3S-160GPa/EPW_anh/D3S.a2f",
}


class Spec:
    """A fixed spectral SHAPE P(dw) = 2 a2F(w) dw / (lambda w), with lambda freed."""

    def __init__(self, name, path):
        om, a2f, cum, lam_sum, lam_int = parse_a2f(path)
        m = om > 1e-8
        self.name, self.om, self.a2f = name, om[m], a2f[m]
        self.lam_sum, self.lam_int = lam_sum, lam_int
        w, a = self.om, self.a2f
        self.lam_t = 2.0 * np.trapezoid(a / w, w)
        self.lam_s = 2.0 * simpson_nonuniform(a / w, w)
        self.lam = self.lam_t
        # normalised shape weights: P(dw) = (2/lam) a2F/w dw  (integrates to 1)
        self.w = w
        self.dens = (2.0 / self.lam) * (a / w)          # dP/dw
        self.norm_t = np.trapezoid(self.dens, w)
        self.w2 = np.trapezoid(self.dens * w * w, w)     # <w^2>
        self.w2_s = simpson_nonuniform(self.dens * w * w, w)
        self.wlog = np.exp(np.trapezoid(self.dens * np.log(w), w))
        self.w2av = np.sqrt(self.w2)
        self.r = self.w2av / self.wlog

    def P_of(self, f_of_x, u):
        """integral f(w/u) P(dw)"""
        return np.trapezoid(self.dens * f_of_x(self.w / u), self.w)

    # ---- Semenok stability threshold, exact g at every node (no spline, no table)
    def lambda_star(self, n_u=400, refine=3):
        def G(u):
            return np.trapezoid(self.dens * g_vec(self.w / u), self.w)

        lo, hi = np.log(self.w.min() * 0.02), np.log(self.w.max() * 60.0)
        edge = False
        for _ in range(refine):
            us = np.exp(np.linspace(lo, hi, n_u))
            vals = np.array([G(u) for u in us])
            i = int(np.argmax(vals))
            edge = (i == 0 or i == n_u - 1)
            lo, hi = np.log(us[max(i - 1, 0)]), np.log(us[min(i + 1, n_u - 1)])
            n_u = 60
        u_s, Gmax = golden_max(G, np.exp(lo), np.exp(hi), tol=1e-13, it=200)
        return 1.0 / Gmax, u_s, edge

    # ---- N=1 KAY criticality:  lambda * int h1(w/u) P(dw) = 1
    #      Note int h1 -> 1 as u->0 and -> 0 as u->infinity, so a root exists
    #      IFF lambda > 1.  (h1 < 1 pointwise => int h1 P < 1 for every u.)
    def u_tc(self, lam):
        if lam <= 1.0:
            return np.nan
        f = lambda u: lam * self.P_of(h1, u) - 1.0
        lo = self.w.min() * 1e-6
        hi = self.w.max() * max(1e3, 10.0 * np.sqrt(lam))
        for _ in range(80):
            if f(hi) < 0:
                break
            hi *= 4.0
        for _ in range(80):
            if f(lo) > 0:
                break
            lo *= 0.25
        if f(lo) * f(hi) > 0:
            return np.nan
        for _ in range(300):
            mid = np.sqrt(lo * hi)
            if f(mid) > 0:
                lo = mid
            else:
                hi = mid
            if hi / lo - 1 < 1e-15:
                break
        return np.sqrt(lo * hi)

    def tc_K(self, lam):
        return self.u_tc(lam) / (2 * np.pi) / KB_MEV_PER_K

    def sat(self, lam):
        """s^(1) = 2 pi Tc / sqrt(lambda <w^2>) = u_Tc / sqrt(lambda <w^2>)."""
        return self.u_tc(lam) / np.sqrt(lam * self.w2)

    def margin_oneliner(self, lam):
        """C_inf sqrt(lambda<w^2>) / Tc^(1)(lambda) = sqrt(g(2)) / s^(1)."""
        return SQRT_G2_HI / self.sat(lam)


print("\n" + "=" * 78)
print("STEP 1 -- independent recomputation of the 6 Mishra spectra")
print("=" * 78)
specs = {}
print(f"{'case':17s} {'lam(trap)':>10} {'lam(simp)':>10} {'lam:Integ':>10} {'lam:Summed':>11} "
      f"{'relerr(Int)':>11} {'norm':>8} {'w2av':>9} {'wlog':>9} {'r':>7}")
for nm, p in FILES.items():
    s = Spec(nm, p)
    specs[nm] = s
    print(f"{nm:17s} {s.lam_t:10.5f} {s.lam_s:10.5f} {s.lam_int:10.5f} {s.lam_sum:11.5f} "
          f"{abs(s.lam_t-s.lam_int)/s.lam_int:11.4%} {s.norm_t:8.5f} "
          f"{s.w2av:9.3f} {s.wlog:9.3f} {s.r:7.4f}")
print("  (trap vs simpson = my own two quadratures; 'Integ' is the file's own cumulative")
print("   column endpoint = the same integral, 'Summed' is EPW's mode sum, a different object.)")
print(f"  max |trap-simp|/trap over the 6 = "
      f"{max(abs(s.lam_t-s.lam_s)/s.lam_t for s in specs.values()):.3e}")

print("\nlambda_*[P]  (exact g at every quadrature node; no spline):")
print(f"{'case':17s} {'lambda_*':>10} {'u_*(meV)':>10} {'edge?':>6} {'headroom=l*/l':>14} "
      f"{'>= 1/g_max?':>12}")
lstar = {}
for nm, s in specs.items():
    ls, us, edge = s.lambda_star()
    lstar[nm] = ls
    print(f"{nm:17s} {ls:10.5f} {us:10.3f} {str(edge):>6} {ls/s.lam:14.4f} "
          f"{str(ls >= lam_star_einstein):>12}")


print("\n" + "=" * 78)
print("STEP 2 -- the two DIFFERENT margins (the prior note reported only the first)")
print("=" * 78)
print("  margin_lambda  = lambda_*[P] / lambda          (a COUPLING ratio)")
print("  margin_joint   = Tc^(1)(lambda_*)/Tc^(1)(lambda)  (a Tc ratio, joint ceiling)")
print("  margin_1L      = C_inf sqrt(lambda<w^2>) / Tc^(1)(lambda)   (one-liner ceiling)")
print()
print(f"{'case':17s} {'lambda':>8} {'Tc^(1)(K)':>10} {'Tc_joint(K)':>12} {'Tc_1L(K)':>10} "
      f"{'m_lambda':>9} {'m_joint':>9} {'m_1L':>8} {'s^(1)':>8}")
rows = {}
for nm, s in specs.items():
    lam, ls = s.lam, lstar[nm]
    tc = s.tc_K(lam)
    tcj = s.tc_K(ls)
    tc1l = C_INF_HI * np.sqrt(lam * s.w2) / KB_MEV_PER_K
    rows[nm] = dict(lam=lam, lstar=ls, tc=tc, tcj=tcj, tc1l=tc1l,
                    m_lam=ls / lam, m_joint=tcj / tc, m_1l=tc1l / tc, sat=s.sat(lam))
    print(f"{nm:17s} {lam:8.4f} {tc:10.2f} {tcj:12.2f} {tc1l:10.2f} "
          f"{ls/lam:9.3f} {tcj/tc:9.3f} {tc1l/tc:8.3f} {s.sat(lam):8.5f}")

print(f"\n  sqrt(g(2)) = {SQRT_G2_HI:.9f}   <-- claimed universal floor for m_1L")
print(f"  min m_1L over the 6 spectra at their actual lambda: "
      f"{min(r['m_1l'] for r in rows.values()):.4f}")


print("\n" + "=" * 78)
print("STEP 3 -- margin = 1 parameter point, ceiling by ceiling")
print("=" * 78)

print("\n(3a) JOINT (stability) ceiling: margin_joint = 1  <=>  lambda = lambda_*[P] identically")
print("     (lambda_*[P] depends only on the SHAPE, so scaling lambda leaves it fixed;")
print("      this is an identity with ZERO degrees of freedom -- its only content is the")
print("      NUMBER lambda_*[P] and the Tc it implies.)")
print(f"{'case':17s} {'lambda now':>10} {'lambda needed':>13} {'x factor':>9} "
      f"{'Tc now(K)':>10} {'Tc at margin=1':>15} {'x factor':>9}")
for nm, r in rows.items():
    print(f"{nm:17s} {r['lam']:10.4f} {r['lstar']:13.4f} {r['lstar']/r['lam']:9.3f} "
          f"{r['tc']:10.2f} {r['tcj']:15.2f} {r['tcj']/r['tc']:9.3f}")

print("\n(3b) ONE-LINER ceiling: scan lambda at FIXED shape, look for margin_1L -> 1")
lam_scan = [1.0, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 1e3, 1e4, 1e6, 1e9]
print(f"{'case':17s} " + " ".join(f"{('L=%g' % L):>9}" for L in lam_scan))
for nm, s in specs.items():
    print(f"{nm:17s} " + " ".join(f"{s.margin_oneliner(L):9.5f}" for L in lam_scan))
print(f"{'(floor)':17s} " + " ".join(f"{SQRT_G2_HI:9.5f}" for _ in lam_scan))

print("\n     -> margin_1L decreases in lambda and converges to sqrt(g(2)) from ABOVE.")
print("        Monotonicity check on a dense log grid lambda in (1, 1e10]:")
for nm, s in specs.items():
    Ls = np.exp(np.linspace(np.log(1.0001), np.log(1e10), 400))
    ms = np.array([s.margin_oneliner(L) for L in Ls])
    ok = np.all(np.diff(ms) <= 1e-12)
    print(f"        {nm:17s} monotone_decreasing={ok}   m(1.0001)={ms[0]:.4f}"
          f"   m(1e10)={ms[-1]:.9f}   min={ms.min():.9f}")
print(f"        sqrt(g(2)) lower endpoint (CERT-D low) = {SQRT_G2_LO:.9f}")
print("        -> the infimum equals sqrt(g(2)) and is NOT attained at any finite lambda.")

print("\n(3c) check of the analytic inequality  s^(1)(lambda,P) <= 1  (=> m_1L >= sqrt(g(2)))")
viol = 0; tot = 0; worst = 0.0
for nm, s in specs.items():
    for L in [0.01, 0.1, 0.5, 1, 2, 5, 20, 100, 1e4, 1e8]:
        v = s.sat(L); tot += 1
        worst = max(worst, v)
        if v > 1.0 + 1e-12:
            viol += 1
            print(f"    VIOLATION {nm} lambda={L}: s={v}")
print(f"    {tot-viol}/{tot} satisfy s^(1) <= 1 ; sup observed s^(1) = {worst:.10f}")

print("\n(3d) shape deformation at FIXED lambda: can a different SHAPE close margin_joint?")
print("     margin_joint = 1 needs lambda >= lambda_*[P], and lambda_*[P] >= 1/g_max for ALL P.")
print(f"     1/g_max = lambda_*[Einstein] = {lam_star_einstein:.6f}")
print(f"{'case':17s} {'lambda':>8} {'best-case headroom = (1/g_max)/lambda':>40}")
for nm, s in specs.items():
    print(f"{nm:17s} {s.lam:8.4f} {lam_star_einstein/s.lam:40.4f}")

print("\n(3e) numerical stress test of  lambda_*[P] >= 1/g_max  on synthetic shapes")
rng = np.random.default_rng(20260919)
wgrid = np.linspace(1.0, 400.0, 1200)


class SynthSpec(Spec):
    def __init__(self, w, dens):
        self.name = "synth"
        self.w = w
        self.dens = dens / np.trapezoid(dens, w)
        self.lam = 1.0
        self.w2 = np.trapezoid(self.dens * w * w, w)
        self.w2av = np.sqrt(self.w2)
        self.wlog = np.exp(np.trapezoid(self.dens * np.log(w), w))
        self.r = self.w2av / self.wlog


nviol = 0
mins = []
for k in range(40):
    npk = rng.integers(1, 5)
    dens = np.zeros_like(wgrid)
    for _ in range(npk):
        c = rng.uniform(5, 380); s_ = rng.uniform(0.7, 60.0); a_ = rng.uniform(0.1, 1.0)
        dens += a_ * np.exp(-0.5 * ((wgrid - c) / s_) ** 2)
    ss = SynthSpec(wgrid, dens + 1e-12)
    ls, _, _ = ss.lambda_star(n_u=200, refine=2)
    mins.append(ls)
    if ls < lam_star_einstein - 1e-9:
        nviol += 1
        print(f"    VIOLATION synth #{k}: lambda_* = {ls:.8f} < {lam_star_einstein:.8f}")
# deliberately near-Einstein (narrow) case: should approach the floor from above
for sig in [20.0, 5.0, 2.0, 1.0, 0.5]:
    dens = np.exp(-0.5 * ((wgrid - 200.0) / sig) ** 2) + 1e-300
    ss = SynthSpec(wgrid, dens)
    ls, _, _ = ss.lambda_star(n_u=200, refine=2)
    print(f"    near-Einstein sigma={sig:5.2f} meV (r={ss.r:.5f}):  lambda_* = {ls:.6f}"
          f"   excess over 1/g_max = {ls-lam_star_einstein:+.3e}")
print(f"    {40-nviol}/40 random multi-peak shapes satisfy lambda_*[P] >= 1/g_max ; "
      f"min over the 40 = {min(mins):.6f}")


print("\n" + "=" * 78)
print("STEP 4 -- decomposition of the one-liner margin at the STABILITY boundary")
print("=" * 78)
print("  m_1L(lambda_*) = sqrt(g(2)) / s^(1)(lambda_*,P)")
print("                 = [ sqrt(g(2)) ]  x  [ 1/s^(1)(lambda_*,P) ]")
print("        factor 1 = N=1-vs-full-operator price (irreducible, analytic)")
print("        factor 2 = finite-lambda price of being capped at lambda_* (< infinity)")
print(f"{'case':17s} {'lambda_*':>9} {'s^(1)(l_*)':>11} {'m_1L(l_*)':>10} "
      f"{'factor1':>8} {'factor2':>8}")
for nm, s in specs.items():
    ls = lstar[nm]
    sa = s.sat(ls)
    print(f"{nm:17s} {ls:9.4f} {sa:11.6f} {SQRT_G2_HI/sa:10.4f} "
          f"{SQRT_G2_HI:8.4f} {1.0/sa:8.4f}")

# Einstein reference at its own stability boundary -- fully closed form
xbar = np.sqrt((1.0 / lam_star_einstein) / (1.0 - 1.0 / lam_star_einstein))
sat_E = (1.0 / xbar) / np.sqrt(lam_star_einstein)
print(f"\n  Einstein reference (closed form): h1(xbar) = 1/lambda_* = {1/lam_star_einstein:.6f}")
print(f"    xbar = w0/u = {xbar:.6f}   (toy note baseline varpi_0 = 0.609542)")
print(f"    u/w0 = {1/xbar:.6f}         (toy note baseline u_0     = 1.640577)")
print(f"    s^(1) = {sat_E:.6f}   m_1L = {SQRT_G2_HI/sat_E:.6f}")


print("\n" + "=" * 78)
print("STEP 5 -- comparison against EXPERIMENTAL Tc (what the margin really is)")
print("=" * 78)
EXP = {"H3S-200GPa-har": 184.04, "H3S-200GPa-anh": 184.04, "H3S-160GPa-anh": 203.0,
       "D3S-200GPa-har": 150.0, "D3S-200GPa-anh": 150.0, "D3S-160GPa-anh": 150.0}
print("  (H3S 200 GPa annealed 184.04 K from the repo's own vector digitisation of")
print("   Drozdov 2015 Fig. 2c; H3S 160 GPa ~ 203 K peak; D3S ~ 150 K. See note for caveats.)")
print(f"{'case':17s} {'Tc_exp(K)':>10} {'Tc^(1)(K)':>10} {'Tc_1L(K)':>10} "
      f"{'frak f':>8} {'m_1L/exp':>9} {'Tc^(1)/exp':>11}")
for nm, s in specs.items():
    r = rows[nm]; te = EXP[nm]
    print(f"{nm:17s} {te:10.2f} {r['tc']:10.2f} {r['tc1l']:10.2f} "
          f"{te/r['tc1l']:8.4f} {r['tc1l']/te:9.3f} {r['tc']/te:11.3f}")

print("\nDONE (steps 0-5).")


print("\n" + "=" * 78)
print("STEP 6 -- DIRECTION AUDIT + the bracket that is actually sound")
print("=" * 78)
print("""
  KAY: frak_k^(1) <= frak_k  =>  Tc^(1) <= Tc_true.  Therefore

      margin_true = C_inf sqrt(l<w^2>)/Tc_true  <=  C_inf sqrt(l<w^2>)/Tc^(1) = margin^(1).

  margin^(1) is an UPPER estimate of the true margin.  Consequences for what may
  and may not be claimed (one-sided-bound direction rule):

   ALLOWED : "margin_true <= margin^(1)"  -- an upper bound on how far the ceiling
             sits above the truth.  This is the direction that produces a USEFUL
             statement here, because it says the ceiling is CLOSE, not far.
   FORBIDDEN: "margin_true >= sqrt(g(2)) because s^(1)<=1".  The inequality
             s^(1) <= 1 constrains only the N=1 LOWER bound, so it says the
             framework's own lower bound can never approach the ceiling.  It says
             NOTHING about whether the true Tc can.  The only sound lower bound on
             margin_true is the one-liner theorem itself: margin_true >= 1.
""")
print(f"{'case':17s} {'margin_true bracket at lambda_actual':>38} "
      f"{'margin_true bracket at lambda_*':>34}")
for nm, s in specs.items():
    ls = lstar[nm]
    print(f"{nm:17s} {('[ 1 , %.4f ]' % s.margin_oneliner(s.lam)):>38} "
          f"{('[ 1 , %.4f ]' % s.margin_oneliner(ls)):>34}")
xbar_E = np.sqrt((1.0/lam_star_einstein)/(1.0 - 1.0/lam_star_einstein))
m_E = SQRT_G2_HI / ((1.0/xbar_E)/np.sqrt(lam_star_einstein))
print(f"{'EINSTEIN@l_*':17s} {'':>38} {('[ 1 , %.4f ]' % m_E):>34}")

print("\n(6a) is s^(1)(lambda_*[P],P) really nearly shape-independent? stress test")
print(f"{'shape':38s} {'r':>8} {'lambda_*':>9} {'s^(1)(l_*)':>11} {'m^(1)(l_*)':>11}")
tests = []
tests.append(("Einstein (closed form)", 1.0, lam_star_einstein,
              (1.0/xbar_E)/np.sqrt(lam_star_einstein)))
for sig, lo_, hi_ in [(40.0, 1.0, 400.0), (90.0, 1.0, 400.0)]:
    dens = np.exp(-0.5*((wgrid-200.0)/sig)**2)
    ss = SynthSpec(wgrid, dens); ls, _, _ = ss.lambda_star(n_u=300, refine=2)
    tests.append((f"gaussian sigma={sig:.0f} meV", ss.r, ls, ss.sat(ls)))
for name, dens in [
        ("flat on [1,400]", np.ones_like(wgrid)),
        ("1/w on [1,400] (very broad)", 1.0/wgrid),
        ("bimodal 30 & 350 meV", np.exp(-0.5*((wgrid-30.)/6)**2)+np.exp(-0.5*((wgrid-350.)/6)**2)),
        ("bimodal 5 & 390 meV (extreme)",
         np.exp(-0.5*((wgrid-5.)/2)**2)+np.exp(-0.5*((wgrid-390.)/2)**2)),
        ("90/10 split 20 & 300 meV",
         9*np.exp(-0.5*((wgrid-20.)/3)**2)+np.exp(-0.5*((wgrid-300.)/3)**2))]:
    ss = SynthSpec(wgrid, dens + 1e-300)
    ls, _, _ = ss.lambda_star(n_u=300, refine=2)
    tests.append((name, ss.r, ls, ss.sat(ls)))
for nm, s in specs.items():
    tests.append((f"REAL {nm}", s.r, lstar[nm], s.sat(lstar[nm])))
for name, rr, ls, sa in tests:
    print(f"{name:38s} {rr:8.4f} {ls:9.4f} {sa:11.6f} {SQRT_G2_HI/sa:11.5f}")
sats = [t[3] for t in tests]
print(f"\n    spread of s^(1)(lambda_*) over ALL shapes above: "
      f"min={min(sats):.6f}  max={max(sats):.6f}  "
      f"(max-min)/min={(max(sats)-min(sats))/min(sats):.3%}")
print(f"    => m^(1)(lambda_*) in [{SQRT_G2_HI/max(sats):.4f}, {SQRT_G2_HI/min(sats):.4f}]")

print("\n(6b) three-factor decomposition of the observed margin (at fixed shape)")
print("     m^(1)(lambda_act) = [sqrt(g(2))] x [1/s^(1)(l_*)] x [s^(1)(l_*)/s^(1)(l_act)]")
print("                        operator      finite-lambda      distance from the")
print("                        price         at the boundary    stability boundary")
print(f"{'case':17s} {'m^(1)(l_act)':>13} {'F1 operator':>12} {'F2 finite-l':>12} "
      f"{'F3 dist':>9} {'product':>9}")
for nm, s in specs.items():
    ls = lstar[nm]
    sa_s, sa_a = s.sat(ls), s.sat(s.lam)
    F1, F2, F3 = SQRT_G2_HI, 1.0/sa_s, sa_s/sa_a
    print(f"{nm:17s} {s.margin_oneliner(s.lam):13.4f} {F1:12.4f} {F2:12.4f} "
          f"{F3:9.4f} {F1*F2*F3:9.4f}")

print("\n(6c) lambda required to reach a given margin^(1) target (fixed shape)")
targets = [2.0, 1.75, 1.5, 1.4, 1.3446, 1.25, 1.2, 1.15, 1.1481, 1.1]
print(f"{'case':17s} " + " ".join(f"{('m=%.4g' % t):>10}" for t in targets))
for nm, s in specs.items():
    out = []
    for t in targets:
        if t <= SQRT_G2_HI:
            out.append(f"{'IMPOSS':>10}")
            continue
        lo_, hi_ = 1.0 + 1e-9, 1e14
        if s.margin_oneliner(hi_) > t:
            out.append(f"{'>1e14':>10}"); continue
        for _ in range(300):
            mid = np.sqrt(lo_ * hi_)
            if s.margin_oneliner(mid) > t:
                lo_ = mid
            else:
                hi_ = mid
            if hi_/lo_ - 1 < 1e-13:
                break
        out.append(f"{np.sqrt(lo_*hi_):10.4g}")
    print(f"{nm:17s} " + " ".join(out))
print(f"    (lambda_* for these shapes is 4.82-5.81, so every target below ~1.34")
print(f"     demands a lambda far beyond the Semenok stability threshold.)")

print("\nDONE (step 6).")


print("\n" + "=" * 78)
print("STEP 7 -- (a) how far can m^(1)(lambda_*) be pushed? random shape search")
print("=" * 78)
print("  The 6 real spectra + Einstein all landed at 1.3446-1.3465, but the 90/10")
print("  bimodal gave 3.18 -- so the near-constancy is NOT universal and the benign")
print("  sampling family was hiding it.  Search properly.")
rng2 = np.random.default_rng(7771)
best = (1e9, None); worst = (0.0, None)
vals = []
for k in range(500):
    npk = rng2.integers(1, 4)
    dens = np.zeros_like(wgrid)
    desc = []
    for _ in range(npk):
        c = rng2.uniform(2, 395)
        sg = np.exp(rng2.uniform(np.log(0.8), np.log(80.0)))
        am = np.exp(rng2.uniform(np.log(0.02), np.log(1.0)))
        dens += am * np.exp(-0.5 * ((wgrid - c) / sg) ** 2)
        desc.append(f"({c:.0f},{sg:.1f},{am:.2f})")
    ss = SynthSpec(wgrid, dens + 1e-300)
    ls, _, _ = ss.lambda_star(n_u=250, refine=2)
    sa = ss.sat(ls)
    if not np.isfinite(sa):
        continue
    m = SQRT_G2_HI / sa
    vals.append(m)
    if m < best[0]:
        best = (m, (ss.r, ls, sa, " ".join(desc)))
    if m > worst[0]:
        worst = (m, (ss.r, ls, sa, " ".join(desc)))
vals = np.array(vals)
print(f"  n = {len(vals)} random shapes (1-3 gaussian peaks, widths 0.8-80 meV,")
print(f"                                amplitudes spanning 50x)")
print(f"  m^(1)(lambda_*):  min = {vals.min():.5f}   median = {np.median(vals):.5f}   "
      f"max = {vals.max():.5f}")
print(f"  Einstein closed-form value                       = {m_E:.5f}")
print(f"  fraction of random shapes BELOW the Einstein value: "
      f"{np.mean(vals < m_E - 1e-9):.4f}  ({int(np.sum(vals < m_E - 1e-9))}/{len(vals)})")
print(f"  argmin shape (r, lambda_*, s^(1), peaks): {best[1]}")
print(f"  argmax shape (r, lambda_*, s^(1), peaks): {worst[1]}")
print("  => CONJECTURE (numerical-evidence only): Einstein MINIMISES m^(1)(lambda_*[P],P);")
print("     the quantity is unbounded above, so 1.345 is NOT a universal constant.")

print("\n" + "=" * 78)
print("STEP 7b -- is lambda<w^2> invariant under the anharmonic renormalisation?")
print("=" * 78)
print("  lambda<w^2> = 2 int w a2F(w) dw = eta/M is the ONLY spectral input to the")
print("  one-liner ceiling.  The repo claims softening cannot move it (Theorem 2).")
print("  Direct test on the harmonic/anharmonic DFT pairs:")
print(f"{'pair':28s} {'lam<w^2>_har':>14} {'lam<w^2>_anh':>14} {'rel diff':>10} "
      f"{'Tc_1L_har':>10} {'Tc_1L_anh':>10}")
for a, b, lbl in [("H3S-200GPa-har", "H3S-200GPa-anh", "H3S 200GPa har->anh"),
                  ("D3S-200GPa-har", "D3S-200GPa-anh", "D3S 200GPa har->anh")]:
    sa_, sb_ = specs[a], specs[b]
    va, vb = sa_.lam * sa_.w2, sb_.lam * sb_.w2
    ta = C_INF_HI * np.sqrt(va) / KB_MEV_PER_K
    tb = C_INF_HI * np.sqrt(vb) / KB_MEV_PER_K
    print(f"{lbl:28s} {va:14.4f} {vb:14.4f} {abs(va-vb)/va:10.4%} {ta:10.3f} {tb:10.3f}")
print("  (for reference, lambda itself moves 2.293->1.721 = -25%, and w2av 129->149 = +15%)")
print("  direct first moment 2*int w a2F dw, recomputed both quadratures:")
for nm, s in specs.items():
    m1_t = 2.0 * np.trapezoid(s.a2f * s.om, s.om)
    m1_s = 2.0 * simpson_nonuniform(s.a2f * s.om, s.om)
    print(f"    {nm:17s} 2*int w a2F = {m1_t:12.5f} (trap)  {m1_s:12.5f} (simp)  "
          f"lam*<w^2> = {s.lam*s.w2:12.5f}")

print("\nDONE (step 7).")
