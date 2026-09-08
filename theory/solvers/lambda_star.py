r"""
theory/solvers/lambda_star.py

Productized implementation of the Semenok-Altshuler-Yuzbashyan stability
upper bound functional lambda_*[alpha2F] (Eq. S12), plus the P3 Sec.4
sigma/skewness expansion around the Einstein point.

    lambda_*[P] = 1 / max_s int_0^inf g(s*w) P(w) dw
    g(x) = 6x^2 + 12x^3 Im psi'(ix) + 6x^4 Re psi''(ix)

P(w) is the normalized Eliashberg shape weight P(w) = (1/lambda) * 2*alpha2F(w)/w
(integrates to 1 over w>0, i.e. P ~ w for a Debye alpha2F~w^2 spectrum);
lambda_* depends only on the *shape* of P, not its overall scale.

Provenance: derived in theory/derivations/p3_lambda_star_kappa.py, independently
re-derived/cross-checked in theory/derivations/_xcheck_p3_kappa_A.py / _B.py
(same closed forms, same x_max/kappa numbers). This module productizes that
work behind a clean API with unit tests (theory/solvers/tests/test_lambda_star.py).

Hard anchors (P3_tc_max_pressure.md Sec.4, Semenok 2407.12922 Eq.9):
    Einstein delta-function spectrum         -> lambda_* = 3.6915 (= 1/g_max, lower bound)
    Debye spectrum (alpha2F ~ w^2, P(w) ~ w) -> lambda_* = 4.7186
"""
import mpmath as mp

mp.mp.dps = 25


# ---------------------------------------------------------------------------
# g(x) and its derivatives (closed form, P3_tc_max_pressure.md Sec.4)
# ---------------------------------------------------------------------------

def _Ipsi(n, x):
    """Im psi^(n)(i x)."""
    return mp.im(mp.polygamma(n, 1j * mp.mpf(x)))


def _Rpsi(n, x):
    """Re psi^(n)(i x)."""
    return mp.re(mp.polygamma(n, 1j * mp.mpf(x)))


def g(x):
    x = mp.mpf(x)
    return 6*x**2 + 12*x**3*_Ipsi(1, x) + 6*x**4*_Rpsi(2, x)


def g1(x):
    """g'(x)."""
    x = mp.mpf(x)
    return 12*x + 36*x**2*_Ipsi(1, x) + 36*x**3*_Rpsi(2, x) - 6*x**4*_Ipsi(3, x)


def g2(x):
    """g''(x)."""
    x = mp.mpf(x)
    return 12 + 72*x*_Ipsi(1, x) + 144*x**2*_Rpsi(2, x) - 60*x**3*_Ipsi(3, x) - 6*x**4*_Rpsi(4, x)


def g3(x):
    """g'''(x)."""
    x = mp.mpf(x)
    return (72*_Ipsi(1, x) + 360*x*_Rpsi(2, x) - 324*x**2*_Ipsi(3, x)
            - 84*x**3*_Rpsi(4, x) + 6*x**4*_Ipsi(5, x))


# ---------------------------------------------------------------------------
# x_max, g_max, lambda0, kappa, c3 -- computed once at import time
# ---------------------------------------------------------------------------

X_MAX = mp.findroot(g1, mp.mpf('0.33'))
G_MAX = g(X_MAX)
LAMBDA0 = 1 / G_MAX                       # = lambda_*(Einstein) = 3.6915 (lower bound)
KAPPA = -X_MAX**2 * g2(X_MAX) / G_MAX     # sigma^2 coefficient of lambda_*/lambda0 is KAPPA/2
C3 = -X_MAX**3 * g3(X_MAX) / (6 * G_MAX)  # sigma^3 (skewness) coefficient


# ---------------------------------------------------------------------------
# Exact lambda_*[P] for a discrete or continuous spectrum shape
# ---------------------------------------------------------------------------

def _maximize_G(Gf, s_guess):
    """max_s Gf(s), found via mp.findroot on the stationarity condition
    d(Gf)/ds = 0, starting near s_guess (typically X_MAX/mean(w))."""
    s0 = mp.findroot(lambda s: mp.diff(Gf, s), mp.mpf(s_guess))
    return s0, Gf(s0)


def lambda_star_discrete(ws, ps):
    """Exact lambda_*[P] for a discrete spectrum given as parallel lists/arrays
    of frequencies ws (>0) and (possibly unnormalized) weights ps (>=0), i.e.
    P(w) = sum_i p_i delta(w-w_i) after normalizing sum(p_i)=1.
    Returns lambda_* as an mpmath mpf."""
    pts = [(mp.mpf(w), mp.mpf(p)) for w, p in zip(ws, ps) if p > 0 and w > 0]
    if not pts:
        raise ValueError("lambda_star_discrete: no positive-weight points")
    tot = sum(p for _, p in pts)
    pts = [(w, p / tot) for w, p in pts]
    mean = sum(w * p for w, p in pts)
    Gf = lambda s: sum(p * g(mp.mpf(s) * w) for w, p in pts)
    _, Gval = _maximize_G(Gf, X_MAX / mean)
    return 1 / Gval


def lambda_star_continuous(P, domain):
    """Exact lambda_*[P] for a continuous shape function P(w) (need not be
    pre-normalized) on the union of intervals `domain` (a list of breakpoints
    passed straight to mp.quad, e.g. [0, w1, w2, ...] for a piecewise-defined
    integrand -- see mpmath.quad docs). Returns lambda_* as an mpmath mpf."""
    Zn = mp.quad(P, domain)
    Pn = lambda w: P(w) / Zn
    mean = mp.quad(lambda w: w * Pn(w), domain)
    Gf = lambda s: mp.quad(lambda w: Pn(w) * g(mp.mpf(s) * w), domain)
    _, Gval = _maximize_G(Gf, X_MAX / mean)
    return 1 / Gval


# ---------------------------------------------------------------------------
# Moments: sigma (relative width), skewness, r = omega_2/omega_log
# ---------------------------------------------------------------------------

def moments_discrete(ws, ps):
    """Return (mean, sigma, skewness, r) for a discrete P(w) weight.
    sigma = sqrt(<w^2>-<w>^2)/<w>        (relative std dev)
    skewness = third standardized moment (gamma_3)
    r = omega_2/omega_log = sqrt(<w^2>)/exp(<ln w>)   (both moments computed
        with the SAME P(w) weight -- P3 Sec.4/Sec.7 sigma<->r convention,
        r-1 = sigma^2 - (1/3)*skew*sigma^3 + O(sigma^4))."""
    pts = [(mp.mpf(w), mp.mpf(p)) for w, p in zip(ws, ps) if p > 0 and w > 0]
    tot = sum(p for _, p in pts)
    pts = [(w, p / tot) for w, p in pts]
    m1 = sum(p * w for w, p in pts)
    m2 = sum(p * w**2 for w, p in pts)
    m3 = sum(p * w**3 for w, p in pts)
    lnm = sum(p * mp.log(w) for w, p in pts)
    var = m2 - m1**2
    sigma = mp.sqrt(var) / m1
    mu3 = m3 - 3*m1*m2 + 2*m1**3
    skew = mu3 / (m1*sigma)**3 if sigma > 0 else mp.mpf(0)
    r = mp.sqrt(m2) / mp.e**lnm
    return m1, sigma, skew, r


# ---------------------------------------------------------------------------
# sigma/skewness expansion (P3_tc_max_pressure.md Sec.4, v0.3.2)
# ---------------------------------------------------------------------------

def lambda_star_expansion(sigma, skew=0.0):
    """lambda_*/lambda0 = 1 + (kappa/2) sigma^2 + c3*skew*sigma^3 + O(sigma^4).

    Valid to <0.4% for symmetric spectra (skew=0) at sigma<~0.2; for spectra
    with a low-frequency tail (skew<0, e.g. Debye/acoustic), only accurate to
    ~1-3% at sigma~0.2 and degrades further for sigma>~0.3 -- use
    lambda_star_discrete/_continuous directly for such spectra
    (see P3_tc_max_pressure.md Sec.4/Sec.7 for the full caveat)."""
    sigma = mp.mpf(sigma)
    skew = mp.mpf(skew)
    return LAMBDA0 * (1 + KAPPA/2 * sigma**2 + C3 * skew * sigma**3)


if __name__ == "__main__":
    print(f"x_max   = {mp.nstr(X_MAX, 10)}")
    print(f"g_max   = {mp.nstr(G_MAX, 10)}")
    print(f"lambda0 = {mp.nstr(LAMBDA0, 10)}   (Einstein anchor, lit 3.6915)")
    print(f"kappa   = {mp.nstr(KAPPA, 10)}   (lit 5.6035)")
    print(f"c3      = {mp.nstr(C3, 10)}   (lit 0.8959)")
    debye = lambda_star_continuous(lambda w: w, [0, 1])
    print(f"Debye (P~w on (0,1]) lambda_* = {mp.nstr(debye, 8)}   (lit 4.7186)")
