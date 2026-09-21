"""Fill factor f := Tc / (C_inf sqrt(lambda <om^2>)) for the Einstein spectrum,
mu* = 0 and mu* > 0, from the exact KAY operator.

One eigh(A) per grid point in varpi (A = O[[[.]]](varpi)); lam_1(M) for M =
A - 2 mu~ u_c u_c^T then comes from the SECULAR ROOT (O(N) per mu~), same device
as Theorem CE' / hpp_fast.  Tc <-> smallest varpi with lambda*lam_1(M) >= 1
(largest T), Nc = floor(r*varpi).
numerical-evidence (float64; N and grid convergence printed)."""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concavity_exact_resolvent import build

CINF = 0.1827262477            # (1/2pi) sqrt(g2), analytic
_TAB = {}


def table(N, VP):
    """for each varpi: (n, u, mu, c-coeffs not yet -- need Nc-dependent uc)"""
    key = (N, len(VP), VP[0], VP[-1])
    if key not in _TAB:
        out = []
        for vp in VP:
            n, u, A = build(N, vp * vp)[:3]
            mu, PSI = eigh(A)
            out.append((n, u, mu, PSI))
        _TAB[key] = out
    return _TAB[key]


def Q_at(rec, vp, mus, lam, r):
    n, u, mu, PSI = rec
    if mus == 0.0:
        return mu[-1]
    Nc = int(math.floor(r * vp))
    if Nc <= 0:
        return mu[-1]
    mut = mus / lam
    c = PSI.T @ np.where(n < Nc, u, 0.0)
    m1, m2 = mu[-1], mu[-2]
    w = (m1 - m2) * 1e-13 + 1e-300
    phi = lambda l: 2 * mut * float(np.sum(c ** 2 / (mu - l))) - 1.0
    return brentq(phi, m2 + w, m1 - w, xtol=1e-16, rtol=1e-15)


def varpi_c(lam, mus, r, N=400, M=900):
    """Q is increasing in varpi; Tc <-> smallest varpi with lam*Q >= 1."""
    lo, hi = 0.15 / math.sqrt(lam), 0.75 * math.exp(1.0 / lam) + 8.0 / math.sqrt(lam)
    VP = np.geomspace(lo, hi, M)
    recs = table(N, VP)
    for i, vp in enumerate(VP):
        if lam * Q_at(recs[i], vp, mus, lam, r) - 1.0 >= 0:
            if i == 0:
                return None
            # refine within [VP[i-1], VP[i]] using the finer of the two records
            a, b = VP[i - 1], VP[i]
            for _ in range(28):
                m = 0.5 * (a + b)
                n_, u_, A_ = build(N, m * m)[:3]
                mu_, PSI_ = eigh(A_)
                q = Q_at((n_, u_, mu_, PSI_), m, mus, lam, r)
                if lam * q - 1.0 >= 0:
                    b = m
                else:
                    a = m
            return b
    return None


def fill(lam, mus=0.0, r=10.0, N=400):
    vp = varpi_c(lam, mus, r, N)
    if vp is None:
        return None, None
    tc = 1.0 / (2 * math.pi * vp)                 # Tc / omega_E
    return tc, tc / (CINF * math.sqrt(lam))


if __name__ == "__main__":
    LAMS = (0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0, 30.0, 100.0)
    print(f"C_inf = {CINF}   ceiling/omega_E = C_inf*sqrt(lambda)   N=400, r=10")
    print(f"{'lam':>6} {'varpi_c':>9} {'Tc/omE':>9} {'f(mu*=0)':>9} | "
          f"{'Tc(.13)':>9} {'f':>7} {'defic':>7} | {'Tc(.10)':>9} {'f':>7} {'defic':>7}")
    for lam in LAMS:
        t0, f0 = fill(lam, 0.0)
        t1, f1 = fill(lam, 0.13, 10.0)
        t2, f2 = fill(lam, 0.10, 10.0)
        d1 = (1 - t1 / t0) * 100
        d2 = (1 - t2 / t0) * 100
        print(f"{lam:>6.1f} {1/(2*math.pi*t0):>9.4f} {t0:>9.5f} {f0:>9.4f} | "
              f"{t1:>9.5f} {f1:>7.4f} {d1:>6.1f}% | {t2:>9.5f} {f2:>7.4f} {d2:>6.1f}%",
              flush=True)

    print("\nasymptotic check: varpi_c -> 1/sqrt(g2*lam) = 0.87088/sqrt(lam), f -> 1")
    for lam in (10.0, 100.0, 1000.0):
        vp = varpi_c(lam, 0.0, 10.0)
        print(f"  lam={lam:>7.0f}  varpi_c={vp:.6f}  pred={0.87088/math.sqrt(lam):.6f}"
              f"   f={1/(2*math.pi*vp)/(CINF*math.sqrt(lam)):.6f}")

    print("\nweak-coupling check: Tc/omE -> 0.2530*exp(-1/lam)*(1+0.3584 lam)")
    for lam in (0.3, 0.4, 0.5, 0.7):
        t0, f0 = fill(lam, 0.0)
        w = 0.2529996837 * math.exp(-1.0 / lam)
        print(f"  lam={lam:>4.2f}  Tc/omE={t0:.6f}  W(lead)={w:.6f}  "
              f"W(1+b1 lam)={w*(1+0.3584*lam):.6f}  f={f0:.4f}")

    print("\nN / grid convergence at lam=2.0")
    for N in (300, 400, 600):
        t0, f0 = fill(2.0, 0.0, 10.0, N)
        t1, f1 = fill(2.0, 0.13, 10.0, N)
        print(f"  N={N}: f(0)={f0:.6f}  f(.13)={f1:.6f}  defic={(1-t1/t0)*100:.2f}%")
