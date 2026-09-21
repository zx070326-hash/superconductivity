"""Explore f(lambda, mu*) in strong coupling limit to extract mu*-dependent corrections.

Key question: Does lim_{lambda->inf} Tc/sqrt(lambda<om^2>) depend on mu*?
If not (expected), extract the 1/lambda correction: f ~ 1 - a(mu*)/lambda^p

Strategy:
1. Distinguish physical mu~ = mu*/lambda -> 0 from N_c = floor(r*varpi) -> 0 artifact
2. Use larger r to separate these effects
3. Extract asymptotic form via Richardson extrapolation or fitting

numerical-evidence (float64)."""
import os, math, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from concavity_exact_resolvent import build

CINF = 0.1827262477


def table(N, VP):
    out = []
    for vp in VP:
        n, u, A = build(N, vp * vp)[:3]
        mu, PSI = eigh(A)
        out.append((n, u, mu, PSI))
    return out


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
    lo, hi = 0.15 / math.sqrt(lam), 0.75 * math.exp(1.0 / lam) + 8.0 / math.sqrt(lam)
    VP = np.geomspace(lo, hi, M)
    recs = table(N, VP)
    for i, vp in enumerate(VP):
        if lam * Q_at(recs[i], vp, mus, lam, r) - 1.0 >= 0:
            if i == 0:
                return None
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
        return None, None, None
    tc = 1.0 / (2 * math.pi * vp)
    Nc = int(math.floor(r * vp))
    return tc, tc / (CINF * math.sqrt(lam)), Nc


if __name__ == "__main__":
    print("=" * 80)
    print("PART 1: Strong coupling behavior for various r (separate physical vs artifact)")
    print("=" * 80)

    LAMS_STRONG = [10.0, 20.0, 30.0, 50.0, 100.0]
    R_VALUES = [10.0, 20.0, 30.0, 50.0]
    MUS = [0.0, 0.10, 0.13]

    print("\nFor each (lambda, r): f(mu*) and N_c")
    print(f"{'lam':>6} {'r':>4} | {'f(0)':>8} {'Nc(0)':>6} | " +
          " | ".join(f"{'f('+str(m)+')':>8} {'Nc':>4} {'def%':>6}" for m in MUS[1:]))

    for lam in LAMS_STRONG:
        for r in R_VALUES:
            tc0, f0, nc0 = fill(lam, 0.0, r, N=400)
            results = [(tc0, f0, nc0)]
            for mus in MUS[1:]:
                tc, f, nc = fill(lam, mus, r, N=400)
                results.append((tc, f, nc))

            line = f"{lam:>6.1f} {r:>4.0f} | {f0:>8.5f} {nc0:>6d} |"
            for i, mus in enumerate(MUS[1:], 1):
                tc, f, nc = results[i]
                deficit = (1 - tc / tc0) * 100
                line += f" {f:>8.5f} {nc:>4d} {deficit:>6.1f}% |"
            print(line)

    print("\n" + "=" * 80)
    print("PART 2: Extract 1/lambda correction via fitting")
    print("=" * 80)

    # For fixed mu*, fit f(lam) ~ 1 - a/lam^p in strong coupling region
    LAMS_FIT = [10.0, 15.0, 20.0, 30.0, 50.0, 100.0]

    for mus in [0.10, 0.13]:
        print(f"\nmu* = {mus}, r = 10")
        print(f"{'lam':>6} {'f':>9} {'1-f':>9} {'Nc':>4} | {'lam*(1-f)':>11} {'lam^2*(1-f)':>13}")

        data = []
        for lam in LAMS_FIT:
            tc, f, nc = fill(lam, mus, 10.0, N=400)
            data.append((lam, f, nc))
            gap = 1 - f
            print(f"{lam:>6.1f} {f:>9.6f} {gap:>9.6f} {nc:>4d} | "
                  f"{lam * gap:>11.6f} {lam * lam * gap:>13.6f}")

        # Simple power law fit: 1 - f ~ a / lam^p
        # Try p=1 and p=2
        from scipy.optimize import curve_fit

        lams = np.array([d[0] for d in data])
        gaps = np.array([1 - d[1] for d in data])

        # p=1: 1-f = a/lam
        def model_p1(lam, a):
            return a / lam

        # p=2: 1-f = a/lam^2
        def model_p2(lam, a):
            return a / lam**2

        # Mixed: 1-f = a1/lam + a2/lam^2
        def model_mixed(lam, a1, a2):
            return a1 / lam + a2 / lam**2

        try:
            popt_p1, _ = curve_fit(model_p1, lams, gaps)
            popt_p2, _ = curve_fit(model_p2, lams, gaps)
            popt_mixed, _ = curve_fit(model_mixed, lams, gaps)

            print(f"\n  Fit 1-f = a/λ:           a = {popt_p1[0]:.6f}")
            print(f"  Fit 1-f = a/λ²:          a = {popt_p2[0]:.6f}")
            print(f"  Fit 1-f = a₁/λ + a₂/λ²:  a₁ = {popt_mixed[0]:.6f}, a₂ = {popt_mixed[1]:.6f}")

            # Residuals
            res_p1 = np.abs(gaps - model_p1(lams, *popt_p1))
            res_p2 = np.abs(gaps - model_p2(lams, *popt_p2))
            res_mixed = np.abs(gaps - model_mixed(lams, *popt_mixed))

            print(f"  Max residual p=1: {np.max(res_p1):.2e}")
            print(f"  Max residual p=2: {np.max(res_p2):.2e}")
            print(f"  Max residual mixed: {np.max(res_mixed):.2e}")
        except Exception as e:
            print(f"  Fit failed: {e}")

    print("\n" + "=" * 80)
    print("PART 3: Compare mu*=0.10 and mu*=0.13 directly")
    print("=" * 80)

    print(f"\n{'lam':>6} | {'f(0.10)':>9} {'f(0.13)':>9} {'ratio':>9} | "
          f"{'def(0.10)%':>11} {'def(0.13)%':>11}")

    for lam in LAMS_FIT:
        tc0, f0, _ = fill(lam, 0.0, 10.0, N=400)
        tc1, f1, _ = fill(lam, 0.10, 10.0, N=400)
        tc2, f2, _ = fill(lam, 0.13, 10.0, N=400)

        ratio = f1 / f2
        def1 = (1 - tc1 / tc0) * 100
        def2 = (1 - tc2 / tc0) * 100

        print(f"{lam:>6.1f} | {f1:>9.6f} {f2:>9.6f} {ratio:>9.6f} | "
              f"{def1:>11.2f} {def2:>11.2f}")
