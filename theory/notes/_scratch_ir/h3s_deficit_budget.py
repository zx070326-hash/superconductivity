"""THE deficit budget: solve the KAY Tc criterion for the ACTUAL H3S/D3S spectral
measure P (from first-hand EPW a2F) and compare to the Einstein spectrum with the
same (lambda, <om^2>).  Zero free parameters except mu*.

F_P(k) = int [[k]](om/2piT) P(dom) = sum_i w_i x_i/(k^2+x_i),  x_i=(om_i/2piT)^2,
P(dom) = (2/lambda) a2F/om dom.  Tc: lambda*lam_1(O[F_P] - 2(mu*/lambda)u_c u_c^T) = 1,
Nc = floor(om_c/2piT) = floor(r*varpi_rms).

Theorem 2 says Q(P) <= Q(Einstein) => Tc(P) <= Tc(E); the GAP is the broadening
deficit, here computed exactly instead of bounded.
numerical-evidence (float64; N and quadrature convergence printed)."""
import os, sys, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from h3s_fill_factor import moments, MEV_K, ROOT
from fill_factor_probe import CINF

CINF_ = CINF


def O_of(fker, N):
    """KAY kernel-generic operator O[f]; fker is an array of length >= 2N+1."""
    n = np.arange(N); u = 1.0 / np.sqrt(2.0 * n + 1.0)
    absd = np.abs(n[:, None] - n[None, :]); sumd = n[:, None] + n[None, :] + 1
    uu = u[:, None] * u[None, :]
    Am = (fker[absd] + fker[sumd]) * uu
    L1 = np.zeros(N); L1[1:] = np.cumsum(2.0 * fker[1:N])[:N - 1] / (2 * n[1:] + 1)
    np.fill_diagonal(Am, -L1 + fker[2 * n + 1] * u ** 2)
    return n, u, Am


def FP(ks2, xs, ws):
    """F_P(k) for k^2 = ks2 (vector); xs = (om_i/2piT)^2, ws = P weights."""
    out = np.einsum('i,ki->k', ws, xs[None, :] / (ks2[:, None] + xs[None, :]))
    return out


def Q_P(T_mev, w, a2f, lam, mus, r, N=400, sqom2=None):
    """lam_1 of O[F_P] - 2 mu~ u_c u_c^T at temperature T (meV)."""
    ws = (2.0 / lam) * (a2f / w)                       # unnormalised P density
    ws = ws / np.trapezoid(ws, w) if False else ws     # keep trapezoid weights below
    # trapezoid weights for int . dom
    if len(w) == 1:
        p = np.array([1.0])                            # single atom = Einstein
    else:
        tw = np.zeros_like(w); tw[1:-1] = 0.5 * (w[2:] - w[:-2])
        tw[0] = 0.5 * (w[1] - w[0]); tw[-1] = 0.5 * (w[-1] - w[-2])
        p = ws * tw
        p = p / p.sum()                                # exact probability weights
    xs = (w / (2 * math.pi * T_mev)) ** 2
    ks = np.arange(2 * N + 1, dtype=float)
    fker = FP(ks ** 2, xs, p); fker[0] = 1.0
    n, u, A = O_of(fker, N)
    if mus == 0.0:
        return float(eigh(A, eigvals_only=True)[-1])
    vp = (sqom2 if sqom2 else math.sqrt(np.sum(p * w ** 2))) / (2 * math.pi * T_mev)
    Nc = int(math.floor(r * vp))
    if Nc <= 0:
        return float(eigh(A, eigvals_only=True)[-1])
    mu, PSI = eigh(A)
    c = PSI.T @ np.where(n < Nc, u, 0.0)
    m1, m2 = mu[-1], mu[-2]
    wdg = (m1 - m2) * 1e-13 + 1e-300
    mut = mus / lam
    phi = lambda l: 2 * mut * float(np.sum(c ** 2 / (mu - l))) - 1.0
    return brentq(phi, m2 + wdg, m1 - wdg, xtol=1e-16, rtol=1e-15)


def Tc_P(w, a2f, lam, mus, r, N=400, sqom2=None):
    """largest T (meV) with lambda*Q >= 1; Q decreasing in T."""
    g = np.geomspace(0.05, 2.0 * math.sqrt(np.sum(w ** 2) / len(w)), 160)
    prev = None
    for T in g:
        v = lam * Q_P(T, w, a2f, lam, mus, r, N, sqom2) - 1.0
        if v < 0:
            if prev is None:
                return None
            a, b = prev, T
            for _ in range(40):
                m = 0.5 * (a + b)
                if lam * Q_P(m, w, a2f, lam, mus, r, N, sqom2) - 1.0 >= 0:
                    a = m
                else:
                    b = m
            return a
        prev = T
    return None


if __name__ == "__main__":
    # ---- validation: general-P code on a single atom must reproduce fill_factor_probe
    from fill_factor_probe import fill as _fill
    print("validation: single-atom P via general-P code  vs  fill_factor_probe (Einstein)")
    for L in (1.0, 1.721, 2.293):
        tA = Tc_P(np.array([1.0]), np.array([L / 2.0]), L, 0.0, 10.0, 400, 1.0)
        tB, _ = _fill(L, 0.0)
        print(f"  lam={L:>6.3f}  Tc/omE: general-P={tA:.8f}  probe={tB:.8f}  "
              f"rel={abs(tA-tB)/tB:.2e}")
    print()
    print("Tc from the ACTUAL spectral measure vs the Einstein spectrum, same (lam,<om2>)")
    print(f"{'dataset':>26} {'lam':>6} {'sq<om2>':>8} | {'ceiling':>8} "
          f"{'Tc_E(0)':>8} {'Tc_P(0)':>8} {'broad':>7} | {'Tc_P(.13)':>9} {'Tc_P(.20)':>9} "
          f"| {'Tc_exp':>7} {'f_exp':>6}")
    EXP = {"H3S-200GPa": 203.0, "D3S-200GPa": 150.0, "H3S-160GPa": 203.0}
    for tag in ("H3S-200GPa", "D3S-200GPa", "H3S-160GPa"):
        for sub in ("EPW_anharmonic", "EPW_harmonic"):
            import glob
            g = glob.glob(f"{ROOT}/{tag}/{sub}/*.a2f")
            g = [x for x in g if not x.endswith(".a2fv")]
            if not g:
                continue
            lam, om2, olog, w, a2f = moments(g[0])
            sq = math.sqrt(om2)
            ceil_K = CINF_ * math.sqrt(lam * om2) * MEV_K
            # Einstein at the same (lam, <om2>): delta at sqrt<om2>
            tE = Tc_P(np.array([sq]), np.array([lam * sq / 2.0]), lam, 0.0, 10.0, 400, sq)
            tP0 = Tc_P(w, a2f, lam, 0.0, 10.0, 400, sq)
            tP1 = Tc_P(w, a2f, lam, 0.13, 10.0, 400, sq)
            tP2 = Tc_P(w, a2f, lam, 0.20, 10.0, 400, sq)
            br = (1 - tP0 / tE) * 100 if (tE and tP0) else float("nan")
            te = EXP.get(tag, float("nan"))
            print(f"{tag+'/'+sub:>26} {lam:>6.3f} {sq:>8.2f} | {ceil_K:>8.1f} "
                  f"{tE*MEV_K:>8.1f} {tP0*MEV_K:>8.1f} {br:>6.1f}% | "
                  f"{tP1*MEV_K:>9.1f} {tP2*MEV_K:>9.1f} | {te:>7.1f} {te/ceil_K:>6.3f}",
                  flush=True)

    print("\nN convergence (H3S-200 anharmonic, mu*=0.13)")
    import glob
    lam, om2, olog, w, a2f = moments(glob.glob(f"{ROOT}/H3S-200GPa/EPW_anharmonic/*.a2f")[0])
    sq = math.sqrt(om2)
    for N in (200, 400, 800):
        t = Tc_P(w, a2f, lam, 0.13, 10.0, N, sq)
        print(f"  N={N}: Tc = {t*MEV_K:.3f} K")
    print("quadrature: drop every 2nd a2F point (grid halving)")
    t = Tc_P(w[::2], a2f[::2], lam, 0.13, 10.0, 400, sq)
    print(f"  half grid: Tc = {t*MEV_K:.3f} K")
