r"""
Numerical probe: non-Lorentz kernels in the generic operator O[f].
NUMERICAL-EVIDENCE ONLY (float64 eigvalsh, finite N=200, finite w-grid).
No certificate claims here -- this is exploratory, per team-lead request
(explore-non-lorentz support task).

Framework (P4_rigorous_tc_bound.md Sec.0, verbatim):
  O[f] := -L1[f] + L2[f] + L3[f]  on l^2(N_0), u_n := (2n+1)^{-1/2}
    L1_{nn} = (2n+1)^{-1} sum_{k=1}^n 2 f(k)              (diagonal)
    L2_{nm} = (1-delta_nm) f(|n-m|) u_n u_m
    L3_{nm} = f(n+m+1) u_n u_m
  Standard (KAY) Lorentz kernel: [[k]](w) = w^2/(k^2+w^2).
  H(w) := O[[[.]](w)],  h(w) := lambda_1(H(w)).
  H(x) := h(sqrt(x)) is proved concave on (0,infty) for the Lorentz kernel
  (Lemma C / "the concavity theorem", certificate-level in the repo).

This script asks: for OTHER single-argument kernel families f(k;w,params)
(still plugged into the SAME O[f] machinery), does h(w) keep the same
qualitative shape (log-type large-w growth, quadratic small-w growth,
H(x)=h(sqrt(x)) concave)?  All kernels reduce to the standard Lorentzian
at their "trivial" parameter value, which is used as an in-script sanity
check against the KNOWN closed-form asymptotics quoted in the theory notes:
  r(0+) = h(w)/w^2 -> g(2) = 1.3181402... as w->0   (recorded constant)
  h(w) ~ ln(w) + gamma_E + ln4 - 3/2 - 1/(12w) - 3/(16w^2)   as w->infty
  (P4_rigorous_tc_bound.md Sec.2 recap (iv)/(iii))
"""
import numpy as np

np.seterr(all="ignore")

EULER_GAMMA = 0.5772156649015329


# ============================================================================
# Generic O[f] builder (structure copied from p4_interchange_uniformity.H_mat,
# generalized: the kernel is passed in as a plain array function of the
# integer argument k, already bound to whatever parameters (w, alpha, ...)
# the caller wants).
# ============================================================================

def build_O(N, kfun):
    """O[f] as an N x N float64 matrix, f given by kfun(k_array) -> values."""
    n = np.arange(N, dtype=np.float64)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    nrm = np.outer(u, u)
    Dm = np.abs(n[:, None] - n[None, :])
    Sm = n[:, None] + n[None, :] + 1.0
    off = kfun(Dm)
    np.fill_diagonal(off, 0.0)
    M = (off + kfun(Sm)) * nrm
    ks = np.arange(1, N, dtype=np.float64)
    cum = np.concatenate(([0.0], np.cumsum(2.0 * kfun(ks))))
    M[np.diag_indices(N)] -= cum / (2.0 * n + 1.0)
    return M


def top_eig(N, kfun):
    return float(np.linalg.eigvalsh(build_O(N, kfun))[-1])


# ============================================================================
# Kernel families.  All reduce to the standard Lorentzian [[k]](w)=w^2/(k^2+w^2)
# at the "trivial" parameter (alpha=0, Ncut=infty, gamma=2).
# ============================================================================

def k_lorentz(k, w):
    return w * w / (k * k + w * w)


def k_yukawa(k, w, alpha):
    """Shifted-Lorentzian ("Yukawa-like" detuning/screening of the mode index):
    kk = k+alpha in the same w^2/(kk^2+w^2) form.  alpha=0 -> standard KAY."""
    kk = k + alpha
    return w * w / (kk * kk + w * w)


def k_expcut(k, w, Ncut):
    """Standard Lorentzian x exp(-k/Ncut) (finite-bandwidth cutoff).
    Ncut->infty -> standard KAY."""
    return (w * w / (k * k + w * w)) * np.exp(-k / Ncut)


def k_pow(k, w, gamma):
    """Power-law generalization of the Lorentzian: w^gamma/(k^gamma+w^gamma).
    gamma=2 -> standard KAY (Lorentzian); gamma=1 Ohmic-like, gamma=3 super-Ohmic-like."""
    return w ** gamma / (k ** gamma + w ** gamma)


KERNELS = {}
for a in (0.0, 0.5, 1.0, 2.0):
    KERNELS[f"yukawa(alpha={a})"] = (lambda w, a=a: (lambda k: k_yukawa(k, w, a)))
for Nc in (50.0, 100.0, 500.0):
    KERNELS[f"expcut(Ncut={Nc:.0f})"] = (lambda w, Nc=Nc: (lambda k: k_expcut(k, w, Nc)))
for g in (1.0, 1.5, 2.0, 3.0):
    KERNELS[f"pow(gamma={g})"] = (lambda w, g=g: (lambda k: k_pow(k, w, g)))
KERNELS["lorentz(baseline)"] = (lambda w: (lambda k: k_lorentz(k, w)))


# ============================================================================
# h(w) sampling + asymptotics + H(x)=h(sqrt(x)) concavity test
# ============================================================================

N = 200


def h_curve(kernel_maker, wgrid):
    return np.array([top_eig(N, kernel_maker(w)) for w in wgrid])


def fit_loglog(w, h, lo, hi):
    """Least-squares log(h) ~ alpha*log(w)+c on w in [lo,hi]; returns (alpha, resid)."""
    m = (w >= lo) & (w <= hi) & (h > 0)
    if m.sum() < 4:
        return None, None
    x = np.log(w[m]); y = np.log(h[m])
    A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    resid = float(np.sqrt(np.mean((pred - y) ** 2)))
    return float(coef[0]), resid


def fit_log_of_w(w, h, lo, hi):
    """Least-squares h ~ A*ln(w)+B on w in [lo,hi]; returns (A,B,resid)."""
    m = (w >= lo) & (w <= hi)
    x = np.log(w[m]); y = h[m]
    A = np.vstack([x, np.ones_like(x)]).T
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    resid = float(np.sqrt(np.mean((pred - y) ** 2)))
    return float(coef[0]), float(coef[1]), resid


def concavity_report(kernel_maker, xgrid):
    """H(x)=h(sqrt(x)); 2nd difference via the standard 3-point formula for a
    NON-uniform grid (xgrid is expected log-spaced so it can resolve x->0):
        H''(x1) ~= 2*[ (h2-h1)/(x2-x1) - (h1-h0)/(x1-x0) ] / (x2-x0)
    which reduces to the usual central difference when the spacing is uniform."""
    w = np.sqrt(xgrid)
    h = h_curve(kernel_maker, w)
    x0, x1, x2 = xgrid[:-2], xgrid[1:-1], xgrid[2:]
    h0, h1, h2 = h[:-2], h[1:-1], h[2:]
    Hpp = 2.0 * ((h2 - h1) / (x2 - x1) - (h1 - h0) / (x1 - x0)) / (x2 - x0)
    return h, Hpp


def run_family(name, kernel_maker, wgrid_log, xgrid_lin):
    print(f"\n=== {name} ===")
    h = h_curve(kernel_maker, wgrid_log)
    # asymptotic characterization: try log-fit vs power-fit on the two largest
    # decades of the grid
    lo, hi = wgrid_log[-1] / 100.0, wgrid_log[-1]
    alpha_pow, resid_pow = fit_loglog(wgrid_log, h, lo, hi)
    A_log, B_log, resid_log = fit_log_of_w(wgrid_log, h, lo, hi)
    # small-w check: h/w^2 -> const?
    small = wgrid_log <= 0.3
    r_small = h[small] / wgrid_log[small] ** 2 if small.any() else np.array([])
    idxs = [0, len(wgrid_log) // 4, len(wgrid_log) // 2, -1]
    print(f"  h(w) at w={wgrid_log[idxs]} = {h[idxs]}")
    print(f"  power-law fit on w in [{lo:.3g},{hi:.3g}]: h~A w^alpha, alpha={alpha_pow:.4f}"
          f" (log-log rmse={resid_pow:.2e})")
    print(f"  log fit         on w in [{lo:.3g},{hi:.3g}]: h~A ln w + B, A={A_log:.4f} B={B_log:.4f}"
          f" (rmse={resid_log:.2e})")
    if r_small.size:
        # report at the SMALLEST w in the window (closest to the w->0 limit);
        # r(w) is monotonically decreasing so this is the tightest estimate.
        print(f"  small-w r(w)=h/w^2 -> {r_small[0]:.6f} (w={wgrid_log[small][0]:.3f}),"
              f" spread over small-w window = {r_small.max()-r_small.min():.2e}")
    # concavity of H(x)=h(sqrt(x))
    hH, Hpp = concavity_report(kernel_maker, xgrid_lin)
    worst = float(Hpp.max())
    worst_x = float(xgrid_lin[1:-1][np.argmax(Hpp)])
    frac_pos = float(np.mean(Hpp > 1e-9))
    print(f"  H''(x) over x in [{xgrid_lin[0]:.3g},{xgrid_lin[-1]:.3g}]:"
          f" max={worst:.4e} at x={worst_x:.4g}, min={Hpp.min():.4e},"
          f" frac(H''>1e-9)={frac_pos:.3f}")
    verdict = "CONCAVE (numerical evidence)" if worst <= 1e-6 else "NOT CONCAVE (violation found)"
    print(f"  VERDICT: {verdict}")
    return dict(name=name, alpha_pow=alpha_pow, resid_pow=resid_pow, A_log=A_log,
                B_log=B_log, resid_log=resid_log, r_small=(r_small[0] if r_small.size else None),
                Hpp_max=worst, Hpp_max_x=worst_x, Hpp_min=float(Hpp.min()),
                frac_pos=frac_pos, verdict=verdict)


if __name__ == "__main__":
    wgrid_log = np.geomspace(0.001, 100.0, 80)
    # log-spaced x-grid: must resolve x->0 (small-w power-law regime, where the
    # gamma>2 concavity violation lives) AND cover large x (log-growth regime).
    xgrid_lin = np.geomspace(1e-6, 20000.0, 1200)  # x=w^2, w in [1e-3, ~141]

    print("N =", N, " w-grid:", wgrid_log[0], "...", wgrid_log[-1], f"({len(wgrid_log)} pts)")
    print("x-grid (for H'', log-spaced):", xgrid_lin[0], "...", xgrid_lin[-1], f"({len(xgrid_lin)} pts)")

    rows = []
    # baseline sanity check FIRST
    rows.append(run_family("lorentz(baseline)", KERNELS["lorentz(baseline)"], wgrid_log, xgrid_lin))
    print("\n  [baseline sanity] known closed forms (P4_rigorous_tc_bound.md):")
    print(f"    r(0+) should -> g(2)=1.3181402...; got {rows[0]['r_small']:.6f}")
    print(f"    large-w log-fit A should -> 1.0 (h ~ ln w + gamma_E + ln4 - 3/2 + ...); got A={rows[0]['A_log']:.4f}")
    print(f"    (gamma_E + ln4 - 3/2 = {EULER_GAMMA + np.log(4) - 1.5:.4f}; fitted B={rows[0]['B_log']:.4f})")

    for name in KERNELS:
        if name == "lorentz(baseline)":
            continue
        rows.append(run_family(name, KERNELS[name], wgrid_log, xgrid_lin))

    print("\n\n=== SUMMARY TABLE ===")
    hdr = "%22s %10s %8s %8s %9s %11s %40s" % ("kernel", "pow-alpha", "log-A", "log-B", "r(0+)", "H'' max", "verdict")
    print(hdr)
    for r in rows:
        print("%22s %10.4f %8.4f %8.4f %9.5f %11.4e %40s" % (
            r["name"], r["alpha_pow"], r["A_log"], r["B_log"],
            (r["r_small"] if r["r_small"] is not None else float("nan")),
            r["Hpp_max"], r["verdict"]))
