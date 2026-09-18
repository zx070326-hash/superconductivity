"""Acceptance-gate test for the expansion-based (non-triangle) third-term fix,
per team-lead's round-3 correction:

  ||rho||^2 = ||r||^2 - 2<r,Ty> + ||Ty||^2   (EXACT expansion, no triangle ineq)

For a trial-vector-only r(v) (v unknown, only sinTheta(tau,v)<=s known), replace
<z,r(v)>  (z = y or z = Ty, arbitrary, not necessarily perp tau) by its
sup/inf over the cone v = sqrt(1-s^2)*tau + s*w, ||w||=1, w perp tau:

  decompose z = a*tau + z_perp  (a = <z,tau>, z_perp perp tau)
  <z_perp, r(v)> in [ sqrt(1-s^2)<z_perp,Opsi tau> - 2 s npsi ||z_perp||,
                       sqrt(1-s^2)<z_perp,Opsi tau> + 2 s npsi ||z_perp|| ]
  <tau, r(v)>   in [ sqrt(1-s^2)s^2 qpsi(tau) - (s|2s^2-1| e_psi + sqrt(1-s^2)s^2 npsi),
                      sqrt(1-s^2)s^2 qpsi(tau) + (s|2s^2-1| e_psi + sqrt(1-s^2)s^2 npsi) ]

Both derived by hand from r(v)=Opsi v - v qpsi(v), v=sqrt(1-s^2)tau+s w, bounding
<tau,Opsi w> by e_psi (TIGHT: = <Opsi tau - qpsi(tau) tau, w> since w perp tau)
and <w,Opsi w> by npsi (Rayleigh quotient on unit w).

At s=0 this must return the EXACT <z, r(tau)> with zero slack (single-point cone).
ACCEPTANCE GATE: feed the TRUE v as "tau" (s=0) and confirm UB_expansion(y) at
m=2 reproduces the true-v UB(y)/R ratios (team's 0.835/0.847/0.899-type numbers,
i.e. matches margin_factor from lemmaC_Rprime_M1.py Section 1), before testing s>0.
Do NOT proceed to CERT-B until this gate passes.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(__file__))
from lemmaC_compact_probe import O_of, top, THETA_STAR


def build(w, N=600):
    x = w * w
    f = lambda k: x / (k * k + x)
    psi = lambda k: 1.0 / (k * k + x) ** 2
    A = O_of(f, N)
    Opsi = O_of(psi, N)
    h, v, l2, ev, EV = top(A)
    g = h - l2
    npsi = np.abs(np.linalg.eigvalsh(Opsi)).max()
    return dict(w=w, x=x, N=N, A=A, Opsi=Opsi, h=h, v=v, l2=l2, g=g, npsi=npsi)


def gap_an(w, h):
    return 0.4 * h if w <= 2 else h - THETA_STAR


def r_bounds(pr, t, s, z):
    """sup/inf of <z, r(v)> over the sinTheta<=s cone around unit t, for
    ARBITRARY z (not necessarily perp t). r(v) = Opsi v - v <v,Opsi v>."""
    Opsi, npsi = pr["Opsi"], pr["npsi"]
    a = float(z @ t)
    z_perp = z - a * t
    nzperp = np.linalg.norm(z_perp)
    Opsi_t = Opsi @ t
    qpsi_t = float(t @ Opsi_t)
    e_psi_vec = Opsi_t - qpsi_t * t
    e_psi = np.linalg.norm(e_psi_vec)
    sq = np.sqrt(max(1.0 - s * s, 0.0))

    cross_perp = float(z_perp @ Opsi_t)
    hi_perp = sq * cross_perp + 2 * s * npsi * nzperp
    lo_perp = sq * cross_perp - 2 * s * npsi * nzperp

    coef = s * abs(2 * s * s - 1.0)
    tau_center = sq * s * s * qpsi_t
    tau_rad = coef * e_psi + sq * s * s * npsi
    hi_tau = tau_center + tau_rad
    lo_tau = tau_center - tau_rad

    if a >= 0:
        hi = a * hi_tau + hi_perp
        lo = a * lo_tau + lo_perp
    else:
        hi = a * lo_tau + hi_perp
        lo = a * hi_tau + lo_perp
    return lo, hi


def ub_expansion(pr, t, s, g_an, y):
    """UB via exact expansion of ||rho||^2 (no triangle), y arbitrary."""
    A, h = pr["A"], pr["h"]
    Ty = h * y - A @ y
    yTy = float(y @ Ty)
    TyTy = float(Ty @ Ty)
    _, B_r_y = r_bounds(pr, t, s, y)       # sup <y, r(v)>  (safe upper)
    lb_rTy, _ = r_bounds(pr, t, s, Ty)     # inf <Ty, r(v)> (safe lower)
    Opsi_t = pr["Opsi"] @ t
    qpsi_t = float(t @ Opsi_t)
    e_psi = np.linalg.norm(Opsi_t - qpsi_t * t)
    E_hi = e_psi if s == 0 else e_psi + s * (abs(qpsi_t) + pr["npsi"])
    return -yTy + 2.0 * B_r_y + E_hi ** 2 / g_an - 2.0 * lb_rTy / g_an + TyTy / g_an


print(f"{'w':>6} {'h':>9} {'g_an':>9} {'margin_factor':>14} {'ratio_dk':>10} "
      f"{'ratio_exp_m1':>13} {'ratio_exp_m2':>13} {'exp_m2/target':>13}")
for w in (0.3, 1.0, 2.0, 3.0, 8.5, 20.0):
    pr = build(w)
    h, v, g = pr["h"], pr["v"], pr["g"]
    g_an = gap_an(w, h)
    t = v.copy()   # s=0 acceptance gate: feed TRUE v as "tau"
    s = 0.0

    Opsi_t = pr["Opsi"] @ t
    qpsi_t = float(t @ Opsi_t)
    r = Opsi_t - qpsi_t * t
    E = np.linalg.norm(r)
    F = None  # not needed for this ratio check; use E^2/g_an vs UB directly
    ratio_dk = E * E / g_an

    y0 = r / E
    ub_m1 = ub_expansion(pr, t, s, g_an, y0)
    # grid/local search over scale c for y=c*y0 (m=1)
    from scipy.optimize import minimize_scalar
    f1 = lambda c: ub_expansion(pr, t, s, g_an, c * y0)
    cs = np.concatenate(([0.0, 1.0], np.geomspace(1e-4, 20.0, 100)))
    vals = [f1(c) for c in cs]
    j = int(np.argmin(vals))
    lo_c, hi_c = cs[max(j - 1, 0)], cs[min(j + 1, len(cs) - 1)]
    if hi_c <= lo_c:
        hi_c = lo_c + 1e-6
    r1 = minimize_scalar(f1, bounds=(lo_c, hi_c + 1e-9), method="bounded")
    ub_m1_best = min(vals[0], vals[j], float(r1.fun))
    ratio_exp_m1 = ub_m1_best / g_an  # this is UB itself already includes /g terms... wait

    # NOTE: ub_expansion already returns a quantity comparable to E_hi^2/g_an
    # (i.e. it IS already-normalized "UB", not UB*g). So ratio vs DK:
    ratio_m1_val = ub_m1_best

    # m=2: second Krylov vector Ty0 Gram-Schmidt against t,y0
    Ty0 = pr["h"] * y0 - pr["A"] @ y0
    y1raw = Ty0 - (Ty0 @ t) * t - (Ty0 @ y0) * y0
    n1 = np.linalg.norm(y1raw)
    ub_m2_best = ub_m1_best
    if n1 > 1e-10:
        y1 = y1raw / n1
        from scipy.optimize import minimize
        f2 = lambda p: ub_expansion(pr, t, s, g_an, p[0] * y0 + p[1] * y1)
        best = ub_m1_best
        bestpt = (float(cs[j]), 0.0)
        for c0 in np.concatenate(([0.0, cs[j]], np.geomspace(1e-4, 10.0, 20))):
            for c1 in np.concatenate(([0.0], np.geomspace(1e-4, 10.0, 15), -np.geomspace(1e-4, 10.0, 15))):
                val = f2((c0, c1))
                if val < best:
                    best = val
                    bestpt = (c0, c1)
        r2 = minimize(f2, x0=np.array(bestpt), method="Nelder-Mead",
                       options=dict(xatol=1e-10, fatol=1e-14, maxiter=1000))
        ub_m2_best = min(best, float(r2.fun))

    # margin_factor from exact top-eigenvector computation (Section-1 style)
    A, x = pr["A"], pr["x"]
    ev_all, EV_all = np.linalg.eigh(A)
    hh = ev_all[-1]
    c = EV_all.T @ r
    c[-1] = 0.0
    R_exact = np.sum(c[:-1] ** 2 / (hh - ev_all[:-1]))
    dk_exact = E * E / g
    margin_factor = dk_exact / R_exact
    target_ratio_true = R_exact  # what UB should approach (>= this, ideally close)

    print(f"{w:6.2f} {h:9.4f} {g_an:9.4f} {margin_factor:14.4f} {ratio_dk:10.4e} "
          f"{ratio_m1_val:13.4e} {ub_m2_best:13.4e} {ub_m2_best/target_ratio_true:13.4f}")
