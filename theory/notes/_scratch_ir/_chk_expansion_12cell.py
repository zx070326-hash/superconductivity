"""s>0 test: 12-cell freed-width trial family, expansion-based (non-triangle)
upper bound on R, per team-lead's round-3 fix. Builds on
_chk_expansion_gate.py's r_bounds()/ub_expansion() (verified there at s=0,
exp_m2/target in [1.0000,1.0015] across w=0.3..20 -- gate PASSED).

Does NOT touch CERT-B (per team-lead's explicit instruction: only look at how
much the leakage term eats into the margin, nothing about cell counts yet).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from lemmaC_compact_probe import O_of, top, THETA_STAR
from scipy.optimize import minimize_scalar, minimize

N_DEFAULT = 600


def build(w, N=N_DEFAULT):
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


def trial_vector(pr, alpha, c=0.5):
    N, w = pr["N"], pr["w"]
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)
    a = alpha * w
    t = u * a * a / (a * a + (n + c) ** 2)
    return t / np.linalg.norm(t)


def sinTheta_bound(pr, t, g_an):
    A = pr["A"]
    rho = t @ (A @ t)
    res = np.linalg.norm(A @ t - rho * t)
    return min(res / g_an, 1.0) if g_an > 0 else np.inf


def r_bounds(pr, t, s, z):
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


def E_hi_of(pr, t, s):
    Opsi_t = pr["Opsi"] @ t
    qpsi_t = float(t @ Opsi_t)
    e_psi = np.linalg.norm(Opsi_t - qpsi_t * t)
    return e_psi if s == 0 else e_psi + s * (abs(qpsi_t) + pr["npsi"])


def F_lo_of(pr, t, s, Ophi, nphi):
    Ft = float(t @ (Ophi @ t))
    e_phi = np.linalg.norm(Ophi @ t - Ft * t)
    return (1 - s * s) * Ft - 2 * s * e_phi - s * s * nphi


def ub_expansion(pr, t, s, g_an, y):
    A, h = pr["A"], pr["h"]
    Ty = h * y - A @ y
    yTy = float(y @ Ty)
    TyTy = float(Ty @ Ty)
    _, B_r_y = r_bounds(pr, t, s, y)
    lb_rTy, _ = r_bounds(pr, t, s, Ty)
    E_hi = E_hi_of(pr, t, s)
    return -yTy + 2.0 * B_r_y + E_hi ** 2 / g_an - 2.0 * lb_rTy / g_an + TyTy / g_an


def best_ub(pr, t, s, g_an, m2=True):
    Opsi_t = pr["Opsi"] @ t
    qpsi_t = float(t @ Opsi_t)
    e_psi_vec = Opsi_t - qpsi_t * t
    e_psi = np.linalg.norm(e_psi_vec)
    E_hi = E_hi_of(pr, t, s)
    dk = E_hi ** 2 / g_an
    if e_psi < 1e-13:
        return dk, dk
    y0 = e_psi_vec / e_psi
    f1 = lambda c: ub_expansion(pr, t, s, g_an, c * y0)
    cs = np.concatenate(([0.0, 1.0], np.geomspace(1e-4, 20.0, 80)))
    vals = [f1(c) for c in cs]
    j = int(np.argmin(vals))
    lo_c, hi_c = cs[max(j - 1, 0)], cs[min(j + 1, len(cs) - 1)]
    if hi_c <= lo_c:
        hi_c = lo_c + 1e-6
    r1 = minimize_scalar(f1, bounds=(lo_c, hi_c + 1e-9), method="bounded")
    ub_m1 = min(dk, vals[0], vals[j], float(r1.fun))
    if not m2:
        return dk, ub_m1
    h = pr["h"]
    Ty0 = h * y0 - pr["A"] @ y0
    y1raw = Ty0 - (Ty0 @ t) * t - (Ty0 @ y0) * y0
    n1 = np.linalg.norm(y1raw)
    ub_m2 = ub_m1
    if n1 > 1e-10:
        y1 = y1raw / n1
        f2 = lambda p: ub_expansion(pr, t, s, g_an, p[0] * y0 + p[1] * y1)
        best = ub_m1
        bestpt = (float(cs[j]), 0.0)
        for c0 in np.concatenate(([0.0, cs[j]], np.geomspace(1e-4, 10.0, 15))):
            for c1 in np.concatenate(([0.0], np.geomspace(1e-4, 10.0, 10), -np.geomspace(1e-4, 10.0, 10))):
                val = f2((c0, c1))
                if val < best:
                    best = val
                    bestpt = (c0, c1)
        r2 = minimize(f2, x0=np.array(bestpt), method="Nelder-Mead",
                       options=dict(xatol=1e-9, fatol=1e-13, maxiter=600))
        ub_m2 = min(best, float(r2.fun))
    return dk, ub_m2


CELLS = [
    (0.30, 0.42, 5.0), (0.42, 0.60, 3.6), (0.60, 0.85, 2.7), (0.85, 1.20, 2.1),
    (1.20, 1.70, 1.8), (1.70, 2.40, 1.6), (2.40, 3.40, 1.5), (3.40, 4.80, 1.42),
    (4.80, 6.80, 1.35), (6.80, 9.60, 1.30), (9.60, 13.6, 1.27), (13.6, 20.0, 1.24),
]

print(f"{'cell':>15} {'alpha':>6} {'worst w':>8} {'s':>7} {'F_lo':>10} {'ratio_dk':>9} {'ratio_exp_m2':>12} status")
worst_all = 0.0

for wa, wb, alpha in CELLS:
    ws = np.geomspace(wa, wb, 11)
    worst = None
    for w in ws:
        pr = build(w)
        g_an = gap_an(w, pr["h"])
        t = trial_vector(pr, alpha)
        s = sinTheta_bound(pr, t, g_an)
        x = pr["x"]
        Ophi = O_of(lambda k: k * k / (k * k + x) ** 3, pr["N"])
        nphi = np.abs(np.linalg.eigvalsh(Ophi)).max()
        F_lo = F_lo_of(pr, t, s, Ophi, nphi)
        if F_lo <= 0:
            continue
        dk_num, ub_num = best_ub(pr, t, s, g_an, m2=True)
        ratio_dk = x * x * dk_num / F_lo
        ratio_exp = x * x * ub_num / F_lo
        if worst is None or ratio_exp > worst[1]:
            worst = (w, ratio_exp, ratio_dk, s, F_lo)
    w_worst, r_exp, r_dk, s_w, Flo_w = worst
    worst_all = max(worst_all, r_exp)
    status = "PASS" if r_exp < 1 else "FAIL"
    print(f"[{wa:5.2f},{wb:5.2f}] {alpha:6.2f} {w_worst:8.3f} {s_w:7.4f} {Flo_w:10.4e} {r_dk:9.4f} {r_exp:12.4f} {status:>6}")

print(f"\nWorst REAL (expansion, s>0) ratio over all 12 cells: {worst_all:.4f}")
print("(reference: DK-only worst was 0.5751; s=0 acceptance gate passed at m=2, "
      "exp/target in [1.0000,1.0015])")
