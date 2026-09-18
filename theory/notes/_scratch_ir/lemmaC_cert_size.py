"""
Size of the smallest 'exact-rational' (CERT-B style) certificate for Lemma C on
(0, 8.5]: for each w, what residual must a trial vector achieve, what support
N_s does that need, and how wide a cell can one trial vector cover?

Assumptions (all optimistic for the certificate, hence a LOWER bound on size):
  * analytic gap lower bound g_an = 0.4 h (w<=2, trace-bound route) or
    h - Theta* (w>=2); h taken at its true value (a real proof would use h_lo<h);
  * all trial-vector quantities (F_t, e_phi, e_psi, norms) at their float64
    values (a real proof would need tail majorants on top);
  * DK propagation exactly as in Sec 1.F.3(e) of P4_einstein_extremal.md.
float64, numerical-evidence.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from lemmaC_compact_probe import O_of, top, THETA_STAR


def pieces(x, N):
    f = lambda k: x / (k * k + x)
    phi = lambda k: k * k / (k * k + x) ** 3
    psi = lambda k: 1.0 / (k * k + x) ** 2
    A, Op, Os = O_of(f, N), O_of(phi, N), O_of(psi, N)
    h, v, l2, _, _ = top(A)
    nphi = np.abs(np.linalg.eigvalsh(Op)).max()
    npsi = np.abs(np.linalg.eigvalsh(Os)).max()
    return A, Op, Os, h, v, nphi, npsi


def dk_ratio(t, A, Op, Os, x, glo, nphi, npsi):
    t = t / np.linalg.norm(t)
    rho = t @ A @ t
    res = np.linalg.norm(A @ t - rho * t)
    s = min(res / glo, 1.0) if glo > 0 else 1.0
    Ft = t @ Op @ t
    e_phi = np.linalg.norm(Op @ t - Ft * t)
    Kt = Os @ t
    qpsi = t @ Kt
    e_psi = np.linalg.norm(Kt - qpsi * t)
    F_lo = (1 - s * s) * Ft - 2 * s * e_phi - s * s * nphi
    E_hi = e_psi + s * (abs(qpsi) + npsi)
    if F_lo <= 0 or glo <= 0:
        return np.inf, res, s
    return x * x * E_hi ** 2 / (glo * F_lo), res, s


if __name__ == "__main__":
    N = 600
    ws = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.5]
    print(f"{'w':>5} {'h':>7} {'g_an':>7} {'g_an/g':>7} | {'Ns':>4} {'res':>8} {'s_bnd':>7} {'ratio':>7} | {'cell half-width (rel x)':>24} {'cells/e-fold':>13}")
    for w in ws:
        x = w * w
        A, Op, Os, h, v, nphi, npsi = pieces(x, N)
        g_an = 0.4 * h if w <= 2 else h - THETA_STAR
        # minimal support N_s with ratio < 0.5
        best = None
        for Ns in [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256]:
            t = np.zeros(N); t[:Ns] = v[:Ns]
            r, res, s = dk_ratio(t, A, Op, Os, x, g_an, nphi, npsi)
            if r < 0.5:
                best = (Ns, res, s, r); break
        if best is None:
            t = v.copy()
            r, res, s = dk_ratio(t, A, Op, Os, x, g_an, nphi, npsi)
            best = (N, res, s, r)
        Ns, res, s, r = best
        t = np.zeros(N); t[:Ns] = v[:Ns]
        # cell half-width: reuse t at x' = x(1+d); find largest d with ratio<1 both sides
        dmax = 0.0
        for d in [0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5]:
            ok = True
            for xp in (x * (1 + d), x * (1 - d)):
                A2, Op2, Os2, h2, v2, nphi2, npsi2 = pieces(xp, N)
                g2 = 0.4 * h2 if np.sqrt(xp) <= 2 else h2 - THETA_STAR
                r2, _, _ = dk_ratio(t, A2, Op2, Os2, xp, g2, nphi2, npsi2)
                if not (r2 < 1.0):
                    ok = False; break
            if ok:
                dmax = d
            else:
                break
        cells_per_efold = (1.0 / (2 * dmax)) if dmax > 0 else np.inf
        print(f"{w:5.2f} {h:7.4f} {g_an:7.4f} {g_an/h:7.3f} | {Ns:4d} {res:8.2e} {s:7.4f} {r:7.4f} | {dmax:24.3f} {cells_per_efold:13.1f}")
