"""
Greedy cell partition of x in [x_lo, x_hi] for a minimal exact-rational
Lemma C certificate:  one rational trial vector (truncated top eigenvector,
support N_s, components rounded to multiples of 2^-12) per cell, analytic gap
lower bound  g_an = 0.4*rho(t)  (w<=2, trace-bound route) or rho(t) - Theta*
(w>=2), DK propagation as in Sec 1.F.3(e).  Cell accepted if the DK ratio
x^2 E_hi^2/(g_an F_lo) < RATIO_MAX at 9 interior sample points (uniformity is
only *sampled* here -- a real proof needs monotone/Lipschitz envelopes, which
cost maybe another factor 2 in cell count).

float64, numerical-evidence.  Purpose: size estimate only.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import numpy as np
from lemmaC_compact_probe import O_of, top, THETA_STAR
from lemmaC_cert_size import pieces

RATIO_MAX = 0.9
N = 600
DEN = 2.0 ** 12


def gap_an(rho, w):
    return 0.4 * rho if w <= 2 else rho - THETA_STAR


def dk_ratio_t(t, x):
    A, Op, Os, h, v, nphi, npsi = pieces(x, N)
    t = t / np.linalg.norm(t)
    rho = t @ A @ t
    res = np.linalg.norm(A @ t - rho * t)
    glo = gap_an(rho, np.sqrt(x))
    if glo <= 0:
        return np.inf
    s = min(res / glo, 1.0)
    Ft = t @ Op @ t
    e_phi = np.linalg.norm(Op @ t - Ft * t)
    Kt = Os @ t
    qpsi = t @ Kt
    e_psi = np.linalg.norm(Kt - qpsi * t)
    F_lo = (1 - s * s) * Ft - 2 * s * e_phi - s * s * nphi
    E_hi = e_psi + s * (abs(qpsi) + npsi)
    if F_lo <= 0:
        return np.inf
    return x * x * E_hi ** 2 / (glo * F_lo)


def trial_at(xc, Ns):
    A, Op, Os, h, v, nphi, npsi = pieces(xc, N)
    t = np.zeros(N)
    t[:Ns] = np.round(v[:Ns] * DEN) / DEN
    return t


def cell_ok(t, xa, xb, npts=9):
    xs = np.exp(np.linspace(np.log(xa), np.log(xb), npts))
    return all(dk_ratio_t(t, x) < RATIO_MAX for x in xs)


if __name__ == "__main__":
    x_lo, x_hi = 0.09, 8.5 ** 2
    xa = x_lo
    cells = []
    while xa < x_hi:
        # try growing cell widths (geometric), pick smallest Ns that works for the widest cell
        best = None
        for fac in [4.0, 3.0, 2.5, 2.0, 1.7, 1.5, 1.3, 1.2, 1.1]:
            xb = min(xa * fac, x_hi)
            xc = np.sqrt(xa * xb)
            for Ns in [4, 6, 8, 12, 16, 24, 32, 48, 64]:
                t = trial_at(xc, Ns)
                if cell_ok(t, xa, xb):
                    best = (xb, Ns, fac)
                    break
            if best:
                break
        if best is None:
            print(f"FAILED to close a cell starting at x={xa:.4g}")
            sys.exit(1)
        xb, Ns, fac = best
        cells.append((xa, xb, Ns))
        print(f"cell [{xa:8.4f}, {xb:8.4f}]  (w in [{np.sqrt(xa):.3f},{np.sqrt(xb):.3f}])  Ns={Ns:3d}  fac={fac}")
        xa = xb
    print(f"\nTOTAL cells on x in [{x_lo},{x_hi}] (w in [0.3, 8.5]): {len(cells)};  max Ns = {max(c[2] for c in cells)};  total rationals = {sum(c[2] for c in cells)}")
