"""Scratch: independent check of Corollary S.1  T_c(P)/T_c^E(<w>) = 1 + a2 s^2 + a3 g3 s^3 + O(s^4)
with a2 = kappa/(2L), a3 = c3/L, L = wbar h'/h, on the true (N=100) operator, wbar_E in {0.5,1,2}."""
import os, sys
import numpy as np
from scipy.optimize import brentq
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "derivations")
sys.path.insert(0, os.path.abspath(ROOT))
from p4_interchange_uniformity import H_mat  # noqa: E402

N = 100


def h(w):
    return np.linalg.eigvalsh(H_mat(N, w, 0))[-1]


def kappa_c3(w):
    E, V = np.linalg.eigh(H_mat(N, w, 0)); v = V[:, -1]; hv = E[-1]
    q2 = v @ H_mat(N, w, 2) @ v; q3 = v @ H_mat(N, w, 3) @ v
    # first-order perturbation: k(P) = h + <v, sum_j H^(j) w^j <xi^j>/j! v> + ...
    # kappa := w^2 q2 / h  (so that h[1 + 1/2 kappa s^2]);  c3 := w^3 q3 / (6 h)
    return w * w * q2 / hv, w ** 3 * q3 / (6 * hv), hv


def kP(P, T):  # P = [(weight, omega)], T in units where wbar = <w>/(2 pi T)
    return np.linalg.eigvalsh(sum(p * H_mat(N, om / (2 * np.pi * T), 0) for p, om in P))[-1]


def Tc(P, lam):
    return brentq(lambda T: lam * kP(P, T) - 1.0, 1e-3, 1e3, xtol=1e-13)


print(f"{'wbarE':>6} {'kappa':>8} {'c3':>8} {'L':>7} {'a2':>8} {'a3':>8} | {'sigma':>6} {'g3':>5} {'(Tc/TcE-1)/s^2':>15} {'pred a2+a3 g3 s':>16} {'resid/s^4':>10}")
for wbarE in (0.5, 1.0, 2.0):
    om_mean = 1.0
    T_E = om_mean / (2 * np.pi * wbarE)
    lam = 1.0 / h(wbarE)
    k2, c3, hv = kappa_c3(wbarE)
    dw = 1e-4
    L = wbarE * (h(wbarE + dw) - h(wbarE - dw)) / (2 * dw) / hv
    a2, a3 = k2 / (2 * L), c3 / L
    for sig in (0.05, 0.1, 0.2):
        for g3 in (0.0, 1.0):
            # 2-atom spectrum with mean 1, rel. std sig, skewness g3: weights p,1-p at xi1<0<xi2
            if g3 == 0.0:
                P = [(0.5, 1 - sig), (0.5, 1 + sig)]
            else:
                # xi = sig*x, two-point with mean0 var1 skew g3: x1=-(sqrt(g3^2+4)+g3)/2... use standard construction
                r = (g3 + np.sqrt(g3 * g3 + 4)) / 2  # x2 = r, x1 = -1/r, p1 = r^2/(1+r^2)
                x1, x2 = -1 / r, r
                p1 = r * r / (1 + r * r)
                P = [(p1, 1 + sig * x1), (1 - p1, 1 + sig * x2)]
            tc = Tc(P, lam)
            ratio = (tc / T_E - 1) / sig ** 2
            pred = a2 + a3 * g3 * sig
            resid = (tc / T_E - 1 - a2 * sig ** 2 - a3 * g3 * sig ** 3) / sig ** 4
            print(f"{wbarE:>6.2f} {k2:>8.4f} {c3:>8.4f} {L:>7.4f} {a2:>8.4f} {a3:>8.4f} | {sig:>6.2f} {g3:>5.1f} {ratio:>15.5f} {pred:>16.5f} {resid:>10.3f}")
