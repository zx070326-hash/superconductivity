"""Probe: sharp (floored) vs smooth Coulomb cutoff -> monotonicity of lam1(M) in w.
M = K_ph(w) - mu*K_c,  K_c = 2 u u^T,  u_n = chi(t_n)/sqrt(2n+1), t_n=(2n+1)/(2 r w).
Sharp chi = indicator reproduces the step model N_c=floor(r*w)-ish.
NOT COMMITTED - scratch probe."""
import numpy as np, scipy.linalg as sla, nx_core as C


def chi_sharp(t):
    return (t < 1.0).astype(float)


def chi_p(p):
    def f(t):
        return 1.0 / (1.0 + t ** (2 * p))
    return f


def dchi_p(p):
    def f(t):
        return -2 * p * t ** (2 * p - 1) / (1.0 + t ** (2 * p)) ** 2
    return f


def Nsize(w, r, cap=1200):
    return int(min(cap, max(300, 4 * max(w, r * w))))


def lam1(w, r, mu, chi, N=None, vec=False):
    N = N or Nsize(w, r)
    n = np.arange(N)
    A = C.Op(C.kern(w * w, N, 0), N)
    t = (2 * n + 1.0) / (2 * r * w)
    u = chi(t) / np.sqrt(2 * n + 1.0)
    A = A - mu * 2 * np.outer(u, u)
    if not vec:
        return sla.eigh(A, eigvals_only=True, subset_by_index=[N - 1, N - 1])[0]
    ev, V = sla.eigh(A, subset_by_index=[N - 2, N - 1])
    v = V[:, 1]
    if v.sum() < 0:
        v = -v
    return ev[1], ev[0], v, u, t, n, N


def margin(w, r, mu, chi, dchi):
    """d lam1/dw = <v,dO/dw v> - mu*4<v,u><v,du/dw>. Both terms >= 0 for chi decreasing."""
    l1, l2, v, u, t, n, N = lam1(w, r, mu, chi, vec=True)
    Ap = C.Op(C.kern(w * w, N, 1), N) * (2 * w)
    Dph = v @ (Ap @ v)
    up = -(t / w) * dchi(t) / np.sqrt(2 * n + 1.0)
    Dc = 4 * (v @ u) * (v @ up)
    return dict(w=w, l1=l1, gap=l1 - l2, Dph=Dph, Dc=Dc, dl=Dph - mu * Dc,
                cone=float(np.min(v[:-1] - v[1:])), vmin=float(v.min()),
                ratio=(mu * Dc / Dph if Dph > 0 else np.inf))


def scan_mono(r, mu, chi, ws):
    vals = np.array([lam1(w, r, mu, chi) for w in ws])
    d = np.diff(vals)
    bad = np.where(d < 0)[0]
    return vals, bad, (d[bad].min() if len(bad) else 0.0)
