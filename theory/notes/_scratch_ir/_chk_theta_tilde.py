"""[numerical-evidence, with an analytic target] The constant that replaces Theta*.

Repo 1.F.3(c) bounds  lam2(H) <= lam1(R - D) <= ||R|| <= Theta* = 0.63493,
i.e. it throws D away before Schur.  Keeping D:

    lam1(R - D) <= sup_n [ sum_m w_m Rk(s_n,s_m) - D_nn ]          (Lemma Lambda, k=0)

Continuum limits of the two pieces (both are midpoint sums of spacing D=1/varpi):
    sum_m w_m Rk(s_n,s_m)  ->  Theta(s) = [ln(1+s^2) + 2 s arctan s] / (2(1+s^2))
    D_nn = (1 + 2 sum_{k=1}^n G(k/varpi))/(2n+1)  ->  arctan(s)/s
(the latter because 1 + 2 sum_{k=1}^n G(kD) = sum_{k=-n}^{n} G(kD) is the midpoint
 sum of int_{-s_n}^{s_n} G = 2 arctan(s_n), over a cell of total width 2 s_n)

Both behave like pi/(2s) as s -> infinity, so they cancel at leading order:

    Thetatil(s) := Theta(s) - arctan(s)/s

This script: locates sup_s Thetatil(s), checks Thetatil(0+) = -1 and the large-s
tail ~ ln(s)/s^2, and checks the discrete sup_n against it over varpi.

Run: OMP_NUM_THREADS=1 python _chk_theta_tilde.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import mpmath as mp

MT, NROW = 200000, 400


def G(y):
    return 1.0 / (1.0 + y * y)


def Rk(s, sg):
    num = 2.0 * s * s * sg * sg * (5.0 + s * s + sg * sg)
    den = ((1.0 + (s - sg) ** 2) * (1.0 + (s + sg) ** 2)
           * (1.0 + s * s) * (1.0 + sg * sg))
    return num / den


def theta_tilde_mp(s):
    s = mp.mpf(s)
    return ((mp.log(1 + s ** 2) + 2 * s * mp.atan(s)) / (2 * (1 + s ** 2))
            - mp.atan(s) / s)


def discrete_sup(w):
    """sup_n [ sum_{m>=0} w_m Rk(s_n,s_m) - D_nn ] with an upward tail pad."""
    x = w * w
    m = np.arange(MT, dtype=float)
    sm = (m + 0.5) / w
    wm = 1.0 / (2.0 * m + 1.0)
    cumG = np.concatenate(([0.0], np.cumsum(G(np.arange(1, NROW + 1) / w))))
    best, arg = -np.inf, -1
    for a in range(0, NROW, 20):
        b = min(a + 20, NROW)
        nn = np.arange(a, b, dtype=float)
        sn = (nn + 0.5) / w
        row = Rk(sn[:, None], sm[None, :]) @ wm
        pad = 2.0 * x / ((2.0 * MT + 1.0) * (MT - nn - 1.0))
        Dnn = (1.0 + 2.0 * cumG[a:b]) / (2.0 * nn + 1.0)
        v = row + pad - Dnn
        if v.max() > best:
            best, arg = float(v.max()), a + int(v.argmax())
    return best, arg


if __name__ == "__main__":
    mp.mp.dps = 40
    print("== Thetatil(s) = Theta(s) - arctan(s)/s ==")
    for s in (1e-3, 0.1, 0.5, 1.0, 1.500904, 2.0, 3.0, 5.0, 10.0, 100.0, 1e4):
        print(f"   s={s:10.4g}   Thetatil={float(theta_tilde_mp(s)):+.8f}")

    # locate the max
    g = mp.findroot(lambda s: mp.diff(theta_tilde_mp, s), mp.mpf('3.0'))
    v = theta_tilde_mp(g)
    print(f"\n   stationary point  s~ = {mp.nstr(g, 12)}")
    print(f"   Thetatil(s~)         = {mp.nstr(v, 12)}")
    # brute-force confirm it is the global max
    ss = np.geomspace(1e-4, 1e6, 20001)
    vv = np.array([float(theta_tilde_mp(s)) for s in ss[::40]])
    print(f"   coarse global scan over s in [1e-4,1e6]: max = {vv.max():.10f}")

    print(f"\n== discrete sup_n vs Thetatil* = {float(v):.8f} ==")
    print(f"{'varpi':>7} {'discrete sup':>13} {'argmax n':>9} {'s_argmax':>9}"
          f" {'/Thetatil*':>11}")
    for w in (0.30, 0.42, 0.60, 0.85, 1.20, 1.70, 2.40, 3.40, 4.80,
              6.80, 9.60, 13.60, 20.00, 45.0):
        d, n = discrete_sup(w)
        print(f"{w:7.2f} {d:13.6f} {n:9d} {(n + 0.5) / w:9.4f}"
              f" {d / float(v):11.4f}")
