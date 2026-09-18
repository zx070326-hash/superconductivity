"""[numerical-evidence] Build the ROW-SUM upward pad from zero (do NOT inherit repo M2).

The row sum in Lemma-Lambda's (Stil) is EXACTLY a midpoint sum:

    w_m = 1/(2m+1) = D/(2 s_m),   s_m = (m+1/2) D,  D = 1/varpi
    =>  sum_{m>=0} w_m Rk(s,s_m) = D * sum_{m>=0} Phi_s(s_m),
        Phi_s(sigma) := Rk(s,sigma) / (2 sigma)

with nodes s_m the midpoints of the cells [mD,(m+1)D] tiling (0,inf).  Hence

    | sum_m w_m Rk(s,s_m) - Theta(s) | <= (D^2/8) * V2[Phi_s],
    V2[Phi_s] := int_0^inf |d^2 Phi_s/dsigma^2| dsigma       (= TV of Phi_s')

using the repo's PROVABLE midpoint constant D^2/8 (sec 1.F.R correction), not 1/24.
Theta(s) = int_0^inf Phi_s(sigma) dsigma is checked below, not assumed.

DIRECTION.  The row sum is ADDED, so it needs an UPPER bound -> pad_R is added.
D_nn is SUBTRACTED, so it needs a LOWER bound -> pad_D is also added.  Both pads
push the final bound UP.  Therefore V2[Phi_s] must be UPPER-estimated: a grid TV
converges to the true TV from BELOW, so a naive grid sum would fail open.  We
instead locate every sign change of Phi_s'' and sum |Delta Phi_s'| between
consecutive critical points -- exact TV up to root-location error, with the
sign-change COUNT checked for stability under grid refinement.

Closed form used (elementary algebra, verified against the repo form below):
    (1+(s-t)^2)(1+(s+t)^2) = (1+s^2+t^2)^2 - 4 s^2 t^2
    Phi_s(t) = s^2 t (5+s^2+t^2) / [ ((1+s^2+t^2)^2 - 4 s^2 t^2)(1+s^2)(1+t^2) ]

Run: OMP_NUM_THREADS=1 python _chk_padR_rowsum.py
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import mpmath as mp

mp.mp.dps = 30
V2G = 3.0 * np.sqrt(3.0) / 2.0          # int |G''| dt, exact


def Rk_repo(s, t):
    """Repo form of R(s,sigma), verbatim from P4_einstein_extremal.md:455."""
    num = 2.0 * s * s * t * t * (5.0 + s * s + t * t)
    den = ((1.0 + (s - t) ** 2) * (1.0 + (s + t) ** 2)
           * (1.0 + s * s) * (1.0 + t * t))
    return num / den


def Phi(s, t):
    """Phi_s(t) = Rk(s,t)/(2t), factored denominator (mp or float)."""
    q = 1 + s * s + t * t
    return s * s * t * (5 + s * s + t * t) / ((q * q - 4 * s * s * t * t)
                                              * (1 + s * s) * (1 + t * t))


def theta(s):
    return (np.log1p(s * s) + 2.0 * s * np.arctan(s)) / (2.0 * (1.0 + s * s))


def d1(s, t):
    return mp.diff(lambda x: Phi(mp.mpf(s), x), t)


def d2(s, t):
    return mp.diff(lambda x: Phi(mp.mpf(s), x), t, 2)


def V2_Phi(s, ngrid=4001, tmax=400.0, verbose=False):
    """TV of Phi_s' = int|Phi_s''|, via exact-ish critical points of Phi_s'."""
    ts = np.geomspace(1e-6, tmax, ngrid)
    d2v = np.array([float(d2(s, mp.mpf(float(t)))) for t in ts])
    idx = np.where(np.sign(d2v[:-1]) != np.sign(d2v[1:]))[0]
    roots = []
    for i in idx:
        r = mp.findroot(lambda x: d2(s, x), (mp.mpf(float(ts[i])),
                                             mp.mpf(float(ts[i + 1]))),
                        solver='bisect', tol=mp.mpf(10) ** (-25))
        roots.append(r)
    # Phi_s'(0+) is finite; Phi_s'(inf)=0.  Endpoints: 0 and tmax (tail added).
    nodes = [mp.mpf(10) ** (-12)] + roots + [mp.mpf(tmax)]
    tv = mp.mpf(0)
    for a, b in zip(nodes[:-1], nodes[1:]):
        tv += abs(d1(s, b) - d1(s, a))
    tv += abs(d1(s, mp.mpf(tmax)))          # tail |Phi'(tmax)-Phi'(inf)|, monotone
    if verbose:
        print(f"      s={s:g}: {len(roots)} critical pts of Phi', TV={float(tv):.6g}")
    return float(tv), len(roots)


if __name__ == "__main__":
    print("== 0. sanity: factored Phi == repo Rk/(2t), and int Phi == Theta ==")
    for s in (0.3, 1.0, 2.86906499489, 10.0):
        t = 1.7
        a = Rk_repo(s, t) / (2 * t)
        b = float(Phi(s, t))
        I = float(mp.quad(lambda x: Phi(mp.mpf(s), x), [0, 1, 10, mp.inf]))
        print(f"   s={s:<14g} Phi rel.diff={abs(a - b) / abs(b):.2e}"
              f"   int Phi={I:.12f}  Theta={theta(s):.12f}"
              f"   diff={I - theta(s):+.2e}")

    print("\n== 1. V2[Phi_s] and grid-refinement stability of the root count ==")
    print(f"{'s':>10} {'V2(n=2001)':>12} {'V2(n=4001)':>12} {'V2(n=8001)':>12}"
          f" {'#roots':>7}")
    SS = [0.5, 1.0, 1.5, 2.0, 2.86906499489, 4.0, 6.0, 10.0, 20.0]
    V2 = {}
    for s in SS:
        vs, ks = [], []
        for ng in (2001, 4001, 8001):
            v, k = V2_Phi(s, ngrid=ng)
            vs.append(v); ks.append(k)
        V2[s] = vs[-1]
        ok = "ok" if len(set(ks)) == 1 else f"UNSTABLE {ks}"
        print(f"{s:10.4f} {vs[0]:12.6f} {vs[1]:12.6f} {vs[2]:12.6f} {ks[-1]:5d} {ok}")

    print("\n   scaling check  V2[Phi_s] * s^k  (looking for the large-s power):")
    for s in SS:
        print(f"      s={s:<10.4f} V2={V2[s]:12.6f}  V2*s^2={V2[s] * s * s:10.5f}"
              f"  V2*s^3={V2[s] * s ** 3:10.5f}")
