"""Verify the analytic chain for Lemma (G) and measure the Perron-family barrier.

Chain (all steps claimed analytic):
  (a) |<y,Ay>| <= <|y|,|A| |y|>  for symmetric A   => |lmin(O)| <= lmax(|O|)
  (b) |O[f]| <= P[f] := O[f] + 2 L1[f]  elementwise, P[f] >= 0 elementwise, linear in f
  (c) lmax(P[f]) <= lam1(O[f]) + 2 lmax(L1[f])
  (d) lmax(L1[f]) = sup_n 2S(n)/(2n+1) < f(1)
  (e) lam1(O[f]) >= O[f]_00 = f(1)
  => |lmin| <= lam1 + 2 f(1) <= 3 lam1.

Barrier measurement: sup_f lmax(|O[f]|)/lam1(O[f])  -- if > 1, no elementwise-majorant
(Perron / absolute-value) argument can ever give c <= 1.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _lemG_reduce import build, ev, step


def chain(f, N=900):
    B = build(f, N)
    O = B["O"]
    l1, lmin = ev(O)
    Oabs = np.abs(O)
    labs = ev(Oabs)[0]
    L1 = np.diag(B["L1"])
    P = O + 2.0 * np.diag(L1)
    lP = ev(P)[0]
    lL1 = L1.max()
    return dict(l1=l1, r=-lmin, labs=labs, lP=lP, lL1=lL1, f1=f(1))


def show(name, f, N=900):
    c = chain(f, N)
    ok = (c["r"] <= c["labs"] + 1e-9 <= c["lP"] + 1e-9
          and c["lP"] <= c["l1"] + 2 * c["lL1"] + 1e-9
          and c["lL1"] < c["f1"] + 1e-12
          and c["l1"] >= c["f1"] - 1e-12)
    print(f"{name:22s} r={c['r']:9.6f} <= |O|={c['labs']:9.6f} <= P={c['lP']:9.6f} "
          f"<= l1+2L1={c['l1']+2*c['lL1']:9.6f} <= 3l1={3*c['l1']:9.6f} | "
          f"r/l1={c['r']/c['l1']:7.5f} |O|/l1={c['labs']/c['l1']:7.5f} "
          f"P/l1={c['lP']/c['l1']:7.5f}  CHAIN_OK={ok}")
    return c


if __name__ == "__main__":
    fam = [(f"step K={K}", step(K)) for K in [1, 2, 3, 5, 10, 30, 100, 300]]
    fam += [("1/k^2", lambda k: 1.0 / k ** 2), ("1/k^3", lambda k: 1.0 / k ** 3),
            ("1/k^1.1", lambda k: 1.0 / k ** 1.1),
            ("E w=0.1", lambda k: 0.01 / (k * k + 0.01)),
            ("E w=1", lambda k: 1.0 / (k * k + 1.0)),
            ("E w=7", lambda k: 49.0 / (k * k + 49.0)),
            ("E w=100", lambda k: 1e4 / (k * k + 1e4)),
            ("exp(-k)", lambda k: float(np.exp(-k))),
            ("exp(-k/20)", lambda k: float(np.exp(-k / 20.0)))]
    for name, f in fam:
        show(name, f)

    print("\n--- random nonneg nonincreasing kernels: barrier sup |O|/l1 ---")
    rng = np.random.default_rng(20260918)
    best = (0, None)
    for trial in range(400):
        p = rng.integers(1, 40)
        a = rng.random(p) ** rng.uniform(0.3, 4.0)
        Ks = np.sort(rng.integers(1, 200, size=p))
        fk = np.zeros(2100)
        for aa, K in zip(a, Ks):
            fk[1:K + 1] += aa
        f = lambda k, fk=fk: fk[k] if k < len(fk) else 0.0
        c = chain(f, N=700)
        rat = c["labs"] / c["l1"]
        if rat > best[0]:
            best = (rat, (a.round(4).tolist(), Ks.tolist()))
    print(f"max |O|/l1 over 400 random kernels: {best[0]:.6f}")
