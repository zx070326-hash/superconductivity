"""EXACT rational evaluation of the explicit analytic bound

    rho = |lmin(O[f])| / lam1(O[f])  <=  1 + 2*Gamma,
    Gamma := sup_K  theta_K / q_K(v),   theta_K = 2K/(2K+1),
    q_K(v) := <v, O[1_{k<=K}] v>/||v||^2  with  v_n = u_n t_n (t rational, finite support,
    v nonneg nonincreasing).

q_K is rational; q_K nondecreasing in K (Cor. A'); theta_K < 1  =>
    Gamma <= max( max_{K<=K0} theta_K/q_K ,  1/q_{K0+1} )   -- a FINITE check.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
from fractions import Fraction as F


def qK_exact(t, K):
    """t: list of Fractions, v_n = u_n t_n, n < p.  Returns q_K(v) as Fraction."""
    p = len(t)
    us = [F(1, 2 * n + 1) for n in range(p)]        # u_n^2
    nrm = sum(us[n] * t[n] * t[n] for n in range(p))
    tot = F(0)
    for n in range(p):
        S = min(n, K)                                # sum_{k<=n} 1_{k<=K}
        diag = -2 * S * us[n]                        # L1 only; L3 diagonal is in the
        tot += diag * us[n] * t[n] * t[n]            # third loop below (n+m+1=2n+1)
        for m in range(p):
            if m == n:
                continue
            if abs(n - m) <= K:
                tot += us[n] * us[m] * t[n] * t[m]
        for m in range(p):
            if n + m + 1 <= K:
                tot += us[n] * us[m] * t[n] * t[m]
    return tot / nrm


def qinf_exact(t):
    p = len(t)
    us = [F(1, 2 * n + 1) for n in range(p)]
    nrm = sum(us[n] * t[n] * t[n] for n in range(p))
    s = sum(us[n] * t[n] for n in range(p))          # <u,v> = sum u_n^2 t_n
    return 2 * s * s / nrm - 1


def gamma_bound(t, K0=400):
    qs = {K: qK_exact(t, K) for K in range(1, K0 + 2)}
    # monotonicity sanity (Cor A'): q_K nondecreasing
    mono = all(qs[K + 1] >= qs[K] for K in range(1, K0 + 1))
    worst, argw = F(0), None
    for K in range(1, K0 + 1):
        th = F(2 * K, 2 * K + 1)
        val = th / qs[K]
        if val > worst:
            worst, argw = val, K
    tail = 1 / qs[K0 + 1]
    return worst, argw, tail, max(worst, tail), mono, qs


if __name__ == "__main__":
    # cross-check q_K against the float operator from _lemG_reduce
    import sys, numpy as np
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from _lemG_reduce import build, step
    for p in (1, 2, 5):
        t = [F(1)] * p
        for K in (1, 2, 5):
            B = build(step(K), 400)
            u = B["u"][:p]
            v = np.array([float(u[n]) * 1.0 for n in range(p)])
            vf = np.zeros(400); vf[:p] = v
            q_float = float(vf @ B["O"] @ vf) / float(vf @ vf)
            q_rat = float(qK_exact(t, K))
            assert abs(q_float - q_rat) < 1e-12, (p, K, q_float, q_rat)
    print("q_K exact-vs-float cross-check PASSED\n")

    best = None
    for p in range(1, 15):
        t = [F(1)] * p
        worst, argw, tail, G, mono, qs = gamma_bound(t, K0=300)
        c = 1 + 2 * G
        if best is None or c < best[0]:
            best = (c, p, G, argw)
        print(f"v=u^({p:2d})  q_1={float(qs[1]):.6f} q_301={float(qs[301]):.6f} "
              f"q_inf={float(qinf_exact(t)):.6f}  maxK<=300: {float(worst):.6f}@K={argw} "
              f" tail={float(tail):.6f}  A<={float(G):.6f}  c = 1+2A <= {float(c):.6f} "
              f" mono={mono}")
    print(f"\nbest flat-t: c <= {float(best[0]):.6f} at p={best[1]}")

    # non-flat trial vectors t (integers), local search
    import itertools
    print("\ninteger t search (p<=6, entries 0..9):")
    bestc, bestt = best[0], [F(1)] * best[1]
    for p in range(2, 7):
        for tt in itertools.product(range(0, 10), repeat=p - 1):
            t = [F(9)] + [F(x) for x in tt]
            if any(t[i] == 0 for i in range(p)):
                continue
            # v must be nonneg nonincreasing: u_n t_n nonincreasing
            vv = [F(t[n], 1) for n in range(p)]
            ok = all((vv[n] ** 2) * F(1, 2 * n + 1) >= (vv[n + 1] ** 2) * F(1, 2 * n + 3)
                     for n in range(p - 1))
            if not ok:
                continue
            w, aw, tl, G, mono, qs = gamma_bound(t, K0=120)
            c = 1 + 2 * G
            if c < bestc:
                bestc, bestt = c, t
    print(f"best integer t = {[str(x) for x in bestt]}  ->  c <= {float(bestc):.6f}")
    w, aw, tl, G, mono, qs = gamma_bound(bestt, K0=300)
    print(f"  recheck K0=300: A <= {G} = {float(G):.9f},  c <= {float(1+2*G):.9f}, mono={mono}")
    print(f"  worst finite K = {aw}, tail 1/q_301 = {float(tl):.9f}")
