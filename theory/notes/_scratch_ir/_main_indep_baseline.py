"""MAIN-SESSION INDEPENDENT BASELINE — written from the derivation text only.

Source of truth for the operator: theory/derivations/P4_main_theorems.md line 17 / §0
    u_n = (2n+1)^{-1/2}
    O[f] = -L1[f] + L2[f] + L3[f]
    L1  diagonal   : (2n+1)^{-1} * sum_{k=1}^{n} 2 f(k)        <-- upper limit INCLUDES k=n
    L2  off-diag   : f(|n-m|) u_n u_m
    L3  (all n,m)  : f(n+m+1) u_n u_m
    F_T(k) = int [[k]](om/2piT) P(dom),  [[k]](vp) = vp^2/(k^2+vp^2),  F_T(0)=1
    P(dom) = (2/lam) a2F(om)/om dom,  lam = 2 int a2F/om dom
    Tc criterion: lam * lam1(O[F_T]) = 1 ;  Coulomb: M = O[F_T] - 2(mu*/lam) u_c u_c^T,
    u_c[n] = u_n * 1[n < N_c],  N_c = floor(r * vp),  vp = sqrt(<om^2>)/(2 pi T)

Deliberately NOT importing anything from the existing scratch scripts: this exists to be
an EXTERNAL baseline for h3s_deficit_cert_v2.py, whose own "alignment check" compares
interval arithmetic against float64 code living in the same file (that check cannot
detect a wrong operator -- 2026-09-17 lesson).

[numerical-evidence] float64 throughout.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import math
import numpy as np

MEV_K = 11.604518
A2F = ("D:/physics/superconductivity/theory/data/external/mishra_h3s_d3s/"
       "Data-H3S-D3S-PH-EPW-Anharmonic/H3S-200GPa/EPW_anharmonic/H3S.a2f")


def read_a2f(path):
    """Return w[meV], a2f, and the file's OWN cumulative column 3 (= int 2 a2F/w)."""
    w, a2f, cum = [], [], []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                x, y, c = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue          # header line
            if x > 0.0:
                w.append(x); a2f.append(y); cum.append(c)
    return np.array(w), np.array(a2f), np.array(cum)


def trapz_weights(x):
    """Plain trapezoid weights, written out rather than calling np.trapezoid."""
    tw = np.empty_like(x)
    tw[0] = 0.5 * (x[1] - x[0])
    tw[-1] = 0.5 * (x[-1] - x[-2])
    tw[1:-1] = 0.5 * (x[2:] - x[:-2])
    return tw


def P_measure(w, a2f):
    tw = trapz_weights(w)
    lam = 2.0 * np.sum((a2f / w) * tw)
    p = (2.0 / lam) * (a2f / w) * tw          # already sums to 1 by construction
    om2 = np.sum(p * w ** 2)
    return lam, p, om2


def F_T(kmax, w, p, T_mev):
    """F_T(k) for k = 0..kmax, F_T(0) := 1."""
    vp = w / (2.0 * math.pi * T_mev)
    k = np.arange(kmax + 1, dtype=float)
    # [[k]] = vp^2/(k^2+vp^2), integrated against p
    F = (p[None, :] * (vp[None, :] ** 2 / (k[:, None] ** 2 + vp[None, :] ** 2))).sum(axis=1)
    F[0] = 1.0
    return F


def O_of_f(F, N):
    """Build O[f] as an explicit N x N matrix, straight from the definition."""
    n = np.arange(N)
    u = 1.0 / np.sqrt(2.0 * n + 1.0)

    # L1[n] = (1/(2n+1)) * sum_{k=1}^{n} 2 F(k)   -- explicit loop, no cumsum tricks
    L1 = np.zeros(N)
    acc = 0.0
    for nn in range(N):
        if nn >= 1:
            acc += 2.0 * F[nn]                 # adds k = nn, so acc = sum_{k=1}^{nn} 2F(k)
        L1[nn] = acc / (2.0 * nn + 1.0)

    M = np.empty((N, N))
    for i in range(N):
        for j in range(N):
            val = F[abs(i - j)] * u[i] * u[j] if i != j else 0.0   # L2: off-diagonal only
            val += F[i + j + 1] * u[i] * u[j]                      # L3
            if i == j:
                val -= L1[i]
            M[i, j] = val
    return n, u, L1, M


def lam1_M(T_mev, w, p, lam, om2, mus, r, N):
    F = F_T(2 * N, w, p, T_mev)
    n, u, L1, A = O_of_f(F, N)
    if mus == 0.0:
        return float(np.linalg.eigvalsh(A)[-1])
    vp = math.sqrt(om2) / (2.0 * math.pi * T_mev)
    Nc = int(math.floor(r * vp))
    uc = np.where(n < Nc, u, 0.0)
    M = A - 2.0 * (mus / lam) * np.outer(uc, uc)
    return float(np.linalg.eigvalsh(M)[-1])


def Tc(w, p, lam, om2, mus, r, N, lo=0.1, hi=60.0, tol=1e-10):
    """Bisect lam*lam1(T) - 1 = 0 in T [meV]; g decreasing in T."""
    g = lambda T: lam * lam1_M(T, w, p, lam, om2, mus, r, N) - 1.0
    glo, ghi = g(lo), g(hi)
    if glo < 0 or ghi > 0:
        raise RuntimeError(f"bracket fails: g({lo})={glo:.3e}, g({hi})={ghi:.3e}")
    while hi - lo > tol * hi:
        mid = 0.5 * (lo + hi)
        if g(mid) >= 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    w, a2f, cum = read_a2f(A2F)
    lam, p, om2 = P_measure(w, a2f)
    print("=== independent baseline (main session), H3S-200GPa anharmonic ===")
    print(f"grid points        : {len(w)}   w in [{w[0]:.4f}, {w[-1]:.4f}] meV")
    print(f"lambda (my trapz)  : {lam:.12f}")
    print(f"lambda (file col3) : {cum[-1]:.12f}   <-- file's own cumulative integral")
    print(f"  reldiff          : {abs(lam - cum[-1]) / lam:.3e}")
    print(f"sqrt<om^2>         : {math.sqrt(om2):.6f} meV")
    print(f"sum p              : {p.sum():.16f}")
    print()

    # Element-level values for cross-checking cert_v2 (T = 1.0 meV, the script's test point)
    N = 400
    F = F_T(2 * N, w, p, 1.0)
    n, u, L1, A = O_of_f(F, N)
    print("T = 1.0 meV, N = 400 -- element-level baseline for cert_v2 comparison")
    for idx in (0, 1, 20, 2 * N):
        print(f"  F[{idx:3d}]   = {F[idx]:.16e}")
    for idx in (0, 1, 20, N - 1):
        print(f"  L1[{idx:3d}]  = {L1[idx]:.16e}")
    for (i, j) in ((0, 0), (0, 1), (20, 20), (20, 21)):
        print(f"  A[{i:2d},{j:2d}] = {A[i, j]:.16e}")
    print()

    for N in (400, 800):
        t0 = Tc(w, p, lam, om2, 0.0, 10.0, N)
        t13 = Tc(w, p, lam, om2, 0.13, 10.0, N)
        vp = math.sqrt(om2) / (2.0 * math.pi * t13)
        print(f"N={N:4d}  Tc(mu*=0)   = {t0 * MEV_K:9.4f} K   ({t0:.6f} meV)")
        print(f"N={N:4d}  Tc(mu*=.13) = {t13 * MEV_K:9.4f} K   ({t13:.6f} meV)"
              f"   varpi={vp:.4f}  Nc=floor(10*varpi)={int(math.floor(10 * vp))}"
              f"   deficit={100 * (1 - t13 / t0):.2f}%", flush=True)
