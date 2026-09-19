"""Theorem SC (strong-coupling window closure) + Theorem TA (two-atom reduction):
numerical verification of the arithmetic and of the structural claims.

Theorem SC.  N_c(T) = floor(omega_c / 2 pi T).  If  lambda * h(1/r) > 1  with
r = omega_c / omega_rms, then Tc(lambda,P;mu*) <= Tc(lambda,E;mu*) for every P
with the given <omega^2> and every mu* >= 0.  Sufficient explicit form:
lambda > 1 + r^2, because h(w) >= H(w)_{00} = w^2/(1+w^2).

Theorem TA.  sup over all P with <w^2>=xbar of lam1(M_P) equals the sup over
two-atom P with the same second moment.  (Proof is elementary; this script only
stress-tests it: 3- and 4-atom spectra must never beat the best two-atom one.)

float64, finite N  ==> numerical-evidence for every number printed here.
"""

import importlib.util
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")

import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "cr_op", os.path.join(_HERE, "_chk_mustar_conerepair_op_2026-09-19.py")
)
op = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(op)

CHK = [0, 0]


def chk(name, cond, extra=""):
    CHK[0] += 1
    if not cond:
        CHK[1] += 1
    print(f"  [{'ok ' if cond else 'FAIL'}] {name} {extra}")


def h_einstein(w, N=None):
    N = N or int(max(400, 60 * w))
    return np.linalg.eigvalsh(op.O_of_kernel(N, op.kernel_einstein(N, w)))[-1]


# ------------------------------------------------- Theorem SC: the constant ---
def sc_table():
    print("\n### Theorem SC: lambda threshold  lam*(r) = 1/h(1/r)  vs the explicit bound 1+r^2")
    print(f"{'r':>6} {'w=1/r':>8} {'h(1/r)':>14} {'lam*(r)':>12} {'1+r^2':>10} {'ratio':>8} {'g2/r^2':>10}")
    g2 = 1.3181402117626  # repo CERT-D lower end, used only as a reference value
    for r in (2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 20.0):
        w = 1.0 / r
        hv = h_einstein(w)
        lam = 1.0 / hv
        bound = 1.0 + r * r
        chk(f"h(1/r) >= w^2/(1+w^2) at r={r}", hv >= w * w / (1 + w * w) - 1e-14,
            f"({hv:.6e} vs {w*w/(1+w*w):.6e})")
        chk(f"1+r^2 >= lam*(r) at r={r}", bound >= lam, f"({bound:.4f} vs {lam:.4f})")
        print(f"{r:6.1f} {w:8.4f} {hv:14.8e} {lam:12.4f} {bound:10.2f} {bound/lam:8.4f} {r*r/g2:10.4f}")
    print("  (g2/r^2 column: the asymptotic lam*(r) -> r^2/g2 since h(w) ~ g2 w^2 as w->0.)")


def sc_case_check():
    """Direct check of the two cases in the proof of Theorem SC, at N_c = 0."""
    print("\n### Theorem SC: case (ii) ingredient -- Theorem 2 at N_c=0 (mu* invisible)")
    rng = np.random.default_rng(7)
    worst = -np.inf
    for _ in range(40):
        w_rms = float(rng.uniform(0.2, 3.0))
        k = int(rng.integers(2, 5))
        xs = rng.uniform(0.05, 25.0, size=k)
        ps = rng.dirichlet(np.ones(k))
        s = np.sqrt(w_rms**2 / float(ps @ (xs**2)))
        ws = xs * s
        N = int(max(300, 24 * ws.max()))
        A = op.O_of_kernel(N, op.kernel_measure(N, ws, ps))
        H = op.O_of_kernel(N, op.kernel_einstein(N, w_rms))
        gap = np.linalg.eigvalsh(H)[-1] - np.linalg.eigvalsh(A)[-1]
        worst = max(worst, -gap)
    chk("Theorem 2 (mu*=0): lam1(O[F_P]) <= h(w_rms), 40 random P", worst <= 1e-12,
        f"(worst violation {worst:.3e})")


# ------------------------------------------ Theorem TA: constructive 2-atom ---
def struct_coeffs(v, kmax):
    """c_k with  <v, O[f] v> = sum_{k>=1} c_k f(k)   (O linear in the kernel)."""
    N = len(v)
    u = op.u_vec(N)
    y = v * u
    c = np.zeros(kmax + 1)
    # L2: off-diagonal Toeplitz, |n-m| = k >= 1, counted twice (n<m and n>m)
    for k in range(1, min(kmax, N - 1) + 1):
        c[k] += 2.0 * float(y[:-k] @ y[k:])
    # L3: Hankel, n+m+1 = k
    idx = np.add.outer(np.arange(N), np.arange(N)) + 1
    outer = np.outer(y, y)
    np.add.at(c, np.clip(idx.ravel(), 0, kmax), outer.ravel() * (idx.ravel() <= kmax))
    # -L1: diagonal  -2 u_n^2 sum_{k=1}^{n} f(k)  ==>  coefficient of f(k) is
    #      -2 sum_{n >= k} v_n^2 u_n^2
    tail = np.cumsum((v * u)[::-1] ** 2)[::-1]  # tail[n] = sum_{m>=n} v_m^2 u_m^2
    for k in range(1, min(kmax, N - 1) + 1):
        c[k] -= 2.0 * tail[k]
    return c


def psi_from_coeffs(c, x):
    """psi(x) = <v, O[f_x] v>,  f_x(k) = x/(k^2+x)  (x = varpi^2)."""
    k = np.arange(1, len(c))
    x = np.atleast_1d(np.asarray(x, dtype=float))
    return (c[1:] * (x[:, None] / (k * k + x[:, None]))).sum(axis=1)


def concave_hull_at(xs, ys, xbar):
    """max over pairs (i,j) with x_i <= xbar <= x_j of the chord value at xbar."""
    best, arg = -np.inf, None
    lo = np.flatnonzero(xs <= xbar)
    hi = np.flatnonzero(xs >= xbar)
    for i in lo:
        for j in hi:
            if xs[j] == xs[i]:
                val, th = ys[i], 1.0
            else:
                th = (xs[j] - xbar) / (xs[j] - xs[i])
                val = th * ys[i] + (1 - th) * ys[j]
            if val > best:
                best, arg = val, (xs[i], xs[j], th)
    return best, arg


def ta_constructive(w_rms, Nc, mut, natoms, seed, ngrid=260):
    """For a random natoms-spectrum P build the two-atom P' from the proof and
    check  lam1(M_P') >= lam1(M_P)."""
    rng = np.random.default_rng(seed)
    xs = np.exp(rng.uniform(np.log(0.02), np.log(40.0), size=natoms))
    ps = rng.dirichlet(0.5 * np.ones(natoms))
    s = np.sqrt(w_rms**2 / float(ps @ (xs**2)))
    ws = xs * s
    N = int(max(300, 24 * ws.max(), 5 * Nc + 40))
    A = op.O_of_kernel(N, op.kernel_measure(N, ws, ps))
    lamP, vP, _ = op.top_pair(op.pencil(A, N, Nc, mut))

    c = struct_coeffs(vP, 2 * N)
    # cross-check the structure coefficients against a direct quadratic form
    xchk = 3.7
    direct = float(vP @ op.O_of_kernel(N, op.kernel_einstein(N, np.sqrt(xchk))) @ vP)
    chk(f"struct coeffs vs direct quad form (natoms={natoms},seed={seed})",
        abs(psi_from_coeffs(c, xchk)[0] - direct) < 1e-11 * max(1.0, abs(direct)),
        f"({psi_from_coeffs(c, xchk)[0]:.10e} vs {direct:.10e})")

    xbar = w_rms**2
    grid = np.unique(np.concatenate([
        np.exp(np.linspace(np.log(xbar * 1e-4), np.log(xbar * 1e4), ngrid)), [xbar]]))
    psi = psi_from_coeffs(c, grid)
    hull, arg = concave_hull_at(grid, psi, xbar)
    x1, x2, th = arg
    w2 = (np.sqrt(x1), np.sqrt(x2))
    p2 = (th, 1.0 - th)
    # ARROW 1 (the concave-hull / Jensen step the proof actually uses):
    #   <v_P, M_{P'} v_P>  >=  <v_P, M_P v_P> = lam1(M_P)
    # Both sides differ from psi only by the P-independent Coulomb term, so it
    # suffices to compare psi-values.
    psi_P = float(ps @ psi_from_coeffs(c, ws**2))
    arrow1 = hull >= psi_P - 1e-11 * max(1.0, abs(psi_P))
    N2 = int(max(300, 24 * max(w2), 5 * Nc + 40))
    A2 = op.O_of_kernel(N2, op.kernel_measure(N2, w2, p2))
    lamQ = np.linalg.eigvalsh(op.pencil(A2, N2, Nc, mut))[-1]
    return lamP, lamQ, (x1, x2, th), arrow1, (hull, psi_P)


def ta_table():
    print("\n### Theorem TA, constructive test: for a random k-atom P, build the")
    print("    two-atom P' that maximises the *same* Rayleigh quotient (concave hull")
    print("    of psi_v at xbar).  ARROW 1 = the hull/Jensen step; ARROW 2 = Rayleigh.")
    print(f"{'wrms':>6} {'Nc':>4} {'mut':>6} {'k':>3} {'lam1(M_P)':>13} {'lam1(M_P2)':>13} "
          f"{'gain':>11} {'x1':>9} {'x2':>10} {'theta':>8} {'2atom?':>7}")
    ndeg = 0
    for (w_rms, Nc, mut) in ((0.6, 1, 0.13), (0.6, 3, 0.13), (1.0, 6, 0.13),
                             (1.0, 3, 0.60), (2.0, 10, 0.30), (2.0, 6, 1.00)):
        for k in (3, 4, 6):
            for seed in (11, 12):
                lamP, lamQ, (x1, x2, th), arrow1, (hull, psi_P) = ta_constructive(
                    w_rms, Nc, mut, k, seed)
                tag = f"(w={w_rms},Nc={Nc},mut={mut},k={k},seed={seed})"
                chk("ARROW 1  hull(xbar) >= int psi dP  " + tag, arrow1,
                    f"({hull:.8e} vs {psi_P:.8e})")
                chk("ARROW 2  lam1(M_P') >= lam1(M_P)   " + tag,
                    lamQ >= lamP - 1e-10 * max(1.0, abs(lamP)),
                    f"(gain {lamQ - lamP:+.3e})")
                nd = 1e-9 < th < 1 - 1e-9
                ndeg += nd
                print(f"{w_rms:6.1f} {Nc:4d} {mut:6.2f} {k:3d} {lamP:13.8f} {lamQ:13.8f} "
                      f"{lamQ - lamP:+11.3e} {x1:9.4f} {x2:10.4f} {th:8.5f} "
                      f"{'yes' if nd else 'degen':>7}")
    print(f"  ({ndeg}/36 constructions are genuinely two-atom; the rest collapse to the")
    print("   single Einstein atom, which is still a legitimate <=2-atom competitor.)")


if __name__ == "__main__":
    sc_table()
    sc_case_check()
    ta_table()
    print(f"\nchecks: {CHK[0]-CHK[1]}/{CHK[0]} passed, {CHK[1]} failed")
    sys.exit(0 if CHK[1] == 0 else 1)
