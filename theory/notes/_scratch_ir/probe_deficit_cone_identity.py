"""Explicit characterisation of the PHYSICAL Jensen-deficit cone, and the
second-order structure of Theorem E at mu*>0.

(1) DEFICIT IDENTITY.  For P = w d_{x1} + (1-w) d_{x2} with mean xbar,
    write s = k^2, V = w(1-w)(x1-x2)^2 = Var(P).  Then

        delta(s) := f_xbar(s) - F_P(s)  =  s * V / [(s+x1)(s+x2)(s+xbar)]

    All s^1 and s^0 terms in the numerator cancel identically; only the
    VARIANCE survives.  Combined with Bauer (delta_P is AFFINE in P on the
    set {mass 1, mean xbar}, whose extreme points are 2-atom measures), this
    gives the FULL explicit description of the physical cone:

        D_xbar = conv { V * psi_{x1,x2} },   psi = k^2/[(k^2+x1)(k^2+x2)(k^2+xbar)]

    This replaces "delta >= 0" (for which E4 is FALSE, see probe_E4_delta_k)
    by a 2-parameter family of explicit rational profiles.

(2) SECOND-ORDER STRUCTURE.  Let M(x) = O[f_x] - 2 mu_t |u_c><u_c|, with top
    pair (lam1, v).  Since f_x(0)=1 is x-independent,

        dM/dx  = O[k^2/(k^2+x)^2],      d2M/dx2 = -2 O[psi_0],
        psi_0(k) = k^2/(k^2+x)^3   (= the x1=x2=xbar member of the family).

    Define the EXTREMALITY FUNCTIONAL   Psi := <v, O[psi_0] v>.
    Second-order extremality of the Einstein point at mu*>0 is exactly Psi >= 0.

    Standard perturbation theory gives
        H''(x) = <v, M'' v> + 2 sum_{j!=1} |<v_j, M' v>|^2 / (lam1 - lam_j)
               = -2 Psi + R,        R >= 0  (level repulsion, always >= 0).

    So H'' > 0 does NOT contradict extremality: the small-varpi "convexity
    window" can come entirely from R.  This script measures the split.

Strictness: numerical-evidence (float64, finite basis).
"""
import os, math
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh

N = 300
_n = np.arange(N)
_u = 1.0 / np.sqrt(2.0 * _n + 1.0)
_ks = np.arange(2 * N + 1, dtype=float)
_absd = np.abs(_n[:, None] - _n[None, :])
_sumd = _n[:, None] + _n[None, :] + 1
_uu = _u[:, None] * _u[None, :]


def O(f):
    """Operator O[f] for an arbitrary channel profile f (len >= 2N+1)."""
    A = (f[_absd] + f[_sumd]) * _uu
    L1 = np.zeros(N)
    L1[1:] = np.cumsum(2.0 * f[1:N])[:N - 1] / (2 * _n[1:] + 1)
    np.fill_diagonal(A, -L1 + f[2 * _n + 1] * _u ** 2)
    return A


def f_ein(x):
    f = x / (_ks ** 2 + x)
    f[0] = 1.0
    return f


def F_meas(atoms):
    f = np.zeros(2 * N + 1)
    for w, x in atoms:
        f += w * x / (_ks ** 2 + x)
    f[0] = sum(w for w, _ in atoms)
    return f


def coulomb(mu_t, Nc):
    uc = np.where(_n < (Nc if Nc else N), _u, 0.0)
    return 2.0 * mu_t * np.outer(uc, uc)


# ---------------------------------------------------------------- part 1
print("=" * 74)
print("(1) DEFICIT IDENTITY   delta(k) = V k^2 / [(k^2+x1)(k^2+x2)(k^2+xbar)]")
print("=" * 74)
print(f"{'w':>6} {'x1':>9} {'x2':>9} {'xbar':>9} {'V':>11} {'max|err|':>12}")
rng = np.random.default_rng(20260917)
worst = 0.0
for _ in range(8):
    w = float(rng.uniform(0.05, 0.95))
    x1 = float(np.exp(rng.uniform(math.log(1e-3), math.log(300.0))))
    x2 = float(np.exp(rng.uniform(math.log(1e-3), math.log(300.0))))
    xbar = w * x1 + (1 - w) * x2
    V = w * (1 - w) * (x1 - x2) ** 2
    direct = f_ein(xbar) - F_meas([(w, x1), (1 - w, x2)])
    s = _ks ** 2
    formula = s * V / ((s + x1) * (s + x2) * (s + xbar))
    err = float(np.max(np.abs(direct - formula)))
    worst = max(worst, err)
    print(f"{w:>6.3f} {x1:>9.4f} {x2:>9.4f} {xbar:>9.4f} {V:>11.4g} {err:>12.3e}")
print(f"\nworst deviation over 8 random 2-atom measures: {worst:.3e}"
      f"   ({'IDENTITY HOLDS' if worst < 1e-12 else 'IDENTITY FAILS'})")

# ---------------------------------------------------------------- part 2
print()
print("=" * 74)
print("(2) H''(x) = -2*Psi + R    [Psi = <v,O[psi_0]v>,  R = level repulsion]")
print("    Psi >= 0  <=>  Einstein point is a 2nd-order max (Theorem E local)")
print("=" * 74)
print(f"{'x':>8} {'varpi':>7} {'mu_t':>6} {'Nc':>5} {'gap l1-l2':>10} "
      f"{'Psi':>11} {'-2Psi':>11} {'R':>11} {'-2Psi+R':>11} {'H\" (FD)':>11}")

for mu_t, Nc in ((0.0, None), (0.05, 6), (0.13, 6), (0.13, 45), (0.26, 6)):
    C = coulomb(mu_t, Nc) if mu_t else 0.0
    for x in (0.01, 0.09, 0.1089, 0.16, 0.25, 1.0, 9.0, 100.0):
        s = _ks ** 2
        M = O(f_ein(x)) - C
        lam, V_ = eigh(M)
        v = V_[:, -1]
        Mp = O(s / (s + x) ** 2)                 # dM/dx
        psi0 = s / (s + x) ** 3
        Psi = float(v @ (O(psi0) @ v))           # <v, O[psi_0] v>
        # level-repulsion term
        c = V_[:, :-1].T @ (Mp @ v)              # <v_j, M' v>, j != top
        R = float(2.0 * np.sum(c ** 2 / (lam[-1] - lam[:-1])))
        # independent finite-difference H''
        dx = x / 50.0
        vals = []
        for i in (-2, -1, 0, 1, 2):
            Mi = O(f_ein(max(x + i * dx, 1e-12))) - C
            vals.append(float(eigh(Mi, eigvals_only=True)[-1]))
        Hpp = (-vals[0] + 16 * vals[1] - 30 * vals[2]
               + 16 * vals[3] - vals[4]) / (12 * dx ** 2)
        flag = "  <<< Psi<0" if Psi < 0 else ""
        print(f"{x:>8.4f} {math.sqrt(x):>7.4f} {mu_t:>6.2f} {str(Nc):>5} "
              f"{lam[-1]-lam[-2]:>10.4f} {Psi:>11.4e} {-2*Psi:>11.4e} "
              f"{R:>11.4e} {-2*Psi+R:>11.4e} {Hpp:>11.4e}{flag}")
    print()
