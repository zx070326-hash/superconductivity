"""Theorem CR ("cone repair") for mu*>0 Einstein extremality -- measurement.

Setting (P4D_coulomb_operator.md sec.1, P4_main_theorems.md sec.0/sec.2):
  A_P  = O[F_P],  M_P = A_P - 2 mut u_c u_c^T,  mut = mu*/lambda,  u_c = P_{Nc} u
  d    = [[.]](w_rms) - F_P >= 0   (kernel-level Jensen, analytic)
  O[d] = H(w_rms) - A_P            (O is linear in the kernel)
  Extremality  <==  <v_P, O[d] v_P> >= 0   (Rayleigh step; Coulomb cancels)

Abel sign split of the top eigenvector v_P = w - xi with w, xi in the cone C:
  a_L := (v_P)_{L-1} - (v_P)_L,   xi := sum_{a_L<0} |a_L| e_L,  w := v_P + xi
Corollary A' / Lemma A give  G11 := <w,O[d]w> >= 0,  G22 := <xi,O[d]xi> >= 0,
G12 := <xi,O[d]w> >= 0.  Hence

  (CS)   G12^2 <= G11 * G22   ==>   <v_P,O[d]v_P> = G11 - 2 G12 + G22 >= 0.

Measures the mut-threshold of: cone / CS / exact Rayleigh / true extremality.
float64, finite N  ==>  numerical-evidence.
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


# ------------------------------------------------------------ model pieces ---
def two_atom(w_rms, weight, spread):
    """P = weight*delta_{w1} + (1-weight)*delta_{w2}, w2/w1 = spread, <w^2> = w_rms^2."""
    r = spread
    # weight*w1^2 + (1-weight)*(r w1)^2 = w_rms^2
    w1 = w_rms / np.sqrt(weight + (1.0 - weight) * r * r)
    return (w1, r * w1), (weight, 1.0 - weight)


def nsize(w_rms, Nc, ws):
    return int(max(300, 24 * max(ws), 5 * Nc + 40))


def build(w_rms, ws, ps, Nc, mut, N):
    A = op.O_of_kernel(N, op.kernel_measure(N, ws, ps))
    H = op.O_of_kernel(N, op.kernel_einstein(N, w_rms))
    Od = H - A
    MP = op.pencil(A, N, Nc, mut)
    ME = op.pencil(H, N, Nc, mut)
    return A, H, Od, MP, ME


def diagnose(w_rms, ws, ps, Nc, mut, N):
    A, H, Od, MP, ME = build(w_rms, ws, ps, Nc, mut, N)
    lamP, vP, _ = op.top_pair(MP)
    lamE = np.linalg.eigvalsh(ME)[-1]
    a = op.abel(vP)
    neg = a < 0.0
    xi_a = np.where(neg, -a, 0.0)
    xi = op.from_abel(xi_a)
    w = vP + xi
    G11 = float(w @ Od @ w)
    G12 = float(xi @ Od @ w)
    G22 = float(xi @ Od @ xi)
    exact = float(vP @ Od @ vP)
    return dict(
        lamP=lamP,
        lamE=lamE,
        extr_ok=(lamP <= lamE),
        in_cone=(not neg.any()),
        n_neg=int(neg.sum()),
        neg_idx=np.flatnonzero(neg)[:5] + 1,  # Abel labels L
        G11=G11,
        G12=G12,
        G22=G22,
        CS_ok=(G12 * G12 <= G11 * G22),
        exact_ok=(exact >= 0.0),
        exact=exact,
        vP=vP,
    )


# ------------------------------------------------------ threshold scanning ---
MUGRID = np.geomspace(1e-5, 8.0, 121)


PREDS = (
    ("cone", lambda r: r["in_cone"]),
    ("rayl", lambda r: r["exact_ok"]),
    ("extr", lambda r: r["extr_ok"]),
)


def thresholds(w_rms, weight, spread, Nc, mugrid=None):
    """First mut on the log grid at which each criterion fails (inf = never)."""
    ws, ps = two_atom(w_rms, weight, spread)
    N = nsize(w_rms, Nc, ws)
    out = {name: np.inf for name, _ in PREDS}
    pending = {name for name, _ in PREDS}
    for mut in (MUGRID if mugrid is None else mugrid):
        r = diagnose(w_rms, ws, ps, Nc, mut, N)
        for name, pred in PREDS:
            if name in pending and not pred(r):
                out[name] = mut
                pending.discard(name)
        if not pending:
            break
    out["N"] = N
    return out


SHAPES = ((0.5, 1.5), (0.9, 3.0), (0.99, 10.0), (0.99, 30.0), (0.998, 100.0))


def main():
    print("=" * 104)
    print("Three nested mut-thresholds for mu*>0 Einstein extremality (two-atom spectra)")
    print("  cone : v_P leaves the cone C          -> repo Proposition P3 hypothesis dies")
    print("  rayl : <v_P,O[d]v_P> < 0              -> the Rayleigh step itself dies")
    print("  extr : lam1(M_P) > lam1(M_E)          -> extremality actually fails")
    print("Always cone <= rayl <= extr.  Reported per (wrms,Nc) as the inf over shapes")
    print("(the binding shape), since a criterion must hold for EVERY P.")
    print("Log grid, 121 points, 1e-5..8; '>8' = never on the grid.  numerical-evidence.")
    print("=" * 104)
    print(f"{'wrms':>5} {'Nc':>4} {'N':>5} {'cone':>10} {'rayl':>10} {'extr':>10} "
          f"{'rayl/cone':>10} {'extr/rayl':>10} {'extr/cone':>10}")
    rows = []
    for w_rms in (0.4, 0.6, 1.0, 2.0):
        for Nc in (1, 3, 6, 10, 20):
            agg = {k: np.inf for k in ("cone", "rayl", "extr")}
            N = 0
            for (weight, spread) in SHAPES:
                t = thresholds(w_rms, weight, spread, Nc)
                N = max(N, t["N"])
                for k in agg:
                    agg[k] = min(agg[k], t[k])
            rows.append((w_rms, Nc, agg))

            def f(x):
                return "      >8" if not np.isfinite(x) else f"{x:10.3e}"

            def g(a, b):
                return "       -" if not (np.isfinite(a) and np.isfinite(b)) else f"{a/b:10.1f}"
            print(f"{w_rms:5.1f} {Nc:4d} {N:5d} {f(agg['cone'])} {f(agg['rayl'])} "
                  f"{f(agg['extr'])} {g(agg['rayl'], agg['cone'])} "
                  f"{g(agg['extr'], agg['rayl'])} {g(agg['extr'], agg['cone'])}")
    return rows


if __name__ == "__main__":
    main()
